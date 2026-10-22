#!/usr/bin/env python3
"""Prepare held-field 1 m aperture/CPML/dz checks from one existing HDF5 bunch.

The reference held-field result is supplied separately to the analysis script.
This command only writes APRL cards; it neither copies nor regenerates particles.
"""

import argparse
import json
import math
from pathlib import Path

from prepare import card


# The physical half aperture is (xy_cells / 2 - xy_cpml_cells) * 5 um.
# Doubling dz resolution preserves physical z-CPML/guard distances.
VARIANTS = {
    "base": {},
    "wide": {"xy_cells": 128},
    "thick-cpml": {"xy_cells": 124, "xy_cpml_cells": 20},
    "long-guard": {"longitudinal_guard_cells": 768},
    "fine-z": {
        "dz_um": 2.5,
        "z_cpml_cells": 24,
        "longitudinal_guard_cells": 768,
        "longitudinal_alignment_cells": 512,
    },
}
DEFAULT_VARIANTS = "wide,thick-cpml,fine-z"


def parse_variants(value):
    names = [part.strip() for part in value.split(",")]
    if not names or any(not name for name in names):
        raise ValueError("--variants must be a nonempty comma-separated list")
    if len(set(names)) != len(names):
        raise ValueError("--variants contains duplicates")
    invalid = [name for name in names if name not in VARIANTS]
    if invalid:
        raise ValueError(
            "unknown variant(s): " + ", ".join(invalid)
            + "; available: " + ", ".join(VARIANTS)
        )
    return names


def variant_geometry(name):
    overrides = VARIANTS[name]
    xy_cells = overrides.get("xy_cells", 104)
    xy_cpml = overrides.get("xy_cpml_cells", 10)
    return {
        "xy_cells": xy_cells,
        "xy_cpml_cells": xy_cpml,
        "physical_half_aperture_um": (xy_cells / 2 - xy_cpml) * 5,
        "dz_um": overrides.get("dz_um", 5.0),
        "z_cpml_cells": overrides.get("z_cpml_cells", 12),
        "longitudinal_guard_cells": overrides.get("longitudinal_guard_cells", 384),
        "longitudinal_alignment_cells": overrides.get("longitudinal_alignment_cells", 256),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path,
                        help="existing 10 um HDF5-v4 particle file; not copied")
    parser.add_argument("--output-dir", required=True, type=Path,
                        help="suite directory; each variant gets its own subdirectory")
    parser.add_argument("--electrons", required=True, type=float,
                        help="represented physical electrons, e.g. 10409720 at weight 1")
    parser.add_argument("--variants", default=DEFAULT_VARIANTS,
                        help="comma-separated: " + ", ".join(VARIANTS))
    args = parser.parse_args()

    input_file = args.input.resolve()
    output_dir = args.output_dir.resolve()
    if not input_file.is_file():
        parser.error("--input is not a file: " + str(input_file))
    if not math.isfinite(args.electrons) or args.electrons <= 0:
        parser.error("--electrons must be positive and finite")
    try:
        names = parse_variants(args.variants)
    except ValueError as exc:
        parser.error(str(exc))
    targets = [output_dir / name / "input.yaml" for name in names]
    existing = [target for target in targets if target.exists()]
    if existing:
        parser.error("refusing to replace existing card: " + str(existing[0]))
    manifest_path = output_dir / "held-convergence-plan.json"
    if manifest_path.exists():
        parser.error("refusing to replace existing plan: " + str(manifest_path))

    plan = {
        "purpose": "held-field startup/field-transient sensitivity; not a physical-truth comparison",
        "input": str(input_file),
        "electrons": args.electrons,
        "profile": "g929-1m",
        "case": "held-field",
        "variants": {},
        "notes": [
            "All variants retain the same 1 m laboratory observation plane and particle source.",
            "Use a fixed central spatial ROI and common laboratory-time interval in comparisons.",
            "after-last-element stops after particle passage; check whether the field tail is truncated.",
        ],
    }
    for name, target in zip(names, targets):
        run_dir = target.parent / "runs" / "g929-1m-held-field"
        plan["variants"][name] = {
            "card": str(target),
            "run_dir": str(run_dir),
            **variant_geometry(name),
        }

    for name, target in zip(names, targets):
        target.parent.mkdir(parents=True, exist_ok=True)
        contents = card(
            "g929-1m", "held-field", input_file, args.electrons,
            target.parent, True, **VARIANTS[name]
        )
        target.write_text(
            f"# Held convergence variant: {name}\n" + contents,
            encoding="utf-8",
        )
        print(f"{name}: {target} -> {plan['variants'][name]['run_dir']}")
    manifest_path.write_text(json.dumps(plan, indent=2, ensure_ascii=False) + "\n",
                             encoding="utf-8")
    print("Plan: " + str(manifest_path))


if __name__ == "__main__":
    main()
