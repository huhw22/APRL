#!/usr/bin/env python3
"""Check output transactions of the single retirement production run."""
import json
from pathlib import Path
import sys
import h5py
import yaml

root = Path(sys.argv[1]).resolve()
card = yaml.safe_load((root / "run.yaml").read_text())
log = (root / "aprl.log").read_text(errors="replace")
if "Configured stop reached:" not in log:
    raise SystemExit("No configured physical stop in APRL log")
products = []
for kind, group in (("field_planes", "field_plane"), ("particle_planes", "particle_plane")):
    for plane in card["detectors"][kind]:
        products.append((root / "detectors" / (plane["name"] + ".h5"), group))
products.append((root / "ledger" / "energy-ledger.h5", "energy_ledger"))
checks = []
for path, name in products:
    with h5py.File(path, "r") as file:
        g = file[name]
        field = name == "field_plane"
        n = int(g["committed_samples" if field else "committed_records"][()])
        data = g["fields" if field else "records"]
        if int(g["complete"][()]) != 1 or n <= 0 or len(data) != n:
            raise SystemExit(f"Output incomplete: {path}")
        if field and len(g["time_s"]) != n:
            raise SystemExit(f"Field time/data length mismatch: {path}")
        checks.append(dict(file=str(path), group=name, records=n, complete=True))
(root / "completion-check.json").write_text(json.dumps(
    dict(complete=True, scientific_validation=False, products=checks), indent=2) + "\n")
