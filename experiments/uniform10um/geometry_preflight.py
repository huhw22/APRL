#!/usr/bin/env python3
"""Prepare or run a bounded high-boost t=0 aperture check.

This is geometry only. It turns off initial fields, feedback, detectors, and
the ledger. With --run, APRL runs directly in interactive mode and receives
SIGINT as soon as all startup geometry/field initialization has completed.
The duration guard is one Maxwell step, so it cannot accidentally become a
full high-boost simulation. Do not run the generated cards directly.
"""

import argparse
import json
import re
import signal
import subprocess
from pathlib import Path


SOURCE = Path(__file__).resolve().parent
EXACT_FRACTION = re.compile(
    r"boosted t=0 transverse ballistic positions outside the physical "
    r"aperture: (\d+)/(\d+) \(([^%]+)%\); outside the full mesh: "
    r"(\d+)/(\d+) \(([^%]+)%\)"
)
BOOSTED_RANGE = re.compile(r"centred boosted z range \[m\]=\[([^]]+)\]")


def set_scalar(text, key, replacement):
    pattern = re.compile(r"(?m)^(\s*" + re.escape(key) + r":\s*).*$")
    changed, count = pattern.subn(lambda match: match.group(1) + replacement,
                                  text, count=1)
    if count != 1:
        raise ValueError(f"template must contain exactly one {key} key")
    return changed


def prepare(input_file, electrons, output_dir, modes, overwrite):
    card_dir = output_dir / "cards"
    card_dir.mkdir(parents=True, exist_ok=True)
    result = []
    for mode in modes:
        source = SOURCE / f"preflight-g929-{mode}-w1.yaml"
        card = card_dir / f"g929-{mode}-geometry.yaml"
        run_dir = output_dir / "runs" / f"g929-{mode}-geometry"
        if card.exists() and not overwrite:
            raise FileExistsError(f"refusing to replace {card}")
        text = source.read_text(encoding="utf-8")
        text = set_scalar(text, "manifest", json.dumps(
            str(run_dir / "run-manifest.yaml")))
        text = set_scalar(text, "file", json.dumps(str(input_file)))
        text = set_scalar(text, "electrons", f"{electrons:.17g}")
        card.write_text(text, encoding="utf-8")
        result.append((mode, card, run_dir))
        print(card)
    return result


def run_one(executable, mode, card, run_dir):
    run_dir.mkdir(parents=True, exist_ok=True)
    log_path = run_dir / "run.log"
    if log_path.exists():
        raise FileExistsError(f"refusing to replace {log_path}")
    count = None
    boosted_range = None
    marker_seen = False
    stop_sent = False
    with log_path.open("w", encoding="utf-8") as log:
        with subprocess.Popen([str(executable), str(card)],
                              stdout=subprocess.PIPE,
                              stderr=subprocess.STDOUT,
                              text=True, bufsize=1) as process:
            assert process.stdout is not None
            for line in process.stdout:
                log.write(line)
                log.flush()
                match = EXACT_FRACTION.search(line)
                if match:
                    count = {
                        "outside_physical": int(match.group(1)),
                        "total_macroparticles": int(match.group(2)),
                        "outside_physical_percent": float(match.group(3)),
                        "outside_mesh": int(match.group(4)),
                        "outside_mesh_percent": float(match.group(6)),
                    }
                match = BOOSTED_RANGE.search(line)
                if match:
                    boosted_range = match.group(1)
                if ("EXPERIMENTAL frozen lab-plane injection:" in line and
                        not stop_sent):
                    marker_seen = True
                    process.send_signal(signal.SIGINT)
                    stop_sent = True
            code = process.wait()
    output = {"mode": mode, "card": str(card), "run_log": str(log_path),
              "exit_code": code, "initialized": marker_seen,
              "boosted_z_range_m": boosted_range,
              "exact_aperture": count}
    print(json.dumps(output, indent=2))
    if code != 0 or not marker_seen or count is None:
        raise RuntimeError(f"geometry preflight failed; inspect {log_path}")
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True,
                        help="HDF5-v4 input; server may regenerate it")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--electrons", type=float, default=10409720.0)
    parser.add_argument("--modes", default="held,ballistic",
                        help="comma-separated subset of held,ballistic")
    parser.add_argument("--overwrite-cards", action="store_true")
    parser.add_argument("--run", action="store_true",
                        help="run single-rank APRL and stop at initialization")
    parser.add_argument("--aprl", type=Path,
                        default=SOURCE.parents[1] / "build/aprl")
    args = parser.parse_args()
    if not args.input.is_file():
        parser.error("--input HDF5 does not exist")
    if args.electrons <= 0:
        parser.error("--electrons must be positive")
    modes = [item.strip() for item in args.modes.split(",") if item.strip()]
    if not modes or len(set(modes)) != len(modes) or any(
            item not in ("held", "ballistic") for item in modes):
        parser.error("--modes must be held,ballistic or a unique subset")
    if args.run and not args.aprl.is_file():
        parser.error("APRL executable not found; pass --aprl")
    source_cpp = SOURCE.parents[1] / "src/simulation.cpp"
    if (args.run and source_cpp.is_file() and
            args.aprl.stat().st_mtime < source_cpp.stat().st_mtime):
        parser.error("APRL executable predates src/simulation.cpp; rebuild before preflight")
    output_dir = args.output_dir.resolve()
    prepared = prepare(args.input.resolve(), args.electrons, output_dir,
                       modes, args.overwrite_cards)
    if args.run:
        for mode, card, run_dir in prepared:
            run_one(args.aprl.resolve(), mode, card, run_dir)
    else:
        print("Cards only. Use --run to stop APRL cleanly after startup; "
              "direct execution would hit its one-step duration guard.")


if __name__ == "__main__":
    main()
