"""Verify the APRL 100k synchrotron inputs and prepare a portable production card."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import tempfile

HERE = Path(__file__).resolve().parent
METHODS = ("osiris", "left-20056")
FIELDS = ("Ex", "Ey", "Ez", "Bx", "By", "Bz")
MEMORIES = {"osiris": (0, 3, 4, 5, 6, 9, 10, 11), "left-20056": (3, 4, 9, 10)}
BEAM_SHA256 = "4a7e3b4ba16078adfaf76000e82024e9a73182ab1b71340017d55b718b89c687"


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as reader:
        for block in iter(lambda: reader.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def components(method):
    return (*FIELDS, *(f"psi{index}" for index in MEMORIES[method]))


def verify(data, methods):
    checks = {}
    for row in (data / "SHA256SUMS").read_text().splitlines():
        digest, name = row.split("  ", 1)
        if name in checks or Path(name).is_absolute() or ".." in Path(name).parts:
            raise ValueError(f"Invalid SHA256SUMS member: {name}")
        checks[name] = digest
    required = {"beam.h5", "run-20104.yaml.in", "prepare.py", "README.zh-CN.md", "DATA.json"}
    for method in methods:
        files = ["startup.yaml", "READY.json", *(name + ".bin" for name in components(method))]
        required.update(method + "/" + name for name in files)
    if checks.get("beam.h5") != BEAM_SHA256:
        raise ValueError("This example requires the frozen 100000-electron beam")
    for name in sorted(required):
        path = data / name
        if not path.is_file() or sha256(path) != checks.get(name):
            raise ValueError(f"Startup package missing or checksum mismatch: {name}")
    if sha256(HERE / "run-20104.yaml.in") != checks["run-20104.yaml.in"]:
        raise ValueError("The installed example card and supplied data package have different revisions")
    return len(required)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, required=True, help="Extracted beam-cpml-100k directory")
    parser.add_argument("--method", choices=METHODS, default="osiris")
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument("--verify", action="store_true", help="Check both methods and the beam without creating a run")
    action.add_argument("--output", type=Path, help="New directory for the 20104 physics run")
    args = parser.parse_args()
    data = args.data.resolve()
    count = verify(data, METHODS if args.verify else (args.method,))
    if args.verify:
        print(f"Verified {count} files for osiris and left-20056.")
        return
    output = args.output.resolve()
    if output.exists():
        parser.error(f"Run directory already exists: {output}")
    card = (HERE / "run-20104.yaml.in").read_text()
    for token, value in (("METHOD", args.method), ("ARRAY_DIR", os.path.relpath(data / args.method, output)),
                         ("BEAM_FILE", os.path.relpath(data / "beam.h5", output))):
        marker = "'@" + token + "@'"
        if card.count(marker) != 1:
            raise ValueError(f"Invalid template substitution: {token}")
        card = card.replace(marker, json.dumps(value))
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=output.parent) as staging:
        staged = Path(staging) / "run"
        staged.mkdir()
        (staged / "run.yaml").write_text(card)
        staged.rename(output)
    print(f"Verified {count} input files. Geometry and beam are set to the 20104 benchmark.")
    print(f"Run from {output} with the installed APRL 2.1 executable and run.yaml.")


if __name__ == "__main__":
    main()
