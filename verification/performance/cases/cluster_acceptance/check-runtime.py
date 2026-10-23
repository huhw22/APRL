#!/usr/bin/env python3
"""Check immutable runtime dependencies and audit node-local installation differences."""
import hashlib
import json
import os
from pathlib import Path
import platform
import re
import shutil
import subprocess
import sys


def digest(path):
    value = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            value.update(chunk)
    return value.hexdigest()


def command(arguments):
    result = subprocess.run(arguments, text=True, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, check=False)
    return result.returncode, result.stdout


def version_line(name, text):
    pattern = r"^slurm\s+\S+\s*$" if name in ("srun", "scontrol") else r"^mpirun \(Open MPI\)\s+\S+\s*$"
    matches = [line.strip() for line in text.splitlines() if re.match(pattern, line.strip())]
    if len(matches) != 1:
        raise ValueError(f"Cannot identify {name} version: {text.strip()}")
    return matches[0]


def libraries(text):
    result = {}
    for line in text.splitlines():
        match = re.match(r"\s*(\S+)\s+=>\s+(/\S+)", line)
        if match:
            result[match[1]] = match[2]
        else:
            match = re.match(r"\s*(/\S+)\s+\(", line)
            if match:
                result[Path(match[1]).name] = match[1]
    return result


def check(record, run=command, hash_file=digest):
    audit = {"host": platform.node(), "node_local_files": {}, "tool_versions": {},
             "executables": {}, "errors": [], "pass": False}
    errors = audit["errors"]
    strict = record["strict_files"]
    cache = {}

    def fingerprint(path):
        if path not in cache:
            cache[path] = {"realpath": str(Path(path).resolve()), "sha256": hash_file(path)}
        return cache[path]

    for path, reference in record["node_local_files"].items():
        try:
            current = fingerprint(path)
            audit["node_local_files"][path] = {**current, "reference": reference,
                "changed": current["sha256"] != reference["sha256"] or current["realpath"] != reference["realpath"]}
        except OSError as exc:
            # A node may provide the same SONAME under another system path.
            audit["node_local_files"][path] = {"reference": reference, "changed": True, "error": str(exc)}
    if str(Path(os.environ.get("APRL_DEPENDENCY_PREFIX", "")).resolve()) != record["dependency_prefix"]:
        errors.append("Dependency prefix differs from the prepared build")
    if Path(sys.executable).resolve() != Path(record["python"]).resolve():
        errors.append("Python differs from the prepared analysis environment")
    for name, expected in record["tools"].items():
        actual = shutil.which(name)
        if actual is None or Path(actual).absolute() != Path(expected):
            errors.append(f"Tool selection changed: {name}: {actual} != {expected}")
    for name, expected in record["tool_versions"].items():
        try:
            status, text = run([record["tools"][name], "--version"])
            actual = version_line(name, text)
            audit["tool_versions"][name] = {"version": actual, "output": text, "returncode": status}
            if status or actual != version_line(name, expected):
                errors.append(f"Incompatible tool version: {name}: {actual}")
        except (OSError, ValueError) as exc:
            errors.append(str(exc))
    if all(name in audit["tool_versions"] for name in ("srun", "scontrol")):
        if audit["tool_versions"]["srun"]["version"] != audit["tool_versions"]["scontrol"]["version"]:
            errors.append("srun and scontrol versions differ")

    for executable, reference in record["executables"].items():
        details = {}
        audit["executables"][executable] = details
        try:
            status, text = run(["ldd", "-r", executable])
            actual = libraries(text)
            details.update(returncode=status, ldd=text, libraries=actual)
            if status or re.search(r"not found|undefined symbol|error while loading|cannot open shared object", text, re.I):
                errors.append(f"Unresolved runtime dependencies: {executable}")
            for soname, expected_path in libraries(reference["ldd"]).items():
                # MPI and shared/private dependencies are pinned by the strict manifest.
                if expected_path not in strict:
                    continue
                current_path = actual.get(soname)
                if current_path is None:
                    errors.append(f"Missing pinned dependency: {executable}: {soname}")
                elif Path(current_path).resolve() != Path(expected_path).resolve() or fingerprint(current_path)["sha256"] != strict[expected_path]:
                    errors.append(f"Pinned dependency changed: {executable}: {soname}")
            details["library_fingerprints"] = {name: fingerprint(path) for name, path in actual.items()}
        except (OSError, ValueError) as exc:
            errors.append(f"Cannot inspect runtime for {executable}: {exc}")

    # Each scientific executable must use the same actual MPI and parallel HDF5.
    workspace = Path(record["workspace"])
    scientific = [str(workspace / path) for path in (
        "build/baseline/aprl", "build/candidate/aprl", "build/candidate/sase_bunch_generator")]
    for family, pattern in (("MPI", r"^libmpi\.so(?:\.|$)"), ("HDF5", r"^libhdf5\.so(?:\.|$)")):
        identities = []
        for executable in scientific:
            details = audit["executables"].get(executable, {})
            selected = [value for name, value in details.get("library_fingerprints", {}).items() if re.match(pattern, name)]
            if len(selected) != 1:
                errors.append(f"Expected exactly one {family} runtime: {executable}")
            else:
                identities.append((selected[0]["realpath"], selected[0]["sha256"]))
        if len(identities) == 3 and len(set(identities)) != 1:
            errors.append(f"Baseline, candidate and generator use different {family} runtimes")
    audit["pass"] = not errors
    return audit


def main():
    if len(sys.argv) != 3:
        raise SystemExit("Usage: check-runtime.py BUILD_RECORD_JSON AUDIT_JSON")
    try:
        audit = check(json.loads(Path(sys.argv[1]).read_text()))
    except Exception as exc:
        audit = {"host": platform.node(), "pass": False, "errors": [f"Runtime preflight failed: {exc}"]}
    with Path(sys.argv[2]).open("x") as stream:
        json.dump(audit, stream, indent=2)
        stream.write("\n")
    for error in audit["errors"]:
        print("ERROR:", error, file=sys.stderr)
    print("Runtime compatibility:", "PASS" if audit["pass"] else "FAIL")
    return 0 if audit["pass"] else 1


if __name__ == "__main__":
    sys.exit(main())
