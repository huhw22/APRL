#!/usr/bin/env python3
"""Record an already built campaign, without running binaries or submitting jobs."""
import hashlib
import json
import os
from pathlib import Path
import platform
import re
import shutil
import subprocess
import sys
from datetime import datetime, timezone


def output(*command):
    return subprocess.check_output(command, text=True, stderr=subprocess.STDOUT).strip()


def digest(path):
    value = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            value.update(chunk)
    return value.hexdigest()


def main():
    workspace = Path(sys.argv[1]).resolve()
    prepared = workspace / "tests/cluster_acceptance"
    record_path = prepared / "build-record.json"
    manifest_path = prepared / "manifest.sha256"
    if record_path.exists() or manifest_path.exists():
        raise SystemExit("Existing campaign manifest: preserve it; use a new prepared campaign for rebuilding.")
    paths = {workspace / "cluster.env"}
    record = {"created_utc": datetime.now(timezone.utc).isoformat(), "host": platform.node(),
              "workspace": str(workspace), "builds": {}, "executables": {},
              "cmake": output("cmake", "--version"), "compiler": output("/usr/bin/g++", "--version"),
              "mpi": output("/usr/bin/mpirun", "--version"),
              "mpi_wrapper": output("/usr/bin/mpic++", "--showme"),
              "python": sys.executable,
              "status": "compiled and statically checked; no Slurm job submitted or scientific test run"}
    prefix = Path(os.environ["APRL_DEPENDENCY_PREFIX"]).resolve()
    record["dependency_prefix"] = str(prefix)
    paths.add(prefix / "share/aprl/runtime-env.sh")
    record["tools"] = {}
    for name in ("srun", "scontrol", "ctest", "cmake", "taskset", "mpirun"):
        selected = shutil.which(name)
        if not selected:
            raise SystemExit(f"Missing tool: {name}")
        record["tools"][name] = str(Path(selected).absolute())
        paths.add(Path(selected).absolute())
    import h5py
    import numpy
    record["python_libraries"] = {"numpy": numpy.__version__, "h5py": h5py.__version__,
                                  "h5py_hdf5": h5py.version.hdf5_version}
    for name, source_name in (("baseline", "aprl-baseline"), ("candidate", "aprl"), ("verification", "aprl")):
        build = workspace / "build" / name
        source = workspace / "source" / source_name
        flags = build / "CMakeFiles/direct_eb_core.dir/flags.make"
        text = flags.read_text()
        match = re.search(r"APRL_SOURCE_REVISION=.*?([0-9a-f]{40})(-dirty)?", text)
        if not match or match.group(2):
            raise SystemExit(f"Missing clean core revision in {flags}")
        revision = match.group(1)
        # Check the numerical implementation still matches the captured build.
        subprocess.run(["git", "-C", str(source), "diff", "--exit-code", revision, "--",
                        "src", "utilities", "CMakeLists.txt", "postprocess"], check=True)
        record["builds"][name] = {
            "source": str(source), "compiled_core_revision": revision,
            "current_source_head": output("git", "-C", str(source), "rev-parse", "HEAD"),
            "source_status": output("git", "-C", str(source), "status", "--short"),
            "flags": text, "cache": str(build / "CMakeCache.txt")}
        paths.update((flags, build / "CMakeCache.txt"))
        paths.update(build.rglob("CTestTestfile.cmake"))
        for path in build.rglob("*"):
            if "CMakeFiles" in path.parts or not path.is_file() or not os.access(path, os.X_OK):
                continue
            with path.open("rb") as stream:
                if stream.read(4) == b"\x7fELF":
                    paths.add(path)
        if name == "verification":
            paths.update(p for p in (build / "verification").glob("*.yaml"))
            paths.update(p for p in (build / "verification").glob("*.cmake"))
    paths.add(workspace / "build/cluster-tools/binding_probe")
    for required in ("build/baseline/aprl", "build/candidate/aprl", "build/candidate/sase_bunch_generator",
                     "build/verification/aprl", "build/verification/verification/aprl_hybrid_output_checks",
                     "build/cluster-tools/binding_probe"):
        path = workspace / required
        if not path.is_file() or not os.access(path, os.X_OK):
            raise SystemExit(f"Missing executable: {path}")
    for path in sorted(paths):
        with path.open("rb") as stream:
            if stream.read(4) != b"\x7fELF":
                continue
        dependencies = output("ldd", str(path))
        if "not found" in dependencies:
            raise SystemExit(f"Unresolved dependency: {path}")
        record["executables"][str(path)] = {"sha256": digest(path), "ldd": dependencies}
    for details in record["executables"].values():
        for dependency in re.findall(r"(?:=>\s+|^\s*)(/\S+)", details["ldd"], re.MULTILINE):
            paths.add(Path(dependency))
    source = workspace / "source/aprl"
    tracked = subprocess.check_output(["git", "-C", str(source), "ls-files", "-z", "--",
                                       "verification", "src", "utilities", "postprocess", "CMakeLists.txt"])
    paths.update(source / os.fsdecode(p) for p in tracked.split(b"\0") if p)
    paths.update(p for p in prepared.iterdir() if p.is_file())
    input_path = workspace / "inputs/elegant-particles.h5"
    stat = input_path.stat()
    with h5py.File(input_path, "r") as file:
        record["elegant_input"] = {"link": str(input_path), "target": str(input_path.resolve()),
                                   "bytes": stat.st_size, "mtime_ns": stat.st_mtime_ns,
                                   "records": len(file["particles/records"]),
                                   "note": "metadata inspected read-only; not used by the small SASE harness"}
    if digest(workspace / "build/candidate/aprl") != digest(workspace / "build/verification/aprl"):
        raise SystemExit("Verification APRL differs from candidate: review before freezing.")
    with record_path.open("x") as stream:
        json.dump(record, stream, indent=2)
        stream.write("\n")
    paths.add(record_path)
    with manifest_path.open("x") as stream:
        for path in sorted(paths):
            stream.write(f"{digest(path)}  {path}\n")
    print(f"Recorded {len(paths)} artifacts: {manifest_path}")


if __name__ == "__main__":
    main()
