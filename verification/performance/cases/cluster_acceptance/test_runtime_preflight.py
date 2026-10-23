#!/usr/bin/env python3
"""Exercise runtime-policy decisions without MPI launches or scheduler calls."""
import copy
import importlib.util
import json
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

spec = importlib.util.spec_from_file_location("runtime_check", Path(__file__).with_name("check-runtime.py"))
runtime = importlib.util.module_from_spec(spec)
spec.loader.exec_module(runtime)


class RuntimeChecks(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.prefix = self.root / "dependencies"
        self.paths = {name: str(self.root / name) for name in ("libmpi.so.40", "libhdf5.so.310", "libstdc++.so.6")}
        self.hashes = {path: "prepared-" + name for name, path in self.paths.items()}
        self.executables = [str(self.root / path) for path in (
            "build/baseline/aprl", "build/candidate/aprl", "build/candidate/sase_bunch_generator")]
        self.ldd = "\n".join(f"{name} => {path} (0x123)" for name, path in self.paths.items())
        self.outputs = {executable: self.ldd for executable in self.executables}
        self.tools = {name: "/usr/bin/" + name for name in ("srun", "scontrol", "mpirun")}
        self.versions = {"srun": "slurm 23.02.6", "scontrol": "slurm 23.02.6", "mpirun": "mpirun (Open MPI) 4.0.3"}
        self.record = {
            "workspace": str(self.root), "dependency_prefix": str(self.prefix), "python": sys.executable,
            "tools": self.tools, "tool_versions": self.versions,
            "strict_files": {self.paths[name]: self.hashes[self.paths[name]] for name in ("libmpi.so.40", "libhdf5.so.310")},
            "node_local_files": {self.paths["libstdc++.so.6"]: {"sha256": self.hashes[self.paths["libstdc++.so.6"]], "realpath": self.paths["libstdc++.so.6"]}},
            "executables": {path: {"sha256": "unused", "ldd": self.ldd} for path in self.executables},
        }

    def run_command(self, arguments):
        if arguments[:2] == ["ldd", "-r"]:
            return 0, self.outputs[arguments[2]]
        name = Path(arguments[0]).name
        return 0, "[admin:00001] opal_ifinit: socket() failed with errno=1\n" + self.versions[name]

    def check(self):
        with patch.dict(os.environ, APRL_DEPENDENCY_PREFIX=str(self.prefix)), \
             patch.object(runtime.shutil, "which", side_effect=self.tools.get):
            return runtime.check(self.record, run=self.run_command, hash_file=self.hashes.__getitem__)

    def test_node_local_hash_drift_is_audited_and_allowed(self):
        self.hashes[self.paths["libstdc++.so.6"]] = "node-local-updated-library"
        result = self.check()
        self.assertTrue(result["pass"], result["errors"])
        self.assertTrue(result["node_local_files"][self.paths["libstdc++.so.6"]]["changed"])

    def test_mpi_hash_drift_is_rejected(self):
        self.hashes[self.paths["libmpi.so.40"]] = "other-mpi-build"
        result = self.check()
        self.assertFalse(result["pass"])
        self.assertTrue(any("Pinned dependency changed" in error for error in result["errors"]))

    def test_private_hdf5_substitution_is_rejected(self):
        other = str(self.root / "other/libhdf5.so.310")
        self.hashes[other] = self.hashes[self.paths["libhdf5.so.310"]]
        self.outputs[self.executables[1]] = self.ldd.replace(self.paths["libhdf5.so.310"], other)
        result = self.check()
        self.assertFalse(result["pass"])
        self.assertTrue(any("Pinned dependency changed" in error for error in result["errors"]))

    def test_unresolved_relocation_is_rejected_even_when_ldd_returns_zero(self):
        self.outputs[self.executables[0]] += "\nundefined symbol: missing_function\t(/some/library.so)"
        result = self.check()
        self.assertFalse(result["pass"])
        self.assertTrue(any("Unresolved runtime dependencies" in error for error in result["errors"]))

    def test_system_library_alias_path_is_allowed(self):
        alias = str(self.root / "other-system/libstdc++.so.6")
        self.hashes[alias] = "compatible-other-node"
        self.outputs = {path: text.replace(self.paths["libstdc++.so.6"], alias) for path, text in self.outputs.items()}
        result = self.check()
        self.assertTrue(result["pass"], result["errors"])

    def test_slurm_version_mismatch_is_rejected(self):
        self.record["tool_versions"] = copy.deepcopy(self.versions)
        self.versions["srun"] = "slurm 24.05.1"
        result = self.check()
        self.assertFalse(result["pass"])
        self.assertTrue(any("Incompatible tool version" in error for error in result["errors"]))

    def test_missing_library_is_rejected(self):
        self.outputs[self.executables[0]] = self.ldd.replace(self.paths["libmpi.so.40"], "not found")
        result = self.check()
        self.assertFalse(result["pass"])
        self.assertTrue(any("Unresolved runtime dependencies" in error for error in result["errors"]))

    def test_cli_preserves_a_failure_audit(self):
        source = self.root / "record.json"
        destination = self.root / "audit.json"
        source.write_text(json.dumps(self.record))
        with patch.object(sys, "argv", ["check-runtime.py", str(source), str(destination)]), \
             patch.object(runtime, "check", side_effect=ValueError("malformed runtime record")), \
             patch("builtins.print"):
            self.assertEqual(runtime.main(), 1)
        result = json.loads(destination.read_text())
        self.assertFalse(result["pass"])
        self.assertIn("malformed runtime record", result["errors"][0])


if __name__ == "__main__":
    unittest.main()
