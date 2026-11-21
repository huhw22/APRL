"""Assemble the versioned 20104 input dataset from its verified native arrays.

Requires Python 3, PyYAML and a C compiler. The installed run helper only needs Python 3.
"""
import argparse
import json
import os
import pathlib
import shutil
import struct
import subprocess
import tarfile
import tempfile

import yaml

from prepare import BEAM_SHA256, HERE, METHODS, components, sha256, verify

READY = {
    "osiris": "2a7a3c54a89eb37a19d8238ca5509503bfed79b734c83c6dc56bdf22f0478a47",
    "left-20056": "c7666d310317732a2ec8f4b1c86a0d79c3839a050f80edb173edebda530f777f",
}
RHO = "9c684e8f850c233b4cb14139acfe2d67c13ac4aa1430d4df23d035a83c01b4e7"
METADATA = "ae894cc88262685b9631e205abfa6cf9af3d9ddba536296e1907245cb423c367"


def checked(path, expected):
    if sha256(path) != expected:
        raise ValueError(f"Archived input does not match the validated 20104 setup: {path}")


def native_manifest(reference, template):
    meta = dict((w[0], w[1:]) for line in (reference / "metadata.txt").read_text().splitlines() if (w := line.split()))
    def doubles(key):
        return [struct.unpack(">d", bytes.fromhex(x))[0] for x in meta[key]]
    mesh = template["mesh"]
    window = mesh["moving_window"]
    boundary = template["boundary"]
    return dict(version=1, reference_rho_sha256=RHO, reference_metadata_sha256=METADATA,
                reference_cells=list(map(int, meta["cells"])), reference_origin_m=doubles("origin_bits"),
                spacing_m=doubles("spacing_bits"), dt_s=doubles("dt_bits")[0],
                reference_beta=doubles("beta_bits")[0], frame_gamma=mesh["boost_gamma"],
                window_cells=mesh["cells"][:2] + [window["cells"]], reference_offset=623,
                cpml_cells=boundary["cells"], target_reflection=boundary["target_reflection"],
                polynomial_order=boundary["polynomial_order"], kappa_max=boundary["kappa_max"],
                alpha_fraction=boundary["alpha_fraction"], beam_sha256=BEAM_SHA256,
                beam_position_offset_m=template["beam"]["input"]["position_offset"],
                beam_electrons=template["beam"]["input"]["electrons"], beam_reference=template["beam"]["reference"],
                transverse_x_binomial3=True, longitudinal_shape=mesh["longitudinal_shape"])


def compress(tree, output, threads):
    compressor = ["pigz", "-n", "-3", "-p", str(threads)] if shutil.which("pigz") else ["gzip", "-n", "-3"]
    with output.open("xb") as file:
        with subprocess.Popen(compressor, stdin=subprocess.PIPE, stdout=file) as proc:
            with tarfile.open(fileobj=proc.stdin, mode="w|", format=tarfile.USTAR_FORMAT) as tar:
                for path in sorted(tree.rglob("*")):
                    if not path.is_file():
                        continue
                    name = str(path.relative_to(tree.parent))
                    info = tarfile.TarInfo(name)
                    info.size = path.stat().st_size
                    info.mode = 0o644
                    with path.open("rb") as src:
                        tar.addfile(info, src)
            proc.stdin.close()
            if proc.wait():
                raise RuntimeError(f"Compression failed: {compressor[0]}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("beam", "osiris", "left-20056", "reference", "output"):
        parser.add_argument("--" + name, type=pathlib.Path, required=True)
    parser.add_argument("--threads", type=int, default=4)
    args = parser.parse_args()
    if args.threads < 1:
        parser.error("--threads must be positive")
    dest = args.output.resolve()
    dest.mkdir(parents=True, exist_ok=True)
    archive = "APRL-2.1.0-beam-cpml-100k.tar.gz"
    artifacts = ["beam-cpml-100k", archive, archive + ".sha256"]
    if any((dest / name).exists() for name in artifacts):
        parser.error("A dataset with this version is already present")
    checked(args.reference / "rho.bin", RHO)
    checked(args.reference / "metadata.txt", METADATA)
    directories = {"osiris": args.osiris, "left-20056": args.left_20056}
    template = yaml.safe_load((HERE / "run-20104.yaml.in").read_text())
    common = native_manifest(args.reference, template)
    with tempfile.TemporaryDirectory(dir=dest, prefix=".input-") as tmp:
        stage = pathlib.Path(tmp)
        tree = stage / artifacts[0]
        tree.mkdir()
        fnv = stage / "fnv1a64"
        subprocess.run([os.environ.get("CC", "cc"), "-std=c99", "-O3", str(HERE / "fnv1a64.c"), "-o", str(fnv)], check=True)
        for name in ("run-20104.yaml.in", "prepare.py", "README.zh-CN.md"):
            shutil.copyfile(HERE / name, tree / name)
        beam = tree / "beam.h5"
        shutil.copyfile(args.beam, beam)
        checked(beam, BEAM_SHA256)
        common["beam_fnv1a64"] = subprocess.check_output([fnv, beam], text=True).strip()
        for method in METHODS:
            source = directories[method]
            target = tree / method
            target.mkdir()
            checked(source / "READY.json", READY[method])
            original = json.loads((source / "READY.json").read_text())
            # The historical exporter recorded an absolute author-workspace path. Keep its digest and physics fields.
            original.pop("archive", None)
            original["archived_ready_sha256"] = READY[method]
            (target / "READY.json").write_text(json.dumps(original, indent=2) + "\n")
            if method == "osiris":
                expected = {k: v["sha256"] for k, v in original["files"].items()}
            else:
                expected = original["sha256"]
            names = components(method)
            if set(expected) != set(names):
                raise ValueError(f"Incomplete native field set for {method}")
            files = [target / (name + ".bin") for name in names]
            for name, path in zip(names, files):
                shutil.copyfile(source / path.name, path)
                checked(path, expected[name])
            fnv_values = subprocess.check_output([fnv, *files], text=True).splitlines()
            if len(fnv_values) != len(names):
                raise RuntimeError("Missing native field fingerprints")
            manifest = dict(common, method=method, arrays_sha256=expected, arrays_fnv1a64=dict(zip(names, fnv_values)))
            (target / "startup.yaml").write_text(yaml.safe_dump(manifest, sort_keys=False))
        (tree / "DATA.json").write_text(json.dumps(dict(version="2.1.0", example="20104",
            methods=list(METHODS), particles=100000, byte_order="little-endian IEEE754 binary64",
            archived_ready_sha256=READY, reference_metadata_sha256=METADATA, reference_rho_sha256=RHO), indent=2) + "\n")
        entries = sorted(p for p in tree.rglob("*") if p.is_file())
        (tree / "SHA256SUMS").write_text("".join(f"{sha256(p)}  {p.relative_to(tree)}\n" for p in entries))
        verify(tree, METHODS)
        compress(tree, stage / archive, args.threads)
        archive_digest = sha256(stage / archive)
        (stage / artifacts[2]).write_text(f"{archive_digest}  {archive}\n")
        for name in artifacts:
            (stage / name).rename(dest / name)
        print(f"{archive_digest}  {archive}")


if __name__ == "__main__":
    main()
