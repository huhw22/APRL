#!/usr/bin/env python3
"""Validate the immutable inputs recorded before an APRL case was launched."""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def reject(message: str) -> int:
    print(f"case provenance mismatch: {message}", file=sys.stderr)
    return 1


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("run_dir", type=Path)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--plan", type=Path, required=True)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--fingerprint", required=True)
    parser.add_argument("--input-sha256", required=True)
    parser.add_argument("--generator-sha256", required=True)
    args = parser.parse_args()

    required = {
        "RUN_FINGERPRINT": args.run_dir / "RUN_FINGERPRINT",
        "submitted-config.yaml": args.run_dir / "submitted-config.yaml",
        "submitted-plan.json": args.run_dir / "submitted-plan.json",
        "input-noise.json": args.run_dir / "input-noise.json",
        "run-manifest.yaml": args.run_dir / "run-manifest.yaml",
    }
    for label, path in required.items():
        if not path.is_file() or path.stat().st_size == 0:
            return reject(f"missing or empty {label}: {path}")
    if not args.input.is_file() or args.input.stat().st_size == 0:
        return reject(f"missing or empty current particle input: {args.input}")
    if sha256(args.input) != args.input_sha256:
        return reject("current particle input does not match --input-sha256")

    saved_fingerprint = required["RUN_FINGERPRINT"].read_text(
        encoding="utf-8"
    ).strip()
    if saved_fingerprint != args.fingerprint:
        return reject(
            "RUN_FINGERPRINT does not match the freshly computed case fingerprint"
        )

    if sha256(required["submitted-config.yaml"]) != sha256(args.config):
        return reject("submitted-config.yaml differs from the current rendered card")
    if sha256(required["submitted-plan.json"]) != sha256(args.plan):
        return reject("submitted-plan.json differs from the current campaign plan")

    try:
        metadata = json.loads(
            required["input-noise.json"].read_text(encoding="utf-8")
        )
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        return reject(f"cannot read input-noise.json: {error}")

    expected = {
        "schema": "aprl-sase-input-metadata-v1",
        "input_file_name": args.input.name,
        "input_sha256": args.input_sha256,
        "generator_sha256": args.generator_sha256,
    }
    for key, value in expected.items():
        if metadata.get(key) != value:
            return reject(
                f"input-noise.json field {key!r} is {metadata.get(key)!r}, "
                f"expected {value!r}"
            )
    if not isinstance(metadata.get("attributes"), dict):
        return reject("input-noise.json has no attributes object")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
