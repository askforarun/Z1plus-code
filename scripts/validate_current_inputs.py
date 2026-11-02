#!/usr/bin/env python3

"""Check the bundled PVA--GLU Z1 inputs and chain-to-atom mappings."""

from __future__ import annotations

import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1] / "examples" / "pva-glu-current"
STRANDS = {17: 2700, 21: 2100, 25: 1800}
TEMPERATURES = (300, 625)
REPLICAS = (1, 2, 3)


def check_one(n: int, replica: int, temperature: int) -> None:
    directory = ROOT / f"n{n}" / f"replica-{replica}" / f"{temperature}K"
    config = (directory / "config.Z1").read_text().splitlines()
    if len(config) < 3:
        raise ValueError(f"{directory}: incomplete config.Z1")

    chains = int(config[0])
    box = [float(value) for value in config[1].split()[:3]]
    lengths = [int(value) for value in config[2].split()]
    expected_length = 2 * n + 1
    expected_chains = STRANDS[n]
    if chains != expected_chains or len(lengths) != chains:
        raise ValueError(f"{directory}: unexpected chain count")
    if any(length != expected_length for length in lengths):
        raise ValueError(f"{directory}: expected {expected_length} carbons per strand")
    if len(box) != 3 or any(not math.isfinite(value) or value <= 0 for value in box):
        raise ValueError(f"{directory}: invalid orthorhombic box")

    coordinate_lines = config[3:]
    if len(coordinate_lines) != chains * expected_length:
        raise ValueError(f"{directory}: coordinate count does not match chain lengths")
    for line in coordinate_lines:
        xyz = [float(value) for value in line.split()[:3]]
        if len(xyz) != 3 or not all(math.isfinite(value) for value in xyz):
            raise ValueError(f"{directory}: invalid coordinate line")

    mapping = (directory / "backbone-info.txt").read_text().splitlines()
    records = [line.strip() for line in mapping if line.strip() and not line.startswith("#")]
    offset = 0
    for chain_id in range(1, chains + 1):
        if offset >= len(records):
            raise ValueError(f"{directory}: mapping ends before chain {chain_id}")
        header = records[offset].split()
        offset += 1
        if len(header) != 2 or tuple(map(int, header)) != (chain_id, expected_length):
            raise ValueError(f"{directory}: invalid mapping header for chain {chain_id}")
        for _ in range(expected_length):
            if offset >= len(records) or len(records[offset].split()) != 1:
                raise ValueError(f"{directory}: incomplete atom mapping for chain {chain_id}")
            int(records[offset])
            offset += 1
    if offset != len(records):
        raise ValueError(f"{directory}: unexpected trailing mapping records")


def check_reference_outputs() -> None:
    directory = ROOT / "n17" / "replica-1" / "300K" / "outputs"
    required = (
        "Lpp_values.dat",
        "N_values.dat",
        "Ree_values.dat",
        "Z1+SP.dat",
        "Z1+initconfig.dat",
        "Z1+summary.dat",
        "Z1+summary.html",
        "Z_values.dat",
    )
    for name in required:
        path = directory / name
        if not path.is_file() or path.stat().st_size == 0:
            raise ValueError(f"Missing or empty reference output: {path}")

    summary = [float(value) for value in (directory / "Z1+summary.dat").read_text().split()]
    if len(summary) < 3 or summary[0] != 1 or summary[1] != STRANDS[17]:
        raise ValueError("Reference Z1+ summary does not match the n=17 example")
    if not math.isclose(summary[2], 2 * 17 + 1):
        raise ValueError("Reference Z1+ summary has an unexpected strand length")


def main() -> None:
    count = 0
    for n in STRANDS:
        for replica in REPLICAS:
            for temperature in TEMPERATURES:
                check_one(n, replica, temperature)
                count += 1
    check_reference_outputs()
    print(f"Validated {count} PVA--GLU Z1 inputs, mappings, and reference outputs.")


if __name__ == "__main__":
    main()
