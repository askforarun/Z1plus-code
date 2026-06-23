#!/usr/bin/env python3

from __future__ import annotations

import argparse
import math
import re
from dataclasses import dataclass
from pathlib import Path


CARBON_NAME = re.compile(r"^C(\d+)$")


@dataclass(slots=True)
class GroAtom:
    residue_id: int
    residue_name: str
    atom_name: str
    atom_id: int
    x: float
    y: float
    z: float


def round_away_from_zero(value: float) -> int:
    if value == 0:
        return 0
    return int(math.copysign(math.floor(abs(value) + 0.5), value))


def parse_gro(path: str | Path, *, scale: float) -> tuple[str, list[GroAtom], tuple[float, float, float]]:
    path = Path(path)
    lines = path.read_text().splitlines()
    if len(lines) < 3:
        raise ValueError(f"{path} is too short to be a valid .gro file")

    title = lines[0].strip()
    atom_count = int(lines[1].strip())
    atom_lines = lines[2:2 + atom_count]
    if len(atom_lines) != atom_count:
        raise ValueError(f"{path} declares {atom_count} atoms but contains {len(atom_lines)} atom lines")
    if len(lines) < atom_count + 3:
        raise ValueError(f"{path} is missing its box line")

    atoms: list[GroAtom] = []
    for line_number, line in enumerate(atom_lines, start=3):
        if len(line) < 44:
            raise ValueError(f"{path}:{line_number} is too short to be a valid .gro atom record")
        try:
            residue_id = int(line[0:5])
            residue_name = line[5:10].strip()
            atom_name = line[10:15].strip()
            atom_id = int(line[15:20])
            x = float(line[20:28]) * scale
            y = float(line[28:36]) * scale
            z = float(line[36:44]) * scale
        except ValueError as exc:
            raise ValueError(f"{path}:{line_number} contains an invalid .gro atom record") from exc
        atoms.append(
            GroAtom(
                residue_id=residue_id,
                residue_name=residue_name,
                atom_name=atom_name,
                atom_id=atom_id,
                x=x,
                y=y,
                z=z,
            )
        )

    box_tokens = lines[2 + atom_count].split()
    if len(box_tokens) < 3:
        raise ValueError(f"{path} box line must contain at least three values")
    box = tuple(float(token) * scale for token in box_tokens[:3])
    return title, atoms, box


def group_residue_blocks(atoms: list[GroAtom]) -> list[list[GroAtom]]:
    if not atoms:
        return []
    blocks: list[list[GroAtom]] = []
    current = [atoms[0]]
    current_key = (atoms[0].residue_id, atoms[0].residue_name)
    for atom in atoms[1:]:
        key = (atom.residue_id, atom.residue_name)
        if key == current_key:
            current.append(atom)
            continue
        blocks.append(current)
        current = [atom]
        current_key = key
    blocks.append(current)
    return blocks


def backbone_atoms_from_pva_block(block: list[GroAtom], residue_name: str) -> list[GroAtom]:
    residue_ids = {atom.residue_id for atom in block}
    residue_names = {atom.residue_name for atom in block}
    if len(residue_ids) != 1 or residue_names != {residue_name}:
        raise ValueError("A residue block does not match the expected polymer residue")

    selected: list[tuple[int, GroAtom]] = []
    for atom in block:
        match = CARBON_NAME.match(atom.atom_name)
        if match is None:
            continue
        selected.append((int(match.group(1)), atom))
    if not selected:
        raise ValueError(f"Residue {block[0].residue_id} {residue_name} does not contain any backbone carbon atoms")

    selected.sort(key=lambda item: item[0])
    expected = list(range(1, len(selected) + 1))
    actual = [index for index, _atom in selected]
    if actual != expected:
        raise ValueError(
            f"Residue {block[0].residue_id} {residue_name} has non-consecutive carbon labels: {actual[:10]}..."
        )
    return [atom for _index, atom in selected]


def unwrap_chain(chain: list[GroAtom], box: tuple[float, float, float]) -> list[tuple[float, float, float]]:
    if not chain:
        return []
    boxx, boxy, boxz = box
    coords = [(chain[0].x, chain[0].y, chain[0].z)]
    prev_raw = chain[0]
    prev_unwrapped = coords[0]
    for atom in chain[1:]:
        dx = atom.x - prev_raw.x
        dy = atom.y - prev_raw.y
        dz = atom.z - prev_raw.z
        dx -= boxx * round_away_from_zero(dx / boxx)
        dy -= boxy * round_away_from_zero(dy / boxy)
        dz -= boxz * round_away_from_zero(dz / boxz)
        current = (
            prev_unwrapped[0] + dx,
            prev_unwrapped[1] + dy,
            prev_unwrapped[2] + dz,
        )
        coords.append(current)
        prev_raw = atom
        prev_unwrapped = current
    return coords


def convert_gro_to_z1(
    *,
    gro_path: str | Path,
    output_path: str | Path,
    info_path: str | Path,
    polymer_residue: str,
    scale: float,
) -> tuple[int, int]:
    title, atoms, box = parse_gro(gro_path, scale=scale)
    blocks = group_residue_blocks(atoms)

    chains: list[list[GroAtom]] = []
    for block in blocks:
        if block[0].residue_name != polymer_residue:
            continue
        chains.append(backbone_atoms_from_pva_block(block, polymer_residue))

    if not chains:
        raise ValueError(f"No residue blocks named {polymer_residue!r} were found in {gro_path}")

    config_lines = [str(len(chains)), f"{box[0]:.15g} {box[1]:.15g} {box[2]:.15g}"]
    config_lines.append(" ".join(str(len(chain)) for chain in chains))

    info_lines = [
        f"# source {title}",
        f"# gro {Path(gro_path).resolve()}",
        f"# polymer_residue {polymer_residue}",
        "# each block below is: chain_id chain_length followed by original .gro atom ids",
    ]

    total_atoms = 0
    for chain_id, chain in enumerate(chains, start=1):
        total_atoms += len(chain)
        unwrapped = unwrap_chain(chain, box)
        for x, y, z in unwrapped:
            config_lines.append(f"{x:.15g} {y:.15g} {z:.15g}")
        info_lines.append(f"{chain_id} {len(chain)}")
        info_lines.extend(str(atom.atom_id) for atom in chain)

    output_path = Path(output_path)
    output_path.write_text("\n".join(config_lines) + "\n")

    info_path = Path(info_path)
    info_path.write_text("\n".join(info_lines) + "\n")

    return len(chains), total_atoms


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Convert a .gro file directly to Z1 format using only polymer backbone carbons from selected residue blocks."
    )
    parser.add_argument("gro_file", help="Input GROMACS .gro file.")
    parser.add_argument(
        "-o",
        "--output",
        default="config.Z1",
        help="Output Z1 file path. Default: config.Z1",
    )
    parser.add_argument(
        "--info",
        default="backbone-info.txt",
        help="Output mapping/info file path. Default: backbone-info.txt",
    )
    parser.add_argument(
        "--polymer-residue",
        default="PVA",
        help="Residue name to treat as polymer chains. Default: PVA",
    )
    parser.add_argument(
        "--scale",
        type=float,
        default=10.0,
        help="Coordinate scaling factor applied to .gro coordinates and box lengths. Default: 10.0",
    )
    args = parser.parse_args(argv)

    chains, backbone_atoms = convert_gro_to_z1(
        gro_path=args.gro_file,
        output_path=args.output,
        info_path=args.info,
        polymer_residue=args.polymer_residue,
        scale=args.scale,
    )
    print(f"created {args.output}")
    print(f"created {args.info}")
    print(f"selected {backbone_atoms} backbone atoms across {chains} chain(s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
