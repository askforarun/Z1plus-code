# Z1plus-code

This repository includes the Z1+ workflow for PVA (polyvinyl alcohol) systems:

1. install Z1+
2. convert a GROMACS `.gro` file to `config.Z1`
3. run Z1+ on that file

The repository includes a PVA example (polymer):

- input `.gro` file: `examples/pva-n51/isotropic_equilibration_gromacs_N51.gro`
- converted Z1 input: `examples/pva-n51/pva-backbone-N51-config.Z1`
- backbone mapping file: `examples/pva-n51/pva-backbone-N51-info.txt`
- example Z1+ outputs: `examples/pva-n51/outputs/`

For this bundled example:

- the original `.gro` file contains `300` PVA chains in a periodic box, along with additional `GLU` residues
- the original atomistic configuration contains `55650` atoms in total
- each chain contributes `51` backbone carbon atoms to the Z1 input
- the converted Z1 file therefore contains `15300` backbone coordinates in total

In this workflow, the Z1 file contains only the polymer backbone used for Z1+ analysis. It does not keep the full atomistic coordinates of hydrogens, side-group atoms, or the additional `GLU` residues. The coordinates are written chain-by-chain, and periodic crossings are unwrapped along each chain before the Z1 file is written.

The Z1 format used here is:

- line 1: number of chains
- line 2: box lengths `boxx boxy boxz`
- line 3: chain lengths for all chains
- remaining lines: `x y z` coordinates for each chain, written in order

The related Z1+ paper is available here:
https://www.sciencedirect.com/science/article/pii/S0010465522002867

## Install Z1+

Z1+ needs `perl` and a Fortran compiler such as `gfortran` or `ifort`.

From the repository root:

```bash
mkdir -p downloads/z1plus
tar -xzf downloads/z1plus.tar.gz -C downloads/z1plus
cd downloads/z1plus
perl Z1+install.pl
cd ../..
```

After that, run the generated launcher from outside the installation directory:

```bash
/absolute/path/to/downloads/z1plus/Z1+ config.Z1
```

## Convert A `.gro` File To Z1 Format

The example below is for PVA data where:

- the polymer residue name is `PVA`
- backbone atoms are named `C1`, `C2`, `C3`, ...
- coordinates in the `.gro` file are in nm and are scaled by `10.0`

Use this Python code:

```python
#!/usr/bin/env python3

from __future__ import annotations

import re
import sys
from pathlib import Path


CARBON_NAME = re.compile(r"^C(\d+)$")


def round_away_from_zero(value: float) -> int:
    if value == 0:
        return 0
    return int(value / abs(value) * int(abs(value) + 0.5))


def parse_gro(path: Path, scale: float) -> tuple[str, list[dict[str, object]], tuple[float, float, float]]:
    lines = path.read_text().splitlines()
    if len(lines) < 3:
        raise ValueError(f"{path} is too short to be a valid .gro file")

    title = lines[0].strip()
    atom_count = int(lines[1].strip())
    atom_lines = lines[2:2 + atom_count]
    if len(atom_lines) != atom_count:
        raise ValueError(f"{path} declares {atom_count} atoms but contains {len(atom_lines)} atom lines")

    atoms: list[dict[str, object]] = []
    for line in atom_lines:
        atoms.append(
            {
                "residue_id": int(line[0:5]),
                "residue_name": line[5:10].strip(),
                "atom_name": line[10:15].strip(),
                "atom_id": int(line[15:20]),
                "x": float(line[20:28]) * scale,
                "y": float(line[28:36]) * scale,
                "z": float(line[36:44]) * scale,
            }
        )

    box_tokens = lines[2 + atom_count].split()
    box = tuple(float(token) * scale for token in box_tokens[:3])
    return title, atoms, box


def group_residue_blocks(atoms: list[dict[str, object]]) -> list[list[dict[str, object]]]:
    if not atoms:
        return []
    blocks = []
    current = [atoms[0]]
    current_key = (atoms[0]["residue_id"], atoms[0]["residue_name"])
    for atom in atoms[1:]:
        key = (atom["residue_id"], atom["residue_name"])
        if key == current_key:
            current.append(atom)
        else:
            blocks.append(current)
            current = [atom]
            current_key = key
    blocks.append(current)
    return blocks


def backbone_atoms_from_pva_block(block: list[dict[str, object]], residue_name: str) -> list[dict[str, object]]:
    selected = []
    for atom in block:
        match = CARBON_NAME.match(str(atom["atom_name"]))
        if match is not None:
            selected.append((int(match.group(1)), atom))
    if not selected:
        raise ValueError(f"Residue {block[0]['residue_id']} {residue_name} has no backbone carbon atoms")
    selected.sort(key=lambda item: item[0])
    return [atom for _, atom in selected]


def unwrap_chain(chain: list[dict[str, object]], box: tuple[float, float, float]) -> list[tuple[float, float, float]]:
    boxx, boxy, boxz = box
    coords = [(float(chain[0]["x"]), float(chain[0]["y"]), float(chain[0]["z"]))]
    prev_raw = chain[0]
    prev_unwrapped = coords[0]
    for atom in chain[1:]:
        dx = float(atom["x"]) - float(prev_raw["x"])
        dy = float(atom["y"]) - float(prev_raw["y"])
        dz = float(atom["z"]) - float(prev_raw["z"])
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
    gro_path: str,
    output_path: str = "config.Z1",
    info_path: str = "backbone-info.txt",
    polymer_residue: str = "PVA",
    scale: float = 10.0,
) -> None:
    title, atoms, box = parse_gro(Path(gro_path), scale)
    blocks = group_residue_blocks(atoms)

    chains = []
    for block in blocks:
        if block[0]["residue_name"] == polymer_residue:
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

    for chain_id, chain in enumerate(chains, start=1):
        for x, y, z in unwrap_chain(chain, box):
            config_lines.append(f"{x:.15g} {y:.15g} {z:.15g}")
        info_lines.append(f"{chain_id} {len(chain)}")
        info_lines.extend(str(atom["atom_id"]) for atom in chain)

    Path(output_path).write_text("\n".join(config_lines) + "\n")
    Path(info_path).write_text("\n".join(info_lines) + "\n")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        raise SystemExit("usage: python gro_to_z1.py input.gro [output.Z1] [info.txt]")
    gro_file = sys.argv[1]
    output_file = sys.argv[2] if len(sys.argv) > 2 else "config.Z1"
    info_file = sys.argv[3] if len(sys.argv) > 3 else "backbone-info.txt"
    convert_gro_to_z1(gro_file, output_file, info_file)
    print(f"created {output_file}")
    print(f"created {info_file}")
```

Example:

```bash
python gro_to_z1.py examples/pva-n51/isotropic_equilibration_gromacs_N51.gro config.Z1 backbone-info.txt
```

This creates:

- `config.Z1` for Z1+
- `backbone-info.txt` with the mapping from Z1 backbone atoms back to the original `.gro` atom ids

## Run Z1+

Once `config.Z1` has been created, run:

```bash
/absolute/path/to/downloads/z1plus/Z1+ config.Z1
```

Z1+ will then generate its usual output files such as `Z1+summary.dat`, `Z1+SP.dat`, `Z1+initconfig.dat`, and related analysis files.
You will find these files in the working directory where you ran the `Z1+ config.Z1` command.

This repository also includes one ready-made output set in `examples/pva-n51/outputs/` so users can compare their own run against a known PVA example.

## Citation

```text
M. Kröger, J. D. Dietz, R. S. Hoy and C. Luap,
The Z1+ package: Shortest multiple disconnected path for the analysis of entanglements in macromolecular systems,
Comput. Phys. Commun. 283 (2023) 108567
https://doi.org/10.1016/j.cpc.2022.108567
```
