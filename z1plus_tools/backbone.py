from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from collections import deque

from .lammps import choose_coordinate_fields, iter_lammps_dump_frames, parse_lammps_data


@dataclass(slots=True)
class BackboneChain:
    mol_id: int
    original_atom_ids: list[int]


@dataclass(slots=True)
class BackboneSummary:
    molecules: int
    linear_atoms: int
    config_path: Path
    info_path: Path
    frames_written: int


class BackboneError(RuntimeError):
    pass


def _parse_dump_box(frame) -> tuple[float, float, float, float]:
    tokens = frame.box_header.split()
    numeric_lines = [tuple(float(value) for value in line.split()) for line in frame.box_lines]
    xlo = numeric_lines[0][0]
    xhi = numeric_lines[0][1]
    ylo = numeric_lines[1][0]
    yhi = numeric_lines[1][1]
    zlo = numeric_lines[2][0]
    zhi = numeric_lines[2][1]
    xy = xz = yz = 0.0
    dcol = (0, 0, 0)
    if 'xy' in tokens:
        xlobound, xhibound, xy = numeric_lines[0]
        ylobound, yhibound, xz = numeric_lines[1]
        zlobound, zhibound, yz = numeric_lines[2]
        xlo, xhi = xlobound, xhibound
        ylo, yhi = ylobound, yhibound
        zlo, zhi = zlobound, zhibound
        if xy:
            epc = xy
            if xy < 0.0:
                xlo -= xy
            if xy > 0.0:
                xhi -= xy
            boxx, boxy, boxz = xhi - xlo, yhi - ylo, zhi - zlo
            dcol = (0, 0, 0)
        elif xz:
            epc = xz
            if xz < 0.0:
                xlo -= xz
            if xz > 0.0:
                xhi -= xz
            boxx, boxy, boxz = xhi - xlo, zhi - zlo, yhi - ylo
            dcol = (0, 1, -1)
        elif yz:
            epc = yz
            if yz < 0.0:
                ylo -= yz
            if yz > 0.0:
                yhi -= yz
            boxx, boxy, boxz = yhi - ylo, zhi - zlo, xhi - xlo
            dcol = (1, 1, -2)
        else:
            epc = 0.0
            boxx, boxy, boxz = xhi - xlo, yhi - ylo, zhi - zlo
    else:
        boxx, boxy, boxz = xhi - xlo, yhi - ylo, zhi - zlo
        epc = 0.0
    return boxx, boxy, boxz, epc, dcol


def _bfs_distances(start: int, adjacency: dict[int, list[int]], members: set[int]) -> dict[int, int]:
    distances = {member: -1 for member in members}
    distances[start] = 0
    queue = deque([start])
    while queue:
        current = queue.popleft()
        for neighbor in adjacency[current]:
            if neighbor in members and distances[neighbor] == -1:
                distances[neighbor] = distances[current] + 1
                queue.append(neighbor)
    return distances


def _longest_path_for_molecule(members: list[int], adjacency: dict[int, list[int]], anchor: int) -> list[int]:
    member_set = set(members)
    distances = _bfs_distances(anchor, adjacency, member_set)
    chainbeg = max(members, key=lambda atom_id: distances[atom_id])
    distances = _bfs_distances(chainbeg, adjacency, member_set)
    chainend = max(members, key=lambda atom_id: distances[atom_id])
    dist = distances[chainend]
    path = [chainend]
    seen = {chainend}
    current = chainend
    while dist > 0:
        dist -= 1
        next_atom = None
        for neighbor in adjacency[current]:
            if neighbor in member_set and distances.get(neighbor) == dist and neighbor not in seen:
                next_atom = neighbor
                break
        if next_atom is None:
            raise BackboneError(f'Failed to reconstruct the longest path ending at atom {chainend}')
        path.append(next_atom)
        seen.add(next_atom)
        current = next_atom
    return path


def extract_backbone(
    *,
    data_path: str | Path,
    dump_path: str | Path | None = None,
    ignore_types: tuple[int, ...] = (),
    config_path: str | Path = 'config.Z1',
    info_path: str | Path = 'backbone-info.txt',
) -> BackboneSummary:
    data = parse_lammps_data(data_path)
    ignore_type_set = set(ignore_types)
    kept_atoms = [atom for atom in data.atoms if atom.atom_type not in ignore_type_set]
    row_by_original_id: dict[int, int] = {}
    original_id_by_row: dict[int, int] = {}
    mol_by_row: dict[int, int] = {}
    type_by_row: dict[int, int] = {}
    x_by_row: dict[int, float] = {}
    y_by_row: dict[int, float] = {}
    z_by_row: dict[int, float] = {}
    anchor_by_mol: dict[int, int] = {}
    for row, atom in enumerate(kept_atoms, start=1):
        row_by_original_id[atom.atom_id] = row
        original_id_by_row[row] = atom.atom_id
        mol_by_row[row] = atom.mol_id
        type_by_row[row] = atom.atom_type
        x_by_row[row] = atom.x
        y_by_row[row] = atom.y
        z_by_row[row] = atom.z
        anchor_by_mol[atom.mol_id] = row

    adjacency = {row: [] for row in range(1, len(kept_atoms) + 1)}
    for bond in data.bonds:
        row1 = row_by_original_id.get(bond.atom1)
        row2 = row_by_original_id.get(bond.atom2)
        if row1 is None or row2 is None:
            continue
        adjacency[row1].append(row2)
        adjacency[row2].append(row1)

    functionalities = [len(adjacency[row]) for row in adjacency]
    if functionalities:
        print(f'functionalities range between {min(functionalities)} and {max(functionalities)}')

    members_by_mol: dict[int, list[int]] = {}
    for row, mol_id in mol_by_row.items():
        members_by_mol.setdefault(mol_id, []).append(row)

    backbone_chains: list[BackboneChain] = []
    linear_atom_id_by_row: dict[int, int] = {}
    linear_mol_by_row: dict[int, int] = {}
    xyz_lines: list[str] = []
    next_linear_id = 1
    info_lines: list[str] = []
    chain_lengths: list[int] = []
    for mol_id in sorted(members_by_mol):
        members = sorted(members_by_mol[mol_id])
        print(f'mol {mol_id} has {len(members)} members')
        path = _longest_path_for_molecule(members, adjacency, anchor_by_mol[mol_id])
        print(
            f'longest path with {len(path) - 1} bonds from id {original_id_by_row[path[-1]]} '
            f'[type {type_by_row[path[-1]]}] to {original_id_by_row[path[0]]} [type {type_by_row[path[0]]}]'
        )
        original_ids = [original_id_by_row[row] for row in path]
        backbone_chains.append(BackboneChain(mol_id=mol_id, original_atom_ids=original_ids))
        chain_lengths.append(len(path))
        info_lines.append(f'{mol_id} {len(path)}')
        for row in path:
            xyz_lines.append(f'{x_by_row[row]} {y_by_row[row]} {z_by_row[row]}')
            linear_mol_by_row[row] = mol_id
            linear_atom_id_by_row[row] = next_linear_id
            next_linear_id += 1
            info_lines.append(str(original_id_by_row[row]))

    config_path = Path(config_path)
    box = data.box
    config_lines = [
        str(len(backbone_chains)),
        f'{box.xhi - box.xlo} {box.yhi - box.ylo} {box.zhi - box.zlo}',
        ' '.join(str(length) for length in chain_lengths),
        *xyz_lines,
    ]
    config_path.write_text('\n'.join(config_lines) + '\n')
    print('created config.Z1')

    info_path = Path(info_path)
    info_path.write_text('\n'.join(info_lines) + ('\n' if info_lines else ''))
    print('created backbone-info.txt')

    frames_written = 1
    if dump_path is not None:
        trajectory_lines: list[str] = []
        for raw_frame in iter_lammps_dump_frames(dump_path):
            if raw_frame.atom_count != len(data.atoms):
                raise BackboneError(
                    f'conflicting data and dump files [{len(data.atoms)} versus {raw_frame.atom_count} atoms]'
                )
            atom_fields = raw_frame.atom_fields
            coord_fields = choose_coordinate_fields(atom_fields)
            boxx, boxy, boxz, epc, dcol = _parse_dump_box(raw_frame)
            per_linear_id: dict[int, str] = {}
            for atom_values in raw_frame.atoms:
                original_id = int(atom_values['id'])
                row = row_by_original_id.get(original_id)
                if row is None or row not in linear_atom_id_by_row:
                    continue
                coords = [float(atom_values[coord_fields[0]]), float(atom_values[coord_fields[1]]), float(atom_values[coord_fields[2]])]
                x = coords[0 + dcol[0]]
                y = coords[1 + dcol[1]]
                z = coords[2 + dcol[2]]
                per_linear_id[linear_atom_id_by_row[row]] = f'{x} {y} {z}'
            timestep = int(raw_frame.timestep)
            print(f'[{dump_path}] processing time step {timestep}')
            trajectory_lines.extend([
                str(len(backbone_chains)),
                f'{boxx} {boxy} {boxz}',
                ' '.join(str(length) for length in chain_lengths),
            ])
            for linear_id in range(1, next_linear_id):
                trajectory_lines.append(per_linear_id[linear_id])
            trajectory_lines.extend(['-1', str(epc)])
        config_path.write_text('\n'.join(trajectory_lines) + ('\n' if trajectory_lines else ''))
        print('created config.Z1 (trajectory file)')
        frames_written = len(list(iter_lammps_dump_frames(dump_path)))

    return BackboneSummary(
        molecules=len(backbone_chains),
        linear_atoms=next_linear_id - 1,
        config_path=config_path,
        info_path=info_path,
        frames_written=frames_written,
    )
