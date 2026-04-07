from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from .lammps import choose_coordinate_fields, iter_lammps_dump_frames, parse_lammps_data


@dataclass(slots=True)
class CorrectedVmdSummary:
    atoms: int
    bonds: int
    max_atom_id: int
    chains: int
    corrected_data_path: Path
    corrected_dump_path: Path | None = None
    frames: int = 0


class _UnionFind:
    def __init__(self, members: list[int]) -> None:
        self.parent = {member: member for member in members}
        self.rank = {member: 0 for member in members}

    def find(self, member: int) -> int:
        parent = self.parent[member]
        if parent != member:
            self.parent[member] = self.find(parent)
        return self.parent[member]

    def union(self, left: int, right: int) -> None:
        root_left = self.find(left)
        root_right = self.find(right)
        if root_left == root_right:
            return
        if self.rank[root_left] < self.rank[root_right]:
            root_left, root_right = root_right, root_left
        self.parent[root_right] = root_left
        if self.rank[root_left] == self.rank[root_right]:
            self.rank[root_left] += 1


def assign_molecule_ids(atom_ids: list[int], bonds: list[tuple[int, int]]) -> dict[int, int]:
    union_find = _UnionFind(atom_ids)
    for atom1, atom2 in bonds:
        union_find.union(atom1, atom2)

    components: dict[int, list[int]] = {}
    for atom_id in atom_ids:
        root = union_find.find(atom_id)
        components.setdefault(root, []).append(atom_id)

    ordered_components = sorted((sorted(component) for component in components.values()), key=lambda component: component[0])
    atom_to_mol: dict[int, int] = {}
    for mol_id, component in enumerate(ordered_components, start=1):
        for atom_id in component:
            atom_to_mol[atom_id] = mol_id
    return atom_to_mol


def correct_vmd_lammps_files(
    *,
    data_path: str | Path,
    dump_path: str | Path | None = None,
    angle: bool = False,
) -> CorrectedVmdSummary:
    del angle  # kept only for CLI compatibility; atom records are auto-detected
    data = parse_lammps_data(data_path)
    atom_ids = sorted(atom.atom_id for atom in data.atoms)
    atom_types = {atom.atom_id: atom.atom_type for atom in data.atoms}
    atom_to_mol = assign_molecule_ids(atom_ids, [(bond.atom1, bond.atom2) for bond in data.bonds])

    corrected_atoms = []
    for atom in sorted(data.atoms, key=lambda item: item.atom_id):
        corrected_atoms.append(
            atom.__class__(
                atom_id=atom.atom_id,
                mol_id=atom_to_mol[atom.atom_id],
                atom_type=atom.atom_type,
                charge=atom.charge,
                x=atom.x,
                y=atom.y,
                z=atom.z,
                ix=atom.ix,
                iy=atom.iy,
                iz=atom.iz,
                style=atom.style,
                comment=atom.comment,
            )
        )

    data_path = Path(data_path)
    corrected_data_path = data_path.with_name(f'{data_path.name}-corrected')
    corrected_data_path.write_text(data.render(corrected_atoms))

    corrected_dump_path: Path | None = None
    frames = 0
    if dump_path is not None:
        dump_path = Path(dump_path)
        corrected_dump_path = dump_path.with_name(f'{dump_path.name}-corrected')
        frame_chunks: list[str] = []
        for frame in iter_lammps_dump_frames(dump_path):
            frames += 1
            if frame.atom_count != len(atom_ids):
                raise ValueError(
                    f'Incompatible data and dump files: data has {len(atom_ids)} atoms, dump frame has {frame.atom_count}'
                )
            coord_fields = choose_coordinate_fields(frame.atom_fields)
            frame_chunks.extend([
                'ITEM: TIMESTEP',
                frame.timestep,
                'ITEM: NUMBER OF ATOMS',
                str(frame.atom_count),
                frame.box_header,
                *frame.box_lines,
                f"ITEM: ATOMS id mol type {' '.join(coord_fields)}",
            ])
            for atom_values in frame.atoms:
                atom_id = int(atom_values['id'])
                if atom_id not in atom_to_mol:
                    raise ValueError(f'Dump frame references atom id {atom_id} that is missing from the data file')
                coords = ' '.join(atom_values[field] for field in coord_fields)
                frame_chunks.append(f'{atom_id} {atom_to_mol[atom_id]} {atom_types[atom_id]} {coords}')
        corrected_dump_path.write_text('\n'.join(frame_chunks) + ('\n' if frame_chunks else ''))

    return CorrectedVmdSummary(
        atoms=len(atom_ids),
        bonds=len(data.bonds),
        max_atom_id=max(atom_ids) if atom_ids else 0,
        chains=len(set(atom_to_mol.values())),
        corrected_data_path=corrected_data_path,
        corrected_dump_path=corrected_dump_path,
        frames=frames,
    )
