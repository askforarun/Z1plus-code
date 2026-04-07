from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import xml.etree.ElementTree as ET

from .lammps import LammpsDataAtom, LammpsDumpFrame, choose_coordinate_fields, iter_lammps_dump_frames, parse_lammps_data


@dataclass(slots=True)
class ImportedChains:
    chain_paths: list[list[int]]
    atom_to_newid: dict[int, int]
    atom_to_newmol: dict[int, int]
    chain_lengths: list[int]
    lammps_mol: dict[int, int]
    atom_type: dict[int, int]
    atom_charge: dict[int, float | None]


@dataclass(slots=True)
class DumpFrameData:
    timestep: int
    boxx: float
    boxy: float
    boxz: float
    xy: float
    xz: float
    yz: float
    coords_by_atom_id: dict[int, tuple[float, float, float]]
    mol_by_atom_id: dict[int, int]
    type_by_atom_id: dict[int, int]
    order: list[int]


@dataclass(slots=True)
class ImportSummary:
    frames_written: int
    chains: int
    atoms: int
    config_path: Path
    table_path: Path | None = None
    out_dump_path: Path | None = None


class ImportError(RuntimeError):
    pass


def _mass_is_hydrogen(mass: float | None) -> bool:
    return mass is not None and abs(round(mass) - 1.0) < 1e-9


def _build_adjacency(atom_ids: list[int], bonds: list[tuple[int, int]]) -> dict[int, list[int]]:
    adjacency = {atom_id: [] for atom_id in atom_ids}
    for atom1, atom2 in bonds:
        if atom1 not in adjacency or atom2 not in adjacency:
            continue
        adjacency[atom1].append(atom2)
        adjacency[atom2].append(atom1)
    return adjacency


def _cut_branched_connections(adjacency: dict[int, list[int]]) -> int:
    deleted = 0
    for atom_id in sorted(adjacency):
        while len(adjacency[atom_id]) > 2:
            neighbor = adjacency[atom_id].pop()
            adjacency[neighbor] = [candidate for candidate in adjacency[neighbor] if candidate != atom_id]
            deleted += 1
    return deleted


def _components(adjacency: dict[int, list[int]]) -> list[list[int]]:
    seen: set[int] = set()
    components: list[list[int]] = []
    for atom_id in sorted(adjacency):
        if atom_id in seen or not adjacency[atom_id]:
            continue
        stack = [atom_id]
        component: list[int] = []
        seen.add(atom_id)
        while stack:
            current = stack.pop()
            component.append(current)
            for neighbor in adjacency[current]:
                if neighbor not in seen:
                    seen.add(neighbor)
                    stack.append(neighbor)
        components.append(sorted(component))
    return components


def _walk_linear_component(adjacency: dict[int, list[int]], component: list[int], start: int) -> list[int]:
    path = [start]
    prev = 0
    current = start
    while True:
        next_candidates = [neighbor for neighbor in adjacency[current] if neighbor != prev]
        if not next_candidates:
            break
        if len(next_candidates) > 1:
            raise ImportError(f'Non-linear component encountered while walking chain from atom {start}')
        prev, current = current, next_candidates[0]
        path.append(current)
    return path


def _build_linear_chains(
    *,
    atom_ids: list[int],
    bonds: list[tuple[int, int]],
    ignore_dumbbells: bool,
    allow_branched: bool,
) -> list[list[int]]:
    adjacency = _build_adjacency(atom_ids, bonds)
    deleted_bonds = _cut_branched_connections(adjacency) if allow_branched else 0
    if deleted_bonds:
        print(f'-branched option: {deleted_bonds} deleted bonds')
    for atom_id, neighbors in adjacency.items():
        if len(neighbors) > 2:
            raise ImportError(
                'Branched structure detected. Use the Python backbone extractor for atomistic models, '
                'rerun with -ignore_H, or pass -branched to split sidechains into separate chains.'
            )

    chain_paths: list[list[int]] = []
    for component in _components(adjacency):
        endpoints = [atom_id for atom_id in component if len(adjacency[atom_id]) == 1]
        if len(component) == 2 and len(endpoints) == 2 and ignore_dumbbells:
            continue
        if len(endpoints) != 2:
            raise ImportError(
                f'Component containing atoms {component[:6]} has {len(endpoints)} terminal atoms; '
                'Z1+import-lammps expects linear chains after filtering.'
            )
        chain_paths.append(_walk_linear_component(adjacency, component, min(endpoints)))
    chain_paths.sort(key=lambda path: path[0])
    return chain_paths


def _imported_chains_from_data(
    data_path: str | Path,
    *,
    ignore_h: bool,
    ignore_types: set[int],
    ignore_dumbbells: bool,
    allow_branched: bool,
) -> tuple[ImportedChains, dict[int, LammpsDataAtom], tuple[float, float, float, float, float, float, float, float, float]]:
    data = parse_lammps_data(data_path)
    atom_map = {atom.atom_id: atom for atom in data.atoms}
    active_ids = [
        atom.atom_id
        for atom in sorted(data.atoms, key=lambda item: item.atom_id)
        if atom.atom_type not in ignore_types
        and not (ignore_h and _mass_is_hydrogen(data.masses.get(atom.atom_type)))
    ]
    active_set = set(active_ids)
    bonds = [
        (bond.atom1, bond.atom2)
        for bond in data.bonds
        if bond.atom1 in active_set and bond.atom2 in active_set
    ]
    chain_paths = _build_linear_chains(
        atom_ids=active_ids,
        bonds=bonds,
        ignore_dumbbells=ignore_dumbbells,
        allow_branched=allow_branched,
    )
    atom_to_newid: dict[int, int] = {}
    atom_to_newmol: dict[int, int] = {}
    chain_lengths: list[int] = []
    next_id = 1
    for chain_index, path in enumerate(chain_paths, start=1):
        chain_lengths.append(len(path))
        for atom_id in path:
            atom_to_newid[atom_id] = next_id
            atom_to_newmol[atom_id] = chain_index
            next_id += 1
    imported = ImportedChains(
        chain_paths=chain_paths,
        atom_to_newid=atom_to_newid,
        atom_to_newmol=atom_to_newmol,
        chain_lengths=chain_lengths,
        lammps_mol={atom_id: atom_map[atom_id].mol_id for atom_id in atom_to_newid},
        atom_type={atom_id: atom_map[atom_id].atom_type for atom_id in atom_to_newid},
        atom_charge={atom_id: atom_map[atom_id].charge for atom_id in atom_to_newid},
    )
    box = data.box
    return imported, atom_map, (box.xlo, box.xhi, box.ylo, box.yhi, box.zlo, box.zhi, box.xy, box.xz, box.yz)


def _parse_dump_box(frame: LammpsDumpFrame) -> tuple[float, float, float, float, float, float, float, float, float]:
    tokens = frame.box_header.split()
    numeric_lines = [tuple(float(value) for value in line.split()) for line in frame.box_lines]
    if 'xy' in tokens:
        xlobound, xhibound, xy = numeric_lines[0]
        ylobound, yhibound, xz = numeric_lines[1]
        zlobound, zhibound, yz = numeric_lines[2]
        xlo = xlobound - min(0.0, xy)
        xhi = xhibound - max(0.0, xy)
        ylo, yhi = ylobound, yhibound
        zlo, zhi = zlobound, zhibound
        return xlo, xhi, ylo, yhi, zlo, zhi, xy, xz, yz
    xlo, xhi = numeric_lines[0][:2]
    ylo, yhi = numeric_lines[1][:2]
    zlo, zhi = numeric_lines[2][:2]
    return xlo, xhi, ylo, yhi, zlo, zhi, 0.0, 0.0, 0.0


def _read_dump_frame(frame: LammpsDumpFrame) -> DumpFrameData:
    xlo, xhi, ylo, yhi, zlo, zhi, xy, xz, yz = _parse_dump_box(frame)
    coord_fields = choose_coordinate_fields(frame.atom_fields)
    if coord_fields[0] == 'xs':
        scalex, scaley, scalez = xhi - xlo, yhi - ylo, zhi - zlo
    else:
        scalex = scaley = scalez = 1.0
    coords_by_atom_id: dict[int, tuple[float, float, float]] = {}
    mol_by_atom_id: dict[int, int] = {}
    type_by_atom_id: dict[int, int] = {}
    order: list[int] = []
    for atom_values in frame.atoms:
        atom_id = int(atom_values['id'])
        order.append(atom_id)
        mol_by_atom_id[atom_id] = int(atom_values.get('mol', 0))
        type_by_atom_id[atom_id] = int(atom_values.get('type', 0))
        coords_by_atom_id[atom_id] = (
            float(atom_values[coord_fields[0]]) * scalex,
            float(atom_values[coord_fields[1]]) * scaley,
            float(atom_values[coord_fields[2]]) * scalez,
        )
    return DumpFrameData(
        timestep=int(frame.timestep),
        boxx=xhi - xlo,
        boxy=yhi - ylo,
        boxz=zhi - zlo,
        xy=xy,
        xz=xz,
        yz=yz,
        coords_by_atom_id=coords_by_atom_id,
        mol_by_atom_id=mol_by_atom_id,
        type_by_atom_id=type_by_atom_id,
        order=order,
    )


def _render_z1_frame(
    chain_paths: list[list[int]],
    coords_by_atom_id: dict[int, tuple[float, float, float]],
    *,
    boxx: float,
    boxy: float,
    boxz: float,
    xy: float,
    timestep: int,
    z1_mode: bool,
) -> str:
    lines = [str(len(chain_paths)), f'{boxx} {boxy} {boxz}', ' '.join(str(len(path)) for path in chain_paths)]
    for path in chain_paths:
        for atom_id in path:
            x, y, z = coords_by_atom_id[atom_id]
            lines.append(f'{x} {y} {z}')
    if z1_mode:
        if xy:
            lines.append(str(xy))
    else:
        lines.extend(['-2', str(xy), str(timestep)])
    return '\n'.join(lines) + '\n'


def _write_mapping_table(imported: ImportedChains, path: str | Path = 'table-lammps-ids-to-Z1-id-mol-bead.txt') -> Path:
    path = Path(path)
    lines: list[str] = []
    for chain_index, chain in enumerate(imported.chain_paths, start=1):
        for bead_index, atom_id in enumerate(chain, start=1):
            lines.append(
                f'{atom_id} {imported.lammps_mol[atom_id]} {imported.atom_to_newid[atom_id]} {chain_index} {bead_index}'
            )
    path.write_text('\n'.join(lines) + ('\n' if lines else ''))
    return path


def _write_truechain_dump(
    imported: ImportedChains,
    coords_by_atom_id: dict[int, tuple[float, float, float]],
    *,
    xlo: float,
    xhi: float,
    ylo: float,
    yhi: float,
    zlo: float,
    zhi: float,
    xy: float,
    out_path: str | Path,
) -> Path:
    out_path = Path(out_path)
    lines = ['ITEM: TIMESTEP', '0', 'ITEM: NUMBER OF ATOMS', str(sum(imported.chain_lengths))]
    if xy:
        lines.extend(['ITEM: BOX BOUNDS xy', f'{xlo} {xhi} {xy}', f'{ylo} {yhi} 0', f'{zlo} {zhi} 0'])
    else:
        lines.extend(['ITEM: BOX BOUNDS', f'{xlo} {xhi}', f'{ylo} {yhi}', f'{zlo} {zhi}'])
    lines.append('ITEM: ATOMS id mol x y z')
    for chain_index, path in enumerate(imported.chain_paths, start=1):
        for atom_id in path:
            new_id = imported.atom_to_newid[atom_id]
            x, y, z = coords_by_atom_id[atom_id]
            lines.append(f'{new_id} {chain_index} {x} {y} {z}')
    out_path.write_text('\n'.join(lines) + '\n')
    return out_path


def _chains_from_dump_only(frame: DumpFrameData, *, ignore_dumbbells: bool) -> ImportedChains:
    if frame.order != list(range(1, len(frame.order) + 1)):
        raise ImportError('atoms not saved using dump_modify sort id')
    chain_paths: list[list[int]] = []
    current_mol = None
    current_path: list[int] = []
    last_mol = 0
    for atom_id in frame.order:
        mol = frame.mol_by_atom_id[atom_id]
        if mol < last_mol or mol > last_mol + 1:
            raise ImportError('chains not saved with increasing mol number')
        if current_mol is None or mol != current_mol:
            if current_path and not (ignore_dumbbells and len(current_path) == 2):
                chain_paths.append(current_path)
            current_mol = mol
            current_path = [atom_id]
            last_mol = mol
        else:
            current_path.append(atom_id)
    if current_path and not (ignore_dumbbells and len(current_path) == 2):
        chain_paths.append(current_path)
    atom_to_newid: dict[int, int] = {}
    atom_to_newmol: dict[int, int] = {}
    next_id = 1
    for chain_index, path in enumerate(chain_paths, start=1):
        for atom_id in path:
            atom_to_newid[atom_id] = next_id
            atom_to_newmol[atom_id] = chain_index
            next_id += 1
    return ImportedChains(
        chain_paths=chain_paths,
        atom_to_newid=atom_to_newid,
        atom_to_newmol=atom_to_newmol,
        chain_lengths=[len(path) for path in chain_paths],
        lammps_mol={atom_id: frame.mol_by_atom_id[atom_id] for atom_id in atom_to_newid},
        atom_type={atom_id: frame.type_by_atom_id[atom_id] for atom_id in atom_to_newid},
        atom_charge={atom_id: None for atom_id in atom_to_newid},
    )


def _import_xml(xml_path: str | Path) -> tuple[list[list[int]], dict[int, tuple[float, float, float]], tuple[float, float, float]]:
    root = ET.parse(xml_path).getroot()
    box_element = root.find('.//box')
    if box_element is None:
        raise ImportError('XML input is missing a <box> element')
    boxx = float(box_element.attrib['lx'])
    boxy = float(box_element.attrib['ly'])
    boxz = float(box_element.attrib['lz'])
    positions = root.find('.//position')
    bonds = root.find('.//bond')
    if positions is None or bonds is None:
        raise ImportError('XML input must contain <position> and <bond> sections')
    atoms = [tuple(float(value) for value in line.split()) for line in positions.text.splitlines() if line.strip()]
    edge_list: list[tuple[int, int]] = []
    for line in bonds.text.splitlines():
        if not line.strip():
            continue
        _btype, atom1, atom2 = line.split()[:3]
        edge_list.append((int(atom1), int(atom2)))
    chain_paths = _build_linear_chains(
        atom_ids=list(range(len(atoms))),
        bonds=edge_list,
        ignore_dumbbells=False,
        allow_branched=False,
    )
    coords = {atom_id: atoms[atom_id] for atom_id in range(len(atoms))}
    return chain_paths, coords, (boxx, boxy, boxz)


def import_lammps_to_z1(
    *,
    dump_path: str | Path | None = None,
    data_path: str | Path | None = None,
    xml_path: str | Path | None = None,
    from_frame: int = 1,
    to_frame: int = 10**10,
    each: int = 1,
    ignore_h: bool = False,
    ignore_dumbbells: bool = False,
    ignore_types: tuple[int, ...] = (),
    z1_mode: bool = False,
    allow_branched: bool = False,
    out_dump_path: str | Path | None = None,
    verbose: bool = False,
    config_path: str | Path = 'config.Z1',
) -> ImportSummary:
    config_path = Path(config_path)
    if xml_path is not None:
        chain_paths, coords_by_atom_id, (boxx, boxy, boxz) = _import_xml(xml_path)
        config_path.write_text(_render_z1_frame(chain_paths, coords_by_atom_id, boxx=boxx, boxy=boxy, boxz=boxz, xy=0.0, timestep=0, z1_mode=True))
        return ImportSummary(frames_written=1, chains=len(chain_paths), atoms=sum(len(path) for path in chain_paths), config_path=config_path)

    table_path: Path | None = None
    imported: ImportedChains | None = None
    atom_map: dict[int, LammpsDataAtom] = {}
    data_box: tuple[float, float, float, float, float, float, float, float, float] | None = None
    if data_path is not None:
        imported, atom_map, data_box = _imported_chains_from_data(
            data_path,
            ignore_h=ignore_h,
            ignore_types=set(ignore_types),
            ignore_dumbbells=ignore_dumbbells,
            allow_branched=allow_branched,
        )
        table_path = _write_mapping_table(imported)
        if verbose:
            print(f'{sum(imported.chain_lengths)} atoms in {len(imported.chain_lengths)} chains after filtering')

    if dump_path is None and data_path is not None:
        assert imported is not None and data_box is not None
        xlo, xhi, ylo, yhi, zlo, zhi, xy, xz, yz = data_box
        if xz or yz:
            raise ImportError('data-only import currently supports orthorhombic or xy-tilted boxes only')
        coords_by_atom_id = {atom_id: (atom_map[atom_id].x, atom_map[atom_id].y, atom_map[atom_id].z) for atom_id in imported.atom_to_newid}
        config_path.write_text(
            _render_z1_frame(
                imported.chain_paths,
                coords_by_atom_id,
                boxx=xhi - xlo,
                boxy=yhi - ylo,
                boxz=zhi - zlo,
                xy=xy,
                timestep=0,
                z1_mode=z1_mode,
            )
        )
        out_dump_written = None
        if out_dump_path is not None:
            out_dump_written = _write_truechain_dump(
                imported,
                coords_by_atom_id,
                xlo=xlo,
                xhi=xhi,
                ylo=ylo,
                yhi=yhi,
                zlo=zlo,
                zhi=zhi,
                xy=xy,
                out_path=out_dump_path,
            )
        return ImportSummary(
            frames_written=1,
            chains=len(imported.chain_paths),
            atoms=sum(imported.chain_lengths),
            config_path=config_path,
            table_path=table_path,
            out_dump_path=out_dump_written,
        )

    frame_chunks: list[str] = []
    frames_written = 0
    dump_frames = iter_lammps_dump_frames(dump_path) if dump_path is not None else []
    for snapshot_index, raw_frame in enumerate(dump_frames, start=1):
        frame = _read_dump_frame(raw_frame)
        if snapshot_index < from_frame or snapshot_index > to_frame or ((snapshot_index - from_frame) % each) != 0:
            continue
        if imported is None:
            current_imported = _chains_from_dump_only(frame, ignore_dumbbells=ignore_dumbbells)
            current_paths = current_imported.chain_paths
            coords_by_atom_id = frame.coords_by_atom_id
            current_chains = current_imported
        else:
            if frame.xz or frame.yz:
                raise ImportError('dump import currently supports orthorhombic or xy-tilted boxes only')
            current_imported = imported
            current_paths = imported.chain_paths
            coords_by_atom_id = {atom_id: frame.coords_by_atom_id[atom_id] for atom_id in imported.atom_to_newid}
            current_chains = imported
        frame_chunks.append(
            _render_z1_frame(
                current_paths,
                coords_by_atom_id,
                boxx=frame.boxx,
                boxy=frame.boxy,
                boxz=frame.boxz,
                xy=frame.xy,
                timestep=frame.timestep,
                z1_mode=z1_mode,
            )
        )
        frames_written += 1
        if verbose:
            print(f'processed dump frame {snapshot_index} timestep {frame.timestep}')
        if out_dump_path is not None and imported is not None and frames_written == 1:
            _write_truechain_dump(
                current_chains,
                coords_by_atom_id,
                xlo=0.0,
                xhi=frame.boxx,
                ylo=0.0,
                yhi=frame.boxy,
                zlo=0.0,
                zhi=frame.boxz,
                xy=frame.xy,
                out_path=out_dump_path,
            )
    if not frame_chunks:
        raise ImportError('No frames matched the requested snapshot range')
    config_path.write_text(''.join(frame_chunks))
    if imported is None:
        table_path = None
        chains = len(current_imported.chain_paths)
        atoms = sum(current_imported.chain_lengths)
    else:
        chains = len(imported.chain_paths)
        atoms = sum(imported.chain_lengths)
    return ImportSummary(
        frames_written=frames_written,
        chains=chains,
        atoms=atoms,
        config_path=config_path,
        table_path=table_path,
        out_dump_path=Path(out_dump_path) if out_dump_path is not None else None,
    )
