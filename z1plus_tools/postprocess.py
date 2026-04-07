from __future__ import annotations

import math
from pathlib import Path

from .formats import Chain, Frame, LammpsAtom, LammpsBond, Node
from .parsers import parse_dat_frames, parse_z1_frames, round_away_from_zero
from .writers import (
    bounds_from_atoms,
    render_dump_frame,
    render_dump_from_atoms,
    render_lammps_data,
    symmetric_box_bounds,
    wrap_coordinate,
)


def _require_single_frame(path: str | Path, parser) -> Frame:
    path = Path(path)
    frames = parser(path)
    if not frames:
        raise ValueError(f"{path} does not contain any frames")
    if len(frames) != 1:
        raise ValueError(f"{path} contains {len(frames)} frames; this command expects exactly one")
    return frames[0]


def _extra_int(extras: tuple[str, ...], index: int, default: int = 0) -> int:
    if index >= len(extras):
        return default
    token = extras[index]
    if token == "":
        return default
    try:
        return int(float(token))
    except ValueError:
        return default


def _append_chain_records(
    *,
    frame: Frame,
    atoms: list[LammpsAtom],
    bonds: list[LammpsBond],
    first_type: int,
    interior_type: int,
    last_type: int,
    bond_type: int,
    mol_offset: int = 0,
    wrap_coordinates: bool = True,
    include_images: bool = True,
    comment_for=None,
) -> dict[int, int]:
    boxx, boxy, boxz = frame.box
    first_ids: dict[int, int] = {}
    atom_id = len(atoms)
    bond_id = len(bonds)
    for chain_index, chain in enumerate(frame.chains, start=1):
        first_ids[chain_index] = atom_id + 1
        for node_index, node in enumerate(chain.nodes, start=1):
            atom_id += 1
            if node_index == 1:
                atom_type = first_type
            elif node_index == len(chain.nodes):
                atom_type = last_type
            else:
                atom_type = interior_type
            if wrap_coordinates:
                x, ix = wrap_coordinate(node.x, boxx)
                y, iy = wrap_coordinate(node.y, boxy)
                z, iz = wrap_coordinate(node.z, boxz)
            else:
                x, y, z = node.x, node.y, node.z
                ix = iy = iz = 0
            atoms.append(
                LammpsAtom(
                    atom_id=atom_id,
                    mol_id=chain_index + mol_offset,
                    atom_type=atom_type,
                    x=x,
                    y=y,
                    z=z,
                    ix=ix if include_images else None,
                    iy=iy if include_images else None,
                    iz=iz if include_images else None,
                    comment=None if comment_for is None else comment_for(chain_index, node_index, node),
                )
            )
            if node_index > 1:
                bond_id += 1
                bonds.append(LammpsBond(bond_id=bond_id, bond_type=bond_type, atom1=atom_id - 1, atom2=atom_id))
    return first_ids


def _unwrap_chain(chain: Chain, box: tuple[float, float, float]) -> list[Node]:
    if not chain.nodes:
        return []
    boxx, boxy, boxz = box
    unwrapped = [Node(chain.nodes[0].x, chain.nodes[0].y, chain.nodes[0].z, chain.nodes[0].extras)]
    prev_raw = chain.nodes[0]
    prev_unwrapped = unwrapped[0]
    for node in chain.nodes[1:]:
        dx = node.x - prev_raw.x
        dy = node.y - prev_raw.y
        dz = node.z - prev_raw.z
        dx -= boxx * round_away_from_zero(dx / boxx)
        dy -= boxy * round_away_from_zero(dy / boxy)
        dz -= boxz * round_away_from_zero(dz / boxz)
        current = Node(
            prev_unwrapped.x + dx,
            prev_unwrapped.y + dy,
            prev_unwrapped.z + dz,
            node.extras,
        )
        unwrapped.append(current)
        prev_raw = node
        prev_unwrapped = current
    return unwrapped


def _shift_nodes(nodes: list[Node], shift: tuple[float, float, float]) -> list[Node]:
    sx, sy, sz = shift
    return [Node(node.x + sx, node.y + sy, node.z + sz, node.extras) for node in nodes]


def _emit_chain(
    *,
    atoms: list[LammpsAtom],
    bonds: list[LammpsBond],
    mol_id: int,
    atom_type: int,
    bond_type: int,
    nodes: list[Node],
    box: tuple[float, float, float],
    folded: bool,
    dat_lines: list[str] | None = None,
    dat_extra=None,
) -> tuple[int, int]:
    boxx, boxy, boxz = box
    atom_id = len(atoms)
    bond_id = len(bonds)
    first_atom_id = atom_id + 1
    last_atom_id = first_atom_id
    for node_index, node in enumerate(nodes, start=1):
        atom_id += 1
        last_atom_id = atom_id
        if folded:
            x, ix = wrap_coordinate(node.x, boxx)
            y, iy = wrap_coordinate(node.y, boxy)
            z, iz = wrap_coordinate(node.z, boxz)
        else:
            x, y, z = node.x, node.y, node.z
            ix = iy = iz = 0
        atoms.append(
            LammpsAtom(
                atom_id=atom_id,
                mol_id=mol_id,
                atom_type=atom_type,
                x=x,
                y=y,
                z=z,
                ix=ix,
                iy=iy,
                iz=iz,
            )
        )
        if dat_lines is not None:
            line = f"{x:.15g} {y:.15g} {z:.15g}"
            if dat_extra is not None:
                extra_text = dat_extra(node_index, node, len(nodes))
                if extra_text:
                    line += f" {extra_text}"
            dat_lines.append(line)
        if node_index > 1:
            bond_id += 1
            bonds.append(LammpsBond(bond_id=bond_id, bond_type=bond_type, atom1=atom_id - 1, atom2=atom_id))
    return first_atom_id, last_atom_id


def _sp_extra_text(reverse_map: dict[int, int], node_index: int, node: Node, chain_length: int) -> str:
    pos = node.extras[0] if len(node.extras) > 0 else "0"
    ent = node.extras[1] if len(node.extras) > 1 else "0"
    ent_chain = _extra_int(node.extras, 2, 0)
    ent_bead = _extra_int(node.extras, 3, 0)
    mapped = reverse_map.get(ent_chain)
    if mapped is None:
        mapped_text = "-1" if node_index < chain_length else "0"
    else:
        mapped_text = str(mapped)
    return f"{pos} {ent} {mapped_text} {ent_bead}"


def _choose_shift(selected_node: Node, entangled_sp_chain: Chain, entangled_bead: int, box: tuple[float, float, float]) -> tuple[float, float, float]:
    boxx, boxy, boxz = box
    candidates: list[tuple[float, float, float, float]] = []
    for bead_offset in (0, 1):
        bead_index = entangled_bead - 1 + bead_offset
        if not (0 <= bead_index < len(entangled_sp_chain.nodes)):
            continue
        node = entangled_sp_chain.nodes[bead_index]
        dx = node.x - selected_node.x
        dy = node.y - selected_node.y
        dz = node.z - selected_node.z
        shiftx = boxx * round_away_from_zero(dx / boxx)
        shifty = boxy * round_away_from_zero(dy / boxy)
        shiftz = boxz * round_away_from_zero(dz / boxz)
        dx -= shiftx
        dy -= shifty
        dz -= shiftz
        distance = math.sqrt(dx * dx + dy * dy + dz * dz)
        candidates.append((distance, -shiftx, -shifty, -shiftz))
    if not candidates:
        return (0.0, 0.0, 0.0)
    _, sx, sy, sz = min(candidates, key=lambda item: item[0])
    return (sx, sy, sz)


def convert_z1_file_to_dump_text(path: str | Path, *, unfolded: bool = False) -> str:
    frames = parse_z1_frames(path)
    return "".join(render_dump_frame(frame, index, unfolded=unfolded) for index, frame in enumerate(frames, start=1))


def convert_dat_file_to_dump_text(path: str | Path, *, unfolded: bool = False) -> str:
    frames = parse_dat_frames(path)
    return "".join(render_dump_frame(frame, index, unfolded=unfolded) for index, frame in enumerate(frames, start=1))


def create_sp_to_data_outputs(
    *,
    export_initconfig: bool = False,
    export_sp: bool = False,
    export_merge: bool = False,
    initconfig_path: str | Path = "Z1+initconfig.dat",
    sp_path: str | Path = "Z1+SP.dat",
) -> list[Path]:
    created: list[Path] = []
    init_frame = None
    init_first_ids: dict[int, int] | None = None

    if export_initconfig or export_merge:
        init_frame = _require_single_frame(initconfig_path, parse_dat_frames)
        init_atoms: list[LammpsAtom] = []
        init_bonds: list[LammpsBond] = []
        init_first_ids = _append_chain_records(
            frame=init_frame,
            atoms=init_atoms,
            bonds=init_bonds,
            first_type=1,
            interior_type=2,
            last_type=3,
            bond_type=1,
            wrap_coordinates=True,
            include_images=True,
        )
        if export_initconfig:
            init_path = Path(initconfig_path).with_suffix(".data")
            init_path.write_text(
                render_lammps_data(
                    title="lammps data file created by Z1+SP-to-data.py (software at https://github.com/mkmat/Z1plus-code)",
                    bounds=symmetric_box_bounds(init_frame.box),
                    atoms=init_atoms,
                    bonds=init_bonds,
                    atom_types=6,
                    bond_types=2,
                    atoms_header="Atoms",
                    include_images=True,
                )
            )
            created.append(init_path)

    sp_frame = None
    if export_sp or export_merge:
        sp_frame = _require_single_frame(sp_path, parse_dat_frames)
        if init_frame is not None and len(init_frame.chains) != len(sp_frame.chains):
            raise ValueError("Z1+initconfig.dat and Z1+SP.dat contain different numbers of chains")
        sp_atoms: list[LammpsAtom] = []
        sp_bonds: list[LammpsBond] = []
        _append_chain_records(
            frame=sp_frame,
            atoms=sp_atoms,
            bonds=sp_bonds,
            first_type=4,
            interior_type=5,
            last_type=6,
            bond_type=1,
            wrap_coordinates=True,
            include_images=True,
        )
        if export_sp:
            sp_out = Path(sp_path).with_suffix(".data")
            sp_out.write_text(
                render_lammps_data(
                    title="lammps data file created by Z1+SP-to-data.py (software at https://github.com/mkmat/Z1plus-code)",
                    bounds=symmetric_box_bounds(sp_frame.box),
                    atoms=sp_atoms,
                    bonds=sp_bonds,
                    atom_types=6,
                    bond_types=2,
                    atoms_header="Atoms",
                    include_images=True,
                )
            )
            created.append(sp_out)

    if export_merge:
        if init_frame is None or sp_frame is None or init_first_ids is None:
            raise ValueError("Merge export requires both Z1+initconfig.dat and Z1+SP.dat")

        merge_atoms: list[LammpsAtom] = []
        merge_bonds: list[LammpsBond] = []
        _append_chain_records(
            frame=init_frame,
            atoms=merge_atoms,
            bonds=merge_bonds,
            first_type=1,
            interior_type=2,
            last_type=3,
            bond_type=1,
            wrap_coordinates=True,
            include_images=True,
        )

        def comment_for(chain_index: int, node_index: int, node: Node) -> str | None:
            ent_chain = _extra_int(node.extras, 2, 0)
            ent_bead = _extra_int(node.extras, 3, 0)
            if ent_chain and ent_bead and ent_chain in init_first_ids:
                return str(init_first_ids[ent_chain] + ent_bead - 1)
            return None

        _append_chain_records(
            frame=sp_frame,
            atoms=merge_atoms,
            bonds=merge_bonds,
            first_type=4,
            interior_type=5,
            last_type=6,
            bond_type=2,
            mol_offset=len(init_frame.chains),
            wrap_coordinates=True,
            include_images=True,
            comment_for=comment_for,
        )

        merge_path = Path("Z1+merge.data")
        merge_path.write_text(
            render_lammps_data(
                title="lammps data file created by Z1+SP-to-data.py (software at https://github.com/mkmat/Z1plus-code)",
                bounds=symmetric_box_bounds(init_frame.box),
                atoms=merge_atoms,
                bonds=merge_bonds,
                atom_types=6,
                bond_types=2,
                atoms_header="Atoms",
                include_images=True,
            )
        )
        created.append(merge_path)

    return created


def export_merged_outputs(
    *,
    dump: bool = True,
    data: bool = False,
    from_frame: int = 1,
    to_frame: int = 1,
    export_all: bool = False,
    initconfig_path: str | Path = "Z1+initconfig.dat",
    sp_path: str | Path = "Z1+SP.dat",
) -> list[Path]:
    config_frames = parse_dat_frames(initconfig_path)
    sp_frames = parse_dat_frames(sp_path)
    if len(config_frames) != len(sp_frames):
        raise ValueError("Z1+initconfig.dat and Z1+SP.dat contain different numbers of frames")

    if export_all:
        from_frame = 1
        to_frame = len(config_frames)

    created: list[Path] = []
    dump_chunks: list[str] = []

    for frame_index, (config_frame, sp_frame) in enumerate(zip(config_frames, sp_frames), start=1):
        if frame_index < from_frame or frame_index > to_frame:
            continue
        atoms: list[LammpsAtom] = []
        bonds: list[LammpsBond] = []
        _append_chain_records(
            frame=config_frame,
            atoms=atoms,
            bonds=bonds,
            first_type=1,
            interior_type=2,
            last_type=3,
            bond_type=1,
            wrap_coordinates=False,
            include_images=False,
        )
        _append_chain_records(
            frame=sp_frame,
            atoms=atoms,
            bonds=bonds,
            first_type=4,
            interior_type=5,
            last_type=6,
            bond_type=2,
            wrap_coordinates=False,
            include_images=False,
        )
        bounds = symmetric_box_bounds(config_frame.box)
        if data:
            data_path = Path(f"Z1+merged-{frame_index}.data")
            data_path.write_text(
                render_lammps_data(
                    title="LAMMPS data file via Z1+export.py (mk@mat.ethz.ch)",
                    bounds=bounds,
                    atoms=atoms,
                    bonds=bonds,
                    atom_types=6,
                    bond_types=2,
                    atoms_header="Atoms # angle",
                    include_images=True,
                )
            )
            created.append(data_path)
        if dump:
            dump_chunks.append(render_dump_from_atoms(timestep=frame_index, bounds=bounds, atoms=atoms))

    if dump and dump_chunks:
        dump_path = Path("Z1+merged.dump")
        dump_path.write_text("".join(dump_chunks))
        created.append(dump_path)

    return created


def create_single_chain_entanglement_outputs(
    *,
    chain_id: int,
    folded: bool = False,
    txt: bool = False,
    dat: bool = False,
    add_sp: bool = False,
    add_ee: bool = False,
    output_path: str | Path | None = None,
    initconfig_path: str | Path = "Z1+initconfig.dat",
    sp_path: str | Path = "Z1+SP.dat",
) -> list[Path]:
    if dat:
        add_sp = True

    init_frame = _require_single_frame(initconfig_path, parse_dat_frames)
    sp_frame = _require_single_frame(sp_path, parse_dat_frames)
    if len(init_frame.chains) != len(sp_frame.chains):
        raise ValueError("Z1+initconfig.dat and Z1+SP.dat contain different numbers of chains")

    chains = len(init_frame.chains)
    if not (1 <= chain_id <= chains):
        raise ValueError(f"Chain id {chain_id} is out of range 1..{chains}")

    selected_index = chain_id - 1
    selected_sp = sp_frame.chains[selected_index]
    total_entanglements = 0
    for node in selected_sp.nodes[1:-1]:
        ent_flag = _extra_int(node.extras, 1, 0)
        ent_chain = _extra_int(node.extras, 2, 0)
        total_entanglements += ent_flag
        if ent_flag and not ent_chain:
            raise ValueError("Z1+SP.dat is missing entanglement chain ids; run Z1+ with -SP+ to generate the required metadata")
    if total_entanglements == 0:
        raise ValueError(f"Chain id {chain_id} does not have any entanglements")

    atoms: list[LammpsAtom] = []
    bonds: list[LammpsBond] = []
    ee_pairs: list[tuple[int, int]] = []
    original_blocks: list[str] = []
    sp_blocks: list[str] = []
    reverse_map: dict[int, int] = {chain_id: 1}
    box = init_frame.box

    selected_original_nodes = _unwrap_chain(init_frame.chains[selected_index], box)
    selected_dat_lines: list[str] = []
    first_atom, last_atom = _emit_chain(
        atoms=atoms,
        bonds=bonds,
        mol_id=chain_id,
        atom_type=1,
        bond_type=1,
        nodes=selected_original_nodes,
        box=box,
        folded=folded,
        dat_lines=selected_dat_lines,
    )
    ee_pairs.append((first_atom, last_atom))
    original_blocks.append("\n".join([str(len(selected_original_nodes)), *selected_dat_lines]))

    entangled_instances: list[tuple[int, tuple[float, float, float]]] = []
    for node in selected_sp.nodes[1:-1]:
        entangled_chain = _extra_int(node.extras, 2, 0)
        entangled_bead = _extra_int(node.extras, 3, 0)
        if not entangled_chain:
            continue
        shift = _choose_shift(node, sp_frame.chains[entangled_chain - 1], entangled_bead, sp_frame.box)
        entangled_instances.append((entangled_chain, shift))
        reverse_map[entangled_chain] = len(entangled_instances) + 1
        entangled_original_nodes = _shift_nodes(_unwrap_chain(init_frame.chains[entangled_chain - 1], box), shift)
        entangled_dat_lines: list[str] = []
        first_atom, last_atom = _emit_chain(
            atoms=atoms,
            bonds=bonds,
            mol_id=entangled_chain,
            atom_type=1,
            bond_type=1,
            nodes=entangled_original_nodes,
            box=box,
            folded=folded,
            dat_lines=entangled_dat_lines,
        )
        ee_pairs.append((first_atom, last_atom))
        original_blocks.append("\n".join([str(len(entangled_original_nodes)), *entangled_dat_lines]))

    if add_sp:
        selected_sp_dat_lines: list[str] = []
        _emit_chain(
            atoms=atoms,
            bonds=bonds,
            mol_id=chain_id + chains,
            atom_type=2,
            bond_type=2,
            nodes=selected_sp.nodes,
            box=box,
            folded=folded,
            dat_lines=selected_sp_dat_lines,
            dat_extra=lambda index, node, length: _sp_extra_text(reverse_map, index, node, length),
        )
        sp_blocks.append("\n".join([str(len(selected_sp.nodes)), *selected_sp_dat_lines]))
        for entangled_chain, shift in entangled_instances:
            entangled_sp_nodes = _shift_nodes(sp_frame.chains[entangled_chain - 1].nodes, shift)
            entangled_sp_dat_lines: list[str] = []
            _emit_chain(
                atoms=atoms,
                bonds=bonds,
                mol_id=entangled_chain + chains,
                atom_type=2,
                bond_type=2,
                nodes=entangled_sp_nodes,
                box=box,
                folded=folded,
                dat_lines=entangled_sp_dat_lines,
                dat_extra=lambda index, node, length: _sp_extra_text(reverse_map, index, node, length),
            )
            sp_blocks.append("\n".join([str(len(entangled_sp_nodes)), *entangled_sp_dat_lines]))

    bond_types = 2 if add_sp else 1
    if add_ee:
        bond_types = 3
        for atom1, atom2 in ee_pairs:
            bonds.append(LammpsBond(bond_id=len(bonds) + 1, bond_type=3, atom1=atom1, atom2=atom2))

    created: list[Path] = []
    default_data_path = Path(f"entangled-with-chain-{chain_id}.data") if output_path is None else Path(output_path)

    if not txt and not dat:
        bounds = symmetric_box_bounds(box) if folded else bounds_from_atoms(atoms)
        default_data_path.write_text(
            render_lammps_data(
                title=f"lammps data file generated via extract-single-chain-entanglements.py {chain_id}",
                bounds=bounds,
                atoms=atoms,
                bonds=bonds,
                atom_types=2 if add_sp else 1,
                bond_types=bond_types,
                atoms_header="Atoms",
                include_images=True,
            )
        )
        created.append(default_data_path)

    if dat:
        original_path = Path(f"Z1+initconfig-chain={chain_id}.dat")
        original_path.write_text(
            "\n".join([str(len(original_blocks)), f"{box[0]:.15g} {box[1]:.15g} {box[2]:.15g}", *original_blocks, ""])
        )
        sp_out_path = Path(f"Z1+SP-chain={chain_id}.dat")
        sp_out_path.write_text(
            "\n".join([str(len(sp_blocks)), f"{box[0]:.15g} {box[1]:.15g} {box[2]:.15g}", *sp_blocks, ""])
        )
        created.extend([original_path, sp_out_path])

    if txt:
        txt_path = Path(f"entangled-with-chain-{chain_id}.txt")
        atom_lines = [
            f"{atom.atom_id} {atom.mol_id} {atom.atom_type} {atom.x:.15g} {atom.y:.15g} {atom.z:.15g} "
            f"{0 if atom.ix is None else atom.ix} {0 if atom.iy is None else atom.iy} {0 if atom.iz is None else atom.iz}"
            for atom in atoms
        ]
        txt_path.write_text("\n".join(atom_lines) + "\n")
        created.append(txt_path)

    return created
