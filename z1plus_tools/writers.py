from __future__ import annotations

from collections.abc import Iterable

from .formats import Frame, LammpsAtom, LammpsBond
from .parsers import round_away_from_zero


def format_number(value: float | int) -> str:
    if isinstance(value, int):
        return str(value)
    return f"{value:.15g}"


def wrap_coordinate(value: float, box_length: float) -> tuple[float, int]:
    image = round_away_from_zero(value / box_length)
    wrapped = value - box_length * image
    return wrapped, image


def symmetric_box_bounds(box: tuple[float, float, float]) -> tuple[float, float, float, float, float, float]:
    hx, hy, hz = (dimension / 2.0 for dimension in box)
    return (-hx, hx, -hy, hy, -hz, hz)


def bounds_from_atoms(atoms: Iterable[LammpsAtom]) -> tuple[float, float, float, float, float, float]:
    atoms = list(atoms)
    if not atoms:
        raise ValueError("Cannot compute bounds for an empty atom list")
    xs = [atom.x for atom in atoms]
    ys = [atom.y for atom in atoms]
    zs = [atom.z for atom in atoms]
    return (min(xs), max(xs), min(ys), max(ys), min(zs), max(zs))


def render_dump_frame(frame: Frame, timestep: int, unfolded: bool = False) -> str:
    boxx, boxy, boxz = frame.box
    xlo, xhi, ylo, yhi, zlo, zhi = symmetric_box_bounds(frame.box)
    atom_lines: list[str] = []
    atom_id = 0
    for mol_id, chain in enumerate(frame.chains, start=1):
        for node_index, node in enumerate(chain.nodes, start=1):
            atom_id += 1
            atom_type = 2 if node_index in (1, len(chain.nodes)) else 1
            if unfolded:
                x, y, z = node.x, node.y, node.z
            else:
                x, _ = wrap_coordinate(node.x, boxx)
                y, _ = wrap_coordinate(node.y, boxy)
                z, _ = wrap_coordinate(node.z, boxz)
            atom_lines.append(
                f"{atom_id} {mol_id} {atom_type} "
                f"{format_number(x)} {format_number(y)} {format_number(z)}"
            )
    coord_fields = "xu yu zu" if unfolded else "x y z"
    return "\n".join(
        [
            "ITEM: TIMESTEP",
            str(timestep),
            "ITEM: NUMBER OF ATOMS",
            str(atom_id),
            "ITEM: BOX BOUNDS pp pp pp",
            f"{format_number(xlo)} {format_number(xhi)}",
            f"{format_number(ylo)} {format_number(yhi)}",
            f"{format_number(zlo)} {format_number(zhi)}",
            f"ITEM: ATOMS id mol type {coord_fields}",
            *atom_lines,
        ]
    ) + "\n"


def render_dump_from_atoms(
    *,
    timestep: int,
    bounds: tuple[float, float, float, float, float, float],
    atoms: list[LammpsAtom],
) -> str:
    xlo, xhi, ylo, yhi, zlo, zhi = bounds
    atom_lines = [
        f"{atom.atom_id} {atom.mol_id} {atom.atom_type} "
        f"{format_number(atom.x)} {format_number(atom.y)} {format_number(atom.z)}"
        for atom in atoms
    ]
    return "\n".join(
        [
            "ITEM: TIMESTEP",
            str(timestep),
            "ITEM: NUMBER OF ATOMS",
            str(len(atoms)),
            "ITEM: BOX BOUNDS pp pp pp",
            f"{format_number(xlo)} {format_number(xhi)}",
            f"{format_number(ylo)} {format_number(yhi)}",
            f"{format_number(zlo)} {format_number(zhi)}",
            "ITEM: ATOMS id mol type x y z",
            *atom_lines,
        ]
    ) + "\n"


def render_lammps_data(
    *,
    title: str,
    bounds: tuple[float, float, float, float, float, float],
    atoms: list[LammpsAtom],
    bonds: list[LammpsBond],
    atom_types: int,
    bond_types: int,
    atoms_header: str = "Atoms",
    include_images: bool = True,
) -> str:
    xlo, xhi, ylo, yhi, zlo, zhi = bounds
    lines = [
        title,
        "",
        f"{len(atoms)} atoms",
        f"{len(bonds)} bonds",
        f"{atom_types} atom types",
        f"{bond_types} bond types",
        "",
        f"{format_number(xlo)} {format_number(xhi)} xlo xhi",
        f"{format_number(ylo)} {format_number(yhi)} ylo yhi",
        f"{format_number(zlo)} {format_number(zhi)} zlo zhi",
        "",
        atoms_header,
        "",
    ]
    for atom in atoms:
        fields = [
            str(atom.atom_id),
            str(atom.mol_id),
            str(atom.atom_type),
            format_number(atom.x),
            format_number(atom.y),
            format_number(atom.z),
        ]
        if include_images:
            fields.extend(
                [
                    str(0 if atom.ix is None else atom.ix),
                    str(0 if atom.iy is None else atom.iy),
                    str(0 if atom.iz is None else atom.iz),
                ]
            )
        line = " ".join(fields)
        if atom.comment:
            line += f" # {atom.comment}"
        lines.append(line)
    lines.extend(["", "Bonds", ""])
    for bond in bonds:
        lines.append(f"{bond.bond_id} {bond.bond_type} {bond.atom1} {bond.atom2}")
    lines.append("")
    return "\n".join(lines)
