from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re

_SECTION_NAMES = [
    'PairIJ Coeffs',
    'BondBond13 Coeffs',
    'MiddleBondTorsion Coeffs',
    'EndBondTorsion Coeffs',
    'AngleAngleTorsion Coeffs',
    'AngleTorsion Coeffs',
    'BondBond Coeffs',
    'BondAngle Coeffs',
    'Improper Coeffs',
    'Dihedral Coeffs',
    'Angle Coeffs',
    'Bond Coeffs',
    'Pair Coeffs',
    'Masses',
    'Atoms',
    'Velocities',
    'Ellipsoids',
    'Lines',
    'Triangles',
    'Bodies',
    'Bonds',
    'Angles',
    'Dihedrals',
    'Impropers',
]
_SECTION_NAMES = sorted(_SECTION_NAMES, key=len, reverse=True)

_COUNT_PATTERNS = {
    'atoms': re.compile(r'^(\d+)\s+atoms\b'),
    'bonds': re.compile(r'^(\d+)\s+bonds\b'),
}


@dataclass(slots=True)
class LammpsBox:
    xlo: float
    xhi: float
    ylo: float
    yhi: float
    zlo: float
    zhi: float
    xy: float = 0.0
    xz: float = 0.0
    yz: float = 0.0


@dataclass(slots=True)
class LammpsDataAtom:
    atom_id: int
    mol_id: int
    atom_type: int
    x: float
    y: float
    z: float
    charge: float | None = None
    ix: int | None = None
    iy: int | None = None
    iz: int | None = None
    style: str = 'full'
    comment: str | None = None

    def render(self) -> str:
        fields = [str(self.atom_id), str(self.mol_id), str(self.atom_type)]
        if self.charge is not None:
            fields.append(_format_number(self.charge))
        fields.extend([
            _format_number(self.x),
            _format_number(self.y),
            _format_number(self.z),
        ])
        if self.ix is not None and self.iy is not None and self.iz is not None:
            fields.extend([str(self.ix), str(self.iy), str(self.iz)])
        line = ' '.join(fields)
        if self.comment:
            line += f' # {self.comment}'
        return line


@dataclass(slots=True)
class LammpsDataBond:
    bond_id: int
    bond_type: int
    atom1: int
    atom2: int
    comment: str | None = None


@dataclass(slots=True)
class LammpsDataSection:
    name: str
    header_line: str
    body_lines: list[str]


@dataclass(slots=True)
class LammpsDataFile:
    path: Path
    title: str
    preamble_lines: list[str]
    sections: list[LammpsDataSection]
    box: LammpsBox
    declared_atoms: int | None
    declared_bonds: int | None
    masses: dict[int, float]
    atoms: list[LammpsDataAtom]
    bonds: list[LammpsDataBond]

    def atom_map(self) -> dict[int, LammpsDataAtom]:
        return {atom.atom_id: atom for atom in self.atoms}

    def render(self, atoms: list[LammpsDataAtom] | None = None) -> str:
        atom_records = sorted(self.atoms if atoms is None else atoms, key=lambda atom: atom.atom_id)
        lines = [self.title, *self.preamble_lines]
        for section in self.sections:
            lines.append(section.header_line)
            if section.name == 'Atoms':
                lines.append('')
                lines.extend(atom.render() for atom in atom_records)
                lines.append('')
            else:
                lines.extend(section.body_lines)
        if lines and lines[-1] != '':
            lines.append('')
        return '\n'.join(lines)


@dataclass(slots=True)
class LammpsDumpFrame:
    timestep: str
    atom_count: int
    box_header: str
    box_lines: tuple[str, str, str]
    atom_fields: tuple[str, ...]
    atoms: list[dict[str, str]]


def _format_number(value: float) -> str:
    return f'{value:.15g}'


def _split_comment(line: str) -> tuple[str, str | None]:
    if '#' not in line:
        return line.rstrip('\n'), None
    content, comment = line.split('#', 1)
    return content.rstrip('\n'), comment.strip() or None


def _strip_content(line: str) -> str:
    content, _ = _split_comment(line)
    return ' '.join(content.strip().split())


def _section_name_from_header(line: str) -> str | None:
    stripped = _strip_content(line)
    for name in _SECTION_NAMES:
        if stripped == name or stripped.startswith(f'{name} '):
            return name
    return None


def parse_lammps_data(path: str | Path) -> LammpsDataFile:
    path = Path(path)
    lines = path.read_text().splitlines()
    if not lines:
        raise ValueError(f'{path} is empty')

    title = lines[0]
    section_indices: list[tuple[int, str]] = []
    for index, line in enumerate(lines[1:], start=1):
        name = _section_name_from_header(line)
        if name is not None:
            section_indices.append((index, name))

    first_section = section_indices[0][0] if section_indices else len(lines)
    preamble_lines = lines[1:first_section]
    declared_atoms = _parse_declared_count(preamble_lines, 'atoms')
    declared_bonds = _parse_declared_count(preamble_lines, 'bonds')
    box = _parse_box(preamble_lines)

    sections: list[LammpsDataSection] = []
    for offset, (start, name) in enumerate(section_indices):
        end = section_indices[offset + 1][0] if offset + 1 < len(section_indices) else len(lines)
        sections.append(
            LammpsDataSection(
                name=name,
                header_line=lines[start],
                body_lines=lines[start + 1:end],
            )
        )

    section_map = {section.name: section for section in sections}
    masses = _parse_masses(section_map.get('Masses'))
    atoms = _parse_atoms(section_map.get('Atoms'))
    bonds = _parse_bonds(section_map.get('Bonds'))

    if declared_atoms is not None and declared_atoms != len(atoms):
        raise ValueError(f'{path} declares {declared_atoms} atoms but contains {len(atoms)} atom records')
    if declared_bonds is not None and declared_bonds != len(bonds):
        raise ValueError(f'{path} declares {declared_bonds} bonds but contains {len(bonds)} bond records')

    return LammpsDataFile(
        path=path,
        title=title,
        preamble_lines=preamble_lines,
        sections=sections,
        box=box,
        declared_atoms=declared_atoms,
        declared_bonds=declared_bonds,
        masses=masses,
        atoms=atoms,
        bonds=bonds,
    )


def _parse_declared_count(lines: list[str], key: str) -> int | None:
    pattern = _COUNT_PATTERNS[key]
    for line in lines:
        stripped = _strip_content(line)
        match = pattern.match(stripped)
        if match:
            return int(match.group(1))
    return None


def _parse_box(lines: list[str]) -> LammpsBox:
    bounds: dict[str, tuple[float, float]] = {}
    tilt = (0.0, 0.0, 0.0)
    for line in lines:
        stripped = _strip_content(line)
        if not stripped:
            continue
        tokens = stripped.split()
        if len(tokens) >= 4 and tokens[2:] == ['xlo', 'xhi']:
            bounds['x'] = (float(tokens[0]), float(tokens[1]))
        elif len(tokens) >= 4 and tokens[2:] == ['ylo', 'yhi']:
            bounds['y'] = (float(tokens[0]), float(tokens[1]))
        elif len(tokens) >= 4 and tokens[2:] == ['zlo', 'zhi']:
            bounds['z'] = (float(tokens[0]), float(tokens[1]))
        elif len(tokens) >= 6 and tokens[3:] == ['xy', 'xz', 'yz']:
            tilt = (float(tokens[0]), float(tokens[1]), float(tokens[2]))
    if set(bounds) != {'x', 'y', 'z'}:
        raise ValueError('LAMMPS data file is missing one or more box-bound lines')
    return LammpsBox(
        xlo=bounds['x'][0],
        xhi=bounds['x'][1],
        ylo=bounds['y'][0],
        yhi=bounds['y'][1],
        zlo=bounds['z'][0],
        zhi=bounds['z'][1],
        xy=tilt[0],
        xz=tilt[1],
        yz=tilt[2],
    )


def _iter_section_records(section: LammpsDataSection | None) -> list[tuple[str, str | None]]:
    if section is None:
        return []
    records: list[tuple[str, str | None]] = []
    for line in section.body_lines:
        stripped = _strip_content(line)
        if not stripped:
            continue
        content, comment = _split_comment(line)
        content = ' '.join(content.strip().split())
        if not content:
            continue
        records.append((content, comment))
    return records


def _parse_masses(section: LammpsDataSection | None) -> dict[int, float]:
    masses: dict[int, float] = {}
    for content, _comment in _iter_section_records(section):
        tokens = content.split()
        if len(tokens) < 2:
            continue
        masses[int(tokens[0])] = float(tokens[1])
    return masses


def _parse_atoms(section: LammpsDataSection | None) -> list[LammpsDataAtom]:
    atoms: list[LammpsDataAtom] = []
    for content, comment in _iter_section_records(section):
        tokens = content.split()
        if len(tokens) == 6:
            atom_id, mol_id, atom_type = map(int, tokens[:3])
            x, y, z = map(float, tokens[3:6])
            atoms.append(LammpsDataAtom(atom_id=atom_id, mol_id=mol_id, atom_type=atom_type, x=x, y=y, z=z, style='angle', comment=comment))
        elif len(tokens) == 7:
            atom_id, mol_id, atom_type = map(int, tokens[:3])
            charge = float(tokens[3])
            x, y, z = map(float, tokens[4:7])
            atoms.append(LammpsDataAtom(atom_id=atom_id, mol_id=mol_id, atom_type=atom_type, charge=charge, x=x, y=y, z=z, style='full', comment=comment))
        elif len(tokens) == 9:
            atom_id, mol_id, atom_type = map(int, tokens[:3])
            x, y, z = map(float, tokens[3:6])
            ix, iy, iz = map(int, tokens[6:9])
            atoms.append(LammpsDataAtom(atom_id=atom_id, mol_id=mol_id, atom_type=atom_type, x=x, y=y, z=z, ix=ix, iy=iy, iz=iz, style='angle', comment=comment))
        elif len(tokens) == 10:
            atom_id, mol_id, atom_type = map(int, tokens[:3])
            charge = float(tokens[3])
            x, y, z = map(float, tokens[4:7])
            ix, iy, iz = map(int, tokens[7:10])
            atoms.append(LammpsDataAtom(atom_id=atom_id, mol_id=mol_id, atom_type=atom_type, charge=charge, x=x, y=y, z=z, ix=ix, iy=iy, iz=iz, style='full', comment=comment))
        else:
            raise ValueError(f'Unsupported Atoms record with {len(tokens)} columns: {content!r}')
    return atoms


def _parse_bonds(section: LammpsDataSection | None) -> list[LammpsDataBond]:
    bonds: list[LammpsDataBond] = []
    for content, comment in _iter_section_records(section):
        tokens = content.split()
        if len(tokens) < 4:
            raise ValueError(f'Invalid bond record: {content!r}')
        bonds.append(
            LammpsDataBond(
                bond_id=int(tokens[0]),
                bond_type=int(tokens[1]),
                atom1=int(tokens[2]),
                atom2=int(tokens[3]),
                comment=comment,
            )
        )
    return bonds


def iter_lammps_dump_frames(path: str | Path):
    path = Path(path)
    with path.open() as handle:
        while True:
            line = handle.readline()
            if line == '':
                break
            if line.strip() == '':
                continue
            if line.strip() != 'ITEM: TIMESTEP':
                raise ValueError(f'{path} is not a supported LAMMPS dump file: expected ITEM: TIMESTEP, found {line.strip()!r}')
            timestep = handle.readline().strip()
            if handle.readline().strip() != 'ITEM: NUMBER OF ATOMS':
                raise ValueError(f'{path} is missing ITEM: NUMBER OF ATOMS after timestep {timestep}')
            atom_count = int(handle.readline().strip())
            box_header = handle.readline().rstrip('\n')
            if not box_header.startswith('ITEM: BOX BOUNDS'):
                raise ValueError(f'{path} is missing ITEM: BOX BOUNDS after timestep {timestep}')
            box_lines = tuple(handle.readline().rstrip('\n') for _ in range(3))
            atoms_header = handle.readline().rstrip('\n')
            if not atoms_header.startswith('ITEM: ATOMS '):
                raise ValueError(f'{path} is missing ITEM: ATOMS after timestep {timestep}')
            atom_fields = tuple(atoms_header.split()[2:])
            atoms: list[dict[str, str]] = []
            for _ in range(atom_count):
                row = handle.readline()
                if row == '':
                    raise ValueError(f'{path} ended mid-frame at timestep {timestep}')
                values = row.strip().split()
                if len(values) != len(atom_fields):
                    raise ValueError(
                        f'{path} has {len(values)} values for timestep {timestep}, expected {len(atom_fields)} from {atoms_header!r}'
                    )
                atoms.append(dict(zip(atom_fields, values)))
            yield LammpsDumpFrame(
                timestep=timestep,
                atom_count=atom_count,
                box_header=box_header,
                box_lines=box_lines,
                atom_fields=atom_fields,
                atoms=atoms,
            )


def choose_coordinate_fields(atom_fields: tuple[str, ...]) -> tuple[str, str, str]:
    for triplet in (('x', 'y', 'z'), ('xu', 'yu', 'zu'), ('xs', 'ys', 'zs')):
        if all(field in atom_fields for field in triplet):
            return triplet
    raise ValueError(f'Dump frame does not contain x/y/z, xu/yu/zu, or xs/ys/zs columns: {atom_fields!r}')
