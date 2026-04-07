from __future__ import annotations

from dataclasses import dataclass


@dataclass(slots=True)
class Node:
    x: float
    y: float
    z: float
    extras: tuple[str, ...] = ()


@dataclass(slots=True)
class Chain:
    nodes: list[Node]


@dataclass(slots=True)
class Frame:
    box: tuple[float, float, float]
    chains: list[Chain]


@dataclass(slots=True)
class LammpsAtom:
    atom_id: int
    mol_id: int
    atom_type: int
    x: float
    y: float
    z: float
    ix: int | None = None
    iy: int | None = None
    iz: int | None = None
    comment: str | None = None


@dataclass(slots=True)
class LammpsBond:
    bond_id: int
    bond_type: int
    atom1: int
    atom2: int
