from __future__ import annotations

import math
from pathlib import Path
from typing import TextIO

from .formats import Chain, Frame, Node


def strip_line(line: str) -> str:
    return ' '.join(line.strip().split())


def round_away_from_zero(value: float) -> int:
    if value == 0:
        return 0
    return int(math.copysign(math.floor(abs(value) + 0.5), value))


def _next_content_line(handle: TextIO) -> str | None:
    while True:
        line = handle.readline()
        if line == '':
            return None
        line = strip_line(line)
        if line:
            return line


def _parse_lengths_line(line: str, chains: int) -> list[int]:
    if '*' in line:
        parts = [part.strip() for part in line.split('*', 1)]
        if len(parts) != 2:
            raise ValueError(f'Invalid compact chain-length line: {line!r}')
        repeat = int(parts[0])
        length = int(parts[1])
        if repeat != chains:
            raise ValueError(
                f'Compact chain-length line says {repeat} chains but frame declares {chains}'
            )
        return [length] * repeat
    lengths = [int(token) for token in line.split()]
    if len(lengths) != chains:
        raise ValueError(
            f'Expected {chains} chain lengths, found {len(lengths)} in line: {line!r}'
        )
    return lengths


def parse_z1_frames(path: str | Path) -> list[Frame]:
    path = Path(path)
    frames: list[Frame] = []
    with path.open() as handle:
        while True:
            chain_line = _next_content_line(handle)
            if chain_line is None:
                break
            chains = int(chain_line)
            box_line = _next_content_line(handle)
            if box_line is None:
                raise ValueError(f'{path} ended before the box line')
            box_tokens = box_line.split()
            if len(box_tokens) < 3:
                raise ValueError(f'{path} has an invalid box line: {box_line!r}')
            box = (float(box_tokens[0]), float(box_tokens[1]), float(box_tokens[2]))
            lengths_line = _next_content_line(handle)
            if lengths_line is None:
                raise ValueError(f'{path} ended before the chain-length line')
            lengths = _parse_lengths_line(lengths_line, chains)
            frame_chains: list[Chain] = []
            for chain_length in lengths:
                nodes: list[Node] = []
                for _ in range(chain_length):
                    coord_line = _next_content_line(handle)
                    if coord_line is None:
                        raise ValueError(f'{path} ended mid-frame')
                    tokens = coord_line.split()
                    if len(tokens) < 3:
                        raise ValueError(f'Invalid coordinate line in {path}: {coord_line!r}')
                    nodes.append(
                        Node(
                            x=float(tokens[0]),
                            y=float(tokens[1]),
                            z=float(tokens[2]),
                            extras=tuple(tokens[3:]),
                        )
                    )
                frame_chains.append(Chain(nodes))
            frames.append(Frame(box=box, chains=frame_chains))
    return frames


def parse_dat_frames(path: str | Path) -> list[Frame]:
    path = Path(path)
    frames: list[Frame] = []
    with path.open() as handle:
        while True:
            chain_line = _next_content_line(handle)
            if chain_line is None:
                break
            chains = int(chain_line)
            box_line = _next_content_line(handle)
            if box_line is None:
                raise ValueError(f'{path} ended before the box line')
            box_tokens = box_line.split()
            if len(box_tokens) < 3:
                raise ValueError(f'{path} has an invalid box line: {box_line!r}')
            box = (float(box_tokens[0]), float(box_tokens[1]), float(box_tokens[2]))
            frame_chains: list[Chain] = []
            for _ in range(chains):
                length_line = _next_content_line(handle)
                if length_line is None:
                    raise ValueError(f'{path} ended before a chain-length line')
                chain_length = int(length_line)
                nodes: list[Node] = []
                for _ in range(chain_length):
                    coord_line = _next_content_line(handle)
                    if coord_line is None:
                        raise ValueError(f'{path} ended mid-frame')
                    tokens = coord_line.split()
                    if len(tokens) < 3:
                        raise ValueError(f'Invalid coordinate line in {path}: {coord_line!r}')
                    nodes.append(
                        Node(
                            x=float(tokens[0]),
                            y=float(tokens[1]),
                            z=float(tokens[2]),
                            extras=tuple(tokens[3:]),
                        )
                    )
                frame_chains.append(Chain(nodes))
            frames.append(Frame(box=box, chains=frame_chains))
    return frames
