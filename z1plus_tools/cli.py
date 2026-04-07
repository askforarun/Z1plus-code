from __future__ import annotations

import argparse
import sys

from .backbone import BackboneError, extract_backbone
from .import_lammps import ImportError, import_lammps_to_z1
from .postprocess import (
    convert_dat_file_to_dump_text,
    convert_z1_file_to_dump_text,
    create_single_chain_entanglement_outputs,
    create_sp_to_data_outputs,
    export_merged_outputs,
)
from .preprocess import correct_vmd_lammps_files


def _comma_ints(value: str) -> tuple[int, ...]:
    if not value:
        return ()
    return tuple(int(token) for token in value.split(',') if token)


def main_z1_dump(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description='Convert a Z1-formatted file to LAMMPS dump format.')
    parser.add_argument('-unfolded', action='store_true', help='Write xu/yu/zu instead of wrapped x/y/z.')
    parser.add_argument('z1_file', help='Z1-formatted configuration or trajectory file.')
    args = parser.parse_args(argv)
    sys.stdout.write(convert_z1_file_to_dump_text(args.z1_file, unfolded=args.unfolded))
    return 0


def main_z1_dat2dump(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description='Convert a Z1+ .dat file to LAMMPS dump format.')
    parser.add_argument('-unfolded', action='store_true', help='Write xu/yu/zu instead of wrapped x/y/z.')
    parser.add_argument('dat_file', help='Z1+initconfig.dat, Z1+SP.dat, or Z1+PPA.dat.')
    args = parser.parse_args(argv)
    sys.stdout.write(convert_dat_file_to_dump_text(args.dat_file, unfolded=args.unfolded))
    return 0


def main_sp_to_data(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description='Convert Z1+initconfig.dat and/or Z1+SP.dat to LAMMPS data format.')
    parser.add_argument('-i', action='store_true', help='Convert Z1+initconfig.dat to Z1+initconfig.data.')
    parser.add_argument('-s', action='store_true', help='Convert Z1+SP.dat to Z1+SP.data.')
    parser.add_argument('-si', action='store_true', help='Create a merged Z1+merge.data file from both inputs.')
    args = parser.parse_args(argv)
    if not (args.i or args.s or args.si):
        parser.error('At least one of -i, -s, or -si is required.')
    created = create_sp_to_data_outputs(
        export_initconfig=args.i or args.si,
        export_sp=args.s or args.si,
        export_merge=args.si,
    )
    for path in created:
        print(f'created {path}')
    return 0


def main_export(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description='Merge Z1+initconfig.dat and Z1+SP.dat into dump or data outputs.')
    parser.add_argument('-data', action='store_true', help='Write Z1+merged-<frame>.data files.')
    parser.add_argument('-dump', action='store_true', help='Write Z1+merged.dump (default if -data is not given).')
    parser.add_argument('-from', dest='from_frame', type=int, default=1, help='First frame to export (1-based).')
    parser.add_argument('-to', dest='to_frame', type=int, default=1, help='Last frame to export (1-based).')
    parser.add_argument('-all', action='store_true', help='Export every frame.')
    args = parser.parse_args(argv)
    dump = args.dump or not args.data
    created = export_merged_outputs(
        dump=dump,
        data=args.data,
        from_frame=args.from_frame,
        to_frame=args.to_frame,
        export_all=args.all,
    )
    for path in created:
        print(f'created {path}')
    return 0


def main_extract_single_chain_entanglements(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description='Extract one chain and all chains entangled with it from Z1+ outputs.')
    parser.add_argument('chain_id', type=int, help='1-based chain id to extract.')
    parser.add_argument('-folded', action='store_true', help='Write folded coordinates.')
    parser.add_argument('-txt', action='store_true', help='Write a plain text atom listing instead of a data file.')
    parser.add_argument('-dat', action='store_true', help='Write Z1+initconfig-chain=<id>.dat and Z1+SP-chain=<id>.dat.')
    parser.add_argument('-SP', dest='add_sp', action='store_true', help='Include shortest paths in the data output.')
    parser.add_argument('-ee', dest='add_ee', action='store_true', help='Add end-to-end bonds.')
    parser.add_argument('-o', dest='output_path', default=None, help='Override the default data-file name.')
    args = parser.parse_args(argv)
    created = create_single_chain_entanglement_outputs(
        chain_id=args.chain_id,
        folded=args.folded,
        txt=args.txt,
        dat=args.dat,
        add_sp=args.add_sp,
        add_ee=args.add_ee,
        output_path=args.output_path,
    )
    for path in created:
        print(f'created {path}')
    return 0


def main_convert_vmd_data_to_proper_data(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description='Repair LAMMPS data and optional dump files by rebuilding molecule ids from bond connectivity.'
    )
    parser.add_argument('-data', dest='data_file', required=True, help='LAMMPS data file to repair.')
    parser.add_argument('-dump', dest='dump_file', default=None, help='Optional LAMMPS dump file to repair.')
    parser.add_argument(
        '-angle',
        action='store_true',
        help='Accepted for CLI compatibility. Atom-section parsing is auto-detected from the record width.',
    )
    args = parser.parse_args(argv)
    summary = correct_vmd_lammps_files(data_path=args.data_file, dump_path=args.dump_file, angle=args.angle)
    print(
        f'data file has {summary.atoms} atoms, {summary.bonds} bonds, max ID {summary.max_atom_id}, '
        f'{summary.chains} chains.'
    )
    print(f'created {summary.corrected_data_path}')
    if summary.corrected_dump_path is not None:
        print(f'created {summary.corrected_dump_path} with {summary.frames} frames')
    return 0


def main_import_lammps(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description='Convert LAMMPS or XML inputs into config.Z1.')
    parser.add_argument('-dump', dest='dump_file', default=None, help='LAMMPS dump trajectory or snapshot.')
    parser.add_argument('-data', dest='data_file', default=None, help='LAMMPS data file with connectivity.')
    parser.add_argument('-from', dest='from_frame', type=int, default=1, help='First snapshot index to export.')
    parser.add_argument('-to', dest='to_frame', type=int, default=10**10, help='Last snapshot index to export.')
    parser.add_argument('-each', dest='each', type=int, default=1, help='Stride between exported snapshots.')
    parser.add_argument('-ignore_H', action='store_true', help='Ignore atoms whose type mass rounds to 1.')
    parser.add_argument('-ignore_dumbbells', action='store_true', help='Ignore chains of length 2.')
    parser.add_argument('-ignore_types', type=_comma_ints, default=(), help='Comma-separated list of atom types to ignore.')
    parser.add_argument('-Z1', dest='z1_mode', action='store_true', help='Write plain Z1 format without the trajectory footer.')
    parser.add_argument('-self-entanglements', action='store_true', help='Accepted for compatibility; currently ignored.')
    parser.add_argument('-branched', dest='allow_branched', action='store_true', help='Delete excess sidechain bonds so every node has degree <= 2.')
    parser.add_argument('-xml', dest='xml_file', default=None, help='HOOMD-style XML input.')
    parser.add_argument('-out-dump', dest='out_dump', default=None, help='Optional dump of the reconstructed true chains.')
    parser.add_argument('-verbose', action='store_true', help='Print additional progress information.')
    args = parser.parse_args(argv)
    if args.dump_file is None and args.data_file is None and args.xml_file is None:
        parser.error('at least one of -dump, -data, or -xml is required')
    try:
        summary = import_lammps_to_z1(
            dump_path=args.dump_file,
            data_path=args.data_file,
            xml_path=args.xml_file,
            from_frame=args.from_frame,
            to_frame=args.to_frame,
            each=args.each,
            ignore_h=args.ignore_H,
            ignore_dumbbells=args.ignore_dumbbells,
            ignore_types=args.ignore_types,
            z1_mode=args.z1_mode,
            allow_branched=args.allow_branched,
            out_dump_path=args.out_dump,
            verbose=args.verbose,
        )
    except ImportError as exc:
        parser.exit(2, f'ERROR: {exc}\n')
    print(f'created {summary.config_path} with {summary.frames_written} frame(s)')
    if summary.table_path is not None:
        print(f'created {summary.table_path}')
    if summary.out_dump_path is not None:
        print(f'created {summary.out_dump_path}')
    return 0


def main_extract_backbone(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description='Extract one linear backbone path per molecule and write config.Z1.')
    parser.add_argument('data_file', help='LAMMPS data file.')
    parser.add_argument('dump_file', nargs='?', default=None, help='Optional LAMMPS dump trajectory.')
    parser.add_argument('-ignore-types', dest='ignore_types', type=_comma_ints, default=(), help='Comma-separated atom types to omit entirely.')
    args = parser.parse_args(argv)
    try:
        summary = extract_backbone(data_path=args.data_file, dump_path=args.dump_file, ignore_types=args.ignore_types)
    except BackboneError as exc:
        parser.exit(2, f'ERROR: {exc}\n')
    print(f'linear backbone contains {summary.linear_atoms} atoms across {summary.molecules} molecule(s)')
    return 0
