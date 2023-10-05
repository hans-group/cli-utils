#!/usr/bin/env python
import argparse
import ast
from pymatgen.core.structure import Structure
from pymatgen.core import surface


def parse_literal(s):
    try:
        return ast.literal_eval(s)
    except (SyntaxError, ValueError):
        raise argparse.ArgumentTypeError(f"Invalid tuple: {s}")


def parse_cli_args():
    parser = argparse.ArgumentParser(
        description="Generate a slab from bulk structure. Requires a POSCAR file in current directory."
    )
    parser.add_argument(
        "-i",
        "--index",
        default=[0, 0, 1],
        dest="miller_index",
        type=parse_literal,
        help="The miller index of chosen plane default to [0,0,1].",
    )
    parser.add_argument(
        "-d",
        "--depth",
        default=10.0,
        dest="slab_depth",
        type=float,
        help="The depth of slab, default to 10 angstrom.",
    )
    parser.add_argument(
        "-v",
        "--vaccum",
        default=20.0,
        dest="minimum_vaccum_size",
        type=float,
        help="The vaccum size of surface. Default to 20 angstrom.",
    )
    parser.add_argument(
        "-s",
        "--shift",
        default=0.0,
        dest="shift_factor",
        type=float,
        help="The shift factor of slab. (0 to 1). Default to 0.",
    )
    parser.add_argument(
        "-S",
        "--size",
        default=(1, 1, 1),
        dest="supercell_size",
        type=parse_literal,
        help="The superecell size. Default to 1x1x1 (Enter 1,1,1).",
    )

    parser.add_argument(
        "-p",
        "--poscar",
        default="./POSCAR",
        type=str,
        help="Path to bulk POSCAR file. Defaults to POSCAR in current directory",
    )

    args = parser.parse_args()

    # Sanity check
    if any(not isinstance(i, int) for i in args.miller_index):
        raise ValueError("Miller indices should be integers")

    if args.slab_depth <= 0:
        raise ValueError("depth should be greater than 0")

    if args.minimum_vaccum_size <= 0:
        raise ValueError("minimum_vaccum_size should be greater than 0")

    if args.shift_factor > 1 or args.shift_factor < 0:
        raise ValueError("Shift factor must lie in 0 ~ 1")

    if any(not isinstance(i, int) for i in args.size) or any(i < 1 for i in args.size):
        raise ValueError("Supercell size should be positive integers")
    return args


if __name__ == "__main__":
    args = parse_cli_args()
    poscar = Structure.from_file(f"{args.poscar}")

    slab_gen = surface.SlabGenerator(
        initial_structure=poscar,
        miller_index=args.miller_index,
        center_slab=False,
        min_slab_size=args.slab_depth,
        min_vacuum_size=args.minimum_vaccum_size,
    )
    slab_pos = slab_gen.get_slab(shift=args.shift_factor) * args.supercell_size
    # shift  : relative value (0 to 1)

    slab_pos.to(filename=f"POSCAR_{''.join(map(str, args.miller_index))}", fmt="poscar")
