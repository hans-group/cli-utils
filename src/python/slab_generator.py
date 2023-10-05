#!/usr/bin/env python
import argparse
import ast
from pymatgen.core.structure import Structure
from pymatgen.core import surface

def parse_tuple(s):
    try:
        return ast.literal_eval(s)
    except (SyntaxError, ValueError):
        raise argparse.ArgumentTypeError(f"Invalid tuple: {s}")


parser = argparse.ArgumentParser(description='Generate a slab from bulk structure. Requires a POSCAR file in current directory.')
parser.add_argument("-i", "--index", default=[0,0,1], dest="miller_index", type=parse_tuple, help="The miller index of chosen plane default to [0,0,1].")
parser.add_argument("-d", "--depth",  default=10, dest="slab_depth", type=parse_tuple, help="The depth of slab, default to 10 angstrom.")
parser.add_argument("-v", "--vaccum", default=20, dest="minimum_vaccum_size", type=parse_tuple, help="The vaccum size of surface. Default to 20 angstrom.")
parser.add_argument("-s", "--shift", default=0, dest="shift_factor", type=parse_tuple, help="The shift factor of slab. (0 to 1). Default to 0.")
parser.add_argument("-S", "--size", default=(1,1,1), dest="supercell_size", type=parse_tuple, help="The superecell size. Default to 1x1x1 (Enter 1,1,1).")

args = parser.parse_args()



poscar_path = './POSCAR'
pos = Structure.from_file(f'{poscar_path}')

slab_pos = surface.SlabGenerator(initial_structure=pos, miller_index=args.miller_index, center_slab=False, min_slab_size=args.slab_depth,
                                 min_vacuum_size=args.minimum_vaccum_size ).get_slab(shift=args.shift_factor) * args.supercell_size  #shift  : relative value (0 to 1)


slab_pos.to(filename=f"POSCAR_{''.join(map(str, args.miller_index))}", fmt='poscar')

