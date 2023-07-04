#!/usr/bin/env python
import shutil

import ase.constraints
import ase.io
import click


@click.command()
@click.argument("input_file", type=click.Path(exists=True))
@click.option("-h", "--height", type=float, help="Height of the surface")
@click.option("-d", "--direct", is_flag=True, help="Use direct coordinates")
def main(input_file, height, direct):
    atoms = ase.io.read(input_file)
    z_coords = (
        atoms.get_scaled_positions()[:, 2] if direct else atoms.get_positions()[:, 2]
    )
    fix_mask = z_coords <= height
    atoms.set_constraint(ase.constraints.FixAtoms(mask=fix_mask))
    shutil.copy(input_file, input_file + ".orig")
    ase.io.write(input_file, atoms, format="vasp", direct=direct, vasp5=True)


if __name__ == "__main__":
    main()
