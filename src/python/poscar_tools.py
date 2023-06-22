#!/usr/bin/env python
from pathlib import Path

import ase.io
import click
import ase.build
import ase.constraints
import numpy as np


# Helper functions
def detect_format(filepath: str) -> str:
    filepath = Path(filepath)
    stem = filepath.stem
    suffix = filepath.suffix

    if suffix == ".vasp":
        return "vasp"
    elif any(s in stem for s in ["POSCAR", "CONTCAR", "vasprun.xml", "XDATCAR"]):
        return "vasp"
    else:
        return suffix[1:]


# Make a cli app using click, with possible subcommands
@click.group()
def app():
    pass


# Subcommand to convert a POSCAR to the other formats
@app.command(help="Convert the format of structure file")
@click.argument("input_file", type=click.Path(exists=True))
@click.argument("output_file", type=click.Path())
@click.option("-i", "--index", type=str, default=":")
@click.option("--from", "from_", type=str, default=None)
@click.option("--to", "to_", type=str, default=None)
def convert_format(input_file, output_file, index, from_, to_):
    atoms = ase.io.read(input_file, index, format=from_)
    ase.io.write(output_file, atoms, format=to_)


# Subcommand to make supercell
@app.command(help="Make a supercell")
@click.argument("input_file", type=click.Path(exists=True), default="POSCAR")
@click.option("--size", "-s", type=int, nargs=3, default=[1, 1, 1])
@click.option("--sort", "-S", is_flag=True, default=True)
@click.option("--output", "-o", type=click.Path(), default=None)
def make_supercell(input_file, size, sort, output):
    atoms = ase.io.read(input_file)
    M = np.array(
        [
            [size[0], 0, 0],
            [0, size[1], 0],
            [0, 0, size[2]],
        ]
    )
    scel = ase.build.make_supercell(atoms, M)
    if sort:
        scel = ase.build.sort(scel)
    if output is None:
        input_file = Path(input_file)
        output = input_file.stem + "_scel" + input_file.suffix
    ase.io.write(output, scel, format=detect_format(output))


# Subcommand to fix the atoms below a certain height
@app.command(help="Fix the atoms in POSCAR below a certain height")
@click.argument("input_file", type=click.Path(exists=True), default="POSCAR")
@click.option("--height", "-h", type=float)
@click.option("--output", "-o", type=click.Path(), default="POSCAR_fixed")
def fix_atoms(input_file, height, output):
    atoms = ase.io.read(input_file)
    atoms.set_constraint(ase.constraints.FixAtoms(mask=atoms.positions[:, 2] < height))
    ase.io.write(output, atoms, format="vasp")


# Subcommand to sort the atoms by specified order
# option --order should be a list of chemical symbols, e.g. --order "H" "O" "C" "N"
@app.command(help="Sort the atoms by specified order")
@click.argument("input_file", type=click.Path(exists=True), default="POSCAR")
@click.option("--order", "-o", type=str, default=None)
def sort_atoms(input_file, order):
    atoms = ase.io.read(input_file)
    if order is None:
        sort_keys = None
    else:
        symbols = order.split(",")
        index_map = {symbol: i for i, symbol in enumerate(symbols)}
        sort_keys = [
            index_map.get(symbol, len(order)) for symbol in atoms.get_chemical_symbols()
        ]
    atoms = ase.build.sort(atoms, tags=sort_keys)
    ase.io.write(input_file, atoms, format=detect_format(input_file))


if __name__ == "__main__":
    app()
