#!/usr/bin/env python
""" Script to combine ESP and charge density into a chgcar file
"""
import click
import numpy as np
from ase import units

from lib.formats import Chgcar, CubeData


@click.command()
@click.option("--chgden", type=click.Path(exists=True), required=True)
@click.option("--totesp", type=click.Path(exists=True), required=True)
def main(chgden, totesp):
    if ".cub" not in chgden:
        raise ValueError("chgden file must be a .cub file")
    if ".cub" not in totesp:
        raise ValueError("totesp file must be a .cub file")

    chgden = CubeData.from_file(chgden)
    totesp = CubeData.from_file(totesp)
    print("Loaded files")

    # Check if the two files contains the same atoms
    identical_cell = np.allclose(chgden.atoms.cell.array, totesp.atoms.cell.array)
    identical_species = np.allclose(chgden.atoms.numbers, totesp.atoms.numbers)
    identical_coords = np.allclose(chgden.atoms.positions, totesp.atoms.positions)
    if not (identical_cell and identical_species and identical_coords):
        raise ValueError("chgden and totesp files do not contain the same atoms")

    # Map the ESP to the atoms

    vol = chgden.atoms.get_volume() / (units.Bohr**3)
    mapped_esp = Chgcar(
        chgden.atoms.cell.array,
        chgden.atoms.symbols,
        chgden.atoms.positions,
        [chgden.data * vol, totesp.data * vol * units.Ha],
    )
    mapped_esp.write("chgden_with_esp.vasp")
    print("Wrote charge density and  ESP to chgden_with_esp.vasp")


if __name__ == "__main__":
    main()
