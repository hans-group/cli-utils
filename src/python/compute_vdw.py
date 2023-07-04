#!/usr/bin/env python
import ase.io
import click
from dftd4.ase import DFTD4



@click.command()
@click.argument("filename", type=click.Path())
@click.option("-m", "--method", default="PBE")
def main(filename, method):
    atoms = ase.io.read(filename)
    atoms.calc = DFTD4(method=method)
    vdw_corr = atoms.get_potential_energy()
    print(vdw_corr)

if __name__ == "__main__":
    main()

