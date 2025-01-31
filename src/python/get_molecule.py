#!/usr/bin/env python
from pathlib import Path
from typing import Optional

import ase.build
import click
from ase.data.pubchem import pubchem_atoms_conformer_search
from lib.get_molecule import utils
from lib.get_molecule.console import console
from lib.get_molecule.xtb_optimize import get_coord_str, xtb_optimize

available_types = click.Choice(["sdf", "xyz", "pdb"])


@click.group()
def app():
    """Get a molecule from various sources and write it to a file.

    Currently available sources:

        - ASE (by ase.build.molecule)
        - PubChem (by CID, name, or SMILES)
        - RDKit (by arbitrary SMILES)

    See help for each subcommand for more information.
    """
    pass


@app.command()
@click.argument("formula")
@click.option("-o", "--output", default=None, help="Output filename")
def from_formula(formula: str, output: Optional[str] = None):
    """Generate molecule from chemical formula (ASE)"""
    if output is None:
        output = f"{formula}.xyz"
    atoms = ase.build.molecule(formula)
    atoms.write(output)


@app.command()
@click.option("--cid", type=int, default=None, help="PubChem CID")
@click.option("--name", type=str, default=None, help="PubChem name")
@click.option("--smiles", type=str, default=None, help="SMILES string")
@click.option("-m", "--optimizemethod", default="mmff", help="Optimization method")
@click.option("--show-info-only", is_flag=True, help="Show molecule info and exit.")
@click.option("--all", "all_conformers", is_flag=True, default=False, help="Find all conformers")
@click.option("-o", "--output", default=None, help="Output filename")
@click.option("-f", "--format", "output_format", type=available_types, default="xyz", help="Output format")
def from_pubchem(
    cid: Optional[int] = None,
    name: Optional[str] = None,
    smiles: Optional[str] = None,
    optimizemethod: str = "mmff",
    show_info_only: bool = False,
    all_conformers: bool = False,
    output: Optional[str] = None,
    output_format: Optional[str] = "xyz",
):
    """Search molecule by Pubchem API.
    Exactly one of cid, name, or SMILES must be specified.
    """
    namespace, value = utils.extract_field(cid, name, smiles)
    filename, output_format = utils.determine_output(namespace, value, output, output_format)
    compound = utils.pubchem_get_compound(namespace, value)
    utils.show_compound_info(compound)
    if show_info_only:
        return

    if not all_conformers:
        output = f"{namespace}_{value}.xyz" if output is None else output
        # generate from canonical smiles
        mol = utils.mol_from_smiles(compound.canonical_smiles)
        if optimizemethod == "xtb":
            console.print("Optimizing with xtb", style="info")
            coords = get_coord_str(mol, output_format)
            # This writes file directly
            xtb_optimize(coords, filename, output_format)
            console.print("Done.", style="info")
            console.print(f"Writing molecule to {filename}", style="info")
        else:
            console.print(f"Writing molecule to {filename}", style="info")
            utils.write_mol(mol, filename, output_format)

    else:
        if output_format != "xyz":
            raise ValueError("Currently only xyz output is supported")
        console.print("Finding all optimized conformers...", style="info")
        conformers = pubchem_atoms_conformer_search(cid=compound.cid)
        console.print(f"Found {len(conformers)} conformers", style="info")
        for i, atoms in enumerate(conformers):
            atoms.set_initial_charges(None)
            filename_i = utils.insert_index(filename, i)
            if optimizemethod == "xtb":
                console.print(f"Optimizing conformer {i} with xtb", style="info")
                coords = get_coord_str(atoms, output_format)
                # This writes file directly
                xtb_optimize(coords, filename_i, output_format)
            else:
                atoms.write(filename_i)
            console.print(f"Writing conformer {i} to {filename_i}", style="info")


@app.command()
@click.argument("smiles")
@click.option("-O", "--optimize", default=True, help="Optimize geometry with force field")
@click.option("-m", "--optimizemethod", default="mmff", help="Optimization method")
@click.option("-o", "--output", type=str, default=None, help="Output filename")
@click.option(
    "-f",
    "--format",
    "output_format",
    type=click.Choice(["sdf", "xyz", "pdb"]),
    default="xyz",
    help="Output format",
)
def from_smiles(
    smiles: str,
    optimize: bool,
    optimizemethod: str,
    output: click.File,
    output_format: str,
):
    """Generate molecule from SMILES string (RDKit).
    You can optimize the molecular geometry with a force field,
    but this may be incorrect.
    """
    console.print("Generating molecule from smiles", style="info")
    mol = utils.mol_from_smiles(smiles, optimize)
    filename, output_format = utils.determine_output("smiles", smiles, output, output_format)

    if optimizemethod == "xtb":
        console.print("Optimizing with xtb", style="info")
        coords = get_coord_str(mol, output_format)
        # This writes file directly
        xtb_optimize(coords, filename, output_format)
        console.print("Done.", style="info")
        console.print(f"Writing molecule to {filename}", style="info")
    else:
        console.print(f"Writing molecule to {filename}", style="info")
        utils.write_mol(mol, filename, output_format)


@app.command()
@click.argument("inputfile", type=str)
@click.option("-m", "--method", default="mmff", help="Optimization method")
def optimize(inputfile, method):
    """Optimize molecular geometry"""
    inputfile = Path(inputfile)
    fileformat = inputfile.suffix[1:]
    if fileformat not in ("xyz", "pdb", "sdf"):
        raise ValueError(f"File format {fileformat} not supported")
    outputfile = inputfile.stem + f"_opt.{fileformat}"
    if method == "xtb":
        console.print(f"Optimizing {inputfile} with xtb", style="info")
        coords = inputfile.read_text()
        xtb_optimize(coords, outputfile, fileformat)
        console.print("Done.", style="info")
    elif method == "mmff":
        console.print(f"Optimizing {inputfile} with MMFF", style="info")
        mol = utils.read_mol(str(inputfile))
        mol = utils.mmff_optimize(mol)
        utils.write_mol(mol, str(outputfile), fileformat)
        console.print("Done.", style="info")
    else:
        raise ValueError("Invalid method. Must be 'xtb' or 'mmff'")


if __name__ == "__main__":
    app()
