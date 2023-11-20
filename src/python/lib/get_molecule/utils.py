from typing import Optional, Tuple, Union

from rdkit import Chem
from rdkit.Chem import AllChem

from .console import console
from . import pubchempy as pcp


def show_compound_info(compound: pcp.Compound):
    print()
    console.print("[gray]Compound name[/gray]:", compound.synonyms[0])
    console.print("[gray]CID[/gray]:", compound.cid)
    console.print("[gray]Canonical SMILES[/gray]:", compound.canonical_smiles)
    console.print("[gray]Isomeric SMILES[/gray]:", compound.isomeric_smiles)
    print()


def extract_field(
    cid: Optional[int] = None, name: Optional[str] = None, smiles: Optional[str] = None
) -> Tuple[str, Union[int, str]]:
    if (bool(cid) + bool(name) + bool(smiles)) != 1:
        raise ValueError("Must specify exactly one of cid, name, or smiles")

    fields = {"cid": cid, "name": name, "smiles": smiles}
    namespace, value = next((k, v) for k, v in fields.items() if v is not None)
    return namespace, value


def determine_output(namespace, value, filename=None, output_format=None):
    if filename is None and output_format is None:
        raise ValueError("Must specify either filename or output_format")
    if output_format is None:
        output_format = filename.split(".")[-1]
    if filename is None:
        if namespace == "cid":
            filename = f"cid_{value}.{output_format}"
        if namespace == "name":
            name = value.lower().replace(" ", "_")
            filename = f"{name}.{output_format}"
        else:
            filename = f"{value}.{output_format}"
    # override output format if filename is specified
    output_format = filename.split(".")[-1]
    return filename, output_format


def insert_index(filename, idx):
    *parts, suffix = filename.split(".")
    prefix = ".".join(parts)
    return f"{prefix}_{idx}.{suffix}"


def pubchem_get_compound(namespace, value):
    console.print(f"Searching for {namespace}={value}...", style="info")
    compound = pcp.get_compounds(value, namespace, record_type="3d")[0]
    # Search again with CID because the first search doesn't contain most
    # of the information we want
    cid = compound.cid
    compound = pcp.Compound.from_cid(cid)
    return compound


def mmff_optimize(mol):
    AllChem.MMFFOptimizeMolecule(mol)
    return mol


def mol_from_smiles(smiles: str, optimize: bool = True):
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
    if optimize:
        mol = mmff_optimize(mol)
    console.print("Generated conformer with rdkit", style="info")
    console.print("[warning]Note[/warning]: this may not be correct.")
    return mol


def read_mol(filename):
    if filename.endswith(".sdf"):
        return Chem.SDMolSupplier(filename)[0]
    elif filename.endswith(".pdb"):
        return Chem.MolFromPDBFile(filename)
    else:
        raise ValueError("File format not supported for RDKit.")


def write_mol(mol, filename, output_format):
    if output_format == "sdf":
        with Chem.SDWriter(filename) as w:
            w.write(mol)
    elif output_format == "xyz":
        Chem.MolToXYZFile(mol, filename)
    elif output_format == "pdb":
        Chem.MolToPDBFile(mol, filename)
