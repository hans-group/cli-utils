import io
import os
import shutil
import subprocess
import tempfile

from ase import Atoms
from rdkit import Chem


def _check_xtb_install():
    try:
        subprocess.check_output(["xtb", "--version"])
    except FileNotFoundError:
        raise FileNotFoundError("xtb is not installed. Please install it first.")


def get_coord_str(mol_or_atoms, output_format="xyz") -> str:
    if output_format == "xyz":
        return get_xyz_str(mol_or_atoms)
    elif output_format == "pdb":
        return get_pdb_str(mol_or_atoms)
    else:
        raise ValueError("Only xyz and pdb are supported")


def get_xyz_str(mol_or_atoms):
    if isinstance(mol_or_atoms, Chem.Mol):
        return Chem.MolToXYZBlock(mol_or_atoms)
    elif isinstance(mol_or_atoms, Atoms):
        f = io.StringIO()
        mol_or_atoms.write(f, format="xyz")
        f.seek(0)
        return f.read()


def get_pdb_str(mol_or_atoms):
    if isinstance(mol_or_atoms, Chem.Mol):
        return Chem.MolToPDBBlock(mol_or_atoms)
    elif isinstance(mol_or_atoms, Atoms):
        f = io.StringIO()
        mol_or_atoms.write(f, format="proteindatabank")
        f.seek(0)
        return f.read()


def xtb_optimize(coord_str, output_filename, coord_format="xyz"):
    null_dev = open(os.devnull, "w")
    _check_xtb_install()
    pwd = os.getcwd()
    with tempfile.TemporaryDirectory() as tmpdir:
        with open(f"{tmpdir}/input.{coord_format}", "w") as f:
            f.write(coord_str)
        os.chdir(tmpdir)
        subprocess.run(
            ["xtb", f"{tmpdir}/input.{coord_format}", "--opt", "--gfn", "2"],
            check=True,
            stdout=null_dev,
            stderr=null_dev,
        )
        os.chdir(pwd)
        output_geom_file = f"{tmpdir}/xtbopt.{coord_format}"
        shutil.copy(output_geom_file, output_filename)
