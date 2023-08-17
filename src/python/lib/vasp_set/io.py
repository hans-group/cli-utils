import os
from typing import Any, Dict, Optional, Union

import toml
from ase import Atoms
from ase.formula import Formula
from ase.io import read, write
from pymatgen.core import Structure
from pymatgen.io.ase import AseAtomsAdaptor

from .inputs import Incar, Kpoints, Potcar


class InputGen:
    """
    Generate VASP inputs "POSCAR", "INCAR", "KPOINTS", "POTCAR"
    """

    def __init__(
        self,
        structure: Union[Atoms, Structure, str, os.PathLike],
        incar_params: Dict[str, Any],
        kpoints_params: Optional[Dict] = None,
        potcar_params: Optional[Dict[str, str]] = None,
    ):

        self.structure = structure
        self.incar_params = incar_params
        self.kpoints_params = kpoints_params
        self.potcar_params = potcar_params

    @staticmethod
    def from_file(
        structure: Union[Atoms, Structure, str, os.PathLike],
        filename="vasp_settings.toml",
        path: Union[str, os.PathLike] = "./",
    ):

        file_path = os.path.join(path, filename)
        with open(file_path, "r") as f:
            config = toml.load(f)
        incar_params = config["incar"]
        if "kpoints" in config:
            kpoints_params = config["kpoints"]
        else:
            kpoints_params = None
        if "potcar" in config:
            potcar_params = config["potcar"]
        else:
            potcar_params = None

        return InputGen(
            structure=structure,
            incar_params=incar_params,
            kpoints_params=kpoints_params,
            potcar_params=potcar_params,
        )

    @staticmethod
    def from_dict(
        structure: Union[Atoms, Structure, str, os.PathLike],
        incar_params: Dict,
        kpoints_params: Optional[Dict] = None,
        potcar_params: Optional[Dict[str, str]] = None,
    ):

        return InputGen(
            structure=structure,
            incar_params=incar_params,
            kpoints_params=kpoints_params,
            potcar_params=potcar_params,
        )

    def write(self, path="./", verbose=False):
        # POSCAR
        if isinstance(self.structure, (str, os.PathLike)):
            structure = read(self.structure)
        elif isinstance(self.structure, Atoms):
            structure = self.structure
            write_poscar(structure, path, verbose=verbose)
        elif isinstance(self.structure, Structure):
            structure = AseAtomsAdaptor.get_atoms(self.structure)
            write_poscar(structure, path, verbose=verbose)
        else:
            raise TypeError

        # INCAR
        incar = Incar.from_config(structure, self.incar_params)
        incar.write(path, verbose=verbose)

        # KPOINTS
        if not self.kpoints_params:
            kpoints = Kpoints.from_atoms(structure)
        else:
            kpoints = Kpoints.from_config(self.kpoints_params)
        kpoints.write(path, verbose=verbose)

        # POTCAR
        if not self.potcar_params:
            potcar = Potcar(structure)
        else:
            potcar = Potcar.from_config(structure, self.potcar_params)
        potcar.write(path, verbose=verbose)

        _write_config(
            incar=self.incar_params,
            kpoints=self.kpoints_params,
            potcar=self.potcar_params,
            path=path,
        )


def _write_config(
    incar: Dict,
    kpoints: Optional[Dict] = None,
    potcar: Optional[Dict] = None,
    path="./",
):
    file_path = os.path.join(path, "input_params.toml")
    with open(file_path, "w") as f:
        inputs = {"incar": incar}
        if kpoints:
            inputs.update({"kpoints": kpoints})
            # kpoints = {"kpoints": kpoints}
            # toml.dump(kpoints, f)
        if potcar:
            inputs.update({"potcar": potcar})
            # potcar = {"potcar": potcar}
            # toml.dump(potcar, f)
        toml.dump(inputs, f)


def sort_atoms(image: Atoms) -> Atoms:
    copied = image.copy()
    atoms = Atoms()
    atoms.cell = copied.cell
    w = copied.symbols
    formula = Formula(f"{w}")
    count = list(formula.count())
    for _, symbol in enumerate(count):
        for i in range(len(copied)):
            if copied[i].symbol == symbol:
                atoms.append(copied[i])
    return atoms


def write_poscar(image: Union[Atoms, Structure], path="./", verbose: bool = True):
    filename = os.path.join(path, "POSCAR")

    if isinstance(image, Structure):
        image = AseAtomsAdaptor.get_atoms(image)
    atoms = sort_atoms(image)
    write(filename, atoms, vasp5=True, direct=True)

    if verbose:
        print("POSCAR generated")
