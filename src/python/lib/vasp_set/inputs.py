import os
import string
from abc import ABCMeta, abstractmethod
from typing import Dict, List, Optional, Tuple, Union

import numpy as np
from ase import Atoms
from ase.formula import Formula

from .data import (
    bool_keys,
    dict_keys,
    exp_keys,
    float_keys,
    int_keys,
    list_bool_keys,
    list_float_keys,
    list_int_keys,
    potcar_path,
    reccomended_potcars,
    special_keys,
    string_keys,
)


class VaspInputBase(metaclass=ABCMeta):
    """
    Base class for VASP input files (e.g. INCAR, POTCAR, KPOINTS)
    """

    @staticmethod
    @abstractmethod
    def from_config(atoms: Optional[Atoms], config: Dict):
        """Generate Inputs from dict_params
        Args:
            config (Dict): input parameters
        """

    def write(self, path: Union[str, os.PathLike], verbose: bool):
        """Write input file
        Args:
            path (Union[str, os.PathLike]): path where file is written
            verbose (bool): print log or not
        """


class Incar(VaspInputBase):
    def __init__(self, atoms, **kwargs):
        self.float_params = {}
        self.exp_params = {}
        self.string_params = {}
        self.int_params = {}
        self.bool_params = {}
        self.list_bool_params = {}
        self.list_int_params = {}
        self.list_float_params = {}
        self.special_params = {}
        self.dict_params = {}
        for key in float_keys:
            self.float_params[key] = None
        for key in exp_keys:
            self.exp_params[key] = None
        for key in string_keys:
            self.string_params[key] = None
        for key in int_keys:
            self.int_params[key] = None
        for key in bool_keys:
            self.bool_params[key] = None
        for key in list_bool_keys:
            self.list_bool_params[key] = None
        for key in list_int_keys:
            self.list_int_params[key] = None
        for key in list_float_keys:
            self.list_float_params[key] = None
        for key in special_keys:
            self.special_params[key] = None
        for key in dict_keys:
            self.dict_params[key] = None

        self.atoms = atoms
        self.spinpol = (
            atoms.get_initial_magnetic_moments().any() or self.int_params["ispin"] == 2
        )
        self.all_symbols = atoms.get_chemical_symbols()
        self.natoms = len(atoms)
        self.formula = atoms.get_chemical_formula()
        self.symbols, symbolcount = self._count_symbols(atoms)

        self.sort, self.resort = self._make_sort(atoms)
        self.sorted_atoms = atoms[self.sort]
        self.atoms_sorted = self._sort_atoms(atoms)

        self.symbol_count = []
        for m in self.symbols:
            self.symbol_count.append([m, symbolcount[m]])

        if (
            ("ldauu" in kwargs)
            and ("ldaul" in kwargs)
            and ("ldauj" in kwargs)
            and ("ldau_luj" in kwargs)
        ):
            raise NotImplementedError(
                "You can either specify ldaul, ldauu, and ldauj OR "
                "ldau_luj. ldau_luj is not a VASP keyword. It is a "
                "dictionary that specifies L, U and J for each "
                "chemical species in the atoms object. "
                "For example for a water molecule:"
                """ldau_luj={'H':{'L':2, 'U':4.0, 'J':0.9},
                      'O':{'L':2, 'U':4.0, 'J':0.9}}"""
            )

        for key in kwargs:
            if key in self.float_params:
                self.float_params[key] = kwargs[key]
            elif key in self.exp_params:
                self.exp_params[key] = kwargs[key]
            elif key in self.string_params:
                self.string_params[key] = kwargs[key]
            elif key in self.int_params:
                self.int_params[key] = kwargs[key]
            elif key in self.bool_params:
                self.bool_params[key] = kwargs[key]
            elif key in self.list_bool_params:
                self.list_bool_params[key] = kwargs[key]
            elif key in self.list_int_params:
                self.list_int_params[key] = kwargs[key]
            elif key in self.list_float_params:
                self.list_float_params[key] = kwargs[key]
            elif key in self.special_params:
                self.special_params[key] = kwargs[key]
            elif key in self.dict_params:
                self.dict_params[key] = kwargs[key]
            elif key in self.input_params:
                self.input_params[key] = kwargs[key]
            else:
                raise TypeError("Parameter not defined: " + key)

    @staticmethod
    def from_config(atoms, config):
        incar_keys = filter(
            lambda key: not isinstance(config[key], dict), config.keys()
        )
        incar_params_collected = {}
        for key in incar_keys:
            incar_params_collected[key] = config[key]

        elements = sorted(set(atoms.symbols), key=lambda x: atoms.symbols.index(x))

        element_specific_keys = filter(
            lambda key: isinstance(config[key], dict), config.keys()
        )
        for key in element_specific_keys:
            key_list = []
            for specie in elements:
                param = config[key][str(specie)]
                key_list.append(param)
            incar_params_collected[key] = key_list

        return Incar(atoms, **incar_params_collected)

    def _sort_atoms(
        self,
        atoms: Atoms,
    ) -> Atoms:
        copied = atoms.copy()
        sorted_atoms = Atoms()
        sorted_atoms.cell = copied.cell
        w = copied.symbols
        formula = Formula(f"{w}")
        count = list(formula.count())
        for _, symbol in enumerate(count):
            for i in range(len(copied)):
                if copied[i].symbol == symbol:
                    sorted_atoms.append(copied[i])
        return sorted_atoms

    def _count_symbols(self, atoms: Atoms):
        """Count symbols in atoms object, excluding a set of indices
        Parameters:
            atoms: Atoms object to be grouped
            exclude: List of indices to be excluded from the counting
        Returns:
            Tuple of (symbols, symbolcount)
            symbols: The unique symbols in the included list
            symbolscount: Count of symbols in the included list
        Example:
        >>> from ase.build import bulk
        >>> atoms = bulk('NaCl', crystalstructure='rocksalt', a=4.1, cubic=True)
        >>> count_symbols(atoms)
        (['Na', 'Cl'], {'Na': 4, 'Cl': 4})
        >>> count_symbols(atoms, exclude=(1, 2, 3))
        (['Na', 'Cl'], {'Na': 3, 'Cl': 2})
        """
        symbols = []
        symbolcount = {}
        for _, symbol in enumerate(atoms.symbols):
            if symbol not in symbols:
                symbols.append(symbol)
                symbolcount[symbol] = 1
            else:
                symbolcount[symbol] += 1
        return symbols, symbolcount

    def _make_sort(self, atoms: Atoms) -> Tuple[List[int], List[int]]:
        symbols, _ = self._count_symbols(atoms)

        # Create sorting list
        srt = []  # type: List[int]

        for symbol in symbols:
            for m, atom in enumerate(atoms):
                if atom.symbol == symbol:
                    srt.append(m)
        # Create the resorting list
        resrt = list(range(len(srt)))
        for n in range(len(resrt)):
            resrt[srt[n]] = n
        return srt, resrt

    def write(self, path="./", verbose=True):
        # cwd = Path.cwd()
        # os.chdir(path)
        filename = os.path.join(path, "INCAR")
        magmom_written = False
        with open(filename, "w") as incar:
            incar.write(f"INCAR created for {self.formula}\n")
            for key, val in self.float_params.items():
                if val is not None:
                    incar.write(" %s = %5.6f\n" % (key.upper(), val))

            for key, val in self.exp_params.items():
                if val is not None:
                    incar.write(" %s = %5.2e\n" % (key.upper(), val))

            for key, val in self.string_params.items():
                if val is not None:
                    incar.write(" %s = %s\n" % (key.upper(), val))

            for key, val in self.int_params.items():
                if val is not None:
                    incar.write(" %s = %d\n" % (key.upper(), val))
                    if key == "ichain" and val > 0:
                        incar.write(" IBRION = 3\n POTIM = 0.0\n")
                        for key, val in self.int_params.items():
                            if key == "iopt" and val is None:
                                print(
                                    "WARNING: optimization is "
                                    "set to LFBGS (IOPT = 1)"
                                )
                                incar.write(" IOPT = 1\n")
                        for key, val in self.exp_params.items():
                            if key == "ediffg" and val is None:
                                RuntimeError("Please set EDIFFG < 0")

            for key, val in self.list_bool_params.items():
                if val is None:
                    pass
                else:
                    incar.write(" %s = " % key.upper())
                    [incar.write("%s " % self._to_vasp_bool(x)) for x in val]
                    incar.write("\n")

            for key, val in self.list_int_params.items():
                if val is None:
                    pass
                elif key == "ldaul" and (self.dict_params["ldau_luj"] is not None):
                    pass
                else:
                    incar.write(" %s = " % key.upper())
                    [incar.write("%d " % x) for x in val]
                    incar.write("\n")

            for key, val in self.list_float_params.items():
                if val is None:
                    pass
                elif (key in ("ldauu", "ldauj")) and (
                    self.dict_params["ldau_luj"] is not None
                ):
                    pass
                elif key == "magmom":
                    if not len(val) == len(self.atoms):
                        if not len(val) == len(self.symbols):
                            msg = (
                                f"Expected length of magmom tag to be {len(self.symbols)}, "
                                f"i.e. 1 value per atom, but got {len(val)}"
                            )
                            raise ValueError(msg)
                        val_new = []
                        for i, symbol in enumerate(self.symbol_count):
                            for _ in range(symbol[1]):
                                val_new.append(val[i])
                        val = val_new

                    # Check if user remembered to specify ispin
                    # note: we do not overwrite ispin if ispin=1
                    if not self.int_params["ispin"]:
                        self.spinpol = True
                        incar.write(" ispin = 2\n".upper())

                    incar.write(" %s = " % key.upper())
                    magmom_written = True
                    # Work out compact a*x b*y notation and write in this form
                    # Assume 1 magmom per atom, ordered as our atoms object
                    val = np.array(val)
                    val = val[self.sort]  # Order in VASP format

                    # Compactify the magmom list to symbol order
                    lst = [[1, val[0]]]
                    for n in range(1, len(val)):
                        if val[n] == val[n - 1]:
                            lst[-1][0] += 1
                        else:
                            lst.append([1, val[n]])
                    incar.write(
                        " ".join(["{:d}*{:.2f}".format(mom[0], mom[1]) for mom in lst])
                    )
                    incar.write("\n")
                else:
                    incar.write(" %s = " % key.upper())
                    [incar.write("%.4f " % x) for x in val]
                    incar.write("\n")

            for key, val in self.bool_params.items():
                if val is not None:
                    incar.write(" %s = " % key.upper())
                    if val:
                        incar.write(".TRUE.\n")
                    else:
                        incar.write(".FALSE.\n")
            for key, val in self.special_params.items():
                if val is not None:
                    incar.write(" %s = " % key.upper())
                    if key == "lreal":
                        if isinstance(val, str):
                            incar.write(val + "\n")
                        elif isinstance(val, bool):
                            if val:
                                incar.write(".TRUE.\n")
                            else:
                                incar.write(".FALSE.\n")

            for key, val in self.dict_params.items():
                if val is not None:
                    if key == "ldau_luj":
                        # User didn't turn on LDAU tag.
                        # Only turn on if ldau is unspecified
                        if self.bool_params["ldau"] is None:
                            self.bool_params["ldau"] = True
                            # At this point we have already parsed our bool params
                            incar.write(" LDAU = .TRUE.\n")
                        llist = ulist = jlist = ""
                        for symbol in self.symbol_count:
                            #  default: No +U
                            luj = val.get(symbol[0], {"L": -1, "U": 0.0, "J": 0.0})
                            llist += " %i" % luj["L"]
                            ulist += " %.3f" % luj["U"]
                            jlist += " %.3f" % luj["J"]
                        incar.write(" LDAUL =%s\n" % llist)
                        incar.write(" LDAUU =%s\n" % ulist)
                        incar.write(" LDAUJ =%s\n" % jlist)

            if (
                self.spinpol
                and not magmom_written
                and self.atoms.get_initial_magnetic_moments().any()
            ):
                if not self.int_params["ispin"]:
                    incar.write(" ispin = 2\n".upper())
                # Write out initial magnetic moments
                magmom = self.atoms.get_initial_magnetic_moments()[self.sort]
                # unpack magmom array if three components specified
                if magmom.ndim > 1:
                    magmom = [item for sublist in magmom for item in sublist]
                list = [[1, magmom[0]]]
                for n in range(1, len(magmom)):
                    if magmom[n] == magmom[n - 1]:
                        list[-1][0] += 1
                    else:
                        list.append([1, magmom[n]])
                incar.write(" magmom = ".upper())
                [incar.write("%i*%.4f " % (mom[0], mom[1])) for mom in list]
                incar.write("\n")
        # os.chdir(cwd)
        if verbose:
            print("INCAR generated")

    def _to_vasp_bool(self, x):
        """Convert Python boolean to string for VASP input
        In case the value was modified to a string already, appropriate strings
        will also be accepted and cast to a standard .TRUE. / .FALSE. format.
        """
        if isinstance(x, str):
            if x.lower() in (".true.", "t"):
                x = True
            elif x.lower() in (".false.", "f"):
                x = False
            else:
                raise ValueError('"%s" not recognised as VASP Boolean')
        assert isinstance(x, bool)
        if x:
            return ".TRUE."
        else:
            return ".FALSE."


class Kpoints(VaspInputBase):
    def __init__(self, mesh: list[int], method: str = "auto"):
        _available_methods = ("auto", "gamma", "monkhorst-pack")
        if method not in _available_methods:
            raise ValueError("The method is not available.")
        elif method in _available_methods:
            self.method = string.capwords(method)

        self.mesh = mesh

    @staticmethod
    def from_config(config):
        mesh = config["mesh"]
        method = config["method"]
        return Kpoints(mesh, method)

    @staticmethod
    def from_atoms(atoms: Atoms, grid: int = 45, method: str = "auto"):
        """_summary_
        Args:
            cell (array): (3,3) array
            grid (int, optional): _description_. Defaults to 30.
        """
        cell = atoms.cell
        lengths = np.linalg.norm(cell, axis=1)
        kpoints = [0, 0, 0]
        for i in range(3):
            kpoints[i] = int(grid // lengths[i])
        return Kpoints(kpoints, method)

    def to_string(self) -> str:
        b1, b2, b3 = self.mesh
        s = f"""KPOINTS with kpts:{self.mesh}
0
{self.method}
{b1} {b2} {b3}
0 0 0"""
        return s

    def write(self, path="./", verbose=True):
        file_path = os.path.join(path, "KPOINTS")
        with open(file_path, "w") as f:
            f.write(self.to_string())
        if verbose:
            print("KPOINTS generated")


class Potcar(VaspInputBase):
    def __init__(self, atoms: Atoms, config: Optional[Dict[str, str]] = None):
        self.elements = sorted(set(atoms.symbols), key=lambda x: atoms.symbols.index(x))
        if not config:
            self.potcars = {elem: reccomended_potcars[elem] for elem in self.elements}
        else:
            self.potcars = config

    @staticmethod
    def from_config(atoms, config=None):
        elements = sorted(set(atoms.symbols), key=lambda x: atoms.symbols.index(x))

        potcar_params_collected = {}
        for key in elements:
            potcar_params_collected[key] = config[key]
        return Potcar(atoms, potcar_params_collected)

    def write(self, path="./", verbose=True):
        filename = os.path.join(path, "POTCAR")
        if os.path.exists(filename):
            os.remove(filename)
        for elem in self.elements:
            potacr_file = potcar_path / self.potcars[elem] / "POTCAR"
            with open(filename, "a") as f:
                f.write(potacr_file.read_text())
        if verbose:
            print("POTCAR generated")
