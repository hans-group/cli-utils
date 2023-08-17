from dataclasses import dataclass
from enum import Enum
from typing import Any, Dict, Optional, Union

import matplotlib as mpl
import matplotlib.pyplot as plt
from ase import Atoms
from ase.formula import Formula
from ase import io
from ase.visualize.plot import plot_atoms
from pymatgen.core.structure import Structure
from pymatgen.io.ase import AseAtomsAdaptor


class JobStatus(Enum):
    """
    Enum class for job status
    """

    INITIALIZED = "INITIALIZED"
    # for situation when job is submitted but not tracked by sacct
    SUBMITTED = "SUBMITTED"
    PENDING = "PENDING"
    RUNNING = "RUNNING"
    COMPLETED = "COMPLETED"
    FAILED = "FAILED"


@dataclass
class Job:
    job_name: str
    job_type: str
    slurm_id: Optional[int]
    work_dir: str
    status: JobStatus

    # For JSON serialization
    def to_dict(self):
        return {
            "job_name": self.job_name,
            "job_type": self.job_type,
            "slurm_id": self.slurm_id,
            "work_dir": self.work_dir,
            "status": self.status.value,
        }

    # For deserialization from JSON
    @staticmethod
    def from_dict(d: Dict[str, Any]) -> "Job":
        return Job(
            job_name=d["job_name"],
            job_type=d["job_type"],
            slurm_id=d["slurm_id"],
            work_dir=d["work_dir"],
            status=JobStatus[d["status"]],
        )


def view_atoms(
    structure: Union[str, Atoms, Structure],
    rotation: str = "0x,0y,0z",
    save_to: Optional[str] = None,
) -> None:
    if isinstance(structure, Structure):
        atoms = AseAtomsAdaptor.get_atoms(structure)
    elif isinstance(structure, Atoms):
        pass
    elif isinstance(structure, str):
        atoms = io.read(atoms)
    else:
        raise RuntimeError("Format not recognized")
    fig, ax = plt.subplots(dpi=150)
    plot_atoms(atoms, ax, radii=0.5, rotation=rotation)
    ax.set_xticks([])
    ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)

    if save_to is not None:
        if save_to.endswith(".png"):
            fig.savefig(save_to, dpi=600)
        elif save_to.endswith(".pdf") or save_to.endswith(".svg"):
            fig.savefig(save_to)
    if mpl.get_backend() != "module://matplotlib_inline.backend_inline":
        fig.show()


def sort_atoms(filename):
    image = io.read(filename)
    copied = image.copy()

    atoms = Atoms()
    atoms.cell = copied.cell
    atoms.constraints = copied.constraints

    w = copied.symbols
    formula = Formula(f"{w}")
    count = list(formula.count())
    for _, symbol in enumerate(count):
        for i in range(len(copied)):
            if copied[i].symbol == symbol:
                atoms.append(copied[i])
    return atoms
