#!/usr/bin/env python
import ase.io
import click
import inquirer
import numpy as np
from ase import Atoms

from lib.vibration_analysis import VibAnalysis


# ========== IO utils ==========
def parse_matrix(s: str):
    lines = (line.strip() for line in s.splitlines())
    try:
        n_rows = int(next(lines))
    except StopIteration as e:
        raise ValueError("Empty file") from e
    except ValueError as e:
        raise ValueError("First line is not an integer") from e

    n_processed_cols = 0
    matrix = []
    while n_processed_cols < n_rows:
        n_cols = len(next(lines).split())
        matrix_fragment = []
        for _ in range(n_rows):
            matrix_elems = next(lines).split()[1:]
            matrix_elems = [float(i) for i in matrix_elems]
            matrix_fragment.append(matrix_elems)
        matrix_fragment = np.array(matrix_fragment)
        matrix.append(matrix_fragment)
        n_processed_cols += n_cols
    matrix = np.hstack(matrix)
    return matrix


def read_hessian(hess_file: str, symmetrize: bool = True):
    with open(hess_file, "r") as f:
        s = f.readlines()
        s = [i.strip() for i in s]
    start_line = s.index("$hessian")
    s = "\n".join(s[start_line + 1 :])
    hess = parse_matrix(s)
    if symmetrize:
        # check max error and raise error if larger than 1e-3
        max_error = np.max(np.abs(hess - hess.T))
        if max_error > 1e-2:
            raise ValueError(f"Symmetrization error: max error = {max_error}")
        hess = (hess + hess.T) / 2
    return hess


def animate_mode(coords, disp, scale, n_mesh=10):
    if coords.shape != disp.shape:
        raise ValueError("Coordinates and displacements must have the same shape")
    scale = abs(scale)
    frames = [coords]

    for i in range(1, n_mesh):
        frames.append(coords + scale * disp * (i + 1) / n_mesh)

    for i in range(n_mesh, 0, -1):
        frames.append(coords + scale * disp * (i - 1) / n_mesh)

    for i in range(1, n_mesh):
        frames.append(coords - scale * disp * (i + 1) / n_mesh)

    for i in range(n_mesh, 0, -1):
        frames.append(coords - scale * disp * (i - 1) / n_mesh)

    return frames


@click.command(help="Animate a vibrational mode from ORCA calculation.")
@click.option("-s", "--scale", default=1.0, help="Scale factor for displacements.")
@click.option("-n", "--num_points", default=10, help="Number of points in the mesh.")
@click.option("--xyz_file", default="orca.xyz", help="XYZ file to read.")
@click.option("--hess_file", default="orca.hess", help="Hessian file to read.")
def main(scale, num_points, xyz_file, hess_file):
    atoms = ase.io.read(xyz_file)
    hess = read_hessian(hess_file)
    ana = VibAnalysis(atoms, hess)
    ana._compute_vib_props()
    ana.displacements.shape, ana.vib_freqs_cm.shape, ana.force_constant.shape

    species = atoms.get_chemical_symbols()
    coords = atoms.get_positions()
    freqs = ana.vib_freqs_cm

    answers = inquirer.prompt(
        [
            inquirer.List(
                "freq",
                message="Which frequency do you want to animate?",
                choices=freqs.tolist(),
            )
        ]
    )
    selected_freq = answers["freq"]
    selected_idx = np.argmin(np.abs(freqs - selected_freq))

    # Get the normal modes of the last geometry.
    disp = ana.displacements[selected_idx]

    # Create a list of frames.
    frames = animate_mode(coords, disp, scale=scale, n_mesh=num_points)
    atomslist = [Atoms(species, frame) for frame in frames]

    ase.io.write(f"vibmode_{selected_freq:.2f}cm-1.xyz", atomslist, format="xyz")


if __name__ == "__main__":
    main()
