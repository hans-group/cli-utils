#!/usr/bin/env python
import click
from pathlib import Path
from typing import Dict, Any
from lib.vasp_set.io import InputGen
import toml


@click.command(help="Setting up input files for VASP calculation.")
@click.option("-d", "--directory", default="./", help="calculation directory")
@click.option(
    "-s",
    "--vasp_settings",
    default="input_params.toml",
    help="vasp parameters settings TOML file",
)
def main(directory: str, vasp_settings: str):
    if not directory:
        workdir = Path.cwd()
    else:
        workdir = Path(directory)

    with open(vasp_settings, "r") as f:
        config = toml.load(f)
        incar_params: Dict[str, Any] = config["incar"]
        kpoints_params: Dict = config["kpoints"]
        if "potcar" in config:
            potcar_params = config["potcar"]
        else:
            potcar_params = None

    inputs = InputGen.from_dict(
        workdir / Path("POSCAR"),
        incar_params=incar_params,
        kpoints_params=kpoints_params,
        potcar_params=potcar_params,
    )
    inputs.write(path=workdir, verbose=True)


if __name__ == "__main__":
    main()
