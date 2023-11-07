#!/usr/bin/env python
import os
import sys

import click

AVAILABLE_MODES = ["vasp"]


@click.command(help="Clean up the outputs of the simulation results")
@click.argument("mode", default="vasp")
def main(mode: str):
    mode = mode.strip().lower()
    if mode not in AVAILABLE_MODES:
        click.echo(
            "Error: Unsupported mode. Supported modes are: {}".format(
                ", ".join(AVAILABLE_MODES)
            )
        )
        sys.exit(1)

    if mode == "vasp":
        file_list = [
            "CHG",
            "CHGCAR",
            "CONTCAR",
            "DOSCAR",
            "EIGENVAL",
            "IBZKPT",
            "OSZICAR",
            "OUTCAR",
            "PCDAT",
            "REPORT",
            "stdout.log",
            "vasprun.xml",
            "WAVECAR",
            "XDATCAR",
            "vaspout.h5",
        ]
        for f in file_list:
            if os.path.isfile(f):
                os.remove(f)
            else:
                pass
        os.system("rm slurm-*")
        os.system("rm stderr-*")


if __name__ == "__main__":
    main()
