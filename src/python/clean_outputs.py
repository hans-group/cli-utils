#!/usr/bin/env python
import os
import sys
import glob
import click

AVAILABLE_MODES = ["vasp"]


def remove_matched_files(pattern: str, dry_run: bool = False):
    """Removes the files that match the pattern
    # ex) */file.txt, file*.txt, file.txt, ...
    """

    for f in glob.glob(pattern):
        if dry_run:
            click.echo("Dry run: will remove {}".format(f))
        else:
            os.remove(f)


@click.command(help="Clean up the outputs of the simulation results")
@click.argument("mode", default="vasp")
@click.option("--keep-log", is_flag=True, help="Keep the log files.")
@click.option("--keep-xml", is_flag=True, help="Keep the xml files (only for vasp).")
@click.option("--keep-chgcar", is_flag=True, help="Keep the CHG and CHGCAR files (only for vasp).")
@click.option("--keep-wavecar", is_flag=True, help="Keep the WAVECAR files (only for vasp).")
@click.option(
    "--dry-run", is_flag=True, help="Dry run. (only shows the files to be removed and does not remove them actually)"
)
def main(mode: str, keep_log: bool, keep_xml: bool, keep_chgcar: bool, keep_wavecar: bool, dry_run: bool):
    mode = mode.strip().lower()
    if mode not in AVAILABLE_MODES:
        click.echo("Error: Unsupported mode. Supported modes are: {}".format(", ".join(AVAILABLE_MODES)))
        sys.exit(1)
    
    file_list = []
    if mode == "vasp":
        file_list += [
            "CONTCAR",
            "DOSCAR",
            "EIGENVAL",
            "IBZKPT",
            "OSZICAR",
            "OUTCAR",
            "PCDAT",
            "REPORT",
            "XDATCAR",
            "vaspout.h5",
        ]
        # Vasp specific options
        if not keep_xml:
            file_list.append("./*.xml")
        if not keep_chgcar:
            file_list.append("CHG")
            file_list.append("CHGCAR")
        if not keep_wavecar:
            file_list.append("WAVECAR")

    # Global options
    if not keep_log:
        file_list.append("./*.log")
        file_list.append("./*.out")

    # Remove the files
    for f in file_list:
        remove_matched_files(f, dry_run=dry_run)


if __name__ == "__main__":
    main()
