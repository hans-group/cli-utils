#!/usr/bin/env python
from pathlib import Path
import argparse
import subprocess as sp

parser = argparse.ArgumentParser()
parser.add_argument("-p", type=str, default="g4")
parser.add_argument("-n", type=str, default="16")
parser.add_argument("-w", type=str, default="")
parser.add_argument("-J", type=str, default="orca_job")
args = parser.parse_args()


def check_input_file():
    all_inp_files = list(Path(".").glob("orca.inp"))
    if len(all_inp_files) == 0:
        raise FileNotFoundError(
            "No input file found. "
            "You must have a file named orca.inp in the current directory."
        )
    input_contents = all_inp_files[0].read_text()
    all_xyz_files = list(Path(".").glob("*.xyz"))
    if "xyzfile" in input_contents.lower() and len(all_xyz_files) == 0:
        raise FileNotFoundError(
            "No xyz file found even though xyzfile is specified in orca.inp. "
            "You must have a xyz file in the current directory."
        )


def main():
    source_dir = Path(__file__).parent
    orca_job_script = source_dir / "resources/submit_orca.sh"

    sbatch_args = ["sbatch"]
    sbatch_args.extend(["-N", "1"])
    sbatch_args.extend(["-p", args.p])
    sbatch_args.extend(["-n", args.n])
    if args.w:
        sbatch_args.extend(["-w", args.w])
    sbatch_args.extend(["-J", args.J])
    sbatch_args.append("submit_orca.sh")

    s = orca_job_script.read_text()
    Path("submit_orca.sh").write_text(s)
    sp.run(sbatch_args)


if __name__ == "__main__":
    main()
