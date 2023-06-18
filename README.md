# Command-line utilities

Useful command-line utilities for research in HANS group.


## Installation

First, make sure that `cargo` is installed.

```
git clone https://github.com/hans-group/cli-utils
cd cli-utils
bash install.sh
```

If it fails, reinstall the latest `rustup` and try again.

## List of commands

- `gotojob`: Go to SLURM job submission directory
- `animate_vibmode_orca`: Animate the normal mode of molecule to `.xyz` file from orca frequency calculation.
- `submit_orca`: Submit ORCA job to SLURM in one line
- `submit_vasp`: Submit VASP job to SLURM in one line
- `vasp_mdstat`: Plot potential/kinetic/total energy and temperature profile of VASP MD run.
- `check_ef_vasp`: Check convergence of energy and force from VASP geometry optimization run.
- `kpointsgen`: Generate `KPOINTS` file for VASP calculation, either from explicit grid or from density of k-points.
- `pos2pot`: Generate `POTCAR` file from `POSCAR` file.
- `ndstat`: Prettier version of `pestat`
- `qst`: Prettier version of `qst`
- `symmetrize_poscar`: Symmetrize POSCAR

See wiki for detailed usage.

## For developers & contributors

### Bash scripts

- Write a command as function definition. 
- If it is short, append to `src/bash/cli_utils_bash.sh`. Otherwise, write to separate file and add `source $source_dir/filename` to `cli_utils_bash.sh`.

### Python scripts

- In the directory `src/python`, **only executable scripts are allowed.** Do not put the library code in it.
- The library codes should be placed in `src/python/lib`.
- Check the errors with `flake8` and format the code with `black`.
