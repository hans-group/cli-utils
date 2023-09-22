#!/usr/bin/env python
import ase.io
import click
import yaml


N_A = 6.022140857e23


def get_volume(spec: dict):
    # Ignore the partial volume of components
    total_volume = 0
    for filename in spec:
        if spec[filename].get("ignore", False):
            continue
        atoms = ase.io.read(filename)
        molar_mass = atoms.get_masses().sum()
        mass = spec[filename]["number"] * molar_mass / N_A  # g
        density = spec[filename]["density"]  # g/cm^3
        volume = mass / density  # cm^3
        total_volume += volume
    return total_volume * 10**24  # Å^3


def generate_packmol_input(spec: dict, output: str, tolerance: float, volume: float):
    side_length = volume ** (1 / 3)
    side_length -= 1.5
    side_length = round(side_length, 3)
    s = "filetype xyz\n"
    s += "tolerance {}\n".format(tolerance)
    s += "output {}\n\n".format(output)
    s += f"# original side length: {volume ** (1 / 3):.3f}\n\n"
    for filename in spec:
        s += "structure {}\n".format(filename)
        s += "  number {}\n".format(spec[filename]["number"])
        s += "  inside box 0. 0. 0. {0} {0} {0}\n".format(side_length)
        s += "end structure\n\n"
    return s


@click.command()
@click.argument("spec", type=click.File("r"), default="spec.yaml")
@click.option("--generate-packmol", is_flag=True)
@click.option("--tolerance", default=2.0)
@click.option("-o", "--output", type=click.File("w"), default="pack.inp")
@click.option("--packmol-output", type=str, default="packed_box.xyz")
def main(spec, generate_packmol, tolerance, output, packmol_output):
    spec = yaml.safe_load(spec)
    volume = get_volume(spec)
    print(f"Volume: {volume} Å^3")
    print(f"Side length (cube): {volume**(1/3)} Å^3")
    if generate_packmol:
        packmol_input = generate_packmol_input(spec, packmol_output, tolerance, volume)
        output.write(packmol_input)


if __name__ == "__main__":
    main()

