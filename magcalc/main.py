import yaml
import argparse
import csv
import numpy as np
from magnetocaloric import Magnetocaloric
from util import print_list_of_items


def read_configuration(name, configuration_file):
    if not name:
        print_list_of_items(configuration_file)
        name = input("Choose material from the list: ")

    with open(configuration_file) as file:
        data = yaml.safe_load(file)

    try:
        data = data[name]
    except KeyError:
        raise ValueError(f'Configuration name "{name}" does not exist in list.')

    phase1_data = data[0]
    phase2_data = data[1]
    params = ["phase", "J", "gJ", "Tc", "theta_D", "F0", "Nm", "N"]

    if not all(k in phase1_data for k in params) or not all(
        k in phase2_data for k in params
    ):
        raise KeyError("Missing necessary parameter.")

    return phase1_data, phase2_data


def main(args=None):
    parser = argparse.ArgumentParser(
        description="Compute magnetocaloric properties of the given material."
    )
    parser.add_argument(
        "--name",
        help="Name of configuration to use.",
    )
    parser.add_argument(
        "--configuration_file",
        default="magcalc/configuration.yaml",
        help="Path to a yaml file with configurations of magnetocaloric materials.",
    )
    parser.add_argument(
        "--initial_temperature",
        type=float,
        help="Initial temperature to compute properties.",
    )
    parser.add_argument(
        "--final_temperature",
        type=float,
        help="Final temperature to compute properties.",
    )
    parser.add_argument(
        "--delta_temperature",
        type=float,
        help="Interval between temperature points.",
    )
    parser.add_argument(
        "--initial_magnetic_field",
        type=float,
        help="Initial magnetic field to compute properties.",
    )
    parser.add_argument(
        "--final_magnetic_field",
        type=float,
        help="Final magnetic field to compute properties.",
    )
    parser.add_argument(
        "--delta_magnetic_field",
        type=float,
        help="Interval between magnetic field points.",
    )
    parser.add_argument(
        "--output",
        default="output.csv",
        help="Path to output csv file.",
    )
    args = parser.parse_args(args)

    data1, data2 = read_configuration(args.name, args.configuration_file)

    magcalc1 = Magnetocaloric(**data1)
    magcalc2 = Magnetocaloric(**data2)

    initial_temperature = args.initial_temperature
    if initial_temperature is None:
        initial_temperature = input("Set initial temperature (in K): ")

    final_temperature = args.final_temperature
    if final_temperature is None:
        final_temperature = input("Set final temperature (in K): ")

    delta_temperature = args.delta_temperature
    if delta_temperature is None:
        delta_temperature = input("Set delta temperature (in K): ")

    initial_magnetic_field = args.initial_magnetic_field
    if initial_magnetic_field is None:
        initial_magnetic_field = input("Set initial magnetic field (in T): ")

    final_magnetic_field = args.final_magnetic_field
    if final_magnetic_field is None:
        final_magnetic_field = input("Set final temperature (in T): ")

    delta_magnetic_field = args.delta_magnetic_field
    if delta_magnetic_field is None:
        delta_magnetic_field = input("Set delta temperature (in T): ")

    header = [
        f"Temperature",
        f"Free energy {magcalc1.phase}",
        f"Entropy {magcalc1.phase}",
        f"Energy {magcalc1.phase}",
        f"Free energy {magcalc2.phase}",
        f"Entropy {magcalc2.phase}",
        f"Energy {magcalc2.phase}",
    ]
    temperatures = np.arange(
        initial_temperature, final_temperature + delta_temperature, delta_temperature
    )
    magnetic_fields = np.arange(
        initial_magnetic_field,
        final_magnetic_field + delta_magnetic_field,
        delta_magnetic_field,
    )
    for b in magnetic_fields:
        free_energy_magcalc1 = []
        free_energy_magcalc2 = []
        entropy_magcalc1 = []
        entropy_magcalc2 = []
        energy_magcalc1 = []
        energy_magcalc2 = []
        for t in temperatures:
            free_energy_magcalc1.append(magcalc1.free_energy(t, b))
            free_energy_magcalc2.append(magcalc2.free_energy(t, b))
            entropy_magcalc1.append(magcalc1.entropy(t, b))
            entropy_magcalc2.append(magcalc2.entropy(t, b))
            energy_magcalc1.append(magcalc1.energy(t, b))
            energy_magcalc2.append(magcalc2.energy(t, b))

        list = [
            temperatures,
            free_energy_magcalc1,
            free_energy_magcalc2,
            entropy_magcalc1,
            entropy_magcalc2,
            energy_magcalc1,
            energy_magcalc2,
        ]

        with open(f"{b}T_{args.output}", mode="w", newline="") as output_file:
            output_writer = csv.writer(output_file, delimiter=",")
            output_writer.writerow(header)
            output_writer.writerows(zip(*list))


if __name__ == "__main__":
    main()
