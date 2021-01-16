import yaml
import argparse
from util import print_list_of_items
from magcalc.models.Landau.Landau import Landau


def read_landau_db(arg, configuration_file):
    with open(configuration_file) as file:
        data = yaml.safe_load(file)

    try:
        data = data[arg]
    except KeyError:
        raise ValueError(f'Configuration name "{arg}" does not exist in list.')

    return data


def main(args=None):
    parser = argparse.ArgumentParser(
        description='Compute magnetocaloric properties of the given material.'
    )
    parser.add_argument(
        '--name',
        help='Name of configuration to use.',
    )
    parser.add_argument(
        '--database_file',
        default='magcalc/models/Landau/Landau_db.yaml',
        help='Path to a yaml file with Landau fits of A, B and C coefficients.',
    )
    args = parser.parse_args(args)
    name = args.name

    if not name:
        print_list_of_items(args.database_file)
        name = input('Choose material from the list: ')

    data = read_landau_db(name, args.database_file)

    print(1)


if __name__ == '__main__':
    main()
