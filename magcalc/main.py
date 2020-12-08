import yaml
import argparse
from magnetocaloric import Magnetocaloric


def print_list_of_materials():
    with open('magcalc/configuration.yaml') as file:
        data = yaml.load(file)
    print('\n'.join(data.keys()))


def read_configuration(arg, configuration_file='magcalc/configuration.yaml'):
    with open(configuration_file) as file:
        data = yaml.load(file)

    try:
        data = data[arg]
    except KeyError:
        raise ValueError(f'Configuration name "{arg}" does not exist in list.')

    phase1_data = data[0]
    phase2_data = data[1]
    params = ['phase', 'J', 'gJ', 'Tc', 'theta_D', 'F0', 'Nm', 'N']

    if not all(k in phase1_data for k in params) or not all(k in phase2_data for k in params):
        raise KeyError('Missing necessary parameter.')

    return phase1_data, phase2_data


def main(args=None):
    parser = argparse.ArgumentParser(
        description='Compute magnetocaloric properties of the given material.'
    )
    parser.add_argument(
        '--name',
        help='Name of configuration to use.',
    )
    args = parser.parse_args(args)
    name = args.name

    if not name:
        print_list_of_materials()
        name = input('Choose material from the list: ')

    data1, data2 = read_configuration(name)

    magcalc1 = Magnetocaloric(**data1)
    magcalc2 = Magnetocaloric(**data2)


if __name__ == '__main__':
    main()
