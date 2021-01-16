import argparse
from util import print_list_of_items


def main(args=None):
    parser = argparse.ArgumentParser(
        description='Compute magnetocaloric properties of the given material.'
    )
    parser.add_argument(
        '--model',
        default='Landau',
        help='Name of model to use.',
    )
    parser.add_argument(
        '--name',
        help='Name of configuration to use.',
    )
    args = parser.parse_args(args)

    if args.model == 'Landau':
        from magcalc.models.Landau.Landau import Landau, read_landau_db
        database = 'magcalc/models/Landau/Landau_db.yaml'
        data = read_landau_db(args.name, database)
    else:
        raise ValueError(f'Model chosen, {args.model}, not recognized.')

    print(1)


if __name__ == '__main__':
    main()
