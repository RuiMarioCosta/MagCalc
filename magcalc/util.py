import yaml


def print_list_of_items(yaml_file):
    with open(yaml_file) as file:
        data = yaml.safe_load(file)
    print('\n'.join(data.keys()))
