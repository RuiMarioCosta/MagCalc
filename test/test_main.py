import pytest
from magcalc.main import read_configuration


def test_readConfiguration_correctName_returnData():
    test_data = (
        {
            'phase': 'monoclinic',
            'J': 3.5,
            'gJ': 2,
            'Tc': 251,
            'theta_D': 250,
            'F0': 0.36,
            'Nm': 20,
            'N': 36,
        },
        {
            'phase': 'orthorhombic(I)',
            'J': 3.5,
            'gJ': 2,
            'Tc': 308,
            'theta_D': 278,
            'F0': 0,
            'Nm': 20,
            'N': 36,
        },
    )

    data = read_configuration(
        'sucessfulConfiguration', configuration_file='test/configuration.yaml'
    )

    assert data == test_data


def test_readConfiguration_wrongName_raiseError():
    with pytest.raises(ValueError) as exc_info:
        read_configuration(
            'someRandomNameNotInConfiguration',
            configuration_file='test/configuration.yaml',
        )
    assert isinstance(exc_info.value, ValueError)


def test_readConfiguration_insufficientParams_raiseError():
    with pytest.raises(KeyError) as exc_info:
        read_configuration(
            'insufficientParams', configuration_file='test/configuration.yaml'
        )
    assert isinstance(exc_info.value, KeyError)
