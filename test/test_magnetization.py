from magnetization import Magnetization
import pytest


class TestMagnetization:
    @pytest.mark.parametrize(
        "temperature, magnetic_field, expected",
        [
            (1, 0, 1),
            (1, 1, 1),
            (1, -1, -1),
            (100000, 0, -7.0150551504213795e-09),
            (100000, 100000, 0.41753673678058895),
            (300, 0, 5.097206453848289e-10),
            (300, 100000, 1),
        ],
    )
    def test_reducedMagnetization_atSomeTemperatureAndMagneticField_getMagnetization(
        self, temperature, magnetic_field, expected
    ):
        J = 1
        gJ = 1
        Tc = 123
        Nm = 1
        mag = Magnetization(J, gJ, Tc, Nm)
        mag2 = Magnetization(J, gJ, Tc, 2 * Nm)

        result = mag.reduced_magnetization(temperature, magnetic_field)
        result2 = mag2.reduced_magnetization(temperature, magnetic_field)

        assert result - expected < 1e-10
        assert result2 - expected < 1e-10

    @pytest.mark.parametrize(
        "temperature, magnetic_field, expected",
        [
            (1, 0, 0),
            (1, 1, 0),
            (1, -1, 0),
            (100000, 0, 9.467107270177856e-05),
            (100000, 100000, 8.30051022037348e-05),
            (300, 0, 9.467107270177856e-05),
            (300, 100000, 0),
        ],
    )
    def test_magneticEntropy_atSomeTemperatureAndMagneticField_getEntropy(
        self, temperature, magnetic_field, expected
    ):
        J = 1
        gJ = 1
        Tc = 123
        Nm = 1
        mag = Magnetization(J, gJ, Tc, Nm)

        result = mag.magnetic_entropy(temperature, magnetic_field)

        assert result - expected < 1e-10
