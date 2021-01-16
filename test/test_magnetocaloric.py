from magnetocaloric import Magnetocaloric
from magnetization import Magnetization
from lattice import Lattice
from electronic import Electronic
import pytest


class TestMagnetocaloric:
    def test_entropy_equalsSumOfParts(self):
        J = 1
        gJ = 1
        Tc = 123
        Nm = 1
        theta_D = 123
        N = 1
        F0 = 1
        mc = Magnetocaloric(J, gJ, Tc, Nm, theta_D, N, F0)
        mag = Magnetization(J, gJ, Tc, Nm)
        lat = Lattice(theta_D, N)
        temperature = 300
        magnetic_field = 0
        expected = mag.magnetic_entropy(
            temperature, magnetic_field
        ) + lat.lattice_entropy(temperature)

        result = mc.entropy(temperature, magnetic_field)

        assert result == expected

    def test_energy_equalsSumOfParts(self):
        J = 1
        gJ = 1
        Tc = 123
        Nm = 1
        theta_D = 123
        N = 1
        F0 = 1
        mc = Magnetocaloric(J, gJ, Tc, Nm, theta_D, N, F0)
        mag = Magnetization(J, gJ, Tc, Nm)
        lat = Lattice(theta_D, N)
        elec = Electronic(F0)
        temperature = 300
        magnetic_field = 0
        expected = (
            mag.magnetic_energy(temperature, magnetic_field)
            + lat.lattice_energy(temperature)
            + elec.electronic_energy()
        )

        result = mc.energy(temperature, magnetic_field)

        assert result == expected

    def test_freeEnergy_equalsSumOfParts(self):
        J = 1
        gJ = 1
        Tc = 123
        Nm = 1
        theta_D = 123
        N = 1
        F0 = 1
        mc = Magnetocaloric(J, gJ, Tc, Nm, theta_D, N, F0)
        mag = Magnetization(J, gJ, Tc, Nm)
        lat = Lattice(theta_D, N)
        elec = Electronic(F0)
        temperature = 300
        magnetic_field = 0
        expected = (
            mag.magnetic_free_energy(temperature, magnetic_field)
            + lat.lattice_free_energy(temperature)
            + elec.electronic_free_energy()
        )

        result = mc.free_energy(temperature, magnetic_field)

        assert result == expected

    @pytest.mark.parametrize(
        "temperature, magnetic_field, expected",
        [
            (1, 0, 0.9960252545283329),
            (1, 1, 0.9959673707102669),
            (1, -1, 0.9959673707102669),
            (100000, 0, -190.31197909563747),
            (100000, 100000, -191.56501587560945),
            (300, 0, 0.8769235083464221),
            (300, 100000, -4.898955954721044),
        ],
    )
    def test_freeEnergy_atSomeTemperatureAndMagneticField_getFreeEnergy(
        self, temperature, magnetic_field, expected
    ):
        J = 1
        gJ = 1
        Tc = 123
        Nm = 1
        theta_D = 123
        N = 1
        F0 = 1
        mc = Magnetocaloric(J, gJ, Tc, Nm, theta_D, N, F0)

        result = mc.free_energy(temperature, magnetic_field)

        assert abs(result - expected) < 1e-10
