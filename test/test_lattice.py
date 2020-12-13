from lattice import Lattice


class TestLattice:
    def test_latticeEntropy_temperature_getEntropy(self):
        lattice = Lattice(123, 1)
        lattice2 = Lattice(123, 2)
        expected = 0.00035104320939824554

        result = lattice.lattice_entropy(123)
        result2 = lattice2.lattice_entropy(123)

        assert result == expected
        assert result2 == 2 * expected

    def test_latticeFreeEnergy_temperature_getFreeEnergy(self):
        lattice = Lattice(123, 1)
        lattice2 = Lattice(123, 2)
        expected = -0.009809044240247611

        result = lattice.lattice_free_energy(123)
        result2 = lattice2.lattice_free_energy(123)

        assert result == expected
        assert result2 == 2 * expected

