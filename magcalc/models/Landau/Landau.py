import numpy as np
import matplotlib.pyplot as plt
import yaml
from util import print_list_of_items


class Landau:
    def __init__(self):
        self.a0 = 0
        self.a1 = 0
        self.a2 = 0
        self.a3 = 0
        self.a4 = 0
        self.a5 = 0
        self.b0 = 0
        self.b1 = 0
        self.b2 = 0
        self.b3 = 0
        self.b4 = 0
        self.b5 = 0
        self.c0 = 0
        self.c1 = 0
        self.c2 = 0
        self.c3 = 0
        self.c4 = 0
        self.c5 = 0

    def set_A_coefs(self, a0, a1, a2, a3, a4, a5):
        self.a0 = a0
        self.a1 = a1
        self.a2 = a2
        self.a3 = a3
        self.a4 = a4
        self.a5 = a5

    def set_B_coefs(self, b0, b1, b2, b3, b4, b5):
        self.b0 = b0
        self.b1 = b1
        self.b2 = b2
        self.b3 = b3
        self.b4 = b4
        self.b5 = b5

    def set_C_coefs(self, c0, c1, c2, c3, c4, c5):
        self.c0 = c0
        self.c1 = c1
        self.c2 = c2
        self.c3 = c3
        self.c4 = c4
        self.c5 = c5

    def coef_A(self, T):
        return (
            self.a0
            + self.a1 * T
            + self.a2 * (T ** 2.0)
            + self.a3 * (T ** 3.0)
            + self.a4 * (T ** 4.0)
            + self.a5 * (T ** 5.0)
        )

    def coef_B(self, T):
        return (
            self.b0
            + self.b1 * T
            + self.b2 * (T ** 2.0)
            + self.b3 * (T ** 3.0)
            + self.b4 * (T ** 4.0)
            + self.b5 * (T ** 5.0)
        )

    def coef_C(self, T):
        return (
            self.c0
            + self.c1 * T
            + self.c2 * (T ** 2.0)
            + self.c3 * (T ** 3.0)
            + self.c4 * (T ** 4.0)
            + self.c5 * (T ** 5.0)
        )

    def G(self, M, T, H):
        return (
            1.0 / 2.0 * self.coef_A(T) * (M ** 2.0)
            + 1.0 / 4.0 * self.coef_B(T) * (M ** 4.0)
            + 1.0 / 6.0 * self.coef_C(T) * (M ** 6.0)
            - H * M
        )

    def dGdM(self, M, T, H):
        return (
            self.coef_A(T) * M
            + self.coef_B(T) * (M ** 3.0)
            + self.coef_C(T) * (M ** 5.0)
            - H
        )

    def plot_G(self, M, T, H):
        """Plots the Landau free energy as a function of temperature or magnetic field.

        Parameters
        ----------
        M : array
            Array with the magnetizations.
        T : scalar
            Temperature.
        H : scalar
            Magnetic field.
        """
        plt.figure()

        M = np.asarray(M)
        T = np.asarray(T)
        H = np.asarray(H)
        landau_free_energy = []
        if M.ndim == 1 and T.ndim == 0 and H.ndim == 0:
            for m in M:
                landau_free_energy.append(self.G(m, T, H))
            plt.plot(M, landau_free_energy, label=f'T={T} K, H={H} T')
            plt.xlabel("M (K)")
        else:
            raise ValueError("Temperature and Magnetic Field need to be a combination of array and scalar")

        plt.title("Landau Free Energy, $G$")
        plt.legend(loc=0, fontsize="small")
        plt.ylabel("$G(M, T, H)$")
        plt.show()


def read_landau_db(name, configuration_file):
    if not name:
        print_list_of_items(configuration_file)
        name = input('Choose material from the list: ')

    with open(configuration_file) as file:
        data = yaml.safe_load(file)

    try:
        data = data[name]
    except KeyError:
        raise ValueError(f'Configuration name "{name}" does not exist in database.')

    return data


if __name__ == "__main__":
    data = read_landau_db('LaFe11.6Si1.4', 'magcalc/models/Landau/Landau_db.yaml')
    a = data['a']
    b = data['b']
    c = data['c']

    landau = Landau()
    landau.set_A_coefs(*a)
    landau.set_B_coefs(*b)
    landau.set_C_coefs(*c)

    print(landau.G(0, 300, 0))
    print(landau.G(1, 300, 0))
    print(landau.G(-1, 300, 0))

    magnetization = np.linspace(-400, 400, 1000)
    landau.plot_G(magnetization, 300, 0)
    landau.plot_G(magnetization, 1, 0)
