import numpy as np
import math

from veloxchem.mmdriver import MMDriver


class TestMMDihedralPotential:

    def test_mm_dihedral_potential(self):

        params = []
        for barrier in [0.12, 3.45]:
            for phase in [0.0, 180.0]:
                for periodicity in range(1, 6 + 1):
                    params.append((barrier, phase, periodicity))

        for barrier, phase, periodicity in params:
            phase_in_radian = phase * math.pi / 180.0

            RB_coefs = MMDriver.get_RB_coefficients(barrier, phase, periodicity)

            c0, c1, c2, c3, c4, c5, c6 = RB_coefs

            RB_energies = []
            Fourier_energies = []

            for phi in range(0, 360, 5):
                phi_in_radian = phi * math.pi / 180.0
                cos_phi = math.cos(phi_in_radian)

                RB_ene = (c0 - c1 * cos_phi + c2 * cos_phi**2 -
                          c3 * cos_phi**3 + c4 * cos_phi**4 - c5 * cos_phi**5 +
                          c6 * cos_phi**6)

                Fourier_ene = barrier * (
                    1.0 +
                    math.cos(periodicity * phi_in_radian - phase_in_radian))

                RB_energies.append(RB_ene)
                Fourier_energies.append(Fourier_ene)

            RB_energies = np.array(RB_energies)
            Fourier_energies = np.array(Fourier_energies)

            assert np.max(np.abs(RB_energies - Fourier_energies)) < 1.0e-10
