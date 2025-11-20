import math

from veloxchem.mmdriver import MMDriver


class TestMMDihedralPotential:

    def test_mm_dihedral_potential(self):

        mm_drv = MMDriver()

        # prepare parameters

        barriers = []
        phases = []
        periodicities = []

        for barrier in [0.12, 3.45]:
            for phase in [0.0, 180.0]:
                for periodicity in range(1, 6 + 1):
                    barriers.append(barrier)
                    phases.append(phase)
                    periodicities.append(periodicity)

        # compare RB with Fourier

        for phi in range(0, 360, 5):
            phi_in_radian = phi * math.pi / 180.0
            cos_phi = math.cos(phi_in_radian)

            sum_Fourier_ene = 0.0

            # compare single term

            for barrier, phase, periodicity in zip(barriers, phases,
                                                   periodicities):
                phase_in_radian = phase * math.pi / 180.0

                RB_coefs = mm_drv.get_RB_coefficients_single_term(
                    barrier, phase, periodicity)

                c0, c1, c2, c3, c4, c5, c6 = RB_coefs

                RB_ene = (c0 - c1 * cos_phi + c2 * cos_phi**2 -
                          c3 * cos_phi**3 + c4 * cos_phi**4 - c5 * cos_phi**5 +
                          c6 * cos_phi**6)

                Fourier_ene = barrier * (
                    1.0 +
                    math.cos(periodicity * phi_in_radian - phase_in_radian))

                sum_Fourier_ene += Fourier_ene

                assert abs(RB_ene - Fourier_ene) < 1.0e-10

            # compare multiple terms

            sum_RB_coefs = mm_drv.get_RB_coefficients_multiple_terms(
                barriers, phases, periodicities)

            c0, c1, c2, c3, c4, c5, c6 = sum_RB_coefs

            sum_RB_ene = (c0 - c1 * cos_phi + c2 * cos_phi**2 -
                          c3 * cos_phi**3 + c4 * cos_phi**4 - c5 * cos_phi**5 +
                          c6 * cos_phi**6)

            assert abs(sum_RB_ene - sum_Fourier_ene) < 1.0e-10
