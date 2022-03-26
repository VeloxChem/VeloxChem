import numpy as np
import pytest

from veloxchem.numerovdriver import NumerovDriver
from veloxchem.veloxchemlib import bohr_in_angstroms

class TestNumerov:

    def read_data(self, num_drv):

        # data for HCl (B3LYP/6-31+G)

        # bond_lengths
        bond_lengths = np.arange(-0.7, 2.01, 0.1) * bohr_in_angstroms() + 1.274

        # ground state energies from unrestricted SCF
        gs_energies = [-460.56992716015156, -460.6384446308048,
                       -460.68782504635186, -460.7231200805613,
                       -460.74738848682415, -460.76316652211233,
                       -460.7727995818075,  -460.7776159219715,
                       -460.7787520592223,  -460.7772924833631,
                       -460.773876362175,   -460.76895973907705,
                       -460.7629956534169,  -460.7563464669844,
                       -460.749195786104,   -460.74171325709665,
                       -460.7340626250854,  -460.7263637578346,
                       -460.71868028553394, -460.71107803329744,
                       -460.7036102246686,  -460.69631334996706,
                       -460.68921356046394, -460.68233412001865,
                       -460.6756969071502,  -460.6693081100631,
                       -460.6631656728619,  -460.657274038127]

        # excited state (2) energies from linear eigensolver
        es_energies = [-460.19069874760595, -460.2608580619231,
                       -460.3117947851528,  -460.3486079617552,
                       -460.3744458876339,  -460.3920655532207,
                       -460.40432853489216, -460.4132561295106,
                       -460.4194047263074,  -460.4225198939118,
                       -460.4228194704999,  -460.42085950086914,
                       -460.41728579953207, -460.41262346749215,
                       -460.4071732207166,  -460.3975116567743,
                       -460.3948896316938,  -460.3884317986093,
                       -460.3819047984704,  -460.37539335329643,
                       -460.3689623048443,  -460.3626551120557,
                       -460.3565041832953,  -460.35053887440154,
                       -460.34478533713377, -460.33925192007575,
                       -460.33393789285736, -460.32885039158884]

        # transition dipole moments [x, y, z]
        properties = [[-1.28759568e-15,  2.22459941e-01, -1.43256940e-01],
                      [ 1.13158158e-14,  2.41481132e-01, -1.20318139e-01],
                      [-1.35713022e-14, -2.47360475e-02, -2.78665367e-01],
                      [-3.76303645e-15,  1.59584399e-01, -2.53214738e-01],
                      [-5.75001380e-14,  3.28555800e-01, -7.88746252e-02],
                      [-7.79643013e-14, -2.39090492e-03, -4.14935710e-01],
                      [ 1.94516569e-13,  5.03487289e-01,  2.41961461e-01],
                      [ 1.27750474e-13, -6.53114707e-01,  3.40754360e-01],
                      [ 2.69960488e-14, -4.47108144e-01, -7.16400586e-01],
                      [-1.79971554e-16,  8.73805537e-01,  1.67507207e-01],
                      [-1.34443986e-12,  1.47391285e-01, -8.96353460e-01],
                      [-4.18936256e-15, -5.80167260e-02, -9.13754249e-01],
                      [ 8.39979975e-15, -9.07639464e-01,  1.31932596e-01],
                      [-1.22658857e-14, -1.13042000e-01,  9.08799260e-01],
                      [-2.74680958e-14,  8.04814063e-01,  4.30648152e-01],
                      [ 1.72351893e+00,  2.39232628e-13, -3.34945278e-13],
                      [ 1.00658893e-14,  9.02612082e-01, -5.88757574e-02],
                      [ 2.11558070e-15,  8.30906521e-01, -3.45939988e-01],
                      [-1.00392560e-14,  8.94391145e-01,  4.66713511e-02],
                      [-4.51255308e-15,  2.83894124e-02,  8.90908368e-01],
                      [ 5.50731026e-15, -8.77165252e-01, -1.34326959e-01],
                      [-2.26963916e-15, -2.11237302e-03,  8.83748280e-01],
                      [ 2.57770382e-15, -9.97110967e-02,  8.74809493e-01],
                      [ 3.12763416e-15,  5.55705394e-01, -6.79219066e-01],
                      [-9.80026117e-16, -8.35646830e-01, -2.59705945e-01],
                      [ 2.80861405e-15,  2.52749744e-01, -8.35547938e-01],
                      [ 2.55382733e-15,  1.80494108e-02, -8.70969843e-01],
                      [-5.42165181e-15, -7.40960191e-01, -4.55375942e-01]]

        num_drv.read_pec_data(bond_lengths, properties, gs_energies,
                              es_energies)
        num_drv.set_reduced_mass(0.9795925092006827)

    def run_numerov(self, numerov_dict, data_lines):

        ref_exc = {}
        ref_osc = {}

        for b in data_lines.keys():
            ref_exc[b] = [float(line.split()[0]) for line in data_lines[b]]
            ref_osc[b] = [float(line.split()[2]) for line in data_lines[b]]

        num_drv = NumerovDriver()

        num_drv.update_settings(numerov_dict)

        self.read_data(num_drv)

        results = num_drv.compute()

        assert num_drv.is_converged

        exc = results['excitation_energies']
        osc = results['oscillator_strengths']

        print(exc, ref_exc)

        if len(exc.keys()) == 3:
            for b in 'PQR':
                assert np.max(np.abs(exc[b] - ref_exc[b])) < 1.0e-2
                assert np.max(np.abs(osc[b] - ref_osc[b])) < 1.0e-4

        if len(exc.keys()) == 2:
            for b in ['absorption', 'emission']:
                assert np.max(np.abs(exc[b] - ref_exc[b])) < 1.0e-2
                #assert np.max(np.abs(osc[b] - ref_osc[b])) < 1.0e-4

    def test_IR_spectrum(self):

        #    ----------------------
        #    Rovibronic IR Spectrum                                                                                     
        #    ----------------------                                                                                     
        #            Excitation energy    Oscillator strength
        # P Branch
        raw_p_data = """
                          2711.26 cm-1            4.00736e-09
                          2691.71 cm-1            5.53061e-09
                          2672.16 cm-1            5.83444e-09
                          2652.62 cm-1            5.14370e-09
                          2633.07 cm-1            3.92274e-09
                          2613.52 cm-1            2.63231e-09
                          2593.97 cm-1            1.56933e-09
                          2574.43 cm-1            8.36248e-10
                          2554.88 cm-1            3.99885e-10
                          2535.33 cm-1            1.72080e-10
                          2515.78 cm-1            6.67727e-11
                          2496.24 cm-1            2.33994e-11
                          2476.69 cm-1            7.41405e-12
                          2457.14 cm-1            2.12593e-12
        """
        # Q Branch
        raw_q_data = """
                          2730.81 cm-1            6.02968e-08
        """
        # R Branch
        raw_r_data = """
                          2750.36 cm-1            1.39977e-09
                          2769.90 cm-1            3.82130e-09
                          2789.45 cm-1            5.27383e-09
                          2809.00 cm-1            5.56355e-09
                          2828.55 cm-1            4.90488e-09
                          2848.09 cm-1            3.74061e-09
                          2867.64 cm-1            2.51009e-09
                          2887.19 cm-1            1.49647e-09
                          2906.74 cm-1            7.97422e-10
                          2926.28 cm-1            3.81319e-10
                          2945.83 cm-1            1.64090e-10
                          2965.38 cm-1            6.36725e-11
                          2984.93 cm-1            2.23130e-11
                          3004.47 cm-1            7.06982e-12
                          3024.02 cm-1            2.02723e-12
        """

        data_lines = {}
        data_lines['P'] = raw_p_data.splitlines()[1:-1]
        data_lines['Q'] = raw_q_data.splitlines()[1:-1]
        data_lines['R'] = raw_r_data.splitlines()[1:-1]

        numerov_dict = {
            'el_transition': 'no',
            'n_vib_states': 2,
            'n_rot_states': len(data_lines['R']),
        }

        self.run_numerov(numerov_dict, data_lines)

    def test_UV_spectrum(self):

        #       ----------------------------------
        #       UV/vis Absorption/Emission Spectrum
        #       -----------------------------------
        #                 Excitation energy    Oscillator strength
        # Absorption
        raw_abs_data = """
                              78549.62 cm-1            9.96843e-03
                              80930.00 cm-1            3.35285e-03
                              83217.43 cm-1            2.66584e-03
                              85411.91 cm-1            5.86076e-04
                              87513.43 cm-1            5.64929e-05
                              89522.00 cm-1            1.10832e-05
                              91437.61 cm-1            9.30815e-06
                              93260.27 cm-1            5.18215e-06
                              94989.97 cm-1            2.08770e-06
                              96626.72 cm-1            6.89661e-07
        """
        # Emission
        raw_em_data = """
                              78549.62 cm-1            9.96843e-03
                              75818.81 cm-1            1.44056e-02
                              73191.54 cm-1            1.43757e-03
                              70667.81 cm-1            1.55611e-04
                              68247.61 cm-1            1.36733e-05
                              65930.96 cm-1            1.43414e-06
                              63717.84 cm-1            3.26145e-09
                              61608.26 cm-1            2.74084e-09
                              59602.22 cm-1            1.56605e-10
                              57699.71 cm-1            8.06237e-12
        """

        data_lines = {}
        data_lines['absorption'] = raw_abs_data.splitlines()[1:-1]
        data_lines['emission'] = raw_em_data.splitlines()[1:-1]

        numerov_dict = {
            'el_transition': 'yes',
            'final_el_state': 2,
            'n_vib_states': len(data_lines['absorption']),
        }

        self.run_numerov(numerov_dict, data_lines)


