from pathlib import Path
from unittest.mock import patch
import numpy as np
import tempfile
import sys

from veloxchem.veloxchemlib import is_mpi_master, mpi_barrier
from veloxchem.main import main


class TestInputECD:

    def create_input_file(self, lines, fname):

        if is_mpi_master():
            with open(fname, 'w') as f_inp:
                for line in lines.splitlines():
                    if line:
                        f_inp.write(f'{line}\n')
        mpi_barrier()

    def get_input_lines(self, xcfun, tamm_dancoff):

        input_lines = """
            @jobs
            task: response
            @end

            @method settings
            basis: sto-3g
            @end

            @scf
            restart: no
            @end

            @response
            restart: no
            property: absorption
            tamm_dancoff: yes
            nstates: 10
            @end

            @molecule
            charge: 0
            multiplicity: 1
            units: au
            xyz:
            C    -5.14663     0.73397    -0.06055
            C    -2.26501     0.55256    -0.00370
            H    -5.97431     0.00406     1.73252
            H    -5.88336    -0.47160    -1.61196
            O    -5.91711     3.23056    -0.50121
            H    -1.66302    -1.44978     0.23531
            O    -1.26657     2.04852     1.94120
            H    -1.48472     1.21997    -1.83424
            H    -1.55201     1.11020     3.54065
            H    -5.88106     4.10897     1.15610
            @end
        """

        if xcfun is not None:
            input_lines = input_lines.replace(
                '@method settings', f'@method settings\nxcfun:{xcfun}')

        if tamm_dancoff:
            return input_lines
        else:
            return input_lines.replace('tamm_dancoff: yes', 'tamm_dancoff: no')

    def run_input_ecd(self, capsys, xcfun, tamm_dancoff, ref_data):

        with tempfile.TemporaryDirectory() as temp_dir:
            input_file = Path(temp_dir, 'vlx_ecd')

            input_lines = self.get_input_lines(xcfun, tamm_dancoff)
            self.create_input_file(input_lines, input_file)

            with patch.object(sys, 'argv', ['vlx', str(input_file)]):
                main()
                captured = capsys.readouterr()

                if is_mpi_master():
                    exc_ene = []
                    osc_str = []
                    rot_str = []

                    for line in captured.out.splitlines():
                        if 'Wave Function Model' in line:
                            wf_model = line.split(':')[1].strip()

                        if 'Total Energy' in line:
                            e_scf = float(line.split()[3])

                        if 'SCF converged in' in line:
                            n_scf_iterations = int(line.split()[4])

                        if 'Debye' in line and 'Total' in line:
                            dipole = float(line.split()[2])

                        if 'Number of States' in line:
                            nstates = int(line.split(':')[1].strip())

                        if 'converged in' in line and 'SCF' not in line:
                            n_ecd_iterations = int(
                                line.split('converged in')[1].split()[0])

                        if 'Osc.Str.' in line:
                            content = line.split()
                            exc_ene.append(float(content[3]))
                            osc_str.append(float(content[-1]))

                        if 'Rot.Str.' in line:
                            content = line.split()
                            rot_str.append(float(content[4]))

                    exc_ene = np.array(exc_ene)
                    osc_str = np.array(osc_str)
                    rot_str = np.array(rot_str)

                    assert wf_model == ref_data['wf_model']
                    assert abs(e_scf - ref_data['e_scf']) < 1.0e-6
                    assert n_scf_iterations == ref_data['n_scf_iterations']
                    assert abs(dipole - ref_data['dipole']) < 1.0e-5

                    assert nstates == ref_data['nstates']
                    assert n_ecd_iterations == ref_data['n_ecd_iterations']
                    assert np.max(
                        np.abs(exc_ene - ref_data['exc_ene'])) < 1.0e-5
                    assert np.max(
                        np.abs(osc_str - ref_data['osc_str'])) < 1.0e-4
                    assert np.max(
                        np.abs(rot_str - ref_data['rot_str'])) < 2.0e-4

    def test_input_rhf_cis_ecd(self, capsys):

        # vlxtag: RHF, Absorption, ECD, CIS

        xcfun = None
        tamm_dancoff = True

        ref_data = {
            'wf_model': 'Spin-Restricted Hartree-Fock',
            'e_scf': -225.9526840315,
            'n_scf_iterations': 5,
            'dipole': 0.923876,
            'nstates': 10,
            'n_ecd_iterations': 10,
            'exc_ene': np.array([
                0.44968597, 0.45609008, 0.51132433, 0.52111186, 0.57126084,
                0.57312116, 0.63501338, 0.63805833, 0.64491553, 0.66005954
            ]),
            'osc_str': np.array([
                0.0067, 0.0076, 0.0043, 0.0009, 0.0262, 0.0179, 0.0908, 0.1243,
                0.0770, 0.1144
            ]),
            'rot_str': np.array([
                -0.008829, 0.009472, -0.006779, 0.006590, -0.092286, 0.063381,
                0.083296, -0.060358, -0.007365, -0.008269
            ]),
        }

        self.run_input_ecd(capsys, xcfun, tamm_dancoff, ref_data)

    def test_input_rhf_tdhf_ecd(self, capsys):

        # vlxtag: RHF, Absorption, ECD, TDHF

        xcfun = None
        tamm_dancoff = False

        ref_data = {
            'wf_model': 'Spin-Restricted Hartree-Fock',
            'e_scf': -225.9526840315,
            'n_scf_iterations': 5,
            'dipole': 0.923876,
            'nstates': 10,
            'n_ecd_iterations': 10,
            'exc_ene': np.array([
                0.44794105, 0.45442484, 0.51029301, 0.51999148, 0.56787603,
                0.56968528, 0.63307749, 0.63615325, 0.64390716, 0.65927314
            ]),
            'osc_str': np.array([
                0.0055, 0.0064, 0.0037, 0.0008, 0.0224, 0.0158, 0.0869, 0.1295,
                0.0629, 0.1133
            ]),
            'rot_str': np.array([
                -0.008410, 0.009532, -0.009151, 0.009084, -0.114150, 0.084127,
                0.095483, -0.075793, -0.001168, -0.007853
            ]),
        }

        self.run_input_ecd(capsys, xcfun, tamm_dancoff, ref_data)

    def test_input_rks_tda_ecd(self, capsys):

        # vlxtag: RKS, Absorption, ECD, TDA

        xcfun = 'b3lyp'
        tamm_dancoff = True

        ref_data = {
            'wf_model': 'Spin-Restricted Kohn-Sham',
            'e_scf': -227.1612519218,
            'n_scf_iterations': 8,
            'dipole': 0.859299,
            'nstates': 10,
            'n_ecd_iterations': 6,
            'exc_ene': np.array([
                0.36834945, 0.39181586, 0.41050739, 0.43685048, 0.43805628,
                0.46396967, 0.46892415, 0.47445040, 0.47872523, 0.49959807
            ]),
            'osc_str': np.array([
                0.0143, 0.0114, 0.0047, 0.0009, 0.0108, 0.0044, 0.0068, 0.0149,
                0.0057, 0.0370
            ]),
            'rot_str': np.array([
                -0.004653, 0.000543, -0.002207, 0.008626, 0.006254, 0.002287,
                -0.045490, -0.011791, 0.012514, 0.008681
            ]),
        }

        self.run_input_ecd(capsys, xcfun, tamm_dancoff, ref_data)

    def test_input_rks_tddft_ecd(self, capsys):

        # vlxtag: RKS, Absorption, ECD, TDDFT

        xcfun = 'b3lyp'
        tamm_dancoff = False

        ref_data = {
            'wf_model': 'Spin-Restricted Kohn-Sham',
            'e_scf': -227.1612519218,
            'n_scf_iterations': 8,
            'dipole': 0.859299,
            'nstates': 10,
            'n_ecd_iterations': 7,
            'exc_ene': np.array([
                0.36688036, 0.39023443, 0.41016395, 0.43625363, 0.43776105,
                0.46364009, 0.46676458, 0.47242464, 0.47787214, 0.49859230
            ]),
            'osc_str': np.array([
                0.0115, 0.0093, 0.0045, 0.0009, 0.0094, 0.0037, 0.0063, 0.0138,
                0.0052, 0.0329
            ]),
            'rot_str': np.array([
                -0.000777, 0.000611, -0.004004, 0.013441, 0.002076, 0.004832,
                -0.073971, 0.004679, 0.011049, 0.019718
            ]),
        }

        self.run_input_ecd(capsys, xcfun, tamm_dancoff, ref_data)
