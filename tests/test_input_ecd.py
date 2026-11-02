from mpi4py import MPI
from pathlib import Path
from unittest.mock import patch
import numpy as np
import pytest
import sys

from veloxchem.veloxchemlib import mpi_master
from veloxchem.inputparser import get_random_string_parallel
from veloxchem.main import main


@pytest.mark.solvers
class TestInputECD:

    def create_input_file(self, lines, fname):

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            with open(fname, 'w') as f_inp:
                for line in lines.splitlines():
                    if line:
                        f_inp.write(f'{line}\n')

        MPI.COMM_WORLD.barrier()

    def remove_input_and_h5_files(self, input_file):

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            if input_file.is_file():
                input_file.unlink()

            final_h5 = Path(str(input_file).replace('.inp', '.h5'))
            if final_h5.is_file():
                final_h5.unlink()

            scf_h5 = Path(str(input_file).replace('.inp', '_scf.h5'))
            if scf_h5.is_file():
                scf_h5.unlink()

            rsp_h5 = Path(str(input_file).replace('.inp', '_rsp.h5'))
            if rsp_h5.is_file():
                rsp_h5.unlink()

        MPI.COMM_WORLD.barrier()

    def get_input_lines(self, xcfun, tamm_dancoff, ri):

        input_lines = """
            @jobs
            task: response
            @end

            @method settings
            basis: sto-3g
            ri: false
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

        if ri:
            input_lines = input_lines.replace('ri: false', 'ri: yes')

        if xcfun is not None:
            input_lines = input_lines.replace(
                '@method settings', f'@method settings\nxcfun:{xcfun}')

        if tamm_dancoff:
            return input_lines
        else:
            return input_lines.replace('tamm_dancoff: yes', 'tamm_dancoff: no')

    def run_input_ecd(self, capsys, xcfun, tamm_dancoff, ri, ref_data):

        here = Path(__file__).parent
        random_string = get_random_string_parallel()
        input_file = here / 'data' / f'vlx_ecd_{random_string}.inp'

        input_lines = self.get_input_lines(xcfun, tamm_dancoff, ri)
        self.create_input_file(input_lines, input_file)

        with patch.object(sys, 'argv', ['vlx', str(input_file)]):
            main()
            captured = capsys.readouterr()

            if MPI.COMM_WORLD.Get_rank() == mpi_master():
                # NOTE: printed only where the SCF formed its Fock matrices
                # with RI; one which fell back to the four-center integrals
                # prints a warning instead.
                if ri:
                    assert ('Using the SIMD resolution of the identity'
                            in captured.out)

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
                assert np.max(np.abs(exc_ene - ref_data['exc_ene'])) < 1.0e-5
                assert np.max(np.abs(osc_str - ref_data['osc_str'])) < 1.0e-4
                assert np.max(np.abs(rot_str - ref_data['rot_str'])) < 2.0e-4

            self.remove_input_and_h5_files(input_file)

    @pytest.mark.parametrize('ri', [False, True])
    def test_input_rhf_cis_ecd(self, capsys, ri):

        # vlxtag: RHF, Absorption, ECD, CIS

        xcfun = None
        tamm_dancoff = True

        ref_data = {
            False: {
                'wf_model': 'Spin-Restricted Hartree-Fock',
                'e_scf': -225.9526840315,
                'n_scf_iterations': 5,
                'dipole': 0.923876,
                'nstates': 10,
                'n_ecd_iterations': 8,
                'exc_ene': np.array([
                    0.44968597, 0.45609008, 0.51132433, 0.52111186, 0.57126084,
                    0.57312116, 0.63501338, 0.63805833, 0.64491553, 0.66005954
                ]),
                'osc_str': np.array([
                    0.0067, 0.0076, 0.0043, 0.0009, 0.0262, 0.0179, 0.0908,
                    0.1243, 0.0770, 0.1143
                ]),
                'rot_str': np.array([
                    -0.008829, 0.009472, -0.006779, 0.006590, -0.092286,
                    0.063381, 0.083296, -0.060358, -0.007365, -0.008269
                ]),
            },
            True: {
                'wf_model': 'Spin-Restricted Hartree-Fock',
                'e_scf': -225.953068460,
                'n_scf_iterations': 10,
                'dipole': 0.92367500,
                'nstates': 10,
                'n_ecd_iterations': 8,
                'exc_ene': np.array([
                    0.44969065, 0.45609353, 0.51134184, 0.52112819, 0.57127885,
                    0.57313854, 0.63503481, 0.63807742, 0.64493047, 0.66007690
                ]),
                'osc_str': np.array([
                    0.0067000, 0.0076000, 0.0043000, 0.0009000, 0.0262000,
                    0.0179000, 0.0908000, 0.1242000, 0.0771000, 0.1144000
                ]),
                'rot_str': np.array([
                    -0.0088270, 0.0094720, -0.0067810, 0.0065920, -0.0923000,
                    0.0633940, 0.0833060, -0.0603300, -0.0074290, -0.0082800
                ])
            }
        }[ri]

        self.run_input_ecd(capsys, xcfun, tamm_dancoff, ri, ref_data)

    @pytest.mark.parametrize('ri', [False, True])
    def test_input_rhf_tdhf_ecd(self, capsys, ri):

        # vlxtag: RHF, Absorption, ECD, TDHF

        xcfun = None
        tamm_dancoff = False

        ref_data = {
            False: {
                'wf_model': 'Spin-Restricted Hartree-Fock',
                'e_scf': -225.9526840315,
                'n_scf_iterations': 5,
                'dipole': 0.923876,
                'nstates': 10,
                'n_ecd_iterations': 9,
                'exc_ene': np.array([
                    0.44794105, 0.45442484, 0.51029301, 0.51999148, 0.56787603,
                    0.56968528, 0.63307749, 0.63615325, 0.64390716, 0.65927314
                ]),
                'osc_str': np.array([
                    0.0055, 0.0064, 0.0037, 0.0008, 0.0224, 0.0158, 0.0869,
                    0.1295, 0.0629, 0.1133
                ]),
                'rot_str': np.array([
                    -0.008410, 0.009532, -0.009151, 0.009084, -0.114150,
                    0.084127, 0.095483, -0.075793, -0.001168, -0.007853
                ]),
            },
            True: {
                'wf_model': 'Spin-Restricted Hartree-Fock',
                'e_scf': -225.953068460,
                'n_scf_iterations': 10,
                'dipole': 0.92367500,
                'nstates': 10,
                'n_ecd_iterations': 9,
                'exc_ene': np.array([
                    0.44794631, 0.45442890, 0.51031077, 0.52000815, 0.56789404,
                    0.56970248, 0.63310001, 0.63617385, 0.64392115, 0.65929056
                ]),
                'osc_str': np.array([
                    0.0055000, 0.0064000, 0.0037000, 0.0008000, 0.0224000,
                    0.0158000, 0.0869000, 0.1294000, 0.0630000, 0.1134000
                ]),
                'rot_str': np.array([
                    -0.0084090, 0.0095330, -0.0091520, 0.0090870, -0.1141620,
                    0.0841410, 0.0954940, -0.0757740, -0.0012220, -0.0078570
                ])
            }
        }[ri]

        self.run_input_ecd(capsys, xcfun, tamm_dancoff, ri, ref_data)

    @pytest.mark.parametrize('ri', [False, True])
    def test_input_rks_tda_ecd(self, capsys, ri):

        # vlxtag: RKS, Absorption, ECD, TDA

        xcfun = 'b3lyp'
        tamm_dancoff = True

        ref_data = {
            False: {
                'wf_model': 'Spin-Restricted Kohn-Sham',
                'e_scf': -227.1612519218,
                'n_scf_iterations': 4,
                'dipole': 0.859299,
                'nstates': 10,
                'n_ecd_iterations': 5,
                'exc_ene': np.array([
                    0.36834945, 0.39181586, 0.41050739, 0.43685048, 0.43805628,
                    0.46396967, 0.46892415, 0.47445040, 0.47872523, 0.49959807
                ]),
                'osc_str': np.array([
                    0.0143, 0.0114, 0.0047, 0.0009, 0.0108, 0.0044, 0.0068,
                    0.0149, 0.0057, 0.0370
                ]),
                'rot_str': np.array([
                    -0.004653, 0.000543, -0.002207, 0.008626, 0.006254,
                    0.002287, -0.045490, -0.011791, 0.012514, 0.008681
                ]),
            },
            True: {
                'wf_model': 'Spin-Restricted Kohn-Sham',
                'e_scf': -227.162522904,
                'n_scf_iterations': 10,
                'dipole': 0.85913600,
                'nstates': 10,
                'n_ecd_iterations': 5,
                'exc_ene': np.array([
                    0.36835089, 0.39181022, 0.41050957, 0.43684377, 0.43805982,
                    0.46397385, 0.46892465, 0.47444941, 0.47871727, 0.49959757
                ]),
                'osc_str': np.array([
                    0.0143000, 0.0114000, 0.0047000, 0.0009000, 0.0108000,
                    0.0044000, 0.0068000, 0.0150000, 0.0057000, 0.0371000
                ]),
                'rot_str': np.array([
                    -0.0046530, 0.0005420, -0.0022060, 0.0086690, 0.0062190,
                    0.0022860, -0.0454850, -0.0118330, 0.0125510, 0.0087110
                ])
            }
        }[ri]

        self.run_input_ecd(capsys, xcfun, tamm_dancoff, ri, ref_data)

    @pytest.mark.parametrize('ri', [False, True])
    def test_input_rks_tddft_ecd(self, capsys, ri):

        # vlxtag: RKS, Absorption, ECD, TDDFT

        xcfun = 'b3lyp'
        tamm_dancoff = False

        ref_data = {
            False: {
                'wf_model': 'Spin-Restricted Kohn-Sham',
                'e_scf': -227.1612519218,
                'n_scf_iterations': 4,
                'dipole': 0.859299,
                'nstates': 10,
                'n_ecd_iterations': 5,
                'exc_ene': np.array([
                    0.36688036, 0.39023443, 0.41016395, 0.43625363, 0.43776105,
                    0.46364009, 0.46676458, 0.47242464, 0.47787214, 0.49859230
                ]),
                'osc_str': np.array([
                    0.0115, 0.0093, 0.0045, 0.0009, 0.0094, 0.0037, 0.0063,
                    0.0138, 0.0052, 0.0329
                ]),
                'rot_str': np.array([
                    -0.000777, 0.000611, -0.004004, 0.013441, 0.002076,
                    0.004832, -0.073971, 0.004679, 0.011049, 0.019718
                ]),
            },
            True: {
                'wf_model': 'Spin-Restricted Kohn-Sham',
                'e_scf': -227.162522904,
                'n_scf_iterations': 10,
                'dipole': 0.85913600,
                'nstates': 10,
                'n_ecd_iterations': 5,
                'exc_ene': np.array([
                    0.36688276, 0.39023009, 0.41016640, 0.43624742, 0.43776504,
                    0.46364485, 0.46676575, 0.47242590, 0.47786234, 0.49859259
                ]),
                'osc_str': np.array([
                    0.0115000, 0.0093000, 0.0045000, 0.0009000, 0.0094000,
                    0.0037000, 0.0063000, 0.0138000, 0.0052000, 0.0329000
                ]),
                'rot_str': np.array([
                    -0.0007780, 0.0006100, -0.0040040, 0.0134840, 0.0020450,
                    0.0048310, -0.0739600, 0.0046400, 0.0110830, 0.0197490
                ])
            }
        }[ri]

        self.run_input_ecd(capsys, xcfun, tamm_dancoff, ri, ref_data)
