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
class TestInputUVVis:

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
            basis: aug-cc-pvdz
            ri: false
            @end

            @scf
            restart: no
            @end

            @response
            restart: no
            property: absorption
            tamm_dancoff: yes
            nstates: 5
            @end

            @molecule
            charge: 0
            multiplicity: 1
            units: au
            xyz:
            O   0.0   0.0   0.0
            H   0.0   1.4   1.1
            H   0.0  -1.4   1.1
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

    def run_input_uvvis(self, capsys, xcfun, tamm_dancoff, ri, ref_data):

        here = Path(__file__).parent
        random_string = get_random_string_parallel()
        input_file = here / 'data' / f'vlx_uvvis_{random_string}.inp'

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
                        n_uvvis_iterations = int(
                            line.split('converged in')[1].split()[0])

                    if 'Osc.Str.' in line:
                        content = line.split()
                        exc_ene.append(float(content[3]))
                        osc_str.append(float(content[-1]))

                exc_ene = np.array(exc_ene)
                osc_str = np.array(osc_str)

                assert wf_model == ref_data['wf_model']
                assert abs(e_scf - ref_data['e_scf']) < 1.0e-6
                assert n_scf_iterations == ref_data['n_scf_iterations']
                assert abs(dipole - ref_data['dipole']) < 1.0e-5

                assert nstates == ref_data['nstates']
                assert n_uvvis_iterations == ref_data['n_uvvis_iterations']
                assert np.max(np.abs(exc_ene - ref_data['exc_ene'])) < 1.0e-5
                assert np.max(np.abs(osc_str - ref_data['osc_str'])) < 1.0e-4

            self.remove_input_and_h5_files(input_file)

    @pytest.mark.parametrize('ri', [False, True])
    def test_input_rhf_cis_uvvis(self, capsys, ri):

        # vlxtag: RHF, Absorption, CIS

        xcfun = None
        tamm_dancoff = True

        ref_data = {
            False: {
                'wf_model': 'Spin-Restricted Hartree-Fock',
                'e_scf': -76.0416975498,
                'n_scf_iterations': 9,
                'dipole': 0.786770,
                'nstates': 5,
                'n_uvvis_iterations': 8,
                'exc_ene': np.array([
                    0.32282999, 0.38491251, 0.41024825, 0.44793326, 0.47139996
                ]),
                'osc_str': np.array([0.0530, 0.0000, 0.1070, 0.0045, 0.0248]),
            },
            True: {
                'wf_model': 'Spin-Restricted Hartree-Fock',
                'e_scf': -76.041672470,
                'n_scf_iterations': 10,
                'dipole': 0.78679500,
                'nstates': 5,
                'n_uvvis_iterations': 8,
                'exc_ene': np.array([
                    0.32283670, 0.38496483, 0.41026026, 0.44791818, 0.47149378
                ]),
                'osc_str': np.array(
                    [0.0531000, 0.0000000, 0.1071000, 0.0046000, 0.0249000])
            }
        }[ri]

        self.run_input_uvvis(capsys, xcfun, tamm_dancoff, ri, ref_data)

    @pytest.mark.parametrize('ri', [False, True])
    def test_input_rhf_tdhf_uvvis(self, capsys, ri):

        # vlxtag: RHF, Absorption, TDHF

        xcfun = None
        tamm_dancoff = False

        ref_data = {
            False: {
                'wf_model': 'Spin-Restricted Hartree-Fock',
                'e_scf': -76.0416975498,
                'n_scf_iterations': 9,
                'dipole': 0.786770,
                'nstates': 5,
                'n_uvvis_iterations': 9,
                'exc_ene': np.array([
                    0.32139517, 0.38339622, 0.40934426, 0.44665234, 0.47004797
                ]),
                'osc_str': np.array([0.0519, 0.0000, 0.1018, 0.0048, 0.0233]),
            },
            True: {
                'wf_model': 'Spin-Restricted Hartree-Fock',
                'e_scf': -76.041672470,
                'n_scf_iterations': 10,
                'dipole': 0.78679500,
                'nstates': 5,
                'n_uvvis_iterations': 9,
                'exc_ene': np.array([
                    0.32140222, 0.38344740, 0.40935610, 0.44663890, 0.47013999
                ]),
                'osc_str': np.array(
                    [0.0519000, 0.0000000, 0.1018000, 0.0048000, 0.0233000])
            }
        }[ri]

        self.run_input_uvvis(capsys, xcfun, tamm_dancoff, ri, ref_data)

    @pytest.mark.parametrize('ri', [False, True])
    def test_input_rks_tda_uvvis(self, capsys, ri):

        # vlxtag: RKS, Absorption, TDA

        xcfun = 'b3lyp'
        tamm_dancoff = True

        ref_data = {
            False: {
                'wf_model': 'Spin-Restricted Kohn-Sham',
                'e_scf': -76.4435484737,
                'n_scf_iterations': 7,
                'dipole': 0.731257,
                'nstates': 5,
                'n_uvvis_iterations': 4,
                'exc_ene': np.array([
                    0.25720123, 0.30977192, 0.33859452, 0.37854004, 0.39054751
                ]),
                'osc_str': np.array([0.0537, 0.0000, 0.0906, 0.0000, 0.0127]),
            },
            True: {
                'wf_model': 'Spin-Restricted Kohn-Sham',
                'e_scf': -76.443575505,
                'n_scf_iterations': 7,
                'dipole': 0.73125400,
                'nstates': 5,
                'n_uvvis_iterations': 4,
                'exc_ene': np.array([
                    0.25720088, 0.30981356, 0.33860211, 0.37846444, 0.39064615
                ]),
                'osc_str': np.array(
                    [0.0537000, 0.0000000, 0.0906000, 0.0000000, 0.0127000])
            }
        }[ri]

        self.run_input_uvvis(capsys, xcfun, tamm_dancoff, ri, ref_data)

    @pytest.mark.parametrize('ri', [False, True])
    def test_input_rks_tddft_uvvis(self, capsys, ri):

        # vlxtag: RKS, Absorption, TDDFT

        xcfun = 'b3lyp'
        tamm_dancoff = False

        ref_data = {
            False: {
                'wf_model': 'Spin-Restricted Kohn-Sham',
                'e_scf': -76.4435484737,
                'n_scf_iterations': 7,
                'dipole': 0.731257,
                'nstates': 5,
                'n_uvvis_iterations': 5,
                'exc_ene': np.array([
                    0.25675266, 0.30971169, 0.33794203, 0.37844329, 0.39025586
                ]),
                'osc_str': np.array([0.0520, 0.0000, 0.0848, 0.0000, 0.0116]),
            },
            True: {
                'wf_model': 'Spin-Restricted Kohn-Sham',
                'e_scf': -76.443575505,
                'n_scf_iterations': 7,
                'dipole': 0.73125400,
                'nstates': 5,
                'n_uvvis_iterations': 5,
                'exc_ene': np.array([
                    0.25675235, 0.30975282, 0.33794984, 0.37836769, 0.39035400
                ]),
                'osc_str': np.array(
                    [0.0520000, 0.0000000, 0.0848000, 0.0000000, 0.0116000])
            }
        }[ri]

        self.run_input_uvvis(capsys, xcfun, tamm_dancoff, ri, ref_data)
