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
class TestInputCppECD:

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

    def get_input_lines(self, xcfun):

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
            property: ecd (cpp)
            frequencies: 0.45-0.55(0.01)
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

        return input_lines

    def run_input_cpp_ecd(self, capsys, xcfun, ref_data):

        here = Path(__file__).parent
        random_string = get_random_string_parallel()
        input_file = here / 'data' / f'vlx_cpp_ecd_{random_string}.inp'

        input_lines = self.get_input_lines(xcfun)
        self.create_input_file(input_lines, input_file)

        with patch.object(sys, 'argv', ['vlx', str(input_file)]):
            main()
            captured = capsys.readouterr()

            if MPI.COMM_WORLD.Get_rank() == mpi_master():
                cpp_line_number = None
                output_lines = captured.out.splitlines()

                for line_idx, line in enumerate(output_lines):
                    if 'Wave Function Model' in line:
                        wf_model = line.split(':')[1].strip()

                    if 'Total Energy' in line:
                        e_scf = float(line.split()[3])

                    if 'SCF converged in' in line:
                        n_scf_iterations = int(line.split()[4])

                    if 'Debye' in line and 'Total' in line:
                        dipole = float(line.split()[2])

                    if 'Number of Frequencies' in line:
                        nfreqs = int(line.split(':')[1].strip())

                    if 'converged in' in line and 'SCF' not in line:
                        n_cpp_iterations = int(
                            line.split('converged in')[1].split()[0])

                    if 'Delta_epsilon[L mol^-1 cm^-1]' in line:
                        cpp_line_number = line_idx

                cpp_energies = []
                cpp_ecd_values = []

                for line in output_lines[cpp_line_number + 2:cpp_line_number +
                                         2 + nfreqs]:
                    content = line.split()
                    cpp_energies.append(float(content[0]))
                    cpp_ecd_values.append(float(content[2]))

                cpp_energies = np.array(cpp_energies)
                cpp_ecd_values = np.array(cpp_ecd_values)

                assert wf_model == ref_data['wf_model']
                assert abs(e_scf - ref_data['e_scf']) < 1.0e-6
                assert n_scf_iterations == ref_data['n_scf_iterations']
                assert abs(dipole - ref_data['dipole']) < 1.0e-5

                assert nfreqs == ref_data['nfreqs']
                assert n_cpp_iterations == ref_data['n_cpp_iterations']
                assert np.max(np.abs(cpp_energies -
                                     ref_data['cpp_energies'])) < 1.0e-5
                assert np.max(
                    np.abs(cpp_ecd_values -
                           ref_data['cpp_ecd_values'])) < 1.0e-4

            self.remove_input_and_h5_files(input_file)

    def test_input_rhf_cpp_ecd(self, capsys):

        xcfun = None

        ref_data = {
            'wf_model': 'Spin-Restricted Hartree-Fock',
            'e_scf': -225.9526840315,
            'n_scf_iterations': 5,
            'dipole': 0.923876,
            'nfreqs': 11,
            'n_cpp_iterations': 10,
            'cpp_energies': np.array([
                0.4500, 0.4600, 0.4700, 0.4800, 0.4900, 0.5000, 0.5100, 0.5200,
                0.5300, 0.5400, 0.5500
            ]),
            'cpp_ecd_values': np.array([
                -1.29153530, 1.80271442, 0.24229833, 0.01949196, -0.14441290,
                -0.77377515, -5.55622119, 5.44284829, 0.55350044, -0.46225168,
                -1.87160197
            ]),
        }

        self.run_input_cpp_ecd(capsys, xcfun, ref_data)

    def test_input_rks_cpp_ecd(self, capsys):

        xcfun = 'b3lyp'

        ref_data = {
            'wf_model': 'Spin-Restricted Kohn-Sham',
            'e_scf': -227.1612519218,
            'n_scf_iterations': 8,
            'dipole': 0.859299,
            'nfreqs': 11,
            'n_cpp_iterations': 5,
            'cpp_energies': np.array([
                0.4500, 0.4600, 0.4700, 0.4800, 0.4900, 0.5000, 0.5100, 0.5200,
                0.5300, 0.5400, 0.5500
            ]),
            'cpp_ecd_values': np.array([
                -1.5111483, -11.89993278, -26.60109317, 3.54769273, 4.04898144,
                20.54478343, 3.95404303, -0.75578452, 8.35512241, -9.02677093,
                21.08438121
            ]),
        }

        self.run_input_cpp_ecd(capsys, xcfun, ref_data)
