from pathlib import Path
from unittest.mock import patch
from random import choice
import numpy as np
import random
import sys

from veloxchem.veloxchemlib import is_mpi_master, mpi_barrier
from veloxchem.main import main


class TestInputPolarizability:

    def create_input_file(self, lines, fname):

        if is_mpi_master():
            with open(fname, 'w') as f_inp:
                for line in lines.splitlines():
                    if line:
                        f_inp.write(f'{line}\n')
        mpi_barrier()

    def remove_input_and_h5_files(self, input_file):

        if is_mpi_master():
            if input_file.is_file():
                input_file.unlink()

            scf_h5 = input_file.with_suffix('.scf.h5')
            if scf_h5.is_file():
                scf_h5.unlink()

            scf_final_h5 = scf_h5.with_suffix('.tensors.h5')
            if scf_final_h5.is_file():
                scf_final_h5.unlink()

            rsp_h5 = input_file.with_suffix('.rsp.h5')
            if rsp_h5.is_file():
                rsp_h5.unlink()

            rsp_solutions_h5 = rsp_h5.with_suffix('.solutions.h5')
            if rsp_solutions_h5.is_file():
                rsp_solutions_h5.unlink()

        mpi_barrier()

    def get_input_lines(self, xcfun):

        input_lines = """
            @jobs
            task: response
            @end

            @method settings
            basis: aug-cc-pvdz
            @end

            @scf
            restart: no
            @end

            @response
            restart: no
            property: polarizability
            frequencies: 0.0, 0.05, 0.1
            batch_size: ==batch_size==
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

        input_lines = input_lines.replace('==batch_size==',
                                          str(random.choice([1, 10, 100])))

        if xcfun is not None:
            input_lines = input_lines.replace(
                '@method settings', f'@method settings\nxcfun:{xcfun}')

        return input_lines

    def run_input_polarizability(self, capsys, xcfun, ref_data):

        here = Path(__file__).parent
        random_string = ''.join([choice('abcdef123456') for i in range(8)])
        input_file = here / 'inputs' / f'vlx_polarizability_{random_string}.inp'

        input_lines = self.get_input_lines(xcfun)
        self.create_input_file(input_lines, input_file)

        with patch.object(sys, 'argv', ['vlx', str(input_file)]):
            main()
            captured = capsys.readouterr()

            if is_mpi_master():
                polarizability = []
                polar_flag = False

                for line in captured.out.splitlines():
                    if 'Wave Function Model' in line:
                        wf_model = line.split(':')[1].strip()

                    if 'converged in' in line and 'SCF' not in line:
                        n_rsp_iterations = int(
                            line.split('converged in')[1].split()[0])

                    if line.split() == ['X', 'Y', 'Z']:
                        polar_flag = True

                    if polar_flag:
                        content = line.split()
                        if len(content) != 4:
                            continue
                        if content[0] in ['X', 'Y', 'Z']:
                            polarizability.append(float(content[1]))
                            polarizability.append(float(content[2]))
                            polarizability.append(float(content[3]))
                        if content[0] == 'Z':
                            polar_flag = False

                polarizability = np.array(polarizability)

                assert wf_model == ref_data['wf_model']
                assert n_rsp_iterations == ref_data['n_rsp_iterations']
                assert np.max(
                    np.abs(polarizability -
                           ref_data['polarizability'])) < 1.0e-4

            self.remove_input_and_h5_files(input_file)

    def test_input_rhf_lr_polarizability(self, capsys):

        # vlxtag: RHF, Polarizability, LR

        xcfun = None

        ref_data = {
            'wf_model': 'Spin-Restricted Hartree-Fock',
            'n_rsp_iterations': 7,
            'polarizability': np.array([
                7.251835, 0, 0, 0, 8.724516, 0, 0, 0, 7.880586, 7.311119, 0, 0,
                0, 8.770615, 0, 0, 0, 7.930012, 7.502485, 0, 0, 0, 8.912738, 0,
                0, 0, 8.084815
            ])
        }

        self.run_input_polarizability(capsys, xcfun, ref_data)

    def test_input_rks_lr_polarizability(self, capsys):

        # vlxtag: RKS, Polarizability, LR

        xcfun = 'b3lyp'

        ref_data = {
            'wf_model': 'Spin-Restricted Kohn-Sham',
            'e_scf': -76.4435484737,
            'n_scf_iterations': 7,
            'n_rsp_iterations': 5,
            'polarizability': np.array([
                8.769167, 0, 0, 0, 9.696169, 0, 0, 0, 9.066348, 8.893297, 0, 0,
                0, 9.757652, 0, 0, 0, 9.147367, 9.318413, 0, 0, 0, 9.948375, 0,
                0, 0, 9.407198
            ])
        }

        self.run_input_polarizability(capsys, xcfun, ref_data)
