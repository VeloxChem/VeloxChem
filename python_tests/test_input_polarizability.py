from pathlib import Path
from unittest.mock import patch
import numpy as np
import tempfile
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

        with tempfile.TemporaryDirectory() as temp_dir:
            input_file = Path(temp_dir, 'vlx_polarizability')

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

    def test_input_rhf_lr_polarizability(self, capsys):

        # vlxtag: RHF, LR, Polarizability

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

        # vlxtag: RKS, LR, Polarizability

        xcfun = 'b3lyp'

        ref_data = {
            'wf_model': 'Spin-Restricted Kohn-Sham',
            'e_scf': -76.4435551612,
            'n_scf_iterations': 7,
            'n_rsp_iterations': 5,
            'polarizability': np.array([
                8.769176, 0, 0, 0, 9.696176, 0, 0, 0, 9.066348, 8.893308, 0, 0,
                0, 9.757660, 0, 0, 0, 9.147369, 9.318430, 0, 0, 0, 9.948384, 0,
                0, 0, 9.407202
            ])
        }

        self.run_input_polarizability(capsys, xcfun, ref_data)
