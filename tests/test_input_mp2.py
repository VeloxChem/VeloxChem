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
class TestInputMp2:

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

        MPI.COMM_WORLD.barrier()

    def get_input_lines(self):

        input_lines = """
            @jobs
            task: mp2
            @end

            @method settings
            basis: def2-svp
            @end

            @molecule
            charge: 0
            multiplicity: 1
            units: au
            xyz:
            Se  0.0   0.0   0.0
            H   0.0   0.0   2.8
            H   0.0   2.8   0.0
            @end
        """

        return input_lines

    def run_input_mp2(self, capsys, ref_data):

        here = Path(__file__).parent
        random_string = get_random_string_parallel()
        input_file = here / 'data' / f'vlx_mp2_{random_string}.inp'

        input_lines = self.get_input_lines()
        self.create_input_file(input_lines, input_file)

        with patch.object(sys, 'argv', ['vlx', str(input_file)]):
            main()
            captured = capsys.readouterr()

            if MPI.COMM_WORLD.Get_rank() == mpi_master():

                for line in captured.out.splitlines():

                    if 'MP2 correlation energy:' in line:
                        content = line.split()
                        mp2_energy = float(content[4])

                assert np.max(np.abs(mp2_energy -
                                     ref_data['mp2_energy'])) < 1.0e-8

            self.remove_input_and_h5_files(input_file)

    def test_input_rhf_mp2(self, capsys):

        ref_data = {
            'mp2_energy': -0.28529088,
        }

        self.run_input_mp2(capsys, ref_data)
