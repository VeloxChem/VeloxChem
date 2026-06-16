from mpi4py import MPI
from pathlib import Path
from unittest.mock import patch
import pytest
import sys

from veloxchem.veloxchemlib import mpi_master
from veloxchem.main import main


@pytest.mark.solvers
class TestInputScfPfon:

    def create_input_file(self, lines, fname):

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            with open(fname, 'w') as f_inp:
                for line in lines.splitlines():
                    if line:
                        f_inp.write(f'{line}\n')

        MPI.COMM_WORLD.barrier()

    def test_input_scf_pfon(self, capsys, tmp_path):

        # Note: bad geometry for testing purposes

        input_lines = """
            @jobs
            task: scf
            @end

            @method settings
            basis: def2-svp
            @end

            @scf
            restart: no
            pfon: yes
            conv_thresh: 1e-3
            acc_type: diis
            @end

            @molecule
            charge: 0
            multiplicity: 1
            units: au
            xyz:
            O   0.0   0.0   0.0
            H   0.0   1.1   1.1
            H   0.0  -1.1   1.1
            @end
        """

        input_file = tmp_path / 'vlx_scf_pfon.inp'
        input_file = Path(MPI.COMM_WORLD.bcast(str(input_file), root=mpi_master()))

        self.create_input_file(input_lines, input_file)

        with patch.object(sys, 'argv', ['vlx', str(input_file)]):
            main()
            captured = capsys.readouterr()

            if MPI.COMM_WORLD.Get_rank() == mpi_master():
                pfon_count = 0
                n_iterations = 0

                for line in captured.out.splitlines():
                    if 'Applying pseudo-FON' in line:
                        pfon_count += 1

                    if 'SCF converged in' in line:
                        n_iterations = int(line.split()[4])

                # pFON is applied 4 times over 6 iterations;
                # 4 < 6 confirms pFON disappears before convergence
                assert pfon_count == 4
                assert n_iterations == 6
