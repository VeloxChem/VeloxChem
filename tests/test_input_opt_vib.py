from mpi4py import MPI
from pathlib import Path
from unittest.mock import patch
import numpy as np
import pytest
import sys

from veloxchem.veloxchemlib import mpi_master
from veloxchem.main import main


@pytest.mark.solvers
class TestInputOptVib:

    def create_input_file(self, lines, fname):

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            with open(fname, 'w') as f_inp:
                for line in lines.splitlines():
                    if line:
                        f_inp.write(f'{line}\n')

        MPI.COMM_WORLD.barrier()


    def test_input_opt_vib(self, capsys, tmp_path):

        # Note: bad geometry for testing purposes
        # `task: optimize` with `hessian: last` triggers optimization
        # followed by a vibrational analysis automatically

        input_lines = """
            @jobs
            task: optimize
            @end

            @method settings
            basis: def2-svp
            @end

            @optimize
            hessian: last
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

        input_file = tmp_path / 'vlx_opt_vib.inp'
        input_file = Path(MPI.COMM_WORLD.bcast(str(input_file), root=mpi_master()))

        self.create_input_file(input_lines, input_file)

        with patch.object(sys, 'argv', ['vlx', str(input_file)]):
            main()
            captured = capsys.readouterr()

            if MPI.COMM_WORLD.Get_rank() == mpi_master():
                freqs = []
                red_masses = []
                force_consts = []
                ir_ints = []

                for line in captured.out.splitlines():
                    if 'Harmonic frequency:' in line:
                        freqs.append(float(line.split(':')[1].split()[0]))
                    if 'Reduced mass:' in line:
                        red_masses.append(float(line.split(':')[1].split()[0]))
                    if 'Force constant:' in line:
                        force_consts.append(float(line.split(':')[1].split()[0]))
                    if 'IR intensity:' in line:
                        ir_ints.append(float(line.split(':')[1].split()[0]))

                # water has 3N-6 = 3 vibrational modes
                assert len(freqs) == 3
                assert len(red_masses) == 3
                assert len(force_consts) == 3
                assert len(ir_ints) == 3

                freqs = np.array(freqs)
                red_masses = np.array(red_masses)
                force_consts = np.array(force_consts)
                ir_ints = np.array(ir_ints)

                ref_freqs = np.array([1750.50, 4148.93, 4245.10])
                ref_masses = np.array([1.0816, 1.0461, 1.0825])
                ref_force_consts = np.array([1.9528, 10.6100, 11.4934])
                ref_ir_ints = np.array([82.7062, 26.8087, 74.4989])

                np.testing.assert_allclose(freqs, ref_freqs, atol=0.05)
                np.testing.assert_allclose(red_masses, ref_masses, atol=0.05)
                np.testing.assert_allclose(force_consts, ref_force_consts, atol=0.05)
                np.testing.assert_allclose(ir_ints, ref_ir_ints, atol=0.05)
