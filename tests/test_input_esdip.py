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
class TestInputExcitedStateDipole:

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

            rsp_h5 = Path(str(input_file).replace('.inp', '_orbrsp.h5'))
            if rsp_h5.is_file():
                rsp_h5.unlink()

        MPI.COMM_WORLD.barrier()

    def get_input_lines(self, xcfun, tamm_dancoff):

        input_lines = """
            @jobs
            task: gradient
            @end

            @method settings
            basis: def2-svp
            @end

            @gradient
            state_deriv_index: 3
            do_first_order_prop: yes
            @end

            @response
            restart: no
            property: absorption
            nstates: 5
            @end

            @molecule
            charge: 0
            multiplicity: 1
            xyz:
            O      0.0000        0.0000        0.0000
            H      0.0000        0.7408        0.5821
            H      0.0000       -0.7408        0.5821
            @end
        """

        if xcfun is not None:
            input_lines = input_lines.replace(
                '@method settings', f'@method settings\nxcfun:{xcfun}')

        if tamm_dancoff:
            return input_lines
        else:
            return input_lines.replace('tamm_dancoff: yes', 'tamm_dancoff: no')

    def run_input_esdip(self, capsys, xcfun, tamm_dancoff, ref_data):

        here = Path(__file__).parent
        random_string = get_random_string_parallel()
        input_file = here / 'data' / f'vlx_esdip_{random_string}.inp'

        input_lines = self.get_input_lines(xcfun, tamm_dancoff)
        self.create_input_file(input_lines, input_file)

        with patch.object(sys, 'argv', ['vlx', str(input_file)]):
            main()
            captured = capsys.readouterr()

            if MPI.COMM_WORLD.Get_rank() == mpi_master():
                """
                                RPA Unrelaxed Dipole Moment(s)
                               ---------------------------------


                                        Excited State 3

                      X   :        -0.000000 a.u.        -0.000000 Debye
                      Y   :         0.000000 a.u.         0.000000 Debye
                      Z   :        -0.837253 a.u.        -2.128085 Debye
                    Total :         0.837253 a.u.         2.128085 Debye


                                 RPA Relaxed Dipole Moment(s)
                                -------------------------------


                                        Excited State 3

                      X   :        -0.000000 a.u.        -0.000000 Debye
                      Y   :         0.000000 a.u.         0.000000 Debye
                      Z   :        -0.406896 a.u.        -1.034225 Debye
                    Total :         0.406896 a.u.         1.034225 Debye
                """

                es_dip_unrelaxed = []
                es_dip_relaxed = []

                read_dip_unrelaxed = False
                read_dip_relaxed = False

                for line in captured.out.splitlines():

                    if line.strip() == 'RPA Unrelaxed Dipole Moment(s)':
                        read_dip_unrelaxed = True
                    elif line.strip() == 'RPA Relaxed Dipole Moment(s)':
                        read_dip_relaxed = True

                    if read_dip_unrelaxed and (':' in line):
                        content = line.split()
                        if content[0] in ['X', 'Y', 'Z'] and content[1] == ':':
                            es_dip_unrelaxed.append(float(content[2]))
                        elif content[0] == 'Total' and content[1] == ':':
                            read_dip_unrelaxed = False

                    if read_dip_relaxed and (':' in line):
                        content = line.split()
                        if content[0] in ['X', 'Y', 'Z'] and content[1] == ':':
                            es_dip_relaxed.append(float(content[2]))
                        elif content[0] == 'Total' and content[1] == ':':
                            read_dip_relaxed = False

                es_dip_unrelaxed = np.array(es_dip_unrelaxed)
                es_dip_relaxed = np.array(es_dip_relaxed)

                assert np.max(
                    np.abs(es_dip_unrelaxed -
                           ref_data['es_dip_unrelaxed'])) < 1.0e-4
                assert np.max(
                    np.abs(es_dip_relaxed -
                           ref_data['es_dip_relaxed'])) < 1.0e-4

            self.remove_input_and_h5_files(input_file)

    def test_input_rks_tddft_esdip(self, capsys):

        # vlxtag: RKS, Absorption, TDDFT

        xcfun = 'pbe0'
        tamm_dancoff = False

        ref_data = {
            'es_dip_unrelaxed': np.array([0.0, 0.0, -0.837253]),
            'es_dip_relaxed': np.array([0.0, 0.0, -0.406896]),
        }

        self.run_input_esdip(capsys, xcfun, tamm_dancoff, ref_data)
