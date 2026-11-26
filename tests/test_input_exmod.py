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
class TestInputExcitonModel:

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

            exmod_h5 = Path(str(input_file).replace('.inp', '_exciton.h5'))
            if exmod_h5.is_file():
                exmod_h5.unlink()

            for idx in [1, 2]:

                final_h5 = Path(
                    str(exmod_h5).replace('_exciton.h5',
                                          f'_exciton_monomer_{idx}.h5'))
                if final_h5.is_file():
                    final_h5.unlink()

                scf_h5 = Path(
                    str(exmod_h5).replace('_exciton.h5',
                                          f'_exciton_monomer_{idx}_scf.h5'))
                if scf_h5.is_file():
                    scf_h5.unlink()

                rsp_h5 = Path(
                    str(exmod_h5).replace('_exciton.h5',
                                          f'_exciton_monomer_{idx}_rsp.h5'))
                if rsp_h5.is_file():
                    rsp_h5.unlink()

        MPI.COMM_WORLD.barrier()

    def get_input_lines(self, xcfun, ri):

        input_lines = """
            @jobs
            task: exciton
            @end

            @method settings
            basis: def2-svp
            ri: false
            @end

            @exciton
            fragments: 2
            atoms_per_fragment: 6
            nstates: 5
            ct_nocc: 1
            ct_nvir: 1
            @end

            @molecule
            charge: 0
            multiplicity: 1
            xyz:
            C         -1.37731        1.01769       -0.71611
            C         -0.04211        1.07142       -0.72602
            H         -1.96225        1.74636       -0.16458
            H         -1.90859        0.23094       -1.24174
            H          0.49049        1.84498       -0.18262
            H          0.54315        0.32947       -1.25941
            C         -1.17537       -1.48468        2.37427
            C          0.06813       -1.06658        2.62697
            H         -1.35657       -2.40378        1.82687
            H          0.92893       -1.63558        2.29127
            H         -2.03527       -0.90348        2.69157
            H          0.24803       -0.13578        3.15527
            @end
        """

        if ri:
            input_lines = input_lines.replace('ri: false', 'ri: yes')

        if xcfun is not None:
            input_lines = input_lines.replace(
                '@method settings', f'@method settings\nxcfun:{xcfun}')

        return input_lines

    def run_input_exmod(self, capsys, xcfun, ri, ref_data):

        here = Path(__file__).parent
        random_string = get_random_string_parallel()
        input_file = here / 'data' / f'vlx_exmod_{random_string}.inp'

        input_lines = self.get_input_lines(xcfun, ri)
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

                is_exmod_section = False

                for line in captured.out.splitlines():

                    if 'Adiabatic excited states:' in line:
                        is_exmod_section = True

                    if is_exmod_section and 'Osc.Str.' in line:
                        content = line.split()
                        exc_ene.append(float(content[3]))
                        osc_str.append(float(content[-1]))

                    if is_exmod_section and 'Rot.Str.' in line:
                        content = line.split()
                        rot_str.append(float(content[4]))

                exc_ene = np.array(exc_ene)
                osc_str = np.array(osc_str)
                rot_str = np.array(rot_str)

                assert np.max(np.abs(exc_ene - ref_data['exc_ene'])) < 1.0e-5
                assert np.max(np.abs(osc_str - ref_data['osc_str'])) < 1.0e-4
                assert np.max(np.abs(rot_str - ref_data['rot_str'])) < 2.0e-4

            self.remove_input_and_h5_files(input_file)

    @pytest.mark.parametrize('ri', [False, True])
    def test_input_rhf_exmod(self, capsys, ri):

        xcfun = None

        ref_data = {
            False: {
                'exc_ene': np.array([
                    0.28932589, 0.31669640, 0.34117880, 0.34154634, 0.34448158,
                    0.34474205, 0.35545359, 0.35550180, 0.37509079, 0.37587732,
                    0.40570065, 0.41752625
                ]),
                'osc_str': np.array([
                    0.031064, 1.241580, 0.000108, 0.013926, 0.000704, 0.000126,
                    0.000001, 0.000025, 0.000134, 0.001633, 0.024283, 0.005320
                ]),
                'rot_str': np.array([
                    0.149686, -0.157138, -0.000080, 0.000033, 0.021327,
                    -0.001754, -0.000358, -0.002835, -0.010791, -0.007212,
                    0.001252, 0.016260
                ]),
            },
            True: {
                'exc_ene': np.array([
                    0.28934303, 0.31671347, 0.34118154, 0.34154941, 0.34449241,
                    0.34475279, 0.35546895, 0.35551709, 0.37511313, 0.37589993,
                    0.40570240, 0.41752876
                ]),
                'osc_str': np.array([
                    0.0311000, 1.2416000, 0.0001000, 0.0139000, 0.0007000,
                    0.0001000, 0.0000000, 0.0000000, 0.0001000, 0.0016000,
                    0.0243000, 0.0053000
                ]),
                'rot_str': np.array([
                    0.1496250, -0.1570830, -8.000000e-05, 3.300000e-05,
                    0.0213310, -0.0017560, -0.0003580, -0.0028380, -0.0107960,
                    -0.0072060, 0.0012530, 0.0162620
                ])
            }
        }[ri]

        self.run_input_exmod(capsys, xcfun, ri, ref_data)
