from mpi4py import MPI
import numpy as np
import unittest
import pytest
import sys
import os
from pathlib import Path
try:
    import cppe
except ImportError:
    pass

from veloxchem.veloxchemlib import mpi_master
from veloxchem.veloxchemlib import hartree_in_ev
from veloxchem.mpitask import MpiTask
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.tdaexcidriver import TDAExciDriver
from veloxchem.tdaorbitalresponse import TdaOrbitalResponse


class TestOrbitalResponse(unittest.TestCase):

    # TODO: this will have to be more general once RPA gradients are implemented
    def run_orbitalresponse(self, inpfile, potfile, xcfun_label,
                            ref_lambda_ao, ref_omega_ao):

        task = MpiTask([inpfile, None], MPI.COMM_WORLD)
        task.input_dict['scf']['checkpoint_file'] = None
		# use a smaller basis for lambda test

        if potfile is not None:
            task.input_dict['method_settings']['potfile'] = potfile

        if xcfun_label is not None:
            task.input_dict['method_settings']['xcfun'] = xcfun_label

        scf_drv = ScfRestrictedDriver(task.mpi_comm, task.ostream)
        scf_drv.update_settings(task.input_dict['scf'],
                                task.input_dict['method_settings'])
        scf_drv.compute(task.molecule, task.ao_basis, task.min_basis)

		# Our references: lambda and omega in AO basis

        tda_solver = TDAExciDriver(task.mpi_comm, task.ostream)
        tda_solver.update_settings({'nstates': 3},
                                   task.input_dict['method_settings'])
        tda_results = tda_solver.compute(task.molecule, task.ao_basis,
                                         scf_drv.scf_tensors)
        tda_eig_vecs = tda_results['eigenvectors']

        orb_resp = TdaOrbitalResponse(task.mpi_comm, task.ostream)
        orb_resp.update_settings({'nstates': 3,
                                 'n_state_deriv': 1},
                                 task.input_dict['method_settings'])
        orb_resp.compute(task.molecule, task.ao_basis,
                         scf_drv.scf_tensors, tda_eig_vecs)

        if task.mpi_rank == mpi_master():
            lambda_ao = orb_resp.lambda_ao
            omega_ao = orb_resp.omega_ao

            self.assertTrue(np.max(np.abs(lambda_ao - ref_lambda_ao)) < 5.0e-4)
            self.assertTrue(np.max(np.abs(omega_ao - ref_omega_ao)) < 5.0e-4)

    def test_tda_hf(self):

        here = Path(__file__).parent
        inpfile = str(here / 'inputs' / 'water_orbrsp.inp')
        lambda_ref_file = str(here / 'inputs' / 'lambda_ao_ref.txt')
        omega_ref_file = str(here / 'inputs' / 'omega_ao_ref.txt')

        potfile = None

        xcfun_label = None

        ref_lambda_ao = np.loadtxt(lambda_ref_file, delimiter=' ')
        ref_omega_ao = np.loadtxt(omega_ref_file, delimiter=' ')

        self.run_orbitalresponse(inpfile, potfile, xcfun_label,
                                 ref_lambda_ao, ref_omega_ao)

    @pytest.mark.skipif(True, reason='DFT gradients not available')
    def test_tda_dft(self):

        here = Path(__file__).parent
        inpfile = str(here / 'inputs' / 'water.inp')

        potfile = None

        xcfun_label = 'b3lyp'

        #  State Frequency   Oscillator Strength    Rotatory  Strength
        #          (eV)      Velocity     Length    Velocity    Length
        #  -----------------------------------------------------------
        raw_data = """
            1     6.9987     0.0571     0.0537     0.0000     0.0000
            2     8.4291     0.0000     0.0000    -0.0000    -0.0000
            3     9.2135     0.0611     0.0906     0.0000     0.0000
            4    10.3003     0.0005     0.0000    -0.0000    -0.0000
            5    10.6270     0.0044     0.0127    -0.0000    -0.0000
        """
        data_lines = raw_data.split(os.linesep)[1:-1]

        self.run_orbitalresponse(inpfile, potfile, xcfun_label, data_lines)

    @pytest.mark.skipif(True, reason='DFT gradients not available')
    def test_tda_dft_slda(self):

        here = Path(__file__).parent
        inpfile = str(here / 'inputs' / 'water.inp')

        potfile = None

        xcfun_label = 'slda'

        #  State Frequency   Oscillator Strength    Rotatory  Strength
        #          (eV)      Velocity     Length    Velocity    Length
        #  -----------------------------------------------------------
        raw_data = """
            1     6.7828     0.0588     0.0561    -0.0000     0.0000
            2     8.2221     0.0000     0.0000     0.0000     0.0000
            3     8.9101     0.0603     0.0901    -0.0000    -0.0000
            4    10.1323     0.0014     0.0003     0.0000     0.0000
            5    10.3444     0.0036     0.0115    -0.0000    -0.0000
        """
        data_lines = raw_data.split(os.linesep)[1:-1]

        self.run_orbitalresponse(inpfile, potfile, xcfun_label, data_lines)

    @pytest.mark.skipif('cppe' not in sys.modules, reason='cppe not available')
    def test_tda_hf_pe(self):

        here = Path(__file__).parent
        inpfile = str(here / 'inputs' / 'pe_water.inp')
        potfile = str(here / 'inputs' / 'pe_water.pot')

        xcfun_label = None

        #  State Frequency   Oscillator Strength    Rotatory  Strength
        #          (eV)      Velocity     Length    Velocity    Length
        #  -----------------------------------------------------------
        raw_data = """
            1     9.1467     0.0758     0.0602     0.1440     0.3887
            2    11.1635     0.0113     0.0121    -0.4940    -0.2295
            3    11.4450     0.1062     0.1234    -0.1259     2.8161
            4    11.9792     0.0004     0.0013     0.0916     0.1963
            5    12.8941     0.0007     0.0006     0.0132    -0.4366
        """
        data_lines = raw_data.split(os.linesep)[1:-1]

        self.run_orbitalresponse(inpfile, potfile, xcfun_label, data_lines)

    @pytest.mark.skipif('cppe' not in sys.modules, reason='cppe not available')
    def test_tda_dft_pe(self):

        here = Path(__file__).parent
        inpfile = str(here / 'inputs' / 'pe_water.inp')
        potfile = str(here / 'inputs' / 'pe_water.pot')

        xcfun_label = 'b3lyp'

        #  State Frequency   Oscillator Strength    Rotatory  Strength
        #          (eV)      Velocity     Length    Velocity    Length
        #  -----------------------------------------------------------
        raw_data = """
            1     7.3245     0.0609     0.0582     0.1219     0.3536
            2     9.0048     0.0067     0.0073    -0.0937     0.1897
            3     9.5204     0.0742     0.1000    -0.2837     2.3062
            4    10.1845     0.0003     0.0003    -0.1075    -0.2451
            5    11.0440     0.0076     0.0029     0.1461    -0.5130
        """
        data_lines = raw_data.split(os.linesep)[1:-1]

        self.run_orbitalresponse(inpfile, potfile, xcfun_label, data_lines)


if __name__ == "__main__":
    unittest.main()
