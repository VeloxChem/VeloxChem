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

        #TODO: replace files with appropriate input data and
        # reference values, once (TD)DFT is available
        # (current files are for TDHF/TDA)  
        inpfile = str(here / 'inputs' / 'water.inp')
        lambda_ref_file = str(here / 'inputs' / 'lambda_ao_ref.txt')
        omega_ref_file = str(here / 'inputs' / 'omega_ao_ref.txt')
        potfile = None

        xcfun_label = 'b3lyp'

        ref_lambda_ao = np.loadtxt(lambda_ref_file, delimiter=' ')
        ref_omega_ao = np.loadtxt(omega_ref_file, delimiter=' ')

        self.run_orbitalresponse(inpfile, potfile, xcfun_label,
                                 ref_lambda_ao, ref_omega_ao)
 

    @pytest.mark.skipif(True, reason='DFT gradients not available')
    def test_tda_dft_slda(self):

        here = Path(__file__).parent

        #TODO: replace files with appropriate input data and
        # reference values, once (TD)DFT is available
        # (current files are for TDHF/TDA)  
        inpfile = str(here / 'inputs' / 'water.inp')
        lambda_ref_file = str(here / 'inputs' / 'lambda_ao_ref.txt')
        omega_ref_file = str(here / 'inputs' / 'omega_ao_ref.txt')
        potfile = None

        xcfun_label = 'slda'

        ref_lambda_ao = np.loadtxt(lambda_ref_file, delimiter=' ')
        ref_omega_ao = np.loadtxt(omega_ref_file, delimiter=' ')

        self.run_orbitalresponse(inpfile, potfile, xcfun_label,
                                 ref_lambda_ao, ref_omega_ao)

    #@pytest.mark.skipif('cppe' not in sys.modules, reason='cppe not available')
    @pytest.mark.skipif(True, reason='cppe not available')
    def test_tda_hf_pe(self):

        here = Path(__file__).parent

        #TODO: replace files with appropriate input data and
        # reference values, once (TD)HF with PE is available
        # (current files are for TDHF/TDA) 
        inpfile = str(here / 'inputs' / 'pe_water.inp')
        potfile = str(here / 'inputs' / 'pe_water.pot')

        lambda_ref_file = str(here / 'inputs' / 'lambda_ao_ref.txt')
        omega_ref_file = str(here / 'inputs' / 'omega_ao_ref.txt')

        xcfun_label = None

        ref_lambda_ao = np.loadtxt(lambda_ref_file, delimiter=' ')
        ref_omega_ao = np.loadtxt(omega_ref_file, delimiter=' ')

        self.run_orbitalresponse(inpfile, potfile, xcfun_label,
                                 ref_lambda_ao, ref_omega_ao)

    #@pytest.mark.skipif('cppe' not in sys.modules, reason='cppe not available')
    @pytest.mark.skipif(True, reason='cppe not available')
    def test_tda_dft_pe(self):

        here = Path(__file__).parent

        #TODO: replace files with appropriate input data and
        # reference values, once (TD)DFT with PE is available
        # (current files are for TDHF/TDA) 
        inpfile = str(here / 'inputs' / 'pe_water.inp')
        potfile = str(here / 'inputs' / 'pe_water.pot')

        lambda_ref_file = str(here / 'inputs' / 'lambda_ao_ref.txt')
        omega_ref_file = str(here / 'inputs' / 'omega_ao_ref.txt')

        xcfun_label = 'b3lyp'

        ref_lambda_ao = np.loadtxt(lambda_ref_file, delimiter=' ')
        ref_omega_ao = np.loadtxt(omega_ref_file, delimiter=' ')

        self.run_orbitalresponse(inpfile, potfile, xcfun_label,
                                 ref_lambda_ao, ref_omega_ao)


if __name__ == "__main__":
    unittest.main()
