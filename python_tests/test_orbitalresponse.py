from mpi4py import MPI
import numpy as np
import unittest
from pathlib import Path
try:
    import cppe
except ImportError:
    pass

from veloxchem.veloxchemlib import mpi_master
from veloxchem.mpitask import MpiTask
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.tdaexcidriver import TDAExciDriver
from veloxchem.tdaorbitalresponse import TdaOrbitalResponse
from veloxchem.checkpoint import read_rsp_hdf5


class TestOrbitalResponse(unittest.TestCase):

    def run_tdaorbitalresponse(self, inpfile, potfile, xcfun_label,
                               orbrsp_ref_file):

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

        orb_resp = TdaOrbitalResponse(task.mpi_comm, task.ostream)
        orb_resp.update_settings({
            'nstates': 3,
            'n_state_deriv': 1
        }, task.input_dict['method_settings'])
        orb_resp_result = orb_resp.compute(task.molecule, task.ao_basis,
                                           scf_drv.scf_tensors, tda_results)

        dft_dict = {'dft_func_label': 'HF'}
        pe_dict = {'potfile_text': ''}

        ref_lambda_ao, ref_omega_ao = read_rsp_hdf5(orbrsp_ref_file,
                                                    ['lambda_tda', 'omega_tda'],
                                                    task.molecule,
                                                    task.ao_basis, dft_dict,
                                                    pe_dict, task.ostream)

        if task.mpi_rank == mpi_master():
            lambda_ao = orb_resp_result['lambda_ao']
            omega_ao = orb_resp_result['omega_ao']

            self.assertTrue(np.max(np.abs(lambda_ao - ref_lambda_ao)) < 5.0e-4)
            self.assertTrue(np.max(np.abs(omega_ao - ref_omega_ao)) < 5.0e-4)

    def test_tda_hf(self):

        here = Path(__file__).parent
        inpfile = str(here / 'inputs' / 'water_orbrsp.inp')
        orbrsp_ref_file = str(here / 'inputs' / 'orbital_response_hf_ref.h5')

        potfile = None

        xcfun_label = None

        self.run_tdaorbitalresponse(inpfile, potfile, xcfun_label,
                                    orbrsp_ref_file)


if __name__ == "__main__":
    unittest.main()
