from mpi4py import MPI
import numpy as np
import unittest
import os

from veloxchem.mpitask import MpiTask
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.tdaexcidriver import TDAExciDriver
from veloxchem.veloxchemlib import hartree_in_ev


class TestDFT(unittest.TestCase):

    def run_tda(self, xcfun_label, ref_exc_ene, ref_osc_str):

        inpfile = os.path.join('inputs', 'water.inp')
        if not os.path.isfile(inpfile):
            inpfile = os.path.join('python_tests', inpfile)

        task = MpiTask([inpfile, None], MPI.COMM_WORLD)

        scf_drv = ScfRestrictedDriver(task.mpi_comm, task.ostream)
        scf_drv.update_settings({}, {'xcfun': xcfun_label})
        scf_drv.compute(task.molecule, task.ao_basis, task.min_basis)

        tda_exci = TDAExciDriver(task.mpi_comm, task.ostream)
        tda_exci.update_settings({'nstates': 5}, {'xcfun': xcfun_label})
        tda_results = tda_exci.compute(task.molecule, task.ao_basis,
                                       scf_drv.scf_tensors)

        exc_ene = tda_results['eigenvalues'] * hartree_in_ev()
        osc_str = tda_results['oscillator_strengths']

        self.assertTrue(np.max(np.abs(exc_ene - ref_exc_ene)) < 5.0e-4)
        self.assertTrue(np.max(np.abs(osc_str - ref_osc_str)) < 1.0e-4)

    def test_bhandhlyp(self):

        data = np.array([
            [7.765647, 0.0536],
            [9.349421, 0.0000],
            [10.050452, 0.0976],
            [11.118765, 0.0013],
            [11.611597, 0.0162],
        ])

        xcfun_label = 'bhandhlyp'
        ref_exc_ene = data[:, 0]
        ref_osc_str = data[:, 1]

        self.run_tda(xcfun_label, ref_exc_ene, ref_osc_str)


if __name__ == "__main__":
    unittest.main()
