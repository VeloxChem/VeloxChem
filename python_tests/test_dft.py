from mpi4py import MPI
import numpy as np
import unittest
import os

from veloxchem.mpitask import MpiTask
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.tdaexcidriver import TDAExciDriver
from veloxchem.lreigensolver import LinearResponseEigenSolver
from veloxchem.veloxchemlib import mpi_master
from veloxchem.veloxchemlib import hartree_in_ev


class TestDFT(unittest.TestCase):

    def run_dft(self, xcfun_label, scf_ref, tda_ref, rpa_ref):

        inpfile = os.path.join('inputs', 'water.inp')
        if not os.path.isfile(inpfile):
            inpfile = os.path.join('python_tests', inpfile)

        task = MpiTask([inpfile, None], MPI.COMM_WORLD)

        # SCF

        scf_drv = ScfRestrictedDriver(task.mpi_comm, task.ostream)
        scf_drv.update_settings({}, {'xcfun': xcfun_label})
        scf_drv.compute(task.molecule, task.ao_basis, task.min_basis)

        if task.mpi_rank == mpi_master():
            e_scf = scf_drv.get_scf_energy()
            self.assertTrue(abs(e_scf - scf_ref) < 1.0e-5)

        # TDA

        tda_exci = TDAExciDriver(task.mpi_comm, task.ostream)
        tda_exci.update_settings({'nstates': 3}, {'xcfun': xcfun_label})
        tda_results = tda_exci.compute(task.molecule, task.ao_basis,
                                       scf_drv.scf_tensors)

        if task.mpi_rank == mpi_master():
            exc_ene = tda_results['eigenvalues'] * hartree_in_ev()
            osc_str = tda_results['oscillator_strengths']

            self.assertTrue(np.max(np.abs(exc_ene - tda_ref[:, 0])) < 5.0e-4)
            self.assertTrue(np.max(np.abs(osc_str - tda_ref[:, 1])) < 1.0e-4)

        # RPA

        rpa_exci = LinearResponseEigenSolver(task.mpi_comm, task.ostream)
        rpa_exci.update_settings({'nstates': 3}, {'xcfun': xcfun_label})
        rpa_results = rpa_exci.compute(task.molecule, task.ao_basis,
                                       scf_drv.scf_tensors)

        if task.mpi_rank == mpi_master():
            exc_ene = rpa_results['eigenvalues'] * hartree_in_ev()
            osc_str = rpa_results['oscillator_strengths']

            self.assertTrue(np.max(np.abs(exc_ene - rpa_ref[:, 0])) < 5.0e-4)
            self.assertTrue(np.max(np.abs(osc_str - rpa_ref[:, 1])) < 1.0e-4)

    def test_bhandhlyp(self):

        xcfun_label = 'bhandhlyp'

        scf_ref = -76.40127686

        tda_ref = np.array([
            [7.765647, 0.0536],
            [9.349421, 0.0000],
            [10.050452, 0.0976],
        ])

        rpa_ref = np.array([
            [7.747461, 0.0519],
            [9.341142, 0.0000],
            [10.032598, 0.0916],
        ])

        self.run_dft(xcfun_label, scf_ref, tda_ref, rpa_ref)


if __name__ == "__main__":
    unittest.main()
