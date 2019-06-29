from mpi4py import MPI
import numpy as np
import unittest
import os

from veloxchem.veloxchemlib import mpi_master
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.tdaexcidriver import TDAExciDriver
from veloxchem.lreigensolver import LinearResponseEigenSolver
from veloxchem.lrsolver import LinearResponseSolver
from veloxchem.mpitask import MpiTask


class TestRspDriver(unittest.TestCase):

    def test_h2se_tda(self):

        # scf
        inpfile = os.path.join('inputs', 'h2se.inp')
        if not os.path.isfile(inpfile):
            inpfile = os.path.join('python_tests', inpfile)

        task = MpiTask([inpfile, None], MPI.COMM_WORLD)
        scf_drv = ScfRestrictedDriver(task.mpi_comm, task.ostream)

        scf_drv.compute(task.molecule, task.ao_basis, task.min_basis)
        scf_tensors = scf_drv.scf_tensors

        # TDA
        tda_exci = TDAExciDriver(task.mpi_comm, task.ostream)

        tda_exci.update_settings({
            'nstates': 3,
            'eri_thresh': 1.0e-12,
            'conv_thresh': 1.0e-4,
        })

        tda_results = tda_exci.compute(task.molecule, task.ao_basis,
                                       scf_tensors)

        if task.mpi_rank == mpi_master():

            reigs = tda_results['eigenvalues']
            osc_strs = tda_results['oscillator_strengths']
            trans_dipoles = tda_results['electric_transition_dipoles']

            ref_eigs = np.array([0.207436, 0.257474, 0.368358])
            ref_osc_strs = np.array([0.0000, 0.0003, 0.2797])
            ref_trans_dipoles = [
                np.array([0.0, 0.0, 0.0]),
                np.array([-0.043213, 0.0, 0.0]),
                np.array([0.0, -0.754589, 0.754589]),
            ]

            self.assertTrue(np.max(np.abs(reigs - ref_eigs)) < 1.0e-6)
            self.assertTrue(np.max(np.abs(osc_strs - ref_osc_strs)) < 1.0e-4)

            for td, ref_td in zip(trans_dipoles, ref_trans_dipoles):
                prefac = 1.0 if np.dot(td, ref_td) >= 0.0 else -1.0
                self.assertTrue(np.max(np.abs(td - ref_td * prefac)) < 1.0e-4)

        # RPA
        lreig = LinearResponseEigenSolver(task.mpi_comm, task.ostream)

        lreig.update_settings({
            'nstates': 3,
            'eri_thresh': 1.0e-12,
            'conv_thresh': 1.0e-4,
        })

        rpa_results = lreig.compute(task.molecule, task.ao_basis, scf_tensors)

        if task.mpi_rank == mpi_master():

            reigs = rpa_results['eigenvalues']
            osc_strs = rpa_results['oscillator_strengths']
            trans_dipoles = rpa_results['electric_transition_dipoles']

            ref_eigs = np.array([0.20565979, 0.25474355, 0.36246841])
            ref_osc_strs = np.array([0.0000, 0.0012, 0.2477])
            ref_trans_dipoles = [
                np.array([0.0, 0.0, 0.0]),
                np.array([0.084626, 0.0, 0.0]),
                np.array([0.0, -0.715935, 0.715935])
            ]

            self.assertTrue(np.max(np.abs(reigs - ref_eigs)) < 1.0e-6)
            self.assertTrue(np.max(np.abs(osc_strs - ref_osc_strs)) < 1.0e-4)

            for td, ref_td in zip(trans_dipoles, ref_trans_dipoles):
                prefac = 1.0 if np.dot(td, ref_td) >= 0.0 else -1.0
                self.assertTrue(np.max(np.abs(td - ref_td * prefac)) < 1.0e-4)

        # polarizability
        lr_solver = LinearResponseSolver(task.mpi_comm, task.ostream)

        lr_solver.update_settings({
            'frequencies': '0, 0.1',
            'eri_thresh': '1.0e-12',
            'conv_thresh': '1.0e-4',
        })

        lr_prop = lr_solver.compute(task.molecule, task.ao_basis, scf_tensors)

        if task.mpi_rank == mpi_master():

            ref_polar = {
                ('x', 'x', 0): 15.26732,
                ('x', 'y', 0): 0.0,
                ('x', 'z', 0): 0.0,
                ('y', 'x', 0): 0.0,
                ('y', 'y', 0): 23.22477,
                ('y', 'z', 0): -0.759275,
                ('z', 'x', 0): 0.0,
                ('z', 'y', 0): -0.759275,
                ('z', 'z', 0): 23.22477,
                ('x', 'x', 0.1): 15.63176,
                ('x', 'y', 0.1): 0.0,
                ('x', 'z', 0.1): 0.0,
                ('y', 'x', 0.1): 0.0,
                ('y', 'y', 0.1): 24.27442,
                ('y', 'z', 0.1): -0.85701,
                ('z', 'x', 0.1): 0.0,
                ('z', 'y', 0.1): -0.85701,
                ('z', 'z', 0.1): 24.27442,
            }

            for key, val in lr_prop.items():
                ref = -ref_polar[key]
                diff = abs(val - ref)
                self.assertTrue(diff < 1.0e-5)

        task.finish()


if __name__ == "__main__":
    unittest.main()
