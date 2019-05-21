from mpi4py import MPI
import numpy as np
import unittest
import os

from veloxchem.veloxchemlib import mpi_master
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.tdaexcidriver import TDAExciDriver
from veloxchem.lrsolver import LinearResponseSolver
from veloxchem.mpitask import MpiTask


class TestRspDriver(unittest.TestCase):

    def test_h2se_tda(self):

        # scf
        inpfile = os.path.join('inputs', 'h2se.inp')
        task = MpiTask([inpfile, None], MPI.COMM_WORLD)
        scf_drv = ScfRestrictedDriver(task.mpi_comm, task.ostream)

        scf_drv.compute(task.molecule, task.ao_basis, task.min_basis)
        scf_tensors = scf_drv.scf_tensors

        # TDA
        tda_exci = TDAExciDriver(task.mpi_comm, task.ostream)

        tda_exci.update_settings({
            'nstates': 3,
            'spin': 'singlet',
            'eri_thresh': 1.0e-12,
            'qq_type': 'QQ_DEN',
            'conv_thresh': 1.0e-4,
            'max_iter': 50,
        })

        tda_results = tda_exci.compute(task.molecule, task.ao_basis,
                                       scf_tensors)

        if task.mpi_rank == mpi_master():

            reigs = tda_results['eigenvalues']
            osc_strs = tda_results['oscillator_strengths']
            trans_dipoles = tda_results['transition_dipoles']

            ref_eigs = np.array([0.207436, 0.257474, 0.368358])
            ref_osc_strs = np.array([0.0000, 0.0003, 0.2797])
            ref_trans_dipoles = [
                np.array([0.0000, 0.0000, 0.0000]),
                np.array([0.0432, 0.0000, 0.0000]),
                np.array([0.0000, -0.7546, 0.7546]),
            ]

            self.assertTrue(np.max(np.abs(reigs - ref_eigs)) < 1.0e-6)
            self.assertTrue(np.max(np.abs(osc_strs - ref_osc_strs)) < 1.0e-4)

            for td, ref_td in zip(trans_dipoles, ref_trans_dipoles):
                prefac = 1.0 if np.dot(td, ref_td) >= 0.0 else -1.0
                self.assertTrue(np.max(np.abs(td - ref_td * prefac)) < 1.0e-4)

        # polarizability
        lr_solver = LinearResponseSolver(task.mpi_comm, task.ostream)

        lr_solver.update_settings({
            'a_ops': 'xyz',
            'b_ops': 'xyz',
            'frequencies': (0,),
            'eri_thresh': 1.0e-12,
            'qq_type': 'QQ_DEN',
            'conv_thresh': 1.0e-5,
            'max_iter': 50,
        })

        lr_prop = lr_solver.compute(task.molecule, task.ao_basis, scf_tensors)

        if task.mpi_rank == mpi_master():

            ref_polar = {
                ('x', 'x', 0): 15.2673,
                ('x', 'y', 0): 0.0,
                ('y', 'x', 0): 0.0,
                ('y', 'y', 0): 23.2248,
                ('x', 'z', 0): 0.0,
                ('z', 'x', 0): 0.0,
                ('y', 'z', 0): -0.759273,
                ('z', 'y', 0): -0.759273,
                ('z', 'z', 0): 23.2248,
            }

            for key, val in lr_prop.items():
                ref = -ref_polar[key]
                diff = abs(val - ref)
                if abs(ref) > 0.0:
                    rel_diff = abs(diff / ref)
                    self.assertTrue(rel_diff < 5.0e-6)
                else:
                    self.assertTrue(diff < 1.0e-10)

        task.finish()


if __name__ == "__main__":
    unittest.main()
