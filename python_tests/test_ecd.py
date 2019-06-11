from mpi4py import MPI
import numpy as np
import unittest
import os

from veloxchem.veloxchemlib import mpi_master
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.lreigensolver import LinearResponseEigenSolver
from veloxchem.mpitask import MpiTask


class TestECD(unittest.TestCase):

    def test_glycol_ecd(self):

        # scf
        inpfile = os.path.join('inputs', 'glycol.inp')
        if not os.path.isfile(inpfile):
            inpfile = os.path.join('python_tests', inpfile)

        task = MpiTask([inpfile, None], MPI.COMM_WORLD)
        scf_drv = ScfRestrictedDriver(task.mpi_comm, task.ostream)

        scf_drv.compute(task.molecule, task.ao_basis, task.min_basis)
        scf_tensors = scf_drv.scf_tensors

        # ecd
        lreig_solver = LinearResponseEigenSolver(task.mpi_comm, task.ostream)

        lreig_solver.update_settings({
            'nstates': '10',
            'eri_thresh': '1.0e-12',
            'qq_type': 'QQ_DEN',
            'conv_thresh': '1.0e-6',
            'max_iter': '50',
        })

        results = lreig_solver.compute(task.molecule, task.ao_basis,
                                       scf_tensors)

        if task.mpi_rank == mpi_master():

            reigs = results['eigenvalues']
            osc_strs = results['oscillator_strengths']
            rot_strs = results['rotatory_strengths']

            ref_eigs = np.array([
                0.447941, 0.454425, 0.510293, 0.519992, 0.567876, 0.569685,
                0.633078, 0.636154, 0.643907, 0.659273
            ])
            ref_osc_strs = np.array([
                0.0055, 0.0064, 0.0037, 0.0008, 0.0224, 0.0159, 0.0869, 0.1293,
                0.0629, 0.1133
            ])
            ref_rot_strs = np.array([
                -0.0169, 0.0191, -0.0183, 0.0181, -0.2282, 0.1686, 0.1913,
                -0.1511, -0.0023, -0.0157
            ])

            self.assertTrue(np.max(np.abs(reigs - ref_eigs)) < 1.0e-6)
            self.assertTrue(np.max(np.abs(osc_strs - ref_osc_strs)) < 1.0e-3)
            self.assertTrue(np.max(np.abs(rot_strs - ref_rot_strs)) < 1.0e-3)

        task.finish()


if __name__ == "__main__":
    unittest.main()
