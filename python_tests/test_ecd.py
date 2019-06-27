from mpi4py import MPI
import numpy as np
import unittest
import os

from veloxchem.veloxchemlib import mpi_master
from veloxchem.veloxchemlib import rotatory_strength_in_cgs
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
            'conv_thresh': '1.0e-4',
        })

        results = lreig_solver.compute(task.molecule, task.ao_basis,
                                       scf_tensors)

        if task.mpi_rank == mpi_master():

            reigs = results['eigenvalues']
            osc_strs = results['oscillator_strengths']
            rot_strs = results['rotatory_strengths']
            rot_strs /= rotatory_strength_in_cgs()

            elec_tms = np.array(results['electric_transition_dipoles'])
            velo_tms = np.array(results['velocity_transition_dipoles'])
            magn_tms = np.array(results['magnetic_transition_dipoles'])

            ref_eigs = np.array([
                0.447941, 0.454425, 0.510293, 0.519992, 0.567876, 0.569685,
                0.633077, 0.636153, 0.643907, 0.659273
            ])
            ref_osc_strs = np.array([
                0.0055, 0.0064, 0.0037, 0.0008, 0.0224, 0.0158, 0.0869, 0.1295,
                0.0629, 0.1133
            ])
            ref_rot_strs = np.array([
                -0.0084, 0.0095, -0.0091, 0.0091, -0.1142, 0.0842, 0.0955,
                -0.0758, -0.0012, -0.0079
            ])

            ref_elec_tms = np.array([[-0.002583, 0.135060, -0.017078],
                                     [0.133762, -0.032233, -0.047421],
                                     [0.002016, -0.059787, -0.085703],
                                     [-0.000505, -0.034204, 0.031874],
                                     [-0.120260, -0.116025, 0.176943],
                                     [-0.147290, -0.018351, -0.140257],
                                     [0.148391, -0.358063, 0.235909],
                                     [0.000985, -0.310990, -0.456661],
                                     [-0.176175, -0.329862, -0.081070],
                                     [0.375434, -0.327897, -0.096529]])

            ref_velo_tms = np.array([[0.011785, 0.034739, -0.004039],
                                     [0.170785, -0.012581, -0.015424],
                                     [-0.007698, 0.006435, 0.032807],
                                     [-0.053912, 0.028785, -0.002683],
                                     [0.123161, 0.234646, -0.161814],
                                     [0.134976, -0.007327, 0.226427],
                                     [-0.090727, 0.151253, -0.117001],
                                     [-0.006175, 0.186112, 0.191462],
                                     [0.101400, 0.170126, 0.057920],
                                     [-0.183785, 0.149051, 0.039218]])
            ref_velo_tms[:, 0] /= -ref_eigs
            ref_velo_tms[:, 1] /= -ref_eigs
            ref_velo_tms[:, 2] /= -ref_eigs

            ref_magn_tms = np.array([[0.350275, 0.154203, 0.482833],
                                     [-0.065412, -0.600916, 0.327435],
                                     [0.212374, -0.042437, 0.342707],
                                     [-0.048419, -0.389297, 0.313643],
                                     [1.430297, 0.115747, 0.455256],
                                     [-0.921004, -1.647440, 0.072091],
                                     [-0.110409, -0.131985, 0.948317],
                                     [-0.652499, -0.419945, 0.890720],
                                     [-0.340337, -0.064951, 0.812425],
                                     [-0.162535, -0.305578, 0.664005]])
            ref_magn_tms *= 0.5

            self.assertTrue(np.max(np.abs(reigs - ref_eigs)) < 1.0e-6)
            self.assertTrue(np.max(np.abs(osc_strs - ref_osc_strs)) < 1.0e-4)
            self.assertTrue(np.max(np.abs(rot_strs - ref_rot_strs)) < 1.0e-4)

            for i in range(elec_tms.shape[0]):
                if np.vdot(elec_tms[i, :], ref_elec_tms[i, :]) < 0.0:
                    elec_tms[i, :] *= -1.0
                if np.vdot(velo_tms[i, :], ref_velo_tms[i, :]) < 0.0:
                    velo_tms[i, :] *= -1.0
                if np.vdot(magn_tms[i, :], ref_magn_tms[i, :]) < 0.0:
                    magn_tms[i, :] *= -1.0
            self.assertTrue(np.max(np.abs(elec_tms - ref_elec_tms)) < 1.0e-4)
            self.assertTrue(np.max(np.abs(velo_tms - ref_velo_tms)) < 2.0e-4)
            self.assertTrue(np.max(np.abs(magn_tms - ref_magn_tms)) < 2.0e-4)

        task.finish()


if __name__ == "__main__":
    unittest.main()
