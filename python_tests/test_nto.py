from mpi4py import MPI
import numpy as np
import unittest
from pathlib import Path

from veloxchem.veloxchemlib import mpi_master
from veloxchem.mpitask import MpiTask
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.tdaexcidriver import TDAExciDriver


class TestNTO(unittest.TestCase):

    def run_nto(self, inpfile, ref_eig_vals, ref_nto_vals):

        task = MpiTask([inpfile, None], MPI.COMM_WORLD)
        task.input_dict['scf']['checkpoint_file'] = None

        scf_drv = ScfRestrictedDriver(task.mpi_comm, task.ostream)
        scf_drv.update_settings(task.input_dict['scf'],
                                task.input_dict['method_settings'])
        scf_drv.compute(task.molecule, task.ao_basis, task.min_basis)

        rsp_drv = TDAExciDriver(task.mpi_comm, task.ostream)
        rsp_drv.update_settings({'nstates': ref_eig_vals.shape[0]},
                                task.input_dict['method_settings'])
        rsp_results = rsp_drv.compute(task.molecule, task.ao_basis,
                                      scf_drv.scf_tensors)

        if task.mpi_rank == mpi_master():
            fock = scf_drv.scf_tensors['F'][0]

            mo = scf_drv.scf_tensors['C']
            nocc = task.molecule.number_of_alpha_electrons()
            mo_occ = mo[:, :nocc]
            mo_vir = mo[:, nocc:]

            eig_vals = rsp_results['eigenvalues']
            eig_vecs = rsp_results['eigenvectors']

            nto_vals = []

            for s in range(ref_eig_vals.shape[0]):
                t_mat = eig_vecs[:, s].reshape(mo_occ.shape[1], mo_vir.shape[1])
                u_mat, s_diag, vh_mat = np.linalg.svd(t_mat, full_matrices=True)
                lam_diag = s_diag**2
                nto_occ = np.matmul(mo_occ, u_mat)
                nto_vir = np.matmul(mo_vir, vh_mat.T)

                for i_nto in range(lam_diag.size):
                    if lam_diag[i_nto] < 0.1:
                        continue
                    occ_vec = nto_occ[:, i_nto]
                    e_hole = np.vdot(occ_vec, np.matmul(fock, occ_vec))
                    vir_vec = nto_vir[:, i_nto]
                    e_particle = np.vdot(vir_vec, np.matmul(fock, vir_vec))
                    nto_vals += [lam_diag[i_nto], e_hole, e_particle]

            nto_vals = np.array(nto_vals)

            self.assertTrue(np.max(np.abs(eig_vals - ref_eig_vals)) < 1.0e-6)
            self.assertTrue(np.max(np.abs(nto_vals - ref_nto_vals)) < 1.0e-4)

    def test_nto_tda(self):

        here = Path(__file__).parent
        inpfile = str(here / 'inputs' / 'c2h4.inp')

        raw_vals = """
            0.316434
            0.346741
            0.348663
            0.358313
            0.386238
        """
        eig_vals = np.array([float(x) for x in raw_vals.split()])

        raw_vals = """
            0.9400 -0.3808 0.1698
            0.6581 -0.5038 0.1758
            0.3401 -0.3808 0.2487
            0.9986 -0.3808 0.2132
            0.6602 -0.3808 0.2456
            0.3390 -0.5038 0.1750
            0.6113 -0.5995 0.1846
            0.3873 -0.3808 0.2497
        """
        nto_vals = np.array([float(x) for x in raw_vals.split()])

        self.run_nto(inpfile, eig_vals, nto_vals)


if __name__ == "__main__":
    unittest.main()
