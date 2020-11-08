from mpi4py import MPI
import numpy as np
import unittest
from pathlib import Path

from veloxchem.veloxchemlib import mpi_master
from veloxchem.mpitask import MpiTask
from veloxchem.cubicgrid import CubicGrid
from veloxchem.molecularorbitals import MolecularOrbitals
from veloxchem.visualizationdriver import VisualizationDriver
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.tdaexcidriver import TDAExciDriver


class TestNTO(unittest.TestCase):

    def run_nto(self, inpfile, ref_eig_vals, ref_nto_vals, ref_cube_vals):

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

        nocc = task.molecule.number_of_alpha_electrons()

        if task.mpi_rank == mpi_master():
            fock = scf_drv.scf_tensors['F'][0]

            mo = scf_drv.scf_tensors['C']
            mo_occ = mo[:, :nocc]
            mo_vir = mo[:, nocc:]

            eig_vals = rsp_results['eigenvalues']
            eig_vecs = rsp_results['eigenvectors']

            nto_vals = []
            cube_vals = []

        for s in range(ref_eig_vals.shape[0]):
            if task.mpi_rank == mpi_master():
                lam_diag, nto_mo = rsp_drv.get_nto(s, eig_vecs, mo_occ, mo_vir)
            else:
                lam_diag = None
                nto_mo = MolecularOrbitals()
            lam_diag = task.mpi_comm.bcast(lam_diag, root=mpi_master())
            nto_mo.broadcast(task.mpi_rank, task.mpi_comm)

            grid = CubicGrid([-0.1, -0.1, -0.1], [0.3, 0.3, 0.2], [2, 2, 3])

            vis_drv = VisualizationDriver(task.mpi_comm)

            i_nto = 0

            ind_occ = nocc - i_nto - 1
            vis_drv.compute(grid, task.molecule, task.ao_basis, nto_mo, ind_occ,
                            'alpha')

            if task.mpi_rank == mpi_master():
                occ_vec = nto_mo.alpha_to_numpy()[:, ind_occ]
                e_hole = np.vdot(occ_vec, np.matmul(fock, occ_vec))

                cube_vals.append(grid.values_to_numpy().flatten())

            ind_vir = nocc + i_nto
            vis_drv.compute(grid, task.molecule, task.ao_basis, nto_mo, ind_vir,
                            'alpha')

            if task.mpi_rank == mpi_master():
                vir_vec = nto_mo.alpha_to_numpy()[:, ind_vir]
                e_particle = np.vdot(vir_vec, np.matmul(fock, vir_vec))

                cube_vals.append(grid.values_to_numpy().flatten())

                nto_vals.append([lam_diag[i_nto], e_hole, e_particle])

        if task.mpi_rank == mpi_master():
            nto_vals = np.array(nto_vals)
            cube_vals = np.array(cube_vals)

            self.assertTrue(np.max(np.abs(eig_vals - ref_eig_vals)) < 1.0e-6)
            self.assertTrue(np.max(np.abs(nto_vals - ref_nto_vals)) < 1.0e-4)

            for s in range(cube_vals.shape[0]):
                if np.vdot(cube_vals[s, :], ref_cube_vals[s, :]) < 0.0:
                    cube_vals[s, :] *= -1.0
            self.assertTrue(np.max(np.abs(cube_vals - ref_cube_vals)) < 1.0e-4)

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
        nstates = eig_vals.size

        raw_vals = """
            0.9400 -0.3808 0.1698
            0.6581 -0.5038 0.1758
            0.9986 -0.3808 0.2132
            0.6602 -0.3808 0.2456
            0.6113 -0.5995 0.1846
        """
        nto_vals = np.array([float(x) for x in raw_vals.split()])
        nto_vals = nto_vals.reshape(nstates, -1)

        raw_vals = """
            3.08300E-02  3.08300E-02  3.41697E-02
            3.02083E-02  3.02083E-02  3.33643E-02
           -6.04152E-02 -6.04152E-02 -6.67266E-02
           -5.92120E-02 -5.92120E-02 -6.51781E-02
            4.19373E-03 -4.19373E-03 -1.33649E-02
            4.06161E-03 -4.06161E-03 -1.29179E-02
           -8.12234E-03  8.12234E-03  2.58328E-02
           -7.87085E-03  7.87085E-03  2.49837E-02
            3.11244E-03 -3.11244E-03 -1.00451E-02
           -6.02484E-03  6.02484E-03  1.94061E-02
            3.00892E-03 -3.00892E-03 -9.69226E-03
           -5.82780E-03  5.82780E-03  1.87359E-02
            4.54510E-03 -4.54510E-03 -1.45426E-02
            4.39760E-03 -4.39760E-03 -1.40423E-02
           -8.79436E-03  8.79436E-03  2.80818E-02
           -8.51372E-03  8.51372E-03  2.71320E-02
           -3.08300E-02 -3.08300E-02 -3.41697E-02
           -3.02083E-02 -3.02083E-02 -3.33643E-02
            6.04152E-02  6.04152E-02  6.67266E-02
            5.92120E-02  5.92120E-02  6.51781E-02
           -7.71797E-03 -7.71797E-03 -1.20185E-02
           -7.34438E-03 -7.34438E-03 -1.10426E-02
           -6.58723E-03 -6.58723E-03 -1.02655E-02
           -6.27868E-03 -6.27868E-03 -9.40533E-03
           -3.08300E-02 -3.08300E-02 -3.41697E-02
           -3.02083E-02 -3.02083E-02 -3.33643E-02
            6.04152E-02  6.04152E-02  6.67266E-02
            5.92120E-02  5.92120E-02  6.51781E-02
            2.62157E-02  2.62157E-02  2.84788E-02
           -5.13050E-02 -5.13050E-02 -5.55604E-02
            2.56648E-02  2.56648E-02  2.77928E-02
           -5.02375E-02 -5.02375E-02 -5.42392E-02
           -3.39884E-01 -3.39884E-01 -3.45196E-01
           -3.31899E-01 -3.31899E-01 -3.36173E-01
           -3.32479E-01 -3.32479E-01 -3.36811E-01
           -3.24750E-01 -3.24750E-01 -3.28112E-01
            4.99796E-03 -4.99796E-03 -1.59590E-02
            4.83409E-03 -4.83409E-03 -1.54049E-02
           -9.66732E-03  9.66732E-03  3.08069E-02
           -9.35542E-03  9.35542E-03  2.97545E-02
        """
        cube_vals = np.array([float(x) for x in raw_vals.split()])
        cube_vals = cube_vals.reshape(nstates * 2, -1)

        self.run_nto(inpfile, eig_vals, nto_vals, cube_vals)


if __name__ == "__main__":
    unittest.main()
