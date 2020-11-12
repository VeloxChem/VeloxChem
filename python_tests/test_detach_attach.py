from mpi4py import MPI
import numpy as np
import unittest
from pathlib import Path

from veloxchem.veloxchemlib import AODensityMatrix
from veloxchem.veloxchemlib import denmat
from veloxchem.veloxchemlib import mpi_master
from veloxchem.mpitask import MpiTask
from veloxchem.cubicgrid import CubicGrid
from veloxchem.visualizationdriver import VisualizationDriver
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.lreigensolver import LinearResponseEigenSolver


class TestDetachAttach(unittest.TestCase):

    def run_detach_attach(self, inpfile, xcfun_label, nstates, ref_cube_vals):

        task = MpiTask([inpfile, None], MPI.COMM_WORLD)
        task.input_dict['scf']['checkpoint_file'] = None
        task.input_dict['method_settings']['xcfun'] = xcfun_label

        scf_drv = ScfRestrictedDriver(task.mpi_comm, task.ostream)
        scf_drv.update_settings(task.input_dict['scf'],
                                task.input_dict['method_settings'])
        scf_drv.compute(task.molecule, task.ao_basis, task.min_basis)

        rsp_drv = LinearResponseEigenSolver(task.mpi_comm, task.ostream)
        rsp_drv.update_settings({'nstates': nstates},
                                task.input_dict['method_settings'])
        rsp_results = rsp_drv.compute(task.molecule, task.ao_basis,
                                      scf_drv.scf_tensors)

        nocc = task.molecule.number_of_alpha_electrons()

        if task.mpi_rank == mpi_master():
            mo = scf_drv.scf_tensors['C']
            mo_occ = mo[:, :nocc]
            mo_vir = mo[:, nocc:]

            eig_vecs = rsp_results['eigenvectors']

            cube_vals = []

        for s in range(nstates):

            if task.mpi_rank == mpi_master():
                eigvec = eig_vecs[:, s].copy()
                half_eigvec_size = eigvec.shape[0] // 2

                eigvec_z = eigvec[:half_eigvec_size].reshape(
                    half_eigvec_size, 1) * np.sqrt(2.0)
                eigvec_y = eigvec[half_eigvec_size:].reshape(
                    half_eigvec_size, 1) * np.sqrt(2.0)

                dens_Dz, dens_Az = rsp_drv.get_detach_attach_densities(
                    0, eigvec_z, mo_occ, mo_vir)
                dens_Dy, dens_Ay = rsp_drv.get_detach_attach_densities(
                    0, eigvec_y, mo_occ, mo_vir)

                dens_DA = AODensityMatrix(
                    [dens_Dz + dens_Dy, dens_Az + dens_Ay], denmat.rest)
            else:
                dens_DA = AODensityMatrix()
            dens_DA.broadcast(task.mpi_rank, task.mpi_comm)

            vis_drv = VisualizationDriver(task.mpi_comm)
            grid = CubicGrid([-0.1, -0.1, -0.1], [0.3, 0.3, 0.2], [2, 2, 3])

            vis_drv.compute(grid, task.molecule, task.ao_basis, dens_DA, 0,
                            'alpha')
            if task.mpi_rank == mpi_master():
                cube_vals.append(grid.values_to_numpy().flatten())

            vis_drv.compute(grid, task.molecule, task.ao_basis, dens_DA, 1,
                            'alpha')
            if task.mpi_rank == mpi_master():
                cube_vals.append(grid.values_to_numpy().flatten())

        if task.mpi_rank == mpi_master():
            cube_vals = np.array(cube_vals)

            for s in range(cube_vals.shape[0]):
                if np.vdot(cube_vals[s, :], ref_cube_vals[s, :]) < 0.0:
                    cube_vals[s, :] *= -1.0
            self.assertTrue(np.max(np.abs(cube_vals - ref_cube_vals)) < 1.0e-4)

    def test_detach_attach_rpa(self):

        here = Path(__file__).parent
        inpfile = str(here / 'inputs' / 'c2h4.inp')

        xcfun_label = 'b3lyp'

        nstates = 3

        raw_vals = """
           -3.27476e-05  -3.27476e-05  -1.32379e-04
           -5.99181e-05  -5.99181e-05  -4.18994e-04
           -8.17569e-05  -8.17569e-05  -1.86632e-04
           -1.05149e-04  -1.05149e-04  -4.50826e-04
            3.41319e-05   3.41319e-05   2.42345e-04
            6.59041e-05   6.59041e-05   2.65393e-04
            9.52284e-05   9.52284e-05   8.65482e-04
            1.21892e-04   1.21892e-04   8.45480e-04
           -3.73966e-03  -3.73967e-03  -4.05024e-03
           -3.59327e-03  -3.59329e-03  -3.87652e-03
           -6.30145e-03  -6.30147e-03  -7.23666e-03
           -6.05426e-03  -6.05429e-03  -6.91549e-03
            1.31340e-04   1.31439e-04   1.11607e-03
            1.36796e-04   1.36866e-04   1.06198e-03
            1.78907e-04   1.79094e-04   1.60605e-03
            1.81161e-04   1.81312e-04   1.52018e-03
           -9.51467e-04  -9.51467e-04  -1.18572e-03
           -9.13013e-04  -9.13013e-04  -1.12955e-03
           -3.58911e-03  -3.58911e-03  -4.45325e-03
           -3.44578e-03  -3.44578e-03  -4.24514e-03
            3.21810e-04   3.21810e-04   5.24184e-04
            2.94778e-04   2.94778e-04   4.61132e-04
            2.73616e-04   2.73616e-04   4.34842e-04
            2.51156e-04   2.51156e-04   3.82716e-04
        """
        cube_vals = np.array([float(x) for x in raw_vals.split()])
        cube_vals = cube_vals.reshape(nstates * 2, -1)

        self.run_detach_attach(inpfile, xcfun_label, nstates, cube_vals)


if __name__ == "__main__":
    unittest.main()
