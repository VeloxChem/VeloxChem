from mpi4py import MPI
import numpy as np
import unittest

from veloxchem.veloxchemlib import CubicGrid
from veloxchem.veloxchemlib import mpi_master
from veloxchem.visualizationdriver import VisualizationDriver
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.mpitask import MpiTask


class TestVisualization(unittest.TestCase):

    def test_visualization_driver(self):

        task = MpiTask(['inputs/h2se.inp', None], MPI.COMM_WORLD)
        scf_drv = ScfRestrictedDriver(task.mpi_comm, task.ostream)

        scf_drv.compute_task(task)
        mol_orbs = scf_drv.mol_orbs
        density = scf_drv.density

        if task.mpi_rank == mpi_master():

            vis_drv = VisualizationDriver()

            grid = CubicGrid([0.3, 0.6, 0.9], [1.0, 1.0, 1.0], [2, 3, 3])

            homo = task.molecule.number_of_alpha_electrons() - 1

            vis_drv.compute(grid, task.molecule, task.ao_basis, mol_orbs, homo,
                            'alpha')

            homo_values = grid.values_to_numpy()

            homo_ref = np.array([
                5.09744e-02, 2.82845e-02, 1.20927e-02, 3.22233e-02, 1.53738e-02,
                6.37159e-03, 1.33161e-02, 6.65229e-03, 3.03229e-03, 1.59826e-01,
                7.64720e-02, 3.17590e-02, 8.77249e-02, 4.56304e-02, 1.99448e-02,
                3.70652e-02, 2.14514e-02, 1.03838e-02
            ]).reshape((2, 3, 3))

            homo_diff = np.abs(homo_values) - np.abs(homo_ref)
            homo_diff_rel = homo_diff / np.abs(homo_ref)

            self.assertTrue(np.max(homo_diff_rel) < 1.0e-5)

            vis_drv.compute(grid, task.molecule, task.ao_basis, density, 0,
                            'alpha')

            dens_alpha = grid.values_to_numpy()
            dens_total = dens_alpha * 2.0

            dens_ref = np.array([
                3.26043e-01, 1.18974e-01, 1.00508e-01, 9.37167e-02, 3.30301e-02,
                1.66540e-02, 7.19205e-02, 1.21828e-02, 3.05223e-03, 1.01686e-01,
                4.68902e-02, 2.50173e-02, 4.55977e-02, 1.69215e-02, 7.68553e-03,
                2.35932e-02, 6.38199e-03, 1.90608e-03
            ]).reshape((2, 3, 3))

            dens_diff = dens_total - dens_ref
            dens_diff_rel = dens_diff / np.abs(dens_ref)

            self.assertTrue(np.max(dens_diff_rel) < 1.0e-5)

        task.finish()

    def test_cubic_grid(self):

        origin = [0.1, 0.2, 0.3]
        step_size = [1.0, 1.2, 1.5]
        num_points = [2, 3, 5]

        grid = CubicGrid(origin, step_size, num_points)

        self.assertEqual(
            origin,
            [grid.x_origin(), grid.y_origin(),
             grid.z_origin()])

        self.assertEqual(
            step_size,
            [grid.x_step_size(),
             grid.y_step_size(),
             grid.z_step_size()])

        self.assertEqual(
            num_points,
            [grid.x_num_points(),
             grid.y_num_points(),
             grid.z_num_points()])

        self.assertEqual(grid.values_to_numpy().shape, tuple(num_points))


if __name__ == '__main__':
    unittest.main()
