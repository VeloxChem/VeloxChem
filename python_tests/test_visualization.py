from pathlib import Path
import numpy as np
import tempfile

from veloxchem.veloxchemlib import is_mpi_master
from veloxchem.cubicgrid import CubicGrid
from veloxchem.visualizationdriver import VisualizationDriver
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.mpitask import MpiTask
from veloxchem.molecule import Molecule


class TestVisualization:

    def test_visualization_driver(self):

        here = Path(__file__).parent
        inpfile = str(here / 'inputs' / 'h2se.inp')

        task = MpiTask([inpfile, None])
        scf_drv = ScfRestrictedDriver(task.mpi_comm, task.ostream)

        scf_drv.compute(task.molecule, task.ao_basis, task.min_basis)
        mol_orbs = scf_drv.mol_orbs
        density = scf_drv.density

        mol_orbs.broadcast(task.mpi_rank, task.mpi_comm)
        density.broadcast(task.mpi_rank, task.mpi_comm)

        grid = CubicGrid([0.3, 0.6, 0.9], [1.0, 1.0, 1.0], [2, 3, 3])
        homo = task.molecule.number_of_alpha_electrons() - 1

        vis_drv = VisualizationDriver(task.mpi_comm)
        vis_drv.compute(grid, task.molecule, task.ao_basis, mol_orbs, homo,
                        'alpha')

        points = [[0.3 + 1.0 * ix, 0.6 + 1.0 * iy, 0.9 + 1.0 * iz]
                  for ix in range(2) for iy in range(3) for iz in range(3)]
        mo_val = vis_drv.get_mo(points, task.molecule, task.ao_basis, mol_orbs,
                                homo, 'alpha')

        if is_mpi_master(task.mpi_comm):
            homo_values = grid.values_to_numpy()

            homo_ref = np.array([
                5.09744e-02, 2.82845e-02, 1.20927e-02, 3.22233e-02, 1.53738e-02,
                6.37159e-03, 1.33161e-02, 6.65229e-03, 3.03229e-03, 1.59826e-01,
                7.64720e-02, 3.17590e-02, 8.77249e-02, 4.56304e-02, 1.99448e-02,
                3.70652e-02, 2.14514e-02, 1.03838e-02
            ]).reshape(2, 3, 3)

            homo_diff = np.abs(homo_values) - np.abs(homo_ref)
            homo_diff_rel = homo_diff / np.abs(homo_ref)
            assert np.max(homo_diff_rel) < 1.0e-5

            mo_val = np.array(mo_val).reshape(2, 3, 3)
            assert np.max(np.abs(np.abs(homo_values) - np.abs(mo_val))) < 1.0e-8

        vis_drv.compute(grid, task.molecule, task.ao_basis, density, 0, 'alpha')

        den_val = vis_drv.get_density(points, task.molecule, task.ao_basis,
                                      density, 'alpha')

        if is_mpi_master(task.mpi_comm):
            dens_alpha = grid.values_to_numpy()
            dens_total = dens_alpha * 2.0

            dens_ref = np.array([
                3.26043e-01, 1.18974e-01, 1.00508e-01, 9.37167e-02, 3.30301e-02,
                1.66540e-02, 7.19205e-02, 1.21828e-02, 3.05223e-03, 1.01686e-01,
                4.68902e-02, 2.50173e-02, 4.55977e-02, 1.69215e-02, 7.68553e-03,
                2.35932e-02, 6.38199e-03, 1.90608e-03
            ]).reshape(2, 3, 3)

            dens_diff = dens_total - dens_ref
            dens_diff_rel = dens_diff / np.abs(dens_ref)
            assert np.max(dens_diff_rel) < 1.0e-5

            den_val = np.array(den_val).reshape(2, 3, 3)
            assert np.max(np.abs(dens_alpha - den_val)) < 1.0e-8

        twoe_val_aa = vis_drv.get_two_particle_density(points, points,
                                                       task.molecule,
                                                       task.ao_basis, density,
                                                       'alpha', 'alpha')

        twoe_val_ab = vis_drv.get_two_particle_density(points, points,
                                                       task.molecule,
                                                       task.ao_basis, density,
                                                       'alpha', 'beta')

        if is_mpi_master(task.mpi_comm):
            twoe_val_aa = np.array(twoe_val_aa).reshape(2, 3, 3)
            twoe_val_ab = np.array(twoe_val_ab).reshape(2, 3, 3)

            assert np.max(np.abs(twoe_val_aa) < 1.0e-8)
            assert np.max(np.abs(twoe_val_ab - den_val**2) < 1.0e-8)

        task.finish()

    def test_cubic_grid(self):

        origin = [0.1, 0.2, 0.3]
        step_size = [1.0, 1.2, 1.5]
        num_points = [2, 3, 5]

        grid = CubicGrid(origin, step_size, num_points)
        data = list(np.random.rand(2 * 3 * 5))
        grid.set_values(data)

        assert origin == grid.get_origin()
        assert step_size == grid.get_step_size()
        assert num_points == grid.get_num_points()
        assert data == list(grid.values_to_numpy().flatten())

    def test_gen_cubic_grid(self):

        mol_str = """
            O   0.0   0.0   0.0
            H   0.0   1.4   1.1
            H   0.0  -1.4   1.1
        """
        mol = Molecule.read_str(mol_str, units='bohr')
        num_points = [2, 3, 5]

        vis_drv = VisualizationDriver()
        grid = vis_drv.gen_cubic_grid(mol, num_points)
        assert num_points == grid.get_num_points()

    def test_gen_cubes(self):

        here = Path(__file__).parent
        inpfile = str(here / 'inputs' / 'water.inp')

        task = MpiTask([inpfile, None])
        scf_drv = ScfRestrictedDriver(task.mpi_comm, task.ostream)

        scf_drv.compute(task.molecule, task.ao_basis, task.min_basis)
        mol_orbs = scf_drv.mol_orbs
        density = scf_drv.density

        mol_orbs.broadcast(task.mpi_rank, task.mpi_comm)
        density.broadcast(task.mpi_rank, task.mpi_comm)

        with tempfile.TemporaryDirectory() as temp_dir:
            dens_cube_fname = str(Path(temp_dir, 'density.cube'))
            homo_cube_fname = str(Path(temp_dir, 'homo.cube'))

            cube_dict = {
                'grid': '2, 3, 5',
                'cubes': 'density(alpha), mo(homo)',
                'files': f'{dens_cube_fname}, {homo_cube_fname}',
            }

            vis_drv = VisualizationDriver()
            vis_drv.gen_cubes(cube_dict, task.molecule, task.ao_basis, mol_orbs,
                              density)

            cubic_grid = vis_drv.gen_cubic_grid(task.molecule, [2, 3, 5])

            vis_drv.compute(cubic_grid, task.molecule, task.ao_basis, density,
                            0, 'alpha')
            if is_mpi_master(task.mpi_comm):
                read_grid = CubicGrid.read_cube(dens_cube_fname)
                assert read_grid.compare(cubic_grid) < 1e-6

            vis_drv.compute(cubic_grid, task.molecule, task.ao_basis, mol_orbs,
                            task.molecule.number_of_alpha_electrons() - 1,
                            'alpha')
            if is_mpi_master(task.mpi_comm):
                read_grid = CubicGrid.read_cube(homo_cube_fname)
                assert read_grid.compare(cubic_grid) < 1e-6
