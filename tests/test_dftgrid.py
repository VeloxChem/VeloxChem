from mpi4py import MPI
from pathlib import Path
import numpy as np

from veloxchem.veloxchemlib import GridDriver, XCIntegrator, DenseMatrix
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis


class TestDftGrid:

    def test_dft_grid(self):

        molstr = """
        O   0.0   0.0   0.0
        H   0.0   1.4   1.1
        H   0.0  -1.4   1.1
        """
        mol = Molecule.read_str(molstr, 'au')

        here = Path(__file__).parent
        basis_path = str(here.parent / 'basis')
        bas = MolecularBasis.read(mol, 'def2-svp', basis_path, ostream=None)

        grid_drv = GridDriver()
        grid_drv.set_level(1)
        mol_grid = grid_drv.generate(mol)

        # test grid and GTO values

        gx = mol_grid.x_to_numpy()
        gy = mol_grid.y_to_numpy()
        gz = mol_grid.z_to_numpy()
        gw = mol_grid.w_to_numpy()

        grid_data = np.vstack((gx, gy, gz, gw))

        xc_drv = XCIntegrator()
        gto_data = xc_drv.compute_gto_values(mol, bas, mol_grid)

        comm = MPI.COMM_WORLD
        gto_data = comm.allreduce(gto_data)
        gto_data_slice = gto_data[:, 400:600]

        here = Path(__file__).parent
        npz_data = np.load(str(here / 'data' / 'h2o.dftgrid.npz'))
        ref_grid_data = npz_data['grid_xyzw']
        ref_gto_data_slice = npz_data['gtos_on_grid_400_600']

        assert np.max(np.abs(grid_data - ref_grid_data)) < 1e-10
        assert np.max(np.abs(gto_data_slice - ref_gto_data_slice)) < 1e-10

        # test KS matrix

        dmat = DenseMatrix(np.load('dens.npz')['density'])
        xcmat = xc_drv.integrate_vxc_fock(mol, bas, dmat, mol_grid,
                                          'closedshell')

        xcmat_np = comm.reduce(xcmat.alpha_to_numpy(), root=0)
        xcmat_energy = comm.reduce(xcmat.get_energy(), root=0)
        xcmat_electrons = comm.reduce(xcmat.get_electrons(), root=0)

        if comm.Get_rank() == 0:
            xcref = np.load('dens.npz')['ks']
            energy_ref = -8.10594893202593
            electrons_ref = 9.996583310662455
            assert np.max(np.abs(xcmat_np - xcref)) < 1.0e-10
            assert abs(xcmat_energy - energy_ref) < 1.0e-10
            assert abs(xcmat_electrons - electrons_ref) < 1.0e-10
