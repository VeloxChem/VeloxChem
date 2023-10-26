from mpi4py import MPI
from pathlib import Path
import numpy as np

from veloxchem.veloxchemlib import GridDriver, XCIntegrator, AODensityMatrix
from veloxchem.veloxchemlib import denmat
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

        # test KS matrix

        xc_drv = XCIntegrator()

        npz_data = np.load(str(here / 'data' / 'h2o.dens.npz'))
        dmat = AODensityMatrix([npz_data['density']], denmat.rest)
        xcmat = xc_drv.integrate_vxc_fock(mol, bas, dmat, mol_grid, 'slater')

        comm = MPI.COMM_WORLD

        xcmat_np = comm.reduce(xcmat.alpha_to_numpy(), root=0)
        xcmat_energy = comm.reduce(xcmat.get_energy(), root=0)
        xcmat_electrons = comm.reduce(xcmat.get_electrons(), root=0)

        if comm.Get_rank() == 0:
            xcref = npz_data['ks']
            energy_ref = -8.10594893202593
            electrons_ref = 9.996583310662455
            assert np.max(np.abs(xcmat_np - xcref)) < 1.0e-10
            assert abs(xcmat_energy - energy_ref) < 1.0e-10
            assert abs(xcmat_electrons - electrons_ref) < 1.0e-10
