from pathlib import Path
import numpy as np
import h5py

from veloxchem.veloxchemlib import XCMolecularHessian
from veloxchem.veloxchemlib import mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.griddriver import GridDriver
from veloxchem.scfrestdriver import ScfRestrictedDriver


class TestXCMolHess:

    def run_xc_mol_hess(self,
                        molecule,
                        basis,
                        xcfun_label,
                        ref_exc_deriv_2,
                        ref_vxc_deriv_1=None):

        grid_level = 6
        scf_conv_thresh = 1.0e-8

        scf_drv = ScfRestrictedDriver()
        scf_drv.xcfun = xcfun_label
        scf_drv.grid_level = grid_level
        scf_drv.conv_thresh = scf_conv_thresh
        scf_drv.ostream.mute()
        scf_results = scf_drv.compute(molecule, basis)

        if scf_drv.rank == mpi_master():
            gs_dm = scf_results['D_alpha']
        else:
            gs_dm = None
        gs_dm = scf_drv.comm.bcast(gs_dm, root=mpi_master())

        grid_drv = GridDriver()
        grid_drv.set_level(grid_level)

        xc_mol_hess = XCMolecularHessian()
        mol_grid = grid_drv.generate(molecule)

        exc_deriv_2 = xc_mol_hess.integrate_exc_hessian(molecule, basis,
                                                        [gs_dm], mol_grid,
                                                        xcfun_label)
        exc_deriv_2 = scf_drv.comm.reduce(exc_deriv_2, root=mpi_master())

        atom_list = list(range(molecule.number_of_atoms()))
        vxc_deriv_1 = xc_mol_hess.integrate_vxc_fock_gradient(
            molecule, basis, [gs_dm], mol_grid, xcfun_label, atom_list)
        vxc_deriv_1 = np.array(vxc_deriv_1)
        vxc_deriv_1 = scf_drv.comm.reduce(vxc_deriv_1, root=mpi_master())

        if scf_drv.rank == mpi_master():
            vxc_deriv_1 = vxc_deriv_1.reshape(ref_vxc_deriv_1.shape)
            assert np.max(np.abs(ref_exc_deriv_2 - exc_deriv_2)) < 1.0e-4
            assert np.max(np.abs(ref_vxc_deriv_1 - vxc_deriv_1)) < 1.0e-4

    def test_xc_mol_hess_lda(self):

        molecule_string = """
        O   0.0   0.0   0.0
        H   0.0   1.4   1.1
        H   0.0  -1.4   1.1
        O  10.0   0.0   0.0
        H  10.0   1.4   1.1
        H  10.0  -1.4   1.1
        """
        units = 'au'

        basis_set_label = 'def2-svp'
        xcfun_label = 'slater'

        molecule = Molecule.read_molecule_string(molecule_string, units=units)
        basis = MolecularBasis.read(molecule, basis_set_label, ostream=None)

        here = Path(__file__).parent
        h5file = str(here / 'data' / 'hessian_lda_data.h5')
        hf = h5py.File(h5file, 'r')
        ref_exc_deriv_2 = np.array(hf.get('exc_geom_deriv_2'))
        ref_vxc_deriv_1 = np.array(hf.get('vxc_geom_deriv_1'))
        hf.close()

        self.run_xc_mol_hess(molecule, basis, xcfun_label, ref_exc_deriv_2,
                             ref_vxc_deriv_1)

    def test_xc_mol_hess_gga(self):

        molecule_string = """
        O   0.0   0.0   0.0
        H   0.0   1.4   1.1
        H   0.0  -1.4   1.1
        O  10.0   0.0   0.0
        H  10.0   1.4   1.1
        H  10.0  -1.4   1.1
        """
        units = 'au'

        basis_set_label = 'def2-svp'
        xcfun_label = 'blyp'

        molecule = Molecule.read_molecule_string(molecule_string, units=units)
        basis = MolecularBasis.read(molecule, basis_set_label, ostream=None)

        here = Path(__file__).parent
        h5file = str(here / 'data' / 'hessian_gga_data.h5')
        hf = h5py.File(h5file, 'r')
        ref_exc_deriv_2 = np.array(hf.get('exc_geom_deriv_2'))
        ref_vxc_deriv_1 = np.array(hf.get('vxc_geom_deriv_1'))
        hf.close()

        self.run_xc_mol_hess(molecule, basis, xcfun_label, ref_exc_deriv_2,
                             ref_vxc_deriv_1)
