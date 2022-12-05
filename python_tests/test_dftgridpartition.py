from pathlib import Path
import numpy as np
import h5py

from veloxchem.veloxchemlib import GridDriver, XCNewIntegrator
from veloxchem.veloxchemlib import is_mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver


class TestDftGridPartition:

    def run_dft_grid_partition(self, mol_str, basis_label, xcfun_label,
                               grid_level, tol, ref_h5file):

        molecule = Molecule.read_str(mol_str, units='angstrom')
        basis = MolecularBasis.read(molecule, basis_label, ostream=None)

        scf_drv = ScfRestrictedDriver()
        scf_drv.xcfun = xcfun_label
        scf_drv.grid_level = grid_level
        scf_drv.ostream.state = False
        scf_drv.compute(molecule, basis)
        density = scf_drv.density

        grid_drv = GridDriver()
        grid_drv.set_level(grid_level)
        mol_grid = grid_drv.generate(molecule)

        xc_drv = XCNewIntegrator()
        mol_grid.partition_grid_points()
        mol_grid.distribute_counts_and_displacements(scf_drv.rank,
                                                     scf_drv.nodes,
                                                     scf_drv.comm)
        vxc = xc_drv.integrate_vxc_fock(molecule, basis, density, mol_grid,
                                        xcfun_label)
        vxc.reduce_sum(scf_drv.rank, scf_drv.nodes, scf_drv.comm)

        if is_mpi_master():
            here = Path(__file__).parent
            h5file = here / 'inputs' / ref_h5file
            hf = h5py.File(h5file, 'r')
            vxc_mat_ref = np.array(hf.get('vxc_mat'))
            hf.close()

            vxc_mat = vxc.get_matrix().to_numpy()
            assert np.max(np.abs(vxc_mat - vxc_mat_ref)) < tol

    def test_slater(self):

        mol_str = """
            O  0.0000000000   0.0000000000  -0.0254395383
            H  0.0000000000   0.7695699584   0.5948147012
            H  0.0000000000  -0.7695699584   0.5948147012
            O  4.0000000000   0.0000000000  -0.0254395383
            H  4.0000000000   0.7695699584   0.5948147012
            H  4.0000000000  -0.7695699584   0.5948147012
        """
        basis_label = 'def2-svp'
        xcfun_label = 'slater'
        grid_level = 4
        tol = 1.0e-11
        ref_h5file = 'vxc_lda.h5'

        self.run_dft_grid_partition(mol_str, basis_label, xcfun_label,
                                    grid_level, tol, ref_h5file)

    def test_b3lyp(self):

        mol_str = """
            O  0.0000000000   0.0000000000  -0.0254395383
            H  0.0000000000   0.7695699584   0.5948147012
            H  0.0000000000  -0.7695699584   0.5948147012
            O  4.0000000000   0.0000000000  -0.0254395383
            H  4.0000000000   0.7695699584   0.5948147012
            H  4.0000000000  -0.7695699584   0.5948147012
        """
        basis_label = 'def2-svp'
        xcfun_label = 'b3lyp'
        grid_level = 4
        tol = 1.0e-10
        ref_h5file = 'vxc_gga.h5'

        self.run_dft_grid_partition(mol_str, basis_label, xcfun_label,
                                    grid_level, tol, ref_h5file)
