import pytest

from veloxchem.veloxchemlib import mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver


@pytest.mark.solvers
class TestCpcmSolvation:

    def run_cpcm_solvation(self, xcfun_label, ref_solv_energy,
                           cpcm_custom_vdw_radii, tol):

        xyz_string = """6
        xyz
        O   -0.1858140  -1.1749469   0.7662596
        H   -0.1285513  -0.8984365   1.6808606
        H   -0.0582782  -0.3702550   0.2638279
        O    0.1747051  11.1050002  -0.7244430
        H   -0.5650842  11.3134964  -1.2949455
        H    0.9282185  11.0652990  -1.3134026
        """
        mol = Molecule.read_xyz_string(xyz_string)

        basis_label = 'def2-svp'
        bas = MolecularBasis.read(mol, basis_label, ostream=None)

        scf_drv = ScfRestrictedDriver()
        scf_drv.xcfun = xcfun_label

        scf_drv.solvation_model = 'cpcm'
        scf_drv.cpcm_grid_per_sphere = 110
        scf_drv.cpcm_custom_vdw_radii = cpcm_custom_vdw_radii

        scf_drv.ostream.mute()
        scf_drv.compute(mol, bas)

        if scf_drv.rank == mpi_master():
            assert abs(ref_solv_energy - scf_drv.cpcm_epol) < tol

    def test_hf(self):

        self.run_cpcm_solvation('hf', -0.022994741795899817, None, 1.0e-8)

    def test_b3lyp(self):

        self.run_cpcm_solvation('b3lyp', -0.020249892377218158, None, 1.0e-6)

    def test_b3lyp_custom_radii(self):

        self.run_cpcm_solvation('b3lyp', -0.00214350811184457, ['O', 3.0],
                                1.0e-6)
