import pytest

from veloxchem.veloxchemlib import mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.scfunrestdriver import ScfUnrestrictedDriver


@pytest.mark.solvers
class TestScfDriverWithRIJK:

    def run_scf(self, scf_flag, mol, bas, xcfun_label, ref_scf_energy, tol):

        if scf_flag == 'restricted':
            scf_drv = ScfRestrictedDriver()
        elif scf_flag == 'unrestricted':
            scf_drv = ScfUnrestrictedDriver()

        scf_drv.ostream.mute()
        scf_drv.xcfun = xcfun_label
        scf_drv.ri_jk = True
        scf_results = scf_drv.compute(mol, bas)

        if scf_drv.rank == mpi_master():
            assert abs(ref_scf_energy - scf_results['scf_energy']) < tol

    def test_rijk_hf(self):

        xyz_string = """6
        xyz
        H      1.2001      0.0363      0.8431
        C      0.7031      0.0083     -0.1305
        H      0.9877      0.8943     -0.7114
        H      1.0155     -0.8918     -0.6742
        O     -0.6582     -0.0067      0.1730
        H     -1.1326     -0.0311     -0.6482
        """

        mol = Molecule.read_xyz_string(xyz_string)
        bas = MolecularBasis.read(mol, 'def2-svp', ostream=None)

        self.run_scf('restricted', mol, bas, None, -114.954012228, 1.0e-8)

    def test_rijk_hf_openshell(self):

        xyz_string = """5
        xyz
        H      1.2001      0.0363      0.8431
        C      0.7031      0.0083     -0.1305
        H      0.9877      0.8943     -0.7114
        H      1.0155     -0.8918     -0.6742
        O     -0.6582     -0.0067      0.1730
        """

        mol = Molecule.read_xyz_string(xyz_string)
        mol.set_multiplicity(2)

        bas = MolecularBasis.read(mol, 'def2-svp', ostream=None)

        self.run_scf('unrestricted', mol, bas, None, -114.331645181, 1.0e-8)
