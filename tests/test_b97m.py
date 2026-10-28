import pytest

from veloxchem.dispersionmodel import DispersionModel
from veloxchem.veloxchemlib import mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver


@pytest.mark.skipif(not DispersionModel.is_available(),
                    reason='dftd4-python not available')
@pytest.mark.solvers
class TestB97M:

    def run_scf(self, xcfun_label, ref_scf_energy, tol):

        mol_string = """
        N   -3.710    3.019   -0.037
        H   -3.702    4.942    0.059
        H   -4.704    2.415    1.497
        H   -4.780    2.569   -1.573
        C   -1.621   -5.080    0.444
        H   -0.819   -6.698   -0.465
        H   -3.412   -4.654   -0.393
        H   -0.381   -3.498    0.222
        H   -1.872   -5.468    2.413
        """
        basis_label = 'def2-svp'

        mol = Molecule.read_molecule_string(mol_string, units='bohr')
        bas = MolecularBasis.read(mol, basis_label, ostream=None)

        scf_drv = ScfRestrictedDriver()
        scf_drv.ostream.mute()
        scf_drv.xcfun = xcfun_label
        scf_results = scf_drv.compute(mol, bas)

        if scf_drv.rank == mpi_master():
            assert abs(ref_scf_energy - scf_results['scf_energy']) < tol

    def test_b97m_d4(self):

        self.run_scf('b97m-d4', -97.0598213335, 1.0e-5)

    def test_wb97m_d4(self):

        self.run_scf('wb97m-d4', -96.9982843462, 1.0e-5)
