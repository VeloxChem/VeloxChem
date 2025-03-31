import pytest

from veloxchem.veloxchemlib import mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver


@pytest.mark.solvers
class TestScfRestrictedSemiNumDriver:

    def run_scf_rest(self, xcfun_label, ref_scf_energy, tol):

        xyz_string = """6
        xyz
        O   -0.1858140  -1.1749469   0.7662596
        H   -0.1285513  -0.8984365   1.6808606
        H   -0.0582782  -0.3702550   0.2638279
        O    0.1747051  11.1050002  -0.7244430
        H   -0.5650842  11.3134964  -1.2949455
        H    0.9282185  11.0652990  -1.3134026
        """
        basis_label = 'def2-svp'

        mol = Molecule.read_xyz_string(xyz_string)

        bas = MolecularBasis.read(mol, basis_label, ostream=None)
        min_bas = MolecularBasis.read(mol, 'ao-start-guess', ostream=None)

        scf_drv = ScfRestrictedDriver()
        scf_drv.xcfun = xcfun_label
        #scf_drv.ri_coulomb = True
        #scf_drv.semi_num_exchange = True
        scf_results = scf_drv.compute(mol, bas, min_bas)
        
        assert False

        if scf_drv.rank == mpi_master():
            assert abs(ref_scf_energy - scf_results['scf_energy']) < tol

    def test_b3lyp(self):

        self.run_scf_rest('b3lyp', -152.7162995916, 1.0e-6)
