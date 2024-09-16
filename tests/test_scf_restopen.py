import pytest

from veloxchem import mpi_master
from veloxchem import Molecule, MolecularBasis
from veloxchem import ScfRestrictedOpenDriver


@pytest.mark.solvers
class TestScfRestrictedOpenDriver:

    def run_scf_restopen(self, xcfun_label, charge, mult, ref_scf_energy, tol):

        xyz_string = """3
        xyz
        O   -0.1858140  -1.1749469   0.7662596
        H   -0.1285513  -0.8984365   1.6808606
        H   -0.0582782  -0.3702550   0.2638279
        """
        basis_label = 'def2-svp'

        mol = Molecule.read_xyz_string(xyz_string)
        mol.set_charge(charge)
        mol.set_multiplicity(mult)
        mol.check_multiplicity()

        bas = MolecularBasis.read(mol, basis_label, ostream=None)
        min_bas = MolecularBasis.read(mol, 'ao-start-guess', ostream=None)

        scf_drv = ScfRestrictedOpenDriver()
        scf_drv.ostream.mute()
        #scf_drv.xcfun = xcfun_label
        scf_results = scf_drv.compute(mol, bas, min_bas)

        if scf_drv.rank == mpi_master():
            assert abs(ref_scf_energy - scf_results['scf_energy']) < tol

    def test_hf(self):

        self.run_scf_restopen('hf', 1, 2, -75.5576463523, 1.0e-8)
        self.run_scf_restopen('hf', 0, 3, -75.7069341237, 1.0e-8)

    """
    def test_slda(self):

        self.run_scf_restopen('slda', 1, 2, -75.5020555275, 1.0e-6)
        self.run_scf_restopen('slda', 0, 3, -75.7009836953, 1.0e-6)

    def test_b3lyp(self):

        self.run_scf_restopen('b3lyp', 1, 2, -75.9000931681, 1.0e-6)
        self.run_scf_restopen('b3lyp', 0, 3, -76.0815511111, 1.0e-6)

    def test_tpssh(self):

        self.run_scf_restopen('tpssh', 1, 2, -75.9023473162, 1.0e-6)
        self.run_scf_restopen('tpssh', 0, 3, -76.0759108977, 1.0e-6)
    """
