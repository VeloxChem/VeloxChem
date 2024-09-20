from veloxchem import Molecule
from veloxchem import MolecularBasis
from veloxchem import ScfRestrictedDriver
from veloxchem import C6Driver


class TestC6:

    def run_c6(self, xcfun_label, ref_c6, tol):

        xyz_string = """3
        xyz
        O   -0.1858140  -1.1749469   0.7662596
        H   -0.1285513  -0.8984365   1.6808606
        H   -0.0582782  -0.3702550   0.2638279
        """
        mol = Molecule.read_xyz_string(xyz_string)
        mol.check_multiplicity()

        basis_label = 'def2-svp'
        bas = MolecularBasis.read(mol, basis_label, ostream=None)

        scf_drv = ScfRestrictedDriver()
        scf_drv.ostream.mute()
        scf_drv.xcfun = xcfun_label
        scf_results = scf_drv.compute(mol, bas)

        lr_drv = C6Driver()
        lr_drv.ostream.mute()
        lr_results = lr_drv.compute(mol, bas, scf_results)

        assert abs(lr_results['c6'] - ref_c6) < tol

    def test_hf(self):

        xcfun_label = 'hf'
        ref_c6 = 16.852249

        self.run_c6(xcfun_label, ref_c6, 1.0e-6)

    def test_b3lyp(self):

        xcfun_label = 'b3lyp'
        ref_c6 = 17.421351

        self.run_c6(xcfun_label, ref_c6, 1.0e-4)
