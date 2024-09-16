from veloxchem import mpi_master
from veloxchem import Molecule, MolecularBasis
from veloxchem import ScfRestrictedDriver


class TestScfRestrictedDriver:

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
        mol.check_multiplicity()

        bas = MolecularBasis.read(mol, basis_label, ostream=None)
        min_bas = MolecularBasis.read(mol, 'ao-start-guess', ostream=None)

        scf_drv = ScfRestrictedDriver()
        scf_drv.ostream.mute()
        #scf_drv.xcfun = xcfun_label
        scf_results = scf_drv.compute(mol, bas, min_bas)

        if scf_drv.rank == mpi_master():
            assert abs(ref_scf_energy - scf_results['scf_energy']) < tol

    def test_hf(self):

        self.run_scf_rest('hf', -151.9220860932, 1.0e-8)

    """
    def test_slda(self):

        self.run_scf_rest('slda', -151.9813034527, 1.0e-6)

    def test_b3lyp(self):

        self.run_scf_rest('b3lyp', -152.7162995916, 1.0e-6)

    def test_tpssh(self):

        self.run_scf_rest('tpssh', -152.7062253947, 1.0e-6)
    """
