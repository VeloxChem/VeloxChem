import pytest

from veloxchem import mpi_master
from veloxchem import Molecule, MolecularBasis
from veloxchem import ScfUnrestrictedDriver


@pytest.mark.solvers
class TestScfUnrestrictedDriver:

    def run_scf_unrest(self, xcfun_label, charge, mult, ref_scf_energy, tol):

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

        scf_drv = ScfUnrestrictedDriver()
        scf_drv.ostream.mute()
        scf_drv.xcfun = xcfun_label
        scf_results = scf_drv.compute(mol, bas, min_bas)

        if scf_drv.rank == mpi_master():
            assert abs(ref_scf_energy - scf_results['scf_energy']) < tol

    def test_hf(self):

        self.run_scf_unrest('hf', 1, 2, -75.5622230223, 1.0e-8)
        self.run_scf_unrest('hf', 0, 3, -75.7122812456, 1.0e-8)

    def test_slda(self):

        self.run_scf_unrest('slda', 1, 2, -75.5030096846, 1.0e-6)
        self.run_scf_unrest('slda', 0, 3, -75.7017434080, 1.0e-6)

    def test_b3lyp(self):

        self.run_scf_unrest('b3lyp', 1, 2, -75.9017886760, 1.0e-6)
        self.run_scf_unrest('b3lyp', 0, 3, -76.0832193747, 1.0e-6)

    def test_camb3lyp(self):

        self.run_scf_unrest('cam-b3lyp', 1, 2, -75.8737059513, 1.0e-6)
        self.run_scf_unrest('cam-b3lyp', 0, 3, -76.0523450131, 1.0e-6)

    def test_tpssh(self):

        self.run_scf_unrest('tpssh', 1, 2, -75.9044411187, 1.0e-6)
        self.run_scf_unrest('tpssh', 0, 3, -76.0780524011, 1.0e-6)
