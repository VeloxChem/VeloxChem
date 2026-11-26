import pytest

from veloxchem.veloxchemlib import mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfunrestdriver import ScfUnrestrictedDriver


@pytest.mark.solvers
class TestScfUnrestrictedDriver:

    def run_scf_unrest(self, xcfun_label, charge, mult, ri, ref_scf_energy,
                       tol):

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

        bas = MolecularBasis.read(mol, basis_label, ostream=None)

        scf_drv = ScfUnrestrictedDriver()
        scf_drv.ri = ri
        scf_drv.ostream.mute()
        scf_drv.xcfun = xcfun_label
        scf_results = scf_drv.compute(mol, bas)

        if scf_drv.rank == mpi_master():
            if ri:
                assert scf_results['ri_jk'] or scf_results['ri_coulomb']
                assert not scf_results['ri_fallback']
            assert abs(ref_scf_energy - scf_results['scf_energy']) < tol

    @pytest.mark.parametrize('ri', [False, True])
    def test_hf(self, ri):

        ref_cation, ref_triplet = {
            False: (-75.5622230223, -75.7122812456),
            True: (-75.5621790576, -75.7122289634),
        }[ri]
        self.run_scf_unrest('hf', 1, 2, ri, ref_cation, 1.0e-8)
        self.run_scf_unrest('hf', 0, 3, ri, ref_triplet, 1.0e-8)

    @pytest.mark.parametrize('ri', [False, True])
    def test_slda(self, ri):

        ref_cation, ref_triplet = {
            False: (-75.5030096846, -75.7017434080),
            True: (-75.5031357373, -75.7020123499),
        }[ri]
        self.run_scf_unrest('slda', 1, 2, ri, ref_cation, 1.0e-6)
        self.run_scf_unrest('slda', 0, 3, ri, ref_triplet, 1.0e-6)

    @pytest.mark.parametrize('ri', [False, True])
    def test_b3lyp(self, ri):

        ref_cation, ref_triplet = {
            False: (-75.9017886760, -76.0832193747),
            True: (-75.9018077132, -76.0832413883),
        }[ri]
        self.run_scf_unrest('b3lyp', 1, 2, ri, ref_cation, 1.0e-6)
        self.run_scf_unrest('b3lyp', 0, 3, ri, ref_triplet, 1.0e-6)

    @pytest.mark.parametrize('ri', [False, True])
    def test_camb3lyp(self, ri):

        ref_cation, ref_triplet = {
            False: (-75.8737059513, -76.0523450131),
            True: (-75.8737254924, -76.0523676093),
        }[ri]
        self.run_scf_unrest('cam-b3lyp', 1, 2, ri, ref_cation, 1.0e-6)
        self.run_scf_unrest('cam-b3lyp', 0, 3, ri, ref_triplet, 1.0e-6)

    @pytest.mark.parametrize('ri', [False, True])
    def test_camb3lyp_100(self, ri):

        ref_cation, ref_triplet = {
            False: (-75.8251370791, -76.0032834669),
            True: (-75.8251562567, -76.0033057367),
        }[ri]
        self.run_scf_unrest('cam-b3lyp-100', 1, 2, ri, ref_cation, 1.0e-6)
        self.run_scf_unrest('cam-b3lyp-100', 0, 3, ri, ref_triplet, 1.0e-6)

    @pytest.mark.parametrize('ri', [False, True])
    def test_tpssh(self, ri):

        ref_cation, ref_triplet = {
            False: (-75.9044411187, -76.0780524011),
            True: (-75.9044690766, -76.0780843997),
        }[ri]
        self.run_scf_unrest('tpssh', 1, 2, ri, ref_cation, 1.0e-6)
        self.run_scf_unrest('tpssh', 0, 3, ri, ref_triplet, 1.0e-6)

    def run_scf_unrest_with_ecp(self, xcfun_label, charge, mult, ri,
                                ref_scf_energy, tol):

        xyz_string = """2
        xyz
        Au   0.000   0.000   0.000
        H    0.000   1.550   0.000
        """
        basis_label = 'def2-svp'

        mol = Molecule.read_xyz_string(xyz_string)
        mol.set_charge(charge)
        mol.set_multiplicity(mult)

        bas = MolecularBasis.read(mol, basis_label, ostream=None)

        scf_drv = ScfUnrestrictedDriver()
        scf_drv.ri = ri
        scf_drv.ostream.mute()
        scf_drv.xcfun = xcfun_label
        scf_drv.acc_type = 'diis'
        scf_results = scf_drv.compute(mol, bas)

        if scf_drv.rank == mpi_master():
            if ri:
                assert scf_results['ri_jk'] or scf_results['ri_coulomb']
                assert not scf_results['ri_fallback']
            assert abs(ref_scf_energy - scf_results['scf_energy']) < tol

    @pytest.mark.parametrize('ri', [False, True])
    def test_hf_with_ecp(self, ri):

        ref_scf_energy = {False: -135.0230116586, True: -135.0229944063}[ri]
        self.run_scf_unrest_with_ecp('hf', 1, 2, ri, ref_scf_energy, 1.0e-8)
