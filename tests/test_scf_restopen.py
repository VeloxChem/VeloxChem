import pytest

from veloxchem.veloxchemlib import mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestopendriver import ScfRestrictedOpenDriver


@pytest.mark.solvers
class TestScfRestrictedOpenDriver:

    def run_scf_restopen(self, xcfun_label, charge, mult, ri, ref_scf_energy,
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

        scf_drv = ScfRestrictedOpenDriver()
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
            False: (-75.5576463523, -75.7069341237),
            True: (-75.5576026658, -75.7068811789),
        }[ri]
        self.run_scf_restopen('hf', 1, 2, ri, ref_cation, 1.0e-8)
        self.run_scf_restopen('hf', 0, 3, ri, ref_triplet, 1.0e-8)

    @pytest.mark.parametrize('ri', [False, True])
    def test_slda(self, ri):

        ref_cation, ref_triplet = {
            False: (-75.5020555275, -75.7009836953),
            True: (-75.5021818641, -75.7012523118),
        }[ri]
        self.run_scf_restopen('slda', 1, 2, ri, ref_cation, 1.0e-6)
        self.run_scf_restopen('slda', 0, 3, ri, ref_triplet, 1.0e-6)

    @pytest.mark.parametrize('ri', [False, True])
    def test_b3lyp(self, ri):

        ref_cation, ref_triplet = {
            False: (-75.9000931681, -76.0815511111),
            True: (-75.9001122207, -76.0815730752),
        }[ri]
        self.run_scf_restopen('b3lyp', 1, 2, ri, ref_cation, 1.0e-6)
        self.run_scf_restopen('b3lyp', 0, 3, ri, ref_triplet, 1.0e-6)

    @pytest.mark.parametrize('ri', [False, True])
    def test_camb3lyp(self, ri):

        ref_cation, ref_triplet = {
            False: (-75.8720418597, -76.0506977224),
            True: (-75.8720614175, -76.0507202772),
        }[ri]
        self.run_scf_restopen('cam-b3lyp', 1, 2, ri, ref_cation, 1.0e-6)
        self.run_scf_restopen('cam-b3lyp', 0, 3, ri, ref_triplet, 1.0e-6)

    @pytest.mark.parametrize('ri', [False, True])
    def test_tpssh(self, ri):

        ref_cation, ref_triplet = {
            False: (-75.9023473162, -76.0759108977),
            True: (-75.9023752762, -76.0759429234),
        }[ri]
        self.run_scf_restopen('tpssh', 1, 2, ri, ref_cation, 1.0e-6)
        self.run_scf_restopen('tpssh', 0, 3, ri, ref_triplet, 1.0e-6)

    def run_scf_restopen_with_ecp(self, xcfun_label, charge, mult, ri,
                                  ref_scf_energy, tol):

        xyz_string = """2
        xyz
        H    0.000   1.550   0.000
        Au   0.000   0.000   0.000
        """
        basis_label = 'def2-svp'

        mol = Molecule.read_xyz_string(xyz_string)
        mol.set_charge(charge)
        mol.set_multiplicity(mult)

        bas = MolecularBasis.read(mol, basis_label, ostream=None)

        scf_drv = ScfRestrictedOpenDriver()
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

        ref_scf_energy = {False: -135.0196479978, True: -135.0196316977}[ri]
        self.run_scf_restopen_with_ecp('hf', 1, 2, ri, ref_scf_energy, 1.0e-8)
