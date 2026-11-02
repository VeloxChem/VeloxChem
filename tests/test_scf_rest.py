import pytest

from veloxchem.veloxchemlib import mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver


@pytest.mark.solvers
class TestScfRestrictedDriver:

    def run_scf_rest(self, xcfun_label, ri, ref_scf_energy, tol):

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

        scf_drv = ScfRestrictedDriver()
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

        ref_scf_energy = {False: -151.9220860932, True: -151.9219743919}[ri]
        self.run_scf_rest('hf', ri, ref_scf_energy, 1.0e-8)

    @pytest.mark.parametrize('ri', [False, True])
    def test_slda(self, ri):

        ref_scf_energy = {False: -151.9813034527, True: -151.9814807155}[ri]
        self.run_scf_rest('slda', ri, ref_scf_energy, 1.0e-6)

    @pytest.mark.parametrize('ri', [False, True])
    def test_b3lyp(self, ri):

        ref_scf_energy = {False: -152.7162995916, True: -152.7163199191}[ri]
        self.run_scf_rest('b3lyp', ri, ref_scf_energy, 1.0e-6)

    @pytest.mark.parametrize('ri', [False, True])
    def test_camb3lyp(self, ri):

        ref_scf_energy = {False: -152.6595419840, True: -152.6595638483}[ri]
        self.run_scf_rest('cam-b3lyp', ri, ref_scf_energy, 1.0e-6)

    @pytest.mark.parametrize('ri', [False, True])
    def test_camb3lyp_100(self, ri):

        ref_scf_energy = {False: -152.5671319422, True: -152.5671534603}[ri]
        self.run_scf_rest('cam-b3lyp-100', ri, ref_scf_energy, 1.0e-6)

    @pytest.mark.parametrize('ri', [False, True])
    def test_tpssh(self, ri):

        ref_scf_energy = {False: -152.7062253947, True: -152.7062633670}[ri]
        self.run_scf_rest('tpssh', ri, ref_scf_energy, 1.0e-6)

    def run_scf_rest_with_ecp(self, xcfun_label, ri, ref_scf_energy, tol):

        xyz_string = """2
        xyz
        Ag  0.0  0.0  0.0
        Br  0.0  2.5  0.0
        """
        basis_label = 'def2-svp'

        mol = Molecule.read_xyz_string(xyz_string)
        bas = MolecularBasis.read(mol, basis_label, ostream=None)

        scf_drv = ScfRestrictedDriver()
        scf_drv.ri = ri
        scf_drv.ostream.mute()
        scf_drv.xcfun = xcfun_label
        scf_drv.acc_type = 'l2_diis'
        scf_results = scf_drv.compute(mol, bas)

        if scf_drv.rank == mpi_master():
            if ri:
                assert scf_results['ri_jk'] or scf_results['ri_coulomb']
                assert not scf_results['ri_fallback']
            assert abs(ref_scf_energy - scf_results['scf_energy']) < tol

    @pytest.mark.parametrize('ri', [False, True])
    def test_hf_with_ecp(self, ri):

        ref_scf_energy = {False: -2718.2407866974, True: -2718.2407205568}[ri]
        self.run_scf_rest_with_ecp('hf', ri, ref_scf_energy, 1.0e-8)
