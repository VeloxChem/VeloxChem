import pytest

from veloxchem.veloxchemlib import mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver


@pytest.mark.solvers
class TestCpcmSolvation:

    def run_cpcm_solvation(self, xcfun_label, ref_solv_energy,
                           cpcm_custom_vdw_radii, tol, *, ri):

        xyz_string = """6
        xyz
        O   -0.1858140  -1.1749469   0.7662596
        H   -0.1285513  -0.8984365   1.6808606
        H   -0.0582782  -0.3702550   0.2638279
        O    0.1747051  11.1050002  -0.7244430
        H   -0.5650842  11.3134964  -1.2949455
        H    0.9282185  11.0652990  -1.3134026
        """
        mol = Molecule.read_xyz_string(xyz_string)

        basis_label = 'def2-svp'
        bas = MolecularBasis.read(mol, basis_label, ostream=None)

        scf_drv = ScfRestrictedDriver()
        scf_drv.ri = ri
        scf_drv.xcfun = xcfun_label

        scf_drv.solvation_model = 'cpcm'
        scf_drv.cpcm_grid_per_sphere = (110, 110)
        scf_drv.cpcm_custom_vdw_radii = cpcm_custom_vdw_radii

        scf_drv.ostream.mute()
        scf_results = scf_drv.compute(mol, bas)

        if scf_drv.rank == mpi_master():
            if ri:
                assert scf_results['ri_jk'] or scf_results['ri_coulomb']
                assert not scf_results['ri_fallback']

        if scf_drv.rank == mpi_master():
            assert abs(ref_solv_energy - scf_drv.cpcm_drv.cpcm_epol) < tol

    @pytest.mark.parametrize('ri', [False, True])
    def test_hf(self, ri):

        self.run_cpcm_solvation('hf',
                                {
                                    False: -0.022994741795899817,
                                    True: -0.023000082327860622
                                }[ri],
                                None,
                                1.0e-8,
                                ri=ri)

    @pytest.mark.parametrize('ri', [False, True])
    def test_b3lyp(self, ri):

        self.run_cpcm_solvation('b3lyp',
                                {
                                    False: -0.020249892377218158,
                                    True: -0.020251215480936596
                                }[ri],
                                None,
                                1.0e-6,
                                ri=ri)

    @pytest.mark.parametrize('ri', [False, True])
    def test_b3lyp_custom_radii(self, ri):

        self.run_cpcm_solvation('b3lyp',
                                -0.00214350811184457, ['O', 3.0],
                                1.0e-6,
                                ri=ri)

    def run_cpcm_solvation_with_ecp(self, ref_solv_energy, tol, *, ri):

        xyz_string = """3
        xyz
        Hg   0.00  0.00  0.00
        Cl   0.00  0.00  2.35
        Cl   0.00  0.10 -2.40
        """
        mol = Molecule.read_xyz_string(xyz_string)

        basis_label = 'def2-svp'
        bas = MolecularBasis.read(mol, basis_label, ostream=None)

        scf_drv = ScfRestrictedDriver()
        scf_drv.ri = ri
        scf_drv.solvation_model = 'cpcm'
        scf_drv.ostream.mute()
        scf_results = scf_drv.compute(mol, bas)

        if scf_drv.rank == mpi_master():
            if ri:
                assert scf_results['ri_jk'] or scf_results['ri_coulomb']
                assert not scf_results['ri_fallback']

        if scf_drv.rank == mpi_master():
            assert abs(ref_solv_energy - scf_drv.cpcm_drv.cpcm_epol) < tol

    @pytest.mark.parametrize('ri', [False, True])
    def test_hf_with_ecp(self, ri):

        self.run_cpcm_solvation_with_ecp({
            False: -0.0465435704,
            True: -0.0465441336
        }[ri],
                                         1.0e-8,
                                         ri=ri)
