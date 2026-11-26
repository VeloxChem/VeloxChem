import pytest

from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.veloxchemlib import mpi_master

pytest.importorskip('py3Dmol')


@pytest.mark.solvers
class TestCpcmVisualization:

    @pytest.mark.parametrize('ri', [False, True])
    def test_cpcm_visualize_grid(self, ri):
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

        scf_drv.solvation_model = 'cpcm'
        scf_drv.cpcm_grid_per_sphere = (110, 110)

        scf_drv.ostream.mute()
        scf_results = scf_drv.compute(mol, bas)

        if scf_drv.rank == mpi_master():
            if ri:
                assert scf_results['ri_jk'] or scf_results['ri_coulomb']
                assert not scf_results['ri_fallback']

        assert scf_drv.cpcm_drv.cpcm_grid is not None
        assert scf_drv.cpcm_drv.cpcm_grid.shape[1] == 6

        scf_drv.cpcm_drv.visualize_cpcm_grid(mol)

    @pytest.mark.parametrize('ri', [False, True])
    def test_cpcm_visualize_cpcm_charges(self, ri):
        pytest.importorskip('matplotlib.pyplot')

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

        scf_drv.solvation_model = 'cpcm'
        scf_drv.cpcm_grid_per_sphere = (110, 110)

        scf_drv.ostream.mute()
        scf_results = scf_drv.compute(mol, bas)

        if scf_drv.rank == mpi_master():
            if ri:
                assert scf_results['ri_jk'] or scf_results['ri_coulomb']
                assert not scf_results['ri_fallback']

        assert scf_drv.cpcm_drv.cpcm_grid is not None
        assert scf_drv.cpcm_drv.cpcm_q is not None
        assert (scf_drv.cpcm_drv.cpcm_q.shape[0] ==
                scf_drv.cpcm_drv.cpcm_grid.shape[0])

        scf_drv.cpcm_drv.visualize_cpcm_charges(mol)
