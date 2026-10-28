import numpy as np
import pytest

from veloxchem.veloxchemlib import mpi_master, hartree_in_ev, fine_structure_constant
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.rixsdriver import RixsDriver


@pytest.mark.solvers
class TestRIXS:

    def run_rixs(self, xcfun_label, basis_label, ref_xsection, ncore, nstates,
                 ncorestates=None, nvir=None, nvalence=None, tda=False, cutoff_ene=None):

        xyz_string = """3
        xyz
        O    0.000000000000        0.000000000000        0.000000000000
        H    0.000000000000        0.740848095288        0.582094932012
        H    0.000000000000       -0.740848095288        0.582094932012
        """
        mol = Molecule.read_xyz_string(xyz_string)

        bas = MolecularBasis.read(mol, basis_label, ostream=None)

        scf_drv = ScfRestrictedDriver()
        scf_drv.ostream.mute()
        scf_drv.xcfun = xcfun_label
        scf_drv.conv_thresh = 1e-8
        scf_results = scf_drv.compute(mol, bas)

        rixs_drv = RixsDriver()
        rixs_drv.ostream.mute()
        rixs_drv.final_state_cutoff = cutoff_ene
        rixs_drv.gamma = .16 / hartree_in_ev()

        rixs_drv.conv_thresh = 1e-6
        rixs_drv.tamm_dancoff = tda
        rixs_drv.nstates = nstates

        if nvir is not None:
            rixs_drv.restricted_subspace = True
            rixs_drv.num_core_orbitals = ncore
            rixs_drv.num_valence_orbitals = nvalence
            rixs_drv.num_virtual_orbitals = nvir
        else:
            rixs_drv.restricted_subspace = False
            rixs_drv.num_core_orbitals = ncore
            rixs_drv.num_core_states = ncorestates

        rixs_res = rixs_drv.compute(mol, bas, scf_results)

        if scf_drv.rank == mpi_master():
            ref_xsection = np.array(ref_xsection) * fine_structure_constant()**4
            assert np.allclose(ref_xsection, rixs_res['cross_sections'][:,0],
                               rtol=1e-5, atol=1e-8)

    def test_hf_svp_rpa_rsa(self):

        ref_xsection = np.array([
            17811.632987086, 1229.758519817, 48.681891774, 733.254590340,
            35.062883474, 40.828959002, 0.105766349, 104.758173889,
            0.269747285, 2.187183719, 0.044225456, 0.033239693,
            3.890917736, 75.817596511, 3.899221697, 4.796194695
        ])

        self.run_rixs('hf', 'def2-svp', ref_xsection, 1, 20, ncorestates=2, nvir=16, nvalence=1, tda=False)

    def test_hf_svp_tda_rsa(self):

        ref_xsection = np.array([
            18466.533352055, 1272.812377437, 50.42744885, 768.626586751,
            34.604271141, 42.258008276, 0.129382754, 108.898665032, 0.235516419,
            2.139067077, 0.061402242, 0.034321189, 4.027053888, 80.492246807,
            4.254749669, 5.275604059
        ])

        self.run_rixs('hf', 'def2-svp', ref_xsection, 1, 20, ncorestates=2, nvir=16, nvalence=1, tda=True)

    def test_b3lyp_svp_rpa_rsa(self):

        ref_xsection = np.array([
            5558.473364257, 29.490542895, 0.486049704, 4.629266414, 0.363044685,
            0.004777157, 0.051060100, 0.440354664, 0.051898838, 0.106438494,
            0.133616609, 0.000000149, 0.000335247, 3.649828944, 1.915088149,
            2.465981090
        ])

        self.run_rixs('b3lyp', 'def2-svp', ref_xsection, 1, 20, ncorestates=2, nvir=16, nvalence=1, tda=False)

    def test_b3lyp_svp_tda_rsa(self):

        ref_xsection = np.array([
            5792.840577899, 30.579097754, 0.502282504, 5.272993513, 0.47202898,
            0.004893339, 0.463831817, 0.063764077, 0.067681437, 0.138472312,
            0.153177591, 0.000000204, 0.000344496, 4.067254404, 2.091462752,
            2.702091265
        ])

        self.run_rixs('b3lyp', 'def2-svp', ref_xsection, 1, 20, ncorestates=2, nvir=16, nvalence=1, tda=True)

    def test_hf_svp_rpa_2s(self):

        ref_xsection = np.array([
            18244.555308374, 1293.746642124, 6818.769734878, 1027.656605437,
            11762.205117774, 374.813155636, 50.422314315, 783.286460299,
            65.956114162, 241.283720267, 26.124851954, 58.892406450,
            166.088961687, 603.195221475, 74.821212081, 0.718789606,
            51.308249625, 43.198426438, 58.966168215, 83.452659794
        ])

        self.run_rixs('hf', 'def2-svp', ref_xsection, 1, 20, ncorestates=4, tda=False)

    def test_hf_svp_tda_2s(self):

        ref_xsection = np.array([
            18282.708141936, 1296.440100667, 6785.96485048, 1018.264720384,
            11792.323854124, 365.810941961, 51.791306589, 814.663876348,
            70.87266473, 239.148738953, 26.53241798, 59.123619902,
            170.502237935, 605.356193762, 67.234797706, 1.560585874,
            59.501431469, 44.226334207, 62.956362303, 83.197727132
        ])

        self.run_rixs('hf', 'def2-svp', ref_xsection, 1, 20, ncorestates=4, tda=True)

    def test_b3lyp_svp_rpa_2s(self):

        ref_xsection = np.array([
            5419.408540219, 29.654485681, 2022.737567512, 50.138868853,
            3601.703649966, 15.488104351, 0.513719573, 4.381978513, 1.945222941,
            3.365520699, 134.695043955, 0.540695855, 2.625666460, 0.008474128,
            7.481248883, 4.981706022, 0.690424900, 0.003741214, 0.428471086,
            0.003525553
        ])

        self.run_rixs('b3lyp', 'def2-svp', ref_xsection, 1, 20, ncorestates=4, tda=False)

    def test_b3lyp_svp_tda_2s(self):

        ref_xsection = np.array([
            5451.525464935, 29.795911911, 2015.46273501, 55.197006007,
            3613.003640347, 16.491281722, 0.518048422, 5.17424475, 1.970785575,
            3.969506943, 134.728441148, 3.079707625, 0.401253864, 0.032427692,
            7.815563849, 5.543194369, 0.832842604, 0.003760196, 0.482878793,
            0.004104707
        ])

        self.run_rixs('b3lyp', 'def2-svp', ref_xsection, 1, 20, ncorestates=4, tda=True)

    def test_hf_svp_rpa_fulldiag(self):

        ref_xsection = np.array([
            23060.546141723, 872.799710229, 8798.472214577, 921.838965056,
            14909.152754469, 278.367031877, 6.099015095, 858.837996761,
            106.762584488
        ])

        self.run_rixs('hf', '6-31G*', ref_xsection, 1, 57, nvir=13, nvalence=4, tda=False, cutoff_ene=31.5 / hartree_in_ev())
