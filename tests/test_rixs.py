import numpy as np
import pytest

from veloxchem.veloxchemlib import mpi_master, hartree_in_ev
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.lreigensolver import LinearResponseEigenSolver
from veloxchem.tdaeigensolver import TdaEigenSolver
from veloxchem.rixsdriver import RixsDriver

@pytest.mark.solvers
class TestRIXS:

    def run_rixs(self, xcfun_label, basis_label, ref_xsection, ncore, nstates,
                  ncorestates=None, nvir=None, nvalence=None, cvs_scf=False, tda=False, tol=1e-6, cutoff_ene=None):

        xyz_string = """3
        C2v
        O    0.000000000000        0.000000000000        0.000000000000
        H    0.000000000000        0.740848095288        0.582094932012
        H    0.000000000000       -0.740848095288        0.582094932012
        """
        mol = Molecule.read_xyz_string(xyz_string)

        bas = MolecularBasis.read(mol, basis_label, ostream=None)

        scf_drv = ScfRestrictedDriver()
        scf_drv.ostream.mute()
        scf_drv.xcfun = xcfun_label
        scf_results = scf_drv.compute(mol, bas)

        lr_drv = LinearResponseEigenSolver()
        if tda:
            lr_drv = TdaEigenSolver()
        lr_drv.ostream.mute()

        if nvir is not None:
            lr_drv.restricted_subspace = True
            lr_drv.num_core_orbitals = ncore
            lr_drv.num_valence_orbitals = nvalence
            lr_drv.num_vir_orbitals = nvir
        else:
            cvs_lr_drv = LinearResponseEigenSolver()
            if tda:
                cvs_lr_drv = TdaEigenSolver()
            cvs_lr_drv.ostream.mute()
            cvs_lr_drv.core_excitation = True
            cvs_lr_drv.num_core_orbitals = ncore
            cvs_lr_drv.nstates = ncorestates
            if cvs_scf is not None:
                scf_drv = ScfRestrictedDriver()
                scf_drv.ostream.mute()
                scf_drv.xcfun = xcfun_label
                cvs_scf_res = scf_drv.compute(mol, bas)
                cvs_rsp_res = cvs_lr_drv.compute(mol, bas, cvs_scf_res)
            else:
                cvs_rsp_res = cvs_lr_drv.compute(mol, bas, scf_results)

        lr_drv.nstates = nstates
        valence_rsp = lr_drv.compute(mol, bas, scf_results)

        rixs_drv = RixsDriver()
        rixs_drv.ostream.mute()
        rixs_drv.final_state_cutoff = cutoff_ene
        rixs_drv.gamma = .16 / hartree_in_ev()

        if nvir is not None:
            rixs_res = rixs_drv.compute(mol, bas, scf_results, valence_rsp)
        elif cvs_scf is not None:
            rixs_res = rixs_drv.compute(mol, bas, scf_results, valence_rsp, cvs_rsp_res, 
                                        cvs_scf_tensors=cvs_scf_res)
        else:
            rixs_res = rixs_drv.compute(mol, bas, scf_results, valence_rsp, cvs_rsp_res)

        if scf_drv.rank == mpi_master():
            assert np.allclose(ref_xsection, rixs_res['cross_sections'][:,0], 
                               rtol=tol, atol=1e-8)

            #assert np.max(np.abs((ref_xsection -
            #                     rixs_res['cross_sections'][:,0]) / ref_xsection)) < tol

    def test_hf_svp_rpa_rsa(self):

        ref_xsection = np.array(
            [17793.04602477,  1230.14877966,    48.72583595,   733.62131704,
            35.28746827,    40.83373707,     0.10615352,   103.18690689,
            0.24773989,     2.04570504,     0.04773057,     0.03343511,
            3.89209812,    75.92903686,     3.90556953,     4.827745  ])

        self.run_rixs('hf', 'def2-svp', ref_xsection, 1, 20, ncorestates=2, nvir=16, nvalence=1, tda=False)

    def test_hf_svp_tda_rsa(self):

        ref_xsection = np.array(
            [18466.4904454 ,  1272.83056895,    50.42815731,   768.62422129,
            34.60402595,    42.25861447,     0.12938393,   108.89842261,
            0.23551592,     2.13906708,     0.06140197,     0.03432167,
            4.02711258,    80.49203813,     4.25473883,     5.27559603])

        self.run_rixs('hf', 'def2-svp', ref_xsection, 1, 20, ncorestates=2, nvir=16, nvalence=1, tda=True)

    def test_b3lyp_svp_rpa_rsa(self):

        ref_xsection = np.array(
            [5555.40347227,   29.49558484,    0.48643769,    4.6474659 ,
            0.3588067 ,    0.00478371,    0.05062501,    0.40879286,
            0.05613822,    0.1177419 ,    0.13567176,    0.00000012,
            0.00033592,    3.67154464,    1.92188469,    2.48548861])

        self.run_rixs('b3lyp', 'def2-svp', ref_xsection, 1, 20, ncorestates=2, nvir=16, nvalence=1, tda=False)

    def test_b3lyp_svp_tda_rsa(self):

        ref_xsection = np.array(
            [5792.84021556,   30.57909736,    0.50228254,    5.27299334,
            0.47202892,    0.00489334,    0.46383182,    0.06376404,
            0.06768143,    0.1384723 ,    0.15317757,    0.0000002 ,
            0.0003445 ,    4.06725417,    2.09146254,    2.7020912 ])

        self.run_rixs('b3lyp', 'def2-svp', ref_xsection, 1, 20, ncorestates=2, nvir=16, nvalence=1, tda=True)
    
    def test_hf_svp_rpa_2s(self):

        ref_xsection = np.array(
            [18243.92834424,  1293.71325113,  6818.74977615,  1027.58810983,
            11761.96369905,   374.80028244,    50.41562301,   783.26769218,
            65.95061377,   241.32429406,    26.12203323,    58.8352719 ,
            166.19299612,   603.14558887,    74.82995014,     0.70763516,
            51.23288065,    43.19207675,    58.91889609,    83.48913356])

        self.run_rixs('hf', 'def2-svp', ref_xsection, 1, 20, ncorestates=4, tda=False)

    def test_hf_svp_tda_2s(self):

        ref_xsection = np.array(
            [18282.66549552,  1296.45917462,  6785.95670576,  1018.28013074,
            11792.29673211,   365.81615569,    51.79197945,   814.66136794,
            70.87364587,   239.1480645 ,    26.5327613 ,    59.12325221,
            170.50239741,   605.35459475,    67.23463214,     1.56059   ,
            59.50158123,    44.2218087 ,    62.95686118,    83.19752305])

        self.run_rixs('hf', 'def2-svp', ref_xsection, 1, 20, ncorestates=4, tda=True)

    def test_b3lyp_svp_rpa_2s(self):

        ref_xsection = np.array(
            [5419.19302002,   29.65392705, 2022.62833796,   50.15397842,
            3601.61140574,   15.47585656,    0.51368214,    4.3825346 ,
            1.93996394,    3.36685515,  134.709674  ,    0.53869339,
            2.63183346,    0.00857954,    7.46612329,    4.98026168,
            0.69017856,    0.00373899,    0.4296243 ,    0.00354383])
        
        self.run_rixs('b3lyp', 'def2-svp', ref_xsection, 1, 20, ncorestates=4, tda=False)

    def test_b3lyp_svp_tda_2s(self):

        ref_xsection = np.array(
            [5451.53264374,   29.79591156, 2015.46547064,   55.19705342,
            3613.00828932,   16.49129082,    0.51804846,    5.17425994,
            1.97078415,    3.96950546,  134.72864448,    3.0797107 ,
            0.40125446,    0.03243186,    7.81557628,    5.54322374,
            0.83284497,    0.0037602 ,    0.48287783,    0.0041094 ])

        self.run_rixs('b3lyp', 'def2-svp', ref_xsection, 1, 20, ncorestates=4, tda=True)
        
    def test_hf_svp_rpa_fulldiag(self):

        ref_xsection = np.array(
            [22977.66629182,   868.83776692,  8780.25926577,   920.42698556,
            14842.72142202,   277.51218997,     6.06482292,   846.58668506,
            106.15901828])

        self.run_rixs('hf', '6-31G*', ref_xsection, 1, 57, nvir=13, nvalence=4, tda=False, cutoff_ene=31.5 / hartree_in_ev())

