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
                  ncorestates=None, nvir=None, nvalence=None, cvs_scf=False, tda=False, cutoff_ene=None):

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

        lr_drv = LinearResponseEigenSolver()
        if tda:
            lr_drv = TdaEigenSolver()
        lr_drv.ostream.mute()

        if nvir is not None:
            lr_drv.restricted_subspace = True
            lr_drv.num_core_orbitals = ncore
            lr_drv.num_valence_orbitals = nvalence
            lr_drv.num_virtual_orbitals = nvir
        else:
            cvs_lr_drv = LinearResponseEigenSolver()
            if tda:
                cvs_lr_drv = TdaEigenSolver()
            cvs_lr_drv.ostream.mute()
            cvs_lr_drv.core_excitation = True
            cvs_lr_drv.num_core_orbitals = ncore
            cvs_lr_drv.nstates = ncorestates
            cvs_lr_drv.conv_thresh = 1e-6
            if cvs_scf:
                cvs_scf_drv = ScfRestrictedDriver()
                cvs_scf_drv.ostream.mute()
                cvs_scf_drv.xcfun = xcfun_label
                cvs_scf_drv.conv_thresh = 1e-8
                cvs_scf_res = cvs_scf_drv.compute(mol, bas)
                cvs_rsp_res = cvs_lr_drv.compute(mol, bas, cvs_scf_res)
            else:
                cvs_rsp_res = cvs_lr_drv.compute(mol, bas, scf_results)

        lr_drv.nstates = nstates
        lr_drv.conv_thresh = 1e-6
        valence_rsp = lr_drv.compute(mol, bas, scf_results)

        rixs_drv = RixsDriver()
        rixs_drv.ostream.mute()
        rixs_drv.final_state_cutoff = cutoff_ene
        rixs_drv.gamma = .16 / hartree_in_ev()

        if nvir is not None:
            rixs_res = rixs_drv.compute(mol, bas, scf_results, valence_rsp)
        elif cvs_scf:
            rixs_res = rixs_drv.compute(mol, bas, scf_results, valence_rsp, cvs_rsp_res, 
                                        cvs_scf_tensors=cvs_scf_res)
        else:
            rixs_res = rixs_drv.compute(mol, bas, scf_results, valence_rsp, cvs_rsp_res)

        if scf_drv.rank == mpi_master():
            assert np.allclose(ref_xsection, rixs_res['cross_sections'][:,0], 
                               rtol=1e-5, atol=1e-8)

    def test_hf_svp_rpa_rsa(self):

        ref_xsection = np.array([
            17793.089879968, 1230.117099682, 48.724646624, 733.623643582,
            35.287774695, 40.832683493, 0.106152313, 103.187106313, 0.247740469,
            2.045706565, 0.047730756, 0.033434253, 3.891997487, 75.929243099,
            3.905579709, 4.827752013
        ])

        self.run_rixs('hf', 'def2-svp', ref_xsection, 1, 20, ncorestates=2, nvir=16, nvalence=1, cvs_scf=True, tda=False)

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
            5555.404035122, 29.495579839, 0.486437651, 4.647466241, 0.358806671,
            0.004783713, 0.050625028, 0.408792947, 0.056138218, 0.117741899,
            0.135671787, 0.000000121, 0.00033592, 3.671544966, 1.92188494,
            2.48548873
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
            18244.015365242, 1293.685311641, 6818.799975572, 1027.568339993,
            11762.311658725, 374.78336501, 50.417074025, 783.298480117,
            65.952379858, 241.303347781, 26.121923673, 58.828983154,
            166.19060267, 603.121079685, 74.825211603, 0.708703218,
            51.227812446, 43.190650184, 58.917940084, 83.487470476
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

        self.run_rixs('hf', 'def2-svp', ref_xsection, 1, 20, ncorestates=4, cvs_scf=True, tda=True)

    def test_b3lyp_svp_rpa_2s(self):

        ref_xsection = np.array([
            5419.18617506, 29.653929176, 2022.630087752, 50.15631254,
            3601.586702326, 15.477390187, 0.513680457, 4.382418983, 1.939015588,
            3.367438139, 134.703296992, 0.538659948, 2.631649603, 0.008582586,
            7.465709298, 4.980095871, 0.690127094, 0.003739568, 0.427996885,
            0.003540102
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
            22977.685596678, 868.830525985, 8780.263141343, 920.421671971,
            14842.733521127, 277.510132332, 6.064776563, 846.587889712,
            106.15909942
        ])

        self.run_rixs('hf', '6-31G*', ref_xsection, 1, 57, nvir=13, nvalence=4, tda=False, cutoff_ene=31.5 / hartree_in_ev())
