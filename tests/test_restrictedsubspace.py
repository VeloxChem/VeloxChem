import numpy as np
import pytest

from veloxchem.veloxchemlib import mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.lreigensolver import LinearResponseEigenSolver
from veloxchem.tdaeigensolver import TdaEigenSolver


@pytest.mark.solvers
class TestRestrictedSubspace:

    def run_scf(self, ri):

        water_xyz = """
        O   0.0   0.0   0.0
        H   0.0   1.4   1.1
        H   0.0  -1.4   1.1
        """
        molecule = Molecule.read_molecule_string(water_xyz, units='au')
        basis = MolecularBasis.read(molecule, 'def2-svpd', ostream=None)
        xcfun_label = 'b3lyp'

        # ground state
        scf_drv = ScfRestrictedDriver()
        scf_drv.ri = ri
        scf_drv.xcfun = xcfun_label
        scf_drv.ostream.mute()
        scf_result = scf_drv.compute(molecule, basis)

        if scf_drv.rank == mpi_master():
            if ri:
                assert scf_result['ri_jk'] or scf_result['ri_coulomb']
                assert not scf_result['ri_fallback']

        return molecule, basis, xcfun_label, scf_result

    @pytest.mark.parametrize('ri', [False, True])
    def test_rsa_core_excitation_rpa(self, ri):

        molecule, basis, xcfun_label, scf_result = self.run_scf(ri)

        # verify rsa reproduces core excitation
        rpa_drv = LinearResponseEigenSolver()
        rpa_drv.xcfun = xcfun_label
        rpa_drv.restricted_subspace = True
        rpa_drv.num_core_orbitals = 1
        rpa_drv.num_valence_orbitals = 0
        rpa_drv.num_virtual_orbitals = 34
        rpa_drv.nstates = 10
        rpa_drv.ostream.mute()
        rpa_result = rpa_drv.compute(molecule, basis, scf_result)

        ref_exc_ene = {
            False: np.array([
                19.09409612, 19.15795970, 19.20965781, 19.21608707, 19.25755683,
                19.28712811, 19.33009438, 19.33571862, 19.36632322, 19.44542822
            ]),
            True: np.array([
                19.09410513, 19.15800754, 19.20970148, 19.21606118, 19.25756739,
                19.28711518, 19.33015103, 19.33577924, 19.36633514, 19.44546557
            ])
        }[ri]
        ref_osc_str = {
            False: np.array([
                0.01, 0.0221, 0.0118, 0.0084, 0.0047, 0.0028, 0., 0.0086,
                0.0044, 0.0185
            ]),
            True: np.array([
                0.010018, 0.022154, 0.011829, 0.008422, 0.004724, 0.002782,
                0.000000, 0.008625, 0.004396, 0.018485
            ])
        }[ri]

        if rpa_drv.rank == mpi_master():
            exc_ene = rpa_result['eigenvalues']
            osc_str = rpa_result['oscillator_strengths']
            assert np.max(np.abs(exc_ene - ref_exc_ene)) < 1.0e-4
            assert np.max(np.abs(osc_str - ref_osc_str)) < 1.0e-4

    @pytest.mark.parametrize('ri', [False, True])
    def test_rsa_core_excitation_tda(self, ri):

        molecule, basis, xcfun_label, scf_result = self.run_scf(ri)

        # verify that rsa reproduces core excitation
        tda_drv = TdaEigenSolver()
        tda_drv.xcfun = xcfun_label
        tda_drv.restricted_subspace = True
        tda_drv.num_core_orbitals = 1
        tda_drv.num_valence_orbitals = 0
        tda_drv.num_virtual_orbitals = 34
        tda_drv.nstates = 10
        tda_drv.ostream.mute()
        tda_result = tda_drv.compute(molecule, basis, scf_result)

        ref_exc_ene = {
            False: np.array([
                19.0941462, 19.15797992, 19.20966688, 19.21610156, 19.25756094,
                19.28713597, 19.33009442, 19.33573228, 19.36632619, 19.44544945
            ]),
            True: np.array([
                19.09415281, 19.15802601, 19.20970930, 19.21607423, 19.25757065,
                19.28712239, 19.33015081, 19.33579136, 19.36633746, 19.44548426
            ])
        }[ri]
        ref_osc_str = {
            False: np.array([
                0.0101, 0.0222, 0.0119, 0.0084, 0.0048, 0.0028, 0., 0.0087,
                0.0044, 0.0186
            ]),
            True: np.array([
                0.010061, 0.022234, 0.011863, 0.008441, 0.004742, 0.002792,
                0.000000, 0.008662, 0.004409, 0.018564
            ])
        }[ri]

        if tda_drv.rank == mpi_master():
            exc_ene = tda_result['eigenvalues']
            osc_str = tda_result['oscillator_strengths']
            assert np.max(np.abs(exc_ene - ref_exc_ene)) < 1.0e-4
            assert np.max(np.abs(osc_str - ref_osc_str)) < 1.0e-4

    def run_rsa(self,
                xcfun_label,
                basis_label,
                ri,
                ref_exc_enes,
                ref_osc_str,
                tol,
                ncore=None,
                nvir=None,
                nval=None,
                tda=False):

        xyz_string = """3
        xyz
        O   -0.1858140  -1.1749469   0.7662596
        H   -0.1285513  -0.8984365   1.6808606
        H   -0.0582782  -0.3702550   0.2638279
        """
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

        nocc = mol.number_of_alpha_occupied_orbitals(bas)

        if scf_drv.rank == mpi_master():
            norb = scf_results['C_alpha'].shape[0]
        else:
            norb = None
        norb = scf_drv.comm.bcast(norb, root=mpi_master())

        if nvir is None:
            ncore = 1
            nval = nocc - 1
            nvir = norb - nocc

        if tda:
            lr_drv = TdaEigenSolver()
        else:
            lr_drv = LinearResponseEigenSolver()

        lr_drv.ostream.mute()
        lr_drv.restricted_subspace = True
        lr_drv.num_core_orbitals = ncore
        lr_drv.num_valence_orbitals = nval
        lr_drv.num_virtual_orbitals = nvir
        lr_drv.nstates = 5
        lr_results = lr_drv.compute(mol, bas, scf_results)

        if lr_drv.rank == mpi_master():
            assert np.max(np.abs(ref_exc_enes -
                                 lr_results['eigenvalues'])) < tol
            assert np.max(
                np.abs(ref_osc_str -
                       lr_results['oscillator_strengths'])) < 1.0e-4

    @pytest.mark.parametrize('ri', [False, True])
    def test_hf_svp_rpa(self, ri):
        # verify that rsa reproduces full valence excitations
        ref_exc_enes = {
            False: np.array(
                [0.33973039, 0.40464346, 0.43325535, 0.49805052, 0.55325390]),
            True: np.array(
                [0.33973201, 0.40463916, 0.43325282, 0.49804458, 0.55326636])
        }[ri]

        ref_osc_str = np.array(
            [0.023665, 0.000000, 0.097765, 0.086454, 0.291919])

        self.run_rsa('hf', 'def2-svp', ri, ref_exc_enes, ref_osc_str, 1.0e-6)

    @pytest.mark.parametrize('ri', [False, True])
    def test_camb3lyp_svp_rpa(self, ri):
        # verify that rsa reproduces full valence excitations

        ref_exc_enes = {
            False: np.array(
                [0.28249530, 0.35533592, 0.36664012, 0.44312170, 0.51548914]),
            True: np.array(
                [0.28249542, 0.35532788, 0.36664161, 0.44311429, 0.51549785])
        }[ri]

        ref_osc_str = np.array(
            [0.018219, 0.000000, 0.077391, 0.059649, 0.273690])

        self.run_rsa('cam-b3lyp', 'def2-svp', ri, ref_exc_enes, ref_osc_str,
                     1.0e-5)

    @pytest.mark.parametrize('ri', [False, True])
    def test_hf_svp_tda(self, ri):
        # verify that rsa reproduces full valence excitations

        ref_exc_enes = {
            False: np.array(
                [0.34189725, 0.40717708, 0.43577895, 0.50150019, 0.55485041]),
            True: np.array(
                [0.34189759, 0.40716615, 0.43576509, 0.50149749, 0.55485206])
        }[ri]

        ref_osc_str = np.array([0.0229, 0.0000, 0.1039, 0.0975, 0.3060])

        self.run_rsa('hf',
                     'def2-svp',
                     ri,
                     ref_exc_enes,
                     ref_osc_str,
                     1.0e-6,
                     tda=True)

    @pytest.mark.parametrize('ri', [False, True])
    def test_b3lyp_svp_tda(self, ri):
        # verify that rsa reproduces full valence excitations

        ref_exc_enes = {
            False: np.array(
                [0.28045483, 0.35053288, 0.36504833, 0.43904962, 0.51570229]),
            True: np.array(
                [0.28045263, 0.35052101, 0.36504338, 0.43903883, 0.51570051])
        }[ri]

        ref_osc_str = {
            False: np.array([0.0180, 0.0000, 0.0854, 0.0695, 0.3006]),
            True: np.array([0.018025, 0.000000, 0.085372, 0.069569, 0.300492])
        }[ri]

        self.run_rsa('b3lyp',
                     'def2-svp',
                     ri,
                     ref_exc_enes,
                     ref_osc_str,
                     1.0e-5,
                     tda=True)

    @pytest.mark.parametrize('ri', [False, True])
    def test_b3lyp_svp_rsa_rpa(self, ri):

        ref_exc_enes = {
            False: np.array(
                [0.28013457, 0.35056008, 0.3692739, 0.44042579, 0.77320255]),
            True: np.array(
                [0.28013454, 0.35055189, 0.36927420, 0.44041705, 0.77319473])
        }[ri]

        ref_osc_str = {
            False: np.array([0.0175, 0., 0.1071, 0.0965, 0.]),
            True: np.array([0.017475, 0.000000, 0.107070, 0.096444, 0.000000])
        }[ri]

        self.run_rsa('b3lyp',
                     'def2-svp',
                     ri,
                     ref_exc_enes,
                     ref_osc_str,
                     1.0e-5,
                     ncore=1,
                     nvir=15,
                     nval=2)

    @pytest.mark.parametrize('ri', [False, True])
    def test_b3lyp_svp_rsa_tda(self, ri):

        ref_exc_enes = {
            False: np.array(
                [0.28099976, 0.35063639, 0.37179368, 0.44191394, 0.77325388]),
            True: np.array(
                [0.28099743, 0.35062442, 0.37178645, 0.44190153, 0.77320496])
        }[ri]

        ref_osc_str = {
            False: np.array([0.0175, 0., 0.1166, 0.1079, 0.]),
            True: np.array([0.017482, 0.000000, 0.116598, 0.107836, 0.000000])
        }[ri]

        self.run_rsa('b3lyp',
                     'def2-svp',
                     ri,
                     ref_exc_enes,
                     ref_osc_str,
                     1.0e-5,
                     ncore=1,
                     nvir=15,
                     nval=2,
                     tda=True)
