from mpi4py import MPI
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

    def run_scf(self):

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
        scf_drv.xcfun = xcfun_label
        scf_drv.ostream.mute()
        scf_result = scf_drv.compute(molecule, basis)

        return molecule, basis, xcfun_label, scf_result

    def test_rsa_core_excitation_rpa(self):

        molecule, basis, xcfun_label, scf_result = self.run_scf()

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

        ref_exc_ene = np.array([
            19.09409611, 19.15795969, 19.20965780, 19.21608706, 19.25755683,
            19.28712810, 19.33009437, 19.33571861, 19.36632321, 19.44738600
        ])
        ref_osc_str = np.array([
            0.0100, 0.0221, 0.0118, 0.0084, 0.0047, 0.0028, 0.0000, 0.0086,
            0.0044, 0.0017
        ])

        if rpa_drv.rank == mpi_master():
            exc_ene = rpa_result['eigenvalues']
            osc_str = rpa_result['oscillator_strengths']
            assert np.max(np.abs(exc_ene - ref_exc_ene)) < 1.0e-4
            assert np.max(np.abs(osc_str - ref_osc_str)) < 1.0e-4

    def test_rsa_core_excitation_tda(self):

        molecule, basis, xcfun_label, scf_result = self.run_scf()

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

        ref_exc_ene = np.array([
            19.09414620, 19.15797993, 19.20966688, 19.21610156, 19.25756094,
            19.28713597, 19.33009443, 19.33573228, 19.36632620, 19.44744147
        ])
        ref_osc_str = np.array([
            0.0101, 0.0222, 0.0119, 0.0084, 0.0048, 0.0028, 0.0000, 0.0087,
            0.0044, 0.0017
        ])

        if tda_drv.rank == mpi_master():
            exc_ene = tda_result['eigenvalues']
            osc_str = tda_result['oscillator_strengths']
            assert np.max(np.abs(exc_ene - ref_exc_ene)) < 1.0e-4
            assert np.max(np.abs(osc_str - ref_osc_str)) < 1.0e-4

    def run_rsa(self, xcfun_label, basis_label, ref_exc_enes, ref_osc_str, tol,
                ncore=None, nvir=None, nval=None, tda=False):

        xyz_string = """3
        xyz
        O   -0.1858140  -1.1749469   0.7662596
        H   -0.1285513  -0.8984365   1.6808606
        H   -0.0582782  -0.3702550   0.2638279
        """
        mol = Molecule.read_xyz_string(xyz_string)
        bas = MolecularBasis.read(mol, basis_label, ostream=None)

        scf_drv = ScfRestrictedDriver()
        scf_drv.ostream.mute()
        scf_drv.xcfun = xcfun_label

        scf_results = scf_drv.compute(mol, bas)

        nocc = mol.number_of_alpha_electrons()

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

    def test_hf_svp_rpa(self):
        # verify that rsa reproduces full valence excitations
        ref_exc_enes = np.array(
            [0.33973039, 0.40464346, 0.43325535, 0.49805052, 0.55325390])

        ref_osc_str = np.array(
            [0.023665, 0.000000, 0.097765, 0.086454, 0.291919])

        self.run_rsa('hf', 'def2-svp', ref_exc_enes, ref_osc_str, 1.0e-6)

    def test_camb3lyp_svp_rpa(self):
        # verify that rsa reproduces full valence excitations

        ref_exc_enes = np.array(
            [0.28249530, 0.35533592, 0.36664012, 0.44312170, 0.51548914])

        ref_osc_str = np.array(
            [0.018219, 0.000000, 0.077391, 0.059649, 0.273690])

        self.run_rsa('cam-b3lyp', 'def2-svp', ref_exc_enes, ref_osc_str, 1.0e-5)

    def test_hf_svp_tda(self):
        # verify that rsa reproduces full valence excitations

        ref_exc_enes = np.array(
            [0.34189725, 0.40717708, 0.43577895, 0.50150019, 0.55485041])

        ref_osc_str = np.array([0.0229, 0.0000, 0.1039, 0.0975, 0.3060])

        self.run_rsa('hf', 'def2-svp', ref_exc_enes, ref_osc_str, 1.0e-6, tda=True)

    def test_b3lyp_svp_tda(self):
        # verify that rsa reproduces full valence excitations

        ref_exc_enes = np.array(
            [0.28045483, 0.35053288, 0.36504833, 0.43904962, 0.51570229])

        ref_osc_str = np.array([0.0180, 0.0000, 0.0854, 0.0695, 0.3006])

        self.run_rsa('b3lyp', 'def2-svp', ref_exc_enes, ref_osc_str, 1.0e-5, tda=True)
    
    def test_b3lyp_svp_rsa_rpa(self):

        ref_exc_enes = np.array(
            [0.28013457, 0.35056008, 0.3692739 , 0.44042579, 0.77320255])

        ref_osc_str = np.array([0.0175, 0.    , 0.1071, 0.0965, 0.    ])

        self.run_rsa('b3lyp', 'def2-svp', ref_exc_enes, ref_osc_str, 1.0e-5,
                     ncore=1, nvir=15, nval=2)

    def test_b3lyp_svp_rsa_tda(self):

        ref_exc_enes = np.array(
            [0.28099976, 0.35063639, 0.37179368, 0.44191394, 0.77325388])

        ref_osc_str = np.array([0.0175, 0.    , 0.1166, 0.1079, 0.    ])

        self.run_rsa('b3lyp', 'def2-svp', ref_exc_enes, ref_osc_str, 1.0e-5,
                     ncore=1, nvir=15, nval=2, tda=True)
