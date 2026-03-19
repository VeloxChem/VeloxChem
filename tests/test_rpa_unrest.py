import numpy as np
import pytest

from veloxchem.veloxchemlib import mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfunrestdriver import ScfUnrestrictedDriver
from veloxchem.lreigensolverunrest import LinearResponseUnrestrictedEigenSolver


@pytest.mark.solvers
class TestUnrestrictedRPA:

    def run_rpa(self, xcfun_label, basis_label, ref_exc_enes, ref_osc_str, tol):

        xyz_string = """3
        xyz
        O      0.000000   0.000000   0.117790
        H      0.000000   0.755453  -0.471161
        H      0.000000  -0.755453  -0.471161
        """
        mol = Molecule.read_xyz_string(xyz_string)
        mol.set_multiplicity(3)

        bas = MolecularBasis.read(mol, basis_label, ostream=None)

        scf_drv = ScfUnrestrictedDriver()
        scf_drv.ostream.mute()
        scf_drv.xcfun = xcfun_label
        scf_results = scf_drv.compute(mol, bas)

        lr_drv = LinearResponseUnrestrictedEigenSolver()
        lr_drv.ostream.mute()
        lr_drv.nstates = 10
        lr_results = lr_drv.compute(mol, bas, scf_results)

        if lr_drv.rank == mpi_master():
            assert np.max(np.abs(ref_exc_enes -
                                 lr_results['eigenvalues'])) < tol
            assert np.max(
                np.abs(ref_osc_str -
                       lr_results['oscillator_strengths'])) < 1.0e-4

    def test_hf_svp(self):

        # vlxtag: UHF, Absorption, TDHF

        ref_exc_enes = np.array([
            0.07000019, 0.09444351, 0.23795183, 0.53171161, 0.54518734,
            0.55496039, 0.57527661, 0.58052563, 0.61929121, 0.63384529
        ])

        ref_osc_str = np.array([
            0.1803, 0.0013, 0.0000, 0.2184, 0.0358, 0.0002, 0.0187, 0.0000,
            0.0789, 0.0253
        ])

        self.run_rpa('hf', 'def2-svp', ref_exc_enes, ref_osc_str, 1.0e-6)

    def test_pbe0_svp(self):

        # vlxtag: UKS, Absorption, TDDFT

        ref_exc_enes = np.array([
            0.08832630, 0.09472770, 0.22512802, 0.48348997, 0.48361504,
            0.51691424, 0.52828051, 0.54298284, 0.56229157, 0.60539714
        ])

        ref_osc_str = np.array([
            0.1596, 0.0012, 0.0000, 0.1985, 0.0338, 0.0072, 0.0000, 0.0308,
            0.0600, 0.1210
        ])

        self.run_rpa('pbe0', 'def2-svp', ref_exc_enes, ref_osc_str, 1.0e-5)

    def run_rpa_with_ecp(self, ref_exc_enes, ref_osc_str, tol):

        xyz_string = """2
        xyz
        Au 0 0 0
        H  0 0 1.55
        """
        mol = Molecule.read_xyz_string(xyz_string)
        mol.set_charge(1)
        mol.set_multiplicity(2)

        bas = MolecularBasis.read(mol, 'def2-svp', ostream=None)

        scf_drv = ScfUnrestrictedDriver()
        scf_drv.ostream.mute()
        scf_results = scf_drv.compute(mol, bas)

        lr_drv = LinearResponseUnrestrictedEigenSolver()
        lr_drv.ostream.mute()
        lr_drv.nstates = 5
        lr_results = lr_drv.compute(mol, bas, scf_results)

        if lr_drv.rank == mpi_master():
            assert np.max(np.abs(ref_exc_enes -
                                 lr_results['eigenvalues'])) < tol
            assert np.max(
                np.abs(ref_osc_str -
                       lr_results['oscillator_strengths'])) < 1.0e-4

    def test_hf_with_ecp(self):

        # vlxtag: UHF, Absorption, TDHF

        ref_exc_enes = np.array(
            [0.08148954, 0.08148954, 0.09916297, 0.09916297, 0.12983930])
        ref_osc_str = np.array([0.0000, 0.0000, 0.0001, 0.0001, 0.0119])

        self.run_rpa_with_ecp(ref_exc_enes, ref_osc_str, 1.0e-6)

    def test_checkpoint_and_restart(self, tmp_path):

        xyz_string = """3
        xyz
        O   -0.1858140  -1.1749469   0.7662596
        H   -0.1285513  -0.8984365   1.6808606
        H   -0.0582782  -0.3702550   0.2638279
        """
        mol = Molecule.read_xyz_string(xyz_string)
        mol.set_charge(1.0)
        mol.set_multiplicity(2)
        bas = MolecularBasis.read(mol, 'def2-svp', ostream=None)

        scf_drv = ScfUnrestrictedDriver()
        scf_drv.ostream.mute()

        filename = str(tmp_path / "water_restart")
        # To avoid inconsistency across MPI ranks
        filename = scf_drv.comm.bcast(filename, root=mpi_master())

        scf_drv.filename = filename
        scf_results = scf_drv.compute(mol, bas)

        lr_drv = LinearResponseUnrestrictedEigenSolver()
        lr_drv.filename = filename
        lr_drv.ostream.mute()

        lr_drv.nstates = 5
        lr_results_not_used = lr_drv.compute(mol, bas, scf_results)

        lr_drv.restart = True
        lr_drv.nstates = 10
        lr_results_first = lr_drv.compute(mol, bas, scf_results)

        lr_drv.restart = False
        lr_drv.nstates = 10
        lr_results_second = lr_drv.compute(mol, bas, scf_results)

        if scf_drv.rank == mpi_master():
            for key in [
                    'eigenvalues',
                    'oscillator_strengths',
                    'rotatory_strengths',
            ]:
                assert np.max(
                    np.abs(lr_results_first[key] -
                           lr_results_second[key])) < 1e-10
