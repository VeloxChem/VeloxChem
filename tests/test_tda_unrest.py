import numpy as np
import pytest

from veloxchem.veloxchemlib import mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfunrestdriver import ScfUnrestrictedDriver
from veloxchem.tdaeigensolverunrest import TdaUnrestrictedEigenSolver


@pytest.mark.solvers
class TestUnrestrictedTDA:

    def run_tda(self,
                xcfun_label,
                basis_label,
                ref_exc_enes,
                ref_osc_str,
                tol,
                max_subspace_dim=None):

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

        lr_drv = TdaUnrestrictedEigenSolver()
        lr_drv.ostream.mute()
        lr_drv.nstates = 10
        lr_drv.max_subspace_dim = max_subspace_dim
        lr_results = lr_drv.compute(mol, bas, scf_results)

        if lr_drv.rank == mpi_master():
            assert np.max(np.abs(ref_exc_enes -
                                 lr_results['eigenvalues'])) < tol
            assert np.max(
                np.abs(ref_osc_str -
                       lr_results['oscillator_strengths'])) < 1.0e-4

    def test_camb3lyp_svp(self):

        # vlxtag: UKS, Absorption, TDA

        ref_exc_enes = np.array([
            0.09065604, 0.09500459, 0.21881146, 0.47239371, 0.48842790,
            0.51728556, 0.52413401, 0.55026290, 0.55530527, 0.59463209
        ])

        ref_osc_str = np.array([
            0.0013, 0.2065, 0.0000, 0.0325, 0.1953, 0.0015, 0.0000, 0.0256,
            0.0701, 0.1727
        ])

        self.run_tda('cam-b3lyp',
                     'def2-svp',
                     ref_exc_enes,
                     ref_osc_str,
                     1.0e-5,
                     max_subspace_dim=80)

    def run_tda_with_ecp(self, ref_exc_enes, ref_osc_str, tol):

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

        lr_drv = TdaUnrestrictedEigenSolver()
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

        # vlxtag: UHF, Absorption, TDA

        ref_exc_enes = np.array(
            [0.08351815, 0.08351815, 0.10137196, 0.10137196, 0.15068886])
        ref_osc_str = np.array([0.0000, 0.0000, 0.0002, 0.0002, 0.0089])

        self.run_tda_with_ecp(ref_exc_enes, ref_osc_str, 1.0e-6)

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

        lr_drv = TdaUnrestrictedEigenSolver()
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
