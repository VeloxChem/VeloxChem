from mpi4py import MPI
from pathlib import Path
import numpy as np
import h5py
import pytest

from veloxchem.veloxchemlib import mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.tdaeigensolver import TdaEigenSolver
from veloxchem.lreigensolver import LinearResponseEigenSolver
from veloxchem.tddfthessiandriver import TddftHessianDriver


class TestTddftHessianDriver:

    def run_tddft_hessian(self, xcfun_label, tamm_dancoff, ref_label, *, ri):

        molecule_xyz_string = """
            3

            O       0.0000000000     0.0000000000     0.1240508479
            H      -0.7695699584     0.0000000000    -0.4962033916
            H       0.7695699584    -0.0000000000    -0.4962033916
            """
        basis_set_label = "6-31G"
        molecule = Molecule.from_xyz_string(molecule_xyz_string)
        basis = MolecularBasis.read(molecule, basis_set_label)

        here = Path(__file__).parent
        h5file = str(here / 'data' / {
            False: 'tddft_hessian.h5',
            True: 'tddft_hessian_ri.h5'
        }[ri])

        scf_drv = ScfRestrictedDriver()
        scf_drv.ri = ri
        scf_drv.xcfun = xcfun_label
        scf_drv.grid_level = 6
        scf_drv.conv_thresh = 1e-10
        scf_drv.ostream.mute()
        scf_results = scf_drv.compute(molecule, basis)

        if scf_drv.rank == mpi_master():
            if ri:
                assert scf_results['ri_jk'] or scf_results['ri_coulomb']
                assert not scf_results['ri_fallback']

        if tamm_dancoff:
            rsp_drv = TdaEigenSolver()
        else:
            rsp_drv = LinearResponseEigenSolver()
        rsp_drv.conv_thresh = 1e-8
        rsp_drv.nstates = 5
        rsp_drv.ostream.mute()
        rsp_results_not_used = rsp_drv.compute(molecule, basis, scf_results)

        grad_dict = {}
        cphf_dict = {'conv_thresh': 1e-8, 'use_subspace_solver': 'yes'}

        hessian_drv = TddftHessianDriver(scf_drv,
                                         rsp_drv=rsp_drv)
        hessian_drv.update_settings(method_dict={}, grad_dict=grad_dict,
                                    cphf_dict=cphf_dict)
        hessian_drv.do_dipole_gradient = True
        hessian_drv.ostream.mute()
        hessian_drv.compute(molecule, basis)

        if MPI.COMM_WORLD.Get_rank() == mpi_master():

            hf = h5py.File(h5file)
            ref_hessian = np.array(hf.get('hessian_' + ref_label))
            ref_dipole_gradient = np.array(
                hf.get('dipole_gradient_' + ref_label))
            hf.close()

            diff_hessian = np.max(np.abs(hessian_drv.hessian - ref_hessian))
            diff_dipole_grad = np.max(
                np.abs(hessian_drv.dipole_gradient - ref_dipole_gradient))

            assert diff_hessian < 1.0e-5
            assert diff_dipole_grad < 1.0e-5

    @pytest.mark.solvers
    @pytest.mark.parametrize('ri', [False, True])
    def test_tda(self, ri):
        xcfun_label = "hf"
        tamm_dancoff = True
        ref_label = "tda_hf"
        self.run_tddft_hessian(xcfun_label, tamm_dancoff, ref_label, ri=ri)

    @pytest.mark.timeconsuming
    @pytest.mark.parametrize('ri', [False, True])
    def test_tda_pbe(self, ri):
        xcfun_label = "pbe"
        tamm_dancoff = True
        ref_label = "tda_pbe"
        self.run_tddft_hessian(xcfun_label, tamm_dancoff, ref_label, ri=ri)

    @pytest.mark.timeconsuming
    @pytest.mark.parametrize('ri', [False, True])
    def test_tda_pbe0(self, ri):
        xcfun_label = "pbe0"
        tamm_dancoff = True
        ref_label = "tda_pbe0"
        self.run_tddft_hessian(xcfun_label, tamm_dancoff, ref_label, ri=ri)

    @pytest.mark.timeconsuming
    @pytest.mark.parametrize('ri', [False, True])
    def test_tda_camb3lyp(self, ri):
        xcfun_label = "cam-b3lyp"
        tamm_dancoff = True
        ref_label = "tda_camb3lyp"
        self.run_tddft_hessian(xcfun_label, tamm_dancoff, ref_label, ri=ri)

    @pytest.mark.solvers
    @pytest.mark.parametrize('ri', [False, True])
    def test_rpa(self, ri):
        xcfun_label = "hf"
        tamm_dancoff = False
        ref_label = "rpa_hf"
        self.run_tddft_hessian(xcfun_label, tamm_dancoff, ref_label, ri=ri)

    @pytest.mark.timeconsuming
    @pytest.mark.parametrize('ri', [False, True])
    def test_tddft_pbe(self, ri):
        xcfun_label = "pbe"
        tamm_dancoff = False
        ref_label = "rpa_pbe"
        self.run_tddft_hessian(xcfun_label, tamm_dancoff, ref_label, ri=ri)

    @pytest.mark.timeconsuming
    @pytest.mark.parametrize('ri', [False, True])
    def test_tddft_pbe0(self, ri):
        xcfun_label = "pbe0"
        tamm_dancoff = False
        ref_label = "rpa_pbe0"
        self.run_tddft_hessian(xcfun_label, tamm_dancoff, ref_label, ri=ri)

    @pytest.mark.timeconsuming
    @pytest.mark.parametrize('ri', [False, True])
    def test_tddft_camb3lyp(self, ri):
        xcfun_label = "cam-b3lyp"
        tamm_dancoff = False
        ref_label = "rpa_camb3lyp"
        self.run_tddft_hessian(xcfun_label, tamm_dancoff, ref_label, ri=ri)

    @pytest.mark.solvers
    @pytest.mark.parametrize('ri', [False, True])
    def test_dipole_gradient_of_charged_molecule(self, ri):

        molecule = Molecule.read_xyz_string("""2
        hydroxide
        O    0.000000    0.000000    0.000000
        H    0.000000    0.000000    0.970000
        """)
        molecule.set_charge(-1)
        basis = MolecularBasis.read(molecule, 'sto-3g', ostream=None)

        scf_drv = ScfRestrictedDriver()
        scf_drv.ri = ri
        scf_drv.ostream.mute()
        scf_results = scf_drv.compute(molecule, basis)

        if scf_drv.rank == mpi_master():
            if ri:
                assert scf_results['ri_jk'] or scf_results['ri_coulomb']
                assert not scf_results['ri_fallback']

        rsp_drv = TdaEigenSolver()
        rsp_drv.nstates = 1
        rsp_drv.ostream.mute()
        rsp_drv.compute(molecule, basis, scf_results)

        hessian_drv = TddftHessianDriver(scf_drv, rsp_drv=rsp_drv)
        hessian_drv.do_dipole_gradient = True
        hessian_drv.ostream.mute()
        hessian_drv.compute(molecule, basis)

        # Moving the whole molecule by t changes the dipole moment of a
        # molecule with charge q by q t when the origin is fixed, so the
        # dipole gradient summed over the atoms is q times the unit matrix.
        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            dipole_gradient = hessian_drv.dipole_gradient.reshape(3, -1, 3)
            assert np.max(np.abs(dipole_gradient.sum(axis=1) +
                                 np.eye(3))) < 1.0e-8
