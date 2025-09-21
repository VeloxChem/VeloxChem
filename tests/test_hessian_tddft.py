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
from veloxchem.tddftgradientdriver import TddftGradientDriver
from veloxchem.tddfthessiandriver import TddftHessianDriver
from veloxchem.vibrationalanalysis import VibrationalAnalysis

class TestTddftHessianDriver:

    def run_tddft_hessian(self, xcfun_label, tamm_dancoff, ref_label):

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
        h5file = str(here / 'data' / 'tddft_hessian.h5')

        scf_drv = ScfRestrictedDriver()
        scf_drv.xcfun = xcfun_label
        scf_drv.grid_level = 6
        scf_drv.conv_thresh = 1e-10
        scf_results = scf_drv.compute(molecule, basis)

        if tamm_dancoff:
            rsp_drv = TdaEigenSolver()
        else:
            rsp_drv = LinearResponseEigenSolver()
        rsp_drv.conv_thresh = 1e-8
        rsp_drv.nstates = 5
        rsp_results = rsp_drv.compute(molecule, basis, scf_results)

        grad_drv = TddftGradientDriver(scf_drv)
        grad_dict = {'state_deriv_index': [1], 'do_first_order_prop': 'yes'}
        cphf_dict = {'conv_thresh': 1e-8, 'use_subspace_solver': 'yes'}
        grad_drv.update_settings(grad_dict=grad_dict,
                                 rsp_dict={'conv_thresh':1e-8},
                                 orbrsp_dict=cphf_dict)

        grad_drv.compute(molecule, basis, scf_drv, rsp_drv, rsp_results)

        hessian_drv = TddftHessianDriver(scf_drv, rsp_drv=rsp_drv,
                                             tddft_grad_drv=grad_drv)
        hessian_drv.do_dipole_gradient = True
        hessian_drv.compute(molecule, basis)

        if MPI.COMM_WORLD.Get_rank() == mpi_master():

            hf = h5py.File(h5file)
            ref_hessian = np.array(hf.get('hessian_'+ref_label))
            ref_dipole_gradient = np.array(hf.get('dipole_gradient_'+ref_label)) 
            hf.close()

            diff_hessian = np.max(np.abs(hessian_drv.hessian - ref_hessian))
            diff_dipole_grad = np.max(np.abs(hessian_drv.dipole_gradient -
                                             ref_dipole_gradient))

            assert diff_hessian < 1.0e-5
            assert diff_dipole_grad < 1.0e-5

    def test_tda(self):
        xcfun_label = "hf"
        tamm_dancoff = True
        ref_label = "tda_hf"
        self.run_tddft_hessian(xcfun_label, tamm_dancoff, ref_label)

    @pytest.mark.timeconsuming
    def test_tda_pbe(self):
        xcfun_label = "pbe"
        tamm_dancoff = True
        ref_label = "tda_pbe"
        self.run_tddft_hessian(xcfun_label, tamm_dancoff, ref_label)

    @pytest.mark.timeconsuming
    def test_tda_pbe0(self):
        xcfun_label = "pbe0"
        tamm_dancoff = True
        ref_label = "tda_pbe0"
        self.run_tddft_hessian(xcfun_label, tamm_dancoff, ref_label)

    def test_rpa(self):
        xcfun_label = "hf"
        tamm_dancoff = False
        ref_label = "rpa_hf"
        self.run_tddft_hessian(xcfun_label, tamm_dancoff, ref_label)

    @pytest.mark.timeconsuming
    def test_tddft_pbe(self):
        xcfun_label = "pbe"
        tamm_dancoff = False
        ref_label = "rpa_pbe"
        self.run_tddft_hessian(xcfun_label, tamm_dancoff, ref_label)

    @pytest.mark.timeconsuming
    def test_tddft_pbe0(self):
        xcfun_label = "pbe0"
        tamm_dancoff = False
        ref_label = "rpa_pbe0"
        self.run_tddft_hessian(xcfun_label, tamm_dancoff, ref_label)
