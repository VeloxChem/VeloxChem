import numpy as np
import sys
import pytest

from veloxchem.veloxchemlib import is_mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.tdaeigensolver import TdaEigenSolver
from veloxchem.lreigensolver import LinearResponseEigenSolver
from veloxchem.gradientdriver import GradientDriver
from veloxchem.tddftgradientdriver import TddftGradientDriver

try:
    import pyscf
except ImportError:
    pass

class TestGrad:
    def run_tddft_grad(self, xcfun_label, tamm_dancoff, ref_grad):

        molecule_string = """
        O 0.0000000000    0.0000000000   -0.0254395383
        H 0.0000000000    0.7695699584    0.5948147012
        H 0.0000000000   -0.7695699584    0.5948147012
        """

        basis_set_label = "def2-svp"

        molecule = Molecule.read_str(molecule_string, units='angstrom')
        basis = MolecularBasis.read(molecule, basis_set_label, ostream=None)

        scf_drv = ScfRestrictedDriver()

        scf_dict = {'conv_thresh': 1e-8}
        method_dict = {'xcfun': xcfun_label, 'grid_level': 5}

        scf_drv.update_settings(scf_dict, method_dict)

        scf_drv.compute(molecule, basis)

        if tamm_dancoff:
            rsp_solver = TdaEigenSolver()
        else:
            rsp_solver = LinearResponseEigenSolver()

        if tamm_dancoff:
            tda = 'yes'
        else:
            tda = 'no'
        rsp_dict = {"tamm_dancoff": tda}

        rsp_solver.update_settings(rsp_dict, method_dict)
        rsp_solver.conv_thresh = 1e-5
        rsp_solver.nstates = 3

        rsp_results = rsp_solver.compute(molecule, basis, scf_drv.scf_tensors)

        tddft_grad = TddftGradientDriver(scf_drv)

        orbrsp_dict = {"conv_thresh": 1e-7}
        tddft_grad.update_settings({}, rsp_dict, orbrsp_dict, method_dict)
        tddft_grad.compute(molecule, basis, rsp_solver, rsp_results)

        if is_mpi_master():
            grad = tddft_grad.get_gradient()[0]
            assert np.max(np.abs(grad - ref_grad)) < 1.0e-6

    @pytest.mark.skipif('pyscf' not in sys.modules,
                        reason='pyscf for integral derivatives not available')
    def test_tda_slater(self):
        xcfun_label = 'slater'
        tamm_dancoff = True
        ref_grad = np.array(
            [[-0.000000000000001, -0.000000000000005,  0.096828839746124],
              [ 0.              ,  -0.059512179812975, -0.048405146528462],
              [-0.              ,   0.059512179812977, -0.048405146528465]])

        self.run_tddft_grad(xcfun_label, tamm_dancoff, ref_grad)

    @pytest.mark.skipif('pyscf' not in sys.modules,
                        reason='pyscf for integral derivatives not available')
    def test_tda_blyp(self):
        xcfun_label = 'blyp'
        tamm_dancoff = True
        ref_grad = np.array(
            [[ 0.000000000000002,  0.000000000000003,  0.081702810507772],
              [-0.               , -0.052036797551521, -0.040837361620094],
              [ 0.               ,  0.052036797551518, -0.040837361620095]])

        self.run_tddft_grad(xcfun_label, tamm_dancoff, ref_grad)

    @pytest.mark.skipif('pyscf' not in sys.modules,
                        reason='pyscf for integral derivatives not available')
    def test_rpa_slater(self):
        xcfun_label = 'slater'
        tamm_dancoff = False
        ref_grad = np.array(
            [[ 0.000000000000001, -0.000000000000006,  0.097148642876961],
              [ 0.,                -0.060393461796934, -0.048565041680259],
              [ 0.,                 0.060393461796936, -0.048565041680257]])

        self.run_tddft_grad(xcfun_label, tamm_dancoff, ref_grad)

    @pytest.mark.skipif('pyscf' not in sys.modules,
                        reason='pyscf for integral derivatives not available')
    def test_rpa_blyp(self):
        xcfun_label = 'blyp'
        tamm_dancoff = False
        ref_grad = np.array(
            [[-0.000000000000001, -0.000000000000001,  0.0818110121047  ],
              [-0.               , -0.052800208718226, -0.040891498953149],
              [ 0.               ,  0.052800208718228, -0.040891498953151]])

        self.run_tddft_grad(xcfun_label, tamm_dancoff, ref_grad)
