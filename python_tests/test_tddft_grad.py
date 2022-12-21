import numpy as np

from veloxchem.veloxchemlib import is_mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.tdaeigensolver import TdaEigenSolver
from veloxchem.lreigensolver import LinearResponseEigenSolver
from veloxchem.gradientdriver import GradientDriver
from veloxchem.tdhfgradientdriver import TdhfGradientDriver

class TestGrad:
    pass
    # TODO: this does not work because it relies on pyscf integral derivatives
    # and pyscf is not available on monkey;
    # add back in when integral derivatives are available. 
   # def run_tddft_grad(self, xcfun_label, tamm_dancoff, ref_grad):

   #     molecule_string = """
   #     O 0.0000000000    0.0000000000   -0.0254395383
   #     H 0.0000000000    0.7695699584    0.5948147012
   #     H 0.0000000000   -0.7695699584    0.5948147012
   #     """

   #     basis_set_label = "def2-svp"

   #     molecule = Molecule.read_str(molecule_string, units='angstrom')
   #     basis = MolecularBasis.read(molecule, basis_set_label, ostream=None)

   #     scf_drv = ScfRestrictedDriver()

   #     scf_drv.dft = True
   #     scf_drv.xcfun = xcfun_label
   #     scf_drv.conv_thresh = 1e-8
   #     scf_drv.checkpoint_file = None

   #     scf_drv.compute(molecule, basis)

   #     if tamm_dancoff:
   #         rsp_solver = TdaEigenSolver()
   #     else:
   #         rsp_solver = LinearResponseEigenSolver()

   #     rsp_solver.dft = True
   #     rsp_solver.tamm_dancoff = tamm_dancoff
   #     rsp_solver.xcfun = scf_drv.xcfun
   #     rsp_solver.conv_thresh = 1e-5
   #     rsp_solver.nstates = 3
   #     rsp_solver.checkpoint_file = None

   #     rsp_results = rsp_solver.compute(molecule, basis, scf_drv.scf_tensors)

   #     tddft_grad = TdhfGradientDriver(scf_drv)
   #     if tamm_dancoff:
   #         tda = 'yes'
   #     else:
   #         tda = 'no'
   #     rsp_dict = {"tamm_dancoff": tda}
   #     method_dict = {'xcfun': xcfun_label}
   #     orbrsp_dict = {"conv_thresh": 1e-7}
   #     tddft_grad.update_settings({}, rsp_dict, orbrsp_dict, method_dict)
   #     tddft_grad.dft = True
   #     tddft_grad.xcfun = scf_drv.xcfun
   #     tddft_grad.compute(molecule, basis, rsp_solver, rsp_results)

   #     if is_mpi_master():
   #         grad = tddft_grad.get_gradient()
   #         print("Gradient:\n", grad)
   #         print("\nReference gradient:\n", ref_grad) 
   #         assert np.max(np.abs(grad - ref_grad)) < 1.0e-6


   # def test_tda_slater(self):
   #     xcfun_label = 'slater'
   #     tamm_dancoff = True
   #     ref_grad = np.array(
   #         [[-0. ,  0.000000000000008,  0.096834496224282],
   #          [-0. , -0.059511770311381, -0.048404984179591],
   #          [ 0. ,  0.059511770311380, -0.048404984179587]])

   #     self.run_tddft_grad(xcfun_label, tamm_dancoff, ref_grad)

   # def test_tda_blyp(self):
   #     xcfun_label = 'blyp'
   #     tamm_dancoff = True
   #     ref_grad = np.array(
   #         [[-0. ,  0.000000000000000,  0.081712506517363],
   #          [-0. , -0.052036762539477, -0.040837608701902],
   #          [ 0. ,  0.052036762539478, -0.040837608701899]])

   #     self.run_tddft_grad(xcfun_label, tamm_dancoff, ref_grad)


   # def test_rpa_slater(self):
   #     xcfun_label = 'slater'
   #     tamm_dancoff = False
   #     ref_grad = np.array(
   #         [[ 0. ,  0.000000000000001,  0.097154237192811],
   #          [-0. , -0.060393010060886, -0.048564843803965],
   #          [ 0. ,  0.060393010060885, -0.048564843803961]])

   #     self.run_tddft_grad(xcfun_label, tamm_dancoff, ref_grad)

   # def test_rpa_blyp(self):
   #     xcfun_label = 'blyp'
   #     tamm_dancoff = False
   #     ref_grad = np.array(
   #         [[ 0. ,  0.000000000000008,  0.081819790394095],
   #          [-0. , -0.052800069387996, -0.040891688938228],
   #          [ 0. ,  0.052800069387988, -0.040891688938229]])

   #     self.run_tddft_grad(xcfun_label, tamm_dancoff, ref_grad)
