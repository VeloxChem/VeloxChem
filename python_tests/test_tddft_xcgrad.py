import numpy as np

from veloxchem.veloxchemlib import denmat
from veloxchem.veloxchemlib import is_mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.aodensitymatrix import AODensityMatrix
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.tdaexcidriver import TDAExciDriver
from veloxchem.lreigensolver import LinearResponseEigenSolver
from veloxchem.tdaorbitalresponse import TdaOrbitalResponse
from veloxchem.rpaorbitalresponse import RpaOrbitalResponse
from veloxchem.gradientdriver import GradientDriver


class TestTddftXCgrad:

    def run_tddft_xcgrad(self, xcfun_label, tamm_dancoff, ref_xcgrad):

        molecule_string = """
            O   0.0   0.0   0.0
            H   0.0   1.4   1.1
            H   0.0  -1.4   1.1
        """
        basis_set_label = 'def2-svp'

        molecule = Molecule.read_str(molecule_string, units='au')
        basis = MolecularBasis.read(molecule, basis_set_label)

        # SCF

        scf_drv = ScfRestrictedDriver()

        scf_drv.dft = True
        scf_drv.xcfun = xcfun_label
        scf_drv.conv_thresh = 1e-8
        scf_drv.checkpoint_file = None

        scf_drv.compute(molecule, basis)

        # linear response and orbital response

        if tamm_dancoff:
            rsp_solver = TDAExciDriver()
            orbrsp_solver = TdaOrbitalResponse()
        else:
            rsp_solver = LinearResponseEigenSolver()
            orbrsp_solver = RpaOrbitalResponse()

        rsp_solver.dft = True
        rsp_solver.xcfun = scf_drv.xcfun
        rsp_solver.conv_thresh = 1e-5
        rsp_solver.nstates = 3
        rsp_solver.checkpoint_file = None

        rsp_results = rsp_solver.compute(molecule, basis, scf_drv.scf_tensors)

        orbrsp_solver.dft = True
        orbrsp_solver.xcfun = rsp_solver.xcfun
        orbrsp_solver.nstates = 3

        orbrsp_results = orbrsp_solver.compute(molecule, basis,
                                               scf_drv.scf_tensors, rsp_results)

        # XC contribution to molecular gradient

        grad_drv = GradientDriver(scf_drv)

        if is_mpi_master():
            gs_dm = scf_drv.scf_tensors['D_alpha']
            gs_density = AODensityMatrix([gs_dm], denmat.rest)

            rhow_dm = 0.5 * orbrsp_results['relaxed_density_ao']
            rhow_dm_sym = 0.5 * (rhow_dm + rhow_dm.T)
            rhow_den_sym = AODensityMatrix([rhow_dm_sym], denmat.rest)

            xmy = orbrsp_results['x_minus_y_ao']
            xmy_sym = 0.5 * (xmy + xmy.T)
            xmy_den_sym = AODensityMatrix([xmy_sym], denmat.rest)

        else:
            gs_density = AODensityMatrix()
            rhow_den_sym = AODensityMatrix()
            xmy_den_sym = AODensityMatrix()

        gs_density.broadcast(grad_drv.rank, grad_drv.comm)
        rhow_den_sym.broadcast(grad_drv.rank, grad_drv.comm)
        xmy_den_sym.broadcast(grad_drv.rank, grad_drv.comm)

        vxc_contrib = grad_drv.grad_vxc_contrib(molecule, basis, rhow_den_sym,
                                                gs_density, xcfun_label)

        vxc_contrib_2 = grad_drv.grad_vxc2_contrib(molecule, basis,
                                                   rhow_den_sym, gs_density,
                                                   gs_density, xcfun_label)

        vxc2_contrib = grad_drv.grad_vxc2_contrib(molecule, basis, xmy_den_sym,
                                                  xmy_den_sym, gs_density,
                                                  xcfun_label)

        vxc2_contrib_2 = grad_drv.grad_vxc3_contrib(molecule, basis,
                                                    xmy_den_sym, xmy_den_sym,
                                                    gs_density, xcfun_label)

        xcgrad = grad_drv.grad_tddft_xc_contrib(molecule, basis, rhow_den_sym,
                                                xmy_den_sym, gs_density,
                                                xcfun_label)

        if is_mpi_master():
            assert np.max(np.abs(xcgrad - ref_xcgrad)) < 1.0e-5
            xcgrad2 = vxc_contrib + vxc_contrib_2 + vxc2_contrib + vxc2_contrib_2
            assert np.max(np.abs(xcgrad2 - ref_xcgrad)) < 1.0e-5

    def test_tda_xcgrad_slater(self):

        xcfun_label = 'slater'
        tamm_dancoff = True
        ref_xcgrad = np.array(
            [[-7.35271097e-17, -7.30008259e-15, 1.10204564e-01],
             [2.66183372e-16, -5.78669206e-02, -5.51022819e-02],
             [-3.50947203e-16, 5.78669206e-02, -5.51022819e-02]])

        self.run_tddft_xcgrad(xcfun_label, tamm_dancoff, ref_xcgrad)

    def test_rpa_xcgrad_slater(self):

        xcfun_label = 'slater'
        tamm_dancoff = False
        ref_xcgrad = np.array(
            [[-4.57598797e-17, 1.16411472e-14, 1.08924175e-01],
             [2.37772515e-16, -5.77605000e-02, -5.44620872e-02],
             [-2.69921252e-16, 5.77605000e-02, -5.44620872e-02]])

        self.run_tddft_xcgrad(xcfun_label, tamm_dancoff, ref_xcgrad)

    def test_tda_xcgrad_blyp(self):

        xcfun_label = 'blyp'
        tamm_dancoff = True
        ref_xcgrad = np.array(
            [[-1.48234494e-16, 2.59450523e-14, 1.11233752e-01],
             [-1.60015912e-16, -5.75591566e-02, -5.56168761e-02],
             [2.13936447e-16, 5.75591566e-02, -5.56168761e-02]])

        self.run_tddft_xcgrad(xcfun_label, tamm_dancoff, ref_xcgrad)

    def test_rpa_xcgrad_blyp(self):

        xcfun_label = 'blyp'
        tamm_dancoff = False
        ref_xcgrad = np.array(
            [[-1.27004678e-16, 1.68531512e-14, 1.09730988e-01],
             [-1.45607378e-16, -5.74167800e-02, -5.48654938e-02],
             [2.04927911e-16, 5.74167800e-02, -5.48654938e-02]])

        self.run_tddft_xcgrad(xcfun_label, tamm_dancoff, ref_xcgrad)
