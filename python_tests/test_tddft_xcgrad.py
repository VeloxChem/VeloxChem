import numpy as np

from veloxchem.veloxchemlib import denmat
from veloxchem.veloxchemlib import is_mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.aodensitymatrix import AODensityMatrix
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.gradientdriver import GradientDriver
from veloxchem.tdaexcidriver import TDAExciDriver
from veloxchem.lreigensolver import LinearResponseEigenSolver
from veloxchem.tdaorbitalresponse import TdaOrbitalResponse
from veloxchem.rpaorbitalresponse import RpaOrbitalResponse


class TestOrbitalResponse:

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

        gs_dm = scf_drv.scf_tensors['D_alpha']
        gs_density = AODensityMatrix([gs_dm], denmat.rest)

        rhow_dm = 0.5 * orbrsp_results['relaxed_density_ao']
        rhow_dm_sym = 0.5 * (rhow_dm + rhow_dm.T)
        rhow_den_sym = AODensityMatrix([rhow_dm_sym], denmat.rest)

        vxc_contrib = grad_drv.grad_vxc_contrib(molecule, basis, rhow_den_sym,
                                                gs_density, xcfun_label)
        vxc_contrib_2 = grad_drv.grad_fxc_contrib(molecule, basis, rhow_den_sym,
                                                  gs_density, gs_density,
                                                  xcfun_label)

        xmy = orbrsp_results['x_minus_y_ao']
        xmy_sym = 0.5 * (xmy + xmy.T)
        xmy_den_sym = AODensityMatrix([xmy_sym], denmat.rest)

        fxc_contrib = grad_drv.grad_fxc_contrib(molecule, basis, xmy_den_sym,
                                                xmy_den_sym, gs_density,
                                                xcfun_label)
        fxc_contrib_2 = grad_drv.grad_gxc_contrib(molecule, basis, xmy, xmy,
                                                  gs_density, xcfun_label)

        if is_mpi_master():
            xcgrad = vxc_contrib
            xcgrad += vxc_contrib_2
            xcgrad += fxc_contrib
            xcgrad += fxc_contrib_2
            assert np.max(np.abs(xcgrad - ref_xcgrad)) < 1.0e-5

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
