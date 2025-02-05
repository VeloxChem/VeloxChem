from pathlib import Path
import numpy as np
import h5py

from veloxchem.veloxchemlib import mpi_master
from veloxchem.veloxchemlib import parse_xc_func
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.gradientdriver import GradientDriver


class TestTddftXCgrad:

    def run_tddft_xcgrad(self, xcfun_label, tamm_dancoff, ref_xcgrad):

        molecule_string = """
            O   0.0   0.0   0.0
            H   0.0   1.4   1.1
            H   0.0  -1.4   1.1
        """
        basis_set_label = 'def2-svp'

        molecule = Molecule.read_molecule_string(molecule_string, units='au')
        basis = MolecularBasis.read(molecule, basis_set_label, ostream=None)

        grad_drv = GradientDriver()
        grad_drv.xcfun = parse_xc_func(xcfun_label)
        grad_drv.ostream.mute()

        if grad_drv.rank == mpi_master():
            here = Path(__file__).parent
            h5file = here / 'data' / 'tddft_xcgrad.dens.h5'

            hf = h5py.File(h5file, 'r')

            td_name = 'tda' if tamm_dancoff else 'rpa'
            dset_name = f'{xcfun_label}_{td_name}'

            gs_dm = np.array(hf.get(f'{dset_name}_D_alpha'))

            rhow_dm = 0.5 * np.array(hf.get(f'{dset_name}_relaxed_density_ao'))
            rhow_dm_sym = 0.5 * (rhow_dm + rhow_dm.T)

            xmy = np.array(hf.get(f'{dset_name}_x_minus_y_ao'))
            xmy_sym = 0.5 * (xmy + xmy.T)

            hf.close()

        else:
            gs_dm = None
            rhow_dm_sym = None
            xmy_sym = None

        gs_dm = grad_drv.comm.bcast(gs_dm, root=mpi_master())
        rhow_dm_sym = grad_drv.comm.bcast(rhow_dm_sym, root=mpi_master())
        xmy_sym = grad_drv.comm.bcast(xmy_sym, root=mpi_master())

        vxc_contrib = grad_drv.grad_vxc_contrib(molecule, basis, [rhow_dm_sym],
                                                [gs_dm], xcfun_label)

        vxc_contrib_2 = grad_drv.grad_vxc2_contrib(molecule, basis,
                                                   [rhow_dm_sym], [gs_dm],
                                                   [gs_dm], xcfun_label)

        vxc2_contrib = grad_drv.grad_vxc2_contrib(molecule, basis, [xmy_sym],
                                                  [xmy_sym], [gs_dm],
                                                  xcfun_label)

        vxc2_contrib_2 = grad_drv.grad_vxc3_contrib(molecule, basis, [xmy_sym],
                                                    [xmy_sym], [gs_dm],
                                                    xcfun_label)

        xcgrad = grad_drv.grad_tddft_xc_contrib(molecule, basis, [rhow_dm_sym],
                                                [xmy_sym], [gs_dm], xcfun_label)

        if grad_drv.rank == mpi_master():
            xcgrad2 = vxc_contrib + vxc_contrib_2 + vxc2_contrib + vxc2_contrib_2

            assert np.max(np.abs(xcgrad - ref_xcgrad)) < 1.0e-5
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
