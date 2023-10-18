import numpy as np
import time
import sys
from mpi4py import MPI

from .veloxchemlib import ElectronRepulsionIntegralsDriver
from .veloxchemlib import ElectricDipoleIntegralsDriver
from .veloxchemlib import AODensityMatrix
from .veloxchemlib import AOFockMatrix
from .veloxchemlib import mpi_master
from .veloxchemlib import denmat
from .veloxchemlib import fockmat
from .veloxchemlib import XCIntegrator
from .veloxchemlib import MolecularGrid
from .veloxchemlib import parse_xc_func
from .veloxchemlib import GridDriver, XCMolecularGradient
from .cphfsolver import CphfSolver
from .polorbitalresponse import PolOrbitalResponse
from .outputstream import OutputStream
from .qqscheme import get_qq_scheme
from .errorhandler import assert_msg_critical
from .inputparser import parse_input
from .sanitychecks import dft_sanity_check

# For PySCF integral derivatives
from .import_from_pyscf import overlap_deriv
from .import_from_pyscf import hcore_deriv
#from .import_from_pyscf import fock_deriv
from .import_from_pyscf import eri_deriv
from .import_from_pyscf import dipole_deriv

# TODO: remove commented out code, use flake8 and yapf

class PolarizabilityGradient():
    """
    Implements the (analytical) gradient of the dipole polarizability
    needed for the calculation of vibrational Raman intensities.
    """

    def __init__(self, comm=None, ostream=None):
        """
        Initializes polarizability gradient driver to default setup.

        :param comm:
            The MPI communicator.
        :param ostream:
            The output stream.
        """

        if comm is None:
            comm = MPI.COMM_WORLD

        if ostream is None:
            ostream = OutputStream(sys.stdout)

        # MPI information
        self.comm = comm
        self.rank = self.comm.Get_rank()
        self.nodes = self.comm.Get_size()

        self.ostream = ostream

        self.polgradient = None

        # DFT information
        self._dft = False
        self.grid_level = 4
        self.xcfun = None

        # Polarizability information
        self.flag = 'Polarizability gradient'
        self.frequencies = (0,)
        self.vector_components = 'xyz'

        self._input_keywords = {
            'gradient': {
                'vector_components': ('str_lower', 'Cartesian components of operator'),
                'frequencies': ('seq_range', 'frequencies'),
            },
            'method_settings': {
                'xcfun': ('str_upper', 'exchange-correlation functional'),
                'grid_level': ('int', 'accuracy level of DFT grid'),
            }
        }

    def update_settings(self, grad_dict, orbrsp_dict=None, method_dict=None):
        """
        Updates response and method settings in polarizability gradient
        computation driver.

        :param grad_dict:
            The input dictionary of gradient input.
        :param orbrsp_dict:
            The dictionary of orbital response input.
        :param method_dict:
            The dictionary of method settings.
        """

        if method_dict is None:
            method_dict = {}
        if orbrsp_dict is None:
            orbrsp_dict = {}

        grad_keywords = {
            key: val[0] for key, val in self._input_keywords['gradient'].items()
        }

        parse_input(self, grad_keywords, grad_dict)

        method_keywords = {
            key: val[0]
            for key, val in self._input_keywords['method_settings'].items()
        }
        
        parse_input(self, method_keywords, method_dict)

        dft_sanity_check(self, 'update_settings')
        
        if 'frequencies' not in orbrsp_dict:
            orbrsp_dict['frequencies'] = self.frequencies

        self.method_dict = dict(method_dict)
        self.orbrsp_dict = dict(orbrsp_dict)

    def compute(self, molecule, basis, scf_tensors, lr_results):
        """
        Performs calculation of analytical polarizability gradient.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        :param scf_tensors:
            The tensors from the converged SCF calculation.
        :param lr_results:
            The results of the linear response calculation.

        :return:
            A dictionary containing the polarizability gradient
            with the input frequencies as keys.
            A dictionary containing the results from the orbital
            response calculations.
        """

        # DEBUG
        if self.rank == mpi_master():
            self.ostream.print_info('GRAD DICT FREQS: {}'.format(self.frequencies))
            self.ostream.print_info('ORBRSP DICT FREQS: {}'.format(
                self.orbrsp_dict['frequencies']
            ))

        # convert sequence of frequencies to list []
        frequencylist = []
        for w in self.frequencies:
            frequencylist.append(w)
        self.frequencies = frequencylist
        self.orbrsp_dict['frequencies'] = frequencylist
        
        # compute orbital response
        orbrsp_drv = PolOrbitalResponse(self.comm, self.ostream)
        orbrsp_drv.update_settings(self.orbrsp_dict, self.method_dict)
        orbrsp_drv.compute(molecule, basis, scf_tensors, lr_results)
        orbrsp_drv.compute_omega(molecule, basis, scf_tensors, lr_results)
        all_orbrsp_results = orbrsp_drv.cphf_results
        self.ostream.print_blank()
        self.ostream.print_info('Finished orbital response calculations')
        self.ostream.flush()

        polgrad_results ={}
        n_freqs = len(self.frequencies)
        for f, w in enumerate(self.frequencies):
            self.ostream.print_info(
                'Building gradient for frequency = {:4.3f}'.format(w)
                )
            self.ostream.flush()
            
            orbrsp_results = all_orbrsp_results[w]
            if self.rank == mpi_master():
                mo = scf_tensors['C'] # only alpha part
                ovlp = scf_tensors['S']
                nao = mo.shape[0]
                nocc = molecule.number_of_alpha_electrons()
                mo_occ = mo[:, :nocc].copy()
                mo_vir = mo[:, nocc:].copy()
                nvir = mo_vir.shape[1]
                natm = molecule.number_of_atoms()

                # only alpha part:
                gs_dm = scf_tensors['D_alpha']
                # TODO: MPI should this be done before loop over freqs?
                cphf_ov = all_orbrsp_results['cphf_ov']
                dof = int(cphf_ov.shape[0] / n_freqs)
                lambda_mo = cphf_ov.reshape(n_freqs, dof, nocc, nvir)[f]

                # spin summation already included:
                x_plus_y = orbrsp_results['x_plus_y_ao']
                x_minus_y = orbrsp_results['x_minus_y_ao']
                dof = x_plus_y.shape[0] # Number of vector components
                omega_ao = orbrsp_results['omega_ao'].reshape(dof, dof, nao, nao)

                # transform lambda multipliers to AO basis and
                # calculate relaxed density matrix
                lambda_ao = np.einsum('mi,xia,na->xmn', mo_occ, lambda_mo,
                                      mo_vir).reshape(dof, dof, nao, nao) # occ-vir
                lambda_ao += lambda_ao.transpose(0,1,3,2) # vir-occ
                rel_dm_ao = orbrsp_results['unrel_dm_ao'] + lambda_ao 

                # analytical polarizability gradient
                pol_gradient = np.zeros((dof, dof, natm, 3))
                # dictionary to translate from numbers to operator components 'xyz'
                component_dict = {0: 'x', 1: 'y', 2: 'z'}

                if self._dft:
                    if self.xcfun.is_hybrid():
                        frac_K = self.xcfun.get_frac_exact_exchange()
                    else:
                        frac_K = 0.0
                else:
                    frac_K = 1.0


                # loop over atoms and contract integral derivatives
                # with density matrices
                # add the corresponding contribution to the gradient
                for i in range(natm):
                    # taking integral derivatives from pyscf
                    d_ovlp = overlap_deriv(molecule, basis, i)
                    d_hcore = hcore_deriv(molecule, basis, i)
                    d_eri = eri_deriv(molecule, basis, i)
                    d_dipole = dipole_deriv(molecule, basis, i)

                    # Calculate the analytic polarizability gradient
                    pol_gradient[:, :, i] += ( np.einsum('xymn,amn->xya',
                                                          2.0 * rel_dm_ao, d_hcore)
                        + 1.0 * np.einsum('xymn,amn->xya', 2.0 * omega_ao, d_ovlp)
                        + 2.0 * np.einsum('mt,xynp,amtnp->xya',
                                           gs_dm, 2.0 * rel_dm_ao, d_eri)
                        - 1.0 * frac_K * np.einsum('mt,xynp,amnpt->xya',
                                                    gs_dm, 2.0 * rel_dm_ao, d_eri)
                        + 1.0 * np.einsum('xmn,ypt,atpmn->xya',
                                           x_plus_y, x_plus_y - x_plus_y.transpose(0,2,1), d_eri)
                        - 0.5 * frac_K * np.einsum('xmn,ypt,atnmp->xya',
                                            x_plus_y, x_plus_y - x_plus_y.transpose(0,2,1), d_eri)
                        + 1.0 * np.einsum('xmn,ypt,atpmn->xya',
                                           x_minus_y, x_minus_y + x_minus_y.transpose(0,2,1), d_eri)
                        - 0.5 * frac_K * np.einsum('xmn,ypt,atnmp->xya',
                                            x_minus_y, x_minus_y + x_minus_y.transpose(0,2,1), d_eri)
                        + 1.0 * np.einsum('xmn,ypt,atpmn->yxa',
                                        x_plus_y, x_plus_y - x_plus_y.transpose(0,2,1), d_eri)
                        - 0.5 * frac_K * np.einsum('xmn,ypt,atnmp->yxa',
                                            x_plus_y, x_plus_y - x_plus_y.transpose(0,2,1), d_eri)
                        + 1.0 * np.einsum('xmn,ypt,atpmn->yxa',
                                           x_minus_y, x_minus_y + x_minus_y.transpose(0,2,1), d_eri)
                        - 0.5 * frac_K * np.einsum('xmn,ypt,atnmp->yxa',
                                        x_minus_y, x_minus_y + x_minus_y.transpose(0,2,1), d_eri)
                        - 2.0 * np.einsum('xmn,yamn->xya', x_minus_y, d_dipole)
                        - 2.0 * np.einsum('xmn,yamn->yxa', x_minus_y, d_dipole)
                                    )

                # Add exchange-correlation contributions to the gradient
                # for now by looping over each component of the polarizability
                if self._dft:
                    xcfun_label = self.xcfun.get_func_label()
                    for i in range(dof):
                        # for j in range(dof): 
                        # TODO:  include as soon as non-diagonal terms are available
                        if self.rank == mpi_master():
                            gs_density = AODensityMatrix([gs_dm], denmat.rest)

                            rhow_dm = 1.0 * rel_dm_ao[i,i]
                            rhow_dm_sym = 0.5 * (rhow_dm + rhow_dm.T)
                            rhow_den_sym = AODensityMatrix([rhow_dm_sym],
                                                           denmat.rest)

                            # Takes only one vector type, but two are needed
                            # to account for the different {x,y,z} components
                            # of the response vectors so we are doing only
                            # diagonal components for now.
                            # The sqrt2 takes into account the fact that we need to
                            # symmetrize with respect to the polarizability
                            # components.
                            # (see contraction with two-electron integrals above).
                            x_minus_y_sym = np.sqrt(2) * 0.5 * (x_minus_y[i] + x_minus_y[i].T)
                            x_minus_y_den_sym = AODensityMatrix([x_minus_y_sym], denmat.rest)

                        else:
                            gs_density = AODensityMatrix()
                            rhow_den_sym = AODensityMatrix()
                            x_minus_y_den_sym = AODensityMatrix()

                        gs_density.broadcast(self.rank, self.comm)
                        rhow_den_sym.broadcast(self.rank, self.comm)
                        x_minus_y_den_sym.broadcast(self.rank, self.comm)

                        polgrad_xcgrad = self.grad_polgrad_xc_contrib(molecule,
                                                basis, rhow_den_sym, x_minus_y_den_sym,
                                                         gs_density, xcfun_label)

                        if self.rank == mpi_master():
                            pol_gradient[i,i] += polgrad_xcgrad

                polgrad_results[(w)] = pol_gradient.reshape(dof, dof, 3 * natm)

        self.ostream.print_blank()
        self.ostream.print_info(
            'Finished calculation of polarizability gradient'
            )
        self.ostream.flush()
        self.polgradient = dict(polgrad_results)
        #return polgrad_results, all_orbrsp_results

    def grad_polgrad_xc_contrib(self, molecule, ao_basis, rhow_den, x_minus_y_den,
                              gs_density, xcfun_label):
        """
        Calculates exchange-correlation contribution to polarizability gradient.

        :param molecule:
            The molecule.
        :param ao_basis:
            The AO basis set.
        :param rhow_den:
            The perturbed density.
        :param x_minus_y_den:
            The X-Y density.
        :param gs_density:
            The ground state density.
        :param xcfun_label:
            The label of the xc functional.

        :return:
            The exchange-correlation contribution to polarizability gradient.
        """

        grid_drv = GridDriver(self.comm)
        grid_drv.set_level(self.grid_level)
        mol_grid = grid_drv.generate(molecule)

        xcgrad_drv = XCMolecularGradient(self.comm)
        polgrad_xcgrad = xcgrad_drv.integrate_vxc_gradient(
            molecule, ao_basis, rhow_den, gs_density, mol_grid, xcfun_label)
        polgrad_xcgrad += xcgrad_drv.integrate_fxc_gradient(
            molecule, ao_basis, rhow_den, gs_density, gs_density, mol_grid,
            xcfun_label)
        polgrad_xcgrad += xcgrad_drv.integrate_fxc_gradient(
            molecule, ao_basis, x_minus_y_den, x_minus_y_den, gs_density, mol_grid,
            xcfun_label)
        polgrad_xcgrad += xcgrad_drv.integrate_kxc_gradient(
            molecule, ao_basis, x_minus_y_den, x_minus_y_den, gs_density, mol_grid,
            xcfun_label)
        #polgrad_xcgrad = xcgrad_drv.integrate_tddft_gradient(
        #    rhow_den, x_minus_y_den, gs_density, molecule, ao_basis, mol_grid,
        #    xcfun_label)
        polgrad_xcgrad = self.comm.reduce(polgrad_xcgrad, root=mpi_master())

        return polgrad_xcgrad

    def _init_dft(self, molecule, scf_tensors):
        """
        Initializes DFT.

        :param molecule:
            The molecule.
        :param scf_tensors:
            The dictionary of tensors from converged SCF wavefunction.

        :return:
            The dictionary of DFT information.
        """

        # generate integration grid
        if self._dft:
            grid_drv = GridDriver(self.comm)
            grid_drv.set_level(self.grid_level)

            grid_t0 = tm.time()
            molgrid = grid_drv.generate(molecule)

            n_grid_points = molgrid.number_of_points()
            self.ostream.print_info(
                'Molecular grid with {0:d} points generated in {1:.2f} sec.'.
                format(n_grid_points,
                       tm.time() - grid_t0))
            self.ostream.print_blank()

            if self.rank == mpi_master():
                gs_density = AODensityMatrix([scf_tensors['D_alpha']],
                                             denmat.rest)
            else:
                gs_density = AODensityMatrix()
            gs_density.broadcast(self.rank, self.comm)
            #molgrid.broadcast(self.rank, self.comm) # TODO duble check

            dft_func_label = self.xcfun.get_func_label().upper()
        else:
            molgrid = MolecularGrid()
            gs_density = AODensityMatrix()
            dft_func_label = 'HF'

        return {
            'molgrid': molgrid,
            'gs_density': gs_density,
            'dft_func_label': dft_func_label,
        }

