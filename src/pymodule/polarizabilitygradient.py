import numpy as np
import time as tm
import sys
from mpi4py import MPI

# from .veloxchemlib import ElectronRepulsionIntegralsDriver
# from .veloxchemlib import ElectricDipoleIntegralsDriver
from .veloxchemlib import AODensityMatrix
# from .veloxchemlib import AOFockMatrix
from .veloxchemlib import mpi_master
from .veloxchemlib import denmat
# from .veloxchemlib import fockmat
# from .veloxchemlib import XCIntegrator
from .veloxchemlib import MolecularGrid
# from .veloxchemlib import parse_xc_func
from .veloxchemlib import GridDriver, XCMolecularGradient
from .veloxchemlib import hartree_in_wavenumber
# from .cphfsolver import CphfSolver
from .polorbitalresponse import PolOrbitalResponse
from .lrsolver import LinearResponseSolver
from .cppsolver import ComplexResponse
from .molecule import Molecule
from .outputstream import OutputStream
# from .qqscheme import get_qq_scheme
# from .errorhandler import assert_msg_critical
from .inputparser import parse_input
from .sanitychecks import dft_sanity_check, polgrad_sanity_check

# For PySCF integral derivatives
from .import_from_pyscf import overlap_deriv
from .import_from_pyscf import hcore_deriv
# from .import_from_pyscf import fock_deriv
from .import_from_pyscf import eri_deriv
from .import_from_pyscf import dipole_deriv


class PolarizabilityGradient():
    """
    Implements the gradient of the dipole polarizability
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

        self.comm = comm
        self.rank = self.comm.Get_rank()
        self.nodes = self.comm.Get_size()

        self.ostream = ostream

        self.polgradient = None
        self.delta_h = 0.001

        self.is_complex = False
        self.damping = 1000.0 / hartree_in_wavenumber() # for numerical complex gradient

        self.numerical = False
        self.do_four_point = False

        self._dft = False
        self.grid_level = 4
        self.xcfun = None

        self.flag = 'Polarizability Gradient'
        self.frequencies = (0,)
        self.vector_components = 'xyz'

        self._input_keywords = {
            'polarizabilitygradient': {
                'vector_components': ('str_lower', 'Cartesian components of operator'),
                'frequencies': ('seq_range', 'frequencies'),
                'numerical': ('bool', 'do numerical integration'),
                'do_four_point': ('bool', 'do four-point numerical integration'),
                'delta_h': ('float', 'the displacement for finite difference'),
                'is_complex': ('bool', 'whether the polarizability is complex'),
                'damping': ('float', 'damping parameter for complex numerical'),
            },
            'method_settings': {
                'xcfun': ('str_upper', 'exchange-correlation functional'),
                'grid_level': ('int', 'accuracy level of DFT grid'),
            }
        }

    def update_settings(self, grad_dict, orbrsp_dict=None, method_dict=None,
                        scf_drv=None):
        """
        Updates response and method settings in polarizability gradient
        computation driver.

        :param grad_dict:
            The input dictionary of gradient input.
        :param orbrsp_dict:
            The dictionary of orbital response (CPHF) input.
        :param method_dict:
            The dictionary of method settings.
        :param scf_drv:
            The SCF driver (only for numerical calculations)
        """

        if method_dict is None:
            method_dict = {}
        if orbrsp_dict is None:
            orbrsp_dict = {}

        grad_keywords = {
            key: val[0] for key, val in 
            self._input_keywords['polarizabilitygradient'].items()
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

        # must be set without an if-statement in case the orbital response
        # dictionary has previously been set
        orbrsp_dict['is_complex'] = self.is_complex

        self.method_dict = dict(method_dict)
        self.orbrsp_dict = dict(orbrsp_dict)
        self.scf_drv = scf_drv

    def compute(self, molecule, basis, scf_tensors, lr_results):
        """
        Calls the correct function to perform the calculation of 
        the polarizability gradient.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        :param scf_tensors:
            The tensors from the converged SCF calculation.
        :param lr_results:
            The results of the linear response calculation.
        """

        # sanity checks
        dft_sanity_check(self, 'compute')
        polgrad_sanity_check(self, self.flag, lr_results)

        start_time = tm.time()

        if self.numerical:
            # print only on the master node
            if self.rank == mpi_master():
                valstr = '*** Calculating numerical polarizability gradient ***'
                self.ostream.print_header(valstr)
                self.ostream.print_blank()
                self.ostream.flush()
            if self.scf_drv is None:
                error_message = 'PolarizabilityGradient: missing input SCF driver '
                error_message += 'for numerical calculations'
                raise ValueError(error_message)
            if self.is_complex:
                self.compute_numerical_complex(molecule, basis, self.scf_drv)
            else:
                self.compute_numerical_real(molecule, basis, self.scf_drv)
        else:
            # print only on the master node
            if self.rank == mpi_master():
                valstr = '*** Calculating analytical polarizability gradient ***'
                self.ostream.print_header(valstr)
                self.ostream.print_blank()
                self.ostream.flush()
            if self.is_complex:
                self.compute_analytical_complex(molecule, basis, scf_tensors, lr_results)
            else:
                self.compute_analytical_real(molecule, basis, scf_tensors, lr_results)

        if self.rank == mpi_master():
            self.print_geometry(molecule)
            if self.numerical:
                self.print_polarizability_gradient(molecule)
            else:
                self.print_polarizability_gradient(molecule)
            valstr = '*** Time spent in polarizability gradient calculation: '
            valstr += '{:.6f} sec ***'.format(tm.time() - start_time)
            self.ostream.print_header(valstr)
            self.ostream.print_blank()
            self.ostream.flush()

    def compute_analytical_real(self, molecule, basis, scf_tensors, lr_results):
        """
        Performs calculation of the real analytical polarizability gradient.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        :param scf_tensors:
            The tensors from the converged SCF calculation.
        :param lr_results:
            The results of the linear response calculation.
        """

        # convert sequence of frequencies to list []
        frequencylist = []
        for w in self.frequencies:
            frequencylist.append(w)
        self.frequencies = frequencylist
        self.orbrsp_dict['frequencies'] = frequencylist

        # compute orbital response
        orbrsp_start_time = tm.time()

        orbrsp_drv = PolOrbitalResponse(self.comm, self.ostream)
        orbrsp_drv.update_settings(self.orbrsp_dict, self.method_dict)
        orbrsp_drv.compute(molecule, basis, scf_tensors, lr_results)
        orbrsp_drv.compute_omega(molecule, basis, scf_tensors, lr_results)
        all_orbrsp_results = orbrsp_drv.cphf_results

        valstr = '** Time spent on orbital response for {} frequencies: '.format(
            len(self.frequencies))
        valstr += '{:.6f} sec **'.format(tm.time() - orbrsp_start_time)
        self.ostream.print_header(valstr)
        self.ostream.print_blank()
        self.ostream.flush()

        polgrad_results = {}
        n_freqs = len(self.frequencies)

        loop_start_time = tm.time()

        for f, w in enumerate(self.frequencies):
            info_msg = 'Building gradient for frequency = {:4.3f}'.format(w)
            self.ostream.print_info(info_msg)
            self.ostream.flush()
            self.ostream.print_blank()

            if self.rank == mpi_master():
                orbrsp_results = all_orbrsp_results[w]
                mo = scf_tensors['C']  # only alpha part
                # ovlp = scf_tensors['S']
                nao = mo.shape[0]
                nocc = molecule.number_of_alpha_electrons()
                mo_occ = mo[:, :nocc].copy()
                mo_vir = mo[:, nocc:].copy()
                nvir = mo_vir.shape[1]
                natm = molecule.number_of_atoms()

                # only alpha part:
                gs_dm = scf_tensors['D_alpha']
                # TODO: MPI should this be done before loop over freqs?
                #cphf_ov = all_orbrsp_results['cphf_ov']
                #dof = int(cphf_ov.shape[0] / n_freqs)
                #lambda_mo = cphf_ov.reshape(n_freqs, dof, nocc, nvir)[f]
                lambda_mo = orbrsp_results['lambda_mo']

                # spin summation already included:
                x_plus_y = orbrsp_results['x_plus_y_ao']
                x_minus_y = orbrsp_results['x_minus_y_ao']
                dof = x_plus_y.shape[0]  # Number of vector components
                omega_ao = orbrsp_results['omega_ao'].reshape(dof, dof, nao, nao)

                # transform lambda multipliers to AO basis and
                # calculate relaxed density matrix
                lambda_ao = np.array([
                    np.linalg.multi_dot([mo_occ, lambda_mo[x], mo_vir.T])
                    for x in range(dof**2)                    
                ]).reshape(dof, dof, nao, nao)
                
                lambda_ao += lambda_ao.transpose(0,1,3,2)  # vir-occ
                rel_dm_ao = orbrsp_results['unrel_dm_ao'] + lambda_ao

                if self._dft:
                    if self.xcfun.is_hybrid():
                        frac_K = self.xcfun.get_frac_exact_exchange()
                    else:
                        frac_K = 0.0
                else:
                    frac_K = 1.0

                pol_gradient = np.zeros((dof, dof, natm, 3))
                # loop over atoms and contract integral derivatives
                # with density matrices
                # add the corresponding contribution to the gradient
                for i in range(natm):

                    integral_start_time = tm.time()

                    # importing integral derivatives from pyscf
                    d_ovlp = overlap_deriv(molecule, basis, i)
                    d_hcore = hcore_deriv(molecule, basis, i)
                    d_eri = eri_deriv(molecule, basis, i)
                    d_dipole = dipole_deriv(molecule, basis, i)

                    valstr = ' * Time spent importing integrals for atom #{}: '.format(i+1)
                    valstr += '{:.6f} sec * '.format(tm.time() - integral_start_time)
                    self.ostream.print_header(valstr)
                    self.ostream.print_blank()
                    self.ostream.flush()

                    gradient_start_time = tm.time()
                    # Calculate the analytic polarizability gradient
                    for x in range(dof):
                        for y in range(dof):
                            for a in range(3):
                                pol_gradient[x, y, i, a] += (
                                np.linalg.multi_dot([
                                    2.0 * rel_dm_ao[x,y].reshape(nao**2), 
                                    d_hcore[a].reshape(nao**2)
                                ])
                                + 1.0 * np.linalg.multi_dot([
                                    2.0 * omega_ao[x,y].reshape(nao**2), 
                                    d_ovlp[a].reshape(nao*nao)
                                ])
                                + 2.0 * np.linalg.multi_dot([
                                    gs_dm.reshape(nao**2),
                                    d_eri[a].reshape(nao**2,nao**2), 
                                    2.0 * rel_dm_ao[x,y].reshape(nao**2)
                                ])
                                - 1.0 * frac_K * np.linalg.multi_dot([
                                    gs_dm.reshape(nao**2),
                                    d_eri[a].transpose(0,3,1,2).reshape(nao**2,nao**2),
                                    2.0 * rel_dm_ao[x,y].reshape(nao**2)
                                ])
                                + 1.0* np.linalg.multi_dot([
                                    x_plus_y[x].reshape(nao**2), 
                                    d_eri[a].transpose(2,3,1,0).reshape(nao**2,nao**2), 
                                    (x_plus_y[y] - x_plus_y[y].T).reshape(nao**2)
                                ])
                                - 0.5 * frac_K * np.linalg.multi_dot([
                                    x_plus_y[x].reshape(nao**2),
                                    d_eri[a].transpose(2,1,3,0).reshape(nao**2,nao**2),
                                    (x_plus_y[y] - x_plus_y[y].T).reshape(nao**2)
                                ])
                                + 1.0 * np.linalg.multi_dot([
                                    x_minus_y[x].reshape(nao**2),
                                    d_eri[a].transpose(2,3,1,0).reshape(nao**2,nao**2),
                                    (x_minus_y[y] + x_minus_y[y].T).reshape(nao**2)
                                ])
                                - 0.5 * frac_K * np.linalg.multi_dot([
                                    x_minus_y[x].reshape(nao**2),
                                    d_eri[a].transpose(2,1,3,0).reshape(nao**2,nao**2),
                                    (x_minus_y[y] + x_minus_y[y].T).reshape(nao**2)
                                ])
                                - 2.0 * np.linalg.multi_dot([
                                    x_minus_y[x].reshape(nao**2),
                                    d_dipole[y,a].reshape(nao**2)
                                ])
                                )

                                pol_gradient[y, x, i, a] += (
                                1.0 * np.linalg.multi_dot([
                                    (x_plus_y[y] - x_plus_y[y].T).reshape(nao**2),
                                    d_eri[a].transpose(1,0,2,3).reshape(nao**2,nao**2), x_plus_y[x].reshape(nao**2)
                                ])
                                - 0.5 * frac_K * np.linalg.multi_dot([
                                    (x_plus_y[y] - x_plus_y[y].T).reshape(nao**2),
                                    d_eri[a].transpose(3,0,2,1).reshape(nao**2,nao**2), x_plus_y[x].reshape(nao**2)
                                ])
                                + 1.0 * np.linalg.multi_dot([
                                    (x_minus_y[y] + x_minus_y[y].T).reshape(nao**2),
                                    d_eri[a].transpose(1,0,2,3).reshape(nao**2,nao**2), x_minus_y[x].reshape(nao**2)
                                ])
                                - 0.5 * frac_K * np.linalg.multi_dot([
                                    (x_minus_y[y] + x_minus_y[y].T).reshape(nao**2),
                                    d_eri[a].transpose(3,0,2,1).reshape(nao**2,nao**2), x_minus_y[x].reshape(nao**2)
                                ])
                                - 2.0 * np.linalg.multi_dot([
                                d_dipole[y,a].reshape(nao**2), x_minus_y[x].reshape(nao**2)
                                ])
                                )

                    valstr = ' * Time spent calculating pol. gradient for '
                    valstr += 'atom #{:d}: {:.6f} sec * '.format(
                        (i + 1), tm.time() - gradient_start_time)
                    self.ostream.print_header(valstr)
                    self.ostream.print_blank()
                    self.ostream.flush()

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

                        polgrad_xcgrad = self.grad_polgrad_xc_contrib(
                            molecule, basis, rhow_den_sym, x_minus_y_den_sym,
                            gs_density, xcfun_label)

                        if self.rank == mpi_master():
                            pol_gradient[i,i] += polgrad_xcgrad

                polgrad_results[w] = pol_gradient.reshape(dof, dof, 3 * natm)

        self.polgradient = dict(polgrad_results)

        valstr = '** Time spent on constructing the gradient for '
        valstr += '{:d} frequencies: {:.6f} sec **'.format(
            n_freqs, tm.time() - loop_start_time)
        self.ostream.print_header(valstr)
        self.ostream.print_blank()
        self.ostream.flush()

    def compute_analytical_complex(self, molecule, basis, scf_tensors, lr_results):
        """
        Performs calculation of the complex analytical polarizability gradient.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        :param scf_tensors:
            The tensors from the converged SCF calculation.
        :param lr_results:
            The results of the CPP calculation.
        """

        # convert sequence of frequencies to list []
        frequencylist = []
        for w in self.frequencies:
            frequencylist.append(w)
        self.frequencies = frequencylist
        self.orbrsp_dict['frequencies'] = frequencylist

        # compute orbital response
        orbrsp_start_time = tm.time()

        orbrsp_drv = PolOrbitalResponse(self.comm, self.ostream)
        orbrsp_drv.update_settings(self.orbrsp_dict, self.method_dict)
        orbrsp_drv.compute(molecule, basis, scf_tensors, lr_results)
        orbrsp_drv.compute_omega(molecule, basis, scf_tensors, lr_results)
        all_orbrsp_results = orbrsp_drv.cphf_results

        valstr = '** Time spent on orbital response for {} frequencies: '.format(
            len(self.frequencies))
        valstr += '{:.6f} sec **'.format(tm.time() - orbrsp_start_time)
        self.ostream.print_header(valstr)
        self.ostream.print_blank()
        self.ostream.flush()

        polgrad_results = {}
        n_freqs = len(self.frequencies)

        loop_start_time = tm.time()

        for f, w in enumerate(self.frequencies):
            info_msg = 'Building gradient for frequency = {:4.3f}'.format(w)
            self.ostream.print_info(info_msg)
            self.ostream.flush()
            self.ostream.print_blank()

            if self.rank == mpi_master():
                orbrsp_results = all_orbrsp_results[w]
                mo = scf_tensors['C']  # only alpha part
                # ovlp = scf_tensors['S']
                nao = mo.shape[0]
                nocc = molecule.number_of_alpha_electrons()
                mo_occ = mo[:, :nocc].copy()
                mo_vir = mo[:, nocc:].copy()
                nvir = mo_vir.shape[1]
                natm = molecule.number_of_atoms()

                # only alpha part:
                gs_dm = scf_tensors['D_alpha']

                # spin summation already included:
                x_plus_y = orbrsp_results['x_plus_y_ao']
                x_minus_y = orbrsp_results['x_minus_y_ao']

                dof = x_plus_y.shape[0]  # Number of vector components

                omega_ao = orbrsp_results['omega_ao'].reshape(dof, dof, nao, nao)
                # cphf_ov = all_orbrsp_results['cphf_ov']
                # len_lambdas = int(cphf_ov.shape[0] / n_freqs)

                # select lambda for frequency with index f and reshape to complex
                # lambda_mo = cphf_ov.reshape(n_freqs, len_lambdas, nocc, nvir)[f]
                # lambda_mo = lambda_mo[:dof**2] + 1j * lambda_mo[dof**2:]
                lambda_mo = orbrsp_results['lambda_mo']

                # transform lambda multipliers to AO basis and
                # calculate relaxed density matrix
                lambda_ao = np.array([
                    np.linalg.multi_dot([mo_occ, lambda_mo[x], mo_vir.T])
                    for x in range(dof**2)                    
                ]).reshape(dof, dof, nao, nao)
                
                lambda_ao += lambda_ao.transpose(0,1,3,2)  # vir-occ
                rel_dm_ao = orbrsp_results['unrel_dm_ao'] + lambda_ao

                if self._dft:
                    if self.xcfun.is_hybrid():
                        frac_K = self.xcfun.get_frac_exact_exchange()
                    else:
                        frac_K = 0.0
                else:
                    frac_K = 1.0

                pol_gradient = np.zeros((dof, dof, natm, 3), dtype=np.complex_)
                # loop over atoms and contract integral derivatives
                # with density matrices
                # add the corresponding contribution to the gradient
                for i in range(natm):

                    integral_start_time = tm.time()

                    # importing integral derivatives from pyscf
                    d_ovlp = overlap_deriv(molecule, basis, i)
                    d_hcore = hcore_deriv(molecule, basis, i)
                    d_eri = eri_deriv(molecule, basis, i)
                    d_dipole = dipole_deriv(molecule, basis, i)

                    valstr = ' * Time spent importing integrals for atom #{}: '.format(i+1)
                    valstr += '{:.6f} sec * '.format(tm.time() - integral_start_time)
                    self.ostream.print_header(valstr)
                    self.ostream.print_blank()
                    self.ostream.flush()

                    gradient_start_time = tm.time()
                    # Calculate the analytic polarizability gradient
                    for x in range(dof):
                        for y in range(dof):
                            for a in range(3):
                                pol_gradient[x, y, i, a] += (
                                np.linalg.multi_dot([
                                    2.0 * rel_dm_ao[x,y].reshape(nao**2), 
                                    d_hcore[a].reshape(nao**2)
                                ])
                                + 1.0 * np.linalg.multi_dot([
                                    2.0 * omega_ao[x,y].reshape(nao**2), 
                                    d_ovlp[a].reshape(nao*nao)
                                ])
                                + 2.0 * np.linalg.multi_dot([
                                    gs_dm.reshape(nao**2),
                                    d_eri[a].reshape(nao**2,nao**2), 
                                    2.0 * rel_dm_ao[x,y].reshape(nao**2)
                                ])
                                - 1.0 * frac_K * np.linalg.multi_dot([
                                    gs_dm.reshape(nao**2),
                                    d_eri[a].transpose(0,3,1,2).reshape(nao**2,nao**2),
                                    2.0 * rel_dm_ao[x,y].reshape(nao**2)
                                ])
                                + 1.0* np.linalg.multi_dot([
                                    x_plus_y[x].reshape(nao**2), 
                                    d_eri[a].transpose(2,3,1,0).reshape(nao**2,nao**2), 
                                    (x_plus_y[y] - x_plus_y[y].T).reshape(nao**2)
                                ])
                                - 0.5 * frac_K * np.linalg.multi_dot([
                                    x_plus_y[x].reshape(nao**2),
                                    d_eri[a].transpose(2,1,3,0).reshape(nao**2,nao**2),
                                    (x_plus_y[y] - x_plus_y[y].T).reshape(nao**2)
                                ])
                                + 1.0 * np.linalg.multi_dot([
                                    x_minus_y[x].reshape(nao**2),
                                    d_eri[a].transpose(2,3,1,0).reshape(nao**2,nao**2),
                                    (x_minus_y[y] + x_minus_y[y].T).reshape(nao**2)
                                ])
                                - 0.5 * frac_K * np.linalg.multi_dot([
                                    x_minus_y[x].reshape(nao**2),
                                    d_eri[a].transpose(2,1,3,0).reshape(nao**2,nao**2),
                                    (x_minus_y[y] + x_minus_y[y].T).reshape(nao**2)
                                ])
                                - 2.0 * np.linalg.multi_dot([
                                    x_minus_y[x].reshape(nao**2),
                                    d_dipole[y,a].reshape(nao**2)
                                ])
                                )

                                pol_gradient[y, x, i, a] += (
                                1.0 * np.linalg.multi_dot([
                                    (x_plus_y[y] - x_plus_y[y].T).reshape(nao**2),
                                    d_eri[a].transpose(1,0,2,3).reshape(nao**2,nao**2), x_plus_y[x].reshape(nao**2)
                                ])
                                - 0.5 * frac_K * np.linalg.multi_dot([
                                    (x_plus_y[y] - x_plus_y[y].T).reshape(nao**2),
                                    d_eri[a].transpose(3,0,2,1).reshape(nao**2,nao**2), x_plus_y[x].reshape(nao**2)
                                ])
                                + 1.0 * np.linalg.multi_dot([
                                    (x_minus_y[y] + x_minus_y[y].T).reshape(nao**2),
                                    d_eri[a].transpose(1,0,2,3).reshape(nao**2,nao**2), x_minus_y[x].reshape(nao**2)
                                ])
                                - 0.5 * frac_K * np.linalg.multi_dot([
                                    (x_minus_y[y] + x_minus_y[y].T).reshape(nao**2),
                                    d_eri[a].transpose(3,0,2,1).reshape(nao**2,nao**2), x_minus_y[x].reshape(nao**2)
                                ])
                                - 2.0 * np.linalg.multi_dot([
                                d_dipole[y,a].reshape(nao**2), x_minus_y[x].reshape(nao**2)
                                ])
                                )

                    valstr = ' * Time spent calculating pol. gradient for '
                    valstr += 'atom #{:d}: {:.6f} sec * '.format(
                        (i + 1), tm.time() - gradient_start_time)
                    self.ostream.print_header(valstr)
                    self.ostream.print_blank()
                    self.ostream.flush()

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

                        polgrad_xcgrad = self.grad_polgrad_xc_contrib(
                            molecule, basis, rhow_den_sym, x_minus_y_den_sym,
                            gs_density, xcfun_label)

                        if self.rank == mpi_master():
                            pol_gradient[i,i] += polgrad_xcgrad

                polgrad_results[w] = pol_gradient.reshape(dof, dof, 3 * natm)

        self.polgradient = dict(polgrad_results)

        valstr = '** Time spent on constructing the gradient for '
        valstr += '{:d} frequencies: {:.6f} sec **'.format(
            n_freqs, tm.time() - loop_start_time)
        self.ostream.print_header(valstr)
        self.ostream.print_blank()
        self.ostream.flush()

    def compute_numerical_real(self, molecule, ao_basis, scf_drv):
        """
        Performs calculation of numerical nuclear gradient
        of the real electric dipole polarizability.

        Moved from tddftgradientdriver.py 18/10/2023

        :param molecule:
            The molecule.
        :param ao_basis:
            The AO basis set.
        :param scf_drv:
            The SCF driver.
        """

        # number of atoms
        natm = molecule.number_of_atoms()

        # atom labels
        labels = molecule.get_labels()

        # atom coordinates (nx3)
        coords = molecule.get_coordinates()

        # number of frequencies
        n_freqs = len(self.frequencies)

        # dictionary
        polgrad_results = {}

        # linear response driver for polarizability calculation
        lr_drv = LinearResponseSolver(self.comm, self.ostream)
        lr_drv.frequencies = self.frequencies

        # polarizability: 3 coordinates x 3 coordinates (ignoring frequencies)
        # polarizability gradient: dictionary goes through 3 coordinates
        # x 3 coordinates, each entry having values for
        # no. atoms x 3 coordinates
        num_polgradient = np.zeros((n_freqs, 3, 3, 3 * natm))

        if not self.do_four_point:
            if self.rank == mpi_master():
                self.ostream.print_blank()
                self.ostream.print_info('do_four_point: False')
                self.ostream.flush()

            # mute the outputstreams
            scf_drv.ostream.mute()
            lr_drv.ostream.mute()

            for i in range(natm):
                for d in range(3):
                    coords[i, d] += self.delta_h
                    new_mol = Molecule(labels, coords, units='au')
                    scf_drv.compute(new_mol, ao_basis)
                    lr_drv._is_converged = False
                    lr_results_p = lr_drv.compute(new_mol, ao_basis,
                                                  scf_drv.scf_tensors)

                    coords[i, d] -= 2.0 * self.delta_h
                    new_mol = Molecule(labels, coords, units='au')
                    scf_drv.compute(new_mol, ao_basis)
                    lr_drv._is_converged = False
                    lr_results_m = lr_drv.compute(new_mol, ao_basis,
                                                  scf_drv.scf_tensors)

                    coords[i, d] += self.delta_h
                    for f, w in enumerate(self.frequencies):
                        self.ostream.print_info('{}'.format(w))
                        self.ostream.flush()
                        for aop, acomp in enumerate('xyz'):
                            for bop, bcomp in enumerate('xyz'):
                                key = (acomp, bcomp, w)
                                if self.rank == mpi_master():
                                    num_polgradient[f, aop, bop, 3 * i + d] = (
                                    (lr_results_p['response_functions'][key] -
                                     lr_results_m['response_functions'][key]) /
                                    (2.0 * self.delta_h))
            if self.rank == mpi_master():
                for f, w in enumerate(self.frequencies):
                    polgrad_results[w] = num_polgradient[f]
                self.polgradient = dict(polgrad_results)
            else:
                self.polgradient = {}

        # four-point approximation for debugging of analytical gradient
        else:
            if self.rank == mpi_master():
                self.ostream.print_info('do_four_point: True')
                self.ostream.flush()
            # mute the outputstreams:
            scf_drv.ostream.mute()
            lr_drv.ostream.mute()

            for i in range(natm):
                for d in range(3):
                    coords[i, d] += self.delta_h
                    new_mol = Molecule(labels, coords, units='au')
                    scf_drv.compute(new_mol, ao_basis)
                    lr_drv._is_converged = False
                    lr_results_p1 = lr_drv.compute(new_mol, ao_basis,
                                                   scf_drv.scf_tensors)

                    coords[i, d] += self.delta_h
                    new_mol = Molecule(labels, coords, units='au')
                    scf_drv.compute(new_mol, ao_basis)
                    lr_drv._is_converged = False
                    lr_results_p2 = lr_drv.compute(new_mol, ao_basis,
                                                   scf_drv.scf_tensors)

                    coords[i, d] -= 3.0 * self.delta_h
                    new_mol = Molecule(labels, coords, units='au')
                    scf_drv.compute(new_mol, ao_basis)
                    lr_drv._is_converged = False
                    lr_results_m1 = lr_drv.compute(new_mol, ao_basis,
                                                   scf_drv.scf_tensors)

                    coords[i, d] -= 1.0 * self.delta_h
                    new_mol = Molecule(labels, coords, units='au')
                    scf_drv.compute(new_mol, ao_basis)
                    lr_drv._is_converged = False
                    lr_results_m2 = lr_drv.compute(new_mol, ao_basis,
                                                   scf_drv.scf_tensors)

                    coords[i, d] += 2.0 * self.delta_h
                    for f, w in enumerate(self.frequencies):
                        for aop, acomp in enumerate('xyz'):
                            for bop, bcomp in enumerate('xyz'):
                                # f'(x) ~ [ f(x - 2h) - 8 f(x - h)
                                # + 8 f(x + h) - f(x + 2h) ] / ( 12h )
                                key = (acomp, bcomp, w)
                                if self.rank == mpi_master():
                                    num_polgradient[f, aop, bop, 3 * i + d] = ((
                                    lr_results_m2['response_functions'][key] -
                                    8.0 * lr_results_m1['response_functions'][key] +
                                    8.0 * lr_results_p1['response_functions'][key] -
                                    lr_results_p2['response_functions'][key]) / (
                                        12.0 * self.delta_h))
            if self.rank == mpi_master():
                for f, w in enumerate(self.frequencies):
                    polgrad_results[w] = num_polgradient[f]
                self.polgradient = dict(polgrad_results)
            else:
                self.polgradient = {}

        scf_drv.ostream.unmute()
        lr_drv.ostream.unmute()

    def compute_numerical_complex(self, molecule, ao_basis, scf_drv):
        """
        Performs calculation of numerical nuclear gradient
        of the complex electric dipole polarizability.

        :param molecule:
            The molecule.
        :param ao_basis:
            The AO basis set.
        :param scf_drv:
            The SCF driver.
        """

        # number of atoms
        natm = molecule.number_of_atoms()

        # atom labels
        labels = molecule.get_labels()

        # atom coordinates (nx3)
        coords = molecule.get_coordinates()

        # number of frequencies
        n_freqs = len(self.frequencies)

        # dictionary
        polgrad_results = {}

        # linear response driver for polarizability calculation
        cpp_drv = ComplexResponse(self.comm, self.ostream)
        cpp_drv.frequencies = self.frequencies
        cpp_drv.damping = self.damping

        # polarizability: 3 coordinates x 3 coordinates (ignoring frequencies)
        # polarizability gradient: dictionary goes through 3 coordinates
        # x 3 coordinates, each entry having values for
        # no. atoms x 3 coordinates
        # datatype: complex
        num_polgradient = np.zeros((n_freqs, 3, 3, 3 * natm), dtype = complex)

        if not self.do_four_point:
            if self.rank == mpi_master():
                self.ostream.print_blank()
                self.ostream.print_info('do_four_point: False')
                self.ostream.flush()

            # mute the outputstreams
            scf_drv.ostream.mute()
            cpp_drv.ostream.mute()

            for i in range(natm):
                for d in range(3):
                    coords[i, d] += self.delta_h
                    new_mol = Molecule(labels, coords, units='au')
                    scf_drv.compute(new_mol, ao_basis)
                    cpp_drv._is_converged = False
                    cpp_results_p = cpp_drv.compute(new_mol, ao_basis,
                                                  scf_drv.scf_tensors)

                    coords[i, d] -= 2.0 * self.delta_h
                    new_mol = Molecule(labels, coords, units='au')
                    scf_drv.compute(new_mol, ao_basis)
                    cpp_drv._is_converged = False
                    cpp_results_m = cpp_drv.compute(new_mol, ao_basis,
                                                  scf_drv.scf_tensors)

                    coords[i, d] += self.delta_h
                    for f, w in enumerate(self.frequencies):
                        self.ostream.print_info('{}'.format(w))
                        self.ostream.flush()
                        for aop, acomp in enumerate('xyz'):
                            for bop, bcomp in enumerate('xyz'):
                                key = (acomp, bcomp, w)
                                if self.rank == mpi_master():
                                    num_polgradient[f, aop, bop, 3 * i + d] = (
                                    (cpp_results_p['response_functions'][key] -
                                     cpp_results_m['response_functions'][key]) /
                                    (2.0 * self.delta_h))
            if self.rank == mpi_master():
                for f, w in enumerate(self.frequencies):
                    polgrad_results[w] = num_polgradient[f]
                self.polgradient = dict(polgrad_results)
            else:
                self.polgradient = {}

        # four-point approximation for debugging of analytical gradient
        else:
            if self.rank == mpi_master():
                self.ostream.print_info('do_four_point: True')
                self.ostream.flush()
            # mute the outputstreams:
            scf_drv.ostream.mute()
            cpp_drv.ostream.mute()

            for i in range(natm):
                for d in range(3):
                    coords[i, d] += self.delta_h
                    new_mol = Molecule(labels, coords, units='au')
                    scf_drv.compute(new_mol, ao_basis)
                    cpp_drv._is_converged = False
                    cpp_results_p1 = cpp_drv.compute(new_mol, ao_basis,
                                                   scf_drv.scf_tensors)

                    coords[i, d] += self.delta_h
                    new_mol = Molecule(labels, coords, units='au')
                    scf_drv.compute(new_mol, ao_basis)
                    cpp_drv._is_converged = False
                    cpp_results_p2 = cpp_drv.compute(new_mol, ao_basis,
                                                   scf_drv.scf_tensors)

                    coords[i, d] -= 3.0 * self.delta_h
                    new_mol = Molecule(labels, coords, units='au')
                    scf_drv.compute(new_mol, ao_basis)
                    cpp_drv._is_converged = False
                    cpp_results_m1 = cpp_drv.compute(new_mol, ao_basis,
                                                   scf_drv.scf_tensors)

                    coords[i, d] -= 1.0 * self.delta_h
                    new_mol = Molecule(labels, coords, units='au')
                    scf_drv.compute(new_mol, ao_basis)
                    cpp_drv._is_converged = False
                    cpp_results_m2 = cpp_drv.compute(new_mol, ao_basis,
                                                   scf_drv.scf_tensors)

                    coords[i, d] += 2.0 * self.delta_h
                    for f, w in enumerate(self.frequencies):
                        for aop, acomp in enumerate('xyz'):
                            for bop, bcomp in enumerate('xyz'):
                                # f'(x) ~ [ f(x - 2h) - 8 f(x - h)
                                # + 8 f(x + h) - f(x + 2h) ] / ( 12h )
                                key = (acomp, bcomp, w)
                                if self.rank == mpi_master():
                                    num_polgradient[f, aop, bop, 3 * i + d] = ((
                                    cpp_results_m2['response_functions'][key] -
                                    8.0 * cpp_results_m1['response_functions'][key] +
                                    8.0 * cpp_results_p1['response_functions'][key] -
                                    cpp_results_p2['response_functions'][key]) / (
                                        12.0 * self.delta_h))
            if self.rank == mpi_master():
                for f, w in enumerate(self.frequencies):
                    polgrad_results[w] = num_polgradient[f]
                self.polgradient = dict(polgrad_results)
            else:
                self.polgradient = {}

        scf_drv.ostream.unmute()
        cpp_drv.ostream.unmute()

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
            # molgrid.broadcast(self.rank, self.comm) # TODO double check

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

    def print_polarizability_gradient(self, molecule):
        """
        Prints the polarizability gradient.

        :param molecule:
            The molecule.
        """

        labels = molecule.get_labels()

        natm = molecule.number_of_atoms()

        nfreq = len(self.frequencies)

        if self.numerical:
            title = 'Numerical '
        else:
            title = 'Analytical '

        # title for result output
        header = title + 'Polarizability Gradients'
        line = '=' * 40 

        self.ostream.print_blank()
        self.ostream.print_header(header)
        self.ostream.print_header(line)
        self.ostream.print_blank()
        self.ostream.flush()

        #gradient = self.polgradient
        for i in range(natm):
            atom_info = '** Atom #{0}: {1} **'.format(i+1, labels[i])
            self.ostream.print_header(atom_info)
            self.ostream.print_blank()

            # column headers
            column_headers = '{:<12}{:<16}{:>12}{:>12}{:>12} \n'.format(
                    "frequency", "alpha components", "coord. x", "coord. y", "coord. z")
            gradient_block =  column_headers

            for w in self.frequencies:
                freq = '{:<12}'.format(round(w, 4))
                current_gradient = self.polgradient[w].reshape(3, 3, natm, 3)
                for aop, acomp in enumerate('xyz'):
                    for bop, bcomp in enumerate('xyz'):
                        row = freq + '{:<16}'.format(acomp + bcomp)
                        for cop, ccomp in enumerate('xyz'):
                            grad_element = current_gradient[aop, bop, i, cop]
                            row += '{:>12}'.format(round(grad_element, 4))
                        row += ' \n'
                gradient_block += row

            self.ostream.print_block(gradient_block)
            self.ostream.print_blank()
            self.ostream.print_blank()

    def print_geometry(self, molecule):
        """
        Prints the geometry.

        :param molecule:
            The molecule.
        """

        self.ostream.print_block(molecule.get_string())
