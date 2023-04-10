import numpy as np
import time as tm

from .veloxchemlib import mpi_master
from .veloxchemlib import AODensityMatrix
from .veloxchemlib import AOFockMatrix
from .veloxchemlib import denmat
from .veloxchemlib import fockmat
from .veloxchemlib import XCIntegrator
from .tddftorbitalresponse import TddftOrbitalResponse
from .molecule import Molecule
from .firstorderprop import FirstOrderProperties
from .lrsolver import LinearResponseSolver
from .gradientdriver import GradientDriver
from .scfgradientdriver import ScfGradientDriver
from .errorhandler import assert_msg_critical
from .inputparser import parse_seq_fixed

# For PySCF integral derivatives
from .import_from_pyscf import overlap_deriv
from .import_from_pyscf import hcore_deriv
from .import_from_pyscf import eri_deriv

class TddftGradientDriver(GradientDriver):
    """
    Implements the analytic gradient driver for excited states at the
    Tamm-Dancoff approximation (TDA) and random phase approximation (RPA)
    level based on an SCF ground state (both HF and DFT).

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.

    Instance variables
        - scf_drv: The SCF driver.
        - gradient: The gradient.
        - tamm_dancoff: Flag if Tamm-Dancoff approximation is employed.
        - state_deriv_index: The excited state of interest.
        - do_first_order_prop: Controls the printout of first-order properties.
        - relaxed_dipole_moment: The relaxed excited-state dipole moment.
        - delta_h: The displacement for finite difference.
        - do_four_point: Flag for four-point finite difference.
    """

    def __init__(self, scf_drv, comm=None, ostream=None):
        """
        Initializes gradient driver.
        """

        super().__init__(scf_drv, comm, ostream)

        self.rank = self.comm.Get_rank()

        self.flag = 'RPA Gradient Driver'
        self.gradient = None

        # SCF (and response) drivers
        self.scf_drv = scf_drv

        # flag on whether RPA or TDA is calculated
        self.tamm_dancoff = False

        # excited state information;
        # if it is set to None, all available states
        # will be calculated.
        self.state_deriv_index = None

        # flag on whether to print excited-state properties
        self.do_first_order_prop = False
        self.relaxed_dipole_moment = None

        # for numerical gradient
        self.delta_h = 0.001

    def update_settings(self,
                        grad_dict,
                        rsp_dict,
                        orbrsp_dict=None,
                        method_dict=None):
        """
        Updates settings in gradient driver.

        :param grad_dict:
            The input dictionary of gradient settings group.
        :param orbrsp_dict:
            The input dictionary of orbital response settings group.
        :param rsp_dict:
            The input dictionary of response settings  group.
        :param method_dict:
            The input dicitonary of method settings group.
        """

        if method_dict is None:
            method_dict = {}

        if orbrsp_dict is None:
            orbrsp_dict = {}

        # basic settings from parent class
        super().update_settings(grad_dict, method_dict)

        if 'tamm_dancoff' in rsp_dict:
            key = rsp_dict['tamm_dancoff'].lower()
            self.tamm_dancoff = True if key in ['yes', 'y'] else False

        if self.tamm_dancoff:
            self.flag = 'TDA Gradient Driver'

        if 'do_first_order_prop' in grad_dict:
            key = grad_dict['do_first_order_prop'].lower()
            self.do_first_order_prop = True if key in ['yes', 'y'] else False
            orbrsp_dict['do_first_order_prop'] = (
                                        grad_dict['do_first_order_prop'] )

        # Excited states of interest
        # NOTE: this is a tuple;
        # the indexing starts at 1.
        if 'state_deriv_index' in grad_dict:
            self.state_deriv_index = parse_seq_fixed(
                    grad_dict['state_deriv_index'], flag='int')
            orbrsp_dict['state_deriv_index'] = grad_dict['state_deriv_index']

        self.grad_dict = dict(grad_dict)
        self.rsp_dict = dict(rsp_dict)
        self.method_dict = dict(method_dict)
        self.orbrsp_dict = dict(orbrsp_dict)

    def compute(self, molecule, basis, rsp_drv, rsp_results, min_basis=None):
        """
        Performs calculation of analytical or numerical gradient.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        :param rsp_drv:
            The RPA or TDA driver.
        :param rsp_results:
            The results from the RPA or TDA calculation.
        """

        # sanity check
        if self.rank == mpi_master():
            all_states = list(np.arange(1, len(rsp_results['eigenvalues'])+1))
            if self.state_deriv_index is not None:
                error_message =  'TdscfGradientDriver: some of the '
                error_message += 'selected states have not been calculated.'
                assert_msg_critical(
                    set(self.state_deriv_index).issubset(all_states),
                        error_message)
            else:
                self.state_deriv_index = all_states

            if self.numerical:
                self.print_header(self.state_deriv_index[:1])
            else:
                self.print_header(self.state_deriv_index)

        start_time = tm.time()

        # compute gradient
        # NOTE: the numerical gradient is calculated for
        # the first state only.
        if self.numerical:
            scf_ostream_state = self.scf_drv.ostream.state
            self.scf_drv.ostream.state = False
            self.compute_numerical(molecule, basis, rsp_drv, rsp_results,
                                   min_basis, self.state_deriv_index[0]-1)
            self.scf_drv.ostream.state = scf_ostream_state
        else:
            self.compute_analytical(molecule, basis, rsp_results)

        # print gradient
        if self.rank == mpi_master():
            self.print_geometry(molecule)
            if self.numerical:
                self.print_gradient(molecule, self.state_deriv_index[:1])
            else:
                self.print_gradient(molecule, self.state_deriv_index)

            valstr = '*** Time spent in gradient calculation: '
            valstr += '{:.2f} sec ***'.format(tm.time() - start_time)
            self.ostream.print_header(valstr)
            self.ostream.print_blank()
            self.ostream.flush()

    def compute_analytical(self, molecule, basis, rsp_results):
        """
        Performs calculation of analytical gradient.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        :param rsp_results:
            The results of the RPA or TDA calculation.
        """

        # select orbital response driver
        if self.tamm_dancoff:
            method = 'TDA'
        else:
            method = 'RPA'

        # SCF results
        scf_tensors = self.scf_drv.scf_tensors

        # compute orbital response
        orbrsp_drv = TddftOrbitalResponse(self.comm, self.ostream)
        orbrsp_drv.update_settings(self.orbrsp_dict, self.rsp_dict,
                                   self.method_dict)
        orbrsp_drv.compute(molecule, basis, scf_tensors,
                           rsp_results)

        orbrsp_results = orbrsp_drv.cphf_results
        omega_ao = orbrsp_drv.compute_omega(molecule, basis, scf_tensors)

        if self.rank == mpi_master():
            # only alpha part
            gs_dm = self.scf_drv.scf_tensors['D_alpha']
            nocc = molecule.number_of_alpha_electrons()
            natm = molecule.number_of_atoms()
            mo = scf_tensors['C']
            mo_occ = mo[:, :nocc]
            mo_vir = mo[:, nocc:]
            nocc = mo_occ.shape[1]
            nvir = mo_vir.shape[1]
            nao = mo_occ.shape[0]

            # TODO: check variable names and make sure they are consistent
            # with cphfsolver.
            # spin summation already included
            x_plus_y_ao = orbrsp_results['x_plus_y_ao']
            x_minus_y_ao = orbrsp_results['x_minus_y_ao']

            # CPHF/CPKS coefficients (lambda Lagrange multipliers)
            cphf_ov = orbrsp_results['cphf_ov']
            unrelaxed_density_ao = orbrsp_results['unrelaxed_density_ao']
            cphf_ao = np.einsum('mi,sia,na->smn', mo_occ, cphf_ov, mo_vir)
            relaxed_density_ao = ( unrelaxed_density_ao + 2.0 * cphf_ao
                        + 2.0 * cphf_ao.transpose(0,2,1) )
            dof = x_plus_y_ao.shape[0]
        else:
            dof = None

        dof = self.comm.bcast(dof, root=mpi_master())

        # ground state gradient
        gs_grad_drv = ScfGradientDriver(self.scf_drv)
        gs_grad_drv.update_settings(self.grad_dict, self.method_dict)

        ostream_state = gs_grad_drv.ostream.state
        gs_grad_drv.ostream.state = False
        gs_grad_drv.compute(molecule, basis)
        gs_grad_drv.ostream.state = ostream_state

        if self.rank == mpi_master():
            self.gradient = np.zeros((dof, natm, 3))
            self.gradient += gs_grad_drv.get_gradient()

            # loop over atoms and contract integral derivatives
            # with density matrices
            # add the corresponding contribution to the gradient
            for i in range(natm):
                # taking integral derivatives from pyscf
                d_ovlp = overlap_deriv(molecule, basis, i)
                d_hcore = hcore_deriv(molecule, basis, i)
                d_eri = eri_deriv(molecule, basis, i)

                self.gradient[:,i] += 1.0 * np.einsum('smn,xmn->sx',
                                                    relaxed_density_ao,
                                                    d_hcore)
                self.gradient[:,i] += 1.0 * np.einsum('smn,xmn->sx',
                                                       2.0 * omega_ao, d_ovlp)

                if self._dft:
                    if self.xcfun.is_hybrid():
                        frac_K = self.xcfun.get_frac_exact_exchange()
                    else:
                        frac_K = 0.0
                else:
                    frac_K = 1.0

                self.gradient[:,i] += 2.0 * np.einsum('mt,snp,xmtnp->sx', gs_dm,
                                                    relaxed_density_ao, d_eri)
                self.gradient[:,i] += -1.0 * frac_K * np.einsum(
                    'mt,snp,xmnpt->sx', gs_dm, relaxed_density_ao, d_eri)

                self.gradient[:,i] += 1.0 * np.einsum('smn,spt,xtpmn->sx',
                                    x_plus_y_ao,
                                    x_plus_y_ao - x_plus_y_ao.transpose(0,2,1),
                                    d_eri)
                self.gradient[:,i] += -0.5 * frac_K * np.einsum(
                                'smn,spt,xtnmp->sx', x_plus_y_ao,
                                x_plus_y_ao - x_plus_y_ao.transpose(0,2,1),
                                d_eri)

                self.gradient[:,i] += 1.0 * np.einsum('smn,spt,xtpmn->sx',
                                  x_minus_y_ao,
                                  x_minus_y_ao + x_minus_y_ao.transpose(0,2,1),
                                  d_eri)
                self.gradient[:,i] += -0.5 * frac_K * np.einsum(
                                'smn,spt,xtnmp->sx', x_minus_y_ao,
                                x_minus_y_ao + x_minus_y_ao.transpose(0,2,1),
                                d_eri)

        # TODO: enable multiple DMs for DFT to avoid for-loops.
        if self._dft:
            xcfun_label = self.scf_drv.xcfun.get_func_label()

            for s in range(dof):
                if self.rank == mpi_master():
                    gs_dm = self.scf_drv.scf_tensors['D_alpha']
                    gs_density = AODensityMatrix([gs_dm], denmat.rest)

                    rhow_dm = 0.5 * relaxed_density_ao[s]
                    rhow_dm_sym = 0.5 * (rhow_dm + rhow_dm.T)
                    rhow_den_sym = AODensityMatrix([rhow_dm_sym], denmat.rest)

                    x_minus_y_sym = 0.5 * (x_minus_y_ao[s] + x_minus_y_ao[s].T)
                    x_minus_y_den_sym = AODensityMatrix([x_minus_y_sym],
                                                        denmat.rest)
                else:
                    gs_density = AODensityMatrix()
                    rhow_den_sym = AODensityMatrix()
                    x_minus_y_den_sym = AODensityMatrix()

                gs_density.broadcast(self.rank, self.comm)
                rhow_den_sym.broadcast(self.rank, self.comm)
                x_minus_y_den_sym.broadcast(self.rank, self.comm)

                tddft_xcgrad = self.grad_tddft_xc_contrib(molecule, basis,
                                               rhow_den_sym, x_minus_y_den_sym,
                                               gs_density, xcfun_label)

                if self.rank == mpi_master():
                    self.gradient[s] += tddft_xcgrad

    def compute_energy(self,
                       molecule,
                       basis,
                       rsp_drv,
                       rsp_results=None,
                       min_basis=None,
                       state_deriv_index=0):
        """
        Computes the energy at the current position.

        :param molecule:
            The molecule.
        :param basis:
            The basis set.
        :param rsp_drv:
            The linear response driver.
        :param rsp_results:
            The results of a linear response calculation.
        :param min_basis:
            The reduced basis set.
        :param state_deriv_index:
            The index of the excited state of interest.
        """

        self.scf_drv.ostream.mute()
        self.scf_drv.compute(molecule, basis, min_basis)
        scf_tensors = self.scf_drv.scf_tensors

        rsp_drv._is_converged = False  # needed by RPA
        rsp_drv.ostream.mute()
        rsp_results = rsp_drv.compute(molecule, basis, scf_tensors)

        if self.rank == mpi_master():
            exc_en = rsp_results['eigenvalues'][state_deriv_index]
        else:
            exc_en = None
        exc_en = self.comm.bcast(exc_en, root=mpi_master())

        self.scf_drv.ostream.unmute()
        rsp_drv.ostream.unmute()

        return self.scf_drv.get_scf_energy() + exc_en

    def compute_numerical_dipole(self,
                                 molecule,
                                 ao_basis,
                                 rsp_drv,
                                 field_strength=1e-5,
                                 min_basis=None):
        """
        Performs calculation of numerical dipole moment at RPA or TDA level.

        :param molecule:
            The molecule.
        :param ao_basis:
            The AO basis set.
        :param rsp_drv:
            The RPA or TDA driver.
        :param field_strength:
            The strength of the external electric field.
        :param min_basis:
            The minimal AO basis set.
        :return:
            The electric dipole moment vector.
        """

        self.scf_drv.ostream.mute()

        # numerical dipole moment of the excited states
        n_states = len(self.state_deriv_index)
        dipole_moment = np.zeros((n_states, 3))
        field = [0.0, 0.0, 0.0]

        for s in range(n_states):
            for i in range(3):
                field[i] = field_strength
                self.scf_drv.electric_field = field
                self.scf_drv.compute(molecule, ao_basis, min_basis)
                scf_tensors = self.scf_drv.scf_tensors
                rsp_drv._is_converged = False  # only needed for RPA
                rsp_results = rsp_drv.compute(molecule, ao_basis, scf_tensors)
                exc_en_plus = rsp_results['eigenvalues'][s]
                e_plus = self.scf_drv.get_scf_energy() + exc_en_plus

                field[i] = -field_strength
                self.scf_drv.compute(molecule, ao_basis, min_basis)
                rsp_drv._is_converged = False
                rsp_results = rsp_drv.compute(molecule, ao_basis,
                                              self.scf_drv.scf_tensors)
                exc_en_minus = rsp_results['eigenvalues'][s]
                e_minus = self.scf_drv.get_scf_energy() + exc_en_minus

                field[i] = 0.0
                dipole_moment[s, i] = ( -(e_plus - e_minus)
                                         / (2.0 * field_strength) )

        self.scf_drv.ostream.unmute()

        return dipole_moment

    def compute_polarizability_grad(self, molecule, ao_basis, min_basis=None):
        """
        Performs calculation of numerical nuclear gradient
        of the electric dipole polarizability.

        :param molecule:
            The molecule.
        :param ao_basis:
            The AO basis set.
        """

        scf_ostream_state = self.scf_drv.ostream.state
        self.scf_drv.ostream.state = False

        # number of atoms
        natm = molecule.number_of_atoms()

        # atom labels
        labels = molecule.get_labels()

        # atom coordinates (nx3)
        coords = molecule.get_coordinates()

        # linear response driver for polarizability calculation
        lr_drv = LinearResponseSolver(self.comm, self.ostream)
        lr_drv.ostream.state = False

        # polarizability: 3 coordinates x 3 coordinates (ignoring frequencies)
        # polarizability gradient: dictionary goes through 3 coordinates
        # x 3 coordinates, each entry having values for
        # no. atoms x 3 coordinates
        self.pol_grad = np.zeros((3, 3, 3 * natm))

        if not self.do_four_point:
            for i in range(natm):
                for d in range(3):
                    coords[i, d] += self.delta_h
                    new_mol = Molecule(labels, coords, units='au')
                    self.scf_drv.compute(new_mol, ao_basis, min_basis)
                    lr_drv._is_converged = False
                    lr_results_p = lr_drv.compute(new_mol, ao_basis,
                                                  self.scf_drv.scf_tensors)

                    coords[i, d] -= 2.0 * self.delta_h
                    new_mol = Molecule(labels, coords, units='au')
                    self.scf_drv.compute(new_mol, ao_basis, min_basis)
                    lr_drv._is_converged = False
                    lr_results_m = lr_drv.compute(new_mol, ao_basis,
                                                  self.scf_drv.scf_tensors)

                    coords[i, d] += self.delta_h
                    for aop, acomp in enumerate('xyz'):
                        for bop, bcomp in enumerate('xyz'):
                            key = (acomp, bcomp, 0.0)
                            self.pol_grad[aop, bop, 3 * i + d] = (
                                (lr_results_p['response_functions'][key] -
                                 lr_results_m['response_functions'][key]) /
                                (2.0 * self.delta_h))

        # four-point approximation for debugging of analytical gradient
        else:
            for i in range(natm):
                for d in range(3):
                    coords[i, d] += self.delta_h
                    new_mol = Molecule(labels, coords, units='au')
                    self.scf_drv.compute(new_mol, ao_basis, min_basis)
                    lr_drv._is_converged = False
                    lr_results_p1 = lr_drv.compute(new_mol, ao_basis,
                                                   self.scf_drv.scf_tensors)

                    coords[i, d] += self.delta_h
                    new_mol = Molecule(labels, coords, units='au')
                    self.scf_drv.compute(new_mol, ao_basis, min_basis)
                    lr_drv._is_converged = False
                    lr_results_p2 = lr_drv.compute(new_mol, ao_basis,
                                                   self.scf_drv.scf_tensors)

                    coords[i, d] -= 3.0 * self.delta_h
                    new_mol = Molecule(labels, coords, units='au')
                    self.scf_drv.compute(new_mol, ao_basis, min_basis)
                    lr_drv._is_converged = False
                    lr_results_m1 = lr_drv.compute(new_mol, ao_basis,
                                                   self.scf_drv.scf_tensors)

                    coords[i, d] -= 1.0 * self.delta_h
                    new_mol = Molecule(labels, coords, units='au')
                    self.scf_drv.compute(new_mol, ao_basis, min_basis)
                    lr_drv._is_converged = False
                    lr_results_m2 = lr_drv.compute(new_mol, ao_basis,
                                                   self.scf_drv.scf_tensors)

                    coords[i, d] += 2.0 * self.delta_h
                    for aop, acomp in enumerate('xyz'):
                        for bop, bcomp in enumerate('xyz'):
                            # f'(x) ~ [ f(x - 2h) - 8 f(x - h)
                            # + 8 f(x + h) - f(x + 2h) ] / ( 12h )
                            key = (acomp, bcomp, 0.0)
                            self.pol_grad[aop, bop, 3 * i + d] = ((
                                lr_results_m2['response_functions'][key] -
                                8.0 * lr_results_m1['response_functions'][key] +
                                8.0 * lr_results_p1['response_functions'][key] -
                                lr_results_p2['response_functions'][key]) / (
                                    12.0 * self.delta_h))

        self.scf_drv.ostream.state = scf_ostream_state
