import numpy as np
import time as tm

from .veloxchemlib import denmat
from .veloxchemlib import mpi_master
from .veloxchemlib import ElectronRepulsionIntegralsDriver
from .veloxchemlib import AODensityMatrix
from .veloxchemlib import AOFockMatrix
from .veloxchemlib import denmat
from .veloxchemlib import fockmat
from .veloxchemlib import XCIntegrator, XCNewIntegrator
from .cphfsolver import CphfSolver
from .qqscheme import get_qq_scheme
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
        # TODO: include rsp_drv in energy_drv
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

    # TODO: should grad_dict and orbrsp_dict be unified?
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
            self.print_header(self.state_deriv_index)

        start_time = tm.time()

        # compute gradient
        # TODO: enable numerical gradient of multiple states
        # / do first selected state only?
        if self.numerical:
            scf_ostream_state = self.scf_drv.ostream.state
            self.scf_drv.ostream.state = False
            self.compute_numerical(molecule, basis, rsp_drv, rsp_results,
                                   min_basis)
            self.scf_drv.ostream.state = scf_ostream_state
        else:
            self.compute_analytical(molecule, basis, rsp_results)

        # print gradient

        if self.rank == mpi_master():
            self.print_geometry(molecule)
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

            # TODO: change names to x_plus_y_ao, etc.!
            # spin summation already included
            xpy = orbrsp_results['x_plus_y_ao']
            xmy = orbrsp_results['x_minus_y_ao']
            lambda_ov = orbrsp_results['cphf_ov']
            unrel_dm_ao = orbrsp_results['unrelaxed_density_ao']
            lambda_ao = np.einsum('mi,sia,na->smn', mo_occ, lambda_ov, mo_vir)
            rel_dm_ao = ( unrel_dm_ao + 2.0 * lambda_ao
                        + 2.0 * lambda_ao.transpose(0,2,1) )
            dof = xpy.shape[0]
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

                self.gradient[:,i] += 1.0 * np.einsum('smn,xmn->sx', rel_dm_ao,
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
                                                    rel_dm_ao, d_eri)
                self.gradient[:,i] += -1.0 * frac_K * np.einsum(
                    'mt,snp,xmnpt->sx', gs_dm, rel_dm_ao, d_eri)

                self.gradient[:,i] += 1.0 * np.einsum('smn,spt,xtpmn->sx', xpy,
                                             xpy - xpy.transpose(0,2,1), d_eri)
                self.gradient[:,i] += -0.5 * frac_K * np.einsum(
                    'smn,spt,xtnmp->sx', xpy, xpy - xpy.transpose(0,2,1), d_eri)

                self.gradient[:,i] += 1.0 * np.einsum('smn,spt,xtpmn->sx', xmy,
                                            xmy + xmy.transpose(0,2,1), d_eri)
                self.gradient[:,i] += -0.5 * frac_K * np.einsum(
                    'smn,spt,xtnmp->sx', xmy, xmy + xmy.transpose(0,2,1), d_eri)

        # TODO: use for-loop for now;
        # TODO: ask Xin if he cab enable multiple DMs for DFT.
        if self._dft:
            xcfun_label = self.scf_drv.xcfun.get_func_label()

            for s in range(dof):
                if self.rank == mpi_master():
                    gs_dm = self.scf_drv.scf_tensors['D_alpha']
                    gs_density = AODensityMatrix([gs_dm], denmat.rest)

                    rhow_dm = 0.5 * rel_dm_ao[s]
                    rhow_dm_sym = 0.5 * (rhow_dm + rhow_dm.T)
                    rhow_den_sym = AODensityMatrix([rhow_dm_sym], denmat.rest)

                    xmy_sym = 0.5 * (xmy[s] + xmy[s].T)
                    xmy_den_sym = AODensityMatrix([xmy_sym], denmat.rest)

                else:
                    gs_density = AODensityMatrix()
                    rhow_den_sym = AODensityMatrix()
                    xmy_den_sym = AODensityMatrix()

                gs_density.broadcast(self.rank, self.comm)
                rhow_den_sym.broadcast(self.rank, self.comm)
                xmy_den_sym.broadcast(self.rank, self.comm)

                tddft_xcgrad = self.grad_tddft_xc_contrib(molecule, basis,
                                                    rhow_den_sym, xmy_den_sym,
                                                    gs_density, xcfun_label)

                if self.rank == mpi_master():
                    self.gradient[s] += tddft_xcgrad


    def compute_energy(self,
                       molecule,
                       basis,
                       rsp_drv,
                       rsp_results=None,
                       min_basis=None):
        """
        Computes the energy at the current position.

        :param molecule:
            The molecule.
        :param basis:
            The basis set.
        :param rsp_drv:
            The linear response driver.
        """

        self.scf_drv.compute(molecule, basis, min_basis)
        scf_tensors = self.scf_drv.scf_tensors

        rsp_drv._is_converged = False  # needed by RPA
        rsp_results = rsp_drv.compute(molecule, basis, scf_tensors)

        if self.rank == mpi_master():
            exc_en = rsp_results['eigenvalues'][self.state_deriv_index]
        else:
            exc_en = None
        exc_en = self.comm.bcast(exc_en, root=mpi_master())

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

        scf_ostream_state = self.scf_drv.ostream.state
        self.scf_drv.ostream.state = False

        # numerical gradient
        dipole_moment = np.zeros((3))
        field = [0.0, 0.0, 0.0]

        for i in range(3):
            field[i] = field_strength
            self.scf_drv.electric_field = field
            self.scf_drv.compute(molecule, ao_basis, min_basis)
            scf_tensors = self.scf_drv.scf_tensors
            rsp_drv._is_converged = False  # only needed for RPA
            rsp_results = rsp_drv.compute(molecule, ao_basis, scf_tensors)
            exc_en_plus = rsp_results['eigenvalues'][self.state_deriv_index]
            e_plus = self.scf_drv.get_scf_energy() + exc_en_plus

            field[i] = -field_strength
            self.scf_drv.compute(molecule, ao_basis, min_basis)
            rsp_drv._is_converged = False
            rsp_results = rsp_drv.compute(molecule, ao_basis,
                                          self.scf_drv.scf_tensors)
            exc_en_minus = rsp_results['eigenvalues'][self.state_deriv_index]
            e_minus = self.scf_drv.get_scf_energy() + exc_en_minus

            field[i] = 0.0
            dipole_moment[i] = -(e_plus - e_minus) / (2.0 * field_strength)

        self.scf_drv.ostream.state = scf_ostream_state

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
        #lr_ostream_state = lr_drv.ostream.state
        lr_drv.ostream.state = False
        # lr_drv has self.a_components and self.b_components = 'xyz'
        # as well as self.frequencies = (0,)
        #lr_drv.update_settings(rsp_settings, method_settings)

        # polarizability: 3 coordinates x 3 coordinates (ignoring frequencies)
        # polarizability gradient: dictionary goes through 3 coordinates x 3 coordinates
        # each entry having values for no. atoms x 3 coordinates
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
                    #for aop in lr_drv.a_components:
                    for aop, acomp in enumerate('xyz'):
                        #for bop in lr_drv.b_components:
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
                    #for aop in lr_drv.a_components:
                    for aop, acomp in enumerate('xyz'):
                        #for bop in lr_drv.b_components:
                        for bop, bcomp in enumerate('xyz'):
                            # f'(x) ~ [ f(x - 2h) - 8 f(x - h) + 8 f(x + h) - f(x + 2h) ] / ( 12h )
                            key = (acomp, bcomp, 0.0)
                            self.pol_grad[aop, bop, 3 * i + d] = ((
                                lr_results_m2['response_functions'][key] -
                                8.0 * lr_results_m1['response_functions'][key] +
                                8.0 * lr_results_p1['response_functions'][key] -
                                lr_results_p2['response_functions'][key]) / (
                                    12.0 * self.delta_h))

        self.scf_drv.ostream.state = scf_ostream_state


class TddftOrbitalResponse(CphfSolver):
    """
    Implements orbital response Lagrange multipliers computation using a
    conjugate gradient scheme for the random phase approximation (RPA)
    level of theory.
    """

    def __init__(self, comm=None, ostream=None):
        """
        Initializes orbital response computation driver to default setup.

        :param comm:
            The MPI communicator.
        :param ostream:
            The output stream.
        """

        self.tamm_dancoff = False
        self.state_deriv_index = None
        self.do_first_order_prop = False

        super().__init__(comm, ostream)

    # TODO: are both orbrsp_dict and rsp_dict necessary?
    def update_settings(self, orbrsp_dict, rsp_dict, method_dict=None):
        """
        Updates response and method settings in orbital response computation
        driver.

        :param orbrsp_dict:
            The dictionary of orbital response settings.
        :param rsp_dict:
            The dictionary of response settings.
        :param method_dict:
            The dictionary of method settings.
        """

        # Excited states of interest
        # NOTE: this is a tuple;
        # the indexing starts at 1.
        if 'state_deriv_index' in orbrsp_dict:
            self.state_deriv_index = parse_seq_fixed(
                    orbrsp_dict['state_deriv_index'], flag='int')

        # Use TDA or not
        if 'tamm_dancoff' in rsp_dict:
            key = rsp_dict['tamm_dancoff'].lower()
            self.tamm_dancoff = True if key in ['yes', 'y'] else False

        # First Order Properties
        if 'do_first_order_prop' in orbrsp_dict:
            key = orbrsp_dict['do_first_order_prop'].lower()
            self.do_first_order_prop = True if key in ['yes', 'y'] else False

        super().update_settings(orbrsp_dict, method_dict)

    def compute(self, molecule, basis, scf_tensors, rsp_results):
        """ Computes the lambda orbital response multipliers and
            excited state first order properties.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        :param scf_tensors:
            The dictionary of tensors from converged SCF wavefunction.
        :param rsp_results:
            The dictionary of results from a converged linear response
            calculation.
        """

        super().compute(molecule, basis, scf_tensors, rsp_results)

        if self.do_first_order_prop:
            first_order_prop = FirstOrderProperties(self.comm, self.ostream)

            # unrelaxed density and dipole moment
            if self.rank == mpi_master():
                if self.tamm_dancoff:
                    method = 'TDA'
                else:
                    method = 'RPA'
                nocc = molecule.number_of_alpha_electrons()
                mo = scf_tensors['C']
                mo_occ = mo[:, :nocc]
                mo_vir = mo[:, nocc:]

                orbrsp_results = self.cphf_results

                lambda_ov = orbrsp_results['cphf_ov']
                unrel_dm_ao = orbrsp_results['unrelaxed_density_ao']
                lambda_ao = np.einsum('mi,sia,na->smn',
                        mo_occ, lambda_ov, mo_vir)
                rel_dm_ao = ( unrel_dm_ao + 2.0 * lambda_ao
                                + 2.0 * lambda_ao.transpose(0,2,1) )

                unrel_density = (scf_tensors['D'][0] + scf_tensors['D'][1] +
                                 unrel_dm_ao)
            else:
                unrel_density = None
            first_order_prop.compute(molecule, basis, unrel_density)

            if self.rank == mpi_master():
                title = method + ' Unrelaxed Dipole Moment(s) '
                first_order_prop.print_properties(molecule, title,
                                                  self.state_deriv_index)

            # relaxed density and dipole moment
            if self.rank == mpi_master():
                rel_density = (scf_tensors['D'][0] + scf_tensors['D'][1] +
                               rel_dm_ao)
            else:
                rel_density = None
            first_order_prop.compute(molecule, basis, rel_density)

            if self.rank == mpi_master():
                self.relaxed_dipole_moment = first_order_prop.get_property(
                        'dipole moment')

                title = method + ' Relaxed Dipole Moment(s) '
                first_order_prop.print_properties(molecule, title,
                                                  self.state_deriv_index)

                self.ostream.print_blank()

    def compute_rhs(self, molecule, basis, scf_tensors, rsp_results):
        """
        Computes the right-hand side (RHS) of the RPA orbital response equation
        including the necessary density matrices using molecular data.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        :param scf_tensors:
            The dictionary of tensors from the converged SCF calculation.
        :param rsp_results:
            The results from a converged linear response calculation.

        :return:
            A dictionary containing the orbital-response RHS and
            unrelaxed one-particle density.
        """

        # DFT information
        dft_dict = self._init_dft(molecule, scf_tensors)
        # ERI information
        eri_dict = self._init_eri(molecule, basis)
        # PE information
        pe_dict = self._init_pe(molecule, basis)

        self.profiler.start_timer('RHS')

        # Workflow:
        # 1) Construct the necessary density matrices
        # 2) Construct the RHS
        # 3) Construct the initial guess => in parent class
        # 4) Write the linear operator for matrix-vector product
        #    => in parent class
        # 5) Run the conjugate gradient => in parent class

        if self.rank == mpi_master():

            # 1) Calculate unrelaxed one-particle and transition density matrix
            ovlp = scf_tensors['S']
            mo = scf_tensors['C']

            nocc = molecule.number_of_alpha_electrons()
            mo_occ = mo[:, :nocc].copy()
            mo_vir = mo[:, nocc:].copy()
            nvir = mo_vir.shape[1]
            nao = mo.shape[0]

            # Take vectors of interest and convert to matrix form
            # n_states x nocc x nvir
            if self.state_deriv_index is None:
                # if no states are selected, calculate all
                # first excitd state is S1 (index starts at 1)
                self.state_deriv_index = list(np.arange(1, 
                                            len(rsp_results['eigenvalues'])+1))

            # number of degrees of freedon:
            dof = len(self.state_deriv_index)
            exc_vec = np.zeros((dof, nocc, nvir))
            deexc_vec = np.zeros((dof, nocc, nvir))

            for i in range(dof):
                ivec = self.state_deriv_index[i] - 1
                exc_vec[i] = (
                        rsp_results['eigenvectors'][:nocc * nvir,
                                                    ivec].reshape(nocc, nvir)
                                )

                if not self.tamm_dancoff:
                    deexc_vec[i] = (
                            rsp_results['eigenvectors'][nocc * nvir:,
                                                     ivec].reshape(nocc, nvir)
                                )

            # Construct plus/minus combinations of excitation
            # and de-excitation part
            xpy = exc_vec + deexc_vec
            xmy = exc_vec - deexc_vec

            # Transform the vectors to the AO basis
            xpy_ao = np.einsum('mi,sia,na->smn', mo_occ, xpy, mo_vir)
            xmy_ao = np.einsum('mi,sia,na->smn', mo_occ, xmy, mo_vir)

            # Calcuate the unrelaxed one-particle density matrix in MO basis
            dm_oo = -0.5 * (  np.einsum('sia,sja->sij', xpy, xpy)
                            + np.einsum('sia,sja->sij', xmy, xmy)
                            )
            dm_vv = 0.5 * (  np.einsum('sia,sib->sab', xpy, xpy)
                           + np.einsum('sia,sib->sab', xmy, xmy)
                            )

            # Transform unrelaxed one-particle density matrix to the AO basis
            unrel_dm_ao = ( np.einsum('mi,sij,nj->smn', mo_occ, dm_oo, mo_occ)
                          + np.einsum('ma,sab,nb->smn', mo_vir, dm_vv, mo_vir)
                            )

            # Make a list of unrel. DMs and excittion vectors:
            dm_ao_list = list(unrel_dm_ao) + list(xpy_ao) + list(xmy_ao)
  
            # 2) Construct the right-hand side
            dm_ao_rhs = AODensityMatrix(dm_ao_list, denmat.rest)

            if self._dft:
                # 3) Construct density matrices for E[3] term:
                # XCIntegrator expects a DM with real and imaginary part,
                # so we set the imaginary part to zero.
                perturbed_dm_ao_list = []
                zero_dm_ao_list = []

                # for each vector, we need to create a list with these elements:
                # xmy_ao[i], 0*xmy_ao[i], xmy_ao[i], 0*xmy_ao[i];
                # and a list with 2 list of zeros for each 4 elements above.

                for s in range(dof):
                    perturbed_dm_ao_list.extend([xmy_ao[s], 0*xmy_ao[s],
                                                 xmy_ao[s], 0*xmy_ao[s]])

                    zero_dm_ao_list.extend([0*xmy_ao[s], 0*xmy_ao[s]])

                perturbed_dm_ao = AODensityMatrix(perturbed_dm_ao_list,
                                                  denmat.rest)

                # corresponds to rho^{omega_b,omega_c} in quadratic response,
                # which is zero for TDDFT orbital response
                zero_dm_ao = AODensityMatrix(zero_dm_ao_list, denmat.rest)
        else:
            dof = None
            dm_ao_rhs = AODensityMatrix()
            if self._dft:
                perturbed_dm_ao = AODensityMatrix()
                zero_dm_ao =  AODensityMatrix()

        dof = self.comm.bcast(dof, root=mpi_master())
        dm_ao_rhs.broadcast(self.rank, self.comm)

        molgrid = dft_dict['molgrid']
        gs_density = dft_dict['gs_density']

        # Fock matrices with corresponding type
        fock_ao_rhs = AOFockMatrix(dm_ao_rhs)

       # Set the vector-related components to general Fock matrix
       # (not 1PDM part)
        for ifock in range(dof, 3*dof):
            fock_ao_rhs.set_fock_type(fockmat.rgenjk, ifock)

        if self._dft:
            perturbed_dm_ao.broadcast(self.rank, self.comm)
            zero_dm_ao.broadcast(self.rank, self.comm)
            # Fock matrix for computing gxc
            fock_gxc_ao = AOFockMatrix(zero_dm_ao)
            if self.xcfun.is_hybrid():
                fact_xc = self.xcfun.get_frac_exact_exchange()
                for ifock in range(fock_ao_rhs.number_of_fock_matrices()):
                    fock_ao_rhs.set_scale_factor(fact_xc, ifock)
                    fock_ao_rhs.set_fock_type(fockmat.restjkx, ifock)
                for ifock in range(fock_gxc_ao.number_of_fock_matrices()):
                    fock_gxc_ao.set_scale_factor(fact_xc, ifock)
                    fock_gxc_ao.set_fock_type(fockmat.rgenjkx, ifock)
                for ifock in range(dof):
                    fock_ao_rhs.set_fock_type(fockmat.restjkx, ifock)
                for ifock in range(dof, 3*dof):
                    fock_ao_rhs.set_fock_type(fockmat.rgenjkx, ifock)
            else:
                for ifock in range(dof):
                    fock_ao_rhs.set_fock_type(fockmat.restj, ifock)
                for ifock in range(dof, 3*dof):
                    fock_ao_rhs.set_fock_type(fockmat.rgenj, ifock)
                for ifock in range(fock_gxc_ao.number_of_fock_matrices()):
                    fock_gxc_ao.set_fock_type(fockmat.rgenj, ifock)
        else:
            fock_gxc_ao = None

        if self._dft:
            if not self.xcfun.is_hybrid():
                for ifock in range(fock_gxc_ao.number_of_fock_matrices()):
                    fock_gxc_ao.scale(2.0, ifock)
            # Quadratic response routine for TDDFT E[3] term g^xc
            # TODO: remove commented out code.
            # xc_drv = XCIntegrator(self.comm)
            # molgrid.distribute(self.rank, self.nodes, self.comm)
            # # Quadratic response routine for TDDFT E[3] term g^xc
            # xc_drv.integrate(fock_gxc_ao, perturbed_dm_ao, zero_dm_ao,
            #                  gs_density, molecule, basis, molgrid,
            #                  self.xcfun.get_func_label(), "qrf")
            xc_drv = XCNewIntegrator()
            molgrid.partition_grid_points()
            molgrid.distribute_counts_and_displacements(self.rank,
                                                self.nodes, self.comm)
            xc_drv.integrate_kxc_fock(fock_gxc_ao, molecule, basis, 
                                      perturbed_dm_ao, zero_dm_ao,
                                      gs_density, molgrid,
                                      self.xcfun.get_func_label(), "qrf")

            fock_gxc_ao.reduce_sum(self.rank, self.nodes, self.comm)

        self._comp_lr_fock(fock_ao_rhs, dm_ao_rhs, molecule, basis,
                          eri_dict, dft_dict, pe_dict, self.profiler)

        # Calculate the RHS and transform it to the MO basis
        if self.rank == mpi_master():
            # Extract the 1PDM contributions
            fock_ao_rhs_1pdm = np.zeros((dof, nao, nao))
            for ifock in range(dof):
                fock_ao_rhs_1pdm[ifock] = fock_ao_rhs.alpha_to_numpy(ifock)

            # Extract the excitation vector contributions
            fock_ao_rhs_xpy = np.zeros((dof, nao, nao))
            fock_ao_rhs_xmy = np.zeros((dof, nao, nao))
            for ifock in range(dof):
                fock_ao_rhs_xpy[ifock] = fock_ao_rhs.alpha_to_numpy(dof+ifock)
                fock_ao_rhs_xmy[ifock] = fock_ao_rhs.alpha_to_numpy(2*dof+ifock)

            # Transform to MO basis:
            fmo_rhs_1pdm = np.einsum('mi,smn,na->sia', mo_occ,
                                      0.5 * fock_ao_rhs_1pdm, mo_vir)

            # TODO: factor out overlap matrix?
            sdp_pds = 0.25 * (
                np.einsum('mn,snt,spt->smp', ovlp, xpy_ao, fock_ao_rhs_xpy)
              + np.einsum('mn,snt,spt->smp', ovlp, xmy_ao, fock_ao_rhs_xmy)
              - np.einsum('smn,smt,tp->snp', fock_ao_rhs_xpy, xpy_ao, ovlp)
              - np.einsum('smn,smt,tp->snp', fock_ao_rhs_xmy, xmy_ao, ovlp)
              - np.einsum('mn,snt,stp->smp', ovlp, xpy_ao, fock_ao_rhs_xpy)
              + np.einsum('mn,snt,stp->smp', ovlp, xmy_ao, fock_ao_rhs_xmy)
              + np.einsum('smn,snt,tp->smp', fock_ao_rhs_xpy, xpy_ao, ovlp)
              - np.einsum('smn,snt,tp->smp', fock_ao_rhs_xmy, xmy_ao, ovlp)
            )   

            rhs_mo = ( fmo_rhs_1pdm 
                     + np.einsum('mi,smn,na->sia', mo_occ, sdp_pds, mo_vir)
                        )

            # Add DFT E[3] contribution to the RHS:
            if self._dft:
                gxc_ao = np.zeros((dof, nao, nao))
                for ifock in range(dof):
                    gxc_ao[ifock] = fock_gxc_ao.alpha_to_numpy(2*ifock)
                gxc_mo = np.einsum('mi,smn,na->sia', mo_occ, gxc_ao, mo_vir)
                rhs_mo += 0.25 * gxc_mo

        self.profiler.stop_timer('RHS')

        if self.rank == mpi_master():
            return {
                'cphf_rhs': rhs_mo,
                'density_occ_occ': dm_oo,
                'density_vir_vir': dm_vv,
                'x_plus_y_ao': xpy_ao,
                'x_minus_y_ao': xmy_ao,
                'unrelaxed_density_ao': unrel_dm_ao,
                'fock_ao_rhs': fock_ao_rhs,
                'fock_gxc_ao': fock_gxc_ao, # None if not DFT
            }
        else:
            return {}

    def compute_omega(self, molecule, basis, scf_tensors):
        """
        Calculates the omega Lagrange multipliers for the overlap matrix.

        :param molecule:
            The molecule.
        :param basis.
            The basis set.
        :param scf_tensors.
            The scf tensors.

        :return:
            a numpy array containing the Lagrange multipliers in AO basis.
        """

        # ERI information
        eri_dict = self._init_eri(molecule, basis)
        # DFT information
        dft_dict = self._init_dft(molecule, scf_tensors)
        # PE information
        pe_dict = self._init_pe(molecule, basis)

        if self.rank == mpi_master():

            # Get overlap, MO coefficients from scf_tensors
            ovlp = scf_tensors['S']
            nocc = molecule.number_of_alpha_electrons()
            mo = scf_tensors['C']
            mo_occ = mo[:, :nocc]
            mo_vir = mo[:, nocc:]
            nocc = mo_occ.shape[1]
            nvir = mo_vir.shape[1]
            nao = mo_occ.shape[0]

            mo_energies = scf_tensors['E']
            eocc = mo_energies[:nocc]
            evir = mo_energies[nocc:]
            eo_diag = np.diag(eocc)
            ev_diag = np.diag(evir)

            xpy_ao = self.cphf_results['x_plus_y_ao']
            xmy_ao = self.cphf_results['x_minus_y_ao']
            dof = xpy_ao.shape[0]

            fock_ao_rhs = self.cphf_results['fock_ao_rhs']

            # The density matrix; only alpha block;
            # Only works for the restricted case.
            D_occ = np.matmul(mo_occ, mo_occ.T)
            D_vir = np.matmul(mo_vir, mo_vir.T)

            # Because the excitation vector is not symmetric,
            # we need both the matrix (for omega_OO, and probably VO)
            # and its transpose (for omega_VV and OV).
            # This comes from the transformation of the 2PDM contribution
            # from MO to AO basis
            fock_ao_rhs_xpy = np.zeros((dof, nao, nao))
            fock_ao_rhs_xmy = np.zeros((dof, nao, nao))
            for ifock in range(dof):
                fock_ao_rhs_xpy[ifock] = (
                        fock_ao_rhs.alpha_to_numpy(ifock+dof) )
                fock_ao_rhs_xmy[ifock] = (
                        fock_ao_rhs.alpha_to_numpy(ifock+2*dof) )

            Fp1_vv = 0.5 * np.einsum('smn,smt,pt->snp',
                                      fock_ao_rhs_xpy, xpy_ao, ovlp)
            Fm1_vv = 0.5 * np.einsum('smn,smt,pt->snp',
                                      fock_ao_rhs_xmy, xmy_ao, ovlp)
            Fp2_vv = 0.5 * np.einsum('smn,snt,pt->smp',
                                      fock_ao_rhs_xpy, xpy_ao, ovlp)
            Fm2_vv = 0.5 * np.einsum('smn,snt,pt->smp',
                                      fock_ao_rhs_xmy, xmy_ao, ovlp)
            Fp1_oo = 0.5 * np.einsum('smn,stn,pt->smp',
                                      fock_ao_rhs_xpy, xpy_ao, ovlp)
            Fm1_oo = 0.5 * np.einsum('smn,stn,pt->smp',
                                      fock_ao_rhs_xmy, xmy_ao, ovlp)
            Fp2_oo = 0.5 * np.einsum('snm,stn,pt->smp',
                                      fock_ao_rhs_xpy, xpy_ao, ovlp)
            Fm2_oo = 0.5 * np.einsum('snm,stn,pt->smp',
                                      fock_ao_rhs_xmy, xmy_ao, ovlp)

            # Construct fock_lambda (the lambda multipliers/cphf coefficients
            # contracted with the two-electron integrals)
            cphf_ov = self.cphf_results['cphf_ov']
            lambda_ao = np.einsum('mi,sia,na->smn', mo_occ, cphf_ov, mo_vir)
            lambda_ao_list = list([lambda_ao[s] for s in range(dof)])
            ao_density_lambda = AODensityMatrix(lambda_ao_list, denmat.rest)
        else:
            ao_density_lambda = AODensityMatrix()

        ao_density_lambda.broadcast(self.rank, self.comm)
        fock_lambda = AOFockMatrix(ao_density_lambda)

        self._comp_lr_fock(fock_lambda, ao_density_lambda, molecule, basis,
                           eri_dict, dft_dict, pe_dict, self.profiler)

        if self.rank == mpi_master():
            # Compute the contributions from the relaxed 1PDM
            # to the omega Lagrange multipliers:
            fock_ao_lambda_np = np.zeros((dof, nao, nao))
            fock_ao_rhs_1pdm = np.zeros((dof, nao, nao))
            for ifock in range(dof):
                fock_ao_lambda_np[ifock] = fock_lambda.alpha_to_numpy(ifock)
                fock_ao_rhs_1pdm[ifock] = fock_ao_rhs.alpha_to_numpy(ifock)

            fmat = (  fock_ao_lambda_np
                    + fock_ao_lambda_np.transpose(0,2,1)
                    + 0.5 * fock_ao_rhs_1pdm
                    )

            omega_1pdm_2pdm_contribs = 0.5 * (
                    np.einsum('mn,snt,tp->smp', D_vir,
                                Fp1_vv + Fm1_vv - Fp2_vv + Fm2_vv, D_vir)
                  + np.einsum('mn,snt,tp->smp', D_occ,
                               Fp1_vv + Fm1_vv - Fp2_vv + Fm2_vv, D_vir)
                  + np.einsum('mn,snt,tp->smp', D_occ, Fp1_vv + Fm1_vv - Fp2_vv
                              + Fm2_vv, D_vir).transpose(0,2,1)
                  + np.einsum('mn,snt,tp->smp', D_occ,
                               Fp1_oo + Fm1_oo - Fp2_oo + Fm2_oo, D_occ)
                  + 2.0 * np.einsum('mn,snt,tp->smp', D_occ, fmat, D_occ)
                                               )

            # Construct the energy-weighted one particle density matrix
            dm_oo = 0.5 * self.cphf_results['density_occ_occ']
            dm_vv = 0.5 * self.cphf_results['density_vir_vir']
            epsilon_dm_ao = np.einsum('mi,ii,sij,nj->smn', mo_occ,
                                       eo_diag, dm_oo, mo_occ)
            epsilon_dm_ao += np.einsum('ma,aa,sab,nb->smn', mo_vir,
                                        ev_diag, dm_vv, mo_vir)
            epsilon_lambda_ao = np.einsum('mi,ii,sia,na->smn', mo_occ,
                                           eo_diag, cphf_ov, mo_vir)
            epsilon_dm_ao += ( epsilon_lambda_ao # OV
                             + epsilon_lambda_ao.transpose(0,2,1) ) # VO

            omega = - epsilon_dm_ao - omega_1pdm_2pdm_contribs

            fock_gxc_ao = self.cphf_results['fock_gxc_ao']

            if fock_gxc_ao is not None:
                factor = -0.25
                fock_gxc_ao_np = np.zeros((dof, nao, nao))
                for ifock in range(dof):
                    fock_gxc_ao_np[ifock] = fock_gxc_ao.alpha_to_numpy(2*ifock)
                omega += factor * np.einsum('mn,snt,tp->smp', D_occ,
                                             fock_gxc_ao_np, D_occ)

            return omega
        else:
            return None
