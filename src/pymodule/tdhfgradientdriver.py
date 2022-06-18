from mpi4py import MPI
import numpy as np
import time as tm
import sys

from .veloxchemlib import ElectricDipoleIntegralsDriver
from .veloxchemlib import mpi_master, dipole_in_debye, denmat
from .molecule import Molecule
from .gradientdriver import GradientDriver
from .scfgradientdriver import ScfGradientDriver
from .outputstream import OutputStream
from .tdaorbitalresponse import TdaOrbitalResponse
from .rpaorbitalresponse import RpaOrbitalResponse
from .errorhandler import assert_msg_critical
from .firstorderprop import FirstOrderProperties
from .lrsolver import LinearResponseSolver
from .aodensitymatrix import AODensityMatrix

# For PySCF integral derivatives
from .import_from_pyscf import overlap_deriv
from .import_from_pyscf import fock_deriv
from .import_from_pyscf import hcore_deriv
from .import_from_pyscf import eri_deriv


class TdhfGradientDriver(GradientDriver):
    """
    Implements the analytic gradient driver for excited states at the
    Tamm-Dancoff approximation (TDA) and random phase approximation (RPA)
    level based on a Hartree-Fock ground state.
    DFT references will be implemented in the future.

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.

    Instance variables
        - scf_drv: The SCF driver.
        - rsp_drv: Wrapper object for RPA/TDA calculations and results (?)
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
        # self.rsp_drv = rsp_drv # rsp_prop object from main.py

        # flag on whether RPA or TDA is calculated
        self.tamm_dancoff = False

        # excited state information, default to first excited state
        self.state_deriv_index = 0

        # flag on whether to print excited-state properties
        self.do_first_order_prop = False
        self.relaxed_dipole_moment = None

        # for numerical gradient
        self.delta_h = 0.001

    def update_settings(self, grad_dict, rsp_dict, orbrsp_dict=None, method_dict=None):
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

        if 'state_deriv_index' in grad_dict:
            # user gives '1' for first excited state, but internal index is 0
            self.state_deriv_index = int(grad_dict['state_deriv_index']) - 1
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

        # sanity check for number of state
        if self.rank == mpi_master():
            assert_msg_critical(
                self.state_deriv_index < rsp_results['eigenvalues'].size,
                'TdhfGradientDriver: not enough states calculated')

            self.print_header(self.state_deriv_index)

        start_time = tm.time()

        if self.numerical:
            self.compute_numerical(molecule, basis, rsp_drv, rsp_results, min_basis)
        else:
            self.compute_analytical(molecule, basis, rsp_results)

        if self.rank == mpi_master():
            # print gradient
            self.print_geometry(molecule)
            self.print_gradient(molecule)

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
            orbrsp_drv = TdaOrbitalResponse(self.comm, self.ostream)
        else:
            method = 'RPA'
            orbrsp_drv = RpaOrbitalResponse(self.comm, self.ostream)

        # SCF results
        scf_tensors = self.scf_drv.scf_tensors

        # compute orbital response
        orbrsp_drv.update_settings(self.orbrsp_dict, self.rsp_dict,
                                   self.method_dict)
        orbrsp_results = orbrsp_drv.compute(molecule, basis, scf_tensors,
                                            rsp_results)

        # TODO: remove this member variable after TDDFT gradient is correct.
        self.orbrsp_results = orbrsp_results

        if self.rank == mpi_master():
            natm = molecule.number_of_atoms()
            # TODO: the following is available on the master node only
            # only alpha part:
            gs_dm = self.scf_drv.scf_tensors['D_alpha']
            omega_ao = orbrsp_results['omega_ao']
            # spin summation already included:
            xpy = orbrsp_results['x_plus_y_ao']
            xmy = orbrsp_results['x_minus_y_ao']
            rel_dm_ao = orbrsp_results['relaxed_density_ao']

            # ground state gradient
            gs_grad_drv = ScfGradientDriver(self.scf_drv)
            gs_grad_drv.update_settings(self.grad_dict, self.method_dict)
            gs_grad_drv.compute(molecule, basis)
            gs_gradient = gs_grad_drv.get_gradient()

            self.gradient = gs_gradient.copy()
            
            # loop over atoms and contract integral derivatives with density matrices
            # add the corresponding contribution to the gradient
            for i in range(natm):
                # taking integral derivatives from pyscf
                d_ovlp = overlap_deriv(molecule, basis, i)
                d_hcore = hcore_deriv(molecule, basis, i)
                d_eri = eri_deriv(molecule, basis, i)

                # TODO: delete commented out code, related to gradient using fock matrix derivative
                self.gradient[i] += ( #np.einsum('mn,xmn->x', 2.0 * gs_dm + rel_dm_ao, d_fock)
                                 ##+1.0 * np.einsum('mn,xmn->x', 2.0 * gs_dm + rel_dm_ao, d_hcore)
                                 +1.0 * np.einsum('mn,xmn->x', rel_dm_ao, d_hcore)
                                 +1.0 * np.einsum('mn,xmn->x', 2.0 * omega_ao, d_ovlp)
                                     )

                if self.dft:
                    if self.xcfun.is_hybrid():
                        fact_exact_exchange = self.xcfun.get_frac_exact_exchange()
                    else:
                        fact_exact_exchange = 0
                else:
                    fact_exact_exchange = 1.0

                self.gradient[i] += (
                                 # GS gradient calculated before using ScfGradientDriver
                                 # TODO: remove commented out code or go back to explicitly 
                                 # including GS here, after TDDFT gradient is correct.
                                 #+2.0 * np.einsum('mt,np,xmtnp->x', gs_dm, gs_dm, d_eri)
                                 #-1.0 * fact_exact_exchange * np.einsum('mt,np,xmnpt->x', gs_dm, gs_dm, d_eri)
                                 +2.0 * np.einsum('mt,np,xmtnp->x', gs_dm, rel_dm_ao, d_eri)
                                 -1.0 * fact_exact_exchange * np.einsum('mt,np,xmnpt->x', gs_dm, rel_dm_ao, d_eri)
                                 #-2.0 * np.einsum('mt,np,xmtnp->x', gs_dm, gs_dm, d_eri)
                                 #+1.0 * np.einsum('mt,np,xmnpt->x', gs_dm, gs_dm, d_eri)
                                 +1.0 * np.einsum('mn,pt,xtpmn->x', xpy, xpy - xpy.T, d_eri)
                                 -0.5 * fact_exact_exchange * np.einsum('mn,pt,xtnmp->x', xpy, xpy - xpy.T, d_eri)
                                 +1.0 * np.einsum('mn,pt,xtpmn->x', xmy, xmy + xmy.T, d_eri)
                                 -0.5 * fact_exact_exchange * np.einsum('mn,pt,xtnmp->x', xmy, xmy + xmy.T, d_eri)
                                )

        if self.dft:
            xcfun_label = self.scf_drv.xcfun.get_func_label()

            gs_dm = self.scf_drv.scf_tensors['D_alpha']
            gs_density = AODensityMatrix([gs_dm], denmat.rest)

            rhow_dm = 0.5 * orbrsp_results['relaxed_density_ao']
            rhow_dm_sym = 0.5 * (rhow_dm + rhow_dm.T)
            rhow_den_sym = AODensityMatrix([rhow_dm_sym], denmat.rest)

            vxc_contrib = self.grad_vxc_contrib(molecule, basis, rhow_den_sym, gs_density, xcfun_label)
            vxc_contrib_2 = self.grad_fxc_contrib(molecule, basis, rhow_den_sym, gs_density, gs_density, xcfun_label)

            xpy_sym = 0.5 * (xpy + xpy.T)
            xpy_den_sym = AODensityMatrix([xpy_sym], denmat.rest)

            fxc_contrib = self.grad_fxc_contrib(molecule, basis, xpy_den_sym, xpy_den_sym, gs_density, xcfun_label)
            fxc_contrib_2 = self.grad_gxc_contrib(molecule, basis, xpy, xpy, gs_density, xcfun_label)

            if self.rank == mpi_master():
                self.gradient += vxc_contrib
                self.gradient += vxc_contrib_2
                self.gradient += fxc_contrib
                self.gradient += fxc_contrib_2

        # Calculate the relaxed and unrelaxed excited-state dipole moment
        firstorderprop = FirstOrderProperties(self.comm, self.ostream)

        if self.rank == mpi_master():
            # unrelaxed density and dipole moment
            unrel_density = (scf_tensors['D'][0] + scf_tensors['D'][1] +
                             orbrsp_results['unrelaxed_density_ao'])
        else:
            unrel_density = None

        firstorderprop.compute(molecule, basis, unrel_density)

        if self.rank == mpi_master():
            if self.do_first_order_prop:
                title = method + ' Unrelaxed Dipole Moment for Excited State ' + str(self.state_deriv_index + 1)
                firstorderprop.print_properties(molecule, title)

            # relaxed density and dipole moment
            if 'relaxed_density_ao' in orbrsp_results:
                rel_density = (scf_tensors['D'][0] + scf_tensors['D'][1] +
                               orbrsp_results['relaxed_density_ao'])
                # TODO: consider computing the relaxed dipole moment in parallel?
                firstorderprop.compute(molecule, basis, rel_density)
                self.relaxed_dipole_moment = firstorderprop.get_property('dipole moment')

            if self.rank == mpi_master():
                if self.do_first_order_prop:
                    title = method + ' Relaxed Dipole Moment for Excited State ' + str(self.state_deriv_index + 1)
                    firstorderprop.print_properties(molecule, title)

            self.ostream.print_blank()

    def compute_energy(self, molecule, basis, rsp_drv, rsp_results=None, min_basis=None):
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

        rsp_drv.is_converged = False  # needed by RPA
        rsp_results = rsp_drv.compute(molecule, basis, scf_tensors)

        if self.rank == mpi_master():
            exc_en = rsp_results['eigenvalues'][self.state_deriv_index]
        else:
            exc_en = None
        exc_en = self.comm.bcast(exc_en, root=mpi_master())
    
        return self.scf_drv.get_scf_energy() + exc_en

    def compute_numerical_dipole(self, molecule, ao_basis, rsp_drv,
                                 field_strength=1e-5, min_basis=None):
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

        start_time = tm.time()

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
            rsp_drv.is_converged = False  # only needed for RPA
            rsp_results = rsp_drv.compute(molecule, ao_basis,
                                               scf_tensors)
            exc_en_plus = rsp_results['eigenvalues'][self.state_deriv_index]
            e_plus = self.scf_drv.get_scf_energy() + exc_en_plus

            field[i] = -field_strength
            self.scf_drv.compute(molecule, ao_basis, min_basis)
            rsp_drv.is_converged = False
            rsp_results = rsp_drv.compute(molecule, ao_basis,
                                               self.scf_drv.scf_tensors)
            exc_en_minus = rsp_results['eigenvalues'][self.state_deriv_index]
            e_minus = self.scf_drv.get_scf_energy() + exc_en_minus

            field[i] = 0.0
            dipole_moment[i] = - (e_plus - e_minus) / (2.0 * field_strength)

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
        # dictionary to translate from numbers to operator components 'xyz'
        component_dict = {0: 'x', 1: 'y', 2: 'z'}

        if not self.do_four_point:
            for i in range(natm):
                for d in range(3):
                    coords[i, d] += self.delta_h
                    new_mol = Molecule(labels, coords, units='au')
                    self.scf_drv.compute(new_mol, ao_basis, min_basis)
                    lr_drv.is_converged = False
                    lr_results_p = lr_drv.compute(new_mol, ao_basis,
                                                       self.scf_drv.scf_tensors)

                    coords[i, d] -= 2.0 * self.delta_h
                    new_mol = Molecule(labels, coords, units='au')
                    self.scf_drv.compute(new_mol, ao_basis, min_basis)
                    lr_drv.is_converged = False
                    lr_results_m = lr_drv.compute(new_mol, ao_basis,
                                                       self.scf_drv.scf_tensors)

                    coords[i, d] += self.delta_h
                    #for aop in lr_drv.a_components:
                    for aop in range(3):
                        #for bop in lr_drv.b_components:
                        for bop in range(3):
                            self.pol_grad[aop, bop, 3*i + d] = (
                                ( lr_results_p['response_functions'][component_dict[aop], component_dict[bop], 0.0]
                                - lr_results_m['response_functions'][component_dict[aop], component_dict[bop], 0.0] ) /
                                (2.0 * self.delta_h) )

        # four-point approximation for debugging of analytical gradient
        else:
            for i in range(natm):
                for d in range(3):
                    coords[i, d] += self.delta_h
                    new_mol = Molecule(labels, coords, units='au')
                    self.scf_drv.compute(new_mol, ao_basis, min_basis)
                    lr_drv.is_converged = False
                    lr_results_p1 = lr_drv.compute(new_mol, ao_basis,
                                                       self.scf_drv.scf_tensors)

                    coords[i, d] += self.delta_h
                    new_mol = Molecule(labels, coords, units='au')
                    self.scf_drv.compute(new_mol, ao_basis, min_basis)
                    lr_drv.is_converged = False
                    lr_results_p2 = lr_drv.compute(new_mol, ao_basis,
                                                       self.scf_drv.scf_tensors)

                    coords[i, d] -= 3.0 * self.delta_h
                    new_mol = Molecule(labels, coords, units='au')
                    self.scf_drv.compute(new_mol, ao_basis, min_basis)
                    lr_drv.is_converged = False
                    lr_results_m1 = lr_drv.compute(new_mol, ao_basis,
                                                       self.scf_drv.scf_tensors)

                    coords[i, d] -= 1.0 * self.delta_h
                    new_mol = Molecule(labels, coords, units='au')
                    self.scf_drv.compute(new_mol, ao_basis, min_basis)
                    lr_drv.is_converged = False
                    lr_results_m2 = lr_drv.compute(new_mol, ao_basis,
                                                       self.scf_drv.scf_tensors)

                    coords[i, d] += 2.0 * self.delta_h
                    #for aop in lr_drv.a_components:
                    for aop in range(3):
                        #for bop in lr_drv.b_components:
                        for bop in range(3):
                    # f'(x) ~ [ f(x - 2h) - 8 f(x - h) + 8 f(x + h) - f(x + 2h) ] / ( 12h )
                            self.pol_grad[aop, bop, 3*i + d] = (
                                ( lr_results_m2['response_functions'][component_dict[aop], component_dict[bop], 0.0]
                                - 8 * lr_results_m1['response_functions'][component_dict[aop], component_dict[bop], 0.0]
                                + 8 * lr_results_p1['response_functions'][component_dict[aop], component_dict[bop], 0.0]
                                - lr_results_p2['response_functions'][component_dict[aop], component_dict[bop], 0.0] ) /
                                (12.0 * self.delta_h) )
