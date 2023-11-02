import numpy as np
import time as tm

from .veloxchemlib import mpi_master
from .veloxchemlib import AODensityMatrix
from .veloxchemlib import AOFockMatrix
from .veloxchemlib import denmat
from .veloxchemlib import fockmat
from .veloxchemlib import XCIntegrator
from .cphfsolver import CphfSolver
from .firstorderprop import FirstOrderProperties
from .inputparser import parse_seq_fixed
from .visualizationdriver import VisualizationDriver

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

        super().__init__(comm, ostream)

        self.tamm_dancoff = False
        self.state_deriv_index = None
        self.do_first_order_prop = False

        self._input_keywords['orbitalresponse'].update({
            'tamm_dancoff': ('bool', 'whether RPA or TDA is calculated'),
            'state_deriv_index': ('seq_fixed_int', 'excited state information'),
            'do_first_order_prop': ('bool', 'do first-order property'),
            }
        )

    # TODO: are both orbrsp_dict and rsp_dict necessary?
    # NOTE: not if we use tamm_dancoff keyword in orbrsp_dict
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

        super().update_settings(orbrsp_dict, method_dict)

        orbrsp_keywords = {
            key: val[0] for key, val in self._input_keywords['orbitalresponse'].items()
        }

        parse_input(self, orbrsp_keywords, orbrsp_dict)

        orbrsp_dict['tamm_dancoff'] = rsp_dict['tamm_dancoff']

        # Excited states of interest
        # NOTE: this is a tuple; the indexing starts at 1.
        #if 'state_deriv_index' in orbrsp_dict:
        #    self.state_deriv_index = parse_seq_fixed(
        #            orbrsp_dict['state_deriv_index'], flag='int')

        ## Use TDA or not
        #if 'tamm_dancoff' in rsp_dict:
        #    key = rsp_dict['tamm_dancoff'].lower()
        #    self.tamm_dancoff = True if key in ['yes', 'y'] else False

        ## First Order Properties
        #if 'do_first_order_prop' in orbrsp_dict:
        #    key = orbrsp_dict['do_first_order_prop'].lower()
        #    self.do_first_order_prop = True if key in ['yes', 'y'] else False

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

    @staticmethod
    def get_full_solution_vector(solution):
        """
        Gets a full solution vector from the distributed solution.

        :param solution:
            The distributed solution as a tuple.

        :return:
            The full solution vector.
        """

        x_ger = solution.get_full_vector(0)
        x_ung = solution.get_full_vector(1)

        if solution.rank == mpi_master():
            x_ger_full = np.hstack((x_ger, x_ger))
            x_ung_full = np.hstack((x_ung, -x_ung))
            return x_ger_full + x_ung_full
        else:
            return None

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
                if self.tamm_dancoff:
                    exc_vec[i] = (
                        rsp_results['eigenvectors'][:nocc * nvir,
                                                    ivec].reshape(nocc, nvir)
                                )
                else:
                    eigenvector = self.get_full_solution_vector(
                                rsp_results['eigenvectors_distributed'][ivec])
                    exc_vec[i] = eigenvector[:nocc * nvir].reshape(nocc, nvir)
                    deexc_vec[i] = eigenvector[nocc * nvir:].reshape(nocc, nvir)

            # Construct plus/minus combinations of excitation
            # and de-excitation part
            x_plus_y = exc_vec + deexc_vec
            x_minus_y = exc_vec - deexc_vec

            # Transform the vectors to the AO basis
            x_plus_y_ao = np.zeros((dof, nao, nao))
            x_minus_y_ao = np.zeros((dof, nao, nao))
            dm_oo = np.zeros((dof, nocc, nocc))
            dm_vv = np.zeros((dof, nvir, nvir))
            unrel_dm_ao = np.zeros((dof, nao, nao))
            for i in range(dof):
                x_plus_y_ao[i] = np.linalg.multi_dot([mo_occ, x_plus_y[i],
                                                      mo_vir.T])
                x_minus_y_ao[i] = np.linalg.multi_dot([mo_occ, x_minus_y[i],
                                                       mo_vir.T])
                dm_oo[i] = -0.5 * ( np.linalg.multi_dot([x_plus_y[i],
                                                         x_plus_y[i].T])
                                  + np.linalg.multi_dot([x_minus_y[i],
                                                         x_minus_y[i].T])
                                    )
                dm_vv[i] = 0.5 * ( np.linalg.multi_dot([x_plus_y[i].T,
                                                        x_plus_y[i]])
                                 + np.linalg.multi_dot([x_minus_y[i].T,
                                                        x_minus_y[i]])
                                    )
                unrel_dm_ao[i] = ( np.linalg.multi_dot([mo_occ, dm_oo[i],
                                                        mo_occ.T])
                                 + np.linalg.multi_dot([mo_vir, dm_vv[i],
                                                        mo_vir.T])
                                    )

            # Make a list of unrel. DMs and excittion vectors:
            dm_ao_list = ( list(unrel_dm_ao) + list(x_plus_y_ao)
                         + list(x_minus_y_ao) )
  
            # 2) Construct the right-hand side
            dm_ao_rhs = AODensityMatrix(dm_ao_list, denmat.rest)

            if self._dft:
                # 3) Construct density matrices for E[3] term:
                # XCIntegrator expects a DM with real and imaginary part,
                # so we set the imaginary part to zero.
                perturbed_dm_ao_list = []
                zero_dm_ao_list = []

                # for each vector, we need to create a list with these elements:
                # x_minus_y_ao[i], 0*x_minus_y_ao[i],
                # x_minus_y_ao[i], 0*x_minus_y_ao[i];
                # and a list with 2 list of zeros for each 4 elements above.

                for s in range(dof):
                    perturbed_dm_ao_list.extend([x_minus_y_ao[s],
                                                 0 * x_minus_y_ao[s],
                                                 x_minus_y_ao[s],
                                                 0 * x_minus_y_ao[s]])

                    zero_dm_ao_list.extend([0 * x_minus_y_ao[s],
                                            0 * x_minus_y_ao[s]])

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
            xc_drv = XCIntegrator()
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
            fock_ao_rhs_x_plus_y = np.zeros((dof, nao, nao))
            fock_ao_rhs_x_minus_y = np.zeros((dof, nao, nao))
            for ifock in range(dof):
                fock_ao_rhs_x_plus_y[ifock] = fock_ao_rhs.alpha_to_numpy(dof
                                                                         + ifock)
                fock_ao_rhs_x_minus_y[ifock] = fock_ao_rhs.alpha_to_numpy(2 * dof
                                                                        + ifock)

            # Transform to MO basis:
            rhs_mo = np.zeros((dof, nocc, nvir))
            for i in range(dof):
                fmo_rhs_1pdm = 0.5 * np.linalg.multi_dot([mo_occ.T,
                                                          fock_ao_rhs_1pdm[i], 
                                                          mo_vir])

                sdp_pds = 0.25 * (
                    np.linalg.multi_dot([ovlp, x_plus_y_ao[i],
                                         fock_ao_rhs_x_plus_y[i].T])
                  + np.linalg.multi_dot([ovlp, x_minus_y_ao[i],
                                         fock_ao_rhs_x_minus_y[i].T])
                  - np.linalg.multi_dot([fock_ao_rhs_x_plus_y[i].T,
                                         x_plus_y_ao[i], ovlp])
                  - np.linalg.multi_dot([fock_ao_rhs_x_minus_y[i].T,
                                         x_minus_y_ao[i], ovlp])
                  - np.linalg.multi_dot([ovlp, x_plus_y_ao[i],
                                         fock_ao_rhs_x_plus_y[i]])
                  + np.linalg.multi_dot([ovlp, x_minus_y_ao[i],
                                         fock_ao_rhs_x_minus_y[i]])
                  + np.linalg.multi_dot([fock_ao_rhs_x_plus_y[i],
                                         x_plus_y_ao[i], ovlp])
                  - np.linalg.multi_dot([fock_ao_rhs_x_minus_y[i],
                                         x_minus_y_ao[i], ovlp])
                                  )

                rhs_mo[i] = ( fmo_rhs_1pdm 
                     + np.linalg.multi_dot([mo_occ.T, sdp_pds, mo_vir])
                        )

            # Add DFT E[3] contribution to the RHS:
            if self._dft:
                gxc_ao = np.zeros((dof, nao, nao))
                gxc_mo = np.zeros((dof, nocc, nvir))
                for ifock in range(dof):
                    gxc_ao[ifock] = fock_gxc_ao.alpha_to_numpy(2*ifock)
                    gxc_mo[ifock] = np.linalg.multi_dot([mo_occ.T, gxc_ao[ifock],
                                                         mo_vir])
                rhs_mo += 0.25 * gxc_mo

        self.profiler.stop_timer('RHS')

        if self.rank == mpi_master():
            return {
                'cphf_rhs': rhs_mo,
                'density_occ_occ': dm_oo,
                'density_vir_vir': dm_vv,
                'x_plus_y_ao': x_plus_y_ao,
                'x_minus_y_ao': x_minus_y_ao,
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

            x_plus_y_ao = self.cphf_results['x_plus_y_ao']
            x_minus_y_ao = self.cphf_results['x_minus_y_ao']
            dof = x_plus_y_ao.shape[0]

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
            fock_ao_rhs_x_plus_y = np.zeros((dof, nao, nao))
            fock_ao_rhs_x_minus_y = np.zeros((dof, nao, nao))
            for ifock in range(dof):
                fock_ao_rhs_x_plus_y[ifock] = (
                        fock_ao_rhs.alpha_to_numpy(ifock+dof) )
                fock_ao_rhs_x_minus_y[ifock] = (
                        fock_ao_rhs.alpha_to_numpy(ifock+2*dof) )

            Fp1_vv = 0.5 * np.einsum('smn,smt,pt->snp',
                                      fock_ao_rhs_x_plus_y, x_plus_y_ao, ovlp)
            Fm1_vv = 0.5 * np.einsum('smn,smt,pt->snp',
                                      fock_ao_rhs_x_minus_y, x_minus_y_ao, ovlp)
            Fp2_vv = 0.5 * np.einsum('smn,snt,pt->smp',
                                      fock_ao_rhs_x_plus_y, x_plus_y_ao, ovlp)
            Fm2_vv = 0.5 * np.einsum('smn,snt,pt->smp',
                                      fock_ao_rhs_x_minus_y, x_minus_y_ao, ovlp)
            Fp1_oo = 0.5 * np.einsum('smn,stn,pt->smp',
                                      fock_ao_rhs_x_plus_y, x_plus_y_ao, ovlp)
            Fm1_oo = 0.5 * np.einsum('smn,stn,pt->smp',
                                      fock_ao_rhs_x_minus_y, x_minus_y_ao, ovlp)
            Fp2_oo = 0.5 * np.einsum('snm,stn,pt->smp',
                                      fock_ao_rhs_x_plus_y, x_plus_y_ao, ovlp)
            Fm2_oo = 0.5 * np.einsum('snm,stn,pt->smp',
                                      fock_ao_rhs_x_minus_y, x_minus_y_ao, ovlp)

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

