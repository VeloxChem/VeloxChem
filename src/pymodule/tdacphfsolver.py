import numpy as np

from .veloxchemlib import ElectronRepulsionIntegralsDriver
from .veloxchemlib import AODensityMatrix
from .veloxchemlib import AOFockMatrix
from .veloxchemlib import mpi_master
from .veloxchemlib import denmat
from .veloxchemlib import fockmat
from .veloxchemlib import XCIntegrator
from .cphfsolver import CphfSolver
from .qqscheme import get_qq_scheme
from .inputparser import parse_seq_fixed

# TODO: move TdaCphfSolver in TdhfGradient?
class TdaCphfSolver(CphfSolver):
    """
    Implements orbital response Lagrange multipliers computation using a
    conjugate gradient scheme for the Tamm-Dancoff Approximation (TDA)
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

        self.tamm_dancoff = True
        self.state_deriv_index = None

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
    
        super().update_settings(orbrsp_dict, method_dict)

    def compute_rhs(self, molecule, basis, scf_tensors, tda_results):
        """
        Computes the right-hand side (RHS) of the TDA orbital response equation
        including the necessary density matrices using molecular data.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        :param scf_tensors:
            The dictionary of tensors from the converged SCF calculation.
        :param tda_results:
            The results from the converged TDA calculation.

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
                                            len(tda_results['eigenvalues'])+1))

            # number of degrees of freedon:
            dof = len(self.state_deriv_index)
            exc_vec = np.zeros((dof, nocc, nvir))

            for i in range(dof):
                ivec = self.state_deriv_index[i] - 1
                exc_vec[i] = (
                        tda_results['eigenvectors'][:, ivec].reshape(nocc, nvir)
                                )

            # Transform the excitation vectors to the AO basis
            exc_vec_ao = np.einsum("mi,sia,na->smn", mo_occ, exc_vec, mo_vir)

            # Calcuate the unrelaxed one-particle density matrix in MO basis
            dm_oo = - np.einsum('sia,sja->sij', exc_vec, exc_vec) 
            dm_vv =   np.einsum('sia,sib->sab', exc_vec, exc_vec)

            # Transform unrelaxed one-particle density matrix to the AO basis
            unrel_dm_ao = ( np.einsum('mi,sij,nj->smn', mo_occ, dm_oo, mo_occ)
                          + np.einsum('ma,sab,nb->smn', mo_vir, dm_vv, mo_vir)
                            )

            # Make a list of unrel. DMs and excittion vectors:
            dm_ao_list = list(unrel_dm_ao) + list(exc_vec_ao)
  
            # 2) Construct the right-hand side
            dm_ao_rhs = AODensityMatrix(dm_ao_list, denmat.rest)

            if self._dft:
                # 3) Construct density matrices for E[3] term:
                # XCIntegrator expects a DM with real and imaginary part,
                # so we set the imaginary part to zero.
                perturbed_dm_ao_list = []
                zero_dm_ao_list = []

                # for each vector, we need to create a list with these elements:
                # exc_vec_ao[i], 0*exc_vec_ao[i], ex_vec_ao[i], 0*exc_vec_ao[i];
                # and a list with 2 list of zeros for each 4 elements above.

                for s in range(dof):
                    perturbed_dm_ao_list.extend([exc_vec_ao[s], 0*exc_vec_ao[s],
                                                 exc_vec_ao[s], 0*exc_vec_ao[s]])

                    zero_dm_ao_list.extend([0*exc_vec_ao[s], 0*exc_vec_ao[s]])

                perturbed_dm_ao = AODensityMatrix(perturbed_dm_ao_list,
                                                  denmat.rest)

                # corresponds to rho^{omega_b,omega_c} in quadratic response,
                # which is zero for TDDFT orbital response
                zero_dm_ao = AODensityMatrix(zero_dm_ao_list, denmat.rest)
        else:
            dm_ao_rhs = AODensityMatrix()
            if self._dft:
                perturbed_dm_ao = AODensityMatrix()
                zero_dm_ao =  AODensityMatrix()

        dm_ao_rhs.broadcast(self.rank, self.comm)

        molgrid = dft_dict['molgrid']
        gs_density = dft_dict['gs_density']

        # Fock matrices with corresponding type
        fock_ao_rhs = AOFockMatrix(dm_ao_rhs)

       # Set the vector-related components to general Fock matrix
       # (not 1PDM part)
        for ifock in range(dof, 2*dof):
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
                for ifock in range(dof, 2*dof):
                    fock_ao_rhs.set_fock_type(fockmat.rgenjkx, ifock)
            else:
                for ifock in range(dof):
                    fock_ao_rhs.set_fock_type(fockmat.restj, ifock)
                for ifock in range(dof, 2*dof):
                    fock_ao_rhs.set_fock_type(fockmat.rgenj, ifock)
                for ifock in range(fock_gxc_ao.number_of_fock_matrices()):
                    fock_gxc_ao.set_fock_type(fockmat.rgenj, ifock)
        else:
            fock_gxc_ao = None

        if self._dft:
            if not self.xcfun.is_hybrid():
                for ifock in range(fock_gxc_ao.number_of_fock_matrices()):
                    fock_gxc_ao.scale(2.0, ifock)
            xc_drv = XCIntegrator(self.comm)
            molgrid.distribute(self.rank, self.nodes, self.comm)
            # Quadratic response routine for TDDFT E[3] term g^xc
            xc_drv.integrate(fock_gxc_ao, perturbed_dm_ao, zero_dm_ao,
                             gs_density, molecule, basis, molgrid,
                             self.xcfun.get_func_label(), "qrf") #"quadratic")

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
            fock_ao_rhs_exc_vec = np.zeros((dof, nao, nao))
            for ifock in range(dof):
                fock_ao_rhs_exc_vec[ifock] = fock_ao_rhs.alpha_to_numpy(
                                                                    dof+ifock)
            # Transform to MO basis:
            fmo_rhs_1pdm = np.einsum('mi,smn,na->sia', mo_occ,
                                      0.5 * fock_ao_rhs_1pdm, mo_vir)
            sdp_pds = (
                  0.5 * np.einsum('mn,snt,slt->sml', ovlp,
                                 exc_vec_ao, fock_ao_rhs_exc_vec)
                - 0.5 * np.einsum('smn,smt,tl->snl', fock_ao_rhs_exc_vec,
                                 exc_vec_ao, ovlp)
                    )

            rhs_mo = ( fmo_rhs_1pdm 
                     + np.einsum('mi,smn,na->sia', mo_occ, sdp_pds, mo_vir)
                        )

            # Add DFT E[3] contribution to the RHS:
            if self._dft:
                gxc_ao = np.zeros((dof,nao,nao))
                for ifock in range(dof):
                    gxc_ao[ifock] = fock_gxc_ao.alpha_to_numpy(2*ifock)
                gxc_mo = np.einsum('mi,smn,na->sia', mo_occ, gxc_ao, mo_vir)
                rhs_mo += 0.25 * gxc_mo

        self.profiler.stop_timer('RHS')

        if self.rank == mpi_master():
            return {
                'cphf_rhs': rhs_mo,
                'dm_oo': dm_oo,
                'dm_vv': dm_vv,
                'xpy_ao': exc_vec_ao,
                'xmy_ao': exc_vec_ao,
                'unrel_dm_ao': unrel_dm_ao,
                'fock_ao_rhs': fock_ao_rhs,
                'fock_gxc_ao': fock_gxc_ao, # None if not DFT
            }
        else:
            return {}

    def compute_omega(self, molecule, basis, scf_tensors):
        """
        Calculates the TDA Lagrange multipliers for the overlap matrix.

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

            exc_vec_ao = self.cphf_results['xpy_ao']
            # The number of degrees of freedom (dof) is the number of
            # excited states we want to calculate.
            dof = exc_vec_ao.shape[0]

            fock_ao_rhs = self.cphf_results['fock_ao_rhs']

            # The density matrix; only alpha block;
            # Only works for the restricted case
            D_occ = np.matmul(mo_occ, mo_occ.T)
            D_vir = np.matmul(mo_vir, mo_vir.T)

            # Because the excitation vector is not symmetric,
            # we need both the matrix (for omega_OO, and probably VO)
            # and its transpose (for omega_VV and OV).
            # This comes from the transformation of the 2PDM contribution
            # from MO to AO basis
            fock_ao_rhs_exc_vec = np.zeros((dof, nao, nao))
            for ifock in range(dof):
                fock_ao_rhs_exc_vec[ifock] = (
                        fock_ao_rhs.alpha_to_numpy(ifock+dof) )

            # TODO: can Ft be replaced by F.transpose(0,2,1)?
            Ft = 0.5 * np.einsum('smn,smt,tl->snl', fock_ao_rhs_exc_vec, 
                                  exc_vec_ao, ovlp)
            F = 0.5 * np.einsum('smn,stn,lt->sml', fock_ao_rhs_exc_vec,
                                 exc_vec_ao, ovlp)

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

            omega_1pdm_2pdm_contribs = (
                   np.einsum('mn,snt,tp->smp', D_occ, F, D_occ)
                 + np.einsum('mn,snt,tp->smp', D_occ, Ft, D_vir)
                 + np.einsum('mn,snt,tp->smp', D_occ, Ft, D_vir).transpose(0,2,1)
                 + np.einsum('mn,snt,tp->smp', D_vir, Ft, D_vir)
                 + np.einsum('mn,snt,tp->smp', D_occ, fmat, D_occ)
                                        )

            # Construct the energy-weighted one particle density matrix
            dm_oo = 0.5 * self.cphf_results['dm_oo']
            dm_vv = 0.5 * self.cphf_results['dm_vv']
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
