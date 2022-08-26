import numpy as np
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
from .veloxchemlib import XCFunctional
from .veloxchemlib import MolecularGrid
from .veloxchemlib import parse_xc_func
from .veloxchemlib import GridDriver, XCMolecularGradient
#from .orbitalresponse import OrbitalResponse
from .cphfsolver import CphfSolver
from .outputstream import OutputStream
from .qqscheme import get_qq_scheme
from .errorhandler import assert_msg_critical

# For PySCF integral derivatives
from .import_from_pyscf import overlap_deriv
from .import_from_pyscf import hcore_deriv
#from .import_from_pyscf import fock_deriv
from .import_from_pyscf import eri_deriv
from .import_from_pyscf import dipole_deriv

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

        self.gradient = None
        self.flag = "Polarizability Gradient"
        # self.numerical = False
        # # flag for two-point or four-point approximation
        # self.do_four_point = False

        # DFT information
        self.dft = False
        self.grid_level = 4
        self.xcfun = XCFunctional()

        # Polarizability information
        self.frequency = 0.0
        self.vector_components = 'xyz'

    def update_settings(self, grad_dict, orbrsp_dict=None, method_dict=None):
        """
        Updates response and method settings in orbital response computation
        driver.

        :param grad_dict:
            The input dictionary of gradient settings group.
        :param orbrsp_dict:
            The dictionary of orbital response settings.
        :param method_dict:
            The dictionary of method settings.
        """

        if method_dict is None:
            method_dict = {}
        if orbrsp_dict is None:
            orbrsp_dict = {}

        if 'dft' in method_dict:
            key = method_dict['dft'].lower()
            self.dft = (key in ['yes', 'y'])
        if 'grid_level' in method_dict:
            self.grid_level = int(method_dict['grid_level'])
        if 'xcfun' in method_dict:
            if 'dft' not in method_dict:
                self.dft = True
            self.xcfun = parse_xc_func(method_dict['xcfun'].upper())
            assert_msg_critical(not self.xcfun.is_undefined(),
                                'PolarizabilityGradient: Undefined XC functional')

        if 'frequency' in grad_dict:
            self.frequency = float(grad_dict['frequency'])
            print("Updating frequency...: ", self.frequency)
            if 'frequency' not in orbrsp_dict:
                orbrsp_dict['frequency'] = grad_dict['frequency']

        if 'vector_components' in grad_dict:
            self.vector_components = grad_dict['vector_components']
            if 'vector_components' not in orbrsp_dict:
                orbrsp_dict['vector_components'] = grad_dict['vector_components']

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
        """

        # orbital response driver
        orbrsp_drv = PolOrbitalResponse(self.comm, self.ostream)

        # compute orbital response
        orbrsp_drv.update_settings(self.orbrsp_dict, self.method_dict)
        orbrsp_drv.compute(molecule, basis, scf_tensors, lr_results)
        orbrsp_drv.compute_omega(molecule, basis, scf_tensors, lr_results)
        orbrsp_results = orbrsp_drv.cphf_results

        if self.rank == mpi_master():
            mo = scf_tensors['C'] # only alpha part
            ovlp = scf_tensors['S']
            nao = mo.shape[0]
            nocc = molecule.number_of_alpha_electrons()
            mo_occ = mo[:, :nocc].copy()
            mo_vir = mo[:, nocc:].copy()
            nvir = mo_vir.shape[1]
            natm = molecule.number_of_atoms()

            # TODO: the following is available on the master node only
            # only alpha part:
            gs_dm = scf_tensors['D_alpha']
            lambda_mo = orbrsp_results['cphf_ov']
            # spin summation already included:
            xpy = orbrsp_results['xpy_ao']
            xmy = orbrsp_results['xmy_ao']
            dof = xpy.shape[0] # Number of vector components
            omega_ao = orbrsp_results['omega_ao'].reshape(dof, dof, nao, nao)


            # transform lambda multipliers to AO basis and calculate relaxed density matrix
            lambda_ao = np.einsum('mi,xia,na->xmn', mo_occ, lambda_mo,
                                  mo_vir).reshape(dof, dof, nao, nao) # occ-vir
            lambda_ao += lambda_ao.transpose(0,1,3,2) # vir-occ
            rel_dm_ao = orbrsp_results['unrel_dm_ao'] + lambda_ao # TODO: undo comment

            # analytical polarizability gradient
            pol_gradient = np.zeros((dof, dof, natm, 3))
            # dictionary to translate from numbers to operator components 'xyz'
            component_dict = {0: 'x', 1: 'y', 2: 'z'}

            if self.dft:
                if self.xcfun.is_hybrid():
                    frac_K = self.xcfun.get_frac_exact_exchange()
                else:
                    frac_K = 0.0
            else:
                frac_K = 1.0

            # TODO: remove this member variable once the freq-dependent gradient is correct.
            freq_dep_contrib = np.zeros((dof, dof, natm, 3))

            # loop over atoms and contract integral derivatives with density matrices
            # add the corresponding contribution to the gradient
            for i in range(natm):
                # taking integral derivatives from pyscf
                d_ovlp = overlap_deriv(molecule, basis, i)
                #d_fock = fock_deriv(molecule, basis, gs_dm, i)
                d_hcore = hcore_deriv(molecule, basis, i)
                d_eri = eri_deriv(molecule, basis, i)
                d_dipole = dipole_deriv(molecule, basis, i)

                # TODO: remove freq_dep_contrib
                freq_dep_contrib[:, :, i] += -0.5 * self.frequency * ( np.einsum('amk,nl,xmn,ykl->xya', d_ovlp, ovlp, xpy, xmy)
                                                         + np.einsum('amk,nl,xmn,ykl->xya', d_ovlp, ovlp, xmy, xpy)
                                                         + np.einsum('mk,anl,xmn,ykl->xya', ovlp, d_ovlp, xpy, xmy)
                                                         + np.einsum('mk,anl,xmn,ykl->xya', ovlp, d_ovlp, xmy, xpy)
                                                          )
                # Calculate the analytic polarizability gradient
                pol_gradient[:, :, i] += ( np.einsum('xymn,amn->xya', 2.0 * rel_dm_ao, d_hcore)
                                 +1.0 * np.einsum('xymn,amn->xya', 2.0 * omega_ao, d_ovlp)
                                 #-0.5 * self.frequency * ( np.einsum('amk,nl,xmn,ykl->xya', d_ovlp, ovlp, xpy, xmy)
                                 #                        + np.einsum('amk,nl,xmn,ykl->xya', d_ovlp, ovlp, xmy, xpy)
                                 #                        + np.einsum('mk,anl,xmn,ykl->xya', ovlp, d_ovlp, xpy, xmy)
                                 #                        + np.einsum('mk,anl,xmn,ykl->xya', ovlp, d_ovlp, xmy, xpy)
                                 #                         )
                                 +1.0 * (
                                 +2.0 * np.einsum('mt,xynp,amtnp->xya', gs_dm, 2.0 * rel_dm_ao, d_eri)
                                 -1.0 * frac_K * np.einsum('mt,xynp,amnpt->xya', gs_dm, 2.0 * rel_dm_ao, d_eri)
                                 +1.0 * np.einsum('xmn,ypt,atpmn->xya', xpy, xpy - xpy.transpose(0,2,1), d_eri)
                                 -0.5 * frac_K * np.einsum('xmn,ypt,atnmp->xya', xpy, xpy - xpy.transpose(0,2,1), d_eri)
                                 +1.0 * np.einsum('xmn,ypt,atpmn->xya', xmy, xmy + xmy.transpose(0,2,1), d_eri)
                                 -0.5 * frac_K * np.einsum('xmn,ypt,atnmp->xya', xmy, xmy + xmy.transpose(0,2,1), d_eri)
                                 +1.0 * np.einsum('xmn,ypt,atpmn->yxa', xpy, xpy - xpy.transpose(0,2,1), d_eri)
                                 -0.5 * frac_K * np.einsum('xmn,ypt,atnmp->yxa', xpy, xpy - xpy.transpose(0,2,1), d_eri)
                                 +1.0 * np.einsum('xmn,ypt,atpmn->yxa', xmy, xmy + xmy.transpose(0,2,1), d_eri)
                                 -0.5 * frac_K * np.einsum('xmn,ypt,atnmp->yxa', xmy, xmy + xmy.transpose(0,2,1), d_eri)
                                 )
                                 -2.0 * np.einsum('xmn,yamn->xya', xmy, d_dipole)
                                 -2.0 * np.einsum('xmn,yamn->yxa', xmy, d_dipole)
                                )

            # Add exchange-correlation contributions to the gradient
            # for now by looping over each component of the polarizability
            if self.dft:
                xcfun_label = self.xcfun.get_func_label()
                for i in range(dof):
                    #for j in range(dof): # include as soon as non-diagonal terms are available
                    if self.rank == mpi_master():
                        gs_density = AODensityMatrix([gs_dm], denmat.rest)

                        rhow_dm = 1.0 * rel_dm_ao[i,i]
                        rhow_dm_sym = 0.5 * (rhow_dm + rhow_dm.T)
                        rhow_den_sym = AODensityMatrix([rhow_dm_sym], denmat.rest)

                        # Takes only one vector type, but two are needed to account for
                        # the different {x,y,z} components of the response vectors
                        # doing only diagonal components for now.
                        # The sqrt2 takes into account the fact that we need to
                        # symmetrize with respect to the polarizability components.
                        # (see contraction with two-electron integrals above).
                        xmy_sym = np.sqrt(2) * 0.5 * (xmy[i] + xmy[i].T)
                        xmy_den_sym = AODensityMatrix([xmy_sym], denmat.rest)

                    else:
                        gs_density = AODensityMatrix()
                        rhow_den_sym = AODensityMatrix()
                        xmy_den_sym = AODensityMatrix()

                    gs_density.broadcast(self.rank, self.comm)
                    rhow_den_sym.broadcast(self.rank, self.comm)
                    xmy_den_sym.broadcast(self.rank, self.comm)

                    polgrad_xcgrad = self.grad_polgrad_xc_contrib(molecule, basis,
                                                                  rhow_den_sym, xmy_den_sym,
                                                                  gs_density, xcfun_label)

                    if self.rank == mpi_master():
                        pol_gradient[i,i] += polgrad_xcgrad



            self.pol_gradient = pol_gradient.reshape(dof, dof, 3 * natm)
            self.freq_dep_contrib = freq_dep_contrib.reshape(dof, dof, 3 * natm)
            self.orbrsp_results = dict(orbrsp_results)

    def grad_polgrad_xc_contrib(self, molecule, ao_basis, rhow_den, xmy_den,
                              gs_density, xcfun_label):
        """
        Calculates exchange-correlation contribution to polarizability gradient.

        :param molecule:
            The molecule.
        :param ao_basis:
            The AO basis set.
        :param rhow_den:
            The perturbed density.
        :param xmy_den:
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
        mol_grid.distribute(self.rank, self.nodes, self.comm)

        xcgrad_drv = XCMolecularGradient(self.comm)
        polgrad_xcgrad = xcgrad_drv.integrate_tddft_gradient(
            rhow_den, xmy_den, gs_density, molecule, ao_basis, mol_grid,
            xcfun_label)
        polgrad_xcgrad = self.comm.reduce(polgrad_xcgrad, root=mpi_master())

        return polgrad_xcgrad

    def init_dft(self, molecule, scf_tensors):
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
        if self.dft:
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
            molgrid.broadcast(self.rank, self.comm) # TODO duble check

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


class PolOrbitalResponse(CphfSolver):
    """
    Implements orbital response Lagrange multipliers computation
    for the polarizability gradient.

    Instance variables
        - frequency: The frequency for which the polarizability
            gradient is computed.
        - vector_components: The components of the response vectors
            corresponding to the operator components in linear response.
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

        self.frequency = 0.0
        self.vector_components = 'xyz'
        self.cphf_results = None

    def update_settings(self, orbrsp_dict, method_dict=None):
        """
        Updates response and method settings in orbital response computation
        driver.

        :param orbrsp_dict:
            The dictionary of orbital response settings.
        :param method_dict:
            The dictionary of method settings.
        """

        super().update_settings(orbrsp_dict, method_dict)

        if 'frequency' in orbrsp_dict:
            self.frequency = float(orbrsp_dict['frequency'])

        if 'vector_components' in orbrsp_dict:
            self.vector_components = orbrsp_dict['vector_components']

    def compute_rhs(self, molecule, basis, scf_tensors, lr_results):
        """
        Computes the right-hand side (RHS) of the polarizability orbital response equation
        including the necessary density matrices using molecular data.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        :param scf_tensors:
            The dictionary of tensors from converged SCF calculation.
        :param lr_results:
            The results from converged linear response calculation.

        :return:
            A dictionary containing the orbital-response RHS and
            unrelaxed one-particle density.
        """

        # DFT information
        dft_dict = self.init_dft(molecule, scf_tensors)

        self.profiler.start_timer('RHS')

        # Workflow:
        # 1) Construct the necessary density matrices
        # 2) Construct the RHS
        # 3) Construct the initial guess => in parent class
        # 4) Run the solver => in parent class

        if self.rank == mpi_master():

            # 1) Calculate unrelaxed one-particle and transition density matrix
            ovlp = scf_tensors['S']
            mo = scf_tensors['C'] # only alpha part

            nao = mo.shape[0]
            nocc = molecule.number_of_alpha_electrons()
            mo_occ = mo[:, :nocc].copy()
            mo_vir = mo[:, nocc:].copy()
            nvir = mo_vir.shape[1]

            # TODO: do we keep this factor like that?
            sqrt2 = np.sqrt(2.0)

            # Check if response vectors exist for desired frequency of gradient
            if (self.vector_components[0], self.frequency) not in lr_results['solutions'].keys():
                raise ValueError("Frequency for gradient not found in linear response results.")

            # Take response vectors and convert to matrix form
            exc_vec = 1/sqrt2 * np.array([lr_results['solutions'][x,
                                self.frequency][:nocc * nvir].reshape(nocc, nvir)
                                for x in self.vector_components])
            deexc_vec = 1/sqrt2 * np.array([lr_results['solutions'][x,
                                self.frequency][nocc * nvir:].reshape(nocc, nvir)
                                for x in self.vector_components])

            # Number of vector components
            dof = exc_vec.shape[0]

            # Construct plus/minus combinations of excitation and de-excitation part
            xpy = exc_vec + deexc_vec
            xmy = exc_vec - deexc_vec

            # Transform the vectors to the AO basis
            xpy_ao = np.einsum('mi,xia,na->xmn', mo_occ, xpy, mo_vir)
            xmy_ao = np.einsum('mi,xia,na->xmn', mo_occ, xmy, mo_vir)

            # Turn them into a list (for AODensityMatrix)
            xpmy_ao_list = list(xpy_ao) + list(xmy_ao)

            # Calculate the symmetrized unrelaxed one-particle density matrix in MO basis
            dm_oo = -0.25 * ( np.einsum('xja,yia->xyij', xpy, xpy) + np.einsum('xja,yia->xyij', xmy, xmy)
                             +np.einsum('yja,xia->xyij', xpy, xpy) + np.einsum('yja,xia->xyij', xmy, xmy)
                            )

            dm_vv = 0.25 * ( np.einsum('xib,yia->xyab', xpy, xpy) + np.einsum('xib,yia->xyab', xmy, xmy)
                            +np.einsum('yib,xia->xyab', xpy, xpy) + np.einsum('yib,xia->xyab', xmy, xmy)
                           )

            # Transform unrelaxed one-particle density matrix to the AO basis and to list
            unrel_dm_ao = np.einsum('mi,xyij,nj->xymn', mo_occ, dm_oo, mo_occ) + np.einsum('ma,xyab,nb->xymn', mo_vir, dm_vv, mo_vir)
            dm_ao_list = list(unrel_dm_ao.reshape(dof**2, nao, nao))

            # 2) Construct the right-hand side
            dm_ao_rhs = AODensityMatrix(dm_ao_list + xpmy_ao_list, denmat.rest)

            if self.dft:
                # 3) Construct density matrices for E[3] term:
                # XCIntegrator expects a DM with real and imaginary part,
                # so we set the imaginary part to zero.
                # Create lists with the corresponding vector components
                perturbed_dm_ao_list = []
                zero_dm_ao_list = []
                for x in range(dof):
                    for y in range(dof): # TODO: only upper triangular matrix and transpose?
                        perturbed_dm_ao_list.extend([xmy_ao[x], 0*xmy_ao[x], xmy_ao[y], 0*xmy_ao[y]])
                        zero_dm_ao_list.extend([0*xmy_ao[x], 0*xmy_ao[y]])

                perturbed_dm_ao = AODensityMatrix(perturbed_dm_ao_list, denmat.rest)

                # corresponds to rho^{omega_b,omega_c} in quadratic response,
                # which is zero for orbital response
                zero_dm_ao = AODensityMatrix(zero_dm_ao_list, denmat.rest)
        else:
            dm_ao_rhs = AODensityMatrix()
            if self.dft:
                perturbed_dm_ao = AODensityMatrix()
                zero_dm_ao =  AODensityMatrix()

        dm_ao_rhs.broadcast(self.rank, self.comm)

        molgrid = dft_dict['molgrid']
        gs_density = dft_dict['gs_density']

        # Fock matrices with corresponding type
        fock_ao_rhs = AOFockMatrix(dm_ao_rhs)
        # Set the vector-related components to general Fock matrix (not 1PDM part)
        for ifock in range(dof**2, dof**2 + 2 * dof):
            fock_ao_rhs.set_fock_type(fockmat.rgenjk, ifock)
        if self.dft:
            perturbed_dm_ao.broadcast(self.rank, self.comm)
            zero_dm_ao.broadcast(self.rank, self.comm)
            # Fock matrix for computing the DFT E[3] term g^xc
            fock_gxc_ao = AOFockMatrix(zero_dm_ao)
            if self.xcfun.is_hybrid():
                fact_xc = self.xcfun.get_frac_exact_exchange()
                for ifock in range(fock_ao_rhs.number_of_fock_matrices()):
                    fock_ao_rhs.set_scale_factor(fact_xc, ifock)
                for ifock in range(fock_gxc_ao.number_of_fock_matrices()):
                    fock_gxc_ao.set_scale_factor(fact_xc, ifock)
                    fock_gxc_ao.set_fock_type(fockmat.rgenjkx, ifock)
                for ifock in range(dof**2):
                    fock_ao_rhs.set_fock_type(fockmat.restjkx, ifock)
                for ifock in range(dof**2, dof**2 + 2 * dof):
                    fock_ao_rhs.set_fock_type(fockmat.rgenjkx, ifock)
            else:
                for ifock in range(dof**2):
                    fock_ao_rhs.set_fock_type(fockmat.restj, ifock)
                for ifock in range(dof**2, dof**2 + 2 * dof):
                    fock_ao_rhs.set_fock_type(fockmat.rgenj, ifock)
                for ifock in range(fock_gxc_ao.number_of_fock_matrices()):
                    fock_gxc_ao.set_fock_type(fockmat.rgenj, ifock)
        else:
            fock_gxc_ao = None # None if not DFT

        eri_drv = ElectronRepulsionIntegralsDriver(self.comm)
        screening = eri_drv.compute(get_qq_scheme(self.qq_type),
                                    self.eri_thresh, molecule, basis)

        # TODO: replace by comp_lr_fock?
        eri_drv.compute(fock_ao_rhs, dm_ao_rhs, molecule, basis, screening)

        if self.dft:
            if not self.xcfun.is_hybrid():
                for ifock in range(fock_ao_rhs.number_of_fock_matrices()):
                    fock_ao_rhs.scale(2.0, ifock)
                for ifock in range(fock_gxc_ao.number_of_fock_matrices()):
                    fock_gxc_ao.scale(2.0, ifock)
            xc_drv = XCIntegrator(self.comm)
            molgrid.distribute(self.rank, self.nodes, self.comm)
            # Linear response routine for f^xc
            xc_drv.integrate(fock_ao_rhs, dm_ao_rhs, gs_density,
                             molecule, basis, molgrid, self.xcfun.get_func_label())
            # Quadratic response routine for TDDFT E[3] term g^xc
            xc_drv.integrate(fock_gxc_ao, perturbed_dm_ao, zero_dm_ao,
                             gs_density, molecule, basis, molgrid,
                             self.xcfun.get_func_label(), "qrf")

            fock_gxc_ao.reduce_sum(self.rank, self.nodes, self.comm)

        fock_ao_rhs.reduce_sum(self.rank, self.nodes, self.comm)

        # Calculate the RHS and transform it to the MO basis
        if self.rank == mpi_master():
            # extract the 1PDM contributions
            fock_ao_rhs_1dm = np.zeros((dof**2, nao, nao))
            for i in range(dof**2):
                fock_ao_rhs_1dm[i] = fock_ao_rhs.alpha_to_numpy(i)

            # Transform to MO basis
            fock_mo_rhs_1dm = np.einsum('mi,xmn,na->xia', mo_occ, fock_ao_rhs_1dm, mo_vir)

            # extract the xpy and xmy contributions
            fock_ao_rhs_xpy = np.zeros((dof, nao, nao))
            fock_ao_rhs_xmy = np.zeros((dof, nao, nao))
            for i in range(dof):
                fock_ao_rhs_xpy[i] = fock_ao_rhs.alpha_to_numpy(dof**2 + i)
                fock_ao_rhs_xmy[i] = fock_ao_rhs.alpha_to_numpy(dof**2 + dof + i)

            # Contract with the second vector
            #sdp_pds = 0.5 * (
            #         np.einsum('mr,xrp,ypz->xymz', ovlp, xpy_ao, fock_ao_rhs_xpy.transpose(0,2,1))
            #        +np.einsum('mr,xrp,ypz->xymz', ovlp, xmy_ao, fock_ao_rhs_xmy.transpose(0,2,1))
            #        -np.einsum('rz,xpr,ypm->xymz', ovlp, xpy_ao, fock_ao_rhs_xpy.transpose(0,2,1))
            #        -np.einsum('rz,xpr,ypm->xymz', ovlp, xmy_ao, fock_ao_rhs_xmy.transpose(0,2,1))
            #        -np.einsum('mr,xrp,yzp->xymz', ovlp, xpy_ao, fock_ao_rhs_xpy.transpose(0,2,1))
            #        +np.einsum('mr,xrp,yzp->xymz', ovlp, xmy_ao, fock_ao_rhs_xmy.transpose(0,2,1))
            #        +np.einsum('rz,xpr,ymp->xymz', ovlp, xpy_ao, fock_ao_rhs_xpy.transpose(0,2,1))
            #        -np.einsum('rz,xpr,ymp->xymz', ovlp, xmy_ao, fock_ao_rhs_xmy.transpose(0,2,1))
            #    )

            fock_mo_rhs_2dm = 0.25 * (
                 np.einsum('cl,ti,la,xytc->xyia', ovlp, mo_occ, mo_vir, np.einsum('xmt,ymc->xytc', fock_ao_rhs_xpy.transpose(0,2,1), xpy_ao))
               - np.einsum('cl,ti,la,xytc->xyia', ovlp, mo_occ, mo_vir, np.einsum('xmt,ymc->xytc', fock_ao_rhs_xmy.transpose(0,2,1), xmy_ao))
               + np.einsum('cl,ti,la,xytc->xyia', ovlp, mo_occ, mo_vir, np.einsum('ymt,xmc->xytc', fock_ao_rhs_xpy.transpose(0,2,1), xpy_ao))
               - np.einsum('cl,ti,la,xytc->xyia', ovlp, mo_occ, mo_vir, np.einsum('ymt,xmc->xytc', fock_ao_rhs_xmy.transpose(0,2,1), xmy_ao)) 
                
               - np.einsum('cl,ti,la,xytc->xyia', ovlp, mo_occ, mo_vir, np.einsum('xtm,ymc->xytc', fock_ao_rhs_xpy.transpose(0,2,1), xpy_ao))
               - np.einsum('cl,ti,la,xytc->xyia', ovlp, mo_occ, mo_vir, np.einsum('xtm,ymc->xytc', fock_ao_rhs_xmy.transpose(0,2,1), xmy_ao))
               - np.einsum('cl,ti,la,xytc->xyia', ovlp, mo_occ, mo_vir, np.einsum('ytm,xmc->xytc', fock_ao_rhs_xpy.transpose(0,2,1), xpy_ao))
               - np.einsum('cl,ti,la,xytc->xyia', ovlp, mo_occ, mo_vir, np.einsum('ytm,xmc->xytc', fock_ao_rhs_xmy.transpose(0,2,1), xmy_ao)) 
                
               + np.einsum('cl,li,ta,xyct->xyia', ovlp, mo_occ, mo_vir, np.einsum('xmt,ycm->xyct', fock_ao_rhs_xpy.transpose(0,2,1), xpy_ao))
               + np.einsum('cl,li,ta,xyct->xyia', ovlp, mo_occ, mo_vir, np.einsum('xmt,ycm->xyct', fock_ao_rhs_xmy.transpose(0,2,1), xmy_ao))
               + np.einsum('cl,li,ta,xyct->xyia', ovlp, mo_occ, mo_vir, np.einsum('ymt,xcm->xyct', fock_ao_rhs_xpy.transpose(0,2,1), xpy_ao))
               + np.einsum('cl,li,ta,xyct->xyia', ovlp, mo_occ, mo_vir, np.einsum('ymt,xcm->xyct', fock_ao_rhs_xmy.transpose(0,2,1), xmy_ao)) 
                
               - np.einsum('cl,li,ta,xytc->xyia', ovlp, mo_occ, mo_vir, np.einsum('xtm,ycm->xytc', fock_ao_rhs_xpy.transpose(0,2,1), xpy_ao))
               + np.einsum('cl,li,ta,xytc->xyia', ovlp, mo_occ, mo_vir, np.einsum('xtm,ycm->xytc', fock_ao_rhs_xmy.transpose(0,2,1), xmy_ao))
               - np.einsum('cl,li,ta,xytc->xyia', ovlp, mo_occ, mo_vir, np.einsum('ytm,xcm->xytc', fock_ao_rhs_xpy.transpose(0,2,1), xpy_ao))
               + np.einsum('cl,li,ta,xytc->xyia', ovlp, mo_occ, mo_vir, np.einsum('ytm,xcm->xytc', fock_ao_rhs_xmy.transpose(0,2,1), xmy_ao))
               ).reshape(dof**2, nocc, nvir)

            # Symmetrize wrt. Cartesian components
            #sdp_pds_sym = 0.5 * (sdp_pds + sdp_pds.transpose(1,0,2,3))

            # Transform 2PDM contributions to MO basis
            #fock_mo_rhs_2dm = np.einsum('mi,xymn,na->xyia', mo_occ, sdp_pds_sym, mo_vir).reshape(dof**2, nocc, nvir)

            # Calculate the dipole contributions to the RHS:
            # Dipole integrals in AO basis
            dipole_drv = ElectricDipoleIntegralsDriver(self.comm)
            dipole_mats = dipole_drv.compute(molecule, basis)
            dipole_ints_ao = np.zeros((dof, nao, nao))
            k = 0
            if 'x' in self.vector_components:
                dipole_ints_ao[k] = dipole_mats.x_to_numpy()
                k += 1
            if 'y' in self.vector_components:
                dipole_ints_ao[k] = dipole_mats.y_to_numpy()
                k += 1
            if 'z' in self.vector_components:
                dipole_ints_ao[k] = dipole_mats.z_to_numpy()

            # Transform them to MO basis (oo and vv blocks only)
            dipole_ints_oo = np.array([np.linalg.multi_dot([mo_occ.T, dipole_ints_ao[x], mo_occ]) for x in range(dof)])
            dipole_ints_vv = np.array([np.linalg.multi_dot([mo_vir.T, dipole_ints_ao[x], mo_vir]) for x in range(dof)])

            # Contract with vectors to get dipole contribution to the RHS
            rhs_dipole_contrib = 0.5 * ( np.einsum('xja,yji->xyia', xmy, dipole_ints_oo)
                                        +np.einsum('yja,xji->xyia', xmy, dipole_ints_oo)
                                       ).reshape(dof**2, nocc, nvir)

            rhs_dipole_contrib += 0.5 * (-np.einsum('xib,yab->xyia', xmy, dipole_ints_vv)
                                         -np.einsum('yib,xab->xyia', xmy, dipole_ints_vv)
                                        ).reshape(dof**2, nocc, nvir)

            rhs_mo = fock_mo_rhs_1dm + fock_mo_rhs_2dm + rhs_dipole_contrib

            # Add DFT E[3] contribution to the RHS:
            if self.dft:
                gxc_ao = np.zeros((dof**2, nao, nao))

                for i in range(dof**2):
                    gxc_ao[i] = fock_gxc_ao.alpha_to_numpy(2*i)

                gxc_mo = np.einsum('mi,xmn,na->xia', mo_occ, gxc_ao, mo_vir)
                # different factor compared to TDDFT orbital response because here vectors are scaled by 1/sqrt(2)
                rhs_mo += 0.5 * gxc_mo

        self.profiler.stop_timer('RHS')

        if self.rank == mpi_master():
            return {
                'cphf_rhs': rhs_mo,
                'dm_oo': dm_oo,
                'dm_vv': dm_vv,
                'xpy_ao': xpy_ao,
                'xmy_ao': xmy_ao,
                'unrel_dm_ao': unrel_dm_ao,
                'fock_ao_rhs': fock_ao_rhs,
                'fock_gxc_ao': fock_gxc_ao, # None if not DFT
                # TODO: remove from dictionary once the gradient is correct
                'fock_mo_rhs_1dm': fock_mo_rhs_1dm,
                'fock_mo_rhs_2dm': fock_mo_rhs_2dm,
                'rhs_dipole_contrib': rhs_dipole_contrib,
            }
        else:
            return {}

    # NOTES:
    #   - epsilon_dm_ao not returned from cphfsolver, to be calculated inside compute_omega
    #   - fock_ao_rhs and fock_gxc_ao come from cphfsolver dictionary
    #   - fock_lambda not returned yet, put in dictionary from cphfsolver (otherwise needs to be recalculated)
    def compute_omega(self, molecule, basis, scf_tensors, lr_results):
        """
        Calculates the polarizability Lagrange multipliers for the overlap matrix.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        :param scf_tensors:
            The tensors from the converged SCF calculation.
        :param lr_results:
            The results from the linear response calculation.
        """

        # ERI information
        eri_dict = self.init_eri(molecule, basis)
        # DFT information
        dft_dict = self.init_dft(molecule, scf_tensors)
        # PE information
        pe_dict = self.init_pe(molecule, basis)

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

            # Transform the overlap matrix to MO basis
            # (required for the frequency-dependent contribuion)
            ovlp_oo = np.linalg.multi_dot([mo_occ.T, ovlp, mo_occ])
            ovlp_vv = np.linalg.multi_dot([mo_vir.T, ovlp, mo_vir])
 
            # Get fock matrices from cphf_results
            fock_ao_rhs = self.cphf_results['fock_ao_rhs']
            fock_gxc_ao = self.cphf_results['fock_gxc_ao']
            dm_oo = self.cphf_results['dm_oo']
            dm_vv = self.cphf_results['dm_vv']
            cphf_ov = self.cphf_results['cphf_ov']

            # TODO: do we keep this factor like that?
            sqrt2 = np.sqrt(2.0)

            # Take response vectors and convert to matrix form
            exc_vec = 1/sqrt2 * np.array([lr_results['solutions'][x,
                                    self.frequency][:nocc * nvir].reshape(nocc, nvir)
                                    for x in self.vector_components])
            deexc_vec = 1/sqrt2 * np.array([lr_results['solutions'][x,
                                    self.frequency][nocc * nvir:].reshape(nocc, nvir)
                                    for x in self.vector_components])

            # Number of vector components
            dof = exc_vec.shape[0]

            xpy = exc_vec + deexc_vec
            xmy = exc_vec - deexc_vec


            # Get dipole moment integrals and transform to MO
            dipole_drv = ElectricDipoleIntegralsDriver(self.comm)
            dipole_mats = dipole_drv.compute(molecule, basis)
            dipole_ints_ao = np.zeros((dof, nao, nao))
            k = 0
            if 'x' in self.vector_components:
                dipole_ints_ao[k] = dipole_mats.x_to_numpy()
                k += 1
            if 'y' in self.vector_components:
                dipole_ints_ao[k] = dipole_mats.y_to_numpy()
                k += 1
            if 'z' in self.vector_components:
                dipole_ints_ao[k] = dipole_mats.z_to_numpy()

            # Transform them to MO basis (oo and ov blocks only)
            dipole_ints_oo = np.array([np.linalg.multi_dot([mo_occ.T, dipole_ints_ao[x], mo_occ]) for x in range(dof)])
            dipole_ints_ov = np.array([np.linalg.multi_dot([mo_occ.T, dipole_ints_ao[x], mo_vir]) for x in range(dof)])


            # Calculate the dipole moment integrals' contribution to omega
            dipole_ints_contrib_oo = 0.5 * ( np.einsum('xjc,yic->xyij',
                                                xmy, dipole_ints_ov)
                                            +  np.einsum('yjc,xic->xyij',
                                                xmy, dipole_ints_ov)
                                            )
            dipole_ints_contrib_ov = 0.5 * ( np.einsum('xka,yki->xyia',
                                                xmy, dipole_ints_oo)
                                            +  np.einsum('yka,xki->xyia',
                                                xmy, dipole_ints_oo)
                                            )

            dipole_ints_contrib_vv = 0.5 * ( np.einsum('xkb,yka->xyab',
                                                xmy, dipole_ints_ov)
                                            +  np.einsum('ykb,xka->xyab',
                                                xmy, dipole_ints_ov)
                                            )
            dipole_ints_contrib_ao = ( np.einsum('mi,xyij,nj->xymn', mo_occ,
                                                 dipole_ints_contrib_oo, mo_occ)
                                     + np.einsum('mi,xyia,na->xymn', mo_occ,
                                                 dipole_ints_contrib_ov, mo_vir)
                                     + np.einsum('mi,xyia,na->xymn', mo_occ,
                                                 dipole_ints_contrib_ov, mo_vir).transpose(0,1,3,2)
                                     + np.einsum('ma,xyab,nb->xymn', mo_vir,
                                                 dipole_ints_contrib_vv, mo_vir)
                                     )

            # Compute the frequency-dependent contribution to the 
            # occupied-occupied and virtual-virtual blocks of omega
            #freq_contrib_oo = + 0.0 * self.frequency * ( np.einsum('xia,yja->xyij', xpy, xmy)
            #                                           + np.einsum('xia,yja->xyij', xmy, xpy)
            #                                         #- 1.0 * np.einsum('xia,yjb->xyij', xpy, xmy)
            #                                         #- 1.0 * np.einsum('xia,yjb->xyij', xmy, xpy)
            #                                         #+ np.einsum('xia,ylb,ik->xykl', xpy, xmy, ovlp_oo)
            #                                         #+ np.einsum('xia,ylb,ik->xykl', xmy, xpy, ovlp_oo)
            #                )

            #freq_contrib_vv = + 0.0 * self.frequency * ( np.einsum('xia,yib->xyab', xpy, xmy)
            #                                           + np.einsum('xia,yib->xyab', xmy, xpy)
            #                                         #- 1.0 * np.einsum('xia,yjb->xyab', xpy, xmy)
            #                                         #- 1.0 * np.einsum('xia,yjb->xyab', xmy, xpy)
            #                                         #+ np.einsum('xia,yjd,ac->xycd', xpy, xmy, ovlp_vv)
            #                                         #+ np.einsum('xia,yjd,ac->xycd', xmy, xpy, ovlp_vv)
            #                )

            #freq_contrib_ao = ( np.einsum('mi,xyij,nj->xymn', mo_occ, freq_contrib_oo, mo_occ)
            #                  + np.einsum('ma,xyab,nb->xymn', mo_vir, freq_contrib_vv, mo_vir)
            #                    )

            # Transform the vectors to the AO basis
            xpy_ao = np.einsum('mi,xia,na->xmn', mo_occ, xpy, mo_vir)
            xmy_ao = np.einsum('mi,xia,na->xmn', mo_occ, xmy, mo_vir)

            # The density matrix; only alpha block;
            # Only works for the restricted case
            # (since scf_tensors['C'] only gives alpha block...)
            D_occ = np.matmul(mo_occ, mo_occ.T)
            D_vir = np.matmul(mo_vir, mo_vir.T)

            # Construct fock_lambda (or fock_cphf)
            cphf_ao = np.einsum('mi,xia,na->xmn', mo_occ, cphf_ov, mo_vir)
            cphf_ao_list = list([cphf_ao[x] for x in range(dof**2)])
            ao_density_cphf = AODensityMatrix(cphf_ao_list, denmat.rest)
            fock_cphf = AOFockMatrix(ao_density_cphf)

        # TODO: what has to be on MPI master and what not?
        self.comp_lr_fock(fock_cphf, ao_density_cphf, molecule,
                          basis, eri_dict, dft_dict, pe_dict, self.profiler)

        # TODO: replace for-loops with np.einsum
        # For now we:
        # - loop over indices m and n
        # - select component m or n in xpy, xmy, fock_ao_rhs and fock_lambda
        # - symmetrize with respect to m and n (only 2PDM?)
        # Notes: fock_ao_rhs is a list with dof**2 matrices corresponding to
        # the contraction of the 1PDMs with ERIs; dof matrices corresponding
        # to the contraction of xpy; and other dof matrices corresponding to
        # the contraction of xmy.

        # TODO: what shape should we use: (dof**2, nao, nao) or (dof, dof, nao, nao)?
        omega = np.zeros((dof*dof, nao, nao))

        # Calculate omega (without for-loops, only parts possible for now)
        # Construct epsilon_dm_ao
        epsilon_dm_ao = -np.einsum('mi,ii,xyij,nj->xymn', mo_occ, eo_diag, dm_oo, mo_occ)
        epsilon_dm_ao -= np.einsum('ma,aa,xyab,nb->xymn', mo_vir, ev_diag, dm_vv, mo_vir)
        epsilon_cphf_ao = np.einsum('mi,ii,xyia,na->xymn', mo_occ, eo_diag, cphf_ov.reshape(dof,dof,nocc,nvir), mo_vir)
        epsilon_dm_ao -= (epsilon_cphf_ao + epsilon_cphf_ao.transpose(0,1,3,2)) # OV + VO

        for m in range(dof):
            for n in range(dof):
                # TODO: delete commented out code
                #epsilon_dm_ao = -np.linalg.multi_dot(
                #    [mo_occ,
                #     np.matmul(eo_diag, dm_oo[m,n]), mo_occ.T])
                #epsilon_dm_ao -= np.linalg.multi_dot(
                #    [mo_vir, np.matmul(ev_diag, dm_vv[m,n]), mo_vir.T])

                #epsilon_cphf_ao = np.linalg.multi_dot(
                #    [mo_occ,
                #     np.matmul(eo_diag, cphf_ov[dof*m+n]), mo_vir.T])
                #epsilon_dm_ao -= (epsilon_cphf_ao + epsilon_cphf_ao.T)

                # TODO: move outside for-loop when all Fock matrices can be
                # extracted into a numpy array at the same time.

                # Because the excitation vector is not symmetric,
                # we need both the matrix (OO block in omega, and probably VO)
                # and its transpose (VV, OV blocks)
                # this comes from the transformation of the 2PDM contribution
                # from MO to AO basis
                fock_ao_rhs_1_m = fock_ao_rhs.alpha_to_numpy(dof**2+m)  # xpy
                fock_ao_rhs_2_m = fock_ao_rhs.alpha_to_numpy(dof**2+dof+m)  # xmy

                fock_ao_rhs_1_n = fock_ao_rhs.alpha_to_numpy(dof**2+n)  # xpy
                fock_ao_rhs_2_n = fock_ao_rhs.alpha_to_numpy(dof**2+dof+n)  # xmy

                Fp1_vv = 0.25 * (  np.linalg.multi_dot([ fock_ao_rhs_1_m.T, xpy_ao[n], ovlp.T])
                                 + np.linalg.multi_dot([ fock_ao_rhs_1_n.T, xpy_ao[m], ovlp.T]) )

                Fm1_vv = 0.25 * ( np.linalg.multi_dot([fock_ao_rhs_2_m.T, xmy_ao[n], ovlp.T])
                                + np.linalg.multi_dot([fock_ao_rhs_2_n.T, xmy_ao[m], ovlp.T]) )


                Fp2_vv = 0.25 * ( np.linalg.multi_dot([fock_ao_rhs_1_m, xpy_ao[n], ovlp.T])
                                + np.linalg.multi_dot([fock_ao_rhs_1_n, xpy_ao[m], ovlp.T]) )

                Fm2_vv = 0.25 * ( np.linalg.multi_dot([fock_ao_rhs_2_m, xmy_ao[n], ovlp.T])
                                + np.linalg.multi_dot([fock_ao_rhs_2_n, xmy_ao[m], ovlp.T]) )


                Fp1_oo = 0.25 * ( np.linalg.multi_dot([fock_ao_rhs_1_m, xpy_ao[n].T, ovlp.T])
                                + np.linalg.multi_dot([fock_ao_rhs_1_n, xpy_ao[m].T, ovlp.T]) )

                Fm1_oo = 0.25 * ( np.linalg.multi_dot([fock_ao_rhs_2_m, xmy_ao[n].T, ovlp.T])
                                + np.linalg.multi_dot([fock_ao_rhs_2_n, xmy_ao[m].T, ovlp.T]) )

                Fp2_oo = 0.25 * ( np.linalg.multi_dot([fock_ao_rhs_1_m.T, xpy_ao[n].T, ovlp.T])
                                + np.linalg.multi_dot([fock_ao_rhs_1_n.T, xpy_ao[m].T, ovlp.T]) )

                Fm2_oo = 0.25 * ( np.linalg.multi_dot([fock_ao_rhs_2_m.T, xmy_ao[n].T, ovlp.T])
                                + np.linalg.multi_dot([fock_ao_rhs_2_n.T, xmy_ao[m].T, ovlp.T]) )
                # We see that:
                # Fp1_vv = Fp1_ov and Fm1_vv = Fm1_ov
                # Fp2_vv = Fp2_ov and Fm2_vv = Fm2_ov

                # Compute the contributions from the 2PDM and the relaxed 1PDM
                # to the omega Lagrange multipliers:
                fmat = ( fock_cphf.alpha_to_numpy(m*dof+n) +
                         fock_cphf.alpha_to_numpy(m*dof+n).T +
                         fock_ao_rhs.alpha_to_numpy(m*dof+n)
                       )
                # dof=3  (0,0), (0,1), (0,2); (1,0), (1,1), (1,2), (2,0), (2,1), (2,2) * dof gamma_{zx} =

                omega_1pdm_2pdm_contribs = - (
                    np.linalg.multi_dot([D_vir, Fp1_vv + Fm1_vv - Fp2_vv + Fm2_vv, D_vir])
                  + np.linalg.multi_dot([D_occ, Fp1_vv + Fm1_vv - Fp2_vv + Fm2_vv, D_vir])
                  + np.linalg.multi_dot([D_occ, Fp1_vv + Fm1_vv - Fp2_vv + Fm2_vv, D_vir]).T
                  + np.linalg.multi_dot([D_occ, Fp1_oo + Fm1_oo - Fp2_oo + Fm2_oo, D_occ])
                  + np.linalg.multi_dot([D_occ, fmat, D_occ])
                  )

                omega[m*dof+n] = epsilon_dm_ao[m,n] + omega_1pdm_2pdm_contribs + dipole_ints_contrib_ao[m,n] #+ freq_contrib_ao[m,n]

                if self.dft:
                    factor = -0.5
                    omega[m*dof+n] += factor * np.linalg.multi_dot([D_occ,
                                              fock_gxc_ao.alpha_to_numpy(2*(m*dof+n)), D_occ])

        # add omega multipliers in AO basis to cphf_results dictionary
        self.cphf_results['omega_ao'] = omega

    def print_cphf_header(self, title):
        self.ostream.print_blank()
        self.ostream.print_header('{:s} Setup'.format(title))
        self.ostream.print_header('=' * (len(title) + 8))
        self.ostream.print_blank()

        str_width = 70

        # print general info
        cur_str = 'Solver Type                     : '
        if self.use_subspace_solver:
            cur_str += 'Iterative Subspace Algorithm'
        else:
            cur_str += 'Conjugate Gradient'
        self.ostream.print_header(cur_str.ljust(str_width))

        cur_str = 'Max. Number of Iterations       : ' + str(self.max_iter)
        self.ostream.print_header(cur_str.ljust(str_width))
        cur_str = 'Convergence Threshold           : {:.1e}'.format(
            self.conv_thresh)
        self.ostream.print_header(cur_str.ljust(str_width))

        cur_str = 'Frequency                       : {:.5f}'.format(self.frequency)
        self.ostream.print_header(cur_str.ljust(str_width))

        cur_str = 'Vector components               : ' + self.vector_components
        self.ostream.print_header(cur_str.ljust(str_width))

        self.ostream.print_blank()
        self.ostream.flush()
