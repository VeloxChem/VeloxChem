import numpy as np
import time as tm
import math

from .veloxchemlib import ElectricDipoleIntegralsDriver
from .veloxchemlib import LinearMomentumIntegralsDriver
from .veloxchemlib import AngularMomentumIntegralsDriver
from .veloxchemlib import ExcitationVector
from .veloxchemlib import AODensityMatrix
from .veloxchemlib import AOFockMatrix
from .veloxchemlib import mpi_master
from .veloxchemlib import rotatory_strength_in_cgs
from .veloxchemlib import szblock
from .veloxchemlib import molorb
from .veloxchemlib import denmat
from .veloxchemlib import fockmat
from .profiler import Profiler
from .linearsolver import LinearSolver
from .blockdavidson import BlockDavidsonSolver
from .molecularorbitals import MolecularOrbitals
from .errorhandler import assert_msg_critical
from .checkpoint import read_rsp_hdf5
from .checkpoint import write_rsp_hdf5


class TDAExciDriver(LinearSolver):
    """
    Implements TDA excited states computation schheme for Hartree-Fock/Kohn-Sham
    level of theory.

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.

    Instance variables
        - nstates: The number of excited states determined by driver.
        - solver: The eigenvalues solver.
    """

    def __init__(self, comm, ostream):
        """
        Initializes TDA excited states computation drived to default setup.
        """

        super().__init__(comm, ostream)

        # excited states information
        self.nstates = 3

        # solver setup
        self.solver = None

    def update_settings(self, rsp_dict, method_dict={}):
        """
        Updates response and method settings in TDA excited states computation
        driver.

        :param rsp_dict:
            The dictionary of response dict.
        :param method_dict:
            The dictionary of method settings.
        """

        super().update_settings(rsp_dict, method_dict)

        if 'nstates' in rsp_dict:
            self.nstates = int(rsp_dict['nstates'])

    def compute(self, molecule, basis, scf_tensors):
        """
        Performs TDA excited states calculation using molecular data.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        :param scf_tensors:
            The dictionary of tensors from converged SCF wavefunction.

        :return:
            A dictionary containing eigenvalues, eigenvectors, transition
            dipole moments, oscillator strengths and rotatory strengths.
        """

        profiler = Profiler({
            'timing': self.timing,
            'profiling': self.profiling,
            'memory_profiling': self.memory_profiling,
            'memory_tracing': self.memory_tracing,
        })

        if self.rank == mpi_master():
            self.print_header('TDA Driver', nstates=self.nstates)

        # set start time

        self.start_time = tm.time()

        # sanity check

        nalpha = molecule.number_of_alpha_electrons()
        nbeta = molecule.number_of_beta_electrons()
        assert_msg_critical(
            nalpha == nbeta,
            'TDAExciDriver: not implemented for unrestricted case')

        # prepare molecular orbitals

        if self.rank == mpi_master():
            mol_orbs = MolecularOrbitals([scf_tensors['C']], [scf_tensors['E']],
                                         molorb.rest)
        else:
            mol_orbs = MolecularOrbitals()

        # ERI information
        eri_dict = self.init_eri(molecule, basis)

        # DFT information
        dft_dict = self.init_dft(molecule, scf_tensors)

        # PE information
        pe_dict = self.init_pe(molecule, basis)

        timing_dict = {}

        # set up trial excitation vectors on master node

        diag_mat, trial_vecs = self.gen_trial_vectors(mol_orbs, molecule)

        # block Davidson algorithm setup

        self.solver = BlockDavidsonSolver()

        # read initial guess from restart file

        n_restart_vectors = 0
        n_restart_iterations = 0

        if self.restart:
            if self.rank == mpi_master():
                rst_trial_mat, rst_sig_mat = read_rsp_hdf5(
                    self.checkpoint_file, ['TDA_trials', 'TDA_sigmas'],
                    molecule, basis, dft_dict, pe_dict, self.ostream)
                self.restart = (rst_trial_mat is not None and
                                rst_sig_mat is not None)
                if rst_trial_mat is not None:
                    n_restart_vectors = rst_trial_mat.shape[1]
            self.restart = self.comm.bcast(self.restart, root=mpi_master())
            n_restart_vectors = self.comm.bcast(n_restart_vectors,
                                                root=mpi_master())
            n_restart_iterations = n_restart_vectors // self.nstates
            if n_restart_vectors % self.nstates != 0:
                n_restart_iterations += 1

        profiler.check_memory_usage('Initial guess')

        # start TDA iteration

        for i in range(self.max_iter):

            profiler.start_timer(i, 'FockBuild')

            # perform linear transformation of trial vectors

            if i >= n_restart_iterations:
                fock, tdens, gsdens = self.get_densities(
                    trial_vecs, scf_tensors, molecule)

                self.comp_lr_fock(fock, tdens, molecule, basis, eri_dict,
                                  dft_dict, pe_dict, timing_dict)

            profiler.stop_timer(i, 'FockBuild')
            profiler.start_timer(i, 'ReducedSpace')

            # solve eigenvalues problem on master node

            if self.rank == mpi_master():

                if i >= n_restart_iterations:
                    trial_mat = self.convert_to_trial_matrix(trial_vecs)
                    sig_mat = self.get_sigmas(fock, scf_tensors, molecule,
                                              trial_mat)
                else:
                    istart = i * self.nstates
                    iend = (i + 1) * self.nstates
                    if iend > n_restart_vectors:
                        iend = n_restart_vectors
                    sig_mat = np.copy(rst_sig_mat[:, istart:iend])
                    trial_mat = np.copy(rst_trial_mat[:, istart:iend])

                self.solver.add_iteration_data(sig_mat, trial_mat, i)

                zvecs = self.solver.compute(diag_mat)

                self.print_iter_data(i)

                trial_vecs = self.convert_to_trial_vectors(
                    mol_orbs, molecule, zvecs)

            profiler.stop_timer(i, 'ReducedSpace')
            if self.dft or self.pe:
                profiler.update_timer(i, timing_dict)

            profiler.check_memory_usage('Iteration {:d}'.format(i + 1))

            profiler.print_memory_tracing(self.ostream)

            # check convergence

            self.check_convergence(i)

            # write checkpoint file

            if (self.rank == mpi_master() and i >= n_restart_iterations):
                trials = self.solver.trial_matrices
                sigmas = self.solver.sigma_matrices
                write_rsp_hdf5(self.checkpoint_file, [trials, sigmas],
                               ['TDA_trials', 'TDA_sigmas'], molecule, basis,
                               dft_dict, pe_dict, self.ostream)

            # finish TDA after convergence

            if self.is_converged:
                break

        # converged?
        if self.rank == mpi_master():
            self.print_convergence('{:d} excited states'.format(self.nstates))

        profiler.print_timing(self.ostream)
        profiler.print_profiling_summary(self.ostream)

        profiler.check_memory_usage('End of TDA driver')
        profiler.print_memory_usage(self.ostream)

        # compute 1e dipole integrals

        dipole_ints = self.comp_dipole_ints(molecule, basis)
        linmom_ints = self.comp_linear_momentum_ints(molecule, basis)
        angmom_ints = self.comp_angular_momentum_ints(molecule, basis)

        # print converged excited states

        if self.rank == mpi_master() and self.is_converged:
            nocc = molecule.number_of_alpha_electrons()
            mo_occ = scf_tensors['C'][:, :nocc]
            mo_vir = scf_tensors['C'][:, nocc:]

            eigvals, rnorms = self.solver.get_eigenvalues()
            eigvecs = self.solver.ritz_vectors

            elec_trans_dipoles = self.comp_elec_trans_dipoles(
                dipole_ints, eigvecs, mo_occ, mo_vir)

            velo_trans_dipoles = self.comp_velo_trans_dipoles(
                linmom_ints, eigvals, eigvecs, mo_occ, mo_vir)

            magn_trans_dipoles = self.comp_magn_trans_dipoles(
                angmom_ints, eigvecs, mo_occ, mo_vir)

            oscillator_strengths = self.comp_oscillator_strengths(
                elec_trans_dipoles, eigvals)

            rotatory_strengths = self.comp_rotatory_strengths(
                velo_trans_dipoles, magn_trans_dipoles)

            return {
                'eigenvalues': eigvals,
                'eigenvectors': eigvecs,
                'electric_transition_dipoles': elec_trans_dipoles,
                'velocity_transition_dipoles': velo_trans_dipoles,
                'magnetic_transition_dipoles': magn_trans_dipoles,
                'oscillator_strengths': oscillator_strengths,
                'rotatory_strengths': rotatory_strengths,
            }
        else:
            return {}

    def gen_trial_vectors(self, mol_orbs, molecule):
        """
        Generates set of TDA trial vectors for given number of excited states
        by selecting primitive excitations wirh lowest approximate energies
        E_ai = e_a-e_i.

        :param mol_orbs:
            The molecular orbitals.
        :param molecule:
            The molecule.

        :return:
            tuple (approximate diagonal of symmetric A, set of trial vectors).
        """

        if self.rank == mpi_master():

            nocc = molecule.number_of_alpha_electrons()
            norb = mol_orbs.number_mos()

            zvec = ExcitationVector(szblock.aa, 0, nocc, nocc, norb, True)
            exci_list = zvec.small_energy_identifiers(mol_orbs, self.nstates)

            diag_mat = zvec.diagonal_to_numpy(mol_orbs)

            trial_vecs = []
            for i in exci_list:
                trial_vecs.append(
                    ExcitationVector(szblock.aa, 0, nocc, nocc, norb, True))
                trial_vecs[-1].set_zcoefficient(1.0, i)

            return diag_mat, trial_vecs

        return None, []

    def check_convergence(self, iteration):
        """
        Checks convergence of excitation energies and set convergence flag on
        all processes within MPI communicator.

        :param iteration:
            The current excited states solver iteration.
        """

        self.cur_iter = iteration

        if self.rank == mpi_master():
            self.is_converged = self.solver.check_convergence(self.conv_thresh)
        else:
            self.is_converged = False

        self.is_converged = self.comm.bcast(self.is_converged,
                                            root=mpi_master())

    def convert_to_trial_matrix(self, trial_vecs):
        """
        Converts set of Z vectors from std::vector<CExcitationVector> to numpy
        2D array.

        :param trial_vecs:
            The Z vectors as std::vector<CExcitationVector>.

        :return:
            The 2D Numpy array.
        """

        nvecs = len(trial_vecs)

        if nvecs > 0:

            trial_mat = trial_vecs[0].zvector_to_numpy()
            for i in range(1, nvecs):
                trial_mat = np.hstack(
                    (trial_mat, trial_vecs[i].zvector_to_numpy()))

            return trial_mat

        return None

    def convert_to_trial_vectors(self, mol_orbs, molecule, zvecs):
        """
        Converts set of Z vectors from numpy 2D array to
        std::vector<CExcitationVector>.

        :param mol_orbs:
            The molecular orbitals.
        :param molecule:
            The molecule.
        :param zvecs:
            The Z vectors as 2D Numpy array.

        :return:
            The Z vectors as std::vector<CExcitationVector>.
        """

        nocc = molecule.number_of_alpha_electrons()
        norb = mol_orbs.number_mos()

        trial_vecs = []
        for i in range(zvecs.shape[1]):
            trial_vecs.append(
                ExcitationVector(szblock.aa, 0, nocc, nocc, norb, True))
            for j in range(zvecs.shape[0]):
                trial_vecs[i].set_zcoefficient(zvecs[j, i], j)

        return trial_vecs

    def get_densities(self, trial_vecs, tensors, molecule):
        """
        Computes the ground-state and transition densities, and initializes the
        Fock matrix.

        :param trial_vecs:
            The Z vectors as std::vector<CExcitationVector>.
        :param tensors:
            The dictionary of tensors from converged SCF wavefunction.
        :param molecule:
            The molecule.

        :return:
            The initialized Fock matrix, the transition density matrix, and the
            ground-state density matrix.
        """

        # form transition densities

        if self.rank == mpi_master():
            nocc = molecule.number_of_alpha_electrons()
            norb = tensors['C'].shape[1]
            nvir = norb - nocc
            mo_occ = tensors['C'][:, :nocc]
            mo_vir = tensors['C'][:, nocc:]
            ao_mats = []
            for vec in trial_vecs:
                mat = vec.zvector_to_numpy().reshape(nocc, nvir)
                mat = np.matmul(mo_occ, np.matmul(mat, mo_vir.T))
                ao_mats.append(mat)
            tdens = AODensityMatrix(ao_mats, denmat.rest)
        else:
            tdens = AODensityMatrix()
        tdens.broadcast(self.rank, self.comm)

        # initialize Fock matrices

        fock = AOFockMatrix(tdens)

        if self.dft:
            if self.xcfun.is_hybrid():
                fock_flag = fockmat.rgenjkx
                fact_xc = self.xcfun.get_frac_exact_exchange()
                for i in range(fock.number_of_fock_matrices()):
                    fock.set_scale_factor(fact_xc, i)
            else:
                fock_flag = fockmat.rgenj
        else:
            fock_flag = fockmat.rgenjk

        for i in range(fock.number_of_fock_matrices()):
            fock.set_fock_type(fock_flag, i)

        # broadcast ground state density

        if self.dft:
            if self.rank == mpi_master():
                gsdens = AODensityMatrix([tensors['D'][0]], denmat.rest)
            else:
                gsdens = AODensityMatrix()
            gsdens.broadcast(self.rank, self.comm)
        else:
            gsdens = None

        return fock, tdens, gsdens

    def get_sigmas(self, fock, tensors, molecule, trial_mat):
        """
        Computes the sigma vectors.

        :param fock:
            The Fock matrix.
        :param tensors:
            The dictionary of tensors from converged SCF wavefunction.
        :param molecule:
            The molecule.
        :param trial_mat:
            The trial vectors as 2D Numpy array.

        :return:
            The sigma vectors as 2D Numpy array.
        """

        nocc = molecule.number_of_alpha_electrons()
        norb = tensors['C'].shape[1]
        nvir = norb - nocc
        mo_occ = tensors['C'][:, :nocc]
        mo_vir = tensors['C'][:, nocc:]
        orb_ene = tensors['E']

        sigma_vecs = []
        for fockind in range(fock.number_of_fock_matrices()):
            # 2e contribution
            mat = fock.to_numpy(fockind)
            mat = np.matmul(mo_occ.T, np.matmul(mat, mo_vir))
            # 1e contribution
            cjb = trial_mat[:, fockind].reshape(nocc, nvir)
            mat += np.matmul(cjb, np.diag(orb_ene[nocc:]).T)
            mat -= np.matmul(np.diag(orb_ene[:nocc]), cjb)
            sigma_vecs.append(mat.reshape(nocc * nvir, 1))

        sigma_mat = sigma_vecs[0]
        for vec in sigma_vecs[1:]:
            sigma_mat = np.hstack((sigma_mat, vec))

        return sigma_mat

    def comp_dipole_ints(self, molecule, basis):
        """
        Computes one-electron dipole integrals.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.

        :return:
            The Cartesian components of one-electron dipole integrals.
        """

        dipole_drv = ElectricDipoleIntegralsDriver(self.comm)
        dipole_mats = dipole_drv.compute(molecule, basis)

        if self.rank == mpi_master():
            return (dipole_mats.x_to_numpy(), dipole_mats.y_to_numpy(),
                    dipole_mats.z_to_numpy())
        else:
            return ()

    def comp_linear_momentum_ints(self, molecule, basis):
        """
        Computes one-electron linear momentum integrals.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.

        :return:
            The Cartesian components of one-electron linear momentum integrals.
        """

        linmom_drv = LinearMomentumIntegralsDriver(self.comm)
        linmom_mats = linmom_drv.compute(molecule, basis)

        if self.rank == mpi_master():
            return (linmom_mats.x_to_numpy(), linmom_mats.y_to_numpy(),
                    linmom_mats.z_to_numpy())
        else:
            return ()

    def comp_angular_momentum_ints(self, molecule, basis):
        """
        Computes one-electron angular momentum integrals.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.

        :return:
            The Cartesian components of one-electron angular momentum integrals.
        """

        angmom_drv = AngularMomentumIntegralsDriver(self.comm)
        angmom_mats = angmom_drv.compute(molecule, basis)

        if self.rank == mpi_master():
            return (angmom_mats.x_to_numpy(), angmom_mats.y_to_numpy(),
                    angmom_mats.z_to_numpy())
        else:
            return ()

    def comp_elec_trans_dipoles(self, dipole_ints, eigvecs, mo_occ, mo_vir):
        """
        Computes electric transition dipole moments in length form.

        :param dipole_ints:
            One-electron dipole integrals.
        :param eigvecs:
            The CI vectors.
        :param mo_occ:
            The occupied MO coefficients.
        :param mo_vir:
            The virtual MO coefficients.

        :return:
            The electric transition dipole moments in length form.
        """

        transition_dipoles = []

        for s in range(self.nstates):
            exc_vec = eigvecs[:, s].reshape(mo_occ.shape[1], mo_vir.shape[1])
            trans_dens = np.matmul(mo_occ, np.matmul(exc_vec, mo_vir.T))
            trans_dens *= math.sqrt(2.0)

            trans_dipole = np.array(
                [np.vdot(trans_dens, dipole_ints[d]) for d in range(3)])

            transition_dipoles.append(trans_dipole)

        return transition_dipoles

    def comp_velo_trans_dipoles(self, linmom_ints, eigvals, eigvecs, mo_occ,
                                mo_vir):
        """
        Computes electric transition dipole moments in velocity form.

        :param linmom_ints:
            One-electron linear momentum integrals.
        :param eigvals:
            The excitation energies.
        :param eigvecs:
            The CI vectors.
        :param mo_occ:
            The occupied MO coefficients.
        :param mo_vir:
            The virtual MO coefficients.

        :return:
            The electric transition dipole moments in velocity form.
        """

        transition_dipoles = []

        for s in range(self.nstates):
            exc_vec = eigvecs[:, s].reshape(mo_occ.shape[1], mo_vir.shape[1])
            trans_dens = np.matmul(mo_occ, np.matmul(exc_vec, mo_vir.T))
            trans_dens *= math.sqrt(2.0)

            trans_dipole = -1.0 / eigvals[s] * np.array(
                [np.vdot(trans_dens, linmom_ints[d]) for d in range(3)])

            transition_dipoles.append(trans_dipole)

        return transition_dipoles

    def comp_magn_trans_dipoles(self, angmom_ints, eigvecs, mo_occ, mo_vir):
        """
        Computes magnetic transition dipole moments.

        :param angmom_ints:
            One-electron angular momentum integrals.
        :param eigvecs:
            The CI vectors.
        :param mo_occ:
            The occupied MO coefficients.
        :param mo_vir:
            The virtual MO coefficients.

        :return:
            The magnetic transition dipole moments.
        """

        transition_dipoles = []

        for s in range(self.nstates):
            exc_vec = eigvecs[:, s].reshape(mo_occ.shape[1], mo_vir.shape[1])
            trans_dens = np.matmul(mo_occ, np.matmul(exc_vec, mo_vir.T))
            trans_dens *= math.sqrt(2.0)

            trans_dipole = 0.5 * np.array(
                [np.vdot(trans_dens, angmom_ints[d]) for d in range(3)])

            transition_dipoles.append(trans_dipole)

        return transition_dipoles

    def comp_oscillator_strengths(self, transition_dipoles, eigvals):
        """
        Computes oscillator strengths.

        :param transition_dipoles:
            The electric transition dipole moments in length form.
        :param eigvals:
            The excitation energies.

        :return:
            The oscillator strengths.
        """

        oscillator_strengths = np.zeros((self.nstates,))

        for s in range(self.nstates):
            exc_ene = eigvals[s]
            dipole_strength = np.sum(transition_dipoles[s]**2)
            oscillator_strengths[s] = 2.0 / 3.0 * dipole_strength * exc_ene

        return oscillator_strengths

    def comp_rotatory_strengths(self, velo_trans_dipoles, magn_trans_dipoles):
        """
        Computes rotatory strengths in CGS unit.

        :param velo_trans_dipoles:
            The electric transition dipole moments in velocity form.
        :param magn_trans_dipoles:
            The magnetic transition dipole moments.

        :return:
            The rotatory strengths in CGS unit.
        """

        rotatory_strengths = np.zeros((self.nstates,))

        for s in range(self.nstates):
            rotatory_strengths[s] = np.dot(velo_trans_dipoles[s],
                                           magn_trans_dipoles[s])
            rotatory_strengths[s] *= rotatory_strength_in_cgs()

        return rotatory_strengths

    def print_iter_data(self, iteration):
        """
        Prints excited states solver iteration data to output stream.

        :param iteration:
            The current excited states solver iteration.
        """

        # iteration header

        exec_str = ' *** Iteration: ' + (str(iteration + 1)).rjust(3)
        exec_str += ' * Reduced Space: '
        exec_str += (str(self.solver.reduced_space_size())).rjust(4)
        rmax, rmin = self.solver.max_min_residual_norms()
        exec_str += ' * Residues (Max,Min): {:.2e} and {:.2e}'.format(
            rmax, rmin)
        self.ostream.print_header(exec_str)
        self.ostream.print_blank()

        # excited states information

        reigs, rnorms = self.solver.get_eigenvalues()
        for i in range(reigs.shape[0]):
            exec_str = 'State {:2d}: {:5.8f} '.format(i + 1, reigs[i])
            exec_str += 'au Residual Norm: {:3.8f}'.format(rnorms[i])
            self.ostream.print_header(exec_str.ljust(84))

        # flush output stream
        self.ostream.print_blank()
        self.ostream.flush()
