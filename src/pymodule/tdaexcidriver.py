import numpy as np
import time as tm
import math

from .veloxchemlib import ElectronRepulsionIntegralsDriver
from .veloxchemlib import ElectricDipoleIntegralsDriver
from .veloxchemlib import LinearMomentumIntegralsDriver
from .veloxchemlib import AngularMomentumIntegralsDriver
from .veloxchemlib import ExcitationVector
from .veloxchemlib import TDASigmaVectorDriver
from .veloxchemlib import GridDriver
from .veloxchemlib import MolecularGrid
from .veloxchemlib import XCFunctional
from .veloxchemlib import mpi_master
from .veloxchemlib import szblock
from .veloxchemlib import molorb
from .veloxchemlib import rotatory_strength_in_cgs
from .veloxchemlib import parse_xc_func
from .blockdavidson import BlockDavidsonSolver
from .molecularorbitals import MolecularOrbitals
from .qqscheme import get_qq_scheme
from .qqscheme import get_qq_type
from .lrmatvecdriver import read_rsp_hdf5
from .lrmatvecdriver import write_rsp_hdf5
from .errorhandler import assert_msg_critical


class TDAExciDriver:
    """
    Implements TDA excited states computation schheme for Hartree-Fock/Kohn-Sham
    level of theory.

    :param nstates:
        The number of excited states determined by driver.
    :param triplet:
        The triplet excited states flag.
    :param eri_thresh:
        The electron repulsion integrals screening threshold.
    :param qq_type:
        The electron repulsion integrals screening scheme.
    :param dft:
        The flag for running DFT.
    :param grid_level:
        The accuracy level of DFT grid.
    :param xcfun:
        The XC functional.
    :param molgrid:
        The molecular grid.
    :param conv_thresh:
        The excited states convergence threshold.
    :param max_iter:
        The maximum number of excited states driver iterations.
    :param cur_iter:
        The current number of excited states driver iterations.
    :param solver:
        The eigenvalues solver.
    :param is_converged:
        The flag for excited states convergence.
    :param rank:
        The rank of MPI process.
    :param nodes:
        The number of MPI processes.
    :param ostream:
        The output stream.
    :param restart:
        The flag for restarting from checkpoint file.
    :param checkpoint_file:
        The name of checkpoint file.
    """

    def __init__(self, comm, ostream):
        """
        Initializes TDA excited states computation drived to default setup.

        :param comm:
            The MPI communicator.
        :param ostream:
            The output stream.
        """

        # excited states information
        self.nstates = 3
        self.triplet = False

        # ERI settings
        self.eri_thresh = 1.0e-15
        self.qq_type = 'QQ_DEN'

        # dft
        self.dft = False
        self.grid_level = 4
        self.xcfun = XCFunctional()
        self.molgrid = MolecularGrid()

        # solver setup
        self.conv_thresh = 1.0e-4
        self.max_iter = 50
        self.cur_iter = 0
        self.solver = None
        self.is_converged = None

        # mpi information
        self.comm = comm
        self.rank = self.comm.Get_rank()
        self.nodes = self.comm.Get_size()

        # output stream
        self.ostream = ostream

        # restart information
        self.restart = True
        self.checkpoint_file = None

    def update_settings(self, rsp_dict, method_dict={}):
        """
        Updates response and method settings in TDA excited states computation
        driver.

        :param rsp_dict:
            The dictionary of response dict.
        :param method_dict:
            The dictionary of method settings.
        """

        if 'nstates' in rsp_dict:
            self.nstates = int(rsp_dict['nstates'])
        if 'spin' in rsp_dict:
            self.triplet = (rsp_dict['spin'][0].upper() == 'T')

        if 'eri_thresh' in rsp_dict:
            self.eri_thresh = float(rsp_dict['eri_thresh'])
        if 'qq_type' in rsp_dict:
            self.qq_type = rsp_dict['qq_type'].upper()

        if 'conv_thresh' in rsp_dict:
            self.conv_thresh = float(rsp_dict['conv_thresh'])
        if 'max_iter' in rsp_dict:
            self.max_iter = int(rsp_dict['max_iter'])

        if 'restart' in rsp_dict:
            key = rsp_dict['restart'].lower()
            self.restart = True if key == 'yes' else False
        if 'checkpoint_file' in rsp_dict:
            self.checkpoint_file = rsp_dict['checkpoint_file']

        if 'dft' in method_dict:
            key = method_dict['dft'].lower()
            self.dft = True if key == 'yes' else False
        if 'grid_level' in method_dict:
            self.grid_level = int(method_dict['grid_level'])
        if 'xcfun' in method_dict:
            self.xcfun = parse_xc_func(method_dict['xcfun'].upper())

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

        if self.rank == mpi_master():
            mol_orbs = MolecularOrbitals([scf_tensors['C']], [scf_tensors['E']],
                                         molorb.rest)
            nocc = molecule.number_of_alpha_electrons()
            mo_occ = scf_tensors['C'][:, :nocc]
            mo_vir = scf_tensors['C'][:, nocc:]
        else:
            mol_orbs = MolecularOrbitals()

        if self.rank == mpi_master():
            self.print_header()

        # set start time

        start_time = tm.time()

        # generate integration grid
        if self.dft:
            grid_drv = GridDriver(self.comm)
            grid_drv.set_level(self.grid_level)

            grid_t0 = tm.time()
            self.molgrid = grid_drv.generate(molecule)
            n_grid_points = self.molgrid.number_of_points()
            self.molgrid.distribute(self.rank, self.nodes, self.comm)
            self.ostream.print_info(
                'Molecular grid with {0:d} points generated in {1:.2f} sec.'.
                format(n_grid_points,
                       tm.time() - grid_t0))
            self.ostream.print_blank()

        if self.dft:
            dft_func_label = self.xcfun.get_func_label().upper()
        else:
            dft_func_label = 'HF'

        eri_drv = ElectronRepulsionIntegralsDriver(self.comm)

        qq_data = eri_drv.compute(get_qq_scheme(self.qq_type), self.eri_thresh,
                                  molecule, basis)

        # set up trial excitation vectors on master node

        diag_mat, trial_vecs = self.gen_trial_vectors(mol_orbs, molecule)

        # initalize sigma vectors driver

        a2x_drv = TDASigmaVectorDriver(self.comm)

        # block Davidson algorithm setup

        self.solver = BlockDavidsonSolver()

        # read initial guess from restart file

        n_restart_vectors = 0
        n_restart_iterations = 0

        if self.restart:
            if self.rank == mpi_master():
                rst_trial_mat, rst_sig_mat = read_rsp_hdf5(
                    self.checkpoint_file, ['TDA_trials', 'TDA_sigmas'],
                    molecule.nuclear_repulsion_energy(),
                    molecule.elem_ids_to_numpy(), basis.get_label(),
                    dft_func_label, self.ostream)
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

        # start TDA iteration

        for i in range(self.max_iter):

            # perform linear transformation of trial vectors

            if i >= n_restart_iterations:
                sig_vecs = a2x_drv.compute(trial_vecs, self.triplet, qq_data,
                                           self.molgrid, self.xcfun, mol_orbs,
                                           molecule, basis)

            # solve eigenvalues problem on master node

            if self.rank == mpi_master():

                if i >= n_restart_iterations:
                    sig_mat = self.convert_to_sigma_matrix(sig_vecs)
                    trial_mat = self.convert_to_trial_matrix(trial_vecs)
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

                if i >= n_restart_iterations:
                    write_rsp_hdf5(self.checkpoint_file,
                                   self.solver.trial_matrices,
                                   self.solver.sigma_matrices,
                                   ['TDA_trials', 'TDA_sigmas'],
                                   molecule.nuclear_repulsion_energy(),
                                   molecule.elem_ids_to_numpy(),
                                   basis.get_label(), dft_func_label,
                                   self.ostream)

            # check convergence

            self.check_convergence(i)

            if self.is_converged:
                break

        if self.rank == mpi_master():
            assert_msg_critical(self.is_converged,
                                'TDA driver: failed to converge')

        # compute 1e dipole integrals

        dipole_ints = self.comp_dipole_ints(molecule, basis)
        linmom_ints = self.comp_linear_momentum_ints(molecule, basis)
        angmom_ints = self.comp_angular_momentum_ints(molecule, basis)

        # print converged excited states

        if self.rank == mpi_master():
            self.print_excited_states(start_time)

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

            nocc = molecule.number_of_electrons() // 2
            norb = mol_orbs.number_mos()

            zvec = ExcitationVector(szblock.aa, 0, nocc, nocc, norb, True)
            exci_list = zvec.small_energy_identifiers(mol_orbs, self.nstates)

            diag_mat = zvec.diagonal_to_numpy(mol_orbs)

            trial_vecs = []
            for i in exci_list:
                trial_vecs.append(
                    ExcitationVector(szblock.aa, 0, nocc, nocc, norb, True))
                trial_vecs[-1].set_zcoefficient(1.0, i)

            return (diag_mat, trial_vecs)

        return (None, [])

    def update_trial_vectors(self, trial_vecs, zvecs):
        """
        Updates set of TDA trial vectors by replacing Z vector coefficients
        with values generated by reduced space solver.

        :param trial_vecs:
            The set of trial vectors.
        :param zvecs:
            The set of Z vector coefficients.
        """

        for i in range(zvecs.shape[1]):
            for j in range(zvecs.shape[0]):
                trial_vecs[i].set_zcoefficient(zvecs[j, i], j)

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

    def convert_to_sigma_matrix(self, sig_vecs):
        """
        Converts set of sigma vectors from std::vector<CDenseMatrix> to numpy
        2D array.

        :param sig_vecs:
            The sigma vectors as std::vector<CDenseMatrix>.

        :return:
            The 2D numpy array.
        """

        nvecs = len(sig_vecs)

        if nvecs > 0:

            sig_mat = sig_vecs[0].to_numpy()
            for i in range(1, nvecs):
                sig_mat = np.hstack((sig_mat, sig_vecs[i].to_numpy()))

            return sig_mat

        return None

    def convert_to_trial_matrix(self, trial_vecs):
        """
        Converts set of Z vectors from std::vector<CExcitationVector> to numpy
        2D array.

        :param trial_vecs:
            The Z vectors as std::vector<CExcitationVector>.

        :return:
            The 2D numpy array.
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
        :param trial_vecs:
            The Z vectors as 2D numpy array.

        :return:
            The Z vectors as std::vector<CExcitationVector>.
        """

        nocc = molecule.number_of_electrons() // 2
        norb = mol_orbs.number_mos()

        trial_vecs = []
        for i in range(zvecs.shape[1]):
            trial_vecs.append(
                ExcitationVector(szblock.aa, 0, nocc, nocc, norb, True))
            for j in range(zvecs.shape[0]):
                trial_vecs[i].set_zcoefficient(zvecs[j, i], j)

        return trial_vecs

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

    def print_header(self):
        """
        Prints TDA driver setup header to output stream.
        """

        self.ostream.print_blank()
        self.ostream.print_header("TDA Driver Setup")
        self.ostream.print_header(18 * "=")
        self.ostream.print_blank()

        str_width = 60

        cur_str = "Number of States          : " + str(self.nstates)
        self.ostream.print_header(cur_str.ljust(str_width))

        cur_str = "Max. Number of Iterations : " + str(self.max_iter)
        self.ostream.print_header(cur_str.ljust(str_width))
        cur_str = "Convergence Threshold     : " + \
            "{:.1e}".format(self.conv_thresh)
        self.ostream.print_header(cur_str.ljust(str_width))

        cur_str = "ERI Screening Scheme      : " + get_qq_type(self.qq_type)
        self.ostream.print_header(cur_str.ljust(str_width))
        cur_str = "ERI Screening Threshold   : " + \
            "{:.1e}".format(self.eri_thresh)
        self.ostream.print_header(cur_str.ljust(str_width))
        self.ostream.print_blank()

        self.ostream.flush()

    def print_iter_data(self, iteration):
        """
        Prints excited states solver iteration data to output stream.

        :param iteration:
            The current excited states solver iteration.
        """

        # iteration header

        exec_str = " *** Iteration: " + (str(iteration + 1)).rjust(3)
        exec_str += " * Reduced Space: "
        exec_str += (str(self.solver.reduced_space_size())).rjust(4)
        rmax, rmin = self.solver.max_min_residual_norms()
        exec_str += " * Residues (Max,Min): {:.2e} and {:.2e}".format(
            rmax, rmin)
        self.ostream.print_header(exec_str)
        self.ostream.print_blank()

        # excited states information

        reigs, rnorms = self.solver.get_eigenvalues()
        for i in range(reigs.shape[0]):
            exec_str = "State {:2d}: {:5.8f} ".format(i + 1, reigs[i])
            exec_str += "au Residual Norm: {:3.8f}".format(rnorms[i])
            self.ostream.print_header(exec_str.ljust(84))

        # flush output stream
        self.ostream.print_blank()
        self.ostream.flush()

    def print_excited_states(self, start_time):
        """
        Prints excited states information to output stream.

        :param start_time:
            The start time of SCF calculation.
        """

        self.ostream.print_blank()

        valstr = "*** {:d} excited states ".format(self.nstates)
        if self.is_converged:
            valstr += "converged"
        else:
            valstr += "NOT converged"
        valstr += " in {:d} iterations. ".format(self.cur_iter + 1)
        valstr += "Time: {:.2f}".format(tm.time() - start_time) + " sec."
        self.ostream.print_header(valstr.ljust(92))

        reigs, rnorms = self.solver.get_eigenvalues()

        for i in range(reigs.shape[0]):
            self.print_state_information(i, reigs[i], rnorms[i])

    def print_state_information(self, iteration, eigval, rnorm):
        """
        Prints excited state information to output stream.

        :param iteration:
            The current excited states solver iteration.
        :param eigval:
            The excitation energy.
        :param rnorm:
            The residual norm.
        """

        self.ostream.print_blank()

        valstr = "Excited State No.{:3d}:".format(iteration + 1)
        self.ostream.print_header(valstr.ljust(92))
        valstr = 21 * "-"
        self.ostream.print_header(valstr.ljust(92))

        valstr = "Excitation Type   : "
        if self.triplet:
            valstr += "Triplet"
        else:
            valstr += "Singlet"
        self.ostream.print_header(valstr.ljust(92))

        valstr = "Excitation Energy : {:5.8f} au".format(eigval)
        self.ostream.print_header(valstr.ljust(92))

        valstr = "Residual Norm     : {:3.8f} au".format(rnorm)
        self.ostream.print_header(valstr.ljust(92))

        self.ostream.print_blank()
