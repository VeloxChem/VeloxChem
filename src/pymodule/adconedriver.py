import numpy as np
import time as tm

from .veloxchemlib import ElectronRepulsionIntegralsDriver
from .veloxchemlib import ExcitationVector
from .veloxchemlib import mpi_master
from .veloxchemlib import szblock
from .veloxchemlib import molorb
from .veloxchemlib import denmat
from .veloxchemlib import fockmat
from .aodensitymatrix import AODensityMatrix
from .aofockmatrix import AOFockMatrix
from .qqscheme import get_qq_scheme
from .qqscheme import get_qq_type
from .blockdavidson import BlockDavidsonSolver
from .molecularorbitals import MolecularOrbitals


class AdcOneDriver:
    """
    Implements ADC(1) computation schheme for Hartree-Fock reference.

    :param nstates:
        The number of excited states determined by driver.
    :param triplet:
        The triplet excited states flag.
    :param eri_thresh:
        The electron repulsion integrals screening threshold.
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
    """

    def __init__(self, comm, ostream):
        """
        Initializes ADC(1) computation drived to default setup.

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

    def update_settings(self, settings):
        """
        Updates settings in ADC(1) computation driver.

        :param settings:
            The settings dictionary for the driver.
        """

        # calculation type
        if 'nstates' in settings:
            self.nstates = int(settings['nstates'])
        if 'spin' in settings:
            self.triplet = (settings['spin'][0].upper() == 'T')

        # solver settings
        if 'conv_thresh' in settings:
            self.conv_thresh = float(settings['conv_thresh'])
        if 'max_iter' in settings:
            self.max_iter = int(settings['max_iter'])

        # ERI settings
        if 'eri_thresh' in settings:
            self.eri_thresh = float(settings['eri_thresh'])
        if 'qq_type' in settings:
            self.qq_type = settings['qq_type'].upper()

    def compute(self, molecule, basis, scf_tensors):
        """
        Performs ADC(1) calculation using molecular data.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        :param scf_tensors:
            The dictionary of tensors from converged SCF wavefunction.

        :return:
            The excitation energies (eigenvalues).
        """

        if self.rank == mpi_master():
            mol_orbs = MolecularOrbitals([scf_tensors['C']], [scf_tensors['E']],
                                         molorb.rest)

            nocc = molecule.number_of_alpha_electrons()
            nvir = scf_tensors['C'].shape[1] - nocc

            Cocc = scf_tensors['C'][:, :nocc]
            Cvir = scf_tensors['C'][:, nocc:]

            fmo = np.matmul(scf_tensors['C'].T,
                            np.matmul(scf_tensors['F'][0], scf_tensors['C']))
            fab = fmo[nocc:, nocc:]
            fij = fmo[:nocc, :nocc]

            self.print_header()

        else:
            mol_orbs = MolecularOrbitals()

        # set start time

        start_time = tm.time()

        eri_drv = ElectronRepulsionIntegralsDriver(self.comm)

        screening = eri_drv.compute(get_qq_scheme(self.qq_type),
                                    self.eri_thresh, molecule, basis)

        # set up trial excitation vectors on master node

        diag_mat, trial_vecs = self.gen_trial_vectors(mol_orbs, molecule)

        # block Davidson algorithm setup

        self.solver = BlockDavidsonSolver()

        for i in range(self.max_iter):

            if self.rank == mpi_master():
                ao_mats = []
                for vec in trial_vecs:
                    mat = vec.zvector_to_numpy().reshape(nocc, nvir)
                    mat = np.matmul(Cocc, np.matmul(mat, Cvir.T))
                    ao_mats.append(mat)
                density = AODensityMatrix(ao_mats, denmat.rest)
            else:
                density = AODensityMatrix()
            density.broadcast(self.rank, self.comm)

            fock = AOFockMatrix(density)
            focktype = fockmat.rgenk if self.triplet else fockmat.rgenjk
            for fockind in range(fock.number_of_fock_matrices()):
                fock.set_fock_type(focktype, fockind)

            eri_drv = ElectronRepulsionIntegralsDriver(self.comm)
            eri_drv.compute(fock, density, molecule, basis, screening)
            fock.reduce_sum(self.rank, self.nodes, self.comm)

            # solve eigenvalues problem on master node

            if self.rank == mpi_master():

                trial_mat = trial_vecs[0].zvector_to_numpy()
                for vec in trial_vecs[1:]:
                    trial_mat = np.hstack((trial_mat, vec.zvector_to_numpy()))

                sigma_vecs = []
                for fockind in range(fock.number_of_fock_matrices()):
                    # 2e contribution
                    prefactor = -1.0 if self.triplet else 1.0
                    mat = prefactor * fock.to_numpy(fockind)
                    mat = np.matmul(Cocc.T, np.matmul(mat, Cvir))
                    # 1e contribution
                    cjb = trial_mat[:, fockind].reshape(nocc, nvir)
                    mat += np.matmul(cjb, fab.T)
                    mat -= np.matmul(fij, cjb)
                    sigma_vecs.append(mat.reshape(nocc * nvir, 1))

                sigma_mat = sigma_vecs[0]
                for vec in sigma_vecs[1:]:
                    sigma_mat = np.hstack((sigma_mat, vec))

                self.solver.add_iteration_data(sigma_mat, trial_mat, i)

                zvecs = self.solver.compute(diag_mat)

                self.print_iter_data(i)

                self.update_trial_vectors(trial_vecs, zvecs)

            # check convergence

            self.check_convergence(i)

            if self.is_converged:
                break

        # print converged excited states

        if self.rank == mpi_master():
            self.print_excited_states(trial_vecs, start_time)

            reigs, rnorms = self.solver.get_eigenvalues()
            return reigs
        else:
            return None

    def print_header(self):
        """
        Prints ADC(1) driver setup header to output stream.
        """

        self.ostream.print_blank()
        self.ostream.print_header("ADC(1) Driver Setup")
        self.ostream.print_header(21 * "=")
        self.ostream.print_blank()

        str_width = 60

        cur_str = "Number Of Excited States  : " + str(self.nstates)
        self.ostream.print_header(cur_str.ljust(str_width))

        cur_str = "Max. Number Of Iterations : " + str(self.max_iter)
        self.ostream.print_header(cur_str.ljust(str_width))
        cur_str = "Convergence Threshold     : " + \
            "{:.1e}".format(self.conv_thresh)
        self.ostream.print_header(cur_str.ljust(str_width))

        cur_str = "ERI screening scheme      : " + get_qq_type(self.qq_type)
        self.ostream.print_header(cur_str.ljust(str_width))
        cur_str = "ERI Screening Threshold   : " + \
            "{:.1e}".format(self.eri_thresh)
        self.ostream.print_header(cur_str.ljust(str_width))
        self.ostream.print_blank()

        self.ostream.flush()

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
            tuple (approximate diagonal of A, set of trial vectors).
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
            valstr += "not converged"
        valstr += " in {:d} iterations. ".format(self.cur_iter + 1)
        valstr += "Time: {:.2f}".format(tm.time() - start_time) + " sec."
        self.ostream.print_header(valstr.ljust(92))

        reigs, rnorms = self.solver.get_eigenvalues()

        for i in range(reigs.shape[0]):
            self.print_state_information(i, reigs[i], rnorms[i])

    def print_state_information(self, iteration, eigval, rnorm):
        """
        Prints excited state information to output stream.
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
