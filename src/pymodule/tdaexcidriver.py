import numpy as np
import time as tm

from .veloxchemlib import ElectronRepulsionIntegralsDriver
from .veloxchemlib import ExcitationVector
from .veloxchemlib import TDASigmaVectorDriver
from .veloxchemlib import mpi_master
from .veloxchemlib import szblock
from .veloxchemlib import molorb
from .qqscheme import get_qq_scheme
from .blockdavidson import BlockDavidsonSolver
from .molecularorbitals import MolecularOrbitals


class TDAExciDriver:
    """Implements TDA excited states computation driver.

    Implements TDA excited states computation schheme for Hatree-Fock/Kohn-Sham
    level of theory.

    Attributes
    ----------
    nstates
        The number of excited states determined by driver.
    triplet
        The triplet excited states flag.
    eri_thresh
        The electron repulsion integrals screening threshold.
    conv_thresh
        The excited states convergence threshold.
    max_iter
        The maximum number of excited states driver iterations.
    cur_iter
        The current number of excited states driver iterations.
    solver
        The eigenvalues solver.
    is_converged
        The flag for excited states convergence.
    rank
        The rank of MPI process.
    nodes
        The number of MPI processes.
    """

    def __init__(self, comm, ostream):
        """Initializes TDA excited states computation driver.

        Initializes TDA excited states computation drived to default setup.

        Parameters
        ----------
        comm
            The MPI communicator.
        ostream
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
        """Updates settings in TDA driver.

        Updates settings in TDA excited states computation driver.

        Parameters
        ----------
        settings
            The settings for the driver.
        """

        if 'nstates' in settings:
            self.nstates = settings['nstates']
        if 'spin' in settings:
            self.triplet = (settings['spin'][0].upper() == 'T')

        if 'eri_thresh' in settings:
            self.eri_thresh = settings['eri_thresh']
        if 'qq_type' in settings:
            self.qq_type = settings['qq_type']

        if 'conv_thresh' in settings:
            self.conv_thresh = settings['conv_thresh']
        if 'max_iter' in settings:
            self.max_iter = settings['max_iter']

    def compute(self, molecule, ao_basis, scf_tensors):
        """Performs TDA excited states calculation.

        Performs TDA excited states calculation using molecular data.

        Parameters
        ----------
        molecule
            The molecule.
        ao_basis
            The AO basis set.
        scf_tensors
            The tensors from converged SCF wavefunction.
        """

        if self.rank == mpi_master():
            mol_orbs = MolecularOrbitals([scf_tensors['C']], [scf_tensors['E']],
                                         molorb.rest)
        else:
            mol_orbs = MolecularOrbitals()

        # set start time

        start_time = tm.time()

        eri_drv = ElectronRepulsionIntegralsDriver(self.comm)

        qq_data = eri_drv.compute(get_qq_scheme(self.qq_type), self.eri_thresh,
                                  molecule, ao_basis)

        # set up trial excitation vectors on master node

        diag_mat, trial_vecs = self.gen_trial_vectors(mol_orbs, molecule)

        # initalize sigma vectors driver

        a2x_drv = TDASigmaVectorDriver(self.comm)

        # block Davidson algorithm setup

        self.solver = BlockDavidsonSolver()

        for i in range(self.max_iter):

            # perform linear transformation of trial vectors

            sig_vecs = a2x_drv.compute(trial_vecs, qq_data, mol_orbs, molecule,
                                       ao_basis)

            # solve eigenvalues problem on master node

            if self.rank == mpi_master():

                sig_mat = self.convert_to_sigma_matrix(sig_vecs)
                trial_mat = self.convert_to_trial_matrix(trial_vecs)

                self.solver.add_iteration_data(sig_mat, trial_mat, i)

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

    def gen_trial_vectors(self, mol_orbs, molecule):
        """Generates set of TDA trial vectors.

        Generates set of TDA trial vectors for given number of excited states
        by selecting primitive excitations wirh lowest approximate energies
        E_ai = e_a-e_i.

        Parameters
        ----------
        mol_orbs
            The molecular orbitals.
        molecule
            The molecule.
        Returns
        -------
            The set of trial vectors.
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
        """Updates set of TDA trial vectors.

        Updates set of TDA trial vectors by replacing Z vector coefficients
        with values generated by reduced space solver.

        Parameters
        ----------
        trial_vecs
            The set of trial vectors.
        zvecs
            The set of Z vector coefficients.
        """

        for i in range(zvecs.shape[1]):
            for j in range(zvecs.shape[0]):
                trial_vecs[i].set_zcoefficient(zvecs[j, i], j)

    def check_convergence(self, iter):
        """Checks convergence of excitation energies and set convergence flag.

        Checks convergence of excitation energies and set convergence flag on
        all processes within MPI communicator.

        Parameters
        ----------
        iter
            The current excited states solver iteration.
        """

        self.cur_iter = iter

        if self.rank == mpi_master():
            self.is_converged = self.solver.check_convergence(self.conv_thresh)
        else:
            self.is_converged = False

        self.is_converged = self.comm.bcast(self.is_converged,
                                            root=mpi_master())

    def convert_to_sigma_matrix(self, sig_vecs):
        """Converts set of sigma vectors from C++ data format to numpy data
        format.

        Converts set of sigma vectors from std::vector<CDenseMatrix> to numpy
        2D array.

        Parameters
        ----------
        sig_vecs
            The sigma vectors as std::vector<CDenseMatrix>.
        Returns
        -------
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
        """Converts set of Z vectors from C++ data format to numpy data
        format.

        Converts set of Z vectors from std::vector<CExcitationVector> to numpy
        2D array.

        Parameters
        ----------
        trial_vecs
            The Z vectors as std::vector<CExcitationVector>.
        Returns
        -------
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

    def print_iter_data(self, iter):
        """Prints excited states solver iteration data to output stream.

        Prints excited states solver iteration data to output stream.

        Parameters
        ----------
        iter
            The current excited states solver iteration.
        """

        # iteration header

        exec_str = " *** Iteration: " + (str(iter + 1)).rjust(3)
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

    def print_excited_states(self, trial_vecs, start_time):
        """Prints excited states information to output stream.

        Prints excited states information to output stream.

        Parameters
        ----------
        start_time
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
            self.print_state_information(i, reigs[i], trial_vecs[i], rnorms[i])

    def print_state_information(self, iter, eval, evec, rnorm):

        self.ostream.print_blank()

        valstr = "Excited State No.{:3d}:".format(iter + 1)
        self.ostream.print_header(valstr.ljust(92))
        valstr = 21 * "-"
        self.ostream.print_header(valstr.ljust(92))

        valstr = "Excitation Type   : "
        if self.triplet:
            valstr += "Triplet"
        else:
            valstr += "Singlet"
        self.ostream.print_header(valstr.ljust(92))

        valstr = "Excitation Energy : {:5.8f} au".format(eval)
        self.ostream.print_header(valstr.ljust(92))

        valstr = "Residual Norm     : {:3.8f} au".format(rnorm)
        self.ostream.print_header(valstr.ljust(92))

        self.ostream.print_blank()
