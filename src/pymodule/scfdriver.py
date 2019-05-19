from collections import deque
from os.path import isfile
import numpy as np
import time as tm
import math

from .veloxchemlib import OverlapIntegralsDriver
from .veloxchemlib import KineticEnergyIntegralsDriver
from .veloxchemlib import NuclearPotentialIntegralsDriver
from .veloxchemlib import ElectronRepulsionIntegralsDriver
from .veloxchemlib import ElectricDipoleIntegralsDriver
from .veloxchemlib import mpi_master
from .aofockmatrix import AOFockMatrix
from .aodensitymatrix import AODensityMatrix
from .molecularorbitals import MolecularOrbitals
from .denguess import DensityGuess
from .qqscheme import get_qq_type
from .qqscheme import get_qq_scheme


class ScfDriver:
    """Implements SCF method (base class).

    Implements SCF method with C2-DIIS and two-level C2-DIIS convergence
    accelerators.

    Attributes
    ----------
    den_guess
        The initial density guess driver.
    acc_type
        The type of SCF convergence accelerator.
    max_err_vecs
        The maximum number of error vectors.
    max_iter
        The maximum number of SCF iterations.
    first_step
        The flag for first step in two-level C2-DIIS convergence acceleration.
    qq_type
        The electron repulsion integrals screening scheme.
    qq_dyn
        The flag for enabling dynamic thresholds in electron repulsion
        integrals screening scheme.
    conv_thresh
        The SCF convergence threshold.
    eri_thresh
        The electron repulsion integrals screening threshold.
    ovl_thresh
        The atomic orbitals linear dependency threshold.
    diis_thresh
        The C2-DIIS switch on threshold.
    use_level_shift
        The flag for usage of level shifting in SCF iterations.
    iter_data
        The list of SCF iteration data (electronic energy, electronic energy
        change, gradient, density change).
    is_converged
        The flag for SCF convergence.
    skip_iter
        The flag for SCF iteration data storage.
    old_energy
        The electronic energy of previous SCF iteration.
    num_iter
        The current number of SCF iterations.
    fock_matrices
        The list of stored Fock/Kohn-Sham matrices.
    den_matrices
        The list of stored density matrices.
    density
        The current density matrix.
    mol_orbs
        The current molecular orbitals.
    nuc_energy
        The nuclear repulsion energy of molecule.
    rank
        The rank of MPI process.
    nodes
        The number of MPI processes.
    restart
        The flag for restarting from checkpoint file
    """

    def __init__(self, comm, ostream):
        """Initializes SCF driver.

        Initializes SCF driver to default setup (convergence threshold, initial
        guess, etc).
        """

        # scf accelerator
        self.acc_type = "L2_DIIS"
        self.max_err_vecs = 10
        self.max_iter = 50
        self.first_step = False

        # screening scheme
        self.qq_type = "QQ_DEN"
        self.qq_dyn = True

        # thresholds
        self.conv_thresh = 1.0e-6
        self.ovl_thresh = 1.0e-6
        self.diis_thresh = 1000.0
        self.eri_thresh = 1.0e-12
        self.eri_thresh_tight = 1.0e-15

        # level shifting
        self.use_level_shift = False

        # iterations data
        self.iter_data = []
        self.is_converged = False
        self.skip_iter = False
        self.old_energy = 0.0
        self.num_iter = 0

        # DIIS data lists
        self.fock_matrices = deque()
        self.den_matrices = deque()

        # density matrix
        self.density = AODensityMatrix()

        # molecular orbitals
        self.mol_orbs = MolecularOrbitals()

        # nuclear repulsion energy
        self.nuc_energy = 0.0

        # mpi information
        self.comm = comm
        self.rank = self.comm.Get_rank()
        self.nodes = self.comm.Get_size()

        # output stream
        self.ostream = ostream

        # restart information
        self.restart = True
        self.checkpoint_file = None

    def update_settings(self, scf_dict):

        if 'acc_type' in scf_dict:
            self.acc_type = scf_dict['acc_type'].upper()
        if 'max_iter' in scf_dict:
            self.max_iter = int(scf_dict['max_iter'])
        if 'conv_thresh' in scf_dict:
            self.conv_thresh = float(scf_dict['conv_thresh'])
        if 'qq_type' in scf_dict:
            self.qq_type = scf_dict['qq_type'].upper()
        if 'eri_thresh' in scf_dict:
            self.eri_thresh = float(scf_dict['eri_thresh'])
        if 'restart' in scf_dict:
            key = scf_dict['restart'].lower()
            self.restart = True if key == 'yes' else False
        if 'checkpoint_file' in scf_dict:
            self.checkpoint_file = scf_dict['checkpoint_file']

    def compute(self, molecule, ao_basis, min_basis):
        """Performs SCF calculation.

        Performs SCF calculation using molecular data

        Parameters
        ----------
        molecule
            The molecule.
        ao_basis
            The AO basis set.
        min_basis
            The minimal AO basis set.
        """

        # initial guess
        if self.restart:
            self.den_guess = DensityGuess("RESTART", self.checkpoint_file)
            self.restart = self.den_guess.validate_checkpoint(
                self.rank, self.comm, molecule.elem_ids_to_numpy(),
                ao_basis.get_label())

        if self.restart:
            self.acc_type = "DIIS"
        else:
            self.den_guess = DensityGuess("SAD")

        # nuclear repulsion energy
        self.nuc_energy = molecule.nuclear_repulsion_energy()

        if self.rank == mpi_master():
            self.print_header()

        # C2-DIIS method
        if self.acc_type == "DIIS":
            self.comp_diis(molecule, ao_basis, min_basis)

        # two level C2-DIIS method
        if self.acc_type == "L2_DIIS":

            # first step
            self.first_step = True

            old_thresh = self.conv_thresh
            self.conv_thresh = 1.0e-3

            old_max_iter = self.max_iter
            self.max_iter = 5

            val_basis = ao_basis.get_valence_basis()

            self.comp_diis(molecule, val_basis, min_basis)

            # second step
            self.first_step = False

            self.conv_thresh = old_thresh

            self.max_iter = old_max_iter

            self.den_guess.guess_type = "PRCMO"

            self.comp_diis(molecule, ao_basis, val_basis)

        self.fock_matrices.clear()
        self.den_matrices.clear()

        if self.rank == mpi_master():
            self.print_scf_energy()
            self.print_ground_state(molecule)
            self.mol_orbs.print_orbitals(molecule, ao_basis, False,
                                         self.ostream)

            if (self.checkpoint_file and
                    isinstance(self.checkpoint_file, str) and
                    isfile(self.checkpoint_file)):
                checkpoint_text = "Checkpoint written to file: "
                checkpoint_text += self.checkpoint_file
                self.ostream.print_info(checkpoint_text)
                self.ostream.print_blank()

    def write_checkpoint(self, nuclear_charges, basis_set):
        """Writes molecular orbitals to checkpoint file"""

        if self.rank == mpi_master() and not self.first_step:
            if self.checkpoint_file and isinstance(self.checkpoint_file, str):
                self.mol_orbs.write_hdf5(self.checkpoint_file, nuclear_charges,
                                         basis_set)

    def comp_diis(self, molecule, ao_basis, min_basis):
        """Performs SCF calculation with C2-DIIS acceleration.

        Performs SCF calculation using molecular data

        Parameters
        ----------
        molecule
            The molecule.
        ao_basis
            The AO basis set.
        min_basis
            The minimal AO basis set.
        """

        start_time = tm.time()

        self.fock_matrices.clear()
        self.den_matrices.clear()

        ovl_mat, kin_mat, npot_mat, dipole_mats = self.comp_one_ints(
            molecule, ao_basis)

        linear_dependency = False

        if self.rank == mpi_master():
            t0 = tm.time()

            oao_mat = ovl_mat.get_ortho_matrix(self.ovl_thresh)

            self.ostream.print_info("Orthogonalization matrix computed in" +
                                    " {:.2f} sec.".format(tm.time() - t0))
            self.ostream.print_blank()

            nrow = oao_mat.number_of_rows()
            ncol = oao_mat.number_of_columns()
            linear_dependency = (nrow != ncol)

            if linear_dependency:
                ndim = nrow - ncol
                self.ostream.print_info(
                    "Removed " + str(ndim) + " linearly dependent" +
                    " vector{:s}.".format('' if ndim == 1 else 's'))
                self.ostream.print_blank()
            self.ostream.flush()

        else:
            oao_mat = None

        linear_dependency = self.comm.bcast(linear_dependency,
                                            root=mpi_master())

        if (linear_dependency and self.eri_thresh > self.eri_thresh_tight):
            self.eri_thresh = self.eri_thresh_tight

            if self.rank == mpi_master():
                self.ostream.print_info("ERI screening threshold tightened to" +
                                        " {:.1e}.".format(self.eri_thresh))
                self.ostream.print_blank()

        eri_drv = ElectronRepulsionIntegralsDriver(self.comm)

        qq_data = eri_drv.compute(get_qq_scheme(self.qq_type), self.eri_thresh,
                                  molecule, ao_basis)

        den_mat = self.comp_guess_density(molecule, ao_basis, min_basis,
                                          ovl_mat)

        den_mat.broadcast(self.rank, self.comm)

        self.density = AODensityMatrix(den_mat)

        fock_mat = AOFockMatrix(den_mat)

        if self.rank == mpi_master():
            self.print_scf_title()

        for i in self.get_scf_range():

            eri_drv.compute(fock_mat, den_mat, molecule, ao_basis, qq_data)

            fock_mat.reduce_sum(self.rank, self.nodes, self.comm)

            e_ee, e_kin, e_en = self.comp_energy(fock_mat, kin_mat, npot_mat,
                                                 den_mat)

            self.comp_full_fock(fock_mat, kin_mat, npot_mat)

            e_grad = self.comp_gradient(fock_mat, ovl_mat, den_mat, oao_mat)

            self.set_skip_iter_flag(i, e_grad)

            diff_den = self.comp_density_change(den_mat, self.density)

            self.add_iter_data(e_ee, e_kin, e_en, e_grad, diff_den)

            self.check_convergence()

            self.print_iter_data(i)

            self.store_diis_data(i, fock_mat, den_mat)

            eff_fock_mat = self.get_effective_fock(fock_mat, ovl_mat, oao_mat)

            self.mol_orbs = self.gen_molecular_orbitals(eff_fock_mat, oao_mat)

            self.write_checkpoint(molecule.elem_ids_to_numpy(),
                                  ao_basis.get_label())

            self.density = AODensityMatrix(den_mat)

            den_mat = self.gen_new_density(molecule)

            den_mat.broadcast(self.rank, self.comm)

            if self.qq_dyn:
                qq_data.set_threshold(self.get_dyn_threshold(e_grad))

            if self.is_converged:
                break

        if self.rank == mpi_master():
            self.scf_tensors = {
                'C': self.mol_orbs.alpha_to_numpy(),
                'E': self.mol_orbs.ea_to_numpy(),
                'S': ovl_mat.to_numpy(),
                'Mu': (dipole_mats.x_to_numpy(), dipole_mats.y_to_numpy(),
                       dipole_mats.z_to_numpy()),
                'D': (self.density.alpha_to_numpy(0),
                      self.density.beta_to_numpy(0)),
                'F': (fock_mat.to_numpy(0), fock_mat.to_numpy(0)),
            }
        else:
            self.scf_tensors = {
                'C': None,
                'E': None,
                'S': None,
                'Mu': None,
                'D': None,
                'F': None
            }

        if self.rank == mpi_master():
            self.print_scf_finish(start_time)

    def comp_one_ints(self, molecule, basis):
        """Computes one-electron integrals required for SCF calculation.

        Computes one-electron integrals (overlap, kinetic energy and nuclear
        potential) using molecular data.

        Parameters
        ----------
        molecule
            The molecule.
        ao_basis
            The AO basis set.
        """

        t0 = tm.time()

        ovl_drv = OverlapIntegralsDriver(self.comm)
        ovl_mat = ovl_drv.compute(molecule, basis)

        t1 = tm.time()

        kin_drv = KineticEnergyIntegralsDriver(self.comm)
        kin_mat = kin_drv.compute(molecule, basis)

        t2 = tm.time()

        npot_drv = NuclearPotentialIntegralsDriver(self.comm)
        npot_mat = npot_drv.compute(molecule, basis)

        t3 = tm.time()

        dipole_drv = ElectricDipoleIntegralsDriver(self.comm)
        dipole_mats = dipole_drv.compute(molecule, basis)

        if self.rank == mpi_master():

            self.ostream.print_info("Overlap matrix computed in" +
                                    " {:.2f} sec.".format(t1 - t0))
            self.ostream.print_blank()

            self.ostream.print_info("Kinetic energy matrix computed in" +
                                    " {:.2f} sec.".format(t2 - t1))
            self.ostream.print_blank()

            self.ostream.print_info("Nuclear potential matrix computed in" +
                                    " {:.2f} sec.".format(t3 - t2))
            self.ostream.print_blank()

            self.ostream.flush()

        return (ovl_mat, kin_mat, npot_mat, dipole_mats)

    def comp_guess_density(self, molecule, ao_basis, min_basis, ovl_mat):
        """Computes initial density guess for SCF calculation.

        Computes initial density guess for SCF using superposition of atomic
        densities or molecular orbitals projection methods.

        Parameters
        ----------
        molecule
            The molecule.
        ao_basis
            The AO basis set.
        min_basis
            The minimal AO basis set.
        ovl_mat
            The overlap matrix between minimal and full AO basis.
        """

        # guess: read from checkpoint file
        if self.den_guess.guess_type == "RESTART":

            return self.den_guess.restart_density(molecule, self.comm,
                                                  self.ostream)

        # guess: superposition of atomic densities
        if self.den_guess.guess_type == "SAD":

            return self.den_guess.sad_density(molecule, ao_basis, min_basis,
                                              ovl_mat, self.comm, self.ostream)

        # guess: projection of molecular orbitals from reduced basis
        if self.den_guess.guess_type == "PRCMO":

            if self.rank == mpi_master():
                return self.den_guess.prcmo_density(molecule, ao_basis,
                                                    min_basis, self.mol_orbs)
            else:
                return AODensityMatrix()

        return AODensityMatrix()

    def set_skip_iter_flag(self, i, e_grad):
        """Sets SCF iteration skiping flag.

        Sets SCF iteration skiping flag based on iteration number and C2-DIIS
        switch on threshold.

        Parameters
        ----------
        i
            The number of current SCF iteration.
        e_grad
            The electronic gradient at current SCF iteration.
        """

        self.num_iter = i

        self.use_level_shift = False

        if i == 0:
            self.skip_iter = True
        else:
            if e_grad < self.diis_thresh:
                self.skip_iter = False
            else:
                self.skip_iter = True

    def comp_energy(self, fock_mat, kin_mat, npot_mat, den_mat):
        """Computes SCF energy components.

            Computes SCF energy components: electronic energy, kinetic energy,
            and nuclear potential energy.

        Parameters
        ----------
        fock_mat
            The Fock/Kohn-Sham matrix (only 2e-part).
        kin_mat
            The kinetic energy matrix.
        npot_mat
            The nuclear potential matrix.
        den_mat
            The density matrix.
        Returns
        -------
            The tuple (electronic energy, kinetic energy, nuclear potential
            energy).
        """

        return (0.0, 0.0, 0.0)

    def comp_full_fock(self, fock_mat, kin_mat, npot_mat):
        """Computes full Fock/Kohn-Sham matrix.

        Computes full Fock/Kohn-Sham matrix by adding to 2e-part of
        Fock/Kohn-Sham matrix the kinetic energy and nuclear potential
        matrices.

        Parameters
        ----------
        fock_mat
            The Fock/Kohn-Sham matrix (2e-part).
        kin_mat
            The kinetic energy matrix.
        npot_mat
            The nuclear potential matrix.
        """

        return

    def comp_gradient(self, fock_mat, ovl_mat, den_mat, oao_mat):
        """Computes electronic gradient.

        Computes electronic gradient using Fock/Kohn-Sham matrix.

        Parameters
        ----------
        fock_mat
            The Fock/Kohn-Sham matrix.
        ovl_mat
            The overlap matrix.
        den_mat
            The density matrix.
        oao_mat
            The orthogonalization matrix.
        Returns
        -------
            The electronic gradient.
        """

        return 0.0

    def comp_density_change(self, den_mat, old_den_mat):
        """Computes norm of density change.

        Computes norm of density change between two density matrices.

        Parameters
        ----------
        den_mat
            The current density matrix.
        old_den_mat
            The previous density matrix.
        Returns
        -------
        The norm of change between two density matrices.
        """

        return 0.0

    def store_diis_data(self, i, fock_mat, den_mat):
        """Stores Fock/Kohn-Sham and density matrices for current iteration.

        Stores Fock/Kohn-Sham and density matrices for current iteration.

        Parameters
        ----------
        i
            The number of current SCF iteration.
        fock_mat
            The Fock/Kohn-Sham matrix.
        den_mat
            The density matrix.
        """

        return

    def get_effective_fock(self, fock_mat, ovl_mat, oao_mat):
        """Computes effective Fock/Kohn-Sham matrix in OAO basis.

        Computes effective Fock/Kohn-Sham matrix in OAO basis by applying
        Lowdin or canonical orthogonalization to AO Fock/Kohn-Sham matrix.

        Parameters
        ----------
        fock_mat
            The Fock/Kohn-Sham matrix.
        ovl_mat
            The overlap matrix.
        oao_mat
            The orthogonalization matrix.
        """

        return None

    def gen_molecular_orbitals(self, fock_mat, oao_mat):
        """Generates molecular orbitals.

        Generates molecular orbital by diagonalizing Fock/Kohn-Sham matrix.

        Parameters
        ----------
        fock_mat
            The Fock/Kohn-Sham matrix.
        oao_mat
            The orthogonalization matrix.
        Returns
        -------
            The molecular orbitals.
        """
        return MolecularOrbitals()

    def gen_new_density(self, molecule):
        """Generates density matrix.

        Generates density matrix from current molecular orbitals.

        Parameters
        ----------
        molecule
            The molecule.
        Returns
        -------
            The density matrix.
        """
        return AODensityMatrix()

    def get_dyn_threshold(self, e_grad):
        """Computes screening threshold for electron repulsion integrals.

        Computes screening threshold for electron repulsion integrals based on
        value of electronic gradient.

        Parameters
        ----------
        e_grad
            The electronic gradient.
        Returns
        -------
            The screening threshold.
        """

        if e_grad < 1.0e-6:
            return self.eri_thresh

        nteri = math.pow(10, math.floor(math.log10(e_grad)))

        nteri = 1.0e-10 * nteri

        if nteri > 1.0e-10:
            return 1.0e-10

        if nteri < self.eri_thresh:
            return self.eri_thresh

        return nteri

    def add_iter_data(self, e_ee, e_kin, e_en, e_grad, diff_den):
        """Adds SCF iteration data to SCF iterations list.

        Adds SCF iteration data (electronic energy, electronic energy change,
        electronic gradient, density difference) to SCF iterations list

        Parameters
        ----------
        e_ee
            The electronic energy.
        e_kin
            The kinetic energy.
        e_en
            The nuclear potential energy.
        e_grad
            The electronic energy gradient.
        diff_den
            The density change with respect to previous SCF iteration.
        """

        e_elec = e_ee + e_kin + e_en + self.nuc_energy

        de_elec = e_elec - self.old_energy

        self.iter_data.append((e_elec, de_elec, e_grad, diff_den))

        self.old_energy = e_elec

    def check_convergence(self):
        """Sets SCF convergence flag.

        Sets SCF convergence flag by checking if convergence condition for
        electronic gradient is fullfiled.
        """

        self.is_converged = False

        if len(self.iter_data) > 1:

            e_elec, de_elec, e_grad, diff_den = self.iter_data[-1]

            if e_grad < self.conv_thresh:
                self.is_converged = True

    def get_scf_range(self):
        """Creates range of SCF iterations.

        Creates range of SCF iterations from maximum number of SCF iterations.

        Returns
        -------
            The range of SCF iterations.
        """

        return range(self.max_iter + 1)

    def print_scf_energy(self):
        """Prints SCF energy information to output stream.

        Prints SCF energy information to output stream.

        Parameters
        ----------
        molecule
            The molecule.
        """

        return

    def print_header(self):
        """Prints SCF setup header to output stream.

        Prints SCF calculation setup details to output stream,
        """

        self.ostream.print_blank()
        self.ostream.print_header("Self Consistent Field Driver Setup")
        self.ostream.print_header(36 * "=")
        self.ostream.print_blank()

        str_width = 80
        cur_str = "Wave Function Model          : " + self.get_scf_type()
        self.ostream.print_header(cur_str.ljust(str_width))
        cur_str = "Initial Guess Model          : " + self.get_guess_type()
        self.ostream.print_header(cur_str.ljust(str_width))

        cur_str = "Convergence Accelerator      : " + self.get_acc_type()
        self.ostream.print_header(cur_str.ljust(str_width))
        cur_str = "Max. Number Of Iterations    : " + str(self.max_iter)
        self.ostream.print_header(cur_str.ljust(str_width))
        cur_str = "Max. Number Of Error Vectors : " + str(self.max_err_vecs)
        self.ostream.print_header(cur_str.ljust(str_width))
        cur_str = "Convergence Threshold        : " + \
            "{:.1e}".format(self.conv_thresh)
        self.ostream.print_header(cur_str.ljust(str_width))

        cur_str = "ERI screening scheme         : " + get_qq_type(self.qq_type)
        self.ostream.print_header(cur_str.ljust(str_width))
        cur_str = "ERI screening mode           : " + self.get_qq_dyn()
        self.ostream.print_header(cur_str.ljust(str_width))
        cur_str = "ERI Screening Threshold      : " + \
            "{:.1e}".format(self.eri_thresh)
        self.ostream.print_header(cur_str.ljust(str_width))
        cur_str = "Linear Dependence Threshold  : " + \
            "{:.1e}".format(self.ovl_thresh)
        self.ostream.print_header(cur_str.ljust(str_width))
        self.ostream.print_blank()

    def print_scf_title(self):
        """Prints SCF cycles header to output stream.

        Prints SCF cycles header to output stream,
        """

        if self.first_step:
            self.ostream.print_info("Starting Reduced Basis SCF calculation...")
        else:
            self.ostream.print_blank()
            self.ostream.print_header("Iter. |   Hartree-Fock Energy, au   | "
                                      "Energy Change, au |  Gradient Norm  | "
                                      "Density Change |")
            self.ostream.print_header(92 * "-")

    def print_scf_finish(self, start_time):
        """Prints SCF calculation finish message to output stream.

        Prints SCF calculation finish message to output stream,

        Parameters
        ----------
        start_time
            The start time of SCF calculation.
        """

        if self.first_step:
            valstr = "...done. SCF energy: "
            valstr += "{:.12f}".format(self.old_energy)
            valstr += " au. Time: "
            valstr += "{:.2f}".format(tm.time() - start_time) + " sec."
            self.ostream.print_info(valstr)
            self.ostream.print_blank()

        else:
            valstr = "*** SCF "
            if self.is_converged:
                valstr += "converged in "
            else:
                valstr += "not converged in "
            valstr += str(self.num_iter)
            valstr += " iterations. Time: "
            valstr += "{:.2f}".format(tm.time() - start_time) + " sec."
            self.ostream.print_blank()
            self.ostream.print_header(valstr.ljust(92))
            self.ostream.print_blank()

        self.ostream.flush()

    def print_iter_data(self, i):
        """Prints SCF iteration data to output stream.

        Prints SCF iteration data to output stream,

        Parameters
        ----------
        i
            The current SCF iteration.
        """

        if self.rank == mpi_master():
            # no output for first step in two level DIIS
            if self.first_step:
                return

            # DIIS or second step in two level DIIS
            if i > 0:

                if len(self.iter_data) > 0:
                    te, diff_te, e_grad, diff_den = self.iter_data[-1]

                if i == 1:
                    diff_te = 0.0
                    diff_den = 0.0

                exec_str = " " + (str(i)).rjust(3) + 4 * " "
                exec_str += ("{:7.12f}".format(te)).center(27) + 3 * " "
                exec_str += ("{:5.10f}".format(diff_te)).center(17) + 3 * " "
                exec_str += ("{:5.8f}".format(e_grad)).center(15) + 3 * " "
                exec_str += ("{:5.8f}".format(diff_den)).center(15) + " "

                self.ostream.print_header(exec_str)
                self.ostream.flush()

    def get_scf_energy(self):
        """Gets SCF energy from previous SCF iteration.

        Gets SCF energy from previous SCF iteration.

        Returns
        -------
        The SCF energy.
        """

        return self.old_energy

    def get_scf_type(self):
        """Gets string with type of SCF calculation.

        Gets string with type of SCF calculation (defined in derrived classes).

        Returns
        -------
            The string with type of SCF calculation.
        """
        return "Undefined"

    def get_guess_type(self):
        """Gets string with type of initial guess.

        Gets string with type of initial guess (superposition of atomic
        densities or projection of molecular orbitals).

        Returns
        -------
            The string with type of initial guess.
        """

        if self.den_guess.guess_type == "SAD":
            return "Superposition of Atomic Densities"

        return "Undefined"

    def get_acc_type(self):
        """Gets string type of SCF convergence accelerator.

        Gets string with type of SCF convergence accelerator (DIIS or two level
        DIIS).

        Returns
        -------
        The string with type of SCF convergence accelerator.
        """

        if self.acc_type == "DIIS":
            return "Direct Inversion of Iterative Subspace"

        if self.acc_type == "L2_DIIS":
            return "Two Level Direct Inversion of Iterative Subspace"

        return "Undefined"

    def get_qq_dyn(self):
        """Gets string with application method of electron repulsion integrals
        screening.

        Gets string with application method (static or dynamic) of electron
        repulsion integrals screening.

        Returns
        -------
        The string with application method of electron repulsion integrals
        screening.
        """

        if self.qq_dyn:
            return "Dynamic"

        return "Static"

    def need_min_basis(self):
        """Determines if minimal AO basis is needed in SCF calculation.

        Determines if minimal AO basis is needed in SCF calculation. Usage of
        two level DIIS accelerator or superposition of atomic densities initial
        guess requires minimal AO basis.

        Returns
        -------
        The flag for need of minimal AO basis.
        """

        if self.acc_type == "L2_DIIS":
            return True

        if self.den_guess.guess_type == "SAD":
            return True

        return False

    def delete_mos(self, mol_orbs, mol_eigs):
        """Generates trimmed molecular orbitals.

        Generates trimmed molecular orbital by deleting MOs with coeficients
        exceeding 1.0 / sqrt(ovl_thresh).

        Parameters
        ----------
        mol_orbs
            The molecular orbitals.
        mol_eigs
            The eigenvalues of molecular orbitals.
        Returns
        -------
            The tuple (trimmed molecular orbitals, eigenvalues).
        """

        fmax = 1.0 / math.sqrt(self.ovl_thresh)

        mvec = np.amax(np.abs(mol_orbs), axis=0)

        molist = []
        for i in range(mvec.shape[0]):
            if mvec[i] < fmax:
                molist.append(i)

        return (mol_orbs[:, molist], mol_eigs[molist])

    def print_ground_state(self, molecule):
        """Prints ground state information to output stream.

        Prints ground state information to output stream.

        Parameters
        ----------
        molecule
            The molecule.
        """

        self.ostream.print_blank()

        self.ostream.print_header("Ground State Information".ljust(92))
        self.ostream.print_header("------------------------".ljust(92))

        chg = molecule.get_charge()
        valstr = "Charge of Molecule            :{:5.1f}".format(chg)
        self.ostream.print_header(valstr.ljust(92))

        mult = molecule.get_multiplicity()
        valstr = "Multiplicity (2S+1)           :{:5.1f}".format(mult)
        self.ostream.print_header(valstr.ljust(92))

        sz = 0.5 * (mult - 1.0)
        valstr = "Magnetic Quantum Number (S_z) :{:5.1f}".format(sz)
        self.ostream.print_header(valstr.ljust(92))

        self.ostream.print_blank()

    def print_energy_components(self):
        """Prints SCF energy components to output stream.

        Prints SCF energy components to output stream.
        """

        enuc = self.nuc_energy

        etot = self.iter_data[-1][0]

        e_el = etot - enuc

        valstr = "Total Energy                       :{:20.10f} au".format(etot)
        self.ostream.print_header(valstr.ljust(92))

        valstr = "Electronic Energy                  :{:20.10f} au".format(e_el)
        self.ostream.print_header(valstr.ljust(92))

        valstr = "Nuclear Repulsion Energy           :{:20.10f} au".format(enuc)
        self.ostream.print_header(valstr.ljust(92))

        self.ostream.print_header(
            "------------------------------------".ljust(92))

        grad = self.iter_data[-1][2]
        valstr = "Gradient Norm                      :{:20.10f} au".format(grad)
        self.ostream.print_header(valstr.ljust(92))
