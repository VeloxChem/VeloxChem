from .veloxchemlib import mpi_master
from .lrsolver import LinearResponseSolver
from .lreigensolver import LinearResponseEigenSolver
from .tdaexcidriver import TDAExciDriver
from .errorhandler import assert_msg_critical
from .qqscheme import get_qq_type


class ResponseDriver:
    """
    Implements response driver for molecular property calculations using
    conventional Hartree-Fock/Kohn-Sham response theory.

    :param rank:
        The rank of MPI process.
    :param nodes:
        The number of MPI processes.
    :param rsp_input:
        The settings dictionary.
    """

    def __init__(self, comm, ostream):
        """
        Initializes response driver to default setup.

        :param comm:
            The MPI communicator.
        :param ostream:
            The output stream.
        """

        # default calculation type
        self.prop_type = 'ABSORPTION'
        self.tamm_dancoff = False
        self.triplet = False
        self.rsp_input = {}
        self.method_input = {}

        # mpi information
        self.comm = comm
        self.rank = self.comm.Get_rank()
        self.nodes = self.comm.Get_size()

        # output stream
        self.ostream = ostream

    def update_settings(self, rsp_input, method_input):
        """
        Updates settings in response solver.

        :param rsp_input:
            The settings dictionary.
        """

        # properties
        if rsp_input['property'].lower() == 'absorption':
            self.prop_type = 'ABSORPTION'
            if 'tamm_dancoff' in rsp_input:
                key = rsp_input['tamm_dancoff'].lower()
                self.tamm_dancoff = True if key in ['yes', 'y'] else False
            if 'spin' in rsp_input:
                key = rsp_input['spin'].lower()
                self.triplet = True if key[0] == 't' else False

        elif rsp_input['property'].lower() == 'polarizability':
            self.prop_type = 'POLARIZABILITY'

        self.rsp_input = dict(rsp_input)
        self.method_input = dict(method_input)

    def compute(self, molecule, ao_basis, scf_tensors):
        """
        Performs molecular property calculation using molecular data

        :param molecule:
            The molecule.
        :param ao_basis:
            The AO basis set.
        :param scf_tensors:
            The dictionary of tensors from converged SCF wavefunction.

        :return:
            The results from the actual response solver.
        """

        #if self.rank == mpi_master():
        #    self.print_header()

        # Linear response eigensolver

        if self.prop_type.upper() in ['ABSORPTION']:
            if not self.tamm_dancoff:
                eigensolver = LinearResponseEigenSolver(self.comm, self.ostream)
                assert_msg_critical(
                    not self.triplet,
                    'LR EigenSolver: not yet implemented for triplets')
            else:
                eigensolver = TDAExciDriver(self.comm, self.ostream)

            eigensolver.update_settings(self.rsp_input, self.method_input)

            return eigensolver.compute(molecule, ao_basis, scf_tensors)

        # Linear response solver

        if self.prop_type.upper() in ['POLARIZABILITY']:
            lr_solver = LinearResponseSolver(self.comm, self.ostream)

            lr_solver.update_settings(self.rsp_input)

            return lr_solver.compute(molecule, ao_basis, scf_tensors)

    def print_header(self):
        """
        Prints molecular property calculation setup details to output stream.
        """

        self.ostream.print_blank()
        self.ostream.print_header('Response Driver Setup')
        self.ostream.print_header(23 * '=')
        self.ostream.print_blank()

        str_width = 60

        cur_str = 'Molecular Property Type   : ' + self.prop_str()
        self.ostream.print_header(cur_str.ljust(str_width))

        if self.prop_type in ['ABSORPTION']:
            if self.tamm_dancoff:
                cur_str = "Response Equations Type   : Tamm-Dancoff"
                self.ostream.print_header(cur_str.ljust(str_width))
            cur_str = 'Number of Excited States  : ' + str(self.nstates)
            self.ostream.print_header(cur_str.ljust(str_width))

        cur_str = 'Max. Number of Iterations : ' + str(self.max_iter)
        self.ostream.print_header(cur_str.ljust(str_width))
        cur_str = 'Convergence Threshold     : ' + \
            '{:.1e}'.format(self.conv_thresh)
        self.ostream.print_header(cur_str.ljust(str_width))

        cur_str = 'ERI Screening Scheme      : ' + get_qq_type(self.qq_type)
        self.ostream.print_header(cur_str.ljust(str_width))
        cur_str = 'ERI Screening Threshold   : ' + \
            '{:.1e}'.format(self.eri_thresh)
        self.ostream.print_header(cur_str.ljust(str_width))
        self.ostream.print_blank()

        self.ostream.flush()

    def prop_str(self):
        """
        Gets string with type of molecular property calculation (Excited
        states, linear and non-linear spectroscopies).

        :return:
            The string with type of molecular property calculation.
        """

        if self.prop_type == 'POLARIZABILITY':
            return 'Polarizability'

        if self.prop_type == 'ABSORPTION':
            if not self.triplet:
                return 'Singlet Excited States'
            else:
                return "Triplet Excited States"

        return 'Undefined'
