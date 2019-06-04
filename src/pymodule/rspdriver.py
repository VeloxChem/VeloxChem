from .veloxchemlib import mpi_master
from .lrsolver import LinearResponseSolver
from .lreigensolver import LinearResponseEigenSolver
from .tdaexcidriver import TDAExciDriver
from .errorhandler import assert_msg_critical
from .qqscheme import get_qq_type


class ResponseDriver:
    """Implements response driver.

    Implements response driver for molecular property calculations using
    conventional Hatree-Fock/Kohn-Sham response theory.

    rank
        The rank of MPI process.
    nodes
        The number of MPI processes.
    """

    def __init__(self, comm, ostream):
        """Initializes Response driver.

        Initializes Response driver to default setup.

        Parameters
        ----------
        comm
            The MPI communicator.
        ostream
            The output stream.
        """

        # calculation type
        self.prop_type = 'EXCITATION'
        self.nstates = 3
        self.tamm_dancoff = False
        self.triplet = False

        # solver settings
        self.conv_thresh = 1.0e-4
        self.max_iter = 50

        # ERI settings
        self.eri_thresh = 1.0e-15
        self.qq_type = 'QQ_DEN'

        # mpi information
        self.comm = comm
        self.rank = self.comm.Get_rank()
        self.nodes = self.comm.Get_size()

        # output stream
        self.ostream = ostream

    def update_settings(self, rsp_input):

        # calculation type
        if 'nstates' in rsp_input:
            self.nstates = int(rsp_input['nstates'])
        if 'tamm_dancoff' in rsp_input:
            key = rsp_input['tamm_dancoff'].lower()
            self.tamm_dancoff = True if key == 'yes' else False
        if 'spin' in rsp_input:
            key = rsp_input['spin'].lower()
            self.triplet = True if key[0] == 't' else False

        # solver settings
        if 'conv_thresh' in rsp_input:
            self.conv_thresh = float(rsp_input['conv_thresh'])
        if 'max_iter' in rsp_input:
            self.max_iter = int(rsp_input['max_iter'])

        # ERI settings
        if 'eri_thresh' in rsp_input:
            self.eri_thresh = float(rsp_input['eri_thresh'])
        if 'qq_type' in rsp_input:
            self.qq_type = rsp_input['qq_type'].upper()

        # properties
        if rsp_input['property'].lower() == 'polarizability':
            self.prop_type = 'POLARIZABILITY'
            self.a_ops, self.b_ops = rsp_input['operators']
            self.frequencies = rsp_input['frequencies']

        elif rsp_input['property'].lower() == 'absorption':
            self.prop_type = 'EXCITATION'

    def compute(self, molecule, ao_basis, scf_tensors):
        """Performs molecular property calculation.

        Performs molecular property calculation using molecular data

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
            self.print_header()

        # Linear response eigensolver

        if self.prop_type.upper() in ['EXCITATION']:
            if not self.tamm_dancoff:
                eigensolver = LinearResponseEigenSolver(self.comm, self.ostream)
                assert_msg_critical(
                    not self.triplet,
                    'LR EigenSolver: not yet implemented for triplets')
            else:
                eigensolver = TDAExciDriver(self.comm, self.ostream)

            eigensolver.update_settings({
                'nstates': self.nstates,
                'eri_thresh': self.eri_thresh,
                'qq_type': self.qq_type,
                'conv_thresh': self.conv_thresh,
                'max_iter': self.max_iter
            })

            return eigensolver.compute(molecule, ao_basis, scf_tensors)

        # Linear response solver

        if self.prop_type.upper() in ['POLARIZABILITY']:
            lr_solver = LinearResponseSolver(self.comm, self.ostream)

            lr_solver.update_settings({
                'a_ops': self.a_ops,
                'b_ops': self.b_ops,
                'frequencies': self.frequencies,
                'eri_thresh': self.eri_thresh,
                'qq_type': self.qq_type,
                'conv_thresh': self.conv_thresh,
                'max_iter': self.max_iter
            })

            return lr_solver.compute(molecule, ao_basis, scf_tensors)

    def print_header(self):
        """Prints response driver setup header to output stream.

        Prints molecular property calculation setup details to output stream.
        """

        self.ostream.print_blank()
        self.ostream.print_header('Response Driver Setup')
        self.ostream.print_header(23 * '=')
        self.ostream.print_blank()

        str_width = 60

        cur_str = 'Molecular Property Type   : ' + self.prop_str()
        self.ostream.print_header(cur_str.ljust(str_width))

        if self.prop_type in ['EXCITATION']:
            if self.tamm_dancoff:
                cur_str = "Response Equations Type   : Tamm-Dancoff"
                self.ostream.print_header(cur_str.ljust(str_width))
            cur_str = 'Number Of Excited States  : ' + str(self.nstates)
            self.ostream.print_header(cur_str.ljust(str_width))

        cur_str = 'Max. Number Of Iterations : ' + str(self.max_iter)
        self.ostream.print_header(cur_str.ljust(str_width))
        cur_str = 'Convergence Threshold     : ' + \
            '{:.1e}'.format(self.conv_thresh)
        self.ostream.print_header(cur_str.ljust(str_width))

        cur_str = 'ERI screening scheme      : ' + get_qq_type(self.qq_type)
        self.ostream.print_header(cur_str.ljust(str_width))
        cur_str = 'ERI Screening Threshold   : ' + \
            '{:.1e}'.format(self.eri_thresh)
        self.ostream.print_header(cur_str.ljust(str_width))
        self.ostream.print_blank()

        self.ostream.flush()

    def prop_str(self):
        """Gets string with type of molecular property calculation.

        Gets string with type of molecular property calculation (Excited
        states, linear and non-linear spectroscopies).

        Returns
        -------
        The string with type of molecular property calculation.
        """

        if self.prop_type == 'POLARIZABILITY':
            return 'Polarizability'

        if self.prop_type == 'EXCITATION':
            if not self.triplet:
                return 'Singlet Excited States'
            else:
                return "Triplet Excited States"

        return 'Undefined'
