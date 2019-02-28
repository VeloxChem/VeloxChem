from .veloxchemlib import AtomBasis
from .veloxchemlib import BasisFunction
from .veloxchemlib import ChemicalElement
from .veloxchemlib import mpi_master
from .veloxchemlib import assert_msg_critical
from .veloxchemlib import to_angular_momentum
from .inputparser import InputParser
from .outputstream import OutputStream
from .molecule import Molecule
from .molecularbasis import MolecularBasis

from os.path import isfile


class MpiTask:

    def __init__(self, fname_list, mpi_comm):

        # mpi settings

        self.mpi_comm = mpi_comm
        self.mpi_rank = mpi_comm.Get_rank()
        self.mpi_size = mpi_comm.Get_size()

        # input/output files
        # on master node:  output_fname is string:    ostream is file handle
        #                  output_fname is "" or "-": ostream is sys.stdout
        # on worker nodes: output_fname is None:      ostream is None

        input_fname = None
        output_fname = None

        if self.mpi_rank == mpi_master():

            assert_msg_critical(
                len(fname_list) >= 1,
                "MpiTask: Need input file name")

            input_fname = fname_list[0]

            output_fname = ""
            if len(fname_list) >= 2:
                output_fname = fname_list[1]

            assert_msg_critical(
                isfile(input_fname),
                "MpiTask: input file %s does not exist" % input_fname)

            assert_msg_critical(
                input_fname != output_fname,
                "MpiTask: input/output file cannot be the same")

        # initialize molecule, basis set and output stream

        self.molecule = Molecule()
        self.ao_basis = MolecularBasis()
        self.min_basis = MolecularBasis()
        self.ostream = OutputStream(output_fname)

        # process input file on master node

        self.input_dict = {}

        if self.mpi_rank == mpi_master():

            self.start_time = self.ostream.print_start_header(self.mpi_size)

            self.ostream.print_info("Reading input file %s..." % input_fname)

            # read input file

            self.input_dict = InputParser(input_fname).get_dict()

            self.ostream.print_info(
                "Found %d control groups." % len(self.input_dict.keys()))
            self.ostream.print_info("...done.")
            self.ostream.print_blank()

            # create molecule

            self.ostream.print_info("Parsing @molecule group...")
            self.ostream.print_info("...done.")
            self.ostream.print_blank()

            self.molecule = Molecule.from_dict(self.input_dict['molecule'])

            self.ostream.print_block(self.molecule.get_string())

            # create basis set

            self.ostream.print_info("Parsing @method settings group...")
            self.ostream.print_info("...done.")
            self.ostream.print_blank()

            basis_path = self.input_dict["method_settings"]["basis_path"]
            basis_name = self.input_dict["method_settings"]["basis"].upper()

            self.ao_basis = MolecularBasis.read(
                self.molecule, basis_name, basis_path)

            self.min_basis = MolecularBasis.read(
                self.molecule, "MIN-CC-PVDZ", basis_path)

            self.ostream.print_block(
                self.ao_basis.get_string("Atomic Basis", self.molecule))

            self.ostream.flush()

        # broadcast input dictionary

        self.input_dict = self.mpi_comm.bcast(
            self.input_dict, root=mpi_master())

        # broadcast molecule and basis set

        self.molecule.broadcast(self.mpi_rank, self.mpi_comm)
        self.ao_basis.broadcast(self.mpi_rank, self.mpi_comm)
        self.min_basis.broadcast(self.mpi_rank, self.mpi_comm)

    def finish(self):

        if (self.mpi_rank == mpi_master()):
            self.ostream.print_finish_header(self.start_time)
