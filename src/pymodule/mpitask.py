from .veloxchemlib import Molecule
from .veloxchemlib import MolecularBasis
from .veloxchemlib import AtomBasis
from .veloxchemlib import BasisFunction
from .veloxchemlib import ChemicalElement
from .veloxchemlib import OutputStream
from .veloxchemlib import mpi_master
from .veloxchemlib import assert_msg_critical
from .veloxchemlib import to_angular_momentum
from .inputparser import InputParser

from os.path import isfile


class MpiTask:

    def __init__(self, fname_list, mpi_comm):

        # mpi settings

        self.mpi_comm = mpi_comm
        self.mpi_rank = mpi_comm.Get_rank()
        self.mpi_size = mpi_comm.Get_size()

        # input/output files

        input_fname = ''
        output_fname = ''

        if self.mpi_rank == mpi_master():

            assert_msg_critical(
                len(fname_list) >= 2,
                "MpiTask: Need input and output file names")

            input_fname = fname_list[0]
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

            self.ostream.put_info("Reading input file %s..." % input_fname)

            # read input file

            input_parser = InputParser(input_fname)
            input_dict = input_parser.get_dict()

            self.input_dict = input_dict

            self.ostream.put_info(
                "Found %d control groups." % len(input_dict.keys()))
            self.ostream.put_info("...done.")
            self.ostream.new_line()

            # create molecule

            self.ostream.put_info("Parsing @molecule group...")
            self.ostream.put_info("...done.")
            self.ostream.new_line()

            self.molecule = input_parser.create_molecule()

            self.molecule.check_proximity(0.1, self.ostream)
            self.molecule.print_geometry(self.ostream)

            # create basis set

            self.ostream.put_info("Parsing @method settings group...")
            self.ostream.put_info("...done.")
            self.ostream.new_line()

            basis_path = input_dict["method_settings"]["basis_path"]
            basis_label = input_dict["method_settings"]["basis"].upper()
            basis_fname = basis_path + '/' + basis_label

            basis_parser = InputParser(basis_fname)
            basis_dict = basis_parser.get_dict()

            assert_msg_critical(
                basis_label == basis_dict['basis_set_name'].upper(),
                "basis set name")

            self.ao_basis = basis_parser.create_basis_set(self.molecule)

            self.ao_basis.print_basis("Atomic Basis", self.molecule,
                                      self.ostream)

            min_basis_fname = basis_path + '/MIN-CC-PVDZ'
            min_basis_parser = InputParser(min_basis_fname)

            self.min_basis = min_basis_parser.create_basis_set(self.molecule)

            self.ostream.flush()

        # broadcast molecule and basis set

        self.molecule.broadcast(self.mpi_rank, self.mpi_comm)
        self.ao_basis.broadcast(self.mpi_rank, self.mpi_comm)
        self.min_basis.broadcast(self.mpi_rank, self.mpi_comm)

    def finish(self):

        if (self.mpi_rank == mpi_master()):
            self.ostream.print_finish_header(self.start_time)
