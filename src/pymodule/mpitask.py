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


class MpiTask:

    def __init__(self, input_fname, output_fname, mpi_comm):

        # mpi settings

        self.mpi_comm = mpi_comm
        self.mpi_rank = mpi_comm.Get_rank()
        self.mpi_size = mpi_comm.Get_size()

        # initialize molecule, basis set and output stream

        self.molecule = Molecule()
        self.ao_basis = MolecularBasis()
        self.min_basis = MolecularBasis()
        self.ostream = OutputStream(output_fname)

        # process input file on master node

        if (self.mpi_rank == mpi_master()):

            self.start_time = self.ostream.print_start_header(self.mpi_size)

            self.ostream.put_info("Reading input file %s..." % input_fname)

            # read input file

            input_dict = InputParser(input_fname).parse()

            self.ostream.put_info(
                "Found %d control groups." % len(input_dict.keys()))
            self.ostream.put_info("...done.")
            self.ostream.new_line()

            # create molecule

            self.ostream.put_info("Parsing @self.molecule group...")
            self.ostream.put_info("...done.")
            self.ostream.new_line()

            self.molecule = InputParser.create_molecule(input_dict)

            self.molecule.check_proximity(0.1, self.ostream)
            self.molecule.print_geometry(self.ostream)

            # create basis set (new code)

            basis_path = input_dict["method_settings"]["basis_path"]
            basis_label = input_dict["method_settings"]["basis"].upper()
            basis_fname = basis_path + '/' + basis_label
            basis_dict = InputParser(basis_fname).parse()

            assert_msg_critical(
                basis_label == basis_dict['basis_set_name'].upper(),
                "basis set name")

            mol_basis = InputParser.create_basis_set(self.molecule, basis_dict)

            mol_basis.print_basis("Atomic Basis", self.molecule, self.ostream)

            min_basis_label = 'MIN-CC-PVDZ'
            min_basis_fname = basis_path + '/' + min_basis_label
            min_basis_dict = InputParser(min_basis_fname).parse()

            min_basis = InputParser.create_basis_set(self.molecule, min_basis_dict)

            # create basis set (old code)

            self.ostream.put_info("Parsing @method settings group...")
            self.ostream.put_info("...done.")
            self.ostream.new_line()

            self.ao_basis = MolecularBasis.from_lib(
                input_dict["method_settings"]["basis"],
                input_dict["method_settings"]["basis_path"], self.molecule,
                self.ostream)
            self.ao_basis.print_basis("Atomic Basis", self.molecule,
                                      self.ostream)

            self.min_basis = MolecularBasis.from_lib(
                "MIN-CC-PVDZ", input_dict["method_settings"]["basis_path"],
                self.molecule, self.ostream)

            self.ostream.flush()

            assert (self.ao_basis == mol_basis)

            assert (self.min_basis == min_basis)

        # broadcast molecule and basis set

        self.molecule.broadcast(self.mpi_rank, self.mpi_comm)
        self.ao_basis.broadcast(self.mpi_rank, self.mpi_comm)
        self.min_basis.broadcast(self.mpi_rank, self.mpi_comm)

    def finish(self):

        if (self.mpi_rank == mpi_master()):
            self.ostream.print_finish_header(self.start_time)
