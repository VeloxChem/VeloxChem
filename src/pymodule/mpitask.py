from .veloxchemlib import Molecule
from .veloxchemlib import MolecularBasis
from .veloxchemlib import OutputStream
from .veloxchemlib import mpi_master
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

            input_dict = InputParser(input_fname).parse()

            self.ostream.put_info(
                "Found %d control groups." % len(input_dict.keys()))
            self.ostream.put_info("...done.")
            self.ostream.new_line()

            self.ostream.put_info("Parsing @self.molecule group...")
            self.ostream.put_info("...done.")
            self.ostream.new_line()

            self.molecule = Molecule.from_xyz(
                input_dict["molecule"]["atom_labels"],
                input_dict["molecule"]["x_coords"],
                input_dict["molecule"]["y_coords"],
                input_dict["molecule"]["z_coords"])

            if "charge" in input_dict["molecule"].keys():
                self.molecule.set_charge(int(input_dict["molecule"]["charge"]))

            if "multiplicity" in input_dict["molecule"].keys():
                self.molecule.set_multiplicity(
                    int(input_dict["molecule"]["multiplicity"]))

            self.molecule.check_multiplicity()
            self.molecule.check_proximity(0.1, self.ostream)
            self.molecule.print_geometry(self.ostream)

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

        # broadcast molecule and basis set

        self.molecule.broadcast(self.mpi_rank, self.mpi_comm)
        self.ao_basis.broadcast(self.mpi_rank, self.mpi_comm)
        self.min_basis.broadcast(self.mpi_rank, self.mpi_comm)

    def finish(self):

        if (self.mpi_rank == mpi_master()):
            self.ostream.print_finish_header(self.start_time)
