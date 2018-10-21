from .VeloxChemLib import InputData
from .VeloxChemLib import InputStream
from .VeloxChemLib import OutputStream
from .VeloxChemLib import MolXYZReader
from .VeloxChemLib import EnvironmentReader
from .VeloxChemLib import BasisReader
from .VeloxChemLib import Molecule
from .VeloxChemLib import MolecularBasis
from .VeloxChemLib import mpi_master


class LocalTask:

    def __init__(self, input_file, output_file):

        # input/output filenames

        self.input_file = input_file
        self.output_file = output_file

        # input data, readers, and molecule objects

        self.input_data = InputData()
        self.xyz_reader = MolXYZReader()
        self.env_reader = EnvironmentReader()
        self.basis_reader = BasisReader()
        self.molecule = Molecule()

        # input/output stream objects

        self.ostream = OutputStream(self.output_file)
        self.input_stream = InputStream(self.input_file, self.ostream)

        # read input file and parse input data

        self.input_stream.read(self.input_data, self.ostream)
        self.xyz_reader.parse(self.molecule, self.input_data, self.ostream)
        self.env_reader.parse(self.input_data, self.ostream)
        self.basis_reader.parse(self.input_data, self.ostream)

        # write molecular geometry

        self.molecule.print_geometry(self.ostream)
        self.ostream.flush()

        # basis set objects

        self.path_to_basis_sets = self.env_reader.get_path_to_basis_sets()

        self.ao_basis = self.basis_reader.get_ao_basis(
            self.path_to_basis_sets, self.molecule, self.ostream)

        self.min_basis = self.basis_reader.get_min_basis(
            self.path_to_basis_sets, self.molecule, self.ostream)


class GlobalTask:

    def __init__(self, input_file, output_file, mpi_comm):

        # mpi settings

        self.mpi_comm = mpi_comm
        self.mpi_rank = mpi_comm.Get_rank()
        self.mpi_size = mpi_comm.Get_size()

        # process input file on master node

        if (self.mpi_rank == mpi_master()):

            task = LocalTask(input_file, output_file)
            self.molecule = task.molecule
            self.ao_basis = task.ao_basis
            self.min_basis = task.min_basis
            self.ostream = task.ostream

        else:

            self.molecule = Molecule()
            self.ao_basis = MolecularBasis()
            self.min_basis = MolecularBasis()
            self.ostream = OutputStream("")

        # broadcast molecule and basis

        self.molecule.broadcast(self.mpi_rank, self.mpi_comm)
        self.ao_basis.broadcast(self.mpi_rank, self.mpi_comm)
        self.min_basis.broadcast(self.mpi_rank, self.mpi_comm)
