from veloxchem.VeloxChemLib import InputData
from veloxchem.VeloxChemLib import MolXYZReader
from veloxchem.VeloxChemLib import EnvironmentReader
from veloxchem.VeloxChemLib import BasisReader
from veloxchem.VeloxChemLib import Molecule
from veloxchem.VeloxChemLib import OutputStream
from veloxchem.VeloxChemLib import InputStream


class Task(object):

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
