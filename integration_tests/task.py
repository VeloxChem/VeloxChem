from VeloxChemMP import *

class Task(object):

    def __init__(self, input_file, output_file):

        # input/output filenames

        self.input_file = input_file
        self.output_file = output_file

        # input data, readers, and molecule objects

        self.input_data   = CInputData()
        self.xyz_reader   = CMolXYZReader()
        self.env_reader   = CEnvironmentReader()
        self.basis_reader = CBasisReader()
        self.molecule     = CMolecule()

        # input/output stream objects

        self.output_stream = COutputStream(self.output_file)
        self.input_stream  = CInputStream(self.input_file, self.output_stream)

        # read input file and parse input data

        self.input_stream.read(self.input_data, self.output_stream)
        self.xyz_reader.parse(self.molecule, self.input_data, self.output_stream)
        self.env_reader.parse(self.input_data, self.output_stream)
        self.basis_reader.parse(self.input_data, self.output_stream)

        # write molecular geometry

        self.molecule.print_geometry(self.output_stream)

        # basis set objects

        self.path_to_basis_sets = self.env_reader.get_path_to_basis_sets()

        self.ao_basis = self.basis_reader.get_ao_basis(self.path_to_basis_sets,
                                                       self.molecule,
                                                       self.output_stream)

        self.min_basis = self.basis_reader.get_min_basis(self.path_to_basis_sets,
                                                         self.molecule,
                                                         self.output_stream)

        self.output_stream.flush()

        assert(self.xyz_reader.get_state() == True)
        assert(self.env_reader.get_state() == True)
        assert(self.basis_reader.get_state() == True)
