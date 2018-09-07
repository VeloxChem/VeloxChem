from VeloxChemMP import *
import __future__

# objects from default constructors

input_data   = CInputData()
xyz_reader   = CMolXYZReader()
env_reader   = CEnvironmentReader()
basis_reader = CBasisReader()
molecule     = CMolecule()

# objects from constructors with arguments

output_stream = COutputStream("test.out")
input_stream  = CInputStream("test.inp", output_stream)

# call methods of objects

input_stream.read(input_data, output_stream)
xyz_reader.parse(molecule, input_data, output_stream)
env_reader.parse(input_data, output_stream)
basis_reader.parse(input_data, output_stream)

path_to_basis_sets = env_reader.get_path_to_basis_sets()

molecule.print_geometry(output_stream)
output_stream.flush()

ao_basis = basis_reader.get_ao_basis(path_to_basis_sets, molecule, output_stream)
print("AO basis set:", ao_basis.get_label())

min_basis = basis_reader.get_min_basis(path_to_basis_sets, molecule, output_stream)
print("Mininal basis set:", min_basis.get_label())

# check state

assert xyz_reader.get_state() == True
assert env_reader.get_state() == True
assert basis_reader.get_state() == True
