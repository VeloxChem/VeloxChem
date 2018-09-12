from mpi4py import MPI
from VeloxChemMP import *

# mpi settings

comm = MPI.COMM_WORLD
rank, size = comm.Get_rank(), comm.Get_size()

# input data, readers, and molecule objects

input_data   = CInputData()
xyz_reader   = CMolXYZReader()
env_reader   = CEnvironmentReader()
basis_reader = CBasisReader()
molecule     = CMolecule()

# input/output stream objects

output_stream = COutputStream("dimer.out")
input_stream  = CInputStream("dimer.inp", output_stream)

# read input file and parse input data

input_stream.read(input_data, output_stream)
xyz_reader.parse(molecule, input_data, output_stream)
env_reader.parse(input_data, output_stream)
basis_reader.parse(input_data, output_stream)

# write molecular geometry

molecule.print_geometry(output_stream)
output_stream.flush()

# sub-molecules (fragments/monomers)
# two parameters: starting index and number of atoms

mol_1 = molecule.get_sub_molecule(0,4)
mol_2 = molecule.get_sub_molecule(4,5)

mol_1.print_geometry(output_stream)
mol_2.print_geometry(output_stream)
output_stream.flush()

# basis set objects

path_to_basis_sets = env_reader.get_path_to_basis_sets()

basis = basis_reader.get_min_basis(path_to_basis_sets, molecule, output_stream)

# compute overlap

overlap_driver = COverlapIntegralsDriver.create(rank, size, comm)

S = overlap_driver.compute(molecule, basis, output_stream, comm)

S11 = overlap_driver.compute(mol_1, basis, output_stream, comm)
S22 = overlap_driver.compute(mol_2, basis, output_stream, comm)
S21 = overlap_driver.compute(mol_2, mol_1, basis, output_stream, comm)
S12 = overlap_driver.compute(mol_1, mol_2, basis, output_stream, comm)

S_new = assemble_overlap_matrices(mol_1, mol_2, basis, basis, S11, S22, S12, S21)

print("Overlap \n", S)
print("Overlap new\n", S_new)

print(S == S_new)

assert(S == S_new)

# check state

assert xyz_reader.get_state() == True
assert env_reader.get_state() == True
assert basis_reader.get_state() == True
