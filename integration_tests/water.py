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

output_stream = COutputStream("water.out")
input_stream  = CInputStream("water.inp", output_stream)

# read input file and parse input data

input_stream.read(input_data, output_stream)
xyz_reader.parse(molecule, input_data, output_stream)
env_reader.parse(input_data, output_stream)
basis_reader.parse(input_data, output_stream)

# write molecular geometry

molecule.print_geometry(output_stream)
output_stream.flush()

# basis set objects

path_to_basis_sets = env_reader.get_path_to_basis_sets()

ao_basis  = basis_reader.get_ao_basis (path_to_basis_sets, molecule, output_stream)
min_basis = basis_reader.get_min_basis(path_to_basis_sets, molecule, output_stream)

print("AO basis set:", ao_basis.get_label())
print("Mininal basis set:", min_basis.get_label())
print()

# compute overlap

overlap_driver = COverlapIntegralsDriver.create(rank, size, comm)

overlap_ao  = overlap_driver.compute(molecule, ao_basis, output_stream, comm)
overlap_min = overlap_driver.compute(molecule, min_basis, output_stream, comm)
overlap_mix = overlap_driver.compute(molecule, ao_basis, min_basis, output_stream, comm)

print("Overlap of AO basis\n", overlap_ao)
print("Overlap of minimal basis\n", overlap_min)
print("Overlap of mixed basis\n", overlap_mix)

# check state

assert xyz_reader.get_state() == True
assert env_reader.get_state() == True
assert basis_reader.get_state() == True
