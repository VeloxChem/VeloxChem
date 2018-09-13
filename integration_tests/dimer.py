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

basis = basis_reader.get_ao_basis(path_to_basis_sets, molecule, output_stream)

basis_1 = basis_reader.get_ao_basis(path_to_basis_sets, mol_1, output_stream)
basis_2 = basis_reader.get_ao_basis(path_to_basis_sets, mol_2, output_stream)

# compute overlap

overlap_driver = COverlapIntegralsDriver.create(rank, size, comm)

S = overlap_driver.compute(molecule, basis, output_stream, comm)

S11 = overlap_driver.compute(mol_1, basis, output_stream, comm)
S22 = overlap_driver.compute(mol_2, basis, output_stream, comm)

S12 = overlap_driver.compute(mol_1, mol_2, basis, output_stream, comm)
S21 = overlap_driver.compute(mol_2, mol_1, basis, output_stream, comm)

S12s = overlap_driver.compute(mol_1, mol_2, basis_1, basis_2, output_stream, comm)
S21s = overlap_driver.compute(mol_2, mol_1, basis_2, basis_1, output_stream, comm)
assert(S12s == S12)
assert(S21s == S21)

S_new = assemble_overlap_matrices(mol_1, mol_2, basis, basis, S11, S22, S12, S21)
assert(S == S_new)
print("Done: Overlap integral")

# compute kinetic energy

kinetic_energy_driver = CKineticEnergyIntegralsDriver.create(rank, size, comm)

T = kinetic_energy_driver.compute(molecule, basis, output_stream, comm)

T11 = kinetic_energy_driver.compute(mol_1, basis, output_stream, comm)
T22 = kinetic_energy_driver.compute(mol_2, basis, output_stream, comm)

T12 = kinetic_energy_driver.compute(mol_1, mol_2, basis, output_stream, comm)
T21 = kinetic_energy_driver.compute(mol_2, mol_1, basis, output_stream, comm)

T12s = kinetic_energy_driver.compute(mol_1, mol_2, basis_1, basis_2, output_stream, comm)
T21s = kinetic_energy_driver.compute(mol_2, mol_1, basis_2, basis_1, output_stream, comm)
assert(T12s == T12)
assert(T21s == T21)

T_new = assemble_kinetic_energy_matrices(mol_1, mol_2, basis, basis, T11, T22, T12, T21)
assert(T == T_new)
print("Done: Kinetic energy integral")

# compute nuclear potential

nuclear_potential_driver = CNuclearPotentialIntegralsDriver.create(rank, size, comm)

V = nuclear_potential_driver.compute(molecule, basis, output_stream, comm)

V11 = nuclear_potential_driver.compute(mol_1, basis, molecule, output_stream, comm)
V22 = nuclear_potential_driver.compute(mol_2, basis, molecule, output_stream, comm)

V12 = nuclear_potential_driver.compute(mol_1, mol_2, basis, molecule, output_stream, comm)
V21 = nuclear_potential_driver.compute(mol_2, mol_1, basis, molecule, output_stream, comm)

V12s = nuclear_potential_driver.compute(mol_1, mol_2, basis_1, basis_2, molecule, output_stream, comm)
V21s = nuclear_potential_driver.compute(mol_2, mol_1, basis_2, basis_1, molecule, output_stream, comm)
assert(V12s == V12)
assert(V21s == V21)

V_new = assemble_nuclear_potential_matrices(mol_1, mol_2, basis, basis, V11, V22, V12, V21)
assert(V == V_new)
print("Done: Nuclear potential integral")

# check state

assert xyz_reader.get_state() == True
assert env_reader.get_state() == True
assert basis_reader.get_state() == True

print("Success!")
