from mpi4py import MPI
from VeloxChemMP import *
from task import *

# mpi settings

comm = MPI.COMM_WORLD
rank, size = comm.Get_rank(), comm.Get_size()

# initialize mandatory objects

molecule = CMolecule()
ao_basis = CMolecularBasis()
min_basis = CMolecularBasis()
ostream = COutputStream("dummy.out")

# process input file on master node

if (rank == mpi_master()):
    task = Task("water.inp", "water.out")
    molecule = task.molecule
    ao_basis = task.ao_basis
    min_basis = task.min_basis
    ostream = task.ostream
    ostream.put_info("Molecular basis set: %s" % ao_basis.get_label())
    ostream.put_info("Minimal basis set: %s" % min_basis.get_label())
    ostream.new_line()

# broadcast molecule and basis

molecule.broadcast(rank, comm)
ao_basis.broadcast(rank, comm)
min_basis.broadcast(rank, comm)

# compute overlap

overlap_driver = COverlapIntegralsDriver.create(rank, size, comm)

S12 = overlap_driver.compute(molecule, min_basis, ao_basis, ostream, comm)
S22 = overlap_driver.compute(molecule, ao_basis, ostream, comm)

sad_driver = CSADGuessDriver.create(rank, size, comm)
density_mat = sad_driver.compute(molecule, min_basis, ao_basis, S12, S22, ostream, comm)

if (rank == mpi_master()):
    print("SAD guess\n%s" % density_mat)
    ostream.flush()
