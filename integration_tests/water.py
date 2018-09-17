from mpi4py import MPI
from VeloxChemMP import *
from task import *
import numpy as np

# mpi settings

comm = MPI.COMM_WORLD
rank, size = comm.Get_rank(), comm.Get_size()

# initialize mandatory objects

molecule = Molecule()
ao_basis = MolecularBasis()
min_basis = MolecularBasis()
ostream = OutputStream("dummy.out")

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

overlap_driver = OverlapIntegralsDriver.create(rank, size, comm)

S12 = overlap_driver.compute(molecule, min_basis, ao_basis, ostream, comm)

S22 = overlap_driver.compute(molecule, ao_basis, ostream, comm)

# compute initial guess

sad_driver = SADGuessDriver.create(rank, size, comm)

D = sad_driver.compute(molecule, min_basis, ao_basis, S12, S22, ostream, comm)

# matrix to numpy

overlap = to_numpy(S22)

density = to_numpy(D)

if (rank == mpi_master()):

    # get attributes

    print()
    print("The dimension of the density matrix is:", end=' ')
    for i in range(density.ndim):
        print(density.shape[i], end=', '),
    print('\n')

    # get number of electrons

    DS = density.dot(overlap)

    nelec = 0.0
    for i in range(DS.shape[0]):
        nelec += DS[i][i]
    nelec *= 2.0

    print("The number of electrons is:", nelec)
    print()

# numpy to matrix

S22new = OverlapMatrix.from_numpy(overlap)

D_new = sad_driver.compute(molecule, min_basis, ao_basis, S12, S22new, ostream, comm);

if (rank == mpi_master()):

    print("Difference =", np.max(np.abs(to_numpy(D) - to_numpy(D_new))))
    print()

# flush output stream

if (rank == mpi_master()):

    ostream.flush()
