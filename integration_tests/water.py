from mpi4py import MPI
from VeloxChemMP import *
from task import *
import numpy as np

# mpi settings

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

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

else:

    molecule = Molecule()
    ao_basis = MolecularBasis()
    min_basis = MolecularBasis()
    ostream = OutputStream("")

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

overlap = S22.to_numpy()
density = D.total_to_numpy(0)

if (rank == mpi_master()):

    # get attributes

    print()
    print("Dimension of density matrix:", end=' ')
    for i in range(density.ndim):
        print(density.shape[i], end=', '),
    print('\n')

    # get number of electrons

    DS = density.dot(overlap)

    nelec = 0.0
    for i in range(DS.shape[0]):
        nelec += DS[i][i]
    nelec *= 2.0

    print("Number of electrons:", nelec)
    print()

    # check DSD-D

    DSD = DS.dot(density)

    print("DSD-D MaxDiff =", np.max(np.abs(DSD - density)))
    print()

# numpy to matrix

S22new   = OverlapMatrix.from_numpy(overlap)
D_rest   = AODensityMatrix.from_numpy_list([density], True)
D_unrest = AODensityMatrix.from_numpy_list([density, density], False)

if (rank == mpi_master()):

    assert(S22new == S22)

    assert(D_rest == D)

# flush output stream

if (rank == mpi_master()):

    ostream.flush()
