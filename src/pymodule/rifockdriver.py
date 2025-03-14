from mpi4py import MPI
import numpy as np

from .veloxchemlib import RIFockDriver
from .veloxchemlib import partition_atoms
from .veloxchemlib import mpi_master
from .veloxchemlib import Matrix

def _RIFockDriver_mpi_prepare_buffers(self, comm, molecule, basis, aux_basis):
    """
    Computes distributed three-center electron integrals tensor.

    :param comm:
        The MPI communicator.
    :param molecule:
        The molecule to compute gradient.
    :param basis:
        The basis set to compute gradient.
    :param aux_basis:
        The auxilary basis set to compute gradient.
    """

    # set up rank, node data
    
    rank = comm.Get_rank()
    nodes = comm.Get_size()

    # select three-center integrals computation method
    
    natoms = molecule.number_of_atoms()
    atoms = partition_atoms(natoms, rank, nodes)
    self.prepare_buffers(molecule, basis, aux_basis, atoms)
        
def _RIFockDriver_mpi_compute_bq_vector(self, comm, density):
    """
    Computes distributed Bq vector.

    :param comm:
        The MPI communicator.
    :param density:
        The density to compute Bq vector.
        
    :return:
        The Bq vector.
    """

    # set up rank, node data
    gv = np.array(self.compute_local_bq_vector(density))
    tv = np.zeros(gv.shape)
    comm.Allreduce([gv, MPI.DOUBLE], [tv, MPI.DOUBLE], op=MPI.SUM)
    return list(tv)

def _RIFockDriver_mpi_compute(self, comm, gamma, density, label):
    """
    Computes Fock matrix.

    :param comm:
        The MPI communicator.
    :param gamma:
        The Bq vector.    
    :param density:
        The density to compute Fock matrix.
    :param label:
        The type of Fock matrix.
        
    :return:
        The Fock matrix.
    """

    # set up rank, node data
    fmat = self.local_compute(density, gamma, label)
    rfmat = Matrix.reduce(fmat, comm, mpi_master())
    return rfmat

RIFockDriver.mpi_prepare_buffers = _RIFockDriver_mpi_prepare_buffers
RIFockDriver.mpi_compute_bq_vector = _RIFockDriver_mpi_compute_bq_vector
RIFockDriver.mpi_compute = _RIFockDriver_mpi_compute
