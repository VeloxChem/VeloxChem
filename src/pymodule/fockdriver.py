from mpi4py import MPI

from .veloxchemlib import FockDriver
from .matrix import Matrix
from .matrices import Matrices

def _FockDriver_mpi_compute(self, comm, screener, density, label, exchange_factor, omega, ithreshold):
    """
    Computes Fock matrix for given density matrix.

    :param comm:
        The MPI communicator.
    :param screener:
        The sreener with ERIs data.
    :param density:
        The density matrix to compute Fock matrix.
    :param label:
        The standard label of Fock matrix type.
    :param exchange_factor:
        The exchange scaling factor.
    :param omega:
        The range separation factor.
    :param ithreshold:
        The threshold of integrals.
    
    :return:
        Fock matrix.
    """

    # compute local Fock matrix
    loc_mat = self.compute(screener, comm.Get_rank(), comm.Get_size(), density, label, exchange_factor, omega, ithreshold)
    
    # reduce Fock matrix
    fock_mat = None
    fock_mat = Matrix.reduce(loc_mat, comm, 0)
    return fock_mat
    
def _FockDriver_mpi_multi_compute(self, comm, screener, densities, labels, exchange_factor, omega, ithreshold):
    """
    Computes Fock matrices for given density matrices.

    :param comm:
        The MPI communicator.
    :param screener:
        The sreener with ERIs data.
    :param densities:
        The density matrices to compute Fock matrices.
    :param label:
        The list of standard Fock matrix type labels.
    :param exchange_factor:
        The exchange scaling factor.
    :param omega:
        The range separation factor.
    :param ithreshold:
        The threshold of integrals.
    
    :return:
        Fock matrices.
    """

    # compute local Fock matrices
    loc_mats = self.compute(screener, comm.Get_rank(), comm.Get_size(), densities, labels, exchange_factor, omega, ithreshold)
    
    # reduce Fock matrices
    fock_mats = None
    fock_mats = Matrices.reduce(loc_mats, comm, 0)
    return fock_mats
    
FockDriver.mpi_compute = _FockDriver_mpi_compute
FockDriver.mpi_multi_compute = _FockDriver_mpi_multi_compute
