from mpi4py import MPI

from .veloxchemlib import FockDriver
from .matrix import Matrix

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
    
FockDriver.mpi_compute = _FockDriver_mpi_compute
