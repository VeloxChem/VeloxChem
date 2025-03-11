from mpi4py import MPI
import numpy as np

from .veloxchemlib import RIFockDriver
from .veloxchemlib import partition_atoms
from .veloxchemlib import mpi_master

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
    
    if nodes > 1:
        atoms = partition_atoms(natoms, rank, nodes)
        self.prepare_buffers(molecule, basis, aux_basis, atoms)
    else:
        self.prepare_buffers(molecule, basis, aux_basis)

RIFockDriver.mpi_prepare_buffers = _RIFockDriver_mpi_prepare_buffers
