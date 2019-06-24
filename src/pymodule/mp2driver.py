import numpy as np

from .veloxchemlib import mpi_master
from .mointsdriver import MOIntegralsDriver
from .subcommunicators import SubCommunicators


class Mp2Driver:
    """
    Implements MP2 driver.

    :param e_mp2:
        The MP2 correlation energy.
    :param comm:
        The MPI communicator.
    :param rank:
        The MPI rank.
    :param nodes:
        Number of MPI processes.
    :param ostream:
        The output stream.
    """

    def __init__(self, comm, ostream):
        """
        Initializes MP2 driver.

        :param comm:
            The MPI communicator.
        :param ostream:
            The output stream.
        """

        self.e_mp2 = None

        # mpi information
        self.comm = comm
        self.rank = self.comm.Get_rank()
        self.nodes = self.comm.Get_size()

        # output stream
        self.ostream = ostream

    def compute(self, molecule, ao_basis, mol_orbs):
        """
        Performs MP2 calculation.

        :param molecule:
            The molecule.
        :param ao_basis:
            The AO basis set.
        :param mol_orbs:
            The molecular orbitals.
        """

        # MO integrals
        moints_drv = MOIntegralsDriver(self.comm, self.ostream)

        grps = [p for p in range(self.nodes)]
        oovv = moints_drv.compute(molecule, ao_basis, mol_orbs, 'OOVV', grps)

        # MP2 energy
        subcomm = SubCommunicators(self.comm, grps)
        local_comm = subcomm.local_comm
        cross_comm = subcomm.cross_comm

        local_master = (local_comm.Get_rank() == mpi_master())
        global_master = (self.rank == mpi_master())

        if local_master:
            orb_ene = mol_orbs.ea_to_numpy()
            nocc = molecule.number_of_alpha_electrons()
            evir = orb_ene[nocc:]
            eab = evir.reshape(-1, 1) + evir

            local_e_mp2 = 0.0

            for pair in oovv.get_gen_pairs():
                eij = orb_ene[pair.first()] + orb_ene[pair.second()]
                denom = eij - eab

                ij = oovv.to_numpy(pair)
                ij_asym = ij - ij.T

                local_e_mp2 += np.sum(ij * (ij + ij_asym) / denom)

            local_e_mp2 = cross_comm.gather(local_e_mp2, root=mpi_master())

        if global_master:
            self.e_mp2 = sum(local_e_mp2)
            mp2_str = '*** MP2 correlation energy: %20.12f a.u.' % self.e_mp2
            self.ostream.print_header(mp2_str.ljust(92))
            self.ostream.print_blank()
