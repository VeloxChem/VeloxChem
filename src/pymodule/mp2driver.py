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
    :param conventional:
        The flag for using conventional (in-memory) AO-to-MO integral
        transformation.
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

        # use conventional (in-memory) AO-to-MO integral transformation?
        self.conventional = False

    def update_settings(self, mp2_dict):
        """
        Updates settings in MP2 driver.

        :param mp2_dict:
            The dictionary of MP2 settings.
        """

        if 'conventional' in mp2_dict:
            key = mp2_dict['conventional'].lower()
            self.conventional = True if key in ['yes', 'y'] else False

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

        if self.conventional:
            self.compute_conventional(molecule, ao_basis, mol_orbs)
        else:
            self.compute_distributed(molecule, ao_basis, mol_orbs)

    def compute_conventional(self, molecule, ao_basis, mol_orbs):
        """
        Performs conventional MP2 calculation.

        :param molecule:
            The molecule.
        :param ao_basis:
            The AO basis set.
        :param mol_orbs:
            The molecular orbitals.
        """

        moints_drv = MOIntegralsDriver(self.comm, self.ostream)

        if self.rank == mpi_master():

            orb_ene = mol_orbs.ea_to_numpy()
            nocc = molecule.number_of_alpha_electrons()
            eocc = orb_ene[:nocc]
            evir = orb_ene[nocc:]
            eab = evir.reshape(-1, 1) + evir

            self.e_mp2 = 0.0
            oovv = moints_drv.compute_in_mem(molecule, ao_basis, mol_orbs,
                                             "OOVV")
            for i in range(oovv.shape[0]):
                for j in range(oovv.shape[1]):
                    ij = oovv[i, j, :, :]
                    ij_antisym = ij - ij.T
                    denom = eocc[i] + eocc[j] - eab
                    self.e_mp2 += np.sum(ij * (ij + ij_antisym) / denom)

            mp2_str = '*** MP2 correlation energy: %20.12f a.u.' % self.e_mp2
            self.ostream.print_header(mp2_str.ljust(92))
            self.ostream.print_blank()

    def compute_distributed(self, molecule, ao_basis, mol_orbs):
        """
        Performs MP2 calculation via distributed Fock builds.

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
