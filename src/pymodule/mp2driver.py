import numpy as np
import sys

from .veloxchemlib import mpi_master
from .mointsdriver import MOIntegralsDriver
from .outputstream import OutputStream
from .mpiutils import split_comm


class Mp2Driver:

    def __init__(self):

        self.e_mp2 = None

    def compute_task(self, task, mol_orbs):

        self.compute(task.molecule, task.ao_basis, mol_orbs, task.mpi_comm,
                     task.ostream)

    def compute(self,
                molecule,
                ao_basis,
                mol_orbs,
                comm,
                ostream=OutputStream(sys.stdout)):

        # MO integrals
        moints_drv = MOIntegralsDriver()

        grps = [p for p in range(comm.Get_size())]
        oovv = moints_drv.compute(molecule, ao_basis, mol_orbs, 'OOVV', grps,
                                  comm, ostream)

        # MP2 energy
        local_comm, cross_comm = split_comm(comm, grps)
        local_master = (local_comm.Get_rank() == mpi_master())
        global_master = (comm.Get_rank() == mpi_master())

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
            ostream.print_header(mp2_str.ljust(92))
            ostream.print_blank()

        local_comm.Disconnect()
        cross_comm.Disconnect()
