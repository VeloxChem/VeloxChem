import time as tm
import sys

from .veloxchemlib import ElectronRepulsionIntegralsDriver
from .veloxchemlib import mpi_master
from .veloxchemlib import fockmat
from .veloxchemlib import moints
from .veloxchemlib import MOIntsBatch
from .veloxchemlib import TwoIndexes
from .outputstream import OutputStream
from .aofockmatrix import AOFockMatrix
from .aodensitymatrix import AODensityMatrix
from .subcommunicators import SubCommunicators
from .qqscheme import get_qq_type
from .qqscheme import get_qq_scheme


class MOIntegralsDriver:

    def __init__(self):

        # screening scheme
        self.qq_type = "QQ_DEN"
        self.eri_thresh = 1.0e-12

        # Fock matrices
        self.num_matrices = 0
        self.batch_size = 3000

    def compute_task(self, task, mol_orbs, mints_type, grps):

        return self.compute(task.molecule, task.ao_basis, mol_orbs, mints_type,
                            grps, task.mpi_comm, task.ostream)

    def compute(self,
                molecule,
                ao_basis,
                mol_orbs,
                mints_type,
                grps,
                global_comm,
                ostream=OutputStream(sys.stdout)):

        # start timer
        start_time = tm.time()

        # split communicators
        subcomm = SubCommunicators(global_comm, grps)
        local_comm = subcomm.local_comm
        cross_comm = subcomm.cross_comm

        global_rank = global_comm.Get_rank()
        global_master = (global_rank == mpi_master())

        local_rank = local_comm.Get_rank()
        local_nodes = local_comm.Get_size()
        local_master = (local_rank == mpi_master())

        cross_rank = cross_comm.Get_rank()
        cross_nodes = cross_comm.Get_size()

        # set up screening data
        eri_drv = ElectronRepulsionIntegralsDriver(local_rank, local_nodes,
                                                   local_comm)

        qq_data = eri_drv.compute(
            get_qq_scheme(self.qq_type), self.eri_thresh, molecule, ao_basis)

        # initialize MO integrals batch
        moints_batch = MOIntsBatch()

        # set up properties of MO integrals batch
        moints_batch.set_ext_indexes(self.get_external_indexes(mints_type))
        moints_batch.set_batch_type(self.get_moints_type(mints_type))

        # broadcast molecular orbitals to local master nodes
        if local_master:
            mol_orbs.broadcast(cross_rank, cross_comm)

            # transformation matrices
            xmat, ymat = self.get_transformation_vectors(
                mints_type, mol_orbs, molecule)

            # set up bra and ket orbitals list
            bra_ids, ket_ids = self.set_orbital_pairs(mints_type, mol_orbs,
                                                      molecule)

            # create MO integrals batch
            batch_pos, batch_dim = self.get_batch_dimensions(bra_ids)

            loc_bra_ids = [
                bra_ids[pos + cross_rank:pos + dim:cross_nodes]
                for pos, dim in zip(batch_pos, batch_dim)
            ]
            loc_ket_ids = [
                ket_ids[pos + cross_rank:pos + dim:cross_nodes]
                for pos, dim in zip(batch_pos, batch_dim)
            ]

        else:
            loc_bra_ids = []
            loc_ket_ids = []

        loc_bra_ids = local_comm.bcast(loc_bra_ids, root=mpi_master())
        loc_ket_ids = local_comm.bcast(loc_ket_ids, root=mpi_master())

        # print title
        if global_master:
            self.print_header(ostream)

        # loop over batches
        for cur_bra_ids, cur_ket_ids in zip(loc_bra_ids, loc_ket_ids):

            # compute pair densities on local master nodes
            if local_master:
                pair_den = mol_orbs.get_pair_density(cur_bra_ids, cur_ket_ids)
            else:
                pair_den = AODensityMatrix()

            # broadcast pair densities via local communicators
            pair_den.broadcast(local_rank, local_comm)

            # compute multiple Fock matrices
            fock_mat = AOFockMatrix(pair_den)

            self.set_fock_matrices_type(mints_type, fock_mat)
            eri_drv.compute(fock_mat, pair_den, molecule, ao_basis, qq_data,
                            local_comm)

            fock_mat.reduce_sum(local_rank, local_nodes, local_comm)

            if local_master:
                # transform Fock matrices and add MO integrals to batch
                moints_batch.append(fock_mat, xmat, ymat, cur_bra_ids,
                                    cur_ket_ids)

        if global_master:
            self.print_finish(start_time, ostream)

        return moints_batch

    def collect_moints_batches(self, moints_batch, grps, global_comm):

        subcomm = SubCommunicators(global_comm, grps)
        local_comm = subcomm.local_comm
        cross_comm = subcomm.cross_comm

        global_master = (global_comm.Get_rank() == mpi_master())
        local_master = (local_comm.Get_rank() == mpi_master())

        cross_rank = cross_comm.Get_rank()
        cross_nodes = cross_comm.Get_size()

        if local_master:
            moints_batch.collect_batches(cross_rank, cross_nodes, cross_comm)

        if not global_master:
            moints_batch = MOIntsBatch()

        return moints_batch

    def print_header(self, ostream):

        ostream.print_blank()
        ostream.print_header("Integrals Transformation (AO->MO) Driver Setup")
        ostream.print_header(46 * "=")
        ostream.print_blank()

        str_width = 80
        cur_str = "Number of Fock matrices      : " + str(self.num_matrices)
        ostream.print_header(cur_str.ljust(str_width))
        cur_str = "Size of Fock matrices batch  : " + str(self.batch_size)
        ostream.print_header(cur_str.ljust(str_width))
        cur_str = "ERI screening scheme         : " + get_qq_type(self.qq_type)
        ostream.print_header(cur_str.ljust(str_width))
        cur_str = "ERI Screening Threshold      : " + \
            "{:.1e}".format(self.eri_thresh)
        ostream.print_header(cur_str.ljust(str_width))
        ostream.print_blank()
        ostream.flush()

    def print_finish(self, start_time, ostream):

        valstr = "*** Integrals transformation (AO->MO) done in "
        valstr += "{:.2f}".format(tm.time() - start_time) + " sec."
        ostream.print_blank()
        ostream.print_header(valstr.ljust(92))
        ostream.print_blank()

    def get_batch_dimensions(self, vec_ids):

        vdim = len(vec_ids)

        # determine number of batches
        nbatch = vdim // self.batch_size
        if (vdim % self.batch_size) != 0:
            nbatch = nbatch + 1

        # compute dimensions of batches
        bdim = vdim // nbatch
        brem = vdim % nbatch

        # populate dimensions list
        vec_dim = []
        for i in range(nbatch):
            if i < brem:
                vec_dim.append(bdim + 1)
            else:
                vec_dim.append(bdim)

        # populate positions list
        vec_pos = []
        cur_pos = 0
        for i in range(nbatch):
            vec_pos.append(cur_pos)
            cur_pos = cur_pos + vec_dim[i]

        return (vec_pos, vec_dim)

    def get_num_orbitals(self, mol_orbs, molecule):

        nocc = molecule.number_of_alpha_electrons()
        nvirt = mol_orbs.number_mos() - nocc

        return (nocc, nvirt)

    def set_orbital_pairs(self, mints_type, mol_orbs, molecule):

        nocc, nvirt = self.get_num_orbitals(mol_orbs, molecule)

        # cases: oooo, ooov, oovv, ovov
        bra_dim = (0, nocc)
        ket_dim = (0, nocc)

        # case: ovvv
        if mints_type == "OVVV":
            ket_dim = (nocc, nocc + nvirt)

        # case: vvvv
        if mints_type == "VVVV":
            bra_dim = (nocc, nocc + nvirt)
            ket_dim = (nocc, nocc + nvirt)

        # set up list of orbital pairs
        bra_ids = []
        ket_ids = []
        for i in range(bra_dim[0], bra_dim[1]):
            for j in range(ket_dim[0], ket_dim[1]):
                bra_ids.append(i)
                ket_ids.append(j)

        # set number of Fock matrices
        self.num_matrices = len(bra_ids)

        return (bra_ids, ket_ids)

    def get_transformation_vectors(self, mints_type, mol_orbs, molecule):

        nocc, nvirt = self.get_num_orbitals(mol_orbs, molecule)

        # cases: oovv, ovov, ovvv, vvvv
        xmat = mol_orbs.alpha_orbitals(nocc, nvirt)
        ymat = mol_orbs.alpha_orbitals(nocc, nvirt)

        # case: oooo
        if mints_type == "OOOO":
            xmat = mol_orbs.alpha_orbitals(0, nocc)
            ymat = mol_orbs.alpha_orbitals(0, nocc)

        # case: ooov
        if mints_type == "OOOV":
            xmat = mol_orbs.alpha_orbitals(0, nocc)

        return (xmat, ymat)

    def set_fock_matrices_type(self, mints_type, fock_matrices):

        if mints_type != "OOVV":
            nfock = fock_matrices.number_of_fock_matrices()
            for i in range(nfock):
                fock_matrices.set_fock_type(fockmat.rgenj, i)

    def get_external_indexes(self, mints_type):

        if mints_type == "OOVV":
            return TwoIndexes(2, 3)

        return TwoIndexes(1, 3)

    def get_moints_type(self, mints_type):

        if mints_type == "OOOO":
            return moints.oooo

        if mints_type == "OOOV":
            return moints.ooov

        if mints_type == "OOVV":
            return moints.oovv

        if mints_type == "OVOV":
            return moints.ovov

        if mints_type == "OVVV":
            return moints.ovvv

        if mints_type == "VVVV":
            return moints.vvvv

        return None
