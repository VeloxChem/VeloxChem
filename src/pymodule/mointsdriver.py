from .veloxchemlib import ElectronRepulsionIntegralsDriver
from .veloxchemlib import MolecularOrbitals
from .veloxchemlib import ScreeningContainer
from .veloxchemlib import mpi_master
from .veloxchemlib import ericut
from .veloxchemlib import fockmat
from .veloxchemlib import molorb
from .veloxchemlib import moints
from .veloxchemlib import MOIntsBatch
from .veloxchemlib import DenseMatrix
from .veloxchemlib import TwoIndexes

from .aofockmatrix import AOFockMatrix
from .aodensitymatrix import AODensityMatrix

from .qqscheme import get_qq_type
from .qqscheme import get_qq_scheme

import numpy as np
import time as tm
import math

class MOIntegralsDriver:

    def __init__(self):
        
        # screening scheme
        self.qq_type = "QQ_DEN"
        self.eri_thresh = 1.0e-12
    
        # mpi information
        self.rank = 0
        self.nodes = 1
        self.dist_exec = False
    
        # Fock matrices
        self.num_matrices = 0
        self.batch_size = 3000
    
    def compute_task(self, task, mol_orbs, mints_type):

        return self.compute(task.molecule, task.ao_basis, mol_orbs, mints_type,
                            task.mpi_comm, task.ostream)
    
    def compute(self, molecule, ao_basis, mol_orbs, mints_type, comm, ostream):
        
        # start timer
        start_time = tm.time()
        
        # mpi communicators
        self.set_global_comm(comm)
        
        # set up bra and ket orbitals list
        
        bra_ids, ket_ids = self.set_orbital_pairs(mints_type, mol_orbs,
                                                  molecule)
    
        if self.rank == mpi_master():
            self.print_header(ostream)
        
        # create MO integrals batch
        batch_pos, batch_dim = self.get_batch_dimensions(bra_ids)
        
        # set up screening data
        eri_drv = ElectronRepulsionIntegralsDriver(self.rank, self.nodes, comm)
            
        qq_data = eri_drv.compute(get_qq_scheme(self.qq_type), self.eri_thresh,
                                  molecule, ao_basis)
                               
        # initialize MO integrals batch
        moints_batch = MOIntsBatch()
        
        # set up properties of MO integrals batch
        moints_batch.set_ext_indexes(self.get_external_indexes(mints_type))
        moints_batch.set_batch_type(self.get_moints_type(mints_type))
        
        xmat, ymat = self.get_transformation_vectors(mints_type, mol_orbs,
                                                     molecule)
        
        # loop over batches
        for i in range(len(batch_dim)):
            
            spos = batch_pos[i]
            epos = batch_pos[i] + batch_dim[i]
            
            cur_bra_ids = bra_ids[spos:epos]
            cur_ket_ids = ket_ids[spos:epos]
        
            # compute pair densities
            pair_den = mol_orbs.get_pair_density(cur_bra_ids, cur_ket_ids);
            pair_den.broadcast(self.rank, comm)
        
            # compute multiple Fock matrices
            fock_mat = AOFockMatrix(pair_den)
            self.set_fock_matrices_type(mints_type, fock_mat)
            eri_drv.compute(fock_mat, pair_den, molecule, ao_basis,
                            qq_data, comm)
        
            # transform Fock matrices and add MO integrals to batch
            moints_batch.append(fock_mat, xmat, ymat, cur_bra_ids, cur_ket_ids)
        
        if self.rank == mpi_master():
            self.print_finish(start_time, ostream)

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

    def print_finish(self, start_time, ostream):
        
            valstr = "*** Integrals transformation (AO->MO) done in "
            valstr += "{:.2f}".format(tm.time() - start_time) + " sec."
            ostream.print_blank()
            ostream.print_header(valstr.ljust(92))
            ostream.print_blank()
    
    def set_global_comm(self, comm):
        
        self.rank = comm.Get_rank()
        self.nodes = comm.Get_size()
    
        if (self.nodes > 1):
            self.dist_exec = True

    def get_batch_dimensions(self, vec_ids):
        
        vdim = len(vec_ids)
        
        # determine number of batches
        nbatch = vdim // self.batch_size
        if (vdim % self.batch_size) != 0:
            nbatch = nbatch + 1
        
        # compute dimensions of batches
        bdim = vdim // nbatch
        brem = vdim %  nbatch
        
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
        nvirt = mol_orbs.get_number_mos() - nocc

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




