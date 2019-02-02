from .veloxchemlib import ElectronRepulsionIntegralsDriver
from .veloxchemlib import MolecularOrbitals
from .veloxchemlib import ScreeningContainer
from .veloxchemlib import mpi_master
from .veloxchemlib import ericut
from .veloxchemlib import molorb

from .aofockmatrix import AOFockMatrix
from .aodensitymatrix import AODensityMatrix

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
        self.num_mats = 0
    
    def compute_task(self, task, mol_orbs):

        self.compute(task.molecule, task.ao_basis, mol_orbs, task.mpi_comm,
                     task.ostream)
    
    def compute(self, molecule, ao_basis, mol_orbs, comm, ostream):
        
        # start timer
        start_time = tm.time()
        
        # mpi communicators
        self.set_global_comm(comm)
        
        # set up bra and ket orbitals list
        
        # dummy set up code here
        bra_ids = []
        ket_ids = []
        
        for i in range(50):
            for j in range(50):
                bra_ids.append(i)
                ket_ids.append(j)
    
        self.num_mats = len(bra_ids)
    
        if self.rank == mpi_master():
            self.print_header(ostream)
        
        # compute pair densities
        
        pair_den = mol_orbs.get_pair_density(bra_ids, ket_ids);
        
        pair_den.broadcast(self.rank, comm)
        
        # set up screening data
        
        eri_drv = ElectronRepulsionIntegralsDriver(self.rank,
                                                   self.nodes,
                                                   comm)
            
        qq_data = eri_drv.compute(self.get_qq_scheme(), self.eri_thresh,
                                  molecule, ao_basis)
    
        # compute multiple Fock matrices
        
        fock_mat = AOFockMatrix(pair_den)
        
        eri_drv.compute(fock_mat, pair_den, molecule, ao_basis, qq_data, comm)
        
        if self.rank == mpi_master():
            self.print_finish(start_time, ostream)

    def print_header(self, ostream):
        
        ostream.print_blank()
        ostream.print_header("Integrals Transformation (AO->MO) Driver Setup")
        ostream.print_header(46 * "=")
        ostream.print_blank()
        
        str_width = 80
        cur_str = "Number of Fock matrices      : " + str(self.num_mats)
        ostream.print_header(cur_str.ljust(str_width))
        cur_str = "ERI screening scheme         : " + self.get_qq_type()
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

    def get_qq_type(self):
        
        if self.qq_type == "QQ":
            return "Cauchy Schwarz"
        
        if self.qq_type == "QQR":
            return "Distance Dependent Cauchy Schwarz"
        
        if self.qq_type == "QQ_DEN":
            return "Cauchy Schwarz + Density"
        
        if self.qq_type == "QQR_DEN":
            return "Distance Dependent Cauchy Schwarz + Density"
        
        return "Undefined"

    def get_qq_scheme(self):
        
        if self.qq_type == "QQ":
            return ericut.qq
        
        if self.qq_type == "QQR":
            return ericut.qqr
        
        if self.qq_type == "QQ_DEN":
            return ericut.qqden

        if self.qq_type == "QQR_DEN":
            return ericut.qqrden

        return None;

