from .veloxchemlib import ElectronRepulsionIntegralsDriver
from .veloxchemlib import OverlapMatrix
from .veloxchemlib import KineticEnergyMatrix
from .veloxchemlib import NuclearPotentialMatrix
from .veloxchemlib import OverlapIntegralsDriver
from .veloxchemlib import KineticEnergyIntegralsDriver
from .veloxchemlib import NuclearPotentialIntegralsDriver
from .veloxchemlib import SADGuessDriver
from .veloxchemlib import MolecularOrbitals
from .veloxchemlib import ScreeningContainer
from .veloxchemlib import mpi_master
from .veloxchemlib import ericut
from .veloxchemlib import molorb

from .aofockmatrix import AOFockMatrix
from .aodensitymatrix import AODensityMatrix

from .denguess import DensityGuess

import numpy as np
import time as tm
import math

from collections import deque

class ScfDriver:

    def __init__(self):
        
        # density guess
        self.den_guess = DensityGuess("SAD")
        
        # scf accelerator
        self.acc_type = "L2_DIIS"
        self.max_err_vecs = 10
        self.max_iter = 50
        self.first_step = False
        
        # screening scheme
        self.qq_type = "QQ_DEN"
        self.qq_dyn = True
        
        # thresholds
        self.conv_thresh = 1.0e-6
        self.eri_thresh  = 1.0e-12
        self.ovl_thresh  = 1.0e-12
        self.diis_thresh = 0.2
        
        # iterations data
        self.iter_data = []
        self.is_converged = False
        self.skip_iter = False
        self.old_energy = 0.0
        self.num_iter = 0
        
        # DIIS data lists
        self.fock_matrices = deque()
        self.den_matrices = deque()
    
        # density matrix
        self.density = AODensityMatrix()

        # molecular orbitals
        self.mol_orbs = MolecularOrbitals()
    
        # nuclear-nuclear repulsion energy
        self.nuc_energy = 0.0
    
        # mpi information
        self.rank = 0
        self.nodes = 1
    

    def compute_task(self, task):

        self.compute(task.molecule, task.ao_basis, task.min_basis,
                     task.mpi_comm, task.ostream)
    
    def compute(self, molecule, ao_basis, min_basis, comm, ostream):
       
        self.rank = comm.Get_rank()
        self.nodes = comm.Get_size()
        
        self.nuc_energy = molecule.nuclear_repulsion_energy()
        
        if self.rank == mpi_master():
            self.print_header(ostream)
    
        # DIIS method
        if self.acc_type == "DIIS":
            self.comp_diis(molecule, ao_basis, min_basis, comm, ostream)
       
        # Two level DIIS method
        if self.acc_type == "L2_DIIS":
        
            # first step
            self.first_step = True
            
            old_thresh = self.conv_thresh
            self.conv_thresh = 1.0e-3
            
            old_max_iter = self.max_iter
            self.max_iter = 5
            
            val_basis = ao_basis.get_valence_basis()
            
            self.comp_diis(molecule, val_basis, min_basis, comm, ostream)
            
            # second step
            self.first_step = False

            self.diis_thresh = 1000.0
        
            self.conv_thresh = old_thresh
            
            self.max_iter = old_max_iter

            self.den_guess.guess_type = "PRCMO"

            self.comp_diis(molecule, ao_basis, val_basis, comm, ostream)

    def comp_diis(self, molecule, ao_basis, min_basis, comm, ostream):
        
        start_time = tm.time()

        self.fock_matrices.clear()
        
        self.den_matrices.clear()
    
        ovl_mat, kin_mat, npot_mat = self.comp_one_ints(molecule, ao_basis,
                                                        comm, ostream)
                                                        
        if self.rank == mpi_master():
            oao_mat = ovl_mat.get_ortho_matrix(self.ovl_thresh, ostream)
        else:
            oao_mat = None
    
        eri_drv = ElectronRepulsionIntegralsDriver.create(self.rank,
                                                          self.nodes,
                                                          comm)
                                                          
        qq_data = eri_drv.compute(self.get_qq_scheme(), self.eri_thresh,
                                  molecule, ao_basis, ostream)

        den_mat = self.comp_guess_density(molecule, ao_basis, min_basis,
                                          ovl_mat, comm, ostream)
                                          
        den_mat.broadcast(self.rank, comm)

        self.density = AODensityMatrix(den_mat)

        fock_mat = AOFockMatrix(den_mat)
    
        if self.rank == mpi_master():
            self.print_scf_title(ostream)

        for i in self.get_scf_range():
            
            eri_drv.compute(fock_mat, den_mat, molecule, ao_basis, qq_data,
                            ostream, comm)
            
            fock_mat.reduce_sum(self.rank, self.nodes, comm)
            
            e_ee, e_kin, e_en = self.comp_energy(fock_mat, kin_mat, npot_mat,
                                                 den_mat, comm)
                
            self.comp_full_fock(fock_mat, kin_mat, npot_mat)
                        
            e_grad = self.comp_gradient(fock_mat, ovl_mat, den_mat, comm)
        
            self.set_skip_iter_flag(i, e_grad)
                
            diff_den = self.comp_density_change(den_mat, self.density, comm)

            self.add_iter_data(e_ee, e_kin, e_en, e_grad, diff_den)

            self.check_convergence()
                
            self.print_iter_data(i, ostream)
                
            self.store_diis_data(i, fock_mat, den_mat)
    
            eff_fock_mat = self.get_effective_fock(fock_mat, ovl_mat, oao_mat)

            self.mol_orbs = self.gen_molecular_orbitals(eff_fock_mat, oao_mat,
                                                        ostream)
            
            self.density = AODensityMatrix(den_mat)

            den_mat = self.gen_new_density(molecule)
            
            den_mat.broadcast(self.rank, comm)
            
            if self.qq_dyn:
                qq_data.set_threshold(self.get_dyn_threshold(e_grad))
            
            if self.is_converged:
                break

        if self.rank == mpi_master():
            self.print_scf_finish(start_time, ostream)

    def comp_one_ints(self, molecule, basis, comm, ostream):
        
        ovl_drv = OverlapIntegralsDriver.create(self.rank, self.nodes, comm)
        ovl_mat = ovl_drv.compute(molecule, basis, ostream, comm)
        
        kin_drv = KineticEnergyIntegralsDriver.create(self.rank, self.nodes,
                                                      comm)
        kin_mat = kin_drv.compute(molecule, basis, ostream, comm)
        
        npot_drv = NuclearPotentialIntegralsDriver.create(self.rank, self.nodes,
                                                          comm)
        npot_mat = npot_drv.compute(molecule, basis, ostream, comm)
        
        return (ovl_mat, kin_mat, npot_mat)
    
    def comp_guess_density(self, molecule, ao_basis, min_basis, ovl_mat,
                           comm, ostream):

        # guess: superposition of atomic densities
        if self.den_guess.guess_type == "SAD":
           
            return self.den_guess.sad_density(molecule, ao_basis, min_basis,
                                              ovl_mat, self.rank, self.nodes,
                                              comm, ostream)
        
        # guess: projection of molecular orbitals from reduced basis
        if self.den_guess.guess_type == "PRCMO":
        
            if self.rank == mpi_master():
                return self.den_guess.prcmo_density(molecule, ao_basis,
                                                    min_basis, self.mol_orbs)
            else:
                return AODensityMatrix()

        return AODensityMatrix()

    def set_skip_iter_flag(self, i, e_grad):
        
        self.num_iter = i
        
        self.use_level_shift = False
        
        if i == 0:
            self.skip_iter = True
        else:
            if e_grad < self.diis_thresh:
                self.skip_iter = False
            else:
                self.skip_iter = True
                    
    def comp_energy(self, fock_mat, kin_mat, npot_mat, den_mat, comm):
        return (0.0, 0.0, 0.0)

    def comp_full_fock(self, fock_mat, kin_mat, npot_mat):
        return

    def comp_gradient(self, fock_mat, ovl_mat, den_mat, comm):
        return 0.0

    def comp_density_change(self, den_mat, old_den_mat, comm):
        return 0.0
    
    def store_diis_data(self, i, fock_mat, den_mat):
        return

    def get_effective_fock(self, fock_mat, ovl_mat, oao_mat):
        return None

    def gen_molecular_orbitals(self, fock_mat, oao_mat, ostream):
        return MolecularOrbitals()

    def gen_new_density(self, molecule):
        return AODensityMatrix()

    def get_dyn_threshold(self, e_grad):

        nteri = math.pow(10, math.floor(math.log10(e_grad)));

        nteri = 1.0e-8 * nteri
    
        if nteri > 1.0e-8:
            return 1.0e-8
        
        if nteri < self.eri_thresh:
            return self.eri_thresh

        return nteri
    
    def add_iter_data(self, e_ee, e_kin, e_en, e_grad, diff_den):
        
        e_elec = e_ee + e_kin + e_en + self.nuc_energy

        de_elec = e_elec - self.old_energy
        
        self.iter_data.append((e_elec, de_elec, e_grad, diff_den))
        
        self.old_energy = e_elec

    def check_convergence(self):
        
        self.is_converged = False
        
        if len(self.iter_data) > 0:
            
            e_elec, de_elec, e_grad, diff_den = self.iter_data[-1]
            
            if e_grad < self.conv_thresh:
                self.is_converged = True

    def get_scf_range(self):
        
        return range(self.max_iter + 1)

    def print_header(self, ostream):
        
        ostream.new_line()
        ostream.put_header("Self Consistent Field Driver Setup")
        ostream.put_header(36 * "=")
        ostream.new_line()
        
        str_width = 80
        cur_str = "Wave Function Model          : " + self.get_scf_type()
        ostream.put_header(cur_str.ljust(str_width))
        cur_str = "Initial Guess Model          : " + self.den_guess.guess_type
        ostream.put_header(cur_str.ljust(str_width))
        
        cur_str = "Convergence Accelerator      : " + self.get_acc_type()
        ostream.put_header(cur_str.ljust(str_width))
        cur_str = "Max. Number Of Iterations    : " + str(self.max_iter)
        ostream.put_header(cur_str.ljust(str_width))
        cur_str = "Max. Number Of Error Vectors : " + str(self.max_err_vecs)
        ostream.put_header(cur_str.ljust(str_width))
        cur_str = "Convergence Threshold        : " + \
            "{:.1e}".format(self.conv_thresh)
        ostream.put_header(cur_str.ljust(str_width))
        
        cur_str = "ERI screening scheme         : " + self.get_qq_type()
        ostream.put_header(cur_str.ljust(str_width))
        cur_str = "ERI screening mode           : " + self.get_qq_dyn()
        ostream.put_header(cur_str.ljust(str_width))
        cur_str = "ERI Screening Threshold      : " + \
            "{:.1e}".format(self.eri_thresh)
        ostream.put_header(cur_str.ljust(str_width))
        cur_str = "Linear Dependence Threshold  : " + \
            "{:.1e}".format(self.ovl_thresh)
        ostream.put_header(cur_str.ljust(str_width))
        ostream.new_line()

    def print_scf_title(self, ostream):
        
        if self.first_step:
            ostream.put_info("Starting Reduced Basis SCF calculation...")
        else:
            ostream.new_line()
            ostream.put_header("Iter. |   Hartree-Fock Energy, au   | "
                               "Energy Change, au |  Gradient Norm  | "
                               "Density Change |")
            ostream.put_header(92 * "-")

    def print_scf_finish(self, start_time, ostream):
    
        if self.first_step:
            valstr = "...done. SCF energy: "
            valstr += "{:.12f}".format(self.old_energy)
            valstr += " au. Time: "
            valstr += "{:.2f}".format(tm.time() - start_time) + " sec."
            ostream.put_info(valstr)
            ostream.new_line()
        else:
            valstr = "*** SCF "
            if self.is_converged:
                valstr += "converged in "
            else:
                valstr += "not converged in "
            valstr += str(self.num_iter)
            valstr += " iterations. Time: "
            valstr += "{:.2f}".format(tm.time() - start_time) + " sec."
            ostream.new_line()
            ostream.put_header(valstr.ljust(92))
            ostream.new_line()

    def print_iter_data(self, i, ostream):
        
        if self.rank == mpi_master():
            # no output for first step in two level DIIS
            if self.first_step:
                return
        
            # DIIS or second step in two level DIIS
            if i > 0:

                if len(self.iter_data) > 0:
                    te, diff_te, e_grad, diff_den = self.iter_data[-1]

                if i == 1:
                    diff_te = 0.0
                    diff_den = 0.0
            
                exec_str  = " " + (str(i)).rjust(3) + 4 * " "
                exec_str += ("{:7.12f}".format(te)).center(27) + 3 * " "
                exec_str += ("{:5.10f}".format(diff_te)).center(17) + 3 * " "
                exec_str += ("{:5.8f}".format(e_grad)).center(15) + 3 * " "
                exec_str += ("{:5.8f}".format(diff_den)).center(15) + " "

                ostream.put_header(exec_str)
                ostream.flush()

    def get_scf_energy(self):
        return self.old_energy

    def get_scf_type(self):
        return "Undefined"

    def get_guess_type(self):
        
        if self.den_guess.guess_type == "SAD":
            return "Superposition of Atomic Densities"
        
        return "Undefined"

    def get_acc_type(self):
        
        if self.acc_type == "DIIS":
            return "Direct Inversion of Iterative Subspace"
        
        if self.acc_type == "L2_DIIS":
            return "Two Level Direct Inversion of Iterative Subspace"
        
        return "Undefined"

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

    def get_qq_dyn(self):
        
        if self.qq_dyn:
            return "Dynamic"
        
        return "Static"

    def need_min_basis(self):
        
        if self.acc_type == "L2_DIIS":
            return True
        
        if self.den_guess.guess_type == "SAD":
            return True
        
        return False

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
