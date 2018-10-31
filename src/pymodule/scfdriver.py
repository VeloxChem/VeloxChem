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

import numpy as np
import time as tm

from collections import deque

class ScfDriver:

    def __init__(self):
        
        self.guess_type = "SAD"
        
        # scf accelerator
        self.acc_type = "L2_DIIS"
        self.max_err_vecs = 10
        self.max_iter = 50
        
        # screening scheme
        self.qq_type = "QQ"
        self.qq_dyn = True
        
        # thresholds
        self.conv_thresh = 1.0e-6
        self.eri_thresh = 1.0e-12
        self.ovl_thresh = 1.0e-12
        self.diis_thresh = 0.1
        
        # iterations data
        self.iter_data = []
        self.is_converged = False
        self.skip_iter = False
        self.old_energy = 0.0
        
        # level shifting data
        self.use_level_shift = False
        self.level_shift = 0.2
        
        # DIIS data lists
        self.fock_matrices = deque()
        self.den_matrices = deque()

    def compute_task(self, task, comm):

        self.compute(task.molecule, task.ao_basis, task.min_basis, comm,
                     task.ostream)
    
    def compute(self, molecule, ao_basis, min_basis, comm, ostream):
       
        loc_rank = comm.Get_rank()
        loc_nodes = comm.Get_size()
        
        if loc_rank == mpi_master():
            self.print_header(ostream)
        
        self.fock_matrices.clear()
        self.den_matrices.clear()
        
        ovl_mat, kin_mat, npot_mat = self.comp_one_ints(molecule, ao_basis,
                                                        loc_rank, loc_nodes,
                                                        comm, ostream)
                                                        
        # DIIS or L2_DIIS method
        if self.acc_type == "DIIS" or self.acc_type == "L2_DIIS":
            eri_drv = ElectronRepulsionIntegralsDriver.create(loc_rank,
                                                              loc_nodes,
                                                              comm)
            qq_data = eri_drv.compute(self.get_qq_scheme(), self.eri_thresh,
                                      molecule, ao_basis, ostream, comm)
                                      
            if self.acc_type == "L2_DIIS":
                val_basis = ao_basis.get_valence_basis()
                val_orbs = self.comp_valence_scf(molecule, val_basis, min_basis,
                                                 loc_rank, loc_nodes, comm,
                                                 ostream)
                prj_orbs = val_orbs.insert(molecule, ao_basis, val_basis)
                den_mat = self.gen_new_density(prj_orbs, molecule)
            else:
                den_mat = self.gen_guess_density(molecule, ao_basis, min_basis,
                                                 ovl_mat, loc_rank, loc_nodes,
                                                 comm, ostream)
           
            if loc_rank == mpi_master():
                oao_mat = ovl_mat.get_ortho_matrix(self.ovl_thresh, ostream)
            
            old_den_mat = AODensityMatrix(den_mat)
            fock_mat = AOFockMatrix(den_mat)
            mol_orbs = MolecularOrbitals()
           
            self.print_scf_title(ostream)
            
            for i in self.get_scf_range():
                eri_drv.compute(fock_mat, den_mat, molecule, ao_basis, qq_data,
                                ostream, comm)
                e_ee, e_kin, e_en = self.comp_energy(fock_mat, kin_mat,
                                                     npot_mat, den_mat)
                self.comp_full_fock(fock_mat, kin_mat, npot_mat)
                
                e_grad = self.comp_gradient(fock_mat, ovl_mat, den_mat)
                diff_den = self.comp_density_change(den_mat, old_den_mat)
                
                self.set_skip_iter_flag(i, e_grad)
                self.add_iter_data(e_ee, e_kin, e_en, e_grad, diff_den)
                self.check_convergence()
                self.print_iter_data(i, ostream)
                
                self.store_diis_data(i, fock_mat, den_mat)
                
                eff_fock_mat = self.get_effective_fock(fock_mat, ovl_mat,
                                                       oao_mat)
                self.apply_level_shift(eff_fock_mat, oao_mat, molecule)
                mol_orbs = self.gen_molecular_orbitals(eff_fock_mat, oao_mat,
                                                       ostream)
            
                old_den_mat = AODensityMatrix(den_mat)
                den_mat = self.gen_new_density(mol_orbs, molecule)
                
                if self.is_converged:
                    break

    def comp_one_ints(self, molecule, basis, loc_rank, loc_nodes, comm,
                      ostream):
        
        ovl_drv = OverlapIntegralsDriver.create(loc_rank, loc_nodes, comm)
        ovl_mat = ovl_drv.compute(molecule, basis, ostream, comm)
        
        kin_drv = KineticEnergyIntegralsDriver.create(loc_rank, loc_nodes, comm)
        kin_mat = kin_drv.compute(molecule, basis, ostream, comm)
        
        npot_drv = NuclearPotentialIntegralsDriver.create(loc_rank, loc_nodes,
                                                          comm)
        npot_mat = npot_drv.compute(molecule, basis, ostream, comm)
        
        return (ovl_mat, kin_mat, npot_mat)
    
    def comp_valence_scf(self, molecule, val_basis, min_basis, loc_rank,
                         loc_nodes, comm, ostream):
        
        start_time = tm.time()
        
        ostream.put_info("Starting Valence SCF calculation...")
        ostream.new_line()
        
        ovl_mat, kin_mat, npot_mat = self.comp_one_ints(molecule, val_basis,
                                                        loc_rank, loc_nodes,
                                                        comm, ostream)
        den_mat = self.gen_guess_density(molecule, val_basis, min_basis,
                                         ovl_mat, loc_rank, loc_nodes, comm,
                                         ostream);
        
        eri_drv = ElectronRepulsionIntegralsDriver.create(loc_rank,
                                                          loc_nodes,
                                                          comm)
        qq_data = eri_drv.compute(self.get_qq_scheme(), self.eri_thresh,
                                  molecule, val_basis, ostream, comm)
       
        if loc_rank == mpi_master():
            oao_mat = ovl_mat.get_ortho_matrix(self.ovl_thresh, ostream)
        
        fock_mat = AOFockMatrix(den_mat)
        mol_orbs = MolecularOrbitals()
        
        for i in range(6):
            eri_drv.compute(fock_mat, den_mat, molecule, val_basis, qq_data,
                            ostream, comm)
            e_ee, e_kin, e_en = self.comp_energy(fock_mat, kin_mat, npot_mat,
                                                 den_mat)
            self.comp_full_fock(fock_mat, kin_mat, npot_mat)
            eff_fock_mat = self.get_effective_fock(fock_mat, ovl_mat, oao_mat)
            
            mol_orbs = self.gen_molecular_orbitals(eff_fock_mat, oao_mat,
                                                   ostream)
            den_mat = self.gen_new_density(mol_orbs, molecule)
        
        valstr = "Valence SCF energy: "
        valstr += "{:.12f}".format(e_ee + e_kin + e_en)
        valstr += " au. Time: "
        valstr += "{:.2f}".format(tm.time() - start_time) + " sec."
        ostream.put_info(valstr)
        ostream.new_line()
        
        return mol_orbs

    def gen_guess_density(self, molecule, ao_basis, min_basis, ovl_mat,
                          loc_rank, loc_nodes, comm, ostream):
        
        if self.guess_type == "SAD":
            ovl_driver = OverlapIntegralsDriver.create(loc_rank, loc_nodes,
                                                       comm)
            ovl_mat_sb = ovl_driver.compute(molecule, min_basis, ao_basis,
                                            ostream, comm)
            sad_drv = SADGuessDriver.create(loc_rank, loc_nodes, comm)
            den_mat = sad_drv.compute(molecule, min_basis, ao_basis, ovl_mat_sb,
                                      ovl_mat, ostream, comm)
            return den_mat
        
        # TODO: implemet restart or other types of density matrices

        return AODensityMatrix()

    def set_skip_iter_flag(self, i, e_grad):
        
        if self.acc_type == "L2_DIIS":
            if i == 0:
                self.skip_iter = True
            else:
                self.skip_iter = False
            self.use_level_shift = False
            return
        
        if self.guess_type == "SAD" and i < 2:
            self.skip_iter = True
            self.use_level_shift = True
            return
        
        if e_grad < self.diis_thresh:
            self.skip_iter = False
            self.use_level_shift = False
        else:
            self.skip_iter = True
            self.use_level_shift = True

    def comp_energy(self, fock_mat, kin_mat, npot_mat, den_mat):
        return (0.0, 0.0, 0.0)

    def comp_full_fock(self, fock_mat, kin_mat, npot_mat):
        return

    def comp_gradient(self, fock_mat, ovl_mat, den_mat):
        return 0.0

    def comp_density_change(self, den_mat, old_den_mat):
        return 0.0
    
    def store_diis_data(self, i, fock_mat, den_mat):
        return

    def get_effective_fock(self, fock_mat, ovl_mat, oao_mat):
        return None

    def apply_level_shift(self, fock_mat, oao_mat, molecule):
        return

    def gen_molecular_orbitals(self, fock_mat, oao_mat, ostream):
        return MolecularOrbitals()

    def gen_new_density(self, mol_orbs, molecule):
        return AODensityMatrix()
    
    def add_iter_data(self, e_ee, e_kin, e_en, e_grad, diff_den):
        e_elec = e_ee + e_kin + e_en
        de_elec = e_elec - self.old_energy
        
        self.iter_data.append((e_elec, de_elec, e_grad, diff_den))
        self.old_energy = e_elec

    def check_convergence(self):
        
        self.is_converged = False
        ddim = len(self.iter_data)
        
        if ddim > 0:
            e_elec, de_elec, e_grad, diff_den = self.iter_data[-1]
            if (abs(de_elec)  < self.conv_thresh and
                abs(e_grad)   < self.conv_thresh and
                abs(diff_den) < self.conv_thresh):
                self.is_converged = True

    def get_scf_range(self):
        
        if self.guess_type == "SAD":
            return range(self.max_iter + 1)
        
        return range(self.max_iter)

    def print_header(self, ostream):
        
        ostream.new_line()
        ostream.put_header("Self Consistent Field Driver Setup")
        ostream.put_header(36 * "=")
        ostream.new_line()
        
        str_width = 80
        cur_str = "WaveFunction Model           : " + self.get_scf_type()
        ostream.put_header(cur_str.ljust(str_width))
        cur_str = "Initial Guess Model          : " + self.get_guess_type()
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
        
        ostream.new_line()
        ostream.put_header("Iter. |   Hartree-Fock Energy, au   | "
                           "Energy Change, au |  Gradient Norm  | "
                           "Density Change |")
        ostream.put_header(92 * "-")

    def print_iter_data(self, i, ostream):
        
        if self.guess_type == "SAD" and i == 0:
            return
        
        ddim = len(self.iter_data)
        if ddim > 0:
            energy, denergy, egrad, dden = self.iter_data[-1]
            if self.guess_type == "SAD" and i == 1:
                denergy = 0.0
                dden = 0.0
            
            exec_str =  " " + (str(i)).rjust(3) + 4 * " "
            exec_str += ("{:7.12f}".format(energy)).center(27) + 3 * " "
            exec_str += ("{:5.10f}".format(denergy)).center(17) + 3 * " "
            exec_str += ("{:5.8f}".format(egrad)).center(15) + 3 * " "
            exec_str += ("{:5.8f}".format(dden)).center(15) + " "
            ostream.put_header(exec_str)
            ostream.flush()

    def get_scf_energy(self, enuc):
        return self.old_energy + enuc

    def get_scf_type(self):
        return "Undefined"

    def get_guess_type(self):
        
        if self.guess_type == "SAD":
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
        
        if self.guess_type == "SAD":
            return True
        
        return False

    def get_qq_scheme(self):
        
        if self.qq_type == "QQ":
            return ericut.qq
        
        if self.qq_type == "QQR":
            return ericut.qqr
        
        return None;
