from .VeloxChemLib import ElectronRepulsionIntegralsDriver
from .VeloxChemLib import OverlapMatrix
from .VeloxChemLib import KineticEnergyMatrix
from .VeloxChemLib import NuclearPotentialMatrix
from .VeloxChemLib import OverlapIntegralsDriver
from .VeloxChemLib import KineticEnergyIntegralsDriver
from .VeloxChemLib import NuclearPotentialIntegralsDriver
from .VeloxChemLib import SADGuessDriver
from .VeloxChemLib import MolecularOrbitals
from .VeloxChemLib import mpi_master
from .VeloxChemLib import ericut
from .VeloxChemLib import molorb

from .aofockmatrix import AOFockMatrix
from .aodensitymatrix import AODensityMatrix

import numpy as np

from collections import deque

class ScfDriver:

    def __init__(self):

        self.guess_type = "SAD"
        
        # scf accelerator
        
        self.acc_type = "DIIS"
        self.max_err_vecs = 10
        self.max_iter = 50
        
        # screening scheme
        
        self.qq_type = "QQR"
        self.qq_dyn = True

        # thresholds
        
        self.conv_thresh = 1.0e-6
        self.eri_thresh = 1.0e-12
        self.ovl_thresh = 1.0e-12
    
        # iterations data
        
        self.iter_data = []
        self.is_converged = False
    
        # DIIS data lists
        
        self.fock_mat_list = deque()
        self.den_mat_list = deque()
    
    def compute(self, molecule, ao_basis, min_basis, comm, ostream):

        # MPI communicator data
        
        loc_rank = comm.Get_rank()
        loc_nodes = comm.Get_size()
        
        # print scf driver setup
        
        if loc_rank == mpi_master():
            self.print_header(ostream)
        
        # clear DIIS data
        
        self.fock_mat_list.clear()
        self.den_mat_list.clear()

        # compute one-electron integrals

        ovl_drv = OverlapIntegralsDriver.create(loc_rank, loc_nodes, comm)
        ovl_mat = ovl_drv.compute(molecule, ao_basis, ostream, comm)
    
        kin_drv = KineticEnergyIntegralsDriver.create(loc_rank, loc_nodes, comm)
        kin_mat = kin_drv.compute(molecule, ao_basis, ostream, comm)
        
        npot_drv = NuclearPotentialIntegralsDriver.create(loc_rank, loc_nodes,
                                                          comm)
        npot_mat = npot_drv.compute(molecule, ao_basis, ostream, comm)
        
        # DIIS method

        if self.acc_type == "DIIS":
            
            # initialize ERI driver
            
            eri_drv = ElectronRepulsionIntegralsDriver.create(loc_rank,
                                                              loc_nodes,
                                                              comm)
        
            qq_data = eri_drv.compute(self.get_qq_scheme(), self.eri_thresh,
                                      molecule, ao_basis, ostream, comm)
            
            # generate initial density
            
            den_mat = self.gen_guess_density(molecule, ao_basis, min_basis,
                                             ovl_mat, ovl_drv, loc_rank,
                                             loc_nodes, comm, ostream)
            
            # compute orthogonalization matrix
            
            if loc_rank == mpi_master():
                oao_mat = ovl_mat.get_ortho_matrix(self.ovl_thresh, ostream)
            
            # set up initial scf sata
        
            old_den_mat = AODensityMatrix(den_mat)
            
            fock_mat = AOFockMatrix(den_mat)
            
            mol_orbs = MolecularOrbitals()
            
            # scf iterations
            
            self.print_scf_title(ostream)
            
            for i in self.get_scf_range():
                
                # compute 2e part of AO fock matrix
                
                eri_drv.compute(fock_mat, den_mat, molecule, ao_basis, qq_data,
                                ostream, comm)
            
            
                # compute electronic energy components
            
                e_ee, e_kin, e_en = self.comp_energy(fock_mat, kin_mat,
                                                     npot_mat, den_mat)
                 
                # construct full AO Fock matrix
                
                self.comp_full_fock(fock_mat, kin_mat, npot_mat)

                # compute electronic gradient

                e_grad = self.comp_gradient(fock_mat, ovl_mat, den_mat)

                # compute density change from last iteration

                diff_den = self.comp_density_change(den_mat, old_den_mat)

                # process scf cycle data

                self.add_iter_data(i, e_ee, e_kin, e_en, e_grad, diff_den)
                self.check_convergence(i)
                self.print_iter_data(i, ostream)
                
                # store DIIS data
                
                self.store_diis_data(i, fock_mat, den_mat)
                
                #print(self.fock_mat_list)
                #print(self.den_mat_list)
                
                # compute molecular orbitals
                
                eff_fock_mat = self.get_effective_fock(fock_mat, ovl_mat,
                                                       oao_mat)
                
                self.apply_level_shift(eff_fock_mat, oao_mat, molecule)
                
                mol_orbs = self.gen_molecular_orbitals(eff_fock_mat, oao_mat,
                                                       ostream)
            
                # update density matrix
            
                old_den_mat = AODensityMatrix(den_mat)
                
                den_mat = self.gen_new_density(mol_orbs, molecule)
            
                if self.is_converged:
                    break
            
    def gen_guess_density(self, molecule, ao_basis, min_basis, ovl_mat,
                          ovl_driver, loc_rank, loc_nodes, comm, ostream):
        
        # guess: superposition of atomic densities
        
        if self.guess_type == "SAD":
            ovl_mat_sb = ovl_driver.compute(molecule, min_basis, ao_basis,
                                            ostream, comm)
                
            sad_drv = SADGuessDriver.create(loc_rank, loc_nodes, comm)
            den_mat = sad_drv.compute(molecule, min_basis, ao_basis, ovl_mat_sb,
                                      ovl_mat, ostream, comm)
            
            return den_mat
    
        # TODO: implemet restart or other types of density matrices
        
        return AODensityMatrix()

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
    
    def add_iter_data(self, i, e_ee, e_kin, e_en, e_grad, diff_den):
        # electronic energy

        e_elec = e_ee + e_kin + e_en
        
        # store scf iteration data
        
        if self.guess_type == "SAD":
            if i == 0:
                return
            elif i == 1:
                self.iter_data.append((0, e_elec, 0.0, e_grad, 0.0))
            else:
                de_elec = e_elec - self.iter_data[i-2][1]
                self.iter_data.append((i-1, e_elec, de_elec, e_grad, diff_den))
        else:
            if i == 0:
                self.iter_data.append((i, e_elec, 0.0, e_grad, 0.0))
            else:
                de_elec = e_elec - self.iter_data[i-1][1]
                self.iter_data.append((i, e_elec, de_elec, e_grad, diff_den))

    def check_convergence(self, i):
        # unpack scf cycle data
        
        if self.guess_type == "SAD":
            if i == 0:
                self.is_converged = False
                return
            j, e_elec, de_elec, e_grad, diff_den = self.iter_data[i-1]
        else:
            j, e_elec, de_elec, e_grad, diff_den = self.iter_data[i]

        # check scf convergence
        
        self.is_converged = False
        if (abs(de_elec)  < self.conv_thresh and
            abs(e_grad)   < self.conv_thresh and
            abs(diff_den) < self.conv_thresh):
            self.is_converged = True

    def get_scf_range(self):
        if self.guess_type == "SAD":
            return range(self.max_iter + 1)
        return rnage(self.max_iter)

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
        
        if self.guess_type == "SAD":
            if i == 0:
                return
            j, energy, denergy, egrad, dden = self.iter_data[i-1]
        else:
            j, energy, denergy, egrad, dden = self.iter_data[i]

        exec_str =  " " + (str(j)).rjust(3) + 4 * " "
        exec_str += ("{:7.12f}".format(energy)).center(27) + 3 * " "
        exec_str += ("{:5.10f}".format(denergy)).center(17) + 3 * " "
        exec_str += ("{:5.8f}".format(egrad)).center(15) + 3 * " "
        exec_str += ("{:5.8f}".format(dden)).center(15) + " "
        ostream.put_header(exec_str)

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
            return "Cauchy–Schwarz"
        if self.qq_type == "QQR":
            return "Distance Dependent Cauchy-Schwarz"
        if self.qq_type == "QQ_DEN":
            return "Cauchy–Schwarz + Density"
        if self.qq_type == "QQR_DEN":
            return "Distance Dependent Cauchy-Schwarz + Density"
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
