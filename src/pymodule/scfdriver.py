from .VeloxChemLib import ElectronRepulsionIntegralsDriver
from .VeloxChemLib import OverlapMatrix
from .VeloxChemLib import KineticEnergyMatrix
from .VeloxChemLib import NuclearPotentialMatrix
from .VeloxChemLib import OverlapIntegralsDriver
from .VeloxChemLib import KineticEnergyIntegralsDriver
from .VeloxChemLib import NuclearPotentialIntegralsDriver
from .VeloxChemLib import SADGuessDriver
from .VeloxChemLib import mpi_master
from .VeloxChemLib import ericut

from .aofockmatrix import AOFockMatrix
from .aodensitymatrix import AODensityMatrix

import numpy as np

class ScfDriver:

    def __init__(self):

        self.guess_type = "SAD"
        
        # scf accelerator
        
        self.acc_type = "DIIS"
        self.max_err_vecs = 10
        self.max_iter = 50
        
        # screeninf scheme
        
        self.qq_type = "QQR"
        self.qq_dyn = True

        # thresholds
        
        self.conv_thresh = 1.0e-6
        self.eri_thresh = 1.0e-12
        self.ovl_thresh = 1.0e-12
    
        # iterations data
        
        self.iter_data = []
    
    def compute(self, molecule, ao_basis, min_basis, comm, ostream):

        # MPI communicator data
        
        loc_rank = comm.Get_rank()
        loc_nodes = comm.Get_size()
        
        # print scf driver setup
        
        if loc_rank == mpi_master():
            self.print_header(ostream)

        # compute one-electron integrals

        ovldrv = OverlapIntegralsDriver.create(loc_rank, loc_nodes, comm)
        ovlmat = ovldrv.compute(molecule, ao_basis, ostream, comm)
    
        kindrv = KineticEnergyIntegralsDriver.create(loc_rank, loc_nodes, comm)
        kinmat = kindrv.compute(molecule, ao_basis, ostream, comm)
        
        npotdrv = NuclearPotentialIntegralsDriver.create(loc_rank, loc_nodes,
                                                         comm)
        npotmat = npotdrv.compute(molecule, ao_basis, ostream, comm)
        
        # DIIS method

        if self.acc_type == "DIIS":

            denmat = self.gen_guess_density(molecule, ao_basis, min_basis,
                                            ovlmat, ovldrv, loc_rank, loc_nodes,
                                            comm, ostream)
            
            if loc_rank == mpi_master():
                oaomat = ovlmat.get_ortho_matrix(self.ovl_thresh, ostream)
            
            eridrv = ElectronRepulsionIntegralsDriver.create(loc_rank,
                                                             loc_nodes,
                                                             comm)
        
            qqdata = eridrv.compute(self.get_qq_scheme(), self.eri_thresh,
                                    molecule, ao_basis, ostream, comm)
            
            # set up scf data
        
            refden = AODensityMatrix(denmat)
            
            fockmat = AOFockMatrix(denmat)
            
            # scf iterations
            
            self.print_scf_title(ostream)
            
            for i in range(self.max_iter):
                eridrv.compute(fockmat, denmat, molecule, ao_basis, qqdata,
                               ostream, comm)
            
                eelec, ekin, enpot = self.comp_energy(fockmat, kinmat, npotmat,
                                                      denmat)
            
                egrad = self.comp_gradient(fockmat, denmat)
                    
                dden = self.comp_density_change(denmat, refden)
                
                self.add_iter_data(i, eelec, ekin, enpot, egrad, dden)
                
                self.print_iter_data(i, ostream)
            
                if dden < self.conv_thresh:
                    break
            
    def gen_guess_density(self, molecule, ao_basis, min_basis, ovlmat,
                          ovl_driver, loc_rank, loc_nodes, comm, ostream):
        # Superposition of atomic densities
        
        if self.guess_type == "SAD":
            ovlmat_sb = ovl_driver.compute(molecule, min_basis, ao_basis,
                                           ostream, comm)
                
            saddrv = SADGuessDriver.create(loc_rank, loc_nodes, comm)
            denmat = saddrv.compute(molecule, min_basis, ao_basis, ovlmat_sb,
                                    ovlmat, ostream, comm)
            return denmat
    
        # TODO: implemet restart or other types of density matrices
        
        return AODensityMatrix()

    def comp_energy(self, fockmat, kinmat, npotmat, denmat):
        return (0.0, 0.0, 0.0)

    def comp_gradient(self, fockmat, denmat):
        return 0.0

    def comp_density_change(self, newden, oldden):
        return 0.0
    
    def add_iter_data(self, iterscf, eelec, ekin, enpot, egrad, dden):
        energy = eelec + ekin + enpot
        if iterscf == 0:
            self.iter_data.append((iterscf, energy, 0.0, egrad, dden))
        else:
            denergy = energy- self.iter_data[iterscf - 1][1]
            self.iter_data.append((iterscf, energy, denergy, egrad, dden))

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

    def print_iter_data(self, iterscf, ostream):
        i, energy, denergy, egrad, dden = self.iter_data[iterscf]
        exec_str =  " " + (str(i)).rjust(3) + 4 * " "
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
            return "Distance Dependent Cauchy Schwarz"
        if self.qq_type == "QQ_DEN":
            return "Cauchy–Schwarz and Density"
        if self.qq_type == "QQR_DEN":
            return "Distance Dependent Cauchy Schwarz and Density"
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
