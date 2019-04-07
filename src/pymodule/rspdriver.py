from .veloxchemlib import mpi_master

from .outputstream import OutputStream
from .tdaexcidriver import TDAExciDriver
from .lrsolver import LinearResponseSolver
from .qqscheme import get_qq_type

import sys


class ResponseDriver:
    """Implements response driver.
        
    Implements response driver for molecular property calculations using
    conventional Hatree-Fock/Kohn-Sham response theory.
    
    rank
        The rank of MPI process.
    nodes
        The number of MPI processes.
    """

    def __init__(self, rsp_input=None):
        """Initializes Response driver.
            
        Initializes Response driver to default setup.
        """
        
        # default calculation type
        self.prop_type = "SINGEX_TDA"
        self.nstates = 3

        # default solver settings
        self.conv_thresh = 1.0e-4
        self.max_iter = 50
        
        # default ERI settings
        self.eri_thresh  = 1.0e-15
        self.qq_type = "QQ_DEN"
        
        if rsp_input:

            if 'conv_thresh' in rsp_input:
                self.conv_thresh = rsp_input['conv_thresh']
            if 'max_iter' in rsp_input:
                self.max_iter = rsp_input['max_iter']
            if 'eri_thresh' in rsp_input:
                self.eri_thresh = rsp_input['eri_thresh']
            if 'qq_type' in rsp_input:
                self.qq_type = rsp_input['qq_type']

            if rsp_input['property'].lower() == 'polarizability':
                self.prop_type = 'POLARIZABILITY'
                self.a_ops, self.b_ops = rsp_input['operators']
                self.frequencies = rsp_input['frequencies']

            elif rsp_input['property'].lower() == 'absorption':
                self.prop_type = 'SINGEX_TDA'
                self.nstates = rsp_input['nstates']
        
        # mpi information
        self.rank = 0
        self.nodes = 1
    
    def compute_task(self, mol_orbs, task):
        """Performs molecular property calculation using response theory.
            
        Performs molecular property calculation using data from MPI task.
            
        Parameters
        ----------
        mol_orbs
            The molecular orbitals.
        task
            The response input data as MPI task.
        """

        self.compute(mol_orbs, task.molecule, task.ao_basis, task.mpi_comm,
                     task.ostream)
    
    def compute(self, mol_orbs, molecule, ao_basis, comm,
                ostream=OutputStream(sys.stdout)):
        """Performs molecular property calculation.
            
        Performs molecular property calculation using molecular data, MPI
        communicator and output stream.
        
        Parameters
        ----------
        mol_orbs
            The molecular orbitals.
        molecule
            The molecule.
        ao_basis
            The AO basis set.
        min_basis
            The minimal AO basis set.
        comm
            The MPI communicator.
        ostream
            The output stream.
        """
        
        self.rank = comm.Get_rank()
        self.nodes = comm.Get_size()
        
        if self.rank == mpi_master():
            self.print_header(ostream)

        # TDA singlet/triplet excited states
        
        if self.prop_type.upper() in ["SINGEX_TDA", "TRIPEX_TDA"]:
            tda_exci = TDAExciDriver(self.rank, self.nodes)
            
            tda_exci.set_number_states(self.nstates)
            tda_exci.set_eri(self.eri_thresh, self.qq_type)
            tda_exci.set_solver(self.conv_thresh, self.max_iter)
            
            return tda_exci.compute(mol_orbs, molecule, ao_basis, comm,
                                    ostream)

        # Linear response solver

        if self.prop_type.upper() in ["POLARIZABILITY"]:
            lr_solver = LinearResponseSolver(self.a_ops, self.b_ops,
                                             self.frequencies)

            lr_solver.set_eri(self.eri_thresh, self.qq_type)
            lr_solver.set_solver(self.conv_thresh, self.max_iter)

            return lr_solver.compute(mol_orbs, molecule, ao_basis, comm,
                                     ostream)

    def print_header(self, ostream):
        """Prints response driver setup header to output stream.
            
        Prints molecular property calculation setup details to output stream.
        
        Parameters
        ----------
        ostream
            The output stream.
        """
    
        ostream.print_blank()
        ostream.print_header("Response Driver Setup")
        ostream.print_header(23 * "=")
        ostream.print_blank()
        
        str_width = 60
       
        cur_str = "Molecular Property Type   : " + self.prop_str()
        ostream.print_header(cur_str.ljust(str_width))
        
        if self.prop_type in ['SINGEX_TDA', 'TRIPEX_TDA']:
            cur_str = "Number Of Excited States  : " + str(self.nstates)
            ostream.print_header(cur_str.ljust(str_width))
        
        if self.prop_type in ['SINGEX_TDA', 'TRIPEX_TDA']:
            cur_str = "Response Equations Type   : Tamm-Dancoff"
            ostream.print_header(cur_str.ljust(str_width))
        
        cur_str = "Max. Number Of Iterations : " + str(self.max_iter)
        ostream.print_header(cur_str.ljust(str_width))
        cur_str = "Convergence Threshold     : " + \
            "{:.1e}".format(self.conv_thresh)
        ostream.print_header(cur_str.ljust(str_width))

        cur_str = "ERI screening scheme      : " + get_qq_type(self.qq_type)
        ostream.print_header(cur_str.ljust(str_width))
        cur_str = "ERI Screening Threshold   : " + \
            "{:.1e}".format(self.eri_thresh)
        ostream.print_header(cur_str.ljust(str_width))
        ostream.print_blank()

        ostream.flush()

    def prop_str(self):
        """Gets string with type of molecular property calculation.
            
        Gets string with type of molecular property calculation (Excited states,
        linear and non-linear spectroscopies).
            
        Returns
        -------
        The string with type of molecular property calculation.
        """
        
        if self.prop_type == "SINGEX_TDA":
            return "Singlet Excited States"

        if self.prop_type == "POLARIZABILITY":
            return "Polarizability"

        return "Undefined"
