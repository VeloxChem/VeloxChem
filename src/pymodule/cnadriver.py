#
#                           VELOXCHEM 1.0-RC2
#         ----------------------------------------------------
#                     An Electronic Structure Code
#
#  Copyright © 2018-2021 by VeloxChem developers. All rights reserved.
#  Contact: https://veloxchem.org/contact
#
#  SPDX-License-Identifier: LGPL-3.0-or-later
#
#  This file is part of VeloxChem.
#
#  VeloxChem is free software: you can redistribute it and/or modify it under
#  the terms of the GNU Lesser General Public License as published by the Free
#  Software Foundation, either version 3 of the License, or (at your option)
#  any later version.
#
#  VeloxChem is distributed in the hope that it will be useful, but WITHOUT
#  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
#  License for more details.
#
#  You should have received a copy of the GNU Lesser General Public License
#  along with VeloxChem. If not, see <https://www.gnu.org/licenses/>.

from mpi4py import MPI
import numpy as np
import time as tm
import sys

from .veloxchemlib import bohr_in_angstroms
from .veloxchemlib import mpi_master
from .veloxchemlib import CommonNeighbors
from .outputstream import OutputStream
from .errorhandler import assert_msg_critical
from .inputparser import parse_input
from .molecule import Molecule

class CnaAnalysisDriver:
    """
    Implements CNA analysis driver.

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.

    Instance variable
        - comm: The MPI communicator.
        - rank: The MPI rank.
        - nodes: Number of MPI processes.
        - ostream: The output stream.
        - cna_bond: The cut-off radius for chemical bond in CNA analysis.
        - cna_rcut: The cut-off radius for chemical bonds environment in
          CNA analysis.
        - xyz_files: The list of xyz files for correlation analysis.
        - lhs_files: The lhs list of xyz files for correlation analysis.
        - rhs_files: The rhs list of xyz files for correlation analysis.
    """

    def __init__(self, comm=None, ostream=None):
        """
        Initializes CNA analysis driver.
        """

        if comm is None:
            comm = MPI.COMM_WORLD

        if ostream is None:
            ostream = OutputStream(sys.stdout)

        # mpi information
        self.comm = comm
        self.rank = self.comm.Get_rank()
        self.nodes = self.comm.Get_size()

        # output stream
        self.ostream = ostream
        
        # CNA parameters
        self.cna_bond = None
        self.cna_rcut = None

        # list of xyz files
        self.xyz_files = None
        self.lhs_files = None
        self.rhs_files = None

    def update_settings(self, cna_dict):
        """
        Updates settings in CNA analysis driver.

        :param cna_dict:
            The dictionary of CNA analysis settings.
        """

        cna_keywords = {
            'xyz_files': 'list',
            'lhs_files': 'list',
            'rhs_files': 'list',
            'cna_bond': 'float',
            'cna_rcut': 'float',
        }
        
        parse_input(self, cna_keywords, cna_dict)
        
        # update CNA bond cut-off radius
        if self.cna_bond is None:
            self.cna_bond = 3.0
        else:
            self.cna_bond /= bohr_in_angstroms()
            
        # update CNA bond environment cut-off radius
        if self.cna_rcut is None:
            self.cna_rcut = 4.5
        else:
            self.cna_rcut /= bohr_in_angstroms()
            
    def compute(self):
        """
        Performs CNA correlation analysis.
        """
        
        self.print_header()
        
        if self.xyz_files is not None:
            self.compute_cna_self()
            
        if self.lhs_files is not None and self.rhs_files is not None:
            self.compute_cna_cross()
            
    def compute_cna_self(self):
        """
        Computes CNA self correlation analysis.
        """
        
        molecules = [Molecule.read_xyz(finp) for finp in self.xyz_files]
        
        cnas = []
        for mol in molecules:
            mcna = CommonNeighbors(mol, self.cna_bond)
            mcna.generate(self.cna_rcut)
            cnas.append(mcna)
        
        self.ostream.print_blank()
        self.ostream.print_header(29 * '-')
        self.ostream.print_header('! Self Correlation Function !')
        self.ostream.print_header(29 * '-')
        self.ostream.print_header('! (i) ! (j) !      J_ij     !')
        self.ostream.print_header(29 * '-')
        matidx = np.tril_indices(len(cnas))
        for i,j in zip(matidx[0],matidx[1]):
            jval = cnas[i].comp_cna(cnas[j])
            self.ostream.print_header('! {:^3d} ! {:^3d} !      {:^.2f}     !'.format(
                i, j, jval))
        self.ostream.print_header(29 * '-')
        
    def compute_cna_cross(self):
        """
        Computes CNA cross correlation analysis.
        """
    
        lhs_molecules = [Molecule.read_xyz(finp) for finp in self.lhs_files]
        rhs_molecules = [Molecule.read_xyz(finp) for finp in self.rhs_files]
    
        lhs_cnas = []
        for mol in lhs_molecules:
            mcna = CommonNeighbors(mol, self.cna_bond)
            mcna.generate(self.cna_rcut)
            lhs_cnas.append(mcna)
            
        rhs_cnas = []
        for mol in rhs_molecules:
            mcna = CommonNeighbors(mol, self.cna_bond)
            mcna.generate(self.cna_rcut)
            rhs_cnas.append(mcna)
    
        self.ostream.print_blank()
        self.ostream.print_header(30 * '-')
        self.ostream.print_header('! Cross Correlation Function !')
        self.ostream.print_header(30 * '-')
        self.ostream.print_header('! (i) ! (j) !      J_ij      !')
        self.ostream.print_header(30 * '-')
        for i in range(len(lhs_cnas)):
            for j in range(len(rhs_cnas)):
                jval = lhs_cnas[i].comp_cna(rhs_cnas[j])
                self.ostream.print_header('! {:^3d} ! {:^3d} !      {:^.2f}      !'.format(
                    i, j, jval))
        self.ostream.print_header(30 * '-')
        
    def print_header(self):
        """
        Prints header for the CNA driver.
        """
        
        self.ostream.print_blank()
        self.ostream.print_header('CNA Driver Setup')
        self.ostream.print_header(18 * '=')
        self.ostream.print_blank()
        
        str_width = 40
        cur_str = 'Bond Distance Cut-Off Radius     : {:.2f}'.format(
            self.cna_bond)
        self.ostream.print_header(cur_str.ljust(str_width))
        cur_str = 'Bonds Environment Cut-Off Radius : {:.2f}'.format(
            self.cna_rcut)
        self.ostream.print_header(cur_str.ljust(str_width))
       
        task_str = 'None'
        if self.xyz_files is not None:
            task_str = 'Self'
        if self.lhs_files is not None and self.rhs_files is not None:
            task_str = 'Cross'
        cur_str = 'Correlation type                 : ' + task_str
        self.ostream.print_header(cur_str.ljust(str_width))
       
        self.ostream.print_blank()
        self.ostream.flush()
