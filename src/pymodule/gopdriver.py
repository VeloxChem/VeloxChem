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
import random as rn
import math
import sys

from .veloxchemlib import bohr_in_angstroms
from .veloxchemlib import mpi_master
from .veloxchemlib import hartree_in_kcalpermol
from .outputstream import OutputStream
from .errorhandler import assert_msg_critical
from .inputparser import parse_input
from .veloxchemlib import CommonNeighbors
from .molecule import Molecule
from .treeblock import TreeBlock
from .xtbdriver import XTBDriver
from .xtbgradientdriver import XTBGradientDriver
from .optimizationdriver import OptimizationDriver

class GlobalOptimizationDriver:
    """
    Implements global optimization driver.

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
        - max_conformers: The maximum number of conformers generated by
          tree growth algorithm.
        - conformer_energies: The list of conformer energies.
        - conformer_geometries: The list of conformer geometries.
        - energy_thresh: The energetic acceptance threshold for conformer.
        - similarity_thresh: The similarity acceptance threshold for conformer.
        - growth_steps: The number of steps in tree growth algorithm.
        - popt_steps: The number of partial geometry optimzation steps.
        - block_data: The geometrical data of tree build block.
    """

    def __init__(self, comm=None, ostream=None):
        """
        Initializes global optimization driver.
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
        
        # tree growth parameters
        self.max_conformers = 500
        self.conformer_energies = []
        self.conformer_geometries = []
        self.energy_thresh = 0.1594
        self.similarity_thresh = 0.1594
        self.growth_steps = 0
        self.popt_steps = 0
        self.block_data = None
        self.tree_growth_unit = None
        

    def update_settings(self, gop_dict):
        """
        Updates settings in global optimization driver.

        :param gop_dict:
            The dictionary of global optimization settings.
        """

        gop_keywords = {
            'max_conformers' : 'int',
            'energy_thresh' : 'float',
            'similarity_thresh' : 'float',
            'growth_steps' : 'int',
            'popt_steps' : 'int',
            'block_data' : 'list',
            'cna_bond': 'float',
            'cna_rcut': 'float',
        }
        
        parse_input(self, gop_keywords, gop_dict)
        
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
            
        # update tree growth unit
        self.tree_growth_unit = TreeBlock(self.block_data)
            
    def compute(self, filename, seed):
        """
        Performs global optimization with given molecular seed.
        
        :param filename:
            The name of input file.
        :param seed:
            The molecular seed.
        """
        
        self.print_header()
    
        # generate conformer candidates
        molgeoms = []
        for _ in range(self.max_conformers):
            mol = self.generate_conformer(seed)
            if mol is not None:
                molgeoms.append(mol)
        
        # select unique conformer candidates
        confgeoms = self.remove_duplicates(molgeoms)
        
        # optimize conformers
        opt_dict = {'max_iter': '900', 'conv_energy': 1.0e-5, 
                    'conv_grms': 1.0e-3, 'conv_gmax': 3.0e-3, 
                    'conv_drms': 1.0e-3, 'conv_dmax': 3.0e-3}
        for i, confgeom in enumerate(confgeoms):
        
            self.ostream.print_blank()
            cur_str = 'Full Geometry Optimization: Conformer Candidate {:d}'.format(i)
            self.ostream.print_header(cur_str)
            self.ostream.print_header(len(cur_str) * '=')
            self.ostream.print_blank()
            
            xtb_drv = XTBDriver(self.comm)
            xtb_drv.set_method('gfn2')
            grad_drv = XTBGradientDriver(xtb_drv, self.comm, self.ostream)
            fname = 'conformer.{:d}'.format(i)
            opt_drv = OptimizationDriver(fname, grad_drv, 'XTB')
            opt_drv.update_settings(opt_dict)
            mol, ene = opt_drv.compute(confgeom, None, None)
            self.conformer_energies.append(ene) 
            self.conformer_geometries.append(mol)

        # sort conformers according to energies 
        zpairs = zip(self.conformer_energies, self.conformer_geometries) 
        spairs = sorted(zpairs, key=lambda pair: pair[0])  
        molenes = [x for x,_ in spairs]
        molgeoms = [x for _,x in spairs]
        self.conformer_energies = molenes 
        self.conformer_geometries = molgeoms 
        
        # print energies and geometries 
        self.print_conformers(filename) 
        

    def generate_conformer(self, seed):
        """
        Generates conformer for given molecular seed.
        
        :param seed:
            The molecular seed.
        :return:
            The generated molecular conformer.
        """
        
        if self.growth_steps == 0:
            return None
        else:
            cur_str = 'Generating conformer candidate...'
            self.ostream.print_info(cur_str)
            
            mol = Molecule(seed)
            
            if self.popt_steps > 0:
            
                # add growth units with optimization 
                nsteps = self.growth_steps // self.popt_steps
                if nsteps > 0:
                    for i in range(self.popt_steps): 
                        natoms = mol.number_of_atoms()
                        for _ in range(nsteps):
                            if not self.add_build_block(mol):
                                return None
                        fstr = 'xyz 1-{:d}'.format(natoms)
                        opt_dict = {'max_iter': '900', 'constraints': ['$freeze', fstr]}
                        xtb_drv = XTBDriver(self.comm)
                        xtb_drv.set_method('gfn2')
                        grad_drv = XTBGradientDriver(xtb_drv, self.comm, self.ostream)
                        filename = 'partial.optimization.{:d}'.format(i)
                        opt_drv = OptimizationDriver(filename, grad_drv, 'XTB')
                        opt_drv.update_settings(opt_dict)
                        mol = opt_drv.compute(mol, None, None)
                
                # add remaining growth units
                nsteps = self.growth_steps % self.popt_steps
                if nsteps > 0:
                    for _ in range(nsteps):
                        if not self.add_build_block(mol):
                            return None

            else:
            
                # add growth units
                for _ in range(self.growth_steps):
                    if not self.add_build_block(mol):
                        return None
                
            cur_str = '...done.'
            self.ostream.print_info(cur_str)
            return mol
            
    def add_build_block(self, molecule):
        """
        Adds the building block to given molecule.
        
        :param molecule:
            The initial molecule.
        :return:
            True if build block is added to molecule, False otherwise.
        """
        
        for i, atom in enumerate(self.tree_growth_unit.atoms):
            atmxyz = self.select_new_coordinates(i, molecule)
            if atmxyz is None:
                return False
            atmlbl = self.tree_growth_unit.atoms[i]
            molecule.add_atom(atmlbl, atmxyz[0], atmxyz[1], atmxyz[2])
        return True
    
    def select_ref_atom(self, index, molecule):
        """
        Select reference atom to append another atom in molecule.
        
        :param index:
            The index of atom in tree block.
        :param molecule:
            The initial molecule.
        :return:
            The randomly selected reference atom.
        """
        ref_indexes = []
        mfval = self.tree_growth_unit.max_cnums[index]
        ref_atom = self.tree_growth_unit.ref_atoms[index]
        for i in molecule.atom_indexes(ref_atom):
            cfval = 1 + molecule.coordination_number(i, self.cna_bond)
            if cfval <= mfval:
                ref_indexes.append(i)
        if len(ref_indexes) > 0:
            return rn.choice(ref_indexes)
        else:
            return None
        
    def select_new_coordinates(self, index, molecule):
        """
        Select new coordinates of atom to be added to molecule.
        
        :param index:
            The index of atom in tree block.
        :param molecule:
            The initial molecule.
        :return:
            The randomly selected coordinates of new atom.
        """
        
        while True:
        
            # select reference atom
            refids = self.select_ref_atom(index, molecule)
            if refids is None:
                return None
            refxyz = molecule.get_atom_coordinates(refids)
            radius = self.tree_growth_unit.bond_lengths[index]
            radius /= bohr_in_angstroms()
            
            minrad = 0.9 * radius
        
            # generate spherical mesh
            phi = np.linspace(0, np.pi, 10)
            theta = np.linspace(0, 2 * np.pi, 10)
            nx = refxyz[0] + radius * np.outer(np.sin(theta), np.cos(phi))
            ny = refxyz[1] + radius * np.outer(np.sin(theta), np.sin(phi))
            nz = refxyz[2] + radius * np.outer(np.cos(theta), np.ones_like(phi))
        
            # screen viable coordinates of new atom
            coords = []
            for x, y, z in zip(nx.flat, ny.flat, nz.flat):
                if molecule.min_distance(x, y, z) > minrad:
                    coords.append([x, y, z])
                    
            # select coordinates
            if len(coords) > 0:
                return rn.choice(coords)
                
        return None
        
    def remove_duplicates(self, molecules):
        """
        Removes duplicates from molecules list according to similarity index.
           
        :param molecules:
            The list of molecules.
        :return:
            The list of unique molecules.
        """
           
        if len(molecules) > 0:
            molgeoms = []
            for mol in molecules:
                if len(molgeoms) == 0:
                    molgeoms.append(mol)
                else:
                    self.add_unique(molgeoms, mol, 0.9)
            return molgeoms
        else:
            return None
                
    def add_unique(self, molecules, molecule, jcutoff):
        """
        Adds molecule to list of unique molecules.
       
        :param molecules:
            The list of unique molecules.
        :param molecule:
            The molecule to be added to list of unique molecules.
        :param jcutoff:
            The inclusion cut-off for similarity index.
        """
       
        rcna = CommonNeighbors(molecule, self.cna_bond)
        rcna.generate(self.cna_rcut)
        
        for mol in molecules:
            lcna = CommonNeighbors(mol, self.cna_bond)
            lcna.generate(self.cna_rcut)
            jval = lcna.comp_cna(rcna)
            if jval > jcutoff:
                return
        
        molecules.append(molecule)

    def print_conformers(self, filename): 
        """                                                                                                                                                                                                                                                                                                     
        Prints conformers data to output stream.
        
        :param filename:
            The name of input file.
        """
        
        self.ostream.print_blank()
        self.ostream.print_header('Conformer Energies')
        self.ostream.print_header(24 * '=')
        self.ostream.print_blank()

        line = '{:>10s}{:>20s} {:>25}{:>30s}'.format(
            'Conformer', 'Energy (a.u.)', 'Relative Energy(a.u.)', 
            'Relative Energy (kcal/mol)') 
        self.ostream.print_header(line) 
        self.ostream.print_header('-' * len(line))
        for i, ene in enumerate(self.conformer_energies):
            rel_ene = ene - self.conformer_energies[0]
            line = '{:5d}{:22.12f}{:22.12f}   {:25.10f}     '.format(
                i + 1, ene, rel_ene, rel_ene * hartree_in_kcalpermol())
            self.ostream.print_header(line) 

        for i, mol in enumerate(self.conformer_geometries):
            fname = filename + '.conformer.{:d}.xyz'.format(i) 
            mol.write_xyz(fname)

    def print_header(self):
        """
        Prints header for the global optimization driver.
        """
        
        self.ostream.print_blank()
        self.ostream.print_header('Global Optimization Driver Setup')
        self.ostream.print_header(34 * '=')
        self.ostream.print_blank()
        
        str_width = 60
        cur_str = 'Maximum Number of Conformers     : {:<6d}'.format(
            self.max_conformers)
        self.ostream.print_header(cur_str.ljust(str_width))
        cur_str = 'Energetic Acceptance Threshold   : {:.2e}'.format(
            self.energy_thresh)
        self.ostream.print_header(cur_str.ljust(str_width))
        cur_str = 'Similarity Acceptance Threshold  : {:.2f}'.format(
            self.similarity_thresh)
        self.ostream.print_header(cur_str.ljust(str_width))
        cur_str = 'Number of Tree Growth Steps      : {:<6d}'.format(
            self.growth_steps)
        self.ostream.print_header(cur_str.ljust(str_width))
        cur_str = 'Tree Growth Unit                 : ' + self.tree_growth_unit.name()
        self.ostream.print_header(cur_str.ljust(str_width))
        cur_str = 'Number of Optimization Steps     : {:<6d}'.format(
            self.popt_steps)
        self.ostream.print_header(cur_str.ljust(str_width))
        cur_str = 'Bond Distance Cut-Off Radius     : {:.2f}'.format(
            self.cna_bond)
        self.ostream.print_header(cur_str.ljust(str_width))
        cur_str = 'Bonds Environment Cut-Off Radius : {:.2f}'.format(
            self.cna_rcut)
        self.ostream.print_header(cur_str.ljust(str_width))
        self.ostream.print_blank()
        
        self.ostream.print_block(self.tree_growth_unit.geom_string())
        self.ostream.print_blank()
        
        self.ostream.flush()
