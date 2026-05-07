#
#                                   VELOXCHEM
#              ----------------------------------------------------
#                          An Electronic Structure Code
#
#  SPDX-License-Identifier: BSD-3-Clause
#
#  Copyright 2018-2025 VeloxChem developers
#
#  Redistribution and use in source and binary forms, with or without modification,
#  are permitted provided that the following conditions are met:
#
#  1. Redistributions of source code must retain the above copyright notice, this
#     list of conditions and the following disclaimer.
#  2. Redistributions in binary form must reproduce the above copyright notice,
#     this list of conditions and the following disclaimer in the documentation
#     and/or other materials provided with the distribution.
#  3. Neither the name of the copyright holder nor the names of its contributors
#     may be used to endorse or promote products derived from this software without
#     specific prior written permission.
#
#  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
#  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
#  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
#  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
#  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
#  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
#  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
#  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
#  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
#  OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

from typing import Any, List, Optional

import numpy as np
import networkx as nx
from ...outputstream import OutputStream
from ...veloxchemlib import mpi_master
from ...errorhandler import assert_msg_critical
from mpi4py import MPI
import sys
from .other import safe_copy
from ..utils.environment import get_data_path
from ..utils.fetch import fetch_pdbfile
from pathlib import Path
from .linker import FrameLinker
from .write import MofWriter
from ..io.pdb_reader import PdbReader
from ..md.linkerforcefield import LinkerForceFieldGenerator
from ..md.gmxfilemerge import GromacsForcefieldMerger
from ..md.solvationbuilder import SolvationBuilder
from ..visualization.viewer import Viewer
from ..md.openmmsetup import OpenmmSetup
from .defects import TerminationDefectGenerator


class Framework:
    """Holds the built MOF: graph, cell, merged data, and options for writing, solvation, and MD.

    Set by MetalOrganicFrameworkBuilder.build(). Provides write() for file output,
    replace() for defects, and access to solvationbuilder and linker FF settings.

    Attributes:
        comm: MPI communicator.
        rank: MPI rank of this process.
        nodes: MPI size (number of processes).
        ostream: Output stream for logging.
        mofwriter: MofWriter instance for writing PDB/GRO/XYZ/CIF.
        framework_data: Merged Cartesian data for the whole framework (set in write()).
        framework_fcoords_data: Merged fractional data (set in write()).
        solvationbuilder: SolvationBuilder instance for adding solvents.
        solvents: List of solvent names or xyz files.
        solvents_molecules: List of solvent molecules.
        solvents_proportions: List of solvent proportions.
        solvents_quantities: List of solvent quantities.
        solvated_gro_file: Path to solvated GRO file (set after solvation).
        target_directory: Output directory for files.
        data_path: Path to database directory.
        save_files: If True, save intermediate/output files.
        linker_ff_name: Name for linker force field (default "Linker").
        linker_charge: Linker charge for FF generation.
        linker_multiplicity: Linker multiplicity.
        linker_reconnect_drv: Driver for reconnection (e.g. 'xtb').
        linker_reconnect_opt: Whether to run reconnection optimization.
        reconnected_linker_molecule: Molecule after reconnection (set by FF step).
        provided_linker_itpfile: Optional path to pre-made linker .itp file.
        filename: Base filename for output (set in write()).
        resp_charges: Whether to use RESP charges for linker.
        _debug: If True, print extra debug messages.
        mof_family: MOF family name (set by builder).
        node_metal: Metal type string.
        dummy_atom_node: Whether nodes include dummy atoms.
        net_spacegroup: Space group of the net.
        net_cell_info: Cell parameters [a, b, c, alpha, beta, gamma].
        net_unit_cell: 3x3 unit cell matrix.
        node_connectivity: Node connectivity.
        linker_connectivity: Linker connectivity (topic).
        linker_fragment_length: Length of linker fragment.
        node_data: Node atom data.
        dummy_atom_node_dict: Dict of dummy atom counts per node.
        termination_data: Termination atom data.
        termination_X_data: Termination X atoms.
        termination_Y_data: Termination Y atoms.
        frame_unit_cell: 3x3 frame unit cell matrix.
        frame_cell_info: Frame cell parameters.
        graph: Edge graph (eG) after building.
        cleaved_graph: Graph after cleaving/defects.
        graph_index_name_dict: Mapping from index to node/edge name.
        graph_matched_vnode_xind: Matched (node, xind, edge) list.
        supercell_info: Supercell dimensions info.
        supercell: Supercell (nx, ny, nz).
        termination: Whether to use terminations.
        termination_name: Name of termination group.
        add_virtual_edge: Whether to add virtual edges for bridge nodes.
        virtual_edge_max_neighbor: Max neighbors for virtual edge.
        sc_unit_cell: 3x3 supercell unit cell matrix.
        sc_unit_cell_inv: Inverse of sc_unit_cell.
        periodicity: Whether framework is periodic.
        linker_fake_edge: Whether linker is fake (zero-length) edge.
        clean_unsaturated_linkers: Whether to remove unsaturated linkers.
        update_node_termination: Whether to update node terminations after removal.
        unsaturated_linkers: List of unsaturated linker names.
        unsaturated_nodes: List of unsaturated node names.
        matched_vnode_xind: Matched (node, xind, edge) list.
        xoo_dict: Dict mapping X index to O indices (XOO) per node.
        residues_info: Dict of residue names and counts.
        solvents_dict: Dict of solvent info after solvation.
        mlp_type: MLP type for energy minimization (e.g. 'mace').
        mlp_model_path: Path to MLP model file.
    """

    def __init__(
        self,
        comm: Optional[Any] = None,
        ostream: Optional[Any] = None,
    ) -> None:
        self.comm = comm or MPI.COMM_WORLD
        self.rank = self.comm.Get_rank()
        self.nodes = self.comm.Get_size()
        self.ostream = ostream or OutputStream(sys.stdout if self.rank ==
                                               mpi_master() else None)

        self.mofwriter = MofWriter(comm=self.comm, ostream=self.ostream)

        #MD preparation
        self.framework_data = None  #merged data for the whole framework, generated in write()
        self.framework_fcoords_data = None  #merged fractional coordinates for the whole framework, generated in write()

        self.solvationbuilder = SolvationBuilder(comm=self.comm,
                                                 ostream=self.ostream)
        self.solvents = []  #list of solvent names or xyz files
        self.solvents_molecules = []  #list of solvent molecules
        self.solvents_proportions = []  #list of solvent proportions
        self.solvents_quantities = []  #list of solvent quantities
        #will be set
        self.solvated_gro_file = None
        #MD simulation

        #others for output and saving
        self.target_directory = None  #will be set to default path if None
        self.data_path = None  #will be set to default path if None

        self.save_files = False
        self.linker_ff_name = "Linker"
        self.linker_charge = None
        self.linker_multiplicity = None
        self.linker_reconnect_drv = 'xtb'
        self.linker_reconnect_opt = True
        self.reconnected_linker_molecule = None  #will be set after linker force field generation
        self.provided_linker_itpfile = None  #if provided, will map directly
        self.filename = None  #will be set in write()
        self.resp_charges = True  #whether to use resp charges for linker force field generation

        #debug
        self._debug = False

        #will be set by MofBuilder.build()
        self.mof_family = None
        self.node_metal = None
        self.dummy_atom_node = None
        self.net_spacegroup = None
        self.net_cell_info = None
        self.net_unit_cell = None
        self.node_connectivity = None
        self.linker_connectivity = None
        self.linker_fragment_length = None
        self.node_data = None
        self.dummy_atom_node_dict = None
        self.termination_data = None
        self.termination_X_data = None
        self.termination_Y_data = None
        self.frame_unit_cell = None
        self.frame_cell_info = None
        self.graph = None
        self.cleaved_graph = None
        self.graph_index_name_dict = None
        self.graph_matched_vnode_xind = None
        self.supercell_info = None
        self.supercell = None
        self.termination = None
        self.termination_name = None
        self.add_virtual_edge = None
        self.virtual_edge_max_neighbor = None
        self.sc_unit_cell = None
        self.sc_unit_cell_inv = None
        self.periodicity = False  #whether the framework are periodic

        self.linker_fake_edge = False  #TODO: check

        self.clean_unsaturated_linkers = False
        self.update_node_termination = True
        self.unsaturated_linkers = []
        self.unsaturated_nodes = []
        self.matched_vnode_xind = None
        self.xoo_dict = None

        self.residues_info = None  #dictionary of residue name and quantity
        self.solvents_dict = None  #dictionary of solvents info after solvation

        #MLP energy minimization
        self.mlp_type = 'mace'  #default MLP type
        self.mlp_model_path = None  #path to the MLP model file

    def replace(
        self,
        replace_indices: Optional[List[int]] = None,
        new_node_pdbfile: Optional[str] = None,
        new_linker_pdbfile: Optional[str] = None,
        new_linker_molecule: Optional[Any] = None,
    ) -> "Framework":
        """Replace nodes or linkers at given indices with new geometry; update cleaved_graph and merged data."""
        replace_indices = replace_indices or []
        self.defectgenerator = TerminationDefectGenerator(comm=self.comm,
                                                          ostream=self.ostream)
        self.defectgenerator.use_termination = self.termination
        self.defectgenerator.linker_connectivity = self.linker_connectivity
        self.defectgenerator.node_connectivity = self.node_connectivity + self.virtual_edge_max_neighbor if self.add_virtual_edge else self.node_connectivity
        self.defectgenerator._debug = self._debug
        self.defectgenerator.eG_index_name_dict = self.graph_index_name_dict
        self.defectgenerator.eG_matched_vnode_xind = self.graph_matched_vnode_xind
        self.defectgenerator.sc_unit_cell_inv = self.sc_unit_cell_inv
        self.defectgenerator.clean_unsaturated_linkers = self.clean_unsaturated_linkers
        self.defectgenerator.update_node_termination = self.update_node_termination
        self.defectgenerator.unsaturated_linkers = self.unsaturated_linkers
        self.defectgenerator.unsaturated_nodes = self.unsaturated_nodes
        #replace
        if new_node_pdbfile is not None:
            self.new_node_pdbfile = new_node_pdbfile
            #use pdbreader to read the replace node pdb files
            pdbreader = PdbReader(comm=self.comm, ostream=self.ostream)
            self.defectgenerator.new_node_data = pdbreader.read_pdb(
                filepath=self.new_node_pdbfile)
            self.defectgenerator.new_node_X_data = pdbreader.X_data
        if new_linker_pdbfile is not None:
            self.new_linker_pdbfile = new_linker_pdbfile
            #use pdbreader to read the replace linker pdb files
            pdbreader = PdbReader(comm=self.comm, ostream=self.ostream)
            self.defectgenerator.new_linker_data = pdbreader.read_pdb(
                filepath=self.new_linker_pdbfile)
            self.defectgenerator.new_linker_X_data = pdbreader.X_data
        if new_linker_molecule is not None:
            self.new_linker_molecule = new_linker_molecule
            #use the molecule directly
            fr_new_linker = FrameLinker(comm=self.comm, ostream=self.ostream)
            fr_new_linker.linker_connectivity = self.linker_connectivity
            if self.save_files:  #TODO: check if the target directory is set
                fr_new_linker.target_directory = self.target_directory
            fr_new_linker.create(molecule=self.new_linker_molecule)
            #pass linker data
            new_linker_center_data = fr_new_linker.linker_center_data
            new_linker_center_X_data = fr_new_linker.linker_center_X_data
            if fr_new_linker.linker_connectivity > 2:
                #recenter com of out data
                new_linker_com = np.mean(
                    fr_new_linker.linker_outer_X_data[:, 5:8].astype(float),
                    axis=0)
                new_linker_outer_data = np.hstack(
                    (fr_new_linker.linker_outer_data[:, 0:5],
                     fr_new_linker.linker_outer_data[:, 5:8].astype(float) -
                     new_linker_com, fr_new_linker.linker_outer_data[:, 8:]))
                new_linker_outer_X_data = np.hstack(
                    (fr_new_linker.linker_outer_X_data[:, 0:5],
                     fr_new_linker.linker_outer_X_data[:, 5:8].astype(float) -
                     new_linker_com, fr_new_linker.linker_outer_X_data[:, 8:]))
                new_linker_frag_length = np.linalg.norm(
                    new_linker_outer_X_data[0, 5:8].astype(float) -
                    new_linker_outer_X_data[1, 5:8].astype(float))
            else:
                new_linker_frag_length = np.linalg.norm(
                    new_linker_center_X_data[0, 5:8].astype(float) -
                    new_linker_center_X_data[1, 5:8].astype(float))

            self.defectgenerator.new_linker_data = new_linker_center_data
            self.defectgenerator.new_linker_X_data = new_linker_center_X_data

        rpG = self.defectgenerator.replace_items(replace_indices,
                                                 self.graph.copy())
        # create a new Framework object to hold the replaced graph
        new_framework = Framework(comm=self.comm, ostream=self.ostream)
        #inherit the properties from the current framework all attributes except graph
        for attr, value in self.__dict__.items():
            if attr not in ['graph', 'defectgenerator', 'framework_data']:
                setattr(new_framework, attr, safe_copy(value))
        new_framework.graph = rpG.copy()
        new_framework.get_merged_data()
        return new_framework

    def remove(self, remove_indices=[]):
        self.defectgenerator = TerminationDefectGenerator(comm=self.comm,
                                                          ostream=self.ostream)
        self.defectgenerator.use_termination = self.termination
        self.defectgenerator.termination_data = self.termination_data
        self.defectgenerator.termination_X_data = self.termination_X_data
        self.defectgenerator.termination_Y_data = self.termination_Y_data
        self.defectgenerator.linker_connectivity = self.linker_connectivity
        self.defectgenerator.node_connectivity = self.node_connectivity + self.vir_edge_max_neighbor if self.add_virtual_edge else self.node_connectivity
        self.defectgenerator._debug = self._debug
        self.defectgenerator.eG_index_name_dict = self.graph_index_name_dict
        self.defectgenerator.sc_unit_cell = self.sc_unit_cell
        self.defectgenerator.sc_unit_cell_inv = self.sc_unit_cell_inv
        self.defectgenerator.clean_unsaturated_linkers = self.clean_unsaturated_linkers  #boolean
        self.defectgenerator.update_node_termination = self.update_node_termination  #boolean
        self.defectgenerator.matched_vnode_xind = self.matched_vnode_xind
        self.defectgenerator.xoo_dict = self.xoo_dict
        self.defectgenerator.use_termination = self.termination
        self.defectgenerator.unsaturated_linkers = self.unsaturated_linkers
        self.defectgenerator.unsaturated_nodes = self.unsaturated_nodes

        #remove
        rmG = self.defectgenerator.remove_items_or_terminate(
            remove_indices, self.graph.copy())
        #create a new Framework object to hold the removed graph
        new_framework = Framework(comm=self.comm, ostream=self.ostream)
        #inherit the properties from the current framework all attributes except graph
        for attr, value in self.__dict__.items():
            if attr not in ['graph', 'defectgenerator', 'framework_data']:
                setattr(new_framework, attr, safe_copy(value))
        new_framework.graph = rmG.copy()
        new_framework.matched_vnode_xind = safe_copy(
            self.defectgenerator.updated_matched_vnode_xind)
        new_framework.unsaturated_linkers = safe_copy(
            self.defectgenerator.unsaturated_linkers)
        new_framework.unsaturated_nodes = safe_copy(
            self.defectgenerator.updated_unsaturated_nodes)
        new_framework.clean_unsaturated_linkers = self.clean_unsaturated_linkers
        new_framework.get_merged_data()
        return new_framework

    def get_merged_data(self, extra_graph=None):
        self.mofwriter.G = extra_graph if extra_graph is not None else self.graph
        self.mofwriter.frame_cell_info = self.supercell_info
        self.mofwriter.sc_unit_cell = self.sc_unit_cell
        self.mofwriter.xoo_dict = self.xoo_dict
        self.mofwriter.dummy_atom_node_dict = self.dummy_atom_node_dict
        self.mofwriter.target_directory = self.target_directory
        self.mofwriter.supercell_boundary = self.supercell
        self.mofwriter._debug = self._debug
        self.framework_data, self.framework_fcoords_data = self.mofwriter.only_get_merged_data(
        )
        self.residues_info = self.mofwriter.residues_info
        self.linker_molecule_data = self.mofwriter.edges_data[
            0] if self.mofwriter.edges_data else None

    def _boundary_overlap_indices(self):
        boundary_overlap_list = []
        for i in self.unsaturated_linkers + self.unsaturated_nodes:
            if i not in self.graph.nodes:
                continue

            def mean_list(lst):
                arr = np.vstack(lst)
                mean_arr = np.mean(arr, axis=0)
                return mean_arr

            com = mean_list(self.graph.nodes[i]['fcoords'])

            #every point should be less than 1
            def primitive_cell_check(point):
                for i in point:
                    if i >= 1:
                        return False
                return True

            if not primitive_cell_check(com):
                index = self.graph.nodes[i]['index']
                boundary_overlap_list.append(index)

        return boundary_overlap_list

    def _mlp_omm_minimize(self, maxIterations=None):
        try:
            from openmm.app import Simulation, PDBFile
            from openmm import LangevinIntegrator, unit
            from openmmml import MLPotential
            from pathlib import Path
        except ImportError:
            assert_msg_critical(
                False, "openmmml is required for Framework.")
        #write a pdb file #with a random name
        self.write(format='pdb')
        pdb = PDBFile(f"{self.filename}.pdb")
        #delete the temporary pdb file after loading
        Path(f"{self.filename}.pdb").unlink()
        potential = MLPotential(self.mlp_type, modelPath=self.mlp_model_path)
        system = potential.createSystem(pdb.topology)
        integrator = LangevinIntegrator(300 * unit.kelvin, 1 / unit.picosecond,
                                        1 * unit.femtoseconds)
        simulation = Simulation(pdb.topology, system, integrator)
        simulation.context.setPositions(pdb.positions)
        if maxIterations is None:
            simulation.minimizeEnergy()
        else:
            simulation.minimizeEnergy(maxIterations)
        state = simulation.context.getState(getPositions=True,
                                            enforcePeriodicBox=True)
        positions = state.getPositions()
        em_ccoords = np.vstack(positions.value_in_unit(
            unit.nanometers)) * 10  # Convert nm to Angstroms
        em_fcoords = np.dot(em_ccoords, self.sc_unit_cell_inv)
        return em_ccoords, em_fcoords

    def write(
        self,
        format: Optional[List[str]] = None,
        filename: Optional[str] = None,
        periodicity: bool = False,
        mlp_em: bool = False,
        mlp_maxIterations: Optional[int] = None,
    ) -> None:
        """Write framework to file(s). format: list of 'pdb'|'gro'|'xyz'|'cif'; optional MLP energy minimization."""
        format = format or []

        if periodicity or self.periodicity:
            self.ostream.print_info("Writing periodic system")
            self.ostream.flush()
            boundary_overlap_list = self._boundary_overlap_indices()
            if len(boundary_overlap_list) > 0:
                #remove these linkers/nodes that overlap with boundary, but remove_defects() will
                #created a new Framework object, only update the graph to mofwriter
                new_framework = self.remove_defects(
                    remove_indices=boundary_overlap_list)
                self.get_merged_data(extra_graph=new_framework.graph)

        if mlp_em:
            self.ostream.print_info(
                "Performing MLP energy minimization before writing output files..."
            )
            self.ostream.flush()
            #update ccoords and fcoords in mofwriter
            em_ccoords, em_fcoords = self._mlp_omm_minimize(mlp_maxIterations)
            em_merged_data = np.hstack((self.framework_data[:, 0:5],
                                        em_ccoords, self.framework_data[:,
                                                                        8:]))
            em_merged_fcoords = np.hstack(
                (self.framework_fcoords_data[:, 0:5],
                 em_fcoords[:self.framework_fcoords_data.shape[0], :],
                 self.framework_fcoords_data[:, 8:]))
            self.mofwriter.merged_data = em_merged_data
            self.mofwriter.merged_fcoords_data = em_merged_fcoords
            self.framework_data = em_merged_data
            self.framework_fcoords_data = em_merged_fcoords

        self.filename = str(
            Path(filename).parent / Path(filename).stem
        ) if filename is not None else f"{self.mof_family}_mofbuilder_output"
        self.mofwriter.filename = self.filename  #whole path including directory
        self.ostream.print_info(f"Writing output files to {self.filename}.*")
        self.ostream.flush()
        if "xyz" in format:
            self.mofwriter.write_xyz(skip_merge=True)
        if "cif" in format:
            self.mofwriter.write_cif(skip_merge=True,
                                     supercell_boundary=self.supercell,
                                     frame_cell_info=self.supercell_info)
        if "pdb" in format:
            self.mofwriter.write_pdb(skip_merge=True)
        if "gro" in format:
            self.mofwriter.write_gro(skip_merge=True)

    def solvate(self,
                solvents_files=[],
                solvents_proportions=[],
                solvents_quantities=[],
                padding_angstrom=10):
        
        self.solvationbuilder.solvents_files = solvents_files if solvents_files else self.solvents
        #if not provided, use TIP3P as default solvent
        if not self.solvationbuilder.solvents_files:
            self.ostream.print_info("No solvents provided, using TIP3P as default solvent.")
            self.solvents_file = str(Path(self.data_path, 'solvents_database','TIP3P.xyz'))
            self.solvationbuilder.solvents_files = [self.solvents_file]
            solvents_proportions = [1]
        self.solvationbuilder.solute_data = self.framework_data
        #mof box as preferred region
        self.solvationbuilder.preferred_region_box = np.array([[0, self.supercell_info[0]+10], [0, self.supercell_info[1]+10], [0, self.supercell_info[2]+10]]) 
        self.solvationbuilder.solvents_proportions = solvents_proportions if solvents_proportions else self.solvents_proportions
        self.solvationbuilder.solvents_quantities = solvents_quantities if solvents_quantities else self.solvents_quantities
        self.solvationbuilder.target_directory = self.target_directory
        self.solvationbuilder.box_size = self.supercell_info[0:3] + np.array(
            [padding_angstrom, padding_angstrom, padding_angstrom])
        self.solvents_dict = self.solvationbuilder.solvate()
        self.framework_data, self.solvents_data = self.solvationbuilder._update_datalines(
        )

        #update residue info
        def update_residues_info(solvents_dict, residue_info):
            for k, v in solvents_dict.items():
                residue_info[k] = v['accepted_quantity']
            return residue_info

        if self.solvents_dict is not None:
            self.residues_info = update_residues_info(self.solvents_dict,
                                                      self.residues_info)
        self.solvation_system_data = np.vstack(
            (self.framework_data, self.solvents_data))

        #write solvated system to gro file
        file_name = f"{self.mof_family}_in_solvent"
        self.solvationbuilder.write_output(output_file=file_name,
                                           format=["gro"])
        self.solvated_gro_file = str(
            Path(self.target_directory, file_name + ".gro"))

        if self.solvents is None or len(self.solvents) == 0:
            self.solvents = self.solvationbuilder.solvents_files

    def generate_linker_forcefield(self):
        self.linker_ff_gen = LinkerForceFieldGenerator(comm=self.comm,
                                                       ostream=self.ostream)
        self.linker_ff_gen.linker_fake_edge = self.linker_fake_edge
        self.linker_ff_gen.target_directory = self.target_directory
        self.linker_ff_gen.linker_ff_name = self.linker_ff_name if self.linker_ff_name is not None else f"{self.mof_family}_linker"
        self.linker_ff_gen.save_files = self.save_files
        self.linker_ff_gen._debug = self._debug
        self.linker_ff_gen.resp_charges = self.resp_charges
        if self.provided_linker_itpfile is not None:
            self.ostream.print_info(
                "Linker force field is provided by the user, will map it directly."
            )
            self.ostream.flush()
            self.linker_ff_gen.src_linker_molecule = self.src_linker_molecule
            self.linker_ff_gen.src_linker_forcefield_itpfile = self.provided_linker_itpfile
            self.linker_ff_gen.linker_residue_name = "EDG"
            self.linker_ff_gen.map_existing_forcefield(
                self.mofwriter.edges_data[0])
            return

        self.linker_ff_gen.linker_optimization = self.linker_reconnect_opt
        self.linker_ff_gen.linker_residue_name = "EDG"
        self.linker_ff_gen.optimize_drv = self.linker_reconnect_drv  # xtb or qm
        #self.linker_ff_gen.linker_ff_name = self.linker_ff_name if self.linker_ff_name is not None else f"{self.mof_family}_linker"
        self.linker_ff_gen.linker_charge = self.linker_charge if self.linker_charge is not None else -1 * int(
            self.linker_connectivity)
        self.linker_ff_gen.linker_multiplicity = self.linker_multiplicity if self.linker_multiplicity is not None else 1
        self.ostream.print_info(
            f"linker charge is set to {self.linker_ff_gen.linker_charge}")
        self.ostream.print_info(
            f"linker multiplicity is set to {self.linker_ff_gen.linker_multiplicity}"
        )
        self.ostream.flush()
        if self.mofwriter.edges_data:
            self.linker_ff_gen.generate_reconnected_molecule_forcefield(
                self.mofwriter.edges_data[0])
        self.reconnected_linker_molecule = self.linker_ff_gen.dest_linker_molecule

    def md_prepare(self):
        #write gro file for the framework
        self.generate_linker_forcefield()
        self.gmx_ff = GromacsForcefieldMerger()
        self.gmx_ff._debug = self._debug
        self.gmx_ff.solvents_dict = self.solvents_dict
        self.gmx_ff.database_dir = self.data_path if self.data_path is not None else get_data_path(
        )
        self.gmx_ff.target_dir = self.target_directory if self.target_directory is not None else Path.cwd(
        )
        self.gmx_ff.node_metal_type = self.node_metal
        self.gmx_ff.dummy_atom_node = self.dummy_atom_node
        self.gmx_ff.solvents_name = [str(Path(i).stem) for i in self.solvents]

        self.gmx_ff.termination_name = self.termination_name
        self.gmx_ff.linker_itp_dir = self.target_directory
        self.gmx_ff.linker_name = self.linker_ff_gen.linker_ff_name
        self.gmx_ff.residues_info = self.residues_info
        self.gmx_ff.mof_name = self.mof_family
        self.gmx_ff.generate_MOF_gromacsfile()
        if self.solvated_gro_file is None:
            self.ostream.print_warning(f"MOF system is not solvated!")
            self.ostream.flush()
            if self.filename is None:
                self.filename = str(
                    Path(self.target_directory) /
                    f"{self.mof_family}_mofbuilder_output")
            system_pbc = False
            #write gro file for the framework only
            self.ostream.print_info(
                f"Writing gro file for the framework only.")
            self.ostream.flush()
            self.ostream.print_info(
                f"MD input gro file will be {self.filename}.gro")
            self.ostream.flush()
            self.write(format=["gro"], filename=self.filename)
            grofile = self.filename + ".gro"
        else:
            grofile = self.solvated_gro_file
            system_pbc = True
        self.ostream.print_info(
            f"MD input gro file: {grofile}, top file: {self.gmx_ff.top_path}")
        self.ostream.flush()

        #setup MD driver
        self.md_driver = OpenmmSetup(gro_file=grofile,
                                     top_file=self.gmx_ff.top_path,
                                     comm=self.comm,
                                     ostream=self.ostream)
        self.md_driver.system_pbc = system_pbc
        # Run EM + NVT + NPT with single continuous PDB trajectory

    def show(self, w=800, h=600, residue_indices=False, residue_name=False):
        self.viewer = Viewer()
        self.viewer.eG_dict = self.graph_index_name_dict
        self.viewer.merged_lines = self.framework_data
        self.viewer.lines_show(w, h, residue_indices, residue_name)
