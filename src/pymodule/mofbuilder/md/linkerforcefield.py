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

import numpy as np
import sys
import networkx as nx
import re
from typing import (
    Any,
    Dict,
    List,
    Optional,
    Sequence,
    Tuple,
    Union,
)
from pathlib import Path
from ...molecule import Molecule
from ...mmforcefieldgenerator import MMForceFieldGenerator
from ...xtbdriver import XtbDriver
from ...optimizationdriver import OptimizationDriver
from ...molecularbasis import MolecularBasis
from ...scfrestdriver import ScfRestrictedDriver
from ...scfunrestdriver import ScfUnrestrictedDriver
from ..io.basic import nn
from ...outputstream import OutputStream
from ...veloxchemlib import mpi_master, hartree_in_kcalpermol, hartree_in_kjpermol
from ...errorhandler import assert_msg_critical
from mpi4py import MPI
from networkx.algorithms.isomorphism import GraphMatcher


class LinkerForceFieldGenerator:
    """Forcefield generation and mapping for linker molecules in MOFs.

    Handles optimization (with QM/xtb), reintroduction of missing bonds, and
    mapping of existing forcefields onto new linker instances.

    Attributes:
        comm (Any): MPI communicator (defaults to MPI.COMM_WORLD).
        rank (int): MPI process rank.
        nodes (int): Number of MPI processes.
        ostream (OutputStream): Output stream for logging.
        linker_optimization (bool): Whether to optimize linker geometry.
        optimize_drv (str): Optimization method, 'xtb' or 'qm'.
        linker_ff_name (str): Name for linker forcefield.
        linker_residue_name (str): Residue name for linker in topology.
        resp_charges (bool): Whether to use RESP charges.
        linker_fake_edge (bool): If True, treats X-X as fake edges.
        linker_charge (int): Net charge of linker molecule.
        linker_multiplicity (int): Multiplicity for linker molecule.
        target_directory (Optional[str]): Directory where output files will be written.
        save_files (bool): Whether to save output structures.
        src_linker_forcefield_itpfile (Optional[str]): Path to reference linker .itp.
        src_linker_molecule (Optional[Molecule]): Reference linker molecule object.
        dest_linker_molecule (Optional[Molecule]): Reconnected linker molecule (after-building).
        dest_molecule_connectivity_matrix (Optional[np.ndarray]): Connectivity matrix for MOF-linked linker.
        linker_itp_path (Optional[Path]): Path to generated linker .itp file.
        free_opt_linker_mol (Optional[Molecule]): Free optimized linker molecule (QM or xtb).
        _debug (bool): Print debug output if True.

    Methods:
        _reconnect_linker_molecule: Generates a VeloxChem Molecule from input coordinates, rebuilding X-X bonds.
        generate_reconnected_molecule_forcefield: Generates a GROMACS forcefield (.itp) for the reconnected linker.
        _reconnect: Build missing X-X or x-x bonds between linker ends.
        _xtb_optimize: Optimize a molecule using the xtb backend.
        _dft_optimize: Optimize a molecule with DFT.
        _find_isomorphism_and_mapping: Find atom mapping between two isomorphic molecules.
        map_existing_forcefield: Maps existing forcefield (.itp) from source to a reconnected/new linker instance.
    """

    def __init__(self, comm: Any = None, ostream: Optional[OutputStream] = None):
        """
        Initialize the forcefield generator.

        Args:
            comm (Any): MPI communicator, defaults to MPI.COMM_WORLD.
            ostream (Optional[OutputStream]): Output stream for logging.
        """
        self.comm = comm or MPI.COMM_WORLD
        self.rank = self.comm.Get_rank()
        self.nodes = self.comm.Get_size()
        self.ostream = ostream or OutputStream(sys.stdout if self.rank == mpi_master() else None)
        # User-settable options
        self.linker_optimization: bool = True
        self.optimize_drv: str = "xtb"  # "xtb" or "qm"
        self.linker_ff_name: str = "Linker"
        self.linker_residue_name: str = "EDG"
        self.resp_charges: bool = True
        self.linker_fake_edge: bool = False
        self.linker_charge: int = -2
        self.linker_multiplicity: int = 1
        self.target_directory: Optional[str] = None
        self.save_files: bool = False
        # Linker reference and MOF linkage data
        self.src_linker_forcefield_itpfile: Optional[str] = None
        self.src_linker_molecule: Optional[Molecule] = None
        self.dest_linker_molecule: Optional[Molecule] = None
        self.dest_molecule_connectivity_matrix: Optional[np.ndarray] = None
        self.linker_itp_path: Optional[Path] = None
        self.free_opt_linker_mol: Optional[Molecule] = None
        self._debug: bool = False

    def _reconnect_linker_molecule(self, linker_mol_data: np.ndarray) -> Molecule:
        """
        Reconnects dummy atoms (X/x) of a linker, reconstructing bonds suitable for MOF frameworks.

        Args:
            linker_mol_data (np.ndarray): Array representing atom data of the linker.

        Returns:
            Molecule: VeloxChem molecule object with reconnected bonds.

        Raises:
            AssertionError: If linker_mol_data is None.
        """
        xyz_string = ''
        atom_idx = 0
        X_indices_coords = []
        lower_x_indices_coords = []
        atom_coords = []
        assert_msg_critical(
            linker_mol_data is not None,
            "no linker molecule set when trying to generate forcefield for linker"
        )
        lines = linker_mol_data
        xyz_string += f"{len(lines)}\n\n"
        for i in range(len(lines)):
            atom_type = lines[i, 0]
            atom_label = lines[i, 1]
            x = float(lines[i, 5])
            y = float(lines[i, 6])
            z = float(lines[i, 7])
            atom_coords.append((nn(atom_label), [x, y, z]))
            if nn(atom_type) in ['X']:
                X_indices_coords.append((atom_idx, [x, y, z]))
            elif nn(atom_type) in ['x']:
                lower_x_indices_coords.append((atom_idx, [x, y, z]))
            atom_idx += 1
            xyz_string += f"{nn(atom_label)} {x} {y} {z}\n"
        molecule = Molecule.read_xyz_string(xyz_string)
        connectivity_matrix = molecule.get_connectivity_matrix()
        if not self.linker_fake_edge:
            connectivity_matrix, connect_constraints = self._reconnect(X_indices_coords, connectivity_matrix)
            connectivity_matrix, connect_constraints = self._reconnect(lower_x_indices_coords, connectivity_matrix, connect_constraints)
        else:
            all_x_indices_coords = X_indices_coords + lower_x_indices_coords
            connectivity_matrix, connect_constraints = self._reconnect(all_x_indices_coords, connectivity_matrix)
        if self._debug:
            self.ostream.print_info(f"constaints{connect_constraints}")
            self.ostream.flush()
        self.reconnected_connectivity_matrix = connectivity_matrix
        self.reconnected_constraints = connect_constraints
        return molecule

    def generate_reconnected_molecule_forcefield(self, linker_mol_data: Optional[np.ndarray] = None) -> None:
        """
        Generates the forcefield and .itp file for a reconnected linker molecule.

        Args:
            linker_mol_data (Optional[np.ndarray]): Data array for the linker. If None, does nothing.

        Returns:
            None
        """
        if linker_mol_data is None:
            return
        molecule = self._reconnect_linker_molecule(linker_mol_data)
        molecule.set_charge(self.linker_charge)
        molecule.set_multiplicity(self.linker_multiplicity)
        if not self.linker_optimization:
            ff_gen = MMForceFieldGenerator()
            ff_gen.topology_update_flag = True
            ff_gen.ostream.mute()
            ff_gen.connectivity_matrix = self.reconnected_connectivity_matrix
            ff_gen.create_topology(molecule, resp=True)
            self.linker_itp_path = Path(self.target_directory, self.linker_ff_name).with_suffix('.itp')
            ffname = str(self.linker_itp_path).removesuffix('.itp')
            ff_gen.write_gromacs_files(filename=ffname, mol_name=self.linker_residue_name)
        elif self.linker_optimization and self.optimize_drv == "qm":
            constrained_opt_linker_mol, scf_result = self._dft_optimize(molecule, self.reconnected_constraints)
            constrained_opt_linker_mol.set_charge(self.linker_charge)
            constrained_opt_linker_mol.set_multiplicity(self.linker_multiplicity)
            free_opt_linker_mol, scf_result = self._dft_optimize(constrained_opt_linker_mol, None)
            free_opt_linker_mol.set_charge(self.linker_charge)
            free_opt_linker_mol.set_multiplicity(self.linker_multiplicity)
            ff_gen = MMForceFieldGenerator()
            self.linker_itp_path = Path(self.target_directory, self.linker_ff_name).with_suffix('.itp')
            ffname = str(self.linker_itp_path).removesuffix('.itp')
            ff_gen.create_topology(free_opt_linker_mol, resp=True)
            ff_gen.write_gromacs_files(filename=ffname, mol_name=self.linker_residue_name)
            self.free_opt_linker_mol = free_opt_linker_mol
        elif self.linker_optimization and self.optimize_drv == "xtb":
            self.ostream.print_info("xtb driver is using for linker optimization")
            self.ostream.flush()
            constrained_opt_linker_mol = self._xtb_optimize(molecule, self.reconnected_constraints)
            constrained_opt_linker_mol.set_charge(self.linker_charge)
            constrained_opt_linker_mol.set_multiplicity(self.linker_multiplicity)
            free_opt_linker_mol, scf_result = self._dft_optimize(constrained_opt_linker_mol, None)
            free_opt_linker_mol.set_charge(self.linker_charge)
            free_opt_linker_mol.set_multiplicity(self.linker_multiplicity)
            self.free_opt_linker_mol = free_opt_linker_mol
            ff_gen = MMForceFieldGenerator()
            self.ostream.print_info("generating forcefield of linker molecule for Gromacs")
            self.ostream.flush()
            self.linker_itp_path = Path(self.target_directory, Path(self.linker_ff_name).with_suffix('.itp'))
            ffname = str(self.linker_itp_path).removesuffix('.itp')
            ff_gen.create_topology(free_opt_linker_mol, resp=self.resp_charges)
            Path(self.linker_itp_path).parent.mkdir(parents=True, exist_ok=True)
            ff_gen.write_gromacs_files(filename=ffname, mol_name=self.linker_residue_name)

    def _reconnect(
        self,
        X_indices_coords: List[Tuple[int, Sequence[float]]],
        connectivity_matrix: np.ndarray,
        connect_constraints: Optional[List[str]] = None
    ) -> Tuple[np.ndarray, List[str]]:
        """
        Builds (re)connections between linker edge atoms, typically between dummy (X/x) atoms.

        Args:
            X_indices_coords (List[Tuple[int, Sequence[float]]]): Atom indices and coordinates for X/x atoms.
            connectivity_matrix (np.ndarray): Original connectivity matrix.
            connect_constraints (Optional[List[str]]): List of additional constraint strings.

        Returns:
            Tuple[np.ndarray, List[str]]: Updated connectivity matrix and bond constraints (for geometry optimization).
        """
        if connect_constraints is None:
            connect_constraints = []
        half_len_X_num = len(X_indices_coords) // 2
        X1_ind_coords = X_indices_coords[:half_len_X_num]
        X2_ind_coords = X_indices_coords[half_len_X_num:]
        if self._debug:
            self.ostream.print_info(f"X1_indices: {X1_ind_coords}")
            self.ostream.print_info(f"X2_indices: {X2_ind_coords}")
            self.ostream.flush()
        for i, j in zip(X1_ind_coords, X2_ind_coords):
            dist = np.linalg.norm(np.array(i[1]) - np.array(j[1]))
            if dist < 4.0:
                connectivity_matrix[i[0], j[0]] = 1
                connectivity_matrix[j[0], i[0]] = 1
                if self._debug:
                    self.ostream.print_info(f"X-X distance {dist} is within threshold, bond created between indices {i[0]} and {j[0]}.")
                connect_constraints.append(f"set distance {i[0]+1} {j[0]+1} 1.54")
        return connectivity_matrix, connect_constraints

    def _xtb_optimize(self, molecule: Molecule, constraints: Optional[List[str]]) -> Molecule:
        """
        Optimize a molecule using the xtb method (semi-empirical).

        Args:
            molecule (Molecule): Molecule to optimize.
            constraints (Optional[List[str]]): Geometric constraints.

        Returns:
            Molecule: Geometry-optimized molecule.
        """
        xtb_drv = XtbDriver()
        xtb_drv.ostream.mute()
        xtb_drv.compute(molecule)
        opt_drv = OptimizationDriver(xtb_drv)
        opt_drv.conv_energy = 1e-04
        opt_drv.conv_drms = 1e-02
        opt_drv.conv_dmax = 2e-02
        opt_drv.conv_grms = 4e-03
        opt_drv.conv_gmax = 8e-03
        opt_drv.constraints = constraints
        opt_drv.tmax = 0.02
        opt_drv.max_iter = 100
        opt_drv.conv_maxiter = True
        opt_results = opt_drv.compute(molecule)
        self.ostream.print_info("Optimization of linker molecule completed successfully with xtb driver.")
        self.ostream.flush()
        opt_mol = Molecule.read_xyz_string(opt_results["final_geometry"])
        opt_energy = opt_results['opt_energies'][-1] * hartree_in_kcalpermol()
        self.ostream.print_info(f"Optimization energy is {opt_energy} kcal/mol")
        self.ostream.flush()
        if self.save_files:
            fname = str(Path(self.target_directory, "linker_opt.xyz"))
            opt_mol.write_xyz_file(fname)
        return opt_mol

    def _dft_optimize(self, molecule: Molecule, constraints: Optional[List[str]]) -> Tuple[Molecule, Any]:
        """
        Perform DFT optimization for a molecule using BLYP/def2-SVP.

        Args:
            molecule (Molecule): Molecule to optimize.
            constraints (Optional[List[str]]): Geometric constraints.

        Returns:
            Tuple[Molecule, Any]: (Optimized molecule, SCF result object)
        """
        if molecule.get_multiplicity() == 1:
            mol_scf_drv = ScfRestrictedDriver()
        else:
            mol_scf_drv = ScfUnrestrictedDriver()
        basis = MolecularBasis.read(molecule, "def2-svp")
        mol_scf_drv.conv_thresh = 1e-4
        mol_scf_drv.xcfun = "blyp"
        mol_scf_drv.ri_coulomb = True
        mol_scf_drv.grid_level = 4
        mol_scf_drv.ostream.mute()
        mol_scf_results = mol_scf_drv.compute(molecule, basis)
        mol_opt_drv = OptimizationDriver(mol_scf_drv)
        mol_opt_drv.ostream.mute()
        mol_opt_drv.conv_energy = 1e-04
        mol_opt_drv.conv_drms = 1e-02
        mol_opt_drv.conv_dmax = 2e-02
        mol_opt_drv.conv_grms = 4e-03
        mol_opt_drv.conv_gmax = 8e-03
        mol_opt_drv.constraints = constraints
        mol_opt_drv.tmax = 0.02
        mol_opt_drv.max_iter = 200
        opt_results = mol_opt_drv.compute(molecule, basis, mol_scf_results)
        self.ostream.print_info("Optimization of linker molecule completed successfully with DFT.")
        self.ostream.flush()
        opt_mol = Molecule.read_xyz_string(opt_results["final_geometry"])
        opt_energy = opt_results['opt_energies'][-1] * hartree_in_kcalpermol()
        self.ostream.print_info(f"Optimization energy is {opt_energy} kcal/mol")
        self.ostream.flush()
        if self.save_files:
            fname = str(Path(self.target_directory, "linker_opt.xyz"))
            opt_mol.write_xyz_file(fname)
        return opt_mol, mol_scf_results

    def _find_isomorphism_and_mapping(
        self,
        src_mol: Molecule,
        dest_mol: Molecule,
        dest_molecule_connectivity_matrix: Optional[np.ndarray] = None
    ) -> Tuple[bool, Optional[Dict[int, int]]]:
        """
        Determines whether two molecules are isomorphic and returns atom index mapping.

        Args:
            src_mol (Molecule): Source/reference molecule.
            dest_mol (Molecule): Destination/target molecule.

        Returns:
            Tuple[bool, Optional[Dict[int, int]]]: (isomorphic, atom index mapping)
        """
        def correct_connectivity_for_hydrogens(mol: Molecule, conn: np.ndarray) -> np.ndarray:
            # find hydrogens with more than 1 connection and correct connectivity matrix to treat them as non-hydrogens for isomorphism
            #the real connection should be the closest non-hydrogen atom
            labels = mol.get_labels()
            element_ids = mol.get_element_ids()
            for i in range(len(labels)):
                if labels[i].startswith('H') and np.sum(conn[i]) > 1:
                    non_hydrogen_indices = [j for j in range(len(labels)) if not labels[j].startswith('H') and conn[i][j] == 1]
                    if len(non_hydrogen_indices) > 0:
                        #sort connected non-hydrogen atoms by distance to the hydrogen, keep label
                        sorted_non_hydrogen_indices = sorted(non_hydrogen_indices, key=lambda j: np.linalg.norm(np.array(mol.get_coordinates_in_angstrom()[i]) - np.array(mol.get_coordinates_in_angstrom()[j])))
                        #check if the closest non-hydrogen have more than 3 connections, if so, find the closest non-hydrogen with less than 3 connections
                        closest_non_hydrogen_index = sorted_non_hydrogen_indices[0]
                        if np.sum(conn[closest_non_hydrogen_index]) > 3:
                            for idx in sorted_non_hydrogen_indices[1:]:
                                if np.sum(conn[idx]) <= 3:
                                    closest_non_hydrogen_index = idx
                                    break
                        conn[i] = np.zeros_like(conn[i])
                        conn[i][closest_non_hydrogen_index] = 1
                        conn[closest_non_hydrogen_index][i] = 1
            return conn

        def correct_connectivity_for_carbons(mol: Molecule, conn: np.ndarray) -> np.ndarray:
            # find carbons with 4 connections and correct connectivity matrix to treat them as non-carbons for isomorphism
            #drop the farest connection for carbons with 4 connections
            for i in range(len(conn)):
                if np.sum(conn[i]) >3 and mol.get_labels()[i].startswith('C'):
                    non_carbon_indices = [j for j in range(len(conn)) if conn[i][j] == 1]
                    if len(non_carbon_indices) > 0:
                        farthest_non_carbon_index = max(non_carbon_indices, key=lambda j: np.linalg.norm(np.array(mol.get_coordinates_in_angstrom()[i]) - np.array(mol.get_coordinates_in_angstrom()[j])))
                        conn[i][farthest_non_carbon_index] = 0
                        conn[farthest_non_carbon_index][i] = 0
            return conn
        
        def correct_connectivity_for_oxygens(mol: Molecule, conn: np.ndarray) -> np.ndarray:
            # find oxygens with 2 connections and correct connectivity matrix to treat them as non-oxygens for isomorphism
            #drop the farest connection for oxygens with 2 connections
            labels = mol.get_labels()
            for i in range(len(conn)):
                if labels[i].startswith('O') and np.sum(conn[i]) == 2:
                    non_oxygen_indices = [j for j in range(len(conn)) if conn[i][j] == 1]
                    if len(non_oxygen_indices) > 0:
                        farthest_non_oxygen_index = max(non_oxygen_indices, key=lambda j: np.linalg.norm(np.array(mol.get_coordinates_in_angstrom()[i]) - np.array(mol.get_coordinates_in_angstrom()[j])))
                        conn[i][farthest_non_oxygen_index] = 0
                        conn[farthest_non_oxygen_index][i] = 0
            return conn


        def correct_connectivity(mol: Molecule, conn: np.ndarray) -> np.ndarray:
            conn = correct_connectivity_for_hydrogens(mol, conn)
            conn = correct_connectivity_for_oxygens(mol, conn)
            return conn
        
        src_mol_connectivity = correct_connectivity(src_mol, src_mol.get_connectivity_matrix())
        if dest_molecule_connectivity_matrix is not None:
            dest_mol_connectivity = correct_connectivity(dest_mol, dest_molecule_connectivity_matrix)
        else:
            dest_mol_connectivity = correct_connectivity(dest_mol, dest_mol.get_connectivity_matrix() )
        bond_num_src = np.sum(src_mol_connectivity) // 2
        bond_num_dest = np.sum(dest_mol_connectivity) // 2
        if (len(src_mol.get_labels()) != len(dest_mol.get_labels())) or bond_num_src != bond_num_dest:
            raise ValueError(
                f"The two molecules are not isomorphic because of different number of atoms or bonds."
                f"{len(src_mol.get_labels())} atoms and {bond_num_src} bonds in source molecule vs "
                f"{len(dest_mol.get_labels())} atoms and {bond_num_dest} bonds in destination molecule."
            )
        src_mol.show(atom_indices=True)
        dest_mol.show(atom_indices=True)


        
        def get_graph_from_molecule(mol: Molecule, conn: np.ndarray) -> nx.Graph:
            labels = mol.get_labels()
            element_ids = mol.get_element_ids()
            G = nx.Graph()
            for n in range(len(labels)):
                G.add_node(n, atom_id=n, element_id=element_ids[n], label=labels[n])
            for i in range(len(labels)):
                for j in range(i, len(labels)):
                    if conn[i][j]:
                        G.add_edge(i, j)
            return G

        def node_match(n1: dict, n2: dict) -> bool:
            return n1['element_id'] == n2['element_id'] and n1['label'] == n2['label']

        G1 = get_graph_from_molecule(src_mol, src_mol_connectivity)
        G2 = get_graph_from_molecule(dest_mol, dest_mol_connectivity)
        GM = GraphMatcher(G1, G2, node_match=node_match)
        isomorphic = GM.is_isomorphic()
        mapping = next(GM.isomorphisms_iter(), None)
        return isomorphic, mapping

    def map_existing_forcefield(self, linker_mol_data: Optional[np.ndarray] = None) -> None:
        """
        Maps an existing linker forcefield from reference (.itp) onto the built linker in a MOF.

        Args:
            linker_mol_data (Optional[np.ndarray]): Coordinates array of the reconnected linker.

        Returns:
            None
        """
        save_itp_path = Path(self.target_directory, self.linker_ff_name).with_suffix('.itp')
        mapper = ForceFieldMapper(comm=self.comm, ostream=self.ostream)
        mapper.dest_molecule = self._reconnect_linker_molecule(linker_mol_data)
        mapper.dest_molecule_connectivity_matrix = self.reconnected_connectivity_matrix
        mapper.src_molecule_forcefield_itpfile = self.src_linker_forcefield_itpfile
        mapper.src_molecule = self.src_linker_molecule
        mapped_sections = mapper._map_forcefield_sections(dest_resname=self.linker_residue_name)
        mapper.write_mapped_itp_file(mapped_sections, str(save_itp_path))
        self.linker_itp_path = save_itp_path


class ForceFieldMapper:
    """Mapping utility for transferring GROMACS forcefield sections between isomorphic molecules.

    Attributes:
        comm (Any): MPI communicator (defaults to MPI.COMM_WORLD).
        rank (int): MPI process rank.
        nodes (int): Number of MPI processes.
        ostream (OutputStream): Output stream for logging.
        src_molecule_forcefield_itpfile (Optional[str]): Path to source molecule itp file.
        src_molecule (Optional[Molecule]): Source/reference molecule.
        dest_molecule_forcefield_itpfile (Optional[str]): Optional itp for destination.
        dest_molecule (Optional[Molecule]): Target molecule after remapping.
        target_directory (Optional[str]): Output directory.
        save_files (bool): Whether to write files.
        _debug (bool): Enable debug output.
        dest_molecule_connectivity_matrix (Optional[np.ndarray]): For alternate connectivity in destination molecule.

    Methods:
        _get_mapping_between_two_molecules: Atom index mapping between two isomorphic molecules.
        _parse_itp_file: Read .itp file and split into sections.
        _format_*: Format forcefield sections (atoms, bonds, angles, etc.) using index mapping and new labels.
        _map_forcefield_sections: Map/convert all supported forcefield sections from source to destination mapping.
        write_mapped_itp_file: Write mapped sections to output .itp file.
    """

    def __init__(self, comm: Any = None, ostream: Optional[OutputStream] = None):
        """
        Initialize the forcefield mapper.

        Args:
            comm (Any): MPI communicator.
            ostream (Optional[OutputStream]): Output stream for logging.
        """
        self.comm = comm or MPI.COMM_WORLD
        self.rank = self.comm.Get_rank()
        self.nodes = self.comm.Get_size()
        self.ostream = ostream or OutputStream(sys.stdout if self.rank == mpi_master() else None)
        self.src_molecule_forcefield_itpfile: Optional[str] = None
        self.src_molecule: Optional[Molecule] = None
        self.dest_molecule_forcefield_itpfile: Optional[str] = None
        self.dest_molecule: Optional[Molecule] = None
        self.target_directory: Optional[str] = None
        self.save_files: bool = False
        self._debug: bool = False
        self.dest_molecule_connectivity_matrix: Optional[np.ndarray] = None

    def _get_mapping_between_two_molecules(
        self,
        src_molecule: Molecule,
        dest_molecule: Molecule
    ) -> Dict[int, int]:
        """
        Finds atom index mapping between source and destination molecule.

        Args:
            src_molecule (Molecule): The reference molecule.
            dest_molecule (Molecule): The built (possibly reordered) molecule.

        Returns:
            Dict[int, int]: Mapping from source to destination indices (both 1-based).

        Raises:
            ValueError: If molecules are not isomorphic.
        """
        isomorphic, mapping = LinkerForceFieldGenerator()._find_isomorphism_and_mapping(src_molecule, dest_molecule, self.dest_molecule_connectivity_matrix)
        if not isomorphic or mapping is None:
            raise ValueError(
                "The linker molecule in MOF is not isomorphic to the reference linker molecule."
            )
        mapping = dict(sorted(mapping.items()))
        mapping = {k + 1: v + 1 for k, v in mapping.items()}
        if self._debug:
            self.ostream.print_info(f"Mapping between source and destination linker molecule: {mapping}")
            self.ostream.flush()
        return mapping

    def _parse_itp_file(self, itpfile: str) -> Dict[str, List[str]]:
        """
        Splits a .itp file into its relevant sections (atoms, bonds, etc.).

        Args:
            itpfile (str): Path to .itp file.

        Returns:
            Dict[str, List[str]]: Dict of section name to its raw lines.

        Raises:
            AssertionError: If file path is None or unreadable.
        """
        assert_msg_critical(itpfile is not None, "No source forcefield itp file provided.")
        sections = {}
        with open(itpfile, 'r') as fp:
            number = []
            lineNumber = 1
            keyword = "]"
            lines = fp.readlines()
            for eachline in lines:
                m = re.search(keyword, eachline)
                if m is not None:
                    number.append(lineNumber - 1)
                lineNumber += 1
            number.append(len(lines))
            number = list(set(number))
            number.sort()
            size = int(len(number))
            for i in range(size - 1):
                start = number[i]
                end = number[i + 1]
                middlelines = lines[start:end]
                middlelines = [line for line in middlelines if line.strip()]
                section = re.findall(r'\[(.*?)\]', middlelines[0])
                title = section[0].split()[0]
                if title == 'dihedrals' and len(middlelines) > 1 and 'impropers' in middlelines[1]:
                    title = 'dihedrals_im'
                sections[title] = middlelines
        if self._debug:
            self.ostream.print_info(f"Sections found in itp file: {list(sections.keys())}")
            self.ostream.flush()
        return sections

    def _format_atomtypes(self, lines: List[str]) -> List[str]:
        """
        Returns the atomtypes section, unchanged.

        Args:
            lines (List[str]): Lines from atomtypes section.

        Returns:
            List[str]: Unchanged lines.
        """
        return lines

    def _format_moleculetype(self, lines: List[str], new_resname: str) -> List[str]:
        """
        Formats the moleculetype section to update with a new residue name.

        Args:
            lines (List[str]): Lines from the section.
            new_resname (str): New residue name.

        Returns:
            List[str]: Updated lines for moleculetype section.
        """
        newff = []
        newff.append(lines[0])
        newff.append(lines[1])
        values = lines[2].split()
        values[0] = new_resname
        formatted_line = "%-7s%7s" % (values[0], values[1])
        newff.append(formatted_line + "\n")
        return newff

    def _format_atoms(
        self,
        lines: List[str],
        src_dest_map: Dict[int, int],
        dest_atomlabels: List[str],
        dest_resname: str
    ) -> List[str]:
        """
        Remaps atom indices/names for atoms section to new linker atom order.

        Args:
            lines (List[str]): Source section lines.
            src_dest_map (Dict[int, int]): Mapping from src to dest (1-based indices).
            dest_atomlabels (List[str]): Atom labels in destination molecule.
            dest_resname (str): Target residue name.

        Returns:
            List[str]: Remapped lines for [atoms] section.
        """
        new_atoms_section = [lines[0], lines[1]]
        for i in range(2, len(lines)):
            values = lines[i].split()
            if len(values) == 0:
                continue
            src_index = int(values[0])
            dest_index = src_dest_map[src_index]
            values[4] = f"{dest_atomlabels[dest_index-1]}{dest_index}"
            values[5] = dest_index
            values[0] = dest_index
            values[6] = float(values[6])
            values[7] = float(values[7])
            values[3] = dest_resname
            if len(values) > 8:
                formatted_line = "%7d%7s%7s%7s%7s%7d%15.8f%15.6f%7s" % (
                    values[0], values[1], values[2], values[3], values[4], values[5],
                    values[6], values[7], values[8])
                new_atoms_section.append(formatted_line)
            else:
                formatted_line = "%7s%7s%7s%7s%7s%7s%15.8f%15.6f" % (
                    values[0], values[1], values[2], values[3], values[4], values[5],
                    values[6], values[7])
                new_atoms_section.append(formatted_line)
        header = new_atoms_section[:2]
        body = sorted(new_atoms_section[2:], key=lambda i: int(i.split()[0]))
        new_main = []
        charge_sum = 0.0
        for line in body:
            values = line.split()
            charge_sum += float(values[6])
            charge_sum_line = line + f" ;qtol  {charge_sum:15.6f}\n"
            new_main.append(charge_sum_line)
        return header + new_main

    def _format_bonds(self, lines: List[str], src_dest_map: Dict[int, int]) -> List[str]:
        """
        Formats the bonds section with remapped indices.

        Args:
            lines (List[str]): Bonds section.
            src_dest_map (Dict[int, int]): Mapping from src to dest.

        Returns:
            List[str]: Remapped bonds lines.
        """
        new_bonds_section = [lines[0], lines[1]]
        for i in range(2, len(lines)):
            values = lines[i].split()
            if len(values) == 0:
                continue
            src_bond_i_index = int(values[0])
            src_bond_j_index = int(values[1])
            values[0] = src_dest_map[src_bond_i_index]
            values[1] = src_dest_map[src_bond_j_index]
            values[3] = float(values[3])
            values[4] = float(values[4])
            formatted_line = "%7s%7s%6s%15.7f%15.6f" % (
                values[0], values[1], values[2], values[3], values[4])
            new_bonds_section.append(formatted_line + "\n")
        return new_bonds_section

    def _format_pairs(self, lines: List[str], src_dest_map: Dict[int, int]) -> List[str]:
        """
        Formats the pairs section with remapped indices.

        Args:
            lines (List[str]): Source pairs section.
            src_dest_map (Dict[int, int]): Atom index mapping.

        Returns:
            List[str]: Remapped pairs lines.
        """
        new_pairs_section = [lines[0], lines[1]]
        for i in range(2, len(lines)):
            values = lines[i].split()
            if len(values) == 0:
                continue
            src_pair_i_index = int(values[0])
            src_pair_j_index = int(values[1])
            values[0] = src_dest_map[src_pair_i_index]
            values[1] = src_dest_map[src_pair_j_index]
            formatted_line = "%7s%7s%6s" % (values[0], values[1], values[2])
            new_pairs_section.append(formatted_line + "\n")
        return new_pairs_section

    def _format_angles(self, lines: List[str], src_dest_map: Dict[int, int]) -> List[str]:
        """
        Formats the angles section with remapped indices.

        Args:
            lines (List[str]): Source angles section.
            src_dest_map (Dict[int, int]): Mapping from src to dest.

        Returns:
            List[str]: Remapped angles lines.
        """
        new_angles_section = [lines[0], lines[1]]
        for i in range(2, len(lines)):
            values = lines[i].split()
            if len(values) == 0:
                continue
            src_angle_i_index = int(values[0])
            src_angle_j_index = int(values[1])
            src_angle_k_index = int(values[2])
            values[0] = src_dest_map[src_angle_i_index]
            values[1] = src_dest_map[src_angle_j_index]
            values[2] = src_dest_map[src_angle_k_index]
            values[4] = float(values[4])
            values[5] = float(values[5])
            formatted_line = "%7s%7s%7s%6s%13.7f%12.6f" % (
                values[0], values[1], values[2], values[3], values[4], values[5])
            new_angles_section.append(formatted_line + "\n")
        return new_angles_section

    def _format_dihedrals(self, lines: List[str], src_dest_map: Dict[int, int]) -> List[str]:
        """
        Formats the dihedrals section using remapped indices.

        Args:
            lines (List[str]): Source dihedrals (proper or improper).
            src_dest_map (Dict[int, int]): Mapping from src to dest.

        Returns:
            List[str]: Remapped dihedrals lines.
        """
        new_dihedrals_section = [lines[0], lines[1], lines[2]]
        for i in range(3, len(lines)):
            values = lines[i].split()
            if len(values) == 0:
                continue
            src_dihedral_i_index = int(values[0])
            src_dihedral_j_index = int(values[1])
            src_dihedral_k_index = int(values[2])
            src_dihedral_l_index = int(values[3])
            values[0] = src_dest_map[src_dihedral_i_index]
            values[1] = src_dest_map[src_dihedral_j_index]
            values[2] = src_dest_map[src_dihedral_k_index]
            values[3] = src_dest_map[src_dihedral_l_index]
            values[5] = float(values[5])
            values[6] = float(values[6])
            formatted_line = "%7s%7s%7s%7s%6s%13.7f%12.7f%3s" % (
                values[0], values[1], values[2], values[3], values[4], values[5], values[6], values[7])
            new_dihedrals_section.append(formatted_line + "\n")
        return new_dihedrals_section

    def _map_forcefield_sections(self, dest_resname: str = "MOL") -> Dict[str, List[str]]:
        """
        Remap and convert all relevant sections from a reference .itp onto the destination molecule.

        Args:
            dest_resname (str): Target residue name for mapping.

        Returns:
            Dict[str, List[str]]: New section contents for .itp
        """
        src_dest_map = self._get_mapping_between_two_molecules(self.src_molecule, self.dest_molecule)
        self.ostream.print_info(f"mapping: {src_dest_map}")
        dest_atomlabels = self.dest_molecule.get_labels()
        sections = self._parse_itp_file(self.src_molecule_forcefield_itpfile)
        if self._debug:
            self.ostream.print_info(f"Forcefield sections found in source itp file: {list(sections.keys())}")
            self.ostream.flush()
        dest_sections = {}
        if 'atomtypes' in sections:
            dest_sections['atomtypes'] = self._format_atomtypes(sections['atomtypes'])
        if 'moleculetype' in sections:
            dest_sections['moleculetype'] = self._format_moleculetype(sections['moleculetype'], dest_resname)
        if 'atoms' in sections:
            dest_sections['atoms'] = self._format_atoms(sections['atoms'], src_dest_map, dest_atomlabels, dest_resname)
        if 'bonds' in sections:
            dest_sections['bonds'] = self._format_bonds(sections['bonds'], src_dest_map)
        if 'pairs' in sections:
            dest_sections['pairs'] = self._format_pairs(sections['pairs'], src_dest_map)
        if 'angles' in sections:
            dest_sections['angles'] = self._format_angles(sections['angles'], src_dest_map)
        if 'dihedrals' in sections:
            dest_sections['dihedrals'] = self._format_dihedrals(sections['dihedrals'], src_dest_map)
        if 'dihedrals_im' in sections:
            dest_sections['dihedrals_im'] = self._format_dihedrals(sections['dihedrals_im'], src_dest_map)
        if self._debug:
            self.ostream.print_info(f"Forcefield sections mapped from source to destination molecule.")
            self.ostream.print_info(f"Sections in mapped itp file: {list(dest_sections.keys())}")
            self.ostream.print_info(f"mapping {src_dest_map}")
            self.ostream.print_info(f"destination atom labels {dest_atomlabels}")
            self.ostream.flush()
        self.mapped_sections = dest_sections
        return dest_sections

    def write_mapped_itp_file(self, mapped_sections: Dict[str, List[str]], output_itpfile: Union[str, Path]) -> None:
        """
        Writes a mapped .itp file from provided section line mappings.

        Args:
            mapped_sections (Dict[str, List[str]]): Section lines by section name.
            output_itpfile (Union[str, Path]): Destination path (.itp).

        Returns:
            None
        """
        if not Path(output_itpfile).suffix == '.itp':
            output_itpfile = Path(output_itpfile).with_suffix('.itp')
        self.ostream.print_info(f"Writing mapped forcefield to itp file: {output_itpfile}")
        self.ostream.flush()
        #create parent directory if not exist
        Path(output_itpfile).parent.mkdir(parents=True, exist_ok=True)
        with open(output_itpfile, 'w') as fp:
            for section_name, lines in mapped_sections.items():
                for line in lines:
                    fp.write(line)
                fp.write("\n")
