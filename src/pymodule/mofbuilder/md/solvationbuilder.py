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
from pathlib import Path
try:
    from scipy.spatial.transform import Rotation as R
    from scipy.spatial import cKDTree
except ImportError:
    pass
from ...molecule import Molecule
from ...mmforcefieldgenerator import MMForceFieldGenerator
from ...xtbdriver import XtbDriver
from ...optimizationdriver import OptimizationDriver
from ...molecularbasis import MolecularBasis
from ...scfrestdriver import ScfRestrictedDriver
from ...scfunrestdriver import ScfUnrestrictedDriver
from ..io.basic import nn
from ..core.other import safe_dict_copy
from ...outputstream import OutputStream
from ...veloxchemlib import mpi_master, hartree_in_kcalpermol, hartree_in_kjpermol
from ...errorhandler import assert_msg_critical
from mpi4py import MPI
import sys
import math
from typing import Any, Dict, List, Optional, Tuple, Union

class SolvationBuilder:
    """Builder for the creation of solvated molecular structures in a periodic box.

    This class provides methods to pack solvent molecules around a solute inside
    a simulation box, avoiding atomistic overlaps and managing multi-solvent systems.

    Attributes:
        comm (MPI.Comm): MPI Communicator used for parallel operations.
        rank (int): Rank of the current process in the communicator.
        nodes (int): Number of processes in the communicator.
        ostream (OutputStream): Output stream for logging.
        buffer (float): Minimum distance between atoms of different residues (Å).
        box_size (np.ndarray|list|None): Simulation box size (Å).
        trial_rounds (int): Number of solvation trials to attempt.
        max_fill_rounds (int): Maximum number of cavity-filling iterations per trial.
        scalar (float): Factor to alter number of random solvent candidates per round.
        solute_file (str|None): Filename of solute structure (xyz).
        solute_data (np.ndarray|None): Pre-loaded solute data array.
        solvents_names (List[str]): List of solvent names.
        solvents_files (List[str]): List of solvent structure xyz filenames.
        solvents_proportions (List[float]): List of solvent proportions (sum to 1).
        solvents_quantities (List[int]): List of solvent molecule counts.
        best_solvents_dict (dict|None): Holds best-performing solvent data after solvation.
        custom_solvent_data (dict): Data for user-supplied solvents.
        preferred_region_box (list|None): Sub-box coordinates for focused packing.
        target_directory (str|Path|None): Directory for output files.
        _debug (bool): Toggle debug logging.
        original_solvents_dict (Optional[dict]): Holds original solvents input dict (set in solvate()).
        solute_dict (Optional[dict]): Solute information (labels/coords/n_atoms).
        rc_solute_coords (Optional[np.ndarray]): Solute coords recentered inside box.
        safe_box (Optional[list]): Internal safe box definition used during packing.
        solvents_datalines (Optional[np.ndarray]): Data lines for solvents in output.
        system_datalines (Optional[np.ndarray]): Solute+solvent combined output datalines.
    """

    def __init__(self, comm: Optional[Any] = None, ostream: Optional[OutputStream] = None):
        """Initialize the SolvationBuilder object.

        Args:
            comm (MPI.Comm, optional): MPI communicator.
            ostream (OutputStream, optional): Output stream for logging.
        """
        self.comm = comm or MPI.COMM_WORLD
        self.rank = self.comm.Get_rank()
        self.nodes = self.comm.Get_size()
        self.ostream = ostream or OutputStream(sys.stdout if self.rank ==
                                               mpi_master() else None)

        self.buffer = 1.8  # Å
        self.box_size: Optional[Union[np.ndarray, list]] = None
        self.trial_rounds = 1
        self.max_fill_rounds = 1000  # Maximum number of filling rounds
        self.scalar = 1.0  # Control the number of candidates generated in each round
        self.solute_file: Optional[str] = None
        self.solute_data: Optional[np.ndarray] = None
        self.solvents_names: List[str] = []
        self.solvents_files: List[str] = []
        self.solvents_proportions: List[float] = []
        self.solvents_quantities: List[int] = []
        self.best_solvents_dict: Optional[dict] = None
        self.custom_solvent_data: dict = {}
        self.preferred_region_box: Optional[List[List[float]]] = None
        self.target_directory: Optional[str] = None
        self._debug: bool = False
        # These will be filled after running .solvate()
        self.original_solvents_dict: Optional[dict] = None
        self.solute_dict: Optional[dict] = None
        self.rc_solute_coords: Optional[np.ndarray] = None
        self.safe_box: Optional[list] = None
        self.solvents_datalines: Optional[np.ndarray] = None
        self.system_datalines: Optional[np.ndarray] = None

    def _read_xyz(self, filename: str) -> Tuple[List[str], np.ndarray]:
        """Read atom labels and coordinates from an xyz file.

        Args:
            filename (str): XYZ file path.

        Returns:
            Tuple[List[str], np.ndarray]: Atom labels and coordinates (centered at origin).
        """
        labels, coords = [], []
        with open(filename) as f:
            lines = f.readlines()
        for line in lines[2:]:
            parts = line.split()
            if len(parts) >= 4:
                labels.append(parts[0])
                coords.append(
                    [float(parts[1]),
                     float(parts[2]),
                     float(parts[3])])
        com = np.mean(coords, axis=0)
        coords = np.array(coords) - com  # center at origin
        return labels, coords

    def _generate_candidates_each_solvent(
        self,
        solvent_coords: np.ndarray,
        solvent_labels: List[str],
        solvent_n_atoms: int,
        target_mol_number: int,
        residue_idx_start: int = 0,
        points_template: Optional[np.ndarray] = None,
        box_size: Optional[List[List[float]]] = None,
        rot: bool = True,
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Generate candidate solvent molecule positions, optional rotation/random placement.

        Args:
            solvent_coords (np.ndarray): Coordinates of single solvent molecule.
            solvent_labels (List[str]): Atom labels of solvent.
            solvent_n_atoms (int): Number of atoms per solvent molecule.
            target_mol_number (int): Number of molecules to generate.
            residue_idx_start (int, optional): Starting residue index.
            points_template (np.ndarray, optional): Preferred placement centers.
            box_size (list, optional): Box for candidate placements.
            rot (bool, optional): If true, randomize orientations.

        Returns:
            Tuple: (random_points, candidates, labels, residue_idx)
                - random_points: shape (N, 3) - selected centers
                - candidates: shape (N*solvent_n_atoms, 3) - candidate atom coordinates
                - labels: (N*solvent_n_atoms, 1)
                - residue_idx: (N*solvent_n_atoms, 1)
        """
        if points_template is None:
            random_points = self._box2randompoints(None, box_size, target_mol_number)
        else:
            if points_template.shape[0] < target_mol_number:
                n_additional = target_mol_number - points_template.shape[0]
                random_points = self._box2randompoints(points_template, box_size, n_additional)
            else:
                random_points = points_template
        target_mol_number = random_points.shape[0]
        if target_mol_number == 0:
            return np.empty((0, 3)), np.empty((0, 3)), np.empty((0, 1)), np.empty((0, 1))
        if rot:
            rots = R.random(target_mol_number).as_matrix()
            coords_exp = solvent_coords[np.newaxis, :, :]
            rot_coords = np.matmul(coords_exp, rots.transpose(0, 2, 1))
            candidates = rot_coords.reshape(-1, 3)
        else:
            candidates = np.tile(solvent_coords, (target_mol_number, 1))
        candidates += np.repeat(random_points, solvent_n_atoms, axis=0)

        labels = np.array(list(solvent_labels) * target_mol_number).reshape(-1, 1)
        residue_idx = np.repeat(
            np.arange(residue_idx_start, residue_idx_start + target_mol_number),
            solvent_n_atoms
        ).reshape(-1, 1)

        return random_points, candidates, labels, residue_idx

    def _box2randompoints(
        self,
        points_template: Optional[np.ndarray],
        box_size: List[List[float]],
        n_additional: int
    ) -> np.ndarray:
        """Generate random Cartesian points within a defined box.

        Args:
            points_template (np.ndarray, optional): Existing template points.
            box_size (list): 3D intervals [[xmin, xmax], ...] for axes.
            n_additional (int): Number of new points to generate.

        Returns:
            np.ndarray: Combined (template + new) points (N, 3).
        """
        if points_template is None:
            points_template = np.empty((0, 3))
        additional_points = np.random.rand(n_additional, 3)
        for i in range(3):
            additional_points[:, i] = additional_points[:, i] * (box_size[i][1] - box_size[i][0]) + box_size[i][0]
        random_points = np.vstack((points_template, additional_points))
        return random_points

    def remove_overlaps_kdtree(
        self,
        existing_coords: np.ndarray,
        candidate_coords: np.ndarray,
        candidate_residues: np.ndarray,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Check for overlaps between candidate atoms and existing atoms/residues.

        Args:
            existing_coords (np.ndarray): Atomic coordinates (existing).
            candidate_coords (np.ndarray): Atomic coordinates (candidates).
            candidate_residues (np.ndarray): Residue indices (candidates).

        Returns:
            Tuple[np.ndarray, np.ndarray]: (keep_mask, drop_mask), Booleans for keeping/dropping atoms.
        """
        candidate_residues = candidate_residues.reshape(-1)

        # === Round 1: overlap with existing atoms ===
        assert_msg_critical(
            "scipy" in sys.modules,
            "SciPy is required for MofBuilder.")
        tree_existing = cKDTree(existing_coords)
        dists, _ = tree_existing.query(candidate_coords, k=1, distance_upper_bound=self.buffer)
        mask_overlap_existing = np.isfinite(dists)
        bad_residues_existing = np.unique(candidate_residues[mask_overlap_existing])

        # === Round 2: overlap among candidates ===
        tree_candidates = cKDTree(candidate_coords)
        pairs = np.array(list(tree_candidates.query_pairs(r=self.buffer)))  # (i,j) pairs within radius

        if pairs.size > 0:
            res_i = candidate_residues[pairs[:, 0]]
            res_j = candidate_residues[pairs[:, 1]]
            mask_diff = res_i != res_j
            bad_residues_candidates = np.unique(np.concatenate([res_i[mask_diff], res_j[mask_diff]]))
        else:
            bad_residues_candidates = np.array([], dtype=candidate_residues.dtype)

        bad_residues = np.union1d(bad_residues_existing, bad_residues_candidates)
        keep_mask = ~np.isin(candidate_residues, bad_residues)
        drop_mask = ~keep_mask
        return keep_mask, drop_mask

    def _remove_overlaps_kdtree(
        self,
        existing_coords: np.ndarray,
        candidate_coords: np.ndarray,
        candidate_residues: np.ndarray,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Refined version of remove_overlaps_kdtree using deterministic winner for overlaps.

        Args:
            existing_coords (np.ndarray): Atomic coordinates (existing).
            candidate_coords (np.ndarray): Atomic coordinates (candidates).
            candidate_residues (np.ndarray): Residue indices (candidates).

        Returns:
            Tuple[np.ndarray, np.ndarray]: (keep_mask, drop_mask)
        """
        candidate_residues = candidate_residues.reshape(-1)

        # === Round 1: overlaps with existing atoms ===
        assert_msg_critical(
            "scipy" in sys.modules,
            "SciPy is required for MofBuilder.")
        tree_existing = cKDTree(existing_coords)
        dists, _ = tree_existing.query(candidate_coords, k=1, distance_upper_bound=self.buffer)
        mask_overlap_existing = np.isfinite(dists)
        bad_residues_existing = np.unique(candidate_residues[mask_overlap_existing])

        # === Round 2: residue–residue overlaps (atom-level proximity) ===
        tree_candidates = cKDTree(candidate_coords)
        atom_pairs = np.array(list(tree_candidates.query_pairs(r=self.buffer)))  # (i, j)

        if len(atom_pairs) > 0:
            res_i = candidate_residues[atom_pairs[:, 0]]
            res_j = candidate_residues[atom_pairs[:, 1]]
            diff_mask = res_i != res_j
            residue_pairs = np.unique(np.sort(np.stack([res_i[diff_mask], res_j[diff_mask]], axis=1)), axis=0)
        else:
            residue_pairs = np.empty((0, 2), dtype=candidate_residues.dtype)

        if len(residue_pairs) > 0:
            smaller = np.minimum(residue_pairs[:, 0], residue_pairs[:, 1])
            larger = np.maximum(residue_pairs[:, 0], residue_pairs[:, 1])
            bad_residues_candidates = np.unique(larger)
            keep_residues = np.unique(smaller)
        else:
            bad_residues_candidates = np.array([], dtype=candidate_residues.dtype)
            keep_residues = np.unique(candidate_residues)
        bad_residues = np.union1d(bad_residues_existing, bad_residues_candidates)
        keep_mask = np.isin(candidate_residues, keep_residues) & ~np.isin(candidate_residues, bad_residues)
        drop_mask = ~keep_mask

        return keep_mask, drop_mask

    def _generate_candidates(
        self,
        sol_dict: Dict[str, dict],
        target_number: int,
        res_start: int = 0,
        points_template: Optional[np.ndarray] = None,
        box_size: Optional[List[List[float]]] = None,
        rot: bool = True
    ) -> Tuple[dict, Dict[str, dict], int, np.ndarray]:
        """Generate placement candidates for each solvent, distributed by target counts.

        Args:
            sol_dict (dict): Dictionary of solvent properties.
            target_number (int): Total number of molecules to generate.
            res_start (int, optional): Residue index starting value.
            points_template (np.ndarray, optional): Placement template.
            box_size (list, optional): Placement bounding box.
            rot (bool, optional): Randomize orientions.

        Returns:
            Tuple
                - all_data (dict): coords, labels, residue_idx, atoms_number.
                - sol_dict (updated), res_start (updated), all_res_com_random_points (centers).
        """
        all_data = {}

        solvent_names = list(sol_dict.keys())
        proportions = [sol_dict[name]['proportion'] for name in solvent_names]

        all_sol_mols = self._distribute_by_proportion(target_number, proportions)

        # Guarantee at least 1 molecule if proportion > 0 and total_number >= number of solvents
        for i, name in enumerate(solvent_names):
            if proportions[i] > 0 and all_sol_mols[i] == 0 and target_number >= len(solvent_names):
                all_sol_mols[i] = 1

        diff = target_number - sum(all_sol_mols)
        if diff != 0:
            idx = np.argmax(all_sol_mols)
            all_sol_mols[idx] += diff

        all_sol_atoms_num = [
            n_mol * sol_dict[solvent_name]['n_atoms']
            for n_mol, solvent_name in zip(all_sol_mols, sol_dict)
        ]

        all_data['atoms_number'] = sum(all_sol_atoms_num)
        all_data['coords'] = np.empty((0, 3))
        all_data['labels'] = np.empty((0, 1))
        all_data['residue_idx'] = np.empty((0, 1))
        start_idx = 0
        if target_number == 0:
            return all_data, sol_dict, res_start, np.empty((0, 3))
        if points_template is not None:
            if points_template.shape[0] < target_number:
                n_additional = target_number - points_template.shape[0]
                points_template = self._box2randompoints(
                    points_template, box_size, n_additional)
            if points_template.shape[0] > target_number:
                points_template = points_template[:target_number]
        all_res_com_random_points = np.empty((0, 3))
        for i, solvent_name in enumerate(sol_dict):
            _target_mol_number = all_sol_mols[i]

            com_random_points, candidates, labels, residue_idx = self._generate_candidates_each_solvent(
                sol_dict[solvent_name]['coords'],
                sol_dict[solvent_name]['labels'],
                sol_dict[solvent_name]['n_atoms'],
                _target_mol_number,
                residue_idx_start=res_start,
                points_template=points_template[sum(all_sol_mols[:i])
                                                :sum(all_sol_mols[:i+1])]
                if points_template is not None else None,
                box_size=box_size,
                rot=rot)
            res_com_random_points = np.repeat(
                com_random_points, sol_dict[solvent_name]['n_atoms'], axis=0)
            ex_residue_idx = np.zeros((sum(all_sol_atoms_num), 1), dtype=bool)
            end_idx = start_idx + _target_mol_number * sol_dict[solvent_name]['n_atoms']
            ex_residue_idx[start_idx:end_idx] = True
            start_idx = end_idx

            res_start += _target_mol_number

            sol_dict[solvent_name]['extended_residue_idx'] = np.vstack(
                (sol_dict[solvent_name]['extended_residue_idx'],
                 ex_residue_idx))
            all_res_com_random_points = np.vstack(
                (all_res_com_random_points, res_com_random_points))

            all_data['coords'] = np.vstack((all_data['coords'], candidates))
            all_data['labels'] = np.vstack((all_data['labels'], labels))
            all_data['residue_idx'] = np.vstack(
                (all_data['residue_idx'], residue_idx))

        return all_data, sol_dict, res_start, all_res_com_random_points

    def add_custom_solvent(self, solvent_file: Union[str, List[str]],
                           density: Union[float, List[float]],
                           molar_mass: Union[float, List[float]]) -> None:
        """Add custom solvent definition.

        Args:
            solvent_file (str or List[str]): File name(s) of the solvent xyz structure(s).
            density (float or List[float]): Density/densities in g/cm³.
            molar_mass (float or List[float]): Molecular mass(es) in g/mol.
        """
        if isinstance(solvent_file, str) and Path(solvent_file).is_file():
            self.custom_solvent_data[str(Path(solvent_file).stem)] = {
                "file": solvent_file,
                "density": density,
                "molar_mass": molar_mass
            }
        elif isinstance(solvent_file, list):
            for f, d, m in zip(solvent_file, density, molar_mass):
                self._add_custom_solvent(f, d, m)

    def _get_density_molarmass(self, name: str) -> Tuple[float, float]:
        """Return density (g/cm³) and molar mass (g/mol) for a known solvent name,
        or from custom_solvent_data if provided.

        Args:
            name (str): Solvent name.

        Returns:
            Tuple[float, float]: (density, molar_mass)

        Raises:
            ValueError: If unknown custom solvent and no density/molar mass info.
        """
        name_l = name.lower()
        if name_l in [
                'water', 'tip3p', 'tip4p', 'tip5p', 'tip4pew', 'spc', 'spce'
        ]:
            return 1.0, 18.015  # g/cm³, g/mol
        elif name_l == 'dmso':
            return 1.1, 78.13
        elif name_l == 'methanol':
            return 0.792, 32.04
        elif name_l == 'ethanol':
            return 0.789, 46.07
        elif name_l == 'acetone':
            return 0.784, 58.08
        elif name_l == 'acetonitrile':
            return 0.786, 41.05
        elif name_l == 'chloroform':
            return 1.48, 119.38
        elif name_l == 'dichloromethane':
            return 1.33, 84.93
        elif name_l == 'toluene':
            return 0.866, 92.14
        elif name_l == 'benzene':
            return 0.876, 78.11
        elif name == "CO2":
            return 1.842, 44.01  # supercritical CO2 at 40C and 200 bar
        else:
            self.ostream.print_info(
                f"Unknown solvent {name}, should provide density and molar mass."
            )
            if name in self.custom_solvent_data:
                if self.custom_solvent_data[name][
                        'density'] is not None and self.custom_solvent_data[
                            name]['molar_mass'] is not None:
                    return self.custom_solvent_data[name][
                        'density'], self.custom_solvent_data[name][
                            'molar_mass']
                else:
                    raise ValueError(
                        f"Custom solvent {name} must have both density and molar mass provided."
                    )
            raise ValueError(
                f"Unknown solvent {name} and no custom solvent parameters provided."
            )

    def _xyzfiles2mols(self, solvents_xyz_files: List[str]) -> List[Molecule]:
        """Load a list of xyz files as Molecule objects.

        Args:
            solvents_xyz_files (List[str]): Paths for solvent xyz files.

        Returns:
            List[Molecule]: List of imported Molecule objects.
        """
        solvents_mols = []
        for solvent_file in solvents_xyz_files:
            mol = Molecule.read_xyz_file(solvent_file)
            solvents_mols.append(mol)
        return solvents_mols

    def _initialize_solvents_dict(
        self,
        solvents_files: Optional[List[str]] = [],
        proportion: Optional[List[float]] = [],
        quantities: Optional[List[int]] = [],
    ) -> Optional[Dict[str, dict]]:
        """Prepare solvent input dictionary for later use in the simulation.

        Args:
            solvents_files (List[str], optional): Solvent xyz file paths.
            proportion (List[float], optional): Solvent proportion (will be normalized).
            quantities (List[int], optional): Number of each solvent molecules.

        Returns:
            dict or None: Solvent definition dictionary, or None if invalid input.
        """
        if proportion:
            total_prop = sum(proportion)
            proportion = [p / total_prop for p in proportion]
        elif quantities:
            total_quant = sum(quantities)
            if total_quant == 0:
                return None
            proportion = [q / total_quant for q in quantities]
        else:
            self.ostream.print_warning(
                f"need solvents quantities or proportions")
            return None

        solvents_mols = self._xyzfiles2mols(solvents_files)
        solvents_names = [Path(f).stem for f in solvents_files] if solvents_files else []
        solvents_dict = {}
        for i, solvent_molecule in enumerate(solvents_mols):
            solvent_name = solvents_names[i]
            ds, molm = self._get_density_molarmass(solvent_name)
            if ds is None or molm is None:
                raise ValueError(
                    f"Solvent {solvent_name} must have both density and molar mass provided."
                )
            solvents_dict[solvent_name] = {
                'molecule': solvent_molecule,
                'density': ds,
                'molar_mass': molm,
                'proportion': proportion[i] if proportion else 0.0,
                'target_quantity': quantities[i] if quantities else self._density2number(
                    ds, molm, self.box_size, proportion[i]),
                'labels': solvent_molecule.get_labels(),
                'coords': solvent_molecule.get_coordinates_in_angstrom() - np.mean(
                    solvent_molecule.get_coordinates_in_angstrom(), axis=0),
                'n_atoms': len(solvent_molecule.get_labels()),
                'extended_residue_idx': np.empty((0, 1), dtype=bool),
                'extended_com_points': np.empty((0, 3), dtype=float)
            }
        return solvents_dict

    def molecule_radii(self, coords: np.ndarray) -> float:
        """Estimate molecule radius by half-max distance from origin.

        Args:
            coords (np.ndarray): Coordinates relative to center.

        Returns:
            float: Estimated radius (Å).
        """
        dist_from_center = np.linalg.norm(coords, axis=1)
        max_dist = np.max(dist_from_center) / 2
        return max_dist

    def mols_radii(self, solvents_dict: Dict[str, dict]) -> float:
        """Return representative max radius from all solvents above 0.1 proportion.

        Args:
            solvents_dict (dict): Solvent info dictionary.

        Returns:
            float: Maximum solvent molecular radius (Å).
        """
        radii = []
        for solvent_name in solvents_dict:
            coords = solvents_dict[solvent_name]['coords']
            if solvents_dict[solvent_name]['proportion'] < 0.1:
                continue
            radii.append(self.molecule_radii(coords))
        max_radii = max(radii)
        return max_radii

    def load_solute_info(
        self,
        solute_file: Optional[str] = None,
        solute_data: Optional[np.ndarray] = None
    ) -> dict:
        """Load solute label/coord data from xyz file or pre-loaded array.

        Args:
            solute_file (str, optional): File path of solute xyz.
            solute_data (np.ndarray, optional): Pre-provided solute data.

        Returns:
            dict: Dict with keys 'labels', 'coords', and 'n_atoms'.
        """
        if solute_file is not None and Path(solute_file).is_file():
            solute_labels, solute_coords = self._read_xyz(solute_file)
        elif solute_data is not None:
            solute_data = np.vstack(solute_data).reshape(-1, 11)
            solute_labels = solute_data[:, 1]
            solute_coords = solute_data[:, 5:8].astype(float)
            solute_coords -= np.mean(solute_coords, axis=0)
        else:
            raise ValueError("Either solute_file or solute_data must be provided.")
        solute_info = {
            'labels': solute_labels,
            'coords': solute_coords,
            'n_atoms': len(solute_labels),
        }
        return solute_info

    def grid_points_template(
        self,
        solvents_dict: Dict[str, dict],
        box_size: Union[np.ndarray, list],
        grid_spacing: Optional[float] = None
    ) -> np.ndarray:
        """Generate a cubic lattice of grid points within the simulation box for use as candidate centers.

        Args:
            solvents_dict (dict): Dictionary of solvent properties.
            box_size (np.ndarray | list): Box size (length 3).
            grid_spacing (float, optional): Grid spacing (Å).

        Returns:
            np.ndarray: Grid points shape (N, 3).
        """
        if self.preferred_region_box is not None:
            x_points = np.arange(self.preferred_region_box[0][0] + grid_spacing, self.preferred_region_box[0][1] - grid_spacing,
                                 2 * grid_spacing)
            y_points = np.arange(self.preferred_region_box[1][0] + grid_spacing, self.preferred_region_box[1][1] - grid_spacing,
                                 2 * grid_spacing)
            z_points = np.arange(self.preferred_region_box[2][0] + grid_spacing, self.preferred_region_box[2][1] - grid_spacing,
                                 2 * grid_spacing)
        else:
            x_points = np.arange(0 + grid_spacing, box_size[0] - grid_spacing, 2 * grid_spacing)
            y_points = np.arange(0 + grid_spacing, box_size[1] - grid_spacing, 2 * grid_spacing)
            z_points = np.arange(0 + grid_spacing, box_size[2] - grid_spacing, 2 * grid_spacing)
        xx, yy, zz = np.meshgrid(x_points, y_points, z_points)
        points_template = np.vstack([xx.ravel(), yy.ravel(), zz.ravel()]).T
        if self._debug:
            self.ostream.print_info(
                f"Generated {points_template.shape[0]} template points for solvent placement with grid spacing {grid_spacing:.2f} Å."
            )
            self.ostream.flush()
        return points_template

    def _distribute_by_proportion(self, total_number: int, proportions: List[float]) -> np.ndarray:
        """Distribute integer molecule counts by a list of proportions, rounding but preserving sum.

        Args:
            total_number (int): Total number to distribute.
            proportions (List[float]): Relative proportions (not necessarily normalized).

        Returns:
            np.ndarray: Integer counts per solvent.
        """
        proportions = np.array(proportions, dtype=float)

        if proportions.sum() == 0:
            return np.zeros_like(proportions, dtype=int)

        proportions = proportions / proportions.sum()

        raw = proportions * total_number
        base = np.floor(raw).astype(int)

        remainder = total_number - base.sum()

        if remainder > 0:
            fractional = raw - base
            order = np.argsort(-fractional)
            for i in range(remainder):
                base[order[i]] += 1

        return base

    def solvate(self) -> Optional[dict]:
        """Perform the primary packing/solvation operation.

        Places all solvent molecules into the simulation box around the solute, avoiding overlaps and
        preserving the desired ratios/densities.

        Returns:
            dict or None: best_solvents_dict, with info on final solvent placements/atoms/etc.
        """
        original_solvents_dict = self._initialize_solvents_dict(
            self.solvents_files, self.solvents_proportions,
            self.solvents_quantities)
        if original_solvents_dict is None:
            return
        if self.solute_data is not None:
            solute_dict = self.load_solute_info(solute_data=self.solute_data)
        elif self.solute_file is not None:
            solute_dict = self.load_solute_info(solute_file=self.solute_file)
        else:
            raise ValueError("Either solute_file or solute_data must be provided.")

        self.original_solvents_dict = original_solvents_dict
        self.solute_dict = solute_dict

        total_number = sum([
            original_solvents_dict[solvent_name]['target_quantity']
            for solvent_name in original_solvents_dict
        ])
        self.ostream.print_info(
            f"Total target solvent mols to add: {total_number}")
        self.ostream.flush()
        trial_rounds = max(1, self.trial_rounds)

        if total_number == 0:
            self.ostream.print_info("No solvents to add.")
            return

        if self.box_size is None:
            solute_radius = self.molecule_radii()
            self.box_size = np.array([solute_radius * 2] * 3)

        best_accepted_coords = None
        best_accepted_labels = None
        best_accepted_residues = None
        max_added = 0
        residue_idx = 0

        grid_spacing = self.mols_radii(original_solvents_dict) + self.buffer
        points_template = self.grid_points_template(original_solvents_dict,
                                                    self.box_size,
                                                    grid_spacing=grid_spacing)
        self.safe_box = [[grid_spacing, self.box_size[0] - grid_spacing],
                         [grid_spacing, self.box_size[1] - grid_spacing],
                         [grid_spacing, self.box_size[2] - grid_spacing]]
        self.rc_solute_coords = solute_dict['coords'] + np.array(
            self.box_size) / 2

        # --- Trial loop for random seeds ---
        for trial in range(trial_rounds):
            if self._debug:
                self.ostream.print_info(
                    f"Starting trial {trial+1}/{trial_rounds}.")
                self.ostream.flush()

            candidates_res_idx = np.empty(0)
            np.random.seed(trial)

            solvents_dict = safe_dict_copy(original_solvents_dict)

            all_candidates_data, solvents_dict, res_start_idx, all_res_com_points = self._generate_candidates(
                solvents_dict,
                points_template.shape[0],
                res_start=0,
                points_template=points_template,
                box_size=self.safe_box,
                rot=True)

            all_candidate_coords = all_candidates_data['coords'].astype(float)
            all_candidate_labels = all_candidates_data['labels']
            all_candidate_residues = all_candidates_data['residue_idx'].astype(int)
            candidates_res_idx = np.r_[candidates_res_idx,
                                       all_candidate_residues.flatten()]
            residue_idx += total_number

            keep_masks = np.empty((0), dtype=bool)

            keep_mask, drop_mask = self._remove_overlaps_kdtree(
                self.rc_solute_coords, all_candidate_coords,
                all_candidate_residues)

            accepted_coords = all_candidate_coords[keep_mask]
            accepted_labels = all_candidate_labels[keep_mask]
            accepted_residues = all_candidate_residues[keep_mask]

            if self._debug:
                self.ostream.print_info(
                    f"Trial {trial+1} initial: {len(set(accepted_residues.flatten()))} added, {len(set(all_candidate_residues[drop_mask].flatten()))} left in cavity."
                )
                self.ostream.flush()
            keep_masks = np.r_[keep_masks, keep_mask]

            max_fill_rounds = self.max_fill_rounds
            round_idx = 0
            round_drop_mask = None

            cavity_number = total_number - len(set(
                accepted_residues.flatten()))
            while round_idx < max_fill_rounds and cavity_number > 0:
                round_idx += 1
                if self._debug:
                    self.ostream.print_info(
                        f"Starting fill round {round_idx} with {cavity_number} mols to fill."
                    )
                    self.ostream.flush()

                if cavity_number == 0:
                    break
                new_points_num = cavity_number if cavity_number > 2000 else max(
                    1000, cavity_number)
                round_all_candidates_data, solvents_dict, _, all_res_com_points = self._generate_candidates(
                    solvents_dict,
                    target_number=new_points_num,
                    res_start=res_start_idx,
                    points_template=None,
                    box_size=self.safe_box)

                res_start_idx += new_points_num

                if self._debug:
                    self.ostream.print_info(
                        f"Round {round_idx}: start kdtree overlap removal...")
                    self.ostream.flush()
                round_keep_mask, round_drop_mask = self._remove_overlaps_kdtree(
                    np.vstack((self.rc_solute_coords, accepted_coords)),
                    round_all_candidates_data['coords'],
                    round_all_candidates_data['residue_idx'])
                if self._debug:
                    self.ostream.print_info(
                        f"Round {round_idx}: kdtree overlap removal done.")
                    self.ostream.flush()

                candidates_res_idx = np.r_[
                    candidates_res_idx,
                    round_all_candidates_data['residue_idx'].flatten()]

                keep_masks = np.r_[keep_masks, round_keep_mask]

                round_keep_coords = round_all_candidates_data['coords'][
                    round_keep_mask]
                round_keep_labels = round_all_candidates_data['labels'][
                    round_keep_mask]
                round_keep_residues = round_all_candidates_data['residue_idx'][
                    round_keep_mask]

                round_drop_residues = round_all_candidates_data['residue_idx'][
                    round_drop_mask]

                if self._debug:
                    self.ostream.print_info(
                        f"Round {round_idx}: {len(set(round_keep_residues.flatten()))} added, {len(set(round_drop_residues.flatten()))} left in cavity."
                    )
                self.ostream.flush()

                keep_res_num = len(set(round_keep_residues.flatten()))
                cavity_number -= keep_res_num

                if cavity_number == 0:
                    break

                accepted_coords = np.vstack(
                    (accepted_coords, round_keep_coords))
                accepted_labels = np.r_[accepted_labels, round_keep_labels]
                accepted_residues = np.r_[accepted_residues,
                                          round_keep_residues]

            if accepted_coords.shape[0] > max_added:
                max_added = accepted_coords.shape[0]
                best_accepted_coords = accepted_coords.copy()
                best_accepted_labels = accepted_labels.copy()
                best_accepted_residues = accepted_residues.copy()
                best_accepted_total_number = len(
                    set(best_accepted_residues.flatten()))
                best_solvents_dict = safe_dict_copy(solvents_dict)
                best_keep_masks = keep_masks.copy()
                best_candidates_res_idx = candidates_res_idx.copy()

        if best_accepted_coords is not None:
            self.ostream.print_info(
                f"Best trial added {best_accepted_total_number} mols, atoms: {best_accepted_coords.shape[0]}."
            )
            self.ostream.flush()
            accepted_proportions = []
            target_proportions = []
            proportion_diff = []
            total_number_limit = total_number
            overshoot_flag = False

            for solvent_name in best_solvents_dict:
                accepted_atoms_number = best_solvents_dict[solvent_name][
                    'extended_residue_idx'][best_keep_masks].sum()
                accepted_quantity = accepted_atoms_number // best_solvents_dict[
                    solvent_name]['n_atoms']
                overshoot_flag = (accepted_quantity
                                  > best_solvents_dict[solvent_name]
                                  ['target_quantity']) or overshoot_flag
                best_solvents_dict[solvent_name][
                    'accepted_atoms_number'] = accepted_atoms_number
                best_solvents_dict[solvent_name][
                    'accepted_quantity'] = accepted_quantity
                best_solvents_dict[solvent_name]['accepted_mols_ind'] = (
                    best_solvents_dict[solvent_name]['extended_residue_idx']
                ).ravel() & best_keep_masks.ravel()
                solvent_residues = best_candidates_res_idx[
                    best_solvents_dict[solvent_name]['accepted_mols_ind']]
                resid_mask = np.isin(best_accepted_residues,
                                     np.unique(solvent_residues))
                best_solvents_dict[solvent_name][
                    'accepted_atoms_coords'] = best_accepted_coords[
                        resid_mask.flatten()]
                best_solvents_dict[solvent_name][
                    'accepted_atoms_labels'] = best_accepted_labels[
                        resid_mask.flatten()]
                best_solvents_dict[solvent_name][
                    'accepted_proportion'] = round(accepted_quantity / len(
                        set(best_accepted_residues.flatten())), 8)
                accepted_proportions.append(
                    best_solvents_dict[solvent_name]['accepted_proportion'])
                target_proportions.append(
                    best_solvents_dict[solvent_name]['proportion'])
                proportion_diff.append(
                    best_solvents_dict[solvent_name]['accepted_proportion'] -
                    best_solvents_dict[solvent_name]['proportion'])
                total_number_limit = min(
                    total_number_limit,
                    math.ceil(accepted_quantity /
                        best_solvents_dict[solvent_name]['proportion']
                        if best_solvents_dict[solvent_name]['proportion'] >
                        0 else total_number))
                best_solvents_dict[solvent_name][
                    'accepted_density'] = self._number2density(
                        accepted_quantity,
                        best_solvents_dict[solvent_name]['molar_mass'],
                        self.box_size)
                self.ostream.print_info("*" * 80)
                self.ostream.print_info(
                    f"total number limit after checking each solvent: {total_number_limit}"
                )
                self.ostream.flush()
            self.ostream.print_info("*" * 80)
            self.ostream.print_info(
                f"total number limit after checking each solvent: {total_number_limit}"
            )
            self.ostream.flush()
            if overshoot_flag:
                for solvent_name in best_solvents_dict:
                    limited_quantity = math.ceil(
                        total_number_limit *
                        best_solvents_dict[solvent_name]['proportion'])
                    if best_solvents_dict[solvent_name][
                            'accepted_quantity'] > limited_quantity:
                        overshoot_number = (
                            accepted_atoms_number - limited_quantity *
                            best_solvents_dict[solvent_name]['n_atoms'])
                        if self._debug:
                            self.ostream.print_info(
                                f"Overshoot {solvent_name}: will kick {overshoot_number} atoms."
                            )

                        best_solvents_dict[solvent_name][
                            'accepted_atoms_number'] = limited_quantity * best_solvents_dict[
                                solvent_name]['n_atoms']
                        best_solvents_dict[solvent_name][
                            'accepted_quantity'] = limited_quantity
                        best_solvents_dict[solvent_name][
                            'accepted_atoms_labels'] = best_solvents_dict[
                                solvent_name][
                                    'accepted_atoms_labels'][:best_solvents_dict[
                                        solvent_name]['accepted_atoms_number']]
                        best_solvents_dict[solvent_name][
                            'accepted_atoms_coords'] = best_solvents_dict[
                                solvent_name][
                                    'accepted_atoms_coords'][:best_solvents_dict[
                                        solvent_name]['accepted_atoms_number']]
                        best_solvents_dict[solvent_name][
                            'accepted_density'] = self._number2density(
                                limited_quantity,
                                best_solvents_dict[solvent_name]['molar_mass'],
                                self.box_size)
                        if self._debug:
                            self.ostream.print_info(
                                f"Kicked {overshoot_number} atoms for {solvent_name}."
                            )

            self.ostream.print_info("*" * 80)
            for solvent_name in best_solvents_dict:
                self.ostream.print_info(
                    f"Final accepted {solvent_name}: {best_solvents_dict[solvent_name]['accepted_quantity']} mols, {best_solvents_dict[solvent_name]['accepted_atoms_number']} atoms."
                )
                self.ostream.print_info(
                    f"Final accepted density of {solvent_name}: {best_solvents_dict[solvent_name]['accepted_density']:.4f} g/cm³. target density: {best_solvents_dict[solvent_name]['proportion'] * best_solvents_dict[solvent_name]['density']:.4f} g/cm³"
                )
            self.ostream.print_info("*" * 80)
            self.ostream.flush()

        else:
            self.ostream.print_warning(
                "No solvent mols were added in any trial.")
            self.ostream.flush()

        self.best_solvents_dict = best_solvents_dict

        return best_solvents_dict

    def _update_datalines(self, res_idx_start: int = 1) -> Tuple[np.ndarray, np.ndarray]:
        """Generate/refresh output datalines for solute and all solvents after packing.

        Args:
            res_idx_start (int): Starting residue idx for solvents.

        Returns:
            Tuple[np.ndarray, np.ndarray]: (solute_data, solvents_datalines)
        """
        if self.solute_data is None:
            self.solute_data = np.empty((self.solute_dict['n_atoms'], 11), dtype=object)
            self.solute_data[:,5:8] = self.solute_dict['coords']
            self.solute_data[:,1] = self.solute_dict['labels']
            self.solute_data[:,0] = self.solute_dict['labels']
            atom_numbers = np.arange(1, self.solute_dict['n_atoms'] + 1)
            self.solute_data[:,2] = atom_numbers
            self.solute_data[:,3] = np.array(['SOLUTE'] * self.solute_dict['n_atoms'])
            residue_numbers = np.repeat(np.arange(res_idx_start, res_idx_start + 1),
                                       self.solute_dict['n_atoms'])
            self.solute_data[:,4] = residue_numbers
            self.solute_data[:,8] = np.array([1] * self.solute_dict['n_atoms'])
            self.solute_data[:,9] = np.array([0] * self.solute_dict['n_atoms'])
            self.solute_data[:,10] = np.array([''] * self.solute_dict['n_atoms'])
            self.solute_data[:,5:8] = self.solute_data[:,5:8].astype(float)
            res_idx_start += 1

        self.solute_data[:, 5:8] = self.solute_data[:, 5:8].astype(float) - np.mean(self.solute_data[:, 5:8].astype(float), axis=0) + np.array(self.box_size) / 2

        solvents_datalines = np.empty((0, 11))
        if self.best_solvents_dict is None:
            self.solvents_datalines = solvents_datalines
            return self.solute_data, self.solvents_datalines

        for solvent, data in self.best_solvents_dict.items():
            labels = np.array(data['accepted_atoms_labels']).reshape(-1, 1)
            coords = np.array(data['accepted_atoms_coords']).reshape(
                -1, 3).astype(float)
            x = coords[:, 0].reshape(-1, 1)
            y = coords[:, 1].reshape(-1, 1)
            z = coords[:, 2].reshape(-1, 1)
            spin = np.zeros(len(labels)).reshape(-1, 1)
            charge = np.zeros(len(labels)).reshape(-1, 1)
            note = np.array([''] * len(labels)).reshape(-1, 1)
            residue_name = np.array([solvent] * len(labels)).reshape(-1, 1)
            atom_number = np.arange(1, len(labels) + 1).reshape(-1, 1)
            residue_number = np.repeat(
                np.arange(res_idx_start, data['accepted_quantity'] + res_idx_start),
                data['n_atoms']).reshape(-1, 1)
            res_idx_start += data['accepted_quantity']

            if len(labels) > 0:
                arr = np.hstack(
                    (labels, labels, atom_number, residue_name, residue_number,
                     x, y, z, spin, charge, note)).reshape(-1, 11)
                self.best_solvents_dict[solvent]['data_lines'] = arr
                solvents_datalines = np.vstack((solvents_datalines, arr))
            else:
                self.best_solvents_dict[solvent]['data_lines'] = []
        self.solvents_datalines = solvents_datalines
        return self.solute_data, self.solvents_datalines

    def write_output(self, output_file: str = "solvated_structure", format: Union[str, List[str]] = []) -> None:
        """Write the solvated system to file(s) in various formats.

        Args:
            output_file (str, optional): Base target name (suffix added per format).
            format (str | list of str, optional): Output format(s): xyz, pdb, gro.
        """
        if self.target_directory is not None:
            output_file = Path(self.target_directory) / output_file
        if not format:
            self.ostream.print_warning(
                "No output format specified, defaulting to 'xyz'.")
            format = ['xyz']
        if isinstance(format, str):
            format = [format]

        self.system_datalines = np.vstack(
            (self.solute_data, self.solvents_datalines))
        header = f"Generated by MofBuilder\n"
        if 'xyz' in format:
            from ..io.xyz_writer import XyzWriter
            xyz_writer = XyzWriter(comm=self.comm, ostream=self.ostream)
            xyz_file = Path(output_file).with_suffix('.xyz')
            xyz_writer.write(filepath=xyz_file,
                             header=header,
                             lines=self.system_datalines)
            self.ostream.print_info(f"Wrote output XYZ file: {xyz_file}")
            self.ostream.flush()
        if 'pdb' in format:
            from ..io.pdb_writer import PdbWriter
            pdb_writer = PdbWriter(comm=self.comm, ostream=self.ostream)
            pdb_file = Path(output_file).with_suffix('.pdb')
            pdb_writer.write(filepath=pdb_file,
                             header=header,
                             lines=self.system_datalines)
            self.ostream.print_info(f"Wrote output PDB file: {pdb_file}")
            self.ostream.flush()
        if 'gro' in format:
            from ..io.gro_writer import GroWriter
            gro_file = Path(output_file).with_suffix('.gro')
            gro_writer = GroWriter(comm=self.comm, ostream=self.ostream)
            gro_writer.write(filepath=gro_file,
                             header=header,
                             lines=self.system_datalines,
                             box=self.box_size)
            self.ostream.print_info(f"Wrote output GRO file: {gro_file}")
            self.ostream.flush()

    def _density2number(
        self,
        density: float,
        molar_mass: float,
        box_size: Union[np.ndarray, list],
        proportion: float = 1.0
    ) -> int:
        """Compute number of molecules for a given density and molar mass in a box.

        Args:
            density (float): Density (g/cm³).
            molar_mass (float): Molecular mass (g/mol).
            box_size (array-like): Simulation box (Å or [Lx, Ly, Lz]).
            proportion (float, optional): Proportion relative to total density.

        Returns:
            int: Number of molecules fitting the space.
        """
        volume_A3 = np.prod(box_size)  # Å³
        volume_cm3 = volume_A3 * 1e-24  # cm³
        N_A = 6.022e23  # Avogadro's number
        n_mols = math.ceil(density * volume_cm3 * N_A * proportion / molar_mass)
        return n_mols

    def _number2density(
        self,
        n_mols: int,
        molar_mass: float,
        box_size: Union[np.ndarray, list]
    ) -> float:
        """Calculate system density for box and known number of moles.

        Args:
            n_mols (int): Number of molecules.
            molar_mass (float): g/mol of the solvent.
            box_size (array-like): Box size (Å³).

        Returns:
            float: Density (g/cm³).
        """
        volume_A3 = np.prod(box_size)  # Å³
        volume_cm3 = volume_A3 * 1e-24  # cm³
        mass_g = n_mols * molar_mass / 6.022e23  # g
        density = mass_g / volume_cm3  # g/cm³
        return density
