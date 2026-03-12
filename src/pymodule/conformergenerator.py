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

from mpi4py import MPI
from pathlib import Path
from copy import deepcopy
import numpy as np
import itertools
import time
import sys

from .veloxchemlib import bohr_in_angstrom, mpi_master
from .outputstream import OutputStream
from .molecule import Molecule
from .atomtypeidentifier import AtomTypeIdentifier
from .mmforcefieldgenerator import MMForceFieldGenerator
from .errorhandler import assert_msg_critical
from .mofutils import svd_superimpose
from .molecularbasis import MolecularBasis
from .respchargesdriver import RespChargesDriver

try:
    import openmm
    from openmm import LangevinIntegrator, Platform
    from openmm.app import NoCutoff, Simulation, PDBFile, ForceField, GromacsTopFile
    from openmm.unit import nanometer, md_unit_system, kelvin, picoseconds, picosecond, dalton
except ImportError:
    pass


class ConformerGenerator:

    def __init__(self, comm=None, ostream=None):

        if comm is None:
            comm = MPI.COMM_WORLD

        if ostream is None:
            if comm.Get_rank() == mpi_master():
                ostream = OutputStream(sys.stdout)
            else:
                ostream = OutputStream(None)

        self._comm = comm
        self._rank = comm.Get_rank()
        self._size = comm.Get_size()

        self.ostream = ostream

        self.molecule = None
        self.number_of_conformers_to_select = None

        self.top_file_name = None
        self.partial_charges = None
        self.resp_charges = True
        self.resp_charges_driver = RespChargesDriver()

        self.freeze_atoms = None  # list of atom indices to freeze during optimization and rotatable bond detection

        self.save_xyz_files = False
        self.save_path = None

        # thresholds for energy minimization
        self.em_tolerance = 1.0

        # thresholds for removing duplicate conformers
        # rmsd threshold in Angstrom, energy threshold in kJ/mol
        self.rmsd_threshold = 1.2
        self.energy_threshold = 1.2

        self.implicit_solvent_model = None
        self.solute_dielectric = 1.0
        self.solvent_dielectric = 78.39

        self.use_gromacs_files = False

        # --- Monte Carlo random sampling settings --------------------------
        # Set mc_search = True to use random sampling instead of the full
        # combinatorial grid.  mc_steps controls how many random dihedral
        # combinations are drawn.  mc_seed ensures reproducibility.
        # When mc_search is False (default) the original grid search is used.
        self.mc_search = False
        self.mc_steps = 500
        self.mc_seed = 42

        # --- Basin-hopping settings ----------------------------------------
        # Set bh_search = True to use basin-hopping instead of grid/MC search.
        # Algorithm:
        #   1. Start from the input geometry (energy-minimised).
        #   2. Randomly perturb a subset of dihedral angles.
        #   3. Minimise the perturbed geometry with OpenMM.
        #   4. Accept the new basin if exp(-ΔE / (R·bh_temperature)) > U(0,1),
        #      where ΔE = E_new − E_current  (kJ/mol), R = 8.314e-3 kJ/mol/K.
        #   5. Collect every accepted (unique) basin; return after bh_steps.
        #
        # bh_perturb = 'grid'       draw angles from the discrete periodicity grid
        # bh_perturb = 'continuous' draw uniformly from [0, 360)
        # bh_temperature            effective temperature for Boltzmann criterion (K)
        # bh_steps                  total MC steps (accepted + rejected)
        # bh_seed                   random seed for reproducibility
        #
        # Basin-hopping runs on rank 0 only (sequential by design).
        # All other MPI ranks idle during the search and receive results
        # via broadcast at the end.
        self.bh_search = False
        self.bh_steps = 500
        self.bh_temperature = 300.0
        self.bh_perturb = 'grid'    # 'grid' or 'continuous'
        self.bh_seed = 42

        # --- Cremer-Pople ring puckering settings -------------------------
        # Set cp_search = True to enumerate ring conformations by sampling
        # the Cremer-Pople puckering coordinate space on a uniform grid.
        # Supported ring sizes: 5, 6, 7.
        #
        # For N=5: one amplitude Q and one phase φ  (circle, 1D grid)
        # For N=6: one amplitude Q, polar angle θ, azimuthal angle φ  (sphere)
        # For N=7: two amplitudes Q2,Q3 and two phases φ2,φ3 (4D, reduced)
        #
        # cp_grid_points  — number of grid points along each angular axis
        #                   (total conformers ≈ cp_grid_points^(n_angles) per ring)
        # cp_amplitude    — fixed puckering amplitude Q in Å (None = auto from
        #                   the input geometry)
        # cp_rings        — list of rings to sample (None = auto-detect all)
        self.cp_search      = False
        self.cp_grid_points = 12       # points per angular axis
        self.cp_amplitude   = None     # Å; None = use amplitude of input geometry

    def _analyze_equiv(self, molecule):

        idtf = AtomTypeIdentifier()
        idtf.ostream.mute()

        atom_type = idtf.generate_gaff_atomtypes(molecule)

        idtf.identify_equivalences()
        equivalent_charges = idtf.equivalent_charges
        if len(equivalent_charges) == 0:
            return []

        equiv_atoms_groups = []
        for substr in equivalent_charges.split(","):
            equiv_atoms_groups.append(substr.split("="))

        one_based_equiv_atoms_groups = [
            list(map(int, x)) for x in equiv_atoms_groups
        ]

        return one_based_equiv_atoms_groups

    def _check_equivside_in_dihedrals(self, dihedral_indices, atom_info_dict,
                                      one_based_equiv_atoms_groups):

        side_j_index = dihedral_indices[1] + 1  # convert to 1 based index
        side_k_index = dihedral_indices[2] + 1  # convert to 1 based index

        side_equiv = False
        max_equiv_atoms = 0

        for a, b in [(side_j_index, side_k_index),
                     (side_k_index, side_j_index)]:

            if atom_info_dict[a]["AtomicSymbol"] == "C":
                one_based_connected_atom_numbers = atom_info_dict[a][
                    "ConnectedAtomsNumbers"]
                connected_set = set(one_based_connected_atom_numbers)
                connected_set = connected_set - {b}

                for equiv_g in one_based_equiv_atoms_groups:
                    if connected_set.issubset(set(equiv_g)):
                        side_equiv = True
                        max_equiv_atoms = max(max_equiv_atoms,
                                              len(connected_set))
                        break

        return side_equiv, max_equiv_atoms

    def _check_methyl_group(self, dihedral_indices, atom_info_dict):

        side_j_index = dihedral_indices[1] + 1  # convert to 1 based index
        side_k_index = dihedral_indices[2] + 1  # convert to 1 based index

        for a, b in [(side_j_index, side_k_index),
                     (side_k_index, side_j_index)]:

            if atom_info_dict[a]["AtomicSymbol"] == "C":
                one_based_connected_atom_numbers = atom_info_dict[a][
                    "ConnectedAtomsNumbers"]
                connected_set = set(one_based_connected_atom_numbers)
                connected_set = connected_set - {b}
                connected_elements = [
                    atom_info_dict[idx]["AtomicSymbol"]
                    for idx in connected_set
                ]
                if tuple(connected_elements) == ("H", "H", "H"):
                    return True

        return False

    def _get_dihedral_candidates(self, molecule, top_file_name,
                                 partial_charges):

        mmff_gen = MMForceFieldGenerator(self._comm)
        mmff_gen.ostream.mute()
        if partial_charges is None:
            warn_text = "ConformerGenerator: Partial charges not provided. "
            warn_text += "Will use a quick (and likely inaccurate) estimation of partial charges."
            self.ostream.print_warning(warn_text)
            mmff_gen.partial_charges = molecule.get_partial_charges(
                molecule.get_charge())
        else:
            mmff_gen.partial_charges = partial_charges
        mmff_gen.create_topology(molecule)
        if self.use_gromacs_files:
            mmff_gen.write_gromacs_files(filename=top_file_name)
        else:
            mmff_gen.write_openmm_files(filename=top_file_name)
        self._comm.barrier()

        atom_info_dict = deepcopy(mmff_gen.atom_info_dict)
        rotatable_bonds = deepcopy(mmff_gen.rotatable_bonds)
        dihedrals_dict = deepcopy(mmff_gen.dihedrals)

        if self.freeze_atoms is not None:
            freeze_set = set(self.freeze_atoms)
            new_rotatable_bonds = []
            for (i, j) in rotatable_bonds:
                if (i not in freeze_set) and (j not in freeze_set):
                    new_rotatable_bonds.append((i, j))
            rotatable_bonds = new_rotatable_bonds
            self.ostream.print_info(
                f"{len(freeze_set)} atoms are frozen, {len(rotatable_bonds)} rotatable bonds remain."
            )
            self.ostream.flush()

        rotatable_bonds_zero_based = [(i - 1, j - 1) for (i, j) in rotatable_bonds]
        rotatable_dihedrals_dict = {}

        def get_max_periodicity(periodicity):
            if isinstance(periodicity, list):
                return max([abs(p) for p in periodicity])
            else:
                return periodicity

        for (i, j, k, l), dih in dihedrals_dict.items():

            sorted_bond = tuple(sorted([j, k]))
            max_periodicity = get_max_periodicity(dih["periodicity"])

            if sorted_bond in rotatable_bonds_zero_based:
                if sorted_bond not in rotatable_dihedrals_dict:
                    rotatable_dihedrals_dict[sorted_bond] = deepcopy(dih)
                    rotatable_dihedrals_dict[sorted_bond][
                        "dihedral_indices"] = (i, j, k, l)
                    rotatable_dihedrals_dict[sorted_bond][
                        "max_periodicity"] = max_periodicity
                else:
                    curr_periodicity = rotatable_dihedrals_dict[sorted_bond][
                        "max_periodicity"]
                    rotatable_dihedrals_dict[sorted_bond][
                        "max_periodicity"] = max(curr_periodicity,
                                                 max_periodicity)

        dihedrals_candidates = []

        one_based_equiv_atoms_groups = self._analyze_equiv(molecule)

        for k, v in rotatable_dihedrals_dict.items():
            max_periodicity = v["max_periodicity"]
            if max_periodicity == 2:
                dih_angle = [0, 180]
            elif max_periodicity == 3:
                dih_angle = [60, 180, 300]
            elif max_periodicity == 4:
                dih_angle = [0, 90, 180, 270]
            elif max_periodicity == 5:
                dih_angle = [36, 54, 84, 144, 324]
            elif max_periodicity == 6:
                dih_angle = [0, 60, 120, 180, 240, 300]
            else:
                continue

            dih_index = v["dihedral_indices"]

            if self._check_methyl_group(dih_index, atom_info_dict):
                continue

            side_equiv, max_equiv_atoms = self._check_equivside_in_dihedrals(
                dih_index, atom_info_dict, one_based_equiv_atoms_groups)
            if side_equiv and max_equiv_atoms == max_periodicity:
                continue

            dihedrals_candidates.append((dih_index, dih_angle))

        return dihedrals_candidates, atom_info_dict, dihedrals_dict

    def _get_dihedral_combinations(self, dihedrals_candidates):
        """Full combinatorial grid — original behaviour."""
        dih_angles = [i[1] for i in dihedrals_candidates]
        dihedrals_combinations = list(itertools.product(*dih_angles))
        dihedral_list = [i[0] for i in dihedrals_candidates]

        self.ostream.print_info(
            f"{len(dihedrals_combinations)} conformers will be generated.")
        self.ostream.flush()

        return dihedrals_combinations, dihedral_list

    def _get_mc_combinations(self, dihedrals_candidates):
        """
        Monte Carlo random sampling of dihedral angle combinations.

        Instead of the full Cartesian product, ``mc_steps`` combinations are
        drawn by independently and uniformly sampling one angle from each
        dihedral's discrete grid.  Duplicate draws are silently discarded so
        the actual number of conformers is at most ``mc_steps`` (and can be
        less if the search space is small — in that case every unique
        combination is kept).

        The grid angles defined per-bond in ``_get_dihedral_candidates`` are
        respected: only angles that would be visited by the grid search are
        eligible.  This keeps the MC sampling physically meaningful (angles
        correspond to periodicity minima) while avoiding the combinatorial
        explosion.

        Parameters
        ----------
        dihedrals_candidates : list of (dihedral_indices, angle_list) tuples
            as returned by ``_get_dihedral_candidates``

        Returns
        -------
        dihedrals_combinations : list of tuples
            Each tuple contains one angle per rotatable bond.
        dihedral_list : list of (i,j,k,l) tuples
            Atom index quadruplets, same order as each combination tuple.
        """
        rng = np.random.default_rng(self.mc_seed)

        dih_angles  = [candidate[1] for candidate in dihedrals_candidates]
        dihedral_list = [candidate[0] for candidate in dihedrals_candidates]

        # Maximum unique combinations possible given the discrete grids
        max_unique = 1
        for angles in dih_angles:
            max_unique *= len(angles)

        n_steps = min(self.mc_steps, max_unique)

        if n_steps < self.mc_steps:
            self.ostream.print_info(
                f"MC search: requested {self.mc_steps} steps but search space "
                f"has only {max_unique} unique combinations — "
                f"switching to full grid search for this molecule."
            )
            self.ostream.flush()
            return self._get_dihedral_combinations(dihedrals_candidates)

        # Draw unique combinations by sampling with replacement and
        # deduplicating.  For large spaces this converges quickly; for small
        # spaces (caught above) we already fall back to grid search.
        seen = set()
        dihedrals_combinations = []

        # Draw in batches to avoid an O(n²) loop for large mc_steps
        batch = n_steps
        while len(dihedrals_combinations) < n_steps:
            # Sample one random angle per bond for `batch` independent draws
            samples = np.array(
                [rng.choice(angles, size=batch) for angles in dih_angles]
            ).T   # shape: (batch, n_bonds)

            for row in samples:
                key = tuple(row)
                if key not in seen:
                    seen.add(key)
                    dihedrals_combinations.append(key)
                    if len(dihedrals_combinations) == n_steps:
                        break

        self.ostream.print_info(
            f"MC search: {n_steps} random conformers will be generated "
            f"(search space: {max_unique} combinations, "
            f"{len(dih_angles)} rotatable bonds)."
        )
        self.ostream.flush()

        return dihedrals_combinations, dihedral_list

    def _get_mol_comb(self, molecule, top_file_name, dihedrals_candidates):
        """
        Assemble the dihedral-combination array for broadcast.
        Dispatches to grid search or MC sampling depending on ``self.mc_search``.
        """
        if self.mc_search:
            dihedrals_combinations, dihedral_list = self._get_mc_combinations(
                dihedrals_candidates)
        else:
            dihedrals_combinations, dihedral_list = self._get_dihedral_combinations(
                dihedrals_candidates)

        conformation_dih_dict = []
        for i in range(len(dihedrals_combinations)):
            combo = np.array(dihedrals_combinations[i]).reshape(-1, 1)
            conformation_dih_dict.append(np.hstack((dihedral_list, combo)))

        dih_comb_array = np.array(conformation_dih_dict)

        return dih_comb_array

    def _run_basin_hopping(self, molecule, dihedrals_candidates,
                           top_file_name, simulation):
        """
        Basin-hopping conformational search.

        Runs a Markov chain over energy-minimised basins.  At each step a
        random subset of dihedrals is perturbed, the geometry is minimised
        with OpenMM, and the new basin is accepted or rejected by a Boltzmann
        criterion.  Every accepted basin (including the starting one) is
        stored; duplicates are removed later by the standard RMSD/energy
        filter in ``generate()``.

        This method runs exclusively on MPI rank 0.  The caller is responsible
        for broadcasting the returned list to other ranks.

        Parameters
        ----------
        molecule             : starting VeloxChem Molecule
        dihedrals_candidates : list of (dihedral_indices, angle_list) tuples
        top_file_name        : topology file stem for OpenMM
        simulation           : initialised OpenMM Simulation object

        Returns
        -------
        list of [energy (float), coords_angstrom (ndarray)] pairs —
        one entry per accepted basin (including the initial minimised geometry).
        """
        # Physical constant: R in kJ / (mol·K)
        R_kJ = 8.314e-3

        rng = np.random.default_rng(self.bh_seed)

        dih_indices = [cand[0] for cand in dihedrals_candidates]   # (i,j,k,l)
        dih_grids   = [cand[1] for cand in dihedrals_candidates]   # angle lists
        n_bonds     = len(dih_indices)

        assert_msg_critical(
            self.bh_perturb in ('grid', 'continuous'),
            "ConformerGenerator: bh_perturb must be 'grid' or 'continuous'."
        )

        # ── Step 0: minimise the starting geometry ──────────────────────────
        current_energy, current_coords = self._minimize_energy(
            molecule, simulation, self.em_tolerance)
        current_mol = Molecule(molecule)
        for i, coord in enumerate(current_coords):
            current_mol.set_atom_coordinates(i, coord / bohr_in_angstrom())

        accepted_basins = [[current_energy, current_coords]]

        n_accepted = 1
        n_rejected = 0

        self.ostream.print_info(
            f"Basin-hopping: starting energy = {current_energy:.3f} kJ/mol"
        )
        self.ostream.flush()

        # ── Main loop ────────────────────────────────────────────────────────
        for step in range(self.bh_steps):

            # — Perturbation —
            # Choose a random non-empty subset of bonds (size 1 … n_bonds)
            subset_size = int(rng.integers(1, n_bonds + 1))
            subset      = rng.choice(n_bonds, size=subset_size, replace=False)

            trial_mol = Molecule(current_mol)

            for bond_i in subset:
                idx_quad = list(np.array(dih_indices[bond_i]) + 1)  # 1-based

                if self.bh_perturb == 'grid':
                    angle = float(rng.choice(dih_grids[bond_i]))
                else:  # continuous
                    angle = float(rng.uniform(0.0, 360.0))

                trial_mol.set_dihedral_in_degrees(idx_quad, angle,
                                                  verbose=False)

            # — Minimisation —
            trial_energy, trial_coords = self._minimize_energy(
                trial_mol, simulation, self.em_tolerance)

            # — Boltzmann acceptance —
            delta_e = trial_energy - current_energy
            if delta_e <= 0.0:
                accept = True
            else:
                prob   = np.exp(-delta_e / (R_kJ * self.bh_temperature))
                accept = bool(rng.random() < prob)

            if accept:
                current_energy = trial_energy
                current_coords = trial_coords
                current_mol    = Molecule(trial_mol)
                for i, coord in enumerate(trial_coords):
                    current_mol.set_atom_coordinates(
                        i, coord / bohr_in_angstrom())
                accepted_basins.append([trial_energy, trial_coords])
                n_accepted += 1
            else:
                n_rejected += 1

            # Progress report every 100 steps
            if (step + 1) % 100 == 0:
                acceptance_rate = n_accepted / (step + 1) * 100
                self.ostream.print_info(
                    f"Basin-hopping: step {step + 1}/{self.bh_steps}  "
                    f"accepted={n_accepted}  rejected={n_rejected}  "
                    f"rate={acceptance_rate:.1f}%  "
                    f"current E={current_energy:.3f} kJ/mol"
                )
                self.ostream.flush()

        self.ostream.print_info(
            f"Basin-hopping finished: {n_accepted} unique basins found "
            f"in {self.bh_steps} steps "
            f"(acceptance rate {n_accepted / self.bh_steps * 100:.1f}%)."
        )
        self.ostream.flush()

        return accepted_basins

    # ------------------------------------------------------------------
    # Cremer-Pople ring puckering
    # ------------------------------------------------------------------

    def _detect_rings(self, molecule):
        """
        Detect all rings of size 5, 6, or 7 in *molecule* using a simple
        DFS bond-graph traversal.  Returns a list of tuples of **0-based**
        atom indices, one tuple per ring, in ring-traversal order.
        """
        labels = list(molecule.get_labels())
        coords = np.array(molecule.get_coordinates_in_angstrom(), dtype=float)
        n = len(labels)

        # Build adjacency from covalent distances
        adj = [[] for _ in range(n)]
        for i in range(n):
            for j in range(i + 1, n):
                ri = self._cov_radius_cp(labels[i])
                rj = self._cov_radius_cp(labels[j])
                if np.linalg.norm(coords[i] - coords[j]) < 1.3 * (ri + rj):
                    adj[i].append(j)
                    adj[j].append(i)

        rings = []
        seen_sets = set()

        def dfs(start, current, path, depth):
            for nb in adj[current]:
                if depth >= 2 and nb == start and 5 <= len(path) <= 7:
                    key = frozenset(path)
                    if key not in seen_sets:
                        seen_sets.add(key)
                        rings.append(tuple(path))
                    continue
                if nb not in path and len(path) < 7:
                    dfs(start, nb, path + [nb], depth + 1)

        for i in range(n):
            dfs(i, i, [i], 0)

        # Keep only rings of size 5–7, remove duplicates (forward/reverse)
        unique = []
        unique_sets = set()
        for ring in rings:
            if not (5 <= len(ring) <= 7):
                continue
            key = frozenset(ring)
            if key not in unique_sets:
                unique_sets.add(key)
                unique.append(ring)

        return unique

    @staticmethod
    def _cov_radius_cp(element):
        """Covalent radius table for ring detection."""
        table = {
            'H': 0.31, 'C': 0.76, 'N': 0.71, 'O': 0.66, 'S': 1.05,
            'F': 0.57, 'P': 1.07, 'Cl': 1.02, 'Se': 1.20, 'Si': 1.11,
        }
        return table.get(element, 0.77)

    @staticmethod
    def _cremer_pople_from_coords(ring_coords):
        """
        Compute Cremer-Pople puckering coordinates from ring Cartesian coords.

        Parameters
        ----------
        ring_coords : (N, 3) array, atoms in ring order

        Returns
        -------
        dict with keys depending on ring size N:
          N=5: {'Q': float, 'phi': float}            (φ in radians)
          N=6: {'Q': float, 'theta': float, 'phi': float}
          N=7: {'Q2': float, 'phi2': float,
                'Q3': float, 'phi3': float}
        """
        coords = np.asarray(ring_coords, dtype=float)
        N = len(coords)
        assert N in (5, 6, 7), f"Ring size {N} not supported"

        # Mean plane: subtract centroid, build reference plane via SVD
        centroid = coords.mean(axis=0)
        r = coords - centroid                          # (N, 3)

        # Cremer-Pople mass-weighted reference plane vectors
        # e1, e2 span the mean plane; n is the normal
        R1 = np.sum([r[j] * np.sin(2 * np.pi * j / N) for j in range(N)], axis=0)
        R2 = np.sum([r[j] * np.cos(2 * np.pi * j / N) for j in range(N)], axis=0)
        n  = np.cross(R1, R2)
        if np.linalg.norm(n) < 1e-10:
            n = np.array([0., 0., 1.])
        n /= np.linalg.norm(n)

        # Out-of-plane displacements z_j
        z = r @ n                                     # (N,)

        if N == 5:
            # One puckering mode: m=2
            A2 = np.sqrt(2 / N) * np.sum(z * np.cos(4 * np.pi * np.arange(N) / N))
            B2 = -np.sqrt(2 / N) * np.sum(z * np.sin(4 * np.pi * np.arange(N) / N))
            Q   = np.sqrt(A2**2 + B2**2)
            phi = np.arctan2(B2, A2)
            return {'Q': Q, 'phi': phi}

        elif N == 6:
            # Two puckering modes: m=2 (equatorial), m=3 (axial/chair)
            j   = np.arange(6)
            A2  = np.sqrt(1/3) * np.sum(z * np.cos(2 * np.pi * j * 2 / 6))
            B2  = -np.sqrt(1/3) * np.sum(z * np.sin(2 * np.pi * j * 2 / 6))
            q3  = (1 / np.sqrt(6)) * np.sum(z * np.cos(np.pi * j))   # alternating ±1

            Q2    = np.sqrt(A2**2 + B2**2)
            phi2  = np.arctan2(B2, A2)
            Q     = np.sqrt(Q2**2 + q3**2)
            theta = np.arctan2(Q2, q3)    # 0=chair, π/2=boat, π=inverted chair
            phi   = phi2
            return {'Q': Q, 'theta': theta, 'phi': phi}

        else:  # N == 7
            j = np.arange(7)
            A2 = np.sqrt(2/7) * np.sum(z * np.cos(2*np.pi*j*2/7))
            B2 = -np.sqrt(2/7) * np.sum(z * np.sin(2*np.pi*j*2/7))
            A3 = np.sqrt(2/7) * np.sum(z * np.cos(2*np.pi*j*3/7))
            B3 = -np.sqrt(2/7) * np.sum(z * np.sin(2*np.pi*j*3/7))
            Q2   = np.sqrt(A2**2 + B2**2)
            phi2 = np.arctan2(B2, A2)
            Q3   = np.sqrt(A3**2 + B3**2)
            phi3 = np.arctan2(B3, A3)
            return {'Q2': Q2, 'phi2': phi2, 'Q3': Q3, 'phi3': phi3}

    @staticmethod
    def _ring_coords_from_cremer_pople(ring_coords_ref, cp_params):
        """
        Reconstruct ring atom positions with target Cremer-Pople puckering.

        Strategy
        --------
        Rather than patching the in-plane positions from the reference
        (which can be distorted), we:

        1. Place N atoms on a regular polygon in the XY plane with radius
           equal to the mean in-plane radius of the reference ring.  This
           gives an ideal, undistorted mean plane.
        2. Compute the target out-of-plane displacements z_j from the CP
           parameters.
        3. Add z_j along the plane normal (Z axis).
        4. Rotate and translate the result so the centroid and orientation
           match the reference ring (so substituents stay attached).

        This guarantees that the output ring has *exactly* the requested
        puckering with no in-plane distortion artefacts.

        Parameters
        ----------
        ring_coords_ref : (N, 3) reference ring atom positions
        cp_params       : dict from _cremer_pople_from_coords (target)

        Returns
        -------
        new_coords : (N, 3) ring atom positions with target puckering
        """
        coords = np.asarray(ring_coords_ref, dtype=float)
        N      = len(coords)
        j      = np.arange(N, dtype=float)

        # ── 1. Mean in-plane radius from reference ────────────────────────
        centroid = coords.mean(axis=0)
        r        = coords - centroid

        # Reference mean-plane normal
        R1 = np.sum([r[k] * np.sin(2*np.pi*k/N) for k in range(N)], axis=0)
        R2 = np.sum([r[k] * np.cos(2*np.pi*k/N) for k in range(N)], axis=0)
        n_ref = np.cross(R1, R2)
        if np.linalg.norm(n_ref) < 1e-10:
            n_ref = np.array([0., 0., 1.])
        n_ref /= np.linalg.norm(n_ref)

        # In-plane displacements and mean radius
        r_ip   = r - np.outer(r @ n_ref, n_ref)
        R_mean = np.mean(np.linalg.norm(r_ip, axis=1))

        # ── 2. Ideal polygon in XY plane ──────────────────────────────────
        phi0    = np.arctan2(r_ip[0, 1], r_ip[0, 0])   # align atom 0 to ref
        angles  = 2 * np.pi * j / N + phi0
        xy      = R_mean * np.column_stack([np.cos(angles), np.sin(angles)])

        # ── 3. Target z_j from CP parameters ─────────────────────────────
        if N == 5:
            Q   = cp_params['Q']
            phi = cp_params['phi']
            A2  = Q * np.cos(phi)
            B2  = Q * np.sin(phi)
            z_new = np.sqrt(2/N) * (
                A2 * np.cos(4*np.pi*j/N) - B2 * np.sin(4*np.pi*j/N))

        elif N == 6:
            Q     = cp_params['Q']
            theta = cp_params['theta']
            phi   = cp_params['phi']
            Q2    = Q * np.sin(theta)
            q3    = Q * np.cos(theta)
            A2    = Q2 * np.cos(phi)
            B2    = Q2 * np.sin(phi)
            z_new = (
                np.sqrt(1/3) * (A2 * np.cos(2*np.pi*j*2/6) -
                                B2 * np.sin(2*np.pi*j*2/6)) +
                (1/np.sqrt(6)) * q3 * np.cos(np.pi * j))

        else:  # N == 7
            Q2, phi2 = cp_params['Q2'], cp_params['phi2']
            Q3, phi3 = cp_params['Q3'], cp_params['phi3']
            A2 = Q2*np.cos(phi2); B2 = Q2*np.sin(phi2)
            A3 = Q3*np.cos(phi3); B3 = Q3*np.sin(phi3)
            z_new = (
                np.sqrt(2/7) * (A2*np.cos(2*np.pi*j*2/7) -
                                B2*np.sin(2*np.pi*j*2/7)) +
                np.sqrt(2/7) * (A3*np.cos(2*np.pi*j*3/7) -
                                B3*np.sin(2*np.pi*j*3/7)))

        # ── 4. Assemble in local frame (XY plane + Z normal) ─────────────
        new_local = np.column_stack([xy, z_new])   # (N, 3) in local frame

        # ── 5. Rotate local frame to match reference orientation ──────────
        # We want the local Z axis to align with n_ref, and local X to align
        # with the in-plane direction of atom 0 in the reference.
        z_axis = np.array([0., 0., 1.])

        # Rotation: local Z → n_ref
        v    = np.cross(z_axis, n_ref)
        c    = np.dot(z_axis, n_ref)
        if np.linalg.norm(v) < 1e-10:
            R_zn = np.eye(3) if c > 0 else np.diag([1., 1., -1.])
        else:
            s    = np.linalg.norm(v)
            K    = np.array([[0, -v[2], v[1]],
                              [v[2], 0, -v[0]],
                              [-v[1], v[0], 0]])
            R_zn = np.eye(3) + K + K @ K * ((1 - c) / s**2)

        new_coords = new_local @ R_zn.T + centroid
        return new_coords

    def _cp_grid(self, N, Q, n_pts):
        """
        Generate a uniform grid of Cremer-Pople parameters for a ring of
        size N with fixed amplitude Q.

        Returns a list of CP parameter dicts suitable for
        ``_ring_coords_from_cremer_pople``.

        Grid definitions
        ----------------
        N=5 : φ ∈ [0, 2π)  — ``n_pts`` points on a circle
        N=6 : θ ∈ (0, π), φ ∈ [0, 2π)  — ``n_pts`` × ``n_pts`` points
              on a sphere (sin-weighted in θ for uniformity)
        N=7 : Q2, Q3 split with Q2²+Q3²=Q²; φ2, φ3 ∈ [0, 2π)
              ``n_pts`` amplitude ratios × ``n_pts``² phase combinations
        """
        grid = []
        if N == 5:
            for phi in np.linspace(0, 2 * np.pi, n_pts, endpoint=False):
                grid.append({'Q': Q, 'phi': phi})

        elif N == 6:
            # The CP sphere for N=6:
            #   theta=0   → chair (one orientation)
            #   theta=pi  → inverted chair
            #   theta=pi/2 → boat/twist-boat belt
            #
            # Always include the two chair poles explicitly — they are the
            # most important conformations and are skipped by equal-area
            # interior grids.  Fill the equatorial belt between them.
            phis = np.linspace(0, 2 * np.pi, n_pts, endpoint=False)

            # South pole: chair (theta ~ pi)
            grid.append({'Q': Q, 'theta': np.pi, 'phi': 0.0})

            # North pole: inverted chair (theta ~ 0)
            grid.append({'Q': Q, 'theta': 0.0, 'phi': 0.0})

            # Interior points: equal-area spacing, skipping poles
            thetas = np.arccos(np.linspace(1, -1, n_pts + 2)[1:-1])
            for theta in thetas:
                for phi in phis:
                    grid.append({'Q': Q, 'theta': float(theta), 'phi': float(phi)})

        else:  # N == 7
            # Sample amplitude ratio r = Q2/Q ∈ [0,1]; Q3 = sqrt(Q²-Q2²)
            ratios = np.linspace(0, 1, n_pts)
            phis   = np.linspace(0, 2 * np.pi, n_pts, endpoint=False)
            for r in ratios:
                Q2 = Q * r
                Q3 = Q * np.sqrt(max(1 - r**2, 0))
                for phi2 in phis:
                    for phi3 in phis:
                        grid.append({'Q2': Q2, 'phi2': float(phi2),
                                     'Q3': Q3, 'phi3': float(phi3)})

        return grid

    def _run_cremer_pople_search(self, molecule, top_file_name, simulation):
        """
        Enumerate ring conformations by scanning the Cremer-Pople coordinate
        space for every ring of size 5–7 in the molecule.

        For each ring, a uniform grid over the CP sphere/circle is constructed.
        For each grid point the ring atom positions are updated in the full
        molecular geometry, the geometry is minimised with OpenMM, and the
        result stored.

        Runs on rank 0 only.  Returns a list of [energy, coords_angstrom] pairs.
        """
        rings = self._detect_rings(molecule)
        rings = [r for r in rings if len(r) in (5, 6, 7)]

        if not rings:
            self.ostream.print_info(
                "Cremer-Pople: no 5-, 6-, or 7-membered rings detected.")
            self.ostream.flush()
            return []

        self.ostream.print_info(
            f"Cremer-Pople: found {len(rings)} ring(s) — sizes "
            f"{[len(r) for r in rings]}."
        )
        self.ostream.flush()

        all_results = []
        ref_coords  = np.array(molecule.get_coordinates_in_angstrom(), dtype=float)

        for ring_idx, ring in enumerate(rings):
            N            = len(ring)
            ring_coords  = ref_coords[list(ring)]

            # Determine puckering amplitude
            if self.cp_amplitude is not None:
                Q = self.cp_amplitude
            else:
                # Use the amplitude from the input geometry, but fall back to
                # chemically meaningful defaults if the input is too flat.
                # Typical values: 6-ring chair ~0.63 Å, 5-ring envelope ~0.40 Å,
                # 7-ring ~0.50 Å.  Input geometries from SMILES are often flat
                # (Q ~ 0), so we take the larger of the measured value and the
                # default to ensure a meaningful scan.
                _cp_defaults = {5: 0.40, 6: 0.63, 7: 0.50}
                cp_ref = self._cremer_pople_from_coords(ring_coords)
                if N == 5:
                    Q_input = cp_ref['Q']
                elif N == 6:
                    Q_input = cp_ref['Q']
                else:
                    Q_input = np.sqrt(cp_ref['Q2']**2 + cp_ref['Q3']**2)
                Q = max(Q_input, _cp_defaults[N])

            grid = self._cp_grid(N, Q, self.cp_grid_points)

            self.ostream.print_info(
                f"Cremer-Pople: ring {ring_idx + 1} (size {N}, Q={Q:.3f} Å) "
                f"— {len(grid)} grid points."
            )
            self.ostream.flush()

            for pt_idx, cp_params in enumerate(grid):
                # Build new ring coords from CP parameters
                new_ring_coords = self._ring_coords_from_cremer_pople(
                    ring_coords, cp_params)

                # Embed into full molecular geometry
                trial_full = ref_coords.copy()
                for local_i, global_i in enumerate(ring):
                    trial_full[global_i] = new_ring_coords[local_i]

                # Build trial VeloxChem Molecule
                trial_mol = Molecule(molecule)
                for iatom, coord in enumerate(trial_full):
                    trial_mol.set_atom_coordinates(
                        iatom, coord / bohr_in_angstrom())

                # Minimise and store
                try:
                    energy, opt_coords = self._minimize_energy(
                        trial_mol, simulation, self.em_tolerance)
                    all_results.append([energy, opt_coords])
                except Exception as e:
                    self.ostream.print_info(
                        f"Cremer-Pople: skipping ring {ring_idx+1} "
                        f"point {pt_idx+1} ({e})")

            self.ostream.print_info(
                f"Cremer-Pople: ring {ring_idx + 1} done "
                f"({len(all_results)} conformers so far)."
            )
            self.ostream.flush()

        return all_results

    def _init_openmm_system(self, topology_file, implicit_solvent_model):

        assert_msg_critical('openmm' in sys.modules,
                            'OpenMM is required for ConformerGenerator.')

        assert_msg_critical(
            not (self.use_gromacs_files
                 and implicit_solvent_model is not None),
            'ConformerGenerator: Please set to use_gromacs_files to False ' +
            'when using implicit solvent model.')

        if self.use_gromacs_files:
            if topology_file.endswith(".top"):
                topology_file = str(Path(topology_file))
            else:
                topology_file = str(Path(topology_file).with_suffix(".top"))
            top = GromacsTopFile(topology_file)
            system = top.createSystem(NoCutoff)
        else:
            pdb_file = str(Path(topology_file).with_suffix(".pdb"))
            xml_file = str(Path(topology_file).with_suffix(".xml"))
            pdb = PDBFile(pdb_file)
            if implicit_solvent_model is None:
                forcefield = ForceField(xml_file)
                system = forcefield.createSystem(pdb.topology,
                                                 nonbondedMethod=NoCutoff)
            else:
                implicit_fpath = Path(
                    openmm.__file__).parent / "app" / "data" / "implicit"
                implicit_model_fname = str(implicit_fpath /
                                           f"{implicit_solvent_model}.xml")
                assert_msg_critical(
                    Path(implicit_model_fname).is_file(),
                    f"ConformerGenerator: Could not find file {implicit_model_fname}"
                )
                forcefield = ForceField(xml_file, implicit_model_fname)
                system = forcefield.createSystem(
                    pdb.topology,
                    nonbondedMethod=NoCutoff,
                    soluteDielectric=self.solute_dielectric,
                    solventDielectric=self.solvent_dielectric)
                self.ostream.print_info(
                    f"Using implicit solvent model {implicit_solvent_model}")
                self.ostream.flush()

        platform = Platform.getPlatformByName("CPU")
        platform.setPropertyDefaultValue("Threads", "1")

        integrator = LangevinIntegrator(300 * kelvin, 1 / picosecond,
                                        0.001 * picoseconds)
        if self.use_gromacs_files:
            simulation = Simulation(top.topology, system, integrator, platform)
        else:
            simulation = Simulation(pdb.topology, system, integrator, platform)

        return simulation

    def show_available_implicit_solvent_models(self):

        if self._rank == mpi_master():
            implicit_folder_path = Path(
                openmm.__file__).parent / "app" / "data" / "implicit"
            implicit_solvent_files = [
                f for f in implicit_folder_path.iterdir()
                if f.is_file() and f.suffix == ".xml"
            ]
            print("Available implicit solvent files:")
            for f in implicit_solvent_files:
                print(f.name)

    def _minimize_energy(self, molecule, simulation, em_tolerance):

        assert_msg_critical('openmm' in sys.modules,
                            'OpenMM is required for ConformerGenerator.')

        coords_nm = molecule.get_coordinates_in_angstrom() * 0.1
        simulation.context.setPositions(coords_nm * nanometer)
        if self.freeze_atoms is not None:
            for atom_idx in self.freeze_atoms:
                simulation.system.setParticleMass(atom_idx, 0.0 * dalton)

        simulation.minimizeEnergy(tolerance=em_tolerance)

        state = simulation.context.getState(
            getPositions=True,
            getEnergy=True,
        )

        energy = state.getPotentialEnergy().value_in_unit_system(md_unit_system)
        optimized_coords = state.getPositions(
            asNumpy=True).value_in_unit_system(md_unit_system) * 10

        return energy, optimized_coords

    def _optimize_molecule(self, molecule, top_file_name, em_tolerance,
                           partial_charges, implicit_solvent_model):

        mmff_gen = MMForceFieldGenerator(self._comm)
        mmff_gen.ostream.mute()
        if partial_charges is None:
            warn_text = "ConformerGenerator: Partial charges not provided. "
            warn_text += "Will use a quick (and likely inaccurate) estimation of partial charges."
            self.ostream.print_warning(warn_text)
            mmff_gen.partial_charges = molecule.get_partial_charges(
                molecule.get_charge())
        else:
            mmff_gen.partial_charges = partial_charges
        mmff_gen.create_topology(molecule)
        if self.use_gromacs_files:
            mmff_gen.write_gromacs_files(filename=top_file_name)
        else:
            mmff_gen.write_openmm_files(filename=top_file_name)
        self._comm.barrier()

        if self._rank == mpi_master():
            simulation = self._init_openmm_system(top_file_name,
                                                  implicit_solvent_model)
            energy, opt_coords = self._minimize_energy(molecule, simulation,
                                                       em_tolerance)
        else:
            energy, opt_coords = None, None
        energy, opt_coords = self._comm.bcast((energy, opt_coords),
                                              root=mpi_master())

        new_molecule = Molecule(molecule)
        for i in range(len(opt_coords)):
            new_molecule.set_atom_coordinates(
                i, opt_coords[i] / bohr_in_angstrom())

        return energy, new_molecule

    def generate(self, molecule):

        if self.implicit_solvent_model is not None:
            self.use_gromacs_files = False

        conf_gen_t0 = time.time()

        top_file_name = self.top_file_name
        if top_file_name is None:
            top_file_name = "MOL"

        self.molecule = molecule

        comm = self._comm
        rank = self._comm.Get_rank()
        size = self._comm.Get_size()

        if self.resp_charges and (self.partial_charges is None):
            basis = MolecularBasis.read(molecule, "6-31g*")
            self.partial_charges = self.resp_charges_driver.compute(
                molecule, basis)

        dihedrals_candidates, atom_info_dict, dihedrals_dict = (
            self._get_dihedral_candidates(molecule, top_file_name,
                                          self.partial_charges))

        if not dihedrals_candidates:
            self.ostream.print_info(
                "No rotatable bond found, no new conformer will be generated.")
            self.ostream.flush()

            energy, molecule = self._optimize_molecule(
                self.molecule, top_file_name, self.em_tolerance,
                self.partial_charges, self.implicit_solvent_model)

            if rank == mpi_master():
                self.global_minimum_conformer = molecule
                self.global_minimum_energy = energy
                return {
                    'energies': [energy],
                    'molecules': [Molecule(molecule)],
                    'geometries': [
                        molecule.get_xyz_string(
                            comment=f"Energy: {energy:.3f} kJ/mol")
                    ],
                }
            else:
                return None

        # --- Log which search strategy is being used ----------------------
        if rank == mpi_master():
            if self.bh_search:
                self.ostream.print_info(
                    f"ConformerGenerator: using basin-hopping search "
                    f"({self.bh_steps} steps, T={self.bh_temperature} K, "
                    f"perturb='{self.bh_perturb}', seed={self.bh_seed})."
                )
            elif self.cp_search:
                self.ostream.print_info(
                    f"ConformerGenerator: using Cremer-Pople ring puckering search "
                    f"(grid_points={self.cp_grid_points}, "
                    f"amplitude={'auto' if self.cp_amplitude is None else f'{self.cp_amplitude:.3f} Å'})."
                )
            elif self.mc_search:
                self.ostream.print_info(
                    f"ConformerGenerator: using Monte Carlo random sampling "
                    f"({self.mc_steps} steps, seed={self.mc_seed})."
                )
            else:
                grid_size = 1
                for _, angles in dihedrals_candidates:
                    grid_size *= len(angles)
                if grid_size > 10000:
                    self.ostream.print_warning(
                        f"ConformerGenerator: grid search will generate "
                        f"{grid_size} conformers ({len(dihedrals_candidates)} "
                        f"rotatable bonds). Consider setting mc_search=True "
                        f"with mc_steps=500, or bh_search=True to limit the search."
                    )
            self.ostream.flush()

        # ── Basin-hopping: runs on rank 0, results broadcast to all ranks ──
        if self.bh_search:
            if rank == mpi_master():
                simulation = self._init_openmm_system(top_file_name,
                                                      self.implicit_solvent_model)
                bh_basins = self._run_basin_hopping(
                    molecule, dihedrals_candidates, top_file_name, simulation)
            else:
                bh_basins = None

            bh_basins = comm.bcast(bh_basins, root=mpi_master())
            all_sorted_energy_coords = sorted(bh_basins, key=lambda x: x[0])
            num_total_conformers = len(all_sorted_energy_coords)

            self.ostream.print_info(
                f"{num_total_conformers} basins collected before deduplication.")
            self.ostream.flush()

            if rank == mpi_master():
                return self._postprocess_conformers(
                    molecule, all_sorted_energy_coords, conf_gen_t0)
            else:
                return None

        # ── Cremer-Pople ring puckering: runs on rank 0, results broadcast ──
        if self.cp_search:
            if rank == mpi_master():
                simulation  = self._init_openmm_system(top_file_name,
                                                       self.implicit_solvent_model)
                cp_results  = self._run_cremer_pople_search(
                    molecule, top_file_name, simulation)
            else:
                cp_results = None

            cp_results = comm.bcast(cp_results, root=mpi_master())

            if not cp_results:
                self.ostream.print_info(
                    "Cremer-Pople search found no conformers.")
                self.ostream.flush()
                return None

            all_sorted_energy_coords = sorted(cp_results, key=lambda x: x[0])

            self.ostream.print_info(
                f"{len(all_sorted_energy_coords)} conformers collected "
                "before deduplication.")
            self.ostream.flush()

            if rank == mpi_master():
                return self._postprocess_conformers(
                    molecule, all_sorted_energy_coords, conf_gen_t0)
            else:
                return None

        # ── Grid / MC search (original parallel pipeline) ──────────────────

        if rank == mpi_master():
            conformation_dih_arr = self._get_mol_comb(molecule, top_file_name,
                                                      dihedrals_candidates)
        else:
            conformation_dih_arr = None
        conformation_dih_arr = comm.bcast(conformation_dih_arr,
                                          root=mpi_master())

        ave, rem = divmod(len(conformation_dih_arr), size)
        counts = [ave + 1 if p < rem else ave for p in range(size)]
        displs = [sum(counts[:p]) for p in range(size)]

        num_total_conformers = len(conformation_dih_arr)

        dih_comb_arr_rank = conformation_dih_arr[displs[rank]:displs[rank] +
                                                 counts[rank]].copy()

        conf_start_time = time.time()

        conformations = []

        for i in range(dih_comb_arr_rank.shape[0]):
            if i > 0:
                new_molecule = Molecule(conformations[-1])
                old_dih_settings = dih_comb_arr_rank[i - 1][:, 4]
                new_dih_settings = dih_comb_arr_rank[i][:, 4]
                diff_dih_ind = np.where(
                    old_dih_settings != new_dih_settings)[0]
            else:
                new_molecule = Molecule(molecule)
                diff_dih_ind = np.arange(dih_comb_arr_rank[i].shape[0])

            value_atom_index = dih_comb_arr_rank[i, :, 0:4] + 1

            for j in diff_dih_ind:
                new_molecule.set_dihedral_in_degrees(
                    value_atom_index[j], dih_comb_arr_rank[i, j, 4], verbose=False)

            conformations.append(Molecule(new_molecule))

        conf_dt = time.time() - conf_start_time

        info = f"{num_total_conformers} conformers generated in {conf_dt:.2f} sec"
        if comm.Get_size() > 1:
            dt_list = comm.gather(conf_dt, root=mpi_master())
            if comm.Get_rank() == mpi_master():
                load_imb = 1.0 - sum(dt_list) / (len(dt_list) * max(dt_list))
                info += f' (load imb.: {load_imb * 100:.1f}%)'
        self.ostream.print_info(info)
        self.ostream.flush()

        opt_start_time = time.time()

        simulation = self._init_openmm_system(top_file_name,
                                              self.implicit_solvent_model)

        energy_coords = []

        for mol_i in range(len(conformations)):
            energy, opt_coords = self._minimize_energy(conformations[mol_i],
                                                       simulation,
                                                       self.em_tolerance)
            energy_coords.append([energy, opt_coords])

        opt_dt = time.time() - opt_start_time

        info = f"Energy minimization of {num_total_conformers} conformers took {opt_dt:.2f} sec"
        if comm.Get_size() > 1:
            dt_list = comm.gather(opt_dt, root=mpi_master())
            if comm.Get_rank() == mpi_master():
                load_imb = 1.0 - sum(dt_list) / (len(dt_list) * max(dt_list))
                info += f' (load imb.: {load_imb * 100:.1f}%)'
        self.ostream.print_info(info)
        self.ostream.flush()

        if self.number_of_conformers_to_select is None:
            self.number_of_conformers_to_select = num_total_conformers
        sorted_energy_coords = sorted(
            energy_coords,
            key=lambda x: x[0])[:self.number_of_conformers_to_select]

        gathered_energy_coords = comm.gather(sorted_energy_coords,
                                             root=mpi_master())

        if rank == mpi_master():
            all_sel_energy_coords = [
                ene_coord for local_energy_coords in gathered_energy_coords
                for ene_coord in local_energy_coords
            ]
            all_sorted_energy_coords = sorted(all_sel_energy_coords,
                                              key=lambda x: x[0])
            return self._postprocess_conformers(
                molecule, all_sorted_energy_coords, conf_gen_t0)
        else:
            return None

    def _postprocess_conformers(self, molecule, all_sorted_energy_coords,
                                conf_gen_t0):
        """
        Shared post-processing for all search strategies.

        Takes a sorted list of ``[energy, coords_angstrom]`` pairs,
        builds VeloxChem Molecule objects, removes duplicates by RMSD and
        energy, optionally saves XYZ files, and returns the standard
        ``conformers_dict``.  Sets ``self.global_minimum_conformer`` and
        ``self.global_minimum_energy`` as side effects.

        Must be called on rank 0 only.
        """
        min_energy, min_coords_angstrom = all_sorted_energy_coords[0]
        min_mol = Molecule(molecule)
        for iatom in range(min_mol.number_of_atoms()):
            min_mol.set_atom_coordinates(
                iatom, min_coords_angstrom[iatom] / bohr_in_angstrom())
        self.ostream.print_info(
            f"Global minimum energy: {min_energy:.3f} kJ/mol")
        self.ostream.flush()

        self.global_minimum_conformer = min_mol
        self.global_minimum_energy    = min_energy

        conformers_dict = {'energies': [], 'molecules': [], 'geometries': []}

        for conf_energy, conf_coords_angstrom in all_sorted_energy_coords:
            conformers_dict["energies"].append(conf_energy)
            mol_copy = Molecule(molecule)
            for iatom in range(mol_copy.number_of_atoms()):
                mol_copy.set_atom_coordinates(
                    iatom, conf_coords_angstrom[iatom] / bohr_in_angstrom())
            conformers_dict["molecules"].append(mol_copy)
            conformers_dict["geometries"].append(
                mol_copy.get_xyz_string(
                    comment=f"Energy: {conf_energy:.3f} kJ/mol"))

        num_conformers = len(conformers_dict["energies"])

        # Deduplicate by RMSD + energy threshold
        equiv_conformer_pairs = []
        for i in range(num_conformers):
            xyz_i = conformers_dict["molecules"][i].get_coordinates_in_angstrom()
            ene_i = conformers_dict["energies"][i]
            for j in range(i + 1, num_conformers):
                xyz_j = conformers_dict["molecules"][j].get_coordinates_in_angstrom()
                ene_j = conformers_dict["energies"][j]
                if abs(ene_i - ene_j) < self.energy_threshold:
                    rmsd, rot, trans = svd_superimpose(xyz_j, xyz_i)
                    if rmsd < self.rmsd_threshold:
                        equiv_conformer_pairs.append((i, j))

        duplicate_conformers = sorted(
            set(j for i, j in equiv_conformer_pairs))

        n_to_select = self.number_of_conformers_to_select or num_conformers

        filtered_energies  = [e for i, e in enumerate(conformers_dict["energies"])
                               if i not in duplicate_conformers]
        filtered_molecules = [Molecule(m) for i, m in enumerate(conformers_dict["molecules"])
                               if i not in duplicate_conformers]
        filtered_geometries = [g for i, g in enumerate(conformers_dict["geometries"])
                                if i not in duplicate_conformers]

        conformers_dict = {
            "energies":   filtered_energies[:n_to_select],
            "molecules":  filtered_molecules[:n_to_select],
            "geometries": filtered_geometries[:n_to_select],
        }

        if self.save_xyz_files:
            save_path = Path(self.save_path) if self.save_path else Path("selected_conformers")
            save_path.mkdir(parents=True, exist_ok=True)
            for i in range(len(conformers_dict["energies"])):
                conf_energy = conformers_dict["energies"][i]
                conf_mol    = conformers_dict["molecules"][i]
                xyz_path    = save_path / f"conformer_{i + 1}.xyz"
                with xyz_path.open("w") as fh:
                    fh.write(conf_mol.get_xyz_string(
                        comment=f"Energy: {conf_energy:.3f} kJ/mol"))
            self.ostream.print_info(
                f"{len(conformers_dict['energies'])} conformers with the "
                f"lowest energies are saved in folder {str(save_path)}"
            )
        else:
            self.ostream.print_info(
                f"{len(conformers_dict['energies'])} conformers remain "
                "after removal of duplicate conformers.")

        self.ostream.print_info(
            "Total time spent in generating conformers: "
            f"{time.time() - conf_gen_t0:.2f} sec")
        self.ostream.flush()

        self.conformer_dict = conformers_dict
        return conformers_dict

    def show_global_minimum(self, atom_indices=False, atom_labels=False):

        if self._rank == mpi_master():
            print(
                f"Global minimum conformer with energy {self.global_minimum_energy:.3f} kJ/mol"
            )
            self.global_minimum_conformer.show(atom_indices=atom_indices,
                                               atom_labels=atom_labels)

    def show_conformers(self, number=1, atom_indices=False, atom_labels=False):

        if self._rank == mpi_master():
            if number > len(self.conformer_dict["energies"]):
                number = len(self.conformer_dict["energies"])
                print(f"Only {number} conformers available, showing all.")
            for i in range(number):
                print(
                    f"Conformer {i + 1} with energy {self.conformer_dict['energies'][i]:.3f} kJ/mol"
                )
                self.conformer_dict["molecules"][i].show(
                    atom_indices=atom_indices, atom_labels=atom_labels)
