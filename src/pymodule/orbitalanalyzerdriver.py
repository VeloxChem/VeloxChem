
"""Shared orbital analysis payload for NBO, VB, and related drivers."""

from dataclasses import dataclass, field
import importlib
import numpy as np

from .oneeints import compute_overlap_integrals
from .veloxchemlib import chemical_element_identifier


@dataclass(frozen=True)
class OrbitalAnalyzerOptions:
    """Options for :class:`OrbitalAnalyzer`."""

    threshold: float = 1.0e-12
    include_mo_analysis: bool = True
    include_nbo_candidates: bool = True
    mo_analysis_top: int = 6
    mo_analysis_threshold: float = 1.0e-2
    lone_pair_min_occupation: float = 1.50
    rydberg_max_occupation: float = 0.50
    bond_min_occupation: float = 1.20
    bond_min_atom_weight: float = 0.10
    pi_min_occupation: float = 0.20
    conjugated_pi_max_path: int = 2
    include_metal_ligand_candidates: bool = True

@dataclass(frozen=True)
class OrbitalAnalysisResult:
    """Structured orbital-analysis result shared by NBO and VB drivers."""

    ao_to_atom: np.ndarray = field(default=None)
    ao_to_l: np.ndarray = field(default=None)
    overlap: np.ndarray = field(default=None)
    density: np.ndarray = field(default=None)
    nao_data: object = None
    spin_data: object = None
    mo_analysis: dict = field(default_factory=dict)
    orbital_candidates: list = field(default_factory=list)
    equivalent_atom_groups: tuple = field(default_factory=tuple)
    local_frames: tuple = field(default_factory=tuple)


# Backward-compatible name used by early OrbitalAnalyzer callers.
OrbitalAnalyzerDiagnostics = OrbitalAnalysisResult


class OrbitalAnalyzer:
    """Centralized orbital analysis and classification for VeloxChem."""

    def __init__(self, molecule, basis, mol_orbs=None, options=None, constraints=None):
        self.molecule = molecule
        self.basis = basis
        self.mol_orbs = mol_orbs
        self.options = options or OrbitalAnalyzerOptions()
        self.constraints = constraints
        self.results = None

    def run(self, mol_orbs=None, constraints=None):
        """Run the shared orbital analysis.

        Without molecular orbitals this returns only AO maps. With molecular
        orbitals it also returns NAO/NPA data, optional MO analysis, spin data,
        and shared NBO-like orbital candidates.
        """

        mol_orbs = self.mol_orbs if mol_orbs is None else mol_orbs
        constraints = self.constraints if constraints is None else constraints
        ao_to_atom, ao_to_l = self.ao_shell_map(self.molecule, self.basis)

        if mol_orbs is None:
            self.results = OrbitalAnalysisResult(
                ao_to_atom=ao_to_atom,
                ao_to_l=ao_to_l,
            )
            return self.results

        nbo = self._nbo_helpers()
        overlap = np.array(compute_overlap_integrals(self.molecule, self.basis))
        density = nbo._build_density_matrix(mol_orbs)
        equivalent_groups = nbo._find_equivalent_atom_groups(self.molecule)
        local_frames = nbo._build_local_atom_frames(self.molecule)
        nao_data = nbo._build_orthonormal_nao_data(
            overlap,
            density,
            ao_to_atom,
            ao_to_l,
            self.molecule.number_of_atoms(),
            equivalent_groups,
            local_frames,
        )
        spin_data = nbo._build_spin_nao_data(mol_orbs, overlap, nao_data)

        mo_analysis = {}
        if self.options.include_mo_analysis:
            mo_analysis = nbo._build_mo_nao_analysis(
                self.molecule,
                mol_orbs,
                overlap,
                nao_data,
                top_n=self.options.mo_analysis_top,
                threshold=self.options.mo_analysis_threshold,
            )

        orbital_candidates = []
        if self.options.include_nbo_candidates:
            orbital_candidates = _build_orbital_candidates(
                self.molecule,
                nao_data,
                spin_data=spin_data,
                lone_pair_min_occupation=self.options.lone_pair_min_occupation,
                rydberg_max_occupation=self.options.rydberg_max_occupation,
                bond_min_occupation=self.options.bond_min_occupation,
                bond_min_atom_weight=self.options.bond_min_atom_weight,
                pi_min_occupation=self.options.pi_min_occupation,
                conjugated_pi_max_path=self.options.conjugated_pi_max_path,
                include_metal_ligand_candidates=(
                    self.options.include_metal_ligand_candidates),
                constraints=constraints,
            )

        self.results = OrbitalAnalysisResult(
            ao_to_atom=ao_to_atom,
            ao_to_l=ao_to_l,
            overlap=overlap,
            density=density,
            nao_data=nao_data,
            spin_data=spin_data,
            mo_analysis=mo_analysis,
            orbital_candidates=orbital_candidates,
            equivalent_atom_groups=equivalent_groups,
            local_frames=local_frames,
        )
        return self.results

    analyze = run

    def classify_nao_data(self, nao_data, spin_data=None, constraints=None):
        """Classify an existing NAO payload into shared orbital candidates."""

        constraints = self.constraints if constraints is None else constraints
        orbital_candidates = _build_orbital_candidates(
            self.molecule,
            nao_data,
            spin_data=spin_data,
            lone_pair_min_occupation=self.options.lone_pair_min_occupation,
            rydberg_max_occupation=self.options.rydberg_max_occupation,
            bond_min_occupation=self.options.bond_min_occupation,
            bond_min_atom_weight=self.options.bond_min_atom_weight,
            pi_min_occupation=self.options.pi_min_occupation,
            conjugated_pi_max_path=self.options.conjugated_pi_max_path,
            include_metal_ligand_candidates=(
                self.options.include_metal_ligand_candidates),
            constraints=constraints,
        )
        ao_to_atom, ao_to_l = self.ao_shell_map(self.molecule, self.basis)
        self.results = OrbitalAnalysisResult(
            ao_to_atom=ao_to_atom,
            ao_to_l=ao_to_l,
            nao_data=nao_data,
            spin_data=spin_data,
            orbital_candidates=orbital_candidates,
        )
        return self.results

    @staticmethod
    def ao_shell_map(molecule, basis):
        l_map = {'s': 0, 'p': 1, 'd': 2, 'f': 3, 'g': 4, 'h': 5}
        ao_to_atom = []
        ao_to_l = []
        for entry in basis.get_ao_basis_map(molecule):
            tokens = entry.split()
            atom_index = int(tokens[0]) - 1
            shell_label = tokens[2]
            shell_type = next(char for char in shell_label if char.isalpha())
            ao_to_atom.append(atom_index)
            ao_to_l.append(l_map.get(shell_type, 0))
        return np.array(ao_to_atom, dtype=int), np.array(ao_to_l, dtype=int)

    @staticmethod
    def _nbo_helpers():
        return importlib.import_module(f'{__package__}.nbodriver')


def _build_orbital_candidates(molecule,
                              nao_data,
                              spin_data=None,
                              lone_pair_min_occupation=1.50,
                              rydberg_max_occupation=0.50,
                              bond_min_occupation=1.20,
                              bond_min_atom_weight=0.10,
                              pi_min_occupation=0.20,
                              conjugated_pi_max_path=2,
                              include_metal_ligand_candidates=True,
                              constraints=None):
    """Build shared orbital candidates from the orthonormal NAO density.

    This is an analyzer-owned candidate layer, not a Lewis assignment layer.
    Candidates are deterministic one-center core/lone-pair orbitals and
    two-center natural orbitals obtained by diagonalizing chemically relevant
    density blocks. NBO and VB drivers should consume these records through
    :class:`OrbitalAnalyzer` rather than depending on a separate classifier
    module.
    """

    nbo = OrbitalAnalyzer._nbo_helpers()

    labels = molecule.get_labels()
    natoms = molecule.number_of_atoms()
    connectivity = molecule.get_connectivity_matrix()
    populations = np.array(nao_data.populations, dtype=float)
    spin_populations = (np.array(spin_data.spin_populations, dtype=float)
                        if spin_data is not None else
                        np.zeros_like(populations))
    density = np.array(nao_data.density, dtype=float)
    atom_map = np.array(nao_data.atom_map, dtype=int)
    angular_map = np.array(nao_data.angular_momentum_map, dtype=int)
    norb = density.shape[0]

    core_by_atom = {}
    one_center_lp = set()
    candidates = []
    serial = 1
    requested_pi_pairs = nbo._requested_pi_pairs(constraints, natoms)
    requested_pi_atoms = {
        int(atom)
        for pair in requested_pi_pairs
        for atom in pair
    }
    closed_shell_anion = (
        int(round(float(molecule.get_multiplicity()))) == 1 and
        float(molecule.get_charge()) < -1.0e-12
    )
    resonance_lone_pair_atoms = set(requested_pi_atoms)
    if closed_shell_anion:
        for atom_i in range(natoms):
            if not nbo._pi_capable_atom(atom_i, atom_map, angular_map):
                continue
            for atom_j in range(atom_i + 1, natoms):
                if not nbo._pi_capable_atom(atom_j, atom_map, angular_map):
                    continue
                distance = nbo._graph_distance(connectivity,
                                               atom_i,
                                               atom_j,
                                               max_depth=conjugated_pi_max_path)
                if distance is None:
                    continue
                resonance_lone_pair_atoms.update((int(atom_i), int(atom_j)))
    resonance_lone_pair_atoms.update(
        nbo._polar_resonance_lone_pair_atoms(molecule, atom_map, angular_map))

    for atom in range(natoms):
        indices = np.where(atom_map == atom)[0]
        ordered = indices[np.argsort(populations[indices])[::-1]]
        nuclear_charge = chemical_element_identifier(labels[atom])
        n_core = min(nbo._n_core_orbitals_from_z(nuclear_charge), len(ordered))
        core_indices = nbo._select_core_indices(ordered,
                                                populations,
                                                angular_map,
                                                n_core)
        core_by_atom[atom] = set(core_indices)

        for nao_index in sorted(core_indices,
                                key=lambda idx: populations[idx],
                                reverse=True):
            vector = nbo._normal_orbital_candidate_vector(norb, {nao_index: 1.0})
            candidates.append({
                'index': serial,
                'type': 'CR',
                'subtype': 'core',
                'atoms': (int(atom),),
                'occupation': nbo._candidate_occupation(vector, density),
                'polarization': {int(atom + 1): 1.0},
                'coefficients': nbo._candidate_coefficients(vector),
                'source': 'one-center NAO occupation',
            })
            serial += 1

        if nuclear_charge <= 2:
            continue

        valence = [idx for idx in ordered if idx not in core_indices]
        spin_active_one_center = set()
        atom_one_electron_added = False
        for nao_index in valence:
            if spin_populations[nao_index] < nbo._SPIN_ACTIVE_MIN_OCCUPATION:
                continue
            spin_active_one_center.add(int(nao_index))
            atom_one_electron_added = True
            vector = nbo._normal_orbital_candidate_vector(norb, {nao_index: 1.0})
            candidates.append({
                'index': serial,
                'type': 'SOMO',
                'subtype': 'one-electron',
                'atoms': (int(atom),),
                'electron_count': 1.0,
                'occupation': nbo._candidate_occupation(vector, density),
                'spin_occupation': float(vector.T @ spin_data.spin_density @ vector),
                'polarization': {int(atom + 1): 1.0},
                'coefficients': nbo._candidate_coefficients(vector),
                'source': 'one-center spin population',
            })
            serial += 1
            candidates.append({
                'index': serial,
                'type': 'LP',
                'subtype': 'radical-lone-pair',
                'atoms': (int(atom),),
                'electron_count': 1.0,
                'occupation': nbo._candidate_occupation(vector, density),
                'spin_occupation': float(vector.T @ spin_data.spin_density @ vector),
                'polarization': {int(atom + 1): 1.0},
                'coefficients': nbo._candidate_coefficients(vector),
                'source': 'one-center spin-active lone-pair candidate',
            })
            serial += 1

        if (spin_data is not None and atom in requested_pi_atoms and
                not atom_one_electron_added):
            p_valence = [
                int(idx) for idx in valence
                if angular_map[idx] == 1 and int(idx) not in spin_active_one_center
            ]
            if p_valence:
                nao_index = max(
                    p_valence,
                    key=lambda idx: abs(float(spin_populations[idx])),
                )
                spin_active_one_center.add(int(nao_index))
                vector = nbo._normal_orbital_candidate_vector(norb, {nao_index: 1.0})
                candidates.append({
                    'index': serial,
                    'type': 'SOMO',
                    'subtype': 'resonance-one-electron',
                    'atoms': (int(atom),),
                    'electron_count': 1.0,
                    'occupation': nbo._candidate_occupation(vector, density),
                    'spin_occupation': float(vector.T @ spin_data.spin_density @ vector),
                    'polarization': {int(atom + 1): 1.0},
                    'coefficients': nbo._candidate_coefficients(vector),
                    'source': 'requested pi-resonance one-electron candidate',
                })
                serial += 1

        lone_pair_indices = set()
        for nao_index in valence:
            if populations[nao_index] < lone_pair_min_occupation:
                continue
            one_center_lp.add(int(nao_index))
            lone_pair_indices.add(int(nao_index))
            vector = nbo._normal_orbital_candidate_vector(norb, {nao_index: 1.0})
            candidates.append({
                'index': serial,
                'type': 'LP',
                'subtype': 'lone-pair',
                'atoms': (int(atom),),
                'occupation': nbo._candidate_occupation(vector, density),
                'polarization': {int(atom + 1): 1.0},
                'coefficients': nbo._candidate_coefficients(vector),
                'source': 'one-center NAO occupation',
            })
            serial += 1

        if atom in resonance_lone_pair_atoms and len(lone_pair_indices) < 3:
            p_valence = [
                int(idx) for idx in valence
                if angular_map[idx] == 1 and
                int(idx) not in spin_active_one_center and
                int(idx) not in lone_pair_indices
            ]
            if p_valence:
                nao_index = max(p_valence, key=lambda idx: populations[idx])
            else:
                nao_index = None
            if nao_index is not None and int(nao_index) not in lone_pair_indices:
                lone_pair_indices.add(int(nao_index))
                vector = nbo._normal_orbital_candidate_vector(norb, {nao_index: 1.0})
                candidates.append({
                    'index': serial,
                    'type': 'LP',
                    'subtype': 'resonance-lone-pair',
                    'atoms': (int(atom),),
                    'occupation': nbo._candidate_occupation(vector, density),
                    'polarization': {int(atom + 1): 1.0},
                    'coefficients': nbo._candidate_coefficients(vector),
                    'source': 'pi-resonance lone-pair candidate',
                })
                serial += 1

        for nao_index in valence:
            if int(nao_index) in lone_pair_indices:
                continue
            if int(nao_index) in spin_active_one_center:
                continue
            if populations[nao_index] > rydberg_max_occupation:
                continue
            vector = nbo._normal_orbital_candidate_vector(norb, {nao_index: 1.0})
            candidates.append({
                'index': serial,
                'type': 'RY',
                'subtype': 'one-center',
                'atoms': (int(atom),),
                'electron_count': 0.0,
                'occupation': nbo._candidate_occupation(vector, density),
                'polarization': {int(atom + 1): 1.0},
                'coefficients': nbo._candidate_coefficients(vector),
                'source': 'one-center low-occupation NAO complement',
            })
            serial += 1

    for atom_i in range(natoms):
        for atom_j in range(atom_i + 1, natoms):
            if connectivity[atom_i, atom_j] == 0:
                continue

            left = [idx for idx in np.where(atom_map == atom_i)[0]
                    if idx not in core_by_atom.get(atom_i, set()) and
                    idx not in one_center_lp and populations[idx] > 1.0e-3]
            right = [idx for idx in np.where(atom_map == atom_j)[0]
                     if idx not in core_by_atom.get(atom_j, set()) and
                     idx not in one_center_lp and populations[idx] > 1.0e-3]
            if not left or not right:
                continue

            pair_indices = np.array(left + right, dtype=int)
            pair_density = density[np.ix_(pair_indices, pair_indices)]
            pair_density = 0.5 * (pair_density + pair_density.T)
            occupations, rotations = np.linalg.eigh(pair_density)
            order = np.argsort(occupations)[::-1]
            max_pair_candidates = min(len(left), len(right), 3)
            accepted = 0
            accepted_candidates = []
            occupied_roots = []
            used_complement_roots = set()

            for root in order:
                if accepted >= max_pair_candidates:
                    break

                occupation = float(occupations[root])
                if occupation < bond_min_occupation:
                    break

                local_vector = rotations[:, root]
                left_weight = float(np.sum(local_vector[:len(left)]**2))
                right_weight = float(np.sum(local_vector[len(left):]**2))
                if (left_weight < bond_min_atom_weight or
                        right_weight < bond_min_atom_weight):
                    continue

                vector = nbo._candidate_vector_from_local(pair_indices,
                                                          local_vector,
                                                          norb)
                subtype = 'sigma' if accepted == 0 else 'pi'
                candidate = {
                    'index': serial,
                    'type': 'BD',
                    'subtype': subtype,
                    'atoms': (int(atom_i), int(atom_j)),
                    'non_sigma': False,
                    'occupation': nbo._candidate_occupation(vector, density),
                    'polarization': {
                        int(atom_i + 1): left_weight,
                        int(atom_j + 1): right_weight,
                    },
                    'coefficients': nbo._candidate_coefficients(vector),
                    'source': 'connected atom-pair density diagonalization',
                }
                candidates.append(candidate)
                accepted_candidates.append(candidate)
                occupied_roots.append(int(root))
                serial += 1
                accepted += 1

            used_complement_roots.update(occupied_roots)
            for candidate in accepted_candidates:
                serial = nbo._append_antibonding_complement(
                    candidates,
                    serial,
                    candidate,
                    occupations,
                    rotations,
                    pair_indices,
                    len(left),
                    density,
                    atom_i,
                    atom_j,
                    bond_min_atom_weight,
                    norb,
                    used_complement_roots,
                    'connected atom-pair antibonding complement',
                )

    existing_pi_pairs = {
        tuple(sorted(candidate.get('atoms', ())))
        for candidate in candidates
        if candidate.get('type') == 'BD' and candidate.get('subtype') == 'pi'
    }

    for atom_i in range(natoms):
        if not nbo._pi_capable_atom(atom_i, atom_map, angular_map):
            continue
        for atom_j in range(atom_i + 1, natoms):
            if not nbo._pi_capable_atom(atom_j, atom_map, angular_map):
                continue

            pair = (atom_i, atom_j)
            distance = nbo._graph_distance(connectivity,
                                           atom_i,
                                           atom_j,
                                           max_depth=conjugated_pi_max_path)
            is_requested = pair in requested_pi_pairs
            if distance is None and not is_requested:
                continue

            allow_requested_pair = (
                not is_requested or
                nbo._allow_requested_pi_pair_candidate(molecule, pair, connectivity)
            )
            if is_requested and not allow_requested_pair:
                continue

            source = ('connected p-p density diagonalization'
                      if connectivity[atom_i, atom_j] != 0 else
                      'conjugated non-sigma p-p density diagonalization')
            non_sigma = bool(connectivity[atom_i, atom_j] == 0)
            force_requested_pair = (
                is_requested and non_sigma and allow_requested_pair
            )
            force_anion_resonance_pair = (
                closed_shell_anion and not non_sigma and distance is not None
            )

            new_candidates, serial = nbo._pi_pair_candidates(
                atom_i,
                atom_j,
                density,
                spin_data.spin_density if spin_data is not None else None,
                atom_map,
                angular_map,
                populations,
                core_by_atom,
                bond_min_atom_weight,
                pi_min_occupation,
                norb,
                serial,
                source,
                non_sigma=non_sigma,
                include_pair_candidate=pair not in existing_pi_pairs,
                force_pair_candidate=(force_requested_pair or
                                      force_anion_resonance_pair),
            )
            candidates.extend(new_candidates)
            if any(candidate.get('electron_count', 2.0) > 1.0
                   for candidate in new_candidates):
                existing_pi_pairs.add(pair)

    if include_metal_ligand_candidates:
        metal_ligand_candidates, serial = _build_metal_ligand_candidates(
            molecule,
            density,
            populations,
            atom_map,
            angular_map,
            core_by_atom,
            norb,
            serial,
        )
        candidates.extend(metal_ligand_candidates)

    return candidates


def _build_metal_ligand_candidates(molecule,
                                   density,
                                   populations,
                                   atom_map,
                                   angular_map,
                                   core_by_atom,
                                   norb,
                                   serial):
    """Build metal-ligand sigma-acceptor and pi-donor diagnostics.

    The records are neutral analyzer candidates. They are not selected Lewis
    bonds and they do not imply a VB active space. The sigma-acceptor channel
    describes ligand-to-metal sigma donation into metal acceptor functions. The
    pi-donor channel describes metal-to-ligand pi back-donation from occupied
    metal d functions into ligand pi-acceptor and nonbonding pi-type functions.
    """

    labels = molecule.get_labels()
    natoms = molecule.number_of_atoms()
    connectivity = molecule.get_connectivity_matrix()
    coords = _safe_coordinates_in_angstrom(molecule)
    candidates = []

    for metal_atom, ligand_atom in _metal_ligand_pairs(
            labels, connectivity, coords):
        metal_indices = np.where(atom_map == metal_atom)[0]
        ligand_indices = np.where(atom_map == ligand_atom)[0]
        if len(metal_indices) == 0 or len(ligand_indices) == 0:
            continue

        metal_valence = [
            int(idx) for idx in metal_indices
            if int(idx) not in core_by_atom.get(metal_atom, set())
        ] or [int(idx) for idx in metal_indices]
        ligand_valence = [
            int(idx) for idx in ligand_indices
            if int(idx) not in core_by_atom.get(ligand_atom, set())
        ] or [int(idx) for idx in ligand_indices]

        # For coordination diagnostics we must not remove the chemically
        # important occupied metal d manifold or ligand nonbonding p/d space
        # just because the generic Lewis candidate builder marked high-
        # population NAOs as core-like.  Back-donation is precisely an
        # occupied metal d -> ligand pi/nonbonding-acceptor diagnostic.
        metal_d_space = [
            int(idx) for idx in metal_indices
            if angular_map[idx] == 2 and populations[idx] > 0.10
        ]
        ligand_pi_space = [
            int(idx) for idx in ligand_indices
            if angular_map[idx] in (1, 2, 3) and populations[idx] < 1.98
        ]

        ligand_sigma_donors = [
            idx for idx in ligand_valence
            if angular_map[idx] <= 1 and populations[idx] > 0.60
        ] or [idx for idx in ligand_valence if populations[idx] > 0.60]
        metal_sigma_acceptors = [
            idx for idx in metal_valence
            if populations[idx] < 1.50
        ] or metal_valence

        if ligand_sigma_donors and metal_sigma_acceptors:
            vector = _metal_ligand_channel_vector(
                ligand_sigma_donors, metal_sigma_acceptors, density, norb)
            strength = _density_coupling_norm(
                density, ligand_sigma_donors, metal_sigma_acceptors)
            candidates.append({
                'index': serial,
                'type': 'ML',
                'subtype': 'sigma-acceptor',
                'atoms': (int(metal_atom), int(ligand_atom)),
                'metal_atom': int(metal_atom),
                'ligand_atom': int(ligand_atom),
                'donor_atom': int(ligand_atom),
                'acceptor_atom': int(metal_atom),
                'donor_acceptor_role': 'ligand sigma donor / metal acceptor',
                'coordination_mode': _coordination_mode(connectivity,
                                                        metal_atom,
                                                        ligand_atom,
                                                        labels),
                'channel': 'ligand-to-metal-sigma-donation',
                'electron_count': 0.0,
                'occupation': float(vector.T @ density @ vector),
                'interaction_strength': strength,
                'donation_strength': strength,
                'back_donation_strength': 0.0,
                'polarization': {
                    int(ligand_atom + 1): 0.5,
                    int(metal_atom + 1): 0.5,
                },
                'coefficients': _candidate_coefficients(vector),
                'source': 'metal-ligand sigma donation diagnostic',
            })
            serial += 1

        metal_pi_donors = [
            idx for idx in metal_d_space
            if populations[idx] > 0.50
        ]
        ligand_pi_acceptors = [
            idx for idx in ligand_pi_space
            if populations[idx] < 1.50
        ]
        if metal_pi_donors and ligand_pi_acceptors:
            vector = _metal_ligand_channel_vector(
                metal_pi_donors, ligand_pi_acceptors, density, norb)
            strength = _density_coupling_norm(
                density, metal_pi_donors, ligand_pi_acceptors)
            candidates.append({
                'index': serial,
                'type': 'ML',
                'subtype': 'pi-donor',
                'atoms': (int(metal_atom), int(ligand_atom)),
                'metal_atom': int(metal_atom),
                'ligand_atom': int(ligand_atom),
                'donor_atom': int(metal_atom),
                'acceptor_atom': int(ligand_atom),
                'donor_acceptor_role': (
                    'metal pi donor / ligand pi nonbonding-or-acceptor'),
                'coordination_mode': _coordination_mode(connectivity,
                                                        metal_atom,
                                                        ligand_atom,
                                                        labels),
                'channel': 'metal-to-ligand-pi-back-donation',
                'electron_count': 0.0,
                'occupation': float(vector.T @ density @ vector),
                'interaction_strength': strength,
                'donation_strength': 0.0,
                'back_donation_strength': strength,
                'polarization': {
                    int(metal_atom + 1): 0.5,
                    int(ligand_atom + 1): 0.5,
                },
                'coefficients': _candidate_coefficients(vector),
                'source': (
                    'metal-ligand pi back-donation diagnostic including '
                    'ligand nonbonding pi-type functions'),
            })
            serial += 1

    return candidates, serial


def _safe_coordinates_in_angstrom(molecule):
    try:
        return np.array(molecule.get_coordinates_in_angstrom(), dtype=float)
    except Exception:
        return None


def _metal_ligand_pairs(labels, connectivity, coords=None):
    pairs = []
    natoms = len(labels)
    for atom_i in range(natoms):
        if not _is_metal_label(labels[atom_i]):
            continue
        nearest_ligand_atom = None
        nearest_distance = float('inf')
        for atom_j in range(natoms):
            if atom_i == atom_j or _is_metal_label(labels[atom_j]):
                continue
            if labels[atom_j].upper() == 'H':
                continue
            connected = bool(connectivity[atom_i, atom_j] != 0)
            close_contact = False
            distance = None
            if coords is not None:
                distance = float(np.linalg.norm(coords[atom_i] - coords[atom_j]))
                close_contact = distance <= _metal_ligand_cutoff(labels[atom_i],
                                                                 labels[atom_j])
            if connected or close_contact:
                pairs.append((int(atom_i), int(atom_j)))
            elif (distance is not None and
                  distance <= _metal_ligand_diagnostic_cutoff(labels[atom_i],
                                                              labels[atom_j]) and
                  distance < nearest_distance):
                nearest_ligand_atom = atom_j
                nearest_distance = distance
        if (nearest_ligand_atom is not None and
                not any(pair[0] == int(atom_i) for pair in pairs)):
            pairs.append((int(atom_i), int(nearest_ligand_atom)))
    return pairs


def _is_metal_label(label):
    nuclear_charge = chemical_element_identifier(label)
    return ((21 <= nuclear_charge <= 30) or
            (39 <= nuclear_charge <= 48) or
            (57 <= nuclear_charge <= 80) or
            (89 <= nuclear_charge <= 112))


def _metal_ligand_cutoff(metal_label, ligand_label):
    ligand_cutoffs = {
        'C': 2.45,
        'N': 2.45,
        'O': 2.35,
        'P': 2.70,
        'S': 2.75,
    }
    return ligand_cutoffs.get(ligand_label.capitalize(), 2.55)


def _metal_ligand_diagnostic_cutoff(metal_label, ligand_label):
    return max(_metal_ligand_cutoff(metal_label, ligand_label), 5.25)


def _coordination_mode(connectivity, metal_atom, ligand_atom, labels):
    metal_neighbors = [
        atom for atom in range(len(labels))
        if atom != ligand_atom and connectivity[ligand_atom, atom] != 0 and
        _is_metal_label(labels[atom])
    ]
    return 'bridging' if len(metal_neighbors) > 1 else 'terminal'


def _density_coupling_norm(density, left_indices, right_indices):
    block = density[np.ix_(left_indices, right_indices)]
    return float(np.linalg.norm(block))


def _metal_ligand_channel_vector(left_indices, right_indices, density, norb):
    block = np.abs(density[np.ix_(left_indices, right_indices)])
    if block.size == 0:
        left_index = int(left_indices[0])
        right_index = int(right_indices[0])
    else:
        local_left, local_right = np.unravel_index(np.argmax(block), block.shape)
        left_index = int(left_indices[local_left])
        right_index = int(right_indices[local_right])
    vector = np.zeros(norb)
    sign = 1.0 if density[left_index, right_index] >= 0.0 else -1.0
    vector[left_index] = 1.0 / np.sqrt(2.0)
    vector[right_index] = sign / np.sqrt(2.0)
    return vector


def _candidate_coefficients(vector, threshold=1.0e-8):
    return {
        int(index): float(value)
        for index, value in enumerate(vector)
        if abs(float(value)) > threshold
    }


# Backward-compatible private name used by early development code.
_build_nbo_candidates = _build_orbital_candidates
