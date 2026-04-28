#
#                                   VELOXCHEM
#              ----------------------------------------------------
#                          An Electronic Structure Code
#
#  SPDX-License-Identifier: BSD-3-Clause
#
#  Copyright 2018-2025 VeloxChem developers

"""
Natural Bond Orbital driver.

This module is a clean restart scaffold. It keeps the public VeloxChem API
importable, implements a small deterministic NAO/NPA baseline, and deliberately
removes the previous threshold-driven Lewis/NBO assignment heuristics.

The replacement implementation should be built as three separate layers:

1. NAO/NPA layer: construct a documented orbital basis and populations.
2. Candidate layer: generate one-centre and two-centre orbital candidates.
3. Assignment layer: select a Lewis structure from candidates using an explicit
   objective and documented constraints.
"""

from dataclasses import asdict, dataclass, field, fields
from itertools import combinations
import sys

from mpi4py import MPI
import numpy as np

from .errorhandler import assert_msg_critical
from .oneeints import compute_overlap_integrals
from .outputstream import OutputStream
from .veloxchemlib import chemical_element_identifier, mpi_master


@dataclass(frozen=True)
class NboComputeOptions:
    """Options for the NBO compute API."""

    include_diagnostics: bool = True
    include_mo_analysis: bool = True
    include_nbo_candidates: bool = True
    include_lewis_assignment: bool = True
    mo_analysis_top: int = 6
    mo_analysis_threshold: float = 1.0e-2
    lone_pair_min_occupation: float = 1.50
    bond_min_occupation: float = 1.20
    bond_min_atom_weight: float = 0.10
    pi_min_occupation: float = 0.20
    conjugated_pi_max_path: int = 2
    max_alternatives: int = 12
    lewis_weight_beta: float = 4.0
    include_nra: bool = False
    nra_subspace: str = 'selected'
    nra_fit_metric: str = 'frobenius'


@dataclass(frozen=True)
class NboConstraints:
    """Constraint container reserved for the future assignment layer."""

    required_bonds: tuple = field(default_factory=tuple)
    forbidden_bonds: tuple = field(default_factory=tuple)
    required_pi_bonds: tuple = field(default_factory=tuple)
    allowed_pi_bonds: tuple = field(default_factory=tuple)
    forbidden_pi_bonds: tuple = field(default_factory=tuple)
    fixed_lone_pairs: tuple = field(default_factory=tuple)
    fixed_core: tuple = field(default_factory=tuple)
    fragment_locks: tuple = field(default_factory=tuple)
    formal_charges: dict = field(default_factory=dict)


@dataclass(frozen=True)
class NaoData:
    """Natural-orbital population data."""

    transform: np.ndarray
    density: np.ndarray
    overlap: np.ndarray
    populations: np.ndarray
    atom_map: np.ndarray
    angular_momentum_map: np.ndarray
    equivalent_atom_groups: tuple = field(default_factory=tuple)
    local_frames: tuple = field(default_factory=tuple)


def _symmetric_orthogonalizer(overlap, threshold=1.0e-10):
    """Return an orthogonalizer X such that X.T @ overlap @ X is identity."""

    eigvals, eigvecs = np.linalg.eigh(overlap)
    keep = eigvals > threshold
    assert_msg_critical(np.any(keep), 'NBO: overlap matrix is rank deficient')
    return eigvecs[:, keep] * eigvals[keep]**-0.5


def _lowdin_inverse_sqrt(overlap, threshold=1.0e-10):
    """Return the square Löwdin inverse square root S^(-1/2)."""

    eigvals, eigvecs = np.linalg.eigh(overlap)
    assert_msg_critical(np.all(eigvals > threshold),
                        'NBO: overlap matrix is rank deficient')
    return (eigvecs * eigvals**-0.5) @ eigvecs.T


def _ao_shell_map(molecule, basis):
    """Map AO indices to atom indices and angular-momentum labels."""

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


def _unit_vector(vector, threshold=1.0e-12):
    """Return a normalized vector or a zero vector if the norm is small."""

    norm = float(np.linalg.norm(vector))
    if norm < threshold:
        return np.zeros(3)
    return vector / norm


def _atom_environment_signature(atom, labels, coords, connectivity, tol=1.0e-3):
    """Build a simple permutation-equivalence signature for an atom."""

    neighbours = np.where(connectivity[atom] != 0)[0]
    pieces = []
    for nbr in neighbours:
        distance = float(np.linalg.norm(coords[atom] - coords[nbr]))
        pieces.append((labels[nbr], int(round(distance / tol))))
    return (labels[atom], tuple(sorted(pieces)))


def _find_equivalent_atom_groups(molecule, tol=1.0e-3):
    """Find atom groups equivalent by element and local connectivity geometry."""

    labels = molecule.get_labels()
    coords = np.array(molecule.get_coordinates_in_bohr())
    connectivity = molecule.get_connectivity_matrix()

    buckets = {}
    for atom in range(molecule.number_of_atoms()):
        signature = _atom_environment_signature(atom,
                                                labels,
                                                coords,
                                                connectivity,
                                                tol=tol)
        buckets.setdefault(signature, []).append(atom)

    groups = [tuple(group) for group in buckets.values() if len(group) > 1]
    return tuple(groups)


def _build_local_atom_frames(molecule):
    """Construct local right-handed coordinate frames for each atom.

    The first axis follows the strongest local bond direction.  The second axis
    is chosen in the local bonding plane when possible, and the third axis is
    the corresponding normal.  Frames are diagnostics and a future hook for NHO
    orientation; the current NAO builder uses them only to report stable local
    geometry information.
    """

    coords = np.array(molecule.get_coordinates_in_bohr())
    connectivity = molecule.get_connectivity_matrix()
    natoms = molecule.number_of_atoms()
    frames = []

    for atom in range(natoms):
        neighbours = np.where(connectivity[atom] != 0)[0]
        if len(neighbours) == 0:
            frames.append(np.eye(3))
            continue

        bond_vectors = [_unit_vector(coords[nbr] - coords[atom])
                        for nbr in neighbours]
        x_axis = bond_vectors[0]
        if len(bond_vectors) > 1:
            trial = bond_vectors[1] - x_axis * float(np.dot(x_axis,
                                                            bond_vectors[1]))
            y_axis = _unit_vector(trial)
        else:
            trial = np.array([1.0, 0.0, 0.0])
            if abs(float(np.dot(trial, x_axis))) > 0.9:
                trial = np.array([0.0, 1.0, 0.0])
            y_axis = _unit_vector(trial - x_axis * float(np.dot(x_axis, trial)))

        z_axis = _unit_vector(np.cross(x_axis, y_axis))
        y_axis = _unit_vector(np.cross(z_axis, x_axis))
        frames.append(np.column_stack([x_axis, y_axis, z_axis]))

    return tuple(frames)


def _symmetrize_density_for_equivalent_atoms(density,
                                             overlap,
                                             ao_to_atom,
                                             ao_to_l,
                                             equivalent_groups):
    """Conservatively symmetrize density by equivalent-atom permutations.

    This first pass only symmetrizes groups whose basis functions are all s-type
    and have identical AO layouts.  That safely enforces equivalent H atoms in
    systems such as H2C=O and CH4 without rotating p/d functions between atoms.

    The symmetrization is an average over AO permutation transforms, not a row
    average.  This is important: row/column averaging is not a valid density
    transformation and can destroy the electron count in a nonorthogonal AO
    basis.  The final electron count is checked and the original density is
    retained if the permutation average is not conservative for the supplied
    overlap matrix.
    """

    original_density = np.array(density, copy=True)
    sym_density = np.array(density, copy=True)
    original_electrons = float(np.trace(overlap @ original_density).real)

    for group in equivalent_groups:
        layouts = []
        indices_by_atom = []
        for atom in group:
            indices = np.where(ao_to_atom == atom)[0]
            indices_by_atom.append(indices)
            layouts.append(tuple(int(ao_to_l[idx]) for idx in indices))

        if len(set(layouts)) != 1:
            continue
        if any(l_val != 0 for l_val in layouts[0]):
            continue

        group_density = np.zeros_like(sym_density)
        ngroup = len(group)
        for shift in range(ngroup):
            permutation = np.arange(sym_density.shape[0])
            for target_pos, target_indices in enumerate(indices_by_atom):
                source_indices = indices_by_atom[(target_pos + shift) % ngroup]
                for target_index, source_index in zip(target_indices,
                                                      source_indices):
                    permutation[target_index] = source_index
            group_density += sym_density[np.ix_(permutation, permutation)]

        sym_density = group_density / float(ngroup)

    sym_density = 0.5 * (sym_density + sym_density.T)
    sym_electrons = float(np.trace(overlap @ sym_density).real)
    if abs(sym_electrons - original_electrons) > max(
            1.0e-8, 1.0e-8 * abs(original_electrons)):
        return original_density

    return sym_density


def _canonicalize_degenerate_rotations(occupations,
                                       rotations,
                                       angular_labels,
                                       threshold=1.0e-5):
    """Make near-degenerate atom-block eigenvectors deterministic."""

    rotations = np.array(rotations, copy=True)
    start = 0
    while start < len(occupations):
        stop = start + 1
        while (stop < len(occupations) and
               abs(float(occupations[stop] - occupations[start])) < threshold):
            stop += 1

        if stop - start > 1:
            cols = list(range(start, stop))
            old_rotations = rotations.copy()
            scores = []
            for col in cols:
                weights = rotations[:, col]**2
                dominant_l = int(angular_labels[int(np.argmax(weights))])
                dominant_ao = int(np.argmax(weights))
                scores.append((dominant_l, dominant_ao, col))
            for new_pos, (_, _, old_col) in enumerate(sorted(scores)):
                rotations[:, start + new_pos] = old_rotations[:, old_col]

        start = stop

    return rotations


def _build_density_matrix(mol_orbs):
    """Build the total AO density matrix from molecular orbitals."""

    from .molecularorbitals import molorb

    def spin_density(occupations, coefficients):
        occ = np.asarray(occupations)
        coeff = np.asarray(coefficients)
        selected = occ > 1.0e-12
        return coeff[:, selected] @ np.diag(occ[selected]) @ coeff[:, selected].T

    density = spin_density(mol_orbs.occa_to_numpy(), mol_orbs.alpha_to_numpy())

    if mol_orbs.get_orbitals_type() == molorb.rest:
        return 2.0 * density

    density += spin_density(mol_orbs.occb_to_numpy(), mol_orbs.beta_to_numpy())
    return density


def _mo_spin_blocks(mol_orbs):
    """Return spin blocks with physical occupations for MO analysis."""

    from .molecularorbitals import molorb

    orb_type = mol_orbs.get_orbitals_type()
    alpha_coeff = mol_orbs.alpha_to_numpy()
    alpha_energies = mol_orbs.ea_to_numpy()
    alpha_occupations = mol_orbs.occa_to_numpy()

    if orb_type == molorb.rest:
        return ((
            'restricted',
            alpha_coeff,
            alpha_energies,
            2.0 * alpha_occupations,
        ),)

    beta_coeff = mol_orbs.beta_to_numpy()
    beta_energies = mol_orbs.eb_to_numpy()
    beta_occupations = mol_orbs.occb_to_numpy()

    return ((
        'alpha',
        alpha_coeff,
        alpha_energies,
        alpha_occupations,
    ), (
        'beta',
        beta_coeff,
        beta_energies,
        beta_occupations,
    ))


def _build_mo_nao_analysis(molecule,
                           mol_orbs,
                           overlap,
                           nao_data,
                           top_n=6,
                           threshold=1.0e-2):
    """Analyze canonical molecular orbitals in the orthonormal NAO basis."""

    labels = molecule.get_labels()
    atom_map = np.array(nao_data.atom_map, dtype=int)
    angular_map = np.array(nao_data.angular_momentum_map, dtype=int)
    l_symbols = {0: 's', 1: 'p', 2: 'd', 3: 'f', 4: 'g', 5: 'h'}
    top_n = max(1, int(top_n))
    threshold = max(0.0, float(threshold))

    spin_blocks = []
    max_norm_error = 0.0
    for spin_label, ao_coeff, energies, occupations in _mo_spin_blocks(mol_orbs):
        nao_coeff = nao_data.transform.T @ overlap @ ao_coeff
        orbitals = []

        for mo_index in range(nao_coeff.shape[1]):
            coeff = nao_coeff[:, mo_index]
            weights = np.abs(coeff)**2
            norm = float(np.sum(weights).real)
            max_norm_error = max(max_norm_error, abs(norm - 1.0))
            if norm > 1.0e-14:
                weights = weights / norm

            atom_weights = np.zeros(molecule.number_of_atoms())
            for nao_index, atom in enumerate(atom_map):
                atom_weights[atom] += float(weights[nao_index].real)

            top_naos = []
            order = np.argsort(weights)[::-1]
            for nao_index in order:
                weight = float(weights[nao_index].real)
                if len(top_naos) >= top_n:
                    break
                if weight < threshold and top_naos:
                    break
                atom = int(atom_map[nao_index])
                l_val = int(angular_map[nao_index])
                top_naos.append({
                    'nao_index': int(nao_index + 1),
                    'atom_index': int(atom + 1),
                    'atom_label': labels[atom],
                    'l': l_symbols.get(l_val, '?'),
                    'coefficient': float(coeff[nao_index].real),
                    'weight': weight,
                })

            top_atoms = []
            for atom in np.argsort(atom_weights)[::-1]:
                weight = float(atom_weights[atom].real)
                if weight < threshold and top_atoms:
                    break
                top_atoms.append({
                    'atom_index': int(atom + 1),
                    'atom_label': labels[atom],
                    'weight': weight,
                })

            orbitals.append({
                'mo_index': int(mo_index + 1),
                'spin': spin_label,
                'energy': float(energies[mo_index]),
                'occupation': float(occupations[mo_index]),
                'normalization': norm,
                'top_naos': top_naos,
                'top_atoms': top_atoms,
            })

        spin_blocks.append({
            'spin': spin_label,
            'coefficients': nao_coeff,
            'energies': np.array(energies, copy=True),
            'occupations': np.array(occupations, copy=True),
            'orbitals': orbitals,
        })

    return {
        'spin_blocks': spin_blocks,
        'top_n': top_n,
        'threshold': threshold,
        'max_normalization_error': float(max_norm_error),
    }


def _lowdin_population_matrix(overlap, density):
    """Return the Löwdin AO population matrix S^(1/2) D S^(1/2)."""

    eigvals, eigvecs = np.linalg.eigh(overlap)
    assert_msg_critical(np.all(eigvals > 1.0e-12),
                        'NBO: AO overlap matrix is not positive definite')
    overlap_half = (eigvecs * np.sqrt(eigvals)) @ eigvecs.T
    return overlap_half @ density @ overlap_half


def _build_orthonormal_nao_data(overlap,
                                density,
                                ao_to_atom,
                                ao_to_l,
                                natoms,
                                equivalent_groups=(),
                                local_frames=()):
    """Build a globally orthonormal NAO basis.

    This restart implementation uses a conservative two-step construction:

    1. Build a square Löwdin-orthogonalized AO basis.  The column ordering is
       still the AO ordering, so each orthogonalized AO keeps its original atom
       and angular-momentum label.
    2. Diagonalize the density matrix within each atom block only.  These
       rotations preserve global orthonormality and prevent the inter-atom
       label stealing seen with full-basis Löwdin or weighted Löwdin mixing.
    """

    density = _symmetrize_density_for_equivalent_atoms(density,
                                                       overlap,
                                                       ao_to_atom,
                                                       ao_to_l,
                                                       equivalent_groups)

    nao = overlap.shape[0]
    lowdin_ao = _lowdin_inverse_sqrt(overlap)
    orth_density = lowdin_ao.T @ overlap @ density @ overlap @ lowdin_ao
    orth_density = 0.5 * (orth_density + orth_density.T)

    block_rotation = np.eye(nao)
    final_atom_map = np.array(ao_to_atom, dtype=int).copy()
    final_angular_map = np.array(ao_to_l, dtype=int).copy()

    for atom in range(natoms):
        atom_indices = np.where(ao_to_atom == atom)[0]
        if atom_indices.size == 0:
            continue

        atom_density = orth_density[np.ix_(atom_indices, atom_indices)]
        atom_density = 0.5 * (atom_density + atom_density.T)
        occupations, rotations = np.linalg.eigh(atom_density)
        order = np.argsort(occupations)[::-1]
        occupations = occupations[order]
        rotations = rotations[:, order]
        rotations = _canonicalize_degenerate_rotations(
            occupations,
            rotations,
            ao_to_l[atom_indices],
        )
        block_rotation[np.ix_(atom_indices, atom_indices)] = rotations

        for local_col, nao_index in enumerate(atom_indices):
            weights = rotations[:, local_col]**2
            l_weights = {}
            for l_val in np.unique(ao_to_l[atom_indices]):
                mask = ao_to_l[atom_indices] == l_val
                l_weights[int(l_val)] = float(np.sum(weights[mask]))
            final_angular_map[nao_index] = max(l_weights, key=l_weights.get)

    transform = lowdin_ao @ block_rotation

    nao_overlap = transform.T @ overlap @ transform
    nao_density = transform.T @ overlap @ density @ overlap @ transform
    nao_density = 0.5 * (nao_density + nao_density.T)
    populations = np.diag(nao_density).real.copy()

    return NaoData(transform=transform,
                   density=nao_density,
                   overlap=nao_overlap,
                   populations=populations,
                   atom_map=final_atom_map,
                   angular_momentum_map=final_angular_map,
                   equivalent_atom_groups=tuple(equivalent_groups),
                   local_frames=tuple(local_frames))


def _n_core_orbitals_from_z(z):
    """Return a simple closed-shell core-orbital count."""

    z = int(round(float(z)))
    if z <= 2:
        return 0
    if z <= 10:
        return 1
    if z <= 18:
        return 5
    if z <= 36:
        return 9
    return 18


def _valence_principal_n_from_z(z):
    """Return a simple valence principal quantum number."""

    z = int(round(float(z)))
    if z <= 2:
        return 1
    if z <= 10:
        return 2
    if z <= 18:
        return 3
    if z <= 36:
        return 4
    if z <= 54:
        return 5
    return 6


def _select_core_indices(indices, populations, angular_map, n_core):
    """Select core NAO indices from atom-local candidates.

    The current scaffold does not track principal quantum numbers explicitly.
    For first-row atoms this function prevents unphysical ``Cor(1p)`` labels by
    assigning the one core orbital from the highest-populated s-type NAO.
    Remaining core slots for heavier atoms are then filled by population.
    """

    if n_core <= 0 or len(indices) == 0:
        return set()

    selected = []

    s_candidates = [idx for idx in indices if angular_map[idx] == 0]
    s_candidates = sorted(s_candidates,
                          key=lambda idx: populations[idx],
                          reverse=True)
    if s_candidates:
        selected.append(s_candidates[0])

    remaining = [idx for idx in indices if idx not in selected]
    remaining = sorted(remaining,
                       key=lambda idx: populations[idx],
                       reverse=True)

    for idx in remaining:
        if len(selected) >= n_core:
            break
        selected.append(idx)

    return set(selected)


def _normal_orbital_candidate_vector(size, components):
    """Return a normalized real coefficient vector from sparse components."""

    vector = np.zeros(size)
    for index, coefficient in components.items():
        vector[int(index)] = float(coefficient)
    norm = float(np.linalg.norm(vector))
    if norm > 1.0e-14:
        vector /= norm
    return vector


def _candidate_occupation(vector, density):
    """Return the density occupation of a normalized NAO-space vector."""

    return float(vector.T @ density @ vector)


def _candidate_coefficients(vector, threshold=1.0e-6):
    """Return compact 1-based NAO coefficients for reporting/storage."""

    return [{
        'nao_index': int(index + 1),
        'coefficient': float(vector[index]),
    } for index in np.where(np.abs(vector) > threshold)[0]]


def _graph_distance(connectivity, start, stop, max_depth=3):
    """Return graph distance up to max_depth, or None if not reached."""

    if start == stop:
        return 0
    visited = {int(start)}
    frontier = {int(start)}
    for depth in range(1, int(max_depth) + 1):
        next_frontier = set()
        for atom in frontier:
            neighbours = np.where(connectivity[atom] != 0)[0]
            for neighbour in neighbours:
                neighbour = int(neighbour)
                if neighbour == stop:
                    return depth
                if neighbour not in visited:
                    visited.add(neighbour)
                    next_frontier.add(neighbour)
        frontier = next_frontier
        if not frontier:
            break
    return None


def _pi_capable_atom(atom, atom_map, angular_map):
    """Return whether an atom has at least one p-type NAO."""

    indices = np.where(atom_map == atom)[0]
    return bool(np.any(angular_map[indices] == 1))


def _requested_pi_pairs(constraints, natoms):
    """Return normalized user-requested pi pairs."""

    pairs = set()
    if constraints is None:
        return pairs
    for attr in ('required_pi_bonds', 'allowed_pi_bonds'):
        for pair in getattr(constraints, attr, ()): 
            pairs.add(_normalize_atom_pair(pair, natoms))
    return pairs


def _pi_pair_candidates(atom_i,
                        atom_j,
                        density,
                        atom_map,
                        angular_map,
                        populations,
                        core_by_atom,
                        bond_min_atom_weight,
                        pi_min_occupation,
                        norb,
                        serial,
                        source):
    """Generate p-space two-center pi candidates for an atom pair."""

    left = [idx for idx in np.where(atom_map == atom_i)[0]
            if idx not in core_by_atom.get(atom_i, set()) and
            angular_map[idx] == 1 and populations[idx] > 1.0e-3]
    right = [idx for idx in np.where(atom_map == atom_j)[0]
             if idx not in core_by_atom.get(atom_j, set()) and
             angular_map[idx] == 1 and populations[idx] > 1.0e-3]
    if not left or not right:
        return [], serial

    pair_indices = np.array(left + right, dtype=int)
    pair_density = density[np.ix_(pair_indices, pair_indices)]
    pair_density = 0.5 * (pair_density + pair_density.T)
    occupations, rotations = np.linalg.eigh(pair_density)
    order = np.argsort(occupations)[::-1]
    candidates = []
    accepted = 0

    for root in order:
        if accepted >= 1:
            break
        occupation = float(occupations[root])
        if occupation < pi_min_occupation:
            break
        local_vector = rotations[:, root]
        left_weight = float(np.sum(local_vector[:len(left)]**2))
        right_weight = float(np.sum(local_vector[len(left):]**2))
        if (left_weight < bond_min_atom_weight or
                right_weight < bond_min_atom_weight):
            continue

        components = {
            int(index): float(coefficient)
            for index, coefficient in zip(pair_indices, local_vector)
        }
        vector = _normal_orbital_candidate_vector(norb, components)
        candidates.append({
            'index': serial,
            'type': 'BD',
            'subtype': 'pi',
            'atoms': (int(atom_i), int(atom_j)),
            'non_sigma': True,
            'occupation': _candidate_occupation(vector, density),
            'polarization': {
                int(atom_i + 1): left_weight,
                int(atom_j + 1): right_weight,
            },
            'coefficients': _candidate_coefficients(vector),
            'source': source,
        })
        serial += 1
        accepted += 1

    return candidates, serial


def _build_nbo_candidates(molecule,
                          nao_data,
                          lone_pair_min_occupation=1.50,
                          bond_min_occupation=1.20,
                          bond_min_atom_weight=0.10,
                          pi_min_occupation=0.20,
                          conjugated_pi_max_path=2,
                          constraints=None):
    """Build first-pass NBO candidates from the orthonormal NAO density.

    This is a candidate layer, not a Lewis assignment layer.  Candidates are
    deterministic one-center core/lone-pair orbitals and two-center natural
    orbitals obtained by diagonalizing connected atom-pair density blocks after
    removing core and obvious lone-pair NAOs from the bond search space.
    """

    labels = molecule.get_labels()
    natoms = molecule.number_of_atoms()
    connectivity = molecule.get_connectivity_matrix()
    populations = np.array(nao_data.populations, dtype=float)
    density = np.array(nao_data.density, dtype=float)
    atom_map = np.array(nao_data.atom_map, dtype=int)
    angular_map = np.array(nao_data.angular_momentum_map, dtype=int)
    norb = density.shape[0]

    core_by_atom = {}
    one_center_lp = set()
    candidates = []
    serial = 1

    for atom in range(natoms):
        indices = np.where(atom_map == atom)[0]
        ordered = indices[np.argsort(populations[indices])[::-1]]
        nuclear_charge = chemical_element_identifier(labels[atom])
        n_core = min(_n_core_orbitals_from_z(nuclear_charge), len(ordered))
        core_indices = _select_core_indices(ordered,
                                            populations,
                                            angular_map,
                                            n_core)
        core_by_atom[atom] = set(core_indices)

        for nao_index in sorted(core_indices,
                                key=lambda idx: populations[idx],
                                reverse=True):
            vector = _normal_orbital_candidate_vector(norb, {nao_index: 1.0})
            candidates.append({
                'index': serial,
                'type': 'CR',
                'subtype': 'core',
                'atoms': (int(atom),),
                'occupation': _candidate_occupation(vector, density),
                'polarization': {int(atom + 1): 1.0},
                'coefficients': _candidate_coefficients(vector),
                'source': 'one-center NAO occupation',
            })
            serial += 1

        if nuclear_charge <= 2:
            continue

        valence = [idx for idx in ordered if idx not in core_indices]
        for nao_index in valence:
            if populations[nao_index] < lone_pair_min_occupation:
                continue
            one_center_lp.add(int(nao_index))
            vector = _normal_orbital_candidate_vector(norb, {nao_index: 1.0})
            candidates.append({
                'index': serial,
                'type': 'LP',
                'subtype': 'lone-pair',
                'atoms': (int(atom),),
                'occupation': _candidate_occupation(vector, density),
                'polarization': {int(atom + 1): 1.0},
                'coefficients': _candidate_coefficients(vector),
                'source': 'one-center NAO occupation',
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

                components = {
                    int(index): float(coefficient)
                    for index, coefficient in zip(pair_indices, local_vector)
                }
                vector = _normal_orbital_candidate_vector(norb, components)
                subtype = 'sigma' if accepted == 0 else 'pi'
                candidates.append({
                    'index': serial,
                    'type': 'BD',
                    'subtype': subtype,
                    'atoms': (int(atom_i), int(atom_j)),
                    'non_sigma': False,
                    'occupation': _candidate_occupation(vector, density),
                    'polarization': {
                        int(atom_i + 1): left_weight,
                        int(atom_j + 1): right_weight,
                    },
                    'coefficients': _candidate_coefficients(vector),
                    'source': 'connected atom-pair density diagonalization',
                })
                serial += 1
                accepted += 1

    requested_pi_pairs = _requested_pi_pairs(constraints, natoms)
    existing_pi_pairs = {
        tuple(sorted(candidate.get('atoms', ())))
        for candidate in candidates
        if candidate.get('type') == 'BD' and candidate.get('subtype') == 'pi'
    }

    for atom_i in range(natoms):
        if not _pi_capable_atom(atom_i, atom_map, angular_map):
            continue
        for atom_j in range(atom_i + 1, natoms):
            if connectivity[atom_i, atom_j] != 0:
                continue
            if not _pi_capable_atom(atom_j, atom_map, angular_map):
                continue

            pair = (atom_i, atom_j)
            distance = _graph_distance(connectivity,
                                       atom_i,
                                       atom_j,
                                       max_depth=conjugated_pi_max_path)
            is_requested = pair in requested_pi_pairs
            if distance is None and not is_requested:
                continue
            if pair in existing_pi_pairs:
                continue

            new_candidates, serial = _pi_pair_candidates(
                atom_i,
                atom_j,
                density,
                atom_map,
                angular_map,
                populations,
                core_by_atom,
                bond_min_atom_weight,
                pi_min_occupation,
                norb,
                serial,
                'conjugated non-sigma p-p density diagonalization',
            )
            candidates.extend(new_candidates)
            if new_candidates:
                existing_pi_pairs.add(pair)

    return candidates


def _candidate_type_counts(candidates):
    """Count NBO candidate types."""

    counts = {}
    for candidate in candidates:
        candidate_type = candidate.get('type', '?')
        counts[candidate_type] = counts.get(candidate_type, 0) + 1
    return counts


def _report_atom_text(candidate, labels):
    """Return atom labels ordered from heavier to lighter for reporting."""

    atoms = tuple(candidate.get('atoms', ()))
    ordered_atoms = sorted(
        atoms,
        key=lambda atom: (-chemical_element_identifier(labels[atom]), atom),
    )
    return '-'.join(f'{labels[atom]}{atom + 1}' for atom in ordered_atoms)


def _nbo_report_sort_key(candidate, labels):
    """Sort NBOs for compact chemical reporting."""

    candidate_type = candidate.get('type', '?')
    subtype = candidate.get('subtype', '')
    atoms = tuple(candidate.get('atoms', ()))
    atomic_numbers = sorted(
        (chemical_element_identifier(labels[atom]) for atom in atoms),
        reverse=True,
    )
    padded_atomic_numbers = tuple([-int(z) for z in atomic_numbers] + [0, 0])[:2]

    type_priority = {
        'CR': 0,
        'BD': 1,
        'LP': 2,
        'RY': 3,
        'BD*': 4,
    }
    subtype_priority = {
        'core': 0,
        'sigma': 0,
        'pi': 1,
        'lone-pair': 0,
    }

    return (
        type_priority.get(candidate_type, 9),
        padded_atomic_numbers,
        subtype_priority.get(subtype, 9),
        tuple(sorted(atoms)),
        -float(candidate.get('occupation', 0.0)),
        int(candidate.get('index', 0)),
    )


def _normalize_atom_pair(pair, natoms):
    """Normalize a user atom pair to a sorted zero-based pair."""

    assert_msg_critical(len(pair) == 2,
                        'NBO constraints: atom pairs must have length two')
    i, j = int(pair[0]), int(pair[1])
    if 1 <= i <= natoms and 1 <= j <= natoms:
        i -= 1
        j -= 1
    assert_msg_critical(0 <= i < natoms and 0 <= j < natoms and i != j,
                        f'NBO constraints: invalid atom pair {pair}')
    return tuple(sorted((i, j)))


def _constraint_pair_sets(constraints, natoms):
    """Return normalized required and forbidden bond-pair sets."""

    required = {
        _normalize_atom_pair(pair, natoms)
        for pair in constraints.required_bonds
    }
    forbidden = {
        _normalize_atom_pair(pair, natoms)
        for pair in constraints.forbidden_bonds
    }
    return required, forbidden


def _constraint_pi_pair_sets(constraints, natoms):
    """Return normalized required, allowed, and forbidden pi-bond sets."""

    required = {
        _normalize_atom_pair(pair, natoms)
        for pair in constraints.required_pi_bonds
    }
    allowed = {
        _normalize_atom_pair(pair, natoms)
        for pair in constraints.allowed_pi_bonds
    }
    forbidden = {
        _normalize_atom_pair(pair, natoms)
        for pair in constraints.forbidden_pi_bonds
    }
    return required, allowed, forbidden


def _candidate_sort_key(candidate):
    """Sort key for deterministic primary candidate assignment."""

    type_priority = {'CR': 0, 'LP': 1, 'BD': 2}
    subtype_priority = {'sigma': 0, 'pi': 1}
    return (
        type_priority.get(candidate.get('type', '?'), 9),
        subtype_priority.get(candidate.get('subtype', ''), 9),
        -float(candidate.get('occupation', 0.0)),
        tuple(candidate.get('atoms', ())),
        int(candidate.get('index', 0)),
    )


def _build_primary_assignment(molecule, candidates, constraints):
    """Select a first Lewis-like primary set from generated candidates.

    This is intentionally a small, deterministic assignment layer.  It selects
    occupied candidates up to the electron-pair count, respecting forbidden
    bonds and prioritizing required bonds.  It does not yet optimize resonance
    alternatives or generate anti-bonding/Rydberg complements.
    """

    natoms = molecule.number_of_atoms()
    target_pairs = int(round(0.5 * molecule.number_of_electrons()))
    required_pairs, forbidden_pairs = _constraint_pair_sets(constraints, natoms)
    required_pi_pairs, allowed_pi_pairs, forbidden_pi_pairs = _constraint_pi_pair_sets(
        constraints,
        natoms,
    )

    def allowed(candidate):
        if candidate.get('type') != 'BD':
            return True
        atoms = tuple(sorted(candidate.get('atoms', ())))
        if atoms in forbidden_pairs:
            return False
        if candidate.get('subtype') == 'pi' and atoms in forbidden_pi_pairs:
            return False
        if (candidate.get('subtype') == 'pi' and
                candidate.get('non_sigma', False) and
                atoms not in required_pi_pairs and atoms not in allowed_pi_pairs):
            return False
        return True

    pool = [candidate for candidate in candidates if allowed(candidate)]
    selected = []
    selected_ids = set()
    warnings = []

    for candidate_type in ('CR',):
        for candidate in sorted((cand for cand in pool
                                 if cand.get('type') == candidate_type),
                                key=_candidate_sort_key):
            selected.append(candidate)
            selected_ids.add(candidate.get('index'))

    for pair in sorted(required_pairs):
        pair_candidates = [
            cand for cand in pool
            if cand.get('type') == 'BD' and
            tuple(sorted(cand.get('atoms', ()))) == pair and
            cand.get('index') not in selected_ids
        ]
        if not pair_candidates:
            warnings.append(
                f'Required bond {pair[0] + 1}-{pair[1] + 1} has no candidate.'
            )
            continue
        for candidate in sorted(pair_candidates, key=_candidate_sort_key):
            if len(selected) >= target_pairs:
                break
            selected.append(candidate)
            selected_ids.add(candidate.get('index'))

    for pair in sorted(required_pi_pairs):
        pair_candidates = [
            cand for cand in pool
            if cand.get('type') == 'BD' and
            cand.get('subtype') == 'pi' and
            tuple(sorted(cand.get('atoms', ()))) == pair and
            cand.get('index') not in selected_ids
        ]
        if not pair_candidates:
            warnings.append(
                f'Required pi bond {pair[0] + 1}-{pair[1] + 1} has no candidate.'
            )
            continue
        for candidate in sorted(pair_candidates, key=_candidate_sort_key):
            if len(selected) >= target_pairs:
                break
            selected.append(candidate)
            selected_ids.add(candidate.get('index'))

    remaining = [
        candidate for candidate in pool
        if candidate.get('index') not in selected_ids and
        candidate.get('type') in {'LP', 'BD'}
    ]
    for candidate in sorted(remaining, key=_candidate_sort_key):
        if len(selected) >= target_pairs:
            break
        selected.append(candidate)
        selected_ids.add(candidate.get('index'))

    if len(selected) < target_pairs:
        warnings.append(
            f'Primary assignment selected {len(selected)} electron pairs; '
            f'target is {target_pairs}.'
        )
    if forbidden_pairs:
        warnings.append('Forbidden bonds were excluded from primary assignment.')

    selected = sorted(selected, key=lambda item: int(item.get('index', 0)))
    occupation_sum = float(sum(candidate.get('occupation', 0.0)
                               for candidate in selected))
    score = occupation_sum - 2.0 * abs(len(selected) - target_pairs)

    return {
        'nbo_list': selected,
        'counts': _candidate_type_counts(selected),
        'score': score,
        'electron_pairs': len(selected),
        'target_electron_pairs': target_pairs,
        'occupation_sum': occupation_sum,
        'warnings': warnings,
    }


def _compatible_pi_combo(combo):
    """Return whether pi candidates form a simple non-overlapping matching."""

    used_atoms = set()
    used_pairs = set()
    for candidate in combo:
        pair = tuple(sorted(candidate.get('atoms', ())))
        if pair in used_pairs:
            return False
        used_pairs.add(pair)
        for atom in pair:
            if atom in used_atoms:
                return False
            used_atoms.add(atom)
    return True


def _enumerate_lewis_alternatives(molecule,
                                  candidates,
                                  primary_assignment,
                                  constraints,
                                  max_alternatives=12,
                                  weight_beta=4.0):
    """Enumerate first Lewis alternatives by varying compatible pi bonds."""

    natoms = molecule.number_of_atoms()
    target_pairs = int(primary_assignment.get(
        'target_electron_pairs',
        round(0.5 * molecule.number_of_electrons()),
    ))
    required_pi_pairs, allowed_pi_pairs, forbidden_pi_pairs = (
        _constraint_pi_pair_sets(constraints, natoms))

    primary_nbos = list(primary_assignment.get('nbo_list', ()))
    fixed = [candidate for candidate in primary_nbos
             if not (candidate.get('type') == 'BD' and
                     candidate.get('subtype') == 'pi')]
    pi_slots = target_pairs - len(fixed)
    warnings = []

    if pi_slots < 0:
        warnings.append('Primary assignment has more non-pi NBOs than electron-pair slots.')
        pi_slots = 0

    pi_pool = []
    for candidate in candidates:
        if candidate.get('type') != 'BD' or candidate.get('subtype') != 'pi':
            continue
        pair = tuple(sorted(candidate.get('atoms', ())))
        if pair in forbidden_pi_pairs:
            continue
        if (candidate.get('non_sigma', False) and
                pair not in required_pi_pairs and pair not in allowed_pi_pairs):
            continue
        pi_pool.append(candidate)

    if len(pi_pool) < pi_slots:
        warnings.append(
            f'Only {len(pi_pool)} pi candidates available for {pi_slots} pi slots.'
        )

    alternatives = []
    for combo in combinations(pi_pool, min(pi_slots, len(pi_pool))):
        combo_pairs = {
            tuple(sorted(candidate.get('atoms', ())))
            for candidate in combo
        }
        if not required_pi_pairs.issubset(combo_pairs):
            continue
        if not _compatible_pi_combo(combo):
            continue

        nbo_list = sorted(fixed + list(combo),
                          key=lambda candidate: int(candidate.get('index', 0)))
        occupation_sum = float(sum(candidate.get('occupation', 0.0)
                                   for candidate in nbo_list))
        allowed_bonus = 0.05 * len(combo_pairs & allowed_pi_pairs)
        score = occupation_sum + allowed_bonus - 2.0 * abs(
            len(nbo_list) - target_pairs)
        alternatives.append({
            'rank': 0,
            'nbo_list': nbo_list,
            'counts': _candidate_type_counts(nbo_list),
            'score': score,
            'weight': 0.0,
            'electron_pairs': len(nbo_list),
            'target_electron_pairs': target_pairs,
            'occupation_sum': occupation_sum,
            'pi_bonds': [tuple(int(atom + 1) for atom in pair)
                         for pair in sorted(combo_pairs)],
            'warnings': list(warnings),
        })

    if not alternatives:
        fallback = dict(primary_assignment)
        fallback['rank'] = 1
        fallback['weight'] = 1.0
        fallback['pi_bonds'] = [
            tuple(int(atom + 1) for atom in sorted(candidate.get('atoms', ())))
            for candidate in primary_nbos
            if candidate.get('type') == 'BD' and candidate.get('subtype') == 'pi'
        ]
        fallback['warnings'] = list(primary_assignment.get('warnings', [])) + warnings
        return [fallback]

    alternatives = sorted(alternatives,
                          key=lambda alt: (-alt['score'],
                                           tuple(alt['pi_bonds'])))
    alternatives = alternatives[:max(1, int(max_alternatives))]
    max_score = max(alt['score'] for alt in alternatives)
    weights = np.array([
        np.exp(float(weight_beta) * (alt['score'] - max_score))
        for alt in alternatives
    ])
    weights_sum = float(np.sum(weights))
    if weights_sum < 1.0e-14:
        weights = np.ones(len(alternatives)) / float(len(alternatives))
    else:
        weights = weights / weights_sum

    for rank, (alternative, weight) in enumerate(zip(alternatives, weights),
                                                 start=1):
        alternative['rank'] = rank
        alternative['weight'] = float(weight)

    return alternatives


def _candidate_vector_from_coefficients(candidate, norb):
    """Rebuild a normalized NAO-space candidate vector from report data."""

    vector = np.zeros(norb)
    for item in candidate.get('coefficients', []):
        vector[int(item['nao_index']) - 1] = float(item['coefficient'])
    norm = float(np.linalg.norm(vector))
    if norm > 1.0e-14:
        vector /= norm
    return vector


def _lewis_density_from_nbos(nbo_list, norb):
    """Build a closed-shell ideal Lewis density from selected NBO vectors."""

    density = np.zeros((norb, norb))
    for candidate in nbo_list:
        vector = _candidate_vector_from_coefficients(candidate, norb)
        if np.linalg.norm(vector) < 1.0e-14:
            continue
        density += 2.0 * np.outer(vector, vector)
    return 0.5 * (density + density.T)


def _nra_subspace_indices(subspace, alternatives, atom_map, angular_map, norb):
    """Return NAO indices used for the requested first-pass NRA fit."""

    subspace = str(subspace).lower()
    if subspace == 'full':
        return np.arange(norb, dtype=int), []

    if subspace == 'pi':
        atoms = set()
        for alternative in alternatives:
            for pair in alternative.get('pi_bonds', []):
                atoms.update(int(atom) - 1 for atom in pair)
        indices = [
            idx for idx in range(norb)
            if int(angular_map[idx]) == 1 and int(atom_map[idx]) in atoms
        ]
        if indices:
            return np.array(indices, dtype=int), []
        return _nra_subspace_indices('selected', alternatives, atom_map,
                                     angular_map, norb)

    selected = set()
    for alternative in alternatives:
        for candidate in alternative.get('nbo_list', []):
            if subspace == 'valence' and candidate.get('type') == 'CR':
                continue
            for item in candidate.get('coefficients', []):
                selected.add(int(item['nao_index']) - 1)

    if subspace in {'selected', 'valence'}:
        if selected:
            return np.array(sorted(selected), dtype=int), []
        return np.arange(norb, dtype=int), [
            f'NRA {subspace} subspace was empty; full NAO space was used.'
        ]

    return np.arange(norb, dtype=int), [
        f"Unknown NRA subspace '{subspace}'; full NAO space was used."
    ]


def _vectorize_symmetric_block(matrix, indices):
    """Vectorize the upper triangle of a selected symmetric matrix block."""

    block = matrix[np.ix_(indices, indices)]
    upper = np.triu_indices(len(indices))
    vector = np.array(block[upper], dtype=float)
    off_diagonal = upper[0] != upper[1]
    vector[off_diagonal] *= np.sqrt(2.0)
    return vector


def _equality_constrained_least_squares(matrix, target):
    """Fit simplex weights with a small active-set least-squares solver."""

    ncols = matrix.shape[1]
    if ncols == 0:
        return np.array([], dtype=float)
    if ncols == 1:
        return np.array([1.0], dtype=float)

    active = list(range(ncols))
    weights = np.zeros(ncols)

    while active:
        active_matrix = matrix[:, active]
        gram = active_matrix.T @ active_matrix
        rhs = active_matrix.T @ target
        ones = np.ones(len(active))
        kkt = np.block([
            [gram, ones[:, np.newaxis]],
            [ones[np.newaxis, :], np.zeros((1, 1))],
        ])
        krhs = np.concatenate([rhs, np.array([1.0])])
        solution = np.linalg.lstsq(kkt, krhs, rcond=None)[0][:-1]

        if np.all(solution >= -1.0e-12):
            weights[active] = np.maximum(solution, 0.0)
            break

        remove_local = int(np.argmin(solution))
        del active[remove_local]

    total = float(np.sum(weights))
    if total < 1.0e-14:
        weights[:] = 1.0 / float(ncols)
    else:
        weights /= total
    return weights


def _build_nra_results(molecule, nao_data, alternatives, subspace, fit_metric):
    """Fit ideal Lewis densities to the actual NAO density."""

    fit_metric = str(fit_metric).lower()
    warnings = []
    norb = nao_data.density.shape[0]

    if fit_metric != 'frobenius':
        warnings.append(
            f"Unknown NRA fit metric '{fit_metric}'; Frobenius metric was used."
        )
        fit_metric = 'frobenius'

    if molecule.get_multiplicity() != 1:
        return {
            'subspace': str(subspace).lower(),
            'fit_metric': fit_metric,
            'weights': [],
            'residual_norm': None,
            'relative_residual': None,
            'structures': [],
            'warnings': [
                'NRA is currently implemented only for closed-shell singlet alternatives.'
            ],
        }

    if not alternatives:
        return {
            'subspace': str(subspace).lower(),
            'fit_metric': fit_metric,
            'weights': [],
            'residual_norm': None,
            'relative_residual': None,
            'structures': [],
            'warnings': ['No Lewis alternatives were available for NRA fitting.'],
        }

    indices, subspace_warnings = _nra_subspace_indices(
        subspace,
        alternatives,
        nao_data.atom_map,
        nao_data.angular_momentum_map,
        norb,
    )
    warnings.extend(subspace_warnings)

    target = _vectorize_symmetric_block(nao_data.density, indices)
    lewis_densities = [
        _lewis_density_from_nbos(alternative.get('nbo_list', []), norb)
        for alternative in alternatives
    ]
    fit_matrix = np.column_stack([
        _vectorize_symmetric_block(density, indices)
        for density in lewis_densities
    ])
    weights = _equality_constrained_least_squares(fit_matrix, target)
    model = fit_matrix @ weights if len(weights) else np.zeros_like(target)
    residual = target - model
    residual_norm = float(np.linalg.norm(residual))
    target_norm = float(np.linalg.norm(target))
    relative_residual = (residual_norm / target_norm
                         if target_norm > 1.0e-14 else 0.0)

    structures = []
    for alternative, weight, lewis_density in zip(alternatives,
                                                  weights,
                                                  lewis_densities):
        structure_vector = _vectorize_symmetric_block(lewis_density, indices)
        structures.append({
            'rank': int(alternative.get('rank', 0)),
            'pi_bonds': list(alternative.get('pi_bonds', [])),
            'nra_weight': float(weight),
            'score_weight': float(alternative.get('weight', 0.0)),
            'score': float(alternative.get('score', 0.0)),
            'residual_norm': float(np.linalg.norm(target - structure_vector)),
        })

    return {
        'subspace': str(subspace).lower(),
        'fit_metric': fit_metric,
        'weights': [float(weight) for weight in weights],
        'residual_norm': residual_norm,
        'relative_residual': relative_residual,
        'structures': structures,
        'warnings': warnings,
    }


class NboDriver:
    """Natural Bond Orbital analysis driver restart scaffold."""

    def __init__(self, comm=None, ostream=None):
        if comm is None:
            comm = MPI.COMM_WORLD

        if ostream is None:
            if comm.Get_rank() == mpi_master():
                ostream = OutputStream(sys.stdout)
            else:
                ostream = OutputStream(None)

        self.comm = comm
        self.rank = comm.Get_rank()
        self.nodes = comm.Get_size()
        self.ostream = ostream

        self.npa_report_level = 'summary'
        self.nbo_report_level = 'summary'
        self.verbose = True
        self._last_molecule = None
        self._last_results = None

    @staticmethod
    def _coerce_dataclass(value, cls, name):
        if value is None:
            return cls()
        if isinstance(value, cls):
            return value
        if isinstance(value, dict):
            valid = {field.name for field in fields(cls)}
            invalid = sorted(set(value) - valid)
            assert_msg_critical(not invalid,
                                f"Invalid {name} keys: {invalid}. "
                                f"Allowed keys are {sorted(valid)}")
            return cls(**value)
        raise TypeError(f'{name} must be None, dict, or {cls.__name__}')

    def compute(self,
                molecule,
                basis,
                mol_orbs,
                mode='npa',
                options=None,
                constraints=None):
        """Compute the clean NAO/NPA baseline.

        Full Lewis/NBO assignment is intentionally not implemented in this
        restart scaffold. ``'primary'`` is accepted as a compatibility alias
        for ``'npa'`` and does not perform Lewis/NBO assignment.
        """

        assert_msg_critical(mode in {'npa', 'primary'},
                    "NBO restart scaffold: only mode='npa' or "
                    "mode='primary' is implemented")

        compute_options = self._coerce_dataclass(options,
                                                 NboComputeOptions,
                                                 'NBO options')
        compute_constraints = self._coerce_dataclass(constraints,
                                                     NboConstraints,
                                                     'NBO constraints')
        results = {}

        if self.rank == mpi_master():
            natoms = molecule.number_of_atoms()
            labels = molecule.get_labels()
            nuclear_charges = np.array([
                chemical_element_identifier(labels[atom])
                for atom in range(natoms)
            ], dtype=float)

            overlap = np.array(compute_overlap_integrals(molecule, basis))
            density = _build_density_matrix(mol_orbs)
            ao_to_atom, ao_to_l = _ao_shell_map(molecule, basis)
            equivalent_groups = _find_equivalent_atom_groups(molecule)
            local_frames = _build_local_atom_frames(molecule)
            nao_data = _build_orthonormal_nao_data(overlap,
                                                   density,
                                                   ao_to_atom,
                                                   ao_to_l,
                                                   natoms,
                                                   equivalent_groups,
                                                   local_frames)
            mo_analysis = (_build_mo_nao_analysis(
                molecule,
                mol_orbs,
                overlap,
                nao_data,
                top_n=compute_options.mo_analysis_top,
                threshold=compute_options.mo_analysis_threshold,
            ) if compute_options.include_mo_analysis else {})
            nbo_candidates = (_build_nbo_candidates(
                molecule,
                nao_data,
                lone_pair_min_occupation=compute_options.lone_pair_min_occupation,
                bond_min_occupation=compute_options.bond_min_occupation,
                bond_min_atom_weight=compute_options.bond_min_atom_weight,
                pi_min_occupation=compute_options.pi_min_occupation,
                conjugated_pi_max_path=compute_options.conjugated_pi_max_path,
                constraints=compute_constraints,
            ) if compute_options.include_nbo_candidates else [])
            candidate_counts = _candidate_type_counts(nbo_candidates)
            primary_assignment = (_build_primary_assignment(
                molecule,
                nbo_candidates,
                compute_constraints,
            ) if compute_options.include_lewis_assignment else {
                'nbo_list': nbo_candidates,
                'counts': candidate_counts,
                'score': 0.0,
                'electron_pairs': len(nbo_candidates),
                'target_electron_pairs': int(round(
                    0.5 * molecule.number_of_electrons())),
                'occupation_sum': float(sum(
                    candidate.get('occupation', 0.0)
                    for candidate in nbo_candidates)),
                'warnings': [
                    'NBO candidate layer only: Lewis/NBO assignment was not requested.'
                ],
            })
            alternatives = _enumerate_lewis_alternatives(
                molecule,
                nbo_candidates,
                primary_assignment,
                compute_constraints,
                max_alternatives=compute_options.max_alternatives,
                weight_beta=compute_options.lewis_weight_beta,
            ) if compute_options.include_lewis_assignment else []
            if alternatives:
                primary_ids = {
                    candidate.get('index')
                    for candidate in primary_assignment.get('nbo_list', [])
                }
                for alternative in alternatives:
                    alternative_ids = {
                        candidate.get('index')
                        for candidate in alternative.get('nbo_list', [])
                    }
                    if alternative_ids == primary_ids:
                        primary_assignment['weight'] = alternative.get('weight', 0.0)
                        primary_assignment['rank'] = alternative.get('rank', 0)
                        break
            primary_counts = primary_assignment['counts']
            nra_results = (_build_nra_results(
                molecule,
                nao_data,
                alternatives,
                compute_options.nra_subspace,
                compute_options.nra_fit_metric,
            ) if compute_options.include_nra else None)

            natural_charges = nuclear_charges.copy()
            for nao_index, atom in enumerate(nao_data.atom_map):
                natural_charges[atom] -= nao_data.populations[nao_index]

            results = {
                'nao_transform': nao_data.transform,
                'nao_density_matrix': nao_data.density,
                'nao_overlap_matrix': nao_data.overlap,
                'nao_populations': nao_data.populations,
                'nao_atom_map': nao_data.atom_map,
                'nao_l_map': nao_data.angular_momentum_map,
                'natural_charges': natural_charges,
                'mo_analysis': mo_analysis,
                'nbo_candidates': nbo_candidates,
                'nbo_list': primary_assignment['nbo_list'],
                'primary': primary_assignment,
                'alternatives': alternatives,
                'diagnostics': {
                    'electron_count': float(np.trace(nao_data.density).real),
                    'orthonormality_error': float(np.linalg.norm(
                        nao_data.overlap - np.eye(nao_data.overlap.shape[0]))),
                    'equivalent_atom_groups': [
                        tuple(int(atom + 1) for atom in group)
                        for group in nao_data.equivalent_atom_groups
                    ],
                    'local_frames': [frame.tolist()
                                     for frame in nao_data.local_frames],
                    'mo_nao_max_normalization_error': float(
                        mo_analysis.get('max_normalization_error', 0.0)),
                    'nbo_candidate_counts': candidate_counts,
                    'primary_nbo_counts': primary_counts,
                    'constraint_summary': {
                        key: len(value) if hasattr(value, '__len__') else value
                        for key, value in asdict(compute_constraints).items()
                    },
                } if compute_options.include_diagnostics else {},
                'provenance': {
                    'api_version': 'nbo-restart-v0',
                    'mode': mode,
                    'options': asdict(compute_options),
                    'constraints': asdict(compute_constraints),
                },
            }
            if nra_results is not None:
                results['nra'] = nra_results

            self._last_molecule = molecule
            self._last_results = results

            if self.verbose:
                self.print_npa_report(molecule, results, level=self.npa_report_level)
                self.print_nbo_report(molecule, results, level=self.nbo_report_level)

        return results

    def set_npa_report_level(self, level):
        """Set NPA report level."""

        valid = {'none', 'summary', 'standard', 'full'}
        assert_msg_critical(level in valid,
                            f"Invalid NPA report level '{level}'. "
                            f"Choose one of {sorted(valid)}")
        self.npa_report_level = level

    def set_nbo_report_level(self, level):
        """Set NBO report level."""

        valid = {'none', 'summary', 'full'}
        assert_msg_critical(level in valid,
                            f"Invalid NBO report level '{level}'. "
                            f"Choose one of {sorted(valid)}")
        self.nbo_report_level = level

    def get_npa_data(self, molecule, results):
        """Return structured NPA reporting data."""

        labels = molecule.get_labels()
        natoms = molecule.number_of_atoms()
        charges = np.array(results['natural_charges'], dtype=float)
        populations = np.array(results['nao_populations'], dtype=float)
        atom_map = np.array(results['nao_atom_map'], dtype=int)
        angular_map = np.array(results['nao_l_map'], dtype=int)

        core_pop = np.zeros(natoms)
        valence_pop = np.zeros(natoms)
        rydberg_pop = np.zeros(natoms)
        valence_by_l = np.zeros((natoms, 6))
        rows = []
        l_labels = {0: 's', 1: 'p', 2: 'd', 3: 'f', 4: 'g', 5: 'h'}

        for atom in range(natoms):
            indices = np.where(atom_map == atom)[0]
            ordered = indices[np.argsort(populations[indices])[::-1]]
            nuclear_charge = chemical_element_identifier(labels[atom])
            n_core = min(_n_core_orbitals_from_z(nuclear_charge), len(ordered))
            n_val = _valence_principal_n_from_z(nuclear_charge)
            core_indices = _select_core_indices(ordered,
                                                populations,
                                                angular_map,
                                                n_core)

            for rank, nao_index in enumerate(ordered):
                occ = float(populations[nao_index])
                l_val = int(angular_map[nao_index])
                l_symbol = l_labels.get(l_val, '?')

                if nao_index in core_indices:
                    partition = 'Core'
                    core_pop[atom] += occ
                    orbital_type = f"Cor({max(1, n_val - 1)}{l_symbol})"
                elif occ > 1.0e-3:
                    partition = 'Valence'
                    valence_pop[atom] += occ
                    orbital_type = f"Val({n_val}{l_symbol})"
                    if 0 <= l_val < valence_by_l.shape[1]:
                        valence_by_l[atom, l_val] += occ
                else:
                    partition = 'Rydberg'
                    rydberg_pop[atom] += occ
                    orbital_type = f"Ryd({n_val + 1}{l_symbol})"

                rows.append({
                    'nao_index': int(nao_index + 1),
                    'atom_label': labels[atom],
                    'atom_index': int(atom + 1),
                    'lang': l_symbol,
                    'type_ao': orbital_type,
                    'partition': partition,
                    'occupancy': occ,
                })

        return {
            'labels': labels,
            'charges': charges,
            'core_pop': core_pop,
            'valence_pop': valence_pop,
            'rydberg_pop': rydberg_pop,
            'total_pop': core_pop + valence_pop + rydberg_pop,
            'valence_by_l': valence_by_l,
            'nao_rows': rows,
            'n_elec': float(molecule.number_of_electrons()),
        }

    def format_npa_report(self, molecule, results, level=None):
        """Return an NPA report string."""

        if level is None:
            level = self.npa_report_level
        valid = {'none', 'summary', 'standard', 'full'}
        assert_msg_critical(level in valid,
                            f"Invalid NPA report level '{level}'. "
                            f"Choose one of {sorted(valid)}")
        if level == 'none':
            return ''

        data = self.get_npa_data(molecule, results)
        lines = []
        width = 74

        if level == 'full':
            lines.append('NATURAL POPULATIONS: Natural atomic orbital occupancies')
            lines.append('')
            lines.append(f"{'NAO':>4} {'Atom #':>7}  {'l':>4}   {'Type(AO)':<10}  {'Occupancy':>10}")
            lines.append('-' * width)
            for row in data['nao_rows']:
                lines.append(
                    f"{row['nao_index']:>4} {row['atom_label']:>2} {row['atom_index']:>2}"
                    f"  {row['lang']:>4}   {row['type_ao']:<10}  {row['occupancy']:>10.5f}"
                )
            lines.append('')

        lines.append('A Summary of Natural Population Analysis:')
        lines.append('')
        lines.append(f"{'Atom #':>8}  {'Natural Charge':>14}  {'Core':>10}  {'Valence':>10}  {'Rydberg':>10}  {'Total':>10}")
        lines.append('-' * width)

        for atom, label in enumerate(data['labels']):
            lines.append(
                f"{label:>2} {atom + 1:>2}  {data['charges'][atom]:>+14.5f}"
                f"  {data['core_pop'][atom]:>10.5f}"
                f"  {data['valence_pop'][atom]:>10.5f}"
                f"  {data['rydberg_pop'][atom]:>10.5f}"
                f"  {data['total_pop'][atom]:>10.5f}"
            )

        lines.append('=' * width)
        lines.append(
            f"{'* Total *':>10}  {data['charges'].sum():>+14.5f}"
            f"  {data['core_pop'].sum():>10.5f}"
            f"  {data['valence_pop'].sum():>10.5f}"
            f"  {data['rydberg_pop'].sum():>10.5f}"
            f"  {data['total_pop'].sum():>10.5f}"
        )

        if level in {'standard', 'full'}:
            lines.append('')
            lines.append('Natural electron configuration (effective valence):')
            lines.append('-' * width)
            l_symbols = ['s', 'p', 'd', 'f', 'g', 'h']
            for atom, label in enumerate(data['labels']):
                nuclear_charge = chemical_element_identifier(label)
                n_val = _valence_principal_n_from_z(nuclear_charge)
                pieces = []
                for l_val, symbol in enumerate(l_symbols):
                    pop = data['valence_by_l'][atom, l_val]
                    if pop > 1.0e-4:
                        pieces.append(f"{n_val}{symbol}({pop:5.2f})")
                config = ''.join(pieces) if pieces else 'n/a'
                lines.append(f'{label:>2} {atom + 1:>2}   {config}')

        return '\n'.join(lines)

    def print_npa_report(self, molecule, results, level=None, ostream=None):
        """Print an NPA report."""

        text = self.format_npa_report(molecule, results, level=level)
        if not text:
            return
        if ostream is None:
            ostream = self.ostream
        if hasattr(ostream, 'print_header'):
            ostream.print_blank()
            for line in text.splitlines():
                ostream.print_header(line)
            ostream.print_blank()
        else:
            print(text, file=ostream if ostream is not None else sys.stdout)

    def npa_report(self,
                   level=None,
                   molecule=None,
                   results=None,
                   ostream=None,
                   return_text=False):
        """Convenience NPA reporting API."""

        molecule = self._last_molecule if molecule is None else molecule
        results = self._last_results if results is None else results
        assert_msg_critical(molecule is not None,
                            'npa_report: no molecule provided and no prior compute() context available')
        assert_msg_critical(results is not None,
                            'npa_report: no results provided and no prior compute() context available')
        text = self.format_npa_report(molecule, results, level=level)
        if return_text:
            return text
        if text:
            print(text, file=ostream if ostream is not None else sys.stdout)
        return text


    def get_mo_data(self, molecule, results):
        """Return structured MO-in-NAO analysis data."""

        assert_msg_critical('mo_analysis' in results,
                            'get_mo_data: results do not contain MO analysis')
        mo_analysis = results['mo_analysis']
        assert_msg_critical(bool(mo_analysis),
                            'get_mo_data: MO analysis was not requested')
        return mo_analysis

    def format_mo_report(self, molecule, results, level=None):
        """Return a molecular-orbital composition report in the NAO basis."""

        if level is None:
            level = 'occupied'
        valid = {'none', 'occupied', 'all'}
        assert_msg_critical(level in valid,
                            f"Invalid MO report level '{level}'. "
                            f"Choose one of {sorted(valid)}")
        if level == 'none':
            return ''

        data = self.get_mo_data(molecule, results)
        lines = []
        width = 92
        lines.append('MO ANALYSIS: Molecular orbital composition in the NAO basis')
        lines.append('=' * width)
        lines.append('Weights are squared NAO expansion coefficients and sum to one per MO.')
        lines.append('')

        for block in data['spin_blocks']:
            spin_label = block['spin']
            orbitals = block['orbitals']
            selected = [orb for orb in orbitals if level == 'all' or
                        orb['occupation'] > 1.0e-8]
            if not selected:
                continue

            lines.append(f'Spin block: {spin_label}')
            lines.append('-' * width)
            lines.append(
                f"{'MO':>4} {'Occ':>8} {'Energy/a.u.':>14}  "
                f"{'Dominant atoms':<26} {'Dominant NAOs'}"
            )
            lines.append('-' * width)

            for orbital in selected:
                atom_pieces = []
                for item in orbital['top_atoms'][:3]:
                    atom_pieces.append(
                        f"{item['atom_label']}{item['atom_index']}:{100.0 * item['weight']:5.1f}%"
                    )
                nao_pieces = []
                for item in orbital['top_naos']:
                    nao_pieces.append(
                        f"{item['nao_index']}({item['atom_label']}{item['atom_index']}{item['l']}):"
                        f"{100.0 * item['weight']:4.1f}%"
                    )
                lines.append(
                    f"{orbital['mo_index']:>4} {orbital['occupation']:>8.4f} "
                    f"{orbital['energy']:>14.6f}  "
                    f"{', '.join(atom_pieces):<26} "
                    f"{', '.join(nao_pieces)}"
                )
            lines.append('')

        return '\n'.join(lines).rstrip()

    def print_mo_report(self, molecule, results, level=None, ostream=None):
        """Print a molecular-orbital composition report."""

        text = self.format_mo_report(molecule, results, level=level)
        if not text:
            return
        if ostream is None:
            ostream = self.ostream
        if hasattr(ostream, 'print_header'):
            ostream.print_blank()
            for line in text.splitlines():
                ostream.print_header(line)
            ostream.print_blank()
        else:
            print(text, file=ostream if ostream is not None else sys.stdout)

    def mo_report(self,
                  level=None,
                  molecule=None,
                  results=None,
                  ostream=None,
                  return_text=False):
        """Convenience MO-in-NAO reporting API."""

        molecule = self._last_molecule if molecule is None else molecule
        results = self._last_results if results is None else results
        assert_msg_critical(molecule is not None,
                            'mo_report: no molecule provided and no prior compute() context available')
        assert_msg_critical(results is not None,
                            'mo_report: no results provided and no prior compute() context available')
        text = self.format_mo_report(molecule, results, level=level)
        if return_text:
            return text
        if text:
            print(text, file=ostream if ostream is not None else sys.stdout)
        return text

    def format_nbo_report(self, molecule, results, level=None):
        """Return an NBO report string."""

        if level is None:
            level = self.nbo_report_level
        valid = {'none', 'summary', 'full'}
        assert_msg_critical(level in valid,
                            f"Invalid NBO report level '{level}'. "
                            f"Choose one of {sorted(valid)}")
        if level == 'none':
            return ''

        assigned_nbos = results.get('nbo_list', [])
        candidates = results.get('nbo_candidates', assigned_nbos)
        primary = results.get('primary', {})
        counts = _candidate_type_counts(assigned_nbos)
        candidate_counts = _candidate_type_counts(candidates)
        labels = molecule.get_labels()
        lines = []
        width = 92

        lines.append('Natural Bond Orbital (NBO) Primary Summary')
        lines.append('=' * width)
        lines.append('First-pass Lewis-like primary assignment from generated NBO candidates.')
        lines.append('Implemented now: NAO/NPA, MO-in-NAO analysis, candidates, and primary selection.')
        lines.append('Report order: CR, BD(sigma), BD(pi), LP, RY, then antibonding/other candidates.')
        lines.append('')
        if counts:
            count_text = ', '.join(f'{key}={value}'
                                   for key, value in sorted(counts.items()))
        else:
            count_text = 'none'
        if candidate_counts:
            candidate_count_text = ', '.join(
                f'{key}={value}' for key, value in sorted(candidate_counts.items()))
        else:
            candidate_count_text = 'none'
        lines.append(f'Primary counts: {count_text}')
        lines.append(f'Candidate pool: {candidate_count_text}')
        if primary:
            lines.append(
                f"Selected pairs: {primary.get('electron_pairs', len(assigned_nbos))}/"
                f"{primary.get('target_electron_pairs', '?')}  "
                f"occupation sum={primary.get('occupation_sum', 0.0):.5f}  "
                f"score={primary.get('score', 0.0):.5f}"
            )
            if 'weight' in primary:
                lines.append(
                    f"Primary resonance rank={primary.get('rank', '?')}  "
                    f"weight={primary.get('weight', 0.0):.5f}"
                )
            for warning in primary.get('warnings', []):
                lines.append(f'Warning: {warning}')

        if level == 'summary':
            lines.append('=' * width)
            return '\n'.join(lines)

        lines.append('')
        lines.append(
            f"{'NBO':>4} {'Type':<10} {'Atoms':<14} {'Occ':>10}  "
            f"{'Polarization':<24} {'Source'}"
        )
        lines.append('-' * width)

        for candidate in sorted(assigned_nbos,
                                key=lambda item: _nbo_report_sort_key(item,
                                                                      labels)):
            atom_text = _report_atom_text(candidate, labels)
            candidate_type = candidate.get('type', '?')
            subtype = candidate.get('subtype', '')
            type_text = candidate_type if not subtype else f'{candidate_type}({subtype})'
            polarization = candidate.get('polarization', {})
            pol_text = ', '.join(
                f'{labels[atom - 1]}{atom}:{100.0 * weight:5.1f}%'
                for atom, weight in polarization.items()
            )
            lines.append(
                f"{candidate.get('index', 0):>4} {type_text:<10.10} "
                f"{atom_text:<14} {candidate.get('occupation', 0.0):>10.5f}  "
                f"{pol_text:<24} {candidate.get('source', '')}"
            )

            if level == 'full':
                pieces = []
                for item in candidate.get('coefficients', []):
                    nao_index = item['nao_index']
                    coefficient = item['coefficient']
                    pieces.append(f"{nao_index}:{coefficient:+.4f}")
                if pieces:
                    lines.append(f"{'':>4} {'coef':<10} {'':<14} {'':>10}  "
                                 f"{', '.join(pieces)}")

        alternatives = results.get('alternatives', [])
        if level == 'full' and alternatives:
            lines.append('')
            lines.append('Lewis/resonance alternatives')
            lines.append('-' * width)
            lines.append(
                f"{'Rank':>4} {'Weight':>10} {'Score':>12} "
                f"{'Pairs':>7} {'Pi bonds'}"
            )
            lines.append('-' * width)
            for alternative in alternatives:
                pi_text = ', '.join(
                    f"{pair[0]}-{pair[1]}"
                    for pair in alternative.get('pi_bonds', [])
                )
                lines.append(
                    f"{alternative.get('rank', 0):>4} "
                    f"{alternative.get('weight', 0.0):>10.5f} "
                    f"{alternative.get('score', 0.0):>12.5f} "
                    f"{alternative.get('electron_pairs', 0):>7}/"
                    f"{alternative.get('target_electron_pairs', '?'):<3} "
                    f"{pi_text}"
                )

        lines.append('=' * width)
        return '\n'.join(lines)

    def print_nbo_report(self, molecule, results, level=None, ostream=None):
        """Print an NBO report."""

        text = self.format_nbo_report(molecule, results, level=level)
        if not text:
            return
        if ostream is None:
            ostream = self.ostream
        if hasattr(ostream, 'print_header'):
            for line in text.splitlines():
                ostream.print_header(line)
            ostream.print_blank()
        else:
            print(text, file=ostream if ostream is not None else sys.stdout)

    def nbo_report(self,
                   level=None,
                   molecule=None,
                   results=None,
                   ostream=None,
                   return_text=False):
        """Convenience NBO reporting API."""

        molecule = self._last_molecule if molecule is None else molecule
        results = self._last_results if results is None else results
        assert_msg_critical(molecule is not None,
                            'nbo_report: no molecule provided and no prior compute() context available')
        assert_msg_critical(results is not None,
                            'nbo_report: no results provided and no prior compute() context available')
        text = self.format_nbo_report(molecule, results, level=level)
        if return_text:
            return text
        if text:
            print(text, file=ostream if ostream is not None else sys.stdout)
        return text