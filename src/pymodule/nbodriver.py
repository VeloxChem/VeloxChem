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

This module implements the VeloxChem NBO analysis pipeline: NAO/NPA
construction, MO composition in the NAO basis, NBO candidate generation,
Lewis/resonance assignment, donor-acceptor diagnostics, and optional NRA/NRT
density fitting for closed-shell and open-shell alternatives.
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
from .orbitalanalyzerdriver import OrbitalAnalyzer, OrbitalAnalyzerOptions


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
    rydberg_max_occupation: float = 0.50
    bond_min_occupation: float = 1.20
    bond_min_atom_weight: float = 0.10
    pi_min_occupation: float = 0.20
    conjugated_pi_max_path: int = 2
    max_alternatives: int = 12
    lewis_weight_beta: float = 4.0
    include_nra: bool = False
    nra_subspace: str = 'selected'
    nra_fit_metric: str = 'frobenius'
    nra_prior_weights: object = None
    nra_prior_strength: float = 0.0
    nra_prior_mode: str = 'regularized'
    nra_spin_fit: str = 'total_spin'


@dataclass(frozen=True)
class NboConstraints:
    """Constraint container for Lewis/NBO assignment and resonance search."""

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


@dataclass(frozen=True)
class SpinNaoData:
    """Spin-resolved NAO density data."""

    alpha_density: np.ndarray
    beta_density: np.ndarray
    spin_density: np.ndarray
    spin_populations: np.ndarray
    unpaired_electrons: float


_ELECTRON_PAIR_PENALTY_WEIGHT = 2.0
_FORMAL_CHARGE_PENALTY_WEIGHT = 1.0
_OCTET_PENALTY_WEIGHT = 2.0
_VALENCE_WARNING_PENALTY_WEIGHT = 1.0
_NONCLASSICAL_PENALTY_WEIGHT = 0.05
_ACCEPTOR_TYPES = {'BD*', 'RY'}
_SPIN_ACTIVE_MIN_OCCUPATION = 0.20


def _foundation_invariant_status(value, tolerance):
    """Return a compact pass/fail diagnostic entry for one invariant."""

    value = float(value)
    tolerance = float(tolerance)
    return {
        'value': value,
        'tolerance': tolerance,
        'passed': bool(value <= tolerance),
    }


def _build_foundation_diagnostics(molecule,
                                  overlap,
                                  ao_density,
                                  nao_data,
                                  natural_charges,
                                  mo_analysis,
                                  candidate_counts,
                                  primary_counts,
                                  compute_constraints):
    """Build stable diagnostics for the NAO/NPA foundation invariants."""

    nao_overlap = nao_data.transform.T @ overlap @ nao_data.transform
    nao_density = (nao_data.transform.T @ overlap @ ao_density @ overlap @
                   nao_data.transform)
    nao_density = 0.5 * (nao_density + nao_density.T)
    identity = np.eye(nao_overlap.shape[0])

    electron_count = float(np.trace(nao_data.density).real)
    target_electron_count = float(molecule.number_of_electrons())
    charge_sum = float(np.sum(natural_charges).real)
    molecular_charge = float(molecule.get_charge())

    orthonormality_error = float(np.linalg.norm(nao_data.overlap - identity))
    density_formula_error = float(np.linalg.norm(nao_data.density - nao_density))
    density_symmetry_error = float(np.linalg.norm(nao_data.density -
                                                  nao_data.density.T))
    electron_count_error = abs(electron_count - target_electron_count)
    charge_conservation_error = abs(charge_sum - molecular_charge)
    population_trace_error = abs(float(np.sum(nao_data.populations)) -
                                 electron_count)
    mo_normalization_error = float(
        mo_analysis.get('max_normalization_error', 0.0))

    tolerances = {
        'orthonormality_error': 1.0e-8,
        'density_formula_error': 1.0e-8,
        'density_symmetry_error': 1.0e-10,
        'electron_count_error': 1.0e-8,
        'charge_conservation_error': 1.0e-8,
        'population_trace_error': 1.0e-10,
        'mo_nao_max_normalization_error': 1.0e-8,
    }
    foundation_invariants = {
        'orthonormality': _foundation_invariant_status(
            orthonormality_error, tolerances['orthonormality_error']),
        'density_transform': _foundation_invariant_status(
            density_formula_error, tolerances['density_formula_error']),
        'density_symmetry': _foundation_invariant_status(
            density_symmetry_error, tolerances['density_symmetry_error']),
        'electron_count': _foundation_invariant_status(
            electron_count_error, tolerances['electron_count_error']),
        'charge_conservation': _foundation_invariant_status(
            charge_conservation_error, tolerances['charge_conservation_error']),
        'population_trace': _foundation_invariant_status(
            population_trace_error, tolerances['population_trace_error']),
        'mo_normalization': _foundation_invariant_status(
            mo_normalization_error,
            tolerances['mo_nao_max_normalization_error']),
    }

    return {
        'electron_count': electron_count,
        'target_electron_count': target_electron_count,
        'electron_count_error': electron_count_error,
        'natural_charge_sum': charge_sum,
        'molecular_charge': molecular_charge,
        'charge_conservation_error': charge_conservation_error,
        'orthonormality_error': orthonormality_error,
        'density_formula_error': density_formula_error,
        'density_symmetry_error': density_symmetry_error,
        'population_trace_error': population_trace_error,
        'foundation_tolerances': tolerances,
        'foundation_invariants': foundation_invariants,
        'foundation_invariants_passed': all(
            item['passed'] for item in foundation_invariants.values()),
        'equivalent_atom_groups': [
            tuple(int(atom + 1) for atom in group)
            for group in nao_data.equivalent_atom_groups
        ],
        'local_frames': [frame.tolist() for frame in nao_data.local_frames],
        'mo_nao_max_normalization_error': mo_normalization_error,
        'nbo_candidate_counts': candidate_counts,
        'primary_nbo_counts': primary_counts,
        'constraint_summary': {
            key: len(value) if hasattr(value, '__len__') else value
            for key, value in asdict(compute_constraints).items()
        },
    }


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


# Use OrbitalAnalyzer for AO mapping
def _ao_shell_map(molecule, basis):
    ocd = OrbitalAnalyzer(molecule, basis)
    diag = ocd.run()
    return diag.ao_to_atom, diag.ao_to_l


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
    local_count = rotations.shape[0]
    priority_values = (np.asarray(angular_labels, dtype=float) *
                       float(local_count + 1) +
                       np.arange(local_count, dtype=float))
    priority_operator = np.diag(priority_values)

    def fix_column_sign(column):
        pivot = int(np.argmax(np.abs(column)))
        if column[pivot] < 0.0:
            return -column
        return column

    start = 0
    while start < len(occupations):
        stop = start + 1
        while (stop < len(occupations) and
               abs(float(occupations[stop] - occupations[start])) < threshold):
            stop += 1

        if stop - start > 1:
            subspace = rotations[:, start:stop]
            projected_priority = subspace.T @ priority_operator @ subspace
            tie_values, tie_vectors = np.linalg.eigh(projected_priority)
            order = np.argsort(tie_values)
            canonical_subspace = subspace @ tie_vectors[:, order]
            for offset in range(stop - start):
                rotations[:, start + offset] = fix_column_sign(
                    canonical_subspace[:, offset])
        else:
            rotations[:, start] = fix_column_sign(rotations[:, start])

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


def _spin_ao_density_matrices(mol_orbs):
    """Build alpha and beta AO density matrices from molecular orbitals."""

    from .molecularorbitals import molorb

    def spin_density(occupations, coefficients):
        occ = np.asarray(occupations)
        coeff = np.asarray(coefficients)
        selected = occ > 1.0e-12
        return coeff[:, selected] @ np.diag(occ[selected]) @ coeff[:, selected].T

    alpha_density = spin_density(mol_orbs.occa_to_numpy(),
                                 mol_orbs.alpha_to_numpy())
    if mol_orbs.get_orbitals_type() == molorb.rest:
        return alpha_density, alpha_density.copy()

    beta_density = spin_density(mol_orbs.occb_to_numpy(),
                                mol_orbs.beta_to_numpy())
    return alpha_density, beta_density


def _build_spin_nao_data(mol_orbs, overlap, nao_data):
    """Build spin-resolved densities in the NAO basis."""

    alpha_ao, beta_ao = _spin_ao_density_matrices(mol_orbs)
    alpha_nao = nao_data.transform.T @ overlap @ alpha_ao @ overlap @ nao_data.transform
    beta_nao = nao_data.transform.T @ overlap @ beta_ao @ overlap @ nao_data.transform
    alpha_nao = 0.5 * (alpha_nao + alpha_nao.T)
    beta_nao = 0.5 * (beta_nao + beta_nao.T)
    spin_density = alpha_nao - beta_nao
    spin_populations = np.diag(spin_density)
    unpaired_electrons = float(np.trace(spin_density).real)

    return SpinNaoData(alpha_density=alpha_nao,
                       beta_density=beta_nao,
                       spin_density=spin_density,
                       spin_populations=spin_populations,
                       unpaired_electrons=unpaired_electrons)


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


def _candidate_vector_from_local(pair_indices, local_vector, norb):
    """Return a normalized full NAO vector from pair-local coefficients."""

    components = {
        int(index): float(coefficient)
        for index, coefficient in zip(pair_indices, local_vector)
    }
    return _normal_orbital_candidate_vector(norb, components)


def _atom_side_weights(local_vector, left_size):
    """Return left/right weights for a two-center local vector."""

    left_weight = float(np.sum(local_vector[:left_size]**2))
    right_weight = float(np.sum(local_vector[left_size:]**2))
    return left_weight, right_weight


def _append_antibonding_complement(candidates,
                                   serial,
                                   parent,
                                   occupations,
                                   rotations,
                                   pair_indices,
                                   left_size,
                                   density,
                                   atom_i,
                                   atom_j,
                                   bond_min_atom_weight,
                                   norb,
                                   used_roots,
                                   source):
    """Append a same-subspace BD* complement for an occupied BD candidate."""

    parent_vector = _candidate_vector_from_coefficients(parent, norb)
    for root in np.argsort(occupations):
        root = int(root)
        if root in used_roots:
            continue

        local_vector = rotations[:, root]
        left_weight, right_weight = _atom_side_weights(local_vector,
                                                       left_size)
        if (left_weight < bond_min_atom_weight or
                right_weight < bond_min_atom_weight):
            continue

        vector = _candidate_vector_from_local(pair_indices,
                                              local_vector,
                                              norb)
        if abs(float(np.dot(parent_vector, vector))) > 1.0e-10:
            continue

        candidates.append({
            'index': serial,
            'type': 'BD*',
            'subtype': parent.get('subtype', ''),
            'atoms': (int(atom_i), int(atom_j)),
            'non_sigma': bool(parent.get('non_sigma', False)),
            'electron_count': 0.0,
            'occupation': _candidate_occupation(vector, density),
            'polarization': {
                int(atom_i + 1): left_weight,
                int(atom_j + 1): right_weight,
            },
            'coefficients': _candidate_coefficients(vector),
            'source': source,
            'parent_index': int(parent.get('index', 0)),
            'parent_type': parent.get('type', ''),
            'parent_subtype': parent.get('subtype', ''),
            'parent_overlap': float(np.dot(parent_vector, vector)),
        })
        used_roots.add(root)
        return serial + 1

    return serial


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


def _polar_resonance_lone_pair_atoms(molecule, atom_map, angular_map):
    """Return terminal atoms in generic polar pi-resonance fragments."""

    labels = molecule.get_labels()
    connectivity = molecule.get_connectivity_matrix()
    natoms = molecule.number_of_atoms()
    nuclear_charges = [
        chemical_element_identifier(labels[atom]) for atom in range(natoms)
    ]
    resonance_atoms = set()

    for center in range(natoms):
        center_charge = nuclear_charges[center]
        if center_charge < 5:
            continue
        if not _pi_capable_atom(center, atom_map, angular_map):
            continue

        terminals = []
        for terminal in np.where(connectivity[center] != 0)[0]:
            terminal = int(terminal)
            terminal_charge = nuclear_charges[terminal]
            if terminal_charge < 7 or terminal_charge <= center_charge:
                continue
            if not _pi_capable_atom(terminal, atom_map, angular_map):
                continue
            terminals.append(terminal)

        if len(terminals) >= 2:
            resonance_atoms.update(terminals)

    return resonance_atoms


def _candidate_coefficient_map(candidate):
    """Return a zero-based coefficient map for an NBO candidate."""

    return {
        int(item['nao_index']) - 1: float(item['coefficient'])
        for item in candidate.get('coefficients', [])
    }


def _candidate_abs_overlap(candidate_a, candidate_b):
    """Return the absolute overlap of two stored candidate vectors."""

    coeffs_a = _candidate_coefficient_map(candidate_a)
    coeffs_b = _candidate_coefficient_map(candidate_b)
    if len(coeffs_a) > len(coeffs_b):
        coeffs_a, coeffs_b = coeffs_b, coeffs_a
    return abs(float(sum(
        coefficient * coeffs_b.get(index, 0.0)
        for index, coefficient in coeffs_a.items()
    )))


def _pi_lone_pair_conflict(pi_candidate,
                           lone_pair_candidate,
                           overlap_threshold=1.0e-5):
    """Return whether a one-center LP occupies the same orbital as a pi bond."""

    lone_pair_atoms = tuple(int(atom) for atom in
                            lone_pair_candidate.get('atoms', ()))
    if len(lone_pair_atoms) != 1:
        return False
    if lone_pair_atoms[0] not in set(int(atom) for atom in
                                     pi_candidate.get('atoms', ())):
        return False
    return _candidate_abs_overlap(pi_candidate, lone_pair_candidate) > overlap_threshold


def _requested_pi_pairs(constraints, natoms):
    """Return normalized user-requested pi pairs."""

    pairs = set()
    if constraints is None:
        return pairs
    for attr in ('required_pi_bonds', 'allowed_pi_bonds'):
        for pair in getattr(constraints, attr, ()):
            pairs.add(_normalize_atom_pair(pair, natoms))
    return pairs


def _allow_requested_non_sigma_pi_candidate(molecule):
    """Return whether a requested non-sigma pi pair may become active."""

    multiplicity = int(round(float(molecule.get_multiplicity())))
    return multiplicity == 1


def _allow_requested_pi_pair_candidate(molecule, pair, connectivity):
    """Return whether a requested pi pair may enter the real active pool."""

    atom_i, atom_j = tuple(int(atom) for atom in pair)
    if connectivity[atom_i, atom_j] != 0:
        return True
    return _allow_requested_non_sigma_pi_candidate(molecule)


def _pi_pair_candidates(atom_i,
                        atom_j,
                        density,
                        spin_density,
                        atom_map,
                        angular_map,
                        populations,
                        core_by_atom,
                        bond_min_atom_weight,
                        pi_min_occupation,
                        norb,
                        serial,
                        source,
                        non_sigma=True,
                        include_pair_candidate=True,
                        force_pair_candidate=False):
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
    occupied_roots = []
    used_complement_roots = set()

    if include_pair_candidate:
        for root in order:
            if accepted >= 1:
                break
            occupation = float(occupations[root])
            if occupation < pi_min_occupation and not force_pair_candidate:
                break
            local_vector = rotations[:, root]
            left_weight = float(np.sum(local_vector[:len(left)]**2))
            right_weight = float(np.sum(local_vector[len(left):]**2))
            if (left_weight < bond_min_atom_weight or
                    right_weight < bond_min_atom_weight):
                continue

            vector = _candidate_vector_from_local(pair_indices,
                                                  local_vector,
                                                  norb)
            candidate = {
                'index': serial,
                'type': 'BD',
                'subtype': 'pi',
                'atoms': (int(atom_i), int(atom_j)),
                'non_sigma': bool(non_sigma),
                'occupation': _candidate_occupation(vector, density),
                'polarization': {
                    int(atom_i + 1): left_weight,
                    int(atom_j + 1): right_weight,
                },
                'coefficients': _candidate_coefficients(vector),
                'source': source,
            }
            candidates.append(candidate)
            occupied_roots.append(int(root))
            serial += 1
            accepted += 1

    used_complement_roots.update(occupied_roots)
    for candidate in [cand for cand in candidates if cand.get('type') == 'BD']:
        serial = _append_antibonding_complement(
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
            f'{source} antibonding complement',
        )

    if spin_density is not None:
        pair_spin_density = spin_density[np.ix_(pair_indices, pair_indices)]
        pair_spin_density = 0.5 * (pair_spin_density + pair_spin_density.T)
        spin_occupations, spin_rotations = np.linalg.eigh(pair_spin_density)
        for root in np.argsort(spin_occupations)[::-1]:
            spin_occupation = float(spin_occupations[root])
            if spin_occupation < _SPIN_ACTIVE_MIN_OCCUPATION:
                break

            local_vector = spin_rotations[:, root]
            left_weight = float(np.sum(local_vector[:len(left)]**2))
            right_weight = float(np.sum(local_vector[len(left):]**2))
            if (left_weight < bond_min_atom_weight or
                    right_weight < bond_min_atom_weight):
                continue

            vector = _candidate_vector_from_local(pair_indices,
                                                  local_vector,
                                                  norb)
            candidates.append({
                'index': serial,
                'type': 'BD',
                'subtype': 'pi',
                'atoms': (int(atom_i), int(atom_j)),
                'non_sigma': bool(non_sigma),
                'electron_count': 1.0,
                'occupation': _candidate_occupation(vector, density),
                'spin_occupation': float(vector.T @ spin_density @ vector),
                'polarization': {
                    int(atom_i + 1): left_weight,
                    int(atom_j + 1): right_weight,
                },
                'coefficients': _candidate_coefficients(vector),
                'source': f'{source} spin-active one-electron pi candidate',
            })
            serial += 1
            break

    return candidates, serial
def _candidate_type_counts(candidates):
    """Count NBO candidate types."""

    counts = {}
    for candidate in candidates:
        candidate_type = candidate.get('type', '?')
        counts[candidate_type] = counts.get(candidate_type, 0) + 1
    return counts


def _candidate_electron_count(candidate):
    """Return the ideal Lewis electron count for a candidate."""

    return float(candidate.get('electron_count', 2.0))


def _split_sigma_pi_nbos(nbo_list):
    """Split selected NBOs into fixed sigma framework and pi part."""

    pi_nbos = [candidate for candidate in nbo_list
               if candidate.get('type') == 'BD' and
               candidate.get('subtype') == 'pi']
    sigma_nbos = [candidate for candidate in nbo_list
                  if candidate not in pi_nbos]
    sigma_occupation = float(sum(candidate.get('occupation', 0.0)
                                 for candidate in sigma_nbos))
    pi_occupation = float(sum(candidate.get('occupation', 0.0)
                              for candidate in pi_nbos))
    return {
        'sigma_nbo_list': sigma_nbos,
        'pi_nbo_list': pi_nbos,
        'sigma_counts': _candidate_type_counts(sigma_nbos),
        'pi_counts': _candidate_type_counts(pi_nbos),
        'sigma_electron_pairs': len(sigma_nbos),
        'pi_electron_pairs': len(pi_nbos),
        'sigma_occupation_sum': sigma_occupation,
        'pi_occupation_sum': pi_occupation,
    }


def _neutral_valence_electron_count(nuclear_charge):
    """Return the neutral-atom valence electron count used for Lewis checks."""

    n_core = _n_core_orbitals_from_z(nuclear_charge)
    return max(0.0, float(nuclear_charge) - 2.0 * float(n_core))


def _lewis_valence_target(nuclear_charge):
    """Return the duet/octet target for Lewis diagnostics."""

    if nuclear_charge <= 2:
        return 2.0
    return 8.0


def _build_lewis_accounting(molecule, nbo_list):
    """Compute per-atom Lewis electron ownership diagnostics."""

    labels = molecule.get_labels()
    natoms = molecule.number_of_atoms()
    nuclear_charges = np.array([
        chemical_element_identifier(labels[atom]) for atom in range(natoms)
    ], dtype=float)
    neutral_valence = np.array([
        _neutral_valence_electron_count(charge) for charge in nuclear_charges
    ], dtype=float)
    valence_targets = np.array([
        _lewis_valence_target(charge) for charge in nuclear_charges
    ], dtype=float)

    atom_electron_count = np.zeros(natoms, dtype=float)
    atom_core_electron_count = np.zeros(natoms, dtype=float)
    selected_electron_count = 0.0
    warnings = []

    for candidate in nbo_list:
        atoms = tuple(int(atom) for atom in candidate.get('atoms', ()))
        if not atoms:
            continue

        contribution = _candidate_electron_count(candidate)
        selected_electron_count += contribution
        share = contribution / float(len(atoms))
        for atom in atoms:
            atom_electron_count[atom] += share
            if candidate.get('type') == 'CR':
                atom_core_electron_count[atom] += share

    atom_valence_electron_count = atom_electron_count - atom_core_electron_count
    formal_charges = neutral_valence - atom_valence_electron_count
    octet_excess = np.maximum(atom_valence_electron_count - valence_targets, 0.0)
    electron_ownership_error = abs(
        float(np.sum(atom_electron_count)) - selected_electron_count)

    for atom, excess in enumerate(octet_excess):
        if excess > 1.0e-8:
            warnings.append(
                f'{labels[atom]}{atom + 1} exceeds the Lewis valence target by '
                f'{excess:.6f} electrons.'
            )

    return {
        'selected_lewis_electron_count': float(selected_electron_count),
        'atom_electron_count': [float(value) for value in atom_electron_count],
        'atom_core_electron_count': [
            float(value) for value in atom_core_electron_count
        ],
        'atom_valence_electron_count': [
            float(value) for value in atom_valence_electron_count
        ],
        'neutral_valence_electron_count': [
            float(value) for value in neutral_valence
        ],
        'formal_charges': [float(value) for value in formal_charges],
        'octet_excess': [float(value) for value in octet_excess],
        'electron_ownership_error': float(electron_ownership_error),
        'valence_warnings': warnings,
    }


def _normalize_atom_index(index, natoms):
    """Normalize a user atom index to zero-based form."""

    atom = int(index)
    if 1 <= atom <= natoms:
        return atom - 1
    assert_msg_critical(0 <= atom < natoms,
                        f'NBO constraints: invalid atom index {index}')
    return atom


def _constraint_formal_charge_targets(molecule, constraints):
    """Return target formal charges for Lewis scoring, or None for no target."""

    natoms = molecule.number_of_atoms()
    targets = np.zeros(natoms, dtype=float)
    formal_charge_constraints = dict(constraints.formal_charges)

    if formal_charge_constraints:
        for atom_index, charge in formal_charge_constraints.items():
            atom = _normalize_atom_index(atom_index, natoms)
            targets[atom] = float(charge)
        return targets

    if abs(float(molecule.get_charge())) < 1.0e-12:
        return targets

    return None


def _build_lewis_score_terms(molecule,
                             assignment,
                             constraints,
                             allowed_pi_bonus=0.0):
    """Build transparent score terms for a Lewis assignment."""

    formal_charges = np.array(assignment.get('formal_charges', []), dtype=float)
    octet_excess = np.array(assignment.get('octet_excess', []), dtype=float)
    target_electrons = float(assignment.get('target_electron_count',
                                            2.0 * assignment.get('target_electron_pairs', 0)))
    selected_electrons = float(assignment.get('selected_lewis_electron_count',
                                              2.0 * assignment.get('electron_pairs', 0)))
    target_formal_charges = _constraint_formal_charge_targets(molecule,
                                                              constraints)

    occupation_reward = float(assignment.get('occupation_sum', 0.0))
    electron_pair_penalty = (_ELECTRON_PAIR_PENALTY_WEIGHT *
                             abs(selected_electrons - target_electrons))

    if target_formal_charges is None:
        formal_charge_error = 0.0
        formal_charge_penalty = 0.0
        target_formal_charge_list = None
    else:
        formal_charge_error = float(
            np.sum((formal_charges - target_formal_charges)**2))
        formal_charge_penalty = (_FORMAL_CHARGE_PENALTY_WEIGHT *
                                 formal_charge_error)
        target_formal_charge_list = [
            float(charge) for charge in target_formal_charges
        ]

    octet_error = float(np.sum(octet_excess**2))
    octet_penalty = _OCTET_PENALTY_WEIGHT * octet_error
    valence_penalty = (_VALENCE_WARNING_PENALTY_WEIGHT *
                       len(assignment.get('valence_warnings', [])))
    nonclassical_count = sum(
        1 for candidate in assignment.get('nbo_list', [])
        if candidate.get('non_sigma', False)
    )
    nonclassical_penalty = (_NONCLASSICAL_PENALTY_WEIGHT *
                            float(nonclassical_count))

    total = (occupation_reward + float(allowed_pi_bonus) -
             electron_pair_penalty - formal_charge_penalty - octet_penalty -
             valence_penalty - nonclassical_penalty)

    return {
        'occupation_reward': occupation_reward,
        'allowed_pi_bonus': float(allowed_pi_bonus),
        'electron_pair_penalty': float(electron_pair_penalty),
        'formal_charge_penalty': float(formal_charge_penalty),
        'formal_charge_error': float(formal_charge_error),
        'target_formal_charges': target_formal_charge_list,
        'octet_penalty': float(octet_penalty),
        'octet_error': float(octet_error),
        'valence_penalty': float(valence_penalty),
        'nonclassical_penalty': float(nonclassical_penalty),
        'nonclassical_count': int(nonclassical_count),
        'total': float(total),
    }


def _active_space_fields(nbo_list, active_candidate_ids=None):
    """Return fixed and variable active-space fields for an assignment."""

    if active_candidate_ids is None:
        active_candidate_ids = {
            int(candidate.get('index', 0))
            for candidate in nbo_list
            if candidate.get('type') == 'BD' and
            candidate.get('subtype') == 'pi'
        }
    else:
        active_candidate_ids = {int(index) for index in active_candidate_ids}

    active_nbos = [candidate for candidate in nbo_list
                   if int(candidate.get('index', 0)) in active_candidate_ids]
    fixed_nbos = [candidate for candidate in nbo_list
                  if int(candidate.get('index', 0)) not in active_candidate_ids]
    active_pi_nbos = [candidate for candidate in active_nbos
                      if candidate.get('type') == 'BD' and
                      candidate.get('subtype') == 'pi']
    active_lone_pair_nbos = [candidate for candidate in active_nbos
                             if candidate.get('type') == 'LP' and
                             _candidate_electron_count(candidate) > 1.0]
    active_one_electron_nbos = [candidate for candidate in active_nbos
                                if _candidate_electron_count(candidate) == 1.0]

    return {
        'fixed_nbo_list': fixed_nbos,
        'active_nbo_list': active_nbos,
        'active_pi_nbo_list': active_pi_nbos,
        'active_lone_pair_nbo_list': active_lone_pair_nbos,
        'active_one_electron_nbo_list': active_one_electron_nbos,
        'active_counts': _candidate_type_counts(active_nbos),
        'fixed_counts': _candidate_type_counts(fixed_nbos),
        'active_electron_pairs': len(active_nbos),
        'active_pi_electron_pairs': len(active_pi_nbos),
        'active_lone_pair_electron_pairs': len(active_lone_pair_nbos),
        'active_one_electron_count': len(active_one_electron_nbos),
        'active_occupation_sum': float(sum(candidate.get('occupation', 0.0)
                                           for candidate in active_nbos)),
        'fixed_occupation_sum': float(sum(candidate.get('occupation', 0.0)
                                          for candidate in fixed_nbos)),
        'active_electron_count': _candidate_list_electron_count(active_nbos),
        'fixed_electron_count': _candidate_list_electron_count(fixed_nbos),
    }


def _assignment_payload(molecule,
                        nbo_list,
                        target_pairs,
                        occupation_sum,
                        warnings,
                        constraints,
                        allowed_pi_bonus=0.0,
                        active_candidate_ids=None):
    """Build a Lewis assignment dictionary with explicit sigma/pi parts."""

    payload = {
        'nbo_list': nbo_list,
        'counts': _candidate_type_counts(nbo_list),
        'electron_pairs': len(nbo_list),
        'target_electron_pairs': target_pairs,
        'target_electron_count': float(molecule.number_of_electrons()),
        'occupation_sum': occupation_sum,
        'warnings': warnings,
    }
    payload.update(_split_sigma_pi_nbos(nbo_list))
    payload.update(_active_space_fields(nbo_list, active_candidate_ids))
    payload.update(_build_lewis_accounting(molecule, nbo_list))
    score_terms = _build_lewis_score_terms(molecule,
                                           payload,
                                           constraints,
                                           allowed_pi_bonus=allowed_pi_bonus)
    payload['score_terms'] = score_terms
    payload['score'] = score_terms['total']
    return payload


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

    type_priority = {'CR': 0, 'LP': 1, 'BD': 2, 'SOMO': 3}
    subtype_priority = {'sigma': 0, 'pi': 1}
    return (
        type_priority.get(candidate.get('type', '?'), 9),
        subtype_priority.get(candidate.get('subtype', ''), 9),
        -float(candidate.get('occupation', 0.0)),
        tuple(candidate.get('atoms', ())),
        int(candidate.get('index', 0)),
    )


def _build_primary_assignment(molecule, candidates, constraints):
    """Select the primary Lewis/NBO set from generated candidates.

    The assignment layer selects occupied candidates up to the electron count,
    respecting forbidden bonds and prioritizing required bonds. Resonance
    alternatives and candidate-only acceptors are handled by later layers.
    """

    natoms = molecule.number_of_atoms()
    target_electrons = float(molecule.number_of_electrons())
    target_pairs = int(np.ceil(0.5 * target_electrons))
    required_pairs, forbidden_pairs = _constraint_pair_sets(constraints, natoms)
    required_pi_pairs, allowed_pi_pairs, forbidden_pi_pairs = _constraint_pi_pair_sets(
        constraints,
        natoms,
    )
    connectivity = molecule.get_connectivity_matrix()
    closed_shell_anion = (
        int(round(float(molecule.get_multiplicity()))) == 1 and
        float(molecule.get_charge()) < -1.0e-12
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
        if (candidate.get('subtype') == 'pi' and
                atoms in (required_pi_pairs | allowed_pi_pairs) and
                not _allow_requested_pi_pair_candidate(molecule,
                                                       atoms,
                                                       connectivity)):
            return False
        return True

    pool = [candidate for candidate in candidates if allowed(candidate)]
    selected = []
    selected_ids = set()
    warnings = []

    def selected_electrons():
        return float(sum(_candidate_electron_count(candidate)
                         for candidate in selected))

    def can_add(candidate):
        trial = selected + [candidate]
        trial_pi = [cand for cand in trial
                    if cand.get('type') == 'BD' and
                    cand.get('subtype') == 'pi']
        return (
            selected_electrons() + _candidate_electron_count(candidate) <=
            target_electrons + 1.0e-12 and
            _compatible_pi_combo(trial_pi) and
            _compatible_pi_lone_pair_combo(trial) and
            _compatible_pi_one_electron_combo(trial)
        )

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
            if not can_add(candidate):
                continue
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
            if not can_add(candidate):
                continue
            selected.append(candidate)
            selected_ids.add(candidate.get('index'))

    remaining = [
        candidate for candidate in pool
        if candidate.get('index') not in selected_ids and
        candidate.get('type') in {'LP', 'SOMO', 'BD'} and
           not (candidate.get('type') == 'SOMO' and
               candidate.get('subtype') == 'resonance-one-electron') and
        not (candidate.get('type') == 'LP' and
             candidate.get('subtype') in {
                 'radical-lone-pair',
                 'resonance-lone-pair',
             })
    ]
    for candidate in sorted(remaining, key=_candidate_sort_key):
        if selected_electrons() >= target_electrons - 1.0e-12:
            break
        if not can_add(candidate):
            continue
        selected.append(candidate)
        selected_ids.add(candidate.get('index'))

    if not required_pi_pairs:
        selected_pi_ids = {
            _candidate_index(candidate)
            for candidate in selected
            if candidate.get('type') == 'BD' and candidate.get('subtype') == 'pi'
        }
        fixed_selected = [
            candidate for candidate in selected
            if _candidate_index(candidate) not in selected_pi_ids
        ]
        pi_pool = [
            candidate for candidate in pool
            if candidate.get('type') == 'BD' and
            candidate.get('subtype') == 'pi'
        ]
        best_pi_combo = _best_pi_matching(pi_pool,
                                          fixed_selected,
                                          target_electrons)
        if (_candidate_list_electron_count(fixed_selected + best_pi_combo) >
                selected_electrons() + 1.0e-12):
            selected = fixed_selected + best_pi_combo
            selected_ids = {_candidate_index(candidate) for candidate in selected}

    if (closed_shell_anion and
            selected_electrons() < target_electrons - 1.0e-12 and
            any(candidate.get('type') == 'BD' and
                candidate.get('subtype') == 'pi'
                for candidate in selected)):
        resonance_lone_pairs = [
            candidate for candidate in pool
            if candidate.get('index') not in selected_ids and
            candidate.get('type') == 'LP' and
            candidate.get('subtype') == 'resonance-lone-pair' and
            _candidate_electron_count(candidate) > 1.0
        ]
        for candidate in sorted(resonance_lone_pairs, key=_candidate_sort_key):
            if selected_electrons() >= target_electrons - 1.0e-12:
                break
            if not can_add(candidate):
                continue
            selected.append(candidate)
            selected_ids.add(candidate.get('index'))

    if selected_electrons() < target_electrons - 1.0e-12:
        warnings.append(
            f'Primary assignment selected {selected_electrons():.6f} electrons; '
            f'target is {target_electrons:.6f}.'
        )
    if forbidden_pairs:
        warnings.append('Forbidden bonds were excluded from primary assignment.')

    selected = sorted(selected, key=lambda item: int(item.get('index', 0)))
    occupation_sum = float(sum(candidate.get('occupation', 0.0)
                               for candidate in selected))

    return _assignment_payload(molecule,
                               selected,
                               target_pairs,
                               occupation_sum,
                               warnings,
                               constraints)


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


def _compatible_pi_one_electron_combo(candidates):
    """Return whether pi bonds and one-electron centers do not overlap."""

    pi_atoms = {
        int(atom)
        for candidate in candidates
        if candidate.get('type') == 'BD' and candidate.get('subtype') == 'pi'
        for atom in candidate.get('atoms', ())
    }
    for candidate in candidates:
        if _candidate_electron_count(candidate) != 1.0:
            continue
        atoms = tuple(int(atom) for atom in candidate.get('atoms', ()))
        if len(atoms) == 1 and atoms[0] in pi_atoms:
            return False
    return True


def _compatible_pi_lone_pair_combo(candidates):
    """Return whether pi bonds and one-center lone pairs do not overlap."""

    pi_candidates = [
        candidate for candidate in candidates
        if candidate.get('type') == 'BD' and candidate.get('subtype') == 'pi'
    ]
    for candidate in candidates:
        if candidate.get('type') != 'LP':
            continue
        if _candidate_electron_count(candidate) <= 1.0:
            continue
        if any(_pi_lone_pair_conflict(pi_candidate, candidate)
               for pi_candidate in pi_candidates):
            return False
    return True


def _best_pi_matching(pi_candidates,
                      fixed_candidates,
                      target_electrons,
                      max_combos=200000):
    """Return the best non-overlapping pi matching for a fixed Lewis set."""

    fixed_electrons = _candidate_list_electron_count(fixed_candidates)
    available_electrons = target_electrons - fixed_electrons
    if available_electrons < 1.0e-12:
        return []

    candidates = sorted(pi_candidates, key=_candidate_sort_key)
    best_combo = []
    best_key = None
    explored = 0

    def evaluate(combo):
        nonlocal best_combo, best_key

        combo_electrons = _candidate_list_electron_count(combo)
        if combo_electrons > available_electrons + 1.0e-12:
            return
        trial = list(fixed_candidates) + list(combo)
        if not _compatible_pi_one_electron_combo(trial):
            return
        if not _compatible_pi_lone_pair_combo(trial):
            return
        key = (
            abs(available_electrons - combo_electrons),
            -combo_electrons,
            -len(combo),
            -float(sum(candidate.get('occupation', 0.0)
                       for candidate in combo)),
            tuple(_candidate_index(candidate) for candidate in combo),
        )
        if best_key is None or key < best_key:
            best_key = key
            best_combo = list(combo)

    def search(start, combo, used_atoms):
        nonlocal explored

        if explored >= max_combos:
            return
        explored += 1
        evaluate(combo)
        for index in range(start, len(candidates)):
            candidate = candidates[index]
            atoms = tuple(int(atom) for atom in candidate.get('atoms', ()))
            if len(atoms) != 2:
                continue
            if any(atom in used_atoms for atom in atoms):
                continue
            next_combo = combo + [candidate]
            if (_candidate_list_electron_count(next_combo) >
                    available_electrons + 1.0e-12):
                continue
            search(index + 1, next_combo, used_atoms | set(atoms))

    search(0, [], set())
    return best_combo


def _enumerate_pi_matchings(pi_candidates,
                            target_electrons,
                            max_matchings=64,
                            required_pi_pairs=None):
    """Return high-scoring non-overlapping pi matchings for alternatives."""

    required_pi_pairs = set(required_pi_pairs or [])
    candidates = sorted(pi_candidates, key=_candidate_sort_key)
    target_electrons = float(target_electrons)
    matchings = []
    seen_signatures = set()

    def add_matching(combo):
        combo_pairs = _active_pi_bonds(combo)
        if not required_pi_pairs.issubset(combo_pairs):
            return
        signature = tuple(sorted(combo_pairs))
        if signature in seen_signatures:
            return
        seen_signatures.add(signature)
        key = (
            -float(sum(candidate.get('occupation', 0.0)
                       for candidate in combo)),
            tuple(_candidate_index(candidate) for candidate in combo),
        )
        matchings.append((key, tuple(combo)))
        matchings.sort(key=lambda item: item[0])
        del matchings[max(1, int(max_matchings)):]

    def search(start, combo, used_atoms, electron_count):
        if electron_count > target_electrons + 1.0e-12:
            return
        if abs(electron_count - target_electrons) <= 1.0e-12:
            add_matching(combo)
            return

        for index in range(start, len(candidates)):
            candidate = candidates[index]
            atoms = tuple(int(atom) for atom in candidate.get('atoms', ()))
            if len(atoms) != 2:
                continue
            if any(atom in used_atoms for atom in atoms):
                continue
            candidate_electrons = _candidate_electron_count(candidate)
            search(index + 1,
                   combo + [candidate],
                   used_atoms | set(atoms),
                   electron_count + candidate_electrons)

    search(0, [], set(), 0.0)
    return [combo for _, combo in matchings]


def _polar_pi_resonance_units(molecule, pi_candidates, lone_pair_candidates):
    """Return generic polar units with interchangeable terminal pi bonds."""

    labels = molecule.get_labels()
    connectivity = molecule.get_connectivity_matrix()
    natoms = molecule.number_of_atoms()
    nuclear_charges = [
        chemical_element_identifier(labels[atom]) for atom in range(natoms)
    ]
    pi_by_pair = {
        tuple(sorted(candidate.get('atoms', ()))): candidate
        for candidate in pi_candidates
        if candidate.get('type') == 'BD' and candidate.get('subtype') == 'pi'
    }
    lone_pairs_by_atom = {}
    seen_lone_pair_ids = set()
    for candidate in lone_pair_candidates:
        if candidate.get('type') != 'LP':
            continue
        if _candidate_electron_count(candidate) <= 1.0:
            continue
        candidate_id = _candidate_index(candidate)
        if candidate_id in seen_lone_pair_ids:
            continue
        seen_lone_pair_ids.add(candidate_id)
        atoms = tuple(int(atom) for atom in candidate.get('atoms', ()))
        if len(atoms) != 1:
            continue
        lone_pairs_by_atom.setdefault(atoms[0], []).append(candidate)

    units = []
    for center in range(natoms):
        center_charge = nuclear_charges[center]
        if center_charge < 5:
            continue

        terminals = []
        pi_by_terminal = {}
        lone_pairs_by_terminal = {}
        for terminal in np.where(connectivity[center] != 0)[0]:
            terminal = int(terminal)
            terminal_charge = nuclear_charges[terminal]
            if terminal_charge < 7 or terminal_charge <= center_charge:
                continue
            pair = tuple(sorted((center, terminal)))
            if pair not in pi_by_pair:
                continue
            terminal_lone_pairs = sorted(
                lone_pairs_by_atom.get(terminal, []),
                key=_candidate_sort_key,
            )
            if not terminal_lone_pairs:
                continue
            terminals.append(terminal)
            pi_by_terminal[terminal] = pi_by_pair[pair]
            lone_pairs_by_terminal[terminal] = terminal_lone_pairs

        if len(terminals) < 2:
            continue
        units.append({
            'center': int(center),
            'terminals': tuple(sorted(terminals)),
            'pi_by_terminal': pi_by_terminal,
            'lone_pairs_by_terminal': lone_pairs_by_terminal,
        })

    return units


def _polar_resonance_unit_state(pi_combo, unit):
    """Return the selected terminal pi state for a polar resonance unit."""

    pi_pairs = _active_pi_bonds(pi_combo)
    selected_terminals = [
        terminal for terminal in unit['terminals']
        if tuple(sorted((unit['center'], terminal))) in pi_pairs
    ]
    if len(selected_terminals) != 1:
        return None
    return int(selected_terminals[0])


def _polar_lone_pair_candidates_for_pi_combo(pi_combo, polar_units):
    """Return coupled terminal LPs for a polar pi combo, or None if invalid."""

    lone_pairs_by_id = {}
    for unit in polar_units:
        selected_terminal = _polar_resonance_unit_state(pi_combo, unit)
        if selected_terminal is None:
            return None
        for terminal in unit['terminals']:
            terminal_lone_pairs = list(unit['lone_pairs_by_terminal'][terminal])
            if terminal == selected_terminal:
                terminal_lone_pairs = sorted(
                    terminal_lone_pairs,
                    key=lambda lone_pair: (
                        any(_pi_lone_pair_conflict(pi_candidate, lone_pair)
                            for pi_candidate in pi_combo),
                        -float(lone_pair.get('occupation', 0.0)),
                        _candidate_index(lone_pair),
                    ),
                )[:2]
            else:
                terminal_lone_pairs = sorted(
                    terminal_lone_pairs,
                    key=_candidate_sort_key,
                )[:3]
            if len(terminal_lone_pairs) < (2 if terminal == selected_terminal else 3):
                return None
            for lone_pair in terminal_lone_pairs:
                lone_pairs_by_id[_candidate_index(lone_pair)] = lone_pair

    return sorted(lone_pairs_by_id.values(), key=_candidate_sort_key)


def _polar_resonance_atoms(pi_combo, polar_units):
    """Return positive centers and negative terminals for polar units."""

    positive_atoms = set()
    negative_atoms = set()
    for unit in polar_units:
        selected_terminal = _polar_resonance_unit_state(pi_combo, unit)
        if selected_terminal is None:
            continue
        positive_atoms.add(int(unit['center']))
        negative_atoms.update(
            int(terminal) for terminal in unit['terminals']
            if int(terminal) != selected_terminal
        )
    return positive_atoms, negative_atoms


def _enumerate_polar_pi_lone_pair_combos(pi_candidates,
                                         polar_units,
                                         target_electrons,
                                         max_matchings=64,
                                         required_pi_pairs=None):
    """Return coupled pi/LP combos for polar resonance alternatives."""

    required_pi_pairs = set(required_pi_pairs or [])
    candidates = sorted(pi_candidates, key=_candidate_sort_key)
    target_electrons = float(target_electrons)
    matchings = []
    seen_signatures = set()

    def add_matching(pi_combo):
        combo_pairs = _active_pi_bonds(pi_combo)
        if not required_pi_pairs.issubset(combo_pairs):
            return
        polar_lone_pairs = _polar_lone_pair_candidates_for_pi_combo(
            pi_combo, polar_units)
        if polar_lone_pairs is None:
            return
        combo = tuple(list(pi_combo) + list(polar_lone_pairs))
        if abs(_candidate_list_electron_count(combo) - target_electrons) > 1.0e-12:
            return
        positive_atoms, negative_atoms = _polar_resonance_atoms(
            pi_combo, polar_units)
        signature = (
            tuple(sorted(combo_pairs)),
            tuple(sorted(negative_atoms)),
            tuple(sorted(positive_atoms)),
        )
        if signature in seen_signatures:
            return
        seen_signatures.add(signature)
        key = (
            -float(sum(candidate.get('occupation', 0.0)
                       for candidate in combo)),
            tuple(_candidate_index(candidate) for candidate in combo),
        )
        matchings.append((key, combo))
        matchings.sort(key=lambda item: item[0])
        del matchings[max(1, int(max_matchings)):]

    def search(start, pi_combo, used_atoms, electron_count):
        if electron_count > target_electrons + 1.0e-12:
            return
        add_matching(tuple(pi_combo))

        for index in range(start, len(candidates)):
            candidate = candidates[index]
            atoms = tuple(int(atom) for atom in candidate.get('atoms', ()))
            if len(atoms) != 2:
                continue
            if any(atom in used_atoms for atom in atoms):
                continue
            candidate_electrons = _candidate_electron_count(candidate)
            search(index + 1,
                   pi_combo + [candidate],
                   used_atoms | set(atoms),
                   electron_count + candidate_electrons)

    search(0, [], set(), 0.0)
    return [combo for _, combo in matchings]


def _active_pi_bonds(active_candidates):
    """Return sorted zero-based pi-bond pairs from active candidates."""

    return {
        tuple(sorted(candidate.get('atoms', ())))
        for candidate in active_candidates
        if candidate.get('type') == 'BD' and candidate.get('subtype') == 'pi'
    }


def _active_one_electron_atoms(active_candidates):
    """Return sorted zero-based atoms carrying active one-electron centers."""

    return {
        int(candidate.get('atoms', ())[0])
        for candidate in active_candidates
        if _candidate_electron_count(candidate) == 1.0 and
        len(candidate.get('atoms', ())) == 1
    }


def _alternative_resonance_signature(alternative):
    """Return a visible resonance signature for deduplicating alternatives."""

    pi_bonds = tuple(sorted(
        tuple(sorted(int(atom) for atom in pair))
        for pair in alternative.get('pi_bonds', [])
    ))
    one_electron_atoms = tuple(sorted(
        int(atom) for atom in alternative.get('active_one_electron_atoms', [])
    ))
    lone_pair_atoms = tuple(sorted(
        int(atom) for atom in alternative.get('active_lone_pair_atoms', [])
    ))
    positive_atoms = tuple(sorted(
        int(atom) for atom in alternative.get('active_positive_atoms', [])
    ))
    return pi_bonds, one_electron_atoms, lone_pair_atoms, positive_atoms


def _alternative_pi_signature(alternative):
    """Return the zero-based pi-bond signature visible in a Lewis alternative."""

    return tuple(sorted(
        tuple(sorted(int(atom) - 1 for atom in pair))
        for pair in alternative.get('pi_bonds', [])
    ))


def _missing_requested_pi_alternatives(requested_pi_pairs,
                                      alternatives,
                                      target_pairs,
                                      active_atoms=None,
                                      expected_one_electron_count=0,
                                      include_lone_pair_resonance=False,
                                      include_positive_resonance=False):
    """Return zero-weight placeholders for requested pi patterns not found."""

    active_atoms = set(active_atoms or [])
    represented_pairs = {
        _alternative_pi_signature(alternative)
        for alternative in alternatives
    }
    missing = []
    for pair in sorted(requested_pi_pairs):
        signature = (pair,)
        if signature in represented_pairs:
            continue
        pi_bonds = [tuple(int(atom + 1) for atom in pair)]
        complement_atoms = sorted(active_atoms - set(pair))
        one_electron_atoms = []
        lone_pair_atoms = []
        positive_atoms = []
        if expected_one_electron_count > 0 and complement_atoms:
            one_electron_atoms = [int(complement_atoms[0] + 1)]
        elif include_lone_pair_resonance and complement_atoms:
            lone_pair_atoms = [int(complement_atoms[0] + 1)]
        elif include_positive_resonance and complement_atoms:
            positive_atoms = [int(complement_atoms[0] + 1)]
        missing.append({
            'rank': 0,
            'weight': 0.0,
            'score': 0.0,
            'electron_pairs': target_pairs,
            'target_electron_pairs': target_pairs,
            'sigma_electron_pairs': max(0, target_pairs - 1),
            'pi_electron_pairs': 1,
            'pi_bonds': pi_bonds,
            'active_one_electron_atoms': one_electron_atoms,
            'active_lone_pair_atoms': lone_pair_atoms,
            'active_positive_atoms': positive_atoms,
            'warnings': [
                f"Requested pi bond {pi_bonds[0][0]}-{pi_bonds[0][1]} "
                "has no compatible Lewis alternative."
            ],
        })
    return missing


def _unique_pi_candidates_by_pair(pi_candidates):
    """Return one deterministic pi candidate per atom pair."""

    best_by_pair = {}
    for candidate in sorted(pi_candidates, key=_candidate_sort_key):
        pair = tuple(sorted(candidate.get('atoms', ())))
        best_by_pair.setdefault(pair, candidate)
    return [best_by_pair[pair] for pair in sorted(best_by_pair)]


def _candidate_index(candidate):
    """Return the integer candidate index."""

    return int(candidate.get('index', 0))


def _candidate_list_electron_count(candidates):
    """Return the ideal electron count for a candidate list."""

    return float(sum(_candidate_electron_count(candidate)
                     for candidate in candidates))


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
    target_electrons = float(primary_assignment.get(
        'target_electron_count',
        molecule.number_of_electrons(),
    ))
    required_pi_pairs, allowed_pi_pairs, forbidden_pi_pairs = (
        _constraint_pi_pair_sets(constraints, natoms))
    connectivity = molecule.get_connectivity_matrix()

    primary_nbos = list(primary_assignment.get('nbo_list', ()))
    warnings = []
    sigma_framework_pairs = {
        tuple(sorted(candidate.get('atoms', ())))
        for candidate in primary_nbos
        if candidate.get('type') == 'BD' and
        candidate.get('subtype') == 'sigma'
    }

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
        if (pair in (required_pi_pairs | allowed_pi_pairs) and
                pair not in sigma_framework_pairs and
                not _allow_requested_non_sigma_pi_candidate(molecule)):
            continue
        pi_pool.append(candidate)
    pi_pool = _unique_pi_candidates_by_pair(pi_pool)

    pi_active_atoms = {
        int(atom)
        for candidate in pi_pool
        for atom in candidate.get('atoms', ())
    }
    donor_lone_pairs = [
        candidate for candidate in primary_nbos
        if candidate.get('type') == 'LP' and
        _candidate_electron_count(candidate) > 1.0 and
        int(candidate.get('atoms', (-1,))[0]) in pi_active_atoms
    ]
    candidate_lone_pairs = [
        candidate for candidate in candidates
        if candidate.get('type') == 'LP' and
        _candidate_electron_count(candidate) > 1.0 and
        len(candidate.get('atoms', ())) == 1 and
        int(candidate.get('atoms', (-1,))[0]) in pi_active_atoms
    ]
    polar_units = _polar_pi_resonance_units(
        molecule,
        pi_pool,
        donor_lone_pairs + candidate_lone_pairs,
    )
    polar_lone_pair_ids = {
        _candidate_index(lone_pair)
        for unit in polar_units
        for lone_pairs in unit['lone_pairs_by_terminal'].values()
        for lone_pair in lone_pairs
    }
    active_lone_pair_ids = set(polar_lone_pair_ids)
    if (int(round(float(molecule.get_multiplicity()))) == 1 and
            float(molecule.get_charge()) < -1.0e-12):
        active_lone_pair_ids.update(
            _candidate_index(candidate) for candidate in donor_lone_pairs
        )
    one_electron_primary_ids = {
        int(candidate.get('index', 0)) for candidate in primary_nbos
        if _candidate_electron_count(candidate) == 1.0
    }

    fixed = [
        candidate for candidate in primary_nbos
        if not (candidate.get('type') == 'BD' and
                candidate.get('subtype') == 'pi') and
        int(candidate.get('index', 0)) not in active_lone_pair_ids and
        int(candidate.get('index', 0)) not in one_electron_primary_ids
    ]
    active_target_electrons = target_electrons - _candidate_list_electron_count(fixed)

    if active_target_electrons < -1.0e-12:
        warnings.append(
            'Primary assignment has more fixed electrons than the molecular target.'
        )
        active_target_electrons = 0.0

    one_electron_pool = [
        candidate for candidate in candidates
        if _candidate_electron_count(candidate) == 1.0 and
        candidate.get('type') in {'SOMO', 'BD'}
    ]
    expected_one_electron_count = max(
        int(round(float(molecule.get_multiplicity()) - 1.0)),
        int(round(active_target_electrons)) % 2,
    )
    include_lone_pair_resonance = (
        expected_one_electron_count == 0 and
        float(molecule.get_charge()) < -1.0e-12
    )
    include_positive_resonance = (
        expected_one_electron_count == 0 and
        float(molecule.get_charge()) > 1.0e-12
    )
    resonance_target_electrons = (
        4.0 if include_lone_pair_resonance else active_target_electrons
    )
    if include_lone_pair_resonance:
        active_lone_pair_ids.update(
            _candidate_index(candidate)
            for candidate in donor_lone_pairs + candidate_lone_pairs
        )

    active_pool_by_id = {}
    active_sources = pi_pool + one_electron_pool
    if include_lone_pair_resonance:
        active_sources += donor_lone_pairs + candidate_lone_pairs
    elif polar_units:
        active_sources += [
            candidate for candidate in donor_lone_pairs + candidate_lone_pairs
            if _candidate_index(candidate) in polar_lone_pair_ids
        ]
    for candidate in active_sources:
        active_pool_by_id[_candidate_index(candidate)] = candidate
    active_pool = sorted(active_pool_by_id.values(), key=_candidate_sort_key)

    if (_candidate_list_electron_count(active_pool) <
            resonance_target_electrons - 1.0e-12):
        warnings.append(
            f'Only {_candidate_list_electron_count(active_pool):.6f} active '
            f'electrons available for {resonance_target_electrons:.6f} active electrons.'
        )

    polar_resonance = bool(polar_units) and expected_one_electron_count == 0
    pi_only_resonance = all(
        candidate.get('type') == 'BD' and candidate.get('subtype') == 'pi'
        for candidate in active_pool
    )
    if polar_resonance:
        combo_iterable = _enumerate_polar_pi_lone_pair_combos(
            pi_pool,
            polar_units,
            resonance_target_electrons,
            max_matchings=max(32 * int(max_alternatives),
                              4 * int(max_alternatives),
                              int(max_alternatives),
                              64),
            required_pi_pairs=required_pi_pairs,
        )
    elif pi_only_resonance:
        combo_iterable = _enumerate_pi_matchings(
            active_pool,
            resonance_target_electrons,
            max_matchings=max(4 * int(max_alternatives), int(max_alternatives), 16),
            required_pi_pairs=required_pi_pairs,
        )
    else:
        combo_iterable = (
            combo
            for combo_size in range(len(active_pool) + 1)
            for combo in combinations(active_pool, combo_size)
        )

    alternatives = []
    for combo in combo_iterable:
        combo_electrons = _candidate_list_electron_count(combo)
        if abs(combo_electrons - resonance_target_electrons) > 1.0e-12:
            continue
        one_electron_count = sum(
            1 for candidate in combo
            if _candidate_electron_count(candidate) == 1.0
        )
        if one_electron_count != expected_one_electron_count:
            continue

        combo_ids = {_candidate_index(candidate) for candidate in combo}
        duplicate_one_center = False
        one_center_vectors = {}
        for candidate in combo:
            if _candidate_electron_count(candidate) != 1.0:
                continue
            if len(candidate.get('atoms', ())) != 1:
                continue
            coefficient_key = tuple(
                (int(item['nao_index']), round(float(item['coefficient']), 12))
                for item in candidate.get('coefficients', [])
            )
            if coefficient_key in one_center_vectors:
                duplicate_one_center = True
                break
            one_center_vectors[coefficient_key] = _candidate_index(candidate)
        if duplicate_one_center:
            continue

        combo_pairs = _active_pi_bonds(combo)
        if not required_pi_pairs.issubset(combo_pairs):
            continue
        occupied_active_atoms = {
            int(atom)
            for pair in combo_pairs
            for atom in pair
        }
        positive_atoms = []
        polar_positive_atoms = set()
        polar_lone_pair_atoms = set()
        if polar_resonance:
            pi_combo_for_polar = [
                candidate for candidate in combo
                if candidate.get('type') == 'BD' and
                candidate.get('subtype') == 'pi'
            ]
            polar_positive_atoms, polar_lone_pair_atoms = _polar_resonance_atoms(
                pi_combo_for_polar,
                polar_units,
            )
        if include_positive_resonance:
            positive_atoms = [
                int(atom + 1)
                for atom in sorted(pi_active_atoms - occupied_active_atoms)
            ]
        if polar_positive_atoms:
            positive_atoms = sorted(
                set(positive_atoms) |
                {int(atom + 1) for atom in polar_positive_atoms}
            )
        pi_combo = [candidate for candidate in combo
                    if candidate.get('type') == 'BD' and
                    candidate.get('subtype') == 'pi']
        lone_pair_combo = [candidate for candidate in combo
                           if candidate.get('type') == 'LP' and
                           _candidate_electron_count(candidate) > 1.0]
        if include_lone_pair_resonance:
            if len(pi_combo) != 1 or len(lone_pair_combo) != 1:
                continue
        if not _compatible_pi_combo(pi_combo):
            continue
        if not _compatible_pi_one_electron_combo(combo):
            continue
        if not polar_resonance and not _compatible_pi_lone_pair_combo(combo):
            continue

        nbo_list = sorted(fixed + list(combo),
                          key=lambda candidate: int(candidate.get('index', 0)))
        occupation_sum = float(sum(candidate.get('occupation', 0.0)
                                   for candidate in nbo_list))
        allowed_bonus = 0.05 * len(combo_pairs & allowed_pi_pairs)
        alternative = _assignment_payload(molecule,
                                          nbo_list,
                                          target_pairs,
                                          occupation_sum,
                                          list(warnings),
                                          constraints,
                                          allowed_pi_bonus=allowed_bonus,
                                          active_candidate_ids=combo_ids)
        alternative.update({
            'rank': 0,
            'weight': 0.0,
            'pi_bonds': [tuple(int(atom + 1) for atom in pair)
                         for pair in sorted(combo_pairs)],
            'active_one_electron_atoms': [
                int(atom + 1)
                for atom in sorted(_active_one_electron_atoms(combo))
            ],
            'active_lone_pair_atoms': (
                [int(atom + 1) for atom in sorted(polar_lone_pair_atoms)]
                if polar_lone_pair_atoms else
                [
                    int(candidate.get('atoms', (0,))[0] + 1)
                    for candidate in alternative.get('active_lone_pair_nbo_list', [])
                ]
            ),
            'active_positive_atoms': positive_atoms,
        })
        alternatives.append(alternative)

    if include_lone_pair_resonance:
        existing_signatures = {
            _alternative_resonance_signature(alternative)
            for alternative in alternatives
        }
        lone_pair_pool = sorted(
            {
                _candidate_index(candidate): candidate
                for candidate in donor_lone_pairs + candidate_lone_pairs
            }.values(),
            key=_candidate_sort_key,
        )
        direct_pi_pairs = (required_pi_pairs | allowed_pi_pairs)
        if not direct_pi_pairs and float(molecule.get_charge()) < -1.0e-12:
            direct_pi_pairs = {
                tuple(sorted(candidate.get('atoms', ())))
                for candidate in pi_pool
                if not candidate.get('non_sigma', False)
            }
        for pair in sorted(direct_pi_pairs):
            pair_pi_candidates = sorted(
                (
                    candidate for candidate in pi_pool
                    if tuple(sorted(candidate.get('atoms', ()))) == pair
                ),
                key=_candidate_sort_key,
            )
            if not pair_pi_candidates:
                continue

            pair_lone_pairs = [
                candidate for candidate in lone_pair_pool
                if len(candidate.get('atoms', ())) == 1 and
                int(candidate.get('atoms', (0,))[0]) not in pair
            ]
            if not pair_lone_pairs:
                continue

            combo = (pair_pi_candidates[0], pair_lone_pairs[0])
            if (abs(_candidate_list_electron_count(combo) -
                    resonance_target_electrons) > 1.0e-12):
                continue
            if not _compatible_pi_lone_pair_combo(combo):
                continue

            combo_ids = {_candidate_index(candidate) for candidate in combo}
            nbo_list = sorted(fixed + list(combo),
                              key=lambda candidate: int(candidate.get('index', 0)))
            occupation_sum = float(sum(candidate.get('occupation', 0.0)
                                       for candidate in nbo_list))
            allowed_bonus = 0.05 * len({pair} & allowed_pi_pairs)
            alternative = _assignment_payload(molecule,
                                              nbo_list,
                                              target_pairs,
                                              occupation_sum,
                                              list(warnings),
                                              constraints,
                                              allowed_pi_bonus=allowed_bonus,
                                              active_candidate_ids=combo_ids)
            alternative.update({
                'rank': 0,
                'weight': 0.0,
                'pi_bonds': [tuple(int(atom + 1) for atom in pair)],
                'active_one_electron_atoms': [],
                'active_lone_pair_atoms': [
                    int(pair_lone_pairs[0].get('atoms', (0,))[0] + 1)
                ],
                'active_positive_atoms': [],
            })
            signature = _alternative_resonance_signature(alternative)
            if signature in existing_signatures:
                continue
            existing_signatures.add(signature)
            alternatives.append(alternative)

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
                                           _alternative_resonance_signature(alt)))
    unique_alternatives = []
    seen_signatures = set()
    for alternative in alternatives:
        signature = _alternative_pi_signature(alternative)
        if signature in seen_signatures:
            continue
        seen_signatures.add(signature)
        unique_alternatives.append(alternative)
    alternatives = unique_alternatives
    if polar_resonance:
        polar_groups = {}
        for alternative in alternatives:
            polar_signature = (
                tuple(sorted(alternative.get('active_lone_pair_atoms', []))),
                tuple(sorted(alternative.get('active_positive_atoms', []))),
            )
            polar_groups.setdefault(polar_signature, []).append(alternative)
        for group in polar_groups.values():
            group.sort(key=lambda alt: (-alt['score'],
                                        _alternative_resonance_signature(alt)))

        diverse_alternatives = []
        while len(diverse_alternatives) < max(1, int(max_alternatives)):
            added = False
            for signature in sorted(polar_groups):
                group = polar_groups[signature]
                if not group:
                    continue
                diverse_alternatives.append(group.pop(0))
                added = True
                if len(diverse_alternatives) >= max(1, int(max_alternatives)):
                    break
            if not added:
                break
        alternatives = diverse_alternatives
    requested_signatures = set()
    if (resonance_target_electrons <= 2.0 + 1.0e-12 or
            expected_one_electron_count > 0 or include_lone_pair_resonance):
        requested_signatures = {
            (pair,) for pair in sorted(required_pi_pairs | allowed_pi_pairs)
        }
    selected_alternatives = []
    selected_signatures = set()
    if requested_signatures:
        for alternative in alternatives:
            signature = _alternative_pi_signature(alternative)
            if signature not in requested_signatures:
                continue
            if signature in selected_signatures:
                continue
            selected_alternatives.append(alternative)
            selected_signatures.add(signature)

    for alternative in alternatives:
        if len(selected_alternatives) >= max(1, int(max_alternatives)):
            break
        signature = _alternative_pi_signature(alternative)
        if signature in selected_signatures:
            continue
        selected_alternatives.append(alternative)
        selected_signatures.add(signature)
    alternatives = selected_alternatives[:max(1, int(max_alternatives))]
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

    for alternative, weight in zip(alternatives, weights):
        alternative['weight'] = float(weight)

    alternatives = sorted(
        alternatives,
        key=lambda alt: (-float(alt.get('weight', 0.0)),
                         -float(alt.get('score', 0.0)),
                         _alternative_resonance_signature(alt)),
    )

    for rank, alternative in enumerate(alternatives, start=1):
        alternative['rank'] = rank

    missing_requested = []
    if requested_signatures:
        missing_requested = _missing_requested_pi_alternatives(
            required_pi_pairs | allowed_pi_pairs,
            alternatives,
            target_pairs,
            active_atoms=pi_active_atoms,
            expected_one_electron_count=expected_one_electron_count,
            include_lone_pair_resonance=include_lone_pair_resonance,
            include_positive_resonance=include_positive_resonance,
        )
    for offset, alternative in enumerate(missing_requested,
                                         start=len(alternatives) + 1):
        alternative['rank'] = offset
    alternatives.extend(missing_requested)

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


def _candidate_atom_indices(candidate):
    """Return one-based atom indices for a candidate."""

    return tuple(int(atom + 1) for atom in candidate.get('atoms', ()))


def _build_donor_acceptor_diagnostics(candidates,
                                      primary_assignment,
                                      density,
                                      max_interactions=64,
                                      min_abs_coupling=0.0):
    """Build donor-acceptor density-coupling diagnostics.

    This is a diagnostic layer only. It uses the NAO density matrix as the
    current coupling operator and does not feed back into Lewis selection.
    """

    norb = density.shape[0]
    donors = [
        candidate for candidate in primary_assignment.get('nbo_list', [])
        if _candidate_electron_count(candidate) > 0.0
    ]
    acceptors = [
        candidate for candidate in candidates
        if candidate.get('type') in _ACCEPTOR_TYPES and
        _candidate_electron_count(candidate) == 0.0
    ]
    interactions = []

    for donor in donors:
        donor_vector = _candidate_vector_from_coefficients(donor, norb)
        if np.linalg.norm(donor_vector) < 1.0e-14:
            continue
        for acceptor in acceptors:
            acceptor_vector = _candidate_vector_from_coefficients(acceptor,
                                                                  norb)
            if np.linalg.norm(acceptor_vector) < 1.0e-14:
                continue

            coupling = float(donor_vector.T @ density @ acceptor_vector)
            overlap = float(np.dot(donor_vector, acceptor_vector))
            abs_coupling = abs(coupling)
            if abs_coupling < min_abs_coupling:
                continue

            interactions.append({
                'donor_index': int(donor.get('index', 0)),
                'donor_type': donor.get('type', ''),
                'donor_subtype': donor.get('subtype', ''),
                'donor_atoms': _candidate_atom_indices(donor),
                'donor_occupation': float(donor.get('occupation', 0.0)),
                'acceptor_index': int(acceptor.get('index', 0)),
                'acceptor_type': acceptor.get('type', ''),
                'acceptor_subtype': acceptor.get('subtype', ''),
                'acceptor_atoms': _candidate_atom_indices(acceptor),
                'acceptor_occupation': float(acceptor.get('occupation', 0.0)),
                'acceptor_parent_index': acceptor.get('parent_index'),
                'overlap': overlap,
                'abs_overlap': abs(overlap),
                'density_coupling': coupling,
                'abs_density_coupling': abs_coupling,
                'density_coupling_squared': float(coupling * coupling),
            })

    interactions = sorted(
        interactions,
        key=lambda item: (-item['abs_density_coupling'],
                          item['donor_index'],
                          item['acceptor_index']),
    )[:max(0, int(max_interactions))]

    return {
        'model': 'nao_density_coupling',
        'description': (
            'Density-coupling diagnostics from occupied primary donors to '
            'candidate-only BD*/RY acceptors; not a second-order perturbation energy.'
        ),
        'donor_count': len(donors),
        'acceptor_count': len(acceptors),
        'interaction_count': len(interactions),
        'interactions': interactions,
        'warnings': [
            'No NBO Fock matrix or energy denominator is used in this first diagnostic layer.'
        ],
    }


def _lewis_density_from_nbos(nbo_list, norb):
    """Build an ideal Lewis density from selected NBO vectors."""

    density = np.zeros((norb, norb))
    for candidate in nbo_list:
        vector = _candidate_vector_from_coefficients(candidate, norb)
        if np.linalg.norm(vector) < 1.0e-14:
            continue
        density += _candidate_electron_count(candidate) * np.outer(vector, vector)
    return 0.5 * (density + density.T)


def _lewis_spin_density_components_from_nbos(nbo_list, norb):
    """Build ideal total and spin Lewis densities from selected NBOs."""

    total_density = np.zeros((norb, norb))
    spin_density = np.zeros((norb, norb))
    for candidate in nbo_list:
        vector = _candidate_vector_from_coefficients(candidate, norb)
        if np.linalg.norm(vector) < 1.0e-14:
            continue
        projector = np.outer(vector, vector)
        electron_count = _candidate_electron_count(candidate)
        total_density += electron_count * projector
        if abs(electron_count - 1.0) < 1.0e-12:
            spin_density += projector

    return (0.5 * (total_density + total_density.T),
            0.5 * (spin_density + spin_density.T))


def _nra_fit_target_and_columns(nao_data,
                                spin_data,
                                alternatives,
                                indices,
                                norb,
                                spin_fit):
    """Build NRA vector target and structure columns for total/spin fitting."""

    spin_fit = str(spin_fit).lower()
    if spin_fit not in {'none', 'total_spin'}:
        spin_fit = 'total_spin'

    use_spin = spin_data is not None and spin_fit == 'total_spin'
    target_parts = [_vectorize_symmetric_block(nao_data.density, indices)]
    if use_spin:
        target_parts.append(_vectorize_symmetric_block(spin_data.spin_density,
                                                       indices))
    target = np.concatenate(target_parts)

    total_densities = []
    spin_densities = []
    columns = []
    for alternative in alternatives:
        total_density, spin_density = _lewis_spin_density_components_from_nbos(
            alternative.get('nbo_list', []), norb)
        total_densities.append(total_density)
        spin_densities.append(spin_density)
        column_parts = [_vectorize_symmetric_block(total_density, indices)]
        if use_spin:
            column_parts.append(_vectorize_symmetric_block(spin_density,
                                                           indices))
        columns.append(np.concatenate(column_parts))

    fit_matrix = np.column_stack(columns)
    return target, fit_matrix, total_densities, spin_densities, use_spin


def _nra_subspace_indices(subspace, alternatives, atom_map, angular_map, norb):
    """Return NAO indices used for the requested NRA/NRT fit."""

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


def _normalize_nra_prior_weights(prior_weights, alternatives):
    """Return normalized NRA prior weights aligned with alternatives."""

    count = len(alternatives)
    if prior_weights is None:
        return None, None, []

    warnings = []
    if isinstance(prior_weights, dict):
        raw = np.zeros(count, dtype=float)
        for position, alternative in enumerate(alternatives):
            rank = int(alternative.get('rank', position + 1))
            label = _nra_structure_label(alternative.get('pi_bonds', []))
            pi_key = tuple(tuple(pair) for pair in alternative.get('pi_bonds', []))
            for key in (rank, str(rank), label, pi_key, str(pi_key)):
                if key in prior_weights:
                    raw[position] = float(prior_weights[key])
                    break
    else:
        raw = np.array(prior_weights, dtype=float).reshape(-1)
        assert_msg_critical(
            raw.size == count,
            'NRA prior weights must have one entry per Lewis alternative.'
        )

    raw = np.maximum(raw, 0.0)
    total = float(np.sum(raw))
    if total < 1.0e-14:
        warnings.append(
            'NRA prior weights were all zero; prior regularization was disabled.'
        )
        return [float(value) for value in raw], None, warnings

    normalized = raw / total
    return [float(value) for value in raw], normalized, warnings


def _nra_weights_with_prior(matrix,
                            target,
                            prior_weights=None,
                            prior_strength=0.0,
                            prior_mode='regularized'):
    """Fit NRA simplex weights with optional prior regularization."""

    warnings = []
    mode = str(prior_mode).lower()
    if mode not in {'regularized', 'fixed'}:
        warnings.append(
            f"Unknown NRA prior mode '{prior_mode}'; regularized mode was used."
        )
        mode = 'regularized'

    strength = max(0.0, float(prior_strength))
    if prior_weights is None:
        weights = _equality_constrained_least_squares(matrix, target)
        return weights, {
            'active': False,
            'mode': mode,
            'strength': strength,
            'input_weights': None,
            'normalized_weights': None,
            'fixed': False,
            'warnings': warnings,
        }

    if mode == 'fixed':
        weights = np.array(prior_weights, dtype=float)
        warnings.append(
            'NRA prior mode fixed: density fitting was evaluated at normalized prior weights.'
        )
        return weights, {
            'active': True,
            'mode': mode,
            'strength': strength,
            'input_weights': None,
            'normalized_weights': [float(weight) for weight in weights],
            'fixed': True,
            'warnings': warnings,
        }

    if strength <= 0.0:
        weights = _equality_constrained_least_squares(matrix, target)
        return weights, {
            'active': False,
            'mode': mode,
            'strength': strength,
            'input_weights': None,
            'normalized_weights': [float(weight) for weight in prior_weights],
            'fixed': False,
            'warnings': warnings,
        }

    scale = np.sqrt(strength)
    augmented_matrix = np.vstack([
        matrix,
        scale * np.eye(matrix.shape[1]),
    ])
    augmented_target = np.concatenate([
        target,
        scale * np.array(prior_weights, dtype=float),
    ])
    weights = _equality_constrained_least_squares(augmented_matrix,
                                                  augmented_target)
    return weights, {
        'active': True,
        'mode': mode,
        'strength': strength,
        'input_weights': None,
        'normalized_weights': [float(weight) for weight in prior_weights],
        'fixed': False,
        'warnings': warnings,
    }


def _normalized_pi_bond_set(pi_bonds):
    """Return a stable one-based pi-bond tuple for structure signatures."""

    return tuple(sorted({tuple(sorted(int(atom) for atom in pair))
                         for pair in pi_bonds}))


def _nra_structure_label(pi_bonds):
    """Return a molecule-independent NRA/NRT structure signature."""

    pi_set = _normalized_pi_bond_set(pi_bonds)
    if not pi_set:
        return 'sigma-only'
    return 'pi:' + ','.join(f'{left}-{right}' for left, right in pi_set)


def _build_nra_results(molecule,
                       nao_data,
                       spin_data,
                       alternatives,
                       subspace,
                       fit_metric,
                       prior_weights=None,
                       prior_strength=0.0,
                       prior_mode='regularized',
                       spin_fit='total_spin'):
    """Fit ideal Lewis densities to the actual NAO density."""

    fit_metric = str(fit_metric).lower()
    warnings = []
    norb = nao_data.density.shape[0]

    if fit_metric != 'frobenius':
        warnings.append(
            f"Unknown NRA fit metric '{fit_metric}'; Frobenius metric was used."
        )
        fit_metric = 'frobenius'

    if not alternatives:
        return {
            'subspace': str(subspace).lower(),
            'fit_metric': fit_metric,
            'spin_fit': 'none',
            'weights': [],
            'residual_norm': None,
            'relative_residual': None,
            'spin_residual_norm': None,
            'relative_spin_residual': None,
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

    use_spin_fit = molecule.get_multiplicity() != 1 and spin_data is not None
    active_spin_fit = spin_fit if use_spin_fit else 'none'
    target, fit_matrix, lewis_densities, lewis_spin_densities, used_spin = (
        _nra_fit_target_and_columns(
            nao_data,
            spin_data if use_spin_fit else None,
            alternatives,
            indices,
            norb,
            active_spin_fit,
        )
    )
    if used_spin:
        warnings.append(
            'Open-shell NRA/NRT fit used concatenated total-density and spin-density residuals.'
        )
    raw_prior_weights, normalized_prior_weights, prior_warnings = (
        _normalize_nra_prior_weights(prior_weights, alternatives)
    )
    warnings.extend(prior_warnings)
    weights, prior = _nra_weights_with_prior(
        fit_matrix,
        target,
        prior_weights=normalized_prior_weights,
        prior_strength=prior_strength,
        prior_mode=prior_mode,
    )
    prior['input_weights'] = raw_prior_weights
    warnings.extend(prior.get('warnings', []))
    model = fit_matrix @ weights if len(weights) else np.zeros_like(target)
    residual = target - model
    residual_norm = float(np.linalg.norm(residual))
    target_norm = float(np.linalg.norm(target))
    relative_residual = (residual_norm / target_norm
                         if target_norm > 1.0e-14 else 0.0)
    total_target = _vectorize_symmetric_block(nao_data.density, indices)
    total_model = np.column_stack([
        _vectorize_symmetric_block(density, indices)
        for density in lewis_densities
    ]) @ weights if len(weights) else np.zeros_like(total_target)
    total_residual = total_target - total_model
    total_residual_norm = float(np.linalg.norm(total_residual))
    total_target_norm = float(np.linalg.norm(total_target))
    relative_total_residual = (total_residual_norm / total_target_norm
                               if total_target_norm > 1.0e-14 else 0.0)
    if used_spin:
        spin_target = _vectorize_symmetric_block(spin_data.spin_density,
                                                 indices)
        spin_model = np.column_stack([
            _vectorize_symmetric_block(density, indices)
            for density in lewis_spin_densities
        ]) @ weights if len(weights) else np.zeros_like(spin_target)
        spin_residual = spin_target - spin_model
        spin_residual_norm = float(np.linalg.norm(spin_residual))
        spin_target_norm = float(np.linalg.norm(spin_target))
        relative_spin_residual = (spin_residual_norm / spin_target_norm
                                  if spin_target_norm > 1.0e-14 else 0.0)
    else:
        spin_residual_norm = None
        relative_spin_residual = None

    structures = []
    for position, (alternative, weight, lewis_density, lewis_spin_density) in enumerate(zip(
            alternatives,
            weights,
            lewis_densities,
            lewis_spin_densities)):
        structure_vector = _vectorize_symmetric_block(lewis_density, indices)
        spin_structure_vector = _vectorize_symmetric_block(lewis_spin_density,
                                                           indices)
        structures.append({
            'rank': int(alternative.get('rank', 0)),
            'pi_bonds': list(alternative.get('pi_bonds', [])),
            'label': _nra_structure_label(alternative.get('pi_bonds', [])),
            'sigma_electron_pairs': int(alternative.get('sigma_electron_pairs', 0)),
            'pi_electron_pairs': int(alternative.get('pi_electron_pairs', 0)),
            'nra_weight': float(weight),
            'prior_weight': (None if normalized_prior_weights is None else
                             float(normalized_prior_weights[position])),
            'score_weight': float(alternative.get('weight', 0.0)),
            'score': float(alternative.get('score', 0.0)),
            'residual_norm': float(np.linalg.norm(total_target - structure_vector)),
            'spin_residual_norm': (None if not used_spin else
                                   float(np.linalg.norm(spin_target - spin_structure_vector))),
        })

    return {
        'subspace': str(subspace).lower(),
        'fit_metric': fit_metric,
        'spin_fit': 'total_spin' if used_spin else 'none',
        'weights': [float(weight) for weight in weights],
        'residual_norm': residual_norm,
        'relative_residual': relative_residual,
        'total_residual_norm': total_residual_norm,
        'relative_total_residual': relative_total_residual,
        'spin_residual_norm': spin_residual_norm,
        'relative_spin_residual': relative_spin_residual,
        'prior': prior,
        'structures': structures,
        'warnings': warnings,
    }


class NboDriver:
    """Natural Bond Orbital analysis driver."""

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
        self.nbo_report_level = 'full'
        self.nra_report_level = 'summary'
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
        """Compute NAO/NPA, NBO candidates, Lewis alternatives, and reports data."""

        assert_msg_critical(mode in {'npa', 'primary'},
                            "NboDriver: choose mode='npa' or mode='primary'")

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

            analyzer_options = OrbitalAnalyzerOptions(
                include_mo_analysis=compute_options.include_mo_analysis,
                include_nbo_candidates=compute_options.include_nbo_candidates,
                mo_analysis_top=compute_options.mo_analysis_top,
                mo_analysis_threshold=compute_options.mo_analysis_threshold,
                lone_pair_min_occupation=compute_options.lone_pair_min_occupation,
                rydberg_max_occupation=compute_options.rydberg_max_occupation,
                bond_min_occupation=compute_options.bond_min_occupation,
                bond_min_atom_weight=compute_options.bond_min_atom_weight,
                pi_min_occupation=compute_options.pi_min_occupation,
                conjugated_pi_max_path=compute_options.conjugated_pi_max_path,
            )
            orbital_analysis = OrbitalAnalyzer(
                molecule,
                basis,
                mol_orbs=mol_orbs,
                options=analyzer_options,
                constraints=compute_constraints,
            ).run()
            overlap = orbital_analysis.overlap
            density = orbital_analysis.density
            nao_data = orbital_analysis.nao_data
            spin_data = orbital_analysis.spin_data
            mo_analysis = orbital_analysis.mo_analysis
            nbo_candidates = orbital_analysis.orbital_candidates
            candidate_counts = _candidate_type_counts(nbo_candidates)
            primary_assignment = (_build_primary_assignment(
                molecule,
                nbo_candidates,
                compute_constraints,
            ) if compute_options.include_lewis_assignment else
                                  _assignment_payload(
                                      molecule,
                                      nbo_candidates,
                                      int(round(0.5 * molecule.number_of_electrons())),
                                      float(sum(candidate.get('occupation', 0.0)
                                                for candidate in nbo_candidates)),
                                      [
                                          'NBO candidate layer only: Lewis/NBO assignment was not requested.'
                                      ],
                                      compute_constraints,
                                  ))
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
            donor_acceptor_diagnostics = _build_donor_acceptor_diagnostics(
                nbo_candidates,
                primary_assignment,
                nao_data.density,
            )
            nra_results = (_build_nra_results(
                molecule,
                nao_data,
                spin_data,
                alternatives,
                compute_options.nra_subspace,
                compute_options.nra_fit_metric,
                prior_weights=compute_options.nra_prior_weights,
                prior_strength=compute_options.nra_prior_strength,
                prior_mode=compute_options.nra_prior_mode,
                spin_fit=compute_options.nra_spin_fit,
            ) if compute_options.include_nra else None)

            natural_charges = nuclear_charges.copy()
            for nao_index, atom in enumerate(nao_data.atom_map):
                natural_charges[atom] -= nao_data.populations[nao_index]

            diagnostics = (_build_foundation_diagnostics(
                molecule,
                overlap,
                density,
                nao_data,
                natural_charges,
                mo_analysis,
                candidate_counts,
                primary_counts,
                compute_constraints,
            ) if compute_options.include_diagnostics else {})
            if diagnostics:
                diagnostics['donor_acceptor_summary'] = {
                    'model': donor_acceptor_diagnostics['model'],
                    'donor_count': donor_acceptor_diagnostics['donor_count'],
                    'acceptor_count': donor_acceptor_diagnostics['acceptor_count'],
                    'interaction_count': donor_acceptor_diagnostics['interaction_count'],
                }

            results = {
                'nao_transform': nao_data.transform,
                'nao_density_matrix': nao_data.density,
                'nao_overlap_matrix': nao_data.overlap,
                'nao_populations': nao_data.populations,
                'nao_atom_map': nao_data.atom_map,
                'nao_l_map': nao_data.angular_momentum_map,
                'nao_alpha_density_matrix': spin_data.alpha_density,
                'nao_beta_density_matrix': spin_data.beta_density,
                'nao_spin_density_matrix': spin_data.spin_density,
                'nao_spin_populations': spin_data.spin_populations,
                'unpaired_electrons': spin_data.unpaired_electrons,
                'natural_charges': natural_charges,
                'mo_analysis': mo_analysis,
                'nbo_candidates': nbo_candidates,
                'nbo_list': primary_assignment['nbo_list'],
                'primary': primary_assignment,
                'alternatives': alternatives,
                'donor_acceptor_diagnostics': donor_acceptor_diagnostics,
                'diagnostics': diagnostics,
                'orbital_analysis': orbital_analysis,
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
                if 'nra' in results:
                    self.print_nra_report(molecule, results, level=self.nra_report_level)
                if hasattr(self.ostream, 'flush'):
                    self.ostream.flush()

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

    def set_nra_report_level(self, level):
        """Set NRA report level."""

        valid = {'none', 'summary', 'full'}
        assert_msg_critical(level in valid,
                            f"Invalid NRA report level '{level}'. "
                            f"Choose one of {sorted(valid)}")
        self.nra_report_level = level

    @staticmethod
    def _viewer_point(point):
        """Return a py3Dmol point dictionary."""

        return {
            'x': float(point[0]),
            'y': float(point[1]),
            'z': float(point[2]),
        }

    @staticmethod
    def _viewer_perpendicular(start, end, scale):
        """Return a stable perpendicular offset for drawing pi bonds."""

        axis = np.array(end, dtype=float) - np.array(start, dtype=float)
        norm = float(np.linalg.norm(axis))
        if norm < 1.0e-12:
            return np.array([0.0, 0.0, float(scale)])
        axis /= norm
        offset = np.cross(axis, np.array([0.0, 0.0, 1.0]))
        if float(np.linalg.norm(offset)) < 1.0e-12:
            offset = np.cross(axis, np.array([0.0, 1.0, 0.0]))
        offset_norm = float(np.linalg.norm(offset))
        if offset_norm < 1.0e-12:
            return np.array([0.0, 0.0, float(scale)])
        return float(scale) * offset / offset_norm

    @staticmethod
    def _viewer_marker_frame(coords, centroid, atom, plane_normal):
        """Return normal and tangent directions for atom-centered markers."""

        if plane_normal is not None:
            direction = np.array(plane_normal, dtype=float)
        else:
            direction = np.array(coords[atom], dtype=float) - np.array(centroid,
                                                                       dtype=float)
            direction_norm = float(np.linalg.norm(direction))
            if direction_norm < 1.0e-12:
                direction = np.array([0.0, 0.0, 1.0])
            else:
                direction /= direction_norm

        radial = np.array(coords[atom], dtype=float) - np.array(centroid,
                                                                dtype=float)
        tangent = radial - float(np.dot(radial, direction)) * direction
        tangent_norm = float(np.linalg.norm(tangent))
        if tangent_norm < 1.0e-12:
            tangent = np.cross(direction, np.array([1.0, 0.0, 0.0]))
            tangent_norm = float(np.linalg.norm(tangent))
        if tangent_norm < 1.0e-12:
            tangent = np.cross(direction, np.array([0.0, 1.0, 0.0]))
            tangent_norm = float(np.linalg.norm(tangent))
        if tangent_norm < 1.0e-12:
            tangent = np.array([1.0, 0.0, 0.0])
        else:
            tangent /= tangent_norm

        side = np.cross(direction, tangent)
        side_norm = float(np.linalg.norm(side))
        if side_norm < 1.0e-12:
            side = np.array([0.0, 1.0, 0.0])
        else:
            side /= side_norm

        return direction, tangent, side

    @staticmethod
    def _viewer_plane_normal(coords, labels):
        """Return a stable normal to the molecular plane."""

        def fit_normal(points):
            points = np.array(points, dtype=float)
            if len(points) < 3:
                return None
            centered = points - np.mean(points, axis=0)
            if float(np.linalg.norm(centered)) < 1.0e-12:
                return None
            try:
                _, singular_values, vh = np.linalg.svd(centered,
                                                       full_matrices=False)
            except np.linalg.LinAlgError:
                return None
            if len(singular_values) < 2 or singular_values[-2] < 1.0e-8:
                return None
            normal = np.array(vh[-1], dtype=float)
            normal_norm = float(np.linalg.norm(normal))
            if normal_norm < 1.0e-12:
                return None
            return normal / normal_norm

        heavy_points = [point for point, label in zip(coords, labels)
                        if str(label).strip().upper() != 'H']
        normal = fit_normal(heavy_points)
        if normal is None:
            normal = fit_normal(coords)
        if normal is None:
            return None
        if normal[2] < 0.0:
            normal *= -1.0
        return normal

    @staticmethod
    def _add_viewer_cylinder(viewer,
                             start,
                             end,
                             color,
                             radius,
                             opacity=1.0,
                             viewer_position=None):
        """Add a cylinder to a py3Dmol viewer or viewer-grid cell."""

        spec = {
            'start': NboDriver._viewer_point(start),
            'end': NboDriver._viewer_point(end),
            'radius': float(radius),
            'color': color,
            'opacity': float(opacity),
        }
        if viewer_position is None:
            viewer.addCylinder(spec)
        else:
            viewer.addCylinder(spec, viewer=viewer_position)

    @staticmethod
    def _add_viewer_sphere(viewer,
                           center,
                           color,
                           radius,
                           opacity=1.0,
                           viewer_position=None):
        """Add a sphere to a py3Dmol viewer or viewer-grid cell."""

        spec = {
            'center': NboDriver._viewer_point(center),
            'radius': float(radius),
            'color': color,
            'opacity': float(opacity),
        }
        if viewer_position is None:
            viewer.addSphere(spec)
        else:
            viewer.addSphere(spec, viewer=viewer_position)

    @staticmethod
    def _add_viewer_label(viewer,
                          text,
                          position,
                          color='black',
                          viewer_position=None):
        """Add a transparent label to a py3Dmol viewer or viewer-grid cell."""

        spec = {
            'position': NboDriver._viewer_point(position),
            'alignment': 'center',
            'fontColor': color,
            'backgroundColor': 0xffffff,
            'backgroundOpacity': 0.0,
        }
        if viewer_position is None:
            viewer.addLabel(str(text), spec)
        else:
            viewer.addLabel(str(text), spec, viewer=viewer_position)

    @staticmethod
    def _selected_structure_alternatives(results, ranks):
        """Return Lewis alternatives selected by one-based ranks."""

        alternatives = list(results.get('alternatives', ()))
        assert_msg_critical(bool(alternatives),
                            'show_structures: results contain no Lewis alternatives')
        if ranks is None:
            return alternatives
        if isinstance(ranks, int):
            ranks = [ranks]
        requested = {int(rank) for rank in ranks}
        selected = [alt for alt in alternatives
                    if int(alt.get('rank', 0)) in requested]
        assert_msg_critical(bool(selected),
                            f'show_structures: no alternatives found for ranks {sorted(requested)}')
        return selected

    def show_structures(self,
                        results=None,
                        molecule=None,
                        ranks=None,
                        width=420,
                        height=330,
                        max_columns=3,
                        show_sigma=True,
                        show_pi=True,
                        show_lone_pairs=True,
                        show_one_electron=True,
                        show_positive_centers=True,
                        show_atom_indices=True,
                        show_labels=True,
                        display=True):
        """Visualize NBO Lewis/resonance alternatives with py3Dmol.

        The method overlays the selected Lewis assignment on the molecular
        geometry: pi bonds and electron dots are darkcyan, lone-pair anion
        markers include two dots and a minus sign, and positive centers are
        dark-blue plus markers.  The input results are read only; no NBO
        analysis is recomputed.
        """

        if self.rank != mpi_master():
            return None

        try:
            import py3Dmol
        except ImportError:
            raise ImportError('show_structures requires py3Dmol')

        molecule = self._last_molecule if molecule is None else molecule
        results = self._last_results if results is None else results
        assert_msg_critical(molecule is not None,
                            'show_structures: no molecule provided and no prior compute() context available')
        assert_msg_critical(results is not None,
                            'show_structures: no results provided and no prior compute() context available')

        alternatives = self._selected_structure_alternatives(results, ranks)
        nviews = len(alternatives)
        ncols = min(max(1, int(max_columns)), nviews)
        nrows = int(np.ceil(float(nviews) / float(ncols)))
        viewer_grid = (nrows, ncols) if nviews > 1 else None
        if viewer_grid is None:
            view = py3Dmol.view(width=int(width), height=int(height))
        else:
            view = py3Dmol.view(width=int(width) * ncols,
                                height=int(height) * nrows,
                                viewergrid=viewer_grid)

        labels = molecule.get_labels()
        coords = np.array(molecule.get_coordinates_in_angstrom(), dtype=float)
        centroid = np.mean(coords, axis=0)
        span = float(np.max(np.linalg.norm(coords - centroid, axis=1)))
        label_position = centroid + np.array([0.0, 0.0, max(1.6, span + 0.8)])
        plane_normal = self._viewer_plane_normal(coords, labels)
        xyz = molecule.get_xyz_string()

        for index, alternative in enumerate(alternatives):
            viewer_position = None
            if viewer_grid is not None:
                viewer_position = (index // ncols, index % ncols)

            stick_style = {'radius': 0.13 if show_sigma else 0.08}
            if not show_sigma:
                stick_style['opacity'] = 0.25
            sphere_style = {'scale': 0.24 if show_sigma else 0.18}

            if viewer_position is None:
                view.addModel(xyz, 'xyz')
                view.setViewStyle({'style': 'outline', 'width': 0.05})
                view.setStyle({}, {
                    'stick': stick_style,
                    'sphere': sphere_style,
                })
            else:
                view.addModel(xyz, 'xyz', viewer=viewer_position)
                view.setViewStyle({'style': 'outline', 'width': 0.05},
                                  viewer=viewer_position)
                view.setStyle({}, {
                    'stick': stick_style,
                    'sphere': sphere_style,
                }, viewer=viewer_position)

            if show_pi:
                pi_pairs = [
                    tuple(int(atom) - 1 for atom in pair)
                    for pair in alternative.get('pi_bonds', [])
                    if len(pair) == 2
                ]
                for atoms in pi_pairs:
                    if min(atoms) < 0 or max(atoms) >= len(coords):
                        continue
                    axis = coords[atoms[1]] - coords[atoms[0]]
                    axis_norm = float(np.linalg.norm(axis))
                    if axis_norm < 1.0e-12:
                        continue
                    axis /= axis_norm
                    trim = min(0.28, 0.22 * axis_norm)
                    offset = None
                    if plane_normal is not None:
                        offset = 0.68 * plane_normal
                    if offset is None:
                        offset = self._viewer_perpendicular(coords[atoms[0]],
                                                            coords[atoms[1]],
                                                            0.68)
                    start = coords[atoms[0]] + trim * axis + offset
                    end = coords[atoms[1]] - trim * axis + offset
                    self._add_viewer_cylinder(view,
                                              start,
                                              end,
                                              '#005f5f',
                                              0.105,
                                              opacity=0.70,
                                              viewer_position=viewer_position)
                    for cap in (start, end):
                        self._add_viewer_sphere(view,
                                                cap,
                                                '#005f5f',
                                                0.105,
                                                opacity=0.70,
                                                viewer_position=viewer_position)
                    self._add_viewer_cylinder(view,
                                              start,
                                              end,
                                              '#008b8b',
                                              0.075,
                                              opacity=1.0,
                                              viewer_position=viewer_position)
                    for cap in (start, end):
                        self._add_viewer_sphere(view,
                                                cap,
                                                '#008b8b',
                                                0.075,
                                                opacity=1.0,
                                                viewer_position=viewer_position)

            if show_lone_pairs:
                lone_pair_atoms = list(alternative.get('active_lone_pair_atoms', []))
                if not lone_pair_atoms:
                    lone_pair_atoms = [
                        int(candidate.get('atoms', (0,))[0] + 1)
                        for candidate in alternative.get('active_lone_pair_nbo_list', [])
                    ]
                for atom_index in lone_pair_atoms:
                    atom = int(atom_index) - 1
                    if atom < 0 or atom >= len(coords):
                        continue
                    direction, tangent, _ = self._viewer_marker_frame(
                        coords, centroid, atom, plane_normal)
                    center = coords[atom] + 0.68 * direction
                    for shift in (-0.15, 0.15):
                        self._add_viewer_sphere(view,
                                                center + shift * tangent,
                                                '#008b8b',
                                                0.105,
                                                opacity=0.96,
                                                viewer_position=viewer_position)
                    minus_center = center + 0.28 * direction
                    minus_start = minus_center - 0.13 * tangent
                    minus_end = minus_center + 0.13 * tangent
                    self._add_viewer_cylinder(view,
                                              minus_start,
                                              minus_end,
                                              '#8b0000',
                                              0.030,
                                              opacity=1.0,
                                              viewer_position=viewer_position)
                    for cap in (minus_start, minus_end):
                        self._add_viewer_sphere(view,
                                                cap,
                                                '#8b0000',
                                                0.030,
                                                opacity=1.0,
                                                viewer_position=viewer_position)

            if show_one_electron:
                for atom_index in alternative.get('active_one_electron_atoms', []):
                    atom = int(atom_index) - 1
                    if atom < 0 or atom >= len(coords):
                        continue
                    direction, _, _ = self._viewer_marker_frame(
                        coords, centroid, atom, plane_normal)
                    center = coords[atom] + 0.68 * direction
                    self._add_viewer_sphere(view,
                                            center,
                                            '#008b8b',
                                            0.115,
                                            opacity=0.96,
                                            viewer_position=viewer_position)

            if show_positive_centers:
                for atom_index in alternative.get('active_positive_atoms', []):
                    atom = int(atom_index) - 1
                    if atom < 0 or atom >= len(coords):
                        continue
                    direction, tangent, _ = self._viewer_marker_frame(
                        coords, centroid, atom, plane_normal)
                    center = coords[atom] + 0.76 * direction
                    arm = 0.16
                    self._add_viewer_cylinder(view,
                                              center - arm * tangent,
                                              center + arm * tangent,
                                              '#003f8c',
                                              0.040,
                                              opacity=1.0,
                                              viewer_position=viewer_position)
                    self._add_viewer_cylinder(view,
                                              center - arm * direction,
                                              center + arm * direction,
                                              '#003f8c',
                                              0.040,
                                              opacity=1.0,
                                              viewer_position=viewer_position)

            if show_atom_indices:
                for atom, label in enumerate(labels):
                    if str(label).strip().upper() == 'H':
                        continue
                    atom_label = f'{label}{atom + 1}'
                    self._add_viewer_label(view,
                                           atom_label,
                                           coords[atom] + np.array([0.0, 0.0, 0.20]),
                                           viewer_position=viewer_position)

            if show_labels:
                rank = int(alternative.get('rank', index + 1))
                weight = float(alternative.get('weight', 0.0))
                title = f'Lewis {rank}: w={weight:.3f}'
                self._add_viewer_label(view,
                                       title,
                                       label_position,
                                       viewer_position=viewer_position)
            if viewer_position is None:
                view.zoomTo()
            else:
                view.zoomTo(viewer=viewer_position)

        if display:
            view.show()
            return None
        return view

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
        if hasattr(ostream, 'print_line'):
            ostream.print_blank()
            for line in text.splitlines():
                ostream.print_line(line)
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
        width = 104
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
        if hasattr(ostream, 'print_line'):
            ostream.print_blank()
            for line in text.splitlines():
                ostream.print_line(line)
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
        width = 104

        lines.append('Natural Bond Orbital (NBO) Primary Summary')
        lines.append('=' * width)
        lines.append('Lewis assignment from generated NBO candidates.')
        lines.append('Available layers: NAO/NPA, MO-in-NAO analysis, candidates, alternatives, diagnostics, and NRA/NRT.')
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
                    f"Primary alternative rank={primary.get('rank', '?')}  "
                    f"Score weight={primary.get('weight', 0.0):.5f}"
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
            show_one_electron = any(
                alternative.get('active_one_electron_atoms')
                for alternative in alternatives
            )
            show_lone_pair = any(
                alternative.get('active_lone_pair_atoms')
                for alternative in alternatives
            )
            show_positive = any(
                alternative.get('active_positive_atoms')
                for alternative in alternatives
            )
            extra_headers = []
            if show_one_electron:
                extra_headers.append(f"{'One-e':<12}")
            if show_lone_pair:
                extra_headers.append(f"{'LP':<12}")
            if show_positive:
                extra_headers.append(f"{'Pos':<12}")
            extra_text = ''.join(f'{header} ' for header in extra_headers)
            lines.append(
                f"{'Rank':>4} {'Score weight':>12} {'Score':>12} "
                f"{'Pairs':>7} {'Sigma/Pi':>9} {extra_text}{'Pi bonds'}"
            )
            lines.append('-' * width)
            for alternative in alternatives:
                pi_text = ', '.join(
                    f"{pair[0]}-{pair[1]}"
                    for pair in alternative.get('pi_bonds', [])
                ) or 'none'
                one_electron_text = ', '.join(
                    f"{labels[atom - 1]}{atom}"
                    for atom in alternative.get('active_one_electron_atoms', [])
                ) or 'none'
                lone_pair_text = ', '.join(
                    f"{labels[atom - 1]}{atom}"
                    for atom in alternative.get('active_lone_pair_atoms', [])
                ) or 'none'
                positive_text = ', '.join(
                    f"{labels[atom - 1]}{atom}"
                    for atom in alternative.get('active_positive_atoms', [])
                ) or 'none'
                prefix = (
                    f"{alternative.get('rank', 0):>4} "
                    f"{alternative.get('weight', 0.0):>12.5f} "
                    f"{alternative.get('score', 0.0):>12.5f} "
                    f"{alternative.get('electron_pairs', 0):>7}/"
                    f"{alternative.get('target_electron_pairs', '?'):<3} "
                    f"{alternative.get('sigma_electron_pairs', 0):>4}/"
                    f"{alternative.get('pi_electron_pairs', 0):<4} "
                )
                if show_one_electron:
                    prefix += f"{one_electron_text:<12} "
                if show_lone_pair:
                    prefix += f"{lone_pair_text:<12} "
                if show_positive:
                    prefix += f"{positive_text:<12} "
                lines.append(
                    prefix + pi_text
                )
            alternative_warnings = []
            for alternative in alternatives:
                alternative_warnings.extend(alternative.get('warnings', []))
            if alternative_warnings:
                lines.append('')
                lines.append('Alternative notes')
                for warning in dict.fromkeys(alternative_warnings):
                    lines.append(f'- {warning}')

        lines.append('=' * width)
        return '\n'.join(lines)

    def print_nbo_report(self, molecule, results, level=None, ostream=None):
        """Print an NBO report."""

        text = self.format_nbo_report(molecule, results, level=level)
        if not text:
            return
        if ostream is None:
            ostream = self.ostream
        if hasattr(ostream, 'print_line'):
            for line in text.splitlines():
                ostream.print_line(line)
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

    def get_nra_data(self, molecule, results):
        """Return structured NRA/NRT reporting data."""

        assert_msg_critical('nra' in results,
                            'get_nra_data: results do not contain NRA data; '
                            'run compute() with include_nra=True')
        nra = results['nra']
        assert_msg_critical(bool(nra),
                            'get_nra_data: NRA data are empty')
        return nra

    def format_nra_report(self, molecule, results, level=None):
        """Return an NRA/NRT density-fit report string."""

        if level is None:
            level = self.nra_report_level
        valid = {'none', 'summary', 'full'}
        assert_msg_critical(level in valid,
                            f"Invalid NRA report level '{level}'. "
                            f"Choose one of {sorted(valid)}")
        if level == 'none':
            return ''

        data = self.get_nra_data(molecule, results)
        structures = data.get('structures', [])
        lines = []
        width = 104

        lines.append('Natural Resonance Analysis (NRA/NRT) Density-Fit Summary')
        lines.append('=' * width)
        lines.append(
            'NRA weights are density-fit weights; Score weights are Lewis ranking weights.'
        )
        lines.append('')
        lines.append(
            f"Subspace: {data.get('subspace', 'n/a')}  "
            f"Metric: {data.get('fit_metric', 'n/a')}  "
            f"Residual: {data.get('residual_norm', 0.0):.6e}  "
            f"Relative residual: {data.get('relative_residual', 0.0):.6e}"
        )
        lines.append(
            f"Spin fit: {data.get('spin_fit', 'none')}  "
            f"Total residual: {data.get('total_residual_norm', 0.0):.6e}  "
            f"Spin residual: {data.get('spin_residual_norm', 0.0) if data.get('spin_residual_norm') is not None else 'n/a'}"
        )
        prior = data.get('prior', {})
        lines.append(
            f"Prior mode: {prior.get('mode', 'regularized')}  "
            f"Prior strength: {prior.get('strength', 0.0):.6g}  "
            f"Prior active: {bool(prior.get('active', False))}"
        )

        for warning in data.get('warnings', []):
            lines.append(f'Warning: {warning}')

        lines.append('')
        show_prior = any(structure.get('prior_weight') is not None
                         for structure in structures)
        if show_prior:
            lines.append(
                f"{'Rank':>4} {'Signature':<18} {'NRA weight':>12} "
                f"{'Prior weight':>12} {'Score weight':>13} {'Score':>12} "
                f"{'Residual':>12}  {'Sigma/Pi':>9} {'Pi bonds'}"
            )
        else:
            lines.append(
                f"{'Rank':>4} {'Signature':<18} {'NRA weight':>12} "
                f"{'Score weight':>13} {'Score':>12} {'Residual':>12}  "
                f"{'Sigma/Pi':>9} {'Pi bonds'}"
            )
        lines.append('-' * width)
        if structures:
            for structure in structures:
                pi_text = ', '.join(
                    f"{pair[0]}-{pair[1]}"
                    for pair in structure.get('pi_bonds', [])
                ) or 'none'
                common = (
                    f"{structure.get('rank', 0):>4} "
                    f"{structure.get('label', 'n/a'):<18.18} "
                    f"{structure.get('nra_weight', 0.0):>12.6f} "
                )
                if show_prior:
                    prior_weight = structure.get('prior_weight')
                    prior_text = ('n/a' if prior_weight is None else
                                  f'{prior_weight:12.6f}')
                    common += f"{prior_text:>12} "
                lines.append(
                    common +
                    f"{structure.get('score_weight', 0.0):>13.6f} "
                    f"{structure.get('score', 0.0):>12.5f} "
                    f"{structure.get('residual_norm', 0.0):>12.5e}  "
                    f"{structure.get('sigma_electron_pairs', 0):>4}/"
                    f"{structure.get('pi_electron_pairs', 0):<4} "
                    f"{pi_text}"
                )
        else:
            lines.append(f"{'':>4} {'No NRA structures available.':<60}")

        if level == 'full' and structures:
            lines.append('')
            lines.append('Structure details')
            lines.append('-' * width)
            alternatives = results.get('alternatives', [])
            alternatives_by_rank = {
                int(alternative.get('rank', 0)): alternative
                for alternative in alternatives
            }
            for structure in structures:
                rank = int(structure.get('rank', 0))
                alternative = alternatives_by_rank.get(rank, {})
                counts = _candidate_type_counts(
                    alternative.get('nbo_list', [])
                )
                if counts:
                    count_text = ', '.join(
                        f'{key}={value}' for key, value in sorted(counts.items())
                    )
                else:
                    count_text = 'n/a'
                electron_pairs = alternative.get('electron_pairs', 'n/a')
                target_pairs = alternative.get('target_electron_pairs', 'n/a')
                lines.append(
                    f"Rank {rank}: selected NBO counts: {count_text}; "
                    f"pairs={electron_pairs}/{target_pairs}"
                )

        lines.append('=' * width)
        return '\n'.join(lines)

    def print_nra_report(self, molecule, results, level=None, ostream=None):
        """Print an NRA/NRT density-fit report."""

        text = self.format_nra_report(molecule, results, level=level)
        if not text:
            return
        if ostream is None:
            ostream = self.ostream
        if hasattr(ostream, 'print_line'):
            for line in text.splitlines():
                ostream.print_line(line)
            ostream.print_blank()
        else:
            print(text, file=ostream if ostream is not None else sys.stdout)

    def nra_report(self,
                   level=None,
                   molecule=None,
                   results=None,
                   ostream=None,
                   return_text=False):
        """Convenience NRA/NRT reporting API."""

        molecule = self._last_molecule if molecule is None else molecule
        results = self._last_results if results is None else results
        assert_msg_critical(molecule is not None,
                            'nra_report: no molecule provided and no prior compute() context available')
        assert_msg_critical(results is not None,
                            'nra_report: no results provided and no prior compute() context available')
        text = self.format_nra_report(molecule, results, level=level)
        if return_text:
            return text
        if text:
            print(text, file=ostream if ostream is not None else sys.stdout)
        return text