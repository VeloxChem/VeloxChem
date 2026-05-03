"""Shared NBO-like orbital candidate classification helpers.

This module is an implementation detail used by :mod:`orbitalanalyzerdriver`.
Public NBO and VB workflows should access these records through
``OrbitalAnalyzer`` so both drivers share the same analysis payload.
"""

import importlib

import numpy as np

from .veloxchemlib import chemical_element_identifier


def _nbo_helpers():
    return importlib.import_module(f'{__package__}.nbodriver')


def _build_nbo_candidates(molecule,
                          nao_data,
                          spin_data=None,
                          lone_pair_min_occupation=1.50,
                          rydberg_max_occupation=0.50,
                          bond_min_occupation=1.20,
                          bond_min_atom_weight=0.10,
                          pi_min_occupation=0.20,
                          conjugated_pi_max_path=2,
                          constraints=None):
    """Build NBO candidates from the orthonormal NAO density.

    This is a candidate layer, not a Lewis assignment layer.  Candidates are
    deterministic one-center core/lone-pair orbitals and two-center natural
    orbitals obtained by diagonalizing connected atom-pair density blocks after
    removing core and obvious lone-pair NAOs from the bond search space.
    """

    nbo = _nbo_helpers()

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

    return candidates
