#
#                                   VELOXCHEM
#              ----------------------------------------------------
#                          An Electronic Structure Code
#
#  SPDX-License-Identifier: BSD-3-Clause
#
#  Copyright 2018-2026 VeloxChem developers
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
"""TREXIO readers and writers for VeloxChem data."""

from pathlib import Path
import importlib
import shutil

import numpy as np

from .errorhandler import assert_msg_critical
from .molecule import Molecule
from .molecularbasis import MolecularBasis
from .molecularorbitals import MolecularOrbitals, molorb
from .veloxchemlib import AtomBasis, AtomCorePotential, BaseCorePotential
from .veloxchemlib import BasisFunction

_ANGULAR_LABELS = 'SPDFGHIKLMNOQRTUVWXYZ'


def _trexio_module():
    """Imports the external TREXIO Python module lazily."""

    try:
        return importlib.import_module('trexio')
    except ImportError as exc:
        raise ImportError(
            'TREXIO support requires the external Python package "trexio". '
            'Install it with conda-forge::trexio and, if needed, '
            'pip install trexio.') from exc


def _get_backend(backend):
    """Converts backend label to TREXIO backend constant."""

    trexio = _trexio_module()

    if isinstance(backend, str):
        backend_label = backend.lower()
        if backend_label in ('hdf5', 'h5'):
            return trexio.TREXIO_HDF5
        if backend_label in ('text', 'txt'):
            return trexio.TREXIO_TEXT
        assert_msg_critical(False, f'Unsupported TREXIO backend {backend!r}.')

    return backend


def _remove_existing_trexio_path(fname):
    """Removes an existing TREXIO file or text-backend directory."""

    path = Path(fname)
    if path.is_dir():
        shutil.rmtree(path)
    elif path.exists():
        path.unlink()


def _open_trexio(fname, mode='r', backend='hdf5'):
    """Opens a TREXIO file."""

    trexio = _trexio_module()
    return trexio.File(str(fname), mode=mode, back_end=_get_backend(backend))


def _safe_close(trexio_file):
    """Closes a TREXIO file handle."""

    if trexio_file is not None:
        trexio_file.close()


def _has(trexio_file, name):
    """Calls a TREXIO has_* function if it exists."""

    trexio = _trexio_module()
    has_func = getattr(trexio, f'has_{name}', None)
    return has_func is not None and has_func(trexio_file)


def _read_if_present(trexio_file, name, default=None):
    """Reads a TREXIO dataset if present."""

    if not _has(trexio_file, name):
        return default

    trexio = _trexio_module()
    return getattr(trexio, f'read_{name}')(trexio_file)


def _angular_momentum_from_label(label):
    """Returns angular momentum for a VeloxChem AO label token."""

    for char in str(label).upper():
        if char in _ANGULAR_LABELS:
            return _ANGULAR_LABELS.index(char)
    assert_msg_critical(False,
                        f'Could not parse AO angular momentum {label!r}.')


def _shell_ordinal_from_label(label):
    """Returns one-based shell ordinal from a VeloxChem AO label token."""

    digits = ''
    for char in str(label):
        if char.isdigit():
            digits += char
        elif digits:
            break
    return int(digits) if digits else 1


def _shell_degeneracy(angular_momentum, cartesian=False):
    """Returns number of AOs in a shell."""

    if cartesian:
        return (angular_momentum + 1) * (angular_momentum + 2) // 2
    return 2 * angular_momentum + 1


def _spherical_m_sequence(angular_momentum):
    """Returns the TREXIO canonical real-spherical m ordering."""

    sequence = [0]
    for m_value in range(1, angular_momentum + 1):
        sequence.extend([m_value, -m_value])

    return sequence


def _spherical_m_rank(angular_momentum, m_value):
    """Returns the TREXIO canonical rank of a spherical AO m value."""

    sequence = _spherical_m_sequence(angular_momentum)
    assert_msg_critical(
        m_value in sequence,
        f'Invalid m={m_value} for angular momentum {angular_momentum}.')
    return sequence.index(m_value)


def _m_value_from_label(label):
    """Returns the m value encoded in a VeloxChem AO label token."""

    token = str(label).strip()
    if '+' in token:
        return int(token.split('+')[-1])
    if '-' in token:
        return -int(token.split('-')[-1])
    return 0


def _primitive_normalization_factors(exponents, angular_momentum):
    """Returns primitive spherical-Gaussian normalization factors."""

    exponents = np.array(exponents, dtype=float)
    factors = np.power(2.0 * exponents / np.pi, 0.75)

    if angular_momentum == 0:
        return factors
    if angular_momentum == 1:
        return factors * 2.0 * np.sqrt(exponents)
    if angular_momentum == 2:
        return factors * (4.0 / np.sqrt(3.0)) * exponents
    if angular_momentum == 3:
        return factors * (8.0 / np.sqrt(15.0)) * exponents * np.sqrt(
            exponents)
    if angular_momentum == 4:
        return factors * (16.0 / np.sqrt(105.0)) * exponents * exponents
    if angular_momentum == 5:
        return factors * (32.0 / np.sqrt(945.0)) * exponents * exponents * np.sqrt(
            exponents)
    if angular_momentum == 6:
        return factors * (64.0 / np.sqrt(10395.0)) * exponents * exponents * exponents

    return np.ones_like(exponents)


def _basis_shells(molecule, basis):
    """Returns basis shell records in TREXIO shell order."""

    basis_sets = basis.basis_sets()
    basis_indices = basis.basis_sets_indices()
    shells = []

    for atom_index in range(molecule.number_of_atoms()):
        atom_basis = basis_sets[basis_indices[atom_index]]
        for basis_function in atom_basis.get_basis_functions():
            shells.append({
                'atom_index': atom_index,
                'angular_momentum': basis_function.get_angular_momentum(),
                'exponents': np.array(basis_function.get_exponents(),
                                      dtype=float),
                'coefficients': np.array(
                    basis_function.get_normalization_factors(), dtype=float),
            })

    return shells


def _shell_lookup_by_atom_angular(shells):
    """Builds a lookup for shell index from atom, angular momentum, ordinal."""

    counters = {}
    lookup = {}

    for shell_index, shell in enumerate(shells):
        key = (shell['atom_index'], shell['angular_momentum'])
        counters[key] = counters.get(key, 0) + 1
        lookup[(shell['atom_index'], shell['angular_momentum'],
                counters[key])] = shell_index

    return lookup


def _ao_shell_indices(molecule, basis, shells):
    """Gets TREXIO AO shell indices in VeloxChem AO order."""

    shell_lookup = _shell_lookup_by_atom_angular(shells)
    ao_shell = []

    for ao_label in basis.get_ao_basis_map(molecule):
        fields = ao_label.split()
        assert_msg_critical(
            len(fields) >= 3,
            f'Could not parse VeloxChem AO basis label {ao_label!r}.')
        atom_index = int(fields[0]) - 1
        angular_momentum = _angular_momentum_from_label(fields[2])
        shell_ordinal = _shell_ordinal_from_label(fields[2])
        key = (atom_index, angular_momentum, shell_ordinal)
        assert_msg_critical(
            key in shell_lookup,
            f'Could not locate basis shell for AO {ao_label!r}.')
        ao_shell.append(shell_lookup[key])

    return np.array(ao_shell, dtype=np.int64)


def _ao_permutation_to_trexio_order(molecule, basis, shells):
    """Returns TREXIO-to-VeloxChem AO permutation and AO shell indices."""

    shell_lookup = _shell_lookup_by_atom_angular(shells)
    ao_records = []

    for vlx_index, ao_label in enumerate(basis.get_ao_basis_map(molecule)):
        fields = ao_label.split()
        assert_msg_critical(
            len(fields) >= 3,
            f'Could not parse VeloxChem AO basis label {ao_label!r}.')
        atom_index = int(fields[0]) - 1
        angular_momentum = _angular_momentum_from_label(fields[2])
        shell_ordinal = _shell_ordinal_from_label(fields[2])
        shell_index = shell_lookup[(atom_index, angular_momentum,
                                    shell_ordinal)]
        m_rank = _spherical_m_rank(angular_momentum,
                                   _m_value_from_label(fields[2]))
        ao_records.append((shell_index, m_rank, vlx_index))

    ao_records.sort(key=lambda item: (item[0], item[1]))
    trexio_to_vlx = np.array([item[2] for item in ao_records], dtype=np.int64)
    ao_shell = np.array([item[0] for item in ao_records], dtype=np.int64)

    return trexio_to_vlx, ao_shell


def _inverse_permutation(permutation):
    """Returns inverse permutation."""

    inverse = np.empty_like(permutation)
    inverse[permutation] = np.arange(permutation.size, dtype=permutation.dtype)
    return inverse


def _write_metadata(trexio_file, description=None):
    """Writes TREXIO metadata."""

    trexio = _trexio_module()

    trexio.write_metadata_code_num(trexio_file, 1)
    trexio.write_metadata_code(trexio_file, ['VeloxChem'])

    if description is not None:
        trexio.write_metadata_description(trexio_file, str(description))


def _write_molecule(trexio_file, molecule, basis=None):
    """Writes molecule and electron information."""

    trexio = _trexio_module()

    labels = list(molecule.get_labels())
    charges = np.array(molecule.get_element_ids(), dtype=float)
    coordinates = np.array(molecule.get_coordinates_in_bohr(), dtype=float)

    trexio.write_nucleus_num(trexio_file, molecule.number_of_atoms())
    trexio.write_nucleus_charge(trexio_file, charges)
    trexio.write_nucleus_coord(trexio_file, coordinates)
    trexio.write_nucleus_label(trexio_file, labels)

    if basis is not None:
        trexio.write_nucleus_repulsion(
            trexio_file,
            float(molecule.effective_nuclear_repulsion_energy(basis)))

    nalpha = int(molecule.number_of_alpha_electrons())
    nbeta = int(molecule.number_of_beta_electrons())
    trexio.write_electron_num(trexio_file, nalpha + nbeta)
    trexio.write_electron_up_num(trexio_file, nalpha)
    trexio.write_electron_dn_num(trexio_file, nbeta)


def _write_basis(trexio_file, molecule, basis):
    """Writes Gaussian basis information."""

    trexio = _trexio_module()
    shells = _basis_shells(molecule, basis)
    prim_num = sum(shell['exponents'].size for shell in shells)

    shell_index = []
    exponent = []
    coefficient = []
    prim_factor = []

    for shell_idx, shell in enumerate(shells):
        primitive_factors = _primitive_normalization_factors(
            shell['exponents'], shell['angular_momentum'])
        shell_index.extend([shell_idx] * shell['exponents'].size)
        exponent.extend(shell['exponents'].tolist())
        coefficient.extend((shell['coefficients'] / primitive_factors).tolist())
        prim_factor.extend(primitive_factors.tolist())

    trexio.write_basis_type(trexio_file, 'Gaussian')
    trexio.write_basis_shell_num(trexio_file, len(shells))
    trexio.write_basis_prim_num(trexio_file, prim_num)
    trexio.write_basis_nucleus_index(
        trexio_file,
        np.array([shell['atom_index'] for shell in shells], dtype=np.int64))
    trexio.write_basis_shell_ang_mom(
        trexio_file,
        np.array([shell['angular_momentum'] for shell in shells],
                 dtype=np.int64))
    trexio.write_basis_shell_factor(trexio_file, np.ones(len(shells)))
    trexio.write_basis_r_power(trexio_file, np.zeros(len(shells),
                                                     dtype=np.int64))
    trexio.write_basis_shell_index(trexio_file,
                                   np.array(shell_index, dtype=np.int64))
    trexio.write_basis_exponent(trexio_file, np.array(exponent, dtype=float))
    trexio.write_basis_coefficient(trexio_file,
                                   np.array(coefficient, dtype=float))
    trexio.write_basis_prim_factor(trexio_file,
                                   np.array(prim_factor, dtype=float))

    return shells


def _write_ao(trexio_file, molecule, basis, shells):
    """Writes AO metadata."""

    trexio = _trexio_module()

    nao = int(basis.get_dimensions_of_basis())
    _, ao_shell = _ao_permutation_to_trexio_order(molecule, basis, shells)

    trexio.write_ao_num(trexio_file, nao)
    trexio.write_ao_cartesian(trexio_file, 0)
    trexio.write_ao_shell(trexio_file, ao_shell)
    trexio.write_ao_normalization(trexio_file, np.ones(nao))


def _write_ecp(trexio_file, molecule, basis):
    """Writes ECP information when present."""

    if not basis.has_ecp():
        return

    trexio = _trexio_module()
    basis_sets = basis.basis_sets()
    basis_indices = basis.basis_sets_indices()

    max_ang_mom_plus_1 = []
    z_core = []
    nucleus_index = []
    ang_mom = []
    exponent = []
    coefficient = []
    power = []

    for atom_index in range(molecule.number_of_atoms()):
        atom_basis = basis_sets[basis_indices[atom_index]]

        if not atom_basis.has_ecp():
            max_ang_mom_plus_1.append(0)
            z_core.append(0)
            continue

        atom_ecp = atom_basis.get_ecp_potential()
        projected_angular = list(atom_ecp.get_angular_momentums())
        local_angular = (max(projected_angular) + 1 if projected_angular else 0)
        max_ang_mom_plus_1.append(local_angular)
        z_core.append(int(atom_ecp.number_of_core_electrons()))

        local_potential = atom_ecp.get_local_potential()
        for expn, coef, radial_order in zip(
                local_potential.get_exponents(), local_potential.get_factors(),
                local_potential.get_radial_orders()):
            nucleus_index.append(atom_index)
            ang_mom.append(local_angular)
            exponent.append(expn)
            coefficient.append(coef)
            power.append(radial_order)

        for proj_ang, projected_potential in zip(
                projected_angular, atom_ecp.get_projected_potentials()):
            for expn, coef, radial_order in zip(
                    projected_potential.get_exponents(),
                    projected_potential.get_factors(),
                    projected_potential.get_radial_orders()):
                nucleus_index.append(atom_index)
                ang_mom.append(proj_ang)
                exponent.append(expn)
                coefficient.append(coef)
                power.append(radial_order)

    if not nucleus_index:
        return

    trexio.write_ecp_max_ang_mom_plus_1(
        trexio_file, np.array(max_ang_mom_plus_1, dtype=np.int64))
    trexio.write_ecp_z_core(trexio_file, np.array(z_core, dtype=np.int64))
    trexio.write_ecp_num(trexio_file, len(nucleus_index))
    trexio.write_ecp_nucleus_index(trexio_file,
                                   np.array(nucleus_index, dtype=np.int64))
    trexio.write_ecp_ang_mom(trexio_file, np.array(ang_mom, dtype=np.int64))
    trexio.write_ecp_exponent(trexio_file, np.array(exponent, dtype=float))
    trexio.write_ecp_coefficient(trexio_file, np.array(coefficient,
                                                       dtype=float))
    trexio.write_ecp_power(trexio_file, np.array(power, dtype=np.int64))


def _molecular_orbitals_from_scf_results(scf_results):
    """Creates MolecularOrbitals from an SCF results dictionary."""

    scf_type = scf_results.get('scf_type', 'restricted')

    if scf_type == 'unrestricted':
        return MolecularOrbitals(
            [scf_results['C_alpha'], scf_results['C_beta']],
            [scf_results['E_alpha'], scf_results['E_beta']],
            [scf_results['occ_alpha'], scf_results['occ_beta']], molorb.unrest)

    if scf_type == 'restricted_openshell':
        return MolecularOrbitals(
            [scf_results['C_alpha']], [scf_results['E_alpha']],
            [scf_results['occ_alpha'], scf_results['occ_beta']],
            molorb.restopen)

    return MolecularOrbitals([scf_results['C_alpha']], [scf_results['E_alpha']],
                             [scf_results['occ_alpha']], molorb.rest)


def _write_molecular_orbitals(trexio_file,
                              mol_orbs,
                              mo_type='VeloxChem',
                              ao_permutation=None):
    """Writes molecular orbitals."""

    trexio = _trexio_module()

    orbitals_type = mol_orbs.get_orbitals_type()
    if ao_permutation is None:
        ao_permutation = np.arange(mol_orbs.alpha_to_numpy().shape[0])

    if orbitals_type == molorb.unrest:
        coef = np.vstack(
            (mol_orbs.alpha_to_numpy()[ao_permutation, :].T,
             mol_orbs.beta_to_numpy()[ao_permutation, :].T))
        energy = np.hstack((mol_orbs.ea_to_numpy(), mol_orbs.eb_to_numpy()))
        occupation = np.hstack(
            (mol_orbs.occa_to_numpy(), mol_orbs.occb_to_numpy()))
        spin = np.hstack((np.zeros(mol_orbs.number_of_mos(), dtype=np.int64),
                          np.ones(mol_orbs.number_of_mos(), dtype=np.int64)))
    else:
        coef = mol_orbs.alpha_to_numpy()[ao_permutation, :].T
        energy = mol_orbs.ea_to_numpy()
        occupation = mol_orbs.occa_to_numpy() + mol_orbs.occb_to_numpy()
        spin = np.zeros(mol_orbs.number_of_mos(), dtype=np.int64)

    trexio.write_mo_type(trexio_file, str(mo_type))
    trexio.write_mo_num(trexio_file, coef.shape[0])
    trexio.write_mo_coefficient(trexio_file, np.array(coef, dtype=float))
    trexio.write_mo_energy(trexio_file, np.array(energy, dtype=float))
    trexio.write_mo_occupation(trexio_file, np.array(occupation, dtype=float))
    trexio.write_mo_spin(trexio_file, spin)


def _write_ao_integrals(trexio_file, scf_results, ao_permutation=None):
    """Writes available AO integral-like matrices from SCF results."""

    trexio = _trexio_module()

    if 'S' in scf_results:
        overlap = np.array(scf_results['S'], dtype=float)
        if ao_permutation is not None:
            overlap = overlap[np.ix_(ao_permutation, ao_permutation)]
        trexio.write_ao_1e_int_overlap(trexio_file, overlap)


def _write_state(trexio_file, scf_results):
    """Writes ground-state energy if available."""

    if 'scf_energy' not in scf_results:
        return

    trexio = _trexio_module()
    trexio.write_state_num(trexio_file, 1)
    trexio.write_state_id(trexio_file, 0)
    trexio.write_state_current_label(trexio_file, 'ground state')
    trexio.write_state_energy(trexio_file, float(scf_results['scf_energy']))


def write_trexio(fname,
                 molecule,
                 basis=None,
                 scf_results=None,
                 mol_orbs=None,
                 backend='hdf5',
                 description=None):
    """Writes VeloxChem data to a TREXIO file.

    :param fname:
        The TREXIO file name.
    :param molecule:
        The molecule.
    :param basis:
        The AO basis set.
    :param scf_results:
        Optional SCF results dictionary.
    :param mol_orbs:
        Optional molecular orbitals. If omitted and ``scf_results`` contains
        orbital arrays, orbitals are reconstructed from the SCF results.
    :param backend:
        ``'hdf5'`` or ``'text'``.
    :param description:
        Optional metadata description.
    """

    _remove_existing_trexio_path(fname)

    trexio_file = None
    try:
        trexio_file = _open_trexio(fname, mode='w', backend=backend)
        _write_metadata(trexio_file, description)
        _write_molecule(trexio_file, molecule, basis)

        ao_permutation = None
        if basis is not None:
            shells = _write_basis(trexio_file, molecule, basis)
            ao_permutation, _ = _ao_permutation_to_trexio_order(
                molecule, basis, shells)
            _write_ao(trexio_file, molecule, basis, shells)
            _write_ecp(trexio_file, molecule, basis)

        if mol_orbs is None and scf_results is not None and 'C_alpha' in scf_results:
            mol_orbs = _molecular_orbitals_from_scf_results(scf_results)

        if mol_orbs is not None:
            mo_type = 'VeloxChem'
            if scf_results is not None:
                mo_type = scf_results.get('scf_type', mo_type)
            _write_molecular_orbitals(trexio_file, mol_orbs, mo_type,
                                      ao_permutation)

        if scf_results is not None:
            _write_ao_integrals(trexio_file, scf_results, ao_permutation)
            _write_state(trexio_file, scf_results)

        trexio_file.flush()
    finally:
        _safe_close(trexio_file)


def read_molecule(fname, backend='hdf5'):
    """Reads a molecule from a TREXIO file."""

    trexio_file = None
    try:
        trexio_file = _open_trexio(fname, mode='r', backend=backend)

        charges = np.array(_read_if_present(trexio_file, 'nucleus_charge'),
                           dtype=float)
        coordinates = np.array(_read_if_present(trexio_file, 'nucleus_coord'),
                               dtype=float)
        labels = _read_if_present(trexio_file, 'nucleus_label', None)

        assert_msg_critical(charges.size > 0,
                            'TREXIO file does not contain nucleus charges.')
        assert_msg_critical(
            coordinates.size > 0,
            'TREXIO file does not contain nucleus coordinates.')

        identifiers = np.rint(charges).astype(int)
        if labels is None:
            atom_basis_labels = None
        else:
            atom_basis_labels = [('', str(label)) for label in labels]

        molecule = Molecule(identifiers, coordinates, 'au', atom_basis_labels)

        electron_num = _read_if_present(trexio_file, 'electron_num', None)
        if electron_num is not None:
            molecule.set_charge(int(round(np.sum(identifiers) - electron_num)))

        up_num = _read_if_present(trexio_file, 'electron_up_num', None)
        dn_num = _read_if_present(trexio_file, 'electron_dn_num', None)
        if up_num is not None and dn_num is not None:
            molecule.set_multiplicity(abs(int(up_num) - int(dn_num)) + 1)

        return molecule
    finally:
        _safe_close(trexio_file)


def _read_basis_from_open_trexio(trexio_file, molecule):
    """Reads a MolecularBasis from an open TREXIO file."""

    basis_type = _read_if_present(trexio_file, 'basis_type')
    assert_msg_critical(
        str(basis_type).lower() == 'gaussian',
        f'Unsupported TREXIO basis type {basis_type!r}.')

    shell_num = int(_read_if_present(trexio_file, 'basis_shell_num'))
    shell_ang_mom = np.array(_read_if_present(trexio_file,
                                              'basis_shell_ang_mom'),
                             dtype=int)
    nucleus_index = np.array(_read_if_present(trexio_file,
                                              'basis_nucleus_index'),
                             dtype=int)
    shell_index = np.array(_read_if_present(trexio_file, 'basis_shell_index'),
                           dtype=int)
    exponent = np.array(_read_if_present(trexio_file, 'basis_exponent'),
                        dtype=float)
    coefficient = np.array(_read_if_present(trexio_file, 'basis_coefficient'),
                           dtype=float)
    prim_factor = np.array(_read_if_present(trexio_file, 'basis_prim_factor',
                                            np.ones_like(coefficient)),
                           dtype=float)
    shell_factor = np.array(_read_if_present(trexio_file, 'basis_shell_factor',
                                             np.ones(shell_num)),
                            dtype=float)

    basis = MolecularBasis()

    for atom_index in range(molecule.number_of_atoms()):
        atom_basis = AtomBasis()
        atom_basis.set_identifier(int(molecule.get_identifiers()[atom_index]))
        atom_basis.set_name('TREXIO')

        atom_shell_indices = [
            idx for idx in range(shell_num) if nucleus_index[idx] == atom_index
        ]
        for atom_shell_index in atom_shell_indices:
            primitive_indices = np.where(shell_index == atom_shell_index)[0]
            factors = (coefficient[primitive_indices] *
                       prim_factor[primitive_indices] *
                       shell_factor[atom_shell_index])
            basis_function = BasisFunction(exponent[primitive_indices].tolist(),
                                           factors.tolist(),
                                           int(shell_ang_mom[atom_shell_index]))
            atom_basis.add(basis_function)

        basis.add(atom_basis)

    _read_ecp_into_basis(trexio_file, basis)

    return basis


def _read_ecp_into_basis(trexio_file, basis):
    """Reads ECP data into an existing molecular basis when present."""

    if not _has(trexio_file, 'ecp_num'):
        return

    max_ang_mom_plus_1 = np.array(_read_if_present(trexio_file,
                                                   'ecp_max_ang_mom_plus_1'),
                                  dtype=int)
    z_core = np.array(_read_if_present(trexio_file, 'ecp_z_core'), dtype=int)
    nucleus_index = np.array(_read_if_present(trexio_file, 'ecp_nucleus_index'),
                             dtype=int)
    ang_mom = np.array(_read_if_present(trexio_file, 'ecp_ang_mom'), dtype=int)
    exponent = np.array(_read_if_present(trexio_file, 'ecp_exponent'),
                        dtype=float)
    coefficient = np.array(_read_if_present(trexio_file, 'ecp_coefficient'),
                           dtype=float)
    power = np.array(_read_if_present(trexio_file, 'ecp_power'), dtype=int)

    for atom_index, atom_basis in enumerate(basis.basis_sets()):
        if atom_index >= max_ang_mom_plus_1.size or z_core[atom_index] == 0:
            continue

        local_ang = int(max_ang_mom_plus_1[atom_index])
        atom_items = np.where(nucleus_index == atom_index)[0]
        local_items = atom_items[ang_mom[atom_items] == local_ang]

        local_potential = BaseCorePotential(exponent[local_items].tolist(),
                                            coefficient[local_items].tolist(),
                                            power[local_items].tolist())

        projected_potentials = []
        projected_angular = []
        for proj_ang in sorted(set(ang_mom[atom_items].tolist())):
            if proj_ang == local_ang:
                continue
            proj_items = atom_items[ang_mom[atom_items] == proj_ang]
            projected_angular.append(int(proj_ang))
            projected_potentials.append(
                BaseCorePotential(exponent[proj_items].tolist(),
                                  coefficient[proj_items].tolist(),
                                  power[proj_items].tolist()))

        if projected_potentials:
            atom_ecp = AtomCorePotential(local_potential, projected_potentials,
                                         projected_angular,
                                         int(z_core[atom_index]))
            atom_basis.set_ecp_potential(atom_ecp)


def read_molecule_and_basis(fname, backend='hdf5'):
    """Reads molecule and basis from a TREXIO file."""

    molecule = read_molecule(fname, backend)
    trexio_file = None
    try:
        trexio_file = _open_trexio(fname, mode='r', backend=backend)
        basis = _read_basis_from_open_trexio(trexio_file, molecule)
        return molecule, basis
    finally:
        _safe_close(trexio_file)


def read_molecular_orbitals(fname, backend='hdf5'):
    """Reads molecular orbitals from a TREXIO file."""

    trexio_file = None
    try:
        trexio_file = _open_trexio(fname, mode='r', backend=backend)
        assert_msg_critical(_has(trexio_file, 'mo_num'),
                            'TREXIO file does not contain molecular orbitals.')
        nmo = int(_read_if_present(trexio_file, 'mo_num'))
        nao = int(_read_if_present(trexio_file, 'ao_num'))
        molecule = read_molecule(fname, backend)
        basis = _read_basis_from_open_trexio(trexio_file, molecule)
        shells = _basis_shells(molecule, basis)
        trexio_to_vlx, _ = _ao_permutation_to_trexio_order(
            molecule, basis, shells)
        vlx_to_trexio = _inverse_permutation(trexio_to_vlx)

        coefficient = np.array(_read_if_present(trexio_file, 'mo_coefficient'),
                               dtype=float).reshape(nmo, nao)
        energy = np.array(_read_if_present(trexio_file, 'mo_energy',
                                           np.zeros(nmo)),
                          dtype=float)
        occupation = np.array(_read_if_present(trexio_file, 'mo_occupation',
                                               np.zeros(nmo)),
                              dtype=float)
        spin = np.array(_read_if_present(trexio_file, 'mo_spin',
                                         np.zeros(nmo, dtype=int)),
                        dtype=int)

        has_alpha = np.any(spin == 0)
        has_beta = np.any(spin == 1)

        if has_alpha and has_beta:
            alpha = np.where(spin == 0)[0]
            beta = np.where(spin == 1)[0]
            return MolecularOrbitals(
                [
                    coefficient[alpha, :].T[vlx_to_trexio, :],
                    coefficient[beta, :].T[vlx_to_trexio, :]
                ],
                [energy[alpha], energy[beta]],
                [occupation[alpha], occupation[beta]], molorb.unrest)

        half_occupation = 0.5 * occupation
        return MolecularOrbitals([coefficient.T[vlx_to_trexio, :]], [energy],
                                 [half_occupation], molorb.rest)
    finally:
        _safe_close(trexio_file)


def read_trexio(fname, backend='hdf5'):
    """Reads available VeloxChem objects from a TREXIO file.

    :return:
        A dictionary with keys among ``molecule``, ``basis``,
        ``molecular_orbitals`` and ``scf_energy``.
    """

    data = {'molecule': read_molecule(fname, backend)}

    try:
        molecule, basis = read_molecule_and_basis(fname, backend)
        data['molecule'] = molecule
        data['basis'] = basis
    except Exception:
        pass

    try:
        data['molecular_orbitals'] = read_molecular_orbitals(fname, backend)
    except Exception:
        pass

    trexio_file = None
    try:
        trexio_file = _open_trexio(fname, mode='r', backend=backend)
        if _has(trexio_file, 'state_energy'):
            state_energy = np.asarray(
                _read_if_present(trexio_file, 'state_energy'), dtype=float)
            data['scf_energy'] = float(state_energy.reshape(-1)[0])
    except Exception:
        pass
    finally:
        _safe_close(trexio_file)

    return data
