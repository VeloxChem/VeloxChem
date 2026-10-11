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


"""
Registry of the grouped MSI v2 model family (schema 2.0.0).

This module holds the table of supported policy values, the definition
records of a grouped model (groups, terms, charts, reference, policy), the
canonical keys with integer normalisation, and the model identity hash
(ADR-001, ADR-006 and SPEC S11 of .claude/tasks/grouped-msi). It contains no
numerical code and imports none of the other grouped modules.

Every selector that is stored controls behaviour: a value outside the table
is rejected with assert_msg_critical, never stored as inactive metadata.
"""

from dataclasses import dataclass, field, fields, is_dataclass, asdict
import hashlib
import json

import numpy as np

from .errorhandler import assert_msg_critical

SCHEMA_VERSION = '2.0.0'

# Policy selectors of schema 2.0.0 and their supported values (ADR-001 with
# the amendments of 2026-10-09). The first value of every entry is the
# default.
SUPPORTED_POLICY_VALUES = {
    'schema_version': ('2.0.0',),
    'construction_mode': ('grouped_outer_mobius_v2',),
    'construction_policy': ('outer_core_projected_mobius_v2',),
    'baseline': ('none',),
    'coordinate_convention': ('cartesian_plain',),
    'bond_representation': ('inverse',),
    'reference_policy': ('image_atlas_periodic_rotor_v1',),
    'cut_policy': ('relaxed_rotor_v1', 'rigid_rotor_v1'),
    'response_policy': ('quadratic_v1', 'none'),
    'inner_bank_policy': ('tapered_shepard_taylor_v1',),
    # ADR-004 amendment 12 (2026-10-10, PROMPT-M2 decision D1): the weights of
    # banks with two or more members. tensor_product_v1 is the product of the
    # members' one-dimensional normalised weights (one support radius per
    # member, no anchoring operator); product_metric_v1 is the single
    # product-metric kernel with the anchoring operator, kept for the A/B
    # regression. Singleton banks are the same under both.
    'inner_weight_scheme': ('tensor_product_v1', 'product_metric_v1'),
    'outer_weight_policy': ('legacy_cartesian_core_v1',),
    'uncovered_policy': ('stop',),
}

GROUP_TYPES = ('terminal_rotor', 'fragment_rotor')

# Bank coefficients of a retained term (ADR-004, ADR-007): the scalar
# increment I, the response gradient b, the response curvature C and the
# anchor-relative valley origin w, with the node orders each may carry.
COEFFICIENT_NAMES = ('I', 'b', 'C', 'w')
NODE_ORDER_LIMITS = {'I': (0, 2), 'b': (0, 1), 'C': (0, 0), 'w': (0, 1)}
DEFAULT_NODE_ORDERS = {'I': 2, 'b': 1, 'C': 0, 'w': 1}

# Numerical parameters of the policy record, with their defaults and the
# validation applied to them: (default, kind). kind is one of 'positive',
# 'nonnegative', 'int_ge_2', 'bool', 'degrees_open_60_120', 'fraction'.
NUMERICAL_PARAMETERS = {
    'atlas_delta_degrees': (90.0, 'degrees_open_60_120'),
    'phase_validity_threshold': (0.5, 'fraction'),
    'rotor_deviation_limit_degrees': (90.0, 'positive'),
    'coordinate_seam_margin_degrees': (10.0, 'nonnegative'),
    'inner_taper_exponent': (4, 'int_ge_2'),
    'inner_singular_exponent': (2, 'int_ge_2'),
    'inner_support_factor': (2.5, 'positive'),
    'outer_exponent_p': (6, 'int_ge_2'),
    'outer_exponent_q': (2, 'int_ge_2'),
    'outer_weight_scale_bohr': (0.5, 'positive'),
    'response_domain_energy_max': (0.01, 'positive'),
    'response_domain_validated': (False, 'bool'),
    'relaxation_gradient_threshold': (1.0e-5, 'positive'),
    'relaxation_constraint_threshold': (1.0e-10, 'positive'),
    'relaxation_step_threshold': (1.0e-8, 'positive'),
    'relaxation_max_iterations': (500, 'int_ge_2'),
    'node_continuity_limit_degrees': (20.0, 'positive'),
    'node_continuity_limit_bond_bohr': (0.2, 'positive'),
    # M1-C (PROMPT-M1-CLOSE.md decisions 3 to 5; seeds from the R1 and R3 measurements)
    'split_rank_gap_min': (1.0e6, 'positive'),            # sigma_r / sigma_{r+1} of G = J J^T at every node
    'split_condition_max': (1.0e5, 'positive'),           # sigma_1 / sigma_r of G at every node
    'legacy_cutoff_margin': (100.0, 'positive'),          # sigma_r(BB^T) >= M 1e-6 and sigma_{r+1} <= 1e-6 / M at chart anchors
    'relaxation_curvature_min': (1.0e-5, 'positive'),     # lambda_min(Z^T K Z) at relaxed nodes, hartree / bohr^2
    'linear_angle_margin_degrees': (10.0, 'positive'),    # no bond angle within this margin of 180 degrees
    # M2 (PROMPT-M2 decision D3, ADR-003 amendment 6): branch diagnostics between neighbouring relaxed nodes. They
    # warn and are recorded; they reject nothing (the investigation procedure gives the verdict). Seeds from
    # measurements/m2_branch_diagnostics.py (2026-10-10): the smooth synthetic models of the M1 and M1-C suites (five
    # fixtures at 40/60 deg, ethanol and methyl acetate at 30 and 20 deg on the standard and the mismatched surface)
    # reach 17.3 (secant), 7.1 (energy), 0.154 (mode ratio) and 3.3 (tangent ratio); the ethanol bifurcation of R1.2
    # reaches 10.4 / 4.8 / 0.31 / 6.95 at 20 deg and 2.0 / 5.4 / 0.55 / 8.8 at 30 deg, so only the tangent ratio
    # separates it (factor about 1.5 on each side); the melatonin bifurcation (13.1 / 5.9 / 0.80 / 1.24) is caught
    # by none of them. The seeds keep every smooth model silent
    'branch_secant_tangent_max': (25.0, 'positive'),      # |(w_{j+1} - w_j) - dphi (w'_j + w'_{j+1}) / 2| relative to the trapezoid size
    #                                                       dphi^3 n^2 max|w'| / 12 of a branch with harmonic n (the mismatch of a flip is O(1))
    'branch_energy_multiplier_max': (10.0, 'positive'),   # |(E_{j+1} - E_j) - dphi (lambda_j + lambda_{j+1}) / 2| relative to dphi^3 n^2 max|lambda| / 12
    'branch_tracked_modes': (4, 'int_ge_2'),              # the lowest k eigenpairs of Z^T K Z at the line start, each followed by its best overlap
    'branch_mode_ratio_min': (0.1, 'positive'),           # smallest ratio of a tracked eigenvalue to its value at the line start
    'branch_tangent_ratio_max': (5.0, 'positive'),        # ratio of |w'| at a node to the anchor's |w'| in the same direction (growth of the tangent)
    'branch_sweep_step_degrees': (5.0, 'positive'),       # the step of the bidirectional warm-start sweep of the investigation
    'branch_degenerate_energy_hartree': (1.0e-6, 'positive'),   # mirror branches: the sweeps end in different basins of equal energy
}


def _check_number(key, value, kind):
    if kind == 'bool':
        assert_msg_critical(
            isinstance(value, (bool, np.bool_)),
            f'grouped MSI policy: {key} must be a boolean, got {value!r}.')
        return
    is_number = isinstance(value, (int, float, np.integer, np.floating)) and not isinstance(value, bool)
    assert_msg_critical(is_number, f'grouped MSI policy: {key} must be a number, got {value!r}.')
    assert_msg_critical(np.isfinite(value), f'grouped MSI policy: {key} must be finite, got {value!r}.')
    if kind == 'positive':
        assert_msg_critical(value > 0, f'grouped MSI policy: {key} must be positive, got {value!r}.')
    elif kind == 'nonnegative':
        assert_msg_critical(value >= 0, f'grouped MSI policy: {key} must not be negative, got {value!r}.')
    elif kind == 'int_ge_2':
        assert_msg_critical(
            isinstance(value, (int, np.integer)) and value >= 2,
            f'grouped MSI policy: {key} must be an integer of at least 2, got {value!r}.')
    elif kind == 'degrees_open_60_120':
        assert_msg_critical(
            60.0 < value < 120.0,
            f'grouped MSI policy: {key} must lie strictly between 60 and 120 degrees '
            f'(coverage of the turn and a plateau), got {value!r}.')
    elif kind == 'fraction':
        assert_msg_critical(
            0.0 < value <= 1.0,
            f'grouped MSI policy: {key} must lie in (0, 1], got {value!r}.')
    else:
        assert_msg_critical(False, f'grouped MSI policy: unknown parameter kind {kind}.')


def validate_node_orders(node_orders):
    """
    Checks a mapping coefficient -> node order against NODE_ORDER_LIMITS.
    """

    assert_msg_critical(
        set(node_orders) == set(COEFFICIENT_NAMES),
        'grouped MSI policy: node_orders must have exactly the keys '
        f'{COEFFICIENT_NAMES}, got {sorted(node_orders)}.')
    for name, order in node_orders.items():
        low, high = NODE_ORDER_LIMITS[name]
        assert_msg_critical(
            isinstance(order, (int, np.integer)) and not isinstance(order, bool) and low <= order <= high,
            f'grouped MSI policy: node order of {name} must be an integer in [{low}, {high}], got {order!r}.')


def validate_policy_record(record):
    """
    Validates a policy record (a mapping) against the tables above.

    Unknown keys, unknown selector values and invalid numerical parameters
    are critical errors. Returns the record as a plain dict with normalised
    integers.
    """

    known = set(SUPPORTED_POLICY_VALUES) | set(NUMERICAL_PARAMETERS) | {'node_orders'}
    unknown = sorted(set(record) - known)
    assert_msg_critical(not unknown, f'grouped MSI policy: unknown keys {unknown}.')
    missing = sorted(known - set(record))
    assert_msg_critical(not missing, f'grouped MSI policy: missing keys {missing}.')

    for key, values in SUPPORTED_POLICY_VALUES.items():
        assert_msg_critical(
            record[key] in values,
            f'grouped MSI policy: unsupported value {record[key]!r} for {key}; '
            f'supported: {list(values)}.')
    for key, (_, kind) in NUMERICAL_PARAMETERS.items():
        _check_number(key, record[key], kind)
    validate_node_orders(record['node_orders'])

    out = normalize_for_json(dict(record))
    # a float-kind parameter given as an integer (2 for 2.0) must hash like the float (review of 2026-10-10)
    for key, (_, kind) in NUMERICAL_PARAMETERS.items():
        if kind in ('positive', 'nonnegative', 'degrees_open_60_120', 'fraction'):
            out[key] = float(out[key])
    return out


@dataclass(frozen=True)
class GroupedModelPolicyV2:
    """
    The persisted policy record of a grouped model: every selector of
    SUPPORTED_POLICY_VALUES, the node orders per coefficient, and every
    numerical parameter that affects evaluation. All of it enters the model
    identity hash.
    """

    schema_version: str = SUPPORTED_POLICY_VALUES['schema_version'][0]
    construction_mode: str = SUPPORTED_POLICY_VALUES['construction_mode'][0]
    construction_policy: str = SUPPORTED_POLICY_VALUES['construction_policy'][0]
    baseline: str = SUPPORTED_POLICY_VALUES['baseline'][0]
    coordinate_convention: str = SUPPORTED_POLICY_VALUES['coordinate_convention'][0]
    bond_representation: str = SUPPORTED_POLICY_VALUES['bond_representation'][0]
    reference_policy: str = SUPPORTED_POLICY_VALUES['reference_policy'][0]
    cut_policy: str = SUPPORTED_POLICY_VALUES['cut_policy'][0]
    response_policy: str = SUPPORTED_POLICY_VALUES['response_policy'][0]
    inner_bank_policy: str = SUPPORTED_POLICY_VALUES['inner_bank_policy'][0]
    inner_weight_scheme: str = SUPPORTED_POLICY_VALUES['inner_weight_scheme'][0]
    outer_weight_policy: str = SUPPORTED_POLICY_VALUES['outer_weight_policy'][0]
    uncovered_policy: str = SUPPORTED_POLICY_VALUES['uncovered_policy'][0]
    node_orders: dict = field(default_factory=lambda: dict(DEFAULT_NODE_ORDERS))
    atlas_delta_degrees: float = NUMERICAL_PARAMETERS['atlas_delta_degrees'][0]
    phase_validity_threshold: float = NUMERICAL_PARAMETERS['phase_validity_threshold'][0]
    rotor_deviation_limit_degrees: float = NUMERICAL_PARAMETERS['rotor_deviation_limit_degrees'][0]
    coordinate_seam_margin_degrees: float = NUMERICAL_PARAMETERS['coordinate_seam_margin_degrees'][0]
    inner_taper_exponent: int = NUMERICAL_PARAMETERS['inner_taper_exponent'][0]
    inner_singular_exponent: int = NUMERICAL_PARAMETERS['inner_singular_exponent'][0]
    inner_support_factor: float = NUMERICAL_PARAMETERS['inner_support_factor'][0]
    outer_exponent_p: int = NUMERICAL_PARAMETERS['outer_exponent_p'][0]
    outer_exponent_q: int = NUMERICAL_PARAMETERS['outer_exponent_q'][0]
    outer_weight_scale_bohr: float = NUMERICAL_PARAMETERS['outer_weight_scale_bohr'][0]
    response_domain_energy_max: float = NUMERICAL_PARAMETERS['response_domain_energy_max'][0]
    response_domain_validated: bool = NUMERICAL_PARAMETERS['response_domain_validated'][0]
    relaxation_gradient_threshold: float = NUMERICAL_PARAMETERS['relaxation_gradient_threshold'][0]
    relaxation_constraint_threshold: float = NUMERICAL_PARAMETERS['relaxation_constraint_threshold'][0]
    relaxation_step_threshold: float = NUMERICAL_PARAMETERS['relaxation_step_threshold'][0]
    relaxation_max_iterations: int = NUMERICAL_PARAMETERS['relaxation_max_iterations'][0]
    node_continuity_limit_degrees: float = NUMERICAL_PARAMETERS['node_continuity_limit_degrees'][0]
    node_continuity_limit_bond_bohr: float = NUMERICAL_PARAMETERS['node_continuity_limit_bond_bohr'][0]
    split_rank_gap_min: float = NUMERICAL_PARAMETERS['split_rank_gap_min'][0]
    split_condition_max: float = NUMERICAL_PARAMETERS['split_condition_max'][0]
    legacy_cutoff_margin: float = NUMERICAL_PARAMETERS['legacy_cutoff_margin'][0]
    relaxation_curvature_min: float = NUMERICAL_PARAMETERS['relaxation_curvature_min'][0]
    linear_angle_margin_degrees: float = NUMERICAL_PARAMETERS['linear_angle_margin_degrees'][0]
    branch_secant_tangent_max: float = NUMERICAL_PARAMETERS['branch_secant_tangent_max'][0]
    branch_energy_multiplier_max: float = NUMERICAL_PARAMETERS['branch_energy_multiplier_max'][0]
    branch_tracked_modes: int = NUMERICAL_PARAMETERS['branch_tracked_modes'][0]
    branch_mode_ratio_min: float = NUMERICAL_PARAMETERS['branch_mode_ratio_min'][0]
    branch_tangent_ratio_max: float = NUMERICAL_PARAMETERS['branch_tangent_ratio_max'][0]
    branch_sweep_step_degrees: float = NUMERICAL_PARAMETERS['branch_sweep_step_degrees'][0]
    branch_degenerate_energy_hartree: float = NUMERICAL_PARAMETERS['branch_degenerate_energy_hartree'][0]

    def __post_init__(self):
        validate_policy_record(self.as_record())

    def as_record(self):
        """
        The policy as a plain dict (the persisted form).
        """

        return {f.name: getattr(self, f.name) for f in fields(self)}

    @classmethod
    def from_record(cls, record):
        return cls(**validate_policy_record(record))

    def replaced(self, **changes):
        record = self.as_record()
        record.update(changes)
        return self.from_record(record)


@dataclass(frozen=True)
class GroupDefinition:
    """
    One rotor group (ADR-003 amendment 1). Atom and row indices are 0-based.

    - group_id: position in the model's group list
    - group_type: 'terminal_rotor' or 'fragment_rotor'
    - attachment_atom a (root side), axis_atom c (moving side): the rotor bond
    - moving_atoms: the atoms the cut rotates, sorted
    - images m, harmonic n, sense: ADR-002 (sense is +1 for m == 1)
    - cyclic_order: the moving atoms in the cyclic order that defines the
      permutation generator (empty for m == 1)
    - rotor_rows: indices of the proper-dihedral rows about the rotor bond in
      the flattened Z-matrix; offsets: their o_k in radians
    - domain: (low, high) in radians of the phase relative to the anchor,
      None for the full period
    - relaxation_zone: non-rotor row indices allowed to relax, None for all
    - mm_rotatable: the MM flag of the bond; user_override: the candidate was
      selected although MM does not consider it rotatable, or was excluded
    """

    group_id: int
    group_type: str
    attachment_atom: int
    axis_atom: int
    moving_atoms: tuple
    images: int
    harmonic: int
    sense: int
    cyclic_order: tuple
    rotor_rows: tuple
    offsets: tuple
    domain: tuple = None
    relaxation_zone: tuple = None
    mm_rotatable: bool = True
    user_override: bool = False
    label: str = ''

    def __post_init__(self):
        assert_msg_critical(
            self.group_type in GROUP_TYPES,
            f'grouped MSI: unsupported group type {self.group_type!r}; supported: {list(GROUP_TYPES)}.')
        assert_msg_critical(
            isinstance(self.images, (int, np.integer)) and self.images >= 1,
            f'grouped MSI: images (m) of group {self.group_id} must be a positive integer.')
        assert_msg_critical(
            isinstance(self.harmonic, (int, np.integer)) and self.harmonic >= 1,
            f'grouped MSI: harmonic (n) of group {self.group_id} must be a positive integer.')
        assert_msg_critical(
            self.sense in (1, -1),
            f'grouped MSI: sense of group {self.group_id} must be +1 or -1.')
        assert_msg_critical(
            len(self.rotor_rows) == len(self.offsets) and len(self.rotor_rows) > 0,
            f'grouped MSI: group {self.group_id} needs one offset per rotor row.')
        if self.images > 1:
            assert_msg_critical(
                len(self.cyclic_order) == self.images and set(self.cyclic_order) == set(self.moving_atoms),
                f'grouped MSI: group {self.group_id} with {self.images} images needs its moving '
                'atoms in cyclic order.')
        if self.domain is not None:
            assert_msg_critical(
                len(self.domain) == 2 and self.domain[0] < self.domain[1],
                f'grouped MSI: domain of group {self.group_id} must be (low, high) with low < high.')

    @property
    def rotor_bond(self):
        return (self.attachment_atom, self.axis_atom)

    def as_record(self):
        return normalize_for_json(asdict(self))


@dataclass(frozen=True)
class CoreChartDefinition:
    """
    A core chart: its anchor and the identities of what it was built from.
    """

    chart_id: int
    anchor_payload_hash: str
    level_of_theory_hash: str
    group_ids: tuple
    label: str = ''


@dataclass(frozen=True)
class ReferenceModelDefinition:
    """
    Identifies the exact Q_a of a chart (ADR-002): policy, atlas width and
    the harmonic of every group.
    """

    reference_policy: str
    atlas_delta_degrees: float
    group_harmonics: tuple
    group_images: tuple

    def __post_init__(self):
        assert_msg_critical(
            self.reference_policy in SUPPORTED_POLICY_VALUES['reference_policy'],
            f'grouped MSI: unsupported reference policy {self.reference_policy!r}.')


@dataclass(frozen=True)
class TermKey:
    """
    Canonical key of a retained term: the sorted group ids and the kind.
    """

    groups: tuple
    kind: str = 'increment'

    def __post_init__(self):
        object.__setattr__(self, 'groups', canonical_group_key(self.groups))
        assert_msg_critical(self.kind in ('increment',), f'grouped MSI: unknown term kind {self.kind!r}.')

    @property
    def order(self):
        return len(self.groups)


@dataclass(frozen=True)
class CouplingEdge:
    """
    A screened coupling between groups (SPEC S8), with the persisted
    screening values and the selection decision.
    """

    groups: tuple
    eta_operator: float
    eta_frobenius: float
    selected: bool

    def __post_init__(self):
        object.__setattr__(self, 'groups', canonical_group_key(self.groups))


def canonical_group_key(groups):
    """
    Sorted tuple of plain ints; duplicates are a critical error.
    """

    key = tuple(sorted(int(g) for g in groups))
    assert_msg_critical(len(set(key)) == len(key), f'grouped MSI: duplicate group ids in {key}.')
    return key


def z_matrix_as_int_tuples(z_matrix):
    """
    The Z-matrix dict with every entry as a tuple of plain ints (SPEC S5:
    tuples mix numpy and Python integers).
    """

    return {
        section: tuple(tuple(int(a) for a in row) for row in z_matrix[section])
        for section in ('bonds', 'angles', 'dihedrals', 'impropers')
    }


def normalize_for_json(obj):
    """
    Recursively converts numpy scalars and arrays, tuples, sets and
    dataclasses into plain JSON types with integers normalised to int and
    floats to float. Dict keys become strings (ints by str, tuples by their
    normalised repr).
    """

    if is_dataclass(obj) and not isinstance(obj, type):
        return normalize_for_json(asdict(obj))
    if isinstance(obj, dict):
        out = {}
        for key, value in obj.items():
            nkey = normalize_for_json(key)
            if not isinstance(nkey, str):
                nkey = repr(nkey)
            assert_msg_critical(nkey not in out, f'grouped MSI: duplicate canonical key {nkey!r}.')
            out[nkey] = normalize_for_json(value)
        return out
    if isinstance(obj, (list, tuple)):
        return [normalize_for_json(v) for v in obj]
    if isinstance(obj, (set, frozenset)):
        return sorted(normalize_for_json(v) for v in obj)
    if isinstance(obj, np.ndarray):
        return normalize_for_json(obj.tolist())
    if isinstance(obj, (bool, np.bool_)):
        return bool(obj)
    if isinstance(obj, (int, np.integer)):
        return int(obj)
    if isinstance(obj, (float, np.floating)):
        value = float(obj)
        assert_msg_critical(np.isfinite(value), 'grouped MSI: non-finite value in a canonical record.')
        return value
    if obj is None or isinstance(obj, str):
        return obj
    assert_msg_critical(False, f'grouped MSI: cannot put an object of type {type(obj).__name__} into a canonical record.')


def canonical_json(obj):
    """
    Canonical JSON of a record: sorted keys, no whitespace, integers
    normalised, floats in their shortest round-trip form.
    """

    return json.dumps(normalize_for_json(obj), sort_keys=True, separators=(',', ':'), allow_nan=False)


def identity_hash(obj):
    """
    SHA-256 of the canonical JSON of a record (ADR-006).
    """

    return hashlib.sha256(canonical_json(obj).encode('utf-8')).hexdigest()


def array_payload_hash(array):
    """
    SHA-256 of an array's little-endian bytes together with its shape and
    dtype (ADR-006).
    """

    arr = np.ascontiguousarray(np.asarray(array))
    if arr.dtype.byteorder == '>':
        arr = arr.byteswap().view(arr.dtype.newbyteorder('<'))
    digest = hashlib.sha256()
    digest.update(repr((arr.shape, arr.dtype.str)).encode('utf-8'))
    digest.update(arr.tobytes())
    return digest.hexdigest()


def model_identity(labels, charge, multiplicity, level_of_theory, z_matrix, policy, groups,
                   reference, retained_terms, extra=None):
    """
    The identity record of a model (ADR-006) and its hash. Every argument is
    a plain record or a dataclass; arrays are normalised.
    """

    record = {
        'labels': list(labels),
        'charge': int(charge),
        'multiplicity': int(multiplicity),
        'level_of_theory': level_of_theory,
        'z_matrix': z_matrix_as_int_tuples(z_matrix),
        'policy': policy.as_record() if isinstance(policy, GroupedModelPolicyV2) else policy,
        'groups': [g.as_record() if isinstance(g, GroupDefinition) else g for g in groups],
        'reference': reference,
        'retained_terms': sorted((list(TermKey(t).groups) if not isinstance(t, TermKey) else list(t.groups)
                                  for t in retained_terms), key=lambda t: (len(t), t)),
    }
    if extra is not None:
        record['extra'] = extra
    return normalize_for_json(record), identity_hash(record)
