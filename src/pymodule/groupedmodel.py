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
One-chart assembly of the grouped MSI v2 model family (SPEC S1, S6, S7;
ADR-007; PROMPT-M1.md sections 3.5 to 3.8):

    F_a(x) = Q_a(x) + sum_S [ I_S(phi_S) + delta^T b_S(phi_S) + 1/2 delta^T C_S(phi_S) delta ],
    delta = s_v( q_v(x) - v(phi) ),   v(phi) = q_v^a + sum_S Delta_S w(phi_S),

with energy, Cartesian gradient and Cartesian Hessian, for one core chart
(W = 1; the outer blend comes in M4). GroupedChartModel.build constructs the
node data on the cuts of the chosen policy (rigid or relaxed), forms the
Möbius increments of every coefficient array from exact source cuts, and
creates the banks; GroupedChartModel.evaluate assembles F_a with its
coverage diagnostics.
"""

import time

import numpy as np

from .errorhandler import assert_msg_critical
from .groupedalgebra import downward_closure, mobius_sign, subsets_of
from .groupedbanks import WEIGHT_SCHEMES, CovariantTransport, IncrementBank
from .groupedcoordinates import DIHEDRAL, wrap_angle
from .groupeddefinitions import (DEFAULT_NODE_ORDERS, EQUIVARIANCE_TOLERANCES, SUPPORTED_POLICY_VALUES, SYMMETRY_DEFECT_LEVELS,
                                 SYMMETRY_DEFECT_PHASE_SETS, SYMMETRY_DEFECT_SEED, canonical_group_key, validation_parameters)
from .groupedrelaxation import (PhaseConstrainedMinimizer, assert_anchor_is_minimum, assert_node_continuity, assert_well_posed,
                                assert_zone_invariance, branch_diagnostics, image_record, investigate_branch)
from .veloxchemlib import bohr_in_angstrom, hartree_in_kcalpermol

TWO_PI = 2.0 * np.pi
KCAL = hartree_in_kcalpermol()                  # kcal/mol per hartree
KCAL_PER_ANGSTROM = KCAL / bohr_in_angstrom()   # kcal/mol/A per hartree/bohr

# the covariance kind of every response coefficient under the image permutations (the valley w is a plain row
# permutation here: the anchor shift enters the assembled valley once, through the face of the moved group)
RESPONSE_KINDS = {'b': 'vector', 'C': 'matrix', 'w': 'vector'}


def cut_key(offsets):
    """The key of a cut record: its phase offsets rounded to 1e-8 rad, as a tuple of ints."""
    return tuple(int(np.rint(float(o) * 1.0e8)) for o in offsets)


def term_label(term):
    """'0_1' for the term (0, 1): the key of a term in persisted records."""
    return '_'.join(str(int(g)) for g in term)


class CutBuilder:
    """
    The per-cut work of a construction (M3, ADR-008 item 9): for a vector of
    phase offsets, the cut geometry of the policy (the rigid cut, or its
    relaxation with the phase-constrained minimiser), the validity checks at
    the node (rank, linear angles, well-posedness), the data source call, the
    residual data and the node data of ADR-007, cached by `cut_key`.
    GroupedChartModel.build computes its nodes with it and the coupling
    selection its probe cuts, so that probes and nodes are the same kind of
    record. `records` may seed the cache with records computed before by the
    same policies and data source (the caller's responsibility: the
    construction passes the records of its own campaign); the builder adds
    to that dict. `touched` collects the keys requested since the last reset.

    :param data_source: callable x (N, 3) -> (energy, gradient (N, 3), hessian (3N, 3N)), called once per cut
    :param energy_gradient: callable x (N, 3) -> (energy, gradient (N, 3)) for the relaxations and the anchor check
        (default: the first two results of data_source, which then computes a Hessian at every relaxation step)
    """

    def __init__(self, coordinates, reference, data_source, cut_policy, response_policy, node_orders,
                 energy_gradient=None, relaxation_zone=None, records=None):
        assert_msg_critical(cut_policy in SUPPORTED_POLICY_VALUES['cut_policy'],
                            f'grouped MSI model: unsupported cut policy {cut_policy!r}.')
        self.coordinates = coordinates
        self.reference = reference
        self.data_source = data_source
        self.cut_policy = cut_policy
        self.response_policy = response_policy
        self.node_orders = dict(node_orders)
        self.policy = coordinates.policy
        if energy_gradient is None:

            def energy_gradient(x):
                return data_source(x)[:2]

        self.energy_gradient = energy_gradient
        self.records = {} if records is None else records
        self.touched = set()
        self.minimizer = None
        self.frozen = None
        if cut_policy == 'relaxed_rotor_v1':
            assert_anchor_is_minimum(coordinates, energy_gradient, self.policy)
            zone = list(coordinates.response_rows) if relaxation_zone is None else [int(r) for r in relaxation_zone]
            assert_zone_invariance(coordinates, zone)
            self.frozen = [int(r) for r in coordinates.response_rows if int(r) not in set(zone)]
            self.minimizer = PhaseConstrainedMinimizer(coordinates, energy_gradient, frozen_rows=self.frozen, policy=self.policy)

    def record(self, offsets):
        """The cut record at the phase offsets (radians, one per group), computed once."""

        coordinates, reference, policy = self.coordinates, self.reference, self.policy
        key = cut_key(offsets)
        self.touched.add(key)
        if key in self.records:
            return self.records[key]
        label = f'the node at offsets {np.degrees(offsets).round(2).tolist()} degrees'
        x_rigid = coordinates.rigid_cut(np.asarray(offsets))
        record = {'offsets': np.array(offsets, dtype=np.float64), 'x_rigid': x_rigid}
        targets = coordinates.phi_anchor + np.asarray(offsets)
        if self.cut_policy == 'rigid_rotor_v1' or not np.any(np.abs(offsets) > 1.0e-12):
            x = x_rigid
            record['relaxation'] = None
        else:
            result = self.minimizer.minimize(x_rigid, targets)
            assert_msg_critical(result.converged,
                                f'grouped MSI model: the relaxation of {label} '
                                f'did not converge (projected gradient {result.projected_gradient_norm:.2e}).')
            x = result.coordinates
            record['relaxation'] = result
        # the node geometry must be valid for the Z-matrix (M1-C decisions 3 and 5): the rank identities at
        # every node, no bond angle within the margin of 180 degrees
        coordinates.check_rank(x, label=label)
        coordinates.assert_angles_away_from_linear(x, label=label)
        energy, gradient, hessian = self.data_source(x)
        residual, g_r, h_r, ref = reference.residual_data(x, energy, gradient, hessian)
        w_prime = None
        if self.cut_policy == 'relaxed_rotor_v1':
            # every relaxed node and the anchor must be a strict constrained minimum (decision 4)
            record['well_posedness'] = assert_well_posed(self.minimizer, x, hessian, targets, policy, label=label)
            if self.node_orders.get('w', 0) >= 1:
                # the valley tangent at every cut, the anchor included (it is
                # the constrained minimum at its own phases)
                _, w_prime = self.minimizer.valley_tangent(x, hessian, targets=targets)
        data = coordinates.extract_node_data(x, residual, g_r, h_r, w_prime=w_prime, label=label)
        record.update({'x': x, 'energy': float(energy), 'residual': residual, 'data': data,
                       'phase_ratio': ref['phase_ratio'], 'diagnostics': data['diagnostics']})
        self.records[key] = record
        return record


class GroupedChartModel:
    """
    :param coordinates: the GroupedCoordinateModel
    :param reference: the GroupedReference Q_a
    :param banks: dict canonical term -> IncrementBank with the anchored
        increments of the coefficients ('I', and with the quadratic response
        'b', 'C', 'w')
    :param response_policy: 'quadratic_v1' or 'none'
    :param cut_policy: 'rigid_rotor_v1' or 'relaxed_rotor_v1'
    :param weight_scheme: the inner weight scheme of the banks (ADR-004
        amendment 12); default: the one the banks carry (they must agree)
    :param support_factor: the support factor the banks were built with;
        default: the one the banks carry (they must agree). Both are
        recorded on the model because a persisted model must store the
        values actually used.
    """

    def __init__(self, coordinates, reference, banks, response_policy, cut_policy, node_records=None,
                 weight_scheme=None, support_factor=None, frozen_rows=None, branch_warnings=None, branch_lines=None):
        assert_msg_critical(response_policy in SUPPORTED_POLICY_VALUES['response_policy'],
                            f'grouped MSI model: unsupported response policy {response_policy!r}.')
        assert_msg_critical(cut_policy in SUPPORTED_POLICY_VALUES['cut_policy'],
                            f'grouped MSI model: unsupported cut policy {cut_policy!r}.')
        self.coordinates = coordinates
        self.reference = reference
        self.banks = {canonical_group_key(k): v for k, v in banks.items()}
        self.terms = sorted(self.banks, key=lambda t: (len(t), t))
        closure = set(downward_closure(self.terms)) - {()}
        assert_msg_critical(
            closure == set(self.terms),
            'grouped MSI model: the retained hypergraph must be downward closed (every face of a retained term is retained).')
        self.response_policy = response_policy
        self.cut_policy = cut_policy
        schemes = sorted({bank.weight_scheme for bank in self.banks.values()})
        factors = sorted({float(bank.support_factor) for bank in self.banks.values()})
        assert_msg_critical(len(schemes) == 1 and len(factors) == 1,
                            f'grouped MSI model: the banks disagree on the weight scheme {schemes} or the support factor {factors}.')
        self.weight_scheme = schemes[0] if weight_scheme is None else str(weight_scheme)
        self.support_factor = factors[0] if support_factor is None else float(support_factor)
        assert_msg_critical(self.weight_scheme == schemes[0] and abs(self.support_factor - factors[0]) < 1.0e-12,
                            f'grouped MSI model: the weight scheme {self.weight_scheme!r} and support factor {self.support_factor} '
                            f'do not match the banks ({schemes[0]!r}, {factors[0]}).')
        self.node_records = node_records if node_records is not None else {}
        # relaxed policy: the frozen rows of the relaxation zone, and the branch diagnostics of the build (D3)
        self.frozen_rows = list(frozen_rows) if frozen_rows else []
        self.branch_warnings = list(branch_warnings) if branch_warnings else []
        self.branch_lines = list(branch_lines) if branch_lines else []
        self.branch_investigations = []
        self.node_grids = {}
        self.node_orders = {}
        self.n_atoms = coordinates.n_atoms
        self.response_rows = coordinates.response_rows
        self.q_v_anchor = coordinates.q_anchor[self.response_rows].copy()
        self.periodic = coordinates.ic.kinds[self.response_rows] >= DIHEDRAL
        # M3: the covariant transport of the banks (for the image elements), the response block of the anchor
        # Hessian in the complete coordinates (the energy norm E_v of ADR-005 amendments 1 and 7), and the record
        # of the build-time symmetry defect (ADR-004 amendment 13)
        self.transport = next((bank.transport for bank in self.banks.values() if bank.transport is not None), None)
        if self.transport is None:
            self.transport = CovariantTransport(coordinates)
        self.energy_norm_hessian = None
        if hasattr(reference, 'hessian_anchor') and hasattr(reference, 'gradient_anchor'):
            split = coordinates.complete_transform(coordinates.x_anchor, reference.gradient_anchor, reference.hessian_anchor,
                                                   check=False, label='the anchor')
            self.energy_norm_hessian = split['H_vv']
        self.symmetry_defect_record = None

    # ------------------------------------------------------------------ build

    @classmethod
    def build(cls, coordinates, reference, data_source, retained_terms, node_grids, cut_policy='rigid_rotor_v1',
              response_policy='quadratic_v1', node_orders=None, support_factor=None, relaxation_zone=None,
              ostream=None, energy_gradient=None, weight_scheme=None, cut_records=None, measure_symmetry=True):
        """
        Builds the model from a data source.

        :param data_source: callable x (N, 3) -> (energy, gradient (N, 3), hessian (3N, 3N)),
            called once at every node
        :param energy_gradient: callable x (N, 3) -> (energy, gradient (N, 3)) for the
            relaxation of the cuts and the anchor check (relaxed policy; about a hundred
            calls per node); default: the first two results of data_source, which then
            also computes a Hessian at every step of every relaxation
        :param retained_terms: the retained hypergraph (group id tuples); it
            must be downward closed
        :param node_grids: dict gid -> phase offsets (radians, relative to the
            anchor phase, containing 0) of the singleton grid of every group
            that occurs in a term; pair and higher grids are the products
        :param node_orders: per coefficient (default: the policy table)
        :param support_factor: the bank support as a multiple of the node
            spacing (default: the policy value inner_support_factor); the
            value used is recorded on the model (ADR-004 amendment 12)
        :param relaxation_zone: non-rotor rows allowed to relax (relaxed
            policy; default all)
        :param weight_scheme: the weight scheme of multi-member banks
            (default: the policy value inner_weight_scheme; ADR-004
            amendment 12). With 'tensor_product_v1' the build asserts that
            every multi-member bank vanishes on its anchor lines to rounding.
        :param cut_records: cut records computed before by a CutBuilder with
            the same policies and data source (M3, ADR-008 item 9: the
            singletons and probes of the selection); the nodes among them are
            not recomputed and no relaxation is replayed. The dict is shared
            with the builder, so the caller receives the new records too.
        :param measure_symmetry: measure the symmetry defect of the built
            model (ADR-004 amendment 13, decision D3; default True)
        """

        policy = coordinates.policy
        node_orders = dict(DEFAULT_NODE_ORDERS if node_orders is None else node_orders)
        if response_policy == 'none':
            node_orders = {'I': node_orders['I']}
        support_factor = policy.inner_support_factor if support_factor is None else float(support_factor)
        weight_scheme = policy.inner_weight_scheme if weight_scheme is None else str(weight_scheme)
        assert_msg_critical(weight_scheme in WEIGHT_SCHEMES,
                            f'grouped MSI model: unsupported weight scheme {weight_scheme!r}; supported: {list(WEIGHT_SCHEMES)}.')
        terms = sorted({canonical_group_key(t) for t in retained_terms}, key=lambda t: (len(t), t))
        closure = [t for t in downward_closure(terms) if t]
        for term in closure:
            for gid in term:
                assert_msg_critical(gid in node_grids, f'grouped MSI model: no node grid for group {gid}.')
                grid = np.asarray(node_grids[gid], dtype=np.float64)
                assert_msg_critical(np.any(np.abs(grid) < 1.0e-12), f'grouped MSI model: the grid of group {gid} must contain the anchor (offset 0).')
                # the grid must be uniformly spaced: a bank has one support radius (the factor times its smallest node
                # distance), so an uneven grid shrinks the support and leaves the torus uncovered (M1-C finding). A
                # partial grid (a declared domain, ADR-004 amendment 6) is allowed: its wrap-around gap may exceed
                # the spacing, never fall below it
                period = TWO_PI / coordinates.groups[gid].images
                points = np.sort(np.mod(grid, period))
                if len(points) > 1:
                    gaps = np.diff(points)
                    wrap = points[0] + period - points[-1]
                    assert_msg_critical(
                        np.abs(gaps - gaps[0]).max() < 1.0e-9 and wrap >= gaps[0] - 1.0e-9,
                        f'grouped MSI model: the grid of group {gid} is not uniformly spaced on its circle (gaps '
                        f'{np.degrees(np.append(gaps, wrap)).round(3).tolist()} degrees): the spacing must divide the period '
                        f'{np.degrees(period):g} degrees, or the grid must be a partial grid with a wrap-around gap of at least the spacing.')

        builder = CutBuilder(coordinates, reference, data_source, cut_policy, response_policy, node_orders,
                             energy_gradient=energy_gradient, relaxation_zone=relaxation_zone, records=cut_records)
        frozen = builder.frozen
        n_g = coordinates.n_groups

        # the continuity check of ADR-003 amendment 3 between neighbouring nodes of every singleton grid (the M1-C
        # review found it missing from the build): a soft mode drifting into another basin is caught here
        if cut_policy == 'relaxed_rotor_v1':
            for gid in sorted({g for term in closure for g in term}):
                grid = np.sort(np.mod(np.asarray(node_grids[gid], dtype=np.float64), TWO_PI / coordinates.groups[gid].images))
                for a, b in zip(grid, np.roll(grid, -1)):
                    off_a, off_b = np.zeros(n_g), np.zeros(n_g)
                    off_a[gid], off_b[gid] = a, b
                    assert_node_continuity(coordinates, builder.record(off_a)['x'], builder.record(off_b)['x'], policy,
                                           label=f'of group {gid} at {np.degrees(b):.1f} degrees against its neighbour at {np.degrees(a):.1f}')

        # every cut of every term's product grid is computed here (the faces through the cut cache), then the banks
        # are assembled from the cut records (the same function re-creates banks from stored cut data)
        for term in closure:
            members = list(term)
            grids = [np.asarray(node_grids[gid], dtype=np.float64) for gid in members]
            mesh = np.array(np.meshgrid(*grids, indexing='ij')).reshape(len(members), -1).T
            for node_offsets in mesh:
                full = np.zeros(n_g)
                full[members] = node_offsets
                for face in subsets_of(term):
                    face_offsets = np.zeros(n_g)
                    for gid in face:
                        face_offsets[gid] = full[gid]
                    builder.record(face_offsets)
        # the model's node records are the cuts this build used (seeded records it did not touch, such as the probes
        # of a dropped candidate, stay out of the model and of its branch diagnostics)
        cut_cache = {key: rec for key, rec in builder.records.items() if key in builder.touched}
        banks = cls.banks_from_cuts(coordinates, cut_cache, closure, node_grids, node_orders, policy, support_factor,
                                    weight_scheme)
        # the branch diagnostics of ADR-003 amendment 6 (D3) between neighbouring relaxed nodes along every line:
        # they warn and are recorded in the node records; the investigation procedure gives the verdict
        warnings, lines = [], []
        if cut_policy == 'relaxed_rotor_v1':
            # diagnostics 1 and 4 need the valley tangent, which only a w bank of order 1 computes
            warnings, lines = branch_diagnostics(coordinates, cut_cache, policy,
                                                 with_tangent=(response_policy == 'quadratic_v1' and node_orders.get('w', 0) >= 1))
        if ostream is not None:
            ostream.print_info(f'grouped MSI model: {len(cut_cache)} cut geometries, {len(banks)} banks, '
                               f'cut policy {cut_policy}, response policy {response_policy}, '
                               f'weight scheme {weight_scheme}, support factor {support_factor:g}')
            for warning in warnings:
                ostream.print_warning(warning['message'])
            ostream.flush()
        model = cls(coordinates, reference, banks, response_policy, cut_policy, node_records=cut_cache,
                    weight_scheme=weight_scheme, support_factor=support_factor, frozen_rows=frozen,
                    branch_warnings=warnings, branch_lines=lines)
        model.node_grids = {int(g): np.array(node_grids[g], dtype=np.float64) for g in node_grids}
        model.node_orders = dict(node_orders)
        if measure_symmetry:
            # D3: every build measures its symmetry defect (no QM); closed hypergraphs are asserted exact (D1(a))
            model.symmetry_defect(ostream=ostream)
        return model

    @staticmethod
    def banks_from_cuts(coordinates, cut_records, terms, node_grids, node_orders, policy, support_factor, weight_scheme,
                        transport=None):
        """
        The banks of the retained terms from the cut records (keyed by
        `cut_key` of the offsets, each with the node data of
        `extract_node_data` under 'data'): the Möbius increments of every
        coefficient over the faces of every node of every term's product
        grid, with the faces kept for the covariant coefficients (per-face
        rule), then the IncrementBank of each term with the given scheme and
        support factor. `build` calls it on the cuts it computed; a stored
        model is re-created from its cut records by the same function. Every
        bank carries the magnitude of the cut data its increments were formed
        from (`data_scale`), the scale of its anchor-line check.
        """

        n_g = coordinates.n_groups
        transport = CovariantTransport(coordinates) if transport is None else transport
        terms = sorted({canonical_group_key(t) for t in terms}, key=lambda t: (len(t), t))
        banks = {}
        for term in terms:
            members = list(term)
            grids = [np.asarray(node_grids[gid], dtype=np.float64) for gid in members]
            mesh = np.array(np.meshgrid(*grids, indexing='ij')).reshape(len(members), -1).T
            points = []
            data_list = []
            data_scale = {name: 0.0 for name in node_orders}
            for node_offsets in mesh:
                full = np.zeros(n_g)
                full[members] = node_offsets
                cut_values = {}
                for face in subsets_of(term):
                    face_offsets = np.zeros(n_g)
                    for gid in face:
                        face_offsets[gid] = full[gid]
                    key = cut_key(face_offsets)
                    assert_msg_critical(key in cut_records,
                                        f'grouped MSI model: the cut at offsets {np.degrees(face_offsets).round(2).tolist()} deg '
                                        f'(a face of the term {term}) is missing from the cut records.')
                    cut_values[face] = cut_records[key]['data']
                increment = {}
                for name in node_orders:
                    value = None
                    first = None
                    second = None
                    faces = {}
                    for face in subsets_of(term):
                        sign = mobius_sign(term, face)
                        d = cut_values[face]
                        v, f, s = GroupedChartModel._coefficient_arrays(name, d, members, node_orders[name])
                        data_scale[name] = max(data_scale[name], float(np.abs(v).max()) if np.size(v) else 0.0,
                                               0.0 if f is None else float(np.abs(f).max()),
                                               0.0 if s is None else float(np.abs(s).max()))
                        # the cut function of a face depends on the face's
                        # members only: its derivatives with respect to the
                        # other members of the term are zero
                        in_face = np.array([gid in face for gid in members], dtype=np.float64)
                        if f is not None:
                            f = f * in_face.reshape((len(members),) + (1,) * (f.ndim - 1))
                        if s is not None:
                            s = s * np.outer(in_face, in_face).reshape((len(members), len(members)) + (1,) * (s.ndim - 2))
                        faces[face] = (sign, np.array(v, dtype=np.float64), f, s)
                        value = sign * v if value is None else value + sign * v
                        if f is not None:
                            first = sign * f if first is None else first + sign * f
                        if s is not None:
                            second = sign * s if second is None else second + sign * s
                    if name == 'I' or len(members) == 1:
                        increment[name] = (value, first, second)
                    else:
                        # covariant coefficients of a pair or higher term:
                        # every face transforms under its own members' images
                        increment[name] = {'faces': faces}
                data_list.append(increment)
                points.append(coordinates.phi_anchor[members] + node_offsets)
            member_records = [{'images': coordinates.groups[gid].images, 'sense': coordinates.groups[gid].sense,
                               'anchor': coordinates.phi_anchor[gid], 'domain': coordinates.groups[gid].domain}
                              for gid in members]
            banks[term] = IncrementBank(term, member_records, np.array(points), data_list, node_orders,
                                        transport=transport, support_factor=support_factor,
                                        taper_exponent=policy.inner_taper_exponent,
                                        singular_exponent=policy.inner_singular_exponent,
                                        weight_scheme=weight_scheme, data_scale=data_scale)
            if banks[term].tensor:
                # ADR-004 amendment 12: no anchoring operator; the raw interpolant and its derivatives along the
                # line must vanish on every anchor line to rounding, which holds exactly when the node data on the
                # lines are the increments of the same faces (they are: the cut records are keyed by the offsets)
                banks[term].assert_vanishes_on_anchor_lines()
        return banks

    def investigate_branch_warnings(self, energy_gradient, ostream=None):
        """
        The investigation procedure of ADR-003 amendment 6 for every pair of
        nodes that the branch diagnostics flagged: the bidirectional
        warm-start sweep and its classification. The outcomes are recorded
        on the model (branch_investigations) and in the node records
        ('branch_investigation'). Returns the list of outcomes; the caller
        stops a construction on any 'reject' verdict.
        """

        coords = self.coordinates
        pairs = {}
        for warning in self.branch_warnings:
            key = (warning['gid'], warning['pair'], warning['line'])
            pairs.setdefault(key, {'warning': warning, 'under_resolved': False})
            if warning['diagnostic'] in ('secant_tangent', 'tangent_ratio_anchor'):
                pairs[key]['under_resolved'] = True
        outcomes = []
        for (gid, (off_a, off_b), line), info in pairs.items():
            period = TWO_PI / coords.groups[gid].images
            base_line = np.asarray(info['warning']['base'], dtype=np.float64)
            rec_a = rec_b = None
            for record in self.node_records.values():
                off = np.asarray(record['offsets'], dtype=np.float64)
                base = off.copy()
                base[gid] = 0.0
                if not np.allclose(base, base_line, atol=1.0e-7):
                    continue
                if np.allclose(np.mod(off[gid], period), np.mod(off_a, period), atol=1.0e-7):
                    rec_a = record
                if np.allclose(np.mod(off[gid], period), np.mod(off_b, period), atol=1.0e-7):
                    rec_b = record
            assert_msg_critical(rec_a is not None and rec_b is not None,
                                f'grouped MSI model: the nodes of the branch warning on {line} were not found.')
            if off_b >= period - 1.0e-9:
                # the wrap-around pair: the second node is the first node's virtual image one period up
                rec_b = image_record(coords, rec_b, gid)
            outcome = investigate_branch(coords, energy_gradient, coords.policy, rec_a, rec_b, gid,
                                         under_resolved=info['under_resolved'], frozen_rows=self.frozen_rows)
            outcome['line'] = line
            outcomes.append(outcome)
            for rec in (rec_a, rec_b):
                # on the stored node record (for the wrap-around pair the first node, not its virtual image), as the
                # diagnostics write back (review of 2026-10-12, second pass, SP-4)
                rec.get('image_of', rec).setdefault('branch_investigation', []).append(
                    {'classification': outcome['classification'], 'verdict': outcome['verdict'], 'message': outcome['message']})
            if ostream is not None:
                (ostream.print_warning if outcome['verdict'] == 'reject' else ostream.print_info)(
                    'grouped MSI relaxation: ' + outcome['message'])
                ostream.flush()
        self.branch_investigations.extend(outcomes)
        return outcomes

    def rebanked(self, weight_scheme=None, support_factor=None):
        """
        A model with the same cuts and node data and new banks under another
        weight scheme or support factor (ADR-004 amendment 12): the A/B
        comparison of the schemes without recomputing the cuts, and the way
        banks are re-created from stored node data. Tensor banks are
        asserted to vanish on their anchor lines.
        """

        scheme = self.weight_scheme if weight_scheme is None else str(weight_scheme)
        factor = self.support_factor if support_factor is None else float(support_factor)
        banks = {}
        for term, bank in self.banks.items():
            members = [{'images': int(m), 'sense': int(s), 'anchor': float(a), 'domain': d}
                       for m, s, a, d in zip(bank.images, bank.sense, bank.anchor, bank.domain)]
            new = IncrementBank(term, members, bank.node_points, bank.node_data, bank.node_orders,
                                transport=bank.transport, support_factor=factor,
                                taper_exponent=bank.taper_exponent, singular_exponent=bank.singular_exponent,
                                weight_scheme=scheme, coefficient_kinds=bank.kinds, data_scale=bank.data_scale)
            if new.tensor:
                new.assert_vanishes_on_anchor_lines()
            banks[term] = new
        model = GroupedChartModel(self.coordinates, self.reference, banks, self.response_policy, self.cut_policy,
                                  node_records=self.node_records, weight_scheme=scheme, support_factor=factor,
                                  frozen_rows=self.frozen_rows, branch_warnings=self.branch_warnings,
                                  branch_lines=self.branch_lines)
        model.node_grids = dict(self.node_grids)
        model.node_orders = dict(self.node_orders)
        return model

    @staticmethod
    def _coefficient_arrays(name, data, members, order):
        """
        (value, first, second) of one coefficient at a cut, restricted to the
        member angles of the term and truncated to the node order.
        """

        idx = list(members)
        if name == 'I':
            value = np.array(data['I'], dtype=np.float64)
            first = data['I1'][idx] if order >= 1 else None
            second = data['I2'][np.ix_(idx, idx)] if order >= 2 else None
        elif name == 'b':
            value = data['b']
            first = data['b1'][idx] if order >= 1 else None
            second = None
        elif name == 'C':
            value = data['C']
            first = None
            second = None
        elif name == 'w':
            value = data['w']
            first = data['w1'][idx] if order >= 1 else None
            second = None
        else:
            assert_msg_critical(False, f'grouped MSI model: unknown coefficient {name!r}.')
        return value, first, second

    # --------------------------------------------------------------- evaluate

    def coverage(self, phi, phase_ratio, q=None):
        """
        The coverage checks of one chart: phase validity, the banks'
        domains and supports, and (with q, decision 5 of M1-C) no bond
        angle within the persisted margin of 180 degrees.
        """

        reasons = []
        coords = self.coordinates
        threshold = coords.policy.phase_validity_threshold
        for gid, ratio in enumerate(phase_ratio):
            if ratio < threshold:
                reasons.append(f'group {gid}: |Z|/|K| = {ratio:.3f} below {threshold}')
        if q is not None:
            margin = coords.policy.linear_angle_margin_degrees
            for row, degrees in coords.linear_angle_rows(q):
                reasons.append(f'{coords.describe_row(row)} at {degrees:.2f} deg within {margin:g} deg of 180')
        for term, bank in self.banks.items():
            covered, why = bank.coverage(phi[list(term)])
            reasons.extend(why)
        return not reasons, reasons

    def evaluate(self, x, order=2, enforce_coverage=True):
        """
        F_a at x: a dict with 'energy', 'gradient' (N, 3), 'hessian' (3N, 3N)
        (None above `order`), 'covered', 'coverage_reasons', 'phase' and
        'reference' (the Q_a result); with the quadratic response also
        'valley_origin', 'response_displacement', and the diagnostic energy
        norm 'E_v' = 1/2 delta^T H_vv^a delta of ADR-005 amendments 1 and 7
        with its Cartesian gradient 'E_v_gradient' at order >= 1 (None for
        the response policy 'none'). An uncovered query is a critical error
        unless enforce_coverage is False (analysis use; the result is then
        labelled uncovered).
        """

        coords = self.coordinates
        x = np.asarray(x, dtype=np.float64).reshape(self.n_atoms, 3)
        if order >= 2:
            q, b_mat, b2 = coords.ic.values_b_b2(x)
        else:
            q, b_mat = coords.ic.values_and_b(x)
            b2 = None
        qbb2 = (q, b_mat, b2)
        ref = self.reference.evaluate(x, order, qbb2=qbb2)
        phi, grad_phi, hess_phi, ratio = coords.phases_cartesian(x, order=max(order, 1), qbb2=qbb2)
        covered, reasons = self.coverage(phi, ratio, q)
        if enforce_coverage:
            assert_msg_critical(covered, 'grouped MSI model: the query is not covered: ' + '; '.join(reasons))
        out = {'covered': covered, 'coverage_reasons': reasons, 'phase': phi, 'reference': ref}
        if not covered:
            out.update({'energy': None, 'gradient': None, 'hessian': None, 'E_v': None, 'E_v_gradient': None})
            return out
        outputs = self._bank_outputs(phi, order)
        out.update(self._assemble(x, order, qbb2, ref, (phi, grad_phi, hess_phi), outputs))
        return out

    def _bank_outputs(self, phi, order):
        """
        The bank values and their derivatives with respect to the member
        angles at the phases phi, as a list of (term, dict name -> (value,
        first, second)) in the order of the retained terms.
        """

        return [(term, self.banks[term].evaluate(phi[list(term)], order)) for term in self.terms]

    def _assemble(self, x, order, qbb2, ref, phases, outputs):
        """
        The assembly of F_a from Q_a and coefficient fields given as a list of
        (members, dict name -> (value, first, second)): the retained banks'
        outputs in `evaluate`, or explicit fields in the symmetry attribution
        (ADR-004 amendment 13), where `members` are the groups whose angle
        derivatives `first` and `second` carry (an empty tuple for a constant
        field). Returns the energy, the Cartesian gradient (N, 3) and Hessian
        (3N, 3N) (None above `order`) and, with the quadratic response, the
        valley origin, the response displacement, E_v and its gradient.
        """

        coords = self.coordinates
        q, b_mat, b2 = qbb2
        phi, grad_phi, hess_phi = phases
        n_g = coords.n_groups
        energy = ref['energy']
        gradient = ref['gradient'].reshape(-1).copy() if order >= 1 else None
        hessian = ref['hessian'].copy() if order >= 2 else None
        quadratic = self.response_policy == 'quadratic_v1'

        # scalar increments
        for term, bank_out in outputs:
            idx = list(term)
            value, first, second = bank_out['I']
            energy += float(value)
            if order >= 1:
                for i, g in enumerate(idx):
                    gradient += first[i] * grad_phi[g]
            if order >= 2:
                for i, g in enumerate(idx):
                    hessian += first[i] * hess_phi[g]
                    for j, h in enumerate(idx):
                        hessian += second[i, j] * np.outer(grad_phi[g], grad_phi[h])
        if not quadratic:
            return {'energy': energy, 'gradient': None if gradient is None else gradient.reshape(self.n_atoms, 3),
                    'hessian': hessian, 'E_v': None, 'E_v_gradient': None}

        # valley origin and its phase derivatives
        n_v = len(self.response_rows)
        v_origin = self.q_v_anchor.copy()
        dv = np.zeros((n_g, n_v)) if order >= 1 else None
        d2v = np.zeros((n_g, n_g, n_v)) if order >= 2 else None
        for term, bank_out in outputs:
            idx = list(term)
            value, first, second = bank_out['w']
            v_origin = v_origin + value
            if order >= 1:
                for i, g in enumerate(idx):
                    dv[g] += first[i]
            if order >= 2:
                for i, g in enumerate(idx):
                    for j, h in enumerate(idx):
                        d2v[g, h] += second[i, j]
        d = q[self.response_rows] - v_origin
        d[self.periodic] = wrap_angle(d[self.periodic])
        delta, s1, s2 = coords.ic.displacement(d, rows=self.response_rows)
        b_v = b_mat[self.response_rows]
        if order >= 1:
            jac_d = b_v - np.einsum('gk,gx->kx', dv, grad_phi)   # dd_k/dx
            jac_delta = s1[:, None] * jac_d
        for term, bank_out in outputs:
            idx = list(term)
            b_val, b_first, b_second = bank_out['b']
            c_val, c_first, c_second = bank_out['C']
            f_vec = b_val + c_val @ delta
            energy += float(delta @ b_val + 0.5 * delta @ c_val @ delta)
            if order >= 1:
                gradient += jac_delta.T @ f_vec
                for i, g in enumerate(idx):
                    coeff = delta @ b_first[i] + 0.5 * delta @ c_first[i] @ delta
                    gradient += coeff * grad_phi[g]
            if order >= 2:
                hessian += jac_delta.T @ c_val @ jac_delta
                # sum_k f_k d2 delta_k = sum_k f_k [s2_k Jd_k Jd_k^T + s1_k d2 d_k]
                c_k = f_vec * s1
                hessian += (jac_d.T * (f_vec * s2)) @ jac_d
                hessian += np.einsum('k,kxy->xy', c_k, b2[self.response_rows])
                hessian -= np.einsum('g,gxy->xy', dv @ c_k, hess_phi)
                cross_v = np.einsum('ghk,k->gh', d2v, c_k)
                hessian -= np.einsum('gh,gx,hy->xy', cross_v, grad_phi, grad_phi)
                for i, g in enumerate(idx):
                    vec = jac_delta.T @ (b_first[i] + c_first[i] @ delta)
                    hessian += np.outer(grad_phi[g], vec) + np.outer(vec, grad_phi[g])
                    coeff = delta @ b_first[i] + 0.5 * delta @ c_first[i] @ delta
                    hessian += coeff * hess_phi[g]
                    for j, h in enumerate(idx):
                        coeff2 = delta @ b_second[i, j] + 0.5 * delta @ c_second[i, j] @ delta
                        hessian += coeff2 * np.outer(grad_phi[g], grad_phi[h])
        # the diagnostic energy norm of the response displacement (ADR-005 amendments 1 and 7)
        e_v = None
        e_v_gradient = None
        if self.energy_norm_hessian is not None:
            h_delta = self.energy_norm_hessian @ delta
            e_v = float(0.5 * delta @ h_delta)
            if order >= 1:
                e_v_gradient = (jac_delta.T @ h_delta).reshape(self.n_atoms, 3)
        return {'energy': energy, 'gradient': None if gradient is None else gradient.reshape(self.n_atoms, 3),
                'hessian': hessian, 'valley_origin': v_origin, 'response_displacement': delta,
                'E_v': e_v, 'E_v_gradient': e_v_gradient}

    def get_energy_and_gradient(self, x):
        out = self.evaluate(x, order=1)
        return out['energy'], out['gradient']

    # ------------------------------------------------- M3: hypergraph structure

    def image_groups(self):
        """The ids of the groups with images (m > 1)."""
        return [gid for gid, g in enumerate(self.coordinates.groups) if g.images > 1]

    def mobius_coefficients(self):
        """
        c_T = sum over the retained S containing T of (-1)^(|S| - |T|): the
        weight of face T's interpolant I_T[X_T] in the assembled response
        coefficients under tensor-product weights (Q5 section 2.2; ADR-004
        amendment 13).
        """

        return {t: sum((-1) ** (len(s) - len(t)) for s in self.terms if set(t) <= set(s)) for t in self.terms}

    def is_closed_under_images(self):
        """
        Whether the retained hypergraph is closed under T -> T + {g} for
        every group g with images: the necessary and sufficient condition
        for exact equivariance with any node data (ADR-004 amendment 13,
        decision D1(a)).
        """

        terms = set(self.terms)
        for t in self.terms:
            for g in self.image_groups():
                if g not in t and canonical_group_key(t + (g,)) not in terms:
                    return False
        return True

    def leaking_faces(self):
        """
        The faces T with c_T != 0 that lack a group with images, each with
        the groups whose permutations expose it: [(T, [g, ...])].
        """

        c_t = self.mobius_coefficients()
        out = []
        for t in self.terms:
            if c_t[t] == 0:
                continue
            lacking = [g for g in self.image_groups() if g not in t]
            if lacking:
                out.append((t, lacking))
        return out

    def raw_face_banks(self):
        """
        I_T[X_T] for every retained term T: a bank of the raw cut data of the
        face (no Moebius sum) on T's own node set with the same weights and
        the whole-value covariance of T's members (Q5 section 2.2; ADR-004
        amendment 13). Needs the node records of every node.
        """

        coords = self.coordinates
        out = {}
        for term, bank in self.banks.items():
            members = list(term)
            points, data = [], []
            for z in bank.node_points:
                full = np.zeros(coords.n_groups)
                full[members] = z - coords.phi_anchor[members]
                key = cut_key(full)
                assert_msg_critical(key in self.node_records,
                                    f'grouped MSI model: no cut record for the node of the term {term} at offsets '
                                    f'{np.degrees(full).round(2).tolist()} degrees (raw face banks need the node records).')
                rec = self.node_records[key]['data']
                entry = {}
                for name in bank.names:
                    v, f, s = GroupedChartModel._coefficient_arrays(name, rec, members, bank.node_orders[name])
                    entry[name] = (np.array(v, dtype=np.float64), f, s)
                data.append(entry)
                points.append(np.array(z, dtype=np.float64))
            member_records = [{'images': int(m), 'sense': int(s), 'anchor': float(a), 'domain': d}
                              for m, s, a, d in zip(bank.images, bank.sense, bank.anchor, bank.domain)]
            out[term] = IncrementBank(term, member_records, np.array(points), data, bank.node_orders,
                                      transport=bank.transport, support_factor=bank.support_factor,
                                      taper_exponent=bank.taper_exponent, singular_exponent=bank.singular_exponent,
                                      anchored=False, weight_scheme=bank.weight_scheme, coefficient_kinds=bank.kinds)
        return out

    def element_powers(self, x, permutation):
        """
        For the image element with atom permutation p at x: the data powers
        (sense_g k_g mod m_g, where phi(Px) = phi(x) + k_g 2 pi / m_g), the
        row permutation pv of the response rows (q_v(Px) = q_v(x)[pv]) with
        its valley shift, and the coordinate-layer residual max |phi(Px) -
        phi(x) - k 2 pi / m| in radians.
        """

        coords = self.coordinates
        q, qp = coords.ic.values(x), coords.ic.values(np.asarray(x)[permutation])
        phi, phip = coords.phases(q, 0)[0], coords.phases(qp, 0)[0]
        powers = []
        residual = 0.0
        for g, grp in enumerate(coords.groups):
            period = TWO_PI / grp.images
            dphi = wrap_angle(phip[g] - phi[g])
            k = int(np.rint(dphi / period))
            residual = max(residual, abs(wrap_angle(dphi - k * period)))
            powers.append(int(grp.sense * k) % grp.images)
        pv, shift = self.transport.permutation(tuple(powers))
        return tuple(powers), pv, shift, residual

    # ------------------------------------------- M3: geometries in the domain

    def valley_geometry(self, offsets):
        """
        The model's own valley geometry at the phase offsets: the rigid cut
        realised onto the valley origin v(phi) at fixed phases (the rigid cut
        itself for the response policy 'none'). Returns (x, info) with the
        realisation info of realise_response_displacement (its 'mismatch' is
        the part of the interpolated valley outside the realisable manifold
        of the redundant non-rotor rows), or (None, info) when the phases are
        not covered or the realisation did not converge.
        """

        coords = self.coordinates
        x_rigid = coords.rigid_cut(np.asarray(offsets, dtype=np.float64))
        out = self.evaluate(x_rigid, 0, enforce_coverage=False)
        if not out['covered']:
            return None, {'converged': False, 'reason': 'uncovered'}
        if out.get('valley_origin') is None:
            return x_rigid, {'converged': True, 'step': 0.0, 'mismatch': 0.0, 'phase_mismatch': 0.0, 'substeps': 0}
        q = coords.ic.values(x_rigid)
        d = np.asarray(out['valley_origin'], dtype=np.float64) - q[self.response_rows]
        d[self.periodic] = wrap_angle(d[self.periodic])
        y, info = coords.realise_response_displacement(x_rigid, d)
        return (y if info['converged'] else None), info

    def response_directions(self, levels, e_max, rng):
        """
        One seeded non-rotor displacement per level f > 0 with the energy
        norm 1/2 dq_v^T H_vv^a dq_v = f E_max (ADR-008 item 6): Gaussian in the
        eigenbasis of H_vv^a restricted to its positive eigenvalues. Returns
        a list over the levels (None for f = 0 or without H_vv^a).
        """

        h = self.energy_norm_hessian
        if h is None:
            return [None for _ in levels]
        lam, vec = np.linalg.eigh(0.5 * (h + h.T))
        keep = lam > 1.0e-12 * max(float(lam.max()), 1.0e-300)
        lam, vec = lam[keep], vec[:, keep]
        out = []
        for f in levels:
            if f == 0.0:
                out.append(None)
                continue
            xi = rng.standard_normal(len(lam))
            out.append(vec @ (xi / np.sqrt(lam)) * np.sqrt(2.0 * float(f) * float(e_max) / float(xi @ xi)))
        return out

    # ------------------------------------- M3: the build-time symmetry defect

    @staticmethod
    def _statistic(values):
        """D2's statistic over a set: the 95th percentile (what gates) and the maximum, NaN entries ignored."""
        v = np.asarray(values, dtype=np.float64).reshape(-1)
        v = v[np.isfinite(v)]
        if v.size == 0:
            return {'p95': float('nan'), 'max': float('nan'), 'n': 0}
        return {'p95': float(np.percentile(v, 95)), 'max': float(v.max()), 'n': int(v.size)}

    def symmetry_defect(self, e_max=None, validation=None, phase_sets=None, levels=None, seed=None, ostream=None,
                        hessian_points=2, attribution=True):
        """
        The symmetry defect of the assembled model (ADR-004 amendment 13,
        decisions D1 to D4), measured without QM on a deterministic seeded
        set of geometries: every group at random phases on the model's own
        valley, displaced to E_v = f E_max for f in `levels`, each with all
        its images under the product group of image elements. Records, in
        D2's statistic, |E(Px) - E(x)| in kcal/mol and the largest per-atom
        |g(Px) - P g(x)| in kcal/mol/A per (point, element), their 95th
        percentile and maximum (in total, per level and per element), the
        Moebius coefficients and the leaking faces, and the per-face
        attribution c_T (Pi_{e|T} - Pi_e) I_T[X_T] with its energy and force
        contribution. A hypergraph closed under the image groups is asserted
        exact (D1(a); tensor-product weights). Above the D4 tolerance a
        warning is printed; nothing stops. The record is kept on the model
        (`symmetry_defect_record`) and returned.
        """

        t0 = time.time()
        coords = self.coordinates
        params = validation_parameters() if validation is None else dict(validation)
        tol_e = float(params['symmetry_defect_energy_max_kcal'])
        tol_f = float(params['symmetry_defect_force_max_kcal_per_angstrom'])
        e_max = coords.policy.response_domain_energy_max if e_max is None else float(e_max)
        phase_sets = SYMMETRY_DEFECT_PHASE_SETS if phase_sets is None else int(phase_sets)
        levels = tuple(SYMMETRY_DEFECT_LEVELS) if levels is None else tuple(float(f) for f in levels)
        seed = SYMMETRY_DEFECT_SEED if seed is None else int(seed)
        elements = coords.image_elements()[1:]
        closed = self.is_closed_under_images()
        record = {
            'status': 'measured', 'closed': bool(closed), 'weight_scheme': self.weight_scheme,
            'response_policy': self.response_policy, 'terms': [list(t) for t in self.terms],
            'e_max': float(e_max), 'levels': [float(f) for f in levels], 'phase_sets': int(phase_sets), 'seed': int(seed),
            'elements': [list(int(e) for e in expo) for expo, _ in elements],
            'tolerances': {'energy_kcal': tol_e, 'force_kcal_per_angstrom': tol_f},
            'mobius_coefficients': {term_label(t): int(c) for t, c in self.mobius_coefficients().items()},
            'leaking_faces': [[term_label(t), [int(g) for g in lacking]] for t, lacking in self.leaking_faces()],
            'asserted_exact': False, 'warning': None,
        }
        if not elements:
            record['status'] = 'no image group'
            record['wall_time'] = time.time() - t0
            self.symmetry_defect_record = record
            return record
        rng = np.random.default_rng(seed)
        offsets_all = rng.uniform(0.0, TWO_PI, size=(phase_sets, coords.n_groups))
        directions = self.response_directions(levels, e_max, rng)
        points = []
        dropped = {'uncovered': 0, 'unrealised': 0, 'no_direction': 0}
        valley_mismatch = 0.0
        for s in range(phase_sets):
            x_v, info = self.valley_geometry(offsets_all[s])
            if x_v is None:
                dropped['uncovered' if info.get('reason') == 'uncovered' else 'unrealised'] += 1
                continue
            valley_mismatch = max(valley_mismatch, float(info.get('mismatch', 0.0)))
            for level, f in enumerate(levels):
                if f == 0.0:
                    x = x_v
                else:
                    if directions[level] is None:
                        dropped['no_direction'] += 1
                        continue
                    x, info_d = coords.realise_response_displacement(x_v, directions[level])
                    if not info_d['converged']:
                        dropped['unrealised'] += 1
                        continue
                out = self.evaluate(x, 1, enforce_coverage=False)
                if not out['covered']:
                    dropped['uncovered'] += 1
                    continue
                points.append({'x': x, 'level': float(f), 'E_v': out['E_v'], 'out': out, 'images': []})
        record['dropped'] = dropped
        record['n_points'] = len(points)
        record['valley_mismatch_max'] = valley_mismatch   # the interpolated valley outside the realisable manifold
        if not points:
            record['status'] = 'no points'
            record['wall_time'] = time.time() - t0
            self.symmetry_defect_record = record
            return record
        n_p, n_e = len(points), len(elements)
        de = np.full((n_p, n_e), np.nan)
        df = np.full((n_p, n_e), np.nan)
        for i, pt in enumerate(points):
            for j, (expo, p) in enumerate(elements):
                out_p = self.evaluate(pt['x'][p], 1, enforce_coverage=False)
                pt['images'].append(out_p)
                if not out_p['covered']:
                    dropped['uncovered'] += 1
                    continue
                de[i, j] = abs(out_p['energy'] - pt['out']['energy']) * KCAL
                df[i, j] = float(np.linalg.norm(out_p['gradient'] - pt['out']['gradient'][p], axis=1).max()) * KCAL_PER_ANGSTROM
        level_of_point = np.array([pt['level'] for pt in points])
        record['energy_kcal'] = de
        record['force_kcal_per_angstrom'] = df
        record['level_of_point'] = level_of_point
        record['E_v_of_point'] = np.array([np.nan if pt['E_v'] is None else float(pt['E_v']) for pt in points])
        stat = {'energy': self._statistic(de), 'force': self._statistic(df), 'per_level': {}, 'per_element': {}}
        for f in levels:
            sel = level_of_point == f
            if sel.any():
                stat['per_level'][f'{f:g}'] = {'energy': self._statistic(de[sel]), 'force': self._statistic(df[sel])}
        for j, (expo, _) in enumerate(elements):
            stat['per_element'][','.join(str(int(e)) for e in expo)] = {'energy': self._statistic(de[:, j]),
                                                                        'force': self._statistic(df[:, j])}
        record['statistic'] = stat
        within = bool(stat['energy']['p95'] <= tol_e and stat['force']['p95'] <= tol_f)
        record['within_tolerance'] = within

        # D1(a): a hypergraph closed under the image groups is exactly equivariant (tensor-product weights; the
        # product-metric scheme of the A/B regression is not exact even for the closure, ADR-004 amendment 12)
        if closed and self.weight_scheme == 'tensor_product_v1':
            hess_rel = 0.0
            for pt in points[:int(hessian_points)]:
                o = self.evaluate(pt['x'], 2, enforce_coverage=False)
                scale = max(1.0, float(np.abs(o['hessian']).max()))
                for expo, p in elements:
                    op = self.evaluate(pt['x'][p], 2, enforce_coverage=False)
                    if not op['covered']:
                        continue
                    p3 = (3 * np.asarray(p)[:, None] + np.arange(3)[None, :]).reshape(-1)
                    hess_rel = max(hess_rel, float(np.abs(op['hessian'] - o['hessian'][np.ix_(p3, p3)]).max()) / scale)
            record['hessian_relative_max'] = hess_rel
            record['asserted_exact'] = True
            tol = EQUIVARIANCE_TOLERANCES
            e_abs = float(np.nanmax(de)) / KCAL
            g_abs = float(np.nanmax(df)) / KCAL_PER_ANGSTROM
            assert_msg_critical(
                e_abs <= tol['energy'] and g_abs <= tol['gradient'] and hess_rel <= tol['hessian'],
                f'grouped MSI model: the retained hypergraph {self.terms} is closed under the image groups, but the model is '
                f'not exactly equivariant: energy {e_abs:.1e} hartree, gradient {g_abs:.1e} hartree/bohr, Hessian {hess_rel:.1e} '
                f'relative (D1(a) limits {tol["energy"]:g}, {tol["gradient"]:g}, {tol["hessian"]:g}).')

        # the per-face attribution (Q5 T3; ADR-004 amendment 13), exact for tensor-product weights
        if attribution and not closed and self.weight_scheme == 'tensor_product_v1' and self.response_policy == 'quadratic_v1':
            record['attribution'] = self._symmetry_attribution(points, elements)

        message = None
        if not closed and not within:
            faces = '; '.join(f'{t} under groups {lacking}' for t, lacking in self.leaking_faces())
            message = (f'grouped MSI model: the symmetry defect of the retained hypergraph {self.terms} is '
                       f'{stat["energy"]["p95"]:.3g} kcal/mol and {stat["force"]["p95"]:.3g} kcal/mol/A at the 95th percentile '
                       f'(max {stat["energy"]["max"]:.3g} / {stat["force"]["max"]:.3g}) over {n_p} geometries with E_v <= E_max '
                       f'({e_max:g} hartree) and {n_e} image elements, above the tolerance {tol_e:g} kcal/mol / {tol_f:g} kcal/mol/A '
                       f'(PROMPT-M3 D4); leaking faces: {faces}.')
            if ostream is not None:
                ostream.print_warning(message)
                ostream.flush()
        record['warning'] = message
        record['wall_time'] = time.time() - t0
        self.symmetry_defect_record = record
        return record

    def _symmetry_attribution(self, points, elements):
        """
        The per-face attribution of the defect (ADR-004 amendment 13): for
        every point and element, the coefficient defect D_X = X(Px) - Pi_e X(x)
        of the assembled b, C and w against its prediction sum_T c_T
        (Pi_{e|T} - Pi_e) I_T[X_T] (relative residual per coefficient), the
        energy and force contribution of every leaking (face, element) from
        the response assembled with the permuted fields of x plus that face's
        leak alone, and the remainder against the measured defect.
        """

        coords = self.coordinates
        periodic = self.periodic
        c_t = self.mobius_coefficients()
        raw = self.raw_face_banks()
        names = [n for n in ('b', 'C', 'w') if n in self.banks[self.terms[0]].names]
        coef_res = {n: 0.0 for n in names}
        coef_scale = {n: 1.0e-300 for n in names}
        faces = {}
        perm_check_e, perm_check_g, coord_res = 0.0, 0.0, 0.0
        rem_e, rem_f = 0.0, 0.0
        for pt in points:
            x, out_x = pt['x'], pt['out']
            phi = out_x['phase']
            outputs_x = self._bank_outputs(phi, 1)
            assembled_x = {n: sum(np.asarray(o[n][0], dtype=np.float64) for _, o in outputs_x) for n in names}
            for j, (expo, p) in enumerate(elements):
                out_p = pt['images'][j]
                if not out_p['covered']:
                    continue
                xp = x[p]
                powers, pv, shift, res = self.element_powers(x, p)
                coord_res = max(coord_res, res)
                q_p, b_p = coords.ic.values_and_b(xp)
                qbb2 = (q_p, b_p, None)
                ref_p = self.reference.evaluate(xp, 1, qbb2=qbb2)
                phi_p, grad_p, _, _ = coords.phases_cartesian(xp, order=1, qbb2=qbb2)
                phases_p = (phi_p, grad_p, None)
                assembled_p = {n: sum(np.asarray(o[n][0], dtype=np.float64) for _, o in self._bank_outputs(phi_p, 1)) for n in names}
                measured = {}
                for n in names:
                    if n == 'C':
                        d = assembled_p[n] - assembled_x[n][np.ix_(pv, pv)]
                    elif n == 'w':
                        d = assembled_p[n] - (assembled_x[n][pv] + shift)
                        d[periodic] = wrap_angle(d[periodic])
                    else:
                        d = assembled_p[n] - assembled_x[n][pv]
                    measured[n] = d
                # the permuted fields at Px: every term's outputs at x transformed as a whole by the element, plus
                # the valley shift once (the exactly equivariant counterpart of the model)
                perm_outputs = []
                for term, o in outputs_x:
                    entry = {'I': (o['I'][0], o['I'][1], None)}
                    for n in names:
                        v, f, _ = o[n]
                        tv, tf, _ = self.transport.transform(RESPONSE_KINDS[n], tuple(powers), np.asarray(v, dtype=np.float64),
                                                             None if f is None else np.asarray(f, dtype=np.float64), None)
                        entry[n] = (tv, tf, None)
                    perm_outputs.append((term, entry))
                constant = {'I': (0.0, np.zeros(0), None)}
                for n in names:
                    shape = assembled_x[n].shape
                    constant[n] = (shift.copy() if n == 'w' else np.zeros(shape), np.zeros((0,) + shape), None)
                perm_outputs.append(((), constant))
                base = self._assemble(xp, 1, qbb2, ref_p, phases_p, perm_outputs)
                perm_check_e = max(perm_check_e, abs(base['energy'] - out_x['energy']))
                perm_check_g = max(perm_check_g, float(np.abs(base['gradient'] - out_x['gradient'][p]).max()))
                sum_de = 0.0
                sum_dg = np.zeros_like(base['gradient'])
                predicted = {n: np.zeros_like(measured[n]) for n in names}
                for t in self.terms:
                    rp = tuple(pw if g in t else 0 for g, pw in enumerate(powers))
                    if rp == tuple(powers) or c_t[t] == 0:
                        continue
                    face_out = raw[t].evaluate(phi[list(t)], 1)
                    leak = {'I': (0.0, np.zeros(len(t)), None)}
                    for n in names:
                        v, f, _ = face_out[n]
                        v = np.asarray(v, dtype=np.float64)
                        f = np.asarray(f, dtype=np.float64)
                        mv, mf, _ = self.transport.transform(RESPONSE_KINDS[n], rp, v, f, None)
                        wv, wf, _ = self.transport.transform(RESPONSE_KINDS[n], tuple(powers), v, f, None)
                        leak[n] = (c_t[t] * (mv - wv), c_t[t] * (mf - wf), None)
                        predicted[n] = predicted[n] + leak[n][0]
                    pert = self._assemble(xp, 1, qbb2, ref_p, phases_p, perm_outputs + [(t, leak)])
                    de_t = (pert['energy'] - base['energy']) * KCAL
                    dg_t = pert['gradient'] - base['gradient']
                    key = f'{term_label(t)}|{",".join(str(int(e)) for e in expo)}'
                    rec = faces.setdefault(key, {'face': [int(g) for g in t], 'element': [int(e) for e in expo], 'c_T': int(c_t[t]),
                                                 'energy_kcal_max': 0.0, 'force_kcal_per_angstrom_max': 0.0,
                                                 **{f'leak_{n}_max': 0.0 for n in names}, **{f'face_{n}_max': 0.0 for n in names}})
                    rec['energy_kcal_max'] = max(rec['energy_kcal_max'], abs(de_t))
                    rec['force_kcal_per_angstrom_max'] = max(rec['force_kcal_per_angstrom_max'],
                                                             float(np.linalg.norm(dg_t, axis=1).max()) * KCAL_PER_ANGSTROM)
                    for n in names:
                        rec[f'leak_{n}_max'] = max(rec[f'leak_{n}_max'], float(np.abs(leak[n][0]).max()))
                        rec[f'face_{n}_max'] = max(rec[f'face_{n}_max'], float(np.abs(face_out[n][0]).max()))
                    sum_de += pert['energy'] - base['energy']
                    sum_dg += dg_t
                for n in names:
                    d = measured[n] - predicted[n]
                    if n == 'w':
                        d[periodic] = wrap_angle(d[periodic])
                    coef_res[n] = max(coef_res[n], float(np.abs(d).max()))
                    coef_scale[n] = max(coef_scale[n], float(np.abs(assembled_x[n]).max()))
                rem_e = max(rem_e, abs((out_p['energy'] - out_x['energy']) - sum_de) * KCAL)
                rem_f = max(rem_f, float(np.linalg.norm((out_p['gradient'] - out_x['gradient'][p]) - sum_dg, axis=1).max())
                            * KCAL_PER_ANGSTROM)
        return {'coefficient_residual_relative': {n: coef_res[n] / coef_scale[n] for n in names},
                'coefficient_scale': {n: coef_scale[n] for n in names},
                'permuted_fields_check': {'energy_hartree': perm_check_e, 'gradient_hartree_per_bohr': perm_check_g},
                'coordinate_layer_residual_rad': coord_res, 'faces': faces,
                'remainder_energy_kcal_max': rem_e, 'remainder_force_kcal_per_angstrom_max': rem_f}
