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

import numpy as np

from .errorhandler import assert_msg_critical
from .groupedalgebra import downward_closure, mobius_sign, subsets_of
from .groupedbanks import WEIGHT_SCHEMES, CovariantTransport, IncrementBank
from .groupedcoordinates import DIHEDRAL, wrap_angle
from .groupeddefinitions import DEFAULT_NODE_ORDERS, SUPPORTED_POLICY_VALUES, canonical_group_key
from .groupedrelaxation import (PhaseConstrainedMinimizer, assert_anchor_is_minimum, assert_node_continuity, assert_well_posed,
                                assert_zone_invariance, branch_diagnostics, image_record, investigate_branch)

TWO_PI = 2.0 * np.pi


def cut_key(offsets):
    """The key of a cut record: its phase offsets rounded to 1e-8 rad, as a tuple of ints."""
    return tuple(int(np.rint(float(o) * 1.0e8)) for o in offsets)


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

    # ------------------------------------------------------------------ build

    @classmethod
    def build(cls, coordinates, reference, data_source, retained_terms, node_grids, cut_policy='rigid_rotor_v1',
              response_policy='quadratic_v1', node_orders=None, support_factor=None, relaxation_zone=None,
              ostream=None, energy_gradient=None, weight_scheme=None):
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
        """

        policy = coordinates.policy
        node_orders = dict(DEFAULT_NODE_ORDERS if node_orders is None else node_orders)
        if energy_gradient is None:

            def energy_gradient(x):
                return data_source(x)[:2]

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

        minimizer = None
        frozen = None
        if cut_policy == 'relaxed_rotor_v1':
            assert_anchor_is_minimum(coordinates, energy_gradient, policy)
            zone = list(coordinates.response_rows) if relaxation_zone is None else [int(r) for r in relaxation_zone]
            assert_zone_invariance(coordinates, zone)
            frozen = [int(r) for r in coordinates.response_rows if int(r) not in set(zone)]
            minimizer = PhaseConstrainedMinimizer(coordinates, energy_gradient, frozen_rows=frozen, policy=policy)

        n_g = coordinates.n_groups
        cut_cache = {}

        def node_data_at(offsets):
            key = cut_key(offsets)
            if key in cut_cache:
                return cut_cache[key]
            label = f'the node at offsets {np.degrees(offsets).round(2).tolist()} degrees'
            x_rigid = coordinates.rigid_cut(np.asarray(offsets))
            record = {'offsets': np.array(offsets, dtype=np.float64), 'x_rigid': x_rigid}
            targets = coordinates.phi_anchor + np.asarray(offsets)
            if cut_policy == 'rigid_rotor_v1' or not np.any(np.abs(offsets) > 1.0e-12):
                x = x_rigid
                record['relaxation'] = None
            else:
                result = minimizer.minimize(x_rigid, targets)
                assert_msg_critical(result.converged,
                                    f'grouped MSI model: the relaxation of {label} '
                                    f'did not converge (projected gradient {result.projected_gradient_norm:.2e}).')
                x = result.coordinates
                record['relaxation'] = result
            # the node geometry must be valid for the Z-matrix (M1-C decisions 3 and 5): the rank identities at
            # every node, no bond angle within the margin of 180 degrees
            coordinates.check_rank(x, label=label)
            coordinates.assert_angles_away_from_linear(x, label=label)
            energy, gradient, hessian = data_source(x)
            residual, g_r, h_r, ref = reference.residual_data(x, energy, gradient, hessian)
            w_prime = None
            if cut_policy == 'relaxed_rotor_v1':
                # every relaxed node and the anchor must be a strict constrained minimum (decision 4)
                record['well_posedness'] = assert_well_posed(minimizer, x, hessian, targets, policy, label=label)
                if node_orders.get('w', 0) >= 1:
                    # the valley tangent at every cut, the anchor included (it is
                    # the constrained minimum at its own phases)
                    _, w_prime = minimizer.valley_tangent(x, hessian, targets=targets)
            data = coordinates.extract_node_data(x, residual, g_r, h_r, w_prime=w_prime, label=label)
            record.update({'x': x, 'energy': float(energy), 'residual': residual, 'data': data,
                           'phase_ratio': ref['phase_ratio'], 'diagnostics': data['diagnostics']})
            cut_cache[key] = record
            return record

        # the continuity check of ADR-003 amendment 3 between neighbouring nodes of every singleton grid (the M1-C
        # review found it missing from the build): a soft mode drifting into another basin is caught here
        if cut_policy == 'relaxed_rotor_v1':
            for gid in sorted({g for term in closure for g in term}):
                grid = np.sort(np.mod(np.asarray(node_grids[gid], dtype=np.float64), TWO_PI / coordinates.groups[gid].images))
                for a, b in zip(grid, np.roll(grid, -1)):
                    off_a, off_b = np.zeros(n_g), np.zeros(n_g)
                    off_a[gid], off_b[gid] = a, b
                    assert_node_continuity(coordinates, node_data_at(off_a)['x'], node_data_at(off_b)['x'], policy,
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
                    node_data_at(face_offsets)
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
        model is re-created from its cut records by the same function.
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
                                        weight_scheme=weight_scheme)
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
                rec.setdefault('branch_investigation', []).append(
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
                                weight_scheme=scheme, coefficient_kinds=bank.kinds)
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
        'reference' (the Q_a result). An uncovered query is a critical error
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
            out.update({'energy': None, 'gradient': None, 'hessian': None})
            return out

        n_g = coords.n_groups
        energy = ref['energy']
        gradient = ref['gradient'].reshape(-1).copy() if order >= 1 else None
        hessian = ref['hessian'].copy() if order >= 2 else None
        quadratic = self.response_policy == 'quadratic_v1'

        # bank values; derivatives with respect to the member angles
        bank_out = {}
        for term, bank in self.banks.items():
            bank_out[term] = bank.evaluate(phi[list(term)], order)

        # scalar increments
        for term in self.terms:
            idx = list(term)
            value, first, second = bank_out[term]['I']
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
            out.update({'energy': energy, 'gradient': None if gradient is None else gradient.reshape(self.n_atoms, 3),
                        'hessian': hessian})
            return out

        # valley origin and its phase derivatives
        n_v = len(self.response_rows)
        v_origin = self.q_v_anchor.copy()
        dv = np.zeros((n_g, n_v)) if order >= 1 else None
        d2v = np.zeros((n_g, n_g, n_v)) if order >= 2 else None
        for term in self.terms:
            idx = list(term)
            value, first, second = bank_out[term]['w']
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
        for term in self.terms:
            idx = list(term)
            b_val, b_first, b_second = bank_out[term]['b']
            c_val, c_first, c_second = bank_out[term]['C']
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
        out.update({'energy': energy, 'gradient': None if gradient is None else gradient.reshape(self.n_atoms, 3),
                    'hessian': hessian, 'valley_origin': v_origin, 'response_displacement': delta})
        return out

    def get_energy_and_gradient(self, x):
        out = self.evaluate(x, order=1)
        return out['energy'], out['gradient']
