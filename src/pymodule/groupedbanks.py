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
Inner increment banks of the grouped MSI v2 model family: policy
`tapered_shepard_taylor_v1` (SPEC S7; ADR-004 with its amendments of
2026-10-09 and amendment 12 of 2026-10-10; PROMPT-M1.md section 3.7).

A bank interpolates the Möbius increment of one retained term S over the
symmetrised rotor angles z = (phi_g)_{g in S}. Node j at z_j carries, per
coefficient, a value and its first and second derivatives in z up to the
coefficient's node order: a local Taylor model t_j(z). The bank is the
Shepard blend of the local models, I_S(z) = sum_j W_j(z) t_j(z), with
normalised weights built from the tapered singular kernel

    d^2 = 2 [1 - cos(m Delta)] / m^2 per angle,   u = d^2 / r^2,
    w = (1 - u)^a / u^p  for u < 1 (a = 4, p = 2),

evaluated in the product form W_j ~ (1 - u_j)^a prod_{k != j} u_k^p, which is
regular at a node and gives 1 - W_j = O(d^4) there.

Two weight schemes exist for banks with two or more members (policy value
`inner_weight_scheme`, ADR-004 amendment 12; singleton banks are the same
under both):

- `tensor_product_v1` (default): W_j(z) = prod_i W^(i)_{j_i}(z_i), the
  product of the members' one-dimensional normalised weights over their own
  singleton grids, with one support radius per member, r_i = c x spacing_i
  < 2 / m_i. A function of one member is then interpolated exactly as that
  member's singleton bank interpolates it (the other members' weights sum
  to one), which makes the assembled model exactly equivariant between
  nodes with relaxed cuts. On an anchor line of member i the only nodes
  with a non-zero weight are those on the line, where the increment and
  its derivatives along the line vanish, so the raw interpolant vanishes
  there with its along-line derivatives and no anchoring operator is
  applied; `GroupedChartModel.build` asserts this to rounding
  (anchor_line_defect). Coverage: every member has a node in its
  one-dimensional support and lies in its domain.
- `product_metric_v1`: the single kernel in the product metric
  d_j^2 = sum_i d_i^2 with one support radius (c x the smallest node
  distance), and the anchoring operator of groupedalgebra applied to the
  interpolant (the M1 formula, kept for the A/B regression of M2).

Coefficients are arrays of any shape. The scalar increment I is invariant
under the image permutations; the array coefficients b, C and w are
covariant: the virtual image of a node k steps up in angle carries the
permuted arrays (ADR-004 amendment 3). The covariance rule holds for the
cut values of every face of a term: a face T transforms under the
permutations of its own members only, and the increment of a term with two
or more members is the signed Möbius sum of its faces' transformed cut
values (the valley shift then enters exactly once, through the singleton
face). Node data of covariant coefficients are therefore stored per face,
{'faces': {T: (sign, value, first, second)}}, and an invariant coefficient as
the plain increment (value, first, second). The metric keeps its period
2 pi / m and the image index that reduces Delta_g to (-pi/m, pi/m] selects
the transformed node data, so a node and its images never coincide.

Face closure of the node set of every multi-member bank is asserted at
creation; the tensor scheme also asserts that the node set is the full
product of the member grids (the partition of unity needs it).
"""

import numpy as np

from .errorhandler import assert_msg_critical
from .groupedalgebra import anchoring_operator, assert_face_closed, subsets_of, restricted_point
from .groupeddefinitions import NODE_ORDER_LIMITS, SUPPORTED_POLICY_VALUES, canonical_group_key

TWO_PI = 2.0 * np.pi

# how the virtual image of a node transforms each coefficient kind
COEFFICIENT_KINDS = {'I': 'invariant', 'b': 'vector', 'C': 'matrix', 'w': 'valley'}

WEIGHT_SCHEMES = SUPPORTED_POLICY_VALUES['inner_weight_scheme']


def reduce_angle(delta, images):
    """
    Reduces an angle difference to [-pi/m, pi/m) and returns it with the
    number of image steps (multiples of 2 pi / m) that were removed (the
    boundary itself has measure zero and zero weight).
    """

    period = TWO_PI / images
    reduced = (delta + 0.5 * period) % period - 0.5 * period
    steps = int(np.rint((delta - reduced) / period))
    return reduced, steps


def metric_distance(delta, images):
    """
    The periodic metric distance sqrt(2 (1 - cos m Delta)) / m of one angle
    difference (2 / m at most, at Delta = pi / m).
    """

    return float(np.sqrt(max(2.0 * (1.0 - np.cos(images * delta)), 0.0)) / images)


class CovariantTransport:
    """
    The transformation of covariant node data under the image permutations
    of the groups (ADR-004 amendment 3), built from a GroupedCoordinateModel:
    for the non-rotor rows, q_v(P x) = q_v(x)[pv], a vector carries b[pv], a
    matrix C[pv, pv], and the anchor-relative valley origin w[pv] + shift
    with shift = wrap(q_v^a[pv] - q_v^a).
    """

    def __init__(self, model):
        self.model = model
        rows = model.response_rows
        self.n_rows = len(rows)
        position = {int(r): i for i, r in enumerate(rows)}
        self._generator = []
        for gid, g in enumerate(model.groups):
            if g.images == 1:
                self._generator.append(None)
                continue
            full = model.row_permutation(gid, 1)
            pv = np.array([position[int(full[r])] for r in rows], dtype=int)
            self._generator.append(pv)
        periodic = model.ic.kinds[rows] >= 2
        self._periodic = periodic
        self._q_anchor = model.q_anchor[rows].copy()
        self._cache = {}

    def permutation(self, powers):
        """
        pv for the product of the generators raised to `powers` (one entry
        per group); the valley shift goes with it.
        """

        key = tuple(int(p) % max(self.model.groups[gid].images, 1) for gid, p in enumerate(powers))
        if key in self._cache:
            return self._cache[key]
        pv = np.arange(self.n_rows)
        for gid, power in enumerate(key):
            gen = self._generator[gid]
            if gen is None or power == 0:
                continue
            for _ in range(power):
                pv = pv[gen]
        shift = self._q_anchor[pv] - self._q_anchor
        shift[self._periodic] = (shift[self._periodic] + np.pi) % TWO_PI - np.pi
        self._cache[key] = (pv, shift)
        return pv, shift

    def transform(self, kind, powers, value, first=None, second=None):
        """
        The node data of the virtual image. Derivatives transform like the
        value, without the shift.
        """

        if kind == 'invariant' or not any(powers):
            return value, first, second
        pv, shift = self.permutation(powers)
        if kind == 'vector' or kind == 'valley':
            out_value = value[..., pv]
            if kind == 'valley':
                out_value = out_value + shift
            out_first = None if first is None else first[..., pv]
            out_second = None if second is None else second[..., pv]
            return out_value, out_first, out_second
        if kind == 'matrix':
            out_value = value[..., pv, :][..., :, pv]
            out_first = None if first is None else first[..., pv, :][..., :, pv]
            out_second = None if second is None else second[..., pv, :][..., :, pv]
            return out_value, out_first, out_second
        assert_msg_critical(False, f'grouped MSI banks: unknown coefficient kind {kind!r}.')


class IncrementBank:
    """
    The bank of one retained term.

    :param term: the group ids of S
    :param members: per member, a dict with 'images' m, 'sense', 'anchor'
        (anchor phase), 'domain' ((low, high) relative to the anchor or None)
    :param node_points: (n_nodes, |S|) node phases (absolute)
    :param node_data: list over nodes of dicts name -> (value, first, second),
        first of shape (|S|,) + value.shape, second (|S|, |S|) + value.shape;
        None where not stored; a covariant coefficient of a term with two or
        more members is given per face as {'faces': {T: (sign, value, first,
        second)}} (see the module docstring)
    :param node_orders: name -> order used (0, 1, 2)
    :param transport: a CovariantTransport, or None for invariant data only
    :param support_factor: the support radius as a multiple of the node
        spacing in the metric (per member with the tensor scheme, the
        smallest node distance of the whole grid with the product-metric
        scheme)
    :param anchored: apply the anchoring operator for |S| >= 2 (product-metric
        scheme only; the tensor scheme never applies it)
    :param weight_scheme: 'tensor_product_v1' (default) or
        'product_metric_v1' (ADR-004 amendment 12); without effect for |S| = 1
    """

    def __init__(self, term, members, node_points, node_data, node_orders, transport=None,
                 support_factor=2.5, taper_exponent=4, singular_exponent=2, anchored=True,
                 coefficient_kinds=None, weight_scheme=None, data_scale=None):
        self.term = canonical_group_key(term)
        # the magnitude of the cut data the increments were formed from, per coefficient (M3): the anchor-line
        # defect is relative to it, so that a bank whose increments vanish to rounding (a term the surface does not
        # have, or the anchor-line banks of the coupling selection) is judged against its data, not its noise
        self.data_scale = {} if data_scale is None else {k: float(v) for k, v in data_scale.items()}
        self.dim = len(self.term)
        assert_msg_critical(len(members) == self.dim, 'grouped MSI banks: one member record per group of the term.')
        self.images = np.array([int(m['images']) for m in members], dtype=int)
        self.sense = np.array([int(m.get('sense', 1)) for m in members], dtype=int)
        self.anchor = np.array([float(m['anchor']) for m in members])
        self.domain = [m.get('domain') for m in members]
        self.node_points = np.asarray(node_points, dtype=np.float64).reshape(-1, self.dim)
        self.n_nodes = self.node_points.shape[0]
        assert_msg_critical(self.n_nodes > 0, 'grouped MSI banks: a bank needs at least one node.')
        assert_msg_critical(len(node_data) == self.n_nodes, 'grouped MSI banks: one data record per node.')
        self.node_data = node_data
        self.names = tuple(node_data[0].keys())
        self.kinds = dict(COEFFICIENT_KINDS)
        if coefficient_kinds:
            self.kinds.update(coefficient_kinds)
        for name in self.names:
            assert_msg_critical(name in self.kinds, f'grouped MSI banks: no coefficient kind for {name!r}.')
        self.node_orders = {}
        for name in self.names:
            order = int(node_orders.get(name, 0))
            limits = NODE_ORDER_LIMITS.get(name, (0, 2))
            assert_msg_critical(limits[0] <= order <= limits[1],
                                f'grouped MSI banks: node order {order} of {name!r} outside {limits}.')
            self.node_orders[name] = order
        self.transport = transport
        self.taper_exponent = int(taper_exponent)
        self.singular_exponent = int(singular_exponent)
        assert_msg_critical(self.taper_exponent >= 2 and self.singular_exponent >= 2,
                            'grouped MSI banks: the taper and singular exponents must be at least 2.')
        self.weight_scheme = WEIGHT_SCHEMES[0] if weight_scheme is None else str(weight_scheme)
        assert_msg_critical(self.weight_scheme in WEIGHT_SCHEMES,
                            f'grouped MSI banks: unsupported weight scheme {self.weight_scheme!r}; supported: {list(WEIGHT_SCHEMES)}.')
        self.support_factor = float(support_factor)
        assert_msg_critical(self.support_factor > 0.0, 'grouped MSI banks: the support factor must be positive.')
        # the tensor-product scheme acts on multi-member banks only; a singleton bank is the same under both
        self.tensor = self.weight_scheme == 'tensor_product_v1' and self.dim >= 2
        self.anchored = bool(anchored) and self.dim >= 2 and not self.tensor
        self.periods = [TWO_PI / m for m in self.images]
        if self.dim >= 2:
            assert_face_closed(self.node_points, self.term, self.anchor, self.periods)
        if self.tensor:
            self.support_radius = None
            self._prepare_members()
        else:
            self.support_radius = self._support_radius(self.support_factor)
            self.member_nodes = None
            self.member_spacing = None
            self.member_radius = None
            self.node_index = None

    # ---------------------------------------------------------------- metric

    def _reduced(self, z):
        """Delta (n_nodes, dim) reduced per dimension, and the image steps."""
        delta = np.zeros((self.n_nodes, self.dim))
        steps = np.zeros((self.n_nodes, self.dim), dtype=int)
        for j in range(self.n_nodes):
            for i in range(self.dim):
                delta[j, i], steps[j, i] = reduce_angle(z[i] - self.node_points[j, i], self.images[i])
        return delta, steps

    def _support_radius(self, factor):
        spacing = None
        for j in range(self.n_nodes):
            delta, _ = self._reduced(self.node_points[j])
            d2 = (2.0 * (1.0 - np.cos(self.images * delta)) / self.images ** 2).sum(axis=1)
            d2[j] = np.inf
            if self.n_nodes > 1:
                spacing = min(spacing, np.sqrt(d2.min())) if spacing is not None else float(np.sqrt(d2.min()))
        if spacing is None or not np.isfinite(spacing) or spacing <= 0.0:
            spacing = min(2.0 / m for m in self.images) / factor * 0.999
        radius = factor * spacing
        bound = min(2.0 / m for m in self.images)
        assert_msg_critical(
            radius < bound,
            f'grouped MSI banks: the support radius {radius:.4f} of the term {self.term} reaches the '
            f'wrap of the metric ({bound:.4f}); use a smaller support factor or a finer grid.')
        return radius

    def _prepare_members(self):
        """
        Tensor scheme: the member grids (the distinct node values of every
        member modulo its period, in node order), the member index of every
        node, the member spacings (the smallest metric distance between two
        member nodes) and the per-member support radii c x spacing_i, each
        asserted below the wrap 2 / m_i. The node set must be the full
        product of the member grids.
        """

        self.member_nodes = []
        self.member_spacing = []
        self.member_radius = []
        index = np.zeros((self.n_nodes, self.dim), dtype=int)
        for i in range(self.dim):
            m = int(self.images[i])
            values = []
            for j in range(self.n_nodes):
                z = float(self.node_points[j, i])
                found = None
                for a, v in enumerate(values):
                    if abs(reduce_angle(z - v, m)[0]) < 1.0e-9:
                        found = a
                        break
                if found is None:
                    values.append(z)
                    found = len(values) - 1
                index[j, i] = found
            values = np.array(values, dtype=np.float64)
            if len(values) > 1:
                spacing = min(metric_distance(values[a] - values[b], m)
                              for a in range(len(values)) for b in range(len(values)) if a != b)
            else:
                spacing = (2.0 / m) / self.support_factor * 0.999
            radius = self.support_factor * spacing
            assert_msg_critical(
                radius < 2.0 / m,
                f'grouped MSI banks: the support radius {radius:.4f} of member {i} (group {self.term[i]}) of the term '
                f'{self.term} reaches the wrap of its metric ({2.0 / m:.4f}); use a smaller support factor or a finer grid.')
            self.member_nodes.append(values)
            self.member_spacing.append(float(spacing))
            self.member_radius.append(float(radius))
        counts = [len(v) for v in self.member_nodes]
        distinct = len({tuple(row) for row in index.tolist()})
        assert_msg_critical(
            distinct == self.n_nodes and self.n_nodes == int(np.prod(counts)),
            f'grouped MSI banks: the node set of the term {self.term} is not the full product of its member grids '
            f'({self.n_nodes} nodes, member grids of {counts} points): the tensor-product weights need a product grid.')
        self.node_index = index

    @property
    def radii(self):
        """The support radii actually used: one per member (tensor scheme) or the single radius."""
        return list(self.member_radius) if self.tensor else [self.support_radius]

    def node_distance(self, z):
        """u_j for every node at the phases z (product-metric scheme)."""
        assert_msg_critical(not self.tensor, 'grouped MSI banks: node_distance is the product-metric quantity; '
                            'the tensor scheme has one distance per member (member_weights).')
        delta, _ = self._reduced(np.asarray(z, dtype=np.float64))
        d2 = (2.0 * (1.0 - np.cos(self.images * delta)) / self.images ** 2).sum(axis=1)
        return d2 / self.support_radius ** 2

    # -------------------------------------------------------------- coverage

    def in_domain(self, z):
        for i in range(self.dim):
            if self.domain[i] is None:
                continue
            low, high = self.domain[i]
            rel = (z[i] - self.anchor[i] + np.pi) % TWO_PI - np.pi
            if not (low - 1.0e-12 <= rel <= high + 1.0e-12):
                return False
        return True

    def coverage(self, z):
        """
        Whether the bank covers z: inside every member's domain, and (tensor
        scheme) every member has a node in its one-dimensional support, or
        (product-metric scheme) at least one node in support, and the same
        at every restricted evaluation point of the anchoring operator.
        Returns (covered, reasons).
        """

        z = np.asarray(z, dtype=np.float64)
        reasons = []
        if not self.in_domain(z):
            reasons.append(f'term {self.term}: outside the declared domain')
        if self.tensor:
            for i in range(self.dim):
                if not self.member_weights(i, z[i], 0)[0]:
                    reasons.append(f'term {self.term}: member {i} (group {self.term[i]}) has no node in support at {z[i]:.4f}')
            return not reasons, reasons
        points = [z]
        if self.anchored:
            points = [restricted_point(z, self.term, face, self.anchor) for face in subsets_of(self.term)]
        for point in points:
            if not np.any(self.node_distance(point) < 1.0):
                reasons.append(f'term {self.term}: no node in support at {np.round(point, 4).tolist()}')
        return not reasons, reasons

    # -------------------------------------------------------------- evaluate

    def _taylor(self, j, dl, steps_j, order):
        """
        The local Taylor models of every coefficient of node j at the reduced
        displacement dl (dim,), with the node data of the virtual image given
        by the image steps: dict name -> (t, dt (dim,)+A, d2t (dim, dim)+A).
        """

        powers = tuple(int(self.sense[i] * steps_j[i]) for i in range(self.dim))
        record = {}
        for name in self.names:
            value, first, second = self._node_arrays(j, name, powers)
            value = np.asarray(value, dtype=np.float64)
            n_order = self.node_orders[name]
            t_val = value.copy()
            t_first = np.zeros((self.dim,) + value.shape) if order >= 1 else None
            t_second = np.zeros((self.dim, self.dim) + value.shape) if order >= 2 else None
            if n_order >= 1 and first is not None:
                first = np.asarray(first, dtype=np.float64)
                t_val = t_val + np.tensordot(dl, first, axes=(0, 0))
                if order >= 1:
                    t_first = t_first + first
            if n_order >= 2 and second is not None:
                second = np.asarray(second, dtype=np.float64)
                t_val = t_val + 0.5 * np.tensordot(dl, np.tensordot(dl, second, axes=(0, 0)), axes=(0, 0))
                if order >= 1:
                    t_first = t_first + np.tensordot(dl, second, axes=(0, 0))
                if order >= 2:
                    t_second = t_second + second
            record[name] = (t_val, t_first, t_second)
        return record

    def _local_models(self, z, order):
        """
        Product-metric scheme. For every node in support: W-independent
        local model values and derivatives of every coefficient at z, the
        weight factors and their derivatives. Returns the support list and
        the arrays.
        """

        delta, steps = self._reduced(z)
        d2 = (2.0 * (1.0 - np.cos(self.images * delta)) / self.images ** 2).sum(axis=1)
        u = d2 / self.support_radius ** 2
        support = [j for j in range(self.n_nodes) if u[j] < 1.0]
        if not support:
            return support, None
        r2 = self.support_radius ** 2
        locals_ = {}
        du = {}
        d2u = {}
        for j in support:
            dl = delta[j]
            du[j] = 2.0 * np.sin(self.images * dl) / self.images / r2 if order >= 1 else None
            d2u[j] = np.diag(2.0 * np.cos(self.images * dl) / r2) if order >= 2 else None
            locals_[j] = self._taylor(j, dl, steps[j], order)
        return support, (u, du, d2u, locals_)

    def _full_powers(self, powers, face=None):
        """
        powers over all groups of the model (zero outside the term, and
        zero outside `face` when a face is given).
        """

        if self.transport is None:
            return powers
        full = [0] * self.transport.model.n_groups
        for i, gid in enumerate(self.term):
            if face is None or gid in face:
                full[gid] = powers[i]
        return tuple(full)

    def _node_arrays(self, j, name, powers):
        """
        (value, first, second) of coefficient `name` at the virtual image of
        node j given by `powers`: the plain increment transformed as a whole
        (invariant data, or a single-member term), or the signed sum of the
        faces' cut values, each transformed by its own members' permutations.
        """

        entry = self.node_data[j][name]
        kind = self.kinds[name]
        if isinstance(entry, dict):
            value = None
            first = None
            second = None
            for face, (sign, v, f, s) in entry['faces'].items():
                if self.transport is not None and kind != 'invariant' and any(powers):
                    v, f, s = self.transport.transform(kind, self._full_powers(powers, face), v, f, s)
                v = np.asarray(v, dtype=np.float64)
                value = sign * v if value is None else value + sign * v
                if f is not None:
                    f = np.asarray(f, dtype=np.float64)
                    first = sign * f if first is None else first + sign * f
                if s is not None:
                    s = np.asarray(s, dtype=np.float64)
                    second = sign * s if second is None else second + sign * s
            return value, first, second
        value, first, second = entry
        if self.transport is not None and kind != 'invariant' and any(powers):
            value, first, second = self.transport.transform(kind, self._full_powers(powers), value, first, second)
        return value, first, second

    def _value_shape(self, name):
        entry = self.node_data[0][name]
        if isinstance(entry, dict):
            return np.asarray(next(iter(entry['faces'].values()))[1]).shape
        return np.asarray(entry[0]).shape

    def _product_form_weights(self, support, u, du, d2u, order, dim):
        """
        The product-form weights hat w_j = (1 - u_j)^a prod_{k != j} u_k^p over
        the nodes in support and their first and second derivatives in a
        `dim`-dimensional argument (du[j] of shape (dim,), d2u[j] of shape
        (dim, dim)), normalised: returns W (dict), dW (dict of (dim,)), d2W
        (dict of (dim, dim)). The one kernel of both weight schemes.
        """

        a_exp, p_exp = self.taper_exponent, self.singular_exponent
        factors = {}
        for j in support:
            # factor list for hat w_j: (value, first, second) of (1-u_j)^a and of u_k^p, k != j
            lst = []
            one = 1.0 - u[j]
            f0 = one ** a_exp
            f1 = -a_exp * one ** (a_exp - 1) * du[j] if order >= 1 else None
            f2 = None
            if order >= 2:
                f2 = a_exp * (a_exp - 1) * one ** (a_exp - 2) * np.outer(du[j], du[j]) - a_exp * one ** (a_exp - 1) * d2u[j]
            lst.append((f0, f1, f2))
            for k in support:
                if k == j:
                    continue
                g0 = u[k] ** p_exp
                g1 = p_exp * u[k] ** (p_exp - 1) * du[k] if order >= 1 else None
                g2 = None
                if order >= 2:
                    g2 = p_exp * (p_exp - 1) * u[k] ** (p_exp - 2) * np.outer(du[k], du[k]) + p_exp * u[k] ** (p_exp - 1) * d2u[k]
                lst.append((g0, g1, g2))
            factors[j] = lst
        hat = {}
        dhat = {}
        d2hat = {}
        for j in support:
            lst = factors[j]
            n = len(lst)
            values = [f[0] for f in lst]
            prod_all = float(np.prod(values))
            hat[j] = prod_all
            if order >= 1:
                leave_one = [float(np.prod(values[:i] + values[i + 1:])) if n > 1 else 1.0 for i in range(n)]
                dhat[j] = sum(leave_one[i] * lst[i][1] for i in range(n))
            if order >= 2:
                d2 = np.zeros((dim, dim))
                for i in range(n):
                    d2 += leave_one[i] * lst[i][2]
                    for k in range(n):
                        if k == i:
                            continue
                        leave_two = float(np.prod([values[m] for m in range(n) if m not in (i, k)])) if n > 2 else 1.0
                        d2 += leave_two * np.outer(lst[i][1], lst[k][1])
                d2hat[j] = d2
        total = sum(hat.values())
        assert_msg_critical(total > 0.0, 'grouped MSI banks: zero total weight inside the support.')
        w_norm = {j: hat[j] / total for j in support}
        dw = {}
        d2w = {}
        if order >= 1:
            dtotal = sum(dhat.values())
            dw = {j: (dhat[j] - w_norm[j] * dtotal) / total for j in support}
        if order >= 2:
            d2total = sum(d2hat.values())
            d2w = {j: (d2hat[j] - np.outer(dw[j], dtotal) - np.outer(dtotal, dw[j]) - w_norm[j] * d2total) / total
                   for j in support}
        return w_norm, dw, d2w

    def _weights(self, support, u, du, d2u, order):
        """
        Product-metric scheme: the normalised product-form weights of the
        nodes in support (W, dW of (dim,), d2W of (dim, dim)).
        """

        return self._product_form_weights(support, u, du, d2u, order, self.dim)

    def member_weights(self, i, z, order=0):
        """
        Tensor scheme: the one-dimensional normalised weights of member i at
        the angle z over its own grid. Returns (support: member node indices,
        W dict, dW dict of scalars, d2W dict of scalars, delta (n_i,), image
        steps (n_i,)).
        """

        assert_msg_critical(self.tensor, 'grouped MSI banks: member weights exist for the tensor scheme only.')
        m = int(self.images[i])
        values = self.member_nodes[i]
        r2 = self.member_radius[i] ** 2
        n = len(values)
        delta = np.zeros(n)
        steps = np.zeros(n, dtype=int)
        for a in range(n):
            delta[a], steps[a] = reduce_angle(z - values[a], m)
        u = 2.0 * (1.0 - np.cos(m * delta)) / m ** 2 / r2
        support = [a for a in range(n) if u[a] < 1.0]
        if not support:
            return support, {}, {}, {}, delta, steps
        du = {a: np.array([2.0 * np.sin(m * delta[a]) / m / r2]) for a in support} if order >= 1 else {}
        d2u = {a: np.array([[2.0 * np.cos(m * delta[a]) / r2]]) for a in support} if order >= 2 else {}
        w_norm, dw, d2w = self._product_form_weights(support, u, du, d2u, order, 1)
        return (support, w_norm, {a: float(v[0]) for a, v in dw.items()}, {a: float(v[0, 0]) for a, v in d2w.items()},
                delta, steps)

    def _tensor_weights(self, z, order):
        """
        Tensor scheme: the nodes in support (every member index in its
        one-dimensional support), their weights W_j = prod_i W^(i)_{j_i} with
        the first and second derivatives by the product rule, and the reduced
        displacement and image steps of every node.
        """

        per = [self.member_weights(i, z[i], order) for i in range(self.dim)]
        support = []
        w_norm, dw, d2w = {}, {}, {}
        delta = np.zeros((self.n_nodes, self.dim))
        steps = np.zeros((self.n_nodes, self.dim), dtype=int)
        for j in range(self.n_nodes):
            idx = self.node_index[j]
            if any(idx[i] not in per[i][1] for i in range(self.dim)):
                continue
            support.append(j)
            w1 = [per[i][1][idx[i]] for i in range(self.dim)]
            delta[j] = [per[i][4][idx[i]] for i in range(self.dim)]
            steps[j] = [per[i][5][idx[i]] for i in range(self.dim)]
            w_norm[j] = float(np.prod(w1))
            if order >= 1:
                d1 = [per[i][2][idx[i]] for i in range(self.dim)]
                dw[j] = np.array([d1[i] * float(np.prod([w1[k] for k in range(self.dim) if k != i])) for i in range(self.dim)])
            if order >= 2:
                d2 = [per[i][3][idx[i]] for i in range(self.dim)]
                mat = np.zeros((self.dim, self.dim))
                for i in range(self.dim):
                    for k in range(self.dim):
                        if i == k:
                            mat[i, i] = d2[i] * float(np.prod([w1[q] for q in range(self.dim) if q != i]))
                        else:
                            mat[i, k] = d1[i] * d1[k] * float(np.prod([w1[q] for q in range(self.dim) if q not in (i, k)]))
                d2w[j] = mat
        return support, w_norm, dw, d2w, delta, steps

    def normalised_weights(self, z, order=0):
        """
        The nodes in support at z and their normalised weights W_j (dict),
        with the first derivatives dW_j/dz (dict of (dim,)) for order >= 1.
        """

        z = np.asarray(z, dtype=np.float64).reshape(self.dim)
        if self.tensor:
            support, w_norm, dw, _, _, _ = self._tensor_weights(z, order)
            return support, w_norm, dw
        support, arrays = self._local_models(z, order)
        if not support:
            return support, {}, {}
        u, du, d2u, _ = arrays
        w_norm, dw, _ = self._weights(support, u, du, d2u, order)
        return support, w_norm, dw

    def evaluate_raw(self, z, order=2, names=None):
        """
        The Shepard-Taylor interpolant h(z) of every coefficient (no
        anchoring): dict name -> (value, first, second) with first of shape
        (dim,) + A and second (dim, dim) + A. Critical error outside the
        support (use coverage first).
        """

        z = np.asarray(z, dtype=np.float64).reshape(self.dim)
        names = self.names if names is None else tuple(names)
        if self.tensor:
            support, w_norm, dw, d2w, delta, steps = self._tensor_weights(z, order)
            assert_msg_critical(support, f'grouped MSI banks: the term {self.term} has no node in support at {z.tolist()}.')
            locals_ = {j: self._taylor(j, delta[j], steps[j], order) for j in support}
        else:
            support, arrays = self._local_models(z, order)
            assert_msg_critical(support, f'grouped MSI banks: the term {self.term} has no node in support at {z.tolist()}.')
            u, du, d2u, locals_ = arrays
            w_norm, dw, d2w = self._weights(support, u, du, d2u, order)
        out = {}
        for name in names:
            shape = self._value_shape(name)
            value = np.zeros(shape)
            first = np.zeros((self.dim,) + shape) if order >= 1 else None
            second = np.zeros((self.dim, self.dim) + shape) if order >= 2 else None
            for j in support:
                t_val, t_first, t_second = locals_[j][name]
                value = value + w_norm[j] * t_val
                if order >= 1:
                    first = first + np.multiply.outer(dw[j], t_val) + w_norm[j] * t_first
                if order >= 2:
                    second = second + np.multiply.outer(d2w[j], t_val) + w_norm[j] * t_second
                    cross = np.multiply.outer(dw[j], t_first)  # (dim, dim)+A: dW_i * dt_k
                    second = second + cross + np.swapaxes(cross, 0, 1)
            out[name] = (value, first, second)
        return out

    def evaluate(self, z, order=2, names=None):
        """
        The bank's interpolant: the raw Shepard-Taylor blend for singleton
        banks and for the tensor scheme, the anchored interpolant A_S h for
        the product-metric scheme with |S| >= 2. dict name -> (value, first,
        second).
        """

        z = np.asarray(z, dtype=np.float64).reshape(self.dim)
        if not self.anchored:
            return self.evaluate_raw(z, order, names)
        names = self.names if names is None else tuple(names)
        out = {}
        for name in names:
            def function(point, order_, name=name):
                res = self.evaluate_raw(point, order_, (name,))[name]
                return res
            out[name] = anchoring_operator(function, self.term, z, self.anchor, order)
        return out

    # ---------------------------------------------------------- anchor lines

    def anchor_line_points(self, i):
        """
        Sample points on the anchor line of member i: z_i at the anchor, the
        other members at every node of their grids and at the midpoints
        between neighbouring nodes (the wrap-around midpoint included),
        restricted to the covered points.
        """

        assert_msg_critical(self.tensor, 'grouped MSI banks: anchor lines are sampled for the tensor scheme only.')
        axes = []
        for k in range(self.dim):
            if k == i:
                axes.append(np.array([self.anchor[i]]))
                continue
            m = int(self.images[k])
            period = TWO_PI / m
            values = np.sort(np.mod(self.member_nodes[k] - self.anchor[k], period)) + self.anchor[k]
            mids = 0.5 * (values + np.roll(values, -1))
            mids[-1] = 0.5 * (values[-1] + values[0] + period)
            axes.append(np.concatenate([values, mids]))
        mesh = np.array(np.meshgrid(*axes, indexing='ij')).reshape(self.dim, -1).T
        # a sample point at the very edge of every node's support of a member (u close to 1 for all of them, which a
        # sparse non-uniform grid can produce at a midpoint) carries a total raw weight near zero and normalised weight
        # derivatives of 1e15, so a rounding residual of the Moebius sum becomes an apparent defect there (review of
        # M3, D1); the check samples points where every member has a node well inside its support (u <= 0.9)
        return [z for z in mesh if self.coverage(z)[0] and all(self._member_min_u(k, z[k]) <= 0.9 for k in range(self.dim))]

    def _member_min_u(self, k, z_k):
        """The smallest u = d^2 / r^2 of member k's nodes at the angle z_k (tensor scheme)."""
        m = int(self.images[k])
        r2 = self.member_radius[k] ** 2
        return min(2.0 * (1.0 - np.cos(m * reduce_angle(z_k - v, m)[0])) / m ** 2 / r2 for v in self.member_nodes[k])

    def anchor_line_defect(self):
        """
        Tensor scheme, |S| >= 2: the largest relative defect of the raw
        interpolant on the anchor lines, per coefficient: a dict name ->
        (value, first derivative along the line, second derivative along the
        line), each relative to the largest node datum of that coefficient
        and order (absolute where that scale is zero). Exact data give zero
        to rounding (the nodes on the line carry a vanishing increment with
        vanishing along-line derivatives, and the nodes off the line have
        zero weight and zero weight derivatives there).
        """

        assert_msg_critical(self.tensor and self.dim >= 2,
                            'grouped MSI banks: the anchor-line defect is defined for tensor banks with two or more members.')
        scales = {}
        for name in self.names:
            s_v = s_f = s_s = 0.0
            for j in range(self.n_nodes):
                v, f, s = self._node_arrays(j, name, (0,) * self.dim)
                s_v = max(s_v, float(np.abs(v).max()) if np.size(v) else 0.0)
                if f is not None and self.node_orders[name] >= 1:
                    s_f = max(s_f, float(np.abs(f).max()))
                if s is not None and self.node_orders[name] >= 2:
                    s_s = max(s_s, float(np.abs(s).max()))
            # the cut data behind the increments set the scale where the increments themselves are smaller
            ds = self.data_scale.get(name, 0.0)
            scales[name] = (max(s_v, ds), max(s_f, ds), max(s_s, ds))
        worst = {name: [0.0, 0.0, 0.0] for name in self.names}

        def rel(err, scale):
            # a coefficient whose node data vanish to rounding (every node of an anchor-line bank of the coupling
            # selection lies on the line, M3) has no relative error: compare absolutely below 1e-20
            return err / scale if scale > 1.0e-20 else err

        for i in range(self.dim):
            along = [k for k in range(self.dim) if k != i]
            for z in self.anchor_line_points(i):
                out = self.evaluate_raw(z, 2)
                for name in self.names:
                    value, first, second = out[name]
                    s_v, s_f, s_s = scales[name]
                    worst[name][0] = max(worst[name][0], rel(float(np.abs(value).max()), s_v))
                    worst[name][1] = max(worst[name][1], rel(float(np.abs(first[along]).max()), max(s_v, s_f)))
                    worst[name][2] = max(worst[name][2], rel(float(np.abs(second[np.ix_(along, along)]).max()), max(s_v, s_f, s_s)))
        return {name: tuple(v) for name, v in worst.items()}

    def assert_vanishes_on_anchor_lines(self, tol=1.0e-12):
        """
        Critical error unless the raw interpolant and its derivatives along
        the line vanish on every anchor line to `tol` (relative).
        """

        worst = self.anchor_line_defect()
        bad = {name: v for name, v in worst.items() if max(v) >= tol}
        assert_msg_critical(
            not bad,
            f'grouped MSI banks: the tensor-product bank of the term {self.term} does not vanish on its anchor lines to '
            f'{tol:.0e} (relative defects value / first / second along the line: '
            + ', '.join(f'{name} {v[0]:.1e} / {v[1]:.1e} / {v[2]:.1e}' for name, v in bad.items())
            + '); the node data on the anchor lines are not exact increments of the same faces.')
        return worst

    # -------------------------------------------------------------- transport

    def transported(self, new_term, new_members, transform):
        """
        The same bank for the groups new_term related to self.term by an
        automorphism of the anchor (ADR-004 amendment 7): node points keep
        their anchor-relative phases, and every covariant array is
        transformed by `transform(kind, value, first, second)`.
        """

        new_points = self.node_points - self.anchor + np.array([float(m['anchor']) for m in new_members])
        new_data = []
        for record in self.node_data:
            new_record = {}
            for name, entry in record.items():
                kind = self.kinds[name]
                if kind == 'invariant':
                    new_record[name] = entry
                elif isinstance(entry, dict):
                    new_record[name] = {'faces': {face: (sign, *transform(kind, v, f, s))
                                                  for face, (sign, v, f, s) in entry['faces'].items()}}
                else:
                    new_record[name] = transform(kind, *entry)
            new_data.append(new_record)
        return self._clone(new_term, new_members, new_points, new_data)

    def _clone(self, new_term, new_members, new_points, new_data):
        bank = IncrementBank.__new__(IncrementBank)
        bank.__dict__.update(self.__dict__)
        bank.term = canonical_group_key(new_term)
        bank.images = np.array([int(m['images']) for m in new_members], dtype=int)
        bank.sense = np.array([int(m.get('sense', 1)) for m in new_members], dtype=int)
        bank.anchor = np.array([float(m['anchor']) for m in new_members])
        bank.domain = [m.get('domain') for m in new_members]
        bank.node_points = np.asarray(new_points, dtype=np.float64).reshape(-1, self.dim)
        bank.node_data = new_data
        bank.periods = [TWO_PI / m for m in bank.images]
        if self.tensor:
            # the member grids follow the new anchors; the per-member radii are carried in member order (they
            # depend on the spacings and the images only, so the recomputation must reproduce them)
            bank._prepare_members()
            assert_msg_critical(
                len(bank.member_radius) == len(self.member_radius)
                and max(abs(a - b) for a, b in zip(bank.member_radius, self.member_radius)) < 1.0e-12,
                f'grouped MSI banks: the transported bank of the term {bank.term} does not reproduce the per-member '
                f'support radii of {self.term} ({self.member_radius} against {bank.member_radius}).')
        return bank
