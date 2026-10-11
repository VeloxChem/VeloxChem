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
Anchored Möbius algebra of the grouped MSI v2 model family (SPEC S4, S7;
ADR-004 amendments).

Subsets of a group index set, Möbius signs, the downward closure of a
retained hypergraph (the source-cut dependency graph), Möbius increments and
their reconstruction, the anchoring operator with value, gradient and
Hessian, and the face-closure check of node sets. Everything works for
scalar, vector and matrix valued quantities alike: a cut value is a numpy
array of any shape, and the algebra is linear in it.
"""

from itertools import combinations

import numpy as np

from .errorhandler import assert_msg_critical
from .groupeddefinitions import canonical_group_key


def subsets_of(groups, min_order=0, max_order=None):
    """
    All subsets of the group ids as canonical keys, ordered by order and then
    lexicographically. min_order 0 includes the empty cut ().
    """

    key = canonical_group_key(groups)
    if max_order is None:
        max_order = len(key)
    out = []
    for order in range(min_order, min(max_order, len(key)) + 1):
        out.extend(tuple(c) for c in combinations(key, order))
    return out


def mobius_sign(outer, inner):
    """
    (-1)^(|S| - |T|) for T a subset of S.
    """

    outer = canonical_group_key(outer)
    inner = canonical_group_key(inner)
    assert_msg_critical(set(inner) <= set(outer), f'grouped MSI algebra: {inner} is not a subset of {outer}.')
    return -1 if (len(outer) - len(inner)) % 2 else 1


def downward_closure(retained):
    """
    Every subset of every retained term, the empty cut included: the
    source-cut dependency graph (SPEC S8), as a sorted list of canonical keys.
    """

    closed = set()
    for term in retained:
        closed.update(subsets_of(term))
    return sorted(closed, key=lambda s: (len(s), s))


def source_graph(retained):
    """
    For every retained term, the list of its faces (the cuts its exact
    increment needs), in the order of the Möbius sum.
    """

    return {canonical_group_key(term): subsets_of(term) for term in retained}


def mobius_increment(cut_values, term):
    """
    Delta_S of an array-valued function from its cut values: sum over the
    subsets T of S of (-1)^(|S|-|T|) R_T. cut_values maps canonical keys to
    arrays (or scalars) of one common shape; every face must be present.
    """

    term = canonical_group_key(term)
    total = None
    for face in subsets_of(term):
        assert_msg_critical(face in cut_values, f'grouped MSI algebra: missing source cut {face} for the term {term}.')
        value = np.asarray(cut_values[face], dtype=np.float64)
        contribution = mobius_sign(term, face) * value
        total = contribution if total is None else total + contribution
    return total


def mobius_increments(cut_values, retained):
    """
    The increments of every term in the downward closure of the retained
    hypergraph.
    """

    return {term: mobius_increment(cut_values, term) for term in downward_closure(retained)}


def reconstruct(increments, term):
    """
    The cut value R_S from the increments of all subsets of S.
    """

    term = canonical_group_key(term)
    total = None
    for face in subsets_of(term):
        assert_msg_critical(face in increments, f'grouped MSI algebra: missing increment {face} for the cut {term}.')
        value = np.asarray(increments[face], dtype=np.float64)
        total = value.copy() if total is None else total + value
    return total


def restricted_point(point, term, members, anchor):
    """
    The point with the coordinates of the members outside the subset
    `members` of `term` replaced by their anchor values. point and anchor are
    arrays over the coordinates of term (in canonical order).
    """

    term = canonical_group_key(term)
    members = canonical_group_key(members)
    out = np.array(anchor, dtype=np.float64).copy()
    for i, g in enumerate(term):
        if g in members:
            out[i] = point[i]
    return out


def anchoring_operator(function, term, point, anchor, order=2):
    """
    The anchored interpolant A_S h at a point of the slice of S (SPEC S7):

        A_S h(z) = sum over T subset of S of (-1)^(|S|-|T|) h(z_T, z^a_(S minus T))

    function(z, order) returns (value, gradient, hessian) of h at z with
    gradient of shape (|S|,) + value.shape and hessian (|S|, |S|) + value.shape
    (None for the derivatives above `order`). Returns the same triple for the
    anchored function. It vanishes identically wherever a coordinate of S is
    at its anchor value, and derivatives with respect to a coordinate that a
    restricted evaluation holds fixed are zero.
    """

    term = canonical_group_key(term)
    n = len(term)
    value = None
    gradient = None
    hessian = None
    for face in subsets_of(term):
        sign = mobius_sign(term, face)
        z = restricted_point(point, term, face, anchor)
        val, grad, hess = function(z, order)
        val = np.asarray(val, dtype=np.float64)
        value = sign * val if value is None else value + sign * val
        if order >= 1:
            grad = np.asarray(grad, dtype=np.float64)
            if gradient is None:
                gradient = np.zeros((n,) + val.shape)
            for i, g in enumerate(term):
                if g in face:
                    gradient[i] += sign * grad[i]
        if order >= 2:
            hess = np.asarray(hess, dtype=np.float64)
            if hessian is None:
                hessian = np.zeros((n, n) + val.shape)
            for i, g in enumerate(term):
                if g not in face:
                    continue
                for j, h in enumerate(term):
                    if h in face:
                        hessian[i, j] += sign * hess[i, j]
    return value, gradient, hessian


def periodic_difference(a, b, period):
    """
    a - b reduced to (-period/2, period/2]; period None or 0 means no
    reduction.
    """

    d = np.asarray(a, dtype=np.float64) - np.asarray(b, dtype=np.float64)
    if period is None or period == 0.0:
        return d
    return (d + 0.5 * period) % period - 0.5 * period


def missing_face_nodes(nodes, term, anchor, periods, tol=1.0e-9):
    """
    Face closure of a node set (ADR-004 amendment 5): for every node and
    every proper subset T of S, the projection of the node onto the face
    (anchor values outside T) must itself be a node. Returns the missing
    points as a list of (node index, face, projected point); empty when the
    set is face-closed. nodes is an array (n_nodes, |S|); periods gives the
    period per coordinate (None for a non-periodic one).
    """

    term = canonical_group_key(term)
    nodes = np.asarray(nodes, dtype=np.float64).reshape(-1, len(term))
    anchor = np.asarray(anchor, dtype=np.float64)
    missing = []
    for j, node in enumerate(nodes):
        for face in subsets_of(term, max_order=len(term) - 1):
            projected = restricted_point(node, term, face, anchor)
            found = False
            for other in nodes:
                dist = max(abs(periodic_difference(projected[i], other[i], periods[i])) for i in range(len(term)))
                if dist <= tol:
                    found = True
                    break
            if not found:
                missing.append((j, face, projected))
    return missing


def assert_face_closed(nodes, term, anchor, periods, tol=1.0e-9):
    """
    Critical error unless the node set is face-closed.
    """

    missing = missing_face_nodes(nodes, term, anchor, periods, tol)
    assert_msg_critical(
        not missing,
        f'grouped MSI algebra: the node set of the term {canonical_group_key(term)} is not face-closed; '
        f'first missing face node: node {missing[0][0]} projected onto {missing[0][1]} at {missing[0][2]}.' if missing else '')
