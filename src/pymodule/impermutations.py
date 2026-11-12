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
Feasible permutations of identical atoms for the permutation-symmetric
Shepard interpolation (InterpolationDriver, permutation_symmetry).

The potential energy surface does not change when identical atoms swap
places. A datapoint computed at the geometry X therefore also gives the energy,
gradient and Hessian at every relabelled geometry P X. The interpolation can
use them without new data by relabelling the current geometry instead:

    E(x) = sum_P sum_i w_i(P x) T_i(P x) / sum_P sum_i w_i(P x).

Only permutations that the dynamics can reach without breaking a bond are
used ("feasible" permutations):

    - 'methyl': the three cyclic relabellings of the three terminal atoms
      of a CX3 group (CH3, CF3, NH3+, ...), i.e. its rotations by 120 deg;
    - 'ring': the 180 deg flip of a six-membered ring that hangs on a single
      rotatable bond (phenyl, p-fluorophenyl, 4-pyridyl, ...): the ortho and
      the meta atoms swap, with their terminal substituents.

The permutations of all groups are combined (every product), so the set is
closed and the interpolant is exactly invariant under it.

A permutation is an index array perm: atom a of the relabelled geometry is
atom perm[a] of the original one, i.e. the relabelled coordinates are
x[perm].
"""

import itertools

import numpy as np

from .errorhandler import assert_msg_critical
from .imrowdescription import RowDescription

PERMUTATION_KINDS = ('none', 'methyl', 'ring', 'all')
# the interpolation costs one evaluation of every datapoint per permutation
MAX_PERMUTATIONS = 1024


def _neighbours(bonds, n_atoms):

    nb = {a: set() for a in range(n_atoms)}
    for i, j in bonds:
        nb[i].add(j)
        nb[j].add(i)
    return nb


def methyl_rotations(bonds, labels):
    """
    The rotations of every CX3 group: an atom with exactly four neighbours,
    three of them terminal atoms of one element.

    :param bonds:
        The bonds (0-based atom pairs) of the z-matrix.
    :param labels:
        The element labels.

    :return:
        One list per group with its three permutations (identity first).
    """

    n = len(labels)
    nb = _neighbours(bonds, n)
    groups = []
    for a in range(n):
        if len(nb[a]) != 4:
            continue
        terminal = sorted(b for b in nb[a] if len(nb[b]) == 1)
        for element in sorted(set(labels[t] for t in terminal)):
            same = [t for t in terminal if labels[t] == element]
            if len(same) != 3:
                continue
            ops = []
            for shift in range(3):
                perm = np.arange(n)
                perm[same] = np.roll(same, -shift)
                ops.append(perm)
            groups.append(ops)
    return groups


def ring_flips(bonds, labels):
    """
    The 180 deg flips of every six-membered ring with exactly one substituent
    that is not a single atom (the ipso substituent, through which the ring
    rotates) and mirror-identical ortho and meta positions.

    :param bonds:
        The bonds (0-based atom pairs) of the z-matrix.
    :param labels:
        The element labels.

    :return:
        One list per ring, [identity, flip].
    """

    n = len(labels)
    bonds = [tuple(int(x) for x in b) for b in bonds]
    nb = _neighbours(bonds, n)
    ring_bonds = [b for b, r in zip(bonds, RowDescription._ring_bond_flags(bonds)) if r]
    rnb = _neighbours(ring_bonds, n)

    # the rings: connected sets of ring atoms where every atom has two ring
    # neighbours (fused rings are skipped)
    rings, seen = [], set()
    for a in range(n):
        if not rnb[a] or a in seen:
            continue
        ring, front = {a}, {a}
        while front:
            front = {y for x in front for y in rnb[x]} - ring
            ring |= front
        seen |= ring
        if len(ring) == 6 and all(len(rnb[x]) == 2 for x in ring):
            rings.append(ring)

    def substituents(x, ring):
        return [y for y in nb[x] if y not in ring]

    groups = []
    for ring in rings:
        branched = [x for x in ring if any(len(nb[y]) > 1 for y in substituents(x, ring))]
        if len(branched) != 1:
            continue
        ipso = branched[0]
        if len(substituents(ipso, ring)) != 1:
            continue
        o1, o2 = sorted(rnb[ipso])
        m1 = next(iter(rnb[o1] - {ipso}))
        m2 = next(iter(rnb[o2] - {ipso}))
        pairs = [(o1, o2), (m1, m2)]
        mirror = True
        for x, y in ((o1, o2), (m1, m2)):
            sx, sy = substituents(x, ring), substituents(y, ring)
            mirror &= labels[x] == labels[y] and len(sx) == len(sy) <= 1
            if mirror and sx:
                mirror &= labels[sx[0]] == labels[sy[0]]
                pairs.append((sx[0], sy[0]))
        if not mirror:
            continue
        flip = np.arange(n)
        for x, y in pairs:
            flip[x], flip[y] = y, x
        groups.append([np.arange(n), flip])
    return groups


def feasible_permutations(bonds, labels, kind):
    """
    All products of the permutations of kind ('none', 'methyl', 'ring' or
    'all'), the identity first.

    :param bonds:
        The bonds (0-based atom pairs) of the z-matrix.
    :param labels:
        The element labels.
    :param kind:
        One of PERMUTATION_KINDS.

    :return:
        The list of permutations (index arrays).
    """

    assert_msg_critical(kind in PERMUTATION_KINDS,
                        f'feasible_permutations: kind must be one of {PERMUTATION_KINDS}, not {kind}.')
    groups = []
    if kind in ('methyl', 'all'):
        groups += methyl_rotations(bonds, labels)
    if kind in ('ring', 'all'):
        groups += ring_flips(bonds, labels)

    size = int(np.prod([len(g) for g in groups])) if groups else 1
    assert_msg_critical(
        size <= MAX_PERMUTATIONS,
        f'feasible_permutations: {size} permutations ({len(groups)} groups) '
        f'exceed {MAX_PERMUTATIONS}; every datapoint is evaluated once per '
        "permutation. Use 'methyl' or 'ring' instead of 'all'.")

    perms = []
    for combo in itertools.product(*groups):
        perm = np.arange(len(labels))
        for op in combo:
            perm = perm[op]
        perms.append(perm)
    return perms if perms else [np.arange(len(labels))]
