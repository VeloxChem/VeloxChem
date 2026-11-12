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

import numpy as np

from .errorhandler import assert_msg_critical

# weight factor of the exocyclic proper dihedrals of saturated rings (lever)
LEVER_EXO_SAT = 0.25
# weight of ring bonds under the bond rule (independent of r0)
RING_BOND_WEIGHT = 2.0
ROTOR_FORMS = ('none', 'b1', 'c5')


def _wrap(x):
    """Wraps angle differences to [-pi, pi)."""

    return (x + np.pi) % (2.0 * np.pi) - np.pi


class RowDescription:
    """
    Describes the rows of the internal coordinates for the Taylor models of
    the interpolation, T = E0 + g.d + 1/2 d.H.d.

    Two choices define a description: the displacement d of each row (its
    shape) and the weight of each row in the pseudo-inverse that turns the
    Cartesian gradient and Hessian into g and H. The internal coordinates are
    redundant, so the weights decide how a motion is shared between rows. Every
    shape has unit slope at the datapoint value p0; then the B rows are the
    plain B rows (bonds in r) and the weight alone sets the share of a row.
    VeloxChem's 1/r bond is the shape r0 (1 - r0/r) with the weight 1/r0^2
    (bohr).

    With every switch at its default this is VeloxChem's own description:
    bonds 1/r, angles theta, proper dihedrals 2 sin(Delta/2), impropers
    2 tan(Delta/2). The switches are independent:

        - bond_scale: the weight of every bond row times c. Alone this is the
          inverse bond coordinate scaled by c, q = c/r.
        - bond_rule: bonds between heavy atoms get the shape
          (r0/2)(1 - r0^2/r^2); ring bonds (a bond whose removal keeps both
          atoms connected) get the weight 2. X-H bonds keep the 1/r shape.
        - angle_rule: angles with an H end that are not ring angles get the
          shape (cos t0 - cos t) / max(sin t0, 0.2).
        - exo_sat_lever: the weight of the proper dihedrals around a ring bond
          that are not ring torsions, in saturated rings (the two central
          atoms not both sp2-like), times 1/4.

    The rotor form (an argument of potential) is chosen at evaluation:
        - 'none': every proper dihedral row is 2 sin(Delta/2).
        - 'b1': the dihedral rows of every non-ring central bond carry the
          collective rotation tau (their mean displacement) as the m-fold term
          (2/m) sin(m tau/2) instead of its one-fold part. m is 3 for a
          methyl-like end or two non-sp2 ends, 2 for two sp2-like ends, and 1
          (no change) otherwise.
        - 'c5': b1, continuous through tau = 180 deg: the rows of a group are
          unwrapped onto the branch of their circular mean, and the row
          deviations (and, for odd m, the terms linear in the collective
          term) are tapered as cos^2 from the window start to 0 at 180 deg.

    The classes come from the z-matrix bonds and the element labels only.

    Instance variables
        - key: the switches that change g and H; datapoints keep the key of
          the description their g and H were made with.
        - cls: the class of every row.
        - groups: (rows, m) of every non-ring central bond.
        - n_bonds, n_angles, n_rows.
    """

    def __init__(self, z_matrix, labels, bond_scale=1.0, bond_rule=False,
                 angle_rule=False, exo_sat_lever=False):
        """
        Initializes the description of the rows of z_matrix.

        :param z_matrix:
            The z-matrix dictionary {'bonds', 'angles', 'dihedrals',
            'impropers'} of 0-based atom-index tuples.
        :param labels:
            The element labels of the atoms.
        :param bond_scale:
            The weight factor c of every bond row (> 0).
        :param bond_rule:
            The bond rule (see the class docstring).
        :param angle_rule:
            The angle rule.
        :param exo_sat_lever:
            The lever on the exocyclic dihedrals of saturated rings.
        """

        assert_msg_critical(float(bond_scale) > 0.0,
                            'RowDescription: bond_scale must be positive.')

        labels = [str(x) for x in labels]
        bonds = [tuple(int(x) for x in b) for b in z_matrix['bonds']]
        angles = [tuple(int(x) for x in a) for a in z_matrix['angles']]
        dihedrals = [tuple(int(x) for x in d) for d in z_matrix['dihedrals']]
        impropers = [tuple(int(x) for x in d) for d in z_matrix['impropers']]

        self.key = (float(bond_scale), bool(bond_rule), bool(angle_rule),
                    bool(exo_sat_lever))

        types = self._atom_types(bonds, labels)
        ring = {frozenset(b) for b, r in zip(bonds, self._ring_bond_flags(bonds)) if r}
        is_h = [x == 'H' for x in labels]

        self.n_bonds = len(bonds)
        self.n_angles = len(angles)
        n_proper = len(dihedrals)
        dihedral_start = self.n_bonds + self.n_angles

        cls = []
        for i, j in bonds:
            if is_h[i] or is_h[j]:
                cls.append('bond XH')
            elif frozenset((i, j)) in ring:
                cls.append('bond ring')
            else:
                cls.append('bond chain')
        for i, j, k in angles:
            if frozenset((i, j)) in ring and frozenset((j, k)) in ring:
                cls.append('angle ring')
            elif is_h[i] or is_h[k]:
                cls.append('angle H')
            else:
                cls.append('angle chain')
        exo_sat = []
        groups = {}
        for r, (i, j, k, l) in enumerate(dihedrals):
            central = frozenset((j, k)) in ring
            ends = frozenset((i, j)) in ring and frozenset((k, l)) in ring
            if central and ends:
                cls.append('proper ring')
            elif central:
                cls.append('proper exo')
                if not (types['sp2'][j] and types['sp2'][k]):
                    exo_sat.append(dihedral_start + r)
            else:
                cls.append('proper rotor')
                groups.setdefault(tuple(sorted((j, k))), []).append(dihedral_start + r)
        cls += ['improper'] * len(impropers)

        self.cls = np.array(cls, dtype=object)
        self.n_rows = len(cls)
        bond_rows = np.arange(self.n_bonds)
        ring_bond = np.array([frozenset(b) in ring for b in bonds], dtype=bool)

        # shapes: power (exponent n), linear, cos, half-sine, half-tan
        self.n = np.full(self.n_rows, np.nan)
        self.n[bond_rows] = -1.0
        heavy_bond = (self.cls[:self.n_bonds] != 'bond XH')
        if bond_rule:
            self.n[bond_rows[heavy_bond]] = -2.0
        self._power = np.zeros(self.n_rows, dtype=bool)
        self._power[bond_rows] = True
        self._cos = (self.cls == 'angle H') if angle_rule else np.zeros(self.n_rows, dtype=bool)
        self._half_sine = np.zeros(self.n_rows, dtype=bool)
        self._half_sine[dihedral_start:dihedral_start + n_proper] = True
        self._half_tan = np.zeros(self.n_rows, dtype=bool)
        self._half_tan[dihedral_start + n_proper:] = True

        # weights: const_weight / p0^2 for the rows with native_weight
        # (bonds, 1/r0^2 in bohr), const_weight otherwise
        self.const_weight = np.ones(self.n_rows)
        self.native_weight = np.zeros(self.n_rows, dtype=bool)
        self.native_weight[bond_rows] = True
        if bond_rule:
            self.native_weight[bond_rows[ring_bond]] = False
            self.const_weight[bond_rows[ring_bond]] = RING_BOND_WEIGHT
        self.const_weight[bond_rows] *= float(bond_scale)
        if exo_sat_lever and exo_sat:
            self.const_weight[np.array(exo_sat)] *= LEVER_EXO_SAT

        # rotor groups: every non-ring central bond; m from the end groups
        self.groups = []
        for bond, rows in groups.items():
            m = self._periodicity(bond, types)
            self.groups.append((np.array(rows), 1 if m is None else m))
        n_groups = len(self.groups)
        width = max((len(rows) for rows, _ in self.groups), default=0)
        self._g_idx = np.zeros((n_groups, width), dtype=int)
        self._g_mask = np.zeros((n_groups, width), dtype=bool)
        self._g_sum = np.zeros((n_groups, self.n_rows))
        for g, (rows, _) in enumerate(self.groups):
            self._g_idx[g, :len(rows)] = rows
            self._g_mask[g, :len(rows)] = True
            self._g_sum[g, rows] = 1.0
        self._g_size = self._g_mask.sum(axis=1)
        self._g_m = np.array([m for _, m in self.groups], dtype=float)
        self._g_odd = (self._g_m % 2 == 1)

    @staticmethod
    def _ring_bond_flags(bonds):
        """
        True for the bonds whose removal keeps both atoms connected (through
        the z-matrix bonds).
        """

        flags = []
        for cut in bonds:
            reached, front = {cut[0]}, {cut[0]}
            while front:
                front = {y for x in front for b in bonds if b != cut and x in b
                         for y in b} - reached
                reached |= front
            flags.append(cut[1] in reached)
        return np.array(flags, dtype=bool)

    @staticmethod
    def _atom_types(bonds, labels):
        """
        Per atom: sp2-like (C with <= 3 neighbours, N with <= 2, terminal O
        or S, and N, O, S bonded to such an atom), linear (C with 2
        neighbours, N with 1) and methyl-like (three or more terminal atoms
        of one element).
        """

        neighbours = {a: set() for a in range(len(labels))}
        for i, j in bonds:
            neighbours[i].add(j)
            neighbours[j].add(i)
        degree = {a: len(v) for a, v in neighbours.items()}
        unsat = {a: (labels[a] == 'C' and degree[a] <= 3) or
                    (labels[a] == 'N' and degree[a] <= 2) or
                    (labels[a] in ('O', 'S') and degree[a] == 1)
                 for a in neighbours}
        sp2 = {a: unsat[a] or (labels[a] in ('N', 'O', 'S') and
                               any(unsat[b] for b in neighbours[a]))
               for a in neighbours}
        linear = {a: (labels[a] == 'C' and degree[a] == 2) or
                     (labels[a] == 'N' and degree[a] == 1)
                  for a in neighbours}
        methyl = set()
        for a in neighbours:
            terminal = [b for b in neighbours[a] if degree[b] == 1]
            elements = set(labels[t] for t in terminal)
            if any(sum(labels[t] == e for t in terminal) >= 3 for e in elements):
                methyl.add(a)
        return {'sp2': sp2, 'linear': linear, 'methyl': methyl}

    @staticmethod
    def _periodicity(bond, types):
        """The rotor periodicity m of a central bond; None if ambiguous."""

        a, b = bond
        if types['linear'][a] or types['linear'][b]:
            return None
        if a in types['methyl'] or b in types['methyl']:
            return 3
        if not types['sp2'][a] and not types['sp2'][b]:
            return 3
        if types['sp2'][a] and types['sp2'][b]:
            return 2
        return None

    def weights(self, p0):
        """The pseudo-inverse weights of the rows at the plain values p0."""

        w = self.const_weight.copy()
        w[self.native_weight] = (self.const_weight[self.native_weight] /
                                 p0[self.native_weight]**2)
        return w

    def curvature(self, p0):
        """The second derivatives f''(p0) of the unit-slope shapes."""

        k = np.zeros(self.n_rows)
        k[self._power] = (self.n[self._power] - 1.0) / p0[self._power]
        k[self._cos] = (np.cos(p0[self._cos]) /
                        np.maximum(np.sin(p0[self._cos]), 0.2))
        return k

    def displacement(self, p, p0):
        """
        The displacements d of the rows and dd/dp at the plain values p
        (bonds in bohr, angles and dihedrals in radian).
        """

        d = p - p0
        chain = np.ones(self.n_rows)

        n, r0 = self.n[self._power], p0[self._power]
        x = p[self._power] / r0
        d[self._power] = r0 * (x**n - 1.0) / n
        chain[self._power] = x**(n - 1.0)

        s0 = np.maximum(np.sin(p0[self._cos]), 0.2)
        d[self._cos] = (np.cos(p0[self._cos]) - np.cos(p[self._cos])) / s0
        chain[self._cos] = np.sin(p[self._cos]) / s0

        x = _wrap(p[self._half_sine] - p0[self._half_sine])
        d[self._half_sine] = 2.0 * np.sin(0.5 * x)
        chain[self._half_sine] = np.cos(0.5 * x)

        x = _wrap(p[self._half_tan] - p0[self._half_tan])
        d[self._half_tan] = 2.0 * np.tan(0.5 * x)
        chain[self._half_tan] = 1.0 / np.maximum(np.cos(0.5 * x)**2, 1.0e-12)

        return d, chain

    def internal_derivatives(self, p0, b_matrix, b2_matrix, gradient, hessian):
        """
        Returns g and H of a datapoint in the rows of this description,
        g = L^T g_x and H = L^T [H_x - sum_q g_q B2_q - sum_q g_q f''_q
        b_q b_q^T] L, with L = (S B)^+ S truncated at the physical rank.

        :param p0:
            The plain internal coordinates of the datapoint.
        :param b_matrix:
            The B-matrix with bonds in r, in the metric of gradient and hessian.
        :param b2_matrix:
            The second derivatives of the coordinates, in the same metric.
        :param gradient:
            The Cartesian gradient (mass-weighted if the database is).
        :param hessian:
            The Cartesian Hessian (mass-weighted if the database is).
        """

        n_atoms = b_matrix.shape[1] // 3
        weight = self.weights(p0)

        # rank of the unweighted rows, capped at 3N - 6
        sv = np.linalg.svd(b_matrix, compute_uv=False)
        rank = min(max(3 * n_atoms - 6, 1), int(np.sum(sv > 1.0e-8 * sv[0])))

        u, s, vt = np.linalg.svd(weight[:, None] * b_matrix, full_matrices=False)
        resolved = (s[0] / s[rank - 1] <= 1.0e10 and
                    (rank == len(s) or s[rank] / s[rank - 1] <= 1.0e-6))
        assert_msg_critical(
            resolved,
            f'RowDescription: the rank {rank} of the weighted B-matrix is not '
            f'resolved (condition {s[0] / s[rank - 1]:.1e}); the weights '
            f'(bond_scale = {self.key[0]}) make the transformation ill-defined.')
        inverse = ((vt[:rank].T / s[:rank]) @ u[:, :rank].T) * weight[None, :]

        g = inverse.T @ np.asarray(gradient, dtype=float).reshape(-1)
        k = self.curvature(p0)
        residual = (0.5 * (hessian + hessian.T) -
                    np.einsum('q,qij->ij', g, b2_matrix, optimize=True))
        bent = np.flatnonzero(k != 0.0)
        residual -= np.einsum('q,qi,qj->ij', g[bent] * k[bent], b_matrix[bent],
                              b_matrix[bent], optimize=True)
        h = inverse.T @ residual @ inverse
        return g, 0.5 * (h + h.T)

    def group_sums(self, g, hessian):
        """
        The sums of g over each rotor group and of H over each block of two
        groups (fixed per datapoint).
        """

        return self._g_sum @ g, self._g_sum @ hessian @ self._g_sum.T

    @staticmethod
    def _window(tau, start):
        """
        w(tau) = 1 for |tau| <= start, cos^2 taper to 0 at 180 deg, and
        dw/dtau. start=None gives w = 1.
        """

        if start is None:
            return np.ones_like(tau), np.zeros_like(tau)
        a = np.abs(tau)
        arg = 0.5 * np.pi * np.clip((a - start) / (np.pi - start), 0.0, 1.0)
        c, s = np.cos(arg), np.sin(arg)
        inside = (a <= start)
        dw = -2.0 * c * s * 0.5 * np.pi / (np.pi - start) * np.sign(tau)
        return np.where(inside, 1.0, c * c), np.where(inside, 0.0, dw)

    def potential(self, energy, g, hessian, sums, p, p0, b_matrix,
                  rotor_form='none', window_deg=150.0):
        """
        Returns the Taylor model of one datapoint: the energy, its gradient
        B^T (...) with respect to the coordinates of b_matrix (with the chain
        rule through the collective rotations), and the internal vector
        dT/dd * dd/dp.

        :param energy, g, hessian:
            E0, g and H of the datapoint in the rows of this description.
        :param sums:
            group_sums(g, hessian).
        :param p, p0:
            The plain internal coordinates at the current geometry and at the
            datapoint.
        :param b_matrix:
            The B-matrix (bonds in r) at the current geometry.
        :param rotor_form:
            'none', 'b1' or 'c5'.
        :param window_deg:
            The start of the C5 window, degree.
        """

        d, chain = self.displacement(p, p0)

        if rotor_form == 'none' or not self.groups:
            hd = hessian @ d
            internal = (g + hd) * chain
            return energy + d @ g + 0.5 * d @ hd, b_matrix.T @ internal, internal

        idx, mask, size, m = self._g_idx, self._g_mask, self._g_size, self._g_m
        p_sum = self._g_sum
        g_sum, k_sum = sums

        dphi = _wrap(p[idx] - p0[idx])
        if rotor_form == 'c5':
            # unwrap onto the branch of the circular mean of the group
            tau_c = np.arctan2(np.where(mask, np.sin(dphi), 0.0).sum(1) / size,
                               np.where(mask, np.cos(dphi), 0.0).sum(1) / size)
            dphi = dphi + 2.0 * np.pi * np.round((tau_c[:, None] - dphi) / (2.0 * np.pi))
            w, dw = self._window(_wrap(np.where(mask, dphi, 0.0).sum(1) / size),
                                 np.radians(window_deg))
        tau = np.where(mask, dphi, 0.0).sum(1) / size
        if rotor_form == 'b1':
            w, dw = self._window(tau, None)

        omega = np.where(self._g_odd, w, 1.0)
        domega = np.where(self._g_odd, dw, 0.0)
        f1, f1_tau = 2.0 * np.sin(0.5 * dphi), 2.0 * np.sin(0.5 * tau)
        f_m, df_m = 2.0 / m * np.sin(0.5 * m * tau), np.cos(0.5 * m * tau)

        delta, ch = d.copy(), chain.copy()
        delta[idx[mask]] = (w[:, None] * (f1 - f1_tau[:, None]))[mask]
        ch[idx[mask]] = (w[:, None] * np.cos(0.5 * dphi))[mask]

        h_delta = hessian @ delta
        h_delta_sum = p_sum @ h_delta
        k_diag = np.diag(k_sum)
        of = omega * f_m
        off = k_sum @ of - k_diag * of

        e = (energy + delta @ g + 0.5 * delta @ h_delta + of @ (g_sum + h_delta_sum) +
             0.5 * np.sum(f_m * f_m * k_diag) + 0.5 * of @ off)

        de_ddelta = g + h_delta + hessian @ (p_sum.T @ of)
        de_df = omega * (g_sum + h_delta_sum) + f_m * k_diag + omega * off
        de_domega = f_m * (g_sum + h_delta_sum) + f_m * off
        internal = de_ddelta * ch
        s = p_sum @ de_ddelta
        t = np.where(mask, de_ddelta[idx] * (f1 - f1_tau[:, None]), 0.0).sum(1)
        coef = -w * np.cos(0.5 * tau) * s + dw * t + de_df * df_m + de_domega * domega

        return e, b_matrix.T @ internal + (coef / size) @ (p_sum @ b_matrix), internal

    @staticmethod
    def plain_bonds(values, b_matrix, n_bonds, use_inverse_bond_length):
        """
        The internal coordinates and B-matrix with the bonds as r, from the
        inverse (1/r) or plain values of VeloxChem.
        """

        p = np.array(values, dtype=float)
        b = np.array(b_matrix, dtype=float)
        if use_inverse_bond_length:
            r = 1.0 / p[:n_bonds]
            p[:n_bonds] = r
            b[:n_bonds] *= -(r * r)[:, None]
        return p, b
