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
The reference model Q_a of the grouped MSI v2 model family: policy
`image_atlas_periodic_rotor_v1` (SPEC S2, S3; ADR-002 with its amendments).

Q_a(x) = sum_b lambda_b(x) T_b(x) over the image charts b of the product of
the cyclic groups of all groups with images. Every chart is the legacy local
Taylor model of the permuted anchor data (InterpolationDatapoint transforms
the permuted Cartesian data to internal coordinates), with the periodic
collective displacement on the rotor rows of every group:

    s_k = S_n(Theta_g) + 2 sin(e_k / 2),  S_n = sin(n Theta) / n,
    T_b = U_a + g_q . s + 1/2 s^T H_q s + 1/2 sum_g [t_n(Theta_g) - S_n(Theta_g)^2] sigma_g^T H_q sigma_g,
    omega_b = prod_g (1 - Theta_g^2 / delta^2)^4 for |Theta_g| < delta,  lambda_b = omega_b / sum_c omega_c.

Energy, Cartesian gradient and Cartesian Hessian are analytic:
grad Q = B^T dQ/dq and grad2 Q = B^T (d2Q/dq2) B + sum_k (dQ/dq_k) B2_k, with
the first and second derivatives of s, of the phases and of the atlas
weights. Without groups the chart is InterpolationDriver.compute_potential
(pinned by a test).
"""

import numpy as np

from .errorhandler import assert_msg_critical
from .groupedcoordinates import rigid_body_projector, wrap_angle


def _taper(theta, delta, order=2):
    """
    omega(Theta) = (1 - (Theta/delta)^2)^4 inside |Theta| < delta, zero
    outside, with its first and second derivatives.
    """

    u = theta / delta
    if abs(u) >= 1.0:
        return 0.0, 0.0, 0.0
    one = 1.0 - u * u
    value = one ** 4
    first = -8.0 * u * one ** 3 / delta
    second = (-8.0 * one ** 3 + 48.0 * u * u * one ** 2) / (delta * delta)
    return value, first, second


class GroupedReference:
    """
    Q_a of one core chart.

    :param model: the GroupedCoordinateModel of the anchor
    :param energy: the anchor energy U_a (hartree)
    :param gradient: the anchor Cartesian gradient (N, 3); its net force and
        torque are removed with the rigid-body projector of the anchor
    :param hessian: the anchor Cartesian Hessian (3N, 3N); its symmetric
        part is used
    :param delta_degrees: the half-width of the atlas plateau transition;
        None takes the policy value
    """

    def __init__(self, model, energy, gradient, hessian, delta_degrees=None):
        self.model = model
        self.ic = model.ic
        self.n_atoms = model.n_atoms
        self.delta = np.radians(model.policy.atlas_delta_degrees if delta_degrees is None else delta_degrees)
        assert_msg_critical(np.pi / 3.0 < self.delta < 2.0 * np.pi / 3.0,
                            'grouped MSI reference: delta must lie strictly between 60 and 120 degrees.')
        self.x_anchor = np.array(model.x_anchor, dtype=np.float64)
        gradient = np.asarray(gradient, dtype=np.float64).reshape(-1)
        projector = rigid_body_projector(self.x_anchor)
        projected = projector @ gradient
        self.rigid_body_gradient_norm = float(np.linalg.norm(gradient - projected))
        hessian = np.asarray(hessian, dtype=np.float64)
        self.hessian_asymmetry = float(np.abs(hessian - hessian.T).max())
        hessian = 0.5 * (hessian + hessian.T)
        self.energy_anchor = float(energy)
        self.gradient_anchor = projected.reshape(self.n_atoms, 3)
        self.hessian_anchor = hessian
        self.atlas_groups = [gid for gid, g in enumerate(model.groups) if g.images > 1]
        self.charts = []
        self.legacy_margins = []
        margin = model.policy.legacy_cutoff_margin
        cutoff = 1.0e-6   # the absolute cut-off of InterpolationDatapoint.transform_gradient_and_hessian
        for expo, perm in model.image_elements():
            p3 = (3 * perm[:, None] + np.arange(3)[None, :]).reshape(-1)
            x_chart = self.x_anchor[perm]
            # the legacy cut-off must be unambiguous at every chart anchor (ADR-002 amendment, M1-C decision 3)
            _, b_chart = self.ic.values_and_b(x_chart)
            sv = np.linalg.svd(b_chart @ b_chart.T, compute_uv=False)
            r = 3 * self.n_atoms - 6
            sigma_r1 = float(sv[r]) if len(sv) > r else 0.0
            m_chart = min(sv[r - 1] / cutoff, cutoff / sigma_r1 if sigma_r1 > 0.0 else np.inf)
            self.legacy_margins.append({'chart': expo, 'sigma_r': float(sv[r - 1]), 'sigma_r1': sigma_r1, 'margin': float(m_chart)})
            assert_msg_critical(
                m_chart >= margin,
                f'grouped MSI reference: the legacy singular-value cut-off (1e-6) is not unambiguous at the chart anchor {expo}: '
                f'sigma_r = {sv[r - 1]:.2e}, sigma_r+1 = {sigma_r1:.2e}, margin {m_chart:.1e} below legacy_cutoff_margin = {margin:g}. '
                'The anchor geometry is close to a Z-matrix singularity (a near-linear angle or a degenerate improper): '
                'choose another anchor or Z-matrix.')
            chart = self.ic.chart_data(x_chart, self.energy_anchor, self.gradient_anchor[perm],
                                       hessian[np.ix_(p3, p3)])
            chart['exponents'] = expo
            chart['permutation'] = perm
            chart['phi0'] = model.phases(chart['q'], order=0)[0]
            chart['curvature'] = np.array([chart['H_q'][np.ix_(list(g.rotor_rows), list(g.rotor_rows))].sum()
                                           for g in model.groups])
            self.charts.append(chart)

    # ------------------------------------------------------------- one chart

    def chart_value(self, chart, q, phases, order=2):
        """
        T_b and its derivatives with respect to q at the internal
        coordinates q. phases = (phi, dphi/dq, d2phi/dq2) from
        model.phases(q). Returns (T, dT/dq, d2T/dq2, diagnostics).
        """

        ic = self.ic
        model = self.model
        phi, dphi, d2phi = phases
        nq = ic.n_rows
        d = q - chart['q']
        s, s1, s2 = ic.displacement(d)
        rotor = np.zeros(nq, dtype=bool)
        jac = None
        if order >= 1:
            jac = np.diag(s1.copy())
        diag_second = s2.copy() if order >= 2 else None
        blocks = []  # per group: (rows, data needed for the second derivative)
        extra = 0.0
        dextra = np.zeros(nq) if order >= 1 else None
        d2extra = np.zeros((nq, nq)) if order >= 2 else None
        thetas = np.zeros(model.n_groups)
        max_dev = np.zeros(model.n_groups)
        for gid, g in enumerate(model.groups):
            rows = list(g.rotor_rows)
            rotor[rows] = True
            n = g.harmonic
            theta = wrap_angle(phi[gid] - chart['phi0'][gid])
            thetas[gid] = theta
            e = wrap_angle(d[rows] - theta)
            max_dev[gid] = np.abs(e).max()
            big_s, big_s1, big_s2 = np.sin(n * theta) / n, np.cos(n * theta), -n * np.sin(n * theta)
            t, t1, t2 = 2.0 * (1.0 - np.cos(n * theta)) / n ** 2, 2.0 * np.sin(n * theta) / n, 2.0 * np.cos(n * theta)
            u, u1, u2 = 2.0 * np.sin(0.5 * e), np.cos(0.5 * e), -0.5 * np.sin(0.5 * e)
            s[rows] = big_s + u
            h_g = chart['curvature'][gid]
            extra += 0.5 * h_g * (t - big_s * big_s)
            if order >= 1:
                th1 = dphi[gid, rows]
                # ds_k/dq_l = S1 Theta_l + u1_k (delta_kl - Theta_l)
                jac[np.ix_(rows, rows)] = np.outer(big_s1 - u1, th1) + np.diag(u1)
                dextra[rows] += 0.5 * h_g * (t1 - 2.0 * big_s * big_s1) * th1
            if order >= 2:
                th2 = d2phi[gid][np.ix_(rows, rows)]
                diag_second[rows] = 0.0
                blocks.append((rows, th1, th2, big_s1, big_s2, u1, u2))
                d2extra[np.ix_(rows, rows)] += 0.5 * h_g * (
                    (t2 - 2.0 * big_s1 ** 2 - 2.0 * big_s * big_s2) * np.outer(th1, th1)
                    + (t1 - 2.0 * big_s * big_s1) * th2)
        g_q, h_q = chart['g_q'], chart['H_q']
        a = g_q + h_q @ s
        value = chart['energy'] + g_q @ s + 0.5 * s @ h_q @ s + extra
        diagnostics = {'theta': thetas, 'max_row_deviation': max_dev}
        if order == 0:
            return value, None, None, diagnostics
        first = jac.T @ a + dextra
        if order == 1:
            return value, first, None, diagnostics
        second = jac.T @ h_q @ jac + np.diag(a * diag_second) + d2extra
        for rows, th1, th2, big_s1, big_s2, u1, u2 in blocks:
            a_k = a[rows]
            sum_a = a_k.sum()
            block = (big_s2 * sum_a) * np.outer(th1, th1) + (big_s1 * sum_a - (a_k * u1).sum()) * th2
            # sum_k a_k u2_k (e_k - Theta)(e_k - Theta)^T with e_k the unit vector of row k
            w = a_k * u2
            block += np.diag(w) - np.outer(w, th1) - np.outer(th1, w) + w.sum() * np.outer(th1, th1)
            second[np.ix_(rows, rows)] += block
        return value, first, second, diagnostics

    # ----------------------------------------------------------- atlas weights

    def atlas_weights(self, phases, order=2):
        """
        omega_b, d omega_b / dq and d2 omega_b / dq2 for every chart, from
        the phases of the groups with images. Only the rotor rows of those
        groups carry derivatives.
        """

        phi, dphi, d2phi = phases
        nq = self.ic.n_rows
        n_charts = len(self.charts)
        omega = np.zeros(n_charts)
        domega = np.zeros((n_charts, nq)) if order >= 1 else None
        d2omega = np.zeros((n_charts, nq, nq)) if order >= 2 else None
        for b, chart in enumerate(self.charts):
            factors = []
            for gid in self.atlas_groups:
                theta = wrap_angle(phi[gid] - chart['phi0'][gid])
                factors.append((gid, *_taper(theta, self.delta, order)))
            values = [f[1] for f in factors]
            omega[b] = float(np.prod(values)) if values else 1.0
            if omega[b] == 0.0 or order == 0:
                continue
            for i, (gid, val, first, second) in enumerate(factors):
                rest = float(np.prod(values[:i] + values[i + 1:])) if len(values) > 1 else 1.0
                th1 = dphi[gid]
                domega[b] += rest * first * th1
                if order >= 2:
                    d2omega[b] += rest * (second * np.outer(th1, th1) + first * d2phi[gid])
                    for j, (hid, val_h, first_h, _) in enumerate(factors):
                        if j == i:
                            continue
                        rest_ij = float(np.prod([values[k] for k in range(len(values)) if k not in (i, j)])) \
                            if len(values) > 2 else 1.0
                        d2omega[b] += rest_ij * first * first_h * np.outer(th1, dphi[hid])
        return omega, domega, d2omega

    # --------------------------------------------------------------- evaluate

    def evaluate(self, x, order=2, qbb2=None):
        """
        Q_a at the Cartesian geometry x (N, 3): a dict with 'energy',
        'gradient' (N, 3), 'hessian' (3N, 3N) (None above `order`),
        'weights' (lambda_b), 'active' (chart indices with omega_b > 0),
        'theta' (per chart, per group), 'phase_ratio' (|Z_g|/|K_g|) and
        'max_row_deviation'. qbb2 may pass precomputed (q, B, B2).
        """

        x = np.asarray(x, dtype=np.float64).reshape(self.n_atoms, 3)
        if qbb2 is None:
            if order >= 2:
                q, b_mat, b2 = self.ic.values_b_b2(x)
            else:
                q, b_mat = self.ic.values_and_b(x)
                b2 = None
        else:
            q, b_mat, b2 = qbb2
        result = self.evaluate_internal(q, order)
        out = {k: result[k] for k in ('energy', 'weights', 'active', 'theta', 'phase_ratio', 'max_row_deviation')}
        out['gradient'] = (b_mat.T @ result['dq']).reshape(self.n_atoms, 3) if order >= 1 else None
        out['hessian'] = None
        if order >= 2:
            out['hessian'] = b_mat.T @ result['d2q'] @ b_mat + np.einsum('k,kxy->xy', result['dq'], b2)
        out['dq'] = result['dq']
        out['d2q'] = result['d2q']
        return out

    def evaluate_internal(self, q, order=2):
        """
        Q_a and its derivatives with respect to the internal coordinates q.
        """

        model = self.model
        phases = model.phases(q, order)
        phi, dphi, d2phi, ratio = phases
        omega, domega, d2omega = self.atlas_weights((phi, dphi, d2phi), order)
        total = omega.sum()
        assert_msg_critical(total > 0.0, 'grouped MSI reference: no chart covers this geometry (every atlas weight is zero).')
        active = [b for b in range(len(self.charts)) if omega[b] > 0.0]
        lam = omega / total
        nq = self.ic.n_rows
        values = np.zeros(len(self.charts))
        firsts = {}
        seconds = {}
        thetas = np.zeros((len(self.charts), model.n_groups))
        max_dev = np.zeros(model.n_groups)
        for b in active:
            val, first, second, diag = self.chart_value(self.charts[b], q, (phi, dphi, d2phi), order)
            values[b] = val
            firsts[b] = first
            seconds[b] = second
            thetas[b] = diag['theta']
            max_dev = np.maximum(max_dev, diag['max_row_deviation'])
        energy = float(lam @ values)
        out = {'energy': energy, 'weights': lam, 'active': active, 'theta': thetas, 'phase_ratio': ratio,
               'max_row_deviation': max_dev, 'dq': None, 'd2q': None}
        if order == 0:
            return out
        dtotal = domega.sum(axis=0)
        dlam = {b: (domega[b] - lam[b] * dtotal) / total for b in active}
        dq = np.zeros(nq)
        for b in active:
            dq += lam[b] * firsts[b] + values[b] * dlam[b]
        out['dq'] = dq
        if order == 1:
            return out
        d2total = d2omega.sum(axis=0)
        d2q = np.zeros((nq, nq))
        for b in active:
            d2lam = (d2omega[b] - np.outer(dlam[b], dtotal) - np.outer(dtotal, dlam[b]) - lam[b] * d2total) / total
            d2q += lam[b] * seconds[b] + np.outer(dlam[b], firsts[b]) + np.outer(firsts[b], dlam[b]) + values[b] * d2lam
        out['d2q'] = d2q
        return out

    def residual_data(self, x, energy, gradient, hessian=None, project=True, qbb2=None):
        """
        R = U - Q_a, g_R = P g_U - grad Q_a and H_R = H_U^sym - grad2 Q_a at x
        (ADR-007, node data step 1). H_R is not projected: its rigid-body part
        (H_R r = omega x g_R for a rotation mode) belongs to every invariant
        function with gradient g_R, and the complete-coordinate split is
        insensitive to it.
        """

        x = np.asarray(x, dtype=np.float64).reshape(self.n_atoms, 3)
        ref = self.evaluate(x, order=2 if hessian is not None else 1, qbb2=qbb2)
        g = np.asarray(gradient, dtype=np.float64).reshape(-1)
        if project:
            g = rigid_body_projector(x) @ g
        residual = float(energy) - ref['energy']
        g_r = (g - ref['gradient'].reshape(-1)).reshape(self.n_atoms, 3)
        h_r = None
        if hessian is not None:
            h = np.asarray(hessian, dtype=np.float64)
            h_r = 0.5 * (h + h.T) - ref['hessian']
        return residual, g_r, h_r, ref
