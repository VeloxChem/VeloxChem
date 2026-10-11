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
Phase-constrained relaxation of the grouped MSI v2 model family: the cut
policy `relaxed_rotor_v1` (ADR-003 amendment 3; PROMPT-M1.md section 3.4).

A node is the minimum of U over the relaxation zone subject to the exact
constraints phi_g(x) = target_g on every group, started from the rigid cut.
The minimiser works in Cartesian space on any callable that returns energy
and gradient:

- the constraints are the phases (analytic gradient from the coordinate
  model), the frozen rows outside the relaxation zone (their B rows), and
  the six rigid-body modes;
- every step is a quasi-Newton (BFGS) step in the tangent space of the
  constraints, found from the KKT system of the current curvature estimate,
  with a backtracking line search on the energy;
- after every trial step the constraints are restored by Newton iterations
  with the constraint Jacobian, so that every accepted iterate satisfies the
  phases to the persisted threshold.

The same module gives the valley tangent of a converged node from the
bordered Lagrangian system (ADR-007), the continuity check between
neighbouring nodes, the anchor-minimum assertion and the invariance check of
the relaxation zone under the group permutations.

M2 (ADR-003 amendment 6, PROMPT-M2 decision D3): the branch diagnostics
between neighbouring relaxed nodes (secant against tangent, energy against
multipliers, tracked soft modes, the ratio of the valley tangent norms),
which warn and are recorded, and the investigation procedure (the
bidirectional warm-start sweep between two nodes and its classification),
whose "reject" verdict stops a construction.
"""

from dataclasses import dataclass, field

import numpy as np

from .errorhandler import assert_msg_critical
from .groupedcoordinates import wrap_angle, BOND, ANGLE


@dataclass
class RelaxationResult:
    coordinates: np.ndarray
    energy: float
    gradient: np.ndarray
    targets: np.ndarray
    multipliers: np.ndarray
    converged: bool
    iterations: int
    evaluations: int
    projected_gradient_norm: float
    constraint_residual: float
    step_norm: float
    history: list = field(default_factory=list)


def rigid_body_modes(x):
    """
    An orthonormal basis (6, 3N) of the translations and the infinitesimal
    rotations about the centroid of x.
    """

    x = np.asarray(x, dtype=np.float64).reshape(-1, 3)
    n = x.shape[0]
    centred = x - x.mean(axis=0)
    modes = []
    for k in range(3):
        v = np.zeros((n, 3))
        v[:, k] = 1.0
        modes.append(v.reshape(-1))
    for k in range(3):
        axis = np.zeros(3)
        axis[k] = 1.0
        modes.append(np.cross(axis, centred).reshape(-1))
    q, _ = np.linalg.qr(np.array(modes).T)
    return q.T


class PhaseConstrainedMinimizer:
    """
    :param model: the GroupedCoordinateModel (phases, B rows)
    :param energy_and_gradient: callable x (N, 3) -> (energy, gradient (N, 3))
    :param frozen_rows: internal-coordinate rows held at their starting
        values (the complement of the relaxation zone among the non-rotor
        rows); None or empty for a full relaxation
    :param policy: the GroupedModelPolicyV2 with the thresholds; None takes
        the model's policy
    """

    def __init__(self, model, energy_and_gradient, frozen_rows=None, policy=None):
        self.model = model
        self.energy_and_gradient = energy_and_gradient
        self.frozen_rows = np.array(sorted(set(int(r) for r in (frozen_rows or ()))), dtype=int)
        rotor = set(model.rotor_rows)
        assert_msg_critical(
            not (set(self.frozen_rows.tolist()) & rotor),
            'grouped MSI relaxation: rotor rows cannot be frozen; they are set by the phase constraints.')
        self.policy = policy if policy is not None else model.policy
        self.n_atoms = model.n_atoms
        self.evaluations = 0

    # ----------------------------------------------------------- constraints

    def constraints(self, x, targets, frozen_values):
        """
        c(x) (phases minus targets, wrapped; frozen rows minus their values;
        zeros for the rigid-body modes) and the Jacobian (n_c, 3N).
        """

        q, b_mat = self.model.ic.values_and_b(x)
        phi, dphi, _, _ = self.model.phases(q, order=1)
        rows = [wrap_angle(phi - targets)]
        jac = [dphi @ b_mat]
        if self.frozen_rows.size:
            rows.append(q[self.frozen_rows] - frozen_values)
            jac.append(b_mat[self.frozen_rows])
        rigid = rigid_body_modes(x)
        rows.append(np.zeros(rigid.shape[0]))
        jac.append(rigid)
        c = np.concatenate([np.atleast_1d(r) for r in rows])
        return c, np.vstack(jac), q

    def restore(self, x, targets, frozen_values, tol, max_iter=20):
        """
        Newton restoration of the constraints: x <- x - A^T (A A^T)^-1 c.
        """

        x = np.array(x, dtype=np.float64)
        for _ in range(max_iter):
            c, jac, _ = self.constraints(x, targets, frozen_values)
            residual = np.abs(c).max()
            if residual < tol:
                return x, residual
            step = jac.T @ np.linalg.solve(jac @ jac.T, c)
            x = x - step.reshape(self.n_atoms, 3)
        c, _, _ = self.constraints(x, targets, frozen_values)
        return x, float(np.abs(c).max())

    def _evaluate(self, x):
        self.evaluations += 1
        energy, gradient = self.energy_and_gradient(x)
        return float(energy), np.asarray(gradient, dtype=np.float64).reshape(-1)

    # -------------------------------------------------------------- minimise

    def minimize(self, x_start, targets=None, initial_hessian=None, max_iterations=None, record_history=False):
        """
        Minimises U subject to phi_g = target_g (targets None: the phases of
        the start geometry) and the frozen rows at their start values.
        Returns a RelaxationResult; convergence needs the projected gradient,
        the constraint residual and the step below the policy thresholds.
        """

        policy = self.policy
        tol_g = policy.relaxation_gradient_threshold
        tol_c = policy.relaxation_constraint_threshold
        tol_s = policy.relaxation_step_threshold
        n_iter = policy.relaxation_max_iterations if max_iterations is None else int(max_iterations)
        x = np.array(x_start, dtype=np.float64).reshape(self.n_atoms, 3)
        q0 = self.model.ic.values(x)
        frozen_values = q0[self.frozen_rows] if self.frozen_rows.size else np.zeros(0)
        if targets is None:
            targets = self.model.phases(q0, order=0)[0]
        targets = np.asarray(targets, dtype=np.float64)
        x, residual = self.restore(x, targets, frozen_values, tol_c)
        self.evaluations = 0  # per call
        n = 3 * self.n_atoms
        hess = np.eye(n)
        if initial_hessian is not None:
            hess = np.array(initial_hessian, dtype=np.float64)
            hess = 0.5 * (hess + hess.T)
        energy, grad = self._evaluate(x)
        history = []
        converged = False
        step_norm = 0.0
        proj_norm = np.inf
        iteration = 0
        # the line search tests the energy (Armijo) until that test first fails, which happens where the expected
        # decrease is below the noise of the energy (DFT grids and SCF convergence; found at PBE in M2); from then
        # on a step is accepted when it lowers the projected gradient norm (the criterion of the Newton polish of
        # M1-C), and a stall of that norm ends the minimisation. On noise-free surfaces the energy test never
        # fails before convergence and the behaviour of M1 is unchanged.
        gradient_mode = False
        norms = []
        for iteration in range(1, n_iter + 1):
            c, jac, _ = self.constraints(x, targets, frozen_values)
            # multipliers and the projected gradient: g_t = g - A^T lambda
            lam = np.linalg.lstsq(jac.T, grad, rcond=None)[0]
            g_t = grad - jac.T @ lam
            proj_norm = float(np.abs(g_t).max())
            if record_history:
                history.append((iteration, energy, proj_norm, float(np.abs(c).max()), step_norm))
            if proj_norm < tol_g and np.abs(c).max() < tol_c:
                converged = True
                break
            if gradient_mode:
                norms.append(proj_norm)
                if len(norms) > 10 and min(norms[-10:]) > 0.9 * norms[-11]:
                    break   # stalled at the noise floor of the gradient: not converged
            # KKT step with the current curvature estimate: minimise 1/2 p^T H p + g^T p, A p = 0
            n_c = jac.shape[0]
            kkt = np.zeros((n + n_c, n + n_c))
            kkt[:n, :n] = hess
            kkt[:n, n:] = jac.T
            kkt[n:, :n] = jac
            rhs = np.concatenate([-grad, np.zeros(n_c)])
            try:
                sol = np.linalg.solve(kkt, rhs)
            except np.linalg.LinAlgError:
                sol = np.linalg.lstsq(kkt, rhs, rcond=None)[0]
            p = sol[:n]
            if p @ g_t > 0.0:
                p = -g_t  # not a descent direction: fall back to steepest descent in the tangent space
            # backtracking line search on the restored point
            alpha = 1.0
            accepted = False
            if not gradient_mode:
                for _ in range(40):
                    trial = x + alpha * p.reshape(self.n_atoms, 3)
                    trial, residual = self.restore(trial, targets, frozen_values, tol_c)
                    e_trial, g_trial = self._evaluate(trial)
                    if e_trial <= energy + 1.0e-4 * alpha * (g_t @ p):
                        accepted = True
                        break
                    alpha *= 0.5
                    if alpha < 1.0e-3 and abs(e_trial - energy) < 1.0e-8:
                        break   # the energy no longer resolves the step (below the SCF noise): switch to the gradient criterion
                if not accepted:
                    gradient_mode = True
            if gradient_mode and not accepted:
                alpha = 1.0
                for _ in range(4):
                    trial = x + alpha * p.reshape(self.n_atoms, 3)
                    trial, residual = self.restore(trial, targets, frozen_values, tol_c)
                    e_trial, g_trial = self._evaluate(trial)
                    c_t, jac_t, _ = self.constraints(trial, targets, frozen_values)
                    lam_t = np.linalg.lstsq(jac_t.T, g_trial, rcond=None)[0]
                    if float(np.abs(g_trial - jac_t.T @ lam_t).max()) < proj_norm:
                        accepted = True
                        break
                    alpha *= 0.5
            if not accepted:
                break
            s = (trial - x).reshape(-1)
            y = g_trial - grad
            step_norm = float(np.linalg.norm(s))
            # BFGS update of the curvature estimate (skipped unless s.y > 0)
            sy = s @ y
            if sy > 1.0e-12 * max(1.0, np.linalg.norm(s) * np.linalg.norm(y)):
                if iteration == 1 and initial_hessian is None:
                    hess = np.eye(n) * (sy / (s @ s))
                hs = hess @ s
                hess = hess - np.outer(hs, hs) / (s @ hs) + np.outer(y, y) / sy
            x, energy, grad = trial, e_trial, g_trial
            if step_norm < tol_s:
                break
        # convergence is the KKT condition at the final iterate: projected
        # gradient and constraint residual below their thresholds
        c, jac, _ = self.constraints(x, targets, frozen_values)
        lam = np.linalg.lstsq(jac.T, grad, rcond=None)[0]
        g_t = grad - jac.T @ lam
        proj_norm = float(np.abs(g_t).max())
        converged = proj_norm < tol_g and np.abs(c).max() < tol_c
        return RelaxationResult(
            coordinates=x, energy=energy, gradient=grad.reshape(self.n_atoms, 3), targets=targets,
            multipliers=lam[:len(targets)], converged=converged, iterations=iteration,
            evaluations=self.evaluations, projected_gradient_norm=proj_norm,
            constraint_residual=float(np.abs(c).max()), step_norm=step_norm, history=history)

    # ---------------------------------------------------------- valley data

    def valley_tangent(self, x, hessian, targets=None):
        """
        The tangents t_h = dx*/dphi_h (G, 3N) of the relaxed valley at a
        converged node x, from the bordered Lagrangian system

            K t_h = A^T mu_h,  A t_h = e_h,   K = H_U - sum_g lambda_g grad2 phi_g - sum_k lambda_k B2_k,

        with the multipliers of grad U = A^T lambda at x and the rigid-body
        modes among the constraints (right-hand side zero). Returns (t, w')
        with w' = B_v t the derivative of the non-rotor rows (G, n_v).
        `targets` is accepted for symmetry with minimize and not needed: the
        multipliers are determined by the gradient at x alone.
        """

        model = self.model
        x = np.asarray(x, dtype=np.float64).reshape(self.n_atoms, 3)
        q, b_mat, b2 = model.ic.values_b_b2(x)
        phi, grad_phi, hess_phi, _ = model.phases_cartesian(x, order=2, qbb2=(q, b_mat, b2))
        if targets is None:
            targets = phi
        energy, grad = self._evaluate(x)
        jac = [grad_phi]
        if self.frozen_rows.size:
            jac.append(b_mat[self.frozen_rows])
        rigid = rigid_body_modes(x)
        jac.append(rigid)
        jac = np.vstack(jac)
        lam = np.linalg.lstsq(jac.T, grad, rcond=None)[0]
        n_g = model.n_groups
        k_mat = np.asarray(hessian, dtype=np.float64)
        k_mat = 0.5 * (k_mat + k_mat.T) - np.einsum('g,gxy->xy', lam[:n_g], hess_phi)
        if self.frozen_rows.size:
            k_mat = k_mat - np.einsum('k,kxy->xy', lam[n_g:n_g + self.frozen_rows.size], b2[self.frozen_rows])
        n = 3 * self.n_atoms
        n_c = jac.shape[0]
        kkt = np.zeros((n + n_c, n + n_c))
        kkt[:n, :n] = k_mat
        kkt[:n, n:] = jac.T
        kkt[n:, :n] = jac
        tangents = np.zeros((n_g, n))
        for h in range(n_g):
            rhs = np.zeros(n + n_c)
            rhs[n + h] = 1.0
            sol = np.linalg.solve(kkt, rhs)
            tangents[h] = sol[:n]
        w_prime = tangents @ b_mat[model.response_rows].T
        return tangents, w_prime


def well_posedness(minimizer, x, hessian, targets, n_modes=4):
    """
    The well-posedness diagnostics of a relaxed node (ADR-003 amendment of
    M1-C, decision 4): with A the constraint Jacobian (phases, frozen rows,
    rigid modes), lambda the multipliers of grad U = A^T lambda and
    K = H_U - sum_g lambda_g grad2 phi_g - sum_k lambda_k B2_k the Hessian
    of the Lagrangian, returns a dict with 'rank_A' and 'n_constraints'
    (LICQ holds when equal), 'sigma_min_A', 'lambda_min' of Z^T K Z (SOSC
    holds when positive), 'rows' (the three internal rows with the largest
    weight in its eigenvector, with their values), 'bordered_condition',
    'multiplier_norm', 'projected_gradient' and 'constraint_residual'.

    M2 (D3): also 'phase_multipliers' (the multipliers of the phase
    constraints, one per group) and the lowest `n_modes` eigenpairs of
    Z^T K Z as 'soft_mode_values' (k,) and 'soft_mode_rows' (k, n_rows),
    the eigenvectors expressed in internal rows (B v / |B_k|, normalised),
    for the mode tracking between neighbouring nodes.
    """

    model = minimizer.model
    x = np.asarray(x, dtype=np.float64).reshape(model.n_atoms, 3)
    q, b_mat, b2 = model.ic.values_b_b2(x)
    phi, grad_phi, hess_phi, _ = model.phases_cartesian(x, order=2, qbb2=(q, b_mat, b2))
    frozen = minimizer.frozen_rows
    _, grad = minimizer.energy_and_gradient(x)
    grad = np.asarray(grad, dtype=np.float64).reshape(-1)
    jac = [grad_phi]
    if frozen.size:
        jac.append(b_mat[frozen])
    jac.append(rigid_body_modes(x))
    jac = np.vstack(jac)
    lam = np.linalg.lstsq(jac.T, grad, rcond=None)[0]
    n_g = model.n_groups
    k_mat = 0.5 * (np.asarray(hessian) + np.asarray(hessian).T) - np.einsum('g,gxy->xy', lam[:n_g], hess_phi)
    if frozen.size:
        k_mat = k_mat - np.einsum('k,kxy->xy', lam[n_g:n_g + frozen.size], b2[frozen])
    _, s, vt = np.linalg.svd(jac, full_matrices=True)
    rank = int((s > 1.0e-10 * s[0]).sum())
    z = vt[rank:].T
    reduced = z.T @ k_mat @ z
    reduced = 0.5 * (reduced + reduced.T)
    evals, evecs = np.linalg.eigh(reduced)
    vec = z @ evecs[:, 0]
    row_norm = np.linalg.norm(b_mat, axis=1)
    weight = np.abs(b_mat @ vec) / np.where(row_norm > 0.0, row_norm, 1.0)
    top = np.argsort(weight)[::-1][:3]
    k_modes = min(int(n_modes), len(evals))
    # every eigenvector of Z^T K Z expressed in internal rows (B v / |B_k|, normalised): the lowest k are persisted
    # as the soft modes, the full set serves the tracking between neighbouring nodes (transient)
    all_rows = (b_mat @ (z @ evecs)).T / np.where(row_norm > 0.0, row_norm, 1.0)[None, :]
    norms = np.linalg.norm(all_rows, axis=1)
    all_rows = all_rows / np.where(norms > 0.0, norms, 1.0)[:, None]
    mode_rows = [all_rows[a] for a in range(k_modes)]
    n = 3 * model.n_atoms
    bordered = np.zeros((n + jac.shape[0], n + jac.shape[0]))
    bordered[:n, :n] = k_mat
    bordered[:n, n:] = jac.T
    bordered[n:, :n] = jac
    sv_b = np.linalg.svd(bordered, compute_uv=False)
    return {
        'rank_A': rank, 'n_constraints': int(jac.shape[0]), 'sigma_min_A': float(s[rank - 1]),
        'lambda_min': float(evals[0]), 'rows': [(model.describe_row(int(r), q), float(weight[r])) for r in top],
        'bordered_condition': float(sv_b[0] / sv_b[-1]) if sv_b[-1] > 0.0 else np.inf,
        'multiplier_norm': float(np.linalg.norm(lam)),
        'projected_gradient': float(np.abs(grad - jac.T @ lam).max()),
        'constraint_residual': float(np.abs(wrap_angle(phi - np.asarray(targets))).max()),
        'phase_multipliers': np.array(lam[:n_g], dtype=np.float64),
        'soft_mode_values': np.array(evals[:k_modes], dtype=np.float64),
        'soft_mode_rows': np.array(mode_rows, dtype=np.float64).reshape(k_modes, -1),
        'mode_values_all': np.array(evals, dtype=np.float64),
        'mode_rows_all': all_rows,
    }


def assert_well_posed(minimizer, x, hessian, targets, policy=None, label='the node'):
    """
    Decision 4 of M1-C: a relaxed node (and the anchor) must be a strict
    constrained minimum, LICQ (full row rank of A) and SOSC with the
    persisted margin lambda_min(Z^T K Z) >= relaxation_curvature_min. By the
    implicit-function theorem these are exactly the conditions that make
    x*(phi) a smooth function and the valley tangent well defined. Returns
    the diagnostics; a violation is a critical error naming the node, the
    rows and the remedy.
    """

    policy = policy if policy is not None else minimizer.policy
    diag = well_posedness(minimizer, x, hessian, targets, n_modes=policy.branch_tracked_modes)
    assert_msg_critical(
        diag['rank_A'] == diag['n_constraints'],
        f'grouped MSI relaxation: the constraints at {label} are linearly dependent (rank {diag["rank_A"]} of '
        f'{diag["n_constraints"]}, sigma_min(A) = {diag["sigma_min_A"]:.1e}): the phase gradients, frozen rows and rigid modes '
        'do not define the node (LICQ fails); check the relaxation zone and the phase validity.')
    rows_txt = ', '.join(f'{name} ({weight:.2f})' for name, weight in diag['rows'])
    assert_msg_critical(
        diag['lambda_min'] >= policy.relaxation_curvature_min,
        f'grouped MSI relaxation: {label} is not a strict constrained minimum: the smallest curvature of the Lagrangian '
        f'on the constraint tangent space is {diag["lambda_min"]:+.3e} hartree/bohr^2 (limit relaxation_curvature_min = '
        f'{policy.relaxation_curvature_min:.0e}), softest direction on {rows_txt}. The valley bifurcates or is degenerate there: '
        'refine or shift the grid, declare a domain that excludes it, or freeze the responsible rows.')
    return diag


def node_continuity(model, x_a, x_b):
    """
    The largest change of the non-rotor rows between two node geometries:
    (max bond length change in bohr, max angular change in degrees).
    """

    qa = model.ic.values(x_a)
    qb = model.ic.values(x_b)
    rows = model.response_rows
    kinds = model.ic.kinds[rows]
    d = qb[rows] - qa[rows]
    bond = kinds == BOND
    bond_change = float(np.abs(1.0 / qb[rows][bond] - 1.0 / qa[rows][bond]).max()) if bond.any() else 0.0
    angular = kinds >= ANGLE
    ang_change = float(np.degrees(np.abs(wrap_angle(d[angular]))).max()) if angular.any() else 0.0
    return bond_change, ang_change


def assert_node_continuity(model, x_a, x_b, policy=None, label=''):
    policy = policy if policy is not None else model.policy
    bond_change, ang_change = node_continuity(model, x_a, x_b)
    assert_msg_critical(
        bond_change <= policy.node_continuity_limit_bond_bohr and ang_change <= policy.node_continuity_limit_degrees,
        f'grouped MSI relaxation: the relaxed node {label} is not continuous with its neighbour '
        f'(bond change {bond_change:.3f} bohr, angular change {ang_change:.1f} degrees).')
    return bond_change, ang_change


def assert_anchor_is_minimum(model, energy_and_gradient, policy=None):
    """
    The anchor must be the constrained minimum at its own phases: the
    gradient projected on the tangent space of the phase constraints (and
    the rigid-body modes) must be below the gradient threshold.
    """

    policy = policy if policy is not None else model.policy
    x = model.x_anchor
    _, grad = energy_and_gradient(x)
    grad = np.asarray(grad, dtype=np.float64).reshape(-1)
    q, b_mat = model.ic.values_and_b(x)
    _, dphi, _, _ = model.phases(q, order=1)
    jac = np.vstack([dphi @ b_mat, rigid_body_modes(x)])
    lam = np.linalg.lstsq(jac.T, grad, rcond=None)[0]
    residual = float(np.abs(grad - jac.T @ lam).max())
    assert_msg_critical(
        residual < policy.relaxation_gradient_threshold,
        f'grouped MSI relaxation: the anchor is not a minimum at its own phases (projected gradient '
        f'{residual:.2e} above {policy.relaxation_gradient_threshold:.1e}); the relaxed cut policy needs an optimised anchor.')
    return residual


def check_zone_invariance(model, zone_rows):
    """
    The relaxation zone (non-rotor rows allowed to move) must be mapped onto
    itself by the row permutation of every group with images, so that
    P x*(phi) = x*(phi + sense 2 pi / m) holds exactly. Returns the groups
    that break it.
    """

    zone = set(int(r) for r in zone_rows)
    broken = []
    for gid, g in enumerate(model.groups):
        if g.images == 1:
            continue
        pi = model.row_permutation(gid, 1)
        if {int(pi[r]) for r in zone} != zone:
            broken.append(gid)
    return broken


def assert_zone_invariance(model, zone_rows):
    broken = check_zone_invariance(model, zone_rows)
    assert_msg_critical(
        not broken,
        f'grouped MSI relaxation: the relaxation zone is not invariant under the permutation of groups {broken}.')


# ------------------------------------------------------------ branch diagnostics (M2, D3)

TWO_PI = 2.0 * np.pi


def _wrapped_row_difference(model, qa, qb):
    """q_b - q_a over the non-rotor rows, wrapped on angular rows."""
    rows = model.response_rows
    d = np.asarray(qb, dtype=np.float64) - np.asarray(qa, dtype=np.float64)
    periodic = model.ic.kinds[rows] >= 2
    d[periodic] = wrap_angle(d[periodic])
    return d


def image_record(model, record, gid, power=1):
    """
    The virtual image of a node record one period of group gid up (ADR-004
    amendment 3): the geometry permuted, the offsets shifted by the period,
    the covariant node data transformed (w -> Pi w + shift, w1 -> Pi w1, the
    soft-mode rows permuted); energy, residual and multipliers unchanged.
    Used for the wrap-around pair of a closed line of a group with images.
    """

    g = model.groups[gid]
    if g.images == 1:
        return record
    out = dict(record)
    pi = model.row_permutation(gid, power)
    rows = model.response_rows
    position = {int(r): i for i, r in enumerate(rows)}
    pv = np.array([position[int(pi[r])] for r in rows], dtype=int)
    shift = model.q_anchor[rows][pv] - model.q_anchor[rows]
    periodic = model.ic.kinds[rows] >= 2
    shift[periodic] = wrap_angle(shift[periodic])
    data = dict(record['data'])
    data['w'] = np.asarray(record['data']['w'])[pv] + shift
    data['w1'] = np.asarray(record['data']['w1'])[:, pv]
    out['data'] = data
    offsets = np.array(record['offsets'], dtype=np.float64)
    offsets[gid] += power * TWO_PI / g.images
    out['offsets'] = offsets
    p = model.cyclic_permutation(gid, power)
    out['x'] = np.asarray(record['x'])[p]
    wp = dict(record['well_posedness'])
    for key in ('soft_mode_rows', 'mode_rows_all'):
        if key in wp:
            wp[key] = np.asarray(wp[key])[:, pi]
    out['well_posedness'] = wp
    out['image_of'] = record
    return out


def branch_lines(model, records):
    """
    The lines of neighbouring nodes along every member direction: for every
    group g and every combination of the other groups' offsets that occurs
    in the node records, the nodes sorted by their offset in g, with the
    wrap-around pair when the grid closes the circle (the first node's
    virtual image one period up, `image_record`, for a group with images).
    The singleton path of g is the line through the anchor (all other
    offsets zero); the lines of a multi-member grid start at a face node.
    Returns a list of dicts with 'gid', 'base' (the other offsets), 'nodes'
    [(offset, record)], 'period', 'closed'.
    """

    lines = {}
    for record in records.values():
        off = np.asarray(record['offsets'], dtype=np.float64)
        for gid in range(model.n_groups):
            base = off.copy()
            base[gid] = 0.0
            key = (gid, tuple(int(np.rint(o * 1.0e8)) for o in base))
            lines.setdefault(key, []).append((float(off[gid]), record))
    out = []
    for (gid, base_key), items in sorted(lines.items()):
        if len(items) < 2:
            continue
        period = TWO_PI / model.groups[gid].images
        items = sorted(((o % period, r) for o, r in items), key=lambda t: t[0])
        offsets = np.array([o for o, _ in items])
        gaps = np.diff(offsets)
        wrap = offsets[0] + period - offsets[-1]
        closed = bool(abs(wrap - gaps[0]) < 1.0e-9) if len(gaps) else False
        out.append({'gid': gid, 'base': np.array(base_key, dtype=np.float64) * 1.0e-8, 'nodes': items,
                    'period': period, 'closed': closed})
    return out


def _follow_modes(tracked_rows, rows_b):
    """
    Follows tracked modes (their eigenvectors in internal rows at the
    previous node) to the next node: every tracked mode is matched to the
    eigenvector of the full spectrum at b with the largest overlap. Returns
    the matched indices and the overlaps.
    """

    overlap = np.abs(np.asarray(tracked_rows) @ np.asarray(rows_b).T)
    match = np.argmax(overlap, axis=1)
    return match, overlap[np.arange(len(match)), match]


def branch_diagnostics(model, records, policy=None, with_tangent=True):
    """
    The branch diagnostics of ADR-003 amendment 6 between neighbouring nodes
    along every line of `branch_lines` (relaxed policy; the records carry
    'offsets', 'x', 'energy', 'data' with w and w1, and 'well_posedness'
    with the phase multipliers and the soft modes):

    1. secant against tangent, s = |(w_{j+1} - w_j) - dphi (w1_j + w1_{j+1}) / 2|
       / max(|w_{j+1} - w_j|, floor) (angular rows wrapped; the floor is a
       tenth of the valley's excursion along the line). On a smooth branch
       the mismatch is the trapezoid error dphi^3 w3 / 12, with w1 and w3
       the first and third derivative of the valley position; for a branch
       with harmonic n that is of size x_w = dphi^3 n^2 max|w1| / 12, and
       the warning compares the mismatch with x_w ('secant_tangent_rel');
       a flip gives an O(1) mismatch;
    2. energy against multipliers, e = |(E_{j+1} - E_j) - dphi (lambda_j +
       lambda_{j+1}) / 2| / max(|E_{j+1} - E_j|, floor), lambda the
       multiplier of the moving member's phase constraint (dE*/dphi =
       lambda); the same expectation with x_e = dphi^3 n^2 max|lambda| / 12
       and the warning on the mismatch relative to x_e
       ('energy_multiplier_rel');
    3. tracked soft modes: the lowest k eigenpairs of Z^T K Z at the line
       start (the anchor on a singleton path, the face node on a line of a
       multi-member grid), each followed from node to node by the
       eigenvector with the largest overlap in internal rows over the full
       spectrum; the smallest ratio of a tracked eigenvalue to its value at
       the start, the crossings (a tracked mode's rank differs from its
       start rank) and the smallest overlap of the matching (below 0.5 the
       tracking is ambiguous, reported);
    4. the ratio of |w1| (the derivative along the line) to the anchor's
       value for the same direction (the growth of the valley tangent; the
       ratio between neighbours is reported but not warned on, because |w1|
       passes through zero at the turning points of the valley).

    Every node record gets 'branch': the worst value of every diagnostic
    over its neighbour pairs, and 'warnings'. Returns (warnings, lines):
    the warnings as dicts naming the line, the pair of nodes (offsets in
    degrees), the diagnostic, its value and the threshold, and the per-line
    tables. Thresholds and floors are the persisted branch_* parameters.
    Diagnostics 1 and 4 need the valley tangent w1 of the node data; with
    `with_tangent` False (a model whose w bank has order 0, or without the
    response) they are reported as zero and never warn.
    """

    policy = policy if policy is not None else model.policy
    lines = branch_lines(model, records)
    anchor_key = tuple(0 for _ in range(model.n_groups))
    anchor = records.get(anchor_key)
    warnings = []
    tables = []
    for record in records.values():
        record['branch'] = {'secant_tangent': 0.0, 'energy_multiplier': 0.0, 'mode_ratio_min': 1.0,
                            'mode_crossings': 0, 'tangent_ratio_neighbours': 1.0, 'tangent_ratio_anchor': 1.0,
                            'warnings': []}
    for line in lines:
        gid = line['gid']
        nodes = line['nodes']
        period = line['period']
        pairs = list(zip(nodes[:-1], nodes[1:]))
        if line['closed']:
            pairs.append((nodes[-1], (nodes[0][0] + period, image_record(model, nodes[0][1], gid))))
        label = (f'group {gid} with the other offsets {np.degrees(line["base"]).round(1).tolist()} deg')
        table = {'label': label, 'gid': gid, 'base': line['base'], 'closed': line['closed'], 'rows': []}
        # the scales of the line: the valley's excursion from the line start, the energy range, the largest tangent
        # and multiplier (the trapezoid expectations of diagnostics 1 and 2 are built from the last two)
        w_start = nodes[0][1]['data']['w']
        scale_w = max(float(np.linalg.norm(_wrapped_row_difference(model, w_start, rec['data']['w']))) for _, rec in nodes)
        energies = [float(rec['energy']) for _, rec in nodes]
        scale_e = max(energies) - min(energies)
        scale_t = max(float(np.linalg.norm(rec['data']['w1'][gid])) for _, rec in nodes)
        scale_l = max(abs(float(rec['well_posedness']['phase_multipliers'][gid])) for _, rec in nodes)
        floor_w = 0.1 * scale_w + 1.0e-12
        floor_e = 0.1 * scale_e + 1.0e-12
        harmonic = max(int(model.groups[gid].harmonic), 1)
        # mode tracking from the line start: the lowest k modes, each followed by its best overlap
        start = nodes[0][1]['well_posedness']
        values0 = np.asarray(start['soft_mode_values'])
        tracked_rows = np.array(start['soft_mode_rows'])
        start_rank = np.arange(len(values0))
        anchor_norm = None
        if anchor is not None:
            anchor_norm = float(np.linalg.norm(anchor['data']['w1'][gid]))
        for (off_a, rec_a), (off_b, rec_b) in pairs:
            dphi = off_b - off_a
            x_w = dphi ** 3 * harmonic ** 2 * scale_t / 12.0 + 1.0e-300
            x_e = dphi ** 3 * harmonic ** 2 * scale_l / 12.0 + 1.0e-300
            da, db = rec_a['data'], rec_b['data']
            d_w = _wrapped_row_difference(model, da['w'], db['w'])
            secant = d_w - dphi * 0.5 * (da['w1'][gid] + db['w1'][gid])
            s_abs = float(np.linalg.norm(secant))
            s_val = s_abs / max(np.linalg.norm(d_w), floor_w)
            lam_a = rec_a['well_posedness']['phase_multipliers'][gid]
            lam_b = rec_b['well_posedness']['phase_multipliers'][gid]
            d_e = rec_b['energy'] - rec_a['energy']
            e_abs = float(abs(d_e - dphi * 0.5 * (lam_a + lam_b)))
            e_val = e_abs / max(abs(d_e), floor_e)
            wb = rec_b['well_posedness']
            match, overlaps = _follow_modes(tracked_rows, wb['mode_rows_all'])
            tracked_rows = np.array(wb['mode_rows_all'])[match]
            crossings = int((match != start_rank).sum())
            ratios = np.asarray(wb['mode_values_all'])[match] / values0
            ratio_min = float(ratios.min())
            n_a = float(np.linalg.norm(da['w1'][gid]))
            n_b = float(np.linalg.norm(db['w1'][gid]))
            floor_t = 0.1 * scale_t + 1.0e-12
            t_neigh = float(max(max(n_a, floor_t) / max(n_b, floor_t), max(n_b, floor_t) / max(n_a, floor_t)))
            t_anchor = float(max(n_a, n_b) / max(anchor_norm, floor_t)) if anchor_norm is not None else 1.0
            row = {'pair': (float(off_a), float(off_b)), 'secant_tangent': s_val, 'energy_multiplier': e_val,
                   'expected': (float(x_w), float(x_e)), 'secant_tangent_rel': s_abs / x_w, 'energy_multiplier_rel': e_abs / x_e,
                   'mode_ratio_min': ratio_min, 'mode_crossings': crossings, 'overlap_min': float(overlaps.min()),
                   'tangent_norms': (n_a, n_b), 'tangent_ratio_neighbours': t_neigh, 'tangent_ratio_anchor': t_anchor,
                   'lambda_min_ratio': float(wb['lambda_min'] / max(values0[0], 1.0e-300)),
                   'scales': (scale_w, scale_e, scale_t, scale_l)}
            table['rows'].append(row)
            pair_txt = f'{np.degrees(off_a):.1f} and {np.degrees(off_b):.1f} deg'
            if not with_tangent:
                s_val = s_abs = 0.0
                t_neigh = t_anchor = 1.0
                row.update({'secant_tangent': 0.0, 'secant_tangent_rel': 0.0, 'tangent_ratio_neighbours': 1.0, 'tangent_ratio_anchor': 1.0})
            checks = (('secant_tangent', s_abs / x_w, policy.branch_secant_tangent_max, 'above'),
                      ('energy_multiplier', e_abs / x_e, policy.branch_energy_multiplier_max, 'above'),
                      ('mode_ratio_min', ratio_min, policy.branch_mode_ratio_min, 'below'),
                      ('tangent_ratio_anchor', t_anchor, policy.branch_tangent_ratio_max, 'above'))
            for rec in (rec_a, rec_b.get('image_of', rec_b)):
                br = rec['branch']
                br['secant_tangent'] = max(br['secant_tangent'], s_abs / x_w)
                br['energy_multiplier'] = max(br['energy_multiplier'], e_abs / x_e)
                br['mode_ratio_min'] = min(br['mode_ratio_min'], ratio_min)
                br['mode_crossings'] += crossings
                br['tangent_ratio_neighbours'] = max(br['tangent_ratio_neighbours'], t_neigh)
                br['tangent_ratio_anchor'] = max(br['tangent_ratio_anchor'], t_anchor)
            for name, value, limit, sense in checks:
                bad = value > limit if sense == 'above' else value < limit
                if bad:
                    warning = {'line': label, 'gid': gid, 'base': line['base'], 'pair': (float(off_a), float(off_b)),
                               'diagnostic': name, 'value': value, 'limit': limit,
                               'message': (f'grouped MSI relaxation: branch diagnostic {name} = {value:.3g} '
                                           f'({sense} the limit {limit:g}) between the nodes at {pair_txt} on the path of {label}'
                                           + (f'; the mismatch relative to the trapezoid size {x_w if name == "secant_tangent" else x_e:.3g}'
                                              if name in ('secant_tangent', 'energy_multiplier') else ''))}
                    warnings.append(warning)
                    for rec in (rec_a, rec_b.get('image_of', rec_b)):
                        rec['branch']['warnings'].append(warning['message'])
        tables.append(table)
    return warnings, tables


def sweep_between_nodes(model, energy_gradient, policy, record_a, record_b, gid, step_degrees=None, frozen_rows=None):
    """
    The bidirectional warm-start sweep between two neighbouring nodes along
    member `gid`: from a to b and from b to a in steps of step_degrees, every
    relaxation started from the previous relaxed point rotated rigidly by
    the step. Returns a dict with the sweep angles (offsets of gid), the
    forward and backward results per angle (x, energy, converged,
    evaluations), the disagreement of the two sweeps at every common angle
    (max bond change in bohr, max angular change in degrees, energy
    difference), and the disagreement of each sweep's end point with the
    policy node there.
    """

    minimizer = PhaseConstrainedMinimizer(model, energy_gradient, frozen_rows=frozen_rows, policy=policy)
    step_deg = policy.branch_sweep_step_degrees if step_degrees is None else float(step_degrees)
    off_a = np.asarray(record_a['offsets'], dtype=np.float64)
    off_b = np.asarray(record_b['offsets'], dtype=np.float64)
    period = TWO_PI / model.groups[gid].images
    delta = (off_b[gid] - off_a[gid]) % period
    if delta == 0.0:
        delta = period
    n_steps = max(1, int(np.rint(delta / np.radians(step_deg))))
    step = delta / n_steps
    angles = [off_a[gid] + i * step for i in range(n_steps + 1)]

    def sweep(x_start, off_start, direction):
        out = []
        x_prev = np.array(x_start, dtype=np.float64)
        for i in range(1, n_steps + 1):
            off = off_start.copy()
            off[gid] = off_start[gid] + direction * i * step
            start = model.rigid_cut({gid: direction * step}, start=x_prev)
            res = minimizer.minimize(start, model.phi_anchor + off)
            out.append({'offset': float(off[gid]), 'x': res.coordinates, 'energy': float(res.energy),
                        'converged': bool(res.converged), 'evaluations': int(res.evaluations)})
            x_prev = res.coordinates
        return out

    forward = sweep(record_a['x'], off_a, +1)
    backward = sweep(record_b['x'], off_b, -1)
    disagreement = []
    for i in range(1, n_steps):
        f = forward[i - 1]
        b = backward[n_steps - i - 1]
        bond, ang = node_continuity(model, f['x'], b['x'])
        disagreement.append({'offset': angles[i], 'bond_change': bond, 'angle_change': ang,
                             'energy_difference': f['energy'] - b['energy']})
    end_forward = node_continuity(model, forward[-1]['x'], record_b['x'])
    end_backward = node_continuity(model, backward[-1]['x'], record_a['x'])
    return {'gid': gid, 'angles': angles, 'forward': forward, 'backward': backward, 'disagreement': disagreement,
            'end_forward': end_forward, 'end_backward': end_backward,
            'end_forward_energy': forward[-1]['energy'] - float(record_b['energy']),
            'end_backward_energy': backward[-1]['energy'] - float(record_a['energy']),
            'evaluations': sum(r['evaluations'] for r in forward + backward)}


def investigate_branch(model, energy_gradient, policy, record_a, record_b, gid, under_resolved=False,
                       step_degrees=None, frozen_rows=None):
    """
    The investigation procedure of ADR-003 amendment 6 for a warning
    between the nodes a and b along member gid: the bidirectional
    warm-start sweep, then the classification

    - 'continuous': the sweeps agree with each other and with the policy
      nodes within the continuity limits; verdict 'accept' (with the
      recommendation of a finer grid when `under_resolved`, that is when
      diagnostic 1 or 4 warned);
    - 'branch_flip': the sweeps disagree with each other or with a policy
      node beyond the limits (a flip or hysteresis); verdict 'reject', the
      remedies being a finer grid, freezing the rows of the bifurcating mode
      (relaxation_zone) or the rigid policy;
    - 'degenerate_mirror_branches': the sweeps end in different basins whose
      energies agree within branch_degenerate_energy_hartree; verdict
      'reject' unless the responsible rows are frozen.

    Returns a dict with 'classification', 'verdict', 'message', 'rows' (the
    non-rotor rows with the largest disagreement) and the sweep.
    """

    sweep = sweep_between_nodes(model, energy_gradient, policy, record_a, record_b, gid, step_degrees, frozen_rows)
    lim_b, lim_a = policy.node_continuity_limit_bond_bohr, policy.node_continuity_limit_degrees
    bad = [d for d in sweep['disagreement'] if d['bond_change'] > lim_b or d['angle_change'] > lim_a]
    end_bad = (sweep['end_forward'][0] > lim_b or sweep['end_forward'][1] > lim_a
               or sweep['end_backward'][0] > lim_b or sweep['end_backward'][1] > lim_a)
    label_a = f'{np.degrees(record_a["offsets"][gid]):.1f} deg'
    label_b = f'{np.degrees(record_b["offsets"][gid]):.1f} deg'
    rows = []
    if bad:
        worst = max(bad, key=lambda d: d['angle_change'] / lim_a + d['bond_change'] / lim_b)
        i = [d['offset'] for d in sweep['disagreement']].index(worst['offset'])
        f = sweep['forward'][i]
        b = sweep['backward'][len(sweep['backward']) - i - 2]
        qf, qb = model.ic.values(f['x']), model.ic.values(b['x'])
        diff = _wrapped_row_difference(model, qf[model.response_rows], qb[model.response_rows])
        top = np.argsort(np.abs(diff))[::-1][:3]
        rows = [(model.describe_row(int(model.response_rows[r]), qf), float(diff[r])) for r in top]
    if not bad and not end_bad:
        classification, verdict = 'continuous', 'accept'
        message = (f'one continuous branch between the nodes at {label_a} and {label_b} of group {gid}: the two sweeps agree '
                   f'with each other (largest change {max([d["angle_change"] for d in sweep["disagreement"]] or [0.0]):.2f} deg, '
                   f'{max([d["bond_change"] for d in sweep["disagreement"]] or [0.0]):.4f} bohr) and with the policy nodes; accepted'
                   + ('; the diagnostics indicate under-resolution: a finer grid is recommended' if under_resolved else ''))
    else:
        degenerate = bool(bad) and all(abs(d['energy_difference']) < policy.branch_degenerate_energy_hartree for d in bad)
        rows_txt = ', '.join(f'{name} ({np.degrees(value):+.1f} deg)' if 'bond' not in name else f'{name} ({value:+.4f})'
                             for name, value in rows)
        if degenerate:
            classification, verdict = 'degenerate_mirror_branches', 'reject'
            message = (f'degenerate mirror branches between the nodes at {label_a} and {label_b} of group {gid}: the two sweeps end '
                       f'in different basins of equal energy (rows {rows_txt}); rejected unless those rows are frozen (relaxation_zone)')
        else:
            classification, verdict = 'branch_flip', 'reject'
            message = (f'branch flip or hysteresis between the nodes at {label_a} and {label_b} of group {gid}: the sweeps disagree '
                       f'{"with each other" if bad else "with a policy node"} beyond the continuity limits (rows {rows_txt}; '
                       f'end points against the policy nodes {sweep["end_forward"][1]:.2f} / {sweep["end_backward"][1]:.2f} deg); '
                       'rejected: refine the grid, freeze the rows of the bifurcating mode (relaxation_zone), or use the rigid policy')
    return {'classification': classification, 'verdict': verdict, 'message': message, 'rows': rows, 'sweep': sweep,
            'gid': gid, 'pair': (float(record_a['offsets'][gid]), float(record_b['offsets'][gid]))}
