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
import warnings

from .veloxchemlib import hartree_in_kcalpermol, bohr_in_angstrom
from .interpolationdriver import InterpolationDriver


class IMTrustRadiusOptimizer:
    """
    Objective and gradient for fitting the confidence radii R_j of a set of
    datapoints to reference energies and gradients.

    For every reference structure s the interpolation is evaluated with the
    trial radii, and the loss

        loss_s = e_x L_e + 1/2 (1 - e_x) L_F,
        L_e = (E_interp - E_ref)^2,
        L_F = beta L_perp + (1 - beta) L_par,

    is formed in kcal/mol and kcal/mol/Angstrom over the Cartesian
    components of the atoms in symmetry_information[3] (nsub of them):
    Delta g = g_interp - g_ref, L_par = (u . Delta g)^2 with u the
    direction of the reference gradient (not divided by nsub) and
    L_perp = |Delta g - (u . Delta g) u|^2 / nsub. For a vanishing
    reference gradient (|g_ref| <= 1e-8) L_F = |Delta g|^2 / nsub.
    fun returns the mean loss over the structures and jac its analytic
    gradient with respect to the radii, from the dw/dR and d(grad w)/dR of
    InterpolationDriver. A structure with an invalid weight state gets the
    penalty loss 1e20 with gradient -1e4.

    The radii are written into the datapoint objects dps themselves
    (dp.confidence_radius) at every evaluation, so after fun or jac the
    datapoints hold the last trial radii; callers that need the original
    radii must restore them.

    :param z_matrix:
        The z-matrix dictionary of the datapoints.
    :param impes_dict:
        The interpolation settings (update_settings of InterpolationDriver).
    :param sym_dict:
        The symmetry_information list of the driver.
    :param dps:
        The InterpolationDatapoint objects whose radii are optimised
        (mutated, see above). dps[0] provides inv_sqrt_masses and
        eq_bond_lengths of the driver geometry.
    :param structure_list:
        The reference molecules.
    :param qm_e:
        The reference energies in Hartree.
    :param qm_g:
        The reference gradients in Hartree/bohr, one (natoms, 3) array per
        structure.
    :param exponent_p_q:
        The Shepard exponents (p, q); p > 10 gives a warning.
    :param e_x:
        The weight of the energy term in [0, 1].
    :param beta:
        The weight of the force component perpendicular to the reference
        gradient (default 0.8).
    :param verbose:
        Unused; kept for backwards compatibility.
    """

    def __init__(self, z_matrix, impes_dict, sym_dict, dps, structure_list, qm_e, qm_g, exponent_p_q, e_x, beta=0.8, verbose=False):
        """
        Initializes the optimizer and its InterpolationDriver (see the class
        docstring for the parameters).
        """

        p, q = exponent_p_q

        if p > 10:
            warnings.warn(
                (
                    f"Large Shepard exponent detected: p = {p}, q = {q}. "
                    "Trust-radius optimization may become numerically unstable because "
                    "the two-part weights can become extremely localized, causing "
                    "underflow of raw weights and division by a near-zero sum of weights. "
                    "Consider using smaller exponents such as p=2,q=1 or p=4,q=2, p=6,q=2"
                ),
                RuntimeWarning,
                stacklevel=2,
            )
        self.z_matrix = z_matrix
        self.impes_dict = impes_dict
        self.sym_dict = sym_dict
        self.dps = dps
        self.structures = structure_list
        self.qm_e = np.asarray(qm_e, dtype=np.float64)
        self.qm_g_flat = np.stack([g.reshape(-1) for g in qm_g], axis=0)  # (S,D)
        self.e_x = float(e_x)
        self.beta = float(beta)
        self.M = len(dps)
        self.S = len(structure_list)
        self.idx = np.asarray([o*3 + i for o in sym_dict[3] for i in range(3)], dtype=int)

        self._init_worker_state(exponent_p_q)

        self._cache_key = None
        self._cache_grad = None
        self._cache_metrics = None
        self._closed = False

    def _init_worker_state(self, exponent_p_q):
        """
        Builds the InterpolationDriver (calc_optim_trust_radius on) over
        the datapoints dps and the per-structure reference data.

        :param exponent_p_q:
            The Shepard exponents (p, q).
        """

        driver = InterpolationDriver(self.z_matrix)
        driver.update_settings(self.impes_dict)
        driver.symmetry_information = self.sym_dict

        driver.impes_coordinate.inv_sqrt_masses = self.dps[0].inv_sqrt_masses
        driver.impes_coordinate.eq_bond_lengths = self.dps[0].eq_bond_lengths

        driver.exponent_p = exponent_p_q[0]
        driver.exponent_q = exponent_p_q[1]
        driver.calc_optim_trust_radius = True

        # the driver evaluates the caller's datapoint objects; their radii
        # are set to the trial values in _eval_structure
        driver.qm_data_points = self.dps

        self._worker_state = {
            "driver": driver,
            "idx": np.asarray(self.idx, dtype=int),
            "nsub": len(self.idx),
            "beta": float(self.beta),
            "e_x": float(self.e_x),
            "conv": float(hartree_in_kcalpermol()),
            "structures": tuple(self.structures),
            "qm_e": np.asarray(self.qm_e, dtype=np.float64),
            "qm_g_flat": np.asarray(self.qm_g_flat, dtype=np.float64),
        }

    def _eval_structure(self, payload):
        """
        Computes the loss of one structure and its gradient with respect to
        the radii.

        :param payload:
            A tuple (structure index, radii).

        :return:
            (loss_s, grad_s, L_e, L_F).
        """
        s_idx, alphas = payload

        W = self._worker_state

        drv = W["driver"]
        idx = W["idx"]
        nsub = W["nsub"]
        beta = W["beta"]
        e_x = W["e_x"]
        conv = W["conv"]
        mol = W["structures"][s_idx]
        E_qm = W["qm_e"][s_idx]
        G_qm_flat = W["qm_g_flat"][s_idx]

        # set the trial radii in the datapoints (the caller's objects)
        for j, dp in enumerate(drv.qm_data_points):
            dp.confidence_radius = alphas[j]

        # Run interpolation for this structure
        drv.compute(mol)
        E_interp = drv.get_energy()                  # Hartree
        G_interp = drv.get_gradient().reshape(-1)    # Hartree/Bohr

        # per-datapoint potentials (M,) and gradients (M, D) of the driver
        P = np.asarray(drv.potentials, dtype=np.float64)
        G = np.asarray(drv.gradients, dtype=np.float64).reshape(len(P), -1)

        w_dα_raw = np.asarray(drv.dw_dalpha_list, dtype=np.float64).reshape(-1)
        w_x_alpha_raw = np.asarray(
            drv.dw_dX_dalpha_list,
            dtype=np.float64,
        ).reshape(len(P), -1)

        def weight_penalty():
            penalty = 1.0e20
            grad_penalty = -1.0e4 * np.ones(self.M, dtype=np.float64)
            return float(penalty), grad_penalty, float(penalty), float(penalty)

        if w_dα_raw.size != len(P) or w_x_alpha_raw.shape[0] != len(P):
            return weight_penalty()

        raw_S = float(
            getattr(
                drv,
                "raw_sum_of_weights",
                getattr(drv, "sum_of_weights", np.nan),
            )
        )

        raw_S_x = np.asarray(
            getattr(
                drv,
                "sum_of_weights_grad",
                np.zeros_like(drv.impes_coordinate.gradient),
            ),
            dtype=np.float64,
        ).reshape(-1)

        used_weight_fallback = bool(getattr(drv, "used_weight_fallback", False))
        used_weight_rescale = bool(getattr(drv, "used_weight_rescale", False))
        weight_eps = float(getattr(drv, "weight_eps", 1.0e-12))
        weight_rescale_factor = float(getattr(drv, "weight_rescale_factor", 1.0))

        if used_weight_rescale:
            if (
                not np.isfinite(weight_rescale_factor)
                or weight_rescale_factor <= 0.0
            ):
                return weight_penalty()

            S = float(getattr(drv, "effective_sum_of_weights", np.nan))
            S_x = np.asarray(
                getattr(drv, "effective_sum_of_weights_grad", raw_S_x),
                dtype=np.float64,
            ).reshape(-1)

            w_dα = w_dα_raw / weight_rescale_factor
            w_x_alpha = w_x_alpha_raw / weight_rescale_factor

        else:
            S = raw_S
            S_x = raw_S_x
            w_dα = w_dα_raw
            w_x_alpha = w_x_alpha_raw

        legacy_invalid_weight_state = (
            not used_weight_rescale
            and ((not np.isfinite(raw_S)) or raw_S <= weight_eps)
        )

        if (
            used_weight_fallback
            or legacy_invalid_weight_state
            or (not np.isfinite(S))
            or S <= 0.0
            or not np.all(np.isfinite(S_x))
            or not np.all(np.isfinite(w_dα))
            or not np.all(np.isfinite(w_x_alpha))
        ):
            return weight_penalty()

        # Residuals in kcal/mol
        dE_res = (E_interp - E_qm) * conv

        # ----- loss (anisotropic force term) -----
        g_sub = (
            (G_interp[idx] - G_qm_flat[idx])
            * conv
            / bohr_in_angstrom()
        )

        h_sub = G_qm_flat[idx] * conv / bohr_in_angstrom()
        nh = np.linalg.norm(h_sub)
        L_iso = (g_sub @ g_sub) / nsub

        if nh > 1.0e-8:
            u = h_sub / nh
            tau = 1.0e-8
            gate = (nh * nh) / (nh * nh + tau * tau)
            gpar = (u @ g_sub) * u
            gper = g_sub - gpar
            L_e = dE_res**2
            L_per = (gper @ gper) / nsub
            L_par = (u @ g_sub)**2
            L_F = (
                gate * (beta * L_per + (1.0 - beta) * L_par)
                + (1.0 - gate) * L_iso
            )
        else:
            L_e = dE_res**2
            L_F = L_iso

        loss_s = 1.0 * e_x * L_e + 0.5 * (1.0 - e_x) * L_F

        G_sub = G[:, idx]
        G_interp_sub = G_interp[idx]
        w_x_alpha_sub = w_x_alpha[:, idx]
        S_x_sub = S_x[idx]

        invS = 1.0 / S
        invS2 = invS * invS

        # dE/dalpha_j = dw_j/dalpha_j * (P_j - E) / S
        dE_dα = (w_dα * (P - E_interp)) * (conv * invS)

        term1 = (
            w_x_alpha_sub * (P - E_interp)[:, None]
        ) * (conv / bohr_in_angstrom() * invS)

        term2 = (
            w_dα[:, None] * (G_sub - G_interp_sub[None, :])
        ) * (conv / bohr_in_angstrom() * invS)

        term3 = (
            (w_dα * (P - E_interp))[:, None]
            * (conv / bohr_in_angstrom())
            * S_x_sub[None, :]
            * invS2
        )

        dG_dα_sub = term1 + term2 - term3

        if nh > 1e-8:
            dL_dg_sub = (2.0 * gate) * ((beta/nsub)*(gper) + (1.0 - beta)*(gpar)) + (2.0 * (1.0 - gate) / nsub) * g_sub
        else:
            dL_dg_sub = 2.0 * g_sub / nsub

        grad_force_part = dG_dα_sub @ dL_dg_sub                     # (M,)
        grad_energy_part = 2.0 * e_x * dE_res * dE_dα                # (M,)
        grad_s = grad_energy_part + 0.5 * (1.0 - e_x) * grad_force_part   # (M,)
        return float(loss_s), grad_s, float(L_e), float(L_F)

    @staticmethod
    def _key(x):
        """
        Returns the cache key of a radius vector (rounded to 1e-14).
        """
        return tuple(np.asarray(x, float).round(14))

    def fun(self, alphas):
        """
        Returns the mean loss over the structures for the radii alphas and
        caches its gradient and the metrics.

        :param alphas:
            The confidence radii of all datapoints in bohr, (M,).
        """
        key = self._key(alphas)
        alphas_arr = np.asarray(alphas, dtype=np.float64).reshape(-1)

        if alphas_arr.size != self.M:
            raise ValueError(
                f"IMTrustRadiusOptimizer.fun expects full alpha size {self.M}, got {alphas_arr.size}"
            )
        if not np.all(np.isfinite(alphas_arr)):
            raise ValueError("IMTrustRadiusOptimizer.fun received non-finite alpha values")

        sum_loss = 0.0
        sum_grad = np.zeros(self.M, dtype=np.float64)
        sum_energy = 0.0
        sum_force = 0.0
        for i in range(self.S):
            loss_s, grad_s, energy_s, force_s = self._eval_structure((i, alphas_arr))
            sum_loss += loss_s
            sum_grad += grad_s
            sum_energy += energy_s
            sum_force += force_s

        mean_loss = float(sum_loss / self.S)
        mean_energy = float(sum_energy / self.S)
        mean_force = float(sum_force / self.S)
        self._cache_key = key
        self._cache_grad = sum_grad / self.S
        self._cache_metrics = {
            'loss_total': mean_loss,
            'loss_energy': mean_energy,
            'loss_force': mean_force,
        }

        return mean_loss

    def jac(self, alphas):
        """
        Returns the gradient of fun with respect to the radii (from the
        cache if fun was evaluated at the same radii).

        :param alphas:
            The confidence radii of all datapoints in bohr, (M,).
        """
        key = self._key(alphas)
        if key == self._cache_key and self._cache_grad is not None:
            return self._cache_grad
        # Fallback: compute once (will also cache)
        _ = self.fun(alphas)
        return self._cache_grad

    def _validate_trainable_idx(self, trainable_idx):
        """
        Returns the indices of the trainable radii as an array (checks range
        and duplicates).
        """
        idx = np.asarray(trainable_idx, dtype=np.int64).reshape(-1)

        if idx.size == 0:
            return idx

        if np.any(idx < 0) or np.any(idx >= self.M):
            raise ValueError(
                f"trainable_idx out of bounds for M={self.M}: {idx.tolist()}"
            )

        if np.unique(idx).size != idx.size:
            raise ValueError("trainable_idx contains duplicates")

        return idx

    def _validate_full_alpha(self, alpha_full):
        """
        Returns the full radius vector as an array (checks size and
        finiteness).
        """
        alpha = np.asarray(alpha_full, dtype=np.float64).reshape(-1)

        if alpha.size != self.M:
            raise ValueError(
                f"Expected full alpha size {self.M}, got {alpha.size}"
            )
        if not np.all(np.isfinite(alpha)):
            raise ValueError("alpha contains non-finite values")

        return alpha

    def pack_reduced_alphas(self, alpha_full, trainable_idx):
        """
        Returns the trainable radii alpha_full[trainable_idx].
        """

        alpha = self._validate_full_alpha(alpha_full)
        idx = self._validate_trainable_idx(trainable_idx)
        return alpha[idx].copy()

    def expand_reduced_alphas(self, x_var, trainable_idx, alpha_full_fixed):
        """
        Returns the full radius vector: alpha_full_fixed with the trainable
        entries replaced by x_var.
        """

        idx = self._validate_trainable_idx(trainable_idx)
        alpha_fixed = self._validate_full_alpha(alpha_full_fixed)

        x = np.asarray(x_var, dtype=np.float64).reshape(-1)
        if x.size != idx.size:
            raise ValueError(
                f"Reduced alpha size mismatch: expected {idx.size}, got {x.size}"
            )
        if not np.all(np.isfinite(x)):
            raise ValueError("x_var contains non-finite values")

        alpha_full = alpha_fixed.copy()
        alpha_full[idx] = x
        return alpha_full

    def project_full_gradient_to_reduced(self, grad_full, trainable_idx):
        """
        Returns the gradient entries of the trainable radii.
        """

        g = np.asarray(grad_full, dtype=np.float64).reshape(-1)
        if g.size != self.M:
            raise ValueError(
                f"Expected full gradient size {self.M}, got {g.size}"
            )
        idx = self._validate_trainable_idx(trainable_idx)
        return g[idx].copy()

    def fun_reduced(self, x_var, trainable_idx, alpha_full_fixed):
        """
        fun for the trainable radii x_var, the others fixed at
        alpha_full_fixed.
        """

        alpha_full = self.expand_reduced_alphas(x_var, trainable_idx, alpha_full_fixed)
        return self.fun(alpha_full)

    def jac_reduced(self, x_var, trainable_idx, alpha_full_fixed):
        """
        jac for the trainable radii x_var, the others fixed at
        alpha_full_fixed.
        """

        alpha_full = self.expand_reduced_alphas(x_var, trainable_idx, alpha_full_fixed)
        g_full = self.jac(alpha_full)
        return self.project_full_gradient_to_reduced(g_full, trainable_idx)

    def get_metrics_reduced(self, x_var, trainable_idx, alpha_full_fixed, evaluate_if_missing=True):
        """
        get_metrics for the trainable radii x_var, the others fixed at
        alpha_full_fixed.
        """
        alpha_full = self.expand_reduced_alphas(x_var, trainable_idx, alpha_full_fixed)
        return self.get_metrics(alpha_full, evaluate_if_missing=evaluate_if_missing)

    def get_metrics(self, alphas=None, evaluate_if_missing=True):
        """
        Returns the mean loss and its energy and force parts
        ({'loss_total', 'loss_energy', 'loss_force'}) of the radii alphas,
        or of the last evaluation if alphas is None.

        :param alphas:
            The confidence radii in bohr, or None.
        :param evaluate_if_missing:
            Evaluate fun if the radii are not the cached ones (otherwise
            return None).
        """

        if alphas is None:
            return dict(self._cache_metrics) if self._cache_metrics is not None else None

        key = self._key(alphas)
        if key == self._cache_key and self._cache_metrics is not None:
            return dict(self._cache_metrics)

        if evaluate_if_missing:
            _ = self.fun(alphas)
            return dict(self._cache_metrics) if self._cache_metrics is not None else None

        return None

    def close(self):
        """
        Releases the driver state.
        """
        if not self._closed:
            self._worker_state = None
            self._closed = True

    def __del__(self):
        """
        Closes the optimizer.
        """
        try:
            self.close()
        except Exception:
            pass


class IMTrustRadiusPairTerms:
    """
    Radius-independent terms of the Shepard interpolation at a fixed list of
    reference structures. The interpolated energy and gradient, the loss of
    IMTrustRadiusOptimizer and its radius gradient follow from them in closed
    form for any confidence radii.

    For reference s and datapoint j the terms are those of
    InterpolationDriver.compute:
    - the distance d_sj and the aligned displacement v_sj of
      cartesian_distance;
    - the multiplicative tc gate metric A_sj and its gradient (zero without tc
      weights);
    - the Taylor value T_sj and its gradient grad T_sj.

    None of them depends on a radius. The radius R_j enters only through
    u = d^2 / R_j^2 in the raw weight

        w = exp(-min(lambda A, clip)) / (u^p + u^q),

    with the conventions of the driver:
    - w = 0 for u > 1e6 or w <= 1e-300;
    - the gate gradient in grad w is zero where lambda A is clipped;
    - d(grad w)/dR uses the unclipped lambda (K52), as
      InterpolationDriver.cartesian_distance does.

    The vectors are stored only for the datapoints free_idx, whose radii may
    change. Every other datapoint enters through per-reference sums at its
    fixed radius (sum w, sum w T, sum grad w, sum (w grad T + T grad w)), which
    is exact. All vectors are flattened Cartesian vectors of length 3 natoms.

    :param n_ref:
        The number of reference structures.
    :param natoms:
        The number of atoms.
    :param radii:
        The radii of all datapoints in bohr. The entries of the datapoints not
        in free_idx stay fixed.
    :param free_idx:
        The indices of the datapoints whose radii may change.
    :param exponent_p_q:
        The Shepard exponents (p, q).
    :param gate_lambda:
        tc_imp_gate_lambda of the driver.
    :param gate_clip:
        tc_imp_gate_exp_clip of the driver.
    """

    def __init__(self, n_ref, natoms, radii, free_idx, exponent_p_q, gate_lambda, gate_clip):
        """
        Allocates the terms (see the class docstring for the parameters).
        """

        self.n_ref = int(n_ref)
        self.natoms = int(natoms)
        self.radii = np.asarray(radii, dtype=np.float64).reshape(-1).copy()
        self.n_dp = self.radii.size
        self.free_idx = np.asarray(free_idx, dtype=np.int64).reshape(-1)
        self.fixed_mask = np.ones(self.n_dp, dtype=bool)
        self.fixed_mask[self.free_idx] = False
        self.p = float(exponent_p_q[0])
        self.q = float(exponent_p_q[1])
        self.gate_lambda = float(gate_lambda)
        self.gate_clip = float(gate_clip)

        nfree = self.free_idx.size
        ndof = 3 * self.natoms

        self.dist = np.zeros((self.n_ref, nfree))
        self.A = np.zeros((self.n_ref, nfree))
        self.T = np.zeros((self.n_ref, nfree))
        self.V = np.zeros((self.n_ref, nfree, ndof))
        self.GA = np.zeros((self.n_ref, nfree, ndof))
        self.GT = np.zeros((self.n_ref, nfree, ndof))

        # sums over the fixed datapoints at their radii
        self.fixed_w = np.zeros(self.n_ref)
        self.fixed_w_max = np.zeros(self.n_ref)
        self.fixed_wT = np.zeros(self.n_ref)
        self.fixed_grad_w = np.zeros((self.n_ref, ndof))
        self.fixed_grad_wT = np.zeros((self.n_ref, ndof))

    @staticmethod
    def memory_estimate(n_ref, natoms, n_free):
        """
        Returns the memory of the stored terms in bytes.
        """

        ndof = 3 * int(natoms)
        return 8 * int(n_ref) * (int(n_free) * (3 + 3 * ndof) + 3 + 2 * ndof)

    def weight_terms(self, dist, A, R, radius_derivatives=False):
        """
        Returns, elementwise, the raw weight w and the coefficients cv, cg of
        grad w = cv v + cg grad A. With radius_derivatives it also returns
        dw/dR and the coefficients dv, dg of d(grad w)/dR = dv v + dg grad A.

        :param dist:
            The distances d (bohr).
        :param A:
            The gate metrics.
        :param R:
            The radii (bohr), broadcast against dist.
        :param radius_derivatives:
            Also return the radius derivatives.
        """

        p, q = self.p, self.q
        lam, clip = self.gate_lambda, self.gate_clip

        u = dist**2 / R**2
        far = u > 1.0e6
        u = np.where(far, 1.0, u)
        up, uq = u**p, u**q
        up1, uq1 = u**(p - 1.0), u**(q - 1.0)
        den = np.maximum(up + uq, 1.0e-300)

        x = lam * A
        clipped = x >= clip
        gate = np.exp(-np.where(clipped, clip, x))

        w_base = 1.0 / den
        w = w_base * gate
        zero = far | (w <= 1.0e-300)
        w = np.where(zero, 0.0, w)
        cv = np.where(zero, 0.0, gate * (-2.0 * (p * up1 + q * uq1) / (R**2 * den**2)))
        cg = np.where(zero | clipped, 0.0, -lam * gate * w_base)

        if not radius_derivatives:
            return w, cv, cg

        dw_base = np.where(far, 0.0, (2.0 / R) * (p * up + q * uq) / den**2)
        t1 = (p * p * up1 + q * q * uq1) / den**2
        t2 = 2.0 * ((p * up1 + q * uq1) * (p * up + q * uq)) / den**3
        dv = np.where(far, 0.0, gate * (4.0 / R**3) * (t1 - t2))
        dg = -lam * gate * dw_base

        return w, cv, cg, gate * dw_base, dv, dg

    def interpolate_reference(self, dist, V, A, GA, T, GT, radii):
        """
        Returns the interpolated energy and gradient (flat) of one reference
        from the terms of all its datapoints, (M,) and (M, 3 natoms) arrays,
        at the radii (M,). Without a positive weight the nearest datapoint is
        used, as in the driver.
        """

        w, cv, cg = self.weight_terms(dist, A, radii)

        if not np.any(w > 0.0):
            nearest = int(np.argmin(dist))
            return float(T[nearest]), np.array(GT[nearest], dtype=np.float64)

        scale = np.max(w)
        w, cv, cg = w / scale, cv / scale, cg / scale
        gw = cv[:, None] * V + cg[:, None] * GA

        S = np.sum(w)
        E = np.dot(w, T) / S
        G = (w @ GT + T @ gw) / S - E * np.sum(gw, axis=0) / S

        return float(E), G

    def set_reference(self, s, dist, V, A, GA, T, GT):
        """
        Stores the terms of reference s: dist, A, T (M,) and V, GA, GT
        (M, 3 natoms), for all datapoints.
        """

        free = self.free_idx
        self.dist[s] = dist[free]
        self.A[s] = A[free]
        self.T[s] = T[free]
        self.V[s] = V[free]
        self.GA[s] = GA[free]
        self.GT[s] = GT[free]

        fixed = self.fixed_mask
        if np.any(fixed):
            w, cv, cg = self.weight_terms(dist[fixed], A[fixed], self.radii[fixed])
            gw = cv[:, None] * V[fixed] + cg[:, None] * GA[fixed]
            self.fixed_w[s] = np.sum(w)
            self.fixed_w_max[s] = np.max(w)
            self.fixed_wT[s] = np.dot(w, T[fixed])
            self.fixed_grad_w[s] = np.sum(gw, axis=0)
            self.fixed_grad_wT[s] = w @ GT[fixed] + T[fixed] @ gw

    def _check_radii(self, radii):
        """
        Returns the radii as an array; the radii of the fixed datapoints must
        equal the stored ones.
        """

        R = np.asarray(radii, dtype=np.float64).reshape(-1)
        if R.size != self.n_dp:
            raise ValueError(f"IMTrustRadiusPairTerms: expected {self.n_dp} radii, got {R.size}")
        if not np.array_equal(R[self.fixed_mask], self.radii[self.fixed_mask]):
            raise ValueError("IMTrustRadiusPairTerms: the radii of the fixed datapoints changed")
        return R

    def _evaluate(self, rows, R, train_pos=None):
        """
        Returns E (n,), G (n, 3 natoms), the sum of the weights, its gradient
        and the rows without a positive weight, all with the weights divided
        by a per-row scale. With train_pos (positions in free_idx) it also
        returns the scaled dw/dR (n, t) and d(grad w)/dR (n, t, 3 natoms).
        """

        Rf = R[self.free_idx]
        dist, A, T = self.dist[rows], self.A[rows], self.T[rows]
        V, GA, GT = self.V[rows], self.GA[rows], self.GT[rows]

        if train_pos is None:
            w, cv, cg = self.weight_terms(dist, A, Rf)
        else:
            w, cv, cg, dw, dv, dg = self.weight_terms(dist, A, Rf, radius_derivatives=True)

        w_max = np.max(w, axis=1) if w.shape[1] > 0 else np.zeros(w.shape[0])
        scale = np.maximum(w_max, self.fixed_w_max[rows])
        no_weight = ~(scale > 0.0)
        scale = np.where(no_weight, 1.0, scale)
        inv = 1.0 / scale

        w, cv, cg = w * inv[:, None], cv * inv[:, None], cg * inv[:, None]
        gw = cv[..., None] * V + cg[..., None] * GA

        S = np.sum(w, axis=1) + self.fixed_w[rows] * inv
        S = np.where(no_weight, 1.0, S)
        E = (np.sum(w * T, axis=1) + self.fixed_wT[rows] * inv) / S
        gS = np.sum(gw, axis=1) + self.fixed_grad_w[rows] * inv[:, None]
        G = ((np.einsum('nf,nfd->nd', w, GT) + np.einsum('nf,nfd->nd', T, gw)
              + self.fixed_grad_wT[rows] * inv[:, None]) / S[:, None]
             - E[:, None] * gS / S[:, None])

        if train_pos is None:
            return E, G, S, gS, no_weight

        t = train_pos
        dw_t = dw[:, t] * inv[:, None]
        dgw_t = (dv[:, t] * inv[:, None])[..., None] * V[:, t] + (dg[:, t] * inv[:, None])[..., None] * GA[:, t]

        return E, G, S, gS, no_weight, dw_t, dgw_t

    def energy_gradient(self, rows, radii, chunk=512):
        """
        Returns the interpolated energies (n,) and flat gradients
        (n, 3 natoms) of the references rows at the radii.
        """

        R = self._check_radii(radii)
        rows = np.asarray(rows, dtype=np.int64).reshape(-1)
        E = np.zeros(rows.size)
        G = np.zeros((rows.size, 3 * self.natoms))

        for start in range(0, rows.size, chunk):
            sub = rows[start:start + chunk]
            E_c, G_c, _, _, no_weight = self._evaluate(sub, R)
            if np.any(no_weight):
                raise ValueError("IMTrustRadiusPairTerms: a reference has no positive weight")
            E[start:start + chunk] = E_c
            G[start:start + chunk] = G_c

        return E, G

    def loss_and_gradient(self, rows, radii, qm_e, qm_g_flat, idx, e_x, beta, train_pos, chunk=512):
        """
        Returns the loss terms of IMTrustRadiusOptimizer._eval_structure for
        the references rows at the radii: loss_s, L_e and L_F (n,) and the
        gradient of loss_s with respect to the radii train_pos (positions in
        free_idx), (n, t). A reference without a positive weight gets the
        penalty of the optimizer (loss 1e20, gradient -1e4).

        :param qm_e:
            The reference energies of rows (Hartree), (n,).
        :param qm_g_flat:
            The reference gradients of rows (Hartree/bohr), (n, 3 natoms).
        :param idx:
            The Cartesian components in the force loss.
        """

        R = self._check_radii(radii)
        rows = np.asarray(rows, dtype=np.int64).reshape(-1)
        train_pos = np.asarray(train_pos, dtype=np.int64).reshape(-1)
        qm_e = np.asarray(qm_e, dtype=np.float64).reshape(-1)
        qm_g_flat = np.asarray(qm_g_flat, dtype=np.float64)

        conv = float(hartree_in_kcalpermol())
        conv_g = conv / bohr_in_angstrom()
        nsub = len(idx)

        loss = np.zeros(rows.size)
        L_e = np.zeros(rows.size)
        L_F = np.zeros(rows.size)
        grad = np.zeros((rows.size, train_pos.size))

        for start in range(0, rows.size, chunk):
            sl = slice(start, start + chunk)
            sub = rows[sl]
            E, G, S, gS, no_weight, dw_t, dgw_t = self._evaluate(sub, R, train_pos)

            dE = (E - qm_e[sl]) * conv
            g = (G[:, idx] - qm_g_flat[sl][:, idx]) * conv_g
            h = qm_g_flat[sl][:, idx] * conv_g
            nh = np.linalg.norm(h, axis=1)

            L_iso = np.sum(g * g, axis=1) / nsub
            big = nh > 1.0e-8
            u = np.where(big[:, None], h / np.where(big, nh, 1.0)[:, None], 0.0)
            gate = (nh * nh) / (nh * nh + 1.0e-8 * 1.0e-8)
            par = np.sum(u * g, axis=1)
            gpar = par[:, None] * u
            gper = g - gpar
            L_per = np.sum(gper * gper, axis=1) / nsub
            L_F_c = np.where(big, gate * (beta * L_per + (1.0 - beta) * par**2) + (1.0 - gate) * L_iso, L_iso)
            L_e_c = dE**2
            loss_c = e_x * L_e_c + 0.5 * (1.0 - e_x) * L_F_c

            dL_dg = np.where(big[:, None],
                             (2.0 * gate[:, None]) * ((beta / nsub) * gper + (1.0 - beta) * gpar)
                             + (2.0 * (1.0 - gate[:, None]) / nsub) * g,
                             2.0 * g / nsub)

            Pm = self.T[sub][:, train_pos] - E[:, None]
            dE_dR = dw_t * Pm * (conv / S[:, None])
            GT_t = self.GT[sub][:, train_pos]
            dG_dR = (dgw_t * Pm[..., None]
                     + dw_t[..., None] * (GT_t - G[:, None, :])
                     - (dw_t * Pm)[..., None] * gS[:, None, :] / S[:, None, None]) * (conv_g / S[:, None, None])
            grad_c = (2.0 * e_x * dE[:, None] * dE_dR
                      + 0.5 * (1.0 - e_x) * np.einsum('ntd,nd->nt', dG_dR[..., idx], dL_dg))

            invalid = no_weight | ~np.isfinite(loss_c) | ~np.all(np.isfinite(grad_c), axis=1)
            loss_c = np.where(invalid, 1.0e20, loss_c)
            L_e_c = np.where(invalid, 1.0e20, L_e_c)
            L_F_c = np.where(invalid, 1.0e20, L_F_c)
            grad_c = np.where(invalid[:, None], -1.0e4, grad_c)

            loss[sl], L_e[sl], L_F[sl], grad[sl] = loss_c, L_e_c, L_F_c, grad_c

        return loss, L_e, L_F, grad


class IMTrustRadiusCachedOptimizer(IMTrustRadiusOptimizer):
    """
    IMTrustRadiusOptimizer evaluated on IMTrustRadiusPairTerms: the same loss
    and gradient, to round-off, in closed form instead of one
    InterpolationDriver.compute per structure and evaluation. The datapoint
    objects are not modified.

    The gradient is computed for trainable_idx only; its other entries are
    zero, which is what the reduced methods (fun_reduced, jac_reduced) use.
    The radii of the datapoints outside pair_terms.free_idx must stay at
    pair_terms.radii.

    :param pair_terms:
        The IMTrustRadiusPairTerms of the reference structures.
    :param rows:
        The indices of the structures of the objective in pair_terms.
    :param qm_e:
        The reference energies of rows in Hartree.
    :param qm_g:
        The reference gradients of rows in Hartree/bohr, one (natoms, 3)
        array per structure.
    :param sym_dict:
        The symmetry_information list of the driver ([3]: the atoms in the
        force loss).
    :param trainable_idx:
        The indices of the trainable radii (a subset of pair_terms.free_idx).
    :param e_x:
        The weight of the energy term in [0, 1].
    :param beta:
        The weight of the force component perpendicular to the reference
        gradient (default 0.8).
    """

    def __init__(self, pair_terms, rows, qm_e, qm_g, sym_dict, trainable_idx, e_x, beta=0.8):
        """
        Initializes the optimizer (see the class docstring for the
        parameters).
        """

        self.pair_terms = pair_terms
        self.rows = np.asarray(rows, dtype=np.int64).reshape(-1)
        self.qm_e = np.asarray(qm_e, dtype=np.float64).reshape(-1)
        self.qm_g_flat = np.stack([np.asarray(g, dtype=np.float64).reshape(-1) for g in qm_g], axis=0)
        self.e_x = float(e_x)
        self.beta = float(beta)
        self.M = pair_terms.n_dp
        self.S = self.rows.size
        self.idx = np.asarray([o*3 + i for o in sym_dict[3] for i in range(3)], dtype=int)

        self.trainable_idx = self._validate_trainable_idx(trainable_idx)
        free_pos = {int(j): k for k, j in enumerate(pair_terms.free_idx)}
        missing = [int(j) for j in self.trainable_idx if int(j) not in free_pos]
        if missing:
            raise ValueError(f"IMTrustRadiusCachedOptimizer: trainable datapoints {missing} are not free in pair_terms")
        self._train_pos = np.asarray([free_pos[int(j)] for j in self.trainable_idx], dtype=np.int64)

        self._worker_state = None
        self._cache_key = None
        self._cache_grad = None
        self._cache_metrics = None
        self._closed = False

    def fun(self, alphas):
        """
        Returns the mean loss over the structures for the radii alphas and
        caches its gradient and the metrics.

        :param alphas:
            The confidence radii of all datapoints in bohr, (M,).
        """
        key = self._key(alphas)
        alphas_arr = np.asarray(alphas, dtype=np.float64).reshape(-1)

        if alphas_arr.size != self.M:
            raise ValueError(
                f"IMTrustRadiusCachedOptimizer.fun expects full alpha size {self.M}, got {alphas_arr.size}"
            )
        if not np.all(np.isfinite(alphas_arr)):
            raise ValueError("IMTrustRadiusCachedOptimizer.fun received non-finite alpha values")

        loss, L_e, L_F, grad = self.pair_terms.loss_and_gradient(
            self.rows, alphas_arr, self.qm_e, self.qm_g_flat, self.idx,
            self.e_x, self.beta, self._train_pos)

        mean_loss = float(np.sum(loss) / self.S)
        grad_full = np.zeros(self.M, dtype=np.float64)
        grad_full[self.trainable_idx] = np.sum(grad, axis=0) / self.S

        self._cache_key = key
        self._cache_grad = grad_full
        self._cache_metrics = {
            'loss_total': mean_loss,
            'loss_energy': float(np.sum(L_e) / self.S),
            'loss_force': float(np.sum(L_F) / self.S),
        }

        return mean_loss
