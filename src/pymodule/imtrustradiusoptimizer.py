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


from mpi4py import MPI
import numpy as np
import scipy
import copy
import warnings

from contextlib import redirect_stderr
from io import StringIO

from .veloxchemlib import mpi_master
from.veloxchemlib import hartree_in_kcalpermol, bohr_in_angstrom
from .interpolationdriver import InterpolationDriver


with redirect_stderr(StringIO()) as fg_err:
    import geometric

try:
    import openmm as mm
    import openmm.app as app
    import openmm.unit as unit
except ImportError:
    pass



class IMTrustRadiusOptimizer:
    def __init__(self, z_matrix, impes_dict, sym_dict, dps, structure_list, qm_e, qm_g, exponent_p_q, e_x, beta=0.8, verbose=False):

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
        self.z_matrix       = z_matrix
        self.impes_dict     = impes_dict
        self.sym_dict       = sym_dict
        self.dps            = dps
        self.structures     = structure_list
        self.qm_e           = np.asarray(qm_e, dtype=np.float64)
        self.qm_g_flat      = np.stack([g.reshape(-1) for g in qm_g], axis=0)  # (S,D)
        self.e_x            = float(e_x)
        self.beta           = float(beta)
        self.M              = len(dps)
        self.S              = len(structure_list)
        self.idx            = np.asarray([o*3 + i for o in sym_dict[3] for i in range(3)], dtype=int)
        self.verbose        = bool(verbose)

        self._init_worker_state(exponent_p_q)

        self._cache_key = None
        self._cache_grad = None
        self._cache_metrics = None
        self._closed = False

    def _init_worker_state(self, exponent_p_q):

        driver = InterpolationDriver(self.z_matrix)
        driver.update_settings(self.impes_dict)
        driver.symmetry_information = self.sym_dict

        driver.impes_coordinate.inv_sqrt_masses = self.dps[0].inv_sqrt_masses
        driver.impes_coordinate.eq_bond_lengths = self.dps[0].eq_bond_lengths

        driver.exponent_p = exponent_p_q[0]
        driver.exponent_q = exponent_p_q[1]
        driver.print = False
        driver.calc_optim_trust_radius = True

        # Keep this if datapoints are mutated during evaluation
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
        """Compute loss and grad contrib for one structure index."""
        s_idx, alphas = payload

        W = self._worker_state

        drv    = W["driver"]
        idx    = W["idx"]; nsub = W["nsub"]
        beta   = W["beta"]; e_x = W["e_x"]; conv = W["conv"]
        mol = W["structures"][s_idx]
        E_qm = W["qm_e"][s_idx]
        G_qm_flat = W["qm_g_flat"][s_idx]

        # Update alphas in this worker's private datapoints
        for j, dp in enumerate(drv.qm_data_points):
            dp.confidence_radius = alphas[j]

        # Run interpolation for this structure
        drv.compute(mol)
        E_interp = drv.get_energy()                  # Hartree
        G_interp = drv.get_gradient().reshape(-1)    # Hartree/Bohr

        # Vectorized getters (M, MxD, etc.)
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

        grad_force_part  = dG_dα_sub @ dL_dg_sub                     # (M,)
        grad_energy_part = 2.0 * e_x * dE_res * dE_dα                # (M,)
        grad_s = grad_energy_part + 0.5 * (1.0 - e_x) * grad_force_part   # (M,)
        return float(loss_s), grad_s, float(L_e), float(L_F)

    @staticmethod
    def _key(x):
        return tuple(np.asarray(x, float).round(14))

    def fun(self, alphas):
        key = self._key(alphas)
        alphas_arr = np.asarray(alphas, dtype=np.float64).reshape(-1)

        if alphas_arr.size != self.M:
            raise ValueError(
                f"AlphaOptimizer.fun expects full alpha size {self.M}, got {alphas_arr.size}"
            )
        if not np.all(np.isfinite(alphas_arr)):
            raise ValueError("AlphaOptimizer.fun received non-finite alpha values")

        futs = (self._eval_structure((i, alphas_arr)) for i in range(self.S))
 

        sum_loss = 0.0
        sum_grad = np.zeros(self.M, dtype=np.float64)
        sum_energy = 0.0
        sum_force = 0.0
        counter = 0
        for loss_s, grad_s, energy_s, force_s in futs:
            sum_loss += loss_s
            sum_grad += grad_s
            sum_energy += energy_s
            sum_force += force_s
            counter += 1

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
        key = self._key(alphas)
        if key == self._cache_key and self._cache_grad is not None:
            return self._cache_grad
        # Fallback: compute once (will also cache)
        _ = self.fun(alphas)
        return self._cache_grad

    def _validate_trainable_idx(self, trainable_idx):
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
        alpha = np.asarray(alpha_full, dtype=np.float64).reshape(-1)

        if alpha.size != self.M:
            raise ValueError(
                f"Expected full alpha size {self.M}, got {alpha.size}"
            )
        if not np.all(np.isfinite(alpha)):
            raise ValueError("alpha contains non-finite values")

        return alpha


    def pack_reduced_alphas(self, alpha_full, trainable_idx):

        alpha = self._validate_full_alpha(alpha_full)
        idx = self._validate_trainable_idx(trainable_idx)
        return alpha[idx].copy()


    def expand_reduced_alphas(self, x_var, trainable_idx, alpha_full_fixed):

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

        g = np.asarray(grad_full, dtype=np.float64).reshape(-1)
        if g.size != self.M:
            raise ValueError(
                f"Expected full gradient size {self.M}, got {g.size}"
            )
        idx = self._validate_trainable_idx(trainable_idx)
        return g[idx].copy()


    def fun_reduced(self, x_var, trainable_idx, alpha_full_fixed):

        alpha_full = self.expand_reduced_alphas(x_var, trainable_idx, alpha_full_fixed)
        return self.fun(alpha_full)


    def jac_reduced(self, x_var, trainable_idx, alpha_full_fixed):

        alpha_full = self.expand_reduced_alphas(x_var, trainable_idx, alpha_full_fixed)
        g_full = self.jac(alpha_full)
        return self.project_full_gradient_to_reduced(g_full, trainable_idx)


    def get_metrics_reduced(self, x_var, trainable_idx, alpha_full_fixed, evaluate_if_missing=True):
        alpha_full = self.expand_reduced_alphas(x_var, trainable_idx, alpha_full_fixed)
        return self.get_metrics(alpha_full, evaluate_if_missing=evaluate_if_missing)


    def get_metrics(self, alphas=None, evaluate_if_missing=True):

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
        if not self._closed:
            # self._pool.shutdown(wait=True)
            self._worker_state = None
            self._closed = True

    def __del__(self):
        try:
            self.close()
        except Exception:
            pass
