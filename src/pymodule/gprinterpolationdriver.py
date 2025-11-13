#                           VELOXCHEM 1.0-RC2
#         ----------------------------------------------------
#                     An Electronic Structure Code
#
#  Copyright © 2018-2021 by VeloxChem developers. All rights reserved.
#  Contact: https://veloxchem.org/contact
#
#  SPDX-License-Identifier: LGPL-3.0-or-later
#
#  This file is part of VeloxChem.
#
#  VeloxChem is free software: you can redistribute it and/or modify it under
#  the terms of the GNU Lesser General Public License as published by the Free
#  Software Foundation, either version 3 of the License, or (at your option)
#  any later version.
#
#  VeloxChem is distributed in the hope that it will be useful, but WITHOUT
#  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
#  License for more details.
#
#  You should have received a copy of the GNU Lesser General Public License
#  along with VeloxChem. If not, see <https://www.gnu.org/licenses/>.

from tkinter import E
from mpi4py import MPI
import multiprocessing as mp
from collections import Counter
import os
import numpy as np
import math
import random
import torch
import gpytorch
from gpytorch.constraints import Interval
from gpytorch.priors import LogNormalPrior
from gpytorch.constraints import GreaterThan
from scipy.optimize import linear_sum_assignment
from scipy.linalg import cho_factor, cho_solve
from collections import deque
import sys
from .profiler import Profiler
import h5py
from contextlib import redirect_stderr
from io import StringIO
from .interpolationdatapoint import InterpolationDatapoint
from .atommapper import AtomMapper
from .molecule import Molecule
from .outputstream import OutputStream
from .veloxchemlib import mpi_master
from .veloxchemlib import hartree_in_kcalpermol, bohr_in_angstrom
from .errorhandler import assert_msg_critical
from .inputparser import (parse_input, print_keywords, print_attributes)

with redirect_stderr(StringIO()) as fg_err:
    import geometric


class ExactDeltaE(gpytorch.models.ExactGP):
    def __init__(self, x_train, y_train, likl, base_kernel=None, nu=2.5):
        super().__init__(x_train, y_train, likl)
        self.mean_module = gpytorch.means.ZeroMean()
        if base_kernel is None:
            # Fallback: old behavior (ARD per feature)
            base_kernel = gpytorch.kernels.MaternKernel(nu=nu, ard_num_dims=x_train.shape[-1])
        self.covar_module = gpytorch.kernels.ScaleKernel(base_kernel)

    def forward(self, x):
        return gpytorch.distributions.MultivariateNormal(
            self.mean_module(x), self.covar_module(x)
        )

class GPRInterpolationDriver:

    def __init__(self, x0, groups, y0, *, standardize=True):
        
        dtype = torch.double
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")


        self.standardize = standardize
        self.dtype = dtype
        self.device = device
        self.verbose_diagnostics = False
        self._kernel_diag_cache = {}
        self._diag_cache_version = 0
        self.scale_mode = "robust"   # or "phys" once you wire the feature types
        self.x_floor = 0.01          # minimum per-feature scale (tune 0.02–0.10)
        self._safe_dists = []
        self._sigma_thresh_dyn = 0.0
        self._ensure_sim_state()
        self._ensure_dist_state()
        self.cause = 'none'
        

        self.X_raw = self._ensure_2d(x0)
        self.y_raw = self._ensure_1d(y0)
        self._keep_cols = None 

        self.X_raw = torch.as_tensor(x0, dtype=dtype, device=device)
        self.y_raw = torch.as_tensor(y0, dtype=dtype, device=device)

        self._fit_data()

        X, y = self._transform(self.X_raw, self.y_raw)

        # --- Build grouped Matérn kernels (one ℓ per group via active_dims) ---
        def _active_dims(arr_like):
            # GPyTorch wants a Python list/tuple of ints
            if arr_like is None:
                return []
            a = arr_like.tolist() if hasattr(arr_like, "tolist") else list(arr_like)
            return sorted(int(i) for i in a)

        bond_idx    = _active_dims(groups.get('bond', []))
        angle_idx   = _active_dims(groups.get('angle', []))
        torsion_idx = _active_dims(groups.get('torsion', []))

        self.group_names  = ["bond", "angle", "torsion"]
        self.group_dims   = {
            "bond":    _active_dims(groups.get("bond", [])),      # list[int]
            "angle":   _active_dims(groups.get("angle", [])),
            "torsion": _active_dims(groups.get("torsion", [])),
        }

        self.num_groups = sum(len(v) > 0 for v in self.group_dims.values())

        ls_prior = LogNormalPrior(0.0, 0.5)  # centered ~1.0 after standardization
        ls_constraint = Interval(
            torch.tensor(0.05, dtype=self.dtype, device=self.device),
            torch.tensor(10.0, dtype=self.dtype, device=self.device),
        )

        kernels = []
        if len(bond_idx):
            kernels.append(gpytorch.kernels.MaternKernel(
                nu=1.5, active_dims=bond_idx,
                lengthscale_prior=ls_prior, lengthscale_constraint=ls_constraint))
        if len(angle_idx):
            kernels.append(gpytorch.kernels.MaternKernel(
                nu=1.5, active_dims=angle_idx,
                lengthscale_prior=ls_prior, lengthscale_constraint=ls_constraint))
        if len(torsion_idx):
            kernels.append(gpytorch.kernels.MaternKernel(
                nu=1.5, active_dims=torsion_idx,   # torsion contains both sin & cos cols
                lengthscale_prior=ls_prior, lengthscale_constraint=ls_constraint))

        if not kernels:
            raise ValueError("No non-empty groups in 'groups' dict.")
        
        self.group_kernels = kernels[:]

        base_kernel = kernels[0]
        for k in kernels[1:]:
            base_kernel = base_kernel + k

        self.likelihood = gpytorch.likelihoods.GaussianLikelihood(
            noise_prior=gpytorch.priors.LogNormalPrior(-4.6, 0.5),  # ~0.01
            noise_constraint=gpytorch.constraints.Interval(1e-4, 5e-2),
        ).to(self.device).double()
        scaled_kernel = gpytorch.kernels.ScaleKernel(
            base_kernel,
            outputscale_prior=gpytorch.priors.LogNormalPrior(0.0, 0.5),
            outputscale_constraint=gpytorch.constraints.Interval(1e-4, 1e2),
        ) 
        self.model = ExactDeltaE(X, y, self.likelihood, base_kernel=base_kernel, nu=1.5).to(self.device).double()
        self.model.covar_module = scaled_kernel
        self.model = self.model.to(self.device).double()

        self.mll = gpytorch.mlls.ExactMarginalLogLikelihood(self.likelihood, self.model)
        self._invalidate_kernel_cache()
  

    def _ensure_sim_state(self):
        if not hasattr(self, "_sim_events"):
            self._sim_events = deque(maxlen=4000)  # (s, added?, abs_res, T)
        if not hasattr(self, "_s_thresh_dyn"):
            self._s_thresh_dyn = 0.4  # sensible starting default in [0,1]
        if not hasattr(self, "prev_s_event"):
            self.prev_s_event = None
        if not hasattr(self, "prev_abs_res_event"):
            self.prev_abs_res_event = None
        # if not hasattr(self, "r_nn_med"):
        #     self._r_nn_med = 1.0

    from contextlib import contextmanager
    @contextmanager
    def jitter_ctx(self, val=1e-3):
        try:
            with gpytorch.settings.cholesky_jitter(val):
                yield
        except TypeError:
            # Older GPyTorch: (float, double, half)
            with gpytorch.settings.cholesky_jitter(float(val), float(val), float(val)):
                yield

    def _train_once(
        self,
        X_tr,
        y_tr,
        *,
        adam_steps: int = 80,
        adam_lr: float = 2e-2,
        lbfgs_lr: float = 0.2,
        lbfgs_max_iter: int = 80,
        jitter: float = 1e-3,
        clip_grad: float = 1.0,
        reinit_ls_if_at_bounds: bool = True,
    ):
        """One warm-up (Adam) + polish (LBFGS) pass.

        Expects X_tr, y_tr already standardized via _transform().
        """
        self.model.train(); self.likelihood.train()

        # ---- Small noise floor; freeze during Adam
        with torch.no_grad():
            try:
                self.likelihood.noise = torch.tensor(1e-3, dtype=self.dtype, device=self.device)
            except AttributeError:
                pass
        if hasattr(self.likelihood, "noise_covar"):
            self.likelihood.noise_covar.raw_noise.requires_grad_(False)

        # Optional: if all LS are at/near bounds, nudge them to sensible start
        def _nudge_lengthscales():
            base = self.model.covar_module.base_kernel
            # collect leaf kernels that have .lengthscale
            leaves = []
            def walk(k):
                if hasattr(k, "kernels"):  # Sum/Prod
                    for kk in k.kernels: walk(kk)
                elif hasattr(k, "lengthscale"):
                    leaves.append(k)
            walk(base)
            with torch.no_grad():
                # initialize around 1.0 in standardized space
                for k in leaves:
                    ls = k.lengthscale
                    lo, hi = None, None
                    try:
                        const = k.param_transform.constraints["raw_lengthscale"]
                        lo = float(const.lower_bound.detach().cpu())
                        hi = float(const.upper_bound.detach().cpu())
                    except Exception:
                        pass
                    target = torch.full_like(ls, 1.0)
                    if reinit_ls_if_at_bounds and lo is not None and hi is not None:
                        at_lo = (ls <= lo * 1.001).all()
                        at_hi = (ls >= hi * 0.999).all()
                        if at_lo or at_hi:
                            k.lengthscale = target
                    else:
                        # single nudge if wildly out of range
                        if (ls > 10.0).any() or (ls < 0.05).any():
                            k.lengthscale = target

        _nudge_lengthscales()

        mll = self.mll

        with self.jitter_ctx(jitter):
            # ---------- Adam warm-up (model-only) ----------
            adam = torch.optim.Adam(self.model.parameters(), lr=adam_lr)
            for _ in range(adam_steps):
                adam.zero_grad(set_to_none=True)
                out  = self.model(X_tr)
                loss = -mll(out, y_tr)
                loss.backward()
                if clip_grad is not None:
                    torch.nn.utils.clip_grad_norm_(self.model.parameters(), max_norm=clip_grad)
                adam.step()

            # ---------- LBFGS polish (model + likelihood) ----------
            # Unfreeze noise now so it can move, and include likelihood params
            if hasattr(self.likelihood, "noise_covar"):
                self.likelihood.noise_covar.raw_noise.requires_grad_(True)

            # Build optimizer over BOTH sets of params
            params = list(self.model.parameters()) + list(self.likelihood.parameters())
            lbfgs = torch.optim.LBFGS(
                params, lr=lbfgs_lr, max_iter=lbfgs_max_iter, line_search_fn="strong_wolfe"
            )

            def closure():
                lbfgs.zero_grad(set_to_none=True)
                out  = self.model(X_tr)
                loss = -mll(out, y_tr)
                if not torch.isfinite(loss):
                    # keep graph valid; tiny positive loss
                    dummy = torch.tensor(0.0, dtype=y_tr.dtype, device=y_tr.device, requires_grad=True)
                    return dummy + 1e-6
                loss.backward()
                return loss

            final_loss = lbfgs.step(closure)

        self.model.eval(); self.likelihood.eval()
            

        return float(final_loss.detach().cpu()) if torch.is_tensor(final_loss) else float(final_loss)



    def _assert_and_clean_finite(self):
        # Ensure shapes
        self.X_raw = self._ensure_2d(self.X_raw)
        self.y_raw = self._ensure_1d(self.y_raw)

        # ---------- 1) ROW filter ----------
        row_ok_X = torch.isfinite(self.X_raw).all(dim=1)   # (N,)
        row_ok_y = torch.isfinite(self.y_raw)              # (N,)
        row_ok   = row_ok_X & row_ok_y
        if row_ok.sum() == 0:
            raise RuntimeError("No valid rows in (X,y): all rows contain non-finite values.")
        if row_ok.sum() < self.X_raw.shape[0]:
            self.X_raw = self.X_raw[row_ok]
            self.y_raw = self.y_raw[row_ok]

        N, D = self.X_raw.shape

        # ---------- 2) COLUMN mask (decide BEFORE sanitizing) ----------
        col_finite_ok = torch.isfinite(self.X_raw).all(dim=0)   # (D,)
        if not col_finite_ok.any():
            raise RuntimeError("No finite feature columns in X_raw.")

        if N >= 3:
            med = self.X_raw.median(dim=0).values
            mad = (self.X_raw - med).abs().median(dim=0).values
            col_var_ok = (1.4826 * mad) > 1e-8
            col_ok = col_finite_ok & col_var_ok
            if not col_ok.any():

                col_ok = col_finite_ok
        else:

            col_ok = col_finite_ok


        self.X_raw = torch.nan_to_num(self.X_raw, nan=0.0, posinf=1e12, neginf=-1e12)
        self.y_raw = torch.nan_to_num(self.y_raw, nan=0.0, posinf=1e12, neginf=-1e12)


        if getattr(self, "_keep_cols", None) is None:

            self._keep_cols = col_ok
        else:
            if self._keep_cols.numel() != D:
                raise RuntimeError(
                    f"Frozen mask length {self._keep_cols.numel()} != X_raw width {D}"
                )



    def _fit_data(self):
        self._assert_and_clean_finite()
        X_view = self._apply_keep_cols(self.X_raw)
        self.X_mean = X_view.mean(dim=0)
        self.X_std  = self._compute_feature_scale(X_view)

        self.y_mean = self.y_raw.mean()
        if self.y_raw.numel() < 3:
            # center only; no scaling
            self.y_std = torch.tensor(1.0, dtype=self.dtype, device=self.device)
        else:
            self.y_std = self.y_raw.std(unbiased=False).clamp_min(1e-6)

    def _compute_feature_scale(self, X):
        eps = 1e-12
        if self.scale_mode == "zscore":
            std = X.std(dim=0, unbiased=False)
        elif self.scale_mode == "robust":
            med = X.median(dim=0).values
            mad = (X - med).abs().median(dim=0).values   # median absolute deviation
            std = 1.4826 * mad                           # ~ N(0,1) consistency
        elif self.scale_mode == "phys":
            # Fallback: if you know your feature layout, set per-feature floors
            # e.g., bonds (Δr in bohr) -> 0.2 bohr; angles cosθ/sinθ -> 0.2; dihedral cos/sin -> 0.2.
            # Here, just use MAD as baseline:
            med = X.median(dim=0).values
            mad = (X - med).abs().median(dim=0).values
            std = 1.4826 * mad
        else:
            std = X.std(dim=0, unbiased=False)
        # Clamp by a *meaningful* floor to avoid exploding z-scores
        floor = torch.full_like(std, self.x_floor)
        std = torch.maximum(std, floor)
        return std
    
    def _ensure_2d(self, X):
        X = torch.as_tensor(X, dtype=self.dtype, device=self.device)
        if X.ndim == 1:
            X = X.unsqueeze(0)
        return X

    def _ensure_1d(self, y):
        y = torch.as_tensor(y, dtype=self.dtype, device=self.device)
        return y.reshape(-1)

    def _apply_keep_cols(self, X):
        X = self._ensure_2d(X)
        if getattr(self, "_keep_cols", None) is not None:
            if self._keep_cols.numel() == X.shape[1]:
                X = X[:, self._keep_cols]
            elif X.shape[1] == int(self._keep_cols.sum().item()):
                # Already masked; leave it
                pass
            else:
                raise RuntimeError(
                    f"Mask length {self._keep_cols.numel()} incompatible with X width {X.shape[1]}"
                )
        return X

    
    def _transform(self, X, y=None):
        X = self._apply_keep_cols(X)
        if not self.standardize:
            if y is None:
                return X
            return X, self._ensure_1d(y)

        Xn = (X - self.X_mean) / self.X_std
        if y is None:
            return Xn
        yn = (self._ensure_1d(y) - self.y_mean) / self.y_std
        return Xn, yn

    def replace_data(self, X_new, y_new, *, clone=True):
        """
        Replace the current training data with (X_new, y_new).
        - Ensures shapes
        - Preserves device/dtype
        - Updates GPyTorch ExactGP model if present
        """
        # Shape checks
        X_new = self._ensure_2d(X_new)
        y_new = self._ensure_1d(y_new)
        if X_new.size(0) != y_new.size(0):
            raise ValueError(f"X_new and y_new must have same length: {X_new.size(0)} vs {y_new.size(0)}")

        # Preserve device/dtype from existing tensors if available
        device = getattr(self, "X_raw", X_new).device if hasattr(self, "X_raw") else getattr(self, "device", None)
        if device is None:
            device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        dtype = getattr(self, "X_raw", X_new).dtype if hasattr(self, "X_raw") else torch.get_default_dtype()

        # Convert and (optionally) clone to decouple external views
        X_new = torch.as_tensor(X_new, device=device, dtype=dtype)
        y_new = torch.as_tensor(y_new, device=device, dtype=dtype)
        if clone:
            X_new = X_new.clone()
            y_new = y_new.clone()

        with torch.no_grad():
            self.X_raw = X_new.contiguous()
            self.y_raw = y_new.contiguous()

            # # If you have a GPyTorch ExactGP model attached
            # self._fit_data()
            # X_tr, y_tr = self._transform(self.X_raw, self.y_raw)
            # if hasattr(self, "model"):
            #     self.model.set_train_data(inputs=X_tr, targets=y_tr, strict=False)

        # If you cache kernels/factorizations, invalidate them here
        if hasattr(self, "_cached_chol"):
            self._cached_chol = None
        self._invalidate_kernel_cache()


    def refit(self, steps="auto", verbose=True):
        
        if verbose:
            print("y_std (raw):", self.y_std)
            print("X_std (raw) min/med/max:", self.X_std.min(), self.X_std.median(), self.X_std.max())
            print("lengthscale:", self.model.covar_module.base_kernel.lengthscale)
            print("outputscale:", self.model.covar_module.outputscale)
            print("noise:", self.likelihood.noise)

        try:
            old_keep = int(self._keep_cols.sum()) if getattr(self, "_keep_cols", None) is not None else None
        except Exception:
            old_keep = None
        # self._keep_cols = None   # <- let _fit_data recompute keep mask from current X_raw
        self._fit_data()         # now scalers *and* mask match the full dataset

        if verbose:
            new_keep = int(self._keep_cols.sum()) if self._keep_cols is not None else -1
            print(f"kept dims: {old_keep} -> {new_keep}")
            print("y_std:", float(self.y_std))

        X_tr, y_tr = self._transform(self.X_raw, self.y_raw)
        X_tr = X_tr.to(self.device, dtype=self.dtype)
        y_tr = y_tr.to(self.device, dtype=self.dtype)

        # Track dataset growth
        if not hasattr(self, "_last_N"): self._last_N = 0
        lastN = int(self._last_N); Nnow = int(self.X_raw.shape[0])
        big_growth = (Nnow >= 50 and (Nnow > max(1.2 * lastN, lastN + 10)))

        if big_growth:
            base = self.model.covar_module.base_kernel
            groups = base.kernels if hasattr(base, "kernels") else [base]
            with torch.no_grad():
                for k in groups:
                    dims = k.active_dims
                    if dims is None or len(dims) == 0: continue
                    D = torch.cdist(X_tr[:, dims], X_tr[:, dims], p=2)
                    D.fill_diagonal_(float('inf'))
                    nn_med = D.min(dim=1).values.median().clamp_min(1e-3)
                    ell = nn_med.clamp(0.05, 10.0)   # respect your Interval constraint
                    k.lengthscale = ell              # set via constrained param
                # neutral variance/noise in standardized space (optional)
                self.model.covar_module.outputscale = y_tr.var().clamp_min(1e-8)
                self.likelihood.noise = (0.05 * y_tr.std()).pow(2).clamp_min(1e-9)

        self._last_N = Nnow

        # (re)attach train data and optimize as you already do
        self.model.set_train_data(X_tr, y_tr, strict=False)
        self._invalidate_kernel_cache()


        if steps == "auto":
            self._train_once(X_tr, y_tr, adam_steps=80, lbfgs_max_iter=80)

        with torch.no_grad(), gpytorch.settings.fast_pred_var():
            pred = self.likelihood(self.model(X_tr))
            mu  = pred.mean
            var = pred.variance.clamp_min(1e-12)
            z   = (y_tr - mu) / var.sqrt()

        # Perfect calibration: z ~ N(0,1)
        z_abs = z.abs()
        cover_1sigma = (z_abs <= 1.0).double().mean().item()  # expect ~0.68
        cover_2sigma = (z_abs <= 2.0).double().mean().item()  # expect ~0.95
        nll = 0.5*( (y_tr-mu)**2/var + var.log() + np.log(2*np.pi) ).mean().item()
        print(f"1σ cover={cover_1sigma:.2f}, 2σ cover={cover_2sigma:.2f}, mean NLL={nll:.3f}")

        # self.model.eval(); self.likelihood.eval()
        if verbose:
            print("lengthscale:", self.model.covar_module.base_kernel.lengthscale)
            print("outputscale:", self.model.covar_module.outputscale)
            print("noise:", self.likelihood.noise)
        
        self._r_nn_med = float('nan')
        prev = getattr(self, "_d_thresh_dyn", None)
        self._r_nn_med = float('nan')                 # force recompute of median NN in new metric
        if hasattr(self, "_dist_events"): self._dist_events.clear()
        if prev is None:
            print('prev', prev)
            self._d_thresh_dyn = 1.0 


        print("#kept dims:", int(self._keep_cols.sum()) if self._keep_cols is not None else -1)
        print("X_std min/median/max:", self.X_std.min().item(), self.X_std.median().item(), self.X_std.max().item())

    
    @torch.no_grad()
    def calibration_report(self):
        Xs, ys = self._transform(self.X_raw, self.y_raw)
        self.model.eval(); self.likelihood.eval()
        with gpytorch.settings.fast_pred_var():
            pred = self.likelihood(self.model(Xs))
        mu  = pred.mean
        var = pred.variance.clamp_min(1e-12)
        z   = (ys - mu) / var.sqrt()
        cover1 = (z.abs() <= 1.0).double().mean().item()
        cover2 = (z.abs() <= 2.0).double().mean().item()
        nll = 0.5 * (((ys - mu)**2 / var) + var.log() + math.log(2*math.pi)).mean().item()
        print(f"[Calib] 1σ={cover1:.2f} (→0.68), 2σ={cover2:.2f} (→0.95), mean NLL={nll:.3f}")


    def kernel_summary(self, *, verbose=True, rtol_at_bound=5e-3):
        """
        Inspect group kernels inside self.model.covar_module.base_kernel (a SumKernel)
        and report:
        - per-group lengthscale (in *standardized* feature units)
        - active_dims (feature indices)
        - outputscale and noise
        - whether any lengthscale is close to its bounds
        Returns a dict you can log or assert on.
        """
        base = self.model.covar_module.base_kernel  # likely a SumKernel
        items = []

        def _walk(kern, path=""):
            # Recurse into Sum/Prod kernels; collect leaf kernels with lengthscale
            if hasattr(kern, "kernels"):  # SumKernel/ProdKernel
                for i, kk in enumerate(kern.kernels):
                    _walk(kk, path=f"{path}.{type(kern).__name__}[{i}]")
            else:
                name = type(kern).__name__
                # Some kernels (e.g., Matern/RBF) have .lengthscale, others don't
                ls = getattr(kern, "lengthscale", None)
                ad = getattr(kern, "active_dims", None)
                bounds = None
                at_lo = at_hi = None
                if ls is not None:
                    # lengthscale is constrained; we can read its current value
                    ls_val = ls.detach().cpu().double().view(-1).numpy()  # shape (1,) if shared; (Dg,) if ARD
                    # If constraints were applied, you can find them under param_name "raw_lengthscale"
                    try:
                        const = kern.param_transform.constraints["raw_lengthscale"]
                        lo = float(const.lower_bound.detach().cpu())
                        hi = float(const.upper_bound.detach().cpu())
                        bounds = (lo, hi)
                        at_lo = bool((ls_val <= lo * (1 + rtol_at_bound)).all())
                        at_hi = bool((ls_val >= hi * (1 - rtol_at_bound)).all())
                    except Exception:
                        pass

                    items.append({
                        "path": path.strip(".") or name,
                        "kernel": name,
                        "active_dims": list(ad.tolist()) if ad is not None else None,
                        "lengthscale": ls_val.tolist(),
                        "bounds": bounds,
                        "near_lower_bound": at_lo,
                        "near_upper_bound": at_hi,
                    })

        _walk(base)

        out = {
            "groups": items,
            "outputscale": float(self.model.covar_module.outputscale.detach().cpu()),
            "noise": float(self.likelihood.noise.detach().cpu()),
        }

        if verbose:
            print("=== Kernel hyperparameters (standardized feature units) ===")
            for it in items:
                ad_str = "all" if it["active_dims"] is None else it["active_dims"]
                ls_str = ", ".join(f"{v:.4g}" for v in it["lengthscale"])
                b_str  = f"[{it['bounds'][0]:.3g}, {it['bounds'][1]:.3g}]" if it["bounds"] else "[-, -]"
                flags = []
                if it["near_lower_bound"]: flags.append("near-LOWER-BOUND")
                if it["near_upper_bound"]: flags.append("near-UPPER-BOUND")
                flag_str = ("  <-- " + ", ".join(flags)) if flags else ""
                print(f"{it['kernel']:>12s} dims={ad_str}  ℓ={ls_str}  bounds={b_str}{flag_str}")
            print(f"Outputscale: {out['outputscale']:.4g}")
            print(f"Noise:       {out['noise']:.4g}")

            # Bonus sanity: pairwise dist in the standardized space
            X_tr_std, _ = self._transform(self.X_raw, self.y_raw)  # returns (X_std, y_std)
            with torch.no_grad():
                d = torch.cdist(X_tr_std, X_tr_std, p=2)
                d.fill_diagonal_(float('inf'))
                print(f"Min pairwise dist (std space): {float(torch.min(d)):.4g}")

        return out


    @torch.no_grad()
    def predict(self, X_query, return_std=True, verbose=False):
        Xq = self._transform(X_query)  # masked + standardized
        # Safety check: active dims union must be within Xq width
        Dq = Xq.shape[-1]
        base = self.model.covar_module.base_kernel
        # base may be a SumKernel; collect all active_dims
        def _collect_active_dims(k):
            if hasattr(k, "kernels"):  # Sum/Prod kernels
                dims = []
                for kk in k.kernels:
                    dims += _collect_active_dims(kk)
                return dims
            return list(getattr(k, "active_dims", torch.arange(Dq)).tolist())
        all_dims = _collect_active_dims(base)
        assert max(all_dims, default=-1) < Dq, f"X_query has width {Dq}, but kernel expects dim {max(all_dims)}."

        self.model.eval(); self.likelihood.eval()
        with torch.no_grad(), gpytorch.settings.fast_pred_var():
            pred = self.likelihood(self.model(Xq))

        # Unstandardize
        if self.standardize:
            mean = pred.mean * self.y_std + self.y_mean
            var  = pred.variance * (self.y_std ** 2)
        else:
            mean, var = pred.mean, pred.variance

        mean = mean.detach().cpu().numpy()
        std  = var.sqrt().detach().cpu().numpy() if return_std else None
        return mean, std
    
    @torch.no_grad()
    def predict_latent(self, X_query, return_std=True, diagnostics=False):
        # same preprocessing
        Xq = self._transform(X_query)
        self.model.eval()
        if diagnostics or self.verbose_diagnostics:
            self.kernel_summary()
            self.calibration_report()
        with gpytorch.settings.fast_pred_var():
            fdist = self.model(Xq)  # LATENT: no likelihood -> epistemic only
        if not return_std:
            return fdist.mean.detach().cpu().numpy(), None

        if self.standardize:
            mean = fdist.mean * self.y_std + self.y_mean
            std  = fdist.variance.sqrt() * self.y_std
        else:
            mean, std = fdist.mean, fdist.variance.sqrt()

        return mean.detach().cpu().numpy(), std.detach().cpu().numpy()


    def reset_distance_history(self):
        if hasattr(self, "_safe_dists"): del self._safe_dists
        if hasattr(self, "_d_thresh_dyn"): del self._d_thresh_dyn


    def record_similarity_event(self, s_value: float, added: bool,
                            abs_residual: float | None, T: float | None):
        """Call once per decision after you know whether you added a point."""
        if not hasattr(self, "_sim_events"):
            self._sim_events = deque(maxlen=4000)
        self._sim_events.append({
            "s": float(s_value),
            "added": bool(added),
            "abs_res": float(abs_residual) if abs_residual is not None else np.inf,
            "T": float(T) if T is not None else 1.0,
        })

    def update_similarity_threshold(self,
                                q_good=0.25, q_bad=0.75,
                                min_good=1, min_bad=2,
                                cap_min=0.20, cap_max=0.95,
                                prior_thresh=0.40,   # your conservative prior
                                prior_strength=8,    # how many pseudo-events the prior counts for
                                alpha_base=0.20, alpha_max=0.35,
                                max_drop=0.02):
        if not getattr(self, "_sim_events", None):
            return
        
        # -------- robust small-n quantiles --------
        def robust_q(arr, q):
            arr = np.asarray(arr).ravel()
            n = arr.size
            arr.sort()
            if n <= 4:
                # no trimming; linear interp
                if n == 1:
                    return float(arr[0])
                idx = q * (n - 1)
                lo, hi = int(np.floor(idx)), int(np.ceil(idx))
                return float(arr[lo] + (idx - lo) * (arr[hi] - arr[lo]))
            # light, adaptive trimming only when we have enough samples
            lo_q = 0.0 if n < 20 else 0.05
            hi_q = 1.0 if n < 200 else 0.995
            lo, hi = np.quantile(arr, [lo_q, hi_q])
            arr = arr[(arr >= lo) & (arr <= hi)]
            return float(np.quantile(arr, q))

        # unpack
        s = np.array([e["s"] for e in self._sim_events], float)
        added = np.array([e["added"] for e in self._sim_events], bool)
        res = np.array([e["abs_res"] for e in self._sim_events], float)
        T   = np.array([e["T"] for e in self._sim_events], float)

        good_mask = (~added) & np.isfinite(res) & (res <= 0.25 * T)
        bad_mask  = added

        good = s[good_mask]
        bad  = s[bad_mask]


        n_good, n_bad = good.size, bad.size
        n_tot = n_good + n_bad

        if n_good < min_good or n_bad < min_bad:
            # With zero bads, still allow gentle learning from good-only evidence:
            # keep threshold near prior but nudge toward a higher quantile of "good".
            if n_good >= min_good and n_bad == 0:
                s_lo_good = robust_q(good, 0.25 if n_good >= 5 else 0.50)  # low quantile or median
                if s_lo_good is None:
                    return
                cand_evidence = float(np.clip(s_lo_good, cap_min, cap_max))
                prev = getattr(self, "_s_thresh_dyn", prior_thresh)

                # Do not increase when we only have good evidence
                cand = min(prev, cand_evidence)

                # gentle smoothing + symmetric caps
                cand = np.clip(cand, prev - max_drop, prev + 0.0)  # no upward move
                alpha = min(alpha_max, alpha_base * n_good / (n_good + 10))
                self._s_thresh_dyn = (1 - alpha) * prev + alpha * cand
                return


        # adapt quantiles for tiny samples
        qg = q_good if n_good >= 5 else 0.50         # median when few good
        qb = q_bad  if n_bad  >= 5 else 1.00         # use max(bad) when few bad

        s_lo_good = robust_q(good, qg)
        s_hi_bad  = robust_q(bad,  qb)

        # add a small safety margin when bad’s sample is tiny (pushes boundary up)
        if n_bad <= 4:
            s_hi_bad += 0.02 * (5 - n_bad)  # +0.08, +0.06, +0.04, +0.02

        cand_evidence = max(s_lo_good, s_hi_bad)
        cand_evidence = float(np.clip(cand_evidence, cap_min, cap_max))

        # -------- shrinkage toward a prior when data are scarce --------
        w_prior = prior_strength / (prior_strength + n_tot)
        cand = w_prior * prior_thresh + (1 - w_prior) * cand_evidence
        max_rise = 0.02
        prev = getattr(self, "_s_thresh_dyn", prior_thresh)
        cand = np.clip(cand, prev - max_drop, prev + max_rise)   # e.g. max_rise=0.02
        alpha = min(alpha_max, alpha_base * n_tot / (n_tot + 10))
        self._s_thresh_dyn = (1 - alpha) * prev + alpha * cand

    

    def _ensure_dist_state(self):
        if not hasattr(self, "_dist_events"):
            self._dist_events = deque(maxlen=4000)
        if not hasattr(self, "_d_thresh_dyn"):
            self._d_thresh_dyn = 1.0   # ← normalized prior, not 3.0
        if not hasattr(self, "prev_d_event"):
            self.prev_d_event = None

    def _metric_whiten(self, X_std, groups):
        """Return whitened, weighted features Y so that ||Y_i - Y_j|| = r_eff(x_i,x_j)."""
        base = self.model.covar_module.base_kernel
        group_kerns = getattr(base, "kernels", [base])  # SumKernel → list
        J = X_std.shape[1]
        # weights a_j initialized to zero
        a = torch.zeros(J, dtype=X_std.dtype, device=X_std.device)
        ls = torch.ones(J, dtype=X_std.dtype, device=X_std.device)
        G = 0
        for ker in group_kerns:
            ad = getattr(ker, "active_dims", None)
            if ad is None or len(ad) == 0:
                continue
            idx = torch.as_tensor(list(ad), dtype=torch.long, device=X_std.device)
            ell = ker.lengthscale.detach().to(X_std.device).reshape(-1)  # (|g|,) if ARD, else (1,)
            if ell.numel() == 1:
                ell = ell.repeat(idx.numel())
            ls[idx] = ell
            a[idx] += 1.0 / max(1, idx.numel())
            G += 1
        if G == 0:
            # fallback: plain standardized RMS distance
            a[:] = 1.0 / J
        else:
            a /= G
        scale = torch.sqrt(torch.clamp(a, min=1e-12)) / torch.clamp(ls, min=1e-12)
        
        return X_std * scale  # Y

    def _invalidate_kernel_cache(self):
        self._diag_cache_version = getattr(self, "_diag_cache_version", 0) + 1
        self._kernel_diag_cache = {}

    def _get_group_diag(self, kernel, Xtr_std):
        cache_key = (id(kernel), self._diag_cache_version, Xtr_std.shape[0])
        diag = self._kernel_diag_cache.get(cache_key)
        if diag is None:
            with torch.no_grad():
                diag = torch.diag(kernel(Xtr_std, Xtr_std).evaluate()).clamp_min(1e-12)
            self._kernel_diag_cache[cache_key] = diag
        return diag

    def _nn_r_stats(self, Xtr_std, groups):
        Y = self._metric_whiten(Xtr_std, groups)
        D = torch.cdist(Y, Y)          # r_eff between all pairs
        D.fill_diagonal_(float('inf'))
        r_nn = D.min(dim=1).values.cpu().numpy()
        return float(np.median(r_nn)), float(np.quantile(r_nn, 0.1)), float(np.quantile(r_nn, 0.9))
    
    def _ensure_r_ref(self, Xtr_std, groups):
        if not hasattr(self, "_r_nn_med") or not np.isfinite(self._r_nn_med):
            m, *_ = self._nn_r_stats(Xtr_std, groups)
            self._r_nn_med = max(1e-8, m)   # avoid divide-by-zero

    def record_distance_event(self, d_norm, added, abs_residual, T):
        self._ensure_dist_state()
        d  = float(d_norm) if d_norm is not None else np.inf
        ar = float(abs_residual) if (abs_residual is not None) else np.inf
        t  = float(T) if (T is not None) else 1.0
        self._dist_events.append((d, bool(added), ar, t, str(self.cause)))

    def _compute_d_norm(self, Xq_std: torch.Tensor, Xtr_std: torch.Tensor, groups: dict) -> float:
        """
        Normalized distance:
        1) whiten both Xq and Xtr into Y via _metric_whiten (ARD RMS metric),
        2) r_raw = min_i ||Yq - Yt[i]||_2,
        3) d_norm = r_raw / median_NN(Y_tr)  (cached in self._r_nn_med).
        Returns a dimensionless number where ~1.0 ≈ 'typical training NN'.
        """
        # Make sure the reference scale exists (computed in the SAME metric)
        self._ensure_r_ref(Xtr_std, groups or {})

        # Raw distance in the whitened space
        Yq = self._metric_whiten(Xq_std, groups)
        Yt = self._metric_whiten(Xtr_std, groups)
        r_raw = float(torch.cdist(Yq, Yt).min().cpu())

        # Normalize
        return r_raw / float(self._r_nn_med)

    def update_distance_threshold(self,
        q_good=0.95, q_bad=0.25,
        min_good=2, min_bad=2,
        cap_min=0.9, cap_max=1.8,           # wider band
        prior_thresh=1.2, prior_strength=12, # more realistic prior
        alpha_base=0.20, alpha_max=0.35,
        max_rise=0.10, max_drop=0.10):

        self._ensure_dist_state()
        if not self._dist_events: return

        d, added, ar, t, cause = map(np.asarray, zip(*self._dist_events))
        good_mask = (~added) & np.isfinite(ar) & (ar <= 0.25*t)
        # count only distance-caused adds if you have it; otherwise fall back to 'added'
        bad_mask  = (added) & (np.asarray(cause, object) == "distance")

        good = d[good_mask]; bad = d[bad_mask]
        n_good, n_bad = good.size, bad.size
        n_tot = n_good + n_bad
        prev = getattr(self, "_d_thresh_dyn", prior_thresh)

        def robust_q(arr, q):
            arr = np.asarray(arr).ravel()
            if arr.size == 0: return None
            arr.sort(); n = arr.size
            if n <= 4:
                if n == 1: return float(arr[0])
                idx = q*(n-1); lo, hi = int(np.floor(idx)), int(np.ceil(idx))
                return float(arr[lo] + (idx-lo)*(arr[hi]-arr[lo]))
            lo_q = 0.0 if n < 20 else 0.05
            hi_q = 1.0 if n < 200 else 0.995
            lo, hi = np.quantile(arr, [lo_q, hi_q])
            arr = arr[(arr >= lo) & (arr <= hi)]
            return float(np.quantile(arr, q))

        # Good-only: let the max distance creep up toward the good tail
        print('n good, min good', n_good, min_good)
        if n_good >= min_good and n_bad == 0:
            d_hi_good = robust_q(good, 0.80 if n_good < 5 else q_good)
            if d_hi_good is None: return
            cand_evidence = float(np.clip(d_hi_good + 0.05, cap_min, cap_max))
            w_prior = prior_strength / (prior_strength + max(1, n_good))
            cand = w_prior*prior_thresh + (1 - w_prior)*cand_evidence
            cand = np.clip(cand, prev - max_drop, prev + max_rise)
            alpha = min(alpha_max, alpha_base * n_tot / (n_tot + 10))
            self._d_thresh_dyn = (1 - alpha)*prev + alpha*cand

            print('In first check', self._d_thresh_dyn)
            return

        if n_good < min_good or n_bad < min_bad: return

        d_hi_good = robust_q(good, q_good)
        d_lo_bad  = robust_q(bad,  q_bad)
        if d_hi_good is None or d_lo_bad is None: return
        if n_bad <= 4: d_lo_bad -= 0.05*(5 - n_bad)

        # balanced boundary (less strict than min(...))
        cand_evidence = float(np.clip(0.5*(d_hi_good + d_lo_bad), cap_min, cap_max))
        w_prior = prior_strength / (prior_strength + n_tot)
        cand = w_prior*prior_thresh + (1 - w_prior)*cand_evidence
        cand = np.clip(cand, prev - max_drop, prev + max_rise)
        alpha = min(alpha_max, alpha_base * n_tot / (n_tot + 10))
        self._d_thresh_dyn = (1 - alpha)*prev + alpha*cand
        print('>t the end', self._d_thresh_dyn)
        

    def get_distance_threshold(self, default=1.0) -> float:

        self._ensure_dist_state()
        return float(self._d_thresh_dyn if hasattr(self, "_d_thresh_dyn") else default)

    def get_similarity_threshold(self, default=0.2):
        # default in your target units (kcal/mol per-atom if that's your y)
        if hasattr(self, "_s_thresh_dyn"):
            return float(self._s_thresh_dyn)
        return float(default)
    
    def get_sigma_threshold(self, default=0.20):
        # default in your target units (kcal/mol per-atom if that's your y)
        if hasattr(self, "_sigma_thresh_dyn"):
            return float(self._sigma_thresh_dyn)
        return float(default)
    
    @torch.no_grad()
    def should_qm_check(self, X_query, T=0.03, k=2.0, sigma_high=None,
                    d_policy="percentile",  # or "fixed"
                    d_fixed=2.5,
                    groups=None,
                    use_effective_dist=True):
        # latent prediction (epistemic only)
        Xq_std = self._transform(X_query)
        with gpytorch.settings.fast_pred_var():
            fdist = self.model(Xq_std)

        # Unstandardize
        if self.standardize:
            mu = fdist.mean * self.y_std + self.y_mean
            sig = fdist.variance.sqrt() * self.y_std
        else:
            mu = fdist.mean
            sig = fdist.variance.sqrt()

        mu = float(mu.squeeze().cpu())
        sig = float(sig.squeeze().cpu())

        # distance guard
        Xtr_std, _ = self._transform(self.X_raw, self.y_raw)
        d_norm = self._compute_d_norm(Xq_std, Xtr_std, groups) if use_effective_dist else \
             float(torch.cdist(Xq_std, Xtr_std).min()) 
        
        d_thresh = self.get_distance_threshold(1.0)

        # bound = abs(mu) + k*sig
        # prior_std = float(self.model.covar_module.outputscale.sqrt().detach().cpu()) \
        #     * (float(self.y_std.cpu()) if self.standardize else 1.0)
        # sigma_high = 0.7 * prior_std
        # print(f"[decision] N={len(self.X_raw)}, prior_std={prior_std:.5f}, "
        #     f"sigma_high={sigma_high:.5f}, sigma_lat={sig:.5f}, "
        #     f"d_norm={d_norm:.3f}, d_thresh={getattr(self, '_d_thresh_dyn', None)}")

        # if bound > T:
        #     return True, {"cause":"bound", "mu":mu, "sigma_lat":sig, "bound":bound, "T":T, "d_norm":d_norm, "k":k}

        # 4) distance test (secondary; only if σ is high)
        # pick sigma_high ~ prior std if not provided
        # prior std ≈ sqrt(outputscale) * y_std
        prior_std = float(self.model.covar_module.outputscale.sqrt().detach().cpu()) * (float(self.y_std.cpu()) if self.standardize else 1.0)
        sigma_high = 0.9 * prior_std  # "near prior" threshold (70%)
        self._sigma_thresh_dyn = sigma_high

        sigma_near_prior = (sig >= sigma_high)

        # if d_policy == "fixed":
        #     d_thresh = d_fixed
        # else:
        #     # percentile policy: use a stored running threshold; fall back to fixed if missing
        #     d_thresh = getattr(self, "_d_thresh_dyn", d_fixed)
        
        print('sigma hihg ', sigma_high)
        
        # Secondary: only if uncertainty is near prior AND similarity is low
        # prior_std   = float(self.model.covar_module.outputscale.sqrt()) * float(self.y_std)
        # sigma_high  = 0.7 * prior_std
        
        s_thresh    = getattr(self, "_s_thresh_dyn", 0.2)  # learned percentile; default conservative
        s_max       = self.max_kernel_similarity(Xq_std, Xtr_std, topk=5)
        
        print('INFOOOO', s_thresh, s_max, d_thresh, d_norm)
        
        trigger = (sig >= sigma_high) and (s_max < s_thresh) 
        print('Kernel similarity', trigger, s_thresh, s_max)

        should_add_by_risk, diag = self.decide_add_by_risk(s_max, d_norm, sig,
                                                  R_add=0.60, R_clear=0.45)
        # Optionally combine with your residual rule:

        cause = (
            "distance"     if (sig >= sigma_high and d_norm > d_thresh) else
            "kernel"       if (s_max < s_thresh)                        else
            "uncertainty"  if (sig >= sigma_high)                       else
            "none"
        )

        self.cause = cause
        self.record_similarity_event(s_max, False, abs_residual=abs(s_max - s_thresh), T=T)
        self.record_distance_event(d_norm, False, abs_residual=abs(d_thresh - d_norm), T=T)

        if should_add_by_risk:
            return True, {"cause": cause, "s_max": s_max, "s_thresh": s_thresh,
                  "d_norm": d_norm, "d_thresh": d_thresh}

        
        # if sigma_near_prior and d_eff > d_thresh:
        #     return True, {"cause":"distance", "mu":mu, "sigma_lat":sig, "bound":bound, "T":T, "d_eff":d_eff, "k":k, "d_thresh":d_thresh}

        return False, {"cause":"none", "mu":mu, "sigma_lat":sig, "T":T, "d_norm":d_norm, "k":k, "d_thresh":d_thresh}
    

    def _norm_sim_group(self, ker, Xq_std, Xtr_std):
        with torch.no_grad():
            K_q_tr = ker(Xq_std, Xtr_std).evaluate()
            K_q_q  = ker(Xq_std, Xq_std).evaluate()[0, 0]
            diag   = self._get_group_diag(ker, Xtr_std)
            s = (K_q_tr / torch.sqrt(K_q_q * diag)).cpu().numpy().ravel()
        s = np.clip(s, 1e-9, 1.0)
        # temper by group size to offset high-D shrinkage
        ad = getattr(ker, "active_dims", None)
        gdim = (len(ad) if ad is not None else 1)
        if gdim > 1:
            s = s ** (1.0 / math.sqrt(gdim))
        return s 

    def max_kernel_similarity(self, Xq_std, Xtr_std, *, topk: int = 5,
                          agg: str = "geom", weight_mode: str = "inv_dim"):
        """
        Size-invariant similarity using your bond/angle/torsion kernels.
        agg: 'geom' (default), 'mean', or 'max'
        weight_mode: 'inv_dim' or 'inv_l2' (1/ell^2 if ARD/ell exists)
        """
        sims = []
        weights = []
        
        if not getattr(self, "group_kernels", None):
            # fallback: try to pull from model's sum kernel
            bk = getattr(self.model.covar_module, "base_kernel", None)
            if hasattr(bk, "kernels"):
                self.group_kernels = list(bk.kernels)
            else:
                # last resort: use the whole kernel (will be dimension-dependent)
                with torch.no_grad():
                    K_q_tr = self.model.covar_module.base_kernel(Xq_std, Xtr_std).evaluate()
                    K_q_q  = self.model.covar_module.base_kernel(Xq_std, Xq_std).evaluate()[0, 0]
                    K_trtr = self.model.covar_module.base_kernel(Xtr_std, Xtr_std).evaluate()
                    diag   = torch.diag(K_trtr).clamp_min(1e-12)
                    s_full = (K_q_tr / torch.sqrt(K_q_q * diag)).cpu().numpy().ravel()
                s_full = np.clip(s_full, 1e-12, 1.0)
                # power correction by number of groups to avoid collapse
                s_full = s_full ** (1.0 / max(1, getattr(self, "num_groups", 1)))
                s_full.sort()
                return float(np.mean(s_full[-min(topk, s_full.size):]))

        # per-group normalized similarities + weights
        for ker in self.group_kernels:
            s_g = self._norm_sim_group(ker, Xq_std, Xtr_std)
            sims.append(s_g)

            if weight_mode == "inv_l2" and hasattr(ker, "lengthscale"):
                ell = ker.lengthscale
                try:
                    ell_val = float(ell.mean().detach().cpu().numpy())
                except Exception:
                    ell_val = 1.0
                w = 1.0 / max(ell_val**2, 1e-6)
            else:
                # inverse by group size as a robust default
                ad = getattr(ker, "active_dims", None)
                gdim = int(len(ad)) if ad is not None else 1
                w = 1.0 / max(gdim, 1)
            weights.append(w)

        S = np.vstack(sims)            # (G, Ntr)
        w = np.asarray(weights, float)
        w = w / (w.sum() if w.sum() > 0 else 1.0)

        if S.size == 0:
            return 0.0

        if agg == "max":
            s_vec = S.max(axis=0)
        elif agg == "mean":
            s_vec = (w[:, None] * S).sum(axis=0)
        else:  # 'geom' (default)
            eps = 1e-9
            s_vec = np.exp((w[:, None] * np.log(np.clip(S, eps, 1.0))).sum(axis=0))

        s_vec.sort()
        k = min(topk, s_vec.size)
        return float(np.mean(s_vec[-k:]))
    
    # @torch.no_grad()
    # def max_kernel_similarity(self, Xq_std: torch.Tensor, Xtr_std: torch.Tensor, topk: int = 1) -> float:
    #     """
    #     Returns max (or top-K mean) *normalized* kernel similarity in [0,1].
    #     Uses the *base* kernel (sum of per-group Matérns) and normalizes by k(x,x).
    #     Xq_std: (1,D) standardized query (use the SAME transform as training)
    #     Xtr_std: (N,D) standardized training
    #     """
    #     base = self.model.covar_module.base_kernel  # SumKernel of group Matérns (no global outputscale)

    #     # k(xq, Xtr): shape (1, N)
    #     k_row = base(Xq_std, Xtr_std).evaluate().squeeze(0)  # (N,)

    #     # Normalize by k(xq,xq) so the diagonal equals 1
    #     k_self = base(Xq_std, Xq_std).evaluate().squeeze().item()  # ~= number of groups (e.g., 3.0)
    #     if k_self <= 0:
    #         return 0.0

    #     s = (k_row / k_self).clamp(min=0.0, max=1.0)  # (N,)

    #     if topk is None or topk <= 1:
    #         return float(s.max().cpu())
    #     else:
    #         topk = min(topk, s.numel())
    #         vals, _ = torch.topk(s, k=topk, largest=True)
    #         return float(vals.mean().cpu())

    
    def _collect_group_ls(self):
        """Return dict {'bond': ℓ_bond, 'angle': ℓ_angle, 'torsion': ℓ_torsion} from the summed kernel."""
        base = self.model.covar_module.base_kernel
        ls_vals = []
        def walk(k):
            if hasattr(k, "kernels"):
                for kk in k.kernels: walk(kk)
            elif hasattr(k, "lengthscale"):
                ls_vals.append(float(k.lengthscale.detach().cpu()))
        walk(base)
        # assume order: [bond, angle, torsion] as constructed
        keys = ["bond", "angle", "torsion"]
        return {k: ls_vals[i] for i, k in enumerate(keys) if i < len(ls_vals)}
    
    @torch.no_grad()
    
    def _min_effective_distance(self, Xq_std: torch.Tensor, Xtr_std: torch.Tensor, groups: dict) -> float:
        """
        Size-invariant, ARD-aware RMS effective distance:
        r_eff^2 = sum_g w_g * mean_j in g [ ((xq - xi)_j / ell_{g,j})^2 ]
        Returns min_i r_eff(xq, xi).
        """
        assert Xq_std.ndim == 2 and Xq_std.shape[0] == 1
        assert Xtr_std.ndim == 2 and Xtr_std.shape[1] == Xq_std.shape[1]

        Yq = self._metric_whiten(Xq_std, groups)
        Yt = self._metric_whiten(Xtr_std, groups)
        return float(torch.cdist(Yq, Yt).min().cpu())


    def _sigmoid(self, x):  # numerically stable logistic
        return 1.0 / (1.0 + np.exp(-np.clip(x, -20, 20)))

    def risk_score(self, s_max: float, d_norm: float, sigma: float,
                   w_s: float = 0.6, w_d: float = 0.3, w_u: float = 0.1,
                   tau_s: float | None = None,
                   tau_d: float | None = None,
                   tau_u: float | None = None):
        """
        Convert (similarity, distance, sigma) into a weighted risk in [0,1].
        Larger risk means 'add' is more likely.
        """
        s_thr = self.get_similarity_threshold(0.40)   # learned min similarity
        d_thr = self.get_distance_threshold(1.00)     # learned max distance
        u_thr = self.get_sigma_threshold(0.20)        # fixed or learned

        # Softness (how sharp the step is around each threshold)
        # Defaults scale with threshold to be roughly 5–10% band
        tau_s = tau_s if tau_s is not None else max(1e-6, 0.08 * s_thr)
        tau_d = tau_d if tau_d is not None else max(1e-6, 0.08 * d_thr)
        tau_u = tau_u if tau_u is not None else max(1e-6, 0.15 * u_thr)

        # Map to risks in [0,1]; each ≈ 0.5 when equal to its threshold
        s_risk = self._sigmoid((s_thr - s_max) / tau_s)   # lower s → higher risk
        d_risk = self._sigmoid((d_norm - d_thr) / tau_d)   # higher d → higher risk
        u_risk = self._sigmoid((sigma - u_thr) / tau_u)   # higher σ → higher risk

        # Normalize weights
        w_sum = max(1e-12, w_s + w_d + w_u)
        w_s, w_d, w_u = w_s / w_sum, w_d / w_sum, w_u / w_sum

        R = float(w_s * s_risk + w_d * d_risk + w_u * u_risk)
        parts = {"s_risk": float(s_risk), "d_risk": float(d_risk), "u_risk": float(u_risk),
                 "s_thr": s_thr, "d_thr": d_thr, "u_thr": u_thr}
        print('in risk score function', R, parts)
        return R, parts

    def decide_add_by_risk(self, s_max: float, d_norm: float, sigma: float,
                           R_add: float = 0.60, R_clear: float = 0.45,
                           ema_beta: float = 0.3,
                           hard_guards: dict | None = None):
        """
        Combine the three signals with weights; use hysteresis on an EMA of risk.
        - R_add: EMA risk above this → add
        - R_clear: EMA risk below this → safe
        - hard_guards: immediate add if any extreme rule trips
        """
        R, parts = self.risk_score(s_max, d_norm, sigma)

        # Hysteresis on an EMA to avoid flapping
        if not hasattr(self, "_risk_ema"):
            self._risk_ema = R
        else:
            self._risk_ema = (1 - ema_beta) * self._risk_ema + ema_beta * R

        add = False

        # Hard safety overrides (fail-fast)
        # Example defaults; tune for your scales
        guards = {
            "s_min": 0.10,          # if similarity below this → add
            "d_max": parts["d_thr"] + 0.10,  # if distance way above thr → add
            "u_max": parts["u_thr"] + 200000.0,  # if sigma way above thr → add
        }
        if hard_guards:
            guards.update(hard_guards)
        
        if d_norm < 1e-6:
            d_norm = 1e-6
        
        prop = s_max/d_norm 
        if parts["s_thr"] > 0.85:
            prop = s_max

        if s_max < guards["s_min"] or prop < parts["s_thr"] or d_norm > guards["d_max"] + 0.1: #or ((sigma > guards["u_max"]) and ((d_norm > guards["d_max"]) or (s_max < parts["s_thr"]))):
            add = True
        
        # print('Add before after the check in decider function', add)
        # # EMA-based decision with hysteresis
        # if not add:
        #     if self._risk_ema >= R_add:
        #         add = True
        #     elif self._risk_ema <= R_clear:
        #         add = False
        #     else:
        #         # keep previous state if you maintain one; default to no-add
        #         add = getattr(self, "_risk_armed", False)

        self._risk_armed = add  # remember last state for hysteresis
        diag = {"R": R, "R_ema": float(self._risk_ema)} | parts

        # print('DIAG', diag)
        return add, diag


    def corrected_energy(self, E_interp_query, X_query, threshold=None):
        dE_mean, dE_std = self.predict(X_query, return_std=True)
        E_hat = E_interp_query + dE_mean
        flag = (dE_std > threshold) if threshold is not None else None
        return E_hat, dE_std, flag
