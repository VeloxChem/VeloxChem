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
from scipy.optimize import linear_sum_assignment
from scipy.linalg import cho_factor, cho_solve
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
from.veloxchemlib import hartree_in_kcalpermol, bohr_in_angstrom
from .errorhandler import assert_msg_critical
from .inputparser import (parse_input, print_keywords, print_attributes)

with redirect_stderr(StringIO()) as fg_err:
    import geometric


class ExactDeltaE(gpytorch.models.ExactGP):
    def __init__(self, x_train, y_train, likl):
        super().__init__(x_train, y_train, likl)
        self.mean_module = gpytorch.means.ZeroMean()
        self.covar_module = gpytorch.kernels.ScaleKernel(
            gpytorch.kernels.MaternKernel(nu=2.5, ard_num_dims=x_train.shape[-1])
        )
    def forward(self, x):
        return gpytorch.distributions.MultivariateNormal(
            self.mean_module(x), self.covar_module(x)
        )

class GPRInterpolationDriver:

    def __init__(self, x0, y0, *, standardize=False):
        
        dtype = torch.double
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

        self.standardize = standardize
        self.dtype = dtype
        self.device = device

        self.X_raw = torch.as_tensor(x0, dtype=dtype, device=device)
        self.y_raw = torch.as_tensor(y0, dtype=dtype, device=device)

        self._fit_data()
        X, y = self._transform(self.X_raw, self.y_raw)
        
        self.likelihood = gpytorch.likelihoods.GaussianLikelihood().to(device).double()
        self.model = ExactDeltaE(X, y, self.likelihood).to(device).double()

        self.model.train(); self.likelihood.train()
        self.mll = gpytorch.mlls.ExactMarginalLogLikelihood(self.likelihood, self.model)
        # initial fit
        # fit_gpytorch_mll(self.mll)
    
    def _fit_data(self):
        if not self.standardize:
            self.X_mean = None
            self.y_mean = None

            self.X_std = None
            self.y_std = None
        
        eps = 1e-12
        self.X_mean = self.X_raw.mean(0)
        self.X_std  = self.X_raw.std(0).clamp_min(eps)
        self.y_mean = self.y_raw.mean(0)
        self.y_std  = self.y_raw.std(0).clamp_min(eps)

    def _transform(self, X, y=None):
        
        if not self.standardize:
            return (X, y) if y is not None else X
        
        Xn = (X - self.X_mean) / self.X_std
        if y is None: return Xn
        yn = (y - self.y_mean) / self.y_std
        
        return Xn, yn

    def add_data(self, X_new, y_new):
        """Append new (x, ΔE) pairs and keep model object; no re-init."""
        Xn = torch.as_tensor(X_new, dtype=self.dtype, device=self.device)
        yn = torch.as_tensor(y_new, dtype=self.dtype, device=self.device)
        self.X_raw = torch.cat([self.X_raw, Xn], dim=0)
        self.y_raw = torch.cat([self.y_raw, yn], dim=0)

    def refit(self, steps="auto"):
        """
        Recompute scalers, update train data in-place, and refit hyperparams.
        'steps="auto"' uses BoTorch’s LBFGS fit; or pass an int to do N Adam steps.
        """
        # Re-standardize with ALL data (most stable for small/medium n)
        self._fit_data()
        X_tr, y_tr = self._transform(self.X_raw, self.y_raw)

        # Update the model’s training tensors (no new object)
        self.model.set_train_data(inputs=X_tr, targets=y_tr, strict=False)

        # Warm-started optimization (params are already near a good solution)
        self.model.train(); self.likelihood.train()
        if steps == "auto":
            # fit_gpytorch_mll(self.mll)
            print('hello')
        else:
            opt = torch.optim.Adam(self.model.parameters(), lr=0.05)
            for _ in range(int(steps)):
                opt.zero_grad()
                out = self.model(X_tr)
                loss = -self.mll(out, y_tr)
                loss.backward()
                opt.step()


    @torch.no_grad()
    def predict(self, X_query, return_std=True):
        Xq_raw = torch.as_tensor(X_query, dtype=self.dtype, device=self.device)
        Xq = self._transform(Xq_raw)
        self.model.eval(); self.likelihood.eval()
        with gpytorch.settings.fast_pred_var():
            pred = self.likelihood(self.model(Xq))
        # unstandardize to physical units
        if self.standardize:
            mean = pred.mean * self.y_std + self.y_mean
            var  = pred.variance * (self.y_std ** 2)
        else:
            mean, var = pred.mean, pred.variance
        return (mean.cpu().numpy(),
                var.sqrt().cpu().numpy() if return_std else None)

    def corrected_energy(self, E_interp_query, X_query, threshold=None):
        dE_mean, dE_std = self.predict(X_query, return_std=True)
        E_hat = E_interp_query + dE_mean
        flag = (dE_std > threshold) if threshold is not None else None
        return E_hat, dE_std, flag