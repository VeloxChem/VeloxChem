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

from math import cos

from networkx import node_clique_number
from mpi4py import MPI
import numpy as np
from scipy.optimize import minimize
import scipy
import h5py
import itertools
import re
import os, copy, math
from pathlib import Path
from sys import stdout
import sys
import random
from time import time
import xml.etree.ElementTree as ET
from xml.dom import minidom
from contextlib import redirect_stderr
from io import StringIO

from .molecule import Molecule
from .veloxchemlib import mpi_master
from.veloxchemlib import hartree_in_kcalpermol, bohr_in_angstrom
from .outputstream import OutputStream
from .errorhandler import assert_msg_critical
from .solvationbuilder import SolvationBuilder
from .optimizationdriver import OptimizationDriver

# Drivers
from .scfrestdriver import ScfRestrictedDriver
from .molecularbasis import MolecularBasis
from .scfgradientdriver import ScfGradientDriver
from .scfhessiandriver import ScfHessianDriver
from .externalscfdriver import ExternalScfDriver
from .externalgradientdriver import ExternalGradientDriver
from .externalhessiandriver import ExternalHessianDriver
from .externalexcitedstatedriver import ExternalExcitedStatesScfDriver
from .externalexcitedstategradientdriver import ExternalExcitedStatesGradientDriver
from .externalexcitedstatehessiandriver import ExternalExcitedStatesHessianDriver
from .externaloptimdriver import ExternalOptimDriver
from .xtbdriver import XtbDriver
from .xtbgradientdriver import XtbGradientDriver
from .xtbhessiandriver import XtbHessianDriver
from .interpolationdriver import InterpolationDriver
from .interpolationdatapoint import InterpolationDatapoint
from .localbayesresidual import LocalBayesResidual

with redirect_stderr(StringIO()) as fg_err:
    import geometric

try:
    import openmm as mm
    import openmm.app as app
    import openmm.unit as unit
except ImportError:
    pass

from concurrent.futures import ProcessPoolExecutor

# ---- global worker state (created once per process) ----
_W = {}

def _init_worker(z_matrix, impes_dict, sym_dict, sym_datapoints, dps, idx, exponent_p_q, beta, e_x):
    """Runs once per process. Build a private driver & static constants."""
    # Avoid oversubscription when each process calls BLAS:
    os.environ["OMP_NUM_THREADS"] = "1"
    os.environ["MKL_NUM_THREADS"] = "1"
    os.environ["OPENBLAS_NUM_THREADS"] = "1"

    driver = InterpolationDriver(z_matrix)
    driver.use_symmetry = True
    driver.update_settings(impes_dict)
    driver.symmetry_information = sym_dict
    driver.qm_symmetry_data_points = sym_datapoints
    driver.distance_thrsh = 1000
    driver.exponent_p = exponent_p_q[0]
    driver.exponent_q = exponent_p_q[1]
    driver.print = False
    driver.calc_optim_trust_radius = True
    # Deep copy datapoints so each process can mutate dp.confidence_radius safely:
    driver.qm_data_points = [copy.deepcopy(dp) for dp in dps]

    _W["driver"] = driver
    _W["idx"]    = np.asarray(idx, dtype=int)
    _W["nsub"]   = len(idx)
    _W["beta"]   = float(beta)
    _W["e_x"]    = float(e_x)
    _W["conv"]   = float(hartree_in_kcalpermol())


def _eval_structure(payload):
    """Compute loss and grad contrib for one structure index."""
    s_idx, mol, E_qm, G_qm_flat, alphas = payload
    drv    = _W["driver"]
    idx    = _W["idx"]; nsub = _W["nsub"]
    beta   = _W["beta"]; e_x = _W["e_x"]; conv = _W["conv"]

    # Update alphas in this worker's private datapoints
    for j, dp in enumerate(drv.qm_data_points):
        dp.confidence_radius = alphas[j]

    # Run interpolation for this structure
    drv.compute(mol)
    E_interp = drv.get_energy()                  # Hartree
    G_interp = drv.get_gradient().reshape(-1)    # Hartree/Bohr

    # Vectorized getters (M, MxD, etc.)
    P       = np.asarray(drv.potentials, dtype=np.float64)                # (M,)
    G       = np.asarray(drv.gradients, dtype=np.float64).reshape(len(P), -1)  # (M,D)
    w_dα    = np.asarray(drv.dw_dalpha_list, dtype=np.float64)                  # (M,)
    w_x     = np.asarray(drv.dw_dX_dalpha_list, dtype=np.float64).reshape(len(P), -1) # (M,D)
    S       = float(drv.sum_of_weights)                                         # ()
    S_x     = np.asarray(drv.sum_of_weights_grad, dtype=np.float64).reshape(-1) # (D,)

    # Residuals in kcal/mol
    dE_res  = (E_interp - E_qm) * conv                                          # ()
    dg_full = (G_interp - G_qm_flat) * conv                                     # (D,)

    # ----- loss (your anisotropic force term) -----
    g_sub   = dg_full[idx]
    h_sub   = (G_qm_flat[idx]) * conv
    nh      = np.linalg.norm(h_sub)
    L_iso   = (g_sub @ g_sub) / nsub

    if nh > 1e-8:
        u    = h_sub / nh
        tau  = 1e-8
        gate = (nh*nh) / (nh*nh + tau*tau)
        gpar = (u @ g_sub) * u
        gper = g_sub - gpar
        L_e   = dE_res**2
        L_per = (gper @ gper) / nsub
        L_par = (u @ g_sub)**2
        L_F   = gate * (beta * L_per + (1.0 - beta) * L_par) + (1.0 - gate) * L_iso
    else:
        L_e = dE_res**2
        L_F = L_iso

    loss_s = 1.0 * e_x * L_e + 0.5 * (1.0 - e_x) * L_F

    # ----- gradient wrt alphas (vector length M) -----
    # dE/dα
    dE_dα = (w_dα * (P - E_interp)) * (conv / S)                                 # (M,)

    # dG/dα (matrix M x D), all in kcal/mol
    term1 = (w_x * (P - E_interp)[:, None]) * (conv / S)                         # (M,D)
    term2 = (w_dα[:, None] * (G - G_interp[None, :])) * (conv / S)               # (M,D)
    term3 = ((w_dα * (P - E_interp)) / (S**2))[:, None] * (conv * S_x[None, :])  # (M,D)
    dG_dα = term1 + term2 - term3                                                # (M,D)

    # dL/dg (full D) via subspace projection
    v = np.zeros_like(dg_full)
    if nh > 1e-8:
        dL_dg_sub = (2.0 * gate) * ((beta/nsub)*(gper) + (1.0 - beta)*(gpar)) + (2.0 * (1.0 - gate) / nsub) * g_sub
    else:
        dL_dg_sub = 2.0 * g_sub / nsub
    v[idx] = dL_dg_sub

    grad_force_part  = dG_dα @ v                                 # (M,)
    grad_energy_part = 2.0 * e_x * dE_res * dE_dα                # (M,)
    grad_s = grad_energy_part + 0.5 * (1.0 - e_x) * grad_force_part   # (M,)

    return float(loss_s), grad_s



class AlphaOptimizer:
    def __init__(self, z_matrix, impes_dict, sym_dict, sym_datapoints, dps, structure_list, qm_e, qm_g, exponent_p_q, e_x, beta=0.8, n_workers=None):
        self.z_matrix       = z_matrix
        self.impes_dict     = impes_dict
        self.sym_dict       = sym_dict
        self.sym_datapoints = sym_datapoints
        self.dps            = dps
        self.structures     = structure_list
        self.qm_e           = np.asarray(qm_e, dtype=np.float64)
        self.qm_g_flat      = np.stack([g.reshape(-1) for g in qm_g], axis=0)  # (S,D)
        self.e_x            = float(e_x)
        self.beta           = float(beta)
        self.M              = len(dps)
        self.S              = len(structure_list)
        self.idx            = np.asarray([o*3 + i for o in sym_dict[3] for i in range(3)], dtype=int)
        
        # 🔒 Force each worker to a single core
        os.environ["OMP_NUM_THREADS"] = "1"
        os.environ["MKL_NUM_THREADS"] = "1"
        os.environ["OPENBLAS_NUM_THREADS"] = "1"
        os.environ["NUMEXPR_NUM_THREADS"] = "1"

        self._pool = ProcessPoolExecutor(
            max_workers=(n_workers or os.cpu_count() or 2),
            initializer=_init_worker,
            initargs=(z_matrix, impes_dict, sym_dict, sym_datapoints, dps, self.idx, exponent_p_q, self.beta, self.e_x),
        )
        self._cache_key = None
        self._cache_grad = None

    @staticmethod
    def _key(x):
        return tuple(np.asarray(x, float).round(14))

    def fun(self, alphas):
        key = self._key(alphas)
        # Build per-structure payloads (tiny):
        payloads = [(i, self.structures[i], self.qm_e[i], self.qm_g_flat[i], np.asarray(alphas, float))
                    for i in range(self.S)]

        # Choose a decent chunksize to lower IPC overhead:
        chunksize = max(1, math.ceil(self.S / (4 * (self._pool._max_workers or 1))))

        # Map in parallel; compute BOTH loss and grad and cache the summed grad
        futs = self._pool.map(_eval_structure, payloads, chunksize=chunksize)
        sum_loss = 0.0
        sum_grad = np.zeros(self.M, dtype=np.float64)
        for loss_s, grad_s in futs:
            sum_loss += loss_s
            sum_grad += grad_s

        self._cache_key  = key
        self._cache_grad = sum_grad

        print('Current alpha and loss', sum_loss, key)
        return float(sum_loss)

    def jac(self, alphas):
        key = self._key(alphas)
        if key == self._cache_key and self._cache_grad is not None:
            return self._cache_grad
        # Fallback: compute once (will also cache)
        _ = self.fun(alphas)
        return self._cache_grad
