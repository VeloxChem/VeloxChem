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

from mpi4py import MPI
import multiprocessing as mp
from collections import Counter
import os
import numpy as np
import math
import random
import scipy
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


class LocalBayesResidual():
    """
    Implements the potential energy surface interpolation driver.
    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.
    Instance variables
        - z_matrix: the Z-matrix of indices defining the internal
            coordinates.
        - impes_coordinate: the new molecular geometry at which the energy
            is to be calculated by interpolation.
        - energy: the energy determined by interpolation.
        - gradient: the Cartesian gradient determined by interpolation.
        - exponent_p: exponent required in the interpolation procedure (both
            simple and Shepard interpolation).
        - exponent_q: exponent required in the Shepard interpolation procedure.
        - confidence_radius: confidence radius required by the Shepard
            interpolation procedure.
        - imforcefield_file: File containing the necessary information to construct the interpolation forcefield
    """
    def __init__(self, z_matrix, comm=None, ostream=None):
        """
        Initializes the IMPES driver.
        """
        # if comm is None:
        #     comm = MPI.COMM_WORLD
        # if ostream is None:
        #     if comm.Get_rank() == mpi_master():
        #         ostream = OutputStream(sys.stdout)
        #     else:
        #         ostream = OutputStream(None)
        # # MPI information
        # self.comm = comm
        # self.rank = self.comm.Get_rank()
        # self.nodes = self.comm.Get_size()
        # ostream = OutputStream(sys.stdout)
        # self.ostream = ostream
        self.z_matrix = z_matrix
        m = len(z_matrix)
        self.m = m
        self.d = 1 + m + m*(m+1)//2
        self.internal_coordinates_values = None
        self.interpolation_energy = None
        self.interpolation_gradient = None
        self.symmetry_information = None
        self.symmetry_dihedral_lists = None
        self.qm_energy = None
        self.qm_gradient = None
        self.internal_coordinates = None
        self.observation_queue = []
        self.N_max = 30
        self.define_internal_coordinates()

    
    def init_bayesian(
        self,
        sigma2      = 1e-5,      # assumes ~0.01 Hartree (0.6 kcal/mol) noise 
        # sigma2 is the trust of the given observation so how much noise it is to be assumed to be in the data
        lambda_E    = 1/1e-3,    # ~1.6 kcal/mol^2 → lambda_E ≈ 384
        # lambda_E in this context it means how much I already trust the energy of the shepard interpolation framework
        lambda_g    = 1/1e-2,    # tighter: lambda_g ≈ 2500
        # lambda_g how much the gradient is allowed to be participating in the correction of the Taylor expansion: high lamda_g assumes the interpolation gradient is already good
        lambda_H    = 1/1e-2):    # even tighter: lambda_H ≈ 10,000): # ###-CHG
        # lambda_H how much the Hessian is allowed to be participating in the correction of the Taylor expansion: high lamda_H assumes the interpolation captures the curvature well

        """
        Allocate Bayesian containers that will learn the *model discrepancy*
            r(q) = E_QM − E_Shep
            s(q) = g_QM − g_Shep
        The prior mean is zero, so we insert **no initial observation**.
        """
        # ---------- dimensions ------------------------------------------------
        

        # ---------- prior precision Λ₀⁻¹ (diagonal) ---------------------------
        self.lambda_inv = np.diag(np.concatenate((
            np.full(1,    1.0 / lambda_E),   # c₀
            np.full(self.m,    1.0 / lambda_g),   # linear coeffs
            np.full(self.d-1-self.m,1.0 / lambda_H)))) # quadratic coeffs
        # ---------- containers -------------------------------------------------
        self.prec_matrix = self.lambda_inv.copy()   # A = Λ₀⁻¹   (no data yet)
        self.rhs         = np.zeros(self.d)              # b = 0
        self.sigma2      = sigma2
        # ---------- Cholesky of prior & mean ----------------------------------
        self.low_trig_chol_prec_matrix = np.linalg.cholesky(self.prec_matrix)
        self.post_mean                 = np.zeros(self.d)   # E[r(q)] = 0
    
    # def bayes_update(self, x_obs, dE, dg, w_shep=1.0):
    #     """
    #     Low-rank Bayesian update with one *residual* observation
    #         dE  = E_QM − E_Shep          (scalar)
    #         dg  = g_QM − g_Shep          (1-D array length m)
    #     x_obs  : internal coordinates of the geometry
    #     w_shep : Shepard-style distance weight (default 1.0 → no damping)
    #     """
    #     # ---------- target variance trick (unchanged) -------------------------
    #     v_target = dE**2                       # keep heuristic
    #     w_scalar = self._target_weight(x_obs, v_target)
    #     # ---------- feature rows Φ and target y --------------------------------
    #     dx_phi = self._phi(x_obs)              # shape (d,)
    #     Gmat   = self._G  (x_obs)              # shape (m,d)
    #     Phi = np.vstack((dx_phi, -Gmat))       # (1+m) × d
    #     y   = np.concatenate(([dE], dg.ravel()))
    #     # ---------- weighting W^(1/2) -----------------------------------------
    #     m       = len(self.internal_coordinates_values)
    #     W_diag  = np.full(1+m, w_scalar * w_shep)  # same weight for E & grads
    #     sqrtW   = np.sqrt(W_diag)
    #     Phi_w = sqrtW[:, None] * Phi
    #     y_w   = sqrtW * y
    #     # ---------- rank-k update of A and b -----------------------------------
    #     self.prec_matrix += Phi_w.T @ Phi_w     # A ← A + ΦᵀWΦ
    #     self.rhs         += Phi_w.T @ y_w       # b ← b + ΦᵀWy
    #     # ---------- refresh posterior ------------------------------------------
    #     self.low_trig_chol_prec_matrix = np.linalg.cholesky(self.prec_matrix)
    #     self.post_mean = scipy.linalg.cho_solve(
    #                         (self.low_trig_chol_prec_matrix, True),
    #                         self.rhs)

    # -------------- 2.  add a *new* energy/gradient row block -------------
    def bayes_update(self, x_obs, dE, dg, w_shep=1.0):
        """Low-rank update with a ghost observation or a neighbour’s point."""
        
        # dE = E_obs - E_int
        # v_target   = dE**2
        # w_diag = self._block_precisions(dE, w_shep, kappa=1.0)
        # w_scalar = self._target_weight(x_obs, v_target)
        # print('Scalar weight', w_scalar)
        # --- feature rows --------------------------------------------------
        dx_phi  = self._phi(x_obs)          # helper below
        Gmat    = self._G  (x_obs)          # (m,d)


        Phi     = np.vstack((dx_phi, -Gmat))
        y       = np.concatenate(([dE], dg.ravel()))
        # --- weighting -----------------------------------------------------
        # m          = len(self.internal_coordinates_values)
        # W_diag     = np.concatenate((
        #         np.full(1 + m, w_scalar),          # E + grads
        #         # if you have Hessian rows:
        #         # np.full(n_H, w_scalar)
        #      ))
        # sqrtW      = np.sqrt(W_diag)
        # Phi_w      = sqrtW[:, None] * Phi
        # y_w        = sqrtW * y
        
        
        # sqrt_w  = np.sqrt(w_diag)
        # Phi_w   = sqrt_w[:, None] * Phi
        # y_w     = sqrt_w * y

        w_obs   = w_shep / self.sigma2
        Phi_w   = np.sqrt(w_obs) * Phi
        y_w     = np.sqrt(w_obs) * y

        # --- rank-k downdate/update ---------------------------------------
        self.prec_matrix += Phi_w.T @ Phi_w
        self.rhs         += Phi_w.T @ y_w
        # refresh Cholesky + mean
        self.low_trig_chol_prec_matrix = np.linalg.cholesky(self.prec_matrix)
        self.post_mean = scipy.linalg.cho_solve(
                            (self.low_trig_chol_prec_matrix, True),
                            self.rhs)
        
    # -------------- 3.  predict μ and σ² at any query point ---------------
    def bayes_predict(self, x_query):
        from scipy.linalg import solve_triangular
        dx_phi = self._phi(x_query)

        mu     = dx_phi @ self.post_mean
        
        # print('in predict', dx_phi, self.post_mean, mu)

        # local variance  φᵀΛφ + σ²
        v = scipy.linalg.solve_triangular(
        self.low_trig_chol_prec_matrix, dx_phi, lower=True)
        var    = dx_phi @ v + self.sigma2
        return mu, var
    # ------------- helpers that live inside the class --------------------
    def _phi(self, x):
        """feature vector  [1, Δq, ½ vec(ΔqΔqᵀ)]   in *internal coordinates*"""
        dx   = x - self.internal_coordinates_values           # Δq

        for i, element in enumerate(self.z_matrix[self.symmetry_information[-1][1]:], start=self.symmetry_information[-1][1]): 
            if tuple(sorted(element)) in self.symmetry_information[6][3]: 
                dx[i] = 0.5 * (1.0 + np.cos(3.0 * dx[i] + np.pi))
            elif tuple(sorted(element)) in self.symmetry_information[6][2]:
                dx[i] = 0.5 * (1.0 + np.cos(2.0 * dx[i] + np.pi))
            else:
                dx[i] = 0.5 * (1.0 + np.cos(1.0 * dx[i] + np.pi)) 

        quad = np.outer(dx, dx)[np.triu_indices_from(
                                np.empty((len(dx),)*2))]       # upper tri
        return np.concatenate(([1.0], dx, 0.5*quad))
    def _G(self, x):
        """Jacobian ∂φ/∂x   size (m,d)   in internal coordinates."""
        dx   = x - self.internal_coordinates_values
        for i, element in enumerate(self.z_matrix[self.symmetry_information[-1][1]:], start=self.symmetry_information[-1][1]): 
            if tuple(sorted(element)) in self.symmetry_information[6][3]: 
                dx[i] = 0.5 * (1.0 + np.cos(3.0 * dx[i] + np.pi))
            elif tuple(sorted(element)) in self.symmetry_information[6][2]:
                dx[i] = 0.5 * (1.0 + np.cos(2.0 * dx[i] + np.pi))
            else:
                dx[i] = 0.5 * (1.0 + np.cos(1.0 * dx[i] + np.pi))
        
        m    = len(dx)
        d    = 1 + m + m*(m+1)//2
        Gmat = np.zeros((m, d))
        Gmat[:, 1:1+m] = -np.eye(m)
        col  = 1 + m
        for i in range(m):
            for j in range(i, m):
                Gmat[i, col] += -0.5 * dx[j]
                Gmat[j, col] += -0.5 * dx[i]
                col += 1
        return Gmat
    
    def bayes_update_sliding(self, x_obs, dE, dg, w_shep=1.0):
        # Store new observation
        self.observation_queue.append((x_obs, dE, dg, w_shep))
        
        # Enforce sliding window
        if len(self.observation_queue) > self.N_max:
            self.observation_queue.pop(0)
        
        # Reset Bayesian containers to prior
        self.prec_matrix = self.lambda_inv.copy()
        self.rhs = np.zeros(self.d)
        
        # Rebuild posterior from all current observations
        for (x, de, grad, w) in self.observation_queue:
            dx_phi = self._phi(x)
            Gmat = self._G(x)
            Phi = np.vstack((dx_phi, -Gmat))
            y = np.concatenate(([de], grad.ravel()))
            w_obs = w / self.sigma2
            
            Phi_w = np.sqrt(w_obs) * Phi
            y_w = np.sqrt(w_obs) * y
            
            self.prec_matrix += Phi_w.T @ Phi_w
            self.rhs += Phi_w.T @ y_w
        
        # Recompute posterior
        self.low_trig_chol_prec_matrix = np.linalg.cholesky(self.prec_matrix)
        self.post_mean = scipy.linalg.cho_solve(
                            (self.low_trig_chol_prec_matrix, True),
                            self.rhs)



    def define_internal_coordinates(self):
        """
        Defines the internal coordinates from the Z-matrix.
        """
        assert_msg_critical(self.z_matrix is not None, 'InterpolationDatapoint: No Z-matrix defined.')
        self.internal_coordinates = []
        for z in self.z_matrix:
            
            if len(z) == 2:
                q = geometric.internal.Distance(*z)
            elif len(z) == 3:
                q = geometric.internal.Angle(*z)
            elif len(z) == 4:
                q = geometric.internal.Dihedral(*z)
            else:
                assert_msg_critical(False, 'InterpolationDatapoint: Invalid entry size in Z-matrix.')
            self.internal_coordinates.append(q)

    def compute_internal_coordinates_values(self, coordinates):
        """
        Creates an array with the values of the internal coordinates
        and saves it in self.
        """

        cartesian_coordinates = coordinates

        n_atoms = cartesian_coordinates.shape[0]
        coords = cartesian_coordinates.reshape((n_atoms * 3))

        int_coords = []

        for q in self.internal_coordinates:
            if (isinstance(q, geometric.internal.Distance)):
                int_coords.append(1.0 / q.value(coords))
            else:
                int_coords.append((q.value(coords)))

        self.internal_coordinates_values = np.array(int_coords)
