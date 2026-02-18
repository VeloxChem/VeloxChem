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
import sys

from .veloxchemlib import mpi_master
from .oneeints import compute_electric_dipole_integrals
from .outputstream import OutputStream
from .errorhandler import assert_msg_critical


class LocalizationDriver:

    def __init__(self, comm=None, ostream=None):
        """
        Initializes the localization driver.
        """

        if comm is None:
            comm = MPI.COMM_WORLD

        if ostream is None:
            if comm.Get_rank() == mpi_master():
                ostream = OutputStream(sys.stdout)
            else:
                ostream = OutputStream(None)

        # mpi information
        self.comm = comm
        self.rank = comm.Get_rank()
        self.nodes = comm.Get_size()

        # outputstream
        self.ostream = ostream

        self.max_iter = 50
        self.threshold = 1e-12

        # choose localization routine: "boys" (only boys for now)
        # TODO: implement "pm" and _compute_pipek_mezey
        self.method = 'boys'

    def compute(self, molecule, basis, mo_coefs, mo_list):
        """
        Send to the selected localization routine.
        """

        assert_msg_critical(
            self.method.lower() in ['boys'],
            'LocalizationDriver: Invalid localization method')

        mo_list = self._parse_mo_list(mo_list)

        if self.method.lower() == 'boys':
            return self._compute_foster_boys(molecule, basis, mo_coefs, mo_list)

    def _compute_foster_boys(self, molecule, basis, mo_coefs, mo_list):
        """
        Foster-Boys localization.
        """

        if self.rank == mpi_master():
            dip_mats = np.array(compute_electric_dipole_integrals(molecule, basis))
            C = mo_coefs.copy()
            C_local = C[:, mo_list].copy()
            m = C_local.shape[1]

            r = np.array([
                np.linalg.multi_dot([C_local.T, x, C_local]) for x in dip_mats
            ])

            loc_converged = False

            for iter_idx in range(self.max_iter):
                max_theta = 0.0
                for i in range(m - 1):
                    for j in range(i+1, m):
                        r_i  = r[:,i,i]
                        r_j  = r[:,j,j]
                        d_ij = r[:,i,j]

                        r_ij = r_i - r_j
                        g = 2.0 * np.dot(r_ij, d_ij)
                        h = np.linalg.norm(r_ij)**2 - 4.0 * np.linalg.norm(d_ij)**2

                        theta_opt = 0.25 * np.arctan2(g,h)
                        if abs(theta_opt) < self.threshold:
                            continue

                        cos_val, sin_val = np.cos(theta_opt), np.sin(theta_opt)
                        Rmat = np.array([[cos_val, -sin_val], [sin_val, cos_val]])

                        C_local[:,[i,j]] = np.matmul(C_local[:,[i,j]].copy(), Rmat)

                        for x in range(3):
                            r[x, [i, j], :] = np.matmul(Rmat.T, r[x, [i, j], :].copy())
                            rx = r[x].copy()
                            rx[:, [i, j]] = np.matmul(rx[:, [i, j]].copy(), Rmat)
                            r[x] = rx

                        max_theta = max(max_theta, abs(theta_opt))

                if max_theta < self.threshold:
                    loc_converged = True
                    break

            if loc_converged:
                C_new = C.copy()
                C_new[:, mo_list] = C_local
                return C_new
            else:
                self.ostream.print_warning(
                    f"Foster-Boys did not converge after {self.max_iter} " +
                    "iterations.")
                return None

        else:
            # non-master rank
            return None

    def _parse_mo_list(self, mo_list):
        """
        Returns a clean list.
        """
        if isinstance(mo_list, (list, tuple, np.ndarray)):
            out = [int(x) for x in mo_list]
        elif isinstance(mo_list, str):
            # TODO: use input parser
            s = mo_list.strip()
            if s.startswith('[') and s.endswith(']'):
                s = s[1:-1]
            parts = [p for p in s.replace(';', ',').split(',') if p.strip() != '']
            out = [int(p) for p in parts]
        assert_msg_critical(len(out) > 0, "mo_list must not be empty.")

        return out
