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
import math
import h5py
import sys

from .veloxchemlib import hartree_in_ev, bohr_in_angstrom, mpi_master
from .oneeints import compute_electric_dipole_integrals
from .subcommunicators import SubCommunicators
from .outputstream import OutputStream
from .errorhandler import assert_msg_critical

from .inputparser import (parse_input, print_keywords, print_attributes,
                          get_random_string_parallel)

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

    def compute(self, molecule, basis, scf_res, mo_list):
        """
        Foster-Boys localization.
        """

        if self.rank == mpi_master():
            dip_mats = np.array(compute_electric_dipole_integrals(molecule, basis))
            C = scf_res['C_alpha']
            C_local = C[:, mo_list].copy()
            m = C_local.shape[1]

            r = np.array([C_local.T @ x @ C_local for x in dip_mats])

            for l in range(self.max_iter):
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
                            
                        cos, sin = np.cos(theta_opt), np.sin(theta_opt)
                        R = np.array([[cos, -sin], 
                                    [sin,  cos]])

                        C_local[:,[i,j]] = C_local[:,[i,j]].copy() @ R

                        for x in range(3):
                            r[x, [i, j], :] = R.T @ r[x, [i, j], :].copy()
                            r[x, :, [i, j]] = r[x, :, [i, j]].copy() @ R

                        max_theta = max(max_theta, abs(theta_opt))

                if max_theta < self.threshold:
                    #print(f'Total iterations: {l}')
                    break
            else:
                self.ostream.print_info(f"Foster–Boys did not converge after {self.max_iter} iterations.")
                    
            C_loc = C.copy()
            C_loc[:, mo_list] = C_local

            ortho = np.max(np.abs(C_loc.T @ scf_res['S'] @ C_loc - np.eye(C_loc.shape[1])))
            if ortho > 1e-9:
                print('[WARNING] Transformed MOs not orthonormal!')

            return C_loc
        else:
            return None


