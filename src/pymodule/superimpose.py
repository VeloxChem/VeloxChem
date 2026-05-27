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


def svd_superimpose(src_arr, target_arr):
    """Compute RMSD and rotation/translation for superimposing two point sets via SVD.

    Ref.: "Least-Squares Fitting of Two 3-D Point Sets", IEEE Trans. Pattern
    Anal. Mach. Intell., 1987, PAMI-9(5), 698-700. DOI: 10.1109/TPAMI.1987.4767965

    Args:
        src_arr: Source point set (N, 3).
        target_arr: Target point set (M, 3); N should equal M for meaningful RMSD.

    Returns:
        Tuple of (rmsd, rot_mat, trans): RMSD, 3x3 rotation matrix, translation vector.
    """

    src_arr = np.array(src_arr)
    target_arr = np.array(target_arr)

    com1 = np.sum(src_arr, axis=0) / src_arr.shape[0]
    com2 = np.sum(target_arr, axis=0) / target_arr.shape[0]

    src_arr -= com1
    target_arr -= com2

    cov_mat = np.matmul(src_arr.T, target_arr)
    U, s, Vt = np.linalg.svd(cov_mat)

    rot_mat = np.matmul(U, Vt)
    if np.linalg.det(rot_mat) < 0:
        Vt[-1, :] *= -1.0
        rot_mat = np.matmul(U, Vt)

    diff = target_arr - np.matmul(src_arr, rot_mat)
    rmsd = np.sqrt(np.sum(diff**2) / diff.shape[0])
    trans = com2 - np.dot(com1, rot_mat)

    return rmsd, rot_mat, trans

