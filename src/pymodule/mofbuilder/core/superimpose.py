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

import itertools
from typing import List, Tuple, Union

import numpy as np


def sort_by_distance(arr: np.ndarray) -> List[Tuple[float, int]]:
    """Sort indices by distance from the first point to each point in arr.

    Args:
        arr: (N, 3) array of points.

    Returns:
        List of (distance, index) tuples sorted by ascending distance.
    """
    distances = [(np.linalg.norm(arr[0] - arr[i]), i) for i in range(len(arr))]
    distances.sort(key=lambda x: x[0])
    return distances


def match_vectors(
    arr1: np.ndarray, arr2: np.ndarray, num: int
) -> Tuple[np.ndarray, np.ndarray]:
    """Select num points from each set by distance-from-first ordering and return aligned subsets.

    Picks the num closest points to the first element in each array, then returns
    those subsets in the same distance order for use in superimposition.

    Args:
        arr1: First (N1, 3) array of points.
        arr2: Second (N2, 3) array of points.
        num: Number of points to select from each (e.g. min(6, len(arr1), len(arr2))).

    Returns:
        Tuple (closest_vectors_arr1, closest_vectors_arr2): (num, 3) arrays.
    """
    sorted_distances_arr1 = sort_by_distance(arr1)
    sorted_distances_arr2 = sort_by_distance(arr2)

    # Select the indices by distance matching in limited number

    indices_arr1 = [sorted_distances_arr1[j][1] for j in range(num)]
    indices_arr2 = [sorted_distances_arr2[j][1] for j in range(num)]

    # reorder the matching vectors# which can induce the smallest RMSD
    closest_vectors_arr1 = np.array([arr1[i] for i in indices_arr1])
    closest_vectors_arr2 = np.array([arr2[i] for i in indices_arr2])

    return closest_vectors_arr1, closest_vectors_arr2


def superimpose(
    src_arr: Union[np.ndarray, List],
    target_arr: Union[np.ndarray, List],
    min_rmsd: float = 1e6,
) -> Tuple[float, np.ndarray, np.ndarray]:
    """Find the best rotation and translation that aligns src_arr to target_arr.

    Procedure:
    - Convert inputs to numpy arrays.
    - Select up to 6 matching vectors from each set based on distance patterns
      (using match_vectors). This reduces the search space for correspondences.
    - Try all permutations of the selected vectors from arr1 and compute the
      SVD-based superposition against the selected vectors from arr2.
    - Keep the rotation/translation that yields the smallest RMSD.

    Args:
        src_arr: Source point set (N, 3).
        target_arr: Target point set (M, 3).
        min_rmsd: Initial RMSD threshold; best solution below this is kept.

    Returns:
        Tuple of (min_rmsd, best_rot, best_tran): best RMSD, 3x3 rotation matrix,
        translation vector (length 3).
    """
    # Ensure inputs are numpy arrays
    src_arr = np.asarray(src_arr)
    target_arr = np.asarray(target_arr)

    # Select up to 6 representative vectors from each array to match by distance
    m_src, m_target = match_vectors(src_arr, target_arr,
                                    min(6, len(src_arr), len(target_arr)))

    # Initialize best transformation to identity/no-translation
    best_rot, best_tran = np.eye(3), np.zeros(3)

    # Try every possible correspondence (permutation) of the selected vectors
    for perm in itertools.permutations(m_src):
        # Compute RMSD, rotation and translation for this correspondence
        rmsd, rot, tran = svd_superimpose(np.asarray(perm), m_target)

        # Keep the transform that gives the smallest RMSD
        if rmsd < min_rmsd:
            min_rmsd, best_rot, best_tran = rmsd, rot, tran

    return min_rmsd, best_rot, best_tran


def svd_superimpose(
    src_arr: Union[np.ndarray, List], target_arr: Union[np.ndarray, List]
) -> Tuple[float, np.ndarray, np.ndarray]:
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


def superimpose_rotation_only(
    arr1: Union[np.ndarray, List],
    arr2: Union[np.ndarray, List],
    min_rmsd: float = 1e6,
) -> Tuple[float, np.ndarray, np.ndarray]:
    """Find the best rotation (no translation) that aligns arr1 to arr2 by minimizing RMSD.

    Uses the same permutation search over matched subsets as superimpose, but keeps
    translation fixed (identity). Useful when only orientation matters.

    Args:
        arr1: Source point set (N, 3).
        arr2: Target point set (M, 3).
        min_rmsd: Initial RMSD threshold; best solution below this is kept.

    Returns:
        Tuple (min_rmsd, best_rot, best_tran): best RMSD, 3x3 rotation, translation (often zero).
    """
    arr1 = np.asarray(arr1)
    arr2 = np.asarray(arr2)
    m_arr1, m_arr2 = match_vectors(arr1, arr2, min(6, len(arr1), len(arr2)))
    best_rot, best_tran = np.eye(3), np.zeros(3)
    for perm in itertools.permutations(m_arr1):
        rmsd, rot, tran = svd_superimpose(np.asarray(perm), m_arr2)
        if rmsd < min_rmsd:
            min_rmsd, best_rot, best_tran = rmsd, rot, tran
            if np.allclose(np.dot(best_tran, np.zeros(3)), 1e-2):
                break

    return min_rmsd, best_rot, best_tran
