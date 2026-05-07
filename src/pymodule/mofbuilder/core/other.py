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

import re
from typing import Any, List, Tuple, Union

import numpy as np
import sys
try:
    from scipy.optimize import linear_sum_assignment
except ImportError:
    pass
from ...errorhandler import assert_msg_critical
from ..utils.geometry import fractional_to_cartesian


def fetch_X_atoms_ind_array(
    array: np.ndarray, column: int, X: str
) -> Tuple[List[int], np.ndarray]:
    """Return indices and rows of atoms whose label (after stripping digits) equals X.

    Args:
        array: Input array; rows are atoms, column is checked for label.
        column: Column index to check for the atom label.
        X: Label to search for (e.g. "X" for connection atoms).

    Returns:
        Tuple of (indices, subarray): indices into array and the subarray of matching rows.
    """
    ind = [
        k for k in range(len(array))
        if re.sub(r"\d", "", array[k, column]) == X
    ]
    x_array = array[ind]
    return ind, x_array


def find_pair_x_edge_fc(
    x_matrix: np.ndarray, edge_matrix: np.ndarray, sc_unit_cell: np.ndarray
) -> Tuple[np.ndarray, np.ndarray]:
    """Find optimal assignment between X atoms and edge points using Cartesian distances.

    Converts fractional coordinates to Cartesian via sc_unit_cell, builds a distance
    matrix, and uses linear sum assignment to pair each X to one edge point.

    Args:
        x_matrix: (N, 3+) array of X atom fractional coordinates (cols 2:5 used if larger).
        edge_matrix: (M, 3+) array of edge point fractional coordinates.
        sc_unit_cell: 3x3 supercell unit cell matrix for fractional_to_cartesian.

    Returns:
        Tuple (row_ind, col_ind): row indices (x) and column indices (edge) of the pairing.
    """
    assert_msg_critical(
        "scipy" in sys.modules,
        "SciPy is required for MofBuilder.")
    dist_matrix = np.zeros((len(x_matrix), len(edge_matrix)))
    x_matrix = fractional_to_cartesian(x_matrix, sc_unit_cell)
    edge_matrix = fractional_to_cartesian(edge_matrix, sc_unit_cell)
    for i in range(len(x_matrix)):
        for j in range(len(edge_matrix)):
            dist_matrix[i, j] = np.linalg.norm(x_matrix[i] - edge_matrix[j])
    row_ind, col_ind = linear_sum_assignment(dist_matrix)
    return row_ind, col_ind


def order_edge_array(
    row_ind: np.ndarray, col_ind: np.ndarray, edges_array: np.ndarray
) -> np.ndarray:
    """Reorder edge points so they follow the assignment given by row_ind and col_ind.

    Splits edges_array by columns, reorders by col_ind so that the i-th output block
    corresponds to row_ind[i], then stacks in order of row_ind.

    Args:
        row_ind: Row indices from the assignment (length N).
        col_ind: Column indices from the assignment (length N).
        edges_array: Array of edge points (one block per column in the assignment).

    Returns:
        Reordered array of edge points.
    """
    old_split = np.vsplit(edges_array, len(col_ind))
    old_order = []
    for i in range(len(col_ind)):
        old_order.append((row_ind[i], col_ind[i],
                          old_split[sorted(col_ind).index(col_ind[i])]))
    new_order = sorted(old_order, key=lambda x: x[0])
    ordered_arr = np.vstack([new_order[j][2] for j in range(len(new_order))])
    return ordered_arr


def safe_dict_copy(d: dict) -> dict:
    """Recursively copy a dict; nested dicts and lists/arrays are copied, not referenced.

    Args:
        d: Dictionary to copy (may contain dict, list, numpy.ndarray, or scalars).

    Returns:
        New dictionary with the same structure and copied values.
    """
    new_d = {}
    for k, v in d.items():
        if isinstance(v, dict):
            new_d[k] = safe_dict_copy(v)
        elif isinstance(v, np.ndarray):
            new_d[k] = v.copy()
        elif isinstance(v, list):
            new_d[k] = list(v)
        else:
            new_d[k] = v
    return new_d


def safe_copy(value: Any) -> Union[dict, list, np.ndarray, Any]:
    """Return a deep copy of value: dict (recursive), list, or ndarray; otherwise return as-is.

    Args:
        value: Value to copy (dict, list, numpy.ndarray, or other).

    Returns:
        Copy of value, or value itself if not dict/list/ndarray.
    """
    if isinstance(value, dict):
        return safe_dict_copy(value)
    elif isinstance(value, np.ndarray):
        return value.copy()
    elif isinstance(value, list):
        return list(value)
    else:
        return value
