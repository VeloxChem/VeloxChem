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
from typing import Tuple, Dict, List, Any


def unit_cell_to_cartesian_matrix(
    aL: float, bL: float, cL: float, alpha: float, beta: float, gamma: float
) -> np.ndarray:
    """Convert unit cell parameters to a 3x3 Cartesian transformation matrix.

    Args:
        aL (float): Unit cell vector length a.
        bL (float): Unit cell vector length b.
        cL (float): Unit cell vector length c.
        alpha (float): Angle (in degrees) between b and c.
        beta (float): Angle (in degrees) between a and c.
        gamma (float): Angle (in degrees) between a and b.

    Returns:
        np.ndarray: A 3x3 transformation matrix mapping fractional to Cartesian coordinates.

    Note:
        All angles must be provided in degrees.

    Example:
        >>> unit_cell_to_cartesian_matrix(10, 10, 10, 90, 90, 90)
        array([[10.,  0.,  0.],
               [ 0., 10.,  0.],
               [ 0.,  0., 10.]])
    """
    pi = np.pi
    aL, bL, cL, alpha, beta, gamma = map(float, (aL, bL, cL, alpha, beta, gamma))
    ax = aL
    ay = 0.0
    az = 0.0
    bx = bL * np.cos(gamma * pi / 180.0)
    by = bL * np.sin(gamma * pi / 180.0)
    bz = 0.0
    cx = cL * np.cos(beta * pi / 180.0)
    cy = (cL * bL * np.cos(alpha * pi / 180.0) - bx * cx) / by
    cz = (cL ** 2.0 - cx ** 2.0 - cy ** 2.0) ** 0.5
    unit_cell = np.asarray([[ax, ay, az], [bx, by, bz], [cx, cy, cz]]).T
    return unit_cell


def fractional_to_cartesian(
    fractional_coords: np.ndarray, T: np.ndarray
) -> np.ndarray:
    """Convert fractional coordinates to Cartesian coordinates.

    Args:
        fractional_coords (np.ndarray): Array of shape (n, 3) with fractional coordinates.
        T (np.ndarray): 3x3 unit cell transformation matrix (from `unit_cell_to_cartesian_matrix`).

    Returns:
        np.ndarray: Array of shape (n, 3) of Cartesian coordinates.

    Example:
        >>> fc = np.array([[0.5, 0.5, 0.5]])
        >>> T = unit_cell_to_cartesian_matrix(10, 10, 10, 90, 90, 90)
        >>> fractional_to_cartesian(fc, T)
        array([[5., 5., 5.]])
    """
    T = T.astype(float)
    fractional_coords = fractional_coords.astype(float)
    return np.dot(T, fractional_coords.T).T


def cartesian_to_fractional(
    cartesian_coords: np.ndarray, unit_cell_inv: np.ndarray
) -> np.ndarray:
    """Convert Cartesian coordinates to fractional coordinates.

    Args:
        cartesian_coords (np.ndarray): Array of shape (n, 3) with Cartesian coordinates.
        unit_cell_inv (np.ndarray): 3x3 inverse transformation matrix of the unit cell.

    Returns:
        np.ndarray: Array of shape (n, 3) of fractional coordinates.

    Example:
        >>> cc = np.array([[5., 5., 5.]])
        >>> T = unit_cell_to_cartesian_matrix(10, 10, 10, 90, 90, 90)
        >>> cartesian_to_fractional(cc, np.linalg.inv(T))
        array([[0.5, 0.5, 0.5]])
    """
    cartesian_coords = cartesian_coords.astype(float)
    unit_cell_inv = unit_cell_inv.astype(float)
    return np.dot(unit_cell_inv, cartesian_coords.T).T


def locate_min_idx(a_array: np.ndarray) -> Tuple[int, int]:
    """Locate the index of the minimum value in a 2D array.

    Args:
        a_array (np.ndarray): A 2D numpy array.

    Returns:
        Tuple[int, int]: Row and column indices of the minimum value.

    Example:
        >>> a_array = np.array([[1, 2], [3, 0]])
        >>> locate_min_idx(a_array)
        (1, 1)
    """
    idx = np.argmin(a_array)
    row_idx = idx // a_array.shape[1]
    col_idx = idx % a_array.shape[1]
    return row_idx, col_idx


def reorthogonalize_matrix(matrix: np.ndarray) -> np.ndarray:
    """Ensure the input matrix is a valid rotation matrix with determinant 1.

    Args:
        matrix (np.ndarray): Square matrix to orthogonalize.

    Returns:
        np.ndarray: Closest valid rotation matrix with determinant 1.

    Example:
        >>> mat = np.eye(3)
        >>> reorthogonalize_matrix(mat)
        array([[1., 0., 0.],
               [0., 1., 0.],
               [0., 0., 1.]])
    """
    U, _, Vt = np.linalg.svd(matrix)
    R = np.dot(U, Vt)
    if np.linalg.det(R) < 0:
        U[:, -1] *= -1
        R = np.dot(U, Vt)
    return R


def find_optimal_pairings(
    node_i_positions: np.ndarray, node_j_positions: np.ndarray
) -> List[int]:
    """Find the optimal one-to-one atom pairing between two nodes using a greedy distance approach.

    Args:
        node_i_positions (np.ndarray): Array of shape (n, 4) for node i's atoms [index, x, y, z].
        node_j_positions (np.ndarray): Array of shape (m, 4) for node j's atoms [index, x, y, z].

    Returns:
        List[int]: Indices [i, j] for the best single match (greedy, not Hungarian/global; see note).

    Note:
        This function currently pairs only the two closest atoms, not the full assignment (Hungarian method).
        For small clusters (single pair), this is sufficient; for larger clusters, a general assignment algorithm is preferred.

    Example:
        >>> pos_i = np.array([[0, 0.0, 0.0, 0.0], [1, 1.0, 0.0, 0.0]])
        >>> pos_j = np.array([[3, 0.1, 0.0, 0.0], [2, 1.0, 1.0, 0.0]])
        >>> find_optimal_pairings(pos_i, pos_j)
        [0, 0]
    """
    num_i, num_j = len(node_i_positions), len(node_j_positions)
    cost_matrix = np.zeros((num_i, num_j))
    for i in range(num_i):
        for j in range(num_j):
            cost_matrix[i, j] = np.linalg.norm(
                node_i_positions[i, 1:] - node_j_positions[j, 1:]
            )
    row_ind, col_ind = locate_min_idx(cost_matrix)
    return [row_ind, col_ind]


def find_edge_pairings(
    sorted_nodes: List[Any],
    sorted_edges: List[Tuple[int, int]],
    atom_positions: Dict[int, np.ndarray],
) -> Dict[Tuple[int, int], List[int]]:
    """Identify optimal atom pairings for each edge in a graph.

    Args:
        sorted_nodes (List[Any]): Vertices (node indices) used in the graph.
        sorted_edges (List[Tuple[int, int]]): List of graph edges as tuples of node indices.
        atom_positions (Dict[int, np.ndarray]): Mapping node index to positions array [[idx, x, y, z], ...].

    Returns:
        Dict[Tuple[int, int], List[int]]: Mapping from each edge to node-local indices in atom_positions.
            For each (i, j), value is [i_idx, j_idx] for best pair.

    Example:
        >>> edges = [(0, 1)]
        >>> atom_positions = {0: np.array([[0,0,0,0]]), 1: np.array([[0,1,0,0]])}
        >>> find_edge_pairings([0,1], edges, atom_positions)
        {(0, 1): [0, 0]}
    """
    edge_pairings: Dict[Tuple[int, int], List[int]] = {}
    for i, j in sorted_edges:
        node_i_positions = atom_positions[i]
        node_j_positions = atom_positions[j]
        pairs = find_optimal_pairings(node_i_positions, node_j_positions)
        edge_pairings[(i, j)] = pairs
    return edge_pairings


def Carte_points_generator(xyz_num: Tuple[int, int, int]) -> np.ndarray:
    """Generate a 3D grid of points with integer coordinates, inclusive in each direction.

    Args:
        xyz_num (Tuple[int, int, int]): Number of divisions in the x, y, and z directions.

    Returns:
        np.ndarray: Array of shape (n, 3) of all integer-lattice points in the grid.

    Example:
        >>> Carte_points_generator((1, 1, 1))
        array([[0, 0, 0],
               [0, 0, 1],
               [0, 1, 0],
               [0, 1, 1],
               [1, 0, 0],
               [1, 0, 1],
               [1, 1, 0],
               [1, 1, 1]])
    """
    x_num, y_num, z_num = xyz_num
    # Use meshgrid for efficient point generation
    x = np.arange(x_num + 1)
    y = np.arange(y_num + 1)
    z = np.arange(z_num + 1)

    # Create meshgrid and reshape to get all points
    xx, yy, zz = np.meshgrid(x, y, z, indexing='ij')
    points = np.vstack([xx.ravel(), yy.ravel(), zz.ravel()]).T

    return points
