import numpy as np
from scipy.optimize import minimize

###
def unit_cell_to_cartesian_matrix(aL, bL, cL, alpha, beta, gamma):
    pi = np.pi
    """Convert unit cell parameters to a Cartesian transformation matrix."""
    aL, bL, cL, alpha, beta, gamma = list(map(float, (aL, bL, cL, alpha, beta, gamma)))
    ax = aL
    ay = 0.0
    az = 0.0
    bx = bL * np.cos(gamma * pi / 180.0)
    by = bL * np.sin(gamma * pi / 180.0)
    bz = 0.0
    cx = cL * np.cos(beta * pi / 180.0)
    cy = (cL * bL * np.cos(alpha * pi / 180.0) - bx * cx) / by
    cz = (cL**2.0 - cx**2.0 - cy**2.0) ** 0.5
    unit_cell = np.asarray([[ax, ay, az], [bx, by, bz], [cx, cy, cz]]).T
    return unit_cell

def fractional_to_cartesian(fractional_coords, T):
    T = T.astype(float)
    fractional_coords = fractional_coords.astype(float)
    """Convert fractional coordinates to Cartesian using the transformation matrix."""
    return np.dot(T, fractional_coords.T).T


def cartesian_to_fractional(cartesian_coords, unit_cell_inv):
    cartesian_coords = cartesian_coords.astype(float)
    unit_cell_inv = unit_cell_inv.astype(float)
    """Convert Cartesian coordinates to fractional coordinates using the inverse transformation matrix."""
    return np.dot(unit_cell_inv, cartesian_coords.T).T

######rotation optimization 2step######
def locate_min_idx(a_array):
    # print(a_array,np.min(a_array))
    idx = np.argmin(a_array)
    row_idx = idx // a_array.shape[1]
    col_idx = idx % a_array.shape[1]
    return row_idx, col_idx


def reorthogonalize_matrix(matrix):
    """
    Ensure the matrix is a valid rotation matrix with determinant = 1.
    """
    U, _, Vt = np.linalg.svd(matrix)
    R = np.dot(U, Vt)
    if np.linalg.det(R) < 0:
        U[:, -1] *= -1
        R = np.dot(U, Vt)
    return R


def objective_function_pre(
    params, G, static_atom_positions, sorted_nodes, sorted_edges
):
    """
    Objective function to minimize distances between paired node to paired node_com along edges.

    Parameters:
        params (numpy.ndarray): Flattened array of rotation matrices.
        G (networkx.Graph): Graph structure.
        atom_positions (dict): Original positions of X atoms for each node.


    Returns:
        float: Total distance metric to minimize.
    """
    num_nodes = len(G.nodes())
    rotation_matrices = params.reshape(num_nodes, 3, 3)
    total_distance = 0.0

    for i, j in sorted_edges:
        R_i = reorthogonalize_matrix(rotation_matrices[i])

        com_i = G.nodes[sorted_nodes[i]]["ccoords"]
        com_j = G.nodes[sorted_nodes[j]]["ccoords"]
        # Rotate positions around their mass center
        rotated_i_positions = (
            np.dot(static_atom_positions[i][:, 1:] - com_i, R_i.T) + com_i
        )

        dist_matrix = np.empty((len(rotated_i_positions), 1))
        for idx_i in range(len(rotated_i_positions)):
            dist = np.linalg.norm(rotated_i_positions[idx_i] - com_j)
            dist_matrix[idx_i, 0] = dist
            # total_distance += dist ** 2
        if np.argmin(dist_matrix) > 1:
            total_distance += 1e4  # penalty for the distance difference
        else:
            total_distance += np.min(dist_matrix) ** 2
        #
        for idx_i in range(len(rotated_i_positions)):
            # second min and min distance difference not max
            if len(dist_matrix[idx_i, :]) > 1:
                second_min_dist = np.partition(dist_matrix[idx_i, :], 1)[1]
            else:
                second_min_dist = np.partition(dist_matrix[idx_i, :], 0)[0]
            diff = second_min_dist - np.min(dist_matrix[idx_i, :])

            if diff < 4:
                total_distance += 1e4

        total_distance += 1e3 / (
            np.max(dist_matrix) - np.min(dist_matrix)
        )  # reward for the distance difference

    return total_distance


def objective_function_after(
    params, G, static_atom_positions, sorted_nodes, sorted_edges
):
    """
    Objective function to minimize distances between paired atoms along edges. just use minimum distance

    Parameters:
        params (numpy.ndarray): Flattened array of rotation matrices.
        G (networkx.Graph): Graph structure.
        atom_positions (dict): Original positions of X atoms for each node.
        edge_pairings (dict): Precomputed pairings for each edge.

    Returns:
        float: Total distance metric to minimize.
    """
    num_nodes = len(G.nodes())
    rotation_matrices = params.reshape(num_nodes, 3, 3)
    total_distance = 0.0

    for i, j in sorted_edges:
        R_i = reorthogonalize_matrix(rotation_matrices[i])
        R_j = reorthogonalize_matrix(rotation_matrices[j])

        com_i = G.nodes[sorted_nodes[i]]["ccoords"]
        com_j = G.nodes[sorted_nodes[j]]["ccoords"]

        # Rotate positions around their mass center
        rotated_i_positions = (
            np.dot(static_atom_positions[i][:, 1:] - com_i, R_i.T) + com_i
        )
        rotated_j_positions = (
            np.dot(static_atom_positions[j][:, 1:] - com_j, R_j.T) + com_j
        )

        dist_matrix = np.empty((len(rotated_i_positions), len(rotated_j_positions)))
        for idx_i in range(len(rotated_i_positions)):
            for idx_j in range(len(rotated_j_positions)):
                dist = np.linalg.norm(
                    rotated_i_positions[idx_i] - rotated_j_positions[idx_j]
                )
                dist_matrix[idx_i, idx_j] = dist

        if np.argmin(dist_matrix) > 1:
            total_distance += 1e4  # penalty for the distance difference
        else:
            total_distance += np.min(dist_matrix) ** 2

        for idx_i in range(len(rotated_i_positions)):
            # second min and min distance difference not max
            if len(dist_matrix[idx_i, :]) > 1:
                second_min_dist = np.partition(dist_matrix[idx_i, :], 1)[1]
            else:
                second_min_dist = np.partition(dist_matrix[idx_i, :], 0)[0]
            diff = second_min_dist - np.min(dist_matrix[idx_i, :])
            if diff < 3:
                total_distance += 1e4
        for idx_j in range(len(rotated_j_positions)):
            # second min and min distance difference not max
            if len(dist_matrix[:, idx_j]) > 1:
                second_min_dist = np.partition(dist_matrix[:, idx_j], 1)[1]
            else:
                second_min_dist = np.partition(dist_matrix[:, idx_j], 0)[0]
            diff = second_min_dist - np.min(dist_matrix[:, idx_j])

            if diff < 3:
                total_distance += 1e4

    return total_distance


def optimize_rotations_pre(
    num_nodes,
    G,
    sorted_nodes,
    sorted_edges,
    atom_positions,
    initial_rotations,
    opt_method,
    maxfun,
    maxiter,
    disp,
    eps,
    iprint,
):
    """
    Optimize rotations for all nodes in the graph.

    Parameters:
        G (networkx.Graph): Graph structure with edges between nodes.
        atom_positions (dict): Positions of X atoms for each node.

    Returns:
        list: Optimized rotation matrices for all nodes.
    """
    print("optimize_rotations_step1")
    # initial_rotations = np.tile(np.eye(3), (num_nodes, 1)).flatten()
    # get a better initial guess, use random rotation matrix combination
    # initial_rotations  = np.array([reorthogonalize_matrix(np.random.rand(3,3)) for i in range(num_nodes)]).flatten()
    static_atom_positions = atom_positions.copy()
    # Precompute edge-specific pairings
    # edge_pairings = find_edge_pairings(sorted_edges, atom_positions)

    result = minimize(
        objective_function_pre,
        initial_rotations,
        args=(G, static_atom_positions, sorted_nodes, sorted_edges),
        method=opt_method,
        options={
            "maxfun": maxfun,
            "maxiter": maxiter,
            "disp": disp,
            "eps": eps,
            "iprint": iprint,
        },
    )

    # optimized_rotations = result.x.reshape(num_nodes, 3, 3)
    # optimized_rotations = [reorthogonalize_matrix(R) for R in optimized_rotations]

    optimized_rotations = result.x
    # optimized_rotations = [reorthogonalize_matrix(R) for R in optimized_rotations]
    ## # Print the optimized pairings after optimization
    ## print("Optimized Pairings (after optimization):")
    ## for (i, j), pairs in edge_pairings.items():
    ##     print(f"Node {i} and Node {j}:")
    ##     for idx_i, idx_j in pairs:
    ##         print(f"  node{i}_{idx_i} -- node{j}_{idx_j}")
    ## print()

    return optimized_rotations, static_atom_positions


def optimize_rotations_after(
    num_nodes,
    G,
    sorted_nodes,
    sorted_edges,
    atom_positions,
    initial_rotations,
    opt_method,
    maxfun,
    maxiter,
    disp,
    eps,
    iprint,
):
    """
    Optimize rotations for all nodes in the graph.

    Parameters:
        G (networkx.Graph): Graph structure with edges between nodes.
        atom_positions (dict): Positions of X atoms for each node.

    Returns:
        list: Optimized rotation matrices for all nodes.
    """
    print("optimize_rotations_step2")
    # get a better initial guess, use random rotation matrix combination
    # initial_rotations  = np.array([reorthogonalize_matrix(np.random.rand(3,3)) for i in range(num_nodes)]).flatten()
    static_atom_positions = atom_positions.copy()
    # Precompute edge-specific pairings
    # edge_pairings = find_edge_pairings(sorted_edges, atom_positions)

    result = minimize(
        objective_function_after,
        initial_rotations,
        args=(G, static_atom_positions, sorted_nodes, sorted_edges),
        method=opt_method,
        options={
            "maxfun": maxfun,
            "maxiter": maxiter,
            "disp": disp,
            "eps": eps,
            "iprint": iprint,
        },
    )

    optimized_rotations = result.x.reshape(num_nodes, 3, 3)
    optimized_rotations = [reorthogonalize_matrix(R) for R in optimized_rotations]

    ## # Print the optimized pairings after optimization
    ## print("Optimized Pairings (after optimization):")
    ## for (i, j), pairs in edge_pairings.items():
    ##     print(f"Node {i} and Node {j}:")
    ##     for idx_i, idx_j in pairs:
    ##         print(f"  node{i}_{idx_i} -- node{j}_{idx_j}")
    ## print()

    return optimized_rotations, static_atom_positions


def apply_rotations_to_atom_positions(
    optimized_rotations, G, sorted_nodes, atom_positions
):
    """
    Apply the optimized rotation matrices to the atom positions.

    Parameters:
        optimized_rotations (list): Optimized rotation matrices for each node.
        G (networkx.Graph): Graph structure.
        atom_positions (dict): Original positions of X atoms for each node.

    Returns:
        dict: Rotated positions for each node.
    """
    rotated_positions = {}

    for i, node in enumerate(sorted_nodes):
        # if node type is V
        # if 'DV' in G.nodes[node]['type']:
        # continue
        R = optimized_rotations[i]

        original_positions = atom_positions[i]

        com = G.nodes[node]["ccoords"]

        # Translate, rotate, and translate back to preserve the mass center
        translated_positions = original_positions - com
        rotated_translated_positions = np.dot(translated_positions, R.T)
        rotated_positions[node] = rotated_translated_positions + com

    return rotated_positions


def find_optimal_pairings(node_i_positions, node_j_positions):
    """
    Find the optimal one-to-one pairing between atoms in two nodes using the Hungarian algorithm.
    """
    num_i, num_j = len(node_i_positions), len(node_j_positions)
    cost_matrix = np.zeros((num_i, num_j))
    for i in range(num_i):
        for j in range(num_j):
            cost_matrix[i, j] = np.linalg.norm(
                node_i_positions[i, 1:] - node_j_positions[j, 1:]
            )

    # row_ind, col_ind = linear_sum_assignment(cost_matrix)
    # print(cost_matrix.shape) #DEBUG
    row_ind, col_ind = locate_min_idx(cost_matrix)
    # print(row_ind,col_ind,cost_matrix) #DEBUG

    return [row_ind, col_ind]


######cell_parameters optimization######


# scale optimizer for the cif parameters update
def scale_objective_function(
    params, old_cell_params, old_cartesian_coords, new_cartesian_coords
):
    a_new, b_new, c_new, _, _, _ = params
    a_old, b_old, c_old, alpha_old, beta_old, gamma_old = old_cell_params

    # Compute transformation matrix for the old unit cell, T is the unit cell matrix
    T_old = unit_cell_to_cartesian_matrix(
        a_old, b_old, c_old, alpha_old, beta_old, gamma_old
    )
    T_old_inv = np.linalg.inv(T_old)
    old_fractional_coords = cartesian_to_fractional(old_cartesian_coords, T_old_inv)

    # backup
    # old_fractional_coords = cartesian_to_fractional(old_cartesian_coords,T_old_inv)

    # Compute transformation matrix for the new unit cell
    T_new = unit_cell_to_cartesian_matrix(
        a_new, b_new, c_new, alpha_old, beta_old, gamma_old
    )
    T_new_inv = np.linalg.inv(T_new)

    # Convert the new Cartesian coordinates to fractional coordinate using the old unit cell

    # Recalculate fractional coordinates from updated Cartesian coordinates
    new_fractional_coords = cartesian_to_fractional(new_cartesian_coords, T_new_inv)

    # Compute difference from original fractional coordinates
    diff = new_fractional_coords - old_fractional_coords
    return np.sum(diff**2)  # Sum of squared differences


# Example usage
def optimize_cell_parameters(cell_info, original_ccoords, updated_ccoords):
    # Old cell parameters (example values)
    old_cell_params = cell_info  # [a, b, c, alpha, beta, gamma]

    # Old Cartesian coordinates of points (example values)
    old_cartesian_coords = np.vstack(
        list(original_ccoords.values())
    )  # original_ccoords

    # New Cartesian coordinates of the same points (example values)
    new_cartesian_coords = np.vstack(list(updated_ccoords.values()))  # updated_ccoords
    # Initial guess for new unit cell parameters (e.g., slightly modified cell)
    initial_params = cell_info

    # Bounds: a, b, c > 3; angles [0, 180]
    bounds = [(3, None), (3, None), (3, None)] + [(20, 180)] * 3

    # Optimize using L-BFGS-B to minimize the objective function
    result = minimize(
        scale_objective_function,
        x0=initial_params,
        args=(old_cell_params, old_cartesian_coords, new_cartesian_coords),
        method="L-BFGS-B",
        bounds=bounds,
    )

    # Extract optimized parameters
    optimized_params = np.round(result.x, 5)
    print(
        "Optimized New Cell Parameters:",
        optimized_params,
        "\nTemplate Cell Parameters:",
        cell_info,
    )

    return optimized_params


######others


def find_edge_pairings(sorted_nodes, sorted_edges, atom_positions):
    """
    Identify optimal pairings for each edge in the graph.

    Parameters:
        G (networkx.Graph): Graph structure with edges between nodes.
        atom_positions (dict): Positions of X atoms for each node.

    Returns:
        dict: Mapping of edges to optimal atom pairs.
              Example: {(0, 1): [(0, 3), (1, 2)], ...}
    """

    edge_pairings = {}

    for i, j in sorted_edges:
        node_i_positions = atom_positions[i]  # [index,x,y,z]
        node_j_positions = atom_positions[j]  # [index,x,y,z]

        # Find optimal pairings for this edge

        pairs = find_optimal_pairings(node_i_positions, node_j_positions)
        # print(sorted_nodes[i],sorted_nodes[j],pairs) #DEBUG
        edge_pairings[(i, j)] = pairs  # update_pairs(pairs,atom_positions,i,j)
        # idx_0,idx_1 = pairs[0]
        # x_idx_0 = atom_positions[i][idx_0][0]
        # x_idx_1 = atom_positions[j][idx_1][0]
    #
    # edge_pairings[(i, j)] = update_pairs(pairs,atom_positions,i,j) #but only first pair match
    # atom_positions[i] = np.delete(atom_positions[i], idx_0, axis=0)
    # atom_positions[j] = np.delete(atom_positions[j], idx_1, axis=0)

    return edge_pairings


def apply_rotations_to_Xatoms_positions(
    optimized_rotations,
    G,
    sorted_nodes,
    sorted_edges_of_sortednodeidx,
    Xatoms_positions_dict,
):
    """
    Apply the optimized rotation matrices to the atom positions.

    Parameters:
        optimized_rotations (list): Optimized rotation matrices for each node.
        G (networkx.Graph): Graph structure.
        atom_positions (dict): Original positions of X atoms for each node.

    Returns:
        dict: Rotated positions for each node.
    """
    rotated_positions = Xatoms_positions_dict.copy()

    for i, node in enumerate(sorted_nodes):
        # if node type is V
        # if 'DV' in G.nodes[node]['type']:
        # continue
        R = optimized_rotations[i]

        original_positions = rotated_positions[i][:, 1:]
        com = G.nodes[node]["ccoords"]

        # Translate, rotate, and translate back to preserve the mass center
        translated_positions = original_positions - com
        rotated_translated_positions = np.dot(translated_positions, R.T)
        rotated_positions[i][:, 1:] = rotated_translated_positions + com
    edge_pair = find_edge_pairings(
        sorted_nodes, sorted_edges_of_sortednodeidx, rotated_positions
    )
    # print("Optimized Pairings (after optimization):") #DEBUG

    optimized_pair = {}

    for (i, j), pair in edge_pair.items():
        # print(f"Node {sorted_nodes[i]} and Node {sorted_nodes[j]}:") #DEBUG
        idx_i, idx_j = pair
        # print(f"  node{sorted_nodes[i]}_{int(idx_i)} -- node{sorted_nodes[j]}_{int(idx_j)}") #DEBUG
        optimized_pair[sorted_nodes[i], sorted_nodes[j]] = (int(idx_i), int(idx_j))

    return rotated_positions, optimized_pair


# use optimized_params to update all of nodes ccoords in G, according to the fccoords
def update_ccoords_by_optimized_cell_params(G, optimized_params):
    sG = G.copy()
    a, b, c, alpha, beta, gamma = optimized_params
    T_unitcell = unit_cell_to_cartesian_matrix(a, b, c, alpha, beta, gamma)
    updated_ccoords = {}
    for n in sG.nodes():
        updated_ccoords[n] = fractional_to_cartesian(
            T_unitcell, sG.nodes[n]["fcoords"].T
        ).T
        sG.nodes[n]["ccoords"] = updated_ccoords[n]
    return sG, updated_ccoords
