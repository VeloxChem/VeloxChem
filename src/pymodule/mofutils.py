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

from pathlib import Path
import numpy as np
import networkx as nx
import itertools
import sys
import re

from .molecule import Molecule
from .errorhandler import assert_msg_critical
from .mofoptimizer import (
    optimize_rotations_pre,
    optimize_rotations_after,
    apply_rotations_to_atom_positions,
    apply_rotations_to_Xatoms_positions,
    update_ccoords_by_optimized_cell_params,
    optimize_cell_parameters,
    expand_setrots,
)


try:
    from scipy.optimize import linear_sum_assignment
except ImportError:
    pass


def nn(s):
    return re.sub(r"\d+", "", s)


def nl(s):
    return re.sub(r"\D+", "", s)


def pname(s):
    # return primitive_cell_vertex_node_name
    # extract V1 from 'V1_[-1.  0.  0.]'
    return s.split("_")[0]


def lname(s):
    # return array format of list of super node name
    # extract [-1.  0.  0.] from 'V1_[-1.  0.  0.]'
    if len(s.split("_")) < 2:
        lis = np.array([0.0, 0.0, 0.0])
    else:
        lis = np.asanyarray(s.split("_")[1][1:-1].split(), dtype=float)
    return lis


def arr_dimension(arr):
    if arr.ndim > 1:
        return 2
    else:
        return 1


def is_list_A_in_B(A, B):
    return all([np.allclose(a, b, atol=1e-9) for a, b in zip(A, B)])


######below are from fetchfile.py####################
def copy_file(old_path, new_path):
    src = Path(old_path)
    dest = Path(new_path)

    if (not dest.is_file()) or (not src.samefile(dest)):
        if not dest.parent.is_dir():
            dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_text(src.read_text())

    print(f"File copied from {old_path} to {new_path}")


def fetch_pdbfile(dir_name, keywords, nokeywords):
    candidates = []
    for pdb in Path(dir_name).rglob("*.pdb"):
        name = pdb.name
        if all(i in name for i in keywords) and all(j not in name for j in nokeywords):
            candidates.append(pdb.name)
    if len(candidates) == 0:
        raise ValueError(f"Cannot find a file including '{keywords}' ")
    elif len(candidates) == 1:
        print("found the file including", keywords)
        return candidates
    elif len(candidates) > 1:
        print("found many files including", keywords)
        return candidates


def read_mof_top_dict(data_path):
    print(data_path)
    if Path(data_path, "MOF_topology_dict").exists():
        mof_top_dict_path = str(Path(data_path, "MOF_topology_dict"))
        with open(mof_top_dict_path, "r") as f:
            lines = f.readlines()
        titles = lines[0].split()
        mofs = lines[1:]
    mof_top_dict = {}
    for mof in mofs:
        mof_name = mof.split()[0]
        if mof_name not in mof_top_dict.keys():
            mof_top_dict[mof_name] = {
                "node_connectivity": int(mof.split()[1]),
                "metal": [mof.split()[2]],
                "linker_topic": int(mof.split()[3]),
                "topology": mof.split()[-1],
            }
        else:
            mof_top_dict[mof_name]["metal"].append(mof.split()[2])
    return mof_top_dict


#####below are from _process_cifstr.py####################
def remove_blank_space(value):
    return re.sub(r"\s", "", value)


def remove_empty_lines(lines):
    newlines = []
    for i in range(len(lines)):
        if lines[i].strip() != "":
            newlines.append(lines[i])
    return newlines


def remove_bracket(value):
    value_float = float(re.sub(r"\(.*?\)", "", value))
    return value_float


def remove_tail_number(value):
    return re.sub(r"\d", "", value)


def add_quotes(value):
    return "'" + value + "'"


def remove_note_lines(lines):
    newlines = []
    for i in range(len(lines)):
        m = re.search("_", lines[i])
        if m is None:
            newlines.append(lines[i])
    return newlines


def extract_quote_lines(lines):
    newlines = []
    for i in range(len(lines)):
        if lines[i].strip()[0] == "'":
            newlines.append(lines[i])
    return newlines


def extract_xyz_lines(lines):
    newlines = []
    for i in range(len(lines)):
        if lines[i].strip()[0] != "_":
            quote_value = add_quotes(remove_blank_space(lines[i]).strip())
            newlines.append(quote_value)
        newlines = remove_empty_lines(newlines)
    return newlines


def remove_quotes(value):
    pattern = r"[\"']([^\"']+)[\"']"
    extracted_values = re.findall(pattern, value)
    return extracted_values[0]


def convert_fraction_to_decimal(expression):
    # Function to replace fraction with its decimal equivalent
    # match as object of re.search
    def replace_fraction(match):
        numerator, denominator = map(int, match.groups())
        return str(numerator / denominator)

    # Regular expression to find fractions
    fraction_pattern = r"(-?\d+)/(\d+)"
    # get fraction
    converted_expression = re.sub(fraction_pattern, replace_fraction, expression)

    return converted_expression


def extract_xyz_coefficients_and_constant(expr_str):
    # Initialize coefficients and constant
    coeffs = {"x": 0, "y": 0, "z": 0}
    constant_term = 0

    # Regular expression to match terms with coefficients and variables
    pattern = r"([+-]?\d*\.?\d*)\s*([xyz])"
    matches = re.findall(pattern, expr_str)

    # Update coefficients based on matches
    for match in matches:
        coeff = match[0]
        variable = match[1]
        if coeff == "" or coeff == "+":
            coeff = 1
        elif coeff == "-":
            coeff = -1
        else:
            coeff = float(coeff)
        coeffs[variable] += coeff

    # match if no constant at tail
    if re.search(r"[a-zA-Z]$", expr_str):
        constant_term = 0
    else:
        # extract tail constant term
        constant_match = re.search(r"([a-zA-Z])(\d*.*)$", expr_str)
        if constant_match:
            # constant_term = constant_match.group(2)
            constant_term = convert_fraction_to_decimal(constant_match.group(2))
        else:
            constant_term = 0

    xyz_coeff_array = np.array([coeffs["x"], coeffs["y"], coeffs["z"]])

    return xyz_coeff_array, constant_term


#########below are from readcif_pdb.py#####################
def find_keyword(keyword, s):
    m = re.search(keyword, s)
    if m:
        return True
    else:
        return False


def read_pdb(pdbfile):
    if not Path(pdbfile).exists():
        raise FileNotFoundError(f"pdb file {pdbfile} not found")
    print(f"trying to read pdb file {pdbfile}")

    inputfile = str(pdbfile)
    with open(inputfile, "r") as fp:
        lines = fp.readlines()
    data = []
    for line in lines:
        line = line.strip()
        if len(line) > 0:  # skip blank line
            if line[0:4] == "ATOM" or line[0:6] == "HETATM":
                value_atom = line[12:16].strip()  # atom_label
                value_x = float(line[30:38])  # x
                value_y = float(line[38:46])  # y
                value_z = float(line[46:54])  # z
                atom_type = line[67:80].strip()  # atom_note
                data.append([value_atom, atom_type, value_x, value_y, value_z])
    return np.vstack(data)


def read_cif(cif_file):
    with open(cif_file, "r") as f:
        lines = f.readlines()
        nonempty_lines = [line for line in lines if line.strip()]
    # nonempty_lines=lines
    keyword1 = "loop_"
    keyword2 = r"x,\s*y,\s*z"
    keyword3 = "-x"
    # find the symmetry sector begin with x,y,z, beteen can have space or tab and comma,but just x start, not '-x'
    # keyword2 = "x,\s*y,\s*z"

    loop_key = []
    loop_key.append(0)
    linenumber = 0
    for i in range(len(nonempty_lines)):  # search for keywords and get linenumber
        # m is find keywor1 or (find keyword2 without keyword3)
        m = find_keyword(keyword1, nonempty_lines[i]) or (
            find_keyword(keyword2, nonempty_lines[i])
            and (not find_keyword(keyword3, nonempty_lines[i]))
        )

        a = re.search("_cell_length_a", nonempty_lines[i])
        b = re.search("_cell_length_b", nonempty_lines[i])
        c = re.search("_cell_length_c", nonempty_lines[i])
        alpha = re.search("_cell_angle_alpha", nonempty_lines[i])
        beta = re.search("_cell_angle_beta", nonempty_lines[i])
        gamma = re.search("_cell_angle_gamma", nonempty_lines[i])

        if m:
            loop_key.append(linenumber)
        # if not nonempty_lines[i].strip():
        #    loop_key.append(linenumber)

        else:
            if a:
                cell_length_a = remove_bracket(nonempty_lines[i].split()[1])
            elif b:
                cell_length_b = remove_bracket(nonempty_lines[i].split()[1])
            elif c:
                cell_length_c = remove_bracket(nonempty_lines[i].split()[1])
            elif alpha:
                cell_angle_alpha = remove_bracket(nonempty_lines[i].split()[1])
            elif beta:
                cell_angle_beta = remove_bracket(nonempty_lines[i].split()[1])
            elif gamma:
                cell_angle_gamma = remove_bracket(nonempty_lines[i].split()[1])

        linenumber += 1
    loop_key.append(len(nonempty_lines))
    list(set(loop_key))
    loop_key.sort()

    cell_info = [
        cell_length_a,
        cell_length_b,
        cell_length_c,
        cell_angle_alpha,
        cell_angle_beta,
        cell_angle_gamma,
    ]

    # find symmetry sectors and atom_site_sectors
    cif_sectors = []
    for i in range(len(loop_key) - 1):
        cif_sectors.append(nonempty_lines[loop_key[i] : loop_key[i + 1]])
    for i in range(len(cif_sectors)):  # find '\s*x,\s*y,\s*z' symmetry sector
        if re.search(keyword2, cif_sectors[i][0]):
            symmetry_sector = cif_sectors[i]

        if len(cif_sectors[i]) > 1:
            if re.search(r"_atom_site_label\s+", cif_sectors[i][1]):  # line0 is _loop
                atom_site_sector = cif_sectors[i]

    return cell_info, symmetry_sector, atom_site_sector


def extract_symmetry_operation_from_lines(symmetry_sector):
    symmetry_operation = []
    for i in range(len(symmetry_sector)):
        # Regular expression to match terms with coefficients and variables
        pattern = r"([+-]?\d*\.?\d*)\s*([xyz])"  # at least find a x/-x/y/-y/z/-z
        match = re.search(pattern, symmetry_sector[i])
        if match:
            string = remove_quotes(symmetry_sector[i].strip("\n"))
            no_space_string = string.replace(" ", "")
            symmetry_operation.append(no_space_string)
    if len(symmetry_operation) < 2:
        print(" no symmetry operation")
        symmetry_operation = ["x,y,z"]
    else:
        print(f"apply {len(symmetry_operation)}  symmetry operation")

    return symmetry_operation


def limit_value_0_1(new_array_metal_xyz):
    # use np.mod to limit the value in [0,1]
    new_array_metal_xyz = np.mod(new_array_metal_xyz, 1)
    return new_array_metal_xyz


def apply_sym_operator(symmetry_operations, array_metal_xyz):
    cell_array_metal_xyz = np.empty((0, 3))
    for sym_line in range(len(symmetry_operations)):
        sym_operation = np.empty((3, 3))
        constants_xyz = np.empty((1, 3))
        for i in range(3):  # x,y,z columns for operation
            coeff_xyz, constant_term = extract_xyz_coefficients_and_constant(
                symmetry_operations[sym_line].split(",")[i]
            )
            # print(f"symmetry_operations[sym_line],{symmetry_operations[sym_line]}\ncoeff_xyz, constant_term,{coeff_xyz, constant_term}")
            sym_operation[:, i] = coeff_xyz
            constants_xyz[:, i] = constant_term
        new_xyz = np.dot(array_metal_xyz, sym_operation) + constants_xyz
        cell_array_metal_xyz = np.vstack((cell_array_metal_xyz, new_xyz))

    round_cell_array_metal_xyz = np.round(limit_value_0_1(cell_array_metal_xyz), 3)
    _, unique_indices = np.unique(round_cell_array_metal_xyz, axis=0, return_index=True)
    unique_indices.sort()
    unique_metal_array = round_cell_array_metal_xyz[unique_indices]

    return unique_metal_array, unique_indices


def process_node_pdb(pdbfile, com_target_type):
    data = read_pdb(pdbfile)
    node_atoms = data[:, 0:2]
    node_ccoords = data[:, 2:5]
    node_ccoords = node_ccoords.astype(float)
    com_type_indices = [
        i for i in range(len(node_atoms)) if nn(node_atoms[i, 0]) == com_target_type
    ]
    x_indices = [j for j in range(len(node_atoms)) if nn(node_atoms[j, 0]) == "X"]
    node_x_ccoords = data[x_indices, 2:5]
    node_x_ccoords = node_x_ccoords.astype(float)
    com_type_ccoords = node_ccoords[com_type_indices]
    com_type = np.mean(com_type_ccoords, axis=0)
    node_ccoords = node_ccoords - com_type
    node_x_ccoords = node_x_ccoords - com_type
    return node_atoms, node_ccoords, node_x_ccoords


def extract_atoms_fcoords_from_lines(atom_site_sector):
    atom_site_lines = []
    keyword = "_"
    for line in atom_site_sector:  # search for keywords and get linenumber
        m = re.search(keyword, line)
        if m is None:
            atom_site_lines.append(line)

    array_atom = np.zeros(
        (len(atom_site_lines), 2), dtype=object
    )  # modified to 2 from 1, NOTE:
    array_xyz = np.zeros((len(atom_site_lines), 3))

    for i in range(len(atom_site_lines)):
        for j in [0, 1, 2, 3, 4]:  # NOTE: modified to 0-4 from 0 2 3 4
            if j < 2:
                array_atom[i, j] = remove_tail_number(atom_site_lines[i].split()[j])
            else:
                array_xyz[i, (j - 2)] = remove_bracket(atom_site_lines[i].split()[j])
    # print(f"array_atom{array_atom}") #DEBUG
    return array_atom, array_xyz


def extract_type_atoms_fcoords_in_primitive_cell(cif_file, target_type):
    cell_info, symmetry_sector, atom_site_sector = read_cif(cif_file)
    array_atom, array_xyz = extract_atoms_fcoords_from_lines(atom_site_sector)

    if len(symmetry_sector) > 1:  # need to apply symmetry operations
        # print(f"apply {len(symmetry_sector)} symmetry operation")
        array_metal_xyz = array_xyz[array_atom[:, 0] == target_type]
        array_metal_xyz = np.round(array_metal_xyz, 3)
        symmetry_sector_neat = extract_quote_lines(symmetry_sector)
        if len(symmetry_sector_neat) < 2:  # if no quote, then find x,y,z
            symmetry_sector_neat = extract_xyz_lines(symmetry_sector)
        symmetry_operations = extract_symmetry_operation_from_lines(
            symmetry_sector_neat
        )
        no_sym_array_metal_xyz, no_sym_indices = apply_sym_operator(
            symmetry_operations, array_metal_xyz
        )
        array_metal_xyz_final = no_sym_array_metal_xyz

    else:  # P1
        print("P1 cell")
        array_metal_xyz = array_xyz[array_atom[:, 0] == target_type]

        array_metal_xyz_final = np.round(array_metal_xyz, 3)
    return cell_info, array_xyz, array_metal_xyz_final


###########below are for prepare class####################
#########below are from add_dummy2node.py####################


# process node pdb file
def nodepdb(filename):
    xyzlines = ""
    inputfile = str(filename)
    with open(inputfile, "r") as fp:
        content = fp.readlines()
        # linesnumber = len(content)
    data = []
    for line in content:
        line = line.strip()
        if len(line) > 0:  # skip blank line
            if line[0:6] == "ATOM" or line[0:6] == "HETATM":
                value_atom = line[12:16].strip()  # atom_label
                value_x = float(line[30:38])  # x
                value_y = float(line[38:46])  # y
                value_z = float(line[46:54])  # z
                atom_type = line[67:80].strip()  # atom_note
                formatted_xyz_line = "{:6s}  {:8.3f} {:8.3f} {:8.3f}".format(
                    atom_type, value_x, value_y, value_z
                )
                xyzlines += formatted_xyz_line + "\n"
                data.append([value_atom, atom_type, value_x, value_y, value_z])
    head = str(len(data)) + "\n" + "\n"
    xyz_lines = head + xyzlines
    return xyz_lines, np.vstack(data)


# make node Graph based on connectivity matrix
def nodepdb2G(pdbfile, metal):
    xyzlines, node_data = nodepdb(pdbfile)
    G = nx.Graph()
    for n in node_data:
        G.add_node(
            n[0],
            ccoords=(np.asarray([float(n[2]), float(n[3]), float(n[4])])),
            type=n[1],
        )

    node_mol = Molecule.read_xyz_string(xyzlines)
    con_matrix = node_mol.get_connectivity_matrix()
    for i in range(len(node_data)):
        for j in range(i, len(node_data)):
            if con_matrix[i, j] == 1:
                if (
                    nn(node_data[i][0]) == nn(node_data[j][0])
                    and nn(node_data[i][0]) == metal
                ):  # skip same atom type bond metal-metal
                    continue
                else:
                    G.add_edge(node_data[i][0], node_data[j][0])
    return G


def fetch_template(metal):
    if metal == "Zr":  # Zr-D = 1A
        template_Zr_D = np.array(
            (
                [0.70710678, 0.70710678, 0.0],
                [-0.70710678, 0.70710678, 0.0],
                [-0.70710678, -0.70710678, 0.0],
                [0.70710678, -0.70710678, 0.0],
                [0.0, 0.70710678, 0.70710678],
                [0.0, -0.70710678, 0.70710678],
                [0.0, -0.70710678, -0.70710678],
                [0.0, 0.70710678, -0.70710678],
            )
        )
        return template_Zr_D
    elif metal == "Hf":  # Hf-D = 1A
        template_Hf_D = np.array(
            (
                [0.70710678, 0.70710678, 0.0],
                [-0.70710678, 0.70710678, 0.0],
                [-0.70710678, -0.70710678, 0.0],
                [0.70710678, -0.70710678, 0.0],
                [0.0, 0.70710678, 0.70710678],
                [0.0, -0.70710678, 0.70710678],
                [0.0, -0.70710678, -0.70710678],
                [0.0, 0.70710678, -0.70710678],
            )
        )
        return template_Hf_D
    elif metal == "Al":  # Al-D = 0.74A
        template_Al_D = np.array(
            (
                [0.74, 0.0, 0.0],
                [0.0, 0.74, 0.0],
                [0.0, 0.0, 0.74],
                [-0.74, 0.0, 0.0],
                [0.0, -0.74, 0.0],
                [0.0, 0.0, -0.74],
            )
        )
        return template_Al_D
    elif metal == "Fe":  # Fe-D = 0.86A
        template_Fe_D = np.array(
            (
                [0.86, 0.0, 0.0],
                [0.0, 0.86, 0.0],
                [0.0, 0.0, 0.86],
                [-0.86, 0.0, 0.0],
                [0.0, -0.86, 0.0],
                [0.0, 0.0, -0.86],
            )
        )
        return template_Fe_D
    elif metal == "Cr":  # Cr-D = 0.82A
        template_Cr_D = np.array(
            (
                [0.82, 0.0, 0.0],
                [0.0, 0.82, 0.0],
                [0.0, 0.0, 0.82],
                [-0.82, 0.0, 0.0],
                [0.0, -0.82, 0.0],
                [0.0, 0.0, -0.82],
            )
        )
        return template_Cr_D


def order_ccoords(d_ccoords, template_Zr_D, target_metal_coords):
    d_ccoords = d_ccoords - target_metal_coords  # centering to origin

    # template_center is already at origin
    min_dist, rot, trans = superimpose(template_Zr_D, d_ccoords)
    ordered_ccoords = np.dot(template_Zr_D, rot) + target_metal_coords
    return ordered_ccoords


def lines_of_atoms(subgraph, subgraph_nodes):
    count = 1
    rows = []
    for sn in subgraph_nodes:
        label = subgraph.nodes[sn]["type"]
        coord = subgraph.nodes[sn]["ccoords"]
        name = sn  # label+str(count)
        count += 1
        rows.append([name, label, coord[0], coord[1], coord[2]])

    return rows


def create_pdb(filename, lines):
    newpdb = []
    newpdb.append("Dummy node \n" + str(filename) + "   GENERATED BY MOF_builder\n")
    # check if the file directory exists and create it if it doesn't
    Path(str(filename) + ".pdb").parent.mkdir(parents=True, exist_ok=True)

    with open(str(filename) + ".pdb", "w") as fp:
        # Iterate over each line in the input file
        for i in range(len(lines)):
            # Split the line into individual values (assuming they are separated by spaces)
            values = lines[i]
            # Extract values based on their positions in the format string
            value1 = "ATOM"
            value2 = int(i + 1)
            value3 = values[0]  # label
            value4 = "MOL"  # residue
            value5 = 1  # residue number
            value6 = float(values[2])  # x
            value7 = float(values[3])  # y
            value8 = float(values[4])  # z
            value9 = "1.00"
            value10 = "0.00"
            value11 = values[1]  # note
            # Format the values using the specified format string
            formatted_line = "%-6s%5d%5s%4s%10d%8.3f%8.3f%8.3f%6s%6s%4s" % (
                value1,
                value2,
                value3,
                value4,
                value5,
                value6,
                value7,
                value8,
                value9,
                value10,
                value11,
            )
            newpdb.append(formatted_line + "\n")
        fp.writelines(newpdb)


def dummynode_get_bonds_from_subgraph(subgraph):
    bonds = []
    for e in list(subgraph.edges()):
        atom1 = e[0]
        atom2 = e[1]
        length = 1  # subgraph.edges[e]['weight']
        sym = "."
        if nn(atom1) == "X" or nn(atom2[0]) == "X":
            bond_type = "A"
        else:
            bond_type = "S"
        bonds.append([atom1, atom2, length, sym, bond_type])

    return bonds


def add_dummy_atoms_nodepdb(pdbfile, metal, nodeG):
    if metal in ["Zr", "Hf"]:
        metal_valence = 4
    elif metal in ["Al", "Fe", "Cr"]:
        metal_valence = 3

    dummy_pdbfile = pdbfile.removesuffix(".pdb") + "_dummy.pdb"

    # dummy_ciffile = 'test12zr_dummy.cif'
    template_metal_D = fetch_template(metal)
    sG = nodeG.copy()

    ind_max = max([int(nl(n)) for n in list(sG.nodes())])
    hydrogen_nodes = [n for n in list(sG.nodes()) if nn(n) == "H"]
    oxygen_nodes = [n for n in list(sG.nodes()) if nn(n) == "O"]
    metal_nodes = [n for n in list(sG.nodes()) if nn(n) == metal]
    count = ind_max + 1
    for mn in metal_nodes:
        neighbor_nodes = sG.adj[mn].copy()
        Ocheck = all(nn(i) == "O" for i in neighbor_nodes)
        if len(neighbor_nodes) == 2 * metal_valence and Ocheck:
            # add dummy
            beginning_cc = sG.nodes[mn]["ccoords"]
            d_ccoords = []
            for nO in neighbor_nodes:
                sO = sG.nodes[nO]["ccoords"]
                cnorm_vec = (sO - beginning_cc) / np.linalg.norm(
                    sO - beginning_cc
                )  # 1 angstrom
                d_ccoord = beginning_cc + cnorm_vec
                d_ccoords.append(d_ccoord)

                sG.remove_edge(mn, nO)
            # order ccords based on template order
            ordered_ccoords = order_ccoords(d_ccoords, template_metal_D, beginning_cc)

            for row in range(len(d_ccoords)):
                sG.add_node("D" + str(count), type="D", ccoords=ordered_ccoords[row])
                sG.add_edge(mn, "D" + str(count))
                count += 1
    # nx.draw_networkx(sG,**options)
    for hn in hydrogen_nodes:
        neighbor_nodes = list(nx.neighbors(sG, hn))
        if len(neighbor_nodes) < 1:
            # find a nearest oxygen to add a connection
            min_dist = 1000
            for on in oxygen_nodes:
                dist = np.linalg.norm(sG.nodes[hn]["ccoords"] - sG.nodes[on]["ccoords"])
                if dist < min_dist:
                    min_dist = dist
                    nearest_o = on
            sG.add_edge(hn, nearest_o)

    sG_subparts = [
        c for c in sorted(nx.connected_components(sG), key=len, reverse=True)
    ]

    head = []
    tail = []
    for sub in sG_subparts:
        l = [nn(i) for i in sub]
        if "X" not in l:
            head.append(sorted(sub))  # include D+metal
        else:
            tail.append(sorted(sub))

    sub_headlens = [len(i) for i in head]
    sub_taillens = [len(i) for i in tail]
    # sum(sub_headlens)+sum(sub_taillens),sub_headlens,sub_taillens

    dummy_count = 0
    hho_count = 0
    ho_count = 0
    o_count = 0
    ooc_count = 0

    for i in sub_headlens:
        if i == 3:
            hho_count += 1
        elif i == 2:
            ho_count += 1
        elif i == 1:
            o_count += 1
        else:
            dummy_count += 1
            dummy_res_len = i
    for j in sub_taillens:
        if j == 3:
            ooc_count += 1

    node_split_dict = {}
    node_split_dict["HHO_count"] = hho_count
    node_split_dict["HO_count"] = ho_count
    node_split_dict["O_count"] = o_count
    node_split_dict["METAL_count"] = dummy_count
    node_split_dict["dummy_res_len"] = dummy_res_len
    node_split_dict["OOC_count"] = ooc_count
    node_split_dict["inres_count"] = hho_count + ho_count + o_count + dummy_count

    node_dict_name = dummy_pdbfile.removesuffix(".pdb") + "_dict"
    node_dict_path = node_dict_name
    with open(node_dict_path, "w") as f:
        for key in list(node_split_dict):
            f.write("{:20} {:<4}".format(key, node_split_dict[key]))
            f.write("\n")
    print(node_dict_path, "saved")
    subpart_nodes = head + tail
    all_lines = []
    all_bonds = []

    for count_i in range(len(subpart_nodes)):
        subnodes = subpart_nodes[count_i]
        subgraph = nx.subgraph(sG, subnodes)
        all_lines += lines_of_atoms(subgraph, sorted(subnodes))
        all_bonds += dummynode_get_bonds_from_subgraph(subgraph)

    # create_cif(all_lines,all_bonds,os.path.dirname(dummy_ciffile),os.path.basename(dummy_ciffile)) #if want dummy node cif file
    create_pdb(dummy_pdbfile.removesuffix(".pdb"), all_lines)
    print(dummy_pdbfile.removesuffix(".pdb") + ".pdb", "saved")
    return all_lines, dummy_pdbfile


#############below are from frag_recognizer.py####################


def get_atom_name_in_subgraph(subgraph, n_id, Xs_indices):
    for ind, value in enumerate(list(subgraph.nodes)):
        if value == n_id:
            if value not in Xs_indices:
                return subgraph.nodes[n_id]["label"] + str(ind + 1)
            else:
                return "X" + str(ind + 1)


def get_bonds_from_subgraph(subgraph, Xs_indices):
    bonds = []
    for e in list(subgraph.edges):
        atom1 = get_atom_name_in_subgraph(subgraph, e[0], Xs_indices)
        atom2 = get_atom_name_in_subgraph(subgraph, e[1], Xs_indices)
        length = subgraph.edges[e]["weight"] / 50  # 50 50 50 box
        sym = "."
        if atom1[0] == "X" or atom2[0] == "X":
            bond_type = "A"
        else:
            bond_type = "S"
        bonds.append([atom1, atom2, length, sym, bond_type])

    return bonds


def create_lG(molecule):
    matrix = molecule.get_connectivity_matrix()
    coords = molecule.get_coordinates_in_angstrom()
    labels = molecule.get_labels()
    dist_matrix = molecule.get_distance_matrix_in_angstrom()
    mass_center_bohr = molecule.center_of_mass_in_bohr()
    bohr_to_angstrom = 0.529177
    mass_center_angstrom = np.asarray(mass_center_bohr) * bohr_to_angstrom
    coords = coords - mass_center_angstrom

    metal_elements_list = [
        "Ag",
        "Al",
        "Au",
        "Ba",
        "Be",
        "Bi",
        "Ca",
        "Cd",
        "Ce",
        "Co",
        "Cr",
        "Cs",
        "Cu",
        "Fe",
        "Ga",
        "Gd",
        "Hf",
        "Hg",
        "In",
        "Ir",
        "K",
        "Li",
        "Mg",
        "Mn",
        "Na",
        "Ni",
        "Pb",
        "Pd",
        "Pt",
        "Rb",
        "Rh",
        "Sc",
        "Sn",
        "Sr",
        "Ti",
        "V",
        "W",
        "Y",
        "Zn",
        "Zr",
    ]

    lG = nx.Graph()
    metals = []
    for i in range(len(labels)):
        lG.add_nodes_from([(i, {"label": labels[i], "coords": coords[i]})])
        if labels[i] in metal_elements_list:
            metals.append(i)

    i = None
    for i in range(len(labels)):
        for j in np.where(matrix[i] == 1)[0]:
            if i not in metals and j not in metals:
                lG.add_edge(i, j, weight=dist_matrix[i, j])
    return lG, metals, mass_center_angstrom


def find_center_cycle_nodes(lG):
    # To find center cycle
    target_nodes = set(nx.center(lG))
    cycles = list(nx.simple_cycles(lG, length_bound=80))
    for cycle in cycles:
        if target_nodes < set(cycle):
            return cycle


def distinguish_G_centers(lG):
    centers = nx.barycenter(lG)
    if len(centers) == 1:
        print("center is a point")
        center_class = "onepoint"
        center_nodes = centers
    elif len(centers) == 2:
        if nx.shortest_path_length(lG, centers[0], centers[1]) == 1:
            print("center is two points")
            center_class = "twopoints"
            center_nodes = centers
        else:
            print("center is a cycle")
            center_class = "cycle"
            center_nodes = find_center_cycle_nodes(lG)
    else:
        print("center is a cycle")
        center_class = "cycle"
        center_nodes = find_center_cycle_nodes(lG)
    return center_class, center_nodes


def find_center_nodes_pair(lG, center_nodes):
    if len(center_nodes) > 6:
        centers = nx.center(lG)

    pairs = []
    for i in range(len(centers)):
        for j in range(i, len(centers)):
            l = nx.shortest_path_length(lG, centers[i], centers[j])
            if l == 1:
                pairs.append([centers[i], centers[j]])

    # loop each pair to find center pair
    for p in pairs:
        a = p[0]
        b = p[1]
        ds = []
        for n in centers:
            if n not in p:
                d = nx.shortest_path_length(lG, a, n)
                ds.append(d)
        if len(set(ds)) < len(ds):
            center_pair = p

    return center_pair


def classify_nodes(lG, center_nodes):
    # Step 1: Identify the center node(s) of the graph
    # center_nodes = nx.center(lG)

    for center_ind in range(len(center_nodes)):
        # to classify which node belonging to which center_node, add 'cnode_l' to lG nodes
        # Compute the shortest path length from the center node to all other nodes
        center = center_nodes[center_ind]
        lengths = nx.single_source_shortest_path_length(lG, center)
        if center_ind == 0:
            for k in lengths:
                lG.nodes[k]["cnodes_l"] = (center, lengths[k])
        elif center_ind > 0:
            for k in lengths:
                if lengths[k] < lG.nodes[k]["cnodes_l"][1]:
                    lG.nodes[k]["cnodes_l"] = (center, lengths[k])
                elif lengths[k] == lG.nodes[k]["cnodes_l"][1]:
                    lG.nodes[k]["cnodes_l"] = (
                        -1,
                        lengths[k],
                    )  # if the node is between any two center nodes
    return lG


def get_pairX_outer_frag(connected_pairXs, outer_frag_nodes):
    for x in list(connected_pairXs):
        pairXs = [connected_pairXs[x][1], connected_pairXs[x][3]]
        if set(pairXs) < set(outer_frag_nodes):
            break
    return pairXs


def cleave_outer_frag_subgraph(lG, pairXs, outer_frag_nodes):
    subgraph_outer_frag = lG.subgraph(outer_frag_nodes)
    kick_nodes = []
    for i in list(outer_frag_nodes):
        if nx.shortest_path_length(
            subgraph_outer_frag, pairXs[0], i
        ) > nx.shortest_path_length(subgraph_outer_frag, pairXs[0], pairXs[1]):
            kick_nodes.append(i)
        elif nx.shortest_path_length(
            subgraph_outer_frag, pairXs[1], i
        ) > nx.shortest_path_length(subgraph_outer_frag, pairXs[0], pairXs[1]):
            kick_nodes.append(i)

    subgraph_single_frag = lG.subgraph(outer_frag_nodes - set(kick_nodes))
    return subgraph_single_frag


def lines_of_center_frag(
    subgraph_center_frag, Xs_indices, metals, labels, coords, mass_center_angstrom
):
    count = 1
    lines = []
    Xs = []
    for cn in list(subgraph_center_frag.nodes):
        label = subgraph_center_frag.nodes[cn]["label"]
        coord = subgraph_center_frag.nodes[cn]["coords"]
        if cn not in Xs_indices:
            name = label + str(count)
        else:
            name = "X" + str(count)
            Xs.append(count - 1)
        count += 1
        lines.append([name, label, coord[0], coord[1], coord[2]])
    for cm in metals:
        label = labels[cm]
        coord = coords[cm] - mass_center_angstrom
        name = label + str(count)
        lines.append([name, label, coord[0], coord[1], coord[2]])
    return lines, Xs


def lines_of_single_frag(subgraph_single_frag, Xs_indices):
    count = 1
    rows = []
    Xs = []
    for sn in list(subgraph_single_frag.nodes):
        label = subgraph_single_frag.nodes[sn]["label"]
        coord = subgraph_single_frag.nodes[sn]["coords"]
        if sn not in Xs_indices:
            name = label + str(count)
        else:
            name = "X" + str(count)
            Xs.append(count - 1)
        count += 1
        rows.append([name, label, coord[0], coord[1], coord[2]])
    return rows, Xs


def process_linker_molecule(
    molecule, linker_topic, save_nodes_dir="nodes", save_edges_dir="edges"
):
    coords = molecule.get_coordinates_in_angstrom()
    labels = molecule.get_labels()
    # To remove center metals
    lG, metals, mass_center_angstrom = create_lG(molecule)
    lG.remove_nodes_from(metals)
    center_class, center_nodes = distinguish_G_centers(lG)
    if linker_topic == 2 and len(center_nodes) > 6:
        center_nodes = find_center_nodes_pair(lG, center_nodes)

    lG = classify_nodes(lG, center_nodes)
    print(center_nodes)

    if center_class == "cycle" and linker_topic > 2:
        print("tritopic/tetratopic/multitopic: center is a cycle")
        connected_pairXs = {}
        Xs_indices = []
        innerX_coords = []
        for k in range(len(center_nodes)):
            linker_C_l = []
            l_list = []
            for n in lG.nodes:
                if (
                    lG.nodes[n]["cnodes_l"][0] == center_nodes[k]
                    and lG.nodes[n]["label"] == "C"
                ):
                    linker_C_l.append((n, lG.nodes[n]["cnodes_l"]))
                    l_list.append(lG.nodes[n]["cnodes_l"][1])
            center_connected_C_ind = [
                ind for ind, value in enumerate(l_list) if value == 1
            ]
            outer_connected_C_ind = [
                ind for ind, value in enumerate(l_list) if value == (max(l_list) - 1)
            ]
            if len(center_connected_C_ind) == 1 and len(outer_connected_C_ind) == 1:
                inner_X = linker_C_l[center_connected_C_ind[0]]
                outer_X = linker_C_l[outer_connected_C_ind[0]]
                if center_nodes[k] not in [inner_X[0], outer_X[0]]:
                    print(
                        "find connected X in edge frag",
                        inner_X[0],
                        outer_X[0],
                        center_nodes[k],
                    )
                    lG.remove_edge(inner_X[0], center_nodes[k])
                    connected_pairXs[center_nodes[k]] = (
                        "inner_X",
                        inner_X[0],
                        "outer_X",
                        outer_X[0],
                    )
                    Xs_indices += [center_nodes[k], inner_X[0], outer_X[0]]
                    innerX_coords.append(lG.nodes[inner_X[0]]["coords"])

        if (
            nx.number_connected_components(lG) != linker_topic + 1
        ):  # for check linker_topics+1
            print("wrong fragments")
            raise ValueError

    elif linker_topic == 2:
        if center_class == "twopoints":
            print("ditopic linker: center are two points")
            Xs_indices = []
            for k in range(len(center_nodes)):
                linker_C_l = []
                l_list = []
                for n in lG.nodes:
                    if (
                        lG.nodes[n]["cnodes_l"][0] == center_nodes[k]
                        and lG.nodes[n]["label"] == "C"
                    ):
                        linker_C_l.append((n, lG.nodes[n]["cnodes_l"]))
                        l_list.append(lG.nodes[n]["cnodes_l"][1])

                outer_connected_C_ind = [
                    ind
                    for ind, value in enumerate(l_list)
                    if value == (max(l_list) - 1)
                ]

                if len(outer_connected_C_ind) == 1:
                    outer_X = linker_C_l[outer_connected_C_ind[0]]
                    if center_nodes[k] not in [outer_X[0]]:
                        print("find connected X in edge:  ", outer_X[0])
                        Xs_indices += [outer_X[0]]

        if center_class == "onepoint":
            print("ditopic linker: center is a point")
            Xs_indices = []
            linker_C_l = []
            l_list = []
            for n in lG.nodes:
                if (
                    lG.nodes[n]["cnodes_l"][0] == center_nodes[0]
                    and lG.nodes[n]["label"] == "C"
                ):
                    linker_C_l.append((n, lG.nodes[n]["cnodes_l"]))
                    l_list.append(lG.nodes[n]["cnodes_l"][1])

            outer_connected_C_ind = [
                ind for ind, value in enumerate(l_list) if value == (max(l_list) - 1)
            ]
            for m in range(len(outer_connected_C_ind)):
                outer_X = linker_C_l[outer_connected_C_ind[m]]
                print("find connected X in edge:  ", outer_X[0])
                Xs_indices += [outer_X[0]]
        if center_class == "cycle":
            print("ditopic linker: center is a cycle")
            connected_pairXs = {}
            Xs_indices = []
            for k in range(len(center_nodes)):
                linker_C_l = []
                l_list = []
                for n in lG.nodes:
                    if (
                        lG.nodes[n]["cnodes_l"][0] == center_nodes[k]
                        and lG.nodes[n]["label"] == "C"
                    ):
                        linker_C_l.append((n, lG.nodes[n]["cnodes_l"]))
                        l_list.append(lG.nodes[n]["cnodes_l"][1])
                outer_connected_C_ind = [
                    ind
                    for ind, value in enumerate(l_list)
                    if value == (max(l_list) - 1)
                ]
                if len(outer_connected_C_ind) == 1:
                    outer_X = linker_C_l[outer_connected_C_ind[0]]
                    if center_nodes[k] not in [outer_X[0]]:
                        print("find connected X in edge:  ", outer_X[0])
                        Xs_indices += [outer_X[0]]

            if len(Xs_indices) < 2:
                print("Xs in the center cycle")
                # the linker is a cycle, but no Xs are found by the dist, then the X in in the center cycle:
                # the node whose adjacents(nonH) more than 2 are the Xs
                for n in center_nodes:
                    adj_nonH_num = 0
                    if lG.nodes[n]["label"] == "C":
                        adj_nodes = list(lG.adj[n])
                        for adj in adj_nodes:
                            if lG.nodes[adj]["label"] != "H":
                                adj_nonH_num += 1
                        if adj_nonH_num > 2:
                            Xs_indices.append(n)

    else:
        raise ValueError(
            "failed to recognize a multitopic linker whose center is not a cycle"
        )

    # write cifs

    if linker_topic > 2:  # multitopic
        frag_nodes = list(sorted(nx.connected_components(lG), key=len, reverse=True))
        for f in frag_nodes:
            if set(center_nodes) < set(f):
                center_frag_nodes = f
            else:
                outer_frag_nodes = f

        subgraph_center_frag = lG.subgraph(center_frag_nodes)
        lines, center_Xs = lines_of_center_frag(
            subgraph_center_frag,
            Xs_indices,
            metals,
            labels,
            coords,
            mass_center_angstrom,
        )
        center_frag_bonds = get_bonds_from_subgraph(subgraph_center_frag, Xs_indices)
        subgraph_center_frag_edges = list(subgraph_center_frag.edges)
        # plot2dedge(lG,coords,subgraph_center_frag_edges,True)
        # plot2dedge(lG,coords,subgraph_center_frag_edges,False)
        pairXs = get_pairX_outer_frag(connected_pairXs, outer_frag_nodes)
        subgraph_single_frag = cleave_outer_frag_subgraph(lG, pairXs, outer_frag_nodes)
        rows, frag_Xs = lines_of_single_frag(subgraph_single_frag, Xs_indices)
        single_frag_bonds = get_bonds_from_subgraph(subgraph_single_frag, Xs_indices)
        if linker_topic == 3:
            print(
                "linker_center_frag:", subgraph_center_frag.number_of_nodes(), center_Xs
            )
            print("linker_outer_frag:", subgraph_single_frag.number_of_nodes(), frag_Xs)
            linker_center_node_pdb_name = Path(save_nodes_dir, "tricenter")
            create_pdb(linker_center_node_pdb_name, lines)
            linker_branch_pdb_name = Path(save_edges_dir, "triedge")
            create_pdb(linker_branch_pdb_name, rows)
            # create_cif(lines,center_frag_bonds,save_nodes_dir,'tricenter.cif')
            # create_cif(rows,single_frag_bonds,save_edges_dir,'triedge.cif')
            return (
                subgraph_center_frag.number_of_nodes(),
                center_Xs,
                subgraph_single_frag.number_of_nodes(),
                frag_Xs,
                linker_center_node_pdb_name + ".pdb",
                linker_branch_pdb_name + ".pdb",
            )
        elif linker_topic == 4:
            print("center_frag:", subgraph_center_frag.number_of_nodes(), center_Xs)
            print("outer_frag:", subgraph_single_frag.number_of_nodes(), frag_Xs)
            linker_center_node_pdb_name = str(Path(save_nodes_dir, "tetracenter"))
            create_pdb(linker_center_node_pdb_name, lines)
            linker_branch_pdb_name = str(Path(save_edges_dir, "tetraedge"))
            create_pdb(linker_branch_pdb_name, rows)
            # create_cif(lines,center_frag_bonds,save_nodes_dir,'tetracenter.cif')
            # create_cif(rows,single_frag_bonds,save_edges_dir,'tetraedge.cif')
            return (
                subgraph_center_frag.number_of_nodes(),
                center_Xs,
                subgraph_single_frag.number_of_nodes(),
                frag_Xs,
                linker_center_node_pdb_name + ".pdb",
                linker_branch_pdb_name + ".pdb",
            )
        else:
            linker_center_node_pdb_name = str(Path(save_nodes_dir, "multicenter"))
            create_pdb(linker_center_node_pdb_name, lines)
            linker_branch_pdb_name = str(Path(save_edges_dir, "multiedge"))
            create_pdb(linker_branch_pdb_name, rows)
            # create_cif(lines,center_frag_bonds,'nodes','multicenter.cif')
            # create_cif(rows,single_frag_bonds,'edges','multiedge.cif')
            return (
                subgraph_center_frag.number_of_nodes(),
                center_Xs,
                subgraph_single_frag.number_of_nodes(),
                frag_Xs,
                linker_center_node_pdb_name + ".pdb",
                linker_branch_pdb_name + ".pdb",
            )

    elif linker_topic == 2:  # ditopic
        pairXs = Xs_indices
        print("pairXs:", pairXs)
        subgraph_center_frag = cleave_outer_frag_subgraph(lG, pairXs, lG.nodes)
        subgraph_center_frag_edges = list(subgraph_center_frag.edges)
        # plot2dedge(lG,coords,subgraph_center_frag_edges,True)
        # plot2dedge(subgraph_center_frag,coords,subgraph_center_frag_edges,False)
        lines, center_Xs = lines_of_center_frag(
            subgraph_center_frag,
            Xs_indices,
            metals,
            labels,
            coords,
            mass_center_angstrom,
        )
        center_frag_bonds = get_bonds_from_subgraph(subgraph_center_frag, Xs_indices)
        # create_cif(lines,center_frag_bonds,'edges','diedge.cif')
        edge_pdb_name = str(Path(save_edges_dir, "diedge"))
        create_pdb(edge_pdb_name, lines)
        print("linker_center_frag:", subgraph_center_frag.number_of_nodes(), center_Xs)
        return (
            subgraph_center_frag.number_of_nodes(),
            center_Xs,
            0,
            [],
            None,
            edge_pdb_name + ".pdb",
        )


##############below are from learn_template#####################


def make_supercell_3x3x3(array_xyz):
    array_x1 = array_xyz + np.array([1, 0, 0])
    array_x2 = array_xyz + np.array([-1, 0, 0])
    array_y1 = array_xyz + np.array([0, 1, 0])
    array_y2 = array_xyz + np.array([0, -1, 0])
    array_x1_y1 = array_xyz + np.array([1, 1, 0])
    array_x1_y2 = array_xyz + np.array([1, -1, 0])
    array_x2_y1 = array_xyz + np.array([-1, 1, 0])
    array_x2_y2 = array_xyz + np.array([-1, -1, 0])

    layer_3x3 = np.vstack(
        (
            array_xyz,
            array_x1,
            array_x2,
            array_y1,
            array_y2,
            array_x1_y1,
            array_x1_y2,
            array_x2_y1,
            array_x2_y2,
        )
    )

    layer_3x3_z1 = layer_3x3 + np.array([0, 0, 1])
    layer_3x3_z2 = layer_3x3 + np.array([0, 0, -1])

    supercell_3x3x3 = np.vstack((layer_3x3, layer_3x3_z1, layer_3x3_z2))

    return supercell_3x3x3


def check_inside_unit_cell(point):
    return all([i >= -0.0 and i < 1.0 for i in point])


def extract_unit_cell(cell_info):
    pi = np.pi
    aL, bL, cL, alpha, beta, gamma = cell_info
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


def find_pair_v_e(vvnode333, eenode333):
    G = nx.Graph()
    pair_ve = []
    for e in eenode333:
        dist = np.linalg.norm(vvnode333 - e, axis=1)
        # find two v which are nearest to e, and at least one v is in [0,1] unit cell
        v1 = vvnode333[np.argmin(dist)]
        v1_idx = np.argmin(dist)
        dist[np.argmin(dist)] = 1000
        v2 = vvnode333[np.argmin(dist)]
        v2_idx = np.argmin(dist)
        # find the center of the pair of v
        center = (v1 + v2) / 2
        # check if there is a v in [0,1] unit cell
        if check_inside_unit_cell(v1) or check_inside_unit_cell(v2):
            # check if the center of the pair of v is around e
            if np.linalg.norm(center - e) < 1e-3:
                G.add_node("V" + str(v1_idx), fcoords=v1, note="V", type="V")
                G.add_node("V" + str(v2_idx), fcoords=v2, note="V", type="V")
                (
                    G.add_edge(
                        "V" + str(v1_idx),
                        "V" + str(v2_idx),
                        fcoords=(v1, v2),
                        fc_center=e,
                    ),
                )
                pair_ve.append((v1, v2, e))
    return pair_ve, len(pair_ve), G


def find_pair_v_e_c(
    vvnode333, ecnode333, eenode333, unit_cell
):  # exist center of linker  in mof
    G = nx.Graph()
    pair_ve = []
    for e in eenode333:
        # print(e, "check")
        # dist_v_e = np.linalg.norm(vvnode333 - e, axis=1)
        dist_v_e = np.linalg.norm(np.dot(unit_cell, (vvnode333 - e).T).T, axis=1)
        # find two v which are nearest to e, and at least one v is in [0,1] unit cell
        v1 = vvnode333[np.argmin(dist_v_e)]
        v1_idx = np.argmin(dist_v_e)
        dist_c_e = np.linalg.norm(np.dot(unit_cell, (ecnode333 - e).T).T, axis=1)
        # find two v which are nearest to e, and at least one v is in [0,1] unit cell
        v2 = ecnode333[np.argmin(dist_c_e)]
        v2_idx = np.argmin(dist_c_e)
        # print(v1, v2, "v1,v2")

        # find the center of the pair of v
        center = (v1 + v2) / 2
        # check if there is a v in [0,1] unit cell
        if check_inside_unit_cell(v1) or check_inside_unit_cell(v2):
            # check if the center of the pair of v is around e
            # if abs(np.linalg.norm(v1 - e)+np.linalg.norm(v2 - e) - np.linalg.norm(v1 - v2))< 1e-2: #v1,v2,e are colinear
            if np.linalg.norm(center - e) < 0.1:
                # print(e,v1,v2,'check')
                G.add_node("V" + str(v1_idx), fcoords=v1, note="V")
                G.add_node("CV" + str(v2_idx), fcoords=v2, note="CV")
                (
                    G.add_edge(
                        "V" + str(v1_idx),
                        "CV" + str(v2_idx),
                        fcoords=(v1, v2),
                        fc_center=e,
                    ),
                )
                pair_ve.append(("V" + str(v1_idx), "CV" + str(v2_idx), e))

    return pair_ve, len(pair_ve), G


# add ccoords to the the nodes in the graph
def add_ccoords(G, unit_cell):
    for n in G.nodes():
        G.nodes[n]["ccoords"] = np.dot(unit_cell, G.nodes[n]["fcoords"])
    return G


# check if after np.mod, the fcoords is the same as before
def check_moded_fcoords(point):
    x, y, z = point[0], point[1], point[2]
    if np.mod(x, 1) != x:
        return False
    if np.mod(y, 1) != y:
        return False
    if np.mod(z, 1) != z:
        return False
    return True


def set_DV_V(G):
    for n in G.nodes():
        ##if G.degree(n) == max(dict(G.degree()).values()):
        ##    G.nodes[n]['type'] = 'V'
        ##    #check if the moded ccoords is in the unit cell
        if check_moded_fcoords(G.nodes[n]["fcoords"]):
            G.nodes[n]["type"] = "V"
            G = nx.relabel_nodes(
                G, {n: pname(n) + "_" + str(np.array([0.0, 0.0, 0.0]))}
            )

        else:
            G.nodes[n]["type"] = "DV"
            # rename node name to V_diff
            diff = G.nodes[n]["fcoords"] - np.mod(G.nodes[n]["fcoords"], 1)
            # find the corresponding V (with type V)node by the moded fcoords
            for n1 in G.nodes():
                if np.all(
                    np.isclose(G.nodes[n1]["fcoords"], np.mod(G.nodes[n]["fcoords"], 1))
                ):
                    n1_name = pname(n1)
            # replace node name
            G = nx.relabel_nodes(G, {n: n1_name + "_" + str(diff)})

    max_degree = max(dict(G.degree()).values())
    return G, max_degree


# if make supercell333 of unique_e list cannot find the e_new, then this e_new should be appended to the unique_e list
# check if the e_new is in the unique_e list, if yes, then this e_new is valid, if no, then this e_new is invalid
def find_unitcell_e(all_e):
    cell_e = [all_e[0]]
    count = 1
    while count < len(all_e):
        e_check = all_e[count]
        supercell_e = make_supercell_3x3x3(cell_e)
        if not np.any(np.all(np.isclose(e_check, supercell_e), axis=1)):
            cell_e.append(e_check)
        count += 1
    return cell_e


# check e in G, find e in unit_cell, use np.mod to filter the e in unit_cell, set the valid e with E type, others are DE type
def set_DE_E(G):
    all_e = []
    for e in G.edges():
        all_e.append(G.edges[e]["fc_center"].copy())
    # all_e = np.vstack([limit_to_abs1(edge) for edge in all_e])
    # limit x,y,z of e to [-1,1]
    unique_e = np.vstack(find_unitcell_e(all_e))
    # print(unique_e) #debug
    for e in G.edges():
        if np.any(np.all(np.isclose(G.edges[e]["fc_center"], unique_e), axis=1)):
            G.edges[e]["type"] = "E"
            # print(G.edges[e]['fc_center'],'E') #debug
        else:
            G.edges[e]["type"] = "DE"
            # print(G.edges[e]['fc_center'],'DE') #debug
    return G


##########below are from terminations.py#####################


def termpdb(filename):
    inputfile = str(filename)
    with open(inputfile, "r") as fp:
        content = fp.readlines()
        # linesnumber = len(content)
    data = []
    for line in content:
        line = line.strip()
        if len(line) > 0:  # skip blank line
            if line[0:6] == "ATOM" or line[0:6] == "HETATM":
                value_atom = line[12:16].strip()  # atom_label
                # resname
                # value2 = 'MOL'  # res_name

                value_x = float(line[30:38])  # x
                value_y = float(line[38:46])  # y
                value_z = float(line[46:54])  # z
                value_charge = float(line[61:66])
                value_note = line[67:80].strip()  # atom_note
                # resnumber
                try:
                    value_res_num = int(line[22:26])
                except ValueError:
                    value_res_num = 1
                data.append(
                    [
                        value_atom,
                        value_charge,
                        value_note,
                        value_res_num,
                        "TERM",
                        value_res_num,
                        value_x,
                        value_y,
                        value_z,
                    ]
                )
    return np.vstack(data)


def Xpdb(data, X):
    indices = [i for i in range(len(data)) if nn(data[i, 0]) == X]
    X_term = data[indices]
    return X_term, indices


######below are from _superimpose.py#####################


def sort_by_distance(arr):
    # Calculate distances from the first element to all other elements
    distances = [(np.linalg.norm(arr[0] - arr[i]), i) for i in range(len(arr))]
    # Sort distances in ascending order
    distances.sort(key=lambda x: x[0])
    return distances


def match_vectors(arr1, arr2, num):
    # Get sorted distances
    sorted_distances_arr1 = sort_by_distance(arr1)
    sorted_distances_arr2 = sort_by_distance(arr2)

    # Select the indices by distance matching in limited number

    indices_arr1 = [sorted_distances_arr1[j][1] for j in range(num)]
    indices_arr2 = [sorted_distances_arr2[j][1] for j in range(num)]

    # reorder the matching vectors# which can induce the smallest RMSD
    closest_vectors_arr1 = np.array([arr1[i] for i in indices_arr1])
    closest_vectors_arr2 = np.array([arr2[i] for i in indices_arr2])

    return closest_vectors_arr1, closest_vectors_arr2


def superimpose(arr1, arr2, min_rmsd=1e6):
    arr1 = np.asarray(arr1)
    arr2 = np.asarray(arr2)
    m_arr1, m_arr2 = match_vectors(arr1, arr2, min(6, len(arr1), len(arr2)))
    best_rot, best_tran = np.eye(3), np.zeros(3)

    for perm in itertools.permutations(m_arr1):
        rmsd, rot, tran = svd_superimpose(np.asarray(perm), m_arr2)
        if rmsd < min_rmsd:
            min_rmsd, best_rot, best_tran = rmsd, rot, tran

    return min_rmsd, best_rot, best_tran


def svd_superimpose(inp_arr1, inp_arr2):
    """
    Calculates RMSD and rotation matrix for superimposing two sets of points,
    using SVD. Ref.: "Least-Squares Fitting of Two 3-D Point Sets", IEEE
    Transactions on Pattern Analysis and Machine Intelligence, 1987, PAMI-9(5),
    698-700. DOI: 10.1109/TPAMI.1987.4767965
    """

    arr1 = np.array(inp_arr1)
    arr2 = np.array(inp_arr2)

    com1 = np.sum(arr1, axis=0) / arr1.shape[0]
    com2 = np.sum(arr2, axis=0) / arr2.shape[0]

    arr1 -= com1
    arr2 -= com2

    cov_mat = np.matmul(arr1.T, arr2)
    U, s, Vt = np.linalg.svd(cov_mat)

    rot_mat = np.matmul(U, Vt)
    if np.linalg.det(rot_mat) < 0:
        Vt[-1, :] *= -1.0
        rot_mat = np.matmul(U, Vt)

    diff = arr2 - np.matmul(arr1, rot_mat)
    rmsd = np.sqrt(np.sum(diff**2) / diff.shape[0])
    trans = com2 - np.dot(com1, rot_mat)

    return rmsd, rot_mat, trans


def superimpose_rotation_only(arr1, arr2, min_rmsd=1e6):
    arr1 = np.asarray(arr1)
    arr2 = np.asarray(arr2)
    m_arr1, m_arr2 = match_vectors(arr1, arr2, min(6, len(arr1), len(arr2)))
    best_rot, best_tran = np.eye(3), np.zeros(3)
    for perm in itertools.permutations(m_arr1):
        rmsd, rot, tran = svd_superimpose(np.asarray(perm), m_arr2)
        if rmsd < min_rmsd:
            min_rmsd, best_rot, best_tran = rmsd, rot, tran
            if np.allclose(np.dot(best_tran, np.zeros(3)), 1e-1):
                break

    return min_rmsd, best_rot, best_tran


##########below are from _place_node_edge.py#####################
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


# Add row indices as the first column
def addidx(array):
    row_indices = np.arange(array.shape[0]).reshape(-1, 1).astype(int)
    new_array = np.hstack((row_indices, array))
    return new_array


def get_edge_lengths(G):
    edge_lengths = {}
    lengths = []
    for e in G.edges():
        i, j = e
        length = np.linalg.norm(G.nodes[i]["ccoords"] - G.nodes[j]["ccoords"])
        length = np.round(length, 3)
        edge_lengths[(i, j)] = length
        edge_lengths[(j, i)] = length
        lengths.append(length)
    # print('edge lengths:',set(lengths)) #debug
    if len(set(lengths)) != 1:
        print("more than one type of edge length")
        # if the length are close, which can be shown by std
        if np.std(lengths) < 0.1:
            print("the edge lengths are close")
        else:
            print("the edge lengths are not close")
        print(set(lengths))
    return edge_lengths, set(lengths)


def update_node_ccoords(G, edge_lengths, start_node, new_edge_length):
    updated_ccoords = {}
    original_ccoords = {}
    updated_ccoords[start_node] = G.nodes[start_node]["ccoords"]
    original_ccoords[start_node] = G.nodes[start_node]["ccoords"]
    updated_node = [start_node]
    # to update all the nodes start from the startnode and spread to neighbors
    for i in range(len(G.nodes()) - 1):
        for n in updated_node:
            for nn in G.neighbors(n):
                if nn in updated_node:
                    continue
                edge = (n, nn)
                edge_length = edge_lengths[edge]
                updated_ccoords[nn] = (
                    updated_ccoords[n]
                    + (G.nodes[nn]["ccoords"] - G.nodes[n]["ccoords"])
                    * new_edge_length
                    / edge_length
                )
                original_ccoords[nn] = G.nodes[nn]["ccoords"]
                updated_node.append(nn)

    return updated_ccoords, original_ccoords


###below are from multiedge_bundling.py#######
def find_pair_x_edge_fc(x_matrix, edge_matrix, sc_unit_cell):
    assert_msg_critical(
        "scipy" in sys.modules, "scipy is required for find_pair_x_edge."
    )
    dist_matrix = np.zeros((len(x_matrix), len(edge_matrix)))
    x_matrix = fractional_to_cartesian(x_matrix, sc_unit_cell)
    edge_matrix = fractional_to_cartesian(edge_matrix, sc_unit_cell)
    for i in range(len(x_matrix)):
        for j in range(len(edge_matrix)):
            dist_matrix[i, j] = np.linalg.norm(x_matrix[i] - edge_matrix[j])
    row_ind, col_ind = linear_sum_assignment(dist_matrix)
    ##print(row_ind, col_ind) # debug
    return row_ind, col_ind


def bundle_multiedge(sG):
    multiedge_bundlings = []
    for n in sG.nodes:
        if "CV" in n and sG.nodes[n]["type"] == "V":
            edges = []
            for con_n in list(sG.neighbors(n)):
                edges.append(sG.edges[n, con_n]["coords"])
                # edges = np.vstack(edges, dtype=float)
                ## extract the X atom in the CV node
                # cv_xatoms = np.asarray(
                #    sG.nodes[n]["x_coords"][:, 2:5], dtype=float
                # )  # modified for extra column of atom type
                # if len(cv_xatoms) < len(edges):
                #    # duplicate the cv_xatoms
                #    cv_xatoms = np.vstack([cv_xatoms] * len(edges))
                ##_, order = find_pair_x_edge(cv_xatoms, edges)
                ##con_n_order = [list(sG.neighbors(n))[i] for i in order]
            multiedge_bundlings.append((n, list(sG.neighbors(n))))
            # make dist matrix of each x to each edge and then use hugerian algorithm to find the pair of x and edge
    return multiedge_bundlings


##########below are from make_eG.py#####################


def order_edge_array(row_ind, col_ind, edges_array):
    # according col_ind order  to reorder the connected edge points
    old_split = np.vsplit(edges_array, len(col_ind))
    old_order = []
    for i in range(len(col_ind)):
        old_order.append(
            (row_ind[i], col_ind[i], old_split[sorted(col_ind).index(col_ind[i])])
        )
    new_order = sorted(old_order, key=lambda x: x[0])
    ordered_arr = np.vstack([new_order[j][2] for j in range(len(new_order))])
    return ordered_arr


def make_paired_Xto_x(ec_arr, merged_arr, neighbor_number, sc_unit_cell):
    """
    ec_arr: the array of the linker center node
    merged_arr: the array of the merged edges
    neighbor_number: the number of the neighbor nodes of the linker center node
    """
    ec_indices, ec_fpoints = fetch_X_atoms_ind_array(ec_arr, 0, "X")
    if len(ec_indices) < neighbor_number:
        # duplicate the cv_xatoms
        ec_fpoints = np.vstack([ec_fpoints] * neighbor_number)
    nei_indices, nei_fcpoints = fetch_X_atoms_ind_array(
        merged_arr, 0, "X"
    )  # extract only the X atoms in the neighbor edges but not the cv nodes: [len(ec_arr) :]
    row_ind, col_ind = find_pair_x_edge_fc(
        ec_fpoints[:, 2:5].astype(float),
        nei_fcpoints[:, 2:5].astype(float),
        sc_unit_cell,
    )  # NOTE: modified to skip atom type

    # according col_ind order  to reorder the connected edge points
    # switch the X to x for nei_fcpoints
    for i in col_ind:
        if nn(merged_arr[nei_indices[i], 0]) == "X":
            merged_arr[nei_indices[i], 0] = "x" + nl(merged_arr[nei_indices[i], 0])
    for k in row_ind:
        if ec_indices[k] < len(ec_arr):
            if nn(ec_arr[ec_indices[k], 0]) == "X":
                ec_arr[ec_indices[k], 0] = "x" + nl(ec_arr[ec_indices[k], 0])

    ordered_edges_points_follow_ecXatoms = order_edge_array(
        row_ind, col_ind, merged_arr
    )
    # remove the duplicated cv_xatoms
    ec_merged_arr = np.vstack((ec_arr, ordered_edges_points_follow_ecXatoms))
    return ec_merged_arr


def superG_to_eG_multitopic(superG, sc_unit_cell):
    # CV + neighbor V -> E+index node
    eG = nx.Graph()
    edge_count = 0
    node_count = 0
    for n in superG.nodes():
        if superG.nodes[n]["note"] == "V":
            node_count += 1
            eG.add_node(
                n,
                f_points=superG.nodes[n]["f_points"],
                fcoords=superG.nodes[n]["fcoords"],
                type="V",
                name="NODE",
                note="V",
                index=node_count,
            )
            superG.nodes[n]["index"] = node_count
            # add virtual edge
            for e in superG.edges(n):
                if superG.edges[e]["type"] == "virtual":
                    eG.add_edge(e[0], e[1], type="virtual")

        elif superG.nodes[n]["note"] == "CV":
            edge_count -= 1
            neighbors = list(superG.neighbors(n))
            merged_edges = superG.nodes[n]["f_points"]
            # reorder neighbors according to the order of X atoms in CV node

            for ne in neighbors:
                merged_edges = np.vstack(
                    [superG.edges[n, ne]["f_points"] for ne in neighbors]
                )

                # follow the order of X atoms in CV node

            # for x atoms in merged_edges, use hungarian algorithm to find the nearest X-X atoms in the neighbor nodes and replace the X to x
            ec_merged_edges = make_paired_Xto_x(
                superG.nodes[n]["f_points"], merged_edges, len(neighbors), sc_unit_cell
            )
            eG.add_node(
                "EDGE_" + str(edge_count),
                f_points=ec_merged_edges,
                fcoords=superG.nodes[n]["fcoords"],
                type="Edge",
                name="EDGE",
                note="E",
                index=edge_count,
            )

            for ne in neighbors:
                eG.add_edge(
                    "EDGE_" + str(edge_count),
                    ne,
                    index="E_" + str(edge_count),
                    type="real",
                )

    return eG, superG


def superG_to_eG_ditopic(superG):
    # V + neighbor V -> E+index node
    eG = nx.Graph()
    edge_count = 0
    node_count = 0
    edge_record = []
    for n in superG.nodes():
        if superG.nodes[n]["note"] == "V":
            node_count += 1
            eG.add_node(
                n,
                f_points=superG.nodes[n]["f_points"],
                fcoords=superG.nodes[n]["fcoords"],
                type=superG.nodes[n]["type"],
                note=superG.nodes[n]["note"],
                name="NODE",
                index=node_count,
            )
            # print('add node',n,'type',superG.nodes[n]['type']) #debug
            superG.nodes[n]["index"] = node_count
            # add virtual edge
            for e in superG.edges(n):
                if superG.edges[e]["type"] == "virtual":
                    eG.add_edge(e[0], e[1], type="virtual")

            neighbors = list(superG.neighbors(n))
            for ne in neighbors:
                if sorted([n, ne]) in edge_record:
                    continue
                edge_record.append(sorted([n, ne]))
                edge_count -= 1
                eG.add_node(
                    "EDGE_" + str(edge_count),
                    f_points=superG.edges[n, ne]["f_points"],
                    fcoords=superG.edges[n, ne]["fcoords"],
                    type="Edge",
                    name="EDGE",
                    note="E",
                    index=edge_count,
                )
                eG.add_edge(n, ne, index="E_" + str(edge_count), type="real")
                eG.add_edge(
                    "EDGE_" + str(edge_count),
                    ne,
                    index="E_" + str(edge_count),
                    type="half",
                )
                eG.add_edge(
                    "EDGE_" + str(edge_count),
                    n,
                    index="E_" + str(edge_count),
                    type="half",
                )

    return eG, superG


def remove_node_by_index(eG, remove_node_list, remove_edge_list):
    for n in eG.nodes():
        if pname(n) != "EDGE":
            if eG.nodes[n]["index"] in remove_node_list:
                eG.remove_node(n)
        if pname(n) == "EDGE":
            if -1 * eG.nodes[n]["index"] in remove_edge_list:
                eG.remove_node(n)
    return eG


def find_nearest_neighbor(i, n_n_distance_matrix):
    n_n_min_distance = np.min(n_n_distance_matrix[i : i + 1, :])
    _, n_j = locate_min_idx(n_n_distance_matrix[i : i + 1, :])
    # print('add virtual edge between',nodes_list[i],nodes_list[n_j])
    # n_n_distance_matrix[i, n_j] = 1000
    # set the column to 1000 to avoid the same atom being selected again
    n_n_distance_matrix[:, n_j] = 1000
    return n_j, n_n_min_distance, n_n_distance_matrix


def find_surrounding_points(ind, n_n_distance_matrix, max_number):
    stop = 0  # if while loop is too long, stop it
    nearest_neighbor = {}
    nearest_neighbor[ind] = []
    while len(nearest_neighbor[ind]) < max_number:
        stop += 1
        if stop > 100:
            break
        n_j, _, n_n_distance_matrix = find_nearest_neighbor(ind, n_n_distance_matrix)
        nearest_neighbor[ind].append(n_j)
    return nearest_neighbor


# Function to find 'XOO' pairs for a specific node
def xoo_pair_ind_node(single_node_fc, sc_unit_cell):
    # if the node x is not surrounded by two o atoms,
    # then modify the fetch_X_atoms_ind_array(single_node, 0, 'O') find_surrounding_points(k, xs_os_dist_matrix, 2)
    # this function is to find the XOO pairs in a specific node(by node_id),
    # this xoo pair is the indice of x and nearest two o atoms in the same node
    # return the indice of x and nearest two o atoms in the same node, which can be convert to a dict with x_index as key and o_indices as value
    # the distance is in cartesian coordinates
    # single_node_fc: coordinates of any node in the main fragment
    # sc_unit_cell: supercell unit cell matrix
    single_node = np.hstack(
        (
            single_node_fc[:, 0:1],
            np.dot(sc_unit_cell, single_node_fc[:, 2:5].astype(float).T).T,
        )
    )  # NOTE: modified to skip atom type
    xind, xs_coords = fetch_X_atoms_ind_array(single_node, 0, "X")
    oind, os_coords = fetch_X_atoms_ind_array(single_node, 0, "O")
    xs_os_dist_matrix = np.zeros((len(xs_coords), len(os_coords)))
    for i in range(len(xs_coords)):
        for j in range(len(os_coords)):
            xs_os_dist_matrix[i, j] = np.linalg.norm(
                xs_coords[i, 1:4].astype(float) - os_coords[j, 1:4].astype(float)
            )
    xoo_ind_list = []
    for k in range(len(xind)):
        nearest_dict = find_surrounding_points(k, xs_os_dist_matrix, 2)
        for key in nearest_dict.keys():
            xoo_ind_list.append(
                [xind[key], sorted([oind[m] for m in nearest_dict[key]])]
            )
    return xoo_ind_list


def get_xoo_dict_of_node(eG, sc_unit_cell):
    # quick check the order of xoo in every node are same, select n0 and n1, if xoo_ind_node0 == xoo_ind_node1, then xoo_dict is the same
    # return xoo dict of every node, key is x index, value is o index
    n0 = [i for i in eG.nodes() if pname(i) != "EDGE"][0]
    n1 = [i for i in eG.nodes() if pname(i) != "EDGE"][1]
    xoo_ind_node0 = xoo_pair_ind_node(
        eG.nodes[n0]["f_points"], sc_unit_cell
    )  # pick node one and get xoo_ind pair
    xoo_ind_node1 = xoo_pair_ind_node(
        eG.nodes[n1]["f_points"], sc_unit_cell
    )  # pick node two and get xoo_ind pair
    if xoo_ind_node0 == xoo_ind_node1:
        xoo_dict = {}
        for xoo in xoo_ind_node0:
            xoo_dict[xoo[0]] = xoo[1]
    else:
        print("the order of xoo in every node are not same, please check the input")
        print("xoo_ind_node0", xoo_ind_node0)
        print("xoo_ind_node1", xoo_ind_node1)
    return xoo_dict


def addxoo2edge_multitopic(eG, sc_unit_cell):
    xoo_dict = get_xoo_dict_of_node(eG, sc_unit_cell)
    matched_vnode_X = []
    unsaturated_linker = []
    # for every X atom in the EDGE node, search for the paired(nearest) X atom in the connected V node
    # and then use the xoo_dict of the connected V node to extract the xoos of the connected V node
    # and then add the xoos to the EDGE node
    # all xoo_node for the V node is the same
    EDGE_nodes = [n for n in eG.nodes() if pname(n) == "EDGE"]
    for n in EDGE_nodes:
        Xs_edge_indices, Xs_edge_fpoints = fetch_X_atoms_ind_array(
            eG.nodes[n]["f_points"], 0, "X"
        )
        Xs_edge_ccpoints = np.hstack(
            (
                Xs_edge_fpoints[:, 0:2],
                np.dot(sc_unit_cell, Xs_edge_fpoints[:, 2:5].astype(float).T).T,
            )
        )  # NOTE: modified to skip atom type
        V_nodes = [i for i in eG.neighbors(n) if pname(i) != "EDGE"]
        if len(V_nodes) == 0:
            # unsaturated_linker.append(n)
            # print(
            #    "no V node connected to this edge node, this linker is a isolated linker, will be ignored",
            #    n,
            # ) # debug
            continue
        all_Xs_vnodes_ind = []
        all_Xs_vnodes_ccpoints = np.zeros((0, 5))
        for v in V_nodes:
            # find the connected V node and its X atoms
            Xs_vnode_indices, Xs_vnode_fpoints = fetch_X_atoms_ind_array(
                eG.nodes[v]["f_points"], 0, "X"
            )
            Xs_vnode_ccpoints = np.hstack(
                (
                    Xs_vnode_fpoints[:, 0:2],
                    np.dot(sc_unit_cell, Xs_vnode_fpoints[:, 2:5].astype(float).T).T,
                )
            )  # NOTE: modified to skip atom type
            for ind in Xs_vnode_indices:
                all_Xs_vnodes_ind.append([v, ind, n])
            all_Xs_vnodes_ccpoints = np.vstack(
                (all_Xs_vnodes_ccpoints, Xs_vnode_ccpoints)
            )
        edgeX_vnodeX_dist_matrix = np.zeros(
            (len(Xs_edge_ccpoints), len(all_Xs_vnodes_ccpoints))
        )
        for i in range(len(Xs_edge_ccpoints)):
            for j in range(len(all_Xs_vnodes_ccpoints)):
                edgeX_vnodeX_dist_matrix[i, j] = np.linalg.norm(
                    Xs_edge_ccpoints[i, 2:5].astype(float)
                    - all_Xs_vnodes_ccpoints[j, 2:5].astype(float)
                )

        for k in range(len(Xs_edge_fpoints)):
            n_j, min_dist, edgeX_vnodeX_dist_matrix = find_nearest_neighbor(
                k, edgeX_vnodeX_dist_matrix
            )

            if min_dist > 4:
                unsaturated_linker.append(n)
                # print(
                #    "no xoo for edge node, this linker is a dangling unsaturated linker",
                #    n,
                # ) # debug
                continue
            # add the xoo to the edge node

            nearest_vnode = all_Xs_vnodes_ind[n_j][0]
            nearest_X_ind_in_vnode = all_Xs_vnodes_ind[n_j][1]
            matched_vnode_X.append(all_Xs_vnodes_ind[n_j])
            corresponding_o_indices = xoo_dict[nearest_X_ind_in_vnode]
            xoo_ind_in_vnode = [[nearest_X_ind_in_vnode] + corresponding_o_indices]
            xoo_fpoints_in_vnode = [
                eG.nodes[nearest_vnode]["f_points"][i] for i in xoo_ind_in_vnode
            ]
            xoo_fpoints_in_vnode = np.vstack(xoo_fpoints_in_vnode)
            eG.nodes[n]["f_points"] = np.vstack(
                (eG.nodes[n]["f_points"], xoo_fpoints_in_vnode)
            )
            # print('add xoo to edge node',n) #debug
    return eG, unsaturated_linker, matched_vnode_X, xoo_dict


def addxoo2edge_ditopic(eG, sc_unit_cell):
    xoo_dict = get_xoo_dict_of_node(eG, sc_unit_cell)
    matched_vnode_X = []
    unsaturated_linker = []
    # for every X atom in the EDGE node, search for the paired(nearest) X atom in the connected V node
    # and then use the xoo_dict of the connected V node to extract the xoos of the connected V node
    # and then add the xoos to the EDGE node
    # all xoo_node for the V node is the same
    EDGE_nodes = [n for n in eG.nodes() if pname(n) == "EDGE"]
    for n in EDGE_nodes:
        Xs_edge_indices, Xs_edge_fpoints = fetch_X_atoms_ind_array(
            eG.nodes[n]["f_points"], 0, "X"
        )
        Xs_edge_ccpoints = np.hstack(
            (
                Xs_edge_fpoints[:, 0:2],
                np.dot(sc_unit_cell, Xs_edge_fpoints[:, 2:5].astype(float).T).T,
            )
        )  # NOTE: modified to skip atom type
        V_nodes = [i for i in eG.neighbors(n) if pname(i) != "EDGE"]
        if len(V_nodes) == 0:
            # unsaturated_linker.append(n)
            print(
                "no V node connected to this edge node, this linker is a isolated linker, will be ignored",
                n,
            )
            continue
        all_Xs_vnodes_ind = []
        all_Xs_vnodes_ccpoints = np.zeros((0, 5))
        for v in V_nodes:
            # find the connected V node
            Xs_vnode_indices, Xs_vnode_fpoints = fetch_X_atoms_ind_array(
                eG.nodes[v]["f_points"], 0, "X"
            )
            Xs_vnode_ccpoints = np.hstack(
                (
                    Xs_vnode_fpoints[:, 0:2],
                    np.dot(sc_unit_cell, Xs_vnode_fpoints[:, 2:5].astype(float).T).T,
                )
            )  # NOTE: modified to skip atom type
            for ind in Xs_vnode_indices:
                all_Xs_vnodes_ind.append([v, ind, n])
            all_Xs_vnodes_ccpoints = np.vstack(
                (all_Xs_vnodes_ccpoints, Xs_vnode_ccpoints)
            )
        edgeX_vnodeX_dist_matrix = np.zeros(
            (len(Xs_edge_ccpoints), len(all_Xs_vnodes_ccpoints))
        )
        for i in range(len(Xs_edge_ccpoints)):
            for j in range(len(all_Xs_vnodes_ccpoints)):
                edgeX_vnodeX_dist_matrix[i, j] = np.linalg.norm(
                    Xs_edge_ccpoints[i, 2:54].astype(float)
                    - all_Xs_vnodes_ccpoints[j, 2:5].astype(float)
                )
        for k in range(len(Xs_edge_fpoints)):
            n_j, min_dist, _ = find_nearest_neighbor(k, edgeX_vnodeX_dist_matrix)
            if min_dist > 4:
                unsaturated_linker.append(n)
                # print(
                #   min_dist,
                #   "no xoo for edge node, this linker is a dangling unsaturated linker",
                #   n,
                # ) # debug
                continue
            # add the xoo to the edge node
            nearest_vnode = all_Xs_vnodes_ind[n_j][0]
            nearest_X_ind_in_vnode = all_Xs_vnodes_ind[n_j][1]
            matched_vnode_X.append(all_Xs_vnodes_ind[n_j])
            corresponding_o_indices = xoo_dict[nearest_X_ind_in_vnode]
            xoo_ind_in_vnode = [[nearest_X_ind_in_vnode] + corresponding_o_indices]
            xoo_fpoints_in_vnode = [
                eG.nodes[nearest_vnode]["f_points"][i] for i in xoo_ind_in_vnode
            ]
            xoo_fpoints_in_vnode = np.vstack(xoo_fpoints_in_vnode)
            eG.nodes[n]["f_points"] = np.vstack(
                (eG.nodes[n]["f_points"], xoo_fpoints_in_vnode)
            )
            # print('add xoo to edge node',n) #debug
    return eG, unsaturated_linker, matched_vnode_X, xoo_dict


def recenter_and_norm_vectors(vectors, extra_mass_center=None):
    vectors = np.array(vectors)
    if extra_mass_center is not None:
        mass_center = extra_mass_center
    else:
        mass_center = np.mean(vectors, axis=0)
    vectors = vectors - mass_center
    vectors = vectors / np.linalg.norm(vectors, axis=1)[:, None]
    return vectors, mass_center


def get_connected_nodes_vectors(node, G):
    # use adjacent nodes to get vectors
    vectors = []
    for i in list(G.neighbors(node)):
        vectors.append(G.nodes[i]["ccoords"])
    return vectors, G.nodes[node]["ccoords"]


def get_rot_trans_matrix(node, G, sorted_nodes, Xatoms_positions_dict):
    node_id = sorted_nodes.index(node)
    node_xvecs = Xatoms_positions_dict[node_id][:, 1:]
    vecsA, _ = recenter_and_norm_vectors(node_xvecs, extra_mass_center=None)
    v2, node_center = get_connected_nodes_vectors(node, G)
    vecsB, _ = recenter_and_norm_vectors(v2, extra_mass_center=node_center)
    rmsd, rot, tran = superimpose_rotation_only(vecsA, vecsB)
    return rot, tran


def find_unsaturated_node(eG, node_topics):
    # find unsaturated node V in eG
    unsaturated_node = []
    for n in eG.nodes():
        if pname(n) != "EDGE":
            real_neighbor = []
            for cn in eG.neighbors(n):
                if eG.edges[(n, cn)]["type"] == "real":
                    real_neighbor.append(cn)
            if len(real_neighbor) < node_topics:
                unsaturated_node.append(n)
    return unsaturated_node


def find_unsaturated_linker(eG, linker_topics):
    # find unsaturated linker in eG
    unsaturated_linker = []
    for n in eG.nodes():
        if pname(n) == "EDGE" and len(list(eG.neighbors(n))) < linker_topics:
            unsaturated_linker.append(n)
    return unsaturated_linker


#######below are from _supercell.py####################
def Carte_points_generator(xyz_num):
    x_num, y_num, z_num = xyz_num
    """this function is to generate a group of 3d points(unit=1) defined by user for further grouping points"""
    unit_dx, unit_dy, unit_dz = np.array(
        [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    )
    # add x layer
    points = np.array([0.0, 0.0, 0.0])
    for i in range(0, x_num + 1):
        points = np.vstack((points, i * unit_dx))
    # add y layer
    points_x = points
    for i in range(0, y_num + 1):
        points = np.vstack((points, points_x + i * unit_dy))
    # add z layer
    points_xy = points
    for i in range(0, z_num + 1):
        points = np.vstack((points, points_xy + i * unit_dz))
    points = np.unique(points, axis=0)
    return points


#####below are from _cluster.py####################


def split_diffs(diffs):
    diff_ele = []
    diff_ele_append = diff_ele.append
    dx = np.array([1, 0, 0])
    dy = np.array([0, 1, 0])
    dz = np.array([0, 0, 1])

    for i in diffs:
        x, y, z = i
        if x != 0:
            diffx = (x * dx).tolist()
            if diffx not in diff_ele:
                diff_ele_append(diffx)
        if y != 0:
            diffy = (y * dy).tolist()
            if diffy not in diff_ele:
                diff_ele_append(diffy)
        if z != 0:
            diffz = (z * dz).tolist()
            if diffz not in diff_ele:
                diff_ele_append(diffz)
        if x * y != 0:
            diffxy = (x * dx + y * dy).tolist()
            if diffxy not in diff_ele:
                diff_ele_append(diffxy)
        if x * z != 0:
            diffxz = (x * dx + z * dz).tolist()
            if diffxz not in diff_ele:
                diff_ele_append(diffxz)
        if y * z != 0:
            diffyz = (y * dy + z * dz).tolist()
            if diffyz not in diff_ele:
                diff_ele_append(diffyz)
        if x * y * z != 0:
            diffxyz = (x * dx + y * dy + z * dz).tolist()
            if diffxyz not in diff_ele:
                diff_ele_append(diffxyz)
    return diff_ele


#######below are from makesuperG.py####################

# find DV cooresponding to V and update the sG, remove dvnode and add vnode with diff_e


def replace_sG_dvnode_with_vnode(sG, diff_e, dvnode, vnode):
    # make name

    sG.add_node(
        vnode + "_" + str(diff_e),
        f_points=np.hstack(
            (
                sG.nodes[vnode]["f_points"][:, 0:2],
                sG.nodes[vnode]["f_points"][:, 2:5].astype(float) + diff_e,
            )
        ),  # NOTE:modified because of extra column of atom type
        fcoords=sG.nodes[vnode]["fcoords"] + diff_e,
        type="SV",
        note=sG.nodes[vnode]["note"],
    )
    # process edge linked to dvnode and relink them to vnode_diff_e
    edges = [e for e in sG.edges(dvnode)]
    for e in edges:
        sG.add_edge(
            vnode + "_" + str(diff_e),
            e[1],
            fcoords=sG.edges[e]["fcoords"],
            f_points=sG.edges[e]["f_points"],
            fc_center=sG.edges[e]["fc_center"],
            type=sG.edges[e]["type"],
        )
        sG.remove_edge(e[0], e[1])
    sG.remove_node(dvnode)
    return sG


# find DV cooresponding to V and update the sG, remove dvnode and add vnode with diff_e
def replace_DV_with_corresponding_V(sG):
    # use distance matrix for moded_DV and V
    moded_dv_fc = []
    v_fc = []
    for n in sG.nodes():
        if sG.nodes[n]["type"] == "DV":
            moded_dv_fc.append((n, np.mod(sG.nodes[n]["fcoords"], 1)))
        else:
            v_fc.append((n, sG.nodes[n]["fcoords"]))

    dist_matrix = np.zeros((len(moded_dv_fc), len(v_fc)))
    for i in range(len(moded_dv_fc)):
        for j in range(len(v_fc)):
            dist_matrix[i, j] = np.linalg.norm(moded_dv_fc[i][1] - v_fc[j][1])

    # check the distance is less than 1e-2
    dv_v_pairs = []
    for k in range(len(moded_dv_fc)):
        if np.min(dist_matrix[k, :]) < 1e-2:
            corresponding_v = np.argmin(dist_matrix[k, :])
            dvnode = moded_dv_fc[k][0]
            vnode = v_fc[corresponding_v][0]
            diff_e = sG.nodes[dvnode]["fcoords"] - sG.nodes[vnode]["fcoords"]
            dv_v_pairs.append((dvnode, vnode + "_" + str(diff_e)))
            sG = replace_sG_dvnode_with_vnode(
                sG, diff_e, moded_dv_fc[k][0], v_fc[corresponding_v][0]
            )
        else:
            print(
                "distance is larger than 1e-2",
                np.min(dist_matrix[k, :]),
                moded_dv_fc[k][0],
                v_fc[np.argmin(dist_matrix[k, :])][0],
            )
    return dv_v_pairs, sG


# check if node is at the boundary of the supercell
# if at boundary, then add the diff_e to the node name and add the diff_e to the f_points
def update_supercell_node_fpoints_loose(sG, supercell):
    # boundary_node_res = []
    # incell_node_res = []
    superG = nx.Graph()
    for n in sG.nodes():
        if sG.nodes[n]["type"] == "SV":  # get rid of SV, sv will be pnode+i
            superG.add_node(
                n,
                f_points=sG.nodes[n]["f_points"],
                fcoords=sG.nodes[n]["fcoords"],
                type="SV",
                note=sG.nodes[n]["note"],
            )

            continue

        # add the node to superG, if lname(n) is np.array([0,0,0])
        superG.add_node(
            pname(n) + "_" + str(lname(n)),
            f_points=sG.nodes[n]["f_points"],
            fcoords=sG.nodes[n]["fcoords"],
            type="V",
            note=sG.nodes[n]["note"],
        )

        arr = (
            sG.nodes[n]["f_points"][:, 2:5].astype(float) + supercell
        )  # NOTE:modified because of extra column of atom type
        moded_arr = np.mod(arr, 1)
        arr = arr.astype(float)
        moded_arr = moded_arr.astype(float)
        row_diff = [i for i in range(len(arr)) if not np.allclose(arr[i], moded_arr[i])]
        diff = [arr[i] - moded_arr[i] for i in row_diff]

        diffs = np.round(np.unique(diff, axis=0), 1)
        diff_ele = split_diffs(diffs)
        # if len(diff_ele) > supercell_Carte.shape[0]:
        #    boundary_node_res.append(n)
        # else:
        #    incell_node_res.append(n)

        for diff_e in diff_ele:
            diff_e = np.asarray(diff_e)
            if (pname(n) + "_" + str(lname(n) + diff_e)) in superG.nodes():
                # print("node already in superG", pname(n) + "_" + str(lname(n) + diff_e))
                continue
            superG.add_node(
                (pname(n) + "_" + str(lname(n) + diff_e)),
                f_points=np.hstack(
                    (
                        sG.nodes[n]["f_points"][:, 0:2],
                        sG.nodes[n]["f_points"][:, 2:5].astype(float) + diff_e,
                    )
                ),  # NOTE:modified because of extra column of atom type
                fcoords=sG.nodes[n]["fcoords"] + diff_e,
                type="SV",
                note=sG.nodes[n]["note"],
            )

    return superG  # ,boundary_node_res,incell_node_res


"""
sG_node_note_set:{'CV', 'V'}
sG_node_type_set:{'DV', 'V'}
sG_edge_type_set:{'DE', 'E'}
"""


# check if edge is at the boundary of the supercell
def update_supercell_edge_fpoints(sG, superG, supercell):
    # boundary edge is DE
    # incell edge is E
    supercell_Carte = Carte_points_generator(supercell)
    for e in sG.edges():
        for i in supercell_Carte:
            s_edge = (
                pname(e[0]) + "_" + str(i + lname(e[0])),
                pname(e[1]) + "_" + str(i + lname(e[1])),
            )

            # check if node e[0]+'_'+str(diff_e) and e[1]+'_'+str(diff_e) in superG
            if (s_edge[0] in superG.nodes()) and (s_edge[1] in superG.nodes()):
                superG.add_edge(
                    s_edge[0],
                    s_edge[1],
                    f_points=np.hstack(
                        (
                            sG.edges[e]["f_points"][:, 0:2],
                            sG.edges[e]["f_points"][:, 2:5].astype(float) + i,
                        )
                    ),  # NOTE:modified because of extra column of atom type
                    fcoords=sG.edges[e]["fcoords"] + i,
                    type=sG.edges[e]["type"],
                )

            elif (s_edge[0] in superG.nodes()) or (s_edge[1] in superG.nodes()):
                if s_edge[0] in superG.nodes():
                    superG.add_node(
                        s_edge[1],
                        f_points=np.hstack(
                            (
                                sG.nodes[e[1]]["f_points"][:, 0:2],
                                sG.nodes[e[1]]["f_points"][:, 2:5].astype(float) + i,
                            )
                        ),  # NOTE:modified because of extra column of atom type
                        fcoords=sG.nodes[e[1]]["fcoords"] + i,
                        type="DSV",
                        note=sG.nodes[e[1]]["note"],
                    )

                else:
                    superG.add_node(
                        s_edge[0],
                        f_points=np.hstack(
                            (
                                sG.nodes[e[0]]["f_points"][:, 0:2],
                                sG.nodes[e[0]]["f_points"][:, 2:5].astype(float) + i,
                            )
                        ),  # NOTE:modified because of extra column of atom type
                        fcoords=sG.nodes[e[0]]["fcoords"] + i,
                        type="DSV",
                        note=sG.nodes[e[0]]["note"],
                    )

                superG.add_edge(
                    s_edge[0],
                    s_edge[1],
                    f_points=np.hstack(
                        (
                            sG.edges[e]["f_points"][:, 0:2],
                            sG.edges[e]["f_points"][:, 2:5].astype(float) + i,
                        )
                    ),  # NOTE:modified because of extra column of atom type
                    fcoords=sG.edges[e]["fcoords"] + i,
                    type="DSE",
                )

    return superG


########## the below is to process the multiedge bundling in superG###########
# need to combine with the multiedge_bundling.py
# replace bundle dvnode with vnode+diff_e
def replace_bundle_dvnode_with_vnode(dv_v_pairs, multiedge_bundlings):
    for dv, v in dv_v_pairs:
        for bund in multiedge_bundlings:
            if dv in bund[1]:
                if pname(bund[1][bund[1].index(dv)]) == pname(dv):
                    bund[1][bund[1].index(dv)] = v
            if dv in bund[0]:
                if pname(bund[0][bund[0].index(dv)]) == pname(dv):
                    bund[0][bund[0].index(dv)] = v
    # update v if no list then add [0,0,0]
    # convert tuple to list
    updated_bundlings = []
    for bund in multiedge_bundlings:
        ec_node = pname(bund[0]) + "_" + str(lname(bund[0]))
        con_nodes = [pname(i) + "_" + str(lname(i)) for i in bund[1]]
        updated_bundlings.append((ec_node, con_nodes))
    return updated_bundlings


# loop bundle and check if any element in the bundle is in the superG, if not, add the element to the superG
def make_super_multiedge_bundlings(prim_multiedge_bundlings, supercell):
    super_multiedge_bundlings = {}
    for i in Carte_points_generator(supercell):
        for bund in prim_multiedge_bundlings:
            ec_node = pname(bund[0]) + "_" + str(i + lname(bund[0]))
            con_nodes = [pname(n) + "_" + str(i + lname(n)) for n in bund[1]]
            super_multiedge_bundlings[ec_node] = con_nodes
    return super_multiedge_bundlings


def update_supercell_bundle(superG, super_multiedge_bundlings):
    # print("*" * 50)
    # print("superG nodes name", list(superG.nodes()))
    # print("*" * 50)
    for ec_node in super_multiedge_bundlings.keys():
        con_nodes = super_multiedge_bundlings[ec_node]
        # order the con_nodes by th x-x pair of the ecnode X atoms
        prim_ecname = pname(ec_node) + "_" + str(np.array([0.0, 0.0, 0.0]))
        if ec_node not in superG.nodes():
            trans = lname(ec_node)
            superG.add_node(
                ec_node,
                f_points=np.hstack(
                    (
                        superG.nodes[prim_ecname]["f_points"][
                            :, 0:2
                        ],  # NOTE:modified because of extra column of atom type
                        superG.nodes[prim_ecname]["f_points"][:, 2:5].astype(float)
                        + trans,
                    )
                ),  # NOTE:modified because of extra column of atom type
                fcoords=superG.nodes[prim_ecname]["fcoords"] + trans,
                type="SV",
                note=superG.nodes[prim_ecname]["note"],
            )
        for j in range(len(con_nodes)):
            cn = con_nodes[j]
            prim_cnname = super_multiedge_bundlings[
                prim_ecname
            ][
                j
            ]  # find prim_ecname in super_multiedge_bundlings and then get the corresponding prim_cnname
            trans = lname(cn) - lname(prim_cnname)
            if cn not in superG.nodes():
                superG.add_node(
                    cn,
                    f_points=np.hstack(
                        (
                            superG.nodes[prim_cnname]["f_points"][
                                :, 0:2
                            ],  # NOTE:modified because of extra column of atom type
                            superG.nodes[prim_cnname]["f_points"][:, 2:5].astype(float)
                            + trans,
                        )
                    ),  # NOTE:modified because of extra column of atom type
                    fcoords=superG.nodes[prim_cnname]["fcoords"] + trans,
                    type="SV",
                    note=superG.nodes[prim_cnname]["note"],
                )
            superG.add_edge(
                ec_node,
                cn,
                f_points=np.hstack(
                    (
                        superG.edges[prim_ecname, prim_cnname]["f_points"][
                            :, 0:2
                        ],  # NOTE:modified because of extra column of atom type
                        superG.edges[prim_ecname, prim_cnname]["f_points"][
                            :, 2:5
                        ].astype(float)
                        + trans,
                    )
                ),  # NOTE:modified because of extra column of atom type
                fcoords=superG.edges[prim_ecname, prim_cnname]["fcoords"] + trans,
                type="DSE",
            )

    return superG


def check_multiedge_bundlings_insuperG(super_multiedge_bundlings, superG):
    super_multiedge_bundlings_edges = []
    for ec_node in super_multiedge_bundlings:
        # check is all CV node in superG are in the super_multiedge_bundlings_edges first element
        cvnodes = [n for n in superG.nodes() if superG.nodes[n]["note"] == "CV"]
        # use set to check if all cvnodes are in the super_multiedge_bundlings_edges
        if set(cvnodes) == set([i[0] for i in super_multiedge_bundlings_edges]):
            return superG
        else:
            # print("not all CV nodes in super_multiedge_bundlings_edges")
            diff_element = set(cvnodes).difference(set(list(super_multiedge_bundlings)))
            # print("to remove diff_element", diff_element)
            # remove the diff_element from the superG
            for n in diff_element:
                superG.remove_node(n)
                # remove all edges linked to the node
                edges = [e for e in superG.edges(n)]
                for e in edges:
                    superG.remove_edge(e[0], e[1])

            return superG


########## the above is to process the multiedge bundling in superG###########


def locate_min_idx(matrix):
    min_idx = np.unravel_index(matrix.argmin(), matrix.shape)
    return min_idx[0], min_idx[1]


def add_virtual_edge(superG, bridge_node_distance, max_neighbor=2):
    # add pillar nodes virtual edges
    nodes_list = [n for n in superG.nodes() if superG.nodes[n]["note"] == "V"]
    n_n_distance_matrix = np.zeros((len(nodes_list), len(nodes_list)))

    for i in range(len(nodes_list)):
        for j in range(len(nodes_list)):
            n_n_distance_matrix[i, j] = np.linalg.norm(
                superG.nodes[nodes_list[i]]["fcoords"]
                - superG.nodes[nodes_list[j]]["fcoords"]
            )
        n_n_distance_matrix[i, i] = 1000
    # use hungrain algorithm to find the shortest path between all nodes

    for i in range(len(nodes_list)):
        neighbor_count = 0
        while neighbor_count < max_neighbor:

            def add_virtual_edge(
                i, n_n_distance_matrix, superG, bridge_node_distance, count
            ):
                n_n_min_distance = np.min(n_n_distance_matrix[i : i + 1, :])
                if n_n_min_distance < bridge_node_distance:
                    _, n_j = locate_min_idx(n_n_distance_matrix[i : i + 1, :])
                    superG.add_edge(nodes_list[i], nodes_list[n_j], type="virtual")
                    # print('add virtual edge between',nodes_list[i],nodes_list[n_j])
                    n_n_distance_matrix[i, n_j] = 1000
                    return True, count + 1, n_n_distance_matrix, superG
                else:
                    return False, count, n_n_distance_matrix, superG

            added, neighbor_count, n_n_distance_matrix, superG = add_virtual_edge(
                i, n_n_distance_matrix, superG, bridge_node_distance, neighbor_count
            )
            if not added:
                break

    return superG


######below are from v2_functions.py####################
# Function to fetch indices and coordinates of atoms with a specific label
def fetch_X_atoms_ind_array(array, column, X):
    # array: input array
    # column: column index to check for label
    # X: label to search for

    ind = [k for k in range(len(array)) if re.sub(r"\d", "", array[k, column]) == X]
    x_array = array[ind]
    return ind, x_array


def sort_nodes_by_type_connectivity(G):
    CV_nodes = [n for n in G.nodes() if G.nodes[n]["note"] == "CV"]
    if len(CV_nodes) == 0:  # ditopic linker MOF
        Vnodes = [n for n in G.nodes() if G.nodes[n]["type"] == "V"]
        DVnodes = [n for n in G.nodes() if G.nodes[n]["type"] == "DV"]
        Vnodes = sorted(Vnodes, key=lambda x: G.degree(x), reverse=True)
        DVnodes = sorted(DVnodes, key=lambda x: G.degree(x), reverse=True)
        return Vnodes + DVnodes
    else:
        # CV+V
        # get CV_Vnode
        CV_Vnodes = [
            n
            for n in G.nodes()
            if G.nodes[n]["type"] == "V" and G.nodes[n]["note"] == "CV"
        ]
        CV_DVnodes = [
            n
            for n in G.nodes()
            if G.nodes[n]["type"] == "DV" and G.nodes[n]["note"] == "CV"
        ]
        V_Vnodes = [
            n
            for n in G.nodes()
            if G.nodes[n]["type"] == "V" and G.nodes[n]["note"] == "V"
        ]
        V_DVnodes = [
            n
            for n in G.nodes()
            if G.nodes[n]["type"] == "DV" and G.nodes[n]["note"] == "V"
        ]
        CV_Vnodes = sorted(CV_Vnodes, key=lambda x: G.degree(x), reverse=True)
        CV_DVnodes = sorted(CV_DVnodes, key=lambda x: G.degree(x), reverse=True)
        V_Vnodes = sorted(V_Vnodes, key=lambda x: G.degree(x), reverse=True)
        V_DVnodes = sorted(V_DVnodes, key=lambda x: G.degree(x), reverse=True)

        return V_Vnodes + CV_Vnodes + V_DVnodes + CV_DVnodes


def check_edge_inunitcell(G, e):
    if "DV" in G.nodes[e[0]]["type"] or "DV" in G.nodes[e[1]]["type"]:
        return False
    return True


def put_V_ahead_of_CV(e):
    if "CV" in e[0] and "V" in e[1]:
        return (e[1], e[0])
    else:
        return e


def find_and_sort_edges_bynodeconnectivity(graph, sorted_nodes):
    all_edges = list(graph.edges())

    sorted_edges = []
    # add unit_cell edge first

    ei = 0
    while ei < len(all_edges):
        e = all_edges[ei]
        if check_edge_inunitcell(graph, e):
            sorted_edges.append(put_V_ahead_of_CV(e))
            all_edges.pop(ei)
        ei += 1
    # sort edge by sorted_nodes
    for n in sorted_nodes:
        ei = 0
        while ei < len(all_edges):
            e = all_edges[ei]
            if n in e:
                if n == e[0]:
                    sorted_edges.append(put_V_ahead_of_CV(e))
                else:
                    sorted_edges.append(put_V_ahead_of_CV((e[1], e[0])))
                all_edges.pop(ei)
            else:
                ei += 1

    return sorted_edges


def make_unsaturated_vnode_xoo_dict(
    unsaturated_node, xoo_dict, matched_vnode_xind, eG, sc_unit_cell
):
    """
    make a dictionary of the unsaturated node and the exposed X connected atom index and the corresponding O connected atoms
    """

    # process matched_vnode_xind make it to a dictionary
    matched_vnode_xind_dict = {}
    for [k, v, e] in matched_vnode_xind:
        if k in matched_vnode_xind_dict.keys():
            matched_vnode_xind_dict[k].append(v)
        else:
            matched_vnode_xind_dict[k] = [v]

    unsaturated_vnode_xind_dict = {}
    xoo_keys = list(xoo_dict.keys())
    # for each unsaturated node, get the upmatched x index and xoo atoms
    for unsat_v in unsaturated_node:
        if unsat_v in matched_vnode_xind_dict.keys():
            unsaturated_vnode_xind_dict[unsat_v] = [
                i for i in xoo_keys if i not in matched_vnode_xind_dict[unsat_v]
            ]
            # print(unsaturated_vnode_xind_dict[unsat_v],'unsaturated_vnode_xind_dict[unsat_v]') #DEBUG
        else:
            unsaturated_vnode_xind_dict[unsat_v] = xoo_keys

    # based on the unsaturated_vnode_xind_dict, add termination to the unsaturated node xoo
    # loop over unsaturated nodes, and find all exposed X atoms and use paied xoo atoms to form a termination
    unsaturated_vnode_xoo_dict = {}
    for vnode, exposed_x_indices in unsaturated_vnode_xind_dict.items():
        for xind in exposed_x_indices:
            x_fpoints = eG.nodes[vnode]["f_points"][xind]
            x_cpoints = np.hstack(
                (x_fpoints[0:2], fractional_to_cartesian(x_fpoints[2:5], sc_unit_cell))
            )  # NOTE: modified add the atom type and atom name
            oo_ind_in_vnode = xoo_dict[xind]
            oo_fpoints_in_vnode = [
                eG.nodes[vnode]["f_points"][i] for i in oo_ind_in_vnode
            ]
            oo_fpoints_in_vnode = np.vstack(oo_fpoints_in_vnode)
            oo_cpoints = np.hstack(
                (
                    oo_fpoints_in_vnode[:, 0:2],
                    fractional_to_cartesian(oo_fpoints_in_vnode[:, 2:5], sc_unit_cell),
                )
            )  # NOTE: modified add the atom type and atom name

            unsaturated_vnode_xoo_dict[(vnode, xind)] = {
                "xind": xind,
                "oo_ind": oo_ind_in_vnode,
                "x_fpoints": x_fpoints,
                "x_cpoints": x_cpoints,
                "oo_fpoints": oo_fpoints_in_vnode,
                "oo_cpoints": oo_cpoints,
            }

    return (
        unsaturated_vnode_xind_dict,
        unsaturated_vnode_xoo_dict,
        matched_vnode_xind_dict,
    )


def update_matched_nodes_xind(
    removed_nodes_list, removed_edges_list, matched_vnode_xind
):
    # if linked edge is removed and the connected node is not removed, then remove this line from matched_vnode_xind
    # add remove the middle xind of the node to matched_vnode_xind_dict[node] list
    to_remove_row = []

    for i in range(len(matched_vnode_xind)):
        node, xind, edge = matched_vnode_xind[i]
        if edge in removed_edges_list and node not in removed_nodes_list:
            # print("remove edge", edge, matched_vnode_xind[i], "from matched_vnode_xind")
            to_remove_row.append(i)
        elif node in removed_nodes_list:
            to_remove_row.append(i)
            # print("remove node", node, matched_vnode_xind[i], "from matched_vnode_xind")

    # remove the rows

    update_matched_vnode_xind = [
        i for j, i in enumerate(matched_vnode_xind) if j not in to_remove_row
    ]
    return update_matched_vnode_xind


# functions for write
# write gro file
def extract_node_edge_term(tG, sc_unit_cell):
    nodes_tG = []
    terms_tG = []
    edges_tG = []
    node_res_num = 0
    term_res_num = 0
    edge_res_num = 0
    nodes_check_set = set()
    nodes_name_set = set()
    edges_check_set = set()
    edges_name_set = set()
    terms_check_set = set()
    terms_name_set = set()
    for n in tG.nodes():
        if pname(n) != "EDGE":
            postions = tG.nodes[n]["noxoo_f_points"]
            name = tG.nodes[n]["name"]
            nodes_check_set.add(len(postions))
            nodes_name_set.add(name)
            if len(nodes_check_set) > len(nodes_name_set):
                raise ValueError("node index is not continuous")
            node_res_num += 1
            nodes_tG.append(
                np.hstack(
                    (
                        np.tile(
                            np.array([node_res_num, name]), (len(postions), 1)
                        ),  # residue number and residue name
                        postions[:, 1:2],  # atom type (element)
                        fractional_to_cartesian(
                            postions[:, 2:5], sc_unit_cell
                        ),  # Cartesian coordinates
                        postions[:, 0:1],  # atom name
                        np.tile(np.array([n]), (len(postions), 1)),
                    )
                )
            )  # node name in eG is added to the last column
            if "term_c_points" in tG.nodes[n]:
                for term_ind_key, c_positions in tG.nodes[n]["term_c_points"].items():
                    terms_check_set.add(len(c_positions))
                    name = "T" + tG.nodes[n]["name"]
                    terms_name_set.add(name)
                    if len(terms_check_set) > len(terms_name_set):
                        raise ValueError("term index is not continuous")

                    term_res_num += 1
                    terms_tG.append(
                        np.hstack(
                            (
                                np.tile(
                                    np.array([term_res_num, name]),
                                    (len(c_positions), 1),
                                ),  # residue number and residue name
                                c_positions[:, 1:2],  # atom type (element)
                                c_positions[:, 2:5],  # Cartesian coordinates
                                c_positions[:, 0:1],  # atom name
                                np.tile(
                                    np.array([term_ind_key]), (len(c_positions), 1)
                                ),
                            )
                        )
                    )  # term name in eG is added to the last column

        elif pname(n) == "EDGE":
            postions = tG.nodes[n]["f_points"]
            name = tG.nodes[n]["name"]
            edges_check_set.add(len(postions))
            edges_name_set.add(name)
            if len(edges_check_set) > len(edges_name_set):
                print(edges_check_set)
                # raise ValueError('edge atom number is not continuous')
                print(
                    "edge atom number is not continuous,ERROR edge name:",
                    len(edges_check_set),
                    len(edges_name_set),
                )
            edge_res_num += 1
            edges_tG.append(
                np.hstack(
                    (
                        np.tile(
                            np.array([edge_res_num, name]), (len(postions), 1)
                        ),  # residue number and residue name
                        postions[:, 1:2],  # atom type (element)
                        fractional_to_cartesian(
                            postions[:, 2:5], sc_unit_cell
                        ),  # Cartesian coordinates
                        postions[:, 0:1],  # atom name
                        np.tile(np.array([n]), (len(postions), 1)),
                    )
                )
            )  # edge name in eG is added to the last column

    # nodes_tG = np.vstack(nodes_tG)
    # terms_tG = np.vstack(terms_tG)
    # edges_tG = np.vstack(edges_tG)
    return nodes_tG, edges_tG, terms_tG, node_res_num, edge_res_num, term_res_num


def check_supercell_box_range(point, supercell, buffer_plus, buffer_minus):
    # to cleave eG to supercell box

    supercell_x = supercell[0] + buffer_plus
    supercell_y = supercell[1] + buffer_plus
    supercell_z = supercell[2] + buffer_plus
    if (
        point[0] >= 0 + buffer_minus
        and point[0] <= supercell_x
        and point[1] >= 0 + buffer_minus
        and point[1] <= supercell_y
        and point[2] >= 0 + buffer_minus
        and point[2] <= supercell_z
    ):
        return True
    else:
        # print(point, 'out of supercell box range:  [',supercell_x,supercell_y,supercell_z, '],   will be excluded') #debug
        return False


def replace_edges_by_callname(
    edge_n_list, eG, sc_unit_cell_inv, new_linker_pdb, prefix="R"
):
    new_linker_atoms, new_linker_ccoords, new_linker_x_ccoords = process_node_pdb(
        new_linker_pdb, "X"
    )
    for edge_n in edge_n_list:
        edge_n = edge_n
        edge_f_points = eG.nodes[edge_n]["f_points"]
        x_indices = [
            i for i in range(len(edge_f_points)) if nn(edge_f_points[i][0]) == "X"
        ]
        edge_x_points = edge_f_points[x_indices]
        edge_com = np.mean(edge_x_points[:, 2:5].astype(float), axis=0)
        edge_x_fcoords = edge_x_points[:, 2:5].astype(float) - edge_com

        new_linker_x_fcoords = cartesian_to_fractional(
            new_linker_x_ccoords, sc_unit_cell_inv
        )
        new_linker_fcoords = cartesian_to_fractional(
            new_linker_ccoords, sc_unit_cell_inv
        )

        _, rot, trans = superimpose(new_linker_x_fcoords, edge_x_fcoords)
        replaced_linker_fcoords = np.dot(new_linker_fcoords, rot) + edge_com
        replaced_linker_f_points = np.hstack(
            (new_linker_atoms, replaced_linker_fcoords)
        )

        eG.nodes[edge_n]["f_points"] = replaced_linker_f_points
        eG.nodes[edge_n]["name"] = prefix + eG.nodes[edge_n]["name"]

    return eG


# the following functions are used for the split node to metal, hho,ho,o and update name and residue number


def extract_node_name_from_gro_resindex(res_index, node_array_list):
    node_array = np.vstack(node_array_list)
    nodes_name = set()
    for node_ind in res_index:
        node_name = node_array[node_array[:, 0] == str(node_ind)][:, -1]
        name_set = set(node_name)
        nodes_name = nodes_name.union(name_set)
    return nodes_name


def make_dummy_split_node_dict(dummy_node_name):
    node_split_dict = {}
    dict_path = dummy_node_name.split(".")[0] + "_dict"
    with open(dict_path, "r") as f:
        lines = f.readlines()
    node_res_counts = 0
    for li in lines:
        li = li.strip("\n")
        key = li[:20].strip(" ")
        value = li[-4:].strip(" ")
        node_split_dict[key] = int(value)
    return node_split_dict


def chunk_array(chunk_list, array, chunk_num, chunksize):
    chunk_list.extend(
        array[i * chunksize : (i + 1) * chunksize] for i in range(chunk_num)
    )
    return chunk_list


def rename_node_arr(node_split_dict, node_arr):
    metal_count = node_split_dict["METAL_count"]
    dummy_len = int(node_split_dict["dummy_res_len"])
    metal_num = metal_count * dummy_len
    hho_num = node_split_dict["HHO_count"] * 3
    ho_num = node_split_dict["HO_count"] * 2
    o_num = node_split_dict["O_count"] * 1
    metal_range = metal_num
    hho_range = metal_range + hho_num
    ho_range = hho_range + ho_num
    o_range = ho_range + o_num
    # print(metal_range,hho_range,ho_range,o_range) #debug

    metals_list = []
    hhos_list = []
    hos_list = []
    os_list = []
    for idx in set(node_arr[:, 0]):
        idx_arr = node_arr[node_arr[:, 0] == idx]
        if metal_num > 0:
            metal = idx_arr[0:metal_range].copy()
            metal[:, 1] = "METAL"
            metals_list = chunk_array(
                metals_list, metal, node_split_dict["METAL_count"], dummy_len
            )
        if hho_num > 0:
            hho = idx_arr[metal_range:hho_range].copy()
            hho[:, 1] = "HHO"
            hhos_list = chunk_array(hhos_list, hho, node_split_dict["HHO_count"], 3)
        if ho_num > 0:
            ho = idx_arr[hho_range:ho_range].copy()
            ho[:, 1] = "HO"
            hos_list = chunk_array(hos_list, ho, node_split_dict["HO_count"], 2)
        if o_num > 0:
            o = idx_arr[ho_range:o_range].copy()
            o[:, 1] = "O"
            os_list = chunk_array(os_list, o, node_split_dict["O_count"], 1)

    return metals_list, hhos_list, hos_list, os_list


def merge_metal_list_to_node_array(
    merged_node_edge_term, metals_list, line_num, res_count
):
    if any([len(metal) == 0 for metal in metals_list]):
        return merged_node_edge_term, line_num, res_count
    for i in range(len(metals_list)):
        metal = metals_list[i]
        metal[:, 0] = i + 1
        formatted_gro_lines, line_num = convert_node_array_to_gro_lines(
            metal, line_num, res_count
        )
        merged_node_edge_term += formatted_gro_lines
    res_count += len(metals_list)
    return merged_node_edge_term, line_num, res_count


def convert_node_array_to_gro_lines(array, line_num_start, res_num_start):
    formatted_gro_lines = []

    for i in range(len(array)):
        line = array[i]
        ind_inres = i + 1
        name = line[1]
        value_atom_number_in_gro = int(ind_inres + line_num_start)  # atom_number
        value_label = re.sub("\d", "", line[2]) + str(ind_inres)  # atom_label
        value_resname = str(name)[0:3]  # +str(eG.nodes[n]['index'])  # residue_name
        value_resnumber = int(res_num_start + int(line[0]))  # residue number
        value_x = 0.1 * float(line[3])  # x
        value_y = 0.1 * float(line[4])  # y
        value_z = 0.1 * float(line[5])  # z
        formatted_line = "%5d%-5s%5s%5d%8.3f%8.3f%8.3f" % (
            value_resnumber,
            value_resname,
            value_label,
            value_atom_number_in_gro,
            value_x,
            value_y,
            value_z,
        )
        formatted_gro_lines.append(formatted_line + "\n")
    return formatted_gro_lines, value_atom_number_in_gro


def merge_node_edge_term(nodes_tG, edges_tG, terms_tG, node_res_num, edge_res_num):
    merged_node_edge_term = []
    line_num = 0
    for node in nodes_tG:
        formatted_gro_lines, line_num = convert_node_array_to_gro_lines(
            node, line_num, 0
        )
        merged_node_edge_term += formatted_gro_lines
    for edge in edges_tG:
        formatted_gro_lines, line_num = convert_node_array_to_gro_lines(
            edge, line_num, node_res_num
        )
        merged_node_edge_term += formatted_gro_lines
    for term in terms_tG:
        formatted_gro_lines, line_num = convert_node_array_to_gro_lines(
            term, line_num, node_res_num + edge_res_num
        )
        merged_node_edge_term += formatted_gro_lines
    return merged_node_edge_term


def save_node_edge_term_gro(merged_node_edge_term, gro_name, dir_name="output_gros"):
    Path(dir_name).mkdir(parents=True, exist_ok=True)
    gro_name = str(Path(dir_name, gro_name))
    with open(gro_name + ".gro", "w") as f:
        head = []
        head.append("eG_NET\n")
        head.append(str(len(merged_node_edge_term)) + "\n")
        f.writelines(head)
        f.writelines(merged_node_edge_term)
        tail = ["20 20 20 \n"]
        f.writelines(tail)


#################below are from display.py######################
def gro_show(gro_file, w=800, h=600, res_id=True, res_name=True):
    try:
        import py3Dmol

        viewer = py3Dmol.view(width=w, height=h)
        with open(gro_file, "r") as f:
            lines = f.readlines()

        viewer.addModel("".join(lines), "gro")
        # viewer.setStyle({"stick": {}})

        viewer.setViewStyle({"style": "outline", "width": 0.05})
        viewer.setStyle({"stick": {}, "sphere": {"scale": 0.20}})
        if res_id or res_name:
            for i in range(2, len(lines) - 1):
                if lines[i].strip() == "":
                    continue
                if lines[i - 1][0:5] == lines[i][0:5]:
                    continue

                value_resnumber = int((lines[i])[0:5])
                value_resname = lines[i][5:10]
                if value_resname.strip() == "TNO":
                    continue
                # value_label = lines[i][10:15]
                # value_atom_number = int(lines[i][15:20])
                value_x = float(lines[i][20:28]) * 10  # x
                value_y = float(lines[i][28:36]) * 10  # y
                value_z = float(lines[i][36:44]) * 10  # z

                text = ""
                if res_name:
                    text += str(value_resname)
                if res_id:
                    text += str(value_resnumber)

                viewer.addLabel(
                    text,
                    {
                        "position": {
                            "x": value_x,
                            "y": value_y,
                            "z": value_z,
                        },
                        "alignment": "center",
                        "fontColor": "white",
                        "font": "Arial",
                        "fontSize": 12,
                        "backgroundColor": "black",
                        "backgroundOpacity": 0.5,
                    },
                )
        viewer.render()
        viewer.zoomTo()
        viewer.show()
    except ImportError:
        raise ImportError("Unable to import py3Dmol")


class NetOptimizer:
    """
    NetOptimizer is a class to optimize the node and edge structure of the MOF, add terminations to nodes.

    :param template_cif (str):
        cif file of the template, including only V and E *(EC)nodes info in primitive cell
    Instance variables:
        - node_cif (str):cif file of the node
        - node_target_type (str):metal atom type of the node
        - node_unit_cell (array):unit cell of the node
        - node_atom (array):2 columns, atom_name, atom_type of the node
        - node_x_fcoords (array):fractional coordinates of the X connected atoms of node
        - node_fcoords (array):fractional coordinates of the whole node
        - node_x_ccoords (array):cartesian coordinates of the X connected atoms of node
        - node_coords (array):cartesian coordinates of the whole node
        - linker_cif (str):cif file of the ditopic linker or branch of multitopic linker
        - linker_unit_cell (array):unit cell of the ditopic linker or branch of multitopic linker
        - linker_atom (array):2 columns, atom_name, atom_type of the ditopic linker or branch of multitopic linker
        - linker_x_fcoords (array):fractional coordinates of the X connected atoms of ditopic linker or branch of multitopic linker
        - linker_fcoords (array):fractional coordinates of the whole ditopic linker or branch of multitopic linker
        - linker_x_ccoords (array):cartesian coordinates of the X connected atoms of ditopic linker or branch of multitopic linker
        - linker_length (float):distance between two X-X connected atoms of the ditopic linker or branch of multitopic linker
        - linker_ccoords (array):cartesian coordinates of the whole ditopic linker or branch of multitopic linker
        - linker_center_cif (str):cif file of the center of the multitopic linker
        - ec_unit_cell (array):unit cell of the center of the multitopic linker
        - ec_atom (array):2 columns, atom_name, atom_type of the center of the multitopic linker
        - ec_x_vecs (array):fractional coordinates of the X connected atoms of the center of the multitopic linker
        - ec_fcoords (array):fractional coordinates of the whole center of the multitopic linker
        - ec_xcoords (array):cartesian coordinates of the X connected atoms of the center of the multitopic linker
        - eccoords (array):cartesian coordinates of the whole center of the multitopic linker
        - constant_length (float):constant length to add to the linker length, normally 1.54 for single bond of C-C,
        because C is always used as the connecting atom in the builder
        - maxfun (int):maximum number of function evaluations for the node rotation optimization
        - opt_method (str):optimization method for the node rotation optimization
        - G (networkx graph):graph of the template
        - node_max_degree (int):maximum degree of the node in the template, should be the same as the node topic
        - sorted_nodes (list):sorted nodes in the template by connectivity
        - sorted_edges (list):sorted edges in the template by connectivity
        - sorted_edges_of_sortednodeidx (list):sorted edges in the template by connectivity with the index of the sorted nodes
        - optimized_rotations (dict):optimized rotations for the nodes in the template
        - optimized_params (array):optimized cell parameters for the template topology to fit the target MOF cell
        - new_edge_length (float):new edge length of the ditopic linker or branch of multitopic linker, 2*constant_length+linker_length
        - optimized_pair (dict): pair of connected nodes in the template with the index of the X connected atoms, used for the edge placement
        - scaled_rotated_node_positions (dict):scaled and rotated node positions in the target MOF cell
        - scaled_rotated_Xatoms_positions (dict):scaled and rotated X connected atom positions of nodes in the target MOF cell
        - sc_unit_cell (array):(scaled) unit cell of the target MOF cell
        - sc_unit_cell_inv (array):inverse of the (scaled) unit cell of the target MOF cell
        - sG_node (networkx graph):graph of the target MOF cell
        - nodes_atom (dict):atom name and atom type of the nodes
        - rotated_node_positions (dict):rotated node positions
        - supercell (array): supercell set by user, along x,y,z direction
        - multiedge_bundlings (dict):multiedge bundlings of center and branches of the multitopic linker, used for the supercell
        construction and merging of center and branches to form one EDGE
        - prim_multiedge_bundlings (dict):multiedge bundlings in primitive cell, used for the supercell construction
        - super_multiedge_bundlings (dict):multiedge bundlings in the supercell, used for the supercell construction
        - dv_v_pairs (dict):DV and V pairs in the template, used for the supercell construction
        - super_multiedge_bundlings (dict):multiedge bundlings in the supercell, used for the supercell construction
        - superG (networkx graph):graph of the supercell
        - add_virtual_edge (bool): add virtual edge to the target MOF cell
        - vir_edge_range (float): range to search the virtual edge between two Vnodes directly, should <= 0.5,
        used for the virtual edge addition of bridge type nodes: nodes and nodes can connect directly without linker
        - vir_edge_max_neighbor (int): maximum number of neighbors of the node with virtual edge, used for the virtual edge addition of bridge type nodes
        - remove_node_list (list):list of nodes to remove in the target MOF cell
        - remove_edge_list (list):list of edges to remove in the target MOF cell
        - eG (networkx graph):graph of the target MOF cell with only EDGE and V nodes
        - node_topic (int):maximum degree of the node in the template, should be the same as the node_max_degree
        - unsaturated_node (list):unsaturated nodes in the target MOF cell
        - term_info (array):information of the node terminations
        - term_coords (array):coordinates of the node terminations
        - term_xoovecs (array):X and O vectors (usually carboxylate group) of the node terminations
        - unsaturated_vnode_xind_dict (dict):unsaturated node and the exposed X connected atom index
        - unsaturated_vnode_xoo_dict (dict):unsaturated node and the exposed X connected atom index and the corresponding O connected atoms
    """

    def __init__(self):
        pass

    def analyze_template_multitopic(self, template_cif):
        """
        analyze the template topology of the multitopic linker

        :param vvnode333 (array):
            supercell of V nodes in template topology
        :param ecnode333 (array):
            supercell of EC nodes (Center of multitopic linker) in template topology
        :param eenode333 (array):
            supercell of E nodes(ditopic linker or branch of multitopic linker) in template
        :param unit_cell (array):
            unit cell of the template
        :param cell_info (array):
            cell information of the template
        """
        template_cell_info, _, vvnode = extract_type_atoms_fcoords_in_primitive_cell(
            template_cif, "V"
        )
        _, _, eenode = extract_type_atoms_fcoords_in_primitive_cell(template_cif, "E")
        _, _, ecnode = extract_type_atoms_fcoords_in_primitive_cell(template_cif, "EC")
        unit_cell = extract_unit_cell(template_cell_info)

        vvnode = np.unique(np.array(vvnode, dtype=float), axis=0)
        eenode = np.unique(np.array(eenode, dtype=float), axis=0)
        ecnode = np.unique(np.array(ecnode, dtype=float), axis=0)
        ##loop over super333xxnode and super333yynode to find the pair of x node in unit cell which pass through the yynode
        vvnode333 = make_supercell_3x3x3(vvnode)
        eenode333 = make_supercell_3x3x3(eenode)
        ecnode333 = make_supercell_3x3x3(ecnode)

        _, _, G = find_pair_v_e_c(vvnode333, ecnode333, eenode333, unit_cell)
        G = add_ccoords(G, unit_cell)
        G, self.node_max_degree = set_DV_V(G)
        self.G = set_DE_E(G)
        self.cell_info = template_cell_info
        self.vvnode333 = vvnode333
        self.eenode333 = eenode333
        self.ecnode333 = ecnode333

    def analyze_template_ditopic(self, template_cif):
        """
        analyze the template topology of the ditopic linker, only V and E nodes in the template

        :param template_cif (str):
            cif file of the template, including only V and E nodes info in primitive cell
        """
        template_cell_info, _, vvnode = extract_type_atoms_fcoords_in_primitive_cell(
            template_cif, "V"
        )
        _, _, eenode = extract_type_atoms_fcoords_in_primitive_cell(template_cif, "E")
        unit_cell = extract_unit_cell(template_cell_info)

        vvnode = np.unique(np.array(vvnode, dtype=float), axis=0)
        eenode = np.unique(np.array(eenode, dtype=float), axis=0)
        ##loop over super333xxnode and super333yynode to find the pair of x node in unicell which pass through the yynode
        vvnode333 = make_supercell_3x3x3(vvnode)
        eenode333 = make_supercell_3x3x3(eenode)
        _, _, G = find_pair_v_e(vvnode333, eenode333)
        G = add_ccoords(G, unit_cell)
        G, self.node_max_degree = set_DV_V(G)
        self.G = set_DE_E(G)
        self.cell_info = template_cell_info
        self.vvnode333 = vvnode333
        self.eenode333 = eenode333

    def node_info(self, node_pdb):
        """
        get the node information

        :param node_cif (str):
            cif file of the node
        :param node_target_type (str):
            metal atom type of the node
        """
        self.node_pdb = node_pdb
        self.node_atom, self.node_ccoords, self.node_x_ccoords = process_node_pdb(
            node_pdb, "X"
        )  # com type could be metal in bridge nodes

    def linker_info(self, linker_pdb):
        """
        get the linker information

        :param linker_cif (str):
            cif file of the ditopic linker or branch of multitopic linker
        """
        self.linker_pdb = linker_pdb
        self.linker_atom, self.linker_ccoords, self.linker_x_ccoords = process_node_pdb(
            linker_pdb, "X"
        )
        self.linker_length = np.linalg.norm(
            self.linker_x_ccoords[0] - self.linker_x_ccoords[1]
        )

    def linker_center_info(self, linker_center_pdb):
        """
        get the linker information

        :param linker_cif (str):
            cif file of the ditopic linker or branch of multitopic linker
        """
        self.linker_center_pdb = linker_center_pdb
        self.ec_atom, self.ec_ccoords, self.ec_x_ccoords = process_node_pdb(
            linker_center_pdb, "X"
        )

    def set_constant_length(self, constant_length):
        """
        set the constant length to add to the linker length, normally 1.54 (default setting)for single bond of C-C, because C is always used as the connecting atom in the builder
        """
        self.constant_length = constant_length

    def set_rotation_optimizer_maxfun(self, maxfun):
        """
        set the maximum number of function evaluations for the node rotation optimization
        """
        self.maxfun = maxfun

    def set_rotation_optimizer_maxiter(self, maxiter):
        """
        set the maximum number of iterations for the node rotation optimization
        """
        self.maxiter = maxiter

    def set_rotation_optimizer_method(self, opt_method):
        """
        set the optimization method for the node rotation optimization
        """
        self.opt_method = opt_method

    def set_rotation_optimizer_display(self, display):
        """
        set the display of the optimization process
        """
        self.display = display

    def set_rotation_optimizer_eps(self, eps):
        """
        set the eps of the optimization
        """
        self.eps = eps

    def set_rotation_optimizer_iprint(self, iprint):
        """
        set the iprint of the optimization
        """
        self.iprint = iprint

    def check_node_template_match(self):
        """
        precheck, check if the number of nodes in the template matches the maximum degree of the node in the template
        """
        return len(self.node_x_ccoords) == self.node_max_degree

    def load_saved_optimized_rotations(self, optimized_rotations):
        """
        use the saved optimized rotations from the previous optimization
        """
        self.saved_optimized_rotations = optimized_rotations
        print("load the saved optimized_rotations from the previous optimization")

    def to_save_optimized_rotations(self, filename):
        """
        save the optimized rotations to the file
        """
        self.to_save_optimized_rotations_filename = filename

    def use_saved_rotations_as_initial_guess(
        self, use_saved_rotations_as_initial_guess
    ):
        """
        use the saved optimized rotations as initial guess
        """
        self.use_saved_rotations_as_initial_guess = use_saved_rotations_as_initial_guess

    def optimize(self):
        """
        two optimization steps:
        1. optimize the node rotation
        2. optimize the cell parameters to fit the target MOF cell
        """
        if not self.check_node_template_match():
            print(
                "The number of nodes in the template does not match the maximum degree of the node in the template"
            )
            raise ValueError(
                "The number of nodes in the template does not match the maximum degree of the node in the template"
            )

        if hasattr(self, "ec_x_ccoords"):
            ec_x_ccoords = self.ec_x_ccoords
            ecoords = self.ec_ccoords

        if not hasattr(self, "opt_method"):
            self.opt_method = "L-BFGS-B"

        if not hasattr(self, "constant_length"):
            self.constant_length = 1.54

        if not hasattr(self, "maxfun"):
            self.maxfun = 15000

        if not hasattr(self, "maxiter"):
            self.maxiter = 15000

        if not hasattr(self, "display"):
            self.display = True

        if not hasattr(self, "eps"):
            self.eps = 1e-8

        if not hasattr(self, "iprint"):
            self.iprint = -1

        G = self.G
        node_xcoords = self.node_x_ccoords
        node_coords = self.node_ccoords
        linker_length = self.linker_length
        opt_method = self.opt_method

        constant_length = self.constant_length

        x_com_length = np.mean([np.linalg.norm(i) for i in node_xcoords])
        sorted_nodes = sort_nodes_by_type_connectivity(G)

        # firstly, check if all V nodes have highest connectivity
        # secondly, sort all DV nodes by connectivity

        sorted_edges = find_and_sort_edges_bynodeconnectivity(G, sorted_nodes)

        nodes_atoms = {}
        for n in sorted_nodes:
            if "CV" in n:
                nodes_atoms[n] = self.ec_atom

            else:
                nodes_atoms[n] = self.node_atom

        Xatoms_positions_dict = {}
        node_positions_dict = {}
        # reindex the nodes in the Xatoms_positions with the index in the sorted_nodes, like G has 16 nodes[2,5,7], but the new dictionary should be [0,1,2]
        for n in sorted_nodes:
            if "CV" in n:
                Xatoms_positions_dict[sorted_nodes.index(n)] = addidx(
                    G.nodes[n]["ccoords"] + ec_x_ccoords
                )
            else:
                Xatoms_positions_dict[sorted_nodes.index(n)] = addidx(
                    G.nodes[n]["ccoords"] + node_xcoords
                )

        for n in sorted_nodes:
            if "CV" in n:
                node_positions_dict[sorted_nodes.index(n)] = (
                    G.nodes[n]["ccoords"] + ecoords
                )
            else:
                node_positions_dict[sorted_nodes.index(n)] = (
                    G.nodes[n]["ccoords"] + node_coords
                )

        # reindex the edges in the G with the index in the sorted_nodes
        sorted_edges_of_sortednodeidx = [
            (sorted_nodes.index(e[0]), sorted_nodes.index(e[1])) for e in sorted_edges
        ]

        # Optimize rotations
        num_nodes = G.number_of_nodes()
        pname_list = [pname(n) for n in sorted_nodes]
        pname_set = set(pname_list)
        pname_set_dict = {}
        for node_pname in pname_set:
            pname_set_dict[node_pname] = {
                "ind_ofsortednodes": [],
            }
        for i, node in enumerate(sorted_nodes):
            pname_set_dict[pname(node)]["ind_ofsortednodes"].append(i)
            if len(pname_set_dict[pname(node)]["ind_ofsortednodes"]) == 1:  # first node
                pname_set_dict[pname(node)]["rot_trans"] = get_rot_trans_matrix(
                    node, G, sorted_nodes, Xatoms_positions_dict
                )  # initial guess
        self.pname_set_dict = pname_set_dict

        for p_name in pname_set_dict:
            rot, trans = pname_set_dict[p_name]["rot_trans"]
            for k in pname_set_dict[p_name]["ind_ofsortednodes"]:
                node = sorted_nodes[k]
                Xatoms_positions_dict[k][:, 1:] = (
                    np.dot(
                        Xatoms_positions_dict[k][:, 1:] - G.nodes[node]["ccoords"], rot
                    )
                    + trans
                    + G.nodes[node]["ccoords"]
                )
                node_positions_dict[k] = (
                    np.dot(node_positions_dict[k] - G.nodes[node]["ccoords"], rot)
                    + trans
                    + G.nodes[node]["ccoords"]
                )
        ###3D free rotation
        if not hasattr(self, "saved_optimized_rotations"):
            print("-" * 80)
            print(" " * 20, "start to optimize the rotations", " " * 20)
            print("-" * 80)

            ##initial_rots = []
            ##
            ##for node in sorted_nodes:
            ##    rot = get_rotation_matrix(node, G, sorted_nodes, Xatoms_positions_dict)
            ##    # print(rot)
            ##    initial_rots.append(rot)
            ##initial_guess_rotations = np.array(
            ##    initial_rots
            ##).flatten()  # Initial guess for rotation matrices

            initial_guess_set_rotations = (
                np.eye(3, 3).reshape(1, 3, 3).repeat(len(pname_set), axis=0)
            )

            ##
            # (
            #    optimized_rotations_pre,
            #    _,
            # ) = optimize_rotations_pre(
            #    num_nodes,
            #    G,
            #    sorted_nodes,
            #    sorted_edges_of_sortednodeidx,
            #    Xatoms_positions_dict,
            #    initial_guess_set_rotations,
            #    pname_set_dict,
            #    opt_method=self.opt_method,
            #    maxfun=self.maxfun,
            #    maxiter=self.maxiter,
            #    disp=self.display,
            #    eps=self.eps,
            #    iprint=self.iprint,
            # )

            (
                optimized_set_rotations,
                _,
            ) = optimize_rotations_after(
                num_nodes,
                G,
                sorted_nodes,
                sorted_edges_of_sortednodeidx,
                Xatoms_positions_dict,
                initial_guess_set_rotations,
                # optimized_rotations_pre,
                pname_set_dict,
                opt_method=self.opt_method,
                maxfun=self.maxfun,
                maxiter=self.maxiter,
                disp=self.display,
                eps=self.eps,
                iprint=self.iprint,
            )
            print("-" * 80)
            print(" " * 20, "rotations optimization completed", " " * 20)
            print("-" * 80)
            # to save the optimized rotations as npy
            if hasattr(self, "to_save_optimized_rotations_filename"):
                np.save(
                    self.to_save_optimized_rotations_filename + ".npy",
                    optimized_set_rotations,
                )
                print(
                    "optimized rotations are saved to: ",
                    self.to_save_optimized_rotations_filename + ".npy",
                )

        else:
            if hasattr(self, "use_saved_rotations_as_initial_guess"):
                if self.use_saved_rotations_as_initial_guess:
                    print("use the saved optimized_rotations as initial guess")
                    print("-" * 80)
                    print(" " * 20, "start to optimize the rotations", " " * 20)
                    print("-" * 80)

                    saved_set_rotations = self.saved_optimized_rotations.reshape(
                        -1, 3, 3
                    )
                    # (
                    #    optimized_rotations_pre,
                    #    _,
                    # ) = optimize_rotations_pre(
                    #    num_nodes,
                    #    G,
                    #    sorted_nodes,
                    #    sorted_edges_of_sortednodeidx,
                    #    Xatoms_positions_dict,
                    #    saved_set_rotations,
                    #    pname_set_dict,
                    #    opt_method=self.opt_method,
                    #    maxfun=self.maxfun,
                    #    maxiter=self.maxiter,
                    #    disp=self.display,
                    #    eps=self.eps,
                    #    iprint=self.iprint,
                    # )

                    (
                        optimized_set_rotations,
                        _,
                    ) = optimize_rotations_after(
                        num_nodes,
                        G,
                        sorted_nodes,
                        sorted_edges_of_sortednodeidx,
                        Xatoms_positions_dict,
                        saved_set_rotations,
                        pname_set_dict,
                        opt_method=self.opt_method,
                        maxfun=self.maxfun,
                        maxiter=self.maxiter,
                        disp=self.display,
                        eps=self.eps,
                        iprint=self.iprint,
                    )
                    print("-" * 80)
                    print(" " * 20, "rotations optimization completed", " " * 20)
                    print("-" * 80)
                    # to save the optimized rotations as npy
                    if hasattr(self, "to_save_optimized_rotations_filename"):
                        np.save(
                            self.to_save_optimized_rotations_filename + ".npy",
                            optimized_set_rotations,
                        )
                        print(
                            "optimized rotations are saved to: ",
                            self.to_save_optimized_rotations_filename + ".npy",
                        )

                else:
                    optimized_set_rotations = self.saved_optimized_rotations.reshape(
                        -1, 3, 3
                    )

            else:
                print(
                    "use the loaded optimized_rotations from the previous optimization"
                )
                optimized_set_rotations = self.saved_optimized_rotations.reshape(
                    -1, 3, 3
                )

        optimized_rotations = expand_setrots(
            pname_set_dict, optimized_set_rotations, sorted_nodes
        )
        # Apply rotations
        rotated_node_positions = apply_rotations_to_atom_positions(
            optimized_rotations, G, sorted_nodes, node_positions_dict
        )

        # Save results to XYZ
        # save_xyz("optimized_nodesstructure.xyz", rotated_node_positions) #DEBUG

        rotated_Xatoms_positions_dict, optimized_pair = (
            apply_rotations_to_Xatoms_positions(
                optimized_rotations,
                G,
                sorted_nodes,
                sorted_edges_of_sortednodeidx,
                Xatoms_positions_dict,
            )
        )

        start_node = sorted_edges[0][0]  # find_nearest_node_to_beginning_point(G)
        # loop all of the edges in G and get the lengths of the edges, length is the distance between the two nodes ccoords
        edge_lengths, lengths = get_edge_lengths(G)

        x_com_length = np.mean([np.linalg.norm(i) for i in node_xcoords])
        new_edge_length = linker_length + 2 * constant_length + 2 * x_com_length
        # update the node ccoords in G by loop edge, start from the start_node, and then update the connected node ccoords by the edge length, and update the next node ccords from the updated node

        updated_ccoords, original_ccoords = update_node_ccoords(
            G, edge_lengths, start_node, new_edge_length
        )
        # exclude the start_node in updated_ccoords and original_ccoords
        updated_ccoords = {k: v for k, v in updated_ccoords.items() if k != start_node}
        original_ccoords = {
            k: v for k, v in original_ccoords.items() if k != start_node
        }

        # use optimized_params to update all of nodes ccoords in G, according to the fccoords
        if not hasattr(self, "optimized_params"):
            print("-" * 80)
            print(" " * 20, "start to optimize the cell parameters", " " * 20)
            print("-" * 80)
            optimized_params = optimize_cell_parameters(
                self.cell_info, original_ccoords, updated_ccoords
            )
            print("-" * 80)
            print(" " * 20, "cell parameters optimization completed", " " * 20)
            print("-" * 80)
        else:
            print("use the optimized_params from the previous optimization")
            optimized_params = self.optimized_params

        sc_unit_cell = unit_cell_to_cartesian_matrix(
            optimized_params[0],
            optimized_params[1],
            optimized_params[2],
            optimized_params[3],
            optimized_params[4],
            optimized_params[5],
        )
        sc_unit_cell_inv = np.linalg.inv(sc_unit_cell)
        sG, scaled_ccoords = update_ccoords_by_optimized_cell_params(
            self.G, optimized_params
        )
        scaled_node_positions_dict = {}
        scaled_Xatoms_positions_dict = {}

        for n in sorted_nodes:
            if "CV" in n:
                scaled_Xatoms_positions_dict[sorted_nodes.index(n)] = addidx(
                    sG.nodes[n]["ccoords"] + ec_x_ccoords
                )
            else:
                scaled_Xatoms_positions_dict[sorted_nodes.index(n)] = addidx(
                    sG.nodes[n]["ccoords"] + node_xcoords
                )

        for n in sorted_nodes:
            if "CV" in n:
                scaled_node_positions_dict[sorted_nodes.index(n)] = (
                    sG.nodes[n]["ccoords"] + ecoords
                )
            else:
                scaled_node_positions_dict[sorted_nodes.index(n)] = (
                    sG.nodes[n]["ccoords"] + node_coords
                )

        # Apply rotations
        for p_name in pname_set_dict:
            rot, trans = pname_set_dict[p_name]["rot_trans"]
            for k in pname_set_dict[p_name]["ind_ofsortednodes"]:
                node = sorted_nodes[k]
                scaled_Xatoms_positions_dict[k][:, 1:] = (
                    np.dot(
                        scaled_Xatoms_positions_dict[k][:, 1:]
                        - sG.nodes[node]["ccoords"],
                        rot,
                    )
                    + trans
                    + sG.nodes[node]["ccoords"]
                )

                scaled_node_positions_dict[k] = (
                    np.dot(
                        scaled_node_positions_dict[k] - sG.nodes[node]["ccoords"], rot
                    )
                    + trans
                    + sG.nodes[node]["ccoords"]
                )

        scaled_rotated_node_positions = apply_rotations_to_atom_positions(
            optimized_rotations, sG, sorted_nodes, scaled_node_positions_dict
        )
        scaled_rotated_Xatoms_positions, optimized_pair = (
            apply_rotations_to_Xatoms_positions(
                optimized_rotations,
                sG,
                sorted_nodes,
                sorted_edges_of_sortednodeidx,
                scaled_Xatoms_positions_dict,
            )
        )
        # Save results to XYZ

        self.sorted_nodes = sorted_nodes
        self.sorted_edges = sorted_edges
        self.sorted_edges_of_sortednodeidx = sorted_edges_of_sortednodeidx
        self.optimized_rotations = optimized_rotations
        self.optimized_params = optimized_params
        self.new_edge_length = new_edge_length
        self.optimized_pair = optimized_pair
        self.scaled_rotated_node_positions = scaled_rotated_node_positions
        self.scaled_rotated_Xatoms_positions = scaled_rotated_Xatoms_positions
        self.sc_unit_cell = sc_unit_cell
        self.sc_unit_cell_inv = sc_unit_cell_inv
        self.sG_node = sG
        self.nodes_atom = nodes_atoms
        self.rotated_node_positions = rotated_node_positions
        self.Xatoms_positions_dict = Xatoms_positions_dict
        self.node_positions_dict = node_positions_dict
        # save_xyz("scale_optimized_nodesstructure.xyz", scaled_rotated_node_positions)

    def place_edge_in_net(self):
        """
        based on the optimized rotations and cell parameters, use optimized pair to find connected X-X pair in optimized cell,
        and place the edge in the target MOF cell

        return:
            sG (networkx graph):graph of the target MOF cell, with scaled and rotated node and edge positions
        """
        # linker_middle_point = np.mean(linker_x_vecs,axis=0)
        linker_xx_vec = self.linker_x_ccoords
        linker_length = self.linker_length
        optimized_pair = self.optimized_pair
        scaled_rotated_Xatoms_positions = self.scaled_rotated_Xatoms_positions
        scaled_rotated_node_positions = self.scaled_rotated_node_positions
        sorted_nodes = self.sorted_nodes
        sG_node = self.sG_node
        sc_unit_cell_inv = self.sc_unit_cell_inv
        nodes_atom = self.nodes_atom

        sG = sG_node.copy()
        scalar = (linker_length + 2 * self.constant_length) / linker_length
        extended_linker_xx_vec = [i * scalar for i in linker_xx_vec]
        norm_xx_vector_record = []
        rot_record = []

        # edges = {}
        for (i, j), pair in optimized_pair.items():
            x_idx_i, x_idx_j = pair
            reindex_i = sorted_nodes.index(i)
            reindex_j = sorted_nodes.index(j)
            x_i = scaled_rotated_Xatoms_positions[reindex_i][x_idx_i][1:]
            x_j = scaled_rotated_Xatoms_positions[reindex_j][x_idx_j][1:]
            x_i_x_j_middle_point = np.mean([x_i, x_j], axis=0)
            xx_vector = np.vstack(
                [x_i - x_i_x_j_middle_point, x_j - x_i_x_j_middle_point]
            )
            norm_xx_vector = xx_vector / np.linalg.norm(xx_vector)

            # print(i, j, reindex_i, reindex_j, x_idx_i, x_idx_j)
            # use superimpose to get the rotation matrix
            # use record to record the rotation matrix for get rid of the repeat calculation
            indices = [
                index
                for index, value in enumerate(norm_xx_vector_record)
                if is_list_A_in_B(norm_xx_vector, value)
            ]
            if len(indices) == 1:
                rot = rot_record[indices[0]]
                # rot = reorthogonalize_matrix(rot)
            else:
                _, rot, _ = superimpose_rotation_only(extended_linker_xx_vec, xx_vector)
                # rot = reorthogonalize_matrix(rot)
                norm_xx_vector_record.append(norm_xx_vector)
                # the rot may be opposite, so we need to check the angle between the two vectors
                # if the angle is larger than 90 degree, we need to reverse the rot
                roted_xx = np.dot(extended_linker_xx_vec, rot)

                if np.dot(roted_xx[1] - roted_xx[0], xx_vector[1] - xx_vector[0]) < 0:
                    ##rotate 180 around the axis of the cross product of the two vectors
                    axis = np.cross(
                        roted_xx[1] - roted_xx[0], xx_vector[1] - xx_vector[0]
                    )
                    # if 001 not linear to the two vectors
                    if np.linalg.norm(axis) == 0:
                        check_z_axis = np.cross(roted_xx[1] - roted_xx[0], [0, 0, 1])
                        if np.linalg.norm(check_z_axis) == 0:
                            axis = np.array([1, 0, 0])
                        else:
                            axis = np.array([0, 0, 1])

                    axis = axis / np.linalg.norm(axis)
                    flip_matrix = np.eye(3) - 2 * np.outer(
                        axis, axis
                    )  # Householder matrix for reflection
                    rot = np.dot(rot, flip_matrix)
                # Flip the last column of the rotation matrix if the determinant is negative

                rot_record.append(rot)

            # use the rotation matrix to rotate the linker x coords
            # rotated_xx = np.dot(extended_linker_xx_vec, rot)
            # print(rotated_xx,'rotated_xx',xx_vector) #DEBUG
            placed_edge_ccoords = (
                np.dot(self.linker_ccoords, rot) + x_i_x_j_middle_point
            )

            placed_edge = np.hstack((np.asarray(self.linker_atom), placed_edge_ccoords))
            sG.edges[(i, j)]["coords"] = x_i_x_j_middle_point
            sG.edges[(i, j)]["c_points"] = placed_edge

            sG.edges[(i, j)]["f_points"] = np.hstack(
                (
                    placed_edge[:, 0:2],
                    cartesian_to_fractional(placed_edge[:, 2:5], sc_unit_cell_inv),
                )
            )  # NOTE: modified add the atom type and atom name

            _, sG.edges[(i, j)]["x_coords"] = fetch_X_atoms_ind_array(
                placed_edge, 0, "X"
            )
            # edges[(i,j)]=placed_edge
        # placed_node = {}
        for k, v in scaled_rotated_node_positions.items():
            # print(k,v)
            # placed_node[k] = np.hstack((nodes_atom[k],v))
            sG.nodes[k]["c_points"] = np.hstack((nodes_atom[k], v))
            sG.nodes[k]["f_points"] = np.hstack(
                (nodes_atom[k], cartesian_to_fractional(v, sc_unit_cell_inv))
            )
            # find the atoms starts with "x" and extract the coordinates
            _, sG.nodes[k]["x_coords"] = fetch_X_atoms_ind_array(
                sG.nodes[k]["c_points"], 0, "X"
            )
        self.sG = sG
        return sG

    def set_supercell(self, supercell):
        """
        set the supercell of the target MOF model
        """
        self.supercell = supercell

    def make_supercell_multitopic(self):
        """
        make the supercell of the multitopic linker MOF
        """
        sG = self.sG
        self.multiedge_bundlings = bundle_multiedge(sG)
        # self.dv_v_pairs, sG = replace_DV_with_corresponding_V(sG) #debug
        superG = update_supercell_node_fpoints_loose(sG, self.supercell)
        superG = update_supercell_edge_fpoints(sG, superG, self.supercell)
        # self.prim_multiedge_bundlings = replace_bundle_dvnode_with_vnode(  #debug
        #    self.dv_v_pairs, self.multiedge_bundlings
        # )
        self.prim_multiedge_bundlings = self.multiedge_bundlings
        self.super_multiedge_bundlings = make_super_multiedge_bundlings(
            self.prim_multiedge_bundlings, self.supercell
        )
        superG = update_supercell_bundle(superG, self.super_multiedge_bundlings)
        superG = check_multiedge_bundlings_insuperG(
            self.super_multiedge_bundlings, superG
        )
        self.superG = superG
        return superG

    def make_supercell_ditopic(self):
        """
        make the supercell of the ditopic linker MOF
        """

        sG = self.sG
        # self.dv_v_pairs, sG = replace_DV_with_corresponding_V(sG)
        superG = update_supercell_node_fpoints_loose(sG, self.supercell)
        superG = update_supercell_edge_fpoints(sG, superG, self.supercell)
        self.superG = superG
        return superG

    def set_virtual_edge(self, bool_x=False, range=0.5, max_neighbor=2):
        """
        set the virtual edge addition for the bridge type nodes,
        range is the range to search the virtual edge between two Vnodes directly, should <= 0.5,
        max_neighbor is the maximum number of neighbors of the node with virtual edge
        """

        self.add_virtual_edge = bool(bool_x)
        self.vir_edge_range = range
        self.vir_edge_max_neighbor = max_neighbor

    def add_virtual_edge_for_bridge_node(self):
        """
        after setting the virtual edge search, add the virtual edge to the target supercell superG MOF
        """
        if self.add_virtual_edge:
            superG = add_virtual_edge(
                self.superG, self.vir_edge_range, self.vir_edge_max_neighbor
            )
            print("add virtual edge")
            self.superG = superG
        else:
            pass

    def set_remove_node_list(self, remove_node_list):
        """
        make defect in the target MOF model by removing nodes
        """
        self.remove_node_list = remove_node_list

    def set_remove_edge_list(self, remove_edge_list):
        """
        make defect in the target MOF model by removing edges
        """
        self.remove_edge_list = remove_edge_list

    def make_eG_from_supereG_multitopic(self):
        """
        make the target MOF cell graph with only EDGE and V, link the XOO atoms to the EDGE
        always need to execute with make_supercell_multitopic
        """

        eG, _ = superG_to_eG_multitopic(self.superG, self.sc_unit_cell)
        self.eG = eG
        return eG

    def add_xoo_to_edge_multitopic(self):
        eG = self.eG
        eG, unsaturated_linker, matched_vnode_xind, xoo_dict = addxoo2edge_multitopic(
            eG, self.sc_unit_cell
        )
        self.unsaturated_linker = unsaturated_linker
        self.matched_vnode_xind = matched_vnode_xind
        self.xoo_dict = xoo_dict
        self.eG = eG
        return eG

    def make_eG_from_supereG_ditopic(self):
        """
        make the target MOF cell graph with only EDGE and V, link the XOO atoms to the EDGE
        always execute with make_supercell_ditopic
        """

        eG, _ = superG_to_eG_ditopic(self.superG)
        self.eG = eG
        return eG

    def add_xoo_to_edge_ditopic(self):
        """
        analyze eG and link the XOO atoms to the EDGE, update eG, for ditopic linker MOF
        """
        eG = self.eG
        eG, unsaturated_linker, matched_vnode_xind, xoo_dict = addxoo2edge_ditopic(
            eG, self.sc_unit_cell
        )
        self.unsaturated_linker = unsaturated_linker
        self.matched_vnode_xind = matched_vnode_xind
        self.xoo_dict = xoo_dict
        self.eG = eG
        return eG

    def main_frag_eG(self):
        """
        only keep the main fragment of the target MOF cell, remove the other fragments, to avoid the disconnected fragments
        """
        eG = self.eG
        self.eG = [eG.subgraph(c).copy() for c in nx.connected_components(eG)][0]
        print("main fragment of the MOF cell is kept")  # ,len(self.eG.nodes()),'nodes')
        # print('fragment size list:',[len(c) for c in nx.connected_components(eG)]) #debug
        return self.eG

    def make_supercell_range_cleaved_eG(self, buffer_plus=0, buffer_minus=0):
        supercell = self.supercell
        new_eG = self.eG.copy()
        eG = self.eG
        removed_edges = []
        removed_nodes = []
        for n in eG.nodes():
            if pname(n) != "EDGE":
                if check_supercell_box_range(
                    eG.nodes[n]["fcoords"], supercell, buffer_plus, buffer_minus
                ):
                    pass
                else:
                    new_eG.remove_node(n)
                    removed_nodes.append(n)
            elif pname(n) == "EDGE":
                if (
                    arr_dimension(eG.nodes[n]["fcoords"]) == 2
                ):  # ditopic linker have two points in the fcoords
                    edge_coords = np.mean(eG.nodes[n]["fcoords"], axis=0)
                elif (
                    arr_dimension(eG.nodes[n]["fcoords"]) == 1
                ):  # multitopic linker have one point in the fcoords from EC
                    edge_coords = eG.nodes[n]["fcoords"]

                if check_supercell_box_range(
                    edge_coords, supercell, buffer_plus, buffer_minus
                ):
                    pass
                else:
                    new_eG.remove_node(n)
                    removed_edges.append(n)

        matched_vnode_xind = self.matched_vnode_xind
        self.matched_vnode_xind = update_matched_nodes_xind(
            removed_nodes,
            removed_edges,
            matched_vnode_xind,
        )

        self.eG = new_eG
        return new_eG, removed_edges, removed_nodes

    def set_node_topic(self, node_topic):
        """
        manually set the node topic, normally should be the same as the maximum degree of the node in the template
        """
        self.node_topic = node_topic

    def find_unsaturated_node_eG(self):
        """
        use the eG to find the unsaturated nodes, whose degree is less than the node topic
        """
        eG = self.eG
        if hasattr(self, "node_topic"):
            node_topic = self.node_topic
        else:
            node_topic = self.node_max_degree
        unsaturated_node = find_unsaturated_node(eG, node_topic)
        self.unsaturated_node = unsaturated_node
        return unsaturated_node

    def find_unsaturated_linker_eG(eG, linker_topics):
        """
        use the eG to find the unsaturated linkers, whose degree is less than linker topic
        """
        new_unsaturated_linker = find_unsaturated_linker(eG, linker_topics)
        return new_unsaturated_linker

    def set_node_terminamtion(self, term_file):
        """
        pdb file, set the node termination file, which contains the information of the node terminations, should have X of connected atom (normally C),
        Y of two connected O atoms (if in carboxylate group) to assist the placement of the node terminations
        """

        term_data = termpdb(term_file)
        term_info = term_data[:, :-3]
        term_coords = term_data[:, -3:]
        xterm, _ = Xpdb(term_data, "X")
        oterm, _ = Xpdb(term_data, "Y")
        term_xvecs = xterm[:, -3:]
        term_ovecs = oterm[:, -3:]
        term_coords = term_coords.astype("float")
        term_xvecs = term_xvecs.astype("float")
        term_ovecs = term_ovecs.astype("float")

        term_ovecs_c = np.mean(np.asarray(term_ovecs), axis=0)
        term_coords = term_coords - term_ovecs_c
        term_xoovecs = np.vstack((term_xvecs, term_ovecs))
        term_xoovecs = term_xoovecs - term_ovecs_c

        self.term_info = term_info
        self.term_coords = term_coords
        self.term_xoovecs = term_xoovecs

    # Function to add node terminations
    def add_terminations_to_unsaturated_node(self):
        """
        use the node terminations to add terminations to the unsaturated nodes

        """
        unsaturated_node = [n for n in self.unsaturated_node if n in self.eG.nodes()]
        xoo_dict = self.xoo_dict
        matched_vnode_xind = self.matched_vnode_xind
        eG = self.eG
        sc_unit_cell = self.sc_unit_cell
        (
            unsaturated_vnode_xind_dict,
            unsaturated_vnode_xoo_dict,
            self.matched_vnode_xind_dict,
        ) = make_unsaturated_vnode_xoo_dict(
            unsaturated_node, xoo_dict, matched_vnode_xind, eG, sc_unit_cell
        )
        # term_file: path to the termination file
        # ex_node_cxo_cc: exposed node coordinates

        node_oovecs_record = []
        for n in eG.nodes():
            eG.nodes[n]["term_c_points"] = {}
        for exvnode_xind_key in unsaturated_vnode_xoo_dict.keys():
            exvnode_x_ccoords = unsaturated_vnode_xoo_dict[exvnode_xind_key][
                "x_cpoints"
            ]
            exvnode_oo_ccoords = unsaturated_vnode_xoo_dict[exvnode_xind_key][
                "oo_cpoints"
            ]
            node_xoo_ccoords = np.vstack([exvnode_x_ccoords, exvnode_oo_ccoords])
            # make the beginning point of the termination to the center of the oo atoms
            node_oo_center_cvec = np.mean(
                exvnode_oo_ccoords[:, 2:5].astype(float), axis=0
            )  # NOTE: modified add the atom type and atom name
            node_xoo_cvecs = (
                node_xoo_ccoords[:, 2:5].astype(float) - node_oo_center_cvec
            )  # NOTE: modified add the atom type and atom name
            node_xoo_cvecs = node_xoo_cvecs.astype("float")
            # use record to record the rotation matrix for get rid of the repeat calculation

            indices = [
                index
                for index, value in enumerate(node_oovecs_record)
                if is_list_A_in_B(node_xoo_cvecs, value[0])
            ]
            if len(indices) == 1:
                rot = node_oovecs_record[indices[0]][1]
            else:
                _, rot, _ = superimpose(self.term_xoovecs, node_xoo_cvecs)

                node_oovecs_record.append((node_xoo_cvecs, rot))
            adjusted_term_vecs = np.dot(self.term_coords, rot) + node_oo_center_cvec
            adjusted_term = np.hstack(
                (
                    np.asarray(self.term_info[:, 0:1]),
                    np.asarray(self.term_info[:, 2:3]),
                    adjusted_term_vecs,
                )
            )
            # add the adjusted term to the terms, add index, add the node name
            unsaturated_vnode_xoo_dict[exvnode_xind_key]["node_term_c_points"] = (
                adjusted_term
            )
            eG.nodes[exvnode_xind_key[0]]["term_c_points"][exvnode_xind_key[1]] = (
                adjusted_term
            )

        self.unsaturated_vnode_xoo_dict = unsaturated_vnode_xoo_dict
        self.eG = eG
        return eG

    def remove_xoo_from_node(self):
        """
        remove the XOO atoms from the node after adding the terminations, add ['noxoo_f_points'] to the node in eG
        """
        eG = self.eG
        xoo_dict = self.xoo_dict

        all_xoo_indices = []
        for x_ind, oo_ind in xoo_dict.items():
            all_xoo_indices.append(x_ind)
            all_xoo_indices.extend(oo_ind)

        for n in eG.nodes():
            if pname(n) != "EDGE":
                all_f_points = eG.nodes[n]["f_points"]
                noxoo_f_points = np.delete(all_f_points, all_xoo_indices, axis=0)
                eG.nodes[n]["noxoo_f_points"] = noxoo_f_points
        self.eG = eG

        return eG

    def write_node_edge_node_gro(self, gro_name):
        """
        write the node, edge, node to the gro file
        """

        nodes_eG, edges_eG, terms_eG, node_res_num, edge_res_num, term_res_num = (
            extract_node_edge_term(self.eG, self.sc_unit_cell)
        )
        merged_node_edge_term = merge_node_edge_term(
            nodes_eG, edges_eG, terms_eG, node_res_num, edge_res_num
        )
        dir_name = "output_gros"
        save_node_edge_term_gro(merged_node_edge_term, gro_name, dir_name)
        print(str(gro_name) + ".gro is saved in folder " + str(dir_name))
        print("node_res_num: ", node_res_num)
        print("edge_res_num: ", edge_res_num)
        print("term_res_num: ", term_res_num)

        self.nodes_eG = nodes_eG
        self.edges_eG = edges_eG
        self.terms_eG = terms_eG
        self.node_res_num = node_res_num
        self.edge_res_num = edge_res_num
        self.term_res_num = term_res_num
        self.merged_node_edge_term = merged_node_edge_term
