from pathlib import Path
import numpy as np
import argparse
import hashlib
import sys
import re


def get_vlx_basis_string(basis_name,
                         output_name,
                         use_gc_and_sp=True,
                         optimize_general=True):
    """
    Gets string for VeloxChem basis set.

    :param basis_name:
        The name of basis set.
    :param output_name:
        The name of output file.
    :param use_gc_and_sp:
        The option for using format for general contraction and SP shells.
    :param optimize_general:
        The option for using optimization for general contraction.

    :return:
        The string for VeloxChem basis set.
    """

    # geometric augmentation

    geom_aug = 0
    if basis_name.upper().startswith('DAUG-'):
        geom_aug = 1
    elif basis_name.upper().startswith('TAUG-'):
        geom_aug = 2
    elif basis_name.upper().startswith('QAUG-'):
        geom_aug = 3

    # (aug-)cc-pCVnZ basis sets: H and He are added from (aug-)cc-pVnZ

    is_cc_pcvnz_basis = (re.search(r'^[DTQ]?(AUG-)?CC-PCV.Z$',
                                   basis_name.upper()) is not None)

    # Process basis set name (or json file name)

    if basis_name[-5:] == '.json' and Path(basis_name).is_file():
        import json
        with open(basis_name, 'r') as f_json:
            basis = json.load(f_json)
        basis_name = Path(basis_name).stem
    else:
        if basis_name.upper() == 'SADLEJ-PVTZ':
            basis_name = 'SADLEJ PVTZ'
        import basis_set_exchange as bse
        if geom_aug >= 1:
            basis_name = basis_name[1:]
        basis = bse.get_basis(basis_name,
                              header=False,
                              optimize_general=optimize_general)
        if is_cc_pcvnz_basis:
            h_he_basis_name = basis_name.upper().replace('CC-PCV', 'CC-PV')
            h_he_basis = bse.get_basis(h_he_basis_name,
                                       header=False,
                                       optimize_general=optimize_general)
        if geom_aug >= 1:
            basis = bse.manip.geometric_augmentation(basis, geom_aug)
            basis_name = 'DTQ'[geom_aug - 1] + basis_name
            if is_cc_pcvnz_basis:
                h_he_basis = bse.manip.geometric_augmentation(
                    h_he_basis, geom_aug)
                h_he_basis_name = 'DTQ'[geom_aug - 1] + h_he_basis_name
        if basis_name.upper() == 'SADLEJ PVTZ':
            basis_name = 'SADLEJ-PVTZ'

    basis_title = basis_name
    if output_name is not None:
        basis_title = f'{Path(output_name).name}'
        if basis_title.upper() != basis_name.upper():
            basis_title += f'  !{basis_name.upper()}'
            if is_cc_pcvnz_basis:
                basis_title += f'  (H and He from {h_he_basis_name})'
    vlx_basis_str = f'@BASIS_SET {basis_title}\n'

    # Go through elements (up to Rn)

    elem_labels = [
        '0', 'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg',
        'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr',
        'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr',
        'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd',
        'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd',
        'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf',
        'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po',
        'At', 'Rn'
    ]

    basis_elements = {}
    if is_cc_pcvnz_basis:
        for elem, atombasis in h_he_basis['elements'].items():
            if int(elem) in [1, 2]:
                basis_elements[elem] = h_he_basis['elements'][elem]
    basis_elements.update(basis['elements'])

    for elem, atombasis in basis_elements.items():

        if int(elem) > 86:
            continue

        vlx_basis_str += f'\n@ATOMBASIS {elem_labels[int(elem)].upper()}\n'

        target_ang_mom = ['S', 'SP', 'P', 'D', 'F', 'G', 'H', 'I']
        atom_shells = {ang_mom: [] for ang_mom in target_ang_mom}

        # Go through shells

        for shell in atombasis['electron_shells']:

            # Check angular momentum

            if len(shell['angular_momentum']) == 1:
                ang_mom = 'SPDFGHI'[shell['angular_momentum'][0]]
            elif shell['angular_momentum'] == [0, 1]:
                ang_mom = 'SP'
            elif shell['angular_momentum'] == [0, 1, 2]:
                ang_mom = 'SPD'
            else:
                ang_mom = None

            # Get exponents and coefficients

            expon_vec = np.array(shell['exponents'], dtype=np.float64)
            coeff_mat = np.array(shell['coefficients'], dtype=np.float64)

            # Mark non-zero coefficients as True

            bool_mat = (coeff_mat != 0.0)

            # An example bool_mat:
            #   [[ True  True  True  True  True  True False False]
            #   [[ True  True  True  True  True  True False False]
            #    [False False False  True False  True False False]
            #    [False False False False Fasle False  True False]]
            #    [False False False False Fasle False False  True]]

            # Find out identical rows and save row indices in groups

            row_indices_list = [[row] for row in range(bool_mat.shape[0])]
            while True:
                merged_row_groups = False
                pair_inds = [(a, b)
                             for a in range(len(row_indices_list))
                             for b in range(a + 1, len(row_indices_list))]
                for a, b in pair_inds:
                    row_a = bool_mat[row_indices_list[a][0], :]
                    row_b = bool_mat[row_indices_list[b][0], :]
                    if (row_a == row_b).all():
                        row_indices_list[a] += row_indices_list[b]
                        row_indices_list[b] = None
                        merged_row_groups = True
                        break
                if merged_row_groups:
                    row_indices_list = [
                        row_indices for row_indices in row_indices_list
                        if row_indices is not None
                    ]
                else:
                    break

            # Find out indices of non-zero columns for each row group

            col_indices_list = []
            for row_indices in row_indices_list:
                col_indices = [
                    col for col, val in enumerate(bool_mat[row_indices[0]])
                    if val
                ]
                col_indices_list.append(col_indices)

            # Go through the blocks and write basis set

            for (row_indices, col_indices) in zip(row_indices_list,
                                                  col_indices_list):

                # Form the exponent/coefficient matrix (column-wise)

                expons = expon_vec[col_indices]
                coeffs = coeff_mat[row_indices, :][:, col_indices]
                matrix = np.vstack((expons, coeffs)).T

                # Get number of primitives

                nprims = matrix.shape[0]

                # Print exponents and coefficients

                if len(ang_mom) == 1:
                    if use_gc_and_sp:
                        ngc = matrix.shape[1] - 1
                        shell_str = f'{ang_mom} {nprims:d}   {ngc:d}\n'
                        for i in range(matrix.shape[0]):
                            shell_str += f'{matrix[i, 0]:18.12e}'
                            for j in range(1, matrix.shape[1]):
                                shell_str += f'{matrix[i, j]:20.12e}'
                            shell_str += '\n'
                        atom_shells[ang_mom].append(shell_str)
                    else:
                        for j in range(1, matrix.shape[1]):
                            shell_str = f'{ang_mom} {nprims:d}   1\n'
                            for i in range(matrix.shape[0]):
                                shell_str += f'{matrix[i, 0]:18.12e}'
                                shell_str += f'{matrix[i, j]:20.12e}'
                                shell_str += '\n'
                            atom_shells[ang_mom].append(shell_str)

                elif ang_mom == 'SP':
                    assert matrix.shape[1] == 3
                    if use_gc_and_sp:
                        shell_str = f'SP {nprims:d}   1\n'
                        for i in range(matrix.shape[0]):
                            shell_str += f'{matrix[i, 0]:18.12e}'
                            shell_str += f'{matrix[i, 1]:20.12e}'
                            shell_str += f'{matrix[i, 2]:20.12e}'
                            shell_str += '\n'
                        atom_shells['SP'].append(shell_str)
                    else:
                        shell_str = f'S {nprims:d}   1\n'
                        for i in range(matrix.shape[0]):
                            shell_str += f'{matrix[i, 0]:18.12e}'
                            shell_str += f'{matrix[i, 1]:20.12e}'
                            shell_str += '\n'
                        atom_shells['SP'].append(shell_str)
                        shell_str = f'P {nprims:d}   1\n'
                        for i in range(matrix.shape[0]):
                            shell_str += f'{matrix[i, 0]:18.12e}'
                            shell_str += f'{matrix[i, 2]:20.12e}'
                            shell_str += '\n'
                        atom_shells['SP'].append(shell_str)

                elif ang_mom == 'SPD':
                    assert matrix.shape[1] == 4
                    if use_gc_and_sp:
                        shell_str = f'SP {nprims:d}   1\n'
                        for i in range(matrix.shape[0]):
                            shell_str += f'{matrix[i, 0]:18.12e}'
                            shell_str += f'{matrix[i, 1]:20.12e}'
                            shell_str += f'{matrix[i, 2]:20.12e}'
                            shell_str += '\n'
                        atom_shells['SP'].append(shell_str)
                    else:
                        shell_str = f'S {nprims:d}   1\n'
                        for i in range(matrix.shape[0]):
                            shell_str += f'{matrix[i, 0]:18.12e}'
                            shell_str += f'{matrix[i, 1]:20.12e}'
                            shell_str += '\n'
                        atom_shells['SP'].append(shell_str)
                        shell_str = f'P {nprims:d}   1\n'
                        for i in range(matrix.shape[0]):
                            shell_str += f'{matrix[i, 0]:18.12e}'
                            shell_str += f'{matrix[i, 2]:20.12e}'
                            shell_str += '\n'
                        atom_shells['SP'].append(shell_str)
                    shell_str = f'D {nprims:d}   1\n'
                    for i in range(matrix.shape[0]):
                        shell_str += f'{matrix[i, 0]:18.12e}'
                        shell_str += f'{matrix[i, 3]:20.12e}'
                        shell_str += '\n'
                    atom_shells['D'].append(shell_str)

        for ang_mom in target_ang_mom:
            vlx_basis_str += ''.join(atom_shells[ang_mom])

        vlx_basis_str += '@END\n'

    md5sum = hashlib.md5(vlx_basis_str.encode('utf-8')).hexdigest()

    return vlx_basis_str + md5sum


if __name__ == '__main__':

    # Process command line arguments

    usage = f'python {Path(__file__).name} '
    usage += '--basis <name_of_basis_set> '
    usage += '[--output <name_of_output_file>]'

    parser = argparse.ArgumentParser(usage=usage)

    parser.add_argument('--basis',
                        dest='basis_name',
                        type=str,
                        default=None,
                        help='name of basis set')

    parser.add_argument('--output',
                        dest='output_name',
                        type=str,
                        default=None,
                        help='name of output file')

    args = parser.parse_args()

    if args.basis_name is None:
        parser.print_help()
        sys.exit(0)

    basis_name = args.basis_name
    output_name = args.output_name

    if output_name is not None:
        output_file = Path(output_name)
        output_name = str(output_file.parent / output_file.name.upper())

    vlx_basis_str = get_vlx_basis_string(basis_name,
                                         output_name,
                                         use_gc_and_sp=False,
                                         optimize_general=True)

    if output_name is None:
        print(vlx_basis_str)
    else:
        with Path(output_name).open('w') as f_out:
            f_out.write(vlx_basis_str)
