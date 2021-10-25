from pathlib import Path
import numpy as np
import argparse
import hashlib
import sys


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
        basis = bse.get_basis(basis_name,
                              header=False,
                              optimize_general=optimize_general)
        if basis_name.upper() == 'SADLEJ PVTZ':
            basis_name = 'SADLEJ-PVTZ'

    basis_title = basis_name
    if output_name is not None:
        basis_title = f'{Path(output_name).name}  !{basis_name.upper()}'
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

    for elem, atombasis in basis['elements'].items():

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
            #   [[ True  True  True  True  True  True False]
            #    [False False False  True False  True False]
            #    [False False False False Fasle False  True]]
            # Need to determine the row- and column-indices of the True values.

            row_indices_list = []
            row_start = 0
            for row in range(1, bool_mat.shape[0]):
                if not (bool_mat[row, :] == bool_mat[row - 1, :]).all():
                    row_indices_list.append(list(range(row_start, row)))
                    row_start = row
            if row_start != bool_mat.shape[0]:
                row_indices_list.append(
                    list(range(row_start, bool_mat.shape[0])))

            col_indices_list = []
            for row_indices in row_indices_list:
                col_indices = []
                for col, val in enumerate(bool_mat[row_indices[0]]):
                    if val:
                        col_indices.append(col)
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
                        shell_str = f'{ang_mom} {nprims:d}  {ngc:d}\n'
                        for i in range(matrix.shape[0]):
                            shell_str += f'{matrix[i, 0]:18.12e}'
                            for j in range(1, matrix.shape[1]):
                                shell_str += f'{matrix[i, j]:20.12e}'
                            shell_str += '\n'
                        atom_shells[ang_mom].append(shell_str)
                    else:
                        for j in range(1, matrix.shape[1]):
                            shell_str = f'{ang_mom} {nprims:d}  1\n'
                            for i in range(matrix.shape[0]):
                                shell_str += f'{matrix[i, 0]:18.12e}'
                                shell_str += f'{matrix[i, j]:20.12e}'
                                shell_str += '\n'
                            atom_shells[ang_mom].append(shell_str)

                elif ang_mom == 'SP':
                    assert matrix.shape[1] == 3
                    if use_gc_and_sp:
                        shell_str = f'SP {nprims:d}  1\n'
                        for i in range(matrix.shape[0]):
                            shell_str += f'{matrix[i, 0]:18.12e}'
                            shell_str += f'{matrix[i, 1]:20.12e}'
                            shell_str += f'{matrix[i, 2]:20.12e}'
                            shell_str += '\n'
                        atom_shells['SP'].append(shell_str)
                    else:
                        shell_str = f'S {nprims:d}  1\n'
                        for i in range(matrix.shape[0]):
                            shell_str += f'{matrix[i, 0]:18.12e}'
                            shell_str += f'{matrix[i, 1]:20.12e}'
                            shell_str += '\n'
                        atom_shells['S'].append(shell_str)
                        shell_str = f'P {nprims:d}  1\n'
                        for i in range(matrix.shape[0]):
                            shell_str += f'{matrix[i, 0]:18.12e}'
                            shell_str += f'{matrix[i, 2]:20.12e}'
                            shell_str += '\n'
                        atom_shells['P'].append(shell_str)

                elif ang_mom == 'SPD':
                    assert matrix.shape[1] == 4
                    if use_gc_and_sp:
                        shell_str = f'SP {nprims:d}  1\n'
                        for i in range(matrix.shape[0]):
                            shell_str += f'{matrix[i, 0]:18.12e}'
                            shell_str += f'{matrix[i, 1]:20.12e}'
                            shell_str += f'{matrix[i, 2]:20.12e}'
                            shell_str += '\n'
                        atom_shells['SP'].append(shell_str)
                    else:
                        shell_str = f'S {nprims:d}  1\n'
                        for i in range(matrix.shape[0]):
                            shell_str += f'{matrix[i, 0]:18.12e}'
                            shell_str += f'{matrix[i, 1]:20.12e}'
                            shell_str += '\n'
                        atom_shells['S'].append(shell_str)
                        shell_str = f'P {nprims:d}  1\n'
                        for i in range(matrix.shape[0]):
                            shell_str += f'{matrix[i, 0]:18.12e}'
                            shell_str += f'{matrix[i, 2]:20.12e}'
                            shell_str += '\n'
                        atom_shells['P'].append(shell_str)
                    shell_str = f'D {nprims:d}  1\n'
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
