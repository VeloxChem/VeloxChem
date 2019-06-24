import os

from .veloxchemlib import MolecularBasis
from .veloxchemlib import AtomBasis
from .veloxchemlib import BasisFunction
from .veloxchemlib import ChemicalElement
from .veloxchemlib import to_angular_momentum
from .inputparser import InputParser
from .outputstream import OutputStream
from .errorhandler import assert_msg_critical


@staticmethod
def _MolecularBasis_read(mol,
                         basis_name,
                         basis_path='.',
                         ostream=OutputStream()):
    """
    Reads AO basis set from file.

    :param mol:
        The molecule.
    :param basis_name:
        Name of the basis set.
    :param basis_path:
        Path to the basis set.
    :param ostream:
        The outputstream.

    :return:
        The AO basis set.
    """

    err_gc = "MolcularBasis.read_file: "
    err_gc += "General contraction currently is not supported"

    # searching order:
    # 1. given basis_path
    # 2. current directory
    # 3. VLXBASISPATH

    basis_file = os.path.join(basis_path, basis_name.upper())

    if not os.path.isfile(basis_file) and basis_path != '.':
        basis_file = os.path.join('.', basis_name.upper())

    if not os.path.isfile(basis_file) and 'VLXBASISPATH' in os.environ:
        basis_file = os.path.join(os.environ['VLXBASISPATH'],
                                  basis_name.upper())

    basis_info = "Reading basis set: " + basis_file
    ostream.print_info(basis_info)
    ostream.print_blank()

    basis_dict = InputParser(basis_file).input_dict

    assert_msg_critical(
        basis_name.upper() == basis_dict['basis_set_name'].upper(),
        "MolecularBasis.read_file: Inconsistent basis set name")

    mol_basis = MolecularBasis()

    elem_comp = mol.get_elemental_composition()

    for elem_id in elem_comp:

        elem = ChemicalElement()
        err = elem.set_atom_type(elem_id)
        assert_msg_critical(err, "ChemicalElement.set_atom_type")

        basis_key = 'atombasis_{}'.format(elem.get_name().lower())
        basis_list = [entry for entry in basis_dict[basis_key]]

        atom_basis = AtomBasis()

        while basis_list:
            shell_title = basis_list.pop(0).split()
            assert_msg_critical(
                len(shell_title) == 3,
                "Basis set parser (shell): {}".format(' '.join(shell_title)))

            angl = to_angular_momentum(shell_title[0])
            npgto = int(shell_title[1])
            ncgto = int(shell_title[2])

            assert_msg_critical(ncgto == 1, err_gc)

            expons = [0.0] * npgto
            coeffs = [0.0] * npgto * ncgto

            for i in range(npgto):
                prims = basis_list.pop(0).split()
                assert_msg_critical(
                    len(prims) == ncgto + 1,
                    "Basis set parser (primitive): {}".format(' '.join(prims)))

                expons[i] = float(prims[0])
                for k in range(ncgto):
                    coeffs[k * npgto + i] = float(prims[k + 1])

            bf = BasisFunction(expons, coeffs, angl)
            bf.normalize()

            atom_basis.add_basis_function(bf)

        atom_basis.set_elemental_id(elem_id)

        mol_basis.add_atom_basis(atom_basis)

    basis_label = basis_dict['basis_set_name'].upper()

    mol_basis.set_label(basis_label)

    return mol_basis


MolecularBasis.read = _MolecularBasis_read
