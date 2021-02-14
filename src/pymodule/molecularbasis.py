from os import environ
from os import walk
from pathlib import Path

from .veloxchemlib import MolecularBasis
from .veloxchemlib import AtomBasis
from .veloxchemlib import BasisFunction
from .veloxchemlib import ChemicalElement
from .veloxchemlib import to_angular_momentum
from .inputparser import InputParser
from .outputstream import OutputStream
from .errorhandler import assert_msg_critical


@staticmethod
def _MolecularBasis_read(mol, basis_name, basis_path='.', ostream=None):
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

    if ostream is None:
        ostream = OutputStream(None)

    err_gc = "MolcularBasis.read: "
    err_gc += "General contraction currently is not supported"

    # searching order:
    # 1. given basis_path
    # 2. current directory
    # 3. VLXBASISPATH

    basis_file = Path(basis_path, basis_name.upper())

    if not basis_file.is_file() and basis_path != '.':
        basis_file = Path('.', basis_name.upper())

    if not basis_file.is_file() and 'VLXBASISPATH' in environ:
        basis_file = Path(environ['VLXBASISPATH'], basis_name.upper())

    basis_info = "Reading basis set: " + str(basis_file)
    ostream.print_info(basis_info)
    ostream.print_blank()

    basis_dict = InputParser(str(basis_file)).input_dict

    assert_msg_critical(
        basis_name.upper() == basis_dict['basis_set_name'].upper(),
        "MolecularBasis.read: Inconsistent basis set name")

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


@staticmethod
def _MolecularBasis_get_avail_basis(element_label):
    """
    Gets the names of available basis sets for an element.

    :param element_label:
        The label of the chemical element.

    :return:
        The tuple of basis sets.
    """

    avail_basis = set()

    basis_path = environ['VLXBASISPATH']

    for root, dirs, files in walk(basis_path, topdown=True):
        for filename in files:
            basis_file = Path(basis_path, filename)
            basis_dict = InputParser(str(basis_file)).input_dict
            for key in list(basis_dict.keys()):
                if 'atombasis_' in key:
                    elem = key.replace('atombasis_', '')
                    if element_label.upper() == elem.upper():
                        avail_basis.add(filename)
        # skip subfolder
        break

    return sorted(list(avail_basis))


MolecularBasis.read = _MolecularBasis_read
MolecularBasis.get_avail_basis = _MolecularBasis_get_avail_basis
