import hashlib
from pathlib import Path

import pytest
from veloxchem import VLX_BASIS_PATH
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.molecule import Molecule
from veloxchem.veloxchemlib import AtomBasis, BasisFunction, is_mpi_master


def build_basis_function(expons, coeffs, angl):

    bf = BasisFunction(expons, coeffs, angl)
    bf.normalize()

    return bf


def test_basis_read(tmp_path):

    mol_text = """
        C        0.0000        0.0000       -1.2436
        C       -0.0000        0.0000        1.2436
        H        0.0000        1.7272       -2.3230
        H        0.0000       -1.7272       -2.3230
        H       -0.0000        1.7272        2.3230
        H       -0.0000       -1.7272        2.3230
    """

    basis_text = """@BASIS_SET TESTBASIS
        @ATOMBASIS H
        S 2  1
        3.425250910000e+00  1.543289700000e-01
        6.239137300000e-01  5.353281400000e-01
        @END
        @ATOMBASIS C
        S 1  1
        7.161683700000e+01  1.543289700000e-01
        P 1  1
        2.941249400000e+00  1.559162700000e-01
        @END
    """

    basis_name = "TESTBASIS"

    if is_mpi_master():

        molecule = Molecule.read_str(mol_text, units="angstrom")

        fname = Path(tmp_path, basis_name)
        with fname.open("w") as f_basis:
            f_basis.write(basis_text)
        basis = MolecularBasis.read(molecule, basis_name, tmp_path)

        ref_basis = MolecularBasis()

        atom_basis = AtomBasis()
        atom_basis.add_basis_function(
            build_basis_function(
                [3.425250910000e00, 6.239137300000e-01],
                [1.543289700000e-01, 5.353281400000e-01],
                0,
            )
        )
        atom_basis.set_elemental_id(1)
        ref_basis.add_atom_basis(atom_basis)

        atom_basis = AtomBasis()
        atom_basis.add_basis_function(
            build_basis_function([7.161683700000e01], [1.543289700000e-01], 0)
        )
        atom_basis.add_basis_function(
            build_basis_function([2.941249400000e00], [1.559162700000e-01], 1)
        )
        atom_basis.set_elemental_id(6)
        ref_basis.add_atom_basis(atom_basis)

        ref_basis.set_label(basis_name)

        assert basis == ref_basis


@pytest.mark.parametrize(
    "basis",
    filter(
        lambda x: "RI" not in x.name,
        [
            _
            for _ in VLX_BASIS_PATH.iterdir()
            if _.is_file() and _.name != "MIN-CC-PVDZ"
        ],
    ),
    ids=(lambda x: x.name),
)
def test_basis_md5(basis):
    # get reference md5

    with basis.open("r") as f:
        # read lines
        lines = f.readlines()
        expected_md5 = lines[-1].strip()

    # verify md5

    calculated_md5 = hashlib.md5("".join(lines[:-1]).encode("utf-8")).hexdigest()

    assert expected_md5 == calculated_md5, f"Testing {basis.name}"


@pytest.mark.parametrize(
    "elem_id,basis_sets",
    [
        (
            "H",
            ["6-311++G(2D,2P)", "AUG-CC-PVTZ", "DEF2-SVPD", "MIN-CC-PVDZ", "STO-6G",],
        ),
        (
            "c",
            ["6-31G", "6-31G(2DF,P)", "AUG-CC-PVDZ", "D-AUG-CC-PVQZ", "SADLEJ-PVTZ",],
        ),
        ("Zn", ["CC-PVTZ", "MIN-CC-PVDZ", "STO-3G"]),
    ],
    ids=["H", "C", "Zn"],
)
def test_avail_basis(elem_id, basis_sets):
    for bas in basis_sets:
        assert bas in MolecularBasis.get_avail_basis(elem_id)
