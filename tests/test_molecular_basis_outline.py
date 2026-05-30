from veloxchem.veloxchemlib import newints
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis


class TestMolecularBasisOutline:

    def water_sto3g(self):

        xyz = """3

        O    0.000000    0.000000    0.000000
        H    0.000000    0.000000    0.950000
        H    0.895670    0.000000   -0.316663
        """
        mol = Molecule.read_xyz_string(xyz)
        bas = MolecularBasis.read(mol, "sto-3g", ostream=None)
        return mol, bas

    def test_counts(self):

        _, bas = self.water_sto3g()
        outline = newints.MolecularBasisOutline(bas)

        # water / STO-3G: O has shells [s, s, p] (3), each H has [s] (1)
        assert outline.number_of_atoms() == 3
        assert outline.number_of_basis_functions() == 5
        # angular-expanded dimension: O -> 1+1+3, H -> 1, H -> 1 = 7 AOs
        assert outline.number_of_atomic_orbitals() == 7

    def test_global_arrays(self):

        _, bas = self.water_sto3g()
        outline = newints.MolecularBasisOutline(bas)

        # indices advance by (2l+1) per shell, atom-major
        assert outline.indices() == [0, 1, 2, 5, 6]
        assert outline.angular_momenta() == [0, 0, 1, 0, 0]
        assert outline.atom_indices() == [0, 0, 0, 1, 2]

    def test_per_atom_indices(self):

        _, bas = self.water_sto3g()
        outline = newints.MolecularBasisOutline(bas)

        assert outline.basis_function_indices(0) == [0, 1, 2]  # O
        assert outline.basis_function_indices(1) == [5]        # H_1
        assert outline.basis_function_indices(2) == [6]        # H_2

    def test_default_is_empty(self):

        outline = newints.MolecularBasisOutline()
        assert outline.number_of_atoms() == 0
        assert outline.number_of_basis_functions() == 0
        assert outline.number_of_atomic_orbitals() == 0

    def test_equality(self):

        _, bas = self.water_sto3g()
        a = newints.MolecularBasisOutline(bas)
        b = newints.MolecularBasisOutline(bas)
        assert a == b
        assert not (a == newints.MolecularBasisOutline())
