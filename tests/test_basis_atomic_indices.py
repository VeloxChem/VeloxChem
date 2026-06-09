from veloxchem.molecularbasis import MolecularBasis
from veloxchem.molecule import Molecule


class TestBasisAtomicIndices:

    def assert_atomic_indices_match_basis_functions(self, molecule, basis_label):

        basis = MolecularBasis.read(molecule, basis_label, ostream=None)

        atomic_indices = basis.atomic_indices()
        basis_functions = basis.basis_functions()

        assert len(atomic_indices) == len(basis_functions)

    def test_regular_basis(self):

        molstr = """
        O   0.0   0.0   0.0
        H   0.0   1.4   0.0
        H   1.2  -0.7   0.0
        """
        molecule = Molecule.read_molecule_string(molstr, units="bohr")

        self.assert_atomic_indices_match_basis_functions(molecule, "def2-svp")

    def test_mixed_basis(self):

        molstr = """
        O   0.0   0.0   0.0   def2-svpd
        H   0.0   1.4   0.0
        H   1.2  -0.7   0.0
        """
        molecule = Molecule.read_molecule_string(molstr, units="bohr")

        self.assert_atomic_indices_match_basis_functions(molecule, "sto-3g")

    def test_ghost_atom(self):

        molstr = """
        O      0.0   0.0   0.0
        Bq_H   0.0   1.4   0.0
        H      1.2  -0.7   0.0
        """
        molecule = Molecule.read_molecule_string(molstr, units="bohr")

        self.assert_atomic_indices_match_basis_functions(molecule, "def2-svp")

    def test_ecp(self):

        molstr = """
        Au     0.0   0.0   0.0
        Au     0.0   0.0   5.0
        """
        molecule = Molecule.read_molecule_string(molstr, units="bohr")

        self.assert_atomic_indices_match_basis_functions(molecule, "def2-svp")
