import numpy as np

from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.veloxchemlib import GtoBlock


class TestTabulaScreeningData:
    """Tests for CGtoBlock.screening_data — the per-contracted-GTO screening
    bounds (max contraction-coefficient magnitude, min/max exponent)."""

    def get_block(self, angmom, npgtos):

        mol = Molecule.read_str(
            "O 0.000  0.000 -1.000\n"
            "H 0.000  1.400 -2.100\n"
            "H 0.000 -1.400 -2.100", 'au')
        bas = MolecularBasis.read(mol, 'DEF2-SVP', ostream=None)
        return GtoBlock(bas, mol, angmom, npgtos)

    def check_block(self, angmom, npgtos):

        block = self.get_block(angmom, npgtos)
        nc = block.number_of_basis_functions()
        npg = block.number_of_primitives()

        # primitive-major SoA -> reshape to (primitive, CGTO)
        exps = np.array(block.exponents()).reshape(npg, nc)
        norms = np.array(block.normalization_factors()).reshape(npg, nc)

        # the screening bounds of each contracted GTO are its column reductions
        for c in range(nc):
            data = block.screening_data(c)
            assert np.isclose(data.max_coefficient,
                              np.max(np.abs(norms[:, c])), 0.0, 1.0e-13)
            assert np.isclose(data.min_exponent,
                              np.min(exps[:, c]), 0.0, 1.0e-13)
            assert np.isclose(data.max_exponent,
                              np.max(exps[:, c]), 0.0, 1.0e-13)

    def test_contracted_s_block(self):

        # s functions, 3 primitives each
        self.check_block(0, 3)

    def test_single_primitive_s_block(self):

        # s functions, 1 primitive each -> min exponent == max exponent
        block = self.get_block(0, 1)
        for c in range(block.number_of_basis_functions()):
            data = block.screening_data(c)
            assert data.min_exponent == data.max_exponent
        self.check_block(0, 1)
