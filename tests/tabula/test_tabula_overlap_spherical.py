import numpy as np

from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.veloxchemlib import GtoBlock, GtoPairBlock, OverlapDriver
from veloxchem.tabulalib import tabula_overlap_spherical


class TestTabulaOverlapSpherical:
    """Step (d) of the overlap recursion — the Cartesian-to-spherical
    assembly. With it the recursion (steps a-d) computes the full overlap
    integral; for an all-s system the spherical component order is
    unambiguous, so the result is validated directly against VeloxChem's
    OverlapDriver."""

    def check_all_s(self, xyz):

        mol = Molecule.read_str(xyz, 'au')
        bas = MolecularBasis.read(mol, 'STO-3G', ostream=None)

        # VeloxChem reference overlap matrix
        reference = OverlapDriver().compute(mol, bas).full_matrix().to_numpy()

        # Tabula — the s-block paired with itself
        block = GtoBlock(bas, mol, 0, 3)
        pair_block = GtoPairBlock(block, block)

        ncgtos = block.number_of_basis_functions()
        orbital_indices = block.orbital_indices()

        spherical = tabula_overlap_spherical(pair_block)
        assert spherical.shape == (1, ncgtos * ncgtos)

        # contracted pair ij = i*ncgtos + j -> overlap of CGTO i with CGTO j;
        # orbital_indices[k+1] is CGTO k's global AO index
        for i in range(ncgtos):
            for j in range(ncgtos):
                ao_i = orbital_indices[i + 1]
                ao_j = orbital_indices[j + 1]
                assert np.isclose(spherical[0, i * ncgtos + j],
                                  reference[ao_i, ao_j], 0.0, 1.0e-10)

    def test_h2(self):

        self.check_all_s("H 0.0 0.0 0.0\nH 0.0 0.0 1.4")

    def test_h4_chain(self):

        self.check_all_s("H 0.0 0.0 0.0\nH 0.0 0.0 1.4\n"
                         "H 0.0 0.0 2.8\nH 0.0 0.0 4.2")
