from pathlib import Path
import numpy as np

from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.oneeints import compute_linear_momentum_integrals


class TestOneElecIntsLinearMomentum:

    def get_molecule_and_basis(self):

        xyz_string = """6
        xyz
        O   -0.1858140  -1.1749469   0.7662596
        H   -0.1285513  -0.8984365   1.6808606
        H   -0.0582782  -0.3702550   0.2638279
        O    0.1747051  11.1050002  -0.7244430
        H   -0.5650842  11.3134964  -1.2949455
        H    0.9282185  11.0652990  -1.3134026
        """
        basis_label = 'def2-svp'

        mol = Molecule.read_xyz_string(xyz_string)
        bas = MolecularBasis.read(mol, basis_label, ostream=None)

        return mol, bas

    def test_linear_momentum_integrals(self):

        here = Path(__file__).parent
        ref_linmom = np.load(str(here / 'data' / 'onee-svp-linmom.npy'))

        mol, bas = self.get_molecule_and_basis()
        linmom = compute_linear_momentum_integrals(mol, bas)
        for comp in range(3):
            assert np.max(np.abs(linmom[comp] - ref_linmom[comp])) < 1.0e-12
