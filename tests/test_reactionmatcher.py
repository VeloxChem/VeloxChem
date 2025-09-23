import pytest
import sys

from veloxchem.molecule import Molecule
from veloxchem.evbdriver import EvbDriver

try:
    import openmm
    import rdkit
except ImportError:
    pass


@pytest.mark.skipif(('openmm' not in sys.modules) or
                    ('rdkit' not in sys.modules),
                    reason='openmm or rdkit not available')
class TestReactionMatcher:

    @staticmethod
    def run_graph_matcher(
            rea_smiles,
            pro_smiles,
            breaking_bonds=set(),
    ):
        rea = []
        rea_charges = []
        for smiles in rea_smiles:
            mol = Molecule.read_smiles(smiles)
            rea.append(mol)
            q = [mol.get_charge() / mol.number_of_atoms()
                ] * mol.number_of_atoms()
            rea_charges.append(q)
        pro = []
        pro_charges = []
        for smiles in pro_smiles:
            mol = Molecule.read_smiles(smiles)
            pro.append(mol)
            q = [mol.get_charge() / mol.number_of_atoms()
                ] * mol.number_of_atoms()
            pro_charges.append(q)

        evb = EvbDriver()
        evb.ffbuilder.reparameterize_bonds = False
        evb.ffbuilder.optimize_ff = False
        evb.build_ff_from_molecules(
            rea,
            pro,
            reactant_partial_charges=rea_charges,
            product_partial_charges=pro_charges,
            breaking_bonds=breaking_bonds,
        )
        return evb.breaking_bonds, evb.forming_bonds

    def test_reactionmatcher(self):
        breaking_bonds, forming_bonds = self.run_graph_matcher(
            ['CCO'],
            ['O', 'C=C'],
        )
        assert breaking_bonds == {(1, 2), (0, 4)}
        assert forming_bonds == {(2, 4)}

        breaking_bonds, forming_bonds = self.run_graph_matcher(
            ['[Cl-]', 'CCBr'],
            ['[Br-]', 'CCCl'],
        )
        assert breaking_bonds == {(2, 3)}
        assert forming_bonds == {(0, 2)}

        breaking_bonds, forming_bonds = self.run_graph_matcher(
            ['C#C', 'CN=[N+]=[N-]'],
            ['CN1N=NC=C1'],
        )
        assert breaking_bonds == set()
        assert forming_bonds == {(1, 7), (0, 5)}

        breaking_bonds, forming_bonds = self.run_graph_matcher(
            ['C1=CCCCC1'],
            ['C=C', 'C=CC=C'],
        )
        assert breaking_bonds == {(2, 3), (4, 5)}
        assert forming_bonds == set()

        breaking_bonds, forming_bonds = self.run_graph_matcher(
            ['C1CC=CCC1'],
            ['C=C', 'C=CC=C'],
        )
        assert breaking_bonds == {(0, 1), (4, 5)}
        assert forming_bonds == set()

        breaking_bonds, forming_bonds = self.run_graph_matcher(
            ['O', 'C(C(=O)O)([O-])'],
            ['O', 'C(=C(O)[O-])O'],
            breaking_bonds={(1, 3)},
        )
        assert breaking_bonds == {(3, 8), (0, 2)}
        assert forming_bonds == {(0, 8), (2, 7)}

        breaking_bonds, forming_bonds = self.run_graph_matcher(
            ['CCCCC=N'],
            ['CC3CCCN3'],
        )
        assert breaking_bonds == {(1, 9)}
        assert forming_bonds == {(4, 9), (1, 5)}
