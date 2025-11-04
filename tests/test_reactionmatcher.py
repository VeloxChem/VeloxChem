import pytest
import sys
from pathlib import Path

from veloxchem.molecule import Molecule
from veloxchem.evbdriver import EvbDriver

try:
    import openmm
    import rdkit
except ImportError:
    pass


@pytest.mark.skipif(('openmm' not in sys.modules)
                    or ('rdkit' not in sys.modules),
                    reason='openmm or rdkit not available')
class TestReactionMatcher:

    def run_graph_matcher(
            self,
            rea_smiles=None,
            pro_smiles=None,
            rea=None,
            pro=None,
            breaking_bonds=set(),
            forming_bonds=set(),
    ):
        if rea_smiles is not None:
            rea = []
            for smiles in rea_smiles:
                mol = Molecule.read_smiles(smiles)
                rea.append(mol)

        if pro_smiles is not None:
            pro = []
            for smiles in pro_smiles:
                mol = Molecule.read_smiles(smiles)
                pro.append(mol)

        rea_charges = []
        for mol in rea:
            q = [mol.get_charge() / mol.number_of_atoms()
                 ] * mol.number_of_atoms()
            rea_charges.append(q)

        pro_charges = []
        for mol in pro:
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
            forced_breaking_bonds=breaking_bonds,
            forced_forming_bonds=forming_bonds,
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

        # Needs monomorphism in find forming edges
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

        # Needs monomorphism in find breaking edges
        breaking_bonds, forming_bonds = self.run_graph_matcher(
            ['CCCCC=N'],
            ['CC3CCCN3'],
        )
        assert breaking_bonds == {(1, 9)}
        assert forming_bonds == {(4, 9), (1, 5)}

        # load forcefields

    @pytest.mark.timeconsuming
    def test_copper_complex_1(self):
        data_path = Path(__file__).parent / 'data'
        rea1 = Molecule.read_xyz_file(str(data_path / 'reamatcher_cu-tfe.xyz'))
        rea2 = Molecule.read_xyz_file(str(data_path / 'reamatcher_c-s.xyz'))
        pro1 = Molecule.read_smiles('c1ccccc1[CH]C')
        pro2 = Molecule.read_xyz_file(str(data_path /
                                          'reamatcher_CuOS-tfe.xyz'))

        pro1.set_multiplicity(2)
        pro2.set_multiplicity(2)

        breaking_bonds, forming_bonds = self.run_graph_matcher(
            rea=[rea1, rea2],
            pro=[pro1, pro2],
            breaking_bonds={(72, 77)},
            forming_bonds={(52, 79)},
        )
        # Input is one-indexed, but inner representation is zero-indexed
        assert breaking_bonds == {(71, 76)}
        assert forming_bonds == {(51, 78)}

        breaking_bonds, forming_bonds = self.run_graph_matcher(
            rea=[rea1, rea2],
            pro=[pro1, pro2],
            breaking_bonds={(72, 77)},
        )
        assert breaking_bonds == {(71, 76)}
        assert forming_bonds == {(51, 78)}

        breaking_bonds, forming_bonds = self.run_graph_matcher(
            rea=[rea1, rea2],
            pro=[pro1, pro2],
            forming_bonds={(52, 79)},
        )
        assert breaking_bonds == {(71, 76)}
        assert forming_bonds == {(51, 78)}

        breaking_bonds, forming_bonds = self.run_graph_matcher(
            rea=[rea1, rea2],
            pro=[pro1, pro2],
        )
        assert breaking_bonds == {(71, 76)}
        assert forming_bonds == {(51, 78)}

    @pytest.mark.timeconsuming
    def test_copper_complex_2(self):
        data_path = Path(__file__).parent / 'data'
        rea = Molecule.read_xyz_file(str(data_path / 'reamatcher_3plus.xyz'))

        pro = Molecule.read_xyz_file(str(data_path / 'reamatcher_5plus.xyz'))

        rea.set_charge(1)
        pro.set_charge(1)

        breaking_bonds, forming_bonds = self.run_graph_matcher(
            rea=[rea],
            pro=[pro],
        )
        # Input is one-indexed, but inner representation is zero-indexed
        assert breaking_bonds == {(1, 99)}
        assert forming_bonds == {(42, 99), (36, 101)}
