import pytest
import sys

from veloxchem.veloxchemlib import mpi_master
from veloxchem.veloxchemlib import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.atombdedriver import AtomBdeDriver

try:
    import rdkit
except ImportError:
    pass


class TestAtomBdeDriver:

    def run_atom_bde(self, molecules):
        bde = AtomBdeDriver()
        bde.mute_output = True
        bde.save_files = False
        bde.compute(molecules)
        if bde._rank == mpi_master():
            assert abs(bde.mols_bdes_list[0]["unique_BDEs_hartree"][0] -
                       0.17657796026960781) < 1.0e-4

    @pytest.mark.skipif("rdkit" not in sys.modules,
                        reason="rdkit not available")
    def test_atom_bde_driver(self):
        self.run_atom_bde(Molecule.read_smiles('C'))
