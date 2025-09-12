import pytest
import sys

from veloxchem.veloxchemlib import mpi_master
from veloxchem.veloxchemlib import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.hydrogenbdedriver import HydrogenBdeDriver

try:
    import rdkit
except ImportError:
    pass


@pytest.mark.hydrogen_bde
class TestHydrogenBdeDriver:

    def run_hydrogen_bde(self, molecules):
        hbde = HydrogenBdeDriver()
        hbde.mute_output = True
        hbde.compute(molecules)
        if hbde._rank == mpi_master():
            assert abs(hbde.mols_bdes_list[0]["unique_BDEs_au"][0] - 0.1752835) < 1.0e-4

    @pytest.mark.skipif("rdkit" not in sys.modules,
                        reason="rdkit not available")
    def test_hydrogen_bde_driver(self):
        self.run_hydrogen_bde(Molecule.read_smiles('C'))
