from mpi4py import MPI
import pytest
import sys

from veloxchem.conformergenerator import ConformerGenerator
from veloxchem.molecule import Molecule

try:
    import openmm
    import rdkit
except ImportError:
    pass


class TestConformerGenerator:

    @pytest.mark.skipif("openmm" not in sys.modules or
                        "rdkit" not in sys.modules,
                        reason="openmm or rdkit not available")
    def test_conformergenerator(self):

        molecule = Molecule.read_smiles(
            "CC1([C@@H](N2[C@H](S1)[C@@H](C2=O)NC(=O)CC3=CC=CC=C3)C(=O)O)C")

        conf = ConformerGenerator()
        conf.top_file_name = "mol.top"
        conf.number_of_conformers_to_select = 10
        conf.save_xyz_files = False
        conf.em_tolerance = 1

        conf.generate(molecule)

        if MPI.COMM_WORLD.Get_rank() == 0:
            minimum_potential_energy = conf.global_minimum_energy
            assert (minimum_potential_energy > -69.0 and
                    minimum_potential_energy < -68.0)
