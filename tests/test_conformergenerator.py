from veloxchem.conformergenerator import ConformerGenerator
from veloxchem.molecule import Molecule
from mpi4py import MPI
import sys
import pytest

try:
    import openmm as omm
except ImportError:
    pass


class TestConformerGenerator:
    def get_global_minimum_potentialenergy_refrange(self):
        return [-80, -74]

    @pytest.mark.skipif("openmm" not in sys.modules, reason="openmm not available")
    def test_conformergenerator(self):
        # Test the ConformerGenerator class
        conf = ConformerGenerator()
        conf.molecule = Molecule.read_smiles(
            "CC1([C@@H](N2[C@H](S1)[C@@H](C2=O)NC(=O)CC3=CC=CC=C3)C(=O)O)C"
        )
        conf.top_file_name = "MOL.top"
        conf.number_of_conformers_to_select = 10
        conf.save_xyz = True  # if True then save the xyz file in filtered folder
        conf.em_tolerance_value = 1  # default is 10
        conf.generate()
        conf.show_global_minimum()
        minimum_potential_energy = conf.global_minimum_energy
        if MPI.COMM_WORLD.Get_rank() == 0:
            assert (
                minimum_potential_energy
                > self.get_global_minimum_potentialenergy_refrange()[0]
                and minimum_potential_energy
                < self.get_global_minimum_potentialenergy_range()[1]
            )
