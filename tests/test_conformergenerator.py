from pathlib import Path
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

        molecule = Molecule.read_xyz_string("""
            41
            xyz
            C    -0.651688  -1.180421   1.366168
            C    -0.492406   0.184620   2.059605
            C     1.024042   0.533559   2.179372
            N     1.298890   1.283508   0.981468
            C     0.307079   2.071282   0.313909
            S    -1.292364   1.457903   0.971649
            C     0.757865   1.684514  -1.011102
            C     1.848950   0.905228  -0.340015
            O     2.801623   0.180952  -0.729915
            N    -0.149953   0.815789  -1.769902
            C     0.296225   0.058533  -2.903641
            O     1.467472   0.203488  -3.347465
            C    -0.628929  -0.917240  -3.565457
            C    -1.005084  -2.012109  -2.604000
            C    -2.266343  -2.015350  -1.984035
            C    -2.599612  -3.020958  -1.071411
            C    -1.673457  -4.020235  -0.760538
            C    -0.413010  -4.017407  -1.364005
            C    -0.077903  -3.017763  -2.282481
            C     1.375784   1.374073   3.378748
            O     1.065536   2.595779   3.419386
            O     2.013858   0.786433   4.466440
            C    -1.228432   0.139705   3.412496
            H    -0.147423  -1.180785   0.376546
            H    -0.208076  -1.989203   1.985821
            H    -1.726288  -1.411541   1.206075
            H     1.655288  -0.382594   2.196082
            H     0.457531   3.153567   0.522512
            H     1.152870   2.526896  -1.619602
            H    -1.177159   0.833900  -1.587634
            H    -0.139686  -1.362860  -4.459210
            H    -1.531366  -0.375309  -3.921910
            H    -2.989198  -1.239665  -2.203985
            H    -3.572558  -3.019769  -0.597396
            H    -1.930713  -4.793625  -0.048526
            H     0.305448  -4.788628  -1.118326
            H     0.904314  -3.022273  -2.738774
            H     2.244827   1.330375   5.289160
            H    -2.280751  -0.188888   3.272574
            H    -1.254614   1.145108   3.882911
            H    -0.725427  -0.569734   4.104181""")

        conf = ConformerGenerator()
        conf.ostream.mute()

        here = Path(__file__).parent
        top_fpath = here / 'data' / 'vlx_conf_gen.top'

        conf.top_file_name = str(top_fpath)
        conf.number_of_conformers_to_select = 10
        conf.save_xyz_files = False
        conf.em_tolerance = 1

        conf_results = conf.generate(molecule)

        if MPI.COMM_WORLD.Get_rank() == 0:
            assert conf_results['energies'][0] == conf.global_minimum_energy
            assert (conf.global_minimum_energy > -69.0 and
                    conf.global_minimum_energy < -68.0)

            if top_fpath.with_suffix('.gro').is_file():
                top_fpath.with_suffix('.gro').unlink()
            if top_fpath.with_suffix('.itp').is_file():
                top_fpath.with_suffix('.itp').unlink()
            if top_fpath.is_file():
                top_fpath.unlink()
