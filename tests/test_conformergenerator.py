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
            35
            xyz
            O        0.771538    -1.238771     1.853844
            C       -0.051511    -2.371837     1.838116
            C        0.080738    -3.105573     0.503592
            C       -0.452963    -2.278502    -0.633378
            C        0.210574    -1.289480    -1.376554
            C        1.650943    -0.906542    -1.157606
            N       -0.593126    -0.645719    -2.284726
            C       -1.866671    -1.109323    -2.277966
            S       -2.056889    -2.321565    -1.136560
            C       -0.172170     0.465004    -3.150326
            C       -0.032037     1.742573    -2.348207
            C       -1.172012     2.244535    -1.700591
            N       -1.097569     3.363469    -0.938320
            C        0.076432     4.022228    -0.778290
            C        0.127848     5.270537     0.037519
            N        1.189997     3.580738    -1.415605
            C        1.177828     2.467573    -2.197055
            N        2.400793     2.085549    -2.838517
            H        0.626128    -0.805154     2.734474
            H       -1.112528    -2.075337     2.000667
            H        0.251935    -3.061146     2.657152
            H        1.135881    -3.397326     0.326209
            H       -0.503738    -4.048934     0.560625
            H        2.170602    -0.796876    -2.130301
            H        1.698175     0.052909    -0.602999
            H        2.204903    -1.671502    -0.581176
            H       -2.657849    -0.726148    -2.909669
            H       -0.932640     0.639688    -3.943441
            H        0.739304     0.170832    -3.695975
            H       -2.125234     1.741525    -1.798694
            H       -0.622099     5.224703     0.855274
            H       -0.091738     6.144399    -0.610381
            H        1.136121     5.393202     0.486503
            H        3.272919     2.621633    -2.632610
            H        2.471830     1.340004    -3.560309
        """)
        molecule.set_charge(1.0)

        conf = ConformerGenerator()
        conf.use_gromacs_files = False
        conf.ostream.mute()

        here = Path(__file__).parent
        top_fpath = here / 'data' / 'vlx_conf_gen.top'

        conf.top_file_name = str(top_fpath)
        conf.save_xyz_files = False

        conf_results = conf.generate(molecule)

        if MPI.COMM_WORLD.Get_rank() == 0:
            assert conf_results['energies'][0] == conf.global_minimum_energy
            assert (conf.global_minimum_energy > -63.0 and
                    conf.global_minimum_energy < -62.0)

            if conf.use_gromacs_files:
                for suffix in ['.top', '.itp', '.gro']:
                    top_fpath.with_suffix(suffix).unlink()
            else:
                for suffix in ['.xml', '.pdb']:
                    top_fpath.with_suffix(suffix).unlink()
