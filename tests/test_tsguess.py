from mpi4py import MPI
from pathlib import Path
import pytest
import sys

from veloxchem.veloxchemlib import mpi_master
from veloxchem.tsguesser import TransitionStateGuesser
from veloxchem.evbdriver import EvbDriver
from veloxchem.optimizationdriver import OptimizationDriver
from veloxchem.molecule import Molecule

try:
    import openmm as mm
except ImportError:
    pass


@pytest.mark.skipif(('openmm' not in sys.modules),
                    reason='openmm not available')
@pytest.mark.timeconsuming
class TestTransitionStateGuesser:

    def test_ts_guesser(self):

        rea = Molecule.read_str("""
        C              1.340494446986         1.458194010980         2.002452328500
        Cl            -0.248712060729         0.603810082277         1.884782585397
        H              1.176565798320         2.498243131144         2.310601520444
        H              1.973676790569         0.954059040738         2.742914700148
        H              1.839289207801         1.441849856114         1.025499053148
        Br             4.316943499142         3.074995880494         2.227046197863"""
                                )

        pro = Molecule.read_str("""
        C              1.644879612388        -0.508351390798        -0.616977056715
        Cl             3.923435352318        -2.563707406826        -1.580365047563
        H              2.132647869077        -0.956075482316         0.257650347399
        H              2.384228922941         0.036800182006        -1.216446463266
        H              1.178719053826        -1.293388510743        -1.224930701103
        Br             0.244603286036         0.757235977330        -0.003426182246"""
                                )

        ref_mm_ts = Molecule.read_str("""
        C              0.287401892289         0.189412663028        -0.181389105630
        Cl            -0.714922847473         1.441402502475        -1.633913908687
        H             -0.506578391806        -0.570890866041        -0.274740014262
        H              1.179217362421         0.137862529013        -0.828734078861
        H              0.201866111086         0.985209360229         0.577838564927
        Br             1.320755564351        -1.098139598962         1.310400014989"""
                                      )
        ref_qm_ts = Molecule.read_str("""
        C             -0.300250466161         0.158721109824         0.132425421981
        Cl            -1.874349934839        -1.112845926481         0.335515890867
        H              0.374806873996        -0.387013493518         0.813286305729
        H             -0.184403283782         0.023551529485        -0.956320205450
        H             -0.827026276510         1.053156785986         0.506230455730
        Br             1.523332262274         1.629693555505        -0.102017809159"""
                                      )
        # check icrmsd between reference and result
        rea.set_charge(-1)
        pro.set_charge(-1)
        ts_guesser = TransitionStateGuesser()
        # ts_guesser.ostream.mute()

        ts_guesser.qm_xcfun = "HF"
        ts_guesser.qm_basis = "STO-3G"
        ts_guesser.do_qm_scan = True

        ts_guesser.save_results_file = False
        ts_guesser.mm_steps = 200
        ts_guesser.lambda_vector = [
            0, 0.2, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.8, 1.0
        ]
        results = ts_guesser.find_transition_state(rea, pro)
        mm_icrmsd = OptimizationDriver.get_ic_rmsd(
            Molecule.read_xyz_string(results['max_mm_xyz']), ref_mm_ts)
        qm_icrmsd = OptimizationDriver.get_ic_rmsd(
            Molecule.read_xyz_string(results['max_qm_xyz']), ref_qm_ts)

        # Generous values based on averages from 100 runs
        max_bond_rms = 0.1
        max_bond_max = 0.2
        max_angle_rms = 10
        max_angle_max = 20

        assert results['max_mm_lambda'] == 0.6
        assert results['max_qm_lambda'] == 0.45

        for icrmsd in [mm_icrmsd, qm_icrmsd]:
            assert icrmsd['bonds']['rms'] < max_bond_rms
            assert icrmsd['bonds']['max'] < max_bond_max
            assert icrmsd['angles']['rms'] < max_angle_rms
            assert icrmsd['angles']['max'] < max_angle_max
