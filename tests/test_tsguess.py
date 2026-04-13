from mpi4py import MPI
from pathlib import Path
import pytest
import sys

from veloxchem.veloxchemlib import mpi_master
from veloxchem.tsguesser import TransitionStateGuesser
from veloxchem.evbdriver import EvbDriver
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

        self.clean_up_files()

        rea = Molecule.read_smiles("C1CCC=C1")
        pro = Molecule.read_smiles("C=CC1CC1")
        ts_guesser = TransitionStateGuesser()
        # ts_guesser.ostream.mute()
        # ts_guesser.qm_xcfun = "HF"
        # ts_guesser.qm_basis = "STO-3G"
        ts_guesser.mm_steps = 200
        ts_guesser.lambda_vector = [0, 0.2, 0.4, 0.45, 0.5, 0.55, 0.6, 0.8, 1.0]
        # ts_guesser.do_qm_scan = True
        results = ts_guesser.find_transition_state(rea, pro)

        data_path = Path(__file__).parent / 'data'
        ref_ts_guesser = TransitionStateGuesser()
        reference_results = ref_ts_guesser.load_results(
            data_path / "ts_guesser_reference_results.h5")

        assert results['breaking_bonds'] == reference_results['breaking_bonds']
        assert results['forming_bonds'] == reference_results['forming_bonds']
        assert round(float(results['max_mm_lambda']),
                     2) == round(float(reference_results['max_mm_lambda']), 2)

        self.clean_up_files()

    def clean_up_files(self):

        fnames = [
            'combined_PRO_0',
            'combined_REA_0',
            'ts_results.h5',
        ]

        fpaths = [Path(f) for f in fnames]

        dir_path = Path('ts_data')
        if dir_path.is_dir():
            fpaths += [x for x in dir_path.iterdir() if x.is_file()]

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            for p in fpaths:
                if p.is_file():
                    p.unlink()
        MPI.COMM_WORLD.barrier()
