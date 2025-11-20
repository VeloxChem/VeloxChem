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
        ts_guesser.ostream.mute()
        ts_guesser.scf_xcfun = "HF"
        ts_guesser.scf_basis = "STO-3G"
        results = ts_guesser.find_TS(rea, pro)
        folder = Path(__file__).parent / 'data'

        reference_results = EvbDriver._load_dict_from_h5(
            folder / "ts_guesser_reference_results.h5")

        assert results['breaking_bonds'] == {
            tuple(reference_results['breaking_bonds'][0])
        }
        assert results['forming_bonds'] == {
            tuple(reference_results['forming_bonds'][0])
        }
        assert results['max_scf_lambda'] == reference_results['max_scf_lambda']

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
