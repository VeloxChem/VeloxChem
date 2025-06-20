from pathlib import Path
import numpy as np
import pytest
import sys

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

    @pytest.mark.skipif('openmm' not in sys.modules,
                        reason='openmm not available')
    @pytest.mark.timeconsuming
    def test_ts_guesser(self):
        rea = Molecule.read_smiles("C1CCC=C1")
        pro = Molecule.read_smiles("C=CC1CC1")
        ts_guesser = TransitionStateGuesser()
        ts_guesser.scf_xcfun = "HF"
        ts_guesser.scf_basis = "STO-3G"
        ts_mol, results = ts_guesser.find_TS(rea, pro)
        folder = Path(__file__).parent / 'data'

        reference_results = EvbDriver._load_dict_from_h5(
            folder / "ts_guesser_reference_results.h5")

        assert results['broken_bonds'] == {
            tuple(reference_results['broken_bonds'][0])
        }
        assert results['formed_bonds'] == {
            tuple(reference_results['formed_bonds'][0])
        }
        assert results['final_lambda'] == reference_results['final_lambda']
