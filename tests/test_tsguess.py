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
        evb = EvbDriver()
        evb.build_ff_from_molecules(rea, pro)
        ts_guesser = TransitionStateGuesser()
        ts_guesser.scf_xcfun = "HF"
        ts_guesser.scf_basis = "STO-3G"
        ts_mol, results = ts_guesser.find_TS(evb)
        folder = Path(__file__).parent / 'data'
        reference_results = evb._load_dict_from_h5(
            folder / "ts_guesser_reference_results.h5")
        self._compare_dict(results, reference_results)

    def _compare_dict(self, dict1, dict2, float_tol=1):

        assert sorted(list(dict1.keys())) == sorted(list(dict2.keys()))

        for key in dict1:

            if key == 'comment' or key == 'mm_energies' or key == 'xyz_geometries' or key == 'final_geometry':
                continue

            val1 = dict1[key]
            val2 = dict2[key]

            type1 = type(val1)
            type2 = type(val2)

            # try to convert val2
            if type1 != type2:
                try:
                    val2 = type1(val2)
                except (ValueError, TypeError):
                    print(
                        f"Type mismatch: {type1} != {type(val2)} for key {key}")
                    assert False

            # compare val1 with val2
            if type1 is dict:
                self._compare_dict(val1, val2)
            elif type1 is float or type1 is np.float64:
                assert abs(val1 - val2) < float_tol
            elif type1 is list or type1 is np.ndarray:
                assert np.allclose(val1, val2, atol=float_tol)
            else:
                assert val1 == val2
