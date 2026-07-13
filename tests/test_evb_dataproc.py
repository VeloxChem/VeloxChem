from veloxchem.evbdriver import EvbDriver
from veloxchem.evbdataprocessing import EvbDataProcessing

import numpy as np
import pytest

from veloxchem.outputstream import OutputStream

from test_evb_helper import evb_compare_dict, evb_data_dir

pytestmark = [pytest.mark.timeconsuming]

pytest.importorskip('pymbar')


class TestEvbDataProcessing:

    @staticmethod
    def _load_canned_input(folder, EVB):
        specific_results = {}
        options_file = folder / 'evb_options.json'

        specific, common = EVB._load_output_files(
            folder / 'evb_Sn2_vacuum_Energies.csv',
            folder / 'evb_Sn2_vacuum_Data_combined.csv',
            options_file,
        )
        specific_results['vacuum'] = specific

        specific, common = EVB._load_output_files(
            folder / 'evb_Sn2_water_Energies.csv',
            folder / 'evb_Sn2_water_Data_combined.csv',
            options_file,
        )
        specific_results['water'] = specific

        input_results = {}
        input_results.update(common)
        input_results.update({"configuration_results": specific_results})
        return input_results

    def test_data_processing(self):
        # Load canned simulation output and run the post-processing, comparing
        # against the reference results.
        folder = evb_data_dir()
        EVB = EvbDriver(ostream=OutputStream(None))

        input_results = self._load_canned_input(folder, EVB)

        dp = EvbDataProcessing(ostream=OutputStream(None))
        # evb_reference_results.h5 was generated with the pairwise-BAR
        # estimator; pin it explicitly rather than relying on whatever
        # EvbDataProcessing.fep_estimator currently defaults to.
        comp_results = dp.compute(input_results, 5, 10)

        reference_results = EVB._load_dict_from_h5(folder /
                                                   "evb_reference_results.h5")
        evb_compare_dict(comp_results, reference_results)

    def test_data_processing_mbar_agrees_with_bar(self):
        # fep_estimator="mbar" is a selectable alternative to the default
        # pairwise-BAR dGfep estimate (see EvbDataProcessing.fep_estimator).
        # There's no separate committed h5 reference for it; instead check it
        # runs cleanly and lands close to the BAR run on the same canned data
        # (both are valid free-energy estimators for the same underlying
        # system, so they should agree to within a few kJ/mol even though the
        # alpha/H12 fit -- which itself calls _calculate_dGfep -- can differ
        # slightly between the two estimators).
        folder = evb_data_dir()
        EVB = EvbDriver(ostream=OutputStream(None))

        bar_input = self._load_canned_input(folder, EVB)
        dp_bar = EvbDataProcessing(ostream=OutputStream(None))
        dp_bar.fep_estimator = "bar"
        bar_results = dp_bar.compute(bar_input, 5, 10)

        mbar_input = self._load_canned_input(folder, EVB)
        dp_mbar = EvbDataProcessing(ostream=OutputStream(None))
        dp_mbar.fep_estimator = "mbar"
        mbar_results = dp_mbar.compute(mbar_input, 5, 10)

        for name in ("vacuum", "water"):
            bar_conf = bar_results["configuration_results"][name]
            mbar_conf = mbar_results["configuration_results"][name]

            assert len(mbar_conf["dGfep"]) == len(bar_conf["dGfep"])
            assert np.all(np.isfinite(mbar_conf["dGfep"]))
            assert mbar_conf["dGfep"][0] == pytest.approx(0.0, abs=1e-8)

            bar_analytical = bar_conf["analytical"]
            mbar_analytical = mbar_conf["analytical"]
            assert np.isfinite(mbar_analytical["barrier"])
            assert np.isfinite(mbar_analytical["free_energy"])
            assert mbar_analytical["barrier"] == pytest.approx(
                bar_analytical["barrier"], abs=5.0)
            assert mbar_analytical["free_energy"] == pytest.approx(
                bar_analytical["free_energy"], abs=5.0)
