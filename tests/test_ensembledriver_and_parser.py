from pathlib import Path
import pytest
import numpy as np

from veloxchem.veloxchemlib import mpi_master
from veloxchem.ensembleparser import EnsembleParser
from veloxchem.ensembledriver import EnsembleDriver

pytest.importorskip("MDAnalysis")
pytest.importorskip("pyframe")


@pytest.mark.timeconsuming
class TestEnsembleDriverOptions:

    def test_compute_with_scf_and_property_options(self, tmp_path):
        """
        Regression test for notebook-style ensemble API:
        - compute(..., scf_options=..., property_options=...)
        - property='absorption' + nstates routes to TD-DFT/LR path
        - PE potfile is frame-specific for each snapshot
        - compare SCF energies and LR eigenvalues to validated reference output
        - compare generated PE+NPE .pot files against reference files line-by-line
        """
        data_dir = Path(__file__).parent / "data"
        traj = str(data_dir / "alpha-helix-acetone-water.xtc")
        top = str(data_dir / "alpha-helix-acetone-water.tpr")

        ens_parser = EnsembleParser()
        ens_parser.ostream.mute()

        ensemble = ens_parser.structures(
            trajectory_file=traj,
            topology_file=top,
            qm_region="resname LIG",
            num_snapshots=2,
            pe_cutoff=3.0,
            npe_cutoff=5.0,
        )

        ens_drv = EnsembleDriver()
        ens_drv.ostream.mute()

        ens_drv.set_env_models(
            pe_model=["CP3", "SEP"],
            npe_model=["ff19sb", "tip3p"],
        )

        scf_options = {
            "scf_type": "restricted",
            "conv_thresh": 1.0e-7,
            "max_iter": 150,
            "xcfun": "CAM-B3LYP",
            "grid_level": 4,
        }

        property_options = {
            "property": "absorption",
            "nstates": 6,
            "nto": True,
        }

        potdir = data_dir
        results = ens_drv.compute(
            ensemble,
            basis_set="6-31G",
            scf_options=scf_options,
            property_options=property_options,
            potdir=potdir,
        )

        assert "scf_all" in results
        assert "rsp_all" in results
        assert len(results["scf_all"]) == 2
        assert len(results["rsp_all"]) == 2

        # Reference values from validated output (frames 0 and 100)
        ref_frames = [0, 100]
        ref_scf = np.array([
            -193.055643527,
            -193.066013148,
        ])
        ref_eigs = {
            0: np.array([0.160392, 0.303838, 0.321107, 0.32427, 0.344735, 0.360693]),
            100: np.array([0.158489, 0.289851, 0.313072, 0.320164, 0.351021, 0.362582]),
        }

        got_frames = [int(frame) for frame, _ in results["scf_all"]]
        assert got_frames == ref_frames

        if ens_drv.rank == mpi_master():
            got_scf = np.array([float(scf_res["scf_energy"]) for _, scf_res in results["scf_all"]])
            np.testing.assert_allclose(got_scf, ref_scf, rtol=0.0, atol=2.0e-6)

            for frame, scf_res in results["scf_all"]:
                expected_name = f"pe_frame_{int(frame):06d}.pot"
                expected_path = potdir / expected_name

                # Each snapshot must use its own PE file
                assert expected_path.is_file()
                assert str(scf_res.get("potfile", "")).endswith(expected_name)

                # Reference .pot regression for mixed PE+NPE writer stability
                reference_name = f"ensemble_alphahelix_reference_{expected_name}"
                reference_path = data_dir / reference_name
                assert reference_path.is_file(), f"Missing reference file: {reference_path}"

                generated_lines = expected_path.read_text().splitlines()
                reference_lines = reference_path.read_text().splitlines()
                assert generated_lines == reference_lines, (
                    f"PE+NPE pot mismatch for frame {frame}: {expected_name}"
                )

            for frame, rsp_res in results["rsp_all"]:
                frame = int(frame)
                eig = np.array(rsp_res.get("eigenvalues", []), dtype=float)
                assert int(rsp_res.get("number_of_states", 0)) == 6
                assert eig.shape == (6,)
                np.testing.assert_allclose(
                    eig,
                    ref_eigs[frame],
                    rtol=0.0,
                    atol=2.0e-5,
                )

            # cleanup
            (data_dir / ".alpha-helix-acetone-water.xtc_offsets.npz").unlink(missing_ok=True)
            (data_dir / "pe_frame_000000.pot").unlink(missing_ok=True)
            (data_dir / "pe_frame_000000.json").unlink(missing_ok=True)
            (data_dir / "pe_frame_000100.pot").unlink(missing_ok=True)
            (data_dir / "pe_frame_000100.json").unlink(missing_ok=True)

            # Ensemble-averaged spectra output
            try:
                import matplotlib.pyplot as plt

                generated_csv = tmp_path / "averaged_spectra.csv"
                ax = ens_drv.plot_uv_vis_spectra(
                    results,
                    show_individual=True,
                    xlim_nm=(100, 200),
                    save_averaged_spectra=True,
                    averaged_spectra_filename=generated_csv,
                )

                assert ax is not None
                assert generated_csv.is_file()

                ref_csv = data_dir / "averaged_spectra.csv"
                assert ref_csv.is_file(), f"Missing reference CSV: {ref_csv}"

                got = np.loadtxt(generated_csv, delimiter=",", skiprows=1)
                ref = np.loadtxt(ref_csv, delimiter=",", skiprows=1)
                assert got.shape == ref.shape
                np.testing.assert_allclose(got, ref, rtol=0.0, atol=1.0e-6)

                plt.close(ax.figure)
            except ImportError:
                pass
