from pathlib import Path
import pytest
import numpy as np

from veloxchem.veloxchemlib import mpi_master, hartree_in_ev
from veloxchem.ensembleparser import EnsembleParser
from veloxchem.ensembledriver import EnsembleDriver

pytest.importorskip("MDAnalysis")


@pytest.mark.solvers
class TestTrajectory:

    def test_compute_with_scf_and_property_options_for_bithio(self):
        """
        Regression test for notebook-style ensemble API:
        - compute(..., scf_options=..., property_options=...)
        - property='absorption' + nstates routes to TD-DFT/LR path
        - the fika potfile is frame-specific for each snapshot
        - compare LR eigenvalues and oscillator strengths to fika reference output
        - check that the generated PE+NPE .pot files exist
        """
        data_dir = Path(__file__).parent / "data"
        traj = str(data_dir / "bithio.xtc")
        top = str(data_dir / "bithio.tpr")

        ens_parser = EnsembleParser()
        ens_parser.ostream.mute()

        ensemble = ens_parser.structures(
            trajectory_file=traj,
            topology_file=top,
            qm_region="resname RES",
            start=300,
            end=400,
            num_snapshots=2,
            pe_cutoff=3.0,
            npe_cutoff=4.5,
        )

        ens_drv = EnsembleDriver()
        ens_drv.ostream.mute()

        ens_drv.set_env_models(
            pe_model=["CP3", "SEP"],
            npe_model=["ff19sb", "tip3p"],
        )

        scf_options = {
            "scf_type": "restricted",
            "max_iter": 150,
        }

        property_options = {
            "property": "absorption",
            "nstates": 2,
        }

        potdir = data_dir
        results = ens_drv.compute(
            ensemble,
            basis_set="sto-3g",
            scf_options=scf_options,
            property_options=property_options,
            potdir=potdir,
        )

        assert "scf_all" in results
        assert "rsp_all" in results
        assert len(results["scf_all"]) == 2
        assert len(results["rsp_all"]) == 2

        # Reference values from validated fika output
        ref_frames = [6, 8]
        ref_eigs = {
            6: np.array([5.78875447, 7.17288998]) / hartree_in_ev(),
            8: np.array([5.62513120, 7.19323719]) / hartree_in_ev(),
        }
        ref_osc = {
            6: np.array([0.54164731, 0.01810991]),
            8: np.array([0.55327598, 0.12162315]),
        }

        got_frames = [int(frame) for frame, _ in results["scf_all"]]
        assert got_frames == ref_frames

        if ens_drv.rank == mpi_master():
            # got_scf = np.array([float(scf_res["scf_energy"]) for _, scf_res in results["scf_all"]])

            for frame, scf_res in results["scf_all"]:
                expected_name = f"pe_frame_{int(frame):06d}.pot"
                expected_path = potdir / expected_name

                # Each snapshot must use its own potfile
                assert expected_path.is_file()

            for frame, rsp_res in results["rsp_all"]:
                frame = int(frame)
                eig = np.array(rsp_res.get("eigenvalues", []), dtype=float)
                osc = np.array(rsp_res.get("oscillator_strengths", []),
                               dtype=float)
                assert int(rsp_res.get("number_of_states", 0)) == 2
                assert eig.shape == (2,)
                assert osc.shape == (2,)

                assert np.max(
                    np.abs(np.array(eig) -
                           np.array(ref_eigs[frame]))) < 0.02 / hartree_in_ev()
                assert np.max(np.abs(np.array(osc) -
                                     np.array(ref_osc[frame]))) < 0.01

            # cleanup
            (data_dir / ".bithio.xtc_offsets.npz").unlink(missing_ok=True)
            (data_dir / "pe_frame_000006.pot").unlink(missing_ok=True)
            (data_dir / "pe_frame_000008.pot").unlink(missing_ok=True)
