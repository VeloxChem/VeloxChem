from pathlib import Path
import numpy as np
import pytest

from veloxchem.mpitask import MpiTask
from veloxchem.ensembledriver import EnsembleDriver
from veloxchem.veloxchemlib import mpi_master, hartree_in_ev

pytest.importorskip("MDAnalysis")


@pytest.mark.solvers
class TestEnsembleDriverInput:

    def test_trajectory_bithio(self):

        here = Path(__file__).parent
        data_dir = here / "data"
        inpfile = data_dir / "bithio.inp"

        task = MpiTask([str(inpfile), None])

        ensemble_dict = dict(task.input_dict["ensemble"])
        ensemble_dict["trajectory_file"] = str(data_dir / "bithio.xtc")
        ensemble_dict["topology_file"] = str(data_dir / "bithio.tpr")
        ensemble_dict["output_dir"] = str(data_dir)

        ens_drv = EnsembleDriver(task.mpi_comm, task.ostream)
        ens_drv.update_settings(
            ensemble_dict,
            task.input_dict.get("spectrum_settings"),
            task.input_dict.get("scf"),
            task.input_dict.get("response"),
            task.input_dict.get("method_settings"),
        )

        try:
            results = ens_drv.compute_from_input()

            if task.mpi_rank == mpi_master():
                ref_frames = [6, 8]
                ref_exc_energies = np.array(
                    [5.790027, 7.172051, 5.628371, 7.199163])
                ref_osc_strengths = np.array(
                    [0.5414, 0.0181, 0.5529, 0.1210])

                got_frames = [
                    int(frame) for frame, _ in results["rsp_all"]
                ]
                assert got_frames == ref_frames

                exc_energies = []
                osc_strengths = []
                for _, rsp_res in results["rsp_all"]:
                    exc_energies.extend(
                        np.asarray(rsp_res["eigenvalues"], dtype=float) *
                        hartree_in_ev())
                    osc_strengths.extend(
                        np.asarray(rsp_res["oscillator_strengths"],
                                   dtype=float))

                assert np.max(
                    np.abs(np.array(exc_energies) -
                           ref_exc_energies)) < 0.02
                assert np.max(
                    np.abs(np.array(osc_strengths) -
                           ref_osc_strengths)) < 0.01

        finally:
            if task.mpi_rank == mpi_master():
                (data_dir / ".bithio.xtc_offsets.npz").unlink(
                    missing_ok=True)
                for frame in (6, 8):
                    (data_dir / f"pe_frame_{frame:06d}.pot").unlink(
                        missing_ok=True)
