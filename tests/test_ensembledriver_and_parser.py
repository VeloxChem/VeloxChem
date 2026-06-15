from pathlib import Path
import pytest
import numpy as np

import veloxchem.ensembleparser as ensembleparser_module
from veloxchem.veloxchemlib import mpi_master
from veloxchem.ensembleparser import EnsembleParser
from veloxchem.ensembledriver import EnsembleDriver

pytest.importorskip("MDAnalysis")
pytest.importorskip("pyframe")


def make_simple_tip3_pdb(tmp_path):
    """Writes a simple PDB containing one QM atom and one TIP3 water."""

    def pdb_record(serial, atom_name, residue_name, resid, x, y, z, element):
        return (
            f'{"HETATM":<6s}{serial:5d} '
            f'{atom_name:^4s} '
            f'{residue_name:>4s}'
            f'A'
            f'{resid:4d}'
            f'    '
            f'{x:8.3f}{y:8.3f}{z:8.3f}'
            f'{1.00:6.2f}{0.00:6.2f}'
            f'          '
            f'{element:>2s}'
        )

    pdb_path = tmp_path / "simple_tip3.pdb"
    pdb_lines = [
        pdb_record(1, "C1", "LIG", 1, 0.0, 0.0, 0.0, "C"),
        pdb_record(2, "OH2", "TIP3", 2, 2.0, 0.0, 0.0, "O"),
        pdb_record(3, "H1", "TIP3", 2, 2.7, 0.0, 0.0, "H"),
        pdb_record(4, "H2", "TIP3", 2, 1.7, 0.7, 0.0, "H"),
        "END",
    ]
    pdb_path.write_text("\n".join(pdb_lines) + "\n")
    return pdb_path


class TestEnsemblePdbAndWaterAliases:

    @pytest.mark.parametrize(
        ("mda_version", "expected_kwargs"),
        [
            ("2.7.9", {"guess_bonds": True}),
            (
                "2.8.0",
                {"to_guess": ("bonds", "angles", "dihedrals")},
            ),
        ],
    )
    def test_mdanalysis_pdb_version_options(
        self, monkeypatch, mda_version, expected_kwargs
    ):
        class UniverseCalled(Exception):
            pass

        universe_kwargs = {}

        def fake_version(distribution_name):
            assert distribution_name == "MDAnalysis"
            return mda_version

        def fake_universe(*args, **kwargs):
            universe_kwargs.update(kwargs)
            raise UniverseCalled

        monkeypatch.setattr(ensembleparser_module, "version", fake_version)
        monkeypatch.setattr(
            ensembleparser_module.mda,
            "Universe",
            fake_universe,
        )

        ens_parser = EnsembleParser()
        ens_parser.ostream.mute()

        with pytest.raises(UniverseCalled):
            ens_parser.structures(
                pdb_file="structure.pdb",
                qm_region="resname LIG",
            )

        assert universe_kwargs == expected_kwargs

    def test_pdb_file_and_tip3_pe_pot_files(self, tmp_path):
        simple_tip3_pdb = make_simple_tip3_pdb(tmp_path)
        ens_parser = EnsembleParser()
        ens_parser.ostream.mute()

        snapshots = ens_parser.structures(
            pdb_file=simple_tip3_pdb,
            qm_region="resname LIG",
            pe_cutoff=3.0,
        )

        assert len(snapshots) == 1
        assert snapshots[0]["qm_atom_names"].tolist() == ["C1"]
        assert snapshots[0]["pe_atom_names"].tolist() == ["OH2", "H1", "H2"]
        assert snapshots[0]["pe_resnames"].tolist() == ["TIP3"] * 3

        tip3p_snapshot = dict(snapshots[0])
        tip3p_snapshot["frame"] = 1
        tip3p_snapshot["pe_resnames"] = np.array(["TIP3P"] * 3, dtype=object)

        ens_drv = EnsembleDriver()
        ens_drv.ostream.mute()
        ens_drv.set_env_models(pe_model="SEP")
        ens_drv.write_pot_files(
            [snapshots[0], tip3p_snapshot],
            outdir=tmp_path,
        )

        if ens_drv.rank == mpi_master():
            for frame, resname in ((0, "TIP3"), (1, "TIP3P")):
                pot_text = (
                    tmp_path / f"pe_frame_{frame:06d}.pot"
                ).read_text()
                charge_lines = pot_text.split(
                    "@charges\n", 1
                )[1].split("@end", 1)[0].splitlines()

                assert [line.split() for line in charge_lines] == [
                    ["O", "-0.67444000", f"{resname}_pe"],
                    ["H", "0.33722000", f"{resname}_pe"],
                    ["H", "0.33722000", f"{resname}_pe"],
                ]
                assert "@polarizabilities\n" in pot_text

    def test_tip3_npe_only_pot_files(self, tmp_path):
        simple_tip3_pdb = make_simple_tip3_pdb(tmp_path)
        ens_parser = EnsembleParser()
        ens_parser.ostream.mute()

        snapshots = ens_parser.structures(
            pdb_file=simple_tip3_pdb,
            qm_region="resname LIG",
            npe_cutoff=3.0,
        )

        tip3p_snapshot = dict(snapshots[0])
        tip3p_snapshot["frame"] = 1
        tip3p_snapshot["npe_resnames"] = np.array(["TIP3P"] * 3, dtype=object)

        ens_drv = EnsembleDriver()
        ens_drv.ostream.mute()
        ens_drv.set_env_models(pe_model="SEP", npe_model="tip3p")
        ens_drv.write_pot_files(
            [snapshots[0], tip3p_snapshot],
            outdir=tmp_path,
        )

        if ens_drv.rank == mpi_master():
            for frame, resname in ((0, "TIP3"), (1, "TIP3P")):
                pot_text = (
                    tmp_path / f"pe_frame_{frame:06d}.pot"
                ).read_text()
                charge_lines = pot_text.split(
                    "@charges\n", 1
                )[1].split("@end", 1)[0].splitlines()

                assert [line.split() for line in charge_lines] == [
                    ["O", "-0.83400000", f"{resname}_npe"],
                    ["H", "0.41700000", f"{resname}_npe"],
                    ["H", "0.41700000", f"{resname}_npe"],
                ]
                assert "@polarizabilities\n" not in pot_text


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
