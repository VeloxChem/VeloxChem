from importlib.util import find_spec
from pathlib import Path
import shutil

import numpy as np
import pytest
from mpi4py import MPI

import veloxchem as vlx


pytestmark = [
    pytest.mark.mofbuilder,
]


EXAMPLE_FILES = ("bdc.xyz", "bdc.itp")
GENERATED_NAMES = (
    "uio66_original.gro",
    "UIO-66_in_solvent.gro",
    "UiO66_MD_em.log",
    "UiO66_MD_em_minimized.pdb",
    "UiO66_MD_whole.pdb",
    "Linker.itp",
    "bdc1.xyz",
    "MOL_scf.h5",
)
GENERATED_DIRS = ("MD_run", "md_database", "nodes", "output")
REFERENCE_EM_ENERGY_KJMOL = -803731.0
EM_ENERGY_TOLERANCE_KJMOL = 500.0
COMM = MPI.COMM_WORLD


def _example_dir():
    return Path(__file__).resolve().parents[1] / "database" / "linker4test"


def _database_dir():
    return Path(__file__).resolve().parents[1] / "database"


def _skip_without_merged_mofbuilder():
    try:
        builder = vlx.MofBuilder(comm=MPI.COMM_SELF)
    except TypeError:
        builder = vlx.MofBuilder()

    if not hasattr(builder, "ostream"):
        pytest.skip("MOFBuilder integration test requires the merged VeloxChem API.")


def _wait_for_master_rank_then_skip():
    if COMM.Get_rank() != 0:
        COMM.Barrier()
        pytest.skip("MOFBuilder integration workflow runs on the MPI master rank.")


def _stage_uio66_example(tmp_path):
    for name in EXAMPLE_FILES:
        source = _example_dir() / name
        if not source.is_file():
            pytest.skip(f"UiO-66 example file is missing: {source}")
        shutil.copyfile(source, tmp_path / name)


def _stage_md_database(tmp_path):
    """Create the small database subset needed after UiO-66 is built."""

    source_db = _database_dir()
    md_db = tmp_path / "md_database"
    (md_db / "nodes_itps").mkdir(parents=True)
    (md_db / "terminations_itps").mkdir()
    (md_db / "solvents_database").mkdir()
    (md_db / "mdps").mkdir()
    (md_db / "amber14sb_OL21.ff").mkdir()

    for name in ("template.top", "Zr.itp", "O.itp", "HO.itp", "HHO.itp"):
        shutil.copyfile(source_db / "nodes_itps" / name,
                        md_db / "nodes_itps" / name)
    for name in ("TIP3P.xyz", "TIP3P.itp"):
        shutil.copyfile(source_db / "solvents_database" / name,
                        md_db / "solvents_database" / name)
    for name in ("acetate.itp", "ooc.itp"):
        shutil.copyfile(source_db / "terminations_itps" / name,
                        md_db / "terminations_itps" / name)
    for name in ("em", "nvt", "npt"):
        (md_db / "mdps" / f"{name}.mdp").write_text(
            "integrator = md\n",
            encoding="utf-8",
        )

    return md_db


def _build_uio66_from_example(tmp_path):
    mol = vlx.Molecule.read_xyz_file("bdc.xyz")

    mof = vlx.MofBuilder(comm=MPI.COMM_SELF, mof_family="UIO-66")
    mof.ostream.mute()
    mof.data_path = str(_database_dir())
    mof.target_directory = str(tmp_path)
    mof.linker_molecule = mol
    mof.node_metal = "Zr"

    return mof.build()


def _first_residue_index(framework):
    return next(iter(framework.graph_index_name_dict))


def _assert_built_uio66(framework):
    assert framework.graph.number_of_nodes() > 0
    assert framework.graph.number_of_edges() > 0
    assert framework.framework_data.shape[1] == 11
    assert framework.framework_fcoords_data.shape[1] == 11
    assert np.isfinite(framework.framework_data[:, 5:8].astype(float)).all()
    assert "EDGE" in framework.residues_info


def _assert_remove_returns_framework(framework):
    removed = framework.remove(remove_indices=[_first_residue_index(framework)])
    assert removed is not framework
    assert removed.graph is not framework.graph
    assert removed.framework_data.shape[1] == 11
    assert np.isfinite(removed.framework_data[:, 5:8].astype(float)).all()


def _write_and_solvate_like_example(framework, tmp_path):
    framework.data_path = str(_stage_md_database(tmp_path))
    framework.write(format=["gro"], filename="uio66_original")
    framework.solvate(solvents_quantities=[0], padding_angstrom=50)

    assert Path("uio66_original.gro").is_file()
    assert Path(framework.solvated_gro_file).is_file()
    assert framework.solvents_dict is None


def _run_openmm_minimization_like_example(framework):
    if find_spec("openmm") is None:
        return None
    from openmm import unit

    framework.provided_linker_itpfile = "bdc.itp"
    framework.md_prepare()
    framework.md_driver.run_pipeline(
        steps=["em"],
        output_prefix="UiO66_MD",
    )

    em_state = framework.md_driver.simulation.context.getState(getEnergy=True)
    return em_state.getPotentialEnergy().value_in_unit(
        unit.kilojoule_per_mole)


def _unlink_generated_files(tmp_path):
    for name in GENERATED_NAMES:
        path = tmp_path / name
        if path.exists():
            path.unlink()

    for dirname in GENERATED_DIRS:
        path = tmp_path / dirname
        if path.exists():
            shutil.rmtree(path)


def test_uio66_example_workflow_build_remove_and_optional_openmm(tmp_path,
                                                                 monkeypatch):
    _skip_without_merged_mofbuilder()
    _wait_for_master_rank_then_skip()
    _stage_uio66_example(tmp_path)
    monkeypatch.chdir(tmp_path)

    try:
        uio = _build_uio66_from_example(tmp_path)
        _assert_built_uio66(uio)
        _assert_remove_returns_framework(uio)

        _write_and_solvate_like_example(uio, tmp_path)
        em_energy = _run_openmm_minimization_like_example(uio)

        if find_spec("openmm") is not None:
            assert em_energy is not None
            assert em_energy == pytest.approx(
                REFERENCE_EM_ENERGY_KJMOL,
                abs=EM_ENERGY_TOLERANCE_KJMOL,
            )
            assert Path("UiO66_MD_em_minimized.pdb").is_file()
    finally:
        _unlink_generated_files(tmp_path)
        COMM.Barrier()

    for name in GENERATED_NAMES:
        assert not (tmp_path / name).exists()
    for dirname in GENERATED_DIRS:
        assert not (tmp_path / dirname).exists()
