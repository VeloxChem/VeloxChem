from mpi4py import MPI
from importlib.util import find_spec
from pathlib import Path
import numpy as np
import pytest

from veloxchem import MofBuilder
from veloxchem import Molecule


def _copy_file(src, dest):
    if (not dest.is_file()) or (not src.samefile(dest)):
        if not dest.parent.is_dir():
            dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_text(src.read_text())


def _example_dir():
    return Path(__file__).resolve().parents[1] / "database" / "linker4test"


def _database_dir():
    return Path(__file__).resolve().parents[1] / "database"


def _stage_uio66_example(tmp_path):
    for name in ("bdc.xyz", "bdc.itp"):
        source = _example_dir() / name
        if not source.is_file():
            pytest.skip(f"UiO-66 example file is missing: {source}")
        _copy_file(source, tmp_path / name)


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
        _copy_file(source_db / "nodes_itps" / name, md_db / "nodes_itps" / name)
    for name in ("TIP3P.xyz", "TIP3P.itp"):
        _copy_file(source_db / "solvents_database" / name, md_db / "solvents_database" / name)
    for name in ("acetate.itp", "ooc.itp"):
        _copy_file(source_db / "terminations_itps" / name, md_db / "terminations_itps" / name)
    for name in ("em", "nvt", "npt"):
        (md_db / "mdps" / f"{name}.mdp").write_text("integrator = md\n", encoding="utf-8")

    return md_db


def _build_uio66_from_example(tmp_path):
    mol = Molecule.read_xyz_file("bdc.xyz")

    mof = MofBuilder(mof_family="UIO-66")
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


@pytest.mark.skipif(MPI.COMM_WORLD.Get_size() != 1,
                    reason='mofbuilder tests require a single MPI rank')
def test_uio66_example_workflow_build_remove_and_optional_openmm(tmp_path,
                                                                 monkeypatch):
    _stage_uio66_example(tmp_path)
    monkeypatch.chdir(tmp_path)

    uio = _build_uio66_from_example(tmp_path)
    _assert_built_uio66(uio)
    _assert_remove_returns_framework(uio)

    _write_and_solvate_like_example(uio, tmp_path)
    em_energy = _run_openmm_minimization_like_example(uio)

    reference_em_energy_kjmol = -803731.0
    em_energy_tolerance_kjmol = 500.0

    if find_spec("openmm") is not None:
        assert em_energy is not None
        assert em_energy == pytest.approx(
            reference_em_energy_kjmol,
            abs=em_energy_tolerance_kjmol,
        )
        assert Path("UiO66_MD_em_minimized.pdb").is_file()
