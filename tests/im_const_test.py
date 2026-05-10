from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
import os
import shutil

import numpy as np
import pytest

from veloxchem.imforcefieldgenerator import IMForceFieldGenerator
from veloxchem.interpolationdatapoint import InterpolationDatapoint
from veloxchem.interpolationdriver import InterpolationDriver
from veloxchem.molecule import Molecule
from veloxchem.xtbdriver import XtbDriver


def _h5py():
    return pytest.importorskip("h5py")


def _generated_db_name():
    return "im_database_0.h5"


def _reference_db_file():
    return Path(__file__).parent / "data" / "im_database_test.h5"


def _nve_energy_limits():
    return {
        "max_abs_drift_kj_mol": 1.0,
        "max_rel_drift": 1.0e-4,
    }


def _dataset_tolerances():
    return (
        ("_cartesian_coordinates", 1.0e-5, 0.0),
        ("_internal_coordinates", 1.0e-8, 0.0),
        ("_cartesian_gradient", 2.0e-6, 0.0),
        ("_gradient", 2.0e-6, 0.0),
        ("_cartesian_hessian", 2.0e-5, 0.0),
        ("_hessian", 2.0e-5, 0.0),
        ("_energy", 1.0e-5, 0.0),
        ("_confidence_radius", 1.0e1, 0.0),
        ("/xyz", 1.0e-9, 0.0)
    )
def get_xyz_structure():
    
    xyz_string = """9

C              2.166668000000         0.686818000000         0.513057000000
C              1.154761000000         0.429392000000        -0.592696000000
O              1.350483000000        -0.852819000000        -1.119526000000
H              3.198553000000         0.616151000000         0.108204000000
H              2.013132000000         1.703523000000         0.932376000000
H              2.040982000000        -0.060527000000         1.325068000000
H              1.286075000000         1.191508000000        -1.392806000000
H              0.127073000000         0.513996000000        -0.174440000000
H              0.667645000000        -0.964224000000        -1.831046000000
"""

    return xyz_string

@dataclass
class ConstructionContext:
    workdir: Path
    molecule: Molecule
    ffg: IMForceFieldGenerator
    results: dict
    db_file: Path


@contextmanager
def _temporary_cwd(path: Path):
    old_cwd = Path.cwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(old_cwd)


def _run_construction_once(workdir: Path) -> ConstructionContext:
    workdir.mkdir(parents=True, exist_ok=True)

    with _temporary_cwd(workdir):
        mol = Molecule.from_xyz_string(get_xyz_structure())
        mol.set_charge(0)
        mol.set_multiplicity(1)

        ffg = IMForceFieldGenerator(
            ground_state_driver=XtbDriver(),
            roots_to_follow=[0],
        )

        ffg.add_conformal_structures = False
        ffg.use_minimized_structures = [False, [], []]
        ffg.use_opt_confidence_radius = [True, "multi_grad", 0.5, 0.3]
        ffg.nsteps = 2000
        ffg.snapshots = 100
        # ffg.desired_point_density = 3
        ffg.energy_threshold = 0.009
        ffg.sampling_settings = {"enabled": False}
        ffg.open_mm_platform = "CPU"
        ffg.ensemble = "NVE"
        ffg.imforcefieldfiles = {0: _generated_db_name()}

        results = ffg.compute(mol)

    return ConstructionContext(
        workdir=workdir,
        molecule=mol,
        ffg=ffg,
        results=results,
        db_file=workdir / _generated_db_name(),
    )


def _cleanup_generated_files(workdir: Path):
    for pattern in (
        "im_database*.h5",
        "trajectory*.pdb",
        "traj*.pdb",
        "ref_structures.h5",
        "checkpoint*",
        "*.out",
        "random_structure.xyz",
    ):
        for path in workdir.glob(pattern):
            if path.is_dir():
                shutil.rmtree(path)
            elif path.exists():
                path.unlink()


@pytest.fixture()
def construction_context(tmp_path):
    workdir = tmp_path / "im_const_construction"
    try:
        yield _run_construction_once(workdir)
    finally:
        _cleanup_generated_files(workdir)


def _label_key(label):
    parts = label.split("_")
    point = int(parts[1]) if len(parts) > 1 and parts[0] == "point" else 10**9
    return point, label


def _normalize_z_matrix(z_matrix):
    return {
        key: [tuple(int(x) for x in row) for row in z_matrix.get(key, [])]
        for key in ("bonds", "angles", "dihedrals", "impropers")
    }


def _read_labels_and_z_matrix(db_file: Path, settings: dict):
    driver = InterpolationDriver()
    driver.update_settings(settings)
    driver.imforcefield_file = str(db_file)

    labels, z_matrix = driver.read_labels()
    return sorted(labels, key=_label_key), _normalize_z_matrix(z_matrix)


def _assert_database_is_runtime_readable(db_file: Path, settings: dict, molecule: Molecule):
    labels, z_matrix = _read_labels_and_z_matrix(db_file, settings)
    assert labels, "Generated interpolation database contains no datapoints"

    n_atoms = len(molecule.get_labels())
    n_cart = 3 * n_atoms

    for label in labels:
        dp = InterpolationDatapoint(z_matrix)
        dp.update_settings(settings)
        dp.read_hdf5(str(db_file), label)

        n_internal = len(dp.z_matrix)
        assert dp.cartesian_coordinates.shape == (n_atoms, 3)
        assert dp.gradient.shape == (n_atoms, 3)
        assert dp.hessian.shape == (n_cart, n_cart)
        assert dp.internal_coordinates_values.shape == (n_internal,)
        assert dp.internal_gradient.shape == (n_internal,)
        assert dp.internal_hessian.shape == (n_internal, n_internal)

        for arr in (
            dp.cartesian_coordinates,
            dp.gradient,
            dp.hessian,
            dp.internal_coordinates_values,
            dp.internal_gradient,
            dp.internal_hessian,
        ):
            assert np.all(np.isfinite(arr))


def _read_h5_datasets(db_file: Path):
    h5py = _h5py()
    datasets = {}

    with h5py.File(db_file, "r") as h5f:
        def visit(name, obj):
            if isinstance(obj, h5py.Dataset):
                datasets[name] = np.asarray(obj[()])

        h5f.visititems(visit)

    return datasets


def _tolerance_for_dataset(name):
    for suffix, atol, rtol in _dataset_tolerances():
        if name.endswith(suffix):
            return atol, rtol
    return 1.0e-12, 0.0

def _parse_xyz_bytes(xyz_entry):
    if isinstance(xyz_entry, bytes):
        text = xyz_entry.decode("utf-8")
    else:
        text = str(xyz_entry)

    lines = [line.strip() for line in text.splitlines() if line.strip()]

    n_atoms = int(lines[0])
    atom_lines = lines[1:1 + n_atoms]

    symbols = []
    coords = []

    for line in atom_lines:
        parts = line.split()
        symbols.append(parts[0])
        coords.append([float(parts[1]), float(parts[2]), float(parts[3])])

    return symbols, np.array(coords, dtype=float)


def _assert_dataset_matches(name, generated, reference):
    assert generated.shape == reference.shape

    if np.issubdtype(reference.dtype, np.floating):
        atol, rtol = _tolerance_for_dataset(name)
        diff = np.max(np.abs(generated - reference))
        assert np.allclose(generated, reference, atol=atol, rtol=rtol), (
            f"{name}: max abs diff {diff:.6e} exceeds atol={atol}, rtol={rtol}"
        )

    elif name.endswith("/xyz"):
        atol, rtol = _tolerance_for_dataset(name)

        for i, (g_xyz, r_xyz) in enumerate(zip(generated, reference)):
            g_symbols, g_coords = _parse_xyz_bytes(g_xyz)
            r_symbols, r_coords = _parse_xyz_bytes(r_xyz)

            assert g_symbols == r_symbols, (
                f"{name}[{i}]: atom labels differ\n"
                f"generated={g_symbols}\nreference={r_symbols}"
            )

            diff = np.max(np.abs(g_coords - r_coords))
            assert np.allclose(g_coords, r_coords, atol=atol, rtol=rtol), (
                f"{name}[{i}]: xyz coordinates differ; "
                f"max abs diff {diff:.6e} exceeds atol={atol}, rtol={rtol}"
            )

    else:
        assert np.array_equal(generated, reference), f"{name}: dataset differs"


def _assert_database_matches_reference(generated_db: Path, reference_db: Path, settings: dict):
    assert generated_db.exists()
    assert reference_db.exists()

    generated_labels, generated_z = _read_labels_and_z_matrix(generated_db, settings)
    reference_labels, reference_z = _read_labels_and_z_matrix(reference_db, settings)

    assert generated_labels == reference_labels
    assert generated_z == reference_z

    generated = _read_h5_datasets(generated_db)
    reference = _read_h5_datasets(reference_db)

    assert set(generated) == set(reference)

    for name in sorted(reference):
        _assert_dataset_matches(name, generated[name], reference[name])


def _flatten_float_series(values):
    flat = []

    def collect(value):
        if value is None:
            return
        if isinstance(value, (list, tuple)):
            for item in value:
                collect(item)
            return
        flat.extend(np.asarray(value, dtype=float).reshape(-1).tolist())

    collect(values)
    return np.asarray(flat, dtype=float)


def _assert_nve_total_energy_conserved(ffg: IMForceFieldGenerator):
    total = _flatten_float_series(ffg.total_energies)
    limits = _nve_energy_limits()

    assert total.size >= 2
    assert np.all(np.isfinite(total))

    drift = float(np.max(np.abs(total - total[0])))
    allowed = max(
        limits["max_abs_drift_kj_mol"],
        limits["max_rel_drift"] * max(abs(float(total[0])), 1.0),
    )

    assert drift <= allowed, (
        f"NVE total energy is not conserved: max drift={drift:.6f} kJ/mol, "
        f"allowed={allowed:.6f} kJ/mol"
    )


def test_generated_database_matches_reference_and_conserves_energy(construction_context):
    ctx = construction_context
    root = int(ctx.ffg.roots_to_follow[0])
    settings = dict(ctx.ffg.states_interpolation_settings[root])

    _assert_database_is_runtime_readable(ctx.db_file, settings, ctx.molecule)
    _assert_database_matches_reference(ctx.db_file, _reference_db_file(), settings)
    _assert_nve_total_energy_conserved(ctx.ffg)
