from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pytest

from veloxchem.imforcefieldgenerator import IMForceFieldGenerator
from veloxchem.interpolationdatapoint import InterpolationDatapoint
from veloxchem.interpolationdriver import InterpolationDriver
from veloxchem.molecule import Molecule
from veloxchem.xtbdriver import XtbDriver


MOLECULE_XYZ = """9

C             -0.655611975218         1.123063867023         1.299900212542
C              0.166450587214        -0.077365693786         0.842810693258
O              0.903745662597        -0.664310974466         1.889735012507
H             -1.182307655193         1.557316856177         0.454684492417
H             -0.006368262160         1.877145043389         1.738051941292
H             -1.389155036746         0.820022591654         2.046171562256
H              0.909654941891         0.230639787827         0.104018658833
H             -0.494545585471        -0.824757720464         0.378163799018
H              0.296192323085        -0.898429757353         2.600555627877"""

# JSON-safe representation of the expected symmetry dictionary.
EXPECTED_SYMMETRY = {
    "gs": [
        [0, 1, 2, 3, 4, 5, 6, 7, 8],
        [[3, 4, 5]],
        [[3, 4, 5]],
        [0, 1, 2, 6, 7, 8],
        [3, 4, 5],
        [[1, 2], [0, 1], [1, 2], [0, 1]],
        [],
        {
            "2": {},
            "3": {
                "(0,1)": [
                    [0, 1, 2, 3], [0, 1, 3, 6], [0, 1, 3, 7],
                    [0, 1, 2, 4], [0, 1, 4, 6], [0, 1, 4, 7],
                    [0, 1, 2, 5], [0, 1, 5, 6], [0, 1, 5, 7],
                ]
            },
        },
        {"(3,0,1,2)": []},
        [21, 33],
    ],
    "es": [],
}

REFERENCE_Z_MATRIX = {0: 
                      {'bonds': [(0, np.int32(1)), (0, np.int32(3)), (0, np.int32(4)), (0, np.int32(5)), (1, np.int32(2)), (1, np.int32(6)), (1, np.int32(7)), (2, np.int32(8))], 
                       'angles': [(np.int32(1), 0, np.int32(3)), (np.int32(1), 0, np.int32(4)), (np.int32(1), 0, np.int32(5)), (np.int32(3), 0, np.int32(4)), (np.int32(3), 0, np.int32(5)), (np.int32(4), 0, np.int32(5)), (0, 1, np.int32(2)), (0, 1, np.int32(6)), (0, 1, np.int32(7)), (np.int32(2), 1, np.int32(6)), (np.int32(2), 1, np.int32(7)), (np.int32(6), 1, np.int32(7)), (1, 2, np.int32(8))], 
                       'dihedrals': [(3, 0, 1, 2), (3, 0, 1, 6), (3, 0, 1, 7), (4, 0, 1, 2), (4, 0, 1, 6), (4, 0, 1, 7), (5, 0, 1, 2), (5, 0, 1, 6), (5, 0, 1, 7), (0, 1, 2, 8), (6, 1, 2, 8), (7, 1, 2, 8)], 
                       'impropers': []}}

## Symmetry optimization molecules comparison:

OPT_MOL_1_STRING = """9

C             -0.655580423664         1.123115508625         1.299905581133
C              0.166449320140        -0.077336412521         0.842837685759
O              0.903752371904        -0.664301479281         1.889675759041
H             -1.182201471608         1.557229186477         0.454606835561
H             -0.006418743400         1.877130716881         1.738239109422
H             -1.389233978409         0.819966682549         2.046038237729
H              0.909641419994         0.230551555451         0.103993901290
H             -0.494693912107        -0.824650830610         0.378245627787
H              0.296340417149        -0.898380927570         2.600549262280"""

OPT_MOL_2_STRING = """9

C             -0.647714404113         1.134335134189         1.296447845196
C              0.171149770792        -0.077303548522         0.845107632414
O              0.351921947330        -1.035191705507         1.861936415014
H             -1.177092780815         1.570514407703         0.453849379992
H             -0.021719671233         1.886550072965         1.765959893521
H             -1.397693947831         0.818363172250         2.022160360592
H              1.176045793982         0.221543381083         0.542361736159
H             -0.320660682769        -0.549882764723        -0.019136216845
H             -0.513099199773        -1.308623730321         2.186897137413"""

OPT_MOL_3_STRING = """9

C             -0.655348659857         1.138863627644         1.296867545148
C              0.163683350983        -0.081642488805         0.846682265696
O             -0.185284729664        -1.264020601080         1.529433631083
H             -1.193474188122         1.583501956968         0.463848757939
H             -0.016774546403         1.900156067723         1.734727117012
H             -1.376944678166         0.836832048032         2.053200501291
H              1.223601596492         0.061871721897         1.064863411744
H              0.058339669271        -0.231414660394        -0.237110752236
H             -1.114939767636        -1.455834761968         1.362959680891"""

OPT_MOL_4_STRING = """9

C             -0.661535066309         1.128097375722         1.298956571720
C              0.160830459103        -0.078956848612         0.840928833654
O             -0.545706317697        -1.293534712398         0.947688936392
H             -1.201152069039         1.577440238953         0.467697680596
H             -0.006154222314         1.890329398828         1.712599756591
H             -1.363731740381         0.823716932687         2.070509517949
H              1.028694858953        -0.210814780197         1.491718794709
H              0.521388572669         0.070466580703        -0.186332296889
H             -1.346953592312        -1.233309727020         0.415095469234"""

OPT_MOL_5_STRING = """9

C             -0.655888545091         1.123051448051         1.299969381073
C              0.162472465612        -0.080425730887         0.844024013115
O              0.905947325261        -0.661846170418         1.889548753643
H             -1.149879446832         1.581875670301         0.448103636019
H             -0.006681704809         1.858571145779         1.768976213881
H             -1.414536349453         0.819082902756         2.020093107427
H              0.902422189935         0.224651088670         0.100602256785
H             -0.500464533914        -0.829493363414         0.385119605825
H              0.302176323885        -0.895963610881         2.603517471811"""

OPT_MOL_6_STRING = """9

C             -0.647225050482         1.132076466130         1.297152546775
C              0.164215418467        -0.083592440538         0.844693396176
O              0.352454572640        -1.033974440902         1.865462253270
H             -1.146202760626         1.592301946343         0.448793396862
H             -0.021138190107         1.864250963250         1.797182615542
H             -1.421862178965         0.815114345749         1.996220384765
H              1.169136618722         0.216573475114         0.540096792373
H             -0.323338911811        -0.549622907571        -0.025195796940
H             -0.509630434066        -1.297842941244         2.206196307525"""

OPT_MOL_7_STRING = """9

C             -0.654790771247         1.132174314204         1.295115494626
C              0.165823468566        -0.087277110279         0.843099413861
O             -0.180949268889        -1.270007827414         1.524383391676
H             -1.155566471514         1.603991938232         0.453881625923
H             -0.025677957276         1.873848499797         1.778107002766
H             -1.410556749706         0.821300531136         2.013811163223
H              1.226629809902         0.056751004794         1.055724725311
H              0.053783537038        -0.225731921612        -0.241950168237
H             -1.125708426544        -1.425778262201         1.414715151998"""

OPT_MOL_8_STRING = """9

C             -0.660950650731         1.125442847606         1.295944981818
C              0.164099312044        -0.083747383798         0.844177409183
O             -0.543708045997        -1.297765166689         0.952833718950
H             -1.160014119028         1.600477862412         0.454173569747
H             -0.014535205542         1.870142442125         1.753480125147
H             -1.397681673782         0.815495534403         2.032633532226
H              1.037060158009        -0.213768528507         1.486914081486
H              0.508276143622         0.065364154984        -0.189161661915
H             -1.363403616834        -1.223015660096         0.451209415865
"""


MOL_OPT_ENERGIES = {1: np.array(-11.394330242869191), 2: np.array(-11.39223432073637), 3: np.array(-11.390385942605404), 4: np.array(-11.392612073904237), 5: np.array(-11.394337798009207), 6: np.array(-11.392480708695425), 7: np.array(-11.390327117819387), 8: np.array(-11.392293372650661)}



# Tolerances separated by intent.
# - roundtrip checks are strict (read/write should be bitwise-close),
# - analytical rebuild checks allow small algebraic/SVD noise.
TOL = {
    "event_energy_abs": 1.0e-12,
    "internal_coordinates_abs": 1.0e-8,
    "internal_gradient_abs": 2.0e-6,
    "internal_hessian_abs": 2.0e-5,
    "backtransform_gradient_abs": 2.0e-6,
    "backtransform_hessian_abs": 2.0e-5,
    "roundtrip_numeric_abs": 1.0e-12,
}


def _normalize_key(key):
    if hasattr(key, "tolist"):
        key = key.tolist()
    if isinstance(key, (tuple, list)):
        return ",".join(str(int(x)) for x in key)
    if isinstance(key, (np.integer,)):
        return str(int(key))
    return str(key).replace(" ", "")


def canonicalize(obj):
    """Convert runtime symmetry payloads into deterministic JSON-safe form."""
    if isinstance(obj, dict):
        items = [(_normalize_key(k), canonicalize(v)) for k, v in obj.items()]
        items.sort(key=lambda kv: kv[0])
        return {k: v for k, v in items}
    if isinstance(obj, (list, tuple)):
        return [canonicalize(v) for v in obj]
    if hasattr(obj, "tolist"):
        return canonicalize(obj.tolist())
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating,)):
        return float(obj)
    return obj


@dataclass
class HookLog:
    events: list

    def by_event(self, name):
        return [e for e in self.events if e["event"] == name]


@dataclass
class ConstructionContext:
    workdir: Path
    molecule: Molecule
    ffg: IMForceFieldGenerator
    results: dict
    hooks: HookLog
    db_file: Path


def _make_hooks(log: HookLog):
    """Collect every emitted event from the generator in one list."""

    def wildcard_cb(envelope):
        log.events.append(envelope)

    return {
        "__strict__": True,
        "__run_id__": "full-single-file",
        "*": wildcard_cb,
    }


def _mode_suffix_from_settings(settings):
    """Return the dataset suffix used by InterpolationDatapoint write/read."""
    if settings.get("use_inverse_bond_length", False):
        bond_tag = "_rinv"
    elif settings.get("use_eq_bond_length", False):
        bond_tag = "_eq"
    else:
        bond_tag = "_r"

    dihedral_tag = "_cosine" if settings.get("use_cosine_dihedral", False) else "_dihedral"
    return f"{bond_tag}{dihedral_tag}"


def _run_construction_once(workdir: Path) -> ConstructionContext:
    """
    Executes exactly one full construction run in an isolated temp directory.

    The resulting context is reused by all assertions so the integration test
    exercises a single realistic run end-to-end.
    """
    import os

    os.chdir(workdir)

    mol = Molecule.from_xyz_string(MOLECULE_XYZ)
    mol.set_charge(0)
    mol.set_multiplicity(1)
    
    gs_driver = XtbDriver()
    ffg = IMForceFieldGenerator(ground_state_driver=gs_driver, roots_to_follow=[0])

    # deterministic/smoke construction settings
    ffg.use_symmetry = True
    ffg.add_conformal_structures = False
    ffg.use_minimized_structures = [False, [], []]
    ffg.use_opt_confidence_radius = [False, "multi_grad", 0.5, 0.5]
    ffg.nsteps = 2000
    ffg.snapshots = 200
    ffg.desired_point_density = 2
    ffg.energy_threshold = 0.00001
    ffg.sampling_settings = {"enabled": False}
    ffg.open_mm_platform = "CPU"
    ffg.ensemble = "NVE"

    db_file = workdir / "im_database_0.h5"
    ffg.imforcefieldfiles = {0: str(db_file)}
    ffg.trajectory_file = str(workdir / "traj_construction.pdb")
    ffg.reference_struc_energy_file = str(workdir / "ref_structures.h5")

    log = HookLog(events=[])
    results = ffg.compute(mol, test_hooks=_make_hooks(log))

    return ConstructionContext(
        workdir=workdir,
        molecule=mol,
        ffg=ffg,
        results=results,
        hooks=log,
        db_file=db_file,
    )


def _assert_setup_symmetry(ctx: ConstructionContext):
    """Checks that runtime setup symmetry payload equals hardcoded reference."""
    ev = ctx.hooks.by_event("ffg.setup_complete")
    assert len(ev) == 1

    runtime_sym = canonicalize(ev[0]["payload"]["symmetry_information"])
    assert runtime_sym == canonicalize(EXPECTED_SYMMETRY)

    runtime_z_matrix = canonicalize(ev[0]["payload"]["z_matrix_summary"])

    assert runtime_z_matrix == canonicalize(REFERENCE_Z_MATRIX)


def _collect_datapoint_events(ctx: ConstructionContext):
    """
    Flatten datapoint-written events from both initial add-point and qmmm phases.

    Returns:
        dict[(root, label)] -> energy_hartree
    """
    observed = {}

    for envelope in ctx.hooks.events:
        event = envelope["event"]
        payload = envelope["payload"]

        if event == "datapoint_written":
            root = payload.get("root")
            label = payload.get("label")
            energy = payload.get("energy_hartree")
        elif event == "qmmm.datapoint_written":
            collector = payload.get("collector", {})
            root = collector.get("root")
            label = collector.get("label")
            energy = collector.get("energy_hartree")
        else:
            continue

        if root is None or label is None or energy is None:
            continue

        observed[(int(root), str(label))] = float(energy)

    return observed


def _read_all_datapoints_from_db(ctx: ConstructionContext):
    """
    Reads labels and datapoints through the production InterpolationDriver path.

    This verifies the same API used in normal runtime restart/loading behavior.
    """
    root = int(ctx.ffg.roots_to_follow[0])
    settings = dict(ctx.ffg.states_interpolation_settings[root])
    z_matrix = ctx.ffg.roots_z_matrix[root]

    driver = InterpolationDriver(z_matrix)
    driver.update_settings(settings)
    driver.imforcefield_file = str(ctx.db_file)

    labels, z_matrix_from_file = driver.read_labels()

    datapoints = []
    for label in labels:
        dp = InterpolationDatapoint(z_matrix_from_file)
        dp.update_settings(settings)
        dp.read_hdf5(str(ctx.db_file), label)
        datapoints.append((label, dp))

    return root, settings, z_matrix_from_file, labels, datapoints


def _assert_cartesian_tensor_shapes(dp: InterpolationDatapoint, n_atoms: int):
    """Basic schema sanity for all mandatory arrays on one datapoint."""
    assert dp.cartesian_coordinates.shape == (n_atoms, 3)
    assert dp.gradient.shape == (n_atoms, 3)
    assert dp.hessian.shape == (3 * n_atoms, 3 * n_atoms)

    n_internal = len(dp.z_matrix)
    assert dp.internal_coordinates_values.shape == (n_internal,)
    assert dp.internal_gradient.shape == (n_internal,)
    assert dp.internal_hessian.shape == (n_internal, n_internal)

    assert np.all(np.isfinite(dp.cartesian_coordinates))
    assert np.all(np.isfinite(dp.gradient))
    assert np.all(np.isfinite(dp.hessian))
    assert np.all(np.isfinite(dp.internal_coordinates_values))
    assert np.all(np.isfinite(dp.internal_gradient))
    assert np.all(np.isfinite(dp.internal_hessian))


def _rebuild_from_cartesian_and_compare(
    dp: InterpolationDatapoint,
    z_matrix,
    settings,
    inv_sqrt_masses,
):
    """
    Recompute internal arrays from stored Cartesian arrays and compare.

    This is the analytical core check:
    - internal coordinates from Cartesian geometry,
    - Cartesian->internal gradient/Hessian transform,
    - internal->Cartesian backtransform consistency.
    """
    rebuilt = InterpolationDatapoint(z_matrix)
    rebuilt.update_settings(settings)

    rebuilt.cartesian_coordinates = np.array(dp.cartesian_coordinates, dtype=np.float64, copy=True)
    rebuilt.eq_bond_lengths = None if dp.eq_bond_lengths is None else np.array(dp.eq_bond_lengths, dtype=np.float64, copy=True)
    rebuilt.mapping_masks = None if dp.mapping_masks is None else np.array(dp.mapping_masks, copy=True)
    rebuilt.inv_sqrt_masses = np.array(inv_sqrt_masses, dtype=np.float64, copy=True)

    rebuilt.energy = float(np.asarray(dp.energy).reshape(-1)[0])
    rebuilt.gradient = np.array(dp.gradient, dtype=np.float64, copy=True)
    rebuilt.hessian = np.array(dp.hessian, dtype=np.float64, copy=True)
    rebuilt.confidence_radius = float(np.asarray(dp.confidence_radius).reshape(-1)[0])

    rebuilt.transform_gradient_and_hessian()

    assert np.allclose(
        rebuilt.internal_coordinates_values,
        dp.internal_coordinates_values,
        atol=TOL["internal_coordinates_abs"],
        rtol=0.0,
    )
    assert np.allclose(
        rebuilt.internal_gradient,
        dp.internal_gradient,
        atol=TOL["internal_gradient_abs"],
        rtol=0.0,
    )
    assert np.allclose(
        rebuilt.internal_hessian,
        dp.internal_hessian,
        atol=TOL["internal_hessian_abs"],
        rtol=0.0,
    )

    grad_back = rebuilt.backtransform_internal_gradient_to_cartesian_coordinates()
    hess_back = rebuilt.backtransform_internal_hessian_to_cartesian_coordinates()

    assert np.allclose(
        grad_back,
        dp.gradient,
        atol=TOL["backtransform_gradient_abs"],
        rtol=0.0,
    )
    assert np.allclose(
        hess_back,
        dp.hessian,
        atol=TOL["backtransform_hessian_abs"],
        rtol=0.0,
    )


def _assert_roundtrip_write_read(ctx: ConstructionContext, z_matrix, settings, labels, datapoints):
    """
    Write every datapoint into a fresh DB and verify strict readback equality.

    This validates database I/O integrity independently of the original file.
    """
    roundtrip_file = ctx.workdir / "im_database_roundtrip_check.h5"

    for label, dp in datapoints:
        dp.write_hdf5(str(roundtrip_file), label)

    for label, dp_ref in datapoints:
        dp_new = InterpolationDatapoint(z_matrix)
        dp_new.update_settings(settings)
        dp_new.read_hdf5(str(roundtrip_file), label)

        assert float(dp_new.energy) == pytest.approx(float(dp_ref.energy), abs=TOL["roundtrip_numeric_abs"])
        assert np.allclose(dp_new.cartesian_coordinates, dp_ref.cartesian_coordinates, atol=TOL["roundtrip_numeric_abs"], rtol=0.0)
        assert np.allclose(dp_new.gradient, dp_ref.gradient, atol=TOL["roundtrip_numeric_abs"], rtol=0.0)
        assert np.allclose(dp_new.hessian, dp_ref.hessian, atol=TOL["roundtrip_numeric_abs"], rtol=0.0)
        assert np.allclose(dp_new.internal_coordinates_values, dp_ref.internal_coordinates_values, atol=TOL["roundtrip_numeric_abs"], rtol=0.0)
        assert np.allclose(dp_new.internal_gradient, dp_ref.internal_gradient, atol=TOL["roundtrip_numeric_abs"], rtol=0.0)
        assert np.allclose(dp_new.internal_hessian, dp_ref.internal_hessian, atol=TOL["roundtrip_numeric_abs"], rtol=0.0)

        if dp_ref.eq_bond_lengths is None:
            assert dp_new.eq_bond_lengths is None
        else:
            assert np.allclose(dp_new.eq_bond_lengths, dp_ref.eq_bond_lengths, atol=TOL["roundtrip_numeric_abs"], rtol=0.0)

        if dp_ref.mapping_masks is None:
            assert dp_new.mapping_masks is None
        else:
            assert np.array_equal(dp_new.mapping_masks, dp_ref.mapping_masks)

        assert dp_new.family_label == dp_ref.family_label
        assert dp_new.bank_role == dp_ref.bank_role
        assert dp_new.cluster_type == dp_ref.cluster_type
        assert dp_new.cluster_id == dp_ref.cluster_id
        assert dp_new.cluster_rotor_ids == dp_ref.cluster_rotor_ids
        assert dp_new.cluster_state_id == dp_ref.cluster_state_id
        assert dp_new.dihedrals_to_rotate == dp_ref.dihedrals_to_rotate

        if dp_ref.phase_signature is None:
            assert dp_new.phase_signature is None
        else:
            assert np.allclose(dp_new.phase_signature, dp_ref.phase_signature, atol=TOL["roundtrip_numeric_abs"], rtol=0.0)


def _assert_database_datapoint_integrity(ctx: ConstructionContext):
    root, settings, z_matrix, labels, datapoints = _read_all_datapoints_from_db(ctx)
    assert len(labels) > 0, "Database is empty; expected at least one datapoint"

    mode_suffix = _mode_suffix_from_settings(settings)
    observed = _collect_datapoint_events(ctx)

    n_atoms = len(ctx.molecule.get_labels())
    inv_sqrt_masses = 1.0 / np.sqrt(np.repeat(ctx.molecule.get_masses(), 3))

    # Label-based mapping avoids fragile counter/index assumptions.
    expected_by_label = {
        "point_1": (OPT_MOL_1_STRING, MOL_OPT_ENERGIES[1]),
        "point_1_symmetry_1": (OPT_MOL_2_STRING, MOL_OPT_ENERGIES[2]),
        "point_1_symmetry_2": (OPT_MOL_3_STRING, MOL_OPT_ENERGIES[3]),
        "point_1_symmetry_3": (OPT_MOL_4_STRING, MOL_OPT_ENERGIES[4]),
        "point_2": (OPT_MOL_5_STRING, MOL_OPT_ENERGIES[5]),
        "point_2_symmetry_1": (OPT_MOL_6_STRING, MOL_OPT_ENERGIES[6]),
        "point_2_symmetry_2": (OPT_MOL_7_STRING, MOL_OPT_ENERGIES[7]),
        "point_2_symmetry_3": (OPT_MOL_8_STRING, MOL_OPT_ENERGIES[8]),
    }

    for label, dp in datapoints:
        key = (root, label)
        assert key in observed, f"Missing datapoint_written event for root={root}, label={label}"
        assert float(dp.energy) == pytest.approx(observed[key], abs=TOL["event_energy_abs"])
        assert str(dp.point_label).endswith(mode_suffix)

        assert label in expected_by_label, f"Unexpected datapoint label in DB: {label}"
        xyz_ref, e_ref = expected_by_label[label]
        cur_molecule = Molecule.from_xyz_string(xyz_ref)

        # Use element-wise comparison; mean-difference can hide errors.
        assert np.allclose(
            dp.cartesian_coordinates,
            cur_molecule.get_coordinates_in_bohr(),
            atol=1.0e-5,
            rtol=0.0,
        ), f"Cartesian mismatch for {label}"

        assert float(dp.energy) == pytest.approx(float(e_ref), abs=1.0e-5)

        _assert_cartesian_tensor_shapes(dp, n_atoms=n_atoms)
        _rebuild_from_cartesian_and_compare(
            dp=dp,
            z_matrix=z_matrix,
            settings=settings,
            inv_sqrt_masses=inv_sqrt_masses,
        )

    _assert_roundtrip_write_read(
        ctx=ctx,
        z_matrix=z_matrix,
        settings=settings,
        labels=labels,
        datapoints=datapoints,
    )



def test_full_pipeline_single_file(tmp_path):
    """
    Single-run integration workflow with deep post-run database validation.

    Phase order:
    1) run construction once,
    2) setup/symmetry validation,
    3) database/datapoint integrity checks.
    """
    workdir = tmp_path / "full_pipeline"
    workdir.mkdir(parents=True, exist_ok=True)

    ctx = _run_construction_once(workdir)

    _assert_setup_symmetry(ctx)
    _assert_database_datapoint_integrity(ctx)
