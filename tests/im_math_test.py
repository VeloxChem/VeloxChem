import os

import numpy as np
import pytest

from veloxchem.molecule import Molecule
from veloxchem.interpolationdatapoint import InterpolationDatapoint
from veloxchem.interpolationdriver import InterpolationDriver
import veloxchem.imtrustradiusoptimizer as alphaoptimizer_module


def _reporting_enabled():
    flag = os.environ.get("IM_TEST_REPORT", "1").strip().lower()
    return flag not in {"0", "false", "no", "off"}


def _detailed_arrays_enabled():
    flag = os.environ.get("IM_TEST_REPORT_ARRAYS", "0").strip().lower()
    return flag not in {"0", "false", "no", "off"}


def _debug_enabled():
    flag = os.environ.get("IM_TEST_DEBUG", "0").strip().lower()
    return flag not in {"0", "false", "no", "off", ""}


def _fd_dump_enabled():
    flag = os.environ.get("IM_TEST_PRINT_FD", "0").strip().lower()
    return flag not in {"0", "false", "no", "off", ""}


def _debug_print(message):
    if _debug_enabled():
        print(f"[IM-DEBUG] {message}")


def _array_str(arr, precision=8):
    return np.array2string(
        np.asarray(arr),
        precision=precision,
        suppress_small=False,
        max_line_width=160,
    )


def _array_comparison_lines(title, analytical, numerical):
    a = np.asarray(analytical, dtype=np.float64)
    n = np.asarray(numerical, dtype=np.float64)
    d = a - n
    return [
        f"{title}:",
        f"analytical={_array_str(a)}",
        f"numerical={_array_str(n)}",
        f"delta={_array_str(d)}",
        f"delta_max_abs={np.max(np.abs(d)):.6e}",
        f"delta_l2={np.linalg.norm(d):.6e}",
    ]


def _print_fd_comparison(tag, analytical, numerical):
    a = np.asarray(analytical, dtype=np.float64)
    n = np.asarray(numerical, dtype=np.float64)
    d = a - n

    print(f"[IM-FD] {tag}")
    print(
        "[IM-FD] summary: "
        f"shape={a.shape}, "
        f"||analytical||={np.linalg.norm(a):.12e}, "
        f"||fd||={np.linalg.norm(n):.12e}, "
        f"||delta||={np.linalg.norm(d):.12e}, "
        f"max_abs_delta={np.max(np.abs(d)):.12e}"
    )
    print(f"[IM-FD] analytical={_array_str(a, precision=12)}")
    print(f"[IM-FD] fd={_array_str(n, precision=12)}")
    print(f"[IM-FD] delta={_array_str(d, precision=12)}")


def _emit_report(request, title, objective, checks):
    if not _reporting_enabled():
        return

    terminal = request.config.pluginmanager.get_plugin("terminalreporter")
    lines = [f"Objective: {objective}", "Checks:"]
    lines.extend(f"  - {check}" for check in checks)

    if terminal is not None:
        terminal.write_sep("-", f"IM-Test: {title}")
        for line in lines:
            terminal.write_line(str(line))
    else:
        print(f"\n[IM-Test] {title}")
        for line in lines:
            print(line)


def _wrap_to_pi(delta):
    return (np.asarray(delta, dtype=np.float64) + np.pi) % (2.0 * np.pi) - np.pi


def _base_molecule_system():
    xyz = """9

C             -0.661535066309         1.128097375722         1.298956571720
C              0.160830459103        -0.078956848612         0.840928833654
O             -0.545706317697        -1.293534712398         0.947688936392
H             -1.201152069039         1.577440238953         0.467697680596
H             -0.006154222315         1.890329398828         1.712599756591
H             -1.363731740381         0.823716932687         2.070509517949
H              1.028694858953        -0.210814780197         1.491718794709
H              0.521388572669         0.070466580703        -0.186332296889
H             -1.346953592312        -1.233309727019         0.415095469234"""

    molecule = Molecule.from_xyz_string(xyz)
    molecule.set_charge(0)
    molecule.set_multiplicity(1)

    z_matrix = {
        "bonds": [
            (0, 1), (0, 3), (0, 4), (0, 5),
            (1, 2), (1, 6), (1, 7), (2, 8),
        ],
        "angles": [
            (1, 0, 3), (1, 0, 4), (1, 0, 5),
            (3, 0, 4), (3, 0, 5), (4, 0, 5),
            (0, 1, 2), (0, 1, 6), (0, 1, 7),
            (2, 1, 6), (2, 1, 7), (6, 1, 7),
            (1, 2, 8),
        ],
        "dihedrals": [
            (3, 0, 1, 2), (3, 0, 1, 6), (3, 0, 1, 7),
            (4, 0, 1, 2), (4, 0, 1, 6), (4, 0, 1, 7),
            (5, 0, 1, 2), (5, 0, 1, 6), (5, 0, 1, 7),
            (0, 1, 2, 8), (6, 1, 2, 8), (7, 1, 2, 8),
        ],
        "impropers": [],
    }

    return molecule, z_matrix


def _identity_runtime_metadata(molecule, z_matrix):
    n_atoms = len(molecule.get_labels())
    dihedral_start = len(z_matrix["bonds"]) + len(z_matrix["angles"])
    dihedral_end = dihedral_start + len(z_matrix["dihedrals"])

    return [
        list(range(n_atoms)),
        [],
        [],
        list(range(n_atoms)),
        [],
        [],
        [],
        {2: {}, 3: {}},
        {},
        [dihedral_start, dihedral_end],
    ]


def _settings(weightfunction_type="cartesian", use_tc_weights=False):
    return {
        "interpolation_type": "shepard",
        "weightfunction_type": str(weightfunction_type),
        "exponent_p": 2,
        "exponent_q": 2,
        "confidence_radius": 0.5,
        "use_inverse_bond_length": True,
        "use_cosine_dihedral": False,
        "use_tc_weights": bool(use_tc_weights),
    }


def _settings_with_bond_mode(settings, bond_mode):
    patched = dict(settings)
    if bond_mode == "inverse":
        patched["use_inverse_bond_length"] = True
    elif bond_mode == "plain_r":
        patched["use_inverse_bond_length"] = False
    else:
        raise ValueError(f"Unsupported bond_mode={bond_mode}")
    return patched


def _make_molecule_like(template_molecule, coordinates_bohr):
    mol = Molecule(template_molecule.get_labels(), coordinates_bohr, "bohr")
    mol.set_charge(template_molecule.get_charge())
    mol.set_multiplicity(template_molecule.get_multiplicity())
    return mol


def _coords_with_dihedral_shift(template_mol, dihedral, delta_rad):
    mol = _make_molecule_like(
        template_mol,
        template_mol.get_coordinates_in_bohr().copy(),
    )
    d = tuple(int(i) + 1 for i in dihedral)
    phi0 = mol.get_dihedral(d, "radian")
    mol.set_dihedral(d, phi0 + float(delta_rad), "radian")
    return mol.get_coordinates_in_bohr().copy()


def _surface_values(coords, center):
    x = np.asarray(coords, dtype=np.float64).reshape(-1)
    c = np.asarray(center, dtype=np.float64).reshape(-1)

    force_constants = np.linspace(0.015, 0.035, x.size)
    dx = x - c

    energy = -11.0 + 0.5 * float(np.dot(force_constants * dx, dx))
    gradient = (force_constants * dx).reshape(np.asarray(coords).shape)
    hessian = np.diag(force_constants)

    return energy, gradient, hessian


def _important_coordinates():
    return {
        "bonds": [(0, 1)],
        "angles": [],
        "dihedrals": [(3, 0, 1, 2)],
        "impropers": [],
    }


def _empty_important_coordinates():
    return {"bonds": [], "angles": [], "dihedrals": [], "impropers": []}


def _make_datapoint(
    molecule,
    z_matrix,
    settings,
    coords,
    center,
    label,
    confidence_radius,
    important=True,
):
    inv_sqrt_masses = 1.0 / np.sqrt(np.repeat(molecule.get_masses(), 3))
    energy, gradient, hessian = _surface_values(coords, center)

    dp = InterpolationDatapoint(z_matrix)
    dp.update_settings(settings)
    dp.cartesian_coordinates = np.array(coords, dtype=np.float64, copy=True)
    dp.inv_sqrt_masses = inv_sqrt_masses
    dp.energy = float(energy)
    dp.gradient = (inv_sqrt_masses * gradient.reshape(-1)).reshape(gradient.shape)
    dp.hessian = inv_sqrt_masses[:, None] * hessian * inv_sqrt_masses[None, :]
    dp.confidence_radius = float(confidence_radius)
    dp.point_label = str(label)
    dp.imp_int_coordinates = (
        _important_coordinates() if important else _empty_important_coordinates()
    )
    dp.transform_gradient_and_hessian()

    return dp


def _build_current_framework_case(
    weightfunction_type="cartesian",
    use_tc_weights=False,
    important=True,
):
    base_mol, z_matrix = _base_molecule_system()
    settings = _settings(weightfunction_type, use_tc_weights)
    runtime_metadata = _identity_runtime_metadata(base_mol, z_matrix)

    center = base_mol.get_coordinates_in_bohr().copy()
    coords_1 = center.copy()
    coords_2 = _coords_with_dihedral_shift(base_mol, (3, 0, 1, 2), 0.30)
    coords_3 = _coords_with_dihedral_shift(base_mol, (0, 1, 2, 8), -0.22)
    coords_3[6, 2] += 0.025

    datapoints = [
        _make_datapoint(base_mol, z_matrix, settings, coords_1, center, "point_1", 0.50, important),
        _make_datapoint(base_mol, z_matrix, settings, coords_2, center, "point_2", 0.55, important),
        _make_datapoint(base_mol, z_matrix, settings, coords_3, center, "point_3", 0.60, important),
    ]

    driver = InterpolationDriver(z_matrix)
    driver.update_settings(settings)
    driver.symmetry_information = runtime_metadata
    driver.impes_coordinate.inv_sqrt_masses = datapoints[0].inv_sqrt_masses
    driver.qm_data_points = datapoints
    driver.labels = [dp.point_label for dp in datapoints]

    probe = _coords_with_dihedral_shift(base_mol, (3, 0, 1, 2), 0.16)
    probe[4, 1] += 0.015

    return base_mol, z_matrix, settings, runtime_metadata, datapoints, driver, probe, center


def _periodic_internal_rows(dp):
    if getattr(dp, "use_cosine_dihedral", False):
        return np.empty(0, dtype=np.int64)

    return np.asarray(
        [i for i, z_entry in enumerate(dp.z_matrix) if len(z_entry) == 4],
        dtype=np.int64,
    )


def _build_fd_validation_datapoint(base_molecule, z_matrix, settings, bond_mode, mass_weighted):
    local_settings = _settings_with_bond_mode(settings, bond_mode)

    dp = InterpolationDatapoint(z_matrix)
    dp.update_settings(local_settings)
    dp.cartesian_coordinates = np.array(
        base_molecule.get_coordinates_in_bohr(),
        dtype=np.float64,
        copy=True,
    )
    dp.define_internal_coordinates()

    if mass_weighted:
        dp.inv_sqrt_masses = 1.0 / np.sqrt(np.repeat(base_molecule.get_masses(), 3))

    return dp


def _xprime_to_cartesian_flat(dp, xprime):
    xprime = np.asarray(xprime, dtype=np.float64).reshape(-1)
    if dp.inv_sqrt_masses is None:
        return xprime.copy()
    return xprime * dp.inv_sqrt_masses


def _evaluate_internal_values_at_xprime(dp, xprime):
    x_cart = _xprime_to_cartesian_flat(dp, xprime)
    base_values = np.fromiter(
        (q.value(x_cart) for q in dp.internal_coordinates),
        dtype=np.float64,
        count=len(dp.internal_coordinates),
    )

    values = base_values.copy()
    row_cache = dp._get_b_matrix_row_cache()

    bond_rows = row_cache["bond_rows"]
    if bond_rows.size > 0 and dp.use_inverse_bond_length:
        values[bond_rows] = 1.0 / base_values[bond_rows]

    if dp.use_cosine_dihedral:
        dihedral_rows = row_cache["dihedral_rows"]
        if dihedral_rows.size > 0:
            phi = base_values[dihedral_rows]
            first_mask = row_cache["dihedral_first"]
            values[dihedral_rows[first_mask]] = np.cos(phi[first_mask])
            values[dihedral_rows[~first_mask]] = np.sin(phi[~first_mask])

    return values


def _finite_difference_b_matrix(dp, step=2.0e-6):
    x_cart = np.asarray(dp.cartesian_coordinates, dtype=np.float64).reshape(-1)
    xprime0 = x_cart.copy() if dp.inv_sqrt_masses is None else x_cart / dp.inv_sqrt_masses

    n_internal = len(dp.z_matrix)
    n_cart = x_cart.size
    fd = np.zeros((n_internal, n_cart), dtype=np.float64)

    q0 = _evaluate_internal_values_at_xprime(dp, xprime0)
    periodic_rows = _periodic_internal_rows(dp)

    for a in range(n_cart):
        e = np.zeros(n_cart, dtype=np.float64)
        e[a] = 1.0

        q_plus = _evaluate_internal_values_at_xprime(dp, xprime0 + step * e)
        q_minus = _evaluate_internal_values_at_xprime(dp, xprime0 - step * e)
        dq = q_plus - q_minus

        if periodic_rows.size > 0:
            d_plus = _wrap_to_pi(q_plus[periodic_rows] - q0[periodic_rows])
            d_minus = _wrap_to_pi(q_minus[periodic_rows] - q0[periodic_rows])
            dq[periodic_rows] = d_plus - d_minus

        fd[:, a] = dq / (2.0 * step)

    return fd


def _finite_difference_b2_matrix(dp, step=1.0e-5):
    x_cart = np.asarray(dp.cartesian_coordinates, dtype=np.float64).reshape(-1)
    xprime0 = x_cart.copy() if dp.inv_sqrt_masses is None else x_cart / dp.inv_sqrt_masses

    n_internal = len(dp.z_matrix)
    n_cart = x_cart.size
    fd2 = np.zeros((n_internal, n_cart, n_cart), dtype=np.float64)

    q0 = _evaluate_internal_values_at_xprime(dp, xprime0)
    periodic_rows = _periodic_internal_rows(dp)

    for a in range(n_cart):
        ea = np.zeros(n_cart, dtype=np.float64)
        ea[a] = 1.0

        q_plus = _evaluate_internal_values_at_xprime(dp, xprime0 + step * ea)
        q_minus = _evaluate_internal_values_at_xprime(dp, xprime0 - step * ea)
        diag = (q_plus - 2.0 * q0 + q_minus) / step**2

        if periodic_rows.size > 0:
            d_plus = _wrap_to_pi(q_plus[periodic_rows] - q0[periodic_rows])
            d_minus = _wrap_to_pi(q_minus[periodic_rows] - q0[periodic_rows])
            diag[periodic_rows] = (d_plus + d_minus) / step**2

        fd2[:, a, a] = diag

        for b in range(a + 1, n_cart):
            eb = np.zeros(n_cart, dtype=np.float64)
            eb[b] = 1.0

            q_pp = _evaluate_internal_values_at_xprime(dp, xprime0 + step * ea + step * eb)
            q_pm = _evaluate_internal_values_at_xprime(dp, xprime0 + step * ea - step * eb)
            q_mp = _evaluate_internal_values_at_xprime(dp, xprime0 - step * ea + step * eb)
            q_mm = _evaluate_internal_values_at_xprime(dp, xprime0 - step * ea - step * eb)

            mixed = (q_pp - q_pm - q_mp + q_mm) / (4.0 * step**2)

            if periodic_rows.size > 0:
                d_pp = _wrap_to_pi(q_pp[periodic_rows] - q0[periodic_rows])
                d_pm = _wrap_to_pi(q_pm[periodic_rows] - q0[periodic_rows])
                d_mp = _wrap_to_pi(q_mp[periodic_rows] - q0[periodic_rows])
                d_mm = _wrap_to_pi(q_mm[periodic_rows] - q0[periodic_rows])
                mixed[periodic_rows] = (d_pp - d_pm - d_mp + d_mm) / (4.0 * step**2)

            fd2[:, a, b] = mixed
            fd2[:, b, a] = mixed

    return fd2


def _evaluate_raw_weight_terms(driver, datapoint, weight_mode):
    if weight_mode == "cartesian":
        return driver.cartesian_distance(datapoint)
    if weight_mode == "internal":
        return driver.internal_distance(datapoint)
    raise ValueError(f"Unknown weight_mode={weight_mode}")


def _finite_difference_raw_weight_gradient_mode(driver, datapoint, coordinates, step, weight_mode):
    fd = np.zeros_like(coordinates, dtype=np.float64)
    excluded_atoms = set(driver.symmetry_information[4])

    for atom_idx in range(coordinates.shape[0]):
        if atom_idx in excluded_atoms:
            continue

        for axis_idx in range(3):
            coords_p = np.array(coordinates, copy=True)
            coords_m = np.array(coordinates, copy=True)
            coords_p[atom_idx, axis_idx] += step
            coords_m[atom_idx, axis_idx] -= step

            driver.define_impes_coordinate(coords_p)
            _, _, denom_p, _, _, _, _, _ = _evaluate_raw_weight_terms(
                driver,
                datapoint,
                weight_mode,
            )

            driver.define_impes_coordinate(coords_m)
            _, _, denom_m, _, _, _, _, _ = _evaluate_raw_weight_terms(
                driver,
                datapoint,
                weight_mode,
            )

            fd[atom_idx, axis_idx] = (1.0 / denom_p - 1.0 / denom_m) / (2.0 * step)

    driver.define_impes_coordinate(coordinates)
    return fd


def _finite_difference_energy_gradient(driver, template_molecule, coordinates, step):
    fd = np.zeros_like(coordinates, dtype=np.float64)

    for atom_idx in range(coordinates.shape[0]):
        for axis_idx in range(3):
            coords_p = np.array(coordinates, copy=True)
            coords_m = np.array(coordinates, copy=True)
            coords_p[atom_idx, axis_idx] += step
            coords_m[atom_idx, axis_idx] -= step

            driver.compute(_make_molecule_like(template_molecule, coords_p))
            e_p = float(driver.get_energy())

            driver.compute(_make_molecule_like(template_molecule, coords_m))
            e_m = float(driver.get_energy())

            fd[atom_idx, axis_idx] = (e_p - e_m) / (2.0 * step)

    driver.compute(_make_molecule_like(template_molecule, coordinates))
    return fd


@pytest.mark.parametrize("bond_mode", ["inverse", "plain_r"])
@pytest.mark.parametrize("mass_weighted", [False, True])
def test_interpolationdatapoint_b_matrix_matches_finite_difference(
    request,
    bond_mode,
    mass_weighted,
):
    base_mol, z_matrix = _base_molecule_system()
    settings = _settings()

    dp = _build_fd_validation_datapoint(
        base_mol,
        z_matrix,
        settings,
        bond_mode=bond_mode,
        mass_weighted=mass_weighted,
    )
    dp.calculate_b_matrix()

    b_fd = _finite_difference_b_matrix(dp, step=2.0e-6)

    abs_err = np.max(np.abs(dp.b_matrix - b_fd))
    rel_err = np.linalg.norm(dp.b_matrix - b_fd) / (np.linalg.norm(b_fd) + 1.0e-12)

    _emit_report(
        request,
        "InterpolationDatapoint B-Matrix FD Validation",
        "Analytical Wilson B matrix must match finite-difference internal-coordinate Jacobian.",
        [
            f"bond_mode={bond_mode}",
            f"mass_weighted={mass_weighted}",
            f"max_abs_error={abs_err:.6e}",
            f"relative_error={rel_err:.6e}",
        ] + (
            _array_comparison_lines("B analytical vs FD", dp.b_matrix, b_fd)
            if _detailed_arrays_enabled()
            else []
        ),
    )

    assert abs_err < 2.0e-5
    assert rel_err < 2.0e-5


@pytest.mark.parametrize("bond_mode", ["inverse", "plain_r"])
@pytest.mark.parametrize("mass_weighted", [False, True])
def test_interpolationdatapoint_b2_matrix_matches_finite_difference(
    request,
    bond_mode,
    mass_weighted,
):
    base_mol, z_matrix = _base_molecule_system()
    settings = _settings()

    dp = _build_fd_validation_datapoint(
        base_mol,
        z_matrix,
        settings,
        bond_mode=bond_mode,
        mass_weighted=mass_weighted,
    )
    dp.calculate_b_matrix()
    dp.calculate_b2_matrix()

    b2_fd = _finite_difference_b2_matrix(dp, step=1.0e-5)

    abs_err = np.max(np.abs(dp.b2_matrix - b2_fd))
    rel_err = np.linalg.norm(dp.b2_matrix - b2_fd) / (np.linalg.norm(b2_fd) + 1.0e-12)

    _emit_report(
        request,
        "InterpolationDatapoint B2-Matrix FD Validation",
        "Analytical Wilson B2 tensor must match finite-difference internal-coordinate Hessian.",
        [
            f"bond_mode={bond_mode}",
            f"mass_weighted={mass_weighted}",
            f"max_abs_error={abs_err:.6e}",
            f"relative_error={rel_err:.6e}",
        ] + (
            _array_comparison_lines("B2 analytical vs FD", dp.b2_matrix, b2_fd)
            if _detailed_arrays_enabled()
            else []
        ),
    )

    assert abs_err < 3.0e-4
    assert rel_err < 5.0e-4


@pytest.mark.parametrize(
    "weight_mode,use_tc_weights,threshold",
    [
        ("cartesian", False, 1.0e-5),
        ("cartesian", True, 1.0e-4),
        ("internal", False, 5.0e-5),
    ],
)
def test_shepard_raw_weight_gradient_matches_fd(
    request,
    weight_mode,
    use_tc_weights,
    threshold,
):
    (
        _base_mol,
        _z_matrix,
        _settings_dict,
        _runtime_metadata,
        datapoints,
        driver,
        probe,
        _center,
    ) = _build_current_framework_case(
        weightfunction_type=weight_mode,
        use_tc_weights=use_tc_weights,
    )

    target_dp = datapoints[1]
    driver.define_impes_coordinate(probe)

    _, _, _, analytic_grad, _, _, _, _ = _evaluate_raw_weight_terms(
        driver,
        target_dp,
        weight_mode,
    )
    fd_grad = _finite_difference_raw_weight_gradient_mode(
        driver,
        target_dp,
        probe,
        step=1.0e-5,
        weight_mode=weight_mode,
    )

    rel_err = np.linalg.norm(analytic_grad - fd_grad) / (np.linalg.norm(fd_grad) + 1.0e-12)

    if _fd_dump_enabled() or _debug_enabled() or rel_err >= threshold:
        _print_fd_comparison(
            f"Raw weight-gradient FD comparison ({weight_mode}, tc={use_tc_weights})",
            analytic_grad,
            fd_grad,
        )

    _emit_report(
        request,
        "Raw Shepard Weight Gradient",
        "Analytical raw weight gradient must match central finite differences.",
        [
            f"weight_mode={weight_mode}",
            f"use_tc_weights={use_tc_weights}",
            f"datapoint_label={target_dp.point_label}",
            f"relative_error={rel_err:.6e}",
            f"threshold={threshold:.6e}",
        ],
    )

    assert rel_err < threshold


@pytest.mark.parametrize(
    "weight_mode,use_tc_weights",
    [
        ("cartesian", False),
        ("cartesian", True),
        ("internal", False),
    ],
)
def test_shepard_normalized_weight_energy_and_gradient_identities(
    request,
    weight_mode,
    use_tc_weights,
):
    (
        base_mol,
        _z_matrix,
        _settings_dict,
        _runtime_metadata,
        datapoints,
        driver,
        probe,
        _center,
    ) = _build_current_framework_case(
        weightfunction_type=weight_mode,
        use_tc_weights=use_tc_weights,
    )

    driver.compute(_make_molecule_like(base_mol, probe))
    driver.define_impes_coordinate(probe)

    raw_weights = []
    raw_weight_grads = []

    for dp in datapoints:
        _, _, denominator, weight_gradient, _, _, _, _ = _evaluate_raw_weight_terms(
            driver,
            dp,
            weight_mode,
        )
        raw_weights.append(1.0 / denominator)
        raw_weight_grads.append(weight_gradient)

    raw_weights = np.asarray(raw_weights, dtype=np.float64)
    raw_weight_grads = np.asarray(raw_weight_grads, dtype=np.float64)

    sum_w = raw_weights.sum()
    sum_grad_w = raw_weight_grads.sum(axis=0)

    normalized_w = raw_weights / sum_w
    normalized_grad_w = (
        raw_weight_grads * sum_w
        - raw_weights[:, None, None] * sum_grad_w
    ) / sum_w**2

    potentials = np.asarray(driver.potentials, dtype=np.float64)
    gradients = np.asarray(driver.gradients, dtype=np.float64)

    energy_reconstructed = float(np.dot(normalized_w, potentials))
    gradient_reconstructed = (
        np.tensordot(normalized_w, gradients, axes=1)
        + np.tensordot(potentials, normalized_grad_w, axes=1)
    )

    assert np.isclose(normalized_w.sum(), 1.0, atol=1.0e-12, rtol=0.0)
    assert np.isclose(energy_reconstructed, driver.get_energy(), atol=1.0e-12, rtol=0.0)
    assert np.allclose(
        gradient_reconstructed,
        driver.get_gradient(),
        atol=1.0e-12,
        rtol=0.0,
    )

    _emit_report(
        request,
        "Normalized Weight Identities",
        "Validate sum(W)=1, E=sum(WU), and gradient assembly from normalized weights.",
        [
            f"weight_mode={weight_mode}",
            f"use_tc_weights={use_tc_weights}",
            f"normalized_weights={_array_str(normalized_w)}",
            f"|E_reconstructed-E_model|={abs(energy_reconstructed - driver.get_energy()):.6e}",
            f"||G_reconstructed-G_model||={np.linalg.norm(gradient_reconstructed - driver.get_gradient()):.6e}",
        ],
    )


@pytest.mark.parametrize(
    "weight_mode,use_tc_weights,threshold",
    [
        ("cartesian", False, 2.0e-2),
        ("cartesian", True, 3.0e-2),
        ("internal", False, 3.0e-2),
    ],
)
def test_shepard_energy_gradient_matches_fd_current_framework(
    request,
    weight_mode,
    use_tc_weights,
    threshold,
):
    (
        base_mol,
        _z_matrix,
        _settings_dict,
        _runtime_metadata,
        _datapoints,
        driver,
        probe,
        _center,
    ) = _build_current_framework_case(
        weightfunction_type=weight_mode,
        use_tc_weights=use_tc_weights,
    )

    driver.compute(_make_molecule_like(base_mol, probe))
    model_grad = np.asarray(driver.get_gradient(), dtype=np.float64)

    fd_grad = _finite_difference_energy_gradient(
        driver,
        base_mol,
        probe,
        step=2.0e-5,
    )

    rel_err = np.linalg.norm(model_grad - fd_grad) / (np.linalg.norm(fd_grad) + 1.0e-12)

    if _fd_dump_enabled() or _debug_enabled() or rel_err >= threshold:
        _print_fd_comparison(
            f"Energy-gradient FD comparison ({weight_mode}, tc={use_tc_weights})",
            model_grad,
            fd_grad,
        )

    _emit_report(
        request,
        "Energy->Gradient FD Consistency",
        "Finite-difference energy derivative must match interpolated Cartesian gradient.",
        [
            f"weight_mode={weight_mode}",
            f"use_tc_weights={use_tc_weights}",
            f"relative_error={rel_err:.6e}",
            f"threshold={threshold:.6e}",
        ],
    )

    assert rel_err < threshold


def test_target_customized_cartesian_weight_gate_is_included():
    (
        _base_mol_plain,
        _z_matrix_plain,
        _settings_plain,
        _runtime_plain,
        datapoints_plain,
        driver_plain,
        probe,
        _center_plain,
    ) = _build_current_framework_case(
        weightfunction_type="cartesian",
        use_tc_weights=False,
        important=True,
    )

    (
        _base_mol_tc,
        _z_matrix_tc,
        _settings_tc,
        _runtime_tc,
        datapoints_tc,
        driver_tc,
        _probe_tc,
        _center_tc,
    ) = _build_current_framework_case(
        weightfunction_type="cartesian",
        use_tc_weights=True,
        important=True,
    )

    target_plain = datapoints_plain[1]
    target_tc = datapoints_tc[1]

    driver_plain.define_impes_coordinate(probe)
    driver_tc.define_impes_coordinate(probe)

    _, _, denominator_base, gradient_base, _, _, _, _ = driver_plain.cartesian_distance(target_plain)
    _, _, denominator_tc, gradient_tc, _, _, _, _ = driver_tc.cartesian_distance(target_tc)

    active_dofs = np.arange(probe.size, dtype=np.int64)
    A_imp, grad_A_imp = driver_tc._important_coordinate_gate_metric(target_tc, active_dofs)

    expected_denominator, expected_gradient = driver_tc._apply_imp_coordinate_penalty_gate(
        denominator_base,
        gradient_base,
        A_imp,
        grad_A_imp,
    )

    assert A_imp > 0.0
    assert np.isfinite(denominator_tc)
    assert 1.0 / denominator_tc < 1.0 / denominator_base
    assert np.isclose(denominator_tc, expected_denominator, atol=1.0e-12, rtol=1.0e-12)
    assert np.allclose(gradient_tc, expected_gradient, atol=1.0e-12, rtol=1.0e-12)


def test_target_customized_cartesian_weights_are_noop_without_important_coordinates():
    (
        _base_mol_plain,
        _z_matrix_plain,
        _settings_plain,
        _runtime_plain,
        datapoints_plain,
        driver_plain,
        probe,
        _center_plain,
    ) = _build_current_framework_case(
        weightfunction_type="cartesian",
        use_tc_weights=False,
        important=False,
    )

    (
        _base_mol_tc,
        _z_matrix_tc,
        _settings_tc,
        _runtime_tc,
        datapoints_tc,
        driver_tc,
        _probe_tc,
        _center_tc,
    ) = _build_current_framework_case(
        weightfunction_type="cartesian",
        use_tc_weights=True,
        important=False,
    )

    driver_plain.define_impes_coordinate(probe)
    driver_tc.define_impes_coordinate(probe)

    for dp_plain, dp_tc in zip(datapoints_plain, datapoints_tc):
        _, _, denominator_plain, gradient_plain, _, _, _, _ = driver_plain.cartesian_distance(dp_plain)
        _, _, denominator_tc, gradient_tc, _, _, _, _ = driver_tc.cartesian_distance(dp_tc)

        assert np.isclose(denominator_tc, denominator_plain, atol=1.0e-12, rtol=1.0e-12)
        assert np.allclose(gradient_tc, gradient_plain, atol=1.0e-12, rtol=1.0e-12)


@pytest.mark.parametrize("use_tc_weights", [False, True])
def test_trust_radius_optimizer_jacobian_matches_fd_without_symmetry(
    request,
    use_tc_weights,
):
    (
        base_mol,
        z_matrix,
        settings,
        runtime_metadata,
        datapoints,
        _driver,
        probe,
        center,
    ) = _build_current_framework_case(
        weightfunction_type="cartesian",
        use_tc_weights=use_tc_weights,
        important=True,
    )

    coords_list = [
        probe,
        probe + np.eye(*probe.shape, k=0) * 3.0e-4,
        _coords_with_dihedral_shift(base_mol, (3, 0, 1, 2), 0.20),
    ]

    structures = [_make_molecule_like(base_mol, coords) for coords in coords_list]

    qm_energies = []
    qm_gradients = []
    for coords in coords_list:
        e, g, _ = _surface_values(coords, center)
        qm_energies.append(e)
        qm_gradients.append(g)

    optimizer = alphaoptimizer_module.IMTrustRadiusOptimizer(
        z_matrix,
        settings,
        runtime_metadata,
        datapoints,
        structures,
        qm_energies,
        qm_gradients,
        exponent_p_q=(2, 2),
        e_x=0.5,
        beta=0.8,
        verbose=False,
    )

    try:
        alphas = np.asarray(
            [float(dp.confidence_radius) for dp in datapoints],
            dtype=np.float64,
        )

        analytic_jac = np.asarray(optimizer.jac(alphas), dtype=np.float64)
        fd_jac = np.zeros_like(alphas)

        step = 1.0e-4
        for idx in range(alphas.size):
            plus = alphas.copy()
            minus = alphas.copy()
            plus[idx] += step
            minus[idx] -= step
            fd_jac[idx] = (optimizer.fun(plus) - optimizer.fun(minus)) / (2.0 * step)

        rel_err = np.linalg.norm(analytic_jac - fd_jac) / (np.linalg.norm(fd_jac) + 1.0e-12)

    finally:
        optimizer.close()

    if _fd_dump_enabled() or _debug_enabled() or rel_err >= 5.0e-4:
        _print_fd_comparison(
            f"Trust-radius optimizer Jacobian FD comparison (tc={use_tc_weights})",
            analytic_jac,
            fd_jac,
        )

    _emit_report(
        request,
        "Trust-Radius Optimizer Jacobian",
        "Analytic confidence-radius Jacobian must match finite differences without symmetry banks.",
        [
            f"use_tc_weights={use_tc_weights}",
            f"n_datapoints={len(datapoints)}",
            f"n_alphas={alphas.size}",
            f"relative_error={rel_err:.6e}",
            "threshold=5.000000e-04",
        ],
    )

    assert rel_err < 5.0e-4
