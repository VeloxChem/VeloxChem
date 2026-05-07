import importlib.util
from pathlib import Path
import sys

import numpy as np
import pytest

from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.scfunrestdriver import ScfUnrestrictedDriver

ANALYZER_PATH = (Path(__file__).resolve().parents[1] / 'src' / 'pymodule' /
                 'orbitalanalyzerdriver.py')
ANALYZER_SPEC = importlib.util.spec_from_file_location(
    'veloxchem.orbitalanalyzerdriver', ANALYZER_PATH)
ANALYZER_MODULE = importlib.util.module_from_spec(ANALYZER_SPEC)
sys.modules['veloxchem.orbitalanalyzerdriver'] = ANALYZER_MODULE
ANALYZER_SPEC.loader.exec_module(ANALYZER_MODULE)

VBDRIVER_PATH = (Path(__file__).resolve().parents[1] / 'src' / 'pymodule' /
                 'vbdriver.py')
VBDRIVER_SPEC = importlib.util.spec_from_file_location('veloxchem.vbdriver',
                                                       VBDRIVER_PATH)
VBDRIVER_MODULE = importlib.util.module_from_spec(VBDRIVER_SPEC)
sys.modules['veloxchem.vbdriver'] = VBDRIVER_MODULE
VBDRIVER_SPEC.loader.exec_module(VBDRIVER_MODULE)

VbComputeOptions = VBDRIVER_MODULE.VbComputeOptions
VbDriver = VBDRIVER_MODULE.VbDriver

HARTREE_TO_KJMOL = 2625.499638


def _h2_molecule(distance):
    molecule = Molecule.read_str(
        f"""
H 0.0 0.0 0.0
H {float(distance):.8f} 0.0 0.0
"""
    )
    molecule.set_charge(0)
    molecule.set_multiplicity(1)
    return molecule


def _run_h2_hf_reference(molecule, basis, unrestricted=False):
    driver = ScfUnrestrictedDriver() if unrestricted else ScfRestrictedDriver()
    driver.ostream.mute()
    driver.xcfun = "hf"
    if unrestricted:
        driver.guess_unpaired_electrons = "1(1),2(-1)"
    driver.compute(molecule, basis)
    return float(driver.get_scf_energy())


def _run_hf_reference(molecule, basis, unrestricted=False):
    return _run_h2_hf_reference(molecule, basis, unrestricted=unrestricted)


def _relative_potential_curve(energies):
    energies = np.asarray(energies, dtype=float)
    return (energies - energies[-1]) * HARTREE_TO_KJMOL


def _run_h2_vb(distance, mode="vbci", include_bovb=False, basis_name="sto-3g"):
    molecule = _h2_molecule(distance)
    basis = MolecularBasis.read(molecule, basis_name, ostream=None)
    result = VbDriver().compute(
        molecule,
        basis,
        options=VbComputeOptions(
            mode=mode,
            optimize_orbitals=(mode == "vbscf"),
            include_bovb=include_bovb,
        ),
    )
    return molecule, basis, result


def _ethylene_molecule(distance=1.339):
    half_distance = 0.5 * float(distance)
    molecule = Molecule.read_str(f"""
C {-half_distance:.8f}  0.0000  0.0000
C  {half_distance:.8f}  0.0000  0.0000
H {-half_distance - 0.5626:.8f}  0.9289  0.0000
H {-half_distance - 0.5626:.8f} -0.9289  0.0000
H  {half_distance + 0.5626:.8f}  0.9289  0.0000
H  {half_distance + 0.5626:.8f} -0.9289  0.0000
""")
    molecule.set_charge(0)
    molecule.set_multiplicity(1)
    return molecule


def _ethane_molecule(distance=1.54):
    half_distance = 0.5 * float(distance)
    molecule = Molecule.read_str(f"""
C {-half_distance:.8f}  0.0000  0.0000
C  {half_distance:.8f}  0.0000  0.0000
H {-half_distance - 0.6300:.8f}  0.9000  0.0000
H {-half_distance - 0.6300:.8f} -0.4500  0.77942286
H {-half_distance - 0.6300:.8f} -0.4500 -0.77942286
H  {half_distance + 0.6300:.8f} -0.9000  0.0000
H  {half_distance + 0.6300:.8f}  0.4500  0.77942286
H  {half_distance + 0.6300:.8f}  0.4500 -0.77942286
""")
    molecule.set_charge(0)
    molecule.set_multiplicity(1)
    return molecule


def _run_ethylene_pi_vb(mode="vbci", include_bovb=False):
    molecule = _ethylene_molecule()
    basis = MolecularBasis.read(molecule, "sto-3g", ostream=None)
    result = VbDriver().compute(
        molecule,
        basis,
        options=VbComputeOptions(
            mode=mode,
            optimize_orbitals=False,
            include_bovb=include_bovb,
            active_candidate_subtype="pi",
            include_ionic=True,
            freeze_inactive_orbitals=True,
        ),
    )
    return molecule, basis, result


def _run_one_bond_vb(molecule,
                     basis,
                     mode="vbci",
                     active_subtype="sigma",
                     reference_orbitals=None):
    return VbDriver().compute(
        molecule,
        basis,
        options=VbComputeOptions(
            mode=mode,
            optimize_orbitals=(mode == "vbscf"),
            include_bovb=(mode == "bovb"),
            active_bond=(0, 1),
            active_candidate_subtype=active_subtype,
            active_bond_reference_orbitals=reference_orbitals,
            include_ionic=True,
            freeze_inactive_orbitals=True,
        ),
    )


ALLYL_CATION_CASE = {
    "name": "allyl_cation",
    "xyz": """
C -1.3000  0.0000  0.0000
C  0.0000  0.0000  0.0000
C  1.3000  0.0000  0.0000
H -1.8500  0.9200  0.0000
H -1.8500 -0.9200  0.0000
H  0.0000  1.0800  0.0000
H  1.8500  0.9200  0.0000
H  1.8500 -0.9200  0.0000
""",
    "charge": 1,
    "multiplicity": 1,
    "active_pi_atoms": (0, 1, 2),
    "active_electron_count": 2,
    "active_spin": "singlet",
}


PI_CASES = [
    {
        "name": "allyl_radical",
        "xyz": """
C -1.3000  0.0000  0.0000
C  0.0000  0.0000  0.0000
C  1.3000  0.0000  0.0000
H -1.8500  0.9200  0.0000
H -1.8500 -0.9200  0.0000
H  0.0000  1.0800  0.0000
H  1.8500  0.9200  0.0000
H  1.8500 -0.9200  0.0000
""",
        "charge": 0,
        "multiplicity": 2,
        "active_pi_atoms": (0, 1, 2),
        "active_electron_count": 3,
        "active_spin": "doublet",
        "expected_templates": 3,
        "expected_types": {"allyl_radical"},
        "expected_determinants": 9,
    },
    {
        "name": "allyl_anion",
        "xyz": """
C -1.3000  0.0000  0.0000
C  0.0000  0.0000  0.0000
C  1.3000  0.0000  0.0000
H -1.8500  0.9200  0.0000
H -1.8500 -0.9200  0.0000
H  0.0000  1.0800  0.0000
H  1.8500  0.9200  0.0000
H  1.8500 -0.9200  0.0000
""",
        "charge": -1,
        "multiplicity": 1,
        "active_pi_atoms": (0, 1, 2),
        "active_electron_count": 4,
        "active_spin": "singlet",
        "expected_templates": 3,
        "expected_types": {"allyl_anion"},
        "expected_determinants": 9,
    },
    {
        "name": "butadiene",
        "xyz": """
C -1.9950  0.0000  0.0000
C -0.6650  0.0000  0.0000
C  0.6650  0.0000  0.0000
C  1.9950  0.0000  0.0000
H -2.5500  0.9200  0.0000
H -2.5500 -0.9200  0.0000
H -0.6650  1.0800  0.0000
H  0.6650 -1.0800  0.0000
H  2.5500  0.9200  0.0000
H  2.5500 -0.9200  0.0000
""",
        "charge": 0,
        "multiplicity": 1,
        "active_pi_atoms": (0, 1, 2, 3),
        "active_electron_count": 4,
        "active_spin": "singlet",
        "expected_templates": 3,
        "expected_types": {"butadiene_kekule", "butadiene_long_bond_pairing"},
        "expected_determinants": 36,
        "min_capture": 1.0e-3,
    },
    {
        "name": "benzene",
        "xyz": """
C  1.3970  0.0000  0.0000
C  0.6985  1.2099  0.0000
C -0.6985  1.2099  0.0000
C -1.3970  0.0000  0.0000
C -0.6985 -1.2099  0.0000
C  0.6985 -1.2099  0.0000
H  2.4810  0.0000  0.0000
H  1.2405  2.1487  0.0000
H -1.2405  2.1487  0.0000
H -2.4810  0.0000  0.0000
H -1.2405 -2.1487  0.0000
H  1.2405 -2.1487  0.0000
""",
        "charge": 0,
        "multiplicity": 1,
        "active_pi_atoms": (0, 1, 2, 3, 4, 5),
        "active_electron_count": 6,
        "active_spin": "singlet",
        "expected_templates": 15,
        "expected_types": {"benzene_kekule", "benzene_dewar"},
        "expected_determinants": 400,
    },
]


def _run_pi_case(case):
    molecule = Molecule.read_str(case["xyz"])
    molecule.set_charge(case["charge"])
    molecule.set_multiplicity(case["multiplicity"])
    basis = MolecularBasis.read(molecule, "sto-3g", ostream=None)

    options = VbComputeOptions(
        mode="vbci",
        optimize_orbitals=False,
        active_pi_atoms=case["active_pi_atoms"],
        active_electron_count=case["active_electron_count"],
        active_spin=case["active_spin"],
        include_ionic=True,
        freeze_inactive_orbitals=True,
    )
    return VbDriver().compute(molecule, basis, options=options)


def _run_pi_case_with_mode(case, mode, basis_name="sto-3g", orbital_amplitude_bound=None,
                           orbital_relaxation_symmetry=None):
    molecule = Molecule.read_str(case["xyz"])
    molecule.set_charge(case["charge"])
    molecule.set_multiplicity(case["multiplicity"])
    basis = MolecularBasis.read(molecule, basis_name, ostream=None)
    options = VbComputeOptions(
        mode=mode,
        optimize_orbitals=(mode == "vbscf"),
        include_bovb=(mode == "bovb"),
        active_pi_atoms=case["active_pi_atoms"],
        active_electron_count=case["active_electron_count"],
        active_spin=case["active_spin"],
        orbital_amplitude_bound=orbital_amplitude_bound,
        orbital_relaxation_symmetry=orbital_relaxation_symmetry,
        include_ionic=True,
        freeze_inactive_orbitals=True,
    )
    return VbDriver().compute(molecule, basis, options=options)


def _run_allyl_cation_pi(mode="vbci", basis_name="sto-3g", orbital_amplitude_bound=None,
                         orbital_relaxation_symmetry=None):
    case = ALLYL_CATION_CASE
    molecule = Molecule.read_str(case["xyz"])
    molecule.set_charge(case["charge"])
    molecule.set_multiplicity(case["multiplicity"])
    basis = MolecularBasis.read(molecule, basis_name, ostream=None)
    options = VbComputeOptions(
        mode=mode,
        optimize_orbitals=(mode == "vbscf"),
        active_pi_atoms=case["active_pi_atoms"],
        active_electron_count=case["active_electron_count"],
        active_spin=case["active_spin"],
        orbital_amplitude_bound=orbital_amplitude_bound,
        orbital_relaxation_symmetry=orbital_relaxation_symmetry,
        include_ionic=True,
        freeze_inactive_orbitals=True,
    )
    return VbDriver().compute(molecule, basis, options=options)


def test_h2_stretched_vbci_vbscf_are_benchmarked_against_rhf_uhf():
    molecule, basis, vbci = _run_h2_vb(4.0, mode="vbci")
    _, _, vbscf = _run_h2_vb(4.0, mode="vbscf")
    _, _, bovb = _run_h2_vb(4.0, mode="bovb", include_bovb=True)
    rhf_energy = _run_h2_hf_reference(molecule, basis, unrestricted=False)
    uhf_energy = _run_h2_hf_reference(molecule, basis, unrestricted=True)

    assert np.isfinite(rhf_energy)
    assert np.isfinite(uhf_energy)
    assert np.isfinite(vbci["energy"])
    assert np.isfinite(vbscf["energy"])
    assert np.isfinite(bovb["energy"])

    assert uhf_energy < rhf_energy - 0.1
    assert vbci["energy"] < rhf_energy - 0.1
    assert abs(vbci["energy"] - uhf_energy) < 1.0e-3
    assert vbscf["energy"] <= vbci["energy"] + 1.0e-8
    assert bovb["energy"] <= vbscf["energy"] + 1.0e-8

    diagnostics = vbci["diagnostics"]
    assert diagnostics["active_space_model"] == "one-active-bond"
    assert diagnostics["retained_overlap_rank"] == 3
    assert diagnostics["generated_structure_labels"] == [
        "covalent",
        "ionic_A_minus_B_plus",
        "ionic_A_plus_B_minus",
    ]
    assert diagnostics["localized_template_labels"] == diagnostics[
        "generated_structure_labels"
    ]
    assert np.isfinite(diagnostics["best_localized_template_energy"])
    assert np.isfinite(diagnostics["resonance_energy"])
    assert diagnostics["resonance_energy"] >= -1.0e-8
    assert np.isclose(np.sum(vbci["weights"]), 1.0, atol=1.0e-8)
    assert np.isclose(np.sum(vbci["lowdin_weights"]), 1.0, atol=1.0e-8)

    vbscf_diagnostics = vbscf["diagnostics"]
    assert vbscf_diagnostics["message"].startswith("H2 VB-SCF result")
    assert vbscf_diagnostics["retained_overlap_rank"] == 3
    assert np.isclose(np.sum(vbscf["lowdin_weights"]), 1.0, atol=1.0e-8)

    bovb_diagnostics = bovb["diagnostics"]
    assert bovb_diagnostics["message"].startswith("H2 BOVB result")
    assert bovb_diagnostics["structure_specific_orbitals"] is True
    assert bovb_diagnostics["retained_overlap_rank"] == 3
    assert np.isclose(np.sum(bovb["lowdin_weights"]), 1.0, atol=1.0e-8)


def test_h2_bovb_uses_center_local_breathing_space_in_split_valence_basis():
    _, _, bovb = _run_h2_vb(
        1.4,
        mode="bovb",
        include_bovb=True,
        basis_name="6-31g",
    )
    diagnostics = bovb["diagnostics"]

    assert np.isfinite(bovb["energy"])
    assert diagnostics["bovb_has_external_breathing_space"] is True
    assert bovb["energy"] < diagnostics["bovb_initial_energy"] - 1.0e-4
    assert abs(diagnostics["bovb_covalent_breathing"]) > 1.0e-4
    assert abs(diagnostics["bovb_ionic_breathing"]) > 1.0e-4
    assert diagnostics["bovb_model"] == (
        "h2-center-local-structure-specific-breathing-orbitals"
    )


def test_h2_vbscf_uses_common_breathing_space_in_split_valence_basis():
    _, _, vbscf = _run_h2_vb(
        0.74,
        mode="vbscf",
        basis_name="6-31g",
    )
    diagnostics = vbscf["diagnostics"]

    assert np.isfinite(vbscf["energy"])
    assert diagnostics["vbscf_model"] == "h2-common-center-local-breathing-orbitals"
    assert diagnostics["vbscf_has_external_breathing_space"] is True
    assert diagnostics["vbscf_used_fixed_orbital_limit"] is False
    assert vbscf["energy"] < diagnostics["vbscf_initial_energy"] - 1.0e-4
    assert abs(diagnostics["vbscf_breathing"]) > 1.0e-4
    assert np.isclose(np.sum(vbscf["weights"]), 1.0, atol=1.0e-8)
    assert np.isclose(np.sum(vbscf["lowdin_weights"]), 1.0, atol=1.0e-8)


def test_h2_split_valence_scan_keeps_stable_atom_centered_active_space():
    labels = []
    for distance in (0.3, 0.8, 1.4, 2.5, 4.0):
        _, _, result = _run_h2_vb(distance, mode="vbci", basis_name="6-31g")
        labels.append(result["active_space"].active_candidate_label)
        diagnostics = result["diagnostics"]
        assert diagnostics["h2_stable_atom_centered_active_space"] is True
        assert diagnostics["active_orbital_labels"] == [
            "active_A_atom_1",
            "active_B_atom_2",
        ]

    assert labels == ["H2_atom_centered_sigma"] * len(labels)


def test_h2_split_valence_vb_methods_reach_uhf_asymptote():
    molecule, basis, vbci = _run_h2_vb(5.0, mode="vbci", basis_name="6-31g")
    _, _, vbscf = _run_h2_vb(5.0, mode="vbscf", basis_name="6-31g")
    _, _, bovb = _run_h2_vb(
        5.0,
        mode="bovb",
        include_bovb=True,
        basis_name="6-31g",
    )
    uhf_energy = _run_h2_hf_reference(molecule, basis, unrestricted=True)

    for result in (vbci, vbscf, bovb):
        assert abs(float(result["energy"]) - uhf_energy) * HARTREE_TO_KJMOL < 10.0


def test_h2_dissociation_curve_compares_rhf_uhf_vb_methods():
    distances = (0.74, 4.0)
    curves = {key: [] for key in ("rhf", "uhf", "vbci", "vbscf", "bovb")}
    labels = []

    for distance in distances:
        molecule, basis, vbci = _run_h2_vb(distance, mode="vbci")
        _, _, vbscf = _run_h2_vb(distance, mode="vbscf")
        _, _, bovb = _run_h2_vb(distance, mode="bovb", include_bovb=True)
        curves["rhf"].append(_run_hf_reference(molecule, basis, unrestricted=False))
        curves["uhf"].append(_run_hf_reference(molecule, basis, unrestricted=True))
        curves["vbci"].append(float(vbci["energy"]))
        curves["vbscf"].append(float(vbscf["energy"]))
        curves["bovb"].append(float(bovb["energy"]))
        labels.append(vbci["active_space"].active_candidate_label)

    relative_curves = {
        method: _relative_potential_curve(energies)
        for method, energies in curves.items()
    }

    assert labels == ["H2_atom_centered_sigma"] * len(distances)
    assert relative_curves["rhf"][-1] == pytest.approx(0.0, abs=1.0e-10)
    assert relative_curves["uhf"][-1] == pytest.approx(0.0, abs=1.0e-10)
    assert relative_curves["vbci"][-1] == pytest.approx(0.0, abs=1.0e-10)
    assert relative_curves["rhf"][0] < -100.0
    assert relative_curves["uhf"][0] < -100.0
    assert relative_curves["vbci"][0] < -100.0
    assert curves["uhf"][-1] < curves["rhf"][-1] - 0.1
    assert curves["vbscf"][-1] <= curves["vbci"][-1] + 1.0e-8
    assert curves["bovb"][-1] <= curves["vbscf"][-1] + 1.0e-8


def test_ethane_cc_dissociation_keeps_requested_sigma_active_bond():
    distances = (1.54, 4.0)
    anchor_molecule = _ethane_molecule(distances[0])
    anchor_basis = MolecularBasis.read(anchor_molecule, "sto-3g", ostream=None)
    anchor = _run_one_bond_vb(
        anchor_molecule,
        anchor_basis,
        mode="vbci",
        active_subtype="sigma",
    )
    reference_orbitals = tuple(anchor["active_space"].active_orbitals)
    curves = {key: [] for key in ("rhf", "uhf", "vbci", "vbscf", "bovb")}
    fallback_seen = False

    for distance in distances:
        molecule = _ethane_molecule(distance)
        basis = MolecularBasis.read(molecule, "sto-3g", ostream=None)
        curves["rhf"].append(_run_hf_reference(molecule, basis, unrestricted=False))
        curves["uhf"].append(_run_hf_reference(molecule, basis, unrestricted=True))
        for mode in ("vbci", "vbscf", "bovb"):
            result = _run_one_bond_vb(
                molecule,
                basis,
                mode=mode,
                active_subtype="sigma",
                reference_orbitals=reference_orbitals,
            )
            diagnostics = result["diagnostics"]
            curves[mode].append(float(result["energy"]))
            assert diagnostics["active_space_model"] == "one-active-bond"
            assert diagnostics["active_candidate_subtype"] == "sigma"
            assert diagnostics["active_bond_state_tracked"] is True
            assert diagnostics["active_orbitals_orthogonalized_to_frozen_space"] is True
            if distance >= 3.0:
                assert diagnostics["active_uhf_frontier_orbitals"] is True
                assert diagnostics["active_stretched_orbital_source"] == "uhf_frontier"
            assert np.isclose(np.sum(result["lowdin_weights"]), 1.0, atol=1.0e-8)
            fallback_seen = fallback_seen or diagnostics["active_candidate_fallback"]

    for energies in curves.values():
        assert np.all(np.isfinite(energies))

    assert fallback_seen is True
    assert curves["uhf"][-1] < curves["rhf"][-1] - 0.1
    assert _relative_potential_curve(curves["uhf"])[0] < -100.0
    assert np.isfinite(_relative_potential_curve(curves["vbci"])[0])
    assert _relative_potential_curve(curves["vbci"])[-1] == pytest.approx(0.0, abs=1.0e-10)
    for mode in ("vbci", "vbscf", "bovb"):
        asymptote_error = (curves[mode][-1] - curves["uhf"][-1]) * HARTREE_TO_KJMOL
        assert asymptote_error > -10.0
        assert asymptote_error < 200.0
    assert curves["vbscf"][-1] <= curves["vbci"][-1] + 1.0e-8
    assert curves["bovb"][-1] <= curves["vbscf"][-1] + 1.0e-8


def test_ethane_split_valence_sigma_vb_methods_reach_uhf_asymptote():
    anchor_molecule = _ethane_molecule(1.54)
    anchor_basis = MolecularBasis.read(anchor_molecule, "6-31g", ostream=None)
    anchor = _run_one_bond_vb(
        anchor_molecule,
        anchor_basis,
        mode="vbci",
        active_subtype="sigma",
    )
    reference_orbitals = tuple(anchor["active_space"].active_orbitals)
    molecule = _ethane_molecule(5.0)
    basis = MolecularBasis.read(molecule, "6-31g", ostream=None)
    uhf_energy = _run_hf_reference(molecule, basis, unrestricted=True)

    for mode in ("vbci", "vbscf", "bovb"):
        result = _run_one_bond_vb(
            molecule,
            basis,
            mode=mode,
            active_subtype="sigma",
            reference_orbitals=reference_orbitals,
        )
        diagnostics = result["diagnostics"]
        asymptote_error = (float(result["energy"]) - uhf_energy) * HARTREE_TO_KJMOL
        assert diagnostics["active_uhf_frontier_orbitals"] is True
        assert diagnostics["active_stretched_orbital_source"] == "uhf_frontier"
        if mode == "bovb":
            assert diagnostics["bovb_conservative_non_h2_limit"] is True
            assert diagnostics["bovb_used_fixed_orbital_limit"] is True
        assert -10.0 <= asymptote_error <= 120.0


def test_ethane_split_valence_sigma_curve_has_no_switch_spike():
    anchor_molecule = _ethane_molecule(1.54)
    anchor_basis = MolecularBasis.read(anchor_molecule, "6-31g", ostream=None)
    anchor = _run_one_bond_vb(
        anchor_molecule,
        anchor_basis,
        mode="vbci",
        active_subtype="sigma",
    )
    reference_orbitals = tuple(anchor["active_space"].active_orbitals)
    distances = (2.04, 2.47, 2.89, 3.31, 5.0)
    curves = {key: [] for key in ("uhf", "vbci", "bovb")}
    stretched_sources = []
    bovb_fixed_limits = []

    for distance in distances:
        molecule = _ethane_molecule(distance)
        basis = MolecularBasis.read(molecule, "6-31g", ostream=None)
        curves["uhf"].append(_run_hf_reference(molecule, basis, unrestricted=True))
        result = _run_one_bond_vb(
            molecule,
            basis,
            mode="vbci",
            active_subtype="sigma",
            reference_orbitals=reference_orbitals,
        )
        curves["vbci"].append(float(result["energy"]))
        stretched_sources.append(result["diagnostics"].get("active_stretched_orbital_source"))
        bovb = _run_one_bond_vb(
            molecule,
            basis,
            mode="bovb",
            active_subtype="sigma",
            reference_orbitals=reference_orbitals,
        )
        curves["bovb"].append(float(bovb["energy"]))
        bovb_fixed_limits.append(bovb["diagnostics"].get("bovb_used_fixed_orbital_limit"))

    relative_vbci = _relative_potential_curve(curves["vbci"])
    relative_bovb = _relative_potential_curve(curves["bovb"])
    relative_uhf = _relative_potential_curve(curves["uhf"])
    vbci_jumps = np.abs(np.diff(relative_vbci))
    bovb_jumps = np.abs(np.diff(relative_bovb))

    assert stretched_sources == ["uhf_frontier"] * len(distances)
    assert bovb_fixed_limits == [True] * len(distances)
    assert np.max(vbci_jumps) < 175.0
    assert np.max(bovb_jumps) < 175.0
    assert np.allclose(relative_bovb, relative_vbci, atol=1.0e-8)
    assert relative_vbci[-1] == pytest.approx(0.0, abs=1.0e-10)
    assert relative_uhf[-1] == pytest.approx(0.0, abs=1.0e-10)
    assert relative_vbci[0] < relative_vbci[-1] - 100.0


def test_ethylene_cc_dissociation_keeps_requested_pi_active_bond():
    distances = (1.339, 4.0)
    anchor_molecule = _ethylene_molecule(distances[0])
    anchor_basis = MolecularBasis.read(anchor_molecule, "sto-3g", ostream=None)
    anchor = _run_one_bond_vb(
        anchor_molecule,
        anchor_basis,
        mode="vbci",
        active_subtype="pi",
    )
    reference_orbitals = tuple(anchor["active_space"].active_orbitals)
    curves = {key: [] for key in ("rhf", "uhf", "vbci", "vbscf", "bovb")}
    fallback_seen = False

    for distance in distances:
        molecule = _ethylene_molecule(distance)
        basis = MolecularBasis.read(molecule, "sto-3g", ostream=None)
        curves["rhf"].append(_run_hf_reference(molecule, basis, unrestricted=False))
        curves["uhf"].append(_run_hf_reference(molecule, basis, unrestricted=True))
        for mode in ("vbci", "vbscf", "bovb"):
            result = _run_one_bond_vb(
                molecule,
                basis,
                mode=mode,
                active_subtype="pi",
                reference_orbitals=reference_orbitals,
            )
            diagnostics = result["diagnostics"]
            curves[mode].append(float(result["energy"]))
            assert diagnostics["active_space_model"] == "one-active-pi-bond"
            assert diagnostics["active_candidate_subtype"] == "pi"
            assert diagnostics["active_bond_state_tracked"] is True
            assert diagnostics["active_orbitals_orthogonalized_to_frozen_space"] is True
            assert np.isclose(np.sum(result["lowdin_weights"]), 1.0, atol=1.0e-8)
            fallback_seen = fallback_seen or diagnostics["active_candidate_fallback"]

    for energies in curves.values():
        assert np.all(np.isfinite(energies))

    assert fallback_seen is True
    assert curves["uhf"][-1] < curves["rhf"][-1] - 0.1
    assert _relative_potential_curve(curves["uhf"])[0] < -100.0
    assert np.isfinite(_relative_potential_curve(curves["vbci"])[0])
    assert _relative_potential_curve(curves["vbci"])[-1] == pytest.approx(0.0, abs=1.0e-10)
    # This is a state-tracking/finite-value regression only. A physical ethylene
    # C-C dissociation curve needs a four-electron sigma+pi active space.
    assert curves["vbscf"][-1] <= curves["vbci"][-1] + 1.0e-8
    assert curves["bovb"][-1] <= curves["vbscf"][-1] + 1.0e-8


def test_ethylene_two_orbital_pi_bovb_uses_frozen_sigma_embedding():
    _, _, bovb = _run_ethylene_pi_vb(mode="bovb", include_bovb=True)
    diagnostics = bovb["diagnostics"]

    assert np.isfinite(bovb["energy"])
    assert np.isfinite(diagnostics["bovb_initial_energy"])
    assert bovb["energy"] <= diagnostics["bovb_initial_energy"] + 1.0e-8

    assert diagnostics["message"].startswith("Two-orbital BOVB result")
    assert diagnostics["bovb_model"] == (
        "two-orbital-center-local-structure-specific-breathing-orbitals"
    )
    assert diagnostics["active_space_model"] == "one-active-pi-bond"
    assert diagnostics["active_candidate_subtype"] == "pi"
    assert diagnostics["structure_specific_orbitals"] is True
    assert diagnostics["bovb_conservative_non_h2_limit"] is True
    assert diagnostics["bovb_used_fixed_orbital_limit"] is True
    assert diagnostics["frozen_hf_embedding"] is True
    assert np.isclose(diagnostics["active_reference_electron_count"], 2.0, atol=1.0e-8)
    assert diagnostics["frozen_electron_count"] > 10.0
    assert diagnostics["retained_overlap_rank"] == 3
    assert np.isclose(np.sum(bovb["lowdin_weights"]), 1.0, atol=1.0e-8)


def test_bovb_request_is_restricted_to_two_orbital_active_spaces():
    case = {
        "xyz": """
C -1.3000  0.0000  0.0000
C  0.0000  0.0000  0.0000
C  1.3000  0.0000  0.0000
H -1.8500  0.9200  0.0000
H -1.8500 -0.9200  0.0000
H  0.0000  1.0800  0.0000
H  1.8500  0.9200  0.0000
H  1.8500 -0.9200  0.0000
""",
        "charge": 1,
        "multiplicity": 1,
    }
    molecule = Molecule.read_str(case["xyz"])
    molecule.set_charge(case["charge"])
    molecule.set_multiplicity(case["multiplicity"])
    basis = MolecularBasis.read(molecule, "sto-3g", ostream=None)

    with pytest.raises(NotImplementedError, match="two-orbital two-electron"):
        VbDriver().compute(
            molecule,
            basis,
            options=VbComputeOptions(
                mode="bovb",
                include_bovb=True,
                active_pi_atoms=(0, 1, 2),
                active_electron_count=2,
                active_spin="singlet",
            ),
        )


def test_allyl_cation_compact_csf_hamiltonian_tracks_full_vbci_reference():
    compact = _run_allyl_cation_pi(mode="compact-csf")

    compact_diag = compact["diagnostics"]
    reference_energy = compact_diag["compact_csf_full_reference_energy"]
    energy_error = compact_diag["compact_csf_energy_error_to_full_reference"]

    assert np.isfinite(compact["energy"])
    assert np.isfinite(reference_energy)
    assert np.isfinite(energy_error)
    assert compact["energy"] >= reference_energy - 1.0e-8
    assert energy_error >= -1.0e-8
    assert np.isclose(compact["energy"] - reference_energy, energy_error, atol=1.0e-8)

    assert compact_diag["active_space_model"] == "fixed-orbital-multicenter-pi"
    assert compact_diag["vb_method"] == "compact-csf"
    assert compact_diag["compact_csf_model"] == (
        "two-electron-singlet-pi-adjacent-bond-templates"
    )
    assert compact_diag["compact_csf_count"] == 2
    assert compact_diag["compact_csf_types"] == ["allyl_cation", "allyl_cation"]
    assert compact_diag["localized_template_labels"] == compact_diag[
        "compact_csf_labels"
    ]
    assert np.isfinite(compact_diag["best_localized_template_energy"])
    assert np.isfinite(compact_diag["resonance_energy"])
    assert compact_diag["resonance_energy"] >= -1.0e-8
    assert compact_diag["compact_csf_retained_rank"] == 2
    assert compact_diag["compact_csf_captured_subspace_weight"] > 0.45
    assert np.isclose(np.sum(compact["lowdin_weights"]), 1.0, atol=1.0e-8)
    assert np.all(np.asarray(compact["lowdin_weights"]) > 0.1)
    assert compact_diag["frozen_hf_embedding"] is True
    assert compact_diag["active_orbitals_orthogonalized_to_frozen_space"] is True


def test_allyl_cation_compact_csf_bovb_lowers_split_valence_compact_limit():
    bovb = _run_allyl_cation_pi(mode="compact-csf-bovb", basis_name="6-31g")

    diagnostics = bovb["diagnostics"]

    assert np.isfinite(bovb["energy"])
    assert diagnostics["vb_method"] == "compact-csf-bovb"
    assert diagnostics["compact_csf_bovb_model"] == (
        "allyl-cation-two-electron-pi-center-local-breathing"
    )
    assert diagnostics["compact_csf_bovb_has_external_breathing_space"] is True
    assert diagnostics["active_orbitals_orthogonalized_to_frozen_space"] is True
    assert diagnostics["localized_template_labels"] == diagnostics[
        "compact_csf_labels"
    ]
    assert np.isfinite(diagnostics["best_localized_template_energy"])
    assert np.isfinite(diagnostics["resonance_energy"])
    assert diagnostics["resonance_energy"] >= -1.0e-8
    assert np.isfinite(diagnostics["compact_csf_bovb_initial_energy"])
    assert bovb["energy"] < diagnostics["compact_csf_bovb_initial_energy"] - 1.0e-4
    assert diagnostics["compact_csf_bovb_energy_lowering"] > 1.0e-4
    assert diagnostics["compact_csf_bovb_used_fixed_orbital_limit"] is False
    assert len(diagnostics["compact_csf_bovb_breathing_amplitudes"]) == 2
    assert any(
        abs(value) > 1.0e-4
        for value in diagnostics["compact_csf_bovb_breathing_amplitudes"]
    )
    assert np.isclose(np.sum(bovb["lowdin_weights"]), 1.0, atol=1.0e-8)


@pytest.mark.timeconsuming
def test_allyl_cation_vbscf_uses_common_pi_breathing_space():
    vbci = _run_allyl_cation_pi(mode="vbci", basis_name="6-31g")
    vbscf = _run_allyl_cation_pi(mode="vbscf", basis_name="6-31g")
    diagnostics = vbscf["diagnostics"]

    assert np.isfinite(vbci["energy"])
    assert np.isfinite(vbscf["energy"])
    assert diagnostics["organic_pi_vbscf_model"] == (
        "common-center-local-pi-breathing-orbitals"
    )
    assert diagnostics["organic_pi_vbscf_has_external_relaxation_space"] is True
    assert diagnostics["organic_pi_vbscf_used_fixed_orbital_limit"] is False
    assert diagnostics["organic_pi_vbscf_energy_lowering"] >= -1.0e-8
    assert vbscf["energy"] <= diagnostics["organic_pi_vbscf_initial_energy"] + 1.0e-8
    assert np.isclose(np.sum(vbscf["lowdin_weights"]), 1.0, atol=1.0e-8)


def test_organic_pi_vbscf_reports_configured_amplitude_bound():
    vbscf = _run_allyl_cation_pi(
        mode="vbscf",
        basis_name="6-31g",
        orbital_amplitude_bound=0.05,
    )
    diagnostics = vbscf["diagnostics"]
    max_abs_amplitude = diagnostics["organic_pi_vbscf_max_abs_orbital_amplitude"]

    assert diagnostics["organic_pi_vbscf_orbital_amplitude_bound"] == pytest.approx(0.05)
    assert max_abs_amplitude <= 0.05 + 1.0e-8
    assert diagnostics["organic_pi_vbscf_hit_amplitude_bound"] is (max_abs_amplitude >= 0.05 - 1.0e-8)


def test_organic_pi_vbscf_equivalent_center_amplitude_mode():
    vbscf = _run_allyl_cation_pi(
        mode="vbscf",
        basis_name="6-31g",
        orbital_amplitude_bound=0.05,
        orbital_relaxation_symmetry="equivalent-centers",
    )
    diagnostics = vbscf["diagnostics"]
    amplitudes = np.asarray(diagnostics["organic_pi_vbscf_orbital_amplitudes"])

    assert diagnostics["organic_pi_vbscf_relaxation_symmetry"] == "equivalent-centers"
    assert diagnostics["organic_pi_vbscf_equivalent_center_amplitude"] is True
    assert len(diagnostics["organic_pi_vbscf_optimizer_parameters"]) == 1
    assert amplitudes.size == 3
    assert np.allclose(amplitudes, amplitudes[0])
    assert np.max(np.abs(amplitudes)) <= 0.05 + 1.0e-8


@pytest.mark.timeconsuming
@pytest.mark.parametrize(
    "case",
    [item for item in PI_CASES if item["name"] in {"allyl_radical", "allyl_anion"}],
    ids=[
        item["name"]
        for item in PI_CASES
        if item["name"] in {"allyl_radical", "allyl_anion"}
    ],
)
def test_allyl_radical_anion_generalized_method_routes(case):
    vbci = _run_pi_case_with_mode(case, "vbci", basis_name="6-31g")
    vbscf = _run_pi_case_with_mode(case, "vbscf", basis_name="6-31g")
    bovb = _run_pi_case_with_mode(case, "bovb", basis_name="6-31g")
    compact = _run_pi_case_with_mode(case, "compact-csf", basis_name="6-31g")

    assert np.isfinite(vbci["energy"])
    assert np.isfinite(vbscf["energy"])
    assert np.isfinite(bovb["energy"])
    assert np.isfinite(compact["energy"])

    vbscf_diag = vbscf["diagnostics"]
    bovb_diag = bovb["diagnostics"]
    compact_diag = compact["diagnostics"]

    assert vbscf_diag["organic_pi_vbscf_model"] == (
        "common-center-local-pi-breathing-orbitals"
    )
    assert vbscf_diag["organic_pi_vbscf_has_external_relaxation_space"] is True
    assert vbscf_diag["organic_pi_vbscf_used_fixed_orbital_limit"] is False
    assert vbscf_diag["organic_pi_vbscf_energy_lowering"] >= -1.0e-8
    assert vbscf["energy"] <= vbscf_diag["organic_pi_vbscf_initial_energy"] + 1.0e-8
    assert bovb_diag["organic_pi_bovb_model"] == (
        "determinant-ci-fixed-orbital-zero-amplitude-limit"
    )
    assert bovb_diag["bovb_used_fixed_orbital_limit"] is True
    assert bovb_diag["organic_pi_bovb_energy_lowering"] == 0.0
    assert np.isclose(
        bovb["energy"],
        bovb_diag["organic_pi_bovb_initial_energy"],
        atol=1.0e-10,
    )

    assert compact_diag["compact_csf_model"] == (
        "graph-template-determinant-ci-subspace"
    )
    assert compact_diag["compact_csf_count"] == case["expected_templates"]
    assert np.isfinite(compact_diag["compact_csf_full_reference_energy"])
    assert np.isfinite(compact_diag["compact_csf_energy_error_to_full_reference"])
    assert compact_diag["compact_csf_energy_error_to_full_reference"] >= -1.0e-8
    assert np.isclose(
        compact["energy"] - compact_diag["compact_csf_full_reference_energy"],
        compact_diag["compact_csf_energy_error_to_full_reference"],
        atol=1.0e-8,
    )
    assert compact_diag["compact_csf_captured_subspace_weight"] >= 0.0
    assert compact_diag["compact_csf_captured_subspace_weight"] <= 1.0 + 1.0e-8
    assert np.isfinite(compact_diag["best_localized_template_energy"])
    assert np.isfinite(compact_diag["resonance_energy"])
    assert compact_diag["resonance_energy"] >= -1.0e-8
    assert np.isclose(np.sum(compact["lowdin_weights"]), 1.0, atol=1.0e-8)


@pytest.mark.timeconsuming
def test_benzene_compact_csf_projection_tracks_full_reference():
    case = next(item for item in PI_CASES if item["name"] == "benzene")
    vbci = _run_pi_case_with_mode(case, "vbci", basis_name="sto-3g")
    compact = _run_pi_case_with_mode(case, "compact-csf", basis_name="sto-3g")
    compact_diag = compact["diagnostics"]

    assert np.isfinite(vbci["energy"])
    assert np.isfinite(compact["energy"])
    assert compact_diag["compact_csf_model"] == (
        "graph-template-determinant-ci-subspace"
    )
    assert compact_diag["compact_csf_count"] == 15
    assert compact_diag["compact_csf_types"].count("benzene_kekule") == 2
    assert compact_diag["compact_csf_types"].count("benzene_dewar") == 13
    assert np.isfinite(compact_diag["compact_csf_full_reference_energy"])
    assert np.isfinite(compact_diag["compact_csf_energy_error_to_full_reference"])
    assert compact_diag["compact_csf_energy_error_to_full_reference"] >= -1.0e-8
    assert np.isclose(
        compact["energy"] - compact_diag["compact_csf_full_reference_energy"],
        compact_diag["compact_csf_energy_error_to_full_reference"],
        atol=1.0e-8,
    )
    assert compact_diag["compact_csf_captured_subspace_weight"] >= 0.0
    assert compact_diag["compact_csf_captured_subspace_weight"] <= 1.0 + 1.0e-8
    assert np.isfinite(compact_diag["best_localized_template_energy"])
    assert np.isfinite(compact_diag["resonance_energy"])
    assert compact_diag["resonance_energy"] >= -1.0e-8
    assert np.isclose(np.sum(compact["lowdin_weights"]), 1.0, atol=1.0e-8)


@pytest.mark.timeconsuming
def test_benzene_vbscf_defaults_to_equivalent_center_stabilization():
    case = next(item for item in PI_CASES if item["name"] == "benzene")
    vbscf = _run_pi_case_with_mode(case, "vbscf", basis_name="sto-3g")
    diagnostics = vbscf["diagnostics"]
    amplitudes = np.asarray(diagnostics["organic_pi_vbscf_orbital_amplitudes"])

    assert np.isfinite(vbscf["energy"])
    assert diagnostics["organic_pi_vbscf_model"] == (
        "common-center-local-pi-breathing-orbitals"
    )
    assert diagnostics["organic_pi_vbscf_auto_equivalent_center_relaxation"] is True
    assert diagnostics["organic_pi_vbscf_relaxation_symmetry"] == "equivalent-centers"
    assert diagnostics["organic_pi_vbscf_equivalent_center_amplitude"] is True
    assert len(diagnostics["organic_pi_vbscf_optimizer_parameters"]) == 1
    assert amplitudes.size == 6
    assert np.allclose(amplitudes, amplitudes[0])
    assert diagnostics["organic_pi_vbscf_orbital_amplitude_bound"] == pytest.approx(0.05)
    assert np.max(np.abs(amplitudes)) <= 0.05 + 1.0e-8
    assert diagnostics["organic_pi_vbscf_energy_lowering"] >= -1.0e-8
    assert vbscf["energy"] <= diagnostics["organic_pi_vbscf_initial_energy"] + 1.0e-8
    assert np.isclose(np.sum(vbscf["lowdin_weights"]), 1.0, atol=1.0e-8)


@pytest.mark.timeconsuming
@pytest.mark.parametrize("case", PI_CASES, ids=[case["name"] for case in PI_CASES])
def test_fixed_orbital_pi_chemical_resonance_diagnostics(case):
    result = _run_pi_case(case)
    diagnostics = result["diagnostics"]

    weights = np.asarray(result["weights"], dtype=float)
    chemical_weights = np.asarray(
        diagnostics["chemical_resonance_weights"], dtype=float
    )
    projection_weights = np.asarray(
        diagnostics["chemical_resonance_projection_weights"], dtype=float
    )
    chemical_types = set(diagnostics["chemical_resonance_types"])

    assert diagnostics["active_space_model"] == "fixed-orbital-multicenter-pi"
    assert diagnostics["determinant_count"] == case["expected_determinants"]
    assert diagnostics["chemical_resonance_count"] == case["expected_templates"]
    assert diagnostics["chemical_resonance_retained_rank"] > 0
    assert case["expected_types"].issubset(chemical_types)
    assert diagnostics["localized_template_energy_model"].startswith(
        "single-basis-vector Rayleigh quotient"
    )
    assert len(diagnostics["localized_template_energies"]) == case[
        "expected_determinants"
    ]
    assert np.isfinite(diagnostics["best_localized_template_energy"])
    assert np.isfinite(diagnostics["resonance_energy"])
    assert diagnostics["resonance_energy"] >= -1.0e-8

    assert np.isfinite(result["energy"])
    assert np.all(np.isfinite(weights))
    assert np.isclose(np.sum(weights), 1.0, atol=1.0e-8)
    assert np.all(np.isfinite(chemical_weights))
    assert np.all(np.isfinite(projection_weights))
    assert np.all(chemical_weights >= -1.0e-12)
    assert np.all(projection_weights >= -1.0e-12)
    assert np.isclose(np.sum(chemical_weights), 1.0, atol=1.0e-8)
    assert diagnostics["chemical_resonance_subspace_weight"] >= 0.0
    assert diagnostics["chemical_resonance_subspace_weight"] <= 1.0 + 1.0e-8

    if case["active_spin"] == "singlet":
        assert diagnostics["determinant_ci_root_selection"] == (
            "lowest alpha/beta exchange-symmetric singlet root"
        )
        assert diagnostics["determinant_ci_spin_exchange_parity"] > 0.5

    if "min_capture" in case:
        assert diagnostics["chemical_resonance_subspace_weight"] > case["min_capture"]


def test_butadiene_singlet_root_and_template_phase_regression():
    case = next(item for item in PI_CASES if item["name"] == "butadiene")
    result = _run_pi_case(case)
    diagnostics = result["diagnostics"]

    assert diagnostics["determinant_ci_spin_exchange_parity"] > 0.5
    assert diagnostics["chemical_resonance_subspace_weight"] > 1.0e-3
    assert np.isclose(
        np.sum(diagnostics["chemical_resonance_weights"]),
        1.0,
        atol=1.0e-8,
    )
    assert "butadiene_kekule" in diagnostics["chemical_resonance_types"]
