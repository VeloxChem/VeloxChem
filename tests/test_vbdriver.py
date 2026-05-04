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


def _ethylene_molecule():
    molecule = Molecule.read_str("""
C -0.6695  0.0000  0.0000
C  0.6695  0.0000  0.0000
H -1.2321  0.9289  0.0000
H -1.2321 -0.9289  0.0000
H  1.2321  0.9289  0.0000
H  1.2321 -0.9289  0.0000
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


def _run_allyl_cation_pi(mode="vbci", basis_name="sto-3g"):
    case = ALLYL_CATION_CASE
    molecule = Molecule.read_str(case["xyz"])
    molecule.set_charge(case["charge"])
    molecule.set_multiplicity(case["multiplicity"])
    basis = MolecularBasis.read(molecule, basis_name, ostream=None)
    options = VbComputeOptions(
        mode=mode,
        optimize_orbitals=False,
        active_pi_atoms=case["active_pi_atoms"],
        active_electron_count=case["active_electron_count"],
        active_spin=case["active_spin"],
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
    assert compact_diag["compact_csf_retained_rank"] == 2
    assert compact_diag["compact_csf_captured_subspace_weight"] > 0.45
    assert np.isclose(np.sum(compact["lowdin_weights"]), 1.0, atol=1.0e-8)
    assert np.all(np.asarray(compact["lowdin_weights"]) > 0.1)
    assert compact_diag["frozen_hf_embedding"] is True


def test_allyl_cation_compact_csf_bovb_lowers_split_valence_compact_limit():
    bovb = _run_allyl_cation_pi(mode="compact-csf-bovb", basis_name="6-31g")

    diagnostics = bovb["diagnostics"]

    assert np.isfinite(bovb["energy"])
    assert diagnostics["vb_method"] == "compact-csf-bovb"
    assert diagnostics["compact_csf_bovb_model"] == (
        "allyl-cation-two-electron-pi-center-local-breathing"
    )
    assert diagnostics["compact_csf_bovb_has_external_breathing_space"] is True
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
