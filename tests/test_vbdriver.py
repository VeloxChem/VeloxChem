import numpy as np
import pytest

from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.vbdriver import VbComputeOptions, VbDriver


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
