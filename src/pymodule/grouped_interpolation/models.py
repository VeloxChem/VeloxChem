"""Immutable records used by grouped interpolation construction and runtime."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any


@dataclass(frozen=True)
class PrimitiveRotor:
    """One topology-defined collective terminal rotor motion."""

    rotor_id: str
    axis_atoms: tuple[int, int]
    moving_side_atoms: tuple[int, ...]
    owned_atoms: tuple[int, ...]
    torsion_rows: tuple[int, ...]
    signature_rows: tuple[int, ...]
    signature_row_types: tuple[str, ...]
    signature_row_scales: tuple[float, ...]
    row_orientations: tuple[int, ...]
    symmetry_order: int
    symmetry_operation_ids: tuple[str, ...] = ()
    motion_kind: str = "terminal_methyl"
    equivalent_atom_groups: tuple[tuple[int, ...], ...] = ()
    # Largest Cartesian deviation between the exact atom permutation and a rigid
    # rotation by one period, measured on the anchor geometry.  It is zero for a
    # group whose own internal geometry carries the full cyclic symmetry and
    # grows once the group is internally distorted -- a methoxy methyl, for
    # instance, has one short anti C-H bond and a narrowed O-C-H angle.  A
    # distorted group is still a valid rotor, but its symmetry images no longer
    # lie on the rigid rotation path that the phase grid samples.
    rigid_symmetry_defect_bohr: float = 0.0

    def __post_init__(self) -> None:
        if len(self.axis_atoms) != 2 or self.axis_atoms[0] == self.axis_atoms[1]:
            raise ValueError("A rotor axis must contain two distinct atoms.")
        if self.symmetry_order < 1:
            raise ValueError("Rotor symmetry_order must be at least one.")
        n_signature = len(self.signature_rows)
        if not (
            n_signature == len(self.signature_row_types)
            == len(self.signature_row_scales)
            == len(self.row_orientations)
        ):
            raise ValueError("Rotor signature metadata lengths do not agree.")
        if not set(self.signature_rows).issubset(self.torsion_rows):
            raise ValueError("Rotor signature rows must be torsion rows.")
        if any(scale <= 0.0 for scale in self.signature_row_scales):
            raise ValueError("Rotor signature scales must be positive.")
        if any(sign not in (-1, 1) for sign in self.row_orientations):
            raise ValueError("Rotor row orientations must be -1 or +1.")
        if self.axis_atoms[1] not in self.moving_side_atoms:
            raise ValueError("The moving-side axis atom must be in moving_side_atoms.")
        if self.axis_atoms[0] in self.moving_side_atoms:
            raise ValueError("The stationary axis atom cannot be in moving_side_atoms.")
        if self.motion_kind not in {"terminal_methyl", "group_rotation"}:
            raise ValueError(f"Unsupported rotor motion_kind: {self.motion_kind}")
        equivalent_atoms = {
            atom for group in self.equivalent_atom_groups for atom in group
        }
        if not equivalent_atoms.issubset(set(self.owned_atoms)):
            raise ValueError("Equivalent atoms must be owned by the rotor.")
        if any(len(group) < 2 for group in self.equivalent_atom_groups):
            raise ValueError("Equivalent-atom groups must contain at least two atoms.")


@dataclass(frozen=True)
class LocalSubset:
    """Persisted row and atom ownership for a local factor or interaction."""

    subset_id: str
    canonical_subset_key: tuple[str, ...]
    role: str
    parent_subset_ids: tuple[str, ...]
    active_atoms: tuple[int, ...]
    environment_atoms: tuple[int, ...]
    relaxation_atoms: tuple[int, ...]
    active_rows: tuple[int, ...]
    response_rows: tuple[int, ...]
    projector_rows: tuple[int, ...]
    relaxation_policy_id: str
    projector_policy_id: str
    anchor_policy_id: str

    def __post_init__(self) -> None:
        if self.role not in {"factor", "overlap", "interaction"}:
            raise ValueError(f"Unsupported grouped subset role: {self.role}")
        if tuple(sorted(self.canonical_subset_key)) != self.canonical_subset_key:
            raise ValueError("canonical_subset_key must be sorted.")


@dataclass(frozen=True)
class GroupedCalculationRequest:
    """A deterministic physical or virtual grouped calculation request."""

    request_id: str
    subset_id: str
    purpose: str
    phase_signature: tuple[float, ...]
    source_geometry_id: str
    coordinate_fingerprint: str
    relaxation_policy_id: str
    projector_policy_id: str
    anchor_id: str
    electronic_state_id: str
    method_id: str
    method_version: str
    is_anchor: bool
    is_virtual_image: bool
    source_request_id: str | None
    symmetry_operation_id: str | None

    def __post_init__(self) -> None:
        if self.is_virtual_image and not self.source_request_id:
            raise ValueError("A virtual request must identify its physical source.")
        if self.is_virtual_image and not self.symmetry_operation_id:
            raise ValueError("A virtual request must identify its symmetry operation.")


@dataclass(frozen=True)
class GroupedCalculationResult:
    """Physical-derivative result at the grouped provider boundary."""

    request_id: str
    converged: bool
    cartesian_coordinates: tuple[tuple[float, float, float], ...]
    energy: float
    cartesian_gradient: tuple[tuple[float, float, float], ...]
    cartesian_hessian: tuple[tuple[float, ...], ...] | None
    internal_coordinates: tuple[float, ...]
    internal_gradient: tuple[float, ...]
    internal_hessian: tuple[tuple[float, ...], ...] | None
    measured_phase_signature: tuple[float, ...]
    constraint_errors: tuple[float, ...]
    coordinate_fingerprint: str
    source_point_label: str
    provenance: str


@dataclass(frozen=True)
class CoordinateDefinition:
    """Canonical, serializable definition of the shared internal chart."""

    rows: tuple[tuple[int, ...], ...]
    row_types: tuple[str, ...]
    units: tuple[str, ...]
    periodicities: tuple[float | None, ...]
    torsion_convention: str
    improper_convention: str
    bond_convention: str
    angle_convention: str
    mass_weighting_convention: str
    equilibrium_bond_lengths_bohr: tuple[float, ...]
    version: str
    fingerprint: str
    sections: tuple[tuple[str, int, int], ...]

    def __post_init__(self) -> None:
        n_rows = len(self.rows)
        if not (
            n_rows == len(self.row_types)
            == len(self.units)
            == len(self.periodicities)
        ):
            raise ValueError("Coordinate-definition row metadata lengths do not agree.")
        if (
            self.bond_convention == "equilibrium_scaled_distance"
            and len(self.equilibrium_bond_lengths_bohr)
            != sum(kind == "bond" for kind in self.row_types)
        ):
            raise ValueError(
                "Equilibrium-scaled charts require one persisted length per bond row."
            )


@dataclass(frozen=True)
class SymmetryOperation:
    """Validated affine coordinate/derivative map for one rotor operation."""

    operation_id: str
    source_rotor_id: str
    symmetry_order: int
    atom_permutation: tuple[int, ...]
    row_permutation: tuple[int, ...]
    row_signs: tuple[int, ...]
    row_offsets: tuple[float, ...]
    phase_offset: float
    geometry_validated: bool
    derivative_validated: bool
    max_geometry_error: float
    max_energy_error: float
    max_gradient_error: float
    max_hessian_error: float
    validation_hash: str
    exact_atom_permutation: bool = False
    coverage_validated: bool = False
    # q'_i = row_signs_i * row_scales_i * q_{perm(i)} + row_offsets_i.  The scale
    # is unity for angles, torsions and plain-distance bonds.  An
    # equilibrium-scaled bond row carries its own reference length, so two rows
    # exchanged by the permutation are related by (r_eq_i / r_eq_perm(i))**2.
    row_scales: tuple[float, ...] = ()

    @property
    def validated(self) -> bool:
        if self.exact_atom_permutation:
            return self.derivative_validated
        return self.geometry_validated and self.derivative_validated


@dataclass(frozen=True)
class GroupedModelPolicy:
    """All numerical choices that affect construction or evaluation."""

    training_phases_degrees: tuple[float, ...] = (0.0, 30.0, 60.0, 90.0, 120.0)
    held_out_phases_degrees: tuple[float, ...] = (15.0, 45.0, 75.0, 105.0)
    group_rotation_training_phases_degrees: tuple[float, ...] = (
        0.0, 60.0, 120.0, 180.0, 240.0, 300.0
    )
    group_rotation_held_out_phases_degrees: tuple[float, ...] = (
        30.0, 90.0, 150.0, 210.0, 270.0, 330.0
    )
    enforce_exact_permutation_symmetry: bool = True
    construct_candidate_interaction_banks: bool = False
    interaction_symmetry_mode: str = "exact_pair_state_images_v1"
    shepard_p: float = 3.0
    shepard_q: float = 4.0
    confidence_radius: float = 1.0
    exact_center_tolerance: float = 1.0e-10
    phase_tolerance: float = 2.0e-5
    signature_scale: float = 1.0
    relaxation_policy_id: str = "rigid"
    projector_policy_id: str = "active_torsions_only_v1"
    anchor_policy_id: str = "single_canonical_anchor_v1"
    weight_metric_id: str = "periodic_signature_mean_v1"
    geometry_symmetry_tolerance_bohr: float = 0.35
    # An exact atom permutation of a symmetric rotor is a symmetry of the
    # potential, so its images are trusted over the whole period by default and
    # ``rigid_symmetry_defect_bohr`` is only reported.  Setting a tolerance here
    # makes a rotor whose own geometry is distorted by more than that fall back
    # to a physically sampled 360-degree phase grid, which costs one calculation
    # per extra phase and only pays off for rigidly driven torsions.
    rigid_symmetry_defect_tolerance_bohr: float | None = None
    # The rigidly rotated one-period structure is always calculated, as the
    # reference the symmetry validation is measured against.  Keeping it as a
    # bank state adds a free physical point, but its own permutation images are
    # not in the bank, so the factor stops being exactly closed under the rotor
    # group and the model loses exact permutation symmetry.  Exactness wins by
    # default; enable this when a rigidly driven torsion matters more.
    bank_includes_rigid_coverage_state: bool = False
    energy_symmetry_tolerance_hartree: float = 2.0e-3
    gradient_symmetry_tolerance_hartree_per_bohr: float = 2.0e-2
    hessian_symmetry_tolerance_hartree_per_bohr2: float = 1.0e-1
    coupling_energy_threshold_hartree: float = 2.0e-3
    coupling_gradient_rms_threshold: float = 2.0e-3
    coupling_probe_phases_degrees: tuple[tuple[float, float], ...] = (
        (30.0, 30.0),
        (30.0, 60.0),
        (60.0, 30.0),
        (60.0, 60.0),
        (90.0, 30.0),
        (30.0, 90.0),
    )
    sidechain_basin_id: str = "canonical_input_basin"
    extra: tuple[tuple[str, Any], ...] = field(default_factory=tuple)
