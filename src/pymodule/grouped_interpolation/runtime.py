"""Production grouped evaluator shared by validation and dynamics."""

from __future__ import annotations

import numpy as np

from .coordinates import assert_coordinate_fingerprint, build_coordinate_definition
from .factors import (
    LocalShepardFactor,
    PairInteractionFactor,
    evaluate_projected_taylor,
)
from .registry import GroupedModelRegistry, find_published_model
from .symmetry import transform_internal_taylor_state


class GroupedRuntimeModel:
    """Load one published signed-residual model and evaluate E/gradient."""

    def __init__(self, bundle):
        self.bundle = bundle
        self.coordinate_definition = bundle["coordinate_definition_json"]
        self.metadata = bundle["metadata_json"]
        self.policy = bundle["policy_json"]
        self.points = bundle["points"]
        self.content_hash = bundle["content_hash"]
        self.status = bundle["status"]
        self.last_factor_diagnostics = {}
        expected = self.coordinate_definition["fingerprint"]
        for request_id, point in self.points.items():
            assert_coordinate_fingerprint(
                point["coordinate_fingerprint"], expected, request_id
            )
        anchor_id = self.metadata["core_point_label"]
        self.anchor_state = self.points[anchor_id]
        row_types = tuple(self.coordinate_definition["row_types"])
        subsets = {
            subset["subset_id"]: subset for subset in bundle["subset_registry_json"]
        }
        rotors = {
            rotor["rotor_id"]: rotor for rotor in bundle["motion_registry_json"]
        }
        exact_operations = {}
        for operation in bundle["symmetry_registry_json"]:
            if not (
                operation.get("exact_atom_permutation", False)
                and operation.get("derivative_validated", False)
            ):
                continue
            exact_operations.setdefault(operation["source_rotor_id"], []).append(
                operation
            )
        self.factors = []
        for subset_id, request_ids in sorted(bundle["factor_registry_json"].items()):
            subset = subsets[subset_id]
            rotor_id = subset["canonical_subset_key"][0]
            states = tuple(self.points[request_id] for request_id in request_ids)
            if any(state["internal_hessian"] is None for state in states):
                raise ValueError(f"Factor {subset_id} contains a state without a Hessian.")
            self.factors.append(
                (
                    subset_id,
                    LocalShepardFactor(
                        subset=subset,
                        rotor=rotors[rotor_id],
                        states=states,
                        anchor_state=self.anchor_state,
                        policy=self.policy,
                        row_types=row_types,
                    ),
                )
            )
        self.interaction_factors = []
        for interaction in bundle.get("interaction_registry_json", []):
            subset_id = interaction["subset_id"]
            # An interaction bank stores reference-minus-additive.  That
            # difference is only a valid correction if it is carried on the
            # whole chart, so the projector is normalized here rather than
            # trusted from the bundle: models published before this was fixed
            # persist the coupled torsions alone and would otherwise discard
            # every other component of the residual gradient.
            subset = dict(subsets[subset_id])
            subset["projector_rows"] = tuple(range(len(row_types)))
            pair_rotors = tuple(
                rotors[rotor_id] for rotor_id in interaction["rotor_ids"]
            )
            joint_states = []
            for record in interaction["state_records"]:
                joint = self.points[record["joint_request_id"]]
                state_a = self.points[record["state_a_request_id"]]
                state_b = self.points[record["state_b_request_id"]]
                if any(
                    state["internal_hessian"] is None
                    for state in (joint, state_a, state_b, self.anchor_state)
                ):
                    raise ValueError(
                        f"Interaction {subset_id} contains a state without a Hessian."
                    )
                joint_states.append(joint)
            # Symmetry images are taken of the raw joint reference states, not
            # of the inclusion-exclusion residual.  The assembled additive model
            # is a Taylor/Shepard construction in labelled internal coordinates
            # and is therefore not equivariant under the exact atom permutation,
            # so a residual carried unchanged onto an image centre would be
            # wrong by the additive model's own non-invariance.  Subtracting the
            # additive model again at the image geometry keeps every bank centre
            # -- physical or virtual -- an exact reproduction of the reference.
            physical_joint_states = tuple(joint_states)
            for source_state in physical_joint_states:
                for rotor_id in interaction["rotor_ids"]:
                    for operation in exact_operations.get(rotor_id, ()):
                        joint_states.append(
                            transform_internal_taylor_state(
                                source_state, operation
                            )
                        )
            residual_states = tuple(
                self._interaction_residual_state(joint, subset["projector_rows"])
                for joint in joint_states
            )
            self.interaction_factors.append(
                (
                    subset_id,
                    PairInteractionFactor(
                        subset=subset,
                        rotors=pair_rotors,
                        states=residual_states,
                        policy=self.policy,
                        row_types=row_types,
                    ),
                )
            )

    def _interaction_residual_state(self, joint_state, projector_rows):
        """Subtract the assembled additive model at this joint state's geometry."""

        coordinates = np.asarray(joint_state["internal_coordinates"], dtype=float)
        additive_energy, additive_gradient, _ = self._evaluate_additive_internal(
            coordinates
        )
        additive_hessian = self._finite_difference_additive_hessian(
            coordinates, projector_rows
        )
        return {
            "request_id": joint_state["request_id"],
            "energy": float(joint_state["energy"]) - additive_energy,
            "internal_coordinates": coordinates,
            "internal_gradient": (
                np.asarray(joint_state["internal_gradient"], dtype=float)
                - additive_gradient
            ),
            "internal_hessian": (
                np.asarray(joint_state["internal_hessian"], dtype=float)
                - additive_hessian
            ),
            "is_virtual_image": bool(joint_state.get("is_virtual_image", False)),
        }

    def _evaluate_additive_internal(
        self, current_internal_coordinates, *, include_virtual=True
    ):
        current = np.asarray(current_internal_coordinates, dtype=float)
        all_rows = tuple(range(len(self.coordinate_definition["rows"])))
        core_energy, core_gradient = evaluate_projected_taylor(
            self.anchor_state,
            current,
            tuple(self.coordinate_definition["row_types"]),
            all_rows,
        )
        energy = float(core_energy)
        gradient = np.asarray(core_gradient, dtype=float)
        diagnostics = {}
        for subset_id, factor in self.factors:
            evaluated = factor.evaluate(current, include_virtual=include_virtual)
            energy += float(evaluated["residual_energy"])
            gradient = gradient + evaluated["residual_gradient"]
            diagnostics[subset_id] = evaluated
        return energy, gradient, diagnostics

    def _finite_difference_additive_hessian(
        self, coordinates, projector_rows=None, step=1.0e-5
    ):
        """Differentiate the assembled additive gradient at an interaction center.

        ``projector_rows`` restricts the differentiation to the columns the
        interaction factor actually projects onto.  Every other column stays
        zero: it is never read back by ``evaluate_projected_taylor``, and
        computing it would cost one additive evaluation per unused coordinate.
        """

        center = np.asarray(coordinates, dtype=float)
        hessian = np.zeros((center.size, center.size), dtype=float)
        columns = (
            range(center.size)
            if projector_rows is None
            else sorted({int(row) for row in projector_rows})
        )
        for column in columns:
            plus = center.copy()
            minus = center.copy()
            plus[column] += step
            minus[column] -= step
            gradient_plus = self._evaluate_additive_internal(plus)[1]
            gradient_minus = self._evaluate_additive_internal(minus)[1]
            hessian[:, column] = (gradient_plus - gradient_minus) / (2.0 * step)
        return 0.5 * (hessian + hessian.T)

    @classmethod
    def from_hdf5(
        cls,
        filename,
        *,
        model_id=None,
        z_matrix=None,
        interpolation_settings=None,
        require_published=True,
    ):
        selected_id = find_published_model(filename, model_id) if require_published else model_id
        if selected_id is None:
            raise ValueError("model_id is required when loading a staging grouped model.")
        registry = GroupedModelRegistry(filename, selected_id)
        bundle = registry.load_bundle(require_published=require_published)
        if z_matrix is not None:
            settings = interpolation_settings or {}
            definition = build_coordinate_definition(
                z_matrix,
                use_inverse_bond_length=bool(settings.get("use_inverse_bond_length", True)),
                use_eq_bond_length=bool(settings.get("use_eq_bond_length", False)),
                use_cos_angle=bool(settings.get("use_cos_angle", False)),
                use_mass_weight=bool(settings.get("use_mass_weight", False)),
                eq_bond_lengths=bundle["coordinate_definition_json"].get(
                    "equilibrium_bond_lengths_bohr", ()),
            )
            assert_coordinate_fingerprint(
                definition.fingerprint,
                bundle["coordinate_definition_json"]["fingerprint"],
                f"runtime model {selected_id}",
            )
        return cls(bundle)

    @property
    def model_id(self):
        return self.metadata["model_id"]

    @property
    def eq_bond_lengths(self):
        return np.asarray(self.metadata["eq_bond_lengths_bohr"], dtype=float)

    def evaluate_internal(self, current_internal_coordinates, *, include_virtual=True):
        current = np.asarray(current_internal_coordinates, dtype=float)
        energy, gradient, diagnostics = self._evaluate_additive_internal(
            current, include_virtual=include_virtual
        )
        self.last_factor_diagnostics = diagnostics
        for subset_id, factor in self.interaction_factors:
            evaluated = factor.evaluate(current)
            energy += float(evaluated["residual_energy"])
            gradient = gradient + evaluated["residual_gradient"]
            self.last_factor_diagnostics[subset_id] = evaluated
        return energy, gradient

    def evaluate_coordinate(self, interpolation_coordinate, molecule):
        energy, internal_gradient = self.evaluate_internal(
            interpolation_coordinate.internal_coordinates_values
        )
        cartesian = interpolation_coordinate.b_matrix.T @ internal_gradient
        if self.coordinate_definition["mass_weighting_convention"] == "cartesian_inverse_sqrt_mass":
            if interpolation_coordinate.inv_sqrt_masses is None:
                raise ValueError("Grouped runtime requires mass weights for this model.")
            cartesian = cartesian / interpolation_coordinate.inv_sqrt_masses
        return float(energy), np.asarray(cartesian, dtype=float).reshape(-1, 3)
