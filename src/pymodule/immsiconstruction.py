#
#                                   VELOXCHEM
#              ----------------------------------------------------
#                          An Electronic Structure Code
#
#  SPDX-License-Identifier: BSD-3-Clause
#

"""Ordered schema-4 MSI construction and promotion from computed QM states."""

from __future__ import annotations

from dataclasses import replace
import math

import numpy as np

from .imcorefamilies import (
    CoreFamilySpec,
    CorePointSpec,
    LocalityAuditRecord,
    PromotionRecord,
    assign_core_candidate,
    audit_local_factor_state,
    descriptor_squared_distance,
    evaluate_environment_descriptor,
    select_descriptor_medoids,
    stable_msi_id,
)
from .imlocalfactorlibrary import FactorBindingSpec
from .immsiregistry import (
    MSIModelRegistry,
    adapt_legacy_banks_to_schema4,
    validate_msi_registry,
    write_msi_registry_atomic,
)


SCHEMA4_CONSTRUCTION_STAGES = (
    "detect_topology",
    "generate_provisional_states",
    "compute_locality_environment_audits",
    "classify_states",
    "finalize_projectors_and_coupling",
    "create_or_reuse_factor_definitions",
    "validate_bindings",
    "promote_core_families_and_points",
    "write_schema4_registry",
)


class Schema4MSIConstructor:
    """Build a validated registry after all provisional QM jets exist.

    Legacy signed-factor banks are used as a compatibility representation for
    already-computed provisional states.  No shared binding is inferred: new
    promoted families receive explicit disabled binding records until a target
    probe or complete native library supplies evidence.
    """

    def __init__(
            self,
            locality_threshold_policy=None,
            promotion_policy=None,
            stage_observer=None):
        self.locality_threshold_policy = dict(
            locality_threshold_policy or {}
        )
        self.promotion_policy = {
            "policy_id": "schema4_promotion_v1",
            "family_domain_radius": 0.50,
            "max_new_family_core_points": 3,
            "max_existing_family_core_points": 2,
            "core_point_min_descriptor_separation": 0.15,
            "contact_regime_cutoff": 0.15,
        }
        self.promotion_policy.update(dict(promotion_policy or {}))
        self.stage_observer = stage_observer
        self.stage_trace = []

    def _stage(self, name):
        expected = SCHEMA4_CONSTRUCTION_STAGES[len(self.stage_trace)]
        if name != expected:
            raise RuntimeError(
                f"Schema-4 construction stage {name!r} occurred before "
                f"required stage {expected!r}."
            )
        self.stage_trace.append(name)
        if self.stage_observer is not None:
            self.stage_observer(name)

    def construct(
            self,
            legacy_banks,
            root,
            descriptor_spec,
            output_file=None):
        """Run all construction stages and optionally atomically persist them."""

        self.stage_trace = []
        self._stage("detect_topology")
        topology = self._detect_topology(legacy_banks, root)

        self._stage("generate_provisional_states")
        provisional_states = self._collect_provisional_states(topology)

        self._stage("compute_locality_environment_audits")
        audited_states = self._audit_states(
            provisional_states, descriptor_spec
        )

        self._stage("classify_states")
        # ``audit_local_factor_state`` assigns the class from persisted raw
        # metrics; keeping this as a separate stage makes ordering observable.
        audits = tuple(item["audit"] for item in audited_states)

        self._stage("finalize_projectors_and_coupling")
        registry = adapt_legacy_banks_to_schema4(legacy_banks, root)
        registry = self._install_descriptor(registry, descriptor_spec, topology)

        self._stage("create_or_reuse_factor_definitions")
        # The adapter creates one complete native definition per legacy term.
        # Reuse remains identity-based and is only introduced by validated
        # bindings, never by guessing from family labels.
        definitions = registry.factor_definitions
        if len(definitions) != len({
                item.factor_definition_id for item in definitions}):
            raise RuntimeError("Duplicate schema-4 factor definition identity.")

        self._stage("validate_bindings")
        validate_msi_registry(registry)

        self._stage("promote_core_families_and_points")
        registry, promotions = self._promote_candidates(
            registry, descriptor_spec, audited_states
        )
        validate_msi_registry(registry)

        self._stage("write_schema4_registry")
        if output_file is not None:
            write_msi_registry_atomic(output_file, registry)
        return registry, audits, promotions

    @staticmethod
    def _detect_topology(legacy_banks, root):
        topology = []
        for family_label, bank in sorted(
                (legacy_banks or {}).items(), key=lambda item: str(item[0])):
            core = bank.get("core")
            if core is None:
                raise ValueError(
                    f"Family {family_label!r} has no provisional core point."
                )
            for term_id, term in sorted(
                    bank.get("terms", {}).items(),
                    key=lambda item: str(item[0])):
                class_id = stable_msi_id(
                    "factor_class",
                    {
                        "root": str(root),
                        "legacy_family_label": str(family_label),
                        "legacy_term_id": str(term_id),
                    },
                )
                topology.append({
                    "family_label": str(family_label),
                    "bank": bank,
                    "core": core,
                    "term_id": str(term_id),
                    "term": term,
                    "factor_class_id": class_id,
                })
        return tuple(topology)

    @staticmethod
    def _collect_provisional_states(topology):
        provisional = []
        for entry in topology:
            term = entry["term"]
            anchors = {
                str(value) for value in term.get("anchor_state_ids", ())
            }
            states = sorted(
                term.get("expected_states", {}).items(),
                key=lambda item: str(item[0]),
            )
            if not states:
                raise ValueError(
                    f"Factor term {entry['term_id']!r} has no provisional states."
                )
            if not anchors:
                anchors = {str(states[0][0])}
            anchor_item = next(
                (
                    (state_id, datapoint)
                    for state_id, datapoint in states
                    if str(state_id) in anchors
                ),
                states[0],
            )
            for state_id, datapoint in states:
                if str(state_id) in anchors:
                    continue
                item = dict(entry)
                item.update({
                    "state_id": str(state_id),
                    "datapoint": datapoint,
                    "anchor_datapoint": anchor_item[1],
                    "datapoint_label": _datapoint_label(datapoint),
                })
                provisional.append(item)
        return tuple(provisional)

    def _audit_states(self, provisional_states, descriptor_spec):
        audited = []
        for item in provisional_states:
            state = item["datapoint"]
            anchor = item["anchor_datapoint"]
            core = item["core"]
            state_descriptor = _descriptor_value(descriptor_spec, state)
            core_descriptor = _descriptor_value(descriptor_spec, core)
            distance_squared = descriptor_squared_distance(
                descriptor_spec, state_descriptor, core_descriptor
            )
            environment_distance = math.sqrt(max(distance_squared, 0.0))
            contact_changed = _contact_regime_changed(
                descriptor_spec,
                core_descriptor,
                state_descriptor,
                float(self.promotion_policy["contact_regime_cutoff"]),
            )
            term = item["term"]
            projector_rows = tuple(
                int(row) for row in (
                    term.get("projector_rows", ())
                    or term.get("active_rows", ())
                )
            )
            boundary_rows = tuple(
                int(row) for row in term.get("response_rows", ())
                if int(row) not in projector_rows
            )
            periodic_rows = _periodic_rows(state)
            try:
                audit = audit_local_factor_state(
                    datapoint_label=item["datapoint_label"],
                    factor_class_id=item["factor_class_id"],
                    state_internal_coordinates=state.internal_coordinates_values,
                    anchor_internal_coordinates=(
                        anchor.internal_coordinates_values
                    ),
                    state_internal_gradient=state.internal_gradient,
                    anchor_internal_gradient=anchor.internal_gradient,
                    state_internal_hessian=state.internal_hessian,
                    projector_rows=projector_rows,
                    environment_distance=environment_distance,
                    contact_regime_changed=contact_changed,
                    boundary_rows=boundary_rows,
                    periodic_rows=periodic_rows,
                    mapping_compatible=_mapping_is_compatible(state, anchor),
                    threshold_policy=self.locality_threshold_policy,
                )
            except (TypeError, ValueError, AttributeError) as exc:
                audit = _ambiguous_audit(item, environment_distance, exc)
            audited.append({
                **item,
                "audit": audit,
                "descriptor": state_descriptor,
                "core_descriptor": core_descriptor,
            })
        return tuple(audited)

    def _install_descriptor(self, registry, descriptor_spec, topology):
        core_descriptor_by_family = {}
        for item in topology:
            core_descriptor_by_family.setdefault(
                item["family_label"],
                _descriptor_value(descriptor_spec, item["core"]),
            )
        family_label_by_id = {
            family.core_family_id: str(
                family.provenance.get("legacy_family_label", "")
            )
            for family in registry.core_families
        }
        families = tuple(
            replace(
                family,
                descriptor_spec_id=descriptor_spec.descriptor_spec_id,
                prototype_descriptor=tuple(
                    core_descriptor_by_family[family_label_by_id[
                        family.core_family_id]]
                ),
                descriptor_domain_radius=float(
                    self.promotion_policy["family_domain_radius"]
                ),
                construction_policy_id="schema4_ordered_v1",
            )
            for family in registry.core_families
        )
        family_by_id = {
            family.core_family_id: family for family in families
        }
        core_by_label = {
            _datapoint_label(item["core"]): item["core"] for item in topology
        }
        points = tuple(
            replace(
                point,
                reference_descriptor=tuple(_descriptor_value(
                    descriptor_spec, core_by_label[point.datapoint_label]
                )),
            )
            for point in registry.core_points
            if point.core_family_id in family_by_id
        )
        model_info = dict(registry.model_info)
        model_info["outer_distance_policy"] = "legacy_plus_contact_torsion_v1"
        model_info["created_by"] = "schema4_ordered_constructor_v1"
        model_info["construction_stage_order"] = list(
            SCHEMA4_CONSTRUCTION_STAGES
        )
        return replace(
            registry,
            model_info=model_info,
            descriptor_specs=(descriptor_spec,),
            core_families=families,
            core_points=points,
        )

    def _promote_candidates(self, registry, descriptor_spec, audited_states):
        family_by_legacy_label = {
            str(family.provenance.get("legacy_family_label", "")): family
            for family in registry.core_families
        }
        promotion_pairs = []
        for item in audited_states:
            source_family = family_by_legacy_label[item["family_label"]]
            identity = " ".join((
                str(item["term_id"]),
                str(item["term"].get("relaxation_policy_id", "")),
            )).lower()
            core_generating_motion = any(
                kind in identity for kind in (
                    "anchored_linker", "alkyl_chain", "amide_side_chain"
                )
            )
            if (
                item["audit"].classification == "core_candidate"
                and not item["audit"].contact_regime_changed
                and not core_generating_motion
            ):
                promotion = _no_action_factor_response_promotion(
                    item["audit"],
                    source_family.core_family_id,
                    str(self.promotion_policy["policy_id"]),
                )
            else:
                promotion = assign_core_candidate(
                    item["audit"],
                    item["descriptor"],
                    registry.core_families,
                    descriptor_spec,
                    source_core_family_id=source_family.core_family_id,
                    policy_id=str(self.promotion_policy["policy_id"]),
                )
            promotion_pairs.append((item, promotion))

        families = list(registry.core_families)
        points = list(registry.core_points)
        bindings = list(registry.factor_bindings)
        index = {
            str(label): dict(record)
            for label, record in registry.datapoint_index.items()
        }

        existing_candidates = [
            pair for pair in promotion_pairs
            if pair[1].decision == "add_core_point_to_existing_family"
            and pair[1].descriptor_distance >= float(
                self.promotion_policy[
                    "core_point_min_descriptor_separation"
                ]
            )
        ]
        selected = []
        existing_by_family = {}
        for pair in existing_candidates:
            existing_by_family.setdefault(
                pair[1].target_core_family_id, []
            ).append(pair)
        for family_id, candidates in sorted(existing_by_family.items()):
            descriptor_array = np.asarray([
                pair[0]["descriptor"] for pair in candidates
            ])
            indices = select_descriptor_medoids(
                descriptor_array,
                int(self.promotion_policy["max_existing_family_core_points"]),
                descriptor_spec.feature_scales,
            )
            selected.extend(candidates[index] for index in indices)
        new_family_candidates = [
            pair for pair in promotion_pairs
            if pair[1].decision == "create_core_family"
        ]
        if new_family_candidates:
            descriptors = np.asarray(
                [pair[0]["descriptor"] for pair in new_family_candidates]
            )
            selected_indices = select_descriptor_medoids(
                descriptors,
                int(self.promotion_policy["max_new_family_core_points"]),
                descriptor_spec.feature_scales,
            )
            representatives = [
                new_family_candidates[index] for index in selected_indices
            ]
            prototype = np.asarray(representatives[0][0]["descriptor"])
            family_id = stable_msi_id(
                "family",
                {
                    "root": registry.root,
                    "descriptor_spec_id": descriptor_spec.descriptor_spec_id,
                    "prototype": prototype.tolist(),
                    "policy_id": self.promotion_policy["policy_id"],
                },
            )
            disabled_binding_ids = []
            for factor_class in registry.factor_classes:
                binding_id = stable_msi_id(
                    "factor_binding",
                    {
                        "target_family": family_id,
                        "factor_class": factor_class.factor_class_id,
                        "mode": "disabled",
                    },
                )
                disabled_binding_ids.append(binding_id)
                bindings.append(FactorBindingSpec(
                    factor_binding_id=binding_id,
                    target_core_family_id=family_id,
                    factor_class_id=factor_class.factor_class_id,
                    factor_definition_id="",
                    mode="disabled",
                    applicability_domain={
                        "policy": "requires_target_probe_or_native_library"
                    },
                    provenance={
                        "reason": "sharing_not_inferred",
                        "next_action": "add_factor_validation_probe",
                    },
                ))
            distances = [
                math.sqrt(max(descriptor_squared_distance(
                    descriptor_spec, pair[0]["descriptor"], prototype
                ), 0.0))
                for pair in new_family_candidates
            ]
            radius = max(
                float(self.promotion_policy["family_domain_radius"]),
                max(distances, default=0.0),
            )
            families.append(CoreFamilySpec(
                core_family_id=family_id,
                root=registry.root,
                descriptor_spec_id=descriptor_spec.descriptor_spec_id,
                prototype_descriptor=tuple(prototype),
                descriptor_domain_radius=radius,
                factor_binding_ids=(),
                construction_policy_id=str(self.promotion_policy["policy_id"]),
                provenance={
                    "promotion_candidate_labels": [
                        pair[0]["datapoint_label"]
                        for pair in new_family_candidates
                    ],
                    "selected_medoid_labels": [
                        pair[0]["datapoint_label"] for pair in representatives
                    ],
                    "pending_disabled_binding_ids": disabled_binding_ids,
                },
            ))
            for item, promotion in representatives:
                promotion = replace(
                    promotion, target_core_family_id=family_id
                )
                selected.append((item, promotion))

        selected_by_label = {
            item["datapoint_label"]: (item, promotion)
            for item, promotion in selected
        }
        final_promotions = []
        for item, promotion in promotion_pairs:
            if item["datapoint_label"] in selected_by_label:
                promotion = selected_by_label[item["datapoint_label"]][1]
            final_promotions.append(promotion)
            record = index[item["datapoint_label"]]
            record["locality_audit"] = item["audit"].to_payload()
            record["promotion_record"] = promotion.to_payload()
            if item["datapoint_label"] not in selected_by_label:
                continue
            roles = set(record.get("model_roles", ()))
            roles.add(str(record.get("model_role", "residual_state")))
            roles.update(("residual_state", "outer_core"))
            record["model_roles"] = sorted(roles)
            record["full_pes_cardinal"] = True
            record["promoted_from_residual"] = True
            target_family_id = promotion.target_core_family_id
            point_id = stable_msi_id(
                "core_point",
                {
                    "family": target_family_id,
                    "label": item["datapoint_label"],
                },
            )
            record["core_family_id"] = target_family_id
            record["core_point_id"] = point_id
            points.append(CorePointSpec(
                core_point_id=point_id,
                core_family_id=target_family_id,
                datapoint_label=item["datapoint_label"],
                reference_descriptor=tuple(item["descriptor"]),
                confidence_radius=float(
                    getattr(item["datapoint"], "confidence_radius", 0.0) or 0.0
                ),
                full_pes_cardinal=True,
                second_order_cardinal=True,
                provenance={
                    "promoted_from_residual": True,
                    "locality_audit_id": item["audit"].locality_audit_id,
                },
            ))

        model_info = dict(registry.model_info)
        model_info["locality_threshold_policy"] = dict(
            self.locality_threshold_policy
        )
        model_info["promotion_policy"] = dict(self.promotion_policy)
        return replace(
            registry,
            model_info=model_info,
            core_families=tuple(families),
            core_points=tuple(points),
            factor_bindings=tuple(bindings),
            datapoint_index=index,
        ), tuple(final_promotions)


def _datapoint_label(datapoint):
    label = getattr(datapoint, "point_label", None)
    if label is None or not str(label).strip():
        raise ValueError("A provisional datapoint has no stable point_label.")
    return str(label)


def _descriptor_value(spec, datapoint):
    value, _ = evaluate_environment_descriptor(
        spec, np.asarray(datapoint.cartesian_coordinates, dtype=np.float64)
    )
    return np.asarray(value, dtype=np.float64)


def _contact_regime_changed(spec, anchor, state, cutoff):
    contact_indices = []
    for index, frozen_feature in enumerate(spec.feature_definitions):
        kind = str(frozen_feature.get("kind", ""))
        if kind in ("smooth_contact", "smooth_contact_pool"):
            contact_indices.append(index)
    if not contact_indices:
        return False
    scales = np.asarray(spec.feature_scales, dtype=np.float64)
    anchor_contact = float(np.sum(
        np.asarray(anchor)[contact_indices] * scales[contact_indices]
    ))
    state_contact = float(np.sum(
        np.asarray(state)[contact_indices] * scales[contact_indices]
    ))
    return (anchor_contact >= cutoff) != (state_contact >= cutoff)


def _periodic_rows(datapoint):
    z_matrix = getattr(datapoint, "z_matrix_dict", None) or {}
    bond_count = len(z_matrix.get("bonds", ()))
    angle_count = len(z_matrix.get("angles", ()))
    dihedral_count = len(z_matrix.get("dihedrals", ()))
    improper_count = len(z_matrix.get("impropers", ()))
    start = bond_count + angle_count
    return tuple(range(start, start + dihedral_count + improper_count))


def _mapping_is_compatible(state, anchor):
    state_masks = getattr(state, "mapping_masks", None)
    anchor_masks = getattr(anchor, "mapping_masks", None)
    if state_masks is None or anchor_masks is None:
        return True
    state_masks = np.asarray(state_masks)
    anchor_masks = np.asarray(anchor_masks)
    return (
        state_masks.ndim == 2
        and anchor_masks.ndim == 2
        and state_masks.shape[1:] == anchor_masks.shape[1:]
    )


def _ambiguous_audit(item, environment_distance, exc):
    policy_id = str(
        item.get("threshold_policy_id", "schema4_locality_unavailable_v1")
    )
    audit_id = stable_msi_id(
        "locality_audit",
        {
            "datapoint_label": item["datapoint_label"],
            "factor_class_id": item["factor_class_id"],
            "reason": type(exc).__name__,
        },
    )
    return LocalityAuditRecord(
        locality_audit_id=audit_id,
        datapoint_label=item["datapoint_label"],
        factor_class_id=item["factor_class_id"],
        omitted_displacement=0.0,
        omitted_gradient_fraction=0.0,
        cross_hessian_fraction=0.0,
        environment_distance=float(environment_distance),
        contact_regime_changed=False,
        classification="ambiguous",
        threshold_policy_id=policy_id,
        metrics={
            "mapping_compatible": False,
            "audit_error": f"{type(exc).__name__}: {exc}",
        },
    )


def _no_action_factor_response_promotion(audit, source_family_id, policy_id):
    """Route non-core local-factor failures to probe/native-library handling."""

    promotion_id = stable_msi_id(
        "promotion",
        {
            "datapoint_label": audit.datapoint_label,
            "locality_audit_id": audit.locality_audit_id,
            "decision": "no_action",
            "reason": "local_factor_response_requires_validation",
            "policy_id": policy_id,
        },
    )
    return PromotionRecord(
        promotion_record_id=promotion_id,
        datapoint_label=audit.datapoint_label,
        locality_audit_id=audit.locality_audit_id,
        decision="no_action",
        source_core_family_id=source_family_id,
        reason_codes=(
            "local_factor_response_requires_validation",
            "do_not_promote_reusable_factor_state",
        ),
        descriptor_distance=audit.environment_distance,
        policy_id=policy_id,
    )
