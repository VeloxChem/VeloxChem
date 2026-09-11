"""Policy boundary for rigid now and constrained relaxation in later milestones."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class RelaxationRecord:
    policy_id: str
    converged: bool
    constraint_errors: tuple[float, ...]
    basin_preserved: bool


class RigidRelaxationPolicy:
    """Explicit no-relaxation policy used by the first grouped milestone."""

    policy_id = "rigid"

    def apply(self, generated_molecule):
        return generated_molecule, RelaxationRecord(
            policy_id=self.policy_id,
            converged=True,
            constraint_errors=(),
            basin_preserved=True,
        )
