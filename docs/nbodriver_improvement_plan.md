# NboDriver implementation summary

This document summarizes the implemented `NboDriver` architecture in VeloxChem. The core roadmap described for the clean NBO/NRA/NRT driver has been completed, and remaining items are listed only as future extensions.

## Architecture

The driver is organized as a layered analysis pipeline:

```text
AO/SCF density and overlap
-> NAO/NPA transformation and diagnostics
-> MO composition in the NAO basis
-> NBO candidate generation
-> Lewis/resonance assignment
-> NRA/NRT density fitting
-> reports and structured diagnostics
```

NRA/NRT remains a post-processing density-fit layer on top of generated Lewis/resonance alternatives. It does not feed back into candidate generation, primary Lewis assignment, or score/ranking weights.

## Implemented capabilities

- Conservative NAO/NPA construction with explicit foundation invariants: `T^T S T`, `D_NAO = T^T S P S T`, electron count, density symmetry, and charge conservation.
- Deterministic handling of near-degenerate atom-block NAO rotations.
- MO-in-NAO composition analysis for restricted and unrestricted orbitals.
- NBO candidate generation for `CR`, `LP`, `BD(sigma)`, `BD(pi)`, `SOMO`, radical one-electron `LP`, radical one-electron `BD(pi)`, `BD*`, and `RY`.
- Lewis/NBO primary assignment with explicit candidate electron counts.
- Bond and pi-bond constraints through `NboConstraints` or dictionary input.
- Lewis/resonance alternatives with separate `sigma_nbo_list` and `pi_nbo_list` fields.
- Fixed/active partition fields for resonance analysis: `fixed_nbo_list`, `active_nbo_list`, `active_pi_nbo_list`, `active_lone_pair_nbo_list`, and `active_one_electron_nbo_list`.
- Lewis electron accounting, atom valence electron counts, formal-charge diagnostics, octet diagnostics, and transparent `score_terms`.
- Lone-pair donation in the resonance active space.
- Closed-shell NRA/NRT density fitting with `selected`, `pi`, `valence`, and `full` subspaces.
- `nra_report(level="summary")` and `nra_report(level="full")` report formatters.
- Molecule-independent pi-bond structure signatures for NRA/NRT reports and prior matching.
- Open-shell spin-resolved NAO densities and one-electron radical alternatives.
- Open-shell NRA/NRT fitting through a combined total-density and spin-density residual.
- `BD*` antibonding complements as same-subspace orthogonal partners to occupied `BD` candidates.
- One-center `RY` acceptor complements as candidate-only records.
- Donor-acceptor diagnostics from occupied donors to candidate-only `BD*` and `RY` acceptors using a density-coupling diagnostic.
- Optional prior-regularized NRA/NRT weights with explicit prior metadata.
- Regression coverage for foundation invariants, NBO counts, constraints, sigma/pi partitions, Lewis accounting, score terms, lone-pair donation, structure-pool invariants, open-shell SOMO/radical alternatives, antibonding/Rydberg complements, donor-acceptor diagnostics, NRA/NRT weights, NRA/NRT priors, open-shell spin-resolved NRA/NRT, labels, and reports.

## Weight terminology

The implementation keeps three weight concepts separate:

| Name | Meaning |
| --- | --- |
| Score/ranking weight | Softmax weight derived from the Lewis assignment score. |
| NRA/NRT density weight | Nonnegative density-fit weight from `results["nra"]`. |
| VB/wavefunction weight | Future state-mixing quantity; not implemented. |

User priors are not a fourth physical weight. They are optional guidance used only in the regularized NRA/NRT density fit and are reported separately as prior metadata.

## Current scientific scope

The driver provides a transparent NBO/NRA/NRT analysis layer suitable for method development and interpretation. Several quantities are deliberately reported as diagnostics rather than final physical observables:

- Donor-acceptor diagnostics use the NAO density matrix as a coupling operator. They are not second-order perturbation energies because no NBO Fock matrix or energy denominator is used.
- NRA/NRT weights are density-reconstruction weights. They are not VB amplitudes or state-mixing coefficients.
- Pi-bond signatures are report annotations. The internal data model remains based on sigma/pi and fixed/active partitions, not named structure classes.
- NBO candidates are density-block natural orbitals in NAO space; a separate NHO directionality layer is outside the current implementation.

## Future extensions

The core roadmap is implemented. Useful future work includes:

- Add a true NHO/hybrid-direction layer on top of the current NAO-space candidate vectors.
- Add an NBO Fock-matrix layer and second-order donor-acceptor perturbation energies.
- Expand optional user-facing structure annotations without adding molecule-specific logic to the core driver.
- Add broader validation notebooks and publication examples for additional radicals, zwitterionic systems, amides, nitro compounds, and transition-metal fragments.
- Investigate VB/BOND-style state coupling if wavefunction weights are desired.
