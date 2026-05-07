# Valence Bond Driver Implementation Plan

## Purpose

This document tracks the current implementation strategy for the VeloxChem valence bond (VB) driver. The driver is being built incrementally from small, inspectable fixed-orbital models toward general non-orthogonal VB, BOVB, VBB, and eventually compact user-facing VB specifications.

The guiding principle is:

> Validate chemically transparent fixed-orbital active spaces before adding orbital optimization or breathing-orbital models.

BOVB is still an important target, but it should be introduced only after the active/inactive orbital partition, fixed-orbital Hamiltonian, spin algebra, and diagnostics are reliable for more than a two-orbital test case.

## Current architecture

### Shared orbital analysis

`OrbitalAnalyzer` is the shared owner of orbital recognition. It provides:

- AO-to-atom and AO-to-angular-momentum maps.
- NAO/NPA data from an RHF/UHF reference.
- Shared candidate records for `CR`, `BD(sigma)`, `BD(pi)`, `LP`, `SOMO`, `BD*`, and `RY` candidates.
- Metal-ligand `ML` diagnostic records with separate `sigma-acceptor` and `pi-donor` channels.
- Candidate labels and atom assignments consumed consistently by NBO and VB.

Candidate classification is implemented inside `orbitalanalyzerdriver.py` as analyzer-private helper logic. Public workflows should not depend on a separate classifier module.

### VB driver responsibilities

`VbDriver` consumes the analyzer payload and owns all VB-specific logic:

- active-space selection,
- active/inactive/frozen candidate partitioning,
- VB active orbital construction,
- VB structure generation,
- overlap and Hamiltonian construction,
- generalized eigenvalue solution,
- structure weights and diagnostics.

The analyzer identifies orbitals; the VB driver decides how those orbitals become a wavefunction.

For metal-ligand systems, the analyzer records are intentionally diagnostic: `ML/sigma-acceptor` marks ligand-to-metal sigma donation, while `ML/pi-donor` marks metal-to-ligand pi back-donation into ligand pi-type nonbonding/acceptor space. `VbDriver` can now explicitly select these records as fixed-orbital active-space seeds. The analyzer does not activate them automatically, and NBO still keeps them out of the primary Lewis table.

## Current status: 2026-05-05

`VbDriver` is currently a fixed-orbital VB development driver with shared analyzer input, validated two-electron spin-adapted algebra, frozen-HF embedding, and a determinant-CI fallback for larger π active spaces.

Handoff status for the current `nbo` branch work:

- H2 is currently the only fully accepted VBCI/VBSCF/BOVB method ladder in this stabilization phase.
- Allyl cation/radical/anion VBCI, compact-CSF, and common-orbital VBSCF routes have been implemented or updated, but the latest split-valence rows still need a clean rerun before replacing the `rerun pending` table entries below.
- Benzene fixed-orbital VBCI and compact-CSF gates are present in the notebook. Benzene independent-amplitude VBSCF is explicitly box-limited and must not be promoted to a final reference.
- Benzene VBSCF now defaults to the symmetry-preserving equivalent-center common-amplitude path for the 6e/6pi singlet checkpoint, using a conservative default breathing-amplitude bound. Treat this as the stabilized VBSCF checkpoint to validate before any determinant-space BOVB work, not an accepted final benzene VBSCF reference yet.
- Benzene determinant-space BOVB is not implemented. The current organic π BOVB route for larger determinant spaces is still the conservative fixed-orbital zero-amplitude wrapper and should not be described as final benzene BOVB.

Current organic bond-breaking progress:

- H2 is the strict two-electron validation gate. VBCI, VBSCF, and H2 BOVB now approach the UHF 5.0 Angstrom asymptote, while RHF correctly fails at dissociation.
- Ethane C-C is split into two roles:
   - the rigid scan is an active-space regression test for state tracking, frozen-HF embedding, and stretched-orbital continuity;
   - the UHF-relaxed constrained scan is the chemically meaningful validation path, using fixed C-C distances from 1.30 to 5.00 Angstrom and evaluating RHF/UHF/VB single points on the same relaxed geometries.
- Ethane stretched sigma bonds now use full-system UHF frontier active orbitals from 1.8 Angstrom onward. This removes the previous active-orbital branch spike around 1.9-2.5 Angstrom.
- Ethane VBCI is currently the useful fixed-orbital VB baseline. VBSCF remains diagnostic because it is still numerically identical to VBCI for this model.
- Non-H2 two-orbital BOVB is deliberately held at the fixed-orbital limit for now. This removes the discontinuous optimized/fixed BOVB branch in ethane, but it is not a completed polyatomic BOVB model.
- Ethylene C-C remains diagnostic only: the current pi-only two-electron active space cannot describe full C=C cleavage. A four-electron sigma+pi active space is required before interpreting ethylene as a physical dissociation curve.

Implemented and validated:

- Shared `OrbitalAnalyzer` consumption for active-orbital selection and active/inactive/frozen candidate diagnostics.
- H₂ one-active-bond VB-CI with covalent/ionic structures, spin-adapted two-electron algebra, generalized overlap metric pruning, and Chirgwin-Coulson/Löwdin weights.
- H₂ VB-SCF now uses a common, symmetry-preserving center-local breathing-orbital relaxation in split-valence bases, while fixed-orbital VB-CI remains the baseline.
- H₂ two-orbital BOVB with structure-specific, center-local breathing active orbitals and regression coverage against the fixed-orbital VB-CI baseline.
- Stable automatic H₂ active-space selection: H₂ scans now use an atom-centered H-H active pair instead of independently selected analyzer `BD_sigma_*` candidates, avoiding active-label switches along dissociation curves.
- Allyl-cation compact-CSF BOVB for the two-electron three-center π singlet checkpoint, using center-local π breathing directions in split-valence bases.
- Ethylene one-active-π VB-CI with frozen sigma/core HF embedding.
- Fixed-orbital multi-center π active spaces through allyl cation, allyl radical, allyl anion, butadiene, and benzene.
- Compact spin-adapted CSF Hamiltonian mode for the allyl cation 2e/3π singlet benchmark, using adjacent π-bond graph templates and the full fixed-orbital VB-CI result as the regression oracle.
- Orthonormal determinant-CI fallback for larger and open-shell active spaces.
- Chemically readable determinant labels and grouped spatial-occupation resonance diagnostics.
- Graph-generated spin-adapted chemical resonance/CSF template projection diagnostics for allyl radical, allyl anion, butadiene, and benzene.
- Graph-automorphism averaged displayed resonance weights with raw unsymmetrized diagnostics retained.
- Singlet determinant-CI root selection by alpha/beta exchange symmetry before CSF projection.
- CSF template phase convention matched to the determinant-CI alpha/beta occupation-string representation.
- Explicit metal-ligand fixed-orbital active spaces through `active_metal_ligand_channels`, including sigma-only and combined sigma-plus-pi back-donation determinant-CI models.
- Metal-ligand diagnostics that label the active model as `sigma-only`, `sigma-plus-backdonation`, or a custom channel selection, and expose whether metal-to-ligand back-donation is blocked or enabled.
- Metal-ligand VB-SCF channel scans for both sigma-only and sigma-plus-backdonation determinant spaces, using a common-active-orbital reference with the same donation/back-donation switch metadata as BOVB.
- Metal-ligand BOVB for analyzer-selected determinant-CI channel spaces, using center-local breathing-orbital relaxation and retaining the fixed determinant-CI zero-amplitude limit as the internal baseline.
- ECP-aware VB AO Hamiltonians for Pd/def2-SVP and other ECP basis workflows, using the same effective-charge nuclear attraction plus ECP one-electron terms as the SCF driver.
- Real Pd--NH3 and Pd--PH3 notebook distance scans with B3LYP constrained geometries and HF total-energy references. These are the currently validated metal-ligand dissociation curves. The notebook also computes state-tracked sigma-channel VB-SCF/BOVB metal-ligand diagnostics, but those traces are not yet accepted as production dissociation-energy curves. Pi back-donation is currently tracked as an `ML/pi-donor` diagnostic strength, not yet as a finished four-electron sigma+pi VB total-energy dissociation curve.
- Notebook validation and source-level regression tests for the fixed-orbital π resonance diagnostics, including butadiene captured-subspace weight and singlet exchange parity.

Current boundaries:

- The multi-center π engine is still fixed-orbital; it is not BOVB.
- Larger π systems currently solve a determinant-CI reference and then project chemically compact CSF templates onto that root; only the allyl cation 2e/3π singlet path has been promoted to a compact CSF Hamiltonian so far.
- BOVB is currently reliable only for the validated H2 two-orbital path and the allyl-cation compact-CSF checkpoint. For non-H2 two-orbital organic active spaces, BOVB is forced to the fixed-orbital limit until a constrained, continuous, size-consistent breathing-orbital model is implemented. Metal-ligand BOVB channel scans remain diagnostics. Metal-ligand sigma-channel scans now support state-tracked active orbitals and corrected ECP total-energy reporting, but the resulting VB-SCF/BOVB sigma traces are not validated production potential curves. The combined sigma+pi back-donation determinant space still needs a size-consistent four-electron total-energy formulation before it should be interpreted quantitatively. Larger organic π systems beyond the allyl checkpoint are still fixed-orbital.

Immediate next method work:

- Implement a real VBSCF orbital-optimization path for embedded organic active bonds; the current ethane VBSCF trace is effectively VBCI.
- Implement constrained non-H2 BOVB breathing orbitals that vary continuously along dissociation and do not over-stabilize separated radical limits.
- Add a four-electron sigma+pi active space for ethylene C=C cleavage.

## VB branch stabilization phase: H2, allyls, and benzene

The immediate branch goal is no longer broad feature expansion. The priority is to stabilize a compact, reproducible validation ladder for fixed-orbital VBCI, VBSCF, and BOVB before resuming general sigma/pi active-space generalization.

The working workflow for this phase is:

1. Update `VbDriver` with one narrow method or diagnostic change.
2. Update `docs/test_vbdriver.ipynb` so the change is exercised on the stabilization ladder.
3. Update this implementation plan to record the change, current interpretation, and remaining open issues.
4. Compile and run the notebook/tests outside this document.
5. Use the reported numerical results or failures to drive the next narrow driver/notebook/documentation update.

### Stabilization validation ladder

The branch is considered stable only when the following systems give consistent and reproducible results for all three method levels:

| System | Active space | Spin | Required methods | Purpose |
| --- | --- | --- | --- | --- |
| H2 | 2e/2orb | singlet | VBCI, VBSCF, BOVB | strict two-electron dissociation and ionic/covalent gate |
| allyl cation | 2e/3pi | singlet | VBCI, VBSCF, BOVB | compact closed-shell three-center resonance gate |
| allyl radical | 3e/3pi | doublet | VBCI, VBSCF, BOVB | open-shell/root-selection and radical-delocalization gate |
| allyl anion | 4e/3pi | singlet | VBCI, VBSCF, BOVB | charge-delocalization and excess-electron resonance gate |
| benzene | 6e/6pi | singlet | VBCI, VBSCF, BOVB | aromatic Kekule/Dewar symmetry and larger pi-space gate |

Butadiene remains useful as an intermediate diagnostic between allyl and benzene, but it is not part of the minimal stabilization gate unless benzene failures require a smaller debug target.

### Required energy and resonance reporting

Every stabilization run must report the same energy quantities so VBCI, VBSCF, and BOVB can be compared without changing definitions between systems:

- `E_SCF`: Hartree-Fock reference total energy for the same molecule, geometry, charge, multiplicity, and basis. Use RHF for closed-shell singlets and ROHF/UHF-style unrestricted references for radicals or stretched-bond checks when appropriate.
- `E_VB`: total energy from the selected VB method.
- `E_corr(VB) = E_VB - E_SCF`: VB energy lowering relative to the HF reference.
- `E_best_template`: lowest localized single-template or single-structure energy inside the same active-space and embedding model.
- `E_res = E_best_template - E_VB`: resonance stabilization energy within the same VB model.
- `chemical_resonance_subspace_weight`: fraction of the determinant-CI root captured by the compact chemically labeled CSF/template subspace.
- Normalized `chemical_resonance_weights`, with graph-automorphism averaged displayed weights and raw unsymmetrized weights retained for diagnostics.

Interpretation rules:

- `E_corr(VB)` is not the same as resonance energy. It mixes active-space/static correlation, resonance, and any orbital-relaxation effect included in the chosen VB model.
- For frozen-pi VBCI, `E_corr(VB)` should be described as active-space/static-correlation lowering relative to HF, not as full dynamic correlation.
- `E_res` is the cleaner resonance diagnostic because it compares the full VB wavefunction against the best localized template in the same orbital and embedding model.
- Dynamic correlation outside the selected active space remains absent unless a later method explicitly adds it. This must be stated whenever VB energies are compared to SCF.

### Definition of compact CSF in this stabilization phase

`compact-CSF` is a chemically compressed projection of the same active-space Hamiltonian, not a replacement for the full VBCI reference. The driver first has a full VB/determinant basis with Hamiltonian `H` and overlap `S`. It then builds a small template matrix `T` whose columns are chemically labeled CSFs or graph templates expanded in that full basis, and solves the projected generalized eigenvalue problem:

- `H_compact = T^T H T`
- `S_compact = T^T S T`

For the allyl cation, this compact basis is the two symmetry-related spin-adapted allyl resonance CSFs corresponding to the two terminal C=C pi bond placements. For allyl radical and allyl anion, it is the graph-generated set of chemically recognizable radical or charge-delocalized templates.

The compact-CSF energy is therefore a subspace variational energy. Its strict ordering check is against `compact_csf_full_reference_energy`, the same-orbital full reference used internally to build the compact Hamiltonian. It should satisfy `E(compact-CSF) >= E(full same-orbital reference)` within numerical tolerance. A separately printed VBCI row is still useful for comparison, but it can be generated by a separate driver call with a different analyzer/SCF gauge; it is not the strict variational reference unless it is guaranteed to be the identical Hamiltonian/orbital object.

The purpose of compact-CSF is chemical interpretation: captured-subspace weight, resonance weights, and a compact resonance-template energy. `compact-CSF-BOVB` means the same compact CSF/template model after structure- or template-specific breathing relaxation; at present this relaxed compact BOVB path is implemented for allyl cation, while radical/anion determinant-space BOVB remains fixed-limit.

### Stabilization result table to fill from the notebook

| System | Method | E_SCF | E_VB | E_corr(VB) | E_best_template | E_res | captured weight | Status |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| H2 | VBCI | -1.1267553135 | -1.1053297868 | 2.143e-02 | -1.1040444070 | 1.285e-03 | n/a | finite; above RHF at equilibrium |
| H2 | VBSCF | -1.1267553135 | -1.1454799073 | -1.872e-02 | -1.1237346456 | 2.175e-02 | n/a | passing; common breathing lowers below VBCI/RHF |
| H2 | BOVB | -1.1267553135 | -1.1461071080 | -1.935e-02 | -1.1363736780 | 9.733e-03 | n/a | passing; below RHF and below VBSCF |
| allyl cation | VBCI | -116.0024522142 | -115.5723501546 | 4.301e-01 | -115.3966926015 | 1.757e-01 | n/a | finite fixed-orbital split-valence baseline |
| allyl cation | VBSCF | rerun pending | rerun pending | rerun pending | rerun pending | rerun pending | rerun pending | common center-local pi breathing implemented; rerun required |
| allyl cation | compact-CSF | rerun pending | rerun pending | rerun pending | rerun pending | rerun pending | rerun pending | compact/full-reference orbital mismatch fixed; rerun required |
| allyl cation | compact-CSF-BOVB | rerun pending | rerun pending | rerun pending | rerun pending | rerun pending | rerun pending | compact/full-reference orbital mismatch fixed; rerun required |
| allyl radical | VBCI | -114.9108292016 | -114.2734241819 | 6.374e-01 | -113.7917019908 | 4.817e-01 | 8.640e-02 | finite fixed-orbital determinant-CI baseline |
| allyl radical | VBSCF | rerun pending | rerun pending | rerun pending | rerun pending | rerun pending | rerun pending | common center-local pi breathing implemented; rerun required |
| allyl radical | BOVB-fixed | -114.9108292016 | -114.2734241819 | 6.374e-01 | -113.8219481851 | 4.515e-01 | 1.413e-01 | fixed-limit wrapper; zero BOVB lowering |
| allyl radical | compact-CSF | -114.9108292016 | -114.1647793047 | 7.460e-01 | -113.9857694228 | 1.790e-01 | 7.402e-02 | finite compact determinant-template diagnostic; 1.086e-01 a.u. above VBCI |
| allyl anion | VBCI | -114.7309008412 | -112.6244915370 | 2.106e+00 | -112.3297891023 | 2.947e-01 | 2.069e-01 | finite fixed-orbital determinant-CI baseline |
| allyl anion | VBSCF | rerun pending | rerun pending | rerun pending | rerun pending | rerun pending | rerun pending | common center-local pi breathing implemented; rerun required |
| allyl anion | BOVB-fixed | -114.7309008412 | -112.6244915370 | 2.106e+00 | -112.3368091531 | 2.877e-01 | 1.999e-01 | fixed-limit wrapper; zero BOVB lowering |
| allyl anion | compact-CSF | -114.7309008412 | -111.9887113601 | 2.742e+00 | -111.8441136782 | 1.446e-01 | 2.069e-01 | finite compact determinant-template diagnostic; 6.358e-01 a.u. above VBCI |
| benzene | VBCI | rerun pending | rerun pending | rerun pending | rerun pending | rerun pending | rerun pending | fixed-orbital 6e/6pi reference gate added to notebook |
| benzene | compact-CSF | rerun pending | rerun pending | rerun pending | rerun pending | rerun pending | rerun pending | Kekule/Dewar compact projection gate added to notebook |
| benzene | VBSCF | -230.6235016817 | -229.0714117103 | 1.552e+00 | -228.7087946074 | 3.626e-01 | 1.472e-28 | constrained common-pi-breathing checkpoint; all amplitudes at +/-0.35 bounds |
| benzene | BOVB | TBD | TBD | TBD | TBD | TBD | TBD | open |

### Acceptance criteria for the stabilization ladder

- All reported energies are finite and reproducible after restarting the notebook kernel.
- Structure weights, Löwdin weights, and normalized chemical resonance weights sum to one within numerical tolerance.
- Captured-subspace weight is reported separately from normalized resonance weights.
- Singlet determinant-CI roots use the lowest alpha/beta exchange-symmetric root before CSF-template projection.
- Symmetry-equivalent resonance templates have equal displayed weights after graph-automorphism averaging.
- H2 dissociation keeps a continuous atom-centered H-H active pair and approaches the broken-symmetry UHF asymptote while RHF correctly fails at stretched distance.
- Allyl cation has symmetric terminal resonance weights at symmetric geometry.
- Allyl radical keeps a stable doublet root and chemically sensible terminal radical delocalization.
- Allyl anion keeps symmetric terminal charge-delocalized resonance forms.
- Benzene reports symmetry-consistent Kekule and Dewar template families; Kekule-equivalent structures must be degenerate/equally weighted in the displayed diagnostics.
- VBSCF must not lie above the corresponding fixed-orbital VBCI result unless an explicit constraint explains the increase.
- BOVB must not lie above the corresponding VBCI/VBSCF result when the relevant breathing relaxations are enabled and no stabilizing constraint is disabled.

### Driver diagnostics required for this phase

`VbDriver` should expose the following diagnostics for every method path as soon as the needed internal quantities are available:

- `vb_method`: explicit method label such as `vbci`, `vbscf`, `bovb`, `compact-csf`, or `compact-csf-bovb`.
- `scf_reference_energy`: the HF total energy used for comparison when computed by the driver or notebook helper.
- `vb_correlation_energy`: `E_VB - E_SCF` when `E_SCF` is available.
- `localized_template_energies`: single-template/localized-structure energies in the same active-space and embedding model.
- `best_localized_template_energy`: the minimum of `localized_template_energies`.
- `resonance_energy`: `best_localized_template_energy - E_VB`.
- `chemical_resonance_subspace_weight`, `chemical_resonance_weights`, and `chemical_resonance_unsymmetrized_weights`.
- Method-specific orbital-relaxation metadata for VBSCF and BOVB, including active-orbital rotation/breathing amplitudes, convergence status, and fixed-orbital zero-amplitude baseline energy where applicable.

Current stabilization implementation note:

- First driver update in this phase adds a shared localized-template energy diagnostic built from the existing reported Hamiltonian and overlap matrices. It does not change the wavefunction model or reinvents the Hamiltonian; it reports single-basis-vector Rayleigh quotients `H_ii/S_ii` for the already constructed VB basis.
- The diagnostic keys are `localized_template_energy_model`, `localized_template_labels`, `localized_template_energies`, `best_localized_template_index`, `best_localized_template_label`, `best_localized_template_energy`, and `resonance_energy`.
- The same helper is now attached to two-orbital VBCI, H2/VBSCF, two-orbital BOVB, compact CSF, compact-CSF BOVB, and determinant-CI pi active-space results.
- Organic pi determinant spaces, including allyl radical and allyl anion, now have generalized method routes for `mode="vbscf"`, `mode="bovb"`, and `mode="compact-csf"`.
- Organic pi VBSCF now has a real common-orbital relaxation route for allyl cation, allyl radical, and allyl anion. It optimizes one center-local pi breathing amplitude per active orbital in a common active-orbital set shared by all structures/determinants. This is not BOVB because structures do not get separate orbital sets, but it is no longer a fixed-limit wrapper. The diagnostics include `organic_pi_vbscf_model="common-center-local-pi-breathing-orbitals"`, the initial VBCI-limit energy, optimized energy, relaxation amplitudes, and whether external breathing space was available.
- The radical/anion BOVB route is currently an explicit zero-amplitude fixed-orbital limit around the same determinant-CI solver. Fixed here also means no structure-specific breathing or other orbital relaxation: it reports `organic_pi_bovb_model="determinant-ci-fixed-orbital-zero-amplitude-limit"`, `bovb_used_fixed_orbital_limit=True`, and zero energy lowering. This keeps the method ladder executable while preventing over-interpretation as real BOVB relaxation.
- The radical/anion compact-CSF route now reuses the existing graph-generated chemical resonance templates as a working compact Hamiltonian subspace over the determinant-CI Hamiltonian. It reports the full determinant-CI reference energy, compact-basis error, captured-subspace weight, and localized-template resonance diagnostics.
- The allyl-cation compact-CSF row exposed a real consistency problem before the benzene gate: compact-CSF built its internal full reference with the frozen-HF effective Hamiltonian but did not apply the same frozen-space projection to the active orbitals that the VBCI route applies. That made the compact internal full reference lower than the displayed VBCI row, so compact-CSF could appear below VBCI in the method table.
- The compact-CSF and compact-CSF-BOVB allyl-cation paths now project active orbitals out of the frozen inactive density after building the frozen-HF embedding, matching the VBCI route internally. The allyl notebook no longer requires a separately computed VBCI row to be bitwise identical to the compact internal reference because separate `compute()` calls can regenerate SCF/analyzer gauges. Instead it reports the displayed-VBCI offset and enforces the variational ordering against the compact route's same-orbital internal full reference. The previous split-valence compact-CSF and compact-CSF-BOVB numbers are superseded and must be regenerated after rebuild.
- The allyl method-ladder assertions now apply the same-reference rule consistently: real VBSCF is checked against its same-call fixed-orbital initial energy, compact-CSF is checked against `compact_csf_full_reference_energy`, and radical/anion BOVB-fixed is checked against its same-call zero-amplitude initial energy rather than against a separately generated displayed VBCI row.
- A first benzene gate has been added to the notebook. It runs fixed-orbital benzene VBCI and compact-CSF in the 6e/6pi active space, verifies 400 determinant structures for the full reference, verifies the 15 compact templates split into two Kekule and thirteen Dewar pairings, reports captured-subspace weight, and enforces compact-CSF variational ordering against the internal same-orbital full reference.
- After the fixed-orbital benzene compact-CSF gate passed, a benzene VBSCF gate was added in 6-31G. It reuses the organic pi common-orbital center-local breathing optimizer, keeps one common active-orbital set for all 400 determinants, and checks the optimized VBSCF energy against the same-call fixed-orbital initial energy. This is a real VBSCF relaxation gate, not BOVB, because no structure/template gets its own orbital set.
- Benzene VBSCF reruns confirm that the independent per-center common-pi-breathing optimizer is box-limited at the default +/-0.35 amplitude bound. The latest 6-31G gate gives fixed-orbital VBCI -227.9017419716 a.u. and common-pi-breathing VBSCF -229.0714117103 a.u.; all six breathing amplitudes sit at +/-0.35 and the compact Kekule/Dewar captured weight collapses to 1.472e-28. This is a constrained-relaxation checkpoint, not yet a final benzene VBSCF reference.
- The independent-amplitude benzene bound scan remains active at every tested box. In 6-31G, bounds 0.05, 0.10, 0.20, and 0.35 give E_VB values -227.8892858027, -228.2042369878, -228.4454626817, and -229.2606545225 a.u.; every row reaches max|amplitude| equal to its bound. This is monotonic box-driven lowering and should not be promoted to a final VBSCF reference.
- `VbComputeOptions.orbital_amplitude_bound` now exposes the organic pi VBSCF breathing-amplitude box. `VbComputeOptions.orbital_relaxation_symmetry="equivalent-centers"` provides a one-parameter, symmetry-preserving common-amplitude path for benzene-like equivalent pi centers. For the benzene 6e/6pi singlet checkpoint, unset symmetry now auto-selects this equivalent-center route and defaults to a conservative +/-0.05 amplitude bound; explicit options still override both choices. The notebook includes a benzene symmetry-preserving scan cell; run it before any determinant-space BOVB implementation.
- Source tests were updated to require finite localized-template and resonance-energy diagnostics for H2 VBCI, allyl-cation compact CSF/BOVB, and the fixed-orbital pi resonance cases.
- Source tests were also added for allyl radical/anion VBCI, relaxed VBSCF, fixed-limit BOVB, compact-CSF method routes, and the benzene compact-CSF projection gate.
- The H2 VBSCF route was changed from a pure two-orbital rotation to a common, symmetry-preserving center-local breathing-orbital optimization. This removes the gauge-only behavior in split-valence bases while keeping one common active-orbital pair for all H2 structures.
- The updated H2 method-ladder notebook run gives finite VBCI, VBSCF, and BOVB rows. VBSCF now lowers from VBCI (-1.1053297868) to -1.1454799073 a.u. and below the RHF reference; BOVB remains slightly lower at -1.1461071080 a.u. because it allows structure-specific breathing. This satisfies the intended ordering `E(BOVB) <= E(VBSCF) <= E(VBCI)` for H2.
- H2 VBSCF still reports localized-template and ionic-pair diagnostics in the optimized orbital basis. The stable method diagnostics are total energy, energy lowering from the fixed-orbital baseline, resonance energy in the same optimized basis, and normalized structure-weight sums.

### Notebook requirements for `docs/test_vbdriver.ipynb`

The validation notebook should be reorganized around the stabilization table rather than isolated demonstrations:

1. Define helper functions to compute the SCF reference and summarize `E_SCF`, `E_VB`, `E_corr(VB)`, `E_best_template`, `E_res`, captured-subspace weight, and weight sums.
2. Run the fixed-orbital VBCI baseline for H2, allyl cation, allyl radical, allyl anion, and benzene first.
3. Add VBSCF rows only after the fixed-orbital VBCI row is stable for the same system.
4. Add BOVB rows only after the corresponding VBCI/VBSCF row is stable and the fixed-orbital zero-amplitude baseline is printed.
5. Print concise result rows that can be copied into this document.
6. Keep larger scans or pedagogical plots secondary to the stabilization table.

Current notebook update:

- Added a branch-stabilization table cell for fixed-orbital VBCI on H2, allyl cation, allyl radical, allyl anion, and benzene.
- The table computes an HF reference, `E_VB`, `E_corr(VB) = E_VB - E_SCF`, `E_best_template`, `E_res`, captured chemical-resonance weight where available, and structure/Lowdin weight sums.
- Added a dedicated H2 method-ladder cell for VBCI, VBSCF, and BOVB. It uses the same H2 geometry, basis, SCF reference, localized-template diagnostics, and monotonic energy checks `E(VBSCF) <= E(VBCI)` and `E(BOVB) <= E(VBSCF)`.
- Added an allyl method-ladder cell. It tests allyl cation VBCI, compact-CSF, and compact-CSF-BOVB in a split-valence basis. It now also tests allyl radical/anion VBCI, fixed-limit VBSCF, fixed-limit BOVB, and determinant-template compact-CSF routes.
- The allyl method-ladder cell now includes real common-orbital VBSCF rows for cation, radical, and anion in 6-31G. After rebuild, rerun this cell before any benzene work; it should verify that VBSCF is no higher than VBCI and that compact-CSF is no lower than its same-orbital internal full reference.
- The allyl method-ladder cell now uses same-call/internal reference checks for fixed-limit and compact rows so separate driver-call gauge changes do not create false failures.
- Added a benzene fixed-orbital compact-CSF gate. It first validates benzene VBCI and the Kekule/Dewar resonance projection before any benzene orbital relaxation is attempted.
- Added a benzene VBSCF common-orbital relaxation gate in 6-31G. It should be run after the benzene fixed-orbital compact-CSF gate and before any determinant-space BOVB implementation.
- Added a benzene VBSCF amplitude-bound scan. If every scanned bound remains active and the compact captured weight remains collapsed, the next driver task is a symmetry-preserving or trust-region benzene VBSCF parameterization, not BOVB yet.
- Added a benzene symmetry-preserving VBSCF scan using `orbital_relaxation_symmetry="equivalent-centers"`. This keeps all six center-local pi breathing amplitudes equal and is the next checkpoint to rerun after rebuild.

### Development order for this phase

1. Lock down fixed-orbital VBCI diagnostics for H2, allyl cation, allyl radical, allyl anion, and benzene.
2. Add SCF comparison and resonance-energy reporting to the notebook and then to `VbDriver` diagnostics where appropriate.
3. Stabilize real VBSCF on H2, with reproducible orbital optimization and an energy no higher than fixed-orbital VBCI unless constrained.
4. Extend VBSCF to allyl cation, then allyl radical, allyl anion, and benzene.
5. Treat BOVB as a separate method path, not as an ambiguous VBSCF option.
6. Validate BOVB first on H2, then allyl cation, then allyl radical, allyl anion, and benzene.
7. Only after the full table is green, resume general sigma/pi separation and frozen-versus-relaxed inactive-orbital options.

### Post-stabilization generalization strategy

After the stabilization ladder is green, generalize in controlled layers:

1. Keep pi-only active spaces as the validated core.
2. Use `OrbitalAnalyzer` only to identify sigma and pi candidate classes and atom assignments; `VbDriver` remains responsible for the active/inactive/frozen partition and wavefunction model.
3. Add explicit frozen-HF sigma/core embedding mode first.
4. Add optional relaxed inactive orbitals only after frozen-HF results remain reproducible.
5. Add structure-specific BOVB active-orbital relaxation before relaxing inactive sigma/core orbitals.
6. Expose explicit user options for pi-only active spaces, mixed sigma/pi active spaces, frozen HF sigma/core orbitals, and relaxed inactive orbitals.
7. Keep VBB as a later Linares et al. Valence Bond BOND method layer, separate from generic active-space construction and separate from BOVB.

## Metal-ligand scan checkpoint: 2026-05-04

The current real Pd--NH3/Pd--PH3 notebook state is deliberately conservative:

- The validated dissociation curves are the B3LYP/def2-SVP constrained-scan energies and the HF single-point total energies evaluated on those geometries.
- The curves are plotted as `E(R) - E(5.0 Angstrom)` in kJ/mol, so the potential well is negative and the 5.0 Angstrom point is the common zero.
- The current grid has internal reference minima, not boundary minima: Pd--NH3 has B3LYP/HF minima at about 2.21/2.49 Angstrom, and Pd--PH3 has B3LYP/HF minima at about 2.33/2.33 Angstrom.
- The corresponding grid dissociation energies are about 96.4/24.5 kJ mol^-1 for Pd--NH3 at B3LYP/HF and about 141.5/52.4 kJ mol^-1 for Pd--PH3 at B3LYP/HF.
- State-tracked sigma-only VB-SCF/BOVB metal-ligand values are still printed as diagnostics, but they are not reported as validated dissociation energies because the present metal-ligand VB total-energy model does not yet give stable production potential curves.
- The next VB task is not another notebook plot. It is a method fix: a size-consistent metal-ligand total-energy formulation for sigma donation plus pi back-donation, followed by regression tests that reject branch switches, ECP Hamiltonian mismatches, and boundary-minimum scans.
- `VbDriver` owns wavefunction algebra; `OrbitalAnalyzer` only supplies recognition and candidate metadata.

## Implemented capabilities

### H₂ one-active-bond VB-CI

Status: implemented and validated as the two-electron/two-orbital baseline.

Features:

- Automatic generated active space for H₂ when no explicit structures/orbitals are provided.
- Covalent plus two ionic structures.
- Spin-adapted two-electron singlet algebra.
- Metric-pruned generalized eigenproblem.
- Chirgwin-Coulson and Löwdin weights.
- AO ERIs through `FockDriver.compute_eri(..., eri_thresh=1.0e-12)`.
- H₂ dissociation notebook comparison to RHF and broken-symmetry UHF.

Important interpretation:

- Chirgwin-Coulson weights are useful but can become negative or larger than one in a strongly nonorthogonal structure basis.
- Löwdin weights are better for bounded visual diagnostics.

### Ethylene one-active-π VB-CI

Status: implemented and validated as the first polyatomic fixed-orbital checkpoint.

The current ethylene model is still a **2-electron / 2-orbital** active problem:

- active candidate: analyzer C=C `BD(pi)`, e.g. `BD_pi_4`,
- inactive sigma framework: C=C sigma plus four C-H sigma candidates,
- active orbitals: carbon-projected π components of the analyzer candidate,
- structures: π covalent plus two π ionic structures,
- embedding: inactive/core/sigma environment kept at the RHF reference level.

The embedding model is:

```text
inactive density = total RHF density - selected active pair density
```

The inactive density contributes an effective one-electron field and a constant frozen-HF energy. This keeps the sigma framework at the HF level while the selected π pair is treated by the VB-CI structure algebra.

User-verified ethylene checkpoint:

- active-space model: `one-active-pi-bond`
- active candidate: `BD_pi_4` / `pi`
- inactive sigma candidates: `BD_sigma_3`, `BD_sigma_7`, `BD_sigma_9`, `BD_sigma_11`, `BD_sigma_13`
- frozen embedding: `True`
- frozen electrons: about `14`
- active reference electrons: about `2`
- structure weights are finite and symmetric for the two equivalent ionic π structures.

## Why orbital relaxation is not the immediate next step

BOVB introduces structure-specific orbital relaxation. It should now follow the fixed-orbital π validation work, but only after the determinant/CSF boundary and compact resonance diagnostics remain covered by source-level tests.

Before production BOVB, the driver must be able to:

1. select a multi-center π system,
2. construct more than two active π orbitals,
3. support more than two active electrons,
4. generate chemically meaningful π VB structures,
5. handle closed-shell and open-shell spin cases,
6. keep sigma/core/inactive orbitals frozen at the HF reference level,
7. report stable weights and diagnostics.

Those fixed-orbital capabilities are now in place through benzene. The next bridge is a compact spin-adapted CSF Hamiltonian layer, followed by carefully constrained BOVB and VBB method work.

## Fixed-orbital multi-center π engine

### Scope

The implemented ladder generalizes from the ethylene 2e/2π model to fixed-orbital π active spaces of the form:

```text
n_active_electrons / n_active_pi_orbitals
```

The current implementation remains fixed-orbital and embedded in the frozen HF sigma/core density. Orbital optimization and BOVB come later.

### Required engine extensions

1. **Active π orbital selection**
   - Use analyzer `BD(pi)`, pi-capable one-center candidates, and NAO atom/angular maps to identify the π-active atoms.
   - Build one localized π active orbital per active atom or per selected π center.
   - Preserve traceability from each active orbital to atom index, NAO indices, and analyzer candidate labels.

2. **Flexible active electron count**
   - Ethylene: 2e/2π.
   - Allyl cation: 2e/3π.
   - Allyl radical: 3e/3π.
   - Allyl anion: 4e/3π.
   - Butadiene: 4e/4π.
   - Benzene: 6e/6π.

3. **General active determinant or CSF basis**
   - Start with determinant-based active-space CI for robustness.
   - Add spin-adapted CSFs once determinant algebra is stable.
   - Support singlets first, then doublets for radicals.

4. **General overlap and Hamiltonian construction**
   - Remove the current restriction to `(1,1)`, `(2,0)`, and `(0,2)` two-orbital structures.
   - Support arbitrary active occupation strings over non-orthogonal active orbitals.
   - Keep metric pruning and Löwdin weights.

5. **Frozen HF sigma/core embedding**
   - Continue using the total HF density minus active reference density as the inactive environment.
   - Report frozen electron count, active electron count, frozen constant energy, and active candidate labels.

6. **VB structure generation**
   - Generate chemically labeled π structures from graph/connectivity and active electron count.
   - Keep the first generator conservative; add richer ionic and long-bond structures later.

## Validation ladder

### 1. Allyl cation: 2e/3π singlet

This is the recommended next implementation target.

Why first:

- closed-shell singlet,
- only two active π electrons,
- three π centers,
- tests multi-center active orbital construction without open-shell complications,
- tests sigma/core frozen embedding beyond ethylene.

Expected model:

- active atoms: the three allyl carbon centers,
- active orbitals: one π orbital on each carbon,
- inactive framework: C-C sigma and C-H sigma candidates plus core candidates,
- active electron count: `2`,
- spin: singlet,
- minimal structures: left and right π-bond resonance forms, with optional ionic/polar extensions later.

Acceptance checks:

- three active π orbitals are produced,
- frozen electron count equals total electrons minus two,
- sigma candidates remain inactive/frozen,
- resonance structures are symmetry-equivalent for symmetric allyl cation,
- weights are normalized and finite,
- overlap metric rank is stable.

### 2. Allyl radical: 3e/3π doublet

Why second:

- adds open-shell spin handling,
- requires one unpaired π electron,
- tests `SOMO`/spin-density information from `OrbitalAnalyzer`.

Expected model:

- active electron count: `3`,
- spin: doublet,
- resonance forms with the radical center delocalized over terminal carbons,
- frozen sigma/core environment.

Acceptance checks:

- doublet spin sector is selected correctly,
- radical populations are chemically sensible,
- left/right radical resonance weights are symmetric in symmetric geometry.

### 3. Allyl anion: 4e/3π singlet

Why third:

- tests a filled 3-center π system with excess electron density,
- requires lone-pair-like/ionic π structures,
- remains closed-shell but has more active electrons than ethylene or allyl cation.

Expected model:

- active electron count: `4`,
- spin: singlet,
- resonance forms with negative charge delocalized over terminal carbons,
- optional ionic/lone-pair-like π structures.

Acceptance checks:

- active electron count is four,
- charge-delocalized structures are generated and labeled,
- terminal symmetry is respected.

### 4. Butadiene: 4e/4π singlet

Why before benzene:

- bridge from 3-center allyl to cyclic/aromatic systems,
- tests a linear conjugated system with two π bonds,
- simpler than benzene while requiring more than one π pair.

Expected model:

- active electron count: `4`,
- active orbitals: four carbon π orbitals,
- structures: two covalent π bonds plus charge-separated/resonance alternatives.

Acceptance checks:

- four active π orbitals are found,
- sigma framework remains frozen,
- dominant structure is chemically sensible.

### 5. Benzene: 6e/6π singlet

Why after butadiene:

- tests cyclic conjugation and aromatic resonance,
- requires Kekulé and Dewar-style structure generation,
- larger overlap/Hamiltonian matrix.

Expected model:

- active electron count: `6`,
- active orbitals: six carbon π orbitals,
- structures: at least two Kekulé structures, then Dewar structures,
- sigma framework frozen.

Acceptance checks:

- two Kekulé structures are degenerate by symmetry,
- Dewar structures are present and have lower weight in first fixed-orbital model,
- structure weights are normalized and stable,
- no BOVB or orbital relaxation is needed to pass the fixed-orbital checkpoint.

### User-validated fixed-orbital π ladder output

The notebook has now run the fixed-orbital π ladder through benzene with frozen-HF sigma/core embedding:

| System | Active space | Spin | Active orbitals | Structures/determinants | Active reference electrons | Frozen electrons | Status |
| --- | --- | --- | ---: | ---: | ---: | ---: | --- |
| allyl cation | 2e/3π | singlet | 3 | 6 spin-adapted structures | 2.0 | 20.0 | validated |
| allyl radical | 3e/3π | doublet | 3 | 9 determinants | 3.0 | 20.0 | validated |
| allyl anion | 4e/3π | singlet | 3 | 9 determinants | 4.0 | 20.0 | validated |
| butadiene | 4e/4π | singlet | 4 | 36 determinants | 4.0 | 26.0 | validated |
| benzene | 6e/6π | singlet | 6 | 400 determinants | 6.0 | 36.0 | validated |

Interpretation:

- The active/inactive electron partition is stable for all current checkpoints.
- The 2e/3π allyl cation remains the chemically transparent spin-adapted VB benchmark.
- Larger and open-shell systems currently use an orthonormal determinant-CI fallback. This is a robust fixed-orbital active-space checkpoint, but the determinants are not yet compact chemically labeled VB resonance structures.
- Determinant labels now include chemically readable spatial occupations, such as radical center, lone-pair/doubly occupied center, empty center, and open-shell singlet center lists.
- Determinant-CI diagnostics also report grouped resonance-structure weights by spatial occupation. This is the first interpretability layer before true spin-adapted CSF Hamiltonians.
- A second interpretability layer now generates graph-based spin-adapted chemical resonance templates and projects the determinant-CI wavefunction onto them. The current template set covers allyl radical centers, allyl anion lone-pair/π-bond forms, butadiene paired π-bond forms, and benzene Kekulé/Dewar-style pairings.
- For singlet determinant-CI active spaces, the driver now selects the lowest alpha/beta exchange-symmetric root before building resonance diagnostics. This avoids projecting spin-adapted CSF templates onto the exchange-antisymmetric M_S = 0 root that appeared first for butadiene in the determinant-only spectrum.
- CSF template phases are now generated directly in the determinant-CI alpha/beta occupation-string representation instead of applying second-quantized creation-order signs a second time. This fixes the butadiene paired π-bond template capture and keeps the displayed resonance weights compact enough for validation.
- Chemical-resonance diagnostics now separate the compact CSF-template subspace capture from the normalized interpretation within that subspace: `chemical_resonance_subspace_weight` reports how much of the determinant-CI root is captured by the generated templates, while `chemical_resonance_weights` are normalized Löwdin-style weights over the retained template metric. These reported weights are graph-automorphism averaged for the current allyl, butadiene, and benzene template graphs so symmetry-equivalent chemical resonance structures carry equal displayed weights; the raw template-metric values remain available as `chemical_resonance_unsymmetrized_weights` to diagnose orbital/reference symmetry breaking.
- The next scientific improvement should promote these diagnostic CSF projections into an optional compact spin-adapted CSF Hamiltonian where chemically useful, while retaining the determinant-CI fallback as a validation reference.

## Future improvements: compact CSFs, BOVB, and VBB

The fixed-orbital π engine is now stable enough to define the next VB roadmap. The next work should proceed in layers so that each new orbital-relaxation or user-facing feature can be checked against the existing determinant-CI and CSF-projection diagnostics.

### 1. Compact spin-adapted CSF Hamiltonians

Goal: promote the current diagnostic CSF projections into an optional working Hamiltonian basis where chemically useful.

Current status: implemented for allyl cation 2e/3π singlet as `mode="compact-csf"`. The compact basis contains the two adjacent π-bond CSFs, solves its own generalized Hamiltonian, and reports the full fixed-orbital VB-CI reference energy plus the captured-subspace weight.

Implementation steps:

1. Reuse the graph-generated CSF templates for allyl, butadiene, and benzene. **Status: allyl cation adjacent π-bond templates are implemented as a working compact basis.**
2. Build the Hamiltonian and overlap matrices directly in the retained CSF-template metric. **Status: implemented for the two-electron singlet multicenter π path.**
3. Compare CSF-Hamiltonian roots against determinant-CI or full fixed-orbital VB-CI roots for the same fixed active orbitals. **Status: allyl cation reports the full fixed-orbital VB-CI oracle.**
4. Keep determinant-CI as the reference fallback and as a regression oracle.
5. Report captured-subspace weight whenever a compact CSF basis is used so users can tell whether the compact basis is faithful.

Acceptance checks:

- Butadiene paired π-bond CSFs retain nonzero captured-subspace weight and exchange-symmetric singlet parity.
- Benzene Kekulé/Dewar template counts and symmetry-orbit weights remain stable.
- Compact CSF energies are explicitly labeled as compact-basis approximations when they do not span the full determinant-CI space.

### 2. BOVB: breathing-orbital valence bond

BOVB should be introduced after the fixed-orbital multi-center π engine and the determinant/CSF regression tests remain stable for allyl, butadiene, and benzene checkpoints.

Recommended BOVB order:

1. Revisit ethylene with two π structures and structure-specific breathing orbitals. **Status: implemented for the current covalent plus ionic one-active-π active space.**
2. Apply BOVB to allyl cation after fixed-orbital 2e/3π works and after the compact CSF/structure-specific orbital layer is ready. **Status: implemented for `mode="compact-csf-bovb"` on allyl-cation-like 2e/3π singlet active spaces.**
3. Extend to butadiene and benzene only after the active-space generator and spin handling are stable.

BOVB requirements:

- structure-specific orbital sets,
- constraints to preserve sigma/core frozen orbitals,
- orbital normalization and nonorthogonality controls,
- stable orbital optimization with line search/trust region,
- diagnostics comparing fixed-orbital VB-CI, VB-SCF, and BOVB.

Implementation steps:

1. Start with ethylene and H₂-like two-orbital cases, where the fixed-orbital baseline is already transparent. **Status: implemented for H₂ and ethylene one-active-π VB active spaces.**
2. Represent each VB structure with its own allowed active-orbital set while sharing frozen core/sigma orbitals where requested. **Status: implemented for the two-orbital/two-electron singlet path, including center-local external breathing directions when the AO basis supplies them.**
3. Add constraints that preserve active-center identity and prevent collapse into unrestricted HF-like orbitals.
4. Optimize orbital parameters and structure coefficients with a robust alternating or joint optimizer.
5. Add diagnostics for structure-specific orbital overlaps, breathing magnitude, orbital gradients, line-search progress, and energy lowering relative to fixed-orbital VB-CI.
6. Only then extend to allyl and butadiene π systems.

Acceptance checks:

- BOVB energy is no higher than the corresponding fixed-orbital VB-CI energy for the same structure set.
- Structure-specific orbitals remain normalized and chemically attached to their intended centers.
- Frozen sigma/core electron count remains stable.
- BOVB output clearly distinguishes fixed-orbital weights, optimized-orbital weights, and any compact CSF projection diagnostics.

### 3. VBB: Valence Bond BOND method

VBB refers to **Valence Bond BOND**, where the second `BOND` stands for **Breathing Orbital Naturally Delocalized**. This is the Linares et al. VBB method and should be documented as a VB electronic-structure method, not as a generic builder/block layer.

The name should stay separate from BOVB:

- BOVB relaxes structure-specific orbitals in breathing-orbital valence bond models.
- VBB uses breathing orbitals that are naturally delocalized according to the Linares et al. method definition.
- A future user-facing constructor/helper layer may still be useful, but it should not be called VBB.

Planned VBB work:

1. **Method definition and notation**
   - Add a concise implementation note for the Linares et al. VBB formalism.
   - Define the relation between VBB structures, breathing orbitals, and naturally delocalized orbital constraints.
   - Keep VBB terminology distinct from generic VB input builders.

2. **VBB orbital model**
   - Introduce naturally delocalized breathing-orbital parameterizations for selected VB structures.
   - Preserve active-center identity, spin coupling, and frozen sigma/core constraints.
   - Compare the allowed VBB orbital space to the corresponding BOVB and fixed-orbital VB spaces.

3. **VBB Hamiltonian and optimization**
   - Reuse the compact CSF or determinant reference machinery for fixed-orbital checks.
   - Add VBB-specific orbital-optimization variables and constraints only after BOVB/fixed-orbital baselines are stable.
   - Report orbital-delocalization diagnostics, breathing amplitudes, structure weights, and energy lowering relative to fixed-orbital VB-CI and BOVB where applicable.

4. **Validation and reporting**
   - Start from two-orbital H₂/ethylene-style benchmarks before applying VBB to allyl, butadiene, and benzene.
   - Verify that VBB energies and diagnostics are labeled separately from BOVB and compact CSF results.
   - Refuse ambiguous VBB specifications unless the Linares et al. constraints can be applied deterministically.

Acceptance checks:

- H₂ and ethylene VBB prototypes reproduce the expected Linares et al. limiting behavior for the chosen benchmark definitions.
- Allyl, butadiene, and benzene VBB extensions remain traceable to the current fixed-orbital π ladder and compact CSF templates.
- Output clearly labels fixed-orbital VB-CI, compact CSF, BOVB, and VBB quantities as distinct models.

### 4. Metal-ligand VB active spaces

`OrbitalAnalyzer` now exposes metal-ligand candidate records and `VbDriver` supports optional fixed-orbital VB active spaces involving ligand donation and back-donation.

Current status:

- `active_metal_ligand_channels=("sigma-acceptor",)` builds a sigma-only donor/acceptor active space and reports `metal_ligand_model="sigma-only"` with back-donation blocked.
- `active_metal_ligand_channels=("sigma-acceptor", "pi-donor")` builds a combined sigma-plus-back-donation active space and reports `metal_ligand_model="sigma-plus-backdonation"` with back-donation enabled.
- The sigma-only model is a 2-orbital/2-electron determinant-CI active space for the ligand-to-metal sigma channel.
- The sigma-plus-back-donation model is a 4-orbital/4-electron determinant-CI active space with both the sigma donor/acceptor pair and the pi donor/acceptor pair active.
- Result diagnostics expose the selected metal-ligand records, active orbital count, active electron count, determinant count, requested channels, and the blocked/enabled back-donation switch.
- `mode="bovb"` with `active_metal_ligand_channels` keeps the same channel switches, optimizes center-local breathing directions for the active metal-ligand orbitals, and reports `metal_ligand_bovb_initial_energy`, `metal_ligand_bovb_energy_lowering`, and per-orbital breathing amplitudes.
- The notebook `docs/metal_ligand_recognition.ipynb` now runs real Pd--NH3 and Pd--PH3 SCF/NBO calculations, reuses the real NAO payload, and compares sigma-only against sigma-plus-back-donation VB-CI and BOVB without mock payloads.

Planned steps:

- Treat metal-ligand candidates as analyzer-provided active-orbital suggestions, not as automatic selected VB structures.
- Add a pedagogical dissociation sequence for metal-ligand bonds analogous to H2: RHF, broken-symmetry UHF, fixed-orbital VB-CI, and BOVB/VBB correlation-orbital relaxation.
- Add metal-ligand VBB benchmarks only after the organic VBB method layer is defined and validated.
- Keep determinant-CI fallback available because compact spin coupling for metal centers may be system dependent.
- Add diagnostics separating ligand-field interpretation from organic π resonance interpretation.

## Immediate implementation checklist

Current fixed-orbital π implementation status:

- [x] Document current H₂ and ethylene status.
- [x] Define the fixed-orbital π active-space roadmap.
- [x] Add an `active_pi_atoms` selection option to `VbComputeOptions`.
- [x] Build a multi-center `VbActiveSpace` representation with active electron count, spin, and arbitrary active π orbitals.
- [x] Implement active π orbital construction for more than two atoms.
- [x] Keep the 2e/3π allyl cation as a spin-adapted two-electron VB checkpoint.
- [x] Add an orthonormal determinant-CI fallback for larger/open-shell fixed-orbital π spaces.
- [x] Add chemically readable determinant labels and grouped spatial-occupation resonance weights.
- [x] Add graph-based spin-adapted resonance/CSF projection diagnostics for allyl radical, allyl anion, butadiene, and benzene.
- [x] Extend frozen-HF embedding from one selected pair density to a general active reference density.
- [x] Add notebook cells for allyl cation and the validation ladder through benzene.
- [x] Add source-level smoke tests for allyl cation, allyl radical, allyl anion, butadiene, and benzene.
- [x] Add source-level regression tests for resonance/CSF labels, counts, projection weights, butadiene captured-subspace weight, and singlet exchange parity.
- [x] Add source-level H2 gate comparing stretched RHF, broken-symmetry UHF, fixed-orbital VB-CI, and H2 VB-SCF before extending the workflow to metal-ligand systems.
- [x] Add a narrow H2 BOVB prototype with structure-specific breathing orbitals and regression coverage against the H2 RHF/UHF/VB-CI/VBSCF gate.
- [x] Generalize the BOVB route from H₂-only to two-orbital/two-electron singlet active spaces and add ethylene one-active-π regression coverage with frozen sigma/core embedding.
- [x] Add center-local external breathing directions to two-orbital BOVB and validate that H₂/6-31G lowers below the fixed-orbital limit while minimal-basis H₂ remains unchanged.
- [x] Stabilize automatic H₂ dissociation scans by forcing the default H₂ active space to an atom-centered H-H pair across all bond distances.
- [x] Add allyl compact-CSF BOVB and validate that the same-orbital zero-amplitude compact limit lowers in a split-valence basis.
- [x] Add source-level regression tests for metal-ligand sigma-only versus sigma-plus-back-donation active-space metadata.
- [x] Add metal-ligand BOVB for sigma-only and sigma-plus-back-donation determinant-CI channel spaces with zero-amplitude fixed-limit diagnostics.
- [x] Add compact spin-adapted CSF Hamiltonian mode for allyl cation 2e/3π with the full fixed-orbital VB-CI result as oracle.
- [ ] Extend compact spin-adapted CSF Hamiltonians to selected determinant-CI diagnostics where chemically useful.
- [x] Generalize BOVB beyond two-orbital/two-electron singlet active spaces for analyzer-selected metal-ligand determinant-CI channel models.
- [ ] Generalize BOVB beyond two-orbital/two-electron singlet active spaces for larger organic π determinant-CI models.
- [ ] Add VBB method notes and prototypes for H₂/ethylene-style benchmarks, then extend to allyl, butadiene, and benzene.
- [x] Add fixed-orbital metal-ligand active-space support after analyzer metal-ligand recognition exists.
- [ ] Add H2 and metal-ligand dissociation notebooks comparing RHF, UHF, fixed-orbital VB-CI, and future BOVB correlation/orbital-relaxation effects.

## Tomorrow restart notes

Start from the validated fixed-orbital baseline, not from BOVB/VBB:

1. Keep H2 as the gatekeeper before metal-ligand VB: stretched RHF must fail relative to broken-symmetry UHF, fixed-orbital VB-CI must track the UHF dissociation limit, and H2 VB-SCF/BOVB must remain no higher than fixed-orbital VB-CI.
2. Treat BOVB as implemented for two-orbital/two-electron singlet active spaces, the allyl compact-CSF checkpoint, and metal-ligand determinant-CI channel spaces. The H₂, allyl, and metal-ligand split-valence benchmarks now exercise real center-local breathing relaxation; larger organic π systems still need their own orbital-optimization layer and regression tests.
3. Re-open `docs/metal_ligand_recognition.ipynb` only after the H2 VB-CI/VB-SCF/BOVB ladder is clean enough to use as the organic reference.
4. The metal-ligand entry points to remember are `active_metal_ligand_channels=("sigma-acceptor",)` for the sigma-only model and `active_metal_ligand_channels=("sigma-acceptor", "pi-donor")` for the combined sigma/back-donation model.
5. The best next implementation task is a Pd--ligand distance checkpoint comparing RHF, UHF, fixed VB-CI, and ML-BOVB for sigma-only versus sigma-plus-backdonation.
6. Do not claim VBB or production metal-ligand VB is implemented. The current metal-ligand BOVB model is a determinant-CI channel-relaxation prototype seeded by analyzer diagnostics.

Current safe stopping point: fixed-orbital H2, ethylene, π-ladder, allyl cation compact CSF/BOVB, metal-ligand sigma-only/sigma-plus-back-donation VB-CI/BOVB diagnostics, and two-orbital H₂/ethylene BOVB are in place; larger organic multicenter orbital relaxation, VBB, and production coordination VB remain roadmap items.

## Documentation and notebook policy

- The notebook is for validation and interpretation, not for core logic.
- All active-space construction and Hamiltonian logic should live in `VbDriver`.
- Notebook cells should print concise diagnostics:
  - active model,
  - active electron/orbital counts,
  - active atoms and labels,
  - inactive sigma/core labels,
  - frozen electron count,
  - energy,
  - overlap rank/eigenvalues,
  - weights by structure.

## Current short-term decision

The initial fixed-orbital multi-center π engine is now in place through benzene, the first non-H₂ BOVB step is in place for ethylene-style two-orbital active spaces, and allyl cation now has a compact CSF Hamiltonian benchmark. The immediate next VB step is to extend compact CSF Hamiltonians to selected determinant-CI template spaces; multicenter BOVB and VBB should follow only with explicit regression tests against the fixed-orbital baseline.
