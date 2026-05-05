from dataclasses import dataclass, field
from typing import Any, List, Optional, Dict, Tuple
import numpy as np
from .orbitalanalyzerdriver import OrbitalAnalyzer, OrbitalAnalyzerOptions


# Data classes for VB driver
@dataclass(frozen=True)
class VbOrbital:
	"""
	Represents a single (possibly non-orthogonal) orbital for VB calculations.
	label: str - human-readable label for the orbital
	coefficients: np.ndarray - AO or NAO coefficients
	center: Optional[int] - atom index or None
	kind: str - e.g. 'active', 'inactive', 'core', etc.
	"""
	label: str
	coefficients: np.ndarray
	center: Optional[int] = None
	kind: str = "active"

@dataclass(frozen=True)
class VbStructure:
	"""
	Represents a single VB structure (occupation pattern over orbitals).
	label: str - human-readable label
	occupation: Tuple[Any, ...] - occupation pattern (to be specified)
	spin: str - e.g. 'singlet', 'doublet'
	charge_pattern: Dict[int, int] - atom index to formal charge
	"""
	label: str
	occupation: Tuple[Any, ...]
	spin: str = "singlet"
	charge_pattern: Dict[int, int] = field(default_factory=dict)


@dataclass(frozen=True)
class VbActiveSpace:
	"""Internal active-space description generated from OrbitalAnalyzer records."""

	active_bond: Tuple[int, ...]
	active_candidate_label: str
	active_candidate: Dict[str, Any]
	active_orbitals: Tuple[VbOrbital, ...]
	structures: Tuple[VbStructure, ...]
	frozen_candidates: Tuple[Dict[str, Any], ...] = field(default_factory=tuple)
	inactive_candidates: Tuple[Dict[str, Any], ...] = field(default_factory=tuple)
	excluded_candidates: Tuple[Dict[str, Any], ...] = field(default_factory=tuple)
	metadata: Dict[str, Any] = field(default_factory=dict)


# Options class for VB driver
@dataclass(frozen=True)
class VbComputeOptions:
	"""
	Options for controlling the VB driver.
	mode: str - calculation mode (e.g. 'vbci', 'vbscf', 'bovb', 'compact-csf', 'compact-csf-bovb')
	max_iter: int - maximum optimization iterations
	conv_thresh: float - convergence threshold
	optimize_orbitals: bool - whether to optimize orbitals
	include_ionic: bool - include ionic structures
	include_bovb: bool - enable BOVB (breathing orbital VB)
	active_bond: optional atom-index pair for one-active-bond generation
	active_candidate_label: optional OrbitalAnalyzer candidate label to activate
	active_candidate_subtype: optional bond subtype filter, currently 'sigma' or 'pi'
	active_metal_ligand_channels: optional tuple of metal-ligand channels to activate
	active_metal_ligand_reference_orbitals: optional bound-geometry orbitals used to state-track metal-ligand scans
	active_bond_reference_orbitals: optional bound-geometry active orbitals used to state-track organic bond scans
	active_pi_atoms: optional atom-index tuple for a fixed-orbital multi-center pi active space
	active_electron_count: optional number of active electrons for generated pi active spaces
	active_spin: optional spin label for generated pi active spaces, currently 'singlet'
	freeze_inactive_orbitals: embed the active pair in the frozen RHF density
	use_active_space: whether automatic generated active spaces are enabled
	"""
	mode: str = "vbscf"
	max_iter: int = 50
	conv_thresh: float = 1.0e-8
	optimize_orbitals: bool = True
	include_ionic: bool = True
	include_bovb: bool = False
	active_bond: Optional[Tuple[int, int]] = None
	active_candidate_label: Optional[str] = None
	active_candidate_subtype: Optional[str] = None
	active_metal_ligand_channels: Optional[Tuple[str, ...]] = None
	active_metal_ligand_reference_orbitals: Optional[Tuple[Any, ...]] = None
	active_bond_reference_orbitals: Optional[Tuple[Any, ...]] = None
	active_pi_atoms: Optional[Tuple[int, ...]] = None
	active_electron_count: Optional[int] = None
	active_spin: str = "singlet"
	freeze_inactive_orbitals: bool = True
	use_active_space: bool = True


# Main VB driver class
class VbDriver:
	def compute(self,
		molecule: Any,
		basis: Any,
		structures: Optional[List['VbStructure']] = None,
		orbitals: Optional[List['VbOrbital']] = None,
		reference_orbitals: Optional[Any] = None,
		options: Optional[VbComputeOptions] = None,
		nao_data: Optional[Any] = None
	) -> Dict[str, Any]:
		"""
		Compute the VB energy and related quantities for arbitrary molecules.

		This version uses the shared OrbitalAnalyzer payload used by the NBO driver.
		"""
		user_supplied_orbitals = orbitals is not None
		user_supplied_structures = structures is not None
		options = options or VbComputeOptions()
		n_ao = basis.get_dimensions_of_basis()
		if orbitals is None and self._is_h2_three_structure_case(molecule, structures):
			ao_to_atom, _ = OrbitalAnalyzer.ao_shell_map(molecule, basis)
			orb1_coeffs = np.zeros(n_ao)
			orb2_coeffs = np.zeros(n_ao)
			orb1_coeffs[np.where(ao_to_atom == 0)[0]] = 1.0
			orb2_coeffs[np.where(ao_to_atom == 1)[0]] = 1.0
			orb1_coeffs /= np.linalg.norm(orb1_coeffs)
			orb2_coeffs /= np.linalg.norm(orb2_coeffs)
			orbitals = [
				VbOrbital(label='1s_A', coefficients=orb1_coeffs, center=0),
				VbOrbital(label='1s_B', coefficients=orb2_coeffs, center=1),
			]
			if self._is_minimal_h2_vbci_case(molecule, structures, orbitals):
				if options is not None and (bool(getattr(options, 'include_bovb', False)) or getattr(options, 'mode', 'vbscf') == 'bovb'):
					return self.compute_bovb_two_orbital(molecule, basis, structures, orbitals, reference_orbitals, options)
				if options is not None and getattr(options, 'mode', 'vbscf') == 'vbscf' and getattr(options, 'optimize_orbitals', True):
					return self.compute_vbscf_h2(molecule, basis, structures, orbitals, reference_orbitals, options)
				return self.compute_vbci_h2(molecule, basis, structures, orbitals, reference_orbitals, options)

		analysis = None
		if nao_data is None:
			try:
				import veloxchem as vlx
				multiplicity = int(round(float(molecule.get_multiplicity())))
				scf_drv = vlx.ScfUnrestrictedDriver() if multiplicity != 1 else vlx.ScfRestrictedDriver()
				scf_results = scf_drv.compute(molecule, basis)
				mol_orbs = reference_orbitals if reference_orbitals is not None else scf_results.get('C', scf_drv.mol_orbs)
				analysis = OrbitalAnalyzer(
					molecule,
					basis,
					mol_orbs=mol_orbs,
					options=OrbitalAnalyzerOptions(include_mo_analysis=False),
				)
				analysis = analysis.run()
				nao_data = analysis.nao_data
			except Exception as e:
				raise RuntimeError("VB driver could not compute orbital analysis automatically. Please provide orbitals or nao_data explicitly.\nReason: " + str(e))
			candidates = analysis.orbital_candidates
		else:
			if reference_orbitals is not None:
				analysis = OrbitalAnalyzer(
					molecule,
					basis,
					mol_orbs=reference_orbitals,
					options=OrbitalAnalyzerOptions(include_mo_analysis=False),
				)
				analysis = analysis.run()
				nao_data = analysis.nao_data
			else:
				analysis = OrbitalAnalyzer(molecule, basis).classify_nao_data(nao_data)
			candidates = analysis.orbital_candidates

		# Build VbOrbital objects for each classified orbital
		vb_orbitals = [
			VbOrbital(
				label=f"{c['type']}_{c['subtype']}_{c['index']}",
				coefficients=self._candidate_coefficient_vector(c, n_ao),
				center=c['atoms'][0] if 'atoms' in c and len(c['atoms']) > 0 else None,
				kind=c['type'],
			)
			for c in candidates
		]
		if orbitals is None:
			orbitals = vb_orbitals

		if (self._should_use_one_active_bond_space(molecule,
													  options,
													  user_supplied_structures,
													  user_supplied_orbitals)):
			active_space = self._build_one_active_bond_space(
				molecule,
				basis,
				candidates,
				nao_data,
				options,
			)
			if bool(getattr(options, 'include_bovb', False)) or options.mode in ('bovb', 'metal-ligand-bovb', 'ml-bovb'):
				if active_space.metadata.get("determinant_ci", False):
					if active_space.metadata.get("metal_ligand_model") is None:
						raise NotImplementedError(
							"BOVB determinant-CI relaxation is currently implemented only for metal-ligand active spaces.")
					active_result = self._compute_metal_ligand_bovb(
						molecule,
						basis,
						list(active_space.structures),
						list(active_space.active_orbitals),
						options,
						active_space=active_space,
						analysis=analysis,
						freeze_inactive=bool(options.freeze_inactive_orbitals),
					)
				elif len(active_space.active_orbitals) != 2:
					raise NotImplementedError(
						"BOVB is currently implemented only for two-orbital two-electron singlet and metal-ligand determinant active spaces.")
				else:
					active_result = self.compute_bovb_two_orbital(
						molecule,
						basis,
						list(active_space.structures),
						list(active_space.active_orbitals),
						reference_orbitals,
						options,
						active_space=active_space,
						analysis=analysis,
						freeze_inactive=bool(options.freeze_inactive_orbitals),
					)
			elif options.mode in ('compact-csf-bovb', 'csf-bovb'):
				active_result = self._compute_compact_csf_bovb(
					molecule,
					basis,
					list(active_space.structures),
					list(active_space.active_orbitals),
					options,
					active_space=active_space,
					analysis=analysis,
					freeze_inactive=bool(options.freeze_inactive_orbitals),
				)
			elif options.mode in ('compact-csf', 'csf'):
				active_result = self._compute_compact_csf(
					molecule,
					basis,
					list(active_space.structures),
					list(active_space.active_orbitals),
					active_space=active_space,
					analysis=analysis,
					freeze_inactive=bool(options.freeze_inactive_orbitals),
				)
			elif options.mode == 'vbscf' and options.optimize_orbitals:
				if (active_space.metadata.get("determinant_ci", False) and
						active_space.metadata.get("metal_ligand_model") is not None):
					active_result = self._compute_metal_ligand_vbscf_determinant(
						molecule,
						basis,
						list(active_space.structures),
						list(active_space.active_orbitals),
						active_space=active_space,
						analysis=analysis,
						freeze_inactive=bool(options.freeze_inactive_orbitals),
					)
				else:
					if len(active_space.active_orbitals) != 2:
						raise RuntimeError("VB-SCF orbital optimization is currently only available for two active orbitals or metal-ligand determinant active spaces.")
					vbscf_active_space = active_space
					if active_space.metadata.get("metal_ligand_model") == "sigma-only":
						vbscf_active_space = self._spin_adapted_metal_ligand_sigma_space(
							active_space,
							options,
						)
					active_result = self.compute_vbscf_h2(
						molecule,
						basis,
						list(vbscf_active_space.structures),
						list(vbscf_active_space.active_orbitals),
						reference_orbitals,
						options,
						active_space=vbscf_active_space,
						analysis=analysis,
						freeze_inactive=bool(options.freeze_inactive_orbitals),
					)
					active_space = vbscf_active_space
			else:
				if active_space.metadata.get("determinant_ci", False):
					active_result = self._compute_active_determinant_ci(
						molecule,
						basis,
						list(active_space.structures),
						list(active_space.active_orbitals),
						active_space=active_space,
						analysis=analysis,
						freeze_inactive=bool(options.freeze_inactive_orbitals),
					)
				else:
					active_result = self._compute_two_orbital_vbci(
						molecule,
						basis,
						list(active_space.structures),
						list(active_space.active_orbitals),
						active_space=active_space,
						analysis=analysis,
						freeze_inactive=bool(options.freeze_inactive_orbitals),
					)
			active_result["active_space"] = active_space
			active_result["orbital_analysis"] = analysis
			active_result.setdefault("diagnostics", {})
			method_message = active_result["diagnostics"].get("message")
			active_result["diagnostics"].update(
				self._active_space_diagnostics(active_space))
			active_result["diagnostics"]["vb_method"] = options.mode
			if method_message:
				active_result["diagnostics"]["active_space_message"] = (
					"Generated VB result from OrbitalAnalyzer active-space candidates.")
				active_result["diagnostics"]["message"] = method_message
			return active_result

		if structures is None:
			structures = [VbStructure(label='default', occupation=(tuple([1]*len(orbitals)),), spin='singlet')]

		if self._is_minimal_h2_vbci_case(molecule, structures, orbitals):
			if options is not None and (bool(getattr(options, 'include_bovb', False)) or getattr(options, 'mode', 'vbscf') == 'bovb'):
				return self.compute_bovb_two_orbital(molecule, basis, structures, orbitals, reference_orbitals, options)
			if options is not None and getattr(options, 'mode', 'vbscf') == 'vbscf' and getattr(options, 'optimize_orbitals', True):
				return self.compute_vbscf_h2(molecule, basis, structures, orbitals, reference_orbitals, options)
			else:
				return self.compute_vbci_h2(molecule, basis, structures, orbitals, reference_orbitals, options)

		# Expose all detailed orbital label information for notebook use
		diagnostics = {
			"message": "Generalized VB driver: orbitals classified using shared OrbitalAnalyzer payload.",
			"n_orbitals": len(vb_orbitals),
			"orbital_analysis_source": "OrbitalAnalyzer",
			"orbital_labels": [o.label for o in vb_orbitals],
			"orbital_types": [o.kind for o in vb_orbitals],
			"candidate_partitions": self._candidate_partition_diagnostics(candidates),
			"bond_candidate_summary": self._bond_candidate_summary(candidates),
			"suggested_active_candidates": self._suggest_active_candidates(candidates),
			"orbital_details": [
				{
					"label": o.label,
					"type": c.get("type"),
					"subtype": c.get("subtype"),
					"atoms": c.get("atoms"),
					"occupation": c.get("occupation"),
					"coefficients": c.get("coefficients"),
					"source": c.get("source"),
				}
				for o, c in zip(vb_orbitals, candidates)
			],
		}
		return {
			"energy": 0.0,
			"structure_coefficients": np.zeros(len(structures)),
			"overlap": np.eye(len(structures)),
			"Hamiltonian": np.eye(len(structures)),
			"weights": np.ones(len(structures)) / max(1, len(structures)),
			"orbitals": [orb.coefficients for orb in orbitals],
			"diagnostics": diagnostics,
			"orbital_analysis": analysis,
		}

	def compute_vbscf_h2(self,
					 molecule,
					 basis,
					 structures,
					 orbitals,
					 reference_orbitals,
					 options,
					 active_space=None,
					 analysis=None,
					 freeze_inactive=False):
		"""
		VB-SCF for H2: optimize two active orbitals by a mixing angle.
		"""
		from scipy.optimize import minimize_scalar

		H_ao, S_ao, eri_ao, e_nuc = self._ao_integrals(molecule, basis)
		embedding = None
		if freeze_inactive and active_space is not None and analysis is not None:
			embedding = self._frozen_hf_embedding(
				H_ao,
				S_ao,
				eri_ao,
				e_nuc,
				active_space,
				analysis,
			)
			if embedding is not None:
				H_ao = embedding["h_effective"]
				e_nuc = embedding["constant_energy"]
				orbitals = self._project_active_orbitals_out_of_frozen_density(
					orbitals,
					S_ao,
					embedding.get("frozen_density"),
				)

		# Assume orbitals[0] and orbitals[1] are initial guesses for 1s_A and 1s_B
		ao1 = orbitals[0].coefficients
		ao2 = orbitals[1].coefficients

		def rotated_orbitals(theta):
			c1 = np.cos(theta)
			c2 = np.sin(theta)
			orbA = c1 * ao1 + c2 * ao2
			orbB = -c2 * ao1 + c1 * ao2
			orbA = self._s_normalize(orbA, S_ao)
			orbB = self._s_normalize(orbB, S_ao)
			return orbA, orbB

		def energy_for_theta(theta):
			orbA, orbB = rotated_orbitals(theta)
			C = np.column_stack([orbA, orbB])
			try:
				S, H = self._build_two_orbital_singlet_matrices(
					structures, C, H_ao, S_ao, eri_ao, e_nuc)
				energy, _, _, _ = self._solve_generalized_vb(H, S)
			except Exception:
				energy = 1e6
			return energy

		res = minimize_scalar(energy_for_theta, bounds=(0, np.pi), method='bounded', options={'xatol': options.conv_thresh if options else 1e-8})
		theta_opt = res.x
		energy_opt = res.fun
		orbA, orbB = rotated_orbitals(theta_opt)
		C = np.column_stack([orbA, orbB])
		S, H = self._build_two_orbital_singlet_matrices(
			structures, C, H_ao, S_ao, eri_ao, e_nuc)
		energy, coeffs, weights, kept, lowdin_weights = self._solve_generalized_vb(
			H, S, return_lowdin=True)
		diagnostics = {
			"message": "H2 VB-SCF result (spin-adapted two-electron algebra)",
			"theta_opt": theta_opt,
			"overlap_condition": float(np.linalg.cond(S)),
			"overlap_eigenvalues": np.linalg.eigvalsh(S).tolist(),
			"retained_overlap_rank": int(len(kept)),
			"weight_scheme": "Chirgwin-Coulson",
			"available_weight_schemes": ["Chirgwin-Coulson", "Lowdin"],
			"optimizer": res,
		}
		if embedding is not None:
			diagnostics.update({
				"frozen_hf_embedding": True,
				"active_orbitals_orthogonalized_to_frozen_space": True,
				"frozen_electron_count": embedding["frozen_electron_count"],
				"active_reference_electron_count": embedding["active_electron_count"],
				"frozen_constant_energy": embedding["constant_energy"],
				"embedding_model": "inactive HF density = total RHF density minus active reference density",
			})
		return {
			"energy": float(energy),
			"structure_coefficients": coeffs,
			"overlap": S,
			"Hamiltonian": H,
			"weights": weights,
			"lowdin_weights": lowdin_weights,
			"orbitals": [orbA, orbB],
			"diagnostics": diagnostics,
		}

	def compute_bovb_h2(self,
					 molecule,
					 basis,
					 structures,
					 orbitals,
					 reference_orbitals,
					 options,
					 active_space=None,
					 analysis=None,
					 freeze_inactive=False):
		return self.compute_bovb_two_orbital(
			molecule,
			basis,
			structures,
			orbitals,
			reference_orbitals,
			options,
			active_space=active_space,
			analysis=analysis,
			freeze_inactive=freeze_inactive,
		)

	def compute_bovb_two_orbital(self,
					 molecule,
					 basis,
					 structures,
					 orbitals,
					 reference_orbitals,
					 options,
					 active_space=None,
					 analysis=None,
					 freeze_inactive=False):
		"""Two-orbital/two-electron BOVB with structure-specific active orbitals."""
		from scipy.optimize import minimize

		if len(orbitals) != 2 or len(structures) != 3:
			raise NotImplementedError(
				"Two-orbital BOVB requires two active orbitals and three covalent/ionic structures.")

		H_ao, S_ao, eri_ao, e_nuc = self._ao_integrals(molecule, basis)
		embedding = None
		if freeze_inactive and active_space is not None and analysis is not None:
			embedding = self._frozen_hf_embedding(
				H_ao,
				S_ao,
				eri_ao,
				e_nuc,
				active_space,
				analysis,
			)
			if embedding is not None:
				H_ao = embedding["h_effective"]
				e_nuc = embedding["constant_energy"]
				orbitals = self._project_active_orbitals_out_of_frozen_density(
					orbitals,
					S_ao,
					embedding.get("frozen_density"),
				)
		ao1 = orbitals[0].coefficients
		ao2 = orbitals[1].coefficients
		breath1 = self._center_breathing_vector(
			molecule,
			basis,
			orbitals[0].center,
			ao1,
			(ao1, ao2),
			S_ao,
		)
		breath2 = self._center_breathing_vector(
			molecule,
			basis,
			orbitals[1].center,
			ao2,
			(ao1, ao2),
			S_ao,
		)
		has_external_breathing = breath1 is not None and breath2 is not None
		if not has_external_breathing:
			breath1 = np.zeros_like(ao1)
			breath2 = np.zeros_like(ao2)

		def symmetric_pair(theta, breathing_amplitude):
			left = np.cos(theta) * ao1 + np.sin(theta) * ao2
			right = np.cos(theta) * ao2 + np.sin(theta) * ao1
			left = left + breathing_amplitude * breath1
			right = right + breathing_amplitude * breath2
			left = self._s_normalize(left, S_ao)
			right = self._s_normalize(right, S_ao)
			return np.column_stack([left, right])

		def structure_coefficients(parameters):
			covalent_theta, ionic_theta, covalent_breathing, ionic_breathing = parameters
			return [
				symmetric_pair(covalent_theta, covalent_breathing),
				symmetric_pair(ionic_theta, ionic_breathing),
				symmetric_pair(ionic_theta, ionic_breathing),
			]

		def energy_for_parameters(parameters):
			try:
				S, H = self._build_structure_specific_two_orbital_singlet_matrices(
					structures,
					structure_coefficients(parameters),
					H_ao,
					S_ao,
					eri_ao,
					e_nuc,
				)
				energy, _, _, _ = self._solve_generalized_vb(H, S)
			except Exception:
				energy = 1.0e6
			return energy

		initial_parameters = np.array([0.0, 0.0, 0.0, 0.0])
		initial_energy = energy_for_parameters(initial_parameters)
		conservative_non_h2_limit = not self._is_h2_molecule(molecule)
		if conservative_non_h2_limit:
			class FixedLimitResult:
				success = True
				message = "non-H2 two-orbital BOVB uses fixed-orbital limit pending constrained breathing model"

			result = FixedLimitResult()
			optimizer_energy = initial_energy
		else:
			result = minimize(
				energy_for_parameters,
				x0=initial_parameters,
				method='L-BFGS-B',
				bounds=[(-0.70, 0.70), (-0.70, 0.70), (-1.50, 1.50), (-1.50, 1.50)],
				options={
					'ftol': options.conv_thresh if options else 1.0e-8,
					'gtol': options.conv_thresh if options else 1.0e-8,
					'maxiter': options.max_iter if options else 50,
				},
			)
			optimizer_energy = float(result.fun) if np.isfinite(result.fun) else 1.0e6
		unphysical_non_h2_lowering = (
			not conservative_non_h2_limit and
			not self._is_h2_molecule(molecule) and
			initial_energy - optimizer_energy > 0.50
		)
		used_fixed_limit = (
			conservative_non_h2_limit or
			optimizer_energy > initial_energy + 1.0e-10 or
			unphysical_non_h2_lowering
		)
		parameters = initial_parameters if used_fixed_limit else np.array(result.x, dtype=float)
		structure_orbitals = structure_coefficients(parameters)
		S, H = self._build_structure_specific_two_orbital_singlet_matrices(
			structures,
			structure_orbitals,
			H_ao,
			S_ao,
			eri_ao,
			e_nuc,
		)
		energy, coeffs, weights, kept, lowdin_weights = self._solve_generalized_vb(
			H, S, return_lowdin=True)
		model_prefix = "h2" if self._is_h2_molecule(molecule) else "two-orbital"
		message_prefix = "H2" if self._is_h2_molecule(molecule) else "Two-orbital"
		diagnostics = {
			"message": f"{message_prefix} BOVB result with structure-specific breathing orbitals",
			"bovb_model": f"{model_prefix}-center-local-structure-specific-breathing-orbitals",
			"bovb_covalent_theta": float(parameters[0]),
			"bovb_ionic_theta": float(parameters[1]),
			"bovb_covalent_breathing": float(parameters[2]),
			"bovb_ionic_breathing": float(parameters[3]),
			"bovb_has_external_breathing_space": bool(has_external_breathing),
			"bovb_optimizer_success": bool(result.success),
			"bovb_optimizer_message": str(result.message),
			"bovb_initial_energy": float(initial_energy),
			"bovb_optimizer_energy": float(optimizer_energy),
			"bovb_used_fixed_orbital_limit": bool(used_fixed_limit),
			"bovb_conservative_non_h2_limit": bool(conservative_non_h2_limit),
			"bovb_rejected_unphysical_non_h2_lowering": bool(unphysical_non_h2_lowering),
			"overlap_condition": float(np.linalg.cond(S)),
			"overlap_eigenvalues": np.linalg.eigvalsh(S).tolist(),
			"retained_overlap_rank": int(len(kept)),
			"weight_scheme": "Chirgwin-Coulson",
			"available_weight_schemes": ["Chirgwin-Coulson", "Lowdin"],
			"structure_specific_orbitals": True,
		}
		if embedding is not None:
			diagnostics.update({
				"frozen_hf_embedding": True,
				"active_orbitals_orthogonalized_to_frozen_space": True,
				"frozen_electron_count": embedding["frozen_electron_count"],
				"active_reference_electron_count": embedding["active_electron_count"],
				"frozen_constant_energy": embedding["constant_energy"],
				"embedding_model": "inactive HF density = total RHF density minus active reference density",
			})
		return {
			"energy": float(energy),
			"structure_coefficients": coeffs,
			"overlap": S,
			"Hamiltonian": H,
			"weights": weights,
			"lowdin_weights": lowdin_weights,
			"orbitals": [matrix for matrix in structure_orbitals],
			"diagnostics": diagnostics,
		}

	def _is_minimal_h2_vbci_case(self, molecule, structures, orbitals):
		try:
			return (
				hasattr(molecule, 'number_of_atoms') and molecule.number_of_atoms() == 2
				and len(orbitals) == 2
				and len(structures) == 3
			)
		except Exception:
			return False

	def _is_h2_three_structure_case(self, molecule, structures):
		try:
			return (
				hasattr(molecule, 'number_of_atoms') and molecule.number_of_atoms() == 2
				and tuple(molecule.get_labels()) == ('H', 'H')
				and structures is not None
				and len(structures) == 3
			)
		except Exception:
			return False

	def _is_h2_molecule(self, molecule):
		try:
			return (
				hasattr(molecule, 'number_of_atoms') and
				molecule.number_of_atoms() == 2 and
				tuple(molecule.get_labels()) == ('H', 'H')
			)
		except Exception:
			return False

	def _candidate_label(self, candidate):
		return f"{candidate['type']}_{candidate.get('subtype', '')}_{candidate['index']}".replace("__", "_")

	def _candidate_coefficient_vector(self, candidate, n_ao):
		coefficients = candidate.get("coefficients", {})
		vector = np.zeros(n_ao)
		if isinstance(coefficients, dict):
			for index, value in coefficients.items():
				index = int(index)
				if 0 <= index < n_ao:
					vector[index] = float(value)
		else:
			for item in coefficients:
				if not isinstance(item, dict):
					continue
				index = int(item.get("nao_index", 0)) - 1
				if 0 <= index < n_ao:
					vector[index] = float(item.get("coefficient", 0.0))
		return vector

	def _candidate_summary_record(self, candidate):
		record = {
			"label": self._candidate_label(candidate),
			"type": candidate.get("type"),
			"subtype": candidate.get("subtype"),
			"atoms": tuple(candidate.get("atoms", ())),
			"occupation": candidate.get("occupation"),
			"source": candidate.get("source"),
		}
		for key in (
			"channel",
			"metal_atom",
			"ligand_atom",
			"donor_atom",
			"acceptor_atom",
			"donor_acceptor_role",
			"coordination_mode",
			"interaction_strength",
			"donation_strength",
			"back_donation_strength",
		):
			if key in candidate:
				record[key] = candidate.get(key)
		return record

	def _bond_candidate_summary(self, candidates):
		bond_records = []
		for candidate in candidates:
			if candidate.get("type") != "BD":
				continue
			if candidate.get("subtype") not in ("sigma", "pi"):
				continue
			atoms = tuple(candidate.get("atoms", ()))
			if len(atoms) != 2:
				continue
			bond_records.append(self._candidate_summary_record(candidate))
		return sorted(
			bond_records,
			key=lambda item: (
				tuple(sorted(item.get("atoms", ()))),
				item.get("subtype") or "",
				-(float(item.get("occupation") or 0.0)),
			),
		)

	def _candidate_partition_diagnostics(self, candidates):
		partitions = {
			"core": [],
			"sigma_bonds": [],
			"pi_bonds": [],
			"lone_pairs": [],
			"radicals": [],
			"metal_ligand": [],
			"antibonds": [],
			"rydberg": [],
			"other": [],
		}
		for candidate in candidates:
			record = self._candidate_summary_record(candidate)
			candidate_type = candidate.get("type")
			candidate_subtype = candidate.get("subtype")
			if candidate_type == "CR":
				partitions["core"].append(record)
			elif candidate_type == "BD" and candidate_subtype == "sigma":
				partitions["sigma_bonds"].append(record)
			elif candidate_type == "BD" and candidate_subtype == "pi":
				partitions["pi_bonds"].append(record)
			elif candidate_type == "LP":
				partitions["lone_pairs"].append(record)
			elif candidate_type == "SOMO":
				partitions["radicals"].append(record)
			elif candidate_type == "ML":
				partitions["metal_ligand"].append(record)
			elif candidate_type == "BD*" or candidate_subtype == "antibonding":
				partitions["antibonds"].append(record)
			elif candidate_type == "RY":
				partitions["rydberg"].append(record)
			else:
				partitions["other"].append(record)
		return partitions

	def _suggest_active_candidates(self, candidates):
		partitions = self._candidate_partition_diagnostics(candidates)
		suggestions = []
		for key in ("pi_bonds", "sigma_bonds"):
			if partitions[key]:
				suggestions.append(
					max(partitions[key],
						key=lambda item: float(item.get("occupation") or 0.0)))
		return suggestions

	def _should_use_one_active_bond_space(self,
									  molecule,
									  options,
									  user_supplied_structures,
									  user_supplied_orbitals):
		if not getattr(options, 'use_active_space', True):
			return False
		if user_supplied_structures or user_supplied_orbitals:
			return False
		if getattr(options, 'active_bond', None) is not None:
			return True
		if getattr(options, 'active_candidate_label', None) is not None:
			return True
		if getattr(options, 'active_candidate_subtype', None) is not None:
			return True
		if getattr(options, 'active_metal_ligand_channels', None) is not None:
			return True
		if getattr(options, 'active_pi_atoms', None) is not None:
			return True
		if bool(getattr(options, 'include_bovb', False)) or getattr(options, 'mode', None) == 'bovb':
			return True
		return self._is_h2_molecule(molecule)

	def _build_one_active_bond_space(self, molecule, basis, candidates, nao_data, options):
		if getattr(options, 'active_metal_ligand_channels', None) is not None:
			return self._build_metal_ligand_active_space(
				molecule,
				basis,
				candidates,
				nao_data,
				options,
			)

		if getattr(options, 'active_pi_atoms', None) is not None:
			active_pi_atoms = tuple(int(atom) for atom in options.active_pi_atoms)
			if len(active_pi_atoms) > 2:
				return self._build_multi_center_pi_space(
					molecule,
					basis,
					candidates,
					nao_data,
					options,
					active_pi_atoms,
				)

		explicit_active_selection = any(
			getattr(options, attribute, None) is not None
			for attribute in (
				'active_bond',
				'active_candidate_label',
				'active_candidate_subtype',
			)
		)
		stable_h2_active_space = self._is_h2_molecule(molecule) and not explicit_active_selection
		active_candidate = self._select_active_bond_candidate(candidates, options)
		if stable_h2_active_space:
			selected_index = int(active_candidate.get('index', 0)) if active_candidate is not None else 0
			active_candidate = {
				'index': selected_index,
				'type': 'H2',
				'subtype': 'atom_centered_sigma',
				'atoms': (0, 1),
				'occupation': 2.0,
				'coefficients': [],
				'source': 'H2 stable atom-centered active bond',
			}
		if active_candidate is None:
			requested_bond = getattr(options, 'active_bond', None)
			if requested_bond is not None:
				active_candidate = {
					'index': -1,
					'type': 'BD',
					'subtype': getattr(options, 'active_candidate_subtype', None) or 'sigma',
					'atoms': tuple(int(atom) for atom in requested_bond),
					'occupation': 0.0,
					'coefficients': [],
					'source': 'explicit active_bond fallback for dissociation scan',
				}
			elif not self._is_h2_molecule(molecule):
				raise RuntimeError("VB active-space builder could not find a two-center bond candidate.")
			else:
				active_candidate = {
					'index': 0,
					'type': 'BD',
					'subtype': 'sigma',
					'atoms': (0, 1),
					'occupation': 2.0,
					'coefficients': [],
					'source': 'H2 fallback active bond',
				}

		active_bond = tuple(int(atom) for atom in active_candidate.get('atoms', ()))
		if len(active_bond) != 2:
			raise RuntimeError("VB active-space builder requires a two-center active bond.")

		reference_active_orbitals = self._reference_active_bond_orbitals(
			getattr(options, 'active_bond_reference_orbitals', None),
			active_bond,
		)
		active_index = int(active_candidate.get('index', -1))
		active_subtype = 'sigma' if stable_h2_active_space else active_candidate.get('subtype')
		stretched_active_orbitals = None
		stretched_active_orbital_source = None
		if (active_index < 0 and active_subtype == 'sigma' and
				self._active_bond_distance_angstrom(molecule, active_bond) >= 1.8):
			stretched_active_orbitals = self._uhf_frontier_active_orbitals(
				molecule,
				basis,
				active_bond,
			)
			if stretched_active_orbitals is not None:
				stretched_active_orbital_source = 'uhf_frontier'
			else:
				stretched_active_orbitals = self._fragment_radical_active_orbitals(
					molecule,
					basis,
					active_bond,
				)
				if stretched_active_orbitals is not None:
					stretched_active_orbital_source = 'fragment_radical'
		if stretched_active_orbitals is not None:
			active_orbitals = stretched_active_orbitals
		elif reference_active_orbitals is not None:
			active_orbitals = reference_active_orbitals
		elif stable_h2_active_space:
			active_orbitals = tuple(
				self._h2_fragment_active_orbitals(molecule, basis, active_bond) or
				self._atom_centered_active_orbitals(molecule, basis, active_bond))
		else:
			active_orbitals = tuple(
				self._candidate_centered_active_orbitals(
					molecule,
					basis,
					active_candidate,
					nao_data,
				) or self._atom_centered_active_orbitals(molecule, basis, active_bond))
		structures = tuple(self._default_one_bond_structures(active_bond, options))
		if stable_h2_active_space:
			active_label = "H2_atom_centered_sigma"
		elif active_index < 0:
			active_label = f"explicit_{active_subtype}_bond_{active_bond[0] + 1}_{active_bond[1] + 1}"
		else:
			active_label = self._candidate_label(active_candidate)

		frozen = []
		inactive = []
		excluded = []
		for candidate in candidates:
			if int(candidate.get('index', -1)) == active_index:
				continue
			candidate_type = candidate.get('type')
			candidate_subtype = candidate.get('subtype')
			if candidate_type == 'CR':
				frozen.append(candidate)
			elif candidate_type in ('BD', 'LP', 'SOMO') and candidate_subtype != 'antibonding':
				inactive.append(candidate)
			else:
				excluded.append(candidate)

		return VbActiveSpace(
			active_bond=active_bond,
			active_candidate_label=active_label,
			active_candidate=active_candidate,
			active_orbitals=active_orbitals,
			structures=structures,
			frozen_candidates=tuple(frozen),
			inactive_candidates=tuple(inactive),
			excluded_candidates=tuple(excluded),
			metadata={
				"source": "OrbitalAnalyzer",
				"model": "one-active-pi-bond" if active_subtype == 'pi' else "one-active-bond",
				"electron_count": 2,
				"include_ionic": bool(options.include_ionic),
				"frozen_hf_reference": bool(options.freeze_inactive_orbitals),
				"h2_stable_atom_centered_active_space": bool(stable_h2_active_space),
				"active_bond_state_tracked": reference_active_orbitals is not None,
				"active_candidate_fallback": active_index < 0,
				"active_uhf_frontier_orbitals": stretched_active_orbital_source == 'uhf_frontier',
				"active_fragment_radical_orbitals": stretched_active_orbital_source == 'fragment_radical',
				"active_stretched_orbital_source": stretched_active_orbital_source,
			},
		)

	def _reference_active_bond_orbitals(self, reference_orbitals, active_bond):
		if reference_orbitals is None:
			return None
		if len(reference_orbitals) != 2:
			raise RuntimeError("State-tracked organic bond scans require exactly two reference active orbitals.")
		tracked_orbitals = []
		for index, orbital in enumerate(reference_orbitals):
			coefficients = np.array(getattr(orbital, 'coefficients', orbital), dtype=float)
			center = getattr(orbital, 'center', None)
			if center is not None and int(center) not in set(active_bond):
				raise RuntimeError("Reference active orbital center is not part of the requested active bond.")
			label = getattr(orbital, 'label', f'state_tracked_active_{index + 1}')
			tracked_orbitals.append(
				VbOrbital(
					label=f"tracked_{label}",
					coefficients=coefficients,
					center=None if center is None else int(center),
					kind=getattr(orbital, 'kind', 'active'),
				)
			)
		return tuple(tracked_orbitals)

	def _build_multi_center_pi_space(self,
								   molecule,
								   basis,
								   candidates,
								   nao_data,
								   options,
								   active_pi_atoms):
		if len(active_pi_atoms) < 3:
			raise RuntimeError("Multi-center pi active spaces require at least three active atoms.")
		spin = str(getattr(options, 'active_spin', 'singlet') or 'singlet').lower()
		electron_count = int(getattr(options, 'active_electron_count', None) or 2)
		n_active_orbitals = len(active_pi_atoms)
		n_alpha, n_beta = self._active_spin_occupations(electron_count, spin)
		if electron_count < 1 or electron_count > 2 * n_active_orbitals:
			raise RuntimeError("Active electron count is incompatible with the requested pi active orbitals.")

		active_atom_set = set(active_pi_atoms)
		pi_candidates = [
			candidate for candidate in candidates
			if candidate.get('type') == 'BD' and
			candidate.get('subtype') == 'pi' and
			set(int(atom) for atom in candidate.get('atoms', ())).issubset(active_atom_set) and
			len(candidate.get('atoms', ())) == 2
		]
		if not pi_candidates:
			raise RuntimeError(
				"VB pi active-space builder could not find analyzer pi candidates on the requested atoms.")

		active_candidate = max(
			pi_candidates,
			key=lambda candidate: float(candidate.get('occupation', 0.0)),
		)
		active_orbitals = tuple(
			self._multi_center_pi_active_orbitals(
				molecule,
				basis,
				active_pi_atoms,
				pi_candidates,
				nao_data,
			)
		)
		determinant_ci = not (spin == 'singlet' and electron_count == 2)
		if determinant_ci:
			structures = tuple(
				self._determinant_pi_structures(
					active_pi_atoms,
					n_alpha,
					n_beta,
					spin,
				)
			)
		else:
			structures = tuple(
				self._default_two_electron_pi_structures(
					active_pi_atoms,
					include_ionic=bool(options.include_ionic),
				)
			)
		active_indices = {int(candidate.get('index', -1)) for candidate in pi_candidates}
		frozen = []
		inactive = []
		excluded = []
		for candidate in candidates:
			if int(candidate.get('index', -1)) in active_indices:
				continue
			candidate_type = candidate.get('type')
			candidate_subtype = candidate.get('subtype')
			if candidate_type == 'CR':
				frozen.append(candidate)
			elif candidate_type in ('BD', 'LP', 'SOMO') and candidate_subtype != 'antibonding':
				inactive.append(candidate)
			else:
				excluded.append(candidate)

		return VbActiveSpace(
			active_bond=tuple(active_pi_atoms),
			active_candidate_label="pi_system_" + "_".join(str(atom + 1) for atom in active_pi_atoms),
			active_candidate={
				'index': int(active_candidate.get('index', -1)),
				'type': 'PI_SYSTEM',
				'subtype': 'pi',
				'atoms': tuple(active_pi_atoms),
				'occupation': float(sum(float(c.get('occupation', 0.0)) for c in pi_candidates)),
				'source': 'multi-center pi system assembled from OrbitalAnalyzer pi candidates',
				'component_candidate_labels': [self._candidate_label(c) for c in pi_candidates],
			},
			active_orbitals=active_orbitals,
			structures=structures,
			frozen_candidates=tuple(frozen),
			inactive_candidates=tuple(inactive),
			excluded_candidates=tuple(excluded),
			metadata={
				"source": "OrbitalAnalyzer",
				"model": "fixed-orbital-multicenter-pi",
				"electron_count": electron_count,
				"spin": spin,
				"n_alpha": n_alpha,
				"n_beta": n_beta,
				"determinant_ci": determinant_ci,
				"determinant_count": len(structures),
				"active_pi_atoms": tuple(active_pi_atoms),
				"component_pi_candidate_labels": [self._candidate_label(c) for c in pi_candidates],
				"include_ionic": bool(options.include_ionic),
				"frozen_hf_reference": bool(options.freeze_inactive_orbitals),
			},
		)

	def _build_metal_ligand_active_space(self,
									molecule,
									basis,
									candidates,
									nao_data,
									options):
		requested_channels = tuple(
			str(channel).lower()
			for channel in getattr(options, 'active_metal_ligand_channels', ())
		)
		if not requested_channels:
			raise RuntimeError("No metal-ligand active channels were requested.")

		selected = []
		for channel in requested_channels:
			matches = [
				candidate for candidate in candidates
				if candidate.get('type') == 'ML' and
				str(candidate.get('subtype', '')).lower() == channel
			]
			if not matches:
				raise RuntimeError(
					f"No metal-ligand candidate found for channel '{channel}'.")
			selected.append(
				max(matches,
					key=lambda candidate: float(
						candidate.get('interaction_strength', 0.0))))

		backdonation_enabled = 'pi-donor' in requested_channels
		if requested_channels == ('sigma-acceptor',):
			metal_ligand_model = 'sigma-only'
		elif requested_channels == ('sigma-acceptor', 'pi-donor'):
			metal_ligand_model = 'sigma-plus-backdonation'
		else:
			metal_ligand_model = 'custom-metal-ligand-channels'
		selected_records = [
			{
				'label': self._candidate_label(candidate),
				'subtype': candidate.get('subtype'),
				'channel': candidate.get('channel'),
				'metal_atom': candidate.get('metal_atom'),
				'ligand_atom': candidate.get('ligand_atom'),
				'donor_atom': candidate.get('donor_atom'),
				'acceptor_atom': candidate.get('acceptor_atom'),
				'donor_acceptor_role': candidate.get('donor_acceptor_role'),
				'donation_strength': float(
					candidate.get('donation_strength', 0.0)),
				'back_donation_strength': float(
					candidate.get('back_donation_strength', 0.0)),
				'interaction_strength': float(
					candidate.get('interaction_strength', 0.0)),
				'occupation': float(candidate.get('occupation', 0.0)),
			}
			for candidate in selected
		]

		reference_orbitals = getattr(
			options,
			'active_metal_ligand_reference_orbitals',
			None,
		)
		if reference_orbitals:
			active_orbitals = self._tracked_metal_ligand_active_orbitals(
				reference_orbitals,
				molecule,
				basis,
			)
			orbital_model = "state-tracked-reference-orbitals"
		else:
			active_orbitals = []
			for candidate in selected:
				active_orbitals.extend(
					self._metal_ligand_channel_active_orbitals(
						candidate,
						nao_data,
						molecule,
						basis,
					)
				)
			orbital_model = "candidate-side-projected-orbitals"
		n_orbitals = len(active_orbitals)
		electron_count = int(getattr(options, 'active_electron_count', None) or
						 2 * len(selected))
		spin = str(getattr(options, 'active_spin', 'singlet') or 'singlet').lower()
		n_alpha, n_beta = self._active_spin_occupations(electron_count, spin)
		structures = tuple(
			self._determinant_active_structures(n_orbitals, n_alpha, n_beta))
		selected_indices = {int(candidate.get('index', -1)) for candidate in selected}
		active_atoms = tuple(
			int(atom) for atom in sorted({
				atom for candidate in selected
				for atom in candidate.get('atoms', ())
			}))
		active_centers = tuple(
			int(orbital.center) if orbital.center is not None else index
			for index, orbital in enumerate(active_orbitals)
		)

		frozen = []
		inactive = []
		excluded = []
		for candidate in candidates:
			if int(candidate.get('index', -1)) in selected_indices:
				continue
			candidate_type = candidate.get('type')
			candidate_subtype = candidate.get('subtype')
			if candidate_type == 'CR':
				frozen.append(candidate)
			elif candidate_type in ('BD', 'LP', 'SOMO') and candidate_subtype != 'antibonding':
				inactive.append(candidate)
			else:
				excluded.append(candidate)

		channel_label = "+".join(requested_channels)
		return VbActiveSpace(
			active_bond=active_atoms,
			active_candidate_label="metal_ligand_" + channel_label.replace('-', '_'),
			active_candidate={
				'index': -1,
				'type': 'ML_SYSTEM',
				'subtype': channel_label,
				'atoms': active_atoms,
				'occupation': float(sum(
					float(candidate.get('occupation', 0.0))
					for candidate in selected)),
				'source': 'combined metal-ligand active space',
				'component_candidate_labels': [
					self._candidate_label(candidate) for candidate in selected
				],
			},
			active_orbitals=tuple(active_orbitals),
			structures=structures,
			frozen_candidates=tuple(frozen),
			inactive_candidates=tuple(inactive),
			excluded_candidates=tuple(excluded),
			metadata={
				"source": "OrbitalAnalyzer",
				"model": "metal-ligand-channel-determinant-ci",
				"metal_ligand_active_orbital_model": orbital_model,
				"metal_ligand_state_tracked_orbitals": bool(reference_orbitals),
				"metal_ligand_model": metal_ligand_model,
				"electron_count": electron_count,
				"spin": spin,
				"n_alpha": n_alpha,
				"n_beta": n_beta,
				"determinant_ci": True,
				"determinant_count": len(structures),
				"active_pi_atoms": active_centers,
				"metal_ligand_channels": requested_channels,
				"backdonation_enabled": bool(backdonation_enabled),
				"backdonation_blocked": not bool(backdonation_enabled),
				"selected_metal_ligand_records": selected_records,
				"component_candidate_labels": [
					self._candidate_label(candidate) for candidate in selected
				],
				"include_ionic": bool(options.include_ionic),
				"frozen_hf_reference": bool(options.freeze_inactive_orbitals),
			},
		)

	def _metal_ligand_channel_active_orbitals(self,
									candidate,
									nao_data,
									molecule,
									basis):
		if nao_data is None:
			raise RuntimeError(
				"Metal-ligand active-space construction requires NAO data.")
		candidate_vector = self._candidate_nao_vector(candidate, nao_data)
		if candidate_vector is None:
			raise RuntimeError(
				"Metal-ligand candidate has no usable NAO coefficient vector.")
		import veloxchem as vlx
		S_ao = vlx.OverlapDriver().compute(molecule, basis).to_numpy()
		atom_map = np.array(nao_data.atom_map, dtype=int)
		orbitals = []
		subtype = str(candidate.get('subtype') or 'ml').replace('-', '_')
		for role, atom_key in (('donor', 'donor_atom'), ('acceptor', 'acceptor_atom')):
			atom = int(candidate.get(atom_key))
			indices = np.where(atom_map == atom)[0]
			if len(indices) == 0:
				raise RuntimeError(
					f"No NAO functions found on metal-ligand {role} atom {atom}.")
			nao_side = np.zeros_like(candidate_vector)
			nao_side[indices] = candidate_vector[indices]
			if np.linalg.norm(nao_side) <= 1.0e-14:
				raise RuntimeError(
					f"No {role} contribution in metal-ligand candidate {subtype}.")
			coeffs = nao_data.transform @ nao_side
			coeffs = self._s_normalize(coeffs, S_ao)
			orbitals.append(
				VbOrbital(
					label=f"active_{subtype}_{role}_atom_{atom + 1}",
					coefficients=coeffs,
					center=atom,
					kind="active",
				)
			)
		return orbitals

	def _tracked_metal_ligand_active_orbitals(self,
									 reference_orbitals,
									 molecule,
									 basis):
		import veloxchem as vlx
		S_ao = vlx.OverlapDriver().compute(molecule, basis).to_numpy()
		orbitals = []
		for index, orbital in enumerate(reference_orbitals):
			if isinstance(orbital, VbOrbital):
				label = orbital.label
				center = orbital.center
				kind = orbital.kind
				coefficients = orbital.coefficients
			else:
				label = getattr(orbital, 'label', f'tracked_ml_orbital_{index + 1}')
				center = getattr(orbital, 'center', None)
				kind = getattr(orbital, 'kind', 'active')
				coefficients = getattr(orbital, 'coefficients', orbital)
			coefficients = np.array(coefficients, dtype=float, copy=True)
			if coefficients.ndim != 1 or coefficients.shape[0] != S_ao.shape[0]:
				raise RuntimeError(
					"Tracked metal-ligand reference orbitals must match the AO basis dimension.")
			coefficients = self._s_normalize(coefficients, S_ao)
			orbitals.append(
				VbOrbital(
					label=str(label),
					coefficients=coefficients,
					center=center,
					kind=str(kind),
				)
			)
		if not orbitals:
			raise RuntimeError("No tracked metal-ligand reference orbitals were provided.")
		return orbitals

	def _spin_adapted_metal_ligand_sigma_space(self, active_space, options):
		records = active_space.metadata.get('selected_metal_ligand_records', ())
		if len(records) != 1 or records[0].get('subtype') != 'sigma-acceptor':
			raise RuntimeError(
				"Metal-ligand VB-SCF is currently available only for the sigma-only channel.")
		if len(active_space.active_orbitals) != 2:
			raise RuntimeError(
				"Metal-ligand sigma VB-SCF requires exactly two active orbitals.")
		metadata = dict(active_space.metadata)
		metadata.update({
			"model": "metal-ligand-sigma-spin-adapted-vbscf",
			"determinant_ci": False,
			"determinant_count": 3 if bool(options.include_ionic) else 1,
			"spin_adapted_sigma_vbscf": True,
		})
		return VbActiveSpace(
			active_bond=active_space.active_bond,
			active_candidate_label=active_space.active_candidate_label,
			active_candidate=active_space.active_candidate,
			active_orbitals=active_space.active_orbitals,
			structures=tuple(
				self._default_one_bond_structures(
					active_space.active_bond,
					options,
				)
			),
			frozen_candidates=active_space.frozen_candidates,
			inactive_candidates=active_space.inactive_candidates,
			excluded_candidates=active_space.excluded_candidates,
			metadata=metadata,
		)

	def _determinant_active_structures(self, n_orbitals, n_alpha, n_beta):
		from itertools import combinations

		structures = []
		for alpha_sites in combinations(range(n_orbitals), n_alpha):
			alpha_occ = [0] * n_orbitals
			for site in alpha_sites:
				alpha_occ[site] = 1
			for beta_sites in combinations(range(n_orbitals), n_beta):
				beta_occ = [0] * n_orbitals
				for site in beta_sites:
					beta_occ[site] = 1
				alpha_label = "a" + "".join(str(site + 1) for site in alpha_sites)
				beta_label = "b" + "".join(str(site + 1) for site in beta_sites)
				structures.append(
					VbStructure(
						label=f"det_{alpha_label}_{beta_label}",
						occupation=(tuple(alpha_occ), tuple(beta_occ)),
						spin='determinant',
					)
				)
		return structures

	def _active_spin_occupations(self, electron_count, spin):
		if spin == 'singlet':
			if electron_count % 2 != 0:
				raise RuntimeError("Singlet pi active spaces require an even active electron count.")
			return electron_count // 2, electron_count // 2
		if spin == 'doublet':
			if electron_count % 2 != 1:
				raise RuntimeError("Doublet pi active spaces require an odd active electron count.")
			return (electron_count + 1) // 2, (electron_count - 1) // 2
		raise RuntimeError("Generated pi active spaces currently support singlet and doublet spin sectors.")

	def _multi_center_pi_active_orbitals(self,
								   molecule,
								   basis,
								   active_pi_atoms,
								   pi_candidates,
								   nao_data):
		if nao_data is None:
			raise RuntimeError("Multi-center pi active-space construction requires NAO data.")
		import veloxchem as vlx
		S_ao = vlx.OverlapDriver().compute(molecule, basis).to_numpy()
		atom_map = np.array(nao_data.atom_map, dtype=int)
		orbitals = []
		candidate_vectors = [
			self._candidate_nao_vector(candidate, nao_data)
			for candidate in pi_candidates
		]
		candidate_vectors = [vector for vector in candidate_vectors if vector is not None]
		for ordinal, atom in enumerate(active_pi_atoms, start=1):
			indices = np.where(atom_map == atom)[0]
			if len(indices) == 0:
				raise RuntimeError(f"No NAO functions found for active pi atom {atom}.")
			nao_side = np.zeros(len(nao_data.populations))
			for candidate_vector in candidate_vectors:
				local = np.zeros_like(candidate_vector)
				local[indices] = candidate_vector[indices]
				if np.linalg.norm(local) <= 1.0e-14:
					continue
				if np.dot(nao_side, local) < 0.0:
					local = -local
				nao_side += local
			if np.linalg.norm(nao_side) <= 1.0e-14:
				raise RuntimeError(
					f"No pi contribution from analyzer candidates on active atom {atom}.")
			coeffs = nao_data.transform @ (nao_side / np.linalg.norm(nao_side))
			coeffs = self._s_normalize(coeffs, S_ao)
			orbitals.append(
				VbOrbital(
					label=f"active_pi_center_{ordinal}_atom_{atom + 1}",
					coefficients=coeffs,
					center=int(atom),
					kind="active",
				)
			)
		return orbitals

	def _default_two_electron_pi_structures(self, active_pi_atoms, include_ionic=True):
		structures = []
		n_active = len(active_pi_atoms)
		for i in range(n_active):
			for j in range(i + 1, n_active):
				occupation = [0] * n_active
				occupation[i] = 1
				occupation[j] = 1
				structures.append(
					VbStructure(
						label=f"covalent_atom_{active_pi_atoms[i] + 1}_atom_{active_pi_atoms[j] + 1}",
						occupation=(tuple(occupation),),
						spin='singlet',
						charge_pattern={active_pi_atoms[i]: 0, active_pi_atoms[j]: 0},
					)
				)
		if include_ionic:
			for i, atom in enumerate(active_pi_atoms):
				occupation = [0] * n_active
				occupation[i] = 2
				structures.append(
					VbStructure(
						label=f"ionic_atom_{atom + 1}",
						occupation=(tuple(occupation),),
						spin='singlet',
						charge_pattern={atom: -1},
					)
				)
		return structures

	def _determinant_pi_structures(self, active_pi_atoms, n_alpha, n_beta, spin):
		from itertools import combinations

		structures = []
		n_active = len(active_pi_atoms)
		for alpha_sites in combinations(range(n_active), n_alpha):
			alpha_occ = [0] * n_active
			for site in alpha_sites:
				alpha_occ[site] = 1
			for beta_sites in combinations(range(n_active), n_beta):
				beta_occ = [0] * n_active
				for site in beta_sites:
					beta_occ[site] = 1
				alpha_label = "a" + "".join(str(active_pi_atoms[i] + 1) for i in alpha_sites)
				beta_label = "b" + "".join(str(active_pi_atoms[i] + 1) for i in beta_sites)
				charge_pattern = {}
				for i, atom in enumerate(active_pi_atoms):
					occupation = alpha_occ[i] + beta_occ[i]
					if occupation == 2:
						charge_pattern[atom] = -1
					elif occupation == 0:
						charge_pattern[atom] = +1
				spatial_label = self._spatial_occupation_label(
					active_pi_atoms,
					tuple(alpha_occ[i] + beta_occ[i] for i in range(n_active)),
					spin,
				)
				structures.append(
					VbStructure(
						label=f"{spatial_label}__{alpha_label}_{beta_label}",
						occupation=(tuple(alpha_occ), tuple(beta_occ)),
						spin=spin,
						charge_pattern=charge_pattern,
					)
				)
		return structures

	def _spatial_occupation_label(self, active_pi_atoms, spatial_occupation, spin):
		fragments = []
		empty_atoms = []
		singly_atoms = []
		doubly_atoms = []
		for atom, occupation in zip(active_pi_atoms, spatial_occupation):
			atom_label = f"atom_{atom + 1}"
			if occupation == 0:
				empty_atoms.append(atom_label)
			elif occupation == 1:
				singly_atoms.append(atom_label)
			elif occupation == 2:
				doubly_atoms.append(atom_label)
			else:
				fragments.append(f"occ{occupation}_{atom_label}")
		if spin == 'doublet' and len(singly_atoms) == 1:
			fragments.append("radical_" + singly_atoms[0])
		elif singly_atoms:
			fragments.append("singlet_open_shell_" + "_".join(singly_atoms))
		if doubly_atoms:
			fragments.append("lone_pair_" + "_".join(doubly_atoms))
		if empty_atoms:
			fragments.append("empty_" + "_".join(empty_atoms))
		if not fragments:
			fragments.append("closed_shell")
		return "_".join(fragments)

	def _select_active_bond_candidate(self, candidates, options):
		requested_label = getattr(options, 'active_candidate_label', None)
		requested_bond = getattr(options, 'active_bond', None)
		requested_subtype = getattr(options, 'active_candidate_subtype', None)
		if requested_subtype is not None:
			requested_subtype = requested_subtype.lower()
		if requested_bond is not None:
			requested_bond = tuple(sorted(int(atom) for atom in requested_bond))

		bond_candidates = []
		for candidate in candidates:
			candidate_type = candidate.get('type')
			if candidate_type not in ('BD', 'ML'):
				continue
			if (candidate_type == 'BD' and
					candidate.get('subtype') not in ('sigma', 'pi')):
				continue
			if candidate_type == 'ML' and requested_label is None and requested_subtype is None:
				continue
			atoms = tuple(int(atom) for atom in candidate.get('atoms', ()))
			if len(atoms) != 2:
				continue
			label = self._candidate_label(candidate)
			if requested_label is not None and label == requested_label:
				return candidate
			if requested_subtype is not None and candidate.get('subtype') != requested_subtype:
				continue
			if requested_bond is not None and tuple(sorted(atoms)) == requested_bond:
				bond_candidates.append(candidate)
			elif requested_bond is None and requested_label is None:
				bond_candidates.append(candidate)

		if not bond_candidates:
			return None
		return max(bond_candidates,
				   key=lambda candidate: float(candidate.get('occupation', 0.0)))

	def _atom_centered_active_orbitals(self, molecule, basis, active_bond):
		ao_to_atom, _ = OrbitalAnalyzer.ao_shell_map(molecule, basis)
		import veloxchem as vlx
		S_ao = vlx.OverlapDriver().compute(molecule, basis).to_numpy()
		orbitals = []
		for side, atom in zip(('A', 'B'), active_bond):
			coeffs = np.zeros(basis.get_dimensions_of_basis())
			indices = np.where(ao_to_atom == atom)[0]
			if len(indices) == 0:
				raise RuntimeError(f"No AO functions found for active-bond atom {atom}.")
			coeffs[indices] = 1.0
			try:
				coeffs = self._s_normalize(coeffs, S_ao)
			except RuntimeError:
				raise RuntimeError(f"Zero active orbital generated for atom {atom}.")
			orbitals.append(
				VbOrbital(
					label=f"active_{side}_atom_{atom + 1}",
					coefficients=coeffs,
					center=int(atom),
					kind="active",
				)
			)
		return orbitals

	def _h2_fragment_active_orbitals(self, molecule, basis, active_bond):
		if not self._is_h2_molecule(molecule):
			return None
		try:
			from .molecule import Molecule
			from .molecularbasis import MolecularBasis
			from .scfunrestdriver import ScfUnrestrictedDriver
			import veloxchem as vlx

			basis_label = basis.get_label()
			if not basis_label:
				basis_label = basis.get_main_basis_label()
			h_atom = Molecule.read_str("H 0.0 0.0 0.0")
			h_atom.set_charge(0)
			h_atom.set_multiplicity(2)
			h_basis = MolecularBasis.read(h_atom, basis_label, ostream=None)
			scf_driver = ScfUnrestrictedDriver()
			scf_driver.ostream.mute()
			scf_driver.xcfun = "hf"
			scf_driver.guess_unpaired_electrons = "1(1)"
			scf_driver.compute(h_atom, h_basis)
			fragment_coefficients = np.array(
				scf_driver.molecular_orbitals._orbitals[0][:, 0],
				dtype=float,
			)

			ao_to_atom, _ = OrbitalAnalyzer.ao_shell_map(molecule, basis)
			S_ao = vlx.OverlapDriver().compute(molecule, basis).to_numpy()
			orbitals = []
			for side, atom in zip(('A', 'B'), active_bond):
				indices = np.where(ao_to_atom == atom)[0]
				if len(indices) != len(fragment_coefficients):
					return None
				coeffs = np.zeros(basis.get_dimensions_of_basis())
				coeffs[indices] = fragment_coefficients
				coeffs = self._s_normalize(coeffs, S_ao)
				orbitals.append(
					VbOrbital(
						label=f"active_{side}_atom_{atom + 1}",
						coefficients=coeffs,
						center=int(atom),
						kind="active",
					)
				)
			return orbitals
		except Exception:
			return None

	def _uhf_frontier_active_orbitals(self, molecule, basis, active_bond):
		try:
			from .scfunrestdriver import ScfUnrestrictedDriver
			import veloxchem as vlx

			scf_driver = ScfUnrestrictedDriver()
			scf_driver.ostream.mute()
			scf_driver.xcfun = "hf"
			scf_driver.guess_unpaired_electrons = "1(1),2(-1)"
			scf_driver.compute(molecule, basis)
			n_electrons = int(molecule.number_of_electrons())
			try:
				multiplicity = int(molecule.get_multiplicity())
			except Exception:
				multiplicity = 1
			n_alpha = (n_electrons + multiplicity - 1) // 2
			n_beta = n_electrons - n_alpha
			if n_alpha <= 0 or n_beta <= 0:
				return None

			overlap = vlx.OverlapDriver().compute(molecule, basis).to_numpy()
			alpha_coefficients = self._s_normalize(
				np.array(scf_driver.molecular_orbitals._orbitals[0][:, n_alpha - 1], dtype=float),
				overlap,
			)
			beta_coefficients = self._s_normalize(
				np.array(scf_driver.molecular_orbitals._orbitals[1][:, n_beta - 1], dtype=float),
				overlap,
			)
			ao_to_atom, _ = OrbitalAnalyzer.ao_shell_map(molecule, basis)
			atom_a, atom_b = (int(active_bond[0]), int(active_bond[1]))

			def atom_weight(coefficients, atom):
				indices = np.where(ao_to_atom == atom)[0]
				if len(indices) == 0:
					return 0.0
				projector = np.zeros_like(coefficients)
				projector[indices] = coefficients[indices]
				return abs(float(projector @ overlap @ coefficients))

			if (atom_weight(alpha_coefficients, atom_b) + atom_weight(beta_coefficients, atom_a) >
					atom_weight(alpha_coefficients, atom_a) + atom_weight(beta_coefficients, atom_b)):
				alpha_coefficients, beta_coefficients = beta_coefficients, alpha_coefficients

			return [
				VbOrbital(
					label=f"full_uhf_alpha_frontier_atom_{atom_a + 1}",
					coefficients=alpha_coefficients,
					center=atom_a,
					kind="active",
				),
				VbOrbital(
					label=f"full_uhf_beta_frontier_atom_{atom_b + 1}",
					coefficients=beta_coefficients,
					center=atom_b,
					kind="active",
				),
			]
		except Exception:
			return None

	def _fragment_radical_active_orbitals(self, molecule, basis, active_bond):
		try:
			from .molecule import Molecule
			from .molecularbasis import MolecularBasis
			from .scfunrestdriver import ScfUnrestrictedDriver
			import veloxchem as vlx

			fragments = self._active_bond_fragments(molecule, active_bond)
			if fragments is None or len(fragments) != 2:
				return None
			basis_label = basis.get_label()
			if not basis_label:
				basis_label = basis.get_main_basis_label()
			labels = molecule.get_labels()
			coordinates = molecule.get_coordinates_in_angstrom()
			parent_ao_to_atom, _ = OrbitalAnalyzer.ao_shell_map(molecule, basis)
			overlap = vlx.OverlapDriver().compute(molecule, basis).to_numpy()
			orbitals = []
			for side, anchor_atom, fragment_atoms in zip(('A', 'B'), active_bond, fragments):
				fragment_atoms = tuple(int(atom) for atom in fragment_atoms)
				if int(anchor_atom) not in fragment_atoms:
					return None
				fragment_xyz = "\n".join(
					f"{labels[atom]} {coordinates[atom, 0]:.12f} {coordinates[atom, 1]:.12f} {coordinates[atom, 2]:.12f}"
					for atom in fragment_atoms
				)
				fragment = Molecule.read_str(fragment_xyz)
				fragment.set_charge(0)
				fragment.set_multiplicity(2)
				fragment_basis = MolecularBasis.read(fragment, basis_label, ostream=None)
				scf_driver = ScfUnrestrictedDriver()
				scf_driver.ostream.mute()
				scf_driver.xcfun = "hf"
				scf_driver.guess_unpaired_electrons = "1(1)"
				scf_driver.compute(fragment, fragment_basis)
				n_alpha = int((fragment.number_of_electrons() + 1) // 2)
				if n_alpha <= 0:
					return None
				fragment_coefficients = np.array(
					scf_driver.molecular_orbitals._orbitals[0][:, n_alpha - 1],
					dtype=float,
				)
				fragment_ao_to_atom, _ = OrbitalAnalyzer.ao_shell_map(fragment, fragment_basis)
				coefficients = np.zeros(basis.get_dimensions_of_basis())
				for local_atom, parent_atom in enumerate(fragment_atoms):
					fragment_indices = np.where(fragment_ao_to_atom == local_atom)[0]
					parent_indices = np.where(parent_ao_to_atom == parent_atom)[0]
					if len(fragment_indices) != len(parent_indices):
						return None
					coefficients[parent_indices] = fragment_coefficients[fragment_indices]
				coefficients = self._s_normalize(coefficients, overlap)
				orbitals.append(
					VbOrbital(
						label=f"fragment_{side}_radical_somo_atom_{anchor_atom + 1}",
						coefficients=coefficients,
						center=int(anchor_atom),
						kind="active",
					)
				)
			return orbitals
		except Exception:
			return None

	def _active_bond_fragments(self, molecule, active_bond):
		try:
			connectivity = np.array(molecule.get_connectivity_matrix(), dtype=int)
		except Exception:
			return None
		atom_a, atom_b = (int(active_bond[0]), int(active_bond[1]))
		connectivity[atom_a, atom_b] = 0
		connectivity[atom_b, atom_a] = 0

		def component(start):
			seen = {int(start)}
			stack = [int(start)]
			while stack:
				atom = stack.pop()
				for neighbor in np.where(connectivity[atom] != 0)[0]:
					neighbor = int(neighbor)
					if neighbor not in seen:
						seen.add(neighbor)
						stack.append(neighbor)
			return tuple(sorted(seen))

		fragment_a = component(atom_a)
		fragment_b = component(atom_b)
		if set(fragment_a) & set(fragment_b):
			return None
		return fragment_a, fragment_b

	def _active_bond_distance_angstrom(self, molecule, active_bond):
		try:
			coordinates = molecule.get_coordinates_in_angstrom()
			atom_a, atom_b = int(active_bond[0]), int(active_bond[1])
			return float(np.linalg.norm(coordinates[atom_a] - coordinates[atom_b]))
		except Exception:
			return 0.0

	def _candidate_nao_vector(self, candidate, nao_data):
		if nao_data is None or not candidate.get('coefficients'):
			return None
		vector = np.zeros(len(nao_data.populations))
		coefficients = candidate.get('coefficients', [])
		if isinstance(coefficients, dict):
			for index, value in coefficients.items():
				index = int(index)
				if 0 <= index < len(vector):
					vector[index] = float(value)
		else:
			for item in coefficients:
				if not isinstance(item, dict):
					continue
				index = int(item.get('nao_index', 0)) - 1
				if 0 <= index < len(vector):
					vector[index] = float(item.get('coefficient', 0.0))
		norm = float(np.linalg.norm(vector))
		if norm <= 1.0e-14:
			return None
		return vector / norm

	def _candidate_centered_active_orbitals(self, molecule, basis, active_candidate, nao_data):
		active_bond = tuple(int(atom) for atom in active_candidate.get('atoms', ()))
		if len(active_bond) != 2:
			return None
		candidate_vector = self._candidate_nao_vector(active_candidate, nao_data)
		if candidate_vector is None:
			return None

		import veloxchem as vlx
		S_ao = vlx.OverlapDriver().compute(molecule, basis).to_numpy()
		atom_map = np.array(nao_data.atom_map, dtype=int)
		orbitals = []
		subtype = active_candidate.get('subtype') or 'bond'
		for side, atom in zip(('A', 'B'), active_bond):
			nao_side = np.zeros_like(candidate_vector)
			indices = np.where(atom_map == atom)[0]
			if len(indices) == 0:
				return None
			nao_side[indices] = candidate_vector[indices]
			if np.linalg.norm(nao_side) <= 1.0e-14:
				return None
			coeffs = nao_data.transform @ nao_side
			try:
				coeffs = self._s_normalize(coeffs, S_ao)
			except RuntimeError:
				return None
			orbitals.append(
				VbOrbital(
					label=f"active_{subtype}_{side}_atom_{atom + 1}",
					coefficients=coeffs,
					center=int(atom),
					kind="active",
				)
			)
		return orbitals

	def _default_one_bond_structures(self, active_bond, options):
		atom_a, atom_b = active_bond
		structures = [
			VbStructure(
				label='covalent',
				occupation=((1, 1),),
				spin='singlet',
				charge_pattern={atom_a: 0, atom_b: 0},
			),
		]
		if options.include_ionic:
			structures.extend([
				VbStructure(
					label='ionic_A_minus_B_plus',
					occupation=((2, 0),),
					spin='singlet',
					charge_pattern={atom_a: -1, atom_b: +1},
				),
				VbStructure(
					label='ionic_A_plus_B_minus',
					occupation=((0, 2),),
					spin='singlet',
					charge_pattern={atom_a: +1, atom_b: -1},
				),
			])
		return structures

	def _active_space_diagnostics(self, active_space):
		inactive_sigma = [
			self._candidate_label(c) for c in active_space.inactive_candidates
			if c.get('type') == 'BD' and c.get('subtype') == 'sigma'
		]
		return {
			"message": "Generated VB result from OrbitalAnalyzer active-space candidates.",
			"active_space_model": active_space.metadata.get("model"),
			"active_bond": active_space.active_bond,
			"active_pi_atoms": active_space.metadata.get("active_pi_atoms"),
			"active_orbital_count": len(active_space.active_orbitals),
			"active_electron_count": active_space.metadata.get("electron_count"),
			"active_spin": active_space.metadata.get("spin"),
			"determinant_count": active_space.metadata.get("determinant_count"),
			"active_candidate_label": active_space.active_candidate_label,
			"active_candidate_type": active_space.active_candidate.get("type"),
			"active_candidate_subtype": active_space.active_candidate.get("subtype"),
			"metal_ligand_model": active_space.metadata.get("metal_ligand_model"),
			"metal_ligand_channels": active_space.metadata.get("metal_ligand_channels", ()),
			"metal_ligand_backdonation_enabled": active_space.metadata.get("backdonation_enabled"),
			"metal_ligand_backdonation_blocked": active_space.metadata.get("backdonation_blocked"),
			"selected_metal_ligand_records": active_space.metadata.get("selected_metal_ligand_records", []),
			"component_pi_candidate_labels": active_space.metadata.get("component_pi_candidate_labels", []),
			"active_bond_state_tracked": bool(active_space.metadata.get("active_bond_state_tracked", False)),
			"active_candidate_fallback": bool(active_space.metadata.get("active_candidate_fallback", False)),
			"active_uhf_frontier_orbitals": bool(active_space.metadata.get("active_uhf_frontier_orbitals", False)),
			"active_fragment_radical_orbitals": bool(active_space.metadata.get("active_fragment_radical_orbitals", False)),
			"active_stretched_orbital_source": active_space.metadata.get("active_stretched_orbital_source"),
			"h2_stable_atom_centered_active_space": bool(
				active_space.metadata.get("h2_stable_atom_centered_active_space", False)
			),
			"active_orbital_labels": [orb.label for orb in active_space.active_orbitals],
			"generated_structure_labels": [structure.label for structure in active_space.structures],
			"frozen_candidate_labels": [self._candidate_label(c) for c in active_space.frozen_candidates],
			"inactive_candidate_labels": [self._candidate_label(c) for c in active_space.inactive_candidates],
			"inactive_sigma_candidate_labels": inactive_sigma,
			"excluded_candidate_labels": [self._candidate_label(c) for c in active_space.excluded_candidates],
			"frozen_hf_reference": bool(active_space.metadata.get("frozen_hf_reference", False)),
			"orbital_analysis_source": active_space.metadata.get("source"),
		}

	def _compute_two_orbital_vbci(self,
							   molecule,
							   basis,
							   structures,
							   orbitals,
							   active_space=None,
							   analysis=None,
							   freeze_inactive=False):
		H_ao, S_ao, eri_ao, e_nuc = self._ao_integrals(molecule, basis)
		embedding = None
		if freeze_inactive and active_space is not None and analysis is not None:
			embedding = self._frozen_hf_embedding(
				H_ao,
				S_ao,
				eri_ao,
				e_nuc,
				active_space,
				analysis,
			)
			if embedding is not None:
				H_ao = embedding["h_effective"]
				e_nuc = embedding["constant_energy"]
				orbitals = self._project_active_orbitals_out_of_frozen_density(
					orbitals,
					S_ao,
					embedding.get("frozen_density"),
				)
		C = np.column_stack([orb.coefficients for orb in orbitals])
		S, H = self._build_two_orbital_singlet_matrices(
			structures, C, H_ao, S_ao, eri_ao, e_nuc)
		energy, coeffs, weights, kept, lowdin_weights = self._solve_generalized_vb(
			H, S, return_lowdin=True)
		diagnostics = {
			"message": "Two-electron VB-CI result from spin-adapted singlet structures.",
			"overlap_condition": float(np.linalg.cond(S)),
			"overlap_eigenvalues": np.linalg.eigvalsh(S).tolist(),
			"retained_overlap_rank": int(len(kept)),
			"weight_scheme": "Chirgwin-Coulson",
			"available_weight_schemes": ["Chirgwin-Coulson", "Lowdin"],
		}
		if embedding is not None:
			diagnostics.update({
				"frozen_hf_embedding": True,
				"active_orbitals_orthogonalized_to_frozen_space": True,
				"frozen_electron_count": embedding["frozen_electron_count"],
				"active_reference_electron_count": embedding["active_electron_count"],
				"frozen_constant_energy": embedding["constant_energy"],
				"embedding_model": "inactive HF density = total RHF density minus active reference density",
			})
		return {
			"energy": float(energy),
			"structure_coefficients": coeffs,
			"overlap": S,
			"Hamiltonian": H,
			"weights": weights,
			"lowdin_weights": lowdin_weights,
			"orbitals": [orb.coefficients for orb in orbitals],
			"diagnostics": diagnostics,
		}

	def _compute_compact_csf(self,
						 molecule,
						 basis,
						 structures,
						 orbitals,
						 active_space=None,
						 analysis=None,
						 freeze_inactive=False):
		if active_space is None:
			raise RuntimeError("Compact CSF mode requires a generated active space.")
		electron_count = int(active_space.metadata.get("electron_count") or 0)
		spin = str(active_space.metadata.get("spin") or "singlet")
		active_pi_atoms = tuple(active_space.metadata.get("active_pi_atoms", ()))
		if electron_count != 2 or spin != "singlet" or len(active_pi_atoms) < 3:
			raise NotImplementedError(
				"Compact CSF Hamiltonians are currently implemented only for two-electron singlet multicenter pi active spaces.")

		H_ao, S_ao, eri_ao, e_nuc = self._ao_integrals(molecule, basis)
		embedding = None
		if freeze_inactive and analysis is not None:
			embedding = self._frozen_hf_embedding(
				H_ao,
				S_ao,
				eri_ao,
				e_nuc,
				active_space,
				analysis,
			)
			if embedding is not None:
				H_ao = embedding["h_effective"]
				e_nuc = embedding["constant_energy"]

		coefficients = np.column_stack([orb.coefficients for orb in orbitals])
		S_full, H_full = self._build_two_orbital_singlet_matrices(
			structures,
			coefficients,
			H_ao,
			S_ao,
			eri_ao,
			e_nuc,
		)
		full_energy, full_coeffs, _, full_kept, _ = self._solve_generalized_vb(
			H_full,
			S_full,
			return_lowdin=True,
		)
		template_records = self._compact_csf_template_records(
			structures,
			active_space,
		)
		if not template_records:
			raise RuntimeError("No compact CSF templates were generated for this active space.")

		template_matrix = np.column_stack([record["vector"] for record in template_records])
		S_compact = template_matrix.T @ S_full @ template_matrix
		H_compact = template_matrix.T @ H_full @ template_matrix
		S_compact = 0.5 * (S_compact + S_compact.T)
		H_compact = 0.5 * (H_compact + H_compact.T)
		energy, compact_coeffs, compact_weights, kept, lowdin_weights = self._solve_generalized_vb(
			H_compact,
			S_compact,
			return_lowdin=True,
		)

		full_metric_rhs = template_matrix.T @ S_full @ full_coeffs
		try:
			projection_coeffs = np.linalg.pinv(S_compact) @ full_metric_rhs
		except np.linalg.LinAlgError:
			projection_coeffs = np.zeros(len(template_records))
		projected_full = template_matrix @ projection_coeffs
		captured_weight = float(projected_full.T @ S_full @ projected_full)
		captured_weight = max(0.0, min(captured_weight, 1.0 + 1.0e-8))

		details = []
		for index, record in enumerate(template_records):
			details.append({
				"label": record["label"],
				"type": record["type"],
				"bond_pairs": record.get("bond_pairs", ()),
				"coefficient": float(compact_coeffs[index]),
				"weight": float(compact_weights[index]),
				"lowdin_weight": float(lowdin_weights[index]),
				"determinant_expansion": record["structure_expansion"],
			})

		diagnostics = {
			"message": "Compact spin-adapted CSF Hamiltonian result from graph resonance templates.",
			"compact_csf_model": "two-electron-singlet-pi-adjacent-bond-templates",
			"compact_csf_labels": [record["label"] for record in template_records],
			"compact_csf_types": [record["type"] for record in template_records],
			"compact_csf_count": len(template_records),
			"compact_csf_details": details,
			"compact_csf_overlap_eigenvalues": np.linalg.eigvalsh(S_compact).tolist(),
			"compact_csf_retained_rank": int(len(kept)),
			"compact_csf_captured_subspace_weight": captured_weight,
			"compact_csf_full_reference_energy": float(full_energy),
			"compact_csf_energy_error_to_full_reference": float(energy - full_energy),
			"full_reference_retained_overlap_rank": int(len(full_kept)),
			"overlap_condition": float(np.linalg.cond(S_compact)),
			"overlap_eigenvalues": np.linalg.eigvalsh(S_compact).tolist(),
			"retained_overlap_rank": int(len(kept)),
			"weight_scheme": "Chirgwin-Coulson compact CSF weights",
			"available_weight_schemes": ["Chirgwin-Coulson", "Lowdin"],
		}
		if embedding is not None:
			diagnostics.update({
				"frozen_hf_embedding": True,
				"frozen_electron_count": embedding["frozen_electron_count"],
				"active_reference_electron_count": embedding["active_electron_count"],
				"frozen_constant_energy": embedding["constant_energy"],
				"embedding_model": "inactive HF density = total RHF density minus active reference density",
			})
		return {
			"energy": float(energy),
			"structure_coefficients": compact_coeffs,
			"overlap": S_compact,
			"Hamiltonian": H_compact,
			"weights": compact_weights,
			"lowdin_weights": lowdin_weights,
			"orbitals": [orb.coefficients for orb in orbitals],
			"diagnostics": diagnostics,
		}

	def _compute_compact_csf_bovb(self,
							  molecule,
							  basis,
							  structures,
							  orbitals,
							  options,
							  active_space=None,
							  analysis=None,
							  freeze_inactive=False):
		from scipy.optimize import minimize

		if active_space is None:
			raise RuntimeError("Compact-CSF BOVB mode requires a generated active space.")
		electron_count = int(active_space.metadata.get("electron_count") or 0)
		spin = str(active_space.metadata.get("spin") or "singlet")
		active_pi_atoms = tuple(active_space.metadata.get("active_pi_atoms", ()))
		if electron_count != 2 or spin != "singlet" or len(active_pi_atoms) != 3:
			raise NotImplementedError(
				"Compact-CSF BOVB is currently implemented only for allyl-cation-like two-electron three-center pi singlet spaces.")

		H_ao, S_ao, eri_ao, e_nuc = self._ao_integrals(molecule, basis)
		embedding = None
		if freeze_inactive and analysis is not None:
			embedding = self._frozen_hf_embedding(
				H_ao,
				S_ao,
				eri_ao,
				e_nuc,
				active_space,
				analysis,
			)
			if embedding is not None:
				H_ao = embedding["h_effective"]
				e_nuc = embedding["constant_energy"]

		coefficients = np.column_stack([orb.coefficients for orb in orbitals])
		active_vectors = tuple(coefficients[:, index] for index in range(coefficients.shape[1]))
		breathing_vectors = []
		for orbital in orbitals:
			breathing_vectors.append(
				self._center_breathing_vector(
					molecule,
					basis,
					orbital.center,
					orbital.coefficients,
					active_vectors,
					S_ao,
				)
			)
		has_external_breathing = all(vector is not None for vector in breathing_vectors)
		breathing_vectors = [
			vector if vector is not None else np.zeros(coefficients.shape[0])
			for vector in breathing_vectors
		]

		template_records = self._compact_csf_template_records(
			structures,
			active_space,
		)
		if not template_records:
			raise RuntimeError("No compact CSF templates were generated for this active space.")
		template_matrix = np.column_stack([record["vector"] for record in template_records])
		atom_position = {int(atom): index for index, atom in enumerate(active_pi_atoms)}

		def coefficient_sets(parameters):
			sets = [np.array(coefficients, dtype=float, copy=True) for _ in structures]
			modified = [set() for _ in structures]
			for template_index, record in enumerate(template_records):
				amplitude = float(parameters[template_index])
				if abs(amplitude) <= 1.0e-14:
					continue
				for expansion in record["structure_expansion"]:
					structure_index = int(expansion["structure_index"])
					for atom_i, atom_j in record.get("bond_pairs", ()): 
						for atom in (atom_i, atom_j):
							orbital_index = atom_position[int(atom)]
							sets[structure_index][:, orbital_index] += (
								amplitude * breathing_vectors[orbital_index]
							)
							modified[structure_index].add(orbital_index)
			for structure_index, matrix in enumerate(sets):
				for orbital_index in modified[structure_index]:
					try:
						matrix[:, orbital_index] = self._s_normalize(
							matrix[:, orbital_index],
							S_ao,
						)
					except RuntimeError:
						sets[structure_index] = np.array(coefficients, dtype=float, copy=True)
						break
			return sets

		def compact_matrices(parameters):
			S_full, H_full = self._build_structure_specific_two_orbital_singlet_matrices(
				structures,
				coefficient_sets(parameters),
				H_ao,
				S_ao,
				eri_ao,
				e_nuc,
			)
			S_compact = template_matrix.T @ S_full @ template_matrix
			H_compact = template_matrix.T @ H_full @ template_matrix
			return 0.5 * (S_compact + S_compact.T), 0.5 * (H_compact + H_compact.T)

		def energy_for_parameters(parameters):
			try:
				S_compact, H_compact = compact_matrices(parameters)
				energy, _, _, _ = self._solve_generalized_vb(H_compact, S_compact)
			except Exception:
				energy = 1.0e6
			return energy

		initial_parameters = np.zeros(len(template_records))
		initial_energy = energy_for_parameters(initial_parameters)
		result = minimize(
			energy_for_parameters,
			x0=initial_parameters,
			method='L-BFGS-B',
			bounds=[(-0.35, 0.35)] * len(template_records),
			options={
				'ftol': options.conv_thresh if options else 1.0e-8,
				'gtol': options.conv_thresh if options else 1.0e-8,
				'maxiter': options.max_iter if options else 50,
			},
		)
		optimizer_energy = float(result.fun) if np.isfinite(result.fun) else 1.0e6
		unphysical_non_h2_lowering = (
			not self._is_h2_molecule(molecule) and
			initial_energy - optimizer_energy > 0.50
		)
		used_fixed_limit = optimizer_energy > initial_energy + 1.0e-10 or unphysical_non_h2_lowering
		parameters = initial_parameters if used_fixed_limit else np.array(result.x, dtype=float)
		S_compact, H_compact = compact_matrices(parameters)
		energy, compact_coeffs, compact_weights, kept, lowdin_weights = self._solve_generalized_vb(
			H_compact,
			S_compact,
			return_lowdin=True,
		)

		S_fixed_full, H_fixed_full = self._build_two_orbital_singlet_matrices(
			structures,
			coefficients,
			H_ao,
			S_ao,
			eri_ao,
			e_nuc,
		)
		full_energy, full_coeffs, _, full_kept, _ = self._solve_generalized_vb(
			H_fixed_full,
			S_fixed_full,
			return_lowdin=True,
		)
		S_fixed_compact = template_matrix.T @ S_fixed_full @ template_matrix
		full_metric_rhs = template_matrix.T @ S_fixed_full @ full_coeffs
		try:
			projection_coeffs = np.linalg.pinv(S_fixed_compact) @ full_metric_rhs
		except np.linalg.LinAlgError:
			projection_coeffs = np.zeros(len(template_records))
		projected_full = template_matrix @ projection_coeffs
		captured_weight = float(projected_full.T @ S_fixed_full @ projected_full)
		captured_weight = max(0.0, min(captured_weight, 1.0 + 1.0e-8))

		details = []
		for index, record in enumerate(template_records):
			details.append({
				"label": record["label"],
				"type": record["type"],
				"bond_pairs": record.get("bond_pairs", ()),
				"coefficient": float(compact_coeffs[index]),
				"weight": float(compact_weights[index]),
				"lowdin_weight": float(lowdin_weights[index]),
				"breathing_amplitude": float(parameters[index]),
				"determinant_expansion": record["structure_expansion"],
			})

		diagnostics = {
			"message": "Compact-CSF BOVB result with center-local structure-specific pi breathing orbitals.",
			"compact_csf_bovb_model": "allyl-cation-two-electron-pi-center-local-breathing",
			"compact_csf_model": "two-electron-singlet-pi-adjacent-bond-templates",
			"compact_csf_labels": [record["label"] for record in template_records],
			"compact_csf_types": [record["type"] for record in template_records],
			"compact_csf_count": len(template_records),
			"compact_csf_details": details,
			"compact_csf_bovb_breathing_amplitudes": [float(value) for value in parameters],
			"compact_csf_bovb_has_external_breathing_space": bool(has_external_breathing),
			"compact_csf_bovb_initial_energy": float(initial_energy),
			"compact_csf_bovb_optimizer_energy": float(optimizer_energy),
			"compact_csf_bovb_energy_lowering": float(initial_energy - energy),
			"compact_csf_bovb_used_fixed_orbital_limit": bool(used_fixed_limit),
			"compact_csf_bovb_optimizer_success": bool(result.success),
			"compact_csf_bovb_optimizer_message": str(result.message),
			"compact_csf_overlap_eigenvalues": np.linalg.eigvalsh(S_compact).tolist(),
			"compact_csf_retained_rank": int(len(kept)),
			"compact_csf_captured_subspace_weight": captured_weight,
			"compact_csf_full_reference_energy": float(full_energy),
			"compact_csf_energy_error_to_full_reference": float(energy - full_energy),
			"full_reference_retained_overlap_rank": int(len(full_kept)),
			"overlap_condition": float(np.linalg.cond(S_compact)),
			"overlap_eigenvalues": np.linalg.eigvalsh(S_compact).tolist(),
			"retained_overlap_rank": int(len(kept)),
			"weight_scheme": "Chirgwin-Coulson compact CSF weights",
			"available_weight_schemes": ["Chirgwin-Coulson", "Lowdin"],
			"structure_specific_orbitals": True,
		}
		if embedding is not None:
			diagnostics.update({
				"frozen_hf_embedding": True,
				"frozen_electron_count": embedding["frozen_electron_count"],
				"active_reference_electron_count": embedding["active_electron_count"],
				"frozen_constant_energy": embedding["constant_energy"],
				"embedding_model": "inactive HF density = total RHF density minus active reference density",
			})
		return {
			"energy": float(energy),
			"structure_coefficients": compact_coeffs,
			"overlap": S_compact,
			"Hamiltonian": H_compact,
			"weights": compact_weights,
			"lowdin_weights": lowdin_weights,
			"orbitals": coefficient_sets(parameters),
			"diagnostics": diagnostics,
		}

	def _compact_csf_template_records(self, structures, active_space):
		active_pi_atoms = tuple(active_space.metadata.get("active_pi_atoms", ()))
		electron_count = int(active_space.metadata.get("electron_count") or 0)
		spin = str(active_space.metadata.get("spin") or "singlet")
		if electron_count != 2 or spin != "singlet" or len(active_pi_atoms) < 3:
			return []
		index_by_occupation = {
			tuple(int(value) for value in structure.occupation[0]): index
			for index, structure in enumerate(structures)
			if len(structure.occupation) == 1
		}
		records = []
		for left, right in zip(range(len(active_pi_atoms) - 1), range(1, len(active_pi_atoms))):
			occupation = [0] * len(active_pi_atoms)
			occupation[left] = 1
			occupation[right] = 1
			structure_index = index_by_occupation.get(tuple(occupation))
			if structure_index is None:
				continue
			vector = np.zeros(len(structures))
			vector[structure_index] = 1.0
			atom_left = active_pi_atoms[left]
			atom_right = active_pi_atoms[right]
			records.append({
				"label": f"allyl_cation_pi_bond_atom_{atom_left + 1}_atom_{atom_right + 1}",
				"type": "allyl_cation" if len(active_pi_atoms) == 3 else "two_electron_pi_bond",
				"bond_pairs": ((atom_left, atom_right),),
				"vector": vector,
				"structure_expansion": [{
					"structure_index": int(structure_index),
					"structure_label": structures[structure_index].label,
					"coefficient": 1.0,
				}],
			})
		return records

	def _compute_active_determinant_ci(self,
								 molecule,
								 basis,
								 structures,
								 orbitals,
								 active_space=None,
								 analysis=None,
								 freeze_inactive=False):
		H_ao, S_ao, eri_ao, e_nuc = self._ao_integrals(molecule, basis)
		embedding = None
		if freeze_inactive and active_space is not None and analysis is not None:
			embedding = self._frozen_hf_embedding(
				H_ao,
				S_ao,
				eri_ao,
				e_nuc,
				active_space,
				analysis,
			)
			if embedding is not None:
				H_ao = embedding["h_effective"]
				e_nuc = embedding["constant_energy"]

		coefficients = np.column_stack([orb.coefficients for orb in orbitals])
		orth_coefficients, active_overlap_values = self._orthonormal_active_orbitals(
			coefficients,
			S_ao,
		)
		H_mo = orth_coefficients.T @ H_ao @ orth_coefficients
		eri_mo = np.einsum('up,vq,lr,ms,uvlm->pqrs',
						 orth_coefficients,
						 orth_coefficients,
						 orth_coefficients,
						 orth_coefficients,
						 eri_ao,
						 optimize=True)
		determinants = [self._structure_spin_bitstring(structure) for structure in structures]
		H = self._build_spin_determinant_hamiltonian(
			determinants,
			H_mo,
			eri_mo,
			float(e_nuc),
		)
		S = np.eye(len(determinants))
		energy, coeffs, weights, kept, lowdin_weights, root_info = self._solve_determinant_ci_root(
			H,
			structures,
			active_space,
		)
		resonance_analysis = self._determinant_resonance_analysis(
			structures,
			coeffs,
			weights,
			lowdin_weights,
			active_space,
		)
		chemical_resonance_analysis = self._chemical_resonance_analysis(
			structures,
			coeffs,
			active_space,
		)
		diagnostics = {
			"message": "Fixed-orbital pi determinant-CI result in an orthonormalized active orbital basis.",
			"overlap_condition": float(np.linalg.cond(S)),
			"overlap_eigenvalues": np.linalg.eigvalsh(S).tolist(),
			"active_orbital_overlap_eigenvalues": active_overlap_values.tolist(),
			"retained_overlap_rank": int(len(kept)),
			"determinant_count": len(determinants),
			"n_alpha": active_space.metadata.get("n_alpha") if active_space is not None else None,
			"n_beta": active_space.metadata.get("n_beta") if active_space is not None else None,
			"weight_scheme": "orthonormal determinant squared coefficients",
			"available_weight_schemes": ["orthonormal determinant squared coefficients"],
			"resonance_structure_labels": resonance_analysis["labels"],
			"resonance_structure_weights": resonance_analysis["weights"],
			"resonance_structure_lowdin_weights": resonance_analysis["lowdin_weights"],
			"resonance_structure_count": len(resonance_analysis["labels"]),
			"resonance_structure_details": resonance_analysis["details"],
			"chemical_resonance_model": chemical_resonance_analysis["model"],
			"chemical_resonance_labels": chemical_resonance_analysis["labels"],
			"chemical_resonance_types": chemical_resonance_analysis["types"],
			"chemical_resonance_weights": chemical_resonance_analysis["weights"],
			"chemical_resonance_unsymmetrized_weights": chemical_resonance_analysis["unsymmetrized_weights"],
			"chemical_resonance_projection_weights": chemical_resonance_analysis["projection_weights"],
			"chemical_resonance_unsymmetrized_projection_weights": chemical_resonance_analysis["unsymmetrized_projection_weights"],
			"chemical_resonance_count": len(chemical_resonance_analysis["labels"]),
			"chemical_resonance_weight_sum": chemical_resonance_analysis["weight_sum"],
			"chemical_resonance_projection_weight_sum": chemical_resonance_analysis["projection_weight_sum"],
			"chemical_resonance_subspace_weight": chemical_resonance_analysis["subspace_weight"],
			"chemical_resonance_overlap_eigenvalues": chemical_resonance_analysis["overlap_eigenvalues"],
			"chemical_resonance_retained_rank": chemical_resonance_analysis["retained_rank"],
			"chemical_resonance_symmetry_model": chemical_resonance_analysis["symmetry_model"],
			"chemical_resonance_symmetry_orbits": chemical_resonance_analysis["symmetry_orbits"],
			"chemical_resonance_details": chemical_resonance_analysis["details"],
		}
		if embedding is not None:
			total_energy = (
				float(embedding["reference_total_energy"]) +
				float(energy) -
				float(embedding["active_reference_energy"])
			)
			is_state_tracked = bool(
				active_space is not None and
				active_space.metadata.get("metal_ligand_state_tracked_orbitals", False)
			)
			diagnostics.update({
				"frozen_hf_embedding": True,
				"frozen_electron_count": embedding["frozen_electron_count"],
				"active_reference_electron_count": embedding["active_electron_count"],
				"frozen_constant_energy": embedding["constant_energy"],
				"embedded_active_reference_energy": embedding["active_reference_energy"],
				"reference_total_hf_energy": embedding["reference_total_energy"],
				"embedded_active_space_energy": float(energy),
				"vb_total_energy": float(total_energy),
				"vb_total_energy_model": "reference_total_hf_energy + embedded_active_space_energy - embedded_active_reference_energy",
				"energy_represents_total_dissociation_energy": bool(is_state_tracked),
				"embedding_model": "inactive HF density = total RHF density minus active reference density",
			})
			if not is_state_tracked:
				diagnostics["energy_represents_total_dissociation_energy_reason"] = (
					"Metal-ligand active orbitals were selected independently for this geometry; "
					"use active_metal_ligand_reference_orbitals for a state-tracked scan.")
			else:
				diagnostics["energy"] = float(total_energy)
		diagnostics.update(root_info)
		reported_energy = float(
			diagnostics.get("energy", energy)
		)
		return {
			"energy": reported_energy,
			"structure_coefficients": coeffs,
			"overlap": S,
			"Hamiltonian": H,
			"weights": weights,
			"lowdin_weights": lowdin_weights,
			"orbitals": [orb.coefficients for orb in orbitals],
			"orthonormal_active_orbitals": orth_coefficients,
			"diagnostics": diagnostics,
		}

	def _compute_metal_ligand_vbscf_determinant(self,
									 molecule,
									 basis,
									 structures,
									 orbitals,
									 active_space=None,
									 analysis=None,
									 freeze_inactive=False):
		result = self._compute_active_determinant_ci(
			molecule,
			basis,
			structures,
			orbitals,
			active_space=active_space,
			analysis=analysis,
			freeze_inactive=freeze_inactive,
		)
		diagnostics = result.setdefault("diagnostics", {})
		if active_space is not None:
			diagnostics.update(self._active_space_diagnostics(active_space))
		diagnostics.update({
			"message": "Metal-ligand VB-SCF result with common active orbitals in the selected channel determinant space.",
			"metal_ligand_vbscf_model": "common-active-orbital-determinant-ci-reference",
			"metal_ligand_vbscf_initial_energy": float(result["energy"]),
			"metal_ligand_vbscf_energy_lowering": 0.0,
			"metal_ligand_vbscf_has_external_relaxation_space": False,
			"metal_ligand_vbscf_orbital_amplitudes": [0.0] * len(orbitals),
			"metal_ligand_orbital_relaxation_method": "vbscf",
			"structure_specific_orbitals": False,
		})
		return result

	def _compute_metal_ligand_bovb(self,
								molecule,
								basis,
								structures,
								orbitals,
								options,
								active_space=None,
								analysis=None,
								freeze_inactive=False,
								method_label="bovb"):
		from scipy.optimize import minimize

		if active_space is None:
			raise RuntimeError("Metal-ligand BOVB requires a generated active space.")
		if active_space.metadata.get("metal_ligand_model") is None:
			raise RuntimeError("Metal-ligand BOVB requires metal-ligand active channels.")

		H_ao, S_ao, eri_ao, e_nuc = self._ao_integrals(molecule, basis)
		embedding = None
		if freeze_inactive and analysis is not None:
			embedding = self._frozen_hf_embedding(
				H_ao,
				S_ao,
				eri_ao,
				e_nuc,
				active_space,
				analysis,
			)
			if embedding is not None:
				H_ao = embedding["h_effective"]
				e_nuc = embedding["constant_energy"]

		coefficients = np.column_stack([orb.coefficients for orb in orbitals])
		active_vectors = tuple(coefficients[:, index] for index in range(coefficients.shape[1]))
		breathing_vectors = [
			self._center_breathing_vector(
				molecule,
				basis,
				orbital.center,
				orbital.coefficients,
				active_vectors,
				S_ao,
			)
			for orbital in orbitals
		]
		has_external_breathing = any(vector is not None for vector in breathing_vectors)
		breathing_vectors = [
			vector if vector is not None else np.zeros(coefficients.shape[0])
			for vector in breathing_vectors
		]

		determinants = [self._structure_spin_bitstring(structure) for structure in structures]

		def coefficient_matrix(parameters):
			matrix = np.array(coefficients, dtype=float, copy=True)
			for orbital_index, amplitude in enumerate(parameters):
				if abs(float(amplitude)) <= 1.0e-14:
					continue
				matrix[:, orbital_index] += float(amplitude) * breathing_vectors[orbital_index]
				try:
					matrix[:, orbital_index] = self._s_normalize(
						matrix[:, orbital_index],
						S_ao,
					)
				except RuntimeError:
					return None
			return matrix

		def solve_for_parameters(parameters):
			matrix = coefficient_matrix(parameters)
			if matrix is None:
				raise RuntimeError("Metal-ligand BOVB generated a singular active orbital.")
			orth_coefficients, active_overlap_values = self._orthonormal_active_orbitals(
				matrix,
				S_ao,
			)
			H_mo = orth_coefficients.T @ H_ao @ orth_coefficients
			eri_mo = np.einsum('up,vq,lr,ms,uvlm->pqrs',
							 orth_coefficients,
							 orth_coefficients,
							 orth_coefficients,
							 orth_coefficients,
							 eri_ao,
							 optimize=True)
			H = self._build_spin_determinant_hamiltonian(
				determinants,
				H_mo,
				eri_mo,
				float(e_nuc),
			)
			energy, coeffs, weights, kept, lowdin_weights, root_info = self._solve_determinant_ci_root(
				H,
				structures,
				active_space,
			)
			return {
				"energy": float(energy),
				"coeffs": coeffs,
				"weights": weights,
				"kept": kept,
				"lowdin_weights": lowdin_weights,
				"root_info": root_info,
				"Hamiltonian": H,
				"active_overlap_values": active_overlap_values,
				"orthonormal_active_orbitals": orth_coefficients,
				"coefficient_matrix": matrix,
			}

		def energy_for_parameters(parameters):
			try:
				return solve_for_parameters(parameters)["energy"]
			except Exception:
				return 1.0e6

		initial_parameters = np.zeros(len(orbitals))
		initial_solution = solve_for_parameters(initial_parameters)
		initial_energy = float(initial_solution["energy"])
		if has_external_breathing:
			trial_starts = [initial_parameters]
			for sign in (1.0, -1.0):
				trial = np.zeros(len(orbitals))
				for index, vector in enumerate(breathing_vectors):
					if np.linalg.norm(vector) > 1.0e-14:
						trial[index] = sign * 0.15
				trial_starts.append(trial)
			best_result = None
			best_energy = 1.0e6
			for start in trial_starts:
				result = minimize(
					energy_for_parameters,
					x0=start,
					method='L-BFGS-B',
					bounds=[(-0.35, 0.35)] * len(orbitals),
					options={
						'ftol': options.conv_thresh if options else 1.0e-8,
						'gtol': options.conv_thresh if options else 1.0e-8,
						'maxiter': options.max_iter if options else 50,
					},
				)
				trial_energy = float(result.fun) if np.isfinite(result.fun) else 1.0e6
				if trial_energy < best_energy:
					best_energy = trial_energy
					best_result = result
			optimizer_energy = best_energy
			optimizer_success = bool(best_result.success) if best_result is not None else False
			optimizer_message = str(best_result.message) if best_result is not None else "No optimizer result was produced."
			optimized_parameters = np.array(best_result.x, dtype=float) if best_result is not None else initial_parameters
		else:
			optimizer_energy = initial_energy
			optimizer_success = False
			optimizer_message = "No center-local external breathing directions were available."
			optimized_parameters = initial_parameters

		used_fixed_limit = optimizer_energy > initial_energy + 1.0e-10
		parameters = initial_parameters if used_fixed_limit else optimized_parameters
		solution = initial_solution if used_fixed_limit else solve_for_parameters(parameters)
		energy = float(solution["energy"])
		coeffs = solution["coeffs"]
		weights = solution["weights"]
		lowdin_weights = solution["lowdin_weights"]
		kept = solution["kept"]
		H = solution["Hamiltonian"]
		S = np.eye(len(determinants))

		resonance_analysis = self._determinant_resonance_analysis(
			structures,
			coeffs,
			weights,
			lowdin_weights,
			active_space,
		)
		chemical_resonance_analysis = self._chemical_resonance_analysis(
			structures,
			coeffs,
			active_space,
		)
		breathing_records = []
		for orbital, amplitude, vector in zip(orbitals, parameters, breathing_vectors):
			breathing_records.append({
				"orbital_label": orbital.label,
				"center": orbital.center,
				"amplitude": float(amplitude),
				"has_external_direction": bool(np.linalg.norm(vector) > 1.0e-14),
			})

		diagnostics = {
			"message": "Metal-ligand BOVB result with center-local channel breathing orbital relaxation.",
			"metal_ligand_orbital_relaxation_method": method_label,
			"metal_ligand_bovb_model": "channel-local-determinant-ci-breathing-orbital-relaxation",
			"metal_ligand_bovb_initial_energy": float(initial_energy),
			"metal_ligand_bovb_optimizer_energy": float(optimizer_energy),
			"metal_ligand_bovb_energy_lowering": float(initial_energy - energy),
			"metal_ligand_bovb_used_fixed_orbital_limit": bool(used_fixed_limit),
			"metal_ligand_bovb_has_external_breathing_space": bool(has_external_breathing),
			"metal_ligand_bovb_breathing_amplitudes": [float(value) for value in parameters],
			"metal_ligand_bovb_breathing_records": breathing_records,
			"metal_ligand_bovb_optimizer_success": bool(optimizer_success),
			"metal_ligand_bovb_optimizer_message": optimizer_message,
			"structure_specific_orbitals": True,
			"overlap_condition": float(np.linalg.cond(S)),
			"overlap_eigenvalues": np.linalg.eigvalsh(S).tolist(),
			"active_orbital_overlap_eigenvalues": solution["active_overlap_values"].tolist(),
			"retained_overlap_rank": int(len(kept)),
			"determinant_count": len(determinants),
			"n_alpha": active_space.metadata.get("n_alpha"),
			"n_beta": active_space.metadata.get("n_beta"),
			"weight_scheme": "orthonormal determinant squared coefficients",
			"available_weight_schemes": ["orthonormal determinant squared coefficients"],
			"resonance_structure_labels": resonance_analysis["labels"],
			"resonance_structure_weights": resonance_analysis["weights"],
			"resonance_structure_lowdin_weights": resonance_analysis["lowdin_weights"],
			"resonance_structure_count": len(resonance_analysis["labels"]),
			"resonance_structure_details": resonance_analysis["details"],
			"chemical_resonance_model": chemical_resonance_analysis["model"],
			"chemical_resonance_labels": chemical_resonance_analysis["labels"],
			"chemical_resonance_types": chemical_resonance_analysis["types"],
			"chemical_resonance_weights": chemical_resonance_analysis["weights"],
			"chemical_resonance_unsymmetrized_weights": chemical_resonance_analysis["unsymmetrized_weights"],
			"chemical_resonance_projection_weights": chemical_resonance_analysis["projection_weights"],
			"chemical_resonance_unsymmetrized_projection_weights": chemical_resonance_analysis["unsymmetrized_projection_weights"],
			"chemical_resonance_count": len(chemical_resonance_analysis["labels"]),
			"chemical_resonance_weight_sum": chemical_resonance_analysis["weight_sum"],
			"chemical_resonance_projection_weight_sum": chemical_resonance_analysis["projection_weight_sum"],
			"chemical_resonance_subspace_weight": chemical_resonance_analysis["subspace_weight"],
			"chemical_resonance_overlap_eigenvalues": chemical_resonance_analysis["overlap_eigenvalues"],
			"chemical_resonance_retained_rank": chemical_resonance_analysis["retained_rank"],
			"chemical_resonance_symmetry_model": chemical_resonance_analysis["symmetry_model"],
			"chemical_resonance_symmetry_orbits": chemical_resonance_analysis["symmetry_orbits"],
			"chemical_resonance_details": chemical_resonance_analysis["details"],
		}
		if embedding is not None:
			total_energy = (
				float(embedding["reference_total_energy"]) +
				float(energy) -
				float(embedding["active_reference_energy"])
			)
			is_state_tracked = bool(
				active_space is not None and
				active_space.metadata.get("metal_ligand_state_tracked_orbitals", False)
			)
			diagnostics.update({
				"frozen_hf_embedding": True,
				"frozen_electron_count": embedding["frozen_electron_count"],
				"active_reference_electron_count": embedding["active_electron_count"],
				"frozen_constant_energy": embedding["constant_energy"],
				"embedded_active_reference_energy": embedding["active_reference_energy"],
				"reference_total_hf_energy": embedding["reference_total_energy"],
				"embedded_active_space_energy": float(energy),
				"vb_total_energy": float(total_energy),
				"vb_total_energy_model": "reference_total_hf_energy + embedded_active_space_energy - embedded_active_reference_energy",
				"energy_represents_total_dissociation_energy": bool(is_state_tracked),
				"embedding_model": "inactive HF density = total RHF density minus active reference density",
			})
			if not is_state_tracked:
				diagnostics["energy_represents_total_dissociation_energy_reason"] = (
					"Metal-ligand active orbitals were selected independently for this geometry; "
					"use active_metal_ligand_reference_orbitals for a state-tracked scan.")
			else:
				diagnostics["energy"] = float(total_energy)
		diagnostics.update(solution["root_info"])
		reported_energy = float(
			diagnostics.get("energy", energy)
		)
		return {
			"energy": reported_energy,
			"structure_coefficients": coeffs,
			"overlap": S,
			"Hamiltonian": H,
			"weights": weights,
			"lowdin_weights": lowdin_weights,
			"orbitals": [solution["coefficient_matrix"][:, index] for index in range(solution["coefficient_matrix"].shape[1])],
			"orthonormal_active_orbitals": solution["orthonormal_active_orbitals"],
			"diagnostics": diagnostics,
		}

	def _solve_determinant_ci_root(self, H, structures, active_space):
		H_sym = 0.5 * (H + H.T)
		eigvals, eigvecs = np.linalg.eigh(H_sym)
		root = 0
		spin = active_space.metadata.get("spin") if active_space is not None else None
		n_alpha = active_space.metadata.get("n_alpha") if active_space is not None else None
		n_beta = active_space.metadata.get("n_beta") if active_space is not None else None
		root_selection = "lowest determinant-CI root"
		spin_exchange_parity = None

		if spin == "singlet" and n_alpha == n_beta:
			permutation = self._alpha_beta_exchange_permutation(structures)
			if permutation is not None:
				parities = np.array([
					float(eigvecs[:, index] @ eigvecs[permutation, index])
					for index in range(eigvecs.shape[1])
				])
				candidate_roots = np.where(parities > 0.5)[0]
				if candidate_roots.size:
					root = int(candidate_roots[0])
					root_selection = "lowest alpha/beta exchange-symmetric singlet root"
				spin_exchange_parity = float(parities[root])

		coeffs = eigvecs[:, root].copy()
		if coeffs[np.argmax(np.abs(coeffs))] < 0.0:
			coeffs *= -1.0
		weights = coeffs**2
		lowdin_weights = weights.copy()
		kept = np.arange(len(structures), dtype=int)
		root_info = {
			"determinant_ci_root_index": int(root),
			"determinant_ci_root_energy": float(eigvals[root]),
			"determinant_ci_root_selection": root_selection,
		}
		if spin_exchange_parity is not None:
			root_info["determinant_ci_spin_exchange_parity"] = spin_exchange_parity
		return float(eigvals[root]), coeffs, weights, kept, lowdin_weights, root_info

	def _alpha_beta_exchange_permutation(self, structures):
		index_by_key = {}
		for index, structure in enumerate(structures):
			alpha_occ = tuple(int(value) for value in structure.occupation[0])
			beta_occ = tuple(int(value) for value in structure.occupation[1])
			index_by_key[(alpha_occ, beta_occ)] = index

		permutation = []
		for structure in structures:
			alpha_occ = tuple(int(value) for value in structure.occupation[0])
			beta_occ = tuple(int(value) for value in structure.occupation[1])
			partner = index_by_key.get((beta_occ, alpha_occ))
			if partner is None:
				return None
			permutation.append(partner)
		return np.asarray(permutation, dtype=int)

	def _determinant_resonance_analysis(self,
									 structures,
									 coeffs,
									 weights,
									 lowdin_weights,
									 active_space):
		active_pi_atoms = tuple(active_space.metadata.get("active_pi_atoms", ())) if active_space is not None else tuple()
		spin = active_space.metadata.get("spin", "singlet") if active_space is not None else "singlet"
		groups = {}
		for index, structure in enumerate(structures):
			alpha_occ = tuple(int(value) for value in structure.occupation[0])
			beta_occ = tuple(int(value) for value in structure.occupation[1])
			spatial_occupation = tuple(alpha_occ[i] + beta_occ[i] for i in range(len(alpha_occ)))
			label = self._spatial_occupation_label(active_pi_atoms, spatial_occupation, spin)
			if label not in groups:
				groups[label] = {
					"label": label,
					"spatial_occupation": spatial_occupation,
					"determinant_labels": [],
					"determinant_indices": [],
					"weight": 0.0,
					"lowdin_weight": 0.0,
					"coefficient_norm": 0.0,
				}
			groups[label]["determinant_labels"].append(structure.label)
			groups[label]["determinant_indices"].append(index)
			groups[label]["weight"] += float(weights[index])
			groups[label]["lowdin_weight"] += float(lowdin_weights[index])
			groups[label]["coefficient_norm"] += float(coeffs[index]**2)

		ordered = sorted(groups.values(), key=lambda item: abs(item["weight"]), reverse=True)
		return {
			"labels": [item["label"] for item in ordered],
			"weights": [item["weight"] for item in ordered],
			"lowdin_weights": [item["lowdin_weight"] for item in ordered],
			"details": ordered,
		}

	def _chemical_resonance_analysis(self, structures, coeffs, active_space):
		if active_space is None:
			return self._empty_chemical_resonance_analysis("no active space")
		active_pi_atoms = tuple(active_space.metadata.get("active_pi_atoms", ()))
		electron_count = active_space.metadata.get("electron_count")
		spin = active_space.metadata.get("spin", "singlet")
		try:
			electron_count = int(electron_count)
		except Exception:
			return self._empty_chemical_resonance_analysis("missing active electron count")
		templates = self._chemical_resonance_templates(
			active_pi_atoms,
			electron_count,
			spin,
		)
		if not templates:
			return self._empty_chemical_resonance_analysis(
				"no graph resonance templates for this active space")

		index_by_occupation = {
			self._structure_occupation_key(structure): index
			for index, structure in enumerate(structures)
		}
		template_records = []
		for template in templates:
			vector = np.zeros(len(structures))
			expansion = []
			for alpha_occ, beta_occ, coefficient in template["terms"]:
				key = (tuple(alpha_occ), tuple(beta_occ))
				index = index_by_occupation.get(key)
				if index is None:
					continue
				vector[index] += float(coefficient)
			norm = float(np.linalg.norm(vector))
			if norm <= 1.0e-14:
				continue
			vector /= norm
			for index, coefficient in enumerate(vector):
				if abs(coefficient) > 1.0e-12:
					expansion.append({
						"determinant_index": int(index),
						"determinant_label": structures[index].label,
						"coefficient": float(coefficient),
					})
			template_records.append({
				"label": template["label"],
				"type": template["type"],
				"bond_pairs": template.get("bond_pairs", ()),
				"special_atom": template.get("special_atom"),
				"vector": vector,
				"determinant_expansion": expansion,
			})
		if not template_records:
			return self._empty_chemical_resonance_analysis(
				"no graph resonance templates with determinant support")

		template_matrix = np.column_stack([record["vector"] for record in template_records])
		metric = template_matrix.T @ template_matrix
		metric = 0.5 * (metric + metric.T)
		rhs = template_matrix.T @ coeffs
		try:
			metric_values, metric_vectors = np.linalg.eigh(metric)
		except np.linalg.LinAlgError:
			return self._empty_chemical_resonance_analysis(
				"chemical resonance template metric diagonalization failed")
		keep = metric_values > 1.0e-10
		if not np.any(keep):
			return self._empty_chemical_resonance_analysis(
				"chemical resonance template metric has no positive-norm subspace")

		metric_inverse = (
			(metric_vectors[:, keep] / metric_values[keep]) @
			metric_vectors[:, keep].T
		)
		template_coefficients = metric_inverse @ rhs
		projected_vector = template_matrix @ template_coefficients
		subspace_weight = float(projected_vector.T @ projected_vector)
		lowdin_coefficients = (
			(metric_vectors[:, keep] * np.sqrt(metric_values[keep])) @
			(metric_vectors[:, keep].T @ template_coefficients)
		)
		lowdin_weights = lowdin_coefficients**2
		lowdin_weight_sum = float(np.sum(lowdin_weights))
		if np.isfinite(lowdin_weight_sum) and lowdin_weight_sum > 0.0:
			lowdin_weights = lowdin_weights / lowdin_weight_sum
		raw_projection_coefficients = rhs
		raw_projection_weights = raw_projection_coefficients**2

		details = []
		for index, record in enumerate(template_records):
			details.append({
				"label": record["label"],
				"type": record["type"],
				"bond_pairs": record["bond_pairs"],
				"special_atom": record["special_atom"],
				"coefficient": float(template_coefficients[index]),
				"lowdin_coefficient": float(lowdin_coefficients[index]),
				"weight": float(lowdin_weights[index]),
				"projection_coefficient": float(raw_projection_coefficients[index]),
				"projection_weight": float(raw_projection_weights[index]),
				"determinant_expansion": record["determinant_expansion"],
			})

		symmetry_analysis = self._chemical_resonance_symmetry_analysis(
			details,
			active_pi_atoms,
		)
		for index, detail in enumerate(details):
			detail["unsymmetrized_weight"] = detail["weight"]
			detail["unsymmetrized_projection_weight"] = detail["projection_weight"]
			detail["symmetry_orbit"] = symmetry_analysis["orbit_by_index"][index]
			detail["weight"] = symmetry_analysis["symmetrized_weights"][index]
			detail["projection_weight"] = symmetry_analysis["symmetrized_projection_weights"][index]

		ordered = sorted(details, key=lambda item: abs(item["weight"]), reverse=True)
		return {
			"model": "graph spin-adapted CSF template projection over determinant-CI wavefunction",
			"labels": [item["label"] for item in ordered],
			"types": [item["type"] for item in ordered],
			"weights": [item["weight"] for item in ordered],
			"unsymmetrized_weights": [item["unsymmetrized_weight"] for item in ordered],
			"projection_weights": [item["projection_weight"] for item in ordered],
			"unsymmetrized_projection_weights": [item["unsymmetrized_projection_weight"] for item in ordered],
			"weight_sum": float(sum(item["weight"] for item in ordered)),
			"projection_weight_sum": float(sum(item["projection_weight"] for item in ordered)),
			"subspace_weight": subspace_weight,
			"overlap_eigenvalues": metric_values.tolist(),
			"retained_rank": int(np.count_nonzero(keep)),
			"symmetry_model": symmetry_analysis["model"],
			"symmetry_orbits": symmetry_analysis["orbits"],
			"details": ordered,
		}

	def _chemical_resonance_symmetry_analysis(self, details, active_pi_atoms):
		if not details:
			return {
				"model": "no chemical resonance templates",
				"symmetrized_weights": [],
				"symmetrized_projection_weights": [],
				"orbit_by_index": [],
				"orbits": [],
			}

		automorphisms = self._active_pi_graph_automorphisms(len(active_pi_atoms))
		if len(automorphisms) == 1:
			weights = [float(item["weight"]) for item in details]
			projection_weights = [float(item["projection_weight"]) for item in details]
			orbits = []
			for index, item in enumerate(details):
				orbits.append({
					"label": f"symmetry_orbit_{index + 1}",
					"type": item["type"],
					"member_count": 1,
					"member_labels": [item["label"]],
					"weight": weights[index],
					"projection_weight": projection_weights[index],
				})
			return {
				"model": "identity graph symmetry",
				"symmetrized_weights": weights,
				"symmetrized_projection_weights": projection_weights,
				"orbit_by_index": [record["label"] for record in orbits],
				"orbits": orbits,
			}

		groups = {}
		keys = []
		for index, detail in enumerate(details):
			key = min(
				self._chemical_resonance_template_symmetry_key(
					detail,
					active_pi_atoms,
					automorphism,
				)
				for automorphism in automorphisms
			)
			keys.append(key)
			groups.setdefault(key, []).append(index)

		symmetrized_weights = [0.0] * len(details)
		symmetrized_projection_weights = [0.0] * len(details)
		orbit_by_index = [""] * len(details)
		orbits = []
		for orbit_index, key in enumerate(sorted(groups), start=1):
			indices = groups[key]
			weight = float(sum(details[index]["weight"] for index in indices))
			projection_weight = float(sum(details[index]["projection_weight"] for index in indices))
			mean_weight = weight / len(indices)
			mean_projection_weight = projection_weight / len(indices)
			label = f"symmetry_orbit_{orbit_index}"
			for index in indices:
				symmetrized_weights[index] = mean_weight
				symmetrized_projection_weights[index] = mean_projection_weight
				orbit_by_index[index] = label
			orbits.append({
				"label": label,
				"type": details[indices[0]]["type"],
				"member_count": len(indices),
				"member_labels": [details[index]["label"] for index in indices],
				"weight": weight,
				"projection_weight": projection_weight,
				"mean_weight": mean_weight,
				"mean_projection_weight": mean_projection_weight,
			})

		return {
			"model": "graph automorphism averaged resonance weights",
			"symmetrized_weights": symmetrized_weights,
			"symmetrized_projection_weights": symmetrized_projection_weights,
			"orbit_by_index": orbit_by_index,
			"orbits": orbits,
		}

	def _active_pi_graph_automorphisms(self, n_active):
		if n_active in (3, 4):
			return [
				tuple(range(n_active)),
				tuple(reversed(range(n_active))),
			]
		if n_active == 6:
			automorphisms = []
			for shift in range(n_active):
				automorphisms.append(
					tuple((site + shift) % n_active for site in range(n_active)))
				automorphisms.append(
					tuple((shift - site) % n_active for site in range(n_active)))
			return automorphisms
		return [tuple(range(n_active))]

	def _chemical_resonance_template_symmetry_key(self,
											 detail,
											 active_pi_atoms,
											 automorphism):
		atom_to_position = {
			int(atom): position for position, atom in enumerate(active_pi_atoms)
		}
		special_atom = detail.get("special_atom")
		if special_atom is None:
			special_position = None
		else:
			special_position = automorphism[atom_to_position[int(special_atom)]]
		mapped_pairs = []
		for atom_i, atom_j in detail.get("bond_pairs", ()): 
			pos_i = automorphism[atom_to_position[int(atom_i)]]
			pos_j = automorphism[atom_to_position[int(atom_j)]]
			mapped_pairs.append(tuple(sorted((pos_i, pos_j))))
		return (
			detail.get("type"),
			special_position,
			tuple(sorted(mapped_pairs)),
		)

	def _empty_chemical_resonance_analysis(self, reason):
		return {
			"model": reason,
			"labels": [],
			"types": [],
			"weights": [],
			"unsymmetrized_weights": [],
			"projection_weights": [],
			"unsymmetrized_projection_weights": [],
			"weight_sum": 0.0,
			"projection_weight_sum": 0.0,
			"subspace_weight": 0.0,
			"overlap_eigenvalues": [],
			"retained_rank": 0,
			"symmetry_model": "no chemical resonance templates",
			"symmetry_orbits": [],
			"details": [],
		}

	def _structure_occupation_key(self, structure):
		return (
			tuple(int(value) for value in structure.occupation[0]),
			tuple(int(value) for value in structure.occupation[1]),
		)

	def _chemical_resonance_templates(self, active_pi_atoms, electron_count, spin):
		n_active = len(active_pi_atoms)
		if n_active == 3 and electron_count == 3 and spin == 'doublet':
			return self._allyl_radical_csf_templates(active_pi_atoms)
		if n_active == 3 and electron_count == 4 and spin == 'singlet':
			return self._allyl_anion_csf_templates(active_pi_atoms)
		if spin == 'singlet' and electron_count == n_active and n_active % 2 == 0:
			return self._closed_shell_pairing_csf_templates(active_pi_atoms)
		return []

	def _allyl_radical_csf_templates(self, active_pi_atoms):
		n_active = len(active_pi_atoms)
		templates = []
		for radical_site in range(n_active):
			paired_sites = tuple(site for site in range(n_active) if site != radical_site)
			terms = self._spin_coupled_terms(
				n_active,
				pairs=(paired_sites,),
				unpaired_alpha=(radical_site,),
			)
			templates.append({
				"label": (
					f"allyl_radical_atom_{active_pi_atoms[radical_site] + 1}__"
					f"pi_bond_atom_{active_pi_atoms[paired_sites[0]] + 1}_atom_{active_pi_atoms[paired_sites[1]] + 1}"
				),
				"type": "allyl_radical",
				"bond_pairs": ((active_pi_atoms[paired_sites[0]], active_pi_atoms[paired_sites[1]]),),
				"special_atom": active_pi_atoms[radical_site],
				"terms": terms,
			})
		return templates

	def _allyl_anion_csf_templates(self, active_pi_atoms):
		n_active = len(active_pi_atoms)
		templates = []
		for lone_pair_site in range(n_active):
			paired_sites = tuple(site for site in range(n_active) if site != lone_pair_site)
			terms = self._spin_coupled_terms(
				n_active,
				pairs=(paired_sites,),
				doubly_occupied=(lone_pair_site,),
			)
			templates.append({
				"label": (
					f"allyl_anion_lone_pair_atom_{active_pi_atoms[lone_pair_site] + 1}__"
					f"pi_bond_atom_{active_pi_atoms[paired_sites[0]] + 1}_atom_{active_pi_atoms[paired_sites[1]] + 1}"
				),
				"type": "allyl_anion",
				"bond_pairs": ((active_pi_atoms[paired_sites[0]], active_pi_atoms[paired_sites[1]]),),
				"special_atom": active_pi_atoms[lone_pair_site],
				"terms": terms,
			})
		return templates

	def _closed_shell_pairing_csf_templates(self, active_pi_atoms):
		n_active = len(active_pi_atoms)
		templates = []
		for pairs in self._perfect_matchings(tuple(range(n_active))):
			resonance_type = self._pairing_resonance_type(n_active, pairs)
			pair_label = "__".join(
				f"atom_{active_pi_atoms[i] + 1}_atom_{active_pi_atoms[j] + 1}"
				for i, j in pairs
			)
			templates.append({
				"label": f"{resonance_type}_pi_bonds_{pair_label}",
				"type": resonance_type,
				"bond_pairs": tuple((active_pi_atoms[i], active_pi_atoms[j]) for i, j in pairs),
				"terms": self._spin_coupled_terms(n_active, pairs=tuple(pairs)),
			})
		return templates

	def _perfect_matchings(self, sites):
		if not sites:
			return [tuple()]
		first = sites[0]
		matchings = []
		for offset in range(1, len(sites)):
			second = sites[offset]
			rest = sites[1:offset] + sites[offset + 1:]
			for matching in self._perfect_matchings(rest):
				matchings.append(((first, second),) + matching)
		return matchings

	def _pairing_resonance_type(self, n_active, pairs):
		if n_active == 4:
			path_edges = {tuple(sorted(edge)) for edge in ((0, 1), (1, 2), (2, 3))}
			pair_edges = {tuple(sorted(edge)) for edge in pairs}
			if pair_edges == {tuple(sorted(edge)) for edge in ((0, 1), (2, 3))}:
				return "butadiene_kekule"
			if pair_edges.issubset(path_edges):
				return "butadiene_covalent_pairing"
			return "butadiene_long_bond_pairing"
		if n_active == 6:
			cycle_edges = {tuple(sorted((i, (i + 1) % n_active))) for i in range(n_active)}
			pair_edges = {tuple(sorted(edge)) for edge in pairs}
			if pair_edges.issubset(cycle_edges):
				return "benzene_kekule"
			return "benzene_dewar"
		return "covalent_pairing"

	def _spin_coupled_terms(self,
						 n_active,
						 pairs=(),
						 doubly_occupied=(),
						 unpaired_alpha=()):
		# The determinant-CI basis used by VbDriver is an alpha-occupation string
		# and a beta-occupation string, not an ordered second-quantized creator
		# string. In this representation a singlet pair (i,j) is the symmetric
		# alpha/beta occupation combination
		#     |i_alpha j_beta> + |j_alpha i_beta>
		# divided by sqrt(2). Applying second-quantized creation signs here would
		# double count phase factors already implicit in the determinant ordering
		# and makes the butadiene covalent CSFs orthogonal to the determinant-CI
		# root.
		configurations = [([0] * n_active, [0] * n_active, 1.0)]
		for site in doubly_occupied:
			for alpha_occ, beta_occ, _ in configurations:
				alpha_occ[site] += 1
				beta_occ[site] += 1
		for site in unpaired_alpha:
			for alpha_occ, _, _ in configurations:
				alpha_occ[site] += 1

		factor = 1.0 / np.sqrt(2.0)
		for i, j in pairs:
			next_configurations = []
			for alpha_occ, beta_occ, coefficient in configurations:
				alpha_ij = list(alpha_occ)
				beta_ij = list(beta_occ)
				alpha_ij[i] += 1
				beta_ij[j] += 1
				next_configurations.append((alpha_ij, beta_ij, coefficient * factor))

				alpha_ji = list(alpha_occ)
				beta_ji = list(beta_occ)
				alpha_ji[j] += 1
				beta_ji[i] += 1
				next_configurations.append((alpha_ji, beta_ji, coefficient * factor))
			configurations = next_configurations

		combined = {}
		for alpha_occ, beta_occ, coefficient in configurations:
			if any(value > 1 for value in alpha_occ):
				continue
			if any(value > 1 for value in beta_occ):
				continue
			key = (tuple(alpha_occ), tuple(beta_occ))
			combined[key] = combined.get(key, 0.0) + float(coefficient)
		return [
			(alpha_occ, beta_occ, coefficient)
			for (alpha_occ, beta_occ), coefficient in combined.items()
			if abs(coefficient) > 1.0e-14
		]

	def _orthonormal_active_orbitals(self, coefficients, s_ao):
		S_active = coefficients.T @ s_ao @ coefficients
		values, vectors = np.linalg.eigh(0.5 * (S_active + S_active.T))
		keep = values > 1.0e-10
		if not np.all(keep):
			raise RuntimeError("Active pi orbital overlap matrix is rank deficient.")
		orthogonalizer = vectors[:, keep] / np.sqrt(values[keep])
		return coefficients @ orthogonalizer, values

	def _structure_spin_bitstring(self, structure):
		if len(structure.occupation) != 2:
			raise RuntimeError("Determinant-CI structures require alpha and beta occupation strings.")
		alpha_occ = tuple(int(value) for value in structure.occupation[0])
		beta_occ = tuple(int(value) for value in structure.occupation[1])
		if len(alpha_occ) != len(beta_occ):
			raise RuntimeError("Alpha and beta occupation strings must have the same length.")
		n_active = len(alpha_occ)
		bitstring = 0
		for index, occupation in enumerate(alpha_occ):
			if occupation:
				bitstring |= (1 << index)
		for index, occupation in enumerate(beta_occ):
			if occupation:
				bitstring |= (1 << (n_active + index))
		return bitstring

	def _build_spin_determinant_hamiltonian(self, determinants, h_mo, eri_mo, constant_energy):
		n_spatial = h_mo.shape[0]
		n_spin = 2 * n_spatial
		H = np.zeros((len(determinants), len(determinants)))
		index_by_det = {det: idx for idx, det in enumerate(determinants)}
		for col, determinant in enumerate(determinants):
			for q in range(n_spin):
				annihilated = self._annihilate_spin_orbital(determinant, q)
				if annihilated is None:
					continue
				det_q, sign_q = annihilated
				q_spatial = q % n_spatial
				q_spin = q // n_spatial
				for p in range(n_spin):
					if p // n_spatial != q_spin:
						continue
					created = self._create_spin_orbital(det_q, p)
					if created is None:
						continue
					det_pq, sign_p = created
					row = index_by_det.get(det_pq)
					if row is not None:
						H[row, col] += sign_q * sign_p * h_mo[p % n_spatial, q_spatial]
			for r in range(n_spin):
				annihilated_r = self._annihilate_spin_orbital(determinant, r)
				if annihilated_r is None:
					continue
				det_r, sign_r = annihilated_r
				for s in range(n_spin):
					annihilated_s = self._annihilate_spin_orbital(det_r, s)
					if annihilated_s is None:
						continue
					det_rs, sign_s = annihilated_s
					for q in range(n_spin):
						created_q = self._create_spin_orbital(det_rs, q)
						if created_q is None:
							continue
						det_qrs, sign_q = created_q
						for p in range(n_spin):
							g = self._antisymmetrized_spin_eri(p, q, r, s, eri_mo)
							if abs(g) <= 1.0e-14:
								continue
							created_p = self._create_spin_orbital(det_qrs, p)
							if created_p is None:
								continue
							det_pqrs, sign_p = created_p
							row = index_by_det.get(det_pqrs)
							if row is not None:
								H[row, col] += 0.5 * sign_r * sign_s * sign_q * sign_p * g
		for idx, determinant in enumerate(determinants):
			H[idx, idx] += constant_energy
		return 0.5 * (H + H.T)

	def _annihilate_spin_orbital(self, determinant, orbital):
		mask = 1 << orbital
		if determinant & mask == 0:
			return None
		sign = -1 if ((determinant & (mask - 1)).bit_count() % 2) else 1
		return determinant ^ mask, sign

	def _create_spin_orbital(self, determinant, orbital):
		mask = 1 << orbital
		if determinant & mask:
			return None
		sign = -1 if ((determinant & (mask - 1)).bit_count() % 2) else 1
		return determinant | mask, sign

	def _antisymmetrized_spin_eri(self, p, q, r, s, eri_mo):
		n_spatial = eri_mo.shape[0]
		p_spatial, q_spatial = p % n_spatial, q % n_spatial
		r_spatial, s_spatial = r % n_spatial, s % n_spatial
		p_spin, q_spin = p // n_spatial, q // n_spatial
		r_spin, s_spin = r // n_spatial, s // n_spatial
		value = 0.0
		if p_spin == r_spin and q_spin == s_spin:
			value += eri_mo[p_spatial, r_spatial, q_spatial, s_spatial]
		if p_spin == s_spin and q_spin == r_spin:
			value -= eri_mo[p_spatial, s_spatial, q_spatial, r_spatial]
		return value

	def _ao_integrals(self, molecule, basis):
		import veloxchem as vlx
		from .fockdriver import FockDriver
		from .oneeints import compute_nuclear_potential_integrals
		from .veloxchemlib import EcpDriver

		T_ao = vlx.KineticEnergyDriver().compute(molecule, basis).to_numpy()
		if basis.has_ecp():
			charges = molecule.get_effective_nuclear_charges(basis)
			coordinates = molecule.get_coordinates_in_bohr()
			V_ao = compute_nuclear_potential_integrals(
				molecule,
				basis,
				charges,
				coordinates,
			)
		else:
			V_ao = compute_nuclear_potential_integrals(molecule, basis)
		H_ao = T_ao + V_ao
		if basis.has_ecp():
			core_electrons = basis.get_number_of_ecp_core_electrons()
			ecp_atom_indices = [
				index for index, nelectrons in enumerate(core_electrons)
				if nelectrons > 0
			]
			if ecp_atom_indices:
				H_ao += EcpDriver().compute(
					molecule,
					basis,
					ecp_atom_indices,
				).to_numpy()
		S_ao = vlx.OverlapDriver().compute(molecule, basis).to_numpy()
		eri_ao = FockDriver().compute_eri(molecule, basis, eri_thresh=1.0e-12)
		if np.max(np.abs(eri_ao)) <= 1.0e-14 and molecule.number_of_electrons() > 1:
			raise RuntimeError(
				"AO electron-repulsion integrals are zero; cannot build a physical VB Hamiltonian.")
		e_nuc = molecule.effective_nuclear_repulsion_energy(basis)
		return H_ao, S_ao, eri_ao, e_nuc

	def _frozen_hf_embedding(self,
							 h_ao,
							 s_ao,
							 eri_ao,
							 e_nuc,
							 active_space,
							 analysis):
		nao_data = getattr(analysis, 'nao_data', None)
		total_density = getattr(analysis, 'density', None)
		if nao_data is None or total_density is None:
			return None
		try:
			total_electrons = float(np.einsum(
				'uv,vu->',
				np.array(total_density, dtype=float),
				s_ao,
				optimize=True,
			))
			active_target_electrons = float(
				active_space.metadata.get('electron_count') or 0.0)
			if active_target_electrons >= total_electrons - 1.0e-8:
				return None
		except Exception:
			pass
		active_vector = self._candidate_nao_vector(
			active_space.active_candidate,
			nao_data,
		)
		if active_vector is not None:
			active_density_nao = 2.0 * np.outer(active_vector, active_vector)
			active_density = nao_data.transform @ active_density_nao @ nao_data.transform.T
		else:
			active_density = self._active_subspace_reference_density(
				active_space,
				total_density,
				s_ao,
			)
			if active_density is None:
				return None
		active_density = 0.5 * (active_density + active_density.T)
		total_density = 0.5 * (np.array(total_density, dtype=float) +
							 np.array(total_density, dtype=float).T)
		frozen_density = total_density - active_density
		frozen_density = 0.5 * (frozen_density + frozen_density.T)

		j_frozen = np.einsum('lm,uvlm->uv', frozen_density, eri_ao, optimize=True)
		k_frozen = np.einsum('lm,ulvm->uv', frozen_density, eri_ao, optimize=True)
		g_frozen = j_frozen - 0.5 * k_frozen
		h_effective = h_ao + g_frozen
		constant_energy = (
			float(np.einsum('uv,uv->', frozen_density, h_ao, optimize=True)) +
			0.5 * float(np.einsum('uv,uv->', frozen_density, g_frozen, optimize=True)) +
			float(e_nuc)
		)
		j_active = np.einsum('lm,uvlm->uv', active_density, eri_ao, optimize=True)
		k_active = np.einsum('lm,ulvm->uv', active_density, eri_ao, optimize=True)
		g_active = j_active - 0.5 * k_active
		active_reference_energy = (
			constant_energy +
			float(np.einsum('uv,uv->', active_density, h_effective, optimize=True)) +
			0.5 * float(np.einsum('uv,uv->', active_density, g_active, optimize=True))
		)
		j_total = np.einsum('lm,uvlm->uv', total_density, eri_ao, optimize=True)
		k_total = np.einsum('lm,ulvm->uv', total_density, eri_ao, optimize=True)
		g_total = j_total - 0.5 * k_total
		reference_total_energy = (
			float(np.einsum('uv,uv->', total_density, h_ao, optimize=True)) +
			0.5 * float(np.einsum('uv,uv->', total_density, g_total, optimize=True)) +
			float(e_nuc)
		)
		active_electrons = float(np.einsum('uv,vu->', active_density, s_ao, optimize=True))
		frozen_electrons = float(np.einsum('uv,vu->', frozen_density, s_ao, optimize=True))

		return {
			"h_effective": 0.5 * (h_effective + h_effective.T),
			"constant_energy": constant_energy,
			"active_reference_energy": active_reference_energy,
			"reference_total_energy": reference_total_energy,
			"active_electron_count": active_electrons,
			"frozen_electron_count": frozen_electrons,
			"frozen_density": frozen_density,
		}

	def _project_active_orbitals_out_of_frozen_density(self, orbitals, overlap, frozen_density):
		if frozen_density is None:
			return orbitals
		projector = 0.5 * np.array(frozen_density, dtype=float) @ overlap
		projected_orbitals = []
		for orbital in orbitals:
			coefficients = np.array(orbital.coefficients, dtype=float)
			projected = coefficients - projector @ coefficients
			try:
				projected = self._s_normalize(projected, overlap)
			except RuntimeError:
				projected = coefficients
			projected_orbitals.append(
				VbOrbital(
					label=orbital.label,
					coefficients=projected,
					center=orbital.center,
					kind=orbital.kind,
				)
			)
		return projected_orbitals

	def _active_subspace_reference_density(self, active_space, total_density, s_ao):
		try:
			coefficients = np.column_stack([
				orbital.coefficients for orbital in active_space.active_orbitals
			])
		except Exception:
			return None
		if coefficients.size == 0:
			return None
		S_active = coefficients.T @ s_ao @ coefficients
		try:
			values, vectors = np.linalg.eigh(0.5 * (S_active + S_active.T))
			keep = values > 1.0e-10
			if not np.any(keep):
				return None
			orthogonalizer = vectors[:, keep] / np.sqrt(values[keep])
			Q = coefficients @ orthogonalizer
			D_projected = Q.T @ s_ao @ np.array(total_density, dtype=float) @ s_ao @ Q
			D_projected = 0.5 * (D_projected + D_projected.T)
			active_density = Q @ D_projected @ Q.T
		except np.linalg.LinAlgError:
			return None
		active_electrons = float(np.einsum('uv,vu->', active_density, s_ao, optimize=True))
		target_electrons = float(active_space.metadata.get('electron_count') or active_electrons)
		if active_electrons > 1.0e-10 and target_electrons > 0.0:
			active_density *= target_electrons / active_electrons
		return 0.5 * (active_density + active_density.T)

	def _s_normalize(self, vector, overlap):
		norm2 = float(vector.T @ overlap @ vector)
		if norm2 <= 1.0e-14:
			raise RuntimeError("Cannot normalize orbital with near-zero overlap norm.")
		return vector / np.sqrt(norm2)

	def _center_breathing_vector(self,
							 molecule,
							 basis,
							 center,
							 active_vector,
							 occupied_vectors,
							 overlap):
		if center is None:
			return None
		try:
			center = int(center)
			ao_to_atom, ao_to_l = OrbitalAnalyzer.ao_shell_map(molecule, basis)
			ao_to_angular = self._ao_angular_component_map(molecule, basis)
		except Exception:
			return None
		center_indices = np.where(ao_to_atom == center)[0]
		if len(center_indices) <= 1:
			return None

		active_weights = {}
		for index in center_indices:
			l_value = int(ao_to_l[index])
			active_weights[l_value] = active_weights.get(l_value, 0.0) + float(active_vector[index]**2)
		if not active_weights:
			return None
		dominant_l = max(active_weights, key=active_weights.get)
		angular_weights = {}
		for index in center_indices:
			if int(ao_to_l[index]) != int(dominant_l):
				continue
			angular_label = ao_to_angular[index]
			angular_weights[angular_label] = angular_weights.get(angular_label, 0.0) + float(active_vector[index]**2)
		dominant_angular = max(angular_weights, key=angular_weights.get) if angular_weights else None
		candidate_indices = [
			int(index) for index in center_indices
			if int(ao_to_l[index]) == int(dominant_l) and
			(dominant_angular is None or ao_to_angular[index] == dominant_angular)
		]
		if len(candidate_indices) <= 1:
			return None

		projectors = []
		for vector in occupied_vectors:
			try:
				projectors.append(self._s_normalize(np.array(vector, dtype=float), overlap))
			except RuntimeError:
				continue

		best_vector = None
		best_norm = 0.0
		for index in candidate_indices:
			trial = np.zeros_like(active_vector, dtype=float)
			trial[index] = 1.0
			for projector in projectors:
				trial = trial - projector * float(projector.T @ overlap @ trial)
			trial_norm = float(trial.T @ overlap @ trial)
			if trial_norm > best_norm:
				best_norm = trial_norm
				best_vector = trial

		if best_vector is None or best_norm <= 1.0e-12:
			return None
		try:
			return self._s_normalize(best_vector, overlap)
		except RuntimeError:
			return None

	def _ao_angular_component_map(self, molecule, basis):
		components = []
		for entry in basis.get_ao_basis_map(molecule):
			tokens = entry.split()
			shell_label = tokens[2] if len(tokens) >= 3 else "s"
			component = shell_label.lstrip("0123456789") or shell_label
			components.append(component)
		return np.array(components, dtype=object)

	def _structure_product_terms(self, structure):
		occupation = tuple(x for occ in structure.occupation for x in occ)
		if sum(occupation) != 2:
			raise RuntimeError(
				"Two-electron singlet VB algebra supports only structures with two active electrons.")
		paired = [index for index, occ in enumerate(occupation) if occ == 2]
		singly = [index for index, occ in enumerate(occupation) if occ == 1]
		if len(paired) == 1 and not singly and all(occ in (0, 2) for occ in occupation):
			orbital = paired[0]
			return ((orbital, orbital, 1.0),)
		if len(singly) == 2 and not paired and all(occ in (0, 1) for occ in occupation):
			left, right = singly
			return ((left, right, 1.0 / np.sqrt(2.0)),
					(right, left, 1.0 / np.sqrt(2.0)))
		raise RuntimeError(
			"Two-electron singlet VB algebra supports only covalent and ionic occupation patterns.")

	def _build_two_orbital_singlet_matrices(self,
										 structures,
										 coefficients,
										 h_ao,
										 s_ao,
										 eri_ao,
										 e_nuc):
		S_mo = coefficients.T @ s_ao @ coefficients
		H_mo = coefficients.T @ h_ao @ coefficients
		eri_mo = np.einsum('up,vq,lr,ms,uvlm->pqrs',
							 coefficients,
							 coefficients,
							 coefficients,
							 coefficients,
							 eri_ao)
		terms = [self._structure_product_terms(structure)
				 for structure in structures]
		n_struct = len(structures)
		S = np.zeros((n_struct, n_struct))
		H = np.zeros((n_struct, n_struct))
		for i, left_terms in enumerate(terms):
			for j, right_terms in enumerate(terms):
				for p, q, left_coeff in left_terms:
					for r, s, right_coeff in right_terms:
						overlap = S_mo[p, r] * S_mo[q, s]
						hamiltonian = (
							H_mo[p, r] * S_mo[q, s] +
							S_mo[p, r] * H_mo[q, s] +
							eri_mo[p, r, q, s] +
							e_nuc * overlap
						)
						factor = left_coeff * right_coeff
						S[i, j] += factor * overlap
						H[i, j] += factor * hamiltonian
		return 0.5 * (S + S.T), 0.5 * (H + H.T)

	def _build_structure_specific_two_orbital_singlet_matrices(self,
													 structures,
													 coefficient_sets,
													 h_ao,
													 s_ao,
													 eri_ao,
													 e_nuc):
		terms = [self._structure_product_terms(structure)
				 for structure in structures]
		n_struct = len(structures)
		S = np.zeros((n_struct, n_struct))
		H = np.zeros((n_struct, n_struct))
		for i, left_terms in enumerate(terms):
			left_coefficients = coefficient_sets[i]
			for j, right_terms in enumerate(terms):
				right_coefficients = coefficient_sets[j]
				S_cross = left_coefficients.T @ s_ao @ right_coefficients
				H_cross = left_coefficients.T @ h_ao @ right_coefficients
				for p, q, left_coeff in left_terms:
					for r, s, right_coeff in right_terms:
						overlap = S_cross[p, r] * S_cross[q, s]
						eri = np.einsum(
							'u,v,l,m,uvlm->',
							left_coefficients[:, p],
							right_coefficients[:, r],
							left_coefficients[:, q],
							right_coefficients[:, s],
							eri_ao,
							optimize=True,
						)
						hamiltonian = (
							H_cross[p, r] * S_cross[q, s] +
							S_cross[p, r] * H_cross[q, s] +
							eri +
							e_nuc * overlap
						)
						factor = left_coeff * right_coeff
						S[i, j] += factor * overlap
						H[i, j] += factor * hamiltonian
		return 0.5 * (S + S.T), 0.5 * (H + H.T)

	def _solve_generalized_vb(self,
							 hamiltonian,
							 overlap,
							 threshold=1.0e-10,
							 return_lowdin=False):
		overlap_values, overlap_vectors = np.linalg.eigh(0.5 * (overlap + overlap.T))
		keep = np.where(overlap_values > threshold)[0]
		if len(keep) == 0:
			raise RuntimeError("VB structure overlap matrix has no positive-norm subspace.")
		X = overlap_vectors[:, keep] / np.sqrt(overlap_values[keep])
		H_orth = X.T @ (0.5 * (hamiltonian + hamiltonian.T)) @ X
		eigvals, eigvecs = np.linalg.eigh(0.5 * (H_orth + H_orth.T))
		root = int(np.argmin(eigvals))
		coeffs = X @ eigvecs[:, root]
		norm = float(coeffs.T @ overlap @ coeffs)
		if norm <= 0.0:
			raise RuntimeError("VB root has non-positive overlap norm.")
		coeffs = coeffs / np.sqrt(norm)
		weights = coeffs * (overlap @ coeffs)
		weight_sum = float(np.sum(weights))
		if abs(weight_sum) > 1.0e-14:
			weights = weights / weight_sum
		lowdin_coeffs = (
			(overlap_vectors[:, keep] * np.sqrt(overlap_values[keep])) @
			(overlap_vectors[:, keep].T @ coeffs)
		)
		lowdin_weights = lowdin_coeffs**2
		lowdin_weight_sum = float(np.sum(lowdin_weights))
		if abs(lowdin_weight_sum) > 1.0e-14:
			lowdin_weights = lowdin_weights / lowdin_weight_sum
		if return_lowdin:
			return float(eigvals[root]), coeffs, weights, keep, lowdin_weights
		return float(eigvals[root]), coeffs, weights, keep

	def compute_vbci_h2(self, molecule, basis, structures, orbitals, reference_orbitals, options):
		return self._compute_two_orbital_vbci(molecule, basis, structures, orbitals)
