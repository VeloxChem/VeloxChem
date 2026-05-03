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
	mode: str - calculation mode (e.g. 'vbscf')
	max_iter: int - maximum optimization iterations
	conv_thresh: float - convergence threshold
	optimize_orbitals: bool - whether to optimize orbitals
	include_ionic: bool - include ionic structures
	include_bovb: bool - enable BOVB (breathing orbital VB)
	active_bond: optional atom-index pair for one-active-bond generation
	active_candidate_label: optional OrbitalAnalyzer candidate label to activate
	active_candidate_subtype: optional bond subtype filter, currently 'sigma' or 'pi'
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
			analysis = OrbitalAnalyzer(molecule, basis).classify_nao_data(nao_data)
			candidates = analysis.orbital_candidates

		# Build VbOrbital objects for each classified orbital
		vb_orbitals = [
			VbOrbital(
				label=f"{c['type']}_{c['subtype']}_{c['index']}",
				coefficients=np.array([coef['coefficient'] for coef in sorted(c['coefficients'], key=lambda x: x['nao_index'])]),
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
			if options.mode == 'vbscf' and options.optimize_orbitals:
				if len(active_space.active_orbitals) != 2:
					raise RuntimeError("VB-SCF orbital optimization is currently only available for two active orbitals.")
				active_result = self.compute_vbscf_h2(
					molecule,
					basis,
					list(active_space.structures),
					list(active_space.active_orbitals),
					reference_orbitals,
					options,
				)
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
			active_result["diagnostics"].update(
				self._active_space_diagnostics(active_space))
			return active_result

		if structures is None:
			structures = [VbStructure(label='default', occupation=(tuple([1]*len(orbitals)),), spin='singlet')]

		if self._is_minimal_h2_vbci_case(molecule, structures, orbitals):
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

	def compute_vbscf_h2(self, molecule, basis, structures, orbitals, reference_orbitals, options):
		"""
		VB-SCF for H2: optimize two active orbitals by a mixing angle.
		"""
		from scipy.optimize import minimize_scalar

		H_ao, S_ao, eri_ao, e_nuc = self._ao_integrals(molecule, basis)

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

	def _candidate_summary_record(self, candidate):
		return {
			"label": self._candidate_label(candidate),
			"type": candidate.get("type"),
			"subtype": candidate.get("subtype"),
			"atoms": tuple(candidate.get("atoms", ())),
			"occupation": candidate.get("occupation"),
			"source": candidate.get("source"),
		}

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
		if getattr(options, 'active_pi_atoms', None) is not None:
			return True
		return self._is_h2_molecule(molecule)

	def _build_one_active_bond_space(self, molecule, basis, candidates, nao_data, options):
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

		active_candidate = self._select_active_bond_candidate(candidates, options)
		if active_candidate is None:
			if not self._is_h2_molecule(molecule):
				raise RuntimeError("VB active-space builder could not find a two-center bond candidate.")
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

		active_orbitals = tuple(
			self._candidate_centered_active_orbitals(
				molecule,
				basis,
				active_candidate,
				nao_data,
			) or self._atom_centered_active_orbitals(molecule, basis, active_bond))
		structures = tuple(self._default_one_bond_structures(active_bond, options))
		active_label = self._candidate_label(active_candidate)
		active_index = int(active_candidate.get('index', -1))
		active_subtype = active_candidate.get('subtype')

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
			},
		)

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
			if candidate.get('type') != 'BD':
				continue
			if candidate.get('subtype') not in ('sigma', 'pi'):
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

	def _candidate_nao_vector(self, candidate, nao_data):
		if nao_data is None or not candidate.get('coefficients'):
			return None
		vector = np.zeros(len(nao_data.populations))
		for item in candidate.get('coefficients', []):
			index = int(item['nao_index']) - 1
			if 0 <= index < len(vector):
				vector[index] = float(item['coefficient'])
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
			"active_electron_count": active_space.metadata.get("electron_count"),
			"active_spin": active_space.metadata.get("spin"),
			"active_candidate_label": active_space.active_candidate_label,
			"active_candidate_type": active_space.active_candidate.get("type"),
			"active_candidate_subtype": active_space.active_candidate.get("subtype"),
			"component_pi_candidate_labels": active_space.metadata.get("component_pi_candidate_labels", []),
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
			diagnostics.update({
				"frozen_hf_embedding": True,
				"frozen_electron_count": embedding["frozen_electron_count"],
				"active_reference_electron_count": embedding["active_electron_count"],
				"frozen_constant_energy": embedding["constant_energy"],
				"embedding_model": "inactive HF density = total RHF density minus active reference density",
			})
		diagnostics.update(root_info)
		return {
			"energy": float(energy),
			"structure_coefficients": coeffs,
			"overlap": S,
			"Hamiltonian": H,
			"weights": weights,
			"lowdin_weights": lowdin_weights,
			"orbitals": [orb.coefficients for orb in orbitals],
			"orthonormal_active_orbitals": orth_coefficients,
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

		T_ao = vlx.KineticEnergyDriver().compute(molecule, basis).to_numpy()
		V_ao = vlx.NuclearPotentialDriver().compute(molecule, basis).to_numpy()
		H_ao = T_ao - V_ao
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
		frozen_density = np.array(total_density, dtype=float) - active_density
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
		active_electrons = float(np.einsum('uv,vu->', active_density, s_ao, optimize=True))
		frozen_electrons = float(np.einsum('uv,vu->', frozen_density, s_ao, optimize=True))

		return {
			"h_effective": 0.5 * (h_effective + h_effective.T),
			"constant_energy": constant_energy,
			"active_electron_count": active_electrons,
			"frozen_electron_count": frozen_electrons,
		}

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
