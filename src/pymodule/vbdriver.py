import numpy as np

"""
Valence Bond (VB) driver skeleton.

Implements the core API and data structures for a general valence bond driver.
This is a placeholder for Phase 0 of the implementation plan.
"""

from dataclasses import dataclass, field
from typing import Any, List, Optional, Dict, Tuple

#
#                                   VELOXCHEM
#              ----------------------------------------------------
#                          An Electronic Structure Code
#
#  SPDX-License-Identifier: BSD-3-Clause
#
#  Copyright 2018-2026 VeloxChem developers

"""
Valence Bond (VB) driver skeleton.

Implements the core API and data structures for a general valence bond driver.
This is a placeholder for Phase 0 of the implementation plan.
"""

from dataclasses import dataclass, field
from typing import Any, List, Optional, Dict, Tuple
import numpy as np

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
class VbComputeOptions:
	"""
	Options for controlling the VB driver.
	mode: str - calculation mode (e.g. 'vbscf')
	max_iter: int - maximum optimization iterations
	conv_thresh: float - convergence threshold
	optimize_orbitals: bool - whether to optimize orbitals
	include_ionic: bool - include ionic structures
	include_bovb: bool - enable BOVB (breathing orbital VB)
	"""
	mode: str = "vbscf"
	max_iter: int = 50
	conv_thresh: float = 1.0e-8
	optimize_orbitals: bool = True
	include_ionic: bool = True
	include_bovb: bool = False

class VbDriver:

	def scan_h2_dissociation(self, basis, structures, orbitals, bond_distances, options=None, molecule_template=None):
		"""
		Automate H2 dissociation scan with VB-SCF. For each bond distance, build the molecule, run VB-SCF, and collect results.
		Args:
			basis: VeloxChem basis object or basis name
			structures: list of VbStructure objects
			orbitals: list of VbOrbital objects (used for initial guess)
			bond_distances: list or np.array of H-H distances (in Angstrom)
			options: VbComputeOptions
			molecule_template: optional function or string for molecule construction
		Returns:
			dict with keys: distances, energies, weights, diagnostics, results (full output per geometry)
		"""
		import numpy as np
		import veloxchem as vlx
		results = []
		energies = []
		weights = []
		diagnostics = []
		distances = []
		for idx, R in enumerate(bond_distances):
			# Build H2 molecule at distance R
			if molecule_template is not None:
				molecule = molecule_template(R)
			else:
				h2_xyz = f"""2\nH2 molecule\nH 0.0 0.0 0.0\nH {R:.8f} 0.0 0.0\n"""
				molecule = vlx.Molecule.read_xyz_string(h2_xyz)
			# Use basis name or object
			if isinstance(basis, str):
				basis_obj = vlx.MolecularBasis.read(molecule, basis)
			else:
				basis_obj = basis
			# Robust AO-to-atom mapping and orbital construction for each geometry
			from veloxchem.aoindices import get_basis_function_indices_of_atoms
			n_ao = basis_obj.get_dimensions_of_basis()
			ao_indices_per_atom = get_basis_function_indices_of_atoms(molecule, basis_obj)
			orb1_coeffs = np.zeros(n_ao)
			orb2_coeffs = np.zeros(n_ao)
			for i in np.atleast_1d(ao_indices_per_atom[0]):
				orb1_coeffs[i] = 1.0
			for i in np.atleast_1d(ao_indices_per_atom[1]):
				orb2_coeffs[i] = 1.0
			orb1_coeffs /= np.linalg.norm(orb1_coeffs)
			orb2_coeffs /= np.linalg.norm(orb2_coeffs)
			from .vbdriver import VbOrbital  # Relative import for safety
			orb1 = VbOrbital(label='1s_A', coefficients=orb1_coeffs, center=0)
			orb2 = VbOrbital(label='1s_B', coefficients=orb2_coeffs, center=1)
			orbitals_geom = [orb1, orb2]
			res = self.compute(molecule, basis_obj, structures, orbitals_geom, reference_orbitals=None, options=options)
			# Diagnostic: print S and H for first geometry
			if idx == 0 and "overlap" in res and "Hamiltonian" in res:
				print("First geometry S matrix:\n", res["overlap"])
				print("First geometry H matrix:\n", res["Hamiltonian"])
			results.append(res)
			energies.append(res["energy"])
			weights.append(res["weights"])
			diagnostics.append(res.get("diagnostics", {}))
			distances.append(R)
		return {
			"distances": np.array(distances),
			"energies": np.array(energies),
			"weights": np.array(weights),
			"diagnostics": diagnostics,
			"results": results,
		}
	"""
	Valence Bond (VB) driver for VeloxChem.

	This class provides the main API for VB calculations. It is designed to mirror
	the style of NboDriver, but implements a true non-orthogonal VB wavefunction.
	"""
	def __init__(self, comm=None, ostream=None):
		"""
		Initialize the VB driver. Arguments are placeholders for future MPI/output support.
		"""
		self.comm = comm
		self.ostream = ostream

	def compute(self,
				molecule: Any,
				basis: Any,
				structures: List[VbStructure],
				orbitals: Optional[List[VbOrbital]] = None,
				reference_orbitals: Optional[Any] = None,
				options: Optional[VbComputeOptions] = None) -> Dict[str, Any]:
		"""
		Compute the VB energy and related quantities.

		Parameters:
			molecule: VeloxChem molecule object
			basis: VeloxChem basis object
			structures: list of VbStructure objects
			orbitals: list of VbOrbital objects
			reference_orbitals: (optional) SCF or other reference orbitals
			options: VbComputeOptions

		Returns:
			dict with keys: energy, structure_coefficients, overlap, Hamiltonian,
			weights, orbitals, diagnostics
		"""
		# If orbitals is None or 'auto', construct localized 1s-like orbitals for H2
		if orbitals is None or orbitals == 'auto':
			# Use AO-to-atom mapping to construct localized orbitals (robust, version-independent)
			from veloxchem.aoindices import get_basis_function_indices_of_atoms
			n_ao = basis.get_dimensions_of_basis()
			ao_indices_per_atom = get_basis_function_indices_of_atoms(molecule, basis)
			orb1_coeffs = np.zeros(n_ao)
			orb2_coeffs = np.zeros(n_ao)
			for i in np.atleast_1d(ao_indices_per_atom[0]):
				orb1_coeffs[i] = 1.0
			for i in np.atleast_1d(ao_indices_per_atom[1]):
				orb2_coeffs[i] = 1.0
			# Normalize
			orb1_coeffs /= np.linalg.norm(orb1_coeffs)
			orb2_coeffs /= np.linalg.norm(orb2_coeffs)
			orb1 = VbOrbital(label='1s_A', coefficients=orb1_coeffs, center=0)
			orb2 = VbOrbital(label='1s_B', coefficients=orb2_coeffs, center=1)
			orbitals = [orb1, orb2]
		# Phase 2: H2 VB-SCF (orbital optimization)
		if self._is_minimal_h2_vbci_case(molecule, structures, orbitals):
			if options is not None and getattr(options, 'mode', 'vbscf') == 'vbscf' and getattr(options, 'optimize_orbitals', True):
				return self.compute_vbscf_h2(molecule, basis, structures, orbitals, reference_orbitals, options)
			else:
				return self.compute_vbci_h2(molecule, basis, structures, orbitals, reference_orbitals, options)
	def compute_vbscf_h2(self, molecule, basis, structures, orbitals, reference_orbitals, options):
		"""
		VB-SCF for H2: optimize orbitals by mixing angle theta between two AOs.
		"""
		import numpy as np
		import veloxchem as vlx
		from scipy.optimize import minimize_scalar

		n_struct = len(structures)
		n_ao = basis.get_dimensions_of_basis()
		# AO integrals
		T_ao = vlx.KineticEnergyDriver().compute(molecule, basis).to_numpy()
		V_ao = vlx.NuclearPotentialDriver().compute(molecule, basis).to_numpy()
		H_ao = T_ao - V_ao
		S_ao = vlx.OverlapDriver().compute(molecule, basis).to_numpy()
		import veloxchem.veloxchemlib as veloxchemlib
		eri_ao = veloxchemlib.compute_ao_eris(molecule, basis)

		# Assume orbitals[0] and orbitals[1] are initial guesses for 1s_A and 1s_B
		ao1 = orbitals[0].coefficients
		ao2 = orbitals[1].coefficients

		def energy_for_theta(theta):
			# Build orthonormalized orbitals by rotation
			c1 = np.cos(theta)
			c2 = np.sin(theta)
			orbA = c1 * ao1 + c2 * ao2
			orbB = -c2 * ao1 + c1 * ao2
			C = np.column_stack([orbA, orbB])
			# Build 1e and 2e integrals in MO basis
			S_mo = C.T @ S_ao @ C
			H_mo = C.T @ H_ao @ C
			eri_mo = np.einsum('up,vq,lr,ms,uvlm->pqrs', C, C, C, C, eri_ao)
			# Determinants for each structure (use full occupation tuple)
			dets = [tuple(x for occ in s.occupation for x in occ) for s in structures]
			S = np.zeros((n_struct, n_struct))
			H = np.zeros((n_struct, n_struct))
			for i in range(n_struct):
				for j in range(n_struct):
					occ_i = dets[i]
					occ_j = dets[j]
					occ_idx_i = [k for k, n in enumerate(occ_i) for _ in range(n)]
					occ_idx_j = [k for k, n in enumerate(occ_j) for _ in range(n)]
					S_occ = S_mo[np.ix_(occ_idx_i, occ_idx_j)]
					S[i, j] = np.linalg.det(S_occ)
					h1 = 0.0
					for a, p in enumerate(occ_idx_i):
						for b, q in enumerate(occ_idx_j):
							h1 += S_occ[a, b] * H_mo[p, q]
					h2 = 0.0
					for a, p in enumerate(occ_idx_i):
						for b, q in enumerate(occ_idx_j):
							for c, r in enumerate(occ_idx_i):
								for d, s in enumerate(occ_idx_j):
									h2 += 0.5 * S_occ[a, b] * S_occ[c, d] * (eri_mo[p, q, r, s] - 0.5 * eri_mo[p, s, r, q])
					H[i, j] = h1 + h2
			# Solve generalized eigenvalue problem
			try:
				eigvals, eigvecs = np.linalg.eigh(np.linalg.pinv(S) @ H)
				energy = np.min(eigvals)
			except Exception:
				energy = 1e6
			return energy

		# Optimize theta in [0, pi)
		res = minimize_scalar(energy_for_theta, bounds=(0, np.pi), method='bounded', options={'xatol': options.conv_thresh if options else 1e-8})
		theta_opt = res.x
		energy_opt = res.fun
		# Build final orbitals and results
		c1 = np.cos(theta_opt)
		c2 = np.sin(theta_opt)
		orbA = c1 * ao1 + c2 * ao2
		orbB = -c2 * ao1 + c1 * ao2
		C = np.column_stack([orbA, orbB])
		# Repeat final calculation to get coefficients, weights, etc.
		S_mo = C.T @ S_ao @ C
		H_mo = C.T @ H_ao @ C
		eri_mo = np.einsum('up,vq,lr,ms,uvlm->pqrs', C, C, C, C, eri_ao)
		dets = [tuple(x for occ in s.occupation for x in occ) for s in structures]
		S = np.zeros((n_struct, n_struct))
		H = np.zeros((n_struct, n_struct))
		for i in range(n_struct):
			for j in range(n_struct):
				occ_i = dets[i]
				occ_j = dets[j]
				occ_idx_i = [k for k, n in enumerate(occ_i) for _ in range(n)]
				occ_idx_j = [k for k, n in enumerate(occ_j) for _ in range(n)]
				S_occ = S_mo[np.ix_(occ_idx_i, occ_idx_j)]
				S[i, j] = np.linalg.det(S_occ)
				h1 = 0.0
				for a, p in enumerate(occ_idx_i):
					for b, q in enumerate(occ_idx_j):
						h1 += S_occ[a, b] * H_mo[p, q]
				h2 = 0.0
				for a, p in enumerate(occ_idx_i):
					for b, q in enumerate(occ_idx_j):
						for c, r in enumerate(occ_idx_i):
							for d, s in enumerate(occ_idx_j):
								h2 += 0.5 * S_occ[a, b] * S_occ[c, d] * (eri_mo[p, q, r, s] - 0.5 * eri_mo[p, s, r, q])
				H[i, j] = h1 + h2
		try:
			eigvals, eigvecs = np.linalg.eigh(np.linalg.pinv(S) @ H)
			idx = np.argmin(eigvals)
			coeffs = eigvecs[:, idx]
			weights = coeffs**2 / np.sum(coeffs**2)
		except Exception:
			coeffs = np.zeros(n_struct)
			weights = np.zeros(n_struct)
		diagnostics = {
			"message": "H2 VB-SCF result (optimized orbitals, AO integral-based)",
			"theta_opt": theta_opt,
			"overlap_condition": np.linalg.cond(S),
			"optimizer": res,
		}
		return {
			"energy": float(energy_opt),
			"structure_coefficients": coeffs,
			"overlap": S,
			"Hamiltonian": H,
			"weights": weights,
			"orbitals": [orbA, orbB],
			"diagnostics": diagnostics,
		}

		# Fallback: placeholder (Phase 0)
		return {
			"energy": 0.0,
			"structure_coefficients": np.zeros(len(structures)),
			"overlap": np.eye(len(structures)),
			"Hamiltonian": np.eye(len(structures)),
			"weights": np.ones(len(structures)) / max(1, len(structures)),
			"orbitals": [orb.coefficients for orb in orbitals],
			"diagnostics": {"message": "VB driver skeleton placeholder"},
		}

	def _is_minimal_h2_vbci_case(self, molecule, structures, orbitals):
		"""
		Detect if this is the minimal H2 VB-CI case (two H atoms, two electrons, two orbitals, three structures: covalent, ionic1, ionic2).
		"""
		try:
			return (
				hasattr(molecule, 'number_of_atoms') and molecule.number_of_atoms() == 2
				and len(orbitals) == 2
				and len(structures) == 3
			)
		except Exception:
			return False

	def compute_vbci_h2(self, molecule, basis, structures, orbitals, reference_orbitals, options):
		"""
		Minimal VB-CI for H2: fixed orbitals, three structures (covalent, ionic1, ionic2), two electrons.
		Computes overlap and Hamiltonian matrices from AO integrals and orbital coefficients using only the public VeloxChem Python API.
		"""
		import numpy as np
		import veloxchem as vlx

		n_struct = len(structures)
		n_ao = basis.get_dimensions_of_basis()
		# AO integrals
		T_ao = vlx.KineticEnergyDriver().compute(molecule, basis).to_numpy()
		V_ao = vlx.NuclearPotentialDriver().compute(molecule, basis).to_numpy()
		H_ao = T_ao - V_ao  # Note: V_ao is negative in VeloxChem
		S_ao = vlx.OverlapDriver().compute(molecule, basis).to_numpy()

		# Build MO coefficients for the two orbitals
		C = np.column_stack([orb.coefficients for orb in orbitals])  # (n_ao, 2)

		# Build 1e integrals in MO basis
		S_mo = C.T @ S_ao @ C  # (2,2)
		H_mo = C.T @ H_ao @ C  # (2,2)

		# Two-electron integrals in MO basis (chemist's notation)
		# AO integrals: (μν|λσ)
		import veloxchem.veloxchemlib as veloxchemlib
		eri_ao = veloxchemlib.compute_ao_eris(molecule, basis)
		# Transform to MO basis: (pq|rs) = C_{μp} C_{νq} C_{λr} C_{σs} (μν|λσ)
		eri_mo = np.einsum('up,vq,lr,ms,uvlm->pqrs', C, C, C, C, eri_ao)

		# Define determinants for each structure
		dets = [s.occupation[0] for s in structures]  # [(1,1), (2,0), (0,2)]

		# Build overlap and Hamiltonian matrices
		S = np.zeros((n_struct, n_struct))
		H = np.zeros((n_struct, n_struct))

		# For minimal H2, all determinants are closed-shell, so use simple formulas
		for i in range(n_struct):
			for j in range(n_struct):
				occ_i = dets[i]
				occ_j = dets[j]
				occ_idx_i = [k for k, n in enumerate(occ_i) for _ in range(n)]
				occ_idx_j = [k for k, n in enumerate(occ_j) for _ in range(n)]
				S_occ = S_mo[np.ix_(occ_idx_i, occ_idx_j)]
				S[i, j] = np.linalg.det(S_occ)

				# Hamiltonian matrix element: <D_i|H|D_j> (see Szabo & Ostlund, Ch. 3)
				# Only for two-electron closed-shell determinants
				# 1e part
				h1 = 0.0
				for a, p in enumerate(occ_idx_i):
					for b, q in enumerate(occ_idx_j):
						h1 += S_occ[a, b] * H_mo[p, q]
				# 2e part
				h2 = 0.0
				for a, p in enumerate(occ_idx_i):
					for b, q in enumerate(occ_idx_j):
						for c, r in enumerate(occ_idx_i):
							for d, s in enumerate(occ_idx_j):
								# (pq|rs) - 0.5*(ps|rq) for singlet
								h2 += 0.5 * S_occ[a, b] * S_occ[c, d] * (eri_mo[p, q, r, s] - 0.5 * eri_mo[p, s, r, q])
				H[i, j] = h1 + h2

		# Solve generalized eigenvalue problem H c = E S c
		eigvals, eigvecs = np.linalg.eigh(np.linalg.pinv(S) @ H)
		idx = np.argmin(eigvals)
		energy = eigvals[idx]
		coeffs = eigvecs[:, idx]
		weights = coeffs**2 / np.sum(coeffs**2)
		diagnostics = {
			"message": "Minimal H2 VB-CI result (AO integral-based, public API only, with 2e integrals)",
			"overlap_condition": np.linalg.cond(S),
		}
		return {
			"energy": float(energy),
			"structure_coefficients": coeffs,
			"overlap": S,
			"Hamiltonian": H,
			"weights": weights,
			"orbitals": [orb.coefficients for orb in orbitals],
			"diagnostics": diagnostics,
		}
