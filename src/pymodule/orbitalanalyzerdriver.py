
"""Shared orbital analysis payload for NBO, VB, and related drivers."""

from dataclasses import dataclass, field
import importlib
import numpy as np

from .oneeints import compute_overlap_integrals


@dataclass(frozen=True)
class OrbitalAnalyzerOptions:
    """Options for :class:`OrbitalAnalyzer`."""

    threshold: float = 1.0e-12
    include_mo_analysis: bool = True
    include_nbo_candidates: bool = True
    mo_analysis_top: int = 6
    mo_analysis_threshold: float = 1.0e-2
    lone_pair_min_occupation: float = 1.50
    rydberg_max_occupation: float = 0.50
    bond_min_occupation: float = 1.20
    bond_min_atom_weight: float = 0.10
    pi_min_occupation: float = 0.20
    conjugated_pi_max_path: int = 2

@dataclass(frozen=True)
class OrbitalAnalysisResult:
    """Structured orbital-analysis result shared by NBO and VB drivers."""

    ao_to_atom: np.ndarray = field(default=None)
    ao_to_l: np.ndarray = field(default=None)
    overlap: np.ndarray = field(default=None)
    density: np.ndarray = field(default=None)
    nao_data: object = None
    spin_data: object = None
    mo_analysis: dict = field(default_factory=dict)
    orbital_candidates: list = field(default_factory=list)
    equivalent_atom_groups: tuple = field(default_factory=tuple)
    local_frames: tuple = field(default_factory=tuple)


# Backward-compatible name used by early OrbitalAnalyzer callers.
OrbitalAnalyzerDiagnostics = OrbitalAnalysisResult


class OrbitalAnalyzer:
    """Centralized orbital analysis and classification for VeloxChem."""

    def __init__(self, molecule, basis, mol_orbs=None, options=None, constraints=None):
        self.molecule = molecule
        self.basis = basis
        self.mol_orbs = mol_orbs
        self.options = options or OrbitalAnalyzerOptions()
        self.constraints = constraints
        self.results = None

    def run(self, mol_orbs=None, constraints=None):
        """Run the shared orbital analysis.

        Without molecular orbitals this returns only AO maps. With molecular
        orbitals it also returns NAO/NPA data, optional MO analysis, spin data,
        and shared NBO-like orbital candidates.
        """

        mol_orbs = self.mol_orbs if mol_orbs is None else mol_orbs
        constraints = self.constraints if constraints is None else constraints
        ao_to_atom, ao_to_l = self.ao_shell_map(self.molecule, self.basis)

        if mol_orbs is None:
            self.results = OrbitalAnalysisResult(
                ao_to_atom=ao_to_atom,
                ao_to_l=ao_to_l,
            )
            return self.results

        nbo = self._nbo_helpers()
        overlap = np.array(compute_overlap_integrals(self.molecule, self.basis))
        density = nbo._build_density_matrix(mol_orbs)
        equivalent_groups = nbo._find_equivalent_atom_groups(self.molecule)
        local_frames = nbo._build_local_atom_frames(self.molecule)
        nao_data = nbo._build_orthonormal_nao_data(
            overlap,
            density,
            ao_to_atom,
            ao_to_l,
            self.molecule.number_of_atoms(),
            equivalent_groups,
            local_frames,
        )
        spin_data = nbo._build_spin_nao_data(mol_orbs, overlap, nao_data)

        mo_analysis = {}
        if self.options.include_mo_analysis:
            mo_analysis = nbo._build_mo_nao_analysis(
                self.molecule,
                mol_orbs,
                overlap,
                nao_data,
                top_n=self.options.mo_analysis_top,
                threshold=self.options.mo_analysis_threshold,
            )

        orbital_candidates = []
        if self.options.include_nbo_candidates:
            from .orbitalclassifier import _build_nbo_candidates

            orbital_candidates = _build_nbo_candidates(
                self.molecule,
                nao_data,
                spin_data=spin_data,
                lone_pair_min_occupation=self.options.lone_pair_min_occupation,
                rydberg_max_occupation=self.options.rydberg_max_occupation,
                bond_min_occupation=self.options.bond_min_occupation,
                bond_min_atom_weight=self.options.bond_min_atom_weight,
                pi_min_occupation=self.options.pi_min_occupation,
                conjugated_pi_max_path=self.options.conjugated_pi_max_path,
                constraints=constraints,
            )

        self.results = OrbitalAnalysisResult(
            ao_to_atom=ao_to_atom,
            ao_to_l=ao_to_l,
            overlap=overlap,
            density=density,
            nao_data=nao_data,
            spin_data=spin_data,
            mo_analysis=mo_analysis,
            orbital_candidates=orbital_candidates,
            equivalent_atom_groups=equivalent_groups,
            local_frames=local_frames,
        )
        return self.results

    analyze = run

    def classify_nao_data(self, nao_data, spin_data=None, constraints=None):
        """Classify an existing NAO payload into shared orbital candidates."""

        constraints = self.constraints if constraints is None else constraints
        from .orbitalclassifier import _build_nbo_candidates

        orbital_candidates = _build_nbo_candidates(
            self.molecule,
            nao_data,
            spin_data=spin_data,
            lone_pair_min_occupation=self.options.lone_pair_min_occupation,
            rydberg_max_occupation=self.options.rydberg_max_occupation,
            bond_min_occupation=self.options.bond_min_occupation,
            bond_min_atom_weight=self.options.bond_min_atom_weight,
            pi_min_occupation=self.options.pi_min_occupation,
            conjugated_pi_max_path=self.options.conjugated_pi_max_path,
            constraints=constraints,
        )
        ao_to_atom, ao_to_l = self.ao_shell_map(self.molecule, self.basis)
        self.results = OrbitalAnalysisResult(
            ao_to_atom=ao_to_atom,
            ao_to_l=ao_to_l,
            nao_data=nao_data,
            spin_data=spin_data,
            orbital_candidates=orbital_candidates,
        )
        return self.results

    @staticmethod
    def ao_shell_map(molecule, basis):
        l_map = {'s': 0, 'p': 1, 'd': 2, 'f': 3, 'g': 4, 'h': 5}
        ao_to_atom = []
        ao_to_l = []
        for entry in basis.get_ao_basis_map(molecule):
            tokens = entry.split()
            atom_index = int(tokens[0]) - 1
            shell_label = tokens[2]
            shell_type = next(char for char in shell_label if char.isalpha())
            ao_to_atom.append(atom_index)
            ao_to_l.append(l_map.get(shell_type, 0))
        return np.array(ao_to_atom, dtype=int), np.array(ao_to_l, dtype=int)

    @staticmethod
    def _nbo_helpers():
        return importlib.import_module(f'{__package__}.nbodriver')
