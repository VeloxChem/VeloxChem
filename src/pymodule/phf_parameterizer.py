import numpy as np

try:
    import openmm.app as app
    import openmm.unit as unit
except ImportError:
    pass

from .artificial_mm_hessian_engine import ArtificialMMHessianEngine
from .force_constant_solver import ForceConstantSolver
from .partial_hessian_extractor import PartialHessianExtractor

from .veloxchemlib import hartree_in_kjpermol, bohr_in_angstrom


class PHFParameterizer:
    """
    Partial Hessian Fitting (PHF) for deriving AMBER-type force constant
    parameters from a QM Hessian matrix.

    Reference:
        Wang, Ozhgibesov & Hirao,
        J. Comput. Chem. 2016, 37, 2349-2359. DOI: 10.1002/jcc.24457

    Usage
    -----
    ::

        phf = PHFParameterizer()
        result = phf.compute(ff_gen, hessian)

    Parameters passed to compute()
    --------------------------------
    ff_gen : veloxchem.MMForceFieldGenerator
        Fully initialised generator after create_topology() has been called.
        Must have ff_gen.molecule, ff_gen.bonds, ff_gen.angles,
        ff_gen.dihedrals populated.
    hessian : np.ndarray, shape (3N, 3N)
        Full Cartesian QM Hessian in atomic units (Hartree / Bohr^2),
        the same array produced by VeloxChem frequency / Hessian drivers.
        Conversion to kJ/mol/nm^2 is applied internally.

    Returns
    -------
    dict with keys:

    ``'bonds'``
        ``{(i, j): k_b}``  —  kJ/mol/nm^2

    ``'angles'``
        ``{(i, j, k): k_a}``  —  kJ/mol/rad^2

    ``'dihedrals'``
        ``{(i, j, k, l): k_d}``  —  kJ/mol

    All indices are zero-based, matching the convention used throughout
    MMForceFieldGenerator.
    """

    def __init__(self):
        self._bonds = None
        self._angles = None
        self._dihedrals = None

        self._bond_constants = {}
        self._angle_constants = {}
        self._dihedral_constants = {}

        self._extractor = PartialHessianExtractor()
        self._solver = ForceConstantSolver()
        self._engine = ArtificialMMHessianEngine()

    # ------------------------------------------------------------------
    # Public entry point
    # ------------------------------------------------------------------

    def compute(self, ff_gen, hessian: np.ndarray) -> dict:
        """
        Run the full PHF fitting procedure.

        Parameters
        ----------
        ff_gen : veloxchem.MMForceFieldGenerator
            Fully initialised force field generator.
        hessian : np.ndarray, shape (3N, 3N)
            QM Hessian in atomic units (Hartree / Bohr^2).

        Returns
        -------
        dict with keys 'bonds', 'angles', 'dihedrals'.
        """
        # Convert QM Hessian from Hartree/Bohr^2 to kJ/mol/nm^2.
        # The Seminario implementation in reparameterize() applies
        # hartree_in_kjpermol() / bohr_to_nm^2 for bond force constants,
        # so we use the same conversion factors here.
        bohr_to_nm = bohr_in_angstrom() * 0.1
        hessian_kjmol_nm2 = hessian * hartree_in_kjpermol() / bohr_to_nm**2

        self._extract_topology(ff_gen)
        self._check_for_rings()

        openmm_system, positions_nm = self._build_openmm_system(ff_gen)
        self._engine.set_system(openmm_system, positions_nm)

        self._fit_dihedrals(hessian_kjmol_nm2)
        self._engine.update_dihedral_constants(self._dihedral_constants)

        self._fit_angles(hessian_kjmol_nm2)
        self._engine.update_angle_constants(self._angle_constants)

        self._fit_bonds(hessian_kjmol_nm2)

        return {
            'bonds': dict(self._bond_constants),
            'angles': dict(self._angle_constants),
            'dihedrals': dict(self._dihedral_constants),
        }

    # ------------------------------------------------------------------
    # Topology extraction
    # ------------------------------------------------------------------

    def _extract_topology(self, ff_gen):
        """
        Pull bond, angle, and dihedral lists from MMForceFieldGenerator.

        ff_gen.bonds, ff_gen.angles, ff_gen.dihedrals are dicts keyed by
        atom-index tuples (zero-based) populated by create_topology().
        """
        self._bonds = list(ff_gen.bonds.keys())
        self._angles = list(ff_gen.angles.keys())
        self._dihedrals = list(ff_gen.dihedrals.keys())

    def _check_for_rings(self):
        """
        Detect cases where two dihedrals share the same terminal atom pair
        (i, l), requiring simultaneous solving of a 2x2 linear system
        (Eq. 17 in the paper). Not yet implemented.
        """
        seen = {}
        for dihedral in self._dihedrals:
            i, j, k, l = dihedral
            key = (i, l)
            if key in seen:
                raise NotImplementedError(
                    f"Ring coupling detected: dihedrals {seen[key]} and "
                    f"{dihedral} share terminal atom pair ({i}, {l}). "
                    "Simultaneous fitting for cyclic systems is not yet "
                    "implemented in PHFParameterizer.")
            seen[key] = dihedral

    # ------------------------------------------------------------------
    # OpenMM system construction
    # ------------------------------------------------------------------

    @staticmethod
    def _build_openmm_system(ff_gen):
        """
        Build an OpenMM System and extract positions from ff_gen.

        Replicates the file-based path used in
        OpenMMDynamics.create_system_from_molecule() for the gas-phase
        case, which is all we need here.

        Returns
        -------
        system : openmm.System
        positions_nm : np.ndarray, shape (N, 3), nanometres
        """
        mol_name = 'PHF_MOL'
        xml_file = f'{mol_name}.xml'
        pdb_file = f'{mol_name}.pdb'

        # write_openmm_files writes a residue XML and a PDB file,
        # exactly as OpenMMDynamics.create_system_from_molecule does.
        ff_gen.write_openmm_files(mol_name, mol_name)

        pdb = app.PDBFile(pdb_file)
        forcefield = app.ForceField(xml_file)
        system = forcefield.createSystem(
            pdb.topology,
            nonbondedMethod=app.NoCutoff,
            constraints=None,  # no constraints — we need all DOF for Hessian
        )

        # Positions from the QM-optimised geometry stored in ff_gen.molecule.
        # VeloxChem stores coordinates in Angstrom; convert to nm.
        coords_ang = ff_gen.molecule.get_coordinates_in_angstrom()
        positions_nm = np.array(coords_ang) / 10.0

        return system, positions_nm

    # ------------------------------------------------------------------
    # Fitting stages
    # ------------------------------------------------------------------

    def _fit_dihedrals(self, hessian: np.ndarray):
        """
        Stage 1: fit k_d for each dihedral i-j-k-l.
        Terminal atom pair used: (i, l).
        H'_0 contains only the nonbonded 1-4 contribution between i and l
        (all other bonded constants are zero at this stage).
        """
        for dihedral in self._dihedrals:
            i, _j, _k, l = dihedral
            h_qm = self._extractor.extract(hessian, i, l)
            h_unit = self._engine.compute_h_unit('dihedral', dihedral)
            h0 = self._engine.compute_h0('dihedral', dihedral)
            k_d = self._solver.solve_dihedral(h_qm, h0, h_unit)
            self._dihedral_constants[dihedral] = k_d

    def _fit_angles(self, hessian: np.ndarray):
        """
        Stage 2: fit k_a for each angle i-j-k.
        Terminal atom pair used: (i, k).
        H'_0 includes dihedral contributions via already-fitted k_d values.
        No nonbonded contribution: AMBER has zero 1-3 interactions.
        """
        for angle in self._angles:
            i, _j, k = angle
            h_qm = self._extractor.extract(hessian, i, k)
            h_unit = self._engine.compute_h_unit('angle', angle)
            h0 = self._engine.compute_h0('angle', angle)
            k_a = self._solver.solve_angle(h_qm, h0, h_unit)
            self._angle_constants[angle] = k_a

    def _fit_bonds(self, hessian: np.ndarray):
        """
        Stage 3: fit k_b for each bond i-j.
        Terminal atom pair used: (i, j).
        H'_0 includes both angle and dihedral contributions via
        already-fitted k_a and k_d values.
        """
        for bond in self._bonds:
            i, j = bond
            h_qm = self._extractor.extract(hessian, i, j)
            h_unit = self._engine.compute_h_unit('bond', bond)
            h0 = self._engine.compute_h0('bond', bond)
            k_b = self._solver.solve_bond(h_qm, h0, h_unit)
            self._bond_constants[bond] = k_b
