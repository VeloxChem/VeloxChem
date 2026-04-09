from mpi4py import MPI
import numpy as np
import sys
import os
import tempfile

try:
    import openmm as mm
    import openmm.app as app
    import openmm.unit as unit
except ImportError:
    pass

from .veloxchemlib import mpi_master
from .outputstream import OutputStream
from .mmhessiandriver import MMHessianDriver
from .partial_hessian_extractor import PartialHessianExtractor


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

    ``'impropers'``
        ``{(i, j, k, l): k_imp}``  —  kJ/mol

    ``'bond_equilibria'``
        ``{(i, j): r_eq}``  —  nanometres

    ``'angle_equilibria'``
        ``{(i, j, k): theta_eq}``  —  degrees

    All indices are zero-based, matching the convention used throughout
    MMForceFieldGenerator.
    """

    def get_reference(self):
        """
        Get string for reference paper.
        """

        return 'Y. Wang, M. Ozhgibesov, and K. Hirao, J. Comput. Chem. 37, 2349-2359 (2016)'

    def __init__(self, comm=None, ostream=None):
        """
        Initializes force field generator.
        """

        if comm is None:
            comm = MPI.COMM_WORLD

        if ostream is None:
            if comm.Get_rank() == mpi_master():
                ostream = OutputStream(sys.stdout)
            else:
                ostream = OutputStream(None)

        # mpi information
        self.comm = comm
        self.rank = self.comm.Get_rank()
        self.nodes = self.comm.Get_size()

        # output stream
        self.ostream = ostream

        self._bonds = None
        self._angles = None
        self._dihedrals = None
        self._impropers = None

        self._bond_constants = {}
        self._angle_constants = {}
        self._dihedral_constants = {}
        self._improper_constants = {}
        self._deferred_angles = set()

        self._extractor = PartialHessianExtractor()
        self._engine = MMHessianDriver()

    # ------------------------------------------------------------------
    # Public entry point
    # ------------------------------------------------------------------

    def compute(self,
                ff_gen,
                hessian: np.ndarray,
                equivalent_atoms: list,
                katachi=True) -> dict:
        """
        Run the full PHF fitting procedure.

        Parameters
        ----------
        ff_gen : veloxchem.MMForceFieldGenerator
            Fully initialised force field generator.
        hessian : np.ndarray, shape (3N, 3N)
            QM Hessian in atomic units (Hartree/Bohr^2). Passed directly
            to the fitting stages without unit conversion, since the MM
            partial Hessian blocks from ArtificialMMHessianEngine are also
            returned in Hartree/Bohr^2.

        Returns
        -------
        dict with keys 'bonds', 'angles', 'dihedrals', 'impropers'.
        """
        self._extract_topology(ff_gen)

        openmm_system, positions_nm = self._build_openmm_system(ff_gen)
        self._engine.set_system(openmm_system, positions_nm)

        self._fit_dihedrals(hessian)
        self._engine.update_dihedral_constants(self._dihedral_constants)

        self._fit_angles(hessian)
        if equivalent_atoms is not None:
            self._average_equivalent_angle_constants(equivalent_atoms)
        self._engine.update_angle_constants(self._angle_constants)

        self._fit_impropers(hessian)
        self._engine.update_improper_constants(self._improper_constants)

        self._fit_bonds(hessian)
        if equivalent_atoms is not None:
            self._average_equivalent_bond_constants(equivalent_atoms)
        self._engine.update_bond_constants(self._bond_constants)

        results = {
            'bonds': dict(self._bond_constants),
            'angles': dict(self._angle_constants),
            'dihedrals': dict(self._dihedral_constants),
            'impropers': dict(self._improper_constants),
        }

        if katachi:
            bond_eq, angle_eq = self.apply_katachi_bonds_angles(
                self._bonds, self._angles, positions_nm, openmm_system)
            if equivalent_atoms is not None:
                self._average_equivalent_bond_equilibria(
                    bond_eq, equivalent_atoms)
                self._average_equivalent_angle_equilibria(
                    angle_eq, equivalent_atoms)
            for key in angle_eq.keys():
                angle_eq[key] = angle_eq[key] * 180 / np.pi
            results.update({
                'bonds_equilibria': bond_eq,
                'angles_equilibria': angle_eq,
                'positions_nm': positions_nm,
            })

        return results

    def _extract_topology(self, ff_gen):
        """
        Pull bond, angle, dihedral, and improper lists from MMForceFieldGenerator.

        ff_gen.bonds, ff_gen.angles, ff_gen.dihedrals, ff_gen.impropers are
        dicts keyed by atom-index tuples (zero-based) populated by create_topology().
        """
        self._bonds = list(ff_gen.bonds.keys())
        self._angles = list(ff_gen.angles.keys())
        self._dihedrals = list(ff_gen.dihedrals.keys())
        self._impropers = list(ff_gen.impropers.keys())

    def _group_all_coupled(self) -> tuple:
        """
        Partition dihedrals and angles into independent and coupled groups.

        For each terminal atom pair (i, l) we collect all dihedrals and
        angles that contribute to the H_{il} partial Hessian block:
          - A dihedral i-j-k-l contributes with terminal pair (i, l).
          - An angle i-j-l contributes with terminal pair (i, l).

        Three cases arise:
          - One dihedral, no angle  -> acyclic, scalar solve (Eq. 13).
          - N>1 dihedrals, no angle -> 6-membered ring, NxN solve (Eq. 17).
          - Any dihedrals + angles  -> 5-membered ring, mixed NxN solve.

        Returns
        -------
        acyclic_dihedrals : list of dihedral tuples
        coupled_dihedrals : dict mapping (i, l) -> [dihedral, ...]
        coupled_mixed : dict mapping (i, l) -> {'dihedrals': [...], 'angles': [...]}
        deferred_angles : set of angle tuples absorbed into coupled_mixed
        """
        angle_by_terminal = {}
        for angle in self._angles:
            i, _j, l = angle
            for key in [(i, l), (l, i)]:
                angle_by_terminal.setdefault(key, []).append(angle)

        dihedral_by_terminal = {}
        for dihedral in self._dihedrals:
            i, _j, _k, l = dihedral
            dihedral_by_terminal.setdefault((i, l), []).append(dihedral)

        acyclic_dihedrals = []
        coupled_dihedrals = {}
        coupled_mixed = {}
        deferred_angles = set()

        for key, dihedrals in dihedral_by_terminal.items():
            angles = angle_by_terminal.get(key, [])

            if not angles:
                if len(dihedrals) == 1:
                    acyclic_dihedrals.append(dihedrals[0])
                else:
                    coupled_dihedrals[key] = dihedrals
            else:
                coupled_mixed[key] = {
                    'dihedrals': dihedrals,
                    'angles': angles,
                }
                for angle in angles:
                    deferred_angles.add(angle)

        return acyclic_dihedrals, coupled_dihedrals, coupled_mixed, deferred_angles

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

        with tempfile.TemporaryDirectory() as tmpdir:
            xml_file = os.path.join(tmpdir, f'{mol_name}.xml')
            pdb_file = os.path.join(tmpdir, f'{mol_name}.pdb')

            # write_openmm_files writes a residue XML and a PDB file,
            # exactly as OpenMMDynamics.create_system_from_molecule does.
            ff_gen.write_openmm_files(os.path.join(tmpdir, mol_name), mol_name)

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
        Stage 1: fit k_d for all dihedrals, and simultaneously k_a for any
        angles coupled with a dihedral via a shared terminal pair (5-ring).

        Three dispatch paths:
          - Acyclic: scalar solve (Eq. 13).
          - Coupled dihedrals only: NxN linear solve (Eq. 17).
          - Mixed dihedral+angle: NxN linear solve (Eq. 17 generalised).
            Coupled angles are stored in _deferred_angles and skipped in
            _fit_angles.
        """
        acyclic, coupled_dihedrals, coupled_mixed, deferred = \
            self._group_all_coupled()
        self._deferred_angles = deferred

        # Acyclic dihedrals: scalar path (Eq. 13)
        for dihedral in acyclic:
            i, _j, _k, l = dihedral
            h_qm = self._extractor.extract(hessian, i, l)
            h_unit = self._engine.compute_h_unit('dihedral', dihedral)
            if not np.any(h_unit):  # dihedral absent from OpenMM system
                continue
            h0 = self._engine.compute_h0('dihedral', dihedral)
            self._dihedral_constants[dihedral] = self.solve_dihedral(
                h_qm, h0, h_unit)

        # Coupled dihedrals only: 6-membered ring path (Eq. 17)
        for (i, l), group in coupled_dihedrals.items():
            h_units_all = [
                self._engine.compute_h_unit('dihedral', d) for d in group
            ]
            present = [(d, hu) for d, hu in zip(group, h_units_all)
                       if np.any(hu)]
            if not present:
                continue
            present_dihedrals, h_units = zip(*present)
            h_qm = self._extractor.extract(hessian, i, l)
            zero_targets = [('dihedral', d) for d in present_dihedrals]
            h0 = self._engine.compute_h0_zeroing(zero_targets, (i, l))
            if len(present_dihedrals) == 1:
                self._dihedral_constants[present_dihedrals[0]] = (
                    self.solve_dihedral(h_qm, h0, h_units[0]))
            else:
                k_values = self.solve_coupled(h_qm, h0, list(h_units))
                for dihedral, k_d in zip(present_dihedrals, k_values):
                    self._dihedral_constants[dihedral] = float(k_d)

        # Mixed dihedral+angle: 5-membered ring path
        for (i, l), members in coupled_mixed.items():
            dihedrals = members['dihedrals']
            angles = members['angles']

            d_hunits = [
                self._engine.compute_h_unit('dihedral', d) for d in dihedrals
            ]
            a_hunits = [self._engine.compute_h_unit('angle', a) for a in angles]

            present_d = [(d, hu) for d, hu in zip(dihedrals, d_hunits)
                         if np.any(hu)]
            present_a = [(a, hu) for a, hu in zip(angles, a_hunits)
                         if np.any(hu)]

            # angles whose dihedral partner vanished must not be deferred
            lost_angles = {
                a
                for a, hu in zip(angles, a_hunits) if not np.any(hu)
            }
            self._deferred_angles -= lost_angles

            if not present_d and not present_a:
                continue

            all_present = present_d + present_a
            coords, h_units = zip(*all_present)
            dihedrals_set = set(dihedrals)
            zero_targets = [('dihedral', c) if c in dihedrals_set else
                            ('angle', c) for c in coords]
            h_qm = self._extractor.extract(hessian, i, l)
            h0 = self._engine.compute_h0_zeroing(zero_targets, (i, l))

            if len(coords) == 1:
                coord = coords[0]
                if coord in dihedrals_set:
                    self._dihedral_constants[coord] = self.solve_dihedral(
                        h_qm, h0, h_units[0])
                else:
                    self._angle_constants[coord] = self.solve_angle(
                        h_qm, h0, h_units[0])
            else:
                k_values = self.solve_coupled(h_qm, h0, list(h_units))
                n_d = len(present_d)
                for (d, _), k_d in zip(present_d, k_values[:n_d]):
                    self._dihedral_constants[d] = float(k_d)
                for (a, _), k_a in zip(present_a, k_values[n_d:]):
                    self._angle_constants[a] = float(k_a)

    def _fit_angles(self, hessian: np.ndarray):
        """
        Stage 2: fit k_a for each angle i-j-k, except those already fitted
        simultaneously with a dihedral during Stage 1 (5-ring coupled case).
        H'_0 includes dihedral contributions via already-fitted k_d values.
        No nonbonded contribution: AMBER has zero 1-3 interactions.
        """
        for angle in self._angles:
            if angle in self._deferred_angles:
                continue
            i, _j, k = angle
            h_qm = self._extractor.extract(hessian, i, k)
            h_unit = self._engine.compute_h_unit('angle', angle)
            h0 = self._engine.compute_h0('angle', angle)
            self._angle_constants[angle] = self.solve_angle(h_qm, h0, h_unit)

    def _fit_impropers(self, hessian: np.ndarray):
        """
        Stage 3: fit k_imp for each improper torsion i-j-k-l.
        Terminal atom pair used: (i, l), same convention as proper dihedrals.
        H'_0 includes dihedral and angle contributions via already-fitted
        k_d and k_a values. No nonbonded contribution: AMBER impropers
        do not have 1-4 scaled interactions.
        """
        for improper in self._impropers:
            i, _j, _k, l = improper
            h_qm = self._extractor.extract(hessian, i, l)
            h_unit = self._engine.compute_h_unit('improper', improper)
            if not np.any(h_unit):  # improper absent from OpenMM system
                continue
            h0 = self._engine.compute_h0('improper', improper)
            k_imp = self.solve_improper(h_qm, h0, h_unit)
            self._improper_constants[improper] = k_imp

    def _fit_bonds(self, hessian: np.ndarray):
        """
        Stage 4: fit k_b for each bond i-j.
        Terminal atom pair used: (i, j).
        H'_0 includes angle, dihedral, and improper contributions via
        already-fitted k_a, k_d, and k_imp values.
        """
        for bond in self._bonds:
            i, j = bond
            h_qm = self._extractor.extract(hessian, i, j)
            h_unit = self._engine.compute_h_unit('bond', bond)
            h0 = self._engine.compute_h0('bond', bond)
            k_b = self.solve_bond(h_qm, h0, h_unit)
            self._bond_constants[bond] = k_b

    def _average_equivalent_angle_constants(self, equivalent_atoms: list):
        """
        Average angle force constants across symmetry-equivalent angles.

        Two angles (i, j, k) and (i', j', k') are equivalent when their
        central atoms share a label and their endpoint-label sets are equal:
        equivalent_atoms[j] == equivalent_atoms[j'] and
        {equivalent_atoms[i], equivalent_atoms[k]} == {equivalent_atoms[i'],
        equivalent_atoms[k']}.

        Modifies self._angle_constants in place.

        Parameters
        ----------
        equivalent_atoms : list of str
            One label per atom (zero-based). Atoms sharing a label are equivalent.
        """
        groups: dict = {}
        for angle in self._angle_constants:
            i, j, k = angle
            li, lj, lk = equivalent_atoms[i], equivalent_atoms[
                j], equivalent_atoms[k]
            key = (min(li, lk), lj, max(li, lk))
            groups.setdefault(key, []).append(angle)

        for group in groups.values():
            mean_k = float(np.mean([self._angle_constants[a] for a in group]))
            for a in group:
                self._angle_constants[a] = mean_k

    def _average_equivalent_bond_constants(self, equivalent_atoms: list):
        """
        Average bond force constants across symmetry-equivalent bonds.

        Two bonds (i, j) and (i', j') are equivalent when
        {equivalent_atoms[i], equivalent_atoms[j]} is the same frozenset.

        Modifies self._bond_constants in place.

        Parameters
        ----------
        equivalent_atoms : list of str
            One label per atom (zero-based). Atoms sharing a label are equivalent.
        """
        groups: dict = {}
        for bond in self._bond_constants:
            i, j = bond
            key = frozenset({equivalent_atoms[i], equivalent_atoms[j]})
            groups.setdefault(key, []).append(bond)

        for group in groups.values():
            mean_k = float(np.mean([self._bond_constants[b] for b in group]))
            for b in group:
                self._bond_constants[b] = mean_k

    def _average_equivalent_bond_equilibria(self, bond_eq: dict,
                                            equivalent_atoms: list):
        """
        Average bond equilibrium lengths across symmetry-equivalent bonds.

        Uses the same grouping key as _average_equivalent_bond_constants:
        frozenset of endpoint labels. Modifies bond_eq in place.

        Parameters
        ----------
        bond_eq : dict, {(i, j): r0_nm}
        equivalent_atoms : list of str
        """
        groups: dict = {}
        for bond in bond_eq:
            i, j = bond
            key = frozenset({equivalent_atoms[i], equivalent_atoms[j]})
            groups.setdefault(key, []).append(bond)

        for group in groups.values():
            mean_r = float(np.mean([bond_eq[b] for b in group]))
            for b in group:
                bond_eq[b] = mean_r

    def _average_equivalent_angle_equilibria(self, angle_eq: dict,
                                             equivalent_atoms: list):
        """
        Average angle equilibrium values across symmetry-equivalent angles.

        Uses the same grouping key as _average_equivalent_angle_constants:
        (min_endpoint_label, center_label, max_endpoint_label).
        Modifies angle_eq in place.

        Parameters
        ----------
        angle_eq : dict, {(i, j, k): theta0_rad}
        equivalent_atoms : list of str
        """
        groups: dict = {}
        for angle in angle_eq:
            i, j, k = angle
            li, lj, lk = equivalent_atoms[i], equivalent_atoms[
                j], equivalent_atoms[k]
            key = (min(li, lk), lj, max(li, lk))
            groups.setdefault(key, []).append(angle)

        for group in groups.values():
            mean_theta = float(np.mean([angle_eq[a] for a in group]))
            for a in group:
                angle_eq[a] = mean_theta

    def _solve_fc(
        self,
        h_qm: np.ndarray,
        h0: np.ndarray,
        h_unit: np.ndarray,
    ) -> float:
        """
        Core least-squares solve shared by all three coordinate types.

        Parameters
        ----------
        h_qm : np.ndarray, shape (3, 3)
            QM partial Hessian block for the terminal atom pair.
        h0 : np.ndarray, shape (3, 3)
            Residual MM partial Hessian: contributions from all terms
            except the one being fitted, evaluated with that term set to 0.
        h_unit : np.ndarray, shape (3, 3)
            Normalised MM partial Hessian: the target term evaluated with
            its force constant set to 1.0 and all nonbonded terms zeroed.

        Returns
        -------
        float
            Fitted force constant k.
        """
        numerator = np.sum((h_qm - h0) * h_unit)
        denominator = np.sum(h_unit**2)

        if abs(denominator) < 1e-30:
            raise ValueError(
                "Denominator in PHF solve is effectively zero. "
                "The h_unit block is all-zero for this internal coordinate, "
                "which means the geometry has this term contributing nothing.")

        return float(numerator / denominator)

    def solve_dihedral(
        self,
        h_qm: np.ndarray,
        h0: np.ndarray,
        h_unit: np.ndarray,
    ) -> float:
        """
        Solve for dihedral force constant k_d (Eq. 13).

        Terminal atom pair is (i, l) for dihedral i-j-k-l.
        h0 contains only the nonbonded 1-4 contribution between i and l.
        """
        return self._solve_fc(h_qm, h0, h_unit)

    def solve_angle(
        self,
        h_qm: np.ndarray,
        h0: np.ndarray,
        h_unit: np.ndarray,
    ) -> float:
        """
        Solve for angle-bending force constant k_a (Eq. 19).

        Terminal atom pair is (i, k) for angle i-j-k.
        h0 contains dihedral contributions involving both i and k,
        computed using k_d values already determined in Stage 1.
        Note: no nonbonded contribution because AMBER has zero 1-3 interactions.
        """
        return self._solve_fc(h_qm, h0, h_unit)

    def solve_improper(
        self,
        h_qm: np.ndarray,
        h0: np.ndarray,
        h_unit: np.ndarray,
    ) -> float:
        """
        Solve for improper-dihedral force constant k_imp.

        Terminal atom pair is (i, l) for improper i-j-k-l.
        h0 includes angle and dihedral contributions involving both i and l,
        using k_a and k_d values already determined in Stages 1 and 2.
        No nonbonded contribution: impropers don't have 1-4 interactions.
        """
        return self._solve_fc(h_qm, h0, h_unit)

    def solve_coupled(
        self,
        h_qm: np.ndarray,
        h0: np.ndarray,
        h_units: list,
    ) -> np.ndarray:
        """
        Solve simultaneously for N coupled force constants (Eq. 17 generalised).

        Handles any combination of coordinate types (dihedrals, angles, or
        mixed) that share the same terminal atom pair. Builds and solves:

            A @ k = b

        where  A[m, n] = sum_{p,s} h_unit_m[p,s] * h_unit_n[p,s]
               b[m]    = sum_{p,s} (H_QM - H'_0)[p,s] * h_unit_m[p,s]

        Parameters
        ----------
        h_qm : np.ndarray, shape (3, 3)
            QM partial Hessian block for the shared terminal atom pair.
        h0 : np.ndarray, shape (3, 3)
            Residual H'_0 with all coupled force constants zeroed.
        h_units : list of np.ndarray, each shape (3, 3)
            One h_unit block per coordinate being fitted, in order.
            For a mixed 5-ring case: dihedrals first, then angles.

        Returns
        -------
        np.ndarray, shape (N,)
            Fitted force constants in the same order as h_units.
        """
        n = len(h_units)
        residual = h_qm - h0

        A = np.zeros((n, n))
        b = np.zeros(n)
        for m in range(n):
            b[m] = np.sum(residual * h_units[m])
            for nn in range(n):
                A[m, nn] = np.sum(h_units[m] * h_units[nn])

        if np.linalg.cond(A) > 1e12:
            raise ValueError(
                "Singular matrix in coupled PHF solve. "
                "The h_unit blocks for the coupled coordinates are linearly "
                "dependent, which typically means two internal coordinates "
                "are geometrically indistinguishable at this geometry.")

        return np.linalg.solve(A, b)

    def solve_bond(
        self,
        h_qm: np.ndarray,
        h0: np.ndarray,
        h_unit: np.ndarray,
    ) -> float:
        """
        Solve for bond-stretching force constant k_b (Eq. 23).

        Terminal atom pair is (i, j) for bond i-j.
        h0 contains angle and dihedral contributions involving both i and j,
        using k_a and k_d values already determined in Stages 1 and 2.
        """
        return self._solve_fc(h_qm, h0, h_unit)

    def _run_mm_minimisation(self, positions_nm,
                             system: 'mm.System') -> np.ndarray:
        """
        Minimise MM energy and return optimised positions in nanometres,
        shape (N, 3).
        """
        integrator = mm.VerletIntegrator(0.001)
        platform = mm.Platform.getPlatformByName('Reference')
        context = mm.Context(system, integrator, platform)
        context.setPositions(positions_nm * unit.nanometer)
        mm.LocalEnergyMinimizer.minimize(context)
        state = context.getState(getPositions=True)
        pos = state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
        del context, integrator
        return pos

    def apply_katachi_bonds_angles(
        self,
        bonds: list,
        angles: list,
        positions_nm: np.ndarray,
        system,
        max_iter: int = 50,
    ) -> 'tuple[dict, dict]':
        """
        Iteratively adjust bond and angle equilibrium values so the MM
        geometry minimum matches the QM geometry (Katachi amendment,
        Wang et al. 2016, p. 2354).

        The MM system used for minimisation already contains the force
        constants fitted by the four PHF stages (dihedrals, angles,
        impropers, bonds), because this method is called after all fitting
        is complete and self._system holds those values.

        Convergence criteria (from the paper):
          bonds:  max|D_i| < 0.0001 Å = 1e-5 nm
          angles: max|D_i| < 0.0028 rad

        Parameters
        ----------
        bonds  : list of (i, j) tuples
        angles : list of (i, j, k) tuples

        Returns
        -------
        bond_equilibria  : {(i, j): r0}   nanometres
        angle_equilibria : {(i, j, k): θ₀} radians
        """
        import copy
        from .reactionsystembuilder import ReactionSystemBuilder

        BOND_THR = 1e-5  # nm  (0.0001 Å)
        ANGLE_THR = 0.002  # rad

        # QM reference geometry
        bond_eq_qm = {
            b:
            ReactionSystemBuilder.measure_length(positions_nm[b[0]],
                                                 positions_nm[b[1]])
            for b in bonds
        }
        angle_eq_qm = {
            a:
            ReactionSystemBuilder.measure_angle(positions_nm[a[0]],
                                                positions_nm[a[1]],
                                                positions_nm[a[2]],
                                                angle_unit='radian')
            for a in angles
        }
        # Working system copy — equilibria updated each iteration
        system = copy.deepcopy(system)
        self._engine._apply_current_best_estimates(system)

        for b in bonds:
            self._engine._set_bond_eq(system, b, bond_eq_qm[b])
        for a in angles:
            self._engine._set_angle_eq(system, a, angle_eq_qm[a])

        # Current equilibria (start from QM values already in the system)
        bond_eq = dict(bond_eq_qm)
        angle_eq = dict(angle_eq_qm)

        for i in range(max_iter):
            pos_mm = self._run_mm_minimisation(positions_nm, system)

            bond_mm = {
                b: ReactionSystemBuilder.measure_length(pos_mm[b[0]],
                                                        pos_mm[b[1]])
                for b in bonds
            }
            angle_mm = {
                a:
                ReactionSystemBuilder.measure_angle(pos_mm[a[0]],
                                                    pos_mm[a[1]],
                                                    pos_mm[a[2]],
                                                    angle_unit='radian')
                for a in angles
            }

            bond_d = {b: bond_eq_qm[b] - bond_mm[b] for b in bonds}
            angle_d = {a: angle_eq_qm[a] - angle_mm[a] for a in angles}

            max_bond_d = max((abs(v) for v in bond_d.values()), default=0.0)
            max_angle_d = max((abs(v) for v in angle_d.values()), default=0.0)

            if max_bond_d < BOND_THR and max_angle_d < ANGLE_THR:
                break

            for b in bonds:
                bond_eq[b] += bond_d[b]
                self._engine._set_bond_eq(system, b, bond_eq[b])
            for a in angles:
                angle_eq[a] += angle_d[a]
                self._engine._set_angle_eq(system, a, angle_eq[a])

        return bond_eq, angle_eq
