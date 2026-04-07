from mpi4py import MPI
import numpy as np
import sys

try:
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

    def compute(self, ff_gen, hessian: np.ndarray) -> dict:
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
        self._filter_unrepresented_dihedrals()

        self._fit_dihedrals(hessian)
        self._engine.update_dihedral_constants(self._dihedral_constants)

        self._fit_angles(hessian)
        self._engine.update_angle_constants(self._angle_constants)

        self._fit_impropers(hessian)
        self._engine.update_improper_constants(self._improper_constants)

        self._fit_bonds(hessian)

        return {
            'bonds': dict(self._bond_constants),
            'angles': dict(self._angle_constants),
            'dihedrals': dict(self._dihedral_constants),
            'impropers': dict(self._improper_constants),
        }

    # ------------------------------------------------------------------
    # Topology extraction
    # ------------------------------------------------------------------

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

    def _filter_unrepresented_dihedrals(self):
        """
        Remove from self._dihedrals any entry that has no corresponding
        PeriodicTorsionForce entry in the OpenMM system.

        This occurs when a dihedral is present in ff_gen.dihedrals (because
        the topology enumerates all four-atom paths) but the GAFF force field
        assigns it no torsion term — for example, the X-ca-c3-X methyl
        rotation has zero barrier and is simply omitted from the XML.
        For such dihedrals h_unit would be all-zero and PHF cannot fit them.
        They are silently dropped here; their barrier stays at zero (free
        rotation), which is the correct force-field answer.

        Must be called after _engine.set_system() so the OpenMM system is
        available via self._engine.get_represented_dihedrals().
        """
        represented = self._engine.get_represented_dihedrals()
        before = len(self._dihedrals)
        self._dihedrals = [
            d for d in self._dihedrals
            if d in represented or (d[3], d[2], d[1], d[0]) in represented
        ]
        skipped = before - len(self._dihedrals)
        if skipped > 0:
            self.ostream.print_info(
                f'PHF: skipped {skipped} dihedral(s) with no force field '
                f'representation (zero-barrier terms omitted from XML).')
            self.ostream.flush()

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
        import os
        import tempfile

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
            h0 = self._engine.compute_h0('dihedral', dihedral)
            self._dihedral_constants[dihedral] = self.solve_dihedral(
                h_qm, h0, h_unit)

        # Coupled dihedrals only: 6-membered ring path (Eq. 17)
        for (i, l), group in coupled_dihedrals.items():
            h_qm = self._extractor.extract(hessian, i, l)
            h_units = [
                self._engine.compute_h_unit('dihedral', d) for d in group
            ]
            zero_targets = [('dihedral', d) for d in group]
            h0 = self._engine.compute_h0_zeroing(zero_targets, (i, l))
            k_values = self.solve_coupled(h_qm, h0, h_units)
            for dihedral, k_d in zip(group, k_values):
                self._dihedral_constants[dihedral] = float(k_d)

        # Mixed dihedral+angle: 5-membered ring path
        for (i, l), members in coupled_mixed.items():
            dihedrals = members['dihedrals']
            angles = members['angles']
            h_qm = self._extractor.extract(hessian, i, l)
            h_units = (
                [self._engine.compute_h_unit('dihedral', d) for d in dihedrals]
                + [self._engine.compute_h_unit('angle', a) for a in angles]
            )
            zero_targets = (
                [('dihedral', d) for d in dihedrals]
                + [('angle', a) for a in angles]
            )
            h0 = self._engine.compute_h0_zeroing(zero_targets, (i, l))
            k_values = self.solve_coupled(h_qm, h0, h_units)
            n_d = len(dihedrals)
            for dihedral, k_d in zip(dihedrals, k_values[:n_d]):
                self._dihedral_constants[dihedral] = float(k_d)
            for angle, k_a in zip(angles, k_values[n_d:]):
                self._angle_constants[angle] = float(k_a)

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
            self._angle_constants[angle] = self.solve_angle(
                h_qm, h0, h_unit)

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

        if abs(np.linalg.det(A)) < 1e-30:
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