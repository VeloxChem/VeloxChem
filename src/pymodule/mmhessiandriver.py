import copy

import numpy as np

try:
    import openmm as mm
    import openmm.unit as unit
except ImportError:
    pass


class MMHessianDriver:
    """
    Computes 3x3 partial MM Hessian blocks for a given internal coordinate
    using controlled OpenMM parameter settings, as required by PHF.

    For each internal coordinate two evaluations are needed:

    Evaluation 1 — h_unit  (Eq. 10 / the direction matrix):
        Set the target force constant to 1.0, zero all partial charges
        and all vdW well-depths in the NonbondedForce.
        The resulting 3x3 partial Hessian is h_unit.
        Because energy is linear in k, h_unit = H_d/k evaluated at k=1.

    Evaluation 2 — h0  (the residual H'_0):
        Set the target force constant to 0.0.
        Restore all other force constants to their current best-estimate
        values (updated by the orchestrator after each fitting stage).
        Restore all nonbonded parameters to their original values.
        The resulting 3x3 partial Hessian is H'_0.

    Hessians are computed numerically via central finite differences of
    forces (OpenMM does not expose analytical second derivatives).

    Units throughout: positions in nm, forces in kJ/mol/nm internally,
    all returned Hessian values in Hartree/Bohr^2.
    """

    _DISPLACEMENT = 1e-4  # nm — step size for finite differences

    # Conversion factor: kJ/mol/nm^2 -> Hartree/Bohr^2
    # Computed lazily to avoid importing veloxchemlib at module level.
    @staticmethod
    def _to_hartree_bohr2() -> float:
        from .veloxchemlib import bohr_in_angstrom, hartree_in_kjpermol
        bohr_to_nm = bohr_in_angstrom() * 0.1
        return bohr_to_nm**2 / hartree_in_kjpermol()

    def __init__(self):
        self._system = None
        self._positions_nm = None  # np.ndarray shape (N, 3), nanometres

        # Current best-estimate force constants, updated by PHFParameterizer
        # after each fitting stage so that H'_0 always uses up-to-date values.
        self._dihedral_constants = {}  # {(i,j,k,l): k_d}  kJ/mol
        self._angle_constants = {}  # {(i,j,k):   k_a}  kJ/mol/rad^2
        # Bond constants are not needed as bonds are the final stage.

    # ------------------------------------------------------------------
    # Setup
    # ------------------------------------------------------------------

    def set_system(self, openmm_system: 'mm.System', positions_nm: np.ndarray):
        """
        Store the base OpenMM system and positions.

        Parameters
        ----------
        openmm_system : openmm.System
            The system produced by VeloxChem's write_openmm_files /
            app.ForceField.createSystem pipeline. Never modified in-place.
        positions_nm : np.ndarray, shape (N, 3)
            Atomic positions in nanometres.
        """
        self._system = openmm_system
        self._positions_nm = positions_nm.copy()

    def update_dihedral_constants(self, dihedral_constants: dict):
        """Called by PHFParameterizer after Stage 1 (dihedral fitting)."""
        self._dihedral_constants = dict(dihedral_constants)

    def update_angle_constants(self, angle_constants: dict):
        """Called by PHFParameterizer after Stage 2 (angle fitting)."""
        self._angle_constants = dict(angle_constants)

    # ------------------------------------------------------------------
    # Public interface
    # ------------------------------------------------------------------

    def compute_full_hessian(self) -> np.ndarray:
        """
        Compute the full 3N x 3N numerical Hessian of the stored MM system
        using central finite differences.

        Rather than creating a new OpenMM Context per atom pair (as
        _numerical_partial_hessian does), a single Context is created and
        each atom is displaced once per Cartesian component.  Each
        displacement yields forces on all atoms simultaneously, so the
        entire Hessian is assembled in 6N force evaluations instead of
        6N^2.

        The result is symmetrised to remove numerical noise.

        Returns
        -------
        np.ndarray, shape (3N, 3N), units kJ/mol/nm^2
        """
        n_atoms = self._positions_nm.shape[0]
        n_dof = 3 * n_atoms
        delta = self._DISPLACEMENT

        integrator = mm.VerletIntegrator(0.001)
        platform = mm.Platform.getPlatformByName('Reference')
        context = mm.Context(self._system, integrator, platform)

        pos_vec3 = [mm.Vec3(*self._positions_nm[n]) for n in range(n_atoms)]
        context.setPositions(pos_vec3)

        hessian = np.zeros((n_dof, n_dof))

        for i in range(n_atoms):
            for p in range(3):
                # +delta
                pos_fwd = list(pos_vec3)
                v = list(pos_fwd[i])
                v[p] += delta
                pos_fwd[i] = mm.Vec3(*v)
                context.setPositions(pos_fwd)
                f_fwd = (context.getState(getForces=True).getForces(
                    asNumpy=True))

                # -delta
                pos_bwd = list(pos_vec3)
                v = list(pos_bwd[i])
                v[p] -= delta
                pos_bwd[i] = mm.Vec3(*v)
                context.setPositions(pos_bwd)
                f_bwd = (context.getState(getForces=True).getForces(
                    asNumpy=True))

                # H[3i+p, 3j+s] = -(F_fwd[j,s] - F_bwd[j,s]) / (2*delta)
                for j in range(n_atoms):
                    for s in range(3):
                        fwd_s = f_fwd[j][s].value_in_unit(
                            unit.kilojoules_per_mole / unit.nanometer)
                        bwd_s = f_bwd[j][s].value_in_unit(
                            unit.kilojoules_per_mole / unit.nanometer)
                        hessian[3 * i + p,
                                3 * j + s] = -(fwd_s - bwd_s) / (2.0 * delta)

        del context, integrator

        # Symmetrise then convert kJ/mol/nm^2 -> Hartree/Bohr^2
        return 0.5 * (hessian + hessian.T) * self._to_hartree_bohr2()

    def compute_h_unit(self, coord_type: str,
                       atom_indices: tuple) -> np.ndarray:
        """
        Compute h_unit: the 3x3 partial Hessian with the target force
        constant set to 1.0 and all nonbonded parameters zeroed.

        Parameters
        ----------
        coord_type : str
            'dihedral', 'angle', or 'bond'.
        atom_indices : tuple
            Full atom index tuple — (i,j,k,l), (i,j,k), or (i,j).

        Returns
        -------
        np.ndarray, shape (3, 3), units Hartree/Bohr^2
        """
        system = copy.deepcopy(self._system)
        self._zero_nonbonded(system)
        self._zero_all_bonded(system)
        self._set_target_force_constant(system, coord_type, atom_indices, k=1.0)
        pair = self._terminal_pair(coord_type, atom_indices)
        return self._numerical_partial_hessian(system, pair)

    def compute_h0(self, coord_type: str, atom_indices: tuple) -> np.ndarray:
        """
        Compute H'_0: the 3x3 partial Hessian with the target force
        constant zeroed and all other parameters at current best estimates.

        Parameters
        ----------
        coord_type : str
            'dihedral', 'angle', or 'bond'.
        atom_indices : tuple
            Full atom index tuple.

        Returns
        -------
        np.ndarray, shape (3, 3), units Hartree/Bohr^2
        """
        system = copy.deepcopy(self._system)
        self._apply_current_best_estimates(system)
        self._set_target_force_constant(system, coord_type, atom_indices, k=0.0)
        pair = self._terminal_pair(coord_type, atom_indices)
        return self._numerical_partial_hessian(system, pair)

    # ------------------------------------------------------------------
    # Private helpers — parameter manipulation
    # ------------------------------------------------------------------

    @staticmethod
    def _terminal_pair(coord_type: str, atom_indices: tuple) -> tuple:
        """Return the (a, b) terminal atom pair used for block extraction."""
        if coord_type == 'dihedral':
            return (atom_indices[0], atom_indices[3])
        elif coord_type == 'angle':
            return (atom_indices[0], atom_indices[2])
        elif coord_type == 'bond':
            return (atom_indices[0], atom_indices[1])
        else:
            raise ValueError(f"Unknown coord_type: {coord_type!r}")

    @staticmethod
    def _zero_nonbonded(system: 'mm.System'):
        """
        Set all partial charges to 0.0 and all vdW well-depths (epsilon)
        to 0.0 in the NonbondedForce. sigma is left unchanged to avoid
        division-by-zero inside OpenMM.
        """
        for force in system.getForces():
            if isinstance(force, mm.NonbondedForce):
                for i in range(force.getNumParticles()):
                    _, sigma, _ = force.getParticleParameters(i)
                    force.setParticleParameters(i, 0.0, sigma, 0.0)
                for i in range(force.getNumExceptions()):
                    idx1, idx2, _, sigma, _ = force.getExceptionParameters(i)
                    force.setExceptionParameters(i, idx1, idx2, 0.0, sigma, 0.0)

    @staticmethod
    def _zero_all_bonded(system: 'mm.System'):
        """
        Set all bonded force constants to 0.0.
        Used when computing h_unit so that only the target term contributes.
        Equilibrium values are left unchanged.
        """
        for force in system.getForces():
            if isinstance(force, mm.HarmonicBondForce):
                for i in range(force.getNumBonds()):
                    a1, a2, r0, _ = force.getBondParameters(i)
                    force.setBondParameters(
                        i, a1, a2, r0,
                        0.0 * unit.kilojoules_per_mole / unit.nanometer**2)

            elif isinstance(force, mm.HarmonicAngleForce):
                for i in range(force.getNumAngles()):
                    a1, a2, a3, theta0, _ = force.getAngleParameters(i)
                    force.setAngleParameters(
                        i, a1, a2, a3, theta0,
                        0.0 * unit.kilojoules_per_mole / unit.radian**2)

            elif isinstance(force, mm.PeriodicTorsionForce):
                for i in range(force.getNumTorsions()):
                    a1, a2, a3, a4, per, phase, _ = force.getTorsionParameters(
                        i)
                    force.setTorsionParameters(i, a1, a2, a3, a4, per, phase,
                                               0.0 * unit.kilojoules_per_mole)

    def _set_target_force_constant(
        self,
        system: 'mm.System',
        coord_type: str,
        atom_indices: tuple,
        k: float,
    ):
        """Set the force constant of the target internal coordinate to k."""
        if coord_type == 'dihedral':
            self._set_dihedral_k(system, atom_indices, k)
        elif coord_type == 'angle':
            self._set_angle_k(system, atom_indices, k)
        elif coord_type == 'bond':
            self._set_bond_k(system, atom_indices, k)

    @staticmethod
    def _set_dihedral_k(system: 'mm.System', indices: tuple, k: float):
        i, j, kk, l = indices
        for force in system.getForces():
            if isinstance(force, mm.PeriodicTorsionForce):
                for idx in range(force.getNumTorsions()):
                    a1, a2, a3, a4, per, phase, _ = force.getTorsionParameters(
                        idx)
                    if (a1, a2, a3, a4) in [(i, j, kk, l), (l, kk, j, i)]:
                        force.setTorsionParameters(
                            idx,
                            a1,
                            a2,
                            a3,
                            a4,
                            per,
                            phase,
                            k * unit.kilojoules_per_mole,
                        )

    @staticmethod
    def _set_angle_k(system: 'mm.System', indices: tuple, k: float):
        i, j, kk = indices
        for force in system.getForces():
            if isinstance(force, mm.HarmonicAngleForce):
                for idx in range(force.getNumAngles()):
                    a1, a2, a3, theta0, _ = force.getAngleParameters(idx)
                    if (a1, a2, a3) in [(i, j, kk), (kk, j, i)]:
                        force.setAngleParameters(
                            idx,
                            a1,
                            a2,
                            a3,
                            theta0,
                            k * unit.kilojoules_per_mole / unit.radian**2,
                        )

    @staticmethod
    def _set_bond_k(system: 'mm.System', indices: tuple, k: float):
        i, j = indices
        for force in system.getForces():
            if isinstance(force, mm.HarmonicBondForce):
                for idx in range(force.getNumBonds()):
                    a1, a2, r0, _ = force.getBondParameters(idx)
                    if {a1, a2} == {i, j}:
                        force.setBondParameters(
                            idx,
                            a1,
                            a2,
                            r0,
                            k * unit.kilojoules_per_mole / unit.nanometer**2,
                        )

    def _apply_current_best_estimates(self, system: 'mm.System'):
        """
        Write the best-estimate force constants determined so far into the
        system copy before computing H'_0. Called in compute_h0 only.
        """
        for (i, j, kk, l), k_val in self._dihedral_constants.items():
            self._set_dihedral_k(system, (i, j, kk, l), k_val)
        for (i, j, kk), k_val in self._angle_constants.items():
            self._set_angle_k(system, (i, j, kk), k_val)

    # ------------------------------------------------------------------
    # Private helpers — numerical Hessian
    # ------------------------------------------------------------------

    def _numerical_partial_hessian(
        self,
        system: 'mm.System',
        atom_pair: tuple,
    ) -> np.ndarray:
        """
        Compute the 3x3 partial Hessian d^2E / (d r_{atom_i} d r_{atom_j})
        numerically via central finite differences of forces.

        For each Cartesian component p of atom_i, displace ± delta and
        read forces on atom_j in all three components s:

            H[p, s] = -dF_{j,s} / d r_{i,p}
                    ≈ -(F_{j,s}(+δ) - F_{j,s}(-δ)) / (2δ)

        Parameters
        ----------
        system : openmm.System
            Already prepared system copy (parameters set as needed).
        atom_pair : tuple
            (atom_i, atom_j) — zero-based indices.

        Returns
        -------
        np.ndarray, shape (3, 3), units Hartree/Bohr^2
        """
        integrator = mm.VerletIntegrator(0.001)
        platform = mm.Platform.getPlatformByName('Reference')
        context = mm.Context(system, integrator, platform)

        # Build position list as Vec3 objects (OpenMM requirement)
        pos_vec3 = [
            mm.Vec3(*self._positions_nm[n])
            for n in range(self._positions_nm.shape[0])
        ]
        context.setPositions(pos_vec3)

        atom_i, atom_j = atom_pair
        delta = self._DISPLACEMENT
        h_block = np.zeros((3, 3))

        for p in range(3):
            # +delta displacement of atom_i along component p
            pos_fwd = list(pos_vec3)
            v = list(pos_fwd[atom_i])
            v[p] += delta
            pos_fwd[atom_i] = mm.Vec3(*v)
            context.setPositions(pos_fwd)
            f_fwd = (context.getState(getForces=True).getForces(asNumpy=True))

            # -delta displacement
            pos_bwd = list(pos_vec3)
            v = list(pos_bwd[atom_i])
            v[p] -= delta
            pos_bwd[atom_i] = mm.Vec3(*v)
            context.setPositions(pos_bwd)
            f_bwd = (context.getState(getForces=True).getForces(asNumpy=True))

            # Central difference: H[p,s] = -(F_fwd - F_bwd) / (2*delta)
            # Forces from OpenMM are in kJ/mol/nm; delta is in nm
            # so the result is in kJ/mol/nm^2.
            for s in range(3):
                fwd_s = f_fwd[atom_j][s].value_in_unit(
                    unit.kilojoules_per_mole / unit.nanometer)
                bwd_s = f_bwd[atom_j][s].value_in_unit(
                    unit.kilojoules_per_mole / unit.nanometer)
                h_block[p, s] = -(fwd_s - bwd_s) / (2.0 * delta)

        del context, integrator
        return h_block * self._to_hartree_bohr2()
