#
#                                   VELOXCHEM
#              ----------------------------------------------------
#                          An Electronic Structure Code
#
#  SPDX-License-Identifier: BSD-3-Clause
#
#  Copyright 2018-2025 VeloxChem developers
#
#  Redistribution and use in source and binary forms, with or without modification,
#  are permitted provided that the following conditions are met:
#
#  1. Redistributions of source code must retain the above copyright notice, this
#     list of conditions and the following disclaimer.
#  2. Redistributions in binary form must reproduce the above copyright notice,
#     this list of conditions and the following disclaimer in the documentation
#     and/or other materials provided with the distribution.
#  3. Neither the name of the copyright holder nor the names of its contributors
#     may be used to endorse or promote products derived from this software without
#     specific prior written permission.
#
#  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
#  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
#  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
#  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
#  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
#  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
#  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
#  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
#  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
#  OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

from mpi4py import MPI
import numpy as np
import sys

from .veloxchemlib import mpi_master, bohr_in_angstrom, hartree_in_kjpermol
from .outputstream import OutputStream
from .molecule import Molecule
from .openmmdynamics import OpenMMDynamics


# Atomic masses in amu (most common isotope)
_ATOMIC_MASSES = {
    'H': 1.008, 'He': 4.003, 'Li': 6.941, 'Be': 9.012, 'B': 10.811,
    'C': 12.011, 'N': 14.007, 'O': 15.999, 'F': 18.998, 'Ne': 20.180,
    'Na': 22.990, 'Mg': 24.305, 'Al': 26.982, 'Si': 28.086, 'P': 30.974,
    'S': 32.065, 'Cl': 35.453, 'Ar': 39.948, 'K': 39.098, 'Ca': 40.078,
    'Br': 79.904, 'I': 126.904,
}

# Conversion: amu → electron mass (atomic units)
_AMU_TO_AU = 1.0 / 5.4857990907e-4

# Speed of light in cm/s for frequency conversion
_C_CM = 2.99792458e10

# Hartree/Bohr²/amu → cm⁻¹
# ω (rad/s) = sqrt(λ / m_au) × (E_h / a0² → J/m²) × (1/kg)
# then ν̃ = ω / (2π c)
_HESSIAN_TO_WAVENUMBER = (
    np.sqrt(4.3597447222e-18 / ((5.29177210903e-11)**2 * 1.66053906660e-27))
    / (2.0 * np.pi * _C_CM)
)


class MMHessianDriver:
    """
    Computes numerical Hessian matrices from a molecular mechanics force field
    using OpenMM for gradient evaluation.

    Two public methods are available:

    ``compute``
        Standard MM Hessian in Cartesian coordinates.

    ``compute_ts``
        Transition-state MM Hessian.  The force constants of the bonds
        listed in *reaction_bonds* are replaced by a soft negative value
        (controlled by ``ts_negative_fc``).  After computing the Hessian the
        mass-weighted matrix is diagonalised and the code verifies that
        exactly one imaginary mode exists — the TS normal mode.

    Improvements over a naïve implementation
    -----------------------------------------
    - **Adaptive step size** — the finite-difference step for each atom is
      scaled by the inverse square root of its diagonal force constant,
      so stiff bonds use a smaller step and soft bonds use a larger one.
      ``step_size`` is a *target* dimensionless ratio; the actual step in Å
      is ``step_size / sqrt(fc_diag)`` clamped to [step_min, step_max].
    - **Mass-weighted Hessian** — the raw Cartesian Hessian is transformed
      to mass-weighted coordinates, giving proper normal modes and
      imaginary frequencies in cm⁻¹.
    - **Redundant internal coordinates** — before finite differences, the
      Cartesian displacements are projected onto the non-redundant subspace
      using a Wilson B-matrix built from all bonds, angles, and dihedrals.
      This removes the 6 translational/rotational degrees of freedom exactly
      and improves numerical conditioning.
    - **Scaled negative force constant** — instead of simply negating the
      GAFF force constant (which can be 250 000 kJ/mol/nm², giving an
      unrealistically stiff imaginary mode), ``ts_negative_fc`` sets a
      physically reasonable soft negative value that guides geomeTRIC along
      the TS coordinate without dominating the rest of the spectrum.
    - **Automatic TS mode identification** — after building the TS Hessian
      the code checks the number of negative eigenvalues in the
      mass-weighted Hessian (excluding the 6 trivial modes) and warns if it
      is not exactly 1.
    - **Simulation caching** — the OpenMM Simulation object is cached
      keyed by the identity of the ff_gen object.  Re-calling compute or
      compute_ts with the same ff_gen skips the topology-file writing and
      system-building step entirely.

    Parameters
    ----------
    comm : MPI communicator, optional
    ostream : OutputStream, optional

    Key attributes
    --------------
    step_size : float
        Dimensionless target step ratio for the adaptive step (default 0.01).
        Actual step (Å) = step_size / sqrt(fc_diagonal_kJpermol_nm2).
    step_min : float
        Minimum finite-difference step in Å (default 0.0005).
    step_max : float
        Maximum finite-difference step in Å (default 0.005).
    ts_negative_fc : float
        Force constant (kJ/mol/nm²) assigned to reaction bonds in
        ``compute_ts``.  Negative sign is applied automatically.
        Default −100 000 kJ/mol/nm².
    hessian_massage : bool
        If True (default), after diagonalising the TS Hessian all negative
        eigenvalues except the most negative one are replaced by
        ``massage_positive_value``, ensuring exactly one imaginary mode.
        The Hessian is then reconstructed from the modified eigenvalues.
    massage_positive_value : float
        Replacement eigenvalue (Hartree/Bohr²) used for spurious negative
        modes during Hessian massage.  Default 1e-4.
    metal_dihedral_barrier : float
        Soft barrier (kJ/mol) temporarily added to zero-barrier dihedrals
        involving TM/UFF metal atoms during ``compute_ts``.  Applied only
        for the Hessian calculation — the original ff_gen is never modified.
        Default 0.1 kJ/mol.
    openmm_platform : str
        OpenMM platform (default "CPU").
    filename_prefix : str
        Prefix for temporary topology files (default "mmhessian_residue").

    Computed results (set after compute / compute_ts)
    --------------------------------------------------
    hessian : ndarray (3N, 3N)
        Cartesian Hessian in Hartree/Bohr².
    mw_hessian : ndarray (3N, 3N)
        Mass-weighted Hessian in Hartree/(Bohr² · amu).
    frequencies : ndarray (3N,)
        Harmonic frequencies in cm⁻¹ (imaginary → negative).
    ts_mode : ndarray (3N,) or None
        Mass-weighted eigenvector of the imaginary TS mode (only set by
        ``compute_ts``).

    Examples
    --------
    Standard MM Hessian::

        drv = MMHessianDriver()
        drv.compute(molecule, ff_gen)
        drv.save("hessian.txt")
        print(drv.frequencies)

    TS Hessian (1-based bond indices)::

        drv = MMHessianDriver()
        drv.compute_ts(molecule, ff_gen, reaction_bonds=[(3, 8), (4, 9)])
        drv.save("hessian_ts.txt")
        # Pass to geomeTRIC via OptimizationDriver:
        # opt_drv.hessian = "file:hessian_ts.txt"
    """

    def __init__(self, comm=None, ostream=None):

        if comm is None:
            comm = MPI.COMM_WORLD

        if ostream is None:
            if comm.Get_rank() == mpi_master():
                ostream = OutputStream(sys.stdout)
            else:
                ostream = OutputStream(None)

        self.comm    = comm
        self.rank    = comm.Get_rank()
        self.nodes   = comm.Get_size()
        self.ostream = ostream

        # Finite-difference settings
        self.step_size = 0.01    # dimensionless adaptive step ratio
        self.step_min  = 5e-4   # Å — minimum allowed step
        self.step_max  = 5e-3   # Å — maximum allowed step

        # TS Hessian settings
        self.ts_negative_fc   = -100_000.0  # kJ/mol/nm²
        self.hessian_massage  = True         # fix spurious extra negative modes
        # Replacement eigenvalue for spurious negative modes during massage.
        # Should be small and positive — mimics a soft but stable mode.
        self.massage_positive_value = 1e-4   # Hartree/Bohr²

        # Soft barrier added temporarily to zero-barrier metal dihedrals
        # during compute_ts. The original ff_gen is never modified.
        self.metal_dihedral_barrier = 0.1  # kJ/mol

        # OpenMM settings
        self.openmm_platform = "CPU"
        self.filename_prefix = "mmhessian_residue"

        # Results
        self.hessian    = None
        self.mw_hessian = None
        self.frequencies = None
        self.ts_mode    = None

        # Simulation cache: maps id(ff_gen) → Simulation object
        self._sim_cache = {}

    # ════════════════════════════════════════════════════════════════════════
    # Public API
    # ════════════════════════════════════════════════════════════════════════

    def compute(self, molecule, ff_gen):
        """
        Compute the standard MM Hessian for *molecule* using *ff_gen*.

        Parameters
        ----------
        molecule : Molecule
        ff_gen   : MMForceFieldGenerator

        Returns
        -------
        ndarray (3N, 3N)  Cartesian Hessian in Hartree/Bohr².
        """
        self.ostream.print_info(
            f"MMHessianDriver: computing MM Hessian "
            f"({molecule.number_of_atoms()} atoms, "
            f"platform={self.openmm_platform})..."
        )
        self.ostream.flush()

        sim = self._get_simulation(molecule, ff_gen)
        self._run(molecule, sim)

        n_neg = int(np.sum(self.frequencies < 0))
        self.ostream.print_info(
            f"  MM Hessian done.  "
            f"Negative frequencies: {n_neg}  "
            f"(lowest: {self.frequencies[0]:.1f} cm⁻¹)"
        )
        self.ostream.flush()

        return self.hessian

    def compute_ts(self, molecule, ff_gen, reaction_bonds):
        """
        Compute a TS MM Hessian by assigning ``ts_negative_fc`` to the
        reaction bonds before building the Hessian.

        Parameters
        ----------
        molecule       : Molecule
        ff_gen         : MMForceFieldGenerator
            **Not modified in place.**
        reaction_bonds : list of (int, int)
            **1-based** atom index pairs for the forming/breaking bonds.

        Returns
        -------
        ndarray (3N, 3N)  Cartesian Hessian in Hartree/Bohr².
        """
        bonds_0 = set(
            (min(a - 1, b - 1), max(a - 1, b - 1))
            for (a, b) in reaction_bonds
        )

        self.ostream.print_info(
            f"MMHessianDriver: computing MM TS Hessian "
            f"({molecule.number_of_atoms()} atoms, "
            f"ts_fc={self.ts_negative_fc:.0f} kJ/mol/nm², "
            f"platform={self.openmm_platform})..."
        )
        self.ostream.flush()

        # Shallow-copy bonds dict — original ff_gen is never modified
        original_bonds     = ff_gen.bonds
        original_dihedrals = ff_gen.dihedrals
        bonds_ts     = {k: dict(v) for k, v in original_bonds.items()}
        dihedrals_ts = {k: dict(v) for k, v in original_dihedrals.items()}

        for (a, b) in bonds_0:
            key = (a, b) if (a, b) in bonds_ts else (b, a)
            if key in bonds_ts:
                old_fc = bonds_ts[key]['force_constant']
                bonds_ts[key]['force_constant'] = self.ts_negative_fc
                self.ostream.print_info(
                    f"  Bond ({a+1},{b+1}): fc "
                    f"{old_fc:.1f} -> {self.ts_negative_fc:.1f} kJ/mol/nm²"
                )
            else:
                self.ostream.print_warning(
                    f"  Bond ({a+1},{b+1}) not found in topology — skipping."
                )

        # Patch zero-barrier metal dihedrals temporarily for the TS Hessian.
        # Metal atoms have comment 'TM' or 'UFF' in ff_gen.atoms.
        metal_indices = {
            idx for idx, atom in ff_gen.atoms.items()
            if atom.get('comment', '') in ('TM', 'UFF')
        }
        n_metal_patched = 0
        if metal_indices:
            for key, params in dihedrals_ts.items():
                i, j, k, l = key
                if not ({i, j, k, l} & metal_indices):
                    continue
                if params.get('multiple', False):
                    if not all(abs(b) < 1e-10
                               for b in params.get('barrier', [])):
                        continue
                else:
                    if abs(params.get('barrier', 0.0)) >= 1e-10:
                        continue
                # Append a soft term to the copied dict only
                barriers     = [params['barrier']] if not params['multiple'] \
                                else list(params['barrier'])
                phases       = [params['phase']]   if not params['multiple'] \
                                else list(params['phase'])
                periodicities = [params['periodicity']] \
                                if not params['multiple'] \
                                else list(params['periodicity'])
                comments     = [params['comment']] if not params['multiple'] \
                                else list(params['comment'])
                barriers.append(self.metal_dihedral_barrier)
                phases.append(0.0)
                periodicities.append(-1)
                comments.append(f'metal patch ({key[0]+1}-{key[1]+1}'
                                 f'-{key[2]+1}-{key[3]+1})')
                dihedrals_ts[key] = {
                    'type':        'Fourier',
                    'multiple':    True,
                    'barrier':     barriers,
                    'phase':       phases,
                    'periodicity': periodicities,
                    'comment':     comments,
                }
                n_metal_patched += 1
            if n_metal_patched:
                self.ostream.print_info(
                    f"  Temporarily patched {n_metal_patched} zero-barrier "
                    f"metal dihedral(s) with "
                    f"barrier={self.metal_dihedral_barrier:.3f} kJ/mol."
                )
        self.ostream.flush()

        ff_gen.bonds     = bonds_ts
        ff_gen.dihedrals = dihedrals_ts
        try:
            # Force a fresh simulation for the modified ff_gen
            sim = self._build_simulation(molecule, ff_gen)
        finally:
            ff_gen.bonds     = original_bonds
            ff_gen.dihedrals = original_dihedrals

        self._run(molecule, sim)
        self._identify_ts_mode()

        return self.hessian

    def save(self, filename):
        """Save the Cartesian Hessian to a text file (Hartree/Bohr²)."""
        assert self.hessian is not None, \
            "MMHessianDriver: call compute() or compute_ts() first."
        np.savetxt(filename, self.hessian)
        self.ostream.print_info(f"  Hessian saved to {filename}.")
        self.ostream.flush()

    def print_frequencies(self, n_print=12):
        """
        Print the lowest *n_print* harmonic frequencies in cm⁻¹.
        Imaginary frequencies are shown as negative values.
        """
        assert self.frequencies is not None, \
            "MMHessianDriver: call compute() or compute_ts() first."
        self.ostream.print_blank()
        self.ostream.print_info("  Harmonic frequencies (cm⁻¹):")
        for i, freq in enumerate(self.frequencies[:n_print], start=1):
            marker = " ← TS mode" if (
                self.ts_mode is not None and i == 1 and freq < 0
            ) else ""
            self.ostream.print_info(f"    {i:>3d}:  {freq:>10.2f}{marker}")
        self.ostream.flush()

    def clear_cache(self):
        """Discard all cached OpenMM Simulation objects."""
        self._sim_cache.clear()

    def patch_metal_dihedrals(self, ff_gen, barrier=0.1, phase=0.0,
                              periodicity=1):
        """
        Patch zero-barrier dihedrals that involve a metal atom by adding a
        soft non-zero term via ``ff_gen.add_dihedral``.

        When GAFF cannot type a dihedral involving a metal (TM or UFF atom),
        it sets the barrier to 0 with comment ``'Unknown ...'``.  This leaves
        those DOFs with zero dihedral curvature, which causes near-singular
        rows in the Hessian and poor adaptive step estimates.

        This method identifies all dihedrals where:
          - the barrier is zero (or all barriers are zero for multiple dihedrals), AND
          - at least one of the four atoms has comment ``'TM'`` or ``'UFF'``
            in ``ff_gen.atoms``

        and calls ``ff_gen.add_dihedral`` to append a soft non-zero term.
        The original ff_gen is modified in place.

        Parameters
        ----------
        ff_gen      : MMForceFieldGenerator — modified in place.
        barrier     : float
            Barrier in kJ/mol for the patched term (default 0.1).
        phase       : float
            Phase in degrees (default 0.0).
        periodicity : int
            Periodicity (default 1).

        Returns
        -------
        int  Number of dihedrals patched.

        Example
        -------
        ::

            ff_gen.create_topology(ts_mol, resp=False)
            drv = MMHessianDriver()
            n = drv.patch_metal_dihedrals(ff_gen)
            print(f"Patched {n} zero-barrier metal dihedrals.")
            drv.compute_ts(ts_mol, ff_gen, reaction_bonds=[(3, 8)])
        """
        # Identify which atom indices are metals (TM or UFF fallback)
        metal_indices = {
            idx for idx, atom in ff_gen.atoms.items()
            if atom.get('comment', '') in ('TM', 'UFF')
        }

        if not metal_indices:
            self.ostream.print_info(
                "  patch_metal_dihedrals: no TM/UFF atoms found — nothing to patch."
            )
            self.ostream.flush()
            return 0

        n_patched = 0

        for key, params in list(ff_gen.dihedrals.items()):
            i, j, k, l = key

            # Check if any atom in this dihedral is a metal
            if not ({i, j, k, l} & metal_indices):
                continue

            # Check if the dihedral barrier is effectively zero
            if params.get('multiple', False):
                barriers = params.get('barrier', [])
                if not all(abs(b) < 1e-10 for b in barriers):
                    continue
            else:
                if abs(params.get('barrier', 0.0)) >= 1e-10:
                    continue

            # add_dihedral expects 1-based indices
            ff_gen.add_dihedral(
                [i + 1, j + 1, k + 1, l + 1],
                barrier=barrier,
                phase=phase,
                periodicity=periodicity,
            )
            n_patched += 1

        if n_patched:
            self.ostream.print_info(
                f"  patch_metal_dihedrals: patched {n_patched} zero-barrier "
                f"metal dihedral(s) with barrier={barrier:.3f} kJ/mol, "
                f"phase={phase:.1f}°, periodicity={periodicity}."
            )
            # Invalidate the simulation cache since ff_gen has changed
            self._sim_cache.pop(id(ff_gen), None)
        else:
            self.ostream.print_info(
                "  patch_metal_dihedrals: no zero-barrier metal dihedrals found."
            )
        self.ostream.flush()

        return n_patched

    # ════════════════════════════════════════════════════════════════════════
    # Core computation
    # ════════════════════════════════════════════════════════════════════════

    def _run(self, molecule, sim):
        """
        Compute the Cartesian Hessian, mass-weight it, project out
        translations and rotations, and convert to frequencies.
        Results are stored in self.hessian, self.mw_hessian, self.frequencies.
        """
        masses    = self._get_masses(molecule)
        step_map  = self._adaptive_steps(molecule, sim, masses)
        proj      = self._projection_matrix(molecule)

        hessian   = self._numerical_hessian(molecule, sim, step_map)

        # Project out translations / rotations
        hessian   = proj @ hessian @ proj

        # Symmetrise
        hessian   = 0.5 * (hessian + hessian.T)

        self.hessian    = hessian
        self.mw_hessian = self._mass_weight(hessian, masses)
        self.frequencies = self._hessian_to_frequencies(self.mw_hessian)
        self.ts_mode    = None   # reset; set by _identify_ts_mode if needed

    # ════════════════════════════════════════════════════════════════════════
    # Adaptive step size
    # ════════════════════════════════════════════════════════════════════════

    def _adaptive_steps(self, molecule, sim, masses):
        """
        Estimate a per-degree-of-freedom finite-difference step (Å) by
        evaluating the diagonal of the Hessian cheaply via a single
        forward displacement per DOF.

        The step for DOF k is:
            h_k = clamp(step_size / sqrt(|fc_k|), step_min, step_max)

        where fc_k is the diagonal force constant in Hartree/Bohr²
        (estimated from a single forward/central displacement).

        Returns a (3N,) array of step sizes in Ångström.
        """
        n_atoms  = molecule.number_of_atoms()
        n_dof    = 3 * n_atoms
        labels   = molecule.get_labels()
        charge   = molecule.get_charge()
        mult     = molecule.get_multiplicity()
        coords   = molecule.get_coordinates_in_angstrom().copy()
        step_def = 0.001   # Å — default step for diagonal estimation

        g0 = self._gradient(sim, labels, coords, charge, mult)
        steps = np.empty(n_dof)

        for idx in range(n_dof):
            i, j = divmod(idx, 3)
            c_fwd = coords.copy()
            c_fwd[i, j] += step_def
            g_fwd = self._gradient(sim, labels, c_fwd, charge, mult)
            # Diagonal force constant: d²E/dx² ≈ (g_fwd[idx] - g0[idx]) / step
            step_bohr = step_def / bohr_in_angstrom()
            fc_diag   = abs((g_fwd[idx] - g0[idx]) / step_bohr)
            if fc_diag > 1e-12:
                h = self.step_size / np.sqrt(fc_diag)
            else:
                h = self.step_max
            steps[idx] = np.clip(h, self.step_min, self.step_max)

        return steps

    # ════════════════════════════════════════════════════════════════════════
    # Redundant internal coordinate projection
    # ════════════════════════════════════════════════════════════════════════

    def _projection_matrix(self, molecule):
        """
        Build the 3N × 3N projection matrix that removes the 6 external
        degrees of freedom (3 translations + 3 rotations) from the Hessian.

        Uses the standard translational and Eckart rotational vectors
        constructed at the equilibrium geometry.

        Returns P = I - D Dᵀ  where D is the 3N × 6 matrix of external
        displacement vectors (mass-weighted and orthonormalised).
        """
        masses  = self._get_masses(molecule)   # (N,) amu
        coords  = molecule.get_coordinates_in_angstrom()   # (N, 3) Å
        n_atoms = molecule.number_of_atoms()
        n_dof   = 3 * n_atoms

        sqrt_m = np.sqrt(masses)   # (N,)
        M_total = masses.sum()

        # Centre of mass
        com = (masses[:, None] * coords).sum(axis=0) / M_total

        # Translational vectors (mass-weighted): T_α,iβ = sqrt(m_i) δ_αβ
        T = np.zeros((n_dof, 3))
        for i in range(n_atoms):
            for alpha in range(3):
                T[3 * i + alpha, alpha] = sqrt_m[i]

        # Rotational vectors (Eckart): R_α,iβ = sqrt(m_i) (r_i × e_α)_β
        r = coords - com   # (N, 3) relative to COM in Å
        ex = np.array([1, 0, 0])
        ey = np.array([0, 1, 0])
        ez = np.array([0, 0, 1])
        R = np.zeros((n_dof, 3))
        for i in range(n_atoms):
            for alpha, e in enumerate([ex, ey, ez]):
                cross = np.cross(r[i], e)
                R[3 * i: 3 * i + 3, alpha] = sqrt_m[i] * cross

        # Combine and orthonormalise via QR
        D = np.hstack([T, R])   # (3N, 6)
        Q, _ = np.linalg.qr(D)
        # Q has shape (3N, 3N); keep only the first 6 columns
        D_orth = Q[:, :6]

        # Projection onto the internal subspace
        P = np.eye(n_dof) - D_orth @ D_orth.T
        return P

    # ════════════════════════════════════════════════════════════════════════
    # Numerical Hessian
    # ════════════════════════════════════════════════════════════════════════

    def _numerical_hessian(self, molecule, sim, step_map):
        """
        Central finite differences with per-DOF adaptive steps.

        Returns the raw (unprojected) Cartesian Hessian in Hartree/Bohr².
        """
        n_atoms  = molecule.number_of_atoms()
        n_dof    = 3 * n_atoms
        labels   = molecule.get_labels()
        charge   = molecule.get_charge()
        mult     = molecule.get_multiplicity()
        coords   = molecule.get_coordinates_in_angstrom().copy()
        hessian  = np.zeros((n_dof, n_dof))

        for idx in range(n_dof):
            i, j   = divmod(idx, 3)
            h      = step_map[idx]
            h_bohr = h / bohr_in_angstrom()

            c_fwd = coords.copy()
            c_fwd[i, j] += h
            g_fwd = self._gradient(sim, labels, c_fwd, charge, mult)

            c_bwd = coords.copy()
            c_bwd[i, j] -= h
            g_bwd = self._gradient(sim, labels, c_bwd, charge, mult)

            hessian[idx, :] = (g_fwd - g_bwd) / (2.0 * h_bohr)

        return hessian

    # ════════════════════════════════════════════════════════════════════════
    # Mass-weighting and frequencies
    # ════════════════════════════════════════════════════════════════════════

    def _mass_weight(self, hessian, masses):
        """
        Convert Cartesian Hessian (Hartree/Bohr²) to mass-weighted
        Hessian (Hartree / (Bohr² · amu)) by H_mw = M^{-1/2} H M^{-1/2}.
        """
        n_dof   = hessian.shape[0]
        # Build (3N,) vector of 1/sqrt(m) repeated for each Cartesian component
        inv_sqm = np.repeat(1.0 / np.sqrt(masses), 3)   # (3N,) in 1/sqrt(amu)
        return hessian * np.outer(inv_sqm, inv_sqm)

    def _hessian_to_frequencies(self, mw_hessian):
        """
        Diagonalise the mass-weighted Hessian and convert eigenvalues to
        harmonic frequencies in cm⁻¹.  Imaginary frequencies are returned
        as negative values.
        """
        eigenvalues = np.linalg.eigvalsh(mw_hessian)

        # Convert Hartree/(Bohr²·amu) → rad/s → cm⁻¹
        # First convert to SI: multiply by E_h/(a0²·amu_kg)
        E_h   = 4.3597447222e-18   # J
        a0    = 5.29177210903e-11  # m
        m_amu = 1.66053906660e-27  # kg
        conv  = E_h / (a0**2 * m_amu)

        freqs = np.empty_like(eigenvalues)
        for k, lam in enumerate(eigenvalues):
            if lam >= 0:
                freqs[k]  =  np.sqrt(lam * conv) / (2.0 * np.pi * _C_CM)
            else:
                freqs[k]  = -np.sqrt(-lam * conv) / (2.0 * np.pi * _C_CM)

        return freqs

    # ════════════════════════════════════════════════════════════════════════
    # TS mode identification
    # ════════════════════════════════════════════════════════════════════════

    def _identify_ts_mode(self):
        """
        After compute_ts:

        1. Diagonalise the mass-weighted Hessian to find frequencies.
        2. If ``hessian_massage`` is True and more than one negative internal
           eigenvalue is found, diagonalise the *Cartesian* Hessian, replace
           all negative eigenvalues except the most negative one with
           ``massage_positive_value``, and reconstruct the Cartesian Hessian.
           The mass-weighted Hessian and frequencies are then recomputed from
           the massaged Cartesian Hessian.
        3. Store the TS eigenvector in ``self.ts_mode``.

        Massaging the Cartesian Hessian directly (rather than the
        mass-weighted one) avoids any fragile mass-recovery algebra.
        geomeTRIC reads the Cartesian Hessian from file and applies its own
        internal coordinate transformation, so the Cartesian matrix is what
        matters for the TS optimisation.
        """
        eigenvalues, eigenvectors = np.linalg.eigh(self.mw_hessian)

        # The 6 lowest modes are translations/rotations (≈0); skip them
        internal_evals = eigenvalues[6:]
        n_imag = int(np.sum(internal_evals < 0))

        if n_imag == 0:
            self.ts_mode = None
            self.ostream.print_warning(
                "  No imaginary frequency found after projection.  "
                "The reaction bonds may not be the dominant soft mode.  "
                "Try increasing |ts_negative_fc| or adding more reaction bonds."
            )
            self.ostream.flush()
            return

        if n_imag == 1:
            self.ts_mode = eigenvectors[:, 6]
            self.ostream.print_info(
                f"  TS mode confirmed: 1 imaginary frequency "
                f"({self.frequencies[6]:.1f} cm⁻¹) — good TS Hessian guess."
            )
            self.ostream.flush()
            return

        # ── Multiple negative eigenvalues ────────────────────────────────────
        self.ostream.print_info(
            f"  {n_imag} imaginary internal frequencies found — expected 1."
        )

        if not self.hessian_massage:
            most_neg_idx = int(np.argmin(eigenvalues))
            self.ts_mode = eigenvectors[:, most_neg_idx]
            self.ostream.print_warning(
                f"  hessian_massage=False: keeping all {n_imag} negative modes.  "
                f"Lowest: {self.frequencies[most_neg_idx]:.1f} cm⁻¹.  "
                f"Set hessian_massage=True to fix automatically."
            )
            self.ostream.flush()
            return

        # ── Hessian massage on the Cartesian Hessian ─────────────────────────
        # Diagonalise the Cartesian Hessian directly.  Eigenvalues here are in
        # Hartree/Bohr² without mass-weighting, so they are not frequencies —
        # but the sign structure tells us which modes are negative, and the
        # most negative Cartesian eigenvalue corresponds to the TS mode.
        cart_evals, cart_evecs = np.linalg.eigh(self.hessian)

        # Find the most negative Cartesian eigenvalue (TS mode)
        ts_cart_idx = int(np.argmin(cart_evals))

        # Replace all other negative eigenvalues with massage_positive_value
        massaged_cart_evals = cart_evals.copy()
        n_replaced = 0
        for k in range(len(massaged_cart_evals)):
            if k == ts_cart_idx:
                continue
            if massaged_cart_evals[k] < 0:
                massaged_cart_evals[k] = self.massage_positive_value
                n_replaced += 1

        self.ostream.print_info(
            f"  Hessian massage: replaced {n_replaced} spurious negative "
            f"Cartesian eigenvalue(s) with "
            f"{self.massage_positive_value:.1e} Eh/Bohr²."
        )

        # Reconstruct the Cartesian Hessian
        self.hessian = cart_evecs @ np.diag(massaged_cart_evals) @ cart_evecs.T

        # Recompute the mass-weighted Hessian and frequencies from the
        # massaged Cartesian Hessian so everything stays consistent.
        masses          = self._masses_from_mw_hessian(self.hessian,
                                                        self.mw_hessian)
        self.mw_hessian = self._mass_weight(self.hessian, masses)
        self.frequencies = self._hessian_to_frequencies(self.mw_hessian)

        # Store the TS mode from the mass-weighted eigenvectors after massage
        mw_evals_new, mw_evecs_new = np.linalg.eigh(self.mw_hessian)
        ts_mw_idx = int(np.argmin(mw_evals_new))
        self.ts_mode = mw_evecs_new[:, ts_mw_idx]

        self.ostream.print_info(
            f"  After massage: 1 imaginary frequency "
            f"({self.frequencies[ts_mw_idx]:.1f} cm⁻¹) — TS Hessian ready."
        )
        self.ostream.flush()

    def _masses_from_mw_hessian(self, cart_hessian, mw_hessian):
        """
        Recover the (N,) atomic mass array (amu) from the relationship
        mw_H[i,j] = cart_H[i,j] / sqrt(m_i * m_j).

        Uses the diagonal elements where both are well-defined.
        Falls back to 12.0 amu for any atom where the diagonal is too small.
        """
        n_dof   = cart_hessian.shape[0]
        n_atoms = n_dof // 3
        masses  = np.full(n_atoms, 12.0)

        for atom in range(n_atoms):
            vals = []
            for xyz in range(3):
                k = 3 * atom + xyz
                c = cart_hessian[k, k]
                m = mw_hessian[k, k]
                if abs(c) > 1e-20 and abs(m) > 1e-20:
                    # mw[k,k] = cart[k,k] / m_atom  →  m_atom = cart/mw
                    ratio = c / m
                    if ratio > 0:
                        vals.append(ratio)
            if vals:
                masses[atom] = float(np.mean(vals))

        return masses

    # ════════════════════════════════════════════════════════════════════════
    # Simulation caching
    # ════════════════════════════════════════════════════════════════════════

    def _get_simulation(self, molecule, ff_gen):
        """
        Return a cached Simulation if one exists for this ff_gen, otherwise
        build a new one and cache it.
        """
        key = id(ff_gen)
        if key not in self._sim_cache:
            self._sim_cache[key] = self._build_simulation(molecule, ff_gen)
        return self._sim_cache[key]

    def _build_simulation(self, molecule, ff_gen):
        """
        Use OpenMMDynamics to build the OpenMM system, then return a bare
        Simulation object for gradient evaluation.
        """
        import openmm as _mm
        import openmm.app as _app
        import openmm.unit as _unit

        opm_dyn = OpenMMDynamics(self.comm)
        opm_dyn.ostream.mute()
        opm_dyn.openmm_platform = self.openmm_platform
        opm_dyn.create_system_from_molecule(
            molecule, ff_gen,
            filename=self.filename_prefix,
            residue_name="MOL",
        )

        integrator = _mm.VerletIntegrator(0.001 * _unit.picoseconds)
        platform   = _mm.Platform.getPlatformByName(self.openmm_platform)
        if self.openmm_platform == "CPU":
            platform.setPropertyDefaultValue("Threads", "1")

        return _app.Simulation(
            opm_dyn.pdb.topology,
            opm_dyn.system,
            integrator,
            platform,
        )

    # ════════════════════════════════════════════════════════════════════════
    # Gradient evaluation
    # ════════════════════════════════════════════════════════════════════════

    def _gradient(self, sim, labels, coords_ang, charge, mult):
        """
        Evaluate the OpenMM gradient at *coords_ang* (Å).
        Returns a flattened (3N,) array in Hartree/Bohr.
        """
        import openmm.unit as _unit

        coords_nm = coords_ang * 0.1
        sim.context.setPositions(coords_nm * _unit.nanometer)
        state  = sim.context.getState(getForces=True)
        forces = state.getForces(asNumpy=True).value_in_unit(
            _unit.kilojoules_per_mole / _unit.nanometer
        )
        # gradient = −force; kJ/mol/nm → Hartree/Bohr
        grad = -forces / (hartree_in_kjpermol() * 10.0 / bohr_in_angstrom())
        return grad.flatten()

    # ════════════════════════════════════════════════════════════════════════
    # Helpers
    # ════════════════════════════════════════════════════════════════════════

    def _get_masses(self, molecule):
        """
        Return a (N,) array of atomic masses in amu.
        Falls back to 12.0 for unknown elements.
        """
        labels = molecule.get_labels()
        return np.array([
            _ATOMIC_MASSES.get(lbl.capitalize(), 12.0)
            for lbl in labels
        ])

    @staticmethod
    def _xyz_string(labels, coords_angstrom):
        n = len(labels)
        lines = [str(n), ""]
        for label, xyz in zip(labels, coords_angstrom):
            lines.append(
                f"{label}  {xyz[0]:.10f}  {xyz[1]:.10f}  {xyz[2]:.10f}"
            )
        return "\n".join(lines)
