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


class MMHessianDriver:
    """
    Computes numerical Hessian matrices from a molecular mechanics force field
    using OpenMM for gradient evaluation.

    Two modes are available:

    ``compute``
        Standard MM Hessian.  All force constants are used as-is.

    ``compute_ts``
        Transition-state MM Hessian.  The force constants of the bonds
        specified in *reaction_bonds* are negated before computing the
        Hessian.  This produces a negative eigenvalue along the concerted
        bond-stretch mode, giving geometry optimisers (e.g. geomeTRIC) a
        proper TS Hessian guess at essentially zero cost.

    The Hessian is computed by central finite differences of the OpenMM
    analytical gradients.  Units are Hartree/Bohr² (consistent with
    VeloxChem's other Hessian drivers).

    Parameters
    ----------
    comm : MPI communicator, optional
    ostream : OutputStream, optional

    Attributes
    ----------
    hessian : ndarray, shape (3N, 3N)
        The most recently computed Hessian matrix.
    step_size : float
        Finite-difference step in Ångström (default 0.001).
    openmm_platform : str
        OpenMM platform name (default "CPU").
        Set to "CUDA" or "OpenCL" for GPU acceleration.
    filename_prefix : str
        Prefix for the temporary topology files written by OpenMMDynamics
        (default "mmhessian_residue").

    Examples
    --------
    Standard MM Hessian::

        ff_gen = MMForceFieldGenerator()
        ff_gen.create_topology(molecule)

        drv = MMHessianDriver()
        drv.compute(molecule, ff_gen)
        np.savetxt("hessian.txt", drv.hessian)

    TS Hessian (1-based bond indices)::

        drv = MMHessianDriver()
        drv.compute_ts(
            molecule, ff_gen,
            reaction_bonds=[(3, 8), (4, 9)],   # 1-based
        )
        np.savetxt("hessian.txt", drv.hessian)
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

        self.hessian          = None
        self.step_size        = 0.001   # Ångström
        self.openmm_platform  = "CPU"
        self.filename_prefix  = "mmhessian_residue"

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
            Fully built force field (``create_topology`` already called).

        Returns
        -------
        ndarray, shape (3N, 3N)
            Hessian in Hartree/Bohr².  Also stored in ``self.hessian``.
        """
        self.ostream.print_info(
            f"MMHessianDriver: computing MM Hessian "
            f"({molecule.number_of_atoms()} atoms, "
            f"step={self.step_size} Å, platform={self.openmm_platform})..."
        )
        self.ostream.flush()

        sim = self._build_simulation(molecule, ff_gen)
        self.hessian = self._numerical_hessian(molecule, sim)

        n_neg = int(np.sum(np.linalg.eigvalsh(self.hessian) < 0))
        self.ostream.print_info(
            f"  MM Hessian done.  Negative eigenvalues: {n_neg}"
        )
        self.ostream.flush()

        return self.hessian

    def compute_ts(self, molecule, ff_gen, reaction_bonds):
        """
        Compute a transition-state MM Hessian by negating the force constants
        of the bonds listed in *reaction_bonds* before the Hessian calculation.

        For a concerted reaction (e.g. Diels-Alder) the coupled negative
        springs combine into a single negative normal mode, giving geomeTRIC
        a correct TS Hessian guess at near-zero cost.

        Parameters
        ----------
        molecule       : Molecule
        ff_gen         : MMForceFieldGenerator
            Fully built force field.  **Not modified in place** — a shallow
            copy of the bonds dict is used temporarily.
        reaction_bonds : list of (int, int)
            **1-based** atom index pairs for the bonds whose force constants
            should be negated.  These are typically the forming and breaking
            bonds of the reaction step.

        Returns
        -------
        ndarray, shape (3N, 3N)
            TS Hessian in Hartree/Bohr².  Also stored in ``self.hessian``.

        Examples
        --------
        ::

            drv.compute_ts(ts_mol, ff_gen,
                           reaction_bonds=[(3, 8), (4, 9)])
        """
        # Convert to 0-based
        bonds_0 = set(
            (min(a - 1, b - 1), max(a - 1, b - 1))
            for (a, b) in reaction_bonds
        )

        self.ostream.print_info(
            f"MMHessianDriver: computing MM TS Hessian "
            f"({molecule.number_of_atoms()} atoms, "
            f"step={self.step_size} Å, platform={self.openmm_platform})..."
        )
        self.ostream.flush()

        # Shallow-copy the bonds dict so the original ff_gen is untouched
        original_bonds = ff_gen.bonds
        bonds_ts = {k: dict(v) for k, v in original_bonds.items()}

        for (a, b) in bonds_0:
            key = (a, b) if (a, b) in bonds_ts else (b, a)
            if key in bonds_ts:
                fc = bonds_ts[key]['force_constant']
                bonds_ts[key]['force_constant'] = -abs(fc)
                self.ostream.print_info(
                    f"  Bond ({a+1},{b+1}): fc {fc:.1f} -> {-abs(fc):.1f} kJ/mol/nm²"
                )
            else:
                self.ostream.print_warning(
                    f"  Bond ({a+1},{b+1}) not found in topology — skipping."
                )
        self.ostream.flush()

        # Temporarily swap in the modified bonds dict
        ff_gen.bonds = bonds_ts
        try:
            sim = self._build_simulation(molecule, ff_gen)
        finally:
            ff_gen.bonds = original_bonds   # always restore

        self.hessian = self._numerical_hessian(molecule, sim)

        eigenvalues = np.linalg.eigvalsh(self.hessian)
        n_neg = int(np.sum(eigenvalues < 0))
        self.ostream.print_info(
            f"  MM TS Hessian done.  Negative eigenvalues: {n_neg}  "
            f"(smallest: {eigenvalues[0]:.6f} Eh/Bohr²)"
        )
        self.ostream.flush()

        return self.hessian

    def save(self, filename):
        """
        Save the most recently computed Hessian to a text file.

        Parameters
        ----------
        filename : str
        """
        assert self.hessian is not None, \
            "MMHessianDriver: no Hessian computed yet — call compute() or compute_ts() first."
        np.savetxt(filename, self.hessian)
        self.ostream.print_info(f"  Hessian saved to {filename}.")
        self.ostream.flush()

    # ════════════════════════════════════════════════════════════════════════
    # Private helpers
    # ════════════════════════════════════════════════════════════════════════

    def _build_simulation(self, molecule, ff_gen):
        """
        Use OpenMMDynamics to build the OpenMM system from *ff_gen*, then
        return a bare Simulation object suitable for gradient evaluation.
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

        sim = _app.Simulation(
            opm_dyn.pdb.topology,
            opm_dyn.system,
            integrator,
            platform,
        )
        return sim

    def _get_gradient(self, sim, molecule):
        """
        Evaluate the OpenMM gradient for *molecule* using *sim*.

        Returns a flattened (3N,) array in Hartree/Bohr.
        """
        import openmm.unit as _unit

        coords_nm = molecule.get_coordinates_in_angstrom() * 0.1
        sim.context.setPositions(coords_nm * _unit.nanometer)
        state  = sim.context.getState(getForces=True)
        forces = state.getForces(asNumpy=True).value_in_unit(
            _unit.kilojoules_per_mole / _unit.nanometer
        )
        # gradient = -force; convert kJ/mol/nm → Hartree/Bohr
        grad = -forces / (hartree_in_kjpermol() * 10.0 / bohr_in_angstrom())
        return grad.flatten()

    def _numerical_hessian(self, molecule, sim):
        """
        Compute the 3N × 3N Hessian by central finite differences of the
        OpenMM gradients.  Step size is ``self.step_size`` Ångström.

        Returns a symmetrised Hessian in Hartree/Bohr².
        """
        n_atoms   = molecule.number_of_atoms()
        n_dof     = 3 * n_atoms
        step_ang  = self.step_size
        step_bohr = step_ang / bohr_in_angstrom()
        labels    = molecule.get_labels()
        charge    = molecule.get_charge()
        mult      = molecule.get_multiplicity()
        coords    = molecule.get_coordinates_in_angstrom().copy()
        hessian   = np.zeros((n_dof, n_dof))

        for i in range(n_atoms):
            for j in range(3):
                idx = 3 * i + j

                # Forward displacement
                c_fwd = coords.copy()
                c_fwd[i, j] += step_ang
                mol_fwd = Molecule.read_xyz_string(
                    self._xyz_string(labels, c_fwd)
                )
                mol_fwd.set_charge(charge)
                mol_fwd.set_multiplicity(mult)
                g_fwd = self._get_gradient(sim, mol_fwd)

                # Backward displacement
                c_bwd = coords.copy()
                c_bwd[i, j] -= step_ang
                mol_bwd = Molecule.read_xyz_string(
                    self._xyz_string(labels, c_bwd)
                )
                mol_bwd.set_charge(charge)
                mol_bwd.set_multiplicity(mult)
                g_bwd = self._get_gradient(sim, mol_bwd)

                hessian[idx, :] = (g_fwd - g_bwd) / (2.0 * step_bohr)

        # Symmetrise to remove small numerical asymmetries
        return 0.5 * (hessian + hessian.T)

    @staticmethod
    def _xyz_string(labels, coords_angstrom):
        """Build a minimal XYZ string from labels and coordinates."""
        n = len(labels)
        lines = [str(n), ""]
        for label, xyz in zip(labels, coords_angstrom):
            lines.append(
                f"{label}  {xyz[0]:.10f}  {xyz[1]:.10f}  {xyz[2]:.10f}"
            )
        return "\n".join(lines)
