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
"""Calculate aqueous CDAM properties for one parameter point.

``Calculator.calculate`` plans the jobs (``AQUEOUS_PROPERTIES`` in
``models.py``), prepares a private copy of the solvated system, then runs:

* ``md``  -> IOD, CN, geometry RMSDs   (``md.py`` run, ``geometry.py`` analysis)
* ``fep`` -> HFE                        (``md.py`` FEP; equilibrated by MD first)
* ``qm``  -> QM_ENERGY                  (VeloxChem single point, below)

and returns the values in the units the targets ask for, with metadata
(method, options, spread, uncertainty) and warnings. Comparing with
references is the ``Evaluator``'s job. The first section holds the cached
VeloxChem single-point SCF, keyed by geometry, state and settings.
"""

from __future__ import annotations

import json
import math
import sys
from typing import Any, Iterable
import time
from dataclasses import dataclass, field, replace
from pathlib import Path

import numpy as np
try:
    import openmm as mm
    import openmm.unit as unit
except ImportError:
    pass

from mpi4py import MPI

from ..errorhandler import assert_msg_critical
from ..molecularbasis import MolecularBasis
from ..molecule import Molecule
from ..outputstream import OutputStream
from ..scfrestdriver import ScfRestrictedDriver
from ..scfunrestdriver import ScfUnrestrictedDriver
from ..veloxchemlib import chemical_element_identifier
from .candidate import CDAMCandidate
from .generator import CoordinateGenerator
from .geometry import analyze_trajectory, count_waters
from .md import (check_fep_backend, fep_uncertainty, finite_size_correction,
                 run_fep, default_platform, md_settings, run_md)
from .models import (print_settings, CalculationResult, FitTarget,
                     GeneratedBatch, OpenMMSystems, ParameterSet,
                     default_ostream, aqueous_property, convert_units)

# ==========================================================================
# qm
# ==========================================================================

DEFAULT_BASIS = "def2-SVP"
SINGLE_POINT_SETTINGS = frozenset({
    "basis",
    "xcfun",
    "conv_thresh",
    "max_iter",
    "level_shifting",
    "pfon",
    "ri_coulomb",
    "solvation_model",
})


def real_atom_molecule(batch, charge: float, multiplicity: int) -> Molecule:
    """The first configuration of ``batch`` without its CDAM dummy sites."""
    atoms = [(atom["element"], [float(v)
                                for v in xyz])
             for atom, xyz in zip(batch.system.atoms, batch.configurations[0])
             if atom.get("role") != "dummy" and
             chemical_element_identifier(atom["element"].upper()) > 0]
    if not atoms:
        raise ValueError(
            "QM single-point calculation requires at least one real atom")
    labels, coordinates = zip(*atoms)
    molecule = Molecule(list(labels), list(coordinates), "angstrom")
    molecule.set_charge(charge)
    molecule.set_multiplicity(multiplicity)
    return molecule


def print_setup(ostream, molecule: Molecule, settings: dict[str, Any]) -> None:
    """Print the QM single-point setup (atoms, charge, multiplicity, basis, SCF options)."""
    multiplicity = molecule.get_multiplicity()
    print_settings(ostream, "QM Single-Point Setup", [
        ("Atoms (dummy sites excluded)", " ".join(molecule.get_labels())),
        ("Molecular Charge", molecule.get_charge()),
        ("Spin Multiplicity", multiplicity),
        ("SCF Type", "restricted" if multiplicity == 1 else "unrestricted"),
        ("Basis Set", settings["basis"]),
        *((name.replace("_", " ").title(), value)
          for name, value in sorted(settings.items())
          if name != "basis"),
    ])


class SinglePointCalculator:
    """Converged SCF energies (Hartree), restricted for singlets, unrestricted otherwise."""

    def __init__(self, comm=None, ostream=None):
        """Store the communicator and output stream; start an empty energy cache."""
        self.comm = comm if comm is not None else MPI.COMM_WORLD
        self.ostream = ostream
        self.cache: dict[tuple, float] = {}

    def energy(self,
               molecule: Molecule,
               settings: dict[str, Any] | None = None,
               metadata: dict[str, Any] | None = None) -> float:
        """SCF energy; ``metadata`` (optional) receives ``cache_hit`` and ``scf_converged``."""
        settings = dict(settings or {})
        unknown = set(settings) - SINGLE_POINT_SETTINGS
        if unknown:
            raise ValueError(
                f"unsupported single-point settings: {sorted(unknown)}")
        settings.setdefault("basis", DEFAULT_BASIS)
        metadata = metadata if metadata is not None else {}
        multiplicity = molecule.get_multiplicity()
        key = (tuple(molecule.get_labels()),
               tuple(molecule.get_atom_basis_labels()),
               tuple(map(tuple, molecule.get_coordinates_in_angstrom())),
               float(molecule.get_charge()), multiplicity,
               json.dumps(settings, sort_keys=True, allow_nan=False))
        metadata["cache_hit"] = metadata["scf_converged"] = key in self.cache
        if key in self.cache:
            self.ostream.print_info("Using cached QM single-point energy.")
            return self.cache[key]

        driver = (ScfRestrictedDriver if multiplicity == 1 else
                  ScfUnrestrictedDriver)(comm=self.comm, ostream=self.ostream)
        for name, value in settings.items():
            if name != "basis":
                setattr(driver, name, value)
        basis = MolecularBasis.read(molecule, settings["basis"])
        self.ostream.print_info(
            f"Starting QM single-point SCF with {settings['basis']}...")
        self.ostream.flush()
        energy = float(driver.compute(molecule, basis)["scf_energy"])
        if not driver.is_converged or not math.isfinite(energy):
            raise RuntimeError(
                "QM single-point SCF did not converge to a finite energy")
        self.cache[key] = energy
        metadata["scf_converged"] = True
        return energy


# ==========================================================================
# calculator
# ==========================================================================

DEFAULT_ION_RADIUS_NM = 0.20  # HFE correction radius when no MD measured R1
SPREAD_DEFINITIONS = {
    "IOD": "shell_distance_std_angstrom",
    "CN": "coordination_number_frame_std"
}


@dataclass
class CalculationPlan:
    """
    Requested properties grouped by the job that computes them.

    :param jobs:
        The target names by job (md, fep, qm).
    :param methods:
        The method label by target name.
    """

    jobs: dict[str, list[str]] = field(default_factory=dict)
    methods: dict[str, str] = field(default_factory=dict)


@dataclass
class _Results:
    """Values, per-target metadata and warnings collected during a calculation."""

    targets: list[FitTarget]
    plan: CalculationPlan
    values: dict[str, Any] = field(default_factory=dict)
    metadata: dict[str, dict] = field(init=False)
    warnings: list[str] = field(default_factory=list)

    def __post_init__(self):
        """Start each target's metadata with the method chosen by the plan."""
        self.metadata = {
            t.name: {
                "method": self.plan.methods[t.name]
            } for t in self.targets
        }

    def finish(self, parameters: ParameterSet, **raw) -> CalculationResult:
        """Convert values to the target units in place and wrap everything in a ``CalculationResult``.

        Extra keyword arguments are stored in ``raw_data`` next to the parameters.
        """
        convert_units(self.values, self.metadata, self.targets)
        return CalculationResult(values=self.values,
                                 metadata=self.metadata,
                                 raw_data={
                                     "parameters": parameters.to_dict(),
                                     **raw
                                 },
                                 warnings=self.warnings)


class Calculator:
    """
    Runs the physical property calculations for one generated parameter point.

    Keyword options of ``calculate``: ``md_options``, ``fep_options``,
    ``qm_options`` and ``work_dir``. In-memory systems are copied first, so the
    caller's batch is never modified.

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.
    """

    requires_openmm_systems = True

    def __init__(self, comm=None, ostream: OutputStream | None = None) -> None:
        """
        Creates the calculator with an output stream and an empty QM single-point
        cache; OpenMM is required.
        """

        assert_msg_critical("openmm" in sys.modules,
                            f"{type(self).__name__}: OpenMM is required")
        self.comm = comm if comm is not None else MPI.COMM_WORLD
        self.ostream = default_ostream(ostream, self.comm)
        self.qm = SinglePointCalculator(self.comm, self.ostream)

    # ------------------------------------------------------------ public API

    def plan(self, targets: Iterable[FitTarget]) -> CalculationPlan:
        """
        Groups targets by job (md/fep/qm) and prints the plan.

        Raises ``ValueError`` for a target name with no registered aqueous property.

        :param targets:
            The targets.

        :return:
            The calculation plan.
        """

        plan = CalculationPlan()
        for target in targets:
            prop = aqueous_property(target.name)
            plan.jobs.setdefault(prop.job, []).append(target.name)
            plan.methods[target.name] = prop.method
        print_settings(self.ostream, "AutoCDAM Calculation Plan",
                       list(plan.methods.items()))
        return plan

    def single_point(self, molecule, *, settings=None, metadata=None) -> float:
        """
        Calculates the converged SCF energy of a real molecule (cached single point).

        :param molecule:
            The molecule.
        :param settings:
            The single-point settings (see ``SINGLE_POINT_SETTINGS``).
        :param metadata:
            A dictionary that receives ``cache_hit`` and ``scf_converged``.

        :return:
            The SCF energy in hartree.
        """

        return self.qm.energy(molecule, settings, metadata)

    def calculate(self, batch: GeneratedBatch, parameters: ParameterSet,
                  targets: Iterable[FitTarget],
                  **kwargs: Any) -> CalculationResult:
        """
        Computes every target for one parameter set.

        QM-only plans return before any OpenMM work. MD runs when MD properties are
        requested or to equilibrate the box before FEP (``equilibrate_first``). An MD
        or trajectory-analysis failure returns a failed result with NaN values.
        Raises ``ValueError`` if ``parameters`` differ from the batch's own parameters,
        and ``check_fep_backend`` errors if HFE is requested without a usable backend.

        :param batch:
            The generated batch.
        :param parameters:
            The parameter set of the batch.
        :param targets:
            The targets.
        :param kwargs:
            The options ``md_options``, ``fep_options``, ``qm_options`` and
            ``work_dir``.

        :return:
            The calculation result.
        """

        targets = list(targets)
        plan = self.plan(targets)
        if "fep" in plan.jobs:
            check_fep_backend()
        if batch.parameters is not None and batch.parameters.values != parameters.values:
            raise ValueError(
                "calculation parameters differ from the generated candidate batch"
            )
        self._print_header(batch, parameters, targets, kwargs)
        results = _Results(targets, plan)

        if "qm" in plan.jobs:
            self._calculate_qm(batch, targets, kwargs, results)
            if set(plan.jobs) == {"qm"}:
                return results.finish(parameters)

        self._initialize_platform()
        systems = self._prepare(batch, parameters)
        work_dir = Path(kwargs.get("work_dir", "."))
        md_options = md_settings(kwargs.get("md_options"))
        fep_options = kwargs.get("fep_options", {})

        final_state, md_analysis, r1_nm = None, {}, DEFAULT_ION_RADIUS_NM
        equilibrate = "fep" in plan.jobs and fep_options.get(
            "equilibrate_first", True)
        if "md" in plan.jobs or equilibrate:
            # HFE alone still runs the MD protocol: short FEP windows started
            # from the minimized template biased Al3+ charging by ~400 kJ/mol.
            final_state, md_results = self._execute_md_block(
                systems, plan.jobs.get("md", []), md_options)
            if final_state is None:
                return self._generate_failed_result(targets, parameters)
            if "md" in plan.jobs:
                r1_nm = md_results.pop("R1", r1_nm)
                md_analysis = md_results.pop("_analysis", {})
                self._record_md(md_results, md_analysis, md_options, results)
            else:
                r1_nm = md_results.get("R1", r1_nm)
                results.metadata["HFE"]["equilibration_md_options"] = md_options

        if "fep" in plan.jobs:
            self._calculate_hfe(systems, final_state, r1_nm, work_dir,
                                fep_options, batch, results)
        return results.finish(parameters, md_analysis=md_analysis)

    # ----------------------------------------------------------------- setup

    def _initialize_platform(self) -> None:
        """First OpenMM work of a calculation: report the default platform."""
        platform, _ = default_platform()
        self.ostream.print_info(f"OpenMM platform: {platform.getName()}")

    def _print_header(self, batch, parameters, targets, kwargs) -> None:
        """Print the sample, candidate, parameters, targets and paths of this calculation."""
        print_settings(self.ostream, "CDAM Property Calculation", [
            ("Sample", batch.metadata.get("sample_id")),
            ("Candidate", batch.metadata.get("candidate_id")),
            ("Potential Form", batch.metadata.get("potential_form")),
            ("Water Model", batch.metadata.get("water_model")),
            ("Parameters", parameters.to_dict() or "candidate defaults"),
            ("Targets", ", ".join(
                f"{t.name} [{t.unit}]" if t.unit else t.name for t in targets)),
            ("Working Directory", kwargs.get("work_dir")),
        ])

    def _prepare(self, batch: GeneratedBatch,
                 parameters: ParameterSet) -> OpenMMSystems:
        """The batch's OpenMM systems (built from its CDAM candidate if absent), solvated copied."""
        systems = batch.openmm_systems
        if systems is None:
            if batch.metadata.get("mode") != "cdam_candidate":
                raise ValueError(
                    "batch has no OpenMM systems and is not a CDAM candidate batch"
                )
            candidate = CDAMCandidate.from_generated_batch(
                batch).with_parameters(parameters)
            batch = candidate.to_generated_batch()
            batch.parameters = parameters
            systems = CoordinateGenerator(
                ostream=self.ostream).create_openmm_template(candidate)

        # A copy keeps repeated calculations independent (MD adds a barostat).
        solvated = mm.XmlSerializer.deserialize(
            mm.XmlSerializer.serialize(systems.solvated_system))
        return replace(systems, solvated_system=solvated)

    # ------------------------------------------------------------------- jobs

    def _calculate_qm(self, batch, targets, kwargs, results: _Results) -> None:
        """Compute QM_ENERGY (Hartree) on the first configuration's real atoms.

        Settings are layered, later wins: defaults, system and batch metadata, the
        target's metadata, then ``qm_options``.
        """
        target = next(t for t in targets if t.name == "QM_ENERGY")
        settings = {
            "qm_basis": DEFAULT_BASIS,
            "multiplicity": 1,
            **batch.system.metadata,
            **batch.metadata
        }
        settings.update(target.metadata)
        settings.update(kwargs.get("qm_options", {}))
        charge = settings.get("formal_charge", sum(batch.system.charges))
        multiplicity = settings.get("multiplicity", 1)

        molecule = real_atom_molecule(batch, charge, multiplicity)
        scf = {
            key: settings[key]
            for key in SINGLE_POINT_SETTINGS
            if key in settings
        }
        scf["basis"] = settings.get("qm_basis", DEFAULT_BASIS)
        print_setup(self.ostream, molecule, scf)
        info = {}
        results.values["QM_ENERGY"] = self.single_point(molecule,
                                                        settings=scf,
                                                        metadata=info)
        results.metadata["QM_ENERGY"].update({
            "qm_basis": settings["qm_basis"],
            "formal_charge": float(charge),
            "multiplicity": settings["multiplicity"],
            "scf_type": "restricted"
                        if settings["multiplicity"] == 1 else "unrestricted",
            "scf_converged": True,
            "cache_hit": info["cache_hit"],
            "scope": "first_frame_real_atoms_without_dummy_sites",
        })

    def _execute_md_block(
            self,
            systems: OpenMMSystems,
            md_targets,
            options=None) -> tuple[mm.State | None, dict[str, Any]]:
        """Run MD of the solvated box and analyse ``md_targets``; (None, {}) if either fails."""
        self.ostream.print_info("Starting NPT MD...")
        started = time.time()
        topology, positions = systems.solvated_topology, systems.solvated_positions
        trajectory = run_md(systems.solvated_system,
                            topology,
                            positions,
                            options,
                            ostream=self.ostream)
        if trajectory is None:
            return None, {}
        try:
            self.ostream.print_info("Analyzing trajectory...")
            self.ostream.flush()
            reference = np.array(positions.value_in_unit(unit.nanometer))
            md_results = analyze_trajectory(trajectory.frames, trajectory.boxes,
                                            topology, md_targets, reference)
        except Exception as error:
            self.ostream.print_warning(f"Trajectory analysis failed: {error}")
            return None, {}
        self.ostream.print_info(
            f"MD block completed in {time.time() - started:.1f} s")
        return trajectory.final_state, md_results

    def _record_md(self, md_results, analysis, md_options,
                   results: _Results) -> None:
        """Store MD values in ``results``; ``*_std`` entries become per-target spread metadata."""
        for name in results.plan.jobs["md"]:
            entry = results.metadata[name]
            entry["md_options"] = md_options
            entry["pilot"] = bool(md_options.get("pilot", False))
            entry["fast_mode"] = bool(md_options.get("fast_mode", False))
            if analysis:
                entry["analysis_method"] = analysis["method"]
                if name == "CN":
                    entry["cutoff_angstrom"] = analysis["cn_cutoff_angstrom"]
        if md_options.get("pilot"):
            results.warnings.append(
                "Pilot MD (0.05 ps) only tests the workflow; its values are "
                "not meaningful.")
        if md_options.get("fast_mode"):
            results.warnings.append(
                "Fast-mode MD (50 ps production); check the final parameters "
                "with the production protocol (fast_mode = False).")
        for key, value in md_results.items():
            if key.endswith("_std"):
                name = key[:-len("_std")]
                if name in results.metadata:
                    results.metadata[name]["spread"] = value
                    results.metadata[name][
                        "spread_definition"] = SPREAD_DEFINITIONS[name]
            else:
                results.values[key] = value

    def _calculate_hfe(self, systems: OpenMMSystems, final_state, r1_nm,
                       work_dir, fep_options, batch, results: _Results) -> None:
        """FEP free energy plus the finite-size correction for this box."""
        meta = results.metadata["HFE"]
        # A number alone does not establish overlap or MBAR convergence; a
        # study quality_check must accept live HFE explicitly.
        meta["quality_status"] = "unreviewed"
        if final_state is None:
            results.warnings.append(
                "HFE correction uses the default ion radius of 0.20 nm; "
                "no MD/RDF radius was calculated.")
        if fep_options:
            results.warnings.append(
                "FEP sampling settings were overridden; inspect the "
                "VeloxChem MBAR output before interpreting HFE.")
        fep_output = run_fep(systems,
                             final_state,
                             work_dir,
                             fep_options,
                             batch,
                             ostream=self.ostream)

        if final_state is not None:
            box = final_state.getPeriodicBoxVectors(asNumpy=True)
        else:
            box = systems.solvated_system.getDefaultPeriodicBoxVectors()
        box_nm = box[0][0].value_in_unit(unit.nanometers)
        density = count_waters(systems.solvated_topology) / box_nm**3
        charge, water_model = batch.metadata["formal_charge"], batch.metadata[
            "water_model"].lower()
        correction = finite_size_correction(charge, box_nm, r1_nm, density,
                                            water_model)
        print_settings(self.ostream, "HFE Finite-Size Correction", [
            ("Formal Charge", charge),
            ("Water Model", water_model),
            ("Box Length (x)", f"{box_nm:.4f} nm"),
            ("Ion Radius", f"{r1_nm:.4f} nm (" +
             ("RDF onset from MD"
              if final_state is not None else "default, no MD") + ")"),
            ("Water Number Density", f"{density:.3f} nm^-3"),
            ("Correction", f"{correction:.2f} kJ/mol"),
        ])

        results.values["HFE"] = fep_output["free_energy"] + correction
        meta["fep_output"] = fep_output
        uncertainty = fep_uncertainty(fep_output)
        if uncertainty is not None:
            meta["uncertainty"] = uncertainty
        meta["analytical_correction_kj_mol"] = correction

    def _generate_failed_result(self, targets: list[FitTarget],
                                parameters: ParameterSet) -> CalculationResult:
        """NaN-valued, ``failed=True`` result used when MD or its analysis crashed."""
        self.ostream.print_warning(
            "Short-circuiting calculation due to MD crash. Returning NaN.")
        return CalculationResult(
            values={t.name: float("nan") for t in targets},
            metadata={
                t.name: {
                    "method": "failed",
                    "error": "OpenMM crash during MD"
                } for t in targets
            },
            raw_data={"parameters": parameters.to_dict()},
            failed=True,
            failure_reason="MD simulation or trajectory analysis failed")
