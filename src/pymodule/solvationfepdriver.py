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

import math
import numbers
import itertools
from mpi4py import MPI
from io import StringIO
from contextlib import redirect_stderr
from pathlib import Path
import numpy as np
import sys
import time
import copy

from .veloxchemlib import mpi_master, bohr_in_angstrom, hartree_in_kjpermol
from .solvationbuilder import SolvationBuilder
from .mmforcefieldgenerator import MMForceFieldGenerator
from .outputstream import OutputStream
from .molecule import Molecule
from .molecularbasis import MolecularBasis
from .scfrestdriver import ScfRestrictedDriver
from .scfunrestdriver import ScfUnrestrictedDriver
from .cpcmdriver import CpcmDriver
from .errorhandler import assert_msg_critical

with redirect_stderr(StringIO()) as fg_err:
    try:
        from pymbar import MBAR, timeseries
    except ImportError:
        pass
    try:
        import openmm as mm
        import openmm.app as app
        import openmm.unit as unit
    except ImportError:
        pass


class SolvationFepDriver:
    """
    Computes the solvation free energy using OpenMM.

    analysis_method defaults to 'block_bootstrap': retain equilibrated frames,
    estimate g from sampled energies and adjacent-state energy differences,
    and scan circular block lengths at (1, 2, 4) times g. Each scale uses
    bootstrap_samples=200 MBAR fits and requires min_bootstrap_blocks=20
    per window. No usable scale gives NaN uncertainty, not zero. The optional
    'subsample' mode uses PyMBAR decorrelation followed by frame bootstrap.
    bootstrap_seed controls reproducibility. Per-stage sampling_diagnostics
    records dimensionless (kT) uncertainty, overlap, temporal and weight ESS.
    Stability and overlap thresholds are diagnostics, not convergence proofs.

    ion_corrections selects 'auto' (default), 'none', or 'custom'.
    Auto selects a bulk-water preset for spce, cspce, tip3p or ctip3p
    at 298.15 K and 1 atm with cubic PME/Ewald electrostatics. Non-SPC/E
    presets are provisional; see _WATER_CORRECTION_PRESETS and result metadata.
    For nonaqueous or undeclared solvents, and for the spc and tip4pew water
    models (no preset), auto reports uncorrected energies with an explicit
    diagnostic; use custom for charged solutes. File inputs declare water_model for aqueous
    runs; it is optional otherwise. A declaration assumes homogeneous bulk,
    and does not verify solvent composition or exclude interfaces in the files.
    Neutral solutes skip corrections. Effective water charges require custom.
    Monatomic PBC corrections require an explicit positive ion_radius_nm.

    Custom uses ion_correction_terms ('pbc_cubic', 'galvani', and/or
    'dielectric_legacy') and explicitly supplied correction_environment,
    correction_solvent, ion_reference_convention, charge_model ('unscaled' or
    'effective'), formal_charge, and term-specific parameters. pbc_cubic uses
    solvent_dielectric and the original approximate spherical-cavity formula
    at the input reference cell. galvani requires galvani_coefficient (signed
    additive kJ/mol per unit charge) and galvani_charge_basis ('formal' or
    'interaction'). dielectric_legacy needs target_dielectric, retains the
    initial-density CPCM approximation, and is restricted to unscaled charges.
    No mode implements a complete ECC electronic-solvation correction.

    The old apply_ion_corrections flag is optional: False explicitly selects
    none and True selects custom. Leave it None when using ion_corrections.
    Resolved preset values are reported in the result, not retained as user
    settings between runs. Charges are never modified by the correction modes.

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.

    Instance variables:
        - padding: The padding for the solvation box.
        - solvent_name: The name of the solvent.
        - resname: The residue name for the solute. Only needed when using input files where solute is not at top of (i.e., index 0) .pdb/.gro.
        - temperature: The temperature for the simulation.
        - pressure: The pressure for the simulation.
        - timestep: The timestep for the simulation.
        - num_equil_steps: The number of steps for the equilibration.
        - num_steps: The number of steps for the simulation.
        - num_snapshots: The number of snapshots to save.
        - cutoff: The cutoff for the nonbonded interactions.
        - constraints: The constraints for the simulation.
        - nonbondedMethod: The nonbonded method for the simulation.
        - solute_ff: The solute force field.
        - platform: The platform for the openMM simulation. None selects the fastest available.
        - precision: GPU precision ('mixed' default; 'single' or 'double').
        - platform_properties: Optional extra OpenMM platform properties (e.g. DeviceIndex).
        - vacuum_platform: Optional platform for the gas-phase stages (None default: same as platform).
        - integrator: 'langevin_middle' (default) or 'langevin' (original protocol).
        - friction: The Langevin collision rate.
        - hydrogen_mass: Optional hydrogen mass for repartitioning (e.g. 1.5 * unit.amu with 4 fs).
        - sequential_seeding: If True, start each lambda window from the final state of the previous one.
        - save_trajectory_xtc: If True, save the trajectories in XTC format.
        - save_system_xml: If True, save the system in XML format.
        - save_energies_txt: If True, save the energies from each simulation in txt format.
        - alpha: The alpha parameter for the Gaussian softcore potential.
        - beta: The beta parameter for the Gaussian softcore potential.
        - x: The x parameter for the Gaussian softcore potential.
        - ion_correction_box: Box used by the pbc_cubic correction: 'average' (default, mean box of the coupled solvated window) or 'input' (starting box).
        - linear_coupling: If True (default), stages 1 and 3 use a Hamiltonian linear in lambda (solute epsilon times (1-lambda)^2, explicit linearly scaled intramolecular pairs). False gives the original parameter scaling.
        - gsc_placeholder_sigma: Sigma used in the GSC term for particles with epsilon = 0 (default 0.1 nm, the VeloxChem water convention; None keeps the force-field placeholder).
        - lambdas_stage1: The lambdas for stage 1.
        - lambdas_stage2: The lambdas for stage 2.
        - lambdas_stage3: The lambdas for stage 3.
        - lambdas_stage4: The lambdas for stage 4 (unused with keep_intramolecular_gsc=True).
        - alchemical_scheme: 'annihilation' (default; intramolecular nonbonded terms switched off, vacuum leg) or 'decoupling' (only solute-environment interactions switched off, GSC between solute and environment, no vacuum leg).
        - keep_intramolecular_gsc: If True (default), the solute-solute GSC stays on in the decoupled state; stage 2 removes only the solute-environment GSC and stage 4 (identically zero) is not run. False gives the original cycle.
        - u_kln: The potential energies across stages.
        - stage: The current stage.
        - final_free_energy: The final free energy.
        - delta_f: The free energy calculations for each stage.
    """

    def __init__(self, comm=None, ostream=None):
        """
        Initialize the OpenMMSolvator class.
        """

        assert_msg_critical(
            'pymbar' in sys.modules and 'openmm' in sys.modules,
            'pymbar and OpenMM are required for SolvationFepDriver.')

        if comm is None:
            comm = MPI.COMM_WORLD

        if ostream is None:
            if comm.Get_rank() == mpi_master():
                ostream = OutputStream(sys.stdout)
            else:
                ostream = OutputStream(None)

        # MPI information
        self.comm = comm
        self.rank = comm.Get_rank()
        self.size = comm.Get_size()

        # Output stream
        self.ostream = ostream

        # Create directory for storing all generated data
        self.output_folder = "solvation_fep_output"

        # Options for the SolvationBuilder
        self.padding = 1.0
        self.solvent_name = 'spce'
        self.water_model = 'spce' # Independent attribute for physical corrections (Galvani, etc.)
        self.resname = None

        # Ensemble and MD options
        self.temperature = 298.15
        self.kT = (unit.AVOGADRO_CONSTANT_NA * unit.BOLTZMANN_CONSTANT_kB * (self.temperature * unit.kelvin)).in_units_of(unit.kilojoule_per_mole)
        self.pressure = 1 * unit.atmospheres
        self.timestep = 2.0 * unit.femtoseconds
        self.num_em_steps = 0  # no limit
        self.num_equil_steps = 5000  # 10 ps
        self.num_steps = 500000  # 1 ns
        self.num_snapshots = 500
        self.cutoff = 1.0 * unit.nanometers
        self.constraints = app.HBonds
        self.nonbondedMethod = app.PME

        # Integration and sampling options.
        # 'langevin_middle' (BAOAB-type splitting) samples configurations more
        # accurately than 'langevin' at the same timestep. Select 'langevin'
        # and sequential_seeding=False to reproduce the original protocol.
        self.integrator = 'langevin_middle'
        self.friction = 1.0 / unit.picoseconds
        # Optional hydrogen mass repartitioning (e.g. 1.5 * unit.amu together
        # with timestep = 4 fs). Masses do not affect configurational free
        # energies; water is never repartitioned. Requires HBonds constraints.
        self.hydrogen_mass = None
        # Start every lambda window from the final state (positions, box and
        # velocities) of the previous window of the same phase. Only the first
        # window of each phase is energy-minimized.
        self.sequential_seeding = True

        # Other objects
        self.solute_ff = None
        self.platform = None
        # GPU precision: 'mixed' is recommended for energy reweighting.
        self.precision = 'mixed'
        self.platform_properties = None
        # Optional separate platform for the vacuum stages (e.g. 'CPU' or
        # 'Reference'). None (default) uses the same choice as the solvated
        # stages, as in the original driver. OpenMM's CPU platform has a large
        # per-step overhead for few-atom systems, so benchmark before changing.
        self.vacuum_platform = None
        self.save_trajectory_xtc = False
        self.save_system_xml = False
        self.save_energies_txt = False

        # Alchemical parameters
        self.alpha = 5 * unit.kilocalories_per_mole
        # GSC range for particles without LJ (epsilon = 0). Their NonbondedForce
        # sigma is an arbitrary placeholder (0.1 nm in VeloxChem, 1 nm in OpenMM
        # water XML files, often 0 in GROMACS topologies) that only the GSC
        # would read. All such particles use this sigma in the GSC term, so the
        # intermediate states do not depend on the force-field file convention.
        # None uses the NonbondedForce sigma as is (original behaviour).
        self.gsc_placeholder_sigma = 0.1 * unit.nanometers
        # Linear coupling in stages 1 and 3 (default): solute epsilon is scaled
        # by (1-lambda)^2, so solute-environment LJ is exactly linear, and all
        # intramolecular solute pairs that are not excluded become explicit pair
        # terms scaled by (1-lambda). The hybrid Hamiltonian is then linear in
        # lambda, as in the dual-topology GSC scheme of Li and Nam, apart from a
        # charged solute's periodic self-image energy (quadratic, configuration
        # independent). End states are unchanged. False restores the original
        # parameter scaling (charge and epsilon times (1-lambda)).
        self.linear_coupling = True
        # Decoupled reference state. True (default since 2026-10-22): the
        # solute-solute GSC (pairs that are not excluded) stays on in the
        # decoupled state of both legs and stage 2 removes only the
        # solute-environment GSC. The decoupled solute does not interact with
        # its environment and has the same Hamiltonian in both legs, so the
        # free energy is unchanged, while the stage-3 end state already is the
        # vacuum reference: stage 4 is identically zero and is not run. This
        # avoids the large, slowly converging switch-off of the intramolecular
        # GSC, which cancels between stages 2 and 4 only in exact arithmetic.
        # False restores the original cycle: stages 2 and 4 remove all GSC, so
        # the decoupled solute keeps only its bonded terms.
        self.keep_intramolecular_gsc = True
        # Alchemical scheme. 'annihilation' (default): stages 1 and 3 switch
        # off all solute nonbonded terms, intramolecular ones included, so a
        # vacuum leg is needed. 'decoupling': only solute-environment
        # interactions are switched off; all intramolecular nonbonded terms
        # (1-4 and 1-5+ pairs, as unscaled explicit exceptions) stay at full
        # strength in every state. The GSC is used as in annihilation (on in
        # stage 1, removed in stage 2), but only between solute and
        # environment: the full intramolecular interactions already prevent
        # collapse. The decoupled state is then environment x gas-phase
        # solute, so stages 3 and 4 are zero and not run (dG = S1 + S2). This
        # avoids artificial conformational changes along lambda in molecules
        # whose torsion terms compensate 1-4 electrostatics (e.g. carboxylic
        # acids: anti becomes favoured once the 1-4 terms are scaled down).
        # Requires linear_coupling=True; keep_intramolecular_gsc is not used.
        self.alchemical_scheme = 'annihilation'
        self.beta = 5
        self.x = 4

        # Lambda schedules (defaults since 2026-10-22, from the MNSol lambda
        # benchmark with keep_intramolecular_gsc=True and linear_coupling=True:
        # dense-grid subsets with per-molecule bias <= 1 kJ/mol and pairwise
        # neighbour overlap >= 0.03; confirmed with independent 1 ns runs).
        # Original schedules: stage 1/3 [0, .2, .4, .6, .8, 1], stage 2
        # [1, .8, .6, .4, .3, .2, .15, .1, .05, .03, 0], stage 4 [1, .5, 0].
        # Stage 1 needs the extra windows near the coupled end (lambda -> 0).
        self.lambdas_stage1 = [0.0, 0.05, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0]
        # Stage 2 is denser where the GSC is weak (water fills the solute volume).
        self.lambdas_stage2 = [1.0, 0.7, 0.5, 0.3, 0.2, 0.15, 0.1, 0.07, 0.05, 0.03, 0.01, 0.0]
        # Stage 3 (vacuum) is smooth.
        self.lambdas_stage3 = [0.0, 0.3, 0.7, 1.0]
        # Stage 4 is only run with keep_intramolecular_gsc=False; the original
        # 3 windows are strongly biased for large flexible solutes (profluralin
        # +12 kJ/mol); 1, .5, .25, .1, .05, 0 agreed with a dense 7-window
        # reference to <= 0.2 kJ/mol (and that with 15 x 1 ns runs to <= 1.5).
        self.lambdas_stage4 = [1.0, 0.5, 0.25, 0.1, 0.05, 0.0]

        # Auto selects a documented bulk-water preset from the declared model.
        # Available terms: 'pbc_cubic', 'galvani', 'dielectric_legacy'.
        self.ion_corrections = 'auto'
        self.apply_ion_corrections = None  # optional legacy override
        self._auto_bulk_solvent = None
        self._input_formal_charge = None
        self.ion_correction_terms = ()
        self.correction_environment = None  # must be 'homogeneous_bulk'
        self.correction_solvent = None      # explicit model/solvent label
        self.ion_reference_convention = None  # user-specified target reference
        self.charge_model = None            # 'unscaled' or 'effective'
        self.formal_charge = None           # integer chemical charge, supplied
        self.interaction_charge = None      # measured from input force field
        self.ion_charge = None              # compatibility alias: unrounded FF charge
        self.solute_volume_nm3 = None
        self.ion_radius_nm = None           # required for monatomic PBC correction
        self.solvent_dielectric = None      # dielectric used in correction formula
        self.target_dielectric = None      # required for legacy dielectric term
        self.galvani_coefficient = None     # signed additive kJ/mol per unit charge
        self.galvani_charge_basis = None    # explicitly 'formal' or 'interaction'
        self._ion_correction_box_nm = None
        # Box for the pbc_cubic correction: 'average' (default) uses the mean
        # cubic box length sampled in the fully coupled solvated window
        # (stage 1, lambda = 0) under NPT; 'input' uses the starting box.
        self.ion_correction_box = 'average'

        # Temporal analysis; lengths and g are measured in saved-frame units.
        self.analysis_method = 'block_bootstrap'  # alternative: 'subsample'
        self.bootstrap_samples = 200
        self.bootstrap_seed = None
        self.block_length_factors = (1.0, 2.0, 4.0)
        self.min_bootstrap_blocks = 20
        self.u_kn_matrices = []
        self.stage = 1

        # Final energies
        self.final_free_energy = 0.0
        self.delta_f = None

    def _reset_calculation_state(self, clear_inputs=False):
        """Discard derived state, preserving explicitly configured user options."""
        self.interaction_charge = self.ion_charge = None
        self._ion_correction_box_nm = None
        self._ion_preset_name = None
        self._active_ion_mode = None
        self._ion_preset_note = None
        self._endpoint_signature = None
        self._endpoint_charges = None
        self.dG_PBC = self.dG_phi = self.dG_D = 0.0
        self.u_kn_matrices = []
        self.snapshots = []
        self.N_k = []
        self._equilibration_diagnostics = []
        self._seed_state = None
        self._has_barostat = None
        self._average_box = None
        self._pbc_box_used_nm = None
        self.stage = 1
        self.final_free_energy = None
        self.delta_f = None
        if clear_inputs:
            self._auto_bulk_solvent = self._input_formal_charge = None
            self._supplied_solute_ff = None
            self._generated_solute_ff = None
            self.solute_ff = None
            for name in ('system_pdb', 'system_xml', 'solute_pdb', 'solute_xml',
                         'system_gro', 'system_top', 'solute_gro', 'solute_top',
                         'other_xml_files', 'base_system', 'topology', 'positions',
                         'vacuum_system', 'vacuum_topology', 'vacuum_positions'):
                setattr(self, name, None)

    @staticmethod
    def _validate_solute_force_field(molecule, ff):
        """Require charges indexed in the supplied molecule's atom order.

        No charge fitting, scaling or automatic atom remapping is performed.
        Element order cannot distinguish permutations of equivalent atoms:
        callers must retain the atom indexing used during parametrization.
        """
        n = molecule.number_of_atoms()
        reference = getattr(ff, 'molecule', None)
        if reference is None or list(reference.get_labels()) != list(molecule.get_labels()):
            raise ValueError('The force field must retain its reference molecule '
                             'with the same atom count and element order as the solute.')
        atoms = getattr(ff, 'atoms', {})
        if set(atoms) != set(range(n)):
            raise ValueError('Force-field atom indices must match solute indices 0..N-1.')
        charges = np.asarray(getattr(ff, 'partial_charges', None), dtype=float)
        if charges.shape != (n,) or not np.all(np.isfinite(charges)):
            raise ValueError('partial_charges must contain one finite charge per solute atom.')
        exported = np.asarray([atoms[i]['charge'] for i in range(n)], dtype=float)
        if exported.shape != (n,) or not np.all(np.isfinite(exported)):
            raise ValueError('Exported atom charges must contain one finite value per atom.')
        if not np.allclose(charges, exported, rtol=0, atol=1e-6):
            raise ValueError('partial_charges and atoms[i]["charge"] disagree; '
                             'update both on the input force-field generator.')

    def _validate_endpoint(self, system, topology, positions, region, vacuum=False):
        """Check the unscaled endpoint before applying any lambda changes."""
        atoms = list(topology.atoms())
        n = system.getNumParticles()
        if len(atoms) != n or len(positions) != n:
            raise ValueError('System, topology and coordinate atom counts disagree.')
        if not region or len(set(region)) != len(region) or any(i < 0 or i >= n for i in region):
            raise ValueError('Invalid or empty solute atom mapping.')
        if vacuum and region != list(range(n)):
            raise ValueError('The vacuum endpoint must contain only the selected solute.')
        signature = tuple((atoms[i].name,
                           atoms[i].element.symbol if atoms[i].element else None)
                          for i in region)
        forces = [f for f in system.getForces() if isinstance(f, mm.NonbondedForce)]
        if len(forces) != 1 or forces[0].getNumParticles() != n:
            raise ValueError('Endpoint validation requires one NonbondedForce with one entry per atom.')
        charges = np.array([forces[0].getParticleParameters(i)[0].value_in_unit(
            unit.elementary_charge) for i in region], dtype=float)
        if not np.all(np.isfinite(charges)):
            raise ValueError('Endpoint contains non-finite solute charges.')
        if self._endpoint_signature is None:
            self._endpoint_signature = signature
            self._endpoint_charges = charges.copy()
        elif (signature != self._endpoint_signature or
              charges.shape != self._endpoint_charges.shape or
              not np.allclose(charges, self._endpoint_charges, rtol=0, atol=1e-6)):
            raise ValueError('Solute atom names/order, elements or unscaled charges '
                             'differ between endpoints. Supply matching atom mappings.')

    def compute(self, molecule, ff_gen_solute=None, solvent='spce', solvent_molecule=None, ff_gen_solvent=None, target_density=None):
        """
        Run the solvation free energy calculation using OpenMM.

        Physical-state charges come from the input force field. Only the
        alchemical lambda scaling required by FEP is applied in this driver.

        :param molecule:
            The VeloxChem molecule object to solvate.
        :param ff_gen_solute:
            The prepared MMForceFieldGenerator for the solute. Its charges are
            used as supplied; compute any charge-model corrections or scaling
            externally before passing this object. If None, the standard
            SolvationBuilder force-field generation is used.
        :param solvent:
            The solvent to use for solvation. Default is spce water.
            Available options: 'cspce', 'ctip3p', 'spc', 'spce', 'tip3p', 'tip4pew', 'ethanol', 'methanol', 'acetone', 'chloroform',
            'hexane', 'toluene', 'dcm', 'benzene', 'dmso', 'thf', 'acetonitrile', 'dmf','other' or 'itself'.
        :param solvent_molecule:
            The VeloxChem molecule object for the solvent. Mandatory for 'other'.
        :param ff_gen_solvent:
            The force field generator for the solvent. If None, it will be calculated.
        :param target_density:
            The target density for the solvent. Mandatory for 'other' or 'itself' solvent options.
        :return:
            A dictionary containing the free energy calculations and the uncertainty for each stage, and the final free energy.
        """

        self._reset_calculation_state(clear_inputs=True)
        self.water_model = solvent
        self._auto_bulk_solvent = solvent if ff_gen_solvent is None else None
        self._input_formal_charge = float(molecule.get_charge())

        if ff_gen_solute is not None:
            self._validate_solute_force_field(molecule, ff_gen_solute)
            self._supplied_solute_ff = ff_gen_solute

        sol_builder = SolvationBuilder()
        sol_builder.steps = 50000

        sol_builder.solvate(solute=molecule,
                            solvent=solvent,
                            solvent_molecule=solvent_molecule,
                            padding=self.padding,
                            target_density=target_density,
                            neutralize=False,
                            equilibrate=True)

        self.solvent_name = solvent
        sol_builder.write_openmm_files(ff_gen_solute, ff_gen_solvent)

        if self.save_trajectory_xtc or self.save_system_xml or self.save_energies_txt:
            Path(self.output_folder).mkdir(parents=True, exist_ok=True)

        if ff_gen_solute is None:
            self._generated_solute_ff = sol_builder.solute_ff
            self.solute_ff = self._generated_solute_ff
        else:
            self.solute_ff = ff_gen_solute

        self._validate_solute_force_field(molecule, self.solute_ff)
        return self._run_stages()

    def compute_solvation_free_energy(self, solute_ff, sol_builder, solute_name='solute', water_model=None):

        self._reset_calculation_state(clear_inputs=True)
        self._validate_solute_force_field(sol_builder.solute, solute_ff)
        sol_builder.write_gromacs_files(solute_ff=solute_ff)

        solute_gro = f'{solute_name}.gro'
        solute_top = f'{solute_name}.top'

        return self.compute_solvation_from_gromacs_files("system.gro", "system.top", solute_gro, solute_top,
            water_model=water_model, solvent=getattr(sol_builder, "solvent_name", None))

    def compute_solvation_from_openmm_files(self, solute_pdb, solute_xml, solvent=None, solvent_molecule=None, ff_gen_solvent=None, target_density=None, system_pdb=None, other_xml_files=None, system_xml=None, water_model=None):
        """
        Run the solvation free energy calculation using OpenMM.
        :param water_model:
            For aqueous inputs: spce, cspce, tip3p or ctip3p.
            Omit for nonaqueous inputs; auto then applies no ion correction.
            Explicit solvent is accepted as a compatibility alias.
            Charges are read from the solute NonbondedForce, never from PDB.
        :param solute_pdb:
            The PDB file with the solute coordinates.
        :param solute_xml:
            The XML file with the solute force field.
        :param other_xml_files:
            A list with the XML files needed for the rest of the system.
        :param system_pdb:
            The PDB file with the solvated system coordinates. If not provided, the system will be built.
        :param solvent:
            The solvent to use for solvation. Default is spce water.
            Available options: 'cspce', 'ctip3p', 'spc', 'spce', 'tip3p', 'tip4pew', 'ethanol', 'methanol', 'acetone', 'chloroform',
            'hexane', 'toluene', 'dcm', 'benzene', 'dmso', 'thf', 'acetonitrile', 'dmf', 'other' or 'itself'.
        :param solvent_molecule:
            The VeloxChem molecule object for the solvent. Mandatory for 'other'.
        :param ff_gen_solvent:
            The force field generator for the solvent. If None, it will be calculated.
        :param target_density:
            The target density for the solvent. Mandatory for 'other' or 'itself' solvent options.

        :return:
            A dictionary containing the free energy calculations and the uncertainty for each stage, and the final free energy.
        """

        self._reset_calculation_state(clear_inputs=True)
        solvent = self._declare_input_solvent(water_model, solvent)
        if system_pdb is None and solvent is None:
            raise ValueError('Specify solvent when building a solvated system, '
                             'e.g. solvent="thf" or water_model="cspce".')
        self.solute_pdb = solute_pdb
        self.solute_xml = solute_xml
        self._input_formal_charge = None

        # To run with pre-built OMM system objects
        if system_xml is not None and system_pdb is not None:
            self.system_pdb = system_pdb
            self.system_xml = system_xml
            self.solute_pdb = solute_pdb
            self.solute_xml = solute_xml
            self.solvent_name = 'omm_system_xml'

        elif system_pdb is None:
            self._auto_bulk_solvent = solvent if ff_gen_solvent is None else None
            self.ostream.print_info(f"No system_pdb provided, using SolvationBuilder to solvate the system with {solvent} solvent")
            self.ostream.flush()
            molecule = Molecule.read_pdb_file(solute_pdb)

            ffgen_solute = MMForceFieldGenerator()
            ffgen_solute.ostream.mute()
            ffgen_solute.create_topology(molecule, resp=False)

            sol_builder = SolvationBuilder()
            sol_builder.steps = 50000
            sol_builder.write_pdb_only = True

            sol_builder.solvate(solute=molecule,
                                solvent=solvent,
                                solvent_molecule=solvent_molecule,
                                padding=self.padding,
                                target_density=target_density,
                                neutralize=False,
                                equilibrate=True)

            sol_builder.write_openmm_files(solute_ff=ffgen_solute,solvent_ffs=ff_gen_solvent)
            self.system_pdb = 'system.pdb'
            self.solvent_name = 'omm_files'
            self.other_xml_files = ['solvent_1.xml']
            if other_xml_files:
                self.other_xml_files.extend(other_xml_files)

        else:
            self.system_pdb = system_pdb
            assert_msg_critical(other_xml_files is not None,
                                "SolvationFepDriver: other_xml_files must be provided if system_pdb is provided")
            self.other_xml_files = other_xml_files
            self.solvent_name = 'omm_files'

        if self.save_trajectory_xtc or self.save_system_xml or self.save_energies_txt:
            Path(self.output_folder).mkdir(parents=True, exist_ok=True)

        return self._run_stages()

    def compute_solvation_from_gromacs_files(self, system_gro, system_top, solute_gro, solute_top, water_model=None, solvent=None):
        """
        Run the solvation free energy calculation using OpenMM.

        :param water_model:
            For aqueous inputs: spce, cspce, tip3p or ctip3p.
            Omit for nonaqueous inputs; auto then applies no ion correction.
            Charges are read from the force field, not GRO coordinates.
        :param system_gro:
            The GRO file with the system coordinates.
        :param system_top:
            The TOP file with the system topology.
        :param solute_gro:
            The GRO file with the solute coordinates.
        :param solute_top:
            The TOP file with the solute topology.
        :return:
            A dictionary containing the free energy calculations and the uncertainty for each stage, and the final free energy.
        """

        self._reset_calculation_state(clear_inputs=True)
        self._declare_input_solvent(water_model, solvent)
        # Special name for the solvent
        self.solvent_name = 'gro_files'
        self.system_gro = system_gro
        self.system_top = system_top
        self.solute_gro = solute_gro
        self.solute_top = solute_top

        if self.save_trajectory_xtc or self.save_system_xml or self.save_energies_txt:
            Path(self.output_folder).mkdir(parents=True, exist_ok=True)

        return self._run_stages()

    def compute_solvation_from_memory(self, solvated_system, solvated_topology, solvated_state, vacuum_system, vacuum_positions, vacuum_topology=None, water_model=None, solvent=None):
        """Use caller-owned endpoints; vacuum_topology must describe the solute."""
        self._reset_calculation_state(clear_inputs=True)
        if vacuum_topology is None:
            raise ValueError("Supply vacuum_topology to validate the solute atom mapping.")
        self._declare_input_solvent(water_model, solvent)
        self.vacuum_topology = vacuum_topology
        self.solvent_name = 'omm_memory'

        self.base_system = copy.deepcopy(solvated_system)
        self.topology = solvated_topology
        self.positions = solvated_state.getPositions(asNumpy=True)

        box_vectors = solvated_state.getPeriodicBoxVectors()
        self.base_system.setDefaultPeriodicBoxVectors(*box_vectors)

        self.vacuum_system = vacuum_system
        self.vacuum_positions = vacuum_positions

        if self.save_trajectory_xtc or self.save_system_xml or self.save_energies_txt:
            Path(self.output_folder).mkdir(parents=True, exist_ok=True)

        return self._run_stages()

    def _estimate_solute_volume(self, topology, positions):
        alchemical_region, _ = self._get_alchemical_region(topology)
        labels = []
        coords_angstrom = []
        atoms = list(topology.atoms())

        for idx in alchemical_region:
            atom = atoms[idx]
            if hasattr(positions[idx], 'value_in_unit'):
                coord = positions[idx].value_in_unit(unit.angstroms)
            else:
                coord = positions[idx] * 10.0
            coords_angstrom.append(coord)
            elem = atom.element.symbol if atom.element is not None else 'C'
            labels.append(elem)

        mol = Molecule(labels, coords_angstrom, 'angstrom')
        coords_bohr = mol.get_coordinates_in_bohr()
        radii_bohr = mol.vdw_radii_to_numpy()

        bohr_to_nm = bohr_in_angstrom() / 10.0
        solute_coords = coords_bohr * bohr_to_nm
        radii = radii_bohr * bohr_to_nm

        min_coords = np.min(solute_coords - radii[:, None], axis=0)
        max_coords = np.max(solute_coords + radii[:, None], axis=0)

        spacing = 0.02
        x = np.arange(min_coords[0], max_coords[0], spacing)
        y = np.arange(min_coords[1], max_coords[1], spacing)
        z = np.arange(min_coords[2], max_coords[2], spacing)

        xv, yv, zv = np.meshgrid(x, y, z, indexing='ij')
        grid_points = np.vstack([xv.ravel(), yv.ravel(), zv.ravel()]).T

        inside = np.zeros(grid_points.shape[0], dtype=bool)
        for i, coord in enumerate(solute_coords):
            dist_sq = np.sum((grid_points - coord)**2, axis=1)
            inside |= (dist_sq <= radii[i]**2)

        voxel_volume = spacing**3
        volume_nm3 = float(np.sum(inside) * voxel_volume)

        return volume_nm3

    def _calculate_dielectric_correction(self, eps_s):
        """
        Calculates the dielectric permittivity correction (dG_D) using VeloxChem's
        CPCM driver and the legacy initial-density approximation. This is not
        an ECC electronic-solvation correction. Both dielectrics are explicit.
        """
        alchemical_region, _ = self._get_alchemical_region(self.vacuum_topology)
        labels = []
        coords_angstrom = []
        atoms = list(self.vacuum_topology.atoms())

        for idx in alchemical_region:
            atom = atoms[idx]
            if hasattr(self.vacuum_positions[idx], 'value_in_unit'):
                coord = self.vacuum_positions[idx].value_in_unit(unit.angstroms)
            else:
                coord = self.vacuum_positions[idx] * 10.0
            coords_angstrom.append(coord)
            elem = atom.element.symbol if atom.element is not None else 'C'
            labels.append(elem)

        mol = Molecule(labels, coords_angstrom, 'angstrom')
        mol.set_charge(float(self.formal_charge))

        silent_stream = OutputStream(None)
        basis = MolecularBasis.read(mol, "STO-3G", ostream=silent_stream)

        if mol.get_multiplicity() > 1:
            scf_drv = ScfUnrestrictedDriver(self.comm, silent_stream)
        else:
            scf_drv = ScfRestrictedDriver(self.comm, silent_stream)

        den_mat_tuple = scf_drv.gen_initial_density_sad(mol, basis, basis)
        total_den_mat = scf_drv._get_total_density_matrix(den_mat_tuple)

        # 1. Calculate Heterogeneous Potential Energy (Solvent = eps_s)
        cpcm_het = CpcmDriver(self.comm, silent_stream)
        cpcm_het.epsilon = eps_s
        cpcm_het.init(mol, basis, do_nuclear=True)
        E_het, _ = cpcm_het.compute_gs_fock(mol, basis, total_den_mat, 1.0e-8)

        # 2. Calculate target-permittivity energy (explicit user setting)
        cpcm_exp = CpcmDriver(self.comm, silent_stream)
        cpcm_exp.epsilon = self.target_dielectric
        cpcm_exp.init(mol, basis, do_nuclear=True)
        E_exp, _ = cpcm_exp.compute_gs_fock(mol, basis, total_den_mat, 1.0e-8)

        dg_d = (E_exp - E_het) * hartree_in_kjpermol()

        return float(dg_d)

    # Signed additive coefficients preserve the original driver's convention.
    # These are pragmatic presets, NOT uniformly validated MD corrections.
    # Modified-model Galvani numbers originate in 3D-RISM (Misin et al.,
    # DOI 10.1021/acs.jpcb.5b10809); their transfer to MD is provisional.
    # 78.4 is an experimental-water approximation, not a measured model epsilon.
    _WATER_CORRECTION_PRESETS = {
        'spce': (71.0, -57.7, 'Existing SPC/E manuscript convention.'),
        'cspce': (78.4, -56.2, 'Provisional: 3D-RISM Galvani coefficient; '
                  'experimental-water dielectric approximation for MD.'),
        'tip3p': (78.4, -48.1, 'Provisional: original-driver Galvani coefficient '
                  '(MD provenance unverified); experimental-water dielectric approximation.'),
        'ctip3p': (78.4, -52.5, 'Provisional: 3D-RISM Galvani coefficient; '
                   'experimental-water dielectric approximation for MD.'),
    }

    def _declare_input_solvent(self, water_model=None, solvent=None):
        """A water declaration enables bulk-water presets; omission does not."""
        water = self._normalize_water_model(water_model) if water_model is not None else None
        if solvent is not None:
            if not isinstance(solvent, str) or not solvent.strip():
                raise ValueError('solvent must be a nonempty name or None.')
            solvent = solvent.strip().lower()
            key = solvent.replace('/', '').replace('-', '')
            if key in self._WATER_CORRECTION_PRESETS:
                solvent = key
            if water is not None and solvent != water:
                raise ValueError('water_model and solvent disagree.')
        selected = water or solvent
        self.water_model = selected if selected in self._WATER_CORRECTION_PRESETS else None
        self._auto_bulk_solvent = selected
        return selected

    @staticmethod
    def _normalize_water_model(model):
        if not isinstance(model, str) or not model.strip():
            raise ValueError('Specify water_model for file/memory input, e.g. '
                             'water_model="cspce". This declares bulk aqueous solvent.')
        key = model.strip().lower().replace('/', '').replace('-', '')
        if key not in SolvationFepDriver._WATER_CORRECTION_PRESETS:
            raise ValueError(f'Unsupported water model {model!r}; supported models: '
                             'spce, cspce, tip3p, ctip3p.')
        return key

    def _resolve_auto_ion_corrections(self):
        """Resolve the preset after reading the actual endpoint charges."""
        solvent = self.correction_solvent or self._auto_bulk_solvent
        key = (solvent.strip().lower().replace('/', '').replace('-', '')
               if isinstance(solvent, str) else None)
        if key not in self._WATER_CORRECTION_PRESETS:
            self.apply_ion_corrections = False
            self.ion_correction_terms = ()
            self._ion_preset_name = 'no_automatic_solvent_correction'
            self._ion_preset_note = (
                f'No automatic ion-correction preset for solvent {solvent!r}. '
                'Reported energy is uncorrected. For aqueous file input with spce, '
                'cspce, tip3p or ctip3p, specify water_model; otherwise a custom '
                'scheme may be supplied.')
            self.ostream.print_info(self._ion_preset_note)
            return
        q = self.interaction_charge
        # Even a scaled ion with an accidentally integral net charge must agree
        # with any chemical charge supplied by the user or molecule object.
        if self.charge_model not in (None, 'unscaled'):
            raise ValueError("Auto does not support effective/ECC charges; use "
                             "ion_corrections='custom' or 'none'.")
        if not np.isclose(q, round(q), rtol=0, atol=1e-5):
            raise ValueError("Auto requires an integer-like FF charge; use custom "
                             "with separate formal and interaction charges for ECC.")
        formal = self.formal_charge
        if formal is None:
            formal = self._input_formal_charge
        if formal is None:
            formal = int(round(q))
        if (not np.isfinite(formal)
                or not np.isclose(formal, round(formal), rtol=0, atol=1e-8)
                or not np.isclose(q, formal, rtol=0, atol=1e-5)):
            raise ValueError('Auto requires matching integer formal and FF charges; '
                             'use custom for effective-charge models.')
        self.formal_charge = int(round(formal))
        if self.formal_charge == 0:
            self.apply_ion_corrections = False
            self.ion_correction_terms = ()
            self._ion_preset_name = 'neutral_no_ion_correction'
            return
        environment = self.correction_environment
        if environment is None and self._auto_bulk_solvent is not None:
            environment = 'homogeneous_bulk'
        solvent = self.correction_solvent or self._auto_bulk_solvent
        if environment != 'homogeneous_bulk':
            raise ValueError('Auto requires homogeneous_bulk. For imported systems, '
                             'explicitly declare correction_environment and '
                             'correction_solvent, or select none/custom.')
        model = self._normalize_water_model(solvent)
        if (self._auto_bulk_solvent is not None
                and self._normalize_water_model(self._auto_bulk_solvent) != model):
            raise ValueError('Declared correction solvent conflicts with the built solvent.')
        eps, coefficient, note = self._WATER_CORRECTION_PRESETS[model]
        self._ion_preset_note = note
        if model != 'spce':
            self.ostream.print_info(f'{model} automatic ion corrections: {note}')
        pressure_atm = self.pressure.value_in_unit(unit.atmospheres)
        if (not np.isclose(self.temperature, 298.15, rtol=0, atol=0.01)
                or not np.isclose(pressure_atm, 1.0, rtol=0, atol=0.01)):
            raise ValueError('The auto water presets are for 298.15 K and 1 atm; '
                             'use custom at other state points.')
        preset = {
            'ion_correction_terms': ('pbc_cubic', 'galvani'),
            'correction_environment': 'homogeneous_bulk',
            'correction_solvent': model,
            'ion_reference_convention': ('manuscript_spce_additive_minus57.7_kj_mol_per_e'
                if model == 'spce' else 'original_driver_signed_additive_water_convention'),
            'charge_model': 'unscaled',
            'solvent_dielectric': eps,
            'galvani_coefficient': coefficient,
            'galvani_charge_basis': 'formal',
        }
        # Avoid silently ignoring user customization; select custom for changes.
        for name, value in preset.items():
            current = getattr(self, name)
            if current is not None and current != () and current != value:
                raise ValueError(f'{name} conflicts with the auto preset; '
                                 "select ion_corrections='custom'.")
        if self.target_dielectric is not None:
            raise ValueError('Auto excludes dielectric_legacy; use custom for that term.')
        for name, value in preset.items():
            setattr(self, name, value)
        self._ion_preset_name = f'bulk_{model}_298.15K_1atm_v1'
        self.apply_ion_corrections = True
        self._validate_ion_correction_settings()

    def _validate_ion_correction_settings(self):
        """Validate declared correction assumptions before production sampling.

        These terms do not implement a complete ECC thermodynamic cycle.
        Galvani coefficients are signed additive coefficients, not voltages.
        The user must match their sign, charge basis and reference convention.
        """
        if not self.apply_ion_corrections:
            return
        terms = self.ion_correction_terms
        allowed = {'pbc_cubic', 'galvani', 'dielectric_legacy'}
        if (not isinstance(terms, (tuple, list)) or not terms
                or any(term not in allowed for term in terms)
                or len(set(terms)) != len(terms)):
            raise ValueError('Select unique ion_correction_terms from '
                             'pbc_cubic, galvani, dielectric_legacy.')
        if self.correction_environment != 'homogeneous_bulk':
            raise ValueError('These corrections require explicitly declared '
                             'homogeneous_bulk; interfaces, slabs and binding '
                             'pockets are not supported by this scheme.')
        for name in ('correction_solvent', 'ion_reference_convention'):
            value = getattr(self, name)
            if not isinstance(value, str) or not value.strip():
                raise ValueError(f'{name} must be explicitly specified.')
        if self.charge_model not in ('unscaled', 'effective'):
            raise ValueError("Set charge_model to 'unscaled' or 'effective'.")
        if (self.formal_charge is None or not np.isfinite(self.formal_charge)
                or not np.isclose(self.formal_charge, round(self.formal_charge),
                                  rtol=0.0, atol=1e-8)):
            raise ValueError('formal_charge must be an explicit integer chemical charge.')
        if 'pbc_cubic' in terms or 'dielectric_legacy' in terms:
            if (self.solvent_dielectric is None
                    or not np.isfinite(self.solvent_dielectric)
                    or self.solvent_dielectric < 1.0):
                raise ValueError('Specify solvent_dielectric >= 1 for this scheme.')
        if 'galvani' in terms:
            if (self.galvani_coefficient is None
                    or not np.isfinite(self.galvani_coefficient)):
                raise ValueError('Specify a signed galvani_coefficient in '
                                 'kJ/mol per unit charge; there is no solvent fallback.')
            if self.galvani_charge_basis not in ('formal', 'interaction'):
                raise ValueError("Select galvani_charge_basis: 'formal' or 'interaction'.")
        if 'dielectric_legacy' in terms:
            if self.charge_model != 'unscaled':
                raise ValueError('The legacy dielectric correction is not '
                                 'supported for effective/ECC charge models.')
            if (self.target_dielectric is None
                    or not np.isfinite(self.target_dielectric)
                    or self.target_dielectric < 1.0):
                raise ValueError('Specify target_dielectric >= 1 for the legacy term.')

    @staticmethod
    def _cubic_box_length_nm(box):
        """Accept only axis-aligned cubic cells; no spherical/slab fallback."""
        if box is None:
            raise ValueError('pbc_cubic requires periodic box vectors.')
        if hasattr(box, 'value_in_unit'):
            values = np.asarray(box.value_in_unit(unit.nanometers), dtype=float)
        else:
            values = np.asarray([
                v.value_in_unit(unit.nanometers) if hasattr(v, 'value_in_unit') else v
                for v in box], dtype=float)
        if values.shape != (3, 3) or not np.all(np.isfinite(values)):
            raise ValueError('Invalid periodic box vectors.')
        lengths = np.diag(values)
        if (np.any(lengths <= 0)
                or not np.allclose(values, np.diag(lengths), rtol=0, atol=1e-8)
                or not np.allclose(lengths, lengths[0], rtol=1e-6, atol=1e-8)):
            raise ValueError('pbc_cubic supports only axis-aligned cubic boxes, '
                             'not orthorhombic or triclinic cells or slabs.')
        return float(lengths[0])

    def _record_interaction_charge(self, system, alchemical_region):
        """Read the unscaled endpoint FF charge without changing any charges."""
        if not alchemical_region:
            raise ValueError('The alchemical region is empty.')
        forces = [f for f in system.getForces() if isinstance(f, mm.NonbondedForce)]
        self.interaction_charge = float(sum(
            f.getParticleParameters(i)[0].value_in_unit(unit.elementary_charge)
            for f in forces for i in alchemical_region))
        self.ion_charge = self.interaction_charge
        self.ostream.print_info(f"Solute force-field net charge: {self.interaction_charge:+.8f} e")
        if not np.isfinite(self.interaction_charge):
            raise ValueError('The input force-field charge is not finite.')
        if self._active_ion_mode == 'auto':
            if len(forces) != 1 or forces[0].getNumParticleParameterOffsets():
                raise ValueError('Auto requires explicit charges in one NonbondedForce '
                                 'without particle offsets; use none/custom.')
            self._resolve_auto_ion_corrections()
        if not self.apply_ion_corrections:
            return
        if len(forces) != 1 or forces[0].getNumParticleParameterOffsets():
            raise ValueError('Ion corrections require one NonbondedForce with '
                             'explicit particle charges and no parameter offsets.')
        if self.charge_model == 'unscaled' and not np.isclose(
                self.interaction_charge, self.formal_charge, rtol=0, atol=1e-5):
            raise ValueError('The FF charge differs from formal_charge; declare '
                             "charge_model='effective' for scaled interaction charges.")
        if 'pbc_cubic' in self.ion_correction_terms:
            if forces[0].getNonbondedMethod() not in (mm.NonbondedForce.PME,
                                                     mm.NonbondedForce.Ewald):
                raise ValueError('pbc_cubic requires 3D PME or Ewald electrostatics.')
            box = system.getDefaultPeriodicBoxVectors()
            self._cubic_box_length_nm(box)
            self._ion_correction_box_nm = box
            if len(alchemical_region) == 1:
                if (self.ion_radius_nm is None or not np.isfinite(self.ion_radius_nm)
                        or self.ion_radius_nm <= 0):
                    raise ValueError('Specify a positive ion_radius_nm for a monatomic ion.')
            elif self.solute_volume_nm3 is not None:
                if not np.isfinite(self.solute_volume_nm3) or self.solute_volume_nm3 <= 0:
                    raise ValueError('solute_volume_nm3 must be positive and finite.')

    def _calculate_ion_corrections(self):
        """Return only the selected additive terms in kJ/mol.

        pbc_cubic retains the existing spherical-cavity approximation and uses
        the input reference cell, not an ensemble-averaged NPT correction.
        It uses the actual interaction charge of the simulated Hamiltonian.
        """
        self._validate_ion_correction_settings()
        if not self.apply_ion_corrections:
            return 0.0, 0.0, 0.0
        if self.interaction_charge is None:
            raise ValueError('Interaction charge has not been evaluated.')
        terms = self.ion_correction_terms
        dg_pbc = dg_phi = dg_d = 0.0
        if 'pbc_cubic' in terms:
            L_input = self._cubic_box_length_nm(self._ion_correction_box_nm)
            if self.ion_correction_box == 'average':
                average = self._average_box
                if average is None or not average.get('cubic', False):
                    raise ValueError('No cubic sampled box from the coupled solvated window '
                                     '(stage 1, lambda = 0); use ion_correction_box="input".')
                L = average['L_nm']
            else:
                L = L_input
            self._pbc_box_used_nm = L
            self.ostream.print_info(f'pbc_cubic box: {self.ion_correction_box} L = {L:.4f} nm '
                                    f'(input box L = {L_input:.4f} nm)')
            region, _ = self._get_alchemical_region(self.vacuum_topology)
            if len(region) == 1:
                r_eff = float(self.ion_radius_nm)
            else:
                volume = self.solute_volume_nm3
                if volume is None:
                    volume = self._estimate_solute_volume(
                        self.vacuum_topology, self.vacuum_positions)
                if not np.isfinite(volume) or volume <= 0:
                    raise ValueError('Invalid solute cavity volume.')
                r_eff = (3.0 * volume / (4.0 * math.pi)) ** (1.0 / 3.0)
            if not 0 < r_eff < L / 2:
                raise ValueError('The spherical cavity radius must be smaller than half the box.')
            eps_s = float(self.solvent_dielectric)
            term1 = 2.837297 / eps_s
            term2 = (1.0 - 1.0 / eps_s) * (
                (4.0 * math.pi / 3.0) * (r_eff / L)**2
                - (16.0 * math.pi**2 / 45.0) * (r_eff / L)**5)
            dg_pbc = (138.935456 * self.interaction_charge**2 / (2.0 * L)
                      * (term1 + term2))
        if 'galvani' in terms:
            q = (self.formal_charge if self.galvani_charge_basis == 'formal'
                 else self.interaction_charge)
            dg_phi = float(q) * float(self.galvani_coefficient)
        if 'dielectric_legacy' in terms:
            dg_d = self._calculate_dielectric_correction(float(self.solvent_dielectric))
        return dg_pbc, dg_phi, dg_d

    def _run_stages(self):
        """Resolve correction settings per call and restore user options afterward."""
        self._reset_calculation_state()
        self._validate_analysis_settings()
        self._validate_md_settings()
        # kT must follow the temperature set by the user after __init__;
        # a stale kT would rescale all reduced potentials passed to MBAR.
        self.kT = self._compute_kT()
        if self.ion_corrections not in ('auto', 'none', 'custom'):
            raise ValueError("ion_corrections must be 'auto', 'none', or 'custom'.")
        mode = self.ion_corrections
        if self.apply_ion_corrections is not None:
            if not isinstance(self.apply_ion_corrections, bool):
                raise ValueError('The legacy apply_ion_corrections override must be bool or None.')
            mode = 'custom' if self.apply_ion_corrections else 'none'
        names = ('apply_ion_corrections', 'ion_correction_terms',
                 'correction_environment', 'correction_solvent',
                 'ion_reference_convention', 'charge_model', 'formal_charge',
                 'solvent_dielectric', 'target_dielectric',
                 'galvani_coefficient', 'galvani_charge_basis')
        original = {name: getattr(self, name) for name in names}
        self._active_ion_mode = mode
        self._ion_preset_name = None
        try:
            self.apply_ion_corrections = mode == 'custom'
            if mode == 'none':
                self.ion_correction_terms = ()
            return self._run_stages_impl()
        finally:
            for name, value in original.items():
                setattr(self, name, value)

    def _run_stages_impl(self):
        self._validate_ion_correction_settings()
        header = "VeloxChem/OpenMM Solvation Free Energy Calculation"
        self.ostream.print_header(header)
        self.ostream.print_header("="*len(header))
        self.ostream.print_blank()
        self.ostream.flush()

        # Print the simulation parameters
        self.ostream.print_line("Simulation parameters:")
        self.ostream.print_line("-"*len("Simulation parameters:"))
        self.ostream.print_blank()
        self.ostream.print_line(f"Temperature: {self.temperature}")
        self.ostream.print_line(f"Pressure: {self.pressure}")
        self.ostream.print_line(f"Timestep: {self.timestep}")
        self.ostream.print_line(f"Non-Bonded Cutoff: {self.cutoff}")
        self.ostream.print_line(f"Integrator: {self.integrator} (friction {self.friction})")
        self.ostream.print_line(f"Hydrogen mass repartitioning: "
                                f"{self.hydrogen_mass if self.hydrogen_mass is not None else 'off'}")
        self.ostream.print_line(f"Sequential seeding of lambda windows: {self.sequential_seeding}")
        self.ostream.print_line(f"Linear coupling in stages 1 and 3: {self.linear_coupling}")
        self.ostream.print_line(f"Alchemical scheme: {self.alchemical_scheme}")
        if self._decoupling:
            self.ostream.print_line("Intramolecular interactions kept at full strength; GSC between "
                                    "solute and environment only (no vacuum leg)")
        else:
            self.ostream.print_line(f"Solute-solute GSC kept in the decoupled state: "
                                    f"{self.keep_intramolecular_gsc}"
                                    + (" (stage 4 not needed)" if self.keep_intramolecular_gsc else ""))
        self.ostream.print_line(f"Platform: {self.platform or 'fastest available'} "
                                f"(precision {self.precision}); vacuum platform: "
                                f"{self.vacuum_platform or self.platform or 'fastest available'}")
        self.ostream.print_line("Thermodynamic Ensemble: NPT (solvated, if a barostat is present)")
        simulation_time = self.num_steps * self.timestep.value_in_unit(unit.nanoseconds)
        self.ostream.print_blank()
        self.ostream.print_line("Alchemical Parameters:")
        self.ostream.print_line("-"*len("Alchemical Parameters:"))
        self.ostream.print_line(f"Number of equilibration steps: {self.num_equil_steps}")
        self.ostream.print_line(f"Number of steps per lambda simulation: {self.num_steps}")
        self.ostream.print_line(f"Number of snapshots per lambda simulation: {self.num_snapshots}")
        self.ostream.print_line(f"Simulation time per lambda: {simulation_time:.2f} ns")
        self.ostream.print_line(f"Lambdas in Stage 1: {self.lambdas_stage1}")
        self.ostream.print_line(f"Lambdas in Stage 2: {self.lambdas_stage2}")
        if self._decoupling:
            self.ostream.print_line("Lambdas in Stages 3 and 4: not run (decoupling, gas-phase reference)")
        else:
            self.ostream.print_line(f"Lambdas in Stage 3: {self.lambdas_stage3}")
        if self._decoupling:
            pass
        elif self.keep_intramolecular_gsc:
            self.ostream.print_line("Lambdas in Stage 4: not run (identical end states)")
        else:
            self.ostream.print_line(f"Lambdas in Stage 4: {self.lambdas_stage4}")
        self.ostream.print_blank()

        # Run the simulations
        self.ostream.print_info("Starting solvated simulation (Stage 1)...\n")
        delta_f_1, free_en_s1 = self._run_lambda_simulations(stage=1)
        self.ostream.flush()

        if self.keep_intramolecular_gsc or self._decoupling:
            self.ostream.print_info("Removing solute-environment GSC potential (Stage 2)...\n")
        else:
            self.ostream.print_info("Removing GSC potential (Stage 2)...\n")
        delta_f_2, free_en_s2 = self._run_lambda_simulations(stage=2)
        self.ostream.flush()

        if self._decoupling:
            # Decoupled state = environment x gas-phase solute: no vacuum leg.
            # The vacuum endpoint is still built (not simulated) to validate the
            # solute mapping/charges and for the solute cavity volume.
            self.stage = 3
            _, self.vacuum_topology, self.vacuum_positions = self._create_vacuum_systems([0.0])
            self.ostream.print_info("Stages 3 and 4 not run: alchemical_scheme='decoupling' keeps all "
                                    "intramolecular interactions, the decoupled state is the gas-phase "
                                    "solute (Stage 3 = Stage 4 = 0 exactly).\n")
            delta_f_3, free_en_s3 = self._skipped_stage_result(), 0.0
            delta_f_4, free_en_s4 = self._skipped_stage_result(), 0.0
            self.u_kn_matrices.extend([None, None])
        else:
            self.ostream.print_info("Starting vacuum simulation (Stage 3)...\n")
            delta_f_3, free_en_s3 = self._run_lambda_simulations(stage=3, vacuum=True)
            self.ostream.flush()

        if self._decoupling:
            pass
        elif self.keep_intramolecular_gsc:
            # The last stage-3 state (solute-solute GSC on, nothing else
            # nonbonded) is already the decoupled reference: Stage 4 = 0.
            self.ostream.print_info("Stage 4 not run: with keep_intramolecular_gsc the stage-3 "
                                    "end state is the decoupled reference (Stage 4 = 0 exactly).\n")
            delta_f_4, free_en_s4 = self._skipped_stage_result(), 0.0
            self.u_kn_matrices.append(None)
        else:
            self.ostream.print_info("Removing GSC potential in vacuum (Stage 4)...\n")
            delta_f_4, free_en_s4 = self._run_lambda_simulations(stage=4, vacuum=True)
        self.ostream.flush()

        # Calculate the final free energy
        final_free_energy = free_en_s1 + free_en_s2 - (free_en_s3 + free_en_s4)

        raw_free_energy = final_free_energy
        self.ostream.print_line(f"Raw FEP free energy: {raw_free_energy:.4f} kJ/mol")
        self.ostream.print_line(f"Ion correction mode: {self._active_ion_mode}; "
                               f"preset: {self._ion_preset_name}")
        if self.apply_ion_corrections:
            self.dG_PBC, self.dG_phi, self.dG_D = self._calculate_ion_corrections()
            final_free_energy += self.dG_PBC + self.dG_phi + self.dG_D
            self.ostream.print_line(
                f"Selected ion corrections: {self.ion_correction_terms}; "
                f"solvent: {self.correction_solvent}; "
                f"reference: {self.ion_reference_convention}")
            self.ostream.print_line(
                f"Formal charge: {self.formal_charge}; "
                f"interaction charge: {self.interaction_charge:.8f}")
            self.ostream.print_line(
                f"Corrections (kJ/mol): PBC={self.dG_PBC:.4f}, "
                f"Galvani={self.dG_phi:.4f}, dielectric={self.dG_D:.4f}")

        self.ostream.print_blank()
        self.ostream.print_line(f"Final free energy: {final_free_energy:.4f} kJ/mol")
        self.ostream.flush()

        self.delta_f = [delta_f_1, delta_f_2, delta_f_3, delta_f_4]

        delta_f = {f"Stage {i+1}": {'Delta_f': (self.delta_f[i]['Delta_f'][-1, 0] * self.kT).value_in_unit(unit.kilojoule_per_mole),
                                    'Uncertainty': (self.delta_f[i]['dDelta_f'][-1, 0] * self.kT).value_in_unit(unit.kilojoule_per_mole)}
                   for i in range(4)}
        for i in range(4):
            delta_f[f'Stage {i+1}']['sampling_diagnostics'] = self.delta_f[i]['sampling_diagnostics']
        total_variance = sum(delta_f[f"Stage {i+1}"]['Uncertainty']**2 for i in range(4))

        delta_f['dG_PBC'] = self.dG_PBC
        delta_f['dG_phi'] = self.dG_phi
        delta_f['dG_D'] = self.dG_D

        delta_f['raw_free_energy'] = raw_free_energy
        delta_f['interaction_charge'] = self.interaction_charge
        delta_f['formal_charge'] = self.formal_charge
        delta_f['ion_correction_terms'] = (list(self.ion_correction_terms)
                                           if self.apply_ion_corrections else [])
        delta_f['ion_corrections'] = self._active_ion_mode
        delta_f['ion_correction_metadata'] = {
            'preset': self._ion_preset_name,
            'preset_note': self._ion_preset_note,
            'provisional': bool(self._ion_preset_note and
                                self._ion_preset_note.startswith('Provisional')),
            'environment': self.correction_environment,
            'solvent': self.correction_solvent or self._auto_bulk_solvent,
            'reference_convention': self.ion_reference_convention,
            'charge_model': self.charge_model,
            'galvani_charge_basis': self.galvani_charge_basis,
            'galvani_coefficient_kj_mol_per_e': self.galvani_coefficient,
            'solvent_dielectric': self.solvent_dielectric,
            'target_dielectric': self.target_dielectric,
            'includes_ecc_electronic_solvation': False,
            'pbc_box': {
                'mode': self.ion_correction_box,
                'used_L_nm': self._pbc_box_used_nm,
                'input_L_nm': (None if self._ion_correction_box_nm is None
                               else self._cubic_box_length_nm(self._ion_correction_box_nm)),
                'sampled': self._average_box,
            },
        }
        delta_f['free_energy'] = final_free_energy
        delta_f['uncertainty'] = float(np.sqrt(total_variance))
        delta_f['md_protocol'] = {
            'temperature_K': float(self.temperature),
            'kT_kj_mol': float(self.kT.value_in_unit(unit.kilojoule_per_mole)),
            'integrator': self.integrator,
            'friction_per_ps': float(self.friction.value_in_unit(unit.picoseconds**-1)),
            'timestep_fs': float(self.timestep.value_in_unit(unit.femtoseconds)),
            'hydrogen_mass_amu': (None if self.hydrogen_mass is None
                                  else self._hydrogen_mass_amu()),
            'sequential_seeding': bool(self.sequential_seeding),
            'precision': self.precision,
            'solvated_barostat': getattr(self, '_has_barostat', None),
            'linear_coupling': bool(self.linear_coupling),
            'keep_intramolecular_gsc': bool(self.keep_intramolecular_gsc),
            'alchemical_scheme': self.alchemical_scheme,
            'gsc_placeholder_sigma_nm': (None if self.gsc_placeholder_sigma is None else
                                         float(self.gsc_placeholder_sigma.value_in_unit(unit.nanometers))),
        }

        self.final_free_energy = final_free_energy

        return delta_f

    @staticmethod
    def _skipped_stage_result():
        """Stage result for a stage that is zero by construction (not simulated)."""
        return {'Delta_f': np.zeros((1, 1)), 'dDelta_f': np.zeros((1, 1)),
                'sampling_diagnostics': {
                    'method': 'not_run', 'equilibration': [], 'windows': [],
                    'messages': ['Not simulated: end states identical (keep_intramolecular_gsc=True '
                                 "or alchemical_scheme='decoupling'); free energy is exactly 0."]}}

    def _run_lambda_simulations(self, stage, vacuum=False):
        """
        Runs simulations for a given stage.

        :param stage:
            The stage number. 1 for solvated simulation, 2 for removing GSC potential, 3 for vacuum simulation.
        :param vacuum:
            If True, run the simulation in vacuum.
        """

        total_sim_time = 0  # Track total simulation time
        total_ns_simulated = 0  # Track total simulated time in ns
        self.stage = stage

        if stage == 1:
            lambdas = self.lambdas_stage1
        elif stage == 2:
            lambdas = self.lambdas_stage2
        elif stage == 3:
            lambdas = self.lambdas_stage3
        elif stage == 4:
            lambdas = self.lambdas_stage4

        self.ostream.print_info(f"Generating systems for stage {stage}")
        self.ostream.flush()

        # Create systems, topology and initial positions
        if vacuum:
            systems, topology, positions = self._create_vacuum_systems(lambdas)
            self.vacuum_topology = topology
            self.vacuum_positions = positions
        else:
            systems, topology, positions = self._create_solvated_systems(lambdas)
            self.topology = topology
            self.positions = positions

        self.snapshots = []
        self.N_k = []
        self._equilibration_diagnostics = []

        phase = 'vacuum' if vacuum else 'solvated'
        rounded = [round(float(x), 6) for x in lambdas]
        coupled_idx = rounded.index(0.0) if 0.0 in rounded else 0
        for lam_idx, lam in enumerate(lambdas):
            lam = round(lam, 2)
            self.ostream.print_info(f"Running lambda = {lam}, stage = {stage}...")
            self.ostream.flush()

            # Windows of the same phase share topology; stage 2 starts from the
            # last stage-1 window, which has the identical Hamiltonian.
            seed = None
            if (self.sequential_seeding and self._seed_state is not None
                    and self._seed_state[0] == phase
                    and len(self._seed_state[1]) == systems[lam_idx].getNumParticles()):
                seed = self._seed_state[1:]

            ns_simulated, sim_time, final_state = self._run_simulation(
                systems[lam_idx], topology, positions, lam, seed=seed, vacuum=vacuum)
            self._seed_state = (phase,) + final_state
            if stage == 1 and not vacuum and lam_idx == coupled_idx:
                # Fully coupled solvated state: its sampled box enters the PBC correction.
                self._record_average_box(self.snapshots[-self.N_k[-1]:], lam)
            total_ns_simulated += ns_simulated
            total_sim_time += sim_time

        # Print total performance metrics for the stage
        overall_ns_per_hour = (total_ns_simulated / total_sim_time) * 3600
        self.ostream.print_info(f"Total simulation time for stage {stage}: {total_sim_time:.2f} s. "
                                f"Total simulated time: {total_ns_simulated:.2f} ns.")
        self.ostream.print_info(f"Overall performance for stage {stage}: {overall_ns_per_hour:.2f} ns/hour")
        self.ostream.flush()

        self.ostream.print_info(f"Recalculating energies for stage {stage}...")
        self.ostream.flush()
        recalculation_start = time.time()

        u_kn = self._recalculate_energies(systems, topology, vacuum=vacuum, lambdas=lambdas)
        self.u_kn_matrices.append(u_kn)

        recalculation_end = time.time()
        recalculation_time = recalculation_end - recalculation_start

        self.ostream.print_info(f"Energy recalculation for stage {stage} took {recalculation_time:.2f} seconds.")
        self.ostream.flush()

        self.ostream.print_info(f"Calculating the free energy with MBAR for stage {stage}...")
        self.ostream.flush()

        delta_f = self._calculate_free_energy(u_kn)
        free_energy = (delta_f['Delta_f'][-1, 0] * self.kT).value_in_unit(unit.kilojoule_per_mole)

        self.ostream.print_line(f"Free energy for stage {stage}: {free_energy:.4f} +/- {(delta_f['dDelta_f'][-1, 0] * self.kT).value_in_unit(unit.kilojoule_per_mole):.4f} kJ/mol")
        self.ostream.flush()

        return delta_f, free_energy

    def _create_solvated_systems(self, lambdas):
        """
        Create solvated system with alchemical scaling.
        """
        # Reading from pre-built OpenMM system XML
        if self.solvent_name == 'omm_system_xml':
            pdb = app.PDBFile(self.system_pdb)
            topology = pdb.topology
            positions = pdb.positions
            with open(self.system_xml, 'r') as f:
                base_system = mm.XmlSerializer.deserialize(f.read())

        elif self.solvent_name == 'omm_memory':
            topology = self.topology
            positions = self.positions
            base_system = self.base_system

        elif self.solvent_name == 'omm_files':
            pdb = app.PDBFile(self.system_pdb)
            initial_system_ff = app.ForceField(self.solute_xml, *self.other_xml_files)
            topology = pdb.topology
            positions = pdb.positions

        elif self.solvent_name == 'gro_files':
            gro = app.GromacsGroFile(self.system_gro)
            initial_system_ff = app.GromacsTopFile(self.system_top,
                                                   periodicBoxVectors=gro.getPeriodicBoxVectors())
            topology = initial_system_ff.topology
            positions = gro.positions

        # From molecule and ffgenerator
        else:
            if self.solvent_name != 'itself':
                pdb = app.PDBFile('system.pdb')
                initial_system_ff = app.ForceField('solute.xml', 'solvent_1.xml')
                topology = pdb.topology
                positions = pdb.positions
            else:
                pdb = app.PDBFile('liquid.pdb')
                initial_system_ff = app.ForceField('liquid.xml')
                topology = pdb.topology
                positions = pdb.positions

        # createSystem from a top file does not require the topology as an argument
        if self.solvent_name == 'gro_files':
            sys_arguments = [self.nonbondedMethod, self.cutoff, self.constraints]
        else:
            sys_arguments = [topology, self.nonbondedMethod, self.cutoff, self.constraints]

        # Build (or deserialize) the unscaled endpoint once per stage and clone
        # it for every lambda; createSystem is expensive for large systems.
        if self.solvent_name in ['omm_system_xml', 'omm_memory']:
            template_system = copy.deepcopy(base_system)
        else:
            template_system = initial_system_ff.createSystem(*sys_arguments)
            template_system.addForce(mm.MonteCarloBarostat(self.pressure, self.temperature * unit.kelvin))
        self._report_barostat(template_system)

        # Define alchemical regions
        alchemical_region, chemical_region = self._get_alchemical_region(topology)
        alchemical_set = set(alchemical_region)
        self._validate_endpoint(template_system, topology, positions, alchemical_region)

        if self.interaction_charge is None:
            self._record_interaction_charge(template_system, alchemical_region)

        self._repartition_hydrogen_masses(template_system, topology)
        self._report_gsc_placeholders(template_system)

        solvated_systems = []
        for lambda_val in lambdas:
            solvated_system = copy.deepcopy(template_system)

            # Custom nonbonded force (Gaussian Softcore)
            gsc_force = self._get_gaussian_softcore_force(lambda_val)

            # Scaling the nonbonded interactions
            original_params = {}
            existing_pairs = set()
            gsc_pairs = None
            intra_gsc = None
            for force in solvated_system.getForces():
                if isinstance(force, mm.NonbondedForce):
                    gsc_force.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
                    gsc_force.setCutoffDistance(self.cutoff)

                    for i in range(force.getNumParticles()):
                        charge, sigma, epsilon = force.getParticleParameters(i)
                        gsc_force.addParticle([self._gsc_sigma(sigma, epsilon)])

                        if i in alchemical_set:
                            original_params[i] = (charge, sigma, epsilon)
                            if self.stage == 1:
                                force.setParticleParameters(i, charge * (1 - lambda_val), sigma,
                                                            epsilon * self._epsilon_scale(lambda_val))
                            elif self.stage == 2:
                                force.setParticleParameters(i, 0, 0, 0)

                    # Handle exceptions (interaction exclusions and 1-4 scaling)
                    for i in range(force.getNumExceptions()):
                        p1, p2, ch, si, ep = force.getExceptionParameters(i)
                        gsc_force.addExclusion(p1, p2)
                        existing_pairs.add((min(p1, p2), max(p1, p2)))
                        # Scale the base NonbondedForce exception if it involves the alchemical region
                        # (decoupling: intramolecular exceptions stay unscaled)
                        if self._decoupling and p1 in alchemical_set and p2 in alchemical_set:
                            pass
                        elif p1 in alchemical_set or p2 in alchemical_set:
                            if self.stage == 1:
                                force.setExceptionParameters(i, p1, p2, ch * (1 - lambda_val), si, ep * (1 - lambda_val))
                            elif self.stage == 2:
                                force.setExceptionParameters(i, p1, p2, 0, si, 0)

                    if self._decoupling:
                        self._add_intramolecular_exceptions(force, gsc_force, alchemical_region,
                                                            original_params, existing_pairs)
                    elif self.stage == 1 and self.linear_coupling:
                        gsc_pairs = self._add_linear_pair_terms(force, gsc_force, alchemical_region,
                                                                lambda_val, original_params, existing_pairs,
                                                                periodic=True)
                    elif self.stage == 2 and self.keep_intramolecular_gsc:
                        # Solute-solute GSC at full strength, independent of lambda:
                        # zero exceptions + GSC exclusions for the non-excluded solute
                        # pairs, their GSC carried by a lambda-free pair force.
                        intra_gsc = self._add_linear_pair_terms(force, gsc_force, alchemical_region,
                                                                1.0, original_params, existing_pairs,
                                                                periodic=True, fixed_gsc=True)

            # Add interaction group for GSC potential
                elif isinstance(force, mm.CustomNonbondedForce):
                    for i in range(force.getNumGlobalParameters()):
                        if force.getGlobalParameterName(i) == 'lambda':
                            current_lambda = lambda_val if self.stage == 1 else 1.0
                            force.setGlobalParameterDefaultValue(i, current_lambda)

            gsc_force.addInteractionGroup(alchemical_region, chemical_region)
            if intra_gsc is None and not self._decoupling:
                gsc_force.addInteractionGroup(alchemical_region, alchemical_region)  # To avoid self-collapse in flexible molecules
            gsc_group = self._free_force_group(solvated_system)
            gsc_force.setForceGroup(gsc_group)
            solvated_system.addForce(gsc_force)
            if gsc_pairs is not None:
                gsc_pairs.setForceGroup(gsc_group)
                solvated_system.addForce(gsc_pairs)
            if intra_gsc is not None:
                # Lambda-independent: stays with the other forces (default group 0).
                solvated_system.addForce(intra_gsc)
            solvated_systems.append(solvated_system)

        return solvated_systems, topology, positions

    def _create_vacuum_systems(self, lambdas):
        """
        Create vacuum system with alchemical scaling.
        """
        if self.solvent_name == 'omm_memory':
            topology = self.vacuum_topology
            positions = self.vacuum_positions
            base_vacuum_system = self.vacuum_system

        elif self.solvent_name in ['omm_files', 'omm_system_xml']:
            pdb = app.PDBFile(self.solute_pdb)
            forcefield_solute = app.ForceField(self.solute_xml)
            topology = pdb.topology
            positions = pdb.positions

        elif self.solvent_name == 'gro_files':
            pdb = app.GromacsGroFile(self.solute_gro)
            forcefield_solute = app.GromacsTopFile(self.solute_top)
            topology = forcefield_solute.topology
            positions = pdb.positions

        # From molecule and ffgenerator
        else:
            self.solute_ff.write_openmm_files('solute_vacuum', 'MOL')
            pdb = app.PDBFile('solute_vacuum.pdb')
            forcefield_solute = app.ForceField('solute_vacuum.xml')
            topology = pdb.topology
            positions = pdb.positions

        if self.solvent_name == 'gro_files':
            template_system = forcefield_solute.createSystem(nonbondedMethod=app.NoCutoff, constraints=self.constraints)
        elif self.solvent_name == 'omm_memory':
            template_system = copy.deepcopy(base_vacuum_system)
        else:
            template_system = forcefield_solute.createSystem(topology, nonbondedMethod=app.NoCutoff, constraints=self.constraints)

        alchemical_region, chemical_region = self._get_alchemical_region(topology)
        self._validate_endpoint(template_system, topology, positions, alchemical_region, vacuum=True)
        self._repartition_hydrogen_masses(template_system, topology)
        self._report_gsc_placeholders(template_system)

        vacuum_systems = []
        for lambda_val in lambdas:
            vacuum_system = copy.deepcopy(template_system)

            gsc_force = self._get_gaussian_softcore_force(lambda_val)

            original_params = {}
            existing_pairs = set()
            gsc_pairs = None
            for force in vacuum_system.getForces():
                if isinstance(force, mm.NonbondedForce):
                    for i in range(force.getNumParticles()):
                        charge, sigma, epsilon = force.getParticleParameters(i)
                        gsc_force.addParticle([self._gsc_sigma(sigma, epsilon)])
                        original_params[i] = (charge, sigma, epsilon)

                        if self.stage == 3:
                            force.setParticleParameters(i, charge * (1 - lambda_val), sigma,
                                                        epsilon * self._epsilon_scale(lambda_val))
                        elif self.stage == 4:
                            force.setParticleParameters(i, 0, 0, 0)

                    # Handle exceptions (interaction exclusions and 1-4 scaling)
                    for i in range(force.getNumExceptions()):
                        p1, p2, ch, si, ep = force.getExceptionParameters(i)
                        gsc_force.addExclusion(p1, p2)
                        existing_pairs.add((min(p1, p2), max(p1, p2)))
                        # Scale the base NonbondedForce exception for the vacuum stages

                        if self.stage == 3:
                            force.setExceptionParameters(i, p1, p2, ch * (1 - lambda_val), si, ep * (1 - lambda_val))
                        elif self.stage == 4:
                            force.setExceptionParameters(i, p1, p2, 0, si, 0)

                    if self.stage == 3 and self.linear_coupling:
                        gsc_pairs = self._add_linear_pair_terms(force, gsc_force, alchemical_region,
                                                                lambda_val, original_params, existing_pairs,
                                                                periodic=False)

                elif isinstance(force, mm.CustomNonbondedForce):
                    for i in range(force.getNumGlobalParameters()):
                        if force.getGlobalParameterName(i) == 'lambda':
                            current_lambda = lambda_val if self.stage == 3 else 1.0
                            force.setGlobalParameterDefaultValue(i, current_lambda)

            gsc_group = self._free_force_group(vacuum_system)
            gsc_force.setForceGroup(gsc_group)
            vacuum_system.addForce(gsc_force)
            if gsc_pairs is not None:
                gsc_pairs.setForceGroup(gsc_group)
                vacuum_system.addForce(gsc_pairs)
            vacuum_systems.append(vacuum_system)

        return vacuum_systems, topology, positions

    def _epsilon_scale(self, lambda_val):
        """Solute epsilon factor in stages 1 and 3.

        With Lorentz-Berthelot mixing the solute-environment well depth is
        sqrt(f * eps_i * eps_j); f = (1-lambda)^2 makes it linear in lambda.
        """
        s = 1.0 - lambda_val
        return s * s if self.linear_coupling else s

    def _add_linear_pair_terms(self, force, gsc_force, region, lambda_val, params, existing,
                               periodic, fixed_gsc=False):
        """Explicit, linearly scaled terms for non-excluded intramolecular pairs.

        Each such solute pair becomes a NonbondedForce exception with
        chargeProd = (1-lambda) q_a q_b and epsilon = (1-lambda) sqrt(eps_a eps_b)
        (unscaled parameters, Lorentz-Berthelot mixing). With PME the explicit
        term plus the reciprocal-space part reproduces the same Ewald pair
        interaction; intramolecular LJ is no longer truncated at the cutoff,
        consistent with the NoCutoff vacuum leg.

        OpenMM requires identical exclusions in all nonbonded forces, so these
        pairs are also excluded from the GSC CustomNonbondedForce and their GSC
        term (same expression and parameters) is returned as a CustomBondForce.
        """
        scale = 1.0 - lambda_val
        pairs = self._get_gaussian_softcore_pair_force(lambda_val, periodic, fixed=fixed_gsc)
        for a, b in itertools.combinations(sorted(region), 2):
            if (a, b) in existing:
                continue
            qa, sa, ea = params[a]
            qb, sb, eb = params[b]
            force.addException(a, b, qa * qb * scale, 0.5 * (sa + sb), (ea * eb).sqrt() * scale)
            gsc_force.addExclusion(a, b)
            sigma_pair = 0.5 * (self._gsc_sigma(sa, ea) + self._gsc_sigma(sb, eb))
            pairs.addBond(a, b, [sigma_pair])
        return pairs if pairs.getNumBonds() > 0 else None

    @property
    def _decoupling(self):
        return self.alchemical_scheme == 'decoupling'

    def _add_intramolecular_exceptions(self, force, gsc_force, region, params, existing):
        """Decoupling: every non-excluded solute pair becomes an unscaled exception.

        chargeProd = q_a q_b, sigma = (sigma_a + sigma_b)/2, epsilon =
        sqrt(eps_a eps_b) from the unscaled parameters (Lorentz-Berthelot), so
        the intramolecular Coulomb and LJ energies do not depend on lambda (and
        LJ is not truncated, as in the gas phase). The pairs are excluded from
        the GSC force (identical exclusions); no solute-solute GSC is added.
        """
        for a, b in itertools.combinations(sorted(region), 2):
            if (a, b) in existing:
                continue
            qa, sa, ea = params[a]
            qb, sb, eb = params[b]
            force.addException(a, b, qa * qb, 0.5 * (sa + sb), (ea * eb).sqrt())
            gsc_force.addExclusion(a, b)

    def _get_gaussian_softcore_pair_force(self, lambda_val, periodic, fixed=False):
        """GSC term for explicit solute pairs (same expression as the GSC force).

        fixed=True gives the full-strength term without a lambda parameter
        (solute-solute GSC kept in the decoupled state, stage 2).
        """
        if fixed:
            energy = 'alpha * exp(-beta * (r/Rmin)^x); Rmin = 1.122462 * sigma;'
        else:
            energy = 'lambda * alpha * exp(-beta * (r/Rmin)^x); Rmin = 1.122462 * sigma;'
        pairs = mm.CustomBondForce(energy)
        if not fixed:
            pairs.addGlobalParameter('lambda', lambda_val)
        pairs.addGlobalParameter('alpha', self.alpha)
        pairs.addGlobalParameter('beta', self.beta)
        pairs.addGlobalParameter('x', self.x)
        pairs.addPerBondParameter('sigma')
        pairs.setUsesPeriodicBoundaryConditions(periodic)
        return pairs

    def _gsc_sigma(self, sigma, epsilon):
        """Sigma used by the GSC term for one particle."""
        if (self.gsc_placeholder_sigma is not None
                and epsilon.value_in_unit(unit.kilojoule_per_mole) == 0.0):
            return self.gsc_placeholder_sigma
        return sigma

    def _report_gsc_placeholders(self, system):
        """Count epsilon = 0 particles whose sigma the GSC replaces."""
        if self.gsc_placeholder_sigma is None:
            return
        target = self.gsc_placeholder_sigma.value_in_unit(unit.nanometers)
        n_zero, changed = 0, {}
        for force in system.getForces():
            if isinstance(force, mm.NonbondedForce):
                for i in range(force.getNumParticles()):
                    _, sigma, epsilon = force.getParticleParameters(i)
                    if epsilon.value_in_unit(unit.kilojoule_per_mole) == 0.0:
                        n_zero += 1
                        value = round(sigma.value_in_unit(unit.nanometers), 6)
                        if not np.isclose(value, target):
                            changed[value] = changed.get(value, 0) + 1
        self._gsc_placeholder_report = {'epsilon_zero_particles': n_zero,
                                        'replaced_sigma_counts_nm': changed}
        if changed:
            detail = ', '.join(f'{n} with sigma {v:g} nm' for v, n in sorted(changed.items()))
            self.ostream.print_info(f'GSC: {n_zero} particles have epsilon = 0; using sigma = '
                                    f'{target:g} nm in the GSC term instead of the placeholder '
                                    f'({detail}). LJ and electrostatics are unchanged.')
            self.ostream.flush()

    def _get_gaussian_softcore_force(self, lambda_val):
        """
        Define the Gaussian softcore potential.
        """
        gsc_energy = 'lambda * alpha * exp(-beta * (r/Rmin)^x); Rmin = 1.122462 * sigma; sigma = 0.5 * (sigma1 + sigma2);'
        gsc_force = mm.CustomNonbondedForce(gsc_energy)
        gsc_force.addGlobalParameter('lambda', lambda_val)
        gsc_force.addGlobalParameter('alpha', self.alpha)
        gsc_force.addGlobalParameter('beta', self.beta)
        gsc_force.addGlobalParameter('x', self.x)
        gsc_force.addPerParticleParameter('sigma')
        return gsc_force

    def _run_simulation(self, system, topology, positions, lambda_value, seed=None, vacuum=False):
        """
        Run the simulation using OpenMM.

        If seed = (positions, box_vectors, velocities) from the previous window
        of the same phase is given, the window starts from that state without
        energy minimization. Otherwise the input positions are minimized and
        velocities are drawn from the Maxwell-Boltzmann distribution.

        Returns simulated time in ns, real elapsed time in seconds, and the
        final (positions, box_vectors, velocities) for seeding the next window.
        """

        simulation = self._create_simulation(topology, system, vacuum=vacuum)
        integrator = simulation.integrator
        context = simulation.context
        temperature = self.temperature * unit.kelvin

        # The simulation is passed explicitly since it is deleted at the end
        def start_from_input(sim):
            sim.context.setPeriodicBoxVectors(*system.getDefaultPeriodicBoxVectors())
            sim.context.setPositions(positions)
            em_t0 = time.time()
            sim.minimizeEnergy(maxIterations=self.num_em_steps)
            self.ostream.print_info(f'Time spent in energy minimization: {time.time() - em_t0:.2f} s')
            sim.context.setVelocitiesToTemperature(temperature)

        seeded = False
        if seed is not None:
            seed_positions, seed_box, seed_velocities = seed
            if seed_box is not None:
                context.setPeriodicBoxVectors(*seed_box)
            context.setPositions(seed_positions)
            energy = context.getState(getEnergy=True).getPotentialEnergy()
            if np.isfinite(energy.value_in_unit(unit.kilojoule_per_mole)):
                seeded = True
                if seed_velocities is not None:
                    context.setVelocities(seed_velocities)
                else:
                    context.setVelocitiesToTemperature(temperature)
                self.ostream.print_info('Starting from the final state of the previous window.')
            else:
                self.ostream.print_info('Seed state has non-finite energy; minimizing input structure.')
        if not seeded:
            start_from_input(simulation)
        self.ostream.flush()

        equil_t0 = time.time()
        # Equilibration. A seeded start that fails falls back to minimization.
        try:
            simulation.step(self.num_equil_steps)
            equil_energy = context.getState(getEnergy=True).getPotentialEnergy()
            if not np.isfinite(equil_energy.value_in_unit(unit.kilojoule_per_mole)):
                raise ValueError('Non-finite energy after equilibration.')
        except (ValueError, mm.OpenMMException) as error:
            if not seeded:
                raise
            self.ostream.print_info(f'Seeded equilibration failed ({error}); '
                                    'restarting from the minimized input structure.')
            start_from_input(simulation)
            simulation.step(self.num_equil_steps)

        self.ostream.print_info(f'Time spent in pre-equilibration: {time.time() - equil_t0:.2f} s')
        self.ostream.flush()

        # Calculate production run time and interval
        interval = self.num_steps // self.num_snapshots
        simulated_ns = interval * self.num_snapshots * self.timestep / unit.nanoseconds

        output_folder = Path(self.output_folder)
        # XTC reporter is optional because it affects performance!
        if self.save_trajectory_xtc:
            trajectory_filename = output_folder / f'trajectory_{lambda_value}_stage{self.stage}.xtc'
            simulation.reporters.append(app.XTCReporter(str(trajectory_filename), interval))
            self.ostream.print_info(f"Trajectory will be saved in {str(trajectory_filename)}")
            self.ostream.flush()
        # Save system xml
        if self.save_system_xml:
            system_filename = output_folder / f'system_{lambda_value}_stage{self.stage}.xml'
            with open(system_filename, 'w') as f:
                f.write(mm.XmlSerializer.serialize(system))
            self.ostream.print_info(f"System will be saved in {str(system_filename)}")
            self.ostream.flush()
        # State data reporter
        if self.save_energies_txt:
            energy_filename = output_folder / f'energy_{lambda_value}_stage{self.stage}.txt'
            simulation.reporters.append(app.StateDataReporter(str(energy_filename), interval,
                                                              step=True, potentialEnergy=True, temperature=True))
            self.ostream.print_info(f"Energies will be saved in {str(energy_filename)}")
            self.ostream.flush()

        time_start = time.time()
        conformations = []
        energies = []
        periodic = system.usesPeriodicBoundaryConditions()
        kT = self.kT.value_in_unit(unit.kilojoule_per_mole)

        for _ in range(self.num_snapshots):
            simulation.step(interval)
            state = context.getState(getPositions=True, getEnergy=True)
            # Plain float32 coordinates in nm halve snapshot memory. All states
            # are re-evaluated at these same coordinates, so u_kn is consistent.
            frame = state.getPositions(asNumpy=True).value_in_unit(unit.nanometers).astype(np.float32)
            energy = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole) / kT  # Reduced potential energy
            box_vectors = (np.asarray(state.getPeriodicBoxVectors(asNumpy=True).value_in_unit(unit.nanometers))
                           if periodic else None)
            # Keep the sampled periodic cell with its positions so equilibration
            # removal and subsampling always select the complete configuration.
            conformations.append((frame, box_vectors))
            energies.append(energy)

        time_end = time.time()
        self._detect_equilibration(conformations, energies)

        final = context.getState(getPositions=True, getVelocities=True)
        final_state = (final.getPositions(asNumpy=True),
                       final.getPeriodicBoxVectors() if periodic else None,
                       final.getVelocities(asNumpy=True))

        elapsed_time = time_end - time_start
        ns_per_hour = (simulated_ns / elapsed_time) * 3600
        self.ostream.print_info(f'Lambda = {lambda_value} completed. Elapsed time: {elapsed_time:.2f} s')
        self.ostream.print_info(f'Performance: {ns_per_hour:.2f} ns/hour')
        self.ostream.flush()

        # Release the (GPU) context before the next window is created.
        del simulation, context, integrator

        return simulated_ns, elapsed_time, final_state

    @staticmethod
    def _free_force_group(system):
        """Highest force group not used by the input forces (for the GSC term)."""
        used = {f.getForceGroup() for f in system.getForces()}
        for group in range(31, -1, -1):
            if group not in used:
                return group
        return 0

    @staticmethod
    def _find_gsc_index(system):
        """Index of the driver's GSC force: the last matching CustomNonbondedForce."""
        indices = [i for i, f in enumerate(system.getForces())
                   if isinstance(f, mm.CustomNonbondedForce)
                   and f.getEnergyFunction().startswith('lambda * alpha * exp(')]
        return indices[-1] if indices else None

    @classmethod
    def _find_gsc_force(cls, system):
        index = cls._find_gsc_index(system)
        return None if index is None else system.getForce(index)

    def _gsc_only_lambda_dependence(self, forcefields, lambdas):
        """True if the systems differ only in the GSC lambda (stages 2 and 4).

        Each system is copied with its GSC lambda reset to 1 and serialized;
        identical XML proves that U_k(x) = U_rest(x) + lambda_k * E_gsc(x)
        exactly, because the GSC energy expression is linear in lambda.
        """
        if lambdas is None or len(lambdas) != len(forcefields) or len(forcefields) < 2:
            return False
        reference_xml, group = None, None
        for system in forcefields:
            clone = copy.deepcopy(system)
            index = self._find_gsc_index(clone)
            if index is None:
                return False
            gsc = clone.getForce(index)
            if group is None:
                group = gsc.getForceGroup()
            # The GSC force must be alone in its force group.
            others = [i for i, f in enumerate(clone.getForces())
                      if i != index and f.getForceGroup() == group]
            if gsc.getForceGroup() != group or others:
                return False
            for i in range(gsc.getNumGlobalParameters()):
                if gsc.getGlobalParameterName(i) == 'lambda':
                    gsc.setGlobalParameterDefaultValue(i, 1.0)
            xml = mm.XmlSerializer.serialize(clone)
            if reference_xml is None:
                reference_xml = xml
            elif xml != reference_xml:
                return False
        return True

    def _recalculate_energies_linear(self, forcefields, vacuum, lambdas):
        """Reconstruct u_kn from one full and one GSC-only evaluation per frame.

        Returns None (caller falls back to direct recalculation) if a spot
        check against direct evaluation of another state fails.
        """
        lambdas = np.asarray(lambdas, dtype=float)
        ref = int(np.argmax(np.abs(lambdas)))
        lam_ref = lambdas[ref]
        if lam_ref == 0.0:
            return None
        kT = self.kT.value_in_unit(unit.kilojoule_per_mole)
        group = self._find_gsc_force(forcefields[ref]).getForceGroup()
        context = self._create_energy_context(forcefields[ref], vacuum=vacuum)
        snapshots = self.snapshots
        u_kn = np.zeros((len(forcefields), len(snapshots)), dtype=np.float64)

        def set_frame(ctx, positions, box_vectors):
            if box_vectors is not None:
                ctx.setPeriodicBoxVectors(mm.Vec3(*box_vectors[0]), mm.Vec3(*box_vectors[1]),
                                          mm.Vec3(*box_vectors[2]))
            ctx.setPositions(positions.astype(np.float64))

        for n, (positions, box_vectors) in enumerate(snapshots):
            set_frame(context, positions, box_vectors)
            total = context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
            e_gsc = context.getState(getEnergy=True, groups={group}).getPotentialEnergy().value_in_unit(
                unit.kilojoule_per_mole)
            rest = total - e_gsc
            u_kn[:, n] = (rest + lambdas * (e_gsc / lam_ref)) / kT
        del context

        # Spot check against direct evaluation of the most distant state.
        check = int(np.argmax(np.abs(lambdas - lam_ref)))
        frames = sorted({0, len(snapshots) // 2, len(snapshots) - 1})
        context = self._create_energy_context(forcefields[check], vacuum=vacuum)
        for n in frames:
            set_frame(context, *snapshots[n])
            direct = context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(
                unit.kilojoule_per_mole) / kT
            if not abs(direct - u_kn[check, n]) <= 1e-2:
                self.ostream.print_info(f'Linear GSC reweighting check failed (|du| = '
                                        f'{abs(direct - u_kn[check, n]):.2e} kT); using direct recalculation.')
                return None
        del context
        return u_kn

    def _recalculate_energies_quadratic(self, forcefields, vacuum, lambdas):
        """Reconstruct u_kn from three states when U is quadratic in lambda.

        With linear coupling, every term in stages 1 and 3 is linear in lambda
        except a charged solute's periodic self-image energy, which is
        quadratic. Three evaluations per frame (first, middle and last state)
        therefore determine all states exactly. Returns None (caller falls back
        to direct recalculation) if a spot check against direct evaluation of
        other states fails.
        """
        lam = np.asarray(lambdas, dtype=float)
        k = len(lam)
        if lambdas is None or k < 4:
            return None
        middle = int(np.argmin(np.abs(lam - 0.5 * (lam[0] + lam[-1]))))
        basis = sorted({0, middle, k - 1})
        if len(basis) < 3:
            return None
        weights = np.vander(lam, 3) @ np.linalg.inv(np.vander(lam[basis], 3))  # k x 3
        kT = self.kT.value_in_unit(unit.kilojoule_per_mole)
        snapshots = self.snapshots

        def evaluate(state, frames):
            context = self._create_energy_context(forcefields[state], vacuum=vacuum)
            values = np.zeros(len(frames))
            for m, n in enumerate(frames):
                positions, box_vectors = snapshots[n]
                if box_vectors is not None:
                    context.setPeriodicBoxVectors(mm.Vec3(*box_vectors[0]), mm.Vec3(*box_vectors[1]),
                                                  mm.Vec3(*box_vectors[2]))
                context.setPositions(positions.astype(np.float64))
                values[m] = context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(
                    unit.kilojoule_per_mole) / kT
            del context
            return values

        all_frames = list(range(len(snapshots)))
        basis_u = np.array([evaluate(b, all_frames) for b in basis])   # 3 x N
        u_kn = weights @ basis_u

        # Spot check: the two non-basis states farthest from any basis state.
        others = [i for i in range(k) if i not in basis]
        others.sort(key=lambda i: -min(abs(lam[i] - lam[b]) for b in basis))
        frames = sorted({0, len(snapshots) // 2, len(snapshots) - 1})
        for state in others[:2]:
            direct = evaluate(state, frames)
            error = np.max(np.abs(direct - u_kn[state, frames]))
            if not error <= 1e-2:
                self.ostream.print_info(f'Quadratic reweighting check failed (|du| = {error:.2e} kT); '
                                        'using direct recalculation.')
                return None
        return u_kn

    def _recalculate_energies(self, forcefields, topology, vacuum=False, lambdas=None):
        # With linear coupling, stages 1 and 3 are (at most) quadratic in lambda.
        if self.stage in (1, 3) and self.linear_coupling:
            time_start = time.time()
            u_kn = self._recalculate_energies_quadratic(forcefields, vacuum, lambdas)
            if u_kn is not None:
                self.ostream.print_info(f'Used quadratic reweighting from 3 of {len(forcefields)} states '
                                        f'({time.time() - time_start:.2f} s).')
                self.ostream.flush()
                return u_kn

        # Stages 2 and 4 only scale the GSC term, so all states follow exactly
        # from one full and one GSC-only evaluation per frame.
        if self.stage in (2, 4) and self._gsc_only_lambda_dependence(forcefields, lambdas):
            time_start = time.time()
            u_kn = self._recalculate_energies_linear(forcefields, vacuum, lambdas)
            if u_kn is not None:
                self.ostream.print_info(f'Used linear GSC reweighting for {len(forcefields)} states '
                                        f'({time.time() - time_start:.2f} s).')
                self.ostream.flush()
                return u_kn

        snapshots = self.snapshots
        u_kn = np.zeros((len(forcefields), sum(self.N_k)), dtype=np.float64)
        kT = self.kT.value_in_unit(unit.kilojoule_per_mole)
        for k, forcefield in enumerate(forcefields):
            self.ostream.print_info(f"Recalculating energies for Forcefield {k}...")
            self.ostream.flush()
            time_start_forcefield = time.time()

            # A bare Context avoids the Simulation/reporter overhead.
            context = self._create_energy_context(forcefield, vacuum=vacuum)
            integrator = context.getIntegrator()

            for n, (positions, box_vectors) in enumerate(snapshots):
                if box_vectors is not None:
                    context.setPeriodicBoxVectors(mm.Vec3(*box_vectors[0]),
                                                  mm.Vec3(*box_vectors[1]),
                                                  mm.Vec3(*box_vectors[2]))
                context.setPositions(positions.astype(np.float64))
                state = context.getState(getEnergy=True)
                # All windows in this stage have the same temperature and
                # pressure. The common beta*p*V term can be omitted from MBAR,
                # but U must be evaluated in each snapshot's sampled cell.
                u_kn[k, n] = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole) / kT

            del context, integrator
            elapsed_time_trajectory = time.time() - time_start_forcefield
            self.ostream.print_info(f"Forcefield {k} energy recalculation took {elapsed_time_trajectory:.2f} seconds.")
            self.ostream.flush()

        return u_kn

    def _detect_equilibration(self, conformations, energies):
        """Discard burn-in only; preserve time order for downstream analysis."""
        a = np.asarray(energies, dtype=float)
        if a.ndim != 1 or len(a) != len(conformations) or len(a) < 4 or not np.all(np.isfinite(a)):
            raise ValueError('Need at least four finite energies matching the saved frames.')
        constant = np.ptp(a) == 0
        if constant:
            t0, g, neff = 0, None, None
        else:
            t0, g, neff = timeseries.detect_equilibration(a)
        t0 = int(t0)
        if len(a) - t0 < 4:
            raise ValueError('Too few production frames after equilibration detection; run longer.')
        self.snapshots.extend(conformations[t0:])
        self.N_k.append(len(a) - t0)
        self._equilibration_diagnostics.append({
            'discarded_frames': t0, 'retained_frames': len(a) - t0,
            'potential_energy_g': None if g is None else float(g),
            'potential_energy_neff': None if neff is None else float(neff),
            'constant_energy': constant,
        })
        self.ostream.print_info(f'Burn-in: {t0} frames; retained {len(a)-t0} '
                                'chronological frames for correlation analysis.')
        self.ostream.flush()

    def _validate_analysis_settings(self):
        if self.analysis_method not in ('block_bootstrap', 'subsample'):
            raise ValueError('analysis_method must be block_bootstrap or subsample.')
        if (not isinstance(self.bootstrap_samples, int) or self.bootstrap_samples < 50
                or not isinstance(self.min_bootstrap_blocks, int) or self.min_bootstrap_blocks < 2):
            raise ValueError('Require bootstrap_samples >= 50 and min_bootstrap_blocks >= 2.')
        factors = np.asarray(self.block_length_factors, dtype=float)
        if (factors.ndim != 1 or len(factors) < 2 or not np.all(np.isfinite(factors))
                or np.any(factors < 1) or np.any(np.diff(factors) <= 0)):
            raise ValueError('Use at least two increasing block_length_factors >= 1.')

    @staticmethod
    def _circular_block_indices(n, length, rng):
        """Circular moving blocks, trimmed to n; caller keeps windows separate."""
        starts = rng.integers(0, n, size=int(np.ceil(n / length)))
        return ((starts[:, None] + np.arange(length)) % n).ravel()[:n]

    def _calculate_free_energy(self, u_kn):
        """MBAR with temporal uncertainty and separate overlap diagnostics.

        Block lengths are heuristics based on observed g, not guarantees of
        independence. Inter-window and inter-stage trajectories are assumed
        independent. Burn-in uncertainty and unsampled states are not covered.
        All energy differences/uncertainties below are dimensionless (kT units).
        """
        self._validate_analysis_settings()
        u = np.asarray(u_kn, dtype=float)
        counts = np.asarray(self.N_k, dtype=int)
        if (u.ndim != 2 or u.shape[0] != len(counts) or u.shape[1] != counts.sum()
                or np.any(counts < 4) or not np.all(np.isfinite(u))):
            raise ValueError('Invalid u_kn/N_k or fewer than four production frames per window.')
        offsets = np.concatenate(([0], np.cumsum(counts)))
        windows, gs = [], []
        for k, n in enumerate(counts):
            sl = slice(offsets[k], offsets[k+1])
            signals = {'sampled_state_energy': u[k, sl]}
            for j in (k-1, k+1):
                if 0 <= j < len(counts):
                    signals[f'delta_u_to_{j}'] = u[j, sl] - u[k, sl]
            observables = {}
            valid_g = []
            for name, values in signals.items():
                if np.ptp(values) == 0:
                    observables[name] = {'g': None, 'temporal_neff': None, 'constant': True}
                    continue
                g = float(timeseries.statistical_inefficiency(values, fast=False))
                if not np.isfinite(g) or g < 1:
                    raise ValueError(f'Invalid correlation estimate in window {k}.')
                valid_g.append(g)
                observables[name] = {'g': g, 'temporal_neff': float(n/g), 'constant': False}
            g = max(valid_g, default=1.0)
            gs.append(g)
            windows.append({'window': k, 'frames': int(n), 'g_for_analysis': g,
                            'correlation_identifiable': bool(valid_g),
                            'observables': observables})
        analysis_counts = counts.copy()
        analysis_u = u
        if self.analysis_method == 'subsample':
            selections = []
            for k, n in enumerate(counts):
                # A supplied conservative g determines the spacing; use one
                # observable only as the length-bearing input to PyMBAR.
                local = timeseries.subsample_correlated_data(
                    u[k, offsets[k]:offsets[k+1]], g=gs[k], conservative=True)
                analysis_counts[k] = len(local)
                selections.extend(offsets[k] + np.asarray(local, dtype=int))
            if np.any(analysis_counts < 2):
                raise ValueError('Fewer than two decorrelated frames in a window; run longer.')
            analysis_u = u[:, selections]
        mbar = MBAR(analysis_u, N_k=analysis_counts, solver_protocol='robust')
        result = mbar.compute_free_energy_differences(compute_uncertainty=False)
        # Keep raw state counts for full-data MBAR, never replace N_k with N_eff.
        overlap = mbar.compute_overlap()
        weight_neff = mbar.compute_effective_sample_number()
        rng = np.random.default_rng(None if self.bootstrap_seed is None else
                                    np.random.SeedSequence([self.bootstrap_seed, self.stage]))
        records, estimates = [], []
        factors = self.block_length_factors if self.analysis_method == 'block_bootstrap' else (1,)
        current_offsets = np.concatenate(([0], np.cumsum(analysis_counts)))
        seen = set()
        for factor in factors:
            lengths = (np.maximum(1, np.ceil(np.asarray(gs)*factor).astype(int))
                       if self.analysis_method == 'block_bootstrap' else np.ones(len(counts), dtype=int))
            signature = tuple(int(x) for x in lengths)
            if signature in seen:
                continue
            seen.add(signature)
            block_counts = analysis_counts // lengths
            record = {'factor': float(factor), 'block_lengths_frames': list(signature),
                      'blocks_per_window': block_counts.tolist(),
                      'usable': bool(np.all(block_counts >= self.min_bootstrap_blocks))}
            records.append(record)
            if not record['usable']:
                record['reason'] = 'Too few blocks; longer trajectories needed.'
                continue
            draws = []
            for _ in range(self.bootstrap_samples):
                indices = np.concatenate([
                    current_offsets[k] + self._circular_block_indices(int(n), int(lengths[k]), rng)
                    for k, n in enumerate(analysis_counts)])
                boot = MBAR(analysis_u[:, indices], N_k=analysis_counts,
                            initial_f_k=mbar.f_k, solver_protocol='robust')
                # Fail rather than silently discard failed replicates.
                draws.append(boot.compute_free_energy_differences(
                    compute_uncertainty=False)['Delta_f'])
            sd = np.std(np.asarray(draws), axis=0, ddof=1)
            if not np.all(np.isfinite(sd)):
                raise ValueError('Non-finite bootstrap uncertainty.')
            record['endpoint_standard_error_kT'] = float(sd[-1, 0])
            estimates.append((record, sd))
        # Use the largest endpoint error among adequately populated block scales.
        if estimates:
            chosen, sd = max(estimates, key=lambda item: item[0]['endpoint_standard_error_kT'])
            result['dDelta_f'] = sd
            selected_lengths = chosen['block_lengths_frames']
        else:
            result['dDelta_f'] = np.full((len(counts), len(counts)), np.nan)
            selected_lengths = None
        stable = None
        if self.analysis_method == 'block_bootstrap' and len(estimates) >= 2:
            e1, e2 = [item[0]['endpoint_standard_error_kT'] for item in estimates[-2:]]
            stable = abs(e2-e1) <= 0.2 * max(e1, e2, np.finfo(float).eps)
        messages = []
        if not estimates:
            messages.append('Insufficient blocks: uncertainty unavailable (NaN).')
        if self.analysis_method == 'block_bootstrap' and stable is not True:
            messages.append('Block-length stability not established; uncertainty is provisional.')
        if any(not record['usable'] for record in records):
            messages.append('Some requested block lengths could not be assessed with enough blocks.')
        if any(not w['correlation_identifiable'] for w in windows):
            messages.append('Constant energy signals cannot establish decorrelation.')
        matrix = np.asarray(overlap['matrix'])
        adjacent = [float(min(matrix[k,k+1], matrix[k+1,k])) for k in range(len(counts)-1)]
        if any(x < 0.03 for x in adjacent):
            messages.append('Weak adjacent overlap (<0.03 heuristic); inspect lambda spacing.')
        diagnostics = {
            'method': self.analysis_method, 'windows': windows,
            'energy_units': 'kT', 'correlation_time_units': 'saved frames',
            'equilibration': list(self._equilibration_diagnostics),
            'analysis_counts': analysis_counts.tolist(),
            'overlap_matrix': matrix.tolist(), 'overlap_scalar': float(np.real(overlap['scalar'])),
            'adjacent_overlap_minima': adjacent,
            'weight_effective_sample_numbers': np.asarray(weight_neff).tolist(),
            'weight_neff_is_not_temporal_neff': True,
            'bootstrap_samples': self.bootstrap_samples,
            'bootstrap_seed': self.bootstrap_seed,
            'block_length_scan': records, 'selected_block_lengths': selected_lengths,
            'block_stability_20_percent_heuristic': stable, 'messages': messages,
            'independent_windows_and_stages_assumed': True,
        }
        result['sampling_diagnostics'] = diagnostics
        for w in windows:
            self.ostream.print_info(f"Window {w['window']}: {w['frames']} production frames, "
                                    f"g(max observed)={w['g_for_analysis']:.2f}")
        for record in records:
            self.ostream.print_info(f'Bootstrap block scan: {record}')
        self.ostream.print_info(f'Adjacent overlap minima: {adjacent}')
        for message in messages:
            self.ostream.print_info(message)
        self.ostream.flush()
        return result

    def plot_overlap_matrix(self, stage=None, filename=None, annotate=True):
        """
        Plot the MBAR overlap matrices of a finished calculation.

        The matrices are taken from the per-stage sampling diagnostics.
        Stages that were not simulated (e.g. Stage 4 with
        keep_intramolecular_gsc, Stages 3/4 with the decoupling scheme)
        are skipped.

        :param stage:
            The stage (1-4) to plot. If None, all simulated stages are plotted.
        :param filename:
            If given, the figure is saved to this file.
        :param annotate:
            If True, the matrix elements are written into the cells.

        :return:
            The matplotlib figure.
        """

        import matplotlib.pyplot as plt

        assert_msg_critical(self.delta_f is not None,
                            'SolvationFepDriver.plot_overlap_matrix: run a calculation first')

        if stage is None:
            stages = [1, 2, 3, 4]
        else:
            assert_msg_critical(stage in (1, 2, 3, 4),
                                'SolvationFepDriver.plot_overlap_matrix: stage must be 1, 2, 3 or 4')
            stages = [stage]

        panels = []
        for s in stages:
            diag = self.delta_f[s - 1].get('sampling_diagnostics', {})
            if 'overlap_matrix' not in diag:
                continue
            matrix = np.asarray(diag['overlap_matrix'], dtype=float)
            lambdas = getattr(self, f'lambdas_stage{s}', None)
            if lambdas is None or len(lambdas) != matrix.shape[0]:
                lambdas = list(range(matrix.shape[0]))
            adjacent = diag.get('adjacent_overlap_minima')
            if adjacent is None:
                adjacent = [min(matrix[k, k + 1], matrix[k + 1, k]) for k in range(matrix.shape[0] - 1)]
            panels.append((s, matrix, lambdas, adjacent))

        assert_msg_critical(len(panels) > 0,
                            'SolvationFepDriver.plot_overlap_matrix: no overlap data for the requested stage(s)')

        fig, axes = plt.subplots(1, len(panels), figsize=(4.0 * len(panels) + 1.0, 4.0), squeeze=False,
                                 constrained_layout=True)

        for ax, (s, matrix, lambdas, adjacent) in zip(axes[0], panels):
            n_states = matrix.shape[0]
            image = ax.imshow(matrix, cmap='Blues', vmin=0.0, vmax=1.0)
            if annotate:
                for i in range(n_states):
                    for j in range(n_states):
                        if matrix[i, j] >= 0.005:
                            color = 'white' if matrix[i, j] > 0.6 else 'black'
                            ax.text(j, i, f'{matrix[i, j]:.2f}'.lstrip('0'), ha='center', va='center',
                                    fontsize=7, color=color)
            labels = [f'{lam:g}' for lam in lambdas]
            ax.set_xticks(range(n_states))
            ax.set_yticks(range(n_states))
            ax.set_xticklabels(labels, rotation=90, fontsize=7)
            ax.set_yticklabels(labels, fontsize=7)
            ax.set_xlabel(r'$\lambda$')
            ax.set_ylabel(r'$\lambda$')
            min_adjacent = min(adjacent) if len(adjacent) > 0 else float('nan')
            ax.set_title(f'Stage {s} (min. adjacent {min_adjacent:.2f})', fontsize=9)

        fig.colorbar(image, ax=list(axes[0]), shrink=0.8, label='Overlap')

        if filename is not None:
            fig.savefig(filename, dpi=300, bbox_inches='tight')

        return fig

    def _get_alchemical_region(self, topology):
        """
        Define alchemical (perturbed) and chemical (unperturbed) regions.
        """
        alchemical_region = []
        chemical_region = []

        for res in topology.residues():
            condition = res.name == self.resname if self.resname else res.index == 0
            if condition:
                for atom in res.atoms():
                    alchemical_region.append(atom.index)
            else:
                for atom in res.atoms():
                    chemical_region.append(atom.index)

        if not alchemical_region:
            raise ValueError("No solute atoms match the selected residue.")
        return alchemical_region, chemical_region

    def _precision_properties(self, platform):
        if (self.precision is not None and platform.getName() in ('CUDA', 'OpenCL', 'HIP')
                and 'Precision' in platform.getPropertyNames()):
            return {'Precision': self.precision}
        return {}

    def _platform_attempts(self, vacuum=False):
        """Ordered (platform, properties) candidates for a new Context.

        An explicitly requested platform (self.platform) is used as given so
        that errors surface. Otherwise, like OpenMM's own default choice,
        platforms are tried from fastest to slowest; a registered platform is
        not necessarily usable for a given system (e.g. OpenCL on a Mac).
        The first working choice is remembered and tried first next time.
        """
        if self.platform is not None and not (vacuum and self.vacuum_platform):
            platform = mm.Platform.getPlatformByName(self.platform)
            properties = self._precision_properties(platform)
            properties.update(self.platform_properties or {})
            return [(platform, properties)], False
        names = []
        if vacuum and self.vacuum_platform:
            names.append(self.vacuum_platform)
        working = getattr(self, '_working_platforms', {}).get(vacuum)
        if working is not None:
            names.append(working)
        if self.platform is not None:
            names.append(self.platform)
        available = sorted((mm.Platform.getPlatform(i) for i in range(mm.Platform.getNumPlatforms())),
                           key=lambda p: p.getSpeed(), reverse=True)
        names.extend(p.getName() for p in available)
        attempts, seen = [], set()
        for name in names:
            if name in seen:
                continue
            seen.add(name)
            try:
                platform = mm.Platform.getPlatformByName(name)
            except Exception:
                continue
            base = self._precision_properties(platform)
            if self.platform_properties and name == self.platform:
                attempts.append((platform, dict(base, **self.platform_properties)))
            if base:
                attempts.append((platform, base))
            attempts.append((platform, {}))
        return attempts, True

    def _build_with_platform(self, build, vacuum=False):
        """Call build(platform, properties) with platform fallback."""
        attempts, fallback = self._platform_attempts(vacuum)
        if not fallback:
            return build(*attempts[0])
        errors = []
        for platform, properties in attempts:
            try:
                result = build(platform, properties)
            except Exception as error:
                errors.append(f'{platform.getName()} {properties}: {error}')
                continue
            if not hasattr(self, '_working_platforms'):
                self._working_platforms = {}
            if self._working_platforms.get(vacuum) != platform.getName():
                self._working_platforms[vacuum] = platform.getName()
                phase = 'vacuum' if vacuum else 'solvated'
                self.ostream.print_info(f'Using OpenMM platform {platform.getName()} '
                                        f'{properties} for the {phase} stages.')
                if errors:
                    self.ostream.print_info('Skipped platforms: ' + '; '.join(errors))
                self.ostream.flush()
            return result
        raise RuntimeError('No usable OpenMM platform found: ' + '; '.join(errors))

    def _create_simulation(self, topology, system, vacuum=False):
        # A fresh integrator per attempt: a failed Context may leave it bound.
        return self._build_with_platform(
            lambda platform, properties: app.Simulation(
                topology, system, self._create_integrator(), platform, properties),
            vacuum)

    def _create_energy_context(self, system, vacuum=False):
        return self._build_with_platform(
            lambda platform, properties: mm.Context(
                system, mm.VerletIntegrator(1.0 * unit.femtoseconds), platform, properties),
            vacuum)

    def _compute_kT(self):
        """Molar thermal energy at the current temperature setting."""
        return (unit.AVOGADRO_CONSTANT_NA * unit.BOLTZMANN_CONSTANT_kB *
                (float(self.temperature) * unit.kelvin)).in_units_of(unit.kilojoule_per_mole)

    def _validate_md_settings(self):
        if not np.isfinite(float(self.temperature)) or float(self.temperature) <= 0:
            raise ValueError('temperature must be a positive number in kelvin.')
        if self.ion_correction_box not in ('average', 'input'):
            raise ValueError("ion_correction_box must be 'average' or 'input'.")
        if self.integrator not in ('langevin_middle', 'langevin'):
            raise ValueError("integrator must be 'langevin_middle' or 'langevin'.")
        if not isinstance(self.keep_intramolecular_gsc, (bool, np.bool_)):
            raise ValueError('keep_intramolecular_gsc must be True or False.')
        if self.alchemical_scheme not in ('annihilation', 'decoupling'):
            raise ValueError("alchemical_scheme must be 'annihilation' or 'decoupling'.")
        if self.alchemical_scheme == 'decoupling' and not self.linear_coupling:
            raise ValueError("alchemical_scheme='decoupling' requires linear_coupling=True.")
        if (not isinstance(self.num_snapshots, numbers.Integral) or self.num_snapshots < 4
                or not isinstance(self.num_steps, numbers.Integral)
                or self.num_steps < self.num_snapshots):
            raise ValueError('Require num_snapshots >= 4 and num_steps >= num_snapshots.')
        if self.hydrogen_mass is not None:
            mass = self._hydrogen_mass_amu()
            if not np.isfinite(mass) or mass < 1.0 or mass > 5.0:
                raise ValueError('hydrogen_mass must be between 1 and 5 amu.')
            if self.constraints not in (app.HBonds, app.AllBonds, app.HAngles):
                raise ValueError('Hydrogen mass repartitioning requires constraints=HBonds or stronger.')
        dt = self.timestep.value_in_unit(unit.femtoseconds)
        if dt > 2.0 + 1e-8 and self.hydrogen_mass is None:
            self.ostream.print_info(f'Warning: timestep {dt} fs without hydrogen mass '
                                    'repartitioning may be unstable.')

    def _create_integrator(self):
        temperature = self.temperature * unit.kelvin
        if self.integrator == 'langevin_middle':
            return mm.LangevinMiddleIntegrator(temperature, self.friction, self.timestep)
        return mm.LangevinIntegrator(temperature, self.friction, self.timestep)

    def _hydrogen_mass_amu(self):
        mass = self.hydrogen_mass
        if hasattr(mass, 'value_in_unit'):
            return float(mass.value_in_unit(unit.dalton))
        return float(mass)

    def _repartition_hydrogen_masses(self, system, topology):
        """Move mass from heavy atoms to their bonded hydrogens (not water).

        Total mass and configurational free energies are unchanged; only
        the dynamics (and thus the stable timestep) is affected.
        """
        if self.hydrogen_mass is None:
            return
        target = self._hydrogen_mass_amu()
        hydrogen = app.element.hydrogen
        count = 0
        for atom1, atom2 in topology.bonds():
            if atom1.element is hydrogen:
                atom1, atom2 = atom2, atom1
            if atom2.element is not hydrogen or atom1.element in (None, hydrogen):
                continue
            symbols = sorted(a.element.symbol for a in atom2.residue.atoms() if a.element is not None)
            if symbols == ['H', 'H', 'O']:
                continue  # rigid water keeps its standard masses
            h_mass = system.getParticleMass(atom2.index).value_in_unit(unit.dalton)
            heavy_mass = system.getParticleMass(atom1.index).value_in_unit(unit.dalton)
            if h_mass == 0.0 or heavy_mass == 0.0:
                continue  # virtual sites / massless particles
            transfer = target - h_mass
            if heavy_mass - transfer <= 0.0:
                raise ValueError('Hydrogen mass repartitioning would give a non-positive heavy-atom mass.')
            system.setParticleMass(atom2.index, target * unit.dalton)
            system.setParticleMass(atom1.index, (heavy_mass - transfer) * unit.dalton)
            count += 1
        self.ostream.print_info(f'Hydrogen mass repartitioning: {count} hydrogens set to {target:.3f} amu.')

    def _record_average_box(self, frames, lambda_value):
        """Mean box of the retained frames of the coupled solvated window.

        L is the cube root of the mean volume. The box counts as cubic if
        every sampled box is axis-aligned with equal edges (an isotropic
        barostat keeps a cubic box cubic).
        """
        boxes = [np.asarray(box, dtype=float) for _, box in frames if box is not None]
        if not boxes:
            self._average_box = None
            return
        volumes = np.array([abs(np.linalg.det(b)) for b in boxes])
        lengths = volumes ** (1.0 / 3.0)
        cubic = all(np.allclose(b, np.diag(np.diag(b)), rtol=0, atol=1e-6)
                    and np.allclose(np.diag(b), b[0, 0], rtol=1e-6, atol=1e-8) for b in boxes)
        self._average_box = {
            'lambda': float(lambda_value),
            'frames': int(len(boxes)),
            'L_nm': float(np.mean(volumes) ** (1.0 / 3.0)),
            'L_std_nm': float(np.std(lengths)),
            'volume_nm3': float(np.mean(volumes)),
            'cubic': bool(cubic),
        }
        self.ostream.print_info(f'Coupled solvated window: mean box L = {self._average_box["L_nm"]:.4f} '
                                f'+/- {self._average_box["L_std_nm"]:.4f} nm over {len(boxes)} frames.')
        self.ostream.flush()

    def _report_barostat(self, system):
        """Record whether the solvated system samples NPT; align barostat temperature.

        A barostat already present in a supplied system (system XML or
        in-memory input) is set to the driver temperature, since Monte Carlo
        volume moves must use the same temperature as the thermostat. A
        different pressure is reported but not changed.
        """
        barostat_types = tuple(getattr(mm, name) for name in (
            'MonteCarloBarostat', 'MonteCarloAnisotropicBarostat',
            'MonteCarloMembraneBarostat', 'MonteCarloFlexibleBarostat') if hasattr(mm, name))
        barostats = [f for f in system.getForces() if isinstance(f, barostat_types)]
        self._has_barostat = bool(barostats)
        if not self._has_barostat and self.stage == 1:
            self.ostream.print_info('Warning: the solvated system has no barostat; '
                                    'the solvated stages sample NVT, not NPT.')
        target_T = float(self.temperature)
        for barostat in barostats:
            current_T = barostat.getDefaultTemperature().value_in_unit(unit.kelvin)
            if not np.isclose(current_T, target_T, rtol=0, atol=1e-6):
                barostat.setDefaultTemperature(target_T * unit.kelvin)
                if self.stage == 1:
                    self.ostream.print_info(f'Barostat temperature changed from {current_T:.2f} K to '
                                            f'the driver temperature {target_T:.2f} K.')
            if isinstance(barostat, mm.MonteCarloBarostat):
                current_p = barostat.getDefaultPressure().value_in_unit(unit.bar)
                target_p = self.pressure.value_in_unit(unit.bar)
                if not np.isclose(current_p, target_p, rtol=1e-6, atol=1e-6) and self.stage == 1:
                    self.ostream.print_info(f'Warning: the supplied barostat pressure ({current_p:.4f} bar) '
                                            f'differs from the driver pressure ({target_p:.4f} bar); '
                                            'the supplied value is kept.')
