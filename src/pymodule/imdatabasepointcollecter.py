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
import h5py
import copy
import math
from pathlib import Path
from sys import stdout
import sys
from time import time
import xml.etree.ElementTree as ET
from xml.dom import minidom
from contextlib import redirect_stderr
from io import StringIO

from .veloxchemlib import mpi_master
from .veloxchemlib import hartree_in_kcalpermol, bohr_in_angstrom
from .molecule import Molecule
from .outputstream import OutputStream
from .errorhandler import assert_msg_critical
from .solvationbuilder import SolvationBuilder
from .optimizationdriver import OptimizationDriver
from .scfrestdriver import ScfRestrictedDriver
from .scfunrestdriver import ScfUnrestrictedDriver
from .molecularbasis import MolecularBasis
from .scfgradientdriver import ScfGradientDriver
from .scfhessiandriver import ScfHessianDriver
from .xtbdriver import XtbDriver
from .xtbgradientdriver import XtbGradientDriver
from .xtbhessiandriver import XtbHessianDriver
from .interpolationdriver import InterpolationDriver
from .interpolationdatapoint import InterpolationDatapoint
from .imtrustradiusoptimizer import (IMTrustRadiusOptimizer, IMTrustRadiusPairTerms,
                                     IMTrustRadiusCachedOptimizer)


with redirect_stderr(StringIO()) as fg_err:
    import geometric

try:
    import openmm as mm
    import openmm.app as app
    import openmm.unit as unit
    from openmm.app.metadynamics import Metadynamics, BiasVariable
except ImportError:
    pass


class IMDatabasePointCollecter:
    """
    Collects new datapoints for an interpolated (MSI / IM) potential energy
    surface by running molecular dynamics on the current interpolation model.

    The QM region of the OpenMM system (system_from_molecule) moves on the
    interpolated PES of an InterpolationDriver. Along the trajectory a cheap
    proximity gate decides when the model is checked against QM (energy and
    gradient). If the error exceeds the thresholds, the structure is added
    as a new datapoint (optionally after a constrained re-optimisation) to the
    HDF5 database, the trajectory is reset to its start and continues;
    structures that pass the check are stored as references in an HDF5
    reference file. The run stops when the database holds
    desired_datapoint_density datapoints, when converged_cycle cycles of
    duration steps pass without a new datapoint, or after nsteps steps.

    The class is driven by IMForceFieldGenerator.compute, which sets the
    attributes listed below and calls system_from_molecule, update_settings
    and run_qmmm. Supported scope: one MPI rank and the ground state
    (roots_to_follow == [0]); run_qmmm checks this.

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.

    dynamics_settings keys (update_settings; required unless a default is
    given):
        - drivers: {'gs': (energy, gradient, Hessian driver), 'es': None};
          ScfRestrictedDriver/ScfUnrestrictedDriver with ScfGradientDriver and
          ScfHessianDriver, or XtbDriver with XtbGradientDriver and
          XtbHessianDriver.
        - basis_set_label: {'gs': basis label, 'es': ...}; the label is read
          also for xTB drivers.
        - roots_to_follow: must be [0].
        - ensemble: 'NVE', 'NVT' or 'NPT' (NPT needs a periodic, solvated
          system).
        - temperature: in K (NVT, NPT; the velocities start at 0.5 T).
        - friction: in 1/ps (NVT, NPT).
        - pressure: in atm, NPT only; default 1 atm.
        - timestep: in fs.
        - nsteps: maximum number of MD steps.
        - snapshots: number of trajectory frames, a frame every
          nsteps // snapshots steps; default: one frame per nsteps steps.
        - trajectory_file: PDB trajectory; default 'trajectory.pdb'
          (overwritten by every run).
        - load_system: path of an OpenMM XML state to start from (the
          collector writes 'checkpoint.xml' after every new datapoint);
          default None (start from the input structure).
        - reference_struc_energy_file: HDF5 reference file; default None,
          which means 'QM_ref_along_traj.h5' in the working directory. An
          existing file is preloaded, a missing one is created.
        - platform: OpenMM platform name ('Reference', 'CPU', 'CUDA',
          'OpenCL'); default None (OpenMM's choice).
        - openmm_precision: 'single', 'mixed' or 'double' (CUDA and OpenCL
          only); default None.
        - desired_datapoint_density: stop when the database holds this many
          datapoints.
        - converged_cycle: stop after this many convergence cycles without a
          new datapoint.
        - duration: number of MD steps of one convergence cycle.
        - energy_threshold: in kcal/mol per atom, |E_QM - E_IM| / natoms.
        - grad_rmsd_thrsh: in kcal/mol/angstrom, RMS over all 3N gradient
          components.
        - start_collect: number of steps after the start and after every
          reset before the proximity gate is evaluated; default 0.
        - sampling_drivers, sampling_settings: optional xTB pre-screen of the
          QM checks against a mirror database; only the sampling_settings
          keys 'enabled' and 'e_thrsh_kcal_per_atom' are read; default off.
        - metadynamics: optional well-tempered metadynamics, {'enabled',
          'variables', 'bias_factor' (> 1), 'hill_height_kjmol',
          'hill_frequency' (steps), 'save_frequency', 'bias_dir',
          'force_group' (default 11)}; every variable is {'type': 'distance',
          'angle' or 'torsion', 'atoms', 'min_nm'/'max_nm'/'width_nm' or
          'min_deg'/'max_deg'/'width_deg', 'periodic'}. Torsion atoms are
          1-based, distance and angle atoms 0-based. All hills are cleared at
          every reset of the trajectory.
        interpolation_settings and sampling_interpolation_settings are dicts
        keyed by root, {0: InterpolationDriver settings including
        'imforcefield_file'}.

    Attributes set by IMForceFieldGenerator before run_qmmm:
        - non_core_symmetry_groups: {'gs': symmetry information, 'es': ...};
          the symmetry information is the list [symmetry groups, rotor
          groups, regrouped atoms, core atoms, non-core atoms, rotatable bonds
          (0-based), index list, symmetry dihedral lists, [], [start, end] of
          the dihedrals in the internal-coordinate vector].
        - density_around_data_point: {root: number of datapoints in the
          database}.
        - use_opt_confidence_radius: [enabled, mode, initial confidence
          radius of new datapoints, energy weight e_x of the radius fit];
          only the mode 'multi_grad' is implemented; default
          [False, 'multi_grad', 0.5, 0.3].
        - gate_distance_prefilter: if True, the proximity gate of
          update_forces aligns only the reference structures whose
          interatomic-distance lower bound of the RMSD does not exclude them.
          The decisions and the candidate cache are those of the full scan.
          Default False.
        - trust_radius_fit_cache: if True, determine_trust_radius_gradient
          evaluates the interpolation once per reference and fit and fits
          the radii in closed form (IMTrustRadiusPairTerms). The groups are
          the same, and the loss and gradient equal those of
          IMTrustRadiusOptimizer to round-off. Default False.
        - identify_relevant_int_coordinates: (bool, user constraints); if
          True, new datapoints are re-optimised with their important internal
          coordinates frozen, together with the user constraints (1-based
          geomeTRIC index tuples or geomeTRIC strings); default (True, []).
        - consider_locality: if True, add further constraints when the
          re-optimised structure leaves the region it should describe. This
          may contaminate the database; use it only if the construction does
          not improve otherwise. Default False.
        - platform: OpenMM platform name.

    Other instance variables:
        - parent_ff, water_ff: OpenMM force fields of standard residues and
          water (solvated systems).
        - padding: solvation box padding in nm; cutoff: nonbonded cutoff in nm.
        - integrator: an OpenMM integrator to use instead of the one created
          for the ensemble; default None.
        - system, topology, simulation, pdb: the OpenMM objects.
        - molecule, labels: the input molecule and its atom labels.
        - qm_atoms, mm_subregion, linking_atoms, qm_force_index: the QM/MM
          partition and the index of the QM force in the system.
        - linking_atom_distance: distance of a linking H atom in angstrom.
        - results of run_qmmm: qm_potentials, qm_mm_interaction_energies,
          mm_potentials, total_potentials, kinetic_energies, total_energies
          (kJ/mol, per step), temperatures (K), state_specific_molecules
          ({root: structure of every step}) and allowed_molecules ({root:
          reference structures, QM and interpolated energies, QM gradients}).

    Reference file (HDF5): one entry per QM check that passed the thresholds
    in /states/<root>/: qm_energy (Hartree), im_energy (Hartree), natoms, xyz
    (xyz string, angstrom), qm_grad_flat (Hartree/bohr, natoms*3 values),
    coords_flat (angstrom, natoms*3 values). The file attributes natoms,
    labels, charge and multiplicity identify the QM region; a file whose
    labels, charge or multiplicity differ is refused when it is preloaded.

    Files written to the working directory: <filename>.xml/.pdb and
    <filename>_system.xml/.pdb from system_from_molecule (filename
    'qm_region' with a QM region; a solvated system also writes the
    SolvationBuilder files system.pdb and solvent_1.xml), the trajectory
    file, checkpoint.xml, pre_opt_database_strucs.xyz (appended) and, by
    default, QM_ref_along_traj.h5.
    """

    # reference file used when reference_struc_energies_file is None
    _DEFAULT_REFERENCE_FILE = 'QM_ref_along_traj.h5'

    @property
    def identfy_relevant_int_coordinates(self):
        """
        Deprecated alias of identify_relevant_int_coordinates.
        """

        return self.identify_relevant_int_coordinates

    @identfy_relevant_int_coordinates.setter
    def identfy_relevant_int_coordinates(self, value):
        self.identify_relevant_int_coordinates = value

    def __init__(self, comm=None, ostream=None):
        """
        Initializes the class with default simulation parameters.
        """
        # MPI and output stream
        if comm is None:
            comm = MPI.COMM_WORLD

        if ostream is None:
            if comm.Get_rank() == mpi_master():
                ostream = OutputStream(stdout)
            else:
                ostream = OutputStream(None)

        self.comm = comm
        self.rank = self.comm.Get_rank()
        self.nodes = self.comm.Get_size()

        # output stream
        self.ostream = ostream

        # Simulation parameters
        self.platform = None
        self.openmm_precision = None
        self.ensemble = None
        self.temperature = None
        self.friction = None
        self.timestep = None
        self.nsteps = None
        self.parent_ff = 'amber03.xml'
        self.water_ff = 'spce.xml'
        self.padding = 1.0
        self.cutoff = 1.0
        self.integrator = None

        # meta dynamics settings
        self.metadynamics_settings = None
        self.metadynamics_enabled = None
        self._metadynamics = None
        self._metadynamics_forcegroup = 11
        # index of the bias force that _initialize_metadynamics added to the system
        self._metadynamics_bias_force_index = None

        # OpenMM objects
        self.system = None
        self.topology = None
        self.simulation = None
        self.pdb = None
        self.labels = []

        # VeloxChem objects
        self.molecule = None

        # QM Region parameters
        self.drivers = None
        self.qm_atoms = None
        self.mm_subregion = None
        self.linking_atoms = None
        self.qm_force_index = None

        self.prev_dens_of_points = None

        self.roots_to_follow = []

        self.snapshots = None
        self.load_system = None
        self.pressure = None
        self.out_file = 'trajectory.pdb'

        # database construction settings (set through update_settings)
        self.duration = None
        self.desired_datpoint_density = None
        self.unadded_cycles = None
        self.root_z_matrix = None

        self.density_around_data_point = None
        self.impes_drivers = None
        self.energy_threshold = None
        self.gradient_rmsd_thrsh = None
        self.start_collect = 0
        # state-resolved symmetry metadata excluding the invariant core atoms
        self.non_core_symmetry_groups = []
        self.interpolation_settings = None
        self.allowed_molecules = None
        self.reference_struc_energies_file = None

        self.current_state = None
        self.starting_state = 0
        self.current_gradient = 0
        self.current_energy = 0
        self.point_checker = 1
        self.last_added = 0
        self.last_point_added = None
        self.add_a_point = False
        self.skipping_value = 0
        self.basis_set_label = None
        self.consider_locality = False

        self.sampling_enabled = False
        self.sampling_settings = None
        self.sampling_driver = None
        self.sampling_impes_drivers = None
        self.sampling_interpolation_settings = None

        self.start_velocities = None
        self.coordinates = None
        self.state_specific_molecules = None
        self.inv_sqrt_masses = None
        self.use_mass_weight = False

        # toggle + optional fixed list for important internal-coordinate constraints
        self.identify_relevant_int_coordinates = (True, [])

        # [enabled, mode, initial confidence radius, energy weight e_x];
        # same default as IMForceFieldGenerator.use_opt_confidence_radius
        self.use_opt_confidence_radius = [False, 'multi_grad', 0.5, 0.3]
        self.confidence_radius_optimized = True

        # exact speed-ups of the proximity gate and the radius fit (opt-in,
        # same results; see the class docstring)
        self.gate_distance_prefilter = False
        self.trust_radius_fit_cache = False

        # Default value for the C-H linker distance
        self.linking_atom_distance = 1.0705

        # Cache of the most recent QM molecule passed through interpolation.
        self._latest_qm_molecule = None
        self._latest_qm_positions_nm = None

        # caches used by update_forces; run_qmmm refreshes the system-dependent ones
        self._topology_atom_labels = None
        self._gradient_to_force_factor = None
        self._qm_external_force = None
        self.previous_candidate_indices = []
        # interatomic distances of the reference structures (gate_distance_prefilter)
        self._gate_bound_cache = None
        # raw weights at radius 0.5 of every (reference, datapoint) pair (trust_radius_fit_cache)
        self._trust_radius_grouping_cache = None

    def _create_platform(self):
        """
        Creates the OpenMM platform named by self.platform.

        :return:
            The OpenMM platform, or None to let OpenMM choose.
        """

        if self.platform is None:
            return None

        platform_name = str(self.platform).strip()
        platform = mm.Platform.getPlatformByName(platform_name)

        # Set single thread because for single molecules multithreading
        # often has little benefit and can add overhead.
        if platform_name.upper() == "CPU":
            platform.setPropertyDefaultValue("Threads", "1")

        # Optional high-precision coordinates on GPU backends.
        if self.openmm_precision is not None and platform_name.upper() in ("CUDA", "OPENCL"):
            precision_map = {
                "single": "single",
                "mixed": "mixed",
                "double": "double",
            }
            precision_key = str(self.openmm_precision).strip().lower()
            if precision_key not in precision_map:
                raise ValueError(
                    "Unsupported OpenMM precision setting. "
                    "Use one of: single, mixed, double.")
            platform.setPropertyDefaultValue("Precision", precision_map[precision_key])

        return platform

    # Method to generate OpenMM system from VeloxChem objects
    def system_from_molecule(self,
                             molecule,
                             z_matrix,
                             ff_gen,
                             solvent='gas',
                             qm_atoms=None,
                             filename='residue',
                             trust_radius=False,
                             residue_name='MOL'):
        """
        Builds the OpenMM system of a molecule whose QM region is driven by the
        interpolation, and writes the OpenMM files of the system to the working
        directory (<filename>.xml/.pdb and <filename>_system.xml/.pdb, with
        filename 'qm_region' when a QM region is given).

        :param molecule:
            VeloxChem molecule object.
        :param z_matrix:
            The z-matrix (internal coordinates) of every root, {root: {'bonds',
            'angles', 'dihedrals', 'impropers'}}; stored as root_z_matrix.
        :param ff_gen:
            VeloxChem forcefield generator object (MMForceFieldGenerator).
        :param solvent:
            Available options:'gas', 'spce', 'tip3p', 'ethanol', 'methanol', 'acetone',
            'chloroform', 'hexane', 'toluene', 'dcm', 'benzene', 'dmso', 'thf',
            'acetonitrile', 'other' or 'itself'.
        :param qm_atoms:
            Options: None, 'all', or list of atom indices for QM region.
        :param filename:
            Base name for the generated files without a QM region. Default is 'residue'.
        :param trust_radius:
            If True, 'trust_radius' is appended to the file names of the QM region.
        :param residue_name:
            Name of the residue. Default is 'MOL'.
        """

        assert_msg_critical('openmm' in sys.modules, 'OpenMM is required for IMDatabasePointCollecter.')

        # Store the molecule object and generate OpenMM compatible files
        self.molecule = molecule
        self.positions = molecule.get_coordinates_in_angstrom()
        self.labels = molecule.get_labels()
        self.root_z_matrix = z_matrix

        # Options for the QM region if it's required.
        # TODO: Take this if else tree to a separate method.
        if qm_atoms:
            filename = 'qm_region'
            if trust_radius:
                filename += 'trust_radius'
            # Create a QM region topology template
            if qm_atoms == 'all':
                msg = 'Full molecule as QM region'
                self.ostream.print_info(msg)
                self.ostream.flush()
                qm_atoms = list(range(len(self.labels)))
                self._create_QM_residue(ff_gen,
                                        qm_atoms,
                                        filename)

            elif isinstance(qm_atoms, list):
                if qm_atoms == list(range(len(self.labels))):
                    msg = 'Full molecule as QM region'
                    self.ostream.print_info(msg)
                    self.ostream.flush()
                    self._create_QM_residue(ff_gen,
                                            qm_atoms,
                                            filename)
                else:
                    msg = 'QM/MM partition inside the molecule'
                    self.ostream.print_info(msg)
                    self.ostream.flush()
                    self._create_QM_subregion(ff_gen,
                                              qm_atoms,
                                              molecule,
                                              filename,
                                              )

            ff_gen.write_pdb(f'{filename}.pdb', 'QMR')

        elif qm_atoms is None:
            ff_gen.write_openmm_files(filename, residue_name)

        # TODO: Incorporate the SystemBuilder here to avoid using the Modeller.
        if solvent == 'gas':
            phase = 'gas'
            self.pdb = app.PDBFile(f'{filename}.pdb')
            # Common forcefield loading, modified according to phase specifics
            forcefield_files = [f'{filename}.xml']

        if solvent != 'gas':
            # Solvate the molecule using the SystemBuilder
            phase = 'periodic'
            sys_builder = SolvationBuilder()
            sys_builder.solvate(solute=molecule,
                                solvent=solvent,
                                padding=self.padding,
                                equilibrate=True,
                                constraint_solute=True
                                )
            sys_builder.write_openmm_files(solute_ff=ff_gen)

            if solvent == 'spce':
                self.water_ff = 'spce.xml'
                forcefield_files = [f'{filename}.xml', self.parent_ff, self.water_ff]
            elif solvent == 'tip3p':
                self.water_ff = 'tip3p.xml'
                forcefield_files = [f'{filename}.xml', self.parent_ff, self.water_ff]
            elif solvent == 'itself':
                # Solute and solvent are the same
                forcefield_files = [f'{filename}.xml', self.parent_ff]
            else:
                # Any other solvent, including custom ones
                solvent_ff = 'solvent_1.xml'
                forcefield_files = [f'{filename}.xml', self.parent_ff, solvent_ff]

            # Load the PDB from the SystemBuilder
            self.pdb = app.PDBFile('system.pdb')

        # Create the ForceField object
        forcefield = app.ForceField(*forcefield_files)

        # Create the System object
        if phase == 'gas':
            self.system = forcefield.createSystem(self.pdb.topology, nonbondedMethod=app.NoCutoff, constraints=app.HBonds)

        else:
            self.system = forcefield.createSystem(self.pdb.topology,
                                                  nonbondedMethod=app.PME,
                                                  nonbondedCutoff=self.cutoff * unit.nanometer,
                                                  constraints=app.HBonds)

        # Modify the system to include QM and QM/MM forces.
        if qm_atoms:
            self.set_qm_mm_system(phase,
                                  ff_gen)

        # Write the system to a xml file (for debugging purposes)
        with open(f'{filename}_system.xml', 'w') as f:
            f.write(mm.XmlSerializer.serialize(self.system))
            msg = f'System parameters written to {filename}_system.xml'
            self.ostream.print_info(msg)
            self.ostream.flush()

        # Write the system to a pdb file (for debugging purposes)
        if phase == 'gas':
            with open(f'{filename}_system.pdb', 'w') as pdb_file:
                app.PDBFile.writeFile(self.pdb.topology,
                                      self.positions,
                                      pdb_file)
            msg = f'System coordinates written to {filename}_system.pdb'
            self.ostream.print_info(msg)
            self.ostream.flush()
        elif phase == 'periodic':
            with open(f'{filename}_system.pdb', 'w') as pdb_file:
                app.PDBFile.writeFile(self.pdb.topology,
                                      self.pdb.positions,
                                      pdb_file)
            msg = f'System coordinates written to {filename}_system.pdb'
            self.ostream.print_info(msg)
            self.ostream.flush()

    def _refresh_interpolation_driver_caches(self, root):
        """
        Marks the runtime caches of the interpolation driver of one root (and of
        its sampling mirror driver) as stale. Call it after every change of the
        datapoints or of their confidence radii.

        :param root:
            The root.
        """
        driver = self.impes_drivers[root]
        driver.mark_runtime_data_cache_dirty()
        if self.sampling_enabled:
            sampling_driver = self.sampling_impes_drivers[root]
            sampling_driver.mark_runtime_data_cache_dirty()

    def _build_metadynamics_cv_force(self, spec):
        """
        Builds the OpenMM force that defines one metadynamics collective variable.

        :param spec:
            One entry of metadynamics['variables']: 'type' ('distance', 'angle' or
            'torsion'), 'atoms' (torsion atoms 1-based, distance and angle atoms
            0-based), 'min_nm'/'max_nm'/'width_nm' (distance) or
            'min_deg'/'max_deg'/'width_deg' (angle, torsion) and 'periodic'.

        :return:
            The force, the minimum, maximum and width of the variable (nm or
            radian) and the periodic flag.
        """

        vtype = spec['type'].lower()
        atoms = spec['atoms']

        if vtype == 'distance':
            cv_force = mm.CustomBondForce("r")
            cv_force.addBond(int(atoms[0]), int(atoms[1]), [])
            min_val = float(spec['min_nm'])
            max_val = float(spec['max_nm'])
            width = float(spec['width_nm'])
            periodic = bool(spec.get('periodic', False))

        elif vtype == 'angle':
            cv_force = mm.CustomAngleForce("theta")
            cv_force.addAngle(int(atoms[0]), int(atoms[1]), int(atoms[2]), [])
            min_val = np.deg2rad(float(spec['min_deg']))
            max_val = np.deg2rad(float(spec['max_deg']))
            width = np.deg2rad(float(spec['width_deg']))
            periodic = bool(spec.get('periodic', False))

        elif vtype == 'torsion':

            cv_force = mm.CustomTorsionForce("theta")
            cv_force.addTorsion(int(atoms[0]) - 1, int(atoms[1]) - 1, int(atoms[2]) - 1, int(atoms[3]) - 1)
            min_val = np.deg2rad(float(spec['min_deg']))
            max_val = np.deg2rad(float(spec['max_deg']))
            width = np.deg2rad(float(spec['width_deg']))
            periodic = bool(spec.get('periodic', False))

        cv_force.setForceGroup(int(self._metadynamics_forcegroup))
        return cv_force, min_val, max_val, width, periodic

    def _initialize_metadynamics(self):
        """
        Adds the well-tempered metadynamics bias (openmm.app.Metadynamics) for
        the variables of metadynamics_settings to the system; does nothing if
        metadynamics is disabled.
        """

        if not self.metadynamics_enabled:
            self._metadynamics = None
            return

        cfg = self.metadynamics_settings

        self._metadynamics_forcegroup = int(cfg.get('force_group', 11))

        variables = []

        for spec in cfg['variables']:

            cv_force, min_val, max_val, width, periodic = self._build_metadynamics_cv_force(spec)

            variables.append(BiasVariable(cv_force, min_val, max_val, width, periodic))

        self._metadynamics = Metadynamics(
            self.system,
            variables,
            self.temperature,
            float(cfg['bias_factor']),
            float(cfg['hill_height_kjmol']) * unit.kilojoules_per_mole,
            int(cfg['hill_frequency']),
            saveFrequency=(int(cfg['save_frequency']) if cfg.get('save_frequency') else None),
            biasDir=cfg.get('bias_dir', None),
        )
        # Usually the newly added metadynamics bias is the last force
        bias_force = self.system.getForce(self.system.getNumForces() - 1)
        bias_force.setForceGroup(self._metadynamics_forcegroup)
        self._metadynamics_bias_force_index = self.system.getNumForces() - 1

    def _remove_metadynamics_bias_force(self):
        """
        Removes the bias force that an earlier _initialize_metadynamics added to
        the system, so that a repeated run_qmmm does not keep a frozen bias.
        """

        index = self._metadynamics_bias_force_index
        self._metadynamics_bias_force_index = None
        if index is None or self.system is None:
            return
        if index < self.system.getNumForces() and isinstance(
                self.system.getForce(index), mm.CustomCVForce):
            self.system.removeForce(index)

    def _reset_metadynamics_bias(self, context):
        """
        Soft reset only: clears the accumulated metadynamics hills and restarts the
        bias growth from zero (called at every reset of the trajectory).

        :param context:
            The OpenMM context.

        :return:
            True if a bias was reset.
        """
        if self._metadynamics is None or not self.metadynamics_enabled:
            return False

        mtd = self._metadynamics
        mtd._selfBias.fill(0.0)
        mtd._totalBias.fill(0.0)
        mtd._loadedBiases = {}

        if len(mtd.variables) == 1:
            mtd._table.setFunctionParameters(mtd._totalBias.flatten(), *mtd._limits)
        else:
            mtd._table.setFunctionParameters(
                *mtd._widths,
                mtd._totalBias.flatten(),
                *mtd._limits
            )
        mtd._force.updateParametersInContext(context)

        return True

    def _validate_metadynamics_settings(self):
        """
        Checks metadynamics_settings: required keys, variable types, atom index
        ranges (torsion atoms 1-based in [1, natoms], distance and angle atoms
        0-based in [0, natoms - 1]), range keys and positive bias parameters.
        Raises ValueError for invalid settings.
        """

        cfg = self.metadynamics_settings

        if Metadynamics is None or BiasVariable is None:

            raise RuntimeError('metadynamics requested, but it is currently unavailable!')

        required = ('variables', 'bias_factor', 'hill_height_kjmol', 'hill_frequency')

        missing = [k for k in required if k not in cfg]

        if missing:
            raise ValueError(f"metadynamics is missing the following variables {missing}")

        variables = cfg['variables']
        if not isinstance(variables, list) or len(variables) == 0:
            raise ValueError("metadynamics.variables must be a non-empty list.")

        natoms = self.system.getNumParticles() if self.system is not None else None
        for idx, var in enumerate(variables):
            vtype = var.get('type', '').lower()
            atoms = var.get('atoms', [])
            if vtype not in ('distance', 'angle', 'torsion'):
                raise ValueError(f"metadynamics.variables[{idx}].type must be distance|angle|torsion.")
            expected = {'distance': 2, 'angle': 3, 'torsion': 4}[vtype]
            if len(atoms) != expected:
                raise ValueError(f"metadynamics.variables[{idx}].atoms must have {expected} entries.")
            if natoms is not None:
                # _build_metadynamics_cv_force reads torsion atoms as 1-based
                # and distance/angle atoms as 0-based indices
                first = 1 if vtype == 'torsion' else 0
                last = first + natoms - 1
                for a in atoms:
                    if a < first or a > last:
                        raise ValueError(
                            f"metadynamics.variables[{idx}] atom index {a} "
                            f"out of range [{first}, {last}].")

            # Unit-range keys
            if vtype in ('angle', 'torsion'):
                for k in ('min_deg', 'max_deg', 'width_deg'):
                    if k not in var:
                        raise ValueError(f"metadynamics.variables[{idx}] missing key '{k}'.")
            else:
                for k in ('min_nm', 'max_nm', 'width_nm'):
                    if k not in var:
                        raise ValueError(f"metadynamics.variables[{idx}] missing key '{k}'.")

        if float(cfg['bias_factor']) <= 1.0:
            raise ValueError("metadynamics.bias_factor must be > 1.0 for well-tempered metadynamics.")
        if float(cfg['hill_height_kjmol']) <= 0.0:
            raise ValueError("metadynamics.hill_height_kjmol must be > 0.")
        if int(cfg['hill_frequency']) <= 0:
            raise ValueError("metadynamics.hill_frequency must be > 0.")

    def update_settings(self, dynamics_settings, interpolation_settings=None, sampling_interpolation_settings=None):
        """
        Reads the dynamics and database-construction settings (the keys, units and
        defaults are listed in the class docstring).

        :param dynamics_settings:
            The dictionary of dynamics settings.
        :param interpolation_settings:
            The interpolation settings keyed by root, {root: settings}; each entry
            holds the InterpolationDriver settings and 'imforcefield_file'.
        :param sampling_interpolation_settings:
            The interpolation settings of the sampling mirror database keyed by
            root (only used when sampling is enabled).
        """

        self.drivers = dynamics_settings['drivers']

        self.interpolation_settings = interpolation_settings
        self.sampling_interpolation_settings = sampling_interpolation_settings

        self.use_mass_weight = interpolation_settings[list(interpolation_settings.keys())[0]].get('use_mass_weight', False) if interpolation_settings else False

        if 'duration' in dynamics_settings:
            self.duration = float(dynamics_settings['duration'])

        if 'temperature' in dynamics_settings:
            self.temperature = dynamics_settings['temperature']

        if 'pressure' in dynamics_settings:
            self.pressure = float(dynamics_settings['pressure'])

        if 'friction' in dynamics_settings:
            self.friction = float(dynamics_settings['friction'])

        if 'timestep' in dynamics_settings:
            self.timestep = float(dynamics_settings['timestep'])

        if 'nsteps' in dynamics_settings:
            self.nsteps = int(dynamics_settings['nsteps'])

        if 'snapshots' in dynamics_settings:
            self.snapshots = int(dynamics_settings['snapshots'])

        # Start the simulation from a given OpenMM XML state (e.g. checkpoint.xml)
        if 'load_system' in dynamics_settings:
            self.load_system = dynamics_settings['load_system']

        if 'trajectory_file' in dynamics_settings:
            self.out_file = dynamics_settings['trajectory_file']
        else:
            self.out_file = 'trajectory.pdb'

        if 'reference_struc_energy_file' in dynamics_settings:
            self.reference_struc_energies_file = dynamics_settings['reference_struc_energy_file']

        if 'platform' in dynamics_settings:
            self.platform = dynamics_settings['platform']

        if 'openmm_precision' in dynamics_settings:
            precision_value = str(dynamics_settings['openmm_precision']).strip().lower()
            if precision_value not in ('single', 'mixed', 'double'):
                raise ValueError(
                    "dynamics_settings['openmm_precision'] must be one of: "
                    "single, mixed, double.")
            self.openmm_precision = precision_value

        # Determines the ensemble in order to set the correct simulation set_up
        if 'ensemble' in dynamics_settings:
            self.ensemble = dynamics_settings['ensemble']

        if "metadynamics" in dynamics_settings:
            self.metadynamics_settings = copy.deepcopy(dynamics_settings["metadynamics"])

        self.metadynamics_enabled = bool(self.metadynamics_settings is not None and
                                         self.metadynamics_settings.get('enabled', False))

        if self.metadynamics_enabled:
            self._validate_metadynamics_settings()
        # DATABASE construction inputs

        # stop when the database holds this many datapoints
        if 'desired_datapoint_density' in dynamics_settings:
            self.desired_datpoint_density = int(dynamics_settings['desired_datapoint_density'])

        # stop after this many cycles of 'duration' steps without a new datapoint
        if 'converged_cycle' in dynamics_settings:
            self.unadded_cycles = int(dynamics_settings['converged_cycle'])

        if 'basis_set_label' in dynamics_settings:
            basis_set_label = dynamics_settings['basis_set_label']
            self.basis_set_label = basis_set_label

        # number of steps after the start and after every reset before the
        # proximity gate (and with it the QM checks) is evaluated
        if 'start_collect' in dynamics_settings:
            self.start_collect = dynamics_settings['start_collect']

        # roots to follow; only [0] (ground state) is supported
        if "roots_to_follow" in dynamics_settings:
            self.roots_to_follow = list(dynamics_settings['roots_to_follow'])

        if 'energy_threshold' in dynamics_settings:
            self.energy_threshold = dynamics_settings['energy_threshold']

        if 'grad_rmsd_thrsh' in dynamics_settings:
            self.gradient_rmsd_thrsh = dynamics_settings['grad_rmsd_thrsh']

        if 'sampling_drivers' in dynamics_settings:
            sampling_settings = dynamics_settings.get('sampling_settings', {
                'enabled':False,
                'e_thrsh_kcal_per_atom': 0.1,
                'g_thrsh_kcal_ang_per_atom':2.0,
                'force_orient_cos': 0.0001})
            self.sampling_enabled = bool(sampling_settings.get('enabled'))
            self.sampling_driver = dynamics_settings.get('sampling_drivers', None)
            self.sampling_settings = sampling_settings

    def _build_run_qmmm_runtime_cache(self):
        """
        Builds cached constants used during the QM/MM integration loop.
        """

        energy_unit = unit.kilojoules_per_mole
        gas_constant = unit.MOLAR_GAS_CONSTANT_R.value_in_unit(energy_unit / unit.kelvin)

        return {
            'energy_unit': energy_unit,
            'qm_mm_groups': {1, 2},
            'mm_groups': {3, 4, 5, 6, 7},
            'dof': self._count_degrees_of_freedom(),
            'gas_constant': gas_constant,
            'metadynamics_enabled': bool(self._metadynamics is not None),
            'metadynamics_group':{int(self._metadynamics_forcegroup)} if self._metadynamics is not None else None,
        }

    def _count_degrees_of_freedom(self):
        """
        Counts the degrees of freedom for the reported temperature the way
        OpenMM's StateDataReporter does: 3 per massive particle, minus one per
        constraint, minus 3 if a CMMotionRemover removes the centre-of-mass motion.

        :return:
            The number of degrees of freedom.
        """

        system = self.system
        zero_mass = 0.0 * unit.dalton
        dof = 0
        for i in range(system.getNumParticles()):
            if system.getParticleMass(i) > zero_mass:
                dof += 3
        for i in range(system.getNumConstraints()):
            p1, p2, _ = system.getConstraintParameters(i)
            if system.getParticleMass(p1) > zero_mass or system.getParticleMass(p2) > zero_mass:
                dof -= 1
        if any(isinstance(force, mm.CMMotionRemover) for force in system.getForces()):
            dof -= 3
        return dof

    def _collect_step_observables_optimized(self, context, runtime_cache):
        """
        Collects the energies and the temperature of the current step.

        :param context:
            The OpenMM context.
        :param runtime_cache:
            The constants from _build_run_qmmm_runtime_cache.

        :return:
            A dictionary with the energies in kJ/mol ('qm', 'qm_mm', 'mm',
            'metad_bias', 'potential', 'kinetic', 'total') and 'temperature' in K.
        """

        energy_unit = runtime_cache['energy_unit']

        qm = float(self.current_energy)
        qm_mm = context.getState(
            getEnergy=True,
            groups=runtime_cache['qm_mm_groups']).getPotentialEnergy().value_in_unit(energy_unit)
        mm = context.getState(
            getEnergy=True,
            groups=runtime_cache['mm_groups']).getPotentialEnergy().value_in_unit(energy_unit)
        kinetic = context.getState(getEnergy=True).getKineticEnergy().value_in_unit(energy_unit)

        meta_bias = 0.0
        if runtime_cache['metadynamics_enabled']:
            group = runtime_cache['metadynamics_group']
            meta_bias_q = context.getState(getEnergy=True, groups=group).getPotentialEnergy()
            meta_bias = meta_bias_q.value_in_unit(energy_unit)

        potential = qm + qm_mm + mm + meta_bias

        temperature = (2.0 * kinetic) / (runtime_cache['dof'] * runtime_cache['gas_constant'])

        total = potential + kinetic

        return {
            'qm': qm,
            'qm_mm': qm_mm,
            'mm': mm,
            'metad_bias':float(meta_bias),
            'potential': potential,
            'kinetic': kinetic,
            'temperature': float(temperature),
            'total': total,
        }

    def _append_step_observables(self, observables):
        """
        Appends per-step observables to trajectory history.
        """

        self.qm_potentials.append(observables['qm'])
        self.qm_mm_interaction_energies.append(observables['qm_mm'])
        self.mm_potentials.append(observables['mm'])
        self.total_potentials.append(observables['potential'])
        self.kinetic_energies.append(observables['kinetic'])
        self.temperatures.append(observables['temperature'])
        self.total_energies.append(observables['total'])

    def _print_run_qmmm_step(self, step, timestep, observables, save_freq):
        """
        Prints the energies of a step, every save_freq steps.

        :param step:
            The step number.
        :param timestep:
            The time step in fs.
        :param observables:
            The observables of the step.
        :param save_freq:
            The printing (and trajectory) frequency in steps.
        """

        if step % save_freq != 0:
            return

        total_std = np.std(self.total_energies) if len(self.total_energies) > 0 else 0.0

        self.ostream.print_blank()
        self.ostream.print_line(
            f"Step: {step} / {self.nsteps}  "
            f"Time: {round((step * timestep) / 1000.0, 2)} ps"
        )
        self.ostream.print_line("-" * 60)
        self.ostream.print_line(
            f"Potential Energy QM region: {observables['qm']:.6f} kJ/mol"
        )
        self.ostream.print_line(
            f"Potential Energy MM region: {observables['mm']:.6f} kJ/mol"
        )
        self.ostream.print_line(
            f"QM/MM Interaction Energy: {observables['qm_mm']:.6f} kJ/mol"
        )
        self.ostream.print_line(
            f"Total Potential Energy: {observables['potential']:.6f} kJ/mol"
        )
        self.ostream.print_line(
            f"Kinetic Energy: {observables['kinetic']:.6f} kJ/mol"
        )
        self.ostream.print_line(
            f"Temperature: {observables['temperature']:.3f} K"
        )
        self.ostream.print_line(
            f"Total Energy: {observables['total']:.6f} +/- {total_std:.6f} kJ/mol"
        )

        if 'metad_bias' in observables:
            self.ostream.print_line(
                f"Metadynamics Bias Energy: {observables['metad_bias']:.6f} kJ/mol"
            )

        for root in self.roots_to_follow:
            self.ostream.print_line(
                f"Current Density root {root}: "
                f"{self.density_around_data_point[root]} -> "
                f"{self.desired_datpoint_density} "
                f"(cycles left: {self.unadded_cycles})"
            )

        self.ostream.print_line(f"Current State (PES): {self.current_state}")
        self.ostream.print_line("-" * 60)
        self.ostream.print_blank()
        self.ostream.flush()

    def _initialize_sampling_impes_drivers(self, inv_sqrt_masses):
        """
        Creates the interpolation drivers of the sampling mirror database (one per
        root) and loads its datapoints in label order.

        :param inv_sqrt_masses:
            The inverse square roots of the Cartesian masses, or None.
        """

        self.sampling_impes_drivers = {root: None for root in self.roots_to_follow}
        self.sampling_qm_data_point_dict = {root: [] for root in self.roots_to_follow}

        for root in self.roots_to_follow:
            drv = InterpolationDriver(self.root_z_matrix[root], ostream=self.ostream)
            drv.update_settings(self.sampling_interpolation_settings[root])
            drv.symmetry_information = self.impes_drivers[root].symmetry_information
            drv.impes_coordinate.inv_sqrt_masses = inv_sqrt_masses

            im_labels, _ = drv.read_labels()
            sorted_im_labels = sorted(im_labels, key=lambda x: int(x.split('_')[1]))
            for label in sorted_im_labels[:]:
                dp = InterpolationDatapoint(self.root_z_matrix[root])
                dp.update_settings(self.sampling_interpolation_settings[root])
                dp.read_hdf5(self.sampling_interpolation_settings[root]['imforcefield_file'], label)

                dp.inv_sqrt_masses = inv_sqrt_masses

                self.sampling_qm_data_point_dict[root].append(dp)

            drv.qm_data_points = self.sampling_qm_data_point_dict[root]

            drv.mark_runtime_data_cache_dirty()
            if len(self.sampling_qm_data_point_dict[root]) > 0:
                drv.impes_coordinate.eq_bond_lengths = self.sampling_qm_data_point_dict[root][0].eq_bond_lengths

            self.sampling_impes_drivers[root] = drv

    def _reload_interpolation_root_from_hdf5(self, root, inv_sqrt_masses):
        """
        Loads all datapoints of one root from the database file in label order
        (point_1, point_2, ...) and hands them to the interpolation driver.

        :param root:
            The root.
        :param inv_sqrt_masses:
            The inverse square roots of the Cartesian masses, or None.
        """

        driver_object = self.impes_drivers[root]

        im_labels, _ = driver_object.read_labels()
        sorted_im_labels = sorted(im_labels, key=lambda x: int(x.split('_')[1]))
        self.qm_data_point_dict[root] = []
        self.sorted_state_spec_im_labels[root] = []

        for label in sorted_im_labels[:]:
            qm_data_point = InterpolationDatapoint(self.root_z_matrix[root])
            qm_data_point.update_settings(self.interpolation_settings[root])
            qm_data_point.read_hdf5(self.interpolation_settings[root]['imforcefield_file'], label)
            qm_data_point.inv_sqrt_masses = inv_sqrt_masses

            self.qm_data_point_dict[root].append(qm_data_point)
            self.sorted_state_spec_im_labels[root].append(label)

        driver_object.qm_data_points = self.qm_data_point_dict[root]
        driver_object.labels = self.sorted_state_spec_im_labels[root]

        if len(self.qm_data_point_dict[root]) > 0:
            driver_object.impes_coordinate.eq_bond_lengths = self.qm_data_point_dict[root][0].eq_bond_lengths

        driver_object.mark_runtime_data_cache_dirty()

    def _check_run_qmmm_scope(self):
        """
        Checks, before anything is set up, that the run is inside the supported
        scope (one MPI rank, ground state only) and that every setting needed by
        run_qmmm is present.
        """

        assert_msg_critical(
            self.nodes == 1,
            'IMDatabasePointCollecter: database construction supports a single '
            'MPI rank only.')

        assert_msg_critical(
            list(self.roots_to_follow) == [0],
            'IMDatabasePointCollecter: only the ground state is supported; '
            f'roots_to_follow must be [0] (got {self.roots_to_follow}). '
            'Excited states and multiple roots are not implemented.')

        assert_msg_critical(
            isinstance(self.drivers, dict) and self.drivers.get('gs') is not None,
            'IMDatabasePointCollecter: dynamics_settings[\'drivers\'][\'gs\'] must '
            'hold the (energy, gradient, Hessian) drivers of the ground state.')

        opt_radius = self.use_opt_confidence_radius
        assert_msg_critical(
            isinstance(opt_radius, (list, tuple)) and len(opt_radius) >= 3,
            'IMDatabasePointCollecter: use_opt_confidence_radius must be a list '
            '[enabled, mode, initial confidence radius(, energy weight)].')
        assert_msg_critical(
            not opt_radius[0] or opt_radius[1] == 'multi_grad',
            'IMDatabasePointCollecter: confidence radius optimization mode '
            f"'{opt_radius[1]}' is not implemented; use 'multi_grad'.")

        assert_msg_critical(
            isinstance(self.non_core_symmetry_groups, dict) and
            'gs' in self.non_core_symmetry_groups,
            'IMDatabasePointCollecter: non_core_symmetry_groups must be a dict '
            "with the key 'gs' (IMForceFieldGenerator.symmetry_information).")

        required = {
            'basis_set_label': (isinstance(self.basis_set_label, dict) and
                                'gs' in self.basis_set_label),
            'duration': self.duration is not None,
            'desired_datapoint_density': self.desired_datpoint_density is not None,
            'converged_cycle': self.unadded_cycles is not None,
            'nsteps': self.nsteps is not None,
            'timestep': self.timestep is not None,
            'ensemble': self.ensemble is not None,
            'energy_threshold': self.energy_threshold is not None,
            'grad_rmsd_thrsh': self.gradient_rmsd_thrsh is not None,
        }
        if self.ensemble in ('NVT', 'NPT'):
            required['temperature'] = self.temperature is not None
            required['friction'] = self.friction is not None
        missing = [key for key, present in required.items() if not present]
        assert_msg_critical(
            not missing,
            'IMDatabasePointCollecter: missing dynamics_settings '
            f'{missing}; set them with update_settings.')

        for name in ('interpolation_settings', 'density_around_data_point',
                     'root_z_matrix'):
            value = getattr(self, name)
            assert_msg_critical(
                isinstance(value, dict) and 0 in value,
                f'IMDatabasePointCollecter: {name} must be a dict with an entry '
                'for root 0 (set by update_settings, by the generator and by '
                'system_from_molecule).')

    def _print_run_qmmm_header(self, integrator, timestep):
        """
        Prints the simulation parameters of run_qmmm.

        :param integrator:
            The OpenMM integrator.
        :param timestep:
            The time step in fs.
        """

        self.ostream.print_blank()
        self.ostream.print_header('QM/MM Simulation Parameters')
        self.ostream.print_header('=' * 60)
        self.ostream.print_info(f'QM driver: {type(self.drivers["gs"][0]).__name__}')
        self.ostream.print_info(f'Ensemble: {self.ensemble}')
        self.ostream.print_info(f'Integration method: {type(integrator).__name__}')
        if self.ensemble in ['NVT', 'NPT']:
            self.ostream.print_info(
                f'Temperature: {self.temperature.value_in_unit(unit.kelvin)} K')
            if self.ensemble == 'NPT':
                pressure = 1.0 if self.pressure is None else float(self.pressure)
                self.ostream.print_info(f'Pressure: {pressure} atm')
            self.ostream.print_info(
                f'Friction: {self.friction.value_in_unit(unit.picosecond**-1)} 1/ps')
        self.ostream.print_info(f'Timestep: {timestep} fs')
        self.ostream.print_info(
            f'Maximum simulation time: {self.nsteps * timestep / 1e6} ns '
            f'({self.nsteps} steps)')
        self.ostream.print_header('=' * 60)
        self.ostream.print_blank()
        self.ostream.flush()

    def run_qmmm(self):
        """
        Runs the database construction: molecular dynamics on the interpolated
        PES with QM checks, datapoint additions and resets, until the database
        holds desired_datapoint_density datapoints or converged_cycle cycles of
        duration steps pass without a new datapoint (or nsteps is reached).

        Requires system_from_molecule, update_settings and the attributes set by
        IMForceFieldGenerator (see the class docstring); the supported scope (one
        MPI rank, roots_to_follow == [0]) and the required settings are checked
        first. run_qmmm can be called again on the same object.

        Writes the trajectory file, checkpoint.xml (after every new datapoint),
        the reference file, pre_opt_database_strucs.xyz and the new datapoints
        (to interpolation_settings[0]['imforcefield_file']). The per-step
        energies are stored in qm_potentials, qm_mm_interaction_energies,
        mm_potentials, total_potentials, kinetic_energies, total_energies (kJ/mol)
        and temperatures (K); the structures of every step in
        state_specific_molecules and the reference structures in allowed_molecules.
        """

        assert_msg_critical('openmm' in sys.modules, 'OpenMM is required for IMDatabasePointCollecter.')

        if self.system is None:
            raise RuntimeError('System has not been created!')

        self._check_run_qmmm_scope()

        # Attach units only once, so that run_qmmm can be called again on the
        # same object (update_settings stores plain numbers)
        if self.temperature is not None and not unit.is_quantity(self.temperature):
            self.temperature = self.temperature * unit.kelvin

        if self.friction is not None and not unit.is_quantity(self.friction):
            self.friction = self.friction / unit.picosecond

        if not unit.is_quantity(self.timestep):
            self.timestep = self.timestep * unit.femtoseconds
        timestep = self.timestep.value_in_unit(unit.femtoseconds)

        self.qm_potentials = []
        self.qm_mm_interaction_energies = []
        self.mm_potentials = []
        self.total_potentials = []
        self.kinetic_energies = []
        self.temperatures = []
        self.total_energies = []

        save_freq = max(1, self.nsteps // self.snapshots) if self.snapshots else max(1, self.nsteps)

        self.topology = self.pdb.topology

        self.positions = self.pdb.positions
        self._topology_atom_labels = [atom.element.symbol for atom in self.topology.atoms()]
        # per-run caches of update_forces (the system may have been rebuilt
        # since an earlier run_qmmm on this object)
        self._qm_external_force = None
        self.previous_candidate_indices = []
        self._gate_bound_cache = None
        self._trust_radius_grouping_cache = None

        # Create or update the integrator
        if self.integrator is None:
            new_integrator = self._create_integrator()
        else:
            new_integrator = self.integrator
        # initializing the metadynamics variables
        self._metadynamics = None
        # a repeated run_qmmm must not keep the bias force of the previous call
        self._remove_metadynamics_bias_force()
        if self.metadynamics_enabled:
            self._initialize_metadynamics()
        self.simulation = app.Simulation(self.topology, self.system, new_integrator, platform=self._create_platform())
        # Load the state if a restart file is provided
        if self.load_system is not None:
            self.simulation.loadState(self.load_system)
        else:
            self.simulation.context.setPositions(self.positions)
            # Set initial velocities if the ensemble is NVT or NPT
            if self.ensemble in ['NVT', 'NPT']:
                self.simulation.context.setVelocitiesToTemperature(self.temperature * 0.5)
        self.start_velocities = self.simulation.context.getState(getVelocities=True).getVelocities()
        # Set up reporting
        self.simulation.reporters.clear()
        self.simulation.reporters.append(app.PDBReporter(self.out_file, save_freq))
        self._print_run_qmmm_header(new_integrator, timestep)

        self.last_point_added = 0
        self.cycle_iteration = self.unadded_cycles

        self.allowed_molecules = {root: {'molecules': [], 'qm_energies': [], 'im_energies':[], 'qm_gradients':[], 'distances': []} for root in self.roots_to_follow}

        for root in self.roots_to_follow:
            self.ostream.print_info(
                f'Datapoints in the database of root {root}: '
                f'{self.density_around_data_point[root]} (target {self.desired_datpoint_density})')
        self.ostream.flush()
        openmm_coordinate = self.simulation.context.getState(getPositions=True).getPositions()
        # coordinates[0] is the reset target after every datapoint addition
        self.coordinates = [openmm_coordinate]

        self.state_specific_molecules = {root: [] for root in self.roots_to_follow}
        self.sorted_state_spec_im_labels = {root: [] for root in self.roots_to_follow}

        self.prev_dens_of_points = {root: 0 for root in self.roots_to_follow}
        self.qm_data_point_dict = {root: [] for root in self.roots_to_follow}
        self.impes_drivers = {root: None for root in self.roots_to_follow}

        # QM reference structures along the trajectory: preload the file once if it
        # exists, otherwise create it; the user-given file name is never replaced
        ref_file = self._reference_file_path()
        qm_labels = self._qm_region_labels()
        self.last_added = 0
        if Path(ref_file).is_file():
            self.allowed_molecules = self.extract_reference_structures_h5(
                ref_file, self.roots_to_follow, expected_labels=qm_labels,
                expected_charge=self.molecule.get_charge(),
                expected_multiplicity=self.molecule.get_multiplicity())
            self.last_added = len(self.allowed_molecules[self.roots_to_follow[0]]['molecules'])
            self.ostream.print_info(
                f'Preloaded {self.last_added} QM reference structures from {ref_file}')
        else:
            self.ostream.print_info(
                f'QM reference structures will be written to the new file {ref_file}')
        self.ostream.flush()
        self._write_reference_file_metadata(ref_file, qm_labels,
                                            self.molecule.get_charge(),
                                            self.molecule.get_multiplicity())

        for root in self.roots_to_follow:
            driver_object = InterpolationDriver(self.root_z_matrix[root], ostream=self.ostream)
            driver_object.update_settings(self.interpolation_settings[root])
            if self.interpolation_settings[root]["use_mass_weight"]:
                masses = self.molecule.get_masses().copy()
                masses_cart = np.repeat(masses, 3)
                inv_sqrt_masses = 1.0 / np.sqrt(masses_cart)
                self.inv_sqrt_masses = inv_sqrt_masses

            driver_object.impes_coordinate.inv_sqrt_masses = self.inv_sqrt_masses
            driver_object.symmetry_information = self.non_core_symmetry_groups['gs']

            self.impes_drivers[root] = driver_object
            self._reload_interpolation_root_from_hdf5(root, self.inv_sqrt_masses)
            self.prev_dens_of_points[root] = len(self.qm_data_point_dict[root])

        if self.sampling_enabled:
            self._initialize_sampling_impes_drivers(self.inv_sqrt_masses)

        self.current_state = self.roots_to_follow[0]
        self.step = 0

        start_time = time()
        runtime_cache = self._build_run_qmmm_runtime_cache()

        for step in range(self.nsteps):

            self.update_forces(self.simulation.context)

            observables = self._collect_step_observables_optimized(
                self.simulation.context,
                runtime_cache)

            self._append_step_observables(observables)
            self._print_run_qmmm_step(step, timestep, observables, save_freq)

            self.step += 1

            if self._metadynamics is not None:
                self._metadynamics.step(self.simulation, 1)
            else:
                self.simulation.step(1)

            reached_target_density = False
            for root in self.roots_to_follow:
                if self.density_around_data_point[root] >= self.desired_datpoint_density:
                    reached_target_density = True
                    break
            if reached_target_density:
                if self.use_opt_confidence_radius[0] and self.confidence_radius_optimized:
                    break
                elif not self.use_opt_confidence_radius[0]:
                    break

            if self.unadded_cycles == 0:
                break

        end_time = time()
        elapsed_time = end_time - start_time
        elapsed_time_days = elapsed_time / (24 * 3600)
        # the loop may stop early (target density or convergence reached)
        steps_run = self.step
        performance = (steps_run * timestep / 1e6) / elapsed_time_days

        self.print_qmmm_summary(elapsed_time, performance, nsteps_run=steps_run)

        return None

    def print_qmmm_summary(self, elapsed_time, performance, nsteps_run=None):
        """
        Prints the summary of run_qmmm.

        :param elapsed_time:
            The wall time in seconds.
        :param performance:
            The performance in ns/day.
        :param nsteps_run:
            The number of steps actually run (default: nsteps).
        """

        if nsteps_run is None:
            nsteps_run = self.nsteps
        self.ostream.print_blank()
        self.ostream.print_header('QM/MM Simulation Completed')
        self.ostream.print_header('==========================')
        self.ostream.print_blank()
        self.ostream.print_header(f"{'Number of steps':<36}{nsteps_run} (maximum {self.nsteps})")
        self.ostream.print_header(f"{'Trajectory file':<36}{self.out_file}")
        self.ostream.print_blank()
        self.ostream.print_header('Simulation Averages')
        self.ostream.print_header('-------------------')
        self.ostream.print_header(
            f"{'QM potential energy':<36}{np.mean(self.qm_potentials): .6f} +/- "
            f"{np.std(self.qm_potentials):.6f} kJ/mol"
        )
        self.ostream.print_header(
            f"{'QM/MM interaction energy':<36}{np.mean(self.qm_mm_interaction_energies): .6f} +/- "
            f"{np.std(self.qm_mm_interaction_energies):.6f} kJ/mol"
        )
        self.ostream.print_header(
            f"{'MM potential energy':<36}{np.mean(self.mm_potentials): .6f} +/- "
            f"{np.std(self.mm_potentials):.6f} kJ/mol"
        )
        self.ostream.print_header(
            f"{'Total potential energy':<36}{np.mean(self.total_potentials): .6f} +/- "
            f"{np.std(self.total_potentials):.6f} kJ/mol"
        )
        self.ostream.print_header(
            f"{'Temperature':<36}{np.mean(self.temperatures): .3f} +/- "
            f"{np.std(self.temperatures):.3f} K"
        )
        self.ostream.print_blank()
        self.ostream.print_info(
            f'Elapsed time: {int(elapsed_time // 60)} minutes, {int(elapsed_time % 60)} seconds'
        )
        self.ostream.print_info(f'Performance: {performance:.2f} ns/day')
        self.ostream.print_blank()
        self.ostream.flush()

    def _create_integrator(self):
        """
        Creates the OpenMM integrator for the ensemble: Langevin for NVT and NPT
        (NPT also adds one MonteCarloBarostat at self.pressure atm, default 1 atm),
        Verlet for NVE.

        :return:
            The OpenMM integrator.
        """

        assert_msg_critical('openmm' in sys.modules, 'OpenMM is required for IMDatabasePointCollecter.')

        # Common parameters for Langevin integrators

        if self.ensemble in ['NVT', 'NPT']:
            integrator = mm.LangevinIntegrator(self.temperature, self.friction, self.timestep)
            integrator.setConstraintTolerance(1e-5)
            if self.ensemble == 'NPT':
                # Add a barostat for pressure control in NPT ensemble;
                # self.pressure is in atm, 1 atm if no pressure was given
                pressure = 1.0
                if self.pressure is not None:
                    pressure = float(self.pressure)
                barostat = None
                for force in self.system.getForces():
                    if isinstance(force, mm.MonteCarloBarostat):
                        barostat = force
                        break
                if barostat is None:
                    barostat = mm.MonteCarloBarostat(pressure * unit.atmospheres,
                                                     self.temperature, 25)
                    self.system.addForce(barostat)
                else:
                    # a repeated run_qmmm reuses the barostat of the first call
                    barostat.setDefaultPressure(pressure * unit.atmospheres)
                    barostat.setDefaultTemperature(self.temperature)
        elif self.ensemble == 'NVE':
            integrator = mm.VerletIntegrator(self.timestep)

        else:
            raise ValueError("Unsupported ensemble type. Please choose 'NVE', 'NVT', or 'NPT'.")

        return integrator

    # Methods to create a QM region/subregion in the system
    def _create_QM_residue(self,
                           ff_gen,
                           qm_atoms,
                           filename='qm_region',
                           residue_name='QMR',
                           write_xml=True):
        """
        This method creates an xml file for a QM region.
        The xml file only contains atomtypes, residue, and nonbonded parameters.

        :param ff_gen:
            VeloxChem forcefield generator object
        :param qm_atoms:
            List of atom indices in the QM region.
        :param filename:
            Name of the files to be generated. Default is 'qm_region'
        :param residue_name:
            Name of the residue. Default is 'QMR'
        :param write_xml:
            If False, only self.qm_atoms is set and no file is written.
        """
        self.qm_atoms = qm_atoms
        if not write_xml:
            return

        atoms = ff_gen.atoms
        bonds = ff_gen.bonds

        # Create the root element of the XML file
        ForceField = ET.Element("ForceField")

        # AtomTypes section
        AtomTypes = ET.SubElement(ForceField, "AtomTypes")

        for i, atom in atoms.items():
            element = ''.join([i for i in atom['name'] if not i.isdigit()])
            attributes = {
                # Name is the atom type_molname
                "name": atom['name'] + '_' + residue_name,
                "class": str(i + 1),
                "element": element,
                "mass": str(atom['mass'])
            }
            ET.SubElement(AtomTypes, "Type", **attributes)

        # Residues section
        Residues = ET.SubElement(ForceField, "Residues")
        Residue = ET.SubElement(Residues, "Residue", name=residue_name)
        for atom_id, atom_data in atoms.items():
            ET.SubElement(Residue, "Atom", name=atom_data['name'], type=atom_data['name'] + '_' + residue_name, charge=str(atom_data['charge']))
        for bond_id, bond_data in bonds.items():
            ET.SubElement(Residue, "Bond", atomName1=atoms[bond_id[0]]['name'], atomName2=atoms[bond_id[1]]['name'])

        # NonbondedForce section
        NonbondedForce = ET.SubElement(ForceField, "NonbondedForce", coulomb14scale=str(ff_gen.fudgeQQ), lj14scale=str(ff_gen.fudgeLJ))
        for atom_id, atom_data in atoms.items():
            attributes = {
                "type": atom_data['name'] + '_' + residue_name,
                "charge": str(0.0),
                "sigma": str(0.0),
                "epsilon": str(0.0)
            }
            ET.SubElement(NonbondedForce, "Atom", **attributes)

        # Generate the tree and write to file
        rough_string = ET.tostring(ForceField, 'utf-8')
        reparsed = minidom.parseString(rough_string)
        indented_string = reparsed.toprettyxml(indent="    ")

        with open(filename + '.xml', 'w') as output_file:
            output_file.write(indented_string)

        self.ostream.print_info(f'QM region parameters written to {filename}.xml')
        self.ostream.flush()

    def _create_QM_subregion(self,
                             ff_gen,
                             qm_atoms,
                             molecule,
                             filename='qm_region',
                             residue_name='QMR',
                             write_xml=True):
        """
        Creates the xml file for a molecule with a QM region.

        :param ff_gen:
            VeloxChem forcefield generator object
        :param qm_atoms:
            List of atom indices in the QM region.
        :param molecule:
            VeloxChem molecule object
        :param filename:
            Name of the files to be generated. Default is 'qm_region'
        :param residue_name:
            Name of the residue. Default is 'QMR'
        :param write_xml:
            If False, only the partition (qm_atoms, mm_subregion,
            linking_atoms) is set and no file is written.
        """

        # Atoms belonging to the QM region do not have bonded parameters
        # Atoms belonging to the MM region have bonded parameters

        # work on a copy: the QM region is extended by the linking atoms below
        qm_atoms = list(qm_atoms)

        # List of atom indices in the molecule
        atom_indices = list(range(molecule.number_of_atoms()))

        connectivity_matrix = molecule.get_connectivity_matrix()

        # Get the atoms connected to the QM region
        qm_connected = []
        for i in qm_atoms:
            qm_connected.extend(np.where(connectivity_matrix[i])[0])
        qm_connected = list(set(qm_connected) - set(qm_atoms))
        self.linking_atoms = qm_connected
        if write_xml:
            self.ostream.print_info(f'Linking atoms indices: {[int(a) for a in qm_connected]}')

        # Get the broken bonds between the QM and MM regions
        self.broken_bonds = []
        for i in qm_connected:
            for j in qm_atoms:
                if connectivity_matrix[i, j]:
                    self.broken_bonds.append((i, j))
        # Set with the atoms in broken bonds
        self.broken_bonds_atoms = list(set([atom for bond in self.broken_bonds for atom in bond]))

        # The linked atoms are part of the QM region since a customforce is used.
        mm_atoms = list(set(atom_indices) - set(qm_atoms) - set(qm_connected))
        # QM atoms are the specified QM region + the linking atoms
        qm_atoms.extend(qm_connected)
        # Order the atoms
        qm_atoms.sort()
        self.qm_atoms = qm_atoms
        self.mm_subregion = mm_atoms
        if write_xml:
            self.ostream.print_info(f'QM subregion (qm_atoms): {[int(a) for a in self.qm_atoms]}')
            self.ostream.print_info(f'MM subregion (mm_subregion): {[int(a) for a in self.mm_subregion]}')
            self.ostream.flush()

        if not write_xml:
            return

        atoms = ff_gen.atoms
        bonds = ff_gen.bonds
        angles = ff_gen.angles
        dihedrals = ff_gen.dihedrals
        impropers = ff_gen.impropers

        # Create the root element of the XML file
        ForceField = ET.Element("ForceField")

        # AtomTypes section
        AtomTypes = ET.SubElement(ForceField, "AtomTypes")

        for i, atom in atoms.items():
            element = ''.join([i for i in atom['name'] if not i.isdigit()])
            attributes = {
                # Name is the atom type_molname
                "name": atom['name'] + '_' + residue_name,
                "class": str(i + 1),
                "element": element,
                "mass": str(atom['mass'])
            }
            ET.SubElement(AtomTypes, "Type", **attributes)

        # Residues section
        Residues = ET.SubElement(ForceField, "Residues")
        Residue = ET.SubElement(Residues, "Residue", name=residue_name)
        for atom_id, atom_data in atoms.items():
            ET.SubElement(Residue, "Atom", name=atom_data['name'], type=atom_data['name'] + '_' + residue_name, charge=str(atom_data['charge']))
        for bond_id, bond_data in bonds.items():
            ET.SubElement(Residue, "Bond", atomName1=atoms[bond_id[0]]['name'], atomName2=atoms[bond_id[1]]['name'])

        # NonbondedForce section
        NonbondedForce = ET.SubElement(ForceField, "NonbondedForce", coulomb14scale=str(ff_gen.fudgeQQ), lj14scale=str(ff_gen.fudgeLJ))
        for atom_id, atom_data in atoms.items():
            if atom_id in qm_atoms:
                charge = 0.0
                sigma = 0.0
                epsilon = 0.0
            else:
                charge = atom_data['charge']
                sigma = atom_data['sigma']
                epsilon = atom_data['epsilon']
            attributes = {
                "type": atom_data['name'] + '_' + residue_name,
                "charge": str(charge),
                "sigma": str(sigma),
                "epsilon": str(epsilon)
            }
            ET.SubElement(NonbondedForce, "Atom", **attributes)

        bonded_atoms = mm_atoms + self.broken_bonds_atoms

        # BondForce section
        BondForce = ET.SubElement(ForceField, "HarmonicBondForce")
        for bond_id, bond_data in bonds.items():
            if (bond_id[0] in bonded_atoms and
                    bond_id[1] in bonded_atoms):
                attributes = {
                    "class1": str(bond_id[0] + 1),
                    "class2": str(bond_id[1] + 1),
                    "length": str(bond_data['equilibrium']),
                    "k": str(bond_data['force_constant'])
                }
                ET.SubElement(BondForce, "Bond", **attributes)

        # AngleForce section
        AngleForce = ET.SubElement(ForceField, "HarmonicAngleForce")
        for angle_id, angle_data in angles.items():
            if (angle_id[0] in bonded_atoms and
                    angle_id[1] in bonded_atoms and
                    angle_id[2] in bonded_atoms):
                attributes = {
                    "class1": str(angle_id[0] + 1),
                    "class2": str(angle_id[1] + 1),
                    "class3": str(angle_id[2] + 1),
                    "angle": str(angle_data['equilibrium'] * np.pi / 180),
                    "k": str(angle_data['force_constant'])
                }
                ET.SubElement(AngleForce, "Angle", **attributes)

        # DihedralForce section
        DihedralForce = ET.SubElement(ForceField, "PeriodicTorsionForce")
        for dihedral_id, dihedral_data in dihedrals.items():
            if (dihedral_id[0] in bonded_atoms and
                    dihedral_id[1] in bonded_atoms and
                    dihedral_id[2] in bonded_atoms and
                    dihedral_id[3] in bonded_atoms):
                # Skip RB dihedrals
                if dihedral_data['type'] == 'RB':
                    continue
                attributes = {
                    "class1": str(dihedral_id[0] + 1),
                    "class2": str(dihedral_id[1] + 1),
                    "class3": str(dihedral_id[2] + 1),
                    "class4": str(dihedral_id[3] + 1),
                    "periodicity1": str(dihedral_data['periodicity']),
                    "phase1": str(dihedral_data['phase'] * np.pi / 180),
                    "k1": str(dihedral_data['barrier'])
                }
                ET.SubElement(DihedralForce, "Proper", **attributes)

        # RB Dihedrals section
        RBForce = ET.SubElement(ForceField, "RBTorsionForce")
        for dihedral_id, dihedral_data in dihedrals.items():
            if (dihedral_id[0] in bonded_atoms and
                    dihedral_id[1] in bonded_atoms and
                    dihedral_id[2] in bonded_atoms and
                    dihedral_id[3] in bonded_atoms):
                # Skip Fourier dihedrals
                if dihedral_data['type'] == 'Fourier':
                    continue
                attributes = {
                    "class1": str(dihedral_id[0] + 1),
                    "class2": str(dihedral_id[1] + 1),
                    "class3": str(dihedral_id[2] + 1),
                    "class4": str(dihedral_id[3] + 1),
                    "c0": str(dihedral_data['RB_coefficients'][0]),
                    "c1": str(dihedral_data['RB_coefficients'][1]),
                    "c2": str(dihedral_data['RB_coefficients'][2]),
                    "c3": str(dihedral_data['RB_coefficients'][3]),
                    "c4": str(dihedral_data['RB_coefficients'][4]),
                    "c5": str(dihedral_data['RB_coefficients'][5])
                }
                ET.SubElement(RBForce, "Torsion", **attributes)

        # ImproperForce section
        ImproperForce = ET.SubElement(ForceField, "PeriodicTorsionForce")  # check the ordering of the improper
        for improper_id, improper_data in impropers.items():
            if (improper_id[0] in bonded_atoms and
                    improper_id[1] in bonded_atoms and
                    improper_id[2] in bonded_atoms and
                    improper_id[3] in bonded_atoms):
                attributes = {
                    "class1": str(improper_id[1] + 1),
                    "class2": str(improper_id[0] + 1),
                    "class3": str(improper_id[2] + 1),
                    "class4": str(improper_id[3] + 1),
                    "periodicity1": str(improper_data['periodicity']),
                    "phase1": str(improper_data['phase'] * np.pi / 180),
                    "k1": str(improper_data['barrier'])
                }
                ET.SubElement(ImproperForce, "Torsion", **attributes)

        # Generate the tree and write to file
        rough_string = ET.tostring(ForceField, 'utf-8')
        reparsed = minidom.parseString(rough_string)
        indented_string = reparsed.toprettyxml(indent="    ")

        with open(filename + '.xml', 'w') as output_file:
            output_file.write(indented_string)

        msg = f'QM subregion parameters written to {filename}.xml'
        self.ostream.print_info(msg)
        self.ostream.flush()

    # Method to set a qm_mm system:
    def set_qm_mm_system(self, phase, ff_gen):
        """
        Adds the external force that carries the interpolated QM forces to the
        system (on self.qm_atoms) and, if there are MM atoms, the QM/MM van der
        Waals and Coulomb interactions; sets the force groups.

        :param phase:
            Phase of the system: 'gas', or 'periodic' (also 'water') for a
            periodic, solvated system.
        :param ff_gen:
            MMForceFieldGenerator object from VeloxChem.
        """

        assert_msg_critical('openmm' in sys.modules, 'OpenMM is required for IMDatabasePointCollecter.')

        # Set the QM/MM Interaction Groups
        total_atoms = self.system.getNumParticles()

        # The MM subregion is counted as regular MM atoms
        qm_group = set(self.qm_atoms)

        mm_group = set(range(total_atoms)) - qm_group

        if not mm_group:
            self.ostream.print_info('No external MM atoms found in the system')
            self.ostream.flush()

        # Add custom hessian forces
        force_expression = mm.CustomExternalForce("-fx*x-fy*y-fz*z")
        self.system.addForce(force_expression)
        force_expression.addPerParticleParameter("fx")
        force_expression.addPerParticleParameter("fy")
        force_expression.addPerParticleParameter("fz")

        for i in self.qm_atoms:
            force_expression.addParticle(i, [0, 0, 0])

        # QM Hessian Force Group
        force_expression.setForceGroup(0)

        # If a MM region is present define the interactions
        if mm_group:

            nonbonded_force = None
            for force in self.system.getForces():
                if isinstance(force, mm.NonbondedForce):
                    nonbonded_force = force
                    break

            if nonbonded_force is None:
                raise ValueError("NonbondedForce not found in the system")

            # CustomNonbondedForce for QM/MM interactions
            vdw = mm.CustomNonbondedForce("4*epsilon*((sigma/r)^12-(sigma/r)^6); sigma=0.5*(sigma1+sigma2); epsilon=sqrt(epsilon1*epsilon2)")
            if phase in ['water', 'periodic']:
                vdw.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
                vdw.setCutoffDistance(self.cutoff)
            vdw.addPerParticleParameter("sigma")
            vdw.addPerParticleParameter("epsilon")

            if phase in ['water', 'periodic']:
                # OpenMM uses Reaction Field method for CustomNB with PBC.
                rfDielectric = nonbonded_force.getReactionFieldDielectric()
                krf = (1 / (self.cutoff**3)) * (rfDielectric - 1) / (2 * rfDielectric + 1)
                crf = (1 / self.cutoff) * (3 * rfDielectric) / (2 * rfDielectric + 1)
                coulomb_rf = f"(138.935456*charge1*charge2)*(1/r + {krf}*r*r - {crf});"
                coulomb = mm.CustomNonbondedForce(coulomb_rf)
                coulomb.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
                coulomb.setCutoffDistance(self.cutoff)
                # Disable long-range electrostatics for nonbonded force
                nonbonded_force.setUseDispersionCorrection(False)
            else:
                coulomb = mm.CustomNonbondedForce("138.935456*charge1*charge2/r;")
                coulomb.setNonbondedMethod(mm.CustomNonbondedForce.NoCutoff)

            coulomb.addPerParticleParameter("charge")

            # Apply the same exclusions to the custom forces as in the NonbondedForce
            # This is needed to avoid errors (all forces must have identical exclusions)
            for i in range(nonbonded_force.getNumExceptions()):
                p1, p2, chargeProd, sigma, epsilon = nonbonded_force.getExceptionParameters(i)
                vdw.addExclusion(p1, p2)
                coulomb.addExclusion(p1, p2)

            self.system.addForce(vdw)
            self.system.addForce(coulomb)

            # Add particles to the custom forces
            # QM region
            for i in qm_group:
                vdw.addParticle([ff_gen.atoms[i]['sigma'] * unit.nanometer,
                                 ff_gen.atoms[i]['epsilon']] * unit.kilojoules_per_mole)
                coulomb.addParticle([ff_gen.atoms[i]['charge']] * unit.elementary_charge)

            # MM region
            # Obtain the sigma, epsilon, and charge values from the system
            nonbonded_force = [f for f in self.system.getForces() if isinstance(f, mm.NonbondedForce)][0]

            for i in mm_group:
                # The charges, sigmas, and epsilons are taken from the system
                charge, sigma, epsilon = nonbonded_force.getParticleParameters(i)
                sigma = sigma * unit.nanometer
                epsilon = epsilon * unit.kilojoules_per_mole
                charge = charge * unit.elementary_charge
                vdw.addParticle([sigma, epsilon])
                coulomb.addParticle([charge])

            vdw.addInteractionGroup(qm_group, mm_group)
            coulomb.addInteractionGroup(qm_group, mm_group)

            # Set force groups
            vdw.setForceGroup(1)
            coulomb.setForceGroup(2)

            # Set a force group for the MM region
            for force in self.system.getForces():
                # Non bonded MM
                if isinstance(force, mm.NonbondedForce):
                    force.setForceGroup(3)
                # Bonded
                elif isinstance(force, mm.HarmonicBondForce):
                    force.setForceGroup(4)
                elif isinstance(force, mm.HarmonicAngleForce):
                    force.setForceGroup(5)
                elif isinstance(force, mm.PeriodicTorsionForce):
                    force.setForceGroup(6)
                elif isinstance(force, mm.RBTorsionForce):
                    force.setForceGroup(7)
                elif isinstance(force, mm.CMMotionRemover):
                    force.setForceGroup(8)

        # Determine the force index for the QM region
        # it is an instance of openmm.openmm.CustomExternalForce
        for i, force in enumerate(self.system.getForces()):
            if isinstance(force, mm.CustomExternalForce):
                self.qm_force_index = i
                break

    def cartesian_just_distance(self, coordinate_1, coordinate_2, non_core_atoms):
        """
        Returns the Cartesian distance between two structures after removing the
        given atoms, moving both centroids to the origin and rotating coordinate_1
        onto coordinate_2.

        :param coordinate_1:
            The coordinates of the first structure, shape (natoms, 3).
        :param coordinate_2:
            The coordinates of the second structure, shape (natoms, 3).
        :param non_core_atoms:
            The indices of the atoms that are ignored.

        :return:
            The norm of the difference vector (in the unit of the coordinates).
        """
        target_coordinates_core = np.delete(coordinate_1, non_core_atoms, axis=0)
        reference_coordinates_core = np.delete(coordinate_2, non_core_atoms, axis=0)
        # First, translate the cartesian coordinates to zero
        target_coordinates = self.calculate_translation_coordinates(target_coordinates_core)
        reference_coordinates = (
            self.calculate_translation_coordinates(reference_coordinates_core))
        # Then, determine the rotation matrix which
        # aligns data_point (target_coordinates)
        # to self.impes_coordinate (reference_coordinates)
        rotation_matrix = geometric.rotate.get_rot(target_coordinates,
                                                   reference_coordinates)
        # Rotate the data point
        rotated_coordinates = np.dot(rotation_matrix, target_coordinates.T).T
        # Calculate the Cartesian distance
        distance_vector = (reference_coordinates - rotated_coordinates)

        return np.linalg.norm(distance_vector)

    def _reference_rmsd_lower_bounds(self, molecule_list, new_coords):
        """
        Returns a lower bound of the RMSD (angstrom) of the proximity gate
        between new_coords and every structure of molecule_list:

            RMSD >= sqrt(sum_{a<b} (d_ab - d'_ab)^2) / N,

        with d_ab and d'_ab the interatomic distances of the two structures
        and N the number of atoms. After any rigid superposition with
        displacements delta_a, sum_{a<b} |delta_a - delta_b|^2
        <= N sum_a |delta_a|^2 and |d_ab - d'_ab| <= |delta_a - delta_b|.
        The interatomic distances of the references are cached; the list must
        only grow (a new list or a changed entry resets the cache).

        :param molecule_list:
            The reference structures.
        :param new_coords:
            The coordinates of the new structure in bohr, (natoms, 3).

        :return:
            The lower bounds, one per reference.
        """

        n_ref = len(molecule_list)
        n_atoms = new_coords.shape[0]
        pair_i, pair_j = np.triu_indices(n_atoms, 1)

        def pair_distances(coords):
            diff = coords[pair_i] - coords[pair_j]
            return np.sqrt(np.einsum('ij,ij->i', diff, diff))

        cache = self._gate_bound_cache
        if (cache is None or cache['list'] is not molecule_list or cache['n'] > n_ref
                or (cache['n'] > 0 and cache['last'] is not molecule_list[cache['n'] - 1])):
            cache = {'list': molecule_list, 'n': 0, 'last': None,
                     'distances': np.zeros((max(n_ref, 64), pair_i.size))}
            self._gate_bound_cache = cache

        if n_ref > cache['distances'].shape[0]:
            grown = np.zeros((max(n_ref, 2 * cache['distances'].shape[0]), pair_i.size))
            grown[:cache['n']] = cache['distances'][:cache['n']]
            cache['distances'] = grown
        for idx in range(cache['n'], n_ref):
            cache['distances'][idx] = pair_distances(molecule_list[idx].get_coordinates_in_bohr())
        if n_ref > cache['n']:
            cache['n'] = n_ref
            cache['last'] = molecule_list[n_ref - 1]

        new_distances = pair_distances(new_coords)
        sum_sq = np.zeros(n_ref)
        chunk = max(1, 2**20 // max(pair_i.size, 1))
        for start in range(0, n_ref, chunk):
            delta = cache['distances'][start:min(start + chunk, n_ref)] - new_distances
            sum_sq[start:start + chunk] = np.einsum('ij,ij->i', delta, delta)

        return np.sqrt(sum_sq) / n_atoms * bohr_in_angstrom()

    def _bounded_reference_scan(self, molecule_list, new_coords, n_atoms, threshold, K, closest_indices):
        """
        Exact replacement of the full scan of the proximity gate
        (gate_distance_prefilter). Only references that the lower bound of
        _reference_rmsd_lower_bounds does not exclude are aligned:
        - the first reference within threshold, in list order, is found among
          those with a bound <= threshold;
        - of the scanned references (up to and including that one, or all),
          the ones that can be among the K closest are aligned in order of
          increasing bound.

        closest_indices (the entries of the cached candidates) is extended in
        place by the aligned references in list order, so that its K closest
        entries, after the stable sort of update_forces, are those of the
        full scan.

        :return:
            True if a reference lies within threshold.
        """

        lower = self._reference_rmsd_lower_bounds(molecule_list, new_coords)
        # relative margin far above the rounding of the bound and of the RMSD
        margin = 1.0 + 1.0e-10

        known = {idx: dist for idx, dist in closest_indices}

        def rmsd(idx):
            if idx not in known:
                checked_coords = molecule_list[idx].get_coordinates_in_bohr()
                checked_distance = self.cartesian_just_distance(checked_coords, new_coords, [])
                known[idx] = (np.linalg.norm(checked_distance) / np.sqrt(n_atoms)) * bohr_in_angstrom()
            return known[idx]

        first = None
        for idx in np.flatnonzero(lower <= threshold * margin):
            if rmsd(int(idx)) <= threshold:
                first = int(idx)
                break
        n_scanned = len(molecule_list) if first is None else first + 1

        # the scanned entries of the full scan with a known RMSD, and the
        # K smallest values among all entries of the full scan known so far
        scanned_known = {idx for idx in known if idx < n_scanned}
        values = sorted([dist for _, dist in closest_indices] + [known[idx] for idx in scanned_known])

        for idx in np.argsort(lower[:n_scanned], kind='stable'):
            idx = int(idx)
            if idx in scanned_known:
                continue
            if len(values) >= K and lower[idx] > values[K - 1] * margin:
                break
            scanned_known.add(idx)
            dist = rmsd(idx)
            values.insert(int(np.searchsorted(values, dist, side='right')), dist)

        closest_indices.extend((idx, known[idx]) for idx in sorted(scanned_known))

        return first is not None

    def update_gradient_and_energy(self, new_positions):
        """
        Evaluates the interpolation at the given QM positions.

        :param new_positions:
            The positions of the QM atoms in nm.

        :return:
            The gradient (Hartree/bohr, shape (natoms, 3)) and the energy (kJ/mol).
        """

        new_molecule = self._build_qm_molecule_from_positions(new_positions)
        self._latest_qm_molecule = new_molecule
        self._latest_qm_positions_nm = np.asarray(new_positions, dtype=np.float64).copy()

        self._compute_all_roots_for_molecule(new_molecule)

        # hartree_in_kcalpermol() * 4.184 differs from hartree_in_kjpermol() in the
        # last digits; it is kept so that stored runs are reproduced bitwise
        potential_kjmol = self.impes_drivers[self.current_state].impes_coordinate.energy * hartree_in_kcalpermol() * 4.184
        self.current_gradient = self.impes_drivers[self.current_state].impes_coordinate.gradient
        self.current_energy = potential_kjmol
        return self.current_gradient, potential_kjmol

    def _qm_region_labels(self):
        """
        Returns the element labels of the QM molecule that interpolation and QM
        see: the QM atoms of the topology, with linking atoms replaced by H.

        :return:
            The list of element labels.
        """

        atom_labels = self._topology_atom_labels
        if atom_labels is None:
            atom_labels = [atom.element.symbol for atom in self.topology.atoms()]
            self._topology_atom_labels = atom_labels

        if self.mm_subregion is not None:
            # Change the linking atoms to H
            return ['H' if i in self.linking_atoms else atom_labels[i] for i in self.qm_atoms]
        return [atom_labels[i] for i in self.qm_atoms]

    def _build_qm_molecule_from_positions(self, positions):
        """
        Builds the QM molecule (charge and multiplicity of self.molecule) from QM
        positions; with a QM/MM partition, linking atoms become H atoms at
        linking_atom_distance.

        :param positions:
            The positions of the QM atoms in nm.

        :return:
            The molecule.
        """

        positions_ang = positions * 10
        qm_atom_labels = self._qm_region_labels()

        new_molecule = None
        # Check if there is a QM/MM partition in the system
        if self.mm_subregion is not None:
            # Create a molecule with a link atom (H)
            # The linking atom is defined in self.linking_atoms
            # It need to be changed to a H atom at 1.0 angstrom
            # The QM region is defined in self.qm_atoms

            # Change the positions of the linking atoms to 1.0 angstrom
            for atom1, atom2 in self.broken_bonds:

                current_distance = np.linalg.norm(positions_ang[atom1] - positions_ang[atom2])
                # Change the position of the linking atom to 1.0 angstrom
                # By construction, the first atom is the linking atom
                direction = (positions_ang[atom2] - positions_ang[atom1]) / current_distance
                positions_ang[atom1] = positions_ang[atom2] - direction * self.linking_atom_distance

            new_molecule = Molecule(qm_atom_labels, positions_ang, units="angstrom")
            new_molecule.set_charge(self.molecule.get_charge())
            new_molecule.set_multiplicity(self.molecule.get_multiplicity())

        else:
            new_molecule = Molecule(qm_atom_labels, positions_ang, units="angstrom")
            new_molecule.set_charge(self.molecule.get_charge())
            new_molecule.set_multiplicity(self.molecule.get_multiplicity())

        return new_molecule

    def _extract_qm_positions_nm(self, context):
        """
        Returns QM coordinates (nm) from OpenMM context as float64 array.
        """

        state = context.getState(getPositions=True)
        positions_nm = state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
        return np.asarray(positions_nm[self.qm_atoms], dtype=np.float64)

    def _get_latest_qm_molecule(self, qm_positions_nm=None):
        """
        Returns the most recent QM molecule used in interpolation.
        Falls back to building from qm_positions_nm if cache is unavailable.
        """

        molecule = self._latest_qm_molecule
        if molecule is not None:
            if qm_positions_nm is None:
                return molecule
            cached_positions = self._latest_qm_positions_nm
            current_positions = np.asarray(qm_positions_nm, dtype=np.float64)
            if cached_positions is not None and cached_positions.shape == current_positions.shape:
                if np.allclose(cached_positions, current_positions, atol=1.0e-12, rtol=0.0):
                    return molecule

        if qm_positions_nm is None:
            return None

        return self._build_qm_molecule_from_positions(
            np.asarray(qm_positions_nm, dtype=np.float64))

    def _compute_all_roots_for_molecule(self, molecule):
        """
        Evaluates the interpolation of every followed root at a molecule.

        :param molecule:
            The molecule.
        """

        for root in self.roots_to_follow:
            driver = self.impes_drivers[root]

            if driver.qm_data_points is not self.qm_data_point_dict[root]:
                driver.qm_data_points = self.qm_data_point_dict[root]

            driver.compute(molecule)

    def update_gradient(self, new_positions):
        """
        Updates and returns the gradient of the QM region.

        :param new_positions:
            The new positions of the atoms in the QM region.
        :return:
            The gradient of the QM region.
        """
        gradient, _ = self.update_gradient_and_energy(
            new_positions,
        )

        return gradient

    def update_forces(self, context):
        """
        Evaluates the interpolation at the current positions, applies the forces
        to the QM atoms and decides whether a QM check is needed (proximity gate);
        runs the QM check and, after an addition or a reference collection,
        resets positions and velocities to the start of the trajectory.

        :param context:
            The OpenMM context object.
        """

        assert_msg_critical('openmm' in sys.modules, 'OpenMM is required for IMDatabasePointCollecter.')

        if self._gradient_to_force_factor is None:
            # Hartree/bohr -> kJ/mol/nm (4.184 kept for bitwise reproducibility)
            self._gradient_to_force_factor = (
                (4.184 * hartree_in_kcalpermol() * 10.0 / bohr_in_angstrom())
                * unit.kilojoule_per_mole / unit.nanometer
            )
        conversion_factor = self._gradient_to_force_factor

        qm_force = self._qm_external_force
        if qm_force is None:
            qm_force = self.system.getForce(self.qm_force_index)
            self._qm_external_force = qm_force

        # Update the forces of the QM region.
        qm_positions = self._extract_qm_positions_nm(context)

        gradient = self.update_gradient(qm_positions)

        # Reuse the exact molecule passed through interpolation compute when available.
        new_molecule = self._get_latest_qm_molecule(qm_positions)

        force = -np.array(gradient) * conversion_factor

        # Correlation Check
        self.add_a_point = False

        if self.skipping_value == 0 and self.point_checker + 1 > self.start_collect:
            bond_rmsd_values = np.asarray(self.impes_drivers[self.current_state].bond_rmsd, dtype=np.float64)
            angle_rmsd_values = np.asarray(self.impes_drivers[self.current_state].angle_rmsd, dtype=np.float64)
            dihedral_rmsd_values = np.asarray(self.impes_drivers[self.current_state].dihedral_rmsd, dtype=np.float64)
            mean_bond_rmsd = float(np.mean(bond_rmsd_values)) if bond_rmsd_values.size else 0.0
            mean_angle_rmsd = float(np.mean(angle_rmsd_values)) if angle_rmsd_values.size else 0.0
            mean_dihedral_rmsd = float(np.mean(dihedral_rmsd_values)) if dihedral_rmsd_values.size else 0.0

            if mean_bond_rmsd > 0.02 or mean_angle_rmsd > 0.1 or mean_dihedral_rmsd > 1.5:
                K = 5  # Number of closest previous matches to cache
                threshold = 5e-2  # RMS distance per atom in angstrom
                new_coords = new_molecule.get_coordinates_in_bohr()
                n_atoms = len(new_molecule.get_labels())
                molecule_list = self.allowed_molecules[self.current_state]['molecules']

                scanned = False
                closest_indices = []
                index_added = []
                # --- 1. Try cached candidates first
                for idx in self.previous_candidate_indices:
                    if idx >= len(molecule_list):
                        continue
                    checked_molecule = molecule_list[idx]
                    checked_coords = checked_molecule.get_coordinates_in_bohr()
                    checked_distance = self.cartesian_just_distance(checked_coords, new_coords, [])
                    normed_dist = (np.linalg.norm(checked_distance) / np.sqrt(n_atoms)) * bohr_in_angstrom()
                    if idx not in index_added:
                        closest_indices.append((idx, normed_dist))
                        index_added.append(idx)
                    if normed_dist <= threshold:
                        scanned = True
                        break
                # --- 2. Fallback: full scan only if no match
                if not scanned and self.gate_distance_prefilter:
                    scanned = self._bounded_reference_scan(
                        molecule_list, new_coords, n_atoms, threshold, K, closest_indices)
                elif not scanned:
                    for idx, checked_molecule in enumerate(molecule_list):
                        checked_coords = checked_molecule.get_coordinates_in_bohr()
                        checked_distance = self.cartesian_just_distance(checked_coords, new_coords, [])
                        normed_dist = (np.linalg.norm(checked_distance) / np.sqrt(n_atoms)) * bohr_in_angstrom()
                        closest_indices.append((idx, normed_dist))
                        if normed_dist <= threshold:
                            scanned = True
                            break
                # --- 3. Update cached candidate indices for the next step
                closest_indices.sort(key=lambda x: x[1])  # Sort by distance
                self.previous_candidate_indices = [i for i, _ in closest_indices[:K]]
                if not scanned:
                    for i, qm_data_point in enumerate(self.qm_data_point_dict[self.current_state], start=1):
                        length_vectors = (self.impes_drivers[self.current_state].impes_coordinate.cartesian_distance_vector(qm_data_point))
                        if (np.linalg.norm(length_vectors) / np.sqrt(len(self.molecule.get_labels()))) * bohr_in_angstrom() <= 5e-2:
                            self.add_a_point = False
                            break
                        if i == len(self.qm_data_point_dict[self.current_state]):
                            self.add_a_point = True
        else:
            self.skipping_value -= 1
            if self.skipping_value < 0:
                self.skipping_value = 0

        self.point_checker += 1

        if self.add_a_point:

            self.point_correlation_check(new_molecule)

        if self.point_checker == 0:

            self.current_state = self.starting_state
            context.setPositions(self.coordinates[0])
            if self._metadynamics is not None and self.metadynamics_enabled:
                self._reset_metadynamics_bias(context)
            # Set initial velocities if the ensemble is NVT or NPT
            if self.ensemble in ['NVT', 'NPT']:
                context.setVelocitiesToTemperature(self.temperature * 0.5)
            else:
                context.setVelocities(self.start_velocities)

            qm_positions = self._extract_qm_positions_nm(context)
            gradient_2 = self.update_gradient(qm_positions)
            new_molecule = self._get_latest_qm_molecule(qm_positions)

            force = -np.array(gradient_2) * conversion_factor

        if (self.point_checker + self.last_point_added) % self.duration == 0 and self.point_checker != 0:
            self.last_point_added = 0
            self.unadded_cycles -= 1
        if self.point_checker < 500 and self.cycle_iteration != self.unadded_cycles:
            self.unadded_cycles += 1

        self.state_specific_molecules[self.current_state].append(new_molecule)
        for i, atom_idx in enumerate(self.qm_atoms):

            qm_force.setParticleParameters(i, atom_idx, force[i])
        qm_force.updateParametersInContext(context)

    # Functions to expand the database

    def _print_point_correlation_result(
        self,
        root,
        energy_difference,
        energy_difference_per_atom,
        rmsd_gradient,
    ):
        """
        Prints the interpolation error found by a QM check.

        :param root:
            The root.
        :param energy_difference:
            The absolute energy difference in Hartree.
        :param energy_difference_per_atom:
            The energy difference per atom in kcal/mol.
        :param rmsd_gradient:
            The gradient difference RMS (over all 3N components) in kcal/mol/angstrom.
        """

        self.ostream.print_line(f"State {root}")
        self.ostream.print_line(
            f"  Energy difference:        {energy_difference: .8e} Eh"
        )
        self.ostream.print_line(
            f"  Energy difference:        {energy_difference * hartree_in_kcalpermol(): .6f} kcal/mol"
        )
        self.ostream.print_line(
            f"  Energy difference/atom:   {energy_difference_per_atom: .6f} kcal/mol"
        )
        self.ostream.print_line(
            f"  Gradient RMSD:            {rmsd_gradient: .6f} kcal/mol/angstrom"
        )

    def _print_datapoint_decision(self, add_point, kept_as_reference=False):
        """
        Prints the outcome of a QM check, after the decision has been taken.

        :param add_point:
            True if the structure is added to the database.
        :param kept_as_reference:
            True if the thresholds were exceeded, but the structure was only
            collected as a reference for the confidence radius fit.
        """

        decision = 'yes' if add_point else 'no'
        if kept_as_reference:
            decision += ' (threshold exceeded; kept as reference for the confidence radius fit)'
        self.ostream.print_line(f"  Add datapoint:            {decision}")
        self.ostream.flush()

    def _sampling_screen_check(self, molecule):
        """
        Optional xTB pre-screen of a QM check: compares the xTB energy with the
        interpolation of the sampling mirror database. If the difference per
        atom is within sampling_settings['e_thrsh_kcal_per_atom'], the QM check is
        skipped and the next 1 to 20 steps are not checked. Only the energy is
        compared; other sampling_settings keys are not read.

        :param molecule:
            The current structure.

        :return:
            A dictionary with 'skip_full_qm' (bool) and 'details' (per root).
        """

        if not self.sampling_enabled or self.sampling_driver is None or self.sampling_impes_drivers is None:
            return {"skip_full_qm": False, "details": {}}

        natms = len(molecule.get_labels())

        sampling_qm, sampling_grad, _ = self.sampling_driver['gs']

        xtb_energy, _, _ = self._compute_energy(sampling_qm, molecule, basis=None)

        e_thr = float(self.sampling_settings["e_thrsh_kcal_per_atom"])

        all_ok = True
        details = {}

        for root in self.roots_to_follow:
            self.sampling_impes_drivers[root].compute(molecule)
            e_im = self.sampling_impes_drivers[root].impes_coordinate.energy

            e_diff = abs(float(xtb_energy[0]) - float(e_im)) / natms * hartree_in_kcalpermol()

            details[root] = {"e_diff_kcal_per_atom": e_diff}
            if not (e_diff <= e_thr):
                all_ok = False

        if all_ok:
            # Use the worst root as the decision value.
            max_e_diff = max(
                float(root_details["e_diff_kcal_per_atom"])
                for root_details in details.values()
            )

            # 0 means perfect agreement, 1 means exactly at threshold.
            rel = max(0.0, min(max_e_diff / e_thr, 1.0))

            # Quadratic falloff: conservative near threshold, generous near zero.
            skip = 20 - (20 - 1) * rel**2
            self.skipping_value = skip

            return {"skip_full_qm": True, "details": details}

        return {"skip_full_qm": False, "details": details}

    def point_correlation_check(self, molecule):
        """
        Runs the QM check at the current structure: computes the QM energy and
        gradient and compares them with the interpolation. If the energy
        difference per atom exceeds energy_threshold or the gradient RMSD
        exceeds grad_rmsd_thrsh, the structure is added to the database (after a
        constrained re-optimisation when identify_relevant_int_coordinates[0] is
        True) or, during the first confidence-radius references, kept as a
        reference; otherwise it is stored in the reference file. With sampling
        enabled, an xTB pre-screen may skip the check.

        :param molecule:
            The molecule corresponding to the current conformation
            in the IM dynamics.
        """

        if self.sampling_enabled:
            screen = self._sampling_screen_check(molecule)
            if screen["skip_full_qm"]:
                self.add_a_point = False
                return

        self.ostream.print_header("Performing point correlation check for current conformation with QM.")
        self.ostream.flush()
        self.add_a_point = False

        # one root, ground state (checked at the start of run_qmmm)
        root = self.roots_to_follow[0]
        current_state_difference = {root: [None, None]}
        state_specific_energies = {root: None}
        state_specific_gradients = {root: None}
        addition_of_state_specific_points = []

        drivers = self.drivers['gs']
        natms = len(molecule.get_labels())
        self.ostream.print_blank()
        self.ostream.print_info("Initialize QM calculation")
        self.ostream.flush()
        current_basis = MolecularBasis.read(molecule, self.basis_set_label['gs'])

        qm_energy, scf_results, rsp_results = self._compute_energy(drivers[0], molecule, current_basis)
        gradients = self._compute_gradient(drivers[1], molecule, current_basis, scf_results, rsp_results)

        impes_coordinate = self.impes_drivers[root].impes_coordinate
        grad = gradients[0]
        energy_difference = abs(qm_energy[0] - impes_coordinate.energy)
        energy_difference_per_atom = energy_difference / natms * hartree_in_kcalpermol()

        gradient_difference = (grad - impes_coordinate.gradient) * hartree_in_kcalpermol() / bohr_in_angstrom()
        rmsd_gradient = np.sqrt((gradient_difference**2).mean())

        self._print_point_correlation_result(
            root=root,
            energy_difference=energy_difference,
            energy_difference_per_atom=energy_difference_per_atom,
            rmsd_gradient=rmsd_gradient,
        )

        state_specific_energies[root] = [qm_energy[0], impes_coordinate.energy]
        current_state_difference[root][0] = energy_difference_per_atom
        current_state_difference[root][1] = rmsd_gradient
        state_specific_gradients[root] = [grad, impes_coordinate.gradient]

        energy_exceeded = current_state_difference[self.current_state][0] > self.energy_threshold
        gradient_exceeded = current_state_difference[self.current_state][1] > self.gradient_rmsd_thrsh
        # with the confidence radius fit on, points are added only once the radii are fitted
        radii_ready = not self.use_opt_confidence_radius[0] or self.confidence_radius_optimized
        if (energy_exceeded or gradient_exceeded) and radii_ready:
            addition_of_state_specific_points.append(self.current_state)

        if len(addition_of_state_specific_points) > 0:

            for root in self.roots_to_follow:
                if self.use_opt_confidence_radius[0] and len(self.allowed_molecules[root]['molecules']) < 20:
                    self.allowed_molecules[root]['molecules'].append(molecule)
                    self.allowed_molecules[root]['im_energies'].append(self.impes_drivers[root].impes_coordinate.energy)
                    self.allowed_molecules[root]['qm_energies'].append(state_specific_energies[root][0])
                    self.allowed_molecules[root]['qm_gradients'].append(state_specific_gradients[root][0])
                    self.last_point_added = self.point_checker - 1
                    self.point_checker = 0
                    self.add_a_point = False

                else:
                    self.add_a_point = True

            self._print_datapoint_decision(self.add_a_point,
                                           kept_as_reference=not self.add_a_point)

        else:
            self._print_datapoint_decision(False)
            for root in self.roots_to_follow:
                self.allowed_molecules[root]['molecules'].append(molecule)
                self.allowed_molecules[root]['im_energies'].append(self.impes_drivers[root].impes_coordinate.energy)
                self.allowed_molecules[root]['qm_energies'].append(state_specific_energies[root][0])
                self.allowed_molecules[root]['qm_gradients'].append(state_specific_gradients[root][0])

                self.write_qm_energy_determined_points_h5(
                    self.allowed_molecules[root]['molecules'][self.last_added:],
                    self.allowed_molecules[root]['qm_energies'][self.last_added:],
                    self.allowed_molecules[root]['qm_gradients'][self.last_added:],
                    self.allowed_molecules[root]['im_energies'][self.last_added:],
                    root)
                self.last_added = len(self.allowed_molecules[root]['molecules'])

            if self.use_opt_confidence_radius[0] and len(self.allowed_molecules[self.current_state]['molecules']) >= 10 and self.density_around_data_point[self.current_state] > 1 and self.prev_dens_of_points[self.current_state] != self.density_around_data_point[self.current_state]:
                self.prev_dens_of_points[self.current_state] = self.density_around_data_point[self.current_state]
                sym_dict = self.non_core_symmetry_groups['gs']

                indices = [i for i in range(len(self.allowed_molecules[self.current_state]['molecules']) - 1)]
                chosen_structures = [self.allowed_molecules[self.current_state]['molecules'][idx] for idx in indices]
                chosen_qm_energies = [self.allowed_molecules[self.current_state]['qm_energies'][idx] for idx in indices]
                chosen_im_energies = [self.allowed_molecules[self.current_state]['im_energies'][idx] for idx in indices]
                chosen_qm_gradients = [self.allowed_molecules[self.current_state]['qm_gradients'][idx] for idx in indices]

                # the mode is 'multi_grad' (checked at the start of run_qmmm)
                trust_radii = self.determine_trust_radius_gradient(
                    chosen_structures,
                    chosen_qm_energies,
                    chosen_qm_gradients,
                    chosen_im_energies,
                    self.qm_data_point_dict[self.current_state],
                    self.interpolation_settings[self.current_state],
                    sym_dict,
                    self.root_z_matrix[self.current_state],
                    exponent_p_q=(self.impes_drivers[self.current_state].exponent_p, self.impes_drivers[self.current_state].exponent_q))

                for idx, trust_radius in enumerate(trust_radii):

                    self.qm_data_point_dict[self.current_state][idx].update_confidence_radius(self.interpolation_settings[self.current_state]['imforcefield_file'], self.sorted_state_spec_im_labels[self.current_state][idx], trust_radius)
                    self.qm_data_point_dict[self.current_state][idx].confidence_radius = trust_radius
                    if self.sampling_enabled:
                        self.sampling_qm_data_point_dict[self.current_state][idx].update_confidence_radius(self.sampling_interpolation_settings[self.current_state]['imforcefield_file'], self.sorted_state_spec_im_labels[self.current_state][idx], trust_radius)
                        self.sampling_qm_data_point_dict[self.current_state][idx].confidence_radius = trust_radius

                self._refresh_interpolation_driver_caches(self.current_state)

                self.confidence_radius_optimized = True
                self.add_a_point = False

        # no adaptive skipping after a QM check: the proximity gate is evaluated
        # again at the next step
        self.skipping_value = 0
        if self.add_a_point:
            self.confidence_radius_optimized = False
            self.ostream.print_blank()
            self.ostream.print_header("Database expansion triggered! --> Initialization of expansion framework.")
            self.ostream.print_blank()
            self.ostream.print_block(molecule.get_xyz_string())
            self.ostream.flush()
            with open('pre_opt_database_strucs.xyz', "a", encoding="utf-8") as f:
                f.write(molecule.get_xyz_string())
                f.write("\n")

            # Implement constraint optimization

            state_specific_molecules = []

            if self.identify_relevant_int_coordinates[0]:
                for state_to_optim in addition_of_state_specific_points:
                    drivers = self.drivers['gs']
                    current_basis = MolecularBasis.read(molecule, self.basis_set_label['gs'])

                    self.impes_drivers[state_to_optim].compute(molecule)
                    current_weights = self.impes_drivers[state_to_optim].weights
                    weights = [value for _, value in current_weights.items()]
                    used_labels = [label_idx for label_idx, _ in current_weights.items()]
                    # Sort labels and weights by descending weight
                    sorted_items_org = sorted(zip(used_labels, weights), key=lambda x: x[1], reverse=True)
                    total_weight = sum(weights)
                    cumulative_weight = 0.0
                    internal_coordinate_datapoints = []
                    for label, weight in sorted_items_org:
                        cumulative_weight += weight
                        internal_coordinate_datapoints.append((self.qm_data_point_dict[state_to_optim][label], weight))
                        if cumulative_weight >= 0.8 * total_weight:
                            break
                    _, candidate_constraints, _, fallback_constraints = self.impes_drivers[state_to_optim].determine_important_internal_coordinates(state_specific_gradients[state_to_optim][0], molecule, self.root_z_matrix[state_to_optim], internal_coordinate_datapoints)
                    selected = list(candidate_constraints)
                    selected = [tuple(int(x) for x in c) for c in selected]
                    selected = list(dict.fromkeys(selected))  # preserve order, remove duplicates
                    main_constraint_list = selected

                    imp_coord_constraint = main_constraint_list.copy()
                    constraints = []
                    selected_atoms = {a for coord in selected for a in coord}

                    for dih in self.root_z_matrix[state_to_optim]['dihedrals']:
                        mid_bond = tuple(sorted((dih[1], dih[2])))

                        if mid_bond not in self.impes_drivers[state_to_optim].symmetry_information[5]:
                            continue

                        touches_fault = bool(selected_atoms.intersection(dih))
                        if touches_fault and dih not in main_constraint_list:
                            constraints.append(dih)

                    main_constraint_list.extend(constraints)
                    full_z = (
                        self.root_z_matrix[state_to_optim]["bonds"]
                        + self.root_z_matrix[state_to_optim]["angles"]
                        + self.root_z_matrix[state_to_optim]["dihedrals"]
                        + self.root_z_matrix[state_to_optim]["impropers"]
                    )

                    coord_to_idx = {
                        tuple(int(x) for x in coord): idx
                        for idx, coord in enumerate(full_z)
                    }

                    b_matrix = self.impes_drivers[state_to_optim].impes_coordinate.b_matrix

                    if len(main_constraint_list) > 2:
                        constraint_mask = self._select_independent_constraints(
                            main_constraint_list, coord_to_idx, b_matrix)
                    else:
                        constraint_mask = main_constraint_list.copy()
                    self.ostream.print_blank()
                    self.ostream.print_info("Selected constraints for the optimization: ")
                    self.ostream.print_line(f"Filtered constraints: {constraint_mask}")
                    if self._user_constraints():
                        self.ostream.print_line(f"User constraints (1-based): {self._user_constraints()}")
                    self.ostream.flush()

                    opt_results, opt_gradient = self._optimize_with_constraints(
                        drivers, molecule, current_basis, scf_results, constraint_mask)

                    optimized_molecule = Molecule.from_xyz_string(opt_results['final_geometry'])
                    optimized_molecule.set_charge(molecule.get_charge())
                    optimized_molecule.set_multiplicity(molecule.get_multiplicity())
                    current_basis = MolecularBasis.read(optimized_molecule, current_basis.get_main_basis_label())
                    # determine distance of the new structure:
                    self.impes_drivers[state_to_optim].compute(optimized_molecule)
                    trial_weights = self.impes_drivers[state_to_optim].weights
                    weights = [value for _, value in trial_weights.items()]
                    used_labels = [label_idx for label_idx, _ in trial_weights.items()]
                    # Sort labels and weights by descending weight
                    sorted_items = sorted(zip(used_labels, weights), key=lambda x: x[1], reverse=True)
                    delta_e = abs(
                        opt_results['opt_energies'][-1] - self.impes_drivers[state_to_optim].get_energy()
                    ) * hartree_in_kcalpermol()

                    current_gradient_difference = (opt_gradient - self.impes_drivers[state_to_optim].impes_coordinate.gradient) * hartree_in_kcalpermol() / bohr_in_angstrom()
                    current_rmsd_gradient = np.sqrt((current_gradient_difference**2).mean())

                    self._print_impes_diagnostics(state_to_optim, sorted_items, delta_e,
                                                  current_rmsd_gradient)
                    needs_locality_fallback = (bool(fallback_constraints) and abs(delta_e) < abs(state_specific_energies[state_to_optim][0] - state_specific_energies[state_to_optim][1]) * hartree_in_kcalpermol() * 0.5 and
                                               bool(fallback_constraints) and current_rmsd_gradient < current_state_difference[state_to_optim][1] * 0.5)

                    if needs_locality_fallback and self.consider_locality:
                        self.ostream.print_blank()
                        self.ostream.print_info("Locality fallback triggered! --> Starting to include additional constraints based on the connectivity to the initially identified important coordinates.")
                        self.ostream.flush()
                        converged = False
                        correct_order = (bool(sorted_items_org[0][0] == sorted_items[0][0]) and abs(delta_e) > abs(state_specific_energies[state_to_optim][0] - state_specific_energies[state_to_optim][1]) * hartree_in_kcalpermol() * 0.3)
                        if correct_order:
                            converged = True
                        increasing_DOFs = 0
                        while not converged:
                            increasing_DOFs += 3
                            if increasing_DOFs >= len(fallback_constraints):
                                break

                            guarded_constraints = list(dict.fromkeys(main_constraint_list + fallback_constraints[:increasing_DOFs]))
                            constraint_mask = self._select_independent_constraints(
                                guarded_constraints, coord_to_idx, b_matrix)

                            if isinstance(drivers[0], (ScfRestrictedDriver, ScfUnrestrictedDriver)):
                                current_basis = MolecularBasis.read(molecule, current_basis.get_main_basis_label())
                            opt_results, opt_gradient = self._optimize_with_constraints(
                                drivers, molecule, current_basis, scf_results, constraint_mask)

                            guarded_constraints = constraint_mask.copy()
                            optimized_molecule = Molecule.from_xyz_string(opt_results['final_geometry'])
                            optimized_molecule.set_charge(molecule.get_charge())
                            optimized_molecule.set_multiplicity(molecule.get_multiplicity())

                            self.impes_drivers[state_to_optim].compute(optimized_molecule)

                            current_weights = self.impes_drivers[state_to_optim].weights
                            weights = [value for _, value in current_weights.items()]
                            used_labels = [label_idx for label_idx, _ in current_weights.items()]
                            # Sort labels and weights by descending weight
                            sorted_items = sorted(zip(used_labels, weights), key=lambda x: x[1], reverse=True)
                            delta_e = abs(opt_results['opt_energies'][-1] - self.impes_drivers[state_to_optim].get_energy()) * hartree_in_kcalpermol()

                            current_gradient_difference = (opt_gradient - self.impes_drivers[state_to_optim].impes_coordinate.gradient) * hartree_in_kcalpermol() / bohr_in_angstrom()
                            current_rmsd_gradient = np.sqrt((current_gradient_difference**2).mean())

                            self._print_impes_diagnostics(
                                state_to_optim, sorted_items, delta_e, current_rmsd_gradient,
                                gradient_target=current_state_difference[state_to_optim][1] * 0.5)

                            prev_fall_const_size = len(fallback_constraints)
                            for key, coord_type_entries in self.root_z_matrix[state_to_optim].items():
                                for coord in coord_type_entries:
                                    if coord not in fallback_constraints:
                                        fallback_constraints.append(coord)
                                        break

                                if len(fallback_constraints) != prev_fall_const_size:
                                    break

                            needs_locality_fallback = (bool(fallback_constraints) and abs(delta_e) < abs(state_specific_energies[state_to_optim][0] - state_specific_energies[state_to_optim][1]) * hartree_in_kcalpermol() * 0.5 and
                                                       bool(fallback_constraints) and current_rmsd_gradient < current_state_difference[state_to_optim][1] * 0.5)

                            if needs_locality_fallback:
                                correct_order = (bool(sorted_items_org[0][0] == sorted_items[0][0]) and abs(delta_e) < abs(state_specific_energies[state_to_optim][0] - state_specific_energies[state_to_optim][1]) * hartree_in_kcalpermol() * 0.3 and
                                                 bool(sorted_items_org[0][0] == sorted_items[0][0]) and current_rmsd_gradient < current_state_difference[state_to_optim][1] * 0.3)
                                if correct_order:
                                    converged = True
                                    main_constraint_list = guarded_constraints.copy()
                                else:
                                    optimized_molecule = molecule
                            else:
                                main_constraint_list = guarded_constraints.copy()
                                converged = True

                    imp_int_coord = {'bonds': [], 'angles': [], 'dihedrals': [], 'impropers': []}
                    for element in imp_coord_constraint:

                        if element in self.root_z_matrix[state_to_optim]['bonds']:
                            imp_int_coord['bonds'].append(element)
                        elif element in self.root_z_matrix[state_to_optim]['angles']:
                            imp_int_coord['angles'].append(element)
                        elif element in self.root_z_matrix[state_to_optim]['dihedrals']:
                            imp_int_coord['dihedrals'].append(element)
                        elif element in self.root_z_matrix[state_to_optim]['impropers']:
                            imp_int_coord['impropers'].append(element)

                    state_specific_molecules.append((optimized_molecule, current_basis, [state_to_optim], imp_int_coord))

                self.add_point(state_specific_molecules)
                self.last_point_added = self.point_checker - 1
                self.point_checker = 0
            else:

                current_basis = MolecularBasis.read(molecule, self.basis_set_label['gs'])
                state_specific_molecules.append((molecule, current_basis, addition_of_state_specific_points, {'bonds': [], 'angles': [], 'dihedrals': [], 'impropers': []}))

                self.add_point(state_specific_molecules)
                self.last_point_added = self.point_checker - 1
                self.point_checker = 0

    def _select_independent_constraints(self, candidates, coord_to_idx, b_matrix):
        """
        Keeps the candidate constraints whose B-matrix rows are linearly
        independent of the rows already kept (independence eta > 0.3); the
        first usable candidate is always kept.

        :param candidates:
            The candidate internal coordinates (0-based index tuples).
        :param coord_to_idx:
            Maps a z-matrix coordinate to its row in the B-matrix.
        :param b_matrix:
            The Wilson B-matrix at the current structure.

        :return:
            The selected coordinates as tuples of int.
        """

        constraint_mask = []
        active_rows = []
        for pos, coord in enumerate(candidates):
            coord = tuple(int(x) for x in coord)
            if coord not in coord_to_idx:
                continue
            row = b_matrix[coord_to_idx[coord], :]
            if np.linalg.norm(row) <= 1.0e-12:
                continue
            # Always keep the primary constraint.
            if pos == 0:
                constraint_mask.append(coord)
                active_rows.append(row)
                continue
            B_active = np.asarray(active_rows, dtype=np.float64)
            eta, _, _ = self.candidate_independence_eta(B_active, row)
            if eta > 0.3:
                constraint_mask.append(coord)
                active_rows.append(row)
        return constraint_mask

    def _print_impes_diagnostics(self, state, sorted_items, delta_e, rmsd_gradient,
                                 gradient_target=None):
        """
        Prints the interpolation error at a re-optimised structure.

        :param state:
            The electronic state.
        :param sorted_items:
            (datapoint index, weight) pairs sorted by descending weight.
        :param delta_e:
            The energy difference in kcal/mol.
        :param rmsd_gradient:
            The gradient difference RMSD in kcal/mol/angstrom.
        :param gradient_target:
            The gradient RMSD the locality fallback aims for (printed if given).
        """

        self.ostream.print_blank()
        self.ostream.print_info(f"IM-PES diagnostics for state {state}")
        self.ostream.print_block("=" * 80)
        self.ostream.print_blank()
        self.ostream.print_block("Weights:")
        for label, weight in sorted_items[:10]:
            self.ostream.print_block(f"  datapoint {label}: {weight: .10e}")

        self.ostream.print_blank()
        self.ostream.print_block("Energy difference:")
        self.ostream.print_block(f"  ΔE = {delta_e: .10f} kcal/mol")

        self.ostream.print_block("=" * 80 + "\n")

        self.ostream.print_block("Gradient difference RMSD: ")
        self.ostream.print_block(f"  {rmsd_gradient} kcal/mol/angstrom")
        if gradient_target is not None:
            self.ostream.print_block(f"  {gradient_target} kcal/mol/angstrom")
        self.ostream.flush()

    def normalize_constraint_rows(self, Bc, eps=1e-14):
        """
        Normalizes the rows of a constraint Jacobian for conditioning diagnostics.

        :param Bc:
            The constraint rows, shape (m, 3N).
        :param eps:
            Rows with a norm at or below eps are inactive.

        :return:
            The row-normalized Jacobian, the norms of the original rows and the mask
            of the active rows (norm > eps).
        """
        Bc = np.asarray(Bc, dtype=np.float64)

        if Bc.size == 0:
            return Bc.copy(), np.zeros(0), np.zeros(0, dtype=bool)

        row_norms = np.linalg.norm(Bc, axis=1)
        active_mask = row_norms > eps

        Bn = np.zeros_like(Bc)
        Bn[active_mask] = Bc[active_mask] / row_norms[active_mask, None]

        return Bn, row_norms, active_mask

    def candidate_independence_eta(
        self,
        B_active,
        b_candidate,
        pinv_rcond=1e-10,
        eps=1e-14,
    ):
        """
        Computes the incremental independence score eta of a candidate constraint
        row: the fraction of the normalized row that is not explained by the span
        of the already accepted rows (eta ~ 1: independent, eta ~ 0: redundant).

        :param B_active:
            The accepted constraint rows, shape (m, 3N).
        :param b_candidate:
            The candidate row, shape (3N,).
        :param pinv_rcond:
            The rcond of the pseudo-inverse.
        :param eps:
            Norms below eps count as zero.

        :return:
            eta, the norm of the residual and the norm of the candidate row.
        """
        b = np.asarray(b_candidate, dtype=np.float64).reshape(1, -1)
        B = np.asarray(B_active, dtype=np.float64)

        # Normalize candidate.
        b_norm = np.linalg.norm(b)
        if b_norm < eps:
            return 0.0, 0.0, float(b_norm)

        bn = b / b_norm

        # No active constraints -> fully independent.
        if B.size == 0 or B.shape[0] == 0:
            return 1.0, 1.0, float(b_norm)

        Bn, row_norms, active_mask = self.normalize_constraint_rows(B, eps=eps)
        Bn = Bn[active_mask]

        if Bn.shape[0] == 0:
            return 1.0, 1.0, float(b_norm)

        # Project candidate row onto row space of active B.
        # Projection onto row space: b_parallel = b @ B^+ @ B
        B_pinv = np.linalg.pinv(Bn, rcond=pinv_rcond)
        b_parallel = bn @ B_pinv @ Bn
        b_perp = bn - b_parallel

        residual_norm = float(np.linalg.norm(b_perp))
        eta = residual_norm / max(float(np.linalg.norm(bn)), eps)

        return float(eta), residual_norm, float(b_norm)

    def _user_constraints(self):
        """
        Returns the user's global constraints (use_minimized_structures[1] of the
        generator): 1-based geomeTRIC index tuples or geomeTRIC strings.

        :return:
            The list of user constraints (empty if none were given).
        """

        setting = self.identify_relevant_int_coordinates
        if isinstance(setting, (list, tuple)) and len(setting) > 1 and setting[1]:
            return list(setting[1])
        return []

    def _with_user_constraints(self, constraints):
        """
        Appends the user's global constraints to the automatically selected
        constraints of the re-optimisation of a new datapoint.

        The automatic constraints are 0-based z-matrix tuples that
        _build_opt_constraint_list shifts by index_offset=1; user tuples are
        1-based, so they are shifted to 0-based here. Strings are passed on
        unchanged. Coordinates that are already constrained are not repeated.

        :param constraints:
            The automatically selected 0-based constraints.

        :return:
            A new list with the user constraints appended.
        """

        combined = list(constraints)
        present = {tuple(int(x) for x in coord) for coord in combined
                   if not isinstance(coord, str)}
        for user_constraint in self._user_constraints():
            if isinstance(user_constraint, str):
                if user_constraint not in combined:
                    combined.append(user_constraint)
                continue
            zero_based = tuple(int(x) - 1 for x in user_constraint)
            if zero_based in present or zero_based[::-1] in present:
                continue
            combined.append(zero_based)
            present.add(zero_based)
        return combined

    def _optimize_with_constraints(self, drivers, molecule, basis, scf_results, constraints):
        """
        Runs the constrained re-optimisation of a structure that is added to
        the database.

        :param drivers:
            The (energy, gradient, Hessian) drivers of the ground state.
        :param molecule:
            The structure from the trajectory.
        :param basis:
            The basis set (SCF drivers only).
        :param scf_results:
            The SCF results of the QM check at this structure (SCF drivers only);
            they are reused instead of running the SCF again.
        :param constraints:
            The automatically selected 0-based constraints; the user's global
            constraints are added.

        :return:
            The optimization results and the final gradient.
        """

        constraints = self._with_user_constraints(constraints)
        if isinstance(drivers[0], (ScfRestrictedDriver, ScfUnrestrictedDriver)):
            return self._run_optimization(drivers[0], molecule, constraints=constraints,
                                          index_offset=1, compute_args=(basis, scf_results))
        if isinstance(drivers[0], XtbDriver):
            return self._run_optimization(drivers[0], molecule, constraints=constraints,
                                          index_offset=1)
        assert_msg_critical(
            False,
            'IMDatabasePointCollecter: unsupported QM driver '
            f'{type(drivers[0]).__name__} for the constrained optimization of '
            'new datapoints; use ScfRestrictedDriver, ScfUnrestrictedDriver or XtbDriver.')

    def _build_opt_constraint_list(self, constraints, index_offset=1):
        """
        Converts constraints to geomeTRIC 'freeze' strings.

        :param constraints:
            Index tuples (2: distance, 3: angle, 4: dihedral) or geomeTRIC strings;
            strings are passed on unchanged.
        :param index_offset:
            Added to every index (1 for 0-based tuples, 0 for 1-based tuples).

        :return:
            The list of geomeTRIC constraint strings.
        """

        opt_constraint_list = []
        for constraint in constraints:
            if isinstance(constraint, str):
                opt_constraint_list.append(constraint)
                continue

            shifted = [value + index_offset for value in constraint]
            if len(shifted) == 2:
                opt_constraint = f"freeze distance {shifted[0]} {shifted[1]}"
            elif len(shifted) == 3:
                opt_constraint = f"freeze angle {shifted[0]} {shifted[1]} {shifted[2]}"
            else:
                opt_constraint = f"freeze dihedral {shifted[0]} {shifted[1]} {shifted[2]} {shifted[3]}"
            opt_constraint_list.append(opt_constraint)

        return opt_constraint_list

    def _run_optimization(self, optimization_driver, molecule, constraints=None, transition=False, index_offset=1, compute_args=None):
        """
        Runs a (constrained) geometry optimization with geomeTRIC.

        :param optimization_driver:
            The energy driver (SCF or XtbDriver).
        :param molecule:
            The start structure.
        :param constraints:
            The constraints (see _build_opt_constraint_list), or None.
        :param transition:
            True for a transition-state search.
        :param index_offset:
            The index offset of tuple constraints.
        :param compute_args:
            Extra arguments of the gradient driver (basis, SCF results) or None.

        :return:
            The optimization results and the gradient at the final structure.
        """

        opt_drv = OptimizationDriver(optimization_driver)
        opt_drv.ostream.mute()
        opt_drv.transition = transition

        if constraints is not None:
            opt_drv.constraints = self._build_opt_constraint_list(constraints, index_offset=index_offset)

        if compute_args is None:
            opt_results = opt_drv.compute(molecule)
        else:
            opt_results = opt_drv.compute(molecule, *compute_args)

        return opt_results, opt_drv.grad_drv.gradient

    def add_point(self, state_specific_molecules):
        """
        Computes the QM energy, gradient and Hessian at the given structures and
        adds them as new datapoints to the database of root 0 (in memory and in
        the HDF5 database file), writes the pending reference structures and saves
        the OpenMM state to checkpoint.xml.

        :param state_specific_molecules:
            List of tuples (molecule, basis, states, important internal
            coordinates); states is ignored (ground state only), the important
            coordinates are a dict {'bonds', 'angles', 'dihedrals', 'impropers'}
            of 0-based tuples in z-matrix orientation.
        """

        self.ostream.print_blank()
        self.ostream.print_header("Adding a new point to the database! --> Starting to compute the necessary information for interpolation.")
        self.ostream.print_blank()
        self.ostream.flush()

        # ground state only (checked at the start of run_qmmm): every new
        # point is computed with the ground-state drivers and stored for root 0
        root = 0
        drivers = self.drivers['gs']

        for molecule, basis, _, imp_int_coordinates in state_specific_molecules:

            energies, scf_results, rsp_results = self._compute_energy(drivers[0], molecule, basis)

            gradients = self._compute_gradient(drivers[1], molecule, basis, scf_results, rsp_results)
            hessians = self._compute_hessian(drivers[2], molecule, basis)

            inv_sqrt_masses = None
            if self.use_mass_weight:
                masses = molecule.get_masses().copy()
                masses_cart = np.repeat(masses, 3)
                inv_sqrt_masses = 1.0 / np.sqrt(masses_cart)

            new_label = f"point_{len(self.qm_data_point_dict[root]) + 1}"

            grad_vec = gradients[0].reshape(-1)         # (3N,)
            hess_mat = hessians[0].reshape(grad_vec.size, grad_vec.size)
            mw_grad_vec = grad_vec
            mw_hess_mat = hess_mat
            if self.use_mass_weight:
                mw_grad_vec = inv_sqrt_masses * grad_vec
                mw_hess_mat = (inv_sqrt_masses[:, None] * hess_mat) * inv_sqrt_masses[None, :]

            impes_coordinate = InterpolationDatapoint(self.root_z_matrix[root])
            impes_coordinate.update_settings(self.interpolation_settings[root])
            impes_coordinate.cartesian_coordinates = molecule.get_coordinates_in_bohr()
            impes_coordinate.imp_int_coordinates = imp_int_coordinates
            impes_coordinate.eq_bond_lengths = self.qm_data_point_dict[root][0].eq_bond_lengths
            impes_coordinate.inv_sqrt_masses = inv_sqrt_masses
            impes_coordinate.energy = energies[0]
            impes_coordinate.gradient = mw_grad_vec.reshape(gradients[0].shape)
            impes_coordinate.hessian = mw_hess_mat.reshape(hessians[0].shape)
            impes_coordinate.transform_gradient_and_hessian()

            impes_coordinate.confidence_radius = self.use_opt_confidence_radius[2]

            category_label = new_label
            if impes_coordinate.use_inverse_bond_length:
                category_label += "_rinv"
            elif impes_coordinate.use_eq_bond_length:
                category_label += "_eq"
            else:
                category_label += "_r"

            category_label += "_dihedral"

            self.qm_data_point_dict[root].append(impes_coordinate)
            self.sorted_state_spec_im_labels[root].append(new_label)

            impes_coordinate.point_label = category_label

            self.impes_drivers[root].qm_data_points = self.qm_data_point_dict[root]
            self.impes_drivers[root].labels = self.sorted_state_spec_im_labels[root]
            self._refresh_interpolation_driver_caches(root)

            self.density_around_data_point[root] += 1

            impes_coordinate.write_hdf5(
                self.interpolation_settings[root]['imforcefield_file'],
                new_label
            )

            self.ostream.print_blank()
            self.ostream.print_header(f"Database expansion: Added point {new_label} to the database of root {root} with energy {energies[0]} Hartree.")
            self.ostream.flush()

            if self.sampling_enabled:
                self._write_sampling_point_from_geometry(
                    root=root,
                    molecule=molecule,
                    label=new_label,
                    template_point=impes_coordinate
                )

            self.write_qm_energy_determined_points_h5(
                self.allowed_molecules[root]['molecules'][self.last_added:],
                self.allowed_molecules[root]['qm_energies'][self.last_added:],
                self.allowed_molecules[root]['qm_gradients'][self.last_added:],
                self.allowed_molecules[root]['im_energies'][self.last_added:],
                root)

            self.last_added = len(self.allowed_molecules[root]['molecules'])

        self.simulation.saveState('checkpoint.xml')

    def _write_sampling_point_from_geometry(self, root, molecule, label, template_point):
        """
        Adds the point of a new datapoint to the sampling mirror database:
        computes xTB energy, gradient and Hessian at the geometry, mass-weights
        them if sampling_interpolation_settings[root]['use_mass_weight'] is set,
        writes the point under the same label and reloads the mirror database.

        :param root:
            The root.
        :param molecule:
            The structure of the new datapoint.
        :param label:
            The label of the new datapoint.
        :param template_point:
            The new datapoint of the main database (geometry, important
            coordinates, confidence radius).
        """

        sampling_qm, sampling_grad, sampling_hess = self.sampling_driver['gs']
        e, _, _ = self._compute_energy(sampling_qm, molecule, basis=None)
        g = self._compute_gradient(
            sampling_grad, molecule, basis=None, scf_results=None, rsp_results=None
        )
        h = self._compute_hessian(sampling_hess, molecule, basis=None)

        inv_sqrt = None
        grad_vec = g[0].reshape(-1)
        hess_mat = h[0].reshape(grad_vec.size, grad_vec.size)
        mw_grad = grad_vec
        mw_hess = hess_mat
        if self.sampling_interpolation_settings[root].get('use_mass_weight', False):
            masses = molecule.get_masses().copy()
            inv_sqrt = 1.0 / np.sqrt(np.repeat(masses, 3))
            mw_grad = inv_sqrt * grad_vec
            mw_hess = (inv_sqrt[:, None] * hess_mat) * inv_sqrt[None, :]

        dp = InterpolationDatapoint(self.root_z_matrix[root])
        dp.update_settings(self.sampling_interpolation_settings[root])

        # Geometry + transformed tensors
        dp.cartesian_coordinates = np.array(template_point.cartesian_coordinates, copy=True)
        dp.eq_bond_lengths = None if template_point.eq_bond_lengths is None else np.array(template_point.eq_bond_lengths, copy=True)
        dp.imp_int_coordinates = copy.deepcopy(
            getattr(
                template_point,
                "imp_int_coordinates",
                {"bonds": [], "angles": [], "dihedrals": [], "impropers": []},
            )
        )
        dp.inv_sqrt_masses = inv_sqrt
        dp.energy = e[0]
        dp.gradient = mw_grad.reshape(g[0].shape)
        dp.hessian = mw_hess.reshape(h[0].shape)
        dp.confidence_radius = template_point.confidence_radius

        dp.point_label = label

        dp.transform_gradient_and_hessian()

        dp.write_hdf5(self.sampling_interpolation_settings[root]['imforcefield_file'], label)

        # Refresh in-memory sampling interpolator
        self._reload_sampling_root_from_hdf5(root, inv_sqrt)

    def _reload_sampling_root_from_hdf5(self, root, inv_sqrt_masses):
        """
        Reloads all points of the sampling mirror database of one root in label
        order.

        :param root:
            The root.
        :param inv_sqrt_masses:
            The inverse square roots of the Cartesian masses, or None.
        """

        driver_object = self.sampling_impes_drivers[root]

        im_labels, _ = driver_object.read_labels()
        im_labels = sorted(im_labels, key=lambda x: int(x.split('_')[1]))
        self.sampling_qm_data_point_dict[root] = []

        for label in im_labels:
            dp = InterpolationDatapoint(self.root_z_matrix[root])
            dp.update_settings(self.sampling_interpolation_settings[root])
            dp.read_hdf5(self.sampling_interpolation_settings[root]['imforcefield_file'], label)
            dp.inv_sqrt_masses = inv_sqrt_masses

            self.sampling_qm_data_point_dict[root].append(dp)

        driver_object.qm_data_points = self.sampling_qm_data_point_dict[root]
        driver_object.labels = self.sorted_state_spec_im_labels[root]

        if len(self.sampling_qm_data_point_dict[root]) > 0:
            driver_object.impes_coordinate.eq_bond_lengths = self.sampling_qm_data_point_dict[root][0].eq_bond_lengths

        driver_object.mark_runtime_data_cache_dirty()

    def _clone_datapoints_with_alphas(self, datapoints, alphas):
        """
        Returns shallow copies of datapoints with the given confidence radii.

        :param datapoints:
            The datapoints.
        :param alphas:
            One confidence radius per datapoint.

        :return:
            The list of copies.
        """

        alpha_vec = np.asarray(alphas, dtype=np.float64).reshape(-1)
        if len(datapoints) != alpha_vec.size:
            raise ValueError(
                f"Alpha/datapoint size mismatch: n_dp={len(datapoints)} vs n_alpha={alpha_vec.size}"
            )

        cloned = []
        for dp, alpha in zip(datapoints, alpha_vec):
            dp_clone = copy.copy(dp)
            dp_clone.confidence_radius = float(alpha)
            cloned.append(dp_clone)

        return cloned

    def _build_interp_driver_for_trust_radius_eval(
        self,
        *,
        z_matrix,
        interpolation_setting,
        sym_dict,
        exponent_p_q,
        datapoints_eval,
    ):
        """
        Builds an interpolation driver on a given set of datapoints, for the
        evaluation of the confidence-radius fit.

        :param z_matrix:
            The z-matrix.
        :param interpolation_setting:
            The interpolation settings.
        :param sym_dict:
            The symmetry information.
        :param exponent_p_q:
            The weight exponents (p, q).
        :param datapoints_eval:
            The datapoints.

        :return:
            The interpolation driver.
        """

        drv = InterpolationDriver(z_matrix, ostream=self.ostream)
        drv.update_settings(interpolation_setting)
        drv.symmetry_information = sym_dict
        drv.calc_optim_trust_radius = False
        drv.exponent_p = exponent_p_q[0]
        drv.exponent_q = exponent_p_q[1]

        drv.impes_coordinate.eq_bond_lengths = datapoints_eval[0].eq_bond_lengths

        if len(datapoints_eval) > 0 and hasattr(datapoints_eval[0], "inv_sqrt_masses"):
            drv.impes_coordinate.inv_sqrt_masses = datapoints_eval[0].inv_sqrt_masses
        elif self.inv_sqrt_masses is not None:
            drv.impes_coordinate.inv_sqrt_masses = self.inv_sqrt_masses

        drv.qm_data_points = datapoints_eval

        return drv

    def _evaluate_interpolation_errors_on_reference_set(
        self,
        *,
        scenario_name,
        molecules,
        qm_energies,
        qm_gradients,
        datapoints_eval,
        z_matrix,
        interpolation_setting,
        sym_dict,
        exponent_p_q,
    ):
        """
        Evaluates the interpolation errors of a set of datapoints at the reference
        structures.

        :param scenario_name:
            The name of the scenario (for the printout).
        :param molecules:
            The reference structures.
        :param qm_energies:
            The QM energies (Hartree).
        :param qm_gradients:
            The QM gradients (Hartree/bohr).
        :param datapoints_eval:
            The datapoints.
        :param z_matrix:
            The z-matrix.
        :param interpolation_setting:
            The interpolation settings.
        :param sym_dict:
            The symmetry information.
        :param exponent_p_q:
            The weight exponents (p, q).

        :return:
            A dictionary with the energy errors (kcal/mol), the gradient RMSDs
            (kcal/mol/angstrom), the weights of every datapoint at every reference
            and summary statistics.
        """

        drv = self._build_interp_driver_for_trust_radius_eval(
            z_matrix=z_matrix,
            interpolation_setting=interpolation_setting,
            sym_dict=sym_dict,
            exponent_p_q=exponent_p_q,
            datapoints_eval=datapoints_eval,
        )

        n_dp = len(drv.qm_data_points)
        n_ref = len(molecules)
        w_dpref = np.zeros((n_dp, n_ref), dtype=np.float64)

        im_energies = []
        im_gradients = []

        for i, mol in enumerate(molecules):
            drv.compute(mol)

            for dp_idx, weight in drv.weights.items():
                w_dpref[int(dp_idx), i] = float(weight)

            im_energies.append(float(drv.get_energy()))
            im_gradients.append(np.asarray(drv.get_gradient(), dtype=np.float64))

        return self._interpolation_error_summary(
            scenario_name=scenario_name,
            n_dp=len(datapoints_eval),
            im_energies=im_energies,
            im_gradients=im_gradients,
            qm_energies=qm_energies,
            qm_gradients=qm_gradients,
            weights_refs=w_dpref,
        )

    def _interpolation_error_summary(self, *, scenario_name, n_dp, im_energies, im_gradients,
                                     qm_energies, qm_gradients, weights_refs):
        """
        Returns the result dictionary of
        _evaluate_interpolation_errors_on_reference_set from the interpolated
        energies (Hartree) and gradients (Hartree/bohr) at the references.

        :param scenario_name:
            The name of the scenario (for the printout).
        :param n_dp:
            The number of datapoints.
        :param im_energies:
            The interpolated energies.
        :param im_gradients:
            The interpolated gradients.
        :param qm_energies:
            The QM energies (Hartree).
        :param qm_gradients:
            The QM gradients (Hartree/bohr).
        :param weights_refs:
            The weights of every datapoint at every reference, (n_dp, n_ref),
            or None.

        :return:
            The result dictionary.
        """

        e_conv = hartree_in_kcalpermol()
        g_conv = hartree_in_kcalpermol() / bohr_in_angstrom()

        e_err = []
        g_rmsd = []

        for i in range(len(im_energies)):
            e_im = float(im_energies[i])
            e_qm = float(np.asarray(qm_energies[i], dtype=np.float64).reshape(-1)[0])
            g_qm = np.asarray(qm_gradients[i], dtype=np.float64)
            g_im = np.asarray(im_gradients[i], dtype=np.float64).reshape(g_qm.shape)

            de = (e_im - e_qm) * e_conv
            dg = (g_im - g_qm) * g_conv
            rmsd = float(np.sqrt(np.mean(dg**2)))

            e_err.append(float(de))
            g_rmsd.append(rmsd)

        e_arr = np.asarray(e_err, dtype=np.float64)
        g_arr = np.asarray(g_rmsd, dtype=np.float64)

        return {
            "scenario": scenario_name,
            "n_dp": int(n_dp),
            "n_ref": int(len(im_energies)),
            "energy_errors_kcal": e_arr.tolist(),
            "gradient_rmsd_kcal_per_mol_A": g_arr.tolist(),
            "weights_refs": weights_refs,
            "summary": {
                "energy_mae": float(np.mean(np.abs(e_arr))),
                "energy_rmse": float(np.sqrt(np.mean(e_arr**2))),
                "energy_max_abs": float(np.max(np.abs(e_arr))),
                "gradient_mean_rmsd": float(np.mean(g_arr)),
                "gradient_rmse_rmsd": float(np.sqrt(np.mean(g_arr**2))),
                "gradient_max_rmsd": float(np.max(g_arr)),
            },
        }

    def _trust_radius_cache_unsupported(self, drv, sym_dict, natoms):
        """
        Returns why the cached radius fit (trust_radius_fit_cache) cannot
        reproduce the interpolation of the driver drv, or None if it can.

        :param drv:
            The interpolation driver of the radius fit.
        :param sym_dict:
            The symmetry information.
        :param natoms:
            The number of atoms.
        """

        if drv.interpolation_type != 'shepard' or drv.weightfunction_type != 'cartesian':
            return 'it needs Shepard interpolation with Cartesian weights'
        if drv.use_tc_weights and drv.tc_weight_mode != 'multiplicative':
            return 'it needs plain or multiplicative tc weights'
        active_atoms = np.delete(np.arange(natoms), sym_dict[4])
        if not np.array_equal(np.asarray(sym_dict[3], dtype=int), active_atoms):
            return 'symmetry_information[3] is not the complement of symmetry_information[4]'
        return None

    @staticmethod
    def _trust_radius_driver_key(drv, sym_dict, z_matrix):
        """
        Returns the settings of drv that the distances and gate metrics depend
        on, as a string.
        """

        return repr((
            drv.interpolation_type, drv.weightfunction_type, drv.exponent_p, drv.exponent_q,
            drv.use_tc_weights, drv.tc_weight_mode, drv.tc_imp_gate_lambda, drv.tc_imp_gate_exp_clip,
            drv.tc_imp_bond_sigma_angstrom, drv.tc_imp_angle_sigma_degrees,
            drv.tc_imp_dihedral_sigma_degrees, drv.tc_imp_improper_sigma_degrees, drv.tc_imp_min_sigma,
            drv.use_inverse_bond_length, drv.use_eq_bond_length,
            [int(a) for a in sym_dict[3]], [int(a) for a in sym_dict[4]], z_matrix,
        ))

    @staticmethod
    def _available_memory_bytes():
        """
        Returns MemAvailable of /proc/meminfo in bytes, or None.
        """

        try:
            with open('/proc/meminfo') as meminfo:
                for line in meminfo:
                    if line.startswith('MemAvailable:'):
                        return int(line.split()[1]) * 1024
        except (OSError, ValueError, IndexError):
            pass
        return None

    def _trust_radius_grouping_weights(self, molecules, datapoints, interpolation_setting,
                                       sym_dict, z_matrix, exponent_p_q):
        """
        Returns the normalised weights of every datapoint at every reference
        with all radii at 0.5, (n_dp, n_ref): the grouping input of
        determine_trust_radius_gradient, equal to the weights of its
        'with_new_datapoint_pre_opt' evaluation.

        The raw weights at radius 0.5 depend on the reference and the
        datapoint only. They are kept in _trust_radius_grouping_cache and
        computed only for new references and new datapoints (the lists only
        grow during a run). The raw weights come from
        InterpolationDriver.cartesian_distance and are normalised as in
        shepard_interpolation, so the weights are the same numbers.

        :return:
            The normalised weights.
        """

        n_ref = len(molecules)
        n_dp = len(datapoints)

        drv = self._build_interp_driver_for_trust_radius_eval(
            z_matrix=z_matrix,
            interpolation_setting=interpolation_setting,
            sym_dict=sym_dict,
            exponent_p_q=exponent_p_q,
            datapoints_eval=self._clone_datapoints_with_alphas(datapoints, [0.5] * n_dp),
        )
        clones = drv.qm_data_points
        key = self._trust_radius_driver_key(drv, sym_dict, z_matrix)

        cache = self._trust_radius_grouping_cache
        if (cache is None or cache['key'] != key
                or len(cache['molecules']) > n_ref or len(cache['datapoints']) > n_dp
                or any(a is not b for a, b in zip(cache['molecules'], molecules))
                or any(a is not b for a, b in zip(cache['datapoints'], datapoints))):
            cache = {'key': key, 'molecules': [], 'datapoints': [],
                     'raw_weights': np.zeros((0, 0)), 'distances': np.zeros((0, 0))}

        n_old_ref = len(cache['molecules'])
        n_old_dp = len(cache['datapoints'])
        raw_weights = np.zeros((n_ref, n_dp))
        distances = np.zeros((n_ref, n_dp))
        raw_weights[:n_old_ref, :n_old_dp] = cache['raw_weights']
        distances[:n_old_ref, :n_old_dp] = cache['distances']

        for s, mol in enumerate(molecules):
            first_new = n_old_dp if s < n_old_ref else 0
            if first_new == n_dp:
                continue
            drv.define_impes_coordinate(mol.get_coordinates_in_bohr())
            for j in range(first_new, n_dp):
                drv._check_datapoint_z_matrix(clones[j])
                distance, _, denominator, _, _, _, _, _ = drv.cartesian_distance(clones[j], store_alpha_gradients=False)
                distances[s, j] = distance
                raw_weights[s, j] = 1.0 / (denominator)

        self._trust_radius_grouping_cache = {
            'key': key, 'molecules': list(molecules), 'datapoints': list(datapoints),
            'raw_weights': raw_weights, 'distances': distances,
        }

        # normalisation of InterpolationDriver.shepard_interpolation
        weights = np.zeros((n_dp, n_ref), dtype=np.float64)
        eps = float(getattr(drv, "weight_eps", 1.0e-12))
        for s in range(n_ref):
            w_i = np.array(raw_weights[s], dtype=np.float64)
            S_raw = float(np.sum(w_i))
            finite_positive = np.isfinite(w_i) & (w_i > 0.0)
            if np.any(finite_positive):
                if np.isfinite(S_raw) and S_raw > eps:
                    w_norm = w_i
                else:
                    w_norm = np.where(finite_positive, w_i / float(np.max(w_i[finite_positive])), 0.0)
                weights[:, s] = w_norm / float(np.sum(w_norm))
            else:
                metric = np.where(np.isfinite(distances[s]), np.abs(distances[s]), np.inf)
                weights[int(np.argmin(metric)), s] = 1.0

        return weights

    def _trust_radius_pair_terms(self, *, molecules, qm_energies, qm_gradients, datapoints,
                                 initial_alphas, free_idx, interpolation_setting, sym_dict,
                                 z_matrix, exponent_p_q):
        """
        Evaluates the interpolation once at every reference and keeps its
        radius-independent terms (IMTrustRadiusPairTerms) for the cached radius
        fit. The 'without_last_datapoint' and 'with_new_datapoint_pre_opt'
        errors of determine_trust_radius_gradient follow from the same terms.

        At every reference the closed form at the current radii is compared
        with the driver. A difference above 1e-8 Hartree in the energy or
        1e-7 Hartree/bohr in a gradient component stops the build.

        :param free_idx:
            The datapoints whose radii are fitted.

        :return:
            (pair terms, result without the last datapoint, result with all
            datapoints at radius 0.5, message); on a mismatch, (None, None,
            None, message).
        """

        drv = self._build_interp_driver_for_trust_radius_eval(
            z_matrix=z_matrix,
            interpolation_setting=interpolation_setting,
            sym_dict=sym_dict,
            exponent_p_q=exponent_p_q,
            datapoints_eval=self._clone_datapoints_with_alphas(datapoints, initial_alphas),
        )

        n_dp = len(datapoints)
        natoms = datapoints[0].cartesian_coordinates.shape[0]
        ndof = 3 * natoms
        active_atoms = np.delete(np.arange(natoms), sym_dict[4])
        active_dofs = np.array([3 * a + i for a in active_atoms for i in range(3)], dtype=int)

        pair_terms = IMTrustRadiusPairTerms(
            len(molecules), natoms, initial_alphas, free_idx, exponent_p_q,
            drv.tc_imp_gate_lambda, drv.tc_imp_gate_exp_clip)

        cols_wo_last = np.arange(n_dp - 1 if n_dp >= 2 else n_dp)
        alphas_wo_last = np.asarray(initial_alphas, dtype=np.float64)[cols_wo_last]
        half = np.full(n_dp, 0.5)

        energies_wo_last, gradients_wo_last = [], []
        energies_half, gradients_half = [], []

        for s, mol in enumerate(molecules):
            drv.pair_term_records = []
            try:
                drv.compute(mol)
                records = drv.pair_term_records
            finally:
                drv.pair_term_records = None

            if len(records) != n_dp:
                return None, None, None, f'{len(records)} distance records for {n_dp} datapoints'

            dist = np.array([rec[0] for rec in records], dtype=np.float64)
            V = np.zeros((n_dp, ndof))
            A = np.zeros(n_dp)
            GA = np.zeros((n_dp, ndof))
            for j, (_, v, a_imp, grad_a_imp) in enumerate(records):
                V[j, active_dofs] = v.reshape(-1)
                if a_imp is not None:
                    A[j] = a_imp
                    GA[j, active_dofs] = grad_a_imp.reshape(-1)
            T = np.asarray(drv.potentials, dtype=np.float64)
            GT = np.asarray(drv.gradients, dtype=np.float64).reshape(n_dp, ndof)

            e_cf, g_cf = pair_terms.interpolate_reference(dist, V, A, GA, T, GT, initial_alphas)
            e_diff = abs(e_cf - float(drv.get_energy()))
            g_diff = float(np.max(np.abs(g_cf - np.asarray(drv.get_gradient(), dtype=np.float64).reshape(-1))))
            if not (e_diff <= 1.0e-8 and g_diff <= 1.0e-7):
                return None, None, None, (f'closed form differs from the driver at reference {s} '
                                          f'(|dE| {e_diff:.2e} Hartree, max |dG| {g_diff:.2e} Hartree/bohr)')

            pair_terms.set_reference(s, dist, V, A, GA, T, GT)

            e_s, g_s = pair_terms.interpolate_reference(
                dist[cols_wo_last], V[cols_wo_last], A[cols_wo_last], GA[cols_wo_last],
                T[cols_wo_last], GT[cols_wo_last], alphas_wo_last)
            energies_wo_last.append(e_s)
            gradients_wo_last.append(g_s)

            e_s, g_s = pair_terms.interpolate_reference(dist, V, A, GA, T, GT, half)
            energies_half.append(e_s)
            gradients_half.append(g_s)

        result_wo_last = self._interpolation_error_summary(
            scenario_name="without_last_datapoint", n_dp=cols_wo_last.size,
            im_energies=energies_wo_last, im_gradients=gradients_wo_last,
            qm_energies=qm_energies, qm_gradients=qm_gradients, weights_refs=None)
        result_half = self._interpolation_error_summary(
            scenario_name="with_new_datapoint_pre_opt", n_dp=n_dp,
            im_energies=energies_half, im_gradients=gradients_half,
            qm_energies=qm_energies, qm_gradients=qm_gradients, weights_refs=None)

        return pair_terms, result_wo_last, result_half, None

    def _print_trust_radius_stage_comparison(self, stage_results):
        """
        Prints the interpolation errors before and after the confidence-radius fit.

        :param stage_results:
            The results of _evaluate_interpolation_errors_on_reference_set.
        """

        order = [
            "without_last_datapoint",
            "with_new_datapoint_pre_opt",
            "with_new_datapoint_post_opt",
        ]
        by_name = {row["scenario"]: row for row in stage_results}

        title = "Trust-Radius Interpolation Comparison"
        header = (
            f"{'Scenario':<32}"
            f"{'N_DP':>8}"
            f"{'E_MAE':>12}"
            f"{'E_RMSE':>12}"
            f"{'E_MAX':>12}"
            f"{'G_MEAN':>12}"
            f"{'G_RMSE':>12}"
            f"{'G_MAX':>12}"
        )

        self.ostream.print_blank()
        self.ostream.print_header(title)
        self.ostream.print_header("=" * len(title))
        self.ostream.print_blank()
        self.ostream.print_header("Interpolation comparison on fixed reference set")
        self.ostream.print_blank()
        self.ostream.print_header(header)
        self.ostream.print_header("-" * len(header))

        for key in order:
            row = by_name.get(key)
            if row is None:
                continue

            s = row["summary"]
            self.ostream.print_header(
                f"{key:<32}"
                f"{row['n_dp']:8d}"
                f"{s['energy_mae']:12.5f}"
                f"{s['energy_rmse']:12.5f}"
                f"{s['energy_max_abs']:12.5f}"
                f"{s['gradient_mean_rmsd']:12.5f}"
                f"{s['gradient_rmse_rmsd']:12.5f}"
                f"{s['gradient_max_rmsd']:12.5f}"
            )

        def _delta(a_name, b_name, field):
            a = by_name.get(a_name, {}).get("summary", {}).get(field, math.nan)
            b = by_name.get(b_name, {}).get("summary", {}).get(field, math.nan)
            return b - a

        delta_header = (
            f"{'Comparison':<32}"
            f"{'Delta E_MAE':>18}"
            f"{'Delta G_MEAN':>18}"
        )

        self.ostream.print_blank()
        self.ostream.print_header("Delta Summary")
        self.ostream.print_header("-------------")
        self.ostream.print_header(delta_header)
        self.ostream.print_header("-" * len(delta_header))

        self.ostream.print_header(
            f"{'pre_opt - without_last':<32}"
            f"{_delta('without_last_datapoint', 'with_new_datapoint_pre_opt', 'energy_mae'):18.6f}"
            f"{_delta('without_last_datapoint', 'with_new_datapoint_pre_opt', 'gradient_mean_rmsd'):18.6f}"
        )

        self.ostream.print_header(
            f"{'post_opt - pre_opt':<32}"
            f"{_delta('with_new_datapoint_pre_opt', 'with_new_datapoint_post_opt', 'energy_mae'):18.6f}"
            f"{_delta('with_new_datapoint_pre_opt', 'with_new_datapoint_post_opt', 'gradient_mean_rmsd'):18.6f}"
        )

        self.ostream.print_blank()
        self.ostream.flush()

    def determine_trust_radius_gradient(self, molecules, qm_energies, qm_gradients, im_energies, datapoints, interpolation_setting, sym_dict, z_matrix, exponent_p_q):
        """
        Fits the confidence radii of the datapoints to the reference structures
        (loss on energies and gradients, scipy L-BFGS-B with basin hopping for
        small reference sets), group by group of datapoints that share references.

        :param molecules:
            The reference structures.
        :param qm_energies:
            The QM energies of the references (Hartree).
        :param qm_gradients:
            The QM gradients of the references (Hartree/bohr).
        :param im_energies:
            The interpolated energies of the references; not used by the fit
            (kept for compatibility).
        :param datapoints:
            The datapoints whose radii are fitted.
        :param interpolation_setting:
            The interpolation settings.
        :param sym_dict:
            The symmetry information.
        :param z_matrix:
            The z-matrix.
        :param exponent_p_q:
            The weight exponents (p, q).

        :return:
            The fitted confidence radii, one per datapoint.
        """

        try:
            from scipy.optimize import minimize
            from scipy.optimize import basinhopping
            from scipy.sparse import csr_matrix
            from scipy.sparse.csgraph import connected_components
        except ImportError:
            error_msg = 'Scipy is required for determine_trust_radius_gradient.'
            assert_msg_critical(False, error_msg)

        def optimize_trust_radius(alpha_full_start, trainable_idx, geom_list, E_ref_list, G_ref_list, dps, impes_dict, sym_dict, exponent_p_q, basin_id=None, pair_terms=None, ref_rows=None):
            """
            Perform the gradient-based optimization to find R*
            that minimizes the sum of squared errors to reference QM energies.
            With pair_terms, the objective is evaluated in closed form on the
            rows ref_rows of pair_terms (IMTrustRadiusCachedOptimizer).
            """

            trainable_idx = np.asarray(trainable_idx, dtype=np.int64).reshape(-1)
            alpha_full_start = np.asarray(alpha_full_start, dtype=np.float64).reshape(-1)

            if trainable_idx.size == 0:
                raise ValueError("optimize_trust_radius called with empty trainable_idx")
            if np.any(trainable_idx < 0) or np.any(trainable_idx >= alpha_full_start.size):
                raise ValueError("trainable_idx out of range for alpha_full_start")

            # Reduced optimization variables (only active/group datapoints)
            x0 = alpha_full_start[trainable_idx].copy()

            bounds = [(1e-3, 1.0)] * len(trainable_idx)
            sample_size = len(geom_list)

            if sample_size < 100:
                n_bh_iter = 10
                gtol = 1e-3
                ftol = 1e-4
                maxiter = 300
            elif sample_size < 400:
                n_bh_iter = 1
                gtol = 1e-3
                ftol = 1e-4
                maxiter = 300
            elif sample_size <= 800:
                n_bh_iter = 1
                gtol = 1e-3
                ftol = 1e-4
                maxiter = 220
            else:
                n_bh_iter = 0
                gtol = 1e-3
                ftol = 1e-4
                maxiter = 140

            # Backward compatible interpretation:
            # [enabled, mode, init_confidence_radius, e_x]
            # Older input may provide only three entries; then fall back to idx 2.
            if isinstance(self.use_opt_confidence_radius, (list, tuple)) and len(self.use_opt_confidence_radius) > 3:
                energy_weight = float(self.use_opt_confidence_radius[3])
            else:
                energy_weight = 0.5

            if pair_terms is not None:
                opt = IMTrustRadiusCachedOptimizer(
                    pair_terms, ref_rows, E_ref_list, G_ref_list, sym_dict, trainable_idx,
                    e_x=energy_weight,
                    beta=0.8)
            else:
                opt = IMTrustRadiusOptimizer(
                    z_matrix, impes_dict, sym_dict, dps,
                    geom_list, E_ref_list, G_ref_list, exponent_p_q,
                    e_x=energy_weight,
                    beta=0.8)

            history = []
            progress_header_printed = {"value": False}
            iter_counter = {"value": 0}
            basin_label = "Init"
            if basin_id is not None:
                basin_label = str(int(basin_id))

            def _expand_alpha(x_var):
                x_var = np.asarray(x_var, dtype=np.float64).reshape(-1)
                if x_var.size != trainable_idx.size:
                    raise ValueError(
                        f"reduced alpha size mismatch: expected {trainable_idx.size}, got {x_var.size}"
                    )
                x_full = alpha_full_start.copy()
                x_full[trainable_idx] = x_var
                return x_full

            def _fmt_float(value):
                if value is None:
                    return "nan"
                try:
                    value = float(value)
                except Exception:
                    return "nan"
                if not np.isfinite(value):
                    return "nan"
                return f"{value:.6e}"

            def _append_history(tag, x_var, **extra):
                x_var_arr = np.asarray(x_var, dtype=np.float64).reshape(-1)
                x_full_arr = _expand_alpha(x_var_arr)

                metrics = opt.get_metrics(x_full_arr, evaluate_if_missing=True)
                entry = {
                    "tag": str(tag),
                    "iteration": int(iter_counter["value"]),
                    "alphas_trainable": x_var_arr.tolist(),
                    "loss_total": float(metrics["loss_total"]) if metrics else math.nan,
                    "loss_energy": float(metrics["loss_energy"]) if metrics else math.nan,
                    "loss_force": float(metrics["loss_force"]) if metrics else math.nan,
                }
                try:
                    jac_full = opt.jac(x_full_arr)
                    entry["jac_l2"] = float(np.linalg.norm(jac_full[trainable_idx]))
                except Exception:
                    entry["jac_l2"] = math.nan

                for key, value in extra.items():
                    if isinstance(value, (np.floating, float)):
                        entry[str(key)] = float(value)
                    elif isinstance(value, (np.integer, int)):
                        entry[str(key)] = int(value)
                    elif isinstance(value, (np.bool_, bool)):
                        entry[str(key)] = bool(value)
                    else:
                        entry[str(key)] = value
                history.append(entry)

                if not progress_header_printed["value"]:
                    title = "Trust-Radius Optimization Progress"
                    header = (
                        f"{'Basin':>7}"
                        f"{'Stage':>16}"
                        f"{'Iter':>8}"
                        f"{'Loss':>16}"
                        f"{'Energy Loss':>16}"
                        f"{'Force Loss':>16}"
                        f"{'Grad Norm':>16}"
                    )

                    self.ostream.print_blank()
                    self.ostream.print_header(title)
                    self.ostream.print_header("=" * len(title))
                    self.ostream.print_blank()
                    self.ostream.print_header(header)
                    self.ostream.print_header("-" * len(header))
                    progress_header_printed["value"] = True

                self.ostream.print_header(
                    f"{basin_label:>7}"
                    f"{entry['tag']:>16}"
                    f"{entry['iteration']:8d}"
                    f"{_fmt_float(entry['loss_total']):>16}"
                    f"{_fmt_float(entry['loss_energy']):>16}"
                    f"{_fmt_float(entry['loss_force']):>16}"
                    f"{_fmt_float(entry['jac_l2']):>16}"
                )
                self.ostream.flush()

            def _local_minimizer_callback(xk):
                iter_counter["value"] += 1
                _append_history("local_iter", xk)

            def _basinhop_callback(x, f, accept):
                _append_history("basinhop_step", x, basinhop_fun=float(f), basinhop_accept=bool(accept))
                return False

            # reduced-space objective/gradient
            def fun_var(x):
                return opt.fun_reduced(x, trainable_idx, alpha_full_start)

            def jac_var(x):
                return opt.jac_reduced(x, trainable_idx, alpha_full_start)

            minimizer_kwargs = {
                "method": "L-BFGS-B",
                "jac": jac_var,
                "bounds": bounds,
                "callback": _local_minimizer_callback,
                "options": {"gtol": gtol, "ftol": ftol, "maxls": 10, "maxiter": maxiter},
            }
            saved_radii = [dp.confidence_radius for dp in dps]

            try:
                _append_history("initial", x0)

                if n_bh_iter > 0:
                    res = basinhopping(
                        fun_var,
                        x0=x0,
                        minimizer_kwargs=minimizer_kwargs,
                        niter=n_bh_iter,
                        callback=_basinhop_callback,
                        disp=False
                    )
                else:
                    res = minimize(
                        fun_var,
                        x0=x0,
                        jac=jac_var,
                        method="L-BFGS-B",
                        bounds=bounds,
                        callback=_local_minimizer_callback,
                        options=minimizer_kwargs["options"]
                    )

                x_var_opt = np.asarray(res.x, dtype=np.float64).copy()
                _append_history("final", x_var_opt)

                # convert reduced solution -> full alpha vector
                x_full_opt = _expand_alpha(x_var_opt)

                setattr(res, "x_var", x_var_opt)
                res.x = x_full_opt

                setattr(res, "history", history)
                setattr(res, "component_summary", opt.get_metrics(x_full_opt, evaluate_if_missing=True))

                self.ostream.print_blank()
                self.ostream.print_header('Trust-Radius Optimization Result')
                self.ostream.print_header('================================')
                self.ostream.print_blank()
                self.ostream.print_header(f"{'Final objective':<32}{float(res.fun): .6f}")
                self.ostream.print_header('Optimized variables')
                self.ostream.print_block(np.array2string(x_var_opt, precision=6, max_line_width=100, threshold=sys.maxsize))
                self.ostream.print_blank()
                self.ostream.flush()
                return res
            finally:
                for dp, radius in zip(dps, saved_radii):
                    dp.confidence_radius = radius
                opt.close()

        initial_alphas = np.array(
            [float(np.asarray(dp.confidence_radius).reshape(-1)[0]) for dp in datapoints],
            dtype=float
        )

        has_last = len(datapoints) >= 2
        if has_last:
            dp_wo_last = datapoints[:-1]
            alphas_wo_last = initial_alphas[:-1]
        else:
            dp_wo_last = datapoints
            alphas_wo_last = initial_alphas

        stage_results = []

        def evaluate_stages_a_b():
            dp_eval_a = self._clone_datapoints_with_alphas(dp_wo_last, alphas_wo_last)
            stage_results.append(
                self._evaluate_interpolation_errors_on_reference_set(
                    scenario_name="without_last_datapoint",
                    molecules=molecules,
                    qm_energies=qm_energies,
                    qm_gradients=qm_gradients,
                    datapoints_eval=dp_eval_a,
                    z_matrix=z_matrix,
                    interpolation_setting=interpolation_setting,
                    sym_dict=sym_dict,
                    exponent_p_q=exponent_p_q,
                )
            )

            averaged_alphas = [0.5 for i in range(len(initial_alphas))]

            dp_eval_b = self._clone_datapoints_with_alphas(datapoints, averaged_alphas)
            stage_results.append(
                self._evaluate_interpolation_errors_on_reference_set(
                    scenario_name="with_new_datapoint_pre_opt",
                    molecules=molecules,
                    qm_energies=qm_energies,
                    qm_gradients=qm_gradients,
                    datapoints_eval=dp_eval_b,
                    z_matrix=z_matrix,
                    interpolation_setting=interpolation_setting,
                    sym_dict=sym_dict,
                    exponent_p_q=exponent_p_q,
                )
            )

        # trust_radius_fit_cache: the grouping weights come from the cached raw
        # weights at radius 0.5, the errors and the fits from one evaluation
        # per reference (IMTrustRadiusPairTerms)
        use_cache = bool(self.trust_radius_fit_cache)
        if use_cache:
            reason = self._trust_radius_cache_unsupported(
                self._build_interp_driver_for_trust_radius_eval(
                    z_matrix=z_matrix,
                    interpolation_setting=interpolation_setting,
                    sym_dict=sym_dict,
                    exponent_p_q=exponent_p_q,
                    datapoints_eval=datapoints,
                ),
                sym_dict,
                datapoints[0].cartesian_coordinates.shape[0],
            )
            if reason is not None:
                self.ostream.print_warning(f'trust_radius_fit_cache is not used: {reason}.')
                self.ostream.flush()
                use_cache = False

        if use_cache:
            time_start = time()
            w_dpref = self._trust_radius_grouping_weights(
                molecules, datapoints, interpolation_setting, sym_dict, z_matrix, exponent_p_q)
            self.ostream.print_info(
                f'Grouping weights at radius 0.5 from the cache: {time() - time_start:.1f} s')
            self.ostream.flush()
        else:
            evaluate_stages_a_b()
            w_dpref = stage_results[-1]['weights_refs']
        n_dp, n_ref = w_dpref.shape

        support_mask = np.zeros((n_dp, n_ref), dtype=bool)
        selected_weight_matrix = np.zeros((n_dp, n_ref), dtype=np.float64)
        reference_sample_idx = np.arange(len(molecules), dtype=int)
        coverage_mass = 0.9
        min_supported_references = 4
        min_shared_references = 4
        min_overlap_fraction = 0.40
        min_references_per_group = 10
        min_datapoints_per_group = 2

        for j in range(n_ref):
            weight_col_agg = w_dpref[:, j]

            finite_mask = np.isfinite(weight_col_agg) & (weight_col_agg > 0.0)
            idx = np.where(finite_mask)[0]
            order = idx[np.argsort(weight_col_agg[idx])[::-1]]
            sorted_weight_col_agg = weight_col_agg[order]
            cumulative = np.cumsum(sorted_weight_col_agg)

            target = float(np.clip(coverage_mass, 0.0, 1.0) * cumulative[-1])
            cutoff_pos = int(np.searchsorted(cumulative, target, side="left"))
            selected_dp = order[:cutoff_pos + 1]

            if selected_dp.size == 0:
                continue
            support_mask[selected_dp, j] = True
            selected_weight_matrix[selected_dp, j] = w_dpref[selected_dp, j]
        support_int = support_mask.astype(np.int32)
        shared_reference_counts = (support_int @ support_int.T)
        support_counts = np.sum(support_int, axis=1).astype(np.int32)

        overlap_denominator = np.minimum.outer(support_counts, support_counts).astype(np.float64)
        normalized_overlap = np.divide(
            shared_reference_counts.astype(np.float64),
            overlap_denominator,
            out=np.zeros_like(shared_reference_counts, dtype=np.float64),
            where=overlap_denominator > 0.0,
        )
        active_dp_mask = support_counts >= min_supported_references
        adjacency_dp = (
            (shared_reference_counts >= min_shared_references)
            & (normalized_overlap >= min_overlap_fraction)
            & active_dp_mask[:, None]
            & active_dp_mask[None, :]
        )

        np.fill_diagonal(adjacency_dp, False)

        graph = csr_matrix(adjacency_dp.astype(np.int8))
        _, dp_labels = connected_components(csgraph=graph, directed=False, return_labels=True)

        basin_groups = []
        active_indices = np.where(active_dp_mask)[0]

        for component_id in np.unique(dp_labels[active_indices]):
            dp_idx = np.where((dp_labels == component_id) & active_dp_mask)[0]

            if dp_idx.size < min_datapoints_per_group:
                continue

            # Reference union of this datapoint component:
            # keep every sampled reference supported by at least one datapoint
            # in the component.
            ref_local_idx = np.where(np.any(support_mask[dp_idx], axis=0))[0]

            if ref_local_idx.size < min_references_per_group:
                continue

            basin_groups.append({
                "basin_id": int(len(basin_groups)),
                "datapoint_indices": dp_idx.astype(int).tolist(),
                "reference_indices": reference_sample_idx[ref_local_idx].astype(int).tolist(),
            })

        groups = {
            "groups": basin_groups,
            "diagnostics": {
                "mode": "normalized_support_overlap",
                "coverage_mass": float(coverage_mass),
                "min_supported_references": int(min_supported_references),
                "min_shared_references": int(min_shared_references),
                "min_overlap_fraction": float(min_overlap_fraction),
                "n_datapoints": int(n_dp),
                "n_reference_samples": int(n_ref),
                "n_active_datapoints": int(np.sum(active_dp_mask)),
                "n_groups": int(len(basin_groups)),
            },
            "support_mask": support_mask,
            "selected_weight_matrix": selected_weight_matrix,
            "shared_reference_counts": shared_reference_counts,
            "normalized_overlap": normalized_overlap,
            "support_counts": support_counts,
        }

        final_alphas = initial_alphas.copy()

        self.ostream.print_blank()
        self.ostream.print_header('Trust-Radius Optimization Groups')
        self.ostream.print_header('================================')
        self.ostream.print_blank()
        for group_idx, group in enumerate(groups.get('groups', []), start=1):
            dp_idx = np.asarray(group.get('datapoint_indices', []), dtype=int)
            ref_idx = np.asarray(group.get('reference_indices', []), dtype=int)
            basin_id = group.get('basin_id', group_idx)

            self.ostream.print_header(f'Group {group_idx}')
            self.ostream.print_header('-' * (len(str(group_idx)) + 6))
            self.ostream.print_header(f"{'Basin':<24}{int(basin_id)}")
            self.ostream.print_header(f"{'Datapoints':<24}{dp_idx.size}")
            self.ostream.print_header(f"{'References':<24}{ref_idx.size}")
            self.ostream.print_header('Datapoint indices')
            self.ostream.print_block(np.array2string(dp_idx, separator=', ', max_line_width=100, threshold=sys.maxsize))
            self.ostream.print_header('Reference indices')
            self.ostream.print_block(np.array2string(ref_idx, separator=', ', max_line_width=100, threshold=sys.maxsize))
            self.ostream.print_blank()

        self.ostream.flush()

        pair_terms = None
        if use_cache:
            free_idx = sorted({int(j) for group in groups['groups'] for j in group['datapoint_indices']})
            natoms = datapoints[0].cartesian_coordinates.shape[0]
            needed = IMTrustRadiusPairTerms.memory_estimate(len(molecules), natoms, len(free_idx))
            available = self._available_memory_bytes()
            message = None
            if available is not None and needed > 0.5 * available:
                message = (f'the terms need {needed / 2**30:.1f} GiB, more than half of the '
                           f'available {available / 2**30:.1f} GiB')
            else:
                time_start = time()
                pair_terms, result_wo_last, result_half, message = self._trust_radius_pair_terms(
                    molecules=molecules,
                    qm_energies=qm_energies,
                    qm_gradients=qm_gradients,
                    datapoints=datapoints,
                    initial_alphas=initial_alphas,
                    free_idx=free_idx,
                    interpolation_setting=interpolation_setting,
                    sym_dict=sym_dict,
                    z_matrix=z_matrix,
                    exponent_p_q=exponent_p_q,
                )
            if pair_terms is None:
                self.ostream.print_warning(f'trust_radius_fit_cache is not used: {message}.')
                self.ostream.flush()
                evaluate_stages_a_b()
            else:
                result_half['weights_refs'] = w_dpref
                stage_results.extend([result_wo_last, result_half])
                self.ostream.print_info(
                    f'Cached radius fit: {len(molecules)} references, {len(free_idx)} of {n_dp} '
                    f'datapoints free, {needed / 2**20:.0f} MiB, built in {time() - time_start:.1f} s')
                self.ostream.flush()

        for group in groups['groups']:

            dp_mask = np.asarray(group['datapoint_indices'], dtype=int)
            ref_mask = group['reference_indices']

            # Group decides active set; no extra freeze policy
            trainable_idx = dp_mask
            if trainable_idx.size == 0:
                continue

            molecules_sub = [molecules[i] for i in ref_mask]
            qm_energies_sub = [qm_energies[i] for i in ref_mask]
            qm_gradients_sub = [qm_gradients[i] for i in ref_mask]

            # IMPORTANT: pass full final_alphas + full datapoints
            res = optimize_trust_radius(
                alpha_full_start=final_alphas,
                trainable_idx=trainable_idx,
                geom_list=molecules_sub,
                E_ref_list=qm_energies_sub,
                G_ref_list=qm_gradients_sub,
                dps=datapoints,  # full set for interpolation physics
                impes_dict=interpolation_setting,
                sym_dict=sym_dict,
                exponent_p_q=exponent_p_q,
                basin_id=group["basin_id"],
                pair_terms=pair_terms,
                ref_rows=ref_mask,
            )

            final_alphas = np.asarray(res['x'], dtype=np.float64).copy()

        if pair_terms is not None:
            im_energies_post, im_gradients_post = pair_terms.energy_gradient(np.arange(len(molecules)), final_alphas)
            stage_results.append(
                self._interpolation_error_summary(
                    scenario_name="with_new_datapoint_post_opt",
                    n_dp=n_dp,
                    im_energies=im_energies_post,
                    im_gradients=im_gradients_post,
                    qm_energies=qm_energies,
                    qm_gradients=qm_gradients,
                    weights_refs=None,
                )
            )
        else:
            dp_eval_c = self._clone_datapoints_with_alphas(datapoints, final_alphas)
            stage_results.append(
                self._evaluate_interpolation_errors_on_reference_set(
                    scenario_name="with_new_datapoint_post_opt",
                    molecules=molecules,
                    qm_energies=qm_energies,
                    qm_gradients=qm_gradients,
                    datapoints_eval=dp_eval_c,
                    z_matrix=z_matrix,
                    interpolation_setting=interpolation_setting,
                    sym_dict=sym_dict,
                    exponent_p_q=exponent_p_q,
                )
            )

        self._print_trust_radius_stage_comparison(stage_results)

        return final_alphas

    def _compute_energy(self, qm_driver, molecule, basis=None):
        """
        Computes the QM energy with the given driver.

        :param qm_driver:
            The energy driver (XtbDriver, ScfRestrictedDriver or
            ScfUnrestrictedDriver).
        :param molecule:
            The molecule.
        :param basis:
            The basis set (SCF drivers only).

        :return:
            The energy as an array of one value (Hartree), the SCF results (SCF
            drivers) and the response results (always None).
        """
        # Determine the type of energy driver, to be able to
        # call it correctly.

        qm_energy = None
        scf_results = None
        rsp_results = None

        # XTB
        if isinstance(qm_driver, XtbDriver):
            qm_driver.ostream.mute()
            qm_driver.compute(molecule)
            qm_energy = qm_driver.get_energy()
            qm_driver.ostream.unmute()
            if qm_energy is None:
                raise RuntimeError('XTB energy is None on this rank after MPI synchronization.')
            qm_energy = np.array([qm_energy])

        # restricted SCF
        elif isinstance(qm_driver, ScfRestrictedDriver) or isinstance(qm_driver, ScfUnrestrictedDriver):
            qm_driver.ostream.mute()
            scf_results = qm_driver.compute(molecule, basis)
            qm_energy = np.array([qm_driver.scf_energy])
            qm_driver.ostream.unmute()
            qm_driver.filename = None
            qm_driver.checkpoint_file = None

        if qm_energy is None:
            error_txt = "Could not compute the QM energy. "
            error_txt += "Please define a QM driver."
            raise ValueError(error_txt)

        return qm_energy, scf_results, rsp_results

    def _compute_gradient(self, grad_driver, molecule, basis=None, scf_results=None, rsp_results=None):
        """
        Computes the QM gradient with the given driver.

        :param grad_driver:
            The gradient driver (XtbGradientDriver or ScfGradientDriver).
        :param molecule:
            The molecule.
        :param basis:
            The basis set (SCF drivers only).
        :param scf_results:
            The SCF results (SCF drivers only).
        :param rsp_results:
            Not used (kept for symmetry with the excited-state drivers).

        :return:
            The gradient as an array of shape (1, natoms, 3) (Hartree/bohr).
        """

        qm_gradient = None

        if isinstance(grad_driver, XtbGradientDriver):
            grad_driver.ostream.mute()
            grad_driver.compute(molecule)
            qm_gradient = grad_driver.gradient
            qm_gradient = np.array([qm_gradient])
            grad_driver.ostream.unmute()

        elif isinstance(grad_driver, ScfGradientDriver):
            grad_driver.ostream.mute()
            grad_driver.compute(molecule, basis, scf_results)
            qm_gradient = grad_driver.gradient
            qm_gradient = np.array([qm_gradient])
            grad_driver.ostream.unmute()

        if qm_gradient is None:
            error_txt = "Could not compute the QM gradient. "
            error_txt += "Please define a QM gradient driver."
            raise ValueError(error_txt)

        return qm_gradient

    # TODO: mute outside to save time?
    def _compute_hessian(self, hess_driver, molecule, basis=None):
        """
        Computes the QM Hessian with the given driver.

        :param hess_driver:
            The Hessian driver (XtbHessianDriver or ScfHessianDriver).
        :param molecule:
            The molecule.
        :param basis:
            The basis set (SCF drivers only).

        :return:
            The Hessian as an array of shape (1, 3N, 3N) (Hartree/bohr**2).
        """

        qm_hessians = None

        if isinstance(hess_driver, XtbHessianDriver):
            hess_driver.ostream.mute()
            hess_driver.compute(molecule)
            qm_hessian = hess_driver.hessian
            hess_driver.ostream.unmute()

            if qm_hessian is None:
                raise RuntimeError('XTB Hessian is None on this rank after MPI synchronization.')
            qm_hessians = np.array([qm_hessian])

        elif isinstance(hess_driver, ScfHessianDriver):
            hess_driver.ostream.mute()
            hess_driver.compute(molecule, basis)
            qm_hessian = hess_driver.hessian
            qm_hessians = np.array([qm_hessian])
            hess_driver.ostream.unmute()

        if qm_hessians is None:
            error_txt = "Could not compute the QM Hessian. "
            error_txt += "Please define a QM Hessian driver."
            raise ValueError(error_txt)

        return qm_hessians

    def write_qm_energy_determined_points_h5(
        self,
        molecules,
        qm_energies,
        qm_gradients,
        im_energies,
        state,
        h5_path=None,
    ):
        """
        Appends reference structures to the reference file. Layout, one entry per
        structure in /states/<state>/: qm_energy (Hartree), im_energy (Hartree),
        natoms, xyz (xyz string, angstrom), qm_grad_flat (Hartree/bohr, natoms*3
        values) and coords_flat (angstrom, natoms*3 values); the file attributes
        natoms, labels, charge and multiplicity identify the QM region.

        :param molecules:
            The structures.
        :param qm_energies:
            The QM energies (Hartree).
        :param qm_gradients:
            The QM gradients, shape (natoms, 3) each (Hartree/bohr).
        :param im_energies:
            The interpolated energies (Hartree).
        :param state:
            The state (root).
        :param h5_path:
            The file; default: the reference file in use (_reference_file_path).
        """

        def _ensure_state_group(h5: h5py.File, state: int, compression="gzip", compression_opts=4):
            """
            Ensure /states/<state>/ datasets exist and are appendable (maxshape=None).
            """
            states_grp = h5.require_group("states")
            g = states_grp.require_group(str(int(state)))

            str_dt = h5py.string_dtype(encoding="utf-8")
            vlen_f64 = h5py.vlen_dtype(np.dtype("float64"))

            def req_1d(name, dtype):
                if name in g:
                    return g[name]
                # chunks=True lets HDF5 pick a reasonable chunking for 1D extendable arrays
                return g.create_dataset(
                    name,
                    shape=(0,),
                    maxshape=(None,),
                    dtype=dtype,
                    chunks=True,
                    compression=compression,
                    compression_opts=compression_opts,
                    shuffle=True,
                )

            req_1d("qm_energy", np.float64)
            req_1d("im_energy", np.float64)
            req_1d("natoms", np.int32)
            req_1d("xyz", str_dt)
            req_1d("qm_grad_flat", vlen_f64)
            req_1d("coords_flat", vlen_f64)

            return g

        h5_path = h5_path or self._reference_file_path()
        state = int(state)
        n = len(molecules)
        if not (len(qm_energies) == len(qm_gradients) == len(im_energies) == n):
            raise ValueError("Input lists must all have the same length.")

        with h5py.File(h5_path, "a") as h5:
            if n > 0 and 'labels' not in h5.attrs:
                self._set_reference_metadata(h5, molecules[0].get_labels(),
                                             molecules[0].get_charge(),
                                             molecules[0].get_multiplicity())
            g = _ensure_state_group(h5, state)
            N0 = g["qm_energy"].shape[0]
            N1 = N0 + n

            for key in ["qm_energy", "im_energy", "natoms", "xyz", "qm_grad_flat", "coords_flat"]:
                g[key].resize((N1,))

            for i, mol in enumerate(molecules):

                xyz = mol.get_xyz_string()
                coords = np.asarray(mol.get_coordinates_in_angstrom(), dtype=np.float64)
                natoms = int(mol.number_of_atoms())

                grad = np.asarray(qm_gradients[i], dtype=np.float64)
                if grad.ndim != 2 or grad.shape[1] != 3:
                    raise ValueError(f"qm_gradients[{i}] must have shape (natoms,3), got {grad.shape}.")
                if grad.shape[0] != natoms:
                    raise ValueError(f"Gradient natoms mismatch at i={i}: grad has {grad.shape[0]}, expected {natoms}.")

                row = N0 + i
                g["qm_energy"][row] = float(qm_energies[i])
                g["im_energy"][row] = float(im_energies[i])
                g["natoms"][row] = natoms
                g["xyz"][row] = xyz
                g["qm_grad_flat"][row] = grad.reshape(-1)
                g["coords_flat"][row] = coords.reshape(-1)

    def _reference_file_path(self):
        """
        Returns the reference file in use: the user-given
        reference_struc_energies_file, or QM_ref_along_traj.h5 in the working
        directory when none was given.

        :return:
            The path of the reference file.
        """

        if self.reference_struc_energies_file is not None:
            return self.reference_struc_energies_file
        return self._DEFAULT_REFERENCE_FILE

    @staticmethod
    def _set_reference_metadata(h5, labels, charge, multiplicity):
        """
        Stores the identity of the QM region as attributes of the reference file.

        :param h5:
            The open h5py file.
        :param labels:
            The element labels of the QM region.
        :param charge:
            The molecular charge.
        :param multiplicity:
            The spin multiplicity.
        """

        h5.attrs['natoms'] = int(len(labels))
        h5.attrs['labels'] = np.array([str(label) for label in labels],
                                      dtype=h5py.string_dtype(encoding='utf-8'))
        h5.attrs['charge'] = float(charge)
        h5.attrs['multiplicity'] = int(multiplicity)

    def _write_reference_file_metadata(self, filename, labels, charge, multiplicity):
        """
        Creates the reference file if it does not exist and stores the identity
        of the QM region in it, unless it is already there.

        :param filename:
            The reference file.
        :param labels:
            The element labels of the QM region.
        :param charge:
            The molecular charge.
        :param multiplicity:
            The spin multiplicity.
        """

        with h5py.File(filename, 'a') as h5:
            h5.require_group('states')
            if 'labels' not in h5.attrs:
                self._set_reference_metadata(h5, labels, charge, multiplicity)

    def extract_reference_structures_h5(self, filename, roots, expected_labels=None,
                                        expected_charge=None, expected_multiplicity=None):
        """
        Reads the reference structures of the given roots from a reference file
        (layout in write_qm_energy_determined_points_h5).

        :param filename:
            The reference file.
        :param roots:
            The roots to read.
        :param expected_labels:
            If given, the element labels of the file and of every entry must match
            them, otherwise the file is refused.
        :param expected_charge:
            If given, must match the charge stored in the file (if stored).
        :param expected_multiplicity:
            If given, must match the multiplicity stored in the file (if stored).

        :return:
            {root: {'molecules', 'qm_energies', 'qm_gradients' (natoms, 3 arrays),
            'im_energies'}}.
        """
        reference_molecules_check = {
            int(root): {"molecules": [], "qm_energies": [], "qm_gradients": [], "im_energies": []}
            for root in roots
        }

        def _normalized(labels):
            return [str(label).strip().capitalize() for label in labels]

        with h5py.File(filename, "r") as h5:
            if expected_labels is not None and 'labels' in h5.attrs:
                stored_labels = [label.decode('utf-8') if isinstance(label, bytes) else label
                                 for label in h5.attrs['labels']]
                assert_msg_critical(
                    _normalized(stored_labels) == _normalized(expected_labels),
                    f'IMDatabasePointCollecter: the reference file {filename} '
                    'belongs to a different molecule (element labels differ). '
                    'Remove it or set reference_struc_energy_file to another file.')
            if expected_charge is not None and 'charge' in h5.attrs:
                assert_msg_critical(
                    float(h5.attrs['charge']) == float(expected_charge),
                    f'IMDatabasePointCollecter: the reference file {filename} '
                    'was written for a different molecular charge.')
            if expected_multiplicity is not None and 'multiplicity' in h5.attrs:
                assert_msg_critical(
                    int(h5.attrs['multiplicity']) == int(expected_multiplicity),
                    f'IMDatabasePointCollecter: the reference file {filename} '
                    'was written for a different spin multiplicity.')

            if "states" not in h5:
                return reference_molecules_check

            states_grp = h5["states"]

            for root in roots:
                s = str(int(root))
                if s not in states_grp:
                    continue

                g = states_grp[s]
                qmE = g["qm_energy"][:]
                imE = g["im_energy"][:]
                natoms = g["natoms"][:]
                xyzs = g["xyz"][:]
                grad_flat = g["qm_grad_flat"][:]

                for i in range(len(qmE)):
                    xyz_string = xyzs[i]
                    if isinstance(xyz_string, bytes):
                        xyz_string = xyz_string.decode("utf-8")

                    mol = Molecule.from_xyz_string(xyz_string)
                    if expected_labels is not None:
                        assert_msg_critical(
                            _normalized(mol.get_labels()) == _normalized(expected_labels),
                            f'IMDatabasePointCollecter: entry {i} of state {root} in the '
                            f'reference file {filename} belongs to a different molecule '
                            '(element labels differ). Remove the file or set '
                            'reference_struc_energy_file to another file.')
                    n = int(natoms[i])
                    gf = np.asarray(grad_flat[i], dtype=np.float64)
                    if gf.size != 3 * n:
                        raise ValueError(f"Stored gradient size mismatch at state={root}, i={i}: got {gf.size}, expected {3*n}.")
                    grad = gf.reshape((n, 3))

                    reference_molecules_check[int(root)]["molecules"].append(mol)
                    reference_molecules_check[int(root)]["qm_energies"].append(float(qmE[i]))
                    reference_molecules_check[int(root)]["qm_gradients"].append(grad)
                    reference_molecules_check[int(root)]["im_energies"].append(float(imE[i]))

        return reference_molecules_check

    def calculate_translation_coordinates(self, given_coordinates):
        """
        Moves the geometric center of the coordinates to the origin.

        :param given_coordinates:
            The coordinates, shape (natoms, 3).

        :return:
            The translated coordinates.
        """
        center = np.mean(given_coordinates, axis=0)
        translated_coordinates = given_coordinates - center

        return translated_coordinates
