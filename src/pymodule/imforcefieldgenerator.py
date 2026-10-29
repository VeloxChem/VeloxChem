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

import numpy as np
import numbers
import random
import itertools
import h5py
import json
import time
from mpi4py import MPI
from contextlib import redirect_stderr
from pathlib import Path
from .outputstream import OutputStream
from io import StringIO
from copy import deepcopy

from .xtbdriver import XtbDriver
from .xtbgradientdriver import XtbGradientDriver
from .xtbhessiandriver import XtbHessianDriver
from .scfrestdriver import ScfRestrictedDriver
from .scfunrestdriver import ScfUnrestrictedDriver
from .scfgradientdriver import ScfGradientDriver
from .scfhessiandriver import ScfHessianDriver
from .tdaeigensolver import TdaEigenSolver
from .lreigensolver import LinearResponseEigenSolver
from .molecularbasis import MolecularBasis
from .interpolationdriver import InterpolationDriver
from .interpolationdatapoint import InterpolationDatapoint
from .imdatabasepointcollecter import IMDatabasePointCollecter
from .mmforcefieldgenerator import MMForceFieldGenerator
from .conformergenerator import ConformerGenerator
from .optimizationdriver import OptimizationDriver
from .molecule import Molecule
from .superimpose import svd_superimpose
from .errorhandler import assert_msg_critical
from .veloxchemlib import hartree_in_kcalpermol, bohr_in_angstrom, boltzmann_in_hartreeperkelvin
from ._version import __version__


with redirect_stderr(StringIO()) as fg_err:
    import geometric


class IMForceFieldGenerator:
    """
    Sets up and controls the construction of an interpolation (MSI, IM)
    database: the seed datapoints from optimised structures, then the
    collector dynamics (IMDatabasePointCollecter) on the interpolated PES,
    which adds datapoints where the interpolation deviates from QM, and a
    final quality check. Ground state, single MPI rank only.

    :param ground_state_driver:
        The QM energy driver: ScfRestrictedDriver, ScfUnrestrictedDriver or
        XtbDriver. The gradient and Hessian drivers are created from it.
    :param roots_to_follow:
        The electronic states (only [0] is supported; default [0]).
    :param comm:
        The MPI communicator (a single rank is supported).
    :param ostream:
        The output stream.

    Usage: set the attributes, call set_up_the_system(molecule) (it freezes
    the interpolation settings), then compute(molecule). Assigning an
    attribute the generator does not use (misspelled or removed setting)
    prints a warning at the start of compute.

    QM drivers and basis sets
        - drivers: {'gs': (energy, gradient, Hessian driver), 'es': None}.
        - sampling_driver: {'gs': the xTB triple of the sampling pre-screen}
          (the ground-state triple itself for XtbDriver).
        - gs_basis_set_label: Basis set of the ground state. Default:
          'def2-svp' (ignored by xTB).
        - es_basis_set_label: Basis set of excited states (not used).
        - open_mm_platform: OpenMM platform of the collector dynamics, e.g.
          'CPU' or 'CUDA'. Default: None (OpenMM chooses).

    Interpolation settings (frozen by set_up_the_system; changing them later
    has no effect and is reported)
        - interpolation_type: 'shepard'.
        - weightfunction_type: Weight distance, 'cartesian' (default) or
          'internal'.
        - exponent_p, exponent_q: Exponents of the Shepard weight
          w = 1/(u^p + u^q), u = (d/R)^2. Default: 6 and 2 (int; strings
          such as '8' are cast by the driver).
        - confidence_radius: Default radius R of the driver settings, bohr.
          Default: 0.5 (new datapoints use use_opt_confidence_radius[2]).
        - use_inverse_bond_length: Bonds as 1/r. Default: True.
        - use_eq_bond_length: Bonds relative to equilibrium lengths, taken
          from the first seed datapoint (the bonds added for a reaction from
          the last structure of the seed batch); needs
          use_inverse_bond_length = False, which otherwise takes precedence.
          Default: False.
        - use_tc_weights: Gate the weights with the important internal
          coordinates of the datapoints. Default: True.
        - tc_weight_mode: 'multiplicative' (default) or 'additive_rhee'.
        - use_mass_weight: Mass-weighted Cartesian gradients and Hessians.
          Default: True.
        - exclude_non_core: Exclude H atoms (and non-core atoms) from the
          weight distance. Default: False.
        - consider_locality: Allow more constraints in the constrained
          optimisation of new datapoints, which keeps them closer to the
          sampled structure (forwarded to the collector). Default: False.

    Seed stage (only for a new database; skipped, with a warning, when the
    database file exists)
        - use_minimized_structures: [optimise the seed structures (default
          True), global optimisation constraints (1-based atom-index tuples of
          length 2, 3, 4 or geomeTRIC constraint strings, default []), unused].
          The global constraints hold in every seed optimisation and network
          scan; the tuple constraints of a seed become the important internal
          coordinates of its datapoint.
        - add_conformal_structures: Seed with the conformers of the MM
          conformer generator. Default: True.
        - reaction_structures: Reaction-path structures, or None (default).
          Each entry is a Molecule, an XYZ string, [molecule, mode,
          constraints] or {'molecule', 'mode', 'constraints'}; mode 'normal',
          'transition' or 'constraint'; constraints 1-based.
        - seed_duplicate_rmsd_thrsh: Optimised seeds with an all-atom RMSD
          below this value (Angstrom) are the same minimum and are added only
          once. Default: None (off; required by the network stage).
        - seed_scan_clash_relief: Rotate neighbouring rotatable bonds when a
          scan seed cannot reach its dihedral by a rigid rotation. Default:
          False.
        - drop_improper_h_dihedrals: Remove the proper dihedrals terminated by
          a hydrogen bonded to an improper centre from a z-matrix built from
          the molecule (the improper and the two valence angles carry the
          wag). False (default), True (every centre) or element symbols of the
          centres, e.g. ['N'] ('N' alone is one symbol). Scanned dihedrals and
          4-atom global constraints are kept.
        - seed_conformer_network: Opt-in network seed stage
          (_seed_conformer_network): every conformer is optimised once, the
          distinct minima are datapoints, minima connected by the rotation of
          one scanned bond are joined by relaxed scans from both ends, and scan
          points are added one at a time where the interpolation is worst; with
          reaction_structures the reaction is an edge of the network too.
          Needs use_minimized_structures[0], seed_duplicate_rmsd_thrsh and
          seed structures. Default: False.
        - seed_network_settings: Settings of the network stage, with the
          defaults below (change single entries; unknown or missing keys and
          invalid values abort the stage before the first optimisation):
            - scan_step_deg: Requested torsion scan step, degree. Default 15.0.
            - state_gap_deg: Smallest circular gap that separates two torsion
              states of a scanned bond, degree. Default 30.0.
            - e_thrsh_kcal_per_atom: Energy threshold of the acceptance score,
              kcal/mol per atom. Default None (energy_threshold).
            - g_rmsd_thrsh: Gradient threshold of the acceptance score, RMS
              over the 3N components, kcal/mol/Angstrom. Default None
              (gradient_rmsd_thrsh).
            - max_points_per_path: Most scan points added per path. Default 6.
            - max_rejections_per_path: Consecutive rejected trial points that
              stop a path. Default 3.
            - include_orphan_arcs: Also scan the arcs of minima without a
              partner minimum. Default False.
            - energy_window_kcal: Drop conformer minima and torsion-scan
              points more than this above the lowest minimum, kcal/mol.
              Default None (off).
            - energy_window_kT: The same window in units of k_B T at the
              dynamics temperature (set at most one window). Default None.
            - reaction_scan_step_angstrom: Step of the driven bond of a
              reaction edge, Angstrom. Default 0.05.
        - seed_structures, seed_conformers, seed_reaction_entries: Seed
          structures stored by set_up_the_system ({root: {1-based dihedral or
          None: [(molecule, mode[, constraints])]}}, the unique conformers, the
          reaction entries kept for the network stage).

    Collector dynamics and selection of new datapoints
        - temperature: Temperature, K. Default: 150.15.
        - ensemble: 'NVE' (default), 'NVT' or 'NPT'.
        - pressure: Pressure of NPT, atm. Default: 1.0.
        - friction: Langevin friction, 1/ps. Default: 1.0.
        - timestep: Time step, fs. Default: 0.5.
        - nsteps: Number of MD steps. Default: 1000.
        - snapshots: Number of saved frames (every nsteps/snapshots steps).
          Default: 1000.
        - duration: Length of one cycle, steps; converged_cycle cycles without
          a new datapoint end the run. Default: 2000.
        - converged_cycle: See duration. Default: 5.
        - start_collect: First step at which datapoints are collected.
          Default: 0.
        - energy_threshold: A structure becomes a datapoint if the
          interpolation error exceeds this energy, kcal/mol per atom.
          Default: 1.5.
        - gradient_rmsd_thrsh: ... or this gradient error, RMS over the 3N
          components, kcal/mol/Angstrom. Default: 1.0.
        - desired_point_density: Total number of datapoints at which the
          construction stops. Default: 50.
        - identify_relevant_int_coordinates: Constraint-optimise new
          datapoints along their important internal coordinates. Default:
          True. (identfy_relevant_int_coordinates is a deprecated alias.)
        - use_opt_confidence_radius: [fit the confidence radii to the
          reference structures (default False), fit mode ('multi_grad'),
          initial radius of every new datapoint in bohr (0.5), energy weight
          of the fit objective (0.3)].
        - gate_distance_prefilter: Exact speed-up of the proximity gate of the
          dynamics: a reference structure is aligned only if the lower bound
          sqrt(sum_{a<b} (d_ab - d'_ab)^2) / N of its RMSD (interatomic
          distances d) does not exclude it. Same decisions as the full scan.
          Default: False.
        - trust_radius_fit_cache: Exact speed-up of the confidence-radius fit:
          one interpolation per reference and fit instead of one per
          reference and objective evaluation, the grouping weights from a
          cache of the raw weights at radius 0.5, and the fit in closed form.
          The groups are the same, and the loss and gradient equal those of
          the default fit to round-off. Memory: 8 x references x free
          datapoints x (9 natoms + 3) bytes. Default: False.
        - sampling_settings: xTB pre-screen of the QM checks, {'enabled':
          False, 'e_thrsh_kcal_per_atom': 0.1} (no other key is read).
        - metadynamics_settings: Well-tempered metadynamics of the dynamics,
          or None (default): {'enabled', 'bias_factor', 'hill_height_kjmol',
          'hill_frequency', 'variables': [{'type': 'torsion', 'atoms' (1-based),
          'min_deg', 'max_deg', 'width_deg', 'periodic'}]}.
        - reference_struc_energy_file: HDF5 file of the QM checks along the
          trajectory; an existing file is preloaded by the collector. Default:
          None (the collector uses QM_ref_along_traj.h5).
        - load_system: OpenMM XML state (written by the collector as
          checkpoint.xml) to start the dynamics from. Default: None.
        - trajectory_file: Not used; the trajectory is written to
          trajectory_<root>.pdb.
        - solvent: 'gas' (the construction is in the gas phase).
        - force_constant: Forwarded to the collector, not used there.
        - nstruc_to_confirm_database_quality: Number of trajectory structures
          of the final quality check. Default: 15.

    Results and internal state
        - molecule: The molecule given to set_up_the_system.
        - roots_z_matrix: {root: z-matrix} with 'bonds', 'angles',
          'dihedrals', 'impropers' (0-based tuples), built from the molecule
          or read from the database.
        - symmetry_information: {'gs', 'es'} atom groups and rotatable bonds
          used by the interpolation driver.
        - all_rotatable_bonds: The 0-based rotatable bonds.
        - imforcefieldfiles, sampling_imforcefieldfiles: The database files
          per root.
        - states_interpolation_settings,
          sampling_states_interpolation_settings: The frozen settings dicts.
        - eq_bond_length, eq_bond_length_irc_bonds: Equilibrium bond lengths
          (bohr) of the eq bond mode.
        - dynamics_settings: The settings passed to the collector.
        - states_data_point_density, qmlabels, im_results: The number of
          datapoints per root, the labels, {'n_datapoints'} (returned by
          compute).
        - qm_energies, total_energies, kinetic_energies, unique_molecules,
          state_specific_molecules: Data of the collector runs.

    Files (in the working directory unless stated)
        - im_database_<root>.h5 (or the given name): the database; the HDF5
          file attribute 'msi_construction_history' records the construction
          settings.
        - <database stem>_org.h5: copy of the seed datapoints.
        - im_database_sampling_<root>.h5: xTB mirror (sampling enabled).
        - seed_network_log.json (next to the database): provenance of the
          network seed stage.
        - full_xyz_traj.xyz, random_xyz_structures.xyz: final quality check.
        - MOL.xml, MOL.pdb (conformer generator); trajectory_<root>.pdb,
          checkpoint.xml, pre_opt_database_strucs.xyz, qm_region*.xml/.pdb and
          the reference file (collector).
    """

    # settings of earlier versions that are no longer read; _check_user_settings
    # reports them with the reason (adapted from the K13 proposal)
    _removed_settings = {
        'qm_omp_threads': 'QM calculations use the number of threads given by OMP_NUM_THREADS.',
        'use_cosine_dihedral': 'Proper dihedrals always use the 2*sin(delta/2) displacement.',
        'force_orient_thrsh': 'New datapoints are selected by energy_threshold and gradient_rmsd_thrsh only.',
        'excitation_pulse': 'Only ground-state database construction is supported.',
        'distance_thrsh': 'The collector selects new datapoints by energy_threshold and gradient_rmsd_thrsh.',
        'sampling_structures': 'The seed structures come from the conformer generator (add_conformal_structures).',
        'datafile': 'The database files are given to set_up_the_system (imforcefieldfiles).',
        'density_of_datapoints': 'The number of datapoints per root is states_data_point_density.',
    }

    # keys of sampling_settings read by the xTB pre-screening of the collector
    _sampling_settings_keys = ('enabled', 'e_thrsh_kcal_per_atom')

    def __init__(self, ground_state_driver, roots_to_follow=None, comm=None, ostream=None):

        if comm is None:
            comm = MPI.COMM_WORLD

        if roots_to_follow is None:
            roots_to_follow = [0]

        assert_msg_critical(
            all(root == 0 for root in roots_to_follow),
            'Only ground-state potential construction is currently supported.'
        )

        self.comm = comm
        self.rank = comm.Get_rank()
        self.nodes = comm.Get_size()
        self.ostream = ostream if ostream is not None else OutputStream.create_mpi_ostream(comm)

        self.open_mm_platform = None

        self.qmlabels = None

        self.state_specific_molecules = None

        self.roots_z_matrix = {}
        self.symmetry_information = None
        self.all_rotatable_bonds = None

        self.reaction_structures = None
        self.seed_structures = None

        # rotate neighbouring rotatable bonds when a scan seed cannot reach
        # its target dihedral by a rigid rotation (steric clash)
        self.seed_scan_clash_relief = False

        # optimised conformer seeds with an all-atom RMSD (Angstrom) below this
        # value are the same minimum and are added only once (None: off)
        self.seed_duplicate_rmsd_thrsh = None

        # opt-in seed stage: optimise every conformer once, connect the distinct
        # minima by single-bond rotations and add relaxed-scan points where the
        # interpolation is worst (_seed_conformer_network); settings in the style
        # of sampling_settings, the unique conformers are stored by set_up_the_system.
        # seed_network_settings holds the defaults (None thresholds: the collector's
        # energy_threshold and gradient_rmsd_thrsh are used; the energy window, in
        # kcal/mol or in units of k_B T at the dynamics temperature, drops conformer
        # minima and torsion-scan points above the lowest minimum + window, None: off)
        self.seed_conformer_network = False
        self.seed_network_settings = {
            'scan_step_deg': 15.0,
            'state_gap_deg': 30.0,
            'e_thrsh_kcal_per_atom': None,
            'g_rmsd_thrsh': None,
            'max_points_per_path': 6,
            'max_rejections_per_path': 3,
            'include_orphan_arcs': False,
            'energy_window_kcal': None,
            'energy_window_kT': None,
            'reaction_scan_step_angstrom': 0.05,
        }
        # the keys of seed_network_settings, for the checks of _check_user_settings
        # and _validate_seed_network_settings
        self._seed_network_settings_keys = tuple(self.seed_network_settings)
        self.seed_conformers = None
        self._network_counters = None
        # reaction seeds kept for the network stage, which turns the reaction into network
        # edges (the conformer block of set_up_the_system replaces seed_structures)
        self.seed_reaction_entries = None

        # Here us the driver set up
        self.gs_basis_set_label = 'def2-svp'
        self.es_basis_set_label = '6-31g*'

        self.drivers = {'gs': None, 'es': None}
        self.sampling_driver = {'gs': None, 'es': None}

        assert_msg_critical(
            isinstance(ground_state_driver, (ScfRestrictedDriver, ScfUnrestrictedDriver, XtbDriver)),
            'IMForceFieldGenerator: the ground-state driver must be a ScfRestrictedDriver, '
            f'ScfUnrestrictedDriver or XtbDriver (got {type(ground_state_driver).__name__}).')

        if isinstance(ground_state_driver, (ScfRestrictedDriver, ScfUnrestrictedDriver)):
            qm_grad_driver = ScfGradientDriver(ground_state_driver)
            qm_hess_driver = ScfHessianDriver(ground_state_driver)

            self.drivers['gs'] = (ground_state_driver, qm_grad_driver, qm_hess_driver)

        # the xTB triple of the sampling pre-screen (the ground-state triple itself for xTB)
        if isinstance(ground_state_driver, XtbDriver):
            qm_grad_driver = XtbGradientDriver(ground_state_driver)
            qm_hess_driver = XtbHessianDriver(ground_state_driver)
            self.drivers['gs'] = (ground_state_driver, qm_grad_driver, qm_hess_driver)
            self.sampling_driver['gs'] = self.drivers['gs']

        else:
            qm_sampling_driver = XtbDriver()
            qm_sampling_grad_driver = XtbGradientDriver(qm_sampling_driver)
            qm_sampling_hess_driver = XtbHessianDriver(qm_sampling_driver)
            self.sampling_driver['gs'] = (qm_sampling_driver, qm_sampling_grad_driver, qm_sampling_hess_driver)

        self.states_interpolation_settings = {root: None for root in roots_to_follow}
        self.states_data_point_density = {root: None for root in roots_to_follow}
        self.roots_to_follow = roots_to_follow
        ##########################################################
        # variables for the interpolation
        self.interpolation_type = 'shepard'
        self.weightfunction_type = 'cartesian'
        self.exponent_p = 6
        self.exponent_q = 2
        self.confidence_radius = 0.5
        self.imforcefieldfiles = None
        self.use_inverse_bond_length = True
        self.use_eq_bond_length = False
        self.use_tc_weights = True
        self.tc_weight_mode = "multiplicative"  # "additive_rhee"
        self.use_mass_weight = True  # set True as it is standard in the YM scheme --> small differences
        self.consider_locality = False

        self.eq_bond_length = None
        self.eq_bond_length_irc_bonds = None

        # variables for the forcefield generation and database expansion
        self.dynamics_settings = None
        self.duration = 2000
        self.temperature = 150.15
        self.pressure = 1.0
        self.force_constant = 1.0
        self.ensemble = 'NVE'
        self.timestep = 0.5
        self.friction = 1.0
        self.nsteps = 1000
        self.snapshots = self.nsteps
        self.trajectory_file = 'trajectory.pdb'
        self.reference_struc_energy_file = None
        # OpenMM XML state (Simulation.saveState, written by the collector as
        # checkpoint.xml after every added point) to start the dynamics from
        self.load_system = None
        self.desired_point_density = 50
        self.converged_cycle = 5
        self.energy_threshold = 1.5
        self.gradient_rmsd_thrsh = 1.0
        self.start_collect = 0
        self.solvent = 'gas'

        # sampling settings
        self.sampling_settings = {
            'enabled': False,
            'e_thrsh_kcal_per_atom': 0.1,
        }

        self.sampling_imforcefieldfiles = None
        self.sampling_states_interpolation_settings = {root: None for root in roots_to_follow}
        # well-tempered metadynamics of the construction dynamics (layout in the class docstring)
        self.metadynamics_settings = None

        # results of the collector runs (one entry per compute call)
        self.qm_energies = []
        self.total_energies = []
        self.kinetic_energies = []
        self.unique_molecules = []

        # results of compute: the number of datapoints per root ('RMSD' and '|D|' are not filled)
        self.im_results = {'n_datapoints': None, 'RMSD': None, '|D|': None}

        # number of random trajectory structures of the final quality check
        self.nstruc_to_confirm_database_quality = 15

        # [optimise the seed structures, global constraints (1-based tuples or geomeTRIC strings), unused]
        self.use_minimized_structures = [True, [], []]
        # seed the database with the conformers of the MM conformer generator
        self.add_conformal_structures = True

        # constraint-optimise new datapoints of the collector along their important internal coordinates
        self.identify_relevant_int_coordinates = True
        self._identfy_alias_warned = False
        # [fit the confidence radii, fit mode, initial radius of every new datapoint (bohr), energy weight of the fit]
        self.use_opt_confidence_radius = [False, 'multi_grad', 0.5, 0.3]
        # exact speed-ups of the proximity gate and of the radius fit (same results)
        self.gate_distance_prefilter = False
        self.trust_radius_fit_cache = False
        # exclude H atoms and non-core atoms from the weight distance
        self.exclude_non_core = False

        # Opt-in z-matrix rule: remove the proper dihedrals terminated by a hydrogen that belongs to an improper
        # (the improper's centre is its first atom and must be bonded to that H). The out-of-plane motion of such
        # an H is then carried by the improper and the two valence angles at the centre instead of by proper
        # dihedrals with a 2 sin(delta/2) displacement, which saturate at large wag angles.
        #   False (default): the z-matrix is unchanged.
        #   True: apply the rule at every improper centre.
        #   list of element symbols, e.g. ['N']: apply it only at centres of those elements.
        # The rule acts only when the z-matrix is built from the molecule; a z-matrix read from an existing
        # database is never modified.
        self.drop_improper_h_dihedrals = False

        self._system_is_set_up = False

        # every attribute defined up to here (and self.molecule, assigned in
        # set_up_the_system) is a setting or internal state of the generator;
        # _check_user_settings reports any other attribute as unused
        self._known_attributes = set(vars(self)) | {'molecule', '_known_attributes'}

    def _warn_identfy_alias(self):
        """
        Warns once that identfy_relevant_int_coordinates is deprecated.
        """

        if not getattr(self, '_identfy_alias_warned', False):
            self.ostream.print_warning('IMForceFieldGenerator: identfy_relevant_int_coordinates is deprecated; '
                                       'use identify_relevant_int_coordinates.')
            self.ostream.flush()
            self._identfy_alias_warned = True

    @property
    def identfy_relevant_int_coordinates(self):
        """
        Deprecated alias of identify_relevant_int_coordinates.
        """

        self._warn_identfy_alias()
        return self.identify_relevant_int_coordinates

    @identfy_relevant_int_coordinates.setter
    def identfy_relevant_int_coordinates(self, value):
        """
        Deprecated alias of identify_relevant_int_coordinates.
        """

        self._warn_identfy_alias()
        self.identify_relevant_int_coordinates = value

    def _check_user_settings(self):
        """
        Reports settings that are assigned but have no effect on the
        construction: attributes the generator does not define (misspelled or
        removed settings), unread keys of sampling_settings and
        seed_network_settings, the legacy reaction bias, and interpolation
        settings changed after an explicit set_up_the_system (the
        interpolation uses the values frozen there). Warnings only. Adapted
        from the K13 proposal.

        :return:
            The list of warning messages.
        """

        messages = []

        known = getattr(self, '_known_attributes', None)
        if known is not None:
            for name in sorted(vars(self)):
                # methods replaced on the instance (tests, user hooks) are not settings
                if name in known or hasattr(type(self), name):
                    continue
                msg = f'IMForceFieldGenerator: setting "{name}" is not used and has no effect.'
                if name in self._removed_settings:
                    msg += ' ' + self._removed_settings[name]
                messages.append(msg)

        if isinstance(self.sampling_settings, dict):
            for key in self.sampling_settings:
                if key not in self._sampling_settings_keys:
                    messages.append(
                        f'IMForceFieldGenerator: sampling_settings["{key}"] is not used and has no effect. '
                        'The xTB pre-screening compares energies only.')

        if isinstance(self.seed_network_settings, dict):
            for key in self.seed_network_settings:
                if key not in self._seed_network_settings_keys:
                    messages.append(
                        f'IMForceFieldGenerator: seed_network_settings["{key}"] is not a setting of the '
                        f'network seed stage (known keys: {sorted(self._seed_network_settings_keys)}); the '
                        'network seed stage refuses to run with it.')

        # the interpolation settings are frozen in set_up_the_system; if it was
        # called explicitly, later attribute changes are not used
        if self._system_is_set_up:
            settings = self.states_interpolation_settings.get(self.roots_to_follow[0])
            for name in sorted(settings or {}):
                if name == 'imforcefield_file' or not hasattr(self, name):
                    continue
                if settings[name] != getattr(self, name):
                    messages.append(
                        f'IMForceFieldGenerator: "{name}" was changed after set_up_the_system; the '
                        f'interpolation uses the value at set-up time ({settings[name]}). Assign the '
                        'interpolation settings before calling set_up_the_system.')

        for msg in messages:
            self.ostream.print_warning(msg)
        if messages:
            self.ostream.flush()

        return messages

    @staticmethod
    def _drop_improper_h_dihedrals(zmat, labels, centre_elements=None, protected_dihedrals=None):
        """
        Removes the proper dihedrals whose terminal atom is a hydrogen that belongs to an improper
        (the improper's centre is its first atom and must be bonded to the H).

        :param zmat:
            The z-matrix dict with keys bonds, angles, dihedrals, impropers (0-based tuples).
        :param labels:
            The atom labels (element symbols) of the molecule.
        :param centre_elements:
            None (every improper centre) or a list of element symbols restricting the centres.
        :param protected_dihedrals:
            Optional 0-based dihedral tuples that are never removed (e.g. constraint dihedrals).

        :return:
            The new z-matrix dict, the removed dihedrals and the hydrogen atoms (0-based).
        """

        bonds = {tuple(sorted(int(a) for a in b)) for b in zmat['bonds']}
        protected = set()
        for d in (protected_dihedrals or []):
            d = tuple(int(x) for x in d)
            protected.add(d)
            protected.add(d[::-1])

        h_atoms = []
        for imp in zmat['impropers']:
            imp = [int(a) for a in imp]
            centre = imp[0]
            if centre_elements is not None and str(labels[centre]).capitalize() not in centre_elements:
                continue
            for a in imp[1:]:
                if (str(labels[a]).capitalize() == 'H' and tuple(sorted((centre, a))) in bonds and
                        a not in h_atoms):
                    h_atoms.append(a)

        removed = []
        kept = []
        for d_raw in zmat['dihedrals']:
            d = tuple(int(x) for x in d_raw)
            if (d[0] in h_atoms or d[3] in h_atoms) and d not in protected:
                removed.append(d)
            else:
                kept.append(d)

        new_zmat = {
            key: [tuple(int(x) for x in c) for c in zmat[key]]
            for key in ('bonds', 'angles', 'dihedrals', 'impropers')
        }
        new_zmat['dihedrals'] = kept

        return new_zmat, removed, h_atoms

    def _drop_rule_centre_elements(self):
        """
        Validates drop_improper_h_dihedrals (False/None, True, or element
        symbols; a single symbol may be given as a string).

        :return:
            Whether the rule is on, and the list of centre elements (None: every
            improper centre).
        """

        flag = self.drop_improper_h_dihedrals
        if flag is None or isinstance(flag, (bool, np.bool_)):
            return bool(flag), None
        if isinstance(flag, str):
            # one element symbol: 'Cl' is the element Cl, not the letters C and l
            flag = [flag]
        assert_msg_critical(
            isinstance(flag, (list, tuple, set)) and all(isinstance(e, str) for e in flag),
            'IMForceFieldGenerator: drop_improper_h_dihedrals must be False, True or a list of element '
            f"symbols such as ['N'] (got {self.drop_improper_h_dihedrals!r}).")

        return len(flag) > 0, [str(e).capitalize() for e in flag]

    def define_z_matrix_dict(self, molecule, add_coordinates=None):
        """
        Builds the primitive internal coordinates of a molecule with
        geomeTRIC (bonds from the connectivity, all angles and proper
        dihedrals; no impropers).

        :param molecule:
            The molecule.
        :param add_coordinates:
            Optional dict of extra coordinates, one per key: 'bond', 'angle',
            'dihedral' or 'impropers' (0-based tuple).

        :return:
            The z-matrix dict with 'bonds', 'angles', 'dihedrals' and
            'impropers' (0-based tuples).
        """

        g_molecule = geometric.molecule.Molecule()
        g_molecule.elem = molecule.get_labels()
        g_molecule.xyzs = [molecule.get_coordinates_in_bohr() * geometric.nifty.bohr2ang]
        g_molecule.build_topology()
        g_molecule.build_bonds()

        bonds = [tuple(x) for x in g_molecule.Data['bonds']]
        angles = [tuple(x) for x in g_molecule.find_angles()]
        dihedrals = [tuple(x) for x in g_molecule.find_dihedrals()]

        zmat = {
            "bonds": bonds[:],
            "angles": angles[:],
            "dihedrals": dihedrals[:],
            "impropers": [],
        }

        if add_coordinates is not None:
            for key, coord in add_coordinates.items():
                c = tuple(coord)
                if key == 'bond' and c not in zmat["bonds"]:
                    zmat["bonds"].append(c)
                elif key == 'angle' and c not in zmat["angles"]:
                    zmat["angles"].append(c)
                elif key == 'dihedral' and c not in zmat["dihedrals"]:
                    zmat["dihedrals"].append(c)
                elif key == 'impropers' and c not in zmat["impropers"]:
                    zmat["impropers"].append(c)

        return zmat

    def set_up_the_system(self, molecule, imforcefieldfiles=None, exclude_rot_bonds=None):
        """
        Prepares the construction: the database file names, the z-matrix of
        every root (read from an existing database, otherwise built from the
        molecule with the MM impropers, the reaction coordinates and the
        optional drop_improper_h_dihedrals rule), the symmetry information,
        the interpolation settings (frozen here: later attribute changes are
        not used), and the seed structures (reaction structures and the
        conformers of the MM conformer generator). Single MPI rank only.

        :param molecule:
            The molecule.
        :param imforcefieldfiles:
            Optional dict {root: database file}. Default: im_database_<root>.h5
            in the working directory. The z-matrix is read from a file that
            exists.
        :param exclude_rot_bonds:
            Optional list of 1-based atom pairs removed from the rotatable
            bonds of the MM force field (no conformer scan around them).
        """

        assert_msg_critical(self.nodes == 1,
                            'IMForceFieldGenerator: the database construction supports a single MPI rank only.')

        bond_consistency_checK = True
        if (self.use_eq_bond_length == True and self.use_inverse_bond_length == True):
            bond_consistency_checK = False
        assert_msg_critical(bond_consistency_checK,
                            "The choice of the bond description needs to be unique! eq can not be used with rinv!")
        
        def filter_rotatable_bonds(rotatable_bonds, exclude_rot_bonds):
            """
            Removes user-defined bonds from the rotatable bond list.

            :param rotatable_bonds:
                The 1-based rotatable bond list from MMForceFieldGenerator.

            :return:
                The filtered 1-based rotatable bond list.
            """

            if not exclude_rot_bonds:
                return rotatable_bonds

            excluded_bonds = set()
            for bond in exclude_rot_bonds:
                assert_msg_critical(
                    len(bond) == 2,
                    'IMForceFieldGenerator.filter_rotatable_bonds: '
                    'excluded_rotatable_bonds must contain atom-index pairs.')

                atom_i, atom_j = int(bond[0]), int(bond[1])

                assert_msg_critical(
                    atom_i > 0 and atom_j > 0 and atom_i != atom_j,
                    'IMForceFieldGenerator.filter_rotatable_bonds: '
                    'excluded_rotatable_bonds expects 1-based atom indices.')

                excluded_bonds.add(tuple(sorted((atom_i, atom_j))))

            filtered_rotatable_bonds = []
            excluded_found = set()

            for atom_i, atom_j in rotatable_bonds:
                bond = tuple(sorted((int(atom_i), int(atom_j))))

                if bond in excluded_bonds:
                    excluded_found.add(bond)
                    continue

                filtered_rotatable_bonds.append([int(atom_i), int(atom_j)])

            missing_bonds = excluded_bonds - excluded_found
            if missing_bonds:
                self.ostream.print_warning(
                    'Some user-excluded rotatable bonds were not present in the '
                    f'MM rotatable bond list: {sorted(missing_bonds)}')
                self.ostream.flush()

            return filtered_rotatable_bonds

        def regroup_by_rotatable_connection(molecule, groups, rotatable_bonds, conn):
            """
            Splits symmetry groups (e.g. the H atoms of a CH3 group) into rotor
            groups attached to one rotatable bond. The symmetry groups passed
            here are empty, so the result is empty.

            :return:
                The regrouped atoms per state ('gs', 'es', 'non_rotatable')
                and the rotor groups per state.
            """

            new_groups = {'gs': [], 'es': [], 'non_rotatable': []}
            rot_groups = {'gs': [], 'es': []}
            labels = molecule.get_labels()

            def determine_state(a1, a2):
                neighbors_a1 = sum(conn[a1])
                neighbors_a2 = sum(conn[a2])
                element_a1 = labels[a1][0]
                element_a2 = labels[a2][0]
                # NH2 rule
                if (element_a1 == 'N' and neighbors_a1 == 2) or (element_a2 == 'N' and neighbors_a2 == 2):
                    return 'gs'

                # Oxygen + sp2 carbon → es
                if (element_a1 == 'O' and element_a2 == 'C' and neighbors_a2 == 3) or \
                        (element_a2 == 'O' and element_a1 == 'C' and neighbors_a1 == 3):
                    return 'es'

                # sp2-sp2 → es
                if neighbors_a1 == 3 and neighbors_a2 == 3:
                    return 'es'

                # default → gs
                return 'gs'

            def determine_similar_groups(groups_to_process):
                separated_groups = []
                for group in groups_to_process:
                    # Dictionary to group atoms by their exact neighborhood
                    neighbor_map = {}
                    for atom in group:
                        # Find the indices of all atoms connected to this atom
                        neighbors = tuple(sorted([i for i, is_bonded in enumerate(conn[atom]) if is_bonded]))
                        neighbor_map.setdefault(neighbors, []).append(atom)

                    # Extract the automatically separated groups
                    for subgroup in neighbor_map.values():
                        separated_groups.append(sorted(subgroup))

                return separated_groups

            # Apply the separation logic to fix lumped groups before iterating
            groups = determine_similar_groups(groups)

            for group in groups:
                if len(group) > 3:
                    continue

                connected_subgroups = {}  # key: (state, atom in bond), value: atoms in group connected to it
                connected_rotatable_bonds = set()

                for atom in group:
                    for a1, a2 in rotatable_bonds:
                        state = determine_state(a1, a2)

                        if conn[atom, a1] or conn[atom, a2]:
                            connected_rotatable_bonds.add((a1, a2))

                        if conn[atom, a1]:
                            connected_subgroups.setdefault((state, a1), []).append(atom)
                        if conn[atom, a2]:
                            connected_subgroups.setdefault((state, a2), []).append(atom)

                if len(connected_rotatable_bonds) > 1:
                    new_groups['non_rotatable'].append(group)
                    continue

                if not connected_subgroups:
                    new_groups['non_rotatable'].append(group)
                else:
                    for (state, _), subgroup in connected_subgroups.items():
                        subgroup = list(set(subgroup))
                        if len(subgroup) > 1 and subgroup not in rot_groups[state]:
                            rot_groups[state].append(sorted(subgroup))
                            new_groups[state].append(sorted(subgroup))

            return new_groups, rot_groups

        def _promote_nonrotatable_ring_torsions_to_impropers(zmat, ff_gen, rotatable_bonds_zero_based, impropers=None):
            """
            Sets the impropers of the z-matrix. With impropers given (the only
            use: the MM impropers), the z-matrix gets exactly these impropers
            and its dihedrals are unchanged. Without them (inactive
            alternative), proper ring torsions around non-rotatable bonds are
            replaced by the MM impropers that touch the same bond (matched by
            the middle bond of the torsion, since the full 4-atom set rarely
            matches).

            :return:
                The z-matrix dict.
            """

            if impropers is None:
                rot_set = {frozenset((int(i), int(j))) for i, j in rotatable_bonds_zero_based}

                # Existing impropers from z-matrix (if any)
                impropers = [tuple(int(x) for x in imp) for imp in zmat.get("impropers", [])]
                impropers_seen = set(impropers)

                # Build non-rotatable ring bond -> impropers map from MM impropers.
                impropers_ff = [tuple(int(x) for x in t) for t in ff_gen.impropers.keys()]
                bond_to_impropers = {}
                for imp in impropers_ff:
                    center = int(imp[0])
                    for neigh in imp[1:]:
                        neigh = int(neigh)
                        bond = frozenset((center, neigh))
                        if bond in rot_set:
                            continue
                        if not ff_gen.is_bond_in_ring(center, neigh):
                            continue
                        bond_to_impropers.setdefault(bond, []).append(imp)

                kept_dihedrals = []
                for dih_raw in zmat["dihedrals"]:
                    dih = tuple(int(x) for x in dih_raw)
                    mid = frozenset((dih[1], dih[2]))

                    # Keep proper terms for rotatable or non-ring bonds.
                    if mid in rot_set or ff_gen.is_bond_in_ring(dih[0], dih[1]) and ff_gen.is_bond_in_ring(dih[1], dih[2]) and ff_gen.is_bond_in_ring(dih[2], dih[3]):
                        kept_dihedrals.append(dih)
                        continue

                    # Non-rotatable ring bond: attach any corresponding impropers.
                    promoted = False
                    for imp in bond_to_impropers.get(mid, []):
                        if imp not in impropers_seen:
                            impropers.append(imp)
                            impropers_seen.add(imp)
                        promoted = True

                    # If no improper was found for this bond, keep the dihedral.
                    if not promoted:
                        kept_dihedrals.append(dih)

                zmat["dihedrals"] = kept_dihedrals
                zmat["impropers"] = impropers
            else:
                zmat["impropers"] = impropers
            return zmat

        # define global symmetry information object for gs and es
        # used in all other classes for the construction
        # 1. list of all atoms
        # 2. all rotatable CH3 groups for ground-state
        # 3. all groups that are assigned for rotation in the ground-state
        # 4. all atoms that are not in the symmetry groups (CH3-groups excluded)
        # 5. all CH3 symmetry groups
        # 6. indices 0-based of all rotatable bonds
        # 7. list of indices grouped for all dihedrals attached to a rotatable bond (z_matrix index)
        # 8. actual dihedrals with the indices assigned to each atom
        # 9. start and end of the dihedral section within the z_matrix
        self.symmetry_information = {'gs': (), 'es': ()}

        self.molecule = molecule

        root_extract_z_matrix = None
        if imforcefieldfiles is not None:
            root_extract_z_matrix = {}
            self.sampling_imforcefieldfiles = {}
            self.imforcefieldfiles = {}
            for root in self.roots_to_follow:
                if root not in imforcefieldfiles:
                    self.imforcefieldfiles[root] = f'im_database_{root}.h5'
                else:
                    self.imforcefieldfiles[root] = imforcefieldfiles[root]
                # the z-matrix is read from the database only if the database exists already
                root_extract_z_matrix[root] = Path(self.imforcefieldfiles[root]).exists()

                self.sampling_imforcefieldfiles[root] = f'im_database_sampling_{root}.h5'
        else:
            self.imforcefieldfiles = {}
            self.sampling_imforcefieldfiles = {}
            standard_files: list[str] = [f'im_database_{root}.h5' for root in self.roots_to_follow]
            standard_sampling_files = [f'im_database_sampling_{root}.h5' for root in self.roots_to_follow]
            root_extract_z_matrix = {}
            for root_idx, standard_file in enumerate(standard_files):

                if Path(standard_file).exists():
                    if self.roots_to_follow[root_idx] not in root_extract_z_matrix:
                        root_extract_z_matrix[self.roots_to_follow[root_idx]] = True
                else:
                    root_extract_z_matrix[self.roots_to_follow[root_idx]] = False
                self.imforcefieldfiles[self.roots_to_follow[root_idx]] = standard_file
                self.sampling_imforcefieldfiles[self.roots_to_follow[root_idx]] = standard_sampling_files[root_idx]

        ff_gen = MMForceFieldGenerator()
        ff_gen.ostream.mute()
        ff_gen.partial_charges = molecule.get_partial_charges(molecule.get_charge())
        ff_gen.create_topology(molecule)

        # The AtomMapper class is based on Turtlemap program which is used to
        # determine equivalent atoms within a molecular structure
        symmetry_groups = (list(range(len(molecule.get_labels()))), [], [])
        rotatable_bonds = deepcopy(ff_gen.rotatable_bonds)
        # add an additional filter to remove rotatable bonds around high energy paths
        if exclude_rot_bonds is not None:
            rotatable_bonds = filter_rotatable_bonds(rotatable_bonds, exclude_rot_bonds)
        # Work in zero-based indexing (same convention as z-matrix dihedrals)
        # and remove all symmetry-related rotatable bonds from the scan list.
        rotatable_bonds_zero_based = [tuple(sorted((i - 1, j - 1))) for (i, j) in rotatable_bonds]

        drop_rule, centre_elements = self._drop_rule_centre_elements()

        for root in self.roots_to_follow:
            if root not in self.roots_z_matrix and not root_extract_z_matrix[root]:
                # if no database is provided construct the primitive internal coordinates using geometric
                self.roots_z_matrix[root] = self.define_z_matrix_dict(molecule)
                if self.reaction_structures is not None:
                    merge_info = self.merge_reaction_internal_coordinates(
                        reaction_structures=self.reaction_structures,
                        root=self.roots_to_follow[0],
                        include_existing_root_zmat=True,
                    )
                    if self.use_eq_bond_length:
                        self.eq_bond_length_irc_bonds = merge_info['added_coordinates']['bonds']
                    self.roots_z_matrix[self.roots_to_follow[0]] = merge_info["global_z_matrix"]

                impropers = list(ff_gen.impropers.keys())
                self.roots_z_matrix[root] = _promote_nonrotatable_ring_torsions_to_impropers(self.roots_z_matrix[root], ff_gen, rotatable_bonds_zero_based, impropers=impropers)

                if drop_rule:
                    protected_dihedrals = []
                    if isinstance(self.use_minimized_structures, (list, tuple)) and len(self.use_minimized_structures) > 1:
                        for constraint in (self.use_minimized_structures[1] or []):
                            if not isinstance(constraint, str) and len(constraint) == 4:
                                protected_dihedrals.append([int(x) - 1 for x in constraint])
                    self.roots_z_matrix[root], removed_dihedrals, wag_h_atoms = self._drop_improper_h_dihedrals(
                        self.roots_z_matrix[root], molecule.get_labels(), centre_elements, protected_dihedrals)
                    self.ostream.print_info(
                        'drop_improper_h_dihedrals: removed '
                        f'{len(removed_dihedrals)} proper dihedral(s) terminated by H atom(s) '
                        f'{[a + 1 for a in wag_h_atoms]} (1-based): '
                        f'{[tuple(a + 1 for a in d) for d in removed_dihedrals]}')
                    self.ostream.flush()
            elif root not in self.roots_z_matrix:
                # generate the z-matrix based for the interpolation database provided
                int_driver = InterpolationDriver()
                int_driver.update_settings({
                    'interpolation_type':self.interpolation_type,
                    'weightfunction_type':self.weightfunction_type,
                    'exponent_p':self.exponent_p,
                    'exponent_q':self.exponent_q,
                    'confidence_radius':self.confidence_radius,
                    'imforcefield_file':self.imforcefieldfiles[root],
                    'use_inverse_bond_length':self.use_inverse_bond_length,
                    'use_eq_bond_length':self.use_eq_bond_length,
                    'use_tc_weights':self.use_tc_weights,
                    'tc_weight_mode':self.tc_weight_mode,
                    'use_mass_weight':self.use_mass_weight,
                })

                _, z_matrix = int_driver.read_labels()
                assert_msg_critical(len(z_matrix["bonds"]) != 0,
                                    "The provided database did not proivde a valid z_matrix: number of bonds = 0! Check the correct bond desctiption mode (r, rinv, eq)!")
                self.roots_z_matrix[root] = z_matrix

            dihedral_start = len(self.roots_z_matrix[root]['bonds']) + len(self.roots_z_matrix[root]['angles'])
            dihedral_end = dihedral_start + len(self.roots_z_matrix[root]['dihedrals'])

            self.all_rotatable_bonds = rotatable_bonds_zero_based
            all_exclusion = [element for rot_bond in rotatable_bonds_zero_based for element in rot_bond]

            symmetry_groups_ref = [groups for groups in symmetry_groups[1] if not any(item in all_exclusion for item in groups)]

            # reduce the symmetry to only CH3 or CH2 symmetry groups for the time being
            regrouped, rot_groups = regroup_by_rotatable_connection(molecule, symmetry_groups_ref, rotatable_bonds_zero_based, molecule.get_connectivity_matrix())

            if root == 0:
                non_core_atoms = [element for group in regrouped['gs'] for element in group]
                core_atoms = [element for element in symmetry_groups[0] if element not in non_core_atoms]
                # [6] z-matrix indices of symmetry-related dihedrals and [7] their groups: none here
                self.symmetry_information['gs'] = [symmetry_groups[0], rot_groups['gs'], regrouped['gs'], core_atoms, non_core_atoms, rotatable_bonds_zero_based, [], {}, [], [dihedral_start, dihedral_end]]

            imforcefieldfile = self.imforcefieldfiles[self.roots_to_follow[0]]
            self.states_interpolation_settings[self.roots_to_follow[0]] = {
                'interpolation_type':self.interpolation_type,
                'weightfunction_type':self.weightfunction_type,
                'exponent_p':self.exponent_p,
                'exponent_q':self.exponent_q,
                'confidence_radius':self.confidence_radius,
                'imforcefield_file':imforcefieldfile,
                'use_inverse_bond_length':self.use_inverse_bond_length,
                'use_eq_bond_length': self.use_eq_bond_length,
                'use_tc_weights': self.use_tc_weights,
                'tc_weight_mode': self.tc_weight_mode,
                'use_mass_weight': self.use_mass_weight,
            }
            self.sampling_states_interpolation_settings[self.roots_to_follow[0]] = self.states_interpolation_settings[self.roots_to_follow[0]].copy()
            self.sampling_states_interpolation_settings[self.roots_to_follow[0]]['imforcefield_file'] = self.sampling_imforcefieldfiles[self.roots_to_follow[0]]

        if self.exclude_non_core:
            new_exclusion = {}
            new_inclusion = {}
            for key, entry in self.symmetry_information.items():
                if len(entry) == 0:
                    continue
                if key not in new_exclusion:
                    new_exclusion[key] = None
                    new_inclusion[key] = None

                new_exclusion[key] = [idx for idx, label in enumerate(molecule.get_labels()) if label == 'H' or idx in entry[4]]
                new_inclusion[key] = [idx for idx, _ in enumerate(molecule.get_labels()) if idx not in new_exclusion[key]]

            for new_key in new_exclusion.keys():

                self.symmetry_information[new_key][4] = new_exclusion[new_key]
                self.symmetry_information[new_key][3] = new_inclusion[new_key]

        if self.reaction_structures is not None:
            self.seed_structures = self.build_initial_seed_structures(
                molecule=molecule,
                reaction_structures=self.reaction_structures,
                include_conformers=False,
                reaction_root=self.roots_to_follow[0],
                reaction_key=None,
            )
            self.seed_reaction_entries = list(self.seed_structures[self.roots_to_follow[0]][None])

        if self.add_conformal_structures and 0 in self.roots_to_follow:

            conformers_plus_ts = {0: {}}

            conformer_gen = ConformerGenerator()
            conformer_gen.resp_charges = False

            conformer_results = conformer_gen.generate(molecule)
            dihedral_candidates = list(getattr(conformer_gen, "dihedral_candidates", []) or [])
            conformal_molecules = list(conformer_results.get("molecules", []) or [])

            if len(conformal_molecules) == 0:
                raise RuntimeError('ConformerGenerator returned no conformers.')

            # the unique conformers, stored once for the opt-in conformer-network seed stage
            self.seed_conformers = [self._copy_molecule(conf_mol, molecule) for conf_mol in conformal_molecules]

            rot_bond_set = {
                tuple(sorted((int(a), int(b))))
                for a, b in rotatable_bonds
            }

            selected_dihedrals = []
            for dih_zero, _angle_grid in dihedral_candidates:
                dih_key = tuple(int(x) + 1 for x in dih_zero)

                if tuple(sorted((dih_key[1], dih_key[2]))) not in rot_bond_set:
                    continue

                selected_dihedrals.append(dih_key)

            if self.seed_conformer_network or drop_rule:
                # the scanned dihedrals become the important coordinates of the scan points (network
                # and classic seeds) and must be z-matrix coordinates (drop_improper_h_dihedrals may
                # have removed one); a z-matrix read from an existing database is never changed
                for root in self.roots_to_follow:
                    if root_extract_z_matrix[root]:
                        continue
                    zmat_dihedrals = [tuple(int(x) for x in d) for d in self.roots_z_matrix[root]['dihedrals']]
                    restored = []
                    for dih_key in selected_dihedrals:
                        dih_zero = tuple(int(x) - 1 for x in dih_key)
                        if dih_zero not in zmat_dihedrals and dih_zero[::-1] not in zmat_dihedrals:
                            self.roots_z_matrix[root]['dihedrals'].append(dih_zero)
                            zmat_dihedrals.append(dih_zero)
                            restored.append(dih_key)
                    if len(restored) > 0:
                        dihedral_start = len(self.roots_z_matrix[root]['bonds']) + len(self.roots_z_matrix[root]['angles'])
                        dihedral_end = dihedral_start + len(self.roots_z_matrix[root]['dihedrals'])
                        if root == 0:
                            self.symmetry_information['gs'][-1] = [dihedral_start, dihedral_end]
                        self.ostream.print_info(
                            f'scanned dihedral(s) {restored} (1-based) kept in the z-matrix '
                            '(important coordinates of the scan points).')
                        self.ostream.flush()

            if len(selected_dihedrals) > 0:
                for dih_key in selected_dihedrals:
                    conformers_plus_ts[0][dih_key] = []

                    for conf_mol in conformal_molecules:
                        mol_i = self._copy_molecule(conf_mol, molecule)
                        if (
                            abs(abs(mol_i.get_dihedral_in_degrees(dih_key)) - 180.0) < 5.0
                            or abs(mol_i.get_dihedral_in_degrees(dih_key)) < 3.0
                        ):
                            mol_i.set_dihedral_in_degrees(
                                dih_key,
                                mol_i.get_dihedral_in_degrees(dih_key) + 10.0,
                                verbose=False,
                            )
                        conformers_plus_ts[0][dih_key].append((mol_i, 'normal'))

            else:
                conformers_plus_ts[0][None] = [(self._copy_molecule(conformal_molecules[0], molecule), 'normal')]

            self.seed_structures = conformers_plus_ts

        self._system_is_set_up = True

    def _bootstrap_sampling_db_from_abinitio_db(self, root):
        """
        Builds the xTB mirror of the database for the sampling pre-screen:
        every datapoint of the ab initio database that is not yet in the
        sampling database is recomputed with the xTB triple
        (sampling_driver) at the same geometry and written under the same
        label (ground state only).

        :param root:
            The electronic state.
        """

        # mirror the current interpolation database point by point
        ref_settings = self.states_interpolation_settings[root]
        sampling_settings = self.sampling_states_interpolation_settings[root]

        ref_file = ref_settings['imforcefield_file']
        sampling_file = sampling_settings['imforcefield_file']

        ref_drv = InterpolationDriver(self.roots_z_matrix[root])
        ref_drv.update_settings(ref_settings)
        ref_db_labels, _ = ref_drv.read_labels()

        existing = set()

        # check if database already contains the given datapoint
        if Path(sampling_file).exists():
            sampling_drv = InterpolationDriver(self.roots_z_matrix[root])
            sampling_drv.update_settings(sampling_settings)
            samp_db_labels, _ = sampling_drv.read_labels()
            existing = set(samp_db_labels)

        mol_labels = self.molecule.get_labels()
        sampling_qm, sampling_grad, sampling_hess = self.sampling_driver['gs']

        for label in ref_db_labels:

            if label in existing:
                continue

            ref_dp = InterpolationDatapoint(self.roots_z_matrix[root])
            ref_dp.update_settings(ref_settings)
            ref_dp.read_hdf5(ref_file, label)

            coords = ref_dp.cartesian_coordinates
            mol = Molecule(mol_labels, coords, 'bohr')
            mol.set_charge(self.molecule.get_charge())
            mol.set_multiplicity(self.molecule.get_multiplicity())

            e, _, _ = self._compute_energy(
                sampling_qm,
                mol,
                basis=None,
            )

            g = self._compute_gradient(
                sampling_grad,
                mol,
                basis=None,
                scf_results=None,
                rsp_results=None,
            )

            h = self._compute_hessian(
                sampling_hess,
                mol,
                basis=None,
            )

            # Build datapoint in sampling DB using same label and geometry metadata
            inv_sqrt = None
            mw_grad = g[0].reshape(-1)
            mw_hess = h[0].reshape(g[0].reshape(-1).size, g[0].reshape(-1).size)
            if sampling_settings["use_mass_weight"]:
                masses = mol.get_masses().copy()
                inv_sqrt = 1.0 / np.sqrt(np.repeat(masses, 3))
                grad_vec = g[0].reshape(-1)
                hess_mat = h[0].reshape(grad_vec.size, grad_vec.size)
                mw_grad = inv_sqrt * grad_vec
                mw_hess = (inv_sqrt[:, None] * hess_mat) * inv_sqrt[None, :]

            samp_dp = InterpolationDatapoint(self.roots_z_matrix[root])
            samp_dp.update_settings(sampling_settings)
            samp_dp.cartesian_coordinates = ref_dp.cartesian_coordinates
            samp_dp.eq_bond_lengths = ref_dp.eq_bond_lengths
            samp_dp.imp_int_coordinates = getattr(ref_dp, 'imp_int_coordinates', [])
            samp_dp.inv_sqrt_masses = inv_sqrt
            samp_dp.energy = e[0]
            samp_dp.gradient = mw_grad.reshape(g[0].shape)
            samp_dp.hessian = mw_hess.reshape(h[0].shape)
            samp_dp.confidence_radius = getattr(ref_dp, 'confidence_radius', self.confidence_radius)
            samp_dp.transform_gradient_and_hessian()

            samp_dp.point_label = ref_dp.point_label

            samp_dp.write_hdf5(sampling_file, label)

    def _write_database_metadata(self, datafile, event, states_basis=None):
        """
        Appends the construction settings to the history stored in the HDF5
        file attribute 'msi_construction_history' of the database (a JSON
        list of dicts; the first entry and the last 19 are kept). The file
        stores no other settings, and attributes are not keys of the file,
        so read_labels and read_hdf5 are unaffected. A failure is reported
        as a warning and never stops the construction.

        :param datafile:
            The database file (must exist).
        :param event:
            A short description of the stage, e.g. 'seed stage finished'.
        :param states_basis:
            The basis-set labels dict, or None.
        """

        root = self.roots_to_follow[0]
        settings = dict(self.states_interpolation_settings.get(root) or {})
        settings.pop('imforcefield_file', None)
        if settings.get('use_inverse_bond_length'):
            bond_mode = 'inverse'
        elif settings.get('use_eq_bond_length'):
            bond_mode = 'eq'
        else:
            bond_mode = 'plain'

        qm_driver = self.drivers['gs'][0] if self.drivers.get('gs') is not None else None
        if isinstance(qm_driver, XtbDriver):
            qm_method = getattr(qm_driver, 'xtb_method', None)
        else:
            qm_method = getattr(qm_driver, 'xcfun', None)
            if qm_method is not None and not isinstance(qm_method, str):
                qm_method = getattr(qm_method, 'get_func_label', lambda: str(qm_method))()

        seed_network_settings = self.seed_network_settings
        if isinstance(seed_network_settings, dict):
            seed_network_settings = dict(seed_network_settings)

        entry = {
            'event': event,
            'time': time.strftime('%Y-%m-%dT%H:%M:%S%z'),
            'veloxchem_version': __version__,
            'bond_mode': bond_mode,
            'interpolation_settings': settings,
            'qm_driver': type(qm_driver).__name__ if qm_driver is not None else None,
            'qm_method': qm_method,
            'basis_set': (states_basis or {}).get('gs', self.gs_basis_set_label),
            'energy_threshold_kcal_per_atom': self.energy_threshold,
            'gradient_rmsd_thrsh': self.gradient_rmsd_thrsh,
            'use_opt_confidence_radius': self.use_opt_confidence_radius,
            'gate_distance_prefilter': self.gate_distance_prefilter,
            'trust_radius_fit_cache': self.trust_radius_fit_cache,
            'identify_relevant_int_coordinates': self.identify_relevant_int_coordinates,
            'use_minimized_structures': self.use_minimized_structures,
            'exclude_non_core': self.exclude_non_core,
            'temperature_K': self.temperature,
            'ensemble': self.ensemble,
            'timestep_fs': self.timestep,
            'nsteps': self.nsteps,
            'sampling_settings': self.sampling_settings,
            'metadynamics_enabled': bool((self.metadynamics_settings or {}).get('enabled', False)),
            'seed_conformer_network': self.seed_conformer_network,
            'seed_network_settings': seed_network_settings,
            'seed_duplicate_rmsd_thrsh': self.seed_duplicate_rmsd_thrsh,
            'seed_scan_clash_relief': self.seed_scan_clash_relief,
            'drop_improper_h_dihedrals': self.drop_improper_h_dihedrals,
        }

        try:
            with h5py.File(datafile, 'a') as h5f:
                history = []
                if 'msi_construction_history' in h5f.attrs:
                    raw = h5f.attrs['msi_construction_history']
                    history = json.loads(raw.decode() if isinstance(raw, bytes) else str(raw))
                history.append(entry)
                if len(history) > 20:
                    history = history[:1] + history[-19:]
                h5f.attrs['msi_construction_history'] = json.dumps(history, default=self._json_default)
        except (OSError, ValueError, TypeError, KeyError) as exc:
            self.ostream.print_warning(f'IMForceFieldGenerator: the construction settings could not be written '
                                       f'to {datafile} ({type(exc).__name__}: {exc}).')
            self.ostream.flush()

    def compute(self, molecule, states_basis=None):

        """
        Constructs the interpolation database. Steps: the settings check
        (unused settings are reported), set_up_the_system if it was not
        called; for a new database the seed stage, which is the network seed
        stage (seed_conformer_network), the classic conformer and reaction
        seeds, or the (optimised) input structure, while for an existing
        database the seed stage is skipped and the dynamics starts from the
        input geometry; the xTB mirror (sampling_settings['enabled']); the
        collector dynamics, which adds datapoints until desired_point_density
        datapoints exist or converged_cycle cycles add none; and the final
        quality check. Single MPI rank and ground state only.

        :param molecule:
            The molecule (the start structure of the dynamics on a restart).
        :param states_basis:
            Optional dict of basis-set labels {'gs': ..., 'es': ...}; missing
            keys are taken from gs_basis_set_label and es_basis_set_label.

        :return:
            The im_results dict; 'n_datapoints' holds the number of
            datapoints per root.
        """

        assert_msg_critical(self.nodes == 1,
                            'IMForceFieldGenerator: the database construction supports a single MPI rank only.')

        self._check_user_settings()

        # the basis-set labels per state: the argument, completed by gs_basis_set_label and es_basis_set_label
        default_basis = {'gs': self.gs_basis_set_label, 'es': self.es_basis_set_label}
        if states_basis is None:
            states_basis = default_basis
        else:
            assert_msg_critical(isinstance(states_basis, dict),
                                "IMForceFieldGenerator.compute: states_basis must be a dict such as "
                                "{'gs': 'def2-svp'}.")
            states_basis = {**default_basis, **states_basis}

        # First set up the system for which the database needs to be constructed
        if not self._system_is_set_up:
            self.ostream.print_warning('System is not set up. Automatic set up ')
            self.set_up_the_system(molecule)

        self.ostream.print_blank()
        self.ostream.print_header('IM Database Construction')
        self.ostream.print_header('========================')
        self.ostream.print_blank()
        self.ostream.print_info(
            f'IM force field files: {self.imforcefieldfiles}'
        )
        self.ostream.print_info('System setup completed.')

        if len(self.roots_to_follow) > 1:
            self.ostream.print_warning(
                'Multi-state database construction is not supported in this version.'
            )
        else:
            self.ostream.print_warning(
                'Single-state collection follows one PES only; ensure state separation is sufficient.'
            )

            self.ostream.print_blank()
            self.ostream.flush()

            imforcefieldfile = self.imforcefieldfiles[self.roots_to_follow[0]]

            self.dynamics_settings = {
                'drivers': self.drivers,
                'basis_set_label': states_basis,
                'duration': self.duration, 'temperature': self.temperature, 'solvent': self.solvent,
                'pressure': self.pressure, 'force_constant': self.force_constant, 'ensemble': self.ensemble,
                'timestep': self.timestep, 'nsteps': self.nsteps, 'friction': self.friction,
                'snapshots': self.snapshots, 'trajectory_file': self.trajectory_file, 'reference_struc_energy_file': self.reference_struc_energy_file,
                'desired_datapoint_density': self.desired_point_density, 'converged_cycle': self.converged_cycle,
                'energy_threshold': self.energy_threshold, 'grad_rmsd_thrsh': self.gradient_rmsd_thrsh,
                'load_system': self.load_system, 'start_collect': self.start_collect, 'roots_to_follow': self.roots_to_follow,
                'sampling_drivers': self.sampling_driver, 'sampling_settings': self.sampling_settings,
                'metadynamics': self.metadynamics_settings,
            }

            molecules_to_add_info = []

            database_exists = Path(imforcefieldfile).exists()
            assert_msg_critical(
                database_exists or not self.seed_conformer_network or self.seed_structures is not None,
                'IMForceFieldGenerator: seed_conformer_network needs seed structures; set '
                'add_conformal_structures = True or give reaction_structures.')
            if database_exists:
                self._write_database_metadata(imforcefieldfile, 'construction restarted', states_basis)

            if self.seed_structures is not None and not database_exists and self.seed_conformer_network:

                molecules_to_add_info = self._seed_conformer_network(molecule, states_basis)

            elif self.seed_structures is not None and not database_exists:

                molecules_to_add_info = []
                # optimised conformers already added as datapoints (all rotatable bonds)
                added_normal_minima = []
                for counter, entry in enumerate(self.seed_structures.items()):

                    key, molecules_info = entry

                    for i, mol_entries in enumerate(molecules_info.items()):
                        current_molecule_to_add_info = []
                        dih_key, mol_info = mol_entries

                        optimized_normal_conformers = []
                        expanded_scan_from_optimized_minima = False
                        seed_index = 0

                        while seed_index < len(mol_info):
                            seed_entry = mol_info[seed_index]
                            seed_index += 1

                            if len(seed_entry) == 2:
                                mol, mode = seed_entry
                                seed_constraints = []
                            else:
                                mol, mode, seed_constraints = seed_entry[0], seed_entry[1], list(seed_entry[2] or [])

                            optimized_molecule_for_scan = None

                            if self.use_minimized_structures[0]:
                                transition = False
                                constraints_global = []
                                if len(self.use_minimized_structures[1]) > 0:
                                    constraints_global.extend(self.use_minimized_structures[1])  # extend, not append

                                if mode == "transition":
                                    transition = True
                                    constraints_global = []

                                elif mode == "constraint" and dih_key is not None:
                                    if dih_key not in constraints_global:
                                        constraints_global.append(dih_key)
                                else:
                                    constraints_global = list(seed_constraints or [])

                                optimized_molecule = None
                                opt_results = None
                                scf_results = None
                                if self.roots_to_follow[0] == 0 and isinstance(self.drivers['gs'][0], ScfRestrictedDriver) or self.roots_to_follow[0] == 0 and isinstance(self.drivers['gs'][0], ScfUnrestrictedDriver):

                                    current_basis = MolecularBasis.read(mol, states_basis['gs'])
                                    _, scf_results, _ = self._compute_energy(self.drivers['gs'][0], mol, current_basis)

                                    opt_results = self._run_optimization(
                                        self.drivers['gs'][0],
                                        mol,
                                        constraints=constraints_global,
                                        transition=transition,
                                        index_offset=0,
                                        compute_args=(current_basis, scf_results),
                                    )
                                    optimized_molecule = opt_results['final_molecule']
                                    optimized_molecule_for_scan = optimized_molecule
                                    current_basis = MolecularBasis.read(optimized_molecule, states_basis['gs'])
                                    current_molecule_to_add_info.append((optimized_molecule, current_basis, self.roots_to_follow, constraints_global,transition))
                                    self.ostream.print_blank()
                                    self.ostream.print_header('Optimized Molecule')
                                    self.ostream.print_header('------------------')
                                    self.ostream.print_block(optimized_molecule.get_xyz_string())
                                    self.ostream.print_blank()
                                    self.ostream.flush()

                                elif self.roots_to_follow[0] == 0 and isinstance(self.drivers['gs'][0], XtbDriver):

                                    opt_results = self._run_optimization(
                                        self.drivers['gs'][0],
                                        mol,
                                        constraints=constraints_global,
                                        transition=transition,
                                        index_offset=0,
                                    )
                                    optimized_molecule = opt_results['final_molecule']
                                    optimized_molecule_for_scan = optimized_molecule
                                    current_basis = MolecularBasis.read(optimized_molecule, states_basis['gs'])
                                    current_molecule_to_add_info.append((optimized_molecule, current_basis, self.roots_to_follow, constraints_global,transition))
                                    self.ostream.print_blank()
                                    self.ostream.print_header('Optimized Molecule')
                                    self.ostream.print_header('------------------')
                                    self.ostream.print_block(optimized_molecule.get_xyz_string())
                                    self.ostream.print_blank()
                                    self.ostream.flush()

                            else:
                                if self.roots_to_follow[0] == 0:
                                    current_basis = MolecularBasis.read(mol, states_basis['gs'])
                                    current_molecule_to_add_info.append((mol, current_basis, self.roots_to_follow, []))
                                else:
                                    current_basis = MolecularBasis.read(mol, states_basis['es'])
                                    current_molecule_to_add_info.append((mol, current_basis, self.roots_to_follow, []))

                            if mode == "normal" and optimized_molecule_for_scan is not None:
                                # a conformer that optimised into a minimum found before is not
                                # added again; it stays a scan minimum only if that minimum came
                                # from another rotatable bond
                                if self._is_same_minimum(optimized_molecule_for_scan, optimized_normal_conformers):
                                    current_molecule_to_add_info.pop()
                                    self.ostream.print_info(f"Conformer {seed_index} of {dih_key} is a known minimum (same key).")
                                else:
                                    if self._is_same_minimum(optimized_molecule_for_scan, added_normal_minima):
                                        current_molecule_to_add_info.pop()
                                        self.ostream.print_info(f"Conformer {seed_index} of {dih_key} is a known minimum (earlier key).")
                                    else:
                                        added_normal_minima.append(optimized_molecule_for_scan)
                                    optimized_normal_conformers.append(optimized_molecule_for_scan)
                                self.ostream.flush()

                            if (
                                not expanded_scan_from_optimized_minima
                                and seed_index == len(mol_info)
                                and dih_key is not None
                                and len(optimized_normal_conformers) > 1
                            ):
                                expanded_scan_from_optimized_minima = True

                                def _w360(a):
                                    return float(a) % 360.0

                                def _k(a):
                                    key = round(_w360(a), 6)
                                    return 0.0 if key >= 360.0 else key

                                def _signed_deg(a):
                                    a = _w360(a)
                                    return a - 360.0 if a > 180.0 else a

                                minima_by_angle = {}
                                for opt_mol in optimized_normal_conformers:
                                    ang = _k(opt_mol.get_dihedral_in_degrees(dih_key))
                                    minima_by_angle.setdefault(ang, opt_mol)

                                minima = sorted(minima_by_angle.items())

                                if len(minima) > 1:
                                    scan_entries = []

                                    for idx, (cur_ang, cur_mol) in enumerate(minima):
                                        nxt_ang, nxt_mol = minima[(idx + 1) % len(minima)]
                                        step = (nxt_ang - cur_ang) % 360.0

                                        if step < 1.0e-8:
                                            continue

                                        for frac, scan_mode in (
                                            (0.25, 'constraint'),
                                            (0.50, 'transition'),
                                            (0.75, 'constraint'),
                                        ):
                                            target_deg = _w360(cur_ang + frac * step)
                                            base_mol = cur_mol if frac <= 0.50 else nxt_mol
                                            scan_entries.append((target_deg, scan_mode, base_mol))
                                    self.ostream.print_blank()
                                    self.ostream.print_line(
                                        f"optimized-minima scan schedule for {dih_key}: "
                                        f"{[(round(_signed_deg(a), 1), m) for a, m, _ in scan_entries]}")
                                    self.ostream.flush()

                                    for target_deg, scan_mode, base_mol in scan_entries:
                                        mol_i = self._copy_molecule(base_mol, molecule)
                                        mol_i.set_dihedral_in_degrees(dih_key, target_deg, verbose=False)

                                        actual = float(mol_i.get_dihedral_in_degrees(dih_key)) % 360.0
                                        if self._circular_diff_deg(actual, target_deg) > 2.0 and self.seed_scan_clash_relief:
                                            relieved_mol, aux_changes = self._place_scan_seed(
                                                [base_mol], dih_key, target_deg)
                                            if relieved_mol is not None:
                                                mol_i = relieved_mol
                                                actual = float(mol_i.get_dihedral_in_degrees(dih_key)) % 360.0
                                                self.ostream.print_info(
                                                    f"{scan_mode} seed {dih_key} = {target_deg:.2f} deg placed with "
                                                    f"clash relief: {[(d, round(x, 1)) for d, x in aux_changes]}")
                                                self.ostream.flush()

                                                # relax the placed geometry with the dihedral frozen,
                                                # so that the TS search starts from a relaxed structure
                                                if scan_mode == 'transition':
                                                    driver = self.drivers['gs'][0]
                                                    compute_args = None
                                                    if isinstance(driver, (ScfRestrictedDriver, ScfUnrestrictedDriver)):
                                                        relax_basis = MolecularBasis.read(mol_i, states_basis['gs'])
                                                        _, relax_scf_results, _ = self._compute_energy(driver, mol_i, relax_basis)
                                                        compute_args = (relax_basis, relax_scf_results)
                                                    mol_i = self._run_optimization(
                                                        driver, mol_i, constraints=[dih_key], index_offset=0,
                                                        compute_args=compute_args)['final_molecule']
                                                    actual = float(mol_i.get_dihedral_in_degrees(dih_key)) % 360.0

                                        if self._circular_diff_deg(actual, target_deg) > 2.0:
                                            raise RuntimeError(
                                                f"Failed to set dihedral {dih_key} to {target_deg:.2f} deg "
                                                f"(actual {actual:.2f} deg)."
                                            )

                                        mol_info.append((mol_i, scan_mode))

                        if len(molecules_to_add_info) == 0:

                            molecules_to_add_info.append(current_molecule_to_add_info[0])

                        self.add_point(current_molecule_to_add_info, self.states_interpolation_settings)

            elif self.seed_structures is None and not database_exists:

                molecules_to_add_info = []
                if self.use_minimized_structures[0]:
                    optimized_molecule = None

                    opt_results = None
                    scf_results = None

                    if self.roots_to_follow[0] == 0 and isinstance(self.drivers['gs'][0], ScfRestrictedDriver) or self.roots_to_follow[0] == 0 and isinstance(self.drivers['gs'][0], ScfUnrestrictedDriver):

                        current_basis = MolecularBasis.read(molecule, states_basis['gs'])
                        _,scf_results, _ = self._compute_energy(self.drivers['gs'][0], molecule, current_basis)

                        opt_results = self._run_optimization(
                            self.drivers['gs'][0],
                            molecule,
                            constraints=self.use_minimized_structures[1],
                            index_offset=0,
                            compute_args=(current_basis, scf_results),
                        )

                        optimized_molecule = opt_results['final_molecule']

                        current_basis = MolecularBasis.read(optimized_molecule, states_basis['gs'])
                        molecules_to_add_info.append((optimized_molecule, current_basis, self.roots_to_follow, self.use_minimized_structures[1]))

                        self.ostream.print_blank()
                        self.ostream.print_header('Optimized Molecule')
                        self.ostream.print_header('------------------')
                        self.ostream.print_block(optimized_molecule.get_xyz_string())
                        self.ostream.print_blank()
                        self.ostream.flush()

                    elif self.roots_to_follow[0] == 0 and isinstance(self.drivers['gs'][0], XtbDriver):

                        opt_results = self._run_optimization(
                            self.drivers['gs'][0],
                            molecule,
                            constraints=self.use_minimized_structures[1],
                            index_offset=0,
                        )

                        optimized_molecule = opt_results['final_molecule']

                        current_basis = MolecularBasis.read(optimized_molecule, states_basis['gs'])
                        molecules_to_add_info.append((optimized_molecule, current_basis, self.roots_to_follow, self.use_minimized_structures[1]))

                        self.ostream.print_blank()
                        self.ostream.print_header('Optimized Molecule')
                        self.ostream.print_header('------------------')
                        self.ostream.print_block(optimized_molecule.get_xyz_string())
                        self.ostream.print_blank()
                        self.ostream.flush()

                    self.add_point(molecules_to_add_info, self.states_interpolation_settings)

                else:
                    if self.roots_to_follow[0] == 0:
                        current_basis = MolecularBasis.read(molecule, states_basis['gs'])
                        molecules_to_add_info.append((molecule, current_basis, self.roots_to_follow, []))
                    else:
                        current_basis = MolecularBasis.read(molecule, states_basis['es'])
                        molecules_to_add_info.append((molecule, current_basis, self.roots_to_follow, []))

                    self.add_point(molecules_to_add_info, self.states_interpolation_settings)

            else:

                # restart: the seed stage only runs for a new database
                self.ostream.print_warning(
                    f'The database {imforcefieldfile} exists: the seed stage (conformer network, conformer '
                    'and reaction seeds, optimisation of the input structure) is skipped and the dynamics '
                    'starts from the input geometry.')
                if self.seed_conformer_network:
                    self.ostream.print_warning(
                        'seed_conformer_network: if an earlier run stopped during the network seed stage, '
                        'the database holds only part of the seeds; remove it to rebuild them.')
                self.ostream.flush()

                if self.roots_to_follow[0] == 0:
                    current_basis = MolecularBasis.read(molecule, states_basis['gs'])
                    molecules_to_add_info.append((molecule, current_basis, self.roots_to_follow, []))
                else:
                    current_basis = MolecularBasis.read(molecule, states_basis['es'])
                    molecules_to_add_info.append((molecule, current_basis, self.roots_to_follow, []))

            if not database_exists and Path(imforcefieldfile).exists():
                self._write_database_metadata(imforcefieldfile, 'seed stage finished', states_basis)

            density_of_datapoints = self.determine_datapoint_density(self.states_interpolation_settings)

            self.states_data_point_density = density_of_datapoints

            if self.sampling_settings.get('enabled', False):
                self._bootstrap_sampling_db_from_abinitio_db(self.roots_to_follow[0])

            dynamics_molecule = molecules_to_add_info[self.roots_to_follow[0]][0]
            forcefield_generator = MMForceFieldGenerator(ostream=self.ostream)
            self.dynamics_settings['trajectory_file'] = f'trajectory_{self.roots_to_follow[0]}.pdb'
            forcefield_generator.partial_charges = dynamics_molecule.get_partial_charges(dynamics_molecule.get_charge())

            forcefield_generator.create_topology(dynamics_molecule)
            im_database_driver = IMDatabasePointCollecter(ostream=self.ostream)
            im_database_driver.non_core_symmetry_groups = self.symmetry_information
            im_database_driver.platform = self.open_mm_platform
            im_database_driver.consider_locality = self.consider_locality

            # set optimization features in the construction run
            im_database_driver.identify_relevant_int_coordinates = (self.identify_relevant_int_coordinates, self.use_minimized_structures[1])
            im_database_driver.use_opt_confidence_radius = self.use_opt_confidence_radius
            im_database_driver.gate_distance_prefilter = self.gate_distance_prefilter
            im_database_driver.trust_radius_fit_cache = self.trust_radius_fit_cache

            im_database_driver.system_from_molecule(dynamics_molecule, self.roots_z_matrix, forcefield_generator, solvent=self.solvent, qm_atoms='all')

            density_of_datapoints = self.determine_datapoint_density(self.states_interpolation_settings)
            desired_point_density = int(self.dynamics_settings['desired_datapoint_density'])
            reached_target_density = False
            current_structure_density = {}
            for root in density_of_datapoints.keys():
                value = density_of_datapoints[root]
                current_structure_density[root] = value
                if value >= desired_point_density:
                    reached_target_density = True
            if not reached_target_density:
                im_database_driver.density_around_data_point = current_structure_density

                im_database_driver.update_settings(self.dynamics_settings, self.states_interpolation_settings, self.sampling_states_interpolation_settings)

                im_database_driver.run_qmmm()

                # individual impes run objects
                self.qm_energies.append(im_database_driver.qm_potentials)
                self.total_energies.append(im_database_driver.total_energies)
                self.kinetic_energies.append(im_database_driver.kinetic_energies)
                self.state_specific_molecules = im_database_driver.state_specific_molecules
                self.unique_molecules.append(im_database_driver.allowed_molecules)

                self._confirm_database_quality(molecule, basis=states_basis, im_settings=self.states_interpolation_settings, given_molecular_structures=self.state_specific_molecules)

            density_of_datapoints = self.determine_datapoint_density(self.states_interpolation_settings)
            self.states_data_point_density = density_of_datapoints

            self.ostream.print_blank()
            self.ostream.print_header('Successfully constructed the interpolation database.')
            self.ostream.print_header('------------------')
            self.ostream.print_block(str(self.states_data_point_density))
            self.ostream.print_blank()
            self.ostream.flush()

            self.im_results['n_datapoints'] = self.states_data_point_density

        return self.im_results

    def merge_reaction_internal_coordinates(
        self,
        reaction_structures,
        root=0,
        include_existing_root_zmat=True,
        reference_molecule=None,
        forced_coordinates=None,
        enforce_same_atoms=True,
    ):
        """
        Builds one global z-matrix for a reaction path by merging all unique
        internal coordinates of the reaction structures (bonds as unordered
        pairs, angles and dihedrals in either direction, impropers by centre).

        :param reaction_structures:
            The reaction-path structures: Molecule, XYZ string, a tuple or
            list whose first item is one of these, or a dict with the key
            'molecule'.
        :param root:
            The root whose roots_z_matrix entry is the starting point.
        :param include_existing_root_zmat:
            If True and the root has a z-matrix, start from it and append the
            missing coordinates of the reaction structures.
        :param reference_molecule:
            The reference for charge, multiplicity and the atom checks
            (default: self.molecule, else the first reaction structure).
        :param forced_coordinates:
            Coordinates that must be in the result, as a dict with the keys
            bond(s), angle(s), dihedral(s), improper(s) and 0-based tuples or
            lists of tuples, e.g. {'bonds': [(donor, H), (acceptor, H)]}.
        :param enforce_same_atoms:
            If True, all structures must have the atom count and the atom
            labels (in order) of the reference.

        :return:
            A dict with 'global_z_matrix', 'per_structure_z_matrices',
            'added_coordinates' (the coordinates appended to the starting
            z-matrix) and 'added_counts'.
        """

        assert_msg_critical(
            isinstance(reaction_structures, (list, tuple)) and len(reaction_structures) > 0,
            "merge_reaction_internal_coordinates: reaction_structures must be a non-empty list/tuple.",
        )

        def _extract_molecule(entry):
            mol_like = None
            if isinstance(entry, Molecule):
                mol_like = entry
            elif isinstance(entry, str):
                mol_like = Molecule.from_xyz_string(entry)
            elif isinstance(entry, dict):
                mol_like = entry.get("molecule", None)
            elif isinstance(entry, (list, tuple)) and len(entry) > 0:
                mol_like = entry[0]

            assert_msg_critical(
                mol_like is not None,
                "merge_reaction_internal_coordinates: failed to extract molecule from one reaction entry.",
            )

            if isinstance(mol_like, Molecule):
                mol = Molecule.from_xyz_string(mol_like.get_xyz_string())
            elif isinstance(mol_like, str):
                mol = Molecule.from_xyz_string(mol_like)
            else:
                assert_msg_critical(
                    False,
                    "merge_reaction_internal_coordinates: molecule entry must be Molecule or XYZ string.",
                )
            return mol

        molecules = [_extract_molecule(entry) for entry in reaction_structures]

        if reference_molecule is None:
            reference_molecule = getattr(self, "molecule", None)
            if reference_molecule is None:
                reference_molecule = molecules[0]

        ref_labels = list(reference_molecule.get_labels())
        ref_n_atoms = len(ref_labels)
        ref_charge = int(reference_molecule.get_charge())
        ref_mult = int(reference_molecule.get_multiplicity())

        for i, mol in enumerate(molecules):
            mol.set_charge(ref_charge)
            mol.set_multiplicity(ref_mult)
            if enforce_same_atoms:
                cur_labels = list(mol.get_labels())
                assert_msg_critical(
                    len(cur_labels) == ref_n_atoms,
                    f"merge_reaction_internal_coordinates: structure {i} has different atom count.",
                )
                assert_msg_critical(
                    cur_labels == ref_labels,
                    f"merge_reaction_internal_coordinates: structure {i} has different atom ordering/labels.",
                )

        sections = ("bonds", "angles", "dihedrals", "impropers")

        def _coord_key(section, coord):
            c = tuple(int(x) for x in coord)
            if section == "bonds":
                return tuple(sorted(c))
            if section == "angles":
                rev = (c[2], c[1], c[0])
                return c if c <= rev else rev
            if section == "dihedrals":
                rev = c[::-1]
                return c if c <= rev else rev
            if section == "impropers":
                if len(c) == 4:
                    return (c[0],) + tuple(sorted(c[1:]))
            return c

        def _normalize_forced_coords(coords_dict):
            out = {k: [] for k in sections}
            if coords_dict is None:
                return out

            key_map = {
                "bond": "bonds", "bonds": "bonds",
                "angle": "angles", "angles": "angles",
                "dihedral": "dihedrals", "dihedrals": "dihedrals",
                "improper": "impropers", "impropers": "impropers",
            }

            for in_key, values in coords_dict.items():
                assert_msg_critical(
                    in_key in key_map,
                    f"merge_reaction_internal_coordinates: invalid forced coordinate key '{in_key}'.",
                )
                dst_key = key_map[in_key]

                if values is None:
                    continue

                if isinstance(values, tuple) and len(values) in (2, 3, 4):
                    values = [values]

                assert_msg_critical(
                    isinstance(values, (list, tuple)),
                    f"merge_reaction_internal_coordinates: forced '{in_key}' must be tuple or list of tuples.",
                )

                for coord in values:
                    assert_msg_critical(
                        isinstance(coord, (list, tuple)),
                        "merge_reaction_internal_coordinates: forced coordinate must be list/tuple.",
                    )
                    out[dst_key].append(tuple(int(x) for x in coord))

            return out

        if include_existing_root_zmat and root in self.roots_z_matrix and self.roots_z_matrix[root] is not None:
            global_z_matrix = {
                key: [tuple(int(x) for x in coord) for coord in self.roots_z_matrix[root].get(key, [])]
                for key in sections
            }
        else:
            global_z_matrix = self.define_z_matrix_dict(reference_molecule)
            for key in sections:
                global_z_matrix.setdefault(key, [])

        seen = {key: set() for key in sections}
        for key in sections:
            for coord in global_z_matrix[key]:
                seen[key].add(_coord_key(key, coord))

        added_coordinates = {key: [] for key in sections}

        def _append_unique(section, coord):
            coord_t = tuple(int(x) for x in coord)
            expected_len = {"bonds": 2, "angles": 3, "dihedrals": 4, "impropers": 4}[section]
            assert_msg_critical(
                len(coord_t) == expected_len,
                f"merge_reaction_internal_coordinates: invalid coordinate length for {section}: {coord_t}",
            )

            key = _coord_key(section, coord_t)
            if key in seen[section]:
                return False

            global_z_matrix[section].append(coord_t)
            seen[section].add(key)
            added_coordinates[section].append(coord_t)
            return True

        per_structure_z_matrices = []
        for mol in molecules:
            zmat_i = self.define_z_matrix_dict(mol)
            for key in sections:
                zmat_i.setdefault(key, [])
            per_structure_z_matrices.append(zmat_i)

            for key in sections:
                for coord in zmat_i[key]:
                    _append_unique(key, coord)

        forced_coords = _normalize_forced_coords(forced_coordinates)
        for key in sections:
            for coord in forced_coords[key]:
                _append_unique(key, coord)

        added_counts = {key: len(vals) for key, vals in added_coordinates.items()}

        return {
            "global_z_matrix": global_z_matrix,
            "per_structure_z_matrices": per_structure_z_matrices,
            "added_coordinates": added_coordinates,
            "added_counts": added_counts,
        }

    def build_initial_seed_structures(
        self,
        molecule,
        reaction_structures=None,
        include_conformers=False,
        reaction_root=0,
        reaction_key=None,
    ):
        """
        Returns the seed structures in the layout of the conformer seeds,
        {state: {dih_key: [(molecule, mode, constraints), ...]}}, with mode
        'normal', 'transition' or 'constraint' and constraints as 1-based
        atom-index tuples or geomeTRIC strings (a transition entry has none).

        :param molecule:
            The molecule (charge and multiplicity of every copy).
        :param reaction_structures:
            The reaction-path entries (Molecule, XYZ string, [molecule, mode,
            constraints] or a dict), or None.
        :param include_conformers:
            True to carry over the conformer seeds of self.seed_structures
            (a 'constraint' entry gets its dihedral as constraint).
        :param reaction_root:
            The root of the reaction seeds.
        :param reaction_key:
            The dih_key under which the reaction seeds are stored.

        :return:
            The seed-structure dict.
        """

        assert_msg_critical(
            isinstance(molecule, Molecule),
            "build_initial_seed_structures: molecule must be a Molecule.",
        )

        roots = [int(r) for r in self.roots_to_follow]
        reaction_root = int(reaction_root)

        assert_msg_critical(
            reaction_root in roots,
            f"build_initial_seed_structures: reaction_root={reaction_root} not in roots_to_follow={roots}.",
        )

        def _copy_norm_molecule(mol_like):
            assert_msg_critical(isinstance(mol_like, (Molecule, str)),
                                "build_initial_seed_structures: expected Molecule or XYZ string.")
            return self._copy_molecule(mol_like, molecule)

        def _normalize_mode(mode):
            m = "normal" if mode is None else str(mode).strip().lower()
            assert_msg_critical(
                m in ("normal", "transition", "constraint"),
                f"build_initial_seed_structures: invalid mode '{mode}'.",
            )
            return m

        def _normalize_constraints(raw_constraints):
            if raw_constraints is None:
                return []
            assert_msg_critical(
                isinstance(raw_constraints, (list, tuple)),
                "build_initial_seed_structures: constraints must be list/tuple or None.",
            )
            out = []
            for c in raw_constraints:
                if isinstance(c, str):
                    out.append(c)
                else:
                    assert_msg_critical(
                        isinstance(c, (list, tuple)) and len(c) in (2, 3, 4),
                        "build_initial_seed_structures: tuple constraints must have length 2/3/4.",
                    )
                    out.append(tuple(int(x) for x in c))
            return out

        seed_structures = {root: {} for root in roots}

        # Optional carry-over from existing conformer-style seeds.
        if include_conformers and isinstance(self.seed_structures, dict):
            for state, dih_map in self.seed_structures.items():
                state_i = int(state)
                if state_i not in seed_structures or not isinstance(dih_map, dict):
                    continue

                for dih_key, entries in dih_map.items():
                    seed_structures[state_i].setdefault(dih_key, [])
                    for entry in entries:
                        if not isinstance(entry, (list, tuple)):
                            continue
                        if len(entry) == 2:
                            mol_obj, mode = entry
                            constraints = []
                        elif len(entry) >= 3:
                            mol_obj, mode, constraints = entry[0], entry[1], entry[2]
                        else:
                            continue

                        mode_norm = _normalize_mode(mode)
                        constraints_norm = _normalize_constraints(constraints)

                        # Conformer behavior: in constraint mode, constrain the scanned dihedral.
                        if (
                            mode_norm == "constraint"
                            and dih_key is not None
                            and isinstance(dih_key, (list, tuple))
                            and len(dih_key) == 4
                        ):
                            dih_tuple = tuple(int(x) for x in dih_key)
                            if dih_tuple not in constraints_norm:
                                constraints_norm.append(dih_tuple)
                        if (mode_norm == 'transition'):
                            constraints_norm = []
                        seed_structures[state_i][dih_key].append(
                            (_copy_norm_molecule(mol_obj), mode_norm, constraints_norm)
                        )

        # Reaction seeds
        if reaction_structures is not None:
            assert_msg_critical(
                isinstance(reaction_structures, (list, tuple)) and len(reaction_structures) > 0,
                "build_initial_seed_structures: reaction_structures must be a non-empty list/tuple.",
            )

            seed_structures[reaction_root].setdefault(reaction_key, [])

            for entry in reaction_structures:
                if isinstance(entry, dict):
                    mol_like = entry.get("molecule", None)
                    mode = entry.get("mode", "normal")
                    constraints = entry.get("constraints", [])
                elif isinstance(entry, (list, tuple)):
                    mol_like = entry[0]
                    mode = entry[1] if len(entry) >= 2 else "normal"
                    constraints = entry[2] if len(entry) >= 3 else []
                else:
                    mol_like = entry
                    mode = "normal"
                    constraints = []

                if mode == 'transition':
                    constraints = []
                seed_structures[reaction_root][reaction_key].append(
                    (_copy_norm_molecule(mol_like), _normalize_mode(mode), _normalize_constraints(constraints))
                )

        return seed_structures

    def determine_datapoint_density(self, imforcefieldfile):
        """
        Counts the datapoints of every root (stores the labels in
        self.qmlabels).

        :param imforcefieldfile:
            The interpolation settings dicts per root (with
            'imforcefield_file').

        :return:
            A dict {root: number of datapoints} (0 if the file does not exist).
        """

        point_densities = {state: 0 for state in self.roots_to_follow}

        for state in point_densities.keys():
            qm_datapoints = []
            if Path(imforcefieldfile[state]['imforcefield_file']).exists():
                impes_driver = InterpolationDriver(self.roots_z_matrix[state])
                impes_driver.update_settings(imforcefieldfile[state])
                self.qmlabels, z_matrix = impes_driver.read_labels()

                for label in self.qmlabels:
                    qm_data_point = InterpolationDatapoint(z_matrix)
                    qm_data_point.update_settings(imforcefieldfile[state])
                    qm_data_point.read_hdf5(imforcefieldfile[state]['imforcefield_file'], label)

                    qm_datapoints.append(qm_data_point)
                    point_densities[state] += 1

        return point_densities

    def _calculate_translation_coordinates_analysis(self, given_coordinates):
        """Center the molecule by translating its geometric center to (0, 0, 0).

           :param given_coordinates:
                Coordinate that is translated to the center.
        """
        center = np.mean(given_coordinates, axis=0)
        translated_coordinates = given_coordinates - center

        return translated_coordinates

    def calculate_distance_to_ref(self, current_coordinates, datapoint_coordinate):
        """
        The Cartesian distance between two structures after centring both on
        their geometric centre and rotating the second onto the first.

        :param current_coordinates:
            The coordinates of the first structure, shape (natoms, 3).
        :param datapoint_coordinate:
            The coordinates of the second structure (same unit).

        :return:
            The norm of the coordinate difference (unit of the input).
        """

        # First, translate the cartesian coordinates to zero
        target_coordinates = self._calculate_translation_coordinates_analysis(datapoint_coordinate)
        reference_coordinates = self._calculate_translation_coordinates_analysis(current_coordinates)

        # Then, determine the rotation matrix which
        # aligns data_point (target_coordinates)
        # to self.impes_coordinate (reference_coordinates)
        rotation_matrix_core = geometric.rotate.get_rot(
            target_coordinates,
            reference_coordinates)

        # Rotate the data point
        rotated_coordinates_core = np.dot(rotation_matrix_core, target_coordinates.T).T
        # Calculate the Cartesian distance
        ref_structure_check = reference_coordinates.copy()
        distance_core = (np.linalg.norm(rotated_coordinates_core - ref_structure_check))

        return distance_core

    def database_extracter(self, datafile, mol_labels, im_settings):
        """
        Reads all datapoints of a database file, sorted by label number.

        :param datafile:
            The database file.
        :param mol_labels:
            The atom labels of the molecule.
        :param im_settings:
            The interpolation settings dict of the root (not modified).

        :return:
            The molecules of the datapoints, the InterpolationDatapoint
            objects and the z-matrix of the file.
        """

        im_driver = InterpolationDriver()  # -> implemented Class in VeloxChem that is capable to perform interpolation calculations for a given molecule and provided z_matrix and database
        # a copy: the caller's (frozen) settings dict is not changed
        im_settings = dict(im_settings)
        im_settings['imforcefield_file'] = datafile
        im_driver.update_settings(im_settings)
        labels, z_matrix = im_driver.read_labels()

        sorted_labels = sorted(labels, key=lambda x: int(x.split('_')[1]))

        data_point_molecules = []
        datapoints = []

        for label in sorted_labels:
            impes_coordinate = InterpolationDatapoint(z_matrix)
            impes_coordinate.update_settings(im_settings)
            impes_coordinate.read_hdf5(datafile, label)  # -> read in function from the ImpesDriver object

            coordinates_in_angstrom = impes_coordinate.cartesian_coordinates * bohr_in_angstrom()

            current_molecule = Molecule(mol_labels, coordinates_in_angstrom, 'angstrom')  # -> creates a VeloxChem Molecule object

            datapoints.append(impes_coordinate)
            data_point_molecules.append(current_molecule)

        return data_point_molecules, datapoints, z_matrix

    def _confirm_database_quality(self, molecule, basis, im_settings, given_molecular_structures=None, improve=True):
        """
        Final quality check: compares the interpolated with the QM energy at
        up to nstruc_to_confirm_database_quality random trajectory structures
        (a selection at least 0.3 Angstrom RMSD from every datapoint is
        preferred). The check only reports: no datapoint is added. All
        checked structures are written to random_xyz_structures.xyz and
        the trajectory to full_xyz_traj.xyz.

        :param molecule:
            The molecule (atom labels).
        :param basis:
            The basis-set labels dict.
        :param im_settings:
            The interpolation settings dicts per root.
        :param given_molecular_structures:
            The trajectory structures, a dict {root: list of molecules}.
        :param improve:
            If True, a structure outside energy_threshold is reported.
        """

        if given_molecular_structures is None:
            raise ValueError('confirm_database_quality requires given_molecular_structures.')

        if not isinstance(given_molecular_structures, dict):
            raise ValueError('confirm_database_quality expects given_molecular_structures to be a dict keyed by root.')

        for root in self.roots_to_follow:
            database_quality = False
            drivers = self.drivers['gs'] if root == 0 else self.drivers['es']
            current_datafile = im_settings[root]['imforcefield_file']
            all_structures_root = []
            random_structure_choices_root = []
            last_qm_energies = []
            last_im_energies = []

            while not database_quality:

                all_structures_root = list(given_molecular_structures.get(root, []))
                if len(all_structures_root) == 0:
                    database_quality = True
                    continue
                else:
                    datapoint_molecules_local, _, _ = self.database_extracter(
                        current_datafile,
                        molecule.get_labels(),
                        im_settings[root],
                    )
                    rmsd = -np.inf
                    counter = 0
                    dist_ok = False
                    selected_molecules = []
                    self.ostream.print_blank()
                    self.ostream.print_header("Confirming Database Quality")
                    self.ostream.print_header('------------------')
                    self.ostream.print_blank()
                    self.ostream.print_header(f"Root {root}: Checking if current structures are well separated from the database conformations.")
                    self.ostream.print_blank()
                    self.ostream.flush()
                    while (not dist_ok) and counter <= 100:
                        selected_molecules = random.sample(
                            all_structures_root,
                            min(self.nstruc_to_confirm_database_quality, len(all_structures_root)),
                        )
                        if len(selected_molecules) == 0 and len(all_structures_root) > 0:
                            selected_molecules = random.sample(
                                all_structures_root,
                                min(self.nstruc_to_confirm_database_quality, len(all_structures_root)),
                            )
                        individual_distances = []
                        for datapoint_molecule in datapoint_molecules_local:
                            for random_struc in selected_molecules:
                                distance_norm = self.calculate_distance_to_ref(
                                    random_struc.get_coordinates_in_bohr(),
                                    datapoint_molecule.get_coordinates_in_bohr(),
                                )
                                individual_distances.append(
                                    distance_norm / np.sqrt(len(molecule.get_labels())) * bohr_in_angstrom()
                                )
                        if len(individual_distances) == 0:
                            rmsd = np.inf
                        else:
                            rmsd = min(individual_distances)
                        counter += 1
                        if rmsd >= 0.3:
                            self.ostream.print_header('------------------')
                            self.ostream.print_header(
                                f'The overall RMSD is {rmsd} -> '
                                'The current structures are well separated from the database conformations! '
                                'loop is discontinued')
                            self.ostream.print_blank()
                            self.ostream.flush()

                            dist_ok = True

                    if not dist_ok:
                        self.ostream.print_warning(
                            f'Root {root}: no random selection of {len(selected_molecules)} structures was '
                            f'separated by 0.3 Angstrom (RMSD) from all datapoints after {counter} attempts '
                            f'(closest: {rmsd:.3f} Angstrom). The last selection is used; structures close '
                            'to datapoints can underestimate the interpolation error.')
                        self.ostream.flush()

                    random_structure_choices_root = {"random_struct_info":list(selected_molecules)}

                impes_driver = self._load_interpolation_driver(molecule, root, im_settings, current_datafile)

                qm_energies = []
                im_energies = []
                database_expanded = False

                for i, mol in enumerate(random_structure_choices_root.get('random_struct_info')[:]):
                    if root == 0:
                        current_basis = MolecularBasis.read(mol, basis['gs'])
                    else:
                        current_basis = MolecularBasis.read(mol, basis['es'])
                    impes_driver.compute(mol)

                    reference_energies, _, _ = self._compute_energy(drivers[0], mol, current_basis)

                    qm_energy_val = float(reference_energies[root]) if len(reference_energies) > root else float(reference_energies[0])
                    qm_energies.append(qm_energy_val)
                    im_energies.append(impes_driver.impes_coordinate.energy)

                    self.ostream.print_blank()
                    self.ostream.print_header(f"#### Structure {i+1} ####")
                    self.ostream.print_blank()
                    self.ostream.print_block(
                        f'delta_E:  {abs(qm_energies[-1] - im_energies[-1]) * hartree_in_kcalpermol()} kcal/mol \n {abs(qm_energies[-1] - im_energies[-1]) * hartree_in_kcalpermol() / len(molecule.get_labels())} '
                        'kcal/mol per atom  \n')
                    self.ostream.flush()

                    if (
                        (abs(qm_energies[-1] - im_energies[-1]) / len(molecule.get_labels()))
                        * hartree_in_kcalpermol()
                        > self.energy_threshold
                        and improve
                    ):
                        self.ostream.print_blank()
                        self.ostream.print_block(
                            'The current structure is not within the desired threshold. It is not added to the '
                            'database. \n The structure can be found in the random_xyz_structures.xyz file')
                        self.ostream.flush()

                last_qm_energies = qm_energies
                last_im_energies = im_energies

                # the check only reports (no datapoint is added), so one pass is done
                if not database_expanded:
                    database_quality = True

            if len(all_structures_root) > 0:
                stride = max(1, int(self.nsteps / self.snapshots)) if self.snapshots else 1
                self.structures_to_xyz_file(
                    all_structures_root[::stride],
                    'full_xyz_traj.xyz',
                )
            if len(random_structure_choices_root) > 0:
                self.structures_to_xyz_file(
                    random_structure_choices_root.get("random_struct_info"),
                    'random_xyz_structures.xyz',
                    last_im_energies,
                    last_qm_energies,
                )

    def structures_to_xyz_file(self, molecules_for_xyz, structure_filename, im_energies=None, qm_energies=None):
        """Writes molecular structures to an XYZ file.

        :param molecules_for_xyz:
            A list of VeloxChem molecular objects.

        :param structure_filename:
            The name of the output file where XYZ structures will be stored.

        :param im_energies:
            An optional list of interpolation energies corresponding to each molecule.

        :param qm_energies:
            An optional list of quantum mechanical energies corresponding to each molecule.

        """

        with open(structure_filename, 'w') as file:
            pass

        for i, dyn_mol in enumerate(molecules_for_xyz):

            current_xyz_string = dyn_mol.get_xyz_string()

            xyz_lines = current_xyz_string.splitlines()

            if len(xyz_lines) >= 2 and im_energies is not None or len(xyz_lines) >= 2 and qm_energies is not None:

                if im_energies is not None and qm_energies is None:
                    xyz_lines[1] += f'Energies  IM: {im_energies[i]}'
                elif im_energies is None and qm_energies is not None:
                    xyz_lines[1] += f'Energies  IM: {qm_energies[i]}'
                else:
                    xyz_lines[1] += f'Energies  QM: {qm_energies[i]}  IM: {im_energies[i]}  delta_E: {abs(qm_energies[i] - im_energies[i])}'

            updated_xyz_string = "\n".join(xyz_lines)

            with open(structure_filename, 'a') as file:
                file.write(f"{updated_xyz_string}\n")

    def add_point(self, molecule_specific_information, interpolation_settings, symmetry_information=None):
        """
        Computes the QM energy, gradient and Hessian of every structure and
        writes it as the next datapoint point_<n+1> of the database of its
        root (and of the _org copy).

        :param molecule_specific_information:
            A list of (molecule, basis, roots, constraints[, transition])
            tuples; the 1-based constraint tuples become the important
            internal coordinates of the datapoint.
        :param interpolation_settings:
            The interpolation settings dicts per root.
        :param symmetry_information:
            Not used (kept for backward compatibility).

        :return:
            The list of written InterpolationDatapoint objects.
        """
        self.ostream.print_blank()
        self.ostream.print_header(f"Constructing initial Database {interpolation_settings[0]['imforcefield_file']} with the provided structure(s) and settings.")
        self.ostream.flush()

        written_datapoints = []

        adjusted_molecule = {'gs': [], 'es': []}

        for entries in molecule_specific_information:
            molecule = entries[0]

            if self.eq_bond_length is None:
                self.eq_bond_length = []
                for element in self.roots_z_matrix[0]['bonds']:

                    if len(element) == 2 and self.use_minimized_structures[0] and self.eq_bond_length_irc_bonds is not None and element not in self.eq_bond_length_irc_bonds:
                        self.eq_bond_length.append(molecule.get_distance([element[0] + 1, element[1] + 1], 'bohr'))
                    elif len(element) == 2 and self.use_minimized_structures[0] and self.eq_bond_length_irc_bonds is not None and element in self.eq_bond_length_irc_bonds:
                        self.eq_bond_length.append(molecule_specific_information[-1][0].get_distance([element[0] + 1, element[1] + 1], 'bohr'))
                    elif len(element) == 2 and self.use_minimized_structures[0]:
                        self.eq_bond_length.append(molecule.get_distance([element[0] + 1, element[1] + 1], 'bohr'))
                    elif len(element) == 2:
                        self.eq_bond_length.append(0.0)

            # (molecule, basis, states, 1-based constraints) per electronic-state block
            if 0 in entries[2]:
                adjusted_molecule['gs'].append((entries[0], entries[1], [0], entries[3]))
            if any(x > 0 for x in entries[2]):
                states = [state for state in entries[2] if state > 0]
                adjusted_molecule['es'].append((entries[0], entries[1], states, entries[3]))

        for key, entries in adjusted_molecule.items():
            if len(entries) == 0:
                continue
            drivers = self.drivers[key]
            assert_msg_critical(drivers is not None,
                                f'IMForceFieldGenerator.add_point: no QM drivers for the {key} state.')

            for point_molecule, point_basis, point_states, point_constraints in entries:

                energies, gradients, hessians = self._compute_qm_data(drivers, point_molecule, point_basis,
                                                                      states=point_states)

                # the mass weighting of the database (frozen in set_up_the_system), not the live attribute
                inv_sqrt_masses = None
                if interpolation_settings[point_states[0]]['use_mass_weight']:
                    masses = point_molecule.get_masses().copy()
                    masses_cart = np.repeat(masses, 3)
                    inv_sqrt_masses = 1.0 / np.sqrt(masses_cart)

                for number in range(len(energies)):
                    target_root = point_states[number]
                    target_file = interpolation_settings[target_root]['imforcefield_file']

                    z_matrix = self.roots_z_matrix[target_root]
                    interpolation_driver = InterpolationDriver(z_matrix)
                    interpolation_driver.update_settings(interpolation_settings[target_root])
                    interpolation_driver.imforcefield_file = target_file

                    sorted_labels = []
                    if Path(target_file).exists():
                        org_labels, z_matrix = interpolation_driver.read_labels()
                        labels = [label for label in org_labels if '_symmetry' not in label]
                        sorted_labels = sorted(labels, key=lambda x: int(x.split('_')[1]))
                    label = f'point_{len(sorted_labels) + 1}'
                    impes_coordinate = self._build_interpolation_datapoint(
                        point_molecule, z_matrix, interpolation_settings[target_root], point_constraints,
                        energies[number], gradients[number], hessians[number], inv_sqrt_masses)

                    self._write_interpolation_datapoint(impes_coordinate, target_root, interpolation_settings, label)
                    written_datapoints.append(impes_coordinate)
                    interpolation_driver.imforcefield_file = target_file

                    labels, z_matrix = interpolation_driver.read_labels()

                    self.ostream.print_blank()
                    self.ostream.print_header(f"Database expansion: Added point {label} to the database of root {target_root} with energy {energies[number]} Hartree.")
                    self.ostream.print_block(f"Database expansion with {', '.join(labels)}")
                    self.ostream.flush()

        return written_datapoints

    # ------------------------------------------------------------------
    # Opt-in seed stage: conformer network with error-driven relaxed scans
    # ------------------------------------------------------------------

    @staticmethod
    def _copy_molecule(molecule, reference=None):
        """
        Returns a copy of a molecule (or the molecule of an XYZ string) with
        the charge and multiplicity of a reference molecule.

        :param molecule:
            The molecule, or an XYZ string.
        :param reference:
            The molecule whose charge and multiplicity the copy gets (default:
            the molecule itself; required for an XYZ string).

        :return:
            The copy.
        """

        if reference is None:
            reference = molecule
        xyz = molecule if isinstance(molecule, str) else molecule.get_xyz_string()
        new_mol = Molecule.from_xyz_string(xyz)
        new_mol.set_charge(reference.get_charge())
        new_mol.set_multiplicity(reference.get_multiplicity())

        return new_mol

    @staticmethod
    def _circular_diff_deg(angle_a, angle_b):
        """
        The absolute circular difference of two angles.

        :param angle_a:
            The first angle in degree (a number or an array).
        :param angle_b:
            The second angle in degree.

        :return:
            The difference in degree, in [0, 180].
        """

        return abs(((angle_a - angle_b + 180.0) % 360.0) - 180.0)

    def _kt_kcal(self):
        """
        k_B T at the dynamics temperature (self.temperature).

        :return:
            k_B T in kcal/mol.
        """

        return boltzmann_in_hartreeperkelvin() * float(self.temperature) * hartree_in_kcalpermol()

    @staticmethod
    def _zmat_oriented(coordinate, entries):
        """
        The orientation of a coordinate stored in a z-matrix section.

        :param coordinate:
            The coordinate (1-based atom indices).
        :param entries:
            The z-matrix section (0-based tuples).

        :return:
            The 1-based coordinate in the orientation of the z-matrix entry
            (forward or reversed), or None if it is not an entry.
        """

        coordinate = tuple(int(x) for x in coordinate)
        zero = tuple(x - 1 for x in coordinate)
        entries = {tuple(int(x) for x in entry) for entry in entries}
        if zero in entries:
            return coordinate
        if zero[::-1] in entries:
            return coordinate[::-1]

        return None

    @staticmethod
    def _json_default(obj):
        """
        JSON conversion of the numpy types, sets and tuples in the provenance
        records (anything else as its string).

        :param obj:
            The object json cannot serialise.

        :return:
            The serialisable object.
        """

        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.floating):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, (set, tuple)):
            return list(obj)
        return str(obj)

    @staticmethod
    def _file_signature(path):
        """
        The modification time (ns) and size of a file, to detect whether a
        call has (re)written it.

        :param path:
            The file path.

        :return:
            The (mtime_ns, size) tuple, or None if the file does not exist.
        """

        try:
            stat = Path(path).stat()
        except OSError:
            return None

        return (stat.st_mtime_ns, stat.st_size)

    @staticmethod
    def _circular_states(values_deg, gap_deg):
        """
        Clusters dihedral values on the circle: the values sorted on [0, 360)
        are split wherever the gap between circular neighbours exceeds
        gap_deg (a gap equal to gap_deg does not split).

        :param values_deg:
            The dihedral values in degree (any range).
        :param gap_deg:
            The smallest gap that separates two states, in degree.

        :return:
            The states sorted by centre, each a dict with the circular mean
            'centre' (degree, [0, 360)), the 'members' (indices into
            values_deg) and the largest circular deviation of a member from
            the centre, 'spread' (degree).
        """

        vals = np.asarray(values_deg, dtype=float) % 360.0
        n = int(vals.size)
        if n == 0:
            return []

        order = np.argsort(vals, kind='stable')
        sorted_vals = vals[order]
        gaps = np.diff(np.append(sorted_vals, sorted_vals[0] + 360.0))
        splits = [i for i in range(n) if gaps[i] > float(gap_deg)]

        if len(splits) == 0:
            groups = [list(range(n))]
        else:
            groups = []
            start = (splits[-1] + 1) % n
            for split in splits:
                group = []
                i = start
                while True:
                    group.append(i)
                    if i == split:
                        break
                    i = (i + 1) % n
                groups.append(group)
                start = (split + 1) % n

        states = []
        for group in groups:
            members = [int(order[i]) for i in group]
            ang = np.deg2rad(vals[members])
            centre = float(np.rad2deg(np.arctan2(np.sin(ang).mean(), np.cos(ang).mean()))) % 360.0
            deviation = IMForceFieldGenerator._circular_diff_deg(vals[members], centre)
            states.append({'centre': centre, 'members': members, 'spread': float(deviation.max())})

        states.sort(key=lambda state: state['centre'])

        return states

    @staticmethod
    def _path_scan_spec(angle_a_deg, angle_b_deg, step_deg, full_circle=False, single=False, direction=1):
        """
        The half-scan specifications of a path from minimum A to minimum B
        along the arc that leaves A in the positive direction of the dihedral.
        Half A runs from A to the arc midpoint in the positive direction,
        half B from B to the midpoint in the negative direction. The end
        minima themselves are not scan points.

        :param angle_a_deg:
            The dihedral of minimum A in degree.
        :param angle_b_deg:
            The dihedral of minimum B in degree.
        :param step_deg:
            The requested scan step in degree (the actual step divides the
            half-arc into an integer number of steps).
        :param full_circle:
            True for a self-loop arc (one torsion state): the arc is 360 deg.
        :param single:
            True for a one-sided scan from A over the whole arc (orphan arcs).
        :param direction:
            With single, -1 scans from A in the negative direction of the
            dihedral to B (the mirrored one-sided scan); otherwise ignored.

        :return:
            A dict with the arc 'length', the 'midpoint' (degree, [0, 360)),
            the actual 'step' and, for 'A' and 'B', the signed 'start' and
            'end' angles (degree, geomeTRIC scan convention) and the number
            of scan points 'n'.
        """

        a = float(angle_a_deg) % 360.0
        b = float(angle_b_deg) % 360.0

        if single and direction < 0:
            length = (a - b) % 360.0
            n_half = max(1, int(round(length / float(step_deg))))
            step = length / n_half
            start_a = a - 360.0 if a > 180.0 else a
            return {
                'length': length,
                'midpoint': b,
                'step': step,
                'A': {'start': start_a - step, 'end': start_a - length, 'n': n_half},
                'B': {'start': start_a, 'end': start_a, 'n': 0},
            }

        length = 360.0 if full_circle else (b - a) % 360.0
        half = length if single else 0.5 * length
        n_half = max(1, int(round(half / float(step_deg))))
        step = half / n_half
        midpoint = (a + half) % 360.0
        start_a = a - 360.0 if a > 180.0 else a
        start_b = b - 360.0 if b > 180.0 else b

        spec = {
            'length': length,
            'midpoint': midpoint,
            'step': step,
            'A': {'start': start_a + step, 'end': start_a + half, 'n': n_half},
            'B': {'start': start_b - step, 'end': start_b - half, 'n': n_half},
        }
        if single:
            spec['B'] = {'start': start_b, 'end': start_b, 'n': 0}

        return spec

    @staticmethod
    def _lowest_frequency_cm1(hessian, masses, coordinates, mass_weighted, return_all=False):
        """
        Harmonic frequencies of a Cartesian Hessian whose translations and
        rotations are projected to zero (not removed). The value returned is
        the minimum over all 3N eigenvalues, so for a true minimum it is one
        of the projected zero modes (about 0 cm^-1, not the lowest
        vibration); it is meaningful as an imaginary-mode test (negative
        values below about -20 cm^-1). Negative values are imaginary modes.
        For a linear molecule one vibration is projected out as well.

        :param hessian:
            The Cartesian Hessian (Hartree/bohr^2), mass-weighted or not.
        :param masses:
            The atomic masses (amu).
        :param coordinates:
            The Cartesian coordinates (bohr), shape (natoms, 3).
        :param mass_weighted:
            True if the Hessian is already M^-1/2 H M^-1/2.
        :param return_all:
            True to return all 3N frequencies sorted (the six projected
            translations and rotations are near zero).

        :return:
            The smallest of the 3N values in cm^-1 (or the sorted array).
        """

        masses = np.asarray(masses, dtype=float)
        natoms = masses.size
        hess = np.asarray(hessian, dtype=float).reshape(3 * natoms, 3 * natoms)
        sqrt_m = np.sqrt(np.repeat(masses, 3))
        if not mass_weighted:
            hess = hess / sqrt_m[:, None] / sqrt_m[None, :]

        xyz = np.asarray(coordinates, dtype=float).reshape(natoms, 3)
        com = (masses[:, None] * xyz).sum(axis=0) / masses.sum()
        rel = xyz - com
        vectors = []
        for k in range(3):
            v = np.zeros((natoms, 3))
            v[:, k] = 1.0
            vectors.append((v * np.sqrt(masses)[:, None]).reshape(-1))
        for k in range(3):
            axis = np.zeros(3)
            axis[k] = 1.0
            v = np.cross(axis, rel)
            vectors.append((v * np.sqrt(masses)[:, None]).reshape(-1))
        q, _ = np.linalg.qr(np.array(vectors).T)
        projector = np.eye(3 * natoms) - q @ q.T
        eigenvalues = np.linalg.eigvalsh(projector @ (0.5 * (hess + hess.T)) @ projector)
        frequencies = np.sign(eigenvalues) * np.sqrt(np.abs(eigenvalues)) * 5140.487

        if return_all:
            return np.sort(frequencies)

        return float(frequencies.min())

    def _network_optimize(self, molecule, constraints, states_basis, transition=False):
        """
        Runs one optimisation of the network stage with the ground-state
        driver (SCF drivers get a converged single point first, as the
        existing seed branch does) and counts it.

        :param molecule:
            The start molecule.
        :param constraints:
            The constraint list (1-based tuples or geomeTRIC strings).
        :param states_basis:
            The basis-set labels dict.
        :param transition:
            True for a transition-state search.

        :return:
            The optimisation results dict.
        """

        driver = self.drivers['gs'][0]
        compute_args = None
        if isinstance(driver, (ScfRestrictedDriver, ScfUnrestrictedDriver)):
            basis = MolecularBasis.read(molecule, states_basis['gs'])
            _, scf_results, _ = self._compute_energy(driver, molecule, basis)
            compute_args = (basis, scf_results)

        self._network_counters['optimizations'] += 1

        return self._run_optimization(driver, molecule, constraints=list(constraints), transition=transition,
                                      index_offset=0, compute_args=compute_args)

    def _seed_conformer_network(self, molecule, states_basis):
        """
        Opt-in seed stage (seed_conformer_network): optimises every conformer
        once and adds the distinct minima, builds the network of minima
        connected by the rotation of one scanned bond, samples the arcs
        between connected minima with relaxed scans from both ends and adds
        scan points one at a time where the current interpolation is worst,
        re-evaluating the updated database after every addition. With
        reaction_structures, the normal reaction structures are minima of the
        network too, and every pair of consecutive ones that differ in their
        bonds is a reaction edge: a transition entry between them is searched
        and kept as an exact datapoint if it is a first-order saddle, the
        exchanged bonds are driven step-wise from both ends with constrained
        optimisations, and the same error-driven growth selects the points.
        An optional energy window (energy_window_kcal, or energy_window_kT
        at the dynamics temperature) drops conformer minima and torsion-scan
        points that the dynamics does not reach; reaction structures and
        reaction paths are never dropped. seed_network_settings is validated
        before the first optimisation. Writes the provenance sidecar
        seed_network_log.json in the directory of the database.

        :param molecule:
            The input molecule.
        :param states_basis:
            The basis-set labels dict.

        :return:
            The molecules_to_add_info list of compute (first entry: the first
            minimum, which starts the dynamics).
        """

        assert_msg_critical(self.nodes == 1,
                            'seed_conformer_network: only a single MPI rank is supported.')
        assert_msg_critical(bool(self.use_minimized_structures[0]),
                            'seed_conformer_network: use_minimized_structures[0] must be True '
                            '(the constrained optimisation is the core of the scheme).')
        reaction = list(self.seed_reaction_entries or []) if self.reaction_structures is not None else None
        assert_msg_critical(reaction is None or sum(1 for e in reaction if e[1] == 'normal') >= 2,
                            'seed_conformer_network: reaction_structures need at least two normal structures.')
        assert_msg_critical(list(self.roots_to_follow) == [0],
                            'seed_conformer_network: only roots_to_follow == [0] is supported.')
        assert_msg_critical(self.seed_duplicate_rmsd_thrsh is not None,
                            'seed_conformer_network: seed_duplicate_rmsd_thrsh (Angstrom) must be set.')
        assert_msg_critical((self.seed_conformers is not None and len(self.seed_conformers) > 0) or reaction is not None,
                            'seed_conformer_network: no conformers were stored by set_up_the_system.')

        # all settings are checked here, before the first optimisation
        settings = self._validate_seed_network_settings()

        root = self.roots_to_follow[0]
        # the provenance sidecar is written next to the database
        log_file = Path(self.states_interpolation_settings[root]['imforcefield_file']).parent / 'seed_network_log.json'
        dih_keys = [key for key in self.seed_structures[root].keys() if key is not None]

        # the scanned dihedrals must be z-matrix coordinates (the important coordinate of a
        # scan point is matched in the orientation stored in the z-matrix)
        for key in dih_keys:
            if self._zmat_oriented(key, self.roots_z_matrix[root]['dihedrals']) is None:
                raise ValueError(f'seed_conformer_network: scanned dihedral {key} (1-based) is not in '
                                 'the z-matrix; check drop_improper_h_dihedrals.')

        self._network_counters = {'optimizations': 0, 'single_points': 0, 'hessians': 0}
        timings = {}
        log = {
            'settings': {k: v for k, v in settings.items()},
            'thresholds': {
                'e_thrsh_kcal_per_atom': float(settings['e_thrsh_kcal_per_atom'] if settings['e_thrsh_kcal_per_atom'] is not None else self.energy_threshold),
                'g_rmsd_thrsh': float(settings['g_rmsd_thrsh'] if settings['g_rmsd_thrsh'] is not None else self.gradient_rmsd_thrsh),
                'global_constraints': [c if isinstance(c, str) else list(c) for c in (self.use_minimized_structures[1] or [])],
                'seed_duplicate_rmsd_thrsh': float(self.seed_duplicate_rmsd_thrsh),
            },
            'scan_dihedrals': [list(k) for k in dih_keys],
            'minima': [], 'dropped_conformers': [], 'bonds': [], 'orphans': [], 'paths': [], 'summary': {},
        }
        if reaction is not None:
            log['reactions'] = []

        def _write_log():
            log['summary'] = {
                'qm_calls': dict(self._network_counters),
                'wall_time_s': {k: round(v, 2) for k, v in timings.items()},
                'n_minima': len(log['minima']),
                'n_paths': len(log['paths']),
                'n_points_written': len(log['minima']) + sum(p.get('accepted', 0) for p in log['paths']),
            }
            with open(log_file, 'w') as fh:
                json.dump(log, fh, indent=1, default=self._json_default)

        self.ostream.print_blank()
        self.ostream.print_header('Seed stage: conformer network (opt-in)')
        self.ostream.print_header('--------------------------------------')
        self.ostream.print_info(f'scanned dihedrals (1-based): {dih_keys}')
        self.ostream.print_info(f'settings: {settings}')
        self.ostream.flush()

        # stage 1: minima
        t0 = time.time()
        minima, molecules_to_add_info = self._seed_minima(molecule, states_basis, dih_keys, settings, log,
                                                          reaction=reaction)
        timings['minima'] = time.time() - t0
        _write_log()

        if len(dih_keys) == 0 and reaction is None:
            self.ostream.print_info('seed_conformer_network: no rotatable bond, only the minimum was added.')
            self.ostream.flush()
            return molecules_to_add_info

        # stage 2: the network
        network = self._conformer_network_paths(minima, dih_keys, settings, log)
        if reaction is not None:
            # reaction edges after the torsion paths; transition states are written here, before the
            # in-memory driver is loaded
            network['paths'] += self._reaction_paths(reaction, minima, states_basis, dih_keys, log,
                                                     start_index=len(network['paths']))
        _write_log()

        # stage 3: error-driven growth along every path
        t0 = time.time()
        datafile = self.states_interpolation_settings[root]['imforcefield_file']
        driver = self._load_interpolation_driver(molecule, root, self.states_interpolation_settings, datafile)
        earlier = []
        for path in network['paths']:
            if path.get('kind') == 'reaction':
                path_log = {k: path[k] for k in ('index', 'kind', 'A', 'B', 'waypoint')}
                path_log['formed'] = [[a + 1, b + 1] for a, b in path['formed']]
                path_log['broken'] = [[a + 1, b + 1] for a, b in path['broken']]
                self.ostream.print_blank()
                self.ostream.print_info(
                    f"path {path['index']}: reaction from minimum {path['A']} to minimum {path['B']}; bonds formed "
                    f"{path_log['formed']}, broken {path_log['broken']} (1-based)" +
                    (f"; transition state {path['waypoint']}" if path['waypoint'] is not None else ''))
                self.ostream.flush()
                candidates = self._reaction_path_scan(path, minima, states_basis, settings, path_log)
                self._grow_path(path, candidates, driver, minima, earlier, states_basis, settings, path_log)
                earlier.append(({path['A'], path['B']}, candidates))
                log['paths'].append(path_log)
                _write_log()
                continue
            path_log = {k: path[k] for k in ('index', 'bond', 'arc', 'A', 'B', 'orphan')}
            path_log['dihedral'] = list(dih_keys[path['bond']])
            self.ostream.print_blank()
            self.ostream.print_info(
                f"path {path['index']}: bond {dih_keys[path['bond']]} arc {path['arc']} "
                f"from minimum {path['A']} ({minima[path['A']]['angles'][path['bond']]:.1f} deg) "
                f"to minimum {path['B']} ({minima[path['B']]['angles'][path['bond']]:.1f} deg)")
            self.ostream.flush()
            candidates = self._relaxed_path_scan(path, network, minima, states_basis, settings, path_log)
            if path['orphan']:
                end_minimum = self._network_orphan_end(path, candidates, minima, molecule, states_basis, dih_keys, log, path_log)
                if end_minimum is None:
                    log['paths'].append(path_log)
                    _write_log()
                    continue
                path['B'] = end_minimum
                path_log['B'] = end_minimum
                # the new minimum is a datapoint now: reload it into the in-memory driver
                driver = self._load_interpolation_driver(molecule, root, self.states_interpolation_settings, datafile)
            candidates = self._path_energy_window(path, candidates, minima, settings, path_log)
            self._grow_path(path, candidates, driver, minima, earlier, states_basis, settings, path_log)
            earlier.append(({path['A'], path['B']}, candidates))
            log['paths'].append(path_log)
            _write_log()
        timings['growth'] = time.time() - t0
        _write_log()

        n_written = log['summary']['n_points_written']
        self.ostream.print_blank()
        self.ostream.print_header('Seed stage: conformer network summary')
        self.ostream.print_info(f"minima: {len(minima)}; paths: {len(network['paths'])}; orphan arcs: {len(log['orphans'])}")
        self.ostream.print_info(f"points written: {n_written}; QM calls: {self._network_counters}")
        self.ostream.print_info(f"wall time (s): {log['summary']['wall_time_s']}")
        self.ostream.flush()

        return molecules_to_add_info

    def _seed_minima(self, molecule, states_basis, dih_keys, settings, log, reaction=None):
        """
        Stage 1: optimises every stored conformer once (global constraints
        use_minimized_structures[1] only), keeps the converged and distinct
        minima (_is_same_minimum with seed_duplicate_rmsd_thrsh, optional
        energy window) and adds them through add_point without important
        coordinates. Minima with an imaginary frequency below -20 cm^-1 stay
        datapoints but are no network nodes. With an energy window
        (_energy_window_kcal), conformer minima more than the window above the
        lowest minimum are dropped before their Hessian is computed.

        :param molecule:
            The input molecule.
        :param states_basis:
            The basis-set labels dict.
        :param dih_keys:
            The scanned dihedrals (1-based).
        :param settings:
            The checked seed_network_settings.
        :param log:
            The provenance dict.
        :param reaction:
            The reaction seed entries (mol, mode, constraints) or None. The
            normal ones are optimised after the conformers; one that is the
            same minimum as a kept one is recorded in that minimum's
            'reaction' list. The energy window never drops a reaction minimum.

        :return:
            The list of minima dicts and the molecules_to_add_info list.
        """

        global_constraints = list(self.use_minimized_structures[1] or [])
        starts = [(conformer, 'conformer', ci) for ci, conformer in enumerate(self.seed_conformers or [])]
        if reaction is not None:
            # after the conformers: the last structure (usually the product) sets the eq lengths of the
            # reaction bonds, as in the existing reaction seeding
            starts += [(entry[0], 'reaction', ri) for ri, entry in enumerate(reaction) if entry[1] == 'normal']
        kept = []
        for start_mol, origin, ci in starts:
            name = 'conformer' if origin == 'conformer' else 'reaction structure'
            mol = self._copy_molecule(start_mol)
            try:
                opt_results = self._network_optimize(mol, global_constraints, states_basis)
            except Exception as exc:
                # any failure of geomeTRIC or of the QM driver drops the structure; the type is logged so
                # that a programming error is not mistaken for a failed optimisation
                reason = f'optimisation failed: {type(exc).__name__}: {exc}'
                log['dropped_conformers'].append({origin: ci, 'reason': reason})
                self.ostream.print_info(f'{name} {ci}: {reason}; dropped.')
                self.ostream.flush()
                continue
            opt_mol = opt_results['final_molecule']
            energy = None
            if 'opt_energies' in opt_results and len(opt_results['opt_energies']) > 0:
                energy = float(opt_results['opt_energies'][-1])
            if self._is_same_minimum(opt_mol, [m['molecule'] for m in kept]):
                if origin == 'reaction':
                    match = next(m for m in kept if self._is_same_minimum(opt_mol, [m['molecule']]))
                    match.setdefault('reaction', []).append(ci)
                    log['dropped_conformers'].append({origin: ci, 'reason': 'same minimum as a kept minimum'})
                    self.ostream.print_info(f'{name} {ci}: same minimum as a kept minimum; recorded there.')
                else:
                    log['dropped_conformers'].append({'conformer': ci, 'reason': 'same minimum as a kept conformer'})
                    self.ostream.print_info(f'conformer {ci}: same minimum as a kept conformer; dropped.')
                self.ostream.flush()
                continue
            kept.append({'conformer': ci if origin == 'conformer' else None, 'molecule': opt_mol, 'energy': energy})
            if origin == 'reaction':
                kept[-1]['reaction'] = [ci]
            self.ostream.print_info(f'{name} {ci}: distinct minimum {len(kept) - 1}' +
                                    (f' (E = {energy:.8f} Hartree)' if energy is not None else ''))
            self.ostream.flush()

        assert_msg_critical(len(kept) > 0, 'seed_conformer_network: no conformer optimisation converged.')

        window = self._energy_window_kcal(settings)
        if all(m['energy'] is not None for m in kept):
            kt = self._kt_kcal()
            e_min = min(m['energy'] for m in kept)
            e_max = (max(m['energy'] for m in kept) - e_min) * hartree_in_kcalpermol()
            self.ostream.print_info(f'{len(kept)} minima within {e_max:.2f} kcal/mol ({e_max / kt:.1f} k_B T at '
                                    f'{float(self.temperature):g} K) of the lowest one.')
            if window is not None:
                log['energy_window'] = {'kcal_per_mol': window, 'temperature_K': float(self.temperature),
                                        'kT_kcal_per_mol': kt, 'kT_multiple': window / kt}
                inside = []
                for m in kept:
                    de = (m['energy'] - e_min) * hartree_in_kcalpermol()
                    if not m.get('reaction') and de > window:
                        log['dropped_conformers'].append({'conformer': m['conformer'], 'reason': 'outside the energy window',
                                                          'dE_kcal_per_mol': de})
                        self.ostream.print_info(f"conformer {m['conformer']}: {de:.2f} kcal/mol ({de / kt:.1f} k_B T) above "
                                                f"the lowest minimum, outside the energy window; dropped (not reached at "
                                                f"{float(self.temperature):g} K).")
                    else:
                        inside.append(m)
                n_out = len(kept) - len(inside)
                kept = inside
                window_key = 'energy_window_kT' if settings.get('energy_window_kT') is not None else 'energy_window_kcal'
                self.ostream.print_info(f'energy window {window:.2f} kcal/mol = {window / kt:.1f} k_B T at '
                                        f'{float(self.temperature):g} K: {n_out} conformer minimum(s) dropped. Raise '
                                        f"seed_network_settings['{window_key}'] (or set it to None) to include "
                                        'higher-energy conformers.')
                if (self.metadynamics_settings or {}).get('enabled', False):
                    self.ostream.print_info('metadynamics is enabled: the biased dynamics can reach energies above '
                                            'the thermal window; consider a larger energy_window_kT.')
            self.ostream.flush()

        molecules_to_add_info = []
        for m in kept:
            basis = MolecularBasis.read(m['molecule'], states_basis['gs'])
            molecules_to_add_info.append((m['molecule'], basis, self.roots_to_follow, []))

        datapoints = self.add_point(molecules_to_add_info, self.states_interpolation_settings)
        self._network_counters['single_points'] += len(kept)
        self._network_counters['hessians'] += len(kept)

        masses = np.asarray(molecule.get_masses(), dtype=float)
        for index, (m, dp) in enumerate(zip(kept, datapoints)):
            m['index'] = index
            m['energy'] = float(dp.energy)
            m['label'] = dp.point_label if dp.point_label is not None else f'point_{index + 1}'
            m['angles'] = [float(m['molecule'].get_dihedral_in_degrees(key)) % 360.0 for key in dih_keys]
            m['lowest_freq_cm1'] = self._lowest_frequency_cm1(dp.hessian, masses, dp.cartesian_coordinates,
                                                              dp.inv_sqrt_masses is not None)
            m['node'] = bool(m['lowest_freq_cm1'] > -20.0)
            if not m['node']:
                self.ostream.print_info(f"minimum {index}: imaginary frequency {m['lowest_freq_cm1']:.1f} cm^-1; "
                                        'kept as a datapoint, not as a network node.')
            log['minima'].append({k: m[k] for k in ('index', 'conformer', 'energy', 'label', 'angles',
                                                    'lowest_freq_cm1', 'node')})
            if m.get('reaction'):
                log['minima'][-1]['reaction'] = list(m['reaction'])
        self.ostream.print_info(f'{len(kept)} distinct minima added as point_1 .. point_{len(kept)}.')
        self.ostream.flush()

        return kept, molecules_to_add_info

    def _validate_seed_network_settings(self):
        """
        Checks seed_network_settings before any QM of the network seed stage:
        unknown and missing keys, types and ranges, and that at most one
        energy window is set.

        :return:
            A copy of the checked settings.
        """

        user = self.seed_network_settings if self.seed_network_settings is not None else {}
        assert_msg_critical(isinstance(user, dict), 'seed_conformer_network: seed_network_settings must be a dict.')
        unknown = sorted(set(user) - set(self._seed_network_settings_keys))
        assert_msg_critical(len(unknown) == 0,
                            f'seed_conformer_network: unknown seed_network_settings key(s) {unknown}; '
                            f'the keys are {sorted(self._seed_network_settings_keys)}.')
        # the defaults live in seed_network_settings itself: a replaced dict must keep all keys
        missing = [key for key in self._seed_network_settings_keys if key not in user]
        assert_msg_critical(len(missing) == 0,
                            f'seed_conformer_network: seed_network_settings misses the key(s) {missing}; '
                            'change single entries of seed_network_settings instead of replacing the dict.')

        settings = dict(user)

        def _number(value):
            return (isinstance(value, numbers.Real) and not isinstance(value, (bool, np.bool_))
                    and bool(np.isfinite(value)))

        for key in ('scan_step_deg', 'state_gap_deg', 'reaction_scan_step_angstrom'):
            assert_msg_critical(_number(settings[key]) and settings[key] > 0,
                                f"seed_conformer_network: seed_network_settings['{key}'] must be a positive "
                                f'number (got {settings[key]!r}).')
        for key in ('e_thrsh_kcal_per_atom', 'g_rmsd_thrsh', 'energy_window_kcal', 'energy_window_kT'):
            assert_msg_critical(settings[key] is None or (_number(settings[key]) and settings[key] > 0),
                                f"seed_conformer_network: seed_network_settings['{key}'] must be None or a "
                                f'positive number (got {settings[key]!r}).')
        for key in ('max_points_per_path', 'max_rejections_per_path'):
            assert_msg_critical(_number(settings[key]) and float(settings[key]).is_integer() and settings[key] >= 1,
                                f"seed_conformer_network: seed_network_settings['{key}'] must be an integer "
                                f'>= 1 (got {settings[key]!r}).')
        assert_msg_critical(isinstance(settings['include_orphan_arcs'], (bool, np.bool_)),
                            "seed_conformer_network: seed_network_settings['include_orphan_arcs'] must be True "
                            f"or False (got {settings['include_orphan_arcs']!r}).")

        # setting both windows aborts here, not after the conformer optimisations
        self._energy_window_kcal(settings)

        return settings

    def _energy_window_kcal(self, settings):
        """
        The energy window of the seed stage: energy_window_kcal, or
        energy_window_kT times k_B T at the dynamics temperature
        (self.temperature). A region more than the window above the lowest
        minimum has a Boltzmann factor below exp(-window / k_B T) and is not
        reached by the dynamics at that temperature.

        :param settings:
            The checked seed_network_settings.

        :return:
            The window in kcal/mol, or None (no window).
        """

        window_kcal = settings.get('energy_window_kcal')
        window_kt = settings.get('energy_window_kT')
        assert_msg_critical(window_kcal is None or window_kt is None,
                            'seed_conformer_network: set energy_window_kcal or energy_window_kT, not both.')
        if window_kt is not None:
            return float(window_kt) * self._kt_kcal()
        return None if window_kcal is None else float(window_kcal)

    def _path_energy_window(self, path, candidates, minima, settings, path_log):
        """
        Barrier check of a torsion path before its growth: logs the highest
        relaxed-scan point above the lower end minimum (the barrier estimate)
        and, with an energy window, removes the scan points more than the
        window above the lowest minimum, which the dynamics does not reach at
        its temperature. They get no QM single point, no Hessian and are not
        part of the acceptance test. Reaction paths are not filtered.

        :param path:
            The path dict.
        :param candidates:
            The candidate list of _relaxed_path_scan.
        :param minima:
            The minima dicts.
        :param settings:
            The checked seed_network_settings.
        :param path_log:
            The provenance dict of this path.

        :return:
            The candidates inside the window (all candidates without a window).
        """

        kcal = hartree_in_kcalpermol()
        kt = self._kt_kcal()
        temperature = float(self.temperature)
        energies = [c['energy'] for c in candidates if c.get('energy') is not None]
        ends = [minima[path['A']].get('energy'), minima[path['B']].get('energy')]
        if len(energies) == 0 or any(e is None for e in ends):
            return candidates

        barrier = (max(energies) - min(ends)) * kcal
        path_log['scan_max_above_lower_end_kcal'] = float(barrier)
        self.ostream.print_info(f'highest scan point {barrier:.2f} kcal/mol ({barrier / kt:.1f} k_B T at '
                                f'{temperature:g} K) above the lower end minimum.')

        window = self._energy_window_kcal(settings)
        if window is None:
            self.ostream.flush()
            return candidates

        e_min = min(m['energy'] for m in minima if m.get('energy') is not None)
        inside = [c for c in candidates if c.get('energy') is None or (c['energy'] - e_min) * kcal <= window]
        inside_ids = {id(c) for c in inside}
        path_log['above_energy_window'] = [c['index'] for c in candidates if id(c) not in inside_ids]
        if len(inside) < len(candidates):
            window_key = 'energy_window_kT' if settings.get('energy_window_kT') is not None else 'energy_window_kcal'
            self.ostream.print_info(f'{len(candidates) - len(inside)} of {len(candidates)} scan points more than '
                                    f'{window:.2f} kcal/mol ({window / kt:.1f} k_B T) above the lowest minimum are not '
                                    f'used (not reached at {temperature:g} K). Raise '
                                    f"seed_network_settings['{window_key}'] to include them.")
        self.ostream.flush()

        return inside

    def _conformer_network_paths(self, minima, dih_keys, settings, log):
        """
        Stage 2: clusters the minima's dihedrals into torsion states per
        scanned bond (circular gap state_gap_deg), builds the arcs between
        circularly adjacent states (a one-state bond gets one full-circle
        arc) and the paths (bond, arc, A, B) between minima that are in the
        arc's two states and share the state of every other scanned bond.
        Minima without a partner are logged as orphans; their arcs are
        scanned only with include_orphan_arcs.

        :param minima:
            The minima dicts of _seed_minima.
        :param dih_keys:
            The scanned dihedrals (1-based).
        :param settings:
            The checked seed_network_settings.
        :param log:
            The provenance dict.

        :return:
            A dict with the 'bonds' (states and arcs), the 'state_of' table
            (minimum index -> list of state indices per bond) and the sorted
            'paths'.
        """

        gap = float(settings['state_gap_deg'])
        nodes = [m['index'] for m in minima if m['node']]
        # torsion paths stay within one connectivity (reaction minima may differ in their bonds), and a
        # scanned dihedral is clustered only over the minima in which it is a chain of bonds
        conn = {m['index']: np.asarray(m['molecule'].get_connectivity_matrix()) for m in minima
                if m['index'] in nodes and hasattr(m.get('molecule'), 'get_connectivity_matrix')}

        def _chain(i, key):
            if i not in conn:
                return True
            return all(conn[i][int(key[c]) - 1, int(key[c + 1]) - 1] > 0 for c in range(3))

        bonds = []
        state_of = {i: [] for i in nodes}
        for b, key in enumerate(dih_keys):
            bond_nodes = [i for i in nodes if _chain(i, key)]
            states = self._circular_states([minima[i]['angles'][b] for i in bond_nodes], gap)
            for state in states:
                state['members'] = [bond_nodes[j] for j in state['members']]
            k = len(states)
            arcs = [(i, (i + 1) % k) for i in range(k)] if k > 1 else ([(0, 0)] if k == 1 else [])
            bonds.append({'key': key, 'states': states, 'arcs': arcs})
            state_index = {}
            for si, state in enumerate(states):
                for m in state['members']:
                    state_index[m] = si
            for i in nodes:
                state_of[i].append(state_index.get(i))
            self.ostream.print_info(
                f"bond {key}: {k} state(s) at {[round(s['centre'], 1) for s in states]} deg, "
                f"members per state {[len(s['members']) for s in states]}, spread {[round(s['spread'], 1) for s in states]} deg")
            log['bonds'].append({'key': list(key), 'states': [dict(s) for s in states], 'arcs': arcs})
        self.ostream.flush()

        def _same_other_states(a, b, bond):
            if a in conn and b in conn and not np.array_equal(conn[a], conn[b]):
                return False
            return all(state_of[a][c] == state_of[b][c] for c in range(len(bonds)) if c != bond)

        paths = []
        for b, bond in enumerate(bonds):
            for ai, (s1, s2) in enumerate(bond['arcs']):
                members1 = bond['states'][s1]['members']
                members2 = bond['states'][s2]['members']
                if s1 == s2:
                    for a in members1:
                        paths.append({'bond': b, 'arc': ai, 'A': a, 'B': a, 'orphan': False, 'direction': +1})
                    continue
                for a in members1:
                    partners = [p for p in members2 if _same_other_states(a, p, b)]
                    if len(partners) == 0:
                        log['orphans'].append({'bond': b, 'arc': ai, 'minimum': a, 'direction': +1, 'scanned': False})
                        if settings.get('include_orphan_arcs', False):
                            paths.append({'bond': b, 'arc': ai, 'A': a, 'B': a, 'orphan': True, 'direction': +1})
                    for p in partners:
                        paths.append({'bond': b, 'arc': ai, 'A': a, 'B': p, 'orphan': False, 'direction': +1})
                for p in members2:
                    if not any(_same_other_states(a, p, b) for a in members1):
                        log['orphans'].append({'bond': b, 'arc': ai, 'minimum': p, 'direction': -1, 'scanned': False})
                        if settings.get('include_orphan_arcs', False):
                            paths.append({'bond': b, 'arc': ai, 'A': p, 'B': p, 'orphan': True, 'direction': -1})

        paths.sort(key=lambda p: (p['bond'], p['arc'], min(minima[p['A']]['energy'], minima[p['B']]['energy']),
                                  p['A'], p['B'], p['direction']))
        for index, path in enumerate(paths):
            path['index'] = index

        n_real = sum(1 for p in paths if not p['orphan'])
        self.ostream.print_info(f"network: {n_real} path(s) between existing minima, {len(log['orphans'])} orphan "
                                f"(minimum, arc direction) pair(s)" +
                                (' (scanned)' if settings.get('include_orphan_arcs', False) else ' (not scanned)'))
        self.ostream.flush()

        return {'bonds': bonds, 'state_of': state_of, 'paths': paths}

    def _relaxed_path_scan(self, path, network, minima, states_basis, settings, path_log):
        """
        Stage 3a: relaxed scans of one path from both end minima to the arc
        midpoint (geomeTRIC 'scan dihedral'; besides the scanned dihedral only
        the global constraints use_minimized_structures[1] are frozen at each
        step, as in the optimisation of the minima; step-wise frozen-dihedral
        optimisations as the fallback). Every scan point is checked for unchanged connectivity and
        for the other scanned torsions staying in the state of the start
        minimum (spread plus half the state gap); a half-scan stops at the
        first violation. The RMSD between the last points of the two halves
        is logged as the hysteresis.

        :param path:
            The path dict.
        :param network:
            The network dict of _conformer_network_paths.
        :param minima:
            The minima dicts.
        :param states_basis:
            The basis-set labels dict.
        :param settings:
            The checked seed_network_settings.
        :param path_log:
            The provenance dict of this path.

        :return:
            The candidate list: dicts with 'half', 'target', 'angle',
            'energy' (Hartree, from the scan), 'molecule' and 'index' (the
            position in the list).
        """

        b = path['bond']
        bond = network['bonds'][b]
        key = bond['key']
        s1, s2 = bond['arcs'][path['arc']]
        step = float(settings['scan_step_deg'])
        tolerance = 0.5 * float(settings['state_gap_deg'])

        # the global constraints hold at every scan point, as in the optimisation of the minima; a global
        # constraint on the scanned dihedral itself is replaced by the scan
        key_set = {tuple(int(x) for x in key), tuple(int(x) for x in key)[::-1]}
        global_constraints = [c for c in (self.use_minimized_structures[1] or [])
                              if isinstance(c, str) or tuple(int(x) for x in c) not in key_set]

        a_min = minima[path['A']]
        b_min = minima[path['B']]
        if path['orphan']:
            # scan from the orphan to the centre of the adjacent state in the arc direction
            other = bond['states'][s2 if path['direction'] > 0 else s1]['centre']
            spec = self._path_scan_spec(a_min['angles'][b], other, step, single=True, direction=path['direction'])
        else:
            spec = self._path_scan_spec(a_min['angles'][b], b_min['angles'][b], step, full_circle=(s1 == s2))
        path_log['spec'] = spec
        path_log['scan_points'] = []
        path_log['violations'] = []

        candidates = []
        last_points = {}
        for half in ('A', 'B'):
            h = spec[half]
            if h['n'] == 0:
                continue
            start_min = a_min if half == 'A' else b_min
            conn_ref = start_min['molecule'].get_connectivity_matrix()
            targets = np.linspace(h['start'], h['end'], h['n'])
            geometries, energies, mode = [], [], 'scan'
            # geomeTRIC writes scan-final.xyz into the working directory for every scan with more than one
            # point; the file of this scan is removed below (a file the scan did not write is kept)
            scan_file = Path('scan-final.xyz')
            scan_file_before = self._file_signature(scan_file)
            try:
                constraint = (f"scan dihedral {key[0]} {key[1]} {key[2]} {key[3]} "
                              f"{h['start']:.4f} {h['end']:.4f} {h['n']}")
                opt_results = self._network_optimize(start_min['molecule'], global_constraints + [constraint],
                                                     states_basis)
                self._network_counters['optimizations'] += h['n'] - 1
                for xyz, energy in zip(opt_results['scan_geometries'], opt_results['scan_energies']):
                    geometries.append(self._copy_molecule(xyz, start_min['molecule']))
                    energies.append(float(energy))
            except Exception as exc:
                mode = 'stepwise'
                self.ostream.print_info(f"half {half}: relaxed scan failed ({type(exc).__name__}: {exc}); "
                                        "step-wise fallback.")
                self.ostream.flush()
                geometries, energies = [], []
                previous = start_min['molecule']
                for target in targets:
                    trial = self._copy_molecule(previous)
                    trial.set_dihedral_in_degrees(key, float(target), verbose=False)
                    err = self._circular_diff_deg(trial.get_dihedral_in_degrees(key), target)
                    if err > 2.0 and self.seed_scan_clash_relief:
                        trial, _ = self._place_scan_seed([previous, start_min['molecule']], key, float(target) % 360.0)
                        if trial is not None:
                            err = self._circular_diff_deg(trial.get_dihedral_in_degrees(key), target)
                    if trial is None or err > 2.0:
                        path_log['violations'].append({'half': half, 'target': float(target), 'reason': 'placement blocked'})
                        break
                    try:
                        opt_results = self._network_optimize(trial, global_constraints + [key], states_basis)
                    except Exception as exc2:
                        path_log['violations'].append({'half': half, 'target': float(target),
                                                       'reason': f'optimisation failed: {type(exc2).__name__}: {exc2}'})
                        break
                    geometries.append(opt_results['final_molecule'])
                    energies.append(float(opt_results['opt_energies'][-1]))
                    previous = opt_results['final_molecule']
            if scan_file.exists() and self._file_signature(scan_file) != scan_file_before:
                scan_file.unlink()

            # one geometry per target at most (geomeTRIC returns one per scan point, the fallback stops early)
            for target, mol, energy in zip(targets, geometries, energies):
                target = float(target)
                actual = float(mol.get_dihedral_in_degrees(key))
                dihedral_ok = self._circular_diff_deg(actual, target) <= 5.0
                connectivity_ok = bool(np.array_equal(mol.get_connectivity_matrix(), conn_ref))
                states_ok = True
                drift = {}
                for c, other_bond in enumerate(network['bonds']):
                    # no state where the dihedral is not a bonded chain (another connectivity class)
                    if c == b or network['state_of'][start_min['index']][c] is None:
                        continue
                    state = other_bond['states'][network['state_of'][start_min['index']][c]]
                    value = float(mol.get_dihedral_in_degrees(other_bond['key']))
                    deviation = self._circular_diff_deg(value, state['centre'])
                    drift[str(list(other_bond['key']))] = round(deviation, 2)
                    if deviation > state['spread'] + tolerance:
                        states_ok = False
                record = {'half': half, 'mode': mode, 'target': target % 360.0, 'angle': actual % 360.0,
                          'energy': energy, 'dihedral_ok': dihedral_ok, 'connectivity_ok': connectivity_ok,
                          'states_ok': states_ok, 'other_torsion_drift_deg': drift}
                path_log['scan_points'].append(record)
                if not (dihedral_ok and connectivity_ok and states_ok):
                    path_log['violations'].append({'half': half, 'target': target % 360.0,
                                                   'reason': 'dihedral' if not dihedral_ok else
                                                   ('connectivity' if not connectivity_ok else 'torsion state')})
                    self.ostream.print_info(f"half {half}: continuity violated at {target % 360.0:.1f} deg "
                                            f"({path_log['violations'][-1]['reason']}); half-scan stopped.")
                    self.ostream.flush()
                    break
                candidates.append({'half': half, 'target': target % 360.0, 'angle': actual % 360.0,
                                   'energy': energy, 'molecule': mol, 'index': len(candidates)})
                last_points[half] = mol

        # the hysteresis is defined only when both halves reached the midpoint
        path_log['hysteresis_rmsd_angstrom'] = None
        complete = all(sum(1 for c in candidates if c['half'] == half) == spec[half]['n'] for half in ('A', 'B'))
        if 'A' in last_points and 'B' in last_points and complete:
            rmsd, _, _ = svd_superimpose(last_points['A'].get_coordinates_in_angstrom(),
                                         last_points['B'].get_coordinates_in_angstrom())
            path_log['hysteresis_rmsd_angstrom'] = float(rmsd)
            self.ostream.print_info(f'midpoint hysteresis (RMSD between the last points of the halves): {rmsd:.4f} Angstrom')
        self.ostream.print_info(f"{len(candidates)} scan point(s) on this path "
                                f"({sum(1 for c in candidates if c['half'] == 'A')} from A, "
                                f"{sum(1 for c in candidates if c['half'] == 'B')} from B)")
        self.ostream.flush()

        return candidates

    def _network_orphan_end(self, path, candidates, minima, molecule, states_basis, dih_keys, log, path_log):
        """
        Optimises the last scan point of an orphan arc (only the global
        constraints use_minimized_structures[1]) and screens it like a
        conformer. A new distinct minimum becomes a minimum datapoint and the
        end of the path; its scan point is removed from the candidates.

        :param path:
            The orphan path dict.
        :param candidates:
            The candidate list of _relaxed_path_scan (changed in place).
        :param minima:
            The minima dicts (the new minimum is appended).
        :param molecule:
            The input molecule (for the masses).
        :param states_basis:
            The basis-set labels dict.
        :param dih_keys:
            The scanned dihedrals (1-based).
        :param log:
            The provenance dict (its orphan entry is marked as scanned).
        :param path_log:
            The provenance dict of this path.

        :return:
            The index of the new minimum, or None (the scan is discarded).
        """

        for orphan in log['orphans']:
            if orphan['bond'] == path['bond'] and orphan['arc'] == path['arc'] and orphan['minimum'] == path['A'] \
                    and orphan['direction'] == path['direction']:
                orphan['scanned'] = True
        if len(candidates) == 0:
            path_log['status'] = 'orphan scan produced no point'
            return None
        try:
            opt_results = self._network_optimize(candidates[-1]['molecule'], list(self.use_minimized_structures[1] or []), states_basis)
        except Exception as exc:
            path_log['status'] = f'orphan end optimisation failed: {type(exc).__name__}: {exc}'
            return None
        end_mol = opt_results['final_molecule']
        if self._is_same_minimum(end_mol, [m['molecule'] for m in minima]):
            path_log['status'] = 'orphan end is a known minimum; scan discarded'
            self.ostream.print_info('orphan arc: the end point relaxed into a known minimum; scan discarded.')
            self.ostream.flush()
            return None

        basis = MolecularBasis.read(end_mol, states_basis['gs'])
        datapoints = self.add_point([(end_mol, basis, self.roots_to_follow, [])], self.states_interpolation_settings)
        self._network_counters['single_points'] += 1
        self._network_counters['hessians'] += 1
        dp = datapoints[0]
        masses = np.asarray(molecule.get_masses(), dtype=float)
        new = {'conformer': None, 'molecule': end_mol, 'energy': float(dp.energy), 'index': len(minima),
               'label': dp.point_label if dp.point_label is not None else f'point_{len(minima) + 1}',
               'angles': [float(end_mol.get_dihedral_in_degrees(key)) % 360.0 for key in dih_keys],
               'lowest_freq_cm1': self._lowest_frequency_cm1(dp.hessian, masses, dp.cartesian_coordinates,
                                                             dp.inv_sqrt_masses is not None)}
        new['node'] = bool(new['lowest_freq_cm1'] > -20.0)
        minima.append(new)
        log['minima'].append({k: new[k] for k in ('index', 'conformer', 'energy', 'label', 'angles', 'lowest_freq_cm1', 'node')})
        candidates.pop()
        self.ostream.print_info(f"orphan arc: new minimum {new['index']} added as {new['label']}.")
        self.ostream.flush()

        return new['index']

    @staticmethod
    def _bond_changes(conn_a, conn_b):
        """
        The bonds formed and broken between two connectivity matrices.

        :param conn_a:
            The connectivity matrix of the start structure.
        :param conn_b:
            The connectivity matrix of the end structure.

        :return:
            The formed and the broken bonds, each a list of 0-based (i, j)
            tuples with i < j.
        """

        conn_a = np.asarray(conn_a)
        conn_b = np.asarray(conn_b)
        natoms = conn_a.shape[0]
        formed = [(i, j) for i in range(natoms) for j in range(i + 1, natoms) if conn_b[i, j] and not conn_a[i, j]]
        broken = [(i, j) for i in range(natoms) for j in range(i + 1, natoms) if conn_a[i, j] and not conn_b[i, j]]

        return formed, broken

    def _reaction_paths(self, reaction, minima, states_basis, dih_keys, log, start_index):
        """
        Builds the reaction edges of the network: every pair of consecutive
        normal reaction structures whose minima differ in their bonds (at
        least one formed and one broken). A transition entry between them is
        searched and, if it is a first-order saddle between the two minima,
        written as an exact datapoint and used as the meeting point of the
        two half-scans. Pairs without a bond exchange are conformational and
        left to the torsion network.

        :param reaction:
            The reaction seed entries (mol, mode, constraints).
        :param minima:
            The minima dicts (each with its 'reaction' entry indices).
        :param states_basis:
            The basis-set labels dict.
        :param dih_keys:
            The scanned dihedrals (1-based).
        :param log:
            The provenance dict.
        :param start_index:
            The index of the first reaction path.

        :return:
            The list of reaction path dicts.
        """

        nodes = []
        for ri, entry in enumerate(reaction):
            if entry[1] == 'normal':
                index = next((m['index'] for m in minima if ri in (m.get('reaction') or [])), None)
                nodes.append(('normal', index, ri))
            elif entry[1] == 'transition':
                nodes.append(('transition', entry[0], ri))
            else:
                log['reactions'].append({'entry': ri, 'skipped': f"mode '{entry[1]}' is not used by the network stage"})

        normals = [(pos, node) for pos, node in enumerate(nodes) if node[0] == 'normal']
        paths = []
        for (pa, na), (pb, nb) in zip(normals[:-1], normals[1:]):
            entry_log = {'from_entry': na[2], 'to_entry': nb[2], 'A': na[1], 'B': nb[1]}
            log['reactions'].append(entry_log)
            if na[1] is None or nb[1] is None or na[1] == nb[1]:
                entry_log['skipped'] = 'no distinct minima at both ends'
                self.ostream.print_info(f"reaction entries {na[2]} -> {nb[2]}: no distinct minima at both ends; skipped.")
                continue
            min_a, min_b = minima[na[1]], minima[nb[1]]
            formed, broken = self._bond_changes(min_a['molecule'].get_connectivity_matrix(),
                                                min_b['molecule'].get_connectivity_matrix())
            entry_log['formed'] = [[a + 1, b + 1] for a, b in formed]
            entry_log['broken'] = [[a + 1, b + 1] for a, b in broken]
            if len(formed) == 0 or len(broken) == 0:
                entry_log['skipped'] = 'no bond formed and broken (conformational: torsion network)'
                self.ostream.print_info(f"reaction entries {na[2]} -> {nb[2]}: no bond exchange; left to the torsion network.")
                continue
            waypoint = None
            for kind, ts_guess, ri in nodes[pa + 1:pb]:
                if kind == 'transition' and waypoint is None:
                    waypoint = self._seed_transition_state(ts_guess, min_a, min_b, formed + broken, minima,
                                                           states_basis, dih_keys, log, entry_log)
            entry_log['waypoint'] = waypoint
            paths.append({'kind': 'reaction', 'A': na[1], 'B': nb[1], 'waypoint': waypoint, 'formed': formed,
                          'broken': broken, 'orphan': False, 'index': start_index + len(paths)})
        self.ostream.print_info(f'reaction edges: {len(paths)}')
        self.ostream.flush()

        return paths

    def _seed_transition_state(self, ts_guess, min_a, min_b, changed, minima, states_basis, dih_keys, log, entry_log):
        """
        Transition-state search from a transition reaction entry. The result
        is kept only if it has exactly one imaginary frequency below
        -20 cm^-1 and every exchanged bond lies between its lengths in the
        two end minima; it is then written as a datapoint without important
        coordinates and appended to the minima as a non-node saddle.

        :param ts_guess:
            The transition reaction entry (start of the search).
        :param min_a:
            The minimum dict of the reactant end.
        :param min_b:
            The minimum dict of the product end.
        :param changed:
            The exchanged bonds (0-based pairs).
        :param minima:
            The minima dicts (the saddle is appended).
        :param states_basis:
            The basis-set labels dict.
        :param dih_keys:
            The scanned dihedrals (1-based).
        :param log:
            The provenance dict.
        :param entry_log:
            The provenance dict of this reaction entry.

        :return:
            The index of the saddle in minima, or None.
        """

        root = self.roots_to_follow[0]
        try:
            opt_results = self._network_optimize(self._copy_molecule(ts_guess), [], states_basis, transition=True)
        except Exception as exc:
            entry_log['transition'] = {'status': f'TS search failed: {type(exc).__name__}: {exc}'}
            self.ostream.print_info(f'transition-state search failed ({type(exc).__name__}: {exc}); the half-scans '
                                    'meet at equal bond order.')
            self.ostream.flush()
            return None
        ts_mol = opt_results['final_molecule']
        basis = MolecularBasis.read(ts_mol, states_basis['gs'])
        energies, gradients, hessians = self._compute_qm_data(self.drivers['gs'], ts_mol, basis, with_hessian=True)
        self._network_counters['single_points'] += 1
        self._network_counters['hessians'] += 1
        masses = np.asarray(ts_mol.get_masses(), dtype=float)
        frequencies = self._lowest_frequency_cm1(hessians[0], masses, ts_mol.get_coordinates_in_bohr(), False,
                                                 return_all=True)
        imaginary = [float(f) for f in frequencies if f < -20.0]
        between = True
        for i, j in changed:
            d_ts = float(ts_mol.get_distance([i + 1, j + 1], 'angstrom'))
            d_a = float(min_a['molecule'].get_distance([i + 1, j + 1], 'angstrom'))
            d_b = float(min_b['molecule'].get_distance([i + 1, j + 1], 'angstrom'))
            if not (min(d_a, d_b) - 0.02 < d_ts < max(d_a, d_b) + 0.02):
                between = False
        entry_log['transition'] = {'imaginary_cm1': imaginary, 'between_end_minima': between,
                                   'energy': float(energies[0])}
        if len(imaginary) != 1 or not between:
            entry_log['transition']['status'] = 'rejected'
            self.ostream.print_info(f'transition-state search: {len(imaginary)} imaginary frequency(ies) '
                                    f'{[round(f, 1) for f in imaginary]}, between the end minima: {between}; '
                                    'rejected, the half-scans meet at equal bond order.')
            self.ostream.flush()
            return None

        inv_sqrt_masses = None
        if self.states_interpolation_settings[root]['use_mass_weight']:
            inv_sqrt_masses = 1.0 / np.sqrt(np.repeat(masses, 3))
        natoms = len(ts_mol.get_labels())
        dp = self._build_interpolation_datapoint(ts_mol, self.roots_z_matrix[root], self.states_interpolation_settings[root],
                                                 [], float(energies[0]), np.asarray(gradients[0]).reshape(natoms, 3),
                                                 hessians[0], inv_sqrt_masses)
        label = self._write_interpolation_datapoint(dp, root, self.states_interpolation_settings)
        saddle = {'conformer': None, 'molecule': ts_mol, 'energy': float(dp.energy), 'index': len(minima),
                  'label': label, 'angles': [float(ts_mol.get_dihedral_in_degrees(key)) % 360.0 for key in dih_keys],
                  'lowest_freq_cm1': float(frequencies.min()), 'node': False, 'kind': 'saddle'}
        minima.append(saddle)
        log['minima'].append({k: saddle[k] for k in ('index', 'conformer', 'energy', 'label', 'angles',
                                                     'lowest_freq_cm1', 'node', 'kind')})
        entry_log['transition']['status'] = f'accepted as {label}'
        self.ostream.print_info(f'transition state {saddle["index"]} added as {label} '
                                f'(imaginary frequency {imaginary[0]:.1f} cm^-1).')
        self.ostream.flush()

        return saddle['index']

    def _reaction_path_scan(self, path, minima, states_basis, settings, path_log):
        """
        Step-wise constrained scans of a reaction edge from both end minima:
        from A the first formed bond, from B the first broken bond is set
        reaction_scan_step_angstrom closer to its value in the other minimum
        at every step (geomeTRIC 'set distance', all other coordinates
        relaxed, each step starts from the previous point). A half stops one
        step before the transition state if there is one, otherwise at the
        first point where the driven bond has at least the Pauling bond order
        n = exp((r_bonded - r) / 0.3 Angstrom) of the opposing bond (r_bonded:
        the bond's length in the minimum where it exists), which the two
        halves then share. Only the exchanged bonds may change connectivity.

        :param path:
            The reaction path dict.
        :param minima:
            The minima dicts.
        :param states_basis:
            The basis-set labels dict.
        :param settings:
            The checked seed_network_settings.
        :param path_log:
            The provenance dict of this path.

        :return:
            The candidate list (the format of _relaxed_path_scan with the
            driven bond as 'imp' and the distance in Angstrom as 'angle').
        """

        step = float(settings['reaction_scan_step_angstrom'])
        a_min = minima[path['A']]
        b_min = minima[path['B']]
        waypoint = minima[path['waypoint']]['molecule'] if path.get('waypoint') is not None else None
        allowed = {tuple(sorted((int(i), int(j)))) for i, j in path['formed'] + path['broken']}
        zmat_bonds = self.roots_z_matrix[self.roots_to_follow[0]]['bonds']
        global_constraints = list(self.use_minimized_structures[1] or [])

        def _dist(mol, bond):
            return float(mol.get_distance([bond[0] + 1, bond[1] + 1], 'angstrom'))

        r_bonded = {}
        for bond in path['formed']:
            r_bonded[bond] = _dist(b_min['molecule'], bond)
        for bond in path['broken']:
            r_bonded[bond] = _dist(a_min['molecule'], bond)

        def _order(mol, bond):
            return float(np.exp((r_bonded[bond] - _dist(mol, bond)) / 0.3))

        def _imp(bond):
            return self._zmat_oriented((bond[0] + 1, bond[1] + 1), zmat_bonds)

        path_log['scan_points'] = []
        path_log['violations'] = []
        path_log['meeting'] = 'transition state' if waypoint is not None else 'equal Pauling bond order'
        if len(path['formed']) > 1 or len(path['broken']) > 1:
            self.ostream.print_info('several bonds change: only the first formed / broken bond is driven, '
                                    'the others relax.')

        candidates = []
        for half, start_min, drive, other in (('A', a_min, path['formed'][0], path['broken'][0]),
                                              ('B', b_min, path['broken'][0], path['formed'][0])):
            start = start_min['molecule']
            conn_ref = np.asarray(start.get_connectivity_matrix())
            r0 = _dist(start, drive)
            r_end = _dist(waypoint, drive) if waypoint is not None else r_bonded[drive]
            total = abs(r_end - r0)
            n_max = int(np.floor(total / step + 1.0e-9))
            if waypoint is not None and n_max * step > total - 0.5 * step:
                n_max -= 1
            sign = 1.0 if r_end > r0 else -1.0
            previous = start
            for k in range(1, n_max + 1):
                target = r0 + sign * k * step
                constraint = f'set distance {drive[0] + 1} {drive[1] + 1} {target:.4f}'
                try:
                    opt_results = self._network_optimize(previous, global_constraints + [constraint], states_basis)
                except Exception as exc:
                    path_log['violations'].append({'half': half, 'target': target,
                                                   'reason': f'optimisation failed: {type(exc).__name__}: {exc}'})
                    break
                mol = opt_results['final_molecule']
                energy = float(opt_results['opt_energies'][-1])
                actual = _dist(mol, drive)
                conn = np.asarray(mol.get_connectivity_matrix())
                changed = [(int(i), int(j)) for i, j in zip(*np.nonzero(np.triu(conn != conn_ref)))]
                connectivity_ok = all(tuple(sorted(c)) in allowed for c in changed)
                distance_ok = abs(actual - target) <= 0.02
                crossed = waypoint is None and _order(mol, drive) >= _order(mol, other)
                path_log['scan_points'].append({'half': half, 'driven': [drive[0] + 1, drive[1] + 1], 'target': target,
                                                'value': actual, 'opposing': _dist(mol, other), 'energy': energy,
                                                'connectivity_ok': connectivity_ok, 'distance_ok': distance_ok,
                                                'crossed': bool(crossed)})
                if not (connectivity_ok and distance_ok):
                    reason = 'connectivity' if not connectivity_ok else 'distance'
                    path_log['violations'].append({'half': half, 'target': target, 'reason': reason})
                    self.ostream.print_info(f'half {half}: continuity violated at {target:.3f} A ({reason}); half-scan stopped.')
                    self.ostream.flush()
                    break
                candidates.append({'half': half, 'target': target, 'angle': actual, 'unit': 'A', 'energy': energy,
                                   'molecule': mol, 'imp': _imp(drive), 'index': len(candidates)})
                previous = mol
                if crossed:
                    break

        self.ostream.print_info(f"{len(candidates)} scan point(s) on this reaction path "
                                f"({sum(1 for c in candidates if c['half'] == 'A')} from A, "
                                f"{sum(1 for c in candidates if c['half'] == 'B')} from B); meeting: {path_log['meeting']}")
        self.ostream.flush()

        return candidates

    def _grow_path(self, path, candidates, driver, minima, earlier, states_basis, settings, path_log):
        """
        Stage 3b: QM energy and gradient at every candidate (cached), the
        normalised acceptance score
        s = max(|E_MSI - E_QM| / (N e_thr), gRMS(G_MSI - G_QM) / g_thr)
        with the in-memory driver, then one trial point at a time: the
        untested candidate with the largest s gets a Hessian and becomes a
        trial datapoint (important coordinate: the scanned dihedral, radius
        use_opt_confidence_radius[2]) appended to the in-memory driver only.
        It is written as point_{n+1} if the path's maximum s and the summed s
        over the acceptance set (this path, earlier paths sharing an end
        minimum, the end minima) both decrease; otherwise it is removed from
        memory. The path stops when all s <= 1, after max_points_per_path
        accepted points or max_rejections_per_path consecutive rejections.

        :param path:
            The path dict.
        :param candidates:
            The candidate list of _relaxed_path_scan or _reaction_path_scan
            (QM data are cached in it).
        :param driver:
            The in-memory InterpolationDriver of the current database.
        :param minima:
            The minima dicts.
        :param earlier:
            The (end-minima set, candidates) tuples of earlier paths.
        :param states_basis:
            The basis-set labels dict.
        :param settings:
            The checked seed_network_settings.
        :param path_log:
            The provenance dict of this path.
        """

        root = self.roots_to_follow[0]
        drivers = self.drivers['gs']
        natoms = len(minima[path['A']]['molecule'].get_labels())
        e_thr = float(settings['e_thrsh_kcal_per_atom'] if settings['e_thrsh_kcal_per_atom'] is not None else self.energy_threshold)
        g_thr = float(settings['g_rmsd_thrsh'] if settings['g_rmsd_thrsh'] is not None else self.gradient_rmsd_thrsh)
        key_zmat = None
        if 'dihedral' in path_log:
            key_zmat = self._zmat_oriented(path_log['dihedral'], self.roots_z_matrix[root]['dihedrals'])

        path_log['iterations'] = []
        path_log['accepted'] = 0
        path_log['rejected'] = 0
        if len(candidates) == 0:
            path_log['status'] = 'no candidates'
            return

        for cand in candidates:
            if 'E_qm' not in cand:
                basis = MolecularBasis.read(cand['molecule'], states_basis['gs'])
                energies, gradients, _ = self._compute_qm_data(drivers, cand['molecule'], basis, with_hessian=False)
                cand['E_qm'] = float(energies[0])
                cand['G_qm'] = np.asarray(gradients[0], dtype=float).reshape(-1)
                cand['basis'] = basis
                self._network_counters['single_points'] += 1

        # the end minima must stay exact: their stored QM data are the reference
        sqrt_masses = np.sqrt(np.repeat(np.asarray(minima[path['A']]['molecule'].get_masses(), dtype=float), 3))
        end_entries = []
        for end in (path['A'], path['B']) + ((path['waypoint'],) if path.get('waypoint') is not None else ()):
            label = minima[end]['label']
            dp = None
            for point in driver.qm_data_points:
                if str(point.point_label).startswith(label + '_'):
                    dp = point
                    break
            if dp is None:
                continue
            grad = np.asarray(dp.gradient, dtype=float).reshape(-1)
            if self.states_interpolation_settings[root]['use_mass_weight']:
                grad = grad * sqrt_masses
            end_entries.append({'molecule': minima[end]['molecule'], 'E_qm': float(dp.energy), 'G_qm': grad, 'minimum': end})

        acceptance_set = list(candidates)
        for ends, cands in earlier:
            if ends & {path['A'], path['B']}:
                acceptance_set.extend(cands)
        acceptance_set.extend(end_entries)

        def _score(entry):
            driver.compute(entry['molecule'])
            e_msi = float(driver.impes_coordinate.energy)
            g_msi = np.asarray(driver.impes_coordinate.gradient, dtype=float).reshape(-1)
            de = abs(e_msi - entry['E_qm']) * hartree_in_kcalpermol() / (natoms * e_thr)
            grms = np.sqrt(np.mean((g_msi - entry['G_qm'])**2)) * hartree_in_kcalpermol() / bohr_in_angstrom() / g_thr
            return max(de, grms), de, grms

        def _metrics():
            scores = [_score(c)[0] for c in candidates]
            set_scores = [_score(c)[0] for c in acceptance_set]
            return scores, max(scores), float(np.sum(set_scores))

        tested = set()
        consecutive_rejections = 0
        status = 'stopped'
        while True:
            scores, path_max, set_sum = _metrics()
            if path_max <= 1.0:
                status = 'converged (all s <= 1)'
                break
            if path_log['accepted'] >= int(settings['max_points_per_path']):
                status = 'max_points_per_path reached'
                break
            if consecutive_rejections >= int(settings['max_rejections_per_path']):
                status = 'max_rejections_per_path reached'
                break
            untested = [i for i in range(len(candidates)) if i not in tested]
            if len(untested) == 0:
                status = 'all candidates tested'
                break
            pick = max(untested, key=lambda i: scores[i])
            cand = candidates[pick]
            s_pick, de_pick, grms_pick = _score(cand)

            if isinstance(drivers[0], (ScfRestrictedDriver, ScfUnrestrictedDriver)):
                # ScfHessianDriver.compute reuses scf_driver.scf_results without a geometry check,
                # and the last SCF above was the one of the last candidate: converge the SCF at
                # the picked candidate before its Hessian
                self._compute_energy(drivers[0], cand['molecule'], cand['basis'])
                self._network_counters['single_points'] += 1

            hessians = self._compute_hessian(drivers[2], cand['molecule'], cand['basis'])
            self._network_counters['hessians'] += 1
            inv_sqrt_masses = None
            if self.states_interpolation_settings[root]['use_mass_weight']:
                inv_sqrt_masses = 1.0 / np.sqrt(np.repeat(cand['molecule'].get_masses().copy(), 3))
            trial = self._build_interpolation_datapoint(
                cand['molecule'], self.roots_z_matrix[root], self.states_interpolation_settings[root],
                [c for c in [cand.get('imp', key_zmat)] if c is not None],
                cand['E_qm'], cand['G_qm'].reshape(natoms, 3), hessians[0], inv_sqrt_masses)
            driver.qm_data_points.append(trial)
            tested.add(pick)

            _, path_max_after, set_sum_after = _metrics()
            accepted = bool(path_max_after < path_max and set_sum_after < set_sum)
            # the candidate's index in the scan (as in scan_points), not its position after the window filter
            iteration = {'candidate': cand.get('index', pick), 'half': cand['half'], 'angle': cand['angle'], 's_before': s_pick,
                         'dE_norm': de_pick, 'gRMS_norm': grms_pick, 'path_max_before': path_max,
                         'set_sum_before': set_sum, 'path_max_after': path_max_after, 'set_sum_after': set_sum_after,
                         'accepted': accepted, 'label': None}
            if accepted:
                label = self._write_interpolation_datapoint(trial, root, self.states_interpolation_settings)
                iteration['label'] = label
                path_log['accepted'] += 1
                consecutive_rejections = 0
            else:
                driver.qm_data_points.pop()
                path_log['rejected'] += 1
                consecutive_rejections += 1
            path_log['iterations'].append(iteration)
            where = (f"{cand['angle']:.1f} deg" if cand.get('unit', 'deg') == 'deg'
                     else f"{cand['angle']:.3f} {cand['unit']}")
            self.ostream.print_info(
                f"  trial at {where} (half {cand['half']}): s = {s_pick:.3f}; path max "
                f"{path_max:.3f} -> {path_max_after:.3f}; set sum {set_sum:.3f} -> {set_sum_after:.3f}; "
                f"{'accepted as ' + str(iteration['label']) if accepted else 'rejected'}")
            self.ostream.flush()

        scores, path_max, set_sum = _metrics()
        path_log['status'] = status
        path_log['final_scores'] = [float(s) for s in scores]
        path_log['final_path_max'] = float(path_max)
        self.ostream.print_info(f"  path done: {status}; accepted {path_log['accepted']}, rejected "
                                f"{path_log['rejected']}, final max s = {path_max:.3f}")
        self.ostream.flush()

    def _load_interpolation_driver(self, molecule, root, im_settings, datafile):
        """
        Builds an InterpolationDriver with all datapoints of a database file
        in memory, set up as the confirmation step does (mass weighting,
        symmetry information, eq bond lengths from the first point).

        :param molecule:
            The molecule (for the masses).
        :param root:
            The electronic state.
        :param im_settings:
            The interpolation settings dicts per root.
        :param datafile:
            The database file.

        :return:
            The driver.
        """

        masses = molecule.get_masses().copy()
        masses_cart = np.repeat(masses, 3)
        inv_sqrt_masses = 1.0 / np.sqrt(masses_cart)

        impes_driver = InterpolationDriver(self.roots_z_matrix[root])
        impes_driver.update_settings(im_settings[root])
        impes_driver.impes_coordinate.inv_sqrt_masses = inv_sqrt_masses
        if root == 0:
            impes_driver.symmetry_information = self.symmetry_information['gs']
        else:
            impes_driver.symmetry_information = self.symmetry_information['es']

        im_labels, _ = impes_driver.read_labels()
        impes_driver.qm_data_points = []

        for label in im_labels:
            qm_data_point = InterpolationDatapoint(self.roots_z_matrix[root])
            qm_data_point.update_settings(im_settings[root])
            qm_data_point.read_hdf5(datafile, label)

            if impes_driver.impes_coordinate.eq_bond_lengths is None:
                impes_driver.impes_coordinate.eq_bond_lengths = qm_data_point.eq_bond_lengths

            impes_driver.qm_data_points.append(qm_data_point)

        return impes_driver

    def _compute_qm_data(self, drivers, molecule, basis, states=None, with_hessian=True):
        """
        Computes the QM energies, gradients and (optionally) Hessians of a
        molecule with a driver triple, as add_point does.

        :param drivers:
            The (energy, gradient, Hessian) drivers.
        :param molecule:
            The molecule.
        :param basis:
            The basis set.
        :param states:
            The states selected from a response solver (None: all).
        :param with_hessian:
            False to skip the Hessian (returned as None).

        :return:
            The energies, gradients and Hessians.
        """

        energies, scf_results, rsp_results = self._compute_energy(drivers[0], molecule, basis)

        if isinstance(drivers[0], LinearResponseEigenSolver) or isinstance(drivers[0], TdaEigenSolver):
            energies = energies[states]

        gradients = self._compute_gradient(drivers[1], molecule, basis, scf_results, rsp_results)
        hessians = self._compute_hessian(drivers[2], molecule, basis) if with_hessian else None

        return energies, gradients, hessians

    def _build_interpolation_datapoint(self, molecule, z_matrix, settings, constraints, energy, gradient, hessian,
                                       inv_sqrt_masses):
        """
        Builds the InterpolationDatapoint of a geometry from its QM data
        (mass weighting of the frozen settings, Cartesian to internal
        transformation, important coordinates from the 1-based constraint
        list, shifted to 0-based atom indices), confidence radius
        use_opt_confidence_radius[2]), as add_point does.

        :param molecule:
            The molecule.
        :param z_matrix:
            The z-matrix dict.
        :param settings:
            The interpolation settings dict of the root.
        :param constraints:
            The constraints (1-based atom-index tuples, in the orientation of
            the z-matrix) that become the important coordinates.
        :param energy:
            The QM energy.
        :param gradient:
            The Cartesian QM gradient.
        :param hessian:
            The Cartesian QM Hessian.
        :param inv_sqrt_masses:
            The inverse square roots of the masses (3N), or None.

        :return:
            The datapoint.
        """

        imp_int_constraints = {'bonds': [], 'angles': [], 'dihedrals': [], 'impropers': []}
        for constraint in constraints:
            constraint = [int(x) - 1 for x in constraint]
            if len(constraint) == 2:
                imp_int_constraints["bonds"].append(constraint)
            if len(constraint) == 3:
                imp_int_constraints["angles"].append(constraint)
            if len(constraint) == 4:
                imp_int_constraints["dihedrals"].append(constraint)

        grad = gradient.copy()
        hess = hessian.copy()
        grad_vec = grad.reshape(-1)         # (3N,)
        hess_mat = hess.reshape(grad_vec.size, grad_vec.size)
        mw_grad_vec = grad_vec
        mw_hess_mat = hess_mat
        if settings['use_mass_weight']:
            mw_grad_vec = inv_sqrt_masses * grad_vec
            mw_hess_mat = (inv_sqrt_masses[:, None] * hess_mat) * inv_sqrt_masses[None, :]

        impes_coordinate = InterpolationDatapoint(z_matrix)
        impes_coordinate.update_settings(settings)
        impes_coordinate.cartesian_coordinates = molecule.get_coordinates_in_bohr()
        impes_coordinate.eq_bond_lengths = self.eq_bond_length
        impes_coordinate.imp_int_coordinates = imp_int_constraints
        impes_coordinate.inv_sqrt_masses = inv_sqrt_masses
        impes_coordinate.energy = energy
        impes_coordinate.gradient = mw_grad_vec.reshape(grad.shape)
        impes_coordinate.hessian = mw_hess_mat.reshape(hess.shape)
        impes_coordinate.transform_gradient_and_hessian()

        trust_radius = self.use_opt_confidence_radius[2]
        impes_coordinate.confidence_radius = trust_radius

        return impes_coordinate

    def _write_interpolation_datapoint(self, impes_coordinate, target_root, interpolation_settings, label=None):
        """
        Writes a datapoint to the database of a root and to its _org copy
        (<database stem>_org.h5, a label already there is not written again),
        as add_point does. Without a label the next point_{n+1} label is used.

        :param impes_coordinate:
            The datapoint.
        :param target_root:
            The electronic state.
        :param interpolation_settings:
            The interpolation settings dicts per root.
        :param label:
            The label, or None.

        :return:
            The label.
        """

        target_file = interpolation_settings[target_root]['imforcefield_file']

        if label is None:
            interpolation_driver = InterpolationDriver(self.roots_z_matrix[target_root])
            interpolation_driver.update_settings(interpolation_settings[target_root])
            interpolation_driver.imforcefield_file = target_file
            sorted_labels = []
            if Path(target_file).exists():
                org_labels, _ = interpolation_driver.read_labels()
                labels = [lbl for lbl in org_labels if '_symmetry' not in lbl]
                sorted_labels = sorted(labels, key=lambda x: int(x.split('_')[1]))
            label = f'point_{len(sorted_labels) + 1}'

        impes_coordinate.write_hdf5(target_file, label)

        # copy of the datapoints written by the generator (seed stage), named after the database
        # (im_database_0.h5 -> im_database_0_org.h5); a copy left by an earlier run must not abort
        # the construction after the QM
        org_file = str(Path(target_file).with_name(Path(target_file).stem + '_org.h5'))
        org_labels = []
        if Path(org_file).exists():
            org_driver = InterpolationDriver(self.roots_z_matrix[target_root])
            org_driver.update_settings(interpolation_settings[target_root])
            org_driver.imforcefield_file = org_file
            org_labels, _ = org_driver.read_labels()
        if label in org_labels:
            self.ostream.print_warning(f'{org_file} already contains {label} from an earlier run; the copy '
                                       'of this datapoint is not written.')
            self.ostream.flush()
        else:
            impes_coordinate.write_hdf5(org_file, label)
        impes_coordinate.point_label = label

        return label

    def _build_opt_constraint_list(self, constraints, index_offset=0):
        """
        Converts constraints into geomeTRIC constraint lines: a tuple of 2, 3
        or 4 atom indices becomes 'freeze distance/angle/dihedral', a string
        is passed unchanged.

        :param constraints:
            The constraints (atom-index tuples or geomeTRIC strings).
        :param index_offset:
            Added to every atom index (0: the tuples are 1-based already, as
            everywhere in this class).

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

    def _is_same_minimum(self, molecule, candidates):
        """
        Checks whether an optimised structure is one of the candidates, i.e.
        the all-atom RMSD after superposition is below
        seed_duplicate_rmsd_thrsh (Angstrom). Always False if that is None.

        :param molecule:
            The optimised molecule.
        :param candidates:
            The molecules to compare with.

        :return:
            True if the molecule is the same minimum as a candidate.
        """

        if self.seed_duplicate_rmsd_thrsh is None:
            return False

        xyz = molecule.get_coordinates_in_angstrom()
        for candidate in candidates:
            rmsd, _, _ = svd_superimpose(xyz, candidate.get_coordinates_in_angstrom())
            if rmsd < self.seed_duplicate_rmsd_thrsh:
                return True

        return False

    def _place_scan_seed(self, base_molecules, dih_key, target_deg, tol_deg=2.0,
                         grid_step_deg=30.0, max_aux_bonds=3, max_aux_simultaneous=2):
        """
        Places a scan seed with dihedral dih_key at target_deg when a rigid
        rotation is blocked by a steric clash (Molecule.set_dihedral then
        stops early without error).

        A rigid rotation is tried from every base molecule first. Otherwise
        up to max_aux_simultaneous of the max_aux_bonds rotatable bonds
        closest to the scanned bond are rotated on a grid of grid_step_deg,
        and the geometry with unchanged connectivity and the largest clearance
        (shortest distance between atoms four or more bonds apart) is kept.

        :param base_molecules:
            The molecules to start from, preferred one first.
        :param dih_key:
            The scanned dihedral (1-based).
        :param target_deg:
            The target value of the dihedral in degree.
        :param tol_deg:
            The accepted deviation from the target in degree.
        :param grid_step_deg:
            The grid step of the auxiliary rotations in degree.
        :param max_aux_bonds:
            The number of nearest rotatable bonds used for relief.
        :param max_aux_simultaneous:
            The largest number of auxiliary bonds rotated at once.

        :return:
            The molecule and the list of (1-based auxiliary dihedral, change in
            degree), or (None, None) if no clash-free geometry was found.
        """

        from scipy.sparse.csgraph import shortest_path

        for base_mol in base_molecules:
            trial = self._copy_molecule(base_mol)
            trial.set_dihedral_in_degrees(dih_key, target_deg, verbose=False)
            if self._circular_diff_deg(trial.get_dihedral_in_degrees(dih_key), target_deg) <= tol_deg:
                return trial, []

        labels = base_molecules[0].get_labels()
        conn = base_molecules[0].get_connectivity_matrix()
        bond_dist = shortest_path(conn, unweighted=True)
        far_pairs = (bond_dist >= 4.0)
        central = tuple(sorted((dih_key[1] - 1, dih_key[2] - 1)))

        def _clearance(mol):
            if not np.any(far_pairs):
                return np.inf
            xyz = mol.get_coordinates_in_bohr()
            dist = np.linalg.norm(xyz[:, None, :] - xyz[None, :, :], axis=-1)
            return float(np.min(dist[far_pairs]))

        def _end_atom(center, exclude):
            nbrs = [int(x) for x in np.nonzero(conn[center])[0] if x != exclude]
            heavy = [x for x in nbrs if labels[x] != 'H']
            pool = heavy if heavy else nbrs
            return min(pool) if pool else None

        # auxiliary torsions around the rotatable bonds closest to the scanned bond
        aux_bonds = sorted(
            {tuple(sorted(bond)) for bond in (self.all_rotatable_bonds or [])} - {central},
            key=lambda bond: min(bond_dist[a, b] for a in bond for b in central))
        aux_dihedrals = []
        for b, c in aux_bonds:
            a, d = _end_atom(b, c), _end_atom(c, b)
            if a is not None and d is not None:
                aux_dihedrals.append((a + 1, b + 1, c + 1, d + 1))
        aux_dihedrals = aux_dihedrals[:max_aux_bonds]

        n_grid = int(round(360.0 / grid_step_deg))
        deltas = [k * grid_step_deg for k in range(1, n_grid)]

        for n_aux in range(1, max_aux_simultaneous + 1):
            best = None
            for base_mol in base_molecules:
                for subset in itertools.combinations(aux_dihedrals, n_aux):
                    for combo in itertools.product(deltas, repeat=n_aux):
                        trial = self._copy_molecule(base_mol)
                        applied = []
                        for aux, delta in zip(subset, combo):
                            phi_0 = trial.get_dihedral_in_degrees(aux)
                            trial.set_dihedral_in_degrees(aux, phi_0 + delta, verbose=False)
                            phi_1 = trial.get_dihedral_in_degrees(aux)
                            applied.append((aux, (phi_1 - phi_0 + 180.0) % 360.0 - 180.0))
                        trial.set_dihedral_in_degrees(dih_key, target_deg, verbose=False)
                        if self._circular_diff_deg(trial.get_dihedral_in_degrees(dih_key), target_deg) > tol_deg:
                            continue
                        if np.max(np.abs(trial.get_connectivity_matrix() - conn)) > 1.0e-10:
                            continue
                        score = (_clearance(trial), -sum(abs(x) for _, x in applied))
                        if best is None or score > best[0]:
                            best = (score, trial, applied)
            if best is not None:
                return best[1], best[2]

        return None, None

    def _run_optimization(self, optimization_driver, molecule, constraints=None, transition=False, index_offset=0,
                          compute_args=None):
        """
        Runs a (constrained) geomeTRIC optimisation or transition-state search
        with a muted OptimizationDriver.

        :param optimization_driver:
            The energy driver (SCF or xTB).
        :param molecule:
            The start molecule.
        :param constraints:
            The constraints (see _build_opt_constraint_list), or None.
        :param transition:
            True for a transition-state search.
        :param index_offset:
            The index offset of the tuple constraints.
        :param compute_args:
            Extra arguments of OptimizationDriver.compute (basis, SCF results
            for SCF drivers), or None.

        :return:
            The optimisation results dict.
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

        return opt_results

    def _compute_energy(self, qm_driver, molecule, basis=None):
        """
        Computes the QM energy with an xTB or SCF driver (muted).

        :param qm_driver:
            The energy driver.
        :param molecule:
            The molecule.
        :param basis:
            The basis set (SCF only).

        :return:
            The energies (array of one value), the SCF results (SCF only,
            else None) and None (no response results).
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
        Computes the QM gradient with an xTB or SCF gradient driver (muted).

        :param grad_driver:
            The gradient driver.
        :param molecule:
            The molecule.
        :param basis:
            The basis set (SCF only).
        :param scf_results:
            The SCF results at this geometry (SCF only).
        :param rsp_results:
            Not used.

        :return:
            The gradients, shape (1, natoms, 3).
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
        Computes the QM Hessian with an xTB or SCF Hessian driver (muted).
        The SCF Hessian driver reuses the SCF results stored in its SCF
        driver, so the last SCF must be at this geometry.

        :param hess_driver:
            The Hessian driver.
        :param molecule:
            The molecule.
        :param basis:
            The basis set (SCF only).

        :return:
            The Hessians, shape (1, 3 natoms, 3 natoms).
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
