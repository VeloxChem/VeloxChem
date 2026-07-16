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
import random
import h5py
import itertools
from mpi4py import MPI
from contextlib import redirect_stderr
from pathlib import Path
from .outputstream import OutputStream
from io import StringIO
from copy import deepcopy
from dataclasses import replace

from .xtbdriver import XtbDriver
from .xtbgradientdriver import XtbGradientDriver
from .xtbhessiandriver import XtbHessianDriver
from .gxtbdriver import GxtbDriver, GxtbGradientDriver, GxtbHessianDriver
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
from .imlocalgroups import LocalRotor, LocalCluster, LocalGroupModel
from .imlocalfactorregistry import (
    COMBINATION_SIGNED_RELAXED_RESIDUAL,
    COMBINATION_SIGNED_FULL,
    build_intersection_coefficients,
    load_signed_factor_banks_for_root,
    write_signed_factor_registry_for_family,
)
from .imcorefamilies import build_contact_torsion_v1_spec
from .immsiconstruction import Schema4MSIConstructor
from .mmforcefieldgenerator import MMForceFieldGenerator
from .conformergenerator import ConformerGenerator
from .optimizationdriver import OptimizationDriver
from .molecule import Molecule
from .errorhandler import assert_msg_critical
from .veloxchemlib import hartree_in_kcalpermol, bohr_in_angstrom


with redirect_stderr(StringIO()) as fg_err:
    import geometric


class IMForceFieldGenerator:
    """
    Class to set up and control the construction of the Interpolation Dynamics (IM) database.

    This class handles multiple aspects of database creation for force field generation, including
    molecule sampling, quantum mechanical calculations, and interpolation settings. It manages
    data such as molecular structures, energy calculations, and control parameters for dynamics.

    Instance Variables:

        - density_of_datapoints: Tracks the number of current interpolation data points in the database. Used to monitor
                                 database size and data density.

        - qm_data_points: Stores interpolation data points.

        - qmlabels: Labels assigned to each interpolation data point in the database.

        - molecules_along_rp: Represents molecular structures along a predefined reaction path or internal coordinate
                              pathway.

        - dihedrals: A list of dihedral angles to be rotated or scanned during simulations. Used to determine a pre-
                        defined path that should be sampled within the interpolation database construction.

        - sampling_structures: Specifies how many structures to generate around rotatable dihedral angles for database
                               population.

        - molecule: The initial molecular structure. (This mus be provided by the user)

        - datafile: Represents the database file (interpolation forcefield) used to store molecular structures and data points.
                    Typically initialized as `im_database.h5`.

        - z_matrix: The original Z-matrix (internal coordinate definition) for the molecule, specifying bonds, angles,
                    and dihedrals. It serves as the basis for internal coordinate transformations.

        - angle_threshold: Defines the range within which dihedral angles can vary during sampling and dynamics.
                           Making sure the smapling for 1 structure stays within a certain constained space.

        - interpolation_settings: A dictionary containing settings for the interpolation.

        - interpolation_type: Defines the type of interpolation used in the force field generation. The default value is
                              'shepard', likely referring to Shepard interpolation, a type of distance-weighted interpolation.

        - qm_driver, qm_grad_driver, qm_hess_driver: Instances of drivers for quantum mechanical (QM) calculations, including
                                                     single-point energy (qm_driver), gradient (qm_grad_driver), and Hessian
                                                     (qm_hess_driver) calculations. These drivers are used to perform QM tasks
                                                     during database construction. (Currently given from the user)

        - dynamics_settings: Holds various settings for molecular dynamics simulations, such as temperature, pressure,
                             force constants, and timestep.

        - basis_set_label: Specifies the basis set for QM calculations, with a default value of 'def2-svp'.

        - xcfun: Specifies the exchange-correlation functional for QM calculations. Default is 'b3lyp'.

        - duration: Specifies the total steps that can pass without a point being added to the database to determine
                    early breaks for a current structure.

        - temperature: Temperature for molecular dynamics simulations, defaulting to 150.15 K.

        - pressure: Pressure value for molecular dynamics, defaulting to 1.0 atm.

        - force_constant: The force constant used in simulations or sampling dynamics, with a default value of 1.0.

        - ensemble: Specifies the statistical ensemble used in the dynamics. Default: 'NVE'.

        - timestep: The time increment (in femtoseconds) between simulation steps. Default: 0.5 fs.

        - nsteps: Number of steps to be performed in the dynamics simulation.

        - snapshots: Specifies how many snapshots of the trajectory to record. Defaults: nsteps.

        - trajectory_file: Name of the file to store the trajectory of the molecular simulation. Default: 'trajectory.pdb'.

        - desired_point_density: Defines the desired density of data points for 1 structure in the database. Default: 50.

        - converged_cycle: Defines the number of cycles required for a simulation or database sampling to be considered
                           converged. Default: 4.

        - energy_threshold: Specifies an energy threshold to determine when a structure is necessary to be added into
                            the interpolation database. Default: 1.5 kcal/mol.

        - start_collect: Specifies at which step in the simulation interpolation datapoint collection should begin. Default: 0.

        - solvent: Specifies the solvent environment for the dynamics. Default: 'gas'. (Should always be gas for the construction)

        - qm_energies: A list to store QM energy from individual simulations or calculations (kj/mol).

        - total_energies: A list to store total energy (kj/mol).

        - molecules: A list to store molecular structures sampled during simulations (kj/mol).

        - kinetic_energies: A list to store kinetic energy from simulations (kj/mol).

        - point_added_molecules: A list of molecules for which new data points were added to the database.

        - unique_molecules: A list of unique molecular structures identified during database construction.

        - dynamics_method: Determines the method to generate molecular structures for the database quality conformation.

        - nstruc_to_confirm_database_quality: Number of randomly selected strucutures for the database quality check.
        
        - consider_locality: False -->  This key_word allows the constraint optimization part to allow more constraints to be considered to keep the optimized
                                        molecule closer to the current molecule. This might lead to more contamination of the interpolation database, should be 
                                        used when system is not improving without stricter locality constraints!   
    """

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

        self.density_of_datapoints = None
        self.qm_data_points = None
        self.qmlabels = None

        self.molecules_along_rp = None
        self.conformal_structures = None
        self.sampling_structures = 1

        self.state_specific_molecules = None
        self.datafile = None
        self.dihedrals_dict = None

        self.roots_z_matrix = {}
        self.int_coord_bond_information = None
        self.symmetry_dihedral_lists = {}
        self.symmetry_information = None
        self.local_group_primitive_model = {}
        self.local_group_model = {}
        self.local_group_coupled_model = {}
        self.local_groups = {}
        self.local_group_coupling_matrix = {}
        self.all_rotatable_bonds = None

        self.reaction_structures = None
        self.seed_structures = None

        # Here us the driver set up
        self.gs_basis_set_label = 'def2-svp'
        self.es_basis_set_label = '6-31g*'

        self.drivers = {'gs': None, 'es': None}
        self.sampling_driver = {'gs': None, 'es': None}

        if isinstance(ground_state_driver, ScfRestrictedDriver):
            # should be necessary to initialize

            qm_grad_driver = ScfGradientDriver(ground_state_driver)
            qm_hess_driver = ScfHessianDriver(ground_state_driver)

            self.drivers['gs'] = (ground_state_driver, qm_grad_driver, qm_hess_driver)

        if isinstance(ground_state_driver, ScfUnrestrictedDriver):
            # should be necessary to initialize

            qm_grad_driver = ScfGradientDriver(ground_state_driver)
            qm_hess_driver = ScfHessianDriver(ground_state_driver)

            self.drivers['gs'] = (ground_state_driver, qm_grad_driver, qm_hess_driver)
        # External Settings

        if isinstance(ground_state_driver, (XtbDriver, GxtbDriver)):
            if isinstance(ground_state_driver, GxtbDriver):
                qm_grad_driver = GxtbGradientDriver(ground_state_driver)
                qm_hess_driver = GxtbHessianDriver(ground_state_driver)
            else:
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
        # This depends on the current system: dihedral dominant small exponents wit q > p
        # strongly chaning electronic structure with bonds --> revers behavior with p >> q (8, 4)
        self.exponent_p = 3
        self.exponent_q = 4
        self.confidence_radius = 0.5
        self.imforcefieldfiles = None
        self.use_inverse_bond_length = True
        self.use_eq_bond_length = False
        self.eq_bond_symmetry_mode = "masked_exact"
        self.use_tc_weights = True
        self.tc_weight_mode = "multiplicative"  # "additive_rhee"
        self.use_mass_weight = True # set True as it is standard in the YM scheme --> small differences
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
        self.desired_point_density = 50
        self.converged_cycle = 5
        self.energy_threshold = 1.5
        self.gradient_rmsd_thrsh = 1.0
        self.distance_thrsh = 0.1
        self.start_collect = 0
        self.solvent = 'gas'
        self.add_bias_force = None
        self.bias_force_reaction_idx = None
        self.bias_force_reaction_prop = None  # this is being set by giving a dihedral, force constant and the final theta, steps_when_increased

        # sampling settings
        self.sampling_settings = {
            'enabled': False,
            'e_thrsh_kcal_per_atom': 0.1,
            'g_rmsd_thrsh_kcal_ang_per_atom':2.0,
            'force_orient_cos': 0.0001
        }

        self.sampling_imforcefieldfiles = None
        self.sampling_states_interpolation_settings = {root: None for root in roots_to_follow}
        self.metadynamics_settings = None
        # Example for the set up in the script/notebook
        # {
        #             "enabled": True,
        #             "bias_factor":10.0,
        #             "hill_height_kjmol":1.2,
        #             "hill_frequency":200,
        #             "variables":
        #             [
        #                 {
        #                     "type":"torsion",
        #                     "atoms":[3,4,6,10],
        #                     "min_deg":-180,
        #                     "max_deg":180,
        #                     "width_deg":12,
        #                     "periodic":True,
        #                 }
        #             ]
        # }

        # individual run variables, information used for database confirmation
        self.qm_energies = []
        self.total_energies = []
        self.molecules = None
        self.kinetic_energies = []
        self.point_added_molecules = []
        self.unique_molecules = []

        # In here I want to store Number_of_dp, exponent_p, exponent_q
        self.im_results = {'n_datapoints': None, 'RMSD': None, '|D|': None}

        # confirm database quality
        self.nstruc_to_confirm_database_quality = 15

        # set boolean for the optimization features of the metho
        self.use_minimized_structures = [True, [], []]  # to use minimum structures, specific constraints, root of the constraints
        self.add_conformal_structures = True  # Add all conformal structures which are being generated for the input molecule

        self.identfy_relevant_int_coordinates = True  # This constraint-optimizes new datapoints during construction run to presever optimal smoothness
        self.use_opt_confidence_radius = [False, 'multi_grad', 0.5, 0.3]  # optimize the set of Trust radii of each datapoint given a set of references (energy, gradient))
        self.exclude_non_core = False  # use only core atoms (H exclusion) for the interpolation (seemed to be less stable)

        self.imp_int_coordinates = []

        self._system_is_set_up = False

        self.local_group_specs = None
        self.auto_local_group_kinds = {
            "anchored_linker",
            "alkyl_chain",
            "alcohol",
            "methoxy",
            "primary_amine",
            "ammonium",
        }
        self.use_local_group_database = False
        self.local_group_phase_library = {
            3: (0.0, np.pi / 6.0, np.pi / 3.0, np.pi / 2.0),
        }
        self.local_group_coupling_threshold = 0.15
        self.local_group_coupling_topology = 'direct_factors'
        self.local_group_force_merge_overlapping_groups = True
        self.local_factor_state_weight_mode = "signature"
        self.local_factor_combination_mode = COMBINATION_SIGNED_FULL
        self.methyl_local_group_active_row_mode = "torsion_only"
        self.relax_local_factor_states = True
        self.local_factor_state_opt_max_iter = 120
        # Schema-3 factors are conditional residuals for one canonical rotor
        # subset.  Keep every rotor outside that subset at the family anchor;
        # otherwise a relaxed state can move another factor's signature and
        # activate both residuals when the training state is reconstructed.
        # Non-rotor environment coordinates remain free to respond.
        self.local_factor_anchor_external_rotors = True
        self.local_factor_detect_response_rows = False
        self.local_factor_response_thresholds = {
            "bond": 0.01,
            "angle": np.deg2rad(1.0),
            "dihedral": np.deg2rad(3.0),
            "improper": np.deg2rad(3.0),
        }
        self.relax_local_group_rotors_during_optimization = True
        self.relax_alkyl_local_group_states = True
        self.relax_methoxy_local_group_states = True
        self.methoxy_state_opt_max_iter = 500
        self.alkyl_chain_phase_values_degrees = (0, 60, 120, 180, 240, 300)
        self.anchored_linker_phase_values_degrees = (0, 60, 120, 180, 240, 300)
        self.anchored_linker_signature_angle_scale = 0.35
        self.anchored_linker_state_opt_max_iter = 500
        self.anchored_linker_relax_amide = True
        self.anchored_linker_freeze_boundary_angles = False
        self.local_group_allow_anchored_linker_methyl_coupling = False
        self.methoxy_phase_values_degrees = (0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340)

        # Schema 4 is opt-in.  Leaving this at 3 preserves every existing
        # construction and HDF5 path unless the user explicitly selects MSI.
        self.msi_schema_version = 3
        self.msi_descriptor_policy = {
            "policy_id": "contact_torsion_v1",
            "distance_coefficient": 1.0,
            "contact_r0_angstrom": 4.0,
            "contact_exponent": 6,
        }
        self.msi_locality_threshold_policy = {
            "policy_id": "schema4_locality_v1",
        }
        self.msi_promotion_policy = {
            "policy_id": "schema4_promotion_v1",
            "family_domain_radius": 0.50,
            "max_new_family_core_points": 3,
            "contact_regime_cutoff": 0.15,
        }
        self.msi_transfer_policy = {
            "policy_id": "schema4_transfer_v1",
            "require_target_probe": True,
        }
        self.msi_construction_stage_trace = ()

    @staticmethod
    def _is_xtb_like_driver(driver):
        return isinstance(driver, (XtbDriver, GxtbDriver))

    @staticmethod
    def _normalize_methyl_local_group_active_row_mode(mode):
        mode = "torsion_only" if mode is None else str(mode).strip().lower()
        mode = mode.replace("-", "_")

        aliases = {
            "": "torsion_only",
            "torsion": "torsion_only",
            "torsions": "torsion_only",
            "torsion_only": "torsion_only",
            "legacy_methyl": "torsion_only",
            "old_methyl": "torsion_only",
            "hydrogen": "hydrogen_internal",
            "hydrogens": "hydrogen_internal",
            "hydrogen_internal": "hydrogen_internal",
            "hydrogen_bond_angle": "hydrogen_internal",
            "hydrogen_bond_angles": "hydrogen_internal",
            "bond_angle": "hydrogen_internal",
            "bond_angles": "hydrogen_internal",
            "active_internal": "hydrogen_internal",
            "full": "full_internal",
            "full_internal": "full_internal",
            "full_methyl": "full_internal",
            "full_methyl_internal": "full_internal",
            "carbon_internal": "full_internal",
        }
        if mode not in aliases:
            raise ValueError(
                "Unsupported methyl_local_group_active_row_mode="
                f"'{mode}'. Use 'torsion_only', 'hydrogen_internal', "
                "or 'full_internal'."
            )
        return aliases[mode]

    def _read_basis_for_driver(self, molecule, driver, basis_label):
        if self._is_xtb_like_driver(driver):
            return None
        return MolecularBasis.read(molecule, basis_label)

    def _refresh_basis_for_driver(self, molecule, driver, basis):
        if self._is_xtb_like_driver(driver):
            return None
        if basis is not None and hasattr(basis, 'get_main_basis_label'):
            return MolecularBasis.read(molecule, basis.get_main_basis_label())
        return basis

    @staticmethod
    def _driver_metadata(driver):
        get_metadata = getattr(driver, 'get_metadata', None)
        if callable(get_metadata):
            return get_metadata()
        return getattr(driver, 'metadata', None)

    def _collect_backend_metadata(self, drivers):
        if drivers is None:
            return {}

        roles = ('energy', 'gradient', 'hessian')
        gxtb_metadata = {}
        for role, driver in zip(roles, drivers):
            if isinstance(driver, (GxtbDriver, GxtbGradientDriver,
                                   GxtbHessianDriver)):
                if role == 'energy' and isinstance(driver, GxtbDriver):
                    metadata = getattr(driver, 'energy_metadata', None)
                    if metadata is None:
                        metadata = self._driver_metadata(driver)
                else:
                    metadata = self._driver_metadata(driver)
                if metadata:
                    gxtb_metadata[role] = metadata

        if not gxtb_metadata:
            return {}

        return {
            'reference_backend': 'gxtb',
            'gxtb': gxtb_metadata,
        }

    def define_z_matrix_dict(self, molecule, add_coordinates=None):
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

    def _detect_methyl_and_tertbutyl_local_groups(
            self,
            molecule,
            ff_gen,
            z_matrix,
            rotatable_bonds_zero_based=None):
        """
        Detect methyl rotors and tert-butyl rotor clusters from the force-field
        topology and return a LocalGroupModel.

        The detector intentionally uses only graph/topology information here:
        the atom labels, the force-field connectivity matrix, the force-field
        rotatable-bond list, and the active interpolation z-matrix.

        Detected objects:
            - methyl rotor:
                axis = heavy_anchor - methyl_carbon
                equivalent units = three H atoms
            - tert-butyl parent rotor:
                axis = external_anchor - quaternary_carbon
                equivalent units = three methyl subtrees
            - primitive local clusters:
                one cluster per detected rotor. This intentionally keeps the
                tert-butyl parent rotation and the three child methyl rotations
                separate until the Hessian coupling pass decides whether they
                need to be merged.

        Atom indices are zero-based throughout, matching the z-matrix.
        """

        def _element(atom_index):
            raw_label = str(labels[int(atom_index)]).strip()
            letters = ''.join(ch for ch in raw_label if ch.isalpha())
            if not letters:
                return raw_label.capitalize()
            return letters[0].upper() + letters[1:].lower()

        def _is_h(atom_index):
            return _element(atom_index) == 'H'

        def _is_c(atom_index):
            return _element(atom_index) == 'C'

        def _neighbors(atom_index):
            return tuple(int(i) for i in np.flatnonzero(conn[int(atom_index)]))

        def _canonical_bond(atom_i, atom_j):
            return tuple(sorted((int(atom_i), int(atom_j))))

        def _component_on_side(start_atom, blocked_bond):
            blocked = frozenset(int(x) for x in blocked_bond)
            visited = set()
            stack = [int(start_atom)]
            while stack:
                atom = stack.pop()
                if atom in visited:
                    continue
                visited.add(atom)
                for neigh in _neighbors(atom):
                    if frozenset((atom, neigh)) == blocked:
                        continue
                    if neigh not in visited:
                        stack.append(neigh)
            return tuple(sorted(visited))

        def _axis_is_allowed(axis):
            if not rotatable_bond_set:
                return True
            return _canonical_bond(*axis) in rotatable_bond_set

        def _torsion_rows_for_axis(axis):
            axis_set = set(int(x) for x in axis)
            rows = []
            for row_idx, coord in enumerate(flat_z_matrix):
                if len(coord) != 4:
                    continue
                if set(int(x) for x in coord[1:3]) == axis_set:
                    rows.append(int(row_idx))
            return tuple(rows)

        def _phase_coordinate(torsion_rows):
            if not torsion_rows:
                return None
            coord = flat_z_matrix[int(torsion_rows[0])]
            return tuple(int(x) for x in coord)

        def _rows_for_local_atoms(local_atoms, active_atoms):
            local_set = set(int(x) for x in local_atoms)
            active_set = set(int(x) for x in active_atoms)
            rows = []
            for row_idx, coord in enumerate(flat_z_matrix):
                coord_atoms = set(int(x) for x in coord)
                if coord_atoms <= active_set and coord_atoms & local_set:
                    rows.append(int(row_idx))
            return tuple(rows)

        labels = molecule.get_labels()
        natoms = len(labels)

        if hasattr(ff_gen, 'connectivity_matrix') and ff_gen.connectivity_matrix is not None:
            conn = np.asarray(ff_gen.connectivity_matrix, dtype=int)
        else:
            conn = np.asarray(molecule.get_connectivity_matrix(), dtype=int)

        if rotatable_bonds_zero_based is None:
            rotatable_bonds_zero_based = [
                _canonical_bond(i - 1, j - 1)
                for i, j in getattr(ff_gen, 'rotatable_bonds', [])
            ]
        rotatable_bond_set = {
            _canonical_bond(i, j)
            for i, j in rotatable_bonds_zero_based
        }

        flat_z_matrix = InterpolationDatapoint.flatten_z_matrix(z_matrix)
        methyl_active_row_mode = (
            self._normalize_methyl_local_group_active_row_mode(
                getattr(
                    self,
                    "methyl_local_group_active_row_mode",
                    "torsion_only",
                )
            )
        )

        methyl_by_carbon = {}
        rotors = {}

        for atom_idx in range(natoms):
            if not _is_c(atom_idx):
                continue

            neigh = _neighbors(atom_idx)
            h_neighbors = tuple(sorted(n for n in neigh if _is_h(n)))
            heavy_neighbors = tuple(sorted(n for n in neigh if not _is_h(n)))

            if len(h_neighbors) != 3 or len(heavy_neighbors) != 1:
                continue

            anchor = int(heavy_neighbors[0])
            axis = (anchor, int(atom_idx))
            torsion_rows = _torsion_rows_for_axis(axis)
            if not torsion_rows:
                continue

            rotor_id = f"methyl_{atom_idx}"
            rotor = LocalRotor(
                rotor_id=rotor_id,
                kind='methyl',
                axis=tuple(int(x) for x in axis),
                symmetry_order=3,
                owned_atoms=tuple(sorted(int(atom) for atom in h_neighbors)),
                unit_atom_sets=tuple((int(h_atom),) for h_atom in h_neighbors),
                torsion_rows=torsion_rows,
                phase_coordinate=_phase_coordinate(torsion_rows),
            )

            methyl_by_carbon[int(atom_idx)] = {
                'anchor': anchor,
                'carbon': int(atom_idx),
                'hydrogens': h_neighbors,
                'rotor_id': rotor_id,
            }
            rotors[rotor_id] = rotor

        clusters = {}
        tertbutyl_child_methyl_rotors = set()

        for atom_idx in range(natoms):
            if not _is_c(atom_idx):
                continue

            neigh = _neighbors(atom_idx)
            h_neighbors = tuple(n for n in neigh if _is_h(n))
            heavy_neighbors = tuple(sorted(n for n in neigh if not _is_h(n)))

            if h_neighbors or len(heavy_neighbors) != 4:
                continue

            methyl_neighbors = tuple(
                n for n in heavy_neighbors
                if n in methyl_by_carbon
                and methyl_by_carbon[n]['anchor'] == atom_idx
            )
            external_neighbors = tuple(
                n for n in heavy_neighbors
                if n not in methyl_neighbors
            )

            if len(methyl_neighbors) != 3 or len(external_neighbors) != 1:
                continue

            external_anchor = int(external_neighbors[0])
            axis = (external_anchor, int(atom_idx))
            if not _axis_is_allowed(axis):
                continue

            torsion_rows = _torsion_rows_for_axis(axis)
            if not torsion_rows:
                continue

            unit_atom_sets = tuple(
                _component_on_side(methyl_carbon, (atom_idx, methyl_carbon))
                for methyl_carbon in methyl_neighbors
            )
            owned_atoms = _component_on_side(atom_idx, (external_anchor, atom_idx))

            parent_rotor_id = f"tertbutyl_parent_{atom_idx}"
            parent_rotor = LocalRotor(
                rotor_id=parent_rotor_id,
                kind='tertbutyl_parent',
                axis=tuple(int(x) for x in axis),
                symmetry_order=3,
                owned_atoms=owned_atoms,
                unit_atom_sets=unit_atom_sets,
                torsion_rows=torsion_rows,
                phase_coordinate=_phase_coordinate(torsion_rows),
            )
            rotors[parent_rotor_id] = parent_rotor

            child_rotor_ids = tuple(
                methyl_by_carbon[methyl_carbon]['rotor_id']
                for methyl_carbon in methyl_neighbors
            )
            tertbutyl_child_methyl_rotors.update(child_rotor_ids)

            active_atoms = tuple(sorted(set(owned_atoms) | {external_anchor}))
            active_rows = tuple(sorted(
                set(_rows_for_local_atoms(owned_atoms, active_atoms))
                | set(parent_rotor.torsion_rows)
            ))

            cluster_id = f"tertbutyl_parent_{atom_idx}"
            clusters[cluster_id] = LocalCluster(
                cluster_id=cluster_id,
                family_label='tertbutyl_parent',
                rotor_ids=(parent_rotor_id,),
                owned_atoms=owned_atoms,
                active_atoms=active_atoms,
                active_rows=active_rows,
            )

        for methyl_carbon, methyl_info in sorted(methyl_by_carbon.items()):
            rotor_id = methyl_info['rotor_id']

            rotor = rotors[rotor_id]
            active_seed = set(rotor.owned_atoms) | set(rotor.axis)
            if methyl_active_row_mode == "torsion_only":
                active_rows = set(int(row) for row in rotor.torsion_rows)
            else:
                row_local_atoms = set(rotor.owned_atoms)
                if methyl_active_row_mode == "full_internal":
                    row_local_atoms.add(int(methyl_info['carbon']))
                active_rows = set(
                    _rows_for_local_atoms(row_local_atoms, active_seed))
                active_rows.update(int(row) for row in rotor.torsion_rows)

            active_atom_set = set(active_seed)
            for row_idx in active_rows:
                active_atom_set.update(
                    int(atom) for atom in flat_z_matrix[int(row_idx)])
            active_atoms = tuple(sorted(active_atom_set))
            active_rows = tuple(sorted(active_rows))

            cluster_id = f"methyl_{methyl_carbon}"
            family_label = (
                'tertbutyl_methyl'
                if rotor_id in tertbutyl_child_methyl_rotors else 'methyl'
            )
            clusters[cluster_id] = LocalCluster(
                cluster_id=cluster_id,
                family_label=family_label,
                rotor_ids=(rotor_id,),
                owned_atoms=rotor.owned_atoms,
                active_atoms=active_atoms,
                active_rows=active_rows,
            )

        local_owned_atoms = set()
        local_active_rows = set()
        for cluster in clusters.values():
            local_owned_atoms.update(int(atom) for atom in cluster.owned_atoms)
            local_active_rows.update(int(row) for row in cluster.active_rows)

        all_atoms = set(range(natoms))
        core_atoms = tuple(sorted(all_atoms - local_owned_atoms))
        core_rows = tuple(
            row_idx for row_idx in range(len(flat_z_matrix))
            if row_idx not in local_active_rows
        )

        local_group_model = LocalGroupModel(
            enabled=bool(clusters),
            rotors=rotors,
            clusters=clusters,
            core_atoms=core_atoms,
            core_rows=core_rows,
        )

        return local_group_model

    def _detect_heteroatom_local_groups(
            self,
            molecule,
            ff_gen,
            z_matrix,
            local_group_model,
            rotatable_bonds_zero_based=None):

        enabled_kinds = set(getattr(self, "auto_local_group_kinds", ()))
        if not enabled_kinds:
            return local_group_model

        labels = tuple(str(x).strip() for x in molecule.get_labels())
        natoms = len(labels)
        conn = np.asarray(ff_gen.connectivity_matrix, dtype=int)
        atom_types = tuple(str(x).strip() for x in ff_gen.atom_types)
        equivalent_atoms = tuple(getattr(ff_gen, "equivalent_atoms", ()) or ())

        flat_z_matrix = InterpolationDatapoint.flatten_z_matrix(z_matrix)

        def element(atom):
            letters = "".join(x for x in labels[int(atom)] if x.isalpha())
            return letters.capitalize()

        def neighbors(atom):
            return tuple(int(x) for x in np.flatnonzero(conn[int(atom)]))

        def canonical_axis(axis):
            return tuple(sorted(int(x) for x in axis))

        if rotatable_bonds_zero_based is None:
            rotatable_bonds_zero_based = tuple(
                canonical_axis((int(a) - 1, int(b) - 1))
                for a, b in getattr(ff_gen, "rotatable_bonds", ())
            )

        rotatable_axes = {
            canonical_axis(axis) for axis in rotatable_bonds_zero_based
        }
        registered_axes = {
            canonical_axis(rotor.axis)
            for rotor in local_group_model.rotors.values()
        }

        def component_on_side(start, blocked_axis):
            blocked = frozenset(int(x) for x in blocked_axis)
            visited = set()
            stack = [int(start)]

            while stack:
                atom = stack.pop()
                if atom in visited:
                    continue
                visited.add(atom)

                for neighbor in neighbors(atom):
                    if frozenset((atom, neighbor)) == blocked:
                        continue
                    if neighbor not in visited:
                        stack.append(neighbor)

            return tuple(sorted(visited))

        def torsion_rows_for_axis(axis):
            axis_set = set(int(x) for x in axis)
            return tuple(
                row_idx
                for row_idx, coordinate in enumerate(flat_z_matrix)
                if (
                    len(coordinate) == 4
                    and set(int(x) for x in coordinate[1:3]) == axis_set
                )
            )

        def oriented_coordinate(coordinate, axis):
            coordinate = tuple(int(x) for x in coordinate)
            axis = tuple(int(x) for x in axis)

            if tuple(coordinate[1:3]) == axis:
                return coordinate
            if tuple(coordinate[2:0:-1]) == axis:
                return tuple(reversed(coordinate))
            return None

        def select_phase_coordinate(axis, torsion_rows, hydrogens):
            hydrogen_set = set(int(x) for x in hydrogens)
            fallback = None

            for row in torsion_rows:
                coordinate = oriented_coordinate(flat_z_matrix[row], axis)
                if coordinate is None:
                    continue
                if fallback is None:
                    fallback = coordinate
                if coordinate[3] in hydrogen_set:
                    return coordinate

            return fallback

        def hydrogens_are_equivalent(hydrogens):
            if len(hydrogens) <= 1:
                return True
            if len(equivalent_atoms) != natoms:
                return True
            return len({
                str(equivalent_atoms[int(atom)])
                for atom in hydrogens
            }) == 1

        def phase_values(kind, symmetry_order):
            if symmetry_order == 3:
                return (0.0, np.pi / 6.0, np.pi / 3.0, np.pi / 2.0)
            if symmetry_order == 2:
                return tuple(np.deg2rad(x) for x in range(0, 180, 30))
            return tuple(np.deg2rad(x) for x in range(0, 360, 30))

        for center in range(natoms):
            atom_type = atom_types[center]
            center_element = element(center)
            center_neighbors = neighbors(center)

            hydrogens = tuple(
                sorted(atom for atom in center_neighbors if element(atom) == "H")
            )
            heavy_atoms = tuple(
                sorted(atom for atom in center_neighbors if element(atom) != "H")
            )

            kind = None
            symmetry_order = 1

            if (
                center_element == "O"
                and atom_type == "oh"
                and len(hydrogens) == 1
                and len(heavy_atoms) == 1
                and "alcohol" in enabled_kinds
            ):
                kind = "alcohol"

            elif (
                center_element == "N"
                and atom_type == "n8"
                and len(hydrogens) == 2
                and len(heavy_atoms) == 1
                and "primary_amine" in enabled_kinds
            ):
                kind = "primary_amine"
                if hydrogens_are_equivalent(hydrogens):
                    symmetry_order = 2

            elif (
                center_element == "N"
                and atom_type == "nz"
                and len(hydrogens) == 3
                and len(heavy_atoms) == 1
                and "ammonium" in enabled_kinds
            ):
                kind = "ammonium"
                if hydrogens_are_equivalent(hydrogens):
                    symmetry_order = 3

            if kind is None:
                continue

            anchor = int(heavy_atoms[0])
            axis = (anchor, int(center))
            axis_key = canonical_axis(axis)

            if axis_key not in rotatable_axes:
                continue

            # A methyl and heteroatom rotor on the same bond describe the same
            # physical torsional degree of freedom. Keep the existing rotor.
            if axis_key in registered_axes:
                continue

            owned_atoms = component_on_side(center, axis)

            # Reject cyclic axes or malformed graph partitions.
            if anchor in owned_atoms:
                continue
            if center not in owned_atoms:
                continue
            if not set(hydrogens).issubset(set(owned_atoms)):
                continue

            torsion_rows = torsion_rows_for_axis(axis)
            if not torsion_rows:
                continue

            phase_coordinate = select_phase_coordinate(
                axis, torsion_rows, hydrogens)
            if phase_coordinate is None:
                continue

            if symmetry_order > 1:
                unit_atom_sets = tuple((int(atom),) for atom in hydrogens)
            else:
                unit_atom_sets = ()

            rotor_id = f"{kind}_{center}"

            rotor = LocalRotor(
                rotor_id=rotor_id,
                kind=kind,
                axis=axis,
                symmetry_order=symmetry_order,
                owned_atoms=owned_atoms,
                unit_atom_sets=unit_atom_sets,
                torsion_rows=torsion_rows,
                phase_coordinate=phase_coordinate,
                phase_values=phase_values(kind, symmetry_order),
            )

            active_seed = set(owned_atoms) | {anchor}
            active_rows = {
                row_idx
                for row_idx, coordinate in enumerate(flat_z_matrix)
                if (
                    set(int(x) for x in coordinate) <= active_seed
                    and set(int(x) for x in coordinate) & set(owned_atoms)
                )
            }
            active_rows.update(torsion_rows)

            active_atoms = set(active_seed)
            for row_idx in active_rows:
                active_atoms.update(
                    int(atom) for atom in flat_z_matrix[row_idx]
                )

            cluster = LocalCluster(
                cluster_id=rotor_id,
                family_label=kind,
                rotor_ids=(rotor_id,),
                owned_atoms=owned_atoms,
                active_atoms=tuple(sorted(active_atoms)),
                active_rows=tuple(sorted(active_rows)),
            )

            local_group_model.rotors[rotor_id] = rotor
            local_group_model.clusters[rotor_id] = cluster
            registered_axes.add(axis_key)

        return local_group_model

    def _detect_methoxy_local_groups(
            self,
            molecule,
            ff_gen,
            z_matrix,
            local_group_model,
            rotatable_bonds_zero_based=None):

        enabled_kinds = set(getattr(self, "auto_local_group_kinds", ()))
        if not ({"methoxy", "aryl_ether"} & enabled_kinds):
            return local_group_model

        labels = tuple(str(x).strip() for x in molecule.get_labels())
        natoms = len(labels)
        conn = np.asarray(ff_gen.connectivity_matrix, dtype=int)
        flat_z_matrix = InterpolationDatapoint.flatten_z_matrix(z_matrix)
        phase_values = tuple(
            np.deg2rad(float(x))
            for x in getattr(
                self,
                "methoxy_phase_values_degrees",
                getattr(
                    self,
                    "alkyl_chain_phase_values_degrees",
                    (0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340),
                ),
            )
        )

        def element(atom):
            letters = "".join(x for x in labels[int(atom)] if x.isalpha())
            return letters.capitalize()

        def neighbors(atom):
            return tuple(int(x) for x in np.flatnonzero(conn[int(atom)]))

        def canonical_axis(axis):
            return tuple(sorted(int(x) for x in axis))

        def component_on_side(start, blocked_axis):
            blocked = frozenset(int(x) for x in blocked_axis)
            visited = set()
            stack = [int(start)]

            while stack:
                atom = stack.pop()
                if atom in visited:
                    continue
                visited.add(atom)

                for neighbor in neighbors(atom):
                    if frozenset((atom, neighbor)) == blocked:
                        continue
                    if neighbor not in visited:
                        stack.append(neighbor)

            return tuple(sorted(visited))

        def torsion_rows_for_axis(axis):
            axis_set = set(int(x) for x in axis)
            return tuple(
                row_idx
                for row_idx, coordinate in enumerate(flat_z_matrix)
                if (
                    len(coordinate) == 4
                    and set(int(x) for x in coordinate[1:3]) == axis_set
                )
            )

        def oriented_coordinate(coordinate, axis):
            coordinate = tuple(int(x) for x in coordinate)
            axis = tuple(int(x) for x in axis)

            if tuple(coordinate[1:3]) == axis:
                return coordinate
            if tuple(coordinate[2:0:-1]) == axis:
                return tuple(reversed(coordinate))
            return None

        def select_phase_coordinate(axis, torsion_rows, moving_atoms):
            moving_atoms = set(int(x) for x in moving_atoms)
            fallback = None

            for row in torsion_rows:
                coordinate = oriented_coordinate(flat_z_matrix[row], axis)
                if coordinate is None:
                    continue
                if fallback is None:
                    fallback = coordinate
                if coordinate[3] in moving_atoms:
                    return coordinate

            return fallback

        def rows_for_local_atoms(local_atoms, active_atoms):
            local_set = set(int(x) for x in local_atoms)
            active_set = set(int(x) for x in active_atoms)
            rows = []
            for row_idx, coordinate in enumerate(flat_z_matrix):
                coord_atoms = set(int(x) for x in coordinate)
                if coord_atoms <= active_set and coord_atoms & local_set:
                    rows.append(int(row_idx))
            return tuple(rows)

        def is_methyl_carbon(atom):
            if element(atom) != "C":
                return False
            neigh = neighbors(atom)
            hydrogens = tuple(n for n in neigh if element(n) == "H")
            heavy = tuple(n for n in neigh if element(n) != "H")
            return len(hydrogens) == 3 and len(heavy) == 1

        if rotatable_bonds_zero_based is None:
            rotatable_bonds_zero_based = tuple(
                canonical_axis((int(a) - 1, int(b) - 1))
                for a, b in getattr(ff_gen, "rotatable_bonds", ())
            )

        rotatable_axes = {
            canonical_axis(axis) for axis in rotatable_bonds_zero_based
        }
        registered_axes = {
            canonical_axis(rotor.axis)
            for rotor in local_group_model.rotors.values()
        }

        for oxygen in range(natoms):
            if element(oxygen) != "O":
                continue

            oxygen_neighbors = neighbors(oxygen)
            hydrogens = tuple(
                atom for atom in oxygen_neighbors if element(atom) == "H"
            )
            heavy_neighbors = tuple(
                atom for atom in oxygen_neighbors if element(atom) != "H"
            )
            if hydrogens or len(heavy_neighbors) != 2:
                continue

            methyl_neighbors = tuple(
                atom for atom in heavy_neighbors if is_methyl_carbon(atom)
            )
            if not methyl_neighbors:
                continue

            for methyl_carbon in methyl_neighbors:
                anchors = tuple(
                    atom for atom in heavy_neighbors if atom != methyl_carbon
                )
                if len(anchors) != 1:
                    continue

                anchor = int(anchors[0])
                axis = (anchor, int(oxygen))
                axis_key = canonical_axis(axis)
                if axis_key not in rotatable_axes:
                    continue
                if axis_key in registered_axes:
                    continue

                owned_atoms = component_on_side(oxygen, axis)
                if anchor in owned_atoms:
                    continue
                if oxygen not in owned_atoms or methyl_carbon not in owned_atoms:
                    continue

                methyl_hydrogens = tuple(
                    atom for atom in neighbors(methyl_carbon)
                    if element(atom) == "H"
                )
                if not set(methyl_hydrogens).issubset(set(owned_atoms)):
                    continue

                torsion_rows = torsion_rows_for_axis(axis)
                if not torsion_rows:
                    continue

                phase_coordinate = select_phase_coordinate(
                    axis, torsion_rows, owned_atoms)
                if phase_coordinate is None:
                    continue

                rotor_id = f"methoxy_{oxygen}_{methyl_carbon}"
                rotor = LocalRotor(
                    rotor_id=rotor_id,
                    kind="methoxy",
                    axis=tuple(int(x) for x in axis),
                    symmetry_order=1,
                    owned_atoms=tuple(sorted(int(x) for x in owned_atoms)),
                    unit_atom_sets=(),
                    torsion_rows=torsion_rows,
                    phase_coordinate=phase_coordinate,
                    phase_values=phase_values,
                )

                active_seed = set(owned_atoms) | {anchor}
                active_rows = set(rows_for_local_atoms(owned_atoms, active_seed))
                active_rows.update(int(row) for row in torsion_rows)

                active_atoms = set(active_seed)
                for row_idx in active_rows:
                    active_atoms.update(
                        int(atom) for atom in flat_z_matrix[row_idx]
                    )

                cluster = LocalCluster(
                    cluster_id=rotor_id,
                    family_label="methoxy",
                    rotor_ids=(rotor_id,),
                    owned_atoms=tuple(sorted(int(x) for x in owned_atoms)),
                    active_atoms=tuple(sorted(active_atoms)),
                    active_rows=tuple(sorted(active_rows)),
                )

                local_group_model.rotors[rotor_id] = rotor
                local_group_model.clusters[rotor_id] = cluster
                registered_axes.add(axis_key)

        return local_group_model

    def _detect_amide_side_chain_local_groups(
            self,
            molecule,
            ff_gen,
            z_matrix,
            local_group_model,
            rotatable_bonds_zero_based=None):

        enabled_kinds = set(getattr(self, "auto_local_group_kinds", ()))
        if "amide_side_chain" not in enabled_kinds:
            return local_group_model

        labels = tuple(str(x).strip() for x in molecule.get_labels())
        natoms = len(labels)
        conn = np.asarray(ff_gen.connectivity_matrix, dtype=int)
        atom_types = tuple(str(x).strip().lower() for x in ff_gen.atom_types)
        flat_z_matrix = InterpolationDatapoint.flatten_z_matrix(z_matrix)
        phase_values = tuple(
            np.deg2rad(float(x))
            for x in getattr(
                self,
                "alkyl_chain_phase_values_degrees",
                (0, 60, 120, 180, 240, 300),
            )
        )

        def element(atom):
            letters = "".join(x for x in labels[int(atom)] if x.isalpha())
            return letters.capitalize()

        def neighbors(atom):
            return tuple(int(x) for x in np.flatnonzero(conn[int(atom)]))

        def canonical_axis(axis):
            return tuple(sorted(int(x) for x in axis))

        def is_sp3_carbon(atom):
            return (
                element(atom) == "C"
                and len(atom_types) == natoms
                and atom_types[int(atom)] == "c3"
            )

        def component_on_side(start, blocked_axis):
            blocked = frozenset(int(x) for x in blocked_axis)
            visited = set()
            stack = [int(start)]

            while stack:
                atom = stack.pop()
                if atom in visited:
                    continue
                visited.add(atom)

                for neighbor in neighbors(atom):
                    if frozenset((atom, neighbor)) == blocked:
                        continue
                    if neighbor not in visited:
                        stack.append(neighbor)

            return tuple(sorted(visited))

        def torsion_rows_for_axis(axis):
            axis_set = set(int(x) for x in axis)
            return tuple(
                row_idx
                for row_idx, coordinate in enumerate(flat_z_matrix)
                if (
                    len(coordinate) == 4
                    and set(int(x) for x in coordinate[1:3]) == axis_set
                )
            )

        def oriented_coordinate(coordinate, axis):
            coordinate = tuple(int(x) for x in coordinate)
            axis = tuple(int(x) for x in axis)

            if tuple(coordinate[1:3]) == axis:
                return coordinate
            if tuple(coordinate[2:0:-1]) == axis:
                return tuple(reversed(coordinate))
            return None

        def select_phase_coordinate(axis, torsion_rows, moving_atoms):
            moving_atoms = set(int(x) for x in moving_atoms)
            fallback = None

            for row in torsion_rows:
                coordinate = oriented_coordinate(flat_z_matrix[row], axis)
                if coordinate is None:
                    continue
                if fallback is None:
                    fallback = coordinate
                if coordinate[3] in moving_atoms:
                    return coordinate

            return fallback

        def rows_for_local_atoms(local_atoms, active_atoms):
            local_set = set(int(x) for x in local_atoms)
            active_set = set(int(x) for x in active_atoms)
            rows = []
            for row_idx, coordinate in enumerate(flat_z_matrix):
                coord_atoms = set(int(x) for x in coordinate)
                if coord_atoms <= active_set and coord_atoms & local_set:
                    rows.append(int(row_idx))
            return tuple(rows)

        def sp3_chain_from_amide_n(first_atom, amide_n):
            chain = set()
            stack = [int(first_atom)]
            while stack:
                atom = stack.pop()
                if atom in chain:
                    continue
                if not is_sp3_carbon(atom):
                    continue
                chain.add(atom)

                for neighbor in neighbors(atom):
                    if neighbor == int(amide_n):
                        continue
                    if neighbor in chain:
                        continue
                    if is_sp3_carbon(neighbor):
                        stack.append(neighbor)

            return tuple(sorted(chain))

        def hydrogens_attached_to_atoms(atoms):
            hydrogens = set()
            for atom in atoms:
                hydrogens.update(
                    int(neighbor)
                    for neighbor in neighbors(atom)
                    if element(neighbor) == "H"
                )
            return hydrogens

        if rotatable_bonds_zero_based is None:
            rotatable_bonds_zero_based = tuple(
                canonical_axis((int(a) - 1, int(b) - 1))
                for a, b in getattr(ff_gen, "rotatable_bonds", ())
            )

        rotatable_axes = {
            canonical_axis(axis) for axis in rotatable_bonds_zero_based
        }
        registered_axes = {
            canonical_axis(rotor.axis)
            for rotor in local_group_model.rotors.values()
        }

        amide_records = []
        for carbonyl_c in range(natoms):
            if element(carbonyl_c) != "C":
                continue

            c_neighbors = neighbors(carbonyl_c)
            oxygen_neighbors = tuple(
                atom for atom in c_neighbors
                if element(atom) == "O"
            )
            nitrogen_neighbors = tuple(
                atom for atom in c_neighbors
                if element(atom) == "N"
            )
            if not oxygen_neighbors or not nitrogen_neighbors:
                continue

            for amide_n in nitrogen_neighbors:
                n_substituents = tuple(
                    atom for atom in neighbors(amide_n)
                    if atom != carbonyl_c and element(atom) != "H"
                )
                if not n_substituents:
                    continue

                for first_atom in n_substituents:
                    if not is_sp3_carbon(first_atom):
                        continue

                    chain_atoms = sp3_chain_from_amide_n(first_atom, amide_n)
                    if not chain_atoms:
                        continue

                    chain_set = set(chain_atoms)
                    boundary_axes = []
                    for chain_atom in chain_atoms:
                        for neighbor in neighbors(chain_atom):
                            if neighbor == amide_n:
                                continue
                            if neighbor in chain_set:
                                continue
                            if element(neighbor) == "H":
                                continue
                            boundary_axes.append((int(neighbor), int(chain_atom)))

                    if boundary_axes:
                        boundary_axis = sorted(
                            boundary_axes,
                            key=lambda axis: (
                                canonical_axis(axis) not in rotatable_axes,
                                axis,
                            ),
                        )[0]
                        owned_atoms = component_on_side(
                            boundary_axis[1], boundary_axis)
                    else:
                        boundary_axis = (int(amide_n), int(first_atom))
                        owned_atoms = component_on_side(first_atom, boundary_axis)

                    raw_owned_atoms = set(int(x) for x in owned_atoms)
                    if int(first_atom) not in raw_owned_atoms:
                        continue

                    amide_core_atoms = {
                        int(carbonyl_c),
                        int(amide_n),
                        *(int(atom) for atom in oxygen_neighbors),
                        *(
                            int(atom)
                            for atom in neighbors(amide_n)
                            if element(atom) == "H"
                        ),
                    }
                    tail_owned_atoms = set(int(atom) for atom in chain_atoms)
                    tail_owned_atoms.update(
                        hydrogens_attached_to_atoms(tail_owned_atoms)
                    )
                    tail_owned_atoms &= raw_owned_atoms
                    tail_owned_atoms -= amide_core_atoms
                    if not tail_owned_atoms:
                        continue

                    amide_records.append({
                        "carbonyl_c": int(carbonyl_c),
                        "amide_n": int(amide_n),
                        "first_atom": int(first_atom),
                        "chain_atoms": chain_atoms,
                        "boundary_axis": tuple(int(x) for x in boundary_axis),
                        "amide_core_atoms": tuple(sorted(amide_core_atoms)),
                        "owned_atoms": tuple(sorted(tail_owned_atoms)),
                    })

        for record in sorted(
                amide_records,
                key=lambda item: (-len(item["owned_atoms"]),
                                  item["carbonyl_c"],
                                  item["amide_n"],
                                  item["first_atom"])):
            cluster_id = (
                f"amide_side_chain_{record['carbonyl_c']}_"
                f"{record['amide_n']}_{record['first_atom']}"
            )
            if cluster_id in local_group_model.clusters:
                continue

            owned_atoms = tuple(int(x) for x in record["owned_atoms"])
            owned_set = set(owned_atoms)
            amide_core_atoms = set(
                int(x) for x in record.get("amide_core_atoms", ())
            )
            active_seed = (
                set(owned_atoms)
                | set(record["boundary_axis"])
                | amide_core_atoms
            )
            amide_axis = canonical_axis(
                (record["carbonyl_c"], record["amide_n"]))
            chain_set = set(int(x) for x in record["chain_atoms"]) & owned_set
            if not chain_set:
                continue

            rotor_ids = []
            torsion_rows = set()
            candidate_axes = set()

            n_axis = canonical_axis((record["amide_n"], record["first_atom"]))
            if n_axis in rotatable_axes:
                candidate_axes.add(tuple(record["boundary_axis"]))
                candidate_axes.add((record["amide_n"], record["first_atom"]))

            for axis_key in rotatable_axes:
                if axis_key == amide_axis:
                    continue
                axis_set = set(axis_key)
                if not axis_set <= active_seed:
                    continue
                if not (axis_set & chain_set):
                    continue
                if axis_key in registered_axes:
                    continue
                candidate_axes.add(tuple(axis_key))

            for axis_key in sorted(canonical_axis(axis) for axis in candidate_axes):
                if axis_key == amide_axis:
                    continue
                if axis_key in registered_axes:
                    continue

                atom_a, atom_b = axis_key
                if atom_a in chain_set and atom_b not in chain_set:
                    oriented_axis = (atom_b, atom_a)
                elif atom_b in chain_set and atom_a not in chain_set:
                    oriented_axis = (atom_a, atom_b)
                elif atom_a == record["amide_n"] and atom_b in chain_set:
                    oriented_axis = (atom_a, atom_b)
                elif atom_b == record["amide_n"] and atom_a in chain_set:
                    oriented_axis = (atom_b, atom_a)
                else:
                    oriented_axis = axis_key

                axis_torsion_rows = torsion_rows_for_axis(oriented_axis)
                if not axis_torsion_rows:
                    continue

                moving_atoms = component_on_side(oriented_axis[1], oriented_axis)
                if not set(int(x) for x in moving_atoms) <= owned_set:
                    continue
                phase_coordinate = select_phase_coordinate(
                    oriented_axis, axis_torsion_rows, moving_atoms)
                if phase_coordinate is None:
                    continue

                rotor_id = (
                    f"{cluster_id}_axis_{oriented_axis[0]}_"
                    f"{oriented_axis[1]}"
                )
                rotor = LocalRotor(
                    rotor_id=rotor_id,
                    kind="amide_side_chain",
                    axis=tuple(int(x) for x in oriented_axis),
                    symmetry_order=1,
                    owned_atoms=tuple(sorted(int(x) for x in moving_atoms)),
                    unit_atom_sets=(),
                    torsion_rows=axis_torsion_rows,
                    phase_coordinate=phase_coordinate,
                    phase_values=phase_values,
                )

                local_group_model.rotors[rotor_id] = rotor
                registered_axes.add(canonical_axis(oriented_axis))
                rotor_ids.append(rotor_id)
                torsion_rows.update(int(row) for row in axis_torsion_rows)

            if not rotor_ids:
                continue

            active_rows = set(rows_for_local_atoms(owned_atoms, active_seed))
            active_rows.update(torsion_rows)

            active_atoms = set(active_seed)
            for row_idx in active_rows:
                active_atoms.update(
                    int(atom) for atom in flat_z_matrix[row_idx]
                )

            local_group_model.clusters[cluster_id] = LocalCluster(
                cluster_id=cluster_id,
                family_label="amide_side_chain",
                rotor_ids=tuple(rotor_ids),
                owned_atoms=owned_atoms,
                active_atoms=tuple(sorted(active_atoms)),
                active_rows=tuple(sorted(active_rows)),
            )

        return local_group_model

    def _detect_anchored_linker_local_groups(
            self,
            molecule,
            ff_gen,
            z_matrix,
            local_group_model,
            rotatable_bonds_zero_based=None):

        enabled_kinds = set(getattr(self, "auto_local_group_kinds", ()))
        if "anchored_linker" not in enabled_kinds:
            return local_group_model

        labels = tuple(str(x).strip() for x in molecule.get_labels())
        natoms = len(labels)
        conn = np.asarray(ff_gen.connectivity_matrix, dtype=int)
        atom_types = tuple(str(x).strip().lower() for x in ff_gen.atom_types)
        flat_z_matrix = InterpolationDatapoint.flatten_z_matrix(z_matrix)
        phase_values = tuple(
            np.deg2rad(float(x))
            for x in getattr(
                self,
                "anchored_linker_phase_values_degrees",
                getattr(
                    self,
                    "alkyl_chain_phase_values_degrees",
                    (0, 60, 120, 180, 240, 300),
                ),
            )
        )
        angle_scale = float(
            getattr(self, "anchored_linker_signature_angle_scale", 0.35)
        )

        def element(atom):
            letters = "".join(x for x in labels[int(atom)] if x.isalpha())
            return letters.capitalize()

        def neighbors(atom):
            return tuple(int(x) for x in np.flatnonzero(conn[int(atom)]))

        def canonical_axis(axis):
            return tuple(sorted(int(x) for x in axis))

        def is_sp3_carbon(atom):
            return (
                element(atom) == "C"
                and len(atom_types) == natoms
                and atom_types[int(atom)] == "c3"
            )

        def hydrogens_attached_to_atoms(atoms):
            hydrogens = set()
            for atom in atoms:
                hydrogens.update(
                    int(neighbor)
                    for neighbor in neighbors(atom)
                    if element(neighbor) == "H"
                )
            return hydrogens

        def component_on_side(start_atom, blocked_axis):
            blocked = frozenset(int(atom) for atom in blocked_axis)
            visited = set()
            stack = [int(start_atom)]

            while stack:
                atom = stack.pop()
                if atom in visited:
                    continue
                visited.add(atom)

                for neighbor in neighbors(atom):
                    if frozenset((atom, neighbor)) == blocked:
                        continue
                    if neighbor not in visited:
                        stack.append(neighbor)

            return tuple(sorted(visited))

        def sp3_component_from(first_atom, blocked_atom):
            component = set()
            stack = [int(first_atom)]
            while stack:
                atom = stack.pop()
                if atom in component:
                    continue
                if not is_sp3_carbon(atom):
                    continue
                component.add(atom)

                for neighbor in neighbors(atom):
                    if neighbor == int(blocked_atom):
                        continue
                    if neighbor in component:
                        continue
                    if is_sp3_carbon(neighbor):
                        stack.append(neighbor)

            return tuple(sorted(component))

        def shortest_path_within(start, target, allowed_atoms):
            allowed = set(int(atom) for atom in allowed_atoms)
            start = int(start)
            target = int(target)
            queue = [(start, (start,))]
            seen = {start}

            while queue:
                atom, path = queue.pop(0)
                if atom == target:
                    return path

                for neighbor in neighbors(atom):
                    if neighbor not in allowed or neighbor in seen:
                        continue
                    seen.add(neighbor)
                    queue.append((neighbor, path + (neighbor,)))

            return ()

        def torsion_rows_for_axis(axis):
            axis_set = set(int(x) for x in axis)
            return tuple(
                row_idx
                for row_idx, coordinate in enumerate(flat_z_matrix)
                if (
                    len(coordinate) == 4
                    and set(int(x) for x in coordinate[1:3]) == axis_set
                )
            )

        def angle_rows_for_angle(angle):
            angle = tuple(int(x) for x in angle)
            rev_angle = tuple(reversed(angle))
            return tuple(
                row_idx
                for row_idx, coordinate in enumerate(flat_z_matrix)
                if (
                    len(coordinate) == 3
                    and (
                        tuple(int(x) for x in coordinate) == angle
                        or tuple(int(x) for x in coordinate) == rev_angle
                    )
                )
            )

        def oriented_coordinate(coordinate, axis):
            coordinate = tuple(int(x) for x in coordinate)
            axis = tuple(int(x) for x in axis)

            if tuple(coordinate[1:3]) == axis:
                return coordinate
            if tuple(coordinate[2:0:-1]) == axis:
                return tuple(reversed(coordinate))
            return None

        def select_phase_coordinate(axis, torsion_rows, left_atoms, right_atoms):
            left_atoms = set(int(x) for x in left_atoms)
            right_atoms = set(int(x) for x in right_atoms)
            fallback = None

            for row in torsion_rows:
                coordinate = oriented_coordinate(flat_z_matrix[row], axis)
                if coordinate is None:
                    continue
                if fallback is None:
                    fallback = coordinate
                if coordinate[0] in left_atoms and coordinate[3] in right_atoms:
                    return coordinate

            return fallback

        def rows_for_local_atoms(local_atoms, active_atoms):
            local_set = set(int(x) for x in local_atoms)
            active_set = set(int(x) for x in active_atoms)
            rows = []
            for row_idx, coordinate in enumerate(flat_z_matrix):
                coord_atoms = set(int(x) for x in coordinate)
                if coord_atoms <= active_set and coord_atoms & local_set:
                    rows.append(int(row_idx))
            return tuple(rows)

        def add_signature_rows(rows, row_types, row_scales, new_rows,
                               row_type, row_scale):
            seen = set(rows)
            for row in new_rows:
                row = int(row)
                if row in seen:
                    continue
                seen.add(row)
                rows.append(row)
                row_types.append(str(row_type))
                row_scales.append(float(row_scale))

        if rotatable_bonds_zero_based is None:
            rotatable_bonds_zero_based = tuple(
                canonical_axis((int(a) - 1, int(b) - 1))
                for a, b in getattr(ff_gen, "rotatable_bonds", ())
            )

        rotatable_axes = {
            canonical_axis(axis) for axis in rotatable_bonds_zero_based
        }
        registered_axes = {
            canonical_axis(rotor.axis)
            for rotor in local_group_model.rotors.values()
        }

        records = []
        for carbonyl_c in range(natoms):
            if element(carbonyl_c) != "C":
                continue

            c_neighbors = neighbors(carbonyl_c)
            oxygen_neighbors = tuple(
                atom for atom in c_neighbors
                if element(atom) == "O"
            )
            nitrogen_neighbors = tuple(
                atom for atom in c_neighbors
                if element(atom) == "N"
            )
            if not oxygen_neighbors or not nitrogen_neighbors:
                continue

            for amide_n in nitrogen_neighbors:
                n_substituents = tuple(
                    atom for atom in neighbors(amide_n)
                    if atom != carbonyl_c and element(atom) != "H"
                )
                for first_atom in n_substituents:
                    if not is_sp3_carbon(first_atom):
                        continue

                    chain_atoms = sp3_component_from(first_atom, amide_n)
                    if len(chain_atoms) < 2:
                        continue

                    chain_set = set(chain_atoms)
                    boundary_axes = []
                    for chain_atom in chain_atoms:
                        for neighbor in neighbors(chain_atom):
                            if neighbor == amide_n:
                                continue
                            if neighbor in chain_set:
                                continue
                            if element(neighbor) == "H":
                                continue
                            if is_sp3_carbon(neighbor):
                                continue
                            boundary_axes.append((int(neighbor), int(chain_atom)))

                    for ring_anchor, boundary_chain_atom in boundary_axes:
                        path = shortest_path_within(
                            first_atom, boundary_chain_atom, chain_atoms)
                        if len(path) < 2:
                            continue

                        owned_atoms = set(int(atom) for atom in path)
                        owned_atoms.update(hydrogens_attached_to_atoms(owned_atoms))

                        amide_anchor_atoms = {
                            int(carbonyl_c),
                            int(amide_n),
                            *(int(atom) for atom in oxygen_neighbors),
                            *(
                                int(atom)
                                for atom in neighbors(amide_n)
                                if element(atom) == "H"
                            ),
                        }
                        amide_relaxation_atoms = set(amide_anchor_atoms)
                        for carbonyl_neighbor in c_neighbors:
                            if carbonyl_neighbor == int(amide_n):
                                continue
                            if int(carbonyl_neighbor) in oxygen_neighbors:
                                continue
                            amide_relaxation_atoms.update(
                                component_on_side(
                                    carbonyl_neighbor,
                                    (carbonyl_c, carbonyl_neighbor),
                                )
                            )
                        amide_relaxation_atoms -= owned_atoms
                        ring_anchor_atoms = {int(ring_anchor)}
                        ring_anchor_atoms.update(
                            int(atom)
                            for atom in neighbors(ring_anchor)
                            if atom != int(boundary_chain_atom)
                            and element(atom) != "H"
                        )
                        anchor_atoms = tuple(
                            sorted(amide_anchor_atoms | ring_anchor_atoms)
                        )

                        records.append({
                            "carbonyl_c": int(carbonyl_c),
                            "oxygen_neighbors": tuple(
                                int(atom) for atom in oxygen_neighbors),
                            "amide_n": int(amide_n),
                            "first_atom": int(first_atom),
                            "ring_anchor": int(ring_anchor),
                            "path": tuple(int(atom) for atom in path),
                            "owned_atoms": tuple(sorted(owned_atoms)),
                            "anchor_atoms": anchor_atoms,
                            "relaxation_atoms": tuple(
                                sorted(amide_relaxation_atoms)
                            ),
                        })

        records = sorted(
            records,
            key=lambda item: (
                len(item["path"]),
                item["carbonyl_c"],
                item["amide_n"],
                item["ring_anchor"],
            ),
            reverse=True,
        )

        claimed_owned_atoms = set(
            int(atom)
            for cluster in local_group_model.clusters.values()
            for atom in cluster.owned_atoms
        )

        for record in records:
            owned_atoms = tuple(int(atom) for atom in record["owned_atoms"])
            owned_set = set(owned_atoms)
            if owned_set & claimed_owned_atoms:
                continue

            path = tuple(int(atom) for atom in record["path"])
            cluster_id = (
                f"anchored_linker_{record['amide_n']}_"
                f"{path[0]}_{path[-1]}_{record['ring_anchor']}"
            )
            if cluster_id in local_group_model.clusters:
                continue

            anchor_atoms = tuple(int(atom) for atom in record["anchor_atoms"])
            anchor_set = set(anchor_atoms)
            relaxation_atoms = tuple(
                int(atom) for atom in record.get("relaxation_atoms", ())
            )
            active_seed = owned_set | anchor_set | set(relaxation_atoms)
            rotor_ids = []
            active_signature_rows = set()

            candidate_axes = [(record["amide_n"], path[0])]
            candidate_axes.extend(
                (path[idx], path[idx + 1])
                for idx in range(len(path) - 1)
            )

            for axis_index, axis in enumerate(candidate_axes):
                axis = tuple(int(atom) for atom in axis)
                axis_key = canonical_axis(axis)
                if axis_key not in rotatable_axes:
                    continue
                if axis_key in registered_axes:
                    continue

                axis_torsion_rows = torsion_rows_for_axis(axis)
                if not axis_torsion_rows:
                    continue

                if axis_index == 0:
                    left_atoms = (
                        record["carbonyl_c"],
                        *record["oxygen_neighbors"],
                    )
                    right_atoms = path[1:] or (record["ring_anchor"],)
                else:
                    left_atoms = (record["amide_n"], *path[:axis_index])
                    right_atoms = (
                        *path[axis_index + 1:],
                        record["ring_anchor"],
                    )

                phase_coordinate = select_phase_coordinate(
                    axis, axis_torsion_rows, left_atoms, right_atoms)
                if phase_coordinate is None:
                    continue

                signature_rows = []
                signature_row_types = []
                signature_row_scales = []
                add_signature_rows(
                    signature_rows,
                    signature_row_types,
                    signature_row_scales,
                    axis_torsion_rows,
                    "torsion",
                    1.0,
                )

                if axis_index == 0:
                    angle_candidates = [
                        (record["carbonyl_c"], record["amide_n"], path[0]),
                    ]
                    if len(path) > 1:
                        angle_candidates.append(
                            (record["amide_n"], path[0], path[1]))
                else:
                    path_axis_idx = axis_index - 1
                    left_angle_atom = (
                        record["amide_n"]
                        if path_axis_idx == 0 else path[path_axis_idx - 1]
                    )
                    right_angle_atom = (
                        path[path_axis_idx + 2]
                        if path_axis_idx + 2 < len(path)
                        else record["ring_anchor"]
                    )
                    angle_candidates = [
                        (left_angle_atom, path[path_axis_idx],
                         path[path_axis_idx + 1]),
                        (path[path_axis_idx], path[path_axis_idx + 1],
                         right_angle_atom),
                    ]

                for angle in angle_candidates:
                    add_signature_rows(
                        signature_rows,
                        signature_row_types,
                        signature_row_scales,
                        angle_rows_for_angle(angle),
                        "angle",
                        angle_scale,
                    )

                if axis_index == len(candidate_axes) - 1:
                    boundary_axis = (path[-1], record["ring_anchor"])
                    add_signature_rows(
                        signature_rows,
                        signature_row_types,
                        signature_row_scales,
                        torsion_rows_for_axis(boundary_axis),
                        "torsion",
                        1.0,
                    )

                rotor_id = (
                    f"{cluster_id}_axis_{axis[0]}_{axis[1]}"
                )
                rotor = LocalRotor(
                    rotor_id=rotor_id,
                    kind="anchored_linker",
                    axis=axis,
                    symmetry_order=1,
                    owned_atoms=owned_atoms,
                    unit_atom_sets=(),
                    torsion_rows=axis_torsion_rows,
                    phase_coordinate=phase_coordinate,
                    phase_values=phase_values,
                    signature_rows=tuple(int(row) for row in signature_rows),
                    signature_row_types=tuple(signature_row_types),
                    signature_row_scales=tuple(signature_row_scales),
                )

                local_group_model.rotors[rotor_id] = rotor
                registered_axes.add(axis_key)
                rotor_ids.append(rotor_id)
                active_signature_rows.update(signature_rows)

            if not rotor_ids:
                continue

            active_rows = set(rows_for_local_atoms(owned_atoms, active_seed))
            active_rows.update(active_signature_rows)

            active_atoms = set(active_seed)
            for row_idx in active_rows:
                active_atoms.update(
                    int(atom) for atom in flat_z_matrix[row_idx]
                )

            local_group_model.clusters[cluster_id] = LocalCluster(
                cluster_id=cluster_id,
                family_label="anchored_linker",
                rotor_ids=tuple(rotor_ids),
                owned_atoms=owned_atoms,
                active_atoms=tuple(sorted(active_atoms)),
                active_rows=tuple(sorted(active_rows)),
                anchor_atoms=anchor_atoms,
                relaxation_atoms=relaxation_atoms,
            )
            claimed_owned_atoms.update(owned_set)

        return local_group_model

    def _detect_alkyl_chain_local_groups(
            self,
            molecule,
            ff_gen,
            z_matrix,
            local_group_model,
            rotatable_bonds_zero_based=None):

        enabled_kinds = set(getattr(self, "auto_local_group_kinds", ()))
        if not ({"alkyl", "alkyl_chain"} & enabled_kinds):
            return local_group_model

        labels = tuple(str(x).strip() for x in molecule.get_labels())
        natoms = len(labels)
        conn = np.asarray(ff_gen.connectivity_matrix, dtype=int)
        atom_types = tuple(str(x).strip().lower() for x in ff_gen.atom_types)
        flat_z_matrix = InterpolationDatapoint.flatten_z_matrix(z_matrix)
        phase_values = tuple(
            np.deg2rad(float(x))
            for x in getattr(
                self,
                "alkyl_chain_phase_values_degrees",
                (0, 60, 120, 180, 240, 300),
            )
        )

        def element(atom):
            letters = "".join(x for x in labels[int(atom)] if x.isalpha())
            return letters.capitalize()

        def neighbors(atom):
            return tuple(int(x) for x in np.flatnonzero(conn[int(atom)]))

        def canonical_axis(axis):
            return tuple(sorted(int(x) for x in axis))

        def is_alkyl_carbon(atom):
            return (
                element(atom) == "C"
                and len(atom_types) == natoms
                and atom_types[int(atom)] == "c3"
            )

        def component_on_side(start, blocked_axis):
            blocked = frozenset(int(x) for x in blocked_axis)
            visited = set()
            stack = [int(start)]

            while stack:
                atom = stack.pop()
                if atom in visited:
                    continue
                visited.add(atom)

                for neighbor in neighbors(atom):
                    if frozenset((atom, neighbor)) == blocked:
                        continue
                    if neighbor not in visited:
                        stack.append(neighbor)

            return tuple(sorted(visited))

        def torsion_rows_for_axis(axis):
            axis_set = set(int(x) for x in axis)
            return tuple(
                row_idx
                for row_idx, coordinate in enumerate(flat_z_matrix)
                if (
                    len(coordinate) == 4
                    and set(int(x) for x in coordinate[1:3]) == axis_set
                )
            )

        def oriented_coordinate(coordinate, axis):
            coordinate = tuple(int(x) for x in coordinate)
            axis = tuple(int(x) for x in axis)

            if tuple(coordinate[1:3]) == axis:
                return coordinate
            if tuple(coordinate[2:0:-1]) == axis:
                return tuple(reversed(coordinate))
            return None

        def select_phase_coordinate(axis, torsion_rows, moving_atoms):
            moving_atoms = set(int(x) for x in moving_atoms)
            fallback = None

            for row in torsion_rows:
                coordinate = oriented_coordinate(flat_z_matrix[row], axis)
                if coordinate is None:
                    continue
                if fallback is None:
                    fallback = coordinate
                if coordinate[3] in moving_atoms:
                    return coordinate

            return fallback

        def rows_for_local_atoms(local_atoms, active_atoms):
            local_set = set(int(x) for x in local_atoms)
            active_set = set(int(x) for x in active_atoms)
            rows = []
            for row_idx, coordinate in enumerate(flat_z_matrix):
                coord_atoms = set(int(x) for x in coordinate)
                if coord_atoms <= active_set and coord_atoms & local_set:
                    rows.append(int(row_idx))
            return tuple(rows)

        def side_is_alkyl_substituent(side_atoms):
            heavy_atoms = tuple(
                atom for atom in side_atoms
                if element(atom) != "H"
            )
            carbon_atoms = tuple(
                atom for atom in heavy_atoms
                if element(atom) == "C"
            )

            if len(carbon_atoms) < 2:
                return False
            if len(carbon_atoms) != len(heavy_atoms):
                return False
            return all(is_alkyl_carbon(atom) for atom in carbon_atoms)

        def complement_has_nonalkyl_core(side_atoms):
            side_set = set(int(x) for x in side_atoms)
            for atom in range(natoms):
                if atom in side_set or element(atom) == "H":
                    continue
                if not is_alkyl_carbon(atom):
                    return True
            return False

        def side_distances(start, side_atoms):
            side_set = set(int(x) for x in side_atoms)
            distances = {int(start): 0}
            stack = [int(start)]

            while stack:
                atom = stack.pop(0)
                for neighbor in neighbors(atom):
                    if neighbor not in side_set:
                        continue
                    if neighbor in distances:
                        continue
                    distances[neighbor] = distances[atom] + 1
                    stack.append(neighbor)

            return distances

        if rotatable_bonds_zero_based is None:
            rotatable_bonds_zero_based = tuple(
                canonical_axis((int(a) - 1, int(b) - 1))
                for a, b in getattr(ff_gen, "rotatable_bonds", ())
            )

        rotatable_axes = {
            canonical_axis(axis) for axis in rotatable_bonds_zero_based
        }
        axis_to_rotor_ids = {}
        for rotor_id, rotor in local_group_model.rotors.items():
            axis_to_rotor_ids.setdefault(
                canonical_axis(rotor.axis), []).append(str(rotor_id))

        tertbutyl_parent_axes = set()
        tertbutyl_parent_owned_sets = []
        for cluster in local_group_model.clusters.values():
            if str(getattr(cluster, "family_label", "")) != "tertbutyl_parent":
                continue

            tertbutyl_parent_owned_sets.append(
                set(int(atom) for atom in cluster.owned_atoms)
            )
            for rotor_id in cluster.rotor_ids:
                rotor = local_group_model.rotors.get(str(rotor_id))
                if rotor is None:
                    continue
                if str(getattr(rotor, "kind", "")) != "tertbutyl_parent":
                    continue
                tertbutyl_parent_axes.add(canonical_axis(rotor.axis))

        candidates = []
        for axis in sorted(rotatable_axes):
            atom_a, atom_b = axis
            if element(atom_a) != "C" or element(atom_b) != "C":
                continue

            for anchor, first_atom in ((atom_a, atom_b), (atom_b, atom_a)):
                if canonical_axis((anchor, first_atom)) in tertbutyl_parent_axes:
                    continue

                side_atoms = component_on_side(first_atom, axis)
                if anchor in side_atoms:
                    continue
                side_set = set(int(atom) for atom in side_atoms)
                if any(
                    side_set <= tertbutyl_atoms
                    for tertbutyl_atoms in tertbutyl_parent_owned_sets
                ):
                    continue
                if not side_is_alkyl_substituent(side_atoms):
                    continue
                if not complement_has_nonalkyl_core(side_atoms):
                    continue

                candidates.append({
                    "axis": (int(anchor), int(first_atom)),
                    "anchor": int(anchor),
                    "first_atom": int(first_atom),
                    "owned_atoms": tuple(sorted(side_atoms)),
                })

        candidates = sorted(
            candidates,
            key=lambda item: (-len(item["owned_atoms"]), item["axis"]),
        )
        maximal_candidates = []
        for candidate in candidates:
            owned = set(candidate["owned_atoms"])
            if any(
                owned < set(other["owned_atoms"])
                for other in maximal_candidates
            ):
                continue
            maximal_candidates.append(candidate)

        for candidate in maximal_candidates:
            cluster_id = f"alkyl_chain_{candidate['first_atom']}"
            if cluster_id in local_group_model.clusters:
                continue

            owned_atoms = tuple(int(x) for x in candidate["owned_atoms"])
            owned_set = set(owned_atoms)
            boundary_axis = canonical_axis(candidate["axis"])
            distances = side_distances(candidate["first_atom"], owned_atoms)
            rotor_ids = []
            torsion_rows = set()

            for axis_key in sorted(rotatable_axes):
                axis_set = set(axis_key)
                if axis_key != boundary_axis and not axis_set <= owned_set:
                    continue
                if not all(element(atom) == "C" for atom in axis_key):
                    continue
                if axis_key != boundary_axis and not all(
                        is_alkyl_carbon(atom) for atom in axis_key):
                    continue

                if axis_key == boundary_axis:
                    oriented_axis = tuple(candidate["axis"])
                else:
                    atom_a, atom_b = axis_key
                    dist_a = distances.get(atom_a)
                    dist_b = distances.get(atom_b)
                    if dist_a is None or dist_b is None or dist_a == dist_b:
                        continue
                    oriented_axis = (
                        (atom_a, atom_b)
                        if dist_a < dist_b else (atom_b, atom_a)
                    )

                existing_rotors = axis_to_rotor_ids.get(
                    canonical_axis(oriented_axis), ())
                if existing_rotors:
                    for rotor_id in existing_rotors:
                        if rotor_id not in rotor_ids:
                            rotor_ids.append(rotor_id)
                        torsion_rows.update(
                            int(row)
                            for row in local_group_model.rotors[
                                rotor_id].torsion_rows
                        )
                    continue

                moving_atoms = component_on_side(
                    oriented_axis[1], oriented_axis)
                if not set(moving_atoms) <= owned_set:
                    continue

                axis_torsion_rows = torsion_rows_for_axis(oriented_axis)
                if not axis_torsion_rows:
                    continue

                phase_coordinate = select_phase_coordinate(
                    oriented_axis, axis_torsion_rows, moving_atoms)
                if phase_coordinate is None:
                    continue

                rotor_id = (
                    f"{cluster_id}_axis_{oriented_axis[0]}_"
                    f"{oriented_axis[1]}"
                )
                rotor = LocalRotor(
                    rotor_id=rotor_id,
                    kind="alkyl_chain",
                    axis=tuple(int(x) for x in oriented_axis),
                    symmetry_order=1,
                    owned_atoms=tuple(sorted(int(x) for x in moving_atoms)),
                    unit_atom_sets=(),
                    torsion_rows=axis_torsion_rows,
                    phase_coordinate=phase_coordinate,
                    phase_values=phase_values,
                )

                local_group_model.rotors[rotor_id] = rotor
                axis_to_rotor_ids.setdefault(
                    canonical_axis(oriented_axis), []).append(rotor_id)
                rotor_ids.append(rotor_id)
                torsion_rows.update(int(row) for row in axis_torsion_rows)

            if not rotor_ids:
                continue

            active_seed = set(owned_atoms) | {int(candidate["anchor"])}
            active_rows = set(rows_for_local_atoms(owned_atoms, active_seed))
            active_rows.update(torsion_rows)

            active_atoms = set(active_seed)
            for row_idx in active_rows:
                active_atoms.update(
                    int(atom) for atom in flat_z_matrix[row_idx]
                )

            local_group_model.clusters[cluster_id] = LocalCluster(
                cluster_id=cluster_id,
                family_label="alkyl_chain",
                rotor_ids=tuple(rotor_ids),
                owned_atoms=owned_atoms,
                active_atoms=tuple(sorted(active_atoms)),
                active_rows=tuple(sorted(active_rows)),
            )

        return local_group_model

    def _build_manual_phenyl_local_group(self, z_matrix, spec):
        """Builds a non-symmetric phenyl local group from a user specification.

        Manual local-group atom indices use VeloxChem's one-based convention.
        """

        def _torsion_rows_for_axis(z_matrix, axis):
            axis_set = set(int(x) for x in axis)
            rows = []
            for row_idx, coord in enumerate(z_matrix):
                if len(coord) != 4:
                    continue
                if set(int(x) for x in coord[1:3]) == axis_set:
                    rows.append(int(row_idx))
            return tuple(rows)

        def _rows_for_local_atoms(z_matrix, local_atoms):
            local_set = set(int(x) for x in local_atoms)
            rows = []
            for row_idx, coord in enumerate(z_matrix):
                coord_atoms = set(int(x) for x in coord)
                if coord_atoms & local_set:
                    rows.append(int(row_idx))
            return tuple(rows)

        def _normalize_indices(values, field_name, expected_size=None, unique=False):
            if values is None:
                raise ValueError(
                    f"Manual phenyl group requires '{field_name}'."
                )

            normalized = tuple(int(value) - 1 for value in values)
            if expected_size is not None and len(normalized) != expected_size:
                raise ValueError(
                    f"Manual phenyl '{field_name}' must contain "
                    f"{expected_size} atom indices."
                )
            if any(atom < 0 or atom >= natoms for atom in normalized):
                raise ValueError(
                    f"Manual phenyl '{field_name}' contains an atom outside "
                    f"the valid one-based range 1..{natoms}."
                )
            if unique:
                normalized = tuple(sorted(set(normalized)))
            return normalized

        allowed_manual_rotors = {
            "phenyl": 1,
            "alcohol": 1,
            "amine": 1,
            "primary_amine": 2,
            "ammonium": 3,
        }

        kind = str(spec["kind"]).lower()
        symmetry_order = int(spec.get("symmetry_order", allowed_manual_rotors[kind]))

        flat_z_matrix = InterpolationDatapoint.flatten_z_matrix(z_matrix)
        natoms = 1 + max(
            int(atom)
            for coord in flat_z_matrix
            for atom in coord
        )

        declared_index_base = int(spec.get("index_base", 1))
        if declared_index_base != 1:
            raise ValueError(
                "Manual local-group atom indices must be one-based."
            )

        owned_atoms = _normalize_indices(
            spec.get("owned_atoms"), "owned_atoms", unique=True)
        axis = _normalize_indices(spec.get("axis"), "axis", expected_size=2)
        phase_coordinate = _normalize_indices(
            spec.get("phase_coordinate"),
            "phase_coordinate",
            expected_size=4,
        )
        phase_values = tuple(np.deg2rad(x) for x in spec["phase_values_degrees"])

        if not owned_atoms:
            raise ValueError("Manual phenyl owned_atoms cannot be empty.")
        if set(phase_coordinate[1:3]) != set(axis):
            raise ValueError(
                "Manual phenyl phase_coordinate must use the specified axis "
                "as its central bond."
            )
        if len(set(axis) & set(owned_atoms)) != 1:
            raise ValueError(
                "Manual phenyl axis must contain exactly one owned phenyl atom "
                "and one core-side boundary atom."
            )

        torsion_rows = _torsion_rows_for_axis(flat_z_matrix, axis)
        if not torsion_rows:
            raise ValueError(
                f"No z-matrix torsion rows were found for phenyl axis {axis}."
            )

        rotor = LocalRotor(
            rotor_id=spec["group_id"],
            kind="phenyl",
            axis=axis,
            symmetry_order=1,
            owned_atoms=owned_atoms,
            unit_atom_sets=(),
            torsion_rows=torsion_rows,
            phase_coordinate=phase_coordinate,
            phase_values=phase_values,
        )

        active_rows = tuple(sorted(
            set(_rows_for_local_atoms(flat_z_matrix, owned_atoms))
            | set(torsion_rows)
        ))
        active_atom_set = set(owned_atoms) | set(axis)
        for row_idx in active_rows:
            active_atom_set.update(int(atom) for atom in flat_z_matrix[row_idx])
        active_atoms = tuple(sorted(active_atom_set))

        cluster = LocalCluster(
            cluster_id=spec["group_id"],
            family_label="phenyl",
            rotor_ids=(rotor.rotor_id,),
            owned_atoms=owned_atoms,
            active_atoms=active_atoms,
            active_rows=active_rows,
        )

        return rotor, cluster

    def _refresh_local_group_core_partition(
            self,
            local_group_model,
            z_matrix,
            natoms):
        """Recomputes the core after adding manual local groups."""
        flat_z_matrix = InterpolationDatapoint.flatten_z_matrix(z_matrix)
        local_owned_atoms = {
            int(atom)
            for cluster in local_group_model.clusters.values()
            for atom in cluster.owned_atoms
        }
        local_active_rows = {
            int(row)
            for cluster in local_group_model.clusters.values()
            for row in cluster.active_rows
        }

        local_group_model.core_atoms = tuple(
            sorted(set(range(int(natoms))) - local_owned_atoms)
        )
        local_group_model.core_rows = tuple(
            row_idx
            for row_idx in range(len(flat_z_matrix))
            if row_idx not in local_active_rows
        )
        local_group_model.enabled = bool(local_group_model.clusters)

    def set_up_the_system(self, molecule, imforcefieldfiles=None, exclude_rot_bonds=None):

        """
        Assign the neccessary variables with respected values.

        :param molecule: original molecule

        :param target_dihedrals: is a list of dihedrals that should be scanned during the dynamics

        :param sampling_structures: devides the searchspace around given rotatbale dihedrals

        """
        
        def filter_rotatable_bonds(rotatable_bonds, exclude_rot_bonds=[]):
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
                    'IMForceFieldGenerator.filter_user_excluded_rotatable_bonds: '
                    'excluded_rotatable_bonds must contain atom-index pairs.')

                atom_i, atom_j = int(bond[0]), int(bond[1])

                assert_msg_critical(
                    atom_i > 0 and atom_j > 0 and atom_i != atom_j,
                    'IMForceFieldGenerator.filter_user_excluded_rotatable_bonds: '
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
            Promote proper ring torsions around non-rotatable bonds to impropers.

            Matching by full 4-atom set is too strict (and often fails), so we match
            by the middle bond of the proper torsion and collect impropers touching
            the same bond.
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
        # 7. list of indices grouped for all dihedrals attached to a rotabale bond (z_matrix index)
        # 8. actual dihedrals with the indices assigned to each atom
        # 9. start and end of the dihedral section within the z_matrix
        self.symmetry_information = {'gs': (), 'es': ()}

        self.molecule = molecule

        root_extract_z_matrix = None
        if imforcefieldfiles is not None:
            root_extract_z_matrix = {}
            self.sampling_imforcefieldfiles = {}
            self.imforcefieldfiles = {}
            for i, root in enumerate(self.roots_to_follow):
                if root not in imforcefieldfiles:
                    self.imforcefieldfiles[self.roots_to_follow[i]] = f'im_database_{root}.h5'
                    root_extract_z_matrix[self.roots_to_follow[i]] = False
                else:
                    root_extract_z_matrix[self.roots_to_follow[i]] = True
                    self.imforcefieldfiles[self.roots_to_follow[i]] = imforcefieldfiles[root]

                self.sampling_imforcefieldfiles[self.roots_to_follow[i]] = f'im_database_sampling_{root}.h5'
        else:
            self.imforcefieldfiles = {}
            self.sampling_imforcefieldfiles = {}
            standard_files: list[str] = [f'im_database_{root}.h5' for root in self.roots_to_follow]
            standard_smapling_files = [f'im_database_sampling_{root}.h5' for root in self.roots_to_follow]
            root_extract_z_matrix = {}
            for root_idx, standard_file in enumerate(standard_files):

                if Path(standard_file).exists():
                    if self.roots_to_follow[root_idx] not in root_extract_z_matrix:
                        root_extract_z_matrix[self.roots_to_follow[root_idx]] = True
                else:
                    root_extract_z_matrix[self.roots_to_follow[root_idx]] = False
                self.imforcefieldfiles[self.roots_to_follow[root_idx]] = standard_file
                self.sampling_imforcefieldfiles[self.roots_to_follow[root_idx]] = standard_smapling_files[root_idx]

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
                    'eq_bond_symmetry_mode':self.eq_bond_symmetry_mode,
                    'use_tc_weights':self.use_tc_weights,
                    'tc_weight_mode':self.tc_weight_mode,
                    'local_factor_state_weight_mode':
                        self.local_factor_state_weight_mode,
                    'local_factor_combination_mode':
                        self.local_factor_combination_mode,
                    'use_mass_weight':self.use_mass_weight,
                })

                _, z_matrix = int_driver.read_labels()
                self.roots_z_matrix[root] = z_matrix

            dihedral_start = len(self.roots_z_matrix[root]['bonds']) + len(self.roots_z_matrix[root]['angles'])
            dihedral_end = dihedral_start + len(self.roots_z_matrix[root]['dihedrals'])

            self.all_rotatable_bonds = rotatable_bonds_zero_based
            excluded_local_atoms = set() 
            if self.use_local_group_database:
                self.local_group_primitive_model[root] = self._detect_methyl_and_tertbutyl_local_groups(
                    molecule=molecule,
                    ff_gen=ff_gen,
                    z_matrix=self.roots_z_matrix[root],
                    rotatable_bonds_zero_based=rotatable_bonds_zero_based,
                )

                self.local_group_primitive_model[root] = (
                    self._detect_heteroatom_local_groups(
                        molecule=molecule,
                        ff_gen=ff_gen,
                        z_matrix=self.roots_z_matrix[root],
                        local_group_model=self.local_group_primitive_model[root],
                        rotatable_bonds_zero_based=rotatable_bonds_zero_based,
                    )
                )

                self.local_group_primitive_model[root] = (
                    self._detect_methoxy_local_groups(
                        molecule=molecule,
                        ff_gen=ff_gen,
                        z_matrix=self.roots_z_matrix[root],
                        local_group_model=self.local_group_primitive_model[root],
                        rotatable_bonds_zero_based=rotatable_bonds_zero_based,
                    )
                )

                self.local_group_primitive_model[root] = (
                    self._detect_amide_side_chain_local_groups(
                        molecule=molecule,
                        ff_gen=ff_gen,
                        z_matrix=self.roots_z_matrix[root],
                        local_group_model=self.local_group_primitive_model[root],
                        rotatable_bonds_zero_based=rotatable_bonds_zero_based,
                    )
                )

                self.local_group_primitive_model[root] = (
                    self._detect_anchored_linker_local_groups(
                        molecule=molecule,
                        ff_gen=ff_gen,
                        z_matrix=self.roots_z_matrix[root],
                        local_group_model=self.local_group_primitive_model[root],
                        rotatable_bonds_zero_based=rotatable_bonds_zero_based,
                    )
                )

                self.local_group_primitive_model[root] = (
                    self._detect_alkyl_chain_local_groups(
                        molecule=molecule,
                        ff_gen=ff_gen,
                        z_matrix=self.roots_z_matrix[root],
                        local_group_model=self.local_group_primitive_model[root],
                        rotatable_bonds_zero_based=rotatable_bonds_zero_based,
                    )
                )
    
                print("\n 2 \n", self.local_group_primitive_model)

                manual_specs = self.local_group_specs or ()
                if isinstance(manual_specs, dict):
                    manual_specs = (manual_specs,)

                for spec in manual_specs:
                    kind = str(spec.get("kind", "")).lower()
                    if kind != "phenyl":
                        raise ValueError(
                            f"Unsupported manual local-group kind '{kind}'."
                        )

                    rotor, cluster = self._build_manual_phenyl_local_group(
                        self.roots_z_matrix[root], spec)
                    group_id = str(spec["group_id"])
                    if (
                        group_id in self.local_group_primitive_model[root].rotors
                        or group_id in self.local_group_primitive_model[root].clusters
                    ):
                        raise ValueError(
                            f"Duplicate local-group id '{group_id}'."
                        )
                    self.local_group_primitive_model[root].rotors[group_id] = rotor
                    self.local_group_primitive_model[root].clusters[group_id] = cluster

                self._refresh_local_group_core_partition(
                    self.local_group_primitive_model[root],
                    self.roots_z_matrix[root],
                    len(molecule.get_labels()),
                )

                self.local_group_model[root] = self.local_group_primitive_model[root]
                self.local_groups[root] = self.local_group_model[root].clusters

                excluded_local_atoms = set()
                for cluster in self.local_group_model[root].clusters.values():
                    cluster_owned = set(int(atom) for atom in cluster.owned_atoms)

                    axis_atoms = set()
                    for rotor_id in cluster.rotor_ids:
                        rotor = self.local_group_model[root].rotors[str(rotor_id)]
                        axis_atoms.update(int(atom) for atom in rotor.axis)

                    excluded_local_atoms.update(cluster_owned - axis_atoms)

                excluded_local_atoms = sorted(excluded_local_atoms)
                      
            all_exclision = [element for rot_bond in rotatable_bonds_zero_based for element in rot_bond]

            symmetry_groups_ref = [groups for groups in symmetry_groups[1] if not any(item in all_exclision for item in groups)]

            # reduce the symmetry to only CH3 or CH2 symmetry groups for the time being
            regrouped, rot_groups = regroup_by_rotatable_connection(molecule, symmetry_groups_ref, rotatable_bonds_zero_based, molecule.get_connectivity_matrix())

            if root == 0:
                non_core_atoms = [element for group in regrouped['gs'] for element in group]
                core_atoms = [element for element in symmetry_groups[0] if element not in excluded_local_atoms]
                # angles_to_set, _, _, self.symmetry_dihedral_lists, dih_list = self._adjust_symmetry_dihedrals(molecule, rot_groups['gs'], rotatable_bonds_zero_based, self.roots_z_matrix[root])
                # dihedrals_to_set = {key: [] for key in angles_to_set.keys()}
                indices_list = []
                # for key, dihedral_list in self.symmetry_dihedral_lists.items():

                #     for i, element in enumerate(self.roots_z_matrix[root]['dihedrals']):
                #         if tuple(sorted(element)) in dihedral_list:
                #             indices_list.append(i)
                self.symmetry_information['gs'] = [symmetry_groups[0], rot_groups['gs'], regrouped['gs'], core_atoms, list(excluded_local_atoms), rotatable_bonds_zero_based, indices_list, self.symmetry_dihedral_lists, [], [dihedral_start, dihedral_end]]

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
                'eq_bond_symmetry_mode': self.eq_bond_symmetry_mode,
                'use_tc_weights': self.use_tc_weights,
                'tc_weight_mode': self.tc_weight_mode,
                'local_factor_state_weight_mode':
                    self.local_factor_state_weight_mode,
                'local_factor_combination_mode':
                    self.local_factor_combination_mode,
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

        if self.add_conformal_structures and 0 in self.roots_to_follow:

            conformers_plus_ts = {0: {}}

            conformer_gen = ConformerGenerator()
            conformer_gen.resp_charges = False

            rot_bonds_to_exclude = []
            if self.use_local_group_database:
                grouped_bonds = set()
                for model in self.local_group_primitive_model.values():
                    if not getattr(model, 'enabled', False):
                        continue
                    grouped_bonds.update(
                        tuple(sorted(int(atom) + 1 for atom in rotor.axis))
                        for rotor in model.rotors.values()
                    )
                rot_bonds_to_exclude = sorted(grouped_bonds)

            conformer_results = conformer_gen.generate(molecule)
            dihedral_candidates = list(getattr(conformer_gen, "dihedral_candidates", []) or [])
            conformal_molecules = list(conformer_results.get("molecules", []) or [])

            if len(conformal_molecules) == 0:
                raise RuntimeError('ConformerGenerator returned no conformers.')

            rot_bond_set = {
                tuple(sorted((int(a), int(b))))
                for a, b in rotatable_bonds
            }
            rot_bond_set.difference_update(rot_bonds_to_exclude)

            selected_dihedrals = []
            for dih_zero, _angle_grid in dihedral_candidates:
                dih_key = tuple(int(x) + 1 for x in dih_zero)

                if tuple(sorted((dih_key[1], dih_key[2]))) not in rot_bond_set:
                    continue

                selected_dihedrals.append((dih_key, len(_angle_grid)))

            if len(selected_dihedrals) > 0:
                counter = 0
                for dih_key, nconf in selected_dihedrals:
                    conformers_plus_ts[0][dih_key] = []
                    
                    for conf_mol in conformal_molecules[counter:counter + nconf]:
                        mol_i = Molecule.from_xyz_string(conf_mol.get_xyz_string())
                        mol_i.set_charge(molecule.get_charge())
                        mol_i.set_multiplicity(molecule.get_multiplicity())
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
                    counter+=nconf
            else:
                mol_i = Molecule.from_xyz_string(conformal_molecules[0].get_xyz_string())
                mol_i.set_charge(molecule.get_charge())
                mol_i.set_multiplicity(molecule.get_multiplicity())
                conformers_plus_ts[0][None] = [(mol_i, 'normal')]

            self.seed_structures = conformers_plus_ts

        self._system_is_set_up = True

    def _bootstrap_sampling_db_from_abinito_db(self, root):
        """
        creating a sample database of a cheap tight-binding method
        for an efficient sampling mode allowing to estimate when ab-inito
        reference calculation is required (only for ground-state for now!)
        """

        # compare to the current interpolation database in order to get a
        # identical mirroring of the ab-inito database -- note that this is
        # currently only for the ground-state
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
        sampling_qm, sampling_grad, sampling_hess = self.sampling_driver['gs']  # if rank == root_rank else None, None, None

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
            backend_metadata = self._collect_backend_metadata(
                self.sampling_driver['gs'])

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
            samp_dp.metadata = deepcopy(backend_metadata)
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

    def compute(self, molecule, states_basis=None):

        """
        Construct the interpolation dynamics database by generating molecular structures,
        performing QM calculations, and collecting data points.

        :param molecule: The input molecular structure for which the database is constructed.
                        This molecule serves as a reference for sampling and simulation tasks.

        The method sets up quantum mechanical drivers, initializes interpolation and dynamics settings,
        and iterates over molecular structures along a predefined reaction path. It runs simulations
        to expand/generate the interpolation forcefield with new data points.

        """

        if self.reference_struc_energy_file is not None and not Path(self.reference_struc_energy_file).exists():
            self.reference_struc_energy_file = None

        # First set up the system for which the database needs to be constructed
        states_basis = {'gs':self.gs_basis_set_label, 'es':self.es_basis_set_label}
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
                'pressure': self.force_constant, 'force_constant': self.force_constant, 'ensemble': self.ensemble,
                'timestep': self.timestep, 'nsteps': self.nsteps, 'friction': self.friction,
                'snapshots': self.snapshots, 'trajectory_file': self.trajectory_file, 'reference_struc_energy_file': self.reference_struc_energy_file,
                'desired_datapoint_density': self.desired_point_density, 'converged_cycle': self.converged_cycle,
                'energy_threshold': self.energy_threshold, 'grad_rmsd_thrsh': self.gradient_rmsd_thrsh,
                'load_system': None, 'start_collect': self.start_collect, 'roots_to_follow': self.roots_to_follow,
                'sampling_drivers': self.sampling_driver, 'sampling_settings': self.sampling_settings,
                'metadynamics': self.metadynamics_settings,
            }

            files_to_add_conf = []
            molecules_to_add_info = []

            if not Path(self.imforcefieldfiles[self.roots_to_follow[0]]).exists():
                files_to_add_conf.append(self.roots_to_follow[0])

            if self.seed_structures is not None and not Path(imforcefieldfile).exists():

                molecules_to_add_info = []
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

                                    if mode == "constraint" and dih_key is not None:
                                        if dih_key not in constraints_global:
                                            constraints_global.append(dih_key)

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

                                elif self.roots_to_follow[0] == 0 and self._is_xtb_like_driver(self.drivers['gs'][0]):

                                    opt_results = self._run_optimization(
                                        self.drivers['gs'][0],
                                        mol,
                                        constraints=constraints_global,
                                        transition=transition,
                                        index_offset=0,
                                    )
                                    optimized_molecule = opt_results['final_molecule']
                                    optimized_molecule_for_scan = optimized_molecule
                                    current_basis = None
                                    current_molecule_to_add_info.append((optimized_molecule, current_basis, self.roots_to_follow, constraints_global,transition))
                                    self.ostream.print_blank()
                                    self.ostream.print_header('Optimized Molecule')
                                    self.ostream.print_header('------------------')
                                    self.ostream.print_block(optimized_molecule.get_xyz_string())
                                    self.ostream.print_blank()
                                    self.ostream.flush()

                            else:
                                if self.roots_to_follow[0] == 0:
                                    current_basis = self._read_basis_for_driver(
                                        mol, self.drivers['gs'][0], states_basis['gs'])
                                    current_molecule_to_add_info.append((mol, current_basis, self.roots_to_follow, []))
                                else:
                                    current_basis = self._read_basis_for_driver(
                                        mol, self.drivers['es'][0], states_basis['es'])
                                    current_molecule_to_add_info.append((mol, current_basis, self.roots_to_follow, []))

                            if mode == "normal" and optimized_molecule_for_scan is not None:
                                optimized_normal_conformers.append(optimized_molecule_for_scan)

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

                                def _circ_err_deg(actual, target):
                                    return abs(((actual - target + 180.0) % 360.0) - 180.0)

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
                                        mol_i = Molecule.from_xyz_string(base_mol.get_xyz_string())
                                        mol_i.set_charge(molecule.get_charge())
                                        mol_i.set_multiplicity(molecule.get_multiplicity())
                                        mol_i.set_dihedral_in_degrees(dih_key, target_deg, verbose=False)

                                        actual = float(mol_i.get_dihedral_in_degrees(dih_key)) % 360.0
                                        if _circ_err_deg(actual, target_deg) > 2.0:
                                            raise RuntimeError(
                                                f"Failed to set dihedral {dih_key} to {target_deg:.2f} deg "
                                                f"(actual {actual:.2f} deg)."
                                            )

                                        mol_info.append((mol_i, scan_mode))

                        if len(molecules_to_add_info) == 0:

                            molecules_to_add_info.append(current_molecule_to_add_info[0])

                        self.add_point(current_molecule_to_add_info, self.states_interpolation_settings, symmetry_information=self.symmetry_information)

            elif self.seed_structures is None and not Path(imforcefieldfile).exists():

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

                    elif self.roots_to_follow[0] == 0 and self._is_xtb_like_driver(self.drivers['gs'][0]):

                        opt_results = self._run_optimization(
                            self.drivers['gs'][0],
                            molecule,
                            constraints=self.use_minimized_structures[1],
                            index_offset=0,
                        )

                        optimized_molecule = opt_results['final_molecule']

                        current_basis = None
                        molecules_to_add_info.append((optimized_molecule, current_basis, self.roots_to_follow, self.use_minimized_structures[1]))

                        self.ostream.print_blank()
                        self.ostream.print_header('Optimized Molecule')
                        self.ostream.print_header('------------------')
                        self.ostream.print_block(optimized_molecule.get_xyz_string())
                        self.ostream.print_blank()
                        self.ostream.flush()

                    self.add_point(molecules_to_add_info, self.states_interpolation_settings, symmetry_information=self.symmetry_information)

                else:
                    if self.roots_to_follow[0] == 0:
                        current_basis = self._read_basis_for_driver(
                            molecule, self.drivers['gs'][0], states_basis['gs'])
                        molecules_to_add_info.append((molecule, current_basis, self.roots_to_follow, []))
                    else:
                        current_basis = self._read_basis_for_driver(
                            molecule, self.drivers['es'][0], states_basis['es'])
                        molecules_to_add_info.append((molecule, current_basis, self.roots_to_follow, []))

                    self.add_point(molecules_to_add_info, self.states_interpolation_settings, symmetry_information=self.symmetry_information)

            else:

                if self.roots_to_follow[0] == 0:
                    current_basis = self._read_basis_for_driver(
                        molecule, self.drivers['gs'][0], states_basis['gs'])
                    molecules_to_add_info.append((molecule, current_basis, self.roots_to_follow, []))
                else:
                    current_basis = self._read_basis_for_driver(
                        molecule, self.drivers['es'][0], states_basis['es'])
                    molecules_to_add_info.append((molecule, current_basis, self.roots_to_follow, []))

                if not Path(imforcefieldfile).exists():
                    self.add_point(molecules_to_add_info, self.states_interpolation_settings, symmetry_information=self.symmetry_information)

            density_of_datapoints = self.determine_datapoint_density(self.states_interpolation_settings)

            self.states_data_point_density = density_of_datapoints

            if self.sampling_settings.get('enabled', False):
                self._bootstrap_sampling_db_from_abinito_db(self.roots_to_follow[0])

            dynamics_molecule = molecules_to_add_info[self.roots_to_follow[0]][0]
            forcefield_generator = MMForceFieldGenerator(ostream=self.ostream)
            self.dynamics_settings['trajectory_file'] = f'trajectory_{self.roots_to_follow[0]}.pdb'
            forcefield_generator.partial_charges = dynamics_molecule.get_partial_charges(dynamics_molecule.get_charge())

            forcefield_generator.create_topology(dynamics_molecule)
            im_database_driver = IMDatabasePointCollecter(ostream=self.ostream)
            im_database_driver.distance_thrsh = self.distance_thrsh
            im_database_driver.non_core_symmetry_groups = self.symmetry_information
            im_database_driver.platform = self.open_mm_platform
            im_database_driver.all_rot_bonds = self.all_rotatable_bonds
            im_database_driver.consider_locality = self.consider_locality

            # set optimization features in the construction run
            im_database_driver.identfy_relevant_int_coordinates = (self.identfy_relevant_int_coordinates, self.use_minimized_structures[1])
            im_database_driver.use_opt_confidence_radius = self.use_opt_confidence_radius

            if (
                self.use_local_group_database
                and any(
                    getattr(model, 'enabled', False)
                    for model in self.local_group_primitive_model.values()
                )
            ):
                im_database_driver.local_group_family_builder = (
                    self.add_point_local_groups
                )
                im_database_driver.local_group_primitive_models = (
                    self.local_group_primitive_model
                )
                im_database_driver.relax_local_group_rotors_during_optimization = (
                    self.relax_local_group_rotors_during_optimization
                )
                im_database_driver.relax_alkyl_local_group_states = (
                    self.relax_alkyl_local_group_states
                )
                im_database_driver.relax_methoxy_local_group_states = (
                    self.relax_methoxy_local_group_states
                )
                im_database_driver.methoxy_state_opt_max_iter = (
                    self.methoxy_state_opt_max_iter
                )

            im_database_driver.system_from_molecule(dynamics_molecule, self.roots_z_matrix, forcefield_generator, solvent=self.solvent, qm_atoms='all')
            if self.bias_force_reaction_prop is not None:
                im_database_driver.bias_force_reaction_idx = self.bias_force_reaction_idx
                im_database_driver.bias_force_reaction_prop = self.bias_force_reaction_prop

            density_of_datapoints = self.determine_datapoint_density(self.states_interpolation_settings)
            self.density_of_datapoints = density_of_datapoints
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
                self.point_added_molecules.append(im_database_driver.point_adding_molecule)
                self.unique_molecules.append(im_database_driver.allowed_molecules)

                self._confirm_database_quality(molecule, basis=states_basis, im_settings=self.states_interpolation_settings, given_molecular_strucutres=self.state_specific_molecules)

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
        Build one global Z-matrix for a reaction path by merging all unique
        internal coordinates found in provided reference structures.

        Parameters
        ----------
        reaction_structures : list
            Reaction-path structures. Supported entry types:
            - Molecule
            - XYZ string
            - tuple/list where first item is Molecule or XYZ string
            - dict with key "molecule"

        root : int, default=0
            Root index whose `self.roots_z_matrix[root]` is used/updated.

        include_existing_root_zmat : bool, default=True
            If True and root z-matrix exists, start from it and append missing
            coordinates from reaction structures.

        reference_molecule : Molecule | None, default=None
            Reference for charge/multiplicity and atom consistency checks.
            If None: uses `self.molecule`, else first reaction structure.

        forced_coordinates : dict | None, default=None
            Coordinates that must be present in final z-matrix.
            Accepted keys (singular/plural):
            - bond / bonds
            - angle / angles
            - dihedral / dihedrals
            - improper / impropers

            Example (0-based):
                {
                    "bonds": [(donor, H), (acceptor, H)],
                    "angles": [(donor, H, acceptor)],
                }

        enforce_same_atoms : bool, default=True
            If True, checks same atom count and same atom-label order
            across all reaction structures.

        Returns
        -------
        dict
            {
                "global_z_matrix": dict,
                "per_structure_z_matrices": list[dict],
                "added_coordinates": dict,
                "added_counts": dict,
            }
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
        Return seed structures in conformer-compatible layout:
            {state: {dih_key: [(mol, mode, constraints), ...]}}

        - mode: 'normal' | 'transition' | 'constraint'
        - constraints: list of geomeTRIC-compatible constraints
        (tuple constraints are assumed 1-based, matching current conformer flow).
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
            if isinstance(mol_like, Molecule):
                mol = Molecule.from_xyz_string(mol_like.get_xyz_string())
            elif isinstance(mol_like, str):
                mol = Molecule.from_xyz_string(mol_like)
            else:
                assert_msg_critical(False, "build_initial_seed_structures: expected Molecule or XYZ string.")
            mol.set_charge(molecule.get_charge())
            mol.set_multiplicity(molecule.get_multiplicity())
            return mol

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

        reseted_point_densities_dict = {state: 0 for state in self.roots_to_follow}

        for state in reseted_point_densities_dict.keys():
            qm_datapoints = []
            if Path(imforcefieldfile[state]['imforcefield_file']).exists():
                impes_driver = InterpolationDriver(self.roots_z_matrix[state])
                impes_driver.update_settings(imforcefieldfile[state])
                self.qmlabels, z_matrix = impes_driver.read_labels()

                for label in self.qmlabels:
                    if "core" in label:
                        qm_data_point = InterpolationDatapoint(z_matrix)
                        qm_data_point.update_settings(imforcefieldfile[state])
                        qm_data_point.read_hdf5(imforcefieldfile[state]['imforcefield_file'], label)

                        qm_datapoints.append(qm_data_point)
                        reseted_point_densities_dict[state] += 1

        return reseted_point_densities_dict

    def _calculate_translation_coordinates_analysis(self, given_coordinates):
        """Center the molecule by translating its geometric center to (0, 0, 0).

           :param given_coordinates:
                Coordinate that is translated to the center.
        """
        center = np.mean(given_coordinates, axis=0)
        translated_coordinates = given_coordinates - center

        return translated_coordinates

    def calculate_distance_to_ref(self, current_coordinates, datapoint_coordinate):
        """Calculates and returns the cartesian distance between
           self.coordinates and data_point coordinates.
           Besides the distance, it also returns the weight gradient,
           which requires the distance vector to be computed.

           :param current_coordinates:
                current molecular coordinates.

           :param data_point:
                InterpolationDatapoint object.

           :returns:
              Norm of the distance between 2 structures.
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
        """Extracts molecular structures from a given database file.

        :param datafile:
            Database file containing interpolation data.

        :param mol_labels:
            List of molecular labels.

        :returns:
            A list of VeloxChem Molecule objects extracted from the database.

        """

        im_driver = InterpolationDriver()  # -> implemented Class in VeloxChem that is capable to perform interpolation calculations for a given molecule and provided z_matrix and database
        im_settings['imforcefield_file'] = datafile
        im_driver.update_settings(im_settings)
        # im_driver.imforcefield_file = datafile
        labels, z_matrix = im_driver.read_labels()
        im_driver.qm_local_factor_banks = im_driver._load_local_factor_banks(
            z_matrix)

        sorted_labels = sorted(labels, key=lambda x: int(x.split('_')[1]))

        # impes_coordinate = InterpolationDatapoint(z_matrix)  # -> implemented Class in VeloxChem that handles all transformations and database changes concerning the interpolation
        data_point_molecules = []
        datapoints = []

        for label in sorted_labels:
            impes_coordinate = InterpolationDatapoint(z_matrix)
            impes_coordinate.update_settings(im_settings)
            impes_coordinate.read_hdf5(datafile, label)  # -> read in function from the ImpesDriver object

            if im_driver._is_local_factor_cluster_datapoint(
                    impes_coordinate):
                continue

            coordinates_in_angstrom = impes_coordinate.cartesian_coordinates * bohr_in_angstrom()

            current_molecule = Molecule(mol_labels, coordinates_in_angstrom, 'angstrom')  # -> creates a VeloxChem Molecule object

            datapoints.append(impes_coordinate)
            data_point_molecules.append(current_molecule)

        return data_point_molecules, datapoints, z_matrix

    def _confirm_database_quality(self, molecule, basis, im_settings, given_molecular_strucutres=None, improve=True):
        """Validates the quality of an interpolation database for a given molecule.

       This function assesses the quality of the provided interpolation database
       comparing the interpolated energy with a QM-reference energy.

       :param molecule:
           A VeloxChem molecule object representing the reference molecular system.

       :param im_database_file:
           Interpolation database file.

       :param given_molecular_strucutres:
           An optional list of additional molecular structures that will be used for the validation.

       :returns:
           List of QM-energies, IM-energies.
        """

        def _calculate_translation_coordinates(given_coordinates):
            """Center the molecule by translating its geometric center to (0, 0, 0)."""
            center = np.mean(given_coordinates, axis=0)
            translated_coordinates = given_coordinates - center

            return translated_coordinates

        def cartesian_just_distance(coordinate_1, coordinate_2, non_core_atoms=[]):
            """Calculates and returns the cartesian distance between
                self.coordinates and data_point coordinates.
                Besides the distance, it also returns the weight gradient,
                which requires the distance vector to be computed.
                :param data_point:
                        InterpolationDatapoint object
                """
            target_coordinates_core = np.delete(coordinate_1, non_core_atoms, axis=0)
            reference_coordinates_core = np.delete(coordinate_2, non_core_atoms, axis=0)
            # First, translate the cartesian coordinates to zero
            target_coordinates = _calculate_translation_coordinates(target_coordinates_core)
            reference_coordinates = (
                _calculate_translation_coordinates(reference_coordinates_core))
            # Then, determine the rotation matrix which
            # aligns data_point (target_coordinates)
            # to self.impes_coordinate (reference_coordinates)
            rotation_matrix = geometric.rotate.get_rot(
                target_coordinates,
                reference_coordinates)
            # Rotate the data point
            rotated_coordinates = np.dot(rotation_matrix, target_coordinates.T).T
            # Calculate the Cartesian distance
            distance_vector = (reference_coordinates - rotated_coordinates)

            return np.linalg.norm(distance_vector)

        def database_distance_check(datafile_mols):

            point_dists = {}

            for mol_idx, mol_1 in enumerate(datafile_mols[:-1]):

                single_p_distances = {}
                for mol_2_idx, mol_2 in enumerate(datafile_mols[mol_idx + 1:], start=mol_idx + 1):
                    distance = cartesian_just_distance(mol_1.get_coordinates_in_bohr(), mol_2.get_coordinates_in_bohr())

                    single_p_distances[mol_2_idx] = distance

                point_dists[mol_idx] = single_p_distances

            return point_dists

        def dist_dict_to_edges(db_distances: dict):
            """
            Convert upper-triangular dict-of-dicts distances into:
            N (number of points), pairs (M,2), dists (M,)
            """
            idxs = set(db_distances.keys())
            if len(idxs) < 2:
                return 1, np.array([]), np.array([])
            for i, row in db_distances.items():
                for j in row.keys():
                    idxs.add(j)
            N = max(idxs) + 1

            pairs = []
            dists = []
            for i, row in db_distances.items():
                for j, dij in row.items():
                    pairs.append((int(i), int(j)))
                    dists.append(float(dij))

            pairs = np.asarray(pairs, dtype=np.int32)
            dists = np.asarray(dists, dtype=np.float64)
            return N, pairs, dists

        overall_db_covergage = {}

        if given_molecular_strucutres is None:
            raise ValueError('confirm_database_quality requires given_molecular_strucutres.')

        if not isinstance(given_molecular_strucutres, dict):
            raise ValueError('confirm_database_quality expects given_molecular_strucutres to be a dict keyed by root.')

        for root in self.roots_to_follow:
            database_quality = False
            drivers = self.drivers['gs'] if root == 0 else self.drivers['es']
            current_datafile = im_settings[root]['imforcefield_file']
            all_structures_root = []
            random_structure_choices_root = []
            last_qm_energies = []
            last_im_energies = []

            while not database_quality:
                # payload = None

                all_structures_root = list(given_molecular_strucutres.get(root, []))
                if len(all_structures_root) == 0:
                    database_quality = True
                    continue
                    # payload['skip_root'] = True
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
                    self.ostream.print_header(f"Root {root}: Checking if current structures are well seperated from the database conformations.")
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
                                'The current structures are well seperated from the database conformations! '
                                'loop is discontinued')
                            self.ostream.print_blank()
                            self.ostream.flush()

                            dist_ok = True

                    random_structure_choices_root = {"random_struct_info":list(selected_molecules)}

                if random_structure_choices_root is None:
                    raise RuntimeError('confirm_database_quality: failed to broadcast root payload.')

                datapoint_molecules, _, _ = self.database_extracter(
                    current_datafile,
                    molecule.get_labels(),
                    im_settings[root],
                )

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
                impes_driver.qm_local_factor_banks = (
                    impes_driver._load_local_factor_banks(
                        self.roots_z_matrix[root], root=root)
                )
                impes_driver.qm_data_points = []
                kept_labels = []

                for label in im_labels:
                    qm_data_point = InterpolationDatapoint(self.roots_z_matrix[root])
                    qm_data_point.update_settings(im_settings[root])
                    qm_data_point.read_hdf5(current_datafile, label)
                    qm_data_point.inv_sqrt_masses = inv_sqrt_masses

                    if impes_driver._is_local_factor_cluster_datapoint(
                            qm_data_point):
                        continue

                    impes_driver.qm_data_points.append(qm_data_point)
                    kept_labels.append(label)

                if not impes_driver.qm_data_points:
                    raise RuntimeError(
                        'Database quality evaluation found no outer datapoints.')

                impes_driver.labels = kept_labels
                impes_driver.impes_coordinate.eq_bond_lengths = (
                    impes_driver.qm_data_points[0].eq_bond_lengths)
                impes_driver.install_local_factor_banks(
                    impes_driver.qm_local_factor_banks)

                qm_energies = []
                im_energies = []
                database_expanded = False

                for i, mol in enumerate(random_structure_choices_root.get('random_struct_info')[:]):
                    if root == 0:
                        current_basis = self._read_basis_for_driver(
                            mol, drivers[0], basis['gs'])
                    else:
                        current_basis = self._read_basis_for_driver(
                            mol, drivers[0], basis['es'])
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
                        self.ostream.print_block("The current structure is not within the desired threshold and can be added to the database if desired! \n The structure can be found in the random_structure.xyz file")
                        self.ostream.flush()
                    else:

                        with h5py.File('summary_output.h5', "a") as h5f:
                            self._append_confirm_database_quality_h5(
                                h5f,
                                [mol],
                                np.array([qm_energy_val]),
                                np.array([impes_driver.impes_coordinate.energy]),
                                root,
                            )

                last_qm_energies = qm_energies
                last_im_energies = im_energies

                if not database_expanded:
                    database_quality = True

                    overall_db_covergage[current_datafile] = {}
                    minmum_distances_in_db = database_distance_check(datapoint_molecules)
                    overall_db_covergage[current_datafile]['db_distances'] = minmum_distances_in_db
                    # N, pairs, dists = dist_dict_to_edges(minmum_distances_in_db)
                    # labels = [f"point {i+1}" for i in range(N)]
                    dist_dict_to_edges(minmum_distances_in_db)

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

    def _write_string_dataset(self, h5f, name, value):
        if value is None:
            return
        dt = h5py.string_dtype(encoding="utf-8")
        h5f.create_dataset(name, data=np.array(value, dtype=object), dtype=dt)

    def _local_group_cluster_coupling_rows(self, local_group_model, cluster):
        rows = []
        for rotor_id in cluster.rotor_ids:
            rotor = local_group_model.rotors.get(str(rotor_id))
            if rotor is None:
                continue
            rows.extend(int(row) for row in rotor.torsion_rows)

        if not rows:
            rows.extend(int(row) for row in cluster.active_rows)

        return tuple(sorted(set(rows)))

    def _local_group_coupling_score(
            self,
            internal_hessian,
            cluster_a_rows,
            cluster_b_rows,
            eps=1.0e-12):
        hessian = np.asarray(internal_hessian, dtype=np.float64)
        if hessian.ndim != 2 or hessian.shape[0] != hessian.shape[1]:
            return 0.0

        n_rows = hessian.shape[0]
        rows_a = tuple(
            row for row in sorted(set(int(x) for x in cluster_a_rows))
            if 0 <= row < n_rows
        )
        rows_b = tuple(
            row for row in sorted(set(int(x) for x in cluster_b_rows))
            if 0 <= row < n_rows
        )

        if not rows_a or not rows_b:
            return 0.0

        H_ab = hessian[np.ix_(rows_a, rows_b)]
        H_aa = hessian[np.ix_(rows_a, rows_a)]
        H_bb = hessian[np.ix_(rows_b, rows_b)]

        numerator = np.linalg.norm(H_ab, ord="fro")
        denominator = np.sqrt(
            np.linalg.norm(H_aa, ord="fro")
            * np.linalg.norm(H_bb, ord="fro")
            + eps
        )

        if denominator <= eps:
            return 0.0

        return float(numerator / denominator)

    def _local_group_clusters_are_nested_rotors(
            self,
            local_group_model,
            cluster_a,
            cluster_b):
        if len(cluster_a.rotor_ids) != 1 or len(cluster_b.rotor_ids) != 1:
            return False

        rotor_a = local_group_model.rotors.get(str(cluster_a.rotor_ids[0]))
        rotor_b = local_group_model.rotors.get(str(cluster_b.rotor_ids[0]))
        if rotor_a is None or rotor_b is None:
            return False

        owned_a = set(int(atom) for atom in rotor_a.owned_atoms)
        owned_b = set(int(atom) for atom in rotor_b.owned_atoms)
        if not owned_a or not owned_b:
            return False

        return owned_a < owned_b or owned_b < owned_a

    def _order_local_group_rotor_ids(self, local_group_model, rotor_ids):
        def rotor_key(rotor_id):
            rotor = local_group_model.rotors.get(str(rotor_id))
            if rotor is None:
                return (0, str(rotor_id))
            return (-len(tuple(rotor.owned_atoms)), str(rotor_id))

        return tuple(str(rotor_id) for rotor_id in sorted(rotor_ids, key=rotor_key))

    def _local_group_maximal_direct_factors(self, adjacency):
        maximal_factors = []

        def bron_kerbosch(current, candidates, excluded):
            if not candidates and not excluded:
                maximal_factors.append(tuple(sorted(current)))
                return

            pivot_pool = candidates | excluded
            if pivot_pool:
                pivot = max(
                    pivot_pool,
                    key=lambda node: len(candidates & adjacency[node]),
                )
                search_nodes = candidates - adjacency[pivot]
            else:
                search_nodes = set(candidates)

            for node in sorted(search_nodes):
                bron_kerbosch(
                    current | {node},
                    candidates & adjacency[node],
                    excluded & adjacency[node],
                )
                candidates.remove(node)
                excluded.add(node)

        bron_kerbosch(set(), set(adjacency), set())

        return tuple(sorted(
            maximal_factors,
            key=lambda factor: (min(factor), -len(factor), factor),
        ))

    def _local_group_connected_components(self, adjacency):
        unvisited = set(adjacency)
        components = []

        while unvisited:
            seed = min(unvisited)
            unvisited.remove(seed)
            component = {seed}
            stack = [seed]

            while stack:
                current = stack.pop()
                neighbors = sorted(adjacency[current] & unvisited)
                for neighbor in neighbors:
                    unvisited.remove(neighbor)
                    component.add(neighbor)
                    stack.append(neighbor)

            components.append(tuple(sorted(component)))

        return tuple(components)

    def _local_group_factor_active_rows(
            self,
            local_group_model,
            rotor_ids,
            fallback_clusters):
        active_rows = []
        for cluster in fallback_clusters:
            active_rows.extend(int(row) for row in cluster.active_rows)

        for rotor_id in rotor_ids:
            rotor = local_group_model.rotors.get(str(rotor_id))
            if rotor is not None:
                active_rows.extend(int(row) for row in rotor.torsion_rows)

        return tuple(sorted(set(active_rows)))

    def _local_group_clusters_must_merge(
            self,
            local_group_model,
            cluster_a,
            cluster_b):
        if not self.local_group_force_merge_overlapping_groups:
            return False

        if self._local_group_clusters_are_nested_rotors(
                local_group_model, cluster_a, cluster_b):
            return False

        owned_a = set(int(atom) for atom in cluster_a.owned_atoms)
        owned_b = set(int(atom) for atom in cluster_b.owned_atoms)
        if owned_a & owned_b:
            return True

        rows_a = set(int(row) for row in cluster_a.active_rows)
        rows_b = set(int(row) for row in cluster_b.active_rows)
        if rows_a & rows_b:
            return True

        return False

    def _local_group_clusters_coupling_is_suppressed(
            self,
            local_group_model,
            cluster_a,
            cluster_b):
        if getattr(
                self,
                "local_group_allow_anchored_linker_methyl_coupling",
                False):
            return False

        a_is_linker = self._local_group_cluster_contains_rotor_kind(
            local_group_model, cluster_a, 'anchored_linker')
        b_is_linker = self._local_group_cluster_contains_rotor_kind(
            local_group_model, cluster_b, 'anchored_linker')
        a_is_methyl = self._local_group_cluster_contains_rotor_kind(
            local_group_model, cluster_a, 'methyl')
        b_is_methyl = self._local_group_cluster_contains_rotor_kind(
            local_group_model, cluster_b, 'methyl')

        return (a_is_linker and b_is_methyl) or (b_is_linker and a_is_methyl)

    def _build_coupled_local_group_model(
            self,
            local_group_model,
            internal_hessian,
            threshold=None):
        threshold = (
            self.local_group_coupling_threshold
            if threshold is None else float(threshold)
        )

        cluster_items = sorted(
            local_group_model.clusters.items(),
            key=lambda item: str(item[0]),
        )
        cluster_ids = tuple(str(cluster_id) for cluster_id, _ in cluster_items)
        n_clusters = len(cluster_items)

        coupling_matrix = np.eye(n_clusters, dtype=np.float64)
        forced_merge_matrix = np.zeros((n_clusters, n_clusters), dtype=bool)

        coupling_rows = {
            str(cluster_id): self._local_group_cluster_coupling_rows(
                local_group_model, cluster)
            for cluster_id, cluster in cluster_items
        }

        for i in range(n_clusters):
            for j in range(i + 1, n_clusters):
                cluster_i = cluster_items[i][1]
                cluster_j = cluster_items[j][1]
                score = self._local_group_coupling_score(
                    internal_hessian,
                    coupling_rows[cluster_ids[i]],
                    coupling_rows[cluster_ids[j]],
                )
                coupling_matrix[i, j] = score
                coupling_matrix[j, i] = score

                forced_merge = self._local_group_clusters_must_merge(
                    local_group_model, cluster_i, cluster_j)
                forced_merge_matrix[i, j] = forced_merge
                forced_merge_matrix[j, i] = forced_merge

        adjacency = {idx: set() for idx in range(n_clusters)}
        direct_edges = []
        for i in range(n_clusters):
            for j in range(i + 1, n_clusters):
                is_direct_edge = (
                    (
                        forced_merge_matrix[i, j]
                        or coupling_matrix[i, j] >= threshold
                    )
                    and not self._local_group_clusters_coupling_is_suppressed(
                        local_group_model,
                        cluster_items[i][1],
                        cluster_items[j][1],
                    )
                )
                if not is_direct_edge:
                    continue

                adjacency[i].add(j)
                adjacency[j].add(i)
                direct_edges.append((
                    cluster_ids[i],
                    cluster_ids[j],
                    float(coupling_matrix[i, j]),
                    bool(forced_merge_matrix[i, j]),
                ))

        if self.local_group_coupling_topology == 'connected_components':
            factors = self._local_group_connected_components(adjacency)
        else:
            factors = self._local_group_maximal_direct_factors(adjacency)

        coupled_clusters = {}
        for factor_index, factor in enumerate(factors):
            component_clusters = [cluster_items[idx][1] for idx in factor]

            if len(component_clusters) == 1:
                cluster = component_clusters[0]
                coupled_clusters[str(cluster.cluster_id)] = cluster
                continue

            rotor_ids = []
            rotor_id_set = set()
            owned_atoms = set()
            active_atoms = set()
            anchor_atoms = set()
            relaxation_atoms = set()

            for cluster in component_clusters:
                owned_atoms.update(int(atom) for atom in cluster.owned_atoms)
                active_atoms.update(int(atom) for atom in cluster.active_atoms)
                anchor_atoms.update(
                    int(atom)
                    for atom in getattr(cluster, "anchor_atoms", ())
                )
                relaxation_atoms.update(
                    int(atom)
                    for atom in getattr(cluster, "relaxation_atoms", ())
                )
                for rotor_id in cluster.rotor_ids:
                    rotor_id = str(rotor_id)
                    if rotor_id not in rotor_id_set:
                        rotor_id_set.add(rotor_id)
                        rotor_ids.append(rotor_id)

            rotor_ids = self._order_local_group_rotor_ids(
                local_group_model, rotor_ids)
            active_rows = self._local_group_factor_active_rows(
                local_group_model,
                rotor_ids,
                component_clusters,
            )
            cluster_id = f"coupled_{factor_index}"
            coupled_clusters[cluster_id] = LocalCluster(
                cluster_id=cluster_id,
                family_label='coupled',
                rotor_ids=rotor_ids,
                owned_atoms=tuple(sorted(owned_atoms)),
                active_atoms=tuple(sorted(active_atoms)),
                active_rows=active_rows,
                anchor_atoms=tuple(sorted(anchor_atoms)),
                relaxation_atoms=tuple(sorted(relaxation_atoms)),
            )

        coupled_model = LocalGroupModel(
            version=local_group_model.version,
            enabled=local_group_model.enabled,
            rotors=local_group_model.rotors,
            clusters=coupled_clusters,
            core_atoms=local_group_model.core_atoms,
            core_rows=local_group_model.core_rows,
        )

        coupling_data = {
            'primitive_cluster_ids': cluster_ids,
            'primitive_cluster_types': tuple(
                str(cluster.family_label) for _, cluster in cluster_items
            ),
            'primitive_rotor_ids': tuple(
                tuple(str(rotor_id) for rotor_id in cluster.rotor_ids)
                for _, cluster in cluster_items
            ),
            'coupled_cluster_ids': tuple(coupled_clusters.keys()),
            'coupling_topology': self.local_group_coupling_topology,
            'direct_edges': tuple(direct_edges),
            'coupling_matrix': coupling_matrix,
            'forced_merge_matrix': forced_merge_matrix,
            'threshold': float(threshold),
            'components': tuple(
                tuple(cluster_ids[idx] for idx in factor)
                for factor in factors
            ),
        }

        return coupled_model, coupling_data

    def _schema3_relaxed_residual_enabled(self):
        mode = getattr(
            self,
            "local_factor_combination_mode",
            COMBINATION_SIGNED_FULL,
        )
        return str(mode).strip().lower().replace(
            "-", "_") == COMBINATION_SIGNED_RELAXED_RESIDUAL

    def _schema4_msi_enabled(self):
        version = int(getattr(self, "msi_schema_version", 3))
        if version not in (2, 3, 4):
            raise ValueError("msi_schema_version must be 2, 3, or 4.")
        return version == 4

    def _schema4_environment_descriptor_spec(self, root):
        if not self._schema4_msi_enabled():
            raise RuntimeError(
                "Schema-4 descriptor construction requires msi_schema_version=4."
            )
        molecule = getattr(self, "molecule", None)
        if molecule is None:
            raise RuntimeError(
                "Schema-4 construction requires the setup molecule."
            )
        local_group_model = self.local_group_model.get(
            root, self.local_group_primitive_model.get(root)
        )
        if local_group_model is None:
            raise RuntimeError(
                f"Schema-4 construction has no local-group model for root {root}."
            )
        policy = dict(self.msi_descriptor_policy)
        return build_contact_torsion_v1_spec(
            descriptor_spec_id=f"descriptor.contact_torsion.root_{root}.v1",
            labels=molecule.get_labels(),
            connectivity=molecule.get_connectivity_matrix(),
            local_group_model=local_group_model,
            z_matrix=self.roots_z_matrix[root],
            distance_coefficient=float(policy["distance_coefficient"]),
            contact_r0_angstrom=float(policy["contact_r0_angstrom"]),
            contact_exponent=int(policy["contact_exponent"]),
        )

    def construct_schema4_registry(
            self,
            target_file,
            root,
            interpolation_settings,
            descriptor_spec=None,
            legacy_banks=None):
        """Audit computed factor states, promote cores, and atomically write MSI."""

        if not self._schema4_msi_enabled():
            raise RuntimeError(
                "construct_schema4_registry requires msi_schema_version=4; "
                "schemas 2/3 are intentionally left on their legacy path."
            )
        if descriptor_spec is None:
            descriptor_spec = self._schema4_environment_descriptor_spec(root)
        if legacy_banks is None:
            legacy_banks = load_signed_factor_banks_for_root(
                target_file,
                root,
                self.roots_z_matrix[root],
                interpolation_settings,
            )
        constructor = Schema4MSIConstructor(
            locality_threshold_policy=self.msi_locality_threshold_policy,
            promotion_policy=self.msi_promotion_policy,
        )
        registry, audits, promotions = constructor.construct(
            legacy_banks=legacy_banks,
            root=root,
            descriptor_spec=descriptor_spec,
            output_file=target_file,
        )
        self.msi_construction_stage_trace = tuple(constructor.stage_trace)
        return registry, audits, promotions

    @staticmethod
    def _active_projector_rows(cluster, rotor_map):
        """Cluster-active rows plus every physical signature/torsion row."""

        rows = set(int(row) for row in cluster.active_rows)
        for rotor_id in cluster.rotor_ids:
            rotor = rotor_map.get(str(rotor_id))
            if rotor is None:
                continue
            rows.update(int(row) for row in rotor.signature_rows)
            rows.update(int(row) for row in rotor.torsion_rows)

        return tuple(sorted(rows))

    def _with_schema3_projector_rows(self, local_group_model):
        """Attach explicit support only to direct anchored-linker factors."""

        if not self._schema3_relaxed_residual_enabled():
            return local_group_model

        clusters = {}
        changed = False
        for cluster_id, cluster in local_group_model.clusters.items():
            if str(cluster.family_label) != "anchored_linker":
                clusters[str(cluster_id)] = cluster
                continue

            projector_rows = self._active_projector_rows(
                cluster, local_group_model.rotors)
            cluster = replace(
                cluster,
                projector_rows=projector_rows,
                projector_policy_id="active_rows_plus_signatures",
            )
            clusters[str(cluster_id)] = cluster
            changed = True

        if not changed:
            return local_group_model

        return replace(local_group_model, clusters=clusters)

    def _local_factor_overlap_cluster_id(self, rotor_key, clusters):
        safe_parts = []
        for rotor_id in rotor_key:
            safe = "".join(
                char if char.isalnum() else "_"
                for char in str(rotor_id)
            ).strip("_")
            safe_parts.append(safe or "rotor")
        base = "overlap_" + "_".join(safe_parts)
        cluster_id = base
        counter = 1
        while cluster_id in clusters:
            cluster_id = f"{base}_{counter}"
            counter += 1
        return cluster_id

    def _augment_schema3_relaxed_residual_clusters(self, local_group_model):
        """Adds explicit canonical overlap clusters required by schema 3."""

        if not self._schema3_relaxed_residual_enabled():
            return local_group_model, ()

        direct_items = sorted(
            (
                (str(cluster_id), cluster)
                for cluster_id, cluster in local_group_model.clusters.items()
                if str(getattr(cluster, "role", "factor")).lower()
                != "overlap"
            ),
            key=lambda item: str(item[0]),
        )
        coefficient_by_key, _ = build_intersection_coefficients(
            (cluster_id, cluster.rotor_ids)
            for cluster_id, cluster in direct_items)

        clusters = dict(local_group_model.clusters)
        key_to_cluster_id = {}
        for cluster_id, cluster in clusters.items():
            key = tuple(sorted(
                str(rotor_id)
                for rotor_id in (
                    getattr(cluster, "canonical_subset_key", ())
                    or getattr(cluster, "rotor_ids", ())
                )
            ))
            key_to_cluster_id[key] = str(cluster_id)

        added_cluster_ids = []
        for rotor_key, coefficient in sorted(coefficient_by_key.items()):
            if abs(float(coefficient)) <= 1.0e-12:
                continue
            if rotor_key in key_to_cluster_id:
                continue

            rotor_set = set(str(rotor_id) for rotor_id in rotor_key)
            parents = [
                (cluster_id, cluster)
                for cluster_id, cluster in direct_items
                if rotor_set.issubset(
                    set(str(rotor_id) for rotor_id in cluster.rotor_ids))
            ]
            if not parents:
                continue

            active_row_sets = [
                set(int(row) for row in cluster.active_rows)
                for _, cluster in parents
                if getattr(cluster, "active_rows", ())
            ]
            active_rows = (
                set.intersection(*active_row_sets)
                if active_row_sets else set()
            )
            owned_atoms = set()
            active_atoms = set()
            anchor_atoms = set()
            relaxation_atoms = set()

            for rotor_id in rotor_key:
                rotor = local_group_model.rotors.get(str(rotor_id))
                if rotor is None:
                    continue
                active_rows.update(int(row) for row in rotor.torsion_rows)
                owned_atoms.update(int(atom) for atom in rotor.owned_atoms)

            for _, cluster in parents:
                parent_active_atoms = set(
                    int(atom) for atom in getattr(cluster, "active_atoms", ()))
                if active_atoms:
                    active_atoms.intersection_update(parent_active_atoms)
                else:
                    active_atoms.update(parent_active_atoms)
                anchor_atoms.update(
                    int(atom)
                    for atom in getattr(cluster, "anchor_atoms", ())
                )
                relaxation_atoms.update(
                    int(atom)
                    for atom in getattr(cluster, "relaxation_atoms", ())
                )

            active_atoms.update(owned_atoms)
            cluster_id = self._local_factor_overlap_cluster_id(
                rotor_key, clusters)
            clusters[cluster_id] = LocalCluster(
                cluster_id=cluster_id,
                family_label="canonical_overlap",
                rotor_ids=tuple(str(rotor_id) for rotor_id in rotor_key),
                owned_atoms=tuple(sorted(owned_atoms)),
                active_atoms=tuple(sorted(active_atoms)),
                active_rows=tuple(sorted(active_rows)),
                anchor_atoms=tuple(sorted(anchor_atoms)),
                relaxation_atoms=tuple(sorted(relaxation_atoms)),
                role="overlap",
                canonical_subset_key=tuple(str(rotor_id)
                                           for rotor_id in rotor_key),
                parent_cluster_ids=tuple(
                    str(parent_id) for parent_id, _ in parents),
                relaxation_policy_id="default",
            )
            key_to_cluster_id[rotor_key] = cluster_id
            added_cluster_ids.append(cluster_id)

        if not added_cluster_ids:
            return local_group_model, ()

        return (
            LocalGroupModel(
                version=local_group_model.version,
                enabled=local_group_model.enabled,
                rotors=local_group_model.rotors,
                clusters=clusters,
                core_atoms=local_group_model.core_atoms,
                core_rows=local_group_model.core_rows,
            ),
            tuple(added_cluster_ids),
        )

    def _with_schema3_anchor_state_ids(self, local_group_model):
        """Returns a schema-3 model whose clusters declare anchor state ids."""

        if not self._schema3_relaxed_residual_enabled():
            return local_group_model

        clusters = {}
        changed = False
        for cluster_id, cluster in local_group_model.clusters.items():
            anchor_state_ids = tuple(
                int(state_id)
                for state_id in (
                    getattr(cluster, "anchor_state_ids", ()) or (0,)
                )
            )
            if anchor_state_ids != getattr(cluster, "anchor_state_ids", ()):
                cluster = replace(
                    cluster,
                    anchor_state_ids=anchor_state_ids,
                )
                changed = True
            clusters[cluster_id] = cluster

        if not changed:
            return local_group_model

        return LocalGroupModel(
            version=local_group_model.version,
            enabled=local_group_model.enabled,
            rotors=local_group_model.rotors,
            clusters=clusters,
            core_atoms=local_group_model.core_atoms,
            core_rows=local_group_model.core_rows,
        )

    def _local_group_phase_values(self, rotor):
        if getattr(rotor, "phase_values", None):
            return tuple(float(x) for x in rotor.phase_values)

        symmetry_order = max(int(getattr(rotor, 'symmetry_order', 1)), 1)
        if symmetry_order in self.local_group_phase_library:
            return tuple(float(x) for x in self.local_group_phase_library[symmetry_order])
        return (0.0,)

    def _oriented_local_rotor_dihedral(self, rotor):
        if rotor.phase_coordinate is None:
            return None

        coord = tuple(int(x) for x in rotor.phase_coordinate)
        axis = tuple(int(x) for x in rotor.axis)

        if tuple(coord[1:3]) == axis:
            return coord
        if tuple(coord[2:0:-1]) == axis:
            return tuple(reversed(coord))

        # Last-resort orientation by central-bond set. This keeps construction
        # possible if geomeTRIC returns the same axis in an unexpected ordering.
        if set(coord[1:3]) == set(axis):
            if coord[1] == axis[0]:
                return coord
            return tuple(reversed(coord))

        return coord

    def _clone_molecule_for_local_group(self, molecule, coordinates_bohr=None):
        if coordinates_bohr is None:
            coordinates_bohr = molecule.get_coordinates_in_bohr()

        cloned = Molecule(
            molecule.get_labels(),
            np.array(coordinates_bohr, dtype=np.float64, copy=True),
            'bohr',
        )
        cloned.set_charge(molecule.get_charge())
        cloned.set_multiplicity(molecule.get_multiplicity())
        return cloned

    def _apply_local_group_state(
            self,
            molecule,
            local_group_model,
            rotor_ids,
            phase_signature,
            skip_kinds=None):
        rotated = self._clone_molecule_for_local_group(molecule)
        skip_kinds = set(str(kind) for kind in (skip_kinds or ()))

        for rotor_id, phase in zip(rotor_ids, phase_signature):
            phase = float(phase)
            if abs(phase) < 1.0e-14:
                continue

            rotor = local_group_model.rotors[str(rotor_id)]
            if str(rotor.kind) in skip_kinds:
                continue

            dihedral = self._oriented_local_rotor_dihedral(rotor)
            if dihedral is None:
                continue

            dihedral_one_based = [int(x) + 1 for x in dihedral]
            current_angle = rotated.get_dihedral(dihedral_one_based, 'radian')
            rotated.set_dihedral(
                dihedral_one_based,
                current_angle + phase,
                'radian',
                verbose=False,
            )

        return rotated

    def _local_factor_relaxation_constraints(
            self,
            molecule,
            anchor_molecule,
            z_matrix,
            local_group_model,
            cluster,
            rotor_ids,
            phase_signature):
        """
        Builds schema-3 canonical local-factor relaxation constraints.

        The constraints intentionally do not freeze Cartesian atoms.  They hold
        the sampled subset rotors at their requested phases and anchor every
        other local rotor at its family/core phase, while allowing the rest of
        the molecule, including core atoms, to relax as response coordinates.
        """

        del z_matrix, cluster

        if anchor_molecule is None:
            anchor_molecule = molecule

        selected_phases = {
            str(rotor_id): float(phase)
            for rotor_id, phase in zip(rotor_ids, phase_signature)
        }
        anchor_external = bool(
            getattr(self, "local_factor_anchor_external_rotors", True)
        )

        constraints = []
        seen = set()
        # Selected subset constraints take precedence if two rotor records use
        # equivalent oriented dihedrals.  This also makes duplicate physical
        # constraints deterministic instead of depending on rotor-id sorting.
        rotor_items = sorted(
            local_group_model.rotors.items(),
            key=lambda item: (
                str(item[0]) not in selected_phases,
                str(item[0]),
            ),
        )
        for rotor_id, rotor in rotor_items:
            rotor_id = str(rotor_id)
            if rotor_id not in selected_phases and not anchor_external:
                continue

            dihedral = self._oriented_local_rotor_dihedral(rotor)
            if dihedral is None:
                continue

            key = tuple(int(atom) for atom in dihedral)
            canonical_key = min(key, tuple(reversed(key)))
            if canonical_key in seen:
                continue
            seen.add(canonical_key)

            phase = selected_phases.get(rotor_id, 0.0)
            dihedral_one_based = [int(atom) + 1 for atom in key]
            anchor_angle = anchor_molecule.get_dihedral(
                dihedral_one_based, 'degree')
            target_angle = float(anchor_angle) + np.rad2deg(float(phase))
            constraints.append(
                "set dihedral "
                + " ".join(str(atom) for atom in dihedral_one_based)
                + f" {target_angle:.12f}"
            )

        return constraints

    def _relax_canonical_local_factor_state(
            self,
            drivers,
            molecule,
            basis,
            z_matrix,
            local_group_model,
            cluster,
            rotor_ids,
            phase_signature,
            anchor_molecule=None):
        if not getattr(self, "relax_local_factor_states", True):
            return molecule, basis

        constraints = self._local_factor_relaxation_constraints(
            molecule,
            anchor_molecule,
            z_matrix,
            local_group_model,
            cluster,
            rotor_ids,
            phase_signature,
        )
        if not constraints:
            return molecule, basis

        max_iter = getattr(self, "local_factor_state_opt_max_iter", 120)
        if self._is_xtb_like_driver(drivers[0]):
            opt_results = self._run_optimization(
                drivers[0],
                molecule,
                constraints=constraints,
                index_offset=1,
                max_iter=max_iter,
            )
            return opt_results['final_molecule'], basis

        _, scf_results, _ = self._compute_energy(drivers[0], molecule, basis)
        opt_results = self._run_optimization(
            drivers[0],
            molecule,
            constraints=constraints,
            index_offset=1,
            compute_args=(basis, scf_results),
            max_iter=max_iter,
        )
        optimized_molecule = opt_results['final_molecule']
        basis = self._refresh_basis_for_driver(
            optimized_molecule, drivers[0], basis)

        return optimized_molecule, basis

    def _local_group_cluster_contains_rotor_kind(
            self,
            local_group_model,
            cluster,
            rotor_kind):
        rotor_kind = str(rotor_kind)
        for rotor_id in getattr(cluster, 'rotor_ids', ()):
            rotor = local_group_model.rotors.get(str(rotor_id))
            if rotor is not None and str(rotor.kind) == rotor_kind:
                return True

        return str(getattr(cluster, 'family_label', '')) == rotor_kind

    def _local_group_cluster_owned_atoms_for_rotor_kind(
            self,
            local_group_model,
            cluster,
            rotor_kind):
        rotor_kind = str(rotor_kind)
        owned_atoms = set()

        for rotor_id in getattr(cluster, 'rotor_ids', ()):
            rotor = local_group_model.rotors.get(str(rotor_id))
            if rotor is None or str(rotor.kind) != rotor_kind:
                continue
            owned_atoms.update(int(atom) for atom in rotor.owned_atoms)

        if not owned_atoms and str(getattr(cluster, 'family_label', '')) == rotor_kind:
            owned_atoms.update(int(atom) for atom in cluster.owned_atoms)

        return tuple(sorted(owned_atoms))

    def _alkyl_local_state_relaxation_constraints(
            self,
            z_matrix,
            local_group_model,
            cluster,
            dihedrals_to_rotate):
        alkyl_atoms = set(
            self._local_group_cluster_owned_atoms_for_rotor_kind(
                local_group_model, cluster, 'alkyl_chain')
        )
        if not alkyl_atoms:
            return []

        constraints = []
        seen = set()

        def add_constraint(coordinate):
            coordinate = tuple(int(atom) for atom in coordinate)
            if coordinate not in seen:
                seen.add(coordinate)
                constraints.append(coordinate)

        all_atoms = set()
        for coordinate in InterpolationDatapoint.flatten_z_matrix(z_matrix):
            all_atoms.update(int(atom) for atom in coordinate)

        frozen_atoms = sorted(all_atoms - alkyl_atoms)
        if frozen_atoms:
            frozen_atoms_one_based = ','.join(
                str(int(atom) + 1) for atom in frozen_atoms)
            constraints.append(f'freeze xyz {frozen_atoms_one_based}')

        for rotor_id in getattr(cluster, 'rotor_ids', ()):
            rotor = local_group_model.rotors.get(str(rotor_id))
            if rotor is None or str(rotor.kind) != 'alkyl_chain':
                continue
            dihedral = self._oriented_local_rotor_dihedral(rotor)
            if dihedral is not None:
                add_constraint(dihedral)

        for dihedral in dihedrals_to_rotate or ():
            add_constraint(dihedral)

        return constraints

    def _relax_alkyl_local_group_state(
            self,
            drivers,
            molecule,
            basis,
            z_matrix,
            local_group_model,
            cluster,
            dihedrals_to_rotate):
        if not getattr(self, 'relax_alkyl_local_group_states', False):
            return molecule, basis

        if not self._local_group_cluster_contains_rotor_kind(
                local_group_model, cluster, 'alkyl_chain'):
            return molecule, basis

        constraints = self._alkyl_local_state_relaxation_constraints(
            z_matrix,
            local_group_model,
            cluster,
            dihedrals_to_rotate,
        )
        if not constraints:
            return molecule, basis

        if self._is_xtb_like_driver(drivers[0]):
            opt_results = self._run_optimization(
                drivers[0],
                molecule,
                constraints=constraints,
                index_offset=1,
                max_iter=getattr(
                    self, "anchored_linker_state_opt_max_iter", 60),
            )
            optimized_molecule = opt_results['final_molecule']
            return optimized_molecule, basis

        _, scf_results, _ = self._compute_energy(drivers[0], molecule, basis)
        opt_results = self._run_optimization(
            drivers[0],
            molecule,
            constraints=constraints,
            index_offset=1,
            compute_args=(basis, scf_results),
            max_iter=getattr(
                self, "anchored_linker_state_opt_max_iter", 60),
        )
        optimized_molecule = opt_results['final_molecule']

        basis = self._refresh_basis_for_driver(
            optimized_molecule, drivers[0], basis)

        return optimized_molecule, basis

    def _methoxy_local_state_relaxation_constraints(
            self,
            molecule,
            local_group_model,
            cluster,
            rotor_ids):
        if not self._local_group_cluster_contains_rotor_kind(
                local_group_model, cluster, 'methoxy'):
            return []

        movable_atoms = set(
            int(atom) for atom in getattr(cluster, 'active_atoms', ()))
        methoxy_atoms = set(
            self._local_group_cluster_owned_atoms_for_rotor_kind(
                local_group_model, cluster, 'methoxy')
        )
        movable_atoms.update(methoxy_atoms)

        for rotor_id in getattr(cluster, 'rotor_ids', ()):
            rotor = local_group_model.rotors.get(str(rotor_id))
            if rotor is None or str(rotor.kind) != 'methoxy':
                continue
            movable_atoms.update(int(atom) for atom in rotor.axis)
            if rotor.phase_coordinate is not None:
                movable_atoms.update(
                    int(atom) for atom in rotor.phase_coordinate)

        if not movable_atoms:
            return []

        constraints = []
        all_atoms = set(range(len(molecule.get_labels())))
        frozen_atoms = sorted(all_atoms - movable_atoms)
        if frozen_atoms:
            frozen_atoms_one_based = ','.join(
                str(int(atom) + 1) for atom in frozen_atoms)
            constraints.append(f'freeze xyz {frozen_atoms_one_based}')

        seen_dihedrals = set()
        for rotor_id in rotor_ids:
            rotor = local_group_model.rotors.get(str(rotor_id))
            if rotor is None:
                continue

            dihedral = self._oriented_local_rotor_dihedral(rotor)
            if dihedral is None:
                continue

            key = tuple(int(atom) for atom in dihedral)
            canonical_key = min(key, tuple(reversed(key)))
            if canonical_key in seen_dihedrals:
                continue
            seen_dihedrals.add(canonical_key)

            constraints.append(
                "freeze dihedral "
                + " ".join(str(int(atom) + 1) for atom in key)
            )

        return constraints

    def _relax_methoxy_local_group_state(
            self,
            drivers,
            molecule,
            basis,
            local_group_model,
            cluster,
            rotor_ids):
        if not getattr(self, 'relax_methoxy_local_group_states', False):
            return molecule, basis

        constraints = self._methoxy_local_state_relaxation_constraints(
            molecule,
            local_group_model,
            cluster,
            rotor_ids,
        )
        if not constraints:
            return molecule, basis

        max_iter = getattr(self, "methoxy_state_opt_max_iter", 120)
        if self._is_xtb_like_driver(drivers[0]):
            opt_results = self._run_optimization(
                drivers[0],
                molecule,
                constraints=constraints,
                index_offset=1,
                max_iter=max_iter,
            )
            optimized_molecule = opt_results['final_molecule']
            return optimized_molecule, basis

        _, scf_results, _ = self._compute_energy(drivers[0], molecule, basis)
        opt_results = self._run_optimization(
            drivers[0],
            molecule,
            constraints=constraints,
            index_offset=1,
            compute_args=(basis, scf_results),
            max_iter=max_iter,
        )
        optimized_molecule = opt_results['final_molecule']

        basis = self._refresh_basis_for_driver(
            optimized_molecule, drivers[0], basis)

        return optimized_molecule, basis

    def _anchored_linker_local_group_state_constraints(
            self,
            molecule,
            z_matrix,
            local_group_model,
            cluster,
            rotor_ids,
            phase_signature):
        owned_atoms = set(int(atom) for atom in cluster.owned_atoms)
        if not owned_atoms:
            return []
        relaxation_atoms = set()
        if getattr(self, "anchored_linker_relax_amide", True):
            relaxation_atoms.update(
                int(atom)
                for atom in getattr(cluster, "relaxation_atoms", ())
            )
        movable_atoms = owned_atoms | relaxation_atoms

        constraints = []
        all_atoms = set(range(len(molecule.get_labels())))
        frozen_atoms = sorted(all_atoms - movable_atoms)
        if frozen_atoms:
            frozen_atoms_one_based = ','.join(
                str(int(atom) + 1) for atom in frozen_atoms)
            constraints.append(f'freeze xyz {frozen_atoms_one_based}')

        if getattr(self, "anchored_linker_freeze_boundary_angles", False):
            def element(atom):
                letters = "".join(
                    x for x in molecule.get_labels()[int(atom)]
                    if x.isalpha())
                return letters.capitalize()

            flat_z_matrix = InterpolationDatapoint.flatten_z_matrix(z_matrix)
            anchor_atoms = set(
                int(atom) for atom in getattr(cluster, "anchor_atoms", ()))
            frozen_angle_constraints = set()
            for rotor_id in rotor_ids:
                rotor = local_group_model.rotors.get(str(rotor_id))
                if rotor is None or str(rotor.kind) != 'anchored_linker':
                    continue
                for row_idx, row_type in zip(
                        getattr(rotor, "signature_rows", ()),
                        getattr(rotor, "signature_row_types", ())):
                    if str(row_type).strip().lower() != "angle":
                        continue
                    row_idx = int(row_idx)
                    if row_idx < 0 or row_idx >= len(flat_z_matrix):
                        continue
                    coordinate = tuple(
                        int(atom) for atom in flat_z_matrix[row_idx])
                    if len(coordinate) != 3:
                        continue
                    coord_set = set(coordinate)
                    local_in_angle = coord_set & owned_atoms
                    anchor_in_angle = coord_set & anchor_atoms
                    if len(local_in_angle) != 2 or len(anchor_in_angle) != 1:
                        continue
                    anchor_atom = next(iter(anchor_in_angle))
                    if element(anchor_atom) != "C":
                        continue
                    key = (
                        coordinate
                        if coordinate <= tuple(reversed(coordinate))
                        else tuple(reversed(coordinate))
                    )
                    if key in frozen_angle_constraints:
                        continue
                    frozen_angle_constraints.add(key)
                    constraints.append(
                        "freeze angle "
                        + " ".join(str(int(atom) + 1) for atom in coordinate)
                    )

        for rotor_id, phase in zip(rotor_ids, phase_signature):
            rotor = local_group_model.rotors.get(str(rotor_id))
            if rotor is None:
                continue

            dihedral = self._oriented_local_rotor_dihedral(rotor)
            if dihedral is None:
                continue

            dihedral_one_based = [int(atom) + 1 for atom in dihedral]
            if str(rotor.kind) != 'anchored_linker':
                if abs(float(phase)) > 1.0e-14:
                    constraints.append(
                        "freeze dihedral "
                        + " ".join(str(atom) for atom in dihedral_one_based)
                    )
                continue

            current_angle = molecule.get_dihedral(
                dihedral_one_based, 'degree')
            target_angle = current_angle + np.rad2deg(float(phase))
            constraints.append(
                "set dihedral "
                + " ".join(str(atom) for atom in dihedral_one_based)
                + f" {target_angle:.12f}"
            )

        return constraints

    def _relax_anchored_linker_local_group_state(
            self,
            drivers,
            molecule,
            basis,
            z_matrix,
            local_group_model,
            cluster,
            rotor_ids,
            phase_signature):
        if not self._local_group_cluster_contains_rotor_kind(
                local_group_model, cluster, 'anchored_linker'):
            return molecule, basis

        constraints = self._anchored_linker_local_group_state_constraints(
            molecule,
            z_matrix,
            local_group_model,
            cluster,
            rotor_ids,
            phase_signature,
        )
        if not constraints:
            return molecule, basis

        if self._is_xtb_like_driver(drivers[0]):
            opt_results = self._run_optimization(
                drivers[0],
                molecule,
                constraints=constraints,
                index_offset=1,
                max_iter=getattr(
                    self, "anchored_linker_state_opt_max_iter", 60),
            )
            optimized_molecule = opt_results['final_molecule']
            return optimized_molecule, basis

        _, scf_results, _ = self._compute_energy(drivers[0], molecule, basis)
        opt_results = self._run_optimization(
            drivers[0],
            molecule,
            constraints=constraints,
            index_offset=1,
            compute_args=(basis, scf_results),
            max_iter=getattr(
                self, "anchored_linker_state_opt_max_iter", 60),
        )
        optimized_molecule = opt_results['final_molecule']

        basis = self._refresh_basis_for_driver(
            optimized_molecule, drivers[0], basis)

        return optimized_molecule, basis

    def _map_local_group_unit_atoms(self, labels, source_atoms, target_atoms):
        source_atoms = tuple(int(x) for x in source_atoms)
        target_atoms = tuple(int(x) for x in target_atoms)

        source_by_element = {}
        target_by_element = {}
        for atom in source_atoms:
            source_by_element.setdefault(str(labels[atom]), []).append(atom)
        for atom in target_atoms:
            target_by_element.setdefault(str(labels[atom]), []).append(atom)

        mapping = {}
        if set(source_by_element) != set(target_by_element):
            for source_atom, target_atom in zip(sorted(source_atoms), sorted(target_atoms)):
                mapping[int(source_atom)] = int(target_atom)
            return mapping

        for element in source_by_element:
            source_group = sorted(source_by_element[element])
            target_group = sorted(target_by_element[element])
            if len(source_group) != len(target_group):
                for source_atom, target_atom in zip(sorted(source_atoms), sorted(target_atoms)):
                    mapping[int(source_atom)] = int(target_atom)
                return mapping
            for source_atom, target_atom in zip(source_group, target_group):
                mapping[int(source_atom)] = int(target_atom)

        return mapping

    def _local_rotor_atom_mapping(self, molecule, rotor, shift):
        symmetry_order = max(int(rotor.symmetry_order), 1)
        shift = int(shift) % symmetry_order
        if shift == 0:
            return {}

        labels = molecule.get_labels()
        units = tuple(tuple(int(atom) for atom in unit) for unit in rotor.unit_atom_sets)
        if len(units) != symmetry_order:
            return {}

        mapping = {}
        for unit_idx, source_atoms in enumerate(units):
            target_atoms = units[(unit_idx + shift) % symmetry_order]
            mapping.update(
                self._map_local_group_unit_atoms(labels, source_atoms, target_atoms)
            )

        return mapping

    def _compose_local_atom_mapping(self, natoms, left_mapping, right_mapping):
        composed = {}
        for atom in range(natoms):
            mid = right_mapping.get(atom, atom)
            target = left_mapping.get(mid, mid)
            if target != atom:
                composed[atom] = target
        return composed

    def _build_local_group_mapping_masks(self, molecule, z_matrix, local_group_model, rotor_ids):
        flat_z_matrix = InterpolationDatapoint.flatten_z_matrix(z_matrix)
        ident = list(range(len(flat_z_matrix)))
        if not rotor_ids:
            return [ident]

        z_lookup = {}
        for row_idx, row in enumerate(flat_z_matrix):
            row_key = tuple(int(x) for x in row)
            z_lookup[row_key] = int(row_idx)
            z_lookup[tuple(reversed(row_key))] = int(row_idx)

        rotors = [local_group_model.rotors[str(rotor_id)] for rotor_id in rotor_ids]
        shift_options = [
            tuple(range(max(int(rotor.symmetry_order), 1)))
            for rotor in rotors
        ]

        natoms = len(molecule.get_labels())
        masks = [ident]
        seen = {tuple(ident)}

        for shifts in itertools.product(*shift_options):
            if all(int(shift) == 0 for shift in shifts):
                continue

            atom_mapping = {}
            # Apply child/local permutations before parent/subtree
            # permutations. This matches the atom relabeling induced by setting
            # child rotor phases on the already parent-rotated geometry.
            for rotor, shift in reversed(list(zip(rotors, shifts))):
                rotor_mapping = self._local_rotor_atom_mapping(
                    molecule, rotor, int(shift))
                atom_mapping = self._compose_local_atom_mapping(
                    natoms, rotor_mapping, atom_mapping)

            mask = []
            valid = True
            for row in flat_z_matrix:
                mapped = tuple(atom_mapping.get(int(atom), int(atom)) for atom in row)
                mapped_idx = z_lookup.get(mapped)
                if mapped_idx is None:
                    valid = False
                    break
                mask.append(mapped_idx)

            if not valid:
                continue

            key = tuple(mask)
            if key not in seen:
                seen.add(key)
                masks.append(mask)

        return masks

    def _set_local_group_datapoint_metadata(
            self,
            datapoint,
            *,
            family_label,
            bank_role,
            local_group_model,
            cluster=None,
            cluster_state_id=0,
            rotor_ids=(),
            dihedrals_to_rotate=None,
            phase_signature=None,
            reference_molecule=None,
            is_anchor=False):

        datapoint.family_label = family_label
        datapoint.bank_role = bank_role
        datapoint.cluster_state_id = int(cluster_state_id)
        datapoint.cluster_rotor_ids = tuple(str(rotor_id) for rotor_id in rotor_ids)
        datapoint.dihedrals_to_rotate = dihedrals_to_rotate
        datapoint.phase_signature = phase_signature
        datapoint.is_anchor = bool(is_anchor)
        datapoint.local_factor_combination_mode = getattr(
            self,
            "local_factor_combination_mode",
            COMBINATION_SIGNED_FULL,
        )

        if bank_role == 'core':
            datapoint.cluster_id = None
            datapoint.cluster_type = None
            datapoint.active_atoms = np.array(local_group_model.core_atoms, dtype=np.int64)
            datapoint.active_rows = np.array(local_group_model.core_rows, dtype=np.int64)
            datapoint.canonical_subset_key = ()
            datapoint.response_rows = ()
            datapoint.relaxation_policy_id = None
            datapoint.anchor_state_ids = ()
            datapoint.local_factor_overlap_source = None
        else:
            datapoint.cluster_id = str(cluster.cluster_id)
            datapoint.cluster_type = str(cluster.family_label)
            datapoint.active_atoms = np.array(cluster.active_atoms, dtype=np.int64)
            datapoint.active_rows = np.array(cluster.active_rows, dtype=np.int64)
            datapoint.canonical_subset_key = tuple(
                str(rotor_id)
                for rotor_id in (
                    getattr(cluster, "canonical_subset_key", ())
                    or getattr(cluster, "rotor_ids", ())
                )
            )
            datapoint.response_rows = tuple(
                int(row) for row in getattr(cluster, "response_rows", ()))
            datapoint.relaxation_policy_id = getattr(
                cluster, "relaxation_policy_id", "default")
            datapoint.anchor_state_ids = tuple(
                int(state_id)
                for state_id in getattr(cluster, "anchor_state_ids", ()))
            datapoint.local_factor_overlap_source = (
                ",".join(
                    str(parent_id)
                    for parent_id in getattr(
                        cluster, "parent_cluster_ids", ())
                )
                if str(getattr(cluster, "role", "factor")).lower() == "overlap"
                else None
            )

        if reference_molecule is not None:
            mapping_rotor_ids = tuple(str(rotor_id) for rotor_id in rotor_ids)
            if bank_role == 'core':
                # The schema-3 core is evaluated with coefficient one and must
                # retain the full local-rotor permutation symmetry.  Passing an
                # empty rotor list here leaves only the identity mask, causing
                # a symmetry-equivalent methyl rotation to look like a large
                # displacement of the full core Taylor surface.
                mapping_rotor_ids = tuple(
                    str(rotor_id)
                    for rotor_id in sorted(
                        local_group_model.rotors, key=lambda value: str(value))
                )
            datapoint.mapping_masks = np.array(
                self._build_local_group_mapping_masks(
                    reference_molecule,
                    datapoint.z_matrix_dict,
                    local_group_model,
                    mapping_rotor_ids,
                ),
                dtype=np.int64,
            )
            if (
                datapoint.use_eq_bond_length
                and str(datapoint.eq_bond_symmetry_mode).strip().lower()
                == 'symmetrized'
            ):
                datapoint.symmetrize_eq_bond_lengths_from_masks(
                    datapoint.mapping_masks)
                datapoint.transform_gradient_and_hessian()

    def _create_local_group_datapoint(
            self,
            molecule,
            z_matrix,
            interpolation_settings,
            energy,
            gradient,
            hessian,
            inv_sqrt_masses,
            eq_bond_lengths,
            imp_int_constraints):

        grad_vec = gradient.reshape(-1)
        hess_mat = hessian.reshape(grad_vec.size, grad_vec.size)

        mw_grad_vec = grad_vec
        mw_hess_mat = hess_mat
        if self.use_mass_weight:
            mw_grad_vec = inv_sqrt_masses * grad_vec
            mw_hess_mat = (
                inv_sqrt_masses[:, None] * hess_mat
            ) * inv_sqrt_masses[None, :]

        datapoint = InterpolationDatapoint(z_matrix)
        datapoint.update_settings(interpolation_settings)
        datapoint.cartesian_coordinates = molecule.get_coordinates_in_bohr()
        datapoint.eq_bond_lengths = eq_bond_lengths
        datapoint.imp_int_coordinates = imp_int_constraints
        datapoint.inv_sqrt_masses = inv_sqrt_masses
        datapoint.energy = energy
        datapoint.gradient = mw_grad_vec.reshape(gradient.shape)
        datapoint.hessian = mw_hess_mat.reshape(hessian.shape)
        datapoint.transform_gradient_and_hessian()
        datapoint.confidence_radius = self.use_opt_confidence_radius[2]

        return datapoint

    def _next_local_group_family_label(self, target_file, interpolation_driver):
        labels = []
        if Path(target_file).exists():
            org_labels, _ = interpolation_driver.read_labels()
            labels = list(org_labels)

        max_index = 0
        for label in labels:
            parts = str(label).split('_')
            if len(parts) < 2 or parts[0] != 'point':
                continue
            try:
                max_index = max(max_index, int(parts[1]))
            except ValueError:
                continue

        return f'point_{max_index + 1}'

    def _local_group_state_jobs(self, local_group_model):
        for cluster in local_group_model.clusters.values():
            rotor_ids = tuple(str(rotor_id) for rotor_id in cluster.rotor_ids)

            yield {
                'cluster': cluster,
                'state_id': 0,
                'rotor_ids': rotor_ids,
                'phase_signature': tuple(0.0 for _ in rotor_ids),
                'is_anchor': True,
            }

            phase_lists = [
                self._local_group_phase_values(local_group_model.rotors[rotor_id])
                for rotor_id in rotor_ids
            ]

            state_id = 1
            for phase_signature in itertools.product(*phase_lists):
                if all(abs(float(phase)) < 1.0e-14 for phase in phase_signature):
                    continue

                yield {
                    'cluster': cluster,
                    'state_id': state_id,
                    'rotor_ids': rotor_ids,
                    'phase_signature': tuple(float(x) for x in phase_signature),
                    'is_anchor': False,
                }
                state_id += 1

    def _detect_local_factor_response_rows(
            self, core_datapoint, cluster_datapoints, z_matrix,
            force_detection=False):
        """Detect internal-coordinate response to constrained relaxation."""

        if not cluster_datapoints or (
            not force_detection
            and not bool(getattr(
                self, "local_factor_detect_response_rows", True))
        ):
            return ()

        thresholds = dict(getattr(
            self,
            "local_factor_response_thresholds",
            {},
        ))
        defaults = {
            "bond": 0.01,
            "angle": np.deg2rad(1.0),
            "dihedral": np.deg2rad(3.0),
            "improper": np.deg2rad(3.0),
        }
        defaults.update(thresholds)

        sections = (
            ("bond", len(z_matrix.get("bonds", ()))),
            ("angle", len(z_matrix.get("angles", ()))),
            ("dihedral", len(z_matrix.get("dihedrals", ()))),
            ("improper", len(z_matrix.get("impropers", ()))),
        )
        row_thresholds = []
        periodic_rows = []
        offset = 0
        for section, count in sections:
            row_thresholds.extend([float(defaults[section])] * int(count))
            if section in ("dihedral", "improper"):
                periodic_rows.extend(range(offset, offset + int(count)))
            offset += int(count)

        core_values = np.asarray(
            core_datapoint.internal_coordinates_values, dtype=np.float64)
        if len(row_thresholds) != core_values.size:
            raise RuntimeError(
                "Local-factor response detection found an inconsistent "
                "z-matrix/internal-coordinate size."
            )

        displacements = []
        for datapoint in cluster_datapoints:
            delta = np.asarray(
                datapoint.internal_coordinates_values,
                dtype=np.float64,
            ) - core_values
            if periodic_rows:
                idx = np.asarray(periodic_rows, dtype=np.int64)
                delta[idx] = (delta[idx] + np.pi) % (2.0 * np.pi) - np.pi
            displacements.append(delta)

        rms = np.sqrt(np.mean(np.square(displacements), axis=0))
        return tuple(int(row) for row in np.flatnonzero(
            rms > np.asarray(row_thresholds, dtype=np.float64)))

    def _promote_local_factor_response_rows(
            self, local_group_model, cluster_datapoints, core_datapoint,
            z_matrix, force_detection=False):
        """Return a model whose factor support includes relaxed responses."""

        flat_z_matrix = InterpolationDatapoint.flatten_z_matrix(z_matrix)
        clusters = {}
        for cluster_id, cluster in local_group_model.clusters.items():
            detected = self._detect_local_factor_response_rows(
                core_datapoint,
                cluster_datapoints.get(str(cluster_id), ()),
                z_matrix,
                force_detection=force_detection,
            )
            active_rows = set(int(row) for row in cluster.active_rows)
            response_rows = set(int(row) for row in cluster.response_rows)
            response_rows.update(detected)
            active_rows.update(response_rows)
            for rotor_id in cluster.rotor_ids:
                rotor = local_group_model.rotors.get(str(rotor_id))
                if rotor is not None:
                    active_rows.update(int(row) for row in rotor.torsion_rows)

            active_atoms = set(int(atom) for atom in cluster.active_atoms)
            for row in active_rows:
                active_atoms.update(int(atom) for atom in flat_z_matrix[row])
            policy_id = getattr(cluster, "relaxation_policy_id", None)
            if not policy_id or str(policy_id) == "default":
                kinds = {
                    str(local_group_model.rotors[str(rotor_id)].kind)
                    for rotor_id in cluster.rotor_ids
                    if str(rotor_id) in local_group_model.rotors
                }
                if "anchored_linker" in kinds:
                    policy_id = "anchored_linker_torsion_relaxed"
                elif "alkyl_chain" in kinds:
                    policy_id = "alkyl_chain_torsion_relaxed"
                elif "methoxy" in kinds:
                    policy_id = "methoxy_torsion_relaxed"
                elif kinds & {"alcohol", "primary_amine", "ammonium"}:
                    policy_id = "heteroatom_torsion_relaxed"
                elif "tert_butyl_parent" in kinds:
                    policy_id = "tertbutyl_parent_torsion_relaxed"
                else:
                    policy_id = "methyl_torsion_relaxed"

            clusters[str(cluster_id)] = replace(
                cluster,
                active_rows=tuple(sorted(active_rows)),
                active_atoms=tuple(sorted(active_atoms)),
                response_rows=tuple(sorted(response_rows)),
                relaxation_policy_id=str(policy_id),
            )

        factor_rows = {
            int(row)
            for cluster in clusters.values()
            for row in cluster.active_rows
        }
        core_rows = tuple(
            int(row) for row in local_group_model.core_rows
            if int(row) not in factor_rows
        )
        return replace(
            local_group_model,
            clusters=clusters,
            core_rows=core_rows,
        )

    def add_point_local_groups(
            self,
            molecule_specific_information,
            interpolation_settings,
            symmetry_information=None,
            write_org_database=True):
        """
        Adds local-group database points.

        For each incoming molecular structure this writes:
            - one core datapoint: point_N_core
            - one anchor per local cluster:
                point_N_cluster_<cluster_id>_state_0_anchor
            - one state datapoint per sampled local phase combination:
                point_N_cluster_<cluster_id>_state_M

        The core point uses core atom/row masks. Each cluster point uses the
        detected cluster atom/row masks and local symmetry permutation masks.
        """

        if symmetry_information is None:
            symmetry_information = {}

        self.ostream.print_blank()
        self.ostream.print_header(
            "Constructing local-group decoupled interpolation database.")
        self.ostream.flush()

        if len(self.drivers) == 0:
            raise ValueError("No energy driver defined.")

        adjusted_molecule = {'gs': [], 'es': []}

        for entries in molecule_specific_information:
            molecule = entries[0]
            symmetry_point = False

            if self.eq_bond_length is None:
                self.eq_bond_length = []
                for element in self.roots_z_matrix[0]['bonds']:
                    if (
                        len(element) == 2
                        and self.use_minimized_structures[0]
                        and self.eq_bond_length_irc_bonds is not None
                        and element not in self.eq_bond_length_irc_bonds
                    ):
                        self.eq_bond_length.append(
                            molecule.get_distance(
                                [element[0] + 1, element[1] + 1], 'bohr'))
                    elif (
                        len(element) == 2
                        and self.use_minimized_structures[0]
                        and self.eq_bond_length_irc_bonds is not None
                        and element in self.eq_bond_length_irc_bonds
                    ):
                        self.eq_bond_length.append(
                            molecule_specific_information[-1][0].get_distance(
                                [element[0] + 1, element[1] + 1], 'bohr'))
                    elif len(element) == 2 and self.use_minimized_structures[0]:
                        self.eq_bond_length.append(
                            molecule.get_distance(
                                [element[0] + 1, element[1] + 1], 'bohr'))
                    elif len(element) == 2:
                        self.eq_bond_length.append(0.0)

            if 0 in entries[2]:
                adjusted_molecule['gs'].append(
                    (entries[0], entries[1], 1, None, [0], symmetry_point, entries[3]))
            if any(x > 0 for x in entries[2]):
                states = [state for state in entries[2] if state > 0]
                adjusted_molecule['es'].append(
                    (entries[0], entries[1], 1, None, states, symmetry_point, entries[3]))

        for key, entries in adjusted_molecule.items():
            if len(entries) == 0:
                continue

            drivers = self.drivers[key]
            for mol_basis in entries:
                imp_int_constraints = {
                    'bonds': [],
                    'angles': [],
                    'dihedrals': [],
                    'impropers': [],
                }
                incoming_constraints = mol_basis[-1]
                if isinstance(incoming_constraints, dict):
                    for coordinate_type in imp_int_constraints:
                        imp_int_constraints[coordinate_type].extend(
                            tuple(int(x) for x in coordinate)
                            for coordinate in incoming_constraints.get(
                                coordinate_type, ())
                        )
                else:
                    for constraint in incoming_constraints:
                        constraint = tuple(int(x) for x in constraint)
                        if len(constraint) == 2:
                            imp_int_constraints["bonds"].append(constraint)
                        elif len(constraint) == 3:
                            imp_int_constraints["angles"].append(constraint)
                        elif len(constraint) == 4:
                            imp_int_constraints["dihedrals"].append(constraint)

                energies, scf_results, rsp_results = self._compute_energy(
                    drivers[0], mol_basis[0], mol_basis[1])

                if isinstance(drivers[0], LinearResponseEigenSolver) or isinstance(drivers[0], TdaEigenSolver):
                    energies = energies[mol_basis[4]]

                gradients = self._compute_gradient(
                    drivers[1], mol_basis[0], mol_basis[1], scf_results, rsp_results)
                hessians = self._compute_hessian(
                    drivers[2], mol_basis[0], mol_basis[1])
                backend_metadata = self._collect_backend_metadata(drivers)

                inv_sqrt_masses = None
                if self.use_mass_weight:
                    masses = mol_basis[0].get_masses().copy()
                    masses_cart = np.repeat(masses, 3)
                    inv_sqrt_masses = 1.0 / np.sqrt(masses_cart)

                for number in range(len(energies)):
                    target_root = mol_basis[4][number]
                    target_file = interpolation_settings[target_root]['imforcefield_file']
                    z_matrix = self.roots_z_matrix[target_root]
                    local_group_model = self.local_group_primitive_model.get(
                        target_root,
                        self.local_group_model.get(target_root),
                    )

                    if local_group_model is None or not local_group_model.enabled:
                        raise RuntimeError(
                            "Local-group database construction was requested, "
                            f"but no local groups were detected for root {target_root}."
                        )

                    interpolation_driver = InterpolationDriver(z_matrix)
                    interpolation_driver.update_settings(
                        interpolation_settings[target_root])
                    interpolation_driver.imforcefield_file = target_file

                    family_label = self._next_local_group_family_label(
                        target_file, interpolation_driver)
                    core_label = f'{family_label}_core'

                    core_dp = self._create_local_group_datapoint(
                        molecule=mol_basis[0],
                        z_matrix=z_matrix,
                        interpolation_settings=interpolation_settings[target_root],
                        energy=energies[number],
                        gradient=gradients[number].copy(),
                        hessian=hessians[number].copy(),
                        inv_sqrt_masses=inv_sqrt_masses,
                        eq_bond_lengths=self.eq_bond_length,
                        imp_int_constraints=imp_int_constraints,
                    )
                    core_dp.metadata = deepcopy(backend_metadata)
                    core_dp.point_label = core_label

                    schema4_ordered_construction = self._schema4_msi_enabled()
                    local_group_model, coupling_data = self._build_coupled_local_group_model(
                        local_group_model,
                        core_dp.internal_hessian,
                    )
                    if not schema4_ordered_construction:
                        local_group_model, added_overlap_clusters = (
                            self._augment_schema3_relaxed_residual_clusters(
                                local_group_model)
                        )
                        local_group_model = (
                            self._with_schema3_anchor_state_ids(local_group_model)
                        )
                        local_group_model = self._with_schema3_projector_rows(
                            local_group_model)
                        if added_overlap_clusters:
                            coupling_data = dict(coupling_data)
                            coupling_data[
                                'schema3_canonical_overlap_cluster_ids'
                            ] = added_overlap_clusters
                            coupling_data['coupled_cluster_ids'] = tuple(
                                local_group_model.clusters.keys())
                    else:
                        coupling_data = dict(coupling_data)
                        coupling_data["construction_status"] = (
                            "provisional_topology_before_response_audit"
                        )
                    self.local_group_coupled_model[target_root] = local_group_model
                    self.local_group_model[target_root] = local_group_model
                    self.local_groups[target_root] = local_group_model.clusters
                    self.local_group_coupling_matrix[target_root] = coupling_data
                    self.local_group_coupling_matrix[
                        (target_root, family_label)
                    ] = coupling_data

                    print(coupling_data)

                    point_labels = [core_label]
                    point_labels_by_cluster = {}
                    generated_cluster_records = []
                    anchored_linker_state_cache = {}
                    local_state_continuation_cache = {}
                    anchored_linker_failed_states = set()

                    for job in self._local_group_state_jobs(local_group_model):
                        cluster = job['cluster']
                        state_id = int(job['state_id'])
                        rotor_ids = tuple(job['rotor_ids'])
                        phase_signature = tuple(job['phase_signature'])
                        label = (
                            f"{family_label}_cluster_"
                            f"{cluster.cluster_id}_state_{state_id}"
                        )
                        if job['is_anchor']:
                            label += "_anchor"

                        dihedrals_to_rotate = []
                        for rotor_id, phase in zip(rotor_ids, phase_signature):
                            if abs(float(phase)) < 1.0e-14:
                                continue
                            dihedral = self._oriented_local_rotor_dihedral(
                                local_group_model.rotors[rotor_id])
                            if dihedral is not None:
                                dihedrals_to_rotate.append(dihedral)
                        if job['is_anchor']:
                            cluster_dp: InterpolationDatapoint = self._create_local_group_datapoint(
                                molecule=mol_basis[0],
                                z_matrix=z_matrix,
                                interpolation_settings=interpolation_settings[target_root],
                                energy=energies[number],
                                gradient=gradients[number].copy(),
                                hessian=hessians[number].copy(),
                                inv_sqrt_masses=inv_sqrt_masses,
                                eq_bond_lengths=self.eq_bond_length,
                                imp_int_constraints=imp_int_constraints,
                            )
                            cluster_dp.metadata = deepcopy(backend_metadata)
                            cluster_molecule = mol_basis[0]
                        else:
                            has_anchored_linker = (
                                self._local_group_cluster_contains_rotor_kind(
                                    local_group_model,
                                    cluster,
                                    'anchored_linker',
                                )
                            )
                            use_canonical_relaxation = (
                                self._schema3_relaxed_residual_enabled()
                            )
                            cluster_molecule = self._apply_local_group_state(
                                mol_basis[0],
                                local_group_model,
                                rotor_ids,
                                phase_signature,
                                skip_kinds=(
                                    {'anchored_linker'}
                                    if (
                                        has_anchored_linker
                                        and not use_canonical_relaxation
                                    ) else None
                                ),
                            )

                            cluster_basis = mol_basis[1]
                            if (
                                not self._is_xtb_like_driver(drivers[0])
                                and cluster_basis is not None
                                and hasattr(cluster_basis, 'get_main_basis_label')
                            ):
                                cluster_basis = MolecularBasis.read(
                                    cluster_molecule,
                                    cluster_basis.get_main_basis_label(),
                                )

                            if use_canonical_relaxation:
                                cluster_molecule, cluster_basis = (
                                    self._relax_canonical_local_factor_state(
                                        drivers,
                                        cluster_molecule,
                                        cluster_basis,
                                        z_matrix,
                                        local_group_model,
                                        cluster,
                                        rotor_ids,
                                        phase_signature,
                                        anchor_molecule=mol_basis[0],
                                    )
                                )
                            elif has_anchored_linker:
                                anchored_phase_signature = tuple(
                                    float(phase)
                                    if (
                                        local_group_model.rotors[
                                            str(rotor_id)].kind
                                        == 'anchored_linker'
                                    ) else 0.0
                                    for rotor_id, phase in zip(
                                        rotor_ids, phase_signature)
                                )
                                anchored_cache_key = (
                                    str(cluster.cluster_id),
                                    tuple(
                                        round(float(phase), 12)
                                        for phase in anchored_phase_signature
                                    ),
                                )
                                if anchored_cache_key in anchored_linker_failed_states:
                                    continue

                                if anchored_cache_key in anchored_linker_state_cache:
                                    cached_molecule, cached_basis = (
                                        anchored_linker_state_cache[
                                            anchored_cache_key]
                                    )
                                    cluster_molecule = (
                                        self._clone_molecule_for_local_group(
                                            cached_molecule)
                                    )
                                    cluster_basis = cached_basis
                                else:
                                    target_signature = anchored_cache_key[1]
                                    predecessor_key = None
                                    predecessor_score = -1
                                    target_nonzero = sum(
                                        abs(float(phase)) > 1.0e-12
                                        for phase in target_signature
                                    )
                                    if target_nonzero > 1:
                                        for cache_key in anchored_linker_state_cache:
                                            if cache_key[0] != anchored_cache_key[0]:
                                                continue
                                            seed_signature = cache_key[1]
                                            if len(seed_signature) != len(
                                                    target_signature):
                                                continue
                                            score = 0
                                            compatible = True
                                            for seed_phase, target_phase in zip(
                                                    seed_signature,
                                                    target_signature):
                                                seed_phase = float(seed_phase)
                                                target_phase = float(
                                                    target_phase)
                                                if abs(target_phase) <= 1.0e-12:
                                                    if abs(seed_phase) > 1.0e-12:
                                                        compatible = False
                                                        break
                                                    continue
                                                if abs(seed_phase - target_phase) <= 1.0e-12:
                                                    score += 1
                                                elif abs(seed_phase) > 1.0e-12:
                                                    compatible = False
                                                    break
                                            if (
                                                compatible
                                                and 0 < score < target_nonzero
                                                and score > predecessor_score
                                            ):
                                                predecessor_key = cache_key
                                                predecessor_score = score

                                    relaxation_phase_signature = (
                                        anchored_phase_signature
                                    )
                                    if predecessor_key is not None:
                                        cached_molecule, cached_basis = (
                                            anchored_linker_state_cache[
                                                predecessor_key]
                                        )
                                        cluster_molecule = (
                                            self._clone_molecule_for_local_group(
                                                cached_molecule)
                                        )
                                        cluster_basis = cached_basis
                                        relaxation_phase_signature = tuple(
                                            float(target) - float(seed)
                                            for target, seed in zip(
                                                target_signature,
                                                predecessor_key[1])
                                        )
                                    else:
                                        cluster_molecule = (
                                            self._clone_molecule_for_local_group(
                                                mol_basis[0])
                                        )
                                        cluster_basis = mol_basis[1]
                                    print("Relax the anchored linker")
                                    try:
                                        cluster_molecule, cluster_basis = (
                                            self._relax_anchored_linker_local_group_state(
                                                drivers,
                                                cluster_molecule,
                                                cluster_basis,
                                                z_matrix,
                                                local_group_model,
                                                cluster,
                                                rotor_ids,
                                                relaxation_phase_signature,
                                            )
                                        )
                                    except Exception as exc:
                                        anchored_linker_failed_states.add(
                                            anchored_cache_key)
                                        self.ostream.print_warning(
                                            "Skipping anchored-linker local state "
                                            f"{label}: {exc}"
                                        )
                                        self.ostream.flush()
                                        continue
                                    anchored_linker_state_cache[
                                        anchored_cache_key] = (
                                            self._clone_molecule_for_local_group(
                                                cluster_molecule),
                                            cluster_basis,
                                        )

                                continuation_key = (str(cluster.cluster_id), anchored_cache_key[1])
                                continuation = local_state_continuation_cache.get(continuation_key)

                                application_phase_signature = phase_signature

                                if continuation is not None:
                                    previous_molecule, previous_basis, previous_phase_signature = continuation
                                    cluster_molecule = self._clone_molecule_for_local_group(previous_molecule)
                                    cluster_basis = previous_basis

                                    application_phase_signature = []
                                    for rotor_id, target_phase, previous_phase in zip(rotor_ids, phase_signature, previous_phase_signature):
                                        rotor = local_group_model.rotors.get(str(rotor_id))

                                        if rotor is not None and str(rotor.kind) == 'anchored_linker':
                                            application_phase_signature.append(0.0)
                                            continue

                                        delta = float(target_phase) - float(previous_phase)
                                        delta = (delta + np.pi) % (2.0 * np.pi) - np.pi
                                        application_phase_signature.append(delta)

                                    application_phase_signature = tuple(application_phase_signature)
                                cluster_molecule = self._apply_local_group_state(
                                    cluster_molecule,
                                    local_group_model,
                                    rotor_ids,
                                    application_phase_signature,
                                    skip_kinds={'anchored_linker'},
                                )
                                cluster_molecule, cluster_basis = (
                                    self._relax_methoxy_local_group_state(
                                        drivers,
                                        cluster_molecule,
                                        cluster_basis,
                                        local_group_model,
                                        cluster,
                                        rotor_ids,
                                    )
                                )

                                local_state_continuation_cache[continuation_key] = (
                                    self._clone_molecule_for_local_group(cluster_molecule),
                                    cluster_basis,
                                    tuple(float(x) for x in phase_signature),
                                )
                                print("I am in the addition of continuation cache", continuation_key)
                            else:
                                print("Relax the methoxy group")
                                cluster_molecule, cluster_basis = (
                                    self._relax_methoxy_local_group_state(
                                        drivers,
                                        cluster_molecule,
                                        cluster_basis,
                                        local_group_model,
                                        cluster,
                                        rotor_ids,
                                    )
                                )
                                print("Relax the alkyl chain")
                                cluster_molecule, cluster_basis = (
                                    self._relax_alkyl_local_group_state(
                                        drivers,
                                        cluster_molecule,
                                        cluster_basis,
                                        z_matrix,
                                        local_group_model,
                                        cluster,
                                        dihedrals_to_rotate,
                                    )
                                )

                            print("After the relaxiation")
                            cluster_energies, cluster_scf, cluster_rsp = (
                                self._compute_energy(
                                    drivers[0], cluster_molecule, cluster_basis)
                            )
                            if isinstance(drivers[0], LinearResponseEigenSolver) or isinstance(drivers[0], TdaEigenSolver):
                                cluster_energies = cluster_energies[mol_basis[4]]

                            cluster_gradients = self._compute_gradient(
                                drivers[1],
                                cluster_molecule,
                                cluster_basis,
                                cluster_scf,
                                cluster_rsp,
                            )
                            cluster_hessians = self._compute_hessian(
                                drivers[2], cluster_molecule, cluster_basis)
                            cluster_backend_metadata = (
                                self._collect_backend_metadata(drivers))

                            cluster_dp = self._create_local_group_datapoint(
                                molecule=cluster_molecule,
                                z_matrix=z_matrix,
                                interpolation_settings=interpolation_settings[target_root],
                                energy=cluster_energies[number],
                                gradient=cluster_gradients[number].copy(),
                                hessian=cluster_hessians[number].copy(),
                                inv_sqrt_masses=inv_sqrt_masses,
                                eq_bond_lengths=self.eq_bond_length,
                                imp_int_constraints=imp_int_constraints,
                            )
                            cluster_dp.metadata = deepcopy(
                                cluster_backend_metadata)

                        cluster_dp.point_label = label
                        generated_cluster_records.append({
                            'datapoint': cluster_dp,
                            'molecule': cluster_molecule,
                            'label': label,
                            'cluster_id': str(cluster.cluster_id),
                            'state_id': state_id,
                            'rotor_ids': rotor_ids,
                            'dihedrals_to_rotate': dihedrals_to_rotate,
                            'phase_signature': phase_signature,
                            'is_anchor': bool(job['is_anchor']),
                        })
                        point_labels.append(label)
                        point_labels_by_cluster.setdefault(
                            str(cluster.cluster_id), {})[int(state_id)] = {
                                'state_id': int(state_id),
                                'label': str(label),
                                'rotor_ids': tuple(str(rotor_id) for rotor_id in rotor_ids),
                                'phase_signature': tuple(
                                    float(phase) for phase in phase_signature),
                                'is_anchor': bool(job['is_anchor']),
                            }
                        
                        self.ostream.print_blank()
                        self.ostream.print_header(
                            "Local-group database expansion: Added cluster point "
                            f"{label} for root {target_root}."
                        ) 
                        self.ostream.flush()
                        
                    if (
                        self._schema3_relaxed_residual_enabled()
                        or schema4_ordered_construction
                    ):
                        datapoints_by_cluster = {}
                        for record in generated_cluster_records:
                            datapoints_by_cluster.setdefault(
                                record['cluster_id'], []).append(
                                    record['datapoint'])
                        local_group_model = (
                            self._promote_local_factor_response_rows(
                                local_group_model,
                                datapoints_by_cluster,
                                core_dp,
                                z_matrix,
                                force_detection=schema4_ordered_construction,
                            )
                        )
                        if schema4_ordered_construction:
                            final_model, final_coupling_data = (
                                self._build_coupled_local_group_model(
                                    local_group_model,
                                    core_dp.internal_hessian,
                                )
                            )
                            provisional_ids = set(
                                str(value)
                                for value in local_group_model.clusters
                            )
                            final_ids = set(
                                str(value) for value in final_model.clusters
                            )
                            if final_ids == provisional_ids:
                                local_group_model = final_model
                            else:
                                # A newly discovered coupling requires joint QM
                                # states.  Preserve independent factors and
                                # persist the deferred edge instead of silently
                                # assembling a factor from missing states.
                                final_coupling_data = dict(final_coupling_data)
                                final_coupling_data[
                                    "deferred_coupling_requires_states"
                                ] = tuple(sorted(final_ids - provisional_ids))
                            coupling_data = dict(final_coupling_data)
                            coupling_data["construction_status"] = (
                                "finalized_after_response_audit"
                            )
                            local_group_model, added_overlap_clusters = (
                                self._augment_schema3_relaxed_residual_clusters(
                                    local_group_model
                                )
                            )
                            local_group_model = (
                                self._with_schema3_anchor_state_ids(
                                    local_group_model
                                )
                            )
                            local_group_model = (
                                self._with_schema3_projector_rows(
                                    local_group_model
                                )
                            )
                            if added_overlap_clusters:
                                coupling_data[
                                    "schema3_canonical_overlap_cluster_ids"
                                ] = added_overlap_clusters
                                coupling_data["coupled_cluster_ids"] = tuple(
                                    local_group_model.clusters.keys()
                                )
                        self.local_group_coupled_model[
                            target_root] = local_group_model
                        self.local_group_model[target_root] = local_group_model
                        self.local_groups[
                            target_root] = local_group_model.clusters
                        if schema4_ordered_construction:
                            self.local_group_coupling_matrix[
                                target_root
                            ] = coupling_data
                            self.local_group_coupling_matrix[
                                (target_root, family_label)
                            ] = coupling_data

                    self._set_local_group_datapoint_metadata(
                        core_dp,
                        family_label=family_label,
                        bank_role='core',
                        local_group_model=local_group_model,
                        rotor_ids=(),
                        reference_molecule=mol_basis[0],
                    )
                    core_dp.write_hdf5(target_file, core_label)
                    if write_org_database:
                        core_dp.write_hdf5(
                            f'im_database_{target_root}_org.h5', core_label)

                    for record in generated_cluster_records:
                        cluster = local_group_model.clusters[
                            record['cluster_id']]
                        cluster_dp = record['datapoint']
                        self._set_local_group_datapoint_metadata(
                            cluster_dp,
                            family_label=family_label,
                            bank_role='cluster',
                            local_group_model=local_group_model,
                            cluster=cluster,
                            cluster_state_id=record['state_id'],
                            rotor_ids=record['rotor_ids'],
                            dihedrals_to_rotate=record[
                                'dihedrals_to_rotate'],
                            phase_signature=record['phase_signature'],
                            reference_molecule=record['molecule'],
                            is_anchor=record['is_anchor'],
                        )
                        cluster_dp.write_hdf5(
                            target_file, record['label'])
                        if write_org_database:
                            cluster_dp.write_hdf5(
                                f'im_database_{target_root}_org.h5',
                                record['label'],
                            )

                    write_signed_factor_registry_for_family(
                        target_file,
                        target_root,
                        family_label,
                        local_group_model,
                        core_label,
                        point_labels_by_cluster,
                        local_factor_combination_mode=(
                            self.local_factor_combination_mode),
                    )

                    if self._schema4_msi_enabled():
                        self.construct_schema4_registry(
                            target_file=target_file,
                            root=target_root,
                            interpolation_settings=(
                                interpolation_settings[target_root]
                            ),
                        )

                    labels, _ = interpolation_driver.read_labels()

                    self.ostream.print_blank()
                    self.ostream.print_header(
                        "Local-group database expansion: Added family "
                        f"{family_label} for root {target_root}."
                    )
                    self.ostream.print_header(
                        f"  core: {core_label}"
                    )
                    primitive_count = len(coupling_data['primitive_cluster_ids'])
                    coupled_count = len(local_group_model.clusters)
                    self.ostream.print_header(
                        "  local-group coupling: "
                        f"{primitive_count} primitive clusters -> "
                        f"{coupled_count} active clusters "
                        f"(threshold {coupling_data['threshold']:.3f})"
                    )
                    self.ostream.print_header(
                        f"  local cluster states written: {len(point_labels) - 1}"
                    )
                    self.ostream.print_block(
                        f"Database expansion with {', '.join(labels)}")
                    self.ostream.flush()

    def add_point(self, molecule_specific_information, interpolation_settings, symmetry_information={}, write_org_database=True):
        """ Adds a new point to the database.

            :param molecule:
                the molecule.
            :param imforcefielddatafile:
                Datafile containing the information of the IM forcefield.

        """
        self.ostream.print_blank()
        self.ostream.print_header(f"Constructing initial Database {interpolation_settings[0]['imforcefield_file']} with the provided structure(s) and settings.")
        self.ostream.flush()
        if len(self.drivers) == 0:
            raise ValueError("No energy driver defined.")

        if (
            self.use_local_group_database
            and any(
                getattr(model, 'enabled', False)
                for model in self.local_group_model.values()
            )
        ):
            return self.add_point_local_groups(
                molecule_specific_information,
                interpolation_settings,
                symmetry_information=symmetry_information,
                write_org_database=write_org_database,
            )

        # define impesdriver to determine if stucture should be added:

        # For symmetry groups of periodicty of 3 it is crucial for the interpolation to set the dihedral to the position between 2 extreme points in order to account
        # for the symmetry correclty using only one reference point!

        # create all molecule combinations

        adjusted_molecule = {'gs': [], 'es': []}

        for entries in molecule_specific_information:
            molecule = entries[0]
            states = entries[2]
            symmetry_point = False

            if self.eq_bond_length is None:
                self.eq_bond_length = []
                for idx, element in enumerate(self.roots_z_matrix[0]['bonds']):

                    if len(element) == 2 and self.use_minimized_structures[0] and self.eq_bond_length_irc_bonds is not None and element not in self.eq_bond_length_irc_bonds:
                        self.eq_bond_length.append(molecule.get_distance([element[0] + 1, element[1] + 1], 'bohr'))
                    elif len(element) == 2 and self.use_minimized_structures[0] and self.eq_bond_length_irc_bonds is not None and element in self.eq_bond_length_irc_bonds:
                        self.eq_bond_length.append(molecule_specific_information[-1][0].get_distance([element[0] + 1, element[1] + 1], 'bohr'))
                    elif len(element) == 2 and self.use_minimized_structures[0]:
                        self.eq_bond_length.append(molecule.get_distance([element[0] + 1, element[1] + 1], 'bohr'))
                    elif len(element) == 2:
                        self.eq_bond_length.append(0.0)

            if 0 in entries[2]:
                adjusted_molecule['gs'].append((entries[0], entries[1], 1, None, [0], symmetry_point, entries[3]))
            if any(x > 0 for x in entries[2]):
                states = [state for state in entries[2] if state > 0]
                adjusted_molecule['es'].append((entries[0], entries[1], 1, None, states, symmetry_point, entries[3]))

        for key, entries in adjusted_molecule.items():
            if len(entries) == 0:
                continue
            drivers = self.drivers[key]

            label_counter = 0
            for mol_basis in entries:
                if not mol_basis[5]:
                    label_counter = 0

                imp_int_constraints = {'bonds': [], 'angles': [], 'dihedrals': [], 'impropers': []}
                for constraint in mol_basis[-1]:
                    if len(constraint) == 2:
                        imp_int_constraints["bonds"].append(constraint)
                    if len(constraint) == 3:
                        imp_int_constraints["angles"].append(constraint)
                    if len(constraint) == 4:
                        imp_int_constraints["dihedrals"].append(constraint)

                energies, scf_results, rsp_results = self._compute_energy(drivers[0], mol_basis[0], mol_basis[1])

                if isinstance(drivers[0], LinearResponseEigenSolver) or isinstance(drivers[0], TdaEigenSolver):
                    energies = energies[mol_basis[4]]

                gradients = self._compute_gradient(drivers[1], mol_basis[0], mol_basis[1], scf_results, rsp_results)
                hessians = self._compute_hessian(drivers[2], mol_basis[0], mol_basis[1])
                backend_metadata = self._collect_backend_metadata(drivers)

                inv_sqrt_masses = None
                if self.use_mass_weight:
                    masses = mol_basis[0].get_masses().copy()
                    masses_cart = np.repeat(masses, 3)
                    inv_sqrt_masses = 1.0 / np.sqrt(masses_cart)

                for number in range(len(energies)):
                    target_root = mol_basis[4][number]
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
                    label = None
                    grad = gradients[number].copy()
                    hess = hessians[number].copy()
                    grad_vec = grad.reshape(-1)         # (3N,)
                    hess_mat = hess.reshape(grad_vec.size, grad_vec.size)
                    mw_grad_vec = grad_vec
                    mw_hess_mat = hess_mat
                    if self.use_mass_weight:
                        mw_grad_vec = inv_sqrt_masses * grad_vec
                        mw_hess_mat = (inv_sqrt_masses[:, None] * hess_mat) * inv_sqrt_masses[None, :]

                    if not mol_basis[5]:
                        label = f'point_{len(sorted_labels) + 1}'

                    impes_coordinate = InterpolationDatapoint(z_matrix)
                    impes_coordinate.update_settings(interpolation_settings[target_root])
                    impes_coordinate.cartesian_coordinates = mol_basis[0].get_coordinates_in_bohr()
                    impes_coordinate.eq_bond_lengths = self.eq_bond_length
                    impes_coordinate.imp_int_coordinates = imp_int_constraints
                    impes_coordinate.inv_sqrt_masses = inv_sqrt_masses
                    impes_coordinate.metadata = deepcopy(backend_metadata)
                    impes_coordinate.energy = energies[number]
                    impes_coordinate.gradient = mw_grad_vec.reshape(grad.shape)
                    impes_coordinate.hessian = mw_hess_mat.reshape(hess.shape)
                    impes_coordinate.transform_gradient_and_hessian()

                    trust_radius = self.use_opt_confidence_radius[2]
                    impes_coordinate.confidence_radius = trust_radius

                    impes_coordinate.write_hdf5(target_file, label)

                    if write_org_database:
                        impes_coordinate.write_hdf5(
                            f'im_database_{target_root}_org.h5', label)
                    interpolation_driver.imforcefield_file = target_file

                    labels, z_matrix = interpolation_driver.read_labels()

                    self.ostream.print_blank()
                    self.ostream.print_header(f"Database expansion: Added point {label} to the database of root {target_root} with energy {energies[number]} Hartree.")
                    self.ostream.print_block(f"Database expansion with {', '.join(labels)}")
                    self.ostream.flush()

                label_counter += 1

    def _calculate_translation_coordinates(self, cart_coord):
        """Center the molecule by translating its geometric center to (0, 0, 0)."""
        center = np.mean(cart_coord, axis=0)
        translated_coordinates = cart_coord - center

        return translated_coordinates

    def _build_opt_constraint_list(self, constraints, index_offset=1):

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

    def _run_optimization(
            self,
            optimization_driver,
            molecule,
            constraints=None,
            transition=False,
            index_offset=1,
            compute_args=None,
            source_molecule=None,
            max_iter=None):

        opt_drv = OptimizationDriver(optimization_driver)
        opt_drv.ostream.mute()
        opt_drv.transition = transition
        if max_iter is not None:
            opt_drv.max_iter = int(max_iter)

        if constraints is not None:
            opt_drv.constraints = self._build_opt_constraint_list(constraints, index_offset=index_offset)

        if compute_args is None:
            opt_results = opt_drv.compute(molecule)
        else:
            opt_results = opt_drv.compute(molecule, *compute_args)

        return opt_results

    def _compute_energy(self, qm_driver, molecule, basis=None):
        """ Computes the QM energy using self.qm_driver.

            :param molecule:
                The molecule.
            :param basis:
                The basis set.

            :returns the QM energy.
        """
        # Dtermine the type of energy driver, to be able to
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

        elif isinstance(qm_driver, GxtbDriver):
            qm_driver.ostream.mute()
            result = qm_driver.compute(molecule)
            qm_driver.ostream.unmute()
            qm_driver.energy_metadata = result.metadata
            qm_energy = np.array([result.energy])

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
        """ Computes the QM gradient using self.grad_driver.

            :param molecule:
                The molecule.
            :param basis:
                The basis set.

            :returns the QM gradient.
        """

        qm_gradient = None

        if isinstance(grad_driver, XtbGradientDriver):

            grad_driver.ostream.mute()
            grad_driver.compute(molecule)
            qm_gradient = grad_driver.gradient
            qm_gradient = np.array([qm_gradient])
            grad_driver.ostream.unmute()

        elif isinstance(grad_driver, GxtbGradientDriver):
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
        """ Computes the QM Hessian using self.hess_driver.

            :param molecule:
                The molecule.
            :param basis:
                The basis set.

            :returns the QM Hessian matrix.
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

        elif isinstance(hess_driver, GxtbHessianDriver):
            hess_driver.ostream.mute()
            hess_driver.compute(molecule)
            qm_hessian = hess_driver.hessian
            hess_driver.ostream.unmute()

            if qm_hessian is None:
                raise RuntimeError(
                    'g-XTB Hessian is None on this rank after MPI synchronization.')
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

    def _append_confirm_database_quality_h5(
        self,
        h5f: h5py.File,
        molecules,
        qm_energies,
        im_energies,
        state: int,         # optional; keep None if you don't have it here
        check_meta: dict | None = None,  # optional metadata (step, tag, etc.)
        compression="gzip",
        compression_opts=4,
    ):
        """
        Append a batch of 'confirmed reference structures' into the *already open* HDF5 file `h5f`.

        Layout (example):
        /confirm_database_quality#000012/
            attrs: created_by, state, ...
            /state_<state>/
                qm_energy    (N,) float64
                im_energy    (N,) float64
                distance     (N,) float64   [optional: only created if distances provided]
                natoms       (N,) int32
                xyz          (N,) utf-8 string
                coords_flat  (N,) vlen<float64>  (flattened natoms*3; empty if not available)
        """
        state = int(state)
        n = len(molecules)

        if not (len(qm_energies) == len(im_energies) == n):
            raise ValueError("molecules, qm_energies, im_energies must have identical length.")

        # --------- choose next batch group name: confirm_database_quality#000000, #000001, ...
        prefix = "confirm_database_quality#"
        existing = [k for k in h5f.keys() if isinstance(k, str) and k.startswith(prefix)]
        if existing:
            # parse numeric suffixes safely
            ids = []
            for k in existing:
                suf = k[len(prefix):]
                try:
                    ids.append(int(suf))
                except ValueError:
                    pass
            next_id = (max(ids) + 1) if ids else 0
        else:
            next_id = 0

        batch_group_name = f"{prefix}{next_id:06d}"
        batch_grp = h5f.require_group(batch_group_name)

        # add minimal metadata
        batch_grp.attrs["state"] = state
        if check_meta:
            for kk, vv in check_meta.items():
                # store simple scalar/string metadata
                if isinstance(vv, (str, int, float, np.integer, np.floating)):
                    batch_grp.attrs[str(kk)] = vv

        # store under a per-state subgroup (keeps it consistent with your file conventions)
        state_grp = batch_grp.require_group(f"state_{state}")

        str_dt = h5py.string_dtype(encoding="utf-8")
        vlen_f64 = h5py.vlen_dtype(np.dtype("float64"))

        def req_1d(name, dtype):
            if name in state_grp:
                return state_grp[name]
            return state_grp.create_dataset(
                name,
                shape=(0,),
                maxshape=(None,),
                dtype=dtype,
                chunks=True,
                compression=compression,
                compression_opts=compression_opts,
                shuffle=True,
            )

        qm_ds = req_1d("qm_energy", np.float64)
        im_ds = req_1d("im_energy", np.float64)
        na_ds = req_1d("natoms", np.int32)
        xyz_ds = req_1d("xyz", str_dt)
        cf_ds = req_1d("coords_flat", vlen_f64)

        # --------- resize and append
        N0 = qm_ds.shape[0]
        # N1 = N0 + n

        for i, mol in enumerate(molecules):
            xyz = mol.get_xyz_string()

            # coords optional
            coords = None
            if hasattr(mol, "get_coordinates"):
                coords = np.asarray(mol.get_coordinates_in_angstrom(), dtype=np.float64)

            # natoms robust inference
            if hasattr(mol, "n_atoms"):
                natoms = int(mol.n_atoms)
            elif coords is not None:
                natoms = int(coords.shape[0])
            else:
                natoms = int(xyz.splitlines()[0].strip())

            row = N0 + i
            row = qm_ds.shape[0]  # current length
            qm_ds.resize((row + 1,))  # grow by 1
            im_ds.resize((row + 1,))  # grow by 1
            na_ds.resize((row + 1,))  # grow by 1
            cf_ds.resize((row + 1,))  # grow by 1
            xyz_ds.resize((row + 1,))  # grow by 1

            qm_ds[row] = float(qm_energies[i])
            im_ds[row] = float(im_energies[i])
            na_ds[row] = natoms
            xyz_ds[row] = xyz

            if coords is None:
                cf_ds[row] = np.array([], dtype=np.float64)
            else:
                if coords.shape != (natoms, 3):
                    raise ValueError(f"Coords must be shape (natoms,3), got {coords.shape}.")
                cf_ds[row] = coords.reshape(-1)

        return batch_group_name  # handy for logging/debugging
