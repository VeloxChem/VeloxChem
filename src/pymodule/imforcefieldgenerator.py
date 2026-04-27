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
import itertools
import os
import random
import h5py
import builtins

from mpi4py import MPI

from contextlib import redirect_stderr
from io import StringIO
from copy import deepcopy
import json


from .xtbdriver import XtbDriver
from .xtbgradientdriver import XtbGradientDriver
from .xtbhessiandriver import XtbHessianDriver
from .scfrestdriver import ScfRestrictedDriver
from .scfunrestdriver import ScfUnrestrictedDriver
from .scfgradientdriver import ScfGradientDriver
from .scfhessiandriver import ScfHessianDriver
from .tdaeigensolver import TdaEigenSolver
from .lreigensolver import LinearResponseEigenSolver
from .tddftgradientdriver import TddftGradientDriver
from .tddfthessiandriver import TddftHessianDriver

from .molecularbasis import MolecularBasis
from .interpolationdriver import InterpolationDriver
from .interpolationdatapoint import InterpolationDatapoint
from .imdatabasepointcollecter import IMDatabasePointCollecter
from .mmforcefieldgenerator import MMForceFieldGenerator
from .conformergenerator import ConformerGenerator
from .optimizationdriver import OptimizationDriver
# from .atommapper import AtomMapper
from .tsguesser import TransitionStateGuesser

from .molecule import Molecule
from .errorhandler import assert_msg_critical
from .veloxchemlib import mpi_master, hartree_in_kcalpermol, bohr_in_angstrom

from .imrotorbuilder import (
    RotorDefinition,
    RotorClusterDefinition,
    RotorClusterInformation,
    RotorClusterStateDefinition,
    RotorClusterAngleLibrary,
    build_rotor_clusters,
    rotor_coupling_score,
)

with redirect_stderr(StringIO()) as fg_err:
    import geometric

# def _rank_root_print(*args, **kwargs):
#     """
#     Print only on MPI root rank to avoid delayed/duplicated worker stdout noise.
#     Falls back to normal print when MPI rank detection is unavailable.
#     """
#     try:
#         if MPI.COMM_WORLD.Get_rank() == mpi_master():
#             builtins.print(*args, **kwargs)
#     except Exception:
#         builtins.print(*args, **kwargs)



# def _collective_call_inline_style(comm, rank, root, phase_name, fn, *args, **kwargs):
#             comm.barrier()
#             local_err = None
#             result = None
#             try:
#                 result = fn(*args, **kwargs)
#             except Exception as exc:
#                 local_err = f"rank {rank}: {exc}"
#             errs = [e for e in comm.allgather(local_err) if e is not None]
#             if errs:
#                 if rank == root:
#                     raise RuntimeError(f"{phase_name} failed\n" + "\n".join(errs))
#                 raise RuntimeError(f"{phase_name} failed on another rank")
#             comm.barrier()
#             return result

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

        - allowed_deviation: A threshold for how much a generated or sampled 
                            structure is allowed to deviate from an expected configuration. Ensures structural integrity 
                            during sampling or dynamics.

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
    """
    
    def __init__(self, ground_state_driver=None, excited_state_driver=None, roots_to_follow=[0], rsp_method='tda'):
        
        assert_msg_critical(all(roots_to_follow) == 0, "The current version is restricted to ground-state potentials. Later version will allow multi-state construction!")
        
        self.open_mm_platform = None

        self.density_of_datapoints = None
        self.qm_data_points = None
        self.qmlabels = None

        self.molecules_along_rp = None
        self.conformal_structures = None
        self.sampling_structures = 1
        
        self.excitation_pulse = None
        self.state_specific_molecules = None
        self.datafile = None
        self.dihedrals_dict = None
        self.allowed_deviation = None

        self.roots_z_matrix = {}
        self.int_coord_bond_information = None
        self.symmetry_information = None
        self.symmetry_rotors = None
        self.rotor_corr_threshold = 0.01
        self.symmetry_dihedral_lists = {}
        self.all_rotatable_bonds = None

        self.reaction_structures = None
        self.seed_structures = None
        self.eq_bond_length = None
        self.eq_bond_length_irc_bonds = None

        # Here us the driver set up
        self.gs_basis_set_label = 'def2-svp'
        self.es_basis_set_label = '6-31g*'

        self.drivers = {'gs': None, 'es':None}
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
        ##########################################################
        ################# External Settings ######################
        ##########################################################

        if isinstance(excited_state_driver, ScfRestrictedDriver):
            
            excited_state_response_driver = TdaEigenSolver()
            if rsp_method == 'tddft':
                excited_state_response_driver = LinearResponseEigenSolver()
            excited_state_response_driver.nstates = 10
            excited_state_grad_driver = TddftGradientDriver(excited_state_driver)
            excitated_roots = [root for root in roots_to_follow if root != 0]
            excited_state_grad_driver.state_deriv_index = excitated_roots
            excited_state_hess_driver = TddftHessianDriver(excited_state_driver, excited_state_response_driver, excited_state_grad_driver)

            self.drivers['es'] = (excited_state_response_driver, excited_state_grad_driver, excited_state_hess_driver, excited_state_driver)


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
        self.exponent_p = 2
        self.exponent_q = 2
        self.confidence_radius = 0.5
        self.imforcefieldfiles = None
        self.use_inverse_bond_length = True
        self.use_eq_bond_length = False
        self.eq_bond_symmetry_mode = "masked_exact"  # "masked_exact" | "symmetrized"
        self.use_cosine_dihedral = False
        self.use_tc_weights = True

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
        self.bias_force_reaction_prop = None # this is being set by giving a dihedral, force constant and the final theta, steps_when_increased

        # sampling settings
        self.sampling_settings = {
            'enabled':False,
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
        self.im_results = {'n_datapoints': None, 'RMSD': None, '|D|':None} 


        # confirm database quality
        self.nstruc_to_confirm_database_quality = 50

        
        # set boolean for the optimization features of the metho
        self.use_minimized_structures = [True, [], []] # to use minimum structures, specific constraints, root of the constraints
        self.add_conformal_structures = True # Add all conformal structures which are being generated for the input molecule
        self.use_symmetry = False # symmetry is considerd only for CH3 rotors for now
        self.cluster_run = False # use the cluster formulation in order to allow individual symmetry groups to be independent --> drastic point reduction
        self.identfy_relevant_int_coordinates = True # This constraint-optimizes new datapoints during construction run to presever optimal smoothness
        self.use_opt_confidence_radius = [False, 'multi_grad', 0.5, 0.3] # optimize the set of Trust radii of each datapoint given a set of references (energy, gradient))
        self.exclude_non_core = False # use only core atoms (H exclusion) for the interpolation (seemed to be less stable)

        self.imp_int_coordinates = []

        self.profile_runtime_timing = False
        self.profile_interpolation_timing = False

        # self.use_local_preload = False
        # self.local_preload_workers = 2
        # self.local_preload_min_tasks = 8
        # self.local_preload_omp_threads = 1

        # self.use_outer_parallel = True
        # self.outer_parallel_workers = 10  # or 0 for auto in driver
        # self.outer_parallel_min_labels = 8
        # self.outer_parallel_chunk_size = 0
        
        # # mpi section of the code Xin needs to check as the code differs from normal MPI integration frameworks
        # self.use_mpi_preload = False
        # self.mpi_control_plane_enabled = False
        # self.mpi_root_worker_mode = False
        # self.mpi_reload_from_hdf5 = False
        # self.mpi_debug_sync = False

        self._test_hooks = {}
        self._test_hook_strict = False
        self._test_event_counter = 0
        self._test_run_id = None

    def _to_test_payload(self, obj):
        """
        Convert payload content to JSON-friendly Python-native types.

        Why:
        - tests should not handle numpy scalar/array types directly,
        - payloads must be serializable for debugging and artifact dumps.
        """
        if isinstance(obj, dict):
            return {str(k): self._to_test_payload(v) for k, v in obj.items()}
        if isinstance(obj, (list, tuple)):
            return [self._to_test_payload(v) for v in obj]
        if isinstance(obj, (np.floating, np.integer)):
            return obj.item()
        if hasattr(obj, "tolist"):
            return self._to_test_payload(obj.tolist())
        return obj

    def _event_envelope(self, event_name, payload):
        """
        Build a stable event envelope consumed by integration tests.
        """
        self._test_event_counter += 1
        return {
            "event": str(event_name),
            "seq": int(self._test_event_counter),
            "run_id": self._test_run_id,
            "roots_to_follow": list(self.roots_to_follow),
            "payload": self._to_test_payload(payload),
        }

    def _emit_test_hook(self, event_name, payload):
        """
        Emit one event to test callbacks.

        Behavior:
        - emits only on MPI root rank (avoids duplicate test events),
        - supports exact-name callback and wildcard callback (`*`),
        - optional strict-mode failure.
        """
        comm = MPI.COMM_WORLD
        if comm.Get_rank() != mpi_master():
            return

        hook_map = getattr(self, "_test_hooks", {})
        fn = hook_map.get(event_name)
        wildcard = hook_map.get("*")
        if fn is None and wildcard is None:
            return

        envelope = self._event_envelope(event_name, payload)

        for cb in (fn, wildcard):
            if cb is None:
                continue
            try:
                cb(envelope)
            except Exception as exc:
                if self._test_hook_strict:
                    raise RuntimeError(
                        f"IMForceFieldGenerator test hook failed for event '{event_name}': {exc}"
                    ) from exc
    
    def _mode_signature(self):
        """
        Return compact mode info used to parametrize tests and debug mismatches.
        """
        if self.use_inverse_bond_length:
            bond_mode = "inverse"
        elif self.use_eq_bond_length:
            bond_mode = "equilibrium"
        else:
            bond_mode = "plain_r"

        return {
            "bond_mode": bond_mode,
            "weightfunction_type": str(self.weightfunction_type),
            "use_cosine_dihedral": bool(self.use_cosine_dihedral),
            "use_tc_weights": bool(self.use_tc_weights),
            "eq_bond_symmetry_mode": str(self.eq_bond_symmetry_mode),
        }

    def _label_count_per_root(self):
        """
        Read current number of labels from each root HDF5.
        Used by tests to assert add-point side effects.
        """
        counts = {}
        for root in self.roots_to_follow:
            fpath = self.states_interpolation_settings[root]["imforcefield_file"]
            if not os.path.exists(fpath):
                counts[int(root)] = 0
                continue
            drv = InterpolationDriver(self.roots_z_matrix[root])
            drv.update_settings(self.states_interpolation_settings[root])
            drv.imforcefield_file = fpath
            labels, _ = drv.read_labels()
            counts[int(root)] = int(len(labels))
        return counts
            
    # def _compute_energy_mpi_safe(self, qm_driver, molecule, basis=None, collective=True, phase_name='Energy calcualtion'):
    #     comm = MPI.COMM_WORLD
    #     rank = comm.Get_rank()
    #     root = mpi_master()

    #     if collective:
    #         return _collective_call_inline_style(
    #             comm, rank, root, phase_name, self._compute_energy, qm_driver, molecule, basis)

    #     return self._compute_energy(qm_driver, molecule, basis)

    # def _compute_gradient_mpi_safe(self, grad_driver, molecule, basis=None, scf_results=None, rsp_results=None, collective=True, phase_name='Gradient Calculation'):
    #     comm = MPI.COMM_WORLD
    #     rank = comm.Get_rank()
    #     root = mpi_master()

    #     if collective:
    #         return _collective_call_inline_style(
    #             comm, rank, root, phase_name, self._compute_gradient, grad_driver, molecule, basis, scf_results, rsp_results)

    #     return self._compute_gradient(grad_driver, molecule, basis, scf_results, rsp_results)

    # def _compute_hessian_mpi_safe(self, hess_driver, molecule, basis=None, collective=True, phase_name='hessian calculation'):
    #     comm = MPI.COMM_WORLD
    #     rank = comm.Get_rank()
    #     root = mpi_master()

    #     if collective:
    #         return _collective_call_inline_style(
    #             comm, rank, root, phase_name, self._compute_hessian, hess_driver, molecule, basis)

    #     return self._compute_hessian(hess_driver, molecule, basis)

    # def _run_optimization_mpi_safe(self, optimization_driver, molecule, constraints=None, transition=False, index_offset=1, compute_args=None, source_molecule=None, collective=True, phase_name='optimization tag'):
    #     comm = MPI.COMM_WORLD
    #     rank = comm.Get_rank()
    #     root = mpi_master()

    #     if collective:
    #         return _collective_call_inline_style(
    #             comm,
    #             rank,
    #             root,
    #             phase_name,
    #             self._run_optimization,
    #             optimization_driver,
    #             molecule,
    #             constraints,
    #             transition,
    #             index_offset,
    #             compute_args,
    #             source_molecule,
    #         )

    #     return self._run_optimization(
    #         optimization_driver,
    #         molecule,
    #         constraints=constraints,
    #         transition=transition,
    #         index_offset=index_offset,
    #         compute_args=compute_args,
    #         source_molecule=source_molecule,
    #     )

    # @staticmethod
    # def _serialize_molecule_for_mpi(mol):
    #     return (mol.get_xyz_string(), int(mol.get_charge()), int(mol.get_multiplicity()))

    # @staticmethod
    # def _deserialize_molecule_from_mpi(payload):
    #     mol = Molecule.from_xyz_string(payload[0])
    #     mol.set_charge(int(payload[1]))
    #     mol.set_multiplicity(int(payload[2]))
    #     return mol

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

    def _run_optimization(self, optimization_driver, molecule, constraints=None, transition=False, index_offset=1, compute_args=None, source_molecule=None):

        opt_drv = OptimizationDriver(optimization_driver)
        opt_drv.ostream.mute()
        opt_drv.transition = transition

        if constraints is not None:
            opt_drv.constraints = self._build_opt_constraint_list(constraints, index_offset=index_offset)

        if compute_args is None:
            opt_results = opt_drv.compute(molecule)
        else:
            opt_results = opt_drv.compute(molecule, *compute_args)

        # if source_molecule is None:
        #     source_molecule = molecule

        # optimized_molecule = Molecule.from_xyz_string(opt_results['final_geometry'])
        # optimized_molecule.set_charge(source_molecule.get_charge())
        # optimized_molecule.set_multiplicity(source_molecule.get_multiplicity())

        return opt_results

    # def _generate_conformers_mpi_synchronized(self, molecule):
    #     """
    #     Run ConformerGenerator as an explicit MPI phase where all ranks enter
    #     the same collective flow and then receive one canonical payload.
    #     """
    #     comm = MPI.COMM_WORLD
    #     rank = comm.Get_rank()
    #     root = mpi_master()

    #     original_cwd = os.getcwd()
    #     rank_workdir = os.path.join(original_cwd,
    #                                 f".imff_conformer_rank_{rank}")
    #     os.makedirs(rank_workdir, exist_ok=True)

    #     local_error = None
    #     conformers_root = None
    #     dihedral_candidates = []

    #     try:
    #         # Enter a dedicated conformer phase together.
    #         comm.barrier()
    #         os.chdir(rank_workdir)
    #         conformer_generator = ConformerGenerator(comm=comm)
    #         conformers_root = conformer_generator.generate(molecule)
    #         dihedral_candidates = list(
    #             getattr(conformer_generator, "dihedral_candidates", []) or [])
    #         comm.barrier()
    #     except Exception as exc:
    #         local_error = f"rank {rank}: {exc}"
    #     finally:
    #         os.chdir(original_cwd)

    #     all_errors = comm.allgather(local_error)
    #     all_errors = [err for err in all_errors if err is not None]
    #     if all_errors:
    #         if rank == root:
    #             err_msg = "\n".join(all_errors)
    #             raise RuntimeError(
    #                 "ConformerGenerator MPI phase failed.\n"
    #                 f"Reported rank errors:\n{err_msg}")
    #         raise RuntimeError(
    #             "ConformerGenerator MPI phase failed on another rank.")

    #     payload = None
    #     if rank == root:
    #         molecules = []
    #         if isinstance(conformers_root, dict):
    #             molecules = list(conformers_root.get("molecules", []) or [])
    #         payload = {
    #             "molecules_xyz": [mol.get_xyz_string() for mol in molecules],
    #             "dihedral_candidates": deepcopy(dihedral_candidates),
    #         }

    #     payload = comm.bcast(payload, root=root)
    #     comm.barrier()

    #     conformers = {"molecules": []}
    #     for xyz in payload.get("molecules_xyz", []):
    #         conf_mol = Molecule.from_xyz_string(xyz)
    #         conf_mol.set_charge(molecule.get_charge())
    #         conf_mol.set_multiplicity(molecule.get_multiplicity())
    #         conformers["molecules"].append(conf_mol)

    #     dihedral_candidates = list(payload.get("dihedral_candidates", []) or [])
    #     return conformers, dihedral_candidates

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

    def _build_rotor_defintions(self, z_matrix, symmetry_dihedral_lists, dihedral_start):
        """
        This function creates a list of Rotors (CH3) groups within the molecule to
        account for suymmetry handling withnin the interpolation decreasing the number
        of necessary points drastically.
        """
        rotors = {}
        rotor_id = 0
        
        sym3_groups = symmetry_dihedral_lists[3]

        dihedral_index_map = {
            tuple(int(x) for x in dih): dihedral_start + idx for idx, dih in enumerate(z_matrix['dihedrals']) 
        }
        print(dihedral_index_map, sym3_groups)

        for center, dihedral_list in sym3_groups.items():

            torsion_rows = []
            torsion_coords = []
            atoms = set()

            row_coord_pairs = []
            for dih in dihedral_list:
                key = tuple(int(x) for x in dih)
                row = dihedral_index_map[key]
                row_coord_pairs.append((row, key))
                atoms.update(key)

            row_coord_pairs.sort(key=lambda x: x[0])
            torsion_rows = tuple(row for row, _ in row_coord_pairs)
            torsion_coords = tuple(coord for _, coord in row_coord_pairs)

            rotors[rotor_id] = RotorDefinition(
                rotor_id=rotor_id,
                center=tuple(int(x) for x in center),
                torsion_rows=torsion_rows,
                torsion_coords=torsion_coords,
                symmetry_order=3,
                atom_group=tuple(sorted(atoms)),
            )
            
            rotor_id += 1
        
        return rotors

    
    def set_up_the_system(self, molecule, extract_z_matrix=None):

        """
        Assign the neccessary variables with respected values. 

        :param molecule: original molecule

        :param target_dihedrals: is a list of dihedrals that should be scanned during the dynamics

        :param sampling_structures: devides the searchspace around given rotatbale dihedrals
            
        """
        def _determine_cX_symmetry_groups(molecule):
            """
            Build symmetry groups for methyl rotors only.

            Returns the same tuple shape used from AtomMapper:
                (atom_map, symmetry_groups, rc_groups)
            where:
                atom_map: all atom indices
                symmetry_groups: list of equivalent H triplets in CH3 groups
                rc_groups: empty (not used in this CH3-only path)
            """
            labels = [str(x).strip().capitalize() for x in molecule.get_labels()]
            conn = molecule.get_connectivity_matrix()
            n_atoms = len(labels)

            atom_map = list(range(n_atoms))
            rc_groups = []

            ch3_groups = []
            seen = set()

            allowed_neighbours = ["H", "F", "Cl", "Br"]

            for c_idx, symbol in enumerate(labels):
                if symbol != "C":
                    continue

                neighbors = [j for j, bonded in enumerate(conn[c_idx]) if bonded]
                h_neighbors = [j for j in neighbors if labels[j] in allowed_neighbours]
                heavy_neighbors = [j for j in neighbors if labels[j] not in allowed_neighbours]

                # Rotor-relevant methyl carbon: C(H)3-X
                if len(neighbors) == 4 and len(h_neighbors) == 3 and len(heavy_neighbors) == 1:
                    group = tuple(sorted(h_neighbors))
                    if group not in seen:
                        seen.add(group)
                        ch3_groups.append(list(group))

            return atom_map, ch3_groups, rc_groups


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
                            for atom in subgroup:
                                print(molecule.get_labels()[atom], sum(conn[atom]))
                            rot_groups[state].append(sorted(subgroup))
                            new_groups[state].append(sorted(subgroup))


            return new_groups, rot_groups
        
        def _promote_nonrotatable_ring_torsions_to_impropers(zmat, ff_gen, rotatable_bonds_zero_based):
            """
            Promote proper ring torsions around non-rotatable bonds to impropers.

            Matching by full 4-atom set is too strict (and often fails), so we match
            by the middle bond of the proper torsion and collect impropers touching
            the same bond.
            """
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
        
        for root in self.roots_to_follow:
            if extract_z_matrix is None:
                extract_z_matrix = {}
                extract_z_matrix[root] = False
            if root not in self.roots_z_matrix and not extract_z_matrix[root]:
                # if no database is provided construct the primitive internal coordinates using geometric
                self.roots_z_matrix[root] = self.define_z_matrix_dict(molecule)
                if self.reaction_structures is not None:
                    merge_info = self.merge_reaction_internal_coordinates(
                        reaction_structures=self.reaction_structures,
                        root=self.roots_to_follow[0],
                        include_existing_root_zmat=True,
                        forced_coordinates=self.reaction_forced_coordinates if hasattr(self, "reaction_forced_coordinates") else None,
                    )

                    if self.use_eq_bond_length:
                        self.eq_bond_length_irc_bonds = merge_info['added_coordinates']['bonds']
                    self.roots_z_matrix[self.roots_to_follow[0]] = merge_info["global_z_matrix"]  
            elif root not in self.roots_z_matrix:
                # generate the z-matrix based for the interpolation database provided
                int_driver = InterpolationDriver()
                int_driver.update_settings({ 'interpolation_type':self.interpolation_type,
                                        'weightfunction_type':self.weightfunction_type,
                                        'exponent_p':self.exponent_p,
                                        'exponent_q':self.exponent_q, 
                                        'confidence_radius':self.confidence_radius,
                                        'imforcefield_file':self.imforcefieldfiles[root],
                                        'use_inverse_bond_length':self.use_inverse_bond_length,
                                        'use_eq_bond_length':self.use_eq_bond_length,
                                        'eq_bond_symmetry_mode':self.eq_bond_symmetry_mode,
                                        'use_cosine_dihedral':self.use_cosine_dihedral,
                                        'use_tc_weights':self.use_tc_weights
                                    })
       
                _, z_matrix = int_driver.read_labels()
                self.roots_z_matrix[root] = z_matrix

            # The AtomMapper class is based on Turtlemap program which is used to
            # determine equivalent atoms within a molecular structure
            # atom_mapper = AtomMapper(molecule, molecule)
            # symmetry_groups = atom_mapper.determine_symmetry_group()

            symmetry_groups = _determine_cX_symmetry_groups(molecule)

            if not self.use_symmetry:
                symmetry_groups = (symmetry_groups[0], [], symmetry_groups[2])
            
            ff_gen = MMForceFieldGenerator()
            ff_gen.ostream.mute()
            ff_gen.partial_charges = molecule.get_partial_charges(molecule.get_charge())
            ff_gen.create_topology(molecule)

            rotatable_bonds = deepcopy(ff_gen.rotatable_bonds)
            org_rotatable_bonds = deepcopy(ff_gen.rotatable_bonds)


            # Work in zero-based indexing (same convention as z-matrix dihedrals)
            # and remove all symmetry-related rotatable bonds from the scan list.
            rotatable_bonds_zero_based = [tuple(sorted((i - 1, j - 1))) for (i, j) in rotatable_bonds]

            self.roots_z_matrix[root] = _promote_nonrotatable_ring_torsions_to_impropers(self.roots_z_matrix[root], ff_gen, rotatable_bonds_zero_based)

            dihedral_start = len(self.roots_z_matrix[root]['bonds']) + len(self.roots_z_matrix[root]['angles'])
            dihedral_end = dihedral_start + len(self.roots_z_matrix[root]['dihedrals'])

            self.all_rotatable_bonds = rotatable_bonds_zero_based
            all_exclision = [element for rot_bond in rotatable_bonds_zero_based for element in rot_bond]
            
            symmetry_groups_ref = [groups for groups in symmetry_groups[1] if not any(item in all_exclision for item in groups)]

            # reduce the symmetry to only CH3 or CH2 symmetry groups for the time being
            regrouped, rot_groups = regroup_by_rotatable_connection(molecule, symmetry_groups_ref, rotatable_bonds_zero_based, molecule.get_connectivity_matrix())
            dih_list = []
            if root == 0 or root == 1 and self.drivers['es'] is not None and self.drivers['es'][0].spin_flip:
            
                non_core_atoms = [element for group in regrouped['gs'] for element in group]
                core_atoms = [element for element in symmetry_groups[0] if element not in non_core_atoms]
                angles_to_set, _, _, self.symmetry_dihedral_lists, dih_list = self._adjust_symmetry_dihedrals(molecule, rot_groups['gs'], rotatable_bonds_zero_based, self.roots_z_matrix[root])
                dihedrals_to_set = {key: [] for key in angles_to_set.keys()}
                indices_list = []
                for key, dihedral_list in self.symmetry_dihedral_lists.items():
                
                    for i, element in enumerate(self.roots_z_matrix[root]['dihedrals']):
                        if tuple(sorted(element)) in dihedral_list:
                            indices_list.append(i)
                self.symmetry_information['gs'] = [symmetry_groups[0], rot_groups['gs'], regrouped['gs'], core_atoms, non_core_atoms, rotatable_bonds_zero_based, indices_list, self.symmetry_dihedral_lists, dihedrals_to_set, [dihedral_start, dihedral_end]]
            if root == 1 and self.drivers['es'] is not None and self.drivers['es'][0].spin_flip is False:
                non_core_atoms = [element for group in regrouped['gs'] for element in group]
                core_atoms = [element for element in symmetry_groups[0] if element not in non_core_atoms]
                angles_to_set, _, _, self.symmetry_dihedral_lists, dih_list = self._adjust_symmetry_dihedrals(molecule, rot_groups['gs'], rotatable_bonds_zero_based,  self.roots_z_matrix[root])
                dihedrals_to_set = {key: [] for key in angles_to_set.keys()}
                indices_list = []
                for key, dihedral_list in self.symmetry_dihedral_lists.items():
                
                    for i, element in enumerate(self.roots_z_matrix[root]['dihedrals']):
                        if tuple(sorted(element)) in dihedral_list:
                            indices_list.append(i)
                self.symmetry_information['es'] = [symmetry_groups[0], rot_groups['es'], regrouped['es'], core_atoms, non_core_atoms, rotatable_bonds_zero_based, indices_list, self.symmetry_dihedral_lists, dihedrals_to_set, [dihedral_start, dihedral_end]]
            if root >= 2 and len(self.symmetry_information['es']) == 0 and self.drivers['es']:
                non_core_atoms = [element for group in regrouped['es'] for element in group]
                core_atoms = [element for element in symmetry_groups[0] if element not in non_core_atoms]
                print(rot_groups['es'], rotatable_bonds_zero_based)
                angles_to_set, _, _, self.symmetry_dihedral_lists, dih_list = self._adjust_symmetry_dihedrals(molecule, rot_groups['es'], rotatable_bonds_zero_based,  self.roots_z_matrix[root])
                dihedrals_to_set = {key: [] for key in angles_to_set.keys()}
                indices_list = []
                for key, dihedral_list in self.symmetry_dihedral_lists.items():
                
                    for i, element in enumerate(self.roots_z_matrix[root]['dihedrals']):
                        if tuple(sorted(element)) in dihedral_list:
                            indices_list.append(i)
                self.symmetry_information['es'] = [symmetry_groups[0], rot_groups['es'], regrouped['es'], core_atoms, non_core_atoms, rotatable_bonds_zero_based, indices_list, self.symmetry_dihedral_lists, dihedrals_to_set, [dihedral_start, dihedral_end]]

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

        self.symmetry_rotors = self._build_rotor_defintions(self.roots_z_matrix[root], dih_list, dihedral_start)

        if self.reaction_structures is not None:
            self.seed_structures = self.build_initial_seed_structures(
                                    molecule=molecule,
                                    reaction_structures=self.reaction_structures,
                                    include_conformers=False,
                                    reaction_root=self.roots_to_follow[0],
                                    reaction_key=None,
                                    )
            
        if self.add_conformal_structures and 0 in self.roots_to_follow:

            conformers_plus_ts = {0 : {}}
        
            # conformal_structures, dihedral_candidates = self._generate_conformers_mpi_synchronized(molecule)
            conformer_gen = ConformerGenerator()
            conformer_gen.resp_charges = False
            
            conformer_results = conformer_gen.generate(molecule)
            dihedral_candidates = list(getattr(conformer_gen, "dihedral_candidates", []) or [])

            conformal_structures = {'molecules': list(conformer_results.get("molecules", []) or [])}
        

        # conformers = {"molecules": []}
        # for xyz in molecules.get("molecules_xyz", []):
        #     conf_mol = Molecule.from_xyz_string(xyz)
        #     conf_mol.set_charge(molecule.get_charge())
        #     conf_mol.set_multiplicity(molecule.get_multiplicity())
        #     conformers["molecules"].append(conf_mol)
            
            # conformal_structures, dihedral_candidates

            # comm = MPI.COMM_WORLD
            # comm.barrier()
            # rank = comm.Get_rank()
            # root = mpi_master()

            conformers_payload = None

            # before adding the conformers generate structures along the dihedral path
            # for making the energy path avaliable in the database
            # if rank == root:
            print('dihedral candidates', dihedral_candidates)
            if len(dihedral_candidates) > 0 and len(conformal_structures['molecules']) > 0:
                seed_molecule = conformal_structures['molecules'][0]
                rot_bond_set = {
                    tuple(sorted((int(a), int(b))))
                    for a, b in org_rotatable_bonds
                }
                def _circ_err_deg(actual, target):
                    return abs(((actual - target + 180.0) % 360.0) - 180.0)
                for dih_zero, angle_grid in dihedral_candidates:
                    dih_key = tuple(int(x) + 1 for x in dih_zero)
                    if tuple(sorted((dih_key[1], dih_key[2]))) not in rot_bond_set:
                        print('skipping dihedral', dih_zero)
                        continue
                    scan_angles = sorted({float(a) % 360.0 for a in angle_grid})
                    conformers_plus_ts[0][dih_key] = []
                    def _w360(a):
                        return float(a) % 360.0
                    def _k(a):
                        return round(_w360(a), 6)
                    def _signed_deg(a):
                        a = _w360(a)
                        return a - 360.0 if a > 180.0 else a
                    # if the same angle is generated multiple times, keep the strongest tag
                    tag_priority = {'constraint': 0, 'transition': 1, 'normal': 2}
                    angle_to_tag = {}
                    def _add_angle(angle, tag):
                        key = _k(angle)
                        old = angle_to_tag.get(key)
                        if old is None or tag_priority[tag] > tag_priority[old]:
                            angle_to_tag[key] = tag
                    # minima from conformer candidates
                    for a in scan_angles:
                        _add_angle(a, 'normal')
                    # for each minima pair on the circle:
                    # - TS midpoint (50%)
                    # - minima↔TS midpoint (25%)
                    # - TS↔next-minima midpoint (75%)
                    if len(scan_angles) > 1:
                        for idx, cur in enumerate(scan_angles):
                            nxt = scan_angles[(idx + 1) % len(scan_angles)]
                            step = (nxt - cur) % 360.0
                            if step < 1.0e-8:
                                continue
                            _add_angle(cur + 0.50 * step, 'transition')
                            _add_angle(cur + 0.25 * step, 'constraint')
                            _add_angle(cur + 0.75 * step, 'constraint')
                    schedule = [(ang, angle_to_tag[ang]) for ang in sorted(angle_to_tag.keys())]
                    print(
                        f"scan schedule for {dih_key}: "
                        f"{[(round(_signed_deg(a), 1), t) for a, t in schedule]}"
                    )
                    for target_deg, mode in schedule:
                        mol_i = Molecule.from_xyz_string(seed_molecule.get_xyz_string())
                        mol_i.set_dihedral_in_degrees(dih_key, target_deg, verbose=False)
                        actual = float(mol_i.get_dihedral_in_degrees(dih_key)) % 360.0
                        if _circ_err_deg(actual, target_deg) > 2.0:
                            raise RuntimeError(
                                f"Failed to set dihedral {dih_key} to {target_deg:.2f} deg "
                                f"(actual {actual:.2f} deg)."
                            )
                        conformers_plus_ts[0][dih_key].append((mol_i, mode))
            else:
                if len(conformal_structures['molecules']) == 0:
                    raise RuntimeError('ConformerGenerator returned no conformers.')
                conformers_plus_ts[0][None] = [(conformal_structures['molecules'][0], 'normal')]
            # conformers_payload =  {}

            # for state, dih_map in conformers_plus_ts.items():
            #     conformers_payload[state] = {}
            #     for dih_key, mol_info in dih_map.items():
            #         key_paload = '__NONE__' if dih_key is None else tuple(dih_key)
            #         conformers_payload[state][key_paload] = [
            #             {
            #                 'xyz':mol_obj.get_xyz_string(),
            #                 'tag':tag,
            #             }
            #             for mol_obj, tag in mol_info
                    # ]
            # conformers_payload = comm.bcast(conformers_payload)
            # comm.barrier()
            

            # rebuilt = {}
            # for state, dih_map in conformers_payload.items():
            #     rebuilt[state] = {}
            #     for key_payload, entries in dih_map.items():
            #         dih_key = None if key_payload == '__NONE__' else tuple(key_payload)
            #         rebuilt_entries = []
            #         for item in entries:
            #             mol_obj = Molecule.from_xyz_string(item['xyz'])
            #             mol_obj.set_charge(molecule.get_charge())
            #             mol_obj.set_multiplicity(molecule.get_multiplicity())
            #             rebuilt_entries.append((mol_obj, item['tag']))
            #         rebuilt[state][dih_key] = rebuilt_entries

            self.seed_structures = conformers_plus_ts


    
    def _bootstrap_sampling_db_from_abinito_db(self, root):
        """
        creating a sample database of a cheap tight-binding method
        for an efficient sampling mode allowing to estimate when ab-inito
        reference calculation is required (only for ground-state for now!)
        """
        # comm = MPI.COMM_WORLD
        # rank = comm.Get_rank()
        # root_rank = mpi_master()

        # compare to the current interpolation database in order to get a 
        # identical mirroring of the ab-inito database -- note that this is
        # currently only for the ground-state
        ref_settings = self.states_interpolation_settings[root]
        sampling_settings = self.sampling_states_interpolation_settings[root]

        ref_file = ref_settings['imforcefield_file']
        sampling_file = sampling_settings['imforcefield_file']

        # ref_file_exists = os.path.exists(ref_file) if rank == root_rank else None
        # ref_file_exists = comm.bcast(ref_file_exists, root=root_rank)
        # if not ref_file_exists:
        #     comm.barrier()
        #     return
        
        ref_drv = InterpolationDriver(self.roots_z_matrix[root])
        ref_drv.update_settings(ref_settings)
        ref_db_labels, _ = ref_drv.read_labels()

        existing = set()

        # check if database already contains the given datapoint
        if os.path.exists(sampling_file):
            sampling_drv = InterpolationDriver(self.roots_z_matrix[root])
            sampling_drv.update_settings(sampling_settings)
            samp_db_labels, _ = sampling_drv.read_labels()
            existing = set(samp_db_labels)
        
        mol_labels = self.molecule.get_labels()
        sampling_qm, sampling_grad, sampling_hess = self.sampling_driver['gs'] #if rank == root_rank else None, None, None

        for label in ref_db_labels:
            print('current label', label, label in existing)
            if label in existing:
                continue

            ref_dp = InterpolationDatapoint(self.roots_z_matrix[root])
            ref_dp.update_settings(ref_settings)
            ref_dp.read_hdf5(ref_file, label)

            coords_ang = ref_dp.cartesian_coordinates * bohr_in_angstrom()
            mol = Molecule(mol_labels, coords_ang, 'angstrom')
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
            samp_dp.mapping_masks = getattr(ref_dp, 'mapping_masks', None)
            samp_dp.imp_int_coordinates = getattr(ref_dp, 'imp_int_coordinates', [])
            samp_dp.inv_sqrt_masses = inv_sqrt
            samp_dp.energy = e[0]
            samp_dp.gradient = mw_grad.reshape(g[0].shape)
            samp_dp.hessian = mw_hess.reshape(h[0].shape)
            samp_dp.confidence_radius = getattr(ref_dp, 'confidence_radius', self.confidence_radius)
            samp_dp.transform_gradient_and_hessian()

            samp_dp.family_label = ref_dp.family_label
            samp_dp.point_label = ref_dp.point_label
            samp_dp.bank_role = ref_dp.bank_role

            # if rank == root_rank:
            samp_dp.write_hdf5(sampling_file, label)

        # comm.barrier()

            
    def compute(self, molecule, states_basis=None, test_hooks=None):

        """
        Construct the interpolation dynamics database by generating molecular structures, 
        performing QM calculations, and collecting data points.

        :param molecule: The input molecular structure for which the database is constructed.
                        This molecule serves as a reference for sampling and simulation tasks.

        The method sets up quantum mechanical drivers, initializes interpolation and dynamics settings,
        and iterates over molecular structures along a predefined reaction path. It runs simulations
        to expand/generate the interpolation forcefield with new data points.

        """

        def _bridge(event_name):
            def _cb(payload):
                self._emit_test_hook(
                    f"qmmm.{event_name}",
                    {
                        "state": int(0),
                        "dihedral_key": key,
                        "structure_index": int(i),
                        "collector": payload,
                    },
                )
            return _cb

        # comm = MPI.COMM_WORLD
        # rank = comm.Get_rank()
        # root = mpi_master()
        
        self._test_hooks = test_hooks or {}
        self._test_hook_strict = bool(self._test_hooks.get("__strict__", False))
        self._test_event_counter = 0
        self._test_run_id = self._test_hooks.get("__run_id__", None)

        if self.reference_struc_energy_file is not None and not os.path.exists(self.reference_struc_energy_file):
            self.reference_struc_energy_file = None

        # First set up the system for which the database needs to be constructed
        states_basis = {'gs':self.gs_basis_set_label, 'es':self.es_basis_set_label}
        root_extract_z_matrix = None
        if self.imforcefieldfiles is not None:
            root_extract_z_matrix = {}
            self.sampling_imforcefieldfiles = {}
            for i, root in enumerate(self.roots_to_follow):
                if root not in self.imforcefieldfiles or not os.path.exists(self.imforcefieldfiles[root]):
                    self.imforcefieldfiles[self.roots_to_follow[i]] = f'im_database_{root}.h5'
                    root_extract_z_matrix[self.roots_to_follow[i]] = False
                else:
                    root_extract_z_matrix[self.roots_to_follow[i]] = True

                self.sampling_imforcefieldfiles[self.roots_to_follow[i]] = f'im_database_sampling_{root}.h5'
        else:
            self.imforcefieldfiles = {}
            self.sampling_imforcefieldfiles = {}
            standard_files = [f'im_database_{root}.h5' for root in self.roots_to_follow]
            standard_smapling_files = [f'im_database_sampling_{root}.h5' for root in self.roots_to_follow]
            root_extract_z_matrix = {}
            for root_idx, standard_file in enumerate(standard_files):
                
                if os.path.exists(standard_file):
                    if self.roots_to_follow[root_idx] not in root_extract_z_matrix:
                        root_extract_z_matrix[self.roots_to_follow[root_idx]] = True 
                else:
                    root_extract_z_matrix[self.roots_to_follow[root_idx]] = False     
                self.imforcefieldfiles[self.roots_to_follow[root_idx]] = standard_file
                self.sampling_imforcefieldfiles[self.roots_to_follow[root_idx]] = standard_smapling_files[root_idx]

        self._emit_test_hook("ffg.compute_start", {
            "imforcefieldfiles": dict(self.imforcefieldfiles or {}),
            "mode_signature": self._mode_signature(),
            "dynamics": {
                "ensemble": str(self.ensemble),
                "nsteps": int(self.nsteps),
                "timestep_fs": float(self.timestep),
                "temperature_K": float(self.temperature),
            },
        })

        print(f'IMPORTANT: IM ForceFieldFile is initalized from the current directory as {self.imforcefieldfiles}')
        self.set_up_the_system(molecule, extract_z_matrix=root_extract_z_matrix)

        self._emit_test_hook("ffg.setup_complete", {
            "z_matrix_summary": self.roots_z_matrix,
            "symmetry_information": self.symmetry_information,
            "all_rotatable_bonds": self.all_rotatable_bonds,
        })

        
        print('Set up the system is here')
        if len(self.roots_to_follow) > 1:
            print('The molecule that is used for the database construction is minimized always based on the ground state Energy -> to find the local minimum based on the '
            'potential energy surface of the ground state!')
            assert_msg_critical(all(self.roots_to_follow) == 0, "Only ground-state potential construction is currently supported! Will be updated to multi-state in the near future!")

                    
        else:
            print('The single-state collection considers only 1 state which can get problematic if the states are not sufficiently seperated!')
            assert_msg_critical(len(self.roots_to_follow) == 1, 'ImForceFieldGenerator: The root to follow is not defined!')

            imforcefieldfile = self.imforcefieldfiles[self.roots_to_follow[0]]
            self.states_interpolation_settings[self.roots_to_follow[0]] = { 'interpolation_type':self.interpolation_type,
                                'weightfunction_type':self.weightfunction_type, 
                                'exponent_p':self.exponent_p,
                                'exponent_q':self.exponent_q, 
                                'confidence_radius':self.confidence_radius,
                                'imforcefield_file':imforcefieldfile,
                                'use_inverse_bond_length':self.use_inverse_bond_length,
                                'use_eq_bond_length':self.use_eq_bond_length,
                                'eq_bond_symmetry_mode':self.eq_bond_symmetry_mode,
                                'use_cosine_dihedral':self.use_cosine_dihedral,
                                'use_tc_weights':self.use_tc_weights,
                                # 'use_mpi_preload': self.use_mpi_preload,
                                # 'use_local_preload': self.use_local_preload,
                                # 'local_preload_workers': self.local_preload_workers,
                                # 'local_preload_min_tasks': self.local_preload_min_tasks,
                                # 'local_preload_omp_threads': self.local_preload_omp_threads,
                                # 'use_outer_parallel': self.use_outer_parallel,
                                # 'outer_parallel_workers': self.outer_parallel_workers,
                                # 'outer_parallel_min_labels': self.outer_parallel_min_labels,
                                # 'outer_parallel_chunk_size': self.outer_parallel_chunk_size,
                            }
            self.sampling_states_interpolation_settings[self.roots_to_follow[0]] = self.states_interpolation_settings[self.roots_to_follow[0]].copy()
            self.sampling_states_interpolation_settings[self.roots_to_follow[0]]['imforcefield_file'] = self.sampling_imforcefieldfiles[self.roots_to_follow[0]]


            self.dynamics_settings = {  'drivers':self.drivers,
                                        'basis_set_label': states_basis,
                                        'duration':self.duration, 'temperature':self.temperature, 'solvent':self.solvent,
                                        'pressure':self.force_constant, 'force_constant': self.force_constant, 'ensemble':self.ensemble,
                                        'timestep': self.timestep, 'nsteps': self.nsteps, 'friction':self.friction,
                                        'snapshots':self.snapshots, 'trajectory_file':self.trajectory_file, 'reference_struc_energy_file':self.reference_struc_energy_file,
                                        'desired_datapoint_density':self.desired_point_density, 'converged_cycle': self.converged_cycle, 
                                        'energy_threshold':self.energy_threshold, 'grad_rmsd_thrsh': self.gradient_rmsd_thrsh,
                                        'load_system': None, 'start_collect': self.start_collect, 'roots_to_follow':self.roots_to_follow,
                                        'sampling_drivers': self.sampling_driver, 'sampling_settings':self.sampling_settings,
                                        'metadynamics':self.metadynamics_settings,
                                        # 'profile_runtime_timing': self.profile_runtime_timing, 'profile_interpolation_timing': self.profile_interpolation_timing,
                                        # 'mpi_control_plane_enabled': self.mpi_control_plane_enabled,
                                        # 'mpi_root_worker_mode': self.mpi_root_worker_mode,
                                        # 'mpi_reload_from_hdf5': self.mpi_reload_from_hdf5,
                                        # 'mpi_debug_sync': self.mpi_debug_sync
                                        }
            
            files_to_add_conf = []
            molecules_to_add_info = []             

            if not os.path.exists(self.imforcefieldfiles[self.roots_to_follow[0]]):
                    files_to_add_conf.append(self.roots_to_follow[0])


            if self.seed_structures is not None and not os.path.exists(imforcefieldfile):
    
                molecules_to_add_info = []
                for counter, entry in enumerate(self.seed_structures.items()):
       
                    key, molecules_info = entry
                    
                    for i, mol_entries in enumerate(molecules_info.items()):
                        current_molecule_to_add_info = []
                        dih_key, mol_info = mol_entries

                        for seed_entry in mol_info:
                            if len(seed_entry) == 2:
                                mol, mode = seed_entry
                                seed_constraints = []
                            else:
                                mol, mode, seed_constraints = seed_entry[0], seed_entry[1], list(seed_entry[2] or [])
          
                            if self.use_minimized_structures[0]:
                                transition = False
                                constraints_global = []
                                if len(self.use_minimized_structures[1]) > 0:
                                    constraints_global.extend(self.use_minimized_structures[1])  # extend, not append

                                if mode == "transition":
                                    transition = True
                                    constraints_global = []

                                if mode == "constraint" and dih_key is not None:
                                    if dih_key not in constraints_global:
                                        constraints_global.append(dih_key)

                                for c in seed_constraints:
                                    if c not in constraints_global:
                                        constraints_global.append(c)

                                if self.use_minimized_structures[0]:       
                                    optimized_molecule = None

                                    opt_results = None
                                    scf_results = None

                                    if self.roots_to_follow[0] == 0 and isinstance(self.drivers['gs'][0], ScfRestrictedDriver) or self.roots_to_follow[0] == 0 and isinstance(self.drivers['gs'][0], ScfUnrestrictedDriver):
                                    

                                        current_basis = MolecularBasis.read(mol, states_basis['gs'])
                                        _, scf_results, _ = self._compute_energy(self.drivers['gs'][0], mol, current_basis)

                                        # scf_results_mpi = self._compute_energy_mpi_safe(
                                        #     self.drivers['gs'][0],
                                        #     mol,
                                        #     current_basis,
                                        #     collective=True,
                                        #     phase_name='Energy calcualtion',
                                        # )
                                        # scf_results = scf_results_mpi[1]
                                        # _, scf_results, _ = self._compute_energy(self.drivers['gs'][0], mol, current_basis)
                                        print('original molecule', mol.get_xyz_string(), constraints_global, transition)
                                        
                                        opt_results = self._run_optimization(
                                            self.drivers['gs'][0],
                                            mol,
                                            constraints=constraints_global,
                                            transition=transition,
                                            index_offset=0,
                                            compute_args=(current_basis, scf_results),
                                            # source_molecule=molecule,
                                        )
                                        # opt_results_mpi = self._run_optimization_mpi_safe(
                                        #     self.drivers['gs'][0],
                                        #     mol,
                                        #     constraints=constraints_global,
                                        #     transition=transition,
                                        #     index_offset=0,
                                        #     compute_args=(current_basis, scf_results),
                                        #     source_molecule=molecule,
                                        #     collective=True,
                                        #     phase_name='optimization tag',
                                        # )
                                        
                                        # opt_results = opt_results_mpi[1]

                                        optimized_molecule = opt_results['final_molecule']

                                        current_basis = MolecularBasis.read(optimized_molecule, states_basis['gs'])
                                        current_molecule_to_add_info.append((optimized_molecule, current_basis, self.roots_to_follow, constraints_global, transition))

                                        print(optimized_molecule.get_xyz_string())
                                  

                                    elif  self.roots_to_follow[0] == 0 and isinstance(self.drivers['gs'][0], XtbDriver):
                                        
                                        opt_results = self._run_optimization(
                                            self.drivers['gs'][0],
                                            mol,
                                            constraints=constraints_global,
                                            transition=transition,
                                            index_offset=0,
                                        )

                                        # opt_results_mpi = self._run_optimization_mpi_safe(
                                        #     self.drivers['gs'][0],
                                        #     mol,
                                        #     constraints=constraints_global,
                                        #     transition=transition,
                                        #     index_offset=0,
                                        #     collective=True,
                                        #     phase_name='optimization tag',
                                        # )
                                        
                                        # opt_results = opt_results_mpi[1]
                                        optimized_molecule = opt_results['final_molecule']

                                        print('Optimization enegiers', opt_results['opt_energies'][-1] * hartree_in_kcalpermol())

                                        current_basis = MolecularBasis.read(optimized_molecule, states_basis['gs'])
                                        current_molecule_to_add_info.append((optimized_molecule, current_basis, self.roots_to_follow, constraints_global, transition))

                                        print('Optimized mol', optimized_molecule.get_xyz_string(), transition)

                                    elif self.roots_to_follow[0] >= 1 and isinstance(self.drivers['es'][0], LinearResponseEigenSolver) or self.roots_to_follow[0] >= 1 and isinstance(self.drivers['es'][0], TdaEigenSolver):

                                            current_basis = MolecularBasis.read(mol, states_basis['es'])
                                            _, _, rsp_results  = self._compute_energy_mpi_safe(
                                                self.drivers['es'][0],
                                                mol,
                                                current_basis,
                                                collective=True,
                                                phase_name='Energy calcualtion',
                                            )
                                            self.drivers['es'][3].ostream.mute()
                                            self.drivers['es'][0].ostream.mute()
                                            opt_results_mpi = self._run_optimization_mpi_safe(
                                                self.drivers['es'][1],
                                                molecule,
                                                constraints=self.use_minimized_structures[1],
                                                index_offset=1,
                                                compute_args=(current_basis, self.drivers['es'][3], self.drivers['es'][0], rsp_results),
                                                source_molecule=molecule,
                                                collective=True,
                                                phase_name='optimization tag',
                                            )

                                            opt_results = opt_results_mpi[1]
                                            optimized_molecule = opt_results['final_molecule']
                                            
                                            print('Optimized Molecule', optimized_molecule.get_xyz_string(), '\n\n', molecule.get_xyz_string())
                                            current_basis = MolecularBasis.read(optimized_molecule, states_basis['es'])
                                            current_molecule_to_add_info.append((optimized_molecule, current_basis, self.roots_to_follow, constraints_global))

                            else:
                                if self.roots_to_follow[0] == 0:
                                    current_basis = MolecularBasis.read(mol, states_basis['gs'])
                                    current_molecule_to_add_info.append((mol, current_basis, self.roots_to_follow, []))
                                else:
                                    current_basis = MolecularBasis.read(mol, states_basis['es'])
                                    current_molecule_to_add_info.append((mol, current_basis, self.roots_to_follow, []))
                        if len(molecules_to_add_info) == 0:

                            molecules_to_add_info.append(current_molecule_to_add_info[0])

                        self.add_point(current_molecule_to_add_info, self.states_interpolation_settings, symmetry_information=self.symmetry_information)
  
            elif not self.add_conformal_structures and not os.path.exists(imforcefieldfile):
        
                molecules_to_add_info = []
                if self.use_minimized_structures[0]:       
                    optimized_molecule = None

                    opt_results = None
                    scf_results = None


                    if self.roots_to_follow[0] == 0 and isinstance(self.drivers['gs'][0], ScfRestrictedDriver) or self.roots_to_follow[0] == 0 and isinstance(self.drivers['gs'][0], ScfUnrestrictedDriver):
                    
                        current_basis = MolecularBasis.read(molecule, states_basis['gs'])
                        _,scf_results, _ = self._compute_energy(self.drivers['gs'][0], molecule, current_basis)
                        # scf_results_mpi = self._compute_energy_mpi_safe(
                        #     self.drivers['gs'][0],
                        #     molecule,
                        #     current_basis,
                        #     collective=True,
                        #     phase_name='Energy calcualtion',
                        # )
                        # scf_results = scf_results_mpi[1]

                        opt_results = self._run_optimization(
                                            self.drivers['gs'][0],
                                            molecule,
                                            constraints=self.use_minimized_structures[1],
                                            index_offset=0,
                                            compute_args=(current_basis, scf_results),
                                            # source_molecule=molecule,
                                            )

                        # opt_results_mpi = self._run_optimization_mpi_safe(
                        #     self.drivers['gs'][0],
                        #     molecule,
                        #     constraints=self.use_minimized_structures[1],
                        #     index_offset=0,
                        #     compute_args=(current_basis, scf_results),
                        #     collective=True,
                        #     phase_name='optimization tag',
                        # )
                        
                        # opt_results = opt_results_mpi[1]
                        optimized_molecule = opt_results['final_molecule']

                        current_basis = MolecularBasis.read(optimized_molecule, states_basis['gs'])
                        molecules_to_add_info.append((optimized_molecule, current_basis, self.roots_to_follow, self.use_minimized_structures[1]))

                        print('Optimized Molecule', optimized_molecule.get_xyz_string())
                    
                    elif  self.roots_to_follow[0] == 0 and isinstance(self.drivers['gs'][0], XtbDriver):
                        
                        opt_results = self._run_optimization(
                                            self.drivers['gs'][0],
                                            molecule,
                                            constraints=self.use_minimized_structures[1],
                                            index_offset=0,
                                        )
                        
                        # opt_results_mpi = self._run_optimization_mpi_safe(
                        #     self.drivers['gs'][0],
                        #     molecule,
                        #     constraints=self.use_minimized_structures[1],
                        #     index_offset=0,
                        #     collective=True,
                        #     phase_name='optimization tag',
                        # )
                        
                        # opt_results = opt_results_mpi[1]
                        optimized_molecule = opt_results['final_molecule']


                        current_basis = MolecularBasis.read(optimized_molecule, states_basis['gs'])
                        molecules_to_add_info.append((optimized_molecule, current_basis, self.roots_to_follow, self.use_minimized_structures[1]))

                        print('optimized molecule', optimized_molecule.get_xyz_string())


                    elif self.roots_to_follow[0] >= 1 and isinstance(self.drivers['es'][0], LinearResponseEigenSolver) or self.roots_to_follow[0] >= 1 and isinstance(self.drivers['es'][0], TdaEigenSolver):

                            current_basis = MolecularBasis.read(molecule, states_basis['es'])
                            _, _, rsp_results  = self._compute_energy_mpi_safe(
                                self.drivers['es'][0],
                                molecule,
                                current_basis,
                                collective=True,
                                phase_name='Energy calcualtion',
                            )
                            self.drivers['es'][3].ostream.mute()
                            self.drivers['es'][0].ostream.mute()
                            opt_results_mpi = self._run_optimization_mpi_safe(
                                self.drivers['es'][1],
                                molecule,
                                constraints=self.use_minimized_structures[1],
                                index_offset=0,
                                compute_args=(current_basis, self.drivers['es'][3], self.drivers['es'][0], rsp_results),
                                collective=True,
                                phase_name='optimization tag',
                            )
                            opt_results = opt_results_mpi[1]
                            optimized_molecule = opt_results['final_molecule']

                            print('Optimized Molecule', optimized_molecule.get_xyz_string(), '\n\n', molecule.get_xyz_string())
                            current_basis = MolecularBasis.read(optimized_molecule, states_basis['es'])
                            molecules_to_add_info.append((optimized_molecule, current_basis, self.roots_to_follow, self.use_minimized_structures[1]))
                    
                    self.add_point(molecules_to_add_info, self.states_interpolation_settings, symmetry_information=self.symmetry_information)

                else:
                    if self.roots_to_follow[0] == 0:
                        current_basis = MolecularBasis.read(molecule, states_basis['gs'])
                        molecules_to_add_info.append((molecule, current_basis, self.roots_to_follow, []))
                    else:
                        current_basis = MolecularBasis.read(molecule, states_basis['es'])
                        molecules_to_add_info.append((molecule, current_basis, self.roots_to_follow, []))
                    
                    self.add_point(molecules_to_add_info, self.states_interpolation_settings, symmetry_information=self.symmetry_information)
            
            else:
                
                if self.roots_to_follow[0] == 0:
                    current_basis = MolecularBasis.read(molecule, states_basis['gs'])
                    molecules_to_add_info.append((molecule, current_basis, self.roots_to_follow, []))
                else:
                    current_basis = MolecularBasis.read(molecule, states_basis['es'])
                    molecules_to_add_info.append((molecule, current_basis, self.roots_to_follow, []))

                if not os.path.exists(imforcefieldfile):
                    self.add_point(molecules_to_add_info, self.states_interpolation_settings, symmetry_information=self.symmetry_information)
            
            density_of_datapoints = self.determine_datapoint_density(self.states_interpolation_settings)

            self.states_data_point_density = density_of_datapoints
            
            if self.sampling_settings.get('enabled', False):
                self._bootstrap_sampling_db_from_abinito_db(self.roots_to_follow[0])
    
            # for counter, (state, dihedral_dict) in enumerate(self.molecules_along_rp.items()):

                # for key, mol_info in dihedral_dict.items():
            dynamics_molecule = molecules_to_add_info[self.roots_to_follow[0]][0]
            forcefield_generator = MMForceFieldGenerator()
            self.dynamics_settings['trajectory_file'] = f'trajectory_{self.roots_to_follow[0]}.pdb'
            forcefield_generator.partial_charges = dynamics_molecule.get_partial_charges(dynamics_molecule.get_charge())
            
            forcefield_generator.create_topology(dynamics_molecule)
            im_database_driver = IMDatabasePointCollecter()
            im_database_driver.distance_thrsh = self.distance_thrsh
            im_database_driver.non_core_symmetry_groups = self.symmetry_information
            im_database_driver.platform = self.open_mm_platform
            im_database_driver.all_rot_bonds = self.all_rotatable_bonds
            
            # set optimization features in the construction run
            im_database_driver.identfy_relevant_int_coordinates = (self.identfy_relevant_int_coordinates, self.use_minimized_structures[1])
            im_database_driver.use_symmetry = self.use_symmetry
            im_database_driver.cluster_run = self.cluster_run
            im_database_driver.symmetry_rotors = self.symmetry_rotors
            im_database_driver.rotor_corr_threshold = self.rotor_corr_threshold
            im_database_driver.use_opt_confidence_radius = self.use_opt_confidence_radius
            
            
            im_database_driver.system_from_molecule(dynamics_molecule, self.roots_z_matrix, forcefield_generator, solvent=self.solvent, qm_atoms='all')  
            if self.bias_force_reaction_prop is not None:
                im_database_driver.bias_force_reaction_idx = self.bias_force_reaction_idx
                im_database_driver.bias_force_reaction_prop = self.bias_force_reaction_prop
            
            density_of_datapoints = self.determine_datapoint_density(self.states_interpolation_settings)
            self.density_of_datapoints = density_of_datapoints
            print('density of points', self.states_data_point_density)
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
                # if key is None:
                #     im_database_driver.allowed_molecule_deviation = self.allowed_deviation
                # else:
                #     im_database_driver.allowed_molecule_deviation = self.allowed_deviation
                im_database_driver.update_settings(self.dynamics_settings, self.states_interpolation_settings, self.sampling_states_interpolation_settings)
                
                self._emit_test_hook("ffg.run_qmmm_start", {
                    "state": int(self.roots_to_follow[0]),
                    "structure_index": int(0),
                    "current_structure_density": current_structure_density,
                })
                im_database_driver.run_qmmm(
                    test_hooks={
                        "run_start": _bridge("run_start"),
                        "step": _bridge("step"),
                        "point_correlation_decision": _bridge("point_correlation_decision"),
                        "datapoint_written": _bridge("datapoint_written"),
                        "alpha_optimization_end": _bridge("alpha_optimization_end"),
                        "run_end": _bridge("run_end"),
                    },
                    collect_step_trace=False,
                    strict_test_hooks=self._test_hook_strict,
                )
                self._emit_test_hook("ffg.run_qmmm_end", {
                    "state": int(self.roots_to_follow[0]),
                    "structure_index": int(0),
                })
 
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
            
            print('The construction of the database was sucessfull', self.states_data_point_density)
            self.im_results['n_datapoints'] = self.states_data_point_density

            self._emit_test_hook("ffg.compute_end", {
                "states_data_point_density": self.states_data_point_density,
                "im_results": self.im_results,
                "final_label_count_per_root": self._label_count_per_root(),
            })
            
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
        def dihedral_to_vector(angle):
            """
            Converts a dihedral angle in degrees to a 2D vector on the unit circle.
            :param angle:
                Angle to be transformed into a sinus and cosinus basis.
            
            :returns:
                Angle in vector form.
            """
            rad = np.radians(angle)
            return np.array([np.cos(rad), np.sin(rad)])
        
        def structure_to_vector(dihedrals):
            """
            Converts a list of dihedral angles to a concatenated vector.
            For N dihedrals, returns a vector of length 2N.
            :param dihedrals:
                List of dihedrals which will be transformed into vector form.
            
            :return:
                A concetenate list of the vectors.
            """
            return np.concatenate([dihedral_to_vector(angle) for angle in dihedrals])
        
        def dihedral_distance_vectorized(dihedrals1, dihedrals2):
            """
            Computes the Euclidean distance between two sets of dihedrals by mapping each
            angle to a unit circle vector.
            
            :param dihedrals1:
                Lists or tuples of angles (in degrees).
            
            :param dihedrals2:
                Lists or tuples of angles (in degrees).
            
            :returns:
                Norm of the distance.
            """
            vec1 = structure_to_vector(dihedrals1)
            vec2 = structure_to_vector(dihedrals2)
            return np.linalg.norm(vec1 - vec2)
        
        # reseted_point_densities_dict = {outer_key: {inner_key: {key: 0 for key in point_densities_dict[outer_key][inner_key].keys()} for inner_key in point_densities_dict[outer_key].keys()} for outer_key in point_densities_dict.keys()}
        reseted_point_densities_dict = {state: 0 for state in self.roots_to_follow}

        for state in reseted_point_densities_dict.keys():
            qm_datapoints = []
            if imforcefieldfile[state]['imforcefield_file'] in os.listdir(os.getcwd()):
                impes_driver = InterpolationDriver(self.roots_z_matrix[state])
                impes_driver.update_settings(imforcefieldfile[state])
                self.qmlabels, z_matrix = impes_driver.read_labels()
               
             
                for label in self.qmlabels:
                    qm_data_point = InterpolationDatapoint(z_matrix)
                    qm_data_point.update_settings(imforcefieldfile[state])
                    qm_data_point.read_hdf5(imforcefieldfile[state]['imforcefield_file'], label)
                    if qm_data_point.bank_role == "core":
                        qm_datapoints.append(qm_data_point)
                        reseted_point_densities_dict[state] += 1
            
            # for specific_dihedral in point_densities_dict[state].keys():
            #     for point in qm_datapoints:
            #         if specific_dihedral is None:
            #             reseted_point_densities_dict[state][specific_dihedral][360] += 1
            #         else:
            #             min_distance = np.inf
            #             key = None
            #             for dihedral in point_densities_dict[state][specific_dihedral].keys():
            #                 datapoint_molecule = Molecule(self.molecule.get_labels(), point.cartesian_coordinates, 'bohr')
            #                 dihedrals_of_dp = [datapoint_molecule.get_dihedral_in_degrees([specific_dihedral[0], specific_dihedral[1], specific_dihedral[2], specific_dihedral[3]])]
            #                 distance_vectorized = dihedral_distance_vectorized([dihedral], dihedrals_of_dp)
            #                 if abs(distance_vectorized) < min_distance:
            #                     min_distance = abs(distance_vectorized)
            #                     key = dihedral
                    
            #             reseted_point_densities_dict[state][specific_dihedral][key] += 1

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
        rotation_matrix_core = geometric.rotate.get_rot(target_coordinates,
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

        im_driver = InterpolationDriver() # -> implemented Class in VeloxChem that is capable to perform interpolation calculations for a given molecule and provided z_matrix and database
        im_driver.update_settings(im_settings)
        # im_driver.imforcefield_file = datafile
        labels, z_matrix = im_driver.read_labels()
    
        sorted_labels = sorted(labels, key=lambda x: int(x.split('_')[1]))

        # impes_coordinate = InterpolationDatapoint(z_matrix) # -> implemented Class in VeloxChem that handles all transformations and database changes concerning the interpolation
        data_point_molecules = []
        datapoints = []

        for label in sorted_labels:
            impes_coordinate = InterpolationDatapoint(z_matrix)
            impes_coordinate.update_settings(im_settings)
            impes_coordinate.read_hdf5(datafile, label) # -> read in function from the ImpesDriver object

            coordinates_in_angstrom = impes_coordinate.cartesian_coordinates * bohr_in_angstrom()

           
            current_molecule = Molecule(mol_labels, coordinates_in_angstrom, 'angstrom') # -> creates a VeloxChem Molecule object

            datapoints.append(impes_coordinate)
            data_point_molecules.append(current_molecule)

        return data_point_molecules, datapoints, z_matrix

    def _list_rotor_cluster_families_from_registry(self, root, im_settings):
        imff_file = im_settings["imforcefield_file"]
        prefix = f"rotor_cluster_registry/root_{root}"
        with h5py.File(imff_file, "r") as h5f:
            if prefix not in h5f:
                return []
            return sorted(name.replace("family_", "", 1) for name in h5f[prefix].keys())


    def _read_cluster_registry_for_family(self, imff_file, root, family_label):
        def _read_scalar_string(ds):
            val = ds[()]
            return val.decode("utf-8") if isinstance(val, bytes) else str(val)

        with h5py.File(imff_file, "r") as h5f:
            prefix = f"rotor_cluster_registry/root_{root}/family_{family_label}"
            info_json = _read_scalar_string(h5f[prefix + "/cluster_info_json"])
            library_json = _read_scalar_string(h5f[prefix + "/cluster_angle_library_json"])
            index_json = _read_scalar_string(h5f[prefix + "/point_index_json"])

        info_payload = json.loads(info_json)
        library_payload = json.loads(library_json)
        index_payload = json.loads(index_json)
        
        rotor_to_cluster_map = {
            int(rid): int(cid) for rid, cid in info_payload["rotor_to_cluster"].items()
        }
        cluster_info = RotorClusterInformation(
            dihedral_start=int(info_payload["dihedral_start"]),
            dihedral_end=int(info_payload["dihedral_end"]),
            rotor_to_cluster=rotor_to_cluster_map,
        )
        for rotor_id, rotor_data in info_payload["rotors"].items():
            cluster_info.rotors[int(rotor_id)] = RotorDefinition(
                rotor_id=int(rotor_id),
                center=tuple(rotor_data["center"]),
                torsion_rows=tuple(rotor_data["torsion_rows"]),
                torsion_coords=tuple(tuple(x) for x in rotor_data["torsion_coords"]),
                symmetry_order=int(rotor_data["symmetry_order"]),
                atom_group=tuple(rotor_data["atom_group"]),
            )
        for cid, cdata in info_payload["clusters"].items():
            cluster_info.clusters[int(cid)] = RotorClusterDefinition(
                cluster_id=int(cid),
                rotor_ids=tuple(cdata["rotor_ids"]),
                cluster_type=cdata["cluster_type"],
                torsion_rows=tuple(cdata["torsion_rows"]),
            )

        angle_library = RotorClusterAngleLibrary()
        for cid, state_list in library_payload.items():
            bank = []
            for state in state_list:
                aa = {}
                for key, val in state["angle_assignment"].items():
                    aa[tuple(int(x) for x in key.split(","))] = float(val)
                bank.append(
                    RotorClusterStateDefinition(
                        cluster_id=int(cid),
                        state_id=int(state["state_id"]),
                        cluster_type=state["cluster_type"],
                        rotor_ids=tuple(state["rotor_ids"]),
                        angle_assignment=aa,
                        dihedrals_to_rotate=(
                            None
                            if state["dihedrals_to_rotate"] is None
                            else tuple(tuple(int(x) for x in row) for row in state["dihedrals_to_rotate"])
                        ),
                        phase_signature=(
                            None if state["phase_signature"] is None
                            else np.asarray(state["phase_signature"], dtype=np.float64)
                        ),
                        is_anchor=bool(state.get("is_anchor", False)),
                        label_suffix=str(state.get("label_suffix", "")),
                    )
                )
            angle_library.state_banks[int(cid)] = tuple(bank)

        point_index = {
            "family_label": index_payload["family_label"],
            "core_label": index_payload["core_label"],
            "cluster_state_labels": {
                int(cid): {int(sid): lbl for sid, lbl in smap.items()}
                for cid, smap in index_payload["cluster_state_labels"].items()
            },
        }
        return cluster_info, angle_library, point_index
    
    def _load_rotor_cluster_bank_for_root(self, root, im_settings):
        out = {}
        imff_file = im_settings["imforcefield_file"]
        families = self._list_rotor_cluster_families_from_registry(root, im_settings)

        for family in families:
            cluster_info, angle_library, point_index = self._read_cluster_registry_for_family(
                imff_file, root, family
            )
            fam = {
                "cluster_info": cluster_info,
                "cluster_angle_library": angle_library,
                "point_index": point_index,
                "core": None,
                "clusters": {
                    cid: {
                        "cluster_type": cluster.cluster_type,
                        "rotor_ids": tuple(cluster.rotor_ids),
                        "expected_states": {state.state_id: None for state in angle_library.state_banks[cid]},
                    }
                    for cid, cluster in cluster_info.clusters.items()
                },
            }

            core_dp = InterpolationDatapoint(self.roots_z_matrix[root])
            core_dp.update_settings(im_settings)
            core_dp.read_hdf5(imff_file, point_index["core_label"])
            fam["core"] = core_dp

            for cid, state_map in point_index["cluster_state_labels"].items():
                for sid, label in state_map.items():
                    dp = InterpolationDatapoint(self.roots_z_matrix[root])
                    dp.update_settings(im_settings)
                    dp.read_hdf5(imff_file, label)
                    fam["clusters"][cid]["expected_states"][sid] = dp

            for cid, cbank in fam["clusters"].items():
                expected = cbank["expected_states"]
                if 0 in expected and expected[0] is None:
                    expected[0] = fam["core"]

            out[family] = fam
        return out
    
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
                rotation_matrix = geometric.rotate.get_rot(target_coordinates,
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

        def write_distances_h5(h5_path: str, group_name: str, db_mol: list, N: int, pairs: np.ndarray, dists: np.ndarray, labels=None):
            """
            Store an edge list + metadata + XYZ structures in an HDF5 file.
            """
            with h5py.File(h5_path, "a") as f:
                # 1. Overwrite group if it exists
                if group_name in f:
                    del f[group_name]
                g = f.create_group(group_name)

                g.attrs["n_points"] = int(N)

                # 2. Extract XYZ strings from the molecule objects
                #    Assuming db_mol is a list of objects that have .get_xyz_string()
                xyz_list = [mol.get_xyz_string() for mol in db_mol]

                # 3. Define the variable-length string data type
                dt = h5py.string_dtype(encoding='utf-8')

                print('file 1', f[group_name])

                # 4. Write the XYZ dataset
                #    We use 'dtype=object' for numpy to handle variable length strings before passing to h5py
                g.create_dataset(
                    "xyz", 
                    data=np.asarray(xyz_list, dtype=object), 
                    dtype=dt, 
                    compression="gzip", 
                    compression_opts=4
                )

                print('file', f[group_name])

                # 5. Optional Labels
                if labels is not None:
                    labels = [str(x) for x in labels]
                    g.create_dataset("labels", data=np.asarray(labels, dtype=object), dtype=dt)

                if int(N) > 1:
                # 6. Store Distance Data
                    g.create_dataset("pairs", data=pairs, compression="gzip", compression_opts=4, shuffle=True)
                    g.create_dataset("dists", data=dists, compression="gzip", compression_opts=4, shuffle=True)
                
        overall_db_covergage = {}

        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        root_rank = mpi_master()

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
                if rank == root_rank:
                    all_structures_root = list(given_molecular_strucutres.get(root, []))
                    # payload = {
                    #     'skip_root': False,
                    #     'current_datafile': current_datafile,
                    #     'random_struct_info': [],
                    # }

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
                                print(
                                    f'The overall RMSD is {rmsd} -> '
                                    'The current structures are well seperated from the database conformations! '
                                    'loop is discontinued'
                                )
                                dist_ok = True
                            else:
                                print(
                                    f'The overall RMSD is {rmsd} -> '
                                    'The current structures are not all well seperated from the database conformations! '
                                    'loop is continued'
                                )

                        random_structure_choices_root = {"random_struct_info":list(selected_molecules)}
                #         payload['random_struct_info'] = [
                #             self._serialize_molecule_for_mpi(mol)
                #             for mol in random_structure_choices_root
                #         ]

                # payload = comm.bcast(payload, root=root_rank)
                # comm.barrier()

                if random_structure_choices_root is None:
                    raise RuntimeError('confirm_database_quality: failed to broadcast root payload.')

                # if payload.get('skip_root', False):
                #     database_quality = True
                #     continue

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

                old_label = None
                im_labels, _ = impes_driver.read_labels()
                impes_driver.qm_data_points = []
            
                for label in im_labels:
                    qm_data_point = InterpolationDatapoint(self.roots_z_matrix[root])
                    qm_data_point.update_settings(im_settings[root])
                    qm_data_point.read_hdf5(current_datafile, label)
                    if qm_data_point.bank_role == "core" or "cluster" not in qm_data_point.point_label and "symmetry" not in qm_data_point.point_label:
                        old_label = qm_data_point.point_label
                        impes_driver.qm_symmetry_data_points[old_label] = [qm_data_point]
                        impes_driver.qm_data_points.append(qm_data_point)
                        if impes_driver.impes_coordinate.eq_bond_lengths is None:
                            impes_driver.impes_coordinate.eq_bond_lengths = qm_data_point.eq_bond_lengths
                    elif qm_data_point.bank_role == "symmetry":
                        impes_driver.qm_symmetry_data_points[old_label].append(qm_data_point)
                
                
                qm_rotor_cluster_banks = self._load_rotor_cluster_bank_for_root(root, im_settings[root])
                impes_driver.qm_rotor_cluster_banks = qm_rotor_cluster_banks

                if impes_driver.qm_rotor_cluster_banks:
                    first_family = next(iter(impes_driver.qm_rotor_cluster_banks.values()))
                    impes_driver.rotor_cluster_information = first_family.get("cluster_info")
                else:
                    impes_driver.rotor_cluster_information = None
                # struct_to_check = [
                #     self._deserialize_molecule_from_mpi(mol_info)
                #     for mol_info in random_structure_choices_root.get('random_struct_info', [])
                # ]

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
                    

                    # scf_results_mpi = self._compute_energy_mpi_safe(
                    #     drivers[0],
                    #     mol,
                    #     current_basis,
                    #     collective=True,
                    #     phase_name='Energy calcualtion',
                    # )
                    # reference_energies, scf_results, rsp_results = (
                    #     scf_results_mpi[0],
                    #     scf_results_mpi[1],
                    #     scf_results_mpi[2],
                    # )

                    qm_energy_val = float(reference_energies[root]) if len(reference_energies) > root else float(reference_energies[0])
                    qm_energies.append(qm_energy_val)
                    im_energies.append(impes_driver.impes_coordinate.energy)
                    print(mol.get_xyz_string())

                    print(f'\n\n ########## Step {i} ######### \n')
                    print(
                        f'delta_E:  {abs(qm_energies[-1] - im_energies[-1]) * hartree_in_kcalpermol()} kcal/mol ::: '
                        f'{abs(qm_energies[-1] - im_energies[-1]) * hartree_in_kcalpermol() / len(molecule.get_labels())} '
                        'kcal/mol per atom  \n'
                    )

                    if (
                        (abs(qm_energies[-1] - im_energies[-1]) / len(molecule.get_labels()))
                        * hartree_in_kcalpermol()
                        > self.energy_threshold
                        and improve
                    ):
                        print(f"The current structure is not within the desired threshold and can be added to the database if desired! \n The structure can be found in the random_structure.xyz file")
  
                    else:
                        # if rank == root_rank:
                        with h5py.File('summary_output.h5', "a") as h5f:
                            self._append_confirm_database_quality_h5(
                                h5f,
                                [mol],
                                np.array([qm_energy_val]),
                                np.array([impes_driver.impes_coordinate.energy]),
                                root,
                            )
                    # comm.barrier()

                # database_expanded = bool(comm.allreduce(bool(database_expanded), op=MPI.LOR))
                last_qm_energies = qm_energies
                last_im_energies = im_energies

                if not database_expanded:
                    database_quality = True
                    # comm.barrier()
                    # if rank == root_rank:
                    overall_db_covergage[current_datafile] = {}
                    minmum_distances_in_db = database_distance_check(datapoint_molecules)
                    overall_db_covergage[current_datafile]['db_distances'] = minmum_distances_in_db
                    N, pairs, dists = dist_dict_to_edges(minmum_distances_in_db)
                    labels = [f"point {i+1}" for i in range(N)]
                    write_distances_h5(
                        current_datafile,
                        group_name="internal_dist_MDS",
                        db_mol=datapoint_molecules,
                        N=N,
                        pairs=pairs,
                        dists=dists,
                        labels=labels,
                    )
                # comm.barrier()
            # if rank == root_rank:
            print('Overall coverage', overall_db_covergage)
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

    def _write_cluster_registry_for_family(
        self,
        imff_file,
        root,
        family_label,
        cluster_info,
        cluster_angle_library,
        point_index,
    ):
        payload_info = {
            "dihedral_start": cluster_info.dihedral_start,
            "dihedral_end": cluster_info.dihedral_end,
            "rotors": {
                int(rotor_id): {
                    "center": list(rotor.center),
                    "torsion_rows": list(rotor.torsion_rows),
                    "torsion_coords": [list(t) for t in rotor.torsion_coords],
                    "symmetry_order": rotor.symmetry_order,
                    "atom_group": list(rotor.atom_group),
                }
                for rotor_id, rotor in cluster_info.rotors.items()
            },
            "clusters": {
                int(cluster_id): {
                    "rotor_ids": list(cluster.rotor_ids),
                    "cluster_type": cluster.cluster_type,
                    "torsion_rows": list(cluster.torsion_rows),
                }
                for cluster_id, cluster in cluster_info.clusters.items()
            },
            "rotor_to_cluster": {
                int(rotor_id): int(cluster_id)
                for rotor_id, cluster_id in cluster_info.rotor_to_cluster.items()
            },
        }

        payload_library = {
            int(cluster_id): [
                {
                    "state_id": state.state_id,
                    "cluster_type": state.cluster_type,
                    "rotor_ids": list(state.rotor_ids),
                    "angle_assignment": {
                        ",".join(str(x) for x in dihedral): float(angle)
                        for dihedral, angle in state.angle_assignment.items()
                    },
                    "dihedrals_to_rotate": [list(x) for x in (state.dihedrals_to_rotate or ())],
                    "phase_signature": (
                        None if state.phase_signature is None
                        else state.phase_signature
                    ),
                    "is_anchor": state.is_anchor,
                    "label_suffix": state.label_suffix,
                }
                for state in state_bank
            ]
            for cluster_id, state_bank in cluster_angle_library.state_banks.items()
        }

        payload_index = {
            "family_label": family_label,
            "core_label": point_index["core_label"],
            "cluster_state_labels": {
                int(cluster_id): {
                    int(state_id): label
                    for state_id, label in state_map.items()
                }
                for cluster_id, state_map in point_index["cluster_state_labels"].items()
            },
        }

        with h5py.File(imff_file, "a") as h5f:
            prefix = f"rotor_cluster_registry/root_{root}/family_{family_label}"
            self._write_string_dataset(
                h5f,
                prefix + "/cluster_info_json",
                json.dumps(payload_info, sort_keys=True),
            )
            self._write_string_dataset(
                h5f,
                prefix + "/cluster_angle_library_json",
                json.dumps(payload_library, sort_keys=True),
            )
            self._write_string_dataset(
                h5f,
                prefix + "/point_index_json",
                json.dumps(payload_index, sort_keys=True),
            )
    
    def _finalize_mapping_masks_and_eq_mode(self, impes_coordinate, masks):
        impes_coordinate.mapping_masks = masks

        if not impes_coordinate.use_eq_bond_length:
            return

        mode = str(getattr(impes_coordinate, 'eq_bond_symmetry_mode', 'masked_exact')).strip().lower()
        if mode == 'symmetrized':
            print('It is the symmetrized mode')
            impes_coordinate.symmetrize_eq_bond_lengths_from_masks(masks)
            impes_coordinate.transform_gradient_and_hessian()


    def add_point_rotor(self, molecule_specific_information, interpolation_settings, symmetry_information={}):
        """ Adds a new point to the database.

            :param molecule:
                the molecule.
            :param imforcefielddatafile:
                Datafile containing the information of the IM forcefield.

        """

        if len(self.drivers) == 0:
            raise ValueError("No energy driver defined.")

        self._emit_test_hook("ffg.initial_add_point_rotor_start", {
            "n_molecules": int(len(molecule_specific_information)),
            "cluster_run": bool(self.cluster_run),
        })
        
        def build_rotor_cluster_information(rotors, cluster_components, coupling, dihedral_start, dihedral_end):
            
            clusters = {}
            rotor_to_clusters = {}

            for cluster_id, rotor_ids in enumerate(cluster_components):
                rows = []
                for rotor_id in rotor_ids:
                    rows.extend(rotors[rotor_id].torsion_rows)
                    rotor_to_clusters[rotor_id] = cluster_id

                cluster_type = "independent" if len(rotor_ids) == 1 else "coupled"
                clusters[cluster_id] = RotorClusterDefinition(
                    cluster_id=cluster_id,
                    rotor_ids=tuple(rotor_ids),
                    cluster_type=cluster_type,
                    torsion_rows=tuple(sorted(rows)),
                )

            return RotorClusterInformation(
                dihedral_start=dihedral_start,
                dihedral_end=dihedral_end,
                rotors=rotors,
                clusters=clusters,
                rotor_to_cluster=rotor_to_clusters,
                coupling_matrix=coupling,
            )
        
        def _get_rotor_phase_library(cluster_info):
            out = {}
            for rotor_id, rotor in cluster_info.rotors.items():
                symmetry_order = max(int(rotor.symmetry_order), 1)
                rotor_period = (2.0 * np.pi) / symmetry_order
                print(symmetry_order, rotor_period)
                # state_id = 0 is the anchor state, so only store the nonzero offsets
                # that will become additional cluster-bank points.
                if symmetry_order == 3:
                    out[rotor_id] = tuple([0.0, np.pi/6.0, np.pi/3.0, np.pi/2.0])

            return out
        
        def _build_cluster_angle_library(
            cluster_info,
            reference_internal_values,
            phase_library,
            ):
            out = RotorClusterAngleLibrary()

            for cluster_id, cluster in cluster_info.clusters.items():
                state_bank = []
                rotor_ids = cluster.rotor_ids

                # Build the anchor state first. It always gets state_id = 0.
                anchor_assignment = {}
                
                for rotor_id in rotor_ids:
                    rotor = cluster_info.rotors[rotor_id]
                    for row, torsion in zip(rotor.torsion_rows, rotor.torsion_coords):
                        angle = float(reference_internal_values[row])
                        anchor_assignment[torsion] = angle

                state_bank.append(
                    RotorClusterStateDefinition(
                        cluster_id=cluster_id,
                        state_id=0,
                        cluster_type=cluster.cluster_type,
                        rotor_ids=rotor_ids,
                        angle_assignment=anchor_assignment,
                        dihedrals_to_rotate=None,
                        phase_signature=None,
                        is_anchor=True,
                        label_suffix="anchor",
                    )
                )

                if cluster.cluster_type == "independent":
                    rotor_id = rotor_ids[0]
                    rotor = cluster_info.rotors[rotor_id]
                    sampled_phases = phase_library[rotor_id]
                    dihedrals_to_rotate = [rotor.torsion_coords[0]]

                    state_skipped = 0
                    for phase_index, phase in enumerate(sampled_phases, start=1):
                        if phase == 0.0:
                            state_skipped += 1
                            continue
                        angle_assignment = {}

                        for torsion in rotor.torsion_coords:
                            angle_assignment[torsion] = anchor_assignment[torsion] + float(phase)


                        state_bank.append(
                            RotorClusterStateDefinition(
                                cluster_id=cluster_id,
                                state_id=phase_index - state_skipped,
                                cluster_type=cluster.cluster_type,
                                rotor_ids=rotor_ids,
                                angle_assignment=angle_assignment,
                                dihedrals_to_rotate=dihedrals_to_rotate,
                                phase_signature=[phase],
                                # label_suffix=f"r{rotor_id}_s{phase_index - state_skipped}",
                            )
                        )

                else:
                    joint_phases = []
                    for rotor_id in rotor_ids:
                        joint_phases.append(phase_library[rotor_id])

                    
  
                    state_id = 1

                    for phase_tuple in itertools.product(*joint_phases):
                        if sum(phase_tuple) == 0.0:
                            continue
                        angle_assignment = {}
                        dihedrals_to_rotate = []
                        current_phase_signiture = []
                        for rotor_id, phase in zip(rotor_ids, phase_tuple):
                            current_phase_signiture.append(phase)
                            rotor = cluster_info.rotors[rotor_id]
                            dihedrals_to_rotate.append(rotor.torsion_coords[0])
                            for torsion in rotor.torsion_coords:
                                angle_assignment[torsion] = anchor_assignment[torsion] + float(phase)

                        state_bank.append(
                            RotorClusterStateDefinition(
                                cluster_id=cluster_id,
                                state_id=state_id,
                                cluster_type=cluster.cluster_type,
                                rotor_ids=rotor_ids,
                                angle_assignment=angle_assignment,
                                dihedrals_to_rotate=dihedrals_to_rotate,
                                phase_signature=current_phase_signiture,
                                # label_suffix=f"state_{state_id}",
                            )
                        )
                        state_id += 1

                out.state_banks[cluster_id] = tuple(state_bank)

            return out
        
        def _enumerate_cluster_bank_jobs(cluster_info, cluster_angle_libary):

            jobs = []

            for cluster_id, cluster in cluster_info.clusters.items():
                for state_def in cluster_angle_libary[cluster_id]:

                    jobs.append({
                        "bank_role": "cluster",
                        "cluster_id": cluster_id,
                        "cluster_type": cluster.cluster_type,
                        "state_id": state_def.state_id,
                        "cluster_rotor_ids": state_def.rotor_ids,
                        "angle_assignment": state_def.angle_assignment,
                        "dihedrals_to_rotate":state_def.dihedrals_to_rotate,
                        "phase_signature": state_def.phase_signature,
                        "is_anchor": state_def.is_anchor,
                        "label_suffix": state_def.label_suffix,
                    })

            return jobs
                

        # define impesdriver to determine if stucture should be added:
        
        ## For symmetry groups of periodicty of 3 it is crucial for the interpolation to set the dihedral to the position between 2 extreme points in order to account
        ## for the symmetry correclty using only one reference point!
        
        # create all molecule combinations

        adjusted_molecule = {'gs': [], 'es': []}
        symmetry_mapping_groups = []
        symmetry_exclusion_groups = []

        # comm = MPI.COMM_WORLD
        # rank = comm.Get_rank()
        # root = mpi_master()

        # use the first incomming structure to determine the cluster grouping

        for entries in molecule_specific_information:
            molecule = entries[0]
            basis = entries[1]
            states = entries[2]
            symmetry_point = False

            transition_state = False
            if len(entries) == 5:
                transition_state = entries[4]
            
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
                adjusted_molecule['gs'].append((entries[0], entries[1], 1, None, [0], symmetry_point, entries[3], transition_state)) 
            if any(x > 0 for x in entries[2]): 
                states = [state for state in entries[2] if state > 0]
                adjusted_molecule['es'].append((entries[0], entries[1], 1, None, states, symmetry_point, entries[3], transition_state)) 

        rotor_to_cluster_inf = {'gs': None, 'es': None}

        for state_key, entries in adjusted_molecule.items():
            if len(entries) == 0:
                continue

            drivers = self.drivers[state_key]
    
            org_roots = None
            if isinstance(drivers[0], LinearResponseEigenSolver) or isinstance(drivers[0], TdaEigenSolver):           
                    org_roots = drivers[1].state_deriv_index
            label_counter = 0
            for mol_basis in entries:
                if not mol_basis[5]:
                    label_counter = 0
                if isinstance(drivers[0], LinearResponseEigenSolver) or isinstance(drivers[0], TdaEigenSolver):
                    root_to_follow_calc = mol_basis[4]           
                    drivers[1].state_deriv_index = root_to_follow_calc 
                    # drivers[2].roots_to_follow = root_to_follow_calc
                
                symmetry_mapping_groups = [item for item in range(len(mol_basis[0].get_labels()))]
                symmetry_exclusion_groups = [item for element in symmetry_information[state_key][1] for item in element]

                energies, scf_results, rsp_results = self._compute_energy(drivers[0], mol_basis[0], mol_basis[1])
                # scf_results_mpi = self._compute_energy_mpi_safe(
                #     drivers[0],
                #     mol_basis[0],
                #     mol_basis[1],
                #     collective=True,
                #     phase_name='Energy calcualtion',
                # )

                # energies, scf_results, rsp_results = scf_results_mpi[0], scf_results_mpi[1], scf_results_mpi[2]
                
                if isinstance(drivers[0], LinearResponseEigenSolver) or isinstance(drivers[0], TdaEigenSolver):
                    energies = energies[mol_basis[4]]
                
                gradients = self._compute_gradient(drivers[1], mol_basis[0], mol_basis[1], scf_results, rsp_results)
                
                # gradient_results_mpi = self._compute_gradient_mpi_safe(
                #     drivers[1],
                #     mol_basis[0],
                #     mol_basis[1],
                #     scf_results,
                #     rsp_results,
                #     collective=True,
                #     phase_name='Gradient Calculation',
                # )
                # gradients = gradient_results_mpi
                
                hessians = self._compute_hessian(drivers[2], mol_basis[0], mol_basis[1])
                # hessians_result_mpi = self._compute_hessian_mpi_safe(
                #     drivers[2],
                #     mol_basis[0],
                #     mol_basis[1],
                #     collective=True,
                #     phase_name='hessian calculation',
                # )

                # hessians = hessians_result_mpi
                
                masses = mol_basis[0].get_masses().copy()
                masses_cart = np.repeat(masses, 3)
                inv_sqrt_masses = 1.0 / np.sqrt(masses_cart)
                
                state_spec_dp = {}
                
                # if rank == root:

                for number in range(len(energies)):
                    target_root = mol_basis[4][number]
                    target_file = interpolation_settings[target_root]['imforcefield_file']
                    
                    z_matrix = self.roots_z_matrix[target_root]
                    interpolation_driver = InterpolationDriver(z_matrix)
                    interpolation_driver.update_settings(interpolation_settings[target_root])
                    interpolation_driver.imforcefield_file = target_file
                    
                    sorted_labels = []
                    if target_file in os.listdir(os.getcwd()):
                        org_labels, z_matrix = interpolation_driver.read_labels()
                        labels = [label for label in org_labels if '_cluster' not in label]
                        sorted_labels = sorted(labels, key=lambda x: int(x.split('_')[1]))
                    label = None
                    grad = gradients[number].copy()
                    hess = hessians[number].copy()
                    grad_vec = grad.reshape(-1)         # (3N,)
                    hess_mat = hess.reshape(grad_vec.size, grad_vec.size)
                    mw_grad_vec = inv_sqrt_masses * grad_vec
                    mw_hess_mat = (inv_sqrt_masses[:, None] * hess_mat) * inv_sqrt_masses[None, :]
                    if mol_basis[5] == False:
                        family_label = f'point_{len(sorted_labels) + 1}'
                        label = f'point_{len(sorted_labels) + 1}_core'
                    
                    impes_coordinate = InterpolationDatapoint(z_matrix)
                    impes_coordinate.eq_bond_lengths = self.eq_bond_length
                    impes_coordinate.update_settings(interpolation_settings[target_root])
                    impes_coordinate.cartesian_coordinates = mol_basis[0].get_coordinates_in_bohr()
                    impes_coordinate.imp_int_coordinates = {'bonds': [], 'angles': [], 'dihedrals': [], 'impropers': []}
                    impes_coordinate.inv_sqrt_masses = inv_sqrt_masses
                    impes_coordinate.energy = energies[number]
                    impes_coordinate.gradient =  mw_grad_vec.reshape(grad.shape)
                    impes_coordinate.hessian = mw_hess_mat.reshape(hess.shape)
                    impes_coordinate.transform_gradient_and_hessian()
                    impes_coordinate.family_label = family_label
                    impes_coordinate.point_label = label
                    impes_coordinate.bank_role = "core"
                    trust_radius = self.use_opt_confidence_radius[2]
                    impes_coordinate.confidence_radius = trust_radius
                    
                    if self.use_symmetry and len(self.symmetry_rotors) > 0:
                        rotation_combinations = None
                        dihedral_list = [rotor.torsion_coords[0] for _, rotor in self.symmetry_rotors.items()]
                        
                        from itertools import product
                        rotations = [0.0, 2.0*np.pi/3.0, 4.0*np.pi/3.0]
                        dihedrals = dihedral_list  # your list of rotatable dihedrals
                        rotation_combinations = list(product(rotations, repeat=len(dihedrals)))
                        
                        org_mask = [i for i in range(len(impes_coordinate.z_matrix))]
                        masks = [org_mask]
                        for combo in rotation_combinations:
                            if all(r == 0 for r in combo):
                                continue  # skip the base geometry if desired
                            
                                
                            rot_mol = Molecule(self.molecule.get_labels(), mol_basis[0].get_coordinates_in_bohr(), 'bohr')
                            for angle, dihedral in zip(combo, dihedrals):
                                dihedral_p_1 = (dihedral[0] + 1, dihedral[1] + 1, dihedral[2] + 1, dihedral[3] + 1)
                                rot_mol.set_dihedral(dihedral_p_1, mol_basis[0].get_dihedral(dihedral_p_1, 'radian') - angle, 'radian')
                            target_coordinates = self._calculate_translation_coordinates(rot_mol.get_coordinates_in_bohr())
                            reference_coordinates = self._calculate_translation_coordinates(mol_basis[0].get_coordinates_in_bohr())
                            
                            target_coordinates_core = target_coordinates.copy()
                            reference_coordinates_core = reference_coordinates.copy()
                                
                            target_coordinates_core = np.delete(target_coordinates, symmetry_exclusion_groups, axis=0)
                            reference_coordinates_core = np.delete(reference_coordinates, symmetry_exclusion_groups, axis=0)
                            rotation_matrix = geometric.rotate.get_rot(target_coordinates_core,
                                                                    reference_coordinates_core)
                            # Rotate the data point
                            rotated_coordinates = np.dot(rotation_matrix, target_coordinates.T).T
                            rotated_coordinates = np.ascontiguousarray(rotated_coordinates)
                            
                            mapping_dict_12 = self._perform_symmetry_assignment(symmetry_mapping_groups, symmetry_exclusion_groups, 
                                                                                mol_basis[0].get_coordinates_in_bohr()[symmetry_exclusion_groups], rotated_coordinates[symmetry_exclusion_groups])
                            
                            
                            z_matrix_dict = {tuple(sorted(element)): i 
                                for i, element in enumerate(impes_coordinate.z_matrix)}
                            
                            mapping_dict = [mapping_dict_12]
                            reorded_int_coords = [impes_coordinate.internal_coordinates_values.copy()]

                            for ord in range(len(mapping_dict)):
                                mask = []
                                # inverse_mapping = {v: k for k, v in mapping_dict[ord].items()}
                                reorded_int_coord = np.zeros_like(impes_coordinate.internal_coordinates_values)
                                for i, element in enumerate(impes_coordinate.z_matrix):
                                    # Otherwise, reorder the element
                                    reordered_element = [mapping_dict[ord].get(x, x) for x in element]
                                    key = tuple(sorted(reordered_element))
                                    z_mat_index = z_matrix_dict.get(key)
                                    mask.append(z_mat_index)
                                    reorded_int_coord[i] = (float(reorded_int_coords[0][z_mat_index]))
                                reorded_int_coords.append(reorded_int_coord)
        
                                masks.append(mask)
                        self._finalize_mapping_masks_and_eq_mode(impes_coordinate, masks)
                    else:        
                        org_mask = [i for i in range(len(impes_coordinate.z_matrix))]
                        masks = [org_mask]
                        self._finalize_mapping_masks_and_eq_mode(impes_coordinate, masks)
                    impes_coordinate.write_hdf5(interpolation_settings[target_root]['imforcefield_file'], label)
                    self._emit_test_hook("datapoint_written", {
                        "root": int(target_root),
                        "label": str(label),
                        "energy_hartree": float(impes_coordinate.energy),
                        "cartesian_coordinats": impes_coordinate.cartesian_coordinates,
                        "confidence_radius": float(impes_coordinate.confidence_radius),
                        "bank_role": str(getattr(impes_coordinate, "bank_role", "core")),
                    })
                    # storing the original database
                    impes_coordinate.write_hdf5(f'im_database_{target_root}_org.h5', label)
                    print('Family label part', impes_coordinate.family_label)
                    print('Point label in core part', impes_coordinate.point_label)
                    state_spec_dp[target_root] = impes_coordinate
                
                state_spec_payload = None
                # if rank == root:
                state_spec_payload = {
                    int(state): {
                        "point_label": dp.point_label,
                        "family_label": dp.family_label,
                        "eq_bond_lengths": np.array(dp.eq_bond_lengths, copy=True),
                        "internal_coordinates_values": np.array(dp.internal_coordinates_values, copy=True),
                        "internal_hessian": np.array(dp.internal_hessian, copy=True),
                        "cartesian_coordinates": np.array(dp.cartesian_coordinates, copy=True),
                    }
                    for state, dp in state_spec_dp.items()
                }    
                # state_spec_payload = comm.bcast(
                #    state_spec_payload if rank == root else None,
                #    root,
                # )    
                for state, cur_dp_payload in state_spec_payload.items():
                    cur_dp = type("CorePointMeta", (), {})()
                    cur_dp.point_label = cur_dp_payload["point_label"]
                    cur_dp.family_label = cur_dp_payload["family_label"]
                    cur_dp.eq_bond_lengths = cur_dp_payload["eq_bond_lengths"]
                    cur_dp.internal_coordinates_values = cur_dp_payload["internal_coordinates_values"]
                    cur_dp.internal_hessian = cur_dp_payload["internal_hessian"]
                    cur_cartesian_coordinates = cur_dp_payload["cartesian_coordinates"]



                    point_index = {"core_label": cur_dp.point_label, "cluster_state_labels": {}}
                    
                    nrot = len(self.symmetry_rotors)
                    coupling_map = np.zeros((nrot, nrot), dtype=np.float64)
                    np.fill_diagonal(coupling_map, 1.0)
                    
                    for rotor_i in range(nrot - 1):
                        for rotor_j in range(rotor_i + 1, nrot):
                            a = self.symmetry_rotors[rotor_i]
                            b = self.symmetry_rotors[rotor_j]
                            rotor_coupling = rotor_coupling_score(cur_dp.internal_hessian, a.torsion_rows, b.torsion_rows)
                            coupling_map[a.rotor_id, b.rotor_id] = rotor_coupling
                            coupling_map[b.rotor_id, a.rotor_id] = rotor_coupling

                    clusters = build_rotor_clusters(self.symmetry_rotors, coupling_map, self.rotor_corr_threshold)

                    rotor_to_cluster_inf[state_key] = build_rotor_cluster_information(
                        self.symmetry_rotors,
                        clusters,
                        coupling_map,
                        symmetry_information[state_key][-1][0],
                        symmetry_information[state_key][-1][1],
                    )

                    angle_phases = _get_rotor_phase_library(rotor_to_cluster_inf[state_key])
                    final_clusters_angle_lib = _build_cluster_angle_library(  rotor_to_cluster_inf[state_key], 
                                                                    cur_dp.internal_coordinates_values,
                                                                    angle_phases)
                    cluster_jobs = _enumerate_cluster_bank_jobs(rotor_to_cluster_inf[state_key], final_clusters_angle_lib)
                    

                    for job in cluster_jobs:
                        if job['is_anchor']:
                            continue
         
                        cur_molecule = Molecule(self.molecule.get_labels(), cur_cartesian_coordinates, 'bohr')
                        cur_molecule.set_charge(self.molecule.get_charge())
                        cur_molecule.set_multiplicity(self.molecule.get_multiplicity())
                        current_angles_setting = []
                        constraints = []
                        for dihedral in job["dihedrals_to_rotate"]:
                            angle = job["angle_assignment"].get(dihedral)
                            cur_molecule.set_dihedral(
                                [
                                    dihedral[0] + 1,
                                    dihedral[1] + 1,
                                    dihedral[2] + 1,
                                    dihedral[3] + 1,
                                ],
                                angle,
                                "radian",
                            )
                            current_angles_setting.append(angle)
                            constraints.append(dihedral)
                        
                        current_basis = MolecularBasis.read(cur_molecule, basis.get_main_basis_label())

                        if not mol_basis[7]:
                            if isinstance(drivers[0], ScfRestrictedDriver):

                                energies, scf_results, rsp_results = self._compute_energy(drivers[0], cur_molecule, current_basis)
                            #     _, scf_results, _ = self._compute_energy_mpi_safe(
                            #         drivers[0],
                            #         cur_molecule,
                            #         current_basis,
                            #         collective=True,
                            #         phase_name='Energy calcualtion',
                            #     )

                                opt_results = self._run_optimization(
                                                drivers[0],
                                                cur_molecule,
                                                constraints=constraints,
                                                index_offset=1,
                                                compute_args=(current_basis, scf_results),
                                                # source_molecule=molecule,
                                                )

                            #     opt_results_mpi = self._run_optimization_mpi_safe(
                            #                     drivers[0],
                            #                     cur_molecule,
                            #                     constraints=constraints,
                            #                     index_offset=1,
                            #                     compute_args=(current_basis, scf_results),
                            #                     source_molecule=molecule,
                            #                     collective=True,
                            #                     phase_name='optimization tag',
                            #                 )
                            #     opt_results = opt_results_mpi[1]
                                optimized_molecule = opt_results['final_molecule']
                                cur_molecule = optimized_molecule
                                
                            elif isinstance(drivers[0], XtbDriver):
                                
                                opt_results = self._run_optimization(
                                                drivers[0],
                                                cur_molecule,
                                                constraints=constraints,
                                                index_offset=1,
                                                # compute_args=(current_basis, scf_results),
                                                # source_molecule=molecule,
                                                )
                                
                                # opt_results_mpi = self._run_optimization_mpi_safe(
                                #     drivers[0],
                                #     cur_molecule,
                                #     constraints=constraints,
                                #     index_offset=1,
                                #     source_molecule=molecule,
                                #     collective=True,
                                #     phase_name='optimization tag',
                                # )
                            
                                # opt_results = opt_results_mpi[1]
                                optimized_molecule = opt_results['final_molecule']

                                cur_molecule = optimized_molecule
                        print(cur_molecule.get_xyz_string())
                        current_basis = MolecularBasis.read(cur_molecule, basis.get_main_basis_label())

                        energies, scf_results, rsp_results = self._compute_energy(drivers[0], cur_molecule, current_basis)
                        # scf_results_mpi = self._compute_energy_mpi_safe(
                        #     drivers[0],
                        #     cur_molecule,
                        #     current_basis,
                        #     collective=True,
                        #     phase_name='Energy calcualtion',
                        # )

                        # energies, scf_results, rsp_results = scf_results_mpi[0], scf_results_mpi[1], scf_results_mpi[2]

                        if isinstance(drivers[0], LinearResponseEigenSolver) or isinstance(drivers[0], TdaEigenSolver):
                            energies = energies[mol_basis[4]]

                        
                        gradients = self._compute_gradient(drivers[1], cur_molecule, current_basis, scf_results, rsp_results)
                        hessians = self._compute_hessian(drivers[2], cur_molecule, current_basis)
                        
                        # gradient_results_mpi = self._compute_gradient_mpi_safe(
                        #     drivers[1],
                        #     cur_molecule,
                        #     current_basis,
                        #     scf_results,
                        #     rsp_results,
                        #     collective=True,
                        #     phase_name='Gradient Calculation',
                        # )
                        # gradients = gradient_results_mpi
                        
                        # hessians_result_mpi = self._compute_hessian_mpi_safe(
                        #     drivers[2],
                        #     cur_molecule,
                        #     current_basis,
                        #     collective=True,
                        #     phase_name='hessian calculation',
                        # )
                        # hessians = hessians_result_mpi#

                        state_spec_dp = {}
                        # if rank == root:
                        masks = []
                        if self.use_symmetry and len(self.symmetry_rotors) > 0:
                            rotation_combinations = None
                            dihedral_list = [rotor.torsion_coords[0] for _, rotor in self.symmetry_rotors.items()]
                            
                            from itertools import product
                            rotations = [0.0, 2.0*np.pi/3.0, 4.0*np.pi/3.0]
                            dihedrals = dihedral_list  # your list of rotatable dihedrals
                            rotation_combinations = list(product(rotations, repeat=len(dihedrals)))
                            
                            org_mask = [i for i in range(len(impes_coordinate.z_matrix))]
                            masks = [org_mask]
                            for combo in rotation_combinations:
                                if all(r == 0 for r in combo):
                                    continue  # skip the base geometry if desired
                                
                                    
                                rot_mol = Molecule(self.molecule.get_labels(), cur_molecule.get_coordinates_in_bohr(), 'bohr')
                                for angle, dihedral in zip(combo, dihedrals):
            
                                    dihedral_p_1 = (dihedral[0] + 1, dihedral[1] + 1, dihedral[2] + 1, dihedral[3] + 1)
                                    rot_mol.set_dihedral(dihedral_p_1, cur_molecule.get_dihedral(dihedral_p_1, 'radian') - angle, 'radian')
                                target_coordinates = self._calculate_translation_coordinates(rot_mol.get_coordinates_in_bohr())
                                reference_coordinates = self._calculate_translation_coordinates(cur_molecule.get_coordinates_in_bohr())
                                
                                target_coordinates_core = target_coordinates.copy()
                                reference_coordinates_core = reference_coordinates.copy()
                                    
                                target_coordinates_core = np.delete(target_coordinates, symmetry_exclusion_groups, axis=0)
                                reference_coordinates_core = np.delete(reference_coordinates, symmetry_exclusion_groups, axis=0)
                                rotation_matrix = geometric.rotate.get_rot(target_coordinates_core,
                                                                        reference_coordinates_core)
                                # Rotate the data point
                                rotated_coordinates = np.dot(rotation_matrix, target_coordinates.T).T
                                rotated_coordinates = np.ascontiguousarray(rotated_coordinates)
                                
                                mapping_dict_12 = self._perform_symmetry_assignment(symmetry_mapping_groups, symmetry_exclusion_groups, 
                                                                                    cur_molecule.get_coordinates_in_bohr()[symmetry_exclusion_groups], rotated_coordinates[symmetry_exclusion_groups])
                                
                                
                                z_matrix_dict = {tuple(sorted(element)): i 
                                    for i, element in enumerate(impes_coordinate.z_matrix)}
                                
                                mapping_dict = [mapping_dict_12]
                                reorded_int_coords = [impes_coordinate.internal_coordinates_values.copy()]

                                for ord in range(len(mapping_dict)):
                                    mask = []
                                    reorded_int_coord = np.zeros_like(impes_coordinate.internal_coordinates_values)
                                    for i, element in enumerate(impes_coordinate.z_matrix):
                                        # Otherwise, reorder the element
                                        reordered_element = [mapping_dict[ord].get(x, x) for x in element]
                                        key = tuple(sorted(reordered_element))
                                        z_mat_index = z_matrix_dict.get(key)
                                        mask.append(z_mat_index)
                                        reorded_int_coord[i] = (float(reorded_int_coords[0][z_mat_index]))
                                    reorded_int_coords.append(reorded_int_coord)
        
                                    masks.append(mask)
                            
                        else:        
                            org_mask = [i for i in range(len(impes_coordinate.z_matrix))]
                            masks = [org_mask]
        
                        for number in range(len(energies)):
                            target_root = mol_basis[4][number]
                            print(cur_dp.point_label, job["bank_role"], job["cluster_id"], job["cluster_type"], job["cluster_rotor_ids"], job["state_id"], job["label_suffix"])
                            grad = gradients[number].copy()
                            hess = hessians[number].copy()
                            grad_vec = grad.reshape(-1)         # (3N,)
                            hess_mat = hess.reshape(grad_vec.size, grad_vec.size)
                            mw_grad_vec = inv_sqrt_masses * grad_vec
                            mw_hess_mat = (inv_sqrt_masses[:, None] * hess_mat) * inv_sqrt_masses[None, :]
                            impes_coordinate = InterpolationDatapoint(z_matrix)
                            impes_coordinate.eq_bond_lengths = cur_dp.eq_bond_lengths
                            impes_coordinate.update_settings(interpolation_settings[target_root])
                            impes_coordinate.cartesian_coordinates = cur_molecule.get_coordinates_in_bohr()
                            impes_coordinate.imp_int_coordinates = {'bonds': [], 'angles': [], 'dihedrals': [], 'impropers': []}
                            impes_coordinate.inv_sqrt_masses = inv_sqrt_masses
                            impes_coordinate.energy = energies[number]
                            impes_coordinate.gradient =  mw_grad_vec.reshape(grad.shape)
                            impes_coordinate.hessian = mw_hess_mat.reshape(hess.shape)
                            impes_coordinate.transform_gradient_and_hessian()
                            
                            label = f"{cur_dp.point_label}_cluster_{job["cluster_id"]}_state_{job["state_id"]}"
                            if job["label_suffix"]:
                                label += f"_{job["label_suffix"]}"
                            print('here is the label', label)
                            impes_coordinate.family_label = cur_dp.family_label
                            impes_coordinate.point_label = label
                            impes_coordinate.bank_role = job["bank_role"]
                            impes_coordinate.cluster_id = job["cluster_id"]
                            impes_coordinate.cluster_type = job["cluster_type"]
                            impes_coordinate.dihedrals_to_rotate = job["dihedrals_to_rotate"]
                            impes_coordinate.cluster_rotor_ids = job["cluster_rotor_ids"]
                            impes_coordinate.cluster_state_id = job["state_id"]
                            trust_radius = self.use_opt_confidence_radius[2]
                            
                            for phase in job["phase_signature"]:
                                print('statements', abs(abs(phase) - np.pi/6) <= 1e-3, abs(abs(phase) - np.pi/2) <= 1e-3)
                                if abs(abs(phase) - np.pi/6) <= 1e-3 or abs(abs(phase) - np.pi/2) <= 1e-3:
                                    trust_radius = 0.15
                                    break
                            print('Given trust_radius', trust_radius)
                            impes_coordinate.confidence_radius = trust_radius
                            
                            
                            self._finalize_mapping_masks_and_eq_mode(impes_coordinate, masks)
                            impes_coordinate.write_hdf5(interpolation_settings[target_root]['imforcefield_file'], label)
                            self._emit_test_hook("datapoint_written", {
                                "root": int(target_root),
                                "label": str(label),
                                "energy_hartree": float(impes_coordinate.energy),
                                "confidence_radius": float(impes_coordinate.confidence_radius),
                                "bank_role": str(getattr(impes_coordinate, "bank_role", "core")),
                            })
                            impes_coordinate.write_hdf5(f'im_database_{target_root}_org.h5', label)
                            print('Symmetry rotor label', label)
                            point_index["cluster_state_labels"].setdefault(int(job["cluster_id"]), {})
                            point_index["cluster_state_labels"][int(job["cluster_id"])][int(job["state_id"])] = label

                    # if rank == root:
                    self._write_cluster_registry_for_family(interpolation_settings[target_root]['imforcefield_file'], target_root, cur_dp.family_label, rotor_to_cluster_inf[state_key], final_clusters_angle_lib, point_index)
                    self._write_cluster_registry_for_family(f'im_database_{target_root}_org.h5', target_root, cur_dp.family_label, rotor_to_cluster_inf[state_key], final_clusters_angle_lib, point_index)
                    
 
                    # comm.barrier()
                    label_counter += 1
                    if any(root > 0 for root in self.roots_to_follow): 
                        if isinstance(drivers[0], LinearResponseEigenSolver) or isinstance(drivers[0], TdaEigenSolver):           
                            drivers[1].state_deriv_index = org_roots

    def add_point(self, molecule_specific_information, interpolation_settings, symmetry_information={}):
        """ Adds a new point to the database.

            :param molecule:
                the molecule.
            :param imforcefielddatafile:
                Datafile containing the information of the IM forcefield.

        """

        if len(self.drivers) == 0:
            raise ValueError("No energy driver defined.")
           
        # define impesdriver to determine if stucture should be added:
        
        ## For symmetry groups of periodicty of 3 it is crucial for the interpolation to set the dihedral to the position between 2 extreme points in order to account
        ## for the symmetry correclty using only one reference point!
        
        # create all molecule combinations

        if self.cluster_run:
            return self.add_point_rotor(molecule_specific_information, interpolation_settings, symmetry_information)

        self._emit_test_hook("ffg.initial_add_point_start", {
            "n_molecules": int(len(molecule_specific_information)),
            "roots": list(self.roots_to_follow),
        })
        adjusted_molecule = {'gs': [], 'es': []}
        symmetry_mapping_groups = []
        symmetry_exclusion_groups = []

        # comm = MPI.COMM_WORLD
        # rank = comm.Get_rank()
        # root = mpi_master()

        for entries in molecule_specific_information:
            molecule = entries[0]
            basis = entries[1]
            states = entries[2]
            outside_constraints = entries[3]
            symmetry_point = False
            transition_state = False
            if len(entries) == 5:
                transition_state = entries[4]

            if self.eq_bond_length is None:
                self.eq_bond_length = []
                for idx, element in enumerate(self.roots_z_matrix[0]['bonds']):
                    if len(element) == 2 and self.use_minimized_structures[0] and self.eq_bond_length_irc_bonds is not None and element not in self.eq_bond_length_irc_bonds:
                        self.eq_bond_length.append(molecule.get_distance([element[0] + 1, element[1] + 1], 'bohr'))
                    if len(element) == 2 and self.use_minimized_structures[0] and self.eq_bond_length_irc_bonds is not None and element in self.eq_bond_length_irc_bonds:
                        self.eq_bond_length.append(molecule_specific_information[-1][0].get_distance([element[0] + 1, element[1] + 1], 'bohr'))
                    elif len(element) == 2 and self.use_minimized_structures[0]:
                        self.eq_bond_length.append(molecule.get_distance([element[0] + 1, element[1] + 1], 'bohr'))
                    elif len(element) == 2:
                        self.eq_bond_length.append(0.0)

            if 0 in states and len(symmetry_information['gs']) != 0 and len(symmetry_information['gs'][2]) != 0 or 1 in states and len(symmetry_information['gs']) != 0 and self.drivers['es'][0].spin_flip and len(symmetry_information['gs'][2]) != 0:
                symmetry_mapping_groups = [item for item in range(len(molecule.get_labels()))]
                symmetry_exclusion_groups = [item for element in symmetry_information['gs'][1] for item in element]
                sym_dihedrals, periodicites, _, _, _ = self._adjust_symmetry_dihedrals(molecule, symmetry_information['gs'][1], symmetry_information['gs'][5],  self.roots_z_matrix[states[0]])
                # Generate all combinations
                keys = list(sym_dihedrals.keys())
                values = list(sym_dihedrals.values())
                combinations = list(itertools.product(*values))
                # Convert to list of dictionaries
                molecule_configs = [dict(zip(keys, combo)) for combo in combinations] # if all(element == combo[0] for element in combo)

                for i, molecule_config in enumerate(molecule_configs):
                    cur_molecule = Molecule.from_xyz_string(molecule.get_xyz_string())
                    cur_molecule.set_charge(molecule.get_charge())
                    cur_molecule.set_multiplicity(molecule.get_multiplicity())
                    dihedral_to_change = []
                    constraints = []
                    constraints.extend(outside_constraints)
                    for dihedral, angle in molecule_config.items():
                        
                        dihedral_p_1 = [dihedral[0] + 1, dihedral[1] + 1, dihedral[2] + 1, dihedral[3] + 1]
                        print('dihedral to change', dihedral_p_1)
                        opt_dihedral_angle = molecule.get_dihedral(dihedral_p_1, 'radian')
     
                        print('Dihedral creation', dihedral, opt_dihedral_angle + angle, angle)
                        
                        cur_molecule.set_dihedral(dihedral_p_1, opt_dihedral_angle + angle, 'radian')
                        dihedral_to_change.append(dihedral_p_1)

                        constraints.append(dihedral_p_1)
                    
                        if self.symmetry_information is not None:
                            print('opt_dihedral_angle', opt_dihedral_angle)
                            self.symmetry_information['gs'][8][dihedral].append(opt_dihedral_angle)

                    if i > 0:
                        symmetry_point = True
                    
                    current_basis = MolecularBasis.read(cur_molecule, basis.get_main_basis_label())

                    opt_results = None
                    scf_results = None
                    if transition_state:
                        if isinstance(self.drivers['gs'][0], ScfRestrictedDriver):
                            _, scf_results, _ = self._compute_energy(self.drivers['gs'][0], cur_molecule, current_basis)
                            # _, scf_results, _ = self._compute_energy_mpi_safe(
                            #     self.drivers['gs'][0],
                            #     cur_molecule,
                            #     current_basis,
                            #     collective=True,
                            #     phase_name='Energy calcualtion',
                            # )
                            opt_results = self._run_optimization(self.drivers['gs'][0], 
                                                                 cur_molecule, 
                                                                 constraints, 
                                                                 index_offset=0, 
                                                                 compute_args=(current_basis, scf_results))
                            # opt_results_mpi = self._run_optimization_mpi_safe(
                            #                 self.drivers['gs'][0],
                            #                 cur_molecule,
                            #                 constraints=constraints,
                            #                 index_offset=0,
                            #                 compute_args=(current_basis, scf_results),
                            #                 source_molecule=molecule,
                            #                 collective=True,
                            #                 phase_name='optimization tag',
                            #             )
                            # print(opt_results_mpi[1]['opt_energies'][-1])
                            # opt_results = opt_results_mpi[1]
                            optimized_molecule = opt_results['final_molecule']
                            cur_molecule = optimized_molecule

                        elif isinstance(self.drivers['gs'][0], XtbDriver):

                            opt_results = self._run_optimization(self.drivers['gs'][0], 
                                                                 cur_molecule, 
                                                                 constraints, 
                                                                 index_offset=0)
                            # opt_results_mpi = self._run_optimization_mpi_safe(
                            #     self.drivers['gs'][0],
                            #     cur_molecule,
                            #     constraints=constraints,
                            #     index_offset=0,
                            #     source_molecule=molecule,
                            #     collective=True,
                            #     phase_name='optimization tag',
                            # )

                            # print(opt_results_mpi[1]['opt_energies'][-1])
                            # opt_results = opt_results_mpi[1]
                            optimized_molecule = opt_results['final_molecule']

                            cur_molecule = optimized_molecule
                        print('Optimized molecule', cur_molecule.get_xyz_string())

                    current_basis = MolecularBasis.read(cur_molecule, basis.get_main_basis_label())

                    adjusted_molecule['gs'].append((cur_molecule, current_basis, periodicites[dihedral],  dihedral_to_change, entries[2], symmetry_point, outside_constraints))
            
            if 1 in states and len(symmetry_information['es']) != 0 and len(symmetry_information['es'][2]) != 0 and self.drivers['es'][0].spin_flip is False:
                symmetry_mapping_groups = [item for item in range(len(molecule.get_labels()))]
                symmetry_exclusion_groups = [item for element in symmetry_information['es'][1] for item in element]
                sym_dihedrals, periodicites, _, _, _ = self._adjust_symmetry_dihedrals(molecule, symmetry_information['es'][1], symmetry_information['es'][5], self.roots_z_matrix[1])
                
                # Generate all combinations
                keys = list(sym_dihedrals.keys())
                values = list(sym_dihedrals.values())
                combinations = list(itertools.product(*values))
                # Convert to list of dictionaries
                molecule_configs = [dict(zip(keys, combo)) for combo in combinations] # if all(element == combo[0] for element in combo)
                print('MOlecule configs', molecule_configs)
                for i, molecule_config in enumerate(molecule_configs):
                    cur_molecule = Molecule.from_xyz_string(molecule.get_xyz_string())
                    cur_molecule.set_charge(molecule.get_charge())
                    cur_molecule.set_multiplicity(molecule.get_multiplicity())
                    dihedral_to_change = []
                    for dihedral, angle in molecule_config.items():
                        print('Dihedral', dihedral)
                        opt_dihedral_angle = cur_molecule.get_dihedral([dihedral[0] + 1, dihedral[1] + 1, dihedral[2] + 1, dihedral[3] + 1], 'radian')
                        cur_molecule.set_dihedral([dihedral[0] + 1, dihedral[1] + 1, dihedral[2] + 1, dihedral[3] + 1], opt_dihedral_angle + angle, 'radian')
                        dihedral_to_change.append([dihedral[0] + 1, dihedral[1] + 1, dihedral[2] + 1, dihedral[3] + 1])
                        print('påt_dihedral', opt_dihedral_angle + angle)
                        if self.symmetry_information is not None:
                            self.symmetry_information['es'][8][dihedral].append(opt_dihedral_angle + angle)
                    if i > 0:
                        symmetry_point = True
                    current_basis = MolecularBasis.read(cur_molecule, basis.get_main_basis_label())
                    adjusted_molecule['es'].append((cur_molecule, current_basis, periodicites[dihedral],  dihedral_to_change, entries[2], symmetry_point))
            elif any(x > 1 for x in states) and len(symmetry_information['es']) != 0 and len(symmetry_information['es'][2]) != 0 and 1 in states and self.drivers['es'][0].spin_flip is True:
                symmetry_mapping_groups = [item for item in range(len(molecule.get_labels()))]
                symmetry_exclusion_groups = [item for element in symmetry_information['es'][1] for item in element]
                sym_dihedrals, periodicites, _, _, _ = self._adjust_symmetry_dihedrals(molecule, symmetry_information['es'][1], symmetry_information['es'][5], self.roots_z_matrix[1])
                
                # Generate all combinations
                keys = list(sym_dihedrals.keys())
                values = list(sym_dihedrals.values())
                combinations = list(itertools.product(*values))
                # Convert to list of dictionaries
                molecule_configs = [dict(zip(keys, combo)) for combo in combinations] # if all(element == combo[0] for element in combo)
                print('MOlecule configs', molecule_configs)
                for i, molecule_config in enumerate(molecule_configs):
                    cur_molecule = Molecule.from_xyz_string(molecule.get_xyz_string())
                    cur_molecule.set_charge(molecule.get_charge())
                    cur_molecule.set_multiplicity(molecule.get_multiplicity())
                    dihedral_to_change = []
                    for dihedral, angle in molecule_config.items():
                        print('Dihedral', dihedral)
                        opt_dihedral_angle = cur_molecule.get_dihedral([dihedral[0] + 1, dihedral[1] + 1, dihedral[2] + 1, dihedral[3] + 1], 'radian')
                        cur_molecule.set_dihedral([dihedral[0] + 1, dihedral[1] + 1, dihedral[2] + 1, dihedral[3] + 1], opt_dihedral_angle + angle, 'radian')
                        dihedral_to_change.append([dihedral[0] + 1, dihedral[1] + 1, dihedral[2] + 1, dihedral[3] + 1])
                        print('påt_dihedral', opt_dihedral_angle + angle)
                        if self.symmetry_information is not None:
                            self.symmetry_information['es'][8][dihedral].append(opt_dihedral_angle + angle)
                    if i > 0:
                        symmetry_point = True
                    
                    current_basis = MolecularBasis.read(cur_molecule, basis.get_main_basis_label())
                    adjusted_molecule['es'].append((cur_molecule, current_basis, periodicites[dihedral],  dihedral_to_change, entries[2], symmetry_point))
            
            elif not symmetry_point:
                if 0 in entries[2]:
                    adjusted_molecule['gs'].append((entries[0], entries[1], 1, None, [0], symmetry_point, entries[3])) 
                if any(x > 0 for x in entries[2]): 
                    states = [state for state in entries[2] if state > 0]
                    adjusted_molecule['es'].append((entries[0], entries[1], 1, None, states, symmetry_point, entries[3])) 

        for key, entries in adjusted_molecule.items():
            if len(entries) == 0:
                continue
            drivers = self.drivers[key]
            org_roots = None
            if isinstance(drivers[0], LinearResponseEigenSolver) or isinstance(drivers[0], TdaEigenSolver):           
                    org_roots = drivers[1].state_deriv_index
            label_counter = 0
            for mol_basis in entries:
                if not mol_basis[5]:
                    label_counter = 0
                if isinstance(drivers[0], LinearResponseEigenSolver) or isinstance(drivers[0], TdaEigenSolver):
                    root_to_follow_calc = mol_basis[4]           
                    drivers[1].state_deriv_index = root_to_follow_calc 
                
                energies, scf_results, rsp_results = self._compute_energy(drivers[0], mol_basis[0], mol_basis[1])

                # scf_results_mpi = self._compute_energy_mpi_safe(
                #     drivers[0],
                #     mol_basis[0],
                #     mol_basis[1],
                #     collective=True,
                #     phase_name='Energy calcualtion',
                # )

                # energies, scf_results, rsp_results = scf_results_mpi[0], scf_results_mpi[1], scf_results_mpi[2]
                
                if isinstance(drivers[0], LinearResponseEigenSolver) or isinstance(drivers[0], TdaEigenSolver):
                    energies = energies[mol_basis[4]]

                gradients = self._compute_gradient(drivers[1], mol_basis[0], mol_basis[1], scf_results, rsp_results)
                hessians = self._compute_hessian(drivers[2], mol_basis[0], mol_basis[1])
                # gradient_results_mpi = self._compute_gradient_mpi_safe(
                #     drivers[1],
                #     mol_basis[0],
                #     mol_basis[1],
                #     scf_results,
                #     rsp_results,
                #     collective=True,
                #     phase_name='Gradient Calculation',
                # )
                # gradients = gradient_results_mpi
                
                # hessians_result_mpi = self._compute_hessian_mpi_safe(
                #     drivers[2],
                #     mol_basis[0],
                #     mol_basis[1],
                #     collective=True,
                #     phase_name='hessian calculation',
                # )

                # hessians = hessians_result_mpi

                masses = mol_basis[0].get_masses().copy()
                masses_cart = np.repeat(masses, 3)
                inv_sqrt_masses = 1.0 / np.sqrt(masses_cart)
                # writing_information = None
                # if rank == root:
                    # writing_information = []
                for number in range(len(energies)):
                    target_root = mol_basis[4][number]
                    target_file = interpolation_settings[target_root]['imforcefield_file']
                    
                    z_matrix = self.roots_z_matrix[target_root]
                    interpolation_driver = InterpolationDriver(z_matrix)
                    interpolation_driver.update_settings(interpolation_settings[target_root])
                    interpolation_driver.imforcefield_file = target_file
                    
                    sorted_labels = []
                    if target_file in os.listdir(os.getcwd()):
                        org_labels, z_matrix = interpolation_driver.read_labels()
                        labels = [label for label in org_labels if '_symmetry' not in label]
                        sorted_labels = sorted(labels, key=lambda x: int(x.split('_')[1]))
                    label = None
                    grad = gradients[number].copy()
                    hess = hessians[number].copy()
                    grad_vec = grad.reshape(-1)         # (3N,)
                    hess_mat = hess.reshape(grad_vec.size, grad_vec.size)
                    mw_grad_vec = inv_sqrt_masses * grad_vec
                    mw_hess_mat = (inv_sqrt_masses[:, None] * hess_mat) * inv_sqrt_masses[None, :]
                    
    
                    if mol_basis[5] == False:
                        label = f'point_{len(sorted_labels) + 1}'
                        old_label = f'point_{len(sorted_labels) + 1}'
                    else:
                        label = f'{old_label}_symmetry_{label_counter}'
                    impes_coordinate = InterpolationDatapoint(z_matrix)
                    impes_coordinate.eq_bond_lengths = self.eq_bond_length
                    impes_coordinate.update_settings(interpolation_settings[target_root])
                    impes_coordinate.cartesian_coordinates = mol_basis[0].get_coordinates_in_bohr()
                    impes_coordinate.imp_int_coordinates = {'bonds': [], 'angles': [], 'dihedrals': [], 'impropers': []}
                    impes_coordinate.inv_sqrt_masses = inv_sqrt_masses
                    impes_coordinate.energy = energies[number]
                    impes_coordinate.gradient =  mw_grad_vec.reshape(grad.shape)
                    impes_coordinate.hessian = mw_hess_mat.reshape(hess.shape)
                    impes_coordinate.transform_gradient_and_hessian()
                    if "_symmetry" in label:
                        impes_coordinate.bank_role = "symmetry"
                    
                    if mol_basis[2] > 1:
                        
                        rotation_combinations = None
                        
                        if mol_basis[2] == 3:
                            from itertools import product
                            rotations = [0.0, 2*np.pi/3, 4*np.pi/3]
                            dihedrals = mol_basis[3]  # your list of rotatable dihedrals
                            rotation_combinations = list(product(rotations, repeat=len(dihedrals)))
                        
                        if mol_basis[2] == 2:
                            from itertools import product
                            rotations = [0.0, np.pi]
                            dihedrals = mol_basis[3]  # your list of rotatable dihedrals
                            rotation_combinations = list(product(rotations, repeat=len(dihedrals)))
                        
                        org_mask = [i for i in range(len(impes_coordinate.z_matrix))]
                        masks = [org_mask]
                        for combo in rotation_combinations:
                            if all(r == 0 for r in combo):
                                continue  # skip the base geometry if desired
                            
                                
                            rot_mol = Molecule(self.molecule.get_labels(), mol_basis[0].get_coordinates_in_bohr(), 'bohr')
                            for angle, dihedral in zip(combo, dihedrals):
                                rot_mol.set_dihedral(dihedral, mol_basis[0].get_dihedral(dihedral, 'radian') - angle, 'radian')
                            target_coordinates = self._calculate_translation_coordinates(rot_mol.get_coordinates_in_bohr())
                            reference_coordinates = self._calculate_translation_coordinates(mol_basis[0].get_coordinates_in_bohr())
                            
                            target_coordinates_core = target_coordinates.copy()
                            reference_coordinates_core = reference_coordinates.copy()
                                
                            target_coordinates_core = np.delete(target_coordinates, symmetry_exclusion_groups, axis=0)
                            reference_coordinates_core = np.delete(reference_coordinates, symmetry_exclusion_groups, axis=0)
                            rotation_matrix = geometric.rotate.get_rot(target_coordinates_core,
                                                                    reference_coordinates_core)
                            # Rotate the data point
                            rotated_coordinates = np.dot(rotation_matrix, target_coordinates.T).T
                            rotated_coordinates = np.ascontiguousarray(rotated_coordinates)
                            mapping_dict_12 = self._perform_symmetry_assignment(symmetry_mapping_groups, symmetry_exclusion_groups, 
                                                                                mol_basis[0].get_coordinates_in_bohr()[symmetry_exclusion_groups], rotated_coordinates[symmetry_exclusion_groups])
                            
                            
                            z_matrix_dict = {tuple(sorted(element)): i 
                                for i, element in enumerate(impes_coordinate.z_matrix)}
                            
                            mapping_dict = [mapping_dict_12]
                            reorded_int_coords = [impes_coordinate.internal_coordinates_values.copy()]
        
                            for ord in range(len(mapping_dict)):
                                mask = []
                                inverse_mapping = {v: k for k, v in mapping_dict[ord].items()}
                                reorded_int_coord = np.zeros_like(impes_coordinate.internal_coordinates_values)
                                for i, element in enumerate(impes_coordinate.z_matrix):
                                    # Otherwise, reorder the element
                                    reordered_element = [inverse_mapping.get(x, x) for x in element]
                                    key = tuple(sorted(reordered_element))
                                    z_mat_index = z_matrix_dict.get(key)
                                    mask.append(z_mat_index)
                                    reorded_int_coord[i] = (float(reorded_int_coords[0][z_mat_index]))
                                reorded_int_coords.append(reorded_int_coord)
                                masks.append(mask)
                            
                        self._finalize_mapping_masks_and_eq_mode(impes_coordinate, masks)
                    else:
                        org_mask = [i for i in range(len(impes_coordinate.z_matrix))]
                        masks = [org_mask]
                        self._finalize_mapping_masks_and_eq_mode(impes_coordinate, masks)
                            
                    trust_radius = self.use_opt_confidence_radius[2]
                    impes_coordinate.confidence_radius = trust_radius
                    
                    impes_coordinate.write_hdf5(target_file, label)
                    self._emit_test_hook("datapoint_written", {
                        "root": int(target_root),
                        "label": str(label),
                        "energy_hartree": float(impes_coordinate.energy),
                        "confidence_radius": float(impes_coordinate.confidence_radius),
                        "bank_role": str(getattr(impes_coordinate, "bank_role", "core")),
                    })
                    impes_coordinate.write_hdf5(f'im_database_{target_root}_org.h5', label)
                    interpolation_driver.imforcefield_file = target_file
    
                    labels, z_matrix = interpolation_driver.read_labels()
                    
                    print(f"Database expansion with {', '.join(labels)}")
                    # writing_information.append({
                    #     'root': int()
                    # })  
            
                # writing_information = comm.bcast(writing_information if rank == root else None, root)
                # comm.barrier()
                label_counter += 1
                if any(root > 0 for root in self.roots_to_follow): 
                    if isinstance(drivers[0], LinearResponseEigenSolver) or isinstance(drivers[0], TdaEigenSolver):           
                        drivers[1].state_deriv_index = org_roots
                    
        self._emit_test_hook("ffg.initial_add_point_done", {
                            "label_count_per_root": self._label_count_per_root(),
                        })
        
    def _calculate_translation_coordinates(self, cart_coord):
        """Center the molecule by translating its geometric center to (0, 0, 0)."""
        center = np.mean(cart_coord, axis=0)
        translated_coordinates = cart_coord - center

        return translated_coordinates

    def _perform_symmetry_assignment(self, atom_map, sym_group, reference_group, datapoint_group):
        """ Performs the atom mapping. """
        from scipy.optimize import linear_sum_assignment
        new_map = np.array(atom_map.copy())
        mapping_dict = {}
        # cost = self.get_dihedral_cost(atom_map, sym_group, non_group_atoms)
        cost = np.linalg.norm(datapoint_group[:, np.newaxis, :] - reference_group[np.newaxis, :, :], axis=2)
        row, col = linear_sum_assignment(cost)

        if not np.equal(row, col).all():

            
            # atom_maps = self.linear_assignment_solver(cost)

            reordred_arr = np.array(sym_group)[col]
            new_map[sym_group] = new_map[reordred_arr]

            mapping_dict = {org: new for org, new in zip(np.array(sym_group), reordred_arr)}
        
        return mapping_dict

    def _adjust_symmetry_dihedrals(self, molecule, symmetry_groups, rot_bonds, z_matrix):
        
        def symmetry_group_dihedral(reference_set, dihedrals, rot_bonds):
            rot_bond_set = {frozenset(bond) for bond in rot_bonds}

            filtered_dihedrals = []
            for d in dihedrals:
                middle_bond = frozenset([d[1], d[2]])

                if middle_bond in rot_bond_set:
                    common_elements = [x for x in [d[0], d[3]] if x in reference_set]
                    if len(common_elements) == 1:
                        filtered_dihedrals.append(d)
            return filtered_dihedrals
        
        all_dihedrals = [element for element in z_matrix['dihedrals']]

        symmetry_group_dihedral_dict = {} 

        sym_groups_dih = {3: {}}

        angles_to_set = {}
        periodicities = {}
        dihedral_groups = {2: {}, 3: {}}

        for symmetry_group in symmetry_groups:
            
            symmetry_group_dihedral_list = symmetry_group_dihedral(symmetry_group, all_dihedrals, rot_bonds)

            symmetry_group_dihedral_dict[tuple(symmetry_group)] = symmetry_group_dihedral_list
            
            if len(symmetry_group) == 3:

                # angles_to_set[symmetry_group_dihedral_list[0]] = ([0.0, np.pi/3.0])
                angles_to_set[tuple(symmetry_group_dihedral_list[0])] = ([0.0, np.pi/6.0, np.pi/3.0, np.pi/2.0])

                periodicities[tuple(symmetry_group_dihedral_list[0])] = 3
                dihedral_groups[3][tuple(symmetry_group_dihedral_list[0][1:3])] = [tuple(sorted(element, reverse=False)) for element in symmetry_group_dihedral_list]
                
                sym_groups_dih[3][tuple(symmetry_group_dihedral_list[0][1:3])] = [tuple(element) for element in symmetry_group_dihedral_list]

            elif len(symmetry_group) == 2:

                # angles_to_set[symmetry_group_dihedral_list[0]] = ([0.0, np.pi/3.0])
                angles_to_set[tuple(symmetry_group_dihedral_list[0])] = ([0.0, np.pi/2.0])

                periodicities[tuple(symmetry_group_dihedral_list[0])] = 2
                dihedral_groups[2].extend([tuple(sorted(element, reverse=False)) for element in symmetry_group_dihedral_list])

        return angles_to_set, periodicities, symmetry_group_dihedral_dict, dihedral_groups, sym_groups_dih


    
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
            qm_energy_local = qm_driver.get_energy()
            if hasattr(qm_driver, 'comm') and qm_driver.comm is not None:
                qm_energy = qm_driver.comm.bcast(
                    qm_energy_local if qm_driver.comm.Get_rank() == mpi_master() else None,
                    root=mpi_master(),
                )
            else:
                qm_energy = qm_energy_local

            if qm_energy is None:
                raise RuntimeError('XTB energy is None on this rank after MPI synchronization.')
            qm_energy = np.array([qm_energy])

        # restricted SCF
        elif isinstance(qm_driver, ScfRestrictedDriver):
            qm_driver.ostream.mute()
            scf_results = qm_driver.compute(molecule, basis)
            qm_energy = np.array([qm_driver.scf_energy])
            qm_driver.ostream.unmute()
            qm_driver.filename = None
            qm_driver.checkpoint_file = None

            print('qm_energy in SCF driver', qm_energy)
        
        elif isinstance(qm_driver, ScfUnrestrictedDriver):
            qm_driver.ostream.mute()
            scf_results = qm_driver.compute(molecule, basis)
            qm_energy = np.array([qm_driver.scf_energy])
            qm_driver.ostream.unmute()
            qm_driver.filename = None
            qm_driver.checkpoint_file = None

            print('qm_energy in SCF driver', qm_energy)

        elif isinstance(qm_driver, LinearResponseEigenSolver) or isinstance(qm_driver, TdaEigenSolver):
            self.drivers['es'][0].ostream.mute()
            self.drivers['es'][3].ostream.mute()
            scf_results = self.drivers['es'][3].compute(molecule, basis)
            self.drivers['es'][3].filename = None
            self.drivers['es'][3].checkpoint_file = None
            scf_energy = self.drivers['es'][3].scf_energy
            
            rsp_results = qm_driver.compute(molecule, basis, scf_results)
            qm_driver.ostream.unmute()
            self.drivers['es'][3].ostream.unmute()
            qm_energy = np.insert(rsp_results['eigenvalues'] + scf_energy, 0, scf_energy)


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
            grad_driver.ostream.unmute()
            qm_gradient = grad_driver.gradient
            qm_gradient = np.array([qm_gradient])

        elif isinstance(grad_driver, ScfGradientDriver):
            grad_driver.ostream.mute()
            grad_driver.compute(molecule, basis, scf_results)
            qm_gradient = grad_driver.gradient
            qm_gradient = np.array([qm_gradient])
            grad_driver.ostream.unmute()
        
        elif isinstance(grad_driver, TddftGradientDriver):
            grad_driver.ostream.mute()
            grad_driver._rsp_results = None
            grad_driver.compute(molecule, basis, self.drivers['es'][3], self.drivers['es'][0], rsp_results)
            qm_gradient = grad_driver.gradient
        
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
            qm_hessian_local = hess_driver.hessian
            if hasattr(hess_driver, 'comm') and hess_driver.comm is not None:
                qm_hessian = hess_driver.comm.bcast(
                    qm_hessian_local if hess_driver.comm.Get_rank() == mpi_master() else None,
                    root=mpi_master(),
                )
            else:
                qm_hessian = qm_hessian_local
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

        elif isinstance(hess_driver, TddftHessianDriver):
            roots = self.drivers['es'][1].state_deriv_index

            qm_hessians = []
            for root in roots:
                self.drivers['es'][1].state_deriv_index = root
                hess_driver.compute(molecule, basis)
                qm_hessian = hess_driver.hessian
                qm_hessians.append(qm_hessian)

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
        N1 = N0 + n

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
            row = qm_ds.shape[0]     # current length
            qm_ds.resize((row + 1,)) # grow by 1
            im_ds.resize((row + 1,)) # grow by 1
            na_ds.resize((row + 1,)) # grow by 1
            cf_ds.resize((row + 1,)) # grow by 1
            xyz_ds.resize((row + 1,)) # grow by 1
        


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