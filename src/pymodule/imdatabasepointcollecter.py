#
#                              VELOXCHEM
#         ----------------------------------------------------
#                     An Electronic Structure Code
#
#  Copyright © 2018-2024 by VeloxChem developers. All rights reserved.
#
#  SPDX-License-Identifier: LGPL-3.0-or-later
#
#  This file is part of VeloxChem.
#
#  VeloxChem is free software: you can redistribute it and/or modify it under
#  the terms of the GNU Lesser General Public License as published by the Free
#  Software Foundation, either version 3 of the License, or (at your option)
#  any later version.
#
#  VeloxChem is distributed in the hope that it will be useful, but WITHOUT
#  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
#  License for more details.
#
#  You should have received a copy of the GNU Lesser General Public License
#  along with VeloxChem. If not, see <https://www.gnu.org/licenses/>.

from math import cos

import builtins
from networkx import node_clique_number
from mpi4py import MPI
import numpy as np
from scipy.optimize import minimize
from scipy.optimize import basinhopping


from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components


import scipy
import h5py
import json
import itertools
import re
import os, copy, math

from pathlib import Path
from sys import stdout
import sys
import random
from time import time
import xml.etree.ElementTree as ET
from xml.dom import minidom
from contextlib import redirect_stderr, contextmanager
from io import StringIO

from typing import Callable, Optional, Dict, Any, List
from copy import deepcopy

from .molecule import Molecule
from .veloxchemlib import mpi_master
from. veloxchemlib import hartree_in_kcalpermol, bohr_in_angstrom, hartree_in_kjpermol
from .outputstream import OutputStream
from .errorhandler import assert_msg_critical
from .solvationbuilder import SolvationBuilder
from .optimizationdriver import OptimizationDriver
from .profiler import Profiler

# Drivers
from .scfrestdriver import ScfRestrictedDriver
from .scfunrestdriver import ScfUnrestrictedDriver
from .molecularbasis import MolecularBasis
from .scfgradientdriver import ScfGradientDriver
from .scfhessiandriver import ScfHessianDriver
from .tdaeigensolver import TdaEigenSolver
from .lreigensolver import LinearResponseEigenSolver
from .tddftgradientdriver import TddftGradientDriver
from .tddfthessiandriver import TddftHessianDriver
from .xtbdriver import XtbDriver
from .xtbgradientdriver import XtbGradientDriver
from .xtbhessiandriver import XtbHessianDriver
from .interpolationdriver import InterpolationDriver
from .interpolationdatapoint import InterpolationDatapoint
from .alphaoptimizer import AlphaOptimizer

from .rotorclass import (
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

try:
    import openmm as mm
    import openmm.app as app
    import openmm.unit as unit
    from openmm.app.metadynamics import Metadynamics, BiasVariable
except ImportError:
    pass

from .trust_radius_grouping_analysis import (improved_group_by_connected_components, 
    explain_grouping_result_for_non_experts,
    TrustRadiusOptimizationAnalyzer,
    summarize_report_for_console)

def _rank_root_print(*args, **kwargs):
    """
    Print only on MPI root rank to avoid delayed/duplicated worker stdout noise.
    Falls back to normal print when MPI rank detection is unavailable.
    """
    try:
        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            builtins.print(*args, **kwargs)
    except Exception:
        builtins.print(*args, **kwargs)


print = _rank_root_print

def _collective_call_inline_style(comm, rank, root, phase_name, fn, *args, **kwargs):
    comm.barrier()
    local_err = None
    result = None
    try:
        result = fn(*args, **kwargs)
    except Exception as exc:
        local_err = f"rank {rank}: {exc}"

    errs = [e for e in comm.allgather(local_err) if e is not None]
    if errs:
        if rank == root:
            raise RuntimeError(f"{phase_name} failed\n" + "\n".join(errs))
        raise RuntimeError(f"{phase_name} failed on another rank")

    comm.barrier()
    return result


class IMDatabasePointCollecter:
    """
    Implements the OpenMMDynamics.
    Performs classical MD and QM/MM simulations.

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.

    Instance variables:
        - platform: The platform for OpenMM. Default is 'reference'
                    Options are 'CPU', 'CUDA' and 'OpenCL'.
        - ensemble: The thermodynamic ensemble used in the simulation.
                    Options are: 'NVE', 'NVT' and 'NPT'.
        - temperature: The temperature.
        - friction: The friction coefficient which couples the system to the heat bath.
        - timestep: The timestep for the integrator in fs.
        - nsteps: The number of steps.
        - parent_ff: The Force Field to be used in standard residues.
        - water_ff: The Force Field to be used for water.
        - box_size: The size of the box for periodic simulations.
        - padding: The padding for the box in nm.
        - cutoff: The cutoff for nonbonded interactions in nm.
        - integrator: The integrator object.
        - system: The OpenMM system object.
        - topology: The OpenMM topology object.
        - simulation: The OpenMM simulation object.
        - constraints: The constraints for the system.
        - pdb: The PDB file object.
        - modeller: The modeller object.
        - labels: The atom labels.
        - molecule: The VeloxChem molecule object.
        - unique_residues: The list of unique residues in the system.
        - unique_molecules: The list of unique molecules in the system.
        - qm_driver: The VeloxChem driver object. Options are XtbDriver and InterpolationDriver.
        - grad_driver: The VeloxChem gradient driver object.
        - qm_atoms: The list of atom indices for the QM region.
        - mm_subregion: The list of atom indices for the MM subregion (part of a molecule).
        - linking_atoms: The list of atom indices for the linking atoms.
        - qm_force_index: The index of the QM force in the system.
        - driver_flag: The flag to determine the driver to be used.
        - scaling_factor: The scaling factor for the QM stabilizer.
        - linking_atom_distance: The distance between the QM and MM regions in angstroms.
    """
    
    def __init__(self, comm=None, ostream=None):
        """
        Initializes the class with default simulation parameters.
        """
        np.set_printoptions(threshold=sys.maxsize)
        # MPI and output stream
        if comm is None:
            comm = MPI.COMM_WORLD

        if ostream is None:
            if comm.Get_rank() == mpi_master():
                ostream = OutputStream(stdout)
            else:
                ostream = OutputStream(None)

        # MPI communicators
        # self.comm remains the external communicator used by legacy code paths.
        # self._mpi_ctrl_comm is a dedicated duplicate used by the root-worker
        # control plane (STEP/STOP).
        # self._mpi_interp_comm is a separate duplicate used by interpolation
        # preload collectives. Keeping them separate prevents cross-talk between
        # control messages and interpolation payload broadcasts.
        self.comm = comm
        try:
            self._mpi_ctrl_comm = self.comm.Dup()
        except Exception:
            self._mpi_ctrl_comm = self.comm

        try:
            self._mpi_interp_comm = self.comm.Dup()
        except Exception:
            self._mpi_interp_comm = self._mpi_ctrl_comm

        # MPI information (based on control communicator)
        self.rank = self._mpi_ctrl_comm.Get_rank()
        self.nodes = self._mpi_ctrl_comm.Get_size()

        # output stream
        self.ostream = ostream

        # Instance variables
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
        self.box_size = 2.0 
        self.padding = 1.0
        self.cutoff = 1.0
        self.scaling_factor = 0.0
        self.integrator = None
        self.bias_force_reaction_idx = None
        self.bias_force_reaction_prop = None
        self.global_theta_list = []

        # meta dynamics settings
        self.metadynamics_settings = None
        self.metadynamics_enabled = None
        self._metadynamics = None
        self._metadynamics_forcegroup = 11
        self.metadynamics_variables = []
        self.metadynamics_cv_forces = []

        self.metadynamics_bias_energies = []
        self.metadynamics_cv_history = []

        # OpenMM objects
        self.system = None
        self.topology = None
        self.integrator = None
        self.simulation = None
        self.constraints = None
        self.pdb = None
        self.modeller = None
        self.labels = []

        # VeloxChem objects
        self.molecule = None
        self.energy_gabs = {}
        self.state_energies = {}
        self.unique_residues = []
        self.unique_molecules = []

        self.current_im_choice = None
        
        # QM Region parameters
        self.drivers = None
        self.qm_atoms = None
        self.mm_subregion = None
        self.linking_atoms = None
        self.qm_force_index = None
        self.driver_flag = None
        self.swapped = False
        self.state_swtiched = False
        self.velo_switch = False
        self.excitation_pulse = None
        self.prev_dens_of_points = None

        self.roots_to_follow = []

        self.snapshots = None
        self.load_system = None
        self.pressure = None
        self.output_file = None
        self.adiabatic_basis = False
        self.density_around_data_point = None
        self.impes_drivers = None
        self.im_labels = None
        self.sorted_im_labels = []
        self.qm_energies = None
        self.qm_data_points = None
        self.point_adding_molecule = {}
        self.energy_threshold = None
        self.gradient_rmsd_thrsh = None
        self.force_orient_thrsh = 0.00001
        self.start_collect = 0
        self.previous_energy_list = []
        # state-resolved symmetry metadata excluding the invariant core atoms
        self.non_core_symmetry_groups = []
        self.interpolation_settings = None
        self.allowed_molecules = None
        self.reference_struc_energies_file = None
        # rotatable-bond definitions used to identify methyl-like rotor groups
        self.all_rot_bonds = None

        self.starting_temperature = None

        self.current_state = None
        self.starting_state = 0
        self.distance_thrsh = 0.1
        self.current_im_choice = None
        self.current_gradient = 0
        self.current_energy = 0
        self.point_checker = 1
        self.last_added = 0
        self.allowed_molecule_deviation = None
        self.last_point_added = None
        self.im_labels = None
        self.add_a_point = False
        self.check_a_point = False
        # if enabled, generates and stores explicit rotor-cluster state banks
        self.cluster_run = None
        self.skipping_value = 0
        self.basis_set_label = None
        self.molecule = None
        self.expansion = True
        self.expansion_molecules = []
        self.dynamics_settings_interpolation_run = None
        self.sampled_molecules = None
        self.bayes_models = None

        self.opt_mols_org_mol_swap = {}

        self.sampling_enabled = False
        self.sampling_settings = None
        self.sampling_driver = None
        self.sampling_impes_drivers = None
        self.sampling_interpolation_settings = None
        self.sampling_qm_datapoints_dict = None
        # sampling-mode mirror of the symmetry-expanded datapoint registry
        self.sampling_qm_symmetry_datapoint_dict = None

        # symmetry considerations
        # per-root-state active rotor cluster information used for setup diagnostics
        self.rotor_cluster_information = {"gs": None, "es": None}
        # root -> family_label -> cluster bank payload loaded from registry
        self.qm_rotor_cluster_banks = None
        # rotor_id -> RotorDefinition mapping for currently analyzed symmetry rotors
        self.symmetry_rotors = None
        self.cluster_run = None
        

        # output_file variables that will be written into self.general_variable_output
        self.gradients = None
        self.start_velocities = None
        self.coordinates = None
        self.coordinates_xyz = None
        self.state_specific_molecules = None
        self.all_gradients = []
        self.inv_sqrt_masses = None
        
        # toggle + optional fixed list for important internal-coordinate constraints
        self.identfy_relevant_int_coordinates = (True, [])
        # master switch for generating/using symmetry-equivalent datapoints
        self.use_symmetry = True
        self.use_opt_confidence_radius = True
        self.confidence_radius_optimized = True
        # root-indexed dictionary of symmetry-expanded QM datapoint banks
        self.qm_symmetry_dict = None
        self.add_gpr_model = True
        
        self.summary_output = 'summary_output.h5'
        self.coordinates_xyz = None

        self.velocities = []

        # Default value for the C-H linker distance
        self.linking_atom_distance = 1.0705 

        self.eq_bond_force_constants = None

        # Runtime profiling controls (optional, default off).
        self.profile_runtime_timing = False
        self.profile_runtime_print_interval = 0
        self.profile_interpolation_timing = False
        self.profile_interpolation_print_summary = True
        self._qmmm_runtime_profiler = None
        self._qmmm_runtime_totals = {}
        self._qmmm_runtime_step = {}
        self._qmmm_runtime_steps = 0

        self.mpi_control_plane_enabled = (self.nodes > 1)
        self.mpi_root_worker_mode = (self.nodes > 1)
        self.mpi_reload_from_hdf5 = True
        self.mpi_debug_sync = False

        # command IDs
        self._MPI_CMD_STEP = 1
        self._MPI_CMD_STOP = 2
        self._MPI_CMD_AUX = 3

        # AUX task IDs
        self._MPI_AUX_EG = 11        # energy + gradient
        self._MPI_AUX_EGH = 12       # energy + gradient + hessian
        self._MPI_AUX_OPT = 13       # constrained optimization

        # True only while worker service loop is active
        self._mpi_worker_service_running = False

        # roots requiring dataset reload on worker ranks
        self._mpi_pending_sync_roots = set()

        # Cache of the most recent QM molecule passed through interpolation.
        self._latest_qm_molecule = None
        self._latest_qm_positions_nm = None

        self._test_hooks: Dict[str, Callable[[Dict[str, Any]], None]] = {}
        self._collect_step_trace: bool = False
        self._step_trace: List[Dict[str, Any]] = []
        self._strict_test_hooks: bool = False
        

    def _create_platform(self):
        """
        Creates an OpenMM platform.

        Returns:
            OpenMM Platform.
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
        
    def _set_test_hooks(
            self, 
            hooks: Optional[Dict[str, Callable[[Dict[str, Any]], None]]] = None,
            collect_step_trace: bool = False,
            strict: bool = False, 
        ):

        self._test_hooks = hooks or {}
        self._collect_step_trace = bool(collect_step_trace)
        self._strict_test_hooks = bool(strict)
        self._step_trace = []
    
    def _emit_test_hook(self, event_name: str, payload: Dict[str, Any]):
        fn = self._test_hooks.get(event_name)
        if fn is None:
            if self._collect_step_trace and event_name == "step":
                self._step_trace.append(deepcopy(payload))
            return

        try:
            fn(payload)
        except Exception as exc:
            if self._strict_test_hooks:
                raise RuntimeError(
                    f"Test hook '{event_name}' failed: {exc}"
                ) from exc

        if self._collect_step_trace and event_name == "step":
            self._step_trace.append(deepcopy(payload))

    def get_step_trace(self):
        return list(self._step_trace)

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
        Generates an OpenMM system from a VeloxChem molecule and a forcefield generator.
        
        :param molecule:
            VeloxChem molecule object.
        :param ff_gen:
            VeloxChem forcefield generator object.
        :param solvent:
            Available options:'gas', 'spce', 'tip3p', 'ethanol', 'methanol', 'acetone', 
            'chloroform', 'hexane', 'toluene', 'dcm', 'benzene', 'dmso', 'thf', 
            'acetonitrile', 'other' or 'itself'.
        :param qm_atoms:
            Options: None, 'all', or list of atom indices for QM region.
        :param filename:
            Base name for the generated files. Default is 'residue'.
        :param residue_name:
            Name of the residue. Default is 'MOL'.
        """

        assert_msg_critical('openmm' in sys.modules, 'OpenMM is required for IMDatabasePointCollecter.')

        # Store the molecule object and generate OpenMM compatible files
        self.molecule = molecule
        self.positions = molecule.get_coordinates_in_angstrom()
        self.labels = molecule.get_labels()
        self.root_z_matrix = z_matrix
        has_multiple_ranks = (self.nodes > 1)
        do_io_on_rank = (not has_multiple_ranks) or self._mpi_is_root()

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
                                        filename,
                                        write_xml=do_io_on_rank)

            elif isinstance (qm_atoms, list):
                if qm_atoms == list(range(len(self.labels))):
                    msg = 'Full molecule as QM region'
                    self.ostream.print_info(msg)
                    self.ostream.flush()
                    self._create_QM_residue(ff_gen,
                                            qm_atoms,
                                            filename,
                                            write_xml=do_io_on_rank)
                msg = 'QM/MM partition inside the molecule'
                self.ostream.print_info(msg)
                self.ostream.flush()
                self._create_QM_subregion(ff_gen,
                                          qm_atoms, 
                                          molecule, 
                                          filename,
                                          write_xml=do_io_on_rank)

            if do_io_on_rank:
                ff_gen.write_pdb(f'{filename}.pdb', 'QMR')

        elif qm_atoms is None:
            if do_io_on_rank:
                ff_gen.write_openmm_files(filename, residue_name)

        if has_multiple_ranks:
            self._mpi_ctrl_comm.Barrier()

        # TODO: Incorporate the SystemBuilder here to avoid using the Modeller.
        if solvent == 'gas':
            phase = 'gas'
            self.pdb = app.PDBFile(f'{filename}.pdb')     
            # Common forcefield loading, modified according to phase specifics
            forcefield_files = [f'{filename}.xml']

        if solvent != 'gas':
            # Solvate the molecule using the SystemBuilder
            phase = 'periodic'
            if do_io_on_rank:
                sys_builder = SolvationBuilder()
                sys_builder.solvate(solute=molecule, 
                                    solvent=solvent,
                                    padding=self.padding,
                                    equilibrate=True,
                                    constraint_solute=True
                                    )
                sys_builder.write_openmm_files(solute_ff=ff_gen)

            if has_multiple_ranks:
                self._mpi_ctrl_comm.Barrier()

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
                print('I am going in here', forcefield_files)

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
            # self.qm_stabilizer(ff_gen)
        
        # Write the system to a xml file (for debugging purposes)
        if do_io_on_rank:
            with open(f'{filename}_system.xml', 'w') as f:
                f.write(mm.XmlSerializer.serialize(self.system))
                msg = f'System parameters written to {filename}_system.xml'
                self.ostream.print_info(msg)
                self.ostream.flush()

        # Write the system to a pdb file (for debugging purposes)
        if do_io_on_rank:
            if phase == 'gas':
                app.PDBFile.writeFile(self.pdb.topology, 
                                      self.positions, 
                                      open(f'{filename}_system.pdb', 'w'))
                msg = f'System coordinates written to {filename}_system.pdb'
                self.ostream.print_info(msg)
                self.ostream.flush()

            elif phase == 'periodic':
                app.PDBFile.writeFile(self.pdb.topology, 
                                      self.pdb.positions, 
                                      open(f'{filename}_system.pdb', 'w'))
                msg = f'System coordinates written to {filename}_system.pdb'
                self.ostream.print_info(msg)
                self.ostream.flush()
        
        self.phase = phase

    def _refresh_interpolation_driver_caches(self, root):
        """
        Refreshes both legacy runtime cache and MPI preload cache for one root.
        Call this after any datapoint/symmetry update.
        """
        driver = self.impes_drivers[root]
        driver.mark_runtime_data_cache_dirty()
        driver.prepare_runtime_data_cache(force=True)
        if getattr(driver, 'mpi_engine', None) is not None:
            driver.mpi_engine.mark_dirty()
            # In root-worker mode, defer preload rebuild until synchronized compute
            # to avoid standalone collectives desynchronizing the worker control loop.
            if not (self._mpi_is_active() and self.mpi_root_worker_mode):
                driver.mpi_engine.preload_static_data(force=True)
        
        if self.sampling_enabled:
            samp_driver = self.sampling_impes_drivers[root]
            samp_driver.mark_runtime_data_cache_dirty()
            samp_driver.prepare_runtime_data_cache(force=True)
            if getattr(driver, 'mpi_engine', None) is not None:
                samp_driver.mpi_engine.mark_dirty()
                # In root-worker mode, defer preload rebuild until synchronized compute
                # to avoid standalone collectives desynchronizing the worker control loop.
                if not (self._mpi_is_active() and self.mpi_root_worker_mode):
                    samp_driver.mpi_engine.preload_static_data(force=True)
            
        self._mpi_mark_root_dirty(root)
    
    def _mpi_is_active(self):
        return bool(self.nodes > 1 and self.mpi_control_plane_enabled)

    def _mpi_is_root(self):
        return self.rank == mpi_master()

    def _mpi_mark_root_dirty(self, root):
        if self._mpi_is_active() and self.mpi_root_worker_mode:
            self._mpi_pending_sync_roots.add(int(root))

    def _mpi_bcast_control(self, message):
        if not self._mpi_is_active():
            return message
        return self._mpi_ctrl_comm.bcast(message if self._mpi_is_root() else None, root=mpi_master())

    def _mpi_decode_cmd(self, cmd_obj):
        if isinstance(cmd_obj, np.ndarray):
            arr = np.asarray(cmd_obj)
            if arr.size == 1:
                return int(arr.reshape(-1)[0])
            raise RuntimeError(
                f"MPI control-channel desync on rank {self.rank}: expected scalar command, "
                f"got ndarray shape={arr.shape} dtype={arr.dtype}")

        if isinstance(cmd_obj, (int, np.integer)):
            return int(cmd_obj)

        raise RuntimeError(
            f"MPI control-channel desync on rank {self.rank}: expected int command, "
            f"got {type(cmd_obj).__name__}")

    @staticmethod
    def _serialize_molecule_for_mpi(mol):
        return (mol.get_xyz_string(), int(mol.get_charge()), int(mol.get_multiplicity()))

    @staticmethod
    def _deserialize_molecule_from_mpi(payload):
        mol = Molecule.from_xyz_string(payload[0])
        mol.set_charge(int(payload[1]))
        mol.set_multiplicity(int(payload[2]))
        return mol

    def _mpi_can_dispatch_aux(self):
        return bool(
            self._mpi_is_active()
            and self.mpi_root_worker_mode
            and self._mpi_is_root()
            and self._mpi_worker_service_running
        )

    def _mpi_run_aux(self, aux_cmd, payload):
        # Root-worker mode: root dispatches AUX command to parked workers.
        if self._mpi_can_dispatch_aux():
            self._mpi_bcast_control({
                'cmd': self._MPI_CMD_AUX,
                'aux_cmd': int(aux_cmd),
                'payload': payload,
            })
            return self._mpi_execute_aux_collective(aux_cmd, payload)

        # Full-SPMD mode (mpi_root_worker_mode=False): all ranks enter here directly.
        if self._mpi_is_active() and (not self.mpi_root_worker_mode):
            return self._mpi_execute_aux_collective(aux_cmd, payload)

        # Root-local fallback.
        return self._mpi_execute_aux_local(aux_cmd, payload)

    def _mpi_execute_aux_collective(self, aux_cmd, payload):
        if int(aux_cmd) == self._MPI_AUX_EG:
            return self._mpi_aux_eval_qm_block(payload, include_hessian=False, collective_mode=True)
        if int(aux_cmd) == self._MPI_AUX_EGH:
            return self._mpi_aux_eval_qm_block(payload, include_hessian=True, collective_mode=True)
        if int(aux_cmd) == self._MPI_AUX_OPT:
            return self._mpi_aux_optimize_block(payload, collective_mode=True)
        raise RuntimeError(f"Unknown MPI AUX command {aux_cmd}")

    def _mpi_execute_aux_local(self, aux_cmd, payload):
        if int(aux_cmd) == self._MPI_AUX_EG:
            return self._mpi_aux_eval_qm_block(payload, include_hessian=False, collective_mode=False)
        if int(aux_cmd) == self._MPI_AUX_EGH:
            return self._mpi_aux_eval_qm_block(payload, include_hessian=True, collective_mode=False)
        if int(aux_cmd) == self._MPI_AUX_OPT:
            return self._mpi_aux_optimize_block(payload, collective_mode=False)
        raise RuntimeError(f"Unknown MPI AUX command {aux_cmd}")
    
    def _mpi_aux_eval_qm_block(self, payload, include_hessian=False, collective_mode=False):
        driver_key = str(payload['driver_key'])   # 'gs' or 'es'
        state_mask = [int(x) for x in payload.get('state_mask', [])]
        molecule = self._deserialize_molecule_from_mpi(payload['molecule'])
        basis_label = str(payload['basis_label'])
        basis = MolecularBasis.read(molecule, basis_label)

        drivers = self.drivers[driver_key]
        restore_state_deriv = None

        if isinstance(drivers[0], (LinearResponseEigenSolver, TdaEigenSolver)):
            restore_state_deriv = list(drivers[1].state_deriv_index)
            drivers[1].state_deriv_index = list(state_mask)

        try:
            energies, scf_results, rsp_results = self._compute_energy_mpi_safe(
                drivers[0], molecule, basis,
                collective=collective_mode,
                phase_name='AUX energy',
            )

            if isinstance(drivers[0], (LinearResponseEigenSolver, TdaEigenSolver)):
                energies = np.asarray(energies)[state_mask]

            gradients = self._compute_gradient_mpi_safe(
                drivers[1], molecule, basis, scf_results, rsp_results,
                collective=collective_mode,
                phase_name='AUX gradient',
            )

            out = {
                'energies': energies,
                'gradients': gradients,
            }

            if include_hessian:
                hessians = self._compute_hessian_mpi_safe(
                    drivers[2], molecule, basis,
                    collective=collective_mode,
                    phase_name='AUX hessian',
                )
                out['hessians'] = hessians

            return out
        finally:
            if restore_state_deriv is not None:
                drivers[1].state_deriv_index = restore_state_deriv

    def _mpi_aux_optimize_block(self, payload, collective_mode=False):
        state_to_optim = int(payload['state_to_optim'])
        molecule = self._deserialize_molecule_from_mpi(payload['molecule'])
        constraints = payload.get('constraints', [])
        basis_label = str(payload['basis_label'])
        basis = MolecularBasis.read(molecule, basis_label)

        if state_to_optim == 0:
            if isinstance(self.drivers['gs'][0], (ScfRestrictedDriver, ScfUnrestrictedDriver)):
                _, scf_results, _ = self._compute_energy_mpi_safe(
                    self.drivers['gs'][0], molecule, basis,
                    collective=collective_mode, phase_name='AUX opt pre-energy')
                optimized_molecule, opt_results = self._run_optimization_mpi_safe(
                    self.drivers['gs'][0],
                    molecule,
                    constraints=constraints,
                    index_offset=1,
                    compute_args=(basis, scf_results),
                    source_molecule=molecule,
                    collective=collective_mode,
                    phase_name='AUX optimization',
                )
            elif isinstance(self.drivers['gs'][0], XtbDriver):
                optimized_molecule, opt_results = self._run_optimization_mpi_safe(
                    self.drivers['gs'][0],
                    molecule,
                    constraints=constraints,
                    index_offset=1,
                    source_molecule=molecule,
                    collective=collective_mode,
                    phase_name='AUX optimization',
                )
            else:
                raise RuntimeError("Unsupported GS optimization driver in AUX optimization.")
        else:
            if not isinstance(self.drivers['es'][0], (LinearResponseEigenSolver, TdaEigenSolver)):
                raise RuntimeError("Unsupported ES optimization driver in AUX optimization.")

            old_state_deriv = list(self.drivers['es'][1].state_deriv_index)
            self.drivers['es'][1].state_deriv_index = [state_to_optim]
            try:
                _, _, rsp_results = self._compute_energy_mpi_safe(
                    self.drivers['es'][0], molecule, basis,
                    collective=collective_mode, phase_name='AUX ES pre-energy')

                optimized_molecule, opt_results = self._run_optimization_mpi_safe(
                    self.drivers['es'][1],
                    molecule,
                    constraints=constraints,
                    index_offset=1,
                    compute_args=(basis, self.drivers['es'][3], self.drivers['es'][0], rsp_results),
                    source_molecule=molecule,
                    collective=collective_mode,
                    phase_name='AUX ES optimization',
                )
            finally:
                self.drivers['es'][1].state_deriv_index = old_state_deriv

        return {
            'optimized_molecule': self._serialize_molecule_for_mpi(optimized_molecule),
            'opt_results': opt_results,
        }
    
    def _aux_eval_eg(self, molecule, basis_label, driver_key, state_mask):
        payload = {
            'molecule': self._serialize_molecule_for_mpi(molecule),
            'basis_label': str(basis_label),
            'driver_key': str(driver_key),
            'state_mask': [int(x) for x in state_mask],
        }
        out = self._mpi_run_aux(self._MPI_AUX_EG, payload)
        return out['energies'], out['gradients']

    def _aux_eval_egh(self, molecule, basis_label, driver_key, state_mask):
        payload = {
            'molecule': self._serialize_molecule_for_mpi(molecule),
            'basis_label': str(basis_label),
            'driver_key': str(driver_key),
            'state_mask': [int(x) for x in state_mask],
        }
        out = self._mpi_run_aux(self._MPI_AUX_EGH, payload)
        return out['energies'], out['gradients'], out['hessians']

    def _aux_optimize(self, state_to_optim, molecule, basis_label, constraints):
        payload = {
            'state_to_optim': int(state_to_optim),
            'molecule': self._serialize_molecule_for_mpi(molecule),
            'basis_label': str(basis_label),
            'constraints': constraints,
        }
        out = self._mpi_run_aux(self._MPI_AUX_OPT, payload)
        optimized_molecule = self._deserialize_molecule_from_mpi(out['optimized_molecule'])
        return optimized_molecule, out['opt_results']


    def _run_qmmm_worker_service(self):
        """
        Worker-side MPI service loop.
        Receives STEP messages from rank 0 and executes interpolation compute only.
        """
        expected_step_idx = 0
        while True:
            msg = self._mpi_bcast_control(None)
            if not isinstance(msg, dict):
                raise RuntimeError(
                    f"MPI control-channel desync on worker rank {self.rank}: "
                    f"unexpected message type={type(msg).__name__}")
            if 'cmd' not in msg:
                raise RuntimeError(
                    f"MPI control-channel desync on worker rank {self.rank}: "
                    f"control dict missing 'cmd'; keys={list(msg.keys())}")

            cmd = self._mpi_decode_cmd(msg['cmd'])

            if cmd == self._MPI_CMD_STOP:
                break

            if cmd == self._MPI_CMD_STEP:
                step_idx_obj = msg.get('step_idx', None)
                if step_idx_obj is not None:
                    step_idx = self._mpi_decode_cmd(step_idx_obj)
                    if step_idx != expected_step_idx:
                        raise RuntimeError(
                            f"MPI STEP index mismatch on worker rank {self.rank}: "
                            f"expected {expected_step_idx}, got {step_idx}")
                    expected_step_idx += 1

                if 'qm_positions_nm' not in msg:
                    raise RuntimeError(
                        f"MPI STEP message missing 'qm_positions_nm' on worker rank {self.rank}")

                sync_roots = msg.get('sync_roots', [])
                qm_positions_nm = msg['qm_positions_nm']
                self._mpi_worker_compute_step(qm_positions_nm, sync_roots)
                continue

            if cmd == self._MPI_CMD_AUX:
                aux_cmd_obj = msg.get('aux_cmd', None)
                if aux_cmd_obj is None:
                    raise RuntimeError(
                        f"MPI AUX message missing 'aux_cmd' on worker rank {self.rank}")
                aux_cmd = self._mpi_decode_cmd(aux_cmd_obj)
                payload = msg.get('payload', {})
                self._mpi_execute_aux_collective(aux_cmd, payload)
                continue

            raise RuntimeError(f"Unknown MPI command {cmd} on worker rank {self.rank}")

    def _mpi_collect_driver_objects_for_root_local_work(self):
        objs = []
        drv_dict = getattr(self, 'drivers', None)
        if isinstance(drv_dict, dict):
            for key in ('gs', 'es'):
                drv_val = drv_dict.get(key)
                if isinstance(drv_val, (list, tuple)):
                    for obj in drv_val:
                        if obj is not None:
                            objs.append(obj)
                elif drv_val is not None:
                    objs.append(drv_val)
        return objs

    def _mpi_apply_comm_override(self, obj, comm, saved, seen):
        if obj is None:
            return

        obj_id = id(obj)
        if obj_id in seen:
            return
        seen.add(obj_id)

        def _record_and_set(target, attr_name, new_value):
            try:
                old_value = getattr(target, attr_name)
            except Exception:
                return
            try:
                setattr(target, attr_name, new_value)
            except Exception:
                return
            saved.append((target, attr_name, old_value))

        # Update both public and private communicator/rank/size fields.
        for comm_attr in ('comm', '_comm'):
            if hasattr(obj, comm_attr):
                _record_and_set(obj, comm_attr, comm)

        for rank_attr in ('rank', '_rank'):
            if hasattr(obj, rank_attr):
                _record_and_set(obj, rank_attr, comm.Get_rank())

        for size_attr in ('nodes', '_nodes', 'size', '_size', 'nnodes', '_nnodes'):
            if hasattr(obj, size_attr):
                _record_and_set(obj, size_attr, comm.Get_size())

        # Recurse through nested veloxchem objects and container members.
        if isinstance(obj, dict):
            iterable = obj.values()
        elif isinstance(obj, (list, tuple, set)):
            iterable = obj
        else:
            iterable = getattr(obj, '__dict__', {}).values()

        for child in iterable:
            if child is None:
                continue
            if isinstance(child, (str, bytes, int, float, complex, bool, np.ndarray)):
                continue
            child_mod = getattr(getattr(child, '__class__', None), '__module__', '')
            if child_mod.startswith('veloxchem') or isinstance(child, (dict, list, tuple, set)):
                self._mpi_apply_comm_override(child, comm, saved, seen)

    @staticmethod
    def _mpi_restore_comm_override(saved):
        for target, attr_name, old_value in reversed(saved):
            try:
                setattr(target, attr_name, old_value)
            except Exception:
                pass

    @contextmanager
    def _mpi_root_local_qm_context(self, extra_objects=None):
        """
        Temporarily forces auxiliary QM drivers to run on COMM_SELF.

        This avoids deadlocks when rank 0 performs point-check/add-point
        QM/optimization work while worker ranks are blocked in the
        root-worker control service loop.
        """

        if not (self._mpi_is_active() and self.mpi_root_worker_mode and self._mpi_is_root()):
            yield
            return

        saved = []
        seen = set()
        local_comm = MPI.COMM_SELF

        for obj in self._mpi_collect_driver_objects_for_root_local_work():
            self._mpi_apply_comm_override(obj, local_comm, saved, seen)

        if extra_objects is not None:
            for obj in extra_objects:
                self._mpi_apply_comm_override(obj, local_comm, saved, seen)

        try:
            yield
        finally:
            self._mpi_restore_comm_override(saved)

    def _opt_compute_mpi_safe(self, opt_drv, molecule, *compute_args):
        with self._mpi_root_local_qm_context(extra_objects=[opt_drv]):
            if self._mpi_is_active() and self.mpi_root_worker_mode and self._mpi_is_root() and hasattr(opt_drv, "comm"):
                if opt_drv.comm.Get_size() != 1:
                    raise RuntimeError("Root-local MPI override failed for OptimizationDriver (comm size != 1).")
            return opt_drv.compute(molecule, *compute_args)


    def _ensure_interp_engine_comm(self, driver):
        engine = getattr(driver, 'mpi_engine', None)
        if engine is None:
            return

        if getattr(engine, 'comm', None) is self._mpi_interp_comm:
            return

        engine.comm = self._mpi_interp_comm
        engine.rank = self._mpi_interp_comm.Get_rank()
        engine.size = self._mpi_interp_comm.Get_size()
    
    def _interp_compute_root_local_serial(self, driver, molecule):
        """
        Runs InterpolationDriver.compute in a root-local serial-safe mode.

        In root-worker MPI mode, this is used for auxiliary/root-only logic
        (point checks, constrained point generation) where workers are parked
        in the control loop. MPI preload collectives must be disabled here.
        """
        if not (self._mpi_is_active() and self.mpi_root_worker_mode and self._mpi_is_root()):
            self._ensure_interp_engine_comm(driver)
            driver.compute(molecule)
            return

        self._ensure_interp_engine_comm(driver)
        old_collective = bool(getattr(driver, 'mpi_collective_compute_active', False))
        old_preload_enabled = bool(getattr(driver, 'mpi_preload_enabled', False))

        try:
            driver.mpi_collective_compute_active = False
            driver.mpi_preload_enabled = False
            driver.compute(molecule)
        finally:
            driver.mpi_collective_compute_active = old_collective
            driver.mpi_preload_enabled = old_preload_enabled

    def add_bias_force(self, atoms, force_constant, target=0.109, anchor=False):
        """
        Method to add a biasing force to the system.

        :param atoms:
            List of atom indices to apply the force.
        :param force_constant:
            Force constant for the biasing force.
        :param target:
            Target equilibrium parameter (distance (nm), angle and torsion (deg))
        """

        assert_msg_critical('openmm' in sys.modules, 'OpenMM is required for IMDatabasePointCollecter.')

        if len(atoms) == 2 and anchor:
            msg = f'Adding stretch force between atom {atoms[0]} and anchor {atoms[1]} with force constant {force_constant}.'
            self.ostream.print_info(msg)
            self.ostream.flush()
            for force in self.system.getForces():
                if isinstance(force, mm.NonbondedForce):
                    force.addParticle(0.0, 1.0, 0.0)  

            force = mm.CustomBondForce("0.5*k*(max(0, r-r0))^2")
            force.addGlobalParameter("k", force_constant)
            force.addGlobalParameter("r0", 0.109)  # nm
            force.addBond(atoms[0], atoms[1])
            self.system.addForce(force)
        
        elif len(atoms) == 2:
            msg = f'Adding stretch force between atoms {atoms[0]} and {atoms[1]} with force constant {force_constant}.'
            self.ostream.print_info(msg)
            self.ostream.flush()
            force = mm.CustomBondForce('0.5*k*(r-r0)^2')
            force.addGlobalParameter('k', force_constant)
            force.addGlobalParameter('r0', target)
            force.addBond(atoms[0], atoms[1])
            self.system.addForce(force)
        elif len(atoms) == 3:
            target_rad = target * np.pi / 180
            msg = f'Adding bend force between atoms {atoms[0]}, {atoms[1]}, and {atoms[2]} with force constant {force_constant}.'
            self.ostream.print_info(msg)
            self.ostream.flush()
            force = mm.CustomAngleForce('0.5*k*(theta-theta0)^2')
            force.addGlobalParameter('k', force_constant)
            force.addGlobalParameter('theta0', target_rad)
            force.addAngle(atoms[0], atoms[1], atoms[2])
            self.system.addForce(force)
        elif len(atoms) == 4:
            target_rad = target * np.pi / 180
            msg = f'Adding torsion force between atoms {atoms[0] - 1}, {atoms[1] - 1}, {atoms[2] -1}, and {atoms[3] - 1} with force constant {force_constant}.'        
            self.ostream.print_info(msg)
            self.ostream.flush()
            force = mm.CustomTorsionForce('k*2*(1 - cos(theta - theta0))')
            # force = mm.CustomTorsionForce('0.5*k*(theta-theta0)^2')
            force.addGlobalParameter('k', force_constant)
            force.addGlobalParameter('theta0', target_rad)
            force.addTorsion(atoms[0] - 1, atoms[1] - 1, atoms[2] - 1, atoms[3] - 1)
            self.system.addForce(force)
        else:
            raise ValueError('Invalid number of atoms for the biasing force.')
     
    def _build_metadynamics_cv_force(self, spec):
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

        if not self.metadynamics_enabled:
            self._metadynamics = None
            self.metadynamics_variables = []
            self.metadynamics_cv_forces = []

            return

        cfg = self.metadynamics_settings

        self._metadynamics_forcegroup = int(cfg.get('force_group', 11))

        variables = []
        cv_forces = []

        for spec in cfg['variables']:

            cv_force, min_val, max_val, width, periodic = self._build_metadynamics_cv_force(spec)

            cv_forces.append(cv_force)
            variables.append(BiasVariable(cv_force, min_val, max_val, width, periodic))

        self.metadynamics_cv_forces = cv_forces
        self.metadynamics_variables = variables

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
        for i in range(self.system.getNumForces()):
            f = self.system.getForce(i)
            print(i, type(f).__name__, f.getForceGroup(), f.getName())

        # Usually the newly added metadynamics bias is the last force
        bias_force = self.system.getForce(self.system.getNumForces() - 1)
        bias_force.setForceGroup(self._metadynamics_forcegroup)
    
    def _reset_metadynamics_bias(self, context, reason='datapoint_added'):
        """
        Soft reset only:
        clear accumulated metadynamics hills and restart bias growth from zero.
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

        # Reset local history containers as well.
        self.metad_bias_energies = []
        self.metad_cv_history = []
        self.metadynamics_bias_energies = []
        self.metadynamics_cv_history = []


        print(
            f"[metadynamics] soft bias reset done (reason={reason}, "
        )
        return True


    def conformational_sampling(self, 
                                ensemble='NVT', 
                                temperature=700, 
                                timestep=2.0, 
                                nsteps=10000, 
                                snapshots=10,
                                minimize=True):

        """
        Runs a high-temperature MD simulation to sample conformations and minimize the energy of these conformations.

        :param ensemble:
            Type of ensemble. Options are 'NVE', 'NVT', 'NPT'. Default is 'NVT'.
        :param temperature:
            Temperature of the system in Kelvin. Default is 700 K.
        :param timestep:
            Timestep of the simulation in femtoseconds. Default is 2.0 fs.
        :param nsteps:
            Number of steps in the simulation. Default is 1000.
        :param snapshots:
            The number of snapshots to save. Default is 10.


        :return:
            Tuple containing the minimized potential energies and the XYZ format strings of the relaxed coordinates.

        """

        assert_msg_critical('openmm' in sys.modules, 'OpenMM is required for IMDatabasePointCollecter.')

        if self.system is None:
            raise RuntimeError('System has not been created!')
        if self.molecule is None:
            raise RuntimeError('Molecule object does not exist!')

        self.ensemble = ensemble
        self.temperature = temperature * unit.kelvin

        self.timestep = timestep * unit.femtosecond
        self.nsteps = nsteps

        self.integrator = self._create_integrator()
        topology = self.modeller.topology if self.phase in ['water', 'periodic'] else self.pdb.topology
        self.positions = self.modeller.positions if self.phase in ['water', 'periodic'] else self.pdb.positions

        self.simulation = app.Simulation(topology, self.system, self.integrator,  platform=self._create_platform())
        self.simulation.context.setPositions(self.positions)
        self.simulation.context.setVelocitiesToTemperature(self.temperature)

        save_freq = nsteps // snapshots if snapshots else nsteps
        energies = []
        opt_coordinates = []

        for step in range(nsteps):
            self.simulation.step(1)
            if step % save_freq == 0:
                state = self.simulation.context.getState(getPositions=True, getEnergy=True)
                energy = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
                coordinates = state.getPositions()
                self.simulation.context.setPositions(coordinates) 

                if minimize:

                    print(f'Step: {step}, Potential energy: {energy}')
                    self.simulation.minimizeEnergy()
                    minimized_state = self.simulation.context.getState(getPositions=True, getEnergy=True)
                    minimized_energy = minimized_state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
                    minimized_coordinates = minimized_state.getPositions()

                    energies.append(minimized_energy)
                    print(f'Minimized energy: {minimized_energy}')
                    xyz = f"{len(self.labels)}\n\n"
                    for label, coord in zip(self.labels, minimized_coordinates):
                        xyz += f"{label} {coord.x * 10} {coord.y * 10} {coord.z * 10}\n"  # Convert nm to Angstroms
                    print('Saved coordinates for step', step)
                    opt_coordinates.append(xyz)

                # Create a xyz with the coordinates
                print('Energy of the conformation:', energy)
                xyz = f"{len(self.labels)}\n\n"
                for label, coord in zip(self.labels, coordinates):
                    xyz += f"{label} {coord.x * 10} {coord.y * 10} {coord.z * 10}\n"
                opt_coordinates.append(xyz)
                energies.append(energy)
                print('Saved coordinates for step', step)

        print('Conformational sampling completed!')
        print(f'Number of conformations: {len(opt_coordinates)}')

        return energies, opt_coordinates

    def sort_points_with_association(self, points, associated_list):
        '''
        This function is sorting the points from a given database
        starting from the first entry to the last (point_0,..., point_n)

        :param points:
            list of database points
        :param associated_list:
            list of labels that are associated to the individual datapoints

        '''
        
        # Define a custom key function to extract the numerical part
        def extract_number(point):
            return int(point.split('_')[1])

        # Combine the points and associated list into tuples
        combined = list(zip(points, associated_list))

        # Sort the combined list using the custom key function
        sorted_combined = sorted(combined, key=lambda x: extract_number(x[0]))

        # Unzip the sorted combined list back into two lists
        sorted_points, sorted_associated_list = zip(*sorted_combined)

        return list(sorted_points), list(sorted_associated_list)
    
    def _validate_metadynamics_settings(self):

        cfg = self.metadynamics_settings

        if Metadynamics is None or BiasVariable is None:

            raise RuntimeError('metadynamics requested, but it is currently unavaliable!')

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
                for a in atoms:
                    if a < 0 or a >= natoms:
                        raise ValueError(f"metadynamics.variables[{idx}] atom index {a} out of range [0, {natoms-1}].")

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
        Updates settings in the ImpesDynamicsDriver.
        :param molecule:
            The Molecule object.
        :param impes_dict:
            The input dictionary of settings for IMPES.
        """
        
        self.drivers = dynamics_settings['drivers']

        self.interpolation_settings = interpolation_settings
        self.sampling_interpolation_settings = sampling_interpolation_settings
        self.dynamics_settings_interpolation_run = dynamics_settings

        if 'print_step' in dynamics_settings:
            self.print_step = int(dynamics_settings['print_step'])

        if 'duration' in dynamics_settings:
            self.duration = float(dynamics_settings['duration'])

        if 'temperature' in dynamics_settings:
            self.temperature = dynamics_settings['temperature']
            self.starting_temperature = dynamics_settings['temperature']

        if 'pressure' in dynamics_settings:
            self.pressure = float(dynamics_settings['pressure'])

        if 'force_constant' in dynamics_settings:
            self.force_constant = float(dynamics_settings['force_constant'])

        if 'friction' in dynamics_settings:
            self.friction = float(dynamics_settings['friction'])

        if 'timestep' in dynamics_settings:
            self.timestep = float(dynamics_settings['timestep'])

        if 'nsteps' in dynamics_settings:
            self.nsteps = int(dynamics_settings['nsteps'])

        if 'snapshots' in dynamics_settings:
            self.snapshots = int(dynamics_settings['snapshots'])

        # Start the simulation form a given state
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
        #################################### DATABASE construciton inputs #############################

        # The desired density around a given starting structure/datapoint
        if 'desired_datapoint_density' in dynamics_settings:
            self.desired_datpoint_density = int(dynamics_settings['desired_datapoint_density'])

        # The number of iteration cycles without database expansion after which the description around the reference point has converged
        if 'converged_cycle' in dynamics_settings:
            self.unadded_cycles = int(dynamics_settings['converged_cycle'])
        
        if 'basis_set_label' in dynamics_settings:
            basis_set_label = dynamics_settings['basis_set_label']
            self.basis_set_label = basis_set_label
    
        # Dertermines if non-adiabatic couplings should be calculate

        if 'cluster_run' in dynamics_settings:
            self.cluster_run = dynamics_settings['cluster_run']
        # TODO: these need to be different dictionary entries
        if 'start_collect' in dynamics_settings:
            self.start_collect = dynamics_settings['start_collect']

        # time step at which QM data point collection should stop
        if 'qmc_stop' in dynamics_settings:
            self.qmc_stop = float(dynamics_settings["qmc_stop"])

        # index of the excited state (in case of TDDFT QM data points)
        if "roots_to_follow" in dynamics_settings:
            self.nstates = len(list(dynamics_settings["roots_to_follow"]))
            self.roots_to_follow = list(dynamics_settings['roots_to_follow'])

        if 'excitation_pulse' in dynamics_settings and dynamics_settings['excitation_pulse'] is not None:
            self.excitation_pulse = list(dynamics_settings['excitation_pulse'])
        
        if 'energy_threshold' in dynamics_settings:
            self.energy_threshold = dynamics_settings['energy_threshold']
        
        if 'grad_rmsd_thrsh' in dynamics_settings:
            self.gradient_rmsd_thrsh = dynamics_settings['grad_rmsd_thrsh']

        if 'sampling_drivers' in dynamics_settings:
            sampling_settings = dynamics_settings.get('sampling_settings', {'enabled':False,
            'e_thrsh_kcal_per_atom': 0.1,
            'g_thrsh_kcal_ang_per_atom':2.0,
            'force_orient_cos': 0.0001})
            self.sampling_enabled = bool(sampling_settings.get('enabled'))
            self.sampling_driver = dynamics_settings.get('sampling_drivers', None)
            self.sampling_settings = sampling_settings

        if 'profile_runtime_timing' in dynamics_settings:
            self.profile_runtime_timing = bool(dynamics_settings['profile_runtime_timing'])
        if 'profile_runtime_print_interval' in dynamics_settings:
            self.profile_runtime_print_interval = int(dynamics_settings['profile_runtime_print_interval'])
        if 'profile_interpolation_timing' in dynamics_settings:
            self.profile_interpolation_timing = bool(dynamics_settings['profile_interpolation_timing'])
        if 'profile_interpolation_print_summary' in dynamics_settings:
            self.profile_interpolation_print_summary = bool(
                dynamics_settings['profile_interpolation_print_summary'])
        if 'mpi_control_plane_enabled' in dynamics_settings:
            self.mpi_control_plane_enabled = bool(dynamics_settings['mpi_control_plane_enabled'])
        if 'mpi_root_worker_mode' in dynamics_settings:
            self.mpi_root_worker_mode = bool(dynamics_settings['mpi_root_worker_mode'])
        if 'mpi_reload_from_hdf5' in dynamics_settings:
            self.mpi_reload_from_hdf5 = bool(dynamics_settings['mpi_reload_from_hdf5'])
        if 'mpi_debug_sync' in dynamics_settings:
            self.mpi_debug_sync = bool(dynamics_settings['mpi_debug_sync'])

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
            'dof': self.system.getNumParticles() * 3 - self.system.getNumConstraints(),
            'has_temp_method': hasattr(self.integrator, 'computeSystemTemperature'),
            'gas_constant': gas_constant,
            'metadynamics_enabled': bool(self._metadynamics is not None),
            'metadynamics_group':{int(self._metadynamics_forcegroup)} if self._metadynamics is not None else None,
        }

    def _initialize_runtime_profilers(self):
        """
        Initializes run-time profilers for QM/MM loop and interpolation drivers.
        """

        if not self.profile_runtime_timing:
            self._qmmm_runtime_profiler = None
            self._qmmm_runtime_totals = {}
            self._qmmm_runtime_step = {}
            self._qmmm_runtime_steps = 0
            return

        self._qmmm_runtime_profiler = Profiler({'timing': True})
        self._qmmm_runtime_profiler.set_timing_key('run_qmmm')
        self._qmmm_runtime_totals = {}
        self._qmmm_runtime_step = {}
        self._qmmm_runtime_steps = 0

    def _add_runtime_timing(self, label, dt):
        """
        Adds timing value to QM/MM profiling aggregates.
        """

        if not self.profile_runtime_timing or self._qmmm_runtime_profiler is None:
            return

        self._qmmm_runtime_profiler.add_timing_info(label, dt)
        self._qmmm_runtime_totals[label] = self._qmmm_runtime_totals.get(label, 0.0) + dt
        self._qmmm_runtime_step[label] = self._qmmm_runtime_step.get(label, 0.0) + dt

    def _finalize_runtime_step(self):
        """
        Finalizes one profiled timestep.
        """

        if not self.profile_runtime_timing or self._qmmm_runtime_profiler is None:
            return

        self._qmmm_runtime_steps += 1
        self._qmmm_runtime_step = {}

    def _print_qmmm_runtime_profile_summary(self):
        """
        Prints runtime breakdown for one average QM/MM timestep.
        """

        if not self.profile_runtime_timing or self._qmmm_runtime_steps == 0:
            return

        steps = self._qmmm_runtime_steps
        totals = self._qmmm_runtime_totals
        avg_step = totals.get('step.total', 0.0) / steps
        avg_update_forces = totals.get('step.update_forces', 0.0) / steps
        avg_impes = totals.get('update_gradient_and_energy.impes_compute_total', 0.0) / steps
        avg_integrator = totals.get('step.integrator', 0.0) / steps
        avg_observables = totals.get('step.observables', 0.0) / steps

        def pct(value):
            return (100.0 * value / avg_step) if avg_step > 0.0 else 0.0

        print('QM/MM runtime profiler summary (average per timestep)')
        print(f'  steps_profiled: {steps}')
        print(f'  timestep.total: {avg_step:.6f} s')
        print(f'  timestep.update_forces: {avg_update_forces:.6f} s ({pct(avg_update_forces):.1f}%)')
        print(f'  timestep.integrator_step: {avg_integrator:.6f} s ({pct(avg_integrator):.1f}%)')
        print(f'  timestep.observables: {avg_observables:.6f} s ({pct(avg_observables):.1f}%)')
        print(f'  impes_driver.compute_total: {avg_impes:.6f} s ({pct(avg_impes):.1f}%)')

        root_keys = sorted(
            key for key in totals if key.startswith('update_gradient_and_energy.impes_compute_root_')
        )
        for key in root_keys:
            avg_val = totals[key] / steps
            print(f'  {key.split(".")[-1]}: {avg_val:.6f} s ({pct(avg_val):.1f}%)')

    def _print_interpolation_runtime_profile_summary(self):
        """
        Prints bottleneck summary from InterpolationDriver runtime profiler.
        """

        if not self.profile_interpolation_timing or not self.profile_interpolation_print_summary:
            return
        if not isinstance(self.impes_drivers, dict):
            return

        for root in self.roots_to_follow:
            driver = self.impes_drivers.get(root)
            if driver is None:
                continue
            summary = driver.get_runtime_profile_summary()
            if not summary['enabled'] or summary['calls'] == 0:
                continue

            totals = summary['totals']
            calls = summary['calls']
            compute_total = totals.get('compute.total', 0.0)
            avg_compute = compute_total / calls if calls > 0 else 0.0

            print(f'InterpolationDriver profiler summary (root={root})')
            print(f'  compute.calls: {calls}')
            print(f'  compute.total: {compute_total:.6f} s (avg {avg_compute:.6f} s/call)')

            hot_labels = [
                'compute.prepare_state',
                'compute.define_impes_coordinate',
                'compute.define_impes_coordinate.define_internal_coordinates',
                'compute.define_impes_coordinate.calculate_b_matrix',
                'compute.define_impes_coordinate.compute_internal_coordinates_values',
                'compute.load_qm_data_points',
                'compute.prepare_runtime_cache',
                'compute.shepard_interpolation',
                'shepard.total',
                'shepard.distance_scan',
                'shepard.point_eval_loop',
                'shepard.compute_potential',
                'shepard.assembly',
                'shepard.filter_close_points',
            ]
            for label in hot_labels:
                if label in totals:
                    avg_val = totals[label] / calls
                    pct_compute = (100.0 * totals[label] / compute_total) if compute_total > 0.0 else 0.0
                    print(f'  {label}: total {totals[label]:.6f} s, avg {avg_val:.6f} s/call ({pct_compute:.1f}% of compute.total)')

    def _collect_step_observables_reference(self, context, runtime_cache):
        """
        Reference observable collection using Quantity arithmetic.
        This mirrors the original implementation for verification.
        """

        energy_unit = runtime_cache['energy_unit']
        dof = runtime_cache['dof']
        
        qm = self.get_qm_potential_energy()
        qm_mm = context.getState(
            getEnergy=True,
            groups=runtime_cache['qm_mm_groups']).getPotentialEnergy()
        mm = context.getState(
            getEnergy=True,
            groups=runtime_cache['mm_groups']).getPotentialEnergy()
        
        meta_bias = 0.0

        if runtime_cache['metadynamics_enabled']:
            group = runtime_cache['metadynamics_group']
            meta_bias_q = context.getState(getEnergy=True, groups=group).getPotentialEnergy()
            meta_bias = meta_bias_q.value_in_unit(energy_unit)


        pot = qm * energy_unit + qm_mm + mm + (meta_bias * energy_unit)
        kinetic = context.getState(getEnergy=True).getKineticEnergy()

        if runtime_cache['has_temp_method']:
            temp = self.integrator.computeSystemTemperature().value_in_unit(unit.kelvin)
        else:
            temp = (2 * kinetic / (dof * unit.MOLAR_GAS_CONSTANT_R)).value_in_unit(unit.kelvin)

        total = pot + kinetic

        return {
            'qm': float(qm),
            'qm_mm': qm_mm.value_in_unit(energy_unit),
            'mm': mm.value_in_unit(energy_unit),
            'metad_bias':float(meta_bias),
            'potential': pot.value_in_unit(energy_unit),
            'kinetic': kinetic.value_in_unit(energy_unit),
            'temperature': float(temp),
            'total': total.value_in_unit(energy_unit),
        }

    def _collect_step_observables_optimized(self, context, runtime_cache):
        """
        Optimized observable collection avoiding repeated Quantity arithmetic.
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

        if runtime_cache['has_temp_method']:
            temperature = self.integrator.computeSystemTemperature().value_in_unit(unit.kelvin)
        else:
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
        self.metad_bias_energies.append(observables.get('metad_bias', 0.0))
        self.total_potentials.append(observables['potential'])
        self.kinetic_energies.append(observables['kinetic'])
        self.temperatures.append(observables['temperature'])
        self.total_energies.append(observables['total'])

        self.gloabal_sim_informations['temperatures'].append(observables['temperature'])
        self.gloabal_sim_informations['state'].append(self.current_state)
        if self.coordinates_xyz:
            self.gloabal_sim_informations['coordinates_ang'].append(self.coordinates_xyz[-1])
        else:
            self.gloabal_sim_informations['coordinates_ang'].append(None)

    def _update_qmmm_verification_stats(self, reference, optimized, stats, step, atol, rtol, fail_fast):
        """
        Compares optimized and reference per-step observables.
        """

        keys = ('qm', 'qm_mm', 'mm', 'potential', 'kinetic', 'temperature', 'total')
        mismatch = []

        for key in keys:
            ref = float(reference[key])
            opt = float(optimized[key])
            abs_diff = abs(ref - opt)
            rel_diff = abs_diff / max(abs(ref), 1.0e-15)

            stats['max_abs'][key] = max(stats['max_abs'][key], abs_diff)
            stats['max_rel'][key] = max(stats['max_rel'][key], rel_diff)

            if abs_diff > (atol + rtol * abs(ref)):
                mismatch.append((key, ref, opt, abs_diff, rel_diff))

        if mismatch:
            stats['num_mismatch_steps'] += 1
            if fail_fast:
                mismatch_txt = ', '.join(
                    f"{key}: ref={ref:.12e}, opt={opt:.12e}, abs={abs_diff:.3e}, rel={rel_diff:.3e}"
                    for key, ref, opt, abs_diff, rel_diff in mismatch
                )
                raise RuntimeError(
                    f'run_qmmm verification failed at step {step}: {mismatch_txt}'
                )

            if stats['num_mismatch_steps'] <= 5:
                print(f'Verification warning at step {step}: {len(mismatch)} field(s) differ above tolerance')

    def _print_run_qmmm_step(self, step, timestep, observables, save_freq):
        """
        Prints per-step simulation diagnostics.
        """

        if step % save_freq != 0:
            return

        print(f"Step: {step} / {self.nsteps} Time: {round((step * timestep) / 1000, 2)} ps")
        print('Potential Energy QM region:', observables['qm'], 'kJ/mol')
        print('Potential Energy MM region:', observables['mm'])
        print('QM/MM Interaction Energy:', observables['qm_mm'])
        print('Total Potential Energy:', observables['potential'])
        print('Kinetic Energy:', observables['kinetic'])
        print('Temperature:', observables['temperature'], 'K')
        print('Total Energy:', observables['total'], '±', np.std(self.total_energies), 'kJ/mol')
        if 'metad_bias' in observables:
            print('Metadynamics Bias Energy:', observables['metad_bias'], 'kJ/mol')
        for root in self.roots_to_follow:
            print('Current Density', self.density_around_data_point[0][root], '-->', self.desired_datpoint_density, self.unadded_cycles)
        print('Current State (PES):', self.current_state)
        print('-' * 60)

    def _print_run_qmmm_verification_summary(self, stats, atol, rtol):
        """
        Prints max observed deviations between reference and optimized observables.
        """

        print('run_qmmm verification summary')
        print(f'Tolerances: atol={atol:.3e}, rtol={rtol:.3e}')
        print(f'Steps with mismatch above tolerance: {stats["num_mismatch_steps"]}')
        for key in ('qm', 'qm_mm', 'mm', 'potential', 'kinetic', 'temperature', 'total'):
            print(
                f'  {key}: max_abs={stats["max_abs"][key]:.3e}, '
                f'max_rel={stats["max_rel"][key]:.3e}'
            )

    def _initialize_sampling_impes_drivers(self, inv_sqrt_masses):
        self.sampling_impes_drivers = {root: None for root in self.roots_to_follow}
        self.sampling_qm_data_point_dict = {root: [] for root in self.roots_to_follow}
        self.sampling_qm_symmetry_datapoint_dict = {root: {} for root in self.roots_to_follow}

        for root in self.roots_to_follow:
            drv = InterpolationDriver(self.root_z_matrix[root])
            drv.update_settings(self.sampling_interpolation_settings[root])
            drv.symmetry_information = self.impes_drivers[root].symmetry_information
            drv.use_symmetry = self.use_symmetry
            drv.impes_coordinate.inv_sqrt_masses = inv_sqrt_masses

            labels, _ = drv.read_labels()
            old_label = None
            for label in labels:
                dp = InterpolationDatapoint(self.root_z_matrix[root])
                dp.update_settings(self.sampling_interpolation_settings[root])
                dp.read_hdf5(self.sampling_interpolation_settings[root]['imforcefield_file'], label)
                dp.inv_sqrt_masses = inv_sqrt_masses
                if dp.bank_role == "core" or "cluster" not in dp.point_label and "symmetry" not in dp.point_label:
                    self.sampling_qm_data_point_dict[root].append(dp)
                    old_label = dp.point_label
                    self.sampling_qm_symmetry_datapoint_dict[root][old_label] = [dp]
                elif dp.bank_role == "symmetry":
                    self.sampling_qm_symmetry_datapoint_dict[root][old_label].append(dp)
            
            drv.qm_data_points = self.sampling_qm_data_point_dict[root]
            drv.qm_symmetry_data_points = self.sampling_qm_symmetry_datapoint_dict[root]
            drv.prepare_runtime_data_cache(force=True)
            
            if len(self.sampling_qm_data_point_dict[root]) > 0:
                drv.impes_coordinate.eq_bond_lengths = self.sampling_qm_data_point_dict[root][0].eq_bond_lengths
            drv.qm_rotor_cluster_banks = self.qm_rotor_cluster_banks[root]
            
            if drv.qm_rotor_cluster_banks:
                first_family = next(iter(drv.qm_rotor_cluster_banks.values()))
                drv.rotor_cluster_information = first_family.get("cluster_info")
            else:
                drv.rotor_cluster_information = None
            
            self.sampling_impes_drivers[root] = drv
    
    def _reload_interpolation_root_from_hdf5(self, root, inv_sqrt_masses):
        driver_object = self.impes_drivers[root]

        im_labels, _ = driver_object.read_labels()

        self.qm_data_point_dict[root] = []
        self.qm_symmetry_datapoint_dict[root] = {}
        self.qm_energies_dict[root] = []
        self.sorted_state_spec_im_labels[root] = []

        old_label = None

        for label in im_labels:
            qm_data_point = InterpolationDatapoint(self.root_z_matrix[root])
            qm_data_point.update_settings(self.interpolation_settings[root])
            qm_data_point.read_hdf5(self.interpolation_settings[root]['imforcefield_file'], label)
            qm_data_point.inv_sqrt_masses = inv_sqrt_masses
            
            if qm_data_point.bank_role == "core" or "cluster" not in qm_data_point.point_label and "symmetry" not in qm_data_point.point_label:


                self.qm_data_point_dict[root].append(qm_data_point)
                self.qm_energies_dict[root].append(qm_data_point.energy)
                self.sorted_state_spec_im_labels[root].append(label)

                old_label = qm_data_point.point_label
                self.qm_symmetry_datapoint_dict[root][old_label] = [qm_data_point]
            elif qm_data_point.bank_role == "symmetry":

                self.qm_symmetry_datapoint_dict[root][old_label].append(qm_data_point)

        driver_object.qm_symmetry_data_points = self.qm_symmetry_datapoint_dict[root]
        driver_object.qm_data_points = self.qm_data_point_dict[root]
        driver_object.labels = self.sorted_state_spec_im_labels[root]
        
        print(len(self.qm_data_point_dict[root]), driver_object.impes_coordinate.eq_bond_lengths, self.qm_data_point_dict[root][0].eq_bond_lengths)
        if len(self.qm_data_point_dict[root]) > 0:
            driver_object.impes_coordinate.eq_bond_lengths = self.qm_data_point_dict[root][0].eq_bond_lengths

        driver_object.mark_runtime_data_cache_dirty()
        driver_object.prepare_runtime_data_cache(force=True)

        use_mpi_preload = bool(self.interpolation_settings[root].get('use_mpi_preload', False))
        # In root-worker mode, avoid immediate preload collectives here; rebuild lazily
        # in synchronized compute path on the next step.
        force_rebuild_preload = not (self._mpi_is_active() and self.mpi_root_worker_mode)
        driver_object.set_mpi_preload_engine(
            comm=self._mpi_interp_comm,
            enabled=(use_mpi_preload and self.nodes > 1),
            force_rebuild=force_rebuild_preload,
        )

        driver_object._mpi_preload_enabled_config = bool(getattr(driver_object, 'mpi_preload_enabled', False))
        self._ensure_interp_engine_comm(driver_object)
        # In root-worker mode, keep preload disabled outside synchronized compute
        # phases to prevent root-only branches from emitting collective packets.
        if self._mpi_is_active() and self.mpi_root_worker_mode:
            driver_object.mpi_preload_enabled = False


    def run_qmmm(self,
                verify_optimized=False,
                verification_atol=1.0e-8,
                verification_rtol=1.0e-6,
                verification_fail_fast=False,
                use_optimized_observables=True,
                test_hooks=None,
                collect_step_trace=False,
                strict_test_hooks=False,
                ):
        """
        Runs a QM/MM simulation using OpenMM, storing the trajectory and simulation data.

        :param verify_optimized:
            If True, computes both reference and optimized observables per step and compares them.
        :param verification_atol:
            Absolute tolerance used for reference-vs-optimized checks.
        :param verification_rtol:
            Relative tolerance used for reference-vs-optimized checks.
        :param verification_fail_fast:
            If True, raises RuntimeError on the first mismatch above tolerance.
        :param use_optimized_observables:
            If True, stores observables from the optimized collector.
            If False, stores observables from the legacy/reference collector.
        """

        assert_msg_critical('openmm' in sys.modules, 'OpenMM is required for IMDatabasePointCollecter.')
        
        if self.system is None:
            raise RuntimeError('System has not been created!')
        
        # test_hook section for consistent set up tests
        self._set_test_hooks(hooks=test_hooks, 
                             collect_step_trace=collect_step_trace,
                             strict=strict_test_hooks)
        
        self._emit_test_hook("run_start", {
            "nsteps": int(self.nsteps),
            "ensemble": str(self.ensemble),
            "roots_to_follow": list(self.roots_to_follow),
        })

        # temperature = self.temperature
        self.temperature_number = self.temperature
        self.temperature = self.temperature * unit.kelvin
        
        # friction = self.friction
        self.friction = self.friction / unit.picosecond
        
        timestep = self.timestep
        self.timestep = timestep * unit.femtoseconds

        self.qm_potentials = []
        self.qm_mm_interaction_energies = []
        self.mm_potentials = []
        self.total_potentials = []
        self.kinetic_energies = []
        self.temperatures = []
        self.total_energies = []
        self.coordinates_xyz = []

        save_freq = max(1, self.nsteps // self.snapshots) if self.snapshots else max(1, self.nsteps)

        is_worker_service_rank = (
            self._mpi_is_active() and self.mpi_root_worker_mode and (not self._mpi_is_root())
        )

        self.topology = self.pdb.topology

        self.positions = self.pdb.positions
        self._topology_atom_labels = [atom.element.symbol for atom in self.topology.atoms()]

        if not is_worker_service_rank:
            # Create or update the integrator
            if self.integrator is None:
                new_integrator = self._create_integrator()
            else:
                new_integrator = self.integrator

            # initializing the metadynamics variables

            self.metad_bias_energies = []
            self.metad_cv_history = []
            self._metadynamics = None

            self.metadynamics_bias_energies = []
            self.metadynamics_cv_history = []

            if self.metadynamics_enabled:
                # Optional guard to avoid double-biasing in early rollout
                if self.bias_force_reaction_prop is not None:
                    raise RuntimeError(
                        "Both legacy bias_force_reaction_prop and metadynamics are enabled. "
                        "Disable one to avoid conflicting bias potentials."
                    )
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
                # else:
                #     self.simulation.context.setVelocitiesToTemperature(150 * unit.kelvin)
            self.start_velocities = self.simulation.context.getState(getVelocities=True).getVelocities()

            # Set up reporting
            self.simulation.reporters.clear()
            self.simulation.reporters.append(app.PDBReporter(self.out_file, save_freq))

            # Print header
            print('QM/MM Simulation Parameters')
            print('=' * 60)
            print('QM Driver:', self.driver_flag)
            print('Ensemble:', self.ensemble)
            print('Integration method:', new_integrator.__class__.__name__)
            if self.ensemble in ['NVT', 'NPT']:
                print('Temperature:', self.temperature, 'K')
                if self.ensemble == 'NPT':
                    print('Pressure:', self.pressure, 'atm')
            print('Friction:', self.friction, '1/ps')
            print('Timestep:', self.timestep, 'fs')
            print('Total simulation time in ns:', self.nsteps * timestep / 1e6)
            print('=' * 60)
        else:
            new_integrator = None
            self.simulation = None
            self.start_velocities = None

        self.last_point_added = 0
        self.cycle_iteration = self.unadded_cycles

        self.allowed_molecules = {root: {'molecules': [], 'qm_energies': [], 'im_energies':[], 'qm_gradients':[], 'distances': []} for root in self.roots_to_follow}

        if not is_worker_service_rank:
            print('current datapoints around the given starting structure', self.desired_datpoint_density, self.density_around_data_point[0], self.density_around_data_point[1], '\n allowed derivation from the given structure', self.allowed_molecule_deviation, '\n ---------------------------------------')
            openmm_coordinate = self.simulation.context.getState(getPositions=True).getPositions()
            self.coordinates = [openmm_coordinate]
        else:
            self.coordinates = []

        self.coordinates_xyz = []
        self.gradients = []
        self.velocities = []
        self.velocities_np = []
        if not is_worker_service_rank:
            self.velocities_np.append(self.simulation.context.getState(getVelocities=True).getVelocities(True))

        self.impes_drivers = []

        self.gloabal_sim_informations = {f'state_{root}':{'pot_energies':[], 'gradients':[]} for root in self.roots_to_follow}
        self.gloabal_sim_informations['coordinates_ang'] = []
        self.gloabal_sim_informations['state'] = []
        self.gloabal_sim_informations['temperatures'] = []

        self.state_specific_molecules = {root: [] for root in self.roots_to_follow}
        self.sorted_state_spec_im_labels = {root: [] for root in self.roots_to_follow}

        self.transition_information = []
        self.prev_dens_of_points = {root: 0 for root in self.roots_to_follow}
        self.bayes_models = {root: [] for root in self.roots_to_follow}
        self.qm_data_point_dict = {root: [] for root in self.roots_to_follow}
        self.qm_symmetry_datapoint_dict = {root: {} for root in self.roots_to_follow}
        self.qm_energies_dict = {root: [] for root in self.roots_to_follow}
        self.impes_drivers = {root: None for root in self.roots_to_follow}
        self.sampled_molecules = {root: {'molecules': [], 'im_energies': [], 'qm_energies': [], 'qm_gradients': [], 'distances': []} for root in self.roots_to_follow}
        self.qm_rotor_cluster_banks = {root: {} for root in self.roots_to_follow}

        if self.reference_struc_energies_file is not None:
            print('Extracting reference structures from', self.reference_struc_energies_file)
            self.sampled_molecules = self.extract_reference_structures_h5(self.reference_struc_energies_file, self.roots_to_follow)
            self.allowed_molecules = self.extract_reference_structures_h5(self.reference_struc_energies_file, self.roots_to_follow)
            self.last_added += len(self.sampled_molecules[self.roots_to_follow[0]]['molecules'])

        else:
             self.reference_struc_energies_file = 'QM_ref_along_traj.h5'

        masses = self.molecule.get_masses().copy()
        masses_cart = np.repeat(masses, 3)
        inv_sqrt_masses = 1.0 / np.sqrt(masses_cart)
        self.inv_sqrt_masses = inv_sqrt_masses
        for root in self.roots_to_follow:
            # Dynamically create an attribute name
            attribute_name = f'impes_driver_{root}'
            # Initialize the object
            driver_object = InterpolationDriver(self.root_z_matrix[root])
            driver_object.update_settings(self.interpolation_settings[root])
            driver_object.enable_runtime_profiling(
                enabled=self.profile_interpolation_timing,
                reset=True,
                print_summary=False)
            # print('Interpolation driver settings updated for root', root, self.eq_bond_force_constants)
            
            driver_object.impes_coordinate.inv_sqrt_masses = inv_sqrt_masses
            if root == 0:
                driver_object.symmetry_information = self.non_core_symmetry_groups['gs']
            elif root == 1 and self.drivers['es'] is not None and self.drivers['es'][0].spin_flip:
                driver_object.symmetry_information = self.non_core_symmetry_groups['gs']
            else:
                driver_object.symmetry_information = self.non_core_symmetry_groups['es']
            

            driver_object.use_symmetry = self.use_symmetry
            self.impes_drivers[root] = driver_object
            self._reload_interpolation_root_from_hdf5(root, inv_sqrt_masses)
            self.prev_dens_of_points[root] = len(self.qm_data_point_dict[root])
            self.qm_rotor_cluster_banks[root] = self._load_rotor_cluster_bank_for_root(root)
            driver_object.qm_rotor_cluster_banks = self.qm_rotor_cluster_banks[root]

            if driver_object.qm_rotor_cluster_banks:
                first_family = next(iter(driver_object.qm_rotor_cluster_banks.values()))
                driver_object.rotor_cluster_information = first_family.get("cluster_info")
            else:
                driver_object.rotor_cluster_information = None
            print('In set up', driver_object.rotor_cluster_information)
            # Set the object as an attribute of the instance
            setattr(self, attribute_name, driver_object)
            # Append the object to the list
            self.impes_drivers[root] = driver_object

        if self.sampling_enabled:
            self._initialize_sampling_impes_drivers(inv_sqrt_masses)

        self.current_state = self.roots_to_follow[0]
 
        self.swap_back = True
        self.prev_state = self.current_state
        self.step = 0 
        self.prev_dE_gpr = 0.0
        self.last_gpr_addition = 0

        if self._mpi_is_active() and self.mpi_root_worker_mode:
            self._mpi_pending_sync_roots.clear()
            # Ensure all initialization-time MPI collectives have completed
            # before starting the control-channel STEP/STOP loop.
            self._mpi_ctrl_comm.Barrier()
            if self._mpi_interp_comm is not self._mpi_ctrl_comm:
                self._mpi_interp_comm.Barrier()
            self._mpi_worker_service_running = True
            if not self._mpi_is_root():
                self._run_qmmm_worker_service()
                self._mpi_worker_service_running = False
                return None

        start_time = time()
        runtime_cache = self._build_run_qmmm_runtime_cache()
        self._initialize_runtime_profilers()
        verification_stats = None
        if verify_optimized:
            verification_stats = {
                'num_mismatch_steps': 0,
                'max_abs': {
                    'qm': 0.0,
                    'qm_mm': 0.0,
                    'mm': 0.0,
                    'potential': 0.0,
                    'kinetic': 0.0,
                    'temperature': 0.0,
                    'total': 0.0,
                },
                'max_rel': {
                    'qm': 0.0,
                    'qm_mm': 0.0,
                    'mm': 0.0,
                    'potential': 0.0,
                    'kinetic': 0.0,
                    'temperature': 0.0,
                    'total': 0.0,
                },
            }

        self.last_force_point = 1
        
        for step in range(self.nsteps):

            step_t0 = time()

            force_t0 = time()

            if self._mpi_is_active() and self.mpi_root_worker_mode and self._mpi_is_root():
                qm_positions_nm = np.array([
                    p.value_in_unit(unit.nanometer)
                    for p in self.simulation.context.getState(getPositions=True).getPositions()
                ])[self.qm_atoms]

                sync_roots = sorted(self._mpi_pending_sync_roots)
                self._mpi_pending_sync_roots.clear()

                if self.mpi_reload_from_hdf5 and len(sync_roots) > 0:
                    for root_sync in sync_roots:
                        self._reload_interpolation_root_from_hdf5(root_sync, self.inv_sqrt_masses)

                self._mpi_bcast_control({
                    'cmd': self._MPI_CMD_STEP,
                    'step_idx': int(step),
                    'qm_positions_nm': np.asarray(qm_positions_nm, dtype=np.float64),
                    'sync_roots': sync_roots,
                })


            self.update_forces(self.simulation.context)
            self._add_runtime_timing('step.update_forces', time() - force_t0)

            obs_t0 = time()
            reference_observables = None
            optimized_observables = None

            if verify_optimized or not use_optimized_observables:
                reference_observables = self._collect_step_observables_reference(
                    self.simulation.context,
                    runtime_cache)
            if verify_optimized or use_optimized_observables:
                optimized_observables = self._collect_step_observables_optimized(
                    self.simulation.context,
                    runtime_cache)

            if use_optimized_observables:
                observables = optimized_observables
            else:
                observables = reference_observables

            if verify_optimized:
                self._update_qmmm_verification_stats(
                    reference_observables,
                    optimized_observables,
                    verification_stats,
                    step,
                    verification_atol,
                    verification_rtol,
                    verification_fail_fast)

            self._append_step_observables(observables)
            self._print_run_qmmm_step(step, timestep, observables, save_freq)
            self._add_runtime_timing('step.observables', time() - obs_t0)

            if (
                self.bias_force_reaction_prop is not None
                and len(self.bias_force_reaction_prop) != 0
                and len(self.global_theta_list) > self.last_force_point
                and self.point_checker != 0
                and self.step % self.bias_force_reaction_prop[3] == 0
            ):
                self.simulation.context.setParameter('theta0', self.global_theta_list[self.last_force_point])
                self.last_force_point += 1

            self.step += 1

            step_payload = {
                "step": int(step),
                "current_state": int(self.current_state),
                "point_checker": int(self.point_checker),
                "add_a_point" : bool(self.add_a_point),
                "qm_energy_kjmol": float(observables["qm"]),
                "potential_kjmol": float(observables["potential"]),
                "kinetic_kjmol": float(observables["kinetic"]),
                "total_kjmol": float(observables["total"]),
                "temperature_K": float(observables["temperature"]),
                "n_datapoints_per_root": {
                    int(root): int(len(self.qm_data_point_dict[root]))
                    for root in self.roots_to_follow
                },
            }
            self._emit_test_hook("step", step_payload)

            # self._update_induction_embedding(phase == 'periodic')
            integ_t0 = time()

            if self._metadynamics is not None:
                self._metadynamics.step(self.simulation, 1)
            else:
                self.simulation.step(1)


            self._add_runtime_timing('step.integrator', time() - integ_t0)
            self.prev_state = self.current_state
            self._add_runtime_timing('step.total', time() - step_t0)
            self._finalize_runtime_step()
            if (self.profile_runtime_timing and self.profile_runtime_print_interval > 0
                    and (step + 1) % self.profile_runtime_print_interval == 0):
                self._print_qmmm_runtime_profile_summary()
            
            reached_target_density = False
            for root in self.roots_to_follow:
                if self.density_around_data_point[0][root] >= self.desired_datpoint_density and self.expansion:
                    reached_target_density = True
                    break
            if reached_target_density:
                if self.use_opt_confidence_radius[0] and self.confidence_radius_optimized:
                    break
                elif not self.use_opt_confidence_radius[0]:
                    break
            if not self.expansion:
                if len(self.expansion_molecules) >= 20:
                    break

            if self.unadded_cycles == 0:
                break

        if self._mpi_is_active() and self.mpi_root_worker_mode and self._mpi_is_root():
            self._mpi_bcast_control({'cmd': self._MPI_CMD_STOP})
        
        if self._mpi_is_active() and self.mpi_root_worker_mode:
            self._mpi_worker_service_running = False

        end_time = time()
        elapsed_time = end_time - start_time
        elapsed_time_days = elapsed_time / (24 * 3600)
        performance = (self.nsteps * timestep / 1e6) / elapsed_time_days
        self._emit_test_hook("run_end", {
            "steps_collected": int(len(self.total_energies)),
            "final_state": int(self.current_state),
            "n_datapoints_per_root": {
                int(root): int(len(self.qm_data_point_dict[root]))
                for root in self.roots_to_follow
            },
        })
        if verify_optimized:
            self._print_run_qmmm_verification_summary(
                verification_stats,
                verification_atol,
                verification_rtol)
        self._print_qmmm_runtime_profile_summary()
        self._print_interpolation_runtime_profile_summary()
        self.output_file_writer(self.summary_output)
        print('QM/MM simulation completed!')
        print(f'Number of steps: {self.nsteps}')
        print(f'Trajectory saved as {self.out_file}')
        print('=' * 60)
        print('Simulation Averages:')
        print('=' * 60)
        print('QM Potential Energy:', np.mean(self.qm_potentials), '±', np.std(self.qm_potentials), 'kJ/mol')
        print('QM/MM Interaction Energy:', np.mean(self.qm_mm_interaction_energies), '±', np.std(self.qm_mm_interaction_energies), 'kJ/mol')
        print('MM Potential Energy:', np.mean(self.mm_potentials), '±', np.std(self.mm_potentials), 'kJ/mol')
        print('Total Potential Energy:', np.mean(self.total_potentials), '±', np.std(self.total_potentials), 'kJ/mol')
        print('Kinetic Energy:', np.mean(self.kinetic_energies), '±', np.std(self.kinetic_energies), 'kJ/mol')
        print('Temperature:', np.mean(self.temperatures), '±', np.std(self.temperatures), 'K')
        print('Total Energy:', np.mean(self.total_energies), '±', np.std(self.total_energies), 'kJ/mol')
        print('=' * 60)
        print(f'Elapsed time: {int(elapsed_time // 60)} minutes, {int(elapsed_time % 60)} seconds')
        print(f'Performance: {performance:.2f} ns/day')
        print(f'Trajectory saved as {self.out_file}')

        if verify_optimized:
            return verification_stats
        return None


    # Post-simulation analysis methods
    def plot_energy(self,
                    components=['total'],
                    filename='energy_plot.png'):
        """
        Plots the energy components of the simulation.

        :param components:
            List of energy components to plot. Default is ['total'].
            It can include 'total', 'kinetic', 'potential', 'temperature', 'qm', 'mm', 'qm_mm'.
        :param filename:
            Name of the output file. Default is 'energy_plot.png'.
        """
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(10, 6))

        # Plot also the trend of the energy components as a linear regression
        if 'total' in components:
            ax.plot(self.total_energies, label='Total Energy', marker='o')
            # Linear regression excluding the first 10% of the data
            x = np.arange(len(self.total_energies))
            m, b = np.polyfit(x[int(0.1 * len(x)):], self.total_energies[int(0.1 * len(x)):], 1)
            ax.plot(x, m * x + b, label=f'Total Energy Trend ({m:.4e} kJ/mol/step)', linestyle='--')

        if 'kinetic' in components:
            ax.plot(self.kinetic_energies, label='Kinetic Energy', marker='o')

        if 'potential' in components:
            ax.plot(self.total_potentials, label='Potential Energy', marker='o')
        if 'temperature' in components:
            ax.plot(self.temperatures, label='Temperature', marker='o')
            ax.set_ylabel('Temperature (K)')
        if 'qm' in components:
            ax.plot(self.qm_potentials, label='QM Energy', marker='o')
        if 'mm' in components:
            ax.plot(self.mm_potentials, label='MM Energy', marker='o')
        if 'qm_mm' in components:
            ax.plot(self.qm_mm_interaction_energies, label='QM/MM Interaction Energy', marker='o')

        ax.set_title('Energy Components of the Simulation')

        ax.set_xlabel('Step')
        ax.set_ylabel('Energy (kJ/mol)')
        ax.legend()
        plt.tight_layout()
        plt.savefig(filename)
        plt.show()
        
    def visualize_trajectory(self, 
                             trajectory_file='trajectory.pdb', 
                             interval=1):
        """
        Visualizes the trajectory of the simulation using py3Dmol.

        :param trajectory_file:
            Path to the PDB file containing the trajectory. Default is 'trajectory.pdb'.
        """
        try:
            import py3Dmol

        except ImportError:
            raise ImportError("py3Dmol is not installed. Please install it using `pip install py3Dmol`.")
        
        viewer = py3Dmol.view(width=800, height=600)

        viewer.addModelsAsFrames(open(trajectory_file, 'r').read(),'pdb', {'keepH': True})
        viewer.animate({'interval': interval, 'loop': "forward", 'reps': 10})
        viewer.setStyle({"stick":{},"sphere": {"scale":0.25}})
        viewer.zoomTo()

        viewer.show()

    # Private methods
    def _save_output(self, output_file):
        """
        Save the simulation output to an out file.
        Reports all the energy components for every snapshot.

        :param output_file:
            Path to the output file.
        """

        with open(output_file + '.out', 'w', encoding='utf-8') as out:
            # Write the header
            out.write('VeloxChem/OpenMM Simulation Output\n')
            out.write('=' * 60 + '\n')
            if self.qm_atoms is not None:
                out.write(f'Simulation type: {self.driver_flag}\n')
            else:
                out.write('Simulation type: Classical MD\n')
            out.write(f'Ensemble: {self.ensemble}\n')
            out.write(f'Temperature: {self.temperature}\n')
            out.write(f'Friction: {self.friction}\n')
            out.write(f'Timestep: {self.timestep} \n')
            out.write(f'Number of steps: {self.nsteps}\n')
            out.write('=' * 60 + '\n')

            # Write the energy components from the lists
            for i in range(len(self.total_energies)):
                out.write(f'Step: {i}\n')
                out.write(f'Potential Energy (kJ/mol): {self.total_potentials[i]:.4f}\n')
                out.write(f'Kinetic Energy (kJ/mol): {self.kinetic_energies[i]:.4f}\n')
                out.write(f'Temperature (K): {self.temperatures[i]:.4f}\n')
                if self.qm_atoms is not None:
                    out.write(f'QM Potential Energy (kJ/mol): {self.qm_potentials[i]:.4f}\n')
                    out.write(f'QM/MM Interaction Energy (kJ/mol): {self.qm_mm_interaction_energies[i]:.4f}\n')
                    out.write(f'MM Potential Energy (kJ/mol): {self.mm_potentials[i]:.4f}\n')
                out.write(f'Total Energy (kJ/mol): {self.total_energies[i]:.4f}\n')
                out.write('-' * 60 + '\n')
            # Write the averages
            out.write('Summary\n')
            out.write('=' * 60 + '\n')
            out.write(f'Average Potential Energy (kJ/mol): {np.mean(self.total_potentials):.4f} ± {np.std(self.total_potentials):.4f}\n')
            out.write(f'Average Kinetic Energy (kJ/mol): {np.mean(self.kinetic_energies):.4f} ± {np.std(self.kinetic_energies):.4f}\n')
            out.write(f'Average Temperature (K): {np.mean(self.temperatures):.4f} ± {np.std(self.temperatures):.4f}\n')
            if self.qm_atoms is not None:
                out.write(f'Average QM Potential Energy (kJ/mol): {np.mean(self.qm_potentials):.4f} ± {np.std(self.qm_potentials):.4f}\n')
                out.write(f'Average QM/MM Interaction Energy (kJ/mol): {np.mean(self.qm_mm_interaction_energies):.4f} ± {np.std(self.qm_mm_interaction_energies):.4f}\n')
                out.write(f'Average MM Potential Energy (kJ/mol): {np.mean(self.mm_potentials):.4f} ± {np.std(self.mm_potentials):.4f}\n')
            out.write(f'Average Total Energy (kJ/mol): {np.mean(self.total_energies):.4f} ± {np.std(self.total_energies):.4f}\n')
            out.write('=' * 60 + '\n')      

    def _format_PDB_file(self, filename):
        """
        Check the format of the PDB file and correct it if necessary.

        :param pdb_file:
            Path to the PDB file to be checked.
        :return:
            Path to the corrected PDB file.
        """

        pdb_path = Path(filename)
        if not pdb_path.is_file():
            raise FileNotFoundError(f"{filename} does not exist.")
        
        pdbstr = pdb_path.read_text()

        # Formatting issues flags
        label_guess_warning = False
        conect_warning = False

        for line in pdbstr.strip().splitlines():
            if line.startswith(('ATOM', 'HETATM')):
                atom_label = line[76:78].strip()
                if not atom_label:
                    label_guess_warning = True

            # Check if the CONECT records are present and skip them
            # If they are not present activate the warning flag
            if line.startswith('CONECT'):
                conect_warning = False
            else:
                conect_warning = True

        # Overwrite the PDB file with the guessed atom labels in columns 77-78
        # if the labels are missing
        if label_guess_warning:

            with open(filename, 'w') as f:
                for line in pdbstr.strip().splitlines():
                    if line.startswith(('ATOM', 'HETATM')):
                        atom_name = line[12:16].strip()
                        atom_label = atom_name[0]
                        new_line = line[:76] + f"{atom_label:>2}" + line[78:] + '\n'
                        f.write(new_line)
                    else:
                        f.write(line + '\n')
            
            print('Warning: Atom labels were guessed based on atom names (first character).')
            print(f'Please verify the atom labels in the {filename} PDB file.')

        if conect_warning:

            # Create a molecule from the PDB file
            molecule = Molecule.read_PDB_file(filename)
            connectivity_matrix = molecule.get_connectivity_matrix()
            # Determine all the bonds in the molecule
            with open(filename, 'a') as f:
                for i in range(connectivity_matrix.shape[0]):
                    for j in range(i + 1, connectivity_matrix.shape[1]):
                        if connectivity_matrix[i, j] == 1:
                            # Convert indices to 1-based index for PDB format and ensure proper column alignment
                            i_index = i + 1
                            j_index = j + 1
                            # Align to the right 
                            con_string = "{:6s}{:>5d}{:>5d}".format('CONECT', i_index, j_index)
                            f.write(con_string + '\n')
       
            print('The CONECT records were not found in the PDB file.')
            print('The connectivity matrix was used to determine the bonds.')

    def _create_integrator(self):
        """
        Creates an OpenMM integrator object based on the specified ensemble type.

        Returns:
            OpenMM Integrator: Configured integrator for the simulation.
        """

        assert_msg_critical('openmm' in sys.modules, 'OpenMM is required for IMDatabasePointCollecter.')

        # Common parameters for Langevin integrators

        if self.ensemble in ['NVT', 'NPT']:
            integrator = mm.LangevinIntegrator(self.temperature, self.friction, self.timestep)
            integrator.setConstraintTolerance(1e-5)
            if self.ensemble == 'NPT':
                # Add a barostat for pressure control in NPT ensemble
                barostat = mm.MonteCarloBarostat(1 * unit.atmospheres, self.temperature, 25)
                self.system.addForce(barostat)
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
        tree = ET.ElementTree(ForceField)
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
            Name of the files to be generated. Default is 'qm_subregion'
        :param residue_name:
            Name of the residue. Default is 'QMR'
        """

        # Atoms belonging to the QM region do not have bonded parameters
        # Atoms belonging to the MM region have bonded parameters

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
            print('Linking atoms indices:', qm_connected)

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
        if write_xml:
            print('QM subregion (self.qm_atoms):', self.qm_atoms)
        self.mm_subregion = mm_atoms
        if write_xml:
            print('MM subregion:(self.mm_subregion)', self.mm_subregion)

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
        ImproperForce = ET.SubElement(ForceField, "PeriodicTorsionForce")
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
        # tree = ET.ElementTree(ForceField)
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
        Configures the system for QM/MM calculations.

        :param qm_atoms:
            List of atom indices to be included in the QM region.
        :param phase:
            Phase of the system ('gas', 'water', 'periodic').
        :param ff_gen:
            MMForceFieldGenerator object from VeloxChem.
        """

        assert_msg_critical('openmm' in sys.modules, 'OpenMM is required for IMDatabasePointCollecter.')

        # Set the QM/MM Interaction Groups
        total_atoms = self.system.getNumParticles()
        
        # The MM subregion is counted as regular MM atoms
        qm_group = set(self.qm_atoms)
        print('QM Group:', qm_group)
        mm_group = set(range(total_atoms)) - qm_group
        print('MM Group:', mm_group)
        if not mm_group:
            print('No external MM atoms found in the system')

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

        if not hasattr(self, "force_IM_embed"):
            self.force_IM_embed = mm.CustomExternalForce("-fx*x - fy*y - fz*z")
            self.force_IM_embed.addPerParticleParameter("fx")
            self.force_IM_embed.addPerParticleParameter("fy")
            self.force_IM_embed.addPerParticleParameter("fz")
            for a in qm_group:
                self.force_IM_embed.addParticle(a, (0.0, 0.0, 0.0))
            self.force_IM_embed.setForceGroup(9)
            self.system.addForce(self.force_IM_embed)

        # 2) Constant-force container for MM back-reaction (group 10)
        if mm_group and not hasattr(self, "force_MM_embed"):
            self.force_MM_embed = mm.CustomExternalForce("-fx*x - fy*y - fz*z")
            self.force_MM_embed.addPerParticleParameter("fx")
            self.force_MM_embed.addPerParticleParameter("fy")
            self.force_MM_embed.addPerParticleParameter("fz")
            for b in mm_group:
                self.force_MM_embed.addParticle(b, (0.0, 0.0, 0.0))
            self.force_MM_embed.setForceGroup(10)
            self.system.addForce(self.force_MM_embed)

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
                # DEBUG Only for testing purposes
                print('Using Reaction Field method for long-range electrostatics!')
                rfDielectric = nonbonded_force.getReactionFieldDielectric()
                krf = (1 / (self.cutoff**3)) * (rfDielectric - 1) / (2*rfDielectric + 1)
                crf = (1 /  self.cutoff) * (3*rfDielectric) / (2*rfDielectric + 1)
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
                vdw.addParticle([ff_gen.atoms[i]['sigma']* unit.nanometer, 
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
    
    def qm_stabilizer(self, ff_gen_qm):
        
        """
        Implements a MM potential to stabilize the QM region.
        The forces are 1% of the regular MM forces.

        :param qm_atoms: 
            List of atom indices to be included in the QM region.
        :param ff_gen_qm: 
            MMForceFieldGenerator object from VeloxChem.
        """

        assert_msg_critical('openmm' in sys.modules, 'OpenMM is required for IMDatabasePointCollecter.')

        # Harmonic bond contribution. Parameters are read from ff_gen_qm
        bonds = ff_gen_qm.bonds
        bond_force = mm.HarmonicBondForce()
        for bond, params in bonds.items():
            bond_force.addBond(*bond,
                            params['equilibrium'] * unit.nanometer,
                            params['force_constant'] * unit.kilojoule_per_mole / unit.nanometer**2 * self.scaling_factor)
        self.system.addForce(bond_force)

        # Harmonic angle contribution. Parameters are read from ff_gen_qm
        angles = ff_gen_qm.angles
        angle_force = mm.HarmonicAngleForce()
        for angle, params in angles.items():
            angle_force.addAngle(*angle,
                                params['equilibrium'] * np.pi / 180 * unit.radian,
                                params['force_constant'] * unit.kilojoule_per_mole / unit.radian**2 * self.scaling_factor)
        self.system.addForce(angle_force)

        # Periodic torsion contribution. Parameters are read from ff_gen_qm
        torsions = ff_gen_qm.dihedrals
        torsion_force = mm.PeriodicTorsionForce()
        rb_torsion_force = mm.RBTorsionForce()
        for torsion, params in torsions.items():
            if params['type'] == 'Fourier':
                torsion_force.addTorsion(*torsion,
                                         params['periodicity'],
                                         params['phase'] * np.pi / 180 * unit.radian,
                                         params['barrier'] * unit.kilojoule_per_mole * self.scaling_factor)
            elif params['type'] == 'RB':
                rb_torsion_force.addTorsion(*torsion,
                                            params['RB_coefficients'][0] * unit.kilojoule_per_mole * self.scaling_factor,
                                            params['RB_coefficients'][1] * unit.kilojoule_per_mole * self.scaling_factor,
                                            params['RB_coefficients'][2] * unit.kilojoule_per_mole * self.scaling_factor,
                                            params['RB_coefficients'][3] * unit.kilojoule_per_mole * self.scaling_factor,
                                            params['RB_coefficients'][4] * unit.kilojoule_per_mole * self.scaling_factor,
                                            params['RB_coefficients'][5] * unit.kilojoule_per_mole * self.scaling_factor)
        self.system.addForce(torsion_force)
        self.system.addForce(rb_torsion_force)

        # Improper torsion contribution. Parameters are read from ff_gen_qm
        impropers = ff_gen_qm.impropers
        improper_force = mm.PeriodicTorsionForce()
        for improper, params in impropers.items():
            improper_force.addTorsion(*improper,
                                    params['periodicity'],
                                    params['phase'] * np.pi / 180 * unit.radian,
                                    params['barrier'] * unit.kilojoule_per_mole * self.scaling_factor)
        self.system.addForce(improper_force)
    
    def update_gradient_and_energy(self, new_positions, collective_mode=None):
        ug_t0 = time()
        new_molecule = self._build_qm_molecule_from_positions(new_positions)
        self._latest_qm_molecule = new_molecule
        self._latest_qm_positions_nm = np.asarray(new_positions, dtype=np.float64).copy()
        if collective_mode is None:
            collective_mode = bool(self._mpi_is_active() and self.mpi_root_worker_mode)
        self._compute_all_roots_for_molecule(
            new_molecule,
            collective_mode=bool(collective_mode))
        self._postprocess_root_state_selection()
        self._add_runtime_timing('update_gradient_and_energy.total', time() - ug_t0)
        potential_kjmol = self.impes_drivers[self.current_state].impes_coordinate.energy * hartree_in_kcalpermol() * 4.184
        self.current_gradient = self.impes_drivers[self.current_state].impes_coordinate.gradient
        self.current_energy = potential_kjmol
        return self.current_gradient, potential_kjmol
    
    def _build_qm_molecule_from_positions(self, positions):
        positions_ang = positions * 10 
        atom_labels = getattr(self, '_topology_atom_labels', None)
        if atom_labels is None:
            atom_labels = [atom.element.symbol for atom in self.topology.atoms()]
            self._topology_atom_labels = atom_labels

        new_molecule = None
        # Check if there is a QM/MM partition in the system
        if self.mm_subregion is not None:
            # Create a molecule with a link atom (H)
            # The linking atom is defined in self.linking_atoms
            # It need to be changed to a H atom at 1.0 angstrom
            # The QM region is defined in self.qm_atoms

            qm_atom_labels = []
            for i in self.qm_atoms:
                # Change the linking atom to H
                if i in self.linking_atoms:
                    qm_atom_labels.append('H')
                else:
                    qm_atom_labels.append(atom_labels[i])

            
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
            # Atom labels for the QM region
            qm_atom_labels = [atom_labels[i] for i in self.qm_atoms]
            new_molecule = Molecule(qm_atom_labels, positions_ang, units="angstrom")
            new_molecule.set_charge(self.molecule.get_charge())
            new_molecule.set_multiplicity(self.molecule.get_multiplicity())
            if not (self._mpi_is_active() and self.mpi_root_worker_mode and (not self._mpi_is_root())):
                self.unique_molecules.append(new_molecule)

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

        molecule = getattr(self, '_latest_qm_molecule', None)
        if molecule is not None:
            if qm_positions_nm is None:
                return molecule
            cached_positions = getattr(self, '_latest_qm_positions_nm', None)
            current_positions = np.asarray(qm_positions_nm, dtype=np.float64)
            if cached_positions is not None and cached_positions.shape == current_positions.shape:
                if np.allclose(cached_positions, current_positions, atol=1.0e-12, rtol=0.0):
                    return molecule

        if qm_positions_nm is None:
            return None

        return self._build_qm_molecule_from_positions(
            np.asarray(qm_positions_nm, dtype=np.float64))
        
    
    def _compute_all_roots_for_molecule(self, molecule, collective_mode=None):
        if collective_mode is None:
            collective_mode = bool(self._mpi_is_active() and self.mpi_root_worker_mode)
        collective_mode = bool(collective_mode)

        for root in self.roots_to_follow:
            driver = self.impes_drivers[root]
            self._ensure_interp_engine_comm(driver)
            if driver.qm_data_points is not self.qm_data_point_dict[root]:
                driver.qm_data_points = self.qm_data_point_dict[root]
            driver.prepare_runtime_data_cache()

            preload_config = bool(getattr(driver, '_mpi_preload_enabled_config', getattr(driver, 'mpi_preload_enabled', False)))

            # Enable preload only while root+workers are in the synchronized
            # compute phase. Outside this phase, keep it disabled.
            if self._mpi_is_active() and self.mpi_root_worker_mode:
                driver.mpi_preload_enabled = bool(preload_config and collective_mode)
            else:
                driver.mpi_preload_enabled = preload_config

            driver.mpi_collective_compute_active = collective_mode
            compute_t0 = time()
            try:
                driver.compute(molecule)
               
            finally:
                driver.mpi_collective_compute_active = False
                if self._mpi_is_active() and self.mpi_root_worker_mode:
                    driver.mpi_preload_enabled = False
                else:
                    driver.mpi_preload_enabled = preload_config
      
            compute_dt = time() - compute_t0
            self._add_runtime_timing(
                f'update_gradient_and_energy.impes_compute_root_{root}',
                compute_dt)
            self._add_runtime_timing(
                'update_gradient_and_energy.impes_compute_total',
                compute_dt)
    
    def _postprocess_root_state_selection(self):
    
        self.current_state = self.roots_to_follow[0]

    def _mpi_worker_compute_step(self, qm_positions_nm, sync_roots):
        if self.mpi_reload_from_hdf5:
            for root in sync_roots:
                self._reload_interpolation_root_from_hdf5(root, self.inv_sqrt_masses)

        # qm_positions_nm is shape (n_qm_atoms, 3)
        new_molecule = self._build_qm_molecule_from_positions(np.asarray(qm_positions_nm, dtype=np.float64))
        self._compute_all_roots_for_molecule(new_molecule, collective_mode=True)


    def update_gradient(self, new_positions, collective_mode=None):
        """
        Updates and returns the gradient of the QM region.

        :param new_positions:
            The new positions of the atoms in the QM region.
        :return:
            The gradient of the QM region.
        """
        gradient, _ = self.update_gradient_and_energy(
            new_positions,
            collective_mode=collective_mode)

        return gradient

    def update_forces(self, context):
        """
        Updates the forces in the system based on a new gradient.

        Args:
            context: The OpenMM context object.
        """

        assert_msg_critical('openmm' in sys.modules, 'OpenMM is required for IMDatabasePointCollecter.')

        if not hasattr(self, '_gradient_to_force_factor'):
            self._gradient_to_force_factor = (
                (4.184 * hartree_in_kcalpermol() * 10.0 / bohr_in_angstrom())
                * unit.kilojoule_per_mole / unit.nanometer
            )
        conversion_factor = self._gradient_to_force_factor

        qm_force = getattr(self, '_qm_external_force', None)
        if qm_force is None:
            qm_force = self.system.getForce(self.qm_force_index)
            self._qm_external_force = qm_force

        # Update the forces of the QM region.
        qm_positions = self._extract_qm_positions_nm(context)

        self.velocities_np.append(context.getState(getVelocities=True).getVelocities(True))
        gradient_t0 = time()
        gradient = self.update_gradient(qm_positions)
        self._add_runtime_timing('update_forces.update_gradient', time() - gradient_t0)

        # Reuse the exact molecule passed through interpolation compute when available.
        new_molecule = self._get_latest_qm_molecule(qm_positions)

        force = -np.array(gradient) * conversion_factor

        self.all_gradients.append(gradient)
        ############################################################
        #################### Correlation Check #####################
        ############################################################
        allowed = True 
        self.add_a_point = False

        if self.density_around_data_point[1] is not None and self.current_state == self.density_around_data_point[2]:
            current_dihedral = new_molecule.get_dihedral_in_degrees([self.density_around_data_point[1][0], self.density_around_data_point[1][1], self.density_around_data_point[1][2], self.density_around_data_point[1][3]])%360.0
            lower, upper = self.allowed_molecule_deviation[self.current_state][self.density_around_data_point[1]][self.density_around_data_point[3]]

            # Case 1: If boundaries do not wrap (e.g., [-60, 60])
            if lower < upper:
                allowed = lower <= current_dihedral <= upper
            else:
                # Case 2: If boundaries wrap around (e.g., [120, -120])
                allowed = current_dihedral >= lower or current_dihedral <= upper

        if allowed:

            openmm_coordinate = context.getState(getPositions=True).getPositions()
            self.coordinates.append(openmm_coordinate)
            self.velocities.append(context.getState(getVelocities=True).getVelocities())
            
            if self.skipping_value == 0 and self.point_checker + 1 > self.start_collect:
                bond_rmsd_values = np.asarray(self.impes_drivers[self.current_state].bond_rmsd, dtype=np.float64)
                angle_rmsd_values = np.asarray(self.impes_drivers[self.current_state].angle_rmsd, dtype=np.float64)
                dihedral_rmsd_values = np.asarray(self.impes_drivers[self.current_state].dihedral_rmsd, dtype=np.float64)

                mean_bond_rmsd = float(np.mean(bond_rmsd_values)) if bond_rmsd_values.size else 0.0
                mean_angle_rmsd = float(np.mean(angle_rmsd_values)) if angle_rmsd_values.size else 0.0
                mean_dihedral_rmsd = float(np.mean(dihedral_rmsd_values)) if dihedral_rmsd_values.size else 0.0

                if mean_bond_rmsd > 0.01 or mean_angle_rmsd > 0.01 or mean_dihedral_rmsd > 1.0:

                    K = 5  # Number of closest previous matches to cache
                    threshold = 5e-2 #self.distance_thrsh - 0.05
                    new_coords = new_molecule.get_coordinates_in_bohr()
                    n_atoms = len(new_molecule.get_labels())
                    sym_group = self.non_core_symmetry_groups['gs' if self.current_state == 0 or self.drivers['es'] is not None and self.drivers['es'][0].spin_flip and self.current_state == 1 else 'es'][4]
                    molecule_list = self.allowed_molecules[self.current_state]['molecules']

                
                    if not hasattr(self, "previous_candidate_indices"):
                        self.previous_candidate_indices = []

                    scanned = False
                    closest_indices = []
                    index_added = []

                    # --- 1. Try cached candidates first
                    for idx in self.previous_candidate_indices:
                        if idx >= len(molecule_list):
                            continue
                        checked_molecule = molecule_list[idx]
                        checked_coords = checked_molecule.get_coordinates_in_bohr()
                        checked_distance = self.cartesian_just_distance(checked_coords, new_coords, sym_group)
                        normed_dist = (np.linalg.norm(checked_distance) / np.sqrt(n_atoms)) * bohr_in_angstrom()

                        if idx not in index_added:
                            closest_indices.append((idx, normed_dist))
                            index_added.append(idx)
                        if normed_dist <= threshold:
                            scanned = True
                            break

                    # --- 2. Fallback: full scan only if no match
                    if not scanned:
                        for idx, checked_molecule in enumerate(molecule_list):
                            checked_coords = checked_molecule.get_coordinates_in_bohr()
                            checked_distance = self.cartesian_just_distance(checked_coords, new_coords, sym_group)

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

                            if (np.linalg.norm(length_vectors) / np.sqrt(len(self.molecule.get_labels()))) * bohr_in_angstrom() <= 5e-2:#abs(self.distance_thrsh - 0.5):
                                self.add_a_point = False
                                break

                            if i == len(self.qm_data_point_dict[self.current_state]):
                                self.add_a_point = True

            else:
                self.skipping_value -= 1
                if self.skipping_value < 0:
                    self.skipping_value = 0
            
            if self.step % 50 == 0 and self.step > self.start_collect:
                self.add_a_point = True


            self.point_checker += 1 
            self.last_gpr_addition += 1

            if self.add_a_point == True:
                pcorr_t0 = time()
                self.point_correlation_check(new_molecule)
                self._add_runtime_timing('update_forces.point_correlation_check', time() - pcorr_t0)
                self.last_gpr_addition = 0
            if self.point_checker == 0:            
                
                self.last_force_point = 0
                self.swap_back = True
                self.current_state = self.starting_state
                context.setPositions(self.coordinates[0])
                self.coordinates = [self.coordinates[0]]
                # context.setVelocities(self.velocities[0])

                if self._metadynamics is not None and self.metadynamics_enabled:
                    self._reset_metadynamics_bias(context, reason='datapoint_added')


            # Set initial velocities if the ensemble is NVT or NPT
                if self.ensemble in ['NVT', 'NPT']:

                    context.setVelocitiesToTemperature(self.temperature * 0.5)
                else:
                    context.setVelocities(self.start_velocities)
                self.velocities = [context.getState(getVelocities=True).getVelocities()]
                
                qm_positions = self._extract_qm_positions_nm(context)
                gradient_2 = self.update_gradient(qm_positions, collective_mode=False)
                new_molecule = self._get_latest_qm_molecule(qm_positions)
                
                force = -np.array(gradient_2) * conversion_factor
            
            if (self.point_checker + self.last_point_added) % self.duration == 0 and self.point_checker != 0:
                print('start when last point was added', self.point_checker + self.last_point_added)
                self.last_point_added = 0
                self.unadded_cycles -= 1

            if self.point_checker < 500 and self.cycle_iteration != self.unadded_cycles:
                self.unadded_cycles += 1
            
            self.coordinates_xyz.append(qm_positions * 10)

            ############################################################
            self.state_specific_molecules[self.current_state].append(new_molecule)
            for i, atom_idx in enumerate(self.qm_atoms):
          
                qm_force.setParticleParameters(i, atom_idx, force[i])
            qm_force.updateParametersInContext(context)


            for root in self.roots_to_follow:
                self.gloabal_sim_informations[f'state_{root}']['pot_energies'].append(self.impes_drivers[root].get_energy() * hartree_in_kcalpermol())
                self.gloabal_sim_informations[f'state_{root}']['gradients'].append(self.impes_drivers[root].get_gradient() * hartree_in_kcalpermol() * (1.0 / bohr_in_angstrom()))
                
        
        else:
            context.setPositions(self.coordinates[0])
            self.coordinates = [self.coordinates[0]]
            if self.ensemble in ['NVT', 'NPT']:

                context.setVelocitiesToTemperature(self.temperature * 0.5)
            else:
                context.setVelocities(self.start_velocities)
            self.velocities = [context.getState(getVelocities=True).getVelocities()]
            qm_positions = self._extract_qm_positions_nm(context)
            gradient_2 = self.update_gradient(qm_positions, collective_mode=False)
            new_molecule = self._get_latest_qm_molecule(qm_positions)
            force = -np.array(gradient_2) * conversion_factor
            
            self.state_specific_molecules[self.current_state].append(new_molecule)
            self.coordinates_xyz.append(qm_positions * 10)
            for i, atom_idx in enumerate(self.qm_atoms):
                qm_force.setParticleParameters(i, atom_idx, force[i])
            qm_force.updateParametersInContext(context)
        
    ####################################################################
    ################ Functions to expand the database ##################
    ####################################################################

    def _sampling_screen_check(self, molecule):

        if not self.sampling_enabled or self.sampling_driver is None or self.sampling_impes_drivers is None:
            return {"skip_full_qm": False, "details": {}}
        
        natms = len(molecule.get_labels())

        sampling_qm, sampling_grad, _ = self.sampling_driver['gs']

        xtb_energy, _, _ = self._compute_energy_mpi_safe(sampling_qm, molecule, basis=None)
        xtb_grad = self._compute_gradient_mpi_safe(sampling_grad, molecule, basis=None, scf_results=None, rsp_results=None)

        e_thr = float(self.sampling_settings["e_thrsh_kcal_per_atom"])
        g_thr = float(self.sampling_settings["g_rmsd_thrsh_kcal_ang_per_atom"])
        c_thr = float(self.sampling_settings["force_orient_cos"])

        all_ok = True
        details = {}
        for root in self.roots_to_follow:
            self._interp_compute_root_local_serial(self.sampling_impes_drivers[root], molecule)
            e_im = self.sampling_impes_drivers[root].impes_coordinate.energy
            g_im = self.sampling_impes_drivers[root].impes_coordinate.gradient

            e_diff = abs(float(xtb_energy[0]) - float(e_im)) / natms * hartree_in_kcalpermol()
            g_diff = (xtb_grad[0] - g_im) * hartree_in_kcalpermol() * bohr_in_angstrom()
            g_rmsd = float(np.sqrt((g_diff ** 2).mean()))

            gq = xtb_grad[0].ravel()
            gi = g_im.ravel()
            denom = np.linalg.norm(gq) * np.linalg.norm(gi)
            cos_theta = 1.0 if denom < 1.0e-15 else float(np.dot(gq, gi) / denom)
            print('xtb e_diff', e_diff)
            details[root] = {"e_diff_kcal_per_atom": e_diff, "g_rmsd": g_rmsd, "cos": cos_theta}
            if not (e_diff <= e_thr and g_rmsd <= g_thr and cos_theta >= c_thr):
                all_ok = False

        if all_ok:
            self.previous_energy_list.append(details[root].get('e_diff_kcal_per_atom'))
            curr_state_diff = details[root].get('e_diff_kcal_per_atom')
            if curr_state_diff == 0:
                curr_state_diff = 1e-8
            if len(self.previous_energy_list) > 3:
                # Compute gradient (rate of change of energy difference)
                grad1 = self.previous_energy_list[-2] - self.previous_energy_list[-3]
                grad2 = self.previous_energy_list[-1] - self.previous_energy_list[-2]
                # Base skipping value calculation
                
                base_skip = min(round(abs(self.energy_threshold / (curr_state_diff)**2)), 20) - 1
                # Adjust skipping value based on gradient
                if grad2 > grad1:  # Energy difference is increasing
                    self.skipping_value = max(1, base_skip - 1)  # Reduce skipping for more frequent checks
                else:  # Energy difference is decreasing
                    self.skipping_value = base_skip + 10  # Increase skipping to check less often
                print(f"Energy Difference: {(details[root].get('e_diff_kcal_per_atom')):.6f}, Gradient: {grad2 - grad1:.6f}, Skipping Value: {self.skipping_value}")
            else:
                self.skipping_value = min(round(abs(self.energy_threshold / (curr_state_diff)**2)), 20)
            print('skipping value', self.skipping_value)
            return {"skip_full_qm": True, "details": details}
        return {"skip_full_qm": False, "details": details}

    def point_correlation_check(self, molecule):
        """ Takes the current point on the PES and checks with a QM-energy
            calculation is necessary based on the current difference to the
            interpolation. Based on the difference the step_size of the next
            step is determined.

            :param molecule:
                The molecule corresponding to the current conformation
                in the IM dynamics.
        """
        
        if self.sampling_enabled:
            screen = self._sampling_screen_check(molecule)
            if screen["skip_full_qm"]:
                self.add_a_point = False
                return
        
        current_state_difference = {}
        state_specific_energies = {}
        state_specific_gradients = {}
        current_state_to_state_difference = {}
        self.add_a_point = False

        all_combinations = list(itertools.combinations(self.roots_to_follow, 2))

        for comb in all_combinations:
            current_state_to_state_difference[comb] = [None, None]
        for root in self.roots_to_follow:
            state_specific_energies[root] = None
            state_specific_gradients[root] = None
            current_state_difference[root] = [None, None, None]

        addition_of_state_specific_points = []

        for state, key in enumerate(self.drivers.keys()):

            drivers = None
            identification_state = state
            
            if state == 0 and self.drivers[key] is not None and 0 in self.roots_to_follow:
                drivers = self.drivers[key]

            mask = None
            if state == 0 and self.drivers['gs'] is not None and any(x == 0 for x in self.roots_to_follow):
                mask = [0] 
            elif state == 1 and self.drivers['es'] is not None and any(x > 0 for x in self.roots_to_follow):
                drivers = self.drivers[key]
                mask = [x for x in self.roots_to_follow if x != 0]
                if 0 not in self.roots_to_follow:
                    identification_state = 0
            else:
                continue
                

            natms = len(molecule.get_labels())
            print('current test mol \n', molecule.get_xyz_string())
            print('############# Energy is QM claculated ############')
            print('identifaction state', identification_state, 'current state', self.current_state)
            current_basis = None
            if identification_state == 0:
                current_basis = MolecularBasis.read(molecule, self.basis_set_label['gs'])
            else:
                current_basis = MolecularBasis.read(molecule, self.basis_set_label['es'])
            
   
            # qm_energy, scf_results, rsp_results = self._compute_energy_mpi_safe(drivers[0], molecule, current_basis)
            
            # if isinstance(drivers[0], LinearResponseEigenSolver) or isinstance(drivers[0], TdaEigenSolver):
            #     qm_energy = qm_energy[mask]


            # gradients = self._compute_gradient_mpi_safe(drivers[1], molecule, current_basis, scf_results, rsp_results)

            # print(qm_energy, gradients)

            qm_energy, gradients = self._aux_eval_eg(
                molecule=molecule,
                basis_label=current_basis.get_main_basis_label(),
                driver_key=key,
                state_mask=mask,
            )

            for e_idx in range(len(qm_energy)):

                grad = gradients[e_idx]         # (3N,)
                energy_difference = (abs(qm_energy[e_idx] - self.impes_drivers[self.roots_to_follow[identification_state + e_idx]].impes_coordinate.energy))


                gradient_difference = (grad - self.impes_drivers[self.roots_to_follow[identification_state + e_idx]].impes_coordinate.gradient) * hartree_in_kcalpermol() * bohr_in_angstrom()
                rmsd_gradient    = np.sqrt((gradient_difference**2).mean())
                
                cos_theta = 1.0
                if rmsd_gradient > 0.5 and energy_difference > 0.01:
                    
                    gq = grad.ravel()
                    gi = self.impes_drivers[self.roots_to_follow[identification_state + e_idx]].impes_coordinate.gradient.ravel()

                    cos_theta = 0.0
                    significant_mask = (np.abs(gq) > 1e-4) | (np.abs(gi) > 1e-4)
                    if np.sum(significant_mask) == 0:
                        cos_theta = 1.0
                    # print('significant components', np.sum(significant_mask))
                    gq_clean = gq[significant_mask]
                    gi_clean = gi[significant_mask]

                    norm_gq_clean = np.linalg.norm(gq_clean)
                    norm_gi_clean = np.linalg.norm(gi_clean)
                    
                    if norm_gq_clean == 0 and norm_gi_clean == 0:
                        cos_theta = 1.0
                    
                    cos_theta = np.dot(gq_clean, gi_clean) / (norm_gq_clean * norm_gi_clean)
                
                print('Energy difference', energy_difference, energy_difference * hartree_in_kcalpermol(), 'kcal/mol', 'energy differences rmsd', energy_difference / natms * hartree_in_kcalpermol())
                print('gradients alignment', cos_theta, 'rmsd gradient', rmsd_gradient, 'kcal/mol/angstrom')

                state_specific_energies[self.roots_to_follow[identification_state + e_idx]] = [qm_energy[e_idx], self.impes_drivers[self.roots_to_follow[identification_state + e_idx]].impes_coordinate.energy]
                current_state_difference[self.roots_to_follow[identification_state + e_idx]][0] = energy_difference / natms * hartree_in_kcalpermol()
                current_state_difference[self.roots_to_follow[identification_state + e_idx]][1] = rmsd_gradient
                current_state_difference[self.roots_to_follow[identification_state + e_idx]][2] = cos_theta
                state_specific_gradients[self.roots_to_follow[identification_state + e_idx]] = [grad, self.impes_drivers[self.roots_to_follow[identification_state + e_idx]].impes_coordinate.gradient]
        
        print(current_state_difference[self.current_state][0] > self.energy_threshold and not self.use_opt_confidence_radius[0]
            or current_state_difference[self.current_state][0] > self.energy_threshold and self.use_opt_confidence_radius[0] and self.confidence_radius_optimized
            or current_state_difference[self.current_state][1] > self.gradient_rmsd_thrsh and not self.use_opt_confidence_radius[0]
            or current_state_difference[self.current_state][1] > self.gradient_rmsd_thrsh and self.use_opt_confidence_radius[0] and self.confidence_radius_optimized)
        if (current_state_difference[self.current_state][0] > self.energy_threshold and not self.use_opt_confidence_radius[0]
            or current_state_difference[self.current_state][0] > self.energy_threshold and self.use_opt_confidence_radius[0] and self.confidence_radius_optimized
            or current_state_difference[self.current_state][1] > self.gradient_rmsd_thrsh and not self.use_opt_confidence_radius[0]
            or current_state_difference[self.current_state][1] > self.gradient_rmsd_thrsh and self.use_opt_confidence_radius[0] and self.confidence_radius_optimized):
            print('The point would be added due to the thresholds')
            addition_of_state_specific_points.append(self.current_state)

        for comb in all_combinations:
            current_state_to_state_difference[comb][0] = abs((state_specific_energies[comb[0]][0] - state_specific_energies[comb[1]][0]) * hartree_in_kjpermol())
            current_state_to_state_difference[comb][1] = abs((state_specific_energies[comb[0]][1] - state_specific_energies[comb[1]][1]) * hartree_in_kjpermol())
            print('current state to state', comb, current_state_to_state_difference[comb])
            if current_state_to_state_difference[comb][0] < 30.0 and abs(current_state_to_state_difference[comb][1] - current_state_to_state_difference[comb][0]) > 10.0:
                if comb[0] not in addition_of_state_specific_points:
                    addition_of_state_specific_points.append(comb[0])
                if comb[1] not in addition_of_state_specific_points:
                    addition_of_state_specific_points.append(comb[1])

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
                                              
            if self.add_a_point:    
            
                if not self.expansion:
                    length_vectors = (self.impes_drivers[-1].impes_coordinate.cartesian_distance_vector(self.qm_data_points[0]))
                    rmsd = (np.linalg.norm(length_vectors) / np.sqrt(len(self.molecule.get_labels())) * bohr_in_angstrom())
                    self.expansion_molecules.append((molecule, energy_difference * hartree_in_kcalpermol(), rmsd, (np.linalg.norm(length_vectors))))
                    self.last_point_added = self.point_checker - 1
                    self.point_checker = 0
        else:
            for root in self.roots_to_follow:
                self.allowed_molecules[root]['molecules'].append(molecule)
                self.allowed_molecules[root]['im_energies'].append(self.impes_drivers[root].impes_coordinate.energy)
                self.allowed_molecules[root]['qm_energies'].append(state_specific_energies[root][0])
                self.allowed_molecules[root]['qm_gradients'].append(state_specific_gradients[root][0])

                self.write_qm_energy_determined_points_h5(self.allowed_molecules[root]['molecules'][self.last_added: ],
                                                    self.allowed_molecules[root]['qm_energies'][self.last_added:],
                                                    self.allowed_molecules[root]['qm_gradients'][self.last_added:],
                                                    self.allowed_molecules[root]['im_energies'][self.last_added:],
                                                    root)
                self.last_added = len(self.allowed_molecules[root]['molecules'])

            if self.use_opt_confidence_radius[0] and len(self.allowed_molecules[self.current_state]['molecules']) >= 10 and self.density_around_data_point[0][self.current_state] > 1 and self.density_around_data_point[0][self.current_state] % 1 == 0 and self.prev_dens_of_points[self.current_state] != self.density_around_data_point[0][self.current_state]:
                self.prev_dens_of_points[self.current_state] = self.density_around_data_point[0][self.current_state]        
                trust_radius = None
                sym_dict = self.non_core_symmetry_groups['gs']
                if self.current_state > 0 or root == 1 and self.drivers['es'] is not None and self.drivers['es'][0].spin_flip:
                    sym_dict = self.non_core_symmetry_groups['es']
                 
                for run_through in range(1):
                    indices = [i for i in range(len(self.allowed_molecules[self.current_state]['molecules']) - 1)]
                    chosen_structures = [self.allowed_molecules[self.current_state]['molecules'][idx] for idx in indices]
                    chosen_qm_energies = [self.allowed_molecules[self.current_state]['qm_energies'][idx] for idx in indices]
                    chosen_im_energies = [self.allowed_molecules[self.current_state]['im_energies'][idx] for idx in indices]
                    chosen_qm_gradients = [self.allowed_molecules[self.current_state]['qm_gradients'][idx] for idx in indices]
                    
                        
                    if self.use_opt_confidence_radius[1] == 'multi_grad':
                        
                        trust_radius = self.determine_trust_radius_gradient(chosen_structures, 
                                                                chosen_qm_energies,
                                                                chosen_qm_gradients,
                                                                chosen_im_energies, 
                                                                self.qm_data_point_dict[self.current_state], 
                                                                self.interpolation_settings[self.current_state],
                                                                self.qm_symmetry_datapoint_dict[self.current_state],
                                                                sym_dict,
                                                                self.root_z_matrix[self.current_state],
                                                                exponent_p_q = (self.impes_drivers[self.current_state].exponent_p, self.impes_drivers[self.current_state].exponent_q))
                    
                    
                    elif self.use_opt_confidence_radius[1] == 'bayes':
                        trust_radius = self.determine_beysian_trust_radius(self.allowed_molecules[self.current_state]['molecules'], 
                                                                        self.allowed_molecules[self.current_state]['qm_energies'], 
                                                                        self.qm_data_point_dict[self.current_state], 
                                                                        self.interpolation_settings[self.current_state], 
                                                                        self.qm_symmetry_datapoint_dict[self.current_state],
                                                                        sym_dict)
                        
      
                    for idx, trust_radius in enumerate(trust_radius):
                        print(self.sorted_state_spec_im_labels[self.current_state][idx])
                        self.qm_data_point_dict[self.current_state][idx].update_confidence_radius(self.interpolation_settings[self.current_state]['imforcefield_file'], self.sorted_state_spec_im_labels[self.current_state][idx], trust_radius)
                        self.qm_data_point_dict[self.current_state][idx].confidence_radius = trust_radius
                        if self.sampling_enabled:
                            self.sampling_qm_data_point_dict[self.current_state][idx].update_confidence_radius(self.sampling_interpolation_settings[self.current_state]['imforcefield_file'], self.sorted_state_spec_im_labels[self.current_state][idx], trust_radius)
                            self.sampling_qm_data_point_dict[self.current_state][idx].confidence_radius = trust_radius

                    
                    self._refresh_interpolation_driver_caches(self.current_state)

                self.confidence_radius_optimized = True
                self.add_a_point = False
        
        self._emit_test_hook("point_correlation_decision", {
            "current_state": int(self.current_state),
            "energy_threshold_kcal_per_atom": float(self.energy_threshold),
            "gradient_threshold": float(self.gradient_rmsd_thrsh),
            "force_orient_threshold": float(self.force_orient_thrsh),
            "add_a_point": bool(self.add_a_point),
            "state_diff_kcal_per_atom": float(current_state_difference[self.current_state][0]),
            "state_grad_rmsd": float(current_state_difference[self.current_state][1]),
            "state_force_orient_cos": float(current_state_difference[self.current_state][2]),
        })

        # calcualte energy gradient
        if self.add_a_point is False:
            self.previous_energy_list.append(current_state_difference[self.current_state][0])
            curr_state_diff = (current_state_difference[self.current_state][0])
            if curr_state_diff == 0:
                curr_state_diff = 1e-8
            if len(self.previous_energy_list) > 3:
                # Compute gradient (rate of change of energy difference)
                grad1 = self.previous_energy_list[-2] - self.previous_energy_list[-3]
                grad2 = self.previous_energy_list[-1] - self.previous_energy_list[-2]
                # Base skipping value calculation
                
                base_skip = min(round(abs(self.energy_threshold / (curr_state_diff)**2)), 20) - 1
                # Adjust skipping value based on gradient
                if grad2 > grad1:  # Energy difference is increasing
                    self.skipping_value = max(1, base_skip - 1)  # Reduce skipping for more frequent checks
                else:  # Energy difference is decreasing
                    self.skipping_value = base_skip + 10  # Increase skipping to check less often
                print(f"Energy Difference: {(current_state_difference[self.current_state][0]):.6f}, Gradient: {grad2 - grad1:.6f}, Skipping Value: {self.skipping_value}")
            else:
                self.skipping_value = min(round(abs(self.energy_threshold / (curr_state_diff)**2)), 20)
            print('len of molecules', len(self.allowed_molecules[self.current_state]['molecules']), self.use_opt_confidence_radius)
        
        self.skipping_value = 0
        if self.add_a_point and self.expansion:
            self.confidence_radius_optimized = False
            print('✨ A point is added! ✨', self.point_checker)
            print(molecule.get_xyz_string())

            ############# Implement constraint optimization ############
            
            state_specific_molecules = []
            
            imp_int_coord = None
            opt_results = None
            
            if self.identfy_relevant_int_coordinates[0]:  
                for state_to_optim in addition_of_state_specific_points:   
                    current_basis = None
                    if state_to_optim  == 0:
                        current_basis = MolecularBasis.read(molecule, self.basis_set_label['gs'])
                    else:
                        current_basis = MolecularBasis.read(molecule, self.basis_set_label['es'])
                                
                    optimized_molecule = molecule
                    self._interp_compute_root_local_serial(self.impes_drivers[state_to_optim], molecule)
                    current_weights = self.impes_drivers[state_to_optim].weights
                    weights = [value for _, value in current_weights.items()]
                    used_labels = [label_idx for label_idx, _ in current_weights.items()]
                    # Sort labels and weights by descending weight
                    sorted_items = sorted(zip(used_labels, weights), key=lambda x: x[1], reverse=True)
                    total_weight = sum(weights)
                    cumulative_weight = 0.0
                    internal_coordinate_datapoints = []
                    for label, weight in sorted_items:
                        cumulative_weight += weight
                        internal_coordinate_datapoints.append((self.qm_data_point_dict[state_to_optim][label], weight))
                        if cumulative_weight >= 0.8 * total_weight:
                            break
                    # qm_datapoints_weighted = [qm_datapoint for qm_datapoint in enumerate if ]
                    print('Items', sorted_items, len(internal_coordinate_datapoints), state_specific_energies[state_to_optim][0])
                    _, candidate_constraints, _ = self.impes_drivers[state_to_optim].determine_important_internal_coordinates(state_specific_gradients[state_to_optim][0], molecule, self.root_z_matrix[state_to_optim], internal_coordinate_datapoints)
                    corr_new_dp_distrib = False

                    imp_coord_constraint = candidate_constraints.copy()
                    main_constraint_list = candidate_constraints.copy()
                    opt_results = None
                    added_bond_key = []
                    for dihedral in self.root_z_matrix[state_to_optim]['dihedrals']:
                        for rot_bond in self.impes_drivers[state_to_optim].symmetry_information[5]:
                            bond_key = tuple(sorted([dihedral[1], dihedral[2]]))
                            if bond_key in added_bond_key:
                                continue
                            if (
                                set([dihedral[1], dihedral[2]]) == set(rot_bond)
                                and bond_key not in self.impes_drivers[state_to_optim].symmetry_information[7][3]
                            ):
                                if dihedral not in main_constraint_list:
                                    added_bond_key.append(bond_key)
                                    main_constraint_list.append(dihedral)
                                   
                    print('Main CONSTRAINTS', main_constraint_list)

                    while not corr_new_dp_distrib:

                        # if state_to_optim == 0 and isinstance(self.drivers['gs'][0], ScfRestrictedDriver) or state_to_optim == 0 and isinstance(self.drivers['gs'][0], ScfUnrestrictedDriver):
                            # _, scf_results, rsp_results = self._compute_energy_mpi_safe(self.drivers['gs'][0], molecule, current_basis)
        
                            # optimized_molecule, opt_results = self._run_optimization_mpi_safe(
                            #     self.drivers['gs'][0],
                            #     molecule,
                            #     constraints=main_constraint_list,
                            #     index_offset=1,
                            #     compute_args=(current_basis, scf_results),
                            #     collective=True,
                            #     phase_name='optimization tag',
                            # )

                        optimized_molecule, opt_results = self._aux_optimize(
                            state_to_optim=state_to_optim,
                            molecule=molecule,
                            basis_label=current_basis.get_main_basis_label(),
                            constraints=main_constraint_list,
                        )

                        optimized_molecule = Molecule.from_xyz_string(opt_results['final_geometry'])
                        optimized_molecule.set_charge(molecule.get_charge())
                        optimized_molecule.set_multiplicity(molecule.get_multiplicity())
                        current_basis = MolecularBasis.read(optimized_molecule, current_basis.get_main_basis_label())

                        print('Optimized Molecule', optimized_molecule.get_xyz_string(), '\n\n', molecule.get_xyz_string())
                        
                        # elif state_to_optim == 0 and isinstance(self.drivers['gs'][0], XtbDriver):

                        #     # optimized_molecule, opt_results = self._run_optimization_mpi_safe(
                        #     #     self.drivers['gs'][0],
                        #     #     molecule,
                        #     #     constraints=main_constraint_list,
                        #     #     index_offset=1,
                        #     #     collective=True,
                        #     #     phase_name='optimization tag',
                        #     # )
                        #     optimized_molecule, opt_results = self._aux_optimize(
                        #         state_to_optim=state_to_optim,
                        #         molecule=molecule,
                        #         basis_label=current_basis.get_main_basis_label(),
                        #         constraints=main_constraint_list,
                        #     )

                        #     optimized_molecule = Molecule.from_xyz_string(opt_results['final_geometry'])
                        #     optimized_molecule.set_charge(molecule.get_charge())
                        #     optimized_molecule.set_multiplicity(molecule.get_multiplicity())

                        
                        # elif state_to_optim >= 1 and isinstance(self.drivers['es'][0], LinearResponseEigenSolver) or state_to_optim >= 1 and isinstance(self.drivers['es'][0], TdaEigenSolver):
                                
                        #     self.drivers['es'][1].state_deriv_index = [state_to_optim]
                        #     opt_drv = OptimizationDriver(self.drivers['es'][1])
                        #     _, _, rsp_results = self._compute_energy_mpi_safe(self.drivers['es'][0], molecule, current_basis)
                        #     opt_drv.ostream.mute()
                            
                        #     opt_constraint_list = []

                        #     for constraint in main_constraint_list:
                        #         if len(constraint) == 2:
                        #             opt_constraint = f"freeze distance {constraint[0] + 1} {constraint[1] + 1}"
                        #             opt_constraint_list.append(opt_constraint)
                                
                        #         elif len(constraint) == 3:
                        #             opt_constraint = f"freeze angle {constraint[0] + 1} {constraint[1] + 1} {constraint[2] + 1}"
                        #             opt_constraint_list.append(opt_constraint)
                            
                        #         else:
                        #             opt_constraint = f"freeze dihedral {constraint[0] + 1} {constraint[1] + 1} {constraint[2] + 1} {constraint[3] + 1}"
                        #             opt_constraint_list.append(opt_constraint)
                            
                        #     for constraint in self.identfy_relevant_int_coordinates[1]:
                        #         if constraint in opt_constraint_list:
                        #             continue
                        #         if len(constraint) == 2:
                        #             opt_constraint = f"freeze distance {constraint[0] + 1} {constraint[1] + 1}"
                        #             opt_constraint_list.append(opt_constraint)
                                
                        #         elif len(constraint) == 3:
                        #             opt_constraint = f"freeze angle {constraint[0] + 1} {constraint[1] + 1} {constraint[2] + 1}"
                        #             opt_constraint_list.append(opt_constraint)
                            
                        #         else:
                        #             opt_constraint = f"freeze dihedral {constraint[0] + 1} {constraint[1] + 1} {constraint[2] + 1} {constraint[3] + 1}"
                        #             opt_constraint_list.append(opt_constraint)
                        #     opt_drv.constraints = opt_constraint_list
                        #     self.drivers['es'][3].ostream.mute()
                        #     self.drivers['es'][0].ostream.mute()
                        #     opt_results = self._opt_compute_mpi_safe(opt_drv, molecule, current_basis, self.drivers['es'][3], self.drivers['es'][0], rsp_results)
                        #     excitated_roots = [root for root in self.roots_to_follow if root != 0]
                        #     self.drivers['es'][1].state_deriv_index = excitated_roots
                        #     optimized_molecule = Molecule.from_xyz_string(opt_results['final_geometry'])
                        #     optimized_molecule.set_charge(molecule.get_charge())
                        #     optimized_molecule.set_multiplicity(molecule.get_multiplicity())
                        #     print('Optimized Molecule', optimized_molecule.get_xyz_string(), '\n\n', molecule.get_xyz_string()) 
                        #     current_basis = MolecularBasis.read(optimized_molecule, current_basis.get_main_basis_label())       
                        
                        imp_int_coord = {'bonds': [], 'angles': [], 'dihedrals': [], 'impropers': []}
                        for element in imp_coord_constraint:
                        
                            if element in self.root_z_matrix[self.current_state]['bonds']:
                                imp_int_coord['bonds'].append(element)
                            elif element in self.root_z_matrix[self.current_state]['angles']:
                                imp_int_coord['angles'].append(element)
                            elif element in self.root_z_matrix[self.current_state]['dihedrals']:
                                imp_int_coord['dihedrals'].append(element)
                            elif element in self.root_z_matrix[self.current_state]['impropers']:
                                imp_int_coord['impropers'].append(element)
                        corr_new_dp_distrib = True
                    state_specific_molecules.append((optimized_molecule, current_basis, [state_to_optim], imp_int_coord))
                   
                    print('New optimized molecule \n', optimized_molecule.get_xyz_string())

                self.add_point(state_specific_molecules, self.non_core_symmetry_groups)
                self.last_point_added = self.point_checker - 1
                self.point_checker = 0
            else:
                    
                current_basis = None
                if self.current_state == 0:
                    current_basis = MolecularBasis.read(molecule, self.basis_set_label['gs'])
                else:
                    current_basis = MolecularBasis.read(molecule, self.basis_set_label['es'])
                state_specific_molecules.append((molecule, current_basis, addition_of_state_specific_points, {'bonds': [], 'angles': [], 'dihedrals': [], 'impropers':[]}))
                    
                
                self.add_point(state_specific_molecules, self.non_core_symmetry_groups)
                self.last_point_added = self.point_checker - 1
                self.point_checker = 0

    
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
                        else state.phase_signature.tolist()
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
    
    
    def _register_core_datapoint_runtime(self, root, impes_coordinate, new_label, energy_value):
        self.qm_data_point_dict[root].append(impes_coordinate)
        self.sorted_state_spec_im_labels[root].append(new_label)
        self.qm_energies_dict[root].append(energy_value)

        core_key = impes_coordinate.point_label
        self.qm_symmetry_datapoint_dict[root][core_key] = [impes_coordinate]

        drv = self.impes_drivers[root]
        drv.qm_data_points = self.qm_data_point_dict[root]
        drv.qm_symmetry_data_points = self.qm_symmetry_datapoint_dict[root]
        drv.labels = self.sorted_state_spec_im_labels[root]
        self._refresh_interpolation_driver_caches(root)


    def _refresh_rotor_cluster_runtime_for_root(self, root):
        self.qm_rotor_cluster_banks[root] = self._load_rotor_cluster_bank_for_root(root)
        drv = self.impes_drivers[root]
        drv.qm_rotor_cluster_banks = self.qm_rotor_cluster_banks[root]

        if drv.qm_rotor_cluster_banks:
            first_family = next(iter(drv.qm_rotor_cluster_banks.values()))
            drv.rotor_cluster_information = first_family.get("cluster_info")
        else:
            drv.rotor_cluster_information = None

        if self.sampling_enabled:
            samp_drv = self.sampling_impes_drivers[root]
            samp_drv.qm_rotor_cluster_banks = self.qm_rotor_cluster_banks[root]
            
            if samp_drv.qm_rotor_cluster_banks:
                first_family = next(iter(samp_drv.qm_rotor_cluster_banks.values()))
                samp_drv.rotor_cluster_information = first_family.get("cluster_info")
            else:
                samp_drv.rotor_cluster_information = None

        self._refresh_interpolation_driver_caches(root)
    
    def _finalize_mapping_masks_and_eq_mode(self, impes_coordinate, masks):
        impes_coordinate.mapping_masks = masks

        if not impes_coordinate.use_eq_bond_length:
            return

        mode = str(getattr(impes_coordinate, 'eq_bond_symmetry_mode', 'masked_exact')).strip().lower()
        if mode == 'symmetrized':
            impes_coordinate.symmetrize_eq_bond_lengths_from_masks(masks)
            impes_coordinate.transform_gradient_and_hessian()

    def add_point_rotor(self, molecule_specific_information, symmetry_information={}):
        """ Adds a new point to the database.

            :param molecule:
                the molecule.
            :param imforcefielddatafile:
                Datafile containing the information of the IM forcefield.

        """

        if len(self.drivers) == 0:
            raise ValueError("No energy driver defined.")
        
        self._emit_test_hook("add_point_start", {
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
                                phase_signature=None,
                            )
                        )

                else:
                    joint_phases = []
                    for rotor_id in rotor_ids:
                        joint_phases.append(phase_library[rotor_id])

                    state_id = 1

                    for phase_tuple in itertools.product(*joint_phases):
                        if all(abs(float(phase)) < 1.0e-12 for phase in phase_tuple):
                            continue
                        angle_assignment = {}
                        dihedrals_to_rotate = []
                        for rotor_id, phase in zip(rotor_ids, phase_tuple):
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
                                phase_signature=None,
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

        comm = self._mpi_collective_comm()
        rank = comm.Get_rank()
        root = mpi_master()

        root_worker_service_mode = bool(
            self._mpi_is_active()
            and self.mpi_root_worker_mode
            and self._mpi_worker_service_running
        )

        if root_worker_service_mode and rank != root:
            return

        # use the first incomming structure to determine the cluster grouping

        for entries in molecule_specific_information:
            molecule = entries[0]
            basis = entries[1]
            states = entries[2]
            symmetry_point = False

            if 0 in entries[2]:
                adjusted_molecule['gs'].append((entries[0], entries[1], 1, None, [0], symmetry_point, entries[3])) 
            if any(x > 0 for x in entries[2]): 
                states = [state for state in entries[2] if state > 0]
                adjusted_molecule['es'].append((entries[0], entries[1], 1, None, states, symmetry_point, entries[3])) 

        eq_bond_length = []
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
                outside_constraints = mol_basis[6]
                if not mol_basis[5]:
                    label_counter = 0
                if isinstance(drivers[0], LinearResponseEigenSolver) or isinstance(drivers[0], TdaEigenSolver):
                    root_to_follow_calc = mol_basis[4]           
                    drivers[1].state_deriv_index = root_to_follow_calc 
                    # drivers[2].roots_to_follow = root_to_follow_calc
                
                symmetry_mapping_groups = [item for item in range(len(mol_basis[0].get_labels()))]
                symmetry_exclusion_groups = [item for element in symmetry_information[state_key][1] for item in element]
   
                # scf_results_mpi = self._compute_energy_mpi_safe(
                #     drivers[0],
                #     mol_basis[0],
                #     mol_basis[1],
                #     collective=False,
                #     phase_name='Energy calculation',
                # )

                # energies, scf_results, rsp_results = scf_results_mpi[0], scf_results_mpi[1], scf_results_mpi[2]
                
                # if isinstance(drivers[0], LinearResponseEigenSolver) or isinstance(drivers[0], TdaEigenSolver):
                #     energies = energies[mol_basis[4]]
                # gradient_results_mpi = self._compute_gradient_mpi_safe(
                #     drivers[1],
                #     mol_basis[0],
                #     mol_basis[1],
                #     scf_results,
                #     rsp_results,
                #     collective=False,
                #     phase_name='Gradient calculation',
                # )
                # gradients = gradient_results_mpi
                
                # hessians_result_mpi = self._compute_hessian_mpi_safe(
                #     drivers[2],
                #     mol_basis[0],
                #     mol_basis[1],
                #     collective=False,
                #     phase_name='Hessian calculation',
                # )
                # hessians = hessians_result_mpi#
                energies, gradients, hessians = self._aux_eval_egh(
                    molecule=mol_basis[0],
                    basis_label=mol_basis[1].get_main_basis_label(),
                    driver_key=state_key,
                    state_mask=mol_basis[4],
                )

                masses = mol_basis[0].get_masses().copy()
                masses_cart = np.repeat(masses, 3)
                inv_sqrt_masses = 1.0 / np.sqrt(masses_cart)

                state_spec_dp = {}
                if rank == root:

                    for number in range(len(energies)):
                        target_root = mol_basis[4][number]
                        target_file = self.interpolation_settings[target_root]['imforcefield_file']
                        
                        z_matrix = self.root_z_matrix[target_root]
                        interpolation_driver = InterpolationDriver(z_matrix)
                        interpolation_driver.update_settings(self.interpolation_settings[target_root])
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

                        if len(eq_bond_length) == 0:
                            for idx, element in enumerate(z_matrix['bonds']):
                                if 1 == 1 or len(element) == 2 and self.use_minimized_structures[0]:
                                    eq_bond_length.append(mol_basis[0].get_distance([element[0] + 1, element[1] + 1], 'bohr'))
                                elif len(element) == 2:
                                    eq_bond_length.append(0.0)
                        
                        impes_coordinate = InterpolationDatapoint(z_matrix)
                        impes_coordinate.eq_bond_lengths = self.qm_data_point_dict[mol_basis[4][number]][0].eq_bond_lengths
                        impes_coordinate.update_settings(self.interpolation_settings[target_root])
                        impes_coordinate.cartesian_coordinates = mol_basis[0].get_coordinates_in_bohr()
                        impes_coordinate.imp_int_coordinates = {'bonds': [], 'angles': [], 'dihedrals': [], 'impropers': []}
                        impes_coordinate.inv_sqrt_masses = inv_sqrt_masses
                        impes_coordinate.energy = energies[number]
                        impes_coordinate.gradient =  mw_grad_vec.reshape(grad.shape)
                        impes_coordinate.hessian = mw_hess_mat.reshape(hess.shape)
                        impes_coordinate.transform_gradient_and_hessian()


                        # impes_coordinate.point_label = label
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

                                target_coordinates = self.calculate_translation_coordinates(rot_mol.get_coordinates_in_bohr())
                                reference_coordinates = self.calculate_translation_coordinates(mol_basis[0].get_coordinates_in_bohr())
                                
                                target_coordinates_core = target_coordinates.copy()
                                reference_coordinates_core = reference_coordinates.copy()
                                    
                                target_coordinates_core = np.delete(target_coordinates, symmetry_exclusion_groups, axis=0)
                                reference_coordinates_core = np.delete(reference_coordinates, symmetry_exclusion_groups, axis=0)
                                rotation_matrix = geometric.rotate.get_rot(target_coordinates_core,
                                                                        reference_coordinates_core)
                                # Rotate the data point
                                rotated_coordinates = np.dot(rotation_matrix, target_coordinates.T).T
                                rotated_coordinates = np.ascontiguousarray(rotated_coordinates)

                                
                                mapping_dict_12 = self.perform_symmetry_assignment(symmetry_mapping_groups, symmetry_exclusion_groups, 
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

                        impes_coordinate.write_hdf5(self.interpolation_settings[mol_basis[4][number]]['imforcefield_file'], label)
                        if self.sampling_enabled:
                            self._write_sampling_point_from_geometry(
                                root=mol_basis[4][number],
                                molecule=mol_basis[0],
                                label=label,
                                template_point=impes_coordinate
                            )
                        
                        self._emit_test_hook("datapoint_written", {
                            "root": int(target_root),
                            "label": str(label),
                            "energy_hartree": float(impes_coordinate.energy),
                            "confidence_radius": float(impes_coordinate.confidence_radius),
                            "bank_role": str(getattr(impes_coordinate, "bank_role", "core")),
                        })
                        
                        self.simulation.saveCheckpoint('checkpoint')    

                        self._register_core_datapoint_runtime(
                            target_root,
                            impes_coordinate,
                            label,
                            energies[number],
                        )
                        self.density_around_data_point[0][target_root] += 1
                        # impes_coordinate.write_hdf5(self.interpolation_settings[target_root]['imforcefield_file'], label)
                        print('Family label part', impes_coordinate.family_label)
                        print('Point label in core part', impes_coordinate.point_label)

                        state_spec_dp[target_root] = impes_coordinate
                        
                state_spec_payload = None
                if rank == root:
                    state_spec_payload = {
                        int(state): {
                            "point_label": dp.point_label,
                            "family_label": dp.family_label,
                            "eq_bond_lengths": np.array(dp.eq_bond_lengths, copy=True),
                            "internal_coordinates_values": np.array(dp.internal_coordinates_values, copy=True),
                            "internal_hessian": np.array(dp.internal_hessian, copy=True),
                            "imp_int_coordinates": copy.deepcopy(
                                getattr(
                                    dp,
                                    "imp_int_coordinates",
                                    {"bonds": [], "angles": [], "dihedrals": [], "impropers": []},
                                )
                            ),
                        }
                        for state, dp in state_spec_dp.items()
                    }

                if root_worker_service_mode:
                    # Workers are in _run_qmmm_worker_service, so do not use raw bcast here.
                    state_spec_payload = state_spec_payload if rank == root else {}
                else:
                    state_spec_payload = comm.bcast(
                        state_spec_payload if rank == root else None,
                        root,
                    )

                for state, cur_dp_payload in state_spec_payload.items():
                    cur_dp = type("CorePointMeta", (), {})()
                    cur_dp.point_label = cur_dp_payload["point_label"]
                    cur_dp.family_label = cur_dp_payload["family_label"]
                    cur_dp.eq_bond_lengths = cur_dp_payload["eq_bond_lengths"]
                    cur_dp.internal_coordinates_values = cur_dp_payload["internal_coordinates_values"]
                    cur_dp.internal_hessian = cur_dp_payload["internal_hessian"]
                    cur_dp.imp_int_coordinates = copy.deepcopy(
                        cur_dp_payload.get(
                            "imp_int_coordinates",
                            {"bonds": [], "angles": [], "dihedrals": [], "impropers": []},
                        )
                    )

                    current_constraints = cur_dp.imp_int_coordinates


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

                    clusters = build_rotor_clusters(self.symmetry_rotors, coupling_map, threshold=0.04)

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
                        cur_molecule = Molecule.from_xyz_string(molecule.get_xyz_string())
                        cur_molecule.set_charge(molecule.get_charge())
                        cur_molecule.set_multiplicity(molecule.get_multiplicity())
                        current_angles_setting = []
                        constraints = []
                        for key ,entries in current_constraints.items():
                            for entry in entries:
                                constraints.append(entry)
                        # constraints = outside_constraints.copy()
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

                        optimized_molecule, opt_results = self._aux_optimize(
                            state_to_optim=state,
                            molecule=cur_molecule,
                            basis_label=current_basis.get_main_basis_label(),
                            constraints=constraints,
                        )
                        cur_molecule = optimized_molecule
                        
                      
                        # if isinstance(drivers[0], ScfRestrictedDriver):
                        #     _, scf_results, _ = self._compute_energy_mpi_safe(
                        #         drivers[0],
                        #         cur_molecule,
                        #         current_basis,
                        #         collective=True,
                        #         phase_name='Energy calcualtion',
                        #     )
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
                        #     optimized_molecule = opt_results['final_molecule']
                        #     cur_molecule = optimized_molecule
                            
                        # elif isinstance(drivers[0], XtbDriver):
                            
                        #     opt_results_mpi = self._run_optimization_mpi_safe(
                        #         drivers[0],
                        #         cur_molecule,
                        #         constraints=constraints,
                        #         index_offset=1,
                        #         source_molecule=molecule,
                        #         collective=True,
                        #         phase_name='optimization tag',
                        #     )
                        
                        #     opt_results = opt_results_mpi[1]
                        #     optimized_molecule = opt_results['final_molecule']
                        #     print(optimized_molecule.get_xyz_string())
                        #     cur_molecule = optimized_molecule
                        current_basis = MolecularBasis.read(cur_molecule, basis.get_main_basis_label())

                        energies, gradients, hessians = self._aux_eval_egh(
                            molecule=cur_molecule,
                            basis_label=current_basis.get_main_basis_label(),
                            driver_key=state_key,
                            state_mask=[state],
                        )

                        # scf_results_mpi = self._compute_energy_mpi_safe(
                        #     drivers[0],
                        #     cur_molecule,
                        #     current_basis,
                        #     collective=True,
                        #     phase_name='Energy calcualtion',
                        # )

                        # energies, scf_results, rsp_results = scf_results_mpi[0], scf_results_mpi[1], scf_results_mpi[2]

                        # print('energies', energies)

                        # if isinstance(drivers[0], LinearResponseEigenSolver) or isinstance(drivers[0], TdaEigenSolver):
                        #     energies = energies[mol_basis[4]]
                        
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
                        if rank == root:
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

                                    target_coordinates = self.calculate_translation_coordinates(rot_mol.get_coordinates_in_bohr())
                                    reference_coordinates = self.calculate_translation_coordinates(cur_molecule.get_coordinates_in_bohr())
                                    
                                    target_coordinates_core = target_coordinates.copy()
                                    reference_coordinates_core = reference_coordinates.copy()
                                        
                                    target_coordinates_core = np.delete(target_coordinates, symmetry_exclusion_groups, axis=0)
                                    reference_coordinates_core = np.delete(reference_coordinates, symmetry_exclusion_groups, axis=0)
                                    rotation_matrix = geometric.rotate.get_rot(target_coordinates_core,
                                                                            reference_coordinates_core)
                                    # Rotate the data point
                                    rotated_coordinates = np.dot(rotation_matrix, target_coordinates.T).T
                                    rotated_coordinates = np.ascontiguousarray(rotated_coordinates)

                                    
                                    mapping_dict_12 = self.perform_symmetry_assignment(symmetry_mapping_groups, symmetry_exclusion_groups, 
                                                                                        cur_molecule.get_coordinates_in_bohr()[symmetry_exclusion_groups], rotated_coordinates[symmetry_exclusion_groups])
                                    
                                    
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
                                            print(reorded_int_coords, z_mat_index, z_matrix_dict)
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
                                impes_coordinate.update_settings(self.interpolation_settings[target_root])
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
                                impes_coordinate.confidence_radius = trust_radius
                                
                                self._finalize_mapping_masks_and_eq_mode(impes_coordinate, masks)

                                impes_coordinate.write_hdf5(self.interpolation_settings[target_root]['imforcefield_file'], label)
                                if self.sampling_enabled:
                                        self._write_sampling_point_from_geometry(
                                            root=target_root,
                                            molecule=cur_molecule,
                                            label=label,
                                            template_point=impes_coordinate,
                                        )
                                self._emit_test_hook("datapoint_written", {
                                    "root": int(target_root),
                                    "label": str(label),
                                    "energy_hartree": float(impes_coordinate.energy),
                                    "confidence_radius": float(impes_coordinate.confidence_radius),
                                    "bank_role": str(getattr(impes_coordinate, "bank_role", "core")),
                                })
                                
                                print('Symmetry rotor label', label)

                                point_index["cluster_state_labels"].setdefault(int(job["cluster_id"]), {})
                                point_index["cluster_state_labels"][int(job["cluster_id"])][int(job["state_id"])] = label

                    if rank == root:
                        self._write_cluster_registry_for_family(self.interpolation_settings[target_root]['imforcefield_file'], target_root, cur_dp.family_label, rotor_to_cluster_inf[state_key], final_clusters_angle_lib, point_index)
                        self._refresh_rotor_cluster_runtime_for_root(target_root)

                    if not root_worker_service_mode:
                        comm.barrier()
                    label_counter += 1
                    if any(root > 0 for root in self.roots_to_follow): 
                        if isinstance(drivers[0], LinearResponseEigenSolver) or isinstance(drivers[0], TdaEigenSolver):           
                            drivers[1].state_deriv_index = org_roots

        self._emit_test_hook("add_point_end", {
            "n_datapoints_per_root": {
                int(root): int(len(self.qm_data_point_dict[root]))
                for root in self.roots_to_follow
            }
        })                

    def add_point(self, state_specific_molecules, symmetry_information):
        """ Adds a new point to the database.

            :param molecule:
                the molecule.
            :param label:
                the label for the new point to be added.    
            :param energy:
                the energy of the previous QM calcualtion.
            :param basis:
                the basis set (if required).
            :scf_results:
                the scf_results of previous QM calculation (if required).
        """
        
        if self.cluster_run:
            return self.add_point_rotor(state_specific_molecules, symmetry_information)

        self._emit_test_hook("add_point_start", {
            "n_molecules": int(len(state_specific_molecules)),
            "cluster_run": bool(self.cluster_run),
        })
            
        adjusted_molecule = {'gs': [], 'es': []}
        symmetry_mapping_groups = []
        symmetry_exclusion_groups = []
        category_label = None

        comm = self._mpi_collective_comm()
        rank = comm.Get_rank()
        root = mpi_master()

        root_worker_service_mode = bool(
            self._mpi_is_active()
            and self.mpi_root_worker_mode
            and self._mpi_worker_service_running
        )

        if root_worker_service_mode and rank != root:
            return

        for entries in state_specific_molecules:
            symmetry_point = False
            if 0 in entries[2] and len(symmetry_information['gs']) != 0 and len(symmetry_information['gs'][2]) != 0:
                symmetry_mapping_groups = [item for item in range(len(entries[0].get_labels()))]
                symmetry_exclusion_groups = [item for element in symmetry_information['gs'][1] for item in element]

                sym_dihedrals, periodicites, _, _ = self.adjust_symmetry_dihedrals(symmetry_information['gs'][1], symmetry_information['gs'][5], self.root_z_matrix[0])
                
                # Generate all combinations
                keys = list(sym_dihedrals.keys())
                values = list(sym_dihedrals.values())
                combinations = list(itertools.product(*values))
                # Convert to list of dictionaries
                molecule_configs = [dict(zip(keys, combo)) for combo in combinations] # if all(element == combo[0] for element in combo)

                for i, molecule_config in enumerate(molecule_configs):
                    cur_molecule = Molecule.from_xyz_string(entries[0].get_xyz_string())
                    cur_molecule.set_charge(entries[0].get_charge())
                    cur_molecule.set_multiplicity(entries[0].get_multiplicity())
                    dihedral_to_change = []
                    optim_dih_0based = []
                    for dihedral, angle in molecule_config.items():
             
                        opt_dihedral_val = cur_molecule.get_dihedral([dihedral[0] + 1, dihedral[1] + 1, dihedral[2] + 1, dihedral[3] + 1], 'radian')
                        cur_molecule.set_dihedral([dihedral[0] + 1, dihedral[1] + 1, dihedral[2] + 1, dihedral[3] + 1], opt_dihedral_val + angle, 'radian')
                        dihedral_to_change.append([dihedral[0] + 1, dihedral[1] + 1, dihedral[2] + 1, dihedral[3] + 1])
                        optim_dih_0based.append(dihedral)
                    current_basis = MolecularBasis.read(cur_molecule, entries[1].get_main_basis_label())
                    
                    if i > 0:
                        symmetry_point = True
                    
                    current_basis = MolecularBasis.read(cur_molecule, entries[1].get_main_basis_label())
                    for key ,add_const in entries[3].items():
                            for ac in add_const:
                                optim_dih_0based.append(ac)
                    
                    optimized_molecule, opt_results = self._aux_optimize(
                            state_to_optim=0,
                            molecule=cur_molecule,
                            basis_label=current_basis.get_main_basis_label(),
                            constraints=optim_dih_0based,
                        )
                    
                    print(i, opt_results['opt_energies'][-1], optimized_molecule.get_xyz_string())
                    adjusted_molecule['gs'].append((optimized_molecule, current_basis, periodicites[dihedral], dihedral_to_change, entries[2], symmetry_point, entries[3]))
            
            elif any(x > 0 for x in entries[2]) and len(symmetry_information['es']) != 0 and len(symmetry_information['es'][2]) != 0:
                symmetry_mapping_groups = [item for item in range(len(entries[0].get_labels()))]
                symmetry_exclusion_groups = [item for element in symmetry_information['es'][1] for item in element]
                state_idx_b_0 = 0
                for ent in entries[2]:
                    if ent > 0:
                        state_idx_b_0 = ent
                        break
                sym_dihedrals, periodicites, _, _ = self.adjust_symmetry_dihedrals(symmetry_information['es'][1], symmetry_information['es'][5], self.root_z_matrix[state_idx_b_0])
                
                # Generate all combinations
                keys = list(sym_dihedrals.keys())
                values = list(sym_dihedrals.values())
                combinations = list(itertools.product(*values))
                # Convert to list of dictionaries
                molecule_configs = [dict(zip(keys, combo)) for combo in combinations] # if all(element == combo[0] for element in combo)
                for i, molecule_config in enumerate(molecule_configs):
                    cur_molecule = Molecule.from_xyz_string(entries[0].get_xyz_string())
                    cur_molecule.set_charge(entries[0].get_charge())
                    cur_molecule.set_multiplicity(entries[0].get_multiplicity())
                    dihedral_to_change = []
                    for dihedral, angle in molecule_config.items():
                        cur_molecule.set_dihedral([dihedral[0] + 1, dihedral[1] + 1, dihedral[2] + 1, dihedral[3] + 1], angle, 'radian')
                        dihedral_to_change.append([dihedral[0] + 1, dihedral[1] + 1, dihedral[2] + 1, dihedral[3] + 1])
                    if i > 0:
                        symmetry_point = True   
                    current_basis = MolecularBasis.read(cur_molecule, entries[1].get_main_basis_label())
                    adjusted_molecule['es'].append((cur_molecule, current_basis, periodicites[dihedral],  dihedral_to_change, entries[2], symmetry_point, entries[3]))
            else:
                if 0 in entries[2]:
                    adjusted_molecule['gs'].append((entries[0], entries[1], 1, None, [0], symmetry_point, entries[3])) 
                if any(x > 0 for x in entries[2]): 
                    states = [state for state in entries[2] if state > 0]
                    adjusted_molecule['es'].append((entries[0], entries[1], 1, None, states, symmetry_point, entries[3])) 

        
        for state_key, entries in adjusted_molecule.items():
            if len(entries) == 0:
                continue
            print(state_key, entries)
            drivers = self.drivers[state_key]
            
            org_roots = self.roots_to_follow[0]
            if any(root > 0 for root in self.roots_to_follow): 
                if isinstance(drivers[0], LinearResponseEigenSolver) or isinstance(drivers[0], TdaEigenSolver):           
                    org_roots = drivers[1].state_deriv_index  

            label_counter = 0
            for mol_basis in entries:
                if any(root > 0 for root in self.roots_to_follow):
                    if isinstance(drivers[0], LinearResponseEigenSolver) or isinstance(drivers, TdaEigenSolver):
                        root_to_follow_calc = mol_basis[4]           
                        drivers[1].state_deriv_index = root_to_follow_calc 
                
                # energies, scf_results, rsp_results = self._compute_energy_mpi_safe(drivers[0],
                #                                                                    mol_basis[0], 
                #                                                                    mol_basis[1],
                #                                                                    collective=False,
                #                                                                    phase_name='Energy calculation')
                
                # if isinstance(drivers[0], LinearResponseEigenSolver) or isinstance(drivers[0], TdaEigenSolver):
                #     energies = energies[mol_basis[4]]

                # gradients = self._compute_gradient_mpi_safe(drivers[1], 
                #                                             mol_basis[0], 
                #                                             mol_basis[1], 
                #                                             scf_results, 
                #                                             rsp_results,
                #                                             collective=False,
                #                                             phase_name='Gradient calculation')
                
                # hessians = self._compute_hessian_mpi_safe(drivers[2], 
                #                                           mol_basis[0], 
                #                                           mol_basis[1],
                #                                           collective=False,
                #                                           phase_name='Energy calculation')

                energies, gradients, hessians = self._aux_eval_egh(
                    molecule=mol_basis[0],
                    basis_label=mol_basis[1].get_main_basis_label(),
                    driver_key=state_key,
                    state_mask=mol_basis[4],
                )

                masses = mol_basis[0].get_masses().copy()
                masses_cart = np.repeat(masses, 3)
                inv_sqrt_masses = 1.0 / np.sqrt(masses_cart)
                for number in range(len(energies)):

                    label = f"point_{len(self.qm_data_point_dict[mol_basis[4][number]])}"
                    if not mol_basis[5]:
                        label = f"point_{len(self.qm_data_point_dict[mol_basis[4][number]]) + 1}"
                    new_label = label

                    grad_vec = gradients[number].reshape(-1)         # (3N,)
                    hess_mat = hessians[number].reshape(grad_vec.size, grad_vec.size)
                    mw_grad_vec = inv_sqrt_masses * grad_vec
                    mw_hess_mat = (inv_sqrt_masses[:, None] * hess_mat) * inv_sqrt_masses[None, :]
                    

                    impes_coordinate = InterpolationDatapoint(self.root_z_matrix[mol_basis[4][number]])
                    impes_coordinate.update_settings(self.interpolation_settings[mol_basis[4][number]])
                    impes_coordinate.eq_bond_lengths = self.qm_data_point_dict[mol_basis[4][number]][0].eq_bond_lengths
                    impes_coordinate.cartesian_coordinates = mol_basis[0].get_coordinates_in_bohr()
                    impes_coordinate.imp_int_coordinates = mol_basis[6]
                    impes_coordinate.inv_sqrt_masses = inv_sqrt_masses
                    impes_coordinate.energy = energies[number]
                    impes_coordinate.gradient = mw_grad_vec.reshape(gradients[number].shape)
                    impes_coordinate.hessian = mw_hess_mat.reshape(hessians[number].shape)
                    impes_coordinate.transform_gradient_and_hessian()

                    impes_coordinate.confidence_radius = float(self.use_opt_confidence_radius[2])

                    if mol_basis[5]:
                            new_label = f'{label}_symmetry_{label_counter}'

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
                        print('rot comb', rotation_combinations)
                        
                        org_mask = [i for i in range(len(impes_coordinate.z_matrix))]
                        masks = [org_mask]
                        for combo in rotation_combinations:
                            if all(r == 0 for r in combo):
                                continue  # skip the base geometry if desired
                            
                    
                            rot_mol = Molecule(self.molecule.get_labels(), mol_basis[0].get_coordinates_in_bohr(), 'bohr')
                            for angle, dihedral in zip(combo, dihedrals):
                                
                                rot_mol.set_dihedral(dihedral, mol_basis[0].get_dihedral(dihedral, 'radian') - angle, 'radian')
                            target_coordinates = self.calculate_translation_coordinates(rot_mol.get_coordinates_in_bohr())
                            reference_coordinates = self.calculate_translation_coordinates(mol_basis[0].get_coordinates_in_bohr())
                            
                            target_coordinates_core = target_coordinates.copy()
                            reference_coordinates_core = reference_coordinates.copy()
                                
                            target_coordinates_core = np.delete(target_coordinates, symmetry_exclusion_groups, axis=0)
                            reference_coordinates_core = np.delete(reference_coordinates, symmetry_exclusion_groups, axis=0)
                            rotation_matrix = geometric.rotate.get_rot(target_coordinates_core,
                                                                    reference_coordinates_core)
                            # Rotate the data point
                            rotated_coordinates = np.dot(rotation_matrix, target_coordinates.T).T
                            rotated_coordinates = np.ascontiguousarray(rotated_coordinates)
    
                            mapping_dict_12 = self.perform_symmetry_assignment(symmetry_mapping_groups, symmetry_exclusion_groups, 
                                                                                mol_basis[0].get_coordinates_in_bohr()[symmetry_exclusion_groups], rotated_coordinates[symmetry_exclusion_groups])
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
                        self._finalize_mapping_masks_and_eq_mode(impes_coordinate, masks)
                    else:
                        org_mask = [i for i in range(len(impes_coordinate.z_matrix))]
                        masks = [org_mask]
                        self._finalize_mapping_masks_and_eq_mode(impes_coordinate, masks)

                    print('Label counter', label_counter, len(self.qm_data_point_dict[mol_basis[4][number]]), new_label)
                    if not mol_basis[5]:
                        category_label = label
                        if impes_coordinate.use_inverse_bond_length:
                            category_label += "_rinv"
                        else:
                            category_label += "_r"
                        if impes_coordinate.use_cosine_dihedral:
                            category_label += "_cosine"
                        else:
                            category_label += "_dihedral"
                        
                        print('dict and key', self.impes_drivers[mol_basis[4][number]].qm_symmetry_data_points, category_label)
                        self.qm_data_point_dict[mol_basis[4][number]].append(impes_coordinate)
                        self.sorted_state_spec_im_labels[mol_basis[4][number]].append(new_label)
                        self.qm_symmetry_datapoint_dict[mol_basis[4][number]][category_label] = [impes_coordinate]
                        self.impes_drivers[mol_basis[4][number]].qm_symmetry_data_points[category_label] = [impes_coordinate]
                        impes_coordinate.point_label = category_label
                        impes_coordinate.family_label = category_label
                        impes_coordinate.bank_role = "core"
                        
                        self.qm_energies_dict[mol_basis[4][number]].append(energies[number])
                        
                        self.impes_drivers[mol_basis[4][number]].qm_data_points = self.qm_data_point_dict[mol_basis[4][number]]
                        self.impes_drivers[mol_basis[4][number]].labels = self.sorted_state_spec_im_labels[mol_basis[4][number]]
                        self._refresh_interpolation_driver_caches(mol_basis[4][number])
                        print('I am writing the file', category_label, self.interpolation_settings[mol_basis[4][number]]['imforcefield_file'])
                    
                        self.density_around_data_point[0][mol_basis[4][number]] += 1

                    else:
                        
                        impes_coordinate.bank_role = "symmetry"
                        self.qm_symmetry_datapoint_dict[mol_basis[4][number]][category_label].append(impes_coordinate)
                        self.impes_drivers[mol_basis[4][number]].qm_symmetry_data_points[category_label].append(impes_coordinate)
                        self._refresh_interpolation_driver_caches(mol_basis[4][number])

                    
                    impes_coordinate.confidence_radius = self.use_opt_confidence_radius[2]

                    print('data list', len(self.allowed_molecules[mol_basis[4][number]]['molecules']))

                    if self._mpi_is_root() or (not self._mpi_is_active()):
                        impes_coordinate.write_hdf5(
                            self.interpolation_settings[mol_basis[4][number]]['imforcefield_file'],
                            new_label
                        )

                        self._emit_test_hook("datapoint_written", {
                            "root": int(mol_basis[4][number]),
                            "label": str(new_label),
                            "energy_hartree": float(impes_coordinate.energy),
                            "confidence_radius": float(impes_coordinate.confidence_radius),
                            "bank_role": str(getattr(impes_coordinate, "bank_role", "core")),
                        })


                    # Call on all active ranks in full-SPMD mode so in-memory sampling caches stay consistent.
                    if self.sampling_enabled:
                        self._write_sampling_point_from_geometry(
                            root=mol_basis[4][number],
                            molecule=mol_basis[0],
                            label=new_label,
                            template_point=impes_coordinate
                        )

                    # Only synchronize in full-SPMD mode; NEVER do this in root-worker service mode.
                    if self._mpi_is_active() and (not root_worker_service_mode):
                        comm.barrier()
                    
                    # if self._mpi_is_root() or (not self._mpi_is_active()):
                    #     impes_coordinate.write_hdf5(self.interpolation_settings[mol_basis[4][number]]['imforcefield_file'], new_label)
                    #     if self.sampling_enabled:
                    #         self._write_sampling_point_from_geometry(
                    #             root=mol_basis[4][number],
                    #             molecule=mol_basis[0],
                    #             label=new_label,
                    #             template_point=impes_coordinate
                    #         )
                    if self._mpi_is_root() or (not self._mpi_is_active()):
                        for root_idx in mol_basis[4]:
                            self.write_qm_energy_determined_points_h5(self.allowed_molecules[root_idx]['molecules'][self.last_added: ],
                                                                self.allowed_molecules[root_idx]['qm_energies'][self.last_added:],
                                                                self.allowed_molecules[root_idx]['qm_gradients'][self.last_added:],
                                                                self.allowed_molecules[root_idx]['im_energies'][self.last_added:],
                                                                root_idx)
                
                self.last_added = len(self.allowed_molecules[root]['molecules'])
                label_counter += 1
            if any(root > 0 for root in self.roots_to_follow):
                                   
                if isinstance(drivers[0], LinearResponseEigenSolver) or isinstance(drivers[0], TdaEigenSolver):           
                    
                    drivers[1].state_deriv_index = org_roots

            if self.current_state in self.opt_mols_org_mol_swap:
                self._interp_compute_root_local_serial(self.impes_drivers[self.current_state], self.opt_mols_org_mol_swap[self.current_state])
                print(self.impes_drivers[self.current_state].weights, self.impes_drivers[self.current_state].get_energy(), self.opt_mols_org_mol_swap[self.current_state].get_xyz_string())

        
        if self._mpi_is_root() or (not self._mpi_is_active()):
            self.simulation.saveCheckpoint('checkpoint')
            self._emit_test_hook("add_point_end", {
                "n_datapoints_per_root": {
                    int(root): int(len(self.qm_data_point_dict[root]))
                    for root in self.roots_to_follow
                }
            })

        if self._mpi_is_active() and (not root_worker_service_mode):
            comm.barrier()
    
    def _list_rotor_cluster_families_from_registry(self, root):
        imff_file = self.interpolation_settings[root]["imforcefield_file"]
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


    def _load_rotor_cluster_bank_for_root(self, root):
        out = {}
        imff_file = self.interpolation_settings[root]["imforcefield_file"]
        families = self._list_rotor_cluster_families_from_registry(root)

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

            core_dp = InterpolationDatapoint(self.root_z_matrix[root])
            core_dp.update_settings(self.interpolation_settings[root])
            core_dp.read_hdf5(imff_file, point_index["core_label"])
            fam["core"] = core_dp

            for cid, state_map in point_index["cluster_state_labels"].items():
                for sid, label in state_map.items():
                    dp = InterpolationDatapoint(self.root_z_matrix[root])
                    dp.update_settings(self.interpolation_settings[root])
                    dp.read_hdf5(imff_file, label)
                    fam["clusters"][cid]["expected_states"][sid] = dp

            for cid, cbank in fam["clusters"].items():
                expected = cbank["expected_states"]
                if 0 in expected and expected[0] is None:
                    expected[0] = fam["core"]

            out[family] = fam
        return out
     
    def _write_sampling_point_from_geometry(self, root, molecule, label, template_point):
        sampling_qm, sampling_grad, sampling_hess = self.sampling_driver['gs']
        e, _, _ = self._compute_energy_mpi_safe(sampling_qm, molecule, basis=None)
        g = self._compute_gradient_mpi_safe(
            sampling_grad, molecule, basis=None, scf_results=None, rsp_results=None
        )
        h = self._compute_hessian_mpi_safe(sampling_hess, molecule, basis=None)

        masses = molecule.get_masses().copy()
        inv_sqrt = 1.0 / np.sqrt(np.repeat(masses, 3))
        grad_vec = g[0].reshape(-1)
        hess_mat = h[0].reshape(grad_vec.size, grad_vec.size)
        mw_grad = inv_sqrt * grad_vec
        mw_hess = (inv_sqrt[:, None] * hess_mat) * inv_sqrt[None, :]

        dp = InterpolationDatapoint(self.root_z_matrix[root])
        dp.update_settings(self.sampling_interpolation_settings[root])

        # Geometry + transformed tensors
        dp.cartesian_coordinates = np.array(template_point.cartesian_coordinates, copy=True)
        dp.eq_bond_lengths = None if template_point.eq_bond_lengths is None else np.array(template_point.eq_bond_lengths, copy=True)
        dp.mapping_masks = None if getattr(template_point, "mapping_masks", None) is None else np.array(template_point.mapping_masks, copy=True)
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

        # Carry rotor-cluster metadata into sampling datapoints.
        dp.point_label = label
        dp.family_label = getattr(template_point, "family_label", None)
        dp.bank_role = str(getattr(template_point, "bank_role", "core") or "core")

        cluster_id = getattr(template_point, "cluster_id", None)
        dp.cluster_id = None if cluster_id is None else int(cluster_id)

        dp.cluster_type = getattr(template_point, "cluster_type", None)

        rotor_ids = getattr(template_point, "cluster_rotor_ids", None)
        dp.cluster_rotor_ids = None if rotor_ids is None else tuple(int(x) for x in rotor_ids)

        cluster_state_id = getattr(template_point, "cluster_state_id", 0)
        dp.cluster_state_id = 0 if cluster_state_id is None else int(cluster_state_id)

        dihedrals_to_rotate = getattr(template_point, "dihedrals_to_rotate", None)
        if dihedrals_to_rotate is None:
            dp.dihedrals_to_rotate = None
        else:
            dp.dihedrals_to_rotate = tuple(
                tuple(int(x) for x in row) for row in dihedrals_to_rotate
            )

        phase_signature = getattr(template_point, "phase_signature", None)
        dp.phase_signature = None if phase_signature is None else np.asarray(phase_signature, dtype=np.float64)

        dp.transform_gradient_and_hessian()

        if self._mpi_is_root() or (not self._mpi_is_active()):
            dp.write_hdf5(self.sampling_interpolation_settings[root]['imforcefield_file'], label)

        # Refresh in-memory sampling interpolator
        self._reload_sampling_root_from_hdf5(root, inv_sqrt)




    def _reload_sampling_root_from_hdf5(self, root, inv_sqrt_masses):
        driver_object = self.sampling_impes_drivers[root]

        im_labels, _ = driver_object.read_labels()

        self.sampling_qm_data_point_dict[root] = []
        self.sampling_qm_symmetry_datapoint_dict[root] = {}

        old_label = None
        for label in im_labels:
            dp = InterpolationDatapoint(self.root_z_matrix[root])
            dp.update_settings(self.sampling_interpolation_settings[root])
            dp.read_hdf5(self.sampling_interpolation_settings[root]['imforcefield_file'], label)
            dp.inv_sqrt_masses = inv_sqrt_masses
            if dp.bank_role == "core" or "cluster" not in dp.point_label and "symmetry" not in dp.point_label:
                self.sampling_qm_data_point_dict[root].append(dp)
                old_label = dp.point_label
                self.sampling_qm_symmetry_datapoint_dict[root][old_label] = [dp]
            elif dp.bank_role == "symmetry":
                self.sampling_qm_symmetry_datapoint_dict[root][old_label].append(dp)
            

        driver_object.qm_symmetry_data_points = self.sampling_qm_symmetry_datapoint_dict[root]
        driver_object.qm_data_points = self.sampling_qm_data_point_dict[root]
        driver_object.labels = self.sorted_state_spec_im_labels[root]

        if len(self.sampling_qm_data_point_dict[root]) > 0:
            driver_object.impes_coordinate.eq_bond_lengths = self.sampling_qm_data_point_dict[root][0].eq_bond_lengths

        driver_object.mark_runtime_data_cache_dirty()
        driver_object.prepare_runtime_data_cache(force=True)

        use_mpi_preload = bool(self.sampling_interpolation_settings[root].get('use_mpi_preload', False))
        # In root-worker mode, avoid immediate preload collectives here; rebuild lazily
        # in synchronized compute path on the next step.
        force_rebuild_preload = not (self._mpi_is_active() and self.mpi_root_worker_mode)
        driver_object.set_mpi_preload_engine(
            comm=self._mpi_interp_comm,
            enabled=(use_mpi_preload and self.nodes > 1),
            force_rebuild=force_rebuild_preload,
        )

        driver_object._mpi_preload_enabled_config = bool(getattr(driver_object, 'mpi_preload_enabled', False))
        self._ensure_interp_engine_comm(driver_object)
        # In root-worker mode, keep preload disabled outside synchronized compute
        # phases to prevent root-only branches from emitting collective packets.
        if self._mpi_is_active() and self.mpi_root_worker_mode:
            driver_object.mpi_preload_enabled = False

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

    
    def determine_beysian_trust_radius(self, molecules, qm_energies, current_datapoints, interpolation_setting, sym_datapoints, sym_dict, z_matrix):
    
        trust_radii = []
        for dp in current_datapoints:   
            sum_sq_error = 0.0
            combined_datapoints = [dp]

            for i, mol in enumerate(molecules):
                _, distance, _ = self.calculate_distance_to_ref(mol.get_coordinates_in_bohr(), dp.cartesian_coordinates)

                interpolation_driver = InterpolationDriver(z_matrix)
                interpolation_driver.update_settings(interpolation_setting)
                interpolation_driver.symmetry_information = sym_dict
                interpolation_driver.qm_symmetry_data_points = sym_datapoints
                interpolation_driver.print = False
                interpolation_driver.qm_data_points = combined_datapoints
                
                interpolation_driver.compute(mol)
                new_im_energy = interpolation_driver.get_energy()
                diff = (new_im_energy - qm_energies[i]) * hartree_in_kcalpermol()
                sum_sq_error += (diff)**2 / (0.1**2 * distance**6)

            bey_trust_radius = (1/sum_sq_error)**(1/6)

            trust_radii.append(bey_trust_radius)
            
        return trust_radii


    def _build_dp_ref_distance_matrix(self, datapoints, molecules, symmetry_info):
        """
        D[i, j] = cartesian distance between datapoint i and reference molecule j
        using calculate_distance_to_ref (translation + rotational alignment + symmetry filtering).
        """
        n_dp = len(datapoints)
        n_ref = len(molecules)
        D_db = np.full((n_dp, n_dp), np.inf, dtype=np.float64)
        D_ref = np.full((n_ref, n_ref), np.inf, dtype=np.float64)
        D_dpref = np.full((n_dp, n_ref), np.inf, dtype=np.float64)


        for i, dp in enumerate(datapoints):
            dp_xyz = dp.cartesian_coordinates
            for j, dp_2 in enumerate(datapoints):
                mol_xyz = dp_2.cartesian_coordinates
                _, dist, _ = self.calculate_distance_to_ref(mol_xyz, dp_xyz, symmetry_info)
                D_db[i, j] = float(dist)
        
        for i, mol_1 in enumerate(molecules):
            mol1_xyz = mol_1.get_coordinates_in_bohr()
            for j, mol_2 in enumerate(molecules):
                mol2_xyz = mol_2.get_coordinates_in_bohr()
                _, dist, _ = self.calculate_distance_to_ref(mol1_xyz, mol2_xyz, symmetry_info)
                D_ref[i, j] = float(dist)

        for i, dp in enumerate(datapoints):
            dp_xyz = dp.cartesian_coordinates
            for j, mol_2 in enumerate(molecules):
                mol2_xyz = mol_2.get_coordinates_in_bohr()
                _, dist, _ = self.calculate_distance_to_ref(dp_xyz, mol2_xyz, symmetry_info)
                D_dpref[i, j] = float(dist)


        return D_db, D_ref, D_dpref
    
    def _clone_datapoints_with_alphas(self, datapoints, alphas):
        cloned = [copy.deepcopy(dp) for dp in datapoints]
        alpha_vec = np.asarray(alphas, dtype=np.float64).reshape(-1)
        if len(cloned) != alpha_vec.size:
            raise ValueError(
                f"Alpha/datapoint size mismatch: n_dp={len(cloned)} vs n_alpha={alpha_vec.size}"
            )
        for i, dp in enumerate(cloned):
            dp.confidence_radius = float(alpha_vec[i])
        return cloned
    
    def _build_interp_driver_for_trust_radius_eval(
        self,
        *,
        z_matrix,
        interpolation_setting,
        sym_dict,
        sym_datapoints,
        exponent_p_q,
        datapoints_eval,
        cluster_banks=None,
    ):
        drv = InterpolationDriver(z_matrix)
        drv.update_settings(interpolation_setting)
        drv.use_symmetry = bool(sym_datapoints)
        drv.symmetry_information = sym_dict
        drv.qm_symmetry_data_points = sym_datapoints
        drv.qm_rotor_cluster_banks = cluster_banks or {}
        if drv.qm_rotor_cluster_banks:
            first_family = next(iter(drv.qm_rotor_cluster_banks.values()))
            drv.rotor_cluster_information = first_family.get("cluster_info")
        else:
            drv.rotor_cluster_information = None
        drv.print = False
        drv.calc_optim_trust_radius = False
        drv.exponent_p = exponent_p_q[0]
        drv.exponent_q = exponent_p_q[1]

        drv.impes_coordinate.eq_bond_lengths = datapoints_eval[0].eq_bond_lengths

        if len(datapoints_eval) > 0 and hasattr(datapoints_eval[0], "inv_sqrt_masses"):
            drv.impes_coordinate.inv_sqrt_masses = datapoints_eval[0].inv_sqrt_masses
        elif getattr(self, "inv_sqrt_masses", None) is not None:
            drv.impes_coordinate.inv_sqrt_masses = self.inv_sqrt_masses

        if getattr(self, "eq_bond_force_constants", None) is not None:
            drv.eq_bond_force_constants = self.eq_bond_force_constants
            drv.impes_coordinate.eq_bond_force_constants = self.eq_bond_force_constants

        drv.qm_data_points = datapoints_eval
        drv.prepare_runtime_data_cache(force=True)
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
        sym_datapoints,
        exponent_p_q,
        cluster_banks=None,
    ):
        if len(molecules) == 0:
            return {
                "scenario": scenario_name,
                "n_dp": int(len(datapoints_eval)),
                "n_ref": 0,
                "energy_errors_kcal": [],
                "gradient_rmsd_kcal_per_mol_A": [],
                "summary": {
                    "energy_mae": math.nan,
                    "energy_rmse": math.nan,
                    "energy_max_abs": math.nan,
                    "gradient_mean_rmsd": math.nan,
                    "gradient_rmse_rmsd": math.nan,
                    "gradient_max_rmsd": math.nan,
                },
            }

        drv = self._build_interp_driver_for_trust_radius_eval(
            z_matrix=z_matrix,
            interpolation_setting=interpolation_setting,
            sym_dict=sym_dict,
            sym_datapoints=sym_datapoints,
            exponent_p_q=exponent_p_q,
            datapoints_eval=datapoints_eval,
            cluster_banks=cluster_banks,
        )

        e_conv = hartree_in_kcalpermol()
        g_conv = hartree_in_kcalpermol() * bohr_in_angstrom()

        e_err = []
        g_rmsd = []

        for i, mol in enumerate(molecules):
            self._interp_compute_root_local_serial(drv, mol)
            e_im = float(drv.get_energy())
            g_im = np.asarray(drv.get_gradient(), dtype=np.float64)

            e_qm = float(np.asarray(qm_energies[i], dtype=np.float64).reshape(-1)[0])
            g_qm = np.asarray(qm_gradients[i], dtype=np.float64)

            de = (e_im - e_qm) * e_conv
            dg = (g_im - g_qm) * g_conv
            rmsd = float(np.sqrt(np.mean(dg**2)))

            e_err.append(float(de))
            g_rmsd.append(rmsd)

        e_arr = np.asarray(e_err, dtype=np.float64)
        g_arr = np.asarray(g_rmsd, dtype=np.float64)

        return {
            "scenario": scenario_name,
            "n_dp": int(len(datapoints_eval)),
            "n_ref": int(len(molecules)),
            "energy_errors_kcal": e_arr.tolist(),
            "gradient_rmsd_kcal_per_mol_A": g_arr.tolist(),
            "summary": {
                "energy_mae": float(np.mean(np.abs(e_arr))),
                "energy_rmse": float(np.sqrt(np.mean(e_arr**2))),
                "energy_max_abs": float(np.max(np.abs(e_arr))),
                "gradient_mean_rmsd": float(np.mean(g_arr)),
                "gradient_rmse_rmsd": float(np.sqrt(np.mean(g_arr**2))),
                "gradient_max_rmsd": float(np.max(g_arr)),
            },
        }
    
    def _print_trust_radius_stage_comparison(self, stage_results):
        order = [
            "without_last_datapoint",
            "with_new_datapoint_pre_opt",
            "with_new_datapoint_post_opt",
        ]
        by_name = {row["scenario"]: row for row in stage_results}

        print("[trust-radius-eval] Interpolation comparison on fixed reference set")
        print(
            "scenario".ljust(32),
            "n_dp".rjust(6),
            "E_MAE".rjust(12),
            "E_RMSE".rjust(12),
            "E_MAX".rjust(12),
            "G_MEAN".rjust(12),
            "G_RMSE".rjust(12),
            "G_MAX".rjust(12),
        )

        for key in order:
            row = by_name.get(key)
            if row is None:
                continue
            s = row["summary"]
            print(
                key.ljust(32),
                f"{row['n_dp']:6d}",
                f"{s['energy_mae']:12.5f}",
                f"{s['energy_rmse']:12.5f}",
                f"{s['energy_max_abs']:12.5f}",
                f"{s['gradient_mean_rmsd']:12.5f}",
                f"{s['gradient_rmse_rmsd']:12.5f}",
                f"{s['gradient_max_rmsd']:12.5f}",
            )

        def _delta(a_name, b_name, field):
            a = by_name.get(a_name, {}).get("summary", {}).get(field, math.nan)
            b = by_name.get(b_name, {}).get("summary", {}).get(field, math.nan)
            return b - a

        print(
            "[trust-radius-eval][delta] pre_opt - without_last:",
            f"dE_MAE={_delta('without_last_datapoint', 'with_new_datapoint_pre_opt', 'energy_mae'):.6f}",
            f"dG_MEAN={_delta('without_last_datapoint', 'with_new_datapoint_pre_opt', 'gradient_mean_rmsd'):.6f}",
        )
        print(
            "[trust-radius-eval][delta] post_opt - pre_opt:",
            f"dE_MAE={_delta('with_new_datapoint_pre_opt', 'with_new_datapoint_post_opt', 'energy_mae'):.6f}",
            f"dG_MEAN={_delta('with_new_datapoint_pre_opt', 'with_new_datapoint_post_opt', 'gradient_mean_rmsd'):.6f}",
        )
    
    def _resolve_trust_radius_cluster_bank(self, root):
        if not self.cluster_run:
            return {}
        bank = self.qm_rotor_cluster_banks.get(root)
        if not bank:
            bank = self._load_rotor_cluster_bank_for_root(root)
            self.qm_rotor_cluster_banks[root] = bank
        return bank or {}


    def _subset_rotor_cluster_bank_for_datapoints(self, cluster_banks, datapoints):
        if not isinstance(cluster_banks, dict) or len(cluster_banks) == 0:
            return {}

        subset = {}
        for dp in datapoints:
            family_label = getattr(dp, "family_label", None) or getattr(dp, "point_label", None)
            if family_label in cluster_banks and family_label not in subset:
                subset[family_label] = cluster_banks[family_label]
        return subset

    def determine_trust_radius_gradient(self, molecules, qm_energies, qm_gradients, im_energies, datapoints, interpolation_setting, sym_datapoints, sym_dict, z_matrix, exponent_p_q):
        
        def optimize_trust_radius(alphas, geom_list, E_ref_list, G_ref_list, E_im_list, dps, impes_dict, sym_datapoints, sym_dict, exponent_p_q, cluster_banks, basin_id=None):
            """
            Perform the gradient-based optimization to find R*
            that minimizes the sum of squared errors to reference QM energies.
            """

            bounds = [(1e-4, 10.5)] * len(dps)
            sample_size = len(geom_list)

            if sample_size < 400:
                size_label = 'small'
                n_bh_iter = 10
                gtol = 1e-3
                ftol = 1e-4
                maxiter = 300
            elif sample_size <= 800:
                size_label = 'medium'
                n_bh_iter = 4
                gtol = 2e-3
                ftol = 5e-4
                maxiter = 220
            else:
                size_label = 'large'
                n_bh_iter = 0
                gtol = 5e-3
                ftol = 5e-4
                maxiter = 140

            # Backward compatible interpretation:
            # [enabled, mode, init_confidence_radius, e_x]
            # Older input may provide only three entries; then fall back to idx 2.
            if isinstance(self.use_opt_confidence_radius, (list, tuple)) and len(self.use_opt_confidence_radius) > 3:
                energy_weight = float(self.use_opt_confidence_radius[3])
            else:
                energy_weight = float(self.use_opt_confidence_radius[2])

            opt = AlphaOptimizer(z_matrix, impes_dict, sym_dict, sym_datapoints, dps,
                 geom_list, E_ref_list, G_ref_list, exponent_p_q,
                 e_x=energy_weight,
                 beta=0.8, n_workers=os.cpu_count(), cluster_banks=cluster_banks)

            history = []
            iter_counter = {"value": 0}
            basin_label = "?"
            if basin_id is not None:
                basin_label = str(int(basin_id))

            def _fmt_float(value):
                if value is None:
                    return "nan"
                try:
                    v = float(value)
                except Exception:
                    return "nan"
                if not np.isfinite(v):
                    return "nan"
                return f"{v:.6e}"

            def _append_history(tag, x, **extra):
                x_arr = np.asarray(x, dtype=np.float64).reshape(-1)
                metrics = opt.get_metrics(x_arr, evaluate_if_missing=True)
                entry = {
                    "tag": str(tag),
                    "iteration": int(iter_counter["value"]),
                    "alphas": x_arr.tolist(),
                    "loss_total": float(metrics["loss_total"]) if metrics else math.nan,
                    "loss_energy": float(metrics["loss_energy"]) if metrics else math.nan,
                    "loss_force": float(metrics["loss_force"]) if metrics else math.nan,
                }
                try:
                    jac = opt.jac(x_arr)
                    entry["jac_l2"] = float(np.linalg.norm(jac))
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
                print(
                    f"[trust-opt][basin {basin_label}][{entry['tag']}][iter={entry['iteration']}] "
                    f"L={_fmt_float(entry['loss_total'])} "
                    f"Le={_fmt_float(entry['loss_energy'])} "
                    f"Lf={_fmt_float(entry['loss_force'])} "
                    f"|g|={_fmt_float(entry['jac_l2'])}"
                )

            def _local_minimizer_callback(xk):
                iter_counter["value"] += 1
                _append_history("local_iter", xk)

            def _basinhop_callback(x, f, accept):
                _append_history(
                    "basinhop_step",
                    x,
                    basinhop_fun=float(f),
                    basinhop_accept=bool(accept))
                return False

            minimizer_kwargs = {
                "method": "L-BFGS-B",
                "jac": opt.jac,
                "bounds": bounds,
                "callback": _local_minimizer_callback,
                "options": {"gtol": gtol, "ftol": ftol, "maxls": 10, "maxiter": maxiter, "disp": False}
            }

            print(
                f"Trust-radius optimization regime: {size_label} (S={sample_size}, basinhopping niter={n_bh_iter})"
            )
            try:
                _append_history("initial", alphas)
                if n_bh_iter > 0:
                    res = basinhopping(
                        opt.fun,
                        x0=alphas,
                        minimizer_kwargs=minimizer_kwargs,
                        niter=n_bh_iter,
                        callback=_basinhop_callback,
                        disp=False)
                else:
                    res = minimize(
                        opt.fun,
                        x0=alphas,
                        jac=opt.jac,
                        method="L-BFGS-B",
                        bounds=bounds,
                        callback=_local_minimizer_callback,
                        options=minimizer_kwargs["options"])
                _append_history("final", res.x)
                setattr(res, "history", history)
                setattr(res, "component_summary", opt.get_metrics(res.x, evaluate_if_missing=True))
                print(
                    "Trust-radius optimization result:",
                    f"fun={float(res.fun):.6f}",
                    f"x={np.asarray(res.x, dtype=np.float64)}")
                
                return res
            finally:
                opt.close()

        def _cluster_banks_for(dp_subset):
            return self._subset_rotor_cluster_bank_for_datapoints(cluster_banks_all, dp_subset)

        self._emit_test_hook("alpha_optimization_start", {
            "n_structures": int(len(molecules)),
            "n_datapoints": int(len(datapoints)),
            "mode": str(self.use_opt_confidence_radius[1]),
        })
        
        initial_alphas = np.array(
            [float(np.asarray(dp.confidence_radius).reshape(-1)[0]) for dp in datapoints],
            dtype=float
        )
        print('INPUT Trust radius', initial_alphas)

        if isinstance(self.use_opt_confidence_radius, (list, tuple)) and len(self.use_opt_confidence_radius) > 3:
            energy_weight = float(self.use_opt_confidence_radius[3])
        else:
            energy_weight = float(self.use_opt_confidence_radius[2])

        
        
        has_last = len(datapoints) >= 2
        if has_last:
            dp_wo_last = datapoints[:-1]
            alphas_wo_last = initial_alphas[:-1]
        else:
            dp_wo_last = datapoints
            alphas_wo_last = initial_alphas
        
        
        stage_results = []
        cluster_banks_all = self._resolve_trust_radius_cluster_bank(self.current_state)     

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
                sym_datapoints=sym_datapoints,
                exponent_p_q=exponent_p_q,
                cluster_banks=_cluster_banks_for(dp_eval_a),
            )
        )

        dp_eval_b = self._clone_datapoints_with_alphas(datapoints, initial_alphas)
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
                sym_datapoints=sym_datapoints,
                exponent_p_q=exponent_p_q,
                cluster_banks=_cluster_banks_for(dp_eval_b),
            )
        )  
        
        D_dp, D_ref, D_dpref = self._build_dp_ref_distance_matrix(datapoints, molecules, sym_dict)


        groups = improved_group_by_connected_components(D_dp, D_ref, D_dpref)

        final_alphas = initial_alphas.copy()

        basin_results = []
        basin_payloads = []

        for group in groups['groups']:
            
            dp_mask = group['datapoint_indices']
            ref_mask = group['reference_indices']

            datapoints_sub = [datapoints[i] for i in dp_mask]
            molecules_sub = [molecules[i] for i in ref_mask]
            qm_energies_sub = [qm_energies[i] for i in ref_mask]
            qm_gradients_sub = [qm_gradients[i] for i in ref_mask]
            im_energies_sub = [im_energies[i] for i in ref_mask]
            alphas_sub = [0.5 for i in dp_mask]

            cluster_banks_sub = _cluster_banks_for(datapoints_sub)
                    

            res = optimize_trust_radius(
                alphas_sub,
                molecules_sub,
                qm_energies_sub,
                qm_gradients_sub,
                im_energies_sub,
                datapoints_sub,
                interpolation_setting,
                sym_datapoints,
                sym_dict,
                exponent_p_q,
                cluster_banks=cluster_banks_sub,
                basin_id=group["basin_id"])

            final_alphas[dp_mask] = res['x']
            print('FINAL Trust radius', res['x'])

            self._emit_test_hook("alpha_optimization_end", {
                "optimized_alphas": [float(x) for x in final_alphas],
            })

            basin_results.append({
                "basin_id": group["basin_id"],
                "n_datapoints": len(group["datapoint_indices"]),
                "n_references": len(group["reference_indices"]),
                "result": res,
            })
            basin_payloads.append({
                "basin_id": int(group["basin_id"]),
                "dp_mask": np.asarray(dp_mask, dtype=int),
                "molecules_sub": molecules_sub,
                "qm_energies_sub": qm_energies_sub,
                "qm_gradients_sub": qm_gradients_sub,
                "datapoints_sub": datapoints_sub,
                "cluster_banks_sub": cluster_banks_sub,
            })
    
        sensitivity_optimizers = []
        objective_fn = None
        if len(basin_payloads) > 0:
            n_workers_sens = max(1, min(4, os.cpu_count() or 1))

            for payload in basin_payloads:
                sens_opt = AlphaOptimizer(
                    z_matrix,
                    interpolation_setting,
                    sym_dict,
                    sym_datapoints,
                    payload["datapoints_sub"],
                    payload["molecules_sub"],
                    payload["qm_energies_sub"],
                    payload["qm_gradients_sub"],
                    exponent_p_q,
                    e_x=energy_weight,
                    beta=0.8,
                    n_workers=n_workers_sens,
                    verbose=False,
                    cluster_banks=payload["cluster_banks_sub"],
                )
                sensitivity_optimizers.append({
                    "dp_mask": payload["dp_mask"],
                    "optimizer": sens_opt,
                })

            def objective_fn(alpha_global):
                alpha_vec = np.asarray(alpha_global, dtype=np.float64).reshape(-1)
                if alpha_vec.size != len(datapoints):
                    raise ValueError(
                        f"objective_fn expects alpha size {len(datapoints)}, got {alpha_vec.size}")
                if len(sensitivity_optimizers) == 0:
                    return 0.0
                losses = []
                for item in sensitivity_optimizers:
                    sub_alpha = alpha_vec[item["dp_mask"]]
                    losses.append(float(item["optimizer"].fun(sub_alpha)))
                return float(np.mean(losses))

        dp_labels = None
        if (isinstance(getattr(self, 'sorted_state_spec_im_labels', None), dict)
                and self.current_state in self.sorted_state_spec_im_labels):
            cand = self.sorted_state_spec_im_labels[self.current_state]
            if len(cand) == len(datapoints):
                dp_labels = list(cand)

        analyzer = TrustRadiusOptimizationAnalyzer()
        try:
            report = analyzer.analyze(
                initial_alphas=initial_alphas,
                final_alphas=final_alphas,
                groups=groups['groups'],
                d_dpref=D_dpref,
                grouping_diagnostics=groups["diagnostics"],
                basin_optimizer_results=basin_results,
                datapoint_labels=dp_labels,
                objective_fn=objective_fn,
            )
        finally:
            for item in sensitivity_optimizers:
                try:
                    item["optimizer"].close()
                except Exception:
                    pass

        self.last_trust_radius_report = report

        print(summarize_report_for_console(report, top_k=10))
        print(report["summary"])
        print(report["problematic_datapoints"][:5])
        if "basin_optimizer_results" in report:
            print("Trust-radius basin convergence summary")
            for basin_entry in report["basin_optimizer_results"]:
                basin_id = basin_entry.get("basin_id", "?")
                result_info = basin_entry.get("result", {})
                history = result_info.get("history", [])
                comp = result_info.get("component_summary", {})

                loss_start = math.nan
                loss_end = math.nan
                if len(history) > 0:
                    loss_start = float(history[0].get("loss_total", math.nan))
                    loss_end = float(history[-1].get("loss_total", math.nan))

                print(
                    f"[trust-opt][basin {basin_id}] "
                    f"success={result_info.get('success', False)} "
                    f"nit={result_info.get('nit', -1)} "
                    f"nfev={result_info.get('nfev', -1)} "
                    f"hist={len(history)} "
                    f"L0={loss_start:.6e} "
                    f"Lf={loss_end:.6e} "
                    f"Le_f={float(comp.get('loss_energy', math.nan)):.6e} "
                    f"Lf_f={float(comp.get('loss_force', math.nan)):.6e}"
                )
        
        
        
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
                sym_datapoints=sym_datapoints,
                exponent_p_q=exponent_p_q,
                cluster_banks=_cluster_banks_for(dp_eval_c),
            )
        )
        
        self._print_trust_radius_stage_comparison(stage_results)

        return final_alphas
    
    def perform_symmetry_assignment(self, atom_map, sym_group, reference_group, datapoint_group):
        """ Performs the atom mapping. """
        from scipy.optimize import linear_sum_assignment
        new_map = np.array(atom_map.copy())
        mapping_dict = {}
        # cost = self.get_dihedral_cost(atom_map, sym_group, non_group_atoms)
        cost = np.linalg.norm(datapoint_group[:, np.newaxis, :] - reference_group[np.newaxis, :, :], axis=2)
        row, col = linear_sum_assignment(cost)
        assigned = False
        if not np.equal(row, col).all():
            assigned = True
            
            # atom_maps = self.linear_assignment_solver(cost)

            reordred_arr = np.array(sym_group)[col]
            new_map[sym_group] = new_map[reordred_arr]

            mapping_dict = {org: new for org, new in zip(np.array(sym_group), reordred_arr)}
        
        return mapping_dict

    def adjust_symmetry_dihedrals(self, symmetry_groups, rot_bonds, z_matrix):
        
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
                

            elif len(symmetry_group) == 2:

                # angles_to_set[symmetry_group_dihedral_list[0]] = ([0.0, np.pi/3.0])
                angles_to_set[tuple(symmetry_group_dihedral_list[0])] = ([0.0, np.pi/2.0])

                periodicities[tuple(symmetry_group_dihedral_list[0])] = 2
                dihedral_groups[2].extend([tuple(sorted(element)) for element in symmetry_group_dihedral_list])

        return angles_to_set, periodicities, symmetry_group_dihedral_dict, dihedral_groups       

    def _mpi_collective_comm(self):
        return self._mpi_ctrl_comm if self._mpi_is_active() else self.comm

    def _run_optimization_mpi_safe(
        self,
        optimization_driver,
        molecule,
        constraints=None,
        transition=False,
        index_offset=1,
        compute_args=None,
        source_molecule=None,
        collective=False,
        phase_name='optimization tag',
    ):
        if collective:
            comm = self._mpi_collective_comm()
            rank = comm.Get_rank()
            root = mpi_master()
            return _collective_call_inline_style(
                comm,
                rank,
                root,
                phase_name,
                self._run_optimization,
                optimization_driver,
                molecule,
                constraints,
                transition,
                index_offset,
                compute_args,
                source_molecule,
                True,
            )

        return self._run_optimization(
            optimization_driver,
            molecule,
            constraints=constraints,
            transition=transition,
            index_offset=index_offset,
            compute_args=compute_args,
            source_molecule=source_molecule,
            collective=False,
        )

    def _run_optimization(
        self,
        optimization_driver,
        molecule,
        constraints=None,
        transition=False,
        index_offset=1,
        compute_args=None,
        source_molecule=None,
        collective=False,
    ):
        opt_drv = OptimizationDriver(optimization_driver)
        opt_drv.ostream.mute()
        opt_drv.transition = transition

        if constraints is not None:
            opt_drv.constraints = self._build_opt_constraint_list(constraints, index_offset=index_offset)

        if compute_args is None:
            if collective:
                opt_results = opt_drv.compute(molecule)
            else:
                opt_results = self._opt_compute_mpi_safe(opt_drv, molecule)
        else:
            if collective:
                opt_results = opt_drv.compute(molecule, *compute_args)
            else:
                opt_results = self._opt_compute_mpi_safe(opt_drv, molecule, *compute_args)

        if source_molecule is None:
            source_molecule = molecule

        optimized_molecule = Molecule.from_xyz_string(opt_results['final_geometry'])
        optimized_molecule.set_charge(source_molecule.get_charge())
        optimized_molecule.set_multiplicity(source_molecule.get_multiplicity())

        return optimized_molecule, opt_results

    def _compute_energy_mpi_safe(self, qm_driver, molecule, basis=None, collective=False, phase_name='Energy calcualtion'):
        if collective:
            comm = self._mpi_collective_comm()
            rank = comm.Get_rank()
            root = mpi_master()
            return _collective_call_inline_style(
                comm, rank, root, phase_name, self.compute_energy, qm_driver, molecule, basis
            )

        with self._mpi_root_local_qm_context(extra_objects=[qm_driver]):
            if self._mpi_is_active() and self.mpi_root_worker_mode and self._mpi_is_root() and hasattr(qm_driver, "comm"):
                if qm_driver.comm.Get_size() != 1:
                    raise RuntimeError("Root-local MPI override failed for QM energy driver (comm size != 1).")
            return self.compute_energy(qm_driver, molecule, basis)

    def _compute_gradient_mpi_safe(
        self, grad_driver, molecule, basis=None, scf_results=None, rsp_results=None,
        collective=False, phase_name='Gradient Calculation'
    ):
        if collective:
            comm = self._mpi_collective_comm()
            rank = comm.Get_rank()
            root = mpi_master()
            return _collective_call_inline_style(
                comm, rank, root, phase_name,
                self.compute_gradient, grad_driver, molecule, basis, scf_results, rsp_results
            )

        with self._mpi_root_local_qm_context(extra_objects=[grad_driver]):
            if self._mpi_is_active() and self.mpi_root_worker_mode and self._mpi_is_root() and hasattr(grad_driver, "comm"):
                if grad_driver.comm.Get_size() != 1:
                    raise RuntimeError("Root-local MPI override failed for QM gradient driver (comm size != 1).")
            return self.compute_gradient(grad_driver, molecule, basis, scf_results, rsp_results)

    def _compute_hessian_mpi_safe(self, hess_driver, molecule, basis=None, collective=False, phase_name='hessian calculation'):
        if collective:
            comm = self._mpi_collective_comm()
            rank = comm.Get_rank()
            root = mpi_master()
            return _collective_call_inline_style(
                comm, rank, root, phase_name, self.compute_hessian, hess_driver, molecule, basis
            )

        with self._mpi_root_local_qm_context(extra_objects=[hess_driver]):
            if self._mpi_is_active() and self.mpi_root_worker_mode and self._mpi_is_root() and hasattr(hess_driver, "comm"):
                if hess_driver.comm.Get_size() != 1:
                    raise RuntimeError("Root-local MPI override failed for QM hessian driver (comm size != 1).")
            return self.compute_hessian(hess_driver, molecule, basis)
    
    def compute_energy(self, qm_driver, molecule, basis=None):
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
            qm_energy = qm_driver.scf_energy
            qm_energy = np.array([qm_energy])
            qm_driver.ostream.unmute()
            qm_driver.filename = None
            qm_driver.checkpoint_file = None
        
        elif isinstance(qm_driver, ScfUnrestrictedDriver):
            qm_driver.ostream.mute()
            scf_results = qm_driver.compute(molecule, basis)
            qm_energy = qm_driver.scf_energy
            qm_energy = np.array([qm_energy])
            qm_driver.ostream.unmute()
            qm_driver.filename = None
            qm_driver.checkpoint_file = None
        
        elif isinstance(qm_driver, LinearResponseEigenSolver) or isinstance(qm_driver, TdaEigenSolver):
            self.drivers['es'][0].ostream.mute()

            scf_results = self.drivers['es'][3].compute(molecule, basis)
            self.drivers['es'][3].filename = None
            self.drivers['es'][3].checkpoint_file = None
            scf_energy = self.drivers['es'][3].scf_energy
            qm_driver.ostream.mute()
            rsp_results = qm_driver.compute(molecule, basis, scf_results)
            
            qm_energy = np.insert(rsp_results['eigenvalues'] + scf_energy, 0, scf_energy)

        if qm_energy is None:
            error_txt = "Could not compute the QM energy. "
            error_txt += "Please define a QM driver."
            raise ValueError(error_txt)

        return qm_energy, scf_results, rsp_results

    def compute_gradient(self, grad_driver, molecule, basis=None, scf_results=None, rsp_results=None):
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
    def compute_hessian(self, hess_driver, molecule, basis=None):
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
                print('In hessian clac', root, hess_driver)
                self.drivers['es'][1].state_deriv_index = root
                hess_driver.compute(molecule, basis)
                qm_hessian = hess_driver.hessian
                qm_hessians.append(qm_hessian)

        if qm_hessians is None:
            error_txt = "Could not compute the QM Hessian. "
            error_txt += "Please define a QM Hessian driver."
            raise ValueError(error_txt)


        return qm_hessians



    def get_qm_potential_energy(self):
        """
        Returns the potential energy of the QM region.

        Args:
            context: The OpenMM context object.
        Returns:
            The potential energy of the QM region.
        """

        potential_energy = self.current_energy

        return potential_energy
    
    def output_file_writer(self, outputfile):
        """
        Writes the current simulation summary data (stored in self.gloabal_sim_informations)
        into the HDF5 file and then resets the dictionary.
        """
        
        valid_checkpoint = (outputfile and isinstance(outputfile, str))

        if valid_checkpoint:
            try:
                h5f = h5py.File(outputfile, 'a')
            except IOError:
                h5f = h5py.File(outputfile, 'w')
        

        # 1. Write temperatures
        ds_name = 'mol_labels'
        if ds_name in h5f:
            print('labels already exist')
        else:
            labels = self.molecule.get_labels()
            dt = h5py.string_dtype(encoding='utf-8')
            label_ds = h5f.create_dataset(
                ds_name,
                shape=(0,),
                maxshape=(None,),
                dtype=dt,
                chunks=True
            )
            # old_size = label_ds.shape[0]
            # new_size = old_size + labels.shape[0]
            # label_ds.resize((new_size,))
            # label_ds[old_size:new_size] = labels

        ds_name = "temperatures"
        if ds_name in h5f:
            temp_ds = h5f[ds_name]
        else:
            temp_ds = h5f.create_dataset(
                ds_name,
                shape=(0,),
                maxshape=(None,),
                dtype='float64',
                chunks=True
            )
        # Convert the list to a NumPy array and append
        new_temps = np.array(self.gloabal_sim_informations['temperatures'], dtype='float64')
        old_size = temp_ds.shape[0]
        new_size = old_size + new_temps.shape[0]
        temp_ds.resize((new_size,))
        temp_ds[old_size:new_size] = new_temps
        # 2. Write state labels
        ds_name = "state"
        # Adjust dtype as needed (here assuming integer states; if strings, you'll need a special dtype)
        if ds_name in h5f:
            state_ds = h5f[ds_name]
        else:
            state_ds = h5f.create_dataset(
                ds_name,
                shape=(0,),
                maxshape=(None,),
                dtype='int32',
                chunks=True
            )
        new_states = np.array(self.gloabal_sim_informations['state'], dtype='int32')
        old_size = state_ds.shape[0]
        new_size = old_size + new_states.shape[0]
        state_ds.resize((new_size,))
        state_ds[old_size:new_size] = new_states
        
        # 3. Write coordinates (assumed to be np.array of shape (n_atoms, 3))
        ds_name = "coordinates_ang"
        coord_list = self.gloabal_sim_informations['coordinates_ang']
        if len(coord_list) > 0:
            # Determine the number of atoms from the first snapshot.
            n_atoms = coord_list[0].shape[0]
            new_coords = np.stack(coord_list, axis=0)  # shape: (n_snapshots, n_atoms, 3)
            if ds_name in h5f:
                coord_ds = h5f[ds_name]
            else:
                coord_ds = h5f.create_dataset(
                    ds_name,
                    data=new_coords,
                    maxshape=(None, n_atoms, 3),
                    dtype=new_coords.dtype,
                    chunks=True
                )
                new_coords = None  # Already stored
            if new_coords is not None:
                old_size = coord_ds.shape[0]
                new_size = old_size + new_coords.shape[0]
                coord_ds.resize((new_size, n_atoms, 3))
                coord_ds[old_size:new_size, :, :] = new_coords
        # 4. Write state-specific data (for each state in self.roots_to_follow)
        for root in self.roots_to_follow:
            key = f'state_{root}'
            state_dict = self.gloabal_sim_informations[key]
            state_group = h5f.require_group(key)
            
            # 1. Store potential energies (assumed to be 1D arrays)
            for subkey in ['pot_energies']:
                if subkey in state_group:
                    ds = state_group[subkey]
                else:
                    ds = state_group.create_dataset(
                        subkey,
                        shape=(0,),
                        maxshape=(None,),
                        dtype='float64',
                        chunks=True
                    )
                new_data = np.array(state_dict[subkey], dtype='float64')
                old_size = ds.shape[0]
                new_size = old_size + new_data.shape[0]
                ds.resize((new_size,))
                ds[old_size:new_size] = new_data
            
            # 2. Store gradients (each entry is an array of shape (n_atoms, 3))
            subkey = 'gradients'
            new_gradients = state_dict[subkey]
            if new_gradients:  # Only proceed if there is new data
                # Determine number of atoms from the first gradient array.
                n_atoms = new_gradients[0].shape[0]
                # Stack new gradients along a new axis: shape becomes (n_snapshots, n_atoms, 3)
                grad_data = np.stack(new_gradients, axis=0)
                
                if subkey in state_group:
                    ds = state_group[subkey]
                else:
                    # Create the dataset with initial data
                    ds = state_group.create_dataset(
                        subkey,
                        data=grad_data,
                        maxshape=(None, n_atoms, 3),
                        dtype=grad_data.dtype,
                        chunks=True
                    )
                    grad_data = None  # Already stored
                # If the dataset already existed, append the new data.
                if grad_data is not None:
                    old_size = ds.shape[0]
                    new_size = old_size + grad_data.shape[0]
                    ds.resize((new_size, n_atoms, 3))
                    ds[old_size:new_size, :, :] = grad_data

            # store information about the checked data with reference calculations

            
    
    def read_simulation_summary(self, fname, roots_to_follow):
        """
        Reads a simulation summary from an HDF5 file.

        :param fname: Name of the HDF5 file.
        :param roots_to_follow: List of root states (e.g., [0, 1]) that were written.
        :return: A dictionary containing all read data.
        """
        data = {}

        with h5py.File(fname, 'r') as h5f:
            # 1. Read temperatures
            if 'temperatures' in h5f:
                data['temperatures'] = np.array(h5f['temperatures'])

            # 2. Read states
            if 'state' in h5f:
                data['state'] = np.array(h5f['state'])

            # 3. Read recoordinates
            if 'coordinates_ang' in h5f:
                data['coordinates_ang'] = np.array(h5f['coordinates_ang'])

            # 4. Read state-specific groups
            for root in roots_to_follow:
                state_key = f'state_{root}'
                if state_key in h5f:
                    state_group = h5f[state_key]
                    state_data = {}

                    if 'pot_energies' in state_group:
                        state_data['pot_energies'] = np.array(state_group['pot_energies'])

                    if 'gradients' in state_group:
                        state_data['gradients'] = np.array(state_group['gradients'])  # shape: (n_snapshots, n_atoms, 3)

                    data[state_key] = state_data

        return data
    
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
        Append reference structures + metadata to an HDF5 file.

        Notes:
        - This replaces writing XYZ text blocks.
        - Stores xyz strings + flattened coords and gradients.
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
        
        h5_path = h5_path or self.reference_struc_energies_file
        state = int(state)

        # Basic sanity
        n = len(molecules)
        if not (len(qm_energies) == len(qm_gradients) == len(im_energies) == n):
            raise ValueError("Input lists must all have the same length.")

        with h5py.File(h5_path, "a") as h5:
            g = _ensure_state_group(h5, state)

            # current length
            N0 = g["qm_energy"].shape[0]
            N1 = N0 + n

            # resize all 1D datasets to accommodate new rows
            for key in ["qm_energy", "im_energy", "natoms", "xyz", "qm_grad_flat", "coords_flat"]:
                g[key].resize((N1,))

            # append row-by-row (ragged arrays)
            for i, mol in enumerate(molecules):
                # xyz for easy reconstruction later
                xyz = mol.get_xyz_string()

                # coordinates if available; optional but recommended
                coords = None
                if hasattr(mol, "get_coordinates"):
                    coords = np.asarray(mol.get_coordinates_in_angstrom(), dtype=np.float64)
                else:
                    # Fallback: parse from xyz string if needed (cheap enough, but implement if you want)
                    coords = None

                # natoms: prefer mol attribute, otherwise infer from gradient/coords/xyz
                if hasattr(mol, "n_atoms"):
                    natoms = int(mol.n_atoms)
                elif coords is not None:
                    natoms = int(coords.shape[0])
                else:
                    # xyz format: first line is natoms
                    natoms = int(xyz.splitlines()[0].strip())

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

                if coords is None:
                    # Store empty array if coords are not accessible; reader can still rebuild from xyz
                    g["coords_flat"][row] = np.array([], dtype=np.float64)
                else:
                    if coords.shape != (natoms, 3):
                        raise ValueError(f"Coords must be shape (natoms,3), got {coords.shape}.")
                    g["coords_flat"][row] = coords.reshape(-1)

    def extract_reference_structures_h5(self, filename, roots):
        """
        Read reference structures from HDF5 and return the same dict structure as your
        current extract_reference_structures(...).

        returns:
        { root/state: {
            'molecules': [...],
            'qm_energies': [...],
            'qm_gradients': [...],   # arrays (natoms,3)
            'im_energies': [...]
            }, ...}
        """
        reference_molecules_check = {
            int(root): {"molecules": [], "qm_energies": [], "qm_gradients": [], "im_energies": []}
            for root in roots
        }

        with h5py.File(filename, "r") as h5:
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
                xyzs = g["xyz"][:]               # array of strings
                grad_flat = g["qm_grad_flat"][:] # array of vlen float arrays

                for i in range(len(qmE)):
                    # Rebuild molecule
                    xyz_string = xyzs[i]
                    if isinstance(xyz_string, bytes):
                        xyz_string = xyz_string.decode("utf-8")

                    mol = Molecule.from_xyz_string(xyz_string)

                    # Rebuild gradient
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
        """Center the molecule by translating its geometric center to (0, 0, 0)."""
        center = np.mean(given_coordinates, axis=0)
        translated_coordinates = given_coordinates - center

        return translated_coordinates
    

    def calculate_distance_to_ref(self, current_coordinates, datapoint_coordinate, symmetry_info):
        """Calculates and returns the cartesian distance between
           self.coordinates and data_point coordinates.
           Besides the distance, it also returns the weight gradient,
           which requires the distance vector to be computed.

           :param data_point:
                InterpolationDatapoint object
        """

        active_atoms = np.delete(np.arange(datapoint_coordinate.shape[0]), symmetry_info[4])

        target_coordinates_core = datapoint_coordinate[active_atoms]
        reference_coordinates_core = current_coordinates[active_atoms]

        # First, translate the cartesian coordinates to zero
        target_coordinates = self.calculate_translation_coordinates(target_coordinates_core)
        reference_coordinates = self.calculate_translation_coordinates(reference_coordinates_core)

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

        distance_vector_core = (ref_structure_check - rotated_coordinates_core)
        distance_vector_core_norm = np.zeros(reference_coordinates.shape[0])

        for i in range(len(distance_vector_core_norm)):
            distance_vector_core_norm[i] += np.linalg.norm(distance_vector_core[i])

        return ref_structure_check, distance_core, distance_vector_core
