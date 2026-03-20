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
import itertools
import re
import os, copy, math
import torch
import gpytorch
from pathlib import Path
from sys import stdout
import sys
import random
from time import time
import xml.etree.ElementTree as ET
from xml.dom import minidom
from contextlib import redirect_stderr, contextmanager
from io import StringIO

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
from .gprinterpolationdriver import GPRInterpolationDriver

with redirect_stderr(StringIO()) as fg_err:
    import geometric

try:
    import openmm as mm
    import openmm.app as app
    import openmm.unit as unit
    from openmm.app.metadynamics import Metadynamics, BiasVariable
except ImportError:
    pass


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


# --- put the custom reporter here ---
class GhostSafePDBReporter(app.PDBReporter):
    def report(self, simulation, state):
        positions = state.getPositions(asNumpy=True)
        n_atoms = simulation.topology.getNumAtoms()
        trimmed_positions = positions[:n_atoms]
        app.PDBFile.writeModel(simulation.topology, trimmed_positions, self._out, self._nextModel)
        self._nextModel += 1


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
        self.force_orient_thrsh = None
        self.collect_qm_points = None
        self.previous_energy_list = []
        self.non_core_symmetry_groups = []
        self.interpolation_settings = None
        self.allowed_molecules = None
        self.reference_struc_energies_file = None
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
        self.sampling_qm_symmetry_datapoint_dict = None
        

        # output_file variables that will be written into self.general_variable_output
        self.gradients = None
        self.velocities = None
        self.start_velocities = None
        self.coordinates = None
        self.coordinates_xyz = None
        self.state_specific_molecules = None
        self.all_gradients = []
        self.inv_sqrt_masses = None
        
        self.identfy_relevant_int_coordinates = (True, [])
        self.use_symmetry = True
        self.use_opt_confidence_radius = True
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

        # roots requiring dataset reload on worker ranks
        self._mpi_pending_sync_roots = set()
        

    def _create_platform(self):
        """
        Creates an OpenMM platform.

        Returns:
            OpenMM Platform.
        """

        if self.platform is None:
            return None
        else:
            platform = mm.Platform.getPlatformByName(self.platform)
            if self.openmm_platform == "CPU":
                platform.setPropertyDefaultValue("Threads", "1")
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
            if cmd != self._MPI_CMD_STEP:
                raise RuntimeError(f"Unknown MPI command {cmd} on worker rank {self.rank}")

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

    def _compute_energy_mpi_safe(self, qm_driver, molecule, basis=None):
        with self._mpi_root_local_qm_context(extra_objects=[qm_driver]):
            if self._mpi_is_active() and self.mpi_root_worker_mode and self._mpi_is_root() and hasattr(qm_driver, "comm"):
                if qm_driver.comm.Get_size() != 1:
                    raise RuntimeError("Root-local MPI override failed for QM energy driver (comm size != 1).")
            return self.compute_energy(qm_driver, molecule, basis)

    def _compute_gradient_mpi_safe(self, grad_driver, molecule, basis=None, scf_results=None, rsp_results=None):
        with self._mpi_root_local_qm_context(extra_objects=[grad_driver]):
            if self._mpi_is_active() and self.mpi_root_worker_mode and self._mpi_is_root() and hasattr(grad_driver, "comm"):
                if grad_driver.comm.Get_size() != 1:
                    raise RuntimeError("Root-local MPI override failed for QM gradient driver (comm size != 1).")
            return self.compute_gradient(grad_driver, molecule, basis, scf_results, rsp_results)

    def _compute_hessian_mpi_safe(self, hess_driver, molecule, basis=None):
        with self._mpi_root_local_qm_context(extra_objects=[hess_driver]):
            if self._mpi_is_active() and self.mpi_root_worker_mode and self._mpi_is_root() and hasattr(hess_driver, "comm"):
                if hess_driver.comm.Get_size() != 1:
                    raise RuntimeError("Root-local MPI override failed for QM hessian driver (comm size != 1).")
            return self.compute_hessian(hess_driver, molecule, basis)

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

    def compute_tetrahedral_anchor(self, positions, Cbeta, Calpha, Hbeta1, Hbeta2, bond_length=0.109):
        def normalize(v):
            return v / np.linalg.norm(v)

        rC  = positions[Cbeta].value_in_unit(unit.nanometer)
        rCa = positions[Calpha].value_in_unit(unit.nanometer)
        rH1 = positions[Hbeta1].value_in_unit(unit.nanometer)
        rH2 = positions[Hbeta2].value_in_unit(unit.nanometer)

        v1 = normalize(rCa - rC)
        v2 = normalize(rH1 - rC)
        v3 = normalize(rH2 - rC)

        # Compute plane normal (sp2 plane)
        normal = normalize(np.cross(v2 - v1, v3 - v1))

        # Average substituent direction (points roughly into the sp2 plane)
        avg = normalize(v1 + v2 + v3)

        # The ghost should point mostly along -avg, but tilted out of the plane
        # Combine avg and normal with weights to get ~109.5° separation
        alpha = np.cos(np.deg2rad(109.5))   # ~ -0.333
        beta  = np.sqrt(1 - alpha**2)      # ~ 0.943
        vghost = normalize(alpha * avg + beta * normal)

        anchor = rC + bond_length * vghost
        return anchor
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
    
        # Dertermines if non-adiabatic couplings should be calculated
        if 'NAC' in dynamics_settings:
            self.calc_NAC = dynamics_settings['NAC']
        else:
            self.calc_NAC = False

        if 'ADT' in dynamics_settings:
            self.adiabatic_basis = True

        if 'cluster_run' in dynamics_settings:
            self.cluster_run = dynamics_settings['cluster_run']
        # TODO: these need to be different dictionary entries
        if 'collect_qm_points_from' in dynamics_settings:
            self.collect_qm_points = dynamics_settings['collect_qm_points_from']

        # time step at which QM data point collection should stop
        if 'qmc_stop' in dynamics_settings:
            self.qmc_stop = float(dynamics_settings["qmc_stop"])

        # index of the excited state (in case of TDDFT QM data points)
        if "roots_to_follow" in dynamics_settings:
            self.nstates = len(list(dynamics_settings["roots_to_follow"]))
            self.roots_to_follow = list(dynamics_settings['roots_to_follow'])

        if 'excitation_pulse' in dynamics_settings and dynamics_settings['excitation_pulse'] is not None:
            self.excitation_pulse = list(dynamics_settings['excitation_pulse'])
        
        # keywords used to reverse the dynamics
        # to a previous iteration.
        if 'reverse' in dynamics_settings:
            self.reverse = True
        else:
            self.reverse = False

        # how many iterations to reverse
        if 'reverse_by' in dynamics_settings:
            self.reverse_by = int(dynamics_settings['reverse_by'])
            self.reverse = True
        
        if 'energy_threshold' in dynamics_settings:
            self.energy_threshold = dynamics_settings['energy_threshold']
        
        if 'grad_rmsd_thrsh' in dynamics_settings:
            self.gradient_rmsd_thrsh = dynamics_settings['grad_rmsd_thrsh']
        
        if 'force_orient_thrsh' in dynamics_settings:
            self.force_orient_thrsh = dynamics_settings['force_orient_thrsh']

        if 'sampling_drivers' in dynamics_settings:
            sampling_settings = dynamics_settings.get('sampling_settings', {'enabled':False,
            'e_thrsh_kcal_per_atom': 0.1,
            'g_thrsh_kcal_ang_per_atom':2.0,
            'force_orient_cos': 0.0001})
            self.sampling_enabled = bool(sampling_settings.get('enabled'))
            self.sampling_driver = dynamics_settings.get('sampling_drivers', None)
            self.sampling_settings = sampling_settings
        

        # threshold to determine if a bond is breaking
        if 'bond_threshold' in dynamics_settings:
            self.bond_threshold = float(dynamics_settings['bond_threshold'])

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
                if '_symmetry' not in label:
                    self.sampling_qm_data_point_dict[root].append(dp)
                    old_label = dp.point_label
                    self.sampling_qm_symmetry_datapoint_dict[root][old_label] = [dp]
                else:
                    self.sampling_qm_symmetry_datapoint_dict[root][old_label].append(dp)

            drv.qm_data_points = self.sampling_qm_data_point_dict[root]
            drv.qm_symmetry_data_points = self.sampling_qm_symmetry_datapoint_dict[root]
            drv.prepare_runtime_data_cache(force=True)
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
            if '_symmetry' not in label:
                qm_data_point = InterpolationDatapoint(self.root_z_matrix[root])
                qm_data_point.update_settings(self.interpolation_settings[root])
                qm_data_point.read_hdf5(self.interpolation_settings[root]['imforcefield_file'], label)
                qm_data_point.inv_sqrt_masses = inv_sqrt_masses

                self.qm_data_point_dict[root].append(qm_data_point)
                self.qm_energies_dict[root].append(qm_data_point.energy)
                self.sorted_state_spec_im_labels[root].append(label)

                old_label = qm_data_point.point_label
                self.qm_symmetry_datapoint_dict[root][old_label] = [qm_data_point]
            else:
                sym_dp = InterpolationDatapoint(self.root_z_matrix[root])
                sym_dp.update_settings(self.interpolation_settings[root])
                sym_dp.read_hdf5(self.interpolation_settings[root]['imforcefield_file'], label)
                self.qm_symmetry_datapoint_dict[root][old_label].append(sym_dp)


        driver_object.qm_symmetry_data_points = self.qm_symmetry_datapoint_dict[root]
        driver_object.qm_data_points = self.qm_data_point_dict[root]
        driver_object.labels = self.sorted_state_spec_im_labels[root]

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
                 use_optimized_observables=True):
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

            # Compute anchor from geometry
            if self.bias_force_reaction_prop is not None:
                if self.bias_force_reaction_prop[4]:
                    ghost_index = self.system.addParticle(0.0)
                    self.add_bias_force(
                        (self.bias_force_reaction_prop[0][0], ghost_index),
                        self.bias_force_reaction_prop[1],
                        anchor=self.bias_force_reaction_prop[4])

                    Cbeta, Calpha, H1beta, H2beta = (
                        self.bias_force_reaction_idx[0],
                        self.bias_force_reaction_idx[1],
                        self.bias_force_reaction_idx[2],
                        self.bias_force_reaction_idx[3],
                    )
                    anchor = self.compute_tetrahedral_anchor(
                        self.positions, Cbeta, Calpha, H1beta, H2beta)
                    anchor_pos = np.array(anchor) * unit.nanometer
                    self.positions = list(self.positions) + [anchor_pos]
                else:
                    self.global_theta_list = np.linspace(0, self.bias_force_reaction_prop[2], 15, endpoint=True)
                    self.add_bias_force(
                        self.bias_force_reaction_prop[0],
                        self.bias_force_reaction_prop[1],
                        self.global_theta_list[0])
            
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
            if self.bias_force_reaction_prop is not None and self.bias_force_reaction_prop[4]:
                self.simulation.reporters.append(GhostSafePDBReporter(self.out_file, save_freq))
            else:
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


            if len(self.sampled_molecules[root]['molecules']) != 0 and self.add_gpr_model:


                # Last molecule only
                x_list, y_list = [], []
                x_groups = None
                for idx, mol in enumerate(self.allowed_molecules[root]['molecules'][:]):
                    # 1) ICs + groups (groups should be identical for every mol)
                    X_single, groups, types = driver_object.compute_ic_and_groups(mol.get_coordinates_in_bohr())
                    if x_groups is None:
                        x_groups = groups
                    else:
                        # sanity: same layout across all mols
                        assert all(np.array_equal(x_groups[k], groups[k]) for k in x_groups.keys()), "IC groups changed!"
                    x_list.append(np.asarray(X_single, dtype=np.float64).reshape(1, -1))
                    # 2) energies -> per-atom ΔE (use mol, not molecule!)
                    per_atom = len(mol.get_labels())
                    y_list.append(
                        np.asarray(
                            (self.allowed_molecules[root]['qm_energies'][idx] - self.allowed_molecules[root]['im_energies'][idx])
                            * hartree_in_kcalpermol() / per_atom,
                            dtype=np.float64
                        ).reshape(-1)
                    )
                # 3) stack batch
                X_np = np.vstack(x_list)              # (N,D)
                y_np = np.concatenate(y_list, axis=0) # (N,)
                
                

                gpr_driver = GPRInterpolationDriver(X_np, x_groups, y_np)
                gpr_driver.refit(steps="auto", verbose=False)
                base = gpr_driver.model.covar_module.base_kernel
                groups = base.kernels if hasattr(base, "kernels") else [base]
                ells = [float(k.lengthscale.detach().cpu().view(-1)[0]) for k in groups]
                print("group ℓ:", ells)

                # outputscale & noise
                print("outputscale:", float(gpr_driver.model.covar_module.outputscale.detach().cpu()))
                print("noise:", float(gpr_driver.likelihood.noise.detach().cpu()))
                driver_object.gpr_intdriver = gpr_driver 
            

            self.impes_drivers[root] = driver_object
            self._reload_interpolation_root_from_hdf5(root, inv_sqrt_masses)
            self.prev_dens_of_points[root] = len(self.qm_data_point_dict[root])
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
            if not self._mpi_is_root():
                self._run_qmmm_worker_service()
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
                break
            
            if not self.expansion:
                if len(self.expansion_molecules) >= 20:
                    break

            if self.unadded_cycles == 0:
                break

        if self._mpi_is_active() and self.mpi_root_worker_mode and self._mpi_is_root():
            self._mpi_bcast_control({'cmd': self._MPI_CMD_STOP})

        end_time = time()
        elapsed_time = end_time - start_time
        elapsed_time_days = elapsed_time / (24 * 3600)
        performance = (self.nsteps * timestep / 1e6) / elapsed_time_days
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

    def run_qmmm_legacy(self,
                        verify_against_optimized=False,
                        verification_atol=1.0e-8,
                        verification_rtol=1.0e-6,
                        verification_fail_fast=False):
        """
        Legacy QM/MM integration entry point.
        Uses the reference observable collection path from the original implementation.
        """

        return self.run_qmmm(
            verify_optimized=verify_against_optimized,
            verification_atol=verification_atol,
            verification_rtol=verification_rtol,
            verification_fail_fast=verification_fail_fast,
            use_optimized_observables=False)

    def run_qmmm_and_verify(self,
                            verification_atol=1.0e-8,
                            verification_rtol=1.0e-6,
                            verification_fail_fast=False):
        """
        Convenience wrapper to run QM/MM dynamics and compare optimized vs reference
        per-step observable collection.
        """

        return self.run_qmmm(
            verify_optimized=True,
            verification_atol=verification_atol,
            verification_rtol=verification_rtol,
            verification_fail_fast=verification_fail_fast,
            use_optimized_observables=True)


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
    
    def _update_induction_embedding(self, periodic: bool):
        # --- fetch positions (nm) ---
        state = self.simulation.context.getState(getPositions=True)
        pos = state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)

        # --- get charges for MM atoms (e) once and cache ---
        if not hasattr(self, "_mm_charge_vec"):
            nb = [f for f in self.system.getForces() if isinstance(f, mm.NonbondedForce)][0]
            q = np.zeros(self.system.getNumParticles())
            for i in range(nb.getNumParticles()):
                qi, sigma, eps = nb.getParticleParameters(i)
                q[i] = qi._value   # elementary charge
            self._mm_charge_vec = q

        # --- Coulomb constant in kJ/mol·nm·e^-2 ---
        Kc = 138.935456

        # --- Reaction Field constants if periodic & using your RF coupling ---
        if periodic:
            nb = [f for f in self.system.getForces() if isinstance(f, mm.NonbondedForce)][0]
            rfDielectric = nb.getReactionFieldDielectric()
            rc = float(self.cutoff)  # nm
            krf = (1.0/(rc**3)) * (rfDielectric - 1.0) / (2.0*rfDielectric + 1.0)
            # crf drops out of the field (it’s a constant in potential), see derivation below
        else:
            krf = 0.0

        # --- work arrays ---
        F_im = np.zeros((len(self.qm_atoms), 3))                # forces to add on IM atoms
        F_mm = np.zeros((len(self.mm_subregion), 3))            # forces to add on MM atoms
        U_ind = 0.0

        # --- user-supplied isotropic polarizabilities for IM sites (nm^3 units) ---
        # Suggestion: start with element-wise constants or a single scalar.
        # e.g., self.alpha_iso[a] already defined in nm^3 (1 Å^3 = 1e-3 nm^3)
        alpha = self.alpha_iso   # dict or array indexed by atom id

        # Helper to map mm_subregion index -> row id in F_mm
        mm_index_of = {b: k for k, b in enumerate(self.mm_subregion)}

        # --- loop over IM sites, accumulate field and forces ---
        for ia, a in enumerate(self.qm_atoms):
            ra = np.array(pos[a], dtype=float)
            Ea = np.zeros(3, dtype=float)

            # field from all MM charges
            for b in self.mm_subregion:
                qb = self._mm_charge_vec[b]
                if qb == 0.0:
                    continue
                R = ra - np.array(pos[b], dtype=float)
                r2 = float(np.dot(R, R)) + 1e-24
                r = np.sqrt(r2)
                R3 = r2 * r

                # Electric field (Reaction Field consistent):
                # E = Kc * q * ( R/r^3  -  2*krf * R )
                Ea += Kc * qb * (R / R3 - 2.0 * krf * R)

            # induction energy contribution at site a
            alpha_a = float(alpha[a])  # nm^3
            U_ind += -0.5 * alpha_a * float(np.dot(Ea, Ea))

            # Now build forces from d/dR of U = -1/2 alpha |E|^2
            for b in self.mm_subregion:
                qb = self._mm_charge_vec[b]
                if qb == 0.0:
                    continue
                R = ra - np.array(pos[b], dtype=float)
                r2 = float(np.dot(R, R)) + 1e-24
                r = np.sqrt(r2)
                RR = np.outer(R, R)

                # J(R)·E = ∂E/∂R · E  with Reaction Field:
                # ∂E/∂R = Kc*q * [ (r^2 I - 3 RR^T)/r^5  -  2*krf I ]
                JdotE = Kc * qb * (((r2*np.eye(3) - 3.0*RR) @ Ea) / (r2*r2*r) - 2.0*krf * Ea)

                # forces: F_a += +alpha * J·E ; F_b += -alpha * J·E
                Fa =  alpha_a * JdotE
                Fb = -alpha_a * JdotE

                F_im[ia] += Fa
                F_mm[mm_index_of[b]] += Fb

        # --- push to context (don’t forget to call updateParametersInContext) ---
        for k, a in enumerate(self.qm_atoms):
            self.force_IM_embed.setParticleParameters(k, a, F_im[k])
        for k, b in enumerate(self.mm_subregion):
            self.force_MM_embed.setParticleParameters(k, b, F_mm[k])

        self.force_IM_embed.updateParametersInContext(self.simulation.context)
        self.force_MM_embed.updateParametersInContext(self.simulation.context)

        # store energy for your own total-energy bookkeeping
        self.last_U_induction = U_ind     # kJ/mol
    
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

    # def update_gradient_and_energy(self, new_positions):
    #     """
    #     Updates and returns the gradient and potential energy of the QM region.

    #     :param new_positions:
    #         The new positions of the atoms in the QM region.
    #     :return:
    #         The gradient and potential energy of the QM region.
    #     """
    #     ug_t0 = time()

    #     # self.coordinates_xyz.append(new_positions * 10)
    #     positions_ang = new_positions * 10 

    #     atom_labels = getattr(self, '_topology_atom_labels', None)
    #     if atom_labels is None:
    #         atom_labels = [atom.element.symbol for atom in self.topology.atoms()]
    #         self._topology_atom_labels = atom_labels

    #     new_molecule = None
    #     # Check if there is a QM/MM partition in the system
    #     if self.mm_subregion is not None:
    #         # Create a molecule with a link atom (H)
    #         # The linking atom is defined in self.linking_atoms
    #         # It need to be changed to a H atom at 1.0 angstrom
    #         # The QM region is defined in self.qm_atoms

    #         qm_atom_labels = []
    #         for i in self.qm_atoms:
    #             # Change the linking atom to H
    #             if i in self.linking_atoms:
    #                 qm_atom_labels.append('H')
    #             else:
    #                 qm_atom_labels.append(atom_labels[i])

            
    #         # Change the positions of the linking atoms to 1.0 angstrom
    #         for atom1, atom2 in self.broken_bonds:

    #             current_distance = np.linalg.norm(positions_ang[atom1] - positions_ang[atom2])
    #             # Change the position of the linking atom to 1.0 angstrom
    #             # By construction, the first atom is the linking atom
    #             direction = (positions_ang[atom2] - positions_ang[atom1]) / current_distance
    #             positions_ang[atom1] = positions_ang[atom2] - direction * self.linking_atom_distance

    #         new_molecule = Molecule(qm_atom_labels, positions_ang, units="angstrom")
    #         new_molecule.set_charge(self.molecule.get_charge())
    #         new_molecule.set_multiplicity(self.molecule.get_multiplicity())

    #     else:
    #         # Atom labels for the QM region
    #         qm_atom_labels = [atom_labels[i] for i in self.qm_atoms]
    #         new_molecule = Molecule(qm_atom_labels, positions_ang, units="angstrom")
    #         new_molecule.set_charge(self.molecule.get_charge())
    #         new_molecule.set_multiplicity(self.molecule.get_multiplicity())
    #         self.unique_molecules.append(new_molecule)
        
    #     for root in self.roots_to_follow:
    #         driver = self.impes_drivers[root]
    #         if driver.qm_data_points is not self.qm_data_point_dict[root]:
    #             driver.qm_data_points = self.qm_data_point_dict[root]
    #         driver.prepare_runtime_data_cache()
    #         compute_t0 = time()
    #         driver.compute(new_molecule)
    #         compute_dt = time() - compute_t0
    #         self._add_runtime_timing(
    #             f'update_gradient_and_energy.impes_compute_root_{root}',
    #             compute_dt)
    #         self._add_runtime_timing(
    #             'update_gradient_and_energy.impes_compute_total',
    #             compute_dt)

    #     if self.nstates > 1:
    #         for root_1 in range(0, self.nstates):
    #             for root_2 in range(root_1 + 1, self.nstates):

    #                 potential_kjmol = self.impes_drivers[self.roots_to_follow[root_1]].impes_coordinate.energy * hartree_in_kcalpermol() * 4.184
    #                 potential_kjmol_2 = self.impes_drivers[self.roots_to_follow[root_2]].impes_coordinate.energy * hartree_in_kcalpermol() * 4.184
    #                 self.transition_information.append((self.roots_to_follow[root_1], self.roots_to_follow[root_2], potential_kjmol - potential_kjmol_2))
    #                 print(f'compare the energies between roots: {self.roots_to_follow[root_1]} -> {self.roots_to_follow[root_2]}', potential_kjmol_2 - potential_kjmol)

    #                 if self.calc_NAC == True and np.linalg.norm(self.velocities_np[-1]) > 0.0:
    #                     current_NAC = self.impes_drivers[self.roots_to_follow[root_1]].impes_coordinate.NAC.flatten()
    #                     current_velocitites = self.velocities_np[-1].flatten() * 4.566180e-4

    #                     hopping_potential = np.exp(-abs((np.pi/(4)) * (( self.impes_drivers[self.roots_to_follow[root_2]].impes_coordinate.energy - self.impes_drivers[self.roots_to_follow[root_1]].impes_coordinate.energy ) / np.linalg.multi_dot([current_NAC, current_velocitites]))))
    #                     print('#######################', '\n\n', hopping_potential, potential_kjmol_2 - potential_kjmol, '\n\n', '#######################')

    #                 if potential_kjmol_2 - potential_kjmol < 5 and self.swap_back:
    #                     # Choose a random integer between 0 and 1
    #                     random_integer = random.randint(0, 1)
    #                     self.swap_back = False
    #                     random_integer = 1.0
    #                     if random_integer == 1 and self.current_state == self.roots_to_follow[root_1]:
    #                         self.current_state = self.roots_to_follow[root_2]
    #                     elif random_integer == 1 and self.current_state == self.roots_to_follow[root_2]:
    #                         self.current_state = self.roots_to_follow[root_1]
    #                     break

        
    #     else:
    #         self.current_state = self.roots_to_follow[0]

    #     if self.excitation_pulse is not None and len(self.roots_to_follow) > 1 and self.point_checker == self.excitation_pulse[0] and self.current_state != self.excitation_pulse[1] and self.swap_back:
    #         self.current_state = self.excitation_pulse[1]
            

    #     potential_kjmol = self.impes_drivers[self.current_state].impes_coordinate.energy * hartree_in_kcalpermol() * 4.184
    #     self.current_gradient = self.impes_drivers[self.current_state].impes_coordinate.gradient

    #     self.current_energy = potential_kjmol
                        
    #         # Potential energy is in Hartree, convert to kJ/mol
            
    #     self._add_runtime_timing('update_gradient_and_energy.total', time() - ug_t0)

    #     return self.current_gradient, potential_kjmol
    
    def update_gradient_and_energy(self, new_positions, collective_mode=None):
        ug_t0 = time()
        new_molecule = self._build_qm_molecule_from_positions(new_positions)
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
        if self.nstates > 1:
            for root_1 in range(0, self.nstates):
                for root_2 in range(root_1 + 1, self.nstates):

                    potential_kjmol = self.impes_drivers[self.roots_to_follow[root_1]].impes_coordinate.energy * hartree_in_kcalpermol() * 4.184
                    potential_kjmol_2 = self.impes_drivers[self.roots_to_follow[root_2]].impes_coordinate.energy * hartree_in_kcalpermol() * 4.184
                    self.transition_information.append((self.roots_to_follow[root_1], self.roots_to_follow[root_2], potential_kjmol - potential_kjmol_2))
                    print(f'compare the energies between roots: {self.roots_to_follow[root_1]} -> {self.roots_to_follow[root_2]}', potential_kjmol_2 - potential_kjmol)

                    if self.calc_NAC == True and np.linalg.norm(self.velocities_np[-1]) > 0.0:
                        current_NAC = self.impes_drivers[self.roots_to_follow[root_1]].impes_coordinate.NAC.flatten()
                        current_velocitites = self.velocities_np[-1].flatten() * 4.566180e-4

                        hopping_potential = np.exp(-abs((np.pi/(4)) * (( self.impes_drivers[self.roots_to_follow[root_2]].impes_coordinate.energy - self.impes_drivers[self.roots_to_follow[root_1]].impes_coordinate.energy ) / np.linalg.multi_dot([current_NAC, current_velocitites]))))
                        print('#######################', '\n\n', hopping_potential, potential_kjmol_2 - potential_kjmol, '\n\n', '#######################')

                    if potential_kjmol_2 - potential_kjmol < 5 and self.swap_back:
                        # Choose a random integer between 0 and 1
                        random_integer = random.randint(0, 1)
                        self.swap_back = False
                        random_integer = 1.0
                        if random_integer == 1 and self.current_state == self.roots_to_follow[root_1]:
                            self.current_state = self.roots_to_follow[root_2]
                        elif random_integer == 1 and self.current_state == self.roots_to_follow[root_2]:
                            self.current_state = self.roots_to_follow[root_1]
                        break

        
        else:
            self.current_state = self.roots_to_follow[0]

        if self.excitation_pulse is not None and len(self.roots_to_follow) > 1 and self.point_checker == self.excitation_pulse[0] and self.current_state != self.excitation_pulse[1] and self.swap_back:
            self.current_state = self.excitation_pulse[1]


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
        atom_labels = getattr(self, '_topology_atom_labels', None)
        if atom_labels is None:
            atom_labels = [atom.element.symbol for atom in self.topology.atoms()]
            self._topology_atom_labels = atom_labels

        qm_force = getattr(self, '_qm_external_force', None)
        if qm_force is None:
            qm_force = self.system.getForce(self.qm_force_index)
            self._qm_external_force = qm_force

        new_positions = context.getState(getPositions=True).getPositions()
        
        # Update the forces of the QM region
        qm_positions = np.array([new_positions[i].value_in_unit(unit.nanometer) for i in self.qm_atoms])

        self.velocities_np.append(context.getState(getVelocities=True).getVelocities(True))
        gradient_t0 = time()
        gradient = self.update_gradient(qm_positions)
        self._add_runtime_timing('update_forces.update_gradient', time() - gradient_t0)

        positions_ang = (qm_positions) * 10

        qm_atom_labels = [atom_labels[i] for i in self.qm_atoms]

        new_molecule = Molecule(qm_atom_labels, positions_ang, units="angstrom")
        new_molecule.set_charge(self.molecule.get_charge())
        new_molecule.set_multiplicity(self.molecule.get_multiplicity())
        
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
            self.gradients.append(gradient)
            
            if self.skipping_value == 0 and self.point_checker + 1 > self.collect_qm_points:
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
                # else:
                #     self.add_a_point = True

                if self.impes_drivers[self.current_state].gpr_intdriver is not None and len(self.allowed_molecules[self.current_state]['molecules']) > 5:
                
                    X_single_list, groups, type = self.impes_drivers[self.current_state].compute_ic_and_groups(new_molecule.get_coordinates_in_bohr())
                    
                    X_new = np.asarray(
                                X_single_list,
                                dtype=np.float64
                            ).reshape(1, -1)
                    
                    mean, std = self.impes_drivers[self.current_state].gpr_intdriver.predict(X_new, return_std=True) 
                    mean_latent, std_latent = self.impes_drivers[self.current_state].gpr_intdriver.predict_latent(X_new, return_std=True) 

                    ### second check

                    trigger, information = self.impes_drivers[self.current_state].gpr_intdriver.should_qm_check(X_new, T=self.energy_threshold)
                    
                    error_p_variance = abs(mean) + std
                    
                    error_proxy = float(np.abs(mean).ravel()[0] + std.ravel()[0])
                    # print('Error proxy', error_proxy, trigger, information)
                    Xt, yt = self.impes_drivers[self.current_state].gpr_intdriver._transform(self.impes_drivers[self.current_state].gpr_intdriver.X_raw, self.impes_drivers[self.current_state].gpr_intdriver.y_raw)
                      
                    print('Here is the information', information, trigger)

                    self.add_a_point = False
                    if trigger:
                        self.add_a_point = True 
                    self.prev_dE_gpr = error_p_variance


            else:
                self.skipping_value -= 1
                if self.skipping_value < 0:
                    self.skipping_value = 0
            
            if self.step % 50 == 0 and self.step > self.collect_qm_points:
                self.add_a_point = True
            
            # if self.current_state != self.prev_state:
            #     self.add_a_point = True

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
                
                new_positions = context.getState(getPositions=True).getPositions()
                qm_positions = np.array([new_positions[i].value_in_unit(unit.nanometer) for i in self.qm_atoms])
                positions_ang = (qm_positions) * 10
                qm_atom_labels = [atom_labels[i] for i in self.qm_atoms]
                new_molecule = Molecule(qm_atom_labels, positions_ang, units="angstrom")
                new_molecule.set_charge(self.molecule.get_charge())
                new_molecule.set_multiplicity(self.molecule.get_multiplicity())
                gradient_2 = self.update_gradient(qm_positions, collective_mode=False)
                
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
            new_positions = context.getState(getPositions=True).getPositions()
            qm_positions = np.array([new_positions[i].value_in_unit(unit.nanometer) for i in self.qm_atoms])
            positions_ang = (qm_positions) * 10
            qm_atom_labels = [atom_labels[i] for i in self.qm_atoms]
            new_molecule = Molecule(qm_atom_labels, positions_ang, units="angstrom")
            new_molecule.set_charge(self.molecule.get_charge())
            new_molecule.set_multiplicity(self.molecule.get_multiplicity())
            gradient_2 = self.update_gradient(qm_positions, collective_mode=False)
            force = -np.array(gradient_2) * conversion_factor
            
            self.state_specific_molecules[self.current_state].append(new_molecule)
            self.coordinates_xyz.append(qm_positions * 10)
            for i, atom_idx in enumerate(self.qm_atoms):
                qm_force.setParticleParameters(i, atom_idx, force[i])
            qm_force.updateParametersInContext(context)

        if self.current_state < self.prev_state and 1 == 2:
            
            qm = self.get_qm_potential_energy()
            self.qm_potentials.append(qm)

            # QM/MM interactions
            qm_mm = context.getState(getEnergy=True, groups={1,2}).getPotentialEnergy()

            # MM region 
            mm = context.getState(getEnergy=True, groups={3,4,5,6,7}).getPotentialEnergy()

            # Total potential energy
            pot = qm * unit.kilojoules_per_mole + qm_mm + mm

            # Kinetic energy
            kinetic = context.getState(getEnergy=True).getKineticEnergy()
    
            # Total energy
            total = pot + kinetic
            E_tot_diff = np.sqrt(abs(self.total_energies[-1] - total.value_in_unit(unit.kilojoules_per_mole)) / kinetic.value_in_unit(unit.kilojoules_per_mole))
    
            # Assume you have a Simulation object called 'simulation'
            state = context.getState(getVelocities=True)
            velocities = state.getVelocities(asNumpy=True)  # unit: nanometers/picosecond

            scaling_factor = E_tot_diff  # for example, to reduce velocities by 20%
            velocities *= scaling_factor

            context.setVelocities(velocities)
        
        if len(self.roots_to_follow) > 1 and self.swap_back == False:# and self.temperatures[-1] > 50:
            state = context.getState(getVelocities=True)
            velocities = state.getVelocities(asNumpy=True)  # unit: nanometers/picosecond

         
            velocities *= 0.01

            context.setVelocities(velocities)
            print('I am setting back the velocities')

        elif len(self.roots_to_follow) > 1 and self.swap_back == False and self.temperatures[-1] < self.temperatures[0]:
            self.swap_back = True
            self.point_checker = 0
        
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
            print(e_diff)
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
        
        screen = self._sampling_screen_check(molecule)
        if screen["skip_full_qm"]:
            # fast accept: no point addition, keep allowed_molecules bookkeeping
            # for root in self.roots_to_follow:
            #     self.allowed_molecules[root]['molecules'].append(molecule)
            #     self.allowed_molecules[root]['im_energies'].append(self.impes_drivers[root].impes_coordinate.energy)
            #     self.allowed_molecules[root]['qm_energies'].append(self.impes_drivers[root].impes_coordinate.energy)
            #     self.allowed_molecules[root]['qm_gradients'].append(self.impes_drivers[root].impes_coordinate.gradient)
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
            print('############# Energy is QM claculated ############')
            print('identifaction state', identification_state, 'current state', self.current_state)
            current_basis = None
            if identification_state == 0:
                current_basis = MolecularBasis.read(molecule, self.basis_set_label['gs'])
            else:
                current_basis = MolecularBasis.read(molecule, self.basis_set_label['es'])
            
   
            qm_energy, scf_results, rsp_results = self._compute_energy_mpi_safe(drivers[0], molecule, current_basis)
            
            if isinstance(drivers[0], LinearResponseEigenSolver) or isinstance(drivers[0], TdaEigenSolver):
                qm_energy = qm_energy[mask]


            gradients = self._compute_gradient_mpi_safe(drivers[1], molecule, current_basis, scf_results, rsp_results)

            for e_idx in range(len(qm_energy)):
                masses = molecule.get_masses().copy()
                masses_cart = np.repeat(masses, 3)
                inv_sqrt_masses = 1.0 / np.sqrt(masses_cart)
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

        if current_state_difference[self.current_state][0] > self.energy_threshold or current_state_difference[self.current_state][1] > self.gradient_rmsd_thrsh or current_state_difference[self.current_state][2] < self.force_orient_thrsh:
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

            for state in self.roots_to_follow:
                if self.use_opt_confidence_radius[0] and len(self.allowed_molecules[state]['molecules']) > 20 and 1 == 2:
                    self.add_a_point = True
                    if not self.expansion:
                        length_vectors = (self.impes_drivers[-1].impes_coordinate.cartesian_distance_vector(self.qm_data_points[0]))
                        rmsd = (np.linalg.norm(length_vectors) / np.sqrt(len(self.molecule.get_labels())) * bohr_in_angstrom())
                        self.expansion_molecules.append((molecule, energy_difference / natms * hartree_in_kcalpermol(), rmsd, (np.linalg.norm(length_vectors))))
                        self.last_point_added = self.point_checker - 1
                        self.point_checker = 0
                
                else:
                    self.add_a_point = True
                    #TODO: Check if GPR implementation is correctly implemented and updates correctly considering the GPRinterpolationdriver class
                    if self.add_gpr_model:
                        print('GPR check is performed')
                        for root in self.roots_to_follow:

                            if len(self.allowed_molecules[state]['molecules']) < 2:
                                continue

                            x_list, y_list = [], []
                            interpolation_driver = InterpolationDriver(self.root_z_matrix[root])
                            interpolation_driver.update_settings(self.interpolation_settings[root])
                            interpolation_driver.symmetry_information = self.impes_drivers[root].symmetry_information
                            interpolation_driver.qm_symmetry_data_points = self.impes_drivers[root].qm_symmetry_data_points
                            interpolation_driver.impes_coordinate.inv_sqrt_masses = self.inv_sqrt_masses
                            interpolation_driver.distance_thrsh = 1000
                            interpolation_driver.use_symmetry = self.use_symmetry
                            interpolation_driver.exponent_p = self.impes_drivers[root].exponent_p
                            interpolation_driver.print = False
                            interpolation_driver.qm_data_points = self.impes_drivers[root].qm_data_points
                            interpolation_driver.calc_optim_trust_radius = True

                            x_groups = None

                            for idx, mol in enumerate(self.allowed_molecules[root]['molecules']):
                                X_single, groups, types = self.impes_drivers[root].compute_ic_and_groups(mol.get_coordinates_in_bohr())
                                if x_groups is None:
                                    x_groups = groups
                                else:                           
                                    assert all(np.array_equal(x_groups[k], groups[k]) for k in x_groups.keys()), "IC groups changed!"

                                x_list.append(np.asarray(X_single, dtype=np.float64).reshape(1, -1))

                 
                                interpolation_driver.compute(mol)
                                im_energy = interpolation_driver.get_energy()
                                self.allowed_molecules[root]['im_energies'][idx] = im_energy

                                per_atom = len(mol.get_labels())
                                y_list.append(
                                    np.asarray(
                                        (self.allowed_molecules[root]['qm_energies'][idx] - im_energy)
                                        * hartree_in_kcalpermol() / per_atom,
                                        dtype=np.float64
                                    ).reshape(-1)
                                )

   
                            X_np = np.vstack(x_list)           
                            y_np = np.concatenate(y_list, axis=0) 
                            new_idx = X_np.shape[0] - 1

                            drv = self.impes_drivers[root].gpr_intdriver
                            
                            if drv is not None and len(self.allowed_molecules[root]['molecules']) > 2:
                                X_t = torch.from_numpy(X_np).to(drv.dtype).to(drv.device)
                                y_t = torch.from_numpy(y_np).to(drv.dtype).to(drv.device)
                                # Current training set (exclude the new point)
                                Xtr = X_t[:-1, :]
                                ytr = y_t[:-1]
                                Xq  = X_t[new_idx:new_idx+1, :]    # (1,D)
                                y_true = float(y_t[new_idx].cpu().numpy())
                                # Transform with the CURRENT drv scalers (not refit yet)
                                Xtr_std, ytr_std = drv._transform(Xtr, ytr)
                                Xq_std           = drv._transform(Xq)
                                # Latent prediction out-of-sample (using current drv trained on X_raw/y_raw)
                                mu_lat, sig_lat = drv.predict_latent(Xq, return_std=True)  # add this method if you haven’t
                                mu_lat = float(mu_lat.squeeze()); sig_lat = float(sig_lat.squeeze())
                                # Effective distance of query to current training set
                                d_norm = drv._compute_d_norm(Xq_std, Xtr_std, groups)  # ← normalized distance
                                d_thresh = drv.get_distance_threshold(1.0)
                                s_max       = drv.max_kernel_similarity(Xq_std, Xtr_std, topk=5)
                                # Scalar absolute residual for THIS query (QM label already available)
                                abs_res = abs(y_true - mu_lat)

                                error_good = (abs_res <= 0.1 * self.energy_threshold)
                                necessary = (not error_good)
                                
                                # Update the dynamic threshold with this single pair
                                drv.record_distance_event(d_norm, necessary, abs_res, T=self.energy_threshold)
                                drv.update_distance_threshold()
                                
                                drv.record_similarity_event(s_max, necessary, abs_res, T=self.energy_threshold)
                                drv.update_similarity_threshold()
                                # ---------- 2) NOW add the point and refit (warm start) ----------
                                drv.replace_data(X_t, y_t, clone=False)
                                drv.refit(steps="auto", verbose=False)
                                # (optional) log query-centric diagnostics after refit
                                print(f"[update] d_norm={d_norm:.3f}, |residual|={abs_res:.4g}, "
                                    f"mu_lat={mu_lat:.4g}, sigma_lat={sig_lat:.4g}, T={self.energy_threshold}")
                                
                            # 4) predict for the last structure (shape (1,D))
                            mean, std = self.impes_drivers[root].gpr_intdriver.predict(X_np[-1:,:], return_std=True)
                            print('### Predictions:', mean.ravel().tolist(), std.ravel().tolist())
                            # 5) (optional) inspect kernel hyperparameters to judge reasonableness
                                              
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

                if self.add_gpr_model:
                    # Last molecule only
                    x_list, y_list = [], []
                    interpolation_driver = InterpolationDriver(self.root_z_matrix[root])
                    interpolation_driver.update_settings(self.interpolation_settings[root])
                    interpolation_driver.symmetry_information = self.impes_drivers[root].symmetry_information
                    interpolation_driver.qm_symmetry_data_points = self.impes_drivers[root].qm_symmetry_data_points
                    interpolation_driver.impes_coordinate.inv_sqrt_masses = self.inv_sqrt_masses
                    interpolation_driver.distance_thrsh = 1000
                    interpolation_driver.eq_bond_force_constants = self.eq_bond_force_constants
                    interpolation_driver.impes_coordinate.eq_bond_force_constants = self.eq_bond_force_constants
                    interpolation_driver.exponent_p = self.impes_drivers[root].exponent_p
                    interpolation_driver.print = False
                    interpolation_driver.qm_data_points = self.impes_drivers[root].qm_data_points
                    interpolation_driver.calc_optim_trust_radius = False
                    interpolation_driver = self.use_symmetry
                    x_groups = None
                    for idx, mol in enumerate(self.allowed_molecules[root]['molecules']):
                        # 1) ICs + groups (groups should be identical for every mol)
                        X_single, groups, types = self.impes_drivers[root].compute_ic_and_groups(mol.get_coordinates_in_bohr())
                        if x_groups is None:
                            x_groups = groups
                        else:
                            # sanity: same layout across all mols
                            assert all(np.array_equal(x_groups[k], groups[k]) for k in x_groups.keys()), "IC groups changed!"
                        x_list.append(np.asarray(X_single, dtype=np.float64).reshape(1, -1))
                        # 2) energies -> per-atom ΔE (use mol, not molecule!)
                        interpolation_driver.compute(mol)
                        im_energy = interpolation_driver.get_energy()
                        self.allowed_molecules[root]['im_energies'][idx] = im_energy
                        per_atom = len(mol.get_labels())
                        y_list.append(
                            np.asarray(
                                (self.allowed_molecules[root]['qm_energies'][idx] - im_energy)
                                * hartree_in_kcalpermol() / per_atom,
                                dtype=np.float64
                            ).reshape(-1)
                        )
                    # 3) stack batch
                    X_np = np.vstack(x_list)              # (N,D)
                    y_np = np.concatenate(y_list, axis=0) # (N,)
                    new_idx = X_np.shape[0] - 1
                    drv = self.impes_drivers[root].gpr_intdriver
                    if drv is None:
                        # First time: create once, with grouped kernel
    
                        gpr_driver = GPRInterpolationDriver(X_np, x_groups, y_np)  # your class builds grouped kernel inside __init__
                        # Fit hypers at least once so lengthscales/noise move off defaults
                        gpr_driver.refit()
                        self.impes_drivers[root].gpr_intdriver = gpr_driver
                        
                    else:
                        X_t = torch.from_numpy(X_np).to(drv.dtype).to(drv.device)
                        y_t = torch.from_numpy(y_np).to(drv.dtype).to(drv.device)
                        # Current training set (exclude the new point)
                        Xtr = X_t[:-1, :]
                        ytr = y_t[:-1]
                        Xq  = X_t[new_idx:new_idx+1, :]    # (1,D)
                        y_true = float(y_t[new_idx].cpu().numpy())
                        # Transform with the CURRENT drv scalers (not refit yet)
                        Xtr_std, ytr_std = drv._transform(Xtr, ytr)
                        Xq_std           = drv._transform(Xq)
                        # Latent prediction out-of-sample (using current drv trained on X_raw/y_raw)
                        mu_lat, sig_lat = drv.predict_latent(Xq, return_std=True)  # add this method if you haven’t
                        mu_lat = float(mu_lat.squeeze()); sig_lat = float(sig_lat.squeeze())
                        # Effective distance of query to current training set
                        d_norm = drv._compute_d_norm(Xq_std, Xtr_std, groups)  # ← normalized distance
                        d_thresh = drv.get_distance_threshold(1.0)
                        s_max       = drv.max_kernel_similarity(Xq_std, Xtr_std, topk=5)
                        # Scalar absolute residual for THIS query (QM label already available)
                        abs_res = abs(y_true - mu_lat)

                        error_good = (abs_res <= 0.1 * self.energy_threshold)
                        necessary = (not error_good)


                        
                        # Update the dynamic threshold with this single pair
                        drv.record_distance_event(d_norm, necessary, abs_res, T=self.energy_threshold)
                        drv.update_distance_threshold()
                        
                        drv.record_similarity_event(s_max, necessary, abs_res, T=self.energy_threshold)
                        drv.update_similarity_threshold()
                        # ---------- 2) NOW add the point and refit (warm start) ----------
                        drv.replace_data(X_t, y_t, clone=False)
                        drv.refit(steps="auto", verbose=False)
                        # (optional) log query-centric diagnostics after refit
                        print(f"[update] d_norm={d_norm:.3f}, |residual|={abs_res:.4g}, "
                            f"mu_lat={mu_lat:.4g}, sigma_lat={sig_lat:.4g}, T={self.energy_threshold}")
                        
                    # 4) predict for the last structure (shape (1,D))
                    mean, std = self.impes_drivers[root].gpr_intdriver.predict(X_np[-1:,:], return_std=True)
                    print('### Predictions:', mean.ravel().tolist(), std.ravel().tolist())
                    # 5) (optional) inspect kernel hyperparameters to judge reasonableness

                self.write_qm_energy_determined_points_h5(self.allowed_molecules[root]['molecules'][self.last_added: ],
                                                    self.allowed_molecules[root]['qm_energies'][self.last_added:],
                                                    self.allowed_molecules[root]['qm_gradients'][self.last_added:],
                                                    self.allowed_molecules[root]['im_energies'][self.last_added:],
                                                    root)
                self.last_added = len(self.allowed_molecules[root]['molecules'])

            print(self.use_opt_confidence_radius[0], len(self.allowed_molecules[self.current_state]['molecules']) >= 10, self.density_around_data_point[0][self.current_state] > 1, self.density_around_data_point[0][self.current_state] % 1 == 0, self.prev_dens_of_points[self.current_state] != self.density_around_data_point[0][self.current_state])
            
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


                self.add_a_point = False
        
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
            print('✨ A point is added! ✨', self.point_checker)
            print(molecule.get_xyz_string())
            ############# Implement constraint optimization ############
            
            state_specific_molecules = []
            
            imp_int_coord = None
            opt_results = None
            
            if self.identfy_relevant_int_coordinates[0]:  
                for state_to_optim in addition_of_state_specific_points:                  
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
                    primary_constraint, backup_constraint, ranked_block = self.impes_drivers[state_to_optim].determine_important_internal_coordinates(state_specific_energies[state_to_optim][0], state_specific_gradients[state_to_optim][0], molecule, self.root_z_matrix[state_to_optim], internal_coordinate_datapoints)
                    corr_new_dp_distrib = False
                    # primary_constraint.append([2,3])
                    # primary_constraint = [[2,3]]
                    main_constraint_list = backup_constraint
                    old_const_len = len(main_constraint_list)
                    opt_results = None
                    while not corr_new_dp_distrib:
                        print('CONSTRAINTS', primary_constraint, backup_constraint, ranked_block)
                        if state_to_optim == 0 and isinstance(self.drivers['gs'][0], ScfRestrictedDriver) or state_to_optim == 0 and isinstance(self.drivers['gs'][0], ScfUnrestrictedDriver):
                            _, scf_tensors, rsp_results = self._compute_energy_mpi_safe(self.drivers['gs'][0], molecule, current_basis)
        
                            opt_drv = OptimizationDriver(self.drivers['gs'][0])
                            opt_drv.ostream.mute()
                            
                            opt_constraint_list = []
                            for constraint in main_constraint_list:
                                if len(constraint) == 2:
                                    opt_constraint = f"freeze distance {constraint[0] + 1} {constraint[1] + 1}"
                                    opt_constraint_list.append(opt_constraint)
                                
                                elif len(constraint) == 3:
                                    opt_constraint = f"freeze angle {constraint[0] + 1} {constraint[1] + 1} {constraint[2] + 1}"
                                    opt_constraint_list.append(opt_constraint)
                            
                                else:
                                    opt_constraint = f"freeze dihedral {constraint[0] + 1} {constraint[1] + 1} {constraint[2] + 1} {constraint[3] + 1}"
                                    opt_constraint_list.append(opt_constraint)
                            
                            for constraint in self.identfy_relevant_int_coordinates[1]:
                                if constraint in opt_constraint_list:
                                    continue
                                if len(constraint) == 2:
                                    opt_constraint = f"freeze distance {constraint[0] + 1} {constraint[1] + 1}"
                                    opt_constraint_list.append(opt_constraint)
                                
                                elif len(constraint) == 3:
                                    opt_constraint = f"freeze angle {constraint[0] + 1} {constraint[1] + 1} {constraint[2] + 1}"
                                    opt_constraint_list.append(opt_constraint)
                            
                                else:
                                    opt_constraint = f"freeze dihedral {constraint[0] + 1} {constraint[1] + 1} {constraint[2] + 1} {constraint[3] + 1}"
                                    opt_constraint_list.append(opt_constraint)
                            
                            opt_drv.constraints = opt_constraint_list
                            opt_results = self._opt_compute_mpi_safe(opt_drv, molecule, current_basis, scf_tensors)
                            optimized_molecule = Molecule.from_xyz_string(opt_results['final_geometry'])
                            optimized_molecule.set_charge(molecule.get_charge())
                            optimized_molecule.set_multiplicity(molecule.get_multiplicity())
                            current_basis = MolecularBasis.read(optimized_molecule, current_basis.get_main_basis_label())

                            # qm_energy, scf_tensors = self.compute_energy(drivers[0], optimized_molecule, opt_current_basis)
                            print('Optimized Molecule', optimized_molecule.get_xyz_string(), '\n\n', molecule.get_xyz_string())
                        
                        elif state_to_optim == 0 and isinstance(self.drivers['gs'][0], XtbDriver):

                            opt_drv = OptimizationDriver(self.drivers['gs'][0])
                            opt_drv.ostream.mute()
                            
                            opt_constraint_list = []
                            for constraint in main_constraint_list:
                                if len(constraint) == 2:
                                    opt_constraint = f"freeze distance {constraint[0] + 1} {constraint[1] + 1}"
                                    opt_constraint_list.append(opt_constraint)
                                
                                elif len(constraint) == 3:
                                    opt_constraint = f"freeze angle {constraint[0] + 1} {constraint[1] + 1} {constraint[2] + 1}"
                                    opt_constraint_list.append(opt_constraint)
                            
                                else:
                                    opt_constraint = f"freeze dihedral {constraint[0] + 1} {constraint[1] + 1} {constraint[2] + 1} {constraint[3] + 1}"
                                    opt_constraint_list.append(opt_constraint)
                            
                            for constraint in self.identfy_relevant_int_coordinates[1]:
                                if constraint in opt_constraint_list:
                                    continue
                                if len(constraint) == 2:
                                    opt_constraint = f"freeze distance {constraint[0] + 1} {constraint[1] + 1}"
                                    opt_constraint_list.append(opt_constraint)
                                
                                elif len(constraint) == 3:
                                    opt_constraint = f"freeze angle {constraint[0] + 1} {constraint[1] + 1} {constraint[2] + 1}"
                                    opt_constraint_list.append(opt_constraint)
                            
                                else:
                                    opt_constraint = f"freeze dihedral {constraint[0] + 1} {constraint[1] + 1} {constraint[2] + 1} {constraint[3] + 1}"
                                    opt_constraint_list.append(opt_constraint)
                            opt_drv.constraints = opt_constraint_list
                            
                            opt_results = self._opt_compute_mpi_safe(opt_drv, molecule)
                            optimized_molecule = Molecule.from_xyz_string(opt_results['final_geometry'])
                            optimized_molecule.set_charge(molecule.get_charge())
                            optimized_molecule.set_multiplicity(molecule.get_multiplicity())

                        
                        elif state_to_optim >= 1 and isinstance(self.drivers['es'][0], LinearResponseEigenSolver) or state_to_optim >= 1 and isinstance(self.drivers['es'][0], TdaEigenSolver):
                                
                            self.drivers['es'][1].state_deriv_index = [state_to_optim]
                            opt_drv = OptimizationDriver(self.drivers['es'][1])
                            _, _, rsp_results = self._compute_energy_mpi_safe(self.drivers['es'][0], molecule, current_basis)
                            opt_drv.ostream.mute()
                            
                            opt_constraint_list = []

                            for constraint in main_constraint_list:
                                if len(constraint) == 2:
                                    opt_constraint = f"freeze distance {constraint[0] + 1} {constraint[1] + 1}"
                                    opt_constraint_list.append(opt_constraint)
                                
                                elif len(constraint) == 3:
                                    opt_constraint = f"freeze angle {constraint[0] + 1} {constraint[1] + 1} {constraint[2] + 1}"
                                    opt_constraint_list.append(opt_constraint)
                            
                                else:
                                    opt_constraint = f"freeze dihedral {constraint[0] + 1} {constraint[1] + 1} {constraint[2] + 1} {constraint[3] + 1}"
                                    opt_constraint_list.append(opt_constraint)
                            
                            for constraint in self.identfy_relevant_int_coordinates[1]:
                                if constraint in opt_constraint_list:
                                    continue
                                if len(constraint) == 2:
                                    opt_constraint = f"freeze distance {constraint[0] + 1} {constraint[1] + 1}"
                                    opt_constraint_list.append(opt_constraint)
                                
                                elif len(constraint) == 3:
                                    opt_constraint = f"freeze angle {constraint[0] + 1} {constraint[1] + 1} {constraint[2] + 1}"
                                    opt_constraint_list.append(opt_constraint)
                            
                                else:
                                    opt_constraint = f"freeze dihedral {constraint[0] + 1} {constraint[1] + 1} {constraint[2] + 1} {constraint[3] + 1}"
                                    opt_constraint_list.append(opt_constraint)
                            opt_drv.constraints = opt_constraint_list
                            self.drivers['es'][3].ostream.mute()
                            self.drivers['es'][0].ostream.mute()
                            opt_results = self._opt_compute_mpi_safe(opt_drv, molecule, current_basis, self.drivers['es'][3], self.drivers['es'][0], rsp_results)
                            excitated_roots = [root for root in self.roots_to_follow if root != 0]
                            self.drivers['es'][1].state_deriv_index = excitated_roots
                            optimized_molecule = Molecule.from_xyz_string(opt_results['final_geometry'])
                            optimized_molecule.set_charge(molecule.get_charge())
                            optimized_molecule.set_multiplicity(molecule.get_multiplicity())
                            print('Optimized Molecule', optimized_molecule.get_xyz_string(), '\n\n', molecule.get_xyz_string()) 
                            current_basis = MolecularBasis.read(optimized_molecule, current_basis.get_main_basis_label())       
                        
                        imp_int_coord = main_constraint_list

                        # same = False
                        # sum_of_values = 0
                        # old_sorted_shorted = []
                        # for item, vlaue in sorted_items:
                        #     sum_of_values += vlaue
                        #     old_sorted_shorted.append(item)
                        #     if sum_of_values > 0.8:       
                        #         break

                        # while not same:
                            
                        #     interpolation_driver = InterpolationDriver(self.root_z_matrix[state_to_optim])
                        #     interpolation_driver.update_settings(self.interpolation_settings[state_to_optim])
                        #     interpolation_driver.symmetry_information = self.impes_drivers[state_to_optim].symmetry_information
                        #     interpolation_driver.qm_symmetry_data_points = self.impes_drivers[state_to_optim].qm_symmetry_data_points
                        #     interpolation_driver.impes_coordinate.inv_sqrt_masses = self.inv_sqrt_masses
                        #     interpolation_driver.distance_thrsh = 1000
                        #     interpolation_driver.exponent_p = self.impes_drivers[state_to_optim].exponent_p
                        #     interpolation_driver.print = False
                        #     interpolation_driver.use_symmetry = self.use_symmetry
                        #     interpolation_driver.qm_data_points = self.impes_drivers[state_to_optim].qm_data_points
                        #     interpolation_driver.calc_optim_trust_radius = True

                        #     interpolation_driver.compute(optimized_molecule)

                        #     print('enegy difference', opt_results['opt_energies'][-1], interpolation_driver.get_energy(), abs(opt_results['opt_energies'][-1] - interpolation_driver.get_energy()) * hartree_in_kcalpermol())

                            
                            # _, distance_to_org, _ = self.calculate_distance_to_ref(optimized_molecule.get_coordinates_in_bohr(), 
                            #                                                        molecule.get_coordinates_in_bohr(), self.impes_drivers[state_to_optim].symmetry_information)
                            # print(distance_to_org)
                            # # distances_dp_to_org = []
                            # # distances_dp_to_opt = []
                            # for dp in self.impes_drivers[state_to_optim].qm_data_points:
                                
                            # #     _, dist, _ = self.calculate_distance_to_ref(molecule.get_coordinates_in_bohr(), dp.cartesian_coordinates)
                            # #     _, dist_opt, _ = self.calculate_distance_to_ref(optimized_molecule.get_coordinates_in_bohr(), dp.cartesian_coordinates)
                            # #     distances_dp_to_org.append(dist)
                            # #     distances_dp_to_opt.append(dist_opt)

                            # # if any(dist for dist in distances_dp_to_org) < distance_to_org:
                            # #     min_distance = min(distances_dp_to_org)
                            # #     min_idx = distances_dp_to_org.index(min_distance)

                            # #     cart_coord = self.impes_drivers[state_to_optim].qm_data_points[min_idx].cartesian_coordinates


                            # #     _, dist, distance_vec = self.calculate_distance_to_ref(optimized_molecule.get_coordinates_in_bohr(), cart_coord)

                            #     self.opt_mols_org_mol_swap[state_to_optim] = molecule
                            #     for imp_coord in backup_constraint:
                            #         if imp_coord not in main_constraint_list:
                            #             main_constraint_list.append(imp_coord)
                            #             break
                            #     same = True
                            
                            # if len(main_constraint_list) != old_const_len:
                            #     old_const_len = len(main_constraint_list)
                            # else:
                            #     same = True
                            #     corr_new_dp_distrib = True
                                
                            # newly_weights = interpolation_driver.weights
                            # new_weights = [value for _, value in newly_weights.items()]
                            # new_used_labels = [label_idx for label_idx, _ in newly_weights.items()]
                            # # Sort labels and weights by descending weight
                            # new_sorted_items = sorted(zip(new_used_labels, new_weights), key=lambda x: x[1], reverse=True)
                        

                            # new_order = [lbl for lbl, _ in new_sorted_items[:len(old_sorted_shorted)]]
                            
                            # same = (new_order == old_sorted_shorted)

                            # # quick check: exact same ordering?
                
                            # print(new_sorted_items, sorted_items, same)
                            # int_diff = []
                            # elem_list = []
                            # if not same:
                            #     for element_idx, element in enumerate(self.root_z_matrix[self.current_state]):
                                    
                            #         target = set(element[1:3])   # positions 2 and 3 → indices 1 and 2

                            #         matches = [p for p in self.all_rot_bonds if set(p).issubset(target)]

                            #         if len(element) == 4 and element not in constraints and len(matches) != 0:

                            #             diff = np.sin(interpolation_driver.impes_coordinate.internal_coordinates_values[element_idx] - self.impes_drivers[state_to_optim].impes_coordinate.internal_coordinates_values[element_idx])
                            #             if abs(diff) > 1e-1:
                            #                 elem_list.append(element)
                            #                 int_diff.append(diff)
                                    
                            #     ordered_diif = sorted(zip(elem_list, int_diff), key=lambda x:x[1], reverse=True)
                                
                            #     print('Len ordered list', ordered_diif)
                                
                            #     if len(ordered_diif) > 0:
                            #         for spec_elem in range(counter):
                            #             if spec_elem == len(ordered_diif):
                            #                 break
                                        
                            #             dihedral_val_to_reset = molecule.get_dihedral_in_degrees((ordered_diif[spec_elem][0][0] + 1, ordered_diif[spec_elem][0][1] + 1, ordered_diif[spec_elem][0][2] + 1, ordered_diif[spec_elem][0][3] + 1))
                                        
                            #             print(dihedral_val_to_reset)
                            #             # print(dihedral_val_to_reset, ordered_diif[spec_elem][0])
                            #             optimized_molecule.set_dihedral_in_degrees((ordered_diif[spec_elem][0][0] + 1, ordered_diif[spec_elem][0][1] + 1, ordered_diif[spec_elem][0][2] + 1, ordered_diif[spec_elem][0][3] + 1), dihedral_val_to_reset)
                                        
                            #             print(optimized_molecule.get_dihedral_in_degrees((ordered_diif[spec_elem][0][0] + 1, ordered_diif[spec_elem][0][1] + 1, ordered_diif[spec_elem][0][2] + 1, ordered_diif[spec_elem][0][3] + 1)))      
                            #             print('Orderd diff element', ordered_diif[spec_elem], old_list_change)
                            #             if ordered_diif[spec_elem][0] in old_list_change:
                            #                 same = True
                            #     else:
                            #         same = True
                                
                            #     counter +=1
                            #     if len(old_list_change) == 0:
                            #         old_list_change = elem_list.copy()
                            #         same = True
                        corr_new_dp_distrib = True
                    state_specific_molecules.append((optimized_molecule, current_basis, [state_to_optim], imp_int_coord))
                    
                print('New optimized molecule \n', optimized_molecule.get_xyz_string())
                self.add_point(state_specific_molecules, self.non_core_symmetry_groups)
                self.last_point_added = self.point_checker - 1
                self.point_checker = 0
                # self.point_adding_molecule[self.step] = (molecule, qm_energy, label)

            else:
                    
                current_basis = None
                if self.current_state == 0:
                    current_basis = MolecularBasis.read(molecule, self.basis_set_label['gs'])
                else:
                    current_basis = MolecularBasis.read(molecule, self.basis_set_label['es'])
                state_specific_molecules.append((molecule, current_basis, addition_of_state_specific_points, []))
                    
                
                self.add_point(state_specific_molecules, self.non_core_symmetry_groups)
                self.last_point_added = self.point_checker - 1
                self.point_checker = 0
                # self.point_adding_molecule[self.step] = (molecule, qm_energy, label)


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

        adjusted_molecule = {'gs': [], 'es': []}
        symmetry_mapping_groups = []
        symmetry_exclusion_groups = []
        category_label = None

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
                    optim_dih_0based.extend(entries[3])

                    # if isinstance(self.drivers['gs'][0], ScfRestrictedDriver):
                    #     _, scf_results, _ = self.compute_energy(self.drivers['gs'][0], cur_molecule, current_basis)
                    #     optimized_molecule, opt_results = self._run_optimization(
                    #                     self.drivers['gs'][1],
                    #                     cur_molecule,
                    #                     constraints=optim_dih_0based,
                    #                     index_offset=1,
                    #                     compute_args=(current_basis, scf_results),
                    #                     source_molecule=entries[0]
                    #                 )
                    #     cur_molecule = optimized_molecule
                        
                    # elif isinstance(self.drivers['gs'][0], XtbDriver):
                    #     optimized_molecule, opt_results = self._run_optimization(
                    #                     self.drivers['gs'][1],
                    #                     cur_molecule,
                    #                     constraints=optim_dih_0based,
                    #                     index_offset=1,
                    #                     source_molecule=entries[0]
                    #                 )
                    #     cur_molecule = optimized_molecule
                    # print('Optimized molecule', cur_molecule.get_xyz_string())
                    # current_basis = MolecularBasis.read(cur_molecule, entries[1].get_main_basis_label())
                    
                    adjusted_molecule['gs'].append((cur_molecule, current_basis, periodicites[dihedral], dihedral_to_change, entries[2], symmetry_point, entries[3]))
            
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

        
        for key, entries in adjusted_molecule.items():
            if len(entries) == 0:
                continue
            print(key, entries)
            drivers = self.drivers[key]
            
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
                
                energies, scf_results, rsp_results = self._compute_energy_mpi_safe(drivers[0], mol_basis[0], mol_basis[1])
                
                if isinstance(drivers[0], LinearResponseEigenSolver) or isinstance(drivers[0], TdaEigenSolver):
                    energies = energies[mol_basis[4]]

                gradients = self._compute_gradient_mpi_safe(drivers[1], mol_basis[0], mol_basis[1], scf_results, rsp_results)
                
                hessians = self._compute_hessian_mpi_safe(drivers[2], mol_basis[0], mol_basis[1])
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
                    

                    
                    sampled_space = False
                    sampled_distances = []
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
                        impes_coordinate.mapping_masks = masks
                    else:
                        org_mask = [i for i in range(len(impes_coordinate.z_matrix))]
                        masks = [org_mask]
                        impes_coordinate.mapping_masks = masks
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
                        self.qm_energies_dict[mol_basis[4][number]].append(energies[number])
                        
                        self.impes_drivers[mol_basis[4][number]].qm_data_points = self.qm_data_point_dict[mol_basis[4][number]]
                        self.impes_drivers[mol_basis[4][number]].labels = self.sorted_state_spec_im_labels[mol_basis[4][number]]
                        self._refresh_interpolation_driver_caches(mol_basis[4][number])
                        print('I am writing the file', category_label, self.interpolation_settings[mol_basis[4][number]]['imforcefield_file'])
                    
                        self.density_around_data_point[0][mol_basis[4][number]] += 1

                    else:
                        
                        self.qm_symmetry_datapoint_dict[mol_basis[4][number]][category_label].append(impes_coordinate)
                        self.impes_drivers[mol_basis[4][number]].qm_symmetry_data_points[category_label].append(impes_coordinate)
                        self._refresh_interpolation_driver_caches(mol_basis[4][number])

                    
                    impes_coordinate.confidence_radius = self.use_opt_confidence_radius[2]

                    print('data list', len(self.allowed_molecules[mol_basis[4][number]]['molecules']))
                    
                    if self._mpi_is_root() or (not self._mpi_is_active()):
                        impes_coordinate.write_hdf5(self.interpolation_settings[mol_basis[4][number]]['imforcefield_file'], new_label)
                        if self.sampling_enabled:
                            self._write_sampling_point_from_geometry(
                                root=mol_basis[4][number],
                                molecule=mol_basis[0],
                                label=new_label,
                                template_point=impes_coordinate
                            )
                    # if mol_basis[5]:
                    #     if self.use_opt_confidence_radius[0] and len(self.allowed_molecules[self.current_state]['molecules']) >= 10:
                            
                    #         trust_radius = None
                    #         sym_dict = self.non_core_symmetry_groups['gs']
                    #         if mol_basis[4][number] > 0 or root == 1 and self.drivers['es'] is not None and self.drivers['es'][0].spin_flip:
                    #             sym_dict = self.non_core_symmetry_groups['es']
                            
                            
                    #         for run_through in range(1):
                    #             indices = [i for i in range(len(self.allowed_molecules[mol_basis[4][number]]['molecules']) - 1)]
                    #             chosen_structures = [self.allowed_molecules[mol_basis[4][number]]['molecules'][idx] for idx in indices]
                    #             chosen_qm_energies = [self.allowed_molecules[mol_basis[4][number]]['qm_energies'][idx] for idx in indices]
                    #             chosen_im_energies = [self.allowed_molecules[mol_basis[4][number]]['im_energies'][idx] for idx in indices]
                    #             chosen_qm_gradients = [self.allowed_molecules[mol_basis[4][number]]['qm_gradients'][idx] for idx in indices]
                                
                                    
                    #             if self.use_opt_confidence_radius[1] == 'multi_grad':
                                    
                    #                 trust_radius = self.determine_trust_radius_gradient(chosen_structures, 
                    #                                                         chosen_qm_energies,
                    #                                                         chosen_qm_gradients,
                    #                                                         chosen_im_energies, 
                    #                                                         self.qm_data_point_dict[mol_basis[4][number]], 
                    #                                                         self.interpolation_settings[mol_basis[4][number]],
                    #                                                         self.qm_symmetry_datapoint_dict[mol_basis[4][number]],
                    #                                                         sym_dict,
                    #                                                         self.root_z_matrix[mol_basis[4][number]],
                    #                                                         exponent_p_q = (self.impes_drivers[mol_basis[4][number]].exponent_p, self.impes_drivers[mol_basis[4][number]].exponent_q))
                                                                            
                                
                                
                    #             elif self.use_opt_confidence_radius[1] == 'bayes':
                    #                 trust_radius = self.determine_beysian_trust_radius(self.allowed_molecules[mol_basis[4][number]]['molecules'], 
                    #                                                                 self.allowed_molecules[mol_basis[4][number]]['qm_energies'], 
                    #                                                                 self.qm_data_point_dict[mol_basis[4][number]], 
                    #                                                                 self.interpolation_settings[mol_basis[4][number]], 
                    #                                                                 self.qm_symmetry_datapoint_dict[mol_basis[4][number]],
                    #                                                                 sym_dict,
                    #                                                                 self.root_z_matrix[mol_basis[4][number]])

                    #             for idx, trust_radius in enumerate(trust_radius):

                    #                 print(self.sorted_state_spec_im_labels[mol_basis[4][number]][idx])
                    #                 self.qm_data_point_dict[mol_basis[4][number]][idx].update_confidence_radius(self.interpolation_settings[mol_basis[4][number]]['imforcefield_file'], self.sorted_state_spec_im_labels[mol_basis[4][number]][idx], trust_radius)
                    #                 self.qm_data_point_dict[mol_basis[4][number]][idx].confidence_radius = trust_radius

                    #                 if idx == len(trust_radius) - 1:

                    #                     if trust_radius < 0.09:

                    #                         self.qm_data_point_dict[mol_basis[4][number]][idx].remove_point_from_hdf5(self.interpolation_settings[mol_basis[4][number]]['imforcefield_file'], new_label, 
                    #                                                                                                   use_inverse_bond_length=self.qm_data_point_dict[mol_basis[4][number]][idx].use_inverse_bond_length, 
                    #                                                                                                   use_cosine_dihedral=self.qm_data_point_dict[mol_basis[4][number]][idx].use_cosine_dihedral, 
                    #                                                                                                   use_eq_bond_length=self.qm_data_point_dict[mol_basis[4][number]][idx].use_eq_bond_length)

                                
                    #             self.simulation.saveCheckpoint('checkpoint')
                                
                    #             self.output_file_writer(self.summary_output)
                    #             state_spec_dict = {'pot_energies':[], 'gradients':[]}
                    #             self.gloabal_sim_informations = {f'state_{root}':state_spec_dict for root in self.roots_to_follow}
                    #             self.gloabal_sim_informations['coordinates_ang'] = []
                    #             self.gloabal_sim_informations['state'] = []
                    #             self.gloabal_sim_informations['temperatures'] = []
                    # exit()
                    for root in mol_basis[4]:
                        self.write_qm_energy_determined_points_h5(self.allowed_molecules[root]['molecules'][self.last_added: ],
                                                            self.allowed_molecules[root]['qm_energies'][self.last_added:],
                                                            self.allowed_molecules[root]['qm_gradients'][self.last_added:],
                                                            self.allowed_molecules[root]['im_energies'][self.last_added:],
                                                            root)
                        self.last_added = len(self.allowed_molecules[root]['molecules'])
                label_counter += 1
            if any(root > 0 for root in self.roots_to_follow):
                                   
                if isinstance(drivers[0], LinearResponseEigenSolver) or isinstance(drivers[0], TdaEigenSolver):           
                    
                    drivers[1].state_deriv_index = org_roots

            if self.current_state in self.opt_mols_org_mol_swap:
                self._interp_compute_root_local_serial(self.impes_drivers[self.current_state], self.opt_mols_org_mol_swap[self.current_state])
                print(self.impes_drivers[self.current_state].weights, self.impes_drivers[self.current_state].get_energy(), self.opt_mols_org_mol_swap[self.current_state].get_xyz_string())

        
        self.simulation.saveCheckpoint('checkpoint')    
     
    def _write_sampling_point_from_geometry(self, root, molecule, label, template_point):
        sampling_qm, sampling_grad, sampling_hess = self.sampling_driver['gs']
        e, _, _ = self._compute_energy_mpi_safe(sampling_qm, molecule, basis=None)
        g = self._compute_gradient_mpi_safe(sampling_grad, molecule, basis=None, scf_results=None, rsp_results=None)
        h = self._compute_hessian_mpi_safe(sampling_hess, molecule, basis=None)

        masses = molecule.get_masses().copy()
        inv_sqrt = 1.0 / np.sqrt(np.repeat(masses, 3))
        grad_vec = g[0].reshape(-1)
        hess_mat = h[0].reshape(grad_vec.size, grad_vec.size)
        mw_grad = inv_sqrt * grad_vec
        mw_hess = (inv_sqrt[:, None] * hess_mat) * inv_sqrt[None, :]

        dp = InterpolationDatapoint(self.root_z_matrix[root])
        dp.update_settings(self.sampling_interpolation_settings[root])
        dp.cartesian_coordinates = template_point.cartesian_coordinates
        dp.eq_bond_lengths = template_point.eq_bond_lengths
        dp.mapping_masks = getattr(template_point, "mapping_masks", None)
        dp.imp_int_coordinates = getattr(template_point, "imp_int_coordinates", [])
        dp.inv_sqrt_masses = inv_sqrt
        dp.energy = e[0]
        dp.gradient = mw_grad.reshape(g[0].shape)
        dp.hessian = mw_hess.reshape(h[0].shape)
        dp.confidence_radius = template_point.confidence_radius
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
            if '_symmetry' not in label:
                qm_data_point = InterpolationDatapoint(self.root_z_matrix[root])
                qm_data_point.update_settings(self.sampling_interpolation_settings[root])
                qm_data_point.read_hdf5(self.sampling_interpolation_settings[root]['imforcefield_file'], label)
                qm_data_point.inv_sqrt_masses = inv_sqrt_masses

                self.sampling_qm_data_point_dict[root].append(qm_data_point)

                old_label = qm_data_point.point_label
                self.sampling_qm_symmetry_datapoint_dict[root][old_label] = [qm_data_point]
            else:
                sym_dp = InterpolationDatapoint(self.root_z_matrix[root])
                sym_dp.update_settings(self.sampling_interpolation_settings[root])
                sym_dp.read_hdf5(self.sampling_interpolation_settings[root]['imforcefield_file'], label)
                self.sampling_qm_symmetry_datapoint_dict[root][old_label].append(sym_dp)


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

    def _run_optimization(self, optimization_driver, molecule, constraints=None, index_offset=1, compute_args=None, source_molecule=None):

        opt_drv = OptimizationDriver(optimization_driver)
        opt_drv.ostream.mute()
        if constraints is not None:
            opt_drv.constraints = self._build_opt_constraint_list(constraints, index_offset=index_offset)

        if compute_args is None:
            opt_results = self._opt_compute_mpi_safe(opt_drv, molecule)
        else:
            opt_results = self._opt_compute_mpi_safe(opt_drv, molecule, *compute_args)

        if source_molecule is None:
            source_molecule = molecule

        optimized_molecule = Molecule.from_xyz_string(opt_results['final_geometry'])
        optimized_molecule.set_charge(source_molecule.get_charge())
        optimized_molecule.set_multiplicity(source_molecule.get_multiplicity())

        return optimized_molecule, opt_results
    
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
                interpolation_driver.distance_thrsh = 1000
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


    def group_by_connected_components(self, D_db, D_ref, D_dpref, cutoff_distance=3.0):
        """
        Groups datapoints and reference points into distinct basins based on a physical 
        distance threshold, naturally isolating outliers without forcing a cluster count.
        """
        n_dp = D_db.shape[0]
        n_ref = D_ref.shape[0]
        n_tot = n_dp + n_ref
        
        # 1. Assemble the Global Distance Matrix
        D_global = np.block([
            [D_db,       D_dpref],
            [D_dpref.T,  D_ref  ]
        ])
        
        # 2. Construct the Adjacency Matrix (The Graph Edges)
        # If the distance is less than the cutoff, the nodes share an edge (1). 
        # Otherwise, the edge is severed (0).
        adjacency_matrix = (D_global <= cutoff_distance).astype(int)
        
        # Remove self-loops
        np.fill_diagonal(adjacency_matrix, 0)
        
        # Convert to sparse matrix for the graph solver
        graph = csr_matrix(adjacency_matrix)
        
        # 3. Solve for Connected Components
        # This automatically finds exactly how many distinct topological islands exist
        n_components, labels = connected_components(csgraph=graph, directed=False, return_labels=True)
        
        # 4. Extract the Valid Basins
        ref_global_indices = np.arange(n_dp, n_tot)
        optimization_groups = []
        
        for cluster_id in range(n_components):
            # Find references in this isolated island
            refs_in_cluster = [
                (global_idx - n_dp) for global_idx in ref_global_indices 
                if labels[global_idx] == cluster_id
            ]
            
            # Only keep islands that actually contain reference points to drive the fit
            if len(refs_in_cluster) > 0:
                dps_in_cluster = [
                    i for i in range(n_dp) 
                    if labels[i] == cluster_id
                ]
                
                # Only return groups where uncharacterized datapoints actually exist
                if len(dps_in_cluster) > 0:
                    optimization_groups.append({
                        'basin_id': int(cluster_id), # Cast to standard int
                        'datapoint_indices': dps_in_cluster,
                        'reference_indices': refs_in_cluster
                    })
                    
        return optimization_groups



    def determine_trust_radius_gradient(self, molecules, qm_energies, qm_gradients, im_energies, datapoints, interpolation_setting, sym_datapoints, sym_dict, z_matrix, exponent_p_q):
        
        def optimize_trust_radius(alphas, geom_list, E_ref_list, G_ref_list, E_im_list, dps, impes_dict, sym_datapoints, sym_dict, exponent_p_q):
            """
            Perform the gradient-based optimization to find R*
            that minimizes the sum of squared errors to reference QM energies.
            """

            bounds = [(1e-6, 1.5)] * len(dps)
            sample_size = len(geom_list)

            if sample_size < 50:
                size_label = 'small'
                n_bh_iter = 10
                gtol = 1e-4
                ftol = 1e-5
                maxiter = 300
            elif sample_size <= 300:
                size_label = 'medium'
                n_bh_iter = 4
                gtol = 2e-4
                ftol = 5e-5
                maxiter = 220
            else:
                size_label = 'large'
                n_bh_iter = 0
                gtol = 5e-4
                ftol = 1e-4
                maxiter = 140

             
            opt = AlphaOptimizer(z_matrix, impes_dict, sym_dict, sym_datapoints, dps,
                 geom_list, E_ref_list, G_ref_list, exponent_p_q,
                 e_x=self.use_opt_confidence_radius[2],
                 beta=0.8, n_workers=os.cpu_count())  # pick sensible n_workers

            minimizer_kwargs = {
                "method": "L-BFGS-B",
                "jac": opt.jac,
                "bounds": bounds,
                "options": {"disp": True, "gtol": gtol, "ftol": ftol, "maxls": 10, "maxiter": maxiter}
            }

            print(f"Trust-radius optimization regime: {size_label} (S={sample_size}, basinhopping niter={n_bh_iter})")
            try:
                if n_bh_iter > 0:
                    res = basinhopping(opt.fun, x0=alphas, minimizer_kwargs=minimizer_kwargs, niter=n_bh_iter)
                else:
                    res = minimize(opt.fun, x0=alphas, jac=opt.jac, method="L-BFGS-B", bounds=bounds, options=minimizer_kwargs["options"])
                print(res)
                return res
            finally:
                opt.close()
        
        inital_alphas = [dp.confidence_radius for dp in datapoints]

        print('INPUT Trust radius', inital_alphas)

        D_dp, D_ref, D_dpref = self._build_dp_ref_distance_matrix(datapoints, molecules, sym_dict)

        groups = self.group_by_connected_components(D_dp, D_ref, D_dpref)
        print(groups)

        # additional_molecules = []
        # additional_energies = []
        # additional_gradients = []
        # for dp in datapoints:
        #     additional_molecules.append(Molecule(molecules[0].get_labels(), dp.cartesian_coordinates, 'bohr'))
        #     additional_energies.append(dp.energy)
        #     additional_gradients.append(dp.gradient)
        # molecules.extend(additional_molecules)
        # qm_energies.extend(additional_energies)
        # qm_gradients.extend(additional_gradients)

        final_alphas = np.ones_like(datapoints)

        for i, dp in enumerate(datapoints):
            
            final_alphas[i] = dp.confidence_radius

        for group in groups:

            dp_mask = group['datapoint_indices']
            ref_mask = group['reference_indices']
            # dp_mask = [i for i in range(len(datapoints))]
            # ref_mask = [i for i in range(len(molecules))]
            datapoints_sub = [datapoints[i] for i in dp_mask]
            molecules_sub = [molecules[i] for i in ref_mask]
            qm_energies_sub = [qm_energies[i] for i in ref_mask]
            qm_gradients_sub = [qm_gradients[i] for i in ref_mask]
            im_energies_sub = [im_energies[i] for i in ref_mask]

            trust_radius = optimize_trust_radius(final_alphas[dp_mask], molecules_sub, qm_energies_sub, qm_gradients_sub, im_energies_sub, datapoints_sub, interpolation_setting, sym_datapoints, sym_dict, exponent_p_q)

            final_alphas[dp_mask] = trust_radius['x']
            print('FINAL Trust radius', trust_radius['x'])

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
        
        all_dihedrals = [element for element in z_matrix if len(element) == 4]

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
