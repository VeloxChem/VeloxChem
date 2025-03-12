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

from mpi4py import MPI
from io import StringIO
from contextlib import redirect_stderr
from pathlib import Path
import numpy as np
import sys
import time
import h5py

from .veloxchemlib import mpi_master
from .solvationbuilder import SolvationBuilder
from .outputstream import OutputStream
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

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.

    Instance variables:
        - padding: The padding for the solvation box.
        - solvent_name: The name of the solvent.
        - temperature: The temperature for the simulation.
        - pressure: The pressure for the simulation.
        - timestep: The timestep for the simulation.
        - num_equil_steps: The number of steps for the equilibration.
        - num_steps: The number of steps for the simulation.
        - number_of_snapshots: The number of snapshots to save.
        - cutoff: The cutoff for the nonbonded interactions.
        - constraints: The constraints for the simulation.
        - nonbondedMethod: The nonbonded method for the simulation.
        - solute_ff: The solute force field.
        - platform: The platform for the openMM simulation.
        - save_trajectory_xtc: If True, save the trajectories in XTC format.
        - alpha: The alpha parameter for the Gaussian softcore potential.
        - beta: The beta parameter for the Gaussian softcore potential.
        - x: The x parameter for the Gaussian softcore potential.
        - lambdas_stage1: The lambdas for stage 1.
        - lambdas_stage2: The lambdas for stage 2.
        - lambdas_stage3: The lambdas for stage 3.
        - lambdas_stage4: The lambdas for stage 4.
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
        self.output_folder = Path("solvation_fep_output")

        # Options for the SolvationBuilder
        self.padding = 2.0
        self.solvent_name = 'spce'
        
        # Ensemble and MD options
        self.temperature = 298.15 * unit.kelvin
        self.pressure = 1 * unit.atmospheres
        self.timestep = 1.0 * unit.femtoseconds 
        self.num_equil_steps = 10000 #10 ps
        self.num_steps = 1000000 # 1 ns
        self.number_of_snapshots = 1000
        self.cutoff = 1.0 * unit.nanometers
        self.constraints = app.HBonds  
        self.nonbondedMethod = app.PME  

        # Other objects
        self.solute_ff = None
        self.platform = None
        self.save_trajectory_xtc = False

        # Alchemical parameters
        self.alpha = 5 * unit.kilocalories_per_mole
        self.beta = 5
        self.x = 4
        # Single parameter for lambdas for stage 1
        # Set to 6 to be on the safe side
        self.lambdas_stage1 = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
        # Asymetric lambdas for stage 2
        self.lambdas_stage2 = [1.0, 0.8, 0.6, 0.4, 0.3, 0.2, 0.15, 0.10, 0.05, 0.03, 0.0]
        # Fixed lambda vector for stage 3 with 6 lambdas
        self.lambdas_stage3 = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
        # Fixed lambda vector for stage 4 with 3 lambdas
        self.lambdas_stage4 = [1.0, 0.5, 0.0] 

        # Storage for potential energies across stages
        self.u_kln_matrices = []
        self.stage = 1

        # Final energies
        self.final_free_energy = 0.0
        self.delta_f = []

    def compute_solvation(self, molecule, ff_gen_solute=None, solvent='spce', solvent_molecule=None, ff_gen_solvent=None, target_density=None):
        """
        Run the solvation free energy calculation using OpenMM.

        :param molecule:
            The VeloxChem molecule object to solvate.
        :param ff_gen_solute:
            The force field generator for the solute. If None, it will be calculated.
        :param solvent:
            The solvent to use for solvation. Default is spce water.
            Available options: 'spce', 'tip3p', 'ethanol', 'methanol', 'acetone', 'chloroform', 
            'hexane', 'toluene', 'dcm', 'benzene', 'dmso', 'thf', 'acetonitrile', 'other' or 'itself'.
        :param solvent_molecule:
            The VeloxChem molecule object for the solvent. Mandatory for 'other'.
        :param ff_gen_solvent:
            The force field generator for the solvent. If None, it will be calculated.
        :param target_density:
            The target density for the solvent. Mandatory for 'other' or 'itself' solvent options.
        :return:
            A list containing the free energy calculations for each stage.
        """
        
        sol_builder = SolvationBuilder()

        sol_builder.solvate(solute=molecule, 
                            solvent=solvent,
                            solvent_molecule=solvent_molecule, 
                            padding=self.padding,
                            target_density=target_density, 
                            neutralize=False, 
                            equilibrate=False)
        
        self.solvent_name = solvent
        self.output_folder.mkdir(parents=True, exist_ok=True)
        # Note: GROMACS files will be used instead of OpenMM files.
        # For some reason, the results are more consistent. (?)
        sol_builder.write_gromacs_files(ff_gen_solute, ff_gen_solvent)

        if not ff_gen_solute:
            self.solute_ff = sol_builder.solute_ff
        else:
            self.solute_ff = ff_gen_solute

        delta_f, final_free_energy = self._run_stages()
        u_kln = self.u_kln_matrices

        return delta_f, final_free_energy, u_kln
        
    def compute_solvation_from_omm_files(self, system_pdb, solute_pdb, solute_xml, other_xml_files):
        """
        Run the solvation free energy calculation using OpenMM.
        :param system_pdb:
            The PDB file with the system coordinates.
        :param solute_pdb:
            The PDB file with the solute coordinates.
        :param solute_xml:
            The XML file with the solute force field.
        :param other_xml_files:
            A list with the XML files needed for the rest of the system.
        :return:
            A list containing the free energy calculations for each stage.
        """

        # Special name for the solvent
        self.solvent_name = 'omm_files'

        self.system_pdb = system_pdb
        self.solute_pdb = solute_pdb
        self.solute_xml = solute_xml
        self.other_xml_files = other_xml_files
        self.output_folder.mkdir(parents=True, exist_ok=True)

        delta_f, final_free_energy = self._run_stages()
        u_kln = self.u_kln_matrices

        return delta_f, final_free_energy, u_kln
        
    def compute_solvation_from_gromacs_files(self, system_gro, system_top, solute_gro, solute_top):
        """
        Run the solvation free energy calculation using OpenMM.
        
        :param system_gro:
            The GRO file with the system coordinates.
        :param system_top:
            The TOP file with the system topology.
        :param solute_gro:
            The GRO file with the solute coordinates.
        :param solute_top:
            The TOP file with the solute topology.
        :return:
            A list containing the free energy calculations for each stage.
        """

        # Special name for the solvent
        self.solvent_name = 'gro_files'

        self.system_gro = system_gro
        self.system_top = system_top
        self.solute_gro = solute_gro
        self.solute_top = solute_top
        self.output_folder.mkdir(parents=True, exist_ok=True)

        delta_f, final_free_energy = self._run_stages()
        u_kln = self.u_kln_matrices

        return delta_f, final_free_energy, u_kln

    def _run_stages(self):
        """
        Manage the simulation stages.
        """
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
        self.ostream.print_line("Thermodynamic Ensemble: NPT")
        simulation_time = self.num_steps * self.timestep.value_in_unit(unit.nanoseconds) 
        self.ostream.print_blank()
        self.ostream.print_line("Alchemical Parameters:")
        self.ostream.print_line("-"*len("Alchemical Parameters:"))
        self.ostream.print_line(f"Number of equilibration steps: {self.num_equil_steps}")
        self.ostream.print_line(f"Number of steps per lambda simulation: {self.num_steps}")
        self.ostream.print_line(f"Number of snapshots per lambda simulation: {self.number_of_snapshots}")
        # Simulation time per lambda in ns with 2 decimal places
        self.ostream.print_line(f"Simulation time per lambda: {simulation_time:.2f} ns")
        self.ostream.print_line(f"Lambdas in Stage 1: {self.lambdas_stage1}")
        self.ostream.print_line(f"Lambdas in Stage 2: {self.lambdas_stage2}")
        self.ostream.print_line(f"Lambdas in Stage 3: {self.lambdas_stage3}")
        self.ostream.print_line(f"Lambdas in Stage 4: {self.lambdas_stage4}")
        self.ostream.print_blank()

        # Run the simulations
        self.ostream.print_info("Starting solvated simulation (Stage 1)...\n")
        delta_f_1, free_en_s1 = self._run_lambda_simulations(stage=1)
        self.ostream.flush()
        
        self.ostream.print_info("Starting removing GSC potential (Stage 2)...\n")
        delta_f_2, free_en_s2 = self._run_lambda_simulations(stage=2)
        self.ostream.flush()

        self.ostream.print_info("Starting vacuum simulation (Stage 3)...\n")
        delta_f_3, free_en_s3 = self._run_lambda_simulations(stage=3, vacuum=True)
        self.ostream.flush()

        self.ostream.print_info("Starting removing GSC potential in vacuum (Stage 4)...\n")
        delta_f_4, free_en_s4 = self._run_lambda_simulations(stage=4, vacuum=True)
        self.ostream.flush()

        # Calculate the final free energy
        final_free_energy = free_en_s1 + free_en_s2 - (free_en_s3 + free_en_s4)
        self.ostream.print_line(f"Final free energy: {final_free_energy} kJ/mol")
        self.ostream.flush()

        self.delta_f = [delta_f_1, delta_f_2, delta_f_3, delta_f_4]
        self.final_free_energy = final_free_energy

        return [delta_f_1, delta_f_2, delta_f_3, delta_f_4], final_free_energy

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

        for lam in lambdas:
            lam = round(lam, 2)
            self.ostream.print_info(f"Running lambda = {lam}, stage = {stage}...")
            self.ostream.flush()

            # Run the simulation for the current lambda
            if vacuum:
                system, topology, positions = self._create_vacuum_system(lam)
            else:
                system, topology, positions = self._create_solvated_system(lam)

            ns_simulated, sim_time = self._run_simulation(system, topology, positions, lam)
            
            total_ns_simulated += ns_simulated
            total_sim_time += sim_time

            # Performance calculation in ns/hour
            ns_per_hour = (ns_simulated / sim_time) * 3600
            self.ostream.print_info(f'Lambda = {lam} completed. Performance: {ns_per_hour:.2f} ns/hour')
            self.ostream.flush()

        # Print total performance metrics for the stage
        overall_ns_per_hour = (total_ns_simulated / total_sim_time) * 3600
        self.ostream.print_info(f"Total simulation time for stage {stage}: {total_sim_time:.2f} s. "
                                f"Total simulated time: {total_ns_simulated:.2f} ns.")
        self.ostream.print_info(f"Overall performance for stage {stage}: {overall_ns_per_hour:.2f} ns/hour")
        self.ostream.flush()

        # Time the energy recalculation stage
        self.ostream.print_info(f"Recalculating energies for stage {stage}...")
        self.ostream.flush()
        recalculation_start = time.time()
        trajectories, forcefields = self._fetch_files(lambdas=lambdas, stage=stage)

        u_kln = self._recalculate_energies(trajectories, forcefields, topology)

        self.u_kln_matrices.append(u_kln)

        recalculation_end = time.time()
        recalculation_time = recalculation_end - recalculation_start

        self.ostream.print_info(f"Energy recalculation for stage {stage} took {recalculation_time:.2f} seconds.")
        self.ostream.flush()

        self.ostream.print_info(f"Calculating the free energy with MBAR for stage {stage}...")
        self.ostream.flush()
        delta_f = self._calculate_free_energy(u_kln)
        free_energy = delta_f['Delta_f'][-1, 0]
        self.ostream.print_line(f"Free energy for stage {stage}: {delta_f['Delta_f'][-1, 0]:.4f} +/- {delta_f['dDelta_f'][-1, 0]:.4f} kJ/mol")
        self.ostream.flush()

        return delta_f, free_energy
    
    # This could be deleted
    def _run_single_lambda_simulation(self, args):
        """
        Runs a single lambda simulation. Can be run in parallel.

        """
        lam, stage, vacuum = args
        lam = round(lam, 2)
        self.ostream.print_info(f"Running lambda = {lam}, stage = {stage}...")
        self.ostream.flush()

        if vacuum:
            system, topology, positions = self._create_vacuum_system(lam)
        else:
            system, topology, positions = self._create_solvated_system(lam)

        self._run_simulation(system, topology, positions, lam)

    def _create_solvated_system(self, lambda_val):
        """
        Create solvated system with alchemical scaling.
        """
        # Reading the system from files
        if self.solvent_name == 'omm_files':
            # Under investigation!
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
                gro = app.GromacsGroFile('system.gro')
                initial_system_ff = app.GromacsTopFile('system.top', 
                                                        periodicBoxVectors=gro.getPeriodicBoxVectors())
                topology = initial_system_ff.topology
                positions = gro.positions
            else:
                gro = app.GromacsGroFile('liquid.gro')
                initial_system_ff = app.GromacsTopFile('liquid.top',
                                                    periodicBoxVectors=gro.getPeriodicBoxVectors())
                topology = initial_system_ff.topology
                positions = gro.positions
        
        # createSystem from a top file does not require the topology as an argument
        if self.solvent_name == 'omm_files':
            sys_arguments = [topology, self.nonbondedMethod, self.cutoff, self.constraints]
        else:
            sys_arguments = [self.nonbondedMethod, self.cutoff, self.constraints]

        solvated_system = initial_system_ff.createSystem(*sys_arguments)

        # Add barostat
        solvated_system.addForce(mm.MonteCarloBarostat(self.pressure, self.temperature))

        # Define alchemical regions
        alchemical_region, chemical_region = self._get_alchemical_regions(topology)

        # Custom nonbonded force (Gaussian Softcore)
        gsc_force = self._get_gaussian_softcore_force(lambda_val)

        # Scaling the nonbonded interactions
        for force in solvated_system.getForces():
            if isinstance(force, mm.NonbondedForce):
                gsc_force.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
                gsc_force.setCutoffDistance(self.cutoff)

                for i in range(force.getNumParticles()):
                    charge, sigma, epsilon = force.getParticleParameters(i)
                    gsc_force.addParticle([sigma])

                    if i in alchemical_region:
                        if self.stage == 1:
                            force.setParticleParameters(i, charge * (1 - lambda_val), sigma, epsilon * (1 - lambda_val))
                        elif self.stage == 2:
                            force.setParticleParameters(i, 0, 0, 0)
                
                # Handle exceptions (interaction exclusions)
                for i in range(force.getNumExceptions()):
                    p1, p2, ch, si, ep = force.getExceptionParameters(i)
                    gsc_force.addExclusion(p1, p2)

        # Add interaction group for GSC potential
        gsc_force.addInteractionGroup(chemical_region, alchemical_region)
        solvated_system.addForce(gsc_force)

        return solvated_system, topology, positions

    def _create_vacuum_system(self, lambda_val):
        """
        Create vacuum system with alchemical scaling.
        """
        # Reading the system from files
        if self.solvent_name == 'omm_files':
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
            self.solute_ff.write_gromacs_files('solute_vacuum', 'MOL')
            gro = app.GromacsGroFile('solute_vacuum.gro')
            forcefield_solute = app.GromacsTopFile('solute_vacuum.top')
            topology = forcefield_solute.topology
            positions = gro.positions

        if self.solvent_name == 'omm_files':
            vacuum_system = forcefield_solute.createSystem(topology, nonbondedMethod=app.NoCutoff, constraints=self.constraints)

        else:
            vacuum_system = forcefield_solute.createSystem(nonbondedMethod=app.NoCutoff, constraints=self.constraints)
        
        gsc_force = self._get_gaussian_softcore_force(lambda_val)

        for force in vacuum_system.getForces():
            if isinstance(force, mm.NonbondedForce):
                for i in range(force.getNumParticles()):
                    charge, sigma, epsilon = force.getParticleParameters(i)
                    gsc_force.addParticle([sigma])

                    if self.stage == 3:
                        force.setParticleParameters(i, charge * (1 - lambda_val), sigma, epsilon * (1 - lambda_val))
                    elif self.stage == 4:
                        force.setParticleParameters(i, 0, 0, 0)
                
                # Handle exceptions (interaction exclusions)
                for i in range(force.getNumExceptions()):
                    p1, p2, ch, si, ep = force.getExceptionParameters(i)
                    gsc_force.addExclusion(p1, p2)

        vacuum_system.addForce(gsc_force)

        return vacuum_system, topology, positions

    def _get_gaussian_softcore_force(self, lambda_val):
        """
        Define the Gaussian softcore potential.
        """
        gsc_energy = 'lambda * alpha * exp(-beta * (r/sigma)^x); sigma = 0.5 * (sigma1 + sigma2);'
        gsc_force = mm.CustomNonbondedForce(gsc_energy)
        gsc_force.addGlobalParameter('lambda', lambda_val)
        gsc_force.addGlobalParameter('alpha', self.alpha)
        gsc_force.addGlobalParameter('beta', self.beta)
        gsc_force.addGlobalParameter('x', self.x)
        gsc_force.addPerParticleParameter('sigma')

        return gsc_force

    def _run_simulation(self, system, topology, positions, lambda_value):
        """
        Run the simulation using OpenMM. Return simulated time in ns and real elapsed time in seconds.
        """
        integrator = mm.LangevinIntegrator(self.temperature, 1.0 / unit.picoseconds, self.timestep)

        if self.platform:
            platform = mm.Platform.getPlatformByName(self.platform)
            simulation = app.Simulation(topology, system, integrator, platform)
        else:
            simulation = app.Simulation(topology, system, integrator)
        
        simulation.context.setPositions(positions)

        # Minimize energy
        simulation.minimizeEnergy()

        # Equilibration
        simulation.step(self.num_equil_steps)

        # Calculate production run time and interval
        interval = self.num_steps // self.number_of_snapshots
        simulated_ns = self.num_steps * self.timestep / unit.nanoseconds  # Simulated time in ns

        # Setup reporters
        output_folder = Path(self.output_folder)  
        # XTC reporter is optional because it affects performance!
        if self.save_trajectory_xtc:
            trajectory_filename = output_folder / f'trajectory_{lambda_value}_stage{self.stage}.xtc'
            simulation.reporters.append(app.XTCReporter(str(trajectory_filename), interval))
            self.ostream.print_info(f"Trajectory will be saved in {str(trajectory_filename)}")
            self.ostream.flush()
        # State data reporter
        energy_filename = output_folder / f'energy_{lambda_value}_stage{self.stage}.txt'
        simulation.reporters.append(app.StateDataReporter(str(energy_filename), interval,
                                                        step=True, potentialEnergy=True, temperature=True))

        # Write the system to a file
        if self.stage in [1, 2]:
            system_filename = output_folder / f'solvated_system_{lambda_value}_stage{self.stage}.xml'
            with open(system_filename, 'w') as f:
                f.write(mm.XmlSerializer.serialize(system))
            self.ostream.print_info(f"System XML written to {system_filename}")
            self.ostream.flush()
            # Name for the hdf5 file
            hdf5_filename = output_folder / f'positions_{lambda_value}_stage{self.stage}.h5'

        
        elif self.stage in [3, 4]:
            system_filename = output_folder / f'vacuum_system_{lambda_value}_stage{self.stage}.xml'
            with open(system_filename, 'w') as f:
                f.write(mm.XmlSerializer.serialize(system))
            self.ostream.print_info(f"System XML written to {system_filename}")
            
            # Name for the hdf5 file
            hdf5_filename = output_folder / f'positions_{lambda_value}_stage{self.stage}.h5'
            
        # Start the production run
        time_start = time.time()
        # Create the hdf5 file
        hf = h5py.File(hdf5_filename, 'w')
        
        # Snapshot counter for the hdf5 file labels
        snapshot = 0
        for s in range(self.num_steps):
            simulation.step(1)
            if (s + 1) % interval == 0:
                # Save the positions in the hdf5 file after stepping
                snapshot += 1
                state = simulation.context.getState(getPositions=True)
                positions = state.getPositions().value_in_unit(unit.nanometers)
                label = f'positions_{snapshot}'
                hf.create_dataset(label, data=positions) 
        hf.close()

        time_end = time.time()

        elapsed_time = time_end - time_start  # Elapsed real-world time in seconds

        self.ostream.print_info(f'Lambda {lambda_value} simulation time: {simulated_ns:.2f} ns. Elapsed time: {elapsed_time:.2f} s')
        self.ostream.flush()
        # Performance in ns/hour
        ns_per_hour = (simulated_ns / elapsed_time) * 3600
        self.ostream.print_info(f'Performance: {ns_per_hour:.2f} ns/hour')
        self.ostream.flush()

        return simulated_ns, elapsed_time

    def _fetch_files(self, lambdas, stage):
        """
        Fetch trajectories and forcefields generated during the simulation.
        """
        trajectories = []
        forcefields = []
        for lam in lambdas:
            lam = round(lam, 2)
            positions_path = self.output_folder / f'positions_{lam}_stage{stage}.h5'
            trajectories.append(str(positions_path))

            if stage in [1, 2]:
                forcefield_path = self.output_folder / f'solvated_system_{lam}_stage{stage}.xml'
            elif stage in [3, 4]:
                forcefield_path = self.output_folder / f'vacuum_system_{lam}_stage{stage}.xml'
            
            forcefields.append(str(forcefield_path))

        return trajectories, forcefields

    def _recalculate_energies(self, trajectories, forcefields, topology):

        u_kln = np.zeros((len(trajectories), len(forcefields), self.number_of_snapshots))

        # k loop for U_kln
        for k, trajectory in enumerate(trajectories):
            time_start_trajectory = time.time()
            self.ostream.print_info(f"Recalculating energies for trajectory {trajectory}...")
            self.ostream.flush()
            # Load the positions from the hdf5 file for recalculation
            positions = []
            with h5py.File(trajectory, 'r') as hf:
                positions = [hf[f'positions_{n+1}'][:] for n in range(self.number_of_snapshots)]

            # l loop for U_kln
            for l in range(len(forcefields)):
                time_start_forcefield = time.time()
                self.ostream.print_info(f"Recalculating energies for forcefield {l}...") 
                self.ostream.flush()
                integrator = mm.VerletIntegrator(1.0 * unit.femtoseconds)

                if self.platform:
                    platform = mm.Platform.getPlatformByName(self.platform)
                    simulation = app.Simulation(topology, forcefields[l], integrator, platform)
                else:
                    simulation = app.Simulation(topology, forcefields[l], integrator)

                # n loop for U_kln
                for n in range(self.number_of_snapshots):
                    positions_with_units = positions[n] * unit.nanometers
                    simulation.context.setPositions(positions_with_units)
                    state = simulation.context.getState(getEnergy=True)
                    u_kln[k, l, n] = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)

                time_end_forcefield = time.time()
                elapsed_time_forcefield = time_end_forcefield - time_start_forcefield
                self.ostream.print_info(f"Forcefield {forcefields[l]} energy recalculation took {elapsed_time_forcefield:.2f} seconds.")
                self.ostream.flush()

                
            time_end_trajectory = time.time()
            elapsed_time_trajectory = time_end_trajectory - time_start_trajectory
            self.ostream.print_info(f"Trajectory {trajectory} energy recalculation took {elapsed_time_trajectory:.2f} seconds.")
            self.ostream.flush()

            # Delete the pdb file to save space
            Path.unlink(trajectory)

        return u_kln

    def _calculate_free_energy(self, u_kln):
        n_states = u_kln.shape[0]
        N_k = np.zeros(n_states)
        subsample_indices = []

        for k in range(n_states):
            t0, g, Neff_max = timeseries.detect_equilibration(u_kln[k, k, :])
            u_equil = u_kln[k, k, t0:]
            indices = timeseries.subsample_correlated_data(u_equil, g=g)
            subsample_indices.append(indices + t0)
            N_k[k] = len(indices)

        u_kln_reduced = np.zeros_like(u_kln)
        for k in range(n_states):
            for l in range(n_states):
                u_kln_reduced[k, l, :int(N_k[k])] = u_kln[k, l, subsample_indices[k]]

        mbar = MBAR(u_kln_reduced, N_k, solver_protocol='robust')
        delta_f = mbar.compute_free_energy_differences()

        return delta_f

    def _get_alchemical_regions(self, topology):
        """
        Define alchemical (perturbed) and chemical (unperturbed) regions.
        """
        alchemical_region = []
        chemical_region = []

        for res in topology.residues():
            if res.index == 0:
                for atom in res.atoms():
                    alchemical_region.append(atom.index)
            else:
                for atom in res.atoms():
                    chemical_region.append(atom.index)
        return alchemical_region, chemical_region

