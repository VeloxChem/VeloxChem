# SPDX-License-Identifier: BSD-3-Clause
"""OpenMM bridge for serial, unconstrained, isolated Cartesian potentials.

OpenMM imports are deliberately local so single-point potentials do not
require OpenMM. State XML stores integrator velocities; JSON identifies the
reconstructed potential and dynamics. Live calculators are not serialized.
"""

import json
from pathlib import Path
from time import perf_counter

import numpy as np

from .molecule import Molecule
from .potential import PotentialResult
from .veloxchemlib import bohr_in_angstrom, hartree_in_kjpermol


class OpenMMPotentialBridge:
    """Persistent callable for PythonForce; returns kJ/mol and kJ/mol/nm."""

    def __init__(self, molecule, potential):
        self.labels = list(molecule.get_labels())
        self.numbers = np.asarray(molecule.get_element_ids(), dtype=int).tolist()
        self.charge = molecule.get_charge()
        self.multiplicity = molecule.get_multiplicity()
        self.potential = potential
        self.model_id = potential.model_id
        self.call_count = 0

    def __call__(self, state):
        from openmm import unit
        positions = state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
        molecule = Molecule(self.labels, np.asarray(positions) * 10.0, 'angstrom')
        molecule.set_charge(self.charge)
        molecule.set_multiplicity(self.multiplicity)
        result = self.potential.compute(molecule, properties=('energy', 'gradient'))
        if result.model_id != self.model_id:
            raise RuntimeError('Potential identity changed during OpenMM evaluation.')
        checked = PotentialResult(result.energy, result.gradient, result.model_id)
        if checked.gradient.shape != (len(self.labels), 3):
            raise ValueError('Potential gradient does not match OpenMM particle count.')
        self.call_count += 1
        scale = hartree_in_kjpermol()
        return checked.energy * scale, -checked.gradient * (10.0 * scale / bohr_in_angstrom())


def create_system_from_potential(dynamics, molecule, potential, *, masses=None,
                                 pbc=False, cell=None):
    """Build a system with exactly one physical force and identity mapping.

    Optional masses are in dalton. Otherwise Molecule.get_masses() is used.
    Bonds in the reporting topology introduce no bonded potential terms.
    """
    import openmm as mm
    from openmm import app, unit

    if dynamics.comm.Get_size() != 1:
        raise ValueError('Potential dynamics supports one MPI rank only.')
    if np.any(pbc) or cell is not None:
        raise NotImplementedError('Potential dynamics is nonperiodic only.')
    if not hasattr(mm, 'PythonForce'):
        raise ImportError('This OpenMM installation does not provide PythonForce.')
    if potential.model_role not in ('full', 'composite'):
        raise ValueError('A bare residual model cannot drive standalone dynamics.')
    if not {'energy', 'gradient'}.issubset(potential.capabilities):
        raise ValueError('Potential dynamics requires energy and Cartesian gradient.')
    numbers = np.asarray(molecule.get_element_ids())
    positions = np.asarray(molecule.get_coordinates_in_angstrom(), dtype=float)
    if (numbers.ndim != 1 or not len(numbers) or not np.all(np.isfinite(numbers))
            or not np.all(numbers == np.floor(numbers)) or
            positions.shape != (len(numbers), 3) or not np.all(np.isfinite(positions))):
        raise ValueError('Invalid molecule atomic numbers or coordinates.')
    masses = np.array(molecule.get_masses() if masses is None else masses, dtype=float, copy=True)
    if masses.shape != numbers.shape or not np.all(np.isfinite(masses)) or np.any(masses <= 0):
        raise ValueError('Masses must be positive finite dalton values, one per atom.')
    # Fail on an incompatible target before constructing an OpenMM context.
    result = potential.compute(molecule)
    result = PotentialResult(result.energy, result.gradient, result.model_id)
    if result.model_id != potential.model_id:
        raise ValueError('Potential result identity differs from its configuration.')
    if result.gradient.shape != positions.shape:
        raise ValueError('Potential gradient does not match molecule.')
    topology = app.Topology()
    residue = topology.addResidue('MOL', topology.addChain())
    atoms = [topology.addAtom(f'{label}{i+1}', app.Element.getByAtomicNumber(int(z)), residue)
             for i, (label, z) in enumerate(zip(molecule.get_labels(), numbers))]
    connectivity = molecule.get_connectivity_matrix()
    for i in range(len(atoms)):
        for j in range(i):
            if connectivity[i, j]:
                topology.addBond(atoms[j], atoms[i])
    system = mm.System()
    for mass in masses:
        system.addParticle(float(mass) * unit.dalton)
    bridge = OpenMMPotentialBridge(molecule, potential)
    force = mm.PythonForce(bridge)
    force.setName(f'VeloxChem potential {potential.model_id}')
    force.setUsesPeriodicBoundaryConditions(False)
    system.addForce(force)
    dynamics.system = system
    dynamics.topology = topology
    dynamics.positions = positions * 0.1 * unit.nanometer
    dynamics._potential_bridge = bridge
    masses.setflags(write=False)
    dynamics._potential_masses = masses
    dynamics._potential_manifest = dict(model_id=potential.model_id,
                                        model=potential.manifest,
                                        atomic_numbers=numbers.astype(int).tolist(),
                                        charge=molecule.get_charge(),
                                        multiplicity=molecule.get_multiplicity(),
                                        masses_dalton=masses.tolist(),
                                        atom_to_particle=list(range(len(numbers))),
                                        periodic=False, constraints=0)
    dynamics.qm_atoms = None
    dynamics.driver_flag = type(potential).__name__
    return system


def run_potential_md(dynamics, *, restart_file, ensemble, temperature, pressure,
                     friction, timestep, nsteps, snapshots, traj_file, state_file,
                     save_last_frame, output_file, minimize, initial_velocities,
                     seed, step_callback):
    """Opt-in native run_md branch; legacy force-field behavior is untouched.

    initial_velocities: physical position-time velocities in nm/ps, NVE only.
    They are shifted once to Verlet's stored half-step velocities. Restored
    XML velocities are already stored velocities and are never shifted again.
    step_callback(step, state) sees the initial state and every completed step.
    """
    import openmm as mm
    from openmm import app, unit

    if ensemble not in ('NVE', 'NVT'):
        raise ValueError('Isolated potential dynamics supports NVE and NVT only.')
    if not isinstance(nsteps, int) or nsteps < 1 or not isinstance(snapshots, int) or snapshots < 0:
        raise ValueError('nsteps must be positive and snapshots nonnegative integers.')
    if not np.isfinite(timestep) or timestep <= 0:
        raise ValueError('timestep must be positive and finite (fs).')
    if not np.isfinite(temperature) or temperature <= 0 or not np.isfinite(friction) or friction < 0:
        raise ValueError('Temperature must be positive and friction nonnegative.')
    if initial_velocities is not None and ensemble != 'NVE':
        raise ValueError('Prescribed physical velocities are currently supported for NVE only.')
    if restart_file is not None and (initial_velocities is not None or seed is not None or minimize):
        raise ValueError('Restart cannot minimize, reseed, or replace saved velocities.')
    if minimize and initial_velocities is not None:
        raise ValueError('Relax first, then provide physical initial velocities at that geometry.')
    if ensemble == 'NVE' and restart_file is None and initial_velocities is None and seed is None:
        raise ValueError('NVE needs prescribed velocities, a thermal seed, or a restart.')
    if seed is not None and (not isinstance(seed, int) or not 0 < seed < 2**31):
        raise ValueError('seed must be an integer in [1, 2**31).')
    if dynamics.system.getNumConstraints() or any(
            dynamics.system.isVirtualSite(i) for i in range(dynamics.system.getNumParticles())):
        raise ValueError('Constraints and virtual sites are unsupported.')
    if dynamics.system.getNumForces() != 1 or dynamics.system.usesPeriodicBoundaryConditions():
        raise ValueError('Potential dynamics requires its sole nonperiodic PythonForce.')
    dynamics.ensemble = ensemble
    dynamics.temperature = temperature * unit.kelvin
    dynamics.friction = friction / unit.picosecond
    dynamics.timestep = timestep * unit.femtosecond
    dynamics.nsteps = nsteps
    integrator = (mm.VerletIntegrator(dynamics.timestep) if ensemble == 'NVE' else
                  mm.LangevinIntegrator(dynamics.temperature, dynamics.friction, dynamics.timestep))
    if seed is not None and ensemble == 'NVT':
        integrator.setRandomNumberSeed(seed)
    configuration = dict(dynamics._potential_manifest, ensemble=ensemble,
                         timestep_fs=timestep, temperature_K=temperature,
                         friction_per_ps=friction,
                         integrator=type(integrator).__name__,
                         velocity_convention='OpenMM integrator stored velocities',
                         openmm_version=mm.__version__)
    if restart_file is not None:
        saved = json.loads(Path(str(restart_file) + '.json').read_text(encoding='utf-8'))
        for key in ('model_id', 'atomic_numbers', 'charge', 'multiplicity', 'masses_dalton',
                    'ensemble', 'timestep_fs', 'temperature_K', 'friction_per_ps', 'integrator'):
            if saved.get(key) != configuration.get(key):
                raise ValueError(f'Restart configuration mismatch: {key}.')
    platform = dynamics._create_platform()
    simulation = app.Simulation(dynamics.topology, dynamics.system, integrator, platform=platform)
    dynamics.simulation = simulation
    dynamics.integrator = integrator
    context = simulation.context
    configuration['platform'] = context.getPlatform().getName()
    configuration['seed'] = seed
    if restart_file is not None:
        simulation.loadState(str(restart_file))
    else:
        context.setPositions(dynamics.positions)
        if minimize:
            simulation.minimizeEnergy()
        if initial_velocities is not None:
            velocity = np.asarray(initial_velocities, dtype=float)
            if velocity.shape != (len(dynamics._potential_masses), 3) or not np.all(np.isfinite(velocity)):
                raise ValueError('Initial velocities must have shape (N, 3) in nm/ps.')
            force = context.getState(getForces=True).getForces(asNumpy=True).value_in_unit(
                unit.kilojoule_per_mole / unit.nanometer)
            stored = velocity - 0.5 * timestep * 0.001 * force / dynamics._potential_masses[:, None]
            context.setVelocities(stored * unit.nanometer / unit.picosecond)
        else:
            # OpenMM applies its integrator-specific thermal velocity offset.
            context.setVelocitiesToTemperature(dynamics.temperature, seed if seed is not None else 0)
    save_freq = max(nsteps // snapshots, 1) if snapshots else 0
    if save_freq:
        simulation.reporters.append(app.PDBReporter(str(traj_file), save_freq))
    for name in ('total_potentials', 'kinetic_energies', 'temperatures', 'total_energies'):
        setattr(dynamics, name, [])
    rows = []

    def observe(step):
        state = context.getState(getEnergy=True, getPositions=True,
                                 getVelocities=True, getForces=step_callback is not None)
        potential = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        kinetic = state.getKineticEnergy().value_in_unit(unit.kilojoule_per_mole)
        temp = 2 * kinetic / (3 * len(dynamics._potential_masses) *
                             unit.MOLAR_GAS_CONSTANT_R.value_in_unit(unit.kilojoule_per_mole / unit.kelvin))
        if not np.all(np.isfinite([potential, kinetic, temp])):
            raise RuntimeError('Nonfinite potential dynamics observable.')
        for name, value in [('total_potentials', potential), ('kinetic_energies', kinetic),
                            ('temperatures', temp), ('total_energies', potential + kinetic)]:
            getattr(dynamics, name).append(value)
        rows.append([state.getTime().value_in_unit(unit.picosecond), potential, kinetic, potential + kinetic, temp])
        if step_callback is not None:
            step_callback(step, state)

    start = perf_counter()
    observe(0)
    for step in range(1, nsteps + 1):
        simulation.step(1)
        observe(step)
    dynamics.wall_seconds = perf_counter() - start
    simulation.saveState(str(state_file))
    configuration.update(nsteps=nsteps, wall_seconds=dynamics.wall_seconds,
                         callback_calls=dynamics._potential_bridge.call_count,
                         final_time_ps=rows[-1][0])
    Path(str(state_file) + '.json').write_text(json.dumps(configuration, indent=2) + '\n')
    np.savetxt(str(output_file) + '.csv', rows, delimiter=',',
               header='time_ps,potential_kJ_mol,kinetic_kJ_mol,total_kJ_mol,temperature_K', comments='')
    Path(str(output_file) + '.json').write_text(json.dumps(configuration, indent=2) + '\n')
    if save_last_frame:
        with open(str(output_file) + '_last.pdb', 'w') as stream:
            app.PDBFile.writeFile(dynamics.topology, context.getState(getPositions=True).getPositions(), stream)
    simulation.reporters.clear()
    dynamics.ostream.print_info(f'{dynamics.driver_flag} {configuration["model_id"]}: '
                               f'{nsteps} steps on {configuration["platform"]}; '
                               f'physical energies in {output_file}.csv')
    dynamics.ostream.flush()
    return configuration
