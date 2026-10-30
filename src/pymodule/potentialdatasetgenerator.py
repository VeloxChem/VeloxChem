# SPDX-License-Identifier: BSD-3-Clause
"""Sampling, selection and reference labeling for molecular potentials."""

import json
from pathlib import Path
import sys
from time import perf_counter
from types import SimpleNamespace

import numpy as np
from mpi4py import MPI

from .inputparser import parse_input
from .outputstream import OutputStream
from .potential import PotentialResult
from .potentialdataset import PotentialDataset, _atomic_json, _atomic_npz, _atomic_text
from .referencepotential import ReferencePotential, _identity, _json_value, _versions
from .veloxchemlib import mpi_master


def _dihedral(positions, atoms):
    i, j, k, l = atoms
    axis = positions[k] - positions[j]
    axis /= np.linalg.norm(axis)
    first, last = positions[i] - positions[j], positions[l] - positions[k]
    first -= np.dot(first, axis) * axis
    last -= np.dot(last, axis) * axis
    if min(np.linalg.norm(first), np.linalg.norm(last)) < 1.e-12:
        raise ValueError('Cannot define a torsion from collinear atoms.')
    return np.arctan2(np.dot(np.cross(axis, first), last), np.dot(first, last))


def _rotors(connectivity, labels):
    """Find heavy-atom graph bridges with heavy neighbors on both sides."""
    connectivity = np.asarray(connectivity, dtype=bool)
    count = len(labels)

    def component(j, k):
        visited, pending = {k}, [k]
        while pending:
            atom = pending.pop()
            for neighbor in np.flatnonzero(connectivity[atom]):
                if {atom, neighbor} == {j, k}:
                    continue
                if neighbor not in visited:
                    visited.add(int(neighbor))
                    pending.append(int(neighbor))
        return visited

    candidates = []
    for j in range(count):
        for k in range(j + 1, count):
            if not connectivity[j, k] or 'H' in (labels[j], labels[k]):
                continue
            side = component(j, k)
            first = [i for i in np.flatnonzero(connectivity[j])
                     if i != k and labels[i] != 'H']
            last = [i for i in np.flatnonzero(connectivity[k])
                    if i != j and labels[i] != 'H']
            if j not in side and first and last:
                candidates.append(([int(first[0]), j, k, int(last[0])], sorted(side)))
    return candidates


def _rotor(connectivity, labels, central_bond=None):
    """Select one graph bridge; public bond indices are one-based."""
    candidates = _rotors(connectivity, labels)
    available = [[item[0][1] + 1, item[0][2] + 1] for item in candidates]
    if central_bond is not None:
        bond = sorted(int(i) - 1 for i in central_bond)
        candidates = [item for item in candidates if item[0][1:3] == bond]
    if len(candidates) != 1:
        raise ValueError('torsion_thermal needs one unambiguous rotatable bridge. '
                         f'Heavy-atom bridge candidates (one-based): {available}. '
                         'Set central_bond to one pair, central_bonds to a list of pairs, '
                         'or use recipe="thermal". Graph bridges alone do not exclude amides.')
    return candidates[0]


def _rotate(positions, atoms, side, angle):
    positions = np.array(positions, dtype=float, copy=True)
    j, k = atoms[1:3]
    delta = angle - _dihedral(positions, atoms)
    axis = positions[k] - positions[j]
    axis /= np.linalg.norm(axis)
    relative = positions[side] - positions[j]
    positions[side] = (relative * np.cos(delta) + np.cross(axis, relative) * np.sin(delta)
                       + np.outer(relative @ axis, axis) * (1 - np.cos(delta)) + positions[j])
    return positions


def _descriptor(positions):
    """All labeled pair distances: rotation/translation invariant, in Angstrom."""
    i, j = np.triu_indices(len(positions), 1)
    return np.linalg.norm(positions[i] - positions[j], axis=1)


def _select(candidates, budget, distance, seed, existing_counts=None):
    """Farthest-point selection with split quotas and global near-duplicate removal."""
    existing_counts = existing_counts or dict(train=0, valid=0, test=0)
    total = budget + sum(existing_counts.values())
    if budget < 0 or total < 3:
        raise ValueError('At least three reference evaluations are needed for three splits.')
    valid = max(1, int(round(total * .1)))
    quotas = dict(train=total - 2 * valid, valid=valid, test=valid)
    quotas = {split: max(0, count - existing_counts[split]) for split, count in quotas.items()}
    descriptors = np.array([_descriptor(item['positions']) for item in candidates])
    if not len(descriptors):
        raise ValueError('No valid sampling candidates.')
    nearest = np.full(len(candidates), np.inf)
    available = np.ones(len(candidates), dtype=bool)
    selected = []
    rng = np.random.default_rng(seed)
    # Independently generated test/validation campaigns get first choice.
    for split in ('test', 'valid', 'train'):
        for _ in range(quotas[split]):
            choices = np.array([i for i, frame in enumerate(candidates)
                                if available[i] and frame['split'] == split], dtype=int)
            if not len(choices):
                break
            index = int(rng.choice(choices) if not selected else choices[np.argmax(nearest[choices])])
            selected.append(index)
            delta = np.sqrt(np.mean((descriptors - descriptors[index]) ** 2, axis=1))
            nearest = np.minimum(nearest, delta)
            available &= delta > distance
    if any(existing_counts[split] == 0 and
           not any(candidates[i]['split'] == split for i in selected) for split in quotas):
        raise ValueError('Selection left an empty split; increase sampling or reduce duplicate_distance.')
    return [candidates[i] for i in selected]


class PotentialDatasetGenerator:
    """Generate a reference dataset for one isolated molecule, on one MPI rank.

    ``compute(molecule)`` samples gas-phase GAFF MM trajectories, selects
    diverse geometries, labels with the supplied xTB or SCF driver, then
    writes HDF5 and MACE inputs. ``compute(molecule, candidates=...)`` accepts
    user-sampled dictionaries containing positions (Angstrom), source with
    group_id, and split (train/valid/test), bypassing automatic sampling.

    ``filename`` names an output directory. Settings are frozen on first use;
    rerunning identical settings resumes committed data and reference attempts.
    max_reference_evaluations is a hard target evaluation budget, including
    failures and interrupted calls. MM setup charges are accounted separately
    unless they have the same reference target, in which case they count too.
    """

    _dataset_keywords = dict(filename='str', max_reference_evaluations='int', seed='int',
                             basis_set_label='str', duplicate_distance='float',
                             minimum_distance='float', maximum_bond_ratio='float',
                             max_retries='int')
    _sampling_keywords = dict(method='str_lower', recipe='str_lower', temperatures='seq_fixed',
                              timestep='float', nsteps='int', equilibration_steps='int',
                              sample_interval='int', torsion_spacing='float',
                              torsion_force_constant='float', central_bond='seq_fixed_int',
                              central_bonds='seq_fixed_int',
                              replicas='int', perturbations='int', perturbation_sigma='float',
                              friction='float', platform='str', partial_charges='seq_fixed')

    def __init__(self, reference_driver, comm=None, ostream=None):
        self.comm = MPI.COMM_WORLD if comm is None else comm
        self.rank = self.comm.Get_rank()
        self.nodes = self.comm.Get_size()
        if self.nodes != 1:
            raise ValueError('PotentialDatasetGenerator supports one MPI rank only.')
        self.ostream = ostream if ostream is not None else OutputStream(
            sys.stdout if self.rank == mpi_master() else None)
        self.reference_driver = reference_driver
        self.filename = 'potential_dataset'
        self.max_reference_evaluations = 1000
        self.seed = 123
        self.basis_set_label = None
        self.duplicate_distance = .01
        self.minimum_distance = .55
        self.maximum_bond_ratio = 1.6
        self.max_retries = 0
        self.sampling = SimpleNamespace(method='mm', recipe='torsion_thermal',
                                        temperatures=[300., 500.], timestep=.5,
                                        nsteps=2000, equilibration_steps=200,
                                        sample_interval=20, torsion_spacing=30.,
                                        torsion_force_constant=100., central_bond=None,
                                        central_bonds=None,
                                        replicas=1, perturbations=2, perturbation_sigma=.015,
                                        friction=1., platform='Reference', partial_charges=None)

    def update_settings(self, dataset_dict, sampling_dict=None):
        """Update typed settings; unknown keywords are errors, not ignored input."""
        for obj, keywords, values in ((self, self._dataset_keywords, dataset_dict),
                                      (self.sampling, self._sampling_keywords, sampling_dict or {})):
            unknown = set(values) - set(keywords)
            if unknown:
                raise ValueError(f'Unknown dataset/sampling settings: {sorted(unknown)}')
            parse_input(obj, keywords, values)
        self._validate_settings()

    def _validate_settings(self):
        s = self.sampling
        for name, value, minimum in [('max_reference_evaluations', self.max_reference_evaluations, 3),
                                     ('seed', self.seed, 1), ('max_retries', self.max_retries, 0),
                                     ('nsteps', s.nsteps, 1), ('sample_interval', s.sample_interval, 1),
                                     ('equilibration_steps', s.equilibration_steps, 0),
                                     ('replicas', s.replicas, 1), ('perturbations', s.perturbations, 0)]:
            if not isinstance(value, int) or value < minimum:
                raise ValueError(f'{name} must be an integer >= {minimum}.')
        values = [self.minimum_distance, self.maximum_bond_ratio, s.timestep,
                  s.torsion_spacing, s.torsion_force_constant, s.friction, *s.temperatures]
        if not s.temperatures or not np.all(np.isfinite(values)) or min(values) <= 0:
            raise ValueError('Sampling scales and temperatures must be positive and finite.')
        if any(not np.isfinite(x) or x < 0 for x in (self.duplicate_distance, s.perturbation_sigma)):
            raise ValueError('Duplicate distance and perturbation sigma must be nonnegative.')
        if s.method != 'mm' or s.recipe not in ('torsion_thermal', 'thermal'):
            raise ValueError('Automatic sampling supports MM thermal or torsion_thermal. '
                             'Use candidates/recorders for other sampling drivers.')
        if s.torsion_spacing > 360 or s.sample_interval > s.nsteps:
            raise ValueError('torsion_spacing must be <= 360 and sample_interval <= nsteps.')
        if s.nsteps % s.sample_interval:
            raise ValueError('nsteps must be a multiple of sample_interval.')
        if len(set(s.temperatures)) != len(s.temperatures):
            raise ValueError('List each sampling temperature only once.')
        if s.central_bond is not None and (len(s.central_bond) != 2 or
                                           len(set(s.central_bond)) != 2 or min(s.central_bond) < 1):
            raise ValueError('central_bond needs two distinct positive one-based atom indices.')
        if s.central_bonds is not None:
            if s.central_bond is not None or s.recipe != 'torsion_thermal':
                raise ValueError('central_bonds requires torsion_thermal and no central_bond.')
            bonds = np.asarray(s.central_bonds)
            if (bonds.ndim != 2 or bonds.shape[0] < 1 or bonds.shape[1] != 2 or
                    not np.issubdtype(bonds.dtype, np.integer) or np.any(bonds < 1) or
                    np.any(bonds[:, 0] == bonds[:, 1])):
                raise ValueError('central_bonds needs a list of distinct positive integer pairs.')
            canonical = [tuple(sorted(pair)) for pair in bonds.tolist()]
            if len(set(canonical)) != len(canonical):
                raise ValueError('central_bonds contains a repeated bond.')

    def _sampling_settings(self):
        """Keep old single-rotor campaign identities unchanged when no list is supplied."""
        settings = dict(vars(self.sampling))
        if settings['central_bonds'] is None:
            settings.pop('central_bonds')
        return settings

    def _reference(self, molecule):
        if isinstance(self.reference_driver, ReferencePotential):
            return self.reference_driver
        return ReferencePotential(self.reference_driver, molecule, self.basis_set_label)

    def _freeze(self, path, settings):
        settings = _json_value(settings)
        if path.exists():
            if json.loads(path.read_text()) != settings:
                raise ValueError(f'Campaign settings changed: {path}. Choose a new filename.')
        else:
            _atomic_json(path, settings)

    def compute(self, molecule, candidates=None):
        """Create/resume the dataset and return paths, counts and coverage diagnostics."""
        self._validate_settings()
        start = perf_counter()
        directory = Path(self.filename).expanduser().resolve()
        reference = self._reference(molecule)
        if reference.comm.Get_size() != 1:
            raise ValueError('Reference evaluation supports one MPI rank only.')
        if candidates is not None:
            candidates = list(candidates)
        with PotentialDataset(directory / 'dataset.h5', molecule) as dataset:
            settings = dict(dataset={key: getattr(self, key) for key in self._dataset_keywords},
                            sampling=self._sampling_settings(), target_id=reference.model_id,
                            initial_geometry_id=dataset.geometry_id(molecule.get_coordinates_in_angstrom()),
                            source='provided candidates' if candidates is not None else 'automatic MM',
                            candidates_sha256=None if candidates is None else _identity(candidates),
                            versions=_versions())
            self._freeze(directory / 'resolved_settings.json', settings)
            target_id = dataset.add_target(reference.manifest)
            selection_path = directory / 'selection.json'
            if selection_path.exists():
                selected_ids = json.loads(selection_path.read_text())
            else:
                if candidates is None:
                    candidates = self._sample(molecule, dataset, reference, directory)
                accepted, rejected = self._filter(molecule, candidates)
                # Existing setup labels count toward the target budget and enter
                # the training set. Exclude their exact/near copies in all arms.
                existing = dataset.labeled_ids(target_id)
                for frame_id in existing:
                    descriptor = _descriptor(dataset.frame(frame_id)['positions'])
                    accepted = [item for item in accepted if np.sqrt(np.mean(
                        (_descriptor(item['positions']) - descriptor) ** 2)) > self.duplicate_distance]
                remaining = self.max_reference_evaluations - self._attempt_count(dataset, target_id)
                existing_counts = {split: sum(dataset.frame(i)['split'] == split for i in existing)
                                   for split in ('train', 'valid', 'test')}
                selected = _select(accepted, remaining, self.duplicate_distance, self.seed, existing_counts)
                selected_ids = existing + [dataset.add_candidate(item['positions'], item['source'],
                                                                item['split']) for item in selected]
                _atomic_json(selection_path, selected_ids)
                _atomic_json(directory / 'selection_report.json', dict(accepted_candidates=len(accepted),
                              rejected=rejected, selected=len(selected_ids),
                              duplicate_metric='RMS of labeled pair-distance differences, Angstrom',
                              duplicate_distance=self.duplicate_distance))
            self._label(dataset, reference, selected_ids, directory)
            results = self._finalize(dataset, reference, directory)
        results['wall_seconds'] = perf_counter() - start
        _atomic_json(directory / 'results.json', results)
        self.ostream.print_info(f'Dataset: {results["dataset_file"]}; '
                               f'{results["successful_labels"]} labels, '
                               f'{results["reference_attempts"]} reference attempts.')
        self.ostream.flush()
        return results

    @staticmethod
    def _attempt_count(dataset, target_id):
        return len(list((dataset.records / 'attempts' / target_id).glob('*.json')))

    def _label(self, dataset, reference, selected_ids, directory):
        target_id = reference.model_id
        complete = set(dataset.labeled_ids(target_id))
        failures = []
        for index, frame_id in enumerate(selected_ids):
            if frame_id in complete:
                continue
            previous = [json.loads(p.read_text()) for p in
                        (dataset.records / 'attempts' / target_id).glob('*.json')]
            tries = sum(record['frame_id'] == frame_id for record in previous)
            while tries <= self.max_retries:
                if self._attempt_count(dataset, target_id) >= self.max_reference_evaluations:
                    break
                attempt = dataset.start_attempt(frame_id, target_id, self.max_reference_evaluations)
                started = perf_counter()
                try:
                    frame = dataset.frame(frame_id)
                    result = reference.compute(dataset.molecule(frame['positions']))
                    dataset.add_result(frame_id, result)
                except Exception as error:
                    failure = dict(frame_id=frame_id, error_type=type(error).__name__, error=str(error))
                    failures.append(failure)
                    dataset.finish_attempt(attempt, 'failed', wall_seconds=perf_counter()-started,
                                           **failure)
                    tries += 1
                else:
                    dataset.finish_attempt(attempt, 'complete', wall_seconds=perf_counter()-started)
                    complete.add(frame_id)
                    break
            if index % 20 == 0:
                self.ostream.print_info(f'Reference labels: {len(complete)}/{len(selected_ids)}')
                self.ostream.flush()
        records = [json.loads(p.read_text()) for p in
                   sorted((dataset.records / 'attempts' / target_id).glob('*.json'))]
        _atomic_text(directory / 'failures.jsonl', ''.join(json.dumps(r) + '\n' for r in records
                                                          if r['status'] != 'complete'))

    def _filter(self, molecule, candidates):
        reference_distances = np.linalg.norm(molecule.get_coordinates_in_angstrom()[:, None] -
                                             molecule.get_coordinates_in_angstrom()[None, :], axis=2)
        bonded = np.asarray(molecule.get_connectivity_matrix(), dtype=bool)
        accepted, rejected, seen = [], dict(nonfinite=0, collision=0, stretched_bond=0, exact_duplicate=0), set()
        groups = {}
        for candidate in candidates:
            item = dict(candidate)
            positions = np.array(item['positions'], dtype=float, copy=True)
            if positions.shape != (molecule.number_of_atoms(), 3) or not np.all(np.isfinite(positions)):
                rejected['nonfinite'] += 1
                continue
            distances = np.linalg.norm(positions[:, None] - positions[None, :], axis=2)
            np.fill_diagonal(distances, np.inf)
            if np.min(distances) < self.minimum_distance:
                rejected['collision'] += 1
                continue
            if np.any(distances[bonded] > self.maximum_bond_ratio * reference_distances[bonded]):
                rejected['stretched_bond'] += 1
                continue
            identity = positions.astype('<f8').tobytes()
            if identity in seen:
                rejected['exact_duplicate'] += 1
                continue
            seen.add(identity)
            split = item['split']
            group = item['source']['group_id']
            if split not in ('train', 'valid', 'test') or not isinstance(group, str) or not group:
                raise ValueError('Candidates require source.group_id and train/valid/test split.')
            if group in groups and groups[group] != split:
                raise ValueError('A sampling group cannot cross splits.')
            groups[group] = split
            item['positions'] = positions
            accepted.append(item)
        return accepted, rejected

    def _sample(self, molecule, dataset, reference, directory):
        """Native OpenMM MM windows; restraining energy never enters reference labels."""
        import openmm as mm
        from openmm import app, unit
        from .mmforcefieldgenerator import MMForceFieldGenerator
        from .xtbdriver import XtbDriver

        s = self.sampling
        sample_dir = directory / 'sampling'
        sample_dir.mkdir(exist_ok=True)
        charges_path = sample_dir / 'charges.json'
        if charges_path.exists():
            charges = json.loads(charges_path.read_text())['partial_charges']
        elif s.partial_charges is not None:
            charges = s.partial_charges
            _atomic_json(charges_path, dict(partial_charges=charges, source='user supplied'))
        else:
            charge_driver = XtbDriver(self.comm, self.ostream)
            charge_driver.set_method('gfn2')
            charge_reference = ReferencePotential(charge_driver, molecule)
            same_target = charge_reference.model_id == reference.model_id
            frame_id = dataset.add_candidate(molecule.get_coordinates_in_angstrom(),
                                            dict(group_id='initial', source='MM charge setup'), 'train')
            attempt = (dataset.start_attempt(frame_id, reference.model_id, self.max_reference_evaluations)
                       if same_target else None)
            try:
                values = charge_driver.compute(molecule)
                charges = values['partial_charges'].tolist()
                if same_target:
                    dataset.add_result(frame_id, PotentialResult(values['energy'], values['gradient'],
                                                                reference.model_id))
                    dataset.finish_attempt(attempt, 'complete', purpose='MM setup and reference label')
                _atomic_json(charges_path, dict(partial_charges=charges, source=charge_reference.manifest,
                                               separate_setup_evaluations=0 if same_target else 1))
            except Exception as error:
                if attempt is not None:
                    dataset.finish_attempt(attempt, 'failed', error=str(error))
                raise
        charges = np.asarray(charges, dtype=float)
        if charges.shape != (molecule.number_of_atoms(),) or not np.all(np.isfinite(charges)):
            raise ValueError('MM sampling needs one finite partial charge per atom.')
        if abs(charges.sum() - molecule.get_charge()) > 1.e-4:
            raise ValueError('MM partial charges do not sum to the molecular charge.')
        ff = MMForceFieldGenerator(self.comm, self.ostream)
        ff.partial_charges = charges.copy()
        ff.create_topology(molecule)  # preprovided charges bypass RESP; resp=False would zero them
        ff.write_openmm_files(str(sample_dir / 'molecule'), 'MOL')
        pdb = app.PDBFile(str(sample_dir / 'molecule.pdb'))
        forcefield = app.ForceField(str(sample_dir / 'molecule.xml'))
        system = forcefield.createSystem(pdb.topology, nonbondedMethod=app.NoCutoff,
                                         constraints=None, rigidWater=False, removeCMMotion=False)
        rotors = []
        windows = [(None, 0, None)]
        multiple = s.central_bonds is not None
        if s.recipe == 'torsion_thermal':
            bonds = s.central_bonds if multiple else [s.central_bond]
            rotors = [_rotor(ff.connectivity_matrix, molecule.get_labels(), bond) for bond in bonds]
            if multiple:
                native_bonds = {tuple(sorted(bond)) for bond in ff.rotatable_bonds}
                for atoms, _ in rotors:
                    bond = tuple(i + 1 for i in atoms[1:3])
                    if bond not in native_bonds:
                        raise ValueError(f'Bond {bond} is excluded by the MM rotatable-bond rules.')
            angles = np.deg2rad(np.arange(-180., 180., s.torsion_spacing)).tolist()
            windows = [(index, window, angle) for index in range(len(rotors))
                       for window, angle in enumerate(angles)]
            for index, (atoms, _) in enumerate(rotors):
                strength = f'k_r{index}' if multiple else 'k'
                center = f'theta0_r{index}' if multiple else 'theta0'
                bias = mm.CustomTorsionForce(f'{strength}*(1-cos(theta-{center}))')
                bias.addGlobalParameter(strength, 0. if multiple else s.torsion_force_constant)
                bias.addGlobalParameter(center, 0.)
                bias.addTorsion(*atoms, [])
                system.addForce(bias)
        manifest = dict(method='GAFF with supplied/GFN2-xTB charges',
                     charge_record=str(charges_path), platform=s.platform, settings=self._sampling_settings(),
                     torsion_atoms_one_based=([i+1 for i in rotors[0][0]]
                                              if rotors and not multiple else None),
                     system_sha256=_identity(mm.XmlSerializer.serialize(system)),
                     note='Independent seeds for each split/window; reference labels exclude bias.')
        if multiple:
            manifest.update(torsions_atoms_one_based=[[i+1 for i in atoms] for atoms, _ in rotors],
                            initialization='Independent uniform initial angles for all selected rotors; '
                                           'one restrained rotor per window, others unrestrained.')
        self._freeze(sample_dir / 'manifest.json', manifest)
        platform = mm.Platform.getPlatformByName(s.platform)
        rng = np.random.default_rng(self.seed)
        candidates = []
        initial = np.asarray(molecule.get_coordinates_in_angstrom())
        for split in ('train', 'valid', 'test'):
            for replica in range(s.replicas):
                for rotor_index, window, angle in windows:
                    for temperature_index, temperature in enumerate(s.temperatures):
                        rotor_tag = f'_b{rotor_index}' if multiple else ''
                        group_id = f'{split}_r{replica}{rotor_tag}_w{window}_t{temperature_index}'
                        trajectory_seed = int(rng.integers(1, 2**31-1))
                        trajectory_rng = np.random.default_rng(trajectory_seed)
                        path = sample_dir / (group_id + '.npz')
                        if path.exists():
                            with np.load(path, allow_pickle=False) as data:
                                positions = data['positions'].copy()
                                sources = json.loads(str(data['sources']))
                        else:
                            integrator = mm.LangevinIntegrator(temperature * unit.kelvin,
                                                               s.friction / unit.picosecond,
                                                               s.timestep * unit.femtosecond)
                            integrator.setRandomNumberSeed(trajectory_seed)
                            properties = {'Threads': '1'} if s.platform == 'CPU' else {}
                            simulation = app.Simulation(pdb.topology, system, integrator, platform, properties)
                            initial_angles = None
                            xyz = initial
                            if multiple:
                                initial_angles = trajectory_rng.uniform(-np.pi, np.pi, len(rotors)).tolist()
                                initial_angles[rotor_index] = angle
                                for (atoms, side), initial_angle in zip(rotors, initial_angles):
                                    xyz = _rotate(xyz, atoms, side, initial_angle)
                            elif angle is not None:
                                atoms, side = rotors[rotor_index]
                                xyz = _rotate(initial, atoms, side, angle)
                            simulation.context.setPositions(xyz * .1 * unit.nanometer)
                            if angle is not None:
                                center = f'theta0_r{rotor_index}' if multiple else 'theta0'
                                simulation.context.setParameter(center, angle)
                                if multiple:
                                    simulation.context.setParameter(f'k_r{rotor_index}', s.torsion_force_constant)
                            simulation.minimizeEnergy(maxIterations=500)
                            simulation.context.setVelocitiesToTemperature(temperature, trajectory_seed)
                            if s.equilibration_steps:
                                simulation.step(s.equilibration_steps)
                            positions, sources = [], []
                            for step in range(s.sample_interval, s.nsteps + 1, s.sample_interval):
                                simulation.step(s.sample_interval)
                                state = simulation.context.getState(getPositions=True, getEnergy=True)
                                xyz = state.getPositions(asNumpy=True).value_in_unit(unit.nanometer) * 10
                                for perturbation in range(s.perturbations + 1):
                                    displacement = (np.zeros_like(xyz) if perturbation == 0 else
                                                    trajectory_rng.normal(0., s.perturbation_sigma, xyz.shape))
                                    positions.append(xyz + displacement)
                                    source = dict(group_id=group_id, source='restrained MM' if rotors else 'MM',
                                                        seed=trajectory_seed, step=step + s.equilibration_steps,
                                                        time_ps=(step+s.equilibration_steps)*s.timestep*.001,
                                                        temperature_K=temperature, window_radians=angle,
                                                        perturbation=perturbation)
                                    if multiple:
                                        source.update(torsion_index=rotor_index,
                                                      central_bond_one_based=[i+1 for i in rotors[rotor_index][0][1:3]],
                                                      initial_torsions_radians=initial_angles)
                                    sources.append(source)
                            _atomic_npz(path, positions=np.array(positions), sources=np.array(json.dumps(sources)))
                            del simulation, integrator
                        candidates.extend(dict(positions=xyz, source=source, split=split)
                                          for xyz, source in zip(positions, sources))
            self.ostream.print_info(f'MM sampling completed for {split}; {len(candidates)} candidates so far.')
            self.ostream.flush()
        return candidates

    def _finalize(self, dataset, reference, directory):
        dataset.consolidate()
        target_id = reference.model_id
        ids = dataset.labeled_ids(target_id)
        split_counts = {split: sum(dataset.frame(i)['split'] == split for i in ids)
                        for split in ('train', 'valid', 'test')}
        energies = [dataset.result(i, target_id).energy for i in ids]
        gradient_rms = [float(np.sqrt(np.mean(dataset.result(i, target_id).gradient**2))) for i in ids]
        coverage = dict(split_counts=split_counts, energy_range_hartree=[min(energies), max(energies)] if ids else [],
                        gradient_rms_range_hartree_bohr=[min(gradient_rms), max(gradient_rms)] if ids else [],
                        groups=sorted({dataset.frame(i)['source']['group_id'] for i in ids}),
                        scientific_accuracy='not certified by dataset generation')
        sampling_manifest = directory / 'sampling' / 'manifest.json'
        if sampling_manifest.exists():
            sampling = json.loads(sampling_manifest.read_text())
            torsion = sampling['torsion_atoms_one_based']
            if torsion is not None:
                atoms = [i - 1 for i in torsion]
                coverage['torsion_atoms_one_based'] = torsion
                coverage['torsion_histograms'] = {}
                for split in split_counts:
                    angles = [np.rad2deg(_dihedral(dataset.frame(i)['positions'], atoms))
                              for i in ids if dataset.frame(i)['split'] == split]
                    counts, edges = np.histogram(angles, bins=np.linspace(-180., 180., 13))
                    coverage['torsion_histograms'][split] = dict(counts=counts.tolist(),
                                                                edges_degrees=edges.tolist())
            if sampling.get('torsions_atoms_one_based'):
                coverage['torsion_histograms_by_bond'] = {}
                for torsion in sampling['torsions_atoms_one_based']:
                    atoms = [i - 1 for i in torsion]
                    histograms = dict(torsion_atoms_one_based=torsion)
                    for split in split_counts:
                        angles = [np.rad2deg(_dihedral(dataset.frame(i)['positions'], atoms))
                                  for i in ids if dataset.frame(i)['split'] == split]
                        counts, edges = np.histogram(angles, bins=np.linspace(-180., 180., 13))
                        histograms[split] = dict(counts=counts.tolist(), edges_degrees=edges.tolist())
                    coverage['torsion_histograms_by_bond'][f'{torsion[1]}-{torsion[2]}'] = histograms
        _atomic_json(directory / 'coverage.json', coverage)
        files = dataset.export_mace(target_id, directory)
        config = self.write_mace_config(files, directory, seed=self.seed)
        return dict(dataset_file=str(dataset.filename), target_id=target_id,
                    training_file=files['train'], validation_file=files['valid'], test_file=files['test'],
                    mace_config_file=config, successful_labels=len(ids), split_counts=split_counts,
                    reference_attempts=self._attempt_count(dataset, target_id),
                    coverage_file=str(directory / 'coverage.json'))

    @staticmethod
    def write_mace_config(files, directory, seed=123, epochs=200):
        """Write a MACE YAML input; training is a separate explicit action."""
        directory = Path(directory).resolve()
        training = directory / 'training'
        config = dict(name='molecular_potential', train_file=files['train'], valid_file=files['valid'],
                      test_file=files['test'], energy_key='REF_energy', forces_key='REF_forces',
                      E0s='average', model='ScaleShiftMACE', hidden_irreps='32x0e', num_interactions=2,
                      correlation=3, r_max=5.0, default_dtype='float64', device='cpu', seed=seed,
                      batch_size=10, valid_batch_size=10, max_num_epochs=epochs,
                      ema=True, ema_decay=.99, amsgrad=True, num_workers=0, save_cpu=True,
                      log_dir=str(training / 'logs'), model_dir=str(training / 'models'),
                      checkpoints_dir=str(training / 'checkpoints'), results_dir=str(training / 'results'))
        if epochs >= 4:
            config.update(swa=True, start_swa=int(.75*epochs))
        path = directory / 'mace.yaml'
        _atomic_text(path, '# Fixed composition: average E0s are fitted offsets, not isolated-atom energies.\n' +
                     '\n'.join(f'{key}: {json.dumps(value)}' for key, value in config.items()) + '\n')
        return str(path)

    def label_existing(self, filename):
        """Add a new reference target to existing geometries, preserving all splits/labels."""
        self._validate_settings()
        directory = Path(self.filename).expanduser().resolve()
        if directory == Path(filename).expanduser().resolve().parent:
            raise ValueError('Choose a separate output directory for the new target exports.')
        with PotentialDataset(filename) as dataset:
            ids = dataset.frame_ids()
            if not ids:
                raise ValueError('Existing dataset has no frames.')
            molecule = dataset.molecule(dataset.frame(ids[0])['positions'])
            reference = self._reference(molecule)
            dataset.add_target(reference.manifest)
            # Round-robin across source splits avoids exhausting a small relabel
            # budget in one split just because geometry hashes sort that way.
            grouped = {split: [i for i in ids if dataset.frame(i)['split'] == split]
                       for split in ('train', 'valid', 'test')}
            if sum(map(len, grouped.values())) != len(ids):
                raise ValueError('Assign source splits before relabeling.')
            ids = [grouped[split][index] for index in range(max(map(len, grouped.values())))
                   for split in ('train', 'valid', 'test') if index < len(grouped[split])]
            self._freeze(directory / 'relabel_settings.json', dict(target_id=reference.model_id,
                         source=str(dataset.filename), frames=ids,
                         budget=self.max_reference_evaluations, max_retries=self.max_retries))
            self._label(dataset, reference, ids, directory)
            return self._finalize(dataset, reference, directory)
