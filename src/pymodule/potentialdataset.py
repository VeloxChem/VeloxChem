# SPDX-License-Identifier: BSD-3-Clause
"""Recoverable isolated-molecule energy/gradient datasets and MACE export.

Closed, atomically renamed NumPy shards are the recovery journal. HDF5 is a
portable, consolidated snapshot; it is never opened in truncate-on-error
mode. One writer holds an advisory file lock for the dataset lifetime.
"""

import copy
import hashlib
import json
import os
from pathlib import Path
import tempfile

import numpy as np

from .molecule import Molecule
from .potential import PotentialResult
from .referencepotential import _identity, _json_value
from .veloxchemlib import bohr_in_angstrom, hartree_in_ev


def _atomic_json(path, value):
    _atomic_text(path, json.dumps(_json_value(value), indent=2, sort_keys=True,
                                 allow_nan=False) + '\n')


def _atomic_text(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix='.' + path.name, dir=path.parent)
    try:
        with os.fdopen(fd, 'w', encoding='utf-8') as stream:
            stream.write(value)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def _atomic_npz(path, **arrays):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix='.' + path.name, dir=path.parent)
    try:
        with os.fdopen(fd, 'wb') as stream:
            np.savez_compressed(stream, **arrays)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


class PotentialDataset:
    """Store one fixed molecular system with multiple separate reference targets.

    ``filename`` is the HDF5 snapshot, normally ``run/dataset.h5``. Its sibling
    ``dataset.records`` directory retains recoverable frames, labels and
    attempt records. Keep both when moving a resumable campaign. Open an
    existing HDF5 file without shards for read-only inspection/export.

    Use ``with PotentialDataset(filename, molecule) as dataset:`` for writes.
    Reads do not take the writer lock. Atom order and electronic state are
    invariant; geometry hashes use exact float64 Cartesian Angstrom bytes.
    """

    def __init__(self, filename, molecule=None):
        self.filename = Path(filename).expanduser().resolve()
        self.records = self.filename.with_suffix('.records')
        self._lock = None
        manifest = self.records / 'system.json'
        if manifest.exists():
            self.system = json.loads(manifest.read_text())
        elif self.filename.exists():
            import h5py
            with h5py.File(self.filename, 'r') as handle:
                if handle.attrs['format'] != 'VeloxChem potential dataset 1':
                    raise ValueError('Unsupported potential dataset format.')
                self.system = json.loads(handle['system'].attrs['manifest'])
        elif molecule is not None:
            self.system = self._system(molecule)
        else:
            raise FileNotFoundError(self.filename)
        if molecule is not None and self._system(molecule) != self.system:
            raise ValueError('Dataset system differs in atom order, mass or electronic state.')

    @staticmethod
    def _system(molecule):
        return dict(labels=list(molecule.get_labels()),
                    atomic_numbers=[int(z) for z in molecule.get_element_ids()],
                    masses_dalton=np.asarray(molecule.get_masses()).tolist(),
                    charge=molecule.get_charge(), multiplicity=molecule.get_multiplicity(),
                    periodic=False)

    def __enter__(self):
        import fcntl
        if self._lock is not None:
            raise RuntimeError('Dataset writer is already open.')
        if self.filename.exists() and not (self.records / 'system.json').exists():
            raise ValueError('Portable HDF5 is read-only; restore its .records directory to resume.')
        self.records.mkdir(parents=True, exist_ok=True)
        self._lock = (self.records / '.writer.lock').open('a+')
        try:
            try:
                fcntl.flock(self._lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError as error:
                raise RuntimeError('Another writer holds this dataset.') from error
            path = self.records / 'system.json'
            if path.exists() and json.loads(path.read_text()) != self.system:
                raise ValueError('Concurrent dataset system mismatch.')
            if not path.exists():
                _atomic_json(path, self.system)
        except BaseException:
            self._lock.close()
            self._lock = None
            raise
        return self

    def __exit__(self, *exc):
        self._lock.close()
        self._lock = None

    def _require_writer(self):
        if self._lock is None:
            raise RuntimeError('Dataset mutation requires a with PotentialDataset(...) block.')

    def molecule(self, positions):
        molecule = Molecule(self.system['labels'], np.asarray(positions), 'angstrom')
        molecule.set_charge(self.system['charge'])
        molecule.set_multiplicity(self.system['multiplicity'])
        return molecule

    def geometry_id(self, positions):
        positions = np.array(positions, dtype='<f8', copy=True)
        if (positions.shape != (len(self.system['labels']), 3) or
                not np.all(np.isfinite(positions))):
            raise ValueError('Coordinates must be finite (N, 3) Angstrom values.')
        positions[positions == 0] = 0  # canonicalize signed zero
        return hashlib.sha256(_identity(self.system).encode() + positions.tobytes()).hexdigest()

    def add_candidate(self, positions, source, split=None):
        """Persist full-precision positions and provenance, without a reference call.

        ``source['group_id']`` identifies a trajectory/perturbation family.
        Optional split is train, valid or test and cannot later change.
        """
        self._require_writer()
        source = _json_value(source)
        if not isinstance(source.get('group_id'), str) or not source['group_id']:
            raise ValueError('Every candidate needs a nonempty source group_id.')
        if split not in (None, 'train', 'valid', 'test'):
            raise ValueError('split must be train, valid or test.')
        frame_id = self.geometry_id(positions)
        path = self.records / 'frames' / (frame_id + '.npz')
        if path.exists():
            old = self.frame(frame_id)
            if old['split'] != split:
                raise ValueError('Identical geometry encountered in different splits.')
            return frame_id
        metadata = dict(frame_id=frame_id, source=source, split=split)
        _atomic_npz(path, positions=np.asarray(positions, dtype=float),
                    metadata=np.array(json.dumps(metadata, allow_nan=False)))
        return frame_id

    def frame_ids(self):
        if (self.records / 'system.json').exists():
            return sorted(path.stem for path in (self.records / 'frames').glob('*.npz'))
        import h5py
        with h5py.File(self.filename, 'r') as handle:
            return sorted(handle['frames'])

    def frame(self, frame_id):
        if (not isinstance(frame_id, str) or len(frame_id) != 64 or
                any(c not in '0123456789abcdef' for c in frame_id)):
            raise KeyError(frame_id)
        if (self.records / 'system.json').exists():
            with np.load(self.records / 'frames' / (frame_id + '.npz'), allow_pickle=False) as data:
                frame = json.loads(str(data['metadata']))
                frame['positions'] = data['positions'].copy()
        else:
            import h5py
            with h5py.File(self.filename, 'r') as handle:
                group = handle['frames'][frame_id]
                frame = json.loads(group.attrs['metadata'])
                frame['positions'] = group['positions_angstrom'][()]
        if self.geometry_id(frame['positions']) != frame_id:
            raise ValueError('Stored geometry hash mismatch.')
        return frame

    def add_target(self, manifest):
        self._require_writer()
        target = copy.deepcopy(manifest['target'])
        target_id = _identity(target)
        if manifest.get('model_id', target_id) != target_id:
            raise ValueError('Reference identity must hash its complete target.')
        for key in ('charge', 'multiplicity'):
            if target[key] != self.system[key]:
                raise ValueError('Reference target does not match the molecular system.')
        if target.get('environment') != 'isolated molecule' or target.get('model_role') != 'full':
            raise ValueError('Only full isolated-molecule reference labels are supported.')
        path = self.records / 'targets' / (target_id + '.json')
        if not path.exists():
            _atomic_json(path, target)
        return target_id

    def target(self, target_id):
        if (self.records / 'system.json').exists():
            target = json.loads((self.records / 'targets' / (target_id + '.json')).read_text())
        else:
            import h5py
            with h5py.File(self.filename, 'r') as handle:
                target = json.loads(handle['targets'][target_id].attrs['manifest'])
        if _identity(target) != target_id:
            raise ValueError('Stored reference target hash mismatch.')
        return target

    def add_result(self, frame_id, result):
        """Commit a successful label bound to this exact geometry and target."""
        self._require_writer()
        self.target(result.model_id)
        self.frame(frame_id)
        result = PotentialResult(result.energy, result.gradient, result.model_id)
        if result.gradient.shape != (len(self.system['labels']), 3):
            raise ValueError('Label gradient does not match atom count.')
        path = self.records / 'labels' / result.model_id / (frame_id + '.npz')
        if path.exists():
            old = self.result(frame_id, result.model_id)
            if old.energy != result.energy or not np.array_equal(old.gradient, result.gradient):
                raise ValueError('A completed reference label is immutable.')
            return
        _atomic_npz(path, energy=result.energy, gradient=result.gradient,
                    frame_id=np.array(frame_id), target_id=np.array(result.model_id))

    def labeled_ids(self, target_id):
        self.target(target_id)
        if (self.records / 'system.json').exists():
            return sorted(path.stem for path in
                          (self.records / 'labels' / target_id).glob('*.npz'))
        import h5py
        with h5py.File(self.filename, 'r') as handle:
            return sorted(handle['labels'][target_id])

    def result(self, frame_id, target_id):
        if (self.records / 'system.json').exists():
            with np.load(self.records / 'labels' / target_id / (frame_id + '.npz'),
                         allow_pickle=False) as data:
                if str(data['frame_id']) != frame_id or str(data['target_id']) != target_id:
                    raise ValueError('Stored label identity mismatch.')
                return PotentialResult(float(data['energy']), data['gradient'], target_id)
        import h5py
        with h5py.File(self.filename, 'r') as handle:
            group = handle['labels'][target_id][frame_id]
            return PotentialResult(group['energy_hartree'][()],
                                   group['gradient_hartree_bohr'][()], target_id)

    def start_attempt(self, frame_id, target_id, budget):
        """Persist an attempt *before* evaluation; interrupted attempts consume budget."""
        self._require_writer()
        self.frame(frame_id)
        self.target(target_id)
        directory = self.records / 'attempts' / target_id
        attempts = sorted(directory.glob('*.json'))
        if len(attempts) >= budget:
            raise RuntimeError('Reference evaluation budget exhausted.')
        path = directory / f'{len(attempts):08d}.json'
        _atomic_json(path, dict(frame_id=frame_id, target_id=target_id, status='started'))
        return path

    def finish_attempt(self, path, status, **details):
        self._require_writer()
        record = json.loads(Path(path).read_text())
        record.update(status=status, **details)
        _atomic_json(path, record)

    def consolidate(self):
        """Write a closed HDF5 snapshot, then replace the previous snapshot."""
        import h5py
        self._require_writer()
        fd, temporary = tempfile.mkstemp(prefix='.dataset-', suffix='.h5',
                                         dir=self.filename.parent)
        os.close(fd)
        try:
            with h5py.File(temporary, 'w') as handle:
                handle.attrs['format'] = 'VeloxChem potential dataset 1'
                handle.create_group('system').attrs['manifest'] = json.dumps(self.system)
                frames = handle.create_group('frames')
                for frame_id in self.frame_ids():
                    frame = self.frame(frame_id)
                    group = frames.create_group(frame_id)
                    group.create_dataset('positions_angstrom', data=frame.pop('positions'))
                    group.attrs['metadata'] = json.dumps(frame)
                targets = handle.create_group('targets')
                labels = handle.create_group('labels')
                for path in sorted((self.records / 'targets').glob('*.json')):
                    target_id = path.stem
                    targets.create_group(target_id).attrs['manifest'] = json.dumps(self.target(target_id))
                    target_labels = labels.create_group(target_id)
                    for frame_id in self.labeled_ids(target_id):
                        result = self.result(frame_id, target_id)
                        group = target_labels.create_group(frame_id)
                        group.create_dataset('energy_hartree', data=result.energy)
                        group.create_dataset('gradient_hartree_bohr', data=result.gradient)
                attempts = [json.loads(p.read_text()) for p in
                            sorted((self.records / 'attempts').glob('*/*.json'))]
                handle.create_dataset('attempts_json', data=json.dumps(attempts))
            with open(temporary, 'rb') as stream:
                os.fsync(stream.fileno())
            os.replace(temporary, self.filename)
        finally:
            if os.path.exists(temporary):
                os.unlink(temporary)
        return str(self.filename)

    def export_mace(self, target_id, directory):
        """Export separate train/valid/test sets using REF_energy and REF_forces.

        Energies are full eV values; physical forces are -gradient in eV/A.
        Seventeen significant digits preserve float64 coordinates and forces.
        No stress or virial labels are fabricated for an isolated molecule.
        """
        directory = Path(directory).resolve()
        ids = self.labeled_ids(target_id)
        partitions = {key: [] for key in ('train', 'valid', 'test')}
        group_splits = {}
        for frame_id in ids:
            frame = self.frame(frame_id)
            split = frame['split']
            if split not in partitions:
                raise ValueError('Assign grouped splits before MACE export.')
            group = frame['source']['group_id']
            if group in group_splits and group_splits[group] != split:
                raise ValueError('A source group crosses dataset splits.')
            group_splits[group] = split
            partitions[split].append(frame_id)
        if any(not group for group in partitions.values()):
            raise ValueError('MACE export requires nonempty train, valid and test sets.')
        result = {}
        h, a = hartree_in_ev(), bohr_in_angstrom()
        for split, frame_ids in partitions.items():
            lines = []
            for frame_id in frame_ids:
                frame = self.frame(frame_id)
                label = self.result(frame_id, target_id)
                forces = -label.gradient * h / a
                info = dict(frame_id=frame_id, group_id=frame['source']['group_id'],
                            target_id=target_id, config_type='molecule')
                lines += [str(len(self.system['labels'])),
                          'Properties=species:S:1:pos:R:3:REF_forces:R:3 '
                          f'REF_energy={label.energy * h:.17g} pbc="F F F" ' +
                          ' '.join(f'{key}={json.dumps(value)}' for key, value in info.items())]
                for symbol, xyz, force in zip(self.system['labels'], frame['positions'], forces):
                    lines.append(symbol + ' ' + ' '.join(f'{x:.17g}' for x in (*xyz, *force)))
            path = directory / (split + '.extxyz')
            _atomic_text(path, '\n'.join(lines) + '\n')
            result[split] = str(path)
        _atomic_json(directory / 'export_manifest.json',
                     dict(target_id=target_id, target=self.target(target_id), splits=partitions,
                          units=dict(energy='eV', forces='eV/angstrom', positions='angstrom'),
                          files={key: dict(path=path, sha256=hashlib.sha256(
                              Path(path).read_bytes()).hexdigest()) for key, path in result.items()}))
        return result


class DatasetRecorder:
    """Full-precision candidate recorder for OpenMMDynamics.step_callback.

    This records positions, not the dynamics engine's possibly biased energy.
    Reference MD can instead call ``record_result`` from ReferencePotential's
    evaluation_callback to retain already computed physical E/G pairs.
    """

    def __init__(self, dataset, trajectory_id, interval=20, split=None):
        if not isinstance(interval, int) or interval < 1:
            raise ValueError('Recorder interval must be a positive integer.')
        self.dataset = dataset
        self.trajectory_id = trajectory_id
        self.interval = interval
        self.split = split
        self.evaluations = 0

    def __call__(self, step, state):
        if step % self.interval:
            return
        from openmm import unit
        positions = state.getPositions(asNumpy=True).value_in_unit(unit.nanometer) * 10
        self.dataset.add_candidate(positions, dict(group_id=self.trajectory_id,
                                   step=int(step), time_ps=state.getTime().value_in_unit(
                                       unit.picosecond), source='OpenMM positions'), self.split)

    def record_result(self, molecule, result):
        self.evaluations += 1
        if (self.evaluations - 1) % self.interval:
            return
        if self.dataset._system(molecule) != self.dataset.system:
            raise ValueError('Recorded reference belongs to a different system.')
        frame_id = self.dataset.add_candidate(
            molecule.get_coordinates_in_angstrom(),
            dict(group_id=self.trajectory_id, evaluation=self.evaluations,
                 source='successful physical reference evaluation'), self.split)
        self.dataset.add_result(frame_id, result)
