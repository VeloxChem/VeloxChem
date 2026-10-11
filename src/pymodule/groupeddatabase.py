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


"""
Storage of the grouped MSI v2 model family (ADR-006; SPEC S11): one HDF5
file per model, `im_grouped_<root>.h5`, never a legacy name and never a
legacy file.

Layout (root attributes: schema_version, construction_mode,
construction_policy, model_identity, commit):

    /policy               the policy record as canonical JSON
    /molecule             labels, charge, multiplicity, the Z-matrix as integer arrays
    /level_of_theory      driver class and settings as canonical JSON (restricted to the supported drivers)
    /registry/...         groups, retained terms, node grids, the reference definition, the anchor calc_id
    /raw/<calc_id>        coordinates (bohr), energy, gradient, hessian (if computed); attributes: kind, cut key,
                          level hash, wall time, threads, hessian_kind, rigid-body norm of the gradient, converged
    /derived/cuts/<key>   the node data of every cut (extract_node_data), the relaxation record, the diagnostics
    /derived/model        the model record: cut and response policy, weight scheme, support factor, node orders,
                          terms, node grids, frozen rows, branch warnings
    /validation/<id>      validation records
    /ledger               one row per physical calculation and per reuse, as JSON
    /manifest             every dataset and attribute set of the file with its hash, and the state `complete`

Commit protocol (one writer): lock `<file>.lock` (O_EXCL, host and pid; a
lock whose process no longer exists on this host is reported and removed),
copy the published file to `<file>.staging` (or create it), add the records,
update the manifest and the commit counter, close; reopen read-only and
verify every hash; fsync the staging file, `os.replace` it onto the
published name, fsync the directory; release the lock. A crash before the
replace leaves the published file untouched and a staging file that the
next writer reports and deletes; a crash after it leaves the new file. The
published file is never opened for writing; readers open it read-only.

Identity: `model_identity` is the SHA-256 of the canonical identity record
(groupeddefinitions.model_identity). A calc_id is the SHA-256 of the
level-of-theory hash and the exact coordinate bytes.

`model_family(path)` is the dispatch of ADR-001: 'legacy' for a file without
the schema marker, 'grouped_v2' for a supported schema, a critical error
otherwise. Legacy readers are never given a grouped file: the name never
matches `im_database_*.h5`.
"""

import glob
import hashlib
import json
import os
import shutil
import socket
import time

import h5py
import numpy as np

from .errorhandler import assert_msg_critical
from .groupeddefinitions import (SCHEMA_VERSION, SUPPORTED_POLICY_VALUES, GroupDefinition, GroupedModelPolicyV2,
                                 array_payload_hash, canonical_json, identity_hash, normalize_for_json,
                                 validate_policy_record, z_matrix_as_int_tuples)

LEGACY_PATTERN = 'im_database_*.h5'
GROUPED_PREFIX = 'im_grouped_'
SUPPORTED_LEVEL_DRIVERS = ('XtbDriver', 'ScfRestrictedDriver', 'ScfUnrestrictedDriver')
HESSIAN_KINDS = ('none', 'analytical', 'numerical', 'analytical_numerical_dispersion')
RAW_KINDS = ('pilot', 'anchor', 'singleton', 'probe', 'pair', 'triple', 'response', 'held_out', 'relaxation', 'noise',
             'hessian_check')
TRANSIENT_KEYS = ('mode_rows_all', 'mode_values_all')


# ------------------------------------------------------------------ names and dispatch

def grouped_file_name(root):
    """`im_grouped_<root>.h5`."""
    return f'{GROUPED_PREFIX}{root}.h5'


def is_legacy_name(path):
    base = os.path.basename(path)
    return base.startswith('im_database_') and base.endswith('.h5')


def assert_grouped_name(path):
    base = os.path.basename(path)
    assert_msg_critical(
        not is_legacy_name(path),
        f'grouped MSI storage: {base!r} is a legacy database name (im_database_*.h5); a grouped model is written '
        f'under {GROUPED_PREFIX}<root>.h5 only.')
    assert_msg_critical(
        base.startswith(GROUPED_PREFIX) and base.endswith('.h5'),
        f'grouped MSI storage: {base!r} is not a grouped model name ({GROUPED_PREFIX}<root>.h5).')


def model_family(path):
    """
    'legacy' for an HDF5 file without the root attribute schema_version,
    'grouped_v2' for a supported schema; a critical error for an unknown
    schema or a file that is not HDF5.
    """

    assert_msg_critical(os.path.isfile(path), f'grouped MSI storage: {path} does not exist.')
    assert_msg_critical(h5py.is_hdf5(path), f'grouped MSI storage: {path} is not an HDF5 file.')
    with h5py.File(path, 'r') as h5:
        if 'schema_version' not in h5.attrs:
            return 'legacy'
        version = _attr_str(h5.attrs['schema_version'])
    assert_msg_critical(
        version in SUPPORTED_POLICY_VALUES['schema_version'],
        f'grouped MSI storage: {path} carries the unsupported schema version {version!r}; '
        f'supported: {list(SUPPORTED_POLICY_VALUES["schema_version"])}.')
    return 'grouped_v2'


def assert_no_legacy_database(directory):
    """A construction never runs in a directory that holds a legacy database (KI-19)."""
    existing = sorted(glob.glob(os.path.join(directory, LEGACY_PATTERN)))
    assert_msg_critical(
        not existing,
        f'grouped MSI storage: the directory {directory} contains the legacy database(s) {[os.path.basename(e) for e in existing]}; '
        'run the construction in a new directory.')


def _attr_str(value):
    return value.decode('utf-8') if isinstance(value, bytes) else str(value)


# ------------------------------------------------------------------ hashes and ids

def calc_id(level_hash, coordinates):
    """SHA-256 of the level-of-theory hash and the exact coordinate bytes (float64, (N, 3), bohr)."""
    arr = np.ascontiguousarray(np.asarray(coordinates, dtype=np.float64).reshape(-1, 3))
    digest = hashlib.sha256()
    digest.update(str(level_hash).encode('utf-8'))
    digest.update(repr(arr.shape).encode('utf-8'))
    digest.update(arr.tobytes())
    return digest.hexdigest()


def level_of_theory_record(driver, **settings):
    """
    The restricted level-of-theory record (ADR-006 amendment 1): the driver
    class name must be one of the supported drivers; the settings are the
    values that affect the energy, gradient and Hessian (functional, basis,
    convergence thresholds, grid level, dispersion, Hessian kind).
    """

    assert_msg_critical(driver in SUPPORTED_LEVEL_DRIVERS,
                        f'grouped MSI storage: unsupported level-of-theory driver {driver!r}; supported: {list(SUPPORTED_LEVEL_DRIVERS)}.')
    record = {'driver': driver}
    record.update(normalize_for_json(settings))
    return record


def level_hash(record):
    return identity_hash(record)


def _json_bytes(obj):
    return np.bytes_(canonical_json(obj).encode('utf-8'))


def _dataset_hash(dataset):
    value = dataset[()]
    if isinstance(value, (bytes, np.bytes_)):
        digest = hashlib.sha256()
        digest.update(b'bytes')
        digest.update(bytes(value))
        return digest.hexdigest()
    return array_payload_hash(np.asarray(value))


def _attrs_hash(attrs):
    record = {}
    for key, value in attrs.items():
        if isinstance(value, (bytes, np.bytes_)):
            value = bytes(value).decode('utf-8')
        record[key] = normalize_for_json(value)
    return identity_hash(record)


def file_hashes(h5):
    """Every dataset ('path') and every attribute set ('path@attrs') of the file with its hash, the manifest excluded."""
    records = {'/@attrs': _attrs_hash(h5.attrs)}

    def visit(name, obj):
        path = '/' + name
        if path == '/manifest':
            return
        if isinstance(obj, h5py.Dataset):
            records[path] = _dataset_hash(obj)
        if len(obj.attrs):
            records[path + '@attrs'] = _attrs_hash(obj.attrs)

    h5.visititems(visit)
    return records


# ------------------------------------------------------------------ records

def record_to_group(group, record):
    """
    Writes a record (a dict) into an HDF5 group: numpy arrays with at least
    one dimension as datasets, nested dicts as subgroups, everything else
    (scalars, strings, lists, tuples, None) in the JSON attribute 'meta'.
    Dataclasses are written through their __dict__. Transient keys are
    skipped.
    """

    meta = {}
    for key, value in record.items():
        if key in TRANSIENT_KEYS:
            continue
        name = str(key)
        if hasattr(value, '__dataclass_fields__'):
            value = {k: v for k, v in value.__dict__.items() if k != 'history'}
        if isinstance(value, dict):
            sub = group.create_group(name)
            record_to_group(sub, value)
        elif isinstance(value, np.ndarray) and value.ndim >= 1:
            group.create_dataset(name, data=np.ascontiguousarray(value))
        else:
            meta[name] = normalize_for_json(value)
    group.attrs['meta'] = _json_bytes(meta)


def group_to_record(group):
    """The inverse of record_to_group."""
    record = {}
    if 'meta' in group.attrs:
        record.update(json.loads(bytes(group.attrs['meta']).decode('utf-8')))
    for name, obj in group.items():
        if isinstance(obj, h5py.Group):
            record[name] = group_to_record(obj)
        else:
            record[name] = np.array(obj[()])
    return record


class GroupedInterpolationDatabase:
    """
    One grouped model file. Reads open the published file read-only;
    writes go through `commit`, which applies the copy-on-write protocol.

    :param path: the file (`im_grouped_<root>.h5`)
    :param ostream: for the reports of stale locks and staging files
    """

    def __init__(self, path, ostream=None):
        self.path = os.path.abspath(path)
        assert_grouped_name(self.path)
        self.ostream = ostream
        self.commit_times = []

    # ------------------------------------------------------------ messages

    def _warn(self, message):
        if self.ostream is not None:
            self.ostream.print_warning(message)
            self.ostream.flush()

    def _info(self, message):
        if self.ostream is not None:
            self.ostream.print_info(message)
            self.ostream.flush()

    # ------------------------------------------------------------ state

    @property
    def exists(self):
        return os.path.isfile(self.path)

    @property
    def lock_path(self):
        return self.path + '.lock'

    @property
    def staging_path(self):
        return self.path + '.staging'

    def read(self):
        """The published file, read-only (a context manager)."""
        assert_msg_critical(self.exists, f'grouped MSI storage: {self.path} does not exist.')
        assert_msg_critical(model_family(self.path) == 'grouped_v2',
                            f'grouped MSI storage: {self.path} is not a grouped model file.')
        return h5py.File(self.path, 'r')

    def root_attributes(self):
        with self.read() as h5:
            return {k: (_attr_str(v) if isinstance(v, (bytes, np.bytes_, str)) else normalize_for_json(v))
                    for k, v in h5.attrs.items()}

    def manifest(self):
        with self.read() as h5:
            assert_msg_critical('manifest' in h5, f'grouped MSI storage: {self.path} has no manifest.')
            return json.loads(bytes(h5['manifest'][()]).decode('utf-8'))

    def verify(self):
        """
        Recomputes every hash of the published file and compares it with
        the manifest; a mismatch, a missing or extra record, or a manifest
        that is not complete is a critical error. Returns the manifest.
        """

        with self.read() as h5:
            assert_msg_critical('manifest' in h5, f'grouped MSI storage: {self.path} has no manifest.')
            manifest = json.loads(bytes(h5['manifest'][()]).decode('utf-8'))
            found = file_hashes(h5)
        assert_msg_critical(manifest.get('state') == 'complete',
                            f'grouped MSI storage: the manifest of {self.path} is not complete.')
        expected = manifest['records']
        missing = sorted(set(expected) - set(found))
        extra = sorted(set(found) - set(expected))
        changed = sorted(k for k in set(expected) & set(found) if expected[k] != found[k])
        assert_msg_critical(
            not missing and not extra and not changed,
            f'grouped MSI storage: the manifest of {self.path} does not match the file: missing {missing[:3]}, '
            f'extra {extra[:3]}, changed {changed[:3]} (the file is corrupt or was modified outside the commit protocol).')
        return manifest

    # ------------------------------------------------------------ lock

    def _acquire_lock(self):
        for attempt in range(3):
            try:
                fd = os.open(self.lock_path, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
                os.write(fd, f'{socket.gethostname()} {os.getpid()} {time.time():.3f}\n'.encode('utf-8'))
                os.close(fd)
                return
            except FileExistsError:
                try:
                    with open(self.lock_path) as f:
                        host, pid, _ = f.read().split()
                    pid = int(pid)
                except (OSError, ValueError):
                    host, pid = '?', -1
                stale = False
                if host == socket.gethostname() and pid > 0:
                    try:
                        os.kill(pid, 0)
                    except ProcessLookupError:
                        stale = True
                    except PermissionError:
                        stale = False
                if stale:
                    self._warn(f'grouped MSI storage: removing the stale lock {self.lock_path} of the finished process {pid}.')
                    try:
                        os.remove(self.lock_path)
                    except FileNotFoundError:
                        pass
                    continue
                assert_msg_critical(False, f'grouped MSI storage: {self.path} is locked by {host} pid {pid} ({self.lock_path}).')
        assert_msg_critical(False, f'grouped MSI storage: could not acquire the lock {self.lock_path}.')

    def _release_lock(self):
        try:
            os.remove(self.lock_path)
        except FileNotFoundError:
            pass

    # ------------------------------------------------------------ commit

    @staticmethod
    def _fsync_path(path, directory=False):
        fd = os.open(path, os.O_RDONLY | (os.O_DIRECTORY if directory else 0))
        try:
            os.fsync(fd)
        finally:
            os.close(fd)

    def commit(self, writer, note=''):
        """
        Applies `writer(h5)` to a staging copy of the published file (or to
        a new file) and publishes it atomically with the protocol of
        ADR-006. Returns the new commit number. The wall time of the commit
        is appended to `commit_times`.
        """

        t0 = time.time()
        directory = os.path.dirname(self.path)
        assert_msg_critical(os.path.isdir(directory), f'grouped MSI storage: the directory {directory} does not exist.')
        assert_msg_critical(os.stat(directory).st_dev == os.stat(os.path.dirname(self.staging_path)).st_dev,
                            'grouped MSI storage: the staging file and the published file must lie on one filesystem.')
        self._acquire_lock()
        try:
            if os.path.exists(self.staging_path):
                self._warn(f'grouped MSI storage: removing the staging file {self.staging_path} of an interrupted commit.')
                os.remove(self.staging_path)
            if self.exists:
                shutil.copy2(self.path, self.staging_path)
                previous = self.root_attributes()
                commit_number = int(previous.get('commit', 0)) + 1
            else:
                commit_number = 1
                # the ledger and the manifest are rewritten at every commit; without a persistent free-space manager
                # HDF5 never reclaims the freed space and the file grows by their size per commit (measured: 7.7 MB
                # against 127 kB after 400 commits). The strategy is a property of the file, set at its creation.
                with h5py.File(self.staging_path, 'w', fs_strategy='fsm', fs_persist=True, fs_threshold=1):
                    pass
            with h5py.File(self.staging_path, 'a') as h5:
                writer(h5)
                h5.attrs['commit'] = int(commit_number)
                h5.attrs['schema_version'] = np.bytes_(SCHEMA_VERSION.encode('utf-8'))
                if 'manifest' in h5:
                    del h5['manifest']
                records = file_hashes(h5)
                manifest = {'state': 'complete', 'commit': int(commit_number), 'note': str(note),
                            'written': time.strftime('%Y-%m-%dT%H:%M:%S'), 'records': records}
                h5.create_dataset('manifest', data=_json_bytes(manifest))
                h5.flush()
            # verification of the staging file before it is published
            with h5py.File(self.staging_path, 'r') as h5:
                found = file_hashes(h5)
                stored = json.loads(bytes(h5['manifest'][()]).decode('utf-8'))
            assert_msg_critical(found == stored['records'] and stored['state'] == 'complete',
                                'grouped MSI storage: the staging file failed its verification; nothing was published.')
            self._fsync_path(self.staging_path)
            os.replace(self.staging_path, self.path)
            self._fsync_path(directory, directory=True)
        finally:
            self._release_lock()
        self.commit_times.append(time.time() - t0)
        return commit_number

    # ------------------------------------------------------------ creation

    def create(self, identity_record, model_identity, policy, molecule_record, level_of_theory, registry, note='created'):
        """
        The first commit: root attributes, /policy, /molecule,
        /level_of_theory, /registry and the empty record groups. The file
        must not exist.
        """

        assert_msg_critical(not self.exists, f'grouped MSI storage: {self.path} exists already.')
        assert_no_legacy_database(os.path.dirname(self.path))
        policy_record = validate_policy_record(policy.as_record() if isinstance(policy, GroupedModelPolicyV2) else policy)
        level = level_of_theory_record(level_of_theory['driver'], **{k: v for k, v in level_of_theory.items() if k != 'driver'})

        def writer(h5):
            h5.attrs['schema_version'] = np.bytes_(SCHEMA_VERSION.encode('utf-8'))
            h5.attrs['construction_mode'] = np.bytes_(policy_record['construction_mode'].encode('utf-8'))
            h5.attrs['construction_policy'] = np.bytes_(policy_record['construction_policy'].encode('utf-8'))
            h5.attrs['model_identity'] = np.bytes_(str(model_identity).encode('utf-8'))
            h5.create_dataset('policy', data=_json_bytes(policy_record))
            h5.create_dataset('identity', data=_json_bytes(identity_record))
            h5.create_dataset('level_of_theory', data=_json_bytes(level))
            mol = h5.create_group('molecule')
            mol.create_dataset('labels', data=np.array([str(lb) for lb in molecule_record['labels']], dtype='S4'))
            mol.attrs['charge'] = int(molecule_record['charge'])
            mol.attrs['multiplicity'] = int(molecule_record['multiplicity'])
            zm = z_matrix_as_int_tuples(molecule_record['z_matrix'])
            for section, width in (('bonds', 2), ('angles', 3), ('dihedrals', 4), ('impropers', 4)):
                rows = np.array(zm[section], dtype=np.int64).reshape(-1, width)
                mol.create_dataset(f'z_matrix_{section}', data=rows)
            reg = h5.create_group('registry')
            record_to_group(reg, registry)
            for name in ('raw', 'derived', 'validation'):
                h5.create_group(name)
            h5.create_dataset('ledger', data=_json_bytes([]))

        return self.commit(writer, note=note)

    # ------------------------------------------------------------ reading records

    def policy(self):
        with self.read() as h5:
            return GroupedModelPolicyV2.from_record(json.loads(bytes(h5['policy'][()]).decode('utf-8')))

    def level_of_theory(self):
        with self.read() as h5:
            return json.loads(bytes(h5['level_of_theory'][()]).decode('utf-8'))

    def identity_record(self):
        with self.read() as h5:
            return json.loads(bytes(h5['identity'][()]).decode('utf-8'))

    def molecule_record(self):
        with self.read() as h5:
            mol = h5['molecule']
            return {'labels': [lb.decode('utf-8') for lb in mol['labels'][()]],
                    'charge': int(mol.attrs['charge']), 'multiplicity': int(mol.attrs['multiplicity']),
                    'z_matrix': {section: [tuple(int(a) for a in row) for row in mol[f'z_matrix_{section}'][()]]
                                 for section in ('bonds', 'angles', 'dihedrals', 'impropers')}}

    def registry(self):
        with self.read() as h5:
            return group_to_record(h5['registry'])

    def ledger(self):
        with self.read() as h5:
            return json.loads(bytes(h5['ledger'][()]).decode('utf-8'))

    def raw_ids(self):
        with self.read() as h5:
            return sorted(h5['raw'].keys())

    def has_raw(self, cid):
        """Whether a converged raw record with this calc_id exists (the file's manifest is verified by `verify`, which the construction calls on a restart)."""
        if not self.exists:
            return False
        with self.read() as h5:
            if cid not in h5['raw']:
                return False
            return bool(h5['raw'][cid].attrs.get('converged', True))

    def raw_record(self, cid):
        with self.read() as h5:
            assert_msg_critical(cid in h5['raw'], f'grouped MSI storage: no raw record {cid}.')
            grp = h5['raw'][cid]
            out = {k: (_attr_str(v) if isinstance(v, (bytes, np.bytes_)) else normalize_for_json(v)) for k, v in grp.attrs.items()}
            out['coordinates'] = np.array(grp['coordinates'][()])
            out['energy'] = float(grp['energy'][()])
            out['gradient'] = np.array(grp['gradient'][()]) if 'gradient' in grp else None
            out['hessian'] = np.array(grp['hessian'][()]) if 'hessian' in grp else None
            return out

    def read_record(self, path):
        with self.read() as h5:
            assert_msg_critical(path in h5, f'grouped MSI storage: no record {path}.')
            return group_to_record(h5[path])

    def has_record(self, path):
        if not self.exists:
            return False
        with self.read() as h5:
            return path in h5

    # ------------------------------------------------------------ writers (used inside commit)

    @staticmethod
    def add_raw(h5, cid, kind, coordinates, energy, gradient=None, hessian=None, attrs=None):
        """
        Adds a raw record (inside a commit). attrs: level_hash, wall time,
        threads, hessian_kind, rigid-body norm of the gradient, converged,
        cut key, note.
        """

        assert_msg_critical(kind in RAW_KINDS, f'grouped MSI storage: unknown raw record kind {kind!r}.')
        assert_msg_critical(cid not in h5['raw'], f'grouped MSI storage: the raw record {cid} exists already.')
        grp = h5['raw'].create_group(cid)
        grp.attrs['kind'] = np.bytes_(kind.encode('utf-8'))
        grp.create_dataset('coordinates', data=np.asarray(coordinates, dtype=np.float64).reshape(-1, 3))
        grp.create_dataset('energy', data=np.float64(energy))
        if gradient is not None:
            grp.create_dataset('gradient', data=np.asarray(gradient, dtype=np.float64).reshape(-1, 3))
        if hessian is not None:
            grp.create_dataset('hessian', data=np.asarray(hessian, dtype=np.float64))
        attrs = dict(attrs or {})
        attrs.setdefault('hessian_kind', 'none' if hessian is None else 'numerical')
        attrs.setdefault('converged', True)
        assert_msg_critical(attrs['hessian_kind'] in HESSIAN_KINDS, f'grouped MSI storage: unknown hessian_kind {attrs["hessian_kind"]!r}.')
        for key, value in attrs.items():
            if isinstance(value, str):
                grp.attrs[key] = np.bytes_(value.encode('utf-8'))
            elif isinstance(value, (bool, np.bool_)):
                grp.attrs[key] = bool(value)
            elif isinstance(value, (int, np.integer)):
                grp.attrs[key] = int(value)
            elif isinstance(value, (float, np.floating)):
                grp.attrs[key] = float(value)
            else:
                grp.attrs[key] = np.bytes_(canonical_json(value).encode('utf-8'))

    @staticmethod
    def write_record(h5, path, record, replace=False):
        """Writes a record under `path` (inside a commit); replaces an existing one only when asked."""
        if path in h5:
            assert_msg_critical(replace, f'grouped MSI storage: the record {path} exists already.')
            del h5[path]
        grp = h5.require_group(path)
        record_to_group(grp, record)

    @staticmethod
    def append_ledger(h5, row):
        rows = json.loads(bytes(h5['ledger'][()]).decode('utf-8'))
        rows.append(normalize_for_json(row))
        del h5['ledger']
        h5.create_dataset('ledger', data=_json_bytes(rows))


# ------------------------------------------------------------------ models

def cut_record_for_storage(record):
    """
    The persisted form of a cut record of GroupedChartModel.build: the
    offsets, the rigid start and the node geometry, the energy and residual,
    the node data, the diagnostics, the well-posedness record (soft modes
    included, the full spectrum excluded), the branch diagnostics and the
    relaxation summary (ADR-006 amendment 4).
    """

    out = {}
    for key, value in record.items():
        if key == 'relaxation':
            if value is None:
                out['relaxation'] = {'relaxed': False}
            else:
                out['relaxation'] = {
                    'relaxed': True, 'targets': np.asarray(value.targets, dtype=np.float64),
                    'multipliers': np.asarray(value.multipliers, dtype=np.float64), 'converged': bool(value.converged),
                    'iterations': int(value.iterations), 'evaluations': int(value.evaluations),
                    'projected_gradient_norm': float(value.projected_gradient_norm),
                    'constraint_residual': float(value.constraint_residual), 'step_norm': float(value.step_norm)}
        elif key in ('well_posedness', 'branch', 'data', 'diagnostics', 'phase_ratio', 'x', 'x_rigid', 'offsets',
                     'energy', 'residual', 'calc_ids', 'branch_investigation'):
            out[key] = value
    return out


def store_model(database, model, note='model'):
    """
    Stores a GroupedChartModel: every cut record under /derived/cuts/<key>
    and the model record under /derived/model, in one commit. The anchor
    raw record and the registry are written by the construction.
    """

    cuts = {('_'.join(str(k) for k in key)): cut_record_for_storage(rec) for key, rec in model.node_records.items()}
    model_record = {
        'cut_policy': model.cut_policy, 'response_policy': model.response_policy,
        'weight_scheme': model.weight_scheme, 'support_factor': float(model.support_factor),
        'node_orders': dict(model.node_orders), 'terms': [list(t) for t in model.terms],
        'node_grids': {str(g): np.asarray(grid, dtype=np.float64) for g, grid in model.node_grids.items()},
        'frozen_rows': list(model.frozen_rows),
        'branch_warnings': [w['message'] for w in model.branch_warnings],
        'radii': {('_'.join(str(g) for g in term)): list(bank.radii) for term, bank in model.banks.items()},
    }

    def writer(h5):
        if 'derived/cuts' in h5:
            del h5['derived/cuts']
        grp = h5.require_group('derived/cuts')
        for key, rec in cuts.items():
            record_to_group(grp.create_group(key), rec)
        GroupedInterpolationDatabase.write_record(h5, 'derived/model', model_record, replace=True)

    return database.commit(writer, note=note)


def load_model(database, coordinates, reference):
    """
    Re-creates the GroupedChartModel of a file from its cut records and
    model record with `GroupedChartModel.banks_from_cuts`, on the given
    coordinate model and reference (rebuilt from /molecule, /registry and
    the anchor raw record by the construction module).
    """

    from .groupedmodel import GroupedChartModel

    database.verify()   # ADR-006: loading verifies the manifest
    record = database.read_record('derived/model')
    cut_records = {}
    with database.read() as h5:
        for key, grp in h5['derived/cuts'].items():
            cut_records[tuple(int(k) for k in key.split('_'))] = group_to_record(grp)
    node_grids = {int(g): np.asarray(grid, dtype=np.float64) for g, grid in record['node_grids'].items()}
    node_orders = {k: int(v) for k, v in record['node_orders'].items()}
    terms = [tuple(int(g) for g in t) for t in record['terms']]
    banks = GroupedChartModel.banks_from_cuts(coordinates, cut_records, terms, node_grids, node_orders, coordinates.policy,
                                              float(record['support_factor']), str(record['weight_scheme']))
    model = GroupedChartModel(coordinates, reference, banks, str(record['response_policy']), str(record['cut_policy']),
                              node_records=cut_records, weight_scheme=str(record['weight_scheme']),
                              support_factor=float(record['support_factor']),
                              frozen_rows=[int(r) for r in record.get('frozen_rows', [])])
    model.node_grids = node_grids
    model.node_orders = node_orders
    return model


def group_definitions_from_records(records):
    """GroupDefinition objects from their persisted records."""
    out = []
    for rec in records:
        rec = dict(rec)
        for key in ('moving_atoms', 'cyclic_order', 'rotor_rows', 'offsets'):
            rec[key] = tuple(rec[key])
        rec['domain'] = None if rec.get('domain') is None else tuple(rec['domain'])
        rec['relaxation_zone'] = None if rec.get('relaxation_zone') is None else tuple(rec['relaxation_zone'])
        out.append(GroupDefinition(**rec))
    return out


def molecule_from_record(record, coordinates_bohr):
    """A Molecule with the labels, charge and multiplicity of the record at the given coordinates (bohr)."""
    from .molecule import Molecule

    mol = Molecule(list(record['labels']), np.asarray(coordinates_bohr, dtype=np.float64).reshape(-1, 3), 'bohr')
    mol.set_charge(int(record['charge']))
    mol.set_multiplicity(int(record['multiplicity']))
    return mol


def candidates_from_groups(groups):
    """The RotorCandidate list that re-creates the persisted groups (the cyclic order is the moving order)."""
    from .groupedcoordinates import RotorCandidate

    out = []
    for g in groups:
        moving = tuple(g.cyclic_order) if g.images > 1 else tuple(g.moving_atoms)
        out.append(RotorCandidate(bond=tuple(sorted((int(g.attachment_atom), int(g.axis_atom)))),
                                  attachment_atom=int(g.attachment_atom), axis_atom=int(g.axis_atom),
                                  moving_atoms=moving, group_type=g.group_type, images=int(g.images),
                                  harmonic=int(g.harmonic), mm_rotatable=bool(g.mm_rotatable)))
    return out


def coordinates_from_file(database, x_anchor=None):
    """
    The GroupedCoordinateModel of a file: the molecule and Z-matrix of
    /molecule, the groups of /registry (domains and relaxation zones
    included), the policy of /policy, on the anchor geometry of the anchor
    raw record (or the given one).
    """

    from .groupedcoordinates import GroupedCoordinateModel

    database.verify()   # ADR-006: loading verifies the manifest
    policy = database.policy()
    mol_record = database.molecule_record()
    registry = database.registry()
    groups = group_definitions_from_records(registry['groups'])
    if x_anchor is None:
        anchor = database.raw_record(registry['anchor_calc_id'])
        x_anchor = anchor['coordinates']
    molecule = molecule_from_record(mol_record, x_anchor)
    candidates = candidates_from_groups(groups)
    domains = {c.bond: g.domain for c, g in zip(candidates, groups) if g.domain is not None}
    zones = {c.bond: g.relaxation_zone for c, g in zip(candidates, groups) if g.relaxation_zone is not None}
    model = GroupedCoordinateModel(molecule, mol_record['z_matrix'], x_anchor, candidates, policy=policy,
                                   domains=domains, relaxation_zones=zones)
    stored = [g.as_record() for g in groups]
    rebuilt = [g.as_record() for g in model.groups]
    assert_msg_critical(stored == rebuilt,
                        'grouped MSI storage: the groups re-created from the file differ from the persisted definitions.')
    return model


def reference_from_file(database, coordinates):
    """
    The GroupedReference of a file from the anchor raw record (energy,
    gradient and the symmetric part of the Hessian), as the construction
    built it.
    """

    from .groupedreference import GroupedReference

    registry = database.registry()
    anchor = database.raw_record(registry['anchor_calc_id'])
    assert_msg_critical(anchor['hessian'] is not None, 'grouped MSI storage: the anchor raw record has no Hessian.')
    hessian = np.asarray(anchor['hessian'], dtype=np.float64)
    hessian = 0.5 * (hessian + hessian.T)
    return GroupedReference(coordinates, anchor['energy'], anchor['gradient'], hessian)


def registry_record(coordinates, retained_terms, node_grids, anchor_calc_id, reference=None, extra=None):
    """
    The /registry record: the group definitions, the retained terms, the
    node grids (radians, relative to the anchor), the reference definition
    and the calc_id of the anchor.
    """

    record = {
        'groups': [g.as_record() for g in coordinates.groups],
        'retained_terms': [list(t) for t in retained_terms],
        'node_grids': {str(int(g)): np.asarray(grid, dtype=np.float64) for g, grid in node_grids.items()},
        'anchor_calc_id': str(anchor_calc_id),
        'reference': {'reference_policy': coordinates.policy.reference_policy,
                      'atlas_delta_degrees': float(coordinates.policy.atlas_delta_degrees),
                      'group_harmonics': [int(g.harmonic) for g in coordinates.groups],
                      'group_images': [int(g.images) for g in coordinates.groups]},
        'core_atoms': [int(a) for a in coordinates.core_atoms],
    }
    if extra:
        record.update(extra)
    return record
