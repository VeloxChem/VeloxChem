#
#                                   VELOXCHEM
#              ----------------------------------------------------
#                          An Electronic Structure Code
#
#  SPDX-License-Identifier: BSD-3-Clause
#
#  Copyright 2018-2026 VeloxChem developers

"""Durable, backend-neutral records for geometry-optimization evaluations."""

from copy import deepcopy
from enum import Enum
import json
from pathlib import Path

import numpy as np


_BOHR_IN_ANGSTROM = 0.529177210903


def _json_safe(value):
    """Recursively converts numerical/backend diagnostics to JSON values."""

    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    if isinstance(value, Enum):
        return _json_safe(value.value)
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, dict):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    return str(value)


class OptimizationStepRecorder:
    """Writes every optimizer evaluation and its later accept/reject decision.

    A new ``run_NNNN`` directory is allocated for each optimization, so a
    restart or recalculation cannot overwrite the evidence from an earlier run.
    Each evaluation has a self-contained directory with XYZ, JSON, gradient,
    energy and overlap files.  ``manifest.json`` provides a compact index.
    """

    schema_version = 1

    @classmethod
    def create(cls, root_directory, labels, driver_name):
        """Allocates and returns the next non-destructive run directory."""

        root = Path(root_directory).expanduser().resolve()
        root.mkdir(parents=True, exist_ok=True)
        run_index = 1
        while True:
            run_directory = root / f'run_{run_index:04d}'
            try:
                run_directory.mkdir()
                break
            except FileExistsError:
                run_index += 1
        return cls(run_directory, labels, driver_name)

    def __init__(self, run_directory, labels, driver_name):
        self.run_directory = Path(run_directory)
        self.evaluations_directory = self.run_directory / 'evaluations'
        self.evaluations_directory.mkdir(parents=True, exist_ok=True)
        self.labels = [str(label) for label in labels]
        self.driver_name = str(driver_name)
        self.records = []
        self._active_record = None
        self._write_manifest()

    def record_evaluation(self,
                          evaluation,
                          coordinates_bohr,
                          energy,
                          gradient,
                          electronic,
                          tracking=None,
                          scan_point=None):
        """Writes one successful energy/gradient evaluation."""

        index = int(evaluation)
        coordinates = np.asarray(coordinates_bohr, dtype=float).reshape(-1, 3)
        gradient = np.asarray(gradient, dtype=float).reshape(-1, 3)
        if coordinates.shape != gradient.shape:
            raise ValueError(
                'OptimizationStepRecorder: geometry/gradient shape mismatch.')
        if len(self.labels) != coordinates.shape[0]:
            raise ValueError(
                'OptimizationStepRecorder: atom-label count mismatch.')

        step_directory = self.evaluations_directory / f'step_{index:06d}'
        step_directory.mkdir(parents=False, exist_ok=False)
        coordinates_angstrom = coordinates * _BOHR_IN_ANGSTROM
        energy = float(energy)
        tracking = _json_safe(tracking)
        electronic = _json_safe(electronic)

        record = {
            'schema_version': self.schema_version,
            'evaluation': index,
            'optimizer_status': 'pending',
            'scan_point': None if scan_point is None else int(scan_point),
            'driver': self.driver_name,
            'energy_hartree': energy,
            'coordinates_bohr': coordinates.tolist(),
            'coordinates_angstrom': coordinates_angstrom.tolist(),
            'gradient_hartree_per_bohr': gradient.tolist(),
            'gradient_rms_hartree_per_bohr':
                float(np.sqrt(np.mean(np.sum(gradient**2, axis=1)))),
            'gradient_max_hartree_per_bohr':
                float(np.max(np.linalg.norm(gradient, axis=1))),
            'electronic': electronic,
            'tracking': tracking,
            'files': {
                'geometry': 'geometry.xyz',
                'gradient': 'gradient.csv',
                'record': 'record.json',
            },
        }

        xyz = self._xyz_text(
            coordinates_angstrom,
            f'evaluation={index} energy_hartree={energy:.12f}')
        self._write_text_atomic(step_directory / 'geometry.xyz', xyz)
        self._write_matrix_csv(
            step_directory / 'gradient.csv', gradient,
            header='gx_hartree_per_bohr,gy_hartree_per_bohr,'
                   'gz_hartree_per_bohr')

        self._write_tracking_files(step_directory, record)
        self._write_json_atomic(step_directory / 'record.json', record)
        self._append_text(self.run_directory / 'all_evaluations.xyz', xyz)

        summary = self._summary(record, step_directory)
        self.records.append(summary)
        self._active_record = {
            'record': record,
            'directory': step_directory,
            'summary': summary,
            'xyz': xyz,
        }
        self._write_manifest()

        return step_directory

    def mark_active(self, status, tracking_reference_advanced=None):
        """Marks the latest pending evaluation accepted or rejected."""

        if self._active_record is None:
            return False
        record = self._active_record['record']
        if record['optimizer_status'] != 'pending':
            return False
        if status not in ('accepted', 'rejected'):
            raise ValueError(
                'OptimizationStepRecorder: status must be accepted or rejected.')

        record['optimizer_status'] = status
        if tracking_reference_advanced is not None:
            record['tracking_reference_advanced'] = bool(
                tracking_reference_advanced)
        self._active_record['summary'].update(self._summary(
            record, self._active_record['directory']))
        self._write_json_atomic(
            self._active_record['directory'] / 'record.json', record)
        self._write_manifest()

        trajectory = self.run_directory / f'{status}_trajectory.xyz'
        self._append_text(trajectory, self._active_record['xyz'])
        return True

    def _write_tracking_files(self, step_directory, record):
        tracking = record.get('tracking')
        if not isinstance(tracking, dict):
            return

        state_energies = tracking.get('state_energies')
        if state_energies is not None:
            energies = np.asarray(state_energies, dtype=float).reshape(-1)
            reference = float(np.min(energies))
            rows = np.column_stack((
                np.arange(1, energies.size + 1),
                energies,
                (energies - reference) * 27.211386245988,
            ))
            self._write_matrix_csv(
                step_directory / 'state_energies.csv', rows,
                header='raw_root,total_energy_hartree,relative_to_lowest_ev')
            record['files']['state_energies'] = 'state_energies.csv'
        else:
            # Serenity exposes excitation energies rather than OpenQP's total
            # state-energy array. Persist the same per-root audit file and add
            # total energies when the SCF reference energy is available.
            excitations = tracking.get(
                'state_excitation_energies_hartree',
                tracking.get('state_excitation_energies'))
            if excitations is not None:
                excitations = np.asarray(
                    excitations, dtype=float).reshape(-1)
                reference = (record.get('electronic') or {}).get(
                    'reference_energy')
                totals = (np.full(excitations.size, np.nan)
                          if reference is None else
                          float(reference) + excitations)
                rows = np.column_stack((
                    np.arange(1, excitations.size + 1),
                    excitations,
                    excitations * 27.211386245988,
                    totals,
                ))
                self._write_matrix_csv(
                    step_directory / 'state_energies.csv', rows,
                    header='raw_root,excitation_energy_hartree,'
                           'excitation_energy_ev,total_energy_hartree')
                record['files']['state_energies'] = 'state_energies.csv'

        matrix_fields = {
            'similarity_matrix': 'assignment_similarity.csv',
            'signed_overlap_matrix': 'selected_signed_overlap.csv',
            'native_overlap_matrix': 'native_state_overlap.csv',
            'response_overlap_matrix': 'mrsf_response_overlap.csv',
            'overlap_matrix': 'overlap_matrix.csv',
        }
        for field, filename in matrix_fields.items():
            values = tracking.get(field)
            if values is None:
                continue
            matrix = np.asarray(values, dtype=float)
            if matrix.ndim != 2:
                continue
            self._write_matrix_csv(step_directory / filename, matrix)
            record['files'][field] = filename

    def _summary(self, record, step_directory):
        tracking = record.get('tracking') or {}
        return {
            'evaluation': int(record['evaluation']),
            'optimizer_status': str(record['optimizer_status']),
            'scan_point': record.get('scan_point'),
            'energy_hartree': float(record['energy_hartree']),
            'selected_raw_root': tracking.get(
                'selected_raw_root', tracking.get('new_state')),
            'assignment_reliable': tracking.get(
                'is_reliable', tracking.get('assignment_confident')),
            'overlap_source': tracking.get('overlap_source'),
            'directory': str(step_directory.relative_to(self.run_directory)),
        }

    def _write_manifest(self):
        manifest = {
            'schema_version': self.schema_version,
            'driver': self.driver_name,
            'run_directory': str(self.run_directory),
            'number_of_evaluations': len(self.records),
            'evaluations': deepcopy(self.records),
        }
        self._write_json_atomic(self.run_directory / 'manifest.json', manifest)

    def _xyz_text(self, coordinates_angstrom, comment):
        lines = [str(len(self.labels)), str(comment)]
        for label, xyz in zip(self.labels, coordinates_angstrom):
            lines.append(
                f'{label:<3s}{xyz[0]:20.12f}{xyz[1]:20.12f}'
                f'{xyz[2]:20.12f}')
        return '\n'.join(lines) + '\n'

    @staticmethod
    def _write_json_atomic(filename, data):
        filename = Path(filename)
        temporary = filename.with_suffix(filename.suffix + '.tmp')
        with temporary.open('w', encoding='utf-8') as output:
            json.dump(_json_safe(data), output, indent=2, sort_keys=True)
            output.write('\n')
        temporary.replace(filename)

    @staticmethod
    def _write_text_atomic(filename, text):
        filename = Path(filename)
        temporary = filename.with_suffix(filename.suffix + '.tmp')
        with temporary.open('w', encoding='utf-8') as output:
            output.write(text)
        temporary.replace(filename)

    @staticmethod
    def _write_matrix_csv(filename, matrix, header=None):
        matrix = np.asarray(matrix)
        lines = [] if header is None else [header]
        for row in np.atleast_2d(matrix):
            lines.append(','.join(f'{float(value):.16e}' for value in row))
        OptimizationStepRecorder._write_text_atomic(
            filename, '\n'.join(lines) + '\n')

    @staticmethod
    def _append_text(filename, text):
        with Path(filename).open('a', encoding='utf-8') as output:
            output.write(text)
