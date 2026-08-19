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

from mpi4py import MPI
from copy import deepcopy
from contextlib import redirect_stderr
from io import StringIO
import numpy as np
import time as tm

from .veloxchemlib import mpi_master
from .outputstream import OutputStream
from .scfgradientdriver import ScfGradientDriver
from .molecule import Molecule
from .inputparser import write_unparsed_input_to_hdf5
from .errorhandler import assert_msg_critical
from .optimizationsteprecorder import OptimizationStepRecorder
from .transitiondensitytracker import (StateTrackingError,
                                       StateTrackingStatus)

with redirect_stderr(StringIO()) as fg_err:
    import geometric

# geomeTRIC constraint keywords and the number of atoms each one takes.
_CONSTRAINT_ATOM_COUNT = {
    'distance': 2,
    'angle': 3,
    'dihedral': 4,
}

# A constrained torsion whose central angle approaches 180 degrees is not
# defined, and geomeTRIC aborts with LinearTorsionError at 175 degrees.
_LINEAR_TORSION_WARNING_DEGREES = 165.0


def _internal_coordinate(coordinate, xyz, atoms):
    """
    Evaluates one primitive internal coordinate from Cartesian coordinates.

    :param coordinate:
        ``'distance'``, ``'angle'`` or ``'dihedral'``.
    :param xyz:
        Cartesian coordinates in bohr, shape ``(natoms, 3)``.
    :param atoms:
        Zero-based atom indices.

    :return:
        Distance in Angstrom, angle or dihedral in degrees, or ``None`` when
        the indices are out of range or the coordinate is degenerate.
    """

    if any(not 0 <= index < xyz.shape[0] for index in atoms):
        return None

    if coordinate == 'distance':
        return float(np.linalg.norm(xyz[atoms[1]] - xyz[atoms[0]]) *
                     geometric.nifty.bohr2ang)

    if coordinate == 'angle':
        u = xyz[atoms[0]] - xyz[atoms[1]]
        v = xyz[atoms[2]] - xyz[atoms[1]]
        norms = np.linalg.norm(u) * np.linalg.norm(v)
        if norms < 1.0e-12:
            return None
        return float(np.degrees(
            np.arccos(np.clip(np.dot(u, v) / norms, -1.0, 1.0))))

    if coordinate == 'dihedral':
        b0 = xyz[atoms[0]] - xyz[atoms[1]]
        b1 = xyz[atoms[2]] - xyz[atoms[1]]
        b2 = xyz[atoms[3]] - xyz[atoms[2]]
        b1_norm = np.linalg.norm(b1)
        if b1_norm < 1.0e-12:
            return None
        b1 = b1 / b1_norm
        # Praxeolitic formula: projections of the terminal bonds onto the
        # plane perpendicular to the central bond.
        v = b0 - np.dot(b0, b1) * b1
        w = b2 - np.dot(b2, b1) * b1
        if np.linalg.norm(v) < 1.0e-12 or np.linalg.norm(w) < 1.0e-12:
            return None
        return float(np.degrees(
            np.arctan2(np.dot(np.cross(b1, v), w), np.dot(v, w))))

    return None


class OptimizationEngine(geometric.engine.Engine):
    """
    Implements optimization engine for geomeTRIC.

    :param grad_drv:
        The gradient driver.
    :param molecule:
        The molecule.
    :params **:
        The input parameters of the main energy driver

    Instance variables
        - molecule: The molecule.
        - grad_drv: The gradient driver.
        - comm: The MPI communicator.
        - rank: The rank of MPI process.
    """

    def __init__(self, grad_drv, molecule, *args):
        """
        Initializes optimization engine for geomeTRIC.
        """

        g_molecule = geometric.molecule.Molecule()
        g_molecule.elem = molecule.get_labels()
        g_molecule.xyzs = [
            molecule.get_coordinates_in_bohr() * geometric.nifty.bohr2ang
        ]
        g_molecule.build_topology()

        super().__init__(g_molecule)

        self.molecule = molecule
        self.args = args
        self.grad_drv = grad_drv

        self.comm = grad_drv.comm
        self.rank = grad_drv.comm.Get_rank()

        self.opt_unparsed_input = None
        self.opt_current_step = 0
        self._last_evaluated_coords = None
        self._accepted_tracking_coords = None
        self.intermediate_data_directory = None
        self._step_recorder = None

        # Constraint reporting.  geomeTRIC logs its own "Scan i/N" banner
        # through a logger that VeloxChem captures into a buffer, so nothing
        # about the constrained coordinate reaches the user's output.  The
        # values are recomputed here directly from the Cartesian coordinates
        # of every evaluated geometry.
        self._constraints = []
        self._scan_point_index = 0

    def set_intermediate_data_directory(self, directory):
        """Enables durable per-evaluation optimization records.

        Only the MPI master writes.  Every optimization allocates a new
        ``run_NNNN`` child, so restarting a failed scan point never overwrites
        the structures and overlap decisions from the previous attempt.
        """

        if directory is None or not str(directory).strip():
            self.intermediate_data_directory = None
            self._step_recorder = None
            return None

        if self.rank == mpi_master():
            recorder = OptimizationStepRecorder.create(
                directory,
                self.molecule.get_labels(),
                type(self.grad_drv).__name__)
            run_directory = str(recorder.run_directory)
            self._step_recorder = recorder
        else:
            run_directory = None
        run_directory = self.comm.bcast(run_directory, root=mpi_master())
        self.intermediate_data_directory = run_directory

        if self.rank == mpi_master():
            self.grad_drv.ostream.print_info(
                'Intermediate optimization data: ' + run_directory)
            self.grad_drv.ostream.flush()
        return run_directory

    def lower(self):
        """
        Required in order to get MECI optimization working in geomeTRIC
        """

        return 'custom engine'

    def set_constraints(self, constraints):
        """
        Registers the geomeTRIC constraint specification for reporting.

        :param constraints:
            The ``OptimizationDriver.constraints`` list, e.g.
            ``['scan dihedral 4 1 2 6 0 180 30']``.  Unparseable entries are
            skipped: this is a reporting path and must never be able to break
            an optimization that geomeTRIC itself accepts.
        """

        self._constraints = []
        self._scan_point_index = 0

        for line in (constraints or []):
            fields = str(line).strip().split()
            if len(fields) < 3:
                continue
            kind, coordinate = fields[0].lower(), fields[1].lower()
            n_atoms = _CONSTRAINT_ATOM_COUNT.get(coordinate)
            if n_atoms is None:
                continue
            try:
                atoms = [int(value) for value in fields[2:2 + n_atoms]]
                values = [float(value) for value in fields[2 + n_atoms:]]
            except ValueError:
                continue
            if len(atoms) != n_atoms:
                continue

            targets = None
            if kind == 'scan' and len(values) >= 3:
                targets = np.linspace(values[0], values[1],
                                      int(round(values[2])))
            elif kind == 'set' and values:
                targets = np.array([values[0]])

            self._constraints.append({
                'kind': kind,
                'coordinate': coordinate,
                # geomeTRIC constraint files are one-based.
                'atoms': [index - 1 for index in atoms],
                'targets': targets,
            })

    def _constraint_report(self, coords):
        """
        Builds the per-step constraint lines for the current geometry.

        :param coords:
            Flat Cartesian coordinates in bohr, as handed to ``calc_new``.

        :return:
            A list of formatted lines; empty when nothing is constrained.
        """

        if not self._constraints:
            return []

        xyz = np.asarray(coords, dtype=float).reshape(-1, 3)
        lines = []

        for constraint in self._constraints:
            value = _internal_coordinate(constraint['coordinate'],
                                         xyz, constraint['atoms'])

            unit = ('Angstrom' if constraint['coordinate'] == 'distance'
                    else 'degrees')
            atom_text = '-'.join(str(index + 1)
                                 for index in constraint['atoms'])
            label = f"{constraint['coordinate']} {atom_text}"
            targets = constraint['targets']

            if value is None:
                # A degenerate coordinate is the single most important thing
                # to report, so this path must stay noisy rather than skip the
                # constraint: it is the geometry geomeTRIC is about to reject.
                lines.append(
                    f'  WARNING: {label} is degenerate at this geometry and '
                    'cannot be evaluated')
            elif targets is None or not targets.size:
                lines.append(f'  {label:<26s} = {value:10.3f} {unit} '
                             f'({constraint["kind"]})')
            else:
                # The constraint is enforced, so within a scan point the value
                # sits on its target; nearest-target matching is therefore
                # exact.  The index only ever moves forward, so the transient
                # geometries at the start of a scan point cannot report a step
                # backwards.
                nearest = int(np.argmin(np.abs(targets - value)))
                if constraint['kind'] == 'scan':
                    self._scan_point_index = max(self._scan_point_index,
                                                 nearest)
                    nearest = self._scan_point_index
                target = float(targets[nearest])

                if constraint['kind'] == 'scan':
                    lines.append(
                        f'  Scan point {nearest + 1:d}/{targets.size:d}: '
                        f'{label} target {target:.3f} {unit}')
                else:
                    lines.append(
                        f'  Fixed {label} target {target:.3f} {unit}')
                lines.append(
                    f'  Current value              = {value:10.3f} {unit} '
                    f'(deviation {value - target:+.3e})')

            # A constrained torsion is undefined once either of its two central
            # angles goes linear, and geomeTRIC then aborts with
            # LinearTorsionError.  Checked whatever the dihedral itself did, so
            # the warning survives the degenerate case.  In practice the
            # molecule is escaping the constraint, which on an excited-state
            # scan usually means the gradient is being taken on the wrong
            # surface.
            if constraint['coordinate'] == 'dihedral':
                for triple in (constraint['atoms'][0:3],
                               constraint['atoms'][1:4]):
                    bend = _internal_coordinate('angle', xyz, triple)
                    if (bend is not None and
                            bend >= _LINEAR_TORSION_WARNING_DEGREES):
                        triple_text = '-'.join(str(i + 1) for i in triple)
                        lines.append(
                            f'  WARNING: angle {triple_text} = {bend:.2f} '
                            'degrees is approaching linear; the constrained '
                            'torsion is becoming undefined')

        return lines

    def save_guess_files(self, dirname):
        """
        geomeTRIC hook fired once per *accepted* optimization step.

        Excited-state drivers that follow a state across geometries stage
        their new reference during ``calc_new`` and commit it here, so a
        nuclear trial that geomeTRIC later rejects never becomes the
        comparison point for the next geometry.

        :param dirname:
            The geomeTRIC scratch directory (unused; state is kept in memory).
        """

        commit = getattr(self.grad_drv, 'commit_tracking_reference', None)
        committed = None
        if commit is not None:
            committed = commit()
        if self.rank == mpi_master() and self._step_recorder is not None:
            self._step_recorder.mark_active(
                'accepted', tracking_reference_advanced=committed)
        if self._last_evaluated_coords is not None:
            self._accepted_tracking_coords = self._last_evaluated_coords.copy()

    def load_guess_files(self, dirname):
        """
        geomeTRIC hook fired when a trial step is *rejected*.

        Discards the staged tracking reference so the next trial is still
        compared against the last accepted geometry.

        :param dirname:
            The geomeTRIC scratch directory (unused; state is kept in memory).
        """

        rollback = getattr(self.grad_drv, 'rollback_tracking_reference', None)
        rolled_back = None
        if rollback is not None:
            rolled_back = rollback()
        if self.rank == mpi_master() and self._step_recorder is not None:
            self._step_recorder.mark_active(
                'rejected',
                tracking_reference_advanced=(
                    False if rollback is not None else None))

    def _tracking_record(self):
        """Returns backend-neutral tracking diagnostics for one evaluation."""

        result = getattr(self.grad_drv, 'last_tracking_result', None)
        if result is not None:
            to_dict = getattr(result, 'to_dict', None)
            return to_dict() if to_dict is not None else result
        return getattr(self.grad_drv, 'tracking_info', None)

    def _electronic_record(self):
        """Collects common state-energy/root metadata from the gradient driver."""

        record = {}
        for attribute in (
                'total_energy', 'reference_energy', 'excited_state_energy',
                'selected_excitation_energy', 'target_state_energies',
                'state_deriv_index', 'target_state_index'):
            value = getattr(self.grad_drv, attribute, None)
            if value is not None:
                record[attribute] = value

        tracker = getattr(self.grad_drv, 'state_tracker', None)
        if tracker is not None:
            revision = getattr(tracker, 'reference_revision', None)
            if revision is not None:
                record['tracking_reference_revision'] = int(revision)
        return record

    def _record_intermediate_step(self, coords, energy, gradient):
        """Writes the just-completed optimizer evaluation on the master rank."""

        if self.rank != mpi_master() or self._step_recorder is None:
            return
        scan_point = self._scan_point_index if self._constraints else None
        self._step_recorder.record_evaluation(
            self.opt_current_step,
            coords,
            energy,
            gradient,
            self._electronic_record(),
            tracking=self._tracking_record(),
            scan_point=scan_point)

    def _compute_energy_gradient(self, molecule):
        """Evaluates one molecule using the driver's atomic contract."""

        if (getattr(self.grad_drv, 'computes_energy_with_gradient', False) and
                not getattr(self.grad_drv, 'numerical', False)):
            self.grad_drv.compute(molecule, *self.args)
            energy = getattr(self.grad_drv, 'total_energy', None)
            assert_msg_critical(
                energy is not None,
                'OptimizationEngine.calc_new: the combined energy-gradient '
                'driver did not provide total_energy.')
        else:
            energy = self.grad_drv.compute_energy(molecule, *self.args)
            self.grad_drv.compute(molecule, *self.args)
        return energy, self.grad_drv.get_gradient()

    def _molecule_at_coordinates(self, coords):
        labels = self.molecule.get_labels()
        atom_basis_labels = self.molecule.get_atom_basis_labels()
        if self.rank == mpi_master():
            molecule = Molecule(labels, np.asarray(coords).reshape(-1, 3),
                                'au', atom_basis_labels)
            molecule.set_charge(self.molecule.get_charge())
            molecule.set_multiplicity(self.molecule.get_multiplicity())
        else:
            molecule = Molecule()
        return self.comm.bcast(molecule, root=mpi_master())

    def _retry_low_overlap(self, target_coords, initial_error):
        """Retries LOW_OVERLAP by root expansion and segment subdivision."""

        tracker = getattr(self.grad_drv, 'state_tracker', None)
        begin_retry = getattr(self.grad_drv, 'begin_tracking_retry', None)
        promote = getattr(
            self.grad_drv, 'promote_tracking_reference_for_retry', None)
        rollback = getattr(
            self.grad_drv, 'rollback_tracking_reference', None)
        if (tracker is None or begin_retry is None or promote is None or
                self._accepted_tracking_coords is None):
            raise initial_error

        expansions = int(getattr(tracker, 'max_root_expansions', 0))
        increment = int(getattr(tracker, 'root_window_increment', 0))
        for attempt in range(expansions):
            begin_retry()
            response_driver = getattr(self.grad_drv, 'rsp_driver', None)
            if response_driver is None or increment <= 0:
                break
            response_driver.set_nstates(
                int(response_driver.nstates) + increment)
            self.grad_drv.ostream.print_info(
                'Serenity tracking LOW_OVERLAP: retrying with '
                f'{response_driver.nstates} response roots '
                f'(window retry {attempt + 1}/{expansions}).')
            self.grad_drv.ostream.flush()
            try:
                return self._compute_energy_gradient(
                    self._molecule_at_coordinates(target_coords))
            except StateTrackingError as error:
                if error.status is not StateTrackingStatus.LOW_OVERLAP:
                    if rollback is not None:
                        rollback()
                    raise
                initial_error = error

        max_depth = int(getattr(tracker, 'max_subdivisions', 0))
        if max_depth <= 0:
            if rollback is not None:
                rollback()
            raise initial_error

        begin_retry()
        start = np.asarray(self._accepted_tracking_coords, dtype=float)
        target = np.asarray(target_coords, dtype=float)
        retry_count = {'value': 0}

        def walk_segment(left, right, depth):
            try:
                return self._compute_energy_gradient(
                    self._molecule_at_coordinates(right))
            except StateTrackingError as error:
                if (error.status is not StateTrackingStatus.LOW_OVERLAP or
                        depth <= 0):
                    raise
                midpoint = 0.5 * (left + right)
                retry_count['value'] += 1
                self.grad_drv.ostream.print_info(
                    'Serenity tracking LOW_OVERLAP: bisecting the accepted-to-'
                    f'trial segment (level {max_depth - depth + 1}/'
                    f'{max_depth}).')
                self.grad_drv.ostream.flush()
                walk_segment(left, midpoint, depth - 1)
                assert_msg_critical(
                    promote(),
                    'OptimizationEngine: successful tracking midpoint did '
                    'not stage a reference.')
                return walk_segment(midpoint, right, depth - 1)

        try:
            energy, gradient = walk_segment(start, target, max_depth)
        except Exception:
            if rollback is not None:
                rollback()
            raise

        if getattr(self.grad_drv, 'tracking_info', None) is not None:
            self.grad_drv.tracking_info['adaptive_subdivisions'] = int(
                retry_count['value'])
        return energy, gradient

    def _validate_step_results(self, energy, gradient, natoms):
        """Validates optimizer step results before returning to geomeTRIC.

        :param energy:
            The energy returned by the driver stack.
        :param gradient:
            The gradient returned by the driver stack.
        :param natoms:
            Number of atoms in the current molecule.
        """

        assert_msg_critical(
            np.isscalar(energy),
            'OptimizationEngine.calc_new: expected scalar energy after MPI '
            'synchronization')

        try:
            gradient = np.asarray(gradient, dtype=float)
        except (TypeError, ValueError):
            gradient = None

        assert_msg_critical(
            gradient is not None,
            'OptimizationEngine.calc_new: expected array-like gradient after '
            'MPI synchronization')

        return float(energy), gradient

    def calc_new(self, coords, dirname):
        """
        Implements calc_new method for the engine.

        :param coords:
            The coordinates.
        :param dirname:
            The relative path.

        :return:
            A dictionary containing energy and gradient.
        """

        start_time = tm.time()

        new_mol = self._molecule_at_coordinates(coords)

        title_txt = f'Optimization Step {self.opt_current_step}'
        self.grad_drv.ostream.print_header(title_txt)
        self.grad_drv.ostream.print_header('=' * (len(title_txt) + 2))
        self.grad_drv.ostream.print_blank()

        for line in self._constraint_report(coords):
            self.grad_drv.ostream.print_info(line)
        if self._constraints:
            self.grad_drv.ostream.print_blank()
            # A scan makes many evaluations per grid point, so a tracking
            # history can only be aligned with the converged geometries if the
            # grid point is recorded on each evaluation.
            self.grad_drv.current_scan_point = self._scan_point_index

        self.grad_drv.ostream.print_info('Computing energy and gradient...')
        self.grad_drv.ostream.flush()

        recovery_flag = hasattr(
            self.grad_drv, '_tracking_recovery_active')
        previous_recovery = getattr(
            self.grad_drv, '_tracking_recovery_active', False)
        if recovery_flag:
            self.grad_drv._tracking_recovery_active = True
        try:
            try:
                energy, gradient = self._compute_energy_gradient(new_mol)
            except StateTrackingError as error:
                if error.status is not StateTrackingStatus.LOW_OVERLAP:
                    raise
                try:
                    energy, gradient = self._retry_low_overlap(coords, error)
                except StateTrackingError as exhausted:
                    tracker = getattr(self.grad_drv, 'state_tracker', None)
                    policy = getattr(tracker, 'failure_policy', 'strict')
                    if (exhausted.status is not
                            StateTrackingStatus.LOW_OVERLAP or
                            policy == 'strict'):
                        raise
                    rollback = getattr(
                        self.grad_drv, 'rollback_tracking_reference', None)
                    if rollback is not None:
                        rollback()
                    # Re-evaluate once from the committed reference so the
                    # configured best-effort/adiabatic policy is applied only
                    # after all bounded recovery attempts have failed.
                    if recovery_flag:
                        self.grad_drv._tracking_recovery_active = False
                    energy, gradient = self._compute_energy_gradient(new_mol)
        finally:
            if recovery_flag:
                self.grad_drv._tracking_recovery_active = previous_recovery

        if hasattr(self.grad_drv, "scf_driver") and isinstance(
                self.grad_drv, ScfGradientDriver):
            checkpoint_file = self.grad_drv.scf_driver.get_checkpoint_file()
            if checkpoint_file is not None and self.opt_unparsed_input is not None:
                if self.rank == mpi_master():
                    write_unparsed_input_to_hdf5(checkpoint_file,
                                                 self.opt_unparsed_input,
                                                 group_name='opt_settings')

        energy = self.comm.bcast(energy, root=mpi_master())
        gradient = self.comm.bcast(gradient, root=mpi_master())
        energy, gradient = self._validate_step_results(
            energy, gradient, new_mol.number_of_atoms())

        self._record_intermediate_step(coords, energy, gradient)

        self._last_evaluated_coords = np.asarray(coords, dtype=float).copy()
        if self.opt_current_step == 0:
            # The first electronic evaluation defines the accepted baseline.
            self._accepted_tracking_coords = self._last_evaluated_coords.copy()
        self.opt_current_step += 1

        if self.rank == mpi_master():
            grad2 = np.sum(gradient**2, axis=1)
            rms_grad = np.sqrt(np.mean(grad2))
            max_grad = np.max(np.sqrt(grad2))
            valstr = '  Energy   : {:.10f} a.u.'.format(energy)
            self.grad_drv.ostream.print_info(valstr)
            valstr = '  Gradient : {:.6e} a.u. (RMS)'.format(rms_grad)
            self.grad_drv.ostream.print_info(valstr)
            valstr = '             {:.6e} a.u. (Max)'.format(max_grad)
            self.grad_drv.ostream.print_info(valstr)
            valstr = '  Time     : {:.2f} sec'.format(tm.time() - start_time)
            self.grad_drv.ostream.print_info(valstr)
            self.grad_drv.ostream.print_blank()
            self.grad_drv.ostream.flush()

        return {
            'energy': energy,
            'gradient': gradient.flatten(),
        }

    def copy_scratch(self, src, dest):
        """
        Implements copy_scratch method for the engine.

        :param src:
            The source.
        :param dest:
            The destination.
        """

        return

    def __deepcopy__(self, memo):
        """
        Implements deepcopy.

        :param memo:
            The memo dictionary for deepcopy.

        :return:
            A deepcopy of self.
        """

        new_engine = OptimizationEngine(deepcopy(self.grad_drv),
                                        deepcopy(self.molecule),
                                        *deepcopy(self.args))

        for key, val in vars(self).items():
            if isinstance(val, (MPI.Intracomm, OutputStream)):
                continue
            if key == '_step_recorder':
                # A copied geomeTRIC engine (for example for a Hessian helper)
                # must not race the live optimizer for the same step paths.
                setattr(new_engine, key, None)
                continue
            setattr(new_engine, key, deepcopy(val))

        return new_engine
