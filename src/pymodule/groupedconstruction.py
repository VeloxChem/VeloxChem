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
Construction of a grouped MSI v2 model from quantum-chemical calculations
(ADR-001, ADR-006; PROMPT-M2 section 5): `GroupedIMForceFieldGenerator`.

The generator composes a legacy IMForceFieldGenerator for the QM plumbing
(energy, gradient and Hessian drivers of XtbDriver, ScfRestrictedDriver and
ScfUnrestrictedDriver) and never calls legacy storage. Every physical
calculation is written as a raw record of the grouped database before the
next one starts, with calc_id deduplication (a geometry computed before is
reused from the file, which is the restart), the rigid-body norm of the
gradient, the asymmetry of the Hessian, the wall time and the thread count;
one ledger row per calculation and per reuse.

The construction of one core chart (M2):

1. the anchor: geomeTRIC optimisation at the level used, then a projected
   Newton polish on the QM gradient with the QM Hessian; the anchor's raw
   record holds E, G and H, and the anchor jet of Q_a is verified before
   any bank exists (SPEC S13 step 1);
2. the pre-flight (decision D2): at the anchor and at every node geometry
   (rigid cuts before any QM, relaxed cuts after their relaxation and
   before their Hessian) the condition and the rank gap of the
   complete-coordinate split and the legacy cut-off margin, each with a
   factor of at least 10 to its limit, otherwise the construction stops;
3. `GroupedChartModel.build` with the energy-and-gradient callable for the
   relaxations and the Hessian only at converged nodes;
4. the branch diagnostics (decision D3): a warning starts the investigation
   procedure, and a "reject" verdict stops the construction;
5. the model is stored (`groupeddatabase.store_model`), the ledger and the
   cost summary are reported.

M3B (ADR-009, PROMPT-M3B decisions E1 to E4): `qm_driver` may also be a
grouped-owned QM level of groupedqm.py (`GxtbLevel`, `XtbPythonLevel`), which
bypasses the legacy generator. Every computed Hessian, of any level, is gated
by the validity battery V0 to V5 against independent QM values at its own
geometry; on a failure it is recomputed once (at the level's check SCC
accuracy where the level has one) and a second failure stops the
construction. At relaxed nodes the battery includes the softest constrained
mode, so a negative SOSC curvature counts as a bifurcation only when the
energies confirm it.
"""

import contextlib
import json
import os
import time

import numpy as np

from .errorhandler import assert_msg_critical
from .groupedcoordinates import GroupedCoordinateModel, legacy_topology, propose_groups, rigid_body_projector, select_groups
from .groupeddatabase import (GroupedInterpolationDatabase, assert_no_legacy_database, calc_id, grouped_file_name,
                              level_hash, level_of_theory_record, registry_record, store_model, touch)
from .groupeddefinitions import (GroupedModelPolicyV2, canonical_group_key, canonical_json, hessian_validity_parameters,
                                 model_identity, selection_parameters, validation_parameters)
from .groupedmodel import GroupedChartModel
from .groupedqm import (BATTERY_PROTOCOL, GroupedQMLevel, battery_summary, constrained_soft_mode, hessian_battery, reevaluate_battery,
                        sosc_verdict)
from .groupedreference import GroupedReference
from .groupedselection import select_couplings
from .imforcefieldgenerator import IMForceFieldGenerator
from .molecule import Molecule
from .outputstream import OutputStream
from .scfrestdriver import ScfRestrictedDriver
from .scfunrestdriver import ScfUnrestrictedDriver
from .xtbdriver import XtbDriver

PREFLIGHT_FACTOR = 10.0
# raw-record kinds whose energy+gradient calls are buffered and written in one commit (ADR-006 amendment 8, decision D7):
# the relaxation steps of a node and the relaxation steps of the branch investigation's sweeps (M3)
BUFFERED_KINDS = ('relaxation', 'investigation')
INVESTIGATION_FLUSH_STEPS = 100


class GroupedIMForceFieldGenerator:
    """
    :param qm_driver: XtbDriver, ScfRestrictedDriver or ScfUnrestrictedDriver,
        configured by the caller (functional, convergence, grid level), or a
        grouped QM level of groupedqm.py (GxtbLevel, XtbPythonLevel; ADR-009)
    :param basis_label: the basis set label for SCF drivers
    :param ostream: the output stream (None: silent)
    :param numerical_hessian: force numerical SCF Hessians (default: the
        driver's own setting, analytical for restricted DFT)

    Attributes set by the caller: ``hessian_check`` (default True) runs,
    after every computed Hessian, one central difference of the gradient
    along a random rigid-body-free direction (step ``hessian_check_step``,
    5e-3 bohr) and stops the construction when the relative disagreement
    with H d exceeds ``hessian_check_tolerance`` (0.1); the two gradients
    are raw records of kind ``hessian_check`` and the error is stored on
    the Hessian's record as ``hessian_fd_error`` (M2: the analytic SCF
    Hessian driver uses the SCF driver's last results whatever geometry
    they belong to, so a Hessian computed without an SCF at its own
    geometry is wrong by orders of magnitude; the SCF is now run
    immediately before every analytic Hessian, and the check is the proof).

    M3B (ADR-009, decision E4): with ``hessian_check`` the validity battery
    (`groupedqm.hessian_battery`, tolerances ``hessian_validity``) gates every
    computed Hessian in place of the M2 check's tolerance, whose value is still
    recorded as ``hessian_fd_error``; a failure is recomputed once and a second
    failure stops the construction. ``hessian_test_hook`` (None) is a callable
    (x, hessian, attempt) -> hessian applied to every computed Hessian, for the
    tests. ``precision_contract`` (True) measures, for a grouped level, the
    gradient against central differences of the energy at the anchor and
    stops above ``precision_contract_factor`` (0.1) times
    ``relaxation_gradient_threshold``.
    """

    def __init__(self, qm_driver, basis_label='def2-svp', ostream=None, numerical_hessian=None):
        self.qm_driver = qm_driver
        self.basis_label = str(basis_label)
        self.ostream = ostream
        # M3B (ADR-009): a grouped-owned QM level bypasses the legacy generator, which does not know it
        self.qm_level = qm_driver if isinstance(qm_driver, GroupedQMLevel) else None
        if self.qm_level is not None:
            self.generator = None
            self.drivers = None
        else:
            self.generator = IMForceFieldGenerator(qm_driver, ostream=OutputStream(None))
            self.drivers = self.generator.drivers['gs']
            if numerical_hessian is not None and hasattr(self.drivers[2], 'numerical'):
                self.drivers[2].numerical = bool(numerical_hessian)
        self.level = self._level_record()
        self.level_hash = level_hash(self.level)
        self.threads = int(os.environ.get('OMP_NUM_THREADS', '0') or 0)
        self.database = None
        self.labels = None
        self.charge = 0
        self.multiplicity = 1
        self.counts = {'energy_gradient': 0, 'hessian': 0, 'reused': 0, 'battery_single_points': 0}
        self.wall = {'energy_gradient': 0.0, 'hessian': 0.0}
        self.preflight_records = []
        self._basis = None
        self._kind = 'relaxation'
        # calculations made before the file exists (the anchor optimisation) are kept here and committed as raw
        # records as soon as the file is created
        self._pending = {}
        # every calculation of this generator, by calc_id, so that a second model file of the same campaign (another
        # grid or policy) reuses a geometry computed for the first one: it is then written into the new file as a raw
        # record and counted as reused
        self._memory = {}
        self.hessian_check = True
        self.hessian_check_step = 5.0e-3
        self.hessian_check_tolerance = 0.1
        self.hessian_checks = []
        # the pre-flight margin of PROMPT-M2 section 5.2 (D2): the smallest factor to a limit at the anchor and the node
        # geometries below which the construction stops. The default is the M2 value; a caller that lowers it for a
        # molecule whose split condition is a known, stable property (melatonin's fused rings: a factor 6, ADR-007
        # scope statement) records the value used in the validation record (M3)
        self.preflight_factor = PREFLIGHT_FACTOR
        # M3 (ADR-006 amendment 8, decision D7): the relaxation gradients of one node are buffered here and written in
        # one commit when the node has converged or failed (before its Hessian, and at the end of the construction);
        # the reuse rows of relaxation steps are buffered with them
        self._relaxation_buffer = {}
        self._pending_rows = []
        # M3: the validation parameters (D4) persisted with the validation record, the selection (ADR-008)
        self.validation = validation_parameters()
        self.selection_parameters = None
        self.selection_record = None
        self.node_kind = None
        # M3B (ADR-009, decision E4): the battery tolerances (validation parameters), the per-Hessian log, the test hook,
        # the precision contract of a grouped level at the anchor
        self.hessian_validity = hessian_validity_parameters()
        self.hessian_validity_log = []
        self.hessian_test_hook = None
        self.precision_contract = True
        self.precision_contract_factor = 0.1
        self.precision_contract_record = None
        self._last_hessian_protocol = None
        self._cut_policy = None
        self._frozen_rows = []

    # ------------------------------------------------------------ level of theory

    def _level_record(self):
        drv = self.qm_driver
        if self.qm_level is not None:
            return self.qm_level.level_record()
        if isinstance(drv, XtbDriver):
            return level_of_theory_record('XtbDriver', method=str(drv.xtb_method), hessian='numerical')
        name = 'ScfRestrictedDriver' if isinstance(drv, ScfRestrictedDriver) else (
            'ScfUnrestrictedDriver' if isinstance(drv, ScfUnrestrictedDriver) else type(drv).__name__)
        hess = self.drivers[2]
        numerical = bool(getattr(hess, 'numerical', False))
        dispersion = bool(getattr(drv, 'dispersion', False))
        kind = 'numerical' if numerical else ('analytical_numerical_dispersion' if dispersion else 'analytical')
        return level_of_theory_record(name, xcfun=str(drv.xcfun) if drv.xcfun is not None else 'hf', basis=self.basis_label,
                                      conv_thresh=float(drv.conv_thresh), max_iter=int(drv.max_iter),
                                      grid_level=(None if drv.grid_level is None else int(drv.grid_level)),
                                      dispersion=dispersion, hessian=kind)

    @property
    def hessian_kind(self):
        return self.level['hessian']

    def _print(self, message, warning=False):
        if self.ostream is not None:
            (self.ostream.print_warning if warning else self.ostream.print_info)(message)
            self.ostream.flush()

    # ------------------------------------------------------------ QM calls

    def _molecule(self, x):
        mol = Molecule(list(self.labels), np.asarray(x, dtype=np.float64).reshape(-1, 3), 'bohr')
        mol.set_charge(int(self.charge))
        mol.set_multiplicity(int(self.multiplicity))
        return mol

    def _basis_for(self, mol):
        if self.qm_level is not None or isinstance(self.qm_driver, XtbDriver):
            return None
        from .molecularbasis import MolecularBasis
        return MolecularBasis.read(mol, self.basis_label, ostream=OutputStream(None))

    def _compute_energy_gradient(self, x):
        if self.qm_level is not None:
            energy, gradient = self.qm_level.energy_gradient(self.labels, x, self.charge, self.multiplicity)
            return float(energy), np.asarray(gradient, dtype=np.float64).reshape(-1, 3), True
        mol = self._molecule(x)
        basis = self._basis_for(mol)
        energy, scf_results, _ = self.generator._compute_energy(self.drivers[0], mol, basis)
        gradient = self.generator._compute_gradient(self.drivers[1], mol, basis, scf_results)
        converged = bool(getattr(self.qm_driver, 'is_converged', True)) if not isinstance(self.qm_driver, XtbDriver) else True
        return float(np.asarray(energy).reshape(-1)[0]), np.asarray(gradient, dtype=np.float64).reshape(-1, 3), converged

    def _compute_hessian(self, x):
        if self.qm_level is not None:
            # E3: central differences of fresh single points by parallel workers; the protocol record (with the pair
            # check of V5) is kept for the battery
            hessian, self._last_hessian_protocol = self.qm_level.hessian(self.labels, x, self.charge, self.multiplicity)
            return hessian
        mol = self._molecule(x)
        basis = self._basis_for(mol)
        hess = self.drivers[2]
        scf_driver = getattr(hess, 'scf_driver', None)
        if scf_driver is not None and not isinstance(self.qm_driver, XtbDriver):
            # ScfHessianDriver.compute takes the SCF driver's last results, whatever geometry they belong to, and runs
            # the SCF itself only when there are none (found in M2: the analytic PBE Hessians of rigid nodes whose
            # energy and gradient had been reused from a record were computed from the density of another geometry
            # and had eigenvalues down to -65 hartree/bohr^2). The SCF is therefore run at x immediately before.
            self.generator._compute_energy(scf_driver, mol, basis)
        hessian = self.generator._compute_hessian(hess, mol, basis)
        n = 3 * len(self.labels)
        return np.asarray(hessian, dtype=np.float64).reshape(n, n)

    def _ledger_row(self, cid, kind, computed, wall, reused, extra=None):
        row = {'date': time.strftime('%Y-%m-%d %H:%M:%S'), 'calc_id': cid, 'kind': kind, 'level_hash': self.level_hash,
               'computed': computed, 'hessian_kind': self.hessian_kind if 'H' in computed else 'none',
               'wall_time': float(wall), 'threads': self.threads, 'reused': bool(reused)}
        if extra:
            row.update(extra)
        return row

    def _reuse_row(self, cid, kind, computed, extra=None):
        """A ledger row of a reuse: committed at once, or buffered with the node's relaxation steps (D7)."""
        row = self._ledger_row(cid, kind, computed, 0.0, True, extra)
        if self.database is None:
            return
        if kind in BUFFERED_KINDS:
            self._pending_rows.append(row)
        else:
            self.database.commit(lambda h5: GroupedInterpolationDatabase.append_ledger(h5, row), note='reuse')

    def flush_relaxation_buffer(self, note='relaxation steps of one node'):
        """
        Writes the buffered relaxation gradients and their ledger rows in one
        commit (ADR-006 amendment 8). Called before every Hessian, by the
        node source, and at the end of a construction.
        """

        if self.database is None or (not self._relaxation_buffer and not self._pending_rows):
            return
        entries = list(self._relaxation_buffer.values())
        rows = list(self._pending_rows)

        def writer(h5):
            for entry in entries:
                if entry['cid'] not in h5['raw']:
                    GroupedInterpolationDatabase.add_raw(h5, entry['cid'], entry['kind'], entry['x'], entry['energy'],
                                                         entry['gradient'], entry['hessian'], attrs=entry['attrs'])
            GroupedInterpolationDatabase.append_ledger_rows(h5, [r for entry in entries for r in entry['rows']] + rows)

        self.database.commit(writer, note=note)
        self._relaxation_buffer = {}
        self._pending_rows = []

    def energy_and_gradient(self, x, kind=None, note=''):
        """
        E and the Cartesian gradient at x, from the file when the geometry
        was computed before (a ledger row marks the reuse), otherwise
        computed and committed as a raw record before returning. Relaxation
        steps (kind 'relaxation') are buffered and committed per node (D7).
        """

        kind = self._kind if kind is None else kind
        x = np.asarray(x, dtype=np.float64).reshape(-1, 3)
        cid = calc_id(self.level_hash, x)
        db = self.database
        buffered = kind in BUFFERED_KINDS and db is not None
        if cid in self._relaxation_buffer:
            entry = self._relaxation_buffer[cid]
            self.counts['reused'] += 1
            entry['rows'].append(self._ledger_row(cid, kind, 'EG', 0.0, True))
            return entry['energy'], entry['gradient']
        if db is not None and db.has_raw(cid):
            rec = db.raw_record(cid)
            if rec['gradient'] is not None:
                self.counts['reused'] += 1
                self._reuse_row(cid, kind, 'EG')
                return rec['energy'], rec['gradient']
        if cid in self._memory:
            # computed for another file of this campaign: written into this file as a raw record, counted as reused
            mem = self._memory[cid]
            self.counts['reused'] += 1
            if db is not None:
                row = self._ledger_row(cid, kind, 'EG', 0.0, True, {'note': 'reused from an earlier file of this session'})
                if buffered:
                    self._relaxation_buffer[cid] = {'cid': cid, 'kind': kind, 'x': x, 'energy': mem['energy'], 'gradient': mem['gradient'],
                                                    'hessian': mem['hessian'], 'attrs': dict(mem['attrs']), 'rows': [row]}
                else:
                    def writer(h5):
                        if cid not in h5['raw']:
                            GroupedInterpolationDatabase.add_raw(h5, cid, kind, x, mem['energy'], mem['gradient'], mem['hessian'],
                                                                 attrs=mem['attrs'])
                        GroupedInterpolationDatabase.append_ledger(h5, row)
                    db.commit(writer, note='reuse from memory')
            return mem['energy'], mem['gradient']
        t0 = time.time()
        energy, gradient, converged = self._compute_energy_gradient(x)
        wall = time.time() - t0
        self.counts['energy_gradient'] += 1
        self.wall['energy_gradient'] += wall
        rigid = float(np.linalg.norm(gradient.reshape(-1) - rigid_body_projector(x) @ gradient.reshape(-1)))
        attrs = {'level_hash': self.level_hash, 'wall_time_energy_gradient': wall, 'threads': self.threads,
                 'hessian_kind': 'none', 'rigid_body_gradient_norm': rigid, 'converged': converged, 'note': str(note)}
        self._memory[cid] = {'energy': energy, 'gradient': gradient, 'hessian': None, 'attrs': dict(attrs)}
        row = self._ledger_row(cid, kind, 'EG', wall, False, {'rigid_body_gradient_norm': rigid})
        if buffered:
            self._relaxation_buffer[cid] = {'cid': cid, 'kind': kind, 'x': x, 'energy': energy, 'gradient': gradient,
                                            'hessian': None, 'attrs': attrs, 'rows': [row]}
        elif db is not None:
            def writer(h5):
                GroupedInterpolationDatabase.add_raw(h5, cid, kind, x, energy, gradient, None, attrs=attrs)
                GroupedInterpolationDatabase.append_ledger(h5, row)
            db.commit(writer, note=f'{kind} EG')
        else:
            self._pending[cid] = {'kind': kind, 'x': x, 'energy': energy, 'gradient': gradient, 'hessian': None, 'attrs': attrs,
                                  'rows': [self._ledger_row(cid, kind, 'EG', wall, False, {'rigid_body_gradient_norm': rigid})]}
        assert_msg_critical(converged, f'grouped MSI construction: the QM calculation at calc_id {cid[:12]} did not converge.')
        return energy, gradient

    def hessian(self, x, kind=None, note='', constrained=False):
        """
        The raw (unsymmetrised) QM Hessian at x, from the file when present,
        otherwise computed, validated (the battery of ADR-009, decision E4)
        and added to the geometry's raw record (which the energy-and-gradient
        call created) before returning. The asymmetry max |H - H^T|, the M2
        check value and the battery result are recorded on the raw record and
        in the ledger row. `constrained`: x is a phase-constrained stationary
        point (a relaxed node, a probe or the anchor of the relaxed policy);
        the battery then includes the softest constrained mode, and a negative
        SOSC curvature counts as a bifurcation only when the energies confirm
        it.
        """

        kind = self._kind if kind is None else kind
        x = np.asarray(x, dtype=np.float64).reshape(-1, 3)
        cid = calc_id(self.level_hash, x)
        db = self.database
        # D7: a node's relaxation steps are written before its Hessian (the node geometry's energy and gradient are
        # among them)
        self.flush_relaxation_buffer()
        if db is not None and db.has_raw(cid):
            rec = db.raw_record(cid)
            if rec['hessian'] is not None:
                self.counts['reused'] += 1
                if self.hessian_check and 'hessian_valid' not in rec:
                    # E4: a stored Hessian without a battery result (a file written before M3B) is checked before use
                    self._validate_stored_hessian(cid, x, rec, kind, constrained)
                elif self.hessian_check and 'hessian_validity' in rec:
                    # a stored result is judged again under the current tolerances (no QM); one that fails them (another
                    # tolerance at the time, review of 2026-10-13) is checked again with QM before use
                    stored = reevaluate_battery(json.loads(rec['hessian_validity']), self.hessian_validity)
                    if not stored['valid']:
                        self._validate_stored_hessian(cid, x, rec, kind, constrained)
                db.commit(lambda h5: GroupedInterpolationDatabase.append_ledger(
                    h5, self._ledger_row(cid, kind, 'H', 0.0, True)), note='reuse')
                return rec['hessian']
        if db is None and cid in self._pending and self._pending[cid]['hessian'] is not None:
            self.counts['reused'] += 1
            return self._pending[cid]['hessian']
        if cid in self._memory and self._memory[cid]['hessian'] is not None:
            mem = self._memory[cid]
            self.counts['reused'] += 1
            if db is not None:
                def writer(h5):
                    raw = h5['raw']
                    if cid in raw:
                        if 'hessian' not in raw[cid]:
                            GroupedInterpolationDatabase.add_hessian(
                                h5, cid, mem['hessian'], {key: value for key, value in mem['attrs'].items()
                                                          if key.startswith('hessian_') or key == 'wall_time_hessian'})
                    else:
                        GroupedInterpolationDatabase.add_raw(h5, cid, kind, x, mem['energy'], mem['gradient'], mem['hessian'],
                                                             attrs=mem['attrs'])
                    GroupedInterpolationDatabase.append_ledger(h5, self._ledger_row(cid, kind, 'H', 0.0, True,
                                                                                    {'note': 'reused from an earlier file of this session'}))
                db.commit(writer, note='reuse from memory')
            return mem['hessian']
        t0 = time.time()
        hessian, protocol = self._computed_hessian(x, attempt=0)
        wall = time.time() - t0
        fd_error = self._check_hessian(cid, x, hessian, kind) if self.hessian_check else None
        validity = None
        recomputed = False
        if self.hessian_check:
            # E4: the battery decides (the M2 value above is recorded only); one recompute, then a stop
            record_eg = self._record_energy_gradient(cid, x)
            validity = self._battery(cid, x, hessian, protocol, constrained, record_eg)
            if not validity['valid']:
                first = validity
                accuracy = self.qm_level.check_accuracy if self.qm_level is not None else None
                self._print(f'grouped MSI construction: the Hessian at calc_id {cid[:12]} ({kind}) fails the validity battery '
                            f'({", ".join(first["failed"])}); it is recomputed once'
                            + (f' with every single point at SCC accuracy {accuracy:g}' if accuracy is not None else ''), warning=True)
                t1 = time.time()
                hessian, protocol = self._computed_hessian(x, attempt=1)
                wall += time.time() - t1
                fd_error = self._check_hessian(cid, x, hessian, kind)
                validity = self._battery(cid, x, hessian, protocol, constrained, record_eg)
                validity['first_attempt'] = battery_summary(first)
                recomputed = True
                if not validity['valid']:
                    self._log_validity(cid, kind, validity, recomputed, stopped=True)
                    self._commit_failed_battery(cid, kind, first, validity)
                    worst = battery_summary(validity)
                    assert_msg_critical(
                        False,
                        f'grouped MSI construction: the {self.hessian_kind} Hessian at calc_id {cid[:12]} ({kind}) fails the validity '
                        f'battery twice (first {first["failed"]}, after the recompute {validity["failed"]}; V1 {worst["V1"]:.1e}, '
                        f'V2 {worst["V2"]:.1e}, V3 {worst["V3"]:.1e}, V4 {worst["V4"]}): it does not belong to the surface at this '
                        f'geometry (failed single points, a stale SCF state or a failed response calculation). The Hessian is not '
                        f'stored; the construction stops (PROMPT-M3B decision E4).')
            self._sosc_rule(cid, kind, validity)
            self._log_validity(cid, kind, validity, recomputed)
        self.wall['hessian'] += wall
        n_sp = 0 if protocol is None else int(protocol.get('single_points', 0))
        asymmetry = float(np.abs(hessian - hessian.T).max())
        h_attrs = {'hessian_kind': self.hessian_kind, 'hessian_asymmetry': asymmetry, 'wall_time_hessian': float(wall)}
        row_extra = {'hessian_asymmetry': asymmetry}
        if fd_error is not None:
            # without hessian_check no M2 value exists; the key is omitted (a NaN cannot enter the canonical ledger)
            h_attrs['hessian_fd_error'] = fd_error
            row_extra['hessian_fd_error'] = fd_error
        if protocol is not None:
            h_attrs['hessian_protocol'] = canonical_json({k: protocol[k] for k in ('scheme', 'step_bohr', 'scc_accuracy', 'single_points')})
            row_extra['hessian_single_points'] = n_sp
        if validity is not None:
            h_attrs.update({'hessian_valid': bool(validity['valid']), 'hessian_recomputed': bool(recomputed),
                            'hessian_validity': canonical_json(validity)})
            row_extra.update({'hessian_valid': bool(validity['valid']), 'hessian_recomputed': bool(recomputed),
                              'hessian_validity': battery_summary(validity)})
        if cid not in self._memory:
            energy_m, gradient_m = self._record_energy_gradient(cid, x)
        self._memory[cid]['hessian'] = hessian
        self._memory[cid]['attrs'].update(h_attrs)
        if db is None:
            if cid not in self._pending:
                energy, gradient = self.energy_and_gradient(x, kind=kind, note=note)
            if cid not in self._pending:
                # the energy and gradient came from the session memory (another file of the campaign): a pending entry is
                # made from it (review of 2026-10-13: a KeyError otherwise)
                mem = self._memory[cid]
                self._pending[cid] = {'kind': kind, 'x': x, 'energy': mem['energy'], 'gradient': mem['gradient'], 'hessian': None,
                                      'attrs': {k: v for k, v in mem['attrs'].items() if not k.startswith('hessian_') and k != 'wall_time_hessian'},
                                      'rows': [self._ledger_row(cid, kind, 'EG', 0.0, True, {'note': 'reused from an earlier file of this session'})]}
            entry = self._pending[cid]
            entry['hessian'] = hessian
            entry['kind'] = kind
            entry['attrs'].update(h_attrs)
            entry['rows'].append(self._ledger_row(cid, kind, 'H', wall, False, row_extra))
            return hessian

        def writer(h5):
            raw = h5['raw']
            if cid in raw:
                GroupedInterpolationDatabase.add_hessian(h5, cid, hessian, h_attrs, kind=kind)
            else:
                mem = self._memory[cid]
                attrs = {'level_hash': self.level_hash, 'threads': self.threads, 'converged': bool(mem['attrs'].get('converged', True)),
                         'rigid_body_gradient_norm': float(np.linalg.norm(
                             mem['gradient'].reshape(-1) - rigid_body_projector(x) @ mem['gradient'].reshape(-1)))}
                attrs.update(h_attrs)
                GroupedInterpolationDatabase.add_raw(h5, cid, kind, x, mem['energy'], mem['gradient'], hessian, attrs=attrs)
            GroupedInterpolationDatabase.append_ledger(h5, self._ledger_row(cid, kind, 'H', wall, False, row_extra))
        db.commit(writer, note=f'{kind} H')
        return hessian

    # ------------------------------------------------------------ M3B: data validity (ADR-009, decision E4)

    def _computed_hessian(self, x, attempt=0):
        """
        A newly computed Hessian (3N, 3N) and its protocol record (None for
        the legacy drivers). attempt 1 is the recompute of E4: for a grouped
        level every single point at its check SCC accuracy; a legacy driver
        has no such control and repeats its own protocol.
        """

        self._last_hessian_protocol = None
        if attempt == 0 or self.qm_level is None:
            hessian = self._compute_hessian(x)
        else:
            hessian, self._last_hessian_protocol = self.qm_level.hessian(self.labels, x, self.charge, self.multiplicity,
                                                                         accuracy=self.qm_level.check_accuracy)
        n = 3 * len(self.labels)
        hessian = np.asarray(hessian, dtype=np.float64).reshape(n, n)
        if self.hessian_test_hook is not None:
            hessian = np.asarray(self.hessian_test_hook(x, hessian, attempt), dtype=np.float64).reshape(n, n)
        self.counts['hessian'] += 1
        return hessian, self._last_hessian_protocol

    def _record_energy_gradient(self, cid, x):
        """The record's energy and gradient at x (memory, buffer, pending, file, else computed into the memory), without a ledger row."""
        if cid in self._memory:
            return self._memory[cid]['energy'], self._memory[cid]['gradient']
        for store in (self._relaxation_buffer, self._pending):
            if cid in store:
                entry = store[cid]
                self._memory[cid] = {'energy': entry['energy'], 'gradient': entry['gradient'], 'hessian': None, 'attrs': dict(entry['attrs'])}
                return entry['energy'], entry['gradient']
        if self.database is not None and self.database.has_raw(cid):
            rec = self.database.raw_record(cid)
            if rec['gradient'] is not None:
                self._memory[cid] = {'energy': rec['energy'], 'gradient': rec['gradient'], 'hessian': None,
                                     'attrs': {'level_hash': self.level_hash, 'threads': self.threads,
                                               'converged': bool(rec.get('converged', True)),
                                               'rigid_body_gradient_norm': float(rec.get('rigid_body_gradient_norm', 0.0))}}
                return rec['energy'], rec['gradient']
        energy, gradient, converged = self._compute_energy_gradient(x)
        self.counts['energy_gradient'] += 1
        self._memory[cid] = {'energy': energy, 'gradient': gradient, 'hessian': None,
                             'attrs': {'level_hash': self.level_hash, 'threads': self.threads, 'converged': converged,
                                       'rigid_body_gradient_norm': float(np.linalg.norm(
                                           gradient.reshape(-1) - rigid_body_projector(x) @ gradient.reshape(-1)))}}
        return energy, gradient

    def _check_source(self):
        """The single-point source of the battery: the level at its check SCC accuracy (in parallel), or the legacy driver."""
        if self.qm_level is not None:
            level = self.qm_level
            return lambda xs: level.energy_gradient_many(self.labels, xs, self.charge, self.multiplicity, accuracy=level.check_accuracy)
        return lambda xs: [self._compute_energy_gradient(y)[:2] for y in xs]

    def _battery(self, cid, x, hessian, protocol, constrained, record_eg):
        """The validity battery of one Hessian (groupedqm.hessian_battery) with the construction's tolerances."""
        mode = None
        if constrained and getattr(self, 'coordinates', None) is not None:
            mode = constrained_soft_mode(self.coordinates, x, record_eg[1], hessian, frozen_rows=self._frozen_rows)
        result = hessian_battery(hessian, x, self._check_source(), record=record_eg, constrained=mode,
                                 pair_check=None if protocol is None else protocol.get('pair_check'),
                                 parameters=self.hessian_validity, seed=int(cid[8:16], 16))
        self.counts['battery_single_points'] = self.counts.get('battery_single_points', 0) + int(result['single_points'])
        result['check_accuracy'] = self.qm_level.check_accuracy if self.qm_level is not None else 'the record level (legacy driver)'
        if protocol is not None:
            result['protocol'] = {k: protocol[k] for k in ('step_bohr', 'scc_accuracy', 'single_points', 'wall_time')}
        return result

    def _sosc_rule(self, cid, kind, validity):
        """
        E4: at a constrained stationary point the SOSC verdict is decided on
        the validated Hessian. The battery passed, so along the softest
        constrained mode the curvature from energies agrees with the
        Hessian's in sign: a negative one is confirmed (a genuine
        bifurcation; the build's well-posedness check then stops as before).
        """

        c = validity.get('V3', {}).get('constrained')
        if c is None:
            return
        verdict = sosc_verdict(validity, float(self.coordinates.policy.relaxation_curvature_min))
        validity['sosc'] = {'lambda_min': c['lambda_min'], 'lambda_energy': c['lambda_energy'], 'verdict': verdict}
        if verdict == 'bifurcation':
            self._print(f'grouped MSI construction: at calc_id {cid[:12]} ({kind}) the smallest constrained curvature '
                        f'{c["lambda_min"]:+.3e} is negative and confirmed by energies ({c["lambda_energy"]:+.3e}) on a valid Hessian: a '
                        f'genuine bifurcation of the valley (PROMPT-M3B section 9: stop and report)', warning=True)

    def _log_validity(self, cid, kind, validity, recomputed, stopped=False):
        entry = battery_summary(validity)
        entry.update({'calc_id': cid[:12], 'kind': kind, 'recomputed': bool(recomputed), 'stopped': bool(stopped),
                      'checked_on_reuse': bool(validity.get('checked_on_reuse', False))})
        if 'sosc' in validity:
            entry['sosc'] = validity['sosc']['verdict']
        self.hessian_validity_log.append(entry)

    def _commit_failed_battery(self, cid, kind, first, second):
        """The two failed battery results go into the ledger before the stop (the file records why it stopped)."""
        if self.database is None:
            return
        row = self._ledger_row(cid, kind, 'battery', 0.0, False, {'hessian_valid': False, 'stopped': True,
                                                                  'first_attempt': battery_summary(first), 'recompute': battery_summary(second)})
        self.database.commit(lambda h5: GroupedInterpolationDatabase.append_ledger(h5, row), note='failed Hessian battery')

    def _validate_stored_hessian(self, cid, x, rec, kind, constrained):
        """E4: a reused Hessian without a battery result is checked before use; the result is written onto its record."""
        validity = self._battery(cid, x, rec['hessian'], None, constrained, (rec['energy'], rec['gradient']))
        validity['checked_on_reuse'] = True
        validity['V5'] = 'not available (a stored Hessian: its displaced single points are not kept)'
        if validity['valid']:
            self._sosc_rule(cid, kind, validity)
        self._log_validity(cid, kind, validity, False, stopped=not validity['valid'])
        if not validity['valid']:
            row = self._ledger_row(cid, kind, 'battery', 0.0, True, {'hessian_valid': False, 'stopped': True, 'checked_on_reuse': True,
                                                                     'battery': battery_summary(validity)})
            self.database.commit(lambda h5: GroupedInterpolationDatabase.append_ledger(h5, row), note='failed battery of a stored Hessian')
        assert_msg_critical(
            validity['valid'],
            f'grouped MSI construction: the stored Hessian at calc_id {cid[:12]} ({kind}) fails the validity battery '
            f'({validity["failed"]}): it does not belong to the surface at this geometry. A stored record is not replaced; build the '
            f'model in a new file (PROMPT-M3B decision E4).')
        text = canonical_json(validity)

        def writer(h5):
            touch(h5, f'/raw/{cid}')
            grp = h5['raw'][cid]
            grp.attrs['hessian_valid'] = True
            grp.attrs['hessian_recomputed'] = False
            grp.attrs['hessian_validity'] = np.bytes_(text.encode('utf-8'))
        self.database.commit(writer, note='battery of a stored Hessian')

    def hessian_validity_summary(self):
        """Counts and worst values of the battery over this generator's Hessians (the validation record, the report)."""
        log = self.hessian_validity_log

        def worst(key, fn=max):
            values = [e[key] for e in log if e.get(key) is not None]
            return fn(values) if values else None

        return {'hessians': len(log),
                'valid_at_first_attempt': sum(1 for e in log if e['valid'] and not e['recomputed'] and not e['checked_on_reuse']),
                'recomputed': sum(1 for e in log if e['recomputed']), 'checked_on_reuse': sum(1 for e in log if e['checked_on_reuse']),
                'stopped': sum(1 for e in log if e['stopped']),
                'worst': {'V0_energy': worst('V0_energy'), 'V0_gradient': worst('V0_gradient'), 'V1': worst('V1'), 'V2': worst('V2'),
                          'V3': worst('V3'), 'V4': worst('V4', min), 'V5_ratio': worst('V5_ratio'), 'V5_max': worst('V5_max')},
                'sosc_bifurcation': sum(1 for e in log if e.get('sosc') == 'bifurcation'),
                'battery_single_points': int(self.counts.get('battery_single_points', 0))}

    def file_validity_summary(self):
        """
        The battery over every Hessian record of the file (review of 2026-10-13: the generator's log covers its session,
        not the file): counts of valid, recomputed and checked-on-reuse records, those without a result, the worst values
        and the SOSC verdicts, from the raw records' attributes.
        """

        if self.database is None:
            return None
        rows = []
        with self.database.read() as h5:
            for cid in h5['raw']:
                grp = h5['raw'][cid]
                if 'hessian' not in grp:
                    continue
                text = grp.attrs.get('hessian_validity')
                if text is None:
                    rows.append(None)
                    continue
                result = json.loads(text.decode('utf-8') if isinstance(text, bytes) else str(text))
                entry = battery_summary(result)
                entry.update({'recomputed': bool(grp.attrs.get('hessian_recomputed', False)),
                              'checked_on_reuse': bool(result.get('checked_on_reuse', False)),
                              'sosc': (result.get('sosc') or {}).get('verdict')})
                rows.append(entry)
        done = [r for r in rows if r is not None]

        def worst(key, fn=max):
            values = [r[key] for r in done if r.get(key) is not None]
            return fn(values) if values else None

        return {'hessian_records': len(rows), 'without_result': len(rows) - len(done), 'valid': sum(1 for r in done if r['valid']),
                'recomputed': sum(1 for r in done if r['recomputed']), 'checked_on_reuse': sum(1 for r in done if r['checked_on_reuse']),
                'sosc': {v: sum(1 for r in done if r['sosc'] == v) for v in sorted({str(r['sosc']) for r in done})},
                'worst': {'V0_energy': worst('V0_energy'), 'V0_gradient': worst('V0_gradient'), 'V1': worst('V1'), 'V2': worst('V2'),
                          'V3': worst('V3'), 'V4': worst('V4', min), 'V5_ratio': worst('V5_ratio'), 'V5_max': worst('V5_max')}}

    def _check_hessian(self, cid, x, hessian, kind):
        """
        The M2 check: one central difference of the gradient along a
        rigid-body-free unit direction d (seeded by the calc_id) against the
        symmetrised Hessian times d, the relative error |fd - H d| / |H d|.
        The two gradients are raw records of kind hessian_check. Since M3B
        (decision E4) the value is recorded on the Hessian's record
        (hessian_fd_error) and the validity battery decides; it no longer
        stops the construction by itself.
        """

        rng = np.random.default_rng(int(cid[:8], 16))
        d = rng.normal(size=x.size)
        d = rigid_body_projector(x) @ d
        d /= np.linalg.norm(d)
        h = float(self.hessian_check_step)
        _, g_plus = self.energy_and_gradient(x + h * d.reshape(-1, 3), kind='hessian_check', note=f'Hessian check of {cid[:12]}, +')
        _, g_minus = self.energy_and_gradient(x - h * d.reshape(-1, 3), kind='hessian_check', note=f'Hessian check of {cid[:12]}, -')
        fd = (np.asarray(g_plus, dtype=np.float64) - np.asarray(g_minus, dtype=np.float64)).reshape(-1) / (2.0 * h)
        hd = 0.5 * (hessian + hessian.T) @ d
        error = float(np.linalg.norm(fd - hd) / max(np.linalg.norm(hd), 1.0e-300))
        self.hessian_checks.append({'calc_id': cid, 'kind': kind, 'error': error})
        return error

    def source(self, x):
        """
        The data source of GroupedChartModel.build: E, G and the symmetrised
        Hessian at a converged node, after the pre-flight of the node
        geometry (D2). Under the relaxed policy the node is a constrained
        stationary point, so its battery includes the constrained soft mode.
        """

        self.preflight(x, label='node')
        energy, gradient = self.energy_and_gradient(x)
        hessian = self.hessian(x, constrained=(self._cut_policy == 'relaxed_rotor_v1'))
        return energy, gradient, 0.5 * (hessian + hessian.T)

    # ------------------------------------------------------------ pre-flight

    def preflight(self, x, label='the geometry', stop=True):
        """
        D2: the condition sigma_1/sigma_r and the rank gap of the
        complete-coordinate split and the legacy cut-off margin of B B^T at
        x, with the factors to their limits; a factor below 10 stops the
        construction (critical error) when `stop`. Returns the record.
        """

        coords = self.coordinates
        policy = coords.policy
        x = np.asarray(x, dtype=np.float64).reshape(-1, 3)
        n_free = 3 * coords.n_atoms - 6
        split = coords.complete_transform(x, np.zeros(3 * coords.n_atoms), None, check=False, label=label)
        diag = split['diagnostics']
        _, b_mat = coords.ic.values_and_b(x)
        sv = np.linalg.svd(b_mat @ b_mat.T, compute_uv=False)
        sigma_r = float(sv[n_free - 1])
        sigma_r1 = float(sv[n_free]) if len(sv) > n_free else 0.0
        margin = min(sigma_r / 1.0e-6, (1.0e-6 / sigma_r1) if sigma_r1 > 0.0 else np.inf)
        factors = {'condition': policy.split_condition_max / max(diag['condition'], 1.0e-300),
                   'gap': diag['gap'] / policy.split_rank_gap_min,
                   'legacy_margin': margin / policy.legacy_cutoff_margin}
        record = {'label': str(label), 'condition': float(diag['condition']), 'gap': float(diag['gap']),
                  'rank': int(diag['rank']), 'legacy_margin': float(margin), 'factors': {k: float(v) for k, v in factors.items()},
                  'rows': diag['rows']}
        self.preflight_records.append(record)
        worst = min(factors.values())
        if stop:
            assert_msg_critical(
                worst >= self.preflight_factor,
                f'grouped MSI construction: the pre-flight at {label} leaves a factor {worst:.1f} below the limits '
                f'(condition {diag["condition"]:.2e} against {policy.split_condition_max:.0e}, gap {diag["gap"]:.2e} against '
                f'{policy.split_rank_gap_min:.0e}, legacy margin {margin:.2e} against {policy.legacy_cutoff_margin:g}); '
                f'the construction stops (PROMPT-M2 section 5.2).')
        return record

    # ------------------------------------------------------------ anchor

    @staticmethod
    def _assert_minimum(h_proj, policy, what):
        """Stops when the rigid-body-projected Hessian has an eigenvalue below -relaxation_curvature_min."""
        lam_min = float(np.linalg.eigvalsh(h_proj)[0])
        assert_msg_critical(lam_min > -float(policy.relaxation_curvature_min),
                            f'grouped MSI construction: {what} (projected Hessian eigenvalue {lam_min:.2e} hartree/bohr^2, below '
                            f'-relaxation_curvature_min). geomeTRIC keeps the symmetry of its start and can stop at a saddle (for '
                            f'example methyls staggered by symmetry); start from a geometry near the minimum, or pass an optimised '
                            f'anchor.')

    def optimise_anchor(self, molecule, polish_threshold=1.0e-6, max_polish=12, anchor=None, policy=None, max_step_bohr=0.05):
        """
        geomeTRIC optimisation at the level used (skipped when `anchor` is
        given: an optimised geometry), then a projected Newton polish on the
        QM gradient with the QM Hessian until the projected gradient is below
        the smaller of `polish_threshold` and a tenth of the policy's
        `relaxation_gradient_threshold` (at most `max_polish` steps, the
        Hessian refreshed every four). The polish refines inside the
        minimum's basin: it stops with an error when the projected Hessian has
        an eigenvalue below -`relaxation_curvature_min` (the start is not near
        a minimum), and a step whose largest Cartesian component exceeds
        `max_step_bohr` is scaled down to it. The same curvature check is
        applied to the final anchor, polished or not (geomeTRIC can stop at a
        saddle from a symmetric start). Every calculation is committed.
        Returns (x_anchor, energy, gradient, raw Hessian, record).
        """

        from .optimizationdriver import OptimizationDriver

        t0 = time.time()
        minimisation = None
        if anchor is None and self.qm_level is not None:
            # M3B (ADR-009): a grouped level has no geomeTRIC driver; a Cartesian quasi-Newton minimisation from the
            # molecule's geometry (its single points are not raw records, as geomeTRIC's are not), then the polish below
            x, minimisation = self._minimise_with_level(molecule, policy if policy is not None else GroupedModelPolicyV2())
        elif anchor is None:
            opt = OptimizationDriver(self.qm_driver)
            opt.ostream.mute()
            basis = self._basis_for(molecule)
            result = opt.compute(molecule) if basis is None else opt.compute(molecule, basis)
            x = np.array(result['final_molecule'].get_coordinates_in_bohr(), dtype=np.float64)
        else:
            x = np.array(anchor, dtype=np.float64).reshape(-1, 3)
        self._kind = 'anchor'
        energy, gradient = self.energy_and_gradient(x, kind='anchor', note='geomeTRIC optimum')
        hessian = self.hessian(x, kind='anchor')
        polish = 0
        proj = rigid_body_projector(x)
        g_t = proj @ gradient.reshape(-1)
        # the polish must leave the anchor a strict minimum for the relaxed policy's check (the max-norm of the gradient
        # projected on the phase constraints and the rigid-body modes below relaxation_gradient_threshold, M1-C decision
        # 4), so the target is the smaller of polish_threshold and a tenth of that threshold; with a fixed Hessian the
        # steps converge linearly only, so the Hessian is recomputed when four steps did not converge (M3: methyl
        # acetate's geomeTRIC optimum needed more than four steps)
        # (review of 2026-10-12, second pass, SP-5) the target follows the policy of the construction; the polish has
        # a curvature guard and a step cap: from a start outside the basin an undamped Newton step of 0.2 bohr made
        # xTB's SCC fail
        policy = policy if policy is not None else GroupedModelPolicyV2()
        target = min(float(polish_threshold), 0.1 * float(policy.relaxation_gradient_threshold))
        capped = 0
        while np.abs(g_t).max() > target and polish < int(max_polish):
            if polish and polish % 4 == 0:
                hessian = self.hessian(x, kind='anchor')
            h_proj = proj @ (0.5 * (hessian + hessian.T)) @ proj
            self._assert_minimum(h_proj, policy, 'the anchor polish starts outside a minimum')
            step = -np.linalg.pinv(h_proj, rcond=1.0e-8, hermitian=True) @ g_t
            largest = float(np.abs(step).max())
            if largest > float(max_step_bohr):
                step *= float(max_step_bohr) / largest
                capped += 1
            x = x + step.reshape(x.shape)
            energy, gradient = self.energy_and_gradient(x, kind='anchor', note=f'Newton polish {polish + 1}')
            proj = rigid_body_projector(x)
            g_t = proj @ gradient.reshape(-1)
            polish += 1
        if polish:
            hessian = self.hessian(x, kind='anchor')
        # the anchor must be a minimum whether or not the polish ran (M3 campaign: geomeTRIC kept the symmetry of the
        # fixture's acetone, both methyls staggered against C=O, and stopped at a second-order saddle)
        proj = rigid_body_projector(x)
        self._assert_minimum(proj @ (0.5 * (hessian + hessian.T)) @ proj, policy, 'the anchor is not a minimum')
        record = {'polish_steps': polish, 'capped_steps': capped, 'projected_gradient': float(np.abs(g_t).max()),
                  'hessian_asymmetry': float(np.abs(hessian - hessian.T).max()), 'wall_time': time.time() - t0}
        if minimisation is not None:
            record['minimisation'] = minimisation
        self._print(f'grouped MSI construction: anchor optimised ({polish} Newton polish steps, projected gradient '
                    f'{record["projected_gradient"]:.1e} hartree/bohr, Hessian asymmetry {record["hessian_asymmetry"]:.1e}, '
                    f'{record["wall_time"]:.0f} s)')
        return x, energy, gradient, hessian, record

    def _minimise_with_level(self, molecule, policy):
        """
        The anchor's minimisation at a grouped level (M3B, decision E7): BFGS
        in Cartesian coordinates on the level's energy and gradient from the
        molecule's geometry, to 0.3 x relaxation_gradient_threshold in the
        largest component. Returns (x, record).
        """

        from scipy.optimize import minimize

        x0 = np.array(molecule.get_coordinates_in_bohr(), dtype=np.float64)
        counter = {'n': 0}
        t0 = time.time()

        def fun(v):
            counter['n'] += 1
            energy, gradient, _ = self._compute_energy_gradient(v.reshape(-1, 3))
            return energy, gradient.reshape(-1)

        gtol = 0.3 * float(policy.relaxation_gradient_threshold)
        res = minimize(fun, x0.reshape(-1), jac=True, method='BFGS', options={'gtol': gtol, 'maxiter': 5000})
        record = {'optimiser': 'scipy BFGS in Cartesian coordinates', 'evaluations': counter['n'], 'success': bool(res.success),
                  'message': str(res.message), 'gradient_max': float(np.abs(res.jac).max()), 'gtol': gtol, 'wall_time': time.time() - t0}
        self.counts['anchor_minimisation'] = self.counts.get('anchor_minimisation', 0) + counter['n']
        self._print(f'grouped MSI construction: anchor minimised at {self.level["driver"]} ({counter["n"]} single points, largest '
                    f'gradient component {record["gradient_max"]:.1e}, {record["message"]})', warning=not res.success)
        # an unconverged minimisation is left to the Newton polish and the minimum guard of optimise_anchor, which judge the
        # anchor; the record keeps the message (review of 2026-10-13)
        return res.x.reshape(-1, 3), record

    def measure_precision_contract(self, x, policy):
        """
        ADR-009: for a grouped level, the gradient against central differences
        of the energy at x (the anchor), Richardson-extrapolated, at the
        level's SCC accuracy; above precision_contract_factor x
        relaxation_gradient_threshold in any component the construction stops
        (PROMPT-M3B section 9). Returns the record.
        """

        record = self.qm_level.precision_contract(self.labels, x, self.charge, self.multiplicity)
        record['limit'] = self.precision_contract_factor * float(policy.relaxation_gradient_threshold)
        record['passed'] = bool(record['max'] <= record['limit'])
        self.precision_contract_record = record
        self._print(f'grouped MSI construction: precision contract at the anchor: |g - dE/dx| max {record["max"]:.1e} (median '
                    f'{record["median"]:.1e}) hartree/bohr at SCC accuracy {record["scc_accuracy"]:g}, limit {record["limit"]:.1e}')
        assert_msg_critical(
            record['passed'],
            f'grouped MSI construction: the gradient of {self.level["driver"]} at SCC accuracy {record["scc_accuracy"]:g} differs from '
            f'the derivative of its energy by {record["max"]:.2e} hartree/bohr at the anchor (coordinate {record["argmax"]}), above '
            f'{record["limit"]:.1e} = {self.precision_contract_factor:g} x relaxation_gradient_threshold: the precision contract of '
            f'ADR-009 fails and the construction stops (PROMPT-M3B section 9).')
        return record

    # ------------------------------------------------------------ construction

    def coordinate_model(self, molecule, x_anchor, fragment_bonds=(), excluded_bonds=(), policy=None, domains=None,
                         relaxation_zones=None):
        """
        The GroupedCoordinateModel of the molecule on the anchor geometry:
        the Z-matrix and MM-rotatable bonds of the legacy topology (no QM),
        the group rules of ADR-003 with the user's selection. Also returns
        the Z-matrix. Must run in a directory without im_database_*.h5.
        """

        self.labels = list(molecule.get_labels())
        self.charge = int(molecule.get_charge())
        self.multiplicity = int(molecule.get_multiplicity())
        policy = policy if policy is not None else GroupedModelPolicyV2()
        z_matrix, rotatable = legacy_topology(molecule)
        _, candidates, _ = propose_groups(self.labels, z_matrix['bonds'], rotatable)
        selected = select_groups(candidates, fragment_bonds=fragment_bonds, excluded_bonds=excluded_bonds, ostream=self.ostream)
        coords = GroupedCoordinateModel(self._molecule(x_anchor), z_matrix, x_anchor, selected, policy=policy,
                                        domains=domains, relaxation_zones=relaxation_zones, ostream=self.ostream)
        return coords, z_matrix

    def construct(self, molecule, retained_terms, node_grids, cut_policy='relaxed_rotor_v1', response_policy='quadratic_v1',
                  node_orders=None, support_factor=None, weight_scheme=None, fragment_bonds=(), excluded_bonds=(),
                  policy=None, domains=None, relaxation_zones=None, root='model', directory='.', restart=True,
                  investigate=True, anchor=None, selection=None):
        """
        Builds and stores the model of one core chart. `node_grids` maps a
        group's bond (0-based sorted pair) or its group id to its phase
        offsets in radians. On a restart (the file exists and verifies, with
        the same identity) every raw record with a valid calc_id is reused.
        `anchor`: a known optimised anchor geometry (bohr), which skips the
        geomeTRIC optimisation (the Newton polish still runs).
        `retained_terms=None` selects the hypergraph by the coupling
        selection of ADR-008 (M3): the singletons are built first, the
        evidence and the probes follow (raw records of kind 'probe'), the
        selected hypergraph is built from the same cut records, and the
        selection record is stored under /derived/selection; `selection`
        overrides SELECTION_PARAMETERS. Returns (model, database).
        """

        directory = os.path.abspath(directory)
        assert_no_legacy_database(directory)
        path = os.path.join(directory, grouped_file_name(root))
        policy = policy if policy is not None else GroupedModelPolicyV2()
        self.labels = list(molecule.get_labels())
        self.charge = int(molecule.get_charge())
        self.multiplicity = int(molecule.get_multiplicity())
        t_start = time.time()
        automatic = retained_terms is None
        self.selection_parameters = selection_parameters(**(selection or {})) if automatic else None
        # the battery log covers this construction (review of 2026-10-13); the file's summary is read from its records
        self.hessian_validity_log = []
        self.counts['battery_single_points'] = 0
        try:
            return self._construct(molecule, retained_terms, node_grids, cut_policy, response_policy, node_orders, support_factor,
                                   weight_scheme, fragment_bonds, excluded_bonds, policy, domains, relaxation_zones, path, restart,
                                   investigate, anchor, automatic, t_start)
        finally:
            # D7: the relaxation steps of an unfinished node are written when the construction stops for any reason
            self.flush_relaxation_buffer(note='relaxation steps at the end of the construction')

    def _construct(self, molecule, retained_terms, node_grids, cut_policy, response_policy, node_orders, support_factor,
                   weight_scheme, fragment_bonds, excluded_bonds, policy, domains, relaxation_zones, path, restart, investigate,
                   anchor, automatic, t_start):

        self.database = GroupedInterpolationDatabase(path, ostream=self.ostream)
        existing = self.database.exists
        self._cut_policy = cut_policy
        self._frozen_rows = []
        if existing:
            assert_msg_critical(restart, f'grouped MSI construction: {path} exists; pass restart=True to continue it.')
            self.database.verify()
            self._print(f'grouped MSI construction: continuing {path} (commit {self.database.root_attributes()["commit"]}, '
                        f'{len(self.database.raw_ids())} raw records)')
            registry = self.database.registry()
            anchor_record = self.database.raw_record(registry['anchor_calc_id'])
            x_anchor = anchor_record['coordinates']
            self.coordinates, z_matrix = self.coordinate_model(molecule, x_anchor, fragment_bonds, excluded_bonds, policy,
                                                               domains, relaxation_zones)
            energy_a, gradient_a = anchor_record['energy'], anchor_record['gradient']
            hessian_a = anchor_record['hessian']
            anchor_calc = registry['anchor_calc_id']
            anchor_info = {'polish_steps': None, 'projected_gradient': None, 'hessian_asymmetry': float(np.abs(hessian_a - hessian_a.T).max()),
                           'wall_time': 0.0}
            self.precision_contract_record = registry.get('precision_contract')
        else:
            # the anchor is optimised before the file exists; its calculations are committed once the file is created
            # with the registry (the anchor's calc_id is part of it), so they are cached here first
            self.database = None
            self._kind = 'anchor'
            x_anchor, energy_a, gradient_a, hessian_a, anchor_info = self.optimise_anchor(molecule, anchor=anchor, policy=policy)
            if self.qm_level is not None and self.precision_contract:
                self.measure_precision_contract(x_anchor, policy)
            self.coordinates, z_matrix = self.coordinate_model(molecule, x_anchor, fragment_bonds, excluded_bonds, policy,
                                                               domains, relaxation_zones)
            anchor_calc = calc_id(self.level_hash, x_anchor)
        coords = self.coordinates
        bond_of = {tuple(sorted(g.rotor_bond)): gid for gid, g in enumerate(coords.groups)}
        grids = {}
        for key, grid in node_grids.items():
            gid = bond_of[tuple(sorted(key))] if isinstance(key, (tuple, list)) else int(key)
            grids[gid] = np.asarray(grid, dtype=np.float64)
        singletons = [(g,) for g in range(coords.n_groups)]
        if automatic:
            # ADR-006 amendment 9: the identity carries the selection parameters in place of typed terms; the selected
            # terms are written into the registry once known (a restart reuses them)
            terms = list(singletons)
            identity_terms = 'selected'
        else:
            terms = sorted({canonical_group_key(t) for t in retained_terms}, key=lambda t: (len(t), t))
            identity_terms = terms
        extra = {'node_grids': {str(g): grids[g].tolist() for g in sorted(grids)},
                 'cut_policy': cut_policy, 'response_policy': response_policy,
                 'weight_scheme': policy.inner_weight_scheme if weight_scheme is None else weight_scheme,
                 'support_factor': policy.inner_support_factor if support_factor is None else float(support_factor)}
        if automatic:
            extra['selection'] = dict(self.selection_parameters)
            extra['validation'] = dict(self.validation)
        identity, model_hash = model_identity(self.labels, self.charge, self.multiplicity, self.level, z_matrix, policy,
                                              coords.groups, {'reference_policy': policy.reference_policy},
                                              [] if automatic else terms, extra=extra)
        if automatic:
            identity['retained_terms'] = identity_terms
        if existing:
            assert_msg_critical(self.database.root_attributes()['model_identity'] == model_hash,
                                f'grouped MSI construction: the identity of {path} differs from the requested model; '
                                'a different model needs its own file.')
            if automatic and self.database.has_record('derived/selection'):
                stored = self.database.read_record('derived/selection')
                terms = sorted({canonical_group_key(t) for t in stored['retained_terms']}, key=lambda t: (len(t), t))
                self.selection_record = stored
                self._print(f'grouped MSI construction: the selection of {path} is reused: retained hypergraph {terms}')
        else:
            registry = registry_record(coords, terms if not automatic else [], grids, anchor_calc,
                                       extra={'anchor': anchor_info, 'level_hash': self.level_hash,
                                              'retained_terms_selected': automatic,
                                              'precision_contract': self.precision_contract_record})
            self.database = GroupedInterpolationDatabase(path, ostream=self.ostream)
            self.database.create(identity, model_hash, policy,
                                 {'labels': self.labels, 'charge': self.charge, 'multiplicity': self.multiplicity, 'z_matrix': z_matrix},
                                 self.level, registry, note='created')
            self._flush_pending()
            if not (self.database.has_raw(anchor_calc) and self.database.raw_record(anchor_calc)['hessian'] is not None):
                # the anchor was computed for an earlier file of this campaign: written here from the session memory
                self.energy_and_gradient(x_anchor, kind='anchor', note='optimised anchor')
                self.hessian(x_anchor, kind='anchor')
            assert_msg_critical(self.database.has_raw(anchor_calc) and self.database.raw_record(anchor_calc)['hessian'] is not None,
                                'grouped MSI construction: the anchor raw record was not written.')

        # the reference and its anchor jet (SPEC S13 step 1), before any bank exists
        hessian_sym = 0.5 * (hessian_a + hessian_a.T)
        self.reference = GroupedReference(coords, energy_a, gradient_a, hessian_sym)
        jet = self.reference.evaluate(x_anchor, order=2)
        proj = rigid_body_projector(x_anchor)
        jet_record = {
            'energy': float(abs(jet['energy'] - energy_a)),
            'gradient_projected': float(np.abs(jet['gradient'].reshape(-1) - proj @ gradient_a.reshape(-1)).max()),
            'gradient_rigid_body_part': float(np.linalg.norm(gradient_a.reshape(-1) - proj @ gradient_a.reshape(-1))),
            'hessian_projected': float(np.abs(proj @ (jet['hessian'] - hessian_sym) @ proj).max()),
            'hessian_asymmetry': float(np.abs(hessian_a - hessian_a.T).max()),
        }
        assert_msg_critical(jet_record['energy'] < 1.0e-10 and jet_record['gradient_projected'] < 1.0e-10
                            and jet_record['hessian_projected'] < 1.0e-8,
                            f'grouped MSI construction: the anchor jet of Q_a is not reproduced: {jet_record}.')
        self.anchor_jet = jet_record
        self._print(f'grouped MSI construction: anchor jet reproduced (energy {jet_record["energy"]:.1e}, projected gradient '
                    f'{jet_record["gradient_projected"]:.1e}, projected Hessian {jet_record["hessian_projected"]:.1e}; rigid-body part of '
                    f'the QM gradient {jet_record["gradient_rigid_body_part"]:.1e}, Hessian asymmetry {jet_record["hessian_asymmetry"]:.1e})')

        # pre-flight at the anchor and at every rigid node geometry (relaxed nodes are checked in `source`)
        self.preflight_records = []
        self.preflight(x_anchor, label='the anchor')
        for term in (singletons if automatic else terms):
            for gid in term:
                for off in grids[gid]:
                    if abs(off) < 1.0e-12:
                        continue
                    offsets = np.zeros(coords.n_groups)
                    offsets[gid] = off
                    self.preflight(coords.rigid_cut(offsets), label=f'the rigid cut of group {gid} at {np.degrees(off):.1f} deg')
        worst = min(min(r['factors'].values()) for r in self.preflight_records)
        self._print(f'grouped MSI construction: pre-flight at the anchor and {len(self.preflight_records) - 1} rigid node geometries: '
                    f'smallest factor to a limit {worst:.1f}')

        # the build: the relaxations use energy_and_gradient (committed per node, D7), the Hessian only at converged nodes
        self._kind = 'relaxation' if cut_policy == 'relaxed_rotor_v1' else 'singleton'
        self.node_kind = None
        t_build = time.time()
        cut_records = {}
        selection_ran = automatic and self.selection_record is None
        if selection_ran:
            # ADR-008: the singletons, the evidence and the probes, then the decision; the records are shared with the build
            t_sel = time.time()
            terms, self.selection_record, cut_records = select_couplings(
                coords, self.reference, self._node_source, grids, cut_policy=cut_policy, response_policy=response_policy,
                node_orders=node_orders, support_factor=support_factor, weight_scheme=weight_scheme,
                energy_gradient=self.energy_and_gradient, cut_records=cut_records, parameters=self.selection_parameters,
                validation=self.validation, ostream=self.ostream, probe_context=self._probe_context)
            self.flush_relaxation_buffer()
            self.selection_wall_time = time.time() - t_sel
            registry = registry_record(coords, terms, grids, anchor_calc,
                                       extra={'anchor': anchor_info, 'level_hash': self.level_hash, 'retained_terms_selected': True,
                                              'precision_contract': self.precision_contract_record})
            record = dict(self.selection_record)
            # ADR-008 item 8: the probe cuts' calc_ids (the generator knows the level hash; the geometry is the record's)
            for cand in record['candidates'].values():
                for cut in cand.get('probes', {}).get('cuts', []):
                    key = tuple(int(k) for k in cut['cut_key'])
                    if key in cut_records:
                        cut['calc_id'] = calc_id(self.level_hash, cut_records[key]['x'])

            def writer(h5):
                GroupedInterpolationDatabase.write_record(h5, 'derived/selection', record, replace=True)
                GroupedInterpolationDatabase.write_record(h5, 'registry', registry, replace=True)
            self.database.commit(writer, note='selection')
        model = GroupedChartModel.build(coords, self.reference, self._node_source, terms, grids, cut_policy=cut_policy,
                                        response_policy=response_policy, node_orders=node_orders, support_factor=support_factor,
                                        relaxation_zone=None, ostream=self.ostream, energy_gradient=self.energy_and_gradient,
                                        weight_scheme=weight_scheme, cut_records=cut_records)
        self.flush_relaxation_buffer()
        if selection_ran:
            # ADR-006 amendment 9: a probe reused as a node of a retained grid is counted as reused (one ledger row each).
            # Only when the selection ran here and shared its cut records with the build: on a restart the build reads
            # those nodes from the file and records their reuse itself (review of 2026-10-12, second pass, SP-8)
            rows = []
            for cand in self.selection_record['candidates'].values():
                for cut in cand.get('probes', {}).get('cuts', []):
                    key = tuple(int(k) for k in cut['cut_key'])
                    if key in model.node_records and 'calc_id' in cut:
                        rows.append(self._ledger_row(cut['calc_id'], 'probe', 'EGH', 0.0, True, {'note': 'probe reused as a node of a retained grid'}))
            if rows:
                self.counts['reused'] += len(rows)
                self.database.commit(lambda h5: GroupedInterpolationDatabase.append_ledger_rows(h5, rows), note='probes reused as nodes')
        self.build_wall_time = time.time() - t_build
        # the branch investigation (D3): a reject verdict stops the construction after the model is stored
        outcomes = []
        if model.branch_warnings and investigate:
            # the sweeps of the investigation are relaxations: their steps are buffered (kind 'investigation', M3) and written
            # every INVESTIGATION_FLUSH_STEPS steps and at the end, so a crash loses at most about one sweep (D7)
            steps = {'n': 0}

            def sweep_energy_gradient(x):
                steps['n'] += 1
                if steps['n'] % INVESTIGATION_FLUSH_STEPS == 0:
                    self.flush_relaxation_buffer(note='relaxation steps of the branch investigation')
                return self.energy_and_gradient(x, kind='investigation')

            self._kind = 'investigation'
            outcomes = model.investigate_branch_warnings(sweep_energy_gradient, ostream=self.ostream)
            self.flush_relaxation_buffer(note='relaxation steps of the branch investigation')
        store_model(self.database, model, note='model')
        validation_record = {'preflight': self.preflight_records, 'anchor_jet': jet_record,
                             'branch_investigations': [{k: o[k] for k in ('classification', 'verdict', 'message', 'gid', 'pair')}
                                                       for o in outcomes],
                             'counts': dict(self.counts), 'wall': dict(self.wall),
                             'build_wall_time': self.build_wall_time, 'total_wall_time': time.time() - t_start,
                             'validated_scope': 'unvalidated', 'preflight_factor': float(self.preflight_factor),
                             # M3 (ADR-004 amendment 13, D3 and D4): the build-time symmetry defect and the tolerances used
                             'validation_parameters': dict(self.validation),
                             'symmetry_defect': model.symmetry_defect_record if model.symmetry_defect_record is not None else {'status': 'not measured'},
                             'retained_terms': [list(t) for t in model.terms], 'hypergraph_selected': bool(automatic),
                             # M3B (ADR-009, E4): the battery of every Hessian of this construction and its tolerances, the
                             # level's precision contract at the anchor
                             'hessian_validity': {'parameters': dict(self.hessian_validity), 'battery_protocol': dict(BATTERY_PROTOCOL),
                                                  'enabled': bool(self.hessian_check), 'this_construction': self.hessian_validity_summary(),
                                                  'file': self.file_validity_summary()},
                             'precision_contract': self.precision_contract_record}
        self.database.commit(lambda h5: GroupedInterpolationDatabase.write_record(
            h5, 'validation/construction', validation_record, replace=True), note='construction record')
        # ADR-006 amendment 8: the full verification of every hash at the end of a construction
        self.database.verify()
        self.model = model
        rejected = [o for o in outcomes if o['verdict'] == 'reject']
        assert_msg_critical(not rejected, 'grouped MSI construction: the branch investigation rejected the model: '
                            + '; '.join(o['message'] for o in rejected))
        self._print(self.cost_summary())
        return model, self.database

    def _flush_pending(self):
        """Commits the calculations made before the file existed (one commit, in order)."""

        pending = list(self._pending.items())
        if not pending:
            return

        def writer(h5):
            for cid, entry in pending:
                if cid in h5['raw']:
                    continue
                GroupedInterpolationDatabase.add_raw(h5, cid, entry['kind'], entry['x'], entry['energy'], entry['gradient'],
                                                     entry['hessian'], attrs=entry['attrs'])
                for row in entry['rows']:
                    GroupedInterpolationDatabase.append_ledger(h5, row)

        self.database.commit(writer, note='calculations before the file existed')
        self._pending = {}

    def _kind_of_cut(self, x):
        """
        The raw-record kind of a node geometry from the number of groups
        rotated away from the anchor: anchor, singleton, pair or triple
        (M3); 'probe' inside the probe context of the selection.
        """

        if self.node_kind == 'probe':
            return 'probe'
        shift = np.abs(self.coordinates.phase_shift(np.asarray(x, dtype=np.float64).reshape(-1, 3)))
        order = int((shift > 1.0e-6).sum())
        return {0: 'anchor', 1: 'singleton', 2: 'pair'}.get(order, 'triple')

    @contextlib.contextmanager
    def _probe_context(self):
        """The raw records computed inside are the probe cuts of the selection (kind 'probe')."""
        kind = self.node_kind
        self.node_kind = 'probe'
        try:
            yield
        finally:
            self.node_kind = kind

    def _node_source(self, x):
        """The data source of the build and the selection: the node kind follows the cut's order (or 'probe')."""
        kind = self._kind
        self._kind = self._kind_of_cut(x)
        try:
            out = self.source(x)
        finally:
            self._kind = kind
        return out

    def cost_summary(self):
        ledger = self.database.ledger() if self.database is not None else []
        kinds = {}
        for row in ledger:
            key = (row['kind'], row['computed'], 'reused' if row['reused'] else 'computed')
            kinds[key] = kinds.get(key, 0) + 1
        images = sum((g.images - 1) * (len(self.model.node_grids.get(gid, [])) if hasattr(self, 'model') else 0)
                     for gid, g in enumerate(self.coordinates.groups)) if hasattr(self, 'coordinates') else 0
        lines = [f'grouped MSI construction: {self.counts["energy_gradient"]} energy+gradient and {self.counts["hessian"]} Hessian '
                 f'calculations ({self.hessian_kind}), {self.counts["reused"]} reused from the file, wall {self.wall["energy_gradient"]:.0f} s + '
                 f'{self.wall["hessian"]:.0f} s at {self.threads} threads; virtual images {images}; commits '
                 f'{len(self.database.commit_times) if self.database is not None else 0} '
                 f'(mean {np.mean(self.database.commit_times) * 1000.0 if self.database is not None and self.database.commit_times else 0.0:.1f} ms, '
                 f'max {max(self.database.commit_times) * 1000.0 if self.database is not None and self.database.commit_times else 0.0:.1f} ms)']
        for (kind, computed, status), n in sorted(kinds.items()):
            lines.append(f'   ledger: {kind:10s} {computed:3s} {status:8s} {n}')
        return '\n'.join(lines)
