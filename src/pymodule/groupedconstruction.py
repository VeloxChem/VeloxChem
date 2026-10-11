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
"""

import os
import time

import numpy as np

from .errorhandler import assert_msg_critical
from .groupedcoordinates import GroupedCoordinateModel, legacy_topology, propose_groups, rigid_body_projector, select_groups
from .groupeddatabase import (GroupedInterpolationDatabase, assert_no_legacy_database, calc_id, grouped_file_name,
                              level_hash, level_of_theory_record, registry_record, store_model)
from .groupeddefinitions import GroupedModelPolicyV2, canonical_group_key, model_identity
from .groupedmodel import GroupedChartModel
from .groupedreference import GroupedReference
from .imforcefieldgenerator import IMForceFieldGenerator
from .molecule import Molecule
from .outputstream import OutputStream
from .scfrestdriver import ScfRestrictedDriver
from .scfunrestdriver import ScfUnrestrictedDriver
from .xtbdriver import XtbDriver

PREFLIGHT_FACTOR = 10.0


class GroupedIMForceFieldGenerator:
    """
    :param qm_driver: XtbDriver, ScfRestrictedDriver or ScfUnrestrictedDriver,
        configured by the caller (functional, convergence, grid level)
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
    """

    def __init__(self, qm_driver, basis_label='def2-svp', ostream=None, numerical_hessian=None):
        self.qm_driver = qm_driver
        self.basis_label = str(basis_label)
        self.ostream = ostream
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
        self.counts = {'energy_gradient': 0, 'hessian': 0, 'reused': 0}
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

    # ------------------------------------------------------------ level of theory

    def _level_record(self):
        drv = self.qm_driver
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
        if isinstance(self.qm_driver, XtbDriver):
            return None
        from .molecularbasis import MolecularBasis
        return MolecularBasis.read(mol, self.basis_label, ostream=OutputStream(None))

    def _compute_energy_gradient(self, x):
        mol = self._molecule(x)
        basis = self._basis_for(mol)
        energy, scf_results, _ = self.generator._compute_energy(self.drivers[0], mol, basis)
        gradient = self.generator._compute_gradient(self.drivers[1], mol, basis, scf_results)
        converged = bool(getattr(self.qm_driver, 'is_converged', True)) if not isinstance(self.qm_driver, XtbDriver) else True
        return float(np.asarray(energy).reshape(-1)[0]), np.asarray(gradient, dtype=np.float64).reshape(-1, 3), converged

    def _compute_hessian(self, x):
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

    def energy_and_gradient(self, x, kind=None, note=''):
        """
        E and the Cartesian gradient at x, from the file when the geometry
        was computed before (a ledger row marks the reuse), otherwise
        computed and committed as a raw record before returning.
        """

        kind = self._kind if kind is None else kind
        x = np.asarray(x, dtype=np.float64).reshape(-1, 3)
        cid = calc_id(self.level_hash, x)
        db = self.database
        if db is not None and db.has_raw(cid):
            rec = db.raw_record(cid)
            if rec['gradient'] is not None:
                self.counts['reused'] += 1
                db.commit(lambda h5: GroupedInterpolationDatabase.append_ledger(
                    h5, self._ledger_row(cid, kind, 'EG', 0.0, True)), note='reuse')
                return rec['energy'], rec['gradient']
        if cid in self._memory:
            # computed for another file of this campaign: written into this file as a raw record, counted as reused
            mem = self._memory[cid]
            self.counts['reused'] += 1
            if db is not None:
                def writer(h5):
                    if cid not in h5['raw']:
                        GroupedInterpolationDatabase.add_raw(h5, cid, kind, x, mem['energy'], mem['gradient'], mem['hessian'],
                                                             attrs=mem['attrs'])
                    GroupedInterpolationDatabase.append_ledger(h5, self._ledger_row(cid, kind, 'EG', 0.0, True,
                                                                                    {'note': 'reused from an earlier file of this session'}))
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
        if db is not None:
            def writer(h5):
                GroupedInterpolationDatabase.add_raw(h5, cid, kind, x, energy, gradient, None, attrs=attrs)
                GroupedInterpolationDatabase.append_ledger(h5, self._ledger_row(cid, kind, 'EG', wall, False,
                                                                                {'rigid_body_gradient_norm': rigid}))
            db.commit(writer, note=f'{kind} EG')
        else:
            self._pending[cid] = {'kind': kind, 'x': x, 'energy': energy, 'gradient': gradient, 'hessian': None, 'attrs': attrs,
                                  'rows': [self._ledger_row(cid, kind, 'EG', wall, False, {'rigid_body_gradient_norm': rigid})]}
        assert_msg_critical(converged, f'grouped MSI construction: the QM calculation at calc_id {cid[:12]} did not converge.')
        return energy, gradient

    def hessian(self, x, kind=None, note=''):
        """
        The raw (unsymmetrised) QM Hessian at x, from the file when present,
        otherwise computed and added to the geometry's raw record (which
        the energy-and-gradient call created) before returning. The
        asymmetry max |H - H^T| is recorded.
        """

        kind = self._kind if kind is None else kind
        x = np.asarray(x, dtype=np.float64).reshape(-1, 3)
        cid = calc_id(self.level_hash, x)
        db = self.database
        if db is not None and db.has_raw(cid):
            rec = db.raw_record(cid)
            if rec['hessian'] is not None:
                self.counts['reused'] += 1
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
                            raw[cid].create_dataset('hessian', data=mem['hessian'])
                            for key in ('hessian_kind', 'hessian_asymmetry', 'wall_time_hessian'):
                                value = mem['attrs'][key]
                                raw[cid].attrs[key] = np.bytes_(value.encode('utf-8')) if isinstance(value, str) else value
                    else:
                        GroupedInterpolationDatabase.add_raw(h5, cid, kind, x, mem['energy'], mem['gradient'], mem['hessian'],
                                                             attrs=mem['attrs'])
                    GroupedInterpolationDatabase.append_ledger(h5, self._ledger_row(cid, kind, 'H', 0.0, True,
                                                                                    {'note': 'reused from an earlier file of this session'}))
                db.commit(writer, note='reuse from memory')
            return mem['hessian']
        t0 = time.time()
        hessian = self._compute_hessian(x)
        wall = time.time() - t0
        self.counts['hessian'] += 1
        self.wall['hessian'] += wall
        asymmetry = float(np.abs(hessian - hessian.T).max())
        fd_error = self._check_hessian(cid, x, hessian, kind) if self.hessian_check else float('nan')
        if cid not in self._memory:
            energy_m, gradient_m, converged_m = self._compute_energy_gradient(x)
            self._memory[cid] = {'energy': energy_m, 'gradient': gradient_m, 'hessian': None,
                                 'attrs': {'level_hash': self.level_hash, 'threads': self.threads, 'converged': converged_m,
                                           'rigid_body_gradient_norm': float(np.linalg.norm(
                                               gradient_m.reshape(-1) - rigid_body_projector(x) @ gradient_m.reshape(-1)))}}
        self._memory[cid]['hessian'] = hessian
        self._memory[cid]['attrs'].update({'hessian_kind': self.hessian_kind, 'hessian_asymmetry': asymmetry,
                                           'wall_time_hessian': float(wall), 'hessian_fd_error': fd_error})
        if db is None:
            if cid not in self._pending:
                energy, gradient = self.energy_and_gradient(x, kind=kind, note=note)
            entry = self._pending[cid]
            entry['hessian'] = hessian
            entry['kind'] = kind
            entry['attrs'].update({'hessian_kind': self.hessian_kind, 'hessian_asymmetry': asymmetry, 'wall_time_hessian': float(wall),
                                   'hessian_fd_error': fd_error})
            entry['rows'].append(self._ledger_row(cid, kind, 'H', wall, False,
                                                  {'hessian_asymmetry': asymmetry, 'hessian_fd_error': fd_error}))
            return hessian
        if db is not None:
            def writer(h5):
                raw = h5['raw']
                if cid in raw:
                    grp = raw[cid]
                    assert_msg_critical('hessian' not in grp, f'grouped MSI storage: the raw record {cid} has a Hessian already.')
                    grp.create_dataset('hessian', data=hessian)
                    grp.attrs['hessian_kind'] = np.bytes_(self.hessian_kind.encode('utf-8'))
                    grp.attrs['hessian_asymmetry'] = asymmetry
                    grp.attrs['wall_time_hessian'] = float(wall)
                    grp.attrs['hessian_fd_error'] = fd_error
                    grp.attrs['kind'] = np.bytes_(kind.encode('utf-8'))
                else:
                    energy, gradient, converged = self._compute_energy_gradient(x)
                    GroupedInterpolationDatabase.add_raw(h5, cid, kind, x, energy, gradient, hessian,
                                                         attrs={'level_hash': self.level_hash, 'threads': self.threads,
                                                                'hessian_kind': self.hessian_kind, 'hessian_asymmetry': asymmetry,
                                                                'wall_time_hessian': float(wall), 'converged': converged,
                                                                'hessian_fd_error': fd_error,
                                                                'rigid_body_gradient_norm': float(np.linalg.norm(
                                                                    gradient.reshape(-1) - rigid_body_projector(x) @ gradient.reshape(-1)))})
                GroupedInterpolationDatabase.append_ledger(h5, self._ledger_row(cid, kind, 'H', wall, False,
                                                                                {'hessian_asymmetry': asymmetry,
                                                                                 'hessian_fd_error': fd_error}))
            db.commit(writer, note=f'{kind} H')
        return hessian

    def _check_hessian(self, cid, x, hessian, kind):
        """
        One central difference of the gradient along a rigid-body-free unit
        direction d (seeded by the calc_id) against the symmetrised Hessian
        times d: the relative error |fd - H d| / |H d|. The two gradients
        are raw records of kind hessian_check. An error above
        hessian_check_tolerance stops the construction before the Hessian
        is stored: the matrix does not belong to the surface at x.
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
        assert_msg_critical(
            error <= self.hessian_check_tolerance,
            f'grouped MSI construction: the {self.hessian_kind} Hessian at calc_id {cid[:12]} ({kind}) disagrees with the central '
            f'difference of the gradient along a random direction by {error:.2e} relative (step {h:g} bohr, limit '
            f'{self.hessian_check_tolerance:g}): it does not belong to the surface at this geometry (a stale SCF state or a '
            f'failed response calculation). The Hessian is not stored.')
        return error

    def source(self, x):
        """
        The data source of GroupedChartModel.build: E, G and the symmetrised
        Hessian at a converged node, after the pre-flight of the node
        geometry (D2).
        """

        self.preflight(x, label='node')
        energy, gradient = self.energy_and_gradient(x)
        hessian = self.hessian(x)
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
                worst >= PREFLIGHT_FACTOR,
                f'grouped MSI construction: the pre-flight at {label} leaves a factor {worst:.1f} below the limits '
                f'(condition {diag["condition"]:.2e} against {policy.split_condition_max:.0e}, gap {diag["gap"]:.2e} against '
                f'{policy.split_rank_gap_min:.0e}, legacy margin {margin:.2e} against {policy.legacy_cutoff_margin:g}); '
                f'the construction stops (PROMPT-M2 section 5.2).')
        return record

    # ------------------------------------------------------------ anchor

    def optimise_anchor(self, molecule, polish_threshold=1.0e-6, max_polish=4, anchor=None):
        """
        geomeTRIC optimisation at the level used, then a projected Newton
        polish on the QM gradient with the QM Hessian until the projected
        gradient is below `polish_threshold` (at most `max_polish` steps).
        Every calculation is committed. Returns (x_anchor, energy, gradient,
        raw Hessian, record).
        """

        from .optimizationdriver import OptimizationDriver

        t0 = time.time()
        if anchor is None:
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
        while np.abs(g_t).max() > polish_threshold and polish < max_polish:
            h_sym = 0.5 * (hessian + hessian.T)
            step = -np.linalg.pinv(proj @ h_sym @ proj, rcond=1.0e-8, hermitian=True) @ g_t
            x = x + step.reshape(x.shape)
            energy, gradient = self.energy_and_gradient(x, kind='anchor', note=f'Newton polish {polish + 1}')
            proj = rigid_body_projector(x)
            g_t = proj @ gradient.reshape(-1)
            polish += 1
        if polish:
            hessian = self.hessian(x, kind='anchor')
        record = {'polish_steps': polish, 'projected_gradient': float(np.abs(g_t).max()),
                  'hessian_asymmetry': float(np.abs(hessian - hessian.T).max()), 'wall_time': time.time() - t0}
        self._print(f'grouped MSI construction: anchor optimised ({polish} Newton polish steps, projected gradient '
                    f'{record["projected_gradient"]:.1e} hartree/bohr, Hessian asymmetry {record["hessian_asymmetry"]:.1e}, '
                    f'{record["wall_time"]:.0f} s)')
        return x, energy, gradient, hessian, record

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
                  investigate=True, anchor=None):
        """
        Builds and stores the model of one core chart. `node_grids` maps a
        group's bond (0-based sorted pair) or its group id to its phase
        offsets in radians. On a restart (the file exists and verifies, with
        the same identity) every raw record with a valid calc_id is reused.
        `anchor`: a known optimised anchor geometry (bohr), which skips the
        geomeTRIC optimisation (the Newton polish still runs). Returns
        (model, database).
        """

        directory = os.path.abspath(directory)
        assert_no_legacy_database(directory)
        path = os.path.join(directory, grouped_file_name(root))
        policy = policy if policy is not None else GroupedModelPolicyV2()
        self.labels = list(molecule.get_labels())
        self.charge = int(molecule.get_charge())
        self.multiplicity = int(molecule.get_multiplicity())
        t_start = time.time()

        self.database = GroupedInterpolationDatabase(path, ostream=self.ostream)
        existing = self.database.exists
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
        else:
            # the anchor is optimised before the file exists; its calculations are committed once the file is created
            # with the registry (the anchor's calc_id is part of it), so they are cached here first
            self.database = None
            self._kind = 'anchor'
            x_anchor, energy_a, gradient_a, hessian_a, anchor_info = self.optimise_anchor(molecule, anchor=anchor)
            self.coordinates, z_matrix = self.coordinate_model(molecule, x_anchor, fragment_bonds, excluded_bonds, policy,
                                                               domains, relaxation_zones)
            anchor_calc = calc_id(self.level_hash, x_anchor)
        coords = self.coordinates
        bond_of = {tuple(sorted(g.rotor_bond)): gid for gid, g in enumerate(coords.groups)}
        grids = {}
        for key, grid in node_grids.items():
            gid = bond_of[tuple(sorted(key))] if isinstance(key, (tuple, list)) else int(key)
            grids[gid] = np.asarray(grid, dtype=np.float64)
        terms = sorted({canonical_group_key(t) for t in retained_terms}, key=lambda t: (len(t), t))
        identity, model_hash = model_identity(self.labels, self.charge, self.multiplicity, self.level, z_matrix, policy,
                                              coords.groups, {'reference_policy': policy.reference_policy}, terms,
                                              extra={'node_grids': {str(g): grids[g].tolist() for g in sorted(grids)},
                                                     'cut_policy': cut_policy, 'response_policy': response_policy,
                                                     'weight_scheme': policy.inner_weight_scheme if weight_scheme is None else weight_scheme,
                                                     'support_factor': policy.inner_support_factor if support_factor is None else float(support_factor)})
        if existing:
            assert_msg_critical(self.database.root_attributes()['model_identity'] == model_hash,
                                f'grouped MSI construction: the identity of {path} differs from the requested model; '
                                'a different model needs its own file.')
        else:
            registry = registry_record(coords, terms, grids, anchor_calc,
                                       extra={'anchor': anchor_info, 'level_hash': self.level_hash})
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
        for term in terms:
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

        # the build: the relaxations use energy_and_gradient (every call committed), the Hessian only at converged nodes
        self._kind = 'relaxation' if cut_policy == 'relaxed_rotor_v1' else 'singleton'
        t_build = time.time()
        self.node_kind = 'singleton'
        model = GroupedChartModel.build(coords, self.reference, self._node_source, terms, grids, cut_policy=cut_policy,
                                        response_policy=response_policy, node_orders=node_orders, support_factor=support_factor,
                                        relaxation_zone=None, ostream=self.ostream, energy_gradient=self.energy_and_gradient,
                                        weight_scheme=weight_scheme)
        self.build_wall_time = time.time() - t_build
        # the branch investigation (D3): a reject verdict stops the construction after the model is stored
        outcomes = []
        if model.branch_warnings and investigate:
            self._kind = 'probe'
            outcomes = model.investigate_branch_warnings(self.energy_and_gradient, ostream=self.ostream)
        store_model(self.database, model, note='model')
        self.database.commit(lambda h5: GroupedInterpolationDatabase.write_record(
            h5, 'validation/construction', {'preflight': self.preflight_records, 'anchor_jet': jet_record,
                                            'branch_investigations': [{k: o[k] for k in ('classification', 'verdict', 'message', 'gid', 'pair')}
                                                                      for o in outcomes],
                                            'counts': dict(self.counts), 'wall': dict(self.wall),
                                            'build_wall_time': self.build_wall_time, 'total_wall_time': time.time() - t_start,
                                            'validated_scope': 'unvalidated'}, replace=True), note='construction record')
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

    def _node_source(self, x):
        """The data source of the build for the planned node kinds (singleton nodes; pairs in M3)."""
        kind = self._kind
        self._kind = self.node_kind
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
