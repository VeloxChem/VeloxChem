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
Coupling selection of the grouped MSI v2 model family (SPEC S8; ADR-008 of
.claude/tasks/grouped-msi/decisions; M3).

Which pair and triple increments the retained hypergraph keeps is decided
from measured evidence, every kind of which is converted into one currency:
the predicted omitted error of the candidate S, |F+(x) - F-(x)| in kcal/mol
and the largest per-atom |grad F+(x) - grad F-(x)| in kcal/mol/A (decision
D2), where F+ and F- are the models with and without the increment of S,
built from the same cut records on sparse node sets at whose nodes the
increment is exact (a tensor-product bank gives its node the weight 1).

- Evidence 1, the anchor-Hessian screening of S8 (eta_op, eta_F with units
  and scales): persisted as a signal, decides nothing.
- Evidence 2 and 3, the anchor-line model of S = T + {g} (no QM): the bank
  of S on grid(T) x {anchor of g}; its g-derivative data are the normal
  derivative of the dropped increment on g's anchor line (the missing force
  at T's nodes), and its covariant per-face data give the leak
  (Pi_g^k - 1) Delta_T X on g's image lines (the symmetry defect that
  dropping S causes, Q5 report section 2.6).
- Evidence 4, the probes (QM): the cuts of S at the product of two probe
  phases per member (the singleton nodes with the largest |I|), each an
  exact Moebius increment; the probe bank on the product grid {0, probes}.
- Every evidence point is also displaced along seeded non-rotor directions
  with E_v = E_max, which makes the response parts (b, C, w) visible.
- Rule: a candidate is retained when its predicted energy or force error
  reaches theta x the D4 grouping tolerance (seed theta 0.25); pairs first
  against the singletons, then the triples; the retained hypergraph is the
  downward closure. Deterministic: the same records give the same selection
  record bitwise.
"""

import contextlib
import time

import numpy as np

from .errorhandler import assert_msg_critical
from .groupedalgebra import downward_closure, subsets_of
from .groupedbanks import metric_distance
from .groupedcoordinates import BOND
from .groupeddefinitions import canonical_group_key, normalize_for_json, selection_parameters, validation_parameters
from .groupedmodel import KCAL, KCAL_PER_ANGSTROM, TWO_PI, CutBuilder, GroupedChartModel, cut_key, term_label

KB_HARTREE_PER_KELVIN = 3.166811563e-6


def _max_atom_norm(vectors):
    return float(np.linalg.norm(np.asarray(vectors, dtype=np.float64).reshape(-1, 3), axis=1).max())


class CouplingSelection:
    """
    The selection of ADR-008 on one core chart.

    :param coordinates: the GroupedCoordinateModel
    :param reference: the GroupedReference Q_a (anchor gradient and Hessian)
    :param builder: a CutBuilder (its records hold the singleton cuts, and
        receive the probe cuts)
    :param node_grids: dict gid -> singleton grid (phase offsets, radians,
        containing 0) of every group
    :param support_factor, weight_scheme: of the banks (the values the final
        model uses)
    :param parameters: overrides of SELECTION_PARAMETERS
    :param validation: overrides of VALIDATION_PARAMETERS (the grouping
        tolerances set the thresholds)
    :param probe_context: a callable returning a context manager that is
        entered around the QM of every probe cut (the construction switches
        the raw-record kind to 'probe' in it); None: no context
    """

    def __init__(self, coordinates, reference, builder, node_grids, support_factor, weight_scheme,
                 parameters=None, validation=None, ostream=None, probe_context=None):
        self.coords = coordinates
        self.reference = reference
        self.builder = builder
        self.records = builder.records
        self.node_grids = {int(g): np.asarray(grid, dtype=np.float64) for g, grid in node_grids.items()}
        self.support_factor = float(support_factor)
        self.weight_scheme = str(weight_scheme)
        self.node_orders = dict(builder.node_orders)
        self.cut_policy = builder.cut_policy
        self.response_policy = builder.response_policy
        self.params = selection_parameters(**(parameters or {}))
        self.validation = validation_parameters(**(validation or {}))
        self.ostream = ostream
        self.probe_context = probe_context if probe_context is not None else (lambda: contextlib.nullcontext())
        self.theta = float(self.params['selection_threshold_factor'])
        self.tol_e = float(self.validation['grouping_error_energy_max_kcal'])
        self.tol_f = float(self.validation['grouping_error_force_max_kcal_per_angstrom'])
        self.tau_e = self.theta * self.tol_e
        self.tau_f = self.theta * self.tol_f
        self.e_max = float(coordinates.policy.response_domain_energy_max)
        self.e_tol = self.tol_e / KCAL
        n_g = coordinates.n_groups
        self.groups = list(range(n_g))
        self.singletons = [(g,) for g in self.groups]
        for g in self.groups:
            assert_msg_critical(g in self.node_grids, f'grouped MSI selection: no singleton grid for group {g}.')
        self._directions = None
        self.record = {
            'parameters': dict(self.params),
            'thresholds': {'energy_kcal': self.tau_e, 'force_kcal_per_angstrom': self.tau_f,
                           'grouping_tolerance_energy_kcal': self.tol_e, 'grouping_tolerance_force_kcal_per_angstrom': self.tol_f},
            'e_max': self.e_max, 'hessian_energy_tolerance_hartree': self.e_tol,
            'cut_policy': self.cut_policy, 'response_policy': self.response_policy, 'node_orders': dict(self.node_orders),
            'support_factor': self.support_factor, 'weight_scheme': self.weight_scheme,
            'node_grids': {str(g): self.node_grids[g].tolist() for g in sorted(self.node_grids)},
            'statistic': 'maximum over the evidence points of |dE| (kcal/mol) and of the largest per-atom |dF| (kcal/mol/A)',
            'candidates': {}, 'probe_phases': {}, 'hessian_screening': {}, 'retained_terms': None,
        }

    def _print(self, message, warning=False):
        if self.ostream is not None:
            (self.ostream.print_warning if warning else self.ostream.print_info)(message)
            self.ostream.flush()

    # --------------------------------------------------------- the records

    def _offsets(self, assignment):
        """The full offset vector with the given gid -> offset entries."""
        off = np.zeros(self.coords.n_groups)
        for gid, value in assignment.items():
            off[int(gid)] = float(value)
        return off

    def ensure_singletons(self):
        """Every singleton node (and the anchor) as a cut record."""
        self.builder.record(self._offsets({}))
        for g in self.groups:
            for off in self.node_grids[g]:
                self.builder.record(self._offsets({g: off}))

    def model_from(self, terms, grids_by_term):
        """
        A GroupedChartModel on `terms` from the shared records: every term
        takes its member grids from grids_by_term[term] (a dict gid -> grid),
        singletons their full grids by default. The banks of different
        terms may use different grids (the sparse banks of the evidence);
        every bank passes the anchor-line check of ADR-004 amendment 12
        against the scale of its cut data.
        """

        terms = sorted({canonical_group_key(t) for t in terms}, key=lambda t: (len(t), t))
        grids_of = {}
        for term in terms:
            grids = {g: self.node_grids[g] for g in term}
            grids.update({int(g): np.asarray(v, dtype=np.float64) for g, v in grids_by_term.get(term, {}).items()})
            grids_of[term] = grids
        # the sparse grids of the evidence banks ({0} plus two probe phases) can have a spacing whose product with the
        # model's support factor reaches the wrap of a member's metric; the evidence models therefore use one common
        # factor that keeps every member radius below the wrap. The predictions are taken at bank nodes, where the
        # node's weight is 1 for any factor, so the factor does not enter them
        factor = self.safe_support_factor(grids_of)
        banks = {}
        for term in terms:
            banks.update(GroupedChartModel.banks_from_cuts(self.coords, self.records, [term], grids_of[term], self.node_orders,
                                                           self.coords.policy, factor, self.weight_scheme))
        model = GroupedChartModel(self.coords, self.reference, banks, self.response_policy, self.cut_policy,
                                  node_records=self.records, weight_scheme=self.weight_scheme, support_factor=factor,
                                  frozen_rows=self.builder.frozen)
        model.node_orders = dict(self.node_orders)
        model.node_grids = dict(self.node_grids)
        return model

    def safe_support_factor(self, grids_of):
        """
        The largest factor not above the model's that keeps c x spacing_i
        below 2 / m_i for every member grid of every term in `grids_of`
        (term -> {gid: grid}); the spacing is the smallest metric distance
        between two grid values of the member (as IncrementBank defines it).
        """

        factor = self.support_factor
        for term, grids in grids_of.items():
            for g in term:
                m = int(self.coords.groups[g].images)
                values = np.asarray(grids[g], dtype=np.float64)
                if len(values) < 2:
                    continue
                spacing = min(metric_distance(values[a] - values[b], m) for a in range(len(values)) for b in range(len(values)) if a != b)
                factor = min(factor, 0.999 * (2.0 / m) / spacing)
        return float(factor)

    # ------------------------------------------------- evidence 1: Hessian

    def hessian_screening(self):
        """
        SPEC S8 and ADR-008 item 3: K_gh = Gamma[Y_g, Y_h] with Gamma the anchor
        Hessian of U in the complete coordinates, Y_g the phase of g and its
        local non-rotor rows, the scales C_g = diag(s^2) (the full period for
        the phase, sqrt(k_B T / Gamma_kk) for the rows), eta_op = |Gamma_gh|_2 /
        E_tol, eta_F = |Gamma_gh|_F / (sqrt(d_g d_h) E_tol). A signal only.
        """

        coords = self.coords
        split = coords.complete_transform(coords.x_anchor, self.reference.gradient_anchor, self.reference.hessian_anchor,
                                          check=False, label='the anchor')
        gamma = split['Gamma']
        n_g = coords.n_groups
        position = {int(r): i for i, r in enumerate(coords.response_rows)}
        k_t = KB_HARTREE_PER_KELVIN * float(self.params['selection_temperature_kelvin'])
        coords_of = {}
        for g in self.groups:
            rows = [('phase', g, int(g))]
            for r in coords.local_rows(g):
                rows.append(('row', int(r), n_g + position[int(r)]))
            entries = []
            for kind, ident, index in rows:
                if kind == 'phase':
                    scale = TWO_PI / coords.groups[g].images
                    entries.append({'name': f'phase of group {g}', 'unit': 'rad', 'scale': float(scale), 'index': int(index),
                                    'curvature': float(gamma[index, index]), 'flag': ''})
                else:
                    curv = float(gamma[index, index])
                    flag = '' if curv > 0.0 else 'non-positive diagonal curvature, |value| floored at 1e-8'
                    scale = float(np.sqrt(k_t / max(abs(curv), 1.0e-8)))
                    unit = 'bohr^-1' if coords.ic.kinds[ident] == BOND else 'rad'
                    entries.append({'name': coords.describe_row(ident), 'unit': unit, 'scale': scale, 'index': int(index),
                                    'curvature': curv, 'flag': flag})
            coords_of[g] = entries
        out = {'temperature_kelvin': float(self.params['selection_temperature_kelvin']), 'k_T_hartree': k_t,
               'energy_tolerance_hartree': self.e_tol, 'group_coordinates': {str(g): coords_of[g] for g in self.groups}, 'pairs': {}}
        for g, h in self._pairs():
            ig = [e['index'] for e in coords_of[g]]
            ih = [e['index'] for e in coords_of[h]]
            s_g = np.array([e['scale'] for e in coords_of[g]])
            s_h = np.array([e['scale'] for e in coords_of[h]])
            k_block = gamma[np.ix_(ig, ih)]
            scaled = (s_g[:, None] * k_block) * s_h[None, :]
            sv = np.linalg.svd(scaled, compute_uv=False)
            out['pairs'][term_label((g, h))] = {
                'd_g': len(ig), 'd_h': len(ih), 'K_norm_2': float(np.linalg.svd(k_block, compute_uv=False)[0]),
                'Gamma_norm_2_hartree': float(sv[0]), 'Gamma_norm_F_hartree': float(np.linalg.norm(scaled)),
                'eta_op': float(sv[0] / self.e_tol), 'eta_F': float(np.linalg.norm(scaled) / (np.sqrt(len(ig) * len(ih)) * self.e_tol)),
                'phase_phase_curvature_hartree': float(gamma[g, h]),
            }
        self.record['hessian_screening'] = out
        return out

    # ----------------------------------------------- candidates and probes

    def _pairs(self):
        return [(g, h) for i, g in enumerate(self.groups) for h in self.groups[i + 1:]]

    def _triples(self, retained_pairs):
        n_g = len(self.groups)
        triples = [(g, h, k) for i, g in enumerate(self.groups) for j, h in enumerate(self.groups[i + 1:], i + 1)
                   for k in self.groups[j + 1:]]
        if n_g <= int(self.params['selection_max_groups_all_triples']):
            return triples
        retained = {canonical_group_key(p) for p in retained_pairs}
        return [t for t in triples if all(canonical_group_key(p) in retained for p in subsets_of(t, 2, 2))]

    def probe_phases(self):
        """
        ADR-008 item 5: per group the singleton nodes (offset != 0) with the
        largest |I|, ties broken by the smaller offset.
        """

        n = int(self.params['selection_probes_per_group'])
        out = {}
        for g in self.groups:
            candidates = []
            for off in self.node_grids[g]:
                if abs(off) < 1.0e-12:
                    continue
                rec = self.records[cut_key(self._offsets({g: off}))]
                candidates.append((-abs(float(rec['data']['I'])), float(off)))
            assert_msg_critical(candidates, f'grouped MSI selection: the grid of group {g} has no node besides the anchor.')
            candidates.sort()
            chosen = [off for _, off in candidates[:n]]
            out[g] = sorted(chosen)
            self.record['probe_phases'][str(g)] = [
                {'offset_rad': off, 'offset_deg': float(np.degrees(off)), 'abs_I_hartree': -neg} for neg, off in candidates[:n]]
        self.probes = out
        return out

    def probe_grid(self, g):
        """{0} plus the probe phases of g, as a sorted array."""
        return np.array(sorted([0.0] + list(self.probes[g])), dtype=np.float64)

    def compute_probes(self, term):
        """The probe cuts of a term: the product of its members' probe phases (QM through the builder)."""
        members = list(term)
        mesh = np.array(np.meshgrid(*[self.probes[g] for g in members], indexing='ij')).reshape(len(members), -1).T
        ids = []
        for node in mesh:
            off = self._offsets({g: v for g, v in zip(members, node)})
            with self.probe_context():
                rec = self.builder.record(off)
            ids.append({'offsets_deg': np.degrees(off).round(6).tolist(), 'energy': float(rec['energy']),
                        'relaxed': rec.get('relaxation') is not None, 'cut_key': list(cut_key(off))})
        return ids

    # ---------------------------------------------------- evidence points

    def directions(self):
        """The seeded response directions with E_v = E_max (ADR-008 item 6), the same at every point."""
        if self._directions is None:
            rng = np.random.default_rng(int(self.params['selection_seed']))
            probe_model = self.model_from(self.singletons, {})
            n = int(self.params['selection_response_displacements'])
            self._directions = probe_model.response_directions([1.0] * n, self.e_max, rng)
        return self._directions

    def _with_displacements(self, x, label):
        """The point itself and its displaced variants (ADR-008 item 6), with the realisation residuals."""
        points = [{'x': x, 'label': label, 'displaced': False}]
        for i, dq in enumerate(self.directions()):
            if dq is None:
                continue
            y, info = self.coords.realise_response_displacement(x, dq)
            if not info['converged']:
                points.append({'x': None, 'label': f'{label}, direction {i}', 'displaced': True,
                               'dropped': f'realisation did not converge (step {info["step"]:.1e} after {info["substeps"]} sub-steps)'})
                continue
            points.append({'x': y, 'label': f'{label}, direction {i}', 'displaced': True})
        return points

    def predict(self, model_minus, model_plus, points):
        """
        ADR-008 item 2 at every point: dE = |F+ - F-| (kcal/mol), dF = the
        largest per-atom |grad F+ - grad F-| (kcal/mol/A), with E_v of F-.
        """

        out = []
        for pt in points:
            entry = {'label': pt['label'], 'displaced': bool(pt['displaced'])}
            if pt.get('x') is None:
                entry['dropped'] = pt.get('dropped', 'no geometry')
                out.append(entry)
                continue
            a = model_minus.evaluate(pt['x'], 1, enforce_coverage=False)
            b = model_plus.evaluate(pt['x'], 1, enforce_coverage=False)
            if not (a['covered'] and b['covered']):
                entry['dropped'] = 'uncovered: ' + '; '.join(a['coverage_reasons'] + b['coverage_reasons'])
                out.append(entry)
                continue
            entry['energy_kcal'] = abs(b['energy'] - a['energy']) * KCAL
            entry['force_kcal_per_angstrom'] = _max_atom_norm(b['gradient'] - a['gradient']) * KCAL_PER_ANGSTROM
            entry['E_v_hartree'] = None if a['E_v'] is None else float(a['E_v'])
            out.append(entry)
        return out

    @staticmethod
    def _maxima(points):
        e = [p['energy_kcal'] for p in points if 'energy_kcal' in p]
        f = [p['force_kcal_per_angstrom'] for p in points if 'force_kcal_per_angstrom' in p]
        return {'energy_kcal': max(e) if e else 0.0, 'force_kcal_per_angstrom': max(f) if f else 0.0,
                'n_points': len(e), 'n_dropped': sum(1 for p in points if 'dropped' in p)}

    # -------------------------------------- evidence 2 and 3: anchor lines

    def anchor_line_evidence(self, term, base_terms, grids_by_term):
        """
        ADR-008 item 4 for the candidate `term`: for every face T of order
        |S| - 1 and the missing member g, the anchor-line bank of S on
        grid(T) x {0}, evaluated against the model without S at T's nodes
        (the boundary force), at their images under g (the leak) and at the
        displaced variants of both.
        """

        coords = self.coords
        model_minus = self.model_from(base_terms, grids_by_term)
        boundary, leak = {}, {}
        sizes = {}
        for face in subsets_of(term, len(term) - 1, len(term) - 1):
            g = [h for h in term if h not in face][0]
            grids = {h: (grids_by_term.get(face, {}).get(h, self.node_grids[h])) for h in face}
            grids[g] = np.array([0.0])
            model_plus = self.model_from(list(base_terms) + [term], {**grids_by_term, term: grids})
            # T's nodes: the product of the face's grids, the anchor excluded (zero evidence there)
            mesh = np.array(np.meshgrid(*[np.asarray(grids[h], dtype=np.float64) for h in face], indexing='ij')).reshape(len(face), -1).T
            node_points, image_points = [], []
            for node in mesh:
                if np.all(np.abs(node) < 1.0e-12):
                    continue
                off = self._offsets({h: v for h, v in zip(face, node)})
                rec = self.records[cut_key(off)]
                label = f'node of {face} at {np.degrees(off[list(face)]).round(1).tolist()} deg'
                node_points.extend(self._with_displacements(rec['x'], label))
                grp = coords.groups[g]
                for k in range(1, grp.images):
                    p = coords.cyclic_permutation(g, k)
                    image_points.extend(self._with_displacements(np.asarray(rec['x'])[p], f'{label}, image {k} of group {g}'))
            key = f'{term_label(face)}|{g}'
            boundary[key] = {'face': list(face), 'group': int(g), 'points': self.predict(model_minus, model_plus, node_points)}
            boundary[key].update(self._maxima(boundary[key]['points']))
            if image_points:
                leak[key] = {'face': list(face), 'group': int(g), 'points': self.predict(model_minus, model_plus, image_points)}
                leak[key].update(self._maxima(leak[key]['points']))
                # the coefficient sizes of the leak (Pi_g^k - 1) Delta_T X against |Delta_T X| at T's nodes (Q5 2.6)
                sizes[key] = self._leak_sizes(model_plus.banks[canonical_group_key(term)], face, g)
        return boundary, leak, sizes

    def _leak_sizes(self, bank, face, g):
        """max |(Pi_g^k - 1) Delta_T X| for b, C, w over the nodes of the anchor-line bank, against max |Delta_T X|."""
        names = [n for n in ('b', 'C', 'w') if n in bank.names]
        if not names:
            return {}
        term = bank.term
        i_g = list(term).index(g)
        grp = self.coords.groups[g]
        out = {n: {'leak_max': 0.0, 'face_max': 0.0} for n in names}
        for j in range(bank.n_nodes):
            for k in range(1, grp.images):
                powers_k = [0] * bank.dim
                powers_k[i_g] = int(grp.sense * k) % grp.images
                for n in names:
                    value_k, _, _ = bank._node_arrays(j, n, tuple(powers_k))
                    out[n]['leak_max'] = max(out[n]['leak_max'], float(np.abs(value_k).max()) if np.size(value_k) else 0.0)
            # |Delta_T X| at the node: the signed sum of the faces that do not contain g, i.e. the lower increment
            entry = bank.node_data[j]
            for n in names:
                faces = entry[n]['faces'] if isinstance(entry[n], dict) else None
                if faces is None:
                    continue
                total = None
                for f_key, (sign, v, _, _) in faces.items():
                    if g in f_key:
                        continue
                    total = sign * np.asarray(v) if total is None else total + sign * np.asarray(v)
                if total is not None:
                    out[n]['face_max'] = max(out[n]['face_max'], float(np.abs(total).max()))
        return out

    # -------------------------------------------------- evidence 4: probes

    def probe_evidence(self, term, base_terms, grids_by_term):
        """
        ADR-008 item 5: the probe bank of the candidate on the product of
        {0, probes} per member, against the model without it, at the probe
        geometries and their displaced variants.
        """

        grids = {g: self.probe_grid(g) for g in term}
        model_minus = self.model_from(base_terms, grids_by_term)
        model_plus = self.model_from(list(base_terms) + [term], {**grids_by_term, term: grids})
        members = list(term)
        mesh = np.array(np.meshgrid(*[self.probes[g] for g in members], indexing='ij')).reshape(len(members), -1).T
        points = []
        for node in mesh:
            off = self._offsets({g: v for g, v in zip(members, node)})
            rec = self.records[cut_key(off)]
            points.extend(self._with_displacements(rec['x'], f'probe at {np.degrees(node).round(1).tolist()} deg'))
        out = {'grids_deg': {str(g): np.degrees(grids[g]).round(6).tolist() for g in term},
               'support_factor_used': float(model_plus.support_factor),
               'points': self.predict(model_minus, model_plus, points)}
        out.update(self._maxima(out['points']))
        return out

    # ------------------------------------------------------------- the rule

    def decide(self, candidate_record):
        """The retention rule of ADR-008 item 7: the first evidence point that reaches a threshold, in the fixed order."""
        for kind in ('boundary', 'leak', 'probes'):
            block = candidate_record.get(kind)
            if not block:
                continue
            groups = block.values() if kind != 'probes' else [block]
            for sub in groups:
                for point in sub['points']:
                    if 'energy_kcal' not in point:
                        continue
                    if point['energy_kcal'] >= self.tau_e:
                        return True, f'{kind}: {point["label"]}: energy {point["energy_kcal"]:.4g} kcal/mol >= {self.tau_e:.4g}'
                    if point['force_kcal_per_angstrom'] >= self.tau_f:
                        return True, f'{kind}: {point["label"]}: force {point["force_kcal_per_angstrom"]:.4g} kcal/mol/A >= {self.tau_f:.4g}'
        return False, 'below threshold'

    def _predicted(self, candidate_record):
        e, f = 0.0, 0.0
        for kind in ('boundary', 'leak', 'probes'):
            block = candidate_record.get(kind)
            if not block:
                continue
            for sub in (block.values() if kind != 'probes' else [block]):
                e = max(e, sub['energy_kcal'])
                f = max(f, sub['force_kcal_per_angstrom'])
        return {'energy_kcal': e, 'force_kcal_per_angstrom': f}

    def run(self):
        """
        The whole selection: singletons, the Hessian screening, the probe
        phases, every pair (free evidence, probes, decision), every triple
        candidate, the downward closure. Returns (retained_terms, record).
        """

        t0 = time.time()
        self.ensure_singletons()
        self.hessian_screening()
        self.probe_phases()
        candidates = self.record['candidates']
        retained = list(self.singletons)
        probe_grids = {}
        decided_pairs = []
        for pair in self._pairs():
            pair = canonical_group_key(pair)
            rec = {'order': 2, 'hessian': self.record['hessian_screening']['pairs'][term_label(pair)]}
            boundary, leak, sizes = self.anchor_line_evidence(pair, self.singletons, {})
            rec['boundary'], rec['leak'], rec['leak_sizes'] = boundary, leak, sizes
            cuts = self.compute_probes(pair)
            probe_grids[pair] = {g: self.probe_grid(g) for g in pair}
            rec['probes'] = self.probe_evidence(pair, self.singletons, {})
            rec['probes']['cuts'] = cuts
            rec['predicted'] = self._predicted(rec)
            keep, reason = self.decide(rec)
            rec['decision'] = 'retained' if keep else 'dropped'
            rec['reason'] = reason
            candidates[term_label(pair)] = rec
            if keep:
                retained.append(pair)
                decided_pairs.append(pair)
            self._print(f'grouped MSI selection: pair {pair} {rec["decision"]} ({reason}); predicted error '
                        f'{rec["predicted"]["energy_kcal"]:.4g} kcal/mol, {rec["predicted"]["force_kcal_per_angstrom"]:.4g} kcal/mol/A; '
                        f'eta_op {rec["hessian"]["eta_op"]:.3g}, eta_F {rec["hessian"]["eta_F"]:.3g}')
        retained_pairs = list(decided_pairs)
        for triple in self._triples(retained_pairs):
            triple = canonical_group_key(triple)
            pairs = [canonical_group_key(p) for p in subsets_of(triple, 2, 2)]
            # the triple's models carry H2 = the singletons plus every retained pair, and the probe banks of the
            # triple's three pairs whether retained or not (ADR-008 item 1), present in both models
            base_pairs = sorted({*pairs, *retained_pairs}, key=lambda t: (len(t), t))
            base = list(self.singletons) + base_pairs
            base_grids = {p: probe_grids[p] for p in base_pairs}
            etas = [self.record['hessian_screening']['pairs'][term_label(p)] for p in pairs]
            rec = {'order': 3, 'hessian': {'note': 'the anchor Hessian is bilinear and carries no triple signal; the largest pair values',
                                           'eta_op': max(e['eta_op'] for e in etas), 'eta_F': max(e['eta_F'] for e in etas)}}
            boundary, leak, sizes = self.anchor_line_evidence(triple, base, base_grids)
            rec['boundary'], rec['leak'], rec['leak_sizes'] = boundary, leak, sizes
            cuts = self.compute_probes(triple)
            rec['probes'] = self.probe_evidence(triple, base, base_grids)
            rec['probes']['cuts'] = cuts
            rec['predicted'] = self._predicted(rec)
            keep, reason = self.decide(rec)
            rec['decision'] = 'retained' if keep else 'dropped'
            rec['reason'] = reason
            candidates[term_label(triple)] = rec
            if keep:
                retained.append(triple)
                for p in pairs:
                    if p not in retained:
                        retained.append(p)
                        candidates[term_label(p)]['decision'] = f'retained as a face of {triple}'
            self._print(f'grouped MSI selection: triple {triple} {rec["decision"]} ({reason}); predicted error '
                        f'{rec["predicted"]["energy_kcal"]:.4g} kcal/mol, {rec["predicted"]["force_kcal_per_angstrom"]:.4g} kcal/mol/A')
        terms = [t for t in downward_closure(retained) if t]
        self.record['retained_terms'] = [list(t) for t in terms]
        self.record['wall_time'] = time.time() - t0
        self._print(f'grouped MSI selection: retained hypergraph {terms} (theta {self.theta:g}: thresholds {self.tau_e:.4g} kcal/mol, '
                    f'{self.tau_f:.4g} kcal/mol/A; {self.record["wall_time"]:.0f} s)')
        return terms, normalize_for_json(self.record)


def reapply_rule(record, theta):
    """
    The retention rule of ADR-008 item 7 applied to a stored selection record
    with another threshold factor (no QM: the evidence points are in the
    record). Returns (retained terms, the per-candidate decisions with their
    reasons, the thresholds with the flag `exact`).

    As in `CouplingSelection.run`, whatever the order in which the record was
    stored or read back: the pairs are decided first, then the triples; a
    triple is a candidate where `run` would make it one at this threshold
    (every triple up to `selection_max_groups_all_triples` groups, above that
    only the triples whose three pairs are retained on their own evidence;
    the others are labelled "not a candidate"); a retained triple pulls in
    its pairs (labelled "retained as a face of", each keeping its own
    reason); the result is the downward closure, the singletons taken from
    the record's node grids.

    Exact up to three groups. From four groups on, a triple's evidence was
    measured against a base model holding the pairs retained at the record's
    own threshold (ADR-008 item 1), which another threshold changes, so the
    triple decisions at another threshold are approximate (`exact` False).
    """

    tol_e = float(record['thresholds']['grouping_tolerance_energy_kcal'])
    tol_f = float(record['thresholds']['grouping_tolerance_force_kcal_per_angstrom'])
    tau_e, tau_f = float(theta) * tol_e, float(theta) * tol_f
    groups = sorted({int(g) for g in record.get('node_grids', {})} | {int(g) for key in record['candidates'] for g in key.split('_')})
    max_all = int(record.get('parameters', {}).get('selection_max_groups_all_triples', selection_parameters()['selection_max_groups_all_triples']))
    candidates = sorted(((canonical_group_key(tuple(int(g) for g in key.split('_'))), key, cand) for key, cand in record['candidates'].items()),
                        key=lambda item: (len(item[0]), item[0]))
    decisions = {}
    retained = [(g,) for g in groups]
    pairs_on_own_evidence = set()
    for term, key, cand in candidates:
        if len(term) == 3 and len(groups) > max_all and not all(
                canonical_group_key(p) in pairs_on_own_evidence for p in subsets_of(term, 2, 2)):
            why = f'a pair of it is dropped at this threshold ({len(groups)} groups, more than selection_max_groups_all_triples = {max_all})'
            decisions[key] = {'decision': 'not a candidate', 'reason': why}
            continue
        keep, reason = False, 'below threshold'
        for kind in ('boundary', 'leak', 'probes'):
            block = cand.get(kind)
            if not block or keep:
                continue
            for sub in (block.values() if kind != 'probes' else [block]):
                for point in sub['points']:
                    if 'energy_kcal' not in point:
                        continue
                    if point['energy_kcal'] >= tau_e:
                        keep, reason = True, f'{kind}: {point["label"]}: energy {point["energy_kcal"]:.4g} kcal/mol >= {tau_e:.4g}'
                        break
                    if point['force_kcal_per_angstrom'] >= tau_f:
                        keep, reason = True, f'{kind}: {point["label"]}: force {point["force_kcal_per_angstrom"]:.4g} kcal/mol/A >= {tau_f:.4g}'
                        break
                if keep:
                    break
        decisions[key] = {'decision': 'retained' if keep else 'dropped', 'reason': reason}
        if not keep:
            continue
        retained.append(term)
        if len(term) == 2:
            pairs_on_own_evidence.add(term)
        if len(term) == 3:
            for pair in subsets_of(term, 2, 2):
                pair = canonical_group_key(pair)
                if pair not in retained:
                    retained.append(pair)
                    own = decisions.get(term_label(pair), {})
                    decisions[term_label(pair)] = {'decision': f'retained as a face of {term}', 'reason': own.get('reason', '')}
    terms = [t for t in downward_closure(retained) if t]
    return terms, decisions, {'energy_kcal': tau_e, 'force_kcal_per_angstrom': tau_f, 'theta': float(theta), 'exact': len(groups) <= 3}


def select_couplings(coordinates, reference, data_source, node_grids, cut_policy='relaxed_rotor_v1', response_policy='quadratic_v1',
                     node_orders=None, support_factor=None, weight_scheme=None, energy_gradient=None, relaxation_zone=None,
                     cut_records=None, parameters=None, validation=None, ostream=None, probe_context=None):
    """
    The coupling selection of ADR-008 as one call: builds the singleton cuts
    (or reuses `cut_records`), computes the evidence and the probes with a
    CutBuilder on `data_source`, and returns (retained_terms, record,
    cut_records). The retained terms are downward closed; the final model is
    built from the same records with GroupedChartModel.build(...,
    cut_records=cut_records).
    """

    from .groupeddefinitions import DEFAULT_NODE_ORDERS
    policy = coordinates.policy
    orders = dict(DEFAULT_NODE_ORDERS if node_orders is None else node_orders)
    if response_policy == 'none':
        orders = {'I': orders['I']}
    builder = CutBuilder(coordinates, reference, data_source, cut_policy, response_policy, orders,
                         energy_gradient=energy_gradient, relaxation_zone=relaxation_zone, records=cut_records)
    selection = CouplingSelection(coordinates, reference, builder, node_grids,
                                  policy.inner_support_factor if support_factor is None else support_factor,
                                  policy.inner_weight_scheme if weight_scheme is None else weight_scheme,
                                  parameters=parameters, validation=validation, ostream=ostream, probe_context=probe_context)
    terms, record = selection.run()
    return terms, record, builder.records
