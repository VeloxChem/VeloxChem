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
Groups, coordinates, cut maps and pullbacks of the grouped MSI v2 model
family (SPEC S5, S6; ADR-002 and ADR-003 with their amendments of
2026-10-09; PROMPT-M1.md sections 3.1, 3.2, 3.4 and 3.6).

Contents
- InternalCoordinates: the redundant Z-matrix coordinates q, B and B2 of a
  geometry and the legacy chart transformation, through InterpolationDatapoint
  (plain Cartesian coordinates, inverse bond lengths; nothing is copied)
- legacy_topology: the Z-matrix and the MM-rotatable bonds of a molecule, as
  IMForceFieldGenerator.set_up_the_system builds them
- bonds_from_user_input: the one place where 1-based user bonds become 0-based
- propose_groups / select_groups: the general group rules (root, candidates,
  orientation, types, m, n, selection with the MM flag and the override)
- GroupedCoordinateModel: the persisted group definitions with rotor rows,
  offsets and sense; the validity checks; the phases with first and second
  derivatives; rigid cut maps with tangents and accelerations (nested term
  included); pullbacks; the response rows; the complete-coordinate
  transformation and the node-data extraction of ADR-007

Indices are 0-based throughout. Angles are in radians, lengths in bohr.
"""

from dataclasses import dataclass
from pathlib import Path
import glob
import itertools

import numpy as np

from .errorhandler import assert_msg_critical
from .interpolationdatapoint import InterpolationDatapoint
from .outputstream import OutputStream
from .groupeddefinitions import GroupDefinition, GroupedModelPolicyV2, z_matrix_as_int_tuples

BOND, ANGLE, DIHEDRAL, IMPROPER = 0, 1, 2, 3
TWO_PI = 2.0 * np.pi


def wrap_angle(a):
    """
    Reduces an angle (or array) to (-pi, pi].
    """

    return (np.asarray(a, dtype=np.float64) + np.pi) % TWO_PI - np.pi


def rigid_body_projector(x):
    """
    The orthogonal projector of R^{3N} that removes the three translations
    and the three infinitesimal rotations about the centroid of x.
    """

    x = np.asarray(x, dtype=np.float64).reshape(-1, 3)
    n = x.shape[0]
    centred = x - x.mean(axis=0)
    modes = []
    for k in range(3):
        v = np.zeros((n, 3))
        v[:, k] = 1.0
        modes.append(v.reshape(-1))
    for k in range(3):
        axis = np.zeros(3)
        axis[k] = 1.0
        modes.append(np.cross(axis, centred).reshape(-1))
    u, sv, _ = np.linalg.svd(np.array(modes).T, full_matrices=False)
    u = u[:, sv > 1.0e-10 * sv[0]]
    return np.eye(3 * n) - u @ u.T


def rodrigues_rotation(coordinates, attachment_atom, axis_atom, moving_atoms, angle):
    """
    Rotates moving_atoms rigidly by angle about the axis attachment_atom ->
    axis_atom of the given geometry; everything else stays. Returns a copy.
    """

    x = np.array(coordinates, dtype=np.float64)
    axis = x[axis_atom] - x[attachment_atom]
    axis = axis / np.linalg.norm(axis)
    k = np.array([[0.0, -axis[2], axis[1]], [axis[2], 0.0, -axis[0]], [-axis[1], axis[0], 0.0]])
    rot = np.eye(3) + np.sin(angle) * k + (1.0 - np.cos(angle)) * (k @ k)
    idx = list(moving_atoms)
    x[idx] = (rot @ (x[idx] - x[axis_atom]).T).T + x[axis_atom]
    return x


class InternalCoordinates:
    """
    The redundant internal coordinates of a Z-matrix (bonds as 1/r, angles,
    proper dihedrals, impropers, in the flattened order of
    InterpolationDatapoint.flatten_z_matrix), with B and B2 of a geometry and
    the legacy chart transformation. Plain Cartesian coordinates only.
    """

    BOND_SETTINGS = {'use_inverse_bond_length': True, 'use_eq_bond_length': False}

    def __init__(self, z_matrix):
        self.z_matrix = z_matrix_as_int_tuples(z_matrix)
        self.rows = [tuple(r) for section in ('bonds', 'angles', 'dihedrals', 'impropers')
                     for r in self.z_matrix[section]]
        n_b, n_a, n_d, n_i = (len(self.z_matrix[s]) for s in ('bonds', 'angles', 'dihedrals', 'impropers'))
        self.kinds = np.array([BOND] * n_b + [ANGLE] * n_a + [DIHEDRAL] * n_d + [IMPROPER] * n_i, dtype=int)
        self.n_rows = len(self.rows)
        self.n_bonds = n_b
        self.bond_end = n_b
        self.angle_end = n_b + n_a
        self.dihedral_end = n_b + n_a + n_d
        self.periodic_rows = np.flatnonzero(self.kinds >= DIHEDRAL)
        self._datapoint = self._new_datapoint()

    def _new_datapoint(self):
        # a null stream: a default OutputStream closes its stream on garbage
        # collection, which is not wanted for the many transient datapoints
        dp = InterpolationDatapoint({k: list(v) for k, v in self.z_matrix.items()}, ostream=OutputStream(None))
        dp.update_settings(dict(self.BOND_SETTINGS))
        dp.eq_bond_lengths = np.zeros(self.n_bonds)
        dp.define_internal_coordinates()
        return dp

    @property
    def n_atoms_from(self):
        return 1 + max(a for row in self.rows for a in row)

    def values(self, x):
        dp = self._datapoint
        dp.cartesian_coordinates = np.array(x, dtype=np.float64)
        dp.compute_internal_coordinates_values()
        return np.array(dp.internal_coordinates_values, dtype=np.float64)

    def values_and_b(self, x):
        dp = self._datapoint
        dp.reset_coordinates_impes_driver(np.array(x, dtype=np.float64))
        return (np.array(dp.internal_coordinates_values, dtype=np.float64), np.array(dp.b_matrix, dtype=np.float64))

    def values_b_b2(self, x):
        dp = self._datapoint
        dp.reset_coordinates(np.array(x, dtype=np.float64))
        return (np.array(dp.internal_coordinates_values, dtype=np.float64), np.array(dp.b_matrix, dtype=np.float64),
                np.array(dp.b2_matrix, dtype=np.float64))

    def chart_data(self, x, energy, gradient, hessian):
        """
        The legacy transformation of (E, g, H) at x to internal coordinates
        (InterpolationDatapoint.transform_gradient_and_hessian): the chart
        data (E, q, g_q, H_q) of a datapoint at x.
        """

        dp = self._new_datapoint()
        dp.cartesian_coordinates = np.array(x, dtype=np.float64)
        dp.energy = float(energy)
        dp.gradient = np.array(gradient, dtype=np.float64).reshape(-1, 3)
        dp.hessian = np.array(hessian, dtype=np.float64)
        dp.transform_gradient_and_hessian()
        return {
            'energy': float(energy),
            'q': np.array(dp.internal_coordinates_values, dtype=np.float64),
            'g_q': np.array(dp.internal_gradient, dtype=np.float64),
            'H_q': np.array(dp.internal_hessian, dtype=np.float64),
            'x': np.array(x, dtype=np.float64),
        }

    def displacement(self, d, rows=None):
        """
        The legacy displacement functions s(d) of compute_potential with
        their first and second derivatives, per row: d for bonds and angles,
        2 sin(d/2) for proper dihedrals and 2 tan(d/2) for impropers, with d
        wrapped to (-pi, pi]. rows selects a subset of rows (default all).
        """

        d = np.array(d, dtype=np.float64)
        kinds = self.kinds if rows is None else self.kinds[np.asarray(rows)]
        s = d.copy()
        s1 = np.ones_like(d)
        s2 = np.zeros_like(d)
        dih = kinds == DIHEDRAL
        if dih.any():
            dd = wrap_angle(d[dih])
            s[dih] = 2.0 * np.sin(0.5 * dd)
            s1[dih] = np.cos(0.5 * dd)
            s2[dih] = -0.5 * np.sin(0.5 * dd)
        imp = kinds == IMPROPER
        if imp.any():
            di = wrap_angle(d[imp])
            sec2 = 1.0 / np.cos(0.5 * di) ** 2
            s[imp] = 2.0 * np.tan(0.5 * di)
            s1[imp] = sec2
            s2[imp] = sec2 * np.tan(0.5 * di)
        return s, s1, s2

    def row_permutation(self, atom_permutation):
        """
        For an atom permutation p with (Px)_i = x_p[i]: the row permutation
        Pi with q(Px)[i] = q(x)[Pi[i]], or None if the Z-matrix is not closed
        under p (rows of proper dihedrals, angles and bonds may be reversed;
        impropers must match exactly).
        """

        p = np.asarray(atom_permutation, dtype=int)
        index = {}
        for i, row in enumerate(self.rows):
            index[(self.kinds[i], row)] = i
            if self.kinds[i] != IMPROPER:
                index.setdefault((self.kinds[i], tuple(reversed(row))), i)
        out = np.empty(self.n_rows, dtype=int)
        for i, row in enumerate(self.rows):
            image = tuple(int(p[a]) for a in row)
            j = index.get((self.kinds[i], image))
            if j is None:
                return None
            out[i] = j
        if len(set(out.tolist())) != self.n_rows:
            return None
        return out


def legacy_topology(molecule):
    """
    The Z-matrix dict (bonds, angles, dihedrals, impropers) and the 0-based
    MM-rotatable bonds of a molecule, built exactly as the legacy
    IMForceFieldGenerator.set_up_the_system builds them (geomeTRIC primitives
    plus the MM impropers; partial charges from the molecule, so no SCF
    runs). The generator reads an existing im_database_*.h5 from the working
    directory instead, so this refuses to run in such a directory.
    """

    existing = sorted(glob.glob('im_database_*.h5'))
    assert_msg_critical(
        not existing,
        'grouped MSI: legacy_topology must run in a directory without im_database_*.h5 '
        f'(found {existing} in {Path.cwd()}).')

    from .imforcefieldgenerator import IMForceFieldGenerator
    from .outputstream import OutputStream
    from .scfrestdriver import ScfRestrictedDriver

    qm_driver = ScfRestrictedDriver(ostream=OutputStream(None))
    generator = IMForceFieldGenerator(qm_driver, ostream=OutputStream(None))
    generator.add_conformal_structures = False
    generator.set_up_the_system(molecule)
    z_matrix = {k: [tuple(int(a) for a in row) for row in v] for k, v in generator.roots_z_matrix[0].items()}
    rotatable = sorted(tuple(sorted((int(i), int(j)))) for i, j in generator.all_rotatable_bonds)
    return z_matrix, rotatable


def bonds_from_user_input(bonds, n_atoms):
    """
    Converts user bond input (1-based atom index pairs, as in the legacy
    rotatable-bond settings) to sorted 0-based pairs. The only place where
    the conversion happens.
    """

    out = []
    for bond in bonds:
        assert_msg_critical(
            len(bond) == 2,
            f'grouped MSI: a bond must be a pair of 1-based atom indices, got {bond!r}.')
        i, j = int(bond[0]), int(bond[1])
        assert_msg_critical(
            1 <= i <= n_atoms and 1 <= j <= n_atoms and i != j,
            f'grouped MSI: bond {bond!r} is not a pair of distinct 1-based atom indices in 1..{n_atoms}.')
        out.append((min(i, j) - 1, max(i, j) - 1))
    return out


def _neighbours(bonds, n_atoms):
    nb = {i: set() for i in range(n_atoms)}
    for i, j in bonds:
        nb[int(i)].add(int(j))
        nb[int(j)].add(int(i))
    return nb


def _side(nb, start, blocked_bond):
    """
    Atoms reachable from start without crossing blocked_bond (i, j).
    """

    i, j = blocked_bond
    seen, stack = {start}, [start]
    while stack:
        a = stack.pop()
        for b in nb[a]:
            if {a, b} == {i, j} or b in seen:
                continue
            seen.add(b)
            stack.append(b)
    return seen


def _ring_atoms_and_bonds(nb, bonds):
    in_ring = {b: b[0] in _side(nb, b[1], b) for b in bonds}
    ring_bonds = [b for b in bonds if in_ring[b]]
    systems = []
    atoms = sorted({a for b in ring_bonds for a in b})
    for a in atoms:
        if any(a in s for s in systems):
            continue
        s, stack = {a}, [a]
        while stack:
            u = stack.pop()
            for v in nb[u]:
                if v not in s and in_ring.get(tuple(sorted((u, v))), False):
                    s.add(v)
                    stack.append(v)
        systems.append(s)
    return in_ring, systems


@dataclass(frozen=True)
class RotorCandidate:
    """
    A candidate bond of the general group rules with its orientation, type
    and defaults (ADR-003 amendment 1). Atom indices 0-based.
    """

    bond: tuple
    attachment_atom: int
    axis_atom: int
    moving_atoms: tuple
    group_type: str
    images: int
    harmonic: int
    mm_rotatable: bool
    oriented_towards_root_side: bool = False

    def describe(self, labels=None):
        kind = 'terminal' if self.group_type == 'terminal_rotor' else 'fragment'
        return (f'bond {self.bond}: {kind} rotor, attach {self.attachment_atom}, axis {self.axis_atom}, '
                f'moving {list(self.moving_atoms)}, m {self.images}, n {self.harmonic}, '
                f'MM-rotatable {self.mm_rotatable}')


def propose_groups(labels, bonds, rotatable_bonds):
    """
    The root and every rotor candidate of a molecule (ADR-003 amendment 1).

    labels: element symbols; bonds: 0-based bonded pairs; rotatable_bonds:
    0-based MM-rotatable pairs. Returns (root atoms as a sorted tuple,
    candidates sorted by bond, rejected candidates as a list of (bond,
    reason)).
    """

    n_atoms = len(labels)
    bonds = sorted({tuple(sorted((int(i), int(j)))) for i, j in bonds})
    rotatable = {tuple(sorted((int(i), int(j)))) for i, j in rotatable_bonds}
    nb = _neighbours(bonds, n_atoms)
    in_ring, systems = _ring_atoms_and_bonds(nb, bonds)
    heavy = [i for i in range(n_atoms) if labels[i] != 'H']
    assert_msg_critical(heavy, 'grouped MSI: the molecule has no heavy atom.')

    if systems:
        root = max(systems, key=lambda s: (len(s), -min(s)))
    else:
        def eccentricity(a):
            dist, front, d = {a: 0}, [a], 0
            while front:
                d += 1
                nxt = []
                for u in front:
                    for v in nb[u]:
                        if v in heavy and v not in dist:
                            dist[v] = d
                            nxt.append(v)
                front = nxt
            return max(dist.values())

        best = min(heavy, key=lambda a: (eccentricity(a), -sum(labels[v] != 'H' for v in nb[a]), a))
        root = {best}
    root = frozenset(root)

    def symmetric_terminal_side(side_atoms):
        # k >= 2 identical terminal atoms and nothing else
        if len(side_atoms) < 2:
            return False
        return all(len(nb[a]) == 1 for a in side_atoms) and len({labels[a] for a in side_atoms}) == 1

    candidates = []
    rejected = []
    for bond in bonds:
        i, j = bond
        if in_ring[bond] or len(nb[i]) < 2 or len(nb[j]) < 2:
            continue
        side_i = _side(nb, i, bond)
        side_j = _side(nb, j, bond)
        if root & side_j and not root & side_i:
            attach, axis, moving_side, root_side = j, i, side_i, side_j
        elif root & side_i and not root & side_j:
            attach, axis, moving_side, root_side = i, j, side_j, side_i
        else:
            assert_msg_critical(False, f'grouped MSI: the root is split by the bond {bond}; this cannot happen.')
        moving = sorted(moving_side - {axis})
        root_moving = sorted(root_side - {attach})
        oriented_towards_root = False
        if symmetric_terminal_side(root_moving) and not symmetric_terminal_side(moving):
            # methanol with the root on the carbon: the symmetric group must be the moving one
            attach, axis = axis, attach
            moving = root_moving
            oriented_towards_root = True
        elif symmetric_terminal_side(root_moving) and symmetric_terminal_side(moving):
            rejected.append((bond, 'symmetric terminal groups at both ends of the bond (ethane type): the '
                                   'permutation of the root-side group would not be represented'))
            continue
        terminal = all(len(nb[a]) == 1 for a in moving)
        if terminal and len(moving) >= 2 and len({labels[a] for a in moving}) == 1:
            m = len(moving)
        else:
            m = 1
        n = m if m > 1 else len(nb[attach]) - 1
        candidates.append(RotorCandidate(
            bond=bond, attachment_atom=attach, axis_atom=axis, moving_atoms=tuple(moving),
            group_type='terminal_rotor' if terminal else 'fragment_rotor', images=m, harmonic=n,
            mm_rotatable=bond in rotatable, oriented_towards_root_side=oriented_towards_root))
    return tuple(sorted(root)), candidates, rejected


def select_groups(candidates, fragment_bonds=(), excluded_bonds=(), override_bonds=(), harmonics=None,
                  ostream=None):
    """
    The default selection (every terminal rotor) with the user's changes:
    fragment_bonds adds fragment rotors, excluded_bonds removes terminal
    rotors, override_bonds allows candidates that MM does not consider
    rotatable (with a warning), harmonics overrides n per bond. Bonds are
    0-based sorted pairs (use bonds_from_user_input for user input).
    """

    by_bond = {c.bond: c for c in candidates}
    fragment_bonds = {tuple(sorted(b)) for b in fragment_bonds}
    excluded_bonds = {tuple(sorted(b)) for b in excluded_bonds}
    override_bonds = {tuple(sorted(b)) for b in override_bonds}
    harmonics = {tuple(sorted(b)): int(n) for b, n in (harmonics or {}).items()}
    for name, bonds in (('fragment', fragment_bonds), ('excluded', excluded_bonds), ('override', override_bonds),
                        ('harmonic', set(harmonics))):
        unknown = sorted(bonds - set(by_bond))
        assert_msg_critical(not unknown, f'grouped MSI: {name} bonds {unknown} are not rotor candidates.')
    for bond in fragment_bonds:
        assert_msg_critical(
            by_bond[bond].group_type == 'fragment_rotor',
            f'grouped MSI: bond {bond} is a terminal rotor and is selected by default; it cannot be added as a fragment rotor.')
    for bond in excluded_bonds:
        assert_msg_critical(
            by_bond[bond].group_type == 'terminal_rotor',
            f'grouped MSI: bond {bond} is a fragment rotor; it is not selected unless named, so it cannot be excluded.')

    selected = []
    for cand in candidates:
        wanted = (cand.group_type == 'terminal_rotor' and cand.bond not in excluded_bonds) or cand.bond in fragment_bonds
        if not wanted:
            continue
        if not cand.mm_rotatable:
            assert_msg_critical(
                cand.bond in override_bonds,
                f'grouped MSI: bond {cand.bond} is not rotatable according to MMForceFieldGenerator; '
                'selecting it needs an explicit override.')
            if ostream is not None:
                ostream.print_warning(f'grouped MSI: bond {cand.bond} is selected as a rotor although the MM force '
                                      'field does not consider it rotatable (user override).')
                ostream.flush()
        if cand.bond in harmonics:
            assert_msg_critical(cand.images == 1, f'grouped MSI: the harmonic of bond {cand.bond} is fixed by its images.')
            cand = RotorCandidate(**{**cand.__dict__, 'harmonic': harmonics[cand.bond]})
        selected.append(cand)
    return selected


class GroupedCoordinateModel:
    """
    The groups of a model on one anchor geometry: persisted definitions, the
    validity checks, the phases, the rigid cut maps with tangents and
    accelerations, the pullbacks, the response rows and the node-data split.

    :param molecule: the Molecule (labels and connectivity)
    :param z_matrix: the Z-matrix dict (legacy_topology)
    :param anchor_coordinates: the anchor geometry in bohr, (N, 3)
    :param candidates: the selected RotorCandidate list
    :param policy: a GroupedModelPolicyV2 (phase validity threshold)
    :param harmonics / domains / relaxation_zones: optional per-bond overrides
    """

    def __init__(self, molecule, z_matrix, anchor_coordinates, candidates, policy=None, domains=None,
                 relaxation_zones=None, ostream=None, check=True):
        self.labels = list(molecule.get_labels())
        self.n_atoms = len(self.labels)
        self.policy = policy if policy is not None else GroupedModelPolicyV2()
        self.ostream = ostream
        self.ic = InternalCoordinates(z_matrix)
        self.x_anchor = np.array(anchor_coordinates, dtype=np.float64).reshape(self.n_atoms, 3)
        conn = np.asarray(molecule.get_connectivity_matrix())
        conn_bonds = sorted({tuple(sorted((int(i), int(j)))) for i, j in zip(*np.nonzero(conn)) if i < j})
        z_bonds = sorted(tuple(sorted(b)) for b in self.ic.z_matrix['bonds'])
        assert_msg_critical(
            conn_bonds == z_bonds,
            'grouped MSI: the bonds of the Z-matrix differ from the connectivity of the molecule.')
        self.bonds = z_bonds
        self.neighbours = _neighbours(self.bonds, self.n_atoms)
        self.in_ring, _ = _ring_atoms_and_bonds(self.neighbours, self.bonds)
        domains = {tuple(sorted(b)): d for b, d in (domains or {}).items()}
        relaxation_zones = {tuple(sorted(b)): z for b, z in (relaxation_zones or {}).items()}

        self.q_anchor = self.ic.values(self.x_anchor)
        self.groups = []
        self._cyclic_perm = []
        self._row_perm = []
        for gid, cand in enumerate(candidates):
            assert_msg_critical(
                not self.in_ring.get(cand.bond, False),
                f'grouped MSI: bond {cand.bond} is a ring bond and cannot be a rotor.')
            group, perm, row_perm = self._define_group(gid, cand, domains.get(cand.bond), relaxation_zones.get(cand.bond))
            self.groups.append(group)
            self._cyclic_perm.append(perm)
            self._row_perm.append(row_perm)
        self.n_groups = len(self.groups)
        self.rotor_rows = sorted(r for g in self.groups for r in g.rotor_rows)
        self.response_rows = np.array([i for i in range(self.ic.n_rows) if i not in set(self.rotor_rows)], dtype=int)
        moving = set(a for g in self.groups for a in g.moving_atoms)
        self.core_atoms = tuple(i for i in range(self.n_atoms) if i not in moving)
        self.phi_anchor = self.phases(self.q_anchor)[0]
        if check:
            self.check_validity()
            self.assert_angles_away_from_linear(self.q_anchor, label='the anchor')

    # ------------------------------------------------------------ definitions

    def _rotor_rows_of(self, attach, axis):
        return [i for i, row in enumerate(self.ic.rows)
                if self.ic.kinds[i] == DIHEDRAL and {row[1], row[2]} == {attach, axis}]

    def _define_group(self, gid, cand, domain, relaxation_zone):
        attach, axis, moving, m = cand.attachment_atom, cand.axis_atom, list(cand.moving_atoms), cand.images
        rows = self._rotor_rows_of(attach, axis)
        assert_msg_critical(rows, f'grouped MSI: no proper dihedral about the bond {cand.bond} in the Z-matrix.')
        q_rows = self.q_anchor[rows]
        perm = np.arange(self.n_atoms)
        row_perm = None
        if m == 1:
            offsets = q_rows.copy()
            sense = 1
            cyclic = ()
        else:
            cyclic = tuple(moving)
            far_end, position = [], []
            for i in rows:
                row = self.ic.rows[i]
                ends = (row[0], row[3])
                inside = [a for a in ends if a in cyclic]
                assert_msg_critical(
                    len(inside) == 1,
                    f'grouped MSI: rotor row {row} of bond {cand.bond} must contain exactly one moving atom.')
                position.append(cyclic.index(inside[0]))
                far_end.append(ends[1] if ends[0] == inside[0] else ends[0])
            far_end, position = np.array(far_end), np.array(position)
            signs = []
            for x_atom in sorted(set(far_end.tolist())):
                assert_msg_critical(
                    sorted(position[far_end == x_atom].tolist()) == list(range(m)),
                    f'grouped MSI: the rotor rows of bond {cand.bond} with far-end atom {x_atom} are incomplete.')
                if m > 2:
                    val = {j: q_rows[(far_end == x_atom) & (position == j)][0] for j in range(m)}
                    for j in range(m):
                        signs.append(int(np.sign(wrap_angle(val[(j + 1) % m] - val[j]))))
            if m > 2:
                assert_msg_critical(
                    len(set(signs)) == 1 and signs[0] != 0,
                    f'grouped MSI: the cyclic order of the moving atoms of bond {cand.bond} has no consistent sense.')
                sense = signs[0]
            else:
                sense = 1
            offsets = np.zeros(len(rows))
            for x_atom in set(far_end.tolist()):
                sel = far_end == x_atom
                ideal = sense * TWO_PI * position[sel] / m
                c_x = np.angle(np.exp(1j * (q_rows[sel] - ideal)).sum())
                offsets[sel] = c_x + ideal
            perm[list(cyclic)] = np.roll(cyclic, -1)
            row_perm = self.ic.row_permutation(perm)
            assert_msg_critical(
                row_perm is not None,
                f'grouped MSI: the Z-matrix is not closed under the cyclic permutation of the moving atoms of bond {cand.bond}.')
        group = GroupDefinition(
            group_id=gid, group_type=cand.group_type, attachment_atom=attach, axis_atom=axis,
            moving_atoms=tuple(sorted(moving)), images=m, harmonic=cand.harmonic, sense=sense,
            cyclic_order=cyclic, rotor_rows=tuple(int(r) for r in rows), offsets=tuple(float(o) for o in offsets),
            domain=None if domain is None else (float(domain[0]), float(domain[1])),
            relaxation_zone=None if relaxation_zone is None else tuple(int(r) for r in relaxation_zone),
            mm_rotatable=cand.mm_rotatable, user_override=not cand.mm_rotatable,
            label=f'{self.labels[attach]}{attach}-{self.labels[axis]}{axis}')
        return group, perm, row_perm

    def cyclic_permutation(self, gid, power=1):
        """
        The atom permutation p with (Px)_i = x_p[i] of the cyclic generator of
        group gid raised to `power` (identity for a group without images).
        """

        p = np.arange(self.n_atoms)
        gen = self._cyclic_perm[gid]
        for _ in range(power % max(self.groups[gid].images, 1)):
            p = p[gen]
        return p

    def row_permutation(self, gid, power=1):
        """
        The row permutation Pi of the generator of group gid raised to power:
        q(Px) = q(x)[Pi].
        """

        pi = np.arange(self.ic.n_rows)
        gen = self._row_perm[gid]
        if gen is None:
            return pi
        for _ in range(power % self.groups[gid].images):
            pi = pi[gen]
        return pi

    def image_elements(self):
        """
        Every element of the product of the image groups as (exponents, atom
        permutation); the identity first.
        """

        out = []
        ranges = [range(g.images) for g in self.groups]
        for expo in itertools.product(*ranges):
            p = np.arange(self.n_atoms)
            for gid, e in enumerate(expo):
                for _ in range(e):
                    p = p[self._cyclic_perm[gid]]
            out.append((tuple(expo), p))
        return out

    def describe(self):
        lines = [f'groups: {self.n_groups}, core atoms: {len(self.core_atoms)}']
        for g in self.groups:
            lines.append(f'  group {g.group_id} ({g.label}): {g.group_type}, bond {g.rotor_bond}, moving {list(g.moving_atoms)}, '
                         f'm {g.images}, n {g.harmonic}, sense {g.sense}, rows {list(g.rotor_rows)}')
        return '\n'.join(lines)

    # ----------------------------------------------------------------- phases

    def phases(self, q, order=2):
        """
        The symmetrised rotor angles of all groups from the internal
        coordinates q: phi (G,), dphi/dq (G, n_rows), d2phi/dq2 (G, n_rows,
        n_rows), and the validity ratios |Z_g|/|K_g|. Derivatives above
        `order` are None.
        """

        q = np.asarray(q, dtype=np.float64)
        phi = np.zeros(self.n_groups)
        ratio = np.zeros(self.n_groups)
        dphi = np.zeros((self.n_groups, self.ic.n_rows)) if order >= 1 else None
        d2phi = np.zeros((self.n_groups, self.ic.n_rows, self.ic.n_rows)) if order >= 2 else None
        for gid, g in enumerate(self.groups):
            rows = list(g.rotor_rows)
            a = q[rows] - np.array(g.offsets)
            z = np.exp(1j * a).sum()
            phi[gid] = np.angle(z)
            ratio[gid] = abs(z) / len(rows)
            c = np.cos(a - phi[gid])
            s = np.sin(a - phi[gid])
            big_c = c.sum()
            if order >= 1:
                dphi[gid, rows] = c / big_c
            if order >= 2:
                block = (np.outer(s, c) + np.outer(c, s) - np.diag(s) * big_c) / big_c ** 2
                d2phi[gid][np.ix_(rows, rows)] = block
        return phi, dphi, d2phi, ratio

    def phases_cartesian(self, x, order=2, qbb2=None):
        """
        phi (G,), grad phi (G, 3N) and the Cartesian Hessians (G, 3N, 3N) of
        the phases at x, through B and B2. qbb2 may pass precomputed (q, B,
        B2).
        """

        if qbb2 is None:
            q, b_mat, b2 = self.ic.values_b_b2(x) if order >= 2 else (*self.ic.values_and_b(x), None)
        else:
            q, b_mat, b2 = qbb2
        phi, dphi, d2phi, ratio = self.phases(q, order)
        grad = dphi @ b_mat if order >= 1 else None
        hess = None
        if order >= 2:
            hess = np.einsum('gkl,kx,ly->gxy', d2phi, b_mat, b_mat) + np.einsum('gk,kxy->gxy', dphi, b2)
        return phi, grad, hess, ratio

    def phase_shift(self, x):
        """
        wrap(phi(x) - phi(anchor)) for all groups.
        """

        return wrap_angle(self.phases(self.ic.values(x), order=0)[0] - self.phi_anchor)

    # --------------------------------------------------------------- rigid cuts

    def rigid_cut(self, angles, start=None):
        """
        The rigid cut geometry: each group rotated by its angle about its
        rotor bond, axes taken from the current geometry (nested groups
        commute). angles: a sequence over groups, or a dict gid -> angle.
        """

        x = np.array(self.x_anchor if start is None else start, dtype=np.float64)
        if isinstance(angles, dict):
            items = sorted(angles.items())
        else:
            items = list(enumerate(angles))
        for gid, angle in items:
            if angle != 0.0:
                g = self.groups[gid]
                x = rodrigues_rotation(x, g.attachment_atom, g.axis_atom, g.moving_atoms, float(angle))
        return x

    def tangent(self, x, gid):
        """
        v_g = n_g x (x_i - x_c) on the moving atoms of group gid, zero
        elsewhere: the derivative of the rigid rotation about the current
        axis, (N, 3).
        """

        g = self.groups[gid]
        axis = x[g.axis_atom] - x[g.attachment_atom]
        axis = axis / np.linalg.norm(axis)
        v = np.zeros_like(x)
        idx = list(g.moving_atoms)
        v[idx] = np.cross(axis, x[idx] - x[g.axis_atom])
        return v

    def acceleration(self, x, gid, hid):
        """
        d2 pi / dphi_g dphi_h of the rigid cut map at x: n_g x v_g for g = h,
        n_outer x v_inner on the inner moving atoms for nested groups, zero
        for disjoint ones (ADR-003 amendment 2).
        """

        g, h = self.groups[gid], self.groups[hid]
        set_g, set_h = set(g.moving_atoms), set(h.moving_atoms)
        if gid == hid or set_h < set_g:
            outer, inner = gid, hid
        elif set_g < set_h:
            outer, inner = hid, gid
        else:
            assert_msg_critical(not (set_g & set_h), 'grouped MSI: moving sets are neither nested nor disjoint.')
            return np.zeros_like(x)
        go = self.groups[outer]
        axis = x[go.axis_atom] - x[go.attachment_atom]
        axis = axis / np.linalg.norm(axis)
        return np.cross(axis, self.tangent(x, inner))

    def pullback(self, x, gradient, hessian=None):
        """
        The derivatives of a Cartesian scalar (gradient (N,3) and optionally
        Hessian (3N,3N) of R) along the rigid cut: dR/dphi_g = v_g . g_R and
        d2R/dphi_g dphi_h = v_g^T H_R v_h + g_R . d2pi/dphi_g dphi_h.
        """

        g_flat = np.asarray(gradient, dtype=np.float64).reshape(-1)
        v = [self.tangent(x, gid).reshape(-1) for gid in range(self.n_groups)]
        first = np.array([vg @ g_flat for vg in v])
        if hessian is None:
            return first, None
        h_mat = np.asarray(hessian, dtype=np.float64)
        second = np.zeros((self.n_groups, self.n_groups))
        for gid in range(self.n_groups):
            for hid in range(self.n_groups):
                second[gid, hid] = v[gid] @ h_mat @ v[hid] + g_flat @ self.acceleration(x, gid, hid).reshape(-1)
        return first, second

    # ------------------------------------------------------------------ checks

    def check_validity(self, angle=37.0 * np.pi / 180.0, tol=1.0e-9):
        """
        The numerical validity checks of ADR-003 amendment 1 at the anchor;
        every failure is a critical error.
        """

        rows = [set(g.rotor_rows) for g in self.groups]
        for a, b in itertools.combinations(range(self.n_groups), 2):
            assert_msg_critical(
                not (rows[a] & rows[b]),
                f'grouped MSI: groups {a} and {b} share rotor rows {sorted(rows[a] & rows[b])}.')
            set_a, set_b = set(self.groups[a].moving_atoms), set(self.groups[b].moving_atoms)
            assert_msg_critical(
                (set_a & set_b) in (set(), set_a, set_b),
                f'grouped MSI: the moving sets of groups {a} and {b} overlap without nesting (not laminar).')
        for gid, g in enumerate(self.groups):
            y = self.rigid_cut({gid: angle})
            dq = wrap_angle(self.ic.values(y) - self.q_anchor)
            changed = set(np.flatnonzero(np.abs(dq) > tol).tolist())
            assert_msg_critical(
                changed == rows[gid],
                f'grouped MSI: a rigid rotation of group {gid} changes the rows {sorted(changed)} instead of '
                f'its rotor rows {sorted(rows[gid])}.')
            assert_msg_critical(
                np.allclose(np.abs(dq[list(rows[gid])]), angle, atol=tol),
                f'grouped MSI: the rotor rows of group {gid} do not all change by the rotation angle.')
            shift = self.phase_shift(y)
            others = [h for h in range(self.n_groups) if h != gid]
            assert_msg_critical(
                np.all(np.abs(shift[others]) < tol) if others else True,
                f'grouped MSI: rotating group {gid} moves the phase of another group by up to '
                f'{np.abs(shift[others]).max() if others else 0.0:.2e} rad.')
            assert_msg_critical(
                abs(wrap_angle(shift[gid] - angle)) < tol,
                f'grouped MSI: the phase of group {gid} changes by {shift[gid]:.9f} instead of {angle:.9f} rad.')
        for a, b in itertools.combinations(range(self.n_groups), 2):
            y1 = self.rigid_cut({b: -1.1}, start=self.rigid_cut({a: 0.7}))
            y2 = self.rigid_cut({a: 0.7}, start=self.rigid_cut({b: -1.1}))
            assert_msg_critical(
                np.abs(y1 - y2).max() < 1.0e-10,
                f'grouped MSI: the cut maps of groups {a} and {b} do not commute.')
        self.check_rank(self.x_anchor)

    def describe_row(self, row, q=None):
        """
        'angle C12-O15-C16' with its value (degrees for angular rows, bohr
        for bonds) when q is given.
        """

        kind = int(self.ic.kinds[row])
        name = {BOND: 'bond', ANGLE: 'angle', DIHEDRAL: 'dihedral', IMPROPER: 'improper'}[kind]
        atoms = '-'.join(f'{self.labels[a]}{a}' for a in self.ic.rows[row])
        if q is None:
            return f'{name} {atoms}'
        value = f'{1.0 / q[row]:.3f} bohr' if kind == BOND else f'{np.degrees(q[row]):.2f} deg'
        return f'{name} {atoms} at {value}'

    def linear_angle_rows(self, q):
        """
        The bond-angle rows within `linear_angle_margin_degrees` of 180
        degrees at the internal coordinates q (decision 5 of M1-C), as
        (row, degrees).
        """

        margin = self.policy.linear_angle_margin_degrees
        out = []
        for row in range(self.ic.n_rows):
            if self.ic.kinds[row] == ANGLE and 180.0 - np.degrees(q[row]) < margin:
                out.append((row, float(np.degrees(q[row]))))
        return out

    def assert_angles_away_from_linear(self, q, label='the geometry'):
        """
        Critical error if a bond angle lies within the persisted margin of
        180 degrees (ADR-005 coverage criterion, applied at construction
        from M1-C on). q may be a geometry (N, 3) or internal coordinates.
        """

        q = np.asarray(q, dtype=np.float64)
        if q.ndim == 2:
            q = self.ic.values(q)
        rows = self.linear_angle_rows(q)
        assert_msg_critical(
            not rows,
            f'grouped MSI: {label} has a bond angle within {self.policy.linear_angle_margin_degrees:g} degrees of 180: '
            + ', '.join(self.describe_row(r, q) for r, _ in rows)
            + '. Near-linear angles make the B rows of the dihedrals through them singular (1/sin); '
            'linear-bend coordinates are not supported in schema 2.0.0: choose another anchor or exclude the geometry.')
        return rows

    def check_rank(self, x, tol=1.0e-6, label='the anchor'):
        """
        rank(B_v) = 3N - 6 - G, rank([B_v; grad phi]) = 3N - 6, B_v v_g = 0
        and the rank of the whole B matrix (ADR-002), at x.
        """

        q, b_mat, b2 = self.ic.values_b_b2(x)
        expected = 3 * self.n_atoms - 6
        sv_all = np.linalg.svd(b_mat, compute_uv=False)
        assert_msg_critical(
            int((sv_all > tol * sv_all[0]).sum()) == expected,
            f'grouped MSI: at {label} the Z-matrix has rank {int((sv_all > tol * sv_all[0]).sum())} instead of 3N - 6 = {expected}.')
        b_v = b_mat[self.response_rows]
        sv = np.linalg.svd(b_v, compute_uv=False)
        rank_v = int((sv > tol * sv[0]).sum())
        assert_msg_critical(
            rank_v == expected - self.n_groups,
            f'grouped MSI: at {label} the non-rotor rows have rank {rank_v} instead of 3N - 6 - G = {expected - self.n_groups}.')
        _, grad, _, _ = self.phases_cartesian(x, order=1, qbb2=(q, b_mat, b2))
        sv2 = np.linalg.svd(np.vstack([b_v, grad]), compute_uv=False)
        rank_full = int((sv2 > tol * sv2[0]).sum())
        assert_msg_critical(
            rank_full == expected,
            f'grouped MSI: at {label} the non-rotor rows with the phase gradients have rank {rank_full} instead of {expected}.')
        for gid in range(self.n_groups):
            residual = np.abs(b_v @ self.tangent(x, gid).reshape(-1)).max()
            assert_msg_critical(
                residual < 1.0e-9,
                f'grouped MSI: at {label} B_v v_g is {residual:.2e} for group {gid}; the non-rotor rows are not invariant under its rotation.')

    # ------------------------------------------------------- node-data split

    def complete_transform(self, x, gradient, hessian=None, qbb2=None, label='the node', check=True):
        """
        The transformation of a rigid-body-free Cartesian gradient (N,3) and
        Hessian (3N,3N) to the complete coordinates y = (phi_1..phi_G, q_v)
        (ADR-007, node data step 2, amended by M1-C decision 3): the
        pseudo-inverse of G = J J^T keeps exactly the r = 3N - 6 largest
        singular values (the rank of a valid model is an invariant), and the
        conditioning is asserted: sigma_r / sigma_{r+1} >= split_rank_gap_min
        and sigma_1 / sigma_r <= split_condition_max (check=False only
        records them). gamma = G^- J g, Gamma = G^- J (H - sum_k gamma_k
        grad2 y_k) J^T G^-. Returns a dict with alpha (G,), beta (n_v,), the
        blocks H_phiphi, H_phiv, H_vv (None without a Hessian), J, and
        'diagnostics' (rank, sigma_r, sigma_r1, gap, condition, rows: the
        three internal rows with the largest weight in the singular vector
        of sigma_r, with their values).
        """

        if qbb2 is None:
            qbb2 = self.ic.values_b_b2(x) if hessian is not None else (*self.ic.values_and_b(x), None)
        q, b_mat, b2 = qbb2
        order = 2 if hessian is not None else 1
        phi, grad_phi, hess_phi, _ = self.phases_cartesian(x, order=order, qbb2=qbb2)
        b_v = b_mat[self.response_rows]
        j_mat = np.vstack([grad_phi, b_v])
        g_mat = j_mat @ j_mat.T
        u, s, vt = np.linalg.svd(g_mat)
        rank = 3 * self.n_atoms - 6
        s_inv = np.zeros_like(s)
        s_inv[:rank] = 1.0 / s[:rank]
        g_minus = (u * s_inv) @ vt
        gamma = g_minus @ (j_mat @ np.asarray(gradient, dtype=np.float64).reshape(-1))
        n_g = self.n_groups
        sigma_r1 = float(s[rank]) if len(s) > rank else 0.0
        gap = float(s[rank - 1] / sigma_r1) if sigma_r1 > 0.0 else np.inf
        condition = float(s[0] / s[rank - 1])
        vector = u[:, rank - 1]
        top = np.argsort(np.abs(vector))[::-1][:3]
        rows = []
        for i in top:
            if i < n_g:
                rows.append((f'phase of group {i} at {np.degrees(phi[i]):.2f} deg', float(vector[i])))
            else:
                rows.append((self.describe_row(int(self.response_rows[i - n_g]), q), float(vector[i])))
        diagnostics = {'rank': rank, 'sigma_r': float(s[rank - 1]), 'sigma_r1': sigma_r1, 'gap': gap,
                       'condition': condition, 'rows': rows}
        if check:
            policy = self.policy
            rows_txt = ', '.join(f'{name} ({weight:+.2f})' for name, weight in rows)
            assert_msg_critical(
                gap >= policy.split_rank_gap_min,
                f'grouped MSI: the complete-coordinate split at {label} has no clear rank: sigma_r = {s[rank - 1]:.2e}, '
                f'sigma_r+1 = {sigma_r1:.2e} (gap {gap:.1e} below split_rank_gap_min = {policy.split_rank_gap_min:.0e}); '
                f'rows of the smallest kept singular vector: {rows_txt}. The geometry is nearly singular for the Z-matrix: '
                'exclude the node (a smaller domain or a different grid) or choose another Z-matrix.')
            assert_msg_critical(
                condition <= policy.split_condition_max,
                f'grouped MSI: the complete-coordinate split at {label} is ill-conditioned: sigma_1 / sigma_r = {condition:.1e} '
                f'above split_condition_max = {policy.split_condition_max:.0e} (sigma_r = {s[rank - 1]:.2e}); rows of the smallest '
                f'kept singular vector: {rows_txt}. Noise in the node gradient and Hessian would be amplified into the forces: '
                'exclude the node or move the grid away from this geometry.')
        out = {'phi': phi, 'alpha': gamma[:n_g], 'beta': gamma[n_g:], 'J': j_mat, 'gamma': gamma,
               'H_phiphi': None, 'H_phiv': None, 'H_vv': None, 'diagnostics': diagnostics}
        if hessian is not None:
            h_mat = np.asarray(hessian, dtype=np.float64)
            curv = np.einsum('g,gxy->xy', gamma[:n_g], hess_phi) + np.einsum('k,kxy->xy', gamma[n_g:], b2[self.response_rows])
            big = g_minus @ j_mat @ (h_mat - curv) @ j_mat.T @ g_minus.T
            big = 0.5 * (big + big.T)
            out['Gamma'] = big
            out['H_phiphi'] = big[:n_g, :n_g]
            out['H_phiv'] = big[:n_g, n_g:]
            out['H_vv'] = big[n_g:, n_g:]
            out['grad2_y'] = (hess_phi, b2[self.response_rows])
        return out

    def valley_displacement(self, x):
        """
        w = wrap(q_v(x) - q_v^a): the anchor-relative non-rotor rows of a
        node geometry (ADR-007).
        """

        d = self.ic.values(x)[self.response_rows] - self.q_anchor[self.response_rows]
        periodic = self.ic.kinds[self.response_rows] >= DIHEDRAL
        d[periodic] = wrap_angle(d[periodic])
        return d

    def extract_node_data(self, x, residual, gradient, hessian, w_prime=None, qbb2=None, label='the node'):
        """
        The node data of ADR-007 at a node geometry x from the residual R, its
        rigid-body-free Cartesian gradient and its Hessian:

            C = H_vv, b = beta, b' = H_phiv + C w', I = R, I' = alpha + w'^T beta,
            I'' = H_phiphi + w'^T H_phiv + H_phiv^T w' + w'^T C w', w = valley displacement.

        w_prime (G, n_v) is the valley tangent; None means zero (rigid cuts, or
        a valley bank of order 0).
        """

        split = self.complete_transform(x, gradient, hessian, qbb2=qbb2, label=label)
        n_g, n_v = self.n_groups, len(self.response_rows)
        w1 = np.zeros((n_g, n_v)) if w_prime is None else np.asarray(w_prime, dtype=np.float64).reshape(n_g, n_v)
        beta, alpha = split['beta'], split['alpha']
        c_mat, h_pv, h_pp = split['H_vv'], split['H_phiv'], split['H_phiphi']
        data = {
            'I': float(residual),
            'I1': alpha + w1 @ beta,
            'I2': h_pp + w1 @ h_pv.T + h_pv @ w1.T + w1 @ c_mat @ w1.T,
            'b': beta.copy(),
            'b1': h_pv + w1 @ c_mat,
            'C': c_mat.copy(),
            'w': self.valley_displacement(x),
            'w1': w1.copy(),
            'phi': split['phi'].copy(),
            'diagnostics': split['diagnostics'],
        }
        return data
