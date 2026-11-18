"""Carbon structures: periodic graphite slabs, finite flakes and nanotubes (graphite only;
solute placement, solvation and the full system are built on top, in veloxchem.electrode).

All builders return a CarbonModel (nm). Slabs are Bernal/AA/ABC stacked,
periodic in x and y, and centred on z = 0: layer 0 (the top face) at
z = +t/2, the last layer at z = -t/2, t = (layers - 1) * INTERLAYER. In the
assembled system the periodic cell boundary runs through the slab, so water
and ions never cross it and z-based energies stay continuous.

Geometry is the ideal lattice (C-C 0.142 nm). There is no separate geometry
optimisation: frozen carbons keep it exactly, and the bonded parameters of the
flexible mode have their minimum at this geometry (see forces.add_carbon_mobility).
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np
from scipy.spatial import cKDTree

from .params import (CARBON_MASS, CARBON_SETS, CC_BOND, GAFF2_HA, GRAPHITE_RESNAME,
                     INTERLAYER, check_status)

MOBILITY_MODES = ('frozen', 'restrained', 'flexible')
CH_BOND = 0.1086          # nm, GAFF2 ca-ha
H_MASS = 1.008
ELEMENT_MASS = {'C': CARBON_MASS, 'H': H_MASS, 'O': 15.999}


@dataclass
class CarbonModel:
    kind: str                      # 'slab' | 'flake' | 'nanotube'
    carbon_set: str
    elements: list                 # 'C', 'H' (edge or group H) or 'O' (group) per atom
    xyz: np.ndarray                # nm
    layer: np.ndarray              # layer index per atom (0 = top); H/O: layer of their anchor
    periodic: tuple                # (x, y, z) periodicity used for bonds
    cell: tuple                    # periodic lengths (nm); None where not periodic
    mobility: str = 'frozen'
    restraint_k: float | None = None
    edge_charges: tuple | None = None     # (q on edge C, q on H) for flakes
    info: dict = field(default_factory=dict)
    warning: str | None = None
    # surface groups (surface.py): per-atom charge, LJ overrides (sigma, epsilon; NaN =
    # default for the element), group id (-1 = carbon skeleton), fragment id (-1 = none),
    # group bonded terms and the bonds that tie group atoms into the bond graph (model indices)
    charge: np.ndarray | None = None
    lj: np.ndarray | None = None
    group: np.ndarray | None = None
    fragment: np.ndarray | None = None    # group id whose template set this skeleton atom's charge
    terms: dict | None = None
    group_bonds: list | None = None
    features: list = field(default_factory=list)
    _bonds: np.ndarray | None = field(default=None, repr=False)

    # ---- generic
    @property
    def n_atoms(self):
        return len(self.elements)

    @property
    def n_carbon(self):
        return sum(1 for e in self.elements if e == 'C')

    @property
    def mass(self):
        return float(sum(ELEMENT_MASS[e] for e in self.elements))

    @property
    def skeleton(self):
        """True for carbon-skeleton atoms (lattice C and their edge H), False for group atoms."""
        if self.group is None:
            return np.ones(self.n_atoms, bool)
        return np.asarray(self.group) < 0

    def bonds(self, cutoff=None):
        """Covalent bonds (index pairs): neighbour search over the skeleton (periodic where
        set) plus the group bonds."""
        if self._bonds is None:
            skel = np.flatnonzero(self.skeleton)
            b = _bonds([self.elements[i] for i in skel], self.xyz[skel], self.periodic,
                       self.cell, cutoff)
            b = skel[b] if len(b) else b
            if self.group_bonds:
                b = np.vstack([b, np.asarray(self.group_bonds, int)])
            self._bonds = b
        return self._bonds

    def charges(self):
        """Partial charge per atom (surface groups, or edge_charges for a flake; else zero)."""
        if self.charge is not None:
            return np.array(self.charge, float)
        q = np.zeros(self.n_atoms)
        if self.edge_charges:
            qc, qh = self.edge_charges
            h_of = {}
            for i, j in self.bonds():
                if self.elements[i] == 'H':
                    h_of[j] = i
                elif self.elements[j] == 'H':
                    h_of[i] = j
            for c, h in h_of.items():
                q[c], q[h] = qc, qh
        return q

    def to_vlx_molecule(self):
        import veloxchem as vlx
        mol = vlx.Molecule(list(self.elements), 10.0 * np.asarray(self.xyz), 'angstrom')
        mol.set_charge(float(round(self.charges().sum(), 6)))
        mol.set_multiplicity(1)
        return mol

    def copy(self, **changes):
        """Deep copy (arrays and lists duplicated, bond cache reset), with field changes."""
        import copy as _copy
        out = _copy.deepcopy(self)
        out._bonds = None
        for k, v in changes.items():
            setattr(out, k, v)
        return out

    # ---- slab-specific
    def _slab(self):
        if self.kind != 'slab':
            raise AttributeError(f'only defined for slabs, this is a {self.kind}')

    @property
    def layers(self):
        return self.info['layers']

    @property
    def thickness(self):
        self._slab()
        return (self.layers - 1) * INTERLAYER

    @property
    def top_z(self):
        return self.thickness / 2

    @property
    def box_xy(self):
        self._slab()
        return self.cell[0], self.cell[1]

    @property
    def cells(self):
        return self.info.get('cells')


GraphiteSlab = CarbonModel   # backwards-compatible name


# ------------------------------------------------------------------ helpers

def _check(carbon_set, mobility, allow_unverified):
    if carbon_set not in CARBON_SETS:
        raise ValueError(f'unknown carbon_set {carbon_set!r}; options {list(CARBON_SETS)}')
    if mobility not in MOBILITY_MODES:
        raise ValueError(f'mobility must be one of {MOBILITY_MODES}')
    return check_status('carbon set', carbon_set, CARBON_SETS[carbon_set], allow_unverified)


def _wrap_tree(xyz, periodic, cell):
    """KD-tree with periodic dimensions; non-periodic ones get a box much larger than the data."""
    xyz = np.asarray(xyz, float)
    lo = xyz.min(0)
    box = np.empty(3)
    pts = xyz - lo
    for d in range(3):
        if periodic[d]:
            box[d] = cell[d]
            pts[:, d] = np.mod(xyz[:, d], cell[d])
            pts[pts[:, d] >= cell[d], d] = 0.0
        else:
            box[d] = pts[:, d].max() + 10.0
    return cKDTree(pts, boxsize=box)


def _bonds(elements, xyz, periodic, cell, cutoff=None):
    tree = _wrap_tree(xyz, periodic, cell)
    pairs = tree.query_pairs(cutoff or 0.16, output_type='ndarray')
    if len(pairs) == 0:
        return np.zeros((0, 2), int)
    el = np.asarray(elements)
    d = np.linalg.norm(_min_image(xyz[pairs[:, 0]] - xyz[pairs[:, 1]], periodic, cell), axis=1)
    # C-C up to 0.16 nm; C-H up to 0.12 nm; never H-H
    is_h = (el[pairs[:, 0]] == 'H') | (el[pairs[:, 1]] == 'H')
    hh = (el[pairs[:, 0]] == 'H') & (el[pairs[:, 1]] == 'H')
    keep = ~hh & np.where(is_h, d < 0.12, d < 0.16)
    return pairs[keep]


def _min_image(dv, periodic, cell):
    dv = np.array(dv, float)
    for d in range(3):
        if periodic[d]:
            dv[:, d] -= np.round(dv[:, d] / cell[d]) * cell[d]
    return dv


def _sheet(nx, ny, cc=CC_BOND):
    """One graphene layer in the orthorhombic (3cc x sqrt3 cc) cell, armchair along x."""
    ax, ay = 3 * cc, math.sqrt(3) * cc
    cell = np.array([[0, 0], [cc / 2, ay / 2], [1.5 * cc, ay / 2], [2 * cc, 0]])
    gx, gy = np.meshgrid(np.arange(nx) * ax, np.arange(ny) * ay, indexing='ij')
    return (cell[None, None] + np.stack([gx, gy], -1)[:, :, None]).reshape(-1, 2), (nx * ax, ny * ay)


STACKING_SHIFT = {'AB': lambda k: (k % 2), 'AA': lambda k: 0, 'ABC': lambda k: (k % 3)}


def lattice_cells(length_x, length_y, cc=CC_BOND, orientation='armchair_x'):
    """Orthorhombic cells closest to the requested lateral size (in the final x/y frame)."""
    ax, ay = 3 * cc, math.sqrt(3) * cc
    if orientation == 'zigzag_x':
        length_x, length_y = length_y, length_x
    return max(1, round(length_x / ax)), max(1, round(length_y / ay))


# ------------------------------------------------------------------ builders

def make_slab(carbon_set, layers=3, lateral_nm=(5.0, 5.0), *, stacking='AB',
              orientation='armchair_x', mobility='frozen', restraint_k=None,
              allow_unverified=False):
    """Periodic graphite slab (layers=1: graphene) centred on z = 0.

    stacking: 'AB' (Bernal, default), 'AA' or 'ABC' (rhombohedral); each layer is
    shifted by one C-C bond along the armchair direction.
    orientation: 'armchair_x' (C-C bonds along x) or 'zigzag_x' (zigzag rows along x).
    mobility: 'frozen' (mass 0), 'restrained' (position restraints, restraint_k,
    default 1000 kJ/mol/nm^2) or 'flexible' (GAFF2 ca bonded terms plus a weak anchor,
    restraint_k default 100 kJ/mol/nm^2, so the slab cannot drift through the z = 0
    cell boundary; restraint_k=0 removes the anchor).
    """
    warning = _check(carbon_set, mobility, allow_unverified)
    if layers < 1:
        raise ValueError('layers must be at least 1')
    if stacking not in STACKING_SHIFT:
        raise ValueError(f'stacking must be one of {list(STACKING_SHIFT)}')
    if orientation not in ('armchair_x', 'zigzag_x'):
        raise ValueError("orientation must be 'armchair_x' or 'zigzag_x'")
    nx, ny = lattice_cells(*lateral_nm, orientation=orientation)
    grid, (lx, ly) = _sheet(nx, ny)
    half = 0.5 * (layers - 1) * INTERLAYER
    xyz, layer = [], []
    for k in range(layers):
        xy = grid + [STACKING_SHIFT[stacking](k) * CC_BOND, 0.0]
        xy[:, 0] %= lx
        xyz.append(np.column_stack([xy, np.full(len(xy), half - k * INTERLAYER)]))
        layer.append(np.full(len(xy), k))
    xyz = np.vstack(xyz)
    if orientation == 'zigzag_x':
        xyz = xyz[:, [1, 0, 2]]
        lx, ly = ly, lx
    if restraint_k is None:
        restraint_k = {'restrained': 1000.0, 'flexible': 100.0}.get(mobility)
    return CarbonModel('slab', carbon_set, ['C'] * len(xyz), xyz, np.concatenate(layer),
                       (True, True, False), (lx, ly, None), mobility, restraint_k, None,
                       dict(layers=layers, cells=(nx, ny), stacking=stacking,
                            orientation=orientation), warning)


def make_flake(carbon_set, size_nm=(2.0, 2.0), *, layers=1, hydrogenate=True, edge_charges=None,
               mobility='flexible', restraint_k=None, allow_unverified=False):
    """Finite rectangular graphene flake (or AB-stacked multilayer), H-terminated edges.

    Carbons left with fewer than two C neighbours are removed; every carbon with two
    C neighbours gets one H (C-H 0.1086 nm, GAFF2 ca-ha). edge_charges=(qC, qH)
    puts charges on edge C/H pairs (e.g. OPLS-AA aromatic -0.115/+0.115); default none.
    Centred at the origin, not periodic.
    """
    warning = _check(carbon_set, mobility, allow_unverified)
    nx, ny = lattice_cells(size_nm[0] + 1.0, size_nm[1] + 1.0)
    grid, _ = _sheet(nx, ny)
    grid = grid - grid.mean(0)
    keep = (np.abs(grid[:, 0]) <= size_nm[0] / 2) & (np.abs(grid[:, 1]) <= size_nm[1] / 2)
    sheet = grid[keep]
    xyz, layer = [], []
    half = 0.5 * (layers - 1) * INTERLAYER
    for k in range(layers):
        xyz.append(np.column_stack([sheet + [(k % 2) * CC_BOND, 0], np.full(len(sheet), half - k * INTERLAYER)]))
        layer.append(np.full(len(sheet), k))
    xyz, layer = np.vstack(xyz), np.concatenate(layer)
    for _ in range(10):   # prune dangling carbons until every C has >= 2 C neighbours
        nb = _neighbour_counts(xyz)
        if np.all(nb >= 2):
            break
        xyz, layer = xyz[nb >= 2], layer[nb >= 2]
    elements = ['C'] * len(xyz)
    if hydrogenate:
        tree = cKDTree(xyz)
        hs, hl = [], []
        for i, p in enumerate(xyz):
            nbr = [j for j in tree.query_ball_point(p, 0.16) if j != i and abs(xyz[j, 2] - p[2]) < 0.05]
            if len(nbr) == 2:
                u = 2 * p - xyz[nbr[0]] - xyz[nbr[1]]
                u[2] = 0.0
                hs.append(p + CH_BOND * u / np.linalg.norm(u))
                hl.append(layer[i])
        if hs:
            xyz = np.vstack([xyz, hs])
            layer = np.concatenate([layer, hl])
            elements += ['H'] * len(hs)
    if mobility == 'restrained' and restraint_k is None:
        restraint_k = 1000.0
    return CarbonModel('flake', carbon_set, elements, xyz, layer, (False, False, False),
                       (None, None, None), mobility, restraint_k, edge_charges,
                       dict(layers=layers, size_nm=tuple(size_nm), hydrogenated=hydrogenate), warning)


def _neighbour_counts(xyz):
    tree = cKDTree(xyz)
    return np.array([sum(1 for j in tree.query_ball_point(p, 0.16)
                         if j != i and abs(xyz[j, 2] - p[2]) < 0.05) for i, p in enumerate(xyz)])


def make_nanotube(carbon_set, n, m, cells=4, *, mobility='frozen', restraint_k=None,
                  allow_unverified=False):
    """Single-wall (n, m) nanotube, axis along z, periodic over `cells` translational units.

    Radius |C_h| / 2pi with C_h = n a1 + m a2 (a = sqrt3 * 0.142 nm). Centred on the
    z axis (x = y = 0). Periodic length along z = cells * |T|.
    """
    warning = _check(carbon_set, mobility, allow_unverified)
    if n < 1 or m < 0 or m > n:
        raise ValueError('use 0 <= m <= n, n >= 1')
    a = math.sqrt(3) * CC_BOND
    a1 = np.array([a * math.sqrt(3) / 2, a / 2])
    a2 = np.array([a * math.sqrt(3) / 2, -a / 2])
    basis = [np.zeros(2), (a1 + a2) / 3]
    dr = math.gcd(2 * m + n, 2 * n + m)
    t1, t2 = (2 * m + n) // dr, -(2 * n + m) // dr
    ch = n * a1 + m * a2
    tv = t1 * a1 + t2 * a2
    lch, lt = np.linalg.norm(ch), np.linalg.norm(tv)
    uch, ut = ch / lch, tv / lt
    per_cell = 4 * (n * n + m * m + n * m) // dr
    span = int(abs(n) + abs(m) + abs(t1) * cells + abs(t2) * cells + 4)
    pts = []
    eps = 1e-6
    for i in range(-span, span + 1):
        for j in range(-span, span + 1):
            for b in basis:
                r = i * a1 + j * a2 + b
                u, v = r @ uch / lch, r @ ut / lt
                if -eps <= u < 1 - eps and -eps <= v < cells - eps:
                    pts.append((u, v))
    pts = np.array(pts)
    if len(pts) != per_cell * cells:
        raise RuntimeError(f'nanotube construction gave {len(pts)} atoms, expected {per_cell * cells}')
    radius = lch / (2 * math.pi)
    phi = 2 * math.pi * pts[:, 0]
    xyz = np.column_stack([radius * np.cos(phi), radius * np.sin(phi), pts[:, 1] * lt])
    if mobility == 'restrained' and restraint_k is None:
        restraint_k = 1000.0
    lz = cells * lt
    return CarbonModel('nanotube', carbon_set, ['C'] * len(xyz), xyz, np.zeros(len(xyz), int),
                       (False, False, True), (None, None, lz), mobility, restraint_k, None,
                       dict(n=n, m=m, cells=cells, radius_nm=radius, length_nm=lz,
                            atoms_per_cell=per_cell), warning)


# ------------------------------------------------------------------ force field / placement

def graphite_forcefield_xml(carbon_set, edge_charges=None):
    """OpenMM ForceField XML with one-atom residues: GRA (carbon, LJ of the carbon set),
    GRH (hydrogen, GAFF2 ha LJ) and GRO (oxygen, placeholder). All have charge 0;
    per-atom charges and group LJ are set after createSystem (forces.add_carbon_mobility),
    and bonded terms are added there too. edge_charges is accepted for compatibility."""
    cs = CARBON_SETS[carbon_set]
    return ('<ForceField><AtomTypes>'
            f'<Type name="EB-GRA-C" class="EB-GRA-C" element="C" mass="{CARBON_MASS}"/>'
            f'<Type name="EB-GRH-H" class="EB-GRH-H" element="H" mass="{H_MASS}"/>'
            f'<Type name="EB-GRO-O" class="EB-GRO-O" element="O" mass="{ELEMENT_MASS["O"]}"/>'
            '</AtomTypes><Residues>'
            f'<Residue name="{GRAPHITE_RESNAME}"><Atom name="C" type="EB-GRA-C"/></Residue>'
            '<Residue name="GRH"><Atom name="H" type="EB-GRH-H"/></Residue>'
            '<Residue name="GRO"><Atom name="O" type="EB-GRO-O"/></Residue>'
            '</Residues><NonbondedForce coulomb14scale="0.8333333333333334" lj14scale="0.5">'
            f'<Atom type="EB-GRA-C" charge="0" sigma="{cs["sigma"]!r}" epsilon="{cs["epsilon"]!r}"/>'
            f'<Atom type="EB-GRH-H" charge="0" sigma="{GAFF2_HA[0]!r}" epsilon="{GAFF2_HA[1]!r}"/>'
            '<Atom type="EB-GRO-O" charge="0" sigma="0.3" epsilon="0"/>'
            '</NonbondedForce></ForceField>')


def residue_names(model):
    """GRA / GRH / GRO per atom (by element)."""
    return [{'C': GRAPHITE_RESNAME, 'H': 'GRH', 'O': 'GRO'}[e] for e in model.elements]
