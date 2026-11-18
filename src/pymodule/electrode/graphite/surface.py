"""Surface features on carbon models: step edges, vacancies, Stone-Wales defects and
oxygen groups, the pieces of a glassy-carbon-like electrode surface.

    surf = add_step(slab, edge='armchair', keep_fraction=0.5)      # outer-layer ribbon, H edges
    surf = add_defects(surf, vacancies=2, stone_wales=1, seed=1)
    surf = add_groups(surf, {'carboxylate': 2, 'epoxide': 2}, seed=1)

Every call returns a new CarbonModel and acts on the outer layer of `face`: 'top'
(layer 0, normal +z) or 'bottom' (last layer, normal -z). model.features lists
what was done where.

Geometry. The lattice stays ideal and flat; sp3 anchors are not puckered. Defects
can be relaxed with the flexible GAFF2 terms (relax=True, all other atoms fixed),
the one place where a force-field geometry optimisation is used. Vacancies are not
reconstructed (no 5-9 bond formation) and their dangling carbons stay bare.

Groups (params.SURFACE_GROUPS) are closed-shell, so their RESP charges come from a
closed-shell cluster. A single basal OH or a lone edge C=O leaves an unpaired pi
electron; basal OH therefore comes as a para pair on one ring, and edge carbonyls
as an armchair o-quinone. Each group is parametrised once per site type (basal,
armchair edge, zigzag edge) on an H-capped graphene cluster cut from the same
lattice (surface_templates.py): GAFF2 types and bonded terms and RESP charges
(HF/6-31G*) from VeloxChem's MMForceFieldGenerator, after relaxing the group atoms
with GAFF2 (skeleton fixed). Every cluster atom that also exists in the slab (the
group, and the skeleton within a few bonds of it, real edge H included) takes its
cluster charge; the charge of each H that only caps the cut cluster is added to
its carbon. The fragment then carries exactly the group's formal charge (anions
delocalise into the ring, so a small fragment would not). All other skeleton atoms
stay neutral. Group LJ is GAFF2; skeleton atoms keep the carbon set's LJ. Group
atoms always move, also on a frozen skeleton.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

from .graphite import CH_BOND, _min_image, _wrap_tree
from .params import CC_BOND, SURFACE_GROUPS, check_status

from ...environment import get_data_path  # noqa: E402

TEMPLATE_FILE = get_data_path() / 'electrode' / 'surface_groups.json'
MATCH_TOL = 0.02          # nm, template reference atom -> model atom


# ------------------------------------------------------------------ graph helpers

def _neighbours(model):
    nbrs = [[] for _ in range(model.n_atoms)]
    for i, j in model.bonds():
        nbrs[int(i)].append(int(j))
        nbrs[int(j)].append(int(i))
    return nbrs


def _outer_layer(model, face):
    if face not in ('top', 'bottom'):
        raise ValueError("face must be 'top' or 'bottom'")
    return 0 if face == 'top' else model.info['layers'] - 1


def _normal(face):
    return np.array([0.0, 0.0, 1.0 if face == 'top' else -1.0])


def _vec(model, i, j):
    """Minimum-image vector from atom i to atom j."""
    return _min_image((model.xyz[j] - model.xyz[i])[None], model.periodic, model.cell)[0]


def edge_carbons(model, face='top'):
    """{atom: 'armchair' | 'zigzag'} for skeleton carbons of the outer layer with two
    carbon neighbours (armchair: bonded to another edge carbon)."""
    nbrs = _neighbours(model)
    lay = _outer_layer(model, face)
    skel = model.skeleton
    el = model.elements
    edge = [i for i in range(model.n_atoms) if skel[i] and el[i] == 'C' and model.layer[i] == lay
            and sum(1 for j in nbrs[i] if el[j] == 'C' and skel[j]) == 2]
    es = set(edge)
    return {i: ('armchair' if any(j in es for j in nbrs[i] if el[j] == 'C') else 'zigzag')
            for i in edge}


def _outward(model, a, nbrs):
    """In-plane unit vector pointing away from an edge carbon's two C neighbours."""
    cs = [j for j in nbrs[a] if model.elements[j] == 'C' and model.skeleton[j]]
    u = -(_vec(model, a, cs[0]) + _vec(model, a, cs[1]))
    u[2] = 0.0
    return u / np.linalg.norm(u)


def _site_frame(model, anchors, kind, face, nbrs):
    """origin, (u, v, n): edge sites -> u outward; basal -> u from anchor 0 to anchor 1."""
    n = _normal(face)
    a0 = model.xyz[anchors[0]]
    if kind == 'edge':
        outs = [_outward(model, a, nbrs) for a in anchors]
        if len(anchors) == 1:
            origin, u = a0, outs[0]
        else:
            origin = a0 + 0.5 * _vec(model, anchors[0], anchors[1])
            u = outs[0] + outs[1]
            u /= np.linalg.norm(u)
    else:
        d = _vec(model, anchors[0], anchors[1])
        origin = a0 + 0.5 * d
        u = d.copy()
        u[2] = 0.0
        u /= np.linalg.norm(u)
    return origin, np.array([u, np.cross(n, u), n])


def _shell(nbrs, anchors, depth, skip=()):
    """Atoms within `depth` bonds of the anchors (anchors included), excluding `skip`."""
    seen = {int(a): 0 for a in anchors}
    frontier = list(seen)
    for d in range(1, depth + 1):
        nxt = []
        for a in frontier:
            for b in nbrs[a]:
                if b not in seen and b not in skip:
                    seen[b] = d
                    nxt.append(b)
        frontier = nxt
    return seen


# ------------------------------------------------------------------ steps and defects

def add_step(model, edge='armchair', keep_fraction=0.5, face='top', offset_nm=0.0):
    """Keep a ribbon of the outer layer, exposing two step edges of the given type.

    The ribbon runs along x or y, whichever is the `edge` direction for the slab's
    orientation (armchair edges run along the armchair direction). Carbons left with
    one neighbour are removed, and every edge carbon gets an H (C-H 0.1086 nm, GAFF2
    ha LJ, charge 0) in the same layer.
    """
    if model.kind != 'slab':
        raise ValueError('add_step needs a periodic slab')
    if edge not in ('armchair', 'zigzag'):
        raise ValueError("edge must be 'armchair' or 'zigzag'")
    if not 0.1 <= keep_fraction <= 0.9:
        raise ValueError('keep_fraction must be between 0.1 and 0.9')
    lay = _outer_layer(model, face)
    orient = model.info.get('orientation', 'armchair_x')
    run_axis = 0 if (orient == 'armchair_x') == (edge == 'armchair') else 1
    cut_axis = 1 - run_axis
    L = model.cell[cut_axis]
    x = np.mod(model.xyz[:, cut_axis] - offset_nm, L)
    keep = ~((model.layer == lay) & ~(x < keep_fraction * L))
    out = _take(model, np.flatnonzero(keep))
    out = _prune(out, lay)
    out = _hydrogenate(out, lay)
    out.features.append(dict(kind='step', edge=edge, face=face, layer=lay,
                             keep_fraction=keep_fraction, run_axis='xy'[run_axis],
                             removed_carbons=int(model.n_carbon - out.n_carbon)))
    out.info['removed_layer_fraction'] = model.info.get('removed_layer_fraction', 0.0) + 1 - keep_fraction
    return out


def add_defects(model, vacancies=0, divacancies=0, stone_wales=0, face='top', seed=None,
                min_separation_nm=1.0, relax=True):
    """Point defects in the outer layer, at random interior sites at least
    min_separation_nm apart (and from edges and groups).

    vacancies / divacancies remove one carbon / a bonded pair; stone_wales rotates a
    C-C bond by 90 degrees about its midpoint (55-77 rings). relax=True minimises the
    outer layer with the flexible GAFF2 terms while every other atom stays fixed.
    """
    rng = np.random.default_rng(seed)
    lay = _outer_layer(model, face)
    m = model.copy()
    taken = _feature_points(m, face)
    placed = []
    for kind, count in (('stone_wales', stone_wales), ('divacancy', divacancies),
                        ('vacancy', vacancies)):
        for _ in range(count):
            nbrs = _neighbours(m)
            cand = _interior_carbons(m, lay, nbrs)
            rng.shuffle(cand)
            for a in cand:
                if _too_close(m, m.xyz[a], taken, min_separation_nm):
                    continue
                b = None
                if kind != 'vacancy':
                    opts = [j for j in nbrs[a] if j in set(cand)]
                    if not opts:
                        continue
                    b = int(rng.choice(opts))
                pos = m.xyz[a].copy()
                if kind == 'stone_wales':
                    m = _stone_wales(m, a, b)
                else:
                    drop = {a} | ({b} if b is not None else set())
                    m = _take(m, [i for i in range(m.n_atoms) if i not in drop])
                taken.append(pos)
                placed.append(dict(kind=kind, face=face, layer=lay, position_nm=pos.tolist()))
                break
            else:
                raise RuntimeError(f'no room for another {kind} (min_separation_nm '
                                   f'{min_separation_nm}); use a larger slab or fewer defects')
    if relax and placed:
        m = _relax_layer(m, lay)
    m.features += placed
    return m


def _interior_carbons(m, lay, nbrs):
    skel = m.skeleton
    return [i for i in range(m.n_atoms) if skel[i] and m.elements[i] == 'C' and m.layer[i] == lay
            and sum(1 for j in nbrs[i] if m.elements[j] == 'C') == 3
            and all(sum(1 for k in nbrs[j] if m.elements[k] == 'C') == 3
                    for j in nbrs[i] if m.elements[j] == 'C')]


def _stone_wales(m, a, b):
    """Rotate bond a-b by 90 degrees in plane about its midpoint."""
    m = m.copy()
    d = _vec(m, a, b)
    mid = m.xyz[a] + 0.5 * d
    rot = np.array([-d[1], d[0], 0.0])
    m.xyz[a] = mid - 0.5 * rot
    m.xyz[b] = mid + 0.5 * rot
    _wrap_inplace(m)
    m._bonds = None
    return m


def _relax_layer(m, lay):
    """Minimise layer `lay` (skeleton) with the flexible GAFF2 terms; everything else fixed."""
    import openmm as mm
    from openmm import unit
    from .export import carbon_system
    flex = m.copy(mobility='flexible', restraint_k=0)
    system, _, x0, _ = carbon_system(flex)
    for i in range(m.n_atoms):
        if not (m.skeleton[i] and m.layer[i] == lay):
            system.setParticleMass(i, 0.0)
    ctx = mm.Context(system, mm.VerletIntegrator(0.001), mm.Platform.getPlatformByName('Reference'))
    ctx.setPositions(x0 * unit.nanometer)
    mm.LocalEnergyMinimizer.minimize(ctx, 1.0, 5000)
    x1 = ctx.getState(getPositions=True).getPositions(asNumpy=True).value_in_unit(unit.nanometer)
    out = m.copy()
    out.xyz = out.xyz + (x1 - x0)
    _wrap_inplace(out)
    out.info.setdefault('relaxed', []).append(dict(layer=int(lay), max_shift_nm=float(
        np.abs(x1 - x0).max()), method='GAFF2 flexible terms, other atoms fixed'))
    return out


def _wrap_inplace(m):
    for d in range(3):
        if m.periodic[d]:
            m.xyz[:, d] = np.mod(m.xyz[:, d], m.cell[d])


def _take(model, keep):
    """Model restricted to atoms `keep` (indices), with group data remapped."""
    keep = np.asarray(sorted(int(k) for k in keep))
    new = {int(o): n for n, o in enumerate(keep)}
    m = model.copy()
    m.elements = [model.elements[i] for i in keep]
    m.xyz = model.xyz[keep].copy()
    m.layer = model.layer[keep].copy()
    for name in ('charge', 'lj', 'group', 'fragment'):
        v = getattr(model, name)
        if v is not None:
            setattr(m, name, np.asarray(v)[keep].copy())
    if model.group_bonds:
        m.group_bonds = [[new[i], new[j]] for i, j in model.group_bonds if i in new and j in new]
    if model.terms:
        m.terms = {k: [[new[int(x)] for x in t[:_NATOMS[k]]] + list(t[_NATOMS[k]:])
                       for t in v if all(int(x) in new for x in t[:_NATOMS[k]])]
                   for k, v in model.terms.items()}
    m._bonds = None
    return m


_NATOMS = {'bonds': 2, 'angles': 3, 'torsions': 4, 'rb': 4}


def _prune(m, lay):
    """Remove outer-layer carbons with fewer than two carbon neighbours, repeatedly."""
    for _ in range(20):
        nbrs = _neighbours(m)
        bad = [i for i in range(m.n_atoms) if m.skeleton[i] and m.elements[i] == 'C'
               and m.layer[i] == lay and sum(1 for j in nbrs[i] if m.elements[j] == 'C') < 2]
        if not bad:
            return m
        m = _take(m, [i for i in range(m.n_atoms) if i not in set(bad)])
    raise RuntimeError('pruning did not converge')


def _hydrogenate(m, lay):
    """An H on every outer-layer carbon with two carbon neighbours and no H."""
    nbrs = _neighbours(m)
    hs = []
    for i in range(m.n_atoms):
        if not (m.skeleton[i] and m.elements[i] == 'C' and m.layer[i] == lay):
            continue
        cs = [j for j in nbrs[i] if m.elements[j] == 'C']
        if len(cs) == 2 and not any(m.elements[j] == 'H' for j in nbrs[i]):
            hs.append(m.xyz[i] + CH_BOND * _outward(m, i, nbrs))
    if not hs:
        return m
    out = m.copy()
    out.elements = m.elements + ['H'] * len(hs)
    out.xyz = np.vstack([m.xyz, hs])
    out.layer = np.concatenate([m.layer, np.full(len(hs), lay)])
    for name, fill in (('charge', 0.0), ('group', -1), ('fragment', -1)):
        v = getattr(m, name)
        if v is not None:
            setattr(out, name, np.concatenate([v, np.full(len(hs), fill)]))
    if m.lj is not None:
        out.lj = np.vstack([m.lj, np.full((len(hs), 2), np.nan)])
    _wrap_inplace(out)
    out._bonds = None
    return out


def _feature_points(m, face):
    pts = [np.asarray(f['position_nm']) for f in m.features if 'position_nm' in f]
    return pts


def _too_close(m, p, points, r):
    for q in points:
        if np.linalg.norm(_min_image((p - q)[None], m.periodic, m.cell)[0]) < r:
            return True
    return False


# ------------------------------------------------------------------ oxygen groups

def ideal_group_atoms(group, origin, frame, anchor_xyz, twist_deg=30.0):
    """Ideal group geometry: list of (element, name, xyz) for the atoms the group adds.

    frame = (u, v, n). anchor_xyz: anchor positions (edge: outward is u for one anchor;
    o-quinone uses each anchor's own outward direction, passed as anchor_xyz[k][1]).
    """
    u, v, n = frame
    t = math.radians(twist_deg)
    w = math.cos(t) * v + math.sin(t) * n          # carboxyl plane: u, w

    def dirn(deg, a, b):
        r = math.radians(deg)
        return math.cos(r) * a + math.sin(r) * b

    if group in ('phenol', 'phenolate'):
        a = anchor_xyz[0]
        o = a + (0.1364 if group == 'phenol' else 0.1280) * u
        atoms = [('O', 'O1', o)]
        if group == 'phenol':
            atoms.append(('H', 'HO1', o + 0.0974 * dirn(71.5, u, v)))
        return atoms
    if group in ('carboxyl', 'carboxylate'):
        a = anchor_xyz[0]
        c = a + (0.1480 if group == 'carboxyl' else 0.1510) * u
        if group == 'carboxylate':
            return [('C', 'C1', c), ('O', 'O1', c + 0.1260 * dirn(62.5, u, w)),
                    ('O', 'O2', c + 0.1260 * dirn(62.5, u, -w))]
        o1 = c + 0.1210 * dirn(57.0, u, w)
        o2 = c + 0.1340 * dirn(67.0, u, -w)
        e = (c - o2) / np.linalg.norm(c - o2)
        p = np.cross(np.cross(u, w), e)                  # in-plane perpendicular to O2-C
        h = [o2 + 0.0974 * (math.cos(math.radians(107)) * e + s * math.sin(math.radians(107)) * p)
             for s in (1, -1)]
        h = min(h, key=lambda x: np.linalg.norm(x - o1))   # syn: H on the C=O side
        return [('C', 'C1', c), ('O', 'O1', o1), ('O', 'O2', o2), ('H', 'HO2', h)]
    if group == 'o_quinone':
        return [('O', f'O{k + 1}', a + 0.1220 * out) for k, (a, out) in enumerate(anchor_xyz)]
    if group == 'epoxide':
        h = math.sqrt(0.1435 ** 2 - (0.5 * CC_BOND) ** 2)
        return [('O', 'O1', origin + h * n)]
    if group == 'hydroxyl_pair':
        atoms = []
        for k, a in enumerate(anchor_xyz):
            tdir = a - origin
            tdir[2] = 0.0
            tdir /= np.linalg.norm(tdir)
            o = a + 0.1430 * n
            atoms += [('O', f'O{k + 1}', o), ('H', f'HO{k + 1}', o + 0.0974 * dirn(71.5, n, tdir))]
        return atoms
    raise ValueError(f'unknown group {group!r}')


def candidate_sites(model, group, face, nbrs=None):
    """Anchor tuples where `group` could go, by site type (no environment check yet)."""
    spec = SURFACE_GROUPS[group]
    nbrs = nbrs or _neighbours(model)
    lay = _outer_layer(model, face)
    if spec['site'] == 'basal':
        inner = set(_interior_carbons(model, lay, nbrs))
        if group == 'epoxide':
            return [(i, j) for i in inner for j in nbrs[i] if j in inner and i < j]
        pairs = []
        for i in inner:          # para partner: 3 bonds away, 0.284 nm
            for j in _shell(nbrs, [i], 3):
                if j in inner and i < j and abs(np.linalg.norm(_vec(model, i, j)) - 2 * CC_BOND) < 0.01:
                    pairs.append((i, j))
        return pairs
    edges = edge_carbons(model, face)
    with_h = {i: t for i, t in edges.items() if any(model.elements[j] == 'H' for j in nbrs[i])}
    if group == 'o_quinone':
        return [(i, j, 'armchair') for i in with_h for j in nbrs[i]
                if j in with_h and i < j and with_h[i] == with_h[j] == 'armchair']
    return [(i, t) for i, t in with_h.items()]


def load_templates(path=TEMPLATE_FILE):
    if not Path(path).exists():
        return {}
    return json.loads(Path(path).read_text())['templates']


def add_groups(model, counts, face='top', seed=None, min_separation_nm=0.8, relax=True,
               templates=None, allow_unverified=False):
    """Place oxygen groups (params.SURFACE_GROUPS) at random free sites of the outer layer.

    counts: {group: number}. Edge groups need edges (add_step, or a flake); basal
    groups need interior carbons. A site is used only if the template's whole cluster
    environment is found there unchanged and no atom of it belongs to another group's
    charge fragment; sites are also at least min_separation_nm (site to site) from
    other groups and from defects. relax=True minimises the new group
    atoms (GAFF2 + charges) with the skeleton fixed, which settles them against the
    local surface (e.g. the terrace below a step).
    """
    rng = np.random.default_rng(seed)
    templates = load_templates() if templates is None else templates
    m = _with_group_arrays(model.copy())
    warnings = []
    taken = _feature_points(m, face)
    for group, count in counts.items():
        if group not in SURFACE_GROUPS:
            raise ValueError(f'unknown group {group!r}; options {list(SURFACE_GROUPS)}')
        w = check_status('surface group', group, SURFACE_GROUPS[group], allow_unverified)
        if w and w not in warnings:
            warnings.append(w)
        for _ in range(count):
            nbrs = _neighbours(m)
            cands = candidate_sites(m, group, face, nbrs)
            order = rng.permutation(len(cands))
            for k in order:
                site = cands[k]
                if group in ('o_quinone',):
                    anchors, site_type = site[:2], site[2]
                elif SURFACE_GROUPS[group]['site'] == 'basal':
                    anchors, site_type = site, 'basal'
                else:
                    anchors, site_type = site[:1], site[1]
                key = f'{group}/{site_type}'
                if key not in templates:
                    raise KeyError(f'no template {key!r} in {TEMPLATE_FILE.name}; generate it with '
                                   'surface_templates.generate_templates')
                origin, frame = _site_frame(m, anchors, 'basal' if site_type == 'basal' else 'edge',
                                            face, nbrs)
                if _too_close(m, origin, taken, min_separation_nm):
                    continue
                placed = _apply_template(m, templates[key], anchors, origin, frame, nbrs, face)
                if placed is None:
                    continue
                m = placed
                taken.append(origin)
                gid = int(m.group.max())
                m.features.append(dict(kind='group', group=group, site=site_type, face=face,
                                       group_id=gid, position_nm=origin.tolist(),
                                       formal_charge=SURFACE_GROUPS[group]['charge']))
                break
            else:
                raise RuntimeError(f'no free {group} site left (min_separation_nm '
                                   f'{min_separation_nm}); use a larger slab, a step, or fewer groups')
    if relax and (~m.skeleton).any():
        m = relax_groups(m)
    m.info.setdefault('warnings', [])
    m.info['warnings'] += [w for w in warnings if w not in m.info['warnings']]
    return m


def _with_group_arrays(m):
    n = m.n_atoms
    if m.charge is None:
        m.charge = m.charges()
    if m.lj is None:
        m.lj = np.full((n, 2), np.nan)
    if m.group is None:
        m.group = np.full(n, -1)
    if m.fragment is None:
        m.fragment = np.full(n, -1)
    if m.terms is None:
        m.terms = {'bonds': [], 'angles': [], 'torsions': [], 'rb': []}
    if m.group_bonds is None:
        m.group_bonds = []
    return m


def _apply_template(m, tpl, anchors, origin, frame, nbrs, face):
    """Place one group from its template; None if the local environment does not match."""
    refs = tpl['skeleton_refs']
    ga = tpl['group_atoms']
    replaced = [j for a in anchors for j in nbrs[a] if m.elements[j] == 'H'] \
        if tpl['site'] != 'basal' else []
    rep = set(replaced)
    env = _shell(nbrs, anchors, 2, skip=rep)
    tree = _wrap_tree(m.xyz, m.periodic, m.cell)
    lo = m.xyz.min(0)
    for mirror in (1.0, -1.0):
        u, v, n = frame
        R = np.array([u, mirror * v, n]).T          # frame -> lab
        mapping = []
        for r in refs:
            p = origin + R @ np.asarray(r['pos'])
            q = p - lo
            for d in range(3):
                if m.periodic[d]:
                    q[d] = np.mod(p[d], m.cell[d])
            d, j = tree.query(q)
            if (d > MATCH_TOL or m.elements[j] != r['element'] or not m.skeleton[j]
                    or m.fragment[j] >= 0 or m.charge[j] != 0 or j in rep):
                break               # missing atom, another group or its charge fragment
            mapping.append(int(j))
        else:
            if (len(set(mapping)) == len(refs) and set(env) <= set(mapping)
                    and _same_bonding(m, refs, mapping, nbrs, rep)):
                break
    else:
        return None
    # new atoms
    out = m.copy()
    gid = int(out.group.max()) + 1
    n0 = out.n_atoms
    new_xyz = [origin + R @ np.asarray(g['pos']) for g in ga]
    lay = int(m.layer[anchors[0]])
    out.elements = out.elements + [g['element'] for g in ga]
    out.xyz = np.vstack([out.xyz, new_xyz])
    out.layer = np.concatenate([out.layer, np.full(len(ga), lay)])
    out.charge = np.concatenate([out.charge, [g['charge'] for g in ga]])
    out.lj = np.vstack([out.lj, [[g['sigma'], g['epsilon']] for g in ga]])
    out.group = np.concatenate([out.group, np.full(len(ga), gid)])
    out.fragment = np.concatenate([out.fragment, np.full(len(ga), gid)])
    for r, j in zip(refs, mapping):
        out.charge[j] = r['charge']
        out.fragment[j] = gid
    index = mapping + list(range(n0, n0 + len(ga)))
    for name, terms in tpl['terms'].items():
        k = _NATOMS[name]
        out.terms[name] += [[index[x] for x in t[:k]] + list(t[k:]) for t in terms]
    out.group_bonds += [[index[i], index[j]] for i, j in tpl['graph_bonds']]
    _wrap_inplace(out)
    out._bonds = None
    if replaced:
        out = _take(out, [i for i in range(out.n_atoms) if i not in set(replaced)])
    return out


def _same_bonding(m, refs, mapping, nbrs, replaced):
    """Every matched atom has the cluster's bonded neighbours, and a carbon wherever the
    cluster had a cap H (so the slab is the cluster's chemistry, not just its shape)."""
    for r, j in zip(refs, mapping):
        have = {x for x in nbrs[j] if x not in replaced}
        want = {mapping[k] for k in r['neighbours']}
        extra = have - want
        if not want <= have or len(extra) != r['n_caps'] or any(
                m.elements[x] != 'C' or not m.skeleton[x] for x in extra):
            return False
    return True


def relax_groups(m, tolerance=1.0):
    """Minimise the group atoms (GAFF2 + charges, vacuum) with the skeleton fixed."""
    import openmm as mm
    from openmm import unit
    from .export import carbon_system
    system, _, x0, _ = carbon_system(m.copy(mobility='frozen', restraint_k=None))
    ctx = mm.Context(system, mm.VerletIntegrator(0.001), mm.Platform.getPlatformByName('Reference'))
    ctx.setPositions(x0 * unit.nanometer)
    ctx.applyConstraints(1e-6)
    mm.LocalEnergyMinimizer.minimize(ctx, tolerance, 5000)
    x1 = ctx.getState(getPositions=True).getPositions(asNumpy=True).value_in_unit(unit.nanometer)
    out = m.copy()
    out.xyz = out.xyz + (x1 - x0)
    _wrap_inplace(out)
    out.info['groups_relaxed'] = dict(max_shift_nm=float(np.abs(x1 - x0).max()),
                                      method='GAFF2 + template charges, vacuum, skeleton fixed')
    return out
