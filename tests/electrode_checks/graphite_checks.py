"""Self-checks for the graphite builder (graphite only; CPU). Run: python tests/test_graphite.py

Group checks use data/surface_groups.json if present (production, from the cluster), else the
STO-3G test templates in tests/data.
"""
import math
import sys
import tempfile
from pathlib import Path

import numpy as np
import openmm as mm
from openmm import unit

HERE = Path(__file__).resolve().parent

from veloxchem.electrode.graphite import (GraphiteBuilder, add_defects, add_groups, add_nbfix_correction,  # noqa: E402
                              add_step, carbon_system, make_flake, make_nanotube, make_slab,
                              write_carbon)
from veloxchem.electrode.graphite.params import CARBON_SETS, SURFACE_GROUPS  # noqa: E402
from veloxchem.electrode.graphite.surface import (_apply_template, _neighbours, _site_frame,  # noqa: E402
                                      edge_carbons, load_templates)
from veloxchem.electrode.graphite.surface_templates import TEST_TEMPLATE_FILE  # noqa: E402

KJ = unit.kilojoule_per_mole
NM = unit.nanometer


def check(name, cond, detail=''):
    print(('PASS ' if cond else 'FAIL ') + name + (f'  [{detail}]' if detail else ''))
    if not cond:
        raise SystemExit(1)


class raises:
    def __init__(self, name):
        self.name = name

    def __enter__(self):
        return self

    def __exit__(self, t, e, tb):
        check(self.name, t is not None and issubclass(t, Exception),
              f'{t.__name__}: {str(e)[:70]}' if t else 'no exception')
        return True


def test_refusals():
    with raises('unverified carbon set refused'):
        make_slab('werder')
    check('wu_aluru (read from paper) accepted', make_slab('wu_aluru', 2, (2.6, 2.5)).warning is None)


def test_nbfix():
    cs = CARBON_SETS['wu_aluru']
    s = mm.System()
    for _ in range(2):
        s.addParticle(12.0)
    s.setDefaultPeriodicBoxVectors(mm.Vec3(5, 0, 0), mm.Vec3(0, 5, 0), mm.Vec3(0, 0, 5))
    nb = mm.NonbondedForce()
    nb.setNonbondedMethod(mm.NonbondedForce.CutoffPeriodic)
    nb.setCutoffDistance(2.0)
    nb.setUseDispersionCorrection(False)
    nb.setReactionFieldDielectric(1.0)
    nb.addParticle(0.0, cs['sigma'], cs['epsilon'])           # graphite C
    nb.addParticle(0.0, 0.316557, 0.650194)                   # SPC/E O LJ
    s.addForce(nb)
    add_nbfix_correction(s, [0], [1], [], 'wu_aluru', cutoff=2.0, switch=None)
    ctx = mm.Context(s, mm.VerletIntegrator(0.001), mm.Platform.getPlatformByName('Reference'))
    sco, eco = cs['pairs']['CO']
    for r in (0.32, 0.36, 0.5):
        ctx.setPositions([mm.Vec3(1, 1, 1), mm.Vec3(1, 1, 1 + r)])
        e = ctx.getState(getEnergy=True).getPotentialEnergy().value_in_unit(KJ)
        ref = 4 * eco * ((sco / r) ** 12 - (sco / r) ** 6)
        check(f'NBFIX C-O pair energy at {r} nm', math.isclose(e, ref, rel_tol=1e-6, abs_tol=1e-9),
              f'{e:.6f} vs {ref:.6f}')


def _stacked_fraction(slab, k0, k1):
    """Fraction of layer-k1 atoms with a layer-k0 atom directly above/below (xy < 0.01 nm)."""
    from scipy.spatial import cKDTree
    a, b = slab.xyz[slab.layer == k0, :2], slab.xyz[slab.layer == k1, :2]
    d, _ = cKDTree(a, boxsize=slab.box_xy).query(np.mod(b, slab.box_xy))
    return float((d < 0.01).mean())


def test_geometry():
    for st, want in (('AB', 0.5), ('AA', 1.0), ('ABC', 0.5)):
        sl = make_slab('steele', 4, (2.6, 2.5), stacking=st)
        f01, f02, f03 = (_stacked_fraction(sl, 0, k) for k in (1, 2, 3))
        ok = math.isclose(f01, want) and (f02 == (1.0 if st in ('AB', 'AA') else 0.5)) and \
            (f03 == (1.0 if st == 'ABC' else want))
        check(f'{st} stacking registry', ok, f'layer0 over 1/2/3: {f01:.2f} {f02:.2f} {f03:.2f}')
    for o, ax in (('armchair_x', 0), ('zigzag_x', 1)):
        sl = make_slab('steele', 1, (2.6, 2.5), orientation=o)
        b = sl.bonds()
        dv = sl.xyz[b[:, 0]] - sl.xyz[b[:, 1]]
        dv[:, :2] -= np.round(dv[:, :2] / sl.box_xy) * sl.box_xy
        L = np.linalg.norm(dv, axis=1)
        par = np.abs(dv[:, ax]) / L > 0.999
        check(f'{o}: 3 bonds of 0.142 nm per C, 1/3 of bonds along {"xy"[ax]}',
              np.allclose(L, 0.142) and np.all(np.bincount(b.ravel()) == 3) and
              math.isclose(par.mean(), 1 / 3), f'{len(b)} bonds')
    fl = make_flake('steele', (2.0, 1.6))
    el, b = np.array(fl.elements), fl.bonds()
    deg = np.bincount(b.ravel(), minlength=fl.n_atoms)
    ch = [np.linalg.norm(fl.xyz[i] - fl.xyz[j]) for i, j in b if 'H' in (el[i], el[j])]
    check('flake: every C 3-coordinated, H on every edge C, C-H 0.1086 nm',
          np.all(deg[el == 'C'] == 3) and np.all(deg[el == 'H'] == 1) and np.allclose(ch, 0.1086),
          f'C{fl.n_carbon}H{fl.n_atoms - fl.n_carbon}')
    flq = make_flake('steele', (2.0, 1.6), edge_charges=(-0.115, 0.115))
    check('charged-edge flake neutral', abs(flq.charges().sum()) < 1e-9 and flq.charges().min() == -0.115)
    for n, m, cells, n_at, r in ((8, 8, 2, 64, 0.5425), (10, 0, 2, 80, 0.3915), (6, 4, 1, 152, 0.3413)):
        t = make_nanotube('steele', n, m, cells)
        deg = np.bincount(t.bonds().ravel(), minlength=t.n_atoms)
        check(f'({n},{m}) nanotube: {n_at} atoms, R {r} nm, all C 3-coordinated',
              t.n_atoms == n_at and abs(t.info['radius_nm'] - r) < 2e-4 and np.all(deg == 3),
              f"R {t.info['radius_nm']:.4f}")


def _forces(system, xyz, box=None):
    ctx = mm.Context(system, mm.VerletIntegrator(0.001), mm.Platform.getPlatformByName('Reference'))
    if box is not None:
        ctx.setPeriodicBoxVectors(*[mm.Vec3(*v) for v in np.diag(box)])
    ctx.setPositions(xyz * NM)
    st = ctx.getState(getForces=True, getEnergy=True)
    return st.getForces(asNumpy=True).value_in_unit(KJ / NM), st.getPotentialEnergy().value_in_unit(KJ), ctx


def test_standalone():
    sheet = make_slab('steele', 1, (2.6, 2.5), mobility='flexible', restraint_k=0)
    s, top, xyz, info = carbon_system(sheet)
    f, e, _ = _forces(s, xyz, np.array(info['box_nm']))
    mi = info['mobility']
    check('flexible graphene: GAFF2 terms, 3 bonds/2 C, 3 angles/C, 12 torsions/2 C',
          mi['n_bonds'] == 1.5 * sheet.n_atoms and mi['n_angles'] == 3 * sheet.n_atoms and
          mi['n_torsions'] == 6 * sheet.n_atoms)
    check('flexible graphene: ideal lattice is the minimum (max |F|)', np.abs(f).max() < 1e-3,
          f'{np.abs(f).max():.2e} kJ/mol/nm')
    fl = make_flake('steele', (2.0, 1.6), edge_charges=(-0.115, 0.115))
    s, top, xyz, info = carbon_system(fl)
    f, e, _ = _forces(s, xyz)
    check('flexible flake: C-H constrained, impropers on edge C, finite forces',
          info['mobility']['n_ch_constraints'] == fl.n_atoms - fl.n_carbon ==
          info['mobility']['n_impropers'] and np.all(np.isfinite(f)), f'max |F| {np.abs(f).max():.1f}')
    s2, _, xyz2, info2 = carbon_system(fl, minimize=True)
    _, e2, ctx = _forces(s2, xyz2)
    mm.LocalEnergyMinimizer.minimize(ctx, 0.1, 20000)       # C-H constrained: compare energies
    e3 = ctx.getState(getEnergy=True).getPotentialEnergy().value_in_unit(KJ)
    check('flake minimize=True relaxes the edges', e2 < e and e2 - e3 < 0.05 * fl.n_atoms,
          f"E {e:.0f} -> {e2:.1f} (re-minimised {e3:.1f}) kJ/mol, "
          f"max shift {info2['minimized']['max_shift_nm']:.4f} nm")
    tube = make_nanotube('steele', 8, 8, 9, mobility='restrained')
    s, top, xyz, info = carbon_system(tube)
    box = np.array(info['box_nm'])
    _, e0, ctx = _forces(s, xyz, box)
    xyz2 = xyz.copy()
    xyz2[:5, 2] += box[2]          # what GPU coordinate wrapping does
    ctx.setPositions(xyz2 * NM)
    e1 = ctx.getState(getEnergy=True).getPotentialEnergy().value_in_unit(KJ)
    check('carbon restraint unchanged by a periodic image shift', abs(e1 - e0) < 1e-6, f'{e1 - e0:.2e}')
    tmp = Path(tempfile.mkdtemp(prefix='eb_carbon_'))
    for name, model in (('slab', make_slab('wu_aluru', 3, (2.6, 2.5), stacking='ABC')),
                        ('flake', fl), ('tube', tube)):
        out = write_carbon(model, tmp / name)
        files = {p.name for p in out.iterdir()}
        sysx = mm.XmlSerializer.deserialize((out / 'system.xml').read_text())
        conect = sum(1 for ln in (out / 'structure.pdb').read_text().splitlines() if ln.startswith('CONECT'))
        check(f'write_carbon({name}): files, system round-trip, CONECT records',
              {'structure.pdb', 'structure.xyz', 'forcefield.xml', 'system.xml', 'metadata.json'} <= files
              and sysx.getNumParticles() == model.n_atoms and conect > 0)


def test_steps():
    for edge in ('armchair', 'zigzag'):
        for orient in ('armchair_x', 'zigzag_x'):
            s = make_slab('steele', 3, (3.0, 3.0), orientation=orient)
            t = add_step(s, edge=edge, keep_fraction=0.5)
            el = np.array(t.elements)
            deg = np.bincount(t.bonds().ravel(), minlength=t.n_atoms)
            h = np.flatnonzero(el == 'H')
            check(f'{edge} step on {orient}: only {edge} edges, every C 3-bonded, H in the top layer',
                  set(edge_carbons(t).values()) == {edge} and np.all(deg[el == 'C'] == 3)
                  and np.all(t.layer[h] == 0) and np.all(deg[h] == 1),
                  f"{t.features[-1]['removed_carbons']} C removed, {len(h)} H")
    with raises('add_step refuses a nanotube'):
        from veloxchem.electrode.graphite import make_nanotube
        add_step(make_nanotube('steele', 8, 8, 4))


def test_defects():
    import networkx as nx
    s = make_slab('steele', 2, (4.0, 4.0))
    d = add_defects(s, vacancies=1, divacancies=1, seed=2, relax=False)
    deg = np.bincount(d.bonds().ravel(), minlength=d.n_atoms)
    check('vacancy + divacancy: 3 atoms removed, 7 dangling carbons in the top layer',
          d.n_atoms == s.n_atoms - 3 and int((deg[d.layer == 0] == 2).sum()) == 7)
    w = add_defects(s, stone_wales=1, seed=3)
    g = nx.Graph()
    g.add_edges_from((int(i), int(j)) for i, j in w.bonds() if w.layer[i] == 0 and w.layer[j] == 0)
    sizes = [len(c) for c in nx.minimum_cycle_basis(g)]
    deg = np.bincount(w.bonds().ravel(), minlength=w.n_atoms)
    check('Stone-Wales: 2 pentagons, heptagons, every C 3-bonded, relaxed',
          sizes.count(5) == 2 and sizes.count(7) >= 1 and np.all(deg == 3) and w.info['relaxed'],
          f"max shift {w.info['relaxed'][0]['max_shift_nm']:.3f} nm")
    with raises('no room for too many defects'):
        add_defects(make_slab('steele', 2, (2.0, 2.0)), vacancies=5, min_separation_nm=1.5)


def _fragment_ok(m):
    """Each group's atoms + the skeleton atoms its template charged carry the formal charge."""
    q = m.charges()
    return [(f['group'], float(q[(m.group == f['group_id']) | (m.fragment == f['group_id'])].sum()),
             f['formal_charge']) for f in m.features if f['kind'] == 'group']


def test_groups(templates):
    base = add_step(make_slab('steele', 3, (4.3, 4.3)), edge='zigzag', keep_fraction=0.5)
    m = add_groups(base, {'carboxylate': 1, 'phenolate': 1, 'phenol': 1, 'carboxyl': 1,
                          'epoxide': 1, 'hydroxyl_pair': 1}, seed=5, templates=templates)
    total = sum(SURFACE_GROUPS[g]['charge'] for g in ('carboxylate', 'phenolate'))
    frags = _fragment_ok(m)
    check('group fragments carry their formal charges, model charge = -2',
          all(abs(qs - qf) < 1e-6 for _, qs, qf in frags) and abs(m.charges().sum() - total) < 1e-6,
          ', '.join(f'{g} {qs:+.3f}' for g, qs, _ in frags))
    gz = m.xyz[~m.skeleton, 2]
    check('groups sit on the top face and were relaxed', gz.min() > m.top_z - 0.36 and
          m.info['groups_relaxed']['max_shift_nm'] < 0.2,   # OH rotation, COO tilting off the terrace
          f"group z {gz.min():.3f}..{gz.max():.3f}, top layer {m.top_z:.3f}, "
          f"relax shift {m.info['groups_relaxed']['max_shift_nm']:.3f} nm")
    s, top, xyz, info = carbon_system(m)
    f, e, _ = _forces(s, xyz, np.array(info['box_nm']))
    mi = info['mobility']
    gm = np.array([s.getParticleMass(i)._value for i in np.flatnonzero(~m.skeleton)])
    check('standalone functionalised slab: frozen skeleton, mobile groups, group terms, finite energy',
          np.all(gm > 0) and mi['group_terms'].get('surface group torsions', 0) > 0 and
          mi['group_constraints'] == 4 and np.isfinite(e), f"{mi['group_terms']}, E {e:.1f} kJ/mol")

    aq = add_groups(add_step(make_slab('steele', 3, (3.5, 3.4)), edge='armchair'),
                    {'o_quinone': 1, 'carboxylate': 1}, seed=1, face='top', templates=templates)
    check('armchair o-quinone + carboxylate placed', {f['group'] for f in aq.features if f['kind'] == 'group'}
          == {'o_quinone', 'carboxylate'})
    bot = add_groups(make_slab('steele', 3, (3.0, 3.0)), {'epoxide': 1}, face='bottom', seed=2,
                     templates=templates)
    check('bottom-face epoxide below the bottom layer', bot.xyz[~bot.skeleton, 2].max() < -bot.top_z - 0.1)

    vac = add_defects(make_slab('steele', 2, (3.0, 3.0)), vacancies=1, seed=4, relax=False)
    nbrs = _neighbours(vac)
    dangling = [i for i in range(vac.n_atoms) if vac.layer[i] == 0 and len(nbrs[i]) == 2]
    j = next(k for k in nbrs[dangling[0]])
    origin, frame = _site_frame(vac, (dangling[0], j), 'basal', 'top', nbrs)
    from veloxchem.electrode.graphite.surface import _with_group_arrays
    check('template refuses a site whose environment differs (next to a vacancy)',
          _apply_template(_with_group_arrays(vac.copy()), templates['epoxide/basal'],
                          (dangling[0], j), origin, frame, nbrs, 'top') is None)
    return m


def templates_for_tests():
    t = load_templates()
    if t:
        return t
    print('using the STO-3G test templates (tests/data); production: data/surface_groups.json')
    return load_templates(TEST_TEMPLATE_FILE)


def test_builder_class(templates):
    gb = (GraphiteBuilder('steele', 3, (3.5, 3.4)).add_step(edge='zigzag')
          .add_groups({'carboxylate': 1, 'epoxide': 1}, seed=2, templates=templates))
    r = gb.roles()
    mol = gb.molecule()
    meta = gb.metadata()
    check('GraphiteBuilder: chained step + groups, roles, molecule charge, metadata',
          len(r['carbon']) == sum(1 for e, k in zip(gb.elements, gb.model.skeleton) if k and e == 'C')
          and len(r['surface_groups']) == 4 and r['carbon_h']
          and mol.get_charge() == -1 and gb.charge == -1 and meta['groups'] == {'carboxylate': 1, 'epoxide': 1}
          and gb.has_surface_features, f"{gb.n_atoms} atoms, groups {meta['groups']}")
    system, top, xyz, info = gb.system()
    check('GraphiteBuilder.system(): standalone carbon system with the group terms',
          system.getNumParticles() == gb.n_atoms and info['mobility']['group_constraints'] == 0)


if __name__ == '__main__':
    test_refusals()
    test_geometry()
    test_standalone()
    test_nbfix()
    test_steps()
    test_defects()
    templates = templates_for_tests()
    test_groups(templates)
    test_builder_class(templates)
    print('all graphite checks passed')
