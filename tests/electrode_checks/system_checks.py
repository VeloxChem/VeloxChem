"""Self-checks for the system layer on top of the graphite (CPU). Run:
python tests/test_system.py [--fep] [--require-fep]   (--fep: tiny end-to-end FEP, ~10 min;
    --require-fep: a skipped FEP-driver check is a failure. Production driver:
    ELECTRODE_FEP_DRIVER=../../third_party/veloxchem/solvationfepdriver.py)
"""
import json
import math
import sys
import tempfile
from pathlib import Path

import numpy as np
import openmm as mm
from openmm import unit

HERE = Path(__file__).resolve().parent

from veloxchem.electrode import build_system  # noqa: E402
from veloxchem.environment import get_data_path  # noqa: E402
from veloxchem.electrode.system_tools import (FIELD_FORCE_NAME, GraphiteBuilder, build, build_bulk,  # noqa: E402
                          electrode_fep_driver, fep_window_systems, h2so4_mixture, ion_correction_settings,
                          koh_mixture, make_simulation, nonbonded_force, slot_centres,
                          place_solute, run_fep, write_fep_inputs)
from veloxchem.electrode.graphite import make_slab  # noqa: E402
from veloxchem.electrode.graphite.params import CARBON_SETS  # noqa: E402
from veloxchem.electrode.graphite.surface import load_templates  # noqa: E402
from veloxchem.electrode.graphite.surface_templates import TEST_TEMPLATE_FILE  # noqa: E402

SOLUTE = get_data_path() / 'electrode' / 'examples' / '14a'
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


def small_system(**kw):
    gb = GraphiteBuilder(kw.pop('carbon_set', 'steele'), layers=2, lateral_nm=(2.6, 2.5),
                         mobility=kw.pop('mobility', 'frozen'),
                         allow_unverified=kw.get('allow_unverified', False))
    sp, prop = koh_mixture(1.0)
    return build(gb, SOLUTE / 'solute.pdb', SOLUTE / 'solute.xml', sp, prop,
                 box_z_nm=4.5, seed=3, **kw)


def test_default_water():
    slab = make_slab('steele', layers=2, lateral_nm=(2.6, 2.5))
    es = build(slab, SOLUTE / 'solute.pdb', SOLUTE / 'solute.xml', box_z_nm=4.5, seed=3)
    names = {r.name for r in es.topology.residues()}
    nb = nonbonded_force(es.system)
    q_sol = round(sum(nb.getParticleParameters(i)[0]._value for i in es.indices['solute']))
    n_k = sum(1 for r in es.topology.residues() if r.name == 'K')
    check('default electrolyte is pure water; K+ only to neutralise the solute',
          names == {'MOL', 'GRA', 'HOH', 'K'} and n_k == -q_sol and
          list(es.metadata['species']) == ['water_spce', 'K_jc_spce'] and abs(es.metadata['total_charge']) < 1e-6,
          f'residues {sorted(names)}, K+ {n_k}, solute {q_sol:+d}')
    bulk = build_bulk(box_nm=2.0, seed=4)
    check('build_bulk default: pure water, no ions', {r.name for r in bulk.topology.residues()} == {'HOH'})
    with raises('species without proportions refused'):
        build(slab, SOLUTE / 'solute.pdb', SOLUTE / 'solute.xml', ['water_spce'], box_z_nm=4.5)


def test_build():
    es = small_system()
    s, idx, m = es.system, es.indices, es.metadata
    nb = nonbonded_force(s)
    check('neutral', abs(m['total_charge']) < 1e-6)
    first = next(es.topology.residues())
    check('solute is residue 0 (FEP alchemical region)', first.name == 'MOL' and
          [a.index for a in first.atoms()] == idx['solute'])
    q, sig, eps = nb.getParticleParameters(idx['carbon'][0])
    cs = CARBON_SETS['steele']
    check('graphite LJ in NonbondedForce = steele, q = 0',
          q._value == 0 and abs(sig._value - cs['sigma']) < 1e-12 and abs(eps._value - cs['epsilon']) < 1e-12)
    check('graphite frozen', all(s.getParticleMass(i)._value == 0 for i in idx['carbon']))
    o, h = es.indices['anion_center'][0], es.indices['anion_center'][0] + 1
    qo, so, eo = (v._value for v in nb.getParticleParameters(o))
    qh, _, eh = (v._value for v in nb.getParticleParameters(h))
    check('default OH- = Bonthuis 2016 (O -1, 0.381 nm, 0.05 kJ/mol; H q 0, no LJ)',
          (qo, round(so, 6), eo, qh, eh) == (-1.0, 0.381, 0.05, 0.0, 0.0))
    check('K+/OH- cross-term warning recorded', any('Bonthuis 2016 OH- with cations' in w for w in m['warnings']))
    gap = m['graphite']['nearest_solvent_heavy_atom_from_face_nm']
    check('solvent heavy atoms >= 0.25 nm from the faces', gap >= 0.249, f'{gap:.3f} nm')
    check('no close contacts reported', not any('closer than' in w for w in m['warnings']))
    n_water = len(idx['water_O'])
    check('water rigid (3 constraints each)', s.getNumConstraints() >= 3 * n_water,
          f'{s.getNumConstraints()} constraints, {n_water} waters')
    field = [f for f in s.getForces() if f.getName() == FIELD_FORCE_NAME][0]
    n_charged = sum(1 for i in range(s.getNumParticles()) if nb.getParticleParameters(i)[0]._value != 0)
    check('field covers every charged particle', field.getNumParticles() == n_charged)
    masses = np.array([s.getParticleMass(i)._value for i in idx['solute']])
    com = (es.positions[idx['solute']] * masses[:, None]).sum(0) / masses.sum()
    h = com[2] - es.positions[idx['top_layer'], 2].mean()
    check('solute COM 0.40 nm above the top layer', abs(h - 0.40) < 2e-3, f'{h:.4f}')
    mob = m['graphite']['mobility']
    check('frozen: close C-C pairs excluded', mob['mode'] == 'frozen' and
          mob['excluded_pairs'] > 5 * len(idx['carbon']),
          f"{mob['excluded_pairs']} pairs for {len(idx['carbon'])} C")
    sim = make_simulation(es, platform='CPU')
    st = sim.context.getState(getPositions=True, getEnergy=True)
    e_min = st.getPotentialEnergy().value_in_unit(KJ)
    c0 = st.getPositions(asNumpy=True).value_in_unit(NM)[idx['carbon']]
    built = es.positions[idx['carbon']]
    check('minimisation kept graphite fixed', np.abs(c0 - built).max() < 1e-6)
    sim.step(500)
    st = sim.context.getState(getPositions=True, getEnergy=True)
    e = st.getPotentialEnergy().value_in_unit(KJ)
    moved = np.abs(st.getPositions(asNumpy=True).value_in_unit(NM)[idx['carbon']] - c0).max()
    check('1 ps MD stable', np.isfinite(e) and e < 0 and e - e_min < 0.1 * abs(e_min),
          f'E {e_min:.0f} -> {e:.0f} kJ/mol')
    check('graphite did not move', moved == 0.0)
    with raises('barostat refused for slab systems'):
        make_simulation(es, barostat=True, platform='CPU', minimize=False)
    print('   ', m['liquid'], m['packing'])
    return es


def test_ion_pair():
    from veloxchem.electrode.system_tools import add_ion_pair_sigma_scaling
    s = mm.System()
    for _ in range(2):
        s.addParticle(20.0)
    s.setDefaultPeriodicBoxVectors(mm.Vec3(5, 0, 0), mm.Vec3(0, 5, 0), mm.Vec3(0, 0, 5))
    nb = mm.NonbondedForce()
    nb.setNonbondedMethod(mm.NonbondedForce.CutoffPeriodic)
    nb.setCutoffDistance(2.0)
    nb.setUseDispersionCorrection(False)
    nb.setReactionFieldDielectric(1.0)
    sk, ek, so, eo = 0.2838403315995121, 1.7978873936000002, 0.381, 0.05   # JC K+, Bonthuis OH- O
    nb.addParticle(0.0, sk, ek)
    nb.addParticle(0.0, so, eo)
    s.addForce(nb)
    add_ion_pair_sigma_scaling(s, [0], [1], 1.5, cutoff=2.0, switch=None)
    ctx = mm.Context(s, mm.VerletIntegrator(0.001), mm.Platform.getPlatformByName('Reference'))
    sg, eg = 1.5 * math.sqrt(sk * so), math.sqrt(ek * eo)
    for r in (0.45, 0.5, 0.7):
        ctx.setPositions([mm.Vec3(1, 1, 1), mm.Vec3(1, 1, 1 + r)])
        e = ctx.getState(getEnergy=True).getPotentialEnergy().value_in_unit(KJ)
        ref = 4 * eg * ((sg / r) ** 12 - (sg / r) ** 6)
        check(f'cation-anion LJ = Bonthuis Eq. 22 at {r} nm', math.isclose(e, ref, rel_tol=1e-6),
              f'{e:.5f} vs {ref:.5f}')


def test_wu_aluru_build():
    es = small_system(carbon_set='wu_aluru')
    names = [f.getName() for f in es.system.getForces()]
    check('wu_aluru build adds the C-water NBFIX force', 'graphite-water NBFIX correction' in names)


# ---- FEP-driver checks and skips (review 2026-11-18, v2)
# A check may be skipped ONLY for a known, deliberately unsupported driver: the GSC guard
# (unsafe placeholder sigma) or a driver whose internals differ. Any other RuntimeError is a
# failure. --require-fep turns every skip into a failure (use it with the production driver:
# ELECTRODE_FEP_DRIVER=third_party/veloxchem/solvationfepdriver.py). Skips are counted and
# reported at the end, with the driver file and SHA-256.
SKIPPED = []
_UNSUPPORTED = ('placeholder sigma reached the soft-core term', 'driver internals differ')


def fep_windows_or_skip(d, lambdas, label):
    try:
        return fep_window_systems(d, lambdas)
    except RuntimeError as e:
        if not any(t in str(e) for t in _UNSUPPORTED):
            raise
        if '--require-fep' in sys.argv:
            raise AssertionError(f'{label}: FEP-driver check could not run ({e})') from e
        SKIPPED.append(f'{label}: {e}')
        print(f'SKIP {label}: {e}')
        return None


def test_fep_driver(es):
    tmp = Path(tempfile.mkdtemp(prefix='eb_fep_'))
    write_fep_inputs(es, tmp)
    d = electrode_fep_driver(es.indices['solute'], es.indices['carbon'])
    d.solute_pdb, d.solute_xml = str(tmp / 'solute.pdb'), str(tmp / 'solute.xml')
    d.system_pdb, d.system_xml = str(tmp / 'system.pdb'), str(tmp / 'system.xml')
    d.solvent_name = 'omm_system_xml'
    d.stage = 1
    res = fep_windows_or_skip(d, [0.0, 0.5, 1.0], 'FEP-driver unit checks')
    if res is None:
        return tmp
    systems, top, _ = res
    s_atom, w_atom = es.indices['solute'][0], es.indices['water_O'][0]

    def field_q(system, p):
        f = [f for f in system.getForces() if f.getName() == FIELD_FORCE_NAME][0]
        for k in range(f.getNumParticles()):
            i, prm = f.getParticleParameters(k)
            if i == p:
                return prm[0]
        return 0.0

    q0 = field_q(systems[0], s_atom)
    qs = [field_q(sy, s_atom) for sy in systems]
    check('solute field charge scaled 1, 0.5, 0 with lambda', np.allclose(qs, [q0, 0.5 * q0, 0.0]),
          str(np.round(qs, 4)))
    check('environment field charge unchanged', len({round(field_q(sy, w_atom), 8) for sy in systems}) == 1)
    # the GSC placeholder guard must have inspected a soft-core force in every window of the
    # stages that use one (stage 1 and 2); a count of 0 means the guard checked nothing
    d.stage = 2
    if fep_windows_or_skip(d, [1.0, 0.5, 0.0], 'FEP-driver stage-2 windows') is not None:
        g = getattr(d, 'electrode_gsc_checked', {})
        check('GSC guard inspected a soft-core force in every stage-1 and stage-2 window',
              len(g.get(1, [])) == 3 and len(g.get(2, [])) == 3 and all(n >= 1 for st in (1, 2) for n in g[st]),
              f'forces per window: stage 1 {g.get(1)}, stage 2 {g.get(2)}')
    d.stage = 3
    vac, vtop, _ = d._create_vacuum_systems([0.0])
    check('vacuum legs build from solute.xml', vac[0].getNumParticles() == len(es.indices['solute']))
    return tmp


def test_mobile(mode):
    es = small_system(mobility=mode)
    s, idx, m = es.system, es.indices, es.metadata
    names = [f.getName() for f in s.getForces()]
    check(f'{mode}: carbon has mass and a periodic restraint',
          all(s.getParticleMass(i)._value > 0 for i in idx['carbon']) and
          'carbon position restraint' in names and m['graphite']['mobility']['mode'] == mode)
    if mode == 'flexible':
        check('flexible: carbon bonded forces use PBC', all(
            f.usesPeriodicBoundaryConditions() for f in s.getForces() if f.getName().startswith('carbon ')))
    sim = make_simulation(es, platform='CPU')
    c0 = sim.context.getState(getPositions=True).getPositions(asNumpy=True).value_in_unit(NM)[idx['carbon']]
    sim.step(500)
    st = sim.context.getState(getPositions=True, getEnergy=True)
    c1 = st.getPositions(asNumpy=True).value_in_unit(NM)[idx['carbon']]
    rms = float(np.sqrt(((c1 - c0) ** 2).sum(1).mean()))
    top_z = c1[np.isin(idx['carbon'], idx['top_layer']), 2].mean()
    thermal = math.sqrt(3 * 2.479 / 1000.0)      # 3D rms of a k = 1000 tether at 298 K
    check(f'{mode}: 1 ps MD stable, carbons move but stay put', np.isfinite(
          st.getPotentialEnergy().value_in_unit(KJ)) and 0 < rms < 1.5 * thermal and
          abs(top_z - m['graphite']['top_face_z']) < 0.03,
          f'rms {rms:.3f} nm (k=1000 tether alone: {thermal:.3f}), top layer z {top_z:.3f}')
    tmp = Path(tempfile.mkdtemp(prefix=f'eb_{mode}_'))
    write_fep_inputs(es, tmp)
    d = electrode_fep_driver(idx['solute'], idx['carbon'], mode)
    d.solute_pdb, d.solute_xml = str(tmp / 'solute.pdb'), str(tmp / 'solute.xml')
    d.system_pdb, d.system_xml = str(tmp / 'system.pdb'), str(tmp / 'system.xml')
    d.solvent_name, d.stage = 'omm_system_xml', 1
    res = fep_windows_or_skip(d, [0.0, 1.0], f'{mode}: FEP-driver window check')
    if res is None:
        return
    systems = res[0]
    check(f'{mode}: FEP driver accepts the system, carbon keeps its mass',
          all(sy.getParticleMass(idx['carbon'][0])._value > 0 for sy in systems))


def test_bulk():
    sp, prop = koh_mixture(1.0)
    es = build_bulk(sp, prop, box_nm=2.5, seed=4)
    check('bulk box neutral', abs(es.metadata['total_charge']) < 1e-6)
    sim = make_simulation(es, barostat=True, platform='CPU')
    sim.step(200)
    v = sim.context.getState().getPeriodicBoxVolume().value_in_unit(NM ** 3)
    check('bulk NPT runs', np.isfinite(v), f'V = {v:.2f} nm3')


def test_h2so4():
    """3.4 mol/kg H2SO4 (= 3 M) as H3O+ (Bonthuis 2016, rigid) + HSO4- (GAFF2/RESP); the -3
    solute is neutralised with H3O+, the acid's own cation."""
    gb = GraphiteBuilder('wu_aluru', layers=1, lateral_nm=(2.6, 2.5))
    sp, prop = h2so4_mixture(3.4)
    check('H2SO4 mixture: water 55.51 - m per m H3O+ and HSO4-',
          sp == ['water_spce', 'H3O_bonthuis2016', 'HSO4_gaff2_resp'] and abs(prop[0] - 52.108) < 1e-9)
    es = build(gb, SOLUTE / 'solute.pdb', SOLUTE / 'solute.xml', sp, prop, box_z_nm=4.5, seed=3,
               liquid_density=1178.0)
    s, idx, m = es.system, es.indices, es.metadata
    c = m['liquid']['species_counts']
    n3, n4, nw = c['H3O_bonthuis2016'], c['HSO4_gaff2_resp'], c['water_spce']
    check('neutral with the acid ions only', abs(m['total_charge']) < 1e-6 and n3 == n4 + 3
          and m['liquid']['counterion'] == 'H3O_bonthuis2016', f'H3O {n3}, HSO4 {n4}, water {nw}')
    check('about 15.3 water per acid', abs(nw / n4 - 52.108 / 3.4) < 2.0, f'{nw / n4:.1f}')
    nb = nonbonded_force(s)
    o = idx['cation_center'][0]
    qo, so, eo = (v._value for v in nb.getParticleParameters(o))
    qh = nb.getParticleParameters(o + 1)[0]._value
    check('H3O+ = Bonthuis 2016 (O -1.4, 0.31 nm, 0.8 kJ/mol; H +0.8)',
          (round(qo, 6), round(so, 6), eo, round(qh, 6)) == (-1.4, 0.31, 0.8, 0.8))
    cons = {}
    for k in range(s.getNumConstraints()):
        i, j, d = s.getConstraintParameters(k)
        cons[frozenset((i, j))] = d.value_in_unit(NM)
    rigid = all(abs(cons.get(frozenset((a, b)), 0) - r) < 1e-9
                for h0 in idx['cation_center'] for a, b, r in
                [(h0, h0 + 1, 0.098), (h0, h0 + 2, 0.098), (h0, h0 + 3, 0.098),
                 (h0 + 1, h0 + 2, 0.1619), (h0 + 1, h0 + 3, 0.1619), (h0 + 2, h0 + 3, 0.1619)])
    check('H3O+ rigid: O-H 0.098 and H-H 0.1619 nm constrained', rigid)
    hso = [r for r in es.topology.residues() if r.name == 'HSO']
    q4 = sum(nb.getParticleParameters(a.index)[0]._value for a in hso[0].atoms())
    check('HSO4- from the GAFF2/RESP file, charge -1, S is the anion centre',
          abs(q4 + 1) < 1e-4 and all(next(r.atoms()).element.symbol == 'S' for r in hso)
          and set(idx['anion_center']) == {next(r.atoms()).index for r in hso})
    check('slots = molecules (water + H3O+ + HSO4-)',
          len(slot_centres(es)) == m['liquid']['n_slots'] == nw + n3 + n4)
    sim = make_simulation(es, platform='CPU')
    sim.step(200)
    x = sim.context.getState(getPositions=True).getPositions(asNumpy=True).value_in_unit(NM)
    h0 = idx['cation_center'][0]
    check('MD runs, H3O+ stays rigid', np.isfinite(x).all()
          and abs(np.linalg.norm(x[h0 + 1] - x[h0 + 2]) - 0.1619) < 1e-4)

    class D:
        ion_corrections = 'auto'
    d = D()
    _, why = ion_correction_settings(d, dict(graphite=dict(layers=0), species={'water_spce': {}},
                                             liquid=dict(electrolyte_species=sp)), {})
    check("bulk reference in the acid: ion corrections 'none'", d.ion_corrections == 'none', why[:40])
    d = D()
    kw, _ = ion_correction_settings(d, dict(graphite=dict(layers=0), species={'water_spce': {}},
                                            liquid=dict(electrolyte_species=['water_spce'])), {})
    check("bulk reference in pure water: still 'auto' spce", d.ion_corrections == 'auto'
          and kw == {'water_model': 'spce'})


def test_tiny_fep(folder):
    res = run_fep(folder, num_steps=200, num_snapshots=10, num_equil_steps=50, platform='CPU')
    check('tiny end-to-end FEP finishes', np.isfinite(res['free_energy']),
          f"{res['free_energy']:.1f} kJ/mol (not converged; plumbing test only)")
    info = res.get('fep_driver') or {}
    check('FEP result records the driver file and SHA-256', len(info.get('sha256', '')) == 64,
          f"{info.get('source_file')} {info.get('sha256')}")
    import veloxchem as vlx  # noqa: F401
    d = electrode_fep_driver()
    if hasattr(d, 'lambdas_stage1'):
        check('driver: old 3-stage lambda schedule and GSC placeholder sigma 0.1 nm',
              list(d.lambdas_stage1)[:3] == [0.0, 0.05, 0.1] and abs(d.gsc_placeholder_sigma.value_in_unit(unit.nanometer) - 0.1) < 1e-12,
              f'{d.lambdas_stage1}, {d.gsc_placeholder_sigma}')


def test_build_functionalised(m):
    sp, prop = koh_mixture(1.0)
    xy = _free_spot(m)
    es = build(m, SOLUTE / 'solute.pdb', SOLUTE / 'solute.xml', sp, prop, box_z_nm=4.5, seed=3,
               solute_xy_nm=xy, height_nm=0.40)
    s, idx, meta = es.system, es.indices, es.metadata
    nb = nonbonded_force(s)
    n_k = sum(1 for r in es.topology.residues() if r.name == 'K')
    n_oh = sum(1 for r in es.topology.residues() if r.name == 'OHM')
    q_sol = round(sum(nb.getParticleParameters(i)[0]._value for i in idx['solute']))
    q_surf = round(sum(nb.getParticleParameters(i)[0]._value for i in idx['carbon_model']))
    check('functionalised build neutral, extra K+ = -(solute + surface charge)',
          abs(meta['total_charge']) < 1e-6 and q_surf == -2 and n_k - n_oh == -(q_sol + q_surf),
          f'K+ {n_k}, OH- {n_oh}, solute {q_sol:+d}, surface {q_surf:+d}')
    check('skeleton frozen, groups mobile, NBFIX/field bookkeeping',
          all(s.getParticleMass(i)._value == 0 for i in idx['carbon'] + idx['carbon_h']) and
          all(s.getParticleMass(i)._value > 0 for i in idx['surface_groups']))
    field = [f for f in s.getForces() if f.getName() == FIELD_FORCE_NAME][0]
    n_charged = sum(1 for i in range(s.getNumParticles()) if nb.getParticleParameters(i)[0]._value != 0)
    check('field covers every charged particle incl. the groups', field.getNumParticles() == n_charged)
    dmin = meta['graphite']['nearest_solvent_heavy_atom_to_carbon_model_atom_nm']
    check('packed solvent keeps >= 0.25 nm from every surface atom, water on the lower terrace',
          dmin >= 0.249 and _water_below_top_layer(es), f'{dmin:.3f} nm')
    sim = make_simulation(es, platform='CPU')
    sim.step(500)
    st = sim.context.getState(getPositions=True, getEnergy=True)
    x = st.getPositions(asNumpy=True).value_in_unit(NM)
    stretch = _max_group_bond_dev(s, x, es)
    check('1 ps MD: stable, group bonds intact', np.isfinite(st.getPotentialEnergy().value_in_unit(KJ))
          and stretch < 0.02, f'max group bond deviation {stretch:.4f} nm')
    tmp = Path(tempfile.mkdtemp(prefix='eb_surf_'))
    write_fep_inputs(es, tmp)
    d = electrode_fep_driver(idx['solute'], idx['carbon'], 'frozen')
    d.solute_pdb, d.solute_xml = str(tmp / 'solute.pdb'), str(tmp / 'solute.xml')
    d.system_pdb, d.system_xml = str(tmp / 'system.pdb'), str(tmp / 'system.xml')
    d.solvent_name, d.stage = 'omm_system_xml', 1
    res = fep_windows_or_skip(d, [0.0, 1.0], 'functionalised: FEP-driver window check')
    if res is None:
        return es
    systems = res[0]
    check('FEP driver accepts the functionalised system', len(systems) == 2)
    return es


def _free_spot(m, height=0.40):
    """A solute xy on the top terrace that clears every group atom and edge H."""
    m = getattr(m, 'model', m)
    from openmm import app
    pdb = app.PDBFile(str(SOLUTE / 'solute.pdb'))
    xyz = np.array(pdb.positions.value_in_unit(NM))
    mass = np.array([a.element.mass._value for a in pdb.topology.atoms()])
    lx, ly = m.box_xy
    best = None
    for x in np.linspace(0, lx, 15, endpoint=False):
        for y in np.linspace(0, ly, 15, endpoint=False):
            try:
                out = place_solute(xyz, mass, m, height, xy=(x, y), min_contact_nm=0.30)
            except ValueError:
                continue
            under = [i for i in range(m.n_atoms) if m.layer[i] == 0 and m.skeleton[i]
                     and np.linalg.norm(((m.xyz[i] - out.mean(0))[:2] + 0.5 * np.array([lx, ly]))
                                        % np.array([lx, ly]) - 0.5 * np.array([lx, ly])) < 0.3]
            if under and (best is None or len(under) > best[0]):
                best = (len(under), (float(x), float(y)))
    if best is None:
        raise RuntimeError('no free spot for the solute')
    return best[1]


def _water_below_top_layer(es):
    """Some water O sits lower than the top layer plane (on the lower terrace of the step)."""
    top = es.metadata['graphite']['top_face_z']
    z = es.positions[es.indices['water_O'], 2]
    return bool(((z > top - 0.335 + 0.2) & (z < top)).any())


def _max_group_bond_dev(system, x, es):
    dev = 0.0
    for f in system.getForces():
        if f.getName() == 'surface group bonds':
            for k in range(f.getNumBonds()):
                i, j, r0, _ = f.getBondParameters(k)
                dv = x[j] - x[i]
                box = np.array(es.metadata['box_nm'])
                dv -= np.round(dv / box) * box
                dev = max(dev, abs(np.linalg.norm(dv) - r0.value_in_unit(NM)))
    return dev


def test_adsorption_legs():
    import types
    from veloxchem.electrode.system_tools import (adsorption_free_energy, build_reference, ion_correction_settings,
                              window_length, KB_KJ_MOL_K)
    # window length vs numerical integration
    lo, hi, k, temp = 0.0, 0.8, 2000.0, 298.15
    h = np.linspace(-1.0, 2.0, 300001)
    u = 0.5 * k * (np.maximum(0, h - hi) ** 2 + np.minimum(0, h - lo) ** 2)
    num = float(np.trapezoid(np.exp(-u / (KB_KJ_MOL_K * temp)), h))
    check('window length = flat part + harmonic walls', math.isclose(window_length(lo, hi, k, temp), num, rel_tol=1e-6),
          f'{window_length(lo, hi, k, temp):.6f} vs {num:.6f} nm')
    # both legs through the working script
    legs = {}
    for leg in ('surface', 'mid'):
        out = Path(tempfile.mkdtemp(prefix=f'eb_{leg}_'))
        cfg = build_system.load_config([])
        cfg.update(output=str(out), box_z_nm=4.5, field_V_per_nm=0.2,
                   graphite=dict(carbon_set='wu_aluru', layers=2, lateral_nm=[2.6, 2.5]),
                   restraint={'leg': leg}, md=dict(minimize=True, steps=300, platform='CPU'))
        es = build_system.main(cfg)
        legs[leg] = (out, es)
        r, md = es.metadata['restraint'], es.metadata['pre_fep_md']['solute_height_nm']
        f = [f for f in es.system.getForces() if f.getName() == 'solute height window'][0]
        check(f'{leg} leg: window force on, solute stayed inside during MD',
              f.getForceGroup() == 5 and r['leg'] == leg and md['fraction_outside_window'] == 0.0,
              f"window {r['lo_nm']:.2f}..{r['hi_nm']:.2f}, h {md['min']:.2f}..{md['max']:.2f} nm")
    # window energy vs analytic, by moving the solute to h = hi + 0.1
    out, es = legs['surface']
    sys_ = mm.XmlSerializer.deserialize(mm.XmlSerializer.serialize(es.system))
    ctx = mm.Context(sys_, mm.VerletIntegrator(0.001), mm.Platform.getPlatformByName('Reference'))
    idx, r = es.indices, es.metadata['restraint']
    m = np.array([sys_.getParticleMass(i)._value for i in idx['solute']])
    x = es.positions.copy()
    hcur = (x[idx['solute'], 2] * m).sum() / m.sum() - x[idx['top_layer'], 2].mean()
    x[idx['solute'], 2] += r['hi_nm'] + 0.1 - hcur
    ctx.setPositions(x * NM)
    e = ctx.getState(getEnergy=True, groups={5}).getPotentialEnergy().value_in_unit(KJ)
    check('window energy 0.5*k*d^2 just outside', math.isclose(e, 0.5 * r['k_kj_mol_nm2'] * 0.01, rel_tol=1e-6),
          f'{e:.4f} kJ/mol')
    # adsorption formula on synthetic FEP results
    for leg, value in (('surface', -1500.0), ('mid', -1480.0)):
        (legs[leg][0] / 'fep_results.json').write_text(json.dumps(
            dict(free_energy=value, uncertainty=1.0, electrode_ion_corrections='none (slab system)')))
    a = adsorption_free_energy(legs['surface'][0], legs['mid'][0])
    kt = KB_KJ_MOL_K * 298.15
    l_s = window_length(0.0, 0.8, 2000.0)
    check('K_ads = l_s exp(-(F_s - F_m)/kT); dG(delta = l_s) = F_s - F_m',
          math.isclose(a['K_ads_nm'], l_s * math.exp(20.0 / kt), rel_tol=1e-9)
          and math.isclose(a['dG_ads_window_reference_kj_mol'], -20.0)
          and math.isclose(a['dG_ads_kj_mol'], -kt * math.log(a['K_ads_nm'] / 1.0)),
          f"K {a['K_ads_nm']:.3g} nm, dG(1 nm) {a['dG_ads_kj_mol']:.2f} kJ/mol")
    with raises('legs with different restraints refused'):
        adsorption_free_energy(legs['mid'][0], legs['surface'][0])
    # bulk reference: cubic box, no graphite, no field
    ref = build_reference(SOLUTE / 'solute.pdb', SOLUTE / 'solute.xml', box_nm=3.2, seed=2)
    b = ref.metadata
    check('bulk reference: cubic, neutral, no graphite, field 0',
          len(set(np.round(b['box_nm'], 6))) == 1 and abs(b['total_charge']) < 1e-6
          and b['graphite']['layers'] == 0 and b['field']['V_per_nm'] == 0.0, f"box {b['box_nm']}")
    # ion-correction choice (stand-in drivers)
    new = types.SimpleNamespace(ion_corrections='auto')
    kw, d1 = ion_correction_settings(new, es.metadata, {})
    ok_slab = new.ion_corrections == 'none' and kw == {}
    kw2, d2 = ion_correction_settings(types.SimpleNamespace(ion_corrections='auto'), b, {})
    _, d3 = ion_correction_settings(types.SimpleNamespace(), es.metadata, {})
    check("ion corrections: slab 'none', bulk SPC/E 'auto' + water_model, old build 'unavailable'",
          ok_slab and kw2 == {'water_model': 'spce'} and d3.startswith('unavailable'), f'{d1} | {d2}')


def test_refuses_nanotube():
    with raises('build() refuses a nanotube'):
        build(GraphiteBuilder.nanotube('steele', 8, 8, 4), SOLUTE / 'solute.pdb', SOLUTE / 'solute.xml',
              box_z_nm=4.5)


def test_working_script():
    out = Path(tempfile.mkdtemp(prefix='eb_script_'))
    cfg = build_system.load_config([])
    cfg.update(output=str(out), box_z_nm=4.5, graphite=dict(carbon_set='wu_aluru', layers=2,
                                                            lateral_nm=[2.6, 2.5]),
               md=dict(minimize=True, steps=50, platform='CPU'))
    es = build_system.main(cfg)
    files = {p.name for p in out.iterdir()}
    sysx = mm.XmlSerializer.deserialize((out / 'system.xml').read_text())
    check('build_system.py: graphite -> water -> system -> MD -> FEP files',
          {'system.xml', 'system.pdb', 'solute.pdb', 'solute.xml', 'metadata.json', 'config.json'} <= files
          and sysx.getNumParticles() == es.system.getNumParticles()
          and set(es.metadata['species']) == {'water_spce', 'K_jc_spce'}
          and es.metadata['pre_fep_md']['steps'] == 50, f'{sysx.getNumParticles()} particles')


def functionalised_graphite():
    t = load_templates() or load_templates(TEST_TEMPLATE_FILE)
    return (GraphiteBuilder('steele', 3, (4.3, 4.3)).add_step(edge='zigzag', keep_fraction=0.5)
            .add_groups({'carboxylate': 1, 'phenolate': 1, 'phenol': 1, 'carboxyl': 1,
                         'epoxide': 1, 'hydroxyl_pair': 1}, seed=5, templates=t))


if __name__ == '__main__':
    test_default_water()
    test_refuses_nanotube()
    es = test_build()
    test_ion_pair()
    test_wu_aluru_build()
    folder = test_fep_driver(es)
    test_mobile('restrained')
    test_mobile('flexible')
    test_bulk()
    test_h2so4()
    test_build_functionalised(functionalised_graphite())
    test_working_script()
    test_adsorption_legs()
    if '--fep' in sys.argv:
        test_tiny_fep(folder)
    from veloxchem.electrode.system_tools import fep_driver_base
    _, src, sha = fep_driver_base()
    print(f'FEP driver: {src}\n  SHA-256 {sha}')
    if SKIPPED:
        print(f'all executed system checks passed, BUT {len(SKIPPED)} FEP-driver check group(s) SKIPPED:')
        for x in SKIPPED:
            print(f'  SKIPPED {x}')
    else:
        print('all system checks passed (no skips)')
