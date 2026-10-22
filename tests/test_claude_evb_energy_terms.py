"""WP-C3: the PES systems against an independent force-field evaluator.

FFEnergyReference (test_claude_evb_reference.py) is written from the
textbook formulas in the WP-C3 plan, section 4, and reads the state
force-field JSONs. Each PES system that build_systems writes
(run/reactant_sys.xml, run/product_sys.xml) is evaluated per force on the
Reference platform and compared term by term.

Tolerance (plan section 7): 1e-6 relative, or 1e-6 kJ/mol absolute for terms
below 1 kJ/mol; pair parameters to 1e-12 relative.
"""

import copy
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from veloxchem.evbdriver import EvbDriver
from veloxchem.molecule import Molecule
from veloxchem.mmforcefieldgenerator import MMForceFieldGenerator
from veloxchem.outputstream import OutputStream

from test_evb_helper import EvbTestHelper
from test_claude_evb_reference import (FF_TERMS, FFEnergyReference,
                                      coulomb_pair, lj_pair)

mm = pytest.importorskip('openmm')

DATA = Path(__file__).parent / 'data'

# fixture: one (force-field json, xyz) per state, in path order
FIXTURES = {
    'ethanol': [('evb_ethanol_ff_data.json', 'evb_ethanol.xyz'),
                ('evb_ethene_H2O_ff_data.json', 'evb_ethene_H2O.xyz')],
    'clch3cl': [('evb_clch3cl_reactant_ff_data.json',
                 'evb_clch3cl_reactant.xyz'),
                ('evb_clch3cl_product_ff_data.json',
                 'evb_clch3cl_product.xyz')],
    # three states, five path-global reactive pairs
    'petase_ms': [(f'evb_petase_ms_state_{k}_ff_data.json',
                   f'evb_petase_ms_state_{k}.xyz') for k in range(3)],
}

# A field large enough to give energies well above the tolerance (V/m)
E_FIELD = [2.0e8, -1.0e8, 5.0e8]

# Force name in the PES system -> evaluator term. Forces not listed must be
# empty (boundary terms, constraints) or carry no energy (CMMotionRemover).
NAME_TO_TERM = {
    'Static reaction harmonic bond': 'bond',
    'Reaction morse bond': 'morse',
    'Reaction angle': 'angle',
    'Reaction proper fourier torsion': 'torsion',
    'Reaction improper fourier torsion': 'improper',
    'Reaction internal Coul static': 'coul_static',
    'Reaction internal LJ static': 'lj_static',
    'Reaction internal Coul dynamic': 'coul_dynamic',
    'Reaction internal LJ dynamic': 'lj_dynamic',
    'NonbondedForce': 'nonbonded_force',
    'Electric field': 'efield',
}
ZERO_FORCES = {
    'CMMotionRemover', 'Boundary reaction harmonic bond',
    'Boundary reaction angle', 'Boundary reaction proper fourier torsion',
    'Boundary reaction improper fourier torsion', 'Bond constraint',
    'Constant bond force constraint', 'Angle constraint',
    'Harmonic torsion constraint'
}

N_PERTURBED = {'ethanol': 20, 'clch3cl': 20, 'petase_ms': 3}
PERTURB_SIGMA = 0.01  # nm


def _tol(reference):
    return 1e-6 * max(abs(reference), 1.0)


def _load_state(json_name, xyz_name):
    state = MMForceFieldGenerator.load_forcefield_from_json_file(
        str(DATA / json_name), ostream=OutputStream(None))
    state.molecule = Molecule.read_xyz_file(str(DATA / xyz_name))
    return state


def _build_pes_systems(name, conf):
    """The PES system of every state of a fixture, built by
    EvbDriver.build_systems.

    The driver is given the states and the bond sets that
    build_ffs_from_molecules would produce (EvbTestHelper.evb_seed_states).
    """

    states = [_load_state(*files) for files in FIXTURES[name]]
    evb = EvbDriver(ostream=OutputStream(None))
    EvbTestHelper.evb_seed_states(evb, states)
    if len(states) == 2:
        files = ['reactant_sys.xml', 'product_sys.xml']
    else:
        files = [f'state_{k}_sys.xml' for k in range(len(states))]
    with EvbTestHelper.evb_chdir_tmp():
        evb.build_systems(configurations=[conf], Lambda=[0.0, 1.0])
        run = Path(evb.system_confs[0]['data_folder']) / 'run'
        return [
            mm.XmlSerializer.deserialize((run / f).read_text())
            for f in files
        ]


def _geometries(name):
    """Every stored structure plus N_PERTURBED Gaussian copies of each (nm)."""

    rng = np.random.default_rng(20261022)
    geoms = []
    for _, xyz in FIXTURES[name]:
        base = Molecule.read_xyz_file(
            str(DATA / xyz)).get_coordinates_in_angstrom() * 0.1
        geoms.append(base)
        geoms += [
            base + rng.normal(0, PERTURB_SIGMA, base.shape)
            for _ in range(N_PERTURBED[name])
        ]
    return geoms


class ForceEvaluator:
    """Per-force energies of a system: every force in its own group."""

    def __init__(self, system):
        system = copy.deepcopy(system)
        assert system.getNumForces() <= 32
        self.names = []
        for index, force in enumerate(system.getForces()):
            force.setForceGroup(index)
            self.names.append(force.getName())
        self.context = mm.Context(system, mm.VerletIntegrator(1e-3),
                                  mm.Platform.getPlatformByName('Reference'))

    def terms(self, xyz):
        """{term: energy}, summed over forces of the same term, plus 'zero'
        (forces that must not contribute) and 'total'."""

        self.context.setPositions(np.asarray(xyz) * mm.unit.nanometer)
        out = dict.fromkeys(FF_TERMS, 0.0)
        out['zero'] = 0.0
        for index, name in enumerate(self.names):
            energy = self.context.getState(
                getEnergy=True, groups={index}).getPotentialEnergy(
                ).value_in_unit(mm.unit.kilojoule_per_mole)
            if name in NAME_TO_TERM:
                out[NAME_TO_TERM[name]] += energy
            else:
                assert name in ZERO_FORCES, f'unexpected force {name}'
                out['zero'] += energy
        out['total'] = self.context.getState(
            getEnergy=True).getPotentialEnergy().value_in_unit(
                mm.unit.kilojoule_per_mole)
        return out


@pytest.fixture(scope='module', params=sorted(FIXTURES))
def built(request):
    name = request.param
    states = [DATA / json_name for json_name, _ in FIXTURES[name]]
    elements = Molecule.read_xyz_file(str(DATA /
                                          FIXTURES[name][0][1])).get_labels()
    vacuum = {'name': 'vacuum', 'temperature': 300.0}
    field = {'name': 'vacuum_E', 'temperature': 300.0, 'E_field': E_FIELD}
    return SimpleNamespace(
        name=name,
        ref=FFEnergyReference(states, elements=elements),
        ref_E=FFEnergyReference(states, E_field=E_FIELD, elements=elements),
        systems=_build_pes_systems(name, vacuum),
        systems_E=_build_pes_systems(name, field),
        geoms=_geometries(name),
    )


def _mismatches(ref_terms, mm_terms, terms):
    bad = []
    for term in terms:
        diff = mm_terms[term] - ref_terms[term]
        if abs(diff) > _tol(ref_terms[term]):
            bad.append((term, ref_terms[term], mm_terms[term], diff))
    return bad


def test_q1_force_groups(built):
    """Q1: every force of each PES system equals its plan section 4 term
    (harmonic bond, Morse, angle, proper Fourier, improper Fourier in
    OpenMM's default ordering (WP-C3 D1), static and dynamic Coulomb/LJ with
    soft-core on the reactive pairs), and the total equals their sum, at the
    stored and perturbed geometries. 1e-6 rel."""

    failures = {}
    for k, system in enumerate(built.systems):
        evaluator = ForceEvaluator(system)
        for g, xyz in enumerate(built.geoms):
            mm_terms = evaluator.terms(xyz)
            ref_terms = built.ref.energies(k, xyz)
            assert abs(mm_terms['zero']) <= 1e-12
            bad = _mismatches(ref_terms, mm_terms, FF_TERMS + ('total', ))
            for term, *values in bad:
                failures.setdefault((k, term), []).append((g, *values))
    summary = {
        key: (len(rows), rows[0]) for key, rows in failures.items()
    }
    assert not failures, f'{built.name}: (state, term): (n, first) {summary}'


def test_q2_nonbonded_force_zero(built):
    """Q2: in vacuum the NonbondedForce holds only the solute, with every
    solute pair a zeroed exception, so its energy is 0 (<= 1e-9 kJ/mol)."""

    n = built.ref.n_atoms()
    for k, system in enumerate(built.systems):
        nb = [f for f in system.getForces() if f.getName() == 'NonbondedForce']
        assert len(nb) == 1
        assert nb[0].getNumExceptions() == n * (n - 1) // 2
        evaluator = ForceEvaluator(system)
        for xyz in built.geoms:
            assert abs(evaluator.terms(xyz)['nonbonded_force']) <= 1e-9


def _custom_bond_params(system, name):
    forces = [f for f in system.getForces() if f.getName() == name]
    assert len(forces) == 1, name
    force = forces[0]
    out = {}
    for b in range(force.getNumBonds()):
        i, j, params = force.getBondParameters(b)
        key = tuple(sorted((i, j)))
        assert key not in out, f'{name}: pair {key} twice'
        out[key] = tuple(params)
    return out


def _close(a, b, rtol=1e-12):
    return np.allclose(a, b, rtol=rtol, atol=0.0)


def test_q3_pair_sets(built):
    """Q3: each solute pair is excluded (1-2, 1-3 in this state's bond graph),
    1-4 scaled (coulomb 0.833, LJ 0.5) or full, and sits in exactly the
    expected CustomBondForce: dynamic for reactive pairs, static otherwise,
    with qq = scaled q_i q_j, sigma = (s_i + s_j)/2, eps = scaled
    sqrt(e_i e_j). Pairs with qq = 0 or eps = 0 are absent from that force.
    Set equality; parameters 1e-12 relative. Also the Morse and harmonic
    bond lists (D, a = sqrt(k/2D), r0; r0, k)."""

    ref = built.ref
    for k, system in enumerate(built.systems):
        pairs = ref.pairs(k)
        for kind in ('static', 'dynamic'):
            coul = _custom_bond_params(system,
                                       f'Reaction internal Coul {kind}')
            lj = _custom_bond_params(system, f'Reaction internal LJ {kind}')
            mine = {
                key: p for key, p in pairs.items()
                if p['reactive'] == (kind == 'dynamic')
                and p['cls'] != 'excluded'
            }
            exp_coul = {key for key, p in mine.items() if p['qq'] != 0.0}
            exp_lj = {key for key, p in mine.items() if p['epsilon'] != 0.0}
            assert set(coul) == exp_coul, (k, kind, set(coul) ^ exp_coul)
            assert set(lj) == exp_lj, (k, kind, set(lj) ^ exp_lj)
            for key in exp_coul:
                qq = mine[key]['qq']
                assert _close(coul[key], [qq, qq]), (k, kind, key, coul[key])
            for key in exp_lj:
                p = mine[key]
                expected = [p['sigma'], p['sigma'], p['epsilon'], p['epsilon']]
                assert _close(lj[key], expected), (k, kind, key, lj[key])

        state = ref.states[k]
        morse = _custom_bond_params(system, 'Reaction morse bond')
        harmonic_force = [
            f for f in system.getForces()
            if f.getName() == 'Static reaction harmonic bond'
        ][0]
        harmonic = {}
        for b in range(harmonic_force.getNumBonds()):
            i, j, r0, fc = harmonic_force.getBondParameters(b)
            harmonic[tuple(sorted((i, j)))] = (r0.value_in_unit(
                mm.unit.nanometer), fc.value_in_unit(
                    mm.unit.kilojoule_per_mole / mm.unit.nanometer**2))
        exp_morse, exp_harm = {}, {}
        for key, bond in state['bonds'].items():
            pair = tuple(sorted(int(x) for x in key.strip('()').split(',')))
            r0, fc = bond['equilibrium'], bond['force_constant']
            if pair in ref.reactive_pairs:
                D = bond.get('D', ref.morse_D)
                exp_morse[pair] = (D, np.sqrt(fc / (2 * D)), r0)
            else:
                exp_harm[pair] = (r0, fc)
        assert set(morse) == set(exp_morse), (k, set(morse) ^ set(exp_morse))
        assert set(harmonic) == set(exp_harm), (k,
                                                set(harmonic) ^ set(exp_harm))
        for key in exp_morse:
            assert _close(morse[key], exp_morse[key]), (k, key, morse[key])
        for key in exp_harm:
            assert _close(harmonic[key], exp_harm[key]), (k, key)


def _scan_points(rc):
    near = rc + np.array([-1e-4, -1e-5, -1e-6, 0.0, 1e-6, 1e-5, 1e-4])
    return np.concatenate([np.linspace(0.5 * rc, 2.0 * rc, 40), near])


def test_q4_soft_core_scan(built):
    """Q4: along each non-excluded reactive pair, scanned from 0.5 r_c to
    2 r_c and at r_c +- 1e-4..1e-6 nm, the dynamic Coulomb and LJ forces
    equal the plain pair potential above r_c and its second-order Taylor
    polynomial about r_c below. Coulomb r_c = 0.3 (1 + qq); LJ
    r_c = 0.7 ((26/7) 4 eps sigma^6)^(1/6), from the pair's effective
    (1-4 scaled) parameters. 1e-6 rel."""

    ref = built.ref
    n_scanned = 0
    failures = []
    for k, system in enumerate(built.systems):
        evaluator = ForceEvaluator(system)
        base = built.geoms[0]
        for (i, j), p in ref.pairs(k).items():
            if not p['reactive'] or p['cls'] == 'excluded':
                continue
            axis = base[j] - base[i]
            axis /= np.linalg.norm(axis)
            rcs = {'coul_dynamic': ref.coulomb_rc(p['qq'])}
            if p['epsilon'] != 0.0:
                rcs['lj_dynamic'] = ref.lj_rc(p['sigma'], p['epsilon'])
            for term, rc in rcs.items():
                if rc <= 0:
                    continue
                for r in _scan_points(rc):
                    xyz = base.copy()
                    xyz[j] = xyz[i] + r * axis
                    coul, lj = ref.nonbonded(k, xyz, 'dynamic')
                    want = coul if term == 'coul_dynamic' else lj
                    got = evaluator.terms(xyz)[term]
                    if abs(got - want) > _tol(want):
                        failures.append((k, (i, j), term, r - rc, got, want))
                n_scanned += 1
    assert n_scanned > 0
    assert not failures, (len(failures), failures)


def test_q5_electric_field(built):
    """Q5: the E-field force equals -F 1e-12 sum_i q_i E.x_i (q in e, E in
    V/m, x in nm, F = e N_A), and the other terms are unchanged. 1e-6 rel."""

    for k, system in enumerate(built.systems_E):
        evaluator = ForceEvaluator(system)
        for xyz in built.geoms:
            mm_terms = evaluator.terms(xyz)
            ref_terms = built.ref_E.energies(k, xyz)
            assert abs(ref_terms['efield']) > 0.01
            bad = _mismatches(ref_terms, mm_terms, ('efield', 'total'))
            assert not bad, (built.name, k, bad)


@pytest.mark.parametrize('term', ['coulomb', 'lj'])
def test_soft_core_switch_point_counted_once(term):
    """Defect C3-1 reproducer: at r == r_c exactly the soft-core pair energy
    is f(r_c), once. OpenMM's step(0) = 1, so step(r - rc) * plain +
    step(rc - r) * polynomial counted both branches there (2 f(r_c)).

    Two particles, rc chosen so that it and r = rc are exact in double."""

    from veloxchem.reactionsystembuilder import ReactionSystemBuilder

    builder = ReactionSystemBuilder(ostream=OutputStream(None))
    ref = FFEnergyReference([{'atoms': {}, 'bonds': {}}])
    if term == 'coulomb':
        qq = 0.5
        force = mm.CustomBondForce(builder._get_couloumb_expression(True))
        for p in ('qqA', 'qqB'):
            force.addPerBondParameter(p)
        params = [qq, qq]
        rc = ref.coulomb_rc(qq)
        plain = coulomb_pair(rc, qq)[0]
    else:
        sigma, epsilon = 0.34, 0.45
        force = mm.CustomBondForce(builder._get_lj_expression(True))
        for p in ('sigmaA', 'sigmaB', 'epsilonA', 'epsilonB'):
            force.addPerBondParameter(p)
        params = [sigma, sigma, epsilon, epsilon]
        rc = ref.lj_rc(sigma, epsilon)
        plain = lj_pair(rc, sigma, epsilon)[0]
    force.addGlobalParameter('l', 0.0)
    force.addBond(0, 1, params)
    system = mm.System()
    system.addParticle(1.0)
    system.addParticle(1.0)
    system.addForce(force)
    context = mm.Context(system, mm.VerletIntegrator(1e-3),
                         mm.Platform.getPlatformByName('Reference'))
    energies = []
    for r in (rc * (1 - 1e-12), rc, rc * (1 + 1e-12)):
        context.setPositions([[0.0, 0.0, 0.0], [r, 0.0, 0.0]] *
                             mm.unit.nanometer)
        energies.append(context.getState(getEnergy=True).getPotentialEnergy(
        ).value_in_unit(mm.unit.kilojoule_per_mole))
    assert np.allclose(energies, plain, rtol=1e-9, atol=0.0), (energies, plain)
