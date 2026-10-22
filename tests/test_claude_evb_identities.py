"""WP-C2: identities and switch points of the EVB systems.

The PES and integration systems that EvbDriver.build_systems writes are
checked against identities that hold for any correct build (WP-C2 plan,
section 4): rigid-motion invariance, zero net force and torque, the
reactant/product mirror, identical states, and the lambda endpoints against
the PES systems. C1 is folded in: the lambda switch of the torsions and
impropers, and the integration soft-core (non-default).

Every force of a system is evaluated in its own group on the Reference
platform (double) and mapped by name. Reference values come from
FFEnergyReference (verified in WP-C3) and from the formulas in the plan.
"""

import copy
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from veloxchem.evbdriver import EvbDriver
from veloxchem.molecule import Molecule
from veloxchem.outputstream import OutputStream
from veloxchem.reaffbuilder import ReactionForceFieldBuilder
from veloxchem.reactionsystembuilder import ReactionSystemBuilder

from test_evb_helper import EvbTestHelper
from test_claude_evb_reference import FARADAY, FFEnergyReference, _ff_key
from test_claude_evb_energy_terms import (DATA, E_FIELD, FIXTURES,
                                          _geometries, _load_state,
                                          _scan_points)

mm = pytest.importorskip('openmm')
from scipy.spatial.transform import Rotation  # noqa: E402

# Symmetric about 0.5, so that the mirror (Q3) maps the set onto itself
LAMBDAS = [0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0]
DEFAULT_SWITCH = 0.25
SWITCH_VALUES = [0.0, 0.1, 0.25, 0.5, 0.9]
VACUUM = {'name': 'vacuum', 'temperature': 300.0}
SURROGATE_FC_FACTOR = 0.2  # bonded_integration_bond_fac

# The integration systems name their static nonbonded forces differently
# from the PES systems; the terms are the same at the endpoints.
INT_TO_PES = {
    'Reaction internal coulomb': 'Reaction internal Coul static',
    'Reaction internal lennard-jones': 'Reaction internal LJ static',
}
# By design (WP-C3 table rows 1-3): the integration system holds every
# reactive pair as a dynamic harmonic bond and excludes its nonbonded terms;
# the PES holds Morse bonds and soft-core nonbonded terms for them.
INT_ONLY = {'Dynamic reaction harmonic bond'}
PES_ONLY = {
    'Reaction morse bond', 'Reaction internal Coul dynamic',
    'Reaction internal LJ dynamic'
}
TORSION_GROUPS = {
    'dihedrals': 'Reaction proper fourier torsion',
    'impropers': 'Reaction improper fourier torsion',
}

N_ROTATIONS = 10
N_PERTURBED_USED = 2  # perturbed copies per stored structure in Q1-Q3, Q5


def _tol(reference, rel=1e-9):
    return rel * max(abs(reference), 1.0)


def _elements(name):
    return Molecule.read_xyz_file(str(DATA /
                                      FIXTURES[name][0][1])).get_labels()


def _state_jsons(files):
    return [json.loads((DATA / json_name).read_text()) for json_name, _ in files]


def _seeded_driver(states):
    """A driver seeded as build_ffs_from_molecules would leave it, including
    the bond registry (evb_seed_states leaves it {}, which gives a pair
    bonded in neither state of a step a zero surrogate force constant)."""

    evb = EvbDriver(ostream=OutputStream(None))
    EvbTestHelper.evb_seed_states(evb, states)
    builder = ReactionForceFieldBuilder(ostream=OutputStream(None))
    evb.bond_registry = builder._build_reactive_bond_registry(
        states, evb.reactive_bonds)
    return evb


_BUILD_CACHE = {}


def _build(files, conf, lambdas=LAMBDAS):
    """PES systems [state k] and integration systems [step k][lambda] read
    back from the run folder that EvbDriver.build_systems writes. Cached per
    (files, configuration, lambdas); callers must not modify the systems."""

    key = (tuple(files), json.dumps(conf, sort_keys=True), tuple(lambdas))
    if key not in _BUILD_CACHE:
        _BUILD_CACHE[key] = _build_uncached(files, conf, lambdas)
    return _BUILD_CACHE[key]


def _build_uncached(files, conf, lambdas):
    states = [_load_state(*f) for f in files]
    evb = _seeded_driver(states)
    n = len(states)
    with EvbTestHelper.evb_chdir_tmp():
        evb.build_systems(configurations=[dict(conf)], Lambda=list(lambdas))
        run = Path(evb.system_confs[0]['data_folder']) / 'run'

        def read(filename):
            return mm.XmlSerializer.deserialize((run / filename).read_text())

        if n == 2:
            pes = [read('reactant_sys.xml'), read('product_sys.xml')]
            steps = [{lam: read(f'{lam:.3f}_sys.xml') for lam in lambdas}]
        else:
            pes = [read(f'state_{k}_sys.xml') for k in range(n)]
            steps = [{
                lam: read(f'step{k}_{lam:.3f}_sys.xml') for lam in lambdas
            } for k in range(n - 1)]
    return SimpleNamespace(pes=pes,
                           steps=steps,
                           registry=evb.bond_registry,
                           reactive=set(evb.reactive_bonds))


def _build_direct(files, conf, lambdas):
    """As _build, but through ReactionSystemBuilder.build_multi_step_systems
    with the arguments the driver passes, for lambda values that
    build_systems would round to 3 decimals (the switch-point scans)."""

    states = [_load_state(*f) for f in files]
    evb = _seeded_driver(states)
    builder = ReactionSystemBuilder(ostream=OutputStream(None))
    builder.water_model = evb.water_model
    steps, pes, _, _ = builder.build_multi_step_systems(
        states=states,
        lambda_vec=list(lambdas),
        configuration=dict(conf),
        reactive_bonds=evb.reactive_bonds,
        bond_registry=evb.bond_registry)
    return SimpleNamespace(pes=[pes[k] for k in range(len(states))],
                           steps=steps)


def _geoms(name, n_perturbed=N_PERTURBED_USED):
    """The stored structures plus their first n_perturbed Gaussian copies
    (from the WP-C3 geometry set, nm)."""

    from test_claude_evb_energy_terms import N_PERTURBED
    all_geoms = _geometries(name)
    block = N_PERTURBED[name] + 1
    out = []
    for start in range(0, len(all_geoms), block):
        out += all_geoms[start:start + 1 + n_perturbed]
    return out


class GroupEvaluator:
    """Per-force energies of a system by force name, each force in its own
    group, and the total forces."""

    def __init__(self, system):
        system = copy.deepcopy(system)
        assert system.getNumForces() <= 32
        self.names = []
        for index, force in enumerate(system.getForces()):
            force.setForceGroup(index)
            self.names.append(force.getName())
        self.context = mm.Context(system, mm.VerletIntegrator(1e-3),
                                  mm.Platform.getPlatformByName('Reference'))

    def _set(self, xyz):
        self.context.setPositions(np.asarray(xyz) * mm.unit.nanometer)

    def energies(self, xyz):
        """{force name: energy} (kJ/mol), plus 'total'."""

        self._set(xyz)
        out = {}
        for index, name in enumerate(self.names):
            energy = self.context.getState(
                getEnergy=True, groups={index}).getPotentialEnergy(
                ).value_in_unit(mm.unit.kilojoule_per_mole)
            out[name] = out.get(name, 0.0) + energy
        out['total'] = self.context.getState(
            getEnergy=True).getPotentialEnergy().value_in_unit(
                mm.unit.kilojoule_per_mole)
        return out

    def total(self, xyz):
        self._set(xyz)
        return self.context.getState(getEnergy=True).getPotentialEnergy(
        ).value_in_unit(mm.unit.kilojoule_per_mole)

    def forces(self, xyz):
        self._set(xyz)
        return np.asarray(
            self.context.getState(getForces=True).getForces(
                asNumpy=True).value_in_unit(mm.unit.kilojoule_per_mole /
                                            mm.unit.nanometer))


def _group_mismatches(got, want, names=None, rel=1e-9):
    names = set(got) | set(want) if names is None else names
    bad = []
    for name in sorted(names):
        a, b = got.get(name, 0.0), want.get(name, 0.0)
        if abs(a - b) > _tol(b, rel):
            bad.append((name, a, b, a - b))
    return bad


def _all_systems(built):
    out = [(f'pes{k}', s) for k, s in enumerate(built.pes)]
    for k, step in enumerate(built.steps):
        out += [(f'step{k} l={lam}', s) for lam, s in step.items()]
    return out


# ---------------------------------------------------------------------------
# Reference energies for the lambda switch (plan section 4, C1)
# ---------------------------------------------------------------------------


def switch_scaling(lam, s):
    """(reactant, product) scale of a switched torsion: a key present in one
    state only, or with a zero barrier on one side. The reactant term is off
    by lambda = 1 - s, the product term on from lambda = s (the comment at
    torsion_lambda_switch; WP-C2 decision (b)). Defined for 0 <= s < 1."""

    return max(0.0, 1 - lam / (1 - s)), max(0.0, (lam - s) / (1 - s))


def _is_switched(left, right, section, key):
    if key not in left[section] or key not in right[section]:
        return True
    return (left[section][key]['barrier'] == 0
            or right[section][key]['barrier'] == 0)


def _restricted(state, section, keys):
    """A copy of a state JSON holding only the given torsions of one
    section (no bonds, so FFEnergyReference sees no reactive pairs)."""

    out = {'atoms': state['atoms'], 'bonds': {}, 'angles': {}}
    for sec in TORSION_GROUPS:
        out[sec] = {}
    out[section] = {k: state[section][k] for k in keys}
    return out


def torsion_classes(left, right, xyz, elements):
    """Per section, the energies of the shared (linear) and switched terms
    of each side of a step, at xyz. {section: {cls: E}} with cls in
    shared_A, shared_B, switched_A, switched_B."""

    term = {'dihedrals': 'torsion', 'impropers': 'improper'}
    out = {}
    for section in TORSION_GROUPS:
        keys = set(left[section]) | set(right[section])
        switched = {k for k in keys if _is_switched(left, right, section, k)}
        energies = {}
        for side, state in (('A', left), ('B', right)):
            mine = set(state[section])
            for cls, subset in (('shared', mine - switched),
                                ('switched', mine & switched)):
                ref = FFEnergyReference(
                    [_restricted(state, section, subset)], elements=elements)
                energies[f'{cls}_{side}'] = ref.energies(0, xyz)[term[section]]
        out[section] = energies
    return out


def torsion_reference(classes, lam, s):
    r, p = switch_scaling(lam, s)
    return ((1 - lam) * classes['shared_A'] + lam * classes['shared_B'] +
            r * classes['switched_A'] + p * classes['switched_B'])


# ---------------------------------------------------------------------------
# Shared builds
# ---------------------------------------------------------------------------


@pytest.fixture(scope='module', params=sorted(FIXTURES))
def path(request):
    name = request.param
    return SimpleNamespace(name=name,
                           built=_build(FIXTURES[name], VACUUM),
                           geoms=_geoms(name))


# ---------------------------------------------------------------------------
# Q1, Q2: rigid motions, net force and torque
# ---------------------------------------------------------------------------


def _rigid_motions(seed=20261022):
    rng = np.random.default_rng(seed)
    rotations = Rotation.random(N_ROTATIONS, random_state=seed)
    shifts = rng.uniform(-2.0, 2.0, (N_ROTATIONS, 3))
    return list(zip(rotations, shifts))


def test_q1_rigid_motion_invariance(path):
    """Q1: E(Rx + t) = E(x) for every force group of every PES and
    integration system, 10 random rotations and translations (+-2 nm) of the
    stored and perturbed geometries. Plan section 4; 1e-9 relative."""

    motions = _rigid_motions()
    failures = []
    for label, system in _all_systems(path.built):
        evaluator = GroupEvaluator(system)
        for g, xyz in enumerate(path.geoms):
            centre = xyz.mean(axis=0)
            base = evaluator.energies(xyz)
            for rotation, shift in motions:
                moved = rotation.apply(xyz - centre) + centre + shift
                bad = _group_mismatches(evaluator.energies(moved), base)
                failures += [(label, g, *b) for b in bad]
    assert not failures, (len(failures), failures[:5])


def _net_force_and_torque(forces, xyz):
    centre = xyz.mean(axis=0)
    net = forces.sum(axis=0)
    torque = np.cross(xyz - centre, forces).sum(axis=0)
    scale = np.abs(forces).sum()
    return net, torque, scale


def test_q2_net_force_and_torque(path):
    """Q2: sum_i F_i = 0 and sum_i (x_i - x_c) x F_i = 0 for every system
    and geometry, to 1e-9 x sum |F_i| (x 1 nm for the torque). Plan
    section 4."""

    failures = []
    for label, system in _all_systems(path.built):
        evaluator = GroupEvaluator(system)
        for g, xyz in enumerate(path.geoms):
            net, torque, scale = _net_force_and_torque(
                evaluator.forces(xyz), xyz)
            for what, vec in (('net', net), ('torque', torque)):
                if np.linalg.norm(vec) > 1e-9 * scale:
                    failures.append((label, g, what, np.linalg.norm(vec),
                                     scale))
    assert not failures, (len(failures), failures[:5])


def test_q2_electric_field_net_force():
    """Q2 positive control: with an external field the net force is
    F 1e-12 (sum_i q_i) E (kJ/mol/nm), non-zero for Cl- + CH3Cl (net charge
    -1); the torque about the centroid is no longer zero either. The PES and
    every integration system, 1e-6 relative."""

    name = 'clch3cl'
    field = {'name': 'vacuum_E', 'temperature': 300.0, 'E_field': E_FIELD}
    built = _build(FIXTURES[name], field)
    expected = FARADAY * 1e-12 * -1.0 * np.asarray(E_FIELD)
    assert np.linalg.norm(expected) > 1.0
    for label, system in _all_systems(built):
        evaluator = GroupEvaluator(system)
        for xyz in _geoms(name):
            net = evaluator.forces(xyz).sum(axis=0)
            assert np.allclose(net, expected, rtol=1e-6, atol=0.0), (label,
                                                                     net,
                                                                     expected)


# ---------------------------------------------------------------------------
# Q3: reactant/product mirror
# ---------------------------------------------------------------------------


@pytest.fixture(scope='module', params=sorted(FIXTURES))
def mirror(request):
    """Forward and reversed path, at the default switch and at s = 0.5."""

    name = request.param
    files = FIXTURES[name]
    out = SimpleNamespace(name=name, geoms=_geoms(name),
                          jsons=_state_jsons(files),
                          elements=_elements(name))
    for s in (DEFAULT_SWITCH, 0.5):
        conf = VACUUM if s == DEFAULT_SWITCH else dict(
            VACUUM, torsion_lambda_switch=s)
        out.__dict__[s] = (_build(files, conf), _build(files[::-1], conf))
    return out


def _mirror_differences(mirror, s):
    """Per (step, lambda, geometry): forward minus mirrored group energies,
    and the PES differences."""

    forward, reverse = mirror.__dict__[s]
    n = len(forward.pes)
    pes_bad = []
    for k in range(n):
        a = GroupEvaluator(forward.pes[k])
        b = GroupEvaluator(reverse.pes[n - 1 - k])
        for g, xyz in enumerate(mirror.geoms):
            bad = _group_mismatches(b.energies(xyz), a.energies(xyz))
            pes_bad += [(k, g, *x) for x in bad]
    rows = []
    for k in range(n - 1):
        for lam in LAMBDAS:
            a = GroupEvaluator(forward.steps[k][lam])
            b = GroupEvaluator(reverse.steps[n - 2 - k][round(1 - lam, 3)])
            for g, xyz in enumerate(mirror.geoms):
                ea, eb = a.energies(xyz), b.energies(xyz)
                rows.append((k, lam, g, xyz, ea, eb))
    return pes_bad, rows


@pytest.mark.parametrize('s', [DEFAULT_SWITCH, 0.5])
def test_q3_mirror(mirror, s):
    """Q3: the reversed path reproduces every force group, the switched
    torsions included: integration step k at lambda = reversed step N-2-k
    at 1 - lambda, and PES_k = PES'_{N-1-k}, at the default
    torsion_lambda_switch and at 0.5. Plan section 4; 1e-9 relative."""

    pes_bad, rows = _mirror_differences(mirror, s)
    assert not pes_bad, pes_bad[:5]
    failures = []
    for k, lam, g, _, ea, eb in rows:
        failures += [(k, lam, g, *b) for b in _group_mismatches(eb, ea)]
    assert not failures, (len(failures), failures[:5])


# ---------------------------------------------------------------------------
# Q4: identical states
# ---------------------------------------------------------------------------


@pytest.mark.parametrize('name', ['ethanol', 'clch3cl'])
def test_q4_identical_states(name):
    """Q4: a path S_0 = S_1 (two loads of the reactant) has no reactive
    pairs, PES_0 = PES_1 per group, a gap E1 - E2 of 0, and an integration
    system equal in total energy to PES_0 at every lambda. Plan section 4;
    1e-9 relative, gap <= 1e-9 kJ/mol."""

    files = [FIXTURES[name][0], FIXTURES[name][0]]
    built = _build(files, VACUUM)
    assert built.reactive == set()
    for system in [built.pes[0]] + list(built.steps[0].values()):
        for force in system.getForces():
            if force.getName() in PES_ONLY | INT_ONLY:
                assert force.getNumBonds() == 0, force.getName()
    pes0, pes1 = (GroupEvaluator(s) for s in built.pes)
    ints = {lam: GroupEvaluator(s) for lam, s in built.steps[0].items()}
    geoms = _geoms(name)[:1 + N_PERTURBED_USED]  # the reactant's structures
    for xyz in geoms:
        e0, e1 = pes0.energies(xyz), pes1.energies(xyz)
        assert not _group_mismatches(e1, e0)
        assert abs(e0['total'] - e1['total']) <= 1e-9
        for lam, evaluator in ints.items():
            total = evaluator.total(xyz)
            assert abs(total - e0['total']) <= _tol(e0['total']), (lam, total,
                                                                  e0['total'])


# ---------------------------------------------------------------------------
# Q5: endpoints against the PES systems
# ---------------------------------------------------------------------------


def surrogate_bond_energy(ref, registry, k, xyz):
    """E of the dynamic harmonic bonds at an end of a step, where the end is
    state k (plan section 4): a reactive pair bonded in state k keeps that
    bond; otherwise r0 = (sigma_i + sigma_j)/2 of state k and
    k = 0.2 fc_ref (registry)."""

    state = ref.states[k]
    bonds = {tuple(sorted(_ff_key(key))): b for key, b in state['bonds'].items()}
    energy = 0.0
    for pair in ref.reactive_pairs:
        i, j = pair
        r = np.linalg.norm(xyz[i] - xyz[j])
        if pair in bonds:
            r0, fc = bonds[pair]['equilibrium'], bonds[pair]['force_constant']
        else:
            r0 = 0.5 * (ref.atom(k, i)['sigma'] + ref.atom(k, j)['sigma'])
            fc = SURROGATE_FC_FACTOR * registry[pair]['fc_ref']
        if fc > 0:
            energy += 0.5 * fc * (r - r0)**2
    return energy


def _endpoint_check(built, ref, geoms, int_names=None):
    """Rows (step, end, geometry, problem) for the endpoint identity."""

    failures = []
    for k, step in enumerate(built.steps):
        for lam, state in ((0.0, k), (1.0, k + 1)):
            e_int_sys = GroupEvaluator(step[lam])
            e_pes_sys = GroupEvaluator(built.pes[state])
            for g, xyz in enumerate(geoms):
                e_int = e_int_sys.energies(xyz)
                e_pes = e_pes_sys.energies(xyz)
                mapped = {INT_TO_PES.get(n, n): v for n, v in e_int.items()}
                shared = (set(mapped) | set(e_pes)) - INT_ONLY - PES_ONLY - {
                    'total'
                }
                if int_names is not None:
                    shared &= int_names
                for b in _group_mismatches(mapped, e_pes, names=shared):
                    failures.append((k, lam, g, 'shared', *b))
                if int_names is not None:
                    continue
                terms = ref.energies(state, xyz)
                dynharm = surrogate_bond_energy(ref, built.registry, state,
                                                xyz)
                expected = dynharm - (terms['morse'] + terms['coul_dynamic'] +
                                      terms['lj_dynamic'])
                measured = e_int['total'] - e_pes['total']
                if abs(measured - expected) > _tol(expected, 1e-6):
                    failures.append(
                        (k, lam, g, 'remainder', measured, expected))
                got = e_int.get('Dynamic reaction harmonic bond', 0.0)
                if abs(got - dynharm) > _tol(dynharm, 1e-6):
                    failures.append((k, lam, g, 'dynamic bond', got, dynharm))
    return failures


def test_q5_endpoints(path):
    """Q5: integration step k at lambda 0 (1) equals PES_k (PES_k+1) in
    every shared force group (1e-9 relative), and the rest closes
    E_int - E_PES = E_dynharm - E_morse - E_dynCoul - E_dynLJ with the
    FFEnergyReference terms and the surrogate formula of plan section 4
    (1e-6 relative). Stored and perturbed geometries."""

    name = path.name
    ref = FFEnergyReference([DATA / j for j, _ in FIXTURES[name]],
                            elements=_elements(name))
    failures = _endpoint_check(path.built, ref, path.geoms)
    assert not failures, (len(failures), failures[:5])


# ---------------------------------------------------------------------------
# Q7: lambda switch of torsions and impropers (C1)
# ---------------------------------------------------------------------------

SWITCHED_FIXTURES = ['ethanol', 'petase_ms']  # Cl- + CH3Cl has no torsions


def _scan_lambdas(s):
    grid = set(np.round(np.linspace(0.0, 1.0, 51), 12))
    for point in (s, 1 - s):
        grid |= {point - 1e-6, point, point + 1e-6}
    return sorted(x for x in grid if 0.0 <= x <= 1.0)


@pytest.mark.parametrize('name', SWITCHED_FIXTURES)
def test_q7_torsion_switch_scan(name):
    """Q7: on a lambda scan (0.02 steps plus s, 1 - s and +-1e-6 around
    them, s = 0.25), the torsion and improper groups of every step follow
    plan section 4: shared keys (1 - l, l), switched keys (r(l), p(l)), per
    class from FFEnergyReference, 1e-9 relative; and they are continuous at
    the switch points (jump <= |dE/dl| 2e-6 + 1e-9 kJ/mol)."""

    s = DEFAULT_SWITCH
    lambdas = _scan_lambdas(s)
    files = FIXTURES[name]
    jsons, elements = _state_jsons(files), _elements(name)
    built = _build_direct(files, VACUUM, lambdas)
    failures = []
    for k, step in enumerate(built.steps):
        xyz = _geoms(name, 0)[k]
        classes = torsion_classes(jsons[k], jsons[k + 1], xyz, elements)
        measured = {group: {} for group in TORSION_GROUPS.values()}
        for lam in lambdas:
            energies = GroupEvaluator(step[lam]).energies(xyz)
            for section, group in TORSION_GROUPS.items():
                got = energies.get(group, 0.0)
                want = torsion_reference(classes[section], lam, s)
                measured[group][lam] = got
                if abs(got - want) > _tol(want):
                    failures.append((k, lam, group, got, want))
        for section, group in TORSION_GROUPS.items():
            c = classes[section]
            slope = (abs(c['shared_A']) + abs(c['shared_B']) +
                     abs(c['switched_A']) / (1 - s) +
                     abs(c['switched_B']) / (1 - s))
            for point in (s, 1 - s):
                jump = (measured[group][point + 1e-6] -
                        measured[group][point - 1e-6])
                if abs(jump) > slope * 2e-6 + 1e-9:
                    failures.append((k, 'jump', point, group, jump))
    assert not failures, (len(failures), failures[:5])


def _torsion_endpoint_mismatches(name, s):
    files = FIXTURES[name]
    conf = dict(VACUUM, torsion_lambda_switch=s)
    built = _build_direct(files, conf, [0.0, 1.0])
    built.registry = None
    return _endpoint_check(built, None, _geoms(name, 0),
                           int_names=set(TORSION_GROUPS.values()))


@pytest.mark.parametrize('name', SWITCHED_FIXTURES)
@pytest.mark.parametrize('s', [x for x in SWITCH_VALUES if x != 1.0])
def test_q7_torsion_switch_endpoints(name, s):
    """Q7: at lambda 0 and 1 the torsion and improper groups of the
    integration systems equal those of the PES systems, for every
    torsion_lambda_switch in (0, 0.1, 0.25, 0.5, 0.9). 1e-9 relative."""

    failures = _torsion_endpoint_mismatches(name, s)
    assert not failures, (len(failures), failures[:5])


@pytest.mark.parametrize('s', [1.0, -0.1])
def test_q7_torsion_switch_out_of_range(s):
    """Defect C2-1: torsion_lambda_switch outside [0, 1) used to be accepted;
    s = 1 dropped the switched torsions at both endpoints. It is refused."""

    with pytest.raises(Exception, match='torsion_lambda_switch'):
        _build_direct(FIXTURES['ethanol'],
                      dict(VACUUM, torsion_lambda_switch=s), [0.0, 1.0])


# ---------------------------------------------------------------------------
# Q8: integration soft-core (C1, non-default)
# ---------------------------------------------------------------------------


def test_q8_integration_soft_core_scan():
    """Q8: with soft_core_coulomb_int and soft_core_lj_int on, the static
    (non-reactive) Coulomb and LJ forces of the ethanol integration systems
    at lambda 0, 0.5, 1 equal (1 - l) f_A + l f_B, each side the plain pair
    potential above its own r_c and its second-order Taylor polynomial about
    r_c below (WP-C3 definition), along scans of one full and one 1-4 pair
    across each side's r_c (r_c exactly and +-1e-6..1e-4 nm). 1e-6
    relative, or 1e-6 kJ/mol below 1 kJ/mol."""

    name = 'ethanol'
    lambdas = [0.0, 0.5, 1.0]
    conf = dict(VACUUM, soft_core_coulomb_int=True, soft_core_lj_int=True)
    built = _build(FIXTURES[name], conf, lambdas)
    ref = FFEnergyReference([DATA / j for j, _ in FIXTURES[name]],
                            soft_core_flags={
                                'coulomb_static': True,
                                'lj_static': True
                            },
                            elements=_elements(name))
    base = _geoms(name, 0)[0]
    pairs = ref.pairs(0)
    chosen = [
        next(key for key, p in pairs.items()
             if not p['reactive'] and p['cls'] == cls)
        for cls in ('full', '1-4')
    ]
    evaluators = {lam: GroupEvaluator(s) for lam, s in built.steps[0].items()}
    groups = {
        'coul': 'Reaction internal coulomb soft-core',
        'lj': 'Reaction internal lennard-jones soft-core'
    }
    failures = []
    n_points = 0
    for i, j in chosen:
        axis = base[j] - base[i]
        axis /= np.linalg.norm(axis)
        rcs = []
        for k in (0, 1):
            p = ref.pairs(k)[(i, j)]
            if p['cls'] == 'excluded':
                continue
            rcs.append(ref.coulomb_rc(p['qq']))
            if p['epsilon'] != 0.0:
                rcs.append(ref.lj_rc(p['sigma'], p['epsilon']))
        for rc in rcs:
            for r in _scan_points(rc):
                xyz = base.copy()
                xyz[j] = xyz[i] + r * axis
                sides = [ref.nonbonded(k, xyz, 'static') for k in (0, 1)]
                for lam, evaluator in evaluators.items():
                    energies = evaluator.energies(xyz)
                    for t, (term, group) in enumerate(groups.items()):
                        want = (1 - lam) * sides[0][t] + lam * sides[1][t]
                        got = energies[group]
                        if abs(got - want) > _tol(want, 1e-6):
                            failures.append(
                                ((i, j), lam, term, r - rc, got, want))
                    n_points += 1
    assert n_points > 0
    assert not failures, (len(failures), failures[:5])
