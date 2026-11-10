"""WP-E1 dev tests: residue cutting in the PDB route of ReactionSystemBuilder.

Dev-only (test_claude_*): never to master. Until the DhlA fixture arrives the
cut is exercised on a fragment cut out of the PETase fixture: Asp177 at 'beta'
(formate + H cap) and His208 at ('CB', 'CG') (imidazole + H cap), with a proton
transfer HD1 (His ND1) -> the nearer Asp oxygen. Charges are set by hand (no
QM); the total is -1 in both states.
"""

import json
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np
import pytest

from veloxchem.molecule import Molecule
from veloxchem.mmforcefieldgenerator import MMForceFieldGenerator
from veloxchem.reactionsystembuilder import ReactionSystemBuilder
from veloxchem.outputstream import OutputStream

from test_evb_helper import EvbTestHelper

mm = pytest.importorskip('openmm')
import openmm.app as mmapp  # noqa: E402
import openmm.unit as mmunit  # noqa: E402

DATA = EvbTestHelper.evb_data_dir()
CIF = DATA / 'evb_petase_A60S.cif'
ASP = ['CG', 'OD1', 'OD2']
HIS = ['CG', 'ND1', 'HD1', 'CE1', 'HE1', 'NE2', 'CD2', 'HD2']
CAP_ASP, CAP_HIS = 3, 12
LAMBDAS = [0.0, 0.4, 1.0]
PETASE_RES = [{'chain': 'B', 'residue': 1}, {'chain': 'A', 'residue': 208},
              {'chain': 'A', 'residue': 177, 'reacting_atoms': ['OD2']},
              {'chain': 'A', 'residue': 131}]


# --------------------------------------------------------------------------
# fixture construction
# --------------------------------------------------------------------------
def _petase_topology():
    cif = mmapp.PDBxFile(str(CIF))
    return cif.topology, cif.positions


def _residue(top, chain, resid):
    return [r for c in top.chains() if c.id == chain
            for r in c.residues() if r.id == str(resid)][0]


def _cut_fragments(swap=False):
    """Reactant/product force fields of the cut fragment; swap=True lists the
    Asp oxygens as OD2, OD1 (reproducer of E1-2)."""
    top, pos = _petase_topology()

    def coords(res):
        return {a.name: np.array(pos[a.index].value_in_unit(mmunit.angstrom))
                for a in res.atoms()}

    asp, his = coords(_residue(top, 'A', 177)), coords(_residue(top, 'A', 208))

    def cap(c, length):
        v = c['CB'] - c['CG']
        return c['CG'] + length * v / np.linalg.norm(v)

    asp_names = ['CG', 'OD2', 'OD1'] if swap else ASP
    acceptor = min(('OD1', 'OD2'),
                   key=lambda o: np.linalg.norm(asp[o] - his['HD1']))
    labels = ['C', 'O', 'O', 'H', 'C', 'N', 'H', 'C', 'H', 'N', 'C', 'H', 'H']
    rea = np.array([asp[n] for n in asp_names] + [cap(asp, 1.09)] +
                   [his[n] for n in HIS] + [cap(his, 1.08)])
    pro = rea.copy()
    i_h, i_o = 4 + HIS.index('HD1'), asp_names.index(acceptor)
    v = his['ND1'] - asp[acceptor]
    pro[i_h] = asp[acceptor] + 0.98 * v / np.linalg.norm(v)
    q_rea = [0.70, -0.80, -0.80, -0.10,
             0.05, -0.35, 0.35, 0.15, 0.10, -0.50, -0.05, 0.10, 0.15]
    q_acid = {acceptor: -0.55, ({'OD1', 'OD2'} - {acceptor}).pop(): -0.50}
    q_pro = ([0.55] + [q_acid[n] for n in asp_names[1:]] + [0.05] +
             [0.00, -0.55, 0.45, 0.00, 0.05, -0.55, -0.05, 0.05, 0.05])
    ffs = []
    for xyz, q in ((rea, q_rea), (pro, q_pro)):
        text = f"{len(labels)}\n\n" + "\n".join(
            f"{a} {x:.6f} {y:.6f} {z:.6f}" for a, (x, y, z) in zip(labels, xyz))
        ff = MMForceFieldGenerator(ostream=OutputStream(None))
        ff.partial_charges = np.array(q)
        ff.create_topology(Molecule.read_xyz_string(text))
        ffs.append(ff)
    info = {'acceptor': acceptor, 'acceptor_index': i_o, 'hd1': i_h,
            'asp_names': asp_names, 'q_rea': q_rea, 'q_pro': q_pro}
    return ffs[0], ffs[1], info


def _cut_config(reacting=('OD2', ), **extra):
    """Cut config; reacting is Asp177's 'reacting_atoms' (None: no key).
    OD2 is the oxygen nearest HD1 in the PETase structure."""
    asp = {'chain': 'A', 'residue': 177, 'cut': 'beta'}
    if reacting is not None:
        asp['reacting_atoms'] = (reacting if isinstance(reacting, str) else
                                 list(reacting))
    cfg = {"name": "petase_cut", "temperature": 300.0, "pdb": str(CIF),
           "pdb_active_res": [asp, {'chain': 'A', 'residue': 208,
                                    'cut': ('CB', 'CG')}]}
    cfg.update(extra)
    return cfg


def _build(reactant, product, config, out=None):
    builder = ReactionSystemBuilder(
        ostream=OutputStream(str(out) if out else None))
    builder.water_model = "cspce"
    with EvbTestHelper.evb_chdir_tmp():
        systems, topology, positions = builder.build_systems(
            reactant, product, list(LAMBDAS), config)
    return builder, systems, topology


@pytest.fixture(scope='module')
def cut_build(tmp_path_factory):
    out = tmp_path_factory.mktemp('cut') / 'build.out'
    rea, pro, info = _cut_fragments()
    builder, systems, topology = _build(rea, pro, _cut_config(), out)
    return builder, systems, topology, info, out.read_text()


@pytest.fixture(scope='module')
def petase_build():
    def load(n):
        f = MMForceFieldGenerator.load_forcefield_from_json_file(
            str(DATA / f'evb_petase_{n}_ff_data.json'),
            ostream=OutputStream(None))
        f.molecule = Molecule.read_xyz_file(str(DATA / f'evb_petase_{n}.xyz'))
        return f

    cfg = {"name": "petase", "temperature": 300.0, "pdb": str(CIF),
           "pdb_active_res": PETASE_RES}
    return _build(load('reactant'), load('product'), cfg)


# --------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------
def _nb(system):
    return next(f for f in system.getForces()
                if isinstance(f, mm.NonbondedForce))


def _charges(system):
    nb = _nb(system)
    return np.array([nb.getParticleParameters(i)[0].value_in_unit(
        mmunit.elementary_charge) for i in range(nb.getNumParticles())])


def _ff14sb_charges():
    path = Path(mmapp.__file__).parent / 'data' / 'amber14' / 'protein.ff14SB.xml'
    root = ET.parse(path).getroot()
    return {r.get('name'): {a.get('name'): float(a.get('charge'))
                            for a in r.findall('Atom')}
            for r in root.find('Residues')}


def _bond_distance_3(topology):
    adjacency = {}
    for bond in topology.bonds():
        adjacency.setdefault(bond[0].index, set()).add(bond[1].index)
        adjacency.setdefault(bond[1].index, set()).add(bond[0].index)

    def distances(start):
        dist, shell = {start: 0}, [start]
        for depth in (1, 2, 3):
            shell = [n for a in shell for n in adjacency.get(a, ())
                     if n not in dist]
            for n in shell:
                dist.setdefault(n, depth)
        return dist

    return distances


def _check_cross_14(builder, system, topology, changed, fragment):
    """Every 1-4 exception (PDB topology) with a changed-charge atom, outside
    the fragment, equals the AMBER 1-4 rule on the current charges; every
    1-2/1-3 exception of a changed atom stays excluded."""
    nb = _nb(system)
    distances = _bond_distance_3(topology)
    n14 = 0
    for k in range(nb.getNumExceptions()):
        p1, p2, qq, sig, eps = nb.getExceptionParameters(k)
        if p1 not in changed and p2 not in changed:
            continue
        if p1 in fragment and p2 in fragment:
            continue
        d = distances(p1).get(p2)
        qq = qq.value_in_unit(mmunit.elementary_charge**2)
        eps = eps.value_in_unit(mmunit.kilojoule_per_mole)
        if d in (1, 2):
            assert qq == 0.0 and eps == 0.0, (p1, p2, d)
            continue
        if d != 3:
            continue
        n14 += 1
        q1, s1, e1 = nb.getParticleParameters(p1)
        q2, s2, e2 = nb.getParticleParameters(p2)
        q1 = q1.value_in_unit(mmunit.elementary_charge)
        q2 = q2.value_in_unit(mmunit.elementary_charge)
        e1 = e1.value_in_unit(mmunit.kilojoule_per_mole)
        e2 = e2.value_in_unit(mmunit.kilojoule_per_mole)
        s = 0.5 * (s1 + s2).value_in_unit(mmunit.nanometer)
        assert qq == pytest.approx(q1 * q2 / 1.2, abs=1e-9), (p1, p2)
        assert eps == pytest.approx(0.5 * np.sqrt(e1 * e2), abs=1e-9)
        assert sig.value_in_unit(mmunit.nanometer) == pytest.approx(s, abs=1e-9)
    return n14


def _changed_atoms(builder):
    fragment = {a.index for a in builder.reaction_atoms.values()}
    changed = set(fragment)
    for cut in builder._cuts:
        changed.add(cut['host'].index)
        changed.update(a.index for a in cut['rest_atoms'])
    return changed, fragment


# --------------------------------------------------------------------------
# Q1: specification
# --------------------------------------------------------------------------
@pytest.fixture(scope='module')
def petase_residues():
    top, _ = _petase_topology()
    first = {}
    for res in top.residues():
        if res.chain.id == 'A':
            first.setdefault(res.name, res)
    return first


@pytest.mark.parametrize("name, cut, host, evb", [
    ('ASP', 'beta', 'CB', ['CG', 'OD1', 'OD2']),
    ('ASP', ('CB', 'CG'), 'CB', ['CG', 'OD1', 'OD2']),
    ('ASP', ['CG', 'CB'], 'CB', ['CG', 'OD1', 'OD2']),
    ('ASP', 'alpha', 'CA', ['CB', 'HB2', 'HB3', 'CG', 'OD1', 'OD2']),
    ('GLU', 'gamma', 'CG', ['CD', 'OE1', 'OE2']),
    ('LYS', 'gamma', 'CG',
     ['CD', 'HD2', 'HD3', 'CE', 'HE2', 'HE3', 'NZ', 'HZ1', 'HZ2', 'HZ3']),
    ('SER', 'beta', 'CB', ['OG', 'HG']),
    ('HIS', 'beta', 'CB', None),
])
def test_cut_specification(petase_residues, name, cut, host, evb):
    builder = ReactionSystemBuilder(ostream=OutputStream(None))
    res = petase_residues[name]
    result = builder._resolve_residue_cut(res, cut)
    assert result['host'].name == host
    names = [a.name for a in result['evb_atoms']]
    if evb is not None:
        assert sorted(names) == sorted(evb)
    assert 'CA' not in names
    all_names = {a.name for a in res.atoms()}
    rest = {a.name for a in result['rest_atoms']}
    assert rest | set(names) | {host} == all_names
    assert not rest & set(names) and host not in rest


@pytest.mark.parametrize("name, cut, match", [
    ('THR', 'beta', 'heavy neighbors'),
    ('VAL', 'beta', 'heavy neighbors'),
    ('ILE', 'beta', 'heavy neighbors'),
    ('PHE', 'gamma', 'heavy neighbors'),
    ('GLY', 'alpha', 'no CB'),
    ('PRO', 'alpha', 'ring'),
    ('ASP', ('CA', 'CG'), 'not bonded'),
    ('ASP', ('CB', 'XX'), 'no atom named'),
    ('ASP', 'delta', 'unknown cut'),
    ('ASP', ('CB', ), 'pair of atom names'),
])
def test_cut_refusals(petase_residues, name, cut, match):
    builder = ReactionSystemBuilder(ostream=OutputStream(None))
    with pytest.raises(Exception, match=match):
        builder._resolve_residue_cut(petase_residues[name], cut)


def test_cut_of_unmatched_residue_refused():
    rea, pro, _ = _cut_fragments()
    cfg = _cut_config()
    cfg['pdb_active_res'] = [{'chain': 'B', 'residue': 1, 'cut': 'beta'}]
    with pytest.raises(Exception, match='matches no amber14 template'):
        _build(rea, pro, cfg)


# --------------------------------------------------------------------------
# Q2: matching, caps, warnings
# --------------------------------------------------------------------------
def test_cut_mapping_and_caps(cut_build):
    builder, systems, topology, info, out = cut_build
    names = {i: (a.residue.id, a.name)
             for i, a in builder.reaction_atoms.items() if i < 13}
    for k, n in enumerate(info['asp_names']):
        assert names[k] == ('177', n)
    for k, n in enumerate(HIS):
        assert names[4 + k] == ('208', n)
    hosts = {cap: (h.residue.id, h.name)
             for cap, h in builder._cut_caps.items()}
    assert hosts == {CAP_ASP: ('177', 'CB'), CAP_HIS: ('208', 'CB')}
    assert builder.boundary_atoms == {}
    for cap in (CAP_ASP, CAP_HIS):
        assert all(s.atoms[cap]['pdb'] == 'unmapped'
                   for s in (builder.reactant, builder.product))
    # 3-bond rule: the reacting O and N are 2 bonds from their hosts
    assert out.count('lie within 2 bonds of the host') == 2
    assert 'ASP177:' + info['acceptor'] in out and 'HIS208:ND1' in out


@pytest.mark.parametrize("swap, reacting, expected", [
    (True, 'closest', 'nearest'),
    (True, ('OD2', ), 'OD2'),
    (False, ('OD1', ), 'OD1'),
])
def test_reacting_atoms(swap, reacting, expected):
    """E1-2: the two Asp oxygens are equivalent to the graph match; the proton
    acceptor goes where 'reacting_atoms' says, whatever the atom order of the
    fragment ('closest': the oxygen nearest HD1)."""
    rea, pro, info = _cut_fragments(swap=swap)
    builder, _, _ = _build(rea, pro, _cut_config(reacting))
    want = info['acceptor'] if expected == 'nearest' else expected
    assert builder.reaction_atoms[info['acceptor_index']].name == want


@pytest.mark.parametrize("reacting, match", [
    (None, "can be placed on symmetry-equivalent PDB atoms"),
    (('CG', ), "on atoms of a forming or breaking bond"),
    (('CB', ), "no atom named CB in the part of the residue"),
    ('nearest', "must be 'closest' or a list of atom names"),
])
def test_reacting_atoms_refusals(reacting, match):
    rea, pro, _ = _cut_fragments()
    with pytest.raises(Exception, match=match):
        _build(rea, pro, _cut_config(reacting))


# --------------------------------------------------------------------------
# Q3: charges
# --------------------------------------------------------------------------
def test_cut_charges(cut_build):
    builder, systems, topology, info, _ = cut_build
    q_rea, q_pro = np.array(info['q_rea']), np.array(info['q_pro'])
    amber = _ff14sb_charges()
    template = {'177': amber['ASP'], '208': amber['HID']}
    evb = [i for i in range(13) if i not in (CAP_ASP, CAP_HIS)]
    cases = [(lam, systems[lam], (1 - lam) * q_rea + lam * q_pro)
             for lam in LAMBDAS]
    cases += [('reactant', systems['reactant'], q_rea),
              ('product', systems['product'], q_pro)]
    rest_sums = []
    for label, system, q in cases:
        charges = _charges(system)
        for i in evb:
            assert charges[builder.reaction_atoms[i].index] == pytest.approx(
                q[i], abs=1e-12), (label, i)
        for cap, host in builder._cut_caps.items():
            assert charges[host.index] == pytest.approx(q[cap], abs=1e-12)
        for cut in builder._cuts:
            rest = cut['rest_atoms']
            ref = template[cut['residue'].id]
            shifts = [charges[a.index] - ref[a.name] for a in rest]
            assert np.ptp(shifts) < 1e-12
            assert shifts[0] == pytest.approx(cut['shift'], abs=1e-12)
            rest_sums.append(sum(charges[a.index] for a in rest))
        total = charges.sum()
        assert total == pytest.approx(round(total), abs=1e-9), label
    assert np.allclose(rest_sums, np.round(rest_sums), atol=1e-9)
    asp = next(c for c in builder._cuts if c['residue'].id == '177')
    assert asp['shift'] == pytest.approx(0.1663 / 8, abs=1e-4)


def test_cut_charges_equal_water(cut_build):
    """The same fragment without the protein: every fragment atom carries the
    charge it has in the cut build, and the cap H the charge of the host."""
    builder, systems, _, info, _ = cut_build
    rea, pro, _ = _cut_fragments()
    vac, vac_systems, _ = _build(rea, pro, {"name": "vacuum",
                                            "temperature": 300.0})
    for lam in LAMBDAS:
        cut_q, vac_q = _charges(systems[lam]), _charges(vac_systems[lam])
        for i in range(13):
            index = (builder._cut_caps[i].index if i in builder._cut_caps
                     else builder.reaction_atoms[i].index)
            assert cut_q[index] == pytest.approx(
                vac_q[vac.reaction_atoms[i].index], abs=1e-12), (lam, i)


# --------------------------------------------------------------------------
# Q4: bonded terms
# --------------------------------------------------------------------------
def test_cut_bonded_terms(cut_build):
    builder, systems, _, _, _ = cut_build
    evb = {a.index for a in builder.reaction_atoms.values()}
    hosts = {h.index for h in builder._cut_caps.values()}

    def terms(system, name):
        f = next(f for f in system.getForces() if f.getName() == name)
        if isinstance(f, mm.HarmonicBondForce):
            return [f.getBondParameters(i) for i in range(f.getNumBonds())]
        if isinstance(f, mm.HarmonicAngleForce):
            return [f.getAngleParameters(i) for i in range(f.getNumAngles())]
        return [f.getTorsionParameters(i) for i in range(f.getNumTorsions())]

    for name, n in (('PDB Bond Force', 2), ('PDB Angle Force', 3),
                    ('PDB Torsion Force', 4)):
        ref = terms(systems[0.0], name)
        for lam in LAMBDAS[1:] + ['reactant', 'product']:
            assert [str(t) for t in terms(systems[lam], name)] == [
                str(t) for t in ref]
        crossing = inside = 0
        for t in ref:
            atoms, k = t[:n], t[-1]
            k = k.value_in_unit(k.unit) if hasattr(k, 'unit') else k
            if all(a in evb for a in atoms):
                assert k == 0, (name, atoms)
                inside += 1
            elif any(a in evb for a in atoms) and any(a in hosts
                                                      for a in atoms):
                crossing += k != 0
        assert inside > 0 and crossing > 0, name
    # no reaction term is built on a host (the cap is not a particle)
    for f in systems[0.4].getForces():
        if f.getName().startswith('Reaction') and hasattr(f, 'getNumBonds'):
            for i in range(f.getNumBonds()):
                p = f.getBondParameters(i)
                assert p[0] not in hosts and p[1] not in hosts


# --------------------------------------------------------------------------
# Q5/Q6: 1-4 exceptions
# --------------------------------------------------------------------------
def test_cut_14_exceptions(cut_build):
    builder, systems, topology, _, _ = cut_build
    changed, fragment = _changed_atoms(builder)
    for key in LAMBDAS + ['reactant', 'product']:
        n14 = _check_cross_14(builder, systems[key], topology, changed,
                              fragment)
        assert n14 == len(builder._pdb_14_exceptions) > 0


def test_whole_residue_14_exceptions(petase_build):
    """E1-1: on the whole-residue PETase path every cross 1-4 exception
    follows the current charges (22 of 60 kept AMBER products before)."""
    builder, systems, topology = petase_build
    changed, fragment = _changed_atoms(builder)
    for key in LAMBDAS + ['reactant', 'product']:
        n14 = _check_cross_14(builder, systems[key], topology, changed,
                              fragment)
        assert n14 == 60


# --------------------------------------------------------------------------
# Q8: implicit solvent, options round trip, finite energies
# --------------------------------------------------------------------------
def test_cut_implicit_solvent_roundtrip():
    rea, pro, info = _cut_fragments()
    cfg = json.loads(json.dumps(_cut_config(implicit_solvent_model='gbn2')))
    assert cfg['pdb_active_res'][1]['cut'] == ['CB', 'CG']
    builder, systems, _ = _build(rea, pro, cfg)
    assert sorted(builder._cut_caps) == [CAP_ASP, CAP_HIS]
    platform = mm.Platform.getPlatformByName('CPU')  # finiteness only
    for key in (0.4, 'product'):
        system = systems[key]
        gb = next(f for f in system.getForces()
                  if isinstance(f, mm.CustomGBForce))
        q_gb = [gb.getParticleParameters(i)[0]
                for i in range(gb.getNumParticles())]
        assert np.allclose(q_gb, _charges(system), atol=1e-12)
        context = mm.Context(system, mm.VerletIntegrator(0.001), platform)
        context.setPositions(np.array(builder.positions) * 0.1)
        energy = context.getState(getEnergy=True).getPotentialEnergy()
        assert np.isfinite(energy.value_in_unit(mmunit.kilojoule_per_mole))
        del context
