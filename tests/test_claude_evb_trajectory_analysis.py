"""Dev tests: per-term and per-residue energies recomputed from a trajectory
(EvbTrajectoryAnalysis, EvbDriver.compute_trajectory_analysis).

References are independent of the module where possible:
- geometry against OpenMM's distance/angle/dihedral functions;
- tracked bonded terms against their force-field formulas;
- residue bonded energies against copies of the systems with every term
  outside the residue switched off;
- in vacuum (no cutoff) the residue-rest total of the PME route equals the
  pairwise class sum exactly, and the residue-internal energy equals an OpenMM
  NonbondedForce holding only the residue;
- E_int at lambda 0 / 1 equals E1 / E2 for a residue away from the reaction,
  and the Coulomb energy with the fragment is linear in lambda;
- with PBC every column is invariant to moving a molecule by a box vector,
  and the Lennard-Jones total equals the pairwise sum without the dispersion
  correction.
"""

import copy
import json
from types import SimpleNamespace

import numpy as np
import pytest

mm = pytest.importorskip("openmm")
mmapp = pytest.importorskip("openmm.app")
mmunit = pytest.importorskip("openmm.unit")
pytest.importorskip("MDAnalysis")

from veloxchem.molecule import Molecule
from veloxchem.mmforcefieldgenerator import MMForceFieldGenerator
from veloxchem.reactionsystembuilder import ReactionSystemBuilder, EvbForceGroup
from veloxchem.evbdriver import EvbDriver
from veloxchem.evbfepdriver import EvbFepDriver
from veloxchem.evbtrajectoryanalysis import (EvbTrajectoryAnalysis,
                                             _parse_force_terms,
                                             _filtered_force)
from veloxchem.errorhandler import VeloxChemError
from veloxchem.outputstream import OutputStream

from test_evb_helper import EvbTestHelper

LAMBDAS = [0.0, 0.5, 1.0]
FG_MAP = {fg.name: fg.value for fg in EvbForceGroup}
TERMS = {
    "nuc": ["A:124:OD2", "B:1:C1"],
    "leave": ["B:1:C1", "B:1:CL1"],
    "ang": ["A:124:CG", "A:124:OD2", "B:1:C1"],
    "tor": ["A:124:N", "A:124:CA", "A:124:CB", "A:124:CG"],
    "cacb": ["A:124:CA", "A:124:CB"],
}
RESIDUES = ["A:124", "A:289", "B:1"]


def _platform():
    return EvbTestHelper.evb_platform() or "CPU"


def _dhla_states():
    data = EvbTestHelper.evb_data_dir()
    states = []
    for name in ("reactant", "product"):
        ff = MMForceFieldGenerator.load_forcefield_from_json_file(
            str(data / f"evb_dhla_{name}_ff_data.json"),
            ostream=OutputStream(None))
        ff.molecule = Molecule.read_xyz_file(str(data / f"evb_dhla_{name}.xyz"))
        states.append(ff)
    return states


def _dhla_config(name="dhla_vac"):
    return {
        "name": name,
        "temperature": 300.0,
        "pdb": str(EvbTestHelper.evb_data_dir() / "evb_dhla_2dhd_reactant.cif"),
        "pdb_active_res": [{
            "chain": "A",
            "residue": 124,
            "cut": "beta",
            "reacting_atoms": ["OD2"]
        }, {
            "chain": "B",
            "residue": 1,
            "reacting_atoms": ["C1", "CL1"]
        }],
    }


@pytest.fixture(scope="module")
def dhla():
    reactant, product = _dhla_states()
    builder = ReactionSystemBuilder(ostream=OutputStream(None))
    with EvbTestHelper.evb_chdir_tmp():
        systems, topology, _ = builder.build_systems(reactant, product,
                                                     LAMBDAS, _dhla_config())
    positions = np.asarray(builder.positions) * 0.1
    return SimpleNamespace(systems=systems,
                           topology=topology,
                           positions=positions)


def _analysis(dhla, terms=TERMS, residues=RESIDUES, **kwargs):
    return EvbTrajectoryAnalysis(dhla.topology,
                                 dhla.systems["reactant"],
                                 FG_MAP,
                                 terms=terms,
                                 residues=residues,
                                 platform=_platform(),
                                 ostream=OutputStream(None),
                                 **kwargs)


def _frame(positions, seed, sigma=0.003):
    rng = np.random.default_rng(seed)
    return positions + rng.normal(0.0, sigma, positions.shape)


def _values(analysis, systems, positions, box=None):
    """All columns of one frame; systems maps Hamiltonian -> system."""
    geo = analysis.frame_geometry(positions, box)
    values = dict(geo["values"])
    for ham, system in systems.items():
        ev = analysis.evaluator(system, ham)
        ev.evaluate(geo, positions, box, values)
        ev.close()
    return values


def _context(system, platform=None):
    return mm.Context(system, mm.VerletIntegrator(0.001),
                      mm.Platform.getPlatformByName(platform or _platform()))


def _energy(system, positions, box=None, groups=-1, platform=None):
    context = _context(system, platform)
    if box is not None:
        context.setPeriodicBoxVectors(*[mm.Vec3(*v) for v in box])
    context.setPositions(positions)
    return context.getState(getEnergy=True, groups=groups).getPotentialEnergy(
    ).value_in_unit(mmunit.kilojoule_per_mole)


def _index(topology, chain, resid, name):
    return next(a.index for a in topology.atoms()
                if a.residue.chain.id == chain and a.residue.id == str(resid)
                and a.name == name)


# ---------------------------------------------------------------------- #
# Selections and geometry
# ---------------------------------------------------------------------- #


def test_selections(dhla):
    an = _analysis(dhla)
    od2 = _index(dhla.topology, "A", 124, "OD2")
    assert an._resolve_atom("A:124:OD2") == od2
    assert an._resolve_atom(("A", 124, "OD2")) == od2
    assert an._resolve_atom(od2) == od2
    # Asp124 is cut (CB-CG): tracked as its reacting part and the rest
    assert [r["label"] for r in an.residues] == [
        "A:124:ASP[fragment]", "A:124:ASP[rest]", "A:289:HIS", "B:1:DCE"
    ]
    assert [r["kind"] for r in an.residues] == [
        "fragment part", "rest part", "residue", "residue"
    ]
    asp = next(r for r in dhla.topology.residues() if r.id == "124")
    asp_atoms = {a.name: a.index for a in asp.atoms()}
    part, rest = an.residues[0]["atoms"], an.residues[1]["atoms"]
    assert sorted(list(part) + list(rest)) == sorted(asp_atoms.values())
    assert sorted(part) == sorted(asp_atoms[n] for n in ("CG", "OD1", "OD2"))
    assert asp_atoms["CB"] in rest  # the cut's host atom
    assert [c["label"] for c in an.cuts] == ["A:124:ASP[cut]"]
    assert list(an.cuts[0]["atoms"]) == sorted(asp_atoms.values())
    # fragment: the carboxylate and DCE, in their own PDB residues
    names = {(a.residue.id, a.name) for a in dhla.topology.atoms()
             if a.index in set(an.fragment_atoms)}
    assert {("124", "OD2"), ("124", "CG"), ("1", "C1")} <= names
    assert ("124", "CA") not in names
    # one fragment class per residue holding reacting atoms
    assert an.class_names == ("protein", "fragment:A:124:ASP",
                              "fragment:B:1:DCE", "water", "ions", "other")
    for k, residue in ((1, an.residues[0]), (2, an.residues[3])):
        assert set(np.where(an.atom_classes == k)[0]) == set(
            residue["atoms"].tolist())
    assert an.atom_classes[asp_atoms["CB"]] == 0
    # a residue wholly inside the fragment, and a group, are not split
    whole = _analysis(dhla, terms=None,
                      residues=["B:1", {"label": "asp",
                                        "atoms": list(asp_atoms.values())}])
    assert [r["label"] for r in whole.residues] == ["B:1:DCE", "asp"]
    assert whole.cuts == []
    with pytest.raises(VeloxChemError, match="reserved"):
        _analysis(dhla, terms=None,
                  residues=[{"label": "classes", "atoms": ["B:1:C1"]}])
    with pytest.raises(VeloxChemError, match="matches 0 atoms"):
        _analysis(dhla, terms=[["A:124:XX", "B:1:C1"]], residues=None)
    with pytest.raises(VeloxChemError, match="2, 3 or 4 atoms"):
        _analysis(dhla, terms=[["A:124:OD2"]], residues=None)
    with pytest.raises(VeloxChemError, match="matches 0 residues"):
        _analysis(dhla, terms=None, residues=["A:9999"])


def test_geometry_against_openmm(dhla):
    an = _analysis(dhla, residues=None)
    x = _frame(dhla.positions, 1)
    values = an.frame_geometry(x, None)["values"]
    expressions = {
        2: "distance(p1,p2)",
        3: "angle(p1,p2,p3)",
        4: "dihedral(p1,p2,p3,p4)"
    }
    for term in an.terms:
        n = len(term["atoms"])
        system = mm.System()
        for _ in range(an.n_atoms):
            system.addParticle(1.0)
        force = mm.CustomCompoundBondForce(n, expressions[n])
        force.addBond(term["atoms"], [])
        system.addForce(force)
        ref = _energy(system, x)
        got = values[f"{term['label']}/{term['kind']}"]
        if n > 2:
            ref = np.degrees(ref)
        assert got == pytest.approx(ref, abs=1e-5)


# ---------------------------------------------------------------------- #
# Bonded terms
# ---------------------------------------------------------------------- #


def _numpy_dihedral(x):
    b1, b2, b3 = x[1] - x[0], x[2] - x[1], x[3] - x[2]
    n1, n2 = np.cross(b1, b2), np.cross(b2, b3)
    y = np.dot(np.cross(n1, n2), b2 / np.linalg.norm(b2))
    return np.arctan2(y, np.dot(n1, n2))


def test_tracked_terms_against_formulas(dhla):
    an = _analysis(dhla, residues=None)
    x = _frame(dhla.positions, 2)
    systems = {"E1": dhla.systems["reactant"], "E_int": dhla.systems[0.5]}
    values = _values(an, systems, x)
    ca, cb = (_index(dhla.topology, "A", 124, n) for n in ("CA", "CB"))
    tor_atoms = {_index(dhla.topology, "A", 124, n)
                 for n in ("N", "CA", "CB", "CG")}
    for ham, system in systems.items():
        bond = 0.0
        torsion = 0.0
        for force in system.getForces():
            if isinstance(force, mm.HarmonicBondForce):
                for k in range(force.getNumBonds()):
                    i, j, r0, kb = force.getBondParameters(k)
                    if {i, j} == {ca, cb}:
                        r = np.linalg.norm(x[i] - x[j])
                        bond += 0.5 * kb._value * (r - r0._value)**2
            if isinstance(force, mm.PeriodicTorsionForce):
                for k in range(force.getNumTorsions()):
                    *atoms, n, phase, kt = force.getTorsionParameters(k)
                    if set(atoms) == tor_atoms:
                        theta = _numpy_dihedral(x[atoms])
                        torsion += kt._value * (1 + np.cos(n * theta -
                                                           phase._value))
        assert bond > 0.0 and torsion != 0.0
        assert values[f"cacb/{ham}/bonded"] == pytest.approx(bond, rel=1e-6)
        assert values[f"tor/{ham}/bonded"] == pytest.approx(torsion, rel=1e-6)


def test_tracked_pair_matches_soft_core_forces(dhla):
    an = _analysis(dhla, residues=None)
    names = {}
    for ham, key in (("E1", "reactant"), ("E2", "product"), ("E_int", 0.5)):
        ev = an.evaluator(dhla.systems[key], ham)
        names[ham] = {
            role: sorted({m["force"] for m in matched})
            for role in ("bonded", "coul", "lj")
            for matched in [ev.matched.get(f"nuc/{ham}/{role}", [])]
        }
        ev.close()
    # the forming bond: Morse in the product PES only, soft-core
    # nonbonded in the reactant PES, a harmonic surrogate in E_int
    assert names["E1"]["bonded"] == []
    assert names["E1"]["coul"] == ["Reaction internal Coul dynamic"]
    assert names["E1"]["lj"] == ["Reaction internal LJ dynamic"]
    assert names["E2"]["bonded"] == ["Reaction morse bond"]
    assert names["E_int"]["bonded"] == ["Dynamic reaction harmonic bond"]


def test_filtered_forces_partition(dhla):
    """Filtered copies of every force add up to the force (the machinery
    behind every bonded and soft-core column)."""
    x = _frame(dhla.positions, 3)
    rng = np.random.default_rng(0)
    for key in ("reactant", 0.5):
        system = dhla.systems[key]
        for force in system.getForces():
            parsed = _parse_force_terms(force)
            if parsed is None or not parsed["atoms"]:
                continue
            # the force in group 31, three filtered parts in groups 0-2
            sub = mm.System()
            for i in range(system.getNumParticles()):
                sub.addParticle(1.0)
            single = copy.deepcopy(force)
            single.setForceGroup(31)
            sub.addForce(single)
            labels = rng.integers(0, 3, len(parsed["atoms"]))
            for g in range(3):
                sub.addForce(
                    _filtered_force(parsed, np.where(labels == g)[0], g))
            context = _context(sub)
            context.setPositions(x)

            def group_energy(g):
                state = context.getState(getEnergy=True, groups=1 << g)
                return state.getPotentialEnergy().value_in_unit(
                    mmunit.kilojoule_per_mole)

            ref = group_energy(31)
            total = sum(group_energy(g) for g in range(3))
            assert total == pytest.approx(ref, rel=1e-6, abs=1e-6), (
                key, force.getName())


def _zeroed_bonded(system, inside):
    """Copy of the bonded forces with every term not entirely inside
    'inside' switched off; groups: 1 bonds, 2 angles, 3 torsions."""
    out = copy.deepcopy(system)
    names = {v: k for k, v in FG_MAP.items()}
    for k in reversed(range(out.getNumForces())):
        force = out.getForce(k)
        group = names.get(force.getForceGroup(), "")
        keep_type = isinstance(force,
                               (mm.HarmonicBondForce, mm.HarmonicAngleForce,
                                mm.PeriodicTorsionForce))
        morse = group == "REA_MORSE_BOND"
        if not (keep_type or morse) or group == "CONSTRAINT":
            out.removeForce(k)
            continue
        if isinstance(force, mm.HarmonicBondForce) or morse:
            force.setForceGroup(1)
            for t in range(force.getNumBonds()):
                if isinstance(force, mm.HarmonicBondForce):
                    i, j, r0, kb = force.getBondParameters(t)
                    if not {i, j} <= inside:
                        force.setBondParameters(t, i, j, r0, 0.0)
                else:
                    i, j, (d, a, re) = force.getBondParameters(t)
                    if not {i, j} <= inside:
                        force.setBondParameters(t, i, j, [0.0, a, re])
        elif isinstance(force, mm.HarmonicAngleForce):
            force.setForceGroup(2)
            for t in range(force.getNumAngles()):
                i, j, l, th, ka = force.getAngleParameters(t)
                if not {i, j, l} <= inside:
                    force.setAngleParameters(t, i, j, l, th, 0.0)
        else:
            force.setForceGroup(3)
            for t in range(force.getNumTorsions()):
                i, j, l, m, n, ph, kt = force.getTorsionParameters(t)
                if not {i, j, l, m} <= inside:
                    force.setTorsionParameters(t, i, j, l, m, n, ph, 0.0)
    return out


def test_residue_bonded_against_zeroing(dhla):
    an = _analysis(dhla, terms=None, residues=["A:124", "B:1", "fragment"])
    x = _frame(dhla.positions, 4)
    systems = {"E1": dhla.systems["reactant"], "E_int": dhla.systems[0.5]}
    values = _values(an, systems, x)
    for residue in an.residues:
        inside = set(residue["atoms"].tolist())
        for ham, system in systems.items():
            zeroed = _zeroed_bonded(system, inside)
            base = f"{residue['label']}/{ham}/bonded"
            for group, category in ((1, "bond"), (2, "angle"), (3,
                                                                "torsion")):
                ref = _energy(zeroed, x, groups=1 << group)
                assert values[f"{base}/{category}"] == pytest.approx(
                    ref, rel=1e-6, abs=1e-6), (base, category)
            assert values[f"{base}/other"] == 0.0
    # The cut: terms inside Asp124 minus those inside either part
    cut = an.cuts[0]
    sides = [set(cut[k].tolist()) for k in ("atoms", "fragment", "rest")]
    for ham, system in systems.items():
        base = f"{cut['label']}/{ham}/bonded"
        for group, category in ((1, "bond"), (2, "angle"), (3, "torsion")):
            whole, frag, rest = (_energy(_zeroed_bonded(system, s), x,
                                         groups=1 << group) for s in sides)
            assert abs(values[f"{base}/{category}"]) > 1e-3, (base, category)
            assert values[f"{base}/{category}"] == pytest.approx(
                whole - frag - rest, rel=1e-6, abs=1e-6), (base, category)
        assert values[f"{base}/other"] == 0.0


# ---------------------------------------------------------------------- #
# Nonbonded terms in vacuum
# ---------------------------------------------------------------------- #


def test_vacuum_pme_route_equals_pairwise(dhla):
    # Reference platform: the CPU platform's single precision leaves the
    # difference of whole-system energies ~4e-3 kJ/mol off
    # the fragment overlaps Asp124: separate switch layers
    an = _analysis(dhla, terms=None, residues=["A:124", "A:289", "fragment"])
    an.platform = "Reference"
    x = _frame(dhla.positions, 5)
    systems = {"E1": dhla.systems["reactant"]}
    values = _values(an, systems, x)
    for residue in an.residues:
        for ham in systems:
            base = f"{residue['label']}/{ham}/nb_rest"
            for kind in ("coul", "lj"):
                pairwise = sum(values[f"{base}/{c}/{kind}"]
                               for c in an.class_names)
                assert values[f"{base}/{kind}"] == pytest.approx(
                    pairwise, abs=1e-6), (base, kind)
    # His289 holds no fragment atom: its internal energy is the
    # NonbondedForce of the residue alone
    his = an.residues[2]
    assert his["label"] == "A:289:HIS"
    nb = copy.deepcopy(next(f for f in dhla.systems["reactant"].getForces()
                            if isinstance(f, mm.NonbondedForce)))
    inside = set(his["atoms"].tolist())
    for i in range(nb.getNumParticles()):
        if i not in inside:
            q, s, e = nb.getParticleParameters(i)
            nb.setParticleParameters(i, 0.0, s, 0.0)
    for k in range(nb.getNumExceptions()):
        i, j, qq, s, e = nb.getExceptionParameters(k)
        if not {i, j} <= inside:
            nb.setExceptionParameters(k, i, j, 0.0, s, 0.0)
    alone = mm.System()
    for i in range(an.n_atoms):
        alone.addParticle(1.0)
    alone.addForce(nb)
    internal = (values["A:289:HIS/E1/nb_internal/coul"] +
                values["A:289:HIS/E1/nb_internal/lj"])
    assert internal == pytest.approx(_energy(alone, x), abs=1e-3)
    # The rest of Asp124 sees its own reacting part and the substrate as
    # separate fragment classes (and nothing is water)
    for cls in ("fragment:A:124:ASP", "fragment:B:1:DCE"):
        assert values[f"A:124:ASP[rest]/E1/nb_rest/{cls}/coul"] != 0.0
    assert values["A:124:ASP[fragment]/E1/nb_rest/fragment:A:124:ASP/coul"] \
        == 0.0
    assert values["A:124:ASP[rest]/E1/nb_rest/water/coul"] == 0.0


def _nonbonded_only(system):
    out = mm.System()
    for i in range(system.getNumParticles()):
        out.addParticle(system.getParticleMass(i))
    out.setDefaultPeriodicBoxVectors(*system.getDefaultPeriodicBoxVectors())
    nb = next(f for f in system.getForces() if isinstance(f, mm.NonbondedForce))
    out.addForce(copy.deepcopy(nb))
    return out


def _class_pair_sum(an, values, ham, kind):
    return sum(values[f"classes/{ham}/{a}~{b}/{kind}"]
               for a, b in an.class_pair_list())


def test_vacuum_class_pairs(dhla):
    """Class pairs sum to the NonbondedForce energy; without a cutoff a
    class pair equals the pairwise class part of a residue holding the whole
    class (the substrate), and the parts of the split Asp124."""
    an = _analysis(dhla, terms=None,
                   residues=["A:124", "B:1"])
    an.platform = "Reference"
    x = _frame(dhla.positions, 8)
    systems = {"E1": dhla.systems["reactant"], "E_int": dhla.systems[0.5]}
    values = _values(an, systems, x)
    assert len(an.class_pair_list()) == 21
    for ham, system in systems.items():
        total = (_class_pair_sum(an, values, ham, "coul") +
                 _class_pair_sum(an, values, ham, "lj"))
        assert total == pytest.approx(_energy(_nonbonded_only(system), x,
                                              platform="Reference"),
                                      abs=1e-6)
        pairs = {
            "protein~fragment:B:1:DCE": "B:1:DCE/{}/nb_rest/protein",
            "protein~fragment:A:124:ASP":
                "A:124:ASP[fragment]/{}/nb_rest/protein",
        }
        for pair, column in pairs.items():
            for kind in ("coul", "lj"):
                ref = values[f"{column.format(ham)}/{kind}"]
                assert abs(ref) > 1e-2, (pair, kind)
                assert values[f"classes/{ham}/{pair}/{kind}"] == \
                    pytest.approx(ref, abs=1e-6), (ham, pair, kind)
        # no water or ions in this system; the fragment pair sits in the
        # soft-core forces, not the NonbondedForce
        for pair in ("protein~water", "water~water", "ions~ions"):
            assert values[f"classes/{ham}/{pair}/coul"] == 0.0
        assert values[f"classes/{ham}/fragment:A:124:ASP~fragment:B:1:DCE/"
                      "coul"] == pytest.approx(0.0, abs=1e-9)
    off = _analysis(dhla, terms=None, residues=["B:1"], class_pairs=False)
    assert not any(c.startswith("classes/") for c in off.columns)


def test_atom_group_equals_its_residue(dhla):
    dce = ["B:1:" + n for n in ("C1", "C2", "CL1", "CL2", "H11", "H12", "H21",
                                "H22")]
    an = _analysis(dhla,
                   terms=None,
                   residues=["B:1", {"label": "dce", "atoms": dce}])
    assert list(an.residues[0]["atoms"]) == list(an.residues[1]["atoms"])
    frag = _analysis(dhla, terms=None, residues=["fragment"])
    assert list(frag.residues[0]["atoms"]) == list(frag.fragment_atoms)
    values = _values(an, {"E1": dhla.systems["reactant"]},
                     _frame(dhla.positions, 7))
    for column in an.columns:
        if column.startswith("B:1:DCE/E1/"):
            other = "dce/" + column[len("B:1:DCE/"):]
            assert values[column] == pytest.approx(values[other],
                                                   abs=5e-3), column


def test_e_int_endpoints_and_linear_charges(dhla):
    an = _analysis(dhla, terms={"cacb": TERMS["cacb"]}, residues=["A:289"])
    x = _frame(dhla.positions, 6)
    v0 = _values(an, {"E1": dhla.systems["reactant"],
                      "E_int": dhla.systems[0.0]}, x)
    v1 = _values(an, {"E2": dhla.systems["product"],
                      "E_int": dhla.systems[1.0]}, x)
    vh = _values(an, {"E_int": dhla.systems[0.5]}, x)
    for column in an.columns:
        if "/E1/" in column:
            assert v0[column] == pytest.approx(
                v0[column.replace("/E1/", "/E_int/")], rel=1e-6,
                abs=1e-4), column
        if "/E2/" in column:
            assert v1[column] == pytest.approx(
                v1[column.replace("/E2/", "/E_int/")], rel=1e-6,
                abs=1e-4), column
    # the fragment charges are linear in lambda, so is His289's Coulomb
    # energy with the fragment
    col = "A:289:HIS/E_int/nb_rest/fragment:B:1:DCE/coul"
    assert v0[col] != pytest.approx(v1[col], abs=1e-3)
    assert vh[col] == pytest.approx(0.5 * (v0[col] + v1[col]), abs=1e-6)


# ---------------------------------------------------------------------- #
# Periodic system
# ---------------------------------------------------------------------- #


@pytest.fixture(scope="module")
def water():
    """Solvated ethanol, the lambda 0 system minimized. Flexible cSPC/E:
    minimizing rigid TIP3P takes ~12 s."""
    pair = EvbTestHelper.evb_ff_pair()
    evb = EvbDriver(ostream=OutputStream(None))
    evb.water_model = "cspce"
    EvbTestHelper.evb_seed_states(evb, [pair.reactant, pair.product])
    with EvbTestHelper.evb_chdir_tmp():
        conf = evb.default_system_configurations("water_NVT")
        conf.update({"name": "water", "padding": 1.0})
        evb.build_systems(configurations=[conf], Lambda=LAMBDAS)
        conf = evb.system_confs[0]
        systems, topology = conf["systems"], conf["topology"]
        positions = np.asarray(conf["initial_positions"]) * 0.1
    context = _context(systems[0.0])
    context.setPositions(positions)
    mm.LocalEnergyMinimizer.minimize(context, 10.0, 300)
    state = context.getState(getPositions=True)
    positions = state.getPositions(asNumpy=True).value_in_unit(
        mmunit.nanometer)
    box = np.array(state.getPeriodicBoxVectors().value_in_unit(
        mmunit.nanometer))
    return SimpleNamespace(systems=systems,
                           topology=topology,
                           positions=positions,
                           box=box)


def _water_analysis(water, systems=None):
    solute = next(r for r in water.topology.residues()
                  if r.name not in ("HOH", "SOL", "ION"))
    first = next(solute.atoms())
    # the water whose oxygen is nearest the solute's first atom
    oxygens = [a for a in water.topology.atoms()
               if a.residue.name in ("HOH", "SOL") and a.element.symbol == "O"]
    d = water.positions[[a.index for a in oxygens]] - water.positions[
        first.index]
    d -= np.round(d / np.diag(water.box)) * np.diag(water.box)
    o = oxygens[int(np.argmin(np.linalg.norm(d, axis=1)))]
    waters = [o.residue]
    an = EvbTrajectoryAnalysis(
        water.topology, (systems or water.systems)["reactant"],
        FG_MAP,
        terms={"pair": [first.index, o.index]},
        residues=[(solute.chain.id, solute.id),
                  (waters[0].chain.id, waters[0].id)],
        hamiltonians=("E1", ),
        platform=_platform(),
        ostream=OutputStream(None))
    return an, solute, waters[0]


def test_periodic_image_invariance(water):
    an, solute, wat = _water_analysis(water)
    systems = {"E1": water.systems["reactant"]}
    ref = _values(an, systems, water.positions, water.box)
    moved = water.positions.copy()
    for atom in wat.atoms():
        moved[atom.index] += water.box[0] - water.box[1]
    for atom in solute.atoms():
        moved[atom.index] -= water.box[2]
    got = _values(an, systems, moved, water.box)
    assert ref["pair/distance"] < 0.5 * water.box[0, 0]
    for column in an.columns:
        assert got[column] == pytest.approx(ref[column], abs=5e-3), column
    # a solvated residue interacts with water and the fragment
    assert ref[f"{an.residues[1]['label']}/E1/nb_rest/water/coul"] != 0.0


def test_periodic_lj_route_equals_pairwise_without_dispersion(water):
    systems = {}
    for key in ("reactant", ):
        system = copy.deepcopy(water.systems[key])
        for force in system.getForces():
            if isinstance(force, mm.NonbondedForce):
                force.setUseDispersionCorrection(False)
        systems[key] = system
    an, _, _ = _water_analysis(water, systems)
    values = _values(an, {"E1": systems["reactant"]}, water.positions,
                     water.box)
    for residue in an.residues:
        base = f"{residue['label']}/E1/nb_rest"
        pairwise = sum(values[f"{base}/{c}/lj"] for c in an.class_names)
        assert values[f"{base}/lj"] == pytest.approx(pairwise, abs=2e-3)


def _dispersion_energy(system, positions, box, eps_on):
    """OpenMM's dispersion correction of the Lennard-Jones-only force with
    epsilon kept where eps_on: E(with) - E(without)."""
    nb = copy.deepcopy(next(f for f in system.getForces()
                            if isinstance(f, mm.NonbondedForce)))
    for i in range(nb.getNumParticles()):
        q, sig, eps = nb.getParticleParameters(i)
        nb.setParticleParameters(i, 0.0, sig, eps if eps_on[i] else 0.0)
    for k in range(nb.getNumExceptions()):
        i, j, qq, sig, eps = nb.getExceptionParameters(k)
        nb.setExceptionParameters(k, i, j, 0.0, sig, 0.0)
    energies = []
    for flag in (True, False):
        nb.setUseDispersionCorrection(flag)
        sub = mm.System()
        for i in range(system.getNumParticles()):
            sub.addParticle(1.0)
        sub.setDefaultPeriodicBoxVectors(*system.getDefaultPeriodicBoxVectors())
        sub.addForce(copy.deepcopy(nb))
        energies.append(_energy(sub, positions, box))
    return energies[0] - energies[1]


def test_periodic_lj_route_has_the_dispersion_cross_term(water):
    """With the dispersion correction on, the Lennard-Jones total differs
    from the pairwise sum by exactly the residue-rest part of OpenMM's
    correction, all - residue off - residue alone."""
    an, _, _ = _water_analysis(water)
    system = water.systems["reactant"]
    nb = next(f for f in system.getForces() if isinstance(f, mm.NonbondedForce))
    assert nb.getUseDispersionCorrection()
    values = _values(an, {"E1": system}, water.positions, water.box)
    for residue in an.residues:
        inside = np.zeros(an.n_atoms, dtype=bool)
        inside[residue["atoms"]] = True
        reference = (
            _dispersion_energy(system, water.positions, water.box,
                               np.ones(an.n_atoms, dtype=bool)) -
            _dispersion_energy(system, water.positions, water.box, ~inside) -
            _dispersion_energy(system, water.positions, water.box, inside))
        base = f"{residue['label']}/E1/nb_rest"
        pairwise = sum(values[f"{base}/{c}/lj"] for c in an.class_names)
        assert abs(reference) > 1e-3
        assert values[f"{base}/lj"] - pairwise == pytest.approx(reference,
                                                                abs=2e-3)


def test_periodic_coulomb_route_is_linear_in_residue_charges(water):
    """The PME-route interaction is bilinear in the charges: doubling the
    water's charges doubles it. Self-energy or net-charge terms leaking in
    would scale with the square."""
    an, _, wat = _water_analysis(water)
    inside = {a.index for a in wat.atoms()}
    scaled = copy.deepcopy(water.systems["reactant"])
    nb = next(f for f in scaled.getForces() if isinstance(f, mm.NonbondedForce))
    for i in inside:
        q, sig, eps = nb.getParticleParameters(i)
        nb.setParticleParameters(i, 2 * q, sig, eps)
    for k in range(nb.getNumExceptions()):
        i, j, qq, sig, eps = nb.getExceptionParameters(k)
        factor = 2 ** len({i, j} & inside)
        nb.setExceptionParameters(k, i, j, factor * qq, sig, eps)
    base = f"{an.residues[1]['label']}/E1/nb_rest"
    single = _values(an, {"E1": water.systems["reactant"]}, water.positions,
                     water.box)
    double = _values(an, {"E1": scaled}, water.positions, water.box)
    assert abs(single[f"{base}/coul"]) > 1.0
    assert double[f"{base}/coul"] == pytest.approx(2 * single[f"{base}/coul"],
                                                   abs=5e-3)
    for cls in an.class_names:
        assert double[f"{base}/{cls}/coul"] == pytest.approx(
            2 * single[f"{base}/{cls}/coul"], abs=1e-9)
    assert double[f"{base}/lj"] == pytest.approx(single[f"{base}/lj"],
                                                 abs=5e-3)


def test_periodic_class_pairs(water):
    """With PME and the dispersion correction, the class pairs sum to
    OpenMM's NonbondedForce energy (Coulomb and Lennard-Jones separately),
    and the solute-water pair equals the solute's nb_rest (the system holds
    nothing else). Doubling the water charges doubles the solute-water
    Coulomb pair and quadruples the water-water one."""
    an, solute, _ = _water_analysis(water)
    system = water.systems["reactant"]
    frag = f"fragment:{solute.chain.id}:{solute.id}:{solute.name}"
    assert an.class_names == ("protein", frag, "water", "ions", "other")
    assert set(np.where(an.atom_classes == 1)[0]) == {
        a.index for a in solute.atoms()}
    values = _values(an, {"E1": system}, water.positions, water.box)

    only = _nonbonded_only(system)
    nb = only.getForce(0)
    assert nb.getUseDispersionCorrection()
    total = _energy(only, water.positions, water.box)
    for i in range(nb.getNumParticles()):
        q, sig, eps = nb.getParticleParameters(i)
        nb.setParticleParameters(i, q, sig, 0.0)
    for k in range(nb.getNumExceptions()):
        i, j, qq, sig, eps = nb.getExceptionParameters(k)
        nb.setExceptionParameters(k, i, j, qq, sig, 0.0)
    coul = _energy(only, water.positions, water.box)
    assert _class_pair_sum(an, values, "E1", "coul") == pytest.approx(
        coul, abs=2e-2)
    assert _class_pair_sum(an, values, "E1", "lj") == pytest.approx(
        total - coul, abs=2e-2)

    pair = f"classes/E1/{frag}~water"
    label = f"{solute.chain.id}:{solute.id}:{solute.name}"
    for kind in ("coul", "lj"):
        assert abs(values[f"{pair}/{kind}"]) > 1.0
        assert values[f"{pair}/{kind}"] == pytest.approx(
            values[f"{label}/E1/nb_rest/{kind}"], abs=5e-3)

    scaled = copy.deepcopy(system)
    nb = next(f for f in scaled.getForces() if isinstance(f, mm.NonbondedForce))
    waters = set(np.where(an.atom_classes == 2)[0].tolist())
    for i in waters:
        q, sig, eps = nb.getParticleParameters(i)
        nb.setParticleParameters(i, 2 * q, sig, eps)
    for k in range(nb.getNumExceptions()):
        i, j, qq, sig, eps = nb.getExceptionParameters(k)
        nb.setExceptionParameters(k, i, j, 2**len({i, j} & waters) * qq, sig,
                                  eps)
    double = _values(an, {"E1": scaled}, water.positions, water.box)
    assert double[f"{pair}/coul"] == pytest.approx(
        2 * values[f"{pair}/coul"], abs=5e-3)
    ww = "classes/E1/water~water/coul"
    assert double[ww] == pytest.approx(4 * values[ww], rel=1e-5)


# ---------------------------------------------------------------------- #
# Driver
# ---------------------------------------------------------------------- #


def _write_run(folder, topology, frames, lam, replica, direction):
    """Trajectory, Energies.csv and Data_combined.csv of a fake run."""
    with (folder / "trajectory.dcd").open("wb") as f:
        dcd = mmapp.DCDFile(f, topology, 0.001)
        for x in frames:
            dcd.writeModel(x * mmunit.nanometer)
    n = len(frames)
    rows = np.column_stack([lam] + [np.zeros(n)] * 5 + [replica, direction])
    np.savetxt(folder / "Energies.csv",
               rows,
               delimiter=", ",
               comments="",
               header="Lambda, reactant PES, product PES, reactant "
               "integration, product integration, Em, replica, direction")
    np.savetxt(folder / "Data_combined.csv",
               np.zeros((n, 6)),
               delimiter=", ",
               comments="",
               header="step, Ep, Ek, Temp, Vol, Dens")


def test_driver_end_to_end(dhla, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    reactant, product = _dhla_states()
    evb = EvbDriver(ostream=OutputStream(None))
    EvbTestHelper.evb_seed_states(evb, [reactant, product])
    evb.build_systems(configurations=[_dhla_config("dhla")], Lambda=LAMBDAS)
    folder = tmp_path / evb.system_confs[0]["data_folder"]

    lam = np.repeat(LAMBDAS, 2)
    lam = np.concatenate([lam, lam[::-1]])
    direction = np.repeat([0, 1], 6)
    replica = np.zeros(12, dtype=int)
    frames = [_frame(dhla.positions, 10 + k, sigma=0.002) for k in range(12)]
    _write_run(folder, dhla.topology, frames, lam, replica, direction)

    loaded = EvbDriver(ostream=OutputStream(None))
    loaded.load_initialisation(str(folder), load_systems=True, load_top=True)
    terms = {"nuc": TERMS["nuc"], "cacb": TERMS["cacb"]}
    loaded.compute_trajectory_analysis(terms=terms,
                                       frame_stride=2,
                                       platform=_platform())

    data = np.loadtxt(folder / "TrajectoryAnalysis.csv",
                      skiprows=1,
                      delimiter=",")
    with (folder / "TrajectoryAnalysis.csv").open() as f:
        header = [c.strip() for c in f.readline().split(",")]
    meta = json.loads((folder / "TrajectoryAnalysis.json").read_text())
    assert list(data[:, 0].astype(int)) == list(range(0, 12, 2))
    assert header[4:] == meta["columns"]
    assert meta["n_frames"] == 12 and meta["frame_stride"] == 2
    assert sorted(meta["matched_terms"]["E_int"]) == ["0.000", "0.500", "1.000"]

    # each row equals a direct evaluation with the frame's own lambda system
    conf = loaded.system_confs[0]
    an = EvbTrajectoryAnalysis(conf["topology"],
                               conf["systems"]["reactant"],
                               json.loads((folder / "run" /
                                           "force_group_name.json").read_text()),
                               terms=terms,
                               platform=_platform(),
                               ostream=OutputStream(None))
    with EvbFepDriver._open_trajectory(folder / "trajectory.dcd",
                                       False) as (_, read):
        read_back = [pos for pos, _ in read(0, 12)]
    fixed = [an.evaluator(conf["systems"]["reactant"], "E1"),
             an.evaluator(conf["systems"]["product"], "E2")]
    per_lambda = {L: an.evaluator(conf["systems"][L], "E_int")
                  for L in LAMBDAS}
    for row in data:
        frame = int(row[0])
        x = read_back[frame]
        geo = an.frame_geometry(x, None)
        ref = dict(geo["values"])
        for ev in fixed + [per_lambda[float(lam[frame])]]:
            ev.evaluate(geo, x, None, ref)
        np.testing.assert_allclose(row[4:], [ref[c] for c in an.columns],
                                   rtol=1e-7, atol=1e-6)

    # loading: aligned with the kept frames, NaN where not analysed
    results = loaded._load_output_from_folders(1, False, 1)
    analysis = results["configuration_results"]["dhla"]["trajectory_analysis"]
    values = analysis["values"]
    assert values.shape == (12, len(an.columns))
    assert np.isnan(values[1::2]).all() and np.isfinite(values[0::2]).all()
    col = EvbTrajectoryAnalysis.column(analysis, "cacb/distance")
    np.testing.assert_allclose(col[0::2], data[:, 4 + an.columns.index(
        "cacb/distance")])
    averages = EvbTrajectoryAnalysis.lambda_averages(
        results["configuration_results"]["dhla"])
    assert list(averages["Lambda"]) == LAMBDAS
    k = an.columns.index("cacb/distance")
    assert list(averages["n"][:, k]) == [2, 2, 2]
    assert np.isfinite(averages["sem"][:, k]).all()

    # compute_force_groups on the reloaded run: options.json gives the
    # forming/breaking bonds back as per-step lists of lists (defect E5-1,
    # TypeError: unhashable type: 'list', in the bonded sanity check)
    loaded.compute_force_groups(platform=_platform(), decompose_bonded=True)
    assert (folder / "bonded_E1_decomp.csv").is_file()

    # a trajectory that does not match Energies.csv is refused
    _write_run(folder, dhla.topology, frames[:10], lam[:10], replica[:10],
               direction[:10])
    np.savetxt(folder / "Energies.csv", np.zeros((12, 8)), delimiter=", ",
               comments="", header="Lambda, reactant PES, product PES, "
               "reactant integration, product integration, Em, replica, "
               "direction")
    with pytest.raises(VeloxChemError, match="row-aligned"):
        loaded.compute_trajectory_analysis(terms=terms, platform=_platform())
