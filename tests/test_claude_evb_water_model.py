"""Water model choice, rigid TIP3P and amber14 counter-ions in the EVB systems.

Dev-only (test_claude_*): never to master.

- EvbDriver.water_model None (default) chooses per configuration: TIP3P for a
  protein (pdb) configuration, the configuration's water solvent otherwise
  (the presets name cSPC/E). A configuration's 'water_model' key, then the
  driver's water_model, override it.
- ReactionSystemBuilder warns when a protein is solvated in water other than
  TIP3P (amber14 is parameterized with TIP3P).
- 'tip3p' and 'ctip3p' water is rigid: O-H and H-H constraints, no harmonic
  water terms, in every lambda and PES system. The references are OpenMM's own
  amber14/tip3p.xml with rigidWater=True.
- With 'tip3p', counter-ions take amber14/tip3p.xml parameters (not UFF).
- constrain_H equilibration removes only the H constraints it added (unit
  tests; constrain_H itself NaNs or stalls at lambda 0.5 on the solvated
  ethanol fixture with either water model, before this change too).
"""

import math
import os

import numpy as np
import pytest

from veloxchem.evbdriver import EvbDriver
from veloxchem.evbfepdriver import EvbFepDriver
from veloxchem.reactionsystembuilder import ReactionSystemBuilder
from veloxchem.outputstream import OutputStream
from veloxchem.solvationbuilder import SolvationBuilder
from veloxchem.waterparameters import get_water_parameters

from test_evb_helper import EvbTestHelper
from test_claude_evb_pdb_neutralize import CASES, _load, _fast_volume

mm = pytest.importorskip("openmm")
import openmm.app as mmapp  # noqa: E402
import openmm.unit as mmunit  # noqa: E402

LAMBDAS = [0.0, 0.5, 1.0]
WARNING = "parameterized with TIP3P"


def _driver(water_model=None):
    evb = EvbDriver(ostream=OutputStream(None))
    evb.water_model = water_model
    return evb


# ---------------------------------------------------------------- resolution


def test_auto_non_protein_keeps_the_preset():
    evb = _driver()
    conf = evb.default_system_configurations("water")
    assert (conf["name"], conf["solvent"]) == ("water_cspce_NVT", "cspce")
    assert evb._configuration_water_model(conf) == "cspce"
    assert conf["name"] == "water_cspce_NVT"
    assert conf["water_model"] == "cspce"


@pytest.mark.parametrize("preset, name", [("water", "water_tip3p_NVT"),
                                          ("water_NPT", "water_tip3p_NPT"),
                                          ("E_field", "water_E_field_tip3p")])
def test_auto_protein_uses_tip3p(preset, name):
    evb = _driver()
    conf = evb.default_system_configurations(preset)
    conf["pdb"] = "protein.cif"
    assert evb._configuration_water_model(conf) == "tip3p"
    assert (conf["name"], conf["solvent"], conf["water_model"]) == (name,
                                                                    "tip3p",
                                                                    "tip3p")


def test_auto_protein_keeps_a_custom_name():
    evb = _driver()
    conf = evb.default_system_configurations("water")
    conf.update({"name": "PRO", "pdb": "protein.cif"})
    assert evb._configuration_water_model(conf) == "tip3p"
    assert conf["name"] == "PRO"


def test_water_solvent_names_the_model_when_none_chosen():
    """Before, the solvent placed the molecules and the driver's model gave
    the parameters, so solvent 'tip3p' ran with cSPC/E parameters."""
    conf = {"name": "w", "solvent": "ctip3p"}
    assert _driver()._configuration_water_model(conf) == "ctip3p"


def test_explicit_choices_win():
    conf = {"name": "PRO", "solvent": "cspce", "pdb": "protein.cif"}
    assert _driver("cspce")._configuration_water_model(conf) == "cspce"
    conf = {"name": "PRO", "solvent": "cspce", "pdb": "protein.cif",
            "water_model": "spce"}
    assert _driver("tip3p")._configuration_water_model(conf) == "spce"
    assert conf["solvent"] == "spce"


def test_unknown_model_refused():
    with pytest.raises(Exception, match="unknown water model"):
        _driver("tip4p")._configuration_water_model({"solvent": "cspce"})


@pytest.mark.parametrize("conf", [{"name": "vacuum"},
                                  {"name": "eth", "solvent": "ethanol"}])
def test_other_configurations_untouched(conf):
    before = dict(conf)
    assert _driver()._configuration_water_model(conf) == "cspce"
    assert conf == before


@pytest.mark.parametrize("model, pdb, warns", [("cspce", "p.cif", True),
                                               ("ctip3p", "p.cif", True),
                                               ("tip3p", "p.cif", False),
                                               ("cspce", None, False)])
def test_builder_warns_for_protein_in_other_water(tmp_path, model, pdb, warns):
    out = tmp_path / "builder.out"
    builder = ReactionSystemBuilder(ostream=OutputStream(str(out)))
    builder.pdb = pdb
    builder.water_model = model
    assert builder._solvent_water_model("cspce") == model
    builder.ostream.flush()
    assert (WARNING in out.read_text()) == warns


# ------------------------------------------- references from OpenMM's amber14


def _amber14_tip3p():
    return mmapp.ForceField("amber14/tip3p.xml")


def test_rigid_distances_match_openmm():
    top = mmapp.Topology()
    res = top.addResidue("HOH", top.addChain())
    o = top.addAtom("O", mmapp.Element.getBySymbol("O"), res)
    h1 = top.addAtom("H1", mmapp.Element.getBySymbol("H"), res)
    h2 = top.addAtom("H2", mmapp.Element.getBySymbol("H"), res)
    top.addBond(o, h1)
    top.addBond(o, h2)
    system = _amber14_tip3p().createSystem(top,
                                           nonbondedMethod=mmapp.NoCutoff,
                                           rigidWater=True)
    reference = sorted(
        system.getConstraintParameters(i)[2].value_in_unit(mmunit.nanometer)
        for i in range(system.getNumConstraints()))
    r_oh, r_hh = ReactionSystemBuilder._rigid_water_distances(
        get_water_parameters()["tip3p"])
    assert np.allclose(sorted([r_oh, r_oh, r_hh]), reference, atol=1e-6)


@pytest.mark.parametrize("element, resname", [("Na", "NA"), ("Cl", "CL"),
                                              ("K", "K")])
def test_ion_parameters_match_openmm(element, resname):
    top = mmapp.Topology()
    res = top.addResidue(resname, top.addChain())
    top.addAtom(resname, mmapp.Element.getBySymbol(element), res)
    system = _amber14_tip3p().createSystem(top,
                                           nonbondedMethod=mmapp.NoCutoff)
    nb = next(f for f in system.getForces()
              if isinstance(f, mm.NonbondedForce))
    q, sigma, epsilon = nb.getParticleParameters(0)
    builder = ReactionSystemBuilder(ostream=OutputStream(None))
    assert builder._amber14_tip3p_ion(element) == pytest.approx(
        (q.value_in_unit(mmunit.elementary_charge),
         sigma.value_in_unit(mmunit.nanometer),
         epsilon.value_in_unit(mmunit.kilojoule_per_mole)), abs=1e-12)


# -------------------------------------------------------------------- builds


def _all_systems(systems):
    return [systems[key] for key in list(LAMBDAS) + ["reactant", "product"]]


def _waters(topology):
    return [[a.index for a in res.atoms()] for res in topology.residues()
            if res.name in ("SOL", "HOH") and len(list(res.atoms())) == 3]


def _constraints(system):
    return {
        tuple(sorted(system.getConstraintParameters(i)[:2])):
        system.getConstraintParameters(i)[2].value_in_unit(mmunit.nanometer)
        for i in range(system.getNumConstraints())
    }


def _water_force_names(system):
    return [f.getName() for f in system.getForces()
            if f.getName().startswith("Solvent") and (
                "harmonic bond" in f.getName()
                or "harmonic angle" in f.getName())]


@pytest.fixture(scope="module", params=["tip3p", "cspce"])
def ethanol_build(request, tmp_path_factory):
    pair = EvbTestHelper.evb_ff_pair()
    builder = ReactionSystemBuilder(ostream=OutputStream(None))
    config = {"name": "w", "temperature": 300.0, "solvent": request.param,
              "padding": 0.8}
    with EvbTestHelper.evb_chdir_tmp():
        systems, topology, _ = builder.build_systems(pair.reactant,
                                                     pair.product,
                                                     list(LAMBDAS), config)
    return request.param, systems, topology


def test_water_rigid_only_for_tip3p(ethanol_build):
    model, systems, topology = ethanol_build
    waters = _waters(topology)
    assert len(waters) > 50
    r_oh, r_hh = ReactionSystemBuilder._rigid_water_distances(
        get_water_parameters()[model])
    for system in _all_systems(systems):
        constraints = _constraints(system)
        if model == "tip3p":
            assert len(constraints) == 3 * len(waters)
            assert _water_force_names(system) == []
            for o, h1, h2 in waters:
                assert constraints[(o, h1)] == pytest.approx(r_oh)
                assert constraints[(o, h2)] == pytest.approx(r_oh)
                assert constraints[(h1, h2)] == pytest.approx(r_hh)
        else:
            assert constraints == {}
            assert len(_water_force_names(system)) == 2


@pytest.fixture(scope="module")
def dhla_tip3p(tmp_path_factory):
    """DhlA (Asp124 cut at beta), TIP3P with no water model set: the solvent
    names it; 17 Na+ from amber14/tip3p.xml."""
    prefix, cif, residues, pdb_charge = CASES["dhla_cut"]
    out = tmp_path_factory.mktemp("dhla_tip3p") / "build.out"
    config = {"name": "dhla", "temperature": 300.0, "pdb": str(EvbTestHelper.
                                                             evb_data_dir() /
                                                             cif),
              "pdb_active_res": residues, "solvent": "tip3p",
              "neutralize": True, "padding": 0.2}
    builder = ReactionSystemBuilder(ostream=OutputStream(str(out)))
    with pytest.MonkeyPatch.context() as patch:
        patch.setattr(SolvationBuilder, "_get_volume", _fast_volume)
        with EvbTestHelper.evb_chdir_tmp():
            systems, topology, _ = builder.build_systems(
                _load(prefix, "reactant"), _load(prefix, "product"),
                list(LAMBDAS), config)
    builder.ostream.flush()
    return systems, topology, pdb_charge, out.read_text()


def test_dhla_ions_from_amber14(dhla_tip3p):
    systems, topology, pdb_charge, text = dhla_tip3p
    ions = [a.index for a in topology.atoms() if a.residue.name == "ION"]
    assert len(ions) == abs(pdb_charge)
    q_ref, sigma_ref, eps_ref = ReactionSystemBuilder(
        ostream=OutputStream(None))._amber14_tip3p_ion("Na")
    for system in _all_systems(systems):
        nb = next(f for f in system.getForces()
                  if isinstance(f, mm.NonbondedForce))
        for i in ions:
            q, sigma, epsilon = nb.getParticleParameters(i)
            assert q.value_in_unit(mmunit.elementary_charge) == q_ref
            assert sigma.value_in_unit(mmunit.nanometer) == sigma_ref
            assert epsilon.value_in_unit(
                mmunit.kilojoule_per_mole) == eps_ref
    assert "amber14/tip3p.xml parameters" in text


def test_dhla_rigid_and_no_warning(dhla_tip3p):
    systems, topology, _, text = dhla_tip3p
    waters = _waters(topology)
    for system in _all_systems(systems):
        assert system.getNumConstraints() == 3 * len(waters)
    assert WARNING not in text


# ---------------------------------------------------------------------- run


@pytest.fixture(scope="module")
def tip3p_run(tmp_path_factory):
    """Solvated ethanol in rigid TIP3P, 3 lambda, 1 replica. About 85 s on
    CPU, nearly all of it in the 7 constrained minimizations (cSPC/E: 18 s);
    the window minimizations are needed (without them lambda 0.5 NaNs)."""
    folder = tmp_path_factory.mktemp("tip3p_run")
    pair = EvbTestHelper.evb_ff_pair()
    evb = _driver("tip3p")
    EvbTestHelper.evb_seed_states(evb, [pair.reactant, pair.product])
    conf = evb.default_system_configurations("water_NVT")
    conf.update({
        "name": "water", "temperature": 300.0, "padding": 1.0,
        "sample_steps": 20, "write_step": 5, "equil_NVT_steps": 20,
        "equil_NPT_steps": 0, "initial_equil_NVT_steps": 50,
        "initial_equil_NPT_steps": 0, "minimize_every_lambda": True,
        "n_replicas": 1, "trajectory_format": "dcd",
    })
    previous = os.getcwd()
    os.chdir(folder)
    try:
        evb.build_systems(configurations=[conf], Lambda=LAMBDAS)
        evb.run_FEP(platform=EvbTestHelper.evb_platform())
    finally:
        os.chdir(previous)
    return evb, folder / evb.system_confs[0]["data_folder"]


@pytest.mark.timeconsuming
def test_run_keeps_water_rigid(tip3p_run):
    """Every frame keeps the rigid water geometry, and the energies are
    finite."""
    evb, data_folder = tip3p_run
    topology = mmapp.PDBxFile(str(data_folder / "topology.cif")).topology
    waters = np.array(_waters(topology))
    r_oh, r_hh = ReactionSystemBuilder._rigid_water_distances(
        get_water_parameters()["tip3p"])
    with EvbFepDriver._open_trajectory(data_folder / "trajectory.dcd",
                                       True) as (n_frames, read):
        assert n_frames == len(LAMBDAS) * 2 * 4
        frames = list(read(0, n_frames))
    for pos, _ in frames:
        d_oh = np.linalg.norm(pos[waters[:, 1:]] - pos[waters[:, :1]], axis=2)
        d_hh = np.linalg.norm(pos[waters[:, 1]] - pos[waters[:, 2]], axis=1)
        # float32 Angstrom in the dcd
        assert np.abs(d_oh - r_oh).max() < 2e-5
        assert np.abs(d_hh - r_hh).max() < 2e-5
    energies = np.genfromtxt(data_folder / "Energies.csv", delimiter=",",
                             skip_header=1)
    assert np.all(np.isfinite(energies))


@pytest.mark.timeconsuming
def test_run_records_the_model(tip3p_run):
    evb, _ = tip3p_run
    conf = evb.system_confs[0]
    assert (conf["solvent"], conf["water_model"]) == ("tip3p", "tip3p")


def test_remove_added_constraints():
    system = mm.System()
    for _ in range(4):
        system.addParticle(1.0)
    system.addConstraint(0, 1, 0.1)
    n_keep = system.getNumConstraints()
    system.addConstraint(2, 3, 0.1)
    system.addConstraint(1, 2, 0.1)
    EvbFepDriver._remove_added_constraints(system, n_keep)
    assert system.getNumConstraints() == 1
    assert system.getConstraintParameters(0)[:2] == [0, 1]


def test_constrain_H_bonds_keeps_rigid_water_once():
    """A water already constrained (rigid) is not constrained again, and the
    removal after the equilibration leaves exactly the rigid constraints (it
    used to remove every constraint of the system)."""
    system = mm.System()
    for mass in (15.999, 1.008, 1.008, 12.011, 1.008):
        system.addParticle(mass)
    for pair in ((0, 1), (0, 2), (1, 2)):
        system.addConstraint(*pair, 0.1)
    bonds = mm.HarmonicBondForce()
    bonds.addBond(0, 1, 0.1, 1000.0)  # the water O-H, constrained already
    bonds.addBond(3, 4, 0.109, 1000.0)  # a solute C-H
    system.addForce(bonds)
    fep = EvbFepDriver(ostream=OutputStream(None))
    n_keep = system.getNumConstraints()
    fep._constrain_H_bonds(system)
    pairs = [tuple(system.getConstraintParameters(i)[:2])
             for i in range(system.getNumConstraints())]
    assert pairs == [(0, 1), (0, 2), (1, 2), (3, 4)]
    fep._remove_added_constraints(system, n_keep)
    assert system.getNumConstraints() == 3


def test_rigid_hh_geometry():
    r_oh, r_hh = ReactionSystemBuilder._rigid_water_distances(
        get_water_parameters()["tip3p"])
    assert r_hh == pytest.approx(2 * 0.09572 * math.sin(math.radians(52.26)))
