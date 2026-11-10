"""Residue cutting in the PDB route of the EVB system builder.

DhlA (PDB 2DHD, reactant complex) with Asp124 cut between CB and CG: only the
carboxylate enters the reacting fragment, formate + 1,2-dichloroethane ->
2-chloroethyl formate + Cl-, and the rest of Asp124 stays in amber14. The
formate H is the cap that stands for CB.

Checked at every lambda and in both PES states:
- every fragment atom carries its fragment charge, and CB carries the cap's;
- the rest of Asp124 sums to an integer, and so does the whole system;
- a 1-4 pair between the fragment and the protein follows the current charges.
The two carboxylate oxygens (and the two ends of DCE) are equivalent to the
graph match, so the reacting atoms must be named with 'reacting_atoms'.
"""

import numpy as np
import pytest

from veloxchem.molecule import Molecule
from veloxchem.mmforcefieldgenerator import MMForceFieldGenerator
from veloxchem.reactionsystembuilder import ReactionSystemBuilder
from veloxchem.outputstream import OutputStream

from test_evb_helper import EvbTestHelper

mm = pytest.importorskip("openmm")
unit = pytest.importorskip("openmm.unit")

LAMBDAS = [0.0, 0.5, 1.0]
EVB_ATOMS = ("CG", "OD1", "OD2")


def load_states():
    data = EvbTestHelper.evb_data_dir()
    states = []
    for name in ("reactant", "product"):
        ff = MMForceFieldGenerator.load_forcefield_from_json_file(
            str(data / f"evb_dhla_{name}_ff_data.json"),
            ostream=OutputStream(None))
        ff.molecule = Molecule.read_xyz_file(str(data / f"evb_dhla_{name}.xyz"))
        states.append(ff)
    return states


def config(reacting=True):
    asp = {"chain": "A", "residue": 124, "cut": "beta"}
    dce = {"chain": "B", "residue": 1}
    if reacting:
        asp["reacting_atoms"] = ["OD2"]
        dce["reacting_atoms"] = ["C1", "CL1"]
    return {
        "name": "dhla_vac",
        "temperature": 300.0,
        "pdb": str(EvbTestHelper.evb_data_dir() / "evb_dhla_2dhd_reactant.cif"),
        "pdb_active_res": [asp, dce],
    }


def build(reacting=True):
    reactant, product = load_states()
    builder = ReactionSystemBuilder(ostream=OutputStream(None))
    with EvbTestHelper.evb_chdir_tmp():
        systems, topology, _ = builder.build_systems(reactant, product,
                                                     LAMBDAS, config(reacting))
    return builder, systems, topology, reactant, product


@pytest.fixture(scope="module")
def dhla():
    return build()


def charges(system):
    nb = next(f for f in system.getForces()
              if isinstance(f, mm.NonbondedForce))
    return np.array([
        nb.getParticleParameters(i)[0].value_in_unit(unit.elementary_charge)
        for i in range(nb.getNumParticles())
    ])


def systems_and_lambdas(systems):
    yield from ((systems[lam], lam) for lam in LAMBDAS)
    yield systems["reactant"], 0.0
    yield systems["product"], 1.0


def test_reacting_atoms_must_be_named():
    with pytest.raises(Exception, match="reacting_atoms"):
        build(reacting=False)


def test_reacting_atoms_carry_the_reaction(dhla):
    builder, _, _, reactant, product = dhla
    forming = set(product.bonds) - set(reactant.bonds)
    breaking = set(reactant.bonds) - set(product.bonds)
    names = {i: atom.name for i, atom in builder.reaction_atoms.items()}
    assert [sorted(names[i] for i in pair) for pair in forming] == [["C1", "OD2"]]
    assert [sorted(names[i] for i in pair) for pair in breaking] == [["C1", "CL1"]]
    (cap, host), = builder._cut_caps.items()
    assert (host.residue.id, host.name) == ("124", "CB")


def test_charges(dhla):
    builder, systems, topology, reactant, product = dhla
    asp = next(r for r in topology.residues() if r.id == "124")
    rest = [a.index for a in asp.atoms() if a.name not in EVB_ATOMS + ("CB", )]
    (cap, host), = builder._cut_caps.items()
    for system, lam in systems_and_lambdas(systems):
        q = charges(system)
        for i, atom in reactant.atoms.items():
            expected = ((1 - lam) * atom["charge"] +
                        lam * product.atoms[i]["charge"])
            index = host.index if i == cap else builder.reaction_atoms[i].index
            assert q[index] == pytest.approx(expected, abs=1e-10)
        assert q[rest].sum() == pytest.approx(round(q[rest].sum()), abs=1e-9)
        assert q.sum() == pytest.approx(round(q.sum()), abs=1e-9)


def test_fragment_protein_14_pair(dhla):
    """OD2-CG-CB-CA: a 1-4 pair across the cut, scaled by amber's 1/1.2."""
    builder, systems, topology, _, _ = dhla
    asp = {a.name: a.index for a in topology.atoms() if a.residue.id == "124"
           and a.residue.chain.id == "A"}
    pair = {asp["OD2"], asp["CA"]}
    for system, _ in systems_and_lambdas(systems):
        nb = next(f for f in system.getForces()
                  if isinstance(f, mm.NonbondedForce))
        q = charges(system)
        qq = next(
            nb.getExceptionParameters(k)[2]
            for k in range(nb.getNumExceptions())
            if set(nb.getExceptionParameters(k)[:2]) == pair)
        assert qq.value_in_unit(unit.elementary_charge**2) == pytest.approx(
            q[asp["OD2"]] * q[asp["CA"]] / 1.2, abs=1e-10)
