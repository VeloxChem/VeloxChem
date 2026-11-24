"""Energies of chosen force-field terms and residues recomputed from a frame
(EvbTrajectoryAnalysis), on the DhlA reactant complex in vacuum.

- A tracked bond reports its length and the energy of its harmonic term.
- Asp124, cut at CB-CG, is tracked as its reacting part (CG, OD1, OD2) and
  the rest; the bond across the cut is in the [cut] column.
- Without a cutoff, the residue-rest interaction taken from the system's
  NonbondedForce (all - residue off - residue alone) equals its pairwise
  split over the atom classes.
- The class-pair columns add up to the NonbondedForce energy.
"""

import copy

import numpy as np
import pytest

from veloxchem.molecule import Molecule
from veloxchem.mmforcefieldgenerator import MMForceFieldGenerator
from veloxchem.reactionsystembuilder import ReactionSystemBuilder
from veloxchem.evbtrajectoryanalysis import EvbTrajectoryAnalysis
from veloxchem.outputstream import OutputStream

from test_evb_helper import EvbTestHelper

mm = pytest.importorskip("openmm")
mmunit = pytest.importorskip("openmm.unit")


def build_dhla():
    data = EvbTestHelper.evb_data_dir()
    states = []
    for name in ("reactant", "product"):
        ff = MMForceFieldGenerator.load_forcefield_from_json_file(
            str(data / f"evb_dhla_{name}_ff_data.json"),
            ostream=OutputStream(None))
        ff.molecule = Molecule.read_xyz_file(str(data / f"evb_dhla_{name}.xyz"))
        states.append(ff)
    config = {
        "name": "dhla_vac",
        "temperature": 300.0,
        "pdb": str(data / "evb_dhla_2dhd_reactant.cif"),
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
    builder = ReactionSystemBuilder(ostream=OutputStream(None))
    with EvbTestHelper.evb_chdir_tmp():
        systems, topology, _ = builder.build_systems(*states, [0.0, 1.0],
                                                     config)
    return systems, topology, np.asarray(builder.positions) * 0.1


def test_trajectory_analysis_dhla_vacuum():
    systems, topology, positions = build_dhla()
    x = positions + np.random.default_rng(0).normal(0.0, 0.003,
                                                    positions.shape)

    analysis = EvbTrajectoryAnalysis(topology,
                                     systems["reactant"],
                                     terms={"cacb": ["A:124:CA", "A:124:CB"]},
                                     residues=["A:124"],
                                     hamiltonians=("E1", ),
                                     platform="Reference",
                                     ostream=OutputStream(None))
    frame = analysis.frame_geometry(x, None)
    values = dict(frame["values"])
    evaluator = analysis.evaluator(systems["reactant"], "E1")
    evaluator.evaluate(frame, x, None, values)
    evaluator.close()

    bonds = next(f for f in systems["reactant"].getForces()
                 if f.getName() == "PDB Bond Force")

    def harmonic(i, j):
        r = np.linalg.norm(x[i] - x[j])
        r0, k = next((p[2]._value, p[3]._value)
                     for p in map(bonds.getBondParameters,
                                  range(bonds.getNumBonds()))
                     if {p[0], p[1]} == {i, j})
        return r, 0.5 * k * (r - r0)**2

    ca, cb = analysis.terms[0]["atoms"]
    r, energy = harmonic(ca, cb)
    assert values["cacb/distance"] == pytest.approx(r, abs=1e-6)
    assert values["cacb/E1/bonded"] == pytest.approx(energy, rel=1e-6)

    part, rest = analysis.residues
    assert part["label"] == "A:124:ASP[fragment]"
    assert sorted(analysis._atoms[i].name for i in part["atoms"]) == [
        "CG", "OD1", "OD2"
    ]
    assert cb in rest["atoms"]
    cg = next(i for i in part["atoms"] if analysis._atoms[i].name == "CG")
    assert values["A:124:ASP[cut]/E1/bonded/bond"] == pytest.approx(
        harmonic(cb, cg)[1], rel=1e-6)

    for residue in analysis.residues:
        base = f"{residue['label']}/E1/nb_rest"
        for kind in ("coul", "lj"):
            pairwise = sum(values[f"{base}/{c}/{kind}"]
                           for c in analysis.class_names)
            assert abs(values[f"{base}/{kind}"]) > 1.0
            assert values[f"{base}/{kind}"] == pytest.approx(pairwise,
                                                             abs=1e-6)

    nonbonded = mm.System()
    for i in range(systems["reactant"].getNumParticles()):
        nonbonded.addParticle(1.0)
    nonbonded.addForce(
        copy.deepcopy(next(f for f in systems["reactant"].getForces()
                           if isinstance(f, mm.NonbondedForce))))
    context = mm.Context(nonbonded, mm.VerletIntegrator(0.001),
                         mm.Platform.getPlatformByName("Reference"))
    context.setPositions(x)
    reference = context.getState(getEnergy=True).getPotentialEnergy()
    pairs = sum(values[f"classes/E1/{a}~{b}/{kind}"]
                for a, b in analysis.class_pair_list()
                for kind in ("coul", "lj"))
    assert pairs == pytest.approx(reference.value_in_unit(
        mmunit.kilojoule_per_mole), abs=1e-6)
