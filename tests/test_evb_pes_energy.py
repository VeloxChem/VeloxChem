"""The EVB PES systems against the force field they are built from.

For Cl- + CH3Cl -> ClCH3 + Cl-, each PES system must equal its state's force
field: harmonic bonds and angles, 1-2 and 1-3 pairs excluded by that state's
own bonds, plain Coulomb and Lennard-Jones otherwise. By design, the reacting
C-Cl bond is a Morse bond, and the nonbonded terms of the unbonded reacting
C...Cl pair are soft-core (second-order Taylor polynomial below r_c).
"""

import json

import numpy as np
import pytest

from veloxchem.evbdriver import EvbDriver
from veloxchem.molecule import Molecule
from veloxchem.mmforcefieldgenerator import MMForceFieldGenerator
from veloxchem.outputstream import OutputStream

from test_evb_helper import EvbTestHelper

mm = pytest.importorskip("openmm")

COULOMB_KE = 138.935458  # kJ mol^-1 nm e^-2
MORSE_D = 500.0  # kJ/mol, the builder's default dissociation energy
STATES = ("reactant", "product")


def key(text):
    return tuple(int(i) for i in text.strip("()").split(","))


def coulomb(r, qq):
    return COULOMB_KE * qq / r


def lennard_jones(r, sigma, epsilon):
    return 4 * epsilon * ((sigma / r)**12 - (sigma / r)**6)


def soft_core(f, r, rc, *params):
    """f for r >= rc, else its second-order Taylor polynomial about rc."""
    if r >= rc:
        return f(r, *params)
    h = 1e-5 * rc
    f0, fp, fm = f(rc, *params), f(rc + h, *params), f(rc - h, *params)
    return (f0 + (fp - fm) / (2 * h) * (r - rc) + 0.5 *
            (fp - 2 * f0 + fm) / h**2 * (r - rc)**2)


def reference_energy(ff, reactive, xyz):
    """PES energy (kJ/mol) of one state at xyz (nm), from its force field."""

    atoms = ff["atoms"]
    bonds = {key(k): b for k, b in ff["bonds"].items()}
    energy = 0.0
    for (i, j), b in bonds.items():
        r = np.linalg.norm(xyz[i] - xyz[j])
        k, r0 = b["force_constant"], b["equilibrium"]
        if (i, j) in reactive:
            a = np.sqrt(k / (2 * MORSE_D))
            energy += MORSE_D * (1 - np.exp(-a * (r - r0)))**2
        else:
            energy += 0.5 * k * (r - r0)**2
    for text, angle in ff["angles"].items():
        i, j, k = key(text)
        u, v = xyz[i] - xyz[j], xyz[k] - xyz[j]
        theta = np.arccos(u @ v / (np.linalg.norm(u) * np.linalg.norm(v)))
        energy += 0.5 * angle["force_constant"] * (
            theta - np.radians(angle["equilibrium"]))**2

    neighbours = [set() for _ in atoms]
    for i, j in bonds:
        neighbours[i].add(j)
        neighbours[j].add(i)
    for i in range(len(atoms)):
        for j in range(i + 1, len(atoms)):
            # 1-2 and 1-3 pairs are excluded; this system has no 1-4 pairs
            if j in neighbours[i] or neighbours[i] & neighbours[j]:
                continue
            ai, aj = atoms[str(i)], atoms[str(j)]
            qq = ai["charge"] * aj["charge"]
            sigma = 0.5 * (ai["sigma"] + aj["sigma"])
            epsilon = np.sqrt(ai["epsilon"] * aj["epsilon"])
            r = np.linalg.norm(xyz[i] - xyz[j])
            if (i, j) in reactive:
                rc_q = 0.3 * (1 + qq)
                rc_lj = 0.7 * ((26 / 7) * 4 * epsilon * sigma**6)**(1 / 6)
                energy += soft_core(coulomb, r, rc_q, qq)
                energy += soft_core(lennard_jones, r, rc_lj, sigma, epsilon)
            else:
                energy += coulomb(r, qq) + lennard_jones(r, sigma, epsilon)
    return energy


def test_pes_energy_equals_force_field():
    data = EvbTestHelper.evb_data_dir()
    ffs, states, geometries = [], [], []
    for name in STATES:
        path = data / f"evb_clch3cl_{name}_ff_data.json"
        ffs.append(json.loads(path.read_text()))
        state = MMForceFieldGenerator.load_forcefield_from_json_file(
            str(path), ostream=OutputStream(None))
        state.molecule = Molecule.read_xyz_file(
            str(data / f"evb_clch3cl_{name}.xyz"))
        states.append(state)
        geometries.append(state.molecule.get_coordinates_in_angstrom() * 0.1)

    evb = EvbDriver(ostream=OutputStream(None))
    EvbTestHelper.evb_seed_states(evb, states)

    with EvbTestHelper.evb_chdir_tmp():
        evb.build_systems(configurations=["vacuum"], Lambda=[0.0, 1.0])
        run = evb.system_confs[0]["data_folder"] + "/run"
        systems = [
            mm.XmlSerializer.deserialize(open(f"{run}/{name}_sys.xml").read())
            for name in STATES
        ]

    platform = mm.Platform.getPlatformByName("Reference")
    for ff, system in zip(ffs, systems):
        context = mm.Context(system, mm.VerletIntegrator(1e-3), platform)
        for xyz in geometries:
            context.setPositions(xyz * mm.unit.nanometer)
            energy = context.getState(getEnergy=True).getPotentialEnergy()
            assert energy.value_in_unit(
                mm.unit.kilojoule_per_mole) == pytest.approx(
                    reference_energy(ff, evb.reactive_bonds, xyz), rel=1e-6)
