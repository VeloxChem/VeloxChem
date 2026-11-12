"""E3-1 dev tests: `neutralize` on the PDB route of ReactionSystemBuilder.

Dev-only (test_claude_*): never to master. The solute molecule of the PDB
route comes from Molecule.read_pdb_file, whose charge is 0, so SolvationBuilder
used to add no counter-ions and left the box charged (DhlA: -17 e). The total
charge is now counted from the PDB system (protein and the rest of any cut
residue) plus the fragment's net charge, before the solvent is added.

Two builds, solvated in cspce with a small padding:
- DhlA with Asp124 cut at 'beta' (protein -16 e, fragment -1 e): 17 Na+;
- PETase with whole residues (protein +7 e, fragment -1 e): 6 Cl-.

SolvationBuilder._get_volume loops in Python over every solute atom pair
(DhlA: 11.8 M pairs, ~35 s). The builds use `_fast_volume`, the same formula
over the overlapping pairs only; the volume sets the number of waters, not
the charge. `test_fast_volume_is_exact` checks it against the original.
"""

import numpy as np
import pytest
from scipy.spatial import cKDTree

from veloxchem.molecule import Molecule
from veloxchem.mmforcefieldgenerator import MMForceFieldGenerator
from veloxchem.reactionsystembuilder import ReactionSystemBuilder
from veloxchem.outputstream import OutputStream
from veloxchem.solvationbuilder import SolvationBuilder
from veloxchem.veloxchemlib import bohr_in_angstrom

from test_evb_helper import EvbTestHelper

mm = pytest.importorskip("openmm")
import openmm.unit as mmunit  # noqa: E402

DATA = EvbTestHelper.evb_data_dir()
LAMBDAS = [0.0, 0.5, 1.0]
PADDING = 0.2

CASES = {
    "dhla_cut": ("evb_dhla", "evb_dhla_2dhd_reactant.cif", [
        {"chain": "A", "residue": 124, "cut": "beta",
         "reacting_atoms": ["OD2"]},
        {"chain": "B", "residue": 1, "reacting_atoms": ["C1", "CL1"]},
    ], -17),
    "petase_whole": ("evb_petase", "evb_petase_A60S.cif", [
        {"chain": "B", "residue": 1},
        {"chain": "A", "residue": 208},
        {"chain": "A", "residue": 177, "reacting_atoms": ["OD2"]},
        {"chain": "A", "residue": 131},
    ], 6),
}


def _load(prefix, name):
    ff = MMForceFieldGenerator.load_forcefield_from_json_file(
        str(DATA / f"{prefix}_{name}_ff_data.json"),
        ostream=OutputStream(None))
    ff.molecule = Molecule.read_xyz_file(str(DATA / f"{prefix}_{name}.xyz"))
    return ff


def _fast_volume(self, molecule):
    """SolvationBuilder._get_volume, summed over overlapping pairs only."""
    radii = molecule.vdw_radii_to_numpy() * bohr_in_angstrom()
    coords = molecule.get_coordinates_in_angstrom()
    volume = np.sum(4.0 / 3.0 * np.pi * radii**3)
    pairs = cKDTree(coords).query_pairs(2.0 * radii.max(),
                                        output_type="ndarray")
    i, j = pairs[:, 0], pairs[:, 1]
    d = np.linalg.norm(coords[i] - coords[j], axis=1)
    r1, r2 = radii[i], radii[j]
    keep = (d < r1 + r2) & (d > 0)
    d, r1, r2 = d[keep], r1[keep], r2[keep]
    h1 = (r1 - r2 + d) * (r1 + r2 - d) / (2 * d)
    h2 = (r2 - r1 + d) * (r2 + r1 - d) / (2 * d)
    volume -= np.sum(np.pi * h1**2 * (3 * r1 - h1) +
                     np.pi * h2**2 * (3 * r2 - h2)) / 6
    return round(float(volume), 2)


@pytest.fixture(scope="module", params=list(CASES))
def solvated(request, tmp_path_factory):
    prefix, cif, residues, pdb_charge = CASES[request.param]
    out = tmp_path_factory.mktemp(request.param) / "build.out"
    config = {
        "name": request.param,
        "temperature": 300.0,
        "pdb": str(DATA / cif),
        "pdb_active_res": residues,
        "solvent": "cspce",
        "neutralize": True,
        "padding": PADDING,
    }
    builder = ReactionSystemBuilder(ostream=OutputStream(str(out)))
    builder.water_model = "cspce"
    with pytest.MonkeyPatch.context() as patch:
        patch.setattr(SolvationBuilder, "_get_volume", _fast_volume)
        with EvbTestHelper.evb_chdir_tmp():
            systems, topology, _ = builder.build_systems(
                _load(prefix, "reactant"), _load(prefix, "product"),
                list(LAMBDAS), config)
    return systems, topology, pdb_charge, out.read_text()


def _charges(system):
    nb = next(f for f in system.getForces()
              if isinstance(f, mm.NonbondedForce))
    return np.array([
        nb.getParticleParameters(i)[0].value_in_unit(mmunit.elementary_charge)
        for i in range(nb.getNumParticles())
    ])


def _all_systems(systems):
    keys = list(LAMBDAS) + ["reactant", "product"]
    return [systems[key] for key in keys]


def test_counterions(solvated):
    systems, topology, pdb_charge, _ = solvated
    ions = [a.index for a in topology.atoms() if a.residue.name == "ION"]
    assert len(ions) == abs(pdb_charge)
    for system in _all_systems(systems):
        assert np.allclose(_charges(system)[ions], -np.sign(pdb_charge),
                           atol=1e-10)


def test_systems_are_neutral(solvated):
    systems, _, _, _ = solvated
    for system in _all_systems(systems):
        assert _charges(system).sum() == pytest.approx(0.0, abs=1e-6)


def test_output_names_the_charge(solvated):
    _, _, pdb_charge, text = solvated
    assert f"Total charge of the PDB system: {pdb_charge:d}" in text
    assert f"The solute has a charge of {pdb_charge}" in text
    if pdb_charge < 0:
        assert f"adding {abs(pdb_charge)} Na counterions" in text


def test_non_integer_charge_refused(monkeypatch):
    """A fragment whose charges do not add up with the protein's to an
    integer stops the build instead of being rounded."""
    prefix, cif, residues, _ = CASES["dhla_cut"]
    reactant, product = _load(prefix, "reactant"), _load(prefix, "product")
    for state in (reactant, product):
        state.atoms[0]["charge"] += 0.01
    builder = ReactionSystemBuilder(ostream=OutputStream(None))
    builder.water_model = "cspce"

    def no_solvent(*args, **kwargs):
        raise AssertionError("solvation reached")

    monkeypatch.setattr(builder, "_add_solvent", no_solvent)
    config = {"name": "dhla", "temperature": 300.0, "pdb": str(DATA / cif),
              "pdb_active_res": residues, "solvent": "cspce",
              "neutralize": True, "padding": PADDING}
    with EvbTestHelper.evb_chdir_tmp():
        with pytest.raises(Exception, match="not an integer"):
            builder.build_systems(reactant, product, list(LAMBDAS), config)


def test_fast_volume_is_exact():
    """The test's volume equals SolvationBuilder's on a small molecule with
    overlapping and non-overlapping pairs."""
    molecule = Molecule.read_xyz_file(str(DATA / "evb_petase_reactant.xyz"))
    builder = SolvationBuilder(ostream=OutputStream(None))
    assert _fast_volume(builder, molecule) == pytest.approx(
        builder._get_volume(molecule), abs=0.011)
