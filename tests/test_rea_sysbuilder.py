"""Smoke tests for the reaction system builder.

Two kinds of checks:
  * serialized-XML reference compare for the canonical vacuum / water / implicit
    configurations, and
  * a construct-and-evaluate matrix over the many configuration branches: build
    the systems and assert every lambda window yields a finite potential energy
    (the cheap correctness oracle). No MD is run.
"""

from types import SimpleNamespace

import pytest

from veloxchem.molecule import Molecule
from veloxchem.mmforcefieldgenerator import MMForceFieldGenerator
from veloxchem.reactionsystembuilder import ReactionSystemBuilder

from veloxchem.outputstream import OutputStream

from test_evb_helper import EvbTestHelper

pytestmark = [pytest.mark.timeconsuming]

pytest.importorskip('openmm')


class SystemBuilderFixtures:
    """Shared setup for the system-builder tests."""

    @staticmethod
    def _build(ff_pair, config, lambdas=(0, 0.4, 1), constraints=None):
        builder = ReactionSystemBuilder(ostream=OutputStream(None))
        builder.water_model = "cspce"
        systems, topology, positions = builder.build_systems(
            ff_pair.reactant,
            ff_pair.product,
            list(lambdas),
            config,
            constraints=constraints,
        )
        return systems, topology, positions

    @staticmethod
    def _branch_configs():
        return {
            "vacuum": ({
                "name": "vacuum",
                "temperature": 300.0
            }, None),
            "vacuum_soft_core": ({
                "name": "vacuum",
                "temperature": 300.0,
                "soft_core_coulomb_pes_static": True,
                "soft_core_lj_pes_static": True,
            }, None),
            "water_NVT": ({
                "name": "water_cspce_NVT",
                "solvent": "cspce",
                "temperature": 300.0,
                "pressure": 1,
                "padding": 1.0,
                "ion_count": 0,
                "neutralize": False,
            }, None),
            "water_NPT": ({
                "name": "water_cspce_NPT",
                "solvent": "cspce",
                "temperature": 300.0,
                "pressure": 1,
                "isobaric": True,
                "padding": 1.0,
                "ion_count": 0,
                "neutralize": False,
            }, None),
            "implicit_gbn": ({
                "name": "vacuum",
                "temperature": 300.0,
                "implicit_solvent_model": "gbn",
                "solute_dielectric": 1.0,
                "solvent_dielectric": 78.39,
            }, None),
            "E_field": ({
                "name": "efield_vacuum",
                "temperature": 300.0,
                "E_field": [0, 0, 10],
            }, None),
            "frozen_atoms": ({
                "name": "vacuum",
                "temperature": 300.0,
                "frozen_atoms": [0, 1],
            }, None),
            "CNT": ({
                "name": "cnt",
                "temperature": 300.0,
                "solvent": "cspce",
                "padding": 1.0,
                "CNT": True,
                "CNT_radius_nm": 0.7,
            }, None),
            "graphene": ({
                "name": "graphene",
                "temperature": 300.0,
                "solvent": "cspce",
                "padding": 1.0,
                "graphene": True,
                "graphene_size_nm": 2.5,
            }, None),
            "constraints": ({
                "name": "vacuum",
                "temperature": 300.0,
            }, [{
                (0, 1): {
                    "type": "harmonic",
                    "equilibrium": 1.5,
                    "force_constant": 1000.0,
                }
            }]),
        }

class TestSystemBuilderReference:
    """Serialized-XML compare against committed references (brittle but exact)."""

    def test_vacuum_water_implicit(self):
        data_dir = EvbTestHelper.evb_data_dir()
        vac_conf = {"name": "vacuum", "temperature": 300.0}

        wat_conf = {
            "name": "water_cspce_NVT",
            "solvent": "cspce",
            "temperature": 300.0,
            "pressure": 1,
            "padding": 1.5,
            "ion_count": 0,
            "neutralize": False,
        }

        impl_conf = {
            "name": "vacuum",
            "temperature": 300.0,
            "implicit_solvent_model": "gbn",
            "solute_dielectric": 1.0,
            "solvent_dielectric": 78.39,
        }

        # A fresh builder / fresh force field is needed for each run, because the
        # builder annotates the force field it is handed.
        lambdas = [0, 0.4, 1]
        vac_systems, _, _ = SystemBuilderFixtures._build(EvbTestHelper.evb_ff_pair(), vac_conf, lambdas)
        wat_systems, _, _ = SystemBuilderFixtures._build(EvbTestHelper.evb_ff_pair(), wat_conf, lambdas)
        impl_systems, _, _ = SystemBuilderFixtures._build(EvbTestHelper.evb_ff_pair(), impl_conf, lambdas)

        for lam, tag in [(0, "0.000"), (0.4, "0.400"), (1, "1.000")]:
            EvbTestHelper.evb_compare_serialized_systems(
                vac_systems[lam], data_dir / f"evb_ethanol_vac_{tag}_sys.xml")
            EvbTestHelper.evb_compare_serialized_systems(
                wat_systems[lam], data_dir / f"evb_ethanol_solv_{tag}_sys.xml")
            EvbTestHelper.evb_compare_serialized_systems(
                impl_systems[lam], data_dir / f"evb_ethanol_impl_{tag}_sys.xml")


# Configuration branches exercised by the energy-finite oracle. Each entry is
# (config, constraints). Kept small/cheap; solvent shells use tight padding.
class TestSystemBuilderBranches:
    """Build each configuration branch and assert finite energies at 0/0.5/1."""

    @pytest.mark.parametrize("name", list(SystemBuilderFixtures._branch_configs().keys()))
    def test_branch_energy_finite(self, name):
        config, constraints = SystemBuilderFixtures._branch_configs()[name]
        lambdas = [0, 0.5, 1]
        systems, topology, positions = SystemBuilderFixtures._build(EvbTestHelper.evb_ff_pair(), config, lambdas,
                                              constraints)
        for lam in lambdas:
            EvbTestHelper.evb_assert_energy_finite(systems[lam],
                                     positions,
                                     label=f"{name}@lambda={lam}")


class TestPdbInput:
    """Fully-verified PDB-input path: the PETase A60S enzyme with a
    bound reacting ligand (chain B res 1) and the Ser/His/Asp catalytic triad
    (chain A res 208/177/131), built in vacuum. Exercises amber14 template
    generation for the reacting residues, reactant<->residue atom mapping, force
    deletion on reacting atoms, and reaction-force interpolation. Compared against
    committed gzipped reference systems."""

    def test_petase_pdb_vacuum(self):
        data_dir = EvbTestHelper.evb_data_dir()
        reactant = MMForceFieldGenerator.load_forcefield_from_json_file(
            str(data_dir / "evb_petase_reactant_ff_data.json"),
            ostream=OutputStream(None))
        reactant.molecule = Molecule.read_xyz_file(
            str(data_dir / "evb_petase_reactant.xyz"))
        product = MMForceFieldGenerator.load_forcefield_from_json_file(
            str(data_dir / "evb_petase_product_ff_data.json"),
            ostream=OutputStream(None))
        product.molecule = Molecule.read_xyz_file(
            str(data_dir / "evb_petase_product.xyz"))

        config = {
            "name":
            "petase_A60S_vac",
            "temperature":
            300.0,
            "pdb":
            str(data_dir / "evb_petase_A60S.cif"),
            "pdb_active_res": [
                {
                    'chain': 'B',
                    'residue': 1
                },
                {
                    'chain': 'A',
                    'residue': 208
                },
                {
                    'chain': 'A',
                    'residue': 177
                },
                {
                    'chain': 'A',
                    'residue': 131
                },
            ],
        }

        ff_pair = SimpleNamespace(reactant=reactant, product=product)
        lambdas = [0, 0.4, 1]
        # The builder writes intermediate files into the cwd.
        with EvbTestHelper.evb_chdir_tmp():
            systems, topology, positions = SystemBuilderFixtures._build(ff_pair, config, lambdas)

        for lam, tag in [(0, "0.000"), (0.4, "0.400"), (1, "1.000")]:
            EvbTestHelper.evb_compare_serialized_systems(
                systems[lam], data_dir / f"evb_petase_pdb_vac_{tag}_sys.xml.gz")


class TestSystemBuilderValidation:
    """Configuration type-validation should raise rather than silently pass."""

    def test_wrong_type_raises(self):
        # padding is declared float; a string must be rejected.
        bad_conf = {
            "name": "vacuum",
            "temperature": 300.0,
            "padding": "not-a-number",
        }
        with pytest.raises((TypeError, ValueError)):
            SystemBuilderFixtures._build(EvbTestHelper.evb_ff_pair(), bad_conf)
