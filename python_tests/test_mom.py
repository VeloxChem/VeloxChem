import numpy as np
import pytest

from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.scfrestopendriver import ScfRestrictedOpenDriver
from veloxchem.scfunrestdriver import ScfUnrestrictedDriver


class TestPDFT:

    def test_MOM_ROHF(self):
        water_xyz = """
        O       0.0000000000     0.0000000000     0.1178336003
        H      -0.7595754146    -0.0000000000    -0.4713344012
        H       0.7595754146     0.0000000000    -0.4713344012
        """

        molecule = Molecule.read_str(water_xyz)
        basis = MolecularBasis.read(molecule, "6-31G")

        # Ground state
        scf_gs = ScfRestrictedDriver()
        scf_gs.compute(molecule, basis)

        # Core-hole
        molecule.set_charge(1)
        molecule.set_multiplicity(2)
        scf_ch = ScfRestrictedOpenDriver()

        scf_ch.maximum_overlap(molecule, scf_gs.mol_orbs, [0,1,2,3,4], [1,2,3,4])
        scf_ch.compute(molecule,basis)

        assert abs(scf_ch.get_scf_energy() + 56.070929384622) < 1.0e-6

    def test_MOM_UDFT(self):
        water_xyz = """
        O       0.0000000000     0.0000000000     0.1178336003
        H      -0.7595754146    -0.0000000000    -0.4713344012
        H       0.7595754146     0.0000000000    -0.4713344012
        """

        molecule = Molecule.read_str(water_xyz)
        basis = MolecularBasis.read(molecule, "6-31G")

        # Ground state
        scf_gs = ScfRestrictedDriver()
        scf_gs.xcfun = "B3LYP"
        scf_gs.compute(molecule, basis)

        # Core-hole
        molecule.set_charge(1)
        molecule.set_multiplicity(2)
        scf_ch = ScfUnrestrictedDriver()
        scf_ch.xcfun = "B3LYP"

        scf_ch.maximum_overlap(molecule, scf_gs.mol_orbs, [0,1,2,3,4], [1,2,3,4])
        scf_ch.compute(molecule,basis)

        assert abs(scf_ch.get_scf_energy() + 56.4455719055) < 1.0e-6
