from pathlib import Path

from veloxchem.veloxchemlib import is_mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.scfrestopendriver import ScfRestrictedOpenDriver
from veloxchem.scfunrestdriver import ScfUnrestrictedDriver


class TestMOM:

    def test_MOM_ROHF(self):

        water_xyz = """
        O       0.0000000000     0.0000000000     0.1178336003
        H      -0.7595754146    -0.0000000000    -0.4713344012
        H       0.7595754146     0.0000000000    -0.4713344012
        """

        molecule = Molecule.read_str(water_xyz)
        basis = MolecularBasis.read(molecule, '6-31G', ostream=None)

        # Ground state
        scf_gs = ScfRestrictedDriver()
        scf_gs.ostream.mute()
        scf_gs.compute(molecule, basis)

        # Core-hole
        molecule.set_charge(1)
        molecule.set_multiplicity(2)
        scf_ch = ScfRestrictedOpenDriver()
        scf_ch.ostream.mute()

        scf_ch.maximum_overlap(
            molecule,
            basis,
            scf_gs.mol_orbs,
            [0, 1, 2, 3, 4],
            [1, 2, 3, 4],
        )
        scf_ch.compute(molecule, basis)

        assert abs(scf_ch.get_scf_energy() + 56.070929384622) < 1.0e-6

        if is_mpi_master():
            scf_h5 = Path(scf_ch.checkpoint_file)
            if scf_h5.is_file():
                scf_h5.unlink()
            scf_final_h5 = scf_h5.with_suffix('.tensors.h5')
            if scf_final_h5.is_file():
                scf_final_h5.unlink()

    def test_MOM_UDFT(self):

        water_xyz = """
        O       0.0000000000     0.0000000000     0.1178336003
        H      -0.7595754146    -0.0000000000    -0.4713344012
        H       0.7595754146     0.0000000000    -0.4713344012
        """

        molecule = Molecule.read_str(water_xyz)
        basis = MolecularBasis.read(molecule, '6-31G', ostream=None)

        # Ground state
        scf_gs = ScfRestrictedDriver()
        scf_gs.ostream.mute()
        scf_gs.xcfun = 'B3LYP'
        scf_gs.compute(molecule, basis)

        # Core-hole
        molecule.set_charge(1)
        molecule.set_multiplicity(2)
        scf_ch = ScfUnrestrictedDriver()
        scf_ch.ostream.mute()
        scf_ch.xcfun = 'B3LYP'

        scf_ch.maximum_overlap(
            molecule,
            basis,
            scf_gs.mol_orbs,
            [0, 1, 2, 3, 4],
            [1, 2, 3, 4],
        )
        scf_ch.compute(molecule, basis)

        assert abs(scf_ch.get_scf_energy() + 56.4455719055) < 1.0e-6

        if is_mpi_master():
            scf_h5 = Path(scf_ch.checkpoint_file)
            if scf_h5.is_file():
                scf_h5.unlink()
            scf_final_h5 = scf_h5.with_suffix('.tensors.h5')
            if scf_final_h5.is_file():
                scf_final_h5.unlink()
