import numpy as np
import pytest
from pathlib import Path

from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.orbitallocalization import OrbitalLocalizationDriver


@pytest.mark.solvers
class TestOrbitalLocalization:

    @staticmethod
    def _build_system():
        h2o_xyz = """3
        xyz
        O   0    0     0
        H   0.2774    0.8929     0.2544
        H   0.6068    -0.2383     -0.7169
        """
        molecule = Molecule.from_xyz_string(h2o_xyz)
        basis = MolecularBasis.read(molecule, "def2-svp")

        # Load reference
        C_occ_path = Path(__file__).parent / "data" / "orbloc_h2o_svp_C_occ.npy"
        C_occ = np.load(C_occ_path)

        return molecule, basis, C_occ

    @staticmethod
    def _align_phases(C_ref, C_test):
        # This is indeed necessary, since even though the representation
        # between the MOs is "identical", there can still be a global phase
        # for each individual MO.
        C_aligned = C_test.copy()
        for i in range(C_ref.shape[1]):
            overlap = np.dot(C_ref[:, i], C_test[:, i])
            if overlap < 0:
                C_aligned[:, i] *= -1.0
        return C_aligned

    def test_boys(self):
        molecule, basis, C = self._build_system()

        loc = OrbitalLocalizationDriver()
        loc.silent = True
        C_loc = loc.boys(molecule, basis, C.copy())

        # Load reference
        ref_path = Path(__file__).parent / "data" / "orbloc_boys_C.npy"
        C_ref = np.load(ref_path)

        # Align phases
        C_loc = self._align_phases(C_ref, C_loc)

        # Compare
        np.testing.assert_allclose(C_loc, C_ref, atol=1e-6)

    def test_pipek_mezey_mulliken(self):
        molecule, basis, C = self._build_system()

        loc = OrbitalLocalizationDriver()
        loc.silent = True
        C_loc = loc.pipek_mezey(molecule, basis, C.copy(), projector="mulliken")

        # Load reference
        ref_path = Path(__file__).parent / "data" / "orbloc_pm_mulliken_C.npy"
        C_ref = np.load(ref_path)

        # Align phases
        C_loc = self._align_phases(C_ref, C_loc)

        # Compare
        np.testing.assert_allclose(C_loc, C_ref, atol=1e-6)

    def test_pipek_mezey_lowdin(self):
        molecule, basis, C = self._build_system()

        loc = OrbitalLocalizationDriver()
        loc.silent = True
        C_loc = loc.pipek_mezey(molecule, basis, C.copy(), projector="lowdin")

        # Load reference
        ref_path = Path(__file__).parent / "data" / "orbloc_pm_lowdin_C.npy"
        C_ref = np.load(ref_path)

        # Align phases
        C_loc = self._align_phases(C_ref, C_loc)

        # Compare
        np.testing.assert_allclose(C_loc, C_ref, atol=1e-6)
