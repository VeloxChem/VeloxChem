import numpy as np
import pytest
from pathlib import Path

from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.orbitallocalization import OrbitalLocalizationDriver
from veloxchem.scfrestdriver import ScfRestrictedDriver


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
        scf_drv = ScfRestrictedDriver()
        scf_drv.ostream.mute()
        scf_res = scf_drv.compute(molecule, basis)

        n_occ = molecule.number_of_alpha_occupied_orbitals(basis)

        return molecule, basis, scf_res, n_occ

    @staticmethod
    def _align_phases(C_ref, C_test):
        # This is indeed necessary, since even though the representation
        # between the MOs is "identical", there can still be a global phase
        # for each individual MO.
        C_aligned = C_test.copy()
        mos_to_swap = {}

        for i in range(C_ref.shape[1]):
            overlap = np.dot(C_ref[:, i], C_test[:, i])
            test_squared = np.dot(C_test[:, i], C_test[:, i])
            if abs(overlap - test_squared) < 1e-6:
                continue
            elif abs(overlap + test_squared) < 1e-6:
                C_aligned[:, i] *= -1.0
            else:
                # check if degenerate MOs got swapped
                for j in range(i+1, C_ref.shape[1]):
                    overlap_ij = np.dot(C_ref[:, j], C_test[:, i])
                    if abs(overlap_ij - test_squared) < 1e-6:
                        mos_to_swap[(i, j)] = 1
                        break
                    elif abs(overlap_ij + test_squared) < 1e-6:
                        mos_to_swap[(i, j)] = 1
                        C_aligned[:, i] *= -1.0
                        break
                    else:
                        continue
        
        if mos_to_swap:
            for pair in mos_to_swap:
                C_aligned[:, [pair[0], pair[1]]] = C_aligned[:, [pair[1], pair[0]]]
            # catch remaining phase conventions from swapped MOs
            TestOrbitalLocalization._align_phases(C_ref, C_aligned)
        return C_aligned

    def test_boys(self):
        molecule, basis, scf_res, n_occ = self._build_system()

        loc = OrbitalLocalizationDriver()
        loc.method = "boys"
        loc.silent = True
        C_loc = loc.compute(molecule, basis, scf_res, mo_range=(1, n_occ))
        C_loc = C_loc["loc_orbs"].alpha_to_numpy()

        # Load reference
        ref_path = Path(__file__).parent / "data" / "orbloc_boys_C.npy"
        C_ref = np.load(ref_path)

        # Align phases
        C_loc = self._align_phases(C_ref, C_loc)

        # Compare
        np.testing.assert_allclose(C_loc, C_ref, atol=1e-6)

    def test_pipek_mezey_mulliken(self):
        molecule, basis, scf_res, n_occ = self._build_system()

        loc = OrbitalLocalizationDriver()
        loc.method = "pm"
        loc.pm_projector = "mulliken"
        loc.silent = True
        C_loc = loc.compute(molecule, basis, scf_res, mo_range=(1, n_occ))
        C_loc = C_loc["loc_orbs"].alpha_to_numpy()

        # Load reference
        ref_path = Path(__file__).parent / "data" / "orbloc_pm_mulliken_C.npy"
        C_ref = np.load(ref_path)

        # Align phases
        C_loc = self._align_phases(C_ref, C_loc)

        # Compare
        np.testing.assert_allclose(C_loc, C_ref, atol=1e-6)

    def test_pipek_mezey_lowdin(self):
        molecule, basis, scf_res, n_occ = self._build_system()

        loc = OrbitalLocalizationDriver()
        loc.method = "pm"
        loc.pm_projector = "lowdin"
        loc.silent = True
        C_loc = loc.compute(molecule, basis, scf_res, mo_range=(1, n_occ))
        C_loc = C_loc["loc_orbs"].alpha_to_numpy()

        # Load reference
        ref_path = Path(__file__).parent / "data" / "orbloc_pm_lowdin_C.npy"
        C_ref = np.load(ref_path)

        # Align phases
        C_loc = self._align_phases(C_ref, C_loc)

        # Compare
        np.testing.assert_allclose(C_loc, C_ref, atol=1e-6)
