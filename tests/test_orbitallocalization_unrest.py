from pathlib import Path
import numpy as np
import pytest

from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.orbitallocalization import OrbitalLocalizationDriver
from veloxchem.scfunrestdriver import ScfUnrestrictedDriver


class TestOrbitalLocalizationUnrestricted:

    @staticmethod
    def _build_system():
        h2o_xyz = """3
        xyz
        O   0    0     0
        H   0.2774    0.8929     0.2544
        H   0.6068    -0.2383     -0.7169
        """
        molecule = Molecule.from_xyz_string(h2o_xyz)
        molecule.set_multiplicity(3)
        basis = MolecularBasis.read(molecule, "def2-svp")
        scf_drv = ScfUnrestrictedDriver()
        scf_drv.ostream.mute()
        scf_res = scf_drv.compute(molecule, basis)

        n_occ_a = molecule.number_of_alpha_occupied_orbitals(basis)
        n_occ_b = molecule.number_of_beta_occupied_orbitals(basis)

        return molecule, basis, scf_res, n_occ_a, n_occ_b

    def _align_phases(self, C_ref, C_test):
        # This is indeed necessary, since even though the representation
        # between the MOs is "identical", there can still be a global phase
        # for each individual MO.
        C_aligned = C_test.copy()
        unmatched = set(range(C_ref.shape[1]))

        for i in range(C_ref.shape[1]):
            for j in list(unmatched):
                if np.max(np.abs(C_ref[:, i] - C_test[:, j])) < 1e-6:
                    C_aligned[:, i] = C_test[:, j]
                    unmatched.remove(j)
                    break
                elif np.max(np.abs(C_ref[:, i] + C_test[:, j])) < 1e-6:
                    C_aligned[:, i] = -C_test[:, j]
                    unmatched.remove(j)
                    break

        return C_aligned

    @pytest.mark.parametrize(
        "method, pm_projector, ref_name",
        [
            ("boys", None, "orbloc_boys_C_unrest"),
            ("pm", "mulliken", "orbloc_pm_mulliken_C_unrest"),
            ("pm", "lowdin", "orbloc_pm_lowdin_C_unrest"),
        ],
    )
    def test_localization(self, method, pm_projector, ref_name):
        molecule, basis, scf_res, n_occ_a, n_occ_b = self._build_system()

        loc = OrbitalLocalizationDriver()
        loc.method = method
        if pm_projector is not None:
            loc.pm_projector = pm_projector
        loc.ostream.mute()

        for mo_range in [None, (1, n_occ_a, 1, n_occ_b)]:
            C_loc = loc.compute(molecule, basis, scf_res, mo_range=mo_range)
            C_loc_a = C_loc["loc_orbs"].alpha_to_numpy()
            C_loc_b = C_loc["loc_orbs"].beta_to_numpy()

            # Load reference
            ref_path_a = Path(__file__).parent / "data" / f"{ref_name}_alpha.npy"
            ref_path_b = Path(__file__).parent / "data" / f"{ref_name}_beta.npy"
            C_ref_occ_a = np.load(ref_path_a)
            C_ref_occ_b = np.load(ref_path_b)
            C_ref_a = np.zeros(C_loc_a.shape)
            C_ref_b = np.zeros(C_loc_b.shape)
            C_ref_a[:, :n_occ_a] = C_ref_occ_a[:, :]
            C_ref_b[:, :n_occ_b] = C_ref_occ_b[:, :]

            # Align phases
            C_loc_a = self._align_phases(C_ref_a, C_loc_a)
            C_loc_b = self._align_phases(C_ref_b, C_loc_b)

            # Compare
            np.testing.assert_allclose(C_loc_a, C_ref_a, atol=1e-6)
            np.testing.assert_allclose(C_loc_b, C_ref_b, atol=1e-6)
