import numpy as np
import pytest
import veloxchem as vlx
from pathlib import Path


@pytest.mark.solvers
class TestOrbitalLocalization:

    @staticmethod
    def _build_system():
        h2o_xyz = '''3

        O   0    0     0
        H   0.2774    0.8929     0.2544
        H   0.6068    -0.2383     -0.7169
        '''
        molecule = vlx.Molecule.from_xyz_string(h2o_xyz)
        basis = vlx.MolecularBasis.read(molecule, "def2-svp")

        # Load reference
        C_occ_path = Path(__file__).parent / "data" / "h2o_svp_C_occ.npy"
        C_occ = np.load(C_occ_path)
        S_path = Path(__file__).parent / "data" / "h2o_svp_S.npy"
        S = np.load(S_path)

        return molecule, basis, C_occ, S

    @staticmethod
    def _compute_dipoles(molecule, basis):
        coords = molecule.get_coordinates_in_bohr()
        nuclear_charges = molecule.get_element_ids()

        origin = np.sum(coords.T * nuclear_charges,
                        axis=1) / np.sum(nuclear_charges)

        dip_ints = vlx.compute_electric_dipole_integrals(
            molecule, basis, origin
        )

        return dip_ints

    @staticmethod
    def _get_atom_map(molecule, basis):
        atom_map_raw = basis.get_ao_basis_map(molecule)
        return [int(atom_map_raw[i].split()[0]) for i in range(len(atom_map_raw))]

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
        molecule, basis, C, S = self._build_system()
        dip_ints = self._compute_dipoles(molecule, basis)

        loc = vlx.OrbitalLocalization(C.copy())
        C_loc = loc.boys(dip_ints)

        # Load reference
        ref_path = Path(__file__).parent / "data" / "boys_C.npy"
        C_ref = np.load(ref_path)

        # Align phases
        C_loc = self._align_phases(C_ref, C_loc)

        # Compare
        np.testing.assert_allclose(C_loc, C_ref, atol=1e-6)

        # Orthonormality check
        I = np.eye(C.shape[1])
        np.testing.assert_allclose(C_loc.T @ S @ C_loc, I, atol=1e-8)

    def test_pipek_mezey_mulliken(self):
        molecule, basis, C, S = self._build_system()
        atom_map = self._get_atom_map(molecule, basis)

        loc = vlx.OrbitalLocalization(C.copy())
        C_loc = loc.pipek_mezey(
            S, atom_map,
            projector="mulliken"
        )

        # Load reference
        ref_path = Path(__file__).parent / "data" / "pm_mulliken_C.npy"
        C_ref = np.load(ref_path)

        # Align phases
        C_loc = self._align_phases(C_ref, C_loc)

        # Compare
        np.testing.assert_allclose(C_loc, C_ref, atol=1e-6)

        # Orthonormality check
        I = np.eye(C.shape[1])
        np.testing.assert_allclose(C_loc.T @ S @ C_loc, I, atol=1e-8)

    def test_pipek_mezey_lowdin(self):
        molecule, basis, C, S = self._build_system()
        atom_map = self._get_atom_map(molecule, basis)

        loc = vlx.OrbitalLocalization(C.copy())
        C_loc = loc.pipek_mezey(
            S, atom_map,
            projector="lowdin"
        )

        # Load reference
        ref_path = Path(__file__).parent / "data" / "pm_lowdin_C.npy"
        C_ref = np.load(ref_path)

        # Align phases
        C_loc = self._align_phases(C_ref, C_loc)

        # Compare
        np.testing.assert_allclose(C_loc, C_ref, atol=1e-6)

        # Orthonormality check
        I = np.eye(C.shape[1])
        np.testing.assert_allclose(C_loc.T @ S @ C_loc, I, atol=1e-8)

