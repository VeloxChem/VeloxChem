from pathlib import Path
import numpy as np
import pytest
from mpi4py import MPI

from veloxchem.veloxchemlib import mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.orbitallocalization import OrbitalLocalizationDriver
from veloxchem.outputstream import OutputStream
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.errorhandler import VeloxChemError


class TestOrbitalLocalization:

    @staticmethod
    def _build_system(*, ri):
        h2o_xyz = """3
        xyz
        O   0    0     0
        H   0.2774    0.8929     0.2544
        H   0.6068    -0.2383     -0.7169
        """
        molecule = Molecule.from_xyz_string(h2o_xyz)
        basis = MolecularBasis.read(molecule, "def2-svp")
        scf_drv = ScfRestrictedDriver()
        scf_drv.ri = ri
        scf_drv.ostream.mute()
        scf_res = scf_drv.compute(molecule, basis)

        if scf_drv.rank == mpi_master():
            if ri:
                assert scf_res['ri_jk'] or scf_res['ri_coulomb']
                assert not scf_res['ri_fallback']

        n_occ = molecule.number_of_alpha_occupied_orbitals(basis)

        return molecule, basis, scf_res, n_occ

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
            ("boys", None, "orbloc_boys_C"),
            ("pm", "mulliken", "orbloc_pm_mulliken_C"),
            ("pm", "lowdin", "orbloc_pm_lowdin_C"),
        ],
    )
    @pytest.mark.parametrize('ri', [False, True])
    def test_localization(self, method, pm_projector, ref_name, ri):
        molecule, basis, scf_res, n_occ = self._build_system(ri=ri)

        loc = OrbitalLocalizationDriver()
        loc.method = method
        if pm_projector is not None:
            loc.pm_projector = pm_projector
        loc.ostream.mute()

        for mo_range in [None, (1, n_occ)]:
            C_loc = loc.compute(molecule, basis, scf_res, mo_range=mo_range)
            C_loc = C_loc["loc_orbs"].alpha_to_numpy()

            # Load reference
            ref_path = Path(__file__).parent / "data" / {
                False: f"{ref_name}.npy",
                True: f"{ref_name}_ri.npy"
            }[ri]
            C_ref_occ = np.load(ref_path)
            C_ref = np.zeros(C_loc.shape)
            C_ref[:, :n_occ] = C_ref_occ[:, :]

            # Align phases
            C_loc = self._align_phases(C_ref, C_loc)

            # Compare
            np.testing.assert_allclose(C_loc, C_ref, atol=1e-6)

    @pytest.mark.parametrize("method", ["boys", "pm"])
    @pytest.mark.parametrize('ri', [False, True])
    def test_compute_does_not_mutate_scf_orbitals(self, method, ri):
        molecule, basis, scf_res, n_occ = self._build_system(ri=ri)
        C_alpha = scf_res["C_alpha"].copy()

        loc = OrbitalLocalizationDriver()
        loc.method = method
        loc.ostream.mute()
        loc.compute(molecule, basis, scf_res, mo_range=(1, n_occ))

        np.testing.assert_allclose(scf_res["C_alpha"], C_alpha, atol=0.0)

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason="pytest.raises only valid in serial")
    @pytest.mark.parametrize('ri', [False, True])
    def test_mo_range_requires_uniform_occupations(self, ri):
        molecule, basis, scf_res, n_occ = self._build_system(ri=ri)

        loc = OrbitalLocalizationDriver()
        loc.ostream.mute()

        with pytest.raises(VeloxChemError, match="same occupation number"):
            loc.compute(molecule, basis, scf_res,
                        mo_range=(n_occ, n_occ + 1))

    @pytest.mark.parametrize("method", ["boys", "pm"])
    @pytest.mark.parametrize('ri', [False, True])
    def test_compute_uses_configured_ostream(self, method, capsys, ri):
        molecule, basis, scf_res, n_occ = self._build_system(ri=ri)

        loc = OrbitalLocalizationDriver(OutputStream(None))
        loc.method = method
        loc.compute(molecule, basis, scf_res, mo_range=(1, n_occ))

        captured = capsys.readouterr()
        assert captured.out == ""
