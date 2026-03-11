from pathlib import Path

import numpy as np
import pytest
from mpi4py import MPI

from veloxchem.veloxchemlib import mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.dispersionmodel import DispersionModel


@pytest.mark.solvers
class TestScfDriverMiscellaneous:
    """
    Supplemental real-SCF coverage for compute() branches not already covered
    in the existing solver tests.

    Existing real tests already cover:
    - default min_basis with and without ECP:
      tests/test_scf_rest.py, tests/test_scf_unrest.py
    - RI-J and RI-JK success paths:
      tests/test_ri_scf.py, tests/test_rijk_scf.py
    - CPCM, SMD, PE, and point charges:
      tests/test_cpcm_solvation_energy.py, tests/test_smd_solvation_energy.py,
      tests/test_embedding.py, tests/test_pointcharges.py
    """

    @staticmethod
    def get_water_and_basis():

        xyz_string = """3
        xyz
        O   -0.1858140  -1.1749469   0.7662596
        H   -0.1285513  -0.8984365   1.6808606
        H   -0.0582782  -0.3702550   0.2638279
        """

        molecule = Molecule.read_xyz_string(xyz_string)
        basis = MolecularBasis.read(molecule, "sto-3g", ostream=None)

        return molecule, basis

    @staticmethod
    def run_hf_scf(molecule, basis, configure=None):

        scf_drv = ScfRestrictedDriver()
        scf_drv.ostream.mute()

        if configure is not None:
            configure(scf_drv)

        scf_results = scf_drv.compute(molecule, basis)

        return scf_drv, scf_results

    @staticmethod
    def is_master():

        return MPI.COMM_WORLD.Get_rank() == mpi_master()

    def test_filename_checkpoint_and_restart(self, tmp_path):

        molecule, basis = self.get_water_and_basis()
        filename = str(tmp_path / "water_restart")
        checkpoint_file = Path(f"{filename}_scf.h5")

        first_drv, first_results = self.run_hf_scf(
            molecule, basis, lambda drv: setattr(drv, "filename", filename))

        assert first_results is not None
        assert first_drv.checkpoint_file == str(checkpoint_file)
        if self.is_master():
            assert checkpoint_file.is_file()

        second_drv, second_results = self.run_hf_scf(
            molecule, basis, lambda drv: setattr(drv, "filename", filename))

        assert second_results is not None
        assert second_drv.checkpoint_file == str(checkpoint_file)
        assert second_drv.restart
        if self.is_master():
            assert second_drv._ref_mol_orbs is not None
            assert second_results["scf_energy"] == pytest.approx(
                first_results["scf_energy"], abs=1.0e-10)

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason='skip pytest.raises for multiple MPI processes')
    def test_conflicting_ri_modes_raise(self):

        molecule, basis = self.get_water_and_basis()

        def configure(driver):
            driver.ri_coulomb = True
            driver.ri_jk = True

        with pytest.raises(AssertionError, match="either ri_coulomb or ri_jk"):
            self.run_hf_scf(molecule, basis, configure)

    def test_grid_level_consistency(self):

        molecule, basis = self.get_water_and_basis()

        def configure(driver):
            driver.xcfun = "slda"
            driver.grid_level = 2

        scf_drv, scf_results = self.run_hf_scf(molecule, basis, configure)

        assert scf_results is not None
        assert scf_drv._mol_grid is not None
        if self.is_master():
            assert scf_results["grid_level"] == 2
            assert scf_results["xcfun"].lower() == "slda"
            assert np.isfinite(scf_results["scf_energy"])

    @pytest.mark.skipif(not DispersionModel.is_available(),
                        reason='dftd4-python not available')
    def test_d4_correction_consistency(self):

        molecule, basis = self.get_water_and_basis()

        ref_drv, ref_results = self.run_hf_scf(molecule, basis)

        def configure(driver):
            driver.dispersion = True

        d4_drv, d4_results = self.run_hf_scf(molecule, basis, configure)

        assert ref_results is not None
        assert d4_results is not None
        assert d4_drv._d4_energy != pytest.approx(0.0, abs=1.0e-12)
        if self.is_master():
            assert (d4_results["scf_energy"] -
                    ref_results["scf_energy"]) == pytest.approx(
                        d4_drv._d4_energy, abs=1.0e-10)

    @pytest.mark.parametrize("requested,expected", [(0, 1), (4, 3)])
    def test_print_level_consistency(self, requested, expected):

        molecule, basis = self.get_water_and_basis()

        def configure(driver):
            driver.print_level = requested

        scf_drv, scf_results = self.run_hf_scf(molecule, basis, configure)

        assert scf_results is not None
        assert scf_drv.print_level == expected

    def test_compute_return_none_not_converged(self):

        molecule, basis = self.get_water_and_basis()

        def configure(driver):
            driver.max_iter = 1
            driver.conv_thresh = 1.0e-12

        scf_drv, scf_results = self.run_hf_scf(molecule, basis, configure)

        assert scf_results is None
        assert not scf_drv.is_converged
