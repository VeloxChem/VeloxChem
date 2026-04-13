from pathlib import Path
from copy import deepcopy

import numpy as np
import pytest
from mpi4py import MPI

from veloxchem.veloxchemlib import mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.molecularorbitals import molorb
from veloxchem.outputstream import OutputStream
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.scfrestopendriver import ScfRestrictedOpenDriver
from veloxchem.scfunrestdriver import ScfUnrestrictedDriver
from veloxchem.dispersionmodel import DispersionModel
from veloxchem.checkpoint import read_molecule_and_basis
from veloxchem.inputparser import unparse_input, read_unparsed_input_from_hdf5


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
    def get_open_shell_water_and_basis(charge=1, multiplicity=2):

        xyz_string = """3
        xyz
        O   -0.1858140  -1.1749469   0.7662596
        H   -0.1285513  -0.8984365   1.6808606
        H   -0.0582782  -0.3702550   0.2638279
        """

        molecule = Molecule.read_xyz_string(xyz_string)
        molecule.set_charge(charge)
        molecule.set_multiplicity(multiplicity)
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

        # To avoid inconsistency across MPI ranks
        comm = MPI.COMM_WORLD
        filename = comm.bcast(filename, root=mpi_master())
        checkpoint_file = comm.bcast(checkpoint_file, root=mpi_master())

        first_drv, first_results = self.run_hf_scf(
            molecule, basis, lambda drv: setattr(drv, "filename", filename))

        assert first_results is not None
        assert first_drv.checkpoint_file is None
        if self.is_master():
            assert checkpoint_file.is_file()

        second_drv, second_results = self.run_hf_scf(
            molecule, basis, lambda drv: setattr(drv, "filename", filename))

        assert second_results is not None
        assert second_drv.checkpoint_file is None
        assert second_drv.restart
        if self.is_master():
            assert second_drv._ref_mol_orbs is not None
            assert second_results["scf_energy"] == pytest.approx(
                first_results["scf_energy"], abs=1.0e-10)

        third_drv = ScfRestrictedDriver()
        third_drv.ostream.mute()

        # reconstruct molecule and basis
        if self.is_master():
            new_molecule, new_basis = read_molecule_and_basis(
                str(checkpoint_file))
        else:
            new_molecule, new_basis = None, None
        new_molecule = third_drv.comm.bcast(new_molecule, root=mpi_master())
        new_basis = third_drv.comm.bcast(new_basis, root=mpi_master())
        third_drv, third_results = self.run_hf_scf(
            new_molecule, new_basis,
            lambda drv: setattr(drv, "filename", filename))

        assert third_results is not None
        assert third_drv.checkpoint_file is None
        assert third_drv.restart
        if self.is_master():
            assert third_drv._ref_mol_orbs is not None
            assert third_results["scf_energy"] == pytest.approx(
                first_results["scf_energy"], abs=1.0e-10)

    def test_explicit_checkpoint_file_overrides_filename(self, tmp_path):

        molecule, basis = self.get_water_and_basis()
        filename = str(tmp_path / "water_named_output")
        checkpoint_file = str(tmp_path / "custom_restart_file.h5")

        comm = MPI.COMM_WORLD
        filename = comm.bcast(filename, root=mpi_master())
        checkpoint_file = comm.bcast(checkpoint_file, root=mpi_master())

        def configure(driver):
            driver.filename = filename
            driver.checkpoint_file = checkpoint_file

        first_drv, first_results = self.run_hf_scf(molecule, basis, configure)

        assert first_results is not None
        assert first_drv.checkpoint_file == checkpoint_file
        if self.is_master():
            assert Path(checkpoint_file).is_file()
            assert not Path(f"{filename}_scf.h5").exists()

        second_drv, second_results = self.run_hf_scf(molecule, basis, configure)

        assert second_results is not None
        assert second_drv.checkpoint_file == checkpoint_file
        assert second_drv.restart
        if self.is_master():
            assert second_drv._ref_mol_orbs is not None
            assert second_results["scf_energy"] == pytest.approx(
                first_results["scf_energy"], abs=1.0e-10)

    def test_no_filename_or_checkpoint_file_disables_checkpoint_output(self):

        molecule, basis = self.get_water_and_basis()

        scf_drv, scf_results = self.run_hf_scf(molecule, basis)

        assert scf_results is not None
        assert scf_drv.checkpoint_file is None
        assert scf_drv.restart is False

    def test_checkpoint_writes_input_groups(self, tmp_path):

        molecule, basis = self.get_water_and_basis()
        filename = str(tmp_path / "water_keywords")
        checkpoint_file = Path(f"{filename}_scf.h5")

        comm = MPI.COMM_WORLD
        filename = comm.bcast(filename, root=mpi_master())
        checkpoint_file = comm.bcast(checkpoint_file, root=mpi_master())

        def configure(driver):
            driver.filename = filename
            driver.max_iter = 37
            driver.density_damping = True
            driver.xcfun = 'blyp'
            driver.ri_coulomb = True
            driver.ri_auxiliary_basis = 'def2-universal-jkfit'
            driver.ri_metric_threshold = 1.0e-10
            driver.solvation_model = 'cpcm'
            driver.cpcm_custom_vdw_radii = ('Mg', '2.0', '1', '1.0')

        scf_drv, scf_results = self.run_hf_scf(molecule, basis, configure)

        assert scf_results is not None
        assert scf_drv.checkpoint_file is None

        scf_keywords = {
            key: val[0] for key, val in scf_drv._input_keywords["scf"].items()
        }
        method_keywords = {
            key: val[0]
            for key, val in scf_drv._input_keywords["method_settings"].items()
        }

        expected_scf = unparse_input(scf_drv, scf_keywords)
        expected_method = unparse_input(scf_drv, method_keywords)

        # test read_unparsed_input_from_hdf5

        if self.is_master():
            assert checkpoint_file.is_file()
            checkpoint_scf_input = read_unparsed_input_from_hdf5(
                str(checkpoint_file), group_name="scf_settings")
            checkpoint_method_input = read_unparsed_input_from_hdf5(
                str(checkpoint_file), group_name="method_settings")
        else:
            checkpoint_scf_input = None
            checkpoint_method_input = None

        checkpoint_scf_input = scf_drv.comm.bcast(checkpoint_scf_input,
                                                  root=mpi_master())
        checkpoint_method_input = scf_drv.comm.bcast(checkpoint_method_input,
                                                     root=mpi_master())

        assert checkpoint_scf_input == expected_scf
        assert checkpoint_method_input == expected_method

        # test update_settings

        second_drv = ScfRestrictedDriver()
        second_drv.ostream.mute()
        second_drv.update_settings(checkpoint_scf_input,
                                   checkpoint_method_input)

        for key in scf_keywords:
            assert getattr(second_drv, key) == getattr(scf_drv, key)
        for key in method_keywords:
            assert getattr(second_drv, key) == getattr(scf_drv, key)

        # test read_settings

        third_drv = ScfRestrictedDriver()
        third_drv.ostream.mute()
        third_drv.checkpoint_file = str(checkpoint_file)
        third_drv.filename = str(tmp_path / "copied_settings_target")

        new_molecule, new_basis = read_molecule_and_basis(str(checkpoint_file))

        third_drv.read_settings(str(checkpoint_file))
        third_results = third_drv.compute(new_molecule, new_basis)

        for key in scf_keywords:
            if key in ('restart', 'filename', 'checkpoint_file'):
                continue
            assert getattr(third_drv, key) == getattr(scf_drv, key)
        for key in method_keywords:
            assert getattr(third_drv, key) == getattr(scf_drv, key)

        assert third_drv.restart
        assert third_drv.filename == str(tmp_path / "copied_settings_target")
        assert third_drv.checkpoint_file == str(checkpoint_file)

        if self.is_master():
            assert third_results['scf_type'] == scf_results['scf_type']
            assert abs(third_results['scf_energy'] -
                       scf_results['scf_energy']) < 1e-12
            assert np.max(
                np.abs(third_results['D_alpha'] -
                       scf_results['D_alpha'])) < 1e-8

        # test read_settings with a different scf_type and invalid restart

        fourth_drv = ScfUnrestrictedDriver()
        fourth_drv.ostream.mute()
        fourth_drv.checkpoint_file = str(checkpoint_file)
        fourth_drv.read_settings(str(checkpoint_file))

        new_molecule, new_basis = read_molecule_and_basis(str(checkpoint_file))
        new_molecule.set_charge(1)
        new_molecule.set_multiplicity(2)

        fourth_results_not_used = fourth_drv.compute(new_molecule, new_basis)

        assert not fourth_drv.restart

    def test_read_settings_imports_only_configuration(self, tmp_path):

        molecule, basis = self.get_water_and_basis()
        filename = str(tmp_path / "water_import_settings")
        checkpoint_file = Path(f"{filename}_scf.h5")

        comm = MPI.COMM_WORLD
        filename = comm.bcast(filename, root=mpi_master())
        checkpoint_file = comm.bcast(checkpoint_file, root=mpi_master())

        def configure(driver):
            driver.filename = filename
            driver.max_iter = 37
            driver.density_damping = True
            driver.xcfun = 'blyp'
            driver.ri_coulomb = True
            driver.ri_auxiliary_basis = 'def2-universal-jkfit'
            driver.ri_metric_threshold = 1.0e-10

        scf_drv, _ = self.run_hf_scf(molecule, basis, configure)

        imported_drv = ScfRestrictedDriver()
        imported_drv.ostream.mute()
        imported_drv.restart = False
        imported_drv.filename = str(tmp_path / "do_not_override")
        imported_drv.checkpoint_file = str(tmp_path / "do_not_override_scf.h5")

        imported_drv.read_settings(str(checkpoint_file))

        assert imported_drv.restart is False
        assert imported_drv.filename == str(tmp_path / "do_not_override")
        assert imported_drv.checkpoint_file == str(tmp_path /
                                                   "do_not_override_scf.h5")
        assert imported_drv.max_iter == scf_drv.max_iter
        assert imported_drv.density_damping == scf_drv.density_damping
        assert imported_drv.xcfun == scf_drv.xcfun
        assert imported_drv.ri_coulomb == scf_drv.ri_coulomb
        assert imported_drv.ri_auxiliary_basis == scf_drv.ri_auxiliary_basis
        assert imported_drv.ri_metric_threshold == scf_drv.ri_metric_threshold

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

    def test_scf_results_being_independent(self):

        molecule, basis = self.get_water_and_basis()

        scf_drv, scf_results = self.run_hf_scf(molecule, basis)

        if self.is_master():
            ref_energy = scf_results["scf_energy"]
            scf_results["scf_energy"] = 0.0

            ref_C_alpha = scf_results["C_alpha"].copy()
            scf_results["C_alpha"][:, :] = 0.0

            assert scf_results is not scf_drv.scf_results
            assert scf_results["C_alpha"] is not scf_drv.scf_results["C_alpha"]

            assert ref_energy == scf_drv.scf_results["scf_energy"]
            assert 0.0 == np.max(
                np.abs(ref_C_alpha - scf_drv.scf_results["C_alpha"]))

            assert 0.0 == scf_results["scf_energy"]
            assert 0.0 == np.max(np.abs(scf_results["C_alpha"]))

            assert 0.0 != scf_drv.scf_results["scf_energy"]
            assert 0.0 != np.max(np.abs(scf_drv.scf_results["C_alpha"]))

    def test_mom(self):

        molecule, basis = self.get_water_and_basis()

        scf_drv, scf_results = self.run_hf_scf(molecule, basis)

        uhf_drv = ScfUnrestrictedDriver()
        uhf_drv.ostream.mute()

        occ_beta = list(range(molecule.number_of_beta_occupied_orbitals(basis)))
        occ_alpha = list(occ_beta)
        occ_alpha[-1] += 1  # HOMO->LUMO excitation

        uhf_drv.maximum_overlap(molecule, basis, scf_drv.molecular_orbitals,
                                occ_alpha, occ_beta)

        scf_results = uhf_drv.compute(molecule, basis)

        if self.is_master():
            assert abs(-74.54939063506086 - scf_results["scf_energy"]) < 1.0e-8

    def test_user_supplied_start_orbitals_are_not_checkpoint_restart(
            self, tmp_path):

        molecule, basis = self.get_water_and_basis()
        checkpoint_file = tmp_path / "user_start_should_not_exist.h5"

        ref_drv, ref_results = self.run_hf_scf(molecule, basis)

        start_drv = ScfRestrictedDriver()
        start_drv.ostream.mute()
        start_drv.checkpoint_file = str(checkpoint_file)

        if self.is_master():
            start_orbitals = ref_drv.molecular_orbitals.alpha_to_numpy().copy()
        else:
            start_orbitals = None

        start_drv.set_start_orbitals(molecule, basis, start_orbitals)
        if self.is_master():
            assert not checkpoint_file.exists()
        start_results = start_drv.compute(molecule, basis)

        assert start_drv.restart is False
        assert start_drv._use_start_orbitals is True
        if self.is_master():
            assert start_drv._ref_mol_orbs is not None
            assert checkpoint_file.exists()
            assert start_results["scf_energy"] == pytest.approx(
                ref_results["scf_energy"], abs=1.0e-10)

    def test_clear_start_orbitals_resets_user_start_mode(self):

        molecule, basis = self.get_water_and_basis()

        ref_drv, _ = self.run_hf_scf(molecule, basis)

        start_drv = ScfRestrictedDriver()
        start_drv.ostream.mute()

        if self.is_master():
            start_orbitals = ref_drv.molecular_orbitals.alpha_to_numpy().copy()
        else:
            start_orbitals = None

        start_drv.set_start_orbitals(molecule, basis, start_orbitals)
        assert start_drv.restart is False
        assert start_drv._use_start_orbitals is True

        start_drv._mom = ('alpha', 'beta')
        start_drv.clear_start_orbitals()

        assert start_drv.restart is False
        assert start_drv._use_start_orbitals is False
        assert start_drv._mom is None

    def test_scfdriver_deepcopy(self):

        molecule, basis = self.get_water_and_basis()

        scf_drv, scf_results = self.run_hf_scf(
            molecule, basis, lambda drv: setattr(drv, "xcfun", "pbe"))

        scf_drv_copy = deepcopy(scf_drv)

        assert scf_drv_copy.xcfun == scf_drv.xcfun
        assert scf_drv_copy.scf_energy == scf_drv.scf_energy
        if self.is_master():
            assert np.allclose(scf_drv_copy.scf_results['D_alpha'],
                               scf_drv.scf_results['D_alpha'])

    def test_restricted_driver_helper_branches(self):

        driver = ScfRestrictedDriver()
        driver.ostream.mute()

        ovl = np.eye(2)
        oao = np.eye(2)
        fock_a = np.array([[1.0, 0.2], [0.2, 0.4]])
        fock_b = np.array([[0.8, 0.1], [0.1, 0.3]])
        density_a = np.array([[1.0, 0.0], [0.0, 0.0]])
        density_b = np.array([[0.85, 0.0], [0.0, 0.15]])

        effective = driver._get_effective_fock((fock_a,), ovl, oao)
        if self.is_master():
            assert np.allclose(effective[0], fock_a)
            assert effective[0] is fock_a

        driver.acc_type = 'c2diis'
        driver.diis_thresh = 1.0
        driver.max_err_vecs = 2
        driver._store_diis_data((fock_a,), (density_a,), ovl, 0.2)

        single_effective = driver._get_effective_fock((fock_b,), ovl, oao)
        if self.is_master():
            assert np.allclose(single_effective[0], fock_a)
            assert single_effective[0] is not driver._fock_matrices_alpha[0]

        driver._store_diis_data((fock_b,), (density_b,), ovl, 0.2)

        diis_effective = driver._get_effective_fock((fock_b,), ovl, oao)
        if self.is_master():
            assert diis_effective[0].shape == fock_a.shape
            assert np.allclose(diis_effective[0], diis_effective[0].T)

        driver.embedding = {'settings': {'embedding_method': 'PE'}}
        assert driver.get_scf_type_str() == (
            'Spin-Restricted Hartree-Fock with PE')

        driver._dft = True
        assert driver.get_scf_type_str() == 'Spin-Restricted Kohn-Sham with PE'

    def test_unrestricted_helper_branches_and_natural_orbitals(self):

        molecule, basis = self.get_open_shell_water_and_basis()
        driver = ScfUnrestrictedDriver()
        driver.ostream.mute()

        ovl = np.eye(2)
        oao = np.eye(2)
        fock_a = np.array([[1.0, 0.2], [0.2, 0.5]])
        fock_b = np.array([[0.9, 0.1], [0.1, 0.3]])
        fock_c = np.array([[0.7, 0.05], [0.05, 0.2]])
        density_a = np.array([[1.0, 0.0], [0.0, 0.0]])
        density_b = np.array([[0.0, 0.0], [0.0, 1.0]])
        density_c = np.array([[0.8, 0.0], [0.0, 0.2]])
        density_d = np.array([[0.1, 0.0], [0.0, 0.9]])

        passthrough = driver._get_effective_fock((fock_a, fock_b), ovl, oao)
        if self.is_master():
            assert np.allclose(passthrough[0], fock_a)
            assert np.allclose(passthrough[1], fock_b)

        driver.acc_type = 'c2diis'
        driver.diis_thresh = 1.0
        driver.max_err_vecs = 2
        driver._store_diis_data((fock_a, fock_b), (density_a, density_b), ovl,
                                0.2)

        single_effective = driver._get_effective_fock((fock_c, fock_c), ovl,
                                                      oao)
        if self.is_master():
            assert np.allclose(single_effective[0], fock_a)
            assert np.allclose(single_effective[1], fock_b)

        driver._store_diis_data((fock_c, fock_c), (density_c, density_d), ovl,
                                0.2)
        driver._store_diis_data((fock_b, fock_a), (density_b, density_a), ovl,
                                0.2)
        if self.is_master():
            assert len(driver._fock_matrices_alpha) == 2
            assert np.allclose(driver._fock_matrices_alpha[0], fock_c)
            assert np.allclose(driver._fock_matrices_beta[0], fock_c)

        diis_effective = driver._get_effective_fock((fock_a, fock_b), ovl, oao)
        if self.is_master():
            assert diis_effective[0].shape == fock_a.shape
            assert diis_effective[1].shape == fock_b.shape

        n_orbitals = basis.get_dimension_of_basis()
        pfon_fock_a = np.diag(np.linspace(-1.5, 1.5, n_orbitals))
        pfon_fock_b = np.diag(np.linspace(-1.2, 1.8, n_orbitals))

        driver.pfon = True
        driver.pfon_temperature = 200000.0
        driver.pfon_nocc = 1
        driver.pfon_nvir = 1

        mol_orbs = driver._gen_molecular_orbitals(molecule, basis,
                                                  (pfon_fock_a, pfon_fock_b),
                                                  np.eye(n_orbitals))

        nocc_a = molecule.number_of_alpha_occupied_orbitals(basis)
        nocc_b = molecule.number_of_beta_occupied_orbitals(basis)

        if self.is_master():
            assert mol_orbs.get_orbitals_type() == molorb.unrest
            assert mol_orbs.occa_to_numpy().sum() == pytest.approx(nocc_a)
            assert mol_orbs.occb_to_numpy().sum() == pytest.approx(nocc_b)
            assert 0.0 < mol_orbs.occa_to_numpy()[nocc_a] < 1.0
            assert 0.0 < mol_orbs.occb_to_numpy()[nocc_b] < 1.0

        scf_results = driver.compute(molecule, basis)
        natural_orbitals = driver.natural_orbitals()

        assert scf_results is not None
        if self.is_master():
            assert natural_orbitals.get_orbitals_type() == molorb.rest
            occupations = natural_orbitals.occa_to_numpy()
            assert np.all(occupations[:-1] >= occupations[1:])
            assert np.all(occupations >= -1.0e-12)
            assert np.all(occupations <= 2.0 + 1.0e-12)

        driver.embedding = {'settings': {'embedding_method': 'PE'}}
        assert driver.get_scf_type_str() == (
            'Spin-Unrestricted Hartree-Fock with PE')

        driver._dft = True
        assert driver.get_scf_type_str() == (
            'Spin-Unrestricted Kohn-Sham with PE')

        copied = deepcopy(driver)
        assert copied.ostream is driver.ostream
        assert copied.comm is driver.comm
        assert copied.xcfun == driver.xcfun
        if self.is_master():
            assert np.allclose(copied.scf_results['D_alpha'],
                               driver.scf_results['D_alpha'])

    def test_restricted_open_helper_branches(self):

        molecule, basis = self.get_open_shell_water_and_basis()
        driver = ScfRestrictedOpenDriver()
        driver.ostream.mute()

        fa = np.array([[1.2, 0.3], [0.3, 0.4]])
        fb = np.array([[0.8, 0.1], [0.1, 0.2]])
        da = np.array([[1.0, 0.0], [0.0, 0.0]])
        db = np.array([[0.3, 0.0], [0.0, 0.1]])
        s = np.eye(2)

        projected = driver.get_projected_fock(fa, fb, da, db, s)
        f0 = 0.5 * (fa + fb)
        inactive = np.matmul(s, db)
        active = np.matmul(s, da - db)
        virtual = np.eye(2) - np.matmul(s, da)
        expected = f0 + np.linalg.multi_dot([inactive, fb - f0, active.T])
        expected += np.linalg.multi_dot([active, fa - f0, virtual.T])
        expected += (expected - f0).T

        assert np.allclose(projected, expected)
        assert np.allclose(projected, projected.T)

        passthrough = driver._get_effective_fock((fa, fb), s, s)
        if self.is_master():
            assert np.allclose(passthrough[0], fa)
            assert np.allclose(passthrough[1], fb)

        driver.acc_type = 'c2diis'
        driver.diis_thresh = 1.0
        driver.max_err_vecs = 2
        driver._store_diis_data((fa, fb), (da, db), s, 0.2)

        single_effective = driver._get_effective_fock((fa, fb), s, s)
        if self.is_master():
            assert len(single_effective) == 1
            assert np.allclose(single_effective[0], projected)

        fa_2 = np.array([[1.0, 0.2], [0.2, 0.6]])
        fb_2 = np.array([[0.7, 0.15], [0.15, 0.5]])
        da_2 = np.array([[0.9, 0.0], [0.0, 0.1]])
        db_2 = np.array([[0.2, 0.0], [0.0, 0.2]])
        driver._store_diis_data((fa_2, fb_2), (da_2, db_2), s, 0.2)

        diis_effective = driver._get_effective_fock((fa_2, fb_2), s, s)
        if self.is_master():
            assert len(diis_effective) == 1
            assert diis_effective[0].shape == fa.shape

        n_orbitals = basis.get_dimension_of_basis()
        pfon_fock = np.diag(np.linspace(-1.5, 1.5, n_orbitals))

        driver.pfon = True
        driver.pfon_temperature = 200000.0
        driver.pfon_nocc = 1
        driver.pfon_nvir = 1

        mol_orbs = driver._gen_molecular_orbitals(molecule, basis, (pfon_fock,),
                                                  np.eye(n_orbitals))

        nocc_a = molecule.number_of_alpha_occupied_orbitals(basis)
        nocc_b = molecule.number_of_beta_occupied_orbitals(basis)

        if self.is_master():
            assert mol_orbs.get_orbitals_type() == molorb.restopen
            assert mol_orbs.occa_to_numpy().sum() == pytest.approx(nocc_a)
            assert mol_orbs.occb_to_numpy().sum() == pytest.approx(nocc_b)
            assert 0.0 < mol_orbs.occa_to_numpy()[nocc_a] < 1.0
            assert 0.0 < mol_orbs.occb_to_numpy()[nocc_b] < 1.0

        driver.embedding = {'settings': {'embedding_method': 'PE'}}
        assert driver.get_scf_type_str() == (
            'Spin-Restricted Open-Shell Hartree-Fock with PE')

        driver._dft = True
        assert driver.get_scf_type_str() == (
            'Spin-Restricted Open-Shell Kohn-Sham with PE')

        copied = deepcopy(driver)
        assert copied.ostream is driver.ostream
        assert copied.comm is driver.comm
        assert copied._scf_type == driver._scf_type
        if self.is_master():
            assert np.allclose(copied._fock_matrices_proj[0],
                               driver._fock_matrices_proj[0])

    def test_guess_unpaired_electrons_warning_for_restricted(self, tmp_path):

        molecule, basis = self.get_water_and_basis()

        comm = MPI.COMM_WORLD
        ostream = OutputStream.create_mpi_ostream(
            comm, str(tmp_path / 'guess_restricted.out'))

        scf_drv = ScfRestrictedDriver(comm, ostream)
        scf_drv.guess_unpaired_electrons = '1(1.0)'
        scf_drv.max_iter = 2
        scf_drv.conv_thresh = 1.0e-12

        scf_results = scf_drv.compute(molecule, basis)

        assert scf_results is None
        if self.is_master():
            output = (tmp_path / 'guess_restricted.out').read_text()
            assert 'Ignoring "guess_unpaired_electrons" in spin-restricted SCF calculation.' in output

    def test_guess_unpaired_electrons_are_parsed_for_unrestricted(
            self, tmp_path):

        molecule, basis = self.get_water_and_basis()
        molecule.set_charge(1)
        molecule.set_multiplicity(2)

        comm = MPI.COMM_WORLD
        ostream = OutputStream.create_mpi_ostream(
            comm, str(tmp_path / 'guess_unrestricted.out'))

        scf_drv = ScfUnrestrictedDriver(comm, ostream)
        scf_drv.guess_unpaired_electrons = '1(1.0), 2(-0.5)'
        scf_drv.max_iter = 2
        scf_drv.conv_thresh = 1.0e-12

        scf_results = scf_drv.compute(molecule, basis)

        assert scf_results is None
        if self.is_master():
            output = (tmp_path / 'guess_unrestricted.out').read_text()
            assert 'Generating initial guess with user-provided information...' in output
            assert '1.0 unpaired alpha electrons on atom 1 (O)' in output
            assert '0.5 unpaired beta  electrons on atom 2 (H)' in output

    def test_level_shifting_and_pfon_real_scf(self, tmp_path):

        molecule, basis = self.get_water_and_basis()

        comm = MPI.COMM_WORLD
        ostream = OutputStream.create_mpi_ostream(
            comm, str(tmp_path / 'level_pfon.out'))

        scf_drv = ScfRestrictedDriver(comm, ostream)
        scf_drv.acc_type = 'diis'
        scf_drv.level_shifting = 0.2
        scf_drv.level_shifting_delta = 0.5
        scf_drv.pfon = True
        scf_drv.pfon_temperature = 100
        scf_drv.pfon_delta_temperature = 50
        scf_drv.max_iter = 2
        scf_drv.conv_thresh = 1.0e-12

        scf_results = scf_drv.compute(molecule, basis)

        assert scf_results is None
        assert len(scf_drv.history) == 3
        assert scf_drv.pfon_temperature == 0
        assert scf_drv.level_shifting == 0.0
        if self.is_master():
            output = (tmp_path / 'level_pfon.out').read_text()
            assert 'Applying level-shifting' in output
            assert 'Applying pseudo-FON' in output

    def test_density_damping_changes_real_scf_trajectory(self):

        molecule, basis = self.get_water_and_basis()

        ref_drv = ScfRestrictedDriver()
        ref_drv.ostream.mute()
        ref_drv.acc_type = 'diis'
        ref_drv.max_iter = 2
        ref_drv.conv_thresh = 1.0e-12
        ref_results = ref_drv.compute(molecule, basis)

        damped_drv = ScfRestrictedDriver()
        damped_drv.ostream.mute()
        damped_drv.acc_type = 'diis'
        damped_drv.max_iter = 2
        damped_drv.conv_thresh = 1.0e-12
        damped_drv.density_damping = True
        damped_results = damped_drv.compute(molecule, basis)

        assert ref_results is None
        assert damped_results is None
        assert len(ref_drv.history) == len(damped_drv.history) == 3

        if self.is_master():
            ref_energies = np.array(
                [step['energy'] for step in ref_drv.history])
            damped_energies = np.array(
                [step['energy'] for step in damped_drv.history])
            assert np.max(np.abs(ref_energies - damped_energies)) > 1.0e-4
            assert np.max(np.abs(ref_drv.density[0] -
                                 damped_drv.density[0])) > 1.0e-4

    def test_effective_nuclear_charges_subtract_ecp_core(self):

        class MoleculeStub:

            @staticmethod
            def get_element_ids():
                return np.array([8.0, 1.0])

            @staticmethod
            def number_of_atoms():
                return 2

        class BasisStub:

            @staticmethod
            def get_number_of_ecp_core_electrons():
                return np.array([2.0, 0.0])

        effective_charges = Molecule.get_effective_nuclear_charges(
            MoleculeStub(), BasisStub())

        assert np.allclose(effective_charges, np.array([6.0, 1.0]))
