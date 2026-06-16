import numpy as np
from mpi4py import MPI
import pytest

from veloxchem.veloxchemlib import mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfunrestdriver import ScfUnrestrictedDriver
from veloxchem.cppsolverunrest import ComplexResponseUnrestrictedSolver
from veloxchem.resultsio import read_results
from veloxchem.errorhandler import VeloxChemError


@pytest.mark.solvers
class TestCppUnrestricted:

    @staticmethod
    def _bcast_path_string(comm, path_string):

        if comm.Get_rank() != mpi_master():
            path_string = None

        return comm.bcast(path_string, root=mpi_master())

    def run_restarted_scf(self, molecule, basis, filename):

        scf_drv = ScfUnrestrictedDriver()
        scf_drv.ostream.mute()
        scf_drv.filename = filename
        scf_drv.compute(molecule, basis)
        scf_drv.restart = True

        return scf_drv.compute(molecule, basis)

    def _configure_cpp(self, driver, filename, frequencies):

        driver.ostream.mute()
        driver.filename = filename
        driver.conv_thresh = 1.0e-5
        driver.max_iter = 120
        driver.a_components = 'xz'
        driver.b_components = 'yz'
        driver.frequencies = frequencies
        driver.damping = 0.02
        driver.non_equilibrium_solv = False
        driver.ri_auxiliary_basis = 'def2-universal-jfit'

    def run_cpp(self,
                xcfun_label,
                cpp_property,
                ref_x_data,
                ref_y_data,
                tol,
                use_subcomms=False,
                max_subspace_dim=None):

        xyz_string = """6
        xyz
        O -3.42904  1.55532  0.01546
        C -1.99249  1.74379  0.02665
        H -1.74709  2.74160  0.44749
        H -1.59636  1.67836 -1.00861
        H -1.51398  0.95881  0.64937
        H -3.84726  2.33620 -0.34927
        """
        mol = Molecule.read_xyz_string(xyz_string)

        mol.set_multiplicity(3)

        basis_label = 'def2-svp'
        bas = MolecularBasis.read(mol, basis_label, ostream=None)

        scf_drv = ScfUnrestrictedDriver()
        scf_drv.ostream.mute()
        scf_drv.xcfun = xcfun_label
        scf_drv.acc_type = 'l2_c2diis'
        scf_results = scf_drv.compute(mol, bas)

        lr_drv = ComplexResponseUnrestrictedSolver()
        lr_drv.ostream.mute()
        lr_drv.property = cpp_property
        lr_drv.frequencies = list(ref_x_data)
        lr_drv.use_subcomms = use_subcomms
        lr_drv.max_subspace_dim = max_subspace_dim
        lr_results = lr_drv.compute(mol, bas, scf_results)

        if lr_drv.rank == mpi_master():
            lr_spec = lr_drv.get_spectrum(lr_results, 'au')
            assert np.max(
                np.abs(np.array(lr_spec['y_data']) -
                       np.array(ref_y_data))) < tol

    def test_hf_absorption(self):

        # vlxtag: UHF, Absorption, CPP

        xcfun_label = 'hf'
        cpp_property = 'absorption'
        ref_x_data = [0.33, 0.34, 0.35]
        ref_y_data = [0.05067935, 0.19924613, 0.39344688]

        self.run_cpp(xcfun_label,
                     cpp_property,
                     ref_x_data,
                     ref_y_data,
                     1.0e-6,
                     max_subspace_dim=120)

    def test_hf_ecd(self):

        # vlxtag: UHF, ECD, CPP

        xcfun_label = 'hf'
        cpp_property = 'ecd'
        ref_x_data = [0.33, 0.34, 0.35]
        ref_y_data = [-0.12250733, -0.62669734, -1.23612738]

        self.run_cpp(xcfun_label,
                     cpp_property,
                     ref_x_data,
                     ref_y_data,
                     1.0e-6,
                     max_subspace_dim=1000)

    def run_cpp_with_ecp(self, xcfun_label, cpp_property, ref_x_data,
                         ref_y_data, tol):

        xyz_string = """3
        xyz
        Au 0 0 0
        H  0 0 1.55
        H  0 1.53 0
        """
        mol = Molecule.read_xyz_string(xyz_string)
        mol.set_multiplicity(2)

        basis_label = 'def2-svp'
        bas = MolecularBasis.read(mol, basis_label, ostream=None)

        scf_drv = ScfUnrestrictedDriver()
        scf_drv.ostream.mute()
        scf_drv.xcfun = xcfun_label
        scf_results = scf_drv.compute(mol, bas)

        lr_drv = ComplexResponseUnrestrictedSolver()
        lr_drv.ostream.mute()
        lr_drv.property = cpp_property
        lr_drv.frequencies = list(ref_x_data)
        lr_results = lr_drv.compute(mol, bas, scf_results)

        if lr_drv.rank == mpi_master():
            lr_spec = lr_drv.get_spectrum(lr_results, 'au')
            assert np.max(
                np.abs(np.array(lr_spec['y_data']) -
                       np.array(ref_y_data))) < tol

    def test_hf_absorption_with_ecp(self):

        xcfun_label = 'hf'
        cpp_property = 'absorption'
        ref_x_data = [0.10, 0.11, 0.12]
        ref_y_data = [0.00666692, 0.00907027, 0.01351851]

        self.run_cpp_with_ecp(xcfun_label, cpp_property, ref_x_data, ref_y_data,
                              1.0e-6)

    def test_cpp_prop_dens_absorption(self, tmp_path):

        xyz_string = """3
        xyz
        O   -0.1858140  -1.1749469   0.7662596
        H   -0.1285513  -0.8984365   1.6808606
        H   -0.0582782  -0.3702550   0.2638279
        """
        mol = Molecule.read_xyz_string(xyz_string)
        mol.set_charge(1)
        mol.set_multiplicity(2)

        bas = MolecularBasis.read(mol, 'sto-3g', verbose=False)

        filename = self._bcast_path_string(
            MPI.COMM_WORLD, str(tmp_path / 'cpp_unrest_prop_dens_absorption'))

        scf_drv = ScfUnrestrictedDriver()
        scf_drv.ostream.mute()
        scf_drv.filename = filename
        scf_results = scf_drv.compute(mol, bas)

        cpp_drv = ComplexResponseUnrestrictedSolver()
        cpp_drv.frequencies = [0.40]
        cpp_drv.property = 'absorption'
        cpp_drv.ostream.mute()
        cpp_results_orig = cpp_drv.compute(mol, bas, scf_results)

        if scf_drv.rank == mpi_master():
            cpp_results_read = read_results(filename + '.h5', 'rsp')
        else:
            cpp_results_read = {}

        for cpp_results in [cpp_results_orig, cpp_results_read]:

            raw_density_dict = cpp_drv.get_cpp_property_densities(
                mol, bas, scf_results, cpp_results, 0.40, normalize_densities=False)
            density_dict = cpp_drv.get_cpp_property_densities(
                mol, bas, scf_results, cpp_results, 0.40, normalize_densities=True)

            if scf_drv.rank == mpi_master():
                overlap = scf_results['S']

                raw_detachment = raw_density_dict['property_density_detachment']
                raw_attachment = raw_density_dict['property_density_attachment']

                raw_detachment_int = -np.sum(raw_detachment * overlap)
                raw_attachment_int = np.sum(raw_attachment * overlap)

                assert raw_detachment_int > 0.0
                assert raw_detachment_int == pytest.approx(raw_attachment_int,
                                                           abs=1.0e-8)

                detachment = density_dict['property_density_detachment']
                attachment = density_dict['property_density_attachment']

                assert -np.sum(detachment * overlap) == pytest.approx(1.0,
                                                                      abs=1.0e-8)
                assert np.sum(attachment * overlap) == pytest.approx(1.0,
                                                                     abs=1.0e-8)

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason='skip pytest.raises for multiple MPI processes')
    def test_compute_rejects_nonlinear_rhs(self):

        xyz_string = """3
        xyz
        O   -0.1858140  -1.1749469   0.7662596
        H   -0.1285513  -0.8984365   1.6808606
        H   -0.0582782  -0.3702550   0.2638279
        """
        mol = Molecule.read_xyz_string(xyz_string)
        mol.set_charge(1)
        mol.set_multiplicity(2)

        bas = MolecularBasis.read(mol, 'sto-3g', ostream=None)

        scf_drv = ScfUnrestrictedDriver()
        scf_drv.ostream.mute()
        scf_results = scf_drv.compute(mol, bas)

        lr_drv = ComplexResponseUnrestrictedSolver()
        lr_drv.ostream.mute()
        lr_drv.frequencies = [0.10]

        with pytest.raises(VeloxChemError,
                           match='not implemented for nonlinear'):
            lr_drv.compute(mol, bas, scf_results, v_grad={('x', 0.10): None})

    def test_restart_adds_initial_guesses_for_more_frequencies(self, tmp_path):

        xyz_string = """3
        xyz
        O   -0.1858140  -1.1749469   0.7662596
        H   -0.1285513  -0.8984365   1.6808606
        H   -0.0582782  -0.3702550   0.2638279
        """
        mol = Molecule.read_xyz_string(xyz_string)
        mol.set_charge(1)
        mol.set_multiplicity(2)

        bas = MolecularBasis.read(mol, 'sto-3g', ostream=None)

        filename = self._bcast_path_string(
            MPI.COMM_WORLD, str(tmp_path / 'cpp_unrest_restart_extra_freqs'))
        scf_results = self.run_restarted_scf(mol, bas, filename)

        reference_drv = ComplexResponseUnrestrictedSolver()
        self._configure_cpp(reference_drv, filename, (0.10,))

        reference_scf_results = dict(scf_results)
        reference_scf_results['filename'] = filename
        reference_drv.compute(mol, bas, reference_scf_results)

        restarted_drv = ComplexResponseUnrestrictedSolver()
        self._configure_cpp(restarted_drv, filename, (0.10, 0.15))
        restarted_drv.checkpoint_file = reference_drv.checkpoint_file
        restarted_drv.restart = True

        restarted_scf_results = dict(scf_results)
        restarted_scf_results['filename'] = filename
        restarted_results = restarted_drv.compute(mol, bas,
                                                  restarted_scf_results)

        fresh_drv = ComplexResponseUnrestrictedSolver()
        fresh_filename = self._bcast_path_string(
            MPI.COMM_WORLD,
            str(tmp_path / 'cpp_unrest_restart_extra_freqs_fresh'))
        self._configure_cpp(fresh_drv, fresh_filename, (0.10, 0.15))

        fresh_scf_results = dict(scf_results)
        fresh_scf_results['filename'] = fresh_drv.filename
        fresh_results = fresh_drv.compute(mol, bas, fresh_scf_results)

        assert restarted_drv.restart is True

        if restarted_drv.rank == mpi_master():
            assert restarted_results['response_functions'].keys(
            ) == fresh_results['response_functions'].keys()
            for key, value in fresh_results['response_functions'].items():
                assert restarted_results['response_functions'][
                    key] == pytest.approx(value, abs=1.0e-8)
