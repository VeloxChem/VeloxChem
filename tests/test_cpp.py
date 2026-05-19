import numpy as np
from mpi4py import MPI
from pathlib import Path
import pytest
import h5py

from veloxchem.veloxchemlib import mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.cppsolver import ComplexResponseSolver


@pytest.mark.solvers
class TestCPP:

    @staticmethod
    def _bcast_path_string(comm, path_string):

        if comm.Get_rank() != mpi_master():
            path_string = None

        return comm.bcast(path_string, root=mpi_master())

    def run_restarted_scf(self, molecule, basis, filename):

        scf_drv = ScfRestrictedDriver()
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

        basis_label = 'def2-svp'
        bas = MolecularBasis.read(mol, basis_label, ostream=None)

        scf_drv = ScfRestrictedDriver()
        scf_drv.ostream.mute()
        scf_drv.xcfun = xcfun_label
        scf_drv.acc_type = 'l2_c2diis'
        scf_results = scf_drv.compute(mol, bas)

        lr_drv = ComplexResponseSolver()
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

    def run_cpp_plot(self,
                     xcfun_label,
                     cpp_property,
                     frequencies,
                     monkeypatch,
                     x_unit='nm',
                     plot_scatter=True):

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

        bas = MolecularBasis.read(mol, '6-31g', ostream=None)

        scf_drv = ScfRestrictedDriver()
        scf_drv.ostream.mute()
        scf_drv.xcfun = xcfun_label
        scf_results = scf_drv.compute(mol, bas)

        cpp_drv = ComplexResponseSolver()
        cpp_drv.ostream.mute()
        cpp_drv.property = cpp_property
        cpp_drv.frequencies = list(frequencies)
        cpp_results = cpp_drv.compute(mol, bas, scf_results)

        if cpp_drv.rank == mpi_master():
            plt = pytest.importorskip('matplotlib.pyplot')
            pytest.importorskip('scipy.interpolate')

            monkeypatch.setattr(plt, 'show', lambda: None)

            try:
                cpp_drv.plot(cpp_results,
                             x_unit=x_unit,
                             plot_scatter=plot_scatter)
            finally:
                plt.close('all')

    def test_hf_absorption(self):

        # vlxtag: RHF, Absorption, CPP

        xcfun_label = 'hf'
        cpp_property = 'absorption'
        ref_x_data = [0.39, 0.40, 0.41]
        ref_y_data = [0.11359517, 0.34554765, 1.41940099]

        self.run_cpp(xcfun_label,
                     cpp_property,
                     ref_x_data,
                     ref_y_data,
                     1.0e-6,
                     max_subspace_dim=120)

    def test_hf_ecd(self):

        # vlxtag: RHF, ECD, CPP

        xcfun_label = 'hf'
        cpp_property = 'ecd'
        ref_x_data = [0.39, 0.40, 0.41]
        ref_y_data = [5.74701958, 35.38200618, -15.89867350]

        self.run_cpp(xcfun_label, cpp_property, ref_x_data, ref_y_data, 1.0e-6)

        self.run_cpp(xcfun_label,
                     cpp_property,
                     ref_x_data,
                     ref_y_data,
                     1.0e-6,
                     use_subcomms=True)

    def test_b3lyp_absorption(self):

        # vlxtag: RKS, Absorption, CPP

        xcfun_label = 'b3lyp'
        cpp_property = 'absorption'
        ref_x_data = [0.39, 0.40, 0.41]
        ref_y_data = [0.20856963, 0.21772352, 0.71434376]

        self.run_cpp(xcfun_label, cpp_property, ref_x_data, ref_y_data, 1.0e-4)

    def test_b3lyp_ecd(self):

        # vlxtag: RKS, ECD, CPP

        xcfun_label = 'b3lyp'
        cpp_property = 'ecd'
        ref_x_data = [0.39, 0.40, 0.41]
        ref_y_data = [0.48472627, -1.18394788, -9.81377126]

        self.run_cpp(xcfun_label, cpp_property, ref_x_data, ref_y_data, 1.0e-4)

        self.run_cpp(xcfun_label,
                     cpp_property,
                     ref_x_data,
                     ref_y_data,
                     1.0e-4,
                     use_subcomms=True)

    def test_plot_absorption(self, monkeypatch):

        self.run_cpp_plot('b3lyp', 'absorption', [0.38, 0.39, 0.40, 0.41],
                          monkeypatch)

    def test_plot_ecd(self, monkeypatch):

        self.run_cpp_plot('hf',
                          'ecd', [0.38, 0.39, 0.40, 0.41],
                          monkeypatch,
                          x_unit='ev',
                          plot_scatter=False)

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason='skip pytest.raises for multiple MPI processes')
    def test_compute_rejects_openshell_molecule(self):

        xyz_string = """3
        xyz
        O   -0.1858140  -1.1749469   0.7662596
        H   -0.1285513  -0.8984365   1.6808606
        H   -0.0582782  -0.3702550   0.2638279
        """
        mol = Molecule.read_xyz_string(xyz_string)
        mol.set_multiplicity(3)

        bas = MolecularBasis.read(mol, 'def2-svp', ostream=None)

        # empty scf results just for testing
        scf_results = {}

        lr_drv = ComplexResponseSolver()
        lr_drv.ostream.mute()

        with pytest.raises(
                AssertionError,
                match="ComplexResponseSolver: not implemented for unrestricted case"):
            lr_results_not_used = lr_drv.compute(mol, bas, scf_results)

    def run_cpp_with_ecp(self, xcfun_label, cpp_property, ref_x_data,
                         ref_y_data, tol):

        xyz_string = """3
        xyz
        Au 0 0 0
        H  0 0 1.55
        H  0 1.53 0
        """
        mol = Molecule.read_xyz_string(xyz_string)
        mol.set_charge(1)

        basis_label = 'def2-svp'
        bas = MolecularBasis.read(mol, basis_label, ostream=None)

        scf_drv = ScfRestrictedDriver()
        scf_drv.ostream.mute()
        scf_drv.xcfun = xcfun_label
        scf_results = scf_drv.compute(mol, bas)

        lr_drv = ComplexResponseSolver()
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
        ref_x_data = [0.20, 0.21]
        ref_y_data = [0.03486827, 0.04576025]

        self.run_cpp_with_ecp(xcfun_label, cpp_property, ref_x_data, ref_y_data,
                              1.0e-6)

    def test_hf_absorption_nonlinear_rhs(self):

        xyz_string = """3
        xyz
        O   -0.1858140  -1.1749469   0.7662596
        H   -0.1285513  -0.8984365   1.6808606
        H   -0.0582782  -0.3702550   0.2638279
        """
        mol = Molecule.read_xyz_string(xyz_string)
        bas = MolecularBasis.read(mol, 'sto-3g', ostream=None)

        scf_drv = ScfRestrictedDriver()
        scf_drv.ostream.mute()
        scf_results = scf_drv.compute(mol, bas)

        lr_drv = ComplexResponseSolver()
        lr_drv.ostream.mute()
        lr_drv.frequencies = [0.10]
        lr_drv.damping = 0.02

        b_grad = lr_drv.get_complex_prop_grad(lr_drv.b_operator,
                                              lr_drv.b_components, mol, bas,
                                              scf_results)
        v_grad = {
            (op, 0.10): grad for op, grad in zip(lr_drv.b_components, b_grad)
        }

        lr_results = lr_drv.compute(mol, bas, scf_results, v_grad=v_grad)

        assert lr_drv.nonlinear is True
        if lr_drv.rank == mpi_master():
            assert set(lr_results.keys()) == {'focks', 'solutions'}
            assert set(lr_results['solutions'].keys()) == set(v_grad.keys())
            assert set(lr_results['focks'].keys()) == set(v_grad.keys())

    def test_restart_adds_initial_guesses_for_more_frequencies(self, tmp_path):

        xyz_string = """3
        xyz
        O   -0.1858140  -1.1749469   0.7662596
        H   -0.1285513  -0.8984365   1.6808606
        H   -0.0582782  -0.3702550   0.2638279
        """
        mol = Molecule.read_xyz_string(xyz_string)
        bas = MolecularBasis.read(mol, 'sto-3g', ostream=None)

        filename = self._bcast_path_string(
            MPI.COMM_WORLD, str(tmp_path / 'cpp_restart_extra_freqs'))
        scf_results = self.run_restarted_scf(mol, bas, filename)

        reference_drv = ComplexResponseSolver()
        self._configure_cpp(reference_drv, filename, (0.10,))

        reference_scf_results = dict(scf_results)
        reference_scf_results['filename'] = filename
        reference_drv.compute(mol, bas, reference_scf_results)

        restarted_drv = ComplexResponseSolver()
        self._configure_cpp(restarted_drv, filename, (0.10, 0.15))
        restarted_drv.checkpoint_file = reference_drv.checkpoint_file
        restarted_drv.restart = True

        restarted_scf_results = dict(scf_results)
        restarted_scf_results['filename'] = filename
        restarted_results = restarted_drv.compute(mol, bas,
                                                  restarted_scf_results)

        fresh_drv = ComplexResponseSolver()
        fresh_filename = self._bcast_path_string(
            MPI.COMM_WORLD, str(tmp_path / 'cpp_restart_extra_freqs_fresh'))
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

    def test_restart_allows_missing_frequency_dataset(self, tmp_path):

        xyz_string = """3
        xyz
        O   -0.1858140  -1.1749469   0.7662596
        H   -0.1285513  -0.8984365   1.6808606
        H   -0.0582782  -0.3702550   0.2638279
        """
        mol = Molecule.read_xyz_string(xyz_string)
        bas = MolecularBasis.read(mol, 'sto-3g', ostream=None)

        filename = self._bcast_path_string(
            MPI.COMM_WORLD, str(tmp_path / 'cpp_restart_missing_freqs'))
        scf_results = self.run_restarted_scf(mol, bas, filename)

        reference_drv = ComplexResponseSolver()
        self._configure_cpp(reference_drv, filename, (0.10,))

        reference_scf_results = dict(scf_results)
        reference_scf_results['filename'] = filename
        reference_results = reference_drv.compute(mol, bas,
                                                  reference_scf_results)
        checkpoint_file = self._bcast_path_string(MPI.COMM_WORLD,
                                                  reference_drv.checkpoint_file)

        if reference_drv.rank == mpi_master():
            assert Path(checkpoint_file).is_file()
            with h5py.File(checkpoint_file, 'a') as h5f:
                del h5f['frequencies']

        reference_drv.comm.barrier()

        restarted_drv = ComplexResponseSolver()
        self._configure_cpp(restarted_drv, filename, (0.10,))
        restarted_drv.checkpoint_file = checkpoint_file
        restarted_drv.restart = True

        restarted_scf_results = dict(scf_results)
        restarted_scf_results['filename'] = filename
        restarted_results = restarted_drv.compute(mol, bas,
                                                  restarted_scf_results)

        assert restarted_drv.restart is True

        if restarted_drv.rank == mpi_master():
            assert restarted_results['response_functions'].keys(
            ) == reference_results['response_functions'].keys()
            for key, value in reference_results['response_functions'].items():
                assert restarted_results['response_functions'][
                    key] == pytest.approx(value, abs=1.0e-8)
