from copy import deepcopy

from mpi4py import MPI
import numpy as np
import pytest

from veloxchem.veloxchemlib import mpi_master
from veloxchem.dispersionmodel import DispersionModel
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.outputstream import OutputStream
from veloxchem.scfgradientdriver import ScfGradientDriver
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.scfunrestdriver import ScfUnrestrictedDriver
from veloxchem.scfrestopendriver import ScfRestrictedOpenDriver
from veloxchem.errorhandler import VeloxChemError


@pytest.mark.solvers
class TestScfGradientDriverMiscellaneous:

    @staticmethod
    def get_h2_molecule_and_basis():

        xyz_string = """
        2
        h2
        H    0.000000    0.000000   -0.350000
        H    0.000000    0.000000    0.350000
        """
        molecule = Molecule.read_xyz_string(xyz_string)
        basis = MolecularBasis.read(molecule, 'sto-3g', ostream=None)
        return molecule, basis

    @staticmethod
    def get_h2o_molecule_and_basis():

        xyz_string = """
        3
        h2o
        O    0.000000    0.000000    0.000000
        H    0.000000   -0.757160    0.586260
        H    0.000000    0.757160    0.586260
        """
        molecule = Molecule.read_xyz_string(xyz_string)
        basis = MolecularBasis.read(molecule, 'sto-3g', ostream=None)
        return molecule, basis

    @staticmethod
    def get_ch3_molecule_and_basis():

        xyz_string = """
        4
        ch3
        C   -1.85334300   -0.63945100    1.29623300
        H   -2.40884500   -1.56570200    1.04276400
        H   -2.24160900   -0.22442700    2.24900500
        H   -1.98830700    0.10613700    0.48589200
        """
        molecule = Molecule.read_xyz_string(xyz_string)
        molecule.set_multiplicity(2)
        basis = MolecularBasis.read(molecule, 'sto-3g', ostream=None)
        return molecule, basis

    @staticmethod
    def run_restricted_scf(molecule, basis, xcfun='hf', *, ri):

        scf_drv = ScfRestrictedDriver()
        scf_drv.ri = ri
        scf_drv.ostream.mute()
        scf_drv.xcfun = xcfun
        scf_results = scf_drv.compute(molecule, basis)

        if scf_drv.rank == mpi_master():
            if ri:
                assert scf_results['ri_jk'] or scf_results['ri_coulomb']
                assert not scf_results['ri_fallback']

        return scf_drv, scf_results

    @staticmethod
    def run_unrestricted_scf(molecule, basis, xcfun='hf', filename=None, *, ri):

        scf_drv = ScfUnrestrictedDriver()
        scf_drv.ri = ri
        scf_drv.ostream.mute()
        scf_drv.xcfun = xcfun
        scf_drv.filename = filename
        scf_results = scf_drv.compute(molecule, basis)

        if scf_drv.rank == mpi_master():
            if ri:
                assert scf_results['ri_jk'] or scf_results['ri_coulomb']
                assert not scf_results['ri_fallback']

        return scf_drv, scf_results

    @staticmethod
    def run_restricted_openshell_scf(molecule, basis, xcfun='hf'):

        scf_drv = ScfRestrictedOpenDriver()
        scf_drv.ostream.mute()
        scf_drv.xcfun = xcfun
        scf_results = scf_drv.compute(molecule, basis)
        return scf_drv, scf_results

    @pytest.mark.parametrize('ri', [False, True])
    def test_compute_uses_internal_scf_results_for_numerical_gradient(
            self, monkeypatch, ri):

        molecule, basis = self.get_h2_molecule_and_basis()
        scf_drv, scf_results = self.run_restricted_scf(molecule, basis, ri=ri)
        grad_drv = ScfGradientDriver(scf_drv)
        grad_drv.numerical = True

        seen = {}

        def fake_compute_numerical(self, mol, ao_basis, scf_results):
            seen['molecule'] = mol
            seen['basis'] = ao_basis
            seen['scf_results'] = scf_results
            self.gradient = np.zeros((mol.number_of_atoms(), 3))

        monkeypatch.setattr(ScfGradientDriver, 'compute_numerical',
                            fake_compute_numerical)

        grad_drv.compute(molecule, basis)

        assert seen['molecule'] is molecule
        assert seen['basis'] is basis
        assert seen['scf_results'] is scf_drv.scf_results
        assert np.allclose(grad_drv.get_gradient(), 0.0)

    @pytest.mark.parametrize('ri', [False, True])
    def test_compute_runs_internal_scf_when_results_missing(self, ri):

        molecule, basis = self.get_h2_molecule_and_basis()

        scf_drv = ScfRestrictedDriver()
        scf_drv.ri = ri
        scf_drv.ostream.mute()
        scf_drv.xcfun = 'hf'

        grad_drv = ScfGradientDriver(scf_drv)
        grad_drv.compute(molecule, basis)

        assert scf_drv.scf_results is not None
        if scf_drv.rank == mpi_master():
            assert scf_drv.scf_results['scf_type'] == 'restricted'
            if ri:
                scf_results = scf_drv.scf_results
                assert scf_results['ri_jk'] or scf_results['ri_coulomb']
                assert not scf_results['ri_fallback']
        assert np.max(np.abs(grad_drv.get_gradient())) > 0.0

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason='skip pytest.raises for multiple MPI processes')
    def test_compute_rejects_restricted_openshell_marker(self):

        molecule, basis = self.get_h2_molecule_and_basis()
        scf_drv, scf_results = self.run_restricted_openshell_scf(
            molecule, basis)
        grad_drv = ScfGradientDriver(scf_drv)

        with pytest.raises(VeloxChemError,
                           match='Not implemented for restricted open-shell'):
            grad_drv.compute(molecule, basis)

    @staticmethod
    def run_ri_j_gradient(scf_type, molecule, basis, numerical):

        scf_drv = scf_type()
        scf_drv.ri = True
        scf_drv.ostream.mute()
        scf_drv.xcfun = 'blyp'
        scf_drv.grid_level = 6
        scf_drv.conv_thresh = 1.0e-10
        scf_drv.compute(molecule, basis)

        assert scf_drv._uses_ri_coulomb()

        grad_drv = ScfGradientDriver(scf_drv)
        grad_drv.ostream.mute()
        grad_drv.numerical = numerical
        grad_drv.compute(molecule, basis)

        return grad_drv.get_gradient()

    @staticmethod
    def run_hf_gradient(molecule, basis, ri, numerical):

        scf_drv = ScfRestrictedDriver()
        scf_drv.ri = ri
        scf_drv.ostream.mute()
        scf_drv.conv_thresh = 1.0e-10
        scf_drv.compute(molecule, basis)

        assert scf_drv._uses_ri_jk() == ri

        grad_drv = ScfGradientDriver(scf_drv)
        grad_drv.ostream.mute()
        grad_drv.numerical = numerical
        grad_drv.compute(molecule, basis)

        return grad_drv.get_gradient()

    # NOTE: def2-SVPD has two s and two p shells on hydrogen, so H2 has as
    # many p as s shells (4 each). The number of p shells is then also one of
    # the p orbital indices (4 to 7), which the four-center gradient takes
    # through the same local indices as the Fock build.
    @pytest.mark.parametrize('ri', [False, True])
    def test_hf_gradient_with_as_many_p_as_s_shells(self, ri):

        molecule = Molecule.read_xyz_string("""2
        h2
        H    0.000000    0.000000    0.000000
        H    0.000000    0.000000    0.800000
        """)
        basis = MolecularBasis.read(molecule, 'def2-svpd', ostream=None)

        analytical = self.run_hf_gradient(molecule, basis, ri, False)
        numerical = self.run_hf_gradient(molecule, basis, ri, True)

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            assert np.max(np.abs(analytical)) > 1.0e-2
            assert np.max(np.abs(analytical - numerical)) < 1.0e-6

    # NOTE: STO-3G is not a def2 basis set. The analytical RI-J gradient is
    # compared with the numerical gradient of the RI-J energy.
    @pytest.mark.parametrize('scf_type, label', [
        (ScfRestrictedDriver, 'h2o'),
        (ScfUnrestrictedDriver, 'ch3'),
    ])
    def test_ri_j_gradient_with_non_def2_basis(self, scf_type, label):

        molecule, basis = getattr(self, f'get_{label}_molecule_and_basis')()

        analytical = self.run_ri_j_gradient(scf_type, molecule, basis, False)
        numerical = self.run_ri_j_gradient(scf_type, molecule, basis, True)

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            assert np.max(np.abs(analytical - numerical)) < 1.0e-5

    def test_restricted_ri_jk_gradient(self):

        molecule, basis = self.get_h2o_molecule_and_basis()

        scf_drv = ScfRestrictedDriver()
        scf_drv.ri = True
        scf_drv.ostream.mute()
        scf_drv.xcfun = 'hf'
        scf_results = scf_drv.compute(molecule, basis)

        grad_drv = ScfGradientDriver(scf_drv)
        grad_drv.compute(molecule, basis, scf_results)

        assert np.max(np.abs(grad_drv.get_gradient())) > 0.0

    def test_unrestricted_ri_jk_gradient(self):

        molecule, basis = self.get_ch3_molecule_and_basis()

        scf_drv = ScfUnrestrictedDriver()
        scf_drv.ri = True
        scf_drv.ostream.mute()
        scf_drv.xcfun = 'hf'
        scf_results = scf_drv.compute(molecule, basis)

        grad_drv = ScfGradientDriver(scf_drv)
        grad_drv.compute(molecule, basis, scf_results)

        assert np.max(np.abs(grad_drv.get_gradient())) > 0.0

    def get_h2o_def2_molecule_and_basis(self):

        molecule, _ = self.get_h2o_molecule_and_basis()
        basis = MolecularBasis.read(molecule, 'def2-svp', ostream=None)
        return molecule, basis

    @pytest.mark.parametrize('xcfun', ['hf', 'blyp', 'cam-b3lyp'])
    def test_ri_gradient_falls_back_to_four_center(self, xcfun):
        """A gradient budget below the peak the RI gradient is predicted to reach
        has the two-electron gradient formed from the four-center integrals, of
        the same density. RI-JK, RI-J and the range separated RI-JK each take
        the fallback apart, so each is run."""

        molecule, basis = self.get_h2o_def2_molecule_and_basis()
        scf_drv, scf_results = self.run_restricted_scf(molecule, basis,
                                                       xcfun=xcfun, ri=True)

        # ten kilobytes, which the peak of any molecule exceeds
        cramped = ScfGradientDriver(scf_drv)
        cramped.ri_gradient_memory_budget = 1.0e-5
        cramped.compute(molecule, basis, scf_results)

        assert cramped._ri_fallback

        roomy = ScfGradientDriver(scf_drv)
        roomy.compute(molecule, basis, scf_results)

        assert not roomy._ri_fallback

        # NOTE: the same density, with the gradient told the SCF formed its
        # Fock matrices without RI, which is the four-center gradient of it.
        scf_drv.ri_jk = False
        scf_drv.ri_coulomb = False

        plain = ScfGradientDriver(scf_drv)
        plain.compute(molecule, basis, scf_results)

        diff = np.max(np.abs(cramped.get_gradient() - plain.get_gradient()))
        assert diff < 1.0e-10

        diff = np.max(np.abs(roomy.get_gradient() - plain.get_gradient()))
        assert diff > 1.0e-8
        assert diff < 1.0e-3

    @pytest.mark.parametrize('unrestricted, xcfun', [
        (False, 'hf'),
        (False, 'blyp'),
        (True, 'b3lyp'),
    ])
    def test_gradient_follows_an_scf_fallback(self, unrestricted, xcfun):
        """An SCF which fell back to the four-center integrals holds no RI
        driver, and its gradient is the four-center one."""

        if unrestricted:
            molecule, _ = self.get_ch3_molecule_and_basis()
        else:
            molecule, _ = self.get_h2o_molecule_and_basis()
        basis = MolecularBasis.read(molecule, 'def2-svp', ostream=None)

        def run(ri, budget=None):
            driver_type = (ScfUnrestrictedDriver
                           if unrestricted else ScfRestrictedDriver)
            scf_drv = driver_type()
            scf_drv.ostream.mute()
            scf_drv.ri = ri
            scf_drv.xcfun = xcfun
            scf_drv.conv_thresh = 1.0e-8
            scf_drv.ri_memory_budget = budget
            scf_results = scf_drv.compute(molecule, basis)
            grad_drv = ScfGradientDriver(scf_drv)
            grad_drv.compute(molecule, basis, scf_results)
            return scf_drv, grad_drv

        scf_drv, fallen = run(ri=True, budget=1.0e-5)

        assert scf_drv._ri_fallback
        assert not fallen._ri_fallback

        _, plain = run(ri=False)

        diff = np.max(np.abs(fallen.get_gradient() - plain.get_gradient()))
        assert diff < 1.0e-6

    def test_unrestricted_range_separated_ri_jk_gradient(self):
        """The open shell range separated RI-JK gradient on one rank names the
        auxiliary atoms, which have no default. It once named the atoms alone and
        was refused by the binding."""

        molecule, _ = self.get_ch3_molecule_and_basis()
        basis = MolecularBasis.read(molecule, 'def2-svp', ostream=None)

        scf_drv, scf_results = self.run_unrestricted_scf(molecule, basis,
                                                         xcfun='cam-b3lyp',
                                                         ri=True)

        ri_drv = ScfGradientDriver(scf_drv)
        ri_drv.compute(molecule, basis, scf_results)

        assert not ri_drv._ri_fallback

        scf_drv.ri_jk = False
        scf_drv.ri_coulomb = False

        plain = ScfGradientDriver(scf_drv)
        plain.compute(molecule, basis, scf_results)

        diff = np.max(np.abs(ri_drv.get_gradient() - plain.get_gradient()))
        assert diff > 1.0e-8
        assert diff < 1.0e-3

    @pytest.mark.parametrize('xcfun', ['hf', 'blyp'])
    def test_ri_gradient_memory_profiling(self, xcfun):
        """With memory profiling the gradient keeps what was predicted for the
        RI term and what the peak grew by over it."""

        molecule, basis = self.get_h2o_def2_molecule_and_basis()

        scf_drv = ScfRestrictedDriver()
        scf_drv.ri = True
        scf_drv.ostream.mute()
        scf_drv.xcfun = xcfun
        scf_drv.memory_profiling = True
        scf_results = scf_drv.compute(molecule, basis)

        grad_drv = ScfGradientDriver(scf_drv)
        grad_drv.compute(molecule, basis, scf_results)

        label, predicted, reserve, rss_before, grown = grad_drv._ri_memory_check

        assert label == ('RI-J' if xcfun == 'blyp' else 'RI-JK')
        assert predicted > 0
        assert reserve > 0
        assert rss_before > 0
        assert grown is not None

    def test_ri_gradient_peak_memory(self):
        """The prediction grows with what the gradient forms: a second set of
        fitted densities for a range separated functional, and a second spin."""

        from veloxchem.veloxchemlib import SimdRIJKGradientDriver
        from veloxchem.veloxchemlib import SimdRIJGradientDriver

        molecule, basis = self.get_h2o_def2_molecule_and_basis()
        aux_basis = MolecularBasis.read(molecule,
                                        'def2-universal-jkfit',
                                        ostream=None,
                                        auxiliary_basis=True)
        atoms = list(range(molecule.number_of_atoms()))

        drv = SimdRIJKGradientDriver(1.0e-12)

        def peak(norbitals, range_separated, distributed=False):
            return drv.peak_memory(molecule, basis, aux_basis, atoms,
                                   norbitals, range_separated, distributed,
                                   256 * 1024**2, True)

        closed = peak([5], False)

        assert closed > 0
        assert peak([5], True) > closed
        assert peak([5, 4], False) > closed
        assert peak([5], False, distributed=True) > 0

        # a rank which holds no auxiliary atoms and forms no two-center term
        # holds the squares of the transformation alone
        assert drv.peak_memory(molecule, basis, aux_basis, [], [5], False,
                               True, 256 * 1024**2, False) < closed

        rij = SimdRIJGradientDriver(1.0e-12)

        assert rij.peak_memory(molecule, basis, aux_basis, atoms, True) > 0

    @pytest.mark.parametrize('ri', [False, True])
    def test_restricted_gradient_timing_output_path(self, ri):

        molecule, basis = self.get_h2_molecule_and_basis()
        scf_drv, scf_results = self.run_restricted_scf(molecule, basis, ri=ri)
        scf_drv.timing = True

        grad_drv = ScfGradientDriver(scf_drv)
        grad_drv.compute(molecule, basis, scf_results)

        assert np.max(np.abs(grad_drv.get_gradient())) > 0.0

    @pytest.mark.parametrize('ri', [False, True])
    def test_compute_energy_and_restart_mode(self, ri):

        molecule, basis = self.get_h2_molecule_and_basis()
        scf_drv, scf_results = self.run_restricted_scf(molecule, basis, ri=ri)
        grad_drv = ScfGradientDriver(scf_drv)

        original_compute = scf_drv.compute
        restart_states = []

        def wrapped_compute(*args, **kwargs):
            restart_states.append(scf_drv.restart)
            return original_compute(*args, **kwargs)

        scf_drv.compute = wrapped_compute

        energy = grad_drv.compute_energy(molecule, basis)

        assert restart_states == [True]
        assert np.isclose(energy, scf_drv.get_scf_energy())

        grad_drv.numerical = True
        grad_drv.compute_energy(molecule, basis)

        assert restart_states == [True, False]

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason='skip pytest.raises for multiple MPI processes')
    @pytest.mark.parametrize('ri', [False, True])
    def test_compute_energy_raises_on_unconverged_scf(self, monkeypatch, ri):

        molecule, basis = self.get_h2_molecule_and_basis()
        scf_drv, scf_results = self.run_restricted_scf(molecule, basis, ri=ri)
        grad_drv = ScfGradientDriver(scf_drv)

        def fake_compute(*args, **kwargs):
            scf_drv._is_converged = False
            return {}

        monkeypatch.setattr(scf_drv, 'compute', fake_compute)

        with pytest.raises(VeloxChemError, match='SCF did not converge'):
            grad_drv.compute_energy(molecule, basis)

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason='skip pytest.raises for multiple MPI processes')
    @pytest.mark.parametrize('ri', [False, True])
    def test_unrestricted_cpcm_rejects_smd_in_gradient(self, ri):

        molecule, basis = self.get_ch3_molecule_and_basis()
        scf_drv, scf_results = self.run_unrestricted_scf(molecule,
                                                         basis,
                                                         'b3lyp',
                                                         ri=ri)
        scf_drv._cpcm = True
        scf_drv._smd = True

        grad_drv = ScfGradientDriver(scf_drv)

        with pytest.raises(VeloxChemError,
                           match='Cannot use SMD in gradient calculation'):
            grad_drv.compute_analytical_unrestricted(molecule, basis,
                                                     scf_results)

    @pytest.mark.parametrize('ri', [False, True])
    def test_unrestricted_cpcm_gradient_and_timing_output(self, tmp_path, ri):

        molstr = """
        N         -1.96309        1.59755       -0.01963
        H         -1.95876        2.61528        0.03109
        H         -2.48929        1.27814        0.79244
        H         -2.52930        1.35928       -0.83265
        """
        molecule = Molecule.read_molecule_string(molstr, units='bohr')
        molecule.set_charge(1)
        molecule.set_multiplicity(2)
        basis = MolecularBasis.read(molecule, 'def2-svp', ostream=None)

        output_file = tmp_path / 'unrestricted_cpcm_gradient.out'
        comm = MPI.COMM_WORLD
        output_file = comm.bcast(output_file, root=mpi_master())
        ostream = OutputStream.create_mpi_ostream(comm, str(output_file))

        scf_drv = ScfUnrestrictedDriver(comm, ostream)
        scf_drv.ri = ri
        scf_drv.xcfun = 'hf'
        scf_drv.solvation_model = 'cpcm'
        scf_drv.cpcm_grid_per_sphere = (110, 110)
        scf_drv.timing = True

        scf_results = scf_drv.compute(molecule, basis)

        if scf_drv.rank == mpi_master():
            if ri:
                assert scf_results['ri_jk'] or scf_results['ri_coulomb']
                assert not scf_results['ri_fallback']

        grad_drv = ScfGradientDriver(scf_drv)
        grad_drv.compute(molecule, basis, scf_results)

        ref_grad = np.array([
            [-2.50053996, 1.05721786, 0.11437744],
            [-0.09280072, -2.53390389, -0.12426393],
            [1.24618117, 0.84073729, -2.04587895],
            [1.34715952, 0.63594873, 2.05576544],
        ])

        assert np.max(np.abs(grad_drv.get_gradient() - ref_grad)) < 1.0e-4

        ostream.close()
        if scf_drv.rank == mpi_master():
            output_text = output_file.read_text()
            assert 'Gradient timing decomposition' in output_text
            assert 'CPCM_grad' in output_text

    @pytest.mark.skipif(not DispersionModel.is_available(),
                        reason='dftd4-python not available')
    @pytest.mark.parametrize('ri', [False, True])
    def test_unrestricted_dispersion_gradient_matches_added_d4_component(
            self, ri):

        molecule, basis = self.get_ch3_molecule_and_basis()

        base_scf_drv, base_results = self.run_unrestricted_scf(molecule,
                                                               basis,
                                                               'b3lyp',
                                                               ri=ri)
        base_grad_drv = ScfGradientDriver(base_scf_drv)
        base_grad_drv.compute(molecule, basis, base_results)
        base_gradient = base_grad_drv.get_gradient()

        disp_scf_drv, disp_results = self.run_unrestricted_scf(molecule,
                                                               basis,
                                                               'b3lyp',
                                                               ri=ri)
        disp_scf_drv.dispersion = True
        disp_grad_drv = ScfGradientDriver(disp_scf_drv)
        disp_grad_drv.compute(molecule, basis, disp_results)
        disp_gradient = disp_grad_drv.get_gradient()

        d4_model = DispersionModel()
        d4_model.compute(molecule, 'b3lyp')

        assert np.max(
            np.abs((disp_gradient - base_gradient) -
                   d4_model.get_gradient())) < 1.0e-8

    @pytest.mark.parametrize('ri', [False, True])
    def test_gradient_driver_deepcopy(self, ri):

        molecule, basis = self.get_ch3_molecule_and_basis()

        scf_drv, scf_results = self.run_unrestricted_scf(molecule,
                                                         basis,
                                                         'b3lyp',
                                                         ri=ri)
        grad_drv = ScfGradientDriver(scf_drv)
        grad_drv.compute(molecule, basis, scf_results)

        grad_drv_copy = deepcopy(grad_drv)

        assert grad_drv_copy.xcfun == grad_drv.xcfun
        assert np.allclose(grad_drv_copy.gradient, grad_drv.gradient)

    @pytest.mark.parametrize('ri', [False, True])
    def test_read_settings(self, tmp_path, ri):

        molecule, basis = self.get_ch3_molecule_and_basis()

        filename = str(tmp_path / "grad_import_settings")

        comm = MPI.COMM_WORLD
        filename = comm.bcast(filename, root=mpi_master())

        scf_drv, scf_results = self.run_unrestricted_scf(molecule,
                                                         basis,
                                                         'b3lyp',
                                                         filename=filename,
                                                         ri=ri)
        grad_drv = ScfGradientDriver(scf_drv)
        grad_drv.compute(molecule, basis, scf_results)

        checkpoint_file = f'{filename}_scf.h5'

        new_scf_drv = ScfUnrestrictedDriver()
        new_scf_drv.ri = ri
        new_grad_drv = ScfGradientDriver(new_scf_drv)
        new_grad_drv.read_settings(checkpoint_file)

        # these should be updated by read_settings
        assert (new_grad_drv.scf_driver.xcfun.get_func_label() ==
                scf_drv.xcfun.get_func_label())
        assert (new_grad_drv.xcfun.get_func_label() ==
                scf_drv.xcfun.get_func_label())

        # these should not be updated by read_settings
        assert new_grad_drv.scf_driver.filename is None
        assert new_grad_drv.scf_driver.checkpoint_file is None
        assert new_grad_drv.gradient is None
