from copy import deepcopy

from mpi4py import MPI
import numpy as np
import pytest

from veloxchem.veloxchemlib import mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.scfunrestdriver import ScfUnrestrictedDriver
from veloxchem.scfrestopendriver import ScfRestrictedOpenDriver
from veloxchem.scfgradientdriver import ScfGradientDriver
from veloxchem.tdaeigensolver import TdaEigenSolver
from veloxchem.tdaeigensolverunrest import TdaUnrestrictedEigenSolver
from veloxchem.lreigensolver import LinearResponseEigenSolver
from veloxchem.quadraticresponsedriver import QuadraticResponseDriver
from veloxchem.threepatransitiondriver import ThreePATransitionDriver
from veloxchem.excitedstatemomentdriver import ExcitedStateMomentDriver
from veloxchem.tddftgradientdriver import TddftGradientDriver
from veloxchem.errorhandler import VeloxChemError
from veloxchem import riresponse
from veloxchem.riresponse import lend_ri
from veloxchem.rijkresponse import initialize_ri_jk
from veloxchem.rijkresponse import ri_jk_fock_matrices


@pytest.mark.solvers
class TestResponseRiOwnership:
    """
    The RI-JK B vectors of a response calculation belong to the solver which
    formed them for that calculation: a solver computed again forms new ones,
    a driver lends its own to the solvers it drives for one call, and a deep
    copy holds none.
    """

    @staticmethod
    def get_water(geometry):

        # NOTE: the second geometry is the first with the hydrogens moved,
        # far enough that B vectors of one used for the other change the
        # excitation energies by more than 1e-3 hartree.
        if geometry == 1:
            xyz = """
            O   0.0   0.00   0.00
            H   0.0   0.76   0.59
            H   0.0  -0.76   0.59
            """
        else:
            xyz = """
            O   0.0   0.00   0.00
            H   0.0   0.86   0.49
            H   0.0  -0.80   0.62
            """

        molecule = Molecule.read_str(xyz, units='angstrom')
        basis = MolecularBasis.read(molecule, 'def2-svp', ostream=None)

        return molecule, basis

    @staticmethod
    def get_water_cation(geometry):

        molecule, _ = TestResponseRiOwnership.get_water(geometry)
        molecule.set_charge(1)
        molecule.set_multiplicity(2)
        basis = MolecularBasis.read(molecule, 'def2-svp', ostream=None)

        return molecule, basis

    @staticmethod
    def run_scf(scf_type, molecule, basis, xcfun):

        scf_drv = scf_type()
        scf_drv.ostream.mute()
        scf_drv.xcfun = xcfun
        scf_drv.conv_thresh = 1.0e-8
        scf_results = scf_drv.compute(molecule, basis)

        return scf_drv, scf_results

    @staticmethod
    def make_eigensolver(solver_type):

        solver = solver_type()
        solver.ostream.mute()
        solver.nstates = 3
        solver.conv_thresh = 1.0e-6

        return solver

    @pytest.mark.parametrize('scf_type, solver_type, get_molecule, xcfun', [
        (ScfRestrictedDriver, TdaEigenSolver, 'get_water', 'hf'),
        (ScfRestrictedDriver, TdaEigenSolver, 'get_water', 'b3lyp'),
        (ScfRestrictedDriver, LinearResponseEigenSolver, 'get_water', 'hf'),
        (ScfRestrictedDriver, LinearResponseEigenSolver, 'get_water',
         'cam-b3lyp'),
        (ScfUnrestrictedDriver, TdaUnrestrictedEigenSolver,
         'get_water_cation', 'b3lyp'),
    ])
    def test_solver_computed_again_at_another_geometry(self, scf_type,
                                                       solver_type,
                                                       get_molecule, xcfun):

        mol_1, bas_1 = getattr(self, get_molecule)(1)
        mol_2, bas_2 = getattr(self, get_molecule)(2)

        _, scf_results_1 = self.run_scf(scf_type, mol_1, bas_1, xcfun)
        _, scf_results_2 = self.run_scf(scf_type, mol_2, bas_2, xcfun)

        reused = self.make_eigensolver(solver_type)
        reused.compute(mol_1, bas_1, scf_results_1)
        reused_results = reused.compute(mol_2, bas_2, scf_results_2)

        assert reused._uses_ri_jk()

        fresh = self.make_eigensolver(solver_type)
        fresh_results = fresh.compute(mol_2, bas_2, scf_results_2)

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            assert np.max(
                np.abs(reused_results['eigenvalues'] -
                       fresh_results['eigenvalues'])) < 1.0e-10

    def test_driver_computed_again_at_another_geometry(self):

        mol_1, bas_1 = self.get_water(1)
        mol_2, bas_2 = self.get_water(2)

        _, scf_results_1 = self.run_scf(ScfRestrictedDriver, mol_1, bas_1,
                                        'b3lyp')
        _, scf_results_2 = self.run_scf(ScfRestrictedDriver, mol_2, bas_2,
                                        'b3lyp')

        rsp_settings = {
            'conv_thresh': 1.0e-6,
            'b_frequencies': [0.2],
            'c_frequencies': [0.2],
            'damping': 0.1,
            'a_component': 'z',
            'b_component': 'z',
            'c_component': 'z',
        }

        reused = QuadraticResponseDriver()
        reused.ostream.mute()
        reused.update_settings(rsp_settings)
        reused.compute(mol_1, bas_1, scf_results_1)
        reused_results = reused.compute(mol_2, bas_2, scf_results_2)

        assert reused._uses_ri_jk()
        assert not reused._ri_jk_lent

        fresh = QuadraticResponseDriver()
        fresh.ostream.mute()
        fresh.update_settings(rsp_settings)
        fresh_results = fresh.compute(mol_2, bas_2, scf_results_2)

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            key = ('qrf', 0.2, 0.2)
            assert abs(fresh_results[key]) > 1.0e-3
            assert abs(reused_results[key] - fresh_results[key]) < 1.0e-10

    def run_tda(self, molecule, basis, scf_results):

        solver = self.make_eigensolver(TdaEigenSolver)
        solver.compute(molecule, basis, scf_results)

        return solver

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason='skip pytest.raises for multiple MPI processes')
    def test_lent_b_vectors_are_returned(self):

        molecule, basis = self.get_water(1)
        _, scf_results = self.run_scf(ScfRestrictedDriver, molecule, basis,
                                      'b3lyp')

        driver = self.run_tda(molecule, basis, scf_results)
        solver = self.make_eigensolver(TdaEigenSolver)

        assert driver._ri_jk_drv is not None

        with lend_ri(driver, solver):
            assert driver._ri_jk_lent
            assert solver._ri_jk_borrowing

            with pytest.raises(VeloxChemError, match='are lent'):
                ri_jk_fock_matrices(driver, basis, None, 0.0)

            with pytest.raises(VeloxChemError, match='are lent'):
                initialize_ri_jk(driver, molecule, basis)

        assert not driver._ri_jk_lent
        assert not solver._ri_jk_borrowing
        assert solver._ri_jk_drv is None
        assert driver._ri_jk_drv is not None

        # NOTE: also when the call they were lent for raises.
        with pytest.raises(RuntimeError, match='raised inside'):
            with lend_ri(driver, solver):
                raise RuntimeError('raised inside')

        assert not driver._ri_jk_lent
        assert not solver._ri_jk_borrowing
        assert solver._ri_jk_drv is None

    def test_solver_lets_go_of_its_ri_j_driver(self):
        """A driver which uses RI-J lends nothing, and the solver it drives
        lets go of the RI-J driver it set up for itself on leaving the block,
        also when the call raises."""

        molecule, basis = self.get_water(1)
        _, scf_results = self.run_scf(ScfRestrictedDriver, molecule, basis,
                                      'blyp')

        driver = self.run_tda(molecule, basis, scf_results)

        assert driver._uses_ri_coulomb()
        assert driver._ri_drv is not None

        for raises in [False, True]:
            solver = self.make_eigensolver(TdaEigenSolver)

            try:
                with lend_ri(driver, solver):
                    solver.compute(molecule, basis, scf_results)
                    assert solver._ri_drv is not None
                    assert solver._ri_j_memory is not None
                    if raises:
                        raise RuntimeError('raised inside')
            except RuntimeError:
                assert raises

            assert solver._ri_drv is None
            assert solver._ri_j_memory is None
            assert driver._ri_drv is not None
            assert not driver._ri_jk_lent

    def test_ri_j_drivers_release_their_solvers(self, monkeypatch):
        """The 3PA transition and excited-state moment drivers with RI-J: each
        solver they drive lets go of the RI-J driver it set up when its block
        of lend_ri ends, and the results match the RI-J references."""

        xyzstr = """3
        xyz
        O   0.000000    0.000000    0.000000
        H   0.758602    0.000000   -0.504284
        H   0.758602    0.000000    0.504284
        """
        molecule = Molecule.read_xyz_string(xyzstr)
        basis = MolecularBasis.read(molecule, 'cc-pvdz', ostream=None)

        scf_drv = ScfRestrictedDriver()
        scf_drv.ostream.mute()
        scf_drv.xcfun = 'blyp'
        scf_drv.ri = True
        scf_drv.conv_thresh = 1.0e-8
        scf_results = scf_drv.compute(molecule, basis)

        if scf_drv.rank == mpi_master():
            assert scf_results['ri_coulomb'] and not scf_results['ri_jk']
            assert not scf_results['ri_fallback']

        # NOTE: the solvers a driver computes are not returned; the release
        # is recorded for each solver which holds an RI-J driver when lend_ri
        # lets go of it.
        released = []
        real_release = riresponse.release_ri_j

        def recorded_release(solver):
            if solver._ri_drv is not None:
                released.append(type(solver).__name__)
            real_release(solver)
            assert solver._ri_drv is None

        monkeypatch.setattr(riresponse, 'release_ri_j', recorded_release)

        three_pa_drv = ThreePATransitionDriver()
        three_pa_drv.ostream.mute()
        three_pa_drv.xcfun = 'blyp'
        three_pa_drv.nstates = 2
        three_pa_results = three_pa_drv.compute(molecule, basis, scf_results)

        assert three_pa_drv._uses_ri_coulomb()
        assert released == [
            'LinearResponseEigenSolver', 'ComplexResponseSolver',
            'ComplexResponseSolver'
        ]

        released.clear()

        esm_drv = ExcitedStateMomentDriver()
        esm_drv.ostream.mute()
        esm_drv.xcfun = 'blyp'
        esm_drv.state = 1
        esm_results = esm_drv.compute(molecule, basis, scf_results)

        assert esm_drv._uses_ri_coulomb()
        assert released == [
            'LinearResponseEigenSolver', 'ComplexResponseSolver'
        ]

        if scf_drv.rank == mpi_master():
            lin_3pa_str = np.array([
                val.real for val in three_pa_results['3pa_strengths']
                ['linear'].values()
            ])
            cir_3pa_str = np.array([
                val.real for val in three_pa_results['3pa_strengths']
                ['circular'].values()
            ])
            es_dipole = np.array(esm_results['excited_state_dipole_moment'])

            assert np.max(
                np.abs(lin_3pa_str - np.array([2865.274, 28484.970]))) < 0.1
            assert np.max(
                np.abs(cir_3pa_str - np.array([2472.439, 71212.424]))) < 0.1
            assert np.max(
                np.abs(es_dipole - np.array([-0.395035, 0.0, 0.0]))) < 1.0e-4

    def test_lending_hands_back_the_memory_of_the_builds(self):
        """The driver's memory check takes the largest RI build and batch of
        the solver it lent its B vectors to."""

        molecule, basis = self.get_water(1)
        _, scf_results = self.run_scf(ScfRestrictedDriver, molecule, basis,
                                      'b3lyp')

        driver = self.run_tda(molecule, basis, scf_results)
        driver._ri_jk_predicted_peak = 100
        driver._ri_batch_memory = 10

        for solver_peak, solver_batch, peak, batch in [(50, 5, 100, 10),
                                                       (300, 30, 300, 30)]:
            solver = self.make_eigensolver(TdaEigenSolver)

            with lend_ri(driver, solver):
                # NOTE: what the solver's builds record while it borrows.
                solver._ri_jk_predicted_peak = solver_peak
                solver._ri_batch_memory = solver_batch

            assert driver._ri_jk_predicted_peak == peak
            assert driver._ri_batch_memory == batch

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason='skip pytest.raises for multiple MPI processes')
    def test_reuse_without_b_vectors_is_refused(self):

        molecule, basis = self.get_water(1)

        solver = self.make_eigensolver(TdaEigenSolver)

        with pytest.raises(VeloxChemError, match='no RI-JK B vectors'):
            initialize_ri_jk(solver, molecule, basis, reuse=True)

    @pytest.mark.parametrize('solver_type', [
        TdaEigenSolver,
        LinearResponseEigenSolver,
    ])
    def test_eigensolver_deepcopy(self, solver_type):

        molecule, basis = self.get_water(1)
        _, scf_results = self.run_scf(ScfRestrictedDriver, molecule, basis,
                                      'b3lyp')

        solver = self.make_eigensolver(solver_type)
        results = solver.compute(molecule, basis, scf_results)

        assert solver._ri_jk_drv is not None

        solver_copy = deepcopy(solver)

        assert solver_copy._ri_jk_drv is None
        assert solver._ri_jk_drv is not None

        copy_results = solver_copy.compute(molecule, basis, scf_results)

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            assert np.max(
                np.abs(copy_results['eigenvalues'] -
                       results['eigenvalues'])) < 1.0e-10

    @pytest.mark.parametrize('scf_type, get_molecule', [
        (ScfRestrictedDriver, 'get_water'),
        (ScfUnrestrictedDriver, 'get_water_cation'),
    ])
    def test_scf_deepcopy_and_gradient(self, scf_type, get_molecule):

        molecule, basis = getattr(self, get_molecule)(1)
        scf_drv, _ = self.run_scf(scf_type, molecule, basis, 'b3lyp')

        assert scf_drv._ri_drv is not None

        scf_copy = deepcopy(scf_drv)

        assert scf_copy._ri_drv is None
        assert scf_drv._ri_drv is not None

        grad_drv = ScfGradientDriver(scf_drv)
        grad_drv.ostream.mute()
        grad_drv.compute(molecule, basis)

        scf_copy.compute(molecule, basis)

        copy_grad_drv = ScfGradientDriver(scf_copy)
        copy_grad_drv.ostream.mute()
        copy_grad_drv.compute(molecule, basis)

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            assert np.max(
                np.abs(copy_grad_drv.get_gradient() -
                       grad_drv.get_gradient())) < 1.0e-8

    def test_restricted_openshell_scf_deepcopy(self):

        molecule, basis = self.get_water_cation(1)
        scf_drv, scf_results = self.run_scf(ScfRestrictedOpenDriver, molecule,
                                            basis, 'b3lyp')

        scf_copy = deepcopy(scf_drv)

        assert scf_copy._ri_drv is None

        copy_results = scf_copy.compute(molecule, basis)

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            assert abs(copy_results['scf_energy'] -
                       scf_results['scf_energy']) < 1.0e-8

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason='skip pytest.raises for multiple MPI processes')
    def test_gradient_of_a_copy_which_has_not_run_the_scf(self):

        molecule, basis = self.get_water(1)
        scf_drv, scf_results = self.run_scf(ScfRestrictedDriver, molecule,
                                            basis, 'b3lyp')

        rsp_drv = self.run_tda(molecule, basis, scf_results)
        rsp_results = rsp_drv.compute(molecule, basis, scf_results)

        scf_copy = deepcopy(scf_drv)

        grad_drv = ScfGradientDriver(scf_copy)
        grad_drv.ostream.mute()

        with pytest.raises(VeloxChemError, match='holds no RI driver'):
            grad_drv.compute(molecule, basis)

        tddft_grad_drv = TddftGradientDriver(scf_copy)
        tddft_grad_drv.ostream.mute()

        with pytest.raises(VeloxChemError, match='holds no RI driver'):
            tddft_grad_drv.compute(molecule, basis, scf_copy, rsp_drv,
                                   rsp_results)
