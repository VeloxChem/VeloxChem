import numpy as np
import pytest
from mpi4py import MPI

from veloxchem.veloxchemlib import mpi_master
from veloxchem.veloxchemlib import SimdRIJKResponseDriver
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.scfunrestdriver import ScfUnrestrictedDriver
from veloxchem.lreigensolver import LinearResponseEigenSolver
from veloxchem.lreigensolverunrest import LinearResponseUnrestrictedEigenSolver
from veloxchem.tdaeigensolver import TdaEigenSolver
from veloxchem.quadraticresponsedriver import QuadraticResponseDriver
from veloxchem.thgdriver import ThgDriver
from veloxchem.cppsolver import ComplexResponseSolver
from veloxchem.profiler import Profiler
from veloxchem.nonlinearsolver import NonlinearSolver
from veloxchem.outputstream import OutputStream
from veloxchem import rijkresponse
from veloxchem import rijresponse


@pytest.mark.solvers
class TestResponseRiMemory:
    """The response's own memory check for RI: builds cut into chunks of
    densities, and the fallback to the four-center integrals where the RI
    driver does not fit."""

    @staticmethod
    def get_water():

        molecule = Molecule.read_molecule_string(
            """
            O  0.0           0.0  0.0
            H   .7586020000  0.0  -.5042840000
            H   .7586020000  0.0   .5042840000
            """,
            units='angstrom')
        basis = MolecularBasis.read(molecule, 'def2-svp', ostream=None)

        return molecule, basis

    @staticmethod
    def get_ch3():

        molecule = Molecule.read_molecule_string(
            """
            C   0.0000000   0.0000000   0.0000000
            H   1.0790000   0.0000000   0.0000000
            H  -0.5395000   0.9344500   0.0000000
            H  -0.5395000  -0.9344500   0.0000000
            """,
            units='angstrom')
        molecule.set_multiplicity(2)
        basis = MolecularBasis.read(molecule, 'def2-svp', ostream=None)

        return molecule, basis

    @staticmethod
    def run_scf(molecule, basis, xcfun, unrestricted=False):

        driver_type = (ScfUnrestrictedDriver
                       if unrestricted else ScfRestrictedDriver)
        scf_drv = driver_type()
        scf_drv.ostream.mute()
        scf_drv.xcfun = xcfun
        scf_drv.conv_thresh = 1.0e-8
        scf_results = scf_drv.compute(molecule, basis)
        assert not scf_drv._ri_fallback

        return scf_results

    @staticmethod
    def four_center(scf_results):
        """The SCF results with RI turned off, so that a response formed on
        them takes the four-center integrals on the same orbitals."""

        results = dict(scf_results)
        results['ri'] = False
        results['ri_jk'] = False
        results['ri_coulomb'] = False
        results['ri_fallback'] = False

        return results

    @staticmethod
    def run_eigensolver(solver_type, molecule, basis, scf_results,
                        budget=None):

        solver = solver_type()
        solver.ostream.mute()
        solver.nstates = 4
        solver.conv_thresh = 1.0e-6
        solver.ri_memory_budget = budget
        results = solver.compute(molecule, basis, scf_results)

        return solver, results

    @pytest.mark.parametrize('xcfun, solver_type', [
        ('b3lyp', LinearResponseEigenSolver),
        ('cam-b3lyp', LinearResponseEigenSolver),
        ('b3lyp', TdaEigenSolver),
    ])
    @pytest.mark.parametrize('count', [1, 2])
    def test_chunked_builds_match(self, monkeypatch, xcfun, solver_type,
                                  count):
        """A build cut into chunks of densities gives the matrices of the
        whole batch at once."""

        molecule, basis = self.get_water()
        scf_results = self.run_scf(molecule, basis, xcfun)

        whole_drv, whole = self.run_eigensolver(solver_type, molecule, basis,
                                                scf_results)
        assert whole_drv._uses_ri_jk()

        monkeypatch.setattr(rijkresponse, '_max_densities',
                            lambda *args: count)

        _, chunked = self.run_eigensolver(solver_type, molecule, basis,
                                          scf_results)

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            assert np.max(
                np.abs(whole['eigenvalues'] - chunked['eigenvalues'])) < 1.0e-10

    @pytest.mark.parametrize('count', [1, 2])
    def test_chunked_unrestricted_builds_match(self, monkeypatch, count):

        molecule, basis = self.get_ch3()
        scf_results = self.run_scf(molecule, basis, 'b3lyp', unrestricted=True)

        whole_drv, whole = self.run_eigensolver(
            LinearResponseUnrestrictedEigenSolver, molecule, basis, scf_results)
        assert whole_drv._uses_ri_jk()

        monkeypatch.setattr(rijkresponse, '_max_densities',
                            lambda *args: count)

        _, chunked = self.run_eigensolver(
            LinearResponseUnrestrictedEigenSolver, molecule, basis, scf_results)

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            assert np.max(
                np.abs(whole['eigenvalues'] - chunked['eigenvalues'])) < 1.0e-10

    @pytest.mark.parametrize('unrestricted', [False, True])
    def test_chunked_ri_j_builds_match(self, monkeypatch, unrestricted):

        if unrestricted:
            molecule, basis = self.get_ch3()
            solver_type = LinearResponseUnrestrictedEigenSolver
        else:
            molecule, basis = self.get_water()
            solver_type = LinearResponseEigenSolver

        scf_results = self.run_scf(molecule, basis, 'blyp', unrestricted)

        whole_drv, whole = self.run_eigensolver(solver_type, molecule, basis,
                                                scf_results)
        assert whole_drv._uses_ri_coulomb()

        monkeypatch.setattr(rijresponse, '_max_densities', lambda *args: 1)

        _, chunked = self.run_eigensolver(solver_type, molecule, basis,
                                          scf_results)

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            assert np.max(
                np.abs(whole['eigenvalues'] - chunked['eigenvalues'])) < 1.0e-10

    def test_chunked_block_diagonal_builds_match(self, monkeypatch):
        """Quadratic response hands over densities of the block diagonal
        shape, which a build forms as two calls."""

        molecule, basis = self.get_water()
        scf_results = self.run_scf(molecule, basis, 'hf')

        def run():
            drv = QuadraticResponseDriver()
            drv.ostream.mute()
            drv.update_settings({
                'conv_thresh': 1.0e-6,
                'b_frequencies': [0.2],
                'c_frequencies': [0.2],
                'damping': 0.1,
                'a_component': 'z',
                'b_component': 'z',
                'c_component': 'x',
            })
            return drv, drv.compute(molecule, basis, scf_results)

        whole_drv, whole = run()
        assert whole_drv._uses_ri_jk()

        monkeypatch.setattr(rijkresponse, '_max_densities',
                            lambda *args: 1)

        _, chunked = run()

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            key = ('qrf', 0.2, 0.2)
            assert abs(whole[key] - chunked[key]) < 1.0e-8

    @pytest.mark.parametrize('xcfun, label', [
        ('hf', 'RI-JK'),
        ('b3lyp', 'RI-JK'),
        ('blyp', 'RI-J'),
    ])
    def test_response_falls_back_on_its_own(self, xcfun, label):
        """An RI SCF whose response does not fit falls back to the four-center
        integrals in response alone, and equals a four-center response on the
        same orbitals."""

        molecule, basis = self.get_water()
        scf_results = self.run_scf(molecule, basis, xcfun)

        fallen_drv, fallen = self.run_eigensolver(LinearResponseEigenSolver,
                                                  molecule,
                                                  basis,
                                                  scf_results,
                                                  budget=1.0e-5)

        assert not fallen_drv._ri_fallback

        if label == 'RI-JK':
            assert fallen_drv._ri_jk_fallback
            assert fallen_drv.ri_jk
            assert not fallen_drv._uses_ri_jk()
            assert fallen_drv._ri_jk_drv is None
        else:
            assert fallen_drv._ri_j_fallback
            assert fallen_drv.ri_coulomb
            assert not fallen_drv._uses_ri_coulomb()
            assert fallen_drv._ri_drv is None

        _, plain = self.run_eigensolver(LinearResponseEigenSolver, molecule,
                                        basis, self.four_center(scf_results))

        _, ri = self.run_eigensolver(LinearResponseEigenSolver, molecule,
                                     basis, scf_results)

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            assert np.max(
                np.abs(fallen['eigenvalues'] - plain['eigenvalues'])) < 1.0e-10
            # and differs from the RI response by the fitting error
            diff = np.max(np.abs(fallen['eigenvalues'] - ri['eigenvalues']))
            assert diff > 1.0e-8
            assert diff < 1.0e-3

    def test_fallback_is_handed_to_sub_solvers(self):
        """A driver which falls back hands that to the solvers it drives, which
        do not form the B vectors again."""

        molecule, basis = self.get_water()
        scf_results = self.run_scf(molecule, basis, 'hf')

        def run(results, budget=None):
            drv = QuadraticResponseDriver()
            drv.ostream.mute()
            drv.ri_memory_budget = budget
            drv.update_settings({
                'conv_thresh': 1.0e-6,
                'b_frequencies': [0.2],
                'c_frequencies': [0.2],
                'damping': 0.1,
                'a_component': 'z',
                'b_component': 'z',
                'c_component': 'x',
            })
            return drv, drv.compute(molecule, basis, results)

        fallen_drv, fallen = run(scf_results, budget=1.0e-5)

        assert fallen_drv._ri_jk_fallback
        assert not fallen_drv._uses_ri_jk()

        _, plain = run(self.four_center(scf_results))

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            key = ('qrf', 0.2, 0.2)
            assert abs(fallen[key] - plain[key]) < 1.0e-8

    @staticmethod
    def run_thg(molecule, basis, scf_results, xcfun):

        drv = ThgDriver(MPI.COMM_WORLD, OutputStream(None))
        drv.update_settings(
            {
                'frequencies': [0.05],
                'damping': 0.1,
                'conv_thresh': 1.0e-6,
            }, {'xcfun': xcfun})
        results = drv.compute(molecule, basis, scf_results)

        return drv, results

    def test_nonlinear_reserve_floor(self):
        """A nonlinear driver keeps back from its B vectors at least the
        smallest batch of its Fock builds: for a THG with a functional one
        group of 12 second-order and 6 third-order densities."""

        molecule, basis = self.get_water()
        scf_results = self.run_scf(molecule, basis, 'b3lyp')

        drv, _ = self.run_thg(molecule, basis, scf_results, 'b3lyp')

        assert not drv._ri_jk_fallback
        assert drv._ri_jk_reserve_floor > 0
        assert drv._ri_memory_reserve >= drv._ri_jk_reserve_floor

        n_ao = basis.get_dimensions_of_basis()

        # NOTE: of each of the 18 densities, at least the master's 4 matrices
        # of n_ao by n_ao and its share of the 6 first-order densities.
        assert drv._ri_jk_reserve_floor >= 18 * (4 + 6 / 18) * 8 * n_ao * n_ao

    def test_nonlinear_driver_falls_back_for_its_floor(self, monkeypatch):
        """A nonlinear driver whose smallest Fock batch does not fit beside the
        B vectors falls back to the four-center integrals before any work, and
        equals a four-center driver on the same orbitals."""

        molecule, basis = self.get_water()
        scf_results = self.run_scf(molecule, basis, 'b3lyp')

        # NOTE: a memory of each density no rank has, which sets the floor
        # alone; the batch sizes of a four-center build do not read it.
        monkeypatch.setattr(NonlinearSolver, '_ri_density_memory',
                            lambda *args, **kwargs: 10**15)

        warnings = []
        monkeypatch.setattr(OutputStream, 'print_warning',
                            lambda ostream, line: warnings.append(line))

        fallen_drv, fallen = self.run_thg(molecule, basis, scf_results,
                                          'b3lyp')

        assert fallen_drv._ri_jk_fallback
        assert not fallen_drv._uses_ri_jk()
        assert fallen_drv._ri_jk_drv is None

        # NOTE: the floor is the reserve here, and the warning gives the free
        # memory and the floor.
        assert fallen_drv._ri_memory_reserve == fallen_drv._ri_jk_reserve_floor
        assert any('for the smallest batch of the nonlinear Fock build' in line
                   and 'four-center integrals' in line for line in warnings)

        monkeypatch.undo()

        _, plain = self.run_thg(molecule, basis, self.four_center(scf_results),
                                'b3lyp')

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            assert fallen['THG'].keys() == plain['THG'].keys()
            for key, value in plain['THG'].items():
                assert abs(fallen['THG'][key] - value) < 1.0e-8 * max(
                    abs(value), 1.0)

    def test_response_memory_model(self):
        """The prediction grows with the densities of a build, and the count
        the budget allows is the largest whose prediction fits."""

        molecule, basis = self.get_water()
        aux_basis = MolecularBasis.read(molecule,
                                        'def2-universal-jkfit',
                                        ostream=None,
                                        auxiliary_basis=True)

        nao = basis.get_dimensions_of_basis()
        nocc = molecule.number_of_alpha_occupied_orbitals(basis)
        all_atoms = list(range(molecule.number_of_atoms()))

        drv = SimdRIJKResponseDriver(1.0e-12, 0, 1024**2)

        forming, build, nheld = drv.fixed_memory(molecule, basis, aux_basis,
                                                 1.0e-12, all_atoms, False)

        assert nheld == aux_basis.get_dimensions_of_basis()
        assert forming > 0
        assert build > 0

        shape = (nao, nocc, 2, 1, 1)

        # forming is the larger peak for a molecule this small, so it would
        # hide the build; it is checked on its own below
        peaks = [
            drv.peak_memory(0, build, *shape, ndens, nheld)
            for ndens in range(1, 9)
        ]
        assert all(b >= a for a, b in zip(peaks, peaks[1:]))
        assert peaks[-1] > peaks[0]

        assert drv.peak_memory(forming, build, *shape, 1, nheld) >= forming
        assert drv.max_densities(forming, forming, build, *shape, nheld) == 0

        # the block diagonal shape holds the first call's matrices
        assert (drv.chunk_memory(nao, nocc, 1, 2, 1, 4, nheld) >
                drv.chunk_memory(nao, nocc, 1, 1, 1, 4, nheld))

        # two spins hold twice the matrices of one
        assert (drv.chunk_memory(nao, nocc, 2, 1, 2, 4, nheld) >
                drv.chunk_memory(nao, nocc, 2, 1, 1, 4, nheld))

        budget = peaks[4]
        count = drv.max_densities(budget, 0, build, *shape, nheld)

        assert count >= 5
        assert drv.peak_memory(0, build, *shape, count, nheld) <= budget
        assert drv.peak_memory(0, build, *shape, count + 1, nheld) > budget

        assert drv.max_densities(peaks[0] - 1, 0, build, *shape, nheld) == 0

        # a share of no atoms holds no B vectors and no transformed factors
        _, empty_build, empty_nheld = drv.fixed_memory(molecule, basis,
                                                       aux_basis, 1.0e-12, [],
                                                       False)
        assert empty_nheld == 0
        assert empty_build < build

    @pytest.mark.parametrize('xcfun', ['b3lyp', 'blyp'])
    def test_memory_profiling(self, xcfun):

        molecule, basis = self.get_water()
        scf_results = self.run_scf(molecule, basis, xcfun)

        solver = LinearResponseEigenSolver()
        solver.ostream.mute()
        solver.nstates = 3
        solver.memory_profiling = True
        solver.compute(molecule, basis, scf_results)

        assert solver._ri_rss_at_setup is not None

        if xcfun == 'b3lyp':
            assert solver._ri_jk_predicted_peak > 0
        else:
            assert solver._ri_j_memory is not None

        # NOTE: recorded on the master rank only, where the batch size is
        # decided; at least one density of the 7 matrices of n_ao by n_ao.
        if solver.rank == mpi_master():
            nao = basis.get_dimensions_of_basis()
            assert solver._ri_batch_memory >= 7 * nao * nao * 8
        else:
            assert solver._ri_batch_memory == 0

    @pytest.mark.parametrize('xcfun', ['b3lyp', 'blyp'])
    def test_memory_profiling_of_a_nonlinear_driver(self, xcfun):

        molecule, basis = self.get_water()
        scf_results = self.run_scf(molecule, basis, xcfun)

        drv = QuadraticResponseDriver()
        drv.ostream.mute()
        drv.update_settings({
            'conv_thresh': 1.0e-6,
            'b_frequencies': [0.2],
            'c_frequencies': [0.2],
            'damping': 0.1,
            'a_component': 'z',
            'b_component': 'z',
            'c_component': 'x',
            'memory_profiling': 'yes',
        })
        drv.compute(molecule, basis, scf_results)

        assert drv._ri_rss_at_setup is not None

        # NOTE: recorded by the driver's own Fock builds, on the master rank
        # only; at least one density of the 3 matrices of n_ao by n_ao.
        if drv.rank == mpi_master():
            nao = basis.get_dimensions_of_basis()
            assert drv._ri_batch_memory >= 3 * nao * nao * 8
        else:
            assert drv._ri_batch_memory == 0

    def test_ri_j_check_of_a_driver_leaves_out_its_solvers(self, monkeypatch):
        """A THG with RI-J sets up its RI-J driver before the linear solvers it
        drives, which set up their own and reset the peak. Its memory check
        leaves out what the solvers reached: 400 MB held and freed inside each
        solver's compute is not read as growth of the driver."""

        if not Profiler.reset_peak_memory():
            pytest.skip('the peak memory cannot be reset on this system')

        molecule, basis = self.get_water()
        scf_results = self.run_scf(molecule, basis, 'blyp')

        held = 400 * 1024**2
        compute = ComplexResponseSolver.compute

        def compute_and_hold(solver, *args, **kwargs):
            results = compute(solver, *args, **kwargs)
            block = np.ones(held // 8)
            del block
            return results

        monkeypatch.setattr(ComplexResponseSolver, 'compute', compute_and_hold)

        drv = ThgDriver(MPI.COMM_WORLD, OutputStream(None))
        drv.update_settings(
            {
                'frequencies': [0.05],
                'damping': 0.1,
                'conv_thresh': 1.0e-6,
                'memory_profiling': 'yes',
            }, {'xcfun': 'blyp'})
        drv.compute(molecule, basis, scf_results)

        assert drv._uses_ri_coulomb()

        grown = max(Profiler.get_peak_memory() - drv._ri_rss_at_setup,
                    getattr(drv, '_ri_peak_grown', 0))
        assert grown < held / 2
