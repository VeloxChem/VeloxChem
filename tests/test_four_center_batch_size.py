from mpi4py import MPI
import numpy as np
import platform
import pytest

from veloxchem.veloxchemlib import mpi_master
from veloxchem.batchsize import get_four_center_batch_size
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.scfunrestdriver import ScfUnrestrictedDriver
from veloxchem.lreigensolver import LinearResponseEigenSolver
from veloxchem.lreigensolverunrest import LinearResponseUnrestrictedEigenSolver
from veloxchem.cppsolver import ComplexResponseSolver
from veloxchem.hessianorbitalresponse import HessianOrbitalResponse
from veloxchem.quadraticresponsedriver import QuadraticResponseDriver
import veloxchem.batchsize as batchsize
import veloxchem.linearsolver as linearsolver
import veloxchem.cphfsolver as cphfsolver
import veloxchem.nonlinearsolver as nonlinearsolver


@pytest.mark.solvers
class TestFourCenterBatchSize:
    """Response Fock builds from the four-center integrals take their batch
    size from get_four_center_batch_size, and a build cut into many batches
    gives what one batch gives."""

    comm = MPI.COMM_WORLD

    nao = 100

    nthreads = 4

    def set_memory(self, monkeypatch, available):
        """Sets the memory get_four_center_batch_size reads as free on every
        rank.

        NOTE: it reads it through get_available_memory and the threads
        through get_number_of_threads, the names batchsize holds."""

        monkeypatch.setattr(batchsize, 'get_available_memory',
                            lambda: available)
        monkeypatch.setattr(batchsize, 'get_number_of_threads',
                            lambda: self.nthreads)

        here = platform.node()

        return available // self.comm.allgather(here).count(here)

    def per_density(self):
        """The memory of a density, 8 x nao^2 bytes for each thread and no
        fewer than _FOUR_CENTER_DENSITY_MATRICES_FLOOR such matrices."""

        return (max(batchsize._FOUR_CENTER_DENSITY_MATRICES_FLOOR,
                    self.nthreads) * self.nao**2 * 8)

    def eri_memory(self):
        """The memory of the ERI integration, held once by a build."""

        return (batchsize._FOUR_CENTER_ERI_MATRICES_PER_THREAD * self.nao**2 *
                8 * self.nthreads)

    def test_the_count_follows_the_room(self, monkeypatch):

        room = self.set_memory(monkeypatch, 8 * 1024**3)

        reserve = self.eri_memory()

        assert get_four_center_batch_size(
            None, 10**9, self.nao,
            self.comm) == (room - reserve) // self.per_density()

        call_memory = room // 4

        assert get_four_center_batch_size(
            None, 10**9, self.nao, self.comm,
            call_memory) == ((room - call_memory - reserve) //
                             self.per_density())

    def test_a_density_takes_the_floor_or_a_matrix_for_each_thread(
            self, monkeypatch):
        """With fewer threads than _FOUR_CENTER_DENSITY_MATRICES_FLOOR a
        density is charged that many matrices, with more one matrix for each
        thread; the ERI memory is three matrices for each thread either way."""

        floor = batchsize._FOUR_CENTER_DENSITY_MATRICES_FLOOR
        matrix = self.nao**2 * 8

        for threads, per_density in [(1, floor * matrix),
                                     (2 * floor, 2 * floor * matrix)]:
            self.nthreads = threads

            room = self.set_memory(monkeypatch, 8 * 1024**3)

            assert get_four_center_batch_size(
                None, 10**9, self.nao,
                self.comm) == (room - self.eri_memory()) // per_density

    def test_the_input_batch_size_and_the_items_bound_the_count(
            self, monkeypatch):

        room = self.set_memory(monkeypatch, 8 * 1024**3)

        assert get_four_center_batch_size(9, 10**9, self.nao, self.comm) == 9
        assert get_four_center_batch_size(None, 5, self.nao, self.comm) == 5
        assert get_four_center_batch_size(None, 0, self.nao, self.comm) == 0

        # NOTE: a call that leaves no room for a density still gives one.
        assert get_four_center_batch_size(None, 10**9, self.nao, self.comm,
                                          2 * room) == 1

    def test_no_room_still_gives_one_density(self, monkeypatch):

        # NOTE: a room below the ERI integration's memory, or none, gives one density; the
        # build is attempted, never stopped.
        for available in [0, 1, self.per_density()]:
            self.set_memory(monkeypatch, available)

            assert get_four_center_batch_size(None, 10**9, self.nao,
                                              self.comm) == 1
            assert get_four_center_batch_size(None, 10**9, self.nao,
                                              self.comm, 1024**4) == 1

    def test_the_rank_with_the_least_room_sets_the_count(self, monkeypatch):

        rank = self.comm.Get_rank()
        nodes = self.comm.Get_size()

        # NOTE: the last rank sees half the memory of the others.
        available = 8 * 1024**3 if rank < nodes - 1 else 4 * 1024**3
        self.set_memory(monkeypatch, available)

        here = platform.node()
        smallest = min(
            a // self.comm.allgather(here).count(h)
            for h, a in self.comm.allgather((here, available)))

        reserve = self.eri_memory()

        assert get_four_center_batch_size(
            None, 10**9, self.nao,
            self.comm) == (smallest - reserve) // self.per_density()

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
    def run_scf(molecule, basis, xcfun, ri=False, unrestricted=False):

        driver_type = (ScfUnrestrictedDriver
                       if unrestricted else ScfRestrictedDriver)
        scf_drv = driver_type()
        scf_drv.ostream.mute()
        scf_drv.ri = ri
        scf_drv.xcfun = xcfun
        scf_drv.conv_thresh = 1.0e-8
        scf_results = scf_drv.compute(molecule, basis)

        return scf_results

    @staticmethod
    def spy(monkeypatch, module, little_memory=False):
        """Records each call of get_four_center_batch_size from a module.

        With little_memory, it sees one byte free (through
        get_available_memory, the name batchsize holds), which leaves room
        for one density in each batch.

        NOTE: the free memory is changed only while get_four_center_batch_size
        runs.
        """

        calls = []
        real = module.get_four_center_batch_size

        def recorded(input_batch_size, n_total, n_ao, comm, call_memory=0):
            if little_memory:
                with monkeypatch.context() as patch:
                    patch.setattr(batchsize, 'get_available_memory',
                                  lambda: 1)
                    size = real(input_batch_size, n_total, n_ao, comm,
                                call_memory)
            else:
                size = real(input_batch_size, n_total, n_ao, comm,
                            call_memory)
            calls.append((n_total, call_memory, size))
            return size

        monkeypatch.setattr(module, 'get_four_center_batch_size', recorded)

        return calls

    @staticmethod
    def assert_cut(calls):
        """At least one build was given the memory of its exchange-correlation
        integration and cut into batches of fewer densities than it builds."""

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            assert any(memory is not None and memory > 0 and size < n_total
                       for n_total, memory, size in calls)

    @staticmethod
    def run_eigensolver(solver_type, molecule, basis, scf_results):

        solver = solver_type()
        solver.ostream.mute()
        solver.nstates = 3
        solver.conv_thresh = 1.0e-6
        results = solver.compute(molecule, basis, scf_results)

        return results

    def test_rpa(self, monkeypatch):

        molecule, basis = self.get_water()
        scf_results = self.run_scf(molecule, basis, 'b3lyp')

        whole = self.run_eigensolver(LinearResponseEigenSolver, molecule,
                                     basis, scf_results)

        calls = self.spy(monkeypatch, linearsolver, little_memory=True)

        cut = self.run_eigensolver(LinearResponseEigenSolver, molecule, basis,
                                   scf_results)

        self.assert_cut(calls)

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            assert np.max(np.abs(whole['eigenvalues'] -
                                 cut['eigenvalues'])) < 1.0e-10

    def test_unrestricted_rpa(self, monkeypatch):

        molecule, basis = self.get_ch3()
        scf_results = self.run_scf(molecule,
                                   basis,
                                   'b3lyp',
                                   unrestricted=True)

        whole = self.run_eigensolver(LinearResponseUnrestrictedEigenSolver,
                                     molecule, basis, scf_results)

        calls = self.spy(monkeypatch, linearsolver, little_memory=True)

        cut = self.run_eigensolver(LinearResponseUnrestrictedEigenSolver,
                                   molecule, basis, scf_results)

        self.assert_cut(calls)

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            assert np.max(np.abs(whole['eigenvalues'] -
                                 cut['eigenvalues'])) < 1.0e-10

    def test_cpp(self, monkeypatch):

        molecule, basis = self.get_water()
        scf_results = self.run_scf(molecule, basis, 'b3lyp')

        def run():
            solver = ComplexResponseSolver()
            solver.ostream.mute()
            solver.frequencies = [0.1, 0.2]
            solver.conv_thresh = 1.0e-6
            return solver.compute(molecule, basis, scf_results)

        whole = run()

        calls = self.spy(monkeypatch, linearsolver, little_memory=True)

        cut = run()

        self.assert_cut(calls)

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            for key, value in whole['response_functions'].items():
                assert abs(value - cut['response_functions'][key]) < 1.0e-9

    def test_cphf(self, monkeypatch):

        molecule, basis = self.get_water()
        scf_results = self.run_scf(molecule, basis, 'b3lyp')

        def run():
            solver = HessianOrbitalResponse()
            solver.update_settings({'conv_thresh': 1.0e-7})
            solver.ostream.mute()
            solver.compute(molecule, basis, scf_results)
            coefficients = solver.cphf_results['dist_cphf_ov']
            vectors = [c.get_full_vector(0) for c in coefficients]
            return np.array(vectors) if MPI.COMM_WORLD.Get_rank(
            ) == mpi_master() else None

        whole = run()

        calls = self.spy(monkeypatch, cphfsolver, little_memory=True)

        cut = run()

        self.assert_cut(calls)

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            assert np.max(np.abs(whole - cut)) < 1.0e-10

    def test_quadratic_response(self, monkeypatch):

        molecule, basis = self.get_water()
        scf_results = self.run_scf(molecule, basis, 'b3lyp')

        def run():
            driver = QuadraticResponseDriver()
            driver.ostream.mute()
            driver.update_settings({
                'conv_thresh': 1.0e-6,
                'b_frequencies': [0.1, 0.2],
                'c_frequencies': [0.2, 0.3],
                'damping': 0.1,
                'a_component': 'z',
                'b_component': 'z',
                'c_component': 'x',
            })
            return driver.compute(molecule, basis, scf_results)

        whole = run()

        calls = self.spy(monkeypatch, nonlinearsolver, little_memory=True)

        cut = run()

        self.assert_cut(calls)

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            keys = [k for k in whole if isinstance(k, tuple) and k[0] == 'qrf']
            assert len(keys) == 2
            for key in keys:
                assert abs(whole[key] - cut[key]) < 1.0e-9

    def test_ri_builds_do_not_take_it(self, monkeypatch):

        molecule, basis = self.get_water()
        scf_results = self.run_scf(molecule, basis, 'b3lyp', ri=True)

        calls = self.spy(monkeypatch, linearsolver)

        self.run_eigensolver(LinearResponseEigenSolver, molecule, basis,
                             scf_results)

        assert calls == []
