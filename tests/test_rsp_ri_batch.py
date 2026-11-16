from mpi4py import MPI
import numpy as np
import pytest

from veloxchem.veloxchemlib import mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.scfunrestdriver import ScfUnrestrictedDriver
from veloxchem.lreigensolver import LinearResponseEigenSolver
from veloxchem.lreigensolverunrest import LinearResponseUnrestrictedEigenSolver
from veloxchem.cppsolver import ComplexResponseSolver
from veloxchem.hessianorbitalresponse import HessianOrbitalResponse
from veloxchem.quadraticresponsedriver import QuadraticResponseDriver
from veloxchem.cubicresponsedriver import CubicResponseDriver
from veloxchem.thgdriver import ThgDriver
from veloxchem.outputstream import OutputStream
import veloxchem.batchsize as batchsize
import veloxchem.linearsolver as linearsolver
import veloxchem.cphfsolver as cphfsolver
import veloxchem.nonlinearsolver as nonlinearsolver


@pytest.mark.solvers
class TestResponseRiBatch:
    """RI Fock builds take their batch size from get_ri_batch_size, and a
    build cut into many batches gives what one batch gives."""

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
    def run_scf(molecule, basis, xcfun, ri=True, unrestricted=False):

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
        """Records each call of get_ri_batch_size from a module.

        With little_memory, get_ri_batch_size sees one byte free (through
        get_available_memory, the name batchsize holds), which leaves room
        for one density in each batch.

        NOTE: the free memory is changed only while get_ri_batch_size runs.
        The other readers of it (get_batch_size in the drivers which call it
        themselves, the memory budgets of the RI drivers) see the machine's.
        """

        calls = []
        real = module.get_ri_batch_size

        def recorded(input_batch_size,
                     n_total,
                     n_ao,
                     density_memory,
                     comm,
                     call_memory=0):
            if little_memory:
                with monkeypatch.context() as patch:
                    patch.setattr(batchsize, 'get_available_memory',
                                  lambda: 1)
                    size = real(input_batch_size, n_total, n_ao,
                                density_memory, comm, call_memory)
            else:
                size = real(input_batch_size, n_total, n_ao, density_memory,
                            comm, call_memory)
            calls.append((n_total, density_memory, size))
            return size

        monkeypatch.setattr(module, 'get_ri_batch_size', recorded)

        return calls

    @staticmethod
    def count_batches(monkeypatch, module):
        """Records each number of batches a module works out."""

        counts = []
        real = module.get_number_of_batches

        def recorded(n_total, batch_size, comm):
            count = real(n_total, batch_size, comm)
            counts.append(count)
            return count

        monkeypatch.setattr(module, 'get_number_of_batches', recorded)

        return counts

    @staticmethod
    def assert_cut(calls):
        """At least one RI build was given a memory of its densities and cut
        into batches of fewer densities than it builds."""

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

    @pytest.mark.parametrize('xcfun', ['b3lyp', 'blyp'])
    def test_rpa(self, monkeypatch, xcfun):

        molecule, basis = self.get_water()
        scf_results = self.run_scf(molecule, basis, xcfun)

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
        batches = self.count_batches(monkeypatch, nonlinearsolver)

        cut = run()

        self.assert_cut(calls)

        # NOTE: a batch with a functional is made of whole groups, one for
        # each pair of frequencies; two pairs can be cut into two batches.
        assert max(batches) > 1

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            keys = [k for k in whole if isinstance(k, tuple) and k[0] == 'qrf']
            assert len(keys) == 2
            for key in keys:
                assert abs(whole[key] - cut[key]) < 1.0e-9

    def test_cubic_response(self, monkeypatch):

        molecule, basis = self.get_water()
        scf_results = self.run_scf(molecule, basis, 'b3lyp')


        def run():
            driver = CubicResponseDriver()
            driver.ostream.mute()
            driver.update_settings(
                {
                    'conv_thresh': 1.0e-6,
                    'b_frequencies': [0.11, 0.12],
                    'c_frequencies': [-0.2, -0.21],
                    'd_frequencies': [0.06, 0.07],
                    'a_component': 'z',
                    'b_component': 'z',
                    'c_component': 'z',
                    'd_component': 'z',
                    'damping': 0.1,
                }, {'xcfun': 'b3lyp'})
            return driver.compute(molecule, basis, scf_results)

        whole = run()

        calls = self.spy(monkeypatch, nonlinearsolver, little_memory=True)
        batches = self.count_batches(monkeypatch, nonlinearsolver)

        cut = run()

        self.assert_cut(calls)
        assert max(batches) > 1

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            keys = [k for k in whole if isinstance(k, tuple) and k[0] == 'crf']
            assert len(keys) == 2
            for key in keys:
                assert abs(whole[key] - cut[key]) < 1.0e-9

    def test_factor_sizes_count_a_shared_array_once(self):

        left = np.zeros((4, 3))
        rights = [np.zeros((4, 2)) for _ in range(5)]
        transposed = [np.zeros((4, 1)) for _ in range(5)]

        assert nonlinearsolver._factor_sizes(None) == (0, 0)
        assert nonlinearsolver._factor_sizes((left, rights)) == (left.nbytes,
                                                                 2)
        assert nonlinearsolver._factor_sizes(
            (left, rights, transposed)) == (left.nbytes, 3)

        # NOTE: the two sets of a dictionary share their left factors.
        other = np.zeros((4, 6))
        split = {
            'second': (left, rights[:2], other, rights[:2]),
            'third': (left, rights[2:], other, transposed[2:]),
        }

        assert nonlinearsolver._factor_sizes(split) == (left.nbytes +
                                                        other.nbytes, 4)

    @staticmethod
    def run_thg(molecule, basis, scf_results, xcfun, frequencies):

        driver = ThgDriver(MPI.COMM_WORLD, OutputStream(None))
        driver.ostream.mute()
        driver.update_settings(
            {
                'frequencies': frequencies,
                'damping': 0.1,
                'conv_thresh': 1.0e-6,
            }, {'xcfun': xcfun})

        return driver.compute(molecule, basis, scf_results)

    def test_nonlinear_batch_memory_holds_the_factors(self, monkeypatch):

        molecule, basis = self.get_water()
        scf_results = self.run_scf(molecule, basis, 'b3lyp')

        xc_memory = []
        real_xc = nonlinearsolver.get_response_fock_xc_call_memory

        def recorded_xc(*args, **kwargs):
            memory = real_xc(*args, **kwargs)
            xc_memory.append(memory)
            return memory

        columns = []
        real_density = nonlinearsolver.get_ri_density_memory

        def recorded_density(n_ao, n_matrices, n_factor_columns, n_sliced,
                             xc_density_memory):
            columns.append(n_factor_columns)
            return real_density(n_ao, n_matrices, n_factor_columns, n_sliced,
                                xc_density_memory)

        calls = []
        real = nonlinearsolver.get_ri_batch_size

        def recorded(input_batch_size,
                     n_total,
                     n_ao,
                     density_memory,
                     comm,
                     call_memory=0):
            calls.append((n_total, call_memory))
            return real(input_batch_size, n_total, n_ao, density_memory,
                        comm, call_memory)

        monkeypatch.setattr(nonlinearsolver,
                            'get_response_fock_xc_call_memory', recorded_xc)
        monkeypatch.setattr(nonlinearsolver, 'get_ri_density_memory',
                            recorded_density)
        monkeypatch.setattr(nonlinearsolver, 'get_ri_batch_size', recorded)

        frequencies = [0.05, 0.06]

        self.run_thg(molecule, basis, scf_results, 'b3lyp', frequencies)

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            n_ao = basis.get_dimensions_of_basis()
            n_orb = scf_results['C_alpha'].shape[1]

            # NOTE: the first build of a THG takes 12 second-order and 6
            # third-order densities for each frequency. Each density has a
            # right factor of the occupied half and one of the virtual half,
            # n_orb columns together, and the two sets share the left factors,
            # n_ao by n_orb together.
            n_fock = 18 * len(frequencies)

            first = [i for i, (n_total, _) in enumerate(calls)
                     if n_total == n_fock]
            assert first

            for i in first:
                assert calls[i][1] == xc_memory[i] + 8 * n_ao * n_orb
                assert columns[i] == n_orb

            # NOTE: every build of the driver is handed factors.
            for i, (_, call_memory) in enumerate(calls):
                assert call_memory > xc_memory[i]
                assert columns[i] > 0

    @pytest.mark.parametrize('xcfun', ['b3lyp', 'hf'])
    def test_thg(self, monkeypatch, xcfun):

        molecule, basis = self.get_water()
        scf_results = self.run_scf(molecule, basis, xcfun)

        frequencies = [0.05, 0.06]

        whole = self.run_thg(molecule, basis, scf_results, xcfun, frequencies)

        calls = self.spy(monkeypatch, nonlinearsolver, little_memory=True)
        batches = self.count_batches(monkeypatch, nonlinearsolver)

        cut = self.run_thg(molecule, basis, scf_results, xcfun, frequencies)

        self.assert_cut(calls)

        # NOTE: with a functional a batch is made of whole groups, one for
        # each frequency; two frequencies can be cut into two batches.
        assert max(batches) > 1

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            assert len(whole['THG']) == len(frequencies)
            for key, value in whole['THG'].items():
                assert abs(value - cut['THG'][key]) < 1.0e-9 * max(
                    abs(value), 1.0)

    def test_four_center_builds_do_not_take_it(self, monkeypatch):

        molecule, basis = self.get_water()
        scf_results = self.run_scf(molecule, basis, 'b3lyp', ri=False)

        calls = self.spy(monkeypatch, linearsolver)

        self.run_eigensolver(LinearResponseEigenSolver, molecule, basis,
                             scf_results)

        assert calls == []
