from collections import namedtuple
from mpi4py import MPI
import pytest

from veloxchem.batchsize import get_ri_batch_size, get_batch_size
import veloxchem.batchsize as batchsize


@pytest.mark.solvers
class TestRIBatchSize:

    comm = MPI.COMM_WORLD

    nao = 100

    nthreads = 4

    @staticmethod
    def set_memory(monkeypatch, available, nthreads):

        # NOTE: the batch size reads the node's memory through psutil and the
        # threads from OMP_NUM_THREADS; both are set here, the same on every
        # rank, and SLURM_NTASKS_PER_NODE to 1, so that the memory is not
        # divided among the ranks.
        memory = namedtuple('memory', ['total', 'available'])

        monkeypatch.setattr(batchsize.psutil, 'virtual_memory',
                            lambda: memory(64 * 1024**3, available))
        monkeypatch.setenv('OMP_NUM_THREADS', str(nthreads))
        monkeypatch.setenv('SLURM_NTASKS_PER_NODE', '1')

    def estimate(self):
        """What get_batch_size takes for each density."""

        return self.nao**2 * 8 * self.nthreads

    def test_a_density_memory_below_the_estimate_changes_nothing(
            self, monkeypatch):

        self.set_memory(monkeypatch, 8 * 1024**3, self.nthreads)

        plain = get_batch_size(None, 10**9, self.nao, self.comm)

        assert get_ri_batch_size(None, 10**9, self.nao, None,
                                 self.comm) == plain
        assert get_ri_batch_size(None, 10**9, self.nao,
                                 self.estimate() // 2, self.comm) == plain

    def test_a_density_memory_above_the_estimate_makes_the_batch_smaller(
            self, monkeypatch):

        available = 8 * 1024**3
        self.set_memory(monkeypatch, available, self.nthreads)

        density_memory = 3 * self.estimate()

        batch_size = get_ri_batch_size(None, 10**9, self.nao, density_memory,
                                       self.comm)

        assert batch_size == available // density_memory
        assert batch_size < get_batch_size(None, 10**9, self.nao, self.comm)

    def test_the_input_batch_size_and_the_items_bound_the_count(
            self, monkeypatch):

        self.set_memory(monkeypatch, 8 * 1024**3, self.nthreads)

        density_memory = 3 * self.estimate()

        assert get_ri_batch_size(9, 10**9, self.nao, density_memory,
                                 self.comm) == 9
        assert get_ri_batch_size(None, 5, self.nao, density_memory,
                                 self.comm) == 5
        assert get_ri_batch_size(None, 0, self.nao, density_memory,
                                 self.comm) == 0
        assert get_ri_batch_size(None, 10**9, self.nao, 1024**4,
                                 self.comm) == 1
