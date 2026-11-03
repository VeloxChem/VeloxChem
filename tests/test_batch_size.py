from collections import namedtuple
from mpi4py import MPI
import platform
import pytest

from veloxchem.batchsize import get_batch_size
import veloxchem.batchsize as batchsize


@pytest.mark.solvers
class TestBatchSize:

    comm = MPI.COMM_WORLD

    nao = 100

    @staticmethod
    def memory(monkeypatch, available, nthreads):

        # NOTE: get_batch_size reads the node's memory through psutil and the
        # threads from OMP_NUM_THREADS; both are set here, the same on every
        # rank.
        memory = namedtuple('memory', ['total', 'available'])

        monkeypatch.setattr(batchsize.psutil, 'virtual_memory',
                            lambda: memory(64 * 1024**3, available))
        monkeypatch.setenv('OMP_NUM_THREADS', str(nthreads))

    def per_density(self, nthreads):

        return self.nao**2 * 8 * nthreads

    def test_the_ranks_on_the_node_divide_its_memory(self, monkeypatch):

        monkeypatch.delenv('SLURM_NTASKS_PER_NODE', raising=False)

        available = 8 * 1024**3
        self.memory(monkeypatch, available, 4)

        here = platform.node()
        on_this_node = self.comm.allgather(here).count(here)

        expected = int((available // on_this_node) / self.per_density(4))

        assert get_batch_size(None, 10**9, self.nao, self.comm) == expected

    def test_slurm_ntasks_per_node_takes_priority(self, monkeypatch):

        # NOTE: a count which is not the number of ranks on this node, so that
        # the result tells which of the two was used.
        monkeypatch.setenv('SLURM_NTASKS_PER_NODE', '7')

        available = 8 * 1024**3
        self.memory(monkeypatch, available, 4)

        expected = int((available // 7) / self.per_density(4))

        assert get_batch_size(None, 10**9, self.nao, self.comm) == expected

    def test_the_input_batch_size_and_the_items_bound_the_count(
            self, monkeypatch):

        monkeypatch.delenv('SLURM_NTASKS_PER_NODE', raising=False)

        self.memory(monkeypatch, 8 * 1024**3, 4)

        assert get_batch_size(9, 10**9, self.nao, self.comm) == 9
        assert get_batch_size(None, 5, self.nao, self.comm) == 5
        assert get_batch_size(None, 0, self.nao, self.comm) == 0
