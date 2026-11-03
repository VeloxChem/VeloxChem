from collections import namedtuple
from mpi4py import MPI
import platform
import pytest

from veloxchem.veloxchemlib import mpi_master
from veloxchem.batchsize import get_ri_batch_size, get_batch_size
import veloxchem.batchsize as batchsize


@pytest.mark.solvers
class TestRIBatchSize:

    comm = MPI.COMM_WORLD

    nao = 100

    nthreads = 4

    def set_memory(self, monkeypatch, available):
        """Sets the memory get_ri_batch_size reads as free on every rank.

        NOTE: get_ri_batch_size reads it through get_available_memory and the
        threads through get_number_of_threads, the names batchsize holds."""

        monkeypatch.setattr(batchsize, 'get_available_memory',
                            lambda: available)
        monkeypatch.setattr(batchsize, 'get_number_of_threads',
                            lambda: self.nthreads)

        here = platform.node()

        return available // self.comm.allgather(here).count(here)

    def estimate(self):
        """The memory of a density get_batch_size takes, per thread."""

        return self.nao**2 * 8 * self.nthreads

    def test_the_count_follows_the_room(self, monkeypatch):

        room = self.set_memory(monkeypatch, 8 * 1024**3)

        density_memory = 3 * self.estimate()

        assert get_ri_batch_size(None, 10**9, self.nao, density_memory,
                                 self.comm) == room // density_memory

        # NOTE: below the estimate the estimate is taken.
        assert get_ri_batch_size(None, 10**9, self.nao, self.estimate() // 2,
                                 self.comm) == room // self.estimate()

    def test_the_input_batch_size_and_the_items_bound_the_count(
            self, monkeypatch):

        self.set_memory(monkeypatch, 8 * 1024**3)

        density_memory = 3 * self.estimate()

        assert get_ri_batch_size(9, 10**9, self.nao, density_memory,
                                 self.comm) == 9
        assert get_ri_batch_size(None, 5, self.nao, density_memory,
                                 self.comm) == 5
        assert get_ri_batch_size(None, 0, self.nao, density_memory,
                                 self.comm) == 0
        assert get_ri_batch_size(None, 10**9, self.nao, 1024**4,
                                 self.comm) == 1

    def test_the_rank_with_the_least_room_sets_the_count(self, monkeypatch):

        rank = self.comm.Get_rank()

        # NOTE: every rank sees different free memory; the smallest room,
        # over the ranks on its node, decides.
        room = self.set_memory(monkeypatch, (8 + rank) * 1024**3)

        density_memory = 3 * self.estimate()

        expected = min(self.comm.allgather(room)) // density_memory

        assert get_ri_batch_size(None, 10**9, self.nao, density_memory,
                                 self.comm) == expected

    def test_a_control_group_limit_is_followed(self, monkeypatch):
        """The RI batch follows the memory a control group leaves, which
        get_batch_size does not read."""

        # NOTE: psutil sees a node with 512 GB free, and the control group
        # leaves 8 GB of it.
        memory = namedtuple('memory', ['total', 'available'])
        monkeypatch.setattr(batchsize.psutil, 'virtual_memory',
                            lambda: memory(1024 * 1024**3, 512 * 1024**3))
        monkeypatch.setenv('OMP_NUM_THREADS', str(self.nthreads))
        monkeypatch.setenv('SLURM_NTASKS_PER_NODE', '1')

        room = self.set_memory(monkeypatch, 8 * 1024**3)

        density_memory = self.estimate()

        ri = get_ri_batch_size(None, 10**9, self.nao, density_memory,
                               self.comm)
        plain = get_batch_size(None, 10**9, self.nao, self.comm)

        assert ri == room // density_memory
        assert ri < plain

        if self.comm.Get_rank() == mpi_master():
            assert plain == 512 * 1024**3 // self.estimate()
