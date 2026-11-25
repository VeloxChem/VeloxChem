from mpi4py import MPI
from pathlib import Path
import io
import os
import sys
import types

from veloxchem.veloxchemlib import mpi_master
from veloxchem.outputstream import OutputStream


class TestStreams:

    def run_ostream_test(self, ostream, state):

        ostream.print_line('')
        ostream.print_title('')
        ostream.print_separator()
        ostream.print_blank()
        assert ostream.get_state() is state

    def test_file(self):

        here = Path(__file__).parent
        outfile = here / 'data' / 'dummy.out'

        ostream = OutputStream(str(outfile))
        self.run_ostream_test(ostream, True)

        ostream.close()
        assert not ostream.get_state()

    def test_path(self):

        here = Path(__file__).parent
        outfile = here / 'data' / 'dummy.out'

        ostream = OutputStream(outfile)
        self.run_ostream_test(ostream, True)

        ostream.close()
        assert not ostream.get_state()

    def test_file_uses_utf8(self, tmp_path, monkeypatch):

        outfile = tmp_path / 'unicode.out'
        builtin_open = open
        output_encoding = None

        def tracking_open(*args, **kwargs):
            nonlocal output_encoding
            output_encoding = kwargs.get('encoding')
            return builtin_open(*args, **kwargs)

        monkeypatch.setattr('builtins.open', tracking_open)

        ostream = OutputStream(outfile)
        assert output_encoding == 'utf-8'

        reference = 'SoftwareX 7, 1–5 (2018)'
        ostream.print_line(reference)
        ostream.close()

        assert outfile.read_text(encoding='utf-8').rstrip() == reference

    def test_stdout(self):

        ostream = OutputStream(sys.stdout)
        assert ostream.get_state()

    def test_stdout_with_unencodable_character(self, monkeypatch):

        raw = io.BytesIO()
        ascii_stdout = io.TextIOWrapper(raw, encoding='ascii')
        monkeypatch.setattr(sys, 'stdout', ascii_stdout)

        ostream = OutputStream(sys.stdout)
        ostream.print_line('SoftwareX 7, 1–5 (2018)')
        ostream.flush()

        assert b'SoftwareX 7, 1?5 (2018)' in raw.getvalue()

    def test_default(self):

        ostream = OutputStream()
        assert ostream.get_state()

    def test_none(self):

        ostream = OutputStream(None)
        self.run_ostream_test(ostream, False)

    def test_mpi(self):

        comm = MPI.COMM_WORLD
        ostream = OutputStream.create_mpi_ostream(comm)

        if comm.Get_rank() == mpi_master():
            assert ostream.get_state()
        else:
            assert not ostream.get_state()

        ostream.close()

    @staticmethod
    def get_blas_threads_text(tmp_path, monkeypatch, n_omp, blas_pools,
                              main_cpus):
        """
        Returns what print_blas_threads_info writes with n_omp OpenMP
        threads, a fake threadpoolctl which lists blas_pools (None for
        threadpoolctl not installed) and a main thread which may run on
        main_cpus logical cores.
        """

        monkeypatch.setattr('veloxchem.outputstream.get_number_of_threads',
                            lambda: n_omp)
        monkeypatch.setattr('os.sched_getaffinity',
                            lambda pid: set(range(main_cpus)),
                            raising=False)

        if blas_pools is None:
            # NOTE: None in sys.modules makes the import raise ImportError.
            monkeypatch.setitem(sys.modules, 'threadpoolctl', None)
        else:
            fake = types.ModuleType('threadpoolctl')
            fake.threadpool_info = lambda: list(blas_pools)
            monkeypatch.setitem(sys.modules, 'threadpoolctl', fake)

        outfile = tmp_path / 'blas.out'
        ostream = OutputStream(outfile)
        ostream.print_blas_threads_info()
        ostream.close()

        return outfile.read_text(encoding='utf-8')

    @staticmethod
    def get_blas_pool(num_threads, threading_layer):

        return {
            'user_api': 'blas',
            'internal_api': 'openblas',
            'num_threads': num_threads,
            'filepath': '/lib/libopenblas.so',
            'version': '0.3.27',
            'threading_layer': threading_layer,
        }

    def test_blas_threads_fewer_than_openmp(self, tmp_path, monkeypatch):

        text = self.get_blas_threads_text(
            tmp_path, monkeypatch, 16, [self.get_blas_pool(2, 'pthreads')], 2)

        assert ('BLAS library libopenblas.so: openblas 0.3.27, pthreads, 2 ' +
                'threads.') in text
        assert '* Warning *' in text
        assert 'fewer threads than the 16 OpenMP threads' in text
        assert 'libopenblas.so (2)' in text
        if hasattr(os, 'sched_getaffinity'):
            assert 'The main thread may run on 2 logical cores' in text

    def test_blas_threads_not_threaded(self, tmp_path, monkeypatch):

        text = self.get_blas_threads_text(
            tmp_path, monkeypatch, 16, [self.get_blas_pool(1, 'disabled')],
            32)

        assert 'libopenblas.so (1)' in text
        # NOTE: the main thread may run on more logical cores than OpenMP
        # threads, so the sentence on the binding of the main thread is left
        # out.
        assert 'The main thread may run on' not in text

    def test_blas_threads_as_openmp(self, tmp_path, monkeypatch):

        # NOTE: an OpenMP BLAS reports the OpenMP threads even where the main
        # thread is bound to one core.
        text = self.get_blas_threads_text(
            tmp_path, monkeypatch, 16, [self.get_blas_pool(16, 'openmp')], 2)

        assert 'openblas 0.3.27, openmp, 16 threads.' in text
        assert '* Warning *' not in text
        assert 'logical cores' not in text

    def test_blas_threads_without_threadpoolctl(self, tmp_path, monkeypatch):

        text = self.get_blas_threads_text(tmp_path, monkeypatch, 16, None, 2)

        assert text.strip() == ''
