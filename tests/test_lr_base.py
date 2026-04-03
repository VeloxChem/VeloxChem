from mpi4py import MPI
import numpy as np
import pytest

from veloxchem.veloxchemlib import mpi_master
from veloxchem.distributedarray import DistributedArray
from veloxchem.lrsolver import LinearResponseSolver
from veloxchem.outputstream import OutputStream


class RecordingOutput:

    def __init__(self):
        self.headers = []
        self.infos = []
        self.warnings = []
        self.blank_count = 0
        self.flushed = False

    def print_header(self, text):
        self.headers.append(text)

    def print_info(self, text):
        self.infos.append(text)

    def print_warning(self, text):
        self.warnings.append(text)

    def print_blank(self):
        self.blank_count += 1

    def flush(self):
        self.flushed = True


@pytest.mark.solvers
class TestLRBaseUtilities:

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason='skip pytest.raises for multiple MPI processes')
    def test_set_lr_property_and_reject_invalid_property(self):

        lr_drv = LinearResponseSolver()

        lr_drv.set_lr_property('Polarizability')

        assert lr_drv.property == 'polarizability'
        assert lr_drv.a_operator == 'electric dipole'
        assert lr_drv.b_operator == 'electric dipole'
        assert lr_drv.a_components == 'xyz'
        assert lr_drv.b_components == 'xyz'

        with pytest.raises(AssertionError,
                           match='LinearResponseSolver: invalid LR property'):
            lr_drv.set_lr_property('hyperpolarizability')

    def test_solution_vector_preconditioning_and_trial_selection(self):

        lr_drv = LinearResponseSolver()

        nrows = 16
        x_ger = np.arange(1.0, nrows + 1.0)
        x_ung = np.arange(21.0, 21.0 + nrows)
        solution = DistributedArray(
            np.column_stack((x_ger, x_ung)),
            lr_drv.comm,
        )
        full_vector = lr_drv.get_full_solution_vector(solution)

        if lr_drv.rank == mpi_master():
            expected_full_vector = np.hstack((x_ger + x_ung, x_ger - x_ung))
            assert np.allclose(full_vector, expected_full_vector)

        pa = np.arange(2.0, 2.0 + nrows)
        pb = np.ones(nrows)
        precond = DistributedArray(
            np.column_stack((pa, pb)),
            lr_drv.comm,
        )
        v_rg = np.arange(1.0, nrows + 1.0)
        v_ru = np.arange(31.0, 31.0 + nrows)
        v_in = DistributedArray(
            np.column_stack((v_rg, v_ru)),
            lr_drv.comm,
        )

        v_out = lr_drv._preconditioning(precond, v_in)
        full_v_out = v_out.get_full_matrix()

        if lr_drv.rank == mpi_master():
            expected_v_out = np.column_stack((pa * v_rg + pb * v_ru,
                                              pb * v_rg + pa * v_ru))
            assert np.allclose(full_v_out, expected_v_out)

        lr_drv.norm_thresh = 0.5
        trial_vectors = {
            ('x', 0.0): DistributedArray(
                np.column_stack((np.arange(1.0, nrows + 1.0), np.zeros(nrows))),
                lr_drv.comm,
            ),
            ('y', 0.1): DistributedArray(
                np.column_stack((np.zeros(nrows), np.arange(101.0, 101.0 + nrows))),
                lr_drv.comm,
            ),
            ('z', 0.2): DistributedArray(
                np.zeros((nrows, 2)),
                lr_drv.comm,
            ),
        }
        diag_precond = {
            w: DistributedArray(
                np.column_stack((np.ones(nrows), np.zeros(nrows))),
                lr_drv.comm,
            )
            for w in (0.0, 0.1, 0.2)
        }

        new_ger, new_ung = lr_drv._precond_trials(trial_vectors, diag_precond)
        full_new_ger = new_ger.get_full_matrix()
        full_new_ung = new_ung.get_full_matrix()

        if lr_drv.rank == mpi_master():
            assert np.allclose(full_new_ger,
                               np.arange(1.0, nrows + 1.0).reshape(-1, 1))
            assert np.allclose(full_new_ung,
                               np.arange(101.0, 101.0 + nrows).reshape(-1, 1))

    def test_print_helpers_emit_polarizability_output(self, monkeypatch):

        lr_drv = LinearResponseSolver()
        lr_drv.property = 'polarizability'
        lr_drv.frequencies = [0.0]
        lr_drv.a_components = 'xy'
        lr_drv.b_components = 'xz'
        lr_drv.print_level = 2
        lr_drv._cur_iter = 1

        iteration_stream = RecordingOutput()
        lr_drv._ostream = iteration_stream

        lr_drv._print_iteration({('x', 0.0): 1.0e-4, ('z', 0.0): 2.0e-4},
                                [('x', 0.0, 1.23), ('z', 0.0, -0.45)])

        assert any('Iteration:   2' in header
                   for header in iteration_stream.headers)
        assert any('<<x;x>>_0.0000' in header
                   for header in iteration_stream.headers)
        assert iteration_stream.flushed

        result_stream = RecordingOutput()
        rsp_funcs = {
            ('x', 'x', 0.0): -1.0,
            ('x', 'z', 0.0): -2.0,
            ('y', 'x', 0.0): -3.0,
            ('y', 'z', 0.0): -4.0,
        }

        lr_drv._print_results(rsp_funcs, result_stream)

        assert any('Polarizability (w=0.0000)' in header
                   for header in result_stream.headers)
        assert any('X' in header and 'Z' in header
                   for header in result_stream.headers)

        property_stream = RecordingOutput()
        monkeypatch.setattr(OutputStream, 'create_mpi_ostream',
                            staticmethod(lambda comm: property_stream))

        lr_drv.print_property(rsp_funcs)

        assert property_stream.headers == result_stream.headers
