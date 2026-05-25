import pytest
from mpi4py import MPI

from veloxchem.errorhandler import assert_msg_critical
from veloxchem.errorhandler import VeloxChemError


class TestErrorHandler:

    def test_assert_msg_critical_accepts_true_condition(self):

        assert_msg_critical(True, "unused")

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason="pytest.raises only valid in serial")
    def test_assert_msg_critical_raises_for_false_condition_in_serial(self):

        with pytest.raises(VeloxChemError, match="expected failure"):
            assert_msg_critical(False, "expected failure")
