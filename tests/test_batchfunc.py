from veloxchem.veloxchemlib import get_batch_index
from veloxchem.veloxchemlib import number_of_batches
from veloxchem.veloxchemlib import get_batch_range


class TestBatchFunc:
    """
    Implements tests for src/general/BatchFunc.hpp
    """

    def test_get_batch_index(self):

        # nelements < nbatches
        assert get_batch_index(0, 3, 4) == 0
        assert get_batch_index(1, 3, 4) == 1
        assert get_batch_index(2, 3, 4) == 2
        assert get_batch_index(3, 3, 4) == 3
        assert get_batch_index(4, 3, 4) == 3

        # nelements < nbatches
        assert get_batch_index(0, 14, 4) == 0
        assert get_batch_index(1, 14, 4) == 4
        assert get_batch_index(2, 14, 4) == 8
        assert get_batch_index(3, 14, 4) == 11
        assert get_batch_index(4, 14, 4) == 14
        assert get_batch_index(5, 14, 4) == 14

    def test_number_of_batches(self):

        assert number_of_batches(10, 3) == 4
        assert number_of_batches(10, 6) == 2
        assert number_of_batches(10, 2) == 5
        assert number_of_batches(10, 5) == 2
        assert number_of_batches(10, 12) == 1

    def test_get_batch_range(self):

        assert get_batch_range(0, 10, 3) == (0, 3)
        assert get_batch_range(1, 10, 3) == (3, 6)
        assert get_batch_range(2, 10, 3) == (6, 9)
        assert get_batch_range(3, 10, 3) == (9, 10)
        assert get_batch_range(4, 10, 3) == (10, 10)

        assert get_batch_range(0, 10, 5) == (0, 5)
        assert get_batch_range(1, 10, 5) == (5, 10)
