from veloxchem.veloxchemlib import number_of_batches
from veloxchem.veloxchemlib import batch_range


class TestBatchFunc:
    """
    Implements tests for src/general/BatchFunc.hpp
    """

    def test_number_of_batches(self):

        assert number_of_batches(10, 3) == 4
        assert number_of_batches(10, 6) == 2
        assert number_of_batches(10, 2) == 5
        assert number_of_batches(10, 5) == 2
        assert number_of_batches(10, 12) == 1

    def test_batch_range(self):

        assert batch_range(0, 10, 3) == (0, 3)
        assert batch_range(1, 10, 3) == (3, 6)
        assert batch_range(2, 10, 3) == (6, 9)
        assert batch_range(3, 10, 3) == (9, 10)
        assert batch_range(4, 10, 3) == (10, 10)
        assert batch_range(0, 10, 5) == (0, 5)
        assert batch_range(1, 10, 5) == (5, 10)
        assert batch_range(2, 10, 5) == (10, 10)
        assert batch_range(0, 3, 4) == (0, 3)
        assert batch_range(1, 3, 4) == (3, 3)

    def test_batch_range_with_offset(self):

        assert batch_range(0, 10, 3, 2) == (2, 5)
        assert batch_range(1, 10, 3, 2) == (5, 8)
        assert batch_range(2, 10, 3, 2) == (8, 11)
        assert batch_range(3, 10, 3, 2) == (11, 12)
        assert batch_range(4, 10, 3, 2) == (12, 12)
        assert batch_range(0, 10, 5, 1) == (1, 6)
        assert batch_range(1, 10, 5, 1) == (6, 11)
        assert batch_range(2, 10, 5, 1) == (11, 11)
        assert batch_range(0, 3, 4, 2) == (2, 5)
        assert batch_range(1, 3, 4, 2) == (5, 5)
