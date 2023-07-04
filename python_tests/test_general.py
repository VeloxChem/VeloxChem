from veloxchem.veloxchemlib import mpi_master
from veloxchem.veloxchemlib import mpi_initialized
from veloxchem.veloxchemlib import bohr_in_angstrom
from veloxchem.veloxchemlib import hartree_in_ev
from veloxchem.veloxchemlib import to_angular_momentum
from veloxchem.errorhandler import assert_msg_critical


class TestGeneral:

    def test_mpi_master(self):

        assert mpi_master() == 0

    def test_mpi_initialized(self):

        assert mpi_initialized()

    def test_assert_msg_critical(self):

        assert_msg_critical(True, '')

    def test_constants(self):

        assert abs(bohr_in_angstrom() - 0.529177) < 1.0e-6
        assert abs(hartree_in_ev() - 27.2114) < 1.0e-4

    def test_angular_momentum(self):

        angmoms = 'SPDFGH'

        for ind, ang in enumerate(angmoms):
            assert to_angular_momentum(ind) == ang
            assert to_angular_momentum(ang) == ind
