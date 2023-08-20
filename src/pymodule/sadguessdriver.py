import numpy as np

from .veloxchemlib import OverlapDriver
from .errorhandler import assert_msg_critical


class SadGuessDriver:
    """
    Implements SAD initial guess.

    Instance variables
        - _num_unpaired_electrons_on_atoms: number of unparied electrons on
          each atom
    """

    def __init__(self):
        """
        Initializes the SAD guess driver.
        """

        self._num_unpaired_electrons_on_atoms = []

    def set_number_of_unpaired_electrons_on_atoms(self, num_unpaired_electrons):

        self._num_unpaired_electrons_on_atoms = list(num_unpaired_electrons)

    @staticmethod
    def _get_occ_1s(nocc):
        """
        Gets occupation numbers for 1s elements.

        :param nocc:
            Number of 1s orbitals.

        :return:
            List of occupation numbers.
        """

        occ = max(0.0, nocc)

        #       1s
        return [occ]

    @staticmethod
    def _get_occ_2s(nocc):
        """
        Gets occupation numbers for 2s elements.

        :param nocc:
            Number of 2s orbitals.

        :return:
            List of occupation numbers.
        """

        occ = max(0.0, nocc)

        #       1s   2s
        return [1.0, occ]

    @staticmethod
    def _get_occ_2s2p(nocc):
        """
        Gets occupation numbers for 2s2p elements.

        :param nocc:
            Number of 2s2p orbitals.

        :return:
            List of occupation numbers.
        """

        occ = max(0.0, nocc)

        #       1s   2s   2p-1 2p0  2p+1
        return [1.0, occ, occ, occ, occ]

    @staticmethod
    def _get_occ_3s(nocc):
        """
        Gets occupation numbers for 3s elements.

        :param nocc:
            Number of 3s orbitals.

        :return:
            List of occupation numbers.
        """

        occ = max(0.0, nocc)

        #       1s   2s   3s   2p-1 2p0  2p+1
        return [1.0, 1.0, occ, 1.0, 1.0, 1.0]

    @staticmethod
    def _get_occ_3s3p(nocc):
        """
        Gets occupation numbers for 3s3p elements.

        :param nocc:
            Number of 3s3p orbitals.

        :return:
            List of occupation numbers.
        """

        occ = max(0.0, nocc)

        #       1s   2s   3s   2p-1 3p-1 2p0  3p0  2p+1 3p+1
        return [1.0, 1.0, occ, 1.0, occ, 1.0, occ, 1.0, occ]

    @staticmethod
    def _get_occ_4s(nocc):
        """
        Gets occupation numbers for 4s elements.

        :param nocc:
            Number of 4s orbitals.

        :return:
            List of occupation numbers.
        """

        occ = max(0.0, nocc)

        #       1s   2s   3s   4s   2p-1 3p-1 2p0  3p0  2p+1 3p+1
        return [1.0, 1.0, 1.0, occ, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]

    @staticmethod
    def _get_occ_3d(nocc):
        """
        Gets occupation numbers for 3d elements.

        :param nocc:
            Number of 3d orbitals.

        :return:
            List of occupation numbers.
        """

        occ = max(0.0, nocc)

        #   1s   2s   3s   4s   2p-1 3p-1 2p0  3p0  2p+1 3p+1 3d-2 3d-1 3d0
        #   3d+1 3d+2
        return [
            1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, occ, occ, occ,
            occ, occ
        ]

    @staticmethod
    def _get_occ_4s4p(nocc):
        """
        Gets occupation numbers for 4s4p elements.

        :param nocc:
            Number of 4s4p orbitals.

        :return:
            List of occupation numbers.
        """

        occ = max(0.0, nocc)

        #   1s   2s   3s   4s   2p-1 3p-1 4p-1 2p0  3p0  4p0  2p+1 3p+1 4p+1
        #   3d-2 3d-1 3d0  3d+1 3d+2
        return [
            1.0, 1.0, 1.0, occ, 1.0, 1.0, occ, 1.0, 1.0, occ, 1.0, 1.0, occ,
            1.0, 1.0, 1.0, 1.0, 1.0
        ]

    @staticmethod
    def _get_occ_5s(nocc):
        """
        Gets occupation numbers for 5s elements.

        :param nocc:
            Number of 5s orbitals.

        :return:
            List of occupation numbers.
        """

        occ = max(0.0, nocc)

        #   1s   2s   3s   4s   5s   2p-1 3p-1 4p-1 2p0  3p0  4p0  2p+1 3p+1
        #   4p+1 3d-2 3d-1 3d0  3d+1 3d+2
        return [
            1.0, 1.0, 1.0, 1.0, occ, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
            1.0, 1.0, 1.0, 1.0, 1.0, 1.0
        ]

    @staticmethod
    def _get_occ_4d(nocc):
        """
        Gets occupation numbers for 4d elements.

        :param nocc:
            Number of 4d orbitals.

        :return:
            List of occupation numbers.
        """

        occ = max(0.0, nocc)

        #   1s   2s   3s   4s   5s   2p-1 3p-1 4p-1 2p0  3p0  4p0  2p+1 3p+1
        #   4p+1 3d-2 4d-2 3d-1 4d-1 3d0  4d0  3d+1 4d+1 3d+2 4d+2
        return [
            1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
            1.0, 1.0, occ, 1.0, occ, 1.0, occ, 1.0, occ, 1.0, occ
        ]

    @staticmethod
    def _get_occ_5s5p(nocc):
        """
        Gets occupation numbers for 5s5p elements.

        :param nocc:
            Number of 5s5p orbitals.

        :return:
            List of occupation numbers.
        """

        occ = max(0.0, nocc)

        #   1s   2s   3s   4s   5s   2p-1 3p-1 4p-1 5p-1 2p0  3p0  4p0  5p0
        #   2p+1 3p+1 4p+1 5p+1 3d-2 4d-2 3d-1 4d-1 3d0  4d0  3d+1 4d+1 3d+2 4d+2
        return [
            1.0, 1.0, 1.0, 1.0, occ, 1.0, 1.0, 1.0, occ, 1.0, 1.0, 1.0, occ,
            1.0, 1.0, 1.0, occ, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0
        ]

    @staticmethod
    def _get_occ_6s(nocc):
        """
        Gets occupation numbers for 6s elements.

        :param nocc:
            Number of 6s orbitals.

        :return:
            List of occupation numbers.
        """

        occ = max(0.0, nocc)

        #   1s   2s   3s   4s   5s   6s   2p-1 3p-1 4p-1 5p-1 2p0  3p0  4p0
        #   5p0  2p+1 3p+1 4p+1 5p+1 3d-2 4d-2 3d-1 4d-1 3d0  4d0  3d+1 4d+1
        #   3d+2 4d+2
        return [
            1.0, 1.0, 1.0, 1.0, 1.0, occ, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
            1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
            1.0, 1.0
        ]

    @staticmethod
    def _get_occ_4f(nocc):
        """
        Gets occupation numbers for 4f elements.

        :param nocc:
            Number of 4f orbitals.

        :return:
            List of occupation numbers.
        """

        occ = max(0.0, nocc)

        #   1s   2s   3s   4s   5s   6s   2p-1 3p-1 4p-1 5p-1 2p0  3p0  4p0
        #   5p0  2p+1 3p+1 4p+1 5p+1 3d-2 4d-2 3d-1 4d-1 3d0  4d0  3d+1 4d+1
        #   3d+2 4d+2 4f-3 4f-2 4f-1 4f0  4f+1 4f+2 4f+3
        return [
            1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
            1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
            1.0, 1.0, occ, occ, occ, occ, occ, occ, occ
        ]

    @staticmethod
    def _get_occ_5d(nocc):
        """
        Gets occupation numbers for 5d elements.

        :param nocc:
            Number of 5d orbitals.

        :return:
            List of occupation numbers.
        """

        occ = max(0.0, nocc)

        #   1s   2s   3s   4s   5s   6s   2p-1 3p-1 4p-1 5p-1 2p0  3p0  4p0
        #   5p0  2p+1 3p+1 4p+1 5p+1 3d-2 4d-2 5d-2 3d-1 4d-1 5d-1 3d0  4d0
        #   5d0  3d+1 4d+1 5d+1 3d+2 4d+2 5d+2 4f-3 4f-2 4f-1 4f0  4f+1 4f+2 4f+3
        return [
            1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
            1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, occ, 1.0, 1.0, occ, 1.0, 1.0,
            occ, 1.0, 1.0, occ, 1.0, 1.0, occ, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0
        ]

    @staticmethod
    def _get_occ_6s6p(nocc):
        """
        Gets occupation numbers for 6s6p elements.

        :param nocc:
            Number of 6s6p orbitals.

        :return:
            List of occupation numbers.
        """

        occ = max(0.0, nocc)

        #   1s   2s   3s   4s   5s   6s   2p-1 3p-1 4p-1 5p-1 6s-1 2p0  3p0
        #   4p0  5p0  6s0  2p+1 3p+1 4p+1 5p+1 6s+1 3d-2 4d-2 5d-2 3d-1 4d-1
        #   5d-1 3d0  4d0  5d0  3d+1 4d+1 5d+1 3d+2 4d+2 5d+2 4f-3 4f-2 4f-1
        #   4f0  4f+1 4f+2 4f+3
        return [
            1.0, 1.0, 1.0, 1.0, 1.0, occ, 1.0, 1.0, 1.0, 1.0, occ, 1.0, 1.0,
            1.0, 1.0, occ, 1.0, 1.0, 1.0, 1.0, occ, 1.0, 1.0, 1.0, 1.0, 1.0,
            1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
            1.0, 1.0, 1.0, 1.0
        ]

    @staticmethod
    def _get_occ_7s(nocc):
        """
        Gets occupation numbers for 7s elements.

        :param nocc:
            Number of 7s orbitals.

        :return:
            List of occupation numbers.
        """

        occ = max(0.0, nocc)

        #   1s   2s   3s   4s   5s   6s   7s   2p-1 3p-1 4p-1 5p-1 6s-1 2p0
        #   3p0  4p0  5p0  6s0  2p+1 3p+1 4p+1 5p+1 6s+1 3d-2 4d-2 5d-2 3d-1
        #   4d-1 5d-1 3d0  4d0  5d0  3d+1 4d+1 5d+1 3d+2 4d+2 5d+2 4f-3 4f-2
        #   4f-1 4f0  4f+1 4f+2 4f+3
        return [
            1.0, 1.0, 1.0, 1.0, 1.0, 1.0, occ, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
            1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
            1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
            1.0, 1.0, 1.0, 1.0, 1.0
        ]

    @staticmethod
    def _get_occ_5f(nocc):
        """
        Gets occupation numbers for 5f elements.

        :param nocc:
            Number of 5f orbitals.

        :return:
            List of occupation numbers.
        """

        occ = max(0.0, nocc)

        #   1s   2s   3s   4s   5s   6s   7s   2p-1 3p-1 4p-1 5p-1 6s-1 2p0
        #   3p0  4p0  5p0  6s0  2p+1 3p+1 4p+1 5p+1 6s+1 3d-2 4d-2 5d-2 3d-1
        #   4d-1 5d-1 3d0  4d0  5d0  3d+1 4d+1 5d+1 3d+2 4d+2 5d+2 4f-3 5f-3
        #   4f-2 5f-2 4f-1 5f-1 4f0  5f0  4f+1 5f+1 4f+2 5f+2 4f+3 5f+3
        return [
            1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
            1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
            1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, occ,
            1.0, occ, 1.0, occ, 1.0, occ, 1.0, occ, 1.0, occ, 1.0, occ
        ]

    def get_occupation_numbers_for_element(self, elem_id, nelec):
        """
        Computes occupation numbers for a given element.

        :param elem_id:
            The element id (nuclear charge).
        :param nelec:
            the number of excessive alpha or beta electrons.

        :return:
            A list of occupation numbers.
        """

        # H,He

        if elem_id == 1:
            return self._get_occ_1s(0.5 + nelec)

        elif elem_id == 2:
            return self._get_occ_1s(1.0 + nelec)

        # Li,Be

        elif elem_id == 3:
            return self._get_occ_2s(0.5 + nelec)

        elif elem_id == 4:
            return self._get_occ_2s(1.0 + nelec)

        # B,C,N,O,F,Ne

        elif elem_id == 5:
            return self._get_occ_2s2p(0.375 + nelec / 4.0)

        elif elem_id == 6:
            return self._get_occ_2s2p(0.500 + nelec / 4.0)

        elif elem_id == 7:
            return self._get_occ_2s2p(0.625 + nelec / 4.0)

        elif elem_id == 8:
            return self._get_occ_2s2p(0.750 + nelec / 4.0)

        elif elem_id == 9:
            return self._get_occ_2s2p(0.875 + nelec / 4.0)

        elif elem_id == 10:
            return self._get_occ_2s2p(1.000 + nelec / 4.0)

        # Na,Mg

        elif elem_id == 11:
            return self._get_occ_3s(0.5 + nelec)

        elif elem_id == 12:
            return self._get_occ_3s(1.0 + nelec)

        # Al,Si,P,S,Cl,Ar

        elif elem_id == 13:
            return self._get_occ_3s3p(0.375 + nelec / 4.0)

        elif elem_id == 14:
            return self._get_occ_3s3p(0.500 + nelec / 4.0)

        elif elem_id == 15:
            return self._get_occ_3s3p(0.625 + nelec / 4.0)

        elif elem_id == 16:
            return self._get_occ_3s3p(0.750 + nelec / 4.0)

        elif elem_id == 17:
            return self._get_occ_3s3p(0.875 + nelec / 4.0)

        elif elem_id == 18:
            return self._get_occ_3s3p(1.000 + nelec / 4.0)

        # K,Ca

        elif elem_id == 19:
            return self._get_occ_4s(0.5 + nelec)

        elif elem_id == 20:
            return self._get_occ_4s(1.0 + nelec)

        # Sc,Ti,V,Cr,Mn,Fe,Co,Ni,Cu,Zn

        elif elem_id == 21:
            return self._get_occ_3d(0.1 + nelec / 5.0)

        elif elem_id == 22:
            return self._get_occ_3d(0.2 + nelec / 5.0)

        elif elem_id == 23:
            return self._get_occ_3d(0.3 + nelec / 5.0)

        elif elem_id == 24:
            return self._get_occ_3d(0.4 + nelec / 5.0)

        elif elem_id == 25:
            return self._get_occ_3d(0.5 + nelec / 5.0)

        elif elem_id == 26:
            return self._get_occ_3d(0.6 + nelec / 5.0)

        elif elem_id == 27:
            return self._get_occ_3d(0.7 + nelec / 5.0)

        elif elem_id == 28:
            return self._get_occ_3d(0.8 + nelec / 5.0)

        elif elem_id == 29:
            return self._get_occ_3d(0.9 + nelec / 5.0)

        elif elem_id == 30:
            return self._get_occ_3d(1.0 + nelec / 5.0)

        # Ga,Ge,As,Se,Br,Kr

        elif elem_id == 31:
            return self._get_occ_4s4p(0.375 + nelec / 4.0)

        elif elem_id == 32:
            return self._get_occ_4s4p(0.500 + nelec / 4.0)

        elif elem_id == 33:
            return self._get_occ_4s4p(0.625 + nelec / 4.0)

        elif elem_id == 34:
            return self._get_occ_4s4p(0.750 + nelec / 4.0)

        elif elem_id == 35:
            return self._get_occ_4s4p(0.875 + nelec / 4.0)

        elif elem_id == 36:
            return self._get_occ_4s4p(1.000 + nelec / 4.0)

        # Rb,Sr

        elif elem_id == 37:
            return self._get_occ_5s(0.5 + nelec)

        elif elem_id == 38:
            return self._get_occ_5s(1.0 + nelec)

        # Y,Zr,Nb,Mo,Tc,Ru,Rh,Pd,Ag,Cd (39-48)

        elif elem_id == 39:
            return self._get_occ_4d(0.1 + nelec / 5.0)

        elif elem_id == 40:
            return self._get_occ_4d(0.2 + nelec / 5.0)

        elif elem_id == 41:
            return self._get_occ_4d(0.3 + nelec / 5.0)

        elif elem_id == 42:
            return self._get_occ_4d(0.4 + nelec / 5.0)

        elif elem_id == 43:
            return self._get_occ_4d(0.5 + nelec / 5.0)

        elif elem_id == 44:
            return self._get_occ_4d(0.6 + nelec / 5.0)

        elif elem_id == 45:
            return self._get_occ_4d(0.7 + nelec / 5.0)

        elif elem_id == 46:
            return self._get_occ_4d(0.8 + nelec / 5.0)

        elif elem_id == 47:
            return self._get_occ_4d(0.9 + nelec / 5.0)

        elif elem_id == 48:
            return self._get_occ_4d(1.0 + nelec / 5.0)

        # In,Sn,Sb,Te,I,Xe

        elif elem_id == 49:
            return self._get_occ_5s5p(0.375 + nelec / 4.0)

        elif elem_id == 50:
            return self._get_occ_5s5p(0.500 + nelec / 4.0)

        elif elem_id == 51:
            return self._get_occ_5s5p(0.625 + nelec / 4.0)

        elif elem_id == 52:
            return self._get_occ_5s5p(0.750 + nelec / 4.0)

        elif elem_id == 53:
            return self._get_occ_5s5p(0.875 + nelec / 4.0)

        elif elem_id == 54:
            return self._get_occ_5s5p(1.000 + nelec / 4.0)

        # Cs,Ba (55-56)

        elif elem_id == 55:
            return self._get_occ_6s(0.5 + nelec)

        elif elem_id == 56:
            return self._get_occ_6s(1.0 + nelec)

        # La,Ce,Pr,Nd,Pm,Sm,Eu,Gb,Tb,Dy,Ho,Er,Tm,Yb

        elif elem_id == 57:
            return self._get_occ_4f(1.0 / 14.0 + nelec / 7.0)

        elif elem_id == 58:
            return self._get_occ_4f(2.0 / 14.0 + nelec / 7.0)

        elif elem_id == 59:
            return self._get_occ_4f(3.0 / 14.0 + nelec / 7.0)

        elif elem_id == 60:
            return self._get_occ_4f(4.0 / 14.0 + nelec / 7.0)

        elif elem_id == 61:
            return self._get_occ_4f(5.0 / 14.0 + nelec / 7.0)

        elif elem_id == 62:
            return self._get_occ_4f(6.0 / 14.0 + nelec / 7.0)

        elif elem_id == 63:
            return self._get_occ_4f(7.0 / 14.0 + nelec / 7.0)

        elif elem_id == 64:
            return self._get_occ_4f(8.0 / 14.0 + nelec / 7.0)

        elif elem_id == 65:
            return self._get_occ_4f(9.0 / 14.0 + nelec / 7.0)

        elif elem_id == 66:
            return self._get_occ_4f(10.0 / 14.0 + nelec / 7.0)

        elif elem_id == 67:
            return self._get_occ_4f(11.0 / 14.0 + nelec / 7.0)

        elif elem_id == 68:
            return self._get_occ_4f(12.0 / 14.0 + nelec / 7.0)

        elif elem_id == 69:
            return self._get_occ_4f(13.0 / 14.0 + nelec / 7.0)

        elif elem_id == 70:
            return self._get_occ_4f(14.0 / 14.0 + nelec / 7.0)

        # Lu,Hf,Ta,W,Re,Os,Ir,Pt,Au,Hg

        elif elem_id == 71:
            return self._get_occ_5d(0.1 + nelec / 5.0)

        elif elem_id == 72:
            return self._get_occ_5d(0.2 + nelec / 5.0)

        elif elem_id == 73:
            return self._get_occ_5d(0.3 + nelec / 5.0)

        elif elem_id == 74:
            return self._get_occ_5d(0.4 + nelec / 5.0)

        elif elem_id == 75:
            return self._get_occ_5d(0.5 + nelec / 5.0)

        elif elem_id == 76:
            return self._get_occ_5d(0.6 + nelec / 5.0)

        elif elem_id == 77:
            return self._get_occ_5d(0.7 + nelec / 5.0)

        elif elem_id == 78:
            return self._get_occ_5d(0.8 + nelec / 5.0)

        elif elem_id == 79:
            return self._get_occ_5d(0.9 + nelec / 5.0)

        elif elem_id == 80:
            return self._get_occ_5d(1.0 + nelec / 5.0)

        # Tl,Pb,Bi,Po,At,Rn (81-86)

        elif elem_id == 81:
            return self._get_occ_6s6p(0.375 + nelec / 4.0)

        elif elem_id == 82:
            return self._get_occ_6s6p(0.500 + nelec / 4.0)

        elif elem_id == 83:
            return self._get_occ_6s6p(0.625 + nelec / 4.0)

        elif elem_id == 84:
            return self._get_occ_6s6p(0.750 + nelec / 4.0)

        elif elem_id == 85:
            return self._get_occ_6s6p(0.875 + nelec / 4.0)

        elif elem_id == 86:
            return self._get_occ_6s6p(1.000 + nelec / 4.0)

        # Fr,Ra (87-88)

        elif elem_id == 87:
            return self._get_occ_7s(0.5 + nelec)

        elif elem_id == 88:
            return self._get_occ_7s(1.0 + nelec)

        # Ac,Th,Pa,U,Np,Pu,Am,Cm (89-96)

        elif elem_id == 89:
            return self._get_occ_5f(1.0 / 14.0 + nelec / 7.0)

        elif elem_id == 90:
            return self._get_occ_5f(2.0 / 14.0 + nelec / 7.0)

        elif elem_id == 91:
            return self._get_occ_5f(3.0 / 14.0 + nelec / 7.0)

        elif elem_id == 92:
            return self._get_occ_5f(4.0 / 14.0 + nelec / 7.0)

        elif elem_id == 93:
            return self._get_occ_5f(5.0 / 14.0 + nelec / 7.0)

        elif elem_id == 94:
            return self._get_occ_5f(6.0 / 14.0 + nelec / 7.0)

        elif elem_id == 95:
            return self._get_occ_5f(7.0 / 14.0 + nelec / 7.0)

        elif elem_id == 96:
            return self._get_occ_5f(8.0 / 14.0 + nelec / 7.0)

        # default

        else:
            return []

    def get_alpha_beta_occ_for_molecule(self, molecule, net_charge,
                                        num_unpaired_electrons):

        partial_charges = molecule.get_partial_charges(net_charge)

        elem_ids = molecule.get_element_ids()
        sum_elem_ids = sum(elem_ids)

        natoms = molecule.number_of_atoms()
        use_hint_for_unpaired_electrons = (len(
            self._num_unpaired_electrons_on_atoms) == natoms)

        sum_unpaired_electrons_on_atoms = 0
        if use_hint_for_unpaired_electrons:
            sum_unpaired_electrons_on_atoms = sum(
                self._num_unpaired_electrons_on_atoms)

        alpha_occ_for_atoms = []
        beta_occ_for_atoms = []

        for i in range(natoms):

            weight = elem_ids[i] / sum_elem_ids
            num_excess_elec = (num_unpaired_electrons -
                               sum_unpaired_electrons_on_atoms) * weight

            a_elec = 0.5 * num_excess_elec
            b_elec = -0.5 * num_excess_elec

            if use_hint_for_unpaired_electrons:
                a_elec += 0.5 * self._num_unpaired_electrons_on_atoms[i]
                b_elec -= 0.5 * self._num_unpaired_electrons_on_atoms[i]

            a_occ = self.get_occupation_numbers_for_element(
                elem_ids[i], -partial_charges[i] * 0.5 + a_elec)
            b_occ = self.get_occupation_numbers_for_element(
                elem_ids[i], -partial_charges[i] * 0.5 + b_elec)

            alpha_occ_for_atoms.append(a_occ)
            beta_occ_for_atoms.append(b_occ)

        return alpha_occ_for_atoms, beta_occ_for_atoms

    def get_ao_indices_of_atoms(self, molecule, basis):
        """
        TODO
        """

        natoms = molecule.number_of_atoms()
        aoinds_atoms = [[] for atomidx in range(natoms)]

        max_angl = basis.max_angular_momentum()

        aoidx = 0
        for angl in range(max_angl + 1):
            for s in range(-angl, angl + 1):
                for atomidx in range(natoms):
                    nao = basis.number_of_basis_functions([atomidx], angl)
                    for i in range(nao):
                        aoinds_atoms[atomidx].append(aoidx)
                        aoidx += 1

        return aoinds_atoms

    def compute(self, molecule, basis_1, basis_2, density_type):
        """
        TODO
        """

        ovl_drv = OverlapDriver

        S12 = ovl_drv.compute(molecule, basis_1, basis_2)
        S22 = ovl_drv.compute(molecule, basis_2)

        return self._comp_sad_guess(molecule, basis_1, basis_2, S12, S22,
                                    density_type)

    def _comp_sad_guess(self, molecule, basis_1, basis_2, S12, S22,
                        density_type):
        """
        TODO
        """

        natoms = molecule.number_of_atoms()

        nao_1 = S12.number_of_rows()
        nao_2 = S12.number_of_columns()

        # sanity check

        err_ovl_size = 'SadGuessDriver._comp_sad_guess: '
        err_ovl_size += 'Mismatch between overlap matrices'
        assert_msg_critical(nao_2 == S22.number_of_rows(), err_ovl_size)
        assert_msg_critical(nao_2 == S22.number_of_columns(), err_ovl_size)

        # AO indices for atoms

        aoinds_atoms_1 = self.get_ao_indices_of_atoms(molecule, basis_1)
        aoinds_atoms_2 = self.get_ao_indices_of_atoms(molecule, basis_2)

        # more sanity check

        count_ao_1 = sum([len(x) for x in aoinds_atoms_1])
        count_ao_2 = sum([len(x) for x in aoinds_atoms_2])

        err_bas_size = 'SadGuessDriver._comp_sad_guess: '
        err_bas_size += 'Mismatch between basis set and overlap matrix'
        assert_msg_critical(count_ao_1 == nao_1 and count_ao_2 == nao_2,
                            err_bas_size)

        # alpha and beta occupation numbers for atoms

        net_charge = molecule.get_charge()
        multiplicity = molecule.get_multiplicity()

        alpha_occ, beta_occ = self.get_alpha_beta_occ_for_molecule(
            molecule, net_charge, multiplicity - 1)

        # C_SAD matrix

        for atomidx in range(natoms):

            aoinds_1 = aoinds_atoms_1[atomidx]
            aoinds_2 = aoinds_atoms_2[atomidx]

            naodim_1 = len(aoinds_1)
            naodim_2 = len(aoinds_2)

            err_ao_size = 'SadGuessDriver._comp_sad_guess: '
            err_ao_size += 'Mismatch between basis set and occupation number'
            assert_msg_critical(
                len(alpha_occ[atomidx]) == naodim_1, err_ao_size)
            assert_msg_critical(len(beta_occ[atomidx]) == naodim_1, err_ao_size)

            block_12 = S12[aoinds_1, aoinds_2]
            block_22 = S22[aoinds_2, aoinds_2]

            mat_c1 = np.diag(np.ones(naodim_1))
            mat_a = np.matmul(block_12.T, mat_c1)

            evals, evecs = np.linalg.eigh(block_22)
            block_22_inv = np.linalg.multi_dot(
                [evecs, np.diag(1.0 / evals), evecs.T])

            mat_m = np.linalg.multi_dot([mat_a.T, block_22_inv, mat_a])

            evals, evecs = np.linalg.eigh(mat_m)
            mat_m_invsqrt = np.linalg.multi_dot(
                [evecs, np.diag(1.0 / np.sqrt(evals)), evecs.T])

            mat_c2 = np.linalg.multidot([block_22_inv, mat_a, mat_m_invsqrt])

            if density_type.lower() == 'restricted':

                csad = np.zeros((nao_2, nao_1))

                sqrt_occ = np.sqrt(alpha_occ[atomidx])

                for j in range(naodim_2):
                    c2_j_sqrt_occ = mat_c2[j] * sqrt_occ

                    for i in range(naodim_1):
                        csad[aoinds_2[j], aoinds_1[i]] = c2_j_sqrt_occ[i]

                return np.matmul(csad, csad.T)

            elif density_type.lower() == 'unrestricted':

                csad_a = np.zeros((nao_2, nao_1))
                csad_b = np.zeros((nao_2, nao_1))

                sqrt_a_occ = np.sqrt(alpha_occ[atomidx])
                sqrt_b_occ = np.sqrt(beta_occ[atomidx])

                for j in range(naodim_2):
                    c2_j_sqrt_a_occ = mat_c2[j] * sqrt_a_occ
                    c2_j_sqrt_b_occ = mat_c2[j] * sqrt_b_occ

                    for i in range(naodim_1):
                        csad_a[aoinds_2[j], aoinds_1[i]] = c2_j_sqrt_a_occ[i]
                        csad_b[aoinds_2[j], aoinds_1[i]] = c2_j_sqrt_b_occ[i]

                return (
                    np.matmul(csad_a, csad_a.T),
                    np.matmul(csad_b, csad_b.T),
                )

            else:

                return np.zeros(0)
