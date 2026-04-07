import numpy as np


class ForceConstantSolver:
    """
    Stateless solver for PHF force constants.

    Implements the analytical least-squares equations (Eq. 13 / 19 / 23) from:
        Wang, Ozhgibesov & Hirao, J. Comput. Chem. 2016, 37, 2349-2359.

    All three solve methods share the same formula:

        k = sum_{p,s} (H_QM - H0)_{ps} * h_unit_{ps}
            / sum_{p,s} h_unit_{ps}^2

    They differ only in physical interpretation and which H0 is passed in.
    Units of the returned k match the units of h_qm and h0 divided by
    the units of h_unit squared — the caller is responsible for tracking units.
    """

    def _solve(
        self,
        h_qm: np.ndarray,
        h0: np.ndarray,
        h_unit: np.ndarray,
    ) -> float:
        """
        Core least-squares solve shared by all three coordinate types.

        Parameters
        ----------
        h_qm : np.ndarray, shape (3, 3)
            QM partial Hessian block for the terminal atom pair.
        h0 : np.ndarray, shape (3, 3)
            Residual MM partial Hessian: contributions from all terms
            except the one being fitted, evaluated with that term set to 0.
        h_unit : np.ndarray, shape (3, 3)
            Normalised MM partial Hessian: the target term evaluated with
            its force constant set to 1.0 and all nonbonded terms zeroed.

        Returns
        -------
        float
            Fitted force constant k.
        """
        numerator = np.sum((h_qm - h0) * h_unit)
        denominator = np.sum(h_unit**2)

        if abs(denominator) < 1e-30:
            raise ValueError(
                "Denominator in PHF solve is effectively zero. "
                "The h_unit block is all-zero for this internal coordinate, "
                "which means the geometry has this term contributing nothing.")

        return float(numerator / denominator)

    def solve_dihedral(
        self,
        h_qm: np.ndarray,
        h0: np.ndarray,
        h_unit: np.ndarray,
    ) -> float:
        """
        Solve for dihedral force constant k_d (Eq. 13).

        Terminal atom pair is (i, l) for dihedral i-j-k-l.
        h0 contains only the nonbonded 1-4 contribution between i and l.
        """
        return self._solve(h_qm, h0, h_unit)

    def solve_angle(
        self,
        h_qm: np.ndarray,
        h0: np.ndarray,
        h_unit: np.ndarray,
    ) -> float:
        """
        Solve for angle-bending force constant k_a (Eq. 19).

        Terminal atom pair is (i, k) for angle i-j-k.
        h0 contains dihedral contributions involving both i and k,
        computed using k_d values already determined in Stage 1.
        Note: no nonbonded contribution because AMBER has zero 1-3 interactions.
        """
        return self._solve(h_qm, h0, h_unit)

    def solve_improper(
        self,
        h_qm: np.ndarray,
        h0: np.ndarray,
        h_unit: np.ndarray,
    ) -> float:
        """
        Solve for improper-dihedral force constant k_imp.

        Terminal atom pair is (i, l) for improper i-j-k-l.
        h0 includes angle and dihedral contributions involving both i and l,
        using k_a and k_d values already determined in Stages 1 and 2.
        No nonbonded contribution: impropers don't have 1-4 interactions.
        """
        return self._solve(h_qm, h0, h_unit)

    def solve_bond(
        self,
        h_qm: np.ndarray,
        h0: np.ndarray,
        h_unit: np.ndarray,
    ) -> float:
        """
        Solve for bond-stretching force constant k_b (Eq. 23).

        Terminal atom pair is (i, j) for bond i-j.
        h0 contains angle and dihedral contributions involving both i and j,
        using k_a and k_d values already determined in Stages 1 and 2.
        """
        return self._solve(h_qm, h0, h_unit)
