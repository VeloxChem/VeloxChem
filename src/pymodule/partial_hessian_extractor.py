import numpy as np


class PartialHessianExtractor:
    """
    Extracts 3x3 partial Hessian blocks from a full 3N x 3N Cartesian
    Hessian matrix for a given pair of atom indices.

    No dependencies on VeloxChem or OpenMM.
    """

    def extract(self, hessian: np.ndarray, atom_i: int,
                atom_j: int) -> np.ndarray:
        """
        Extract the 3x3 cross-derivative block d^2E / (d_coords_i d_coords_j)
        from the full Hessian matrix.

        Parameters
        ----------
        hessian : np.ndarray, shape (3N, 3N)
            Full Cartesian Hessian matrix in any consistent unit.
        atom_i : int
            Zero-based index of the first atom.
        atom_j : int
            Zero-based index of the second atom.

        Returns
        -------
        np.ndarray, shape (3, 3)
        """
        row = 3 * atom_i
        col = 3 * atom_j
        return hessian[row:row + 3, col:col + 3].copy()
