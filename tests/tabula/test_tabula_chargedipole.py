import numpy as np

from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem import NuclearPotentialGeom010Driver
from veloxchem.tabulalib import (
    TabulaNuclearAttractionDriver,
    TabulaChargeDipoleDriver,
)


WATER = (
    "O 0.000  0.000 -1.000\n"
    "H 0.000  1.400 -2.100\n"
    "H 0.000 -1.400 -2.100"
)

# off-atom dipole sites with generic moments (au): (moment, position)
SITES = [
    ([0.30, -0.40, 0.55], [1.5, 1.0, 0.5]),
    ([-0.70, 0.20, 0.10], [-1.0, 2.0, -0.5]),
    ([0.10, 0.05, -0.90], [0.3, -1.4, 1.2]),
]


class TestTabulaChargeDipole:
    """Validates ChargeDipoleDriver, the matrix Sum_N d_N . (a|(r-N)/|r-N|^3|c).

    Since (r-N)/|r-N|^3 = grad_N (1/|r-N|), the charge-dipole matrix is
    Sum_N d_N . grad_R V_R with V_R the unit-charge nuclear potential at R, so
    we finite-difference Tabula's own (VeloxChem-validated) nuclear attraction
    over all three axes. We also cross-check the X field component against
    VeloxChem's NuclearPotentialGeom010Driver (its CMatrices Python extraction
    is only reliable for the first key 'X'; the all-axis FD covers Y and Z)."""

    def _setup(self, basis_label="def2-svp"):
        mol = Molecule.read_str(WATER, "au")
        bas = MolecularBasis.read(mol, basis_label, ostream=None)
        return mol, bas

    @staticmethod
    def _vlx_field_axis(mol, bas, R, axis):
        # single-site nuclear-potential gradient (unit charge at R), field
        # component `axis`. A single-site call is used per site because the
        # multi-site geom010 path mis-aggregates the Y/Z matrices, while the
        # single-site matrices extract cleanly for all axes.
        mats = NuclearPotentialGeom010Driver().compute(mol, bas, [1.0], [R])
        cm = mats.matrix("XYZ"[axis])  # keep CMatrix/SubMatrix alive across to_numpy
        sm = cm.full_matrix()
        return np.array(sm.to_numpy(), copy=True)

    def test_finite_difference_vs_nuclear(self):
        # charge-dipole == d . grad_R V_R, finite-differenced over all 3 axes
        mol, bas = self._setup()
        nuc = TabulaNuclearAttractionDriver()
        cdp = TabulaChargeDipoleDriver()

        def V(R):
            return nuc.compute_external(mol, bas, [1.0], [R], 0.0).to_numpy()

        def fd(moment, R, h=2.0e-4):
            g = np.zeros_like(V(R))
            for ax in range(3):
                rp, rm = list(R), list(R)
                rp[ax] += h
                rm[ax] -= h
                g += moment[ax] * (V(rp) - V(rm)) / (2.0 * h)
            return g

        moments = [m for (m, _) in SITES]
        coords = [r for (_, r) in SITES]

        tab = cdp.compute(mol, bas, moments, coords, 0.0).to_numpy()
        ref = sum(fd(m, r) for (m, r) in SITES)

        # FD truncation at h=2e-4 sits ~1e-8; the wrong sign would be ~0.85
        assert np.allclose(tab, ref, 0.0, 1.0e-6)

    def test_vs_veloxchem_geom010(self):
        # independent analytic reference: Sum_N d_N . field_N, built from
        # VeloxChem's NuclearPotentialGeom010Driver one site at a time
        mol, bas = self._setup()
        moments = [m for (m, _) in SITES]
        coords = [r for (_, r) in SITES]

        tab = TabulaChargeDipoleDriver().compute(mol, bas, moments, coords, 0.0).to_numpy()

        ref = np.zeros_like(tab)
        for moment, R in zip(moments, coords):
            for ax in range(3):
                ref += moment[ax] * self._vlx_field_axis(mol, bas, R, ax)

        assert np.allclose(tab, ref, 0.0, 1.0e-10)

    def test_symmetric(self):
        # the operator is a multiplicative function of r, so the matrix is symmetric
        mol, bas = self._setup()
        moments = [m for (m, _) in SITES]
        coords = [r for (_, r) in SITES]
        tab = TabulaChargeDipoleDriver().compute(mol, bas, moments, coords, 0.0).to_numpy()
        assert np.array_equal(tab, tab.T)

    def test_sparse_matches_dense(self):
        # block-sparse storage reproduces the dense matrix
        mol, bas = self._setup()
        moments = [m for (m, _) in SITES]
        coords = [r for (_, r) in SITES]
        drv = TabulaChargeDipoleDriver()
        dense = drv.compute(mol, bas, moments, coords, 0.0).to_numpy()
        sparse = drv.compute_sparse(mol, bas, moments, coords, 1.0e-12).to_dense().to_numpy()
        assert np.allclose(dense, sparse, 0.0, 1.0e-12)
