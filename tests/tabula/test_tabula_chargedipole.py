import numpy as np

from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem import NuclearPotentialGeom010Driver
from veloxchem.tabulalib import (
    TabulaNuclearAttractionDriver,
    TabulaChargeDipoleDriver,
    TabulaDenseMatrix,
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

    def test_screening_sound(self):
        # a stretched H chain: a moderate screen drops far bra-ket blocks, and
        # the screened matrix stays within the (sound, conservative) threshold
        xyz = "\n".join(f"H 0.0 0.0 {i * 2.6:.3f}" for i in range(60))
        mol = Molecule.read_str(xyz, "au")
        bas = MolecularBasis.read(mol, "def2-svp", ostream=None)
        coords = [[0.5, 0.3, 10.0], [-0.4, 0.2, 70.0]]
        moments = [[0.3, -0.4, 0.5], [0.1, 0.2, -0.3]]
        drv = TabulaChargeDipoleDriver()

        exact = drv.compute(mol, bas, moments, coords, 0.0).to_numpy()
        threshold = 1.0e-9
        sparse = drv.compute_sparse(mol, bas, moments, coords, threshold)
        screened = sparse.to_dense().to_numpy()

        # sound: the screened matrix stays within the threshold of exact
        assert np.allclose(screened, exact, 0.0, threshold)
        # and screening actually dropped far blocks
        dim = exact.shape[0]
        assert sparse.stored_element_count() < dim * dim

    def test_auto_exact_for_small_molecule(self):
        # a small molecule (< 100 atoms): the auto default (threshold < 0) is
        # exact dense, bit-for-bit identical to an explicit threshold of 0
        mol, bas = self._setup()
        moments = [m for (m, _) in SITES]
        coords = [r for (_, r) in SITES]
        drv = TabulaChargeDipoleDriver()
        auto = drv.compute(mol, bas, moments, coords).to_numpy()
        exact = drv.compute(mol, bas, moments, coords, 0.0).to_numpy()
        assert np.array_equal(auto, exact)

    def test_auto_screen_for_large_molecule(self):
        # a 100-atom stretched H chain (>= 100 atoms) triggers the auto-screen
        # default; the screened-dense matrix stays within the conservative auto
        # threshold of the exact one, and did drop far blocks
        xyz = "\n".join(f"H 0.0 0.0 {i * 2.6:.3f}" for i in range(100))
        mol = Molecule.read_str(xyz, "au")
        bas = MolecularBasis.read(mol, "def2-svp", ostream=None)
        coords = [[0.5, 0.3, 20.0], [-0.4, 0.2, 200.0]]
        moments = [[0.3, -0.4, 0.5], [0.1, 0.2, -0.3]]
        drv = TabulaChargeDipoleDriver()
        auto = drv.compute(mol, bas, moments, coords).to_numpy()        # default: auto-screen
        exact = drv.compute(mol, bas, moments, coords, 0.0).to_numpy()  # explicit exact dense
        assert np.allclose(auto, exact, 0.0, 1.0e-12)
        assert not np.array_equal(auto, exact)

    @staticmethod
    def _density(n):
        # an arbitrary symmetric "density" (the field contraction is linear in D,
        # so correctness needs no physical density)
        rng = np.random.default_rng(0)
        a = rng.standard_normal((n, n))
        return a + a.T

    def test_field_finite_difference(self):
        # E_i(R) = d/dR_i <D, V(R)>, V(R) = <D, nuclear(unit charge @ R)>
        mol, bas = self._setup()
        nuc = TabulaNuclearAttractionDriver()
        cdp = TabulaChargeDipoleDriver()
        n = nuc.compute(mol, bas, 0.0).to_numpy().shape[0]
        D = self._density(n)
        coords = [r for (_, r) in SITES]

        E = cdp.compute_field(mol, bas, TabulaDenseMatrix.from_numpy(D), coords)
        assert E.shape == (len(coords), 3)

        def V(R):
            return float(np.sum(D * nuc.compute_external(mol, bas, [1.0], [R], 0.0).to_numpy()))

        h = 2.0e-4
        fd = np.zeros_like(E)
        for k, R in enumerate(coords):
            for ax in range(3):
                rp, rm = list(R), list(R)
                rp[ax] += h
                rm[ax] -= h
                fd[k, ax] = (V(rp) - V(rm)) / (2.0 * h)
        assert np.allclose(E, fd, 0.0, 1.0e-6)

    def test_field_vs_veloxchem(self):
        # E_i(R) = <D, geom010 field component i at R>, per point (single-site)
        mol, bas = self._setup()
        cdp = TabulaChargeDipoleDriver()
        n = cdp.compute(mol, bas, [[1.0, 0.0, 0.0]], [[0.0, 0.0, 0.0]], 0.0).to_numpy().shape[0]
        D = self._density(n)
        coords = [r for (_, r) in SITES]

        E = cdp.compute_field(mol, bas, TabulaDenseMatrix.from_numpy(D), coords)

        g010 = NuclearPotentialGeom010Driver()
        ref = np.zeros_like(E)
        for k, R in enumerate(coords):
            mats = g010.compute(mol, bas, [1.0], [R])
            for ax in range(3):
                sm = mats.matrix("XYZ"[ax]).full_matrix()
                ref[k, ax] = np.sum(D * np.array(sm.to_numpy(), copy=True))
        assert np.allclose(E, ref, 0.0, 1.0e-9)

    def test_field_screening_sound(self):
        # a stretched H chain: the density-weighted shell-pair screen drops the
        # bra-ket-distant pairs, and the screened field stays within threshold
        xyz = "\n".join(f"H 0.0 0.0 {i * 2.6:.3f}" for i in range(40))
        mol = Molecule.read_str(xyz, "au")
        bas = MolecularBasis.read(mol, "def2-svp", ostream=None)
        cdp = TabulaChargeDipoleDriver()
        n = cdp.compute(mol, bas, [[1.0, 0.0, 0.0]], [[0.0, 0.0, 0.0]], 0.0).to_numpy().shape[0]
        D = self._density(n)
        Dt = TabulaDenseMatrix.from_numpy(D)
        coords = [[0.5, 0.3, 13.0], [-0.4, 0.2, 50.0], [1.0, 1.0, 90.0]]

        threshold = 1.0e-9
        exact = cdp.compute_field(mol, bas, Dt, coords, 0.0)
        screened = cdp.compute_field(mol, bas, Dt, coords, threshold)

        assert np.allclose(screened, exact, 0.0, threshold)  # sound
        assert not np.array_equal(screened, exact)           # screening dropped shell-pairs

    def test_field_dual_of_matrix(self):
        # the two modes are transposes: E_i(R) = <D, matrix(unit moment axis i @ R)>
        mol, bas = self._setup()
        cdp = TabulaChargeDipoleDriver()
        R = [0.7, -0.5, 1.1]
        n = cdp.compute(mol, bas, [[1.0, 0.0, 0.0]], [R], 0.0).to_numpy().shape[0]
        D = self._density(n)

        E = cdp.compute_field(mol, bas, TabulaDenseMatrix.from_numpy(D), [R])[0]
        for ax in range(3):
            unit = [[1.0 if k == ax else 0.0 for k in range(3)]]
            M = cdp.compute(mol, bas, unit, [R], 0.0).to_numpy()
            assert abs(E[ax] - float(np.sum(D * M))) < 1.0e-11
