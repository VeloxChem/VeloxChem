import numpy as np

from veloxchem.veloxchemlib import BasisFunction
from veloxchem.veloxchemlib import ScreenedBasisFunctionPair
from veloxchem.veloxchemlib import compute_pa, compute_pb, compute_overlap


class TestObaraSaikaFunc:

    def make_pair(self):

        # bra: 2 primitives, ket: 2 primitives, 3 atom pairs
        bra = BasisFunction([1.0, 2.5], [0.5, 0.7], 0)
        ket = BasisFunction([3.0, 0.8], [0.9, 1.1], 0)

        ax, ay, az = [1.5, 0.0, -1.0], [-0.5, 1.0, 0.3], [0.7, 2.0, 0.0]
        bx, by, bz = [0.3, 1.0, 0.5], [0.2, 1.0, -0.4], [-0.1, 0.0, 1.2]
        atoms = [0, 0, 0]

        pair = ScreenedBasisFunctionPair(bra, 0, ket, 0, ax, ay, az, bx, by, bz,
                                         atoms, atoms)

        # AB = A - B, shape (3, npairs)
        ab = np.array([
            np.array(ax) - np.array(bx),
            np.array(ay) - np.array(by),
            np.array(az) - np.array(bz),
        ])

        return pair, bra, ket, ab

    def test_pa(self):

        pair, bra, ket, ab = self.make_pair()
        be, ke = bra.get_exponents(), ket.get_exponents()
        npb, npk = len(be), len(ke)
        npairs = pair.number_of_pairs()
        block = npb * npk

        ref = np.zeros((3 * block, npairs))
        for c in range(3):
            for i in range(npb):
                for j in range(npk):
                    factor = -ke[j] / (be[i] + ke[j])
                    ref[c * block + i * npk + j, :] = factor * ab[c, :]

        pa = compute_pa(pair).to_numpy()
        assert pa.shape == (3 * block, npairs)
        assert np.allclose(pa, ref, rtol=1e-12, atol=1e-12)

    def test_pb(self):

        pair, bra, ket, ab = self.make_pair()
        be, ke = bra.get_exponents(), ket.get_exponents()
        npb, npk = len(be), len(ke)
        npairs = pair.number_of_pairs()
        block = npb * npk

        ref = np.zeros((3 * block, npairs))
        for c in range(3):
            for i in range(npb):
                for j in range(npk):
                    factor = be[i] / (be[i] + ke[j])
                    ref[c * block + i * npk + j, :] = factor * ab[c, :]

        pb = compute_pb(pair).to_numpy()
        assert pb.shape == (3 * block, npairs)
        assert np.allclose(pb, ref, rtol=1e-12, atol=1e-12)

    def test_overlap(self):

        pair, bra, ket, ab = self.make_pair()
        be, ke = bra.get_exponents(), ket.get_exponents()
        bc, kc = bra.get_normalization_factors(), ket.get_normalization_factors()
        npb, npk = len(be), len(ke)
        npairs = pair.number_of_pairs()

        r2 = ab[0, :]**2 + ab[1, :]**2 + ab[2, :]**2

        ref = np.zeros((npb * npk, npairs))
        for i in range(npb):
            for j in range(npk):
                s = be[i] + ke[j]
                mu = be[i] * ke[j] / s
                pref = bc[i] * kc[j] * (np.pi / s)**1.5
                ref[i * npk + j, :] = pref * np.exp(-mu * r2)

        ovl = compute_overlap(pair).to_numpy()
        assert ovl.shape == (npb * npk, npairs)
        assert np.allclose(ovl, ref, rtol=1e-12, atol=1e-12)
