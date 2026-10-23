import numpy as np
import pytest

from veloxchem.veloxchemlib import SimdRIJFockDriver
from veloxchem.veloxchemlib import TwoCenterElectronRepulsionDriver
from veloxchem.veloxchemlib import make_matrix, mat_t, rimode
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.packedmatrix import PackedMatrix

# The Coulomb fitting of both RI-J drivers is checked against an exact solve of the
# metric, E_J = 1/2 gamma^T M^-1 gamma with gamma(P) = sum (P|mn) D(mn). The drivers
# solve with a Cholesky factor of the metric; an inverse of it formed outright loses
# the precision the condition of the metric costs, which on larger molecules is 1e-7
# to 1e-6 Hartree of the Coulomb energy.


@pytest.mark.solvers
class TestRIJFitting:

    def setup_system(self):

        xyz_string = """6
        xyz
        O   -0.1858140  -1.1749469   0.7662596
        H   -0.1285513  -0.8984365   1.6808606
        H   -0.0582782  -0.3702550   0.2638279
        O    0.1747051   2.1050002  -0.7244430
        H   -0.5650842   2.3134964  -1.2949455
        H    0.9282185   2.0652990  -1.3134026
        """

        mol = Molecule.read_xyz_string(xyz_string)
        bas = MolecularBasis.read(mol, 'def2-svp', ostream=None)
        aux = MolecularBasis.read(mol, 'def2-universal-jfit', ostream=None)

        nao = bas.get_dimensions_of_basis()

        rng = np.random.default_rng(7)
        orbitals = rng.standard_normal((nao, 10)) * 0.3
        density = orbitals @ orbitals.T

        return mol, bas, aux, density

    def exact_energy(self, mol, aux, gamma):

        metric = TwoCenterElectronRepulsionDriver().compute(mol, aux).to_numpy()

        return 0.5 * gamma @ np.linalg.solve(metric, gamma)

    def test_simd_driver(self):

        mol, bas, aux, density = self.setup_system()
        nao = density.shape[0]

        drv = SimdRIJFockDriver()
        metric = drv.make_metric(mol, aux, 1.0e-12)
        assert metric.get_type() == mat_t.lower_triangular

        drv.prepare(mol, bas, aux, 1.0e-12, 1024**3, 1.0e-12,
                    rimode.in_memory, metric, 0, 1)

        packed = PackedMatrix(nao, nao, mat_t.symmetric)
        packed.from_numpy(np.ascontiguousarray(density))

        gamma = np.array(drv.compute_gamma(packed, drv.owned_parts()))
        fitted = np.array(drv.solve_fitting(gamma))

        fock = PackedMatrix(nao, nao, mat_t.symmetric)
        fock.zero()
        drv.compute_coulomb(list(fitted), drv.owned_parts(), fock)

        energy = 0.5 * np.sum(density * fock.to_numpy())

        assert abs(energy - self.exact_energy(mol, aux, gamma)) < 1.0e-10

        # an inverted metric (symmetric) handed in is still multiplied by
        inverse = PackedMatrix(metric.number_of_rows(),
                               metric.number_of_columns(), mat_t.symmetric)
        full = TwoCenterElectronRepulsionDriver().compute(mol, aux).to_numpy()
        inverse.from_numpy(np.ascontiguousarray(np.linalg.inv(full)))

        drv_inv = SimdRIJFockDriver()
        drv_inv.prepare(mol, bas, aux, 1.0e-12, 1024**3, 1.0e-12,
                        rimode.in_memory, inverse, 0, 1)

        fitted_inv = np.array(drv_inv.solve_fitting(gamma))

        assert np.max(np.abs(fitted_inv - fitted)) < 1.0e-8
