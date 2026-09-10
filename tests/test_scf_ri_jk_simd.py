import pytest

from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.outputstream import OutputStream
from veloxchem.scfrestdriver import ScfRestrictedDriver

# NOTE: the SIMD RI-JK driver is an alternative path through the closed shell Fock
# build, selected by ri_jk_simd. It is the same approximation as the conventional
# RI-JK path, with the same auxiliary basis and the same metric threshold, so the
# two must converge to the same energy. That one number exercises the whole chain,
# the inverted metric, the B vectors, the Y vector, the Coulomb matrix, the W
# matrices and the exchange, through a calculation rather than through constructed
# input.


class TestScfRiJkSimd:

    @pytest.fixture
    def molecule(self):

        return Molecule.read_str(
            """O  0.000  0.000  0.117
               H  0.000  0.757 -0.467
               H  0.000 -0.757 -0.467""", 'angstrom')

    @pytest.fixture
    def basis(self, molecule):

        return MolecularBasis.read(molecule, 'def2-svp', ostream=None)

    def run_scf(self, molecule, basis, **settings):

        driver = ScfRestrictedDriver(ostream=OutputStream(None))

        driver.conv_thresh = 1.0e-8

        for key, value in settings.items():
            setattr(driver, key, value)

        return driver.compute(molecule, basis)['scf_energy']

    @pytest.mark.parametrize('xcfun', [None, 'PBE0', 'B3LYP'])
    def test_simd_matches_the_conventional_driver(self, molecule, basis, xcfun):
        """Hartree-Fock scales the exchange by one and the hybrids by a fraction of
        it, so the functionals cover the scaling as well as the path."""

        settings = {
            'ri_jk': True,
            'ri_auxiliary_basis': 'def2-universal-jkfit',
        }

        if xcfun is not None:
            settings['xcfun'] = xcfun

        conventional = self.run_scf(molecule, basis, **settings)

        simd = self.run_scf(molecule, basis, ri_jk_simd=True, **settings)

        assert abs(simd - conventional) < 1.0e-10, (
            f"{xcfun or 'HF'}: {simd:.12f} against {conventional:.12f}")

    def test_the_approximation_is_the_one_expected(self, molecule, basis):
        """The energy has to differ from the one without the approximation, or the
        test above would pass with the path never taken."""

        plain = self.run_scf(molecule, basis)

        simd = self.run_scf(molecule,
                            basis,
                            ri_jk=True,
                            ri_jk_simd=True,
                            ri_auxiliary_basis='def2-universal-jkfit')

        assert abs(simd - plain) > 1.0e-6
        assert abs(simd - plain) < 1.0e-3

    def test_the_flag_is_off_by_default(self):

        driver = ScfRestrictedDriver(ostream=OutputStream(None))

        assert not driver.ri_jk_simd
        assert driver.ri_memory_budget is None
