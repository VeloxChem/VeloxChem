import pytest

from veloxchem.veloxchemlib import mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.outputstream import OutputStream
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.scfgradientdriver import ScfGradientDriver
from veloxchem.scfhessiandriver import ScfHessianDriver


@pytest.mark.solvers
class TestScfSetupHeader:
    """
    Verifies that a newly constructed SCF gradient/Hessian driver with
    unset own DFT state inherits the exchange-correlation functional and
    grid level from existing SCF results, and that the setup header
    printed by compute() reflects the inherited settings.
    """

    def make_molecule(self):

        return Molecule.from_xyz_string("""3
            water
            O 0.0 0.0 0.1173
            H 0.0 0.7572 -0.4692
            H 0.0 -0.7572 -0.4692
            """)

    def run_scf(self, xcfun_label, molecule, basis, grid_level=None):

        scf_drv = ScfRestrictedDriver()
        scf_drv.ostream.mute()
        if xcfun_label is not None:
            scf_drv.xcfun = xcfun_label
            scf_drv.grid_level = grid_level
        scf_drv.compute(molecule, basis)
        return scf_drv

    def compute_with_capture(self, drv, molecule, basis, tmp_path):

        if drv.rank == mpi_master():
            outfile = tmp_path / 'output.txt'
            drv.ostream = OutputStream(outfile)
        else:
            drv.ostream = OutputStream(None)

        drv.compute(molecule, basis)

        if drv.rank == mpi_master():
            return outfile.read_text()
        return ''

    def dft_header_lines(self, text):

        return [
            line.strip() for line in text.splitlines()
            if 'Exchange-Correlation Functional' in line
            or 'Molecular Grid Level' in line
        ]

    def header_value_lines(self, lines, label, value):

        # match on the value after the colon so that the assertions are
        # insensitive to the fixed-width formatting of the header lines
        return [
            line for line in lines
            if line.startswith(label)
            and line.split(':')[-1].strip() == str(value)
        ]

    def hessian_driver_section(self, text):

        if 'SCF Hessian Driver Setup' in text:
            # the CPHF solver prints its own setup header afterwards
            return text.split('SCF Hessian Driver Setup')[1].split(
                'Coupled-Perturbed')[0]
        return text

    def test_gradient_header_inherits_dft_from_scf_results(self, tmp_path):

        molecule = self.make_molecule()
        basis = MolecularBasis.read(molecule, 'sto-3g')

        scf_drv = self.run_scf('b3lyp', molecule, basis, grid_level=5)

        # newly constructed driver without method settings
        grad_drv = ScfGradientDriver(scf_drv)
        assert grad_drv.xcfun is None
        assert grad_drv._dft is False

        text = self.compute_with_capture(grad_drv, molecule, basis, tmp_path)

        # DFT state inherited from the SCF results by the sanity checks
        assert grad_drv.xcfun is not None
        assert grad_drv._dft is True
        assert grad_drv.grid_level == 5

        if grad_drv.rank == mpi_master():
            # the setup header is written by the master rank only
            lines = self.dft_header_lines(text)
            assert self.header_value_lines(
                lines, 'Exchange-Correlation Functional', 'B3LYP')
            assert self.header_value_lines(lines, 'Molecular Grid Level', 5)

    def test_gradient_header_with_hf_scf_results(self, tmp_path):

        molecule = self.make_molecule()
        basis = MolecularBasis.read(molecule, 'sto-3g')

        scf_drv = self.run_scf(None, molecule, basis)

        grad_drv = ScfGradientDriver(scf_drv)
        assert grad_drv.xcfun is None
        assert grad_drv._dft is False

        text = self.compute_with_capture(grad_drv, molecule, basis, tmp_path)

        assert grad_drv._dft is False

        if grad_drv.rank == mpi_master():
            assert self.dft_header_lines(text) == []

    def test_hessian_header_inherits_dft_from_scf_results(self, tmp_path):

        molecule = self.make_molecule()
        basis = MolecularBasis.read(molecule, 'sto-3g')

        scf_drv = self.run_scf('b3lyp', molecule, basis, grid_level=5)

        # newly constructed driver without method settings
        hess_drv = ScfHessianDriver(scf_drv)
        assert hess_drv.xcfun is None
        assert hess_drv._dft is False

        text = self.compute_with_capture(hess_drv, molecule, basis, tmp_path)

        # DFT state inherited from the SCF results by the sanity checks
        assert hess_drv.xcfun is not None
        assert hess_drv._dft is True
        assert hess_drv.grid_level == 5

        if hess_drv.rank == mpi_master():
            # the setup header is written by the master rank only
            lines = self.dft_header_lines(self.hessian_driver_section(text))
            assert self.header_value_lines(
                lines, 'Exchange-Correlation Functional', 'B3LYP')
            assert self.header_value_lines(lines, 'Molecular Grid Level', 5)

    def test_hessian_header_with_hf_scf_results(self, tmp_path):

        molecule = self.make_molecule()
        basis = MolecularBasis.read(molecule, 'sto-3g')

        scf_drv = self.run_scf(None, molecule, basis)

        hess_drv = ScfHessianDriver(scf_drv)
        assert hess_drv.xcfun is None
        assert hess_drv._dft is False

        text = self.compute_with_capture(hess_drv, molecule, basis, tmp_path)

        assert hess_drv._dft is False

        if hess_drv.rank == mpi_master():
            assert self.dft_header_lines(self.hessian_driver_section(text)) == []

    @pytest.mark.parametrize('driver_cls',
                             [ScfGradientDriver, ScfHessianDriver])
    def test_header_runs_scf_on_demand_before_inheriting_settings(
            self, driver_cls, tmp_path):

        molecule = self.make_molecule()
        basis = MolecularBasis.read(molecule, 'sto-3g')

        scf_drv = ScfRestrictedDriver()
        scf_drv.ostream.mute()
        scf_drv.xcfun = 'b3lyp'
        scf_drv.grid_level = 5
        assert scf_drv.scf_results is None

        drv = driver_cls(scf_drv)
        text = self.compute_with_capture(drv, molecule, basis, tmp_path)

        assert scf_drv.scf_results is not None
        assert drv.xcfun is not None
        assert drv._dft is True
        assert drv.grid_level == 5

        if drv.rank == mpi_master():
            if isinstance(drv, ScfHessianDriver):
                text = self.hessian_driver_section(text)
            lines = self.dft_header_lines(text)
            assert self.header_value_lines(
                lines, 'Exchange-Correlation Functional', 'B3LYP')
            assert self.header_value_lines(lines, 'Molecular Grid Level', 5)
