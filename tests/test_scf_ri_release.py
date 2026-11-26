from mpi4py import MPI
from pathlib import Path
from unittest.mock import patch
import numpy as np
import pytest
import sys

from veloxchem.veloxchemlib import mpi_master
from veloxchem.veloxchemlib import SimdRIJKFockDriver
from veloxchem.errorhandler import VeloxChemError
from veloxchem.inputparser import get_random_string_parallel
from veloxchem.main import main
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.outputstream import OutputStream
from veloxchem.scfdriver import ScfDriver
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.scfunrestdriver import ScfUnrestrictedDriver
from veloxchem.scfgradientdriver import ScfGradientDriver


@pytest.mark.solvers
class TestScfRiRelease:
    """The SCF releases the storage its RI-JK driver keeps between Fock builds
    when it ends, and release_ri releases the whole RI driver."""

    comm = MPI.COMM_WORLD

    @staticmethod
    def get_water():

        molecule = Molecule.read_str(
            """O  0.000  0.000  0.117
               H  0.000  0.757 -0.467
               H  0.000 -0.757 -0.467""", 'angstrom')
        basis = MolecularBasis.read(molecule, 'def2-svp', ostream=None)

        return molecule, basis

    @staticmethod
    def get_ch3():

        molecule = Molecule.read_str(
            """C   0.0000000   0.0000000   0.0000000
               H   1.0790000   0.0000000   0.0000000
               H  -0.5395000   0.9344500   0.0000000
               H  -0.5395000  -0.9344500   0.0000000""", 'angstrom')
        molecule.set_multiplicity(2)
        basis = MolecularBasis.read(molecule, 'def2-svp', ostream=None)

        return molecule, basis

    @staticmethod
    def run_scf(molecule, basis, unrestricted=False):

        driver_type = (ScfUnrestrictedDriver
                       if unrestricted else ScfRestrictedDriver)
        scf_drv = driver_type(ostream=OutputStream(None))
        scf_drv.ri = True
        scf_drv.xcfun = 'b3lyp'
        scf_drv.conv_thresh = 1.0e-8
        scf_results = scf_drv.compute(molecule, basis)

        return scf_drv, scf_results

    def densities(self, scf_drv, scf_results):
        """The converged alpha and beta densities, on every rank."""

        if self.comm.Get_rank() == mpi_master():
            densities = (scf_results['D_alpha'], scf_results['D_beta'])
        else:
            densities = None

        return self.comm.bcast(densities, root=mpi_master())

    def test_compute_releases_the_work_storage(self, monkeypatch):

        calls = []
        real = ScfDriver._release_ri_work_storage

        def recorded(scf_drv):
            calls.append(type(scf_drv._ri_drv).__name__)
            real(scf_drv)

        monkeypatch.setattr(ScfDriver, '_release_ri_work_storage', recorded)

        molecule, basis = self.get_water()
        scf_drv, _ = self.run_scf(molecule, basis)

        assert calls == ['SimdRIJKFockDriver']
        assert scf_drv._ri_drv.is_prepared()

    def test_fock_after_released_work_storage(self):

        molecule, basis = self.get_water()
        scf_drv, scf_results = self.run_scf(molecule, basis)

        assert isinstance(scf_drv._ri_drv, SimdRIJKFockDriver)

        den_a, _ = self.densities(scf_drv, scf_results)

        # NOTE: each rank compares its own share of the Fock matrix.
        kept = scf_drv._simd_ri_jk_fock(den_a, 0.2, use_density_factor=True)

        scf_drv._ri_drv.release_work_storage()

        formed_again = scf_drv._simd_ri_jk_fock(den_a,
                                                0.2,
                                                use_density_factor=True)

        assert np.max(np.abs(kept - formed_again)) < 1.0e-12

    def test_unrestricted_fock_after_released_work_storage(self):

        molecule, basis = self.get_ch3()
        scf_drv, scf_results = self.run_scf(molecule,
                                            basis,
                                            unrestricted=True)

        assert isinstance(scf_drv._ri_drv, SimdRIJKFockDriver)

        den_a, den_b = self.densities(scf_drv, scf_results)

        kept = scf_drv._simd_ri_jk_fock_unrestricted(den_a,
                                                     den_b,
                                                     0.2,
                                                     use_density_factor=True)

        scf_drv._ri_drv.release_work_storage()

        formed_again = scf_drv._simd_ri_jk_fock_unrestricted(
            den_a, den_b, 0.2, use_density_factor=True)

        for a, b in zip(kept, formed_again):
            assert np.max(np.abs(a - b)) < 1.0e-12

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason='the error aborts the run on more ranks')
    def test_gradient_after_release_ri_stops(self):

        molecule, basis = self.get_water()
        scf_drv, scf_results = self.run_scf(molecule, basis)

        scf_drv.release_ri()

        assert scf_drv._ri_drv is None

        grad_drv = ScfGradientDriver(scf_drv)
        grad_drv.ostream.mute()

        with pytest.raises(VeloxChemError, match='release_ri'):
            grad_drv.compute(molecule, basis, scf_results)

    def test_scf_after_release_ri_forms_a_new_driver(self):

        molecule, basis = self.get_water()
        scf_drv, scf_results = self.run_scf(molecule, basis)

        grad_drv = ScfGradientDriver(scf_drv)
        grad_drv.ostream.mute()
        grad_drv.compute(molecule, basis, scf_results)
        before = grad_drv.get_gradient()

        scf_drv.release_ri()

        scf_results = scf_drv.compute(molecule, basis)

        assert not scf_drv._ri_released
        assert scf_drv._ri_drv.is_prepared()

        grad_drv = ScfGradientDriver(scf_drv)
        grad_drv.ostream.mute()
        grad_drv.compute(molecule, basis, scf_results)
        after = grad_drv.get_gradient()

        if self.comm.Get_rank() == mpi_master():
            assert np.max(np.abs(before - after)) < 1.0e-7


@pytest.mark.solvers
class TestInputResponseReleasesScfRi:
    """A response task from an input file releases the SCF's RI driver unless
    a section that is given the SCF driver follows."""

    input_lines = """
        @jobs
        task: response
        @end

        @method settings
        xcfun: b3lyp
        basis: def2-svp
        ri: yes
        @end

        @scf
        restart: no
        @end

        @response
        restart: no
        property: polarizability
        frequencies: 0.0
        @end

        @molecule
        charge: 0
        multiplicity: 1
        units: au
        xyz:
        O   0.0   0.0   0.0
        H   0.0   1.4   1.1
        H   0.0  -1.4   1.1
        @end
    """

    def run_input(self, monkeypatch, lines):

        calls = []
        real = ScfDriver.release_ri

        def recorded(scf_drv):
            calls.append(scf_drv._ri_drv is not None)
            real(scf_drv)

        monkeypatch.setattr(ScfDriver, 'release_ri', recorded)

        here = Path(__file__).parent
        random_string = get_random_string_parallel()
        input_file = here / 'data' / f'vlx_ri_release_{random_string}.inp'

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            with open(input_file, 'w') as f_inp:
                for line in lines.splitlines():
                    if line.strip():
                        f_inp.write(f'{line.strip()}\n')

        MPI.COMM_WORLD.barrier()

        with patch.object(sys, 'argv', ['vlx', str(input_file)]):
            main()

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            for suffix in ['.inp', '.h5', '_scf.h5', '_rsp.h5']:
                path = Path(str(input_file).replace('.inp', suffix))
                if path.is_file():
                    path.unlink()

        MPI.COMM_WORLD.barrier()

        return calls

    def test_response_releases(self, monkeypatch, capsys):

        calls = self.run_input(monkeypatch, self.input_lines)

        # NOTE: True is an RI driver held when it was released.
        assert calls == [True]

        captured = capsys.readouterr()

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            assert 'Using the SIMD resolution' in captured.out
            assert 'Polarizability' in captured.out

    @pytest.mark.parametrize('section',
                             ['gradient', 'vibrational', 'optimize_excited_state'])
    def test_a_following_section_keeps(self, monkeypatch, capsys, section):

        # NOTE: the property is the polarizability, for which none of these
        # sections runs anything; their presence alone keeps the driver. A
        # section needs a line, as the input parser rejects an empty one.
        lines = self.input_lines.replace(
            '@molecule',
            f'@{section.replace("_", " ")}\nnumerical: no\n@end\n\n@molecule')

        calls = self.run_input(monkeypatch, lines)

        assert calls == []

        captured = capsys.readouterr()

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            assert 'Polarizability' in captured.out
