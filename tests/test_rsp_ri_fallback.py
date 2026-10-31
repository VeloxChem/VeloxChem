import h5py
import numpy as np
import pytest
from mpi4py import MPI

from veloxchem.veloxchemlib import mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.scfunrestdriver import ScfUnrestrictedDriver
from veloxchem.lreigensolver import LinearResponseEigenSolver
from veloxchem.lreigensolverunrest import LinearResponseUnrestrictedEigenSolver
from veloxchem.tdaeigensolver import TdaEigenSolver
from veloxchem.cppsolver import ComplexResponseSolver
from veloxchem.quadraticresponsedriver import QuadraticResponseDriver
from veloxchem.resultsio import write_scf_results_to_hdf5, read_results


@pytest.mark.solvers
class TestResponseRiFallback:
    """Response follows an SCF which fell back from RI to the four-center
    integrals for want of memory."""

    @staticmethod
    def get_water():

        molecule = Molecule.read_molecule_string(
            """
            O  0.0           0.0  0.0
            H   .7586020000  0.0  -.5042840000
            H   .7586020000  0.0   .5042840000
            """,
            units='angstrom')
        basis = MolecularBasis.read(molecule, 'def2-svp', ostream=None)

        return molecule, basis

    @staticmethod
    def get_ch3():

        molecule = Molecule.read_molecule_string(
            """
            C   0.0000000   0.0000000   0.0000000
            H   1.0790000   0.0000000   0.0000000
            H  -0.5395000   0.9344500   0.0000000
            H  -0.5395000  -0.9344500   0.0000000
            """,
            units='angstrom')
        molecule.set_multiplicity(2)
        basis = MolecularBasis.read(molecule, 'def2-svp', ostream=None)

        return molecule, basis

    @staticmethod
    def run_scf(molecule, basis, xcfun, ri, budget=None, unrestricted=False):

        driver_type = (ScfUnrestrictedDriver
                       if unrestricted else ScfRestrictedDriver)
        scf_drv = driver_type()
        scf_drv.ostream.mute()
        scf_drv.ri = ri
        scf_drv.xcfun = xcfun
        scf_drv.conv_thresh = 1.0e-8
        scf_drv.ri_memory_budget = budget
        scf_results = scf_drv.compute(molecule, basis)

        return scf_drv, scf_results

    def run_both(self, molecule, basis, xcfun, make_solver, unrestricted=False):
        """Runs the response on an SCF which fell back and on a plain
        four-center SCF, and returns the two solvers and their results."""

        scf_drv, fallen_scf = self.run_scf(molecule,
                                           basis,
                                           xcfun,
                                           ri=True,
                                           budget=1.0e-5,
                                           unrestricted=unrestricted)
        assert scf_drv._ri_fallback

        fallen_drv = make_solver()
        fallen_drv.ostream.mute()
        fallen = fallen_drv.compute(molecule, basis, fallen_scf)

        _, plain_scf = self.run_scf(molecule,
                                    basis,
                                    xcfun,
                                    ri=False,
                                    unrestricted=unrestricted)

        plain_drv = make_solver()
        plain_drv.ostream.mute()
        plain = plain_drv.compute(molecule, basis, plain_scf)

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            assert fallen_scf['ri_fallback']
            assert abs(fallen_scf['scf_energy'] -
                       plain_scf['scf_energy']) < 1.0e-10

        return fallen_drv, fallen, plain_drv, plain

    @staticmethod
    def check_followed(fallen_drv, plain_drv, label):
        """The solver read the fallback, kept what was asked for, and formed
        no RI driver."""

        assert fallen_drv._ri_fallback
        assert not plain_drv._ri_fallback

        if label == 'RI-J':
            assert fallen_drv.ri_coulomb
            assert not fallen_drv._uses_ri_coulomb()
            assert fallen_drv._ri_drv is None
        else:
            assert fallen_drv.ri_jk
            assert not fallen_drv._uses_ri_jk()
            assert fallen_drv._ri_jk_drv is None

    @pytest.mark.parametrize('xcfun, label', [
        ('hf', 'RI-JK'),
        ('b3lyp', 'RI-JK'),
        ('cam-b3lyp', 'RI-JK'),
        ('blyp', 'RI-J'),
    ])
    def test_rpa(self, xcfun, label):

        molecule, basis = self.get_water()

        def make_solver():
            solver = LinearResponseEigenSolver()
            solver.nstates = 3
            solver.conv_thresh = 1.0e-6
            return solver

        fallen_drv, fallen, plain_drv, plain = self.run_both(
            molecule, basis, xcfun, make_solver)

        self.check_followed(fallen_drv, plain_drv, label)

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            assert np.max(
                np.abs(fallen['eigenvalues'] - plain['eigenvalues'])) < 1.0e-8

    @pytest.mark.parametrize('xcfun, label', [
        ('b3lyp', 'RI-JK'),
        ('blyp', 'RI-J'),
    ])
    def test_tda(self, xcfun, label):

        molecule, basis = self.get_water()

        def make_solver():
            solver = TdaEigenSolver()
            solver.nstates = 3
            solver.conv_thresh = 1.0e-6
            return solver

        fallen_drv, fallen, plain_drv, plain = self.run_both(
            molecule, basis, xcfun, make_solver)

        self.check_followed(fallen_drv, plain_drv, label)

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            assert np.max(
                np.abs(fallen['eigenvalues'] - plain['eigenvalues'])) < 1.0e-8

    @pytest.mark.parametrize('xcfun, label', [
        ('b3lyp', 'RI-JK'),
        ('blyp', 'RI-J'),
    ])
    def test_unrestricted_rpa(self, xcfun, label):

        molecule, basis = self.get_ch3()

        def make_solver():
            solver = LinearResponseUnrestrictedEigenSolver()
            solver.nstates = 3
            solver.conv_thresh = 1.0e-6
            return solver

        fallen_drv, fallen, plain_drv, plain = self.run_both(molecule,
                                                             basis,
                                                             xcfun,
                                                             make_solver,
                                                             unrestricted=True)

        self.check_followed(fallen_drv, plain_drv, label)

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            assert np.max(
                np.abs(fallen['eigenvalues'] - plain['eigenvalues'])) < 1.0e-8

    def test_cpp(self):

        molecule, basis = self.get_water()

        def make_solver():
            solver = ComplexResponseSolver()
            solver.frequencies = [0.0, 0.1]
            solver.conv_thresh = 1.0e-6
            return solver

        fallen_drv, fallen, plain_drv, plain = self.run_both(
            molecule, basis, 'b3lyp', make_solver)

        self.check_followed(fallen_drv, plain_drv, 'RI-JK')

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            assert (fallen['response_functions'].keys() ==
                    plain['response_functions'].keys())
            for key, value in plain['response_functions'].items():
                # the response functions are of order 10, and the two SCFs
                # are converged separately to 1e-8
                assert abs(fallen['response_functions'][key] - value) < 1.0e-7

    def test_quadratic_response(self):

        molecule, basis = self.get_water()

        def make_solver():
            solver = QuadraticResponseDriver()
            solver.update_settings({
                'conv_thresh': 1.0e-6,
                'b_frequencies': [0.2],
                'c_frequencies': [0.2],
                'damping': 0.1,
                'a_component': 'z',
                'b_component': 'z',
                'c_component': 'x',
            })
            return solver

        fallen_drv, fallen, plain_drv, plain = self.run_both(
            molecule, basis, 'hf', make_solver)

        self.check_followed(fallen_drv, plain_drv, 'RI-JK')

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            # the response function is of order 30, and the two SCFs are
            # converged separately to 1e-8
            key = ('qrf', 0.2, 0.2)
            assert abs(fallen[key] - plain[key]) < 1.0e-6

    def test_ri_response_without_fallback(self):
        """An SCF with memory to spare does not fall back, and neither does
        its response."""

        molecule, basis = self.get_water()

        scf_drv, scf_results = self.run_scf(molecule, basis, 'b3lyp', ri=True)
        assert not scf_drv._ri_fallback

        solver = LinearResponseEigenSolver()
        solver.ostream.mute()
        solver.nstates = 3
        solver.compute(molecule, basis, scf_results)

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            assert not scf_results['ri_fallback']

        assert not solver._ri_fallback
        assert solver._uses_ri_jk()
        assert solver._ri_jk_drv is not None
        assert solver._ri_jk_drv.is_prepared()

    def test_ri_fallback_is_kept_in_the_checkpoint(self, tmp_path):

        if MPI.COMM_WORLD.Get_rank() != mpi_master():
            return

        for fallback in (True, False):
            fname = str(tmp_path / f'scf_{fallback}.h5')
            with h5py.File(fname, 'w'):
                pass

            write_scf_results_to_hdf5(
                fname, {
                    'ri': True,
                    'ri_coulomb': False,
                    'ri_jk': True,
                    'ri_fallback': fallback,
                })

            results = read_results(fname, 'scf')

            assert bool(results['ri_fallback']) is fallback
