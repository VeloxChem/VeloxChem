from mpi4py import MPI
from pathlib import Path
from unittest.mock import patch
import numpy as np
import pytest
import sys

from veloxchem.veloxchemlib import mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.linearsolver import LinearSolver
from veloxchem.lreigensolver import LinearResponseEigenSolver
from veloxchem.tdaeigensolver import TdaEigenSolver
from veloxchem.lrsolver import LinearResponseSolver
from veloxchem.tddftgradientdriver import TddftGradientDriver
from veloxchem.tddfthessiandriver import TddftHessianDriver
from veloxchem.riresponse import lend_ri
from veloxchem.inputparser import get_random_string_parallel
from veloxchem.main import main
from veloxchem.errorhandler import VeloxChemError


def holds_ri(solver):

    return solver._ri_drv is not None or solver._ri_jk_drv is not None


@pytest.mark.solvers
class TestResponseReleaseRi:
    """
    A linear solver lets go of its RI set on release_ri, and the input-file
    path releases the response's set before the CPHF that follows it.
    """

    @staticmethod
    def get_water():

        molecule = Molecule.read_str(
            """
            O 0.0000000000    0.0000000000   -0.0254395383
            H 0.0000000000    0.7695699584    0.5948147012
            H 0.0000000000   -0.7695699584    0.5948147012
            """,
            units='angstrom')
        basis = MolecularBasis.read(molecule, 'def2-svp', ostream=None)

        return molecule, basis

    @staticmethod
    def run_scf(molecule, basis, xcfun):

        scf_drv = ScfRestrictedDriver()
        scf_drv.ri = True
        scf_drv.xcfun = xcfun
        scf_drv.ostream.mute()

        return scf_drv.compute(molecule, basis)

    @staticmethod
    def make_solver(solver_type):

        solver = solver_type()
        solver.ostream.mute()
        if solver_type is LinearResponseSolver:
            solver.a_operator = 'electric dipole'
            solver.b_operator = 'electric dipole'
            solver.frequencies = (0.0,)
        else:
            solver.nstates = 2

        return solver

    @staticmethod
    def get_values(solver, results):

        if solver.rank != mpi_master():
            return None
        if isinstance(solver, LinearResponseSolver):
            return np.array(list(results['response_functions'].values()))
        return results['eigenvalues']

    @pytest.mark.parametrize('solver_type', [
        TdaEigenSolver, LinearResponseEigenSolver, LinearResponseSolver
    ])
    @pytest.mark.parametrize('xcfun', ['blyp', 'b3lyp'])
    def test_release_ri(self, solver_type, xcfun):

        molecule, basis = self.get_water()
        scf_results = self.run_scf(molecule, basis, xcfun)

        solver = self.make_solver(solver_type)
        results = solver.compute(molecule, basis, scf_results)
        values = self.get_values(solver, results)

        assert holds_ri(solver)

        solver.release_ri()

        assert not holds_ri(solver)
        assert not solver._ri_j_fallback
        assert not solver._ri_jk_fallback

        # NOTE: computed again, the solver sets up its RI anew.
        results_again = solver.compute(molecule, basis, scf_results)
        values_again = self.get_values(solver, results_again)

        assert holds_ri(solver)

        if solver.rank == mpi_master():
            assert np.max(np.abs(values_again - values)) < 1.0e-8

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason='pytest.raises only valid in serial')
    def test_lent_and_borrowed_sets_are_refused(self):

        molecule, basis = self.get_water()
        scf_results = self.run_scf(molecule, basis, 'b3lyp')

        driver = self.make_solver(TdaEigenSolver)
        driver.compute(molecule, basis, scf_results)
        solver = self.make_solver(TdaEigenSolver)

        assert driver._ri_jk_drv is not None

        with lend_ri(driver, solver):
            with pytest.raises(VeloxChemError, match='borrowed'):
                solver.release_ri()
            with pytest.raises(VeloxChemError, match='are lent'):
                driver.release_ri()

        assert driver._ri_jk_drv is not None

    def record_releases(self, monkeypatch):

        releases = []
        release_ri = LinearSolver.release_ri

        def recorded_release_ri(solver):
            releases.append((type(solver).__name__, holds_ri(solver)))
            release_ri(solver)

        monkeypatch.setattr(LinearSolver, 'release_ri', recorded_release_ri)

        return releases

    def run_scf_and_tda(self):

        molecule, basis = self.get_water()

        scf_drv = ScfRestrictedDriver()
        scf_drv.ri = True
        scf_drv.xcfun = 'blyp'
        scf_drv.ostream.mute()
        scf_drv.compute(molecule, basis)

        rsp_solver = self.make_solver(TdaEigenSolver)
        rsp_results = rsp_solver.compute(molecule, basis, scf_drv.scf_results)

        return molecule, basis, scf_drv, rsp_solver, rsp_results

    def test_compute_energy_releases_the_response(self, monkeypatch):

        molecule, basis, scf_drv, rsp_solver, rsp_results = (
            self.run_scf_and_tda())

        releases = self.record_releases(monkeypatch)

        # NOTE: an optimization step computes the energy, and with it the
        # response, before the gradient.
        grad_drv = TddftGradientDriver(scf_drv)
        grad_drv.ostream.mute()
        grad_drv.state_deriv_index = [1]
        energy = grad_drv.compute_energy(molecule, basis, scf_drv,
                                         rsp_solver, rsp_results)

        assert releases == [('TdaEigenSolver', True)]
        assert not holds_ri(rsp_solver)

        if scf_drv.rank == mpi_master():
            ref_energy = (scf_drv.scf_results['scf_energy'] +
                          rsp_results['eigenvalues'][0])
            assert abs(energy - ref_energy) < 1.0e-8

    def test_tddft_hessian_gradient_releases_the_response(self, monkeypatch):

        molecule, basis, scf_drv, rsp_solver, _ = self.run_scf_and_tda()

        releases = self.record_releases(monkeypatch)

        hess_drv = TddftHessianDriver(scf_drv, rsp_solver)
        hess_drv.ostream.mute()
        hess_drv.state_deriv_index = [1]
        hess_drv.compute_gradient(molecule, basis)

        # NOTE: the response's set before the orbital response, and the
        # orbital response's after its solve
        # (TddftGradientDriver.compute_analytical).
        assert releases == [('TdaEigenSolver', True),
                            ('TddftOrbitalResponse', True)]

    input_lines = {
        'gradient_task': """
            @jobs
            task: gradient
            @end
            @method settings
            xcfun: blyp
            basis: def2-svp
            @end
            @response
            property: absorption
            nstates: 2
            @end
            @gradient
            state_deriv_index: 1
            @end
        """,
        'response_task': """
            @jobs
            task: response
            @end
            @method settings
            xcfun: blyp
            basis: def2-svp
            @end
            @response
            property: absorption
            nstates: 2
            @end
            @gradient
            state_deriv_index: 1
            @end
        """,
        'polgrad_task': """
            @jobs
            task: polarizability_gradient
            @end
            @method settings
            xcfun: blyp
            basis: def2-svp
            @end
            @polarizability_gradient
            frequencies: 0.0
            @end
        """,
    }

    molecule_lines = """
        @molecule
        charge: 0
        multiplicity: 1
        xyz:
        O 0.0000000000    0.0000000000   -0.0254395383
        H 0.0000000000    0.7695699584    0.5948147012
        H 0.0000000000   -0.7695699584    0.5948147012
        @end
    """

    def run_input(self, monkeypatch, lines):

        releases = self.record_releases(monkeypatch)

        here = Path(__file__).parent
        random_string = get_random_string_parallel()
        input_file = (here / 'data' /
                      f'vlx_rsp_ri_release_{random_string}.inp')

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            with open(input_file, 'w') as f_inp:
                for line in (lines + self.molecule_lines).splitlines():
                    if line.strip():
                        f_inp.write(f'{line.strip()}\n')

        MPI.COMM_WORLD.barrier()

        with patch.object(sys, 'argv', ['vlx', str(input_file)]):
            main()

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            for suffix in ['.inp', '.h5', '_scf.h5', '_rsp.h5', '_orbrsp.h5']:
                path = Path(str(input_file).replace('.inp', suffix))
                if path.is_file():
                    path.unlink()

        MPI.COMM_WORLD.barrier()

        return releases

    @pytest.mark.parametrize('task, expected', [
        ('gradient_task', [('LinearResponseEigenSolver', True),
                           ('TddftOrbitalResponse', True)]),
        ('response_task', [('LinearResponseEigenSolver', True),
                           ('TddftOrbitalResponse', True)]),
        ('polgrad_task', [('LinearResponseSolver', True)]),
    ])
    def test_input_releases_the_response(self, monkeypatch, task, expected):

        releases = self.run_input(monkeypatch, self.input_lines[task])

        # NOTE: True is an RI set held when it was released; the response's
        # comes first, before the CPHF that follows it.
        assert releases == expected
