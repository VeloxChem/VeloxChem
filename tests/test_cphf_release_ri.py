from mpi4py import MPI
import pytest

from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.scfunrestdriver import ScfUnrestrictedDriver
from veloxchem.tdaeigensolver import TdaEigenSolver
from veloxchem.cphfsolver import CphfSolver
from veloxchem.tddftorbitalresponse import TddftOrbitalResponse
from veloxchem.tddftgradientdriver import TddftGradientDriver
from veloxchem.scfhessiandriver import ScfHessianDriver
from veloxchem.errorhandler import VeloxChemError


def get_water(charge=0, multiplicity=1):

    molecule = Molecule.read_str(
        """
        O 0.0000000000    0.0000000000   -0.0254395383
        H 0.0000000000    0.7695699584    0.5948147012
        H 0.0000000000   -0.7695699584    0.5948147012
        """,
        units='angstrom')
    molecule.set_charge(charge)
    molecule.set_multiplicity(multiplicity)
    basis = MolecularBasis.read(molecule, 'def2-svp', ostream=None)

    return molecule, basis


def run_scf(scf_type, molecule, basis, xcfun):

    scf_drv = scf_type()
    scf_drv.ri = True
    scf_drv.xcfun = xcfun
    scf_drv.ostream.mute()
    scf_drv.compute(molecule, basis)

    return scf_drv


def run_tda(molecule, basis, scf_drv):

    rsp_solver = TdaEigenSolver()
    rsp_solver.nstates = 2
    rsp_solver.ostream.mute()
    rsp_results = rsp_solver.compute(molecule, basis, scf_drv.scf_results)

    return rsp_solver, rsp_results


def holds_ri(solver):

    return solver._ri_drv is not None or solver._ri_jk_drv is not None


class TestCphfReleaseRi:
    """
    A CPHF solver lets go of its RI set on release_ri, and the drivers that
    run one call it once its results are taken.
    """

    @pytest.mark.parametrize('xcfun, uses_ri', [
        ('blyp', '_uses_ri_coulomb'),
        ('b3lyp', '_uses_ri_jk'),
    ])
    def test_release_ri(self, xcfun, uses_ri):

        molecule, basis = get_water()
        scf_drv = run_scf(ScfRestrictedDriver, molecule, basis, xcfun)
        _, rsp_results = run_tda(molecule, basis, scf_drv)

        orbrsp = TddftOrbitalResponse()
        orbrsp.update_settings({
            'tamm_dancoff': 'yes',
            'state_deriv_index': '1'
        })
        orbrsp.ostream.mute()
        orbrsp.compute(molecule, basis, scf_drv.scf_results, rsp_results)
        orbrsp.compute_omega(molecule, basis, scf_drv.scf_results)

        assert getattr(orbrsp, uses_ri)()
        assert holds_ri(orbrsp)

        orbrsp.release_ri()

        assert not holds_ri(orbrsp)
        assert not orbrsp._ri_j_fallback
        assert not orbrsp._ri_jk_fallback

        if MPI.COMM_WORLD.Get_size() == 1:
            # NOTE: pytest.raises is only valid in serial
            with pytest.raises(VeloxChemError, match='no RI-J'):
                orbrsp.compute_omega(molecule, basis, scf_drv.scf_results)

    @pytest.mark.parametrize('xcfun, fallback', [
        ('blyp', '_ri_j_fallback'),
        ('b3lyp', '_ri_jk_fallback'),
    ])
    def test_release_ri_after_fallback(self, xcfun, fallback):

        molecule, basis = get_water()
        scf_drv = run_scf(ScfRestrictedDriver, molecule, basis, xcfun)
        _, rsp_results = run_tda(molecule, basis, scf_drv)

        orbrsp = TddftOrbitalResponse()
        orbrsp.update_settings({
            'tamm_dancoff': 'yes',
            'state_deriv_index': '1'
        })
        orbrsp.ostream.mute()

        # NOTE: a budget of 1 kB, which the RI set does not fit in.
        orbrsp.ri_memory_budget = 1.0e-6
        orbrsp.compute(molecule, basis, scf_drv.scf_results, rsp_results)
        orbrsp.compute_omega(molecule, basis, scf_drv.scf_results)

        assert getattr(orbrsp, fallback)
        assert not holds_ri(orbrsp)

        orbrsp.release_ri()

        assert not orbrsp._ri_j_fallback
        assert not orbrsp._ri_jk_fallback

        if MPI.COMM_WORLD.Get_size() == 1:
            # NOTE: pytest.raises is only valid in serial
            with pytest.raises(VeloxChemError, match='no RI-J'):
                orbrsp.compute_omega(molecule, basis, scf_drv.scf_results)

    @pytest.mark.parametrize('driver', [
        'hessian', 'hessian_unrestricted', 'tddft_gradient'
    ])
    def test_driver_releases_ri(self, monkeypatch, driver):

        releases = []
        release_ri = CphfSolver.release_ri

        def recorded_release_ri(solver):
            held = holds_ri(solver)
            release_ri(solver)
            releases.append((type(solver).__name__, held, holds_ri(solver)))

        monkeypatch.setattr(CphfSolver, 'release_ri', recorded_release_ri)

        if driver == 'hessian_unrestricted':
            molecule, basis = get_water(charge=1, multiplicity=2)
            scf_drv = run_scf(ScfUnrestrictedDriver, molecule, basis, 'blyp')
        else:
            molecule, basis = get_water()
            scf_drv = run_scf(ScfRestrictedDriver, molecule, basis, 'blyp')

        if driver == 'tddft_gradient':
            rsp_solver, rsp_results = run_tda(molecule, basis, scf_drv)
            grad_drv = TddftGradientDriver(scf_drv)
            grad_drv.update_settings({}, {'tamm_dancoff': 'yes'}, {}, {})
            grad_drv.ostream.mute()
            grad_drv.compute(molecule, basis, scf_drv, rsp_solver,
                             rsp_results)
        else:
            hess_drv = ScfHessianDriver(scf_drv)
            hess_drv.ostream.mute()
            hess_drv.compute(molecule, basis)

        assert len(releases) == 1
        name, held_before, held_after = releases[0]
        assert held_before
        assert not held_after
