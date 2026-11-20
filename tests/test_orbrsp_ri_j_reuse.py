from mpi4py import MPI
import pytest

from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.tdaeigensolver import TdaEigenSolver
from veloxchem.lrsolver import LinearResponseSolver
from veloxchem.tddftorbitalresponse import TddftOrbitalResponse
from veloxchem.polorbitalresponse import PolOrbitalResponse
from veloxchem.errorhandler import VeloxChemError


@pytest.fixture(scope='module')
def ground_state():

    molecule = Molecule.read_str(
        """
        O 0.0000000000    0.0000000000   -0.0254395383
        H 0.0000000000    0.7695699584    0.5948147012
        H 0.0000000000   -0.7695699584    0.5948147012
        """,
        units='angstrom')
    basis = MolecularBasis.read(molecule, 'def2-svp', ostream=None)

    scf_drv = ScfRestrictedDriver()
    scf_drv.ri = True
    scf_drv.xcfun = 'blyp'
    scf_drv.ostream.mute()
    scf_results = scf_drv.compute(molecule, basis)

    assert scf_drv._uses_ri_coulomb()

    return molecule, basis, scf_results


class TestOrbitalResponseRiJReuse:
    """
    The omega of an orbital response uses the RI-J driver its compute set up,
    or the fallback its compute took, and does not decide again.
    """

    def make_solver(self, kind, molecule, basis, scf_results):

        if kind == 'tddft':
            rsp_solver = TdaEigenSolver()
            rsp_solver.nstates = 2
            rsp_solver.ostream.mute()
            rsp_results = rsp_solver.compute(molecule, basis, scf_results)

            solver = TddftOrbitalResponse()
            solver.update_settings({
                'tamm_dancoff': 'yes',
                'state_deriv_index': '1'
            })

            def compute():
                solver.compute(molecule, basis, scf_results, rsp_results)

            def compute_omega():
                solver.compute_omega(molecule, basis, scf_results)

        else:
            lr_drv = LinearResponseSolver()
            lr_drv.a_operator = 'electric dipole'
            lr_drv.b_operator = 'electric dipole'
            lr_drv.update_settings({'frequencies': (0.0, 0.4)})
            lr_drv.ostream.mute()
            lr_results = lr_drv.compute(molecule, basis, scf_results)

            solver = PolOrbitalResponse()
            solver.update_settings({'frequencies': (0.0, 0.4)})

            def compute():
                solver.compute(molecule, basis, scf_results, lr_results)

            def compute_omega():
                solver.compute_omega(molecule, basis, scf_results, lr_results)

        solver.ostream.mute()

        return solver, compute, compute_omega

    @pytest.mark.parametrize('kind', ['tddft', 'polgrad'])
    def test_omega_uses_the_driver_of_compute(self, ground_state, kind):

        solver, compute, compute_omega = self.make_solver(kind, *ground_state)

        compute()

        assert solver._uses_ri_coulomb()
        ri_drv = solver._ri_drv
        assert ri_drv is not None

        compute_omega()

        assert solver._ri_drv is ri_drv

    @pytest.mark.parametrize('kind', ['tddft', 'polgrad'])
    def test_omega_keeps_the_fallback_of_compute(self, ground_state, kind):

        solver, compute, compute_omega = self.make_solver(kind, *ground_state)

        # NOTE: a budget of 1 kB, which the RI-J integrals do not fit in.
        solver.ri_memory_budget = 1.0e-6
        compute()

        assert solver._ri_j_fallback

        # NOTE: with no budget set, a new decision would take RI-J.
        solver.ri_memory_budget = None
        compute_omega()

        assert solver._ri_j_fallback
        assert solver._ri_drv is None

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason='pytest.raises only valid in serial')
    def test_omega_without_compute_is_refused(self, ground_state):

        molecule, basis, scf_results = ground_state

        solver = TddftOrbitalResponse()
        solver.ostream.mute()

        with pytest.raises(VeloxChemError, match='no RI-J driver'):
            solver.compute_omega(molecule, basis, scf_results)
