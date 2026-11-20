import pytest

from veloxchem.veloxchemlib import mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.lreigensolver import LinearResponseEigenSolver
from veloxchem.tdaeigensolver import TdaEigenSolver
from veloxchem.tddftorbitalresponse import TddftOrbitalResponse


@pytest.mark.solvers
class TestTddftOrbitalResponseCpcmSettings:
    """
    A TDDFT orbital response takes the C-PCM settings of the response it
    follows from its results: cpcm_optical_epsilon and non_equilibrium_solv
    are response settings, not in the SCF results.
    """

    @pytest.mark.parametrize('solver_type', [TdaEigenSolver,
                                             LinearResponseEigenSolver])
    @pytest.mark.parametrize('non_equilibrium_solv', [True, False])
    def test_orbital_response_takes_response_cpcm_settings(
            self, solver_type, non_equilibrium_solv):

        molecule = Molecule.read_str(
            """
            O 0.0000000000    0.0000000000   -0.0254395383
            H 0.0000000000    0.7695699584    0.5948147012
            H 0.0000000000   -0.7695699584    0.5948147012
            """,
            units='angstrom')
        basis = MolecularBasis.read(molecule, 'def2-svp', ostream=None)

        scf_drv = ScfRestrictedDriver()
        scf_drv.ri = False
        scf_drv.solvation_model = 'cpcm'
        scf_drv.cpcm_epsilon = 37.5
        scf_drv.ostream.mute()
        scf_results = scf_drv.compute(molecule, basis)

        # NOTE: not the defaults of a response (1.777849 and True)
        optical_epsilon = 1.806

        rsp_solver = solver_type()
        rsp_solver.nstates = 2
        rsp_solver.cpcm_optical_epsilon = optical_epsilon
        rsp_solver.non_equilibrium_solv = non_equilibrium_solv
        rsp_solver.ostream.mute()
        rsp_results = rsp_solver.compute(molecule, basis, scf_results)

        if rsp_solver.rank == mpi_master():
            assert rsp_results['solvation_model'] == 'cpcm'
            assert rsp_results['cpcm_optical_epsilon'] == optical_epsilon
            assert rsp_results['non_equilibrium_solv'] == non_equilibrium_solv

        orbrsp = TddftOrbitalResponse()
        orbrsp.update_settings({'state_deriv_index': '1'})
        orbrsp.ostream.mute()
        orbrsp.compute(molecule, basis, scf_results, rsp_results)

        assert orbrsp._cpcm
        assert orbrsp.cpcm_optical_epsilon == optical_epsilon
        assert orbrsp.non_equilibrium_solv == non_equilibrium_solv
        assert orbrsp.cpcm_drv.optical_epsilon == optical_epsilon
