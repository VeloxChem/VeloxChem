import pytest

import veloxchem.linearsolver as linearsolver
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.tdaeigensolver import TdaEigenSolver
from veloxchem.tddftorbitalresponse import TddftOrbitalResponse


class TestTddftOrbitalResponseRiSetup:
    """
    The compute of a TDDFT orbital response sets up its RI once: the right-hand
    side uses what CphfSolver.compute_subspace_solver set up for the solve.
    """

    @pytest.mark.parametrize('xcfun, setup_name, uses_ri', [
        ('blyp', 'initialize_ri_j', '_uses_ri_coulomb'),
        ('b3lyp', 'initialize_ri_jk', '_uses_ri_jk'),
    ])
    def test_compute_sets_up_ri_once(self, monkeypatch, xcfun, setup_name,
                                     uses_ri):

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
        scf_drv.xcfun = xcfun
        scf_drv.ostream.mute()
        scf_results = scf_drv.compute(molecule, basis)

        rsp_solver = TdaEigenSolver()
        rsp_solver.nstates = 2
        rsp_solver.ostream.mute()
        rsp_results = rsp_solver.compute(molecule, basis, scf_results)

        orbrsp = TddftOrbitalResponse()
        orbrsp.update_settings({
            'tamm_dancoff': 'yes',
            'state_deriv_index': '1'
        })
        orbrsp.ostream.mute()

        setup = getattr(linearsolver, setup_name)
        setups = []

        def counted_setup(solver, molecule, basis, reuse=False):
            if solver is orbrsp:
                setups.append(reuse)
            return setup(solver, molecule, basis, reuse)

        monkeypatch.setattr(linearsolver, setup_name, counted_setup)

        orbrsp.compute(molecule, basis, scf_results, rsp_results)

        assert getattr(orbrsp, uses_ri)()
        assert setups == [False]
