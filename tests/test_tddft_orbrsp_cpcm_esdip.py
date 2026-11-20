import numpy as np
import pytest

from veloxchem.veloxchemlib import mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.lreigensolver import LinearResponseEigenSolver
from veloxchem.tdaeigensolver import TdaEigenSolver
from veloxchem.tddftgradientdriver import TddftGradientDriver
from veloxchem.firstorderprop import FirstOrderProperties


@pytest.mark.solvers
class TestTddftOrbitalResponseCpcmDipole:
    """
    The relaxed dipole of the first excited state with C-PCM, equilibrium and
    non-equilibrium, against finite differences of its energy in a field.
    """

    # NOTE: references: -dE/dF_z of E_SCF + omega_1 by central differences,
    # F_z = +-1e-3 a.u., SCF conv_thresh 1e-10, response conv_thresh 1e-8, the
    # integrals of the row (four-center or RI). An electric field together
    # with C-PCM is refused by environment_compatibility_sanity_check; it was
    # replaced by a function that returns for those runs. The analytic dipoles
    # of this test were within 2e-6 of them.
    @pytest.mark.parametrize(
        'xcfun, ri, solver_type, non_equilibrium_solv, ref_dipole_z', [
            ('hf', False, TdaEigenSolver, False, -0.260472),
            ('hf', False, TdaEigenSolver, True, -0.261823),
            ('hf', False, LinearResponseEigenSolver, True, -0.268567),
            ('b3lyp', False, TdaEigenSolver, False, -0.307216),
            ('b3lyp', False, TdaEigenSolver, True, -0.307370),
            ('b3lyp', True, TdaEigenSolver, True, -0.307308),
            ('b3lyp', True, LinearResponseEigenSolver, True, -0.304020),
        ])
    def test_relaxed_dipole(self, xcfun, ri, solver_type,
                            non_equilibrium_solv, ref_dipole_z):

        molecule = Molecule.read_str(
            """
            O 0.0000000000    0.0000000000   -0.0254395383
            H 0.0000000000    0.7695699584    0.5948147012
            H 0.0000000000   -0.7695699584    0.5948147012
            """,
            units='angstrom')
        basis = MolecularBasis.read(molecule, 'def2-svp', ostream=None)

        scf_drv = ScfRestrictedDriver()
        scf_drv.ri = ri
        scf_drv.xcfun = xcfun
        scf_drv.conv_thresh = 1.0e-10
        scf_drv.solvation_model = 'cpcm'
        scf_drv.cpcm_epsilon = 37.5
        scf_drv.ostream.mute()
        scf_results = scf_drv.compute(molecule, basis)

        rsp_solver = solver_type()
        rsp_solver.nstates = 1
        rsp_solver.conv_thresh = 1.0e-8
        rsp_solver.cpcm_optical_epsilon = 1.806
        rsp_solver.non_equilibrium_solv = non_equilibrium_solv
        rsp_solver.ostream.mute()
        rsp_results = rsp_solver.compute(molecule, basis, scf_results)

        grad_drv = TddftGradientDriver(scf_drv)
        grad_drv.state_deriv_index = [1]
        grad_drv.ostream.mute()
        densities = grad_drv.compute_excited_state_densities(
            molecule, basis, rsp_results)

        if scf_drv.rank == mpi_master():
            relaxed_density = (densities['relaxed_density'][0] +
                               scf_results['D_alpha'] + scf_results['D_beta'])
        else:
            relaxed_density = None

        prop = FirstOrderProperties()
        prop.compute(molecule, basis, relaxed_density)

        if scf_drv.rank == mpi_master():
            dipole = prop.get_property('dipole_moment')
            assert np.max(np.abs(dipole[:2])) < 1.0e-6
            assert abs(dipole[2] - ref_dipole_z) < 1.0e-5
