import numpy as np
import pytest

from veloxchem.lreigensolver import LinearResponseEigenSolver
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.molecule import Molecule
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.tdaeigensolver import TdaEigenSolver
from veloxchem.tddftgradientdriver import TddftGradientDriver
from veloxchem.veloxchemlib import mpi_master


@pytest.mark.solvers
class TestTddftGradientRiJk:

    def run_gradient(self, molecule, basis, xcfun, tamm_dancoff, use_rijk):

        method_dict = {'xcfun': xcfun, 'grid_level': 5}

        if use_rijk:
            method_dict.update({
                'ri_jk': True,
                'ri_jk_simd': True,
                'ri_auxiliary_basis': 'def2-universal-jkfit',
            })

        scf_drv = ScfRestrictedDriver()
        scf_drv.update_settings({'conv_thresh': 1e-8}, method_dict)
        scf_drv.ostream.mute()
        scf_drv.compute(molecule, basis)

        rsp_dict = {'tamm_dancoff': 'yes' if tamm_dancoff else 'no'}

        if tamm_dancoff:
            rsp_drv = TdaEigenSolver()
        else:
            rsp_drv = LinearResponseEigenSolver()

        rsp_drv.update_settings(rsp_dict, method_dict)
        rsp_drv.conv_thresh = 1e-5
        rsp_drv.nstates = 2
        rsp_drv.ostream.mute()
        rsp_results = rsp_drv.compute(molecule, basis, scf_drv.scf_results)

        grad_drv = TddftGradientDriver(scf_drv)
        grad_drv.update_settings({}, rsp_dict, {'conv_thresh': 1e-7},
                                 method_dict)
        grad_drv.ostream.mute()
        grad_drv.compute(molecule, basis, scf_drv, rsp_drv, rsp_results)

        if grad_drv.rank == mpi_master():
            return grad_drv.get_gradient()

        return None

    @pytest.mark.parametrize('xcfun, tamm_dancoff', [
        ('hf', True),
        ('b3lyp', False),
        ('cam-b3lyp', True),
    ])
    def test_rijk_matches_four_center(self, xcfun, tamm_dancoff):

        molecule = Molecule.read_str(
            """
            O 0.0000000000 0.0000000000 -0.0254395383
            H 0.0000000000 0.7695699584 0.5948147012
            H 0.0000000000 -0.7695699584 0.5948147012
            """, units='angstrom')

        basis = MolecularBasis.read(molecule, 'def2-svp', ostream=None)

        ref_grad = self.run_gradient(molecule, basis, xcfun, tamm_dancoff,
                                     False)
        rijk_grad = self.run_gradient(molecule, basis, xcfun, tamm_dancoff,
                                      True)

        if ref_grad is not None:
            max_diff = np.max(np.abs(ref_grad - rijk_grad))

            # The RI-JK approximation must differ from the four-centre result,
            # otherwise the branch could silently have fallen back.
            assert max_diff > 1.0e-8
            assert max_diff < 1.0e-4
