import numpy as np
import pytest

from veloxchem.veloxchemlib import mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.lreigensolver import LinearResponseEigenSolver
from veloxchem.tdaeigensolver import TdaEigenSolver


@pytest.mark.solvers
class TestExcitationCpcm:

    def run_cpcm_exc(self,
                     molecule,
                     xcfun_label,
                     basis_label,
                     ref_exc,
                     ref_osc,
                     cpcm_custom_vdw_radii,
                     tol,
                     tda=False,
                     noneq_solv=False,
                     *,
                     ri):

        basis = MolecularBasis.read(molecule, basis_label)

        scf_drv = ScfRestrictedDriver()
        scf_drv.ri = ri
        scf_drv.xcfun = xcfun_label

        scf_drv.solvation_model = 'cpcm'
        scf_drv.cpcm_grid_per_sphere = (110, 110)
        scf_drv.cpcm_custom_vdw_radii = cpcm_custom_vdw_radii

        scf_drv.ostream.mute()
        scf_results = scf_drv.compute(molecule, basis)

        if scf_drv.rank == mpi_master():
            if ri:
                assert scf_results['ri_jk'] or scf_results['ri_coulomb']
                assert not scf_results['ri_fallback']

        if tda:
            rsp_drv = TdaEigenSolver()
        else:
            rsp_drv = LinearResponseEigenSolver()

        rsp_drv.ostream.mute()
        rsp_drv.nstates = 5
        # solvation model info in scf_results will be used
        # rsp_drv.solvation_model = 'cpcm'
        # rsp_drv.cpcm_grid_per_sphere = (110, 110)
        # rsp_drv.cpcm_custom_vdw_radii = cpcm_custom_vdw_radii
        rsp_drv.non_equilibrium_solv = noneq_solv

        rsp_results = rsp_drv.compute(molecule, basis, scf_results)

        if rsp_drv.rank == mpi_master():
            exc = rsp_results['eigenvalues']
            osc = rsp_results['oscillator_strengths']

            assert np.max(np.abs(exc - ref_exc)) < tol
            assert np.max(np.abs(osc - ref_osc)) < tol

    @pytest.mark.parametrize('ri', [False, True])
    def test_nh3_def2svp(self, ri):

        molstr = """
        N         -1.96309        1.59755       -0.01963
        H         -1.95876        2.61528        0.03109
        H         -2.48929        1.27814        0.79244
        H         -2.52930        1.35928       -0.83265
        """
        mol = Molecule.read_molecule_string(molstr, units='bohr')

        # LR-HF

        ref_exc = {
            False: np.array(
                [0.40521822, 0.4860046, 0.48603203, 0.77782504, 0.77783245]),
            True: np.array(
                [0.40523372, 0.48604612, 0.48607355, 0.77785390, 0.77786131])
        }[ri]
        ref_osc = {
            False: np.array(
                [0.06158212, 0.01578646, 0.01578704, 0.12326755, 0.12328609]),
            True: np.array(
                [0.06159313, 0.01580349, 0.01580410, 0.12328536, 0.12330146])
        }[ri]
        self.run_cpcm_exc(mol,
                          'hf',
                          'def2-svp',
                          ref_exc,
                          ref_osc,
                          None,
                          1.0e-5,
                          ri=ri)

        # TDA-HF

        ref_exc = {
            False: np.array(
                [0.40578896, 0.48668494, 0.48671237, 0.77825217, 0.7782596]),
            True: np.array(
                [0.40580432, 0.48672642, 0.48675384, 0.77828084, 0.77828826])
        }[ri]
        ref_osc = {
            False: np.array(
                [0.06022241, 0.01533443, 0.01533621, 0.11997976, 0.11999964]),
            True: np.array(
                [0.06023306, 0.01535261, 0.01535439, 0.11999324, 0.12001255])
        }[ri]
        self.run_cpcm_exc(mol,
                          'hf',
                          'def2-svp',
                          ref_exc,
                          ref_osc,
                          None,
                          1.0e-5,
                          tda=True,
                          ri=ri)

        # LR-HF + Non-eq. solvation

        ref_exc = {
            False: np.array(
                [0.40759806, 0.48702943, 0.48705686, 0.78026947, 0.78027528]),
            True: np.array(
                [0.40761395, 0.48707187, 0.48709929, 0.78029860, 0.78030440])
        }[ri]
        ref_osc = {
            False: np.array(
                [0.05766621, 0.01458069, 0.0145834, 0.11290953, 0.11292115]),
            True: np.array(
                [0.05767554, 0.01459643, 0.01459916, 0.11292518, 0.11293679])
        }[ri]
        self.run_cpcm_exc(mol,
                          'hf',
                          'def2-svp',
                          ref_exc,
                          ref_osc,
                          None,
                          1.0e-5,
                          noneq_solv=True,
                          ri=ri)

        # LR-DFT

        ref_exc = {
            False: np.array(
                [0.34592948, 0.42604827, 0.42607682, 0.68311411, 0.68312656]),
            True: np.array(
                [0.34594433, 0.42608560, 0.42611415, 0.68315692, 0.68316938])
        }[ri]
        ref_osc = {
            False: np.array(
                [0.04423422, 0.01472548, 0.01472579, 0.09518219, 0.095212]),
            True: np.array(
                [0.04423968, 0.01474088, 0.01474140, 0.09517495, 0.09520548])
        }[ri]
        self.run_cpcm_exc(mol,
                          'b3lyp',
                          'def2-svp',
                          ref_exc,
                          ref_osc,
                          None,
                          1.0e-5,
                          ri=ri)

        # LR-HF + custom vdW radii

        ref_exc = {
            False: np.array(
                [0.40041427, 0.4813658, 0.48136975, 0.7765264, 0.7765328]),
            True: np.array(
                [0.40042899, 0.48140782, 0.48141177, 0.77655462, 0.77656102])
        }[ri]
        ref_osc = {
            False: np.array(
                [0.05829868, 0.0139734, 0.01397378, 0.11872162, 0.11872378]),
            True: np.array(
                [0.05830918, 0.01398858, 0.01398898, 0.11873810, 0.11873836])
        }[ri]
        self.run_cpcm_exc(mol,
                          'hf',
                          'def2-svp',
                          ref_exc,
                          ref_osc, ['N', 1.8],
                          1.0e-5,
                          ri=ri)
