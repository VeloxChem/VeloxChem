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
                     noneq_solv=False):

        basis = MolecularBasis.read(molecule, basis_label)

        scf_drv = ScfRestrictedDriver()
        scf_drv.xcfun = xcfun_label

        scf_drv.solvation_model = 'cpcm'
        scf_drv.cpcm_grid_per_sphere = (110, 110)
        scf_drv.cpcm_custom_vdw_radii = cpcm_custom_vdw_radii

        scf_drv.ostream.mute()
        scf_results = scf_drv.compute(molecule, basis)

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

    def test_nh3_def2svp(self):

        molstr = """
        N         -1.96309        1.59755       -0.01963
        H         -1.95876        2.61528        0.03109
        H         -2.48929        1.27814        0.79244
        H         -2.52930        1.35928       -0.83265
        """
        mol = Molecule.read_molecule_string(molstr, units='bohr')

        # LR-HF

        ref_exc = np.array(
            [0.40521822, 0.4860046, 0.48603203, 0.77782504, 0.77783245])
        ref_osc = np.array(
            [0.06158212, 0.01578646, 0.01578704, 0.12326755, 0.12328609])
        self.run_cpcm_exc(mol, 'hf', 'def2-svp', ref_exc, ref_osc, None, 1.0e-5)

        # TDA-HF

        ref_exc = np.array(
            [0.40578896, 0.48668494, 0.48671237, 0.77825217, 0.7782596])
        ref_osc = np.array(
            [0.06022241, 0.01533443, 0.01533621, 0.11997976, 0.11999964])
        self.run_cpcm_exc(mol,
                          'hf',
                          'def2-svp',
                          ref_exc,
                          ref_osc,
                          None,
                          1.0e-5,
                          tda=True)

        # LR-HF + Non-eq. solvation

        ref_exc = np.array(
            [0.40759806, 0.48702943, 0.48705686, 0.78026947, 0.78027528])
        ref_osc = np.array(
            [0.05766621, 0.01458069, 0.0145834, 0.11290953, 0.11292115])
        self.run_cpcm_exc(mol,
                          'hf',
                          'def2-svp',
                          ref_exc,
                          ref_osc,
                          None,
                          1.0e-5,
                          noneq_solv=True)

        # LR-DFT

        ref_exc = np.array(
            [0.34592948, 0.42604827, 0.42607682, 0.68311411, 0.68312656])
        ref_osc = np.array(
            [0.04423422, 0.01472548, 0.01472579, 0.09518219, 0.095212])
        self.run_cpcm_exc(mol, 'b3lyp', 'def2-svp', ref_exc, ref_osc, None,
                          1.0e-5)

        # LR-HF + custom vdW radii

        ref_exc = np.array(
            [0.40041427, 0.4813658, 0.48136975, 0.7765264, 0.7765328])
        ref_osc = np.array(
            [0.05829868, 0.0139734, 0.01397378, 0.11872162, 0.11872378])
        self.run_cpcm_exc(mol, 'hf', 'def2-svp', ref_exc, ref_osc, ['N', 1.8],
                          1.0e-5)
