import numpy as np
import pytest

from veloxchem.veloxchemlib import mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.cppsolver import ComplexResponseSolver


@pytest.mark.solvers
class TestCPPCpcm:

    def run_cpp_cpcm(self,
                     xcfun_label,
                     cpp_property,
                     ref_x_data,
                     ref_y_data,
                     tol,
                     noneq_solv=False,
                     *,
                     ri):

        xyz_string = """6
        xyz
        O -3.42904  1.55532  0.01546
        C -1.99249  1.74379  0.02665
        H -1.74709  2.74160  0.44749
        H -1.59636  1.67836 -1.00861
        H -1.51398  0.95881  0.64937
        H -3.84726  2.33620 -0.34927
        """
        mol = Molecule.read_xyz_string(xyz_string)

        basis_label = 'def2-svp'
        bas = MolecularBasis.read(mol, basis_label, ostream=None)

        # Note: vacuum-SCF + CPCM-response for testing purposes only

        scf_drv = ScfRestrictedDriver()
        scf_drv.ri = ri
        scf_drv.ostream.mute()
        scf_drv.xcfun = xcfun_label
        scf_results = scf_drv.compute(mol, bas)

        if scf_drv.rank == mpi_master():
            if ri:
                assert scf_results['ri_jk'] or scf_results['ri_coulomb']
                assert not scf_results['ri_fallback']

        lr_drv = ComplexResponseSolver()
        lr_drv.ostream.mute()
        lr_drv.property = cpp_property
        lr_drv.frequencies = list(ref_x_data)
        lr_drv.solvation_model = 'cpcm'
        lr_drv.cpcm_grid_per_sphere = (110, 110)
        lr_drv.non_equilibrium_solv = noneq_solv
        lr_results = lr_drv.compute(mol, bas, scf_results)

        if lr_drv.rank == mpi_master():
            lr_spec = lr_drv.get_spectrum(lr_results, 'au')
            assert np.max(np.abs(np.array(lr_spec['y_data']) -
                                 ref_y_data)) < tol

    @pytest.mark.parametrize('ri', [False, True])
    def test_b3lyp_absorption_cpcm(self, ri):

        # vlxtag: RKS, Absorption, CPP, CPCM

        xcfun_label = 'b3lyp'
        cpp_property = 'absorption'
        ref_x_data = [0.39, 0.40, 0.41]
        ref_y_data = {
            False: [0.32158440, 0.66111297, 2.15027316],
            True: [0.32163522, 0.66152443, 2.14998963]
        }[ri]

        self.run_cpp_cpcm(xcfun_label,
                          cpp_property,
                          ref_x_data,
                          ref_y_data,
                          1.0e-4,
                          ri=ri)

    @pytest.mark.parametrize('ri', [False, True])
    def test_b3lyp_absorption_cpcm_noneq(self, ri):

        # vlxtag: RKS, Absorption, CPP, CPCM

        xcfun_label = 'b3lyp'
        cpp_property = 'absorption'
        ref_x_data = [0.39, 0.40, 0.41]
        ref_y_data = {
            False: [0.23785824, 0.32046951, 1.27698344],
            True: [0.23789533, 0.32059839, 1.27758983]
        }[ri]

        self.run_cpp_cpcm(xcfun_label,
                          cpp_property,
                          ref_x_data,
                          ref_y_data,
                          1.0e-4,
                          noneq_solv=True,
                          ri=ri)
