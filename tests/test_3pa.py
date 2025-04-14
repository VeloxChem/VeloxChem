from mpi4py import MPI
import numpy as np
import pytest

from veloxchem.veloxchemlib import mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.threepatransitiondriver import ThreePATransitionDriver


class Test3PA:

    def run_scf(self, xcfun_label, basis_set_label):

        molecule_string = """
        C                 -1.96204900    0.12101600    0.00005000
        H                 -2.10925700    0.75421400    0.88802600
        H                 -2.10934300    0.75429600   -0.88785300
        H                 -2.75189700   -0.64360900    0.00005200
        C                 -0.56837700   -0.51408800   -0.00005600
        H                 -0.46300600   -1.16863200    0.88068200
        H                 -0.46308600   -1.16854400   -0.88087100
        C                  0.56837700    0.51408800   -0.00005600
        H                  0.46308600    1.16854200   -0.88087200
        H                  0.46300600    1.16863400    0.88068100
        C                  1.96204900   -0.12101600    0.00005100
        H                  2.10933600   -0.75431300   -0.88784200
        H                  2.75189700    0.64360900    0.00003100
        H                  2.10926400   -0.75419700    0.88803700"""

        molecule = Molecule.read_molecule_string(molecule_string,
                                                 units='angstrom')
        molecule.set_charge(0)
        molecule.set_multiplicity(1)

        basis = MolecularBasis.read(molecule, basis_set_label, ostream=None)

        scf_drv = ScfRestrictedDriver()
        scf_drv.ostream.mute()
        scf_drv.conv_thresh = 1.0e-8
        scf_drv.xcfun = xcfun_label

        scf_results = scf_drv.compute(molecule, basis)

        return molecule, basis, scf_results

    def run_3pa(self, xcfun_label, basis_set_label, ref_result):

        molecule, ao_basis, scf_results = self.run_scf(xcfun_label,
                                                       basis_set_label)

        rsp_settings = {'nstates': 2, 'conv_thresh': 1.0e-6}
        method_settings = {'xcfun': xcfun_label, 'grid_level': 3}

        three_pa_drv = ThreePATransitionDriver()
        three_pa_drv.ostream.mute()
        three_pa_drv.update_settings(rsp_settings, method_settings)

        results_ThreePA = three_pa_drv.compute(molecule, ao_basis, scf_results)

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            for key in ref_result:
                assert abs(results_ThreePA['transition_moments'][key] -
                           ref_result[key].real) < 1.0e-4

    @pytest.mark.timeconsuming
    def test_hf_3pa(self):

        ref_result = {
            ('xxx', -0.1467): np.complex128(0.004150342657559513 + 0j),
            ('xxy', -0.1467): np.complex128(0.005305105074934268 + 0j),
            ('xxz', -0.1467): np.complex128(-1.622616122867411e-07 + 0j),
            ('xyx', -0.1467): np.complex128(0.005305105074934268 + 0j),
            ('xyy', -0.1467): np.complex128(0.002908943381136751 + 0j),
            ('xyz', -0.1467): np.complex128(-3.1757116812818264e-09 + 0j),
            ('xzx', -0.1467): np.complex128(-1.622616122867411e-07 + 0j),
            ('xzy', -0.1467): np.complex128(-3.1757116812818264e-09 + 0j),
            ('xzz', -0.1467): np.complex128(0.004618421550471675 + 0j),
            ('yxx', -0.1467): np.complex128(0.005304964313157241 + 0j),
            ('yxy', -0.1467): np.complex128(0.0029089399271093404 + 0j),
            ('yxz', -0.1467): np.complex128(-3.1769556026000925e-09 + 0j),
            ('yyx', -0.1467): np.complex128(0.0029089399271093404 + 0j),
            ('yyy', -0.1467): np.complex128(-0.0028604367799099505 + 0j),
            ('yyz', -0.1467): np.complex128(-1.3021113623979663e-07 + 0j),
            ('yzx', -0.1467): np.complex128(-3.1769556026000925e-09 + 0j),
            ('yzy', -0.1467): np.complex128(-1.3021113623979663e-07 + 0j),
            ('yzz', -0.1467): np.complex128(0.01259096082601029 + 0j),
            ('zxx', -0.1467): np.complex128(-1.621424359970057e-07 + 0j),
            ('zxy', -0.1467): np.complex128(-3.180557771865718e-09 + 0j),
            ('zxz', -0.1467): np.complex128(0.004618395997731831 + 0j),
            ('zyx', -0.1467): np.complex128(-3.180557771865718e-09 + 0j),
            ('zyy', -0.1467): np.complex128(-1.3013450232022752e-07 + 0j),
            ('zyz', -0.1467): np.complex128(0.012591131186355792 + 0j),
            ('zzx', -0.1467): np.complex128(0.004618395997731831 + 0j),
            ('zzy', -0.1467): np.complex128(0.012591131186355792 + 0j),
            ('zzz', -0.1467): np.complex128(-9.191823779195945e-09 + 0j),
            ('xxx', -0.15): np.complex128(-6.45236935791603e-07 + 0j),
            ('xxy', -0.15): np.complex128(3.0359181089205776e-07 + 0j),
            ('xxz', -0.15): np.complex128(0.0037816305513540256 + 0j),
            ('xyx', -0.15): np.complex128(3.0359181089205776e-07 + 0j),
            ('xyy', -0.15): np.complex128(-1.5000387866468187e-08 + 0j),
            ('xyz', -0.15): np.complex128(0.004115527910078495 + 0j),
            ('xzx', -0.15): np.complex128(0.0037816305513540256 + 0j),
            ('xzy', -0.15): np.complex128(0.004115527910078495 + 0j),
            ('xzz', -0.15): np.complex128(-9.356474363678775e-08 + 0j),
            ('yxx', -0.15): np.complex128(2.986295184396112e-07 + 0j),
            ('yxy', -0.15): np.complex128(-1.5504164580806843e-08 + 0j),
            ('yxz', -0.15): np.complex128(0.004115537019178068 + 0j),
            ('yyx', -0.15): np.complex128(-1.5504164580806843e-08 + 0j),
            ('yyy', -0.15): np.complex128(-1.3254286147552143e-07 + 0j),
            ('yyz', -0.15): np.complex128(0.005578826735845969 + 0j),
            ('yzx', -0.15): np.complex128(0.004115537019178068 + 0j),
            ('yzy', -0.15): np.complex128(0.005578826735845969 + 0j),
            ('yzz', -0.15): np.complex128(-5.3222633952660965e-08 + 0j),
            ('zxx', -0.15): np.complex128(0.0037814062528413336 + 0j),
            ('zxy', -0.15): np.complex128(0.004115468138739653 + 0j),
            ('zxz', -0.15): np.complex128(-9.330536369544536e-08 + 0j),
            ('zyx', -0.15): np.complex128(0.004115468138739653 + 0j),
            ('zyy', -0.15): np.complex128(0.005578718282424726 + 0j),
            ('zyz', -0.15): np.complex128(-5.164749573297478e-08 + 0j),
            ('zzx', -0.15): np.complex128(-9.330536369544536e-08 + 0j),
            ('zzy', -0.15): np.complex128(-5.164749573297481e-08 + 0j),
            ('zzz', -0.15): np.complex128(0.0068737448714075566 + 0j)
        }

        self.run_3pa('hf', '6-31G', ref_result)

    @pytest.mark.solvers
    def test_gga_hyb_3pa(self):

        ref_result = {
            ('xxx', -0.1237): np.complex128(0.03223358436012977 + 0j),
            ('xxy', -0.1237): np.complex128(0.01834963889836454 + 0j),
            ('xxz', -0.1237): np.complex128(-6.147944414596743e-07 + 0j),
            ('xyx', -0.1237): np.complex128(0.01834963889836454 + 0j),
            ('xyy', -0.1237): np.complex128(0.004792578814914143 + 0j),
            ('xyz', -0.1237): np.complex128(-3.643217841570422e-07 + 0j),
            ('xzx', -0.1237): np.complex128(-6.147944414596743e-07 + 0j),
            ('xzy', -0.1237): np.complex128(-3.643217841570422e-07 + 0j),
            ('xzz', -0.1237): np.complex128(0.015111059474104979 + 0j),
            ('yxx', -0.1237): np.complex128(0.018349769275753696 + 0j),
            ('yxy', -0.1237): np.complex128(0.004792587368420366 + 0j),
            ('yxz', -0.1237): np.complex128(-3.6443067500127315e-07 + 0j),
            ('yyx', -0.1237): np.complex128(0.004792587368420366 + 0j),
            ('yyy', -0.1237): np.complex128(0.0019581019255517436 + 0j),
            ('yyz', -0.1237): np.complex128(3.480288865255999e-08 + 0j),
            ('yzx', -0.1237): np.complex128(-3.6443067500127315e-07 + 0j),
            ('yzy', -0.1237): np.complex128(3.480288865255999e-08 + 0j),
            ('yzz', -0.1237): np.complex128(0.002560896646880898 + 0j),
            ('zxx', -0.1237): np.complex128(-6.145646242968864e-07 + 0j),
            ('zxy', -0.1237): np.complex128(-3.6441334985787935e-07 + 0j),
            ('zxz', -0.1237): np.complex128(0.015110991083030729 + 0j),
            ('zyx', -0.1237): np.complex128(-3.6441334985787935e-07 + 0j),
            ('zyy', -0.1237): np.complex128(3.507238446686786e-08 + 0j),
            ('zyz', -0.1237): np.complex128(0.002560829732424575 + 0j),
            ('zzx', -0.1237): np.complex128(0.015110991083030729 + 0j),
            ('zzy', -0.1237): np.complex128(0.002560829732424575 + 0j),
            ('zzz', -0.1237): np.complex128(-2.0342758575499336e-07 + 0j),
            ('xxx', -0.1241): np.complex128(-2.7845207774538852e-06 + 0j),
            ('xxy', -0.1241): np.complex128(-1.8703636855634042e-06 + 0j),
            ('xxz', -0.1241): np.complex128(-0.00500069945896032 + 0j),
            ('xyx', -0.1241): np.complex128(-1.8703636855634042e-06 + 0j),
            ('xyy', -0.1241): np.complex128(-3.4826407578817637e-07 + 0j),
            ('xyz', -0.1241): np.complex128(-0.0033521651510168625 + 0j),
            ('xzx', -0.1241): np.complex128(-0.00500069945896032 + 0j),
            ('xzy', -0.1241): np.complex128(-0.0033521651510168625 + 0j),
            ('xzz', -0.1241): np.complex128(-1.1920031280782945e-06 + 0j),
            ('yxx', -0.1241): np.complex128(-1.8715369667956103e-06 + 0j),
            ('yxy', -0.1241): np.complex128(-3.485137888993821e-07 + 0j),
            ('yxz', -0.1241): np.complex128(-0.003352232257107714 + 0j),
            ('yyx', -0.1241): np.complex128(-3.48513788899382e-07 + 0j),
            ('yyy', -0.1241): np.complex128(-1.0439605813302163e-07 + 0j),
            ('yyz', -0.1241): np.complex128(0.0018287887606629495 + 0j),
            ('yzx', -0.1241): np.complex128(-0.003352232257107714 + 0j),
            ('yzy', -0.1241): np.complex128(0.0018287887606629495 + 0j),
            ('yzz', -0.1241): np.complex128(-2.2653585088933467e-07 + 0j),
            ('zxx', -0.1241): np.complex128(-0.004999754893794941 + 0j),
            ('zxy', -0.1241): np.complex128(-0.0033521450996174852 + 0j),
            ('zxz', -0.1241): np.complex128(-1.1919480592910493e-06 + 0j),
            ('zyx', -0.1241): np.complex128(-0.0033521450996174852 + 0j),
            ('zyy', -0.1241): np.complex128(0.0018292921195665953 + 0j),
            ('zyz', -0.1241): np.complex128(-2.2630304913808498e-07 + 0j),
            ('zzx', -0.1241): np.complex128(-1.1919480592910493e-06 + 0j),
            ('zzy', -0.1241): np.complex128(-2.2630304913808498e-07 + 0j),
            ('zzz', -0.1241): np.complex128(-0.003551589140212344 + 0j)
        }
        self.run_3pa('b3lyp', '6-31G', ref_result)

    @pytest.mark.timeconsuming
    def test_mgga_3pa(self):

        ref_result = {
            ('xxx', -0.1324): np.complex128(0.016143992546030717 + 0j),
            ('xxy', -0.1324): np.complex128(0.011433196209894435 + 0j),
            ('xxz', -0.1324): np.complex128(-2.909578049605096e-07 + 0j),
            ('xyx', -0.1324): np.complex128(0.011433196209894435 + 0j),
            ('xyy', -0.1324): np.complex128(0.0035939811536211835 + 0j),
            ('xyz', -0.1324): np.complex128(-9.381586200876496e-08 + 0j),
            ('xzx', -0.1324): np.complex128(-2.909578049605096e-07 + 0j),
            ('xzy', -0.1324): np.complex128(-9.381586200876496e-08 + 0j),
            ('xzz', -0.1324): np.complex128(0.010297992023668177 + 0j),
            ('yxx', -0.1324): np.complex128(0.011433194672089383 + 0j),
            ('yxy', -0.1324): np.complex128(0.0035939956683501977 + 0j),
            ('yxz', -0.1324): np.complex128(-9.380985991963196e-08 + 0j),
            ('yyx', -0.1324): np.complex128(0.0035939956683501977 + 0j),
            ('yyy', -0.1324): np.complex128(-3.772721498869502e-05 + 0j),
            ('yyz', -0.1324): np.complex128(-7.750376553842349e-08 + 0j),
            ('yzx', -0.1324): np.complex128(-9.380985991963196e-08 + 0j),
            ('yzy', -0.1324): np.complex128(-7.750376553842349e-08 + 0j),
            ('yzz', -0.1324): np.complex128(0.005960554313942167 + 0j),
            ('zxx', -0.1324): np.complex128(-2.9068801118978245e-07 + 0j),
            ('zxy', -0.1324): np.complex128(-9.38249450815472e-08 + 0j),
            ('zxz', -0.1324): np.complex128(0.010297978387162277 + 0j),
            ('zyx', -0.1324): np.complex128(-9.382494508154723e-08 + 0j),
            ('zyy', -0.1324): np.complex128(-7.722259532914235e-08 + 0j),
            ('zyz', -0.1324): np.complex128(0.005960542483420256 + 0j),
            ('zzx', -0.1324): np.complex128(0.010297978387162277 + 0j),
            ('zzy', -0.1324): np.complex128(0.005960542483420256 + 0j),
            ('zzz', -0.1324): np.complex128(-1.4274957896669311e-08 + 0j),
            ('xxx', -0.1336): np.complex128(-1.2794112770158308e-05 + 0j),
            ('xxy', -0.1336): np.complex128(-2.632560965098581e-06 + 0j),
            ('xxz', -0.1336): np.complex128(-0.004441980602565827 + 0j),
            ('xyx', -0.1336): np.complex128(-2.632560965098581e-06 + 0j),
            ('xyy', -0.1336): np.complex128(-7.512397669352436e-07 + 0j),
            ('xyz', -0.1336): np.complex128(-0.003253188067700807 + 0j),
            ('xzx', -0.1336): np.complex128(-0.004441980602565827 + 0j),
            ('xzy', -0.1336): np.complex128(-0.003253188067700807 + 0j),
            ('xzz', -0.1336): np.complex128(5.063563294068768e-08 + 0j),
            ('yxx', -0.1336): np.complex128(-2.626636762347925e-06 + 0j),
            ('yxy', -0.1336): np.complex128(-7.527295205411167e-07 + 0j),
            ('yxz', -0.1336): np.complex128(-0.003253242998265797 + 0j),
            ('yyx', -0.1336): np.complex128(-7.527295205411167e-07 + 0j),
            ('yyy', -0.1336): np.complex128(1.6417405312372955e-07 + 0j),
            ('yyz', -0.1336): np.complex128(-0.00012561528980910874 + 0j),
            ('yzx', -0.1336): np.complex128(-0.003253242998265797 + 0j),
            ('yzy', -0.1336): np.complex128(-0.00012561528980910852 + 0j),
            ('yzz', -0.1336): np.complex128(-9.287426544134523e-07 + 0j),
            ('zxx', -0.1336): np.complex128(-0.004441509229306287 + 0j),
            ('zxy', -0.1336): np.complex128(-0.003253154395228086 + 0j),
            ('zxz', -0.1336): np.complex128(5.175853821221177e-08 + 0j),
            ('zyx', -0.1336): np.complex128(-0.003253154395228086 + 0j),
            ('zyy', -0.1336): np.complex128(-0.00012537432676904464 + 0j),
            ('zyz', -0.1336): np.complex128(-9.340065723677616e-07 + 0j),
            ('zzx', -0.1336): np.complex128(5.175853821221177e-08 + 0j),
            ('zzy', -0.1336): np.complex128(-9.340065723677616e-07 + 0j),
            ('zzz', -0.1336): np.complex128(-0.004538927293386664 + 0j)
        }

        self.run_3pa('PWB6K', '6-31G', ref_result)
