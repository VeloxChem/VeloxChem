from mpi4py import MPI
import pytest

from veloxchem.veloxchemlib import mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.quadraticresponsedriver import QuadraticResponseDriver


@pytest.mark.solvers
class TestQrfDFT:

    def run_scf(self, xcfun_label, basis_set_label):

        molecule_string = """
            O   0.0   0.0   0.0
            H   0.0   1.4   1.1
            H   0.0  -1.4   1.1
        """

        molecule = Molecule.read_molecule_string(molecule_string, units='au')
        basis = MolecularBasis.read(molecule, basis_set_label, ostream=None)

        scf_drv = ScfRestrictedDriver()
        scf_drv.ostream.mute()
        scf_drv.xcfun = xcfun_label
        scf_results = scf_drv.compute(molecule, basis)

        return molecule, basis, scf_results

    def run_qrf(self, xcfun_label, basis_set_label, ref_result):

        molecule, basis, scf_results = self.run_scf(xcfun_label,
                                                    basis_set_label)
        method_settings = {'xcfun': xcfun_label}

        qrf_drv = QuadraticResponseDriver()
        qrf_drv.ostream.mute()

        rsp_settings = {
            'conv_thresh': 1.0e-6,
            'b_frequencies': [-0.1],
            'c_frequencies': [0.3],
            'damping': 0.1,
            'a_component': 'x',
            'b_component': 'x',
            'c_component': 'z'
        }
        qrf_drv.update_settings(rsp_settings, method_settings)
        qrf_result_xxz = qrf_drv.compute(molecule, basis, scf_results)

        rsp_settings = {
            'conv_thresh': 1.0e-6,
            'b_frequencies': [-0.1],
            'c_frequencies': [0.3],
            'damping': 0.1,
            'a_component': 'y',
            'b_component': 'z',
            'c_component': 'y'
        }
        qrf_drv.update_settings(rsp_settings, method_settings)
        qrf_result_yzy = qrf_drv.compute(molecule, basis, scf_results)

        rsp_settings = {
            'conv_thresh': 1.0e-6,
            'b_frequencies': [-0.1],
            'c_frequencies': [0.3],
            'damping': 0.1,
            'a_component': 'z',
            'b_component': 'z',
            'c_component': 'z'
        }
        qrf_drv.update_settings(rsp_settings, method_settings)
        qrf_result_zzz = qrf_drv.compute(molecule, basis, scf_results)

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            thresh = 1.0e-5

            # x-component
            assert abs(qrf_result_xxz[('qrf', -0.1, 0.3)].real -
                       ref_result['xxz'].real) < thresh
            assert abs(qrf_result_xxz[('qrf', -0.1, 0.3)].imag -
                       ref_result['xxz'].imag) < thresh

            # y-component
            assert abs(qrf_result_yzy[('qrf', -0.1, 0.3)].real -
                       ref_result['yzy'].real) < thresh
            assert abs(qrf_result_yzy[('qrf', -0.1, 0.3)].imag -
                       ref_result['yzy'].imag) < thresh

            # z-component
            assert abs(qrf_result_zzz[('qrf', -0.1, 0.3)].real -
                       ref_result['zzz'].real) < thresh
            assert abs(qrf_result_zzz[('qrf', -0.1, 0.3)].imag -
                       ref_result['zzz'].imag) < thresh

    def test_lda_svp_qrf(self):

        ref_result = {
            'xxz': -0.998188471061107 - 2.5315931718063815j,
            'yzy': -27.1966181165333 - 13.152729994270668j,
            'zzz': -13.210192641695352 - 24.69010870024949j,
        }

        self.run_qrf('slda', 'def2-svp', ref_result)

    def test_gga_hyb_svp_qrf(self):

        ref_result = {
            'xxz': -1.0566216689359318 - 2.075487685485075j,
            'yzy': -26.532830780765988 - 12.366282985289274j,
            'zzz': -13.99075160995493 - 21.01055466618536j,
        }

        self.run_qrf('b3lyp', 'def2-svp', ref_result)

    def test_mgga_hyb_svp_qrf(self):

        ref_result = {
            'xxz': -1.3205225563995615 - 2.0352301782331748j,
            'yzy': -26.641495799507723 - 11.843200356445866j,
            'zzz': -15.41500249265054 - 19.8822927914182j,
        }

        self.run_qrf('tpssh', 'def2-svp', ref_result)
