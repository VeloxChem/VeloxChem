import pytest

from veloxchem.veloxchemlib import is_mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.quadraticresponsedriver import QuadraticResponseDriver


@pytest.mark.solvers
class TestQrf:

    def run_scf(self):

        molecule_string = """
            O   0.0   0.0   0.0
            H   0.0   1.4   1.1
            H   0.0  -1.4   1.1
        """
        basis_set_label = 'cc-pVDZ'
        xcfun_label = 'B3LYP'

        molecule = Molecule.read_str(molecule_string, units='au')
        basis = MolecularBasis.read(molecule, basis_set_label, ostream=None)

        scf_drv = ScfRestrictedDriver()
        scf_drv.ostream.state = False
        scf_drv.update_settings({}, {'xcfun': xcfun_label})
        scf_drv.compute(molecule, basis)

        return molecule, basis, xcfun_label, scf_drv.scf_tensors

    def run_qrf(self, ref_result):

        molecule, basis, xcfun_label, scf_tensors = self.run_scf()
        method_settings = {'xcfun': xcfun_label}

        qrf_drv = QuadraticResponseDriver()
        qrf_drv.ostream.state = False

        rsp_settings = {
            'conv_thresh': 1.0e-6,
            'b_frequencies': [-0.1],
            'c_frequencies': [0.3],
            'damping': 0.1,
            'a_components': 'x',
            'b_components': 'x',
            'c_components': 'z'
        }
        qrf_drv.update_settings(rsp_settings, method_settings)
        qrf_result_xxz = qrf_drv.compute(molecule, basis, scf_tensors)

        rsp_settings = {
            'conv_thresh': 1.0e-6,
            'b_frequencies': [-0.1],
            'c_frequencies': [0.3],
            'damping': 0.1,
            'a_components': 'y',
            'b_components': 'z',
            'c_components': 'y'
        }
        qrf_drv.update_settings(rsp_settings, method_settings)
        qrf_result_yzy = qrf_drv.compute(molecule, basis, scf_tensors)

        rsp_settings = {
            'conv_thresh': 1.0e-6,
            'b_frequencies': [-0.1],
            'c_frequencies': [0.3],
            'damping': 0.1,
            'a_components': 'z',
            'b_components': 'z',
            'c_components': 'z'
        }
        qrf_drv.update_settings(rsp_settings, method_settings)
        qrf_result_zzz = qrf_drv.compute(molecule, basis, scf_tensors)

        if is_mpi_master():
            thresh = 1.0e-4

            # x-component
            assert abs(qrf_result_xxz[(-0.1, 0.3)].real -
                       ref_result['xxz'].real) < thresh
            assert abs(qrf_result_xxz[(-0.1, 0.3)].imag -
                       ref_result['xxz'].imag) < thresh

            # y-component
            assert abs(qrf_result_yzy[(-0.1, 0.3)].real -
                       ref_result['yzy'].real) < thresh
            assert abs(qrf_result_yzy[(-0.1, 0.3)].imag -
                       ref_result['yzy'].imag) < thresh

            # z-component
            assert abs(qrf_result_zzz[(-0.1, 0.3)].real -
                       ref_result['zzz'].real) < thresh
            assert abs(qrf_result_zzz[(-0.1, 0.3)].imag -
                       ref_result['zzz'].imag) < thresh

    def test_qrf(self):

        ref_result = {
            'xxz': 1.89373381 + 2.46156092j,
            'yzy': 24.52654293 + 11.22800902j,
            'zzz': 13.51752880 + 18.68141118j,
        }

        self.run_qrf(ref_result)
