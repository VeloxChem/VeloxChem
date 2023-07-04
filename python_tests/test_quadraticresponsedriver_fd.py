from mpi4py import MPI
import pytest

from veloxchem.veloxchemlib import is_single_node, is_mpi_master
from veloxchem.outputstream import OutputStream
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.cppsolver import ComplexResponse
from veloxchem.quadraticresponsedriver import QuadraticResponseDriver


class TestQrfFD:

    def run_qrf_fd(self, xcfun_label, basis_set_label, components, freqs):

        comm = MPI.COMM_WORLD
        ostream = OutputStream(None)

        molecule_string = """
        O   0.0   0.0   0.0
        H   0.0   1.4   1.1
        H   0.0  -1.4   1.1
        """
        molecule = Molecule.read_molecule_string(molecule_string, units='au')
        basis = MolecularBasis.read(molecule, basis_set_label, ostream=None)

        a, b, c = components
        wb, wc = freqs

        scf_conv_thresh = 1.0e-8
        rsp_conv_thresh = 1.0e-5

        # SCF

        scf_settings = {'conv_thresh': scf_conv_thresh}
        method_settings = {'xcfun': xcfun_label}

        scfdrv = ScfRestrictedDriver(comm, ostream)
        scfdrv.update_settings(scf_settings, method_settings)
        scf_result = scfdrv.compute(molecule, basis)

        # QRF

        qrf = QuadraticResponseDriver(comm, ostream)

        rsp_settings = {
            'conv_thresh': rsp_conv_thresh,
            'a_components': a,
            'b_components': b,
            'c_components': c,
            'b_frequencies': [wb],
            'c_frequencies': [wc],
        }
        qrf.update_settings(rsp_settings, method_settings)
        quad_result = qrf.compute(molecule, basis, scf_result)

        if is_mpi_master():
            beta = -quad_result[(wb, wc)]

        # permutation test

        rsp_settings = {
            'conv_thresh': rsp_conv_thresh,
            'a_components': a,
            'b_components': c,
            'c_components': b,
            'b_frequencies': [wc],
            'c_frequencies': [wb],
        }
        qrf.update_settings(rsp_settings, method_settings)
        quad_result = qrf.compute(molecule, basis, scf_result)

        if is_mpi_master():
            beta_perm = -quad_result[(wc, wb)]
            assert abs(beta.real - beta_perm.real) < 1.0e-6
            assert abs(beta.imag - beta_perm.imag) < 1.0e-6

        # conjugate test

        rsp_settings = {
            'conv_thresh': rsp_conv_thresh,
            'a_components': a,
            'b_components': b,
            'c_components': c,
            'b_frequencies': [-wb],
            'c_frequencies': [-wc],
        }
        qrf.update_settings(rsp_settings, method_settings)
        quad_result = qrf.compute(molecule, basis, scf_result)

        if is_mpi_master():
            beta_conj = -quad_result[(-wb, -wc)]
            assert abs(beta.real - beta_conj.real) < 1.0e-6
            assert abs(beta.imag + beta_conj.imag) < 1.0e-6

        # finite difference test

        rsp_settings = {
            'conv_thresh': rsp_conv_thresh,
            'a_components': a,
            'b_components': b,
            'c_components': c,
            'b_frequencies': [wb],
            'c_frequencies': [0],
            'damping': 0,
        }
        qrf.update_settings(rsp_settings, method_settings)
        qrf_result = qrf.compute(molecule, basis, scf_result)

        if is_mpi_master():
            beta_0 = -qrf_result[(wb, 0)].real
            assert abs(-qrf_result[(wb, 0)].imag) < 1.0e-6

        cpp_settings = {
            'conv_thresh': rsp_conv_thresh,
            'frequencies': [wb],
            'damping': 0,
        }

        delta_ef = 1.0e-4
        fd_index = {'x': 0, 'y': 1, 'z': 2}[c]

        efield_plus = [0.0, 0.0, 0.0]
        efield_minus = [0.0, 0.0, 0.0]
        efield_plus[fd_index] = delta_ef
        efield_minus[fd_index] = -delta_ef

        method_dict_plus = dict(method_settings)
        method_dict_minus = dict(method_settings)
        method_dict_plus['electric_field'] = efield_plus
        method_dict_minus['electric_field'] = efield_minus

        scf_drv_plus = ScfRestrictedDriver(comm, ostream)
        scf_drv_plus.update_settings(scf_settings, method_dict_plus)
        scf_result_plus = scf_drv_plus.compute(molecule, basis)

        cpp_plus = ComplexResponse(comm, ostream)
        cpp_plus.update_settings(cpp_settings, method_dict_plus)
        cpp_result_plus = cpp_plus.compute(molecule, basis, scf_result_plus)

        scf_drv_minus = ScfRestrictedDriver(comm, ostream)
        scf_drv_minus.update_settings(scf_settings, method_dict_minus)
        scf_result_minus = scf_drv_minus.compute(molecule, basis)

        cpp_minus = ComplexResponse(comm, ostream)
        cpp_minus.update_settings(cpp_settings, method_dict_minus)
        cpp_result_minus = cpp_minus.compute(molecule, basis, scf_result_minus)

        if is_mpi_master():
            alpha_plus = -cpp_result_plus['response_functions'][(a, b, wb)]
            alpha_minus = -cpp_result_minus['response_functions'][(a, b, wb)]
            assert abs(alpha_plus.imag) < 1.0e-6
            assert abs(alpha_minus.imag) < 1.0e-6

            beta_0_fd = (alpha_plus.real - alpha_minus.real) / (2.0 * delta_ef)
            assert abs(beta_0 - beta_0_fd) / abs(beta_0_fd) < 1.0e-5

    @pytest.mark.skipif(is_single_node(), reason="multi-node only")
    def test_lda_qrf_fd(self):

        self.run_qrf_fd('slda', 'def2-svp', 'yyz', [0.3, -0.11])

    @pytest.mark.skipif(is_single_node(), reason="multi-node only")
    def test_gga_qrf_fd(self):

        self.run_qrf_fd('bp86', 'def2-svp', 'yyz', [0.3, -0.11])

    @pytest.mark.skipif(is_single_node(), reason="multi-node only")
    def test_mgga_qrf_fd(self):

        self.run_qrf_fd('tpssh', 'def2-svp', 'yyz', [0.3, -0.11])
