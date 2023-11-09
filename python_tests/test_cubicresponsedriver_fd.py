from mpi4py import MPI
import pytest

from veloxchem.veloxchemlib import is_single_node, is_mpi_master
from veloxchem.outputstream import OutputStream
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.quadraticresponsedriver import QuadraticResponseDriver
from veloxchem.cubicresponsedriver import CubicResponseDriver


class TestCrfFD:

    def run_crf_fd(self, xcfun_label, basis_set_label, components, freqs):

        comm = MPI.COMM_WORLD
        ostream = OutputStream(None)

        molecule_string = """
        O   0.0   0.0   0.0
        H   0.0   1.4   1.1
        H   0.0  -1.4   1.1
        """
        molecule = Molecule.read_molecule_string(molecule_string, units='au')
        basis = MolecularBasis.read(molecule, basis_set_label, ostream=None)

        a, b, c, d = components
        wb, wc, wd = freqs

        scf_conv_thresh = 1.0e-8
        rsp_conv_thresh = 1.0e-5

        # SCF

        scf_settings = {'conv_thresh': scf_conv_thresh}
        method_settings = {'xcfun': xcfun_label}

        scfdrv = ScfRestrictedDriver(comm, ostream)
        scfdrv.update_settings(scf_settings, method_settings)
        scf_result = scfdrv.compute(molecule, basis)

        # CRF

        crf = CubicResponseDriver(comm, ostream)

        rsp_settings = {
            'conv_thresh': rsp_conv_thresh,
            'a_component': a,
            'b_component': b,
            'c_component': c,
            'd_component': d,
            'b_frequencies': [wb],
            'c_frequencies': [wc],
            'd_frequencies': [wd],
        }
        crf.update_settings(rsp_settings, method_settings)
        crf_result = crf.compute(molecule, basis, scf_result)

        if is_mpi_master():
            gamma = -crf_result[('gamma', wb, wc, wd)]

        # permutation test

        rsp_settings = {
            'conv_thresh': rsp_conv_thresh,
            'a_component': a,
            'b_component': b,
            'c_component': d,
            'd_component': c,
            'b_frequencies': [wb],
            'c_frequencies': [wd],
            'd_frequencies': [wc],
        }
        crf.update_settings(rsp_settings, method_settings)
        crf_result = crf.compute(molecule, basis, scf_result)

        if is_mpi_master():
            gamma_perm = -crf_result[('gamma', wb, wd, wc)]
            assert abs(gamma.real - gamma_perm.real) < 1.0e-6
            assert abs(gamma.imag - gamma_perm.imag) < 1.0e-6

        # conjugate test

        rsp_settings = {
            'conv_thresh': rsp_conv_thresh,
            'a_component': a,
            'b_component': b,
            'c_component': c,
            'd_component': d,
            'b_frequencies': [-wb],
            'c_frequencies': [-wc],
            'd_frequencies': [-wd],
        }
        crf.update_settings(rsp_settings, method_settings)
        crf_result = crf.compute(molecule, basis, scf_result)

        if is_mpi_master():
            gamma_conj = -crf_result[('gamma', -wb, -wc, -wd)]
            assert abs(gamma.real - gamma_conj.real) < 1.0e-6
            assert abs(gamma.imag + gamma_conj.imag) < 1.0e-6

        # finite difference test

        rsp_settings = {
            'conv_thresh': rsp_conv_thresh,
            'a_component': a,
            'b_component': b,
            'c_component': c,
            'd_component': d,
            'b_frequencies': [wb],
            'c_frequencies': [wc],
            'd_frequencies': [0],
            'damping': 0,
        }
        crf.update_settings(rsp_settings, method_settings)
        crf_result = crf.compute(molecule, basis, scf_result)

        if is_mpi_master():
            gamma_0 = -crf_result[('gamma', wb, wc, 0)].real
            assert abs(-crf_result[('gamma', wb, wc, 0)].imag) < 1.0e-6

        qrf_settings = {
            'conv_thresh': rsp_conv_thresh,
            'a_component': a,
            'b_component': b,
            'c_component': c,
            'b_frequencies': [wb],
            'c_frequencies': [wc],
            'damping': 0,
        }

        delta_ef = 1.0e-4
        fd_index = {'x': 0, 'y': 1, 'z': 2}[d]

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

        qrf_plus = QuadraticResponseDriver(comm, ostream)
        qrf_plus.update_settings(qrf_settings, method_settings)
        quad_result_plus = qrf_plus.compute(molecule, basis, scf_result_plus)

        scf_drv_minus = ScfRestrictedDriver(comm, ostream)
        scf_drv_minus.update_settings(scf_settings, method_dict_minus)
        scf_result_minus = scf_drv_minus.compute(molecule, basis)

        qrf_minus = QuadraticResponseDriver(comm, ostream)
        qrf_minus.update_settings(qrf_settings, method_settings)
        quad_result_minus = qrf_minus.compute(molecule, basis, scf_result_minus)

        if is_mpi_master():
            beta_plus = -quad_result_plus[(wb, wc)]
            beta_minus = -quad_result_minus[(wb, wc)]
            assert abs(beta_plus.imag) < 1.0e-6
            assert abs(beta_minus.imag) < 1.0e-6

            gamma_0_fd = (beta_plus.real - beta_minus.real) / (2.0 * delta_ef)
            assert abs(gamma_0 - gamma_0_fd) / abs(gamma_0_fd) < 1.0e-5

    @pytest.mark.skipif(is_single_node(), reason="multi-node only")
    def test_lda_crf_fd(self):

        self.run_crf_fd('slda', 'def2-svp', 'zyyz', [0.11, -0.3, 0.05])

    @pytest.mark.skipif(is_single_node(), reason="multi-node only")
    def test_gga_crf_fd(self):

        self.run_crf_fd('pbe0', 'def2-svp', 'zyyz', [0.11, -0.3, 0.05])

    @pytest.mark.skipif(is_single_node(), reason="multi-node only")
    def test_mgga_crf_fd(self):

        self.run_crf_fd('tpssh', 'def2-svp', 'zyyz', [0.11, -0.3, 0.05])
