from mpi4py import MPI
import numpy as np
import pytest

from veloxchem.veloxchemlib import mpi_master
from veloxchem.outputstream import OutputStream
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.quadraticresponsedriver import QuadraticResponseDriver
from veloxchem.cubicresponsedriver import CubicResponseDriver


@pytest.mark.timeconsuming
class TestCrfFD:

    def run_crf_fd(self, xcfun_label, basis_set_label, components, freqs,
                   point_charges=None):

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

        if point_charges is not None:
            # rerun the SCF with the point-charge environment; the response
            # drivers then inherit the environment through the SCF reference
            # (orbitals and Fock matrix), since the static point-charge
            # potential does not enter the response equations directly
            scf_result_vacuum = scf_result
            scfdrv.point_charges = point_charges
            scf_result = scfdrv.compute(molecule, basis)
            if MPI.COMM_WORLD.Get_rank() == mpi_master():
                assert abs(scf_result['scf_energy'] -
                           scf_result_vacuum['scf_energy']) > 1.0e-3

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

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            gamma = -crf_result[('crf', wb, wc, wd)]

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

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            gamma_perm = -crf_result[('crf', wb, wd, wc)]
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

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            gamma_conj = -crf_result[('crf', -wb, -wc, -wd)]
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

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            gamma_0 = -crf_result[('crf', wb, wc, 0)].real
            assert abs(-crf_result[('crf', wb, wc, 0)].imag) < 1.0e-6

        if point_charges is not None:
            # check that the point-charge environment changes the response
            crf_vacuum = CubicResponseDriver(comm, ostream)
            crf_vacuum.update_settings(rsp_settings, method_settings)
            crf_result_vacuum = crf_vacuum.compute(molecule, basis,
                                                   scf_result_vacuum)
            if MPI.COMM_WORLD.Get_rank() == mpi_master():
                gamma_0_vacuum = -crf_result_vacuum[('crf', wb, wc, 0)].real
                assert abs(-crf_result_vacuum[('crf', wb, wc,
                                               0)].imag) < 1.0e-6
                assert abs(gamma_0 - gamma_0_vacuum) > 1.0e-2 * abs(
                    gamma_0_vacuum)

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
        if point_charges is not None:
            scf_drv_plus.point_charges = point_charges
        scf_result_plus = scf_drv_plus.compute(molecule, basis)

        qrf_plus = QuadraticResponseDriver(comm, ostream)
        qrf_plus.update_settings(qrf_settings, method_settings)
        quad_result_plus = qrf_plus.compute(molecule, basis, scf_result_plus)

        scf_drv_minus = ScfRestrictedDriver(comm, ostream)
        scf_drv_minus.update_settings(scf_settings, method_dict_minus)
        if point_charges is not None:
            scf_drv_minus.point_charges = point_charges
        scf_result_minus = scf_drv_minus.compute(molecule, basis)

        qrf_minus = QuadraticResponseDriver(comm, ostream)
        qrf_minus.update_settings(qrf_settings, method_settings)
        quad_result_minus = qrf_minus.compute(molecule, basis, scf_result_minus)

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            beta_plus = quad_result_plus[('qrf', wb, wc)]
            beta_minus = quad_result_minus[('qrf', wb, wc)]
            assert abs(beta_plus.imag) < 1.0e-6
            assert abs(beta_minus.imag) < 1.0e-6

            gamma_0_fd = (beta_plus.real - beta_minus.real) / (2.0 * delta_ef)
            assert abs(gamma_0 - gamma_0_fd) / abs(gamma_0_fd) < 1.0e-5

    def test_hf_crf_fd(self):

        self.run_crf_fd('hf', 'def2-tzvp', 'zyyz', [0.11, -0.3, 0.05])

    def test_lda_crf_fd(self):

        self.run_crf_fd('slda', 'def2-svp', 'zyyz', [0.11, -0.3, 0.05])

    def test_gga_hyb_crf_fd(self):

        self.run_crf_fd('b3lyp', 'def2-tzvp', 'zyyz', [0.11, -0.3, 0.05])

    def test_gga_hyb_pc_crf_fd(self):

        # point-charge environment (positions in bohr, charges in e); strong
        # enough to clearly change the response vs. vacuum, but weak enough
        # that the finite-difference test retains its accuracy
        point_charges = np.array([
            # x, y, z, q, sigma, epsilon
            [3.0, 2.0, 1.0, 0.10, 0.0, 0.0],
            [-3.0, -2.0, -1.0, -0.10, 0.0, 0.0],
            [0.0, 0.0, 3.5, 0.06, 0.0, 0.0],
            [2.5, -2.5, 0.0, -0.06, 0.0, 0.0],
        ]).T

        self.run_crf_fd('b3lyp', 'def2-tzvp', 'zyyz', [0.11, -0.3, 0.05],
                        point_charges)

    def test_gga_rsh_crf_fd(self):

        self.run_crf_fd('cam-b3lyp', 'def2-svp', 'zyyz', [0.11, -0.3, 0.05])

    def test_mgga_crf_fd(self):

        self.run_crf_fd('m06-2x', 'def2-svp', 'zyyz', [0.11, -0.3, 0.05])
