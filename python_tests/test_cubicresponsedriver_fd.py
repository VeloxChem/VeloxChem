from mpi4py import MPI
import pytest

from veloxchem.veloxchemlib import is_mpi_master
from veloxchem.outputstream import OutputStream
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.quadraticresponsedriver import QuadraticResponseDriver
from veloxchem.cubicresponsedriver import CubicResponseDriver


@pytest.mark.solvers
class TestCrfFD:

    def run_crf_fd(self, xcfun_label):

        comm = MPI.COMM_WORLD
        ostream = OutputStream(None)

        molecule_string = """
        O   0.0   0.0   0.0
        H   0.0   1.4   1.1
        H   0.0  -1.4   1.1
        """
        basis_set_label = 'def2-svp'
        scf_conv_thresh = 1.0e-8
        rsp_conv_thresh = 1.0e-5

        molecule = Molecule.read_str(molecule_string, units='au')
        basis = MolecularBasis.read(molecule, basis_set_label, ostream=None)

        # CR driver

        scf_settings = {'conv_thresh': scf_conv_thresh}
        rsp_settings = {
            'conv_thresh': rsp_conv_thresh,
            'damping': 0,
            'a_components': 'z',
            'b_components': 'y',
            'c_components': 'y',
            'd_components': 'z',
            'b_frequencies': [0.1],
            'c_frequencies': [-0.3],
            'd_frequencies': [0]
        }
        method_settings = {'xcfun': xcfun_label}

        scfdrv = ScfRestrictedDriver(comm, ostream)
        scfdrv.update_settings(scf_settings, method_settings)
        scfdrv.compute(molecule, basis)

        crf = CubicResponseDriver(comm, ostream)
        crf.update_settings(rsp_settings, method_settings)
        crf_results = crf.compute(molecule, basis, scfdrv.scf_tensors)

        if is_mpi_master():
            gamma_zyyz = -crf_results[('gamma', 0.1, -0.3, 0)].real

        # Finite difference

        delta_ef = 1.0e-4
        efield_plus = [0.0, 0.0, delta_ef]
        efield_minus = [0.0, 0.0, -delta_ef]

        method_dict_plus = dict(method_settings)
        method_dict_minus = dict(method_settings)
        method_dict_plus['electric_field'] = efield_plus
        method_dict_minus['electric_field'] = efield_minus

        scf_drv_plus = ScfRestrictedDriver(comm, ostream)
        scf_drv_plus.update_settings(scf_settings, method_dict_plus)
        scf_result_plus = scf_drv_plus.compute(molecule, basis)

        quad_solver_plus = QuadraticResponseDriver(comm, ostream)
        quad_solver_plus.update_settings(rsp_settings, method_settings)
        quad_result_plus = quad_solver_plus.compute(molecule, basis,
                                                    scf_result_plus)

        scf_drv_minus = ScfRestrictedDriver(comm, ostream)
        scf_drv_minus.update_settings(scf_settings, method_dict_minus)
        scf_result_minus = scf_drv_minus.compute(molecule, basis)

        quad_solver_minus = QuadraticResponseDriver(comm, ostream)
        quad_solver_minus.update_settings(rsp_settings, method_settings)
        quad_result_minus = quad_solver_minus.compute(molecule, basis,
                                                      scf_result_minus)

        if is_mpi_master():
            beta_zyy_plus = -quad_result_plus[(0.1, -0.3)].real
            beta_zyy_minus = -quad_result_minus[(0.1, -0.3)].real
            gamma_zyyz_fd = (beta_zyy_plus - beta_zyy_minus) / (2.0 * delta_ef)

            rel_diff = abs(gamma_zyyz - gamma_zyyz_fd) / abs(gamma_zyyz_fd)
            assert rel_diff < 1.0e-5

    def test_lda_crf_fd(self):

        self.run_crf_fd('slda')

    def test_gga_crf_fd(self):

        self.run_crf_fd('pbe0')

    def test_mgga_crf_fd(self):

        self.run_crf_fd('tpssh')
