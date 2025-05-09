from mpi4py import MPI
import pytest

from veloxchem.outputstream import OutputStream
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.cppsolver import ComplexResponse
from veloxchem.quadraticresponsedriver import QuadraticResponseDriver


@pytest.mark.timeconsuming
class TestQrfOperatorsFD:

    def test_qrf_op_fd(self):

        comm = MPI.COMM_WORLD
        ostream = OutputStream(None)

        basis_set_label = 'def2-svp'

        molecule = Molecule.read_xyz_string("""4
            Hydrogen peroxide
            O  -0.65564532 -0.06106286 -0.03621403
            O   0.65564532  0.06106286 -0.03621403
            H  -0.97628735  0.65082652  0.57474201
            H   0.97628735 -0.65082652  0.57474201""")

        basis = MolecularBasis.read(molecule, basis_set_label, ostream=None)

        ops = ['linear momentum', 'electric dipole']

        wb, wc = (0.0656, 0.0)
        wa = -(wb + wc)

        a, b, c = 'zzz'

        # SCF

        scfdrv = ScfRestrictedDriver(comm, ostream)
        scf_result = scfdrv.compute(molecule, basis)

        # finite difference prep

        delta_ef = 1.0e-4
        fd_index = 2  # 'z'

        efield_plus = [0.0, 0.0, 0.0]
        efield_minus = [0.0, 0.0, 0.0]
        efield_plus[fd_index] = delta_ef
        efield_minus[fd_index] = -delta_ef

        scf_drv_plus = ScfRestrictedDriver(comm, ostream)
        scf_drv_plus.electric_field = efield_plus
        scf_result_plus = scf_drv_plus.compute(molecule, basis)

        scf_drv_minus = ScfRestrictedDriver(comm, ostream)
        scf_drv_minus.electric_field = efield_minus
        scf_result_minus = scf_drv_minus.compute(molecule, basis)

        # QRF prep

        qrf = QuadraticResponseDriver(comm, ostream)

        tol = 1.0e-4

        for op_a in ops:
            for op_b in ops:

                op_c = 'electric dipole'

                # finite difference

                lrf_plus = ComplexResponse(comm, ostream)
                lrf_plus.frequencies = [wb]
                lrf_plus.a_operator = op_a
                lrf_plus.b_operator = op_b
                lrf_plus.a_components = a
                lrf_plus.b_components = b
                lrf_plus.damping = 0.0
                lrf_plus.electric_field = efield_plus
                lrf_result_plus = lrf_plus.compute(molecule, basis,
                                                   scf_result_plus)

                lrf_minus = ComplexResponse(comm, ostream)
                lrf_minus.frequencies = [wb]
                lrf_minus.a_operator = op_a
                lrf_minus.b_operator = op_b
                lrf_minus.a_components = a
                lrf_minus.b_components = b
                lrf_minus.damping = 0.0
                lrf_minus.electric_field = efield_minus
                lrf_result_minus = lrf_minus.compute(molecule, basis,
                                                     scf_result_minus)

                rsp_func_plus = lrf_result_plus['response_functions'][(a, b,
                                                                       wb)]
                rsp_func_minus = lrf_result_minus['response_functions'][(a, b,
                                                                         wb)]

                rsp_func_fd = (rsp_func_plus - rsp_func_minus) / (2.0 *
                                                                  delta_ef)

                # QRF

                rsp_settings = {
                    'b_frequencies': [wb],
                    'c_frequencies': [wc],
                    'a_operator': op_a,
                    'b_operator': op_b,
                    'c_operator': op_c,
                    'a_component': a,
                    'b_component': b,
                    'c_component': c,
                    'damping': 0.0,
                }
                qrf.update_settings(rsp_settings)
                qrf_result = qrf.compute(molecule, basis, scf_result)

                assert abs(qrf_result[('qrf', wb, wc)].real -
                           rsp_func_fd.real) < tol
                assert abs(qrf_result[('qrf', wb, wc)].imag -
                           rsp_func_fd.imag) < tol

                # QRF swap_ab

                qrf = QuadraticResponseDriver(comm, ostream)

                rsp_settings = {
                    'b_frequencies': [wa],
                    'c_frequencies': [wc],
                    'a_operator': op_b,
                    'b_operator': op_a,
                    'c_operator': op_c,
                    'a_component': b,
                    'b_component': a,
                    'c_component': c,
                    'damping': 0.0,
                }
                qrf.update_settings(rsp_settings)
                qrf_result = qrf.compute(molecule, basis, scf_result)

                assert abs(qrf_result[('qrf', wa, wc)].real -
                           rsp_func_fd.real) < tol
                assert abs(qrf_result[('qrf', wa, wc)].imag -
                           rsp_func_fd.imag) < tol

                # QRF swap_ac

                qrf = QuadraticResponseDriver(comm, ostream)

                rsp_settings = {
                    'b_frequencies': [wb],
                    'c_frequencies': [wa],
                    'a_operator': op_c,
                    'b_operator': op_b,
                    'c_operator': op_a,
                    'a_component': c,
                    'b_component': b,
                    'c_component': a,
                    'damping': 0.0,
                }
                qrf.update_settings(rsp_settings)
                qrf_result = qrf.compute(molecule, basis, scf_result)

                assert abs(qrf_result[('qrf', wb, wa)].real -
                           rsp_func_fd.real) < tol
                assert abs(qrf_result[('qrf', wb, wa)].imag -
                           rsp_func_fd.imag) < tol

                # QRF swap_bc

                qrf = QuadraticResponseDriver(comm, ostream)

                rsp_settings = {
                    'b_frequencies': [wc],
                    'c_frequencies': [wb],
                    'a_operator': op_a,
                    'b_operator': op_c,
                    'c_operator': op_b,
                    'a_component': a,
                    'b_component': c,
                    'c_component': b,
                    'damping': 0.0,
                }
                qrf.update_settings(rsp_settings)
                qrf_result = qrf.compute(molecule, basis, scf_result)

                assert abs(qrf_result[('qrf', wc, wb)].real -
                           rsp_func_fd.real) < tol
                assert abs(qrf_result[('qrf', wc, wb)].imag -
                           rsp_func_fd.imag) < tol
