from mpi4py import MPI
import pytest

from veloxchem.outputstream import OutputStream
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.quadraticresponsedriver import QuadraticResponseDriver
from veloxchem.cubicresponsedriver import CubicResponseDriver


@pytest.mark.timeconsuming
class TestCrfOperatorsFD:

    def test_crf_op_fd(self):

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

        wb, wc, wd = (0.0656, 0.0445, 0.0)
        wa = -(wb + wc + wd)

        a, b, c, d = 'zzzz'

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

        # CRF prep

        crf = CubicResponseDriver(comm, ostream)

        tol = 1.0e-3

        for op_a in ops:
            for op_b in ops:
                for op_c in ops:

                    op_d = 'electric dipole'

                    # finite difference

                    qrf_plus = QuadraticResponseDriver(comm, ostream)
                    qrf_plus.b_frequencies = [wb]
                    qrf_plus.c_frequencies = [wc]
                    qrf_plus.a_operator = op_a
                    qrf_plus.b_operator = op_b
                    qrf_plus.c_operator = op_c
                    qrf_plus.a_component = a
                    qrf_plus.b_component = b
                    qrf_plus.c_component = c
                    qrf_plus.damping = 0
                    qrf_plus.electric_field = efield_plus
                    qrf_result_plus = qrf_plus.compute(molecule, basis,
                                                       scf_result_plus)

                    qrf_minus = QuadraticResponseDriver(comm, ostream)
                    qrf_minus.b_frequencies = [wb]
                    qrf_minus.c_frequencies = [wc]
                    qrf_minus.a_operator = op_a
                    qrf_minus.b_operator = op_b
                    qrf_minus.c_operator = op_c
                    qrf_minus.a_component = a
                    qrf_minus.b_component = b
                    qrf_minus.c_component = c
                    qrf_minus.damping = 0
                    qrf_minus.electric_field = efield_minus
                    qrf_result_minus = qrf_minus.compute(
                        molecule, basis, scf_result_minus)

                    rsp_func_plus = qrf_result_plus[('qrf', wb, wc)]
                    rsp_func_minus = qrf_result_minus[('qrf', wb, wc)]

                    rsp_func_fd = (rsp_func_plus - rsp_func_minus) / (2.0 *
                                                                      delta_ef)

                    # CRF

                    rsp_settings = {
                        'b_frequencies': [wb],
                        'c_frequencies': [wc],
                        'd_frequencies': [wd],
                        'a_operator': op_a,
                        'b_operator': op_b,
                        'c_operator': op_c,
                        'd_operator': op_d,
                        'a_component': a,
                        'b_component': b,
                        'c_component': c,
                        'd_component': d,
                        'damping': 0.0,
                    }
                    crf.update_settings(rsp_settings)
                    crf_result = crf.compute(molecule, basis, scf_result)

                    assert abs(crf_result[('crf', wb, wc, wd)].real -
                               rsp_func_fd.real) < tol
                    assert abs(crf_result[('crf', wb, wc, wd)].imag -
                               rsp_func_fd.imag) < tol

                    # CRF swap_ab

                    crf = CubicResponseDriver(comm, ostream)

                    rsp_settings = {
                        'b_frequencies': [wa],
                        'c_frequencies': [wc],
                        'd_frequencies': [wd],
                        'a_operator': op_b,
                        'b_operator': op_a,
                        'c_operator': op_c,
                        'd_operator': op_d,
                        'a_component': b,
                        'b_component': a,
                        'c_component': c,
                        'd_component': d,
                        'damping': 0.0,
                    }
                    crf.update_settings(rsp_settings)
                    crf_result = crf.compute(molecule, basis, scf_result)

                    assert abs(crf_result[('crf', wa, wc, wd)].real -
                               rsp_func_fd.real) < tol
                    assert abs(crf_result[('crf', wa, wc, wd)].imag -
                               rsp_func_fd.imag) < tol

                    # CRF swap_ac

                    crf = CubicResponseDriver(comm, ostream)

                    rsp_settings = {
                        'b_frequencies': [wb],
                        'c_frequencies': [wa],
                        'd_frequencies': [wd],
                        'a_operator': op_c,
                        'b_operator': op_b,
                        'c_operator': op_a,
                        'd_operator': op_d,
                        'a_component': c,
                        'b_component': b,
                        'c_component': a,
                        'd_component': d,
                        'damping': 0.0,
                    }
                    crf.update_settings(rsp_settings)
                    crf_result = crf.compute(molecule, basis, scf_result)

                    assert abs(crf_result[('crf', wb, wa, wd)].real -
                               rsp_func_fd.real) < tol
                    assert abs(crf_result[('crf', wb, wa, wd)].imag -
                               rsp_func_fd.imag) < tol

                    # CRF swap_ad

                    crf = CubicResponseDriver(comm, ostream)

                    rsp_settings = {
                        'b_frequencies': [wb],
                        'c_frequencies': [wc],
                        'd_frequencies': [wa],
                        'a_operator': op_d,
                        'b_operator': op_b,
                        'c_operator': op_c,
                        'd_operator': op_a,
                        'a_component': d,
                        'b_component': b,
                        'c_component': c,
                        'd_component': a,
                        'damping': 0.0,
                    }
                    crf.update_settings(rsp_settings)
                    crf_result = crf.compute(molecule, basis, scf_result)

                    assert abs(crf_result[('crf', wb, wc, wa)].real -
                               rsp_func_fd.real) < tol
                    assert abs(crf_result[('crf', wb, wc, wa)].imag -
                               rsp_func_fd.imag) < tol

                    # CRF swap_bc

                    crf = CubicResponseDriver(comm, ostream)

                    rsp_settings = {
                        'b_frequencies': [wc],
                        'c_frequencies': [wb],
                        'd_frequencies': [wd],
                        'a_operator': op_a,
                        'b_operator': op_c,
                        'c_operator': op_b,
                        'd_operator': op_d,
                        'a_component': a,
                        'b_component': c,
                        'c_component': b,
                        'd_component': d,
                        'damping': 0.0,
                    }
                    crf.update_settings(rsp_settings)
                    crf_result = crf.compute(molecule, basis, scf_result)

                    assert abs(crf_result[('crf', wc, wb, wd)].real -
                               rsp_func_fd.real) < tol
                    assert abs(crf_result[('crf', wc, wb, wd)].imag -
                               rsp_func_fd.imag) < tol

                    # CRF swap_bd

                    crf = CubicResponseDriver(comm, ostream)

                    rsp_settings = {
                        'b_frequencies': [wd],
                        'c_frequencies': [wc],
                        'd_frequencies': [wb],
                        'a_operator': op_a,
                        'b_operator': op_d,
                        'c_operator': op_c,
                        'd_operator': op_b,
                        'a_component': a,
                        'b_component': d,
                        'c_component': c,
                        'd_component': b,
                        'damping': 0.0,
                    }
                    crf.update_settings(rsp_settings)
                    crf_result = crf.compute(molecule, basis, scf_result)

                    assert abs(crf_result[('crf', wd, wc, wb)].real -
                               rsp_func_fd.real) < tol
                    assert abs(crf_result[('crf', wd, wc, wb)].imag -
                               rsp_func_fd.imag) < tol

                    # CRF swap_cd

                    crf = CubicResponseDriver(comm, ostream)

                    rsp_settings = {
                        'b_frequencies': [wb],
                        'c_frequencies': [wd],
                        'd_frequencies': [wc],
                        'a_operator': op_a,
                        'b_operator': op_b,
                        'c_operator': op_d,
                        'd_operator': op_c,
                        'a_component': a,
                        'b_component': b,
                        'c_component': d,
                        'd_component': c,
                        'damping': 0.0,
                    }
                    crf.update_settings(rsp_settings)
                    crf_result = crf.compute(molecule, basis, scf_result)

                    assert abs(crf_result[('crf', wb, wd, wc)].real -
                               rsp_func_fd.real) < tol
                    assert abs(crf_result[('crf', wb, wd, wc)].imag -
                               rsp_func_fd.imag) < tol
