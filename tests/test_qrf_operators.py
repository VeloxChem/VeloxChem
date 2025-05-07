from mpi4py import MPI
import pytest

from veloxchem.veloxchemlib import mpi_master
from veloxchem.quadraticresponsedriver import QuadraticResponseDriver
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver


@pytest.mark.solvers
class TestQrfOperators:

    def run_qrf(self, basis_set_label, freqs, operators, components,
                ref_result):

        molecule = Molecule.read_xyz_string("""4
            Hydrogen peroxide
            O  -0.65564532 -0.06106286 -0.03621403
            O   0.65564532  0.06106286 -0.03621403
            H  -0.97628735  0.65082652  0.57474201
            H   0.97628735 -0.65082652  0.57474201""")

        basis = MolecularBasis.read(molecule, basis_set_label, ostream=None)

        scf_drv = ScfRestrictedDriver()
        scf_drv.ostream.mute()
        scf_results = scf_drv.compute(molecule, basis)

        qrf_drv = QuadraticResponseDriver()
        qrf_drv.ostream.mute()

        wb, wc = freqs
        op_a, op_b, op_c = operators
        a, b, c = components

        rsp_settings = {
            'b_frequencies': [wb],
            'c_frequencies': [wc],
            'damping': 0.0,
            'a_operator': op_a,
            'b_operator': op_b,
            'c_operator': op_c,
            'a_component': a,
            'b_component': b,
            'c_component': c,
        }
        qrf_drv.update_settings(rsp_settings)
        qrf_result = qrf_drv.compute(molecule, basis, scf_results)

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            assert abs(qrf_result[('qrf', wb, wc)].real -
                       ref_result.real) < 0.01
            assert abs(qrf_result[('qrf', wb, wc)].imag -
                       ref_result.imag) < 0.01

    def test_svp_qrf_eme(self):

        wb_wc = (0.0656, 0.0)

        operators = [
            'electric dipole',
            'angular momentum',
            'electric dipole',
        ]

        self.run_qrf('def2-svp', wb_wc, operators, 'zzz', -1.72680998j)
        self.run_qrf('def2-svp', wb_wc, operators, 'xyz', -3.00108916j)
        self.run_qrf('def2-svp', wb_wc, operators, 'zxy', 1.02519503j)
        self.run_qrf('def2-svp', wb_wc, operators, 'yyx', 0.0j)

    def test_svp_qrf_vmv(self):

        wb_wc = (0.0656, 0.0)

        operators = [
            'linear momentum',
            'angular momentum',
            'linear momentum',
        ]

        self.run_qrf('def2-svp', wb_wc, operators, 'zzz', -0.03225185j)
        self.run_qrf('def2-svp', wb_wc, operators, 'xyz', 0.53669370j)
        self.run_qrf('def2-svp', wb_wc, operators, 'zxy', 0.40995799j)
        self.run_qrf('def2-svp', wb_wc, operators, 'yyx', 0.0j)

    def test_svp_qrf_general(self):

        wb_wc = (0.0656, 0.0445)

        operator_mapping = {
            'e': 'electric dipole',
            'm': 'angular momentum',
        }

        op_str_list = [
            'eee',
            'mee',
            'eme',
            'eem',
            'mem',
            'mme',
            'emm',
            'mmm',
        ]

        ref_val_list = [
            12.74758346,
            -3.05827947,
            -1.82846015,
            -1.24159433,
            20.02376852,
            20.57829175,
            -18.22153615,
            -0.01101933,
        ]

        for op_str, ref_val in zip(op_str_list, ref_val_list):

            operators = [operator_mapping[x] for x in op_str]

            if op_str.count('m') % 2 == 0:
                rsp_ref_val = ref_val
            else:
                rsp_ref_val = ref_val * 1j

            self.run_qrf('def2-svp', wb_wc, operators, 'zzz', rsp_ref_val)
