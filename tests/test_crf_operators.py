from mpi4py import MPI
import pytest

from veloxchem.veloxchemlib import mpi_master
from veloxchem.cubicresponsedriver import CubicResponseDriver
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver


@pytest.mark.solvers
class TestCrfOperators:

    def run_crf(self, basis_set_label, freqs, operators, components,
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

        crf_drv = CubicResponseDriver()
        crf_drv.ostream.mute()

        wb, wc, wd = freqs
        op_a, op_b, op_c, op_d = operators
        a, b, c, d = components

        rsp_settings = {
            'b_frequencies': [wb],
            'c_frequencies': [wc],
            'd_frequencies': [wd],
            'damping': 0.0,
            'a_operator': op_a,
            'b_operator': op_b,
            'c_operator': op_c,
            'd_operator': op_d,
            'a_component': a,
            'b_component': b,
            'c_component': c,
            'd_component': d,
        }
        crf_drv.update_settings(rsp_settings)

        crf_result = crf_drv.compute(molecule, basis, scf_results)

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            for key in ['T4', 'X3', 'A3']:
                assert abs(crf_result[(f'crf_{key}_term', wb, wc, wd)].real -
                           ref_result[key].real) < 0.03
                assert abs(crf_result[(f'crf_{key}_term', wb, wc, wd)].imag -
                           ref_result[key].imag) < 0.03

    def test_svp_crf_general(self):

        wb_wc_wd = (0.0656, 0.0445, 0.0353)

        operator_mapping = {
            'e': 'electric dipole',
            'm': 'angular momentum',
        }

        op_str_list = [
            'eeee',
            'eeem',
            'eeme',
            'eemm',
            'emee',
            'emem',
            'emme',
            'emmm',
            'meee',
            'meem',
            'meme',
            'memm',
            'mmee',
            'mmem',
            'mmme',
            'mmmm',
        ]

        ref_X3_val_list = [
            47.57223563,
            -5.06330192,
            -5.63413078,
            -40.60636029,
            -6.95738248,
            -39.70652823,
            -39.29391410,
            59.20051475,
            -15.64102054,
            49.46594516,
            50.05185396,
            23.91546081,
            51.44459303,
            18.16626025,
            15.67518228,
            -808.48948896,
        ]

        ref_T4_val_list = [
            -15.33001762,
            1.92386955,
            2.07261546,
            7.11835636,
            2.43085424,
            7.08341324,
            7.07262457,
            -9.91835395,
            5.07908441,
            -13.14573255,
            -13.23720436,
            -2.21476691,
            -13.51463078,
            -1.69584455,
            -1.48975542,
            166.41452167,
        ]

        ref_A3_val_list = [
            14.10592854,
            0.32549442,
            0.02774138,
            -14.93636313,
            -0.65995394,
            -14.94644259,
            -14.90023947,
            4.41571969,
            -0.82504380,
            11.77588187,
            11.82374685,
            3.28809770,
            11.86044732,
            0.64417968,
            -0.50271541,
            -234.95640919,
        ]

        for op_str, ref_T4_val, ref_X3_val, ref_A3_val, in zip(
                op_str_list, ref_T4_val_list, ref_X3_val_list, ref_A3_val_list):

            operators = [operator_mapping[x] for x in op_str]

            if op_str.count('m') % 2 == 0:
                rsp_ref_T4_val = ref_T4_val
                rsp_ref_X3_val = ref_X3_val
                rsp_ref_A3_val = ref_A3_val
            else:
                rsp_ref_T4_val = ref_T4_val * 1j
                rsp_ref_X3_val = ref_X3_val * 1j
                rsp_ref_A3_val = ref_A3_val * 1j

            rsp_ref_result = {
                'T4': rsp_ref_T4_val,
                'X3': rsp_ref_X3_val,
                'A3': rsp_ref_A3_val,
            }

            self.run_crf('def2-svp', wb_wc_wd, operators, 'zzzz',
                         rsp_ref_result)
