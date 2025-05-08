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
            crf_terms = crf_result['crf_terms']

            for key in ref_result:
                val = crf_terms[(f'crf_{key}_term', wb, wc, wd)]
                ref = ref_result[key]

                real_diff = abs(val.real - ref.real)
                if ref.real != 0.0:
                    real_diff_rel = abs(val.real / ref.real - 1.0)
                else:
                    real_diff_rel = real_diff
                assert real_diff < 0.01 or real_diff_rel < 0.0075

                imag_diff = abs(val.imag - ref.imag)
                if ref.imag != 0.0:
                    imag_diff_rel = abs(val.imag / ref.imag - 1.0)
                else:
                    imag_diff_rel = imag_diff
                assert imag_diff < 0.01 or imag_diff_rel < 0.0075

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

        ref_E3_val_list = [
            41.67203860,
            -5.06448766,
            -5.74020667,
            -13.02657147,
            -7.27793590,
            -12.89558579,
            -12.64484193,
            -5.44237093,
            -15.53443126,
            37.54122754,
            38.40397450,
            -7.08576100,
            40.58914831,
            -7.14492322,
            -7.24720010,
            111.26824095,
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

        ref_X2_val_list = [
            -114.85879508,
            22.28385162,
            23.54243572,
            18.97427492,
            26.50052407,
            17.71902413,
            17.08431540,
            -44.47939901,
            62.35822578,
            -163.20347830,
            -168.45666061,
            -9.96838888,
            -181.35236975,
            -6.12087456,
            -4.43913877,
            591.37780309,
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

        ref_A2_val_list = [
            -84.72521225,
            -4.74704615,
            -2.13983171,
            112.88885116,
            3.75304949,
            113.01886419,
            111.98356759,
            -7.53991587,
            3.34010323,
            -24.99025713,
            -25.49746772,
            -10.10244665,
            -26.60018023,
            -5.40880979,
            -3.31133549,
            495.69636664,
        ]

        for (op_str, ref_T4_val, ref_E3_val, ref_X3_val, ref_X2_val, ref_A3_val,
             ref_A2_val) in zip(op_str_list, ref_T4_val_list, ref_E3_val_list,
                                ref_X3_val_list, ref_X2_val_list,
                                ref_A3_val_list, ref_A2_val_list):

            operators = [operator_mapping[x] for x in op_str]

            if op_str.count('m') % 2 == 0:
                rsp_ref_T4_val = ref_T4_val
                rsp_ref_E3_val = ref_E3_val
                rsp_ref_X3_val = ref_X3_val
                rsp_ref_X2_val = ref_X2_val
                rsp_ref_A3_val = ref_A3_val
                rsp_ref_A2_val = ref_A2_val
            else:
                rsp_ref_T4_val = ref_T4_val * 1j
                rsp_ref_E3_val = ref_E3_val * 1j
                rsp_ref_X3_val = ref_X3_val * 1j
                rsp_ref_X2_val = ref_X2_val * 1j
                rsp_ref_A3_val = ref_A3_val * 1j
                rsp_ref_A2_val = ref_A2_val * 1j

            rsp_ref_result = {
                'T4': rsp_ref_T4_val,
                'E3': rsp_ref_E3_val,
                'X3': rsp_ref_X3_val,
                'X2': rsp_ref_X2_val,
                'A3': rsp_ref_A3_val,
                'A2': rsp_ref_A2_val,
            }

            self.run_crf('def2-svp', wb_wc_wd, operators, 'zzzz',
                         rsp_ref_result)
