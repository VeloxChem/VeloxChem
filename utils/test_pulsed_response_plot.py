from pathlib import Path
import unittest
import h5py

import pulsed_response_plot as prp


class TestPulsedRspPlot(unittest.TestCase):

    def test_mirror_and_zeropad(self):

        # Simple test using real non-negative integers
        probe_arr_a = [0, 1, 2, 3]
        zp_test_a = [0, 0, 0]
        exp_ans_a = [0, 0, 0, 3, 2, 1, 0, 1, 2, 3, 0, 0, 0]

        self.assertEqual(prp.mirror_and_zeropad(probe_arr_a, zp_test_a),
                         exp_ans_a)

        # More elaborate test using complex numbers
        probe_arr_b = [0.0, 4.0 + 3.0j, 1.5, 1.0 - 5.0j, -2.0j]
        zp_test_b = [0.0, 0.0]
        exp_ans_b = [
            0.0, 0.0, 2.0j, 1.0 + 5.0j, 1.5, 4.0 - 3.0j, 0.0, 4.0 + 3.0j, 1.5,
            1.0 - 5.0j, -2.0j, 0.0, 0.0
        ]

        self.assertEqual(prp.mirror_and_zeropad(probe_arr_b, zp_test_b),
                         exp_ans_b)

    def test_get_txt_exp(self):

        # Simple test
        self.assertEqual(prp.get_txt_exp(3), r'$\times 10^{3}$')

        # 10**0 results should return no exponent string
        self.assertEqual(prp.get_txt_exp(0), '')

        # A couple of extra tests for the -1 and +1 cases
        self.assertEqual(prp.get_txt_exp(-1), r'$\times 10^{-1}$')
        self.assertEqual(prp.get_txt_exp(1), r'$\times 10^{1}$')

        # All of the following should raise errors since argument was not
        # an integer
        with self.assertRaises(ValueError):
            prp.get_txt_exp('abcd')

        with self.assertRaises(ValueError):
            prp.get_txt_exp('0.5')

        with self.assertRaises(ValueError):
            prp.get_txt_exp('-1.6')

        with self.assertRaises(ValueError):
            prp.get_txt_exp('2.0')

        with self.assertRaises(ValueError):
            prp.get_txt_exp('0.0')

    def test_get_plot_ceil(self):

        # 0.05 is 5e-2
        r_t, ex_t, txt_t = prp.get_plot_ceil(0.05)
        self.assertAlmostEqual(r_t, 5.0, places=12)
        self.assertEqual(ex_t, -2)
        self.assertEqual(txt_t, r'$\times 10^{-2}$')

        # 6546 is 6.546e+3
        r_t, ex_t, txt_t = prp.get_plot_ceil(6546)
        self.assertAlmostEqual(r_t, 6.546, places=12)
        self.assertEqual(ex_t, 3)
        self.assertEqual(txt_t, r'$\times 10^{3}$')

        # 3.0 is 3.0e0
        r_t, ex_t, txt_t = prp.get_plot_ceil(3.0)
        self.assertAlmostEqual(r_t, 3.0, places=12)
        self.assertEqual(ex_t, 0)
        self.assertEqual(txt_t, '')

        # 9.9999 is 9.9999e0
        r_t, ex_t, txt_t = prp.get_plot_ceil(9.9999)
        self.assertAlmostEqual(r_t, 9.9999, places=12)
        self.assertEqual(ex_t, 0)
        self.assertEqual(txt_t, '')

        # 10.0 is 1.0e+1
        r_t, ex_t, txt_t = prp.get_plot_ceil(10.0)
        self.assertAlmostEqual(r_t, 1.0, places=12)
        self.assertEqual(ex_t, 1)
        self.assertEqual(txt_t, r'$\times 10^{1}$')

        # 10.00001 is 1.000001e+1
        r_t, ex_t, txt_t = prp.get_plot_ceil(10.00001)
        self.assertAlmostEqual(r_t, 1.000001, places=12)
        self.assertEqual(ex_t, 1)
        self.assertEqual(txt_t, r'$\times 10^{1}$')

        # 0.999 is 9.99e-1
        r_t, ex_t, txt_t = prp.get_plot_ceil(0.999)
        self.assertAlmostEqual(r_t, 9.99, places=12)
        self.assertEqual(ex_t, -1)
        self.assertEqual(txt_t, r'$\times 10^{-1}$')

        # 0.100001 is 1.00001e-1
        r_t, ex_t, txt_t = prp.get_plot_ceil(0.100001)
        self.assertAlmostEqual(r_t, 1.00001, places=12)
        self.assertEqual(ex_t, -1)
        self.assertEqual(txt_t, r'$\times 10^{-1}$')

        # 1.0 is 1.0e0
        r_t, ex_t, txt_t = prp.get_plot_ceil(1.0)
        self.assertAlmostEqual(r_t, 1.0, places=12)
        self.assertEqual(ex_t, 0)
        self.assertEqual(txt_t, '')

    def test_plot_pulsed_response(self):

        here = Path(__file__).parent
        ref_file = str(here / 'test_data' /
                       'water_near_monochr_vlx_outp_ref.h5')
        res_ref_file = str(here / 'test_data' /
                           'water_near_monochr_plot_result_ref.h5')

        # "Bootstrap-style" test using water data under assumption that this
        # data was correct
        (xfreqs, alpha_real_plot, alpha_imag_plot, field_w_real_plot,
         field_w_env_plot, t, dipmom_t_real_plot,
         field_t_real_plot) = prp.plot_pulsed_response(
             ['dummy', ref_file, 'save=n'], is_this_a_test=True)

        hf = h5py.File(res_ref_file, 'r')

        xfreqs_ref = hf.get('xfreqs')[()]
        self.assertEqual(xfreqs.size, xfreqs_ref.size)
        for i in range(xfreqs.size):
            self.assertAlmostEqual(xfreqs[i], xfreqs_ref[i], places=12)

        alpha_real_plot_ref = hf.get('alpha_real_plot')[()]
        self.assertEqual(alpha_real_plot.size, alpha_real_plot_ref.size)
        for i in range(alpha_real_plot.size):
            self.assertAlmostEqual(alpha_real_plot[i],
                                   alpha_real_plot_ref[i],
                                   places=12)

        alpha_imag_plot_ref = hf.get('alpha_imag_plot')[()]
        self.assertEqual(alpha_imag_plot.size, alpha_imag_plot_ref.size)
        for i in range(alpha_imag_plot.size):
            self.assertAlmostEqual(alpha_imag_plot[i],
                                   alpha_imag_plot_ref[i],
                                   places=12)

        field_w_real_plot_ref = hf.get('field_w_real_plot')[()]
        self.assertEqual(field_w_real_plot.size, field_w_real_plot_ref.size)
        for i in range(field_w_real_plot.size):
            self.assertAlmostEqual(field_w_real_plot[i],
                                   field_w_real_plot_ref[i],
                                   places=12)

        field_w_env_plot_ref = hf.get('field_w_env_plot')[()]
        self.assertEqual(field_w_env_plot.size, field_w_env_plot_ref.size)
        for i in range(field_w_env_plot.size):
            self.assertAlmostEqual(field_w_env_plot[i],
                                   field_w_env_plot_ref[i],
                                   places=12)

        t_ref = hf.get('t')[()]
        self.assertEqual(t.size, t_ref.size)
        for i in range(t.size):
            self.assertAlmostEqual(t[i], t_ref[i], places=12)

        dipmom_t_real_plot_ref = hf.get('dipmom_t_real_plot')[()]
        self.assertEqual(dipmom_t_real_plot.size, dipmom_t_real_plot_ref.size)
        for i in range(dipmom_t_real_plot.size):
            self.assertAlmostEqual(dipmom_t_real_plot[i],
                                   dipmom_t_real_plot_ref[i],
                                   places=12)

        field_t_real_plot_ref = hf.get('field_t_real_plot')[()]
        self.assertEqual(field_t_real_plot.size, field_t_real_plot_ref.size)
        for i in range(field_t_real_plot.size):
            self.assertAlmostEqual(field_t_real_plot[i],
                                   field_t_real_plot_ref[i],
                                   places=12)

        # To do:
        # - Add tests, data for various command-line arguments concerning
        # plotting


if __name__ == "__main__":

    unittest.main()
