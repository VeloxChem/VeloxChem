from copy import deepcopy

from veloxchem.veloxchemlib import parse_xc_func


class TestXCFunctional:

    def test_xcfunctional_deepcopy(self):

        f_a = parse_xc_func('cam-b3lyp')
        f_b = deepcopy(f_a)

        assert f_a.get_frac_exact_exchange() == f_b.get_frac_exact_exchange()
        assert f_a.is_hybrid() == f_b.is_hybrid()
        assert f_a.is_range_separated() == f_b.is_range_separated()
        assert f_a.get_rs_alpha() == f_b.get_rs_alpha()
        assert f_a.get_rs_beta() == f_b.get_rs_beta()
        assert f_a.get_rs_omega() == f_b.get_rs_omega()
        assert f_a.get_func_type() == f_b.get_func_type()
        assert f_a.get_func_label() == f_b.get_func_label()

    def test_xcfunctional_parameters(self):

        f_a = parse_xc_func('blyp')
        assert f_a.get_frac_exact_exchange() == 0.0
        assert not f_a.is_hybrid()

        f_a = parse_xc_func('b3lyp')
        assert f_a.get_frac_exact_exchange() == 0.2
        assert f_a.is_hybrid()
        assert not f_a.is_range_separated()

        f_a = parse_xc_func('cam-b3lyp')
        assert f_a.get_frac_exact_exchange() == 0.65
        assert f_a.get_rs_omega() == 0.33
        assert f_a.is_hybrid()
        assert f_a.is_range_separated()

        f_a = parse_xc_func('cam-b3lyp-100')
        assert f_a.get_frac_exact_exchange() == 1.0
        assert f_a.get_rs_omega() == 0.33
        assert f_a.is_hybrid()
        assert f_a.is_range_separated()
