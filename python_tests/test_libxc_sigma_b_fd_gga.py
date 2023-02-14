import numpy as np

from veloxchem.veloxchemlib import XCNewIntegrator


class TestLibxcSigmaB:

    def check_term(self, tol, label, val, val_fd):

        if abs(val_fd) < tol:
            assert abs(val - val_fd) < tol
        else:
            assert abs(val - val_fd) / abs(val_fd) < tol

    def run_sigma_b_fd(self, xcfun_label, tol):

        delta_h = 1.0e-6

        xc_drv = XCNewIntegrator()

        rho_a, rho_b = 0.17, 0.15
        nabla_a, nabla_b = 0.26, 0.29
        sigma_aa = nabla_a * nabla_a
        sigma_ab = nabla_a * nabla_b
        sigma_bb = nabla_b * nabla_b

        inp = {
            'rho': np.array([rho_a, rho_b]),
            'sigma': np.array([sigma_aa, sigma_ab, sigma_bb]),
        }

        lxc = xc_drv.compute_lxc_for_gga(xcfun_label, inp['rho'], inp['sigma'])
        kxc = xc_drv.compute_kxc_for_gga(xcfun_label, inp['rho'], inp['sigma'])
        fxc = xc_drv.compute_fxc_for_gga(xcfun_label, inp['rho'], inp['sigma'])
        vxc = xc_drv.compute_exc_vxc_for_gga(xcfun_label, inp['rho'],
                                             inp['sigma'])

        inp_plus = dict(inp)
        nabla_b_plus = nabla_b + delta_h
        inp_plus['sigma'][1] = nabla_b_plus * nabla_a
        inp_plus['sigma'][2] = nabla_b_plus * nabla_b_plus

        kxc_plus = xc_drv.compute_kxc_for_gga(xcfun_label, inp_plus['rho'],
                                              inp_plus['sigma'])
        fxc_plus = xc_drv.compute_fxc_for_gga(xcfun_label, inp_plus['rho'],
                                              inp_plus['sigma'])
        vxc_plus = xc_drv.compute_exc_vxc_for_gga(xcfun_label, inp_plus['rho'],
                                                  inp_plus['sigma'])
        exc_plus = {'zk': vxc_plus['exc']}

        inp_minus = dict(inp)
        nabla_b_minus = nabla_b - delta_h
        inp_minus['sigma'][1] = nabla_b_minus * nabla_a
        inp_minus['sigma'][2] = nabla_b_minus * nabla_b_minus

        kxc_minus = xc_drv.compute_kxc_for_gga(xcfun_label, inp_minus['rho'],
                                               inp_minus['sigma'])
        fxc_minus = xc_drv.compute_fxc_for_gga(xcfun_label, inp_minus['rho'],
                                               inp_minus['sigma'])
        vxc_minus = xc_drv.compute_exc_vxc_for_gga(xcfun_label,
                                                   inp_minus['rho'],
                                                   inp_minus['sigma'])
        exc_minus = {'zk': vxc_minus['exc']}

        vsigma_b = vxc['vsigma'][0][2]
        vsigma_c = vxc['vsigma'][0][1]
        e_a_plus = exc_plus['zk'][0][0] * (rho_a + rho_b)
        e_a_minus = exc_minus['zk'][0][0] * (rho_a + rho_b)
        vsigma_b_fd = (e_a_plus - e_a_minus) / (delta_h * 2)
        self.check_term(tol, 'vsigma_b',
                        2 * vsigma_b * nabla_b + vsigma_c * nabla_a,
                        vsigma_b_fd)

        v2sigma2_bb = fxc["v2sigma2"][0][5]
        v2sigma2_cb = fxc["v2sigma2"][0][4]
        vsigma_b_plus = vxc_plus["vsigma"][0][2]
        vsigma_b_minus = vxc_minus["vsigma"][0][2]
        v2sigma2_bb_fd = (vsigma_b_plus - vsigma_b_minus) / (delta_h * 2)
        self.check_term(tol, "v2sigma2_bb",
                   2 * v2sigma2_bb * nabla_b + v2sigma2_cb * nabla_a,
                   v2sigma2_bb_fd)

        v3sigma3_bbb = kxc["v3sigma3"][0][9]
        v3sigma3_cbb = kxc["v3sigma3"][0][8]
        v2sigma2_bb_plus = fxc_plus["v2sigma2"][0][5]
        v2sigma2_bb_minus = fxc_minus["v2sigma2"][0][5]
        v3sigma3_bbb_fd = (v2sigma2_bb_plus - v2sigma2_bb_minus) / (delta_h * 2)
        self.check_term(tol, "v3sigma3_bbb",
                   2 * v3sigma3_bbb * nabla_b + v3sigma3_cbb * nabla_a,
                   v3sigma3_bbb_fd)

        v4sigma4_bbbb = lxc["v4sigma4"][0][14]
        v4sigma4_cbbb = lxc["v4sigma4"][0][13]
        v3sigma3_bbb_plus = kxc_plus["v3sigma3"][0][9]
        v3sigma3_bbb_minus = kxc_minus["v3sigma3"][0][9]
        v4sigma4_bbbb_fd = (v3sigma3_bbb_plus - v3sigma3_bbb_minus) / (delta_h *
                                                                       2)
        self.check_term(tol, "v4sigma4_bbbb",
                   2 * v4sigma4_bbbb * nabla_b + v4sigma4_cbbb * nabla_a,
                   v4sigma4_bbbb_fd)

    def test_sigma_b_fd_blyp(self):

        self.run_sigma_b_fd('blyp', 1.0e-7)

    def test_sigma_b_fd_b3lyp(self):

        self.run_sigma_b_fd('b3lyp', 1.0e-7)

    def test_sigma_b_fd_pbe0(self):

        self.run_sigma_b_fd('pbe0', 1.0e-7)

    def test_sigma_b_fd_bp86(self):

        self.run_sigma_b_fd('bp86', 1.0e-7)
