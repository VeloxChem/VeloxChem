import numpy as np

from veloxchem.veloxchemlib import XCNewIntegrator


class TestLibxcSigmaA:

    def check_term(self, tol, label, val, val_fd):

        if abs(val_fd) < tol:
            assert abs(val - val_fd) < tol
        else:
            assert abs(val - val_fd) / abs(val_fd) < tol

    def run_sigma_a_fd(self, xcfun_label, tol):

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
        nabla_a_plus = nabla_a + delta_h
        inp_plus['sigma'][0] = nabla_a_plus * nabla_a_plus
        inp_plus['sigma'][1] = nabla_a_plus * nabla_b

        kxc_plus = xc_drv.compute_kxc_for_gga(xcfun_label, inp_plus['rho'],
                                              inp_plus['sigma'])
        fxc_plus = xc_drv.compute_fxc_for_gga(xcfun_label, inp_plus['rho'],
                                              inp_plus['sigma'])
        vxc_plus = xc_drv.compute_exc_vxc_for_gga(xcfun_label, inp_plus['rho'],
                                                  inp_plus['sigma'])
        exc_plus = {'zk': vxc_plus['exc']}

        inp_minus = dict(inp)
        nabla_a_minus = nabla_a - delta_h
        inp_minus['sigma'][0] = nabla_a_minus * nabla_a_minus
        inp_minus['sigma'][1] = nabla_a_minus * nabla_b

        kxc_minus = xc_drv.compute_kxc_for_gga(xcfun_label, inp_minus['rho'],
                                               inp_minus['sigma'])
        fxc_minus = xc_drv.compute_fxc_for_gga(xcfun_label, inp_minus['rho'],
                                               inp_minus['sigma'])
        vxc_minus = xc_drv.compute_exc_vxc_for_gga(xcfun_label,
                                                   inp_minus['rho'],
                                                   inp_minus['sigma'])
        exc_minus = {'zk': vxc_minus['exc']}

        vsigma_a = vxc['vsigma'][0][0]
        vsigma_c = vxc['vsigma'][0][1]
        e_a_plus = exc_plus['zk'][0][0] * (rho_a + rho_b)
        e_a_minus = exc_minus['zk'][0][0] * (rho_a + rho_b)
        vsigma_a_fd = (e_a_plus - e_a_minus) / (delta_h * 2)
        self.check_term(tol, 'vsigma_a',
                        2 * vsigma_a * nabla_a + vsigma_c * nabla_b,
                        vsigma_a_fd)

        v2sigma2_aa = fxc["v2sigma2"][0][0]
        v2sigma2_ac = fxc["v2sigma2"][0][1]
        vsigma_a_plus = vxc_plus["vsigma"][0][0]
        vsigma_a_minus = vxc_minus["vsigma"][0][0]
        v2sigma2_aa_fd = (vsigma_a_plus - vsigma_a_minus) / (delta_h * 2)
        self.check_term(tol, "v2sigma2_aa",
                        2 * v2sigma2_aa * nabla_a + v2sigma2_ac * nabla_b,
                        v2sigma2_aa_fd)

        v2sigma2_ac = fxc["v2sigma2"][0][1]
        v2sigma2_cc = fxc["v2sigma2"][0][3]
        vsigma_c_plus = vxc_plus["vsigma"][0][1]
        vsigma_c_minus = vxc_minus["vsigma"][0][1]
        v2sigma2_ac_fd = (vsigma_c_plus - vsigma_c_minus) / (delta_h * 2)
        self.check_term(tol, "v2sigma2_ac",
                        2 * v2sigma2_ac * nabla_a + v2sigma2_cc * nabla_b,
                        v2sigma2_ac_fd)

        v2sigma2_ab = fxc["v2sigma2"][0][2]
        v2sigma2_cb = fxc["v2sigma2"][0][4]
        vsigma_b_plus = vxc_plus["vsigma"][0][2]
        vsigma_b_minus = vxc_minus["vsigma"][0][2]
        v2sigma2_ab_fd = (vsigma_b_plus - vsigma_b_minus) / (delta_h * 2)
        self.check_term(tol, "v2sigma2_ab",
                        2 * v2sigma2_ab * nabla_a + v2sigma2_cb * nabla_b,
                        v2sigma2_ab_fd)

        v3sigma3_aaa = kxc["v3sigma3"][0][0]
        v3sigma3_aac = kxc["v3sigma3"][0][1]
        v2sigma2_aa_plus = fxc_plus["v2sigma2"][0][0]
        v2sigma2_aa_minus = fxc_minus["v2sigma2"][0][0]
        v3sigma3_aaa_fd = (v2sigma2_aa_plus - v2sigma2_aa_minus) / (delta_h * 2)
        self.check_term(tol, "v3sigma3_aaa",
                        2 * v3sigma3_aaa * nabla_a + v3sigma3_aac * nabla_b,
                        v3sigma3_aaa_fd)

        v3sigma3_aac = kxc["v3sigma3"][0][1]
        v3sigma3_acc = kxc["v3sigma3"][0][3]
        v2sigma2_ac_plus = fxc_plus["v2sigma2"][0][1]
        v2sigma2_ac_minus = fxc_minus["v2sigma2"][0][1]
        v3sigma3_aac_fd = (v2sigma2_ac_plus - v2sigma2_ac_minus) / (delta_h * 2)
        self.check_term(tol, "v3sigma3_aac",
                        2 * v3sigma3_aac * nabla_a + v3sigma3_acc * nabla_b,
                        v3sigma3_aac_fd)

        v3sigma3_aab = kxc["v3sigma3"][0][2]
        v3sigma3_acb = kxc["v3sigma3"][0][4]
        v2sigma2_ab_plus = fxc_plus["v2sigma2"][0][2]
        v2sigma2_ab_minus = fxc_minus["v2sigma2"][0][2]
        v3sigma3_aab_fd = (v2sigma2_ab_plus - v2sigma2_ab_minus) / (delta_h * 2)
        self.check_term(tol, "v3sigma3_aab",
                        2 * v3sigma3_aab * nabla_a + v3sigma3_acb * nabla_b,
                        v3sigma3_aab_fd)

        v3sigma3_acc = kxc["v3sigma3"][0][3]
        v3sigma3_ccc = kxc["v3sigma3"][0][6]
        v2sigma2_cc_plus = fxc_plus["v2sigma2"][0][3]
        v2sigma2_cc_minus = fxc_minus["v2sigma2"][0][3]
        v3sigma3_acc_fd = (v2sigma2_cc_plus - v2sigma2_cc_minus) / (delta_h * 2)
        self.check_term(tol, "v3sigma3_acc",
                        2 * v3sigma3_acc * nabla_a + v3sigma3_ccc * nabla_b,
                        v3sigma3_acc_fd)

        v3sigma3_acb = kxc["v3sigma3"][0][4]
        v3sigma3_ccb = kxc["v3sigma3"][0][7]
        v2sigma2_cb_plus = fxc_plus["v2sigma2"][0][4]
        v2sigma2_cb_minus = fxc_minus["v2sigma2"][0][4]
        v3sigma3_acb_fd = (v2sigma2_cb_plus - v2sigma2_cb_minus) / (delta_h * 2)
        self.check_term(tol, "v3sigma3_acb",
                        2 * v3sigma3_acb * nabla_a + v3sigma3_ccb * nabla_b,
                        v3sigma3_acb_fd)

        v3sigma3_abb = kxc["v3sigma3"][0][5]
        v3sigma3_cbb = kxc["v3sigma3"][0][8]
        v2sigma2_bb_plus = fxc_plus["v2sigma2"][0][5]
        v2sigma2_bb_minus = fxc_minus["v2sigma2"][0][5]
        v3sigma3_abb_fd = (v2sigma2_bb_plus - v2sigma2_bb_minus) / (delta_h * 2)
        self.check_term(tol, "v3sigma3_abb",
                        2 * v3sigma3_abb * nabla_a + v3sigma3_cbb * nabla_b,
                        v3sigma3_abb_fd)

        v4sigma4_aaaa = lxc["v4sigma4"][0][0]
        v4sigma4_aaac = lxc["v4sigma4"][0][1]
        v3sigma3_aaa_plus = kxc_plus["v3sigma3"][0][0]
        v3sigma3_aaa_minus = kxc_minus["v3sigma3"][0][0]
        v4sigma4_aaaa_fd = (v3sigma3_aaa_plus - v3sigma3_aaa_minus) / (delta_h *
                                                                       2)
        self.check_term(tol, "v4sigma4_aaaa",
                        2 * v4sigma4_aaaa * nabla_a + v4sigma4_aaac * nabla_b,
                        v4sigma4_aaaa_fd)

        v4sigma4_aaac = lxc["v4sigma4"][0][1]
        v4sigma4_aacc = lxc["v4sigma4"][0][3]
        v3sigma3_aac_plus = kxc_plus["v3sigma3"][0][1]
        v3sigma3_aac_minus = kxc_minus["v3sigma3"][0][1]
        v4sigma4_aaac_fd = (v3sigma3_aac_plus - v3sigma3_aac_minus) / (delta_h *
                                                                       2)
        self.check_term(tol, "v4sigma4_aaac",
                        2 * v4sigma4_aaac * nabla_a + v4sigma4_aacc * nabla_b,
                        v4sigma4_aaac_fd)

        v4sigma4_aaab = lxc["v4sigma4"][0][2]
        v4sigma4_aacb = lxc["v4sigma4"][0][4]
        v3sigma3_aab_plus = kxc_plus["v3sigma3"][0][2]
        v3sigma3_aab_minus = kxc_minus["v3sigma3"][0][2]
        v4sigma4_aaab_fd = (v3sigma3_aab_plus - v3sigma3_aab_minus) / (delta_h *
                                                                       2)
        self.check_term(tol, "v4sigma4_aaab",
                        2 * v4sigma4_aaab * nabla_a + v4sigma4_aacb * nabla_b,
                        v4sigma4_aaab_fd)

        v4sigma4_aacc = lxc["v4sigma4"][0][3]
        v4sigma4_accc = lxc["v4sigma4"][0][6]
        v3sigma3_acc_plus = kxc_plus["v3sigma3"][0][3]
        v3sigma3_acc_minus = kxc_minus["v3sigma3"][0][3]
        v4sigma4_aacc_fd = (v3sigma3_acc_plus - v3sigma3_acc_minus) / (delta_h *
                                                                       2)
        self.check_term(tol, "v4sigma4_aacc",
                        2 * v4sigma4_aacc * nabla_a + v4sigma4_accc * nabla_b,
                        v4sigma4_aacc_fd)

        v4sigma4_aacb = lxc["v4sigma4"][0][4]
        v4sigma4_accb = lxc["v4sigma4"][0][7]
        v3sigma3_acb_plus = kxc_plus["v3sigma3"][0][4]
        v3sigma3_acb_minus = kxc_minus["v3sigma3"][0][4]
        v4sigma4_aacb_fd = (v3sigma3_acb_plus - v3sigma3_acb_minus) / (delta_h *
                                                                       2)
        self.check_term(tol, "v4sigma4_aacb",
                        2 * v4sigma4_aacb * nabla_a + v4sigma4_accb * nabla_b,
                        v4sigma4_aacb_fd)

        v4sigma4_aabb = lxc["v4sigma4"][0][5]
        v4sigma4_acbb = lxc["v4sigma4"][0][8]
        v3sigma3_abb_plus = kxc_plus["v3sigma3"][0][5]
        v3sigma3_abb_minus = kxc_minus["v3sigma3"][0][5]
        v4sigma4_aabb_fd = (v3sigma3_abb_plus - v3sigma3_abb_minus) / (delta_h *
                                                                       2)
        self.check_term(tol, "v4sigma4_aabb",
                        2 * v4sigma4_aabb * nabla_a + v4sigma4_acbb * nabla_b,
                        v4sigma4_aabb_fd)

        v4sigma4_accc = lxc["v4sigma4"][0][6]
        v4sigma4_cccc = lxc["v4sigma4"][0][10]
        v3sigma3_ccc_plus = kxc_plus["v3sigma3"][0][6]
        v3sigma3_ccc_minus = kxc_minus["v3sigma3"][0][6]
        v4sigma4_accc_fd = (v3sigma3_ccc_plus - v3sigma3_ccc_minus) / (delta_h *
                                                                       2)
        self.check_term(tol, "v4sigma4_accc",
                        2 * v4sigma4_accc * nabla_a + v4sigma4_cccc * nabla_b,
                        v4sigma4_accc_fd)

        v4sigma4_accb = lxc["v4sigma4"][0][7]
        v4sigma4_cccb = lxc["v4sigma4"][0][11]
        v3sigma3_ccb_plus = kxc_plus["v3sigma3"][0][7]
        v3sigma3_ccb_minus = kxc_minus["v3sigma3"][0][7]
        v4sigma4_accb_fd = (v3sigma3_ccb_plus - v3sigma3_ccb_minus) / (delta_h *
                                                                       2)
        self.check_term(tol, "v4sigma4_accb",
                        2 * v4sigma4_accb * nabla_a + v4sigma4_cccb * nabla_b,
                        v4sigma4_accb_fd)

        v4sigma4_acbb = lxc["v4sigma4"][0][8]
        v4sigma4_ccbb = lxc["v4sigma4"][0][12]
        v3sigma3_cbb_plus = kxc_plus["v3sigma3"][0][8]
        v3sigma3_cbb_minus = kxc_minus["v3sigma3"][0][8]
        v4sigma4_acbb_fd = (v3sigma3_cbb_plus - v3sigma3_cbb_minus) / (delta_h *
                                                                       2)
        self.check_term(tol, "v4sigma4_acbb",
                        2 * v4sigma4_acbb * nabla_a + v4sigma4_ccbb * nabla_b,
                        v4sigma4_acbb_fd)

        v4sigma4_abbb = lxc["v4sigma4"][0][9]
        v4sigma4_cbbb = lxc["v4sigma4"][0][13]
        v3sigma3_bbb_plus = kxc_plus["v3sigma3"][0][9]
        v3sigma3_bbb_minus = kxc_minus["v3sigma3"][0][9]
        v4sigma4_abbb_fd = (v3sigma3_bbb_plus - v3sigma3_bbb_minus) / (delta_h *
                                                                       2)
        self.check_term(tol, "v4sigma4_abbb",
                        2 * v4sigma4_abbb * nabla_a + v4sigma4_cbbb * nabla_b,
                        v4sigma4_abbb_fd)

    def test_sigma_a_fd_blyp(self):

        self.run_sigma_a_fd('blyp', 1.0e-7)

    def test_sigma_a_fd_b3lyp(self):

        self.run_sigma_a_fd('b3lyp', 1.0e-7)

    def test_sigma_a_fd_pbe0(self):

        self.run_sigma_a_fd('pbe0', 1.0e-7)

    def test_sigma_a_fd_bp86(self):

        self.run_sigma_a_fd('bp86', 1.0e-7)
