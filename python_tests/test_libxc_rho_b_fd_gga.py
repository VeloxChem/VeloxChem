import numpy as np

from veloxchem.veloxchemlib import XCNewIntegrator


class TestLibxcRhoB:

    def check_term(self, tol, label, val, val_fd):

        if abs(val_fd) < tol:
            assert abs(val - val_fd) < tol
        else:
            assert abs(val - val_fd) / abs(val_fd) < tol

    def run_rho_b_fd(self, xcfun_label, tol):

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
        inp_plus['rho'][1] = rho_b + delta_h

        kxc_plus = xc_drv.compute_kxc_for_gga(xcfun_label, inp_plus['rho'],
                                              inp_plus['sigma'])
        fxc_plus = xc_drv.compute_fxc_for_gga(xcfun_label, inp_plus['rho'],
                                              inp_plus['sigma'])
        vxc_plus = xc_drv.compute_exc_vxc_for_gga(xcfun_label, inp_plus['rho'],
                                                  inp_plus['sigma'])
        exc_plus = {'zk': vxc_plus['exc']}

        inp_minus = dict(inp)
        inp_minus['rho'][1] = rho_b - delta_h

        kxc_minus = xc_drv.compute_kxc_for_gga(xcfun_label, inp_minus['rho'],
                                               inp_minus['sigma'])
        fxc_minus = xc_drv.compute_fxc_for_gga(xcfun_label, inp_minus['rho'],
                                               inp_minus['sigma'])
        vxc_minus = xc_drv.compute_exc_vxc_for_gga(xcfun_label,
                                                   inp_minus['rho'],
                                                   inp_minus['sigma'])
        exc_minus = {'zk': vxc_minus['exc']}

        vrho_b = vxc['vrho'][0][1]
        e_a_plus = exc_plus['zk'][0][0] * (rho_a + rho_b + delta_h)
        e_a_minus = exc_minus['zk'][0][0] * (rho_a + rho_b - delta_h)
        vrho_b_fd = (e_a_plus - e_a_minus) / (delta_h * 2)
        self.check_term(tol, 'vrho_b', vrho_b, vrho_b_fd)

        v2rho2_bb = fxc["v2rho2"][0][2]
        vrho_b_plus = vxc_plus["vrho"][0][1]
        vrho_b_minus = vxc_minus["vrho"][0][1]
        v2rho2_bb_fd = (vrho_b_plus - vrho_b_minus) / (delta_h * 2)
        self.check_term(tol, "v2rho2_bb", v2rho2_bb, v2rho2_bb_fd)

        v2rhosigma_ba = fxc["v2rhosigma"][0][3]
        vsigma_a_plus = vxc_plus["vsigma"][0][0]
        vsigma_a_minus = vxc_minus["vsigma"][0][0]
        v2rhosigma_ba_fd = (vsigma_a_plus - vsigma_a_minus) / (delta_h * 2)
        self.check_term(tol, "v2rhosigma_ba", v2rhosigma_ba, v2rhosigma_ba_fd)

        v2rhosigma_bc = fxc["v2rhosigma"][0][4]
        vsigma_c_plus = vxc_plus["vsigma"][0][1]
        vsigma_c_minus = vxc_minus["vsigma"][0][1]
        v2rhosigma_bc_fd = (vsigma_c_plus - vsigma_c_minus) / (delta_h * 2)
        self.check_term(tol, "v2rhosigma_bc", v2rhosigma_bc, v2rhosigma_bc_fd)

        v2rhosigma_bb = fxc["v2rhosigma"][0][5]
        vsigma_b_plus = vxc_plus["vsigma"][0][2]
        vsigma_b_minus = vxc_minus["vsigma"][0][2]
        v2rhosigma_bb_fd = (vsigma_b_plus - vsigma_b_minus) / (delta_h * 2)
        self.check_term(tol, "v2rhosigma_bb", v2rhosigma_bb, v2rhosigma_bb_fd)

        v3rho3_bbb = kxc["v3rho3"][0][3]
        v2rho2_bb_plus = fxc_plus["v2rho2"][0][2]
        v2rho2_bb_minus = fxc_minus["v2rho2"][0][2]
        v3rho3_bbb_fd = (v2rho2_bb_plus - v2rho2_bb_minus) / (delta_h * 2)
        self.check_term(tol, "v3rho3_bbb", v3rho3_bbb, v3rho3_bbb_fd)

        v3rho2sigma_bba = kxc["v3rho2sigma"][0][6]
        v2rhosigma_ba_plus = fxc_plus["v2rhosigma"][0][3]
        v2rhosigma_ba_minus = fxc_minus["v2rhosigma"][0][3]
        v3rho2sigma_bba_fd = (v2rhosigma_ba_plus -
                              v2rhosigma_ba_minus) / (delta_h * 2)
        self.check_term(tol, "v3rho2sigma_bba", v3rho2sigma_bba,
                        v3rho2sigma_bba_fd)

        v3rho2sigma_bbc = kxc["v3rho2sigma"][0][7]
        v2rhosigma_bc_plus = fxc_plus["v2rhosigma"][0][4]
        v2rhosigma_bc_minus = fxc_minus["v2rhosigma"][0][4]
        v3rho2sigma_bbc_fd = (v2rhosigma_bc_plus -
                              v2rhosigma_bc_minus) / (delta_h * 2)
        self.check_term(tol, "v3rho2sigma_bbc", v3rho2sigma_bbc,
                        v3rho2sigma_bbc_fd)

        v3rho2sigma_bbb = kxc["v3rho2sigma"][0][8]
        v2rhosigma_bb_plus = fxc_plus["v2rhosigma"][0][5]
        v2rhosigma_bb_minus = fxc_minus["v2rhosigma"][0][5]
        v3rho2sigma_bbb_fd = (v2rhosigma_bb_plus -
                              v2rhosigma_bb_minus) / (delta_h * 2)
        self.check_term(tol, "v3rho2sigma_bbb", v3rho2sigma_bbb,
                        v3rho2sigma_bbb_fd)

        v3rhosigma2_baa = kxc["v3rhosigma2"][0][6]
        v2sigma2_aa_plus = fxc_plus["v2sigma2"][0][0]
        v2sigma2_aa_minus = fxc_minus["v2sigma2"][0][0]
        v3rhosigma2_baa_fd = (v2sigma2_aa_plus - v2sigma2_aa_minus) / (delta_h *
                                                                       2)
        self.check_term(tol, "v3rhosigma2_baa", v3rhosigma2_baa,
                        v3rhosigma2_baa_fd)

        v3rhosigma2_bac = kxc["v3rhosigma2"][0][7]
        v2sigma2_ac_plus = fxc_plus["v2sigma2"][0][1]
        v2sigma2_ac_minus = fxc_minus["v2sigma2"][0][1]
        v3rhosigma2_bac_fd = (v2sigma2_ac_plus - v2sigma2_ac_minus) / (delta_h *
                                                                       2)
        self.check_term(tol, "v3rhosigma2_bac", v3rhosigma2_bac,
                        v3rhosigma2_bac_fd)

        v3rhosigma2_bab = kxc["v3rhosigma2"][0][8]
        v2sigma2_ab_plus = fxc_plus["v2sigma2"][0][2]
        v2sigma2_ab_minus = fxc_minus["v2sigma2"][0][2]
        v3rhosigma2_bab_fd = (v2sigma2_ab_plus - v2sigma2_ab_minus) / (delta_h *
                                                                       2)
        self.check_term(tol, "v3rhosigma2_bab", v3rhosigma2_bab,
                        v3rhosigma2_bab_fd)

        v3rhosigma2_bcc = kxc["v3rhosigma2"][0][9]
        v2sigma2_cc_plus = fxc_plus["v2sigma2"][0][3]
        v2sigma2_cc_minus = fxc_minus["v2sigma2"][0][3]
        v3rhosigma2_bcc_fd = (v2sigma2_cc_plus - v2sigma2_cc_minus) / (delta_h *
                                                                       2)
        self.check_term(tol, "v3rhosigma2_bcc", v3rhosigma2_bcc,
                        v3rhosigma2_bcc_fd)

        v3rhosigma2_bcb = kxc["v3rhosigma2"][0][10]
        v2sigma2_cb_plus = fxc_plus["v2sigma2"][0][4]
        v2sigma2_cb_minus = fxc_minus["v2sigma2"][0][4]
        v3rhosigma2_bcb_fd = (v2sigma2_cb_plus - v2sigma2_cb_minus) / (delta_h *
                                                                       2)
        self.check_term(tol, "v3rhosigma2_bcb", v3rhosigma2_bcb,
                        v3rhosigma2_bcb_fd)

        v3rhosigma2_bbb = kxc["v3rhosigma2"][0][11]
        v2sigma2_bb_plus = fxc_plus["v2sigma2"][0][5]
        v2sigma2_bb_minus = fxc_minus["v2sigma2"][0][5]
        v3rhosigma2_bbb_fd = (v2sigma2_bb_plus - v2sigma2_bb_minus) / (delta_h *
                                                                       2)
        self.check_term(tol, "v3rhosigma2_bbb", v3rhosigma2_bbb,
                        v3rhosigma2_bbb_fd)

        v4rho4_bbbb = lxc["v4rho4"][0][4]
        v3rho3_bbb_plus = kxc_plus["v3rho3"][0][3]
        v3rho3_bbb_minus = kxc_minus["v3rho3"][0][3]
        v4rho4_bbbb_fd = (v3rho3_bbb_plus - v3rho3_bbb_minus) / (delta_h * 2)
        self.check_term(tol, "v4rho4_bbbb", v4rho4_bbbb, v4rho4_bbbb_fd)

        v4rho3sigma_bbba = lxc["v4rho3sigma"][0][9]
        v3rho2sigma_bba_plus = kxc_plus["v3rho2sigma"][0][6]
        v3rho2sigma_bba_minus = kxc_minus["v3rho2sigma"][0][6]
        v4rho3sigma_bbba_fd = (v3rho2sigma_bba_plus -
                               v3rho2sigma_bba_minus) / (delta_h * 2)
        self.check_term(tol, "v4rho3sigma_bbba", v4rho3sigma_bbba,
                        v4rho3sigma_bbba_fd)

        v4rho3sigma_bbbc = lxc["v4rho3sigma"][0][10]
        v3rho2sigma_bbc_plus = kxc_plus["v3rho2sigma"][0][7]
        v3rho2sigma_bbc_minus = kxc_minus["v3rho2sigma"][0][7]
        v4rho3sigma_bbbc_fd = (v3rho2sigma_bbc_plus -
                               v3rho2sigma_bbc_minus) / (delta_h * 2)
        self.check_term(tol, "v4rho3sigma_bbbc", v4rho3sigma_bbbc,
                        v4rho3sigma_bbbc_fd)

        v4rho3sigma_bbbb = lxc["v4rho3sigma"][0][11]
        v3rho2sigma_bbb_plus = kxc_plus["v3rho2sigma"][0][8]
        v3rho2sigma_bbb_minus = kxc_minus["v3rho2sigma"][0][8]
        v4rho3sigma_bbbb_fd = (v3rho2sigma_bbb_plus -
                               v3rho2sigma_bbb_minus) / (delta_h * 2)
        self.check_term(tol, "v4rho3sigma_bbbb", v4rho3sigma_bbbb,
                        v4rho3sigma_bbbb_fd)

        v4rho2sigma2_bbaa = lxc["v4rho2sigma2"][0][12]
        v3rhosigma2_baa_plus = kxc_plus["v3rhosigma2"][0][6]
        v3rhosigma2_baa_minus = kxc_minus["v3rhosigma2"][0][6]
        v4rho2sigma2_bbaa_fd = (v3rhosigma2_baa_plus -
                                v3rhosigma2_baa_minus) / (delta_h * 2)
        self.check_term(tol, "v4rho2sigma2_bbaa", v4rho2sigma2_bbaa,
                        v4rho2sigma2_bbaa_fd)

        v4rho2sigma2_bbac = lxc["v4rho2sigma2"][0][13]
        v3rhosigma2_bac_plus = kxc_plus["v3rhosigma2"][0][7]
        v3rhosigma2_bac_minus = kxc_minus["v3rhosigma2"][0][7]
        v4rho2sigma2_bbac_fd = (v3rhosigma2_bac_plus -
                                v3rhosigma2_bac_minus) / (delta_h * 2)
        self.check_term(tol, "v4rho2sigma2_bbac", v4rho2sigma2_bbac,
                        v4rho2sigma2_bbac_fd)

        v4rho2sigma2_bbab = lxc["v4rho2sigma2"][0][14]
        v3rhosigma2_bab_plus = kxc_plus["v3rhosigma2"][0][8]
        v3rhosigma2_bab_minus = kxc_minus["v3rhosigma2"][0][8]
        v4rho2sigma2_bbab_fd = (v3rhosigma2_bab_plus -
                                v3rhosigma2_bab_minus) / (delta_h * 2)
        self.check_term(tol, "v4rho2sigma2_bbab", v4rho2sigma2_bbab,
                        v4rho2sigma2_bbab_fd)

        v4rho2sigma2_bbcc = lxc["v4rho2sigma2"][0][15]
        v3rhosigma2_bcc_plus = kxc_plus["v3rhosigma2"][0][9]
        v3rhosigma2_bcc_minus = kxc_minus["v3rhosigma2"][0][9]
        v4rho2sigma2_bbcc_fd = (v3rhosigma2_bcc_plus -
                                v3rhosigma2_bcc_minus) / (delta_h * 2)
        self.check_term(tol, "v4rho2sigma2_bbcc", v4rho2sigma2_bbcc,
                        v4rho2sigma2_bbcc_fd)

        v4rho2sigma2_bbcb = lxc["v4rho2sigma2"][0][16]
        v3rhosigma2_bcb_plus = kxc_plus["v3rhosigma2"][0][10]
        v3rhosigma2_bcb_minus = kxc_minus["v3rhosigma2"][0][10]
        v4rho2sigma2_bbcb_fd = (v3rhosigma2_bcb_plus -
                                v3rhosigma2_bcb_minus) / (delta_h * 2)
        self.check_term(tol, "v4rho2sigma2_bbcb", v4rho2sigma2_bbcb,
                        v4rho2sigma2_bbcb_fd)

        v4rho2sigma2_bbbb = lxc["v4rho2sigma2"][0][17]
        v3rhosigma2_bbb_plus = kxc_plus["v3rhosigma2"][0][11]
        v3rhosigma2_bbb_minus = kxc_minus["v3rhosigma2"][0][11]
        v4rho2sigma2_bbbb_fd = (v3rhosigma2_bbb_plus -
                                v3rhosigma2_bbb_minus) / (delta_h * 2)
        self.check_term(tol, "v4rho2sigma2_bbbb", v4rho2sigma2_bbbb,
                        v4rho2sigma2_bbbb_fd)

        v4rhosigma3_baaa = lxc["v4rhosigma3"][0][10]
        v3sigma3_aaa_plus = kxc_plus["v3sigma3"][0][0]
        v3sigma3_aaa_minus = kxc_minus["v3sigma3"][0][0]
        v4rhosigma3_baaa_fd = (v3sigma3_aaa_plus -
                               v3sigma3_aaa_minus) / (delta_h * 2)
        self.check_term(tol, "v4rhosigma3_baaa", v4rhosigma3_baaa,
                        v4rhosigma3_baaa_fd)

        v4rhosigma3_baac = lxc["v4rhosigma3"][0][11]
        v3sigma3_aac_plus = kxc_plus["v3sigma3"][0][1]
        v3sigma3_aac_minus = kxc_minus["v3sigma3"][0][1]
        v4rhosigma3_baac_fd = (v3sigma3_aac_plus -
                               v3sigma3_aac_minus) / (delta_h * 2)
        self.check_term(tol, "v4rhosigma3_baac", v4rhosigma3_baac,
                        v4rhosigma3_baac_fd)

        v4rhosigma3_baab = lxc["v4rhosigma3"][0][12]
        v3sigma3_aab_plus = kxc_plus["v3sigma3"][0][2]
        v3sigma3_aab_minus = kxc_minus["v3sigma3"][0][2]
        v4rhosigma3_baab_fd = (v3sigma3_aab_plus -
                               v3sigma3_aab_minus) / (delta_h * 2)
        self.check_term(tol, "v4rhosigma3_baab", v4rhosigma3_baab,
                        v4rhosigma3_baab_fd)

        v4rhosigma3_bacc = lxc["v4rhosigma3"][0][13]
        v3sigma3_acc_plus = kxc_plus["v3sigma3"][0][3]
        v3sigma3_acc_minus = kxc_minus["v3sigma3"][0][3]
        v4rhosigma3_bacc_fd = (v3sigma3_acc_plus -
                               v3sigma3_acc_minus) / (delta_h * 2)
        self.check_term(tol, "v4rhosigma3_bacc", v4rhosigma3_bacc,
                        v4rhosigma3_bacc_fd)

        v4rhosigma3_bacb = lxc["v4rhosigma3"][0][14]
        v3sigma3_acb_plus = kxc_plus["v3sigma3"][0][4]
        v3sigma3_acb_minus = kxc_minus["v3sigma3"][0][4]
        v4rhosigma3_bacb_fd = (v3sigma3_acb_plus -
                               v3sigma3_acb_minus) / (delta_h * 2)
        self.check_term(tol, "v4rhosigma3_bacb", v4rhosigma3_bacb,
                        v4rhosigma3_bacb_fd)

        v4rhosigma3_babb = lxc["v4rhosigma3"][0][15]
        v3sigma3_abb_plus = kxc_plus["v3sigma3"][0][5]
        v3sigma3_abb_minus = kxc_minus["v3sigma3"][0][5]
        v4rhosigma3_babb_fd = (v3sigma3_abb_plus -
                               v3sigma3_abb_minus) / (delta_h * 2)
        self.check_term(tol, "v4rhosigma3_babb", v4rhosigma3_babb,
                        v4rhosigma3_babb_fd)

        v4rhosigma3_bccc = lxc["v4rhosigma3"][0][16]
        v3sigma3_ccc_plus = kxc_plus["v3sigma3"][0][6]
        v3sigma3_ccc_minus = kxc_minus["v3sigma3"][0][6]
        v4rhosigma3_bccc_fd = (v3sigma3_ccc_plus -
                               v3sigma3_ccc_minus) / (delta_h * 2)
        self.check_term(tol, "v4rhosigma3_bccc", v4rhosigma3_bccc,
                        v4rhosigma3_bccc_fd)

        v4rhosigma3_bccb = lxc["v4rhosigma3"][0][17]
        v3sigma3_ccb_plus = kxc_plus["v3sigma3"][0][7]
        v3sigma3_ccb_minus = kxc_minus["v3sigma3"][0][7]
        v4rhosigma3_bccb_fd = (v3sigma3_ccb_plus -
                               v3sigma3_ccb_minus) / (delta_h * 2)
        self.check_term(tol, "v4rhosigma3_bccb", v4rhosigma3_bccb,
                        v4rhosigma3_bccb_fd)

        v4rhosigma3_bcbb = lxc["v4rhosigma3"][0][18]
        v3sigma3_cbb_plus = kxc_plus["v3sigma3"][0][8]
        v3sigma3_cbb_minus = kxc_minus["v3sigma3"][0][8]
        v4rhosigma3_bcbb_fd = (v3sigma3_cbb_plus -
                               v3sigma3_cbb_minus) / (delta_h * 2)
        self.check_term(tol, "v4rhosigma3_bcbb", v4rhosigma3_bcbb,
                        v4rhosigma3_bcbb_fd)

        v4rhosigma3_bbbb = lxc["v4rhosigma3"][0][19]
        v3sigma3_bbb_plus = kxc_plus["v3sigma3"][0][9]
        v3sigma3_bbb_minus = kxc_minus["v3sigma3"][0][9]
        v4rhosigma3_bbbb_fd = (v3sigma3_bbb_plus -
                               v3sigma3_bbb_minus) / (delta_h * 2)
        self.check_term(tol, "v4rhosigma3_bbbb", v4rhosigma3_bbbb,
                        v4rhosigma3_bbbb_fd)

    def test_rho_b_fd_blyp(self):

        self.run_rho_b_fd('blyp', 1.0e-7)

    def test_rho_b_fd_b3lyp(self):

        self.run_rho_b_fd('b3lyp', 1.0e-7)

    def test_rho_b_fd_pbe0(self):

        self.run_rho_b_fd('pbe0', 1.0e-7)

    def test_rho_b_fd_bp86(self):

        self.run_rho_b_fd('bp86', 1.0e-7)
