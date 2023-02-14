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
        lapl_a, lapl_b = 0.1, 0.11
        tau_a, tau_b = 0.12, 0.125

        inp = {
            'rho': np.array([rho_a, rho_b]),
            'sigma': np.array([sigma_aa, sigma_ab, sigma_bb]),
            'lapl': np.array([lapl_a, lapl_b]),
            'tau': np.array([tau_a, tau_b]),
        }

        lxc = xc_drv.compute_lxc_for_mgga(xcfun_label, inp['rho'], inp['sigma'],
                                          inp['lapl'], inp['tau'])
        kxc = xc_drv.compute_kxc_for_mgga(xcfun_label, inp['rho'], inp['sigma'],
                                          inp['lapl'], inp['tau'])
        fxc = xc_drv.compute_fxc_for_mgga(xcfun_label, inp['rho'], inp['sigma'],
                                          inp['lapl'], inp['tau'])
        vxc = xc_drv.compute_exc_vxc_for_mgga(xcfun_label, inp['rho'],
                                              inp['sigma'], inp['lapl'],
                                              inp['tau'])

        inp_plus = dict(inp)
        nabla_a_plus = nabla_a + delta_h
        inp_plus['sigma'][0] = nabla_a_plus * nabla_a_plus
        inp_plus['sigma'][1] = nabla_a_plus * nabla_b

        kxc_plus = xc_drv.compute_kxc_for_mgga(xcfun_label, inp_plus['rho'],
                                               inp_plus['sigma'],
                                               inp_plus['lapl'],
                                               inp_plus['tau'])
        fxc_plus = xc_drv.compute_fxc_for_mgga(xcfun_label, inp_plus['rho'],
                                               inp_plus['sigma'],
                                               inp_plus['lapl'],
                                               inp_plus['tau'])
        vxc_plus = xc_drv.compute_exc_vxc_for_mgga(xcfun_label, inp_plus['rho'],
                                                   inp_plus['sigma'],
                                                   inp_plus['lapl'],
                                                   inp_plus['tau'])
        exc_plus = {'zk': vxc_plus['exc']}

        inp_minus = dict(inp)
        nabla_a_minus = nabla_a - delta_h
        inp_minus['sigma'][0] = nabla_a_minus * nabla_a_minus
        inp_minus['sigma'][1] = nabla_a_minus * nabla_b

        kxc_minus = xc_drv.compute_kxc_for_mgga(xcfun_label, inp_minus['rho'],
                                                inp_minus['sigma'],
                                                inp_minus['lapl'],
                                                inp_minus['tau'])
        fxc_minus = xc_drv.compute_fxc_for_mgga(xcfun_label, inp_minus['rho'],
                                                inp_minus['sigma'],
                                                inp_minus['lapl'],
                                                inp_minus['tau'])
        vxc_minus = xc_drv.compute_exc_vxc_for_mgga(xcfun_label,
                                                    inp_minus['rho'],
                                                    inp_minus['sigma'],
                                                    inp_minus['lapl'],
                                                    inp_minus['tau'])
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

        v2sigmatau_aa = fxc["v2sigmatau"][0][0]
        v2sigmatau_ca = fxc["v2sigmatau"][0][2]
        vtau_a_plus = vxc_plus["vtau"][0][0]
        vtau_a_minus = vxc_minus["vtau"][0][0]
        v2sigmatau_aa_fd = (vtau_a_plus - vtau_a_minus) / (delta_h * 2)
        self.check_term(tol, "v2sigmatau_aa",
                        2 * v2sigmatau_aa * nabla_a + v2sigmatau_ca * nabla_b,
                        v2sigmatau_aa_fd)

        v2sigmatau_ab = fxc["v2sigmatau"][0][1]
        v2sigmatau_cb = fxc["v2sigmatau"][0][3]
        vtau_b_plus = vxc_plus["vtau"][0][1]
        vtau_b_minus = vxc_minus["vtau"][0][1]
        v2sigmatau_ab_fd = (vtau_b_plus - vtau_b_minus) / (delta_h * 2)
        self.check_term(tol, "v2sigmatau_ab",
                        2 * v2sigmatau_ab * nabla_a + v2sigmatau_cb * nabla_b,
                        v2sigmatau_ab_fd)

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

        v3sigma2tau_aaa = kxc["v3sigma2tau"][0][0]
        v3sigma2tau_aca = kxc["v3sigma2tau"][0][2]
        v2sigmatau_aa_plus = fxc_plus["v2sigmatau"][0][0]
        v2sigmatau_aa_minus = fxc_minus["v2sigmatau"][0][0]
        v3sigma2tau_aaa_fd = (v2sigmatau_aa_plus -
                              v2sigmatau_aa_minus) / (delta_h * 2)
        self.check_term(
            tol, "v3sigma2tau_aaa",
            2 * v3sigma2tau_aaa * nabla_a + v3sigma2tau_aca * nabla_b,
            v3sigma2tau_aaa_fd)

        v3sigma2tau_aab = kxc["v3sigma2tau"][0][1]
        v3sigma2tau_acb = kxc["v3sigma2tau"][0][3]
        v2sigmatau_ab_plus = fxc_plus["v2sigmatau"][0][1]
        v2sigmatau_ab_minus = fxc_minus["v2sigmatau"][0][1]
        v3sigma2tau_aab_fd = (v2sigmatau_ab_plus -
                              v2sigmatau_ab_minus) / (delta_h * 2)
        self.check_term(
            tol, "v3sigma2tau_aab",
            2 * v3sigma2tau_aab * nabla_a + v3sigma2tau_acb * nabla_b,
            v3sigma2tau_aab_fd)

        v3sigma2tau_aca = kxc["v3sigma2tau"][0][2]
        v3sigma2tau_cca = kxc["v3sigma2tau"][0][6]
        v2sigmatau_ca_plus = fxc_plus["v2sigmatau"][0][2]
        v2sigmatau_ca_minus = fxc_minus["v2sigmatau"][0][2]
        v3sigma2tau_aca_fd = (v2sigmatau_ca_plus -
                              v2sigmatau_ca_minus) / (delta_h * 2)
        self.check_term(
            tol, "v3sigma2tau_aca",
            2 * v3sigma2tau_aca * nabla_a + v3sigma2tau_cca * nabla_b,
            v3sigma2tau_aca_fd)

        v3sigma2tau_acb = kxc["v3sigma2tau"][0][3]
        v3sigma2tau_ccb = kxc["v3sigma2tau"][0][7]
        v2sigmatau_cb_plus = fxc_plus["v2sigmatau"][0][3]
        v2sigmatau_cb_minus = fxc_minus["v2sigmatau"][0][3]
        v3sigma2tau_acb_fd = (v2sigmatau_cb_plus -
                              v2sigmatau_cb_minus) / (delta_h * 2)
        self.check_term(
            tol, "v3sigma2tau_acb",
            2 * v3sigma2tau_acb * nabla_a + v3sigma2tau_ccb * nabla_b,
            v3sigma2tau_acb_fd)

        v3sigma2tau_aba = kxc["v3sigma2tau"][0][4]
        v3sigma2tau_cba = kxc["v3sigma2tau"][0][8]
        v2sigmatau_ba_plus = fxc_plus["v2sigmatau"][0][4]
        v2sigmatau_ba_minus = fxc_minus["v2sigmatau"][0][4]
        v3sigma2tau_aba_fd = (v2sigmatau_ba_plus -
                              v2sigmatau_ba_minus) / (delta_h * 2)
        self.check_term(
            tol, "v3sigma2tau_aba",
            2 * v3sigma2tau_aba * nabla_a + v3sigma2tau_cba * nabla_b,
            v3sigma2tau_aba_fd)

        v3sigma2tau_abb = kxc["v3sigma2tau"][0][5]
        v3sigma2tau_cbb = kxc["v3sigma2tau"][0][9]
        v2sigmatau_bb_plus = fxc_plus["v2sigmatau"][0][5]
        v2sigmatau_bb_minus = fxc_minus["v2sigmatau"][0][5]
        v3sigma2tau_abb_fd = (v2sigmatau_bb_plus -
                              v2sigmatau_bb_minus) / (delta_h * 2)
        self.check_term(
            tol, "v3sigma2tau_abb",
            2 * v3sigma2tau_abb * nabla_a + v3sigma2tau_cbb * nabla_b,
            v3sigma2tau_abb_fd)

        v3sigmatau2_aaa = kxc["v3sigmatau2"][0][0]
        v3sigmatau2_caa = kxc["v3sigmatau2"][0][3]
        v2tau2_aa_plus = fxc_plus["v2tau2"][0][0]
        v2tau2_aa_minus = fxc_minus["v2tau2"][0][0]
        v3sigmatau2_aaa_fd = (v2tau2_aa_plus - v2tau2_aa_minus) / (delta_h * 2)
        self.check_term(
            tol, "v3sigmatau2_aaa",
            2 * v3sigmatau2_aaa * nabla_a + v3sigmatau2_caa * nabla_b,
            v3sigmatau2_aaa_fd)

        v3sigmatau2_aab = kxc["v3sigmatau2"][0][1]
        v3sigmatau2_cab = kxc["v3sigmatau2"][0][4]
        v2tau2_ab_plus = fxc_plus["v2tau2"][0][1]
        v2tau2_ab_minus = fxc_minus["v2tau2"][0][1]
        v3sigmatau2_aab_fd = (v2tau2_ab_plus - v2tau2_ab_minus) / (delta_h * 2)
        self.check_term(
            tol, "v3sigmatau2_aab",
            2 * v3sigmatau2_aab * nabla_a + v3sigmatau2_cab * nabla_b,
            v3sigmatau2_aab_fd)

        v3sigmatau2_abb = kxc["v3sigmatau2"][0][2]
        v3sigmatau2_cbb = kxc["v3sigmatau2"][0][5]
        v2tau2_bb_plus = fxc_plus["v2tau2"][0][2]
        v2tau2_bb_minus = fxc_minus["v2tau2"][0][2]
        v3sigmatau2_abb_fd = (v2tau2_bb_plus - v2tau2_bb_minus) / (delta_h * 2)
        self.check_term(
            tol, "v3sigmatau2_abb",
            2 * v3sigmatau2_abb * nabla_a + v3sigmatau2_cbb * nabla_b,
            v3sigmatau2_abb_fd)

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

        v4sigma3tau_aaaa = lxc["v4sigma3tau"][0][0]
        v4sigma3tau_aaca = lxc["v4sigma3tau"][0][2]
        v3sigma2tau_aaa_plus = kxc_plus["v3sigma2tau"][0][0]
        v3sigma2tau_aaa_minus = kxc_minus["v3sigma2tau"][0][0]
        v4sigma3tau_aaaa_fd = (v3sigma2tau_aaa_plus -
                               v3sigma2tau_aaa_minus) / (delta_h * 2)
        self.check_term(
            tol, "v4sigma3tau_aaaa",
            2 * v4sigma3tau_aaaa * nabla_a + v4sigma3tau_aaca * nabla_b,
            v4sigma3tau_aaaa_fd)

        v4sigma3tau_aaab = lxc["v4sigma3tau"][0][1]
        v4sigma3tau_aacb = lxc["v4sigma3tau"][0][3]
        v3sigma2tau_aab_plus = kxc_plus["v3sigma2tau"][0][1]
        v3sigma2tau_aab_minus = kxc_minus["v3sigma2tau"][0][1]
        v4sigma3tau_aaab_fd = (v3sigma2tau_aab_plus -
                               v3sigma2tau_aab_minus) / (delta_h * 2)
        self.check_term(
            tol, "v4sigma3tau_aaab",
            2 * v4sigma3tau_aaab * nabla_a + v4sigma3tau_aacb * nabla_b,
            v4sigma3tau_aaab_fd)

        v4sigma3tau_aaca = lxc["v4sigma3tau"][0][2]
        v4sigma3tau_acca = lxc["v4sigma3tau"][0][6]
        v3sigma2tau_aca_plus = kxc_plus["v3sigma2tau"][0][2]
        v3sigma2tau_aca_minus = kxc_minus["v3sigma2tau"][0][2]
        v4sigma3tau_aaca_fd = (v3sigma2tau_aca_plus -
                               v3sigma2tau_aca_minus) / (delta_h * 2)
        self.check_term(
            tol, "v4sigma3tau_aaca",
            2 * v4sigma3tau_aaca * nabla_a + v4sigma3tau_acca * nabla_b,
            v4sigma3tau_aaca_fd)

        v4sigma3tau_aacb = lxc["v4sigma3tau"][0][3]
        v4sigma3tau_accb = lxc["v4sigma3tau"][0][7]
        v3sigma2tau_acb_plus = kxc_plus["v3sigma2tau"][0][3]
        v3sigma2tau_acb_minus = kxc_minus["v3sigma2tau"][0][3]
        v4sigma3tau_aacb_fd = (v3sigma2tau_acb_plus -
                               v3sigma2tau_acb_minus) / (delta_h * 2)
        self.check_term(
            tol, "v4sigma3tau_aacb",
            2 * v4sigma3tau_aacb * nabla_a + v4sigma3tau_accb * nabla_b,
            v4sigma3tau_aacb_fd)

        v4sigma3tau_aaba = lxc["v4sigma3tau"][0][4]
        v4sigma3tau_acba = lxc["v4sigma3tau"][0][8]
        v3sigma2tau_aba_plus = kxc_plus["v3sigma2tau"][0][4]
        v3sigma2tau_aba_minus = kxc_minus["v3sigma2tau"][0][4]
        v4sigma3tau_aaba_fd = (v3sigma2tau_aba_plus -
                               v3sigma2tau_aba_minus) / (delta_h * 2)
        self.check_term(
            tol, "v4sigma3tau_aaba",
            2 * v4sigma3tau_aaba * nabla_a + v4sigma3tau_acba * nabla_b,
            v4sigma3tau_aaba_fd)

        v4sigma3tau_aabb = lxc["v4sigma3tau"][0][5]
        v4sigma3tau_acbb = lxc["v4sigma3tau"][0][9]
        v3sigma2tau_abb_plus = kxc_plus["v3sigma2tau"][0][5]
        v3sigma2tau_abb_minus = kxc_minus["v3sigma2tau"][0][5]
        v4sigma3tau_aabb_fd = (v3sigma2tau_abb_plus -
                               v3sigma2tau_abb_minus) / (delta_h * 2)
        self.check_term(
            tol, "v4sigma3tau_aabb",
            2 * v4sigma3tau_aabb * nabla_a + v4sigma3tau_acbb * nabla_b,
            v4sigma3tau_aabb_fd)

        v4sigma3tau_acca = lxc["v4sigma3tau"][0][6]
        v4sigma3tau_ccca = lxc["v4sigma3tau"][0][12]
        v3sigma2tau_cca_plus = kxc_plus["v3sigma2tau"][0][6]
        v3sigma2tau_cca_minus = kxc_minus["v3sigma2tau"][0][6]
        v4sigma3tau_acca_fd = (v3sigma2tau_cca_plus -
                               v3sigma2tau_cca_minus) / (delta_h * 2)
        self.check_term(
            tol, "v4sigma3tau_acca",
            2 * v4sigma3tau_acca * nabla_a + v4sigma3tau_ccca * nabla_b,
            v4sigma3tau_acca_fd)

        v4sigma3tau_accb = lxc["v4sigma3tau"][0][7]
        v4sigma3tau_cccb = lxc["v4sigma3tau"][0][13]
        v3sigma2tau_ccb_plus = kxc_plus["v3sigma2tau"][0][7]
        v3sigma2tau_ccb_minus = kxc_minus["v3sigma2tau"][0][7]
        v4sigma3tau_accb_fd = (v3sigma2tau_ccb_plus -
                               v3sigma2tau_ccb_minus) / (delta_h * 2)
        self.check_term(
            tol, "v4sigma3tau_accb",
            2 * v4sigma3tau_accb * nabla_a + v4sigma3tau_cccb * nabla_b,
            v4sigma3tau_accb_fd)

        v4sigma3tau_acba = lxc["v4sigma3tau"][0][8]
        v4sigma3tau_ccba = lxc["v4sigma3tau"][0][14]
        v3sigma2tau_cba_plus = kxc_plus["v3sigma2tau"][0][8]
        v3sigma2tau_cba_minus = kxc_minus["v3sigma2tau"][0][8]
        v4sigma3tau_acba_fd = (v3sigma2tau_cba_plus -
                               v3sigma2tau_cba_minus) / (delta_h * 2)
        self.check_term(
            tol, "v4sigma3tau_acba",
            2 * v4sigma3tau_acba * nabla_a + v4sigma3tau_ccba * nabla_b,
            v4sigma3tau_acba_fd)

        v4sigma3tau_acbb = lxc["v4sigma3tau"][0][9]
        v4sigma3tau_ccbb = lxc["v4sigma3tau"][0][15]
        v3sigma2tau_cbb_plus = kxc_plus["v3sigma2tau"][0][9]
        v3sigma2tau_cbb_minus = kxc_minus["v3sigma2tau"][0][9]
        v4sigma3tau_acbb_fd = (v3sigma2tau_cbb_plus -
                               v3sigma2tau_cbb_minus) / (delta_h * 2)
        self.check_term(
            tol, "v4sigma3tau_acbb",
            2 * v4sigma3tau_acbb * nabla_a + v4sigma3tau_ccbb * nabla_b,
            v4sigma3tau_acbb_fd)

        v4sigma3tau_abba = lxc["v4sigma3tau"][0][10]
        v4sigma3tau_cbba = lxc["v4sigma3tau"][0][16]
        v3sigma2tau_bba_plus = kxc_plus["v3sigma2tau"][0][10]
        v3sigma2tau_bba_minus = kxc_minus["v3sigma2tau"][0][10]
        v4sigma3tau_abba_fd = (v3sigma2tau_bba_plus -
                               v3sigma2tau_bba_minus) / (delta_h * 2)
        self.check_term(
            tol, "v4sigma3tau_abba",
            2 * v4sigma3tau_abba * nabla_a + v4sigma3tau_cbba * nabla_b,
            v4sigma3tau_abba_fd)

        v4sigma3tau_abbb = lxc["v4sigma3tau"][0][11]
        v4sigma3tau_cbbb = lxc["v4sigma3tau"][0][17]
        v3sigma2tau_bbb_plus = kxc_plus["v3sigma2tau"][0][11]
        v3sigma2tau_bbb_minus = kxc_minus["v3sigma2tau"][0][11]
        v4sigma3tau_abbb_fd = (v3sigma2tau_bbb_plus -
                               v3sigma2tau_bbb_minus) / (delta_h * 2)
        self.check_term(
            tol, "v4sigma3tau_abbb",
            2 * v4sigma3tau_abbb * nabla_a + v4sigma3tau_cbbb * nabla_b,
            v4sigma3tau_abbb_fd)

        v4sigma2tau2_aaaa = lxc["v4sigma2tau2"][0][0]
        v4sigma2tau2_acaa = lxc["v4sigma2tau2"][0][3]
        v3sigmatau2_aaa_plus = kxc_plus["v3sigmatau2"][0][0]
        v3sigmatau2_aaa_minus = kxc_minus["v3sigmatau2"][0][0]
        v4sigma2tau2_aaaa_fd = (v3sigmatau2_aaa_plus -
                                v3sigmatau2_aaa_minus) / (delta_h * 2)
        self.check_term(
            tol, "v4sigma2tau2_aaaa",
            2 * v4sigma2tau2_aaaa * nabla_a + v4sigma2tau2_acaa * nabla_b,
            v4sigma2tau2_aaaa_fd)

        v4sigma2tau2_aaab = lxc["v4sigma2tau2"][0][1]
        v4sigma2tau2_acab = lxc["v4sigma2tau2"][0][4]
        v3sigmatau2_aab_plus = kxc_plus["v3sigmatau2"][0][1]
        v3sigmatau2_aab_minus = kxc_minus["v3sigmatau2"][0][1]
        v4sigma2tau2_aaab_fd = (v3sigmatau2_aab_plus -
                                v3sigmatau2_aab_minus) / (delta_h * 2)
        self.check_term(
            tol, "v4sigma2tau2_aaab",
            2 * v4sigma2tau2_aaab * nabla_a + v4sigma2tau2_acab * nabla_b,
            v4sigma2tau2_aaab_fd)

        v4sigma2tau2_aabb = lxc["v4sigma2tau2"][0][2]
        v4sigma2tau2_acbb = lxc["v4sigma2tau2"][0][5]
        v3sigmatau2_abb_plus = kxc_plus["v3sigmatau2"][0][2]
        v3sigmatau2_abb_minus = kxc_minus["v3sigmatau2"][0][2]
        v4sigma2tau2_aabb_fd = (v3sigmatau2_abb_plus -
                                v3sigmatau2_abb_minus) / (delta_h * 2)
        self.check_term(
            tol, "v4sigma2tau2_aabb",
            2 * v4sigma2tau2_aabb * nabla_a + v4sigma2tau2_acbb * nabla_b,
            v4sigma2tau2_aabb_fd)

        v4sigma2tau2_acaa = lxc["v4sigma2tau2"][0][3]
        v4sigma2tau2_ccaa = lxc["v4sigma2tau2"][0][9]
        v3sigmatau2_caa_plus = kxc_plus["v3sigmatau2"][0][3]
        v3sigmatau2_caa_minus = kxc_minus["v3sigmatau2"][0][3]
        v4sigma2tau2_acaa_fd = (v3sigmatau2_caa_plus -
                                v3sigmatau2_caa_minus) / (delta_h * 2)
        self.check_term(
            tol, "v4sigma2tau2_acaa",
            2 * v4sigma2tau2_acaa * nabla_a + v4sigma2tau2_ccaa * nabla_b,
            v4sigma2tau2_acaa_fd)

        v4sigma2tau2_acab = lxc["v4sigma2tau2"][0][4]
        v4sigma2tau2_ccab = lxc["v4sigma2tau2"][0][10]
        v3sigmatau2_cab_plus = kxc_plus["v3sigmatau2"][0][4]
        v3sigmatau2_cab_minus = kxc_minus["v3sigmatau2"][0][4]
        v4sigma2tau2_acab_fd = (v3sigmatau2_cab_plus -
                                v3sigmatau2_cab_minus) / (delta_h * 2)
        self.check_term(
            tol, "v4sigma2tau2_acab",
            2 * v4sigma2tau2_acab * nabla_a + v4sigma2tau2_ccab * nabla_b,
            v4sigma2tau2_acab_fd)

        v4sigma2tau2_acbb = lxc["v4sigma2tau2"][0][5]
        v4sigma2tau2_ccbb = lxc["v4sigma2tau2"][0][11]
        v3sigmatau2_cbb_plus = kxc_plus["v3sigmatau2"][0][5]
        v3sigmatau2_cbb_minus = kxc_minus["v3sigmatau2"][0][5]
        v4sigma2tau2_acbb_fd = (v3sigmatau2_cbb_plus -
                                v3sigmatau2_cbb_minus) / (delta_h * 2)
        self.check_term(
            tol, "v4sigma2tau2_acbb",
            2 * v4sigma2tau2_acbb * nabla_a + v4sigma2tau2_ccbb * nabla_b,
            v4sigma2tau2_acbb_fd)

        v4sigma2tau2_abaa = lxc["v4sigma2tau2"][0][6]
        v4sigma2tau2_cbaa = lxc["v4sigma2tau2"][0][12]
        v3sigmatau2_baa_plus = kxc_plus["v3sigmatau2"][0][6]
        v3sigmatau2_baa_minus = kxc_minus["v3sigmatau2"][0][6]
        v4sigma2tau2_abaa_fd = (v3sigmatau2_baa_plus -
                                v3sigmatau2_baa_minus) / (delta_h * 2)
        self.check_term(
            tol, "v4sigma2tau2_abaa",
            2 * v4sigma2tau2_abaa * nabla_a + v4sigma2tau2_cbaa * nabla_b,
            v4sigma2tau2_abaa_fd)

        v4sigma2tau2_abab = lxc["v4sigma2tau2"][0][7]
        v4sigma2tau2_cbab = lxc["v4sigma2tau2"][0][13]
        v3sigmatau2_bab_plus = kxc_plus["v3sigmatau2"][0][7]
        v3sigmatau2_bab_minus = kxc_minus["v3sigmatau2"][0][7]
        v4sigma2tau2_abab_fd = (v3sigmatau2_bab_plus -
                                v3sigmatau2_bab_minus) / (delta_h * 2)
        self.check_term(
            tol, "v4sigma2tau2_abab",
            2 * v4sigma2tau2_abab * nabla_a + v4sigma2tau2_cbab * nabla_b,
            v4sigma2tau2_abab_fd)

        v4sigma2tau2_abbb = lxc["v4sigma2tau2"][0][8]
        v4sigma2tau2_cbbb = lxc["v4sigma2tau2"][0][14]
        v3sigmatau2_bbb_plus = kxc_plus["v3sigmatau2"][0][8]
        v3sigmatau2_bbb_minus = kxc_minus["v3sigmatau2"][0][8]
        v4sigma2tau2_abbb_fd = (v3sigmatau2_bbb_plus -
                                v3sigmatau2_bbb_minus) / (delta_h * 2)
        self.check_term(
            tol, "v4sigma2tau2_abbb",
            2 * v4sigma2tau2_abbb * nabla_a + v4sigma2tau2_cbbb * nabla_b,
            v4sigma2tau2_abbb_fd)

        v4sigmatau3_aaaa = lxc["v4sigmatau3"][0][0]
        v4sigmatau3_caaa = lxc["v4sigmatau3"][0][4]
        v3tau3_aaa_plus = kxc_plus["v3tau3"][0][0]
        v3tau3_aaa_minus = kxc_minus["v3tau3"][0][0]
        v4sigmatau3_aaaa_fd = (v3tau3_aaa_plus - v3tau3_aaa_minus) / (delta_h *
                                                                      2)
        self.check_term(
            tol, "v4sigmatau3_aaaa",
            2 * v4sigmatau3_aaaa * nabla_a + v4sigmatau3_caaa * nabla_b,
            v4sigmatau3_aaaa_fd)

        v4sigmatau3_aaab = lxc["v4sigmatau3"][0][1]
        v4sigmatau3_caab = lxc["v4sigmatau3"][0][5]
        v3tau3_aab_plus = kxc_plus["v3tau3"][0][1]
        v3tau3_aab_minus = kxc_minus["v3tau3"][0][1]
        v4sigmatau3_aaab_fd = (v3tau3_aab_plus - v3tau3_aab_minus) / (delta_h *
                                                                      2)
        self.check_term(
            tol, "v4sigmatau3_aaab",
            2 * v4sigmatau3_aaab * nabla_a + v4sigmatau3_caab * nabla_b,
            v4sigmatau3_aaab_fd)

        v4sigmatau3_aabb = lxc["v4sigmatau3"][0][2]
        v4sigmatau3_cabb = lxc["v4sigmatau3"][0][6]
        v3tau3_abb_plus = kxc_plus["v3tau3"][0][2]
        v3tau3_abb_minus = kxc_minus["v3tau3"][0][2]
        v4sigmatau3_aabb_fd = (v3tau3_abb_plus - v3tau3_abb_minus) / (delta_h *
                                                                      2)
        self.check_term(
            tol, "v4sigmatau3_aabb",
            2 * v4sigmatau3_aabb * nabla_a + v4sigmatau3_cabb * nabla_b,
            v4sigmatau3_aabb_fd)

        v4sigmatau3_abbb = lxc["v4sigmatau3"][0][3]
        v4sigmatau3_cbbb = lxc["v4sigmatau3"][0][7]
        v3tau3_bbb_plus = kxc_plus["v3tau3"][0][3]
        v3tau3_bbb_minus = kxc_minus["v3tau3"][0][3]
        v4sigmatau3_abbb_fd = (v3tau3_bbb_plus - v3tau3_bbb_minus) / (delta_h *
                                                                      2)
        self.check_term(
            tol, "v4sigmatau3_abbb",
            2 * v4sigmatau3_abbb * nabla_a + v4sigmatau3_cbbb * nabla_b,
            v4sigmatau3_abbb_fd)

    def test_sigma_a_fd_pkzb(self):

        self.run_sigma_a_fd('pkzb', 1.0e-7)

    def test_sigma_a_fd_tpssh(self):

        self.run_sigma_a_fd('tpssh', 1.0e-7)

    def test_sigma_a_fd_m06(self):

        self.run_sigma_a_fd('m06', 1.0e-6)

    def test_sigma_a_fd_scan(self):

        self.run_sigma_a_fd('scan', 1.0e-7)
