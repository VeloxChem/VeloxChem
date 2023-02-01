import numpy as np

from veloxchem.veloxchemlib import XCNewIntegrator


class TestLibxcRhoA:

    def check_term(self, tol, label, val, val_fd):

        if abs(val_fd) < tol:
            assert abs(val - val_fd) < tol
        else:
            assert abs(val - val_fd) / abs(val_fd) < tol

    def run_rho_a_fd(self, xcfun_label, tol):

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
        inp_plus['rho'][0] = rho_a + delta_h

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
        inp_minus['rho'][0] = rho_a - delta_h

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

        vrho_a = vxc['vrho'][0][0]
        e_a_plus = exc_plus['zk'][0][0] * (rho_a + rho_b + delta_h)
        e_a_minus = exc_minus['zk'][0][0] * (rho_a + rho_b - delta_h)
        vrho_a_fd = (e_a_plus - e_a_minus) / (delta_h * 2)
        self.check_term(tol, 'vrho_a', vrho_a, vrho_a_fd)

        v2rho2_aa = fxc["v2rho2"][0][0]
        vrho_a_plus = vxc_plus["vrho"][0][0]
        vrho_a_minus = vxc_minus["vrho"][0][0]
        v2rho2_aa_fd = (vrho_a_plus - vrho_a_minus) / (delta_h * 2)
        self.check_term(tol, "v2rho2_aa", v2rho2_aa, v2rho2_aa_fd)

        v2rho2_ab = fxc["v2rho2"][0][1]
        vrho_b_plus = vxc_plus["vrho"][0][1]
        vrho_b_minus = vxc_minus["vrho"][0][1]
        v2rho2_ab_fd = (vrho_b_plus - vrho_b_minus) / (delta_h * 2)
        self.check_term(tol, "v2rho2_ab", v2rho2_ab, v2rho2_ab_fd)

        v2rhosigma_aa = fxc["v2rhosigma"][0][0]
        vsigma_a_plus = vxc_plus["vsigma"][0][0]
        vsigma_a_minus = vxc_minus["vsigma"][0][0]
        v2rhosigma_aa_fd = (vsigma_a_plus - vsigma_a_minus) / (delta_h * 2)
        self.check_term(tol, "v2rhosigma_aa", v2rhosigma_aa, v2rhosigma_aa_fd)

        v2rhosigma_ac = fxc["v2rhosigma"][0][1]
        vsigma_c_plus = vxc_plus["vsigma"][0][1]
        vsigma_c_minus = vxc_minus["vsigma"][0][1]
        v2rhosigma_ac_fd = (vsigma_c_plus - vsigma_c_minus) / (delta_h * 2)
        self.check_term(tol, "v2rhosigma_ac", v2rhosigma_ac, v2rhosigma_ac_fd)

        v2rhosigma_ab = fxc["v2rhosigma"][0][2]
        vsigma_b_plus = vxc_plus["vsigma"][0][2]
        vsigma_b_minus = vxc_minus["vsigma"][0][2]
        v2rhosigma_ab_fd = (vsigma_b_plus - vsigma_b_minus) / (delta_h * 2)
        self.check_term(tol, "v2rhosigma_ab", v2rhosigma_ab, v2rhosigma_ab_fd)

        v2rhotau_aa = fxc["v2rhotau"][0][0]
        vtau_a_plus = vxc_plus["vtau"][0][0]
        vtau_a_minus = vxc_minus["vtau"][0][0]
        v2rhotau_aa_fd = (vtau_a_plus - vtau_a_minus) / (delta_h * 2)
        self.check_term(tol, "v2rhotau_aa", v2rhotau_aa, v2rhotau_aa_fd)

        v2rhotau_ab = fxc["v2rhotau"][0][1]
        vtau_b_plus = vxc_plus["vtau"][0][1]
        vtau_b_minus = vxc_minus["vtau"][0][1]
        v2rhotau_ab_fd = (vtau_b_plus - vtau_b_minus) / (delta_h * 2)
        self.check_term(tol, "v2rhotau_ab", v2rhotau_ab, v2rhotau_ab_fd)

        v3rho3_aaa = kxc["v3rho3"][0][0]
        v2rho2_aa_plus = fxc_plus["v2rho2"][0][0]
        v2rho2_aa_minus = fxc_minus["v2rho2"][0][0]
        v3rho3_aaa_fd = (v2rho2_aa_plus - v2rho2_aa_minus) / (delta_h * 2)
        self.check_term(tol, "v3rho3_aaa", v3rho3_aaa, v3rho3_aaa_fd)

        v3rho3_aab = kxc["v3rho3"][0][1]
        v2rho2_ab_plus = fxc_plus["v2rho2"][0][1]
        v2rho2_ab_minus = fxc_minus["v2rho2"][0][1]
        v3rho3_aab_fd = (v2rho2_ab_plus - v2rho2_ab_minus) / (delta_h * 2)
        self.check_term(tol, "v3rho3_aab", v3rho3_aab, v3rho3_aab_fd)

        v3rho3_abb = kxc["v3rho3"][0][2]
        v2rho2_bb_plus = fxc_plus["v2rho2"][0][2]
        v2rho2_bb_minus = fxc_minus["v2rho2"][0][2]
        v3rho3_abb_fd = (v2rho2_bb_plus - v2rho2_bb_minus) / (delta_h * 2)
        self.check_term(tol, "v3rho3_abb", v3rho3_abb, v3rho3_abb_fd)

        v3rho2sigma_aaa = kxc["v3rho2sigma"][0][0]
        v2rhosigma_aa_plus = fxc_plus["v2rhosigma"][0][0]
        v2rhosigma_aa_minus = fxc_minus["v2rhosigma"][0][0]
        v3rho2sigma_aaa_fd = (v2rhosigma_aa_plus -
                              v2rhosigma_aa_minus) / (delta_h * 2)
        self.check_term(tol, "v3rho2sigma_aaa", v3rho2sigma_aaa,
                        v3rho2sigma_aaa_fd)

        v3rho2sigma_aac = kxc["v3rho2sigma"][0][1]
        v2rhosigma_ac_plus = fxc_plus["v2rhosigma"][0][1]
        v2rhosigma_ac_minus = fxc_minus["v2rhosigma"][0][1]
        v3rho2sigma_aac_fd = (v2rhosigma_ac_plus -
                              v2rhosigma_ac_minus) / (delta_h * 2)
        self.check_term(tol, "v3rho2sigma_aac", v3rho2sigma_aac,
                        v3rho2sigma_aac_fd)

        v3rho2sigma_aab = kxc["v3rho2sigma"][0][2]
        v2rhosigma_ab_plus = fxc_plus["v2rhosigma"][0][2]
        v2rhosigma_ab_minus = fxc_minus["v2rhosigma"][0][2]
        v3rho2sigma_aab_fd = (v2rhosigma_ab_plus -
                              v2rhosigma_ab_minus) / (delta_h * 2)
        self.check_term(tol, "v3rho2sigma_aab", v3rho2sigma_aab,
                        v3rho2sigma_aab_fd)

        v3rho2sigma_aba = kxc["v3rho2sigma"][0][3]
        v2rhosigma_ba_plus = fxc_plus["v2rhosigma"][0][3]
        v2rhosigma_ba_minus = fxc_minus["v2rhosigma"][0][3]
        v3rho2sigma_aba_fd = (v2rhosigma_ba_plus -
                              v2rhosigma_ba_minus) / (delta_h * 2)
        self.check_term(tol, "v3rho2sigma_aba", v3rho2sigma_aba,
                        v3rho2sigma_aba_fd)

        v3rho2sigma_abc = kxc["v3rho2sigma"][0][4]
        v2rhosigma_bc_plus = fxc_plus["v2rhosigma"][0][4]
        v2rhosigma_bc_minus = fxc_minus["v2rhosigma"][0][4]
        v3rho2sigma_abc_fd = (v2rhosigma_bc_plus -
                              v2rhosigma_bc_minus) / (delta_h * 2)
        self.check_term(tol, "v3rho2sigma_abc", v3rho2sigma_abc,
                        v3rho2sigma_abc_fd)

        v3rho2sigma_abb = kxc["v3rho2sigma"][0][5]
        v2rhosigma_bb_plus = fxc_plus["v2rhosigma"][0][5]
        v2rhosigma_bb_minus = fxc_minus["v2rhosigma"][0][5]
        v3rho2sigma_abb_fd = (v2rhosigma_bb_plus -
                              v2rhosigma_bb_minus) / (delta_h * 2)
        self.check_term(tol, "v3rho2sigma_abb", v3rho2sigma_abb,
                        v3rho2sigma_abb_fd)

        v3rho2tau_aaa = kxc["v3rho2tau"][0][0]
        v2rhotau_aa_plus = fxc_plus["v2rhotau"][0][0]
        v2rhotau_aa_minus = fxc_minus["v2rhotau"][0][0]
        v3rho2tau_aaa_fd = (v2rhotau_aa_plus - v2rhotau_aa_minus) / (delta_h *
                                                                     2)
        self.check_term(tol, "v3rho2tau_aaa", v3rho2tau_aaa, v3rho2tau_aaa_fd)

        v3rho2tau_aab = kxc["v3rho2tau"][0][1]
        v2rhotau_ab_plus = fxc_plus["v2rhotau"][0][1]
        v2rhotau_ab_minus = fxc_minus["v2rhotau"][0][1]
        v3rho2tau_aab_fd = (v2rhotau_ab_plus - v2rhotau_ab_minus) / (delta_h *
                                                                     2)
        self.check_term(tol, "v3rho2tau_aab", v3rho2tau_aab, v3rho2tau_aab_fd)

        v3rho2tau_aba = kxc["v3rho2tau"][0][2]
        v2rhotau_ba_plus = fxc_plus["v2rhotau"][0][2]
        v2rhotau_ba_minus = fxc_minus["v2rhotau"][0][2]
        v3rho2tau_aba_fd = (v2rhotau_ba_plus - v2rhotau_ba_minus) / (delta_h *
                                                                     2)
        self.check_term(tol, "v3rho2tau_aba", v3rho2tau_aba, v3rho2tau_aba_fd)

        v3rho2tau_abb = kxc["v3rho2tau"][0][3]
        v2rhotau_bb_plus = fxc_plus["v2rhotau"][0][3]
        v2rhotau_bb_minus = fxc_minus["v2rhotau"][0][3]
        v3rho2tau_abb_fd = (v2rhotau_bb_plus - v2rhotau_bb_minus) / (delta_h *
                                                                     2)
        self.check_term(tol, "v3rho2tau_abb", v3rho2tau_abb, v3rho2tau_abb_fd)

        v3rhosigma2_aaa = kxc["v3rhosigma2"][0][0]
        v2sigma2_aa_plus = fxc_plus["v2sigma2"][0][0]
        v2sigma2_aa_minus = fxc_minus["v2sigma2"][0][0]
        v3rhosigma2_aaa_fd = (v2sigma2_aa_plus - v2sigma2_aa_minus) / (delta_h *
                                                                       2)
        self.check_term(tol, "v3rhosigma2_aaa", v3rhosigma2_aaa,
                        v3rhosigma2_aaa_fd)

        v3rhosigma2_aac = kxc["v3rhosigma2"][0][1]
        v2sigma2_ac_plus = fxc_plus["v2sigma2"][0][1]
        v2sigma2_ac_minus = fxc_minus["v2sigma2"][0][1]
        v3rhosigma2_aac_fd = (v2sigma2_ac_plus - v2sigma2_ac_minus) / (delta_h *
                                                                       2)
        self.check_term(tol, "v3rhosigma2_aac", v3rhosigma2_aac,
                        v3rhosigma2_aac_fd)

        v3rhosigma2_aab = kxc["v3rhosigma2"][0][2]
        v2sigma2_ab_plus = fxc_plus["v2sigma2"][0][2]
        v2sigma2_ab_minus = fxc_minus["v2sigma2"][0][2]
        v3rhosigma2_aab_fd = (v2sigma2_ab_plus - v2sigma2_ab_minus) / (delta_h *
                                                                       2)
        self.check_term(tol, "v3rhosigma2_aab", v3rhosigma2_aab,
                        v3rhosigma2_aab_fd)

        v3rhosigma2_acc = kxc["v3rhosigma2"][0][3]
        v2sigma2_cc_plus = fxc_plus["v2sigma2"][0][3]
        v2sigma2_cc_minus = fxc_minus["v2sigma2"][0][3]
        v3rhosigma2_acc_fd = (v2sigma2_cc_plus - v2sigma2_cc_minus) / (delta_h *
                                                                       2)
        self.check_term(tol, "v3rhosigma2_acc", v3rhosigma2_acc,
                        v3rhosigma2_acc_fd)

        v3rhosigma2_acb = kxc["v3rhosigma2"][0][4]
        v2sigma2_cb_plus = fxc_plus["v2sigma2"][0][4]
        v2sigma2_cb_minus = fxc_minus["v2sigma2"][0][4]
        v3rhosigma2_acb_fd = (v2sigma2_cb_plus - v2sigma2_cb_minus) / (delta_h *
                                                                       2)
        self.check_term(tol, "v3rhosigma2_acb", v3rhosigma2_acb,
                        v3rhosigma2_acb_fd)

        v3rhosigma2_abb = kxc["v3rhosigma2"][0][5]
        v2sigma2_bb_plus = fxc_plus["v2sigma2"][0][5]
        v2sigma2_bb_minus = fxc_minus["v2sigma2"][0][5]
        v3rhosigma2_abb_fd = (v2sigma2_bb_plus - v2sigma2_bb_minus) / (delta_h *
                                                                       2)
        self.check_term(tol, "v3rhosigma2_abb", v3rhosigma2_abb,
                        v3rhosigma2_abb_fd)

        v3rhosigmatau_aaa = kxc["v3rhosigmatau"][0][0]
        v2sigmatau_aa_plus = fxc_plus["v2sigmatau"][0][0]
        v2sigmatau_aa_minus = fxc_minus["v2sigmatau"][0][0]
        v3rhosigmatau_aaa_fd = (v2sigmatau_aa_plus -
                                v2sigmatau_aa_minus) / (delta_h * 2)
        self.check_term(tol, "v3rhosigmatau_aaa", v3rhosigmatau_aaa,
                        v3rhosigmatau_aaa_fd)

        v3rhosigmatau_aab = kxc["v3rhosigmatau"][0][1]
        v2sigmatau_ab_plus = fxc_plus["v2sigmatau"][0][1]
        v2sigmatau_ab_minus = fxc_minus["v2sigmatau"][0][1]
        v3rhosigmatau_aab_fd = (v2sigmatau_ab_plus -
                                v2sigmatau_ab_minus) / (delta_h * 2)
        self.check_term(tol, "v3rhosigmatau_aab", v3rhosigmatau_aab,
                        v3rhosigmatau_aab_fd)

        v3rhosigmatau_aca = kxc["v3rhosigmatau"][0][2]
        v2sigmatau_ca_plus = fxc_plus["v2sigmatau"][0][2]
        v2sigmatau_ca_minus = fxc_minus["v2sigmatau"][0][2]
        v3rhosigmatau_aca_fd = (v2sigmatau_ca_plus -
                                v2sigmatau_ca_minus) / (delta_h * 2)
        self.check_term(tol, "v3rhosigmatau_aca", v3rhosigmatau_aca,
                        v3rhosigmatau_aca_fd)

        v3rhosigmatau_acb = kxc["v3rhosigmatau"][0][3]
        v2sigmatau_cb_plus = fxc_plus["v2sigmatau"][0][3]
        v2sigmatau_cb_minus = fxc_minus["v2sigmatau"][0][3]
        v3rhosigmatau_acb_fd = (v2sigmatau_cb_plus -
                                v2sigmatau_cb_minus) / (delta_h * 2)
        self.check_term(tol, "v3rhosigmatau_acb", v3rhosigmatau_acb,
                        v3rhosigmatau_acb_fd)

        v3rhosigmatau_aba = kxc["v3rhosigmatau"][0][4]
        v2sigmatau_ba_plus = fxc_plus["v2sigmatau"][0][4]
        v2sigmatau_ba_minus = fxc_minus["v2sigmatau"][0][4]
        v3rhosigmatau_aba_fd = (v2sigmatau_ba_plus -
                                v2sigmatau_ba_minus) / (delta_h * 2)
        self.check_term(tol, "v3rhosigmatau_aba", v3rhosigmatau_aba,
                        v3rhosigmatau_aba_fd)

        v3rhosigmatau_abb = kxc["v3rhosigmatau"][0][5]
        v2sigmatau_bb_plus = fxc_plus["v2sigmatau"][0][5]
        v2sigmatau_bb_minus = fxc_minus["v2sigmatau"][0][5]
        v3rhosigmatau_abb_fd = (v2sigmatau_bb_plus -
                                v2sigmatau_bb_minus) / (delta_h * 2)
        self.check_term(tol, "v3rhosigmatau_abb", v3rhosigmatau_abb,
                        v3rhosigmatau_abb_fd)

        v3rhotau2_aaa = kxc["v3rhotau2"][0][0]
        v2tau2_aa_plus = fxc_plus["v2tau2"][0][0]
        v2tau2_aa_minus = fxc_minus["v2tau2"][0][0]
        v3rhotau2_aaa_fd = (v2tau2_aa_plus - v2tau2_aa_minus) / (delta_h * 2)
        self.check_term(tol, "v3rhotau2_aaa", v3rhotau2_aaa, v3rhotau2_aaa_fd)

        v3rhotau2_aab = kxc["v3rhotau2"][0][1]
        v2tau2_ab_plus = fxc_plus["v2tau2"][0][1]
        v2tau2_ab_minus = fxc_minus["v2tau2"][0][1]
        v3rhotau2_aab_fd = (v2tau2_ab_plus - v2tau2_ab_minus) / (delta_h * 2)
        self.check_term(tol, "v3rhotau2_aab", v3rhotau2_aab, v3rhotau2_aab_fd)

        v3rhotau2_abb = kxc["v3rhotau2"][0][2]
        v2tau2_bb_plus = fxc_plus["v2tau2"][0][2]
        v2tau2_bb_minus = fxc_minus["v2tau2"][0][2]
        v3rhotau2_abb_fd = (v2tau2_bb_plus - v2tau2_bb_minus) / (delta_h * 2)
        self.check_term(tol, "v3rhotau2_abb", v3rhotau2_abb, v3rhotau2_abb_fd)

        v4rho4_aaaa = lxc["v4rho4"][0][0]
        v3rho3_aaa_plus = kxc_plus["v3rho3"][0][0]
        v3rho3_aaa_minus = kxc_minus["v3rho3"][0][0]
        v4rho4_aaaa_fd = (v3rho3_aaa_plus - v3rho3_aaa_minus) / (delta_h * 2)
        self.check_term(tol, "v4rho4_aaaa", v4rho4_aaaa, v4rho4_aaaa_fd)

        v4rho4_aaab = lxc["v4rho4"][0][1]
        v3rho3_aab_plus = kxc_plus["v3rho3"][0][1]
        v3rho3_aab_minus = kxc_minus["v3rho3"][0][1]
        v4rho4_aaab_fd = (v3rho3_aab_plus - v3rho3_aab_minus) / (delta_h * 2)
        self.check_term(tol, "v4rho4_aaab", v4rho4_aaab, v4rho4_aaab_fd)

        v4rho4_aabb = lxc["v4rho4"][0][2]
        v3rho3_abb_plus = kxc_plus["v3rho3"][0][2]
        v3rho3_abb_minus = kxc_minus["v3rho3"][0][2]
        v4rho4_aabb_fd = (v3rho3_abb_plus - v3rho3_abb_minus) / (delta_h * 2)
        self.check_term(tol, "v4rho4_aabb", v4rho4_aabb, v4rho4_aabb_fd)

        v4rho4_abbb = lxc["v4rho4"][0][3]
        v3rho3_bbb_plus = kxc_plus["v3rho3"][0][3]
        v3rho3_bbb_minus = kxc_minus["v3rho3"][0][3]
        v4rho4_abbb_fd = (v3rho3_bbb_plus - v3rho3_bbb_minus) / (delta_h * 2)
        self.check_term(tol, "v4rho4_abbb", v4rho4_abbb, v4rho4_abbb_fd)

        v4rho3sigma_aaaa = lxc["v4rho3sigma"][0][0]
        v3rho2sigma_aaa_plus = kxc_plus["v3rho2sigma"][0][0]
        v3rho2sigma_aaa_minus = kxc_minus["v3rho2sigma"][0][0]
        v4rho3sigma_aaaa_fd = (v3rho2sigma_aaa_plus -
                               v3rho2sigma_aaa_minus) / (delta_h * 2)
        self.check_term(tol, "v4rho3sigma_aaaa", v4rho3sigma_aaaa,
                        v4rho3sigma_aaaa_fd)

        v4rho3sigma_aaac = lxc["v4rho3sigma"][0][1]
        v3rho2sigma_aac_plus = kxc_plus["v3rho2sigma"][0][1]
        v3rho2sigma_aac_minus = kxc_minus["v3rho2sigma"][0][1]
        v4rho3sigma_aaac_fd = (v3rho2sigma_aac_plus -
                               v3rho2sigma_aac_minus) / (delta_h * 2)
        self.check_term(tol, "v4rho3sigma_aaac", v4rho3sigma_aaac,
                        v4rho3sigma_aaac_fd)

        v4rho3sigma_aaab = lxc["v4rho3sigma"][0][2]
        v3rho2sigma_aab_plus = kxc_plus["v3rho2sigma"][0][2]
        v3rho2sigma_aab_minus = kxc_minus["v3rho2sigma"][0][2]
        v4rho3sigma_aaab_fd = (v3rho2sigma_aab_plus -
                               v3rho2sigma_aab_minus) / (delta_h * 2)
        self.check_term(tol, "v4rho3sigma_aaab", v4rho3sigma_aaab,
                        v4rho3sigma_aaab_fd)

        v4rho3sigma_aaba = lxc["v4rho3sigma"][0][3]
        v3rho2sigma_aba_plus = kxc_plus["v3rho2sigma"][0][3]
        v3rho2sigma_aba_minus = kxc_minus["v3rho2sigma"][0][3]
        v4rho3sigma_aaba_fd = (v3rho2sigma_aba_plus -
                               v3rho2sigma_aba_minus) / (delta_h * 2)
        self.check_term(tol, "v4rho3sigma_aaba", v4rho3sigma_aaba,
                        v4rho3sigma_aaba_fd)

        v4rho3sigma_aabc = lxc["v4rho3sigma"][0][4]
        v3rho2sigma_abc_plus = kxc_plus["v3rho2sigma"][0][4]
        v3rho2sigma_abc_minus = kxc_minus["v3rho2sigma"][0][4]
        v4rho3sigma_aabc_fd = (v3rho2sigma_abc_plus -
                               v3rho2sigma_abc_minus) / (delta_h * 2)
        self.check_term(tol, "v4rho3sigma_aabc", v4rho3sigma_aabc,
                        v4rho3sigma_aabc_fd)

        v4rho3sigma_aabb = lxc["v4rho3sigma"][0][5]
        v3rho2sigma_abb_plus = kxc_plus["v3rho2sigma"][0][5]
        v3rho2sigma_abb_minus = kxc_minus["v3rho2sigma"][0][5]
        v4rho3sigma_aabb_fd = (v3rho2sigma_abb_plus -
                               v3rho2sigma_abb_minus) / (delta_h * 2)
        self.check_term(tol, "v4rho3sigma_aabb", v4rho3sigma_aabb,
                        v4rho3sigma_aabb_fd)

        v4rho3sigma_abba = lxc["v4rho3sigma"][0][6]
        v3rho2sigma_bba_plus = kxc_plus["v3rho2sigma"][0][6]
        v3rho2sigma_bba_minus = kxc_minus["v3rho2sigma"][0][6]
        v4rho3sigma_abba_fd = (v3rho2sigma_bba_plus -
                               v3rho2sigma_bba_minus) / (delta_h * 2)
        self.check_term(tol, "v4rho3sigma_abba", v4rho3sigma_abba,
                        v4rho3sigma_abba_fd)

        v4rho3sigma_abbc = lxc["v4rho3sigma"][0][7]
        v3rho2sigma_bbc_plus = kxc_plus["v3rho2sigma"][0][7]
        v3rho2sigma_bbc_minus = kxc_minus["v3rho2sigma"][0][7]
        v4rho3sigma_abbc_fd = (v3rho2sigma_bbc_plus -
                               v3rho2sigma_bbc_minus) / (delta_h * 2)
        self.check_term(tol, "v4rho3sigma_abbc", v4rho3sigma_abbc,
                        v4rho3sigma_abbc_fd)

        v4rho3sigma_abbb = lxc["v4rho3sigma"][0][8]
        v3rho2sigma_bbb_plus = kxc_plus["v3rho2sigma"][0][8]
        v3rho2sigma_bbb_minus = kxc_minus["v3rho2sigma"][0][8]
        v4rho3sigma_abbb_fd = (v3rho2sigma_bbb_plus -
                               v3rho2sigma_bbb_minus) / (delta_h * 2)
        self.check_term(tol, "v4rho3sigma_abbb", v4rho3sigma_abbb,
                        v4rho3sigma_abbb_fd)

        v4rho3tau_aaaa = lxc["v4rho3tau"][0][0]
        v3rho2tau_aaa_plus = kxc_plus["v3rho2tau"][0][0]
        v3rho2tau_aaa_minus = kxc_minus["v3rho2tau"][0][0]
        v4rho3tau_aaaa_fd = (v3rho2tau_aaa_plus -
                             v3rho2tau_aaa_minus) / (delta_h * 2)
        self.check_term(tol, "v4rho3tau_aaaa", v4rho3tau_aaaa,
                        v4rho3tau_aaaa_fd)

        v4rho3tau_aaab = lxc["v4rho3tau"][0][1]
        v3rho2tau_aab_plus = kxc_plus["v3rho2tau"][0][1]
        v3rho2tau_aab_minus = kxc_minus["v3rho2tau"][0][1]
        v4rho3tau_aaab_fd = (v3rho2tau_aab_plus -
                             v3rho2tau_aab_minus) / (delta_h * 2)
        self.check_term(tol, "v4rho3tau_aaab", v4rho3tau_aaab,
                        v4rho3tau_aaab_fd)

        v4rho3tau_aaba = lxc["v4rho3tau"][0][2]
        v3rho2tau_aba_plus = kxc_plus["v3rho2tau"][0][2]
        v3rho2tau_aba_minus = kxc_minus["v3rho2tau"][0][2]
        v4rho3tau_aaba_fd = (v3rho2tau_aba_plus -
                             v3rho2tau_aba_minus) / (delta_h * 2)
        self.check_term(tol, "v4rho3tau_aaba", v4rho3tau_aaba,
                        v4rho3tau_aaba_fd)

        v4rho3tau_aabb = lxc["v4rho3tau"][0][3]
        v3rho2tau_abb_plus = kxc_plus["v3rho2tau"][0][3]
        v3rho2tau_abb_minus = kxc_minus["v3rho2tau"][0][3]
        v4rho3tau_aabb_fd = (v3rho2tau_abb_plus -
                             v3rho2tau_abb_minus) / (delta_h * 2)
        self.check_term(tol, "v4rho3tau_aabb", v4rho3tau_aabb,
                        v4rho3tau_aabb_fd)

        v4rho3tau_abba = lxc["v4rho3tau"][0][4]
        v3rho2tau_bba_plus = kxc_plus["v3rho2tau"][0][4]
        v3rho2tau_bba_minus = kxc_minus["v3rho2tau"][0][4]
        v4rho3tau_abba_fd = (v3rho2tau_bba_plus -
                             v3rho2tau_bba_minus) / (delta_h * 2)
        self.check_term(tol, "v4rho3tau_abba", v4rho3tau_abba,
                        v4rho3tau_abba_fd)

        v4rho3tau_abbb = lxc["v4rho3tau"][0][5]
        v3rho2tau_bbb_plus = kxc_plus["v3rho2tau"][0][5]
        v3rho2tau_bbb_minus = kxc_minus["v3rho2tau"][0][5]
        v4rho3tau_abbb_fd = (v3rho2tau_bbb_plus -
                             v3rho2tau_bbb_minus) / (delta_h * 2)
        self.check_term(tol, "v4rho3tau_abbb", v4rho3tau_abbb,
                        v4rho3tau_abbb_fd)

        v4rho2sigma2_aaaa = lxc["v4rho2sigma2"][0][0]
        v3rhosigma2_aaa_plus = kxc_plus["v3rhosigma2"][0][0]
        v3rhosigma2_aaa_minus = kxc_minus["v3rhosigma2"][0][0]
        v4rho2sigma2_aaaa_fd = (v3rhosigma2_aaa_plus -
                                v3rhosigma2_aaa_minus) / (delta_h * 2)
        self.check_term(tol, "v4rho2sigma2_aaaa", v4rho2sigma2_aaaa,
                        v4rho2sigma2_aaaa_fd)

        v4rho2sigma2_aaac = lxc["v4rho2sigma2"][0][1]
        v3rhosigma2_aac_plus = kxc_plus["v3rhosigma2"][0][1]
        v3rhosigma2_aac_minus = kxc_minus["v3rhosigma2"][0][1]
        v4rho2sigma2_aaac_fd = (v3rhosigma2_aac_plus -
                                v3rhosigma2_aac_minus) / (delta_h * 2)
        self.check_term(tol, "v4rho2sigma2_aaac", v4rho2sigma2_aaac,
                        v4rho2sigma2_aaac_fd)

        v4rho2sigma2_aaab = lxc["v4rho2sigma2"][0][2]
        v3rhosigma2_aab_plus = kxc_plus["v3rhosigma2"][0][2]
        v3rhosigma2_aab_minus = kxc_minus["v3rhosigma2"][0][2]
        v4rho2sigma2_aaab_fd = (v3rhosigma2_aab_plus -
                                v3rhosigma2_aab_minus) / (delta_h * 2)
        self.check_term(tol, "v4rho2sigma2_aaab", v4rho2sigma2_aaab,
                        v4rho2sigma2_aaab_fd)

        v4rho2sigma2_aacc = lxc["v4rho2sigma2"][0][3]
        v3rhosigma2_acc_plus = kxc_plus["v3rhosigma2"][0][3]
        v3rhosigma2_acc_minus = kxc_minus["v3rhosigma2"][0][3]
        v4rho2sigma2_aacc_fd = (v3rhosigma2_acc_plus -
                                v3rhosigma2_acc_minus) / (delta_h * 2)
        self.check_term(tol, "v4rho2sigma2_aacc", v4rho2sigma2_aacc,
                        v4rho2sigma2_aacc_fd)

        v4rho2sigma2_aacb = lxc["v4rho2sigma2"][0][4]
        v3rhosigma2_acb_plus = kxc_plus["v3rhosigma2"][0][4]
        v3rhosigma2_acb_minus = kxc_minus["v3rhosigma2"][0][4]
        v4rho2sigma2_aacb_fd = (v3rhosigma2_acb_plus -
                                v3rhosigma2_acb_minus) / (delta_h * 2)
        self.check_term(tol, "v4rho2sigma2_aacb", v4rho2sigma2_aacb,
                        v4rho2sigma2_aacb_fd)

        v4rho2sigma2_aabb = lxc["v4rho2sigma2"][0][5]
        v3rhosigma2_abb_plus = kxc_plus["v3rhosigma2"][0][5]
        v3rhosigma2_abb_minus = kxc_minus["v3rhosigma2"][0][5]
        v4rho2sigma2_aabb_fd = (v3rhosigma2_abb_plus -
                                v3rhosigma2_abb_minus) / (delta_h * 2)
        self.check_term(tol, "v4rho2sigma2_aabb", v4rho2sigma2_aabb,
                        v4rho2sigma2_aabb_fd)

        v4rho2sigma2_abaa = lxc["v4rho2sigma2"][0][6]
        v3rhosigma2_baa_plus = kxc_plus["v3rhosigma2"][0][6]
        v3rhosigma2_baa_minus = kxc_minus["v3rhosigma2"][0][6]
        v4rho2sigma2_abaa_fd = (v3rhosigma2_baa_plus -
                                v3rhosigma2_baa_minus) / (delta_h * 2)
        self.check_term(tol, "v4rho2sigma2_abaa", v4rho2sigma2_abaa,
                        v4rho2sigma2_abaa_fd)

        v4rho2sigma2_abac = lxc["v4rho2sigma2"][0][7]
        v3rhosigma2_bac_plus = kxc_plus["v3rhosigma2"][0][7]
        v3rhosigma2_bac_minus = kxc_minus["v3rhosigma2"][0][7]
        v4rho2sigma2_abac_fd = (v3rhosigma2_bac_plus -
                                v3rhosigma2_bac_minus) / (delta_h * 2)
        self.check_term(tol, "v4rho2sigma2_abac", v4rho2sigma2_abac,
                        v4rho2sigma2_abac_fd)

        v4rho2sigma2_abab = lxc["v4rho2sigma2"][0][8]
        v3rhosigma2_bab_plus = kxc_plus["v3rhosigma2"][0][8]
        v3rhosigma2_bab_minus = kxc_minus["v3rhosigma2"][0][8]
        v4rho2sigma2_abab_fd = (v3rhosigma2_bab_plus -
                                v3rhosigma2_bab_minus) / (delta_h * 2)
        self.check_term(tol, "v4rho2sigma2_abab", v4rho2sigma2_abab,
                        v4rho2sigma2_abab_fd)

        v4rho2sigma2_abcc = lxc["v4rho2sigma2"][0][9]
        v3rhosigma2_bcc_plus = kxc_plus["v3rhosigma2"][0][9]
        v3rhosigma2_bcc_minus = kxc_minus["v3rhosigma2"][0][9]
        v4rho2sigma2_abcc_fd = (v3rhosigma2_bcc_plus -
                                v3rhosigma2_bcc_minus) / (delta_h * 2)
        self.check_term(tol, "v4rho2sigma2_abcc", v4rho2sigma2_abcc,
                        v4rho2sigma2_abcc_fd)

        v4rho2sigma2_abcb = lxc["v4rho2sigma2"][0][10]
        v3rhosigma2_bcb_plus = kxc_plus["v3rhosigma2"][0][10]
        v3rhosigma2_bcb_minus = kxc_minus["v3rhosigma2"][0][10]
        v4rho2sigma2_abcb_fd = (v3rhosigma2_bcb_plus -
                                v3rhosigma2_bcb_minus) / (delta_h * 2)
        self.check_term(tol, "v4rho2sigma2_abcb", v4rho2sigma2_abcb,
                        v4rho2sigma2_abcb_fd)

        v4rho2sigma2_abbb = lxc["v4rho2sigma2"][0][11]
        v3rhosigma2_bbb_plus = kxc_plus["v3rhosigma2"][0][11]
        v3rhosigma2_bbb_minus = kxc_minus["v3rhosigma2"][0][11]
        v4rho2sigma2_abbb_fd = (v3rhosigma2_bbb_plus -
                                v3rhosigma2_bbb_minus) / (delta_h * 2)
        self.check_term(tol, "v4rho2sigma2_abbb", v4rho2sigma2_abbb,
                        v4rho2sigma2_abbb_fd)

        v4rho2sigmatau_aaaa = lxc["v4rho2sigmatau"][0][0]
        v3rhosigmatau_aaa_plus = kxc_plus["v3rhosigmatau"][0][0]
        v3rhosigmatau_aaa_minus = kxc_minus["v3rhosigmatau"][0][0]
        v4rho2sigmatau_aaaa_fd = (v3rhosigmatau_aaa_plus -
                                  v3rhosigmatau_aaa_minus) / (delta_h * 2)
        self.check_term(tol, "v4rho2sigmatau_aaaa", v4rho2sigmatau_aaaa,
                        v4rho2sigmatau_aaaa_fd)

        v4rho2sigmatau_aaab = lxc["v4rho2sigmatau"][0][1]
        v3rhosigmatau_aab_plus = kxc_plus["v3rhosigmatau"][0][1]
        v3rhosigmatau_aab_minus = kxc_minus["v3rhosigmatau"][0][1]
        v4rho2sigmatau_aaab_fd = (v3rhosigmatau_aab_plus -
                                  v3rhosigmatau_aab_minus) / (delta_h * 2)
        self.check_term(tol, "v4rho2sigmatau_aaab", v4rho2sigmatau_aaab,
                        v4rho2sigmatau_aaab_fd)

        v4rho2sigmatau_aaca = lxc["v4rho2sigmatau"][0][2]
        v3rhosigmatau_aca_plus = kxc_plus["v3rhosigmatau"][0][2]
        v3rhosigmatau_aca_minus = kxc_minus["v3rhosigmatau"][0][2]
        v4rho2sigmatau_aaca_fd = (v3rhosigmatau_aca_plus -
                                  v3rhosigmatau_aca_minus) / (delta_h * 2)
        self.check_term(tol, "v4rho2sigmatau_aaca", v4rho2sigmatau_aaca,
                        v4rho2sigmatau_aaca_fd)

        v4rho2sigmatau_aacb = lxc["v4rho2sigmatau"][0][3]
        v3rhosigmatau_acb_plus = kxc_plus["v3rhosigmatau"][0][3]
        v3rhosigmatau_acb_minus = kxc_minus["v3rhosigmatau"][0][3]
        v4rho2sigmatau_aacb_fd = (v3rhosigmatau_acb_plus -
                                  v3rhosigmatau_acb_minus) / (delta_h * 2)
        self.check_term(tol, "v4rho2sigmatau_aacb", v4rho2sigmatau_aacb,
                        v4rho2sigmatau_aacb_fd)

        v4rho2sigmatau_aaba = lxc["v4rho2sigmatau"][0][4]
        v3rhosigmatau_aba_plus = kxc_plus["v3rhosigmatau"][0][4]
        v3rhosigmatau_aba_minus = kxc_minus["v3rhosigmatau"][0][4]
        v4rho2sigmatau_aaba_fd = (v3rhosigmatau_aba_plus -
                                  v3rhosigmatau_aba_minus) / (delta_h * 2)
        self.check_term(tol, "v4rho2sigmatau_aaba", v4rho2sigmatau_aaba,
                        v4rho2sigmatau_aaba_fd)

        v4rho2sigmatau_aabb = lxc["v4rho2sigmatau"][0][5]
        v3rhosigmatau_abb_plus = kxc_plus["v3rhosigmatau"][0][5]
        v3rhosigmatau_abb_minus = kxc_minus["v3rhosigmatau"][0][5]
        v4rho2sigmatau_aabb_fd = (v3rhosigmatau_abb_plus -
                                  v3rhosigmatau_abb_minus) / (delta_h * 2)
        self.check_term(tol, "v4rho2sigmatau_aabb", v4rho2sigmatau_aabb,
                        v4rho2sigmatau_aabb_fd)

        v4rho2sigmatau_abaa = lxc["v4rho2sigmatau"][0][6]
        v3rhosigmatau_baa_plus = kxc_plus["v3rhosigmatau"][0][6]
        v3rhosigmatau_baa_minus = kxc_minus["v3rhosigmatau"][0][6]
        v4rho2sigmatau_abaa_fd = (v3rhosigmatau_baa_plus -
                                  v3rhosigmatau_baa_minus) / (delta_h * 2)
        self.check_term(tol, "v4rho2sigmatau_abaa", v4rho2sigmatau_abaa,
                        v4rho2sigmatau_abaa_fd)

        v4rho2sigmatau_abab = lxc["v4rho2sigmatau"][0][7]
        v3rhosigmatau_bab_plus = kxc_plus["v3rhosigmatau"][0][7]
        v3rhosigmatau_bab_minus = kxc_minus["v3rhosigmatau"][0][7]
        v4rho2sigmatau_abab_fd = (v3rhosigmatau_bab_plus -
                                  v3rhosigmatau_bab_minus) / (delta_h * 2)
        self.check_term(tol, "v4rho2sigmatau_abab", v4rho2sigmatau_abab,
                        v4rho2sigmatau_abab_fd)

        v4rho2sigmatau_abca = lxc["v4rho2sigmatau"][0][8]
        v3rhosigmatau_bca_plus = kxc_plus["v3rhosigmatau"][0][8]
        v3rhosigmatau_bca_minus = kxc_minus["v3rhosigmatau"][0][8]
        v4rho2sigmatau_abca_fd = (v3rhosigmatau_bca_plus -
                                  v3rhosigmatau_bca_minus) / (delta_h * 2)
        self.check_term(tol, "v4rho2sigmatau_abca", v4rho2sigmatau_abca,
                        v4rho2sigmatau_abca_fd)

        v4rho2sigmatau_abcb = lxc["v4rho2sigmatau"][0][9]
        v3rhosigmatau_bcb_plus = kxc_plus["v3rhosigmatau"][0][9]
        v3rhosigmatau_bcb_minus = kxc_minus["v3rhosigmatau"][0][9]
        v4rho2sigmatau_abcb_fd = (v3rhosigmatau_bcb_plus -
                                  v3rhosigmatau_bcb_minus) / (delta_h * 2)
        self.check_term(tol, "v4rho2sigmatau_abcb", v4rho2sigmatau_abcb,
                        v4rho2sigmatau_abcb_fd)

        v4rho2sigmatau_abba = lxc["v4rho2sigmatau"][0][10]
        v3rhosigmatau_bba_plus = kxc_plus["v3rhosigmatau"][0][10]
        v3rhosigmatau_bba_minus = kxc_minus["v3rhosigmatau"][0][10]
        v4rho2sigmatau_abba_fd = (v3rhosigmatau_bba_plus -
                                  v3rhosigmatau_bba_minus) / (delta_h * 2)
        self.check_term(tol, "v4rho2sigmatau_abba", v4rho2sigmatau_abba,
                        v4rho2sigmatau_abba_fd)

        v4rho2sigmatau_abbb = lxc["v4rho2sigmatau"][0][11]
        v3rhosigmatau_bbb_plus = kxc_plus["v3rhosigmatau"][0][11]
        v3rhosigmatau_bbb_minus = kxc_minus["v3rhosigmatau"][0][11]
        v4rho2sigmatau_abbb_fd = (v3rhosigmatau_bbb_plus -
                                  v3rhosigmatau_bbb_minus) / (delta_h * 2)
        self.check_term(tol, "v4rho2sigmatau_abbb", v4rho2sigmatau_abbb,
                        v4rho2sigmatau_abbb_fd)

        v4rho2tau2_aaaa = lxc["v4rho2tau2"][0][0]
        v3rhotau2_aaa_plus = kxc_plus["v3rhotau2"][0][0]
        v3rhotau2_aaa_minus = kxc_minus["v3rhotau2"][0][0]
        v4rho2tau2_aaaa_fd = (v3rhotau2_aaa_plus -
                              v3rhotau2_aaa_minus) / (delta_h * 2)
        self.check_term(tol, "v4rho2tau2_aaaa", v4rho2tau2_aaaa,
                        v4rho2tau2_aaaa_fd)

        v4rho2tau2_aaab = lxc["v4rho2tau2"][0][1]
        v3rhotau2_aab_plus = kxc_plus["v3rhotau2"][0][1]
        v3rhotau2_aab_minus = kxc_minus["v3rhotau2"][0][1]
        v4rho2tau2_aaab_fd = (v3rhotau2_aab_plus -
                              v3rhotau2_aab_minus) / (delta_h * 2)
        self.check_term(tol, "v4rho2tau2_aaab", v4rho2tau2_aaab,
                        v4rho2tau2_aaab_fd)

        v4rho2tau2_aabb = lxc["v4rho2tau2"][0][2]
        v3rhotau2_abb_plus = kxc_plus["v3rhotau2"][0][2]
        v3rhotau2_abb_minus = kxc_minus["v3rhotau2"][0][2]
        v4rho2tau2_aabb_fd = (v3rhotau2_abb_plus -
                              v3rhotau2_abb_minus) / (delta_h * 2)
        self.check_term(tol, "v4rho2tau2_aabb", v4rho2tau2_aabb,
                        v4rho2tau2_aabb_fd)

        v4rho2tau2_abaa = lxc["v4rho2tau2"][0][3]
        v3rhotau2_baa_plus = kxc_plus["v3rhotau2"][0][3]
        v3rhotau2_baa_minus = kxc_minus["v3rhotau2"][0][3]
        v4rho2tau2_abaa_fd = (v3rhotau2_baa_plus -
                              v3rhotau2_baa_minus) / (delta_h * 2)
        self.check_term(tol, "v4rho2tau2_abaa", v4rho2tau2_abaa,
                        v4rho2tau2_abaa_fd)

        v4rho2tau2_abab = lxc["v4rho2tau2"][0][4]
        v3rhotau2_bab_plus = kxc_plus["v3rhotau2"][0][4]
        v3rhotau2_bab_minus = kxc_minus["v3rhotau2"][0][4]
        v4rho2tau2_abab_fd = (v3rhotau2_bab_plus -
                              v3rhotau2_bab_minus) / (delta_h * 2)
        self.check_term(tol, "v4rho2tau2_abab", v4rho2tau2_abab,
                        v4rho2tau2_abab_fd)

        v4rho2tau2_abbb = lxc["v4rho2tau2"][0][5]
        v3rhotau2_bbb_plus = kxc_plus["v3rhotau2"][0][5]
        v3rhotau2_bbb_minus = kxc_minus["v3rhotau2"][0][5]
        v4rho2tau2_abbb_fd = (v3rhotau2_bbb_plus -
                              v3rhotau2_bbb_minus) / (delta_h * 2)
        self.check_term(tol, "v4rho2tau2_abbb", v4rho2tau2_abbb,
                        v4rho2tau2_abbb_fd)

        v4rhosigma3_aaaa = lxc["v4rhosigma3"][0][0]
        v3sigma3_aaa_plus = kxc_plus["v3sigma3"][0][0]
        v3sigma3_aaa_minus = kxc_minus["v3sigma3"][0][0]
        v4rhosigma3_aaaa_fd = (v3sigma3_aaa_plus -
                               v3sigma3_aaa_minus) / (delta_h * 2)
        self.check_term(tol, "v4rhosigma3_aaaa", v4rhosigma3_aaaa,
                        v4rhosigma3_aaaa_fd)

        v4rhosigma3_aaac = lxc["v4rhosigma3"][0][1]
        v3sigma3_aac_plus = kxc_plus["v3sigma3"][0][1]
        v3sigma3_aac_minus = kxc_minus["v3sigma3"][0][1]
        v4rhosigma3_aaac_fd = (v3sigma3_aac_plus -
                               v3sigma3_aac_minus) / (delta_h * 2)
        self.check_term(tol, "v4rhosigma3_aaac", v4rhosigma3_aaac,
                        v4rhosigma3_aaac_fd)

        v4rhosigma3_aaab = lxc["v4rhosigma3"][0][2]
        v3sigma3_aab_plus = kxc_plus["v3sigma3"][0][2]
        v3sigma3_aab_minus = kxc_minus["v3sigma3"][0][2]
        v4rhosigma3_aaab_fd = (v3sigma3_aab_plus -
                               v3sigma3_aab_minus) / (delta_h * 2)
        self.check_term(tol, "v4rhosigma3_aaab", v4rhosigma3_aaab,
                        v4rhosigma3_aaab_fd)

        v4rhosigma3_aacc = lxc["v4rhosigma3"][0][3]
        v3sigma3_acc_plus = kxc_plus["v3sigma3"][0][3]
        v3sigma3_acc_minus = kxc_minus["v3sigma3"][0][3]
        v4rhosigma3_aacc_fd = (v3sigma3_acc_plus -
                               v3sigma3_acc_minus) / (delta_h * 2)
        self.check_term(tol, "v4rhosigma3_aacc", v4rhosigma3_aacc,
                        v4rhosigma3_aacc_fd)

        v4rhosigma3_aacb = lxc["v4rhosigma3"][0][4]
        v3sigma3_acb_plus = kxc_plus["v3sigma3"][0][4]
        v3sigma3_acb_minus = kxc_minus["v3sigma3"][0][4]
        v4rhosigma3_aacb_fd = (v3sigma3_acb_plus -
                               v3sigma3_acb_minus) / (delta_h * 2)
        self.check_term(tol, "v4rhosigma3_aacb", v4rhosigma3_aacb,
                        v4rhosigma3_aacb_fd)

        v4rhosigma3_aabb = lxc["v4rhosigma3"][0][5]
        v3sigma3_abb_plus = kxc_plus["v3sigma3"][0][5]
        v3sigma3_abb_minus = kxc_minus["v3sigma3"][0][5]
        v4rhosigma3_aabb_fd = (v3sigma3_abb_plus -
                               v3sigma3_abb_minus) / (delta_h * 2)
        self.check_term(tol, "v4rhosigma3_aabb", v4rhosigma3_aabb,
                        v4rhosigma3_aabb_fd)

        v4rhosigma3_accc = lxc["v4rhosigma3"][0][6]
        v3sigma3_ccc_plus = kxc_plus["v3sigma3"][0][6]
        v3sigma3_ccc_minus = kxc_minus["v3sigma3"][0][6]
        v4rhosigma3_accc_fd = (v3sigma3_ccc_plus -
                               v3sigma3_ccc_minus) / (delta_h * 2)
        self.check_term(tol, "v4rhosigma3_accc", v4rhosigma3_accc,
                        v4rhosigma3_accc_fd)

        v4rhosigma3_accb = lxc["v4rhosigma3"][0][7]
        v3sigma3_ccb_plus = kxc_plus["v3sigma3"][0][7]
        v3sigma3_ccb_minus = kxc_minus["v3sigma3"][0][7]
        v4rhosigma3_accb_fd = (v3sigma3_ccb_plus -
                               v3sigma3_ccb_minus) / (delta_h * 2)
        self.check_term(tol, "v4rhosigma3_accb", v4rhosigma3_accb,
                        v4rhosigma3_accb_fd)

        v4rhosigma3_acbb = lxc["v4rhosigma3"][0][8]
        v3sigma3_cbb_plus = kxc_plus["v3sigma3"][0][8]
        v3sigma3_cbb_minus = kxc_minus["v3sigma3"][0][8]
        v4rhosigma3_acbb_fd = (v3sigma3_cbb_plus -
                               v3sigma3_cbb_minus) / (delta_h * 2)
        self.check_term(tol, "v4rhosigma3_acbb", v4rhosigma3_acbb,
                        v4rhosigma3_acbb_fd)

        v4rhosigma3_abbb = lxc["v4rhosigma3"][0][9]
        v3sigma3_bbb_plus = kxc_plus["v3sigma3"][0][9]
        v3sigma3_bbb_minus = kxc_minus["v3sigma3"][0][9]
        v4rhosigma3_abbb_fd = (v3sigma3_bbb_plus -
                               v3sigma3_bbb_minus) / (delta_h * 2)
        self.check_term(tol, "v4rhosigma3_abbb", v4rhosigma3_abbb,
                        v4rhosigma3_abbb_fd)

        v4rhosigma2tau_aaaa = lxc["v4rhosigma2tau"][0][0]
        v3sigma2tau_aaa_plus = kxc_plus["v3sigma2tau"][0][0]
        v3sigma2tau_aaa_minus = kxc_minus["v3sigma2tau"][0][0]
        v4rhosigma2tau_aaaa_fd = (v3sigma2tau_aaa_plus -
                                  v3sigma2tau_aaa_minus) / (delta_h * 2)
        self.check_term(tol, "v4rhosigma2tau_aaaa", v4rhosigma2tau_aaaa,
                        v4rhosigma2tau_aaaa_fd)

        v4rhosigma2tau_aaab = lxc["v4rhosigma2tau"][0][1]
        v3sigma2tau_aab_plus = kxc_plus["v3sigma2tau"][0][1]
        v3sigma2tau_aab_minus = kxc_minus["v3sigma2tau"][0][1]
        v4rhosigma2tau_aaab_fd = (v3sigma2tau_aab_plus -
                                  v3sigma2tau_aab_minus) / (delta_h * 2)
        self.check_term(tol, "v4rhosigma2tau_aaab", v4rhosigma2tau_aaab,
                        v4rhosigma2tau_aaab_fd)

        v4rhosigma2tau_aaca = lxc["v4rhosigma2tau"][0][2]
        v3sigma2tau_aca_plus = kxc_plus["v3sigma2tau"][0][2]
        v3sigma2tau_aca_minus = kxc_minus["v3sigma2tau"][0][2]
        v4rhosigma2tau_aaca_fd = (v3sigma2tau_aca_plus -
                                  v3sigma2tau_aca_minus) / (delta_h * 2)
        self.check_term(tol, "v4rhosigma2tau_aaca", v4rhosigma2tau_aaca,
                        v4rhosigma2tau_aaca_fd)

        v4rhosigma2tau_aacb = lxc["v4rhosigma2tau"][0][3]
        v3sigma2tau_acb_plus = kxc_plus["v3sigma2tau"][0][3]
        v3sigma2tau_acb_minus = kxc_minus["v3sigma2tau"][0][3]
        v4rhosigma2tau_aacb_fd = (v3sigma2tau_acb_plus -
                                  v3sigma2tau_acb_minus) / (delta_h * 2)
        self.check_term(tol, "v4rhosigma2tau_aacb", v4rhosigma2tau_aacb,
                        v4rhosigma2tau_aacb_fd)

        v4rhosigma2tau_aaba = lxc["v4rhosigma2tau"][0][4]
        v3sigma2tau_aba_plus = kxc_plus["v3sigma2tau"][0][4]
        v3sigma2tau_aba_minus = kxc_minus["v3sigma2tau"][0][4]
        v4rhosigma2tau_aaba_fd = (v3sigma2tau_aba_plus -
                                  v3sigma2tau_aba_minus) / (delta_h * 2)
        self.check_term(tol, "v4rhosigma2tau_aaba", v4rhosigma2tau_aaba,
                        v4rhosigma2tau_aaba_fd)

        v4rhosigma2tau_aabb = lxc["v4rhosigma2tau"][0][5]
        v3sigma2tau_abb_plus = kxc_plus["v3sigma2tau"][0][5]
        v3sigma2tau_abb_minus = kxc_minus["v3sigma2tau"][0][5]
        v4rhosigma2tau_aabb_fd = (v3sigma2tau_abb_plus -
                                  v3sigma2tau_abb_minus) / (delta_h * 2)
        self.check_term(tol, "v4rhosigma2tau_aabb", v4rhosigma2tau_aabb,
                        v4rhosigma2tau_aabb_fd)

        v4rhosigma2tau_acca = lxc["v4rhosigma2tau"][0][6]
        v3sigma2tau_cca_plus = kxc_plus["v3sigma2tau"][0][6]
        v3sigma2tau_cca_minus = kxc_minus["v3sigma2tau"][0][6]
        v4rhosigma2tau_acca_fd = (v3sigma2tau_cca_plus -
                                  v3sigma2tau_cca_minus) / (delta_h * 2)
        self.check_term(tol, "v4rhosigma2tau_acca", v4rhosigma2tau_acca,
                        v4rhosigma2tau_acca_fd)

        v4rhosigma2tau_accb = lxc["v4rhosigma2tau"][0][7]
        v3sigma2tau_ccb_plus = kxc_plus["v3sigma2tau"][0][7]
        v3sigma2tau_ccb_minus = kxc_minus["v3sigma2tau"][0][7]
        v4rhosigma2tau_accb_fd = (v3sigma2tau_ccb_plus -
                                  v3sigma2tau_ccb_minus) / (delta_h * 2)
        self.check_term(tol, "v4rhosigma2tau_accb", v4rhosigma2tau_accb,
                        v4rhosigma2tau_accb_fd)

        v4rhosigma2tau_acba = lxc["v4rhosigma2tau"][0][8]
        v3sigma2tau_cba_plus = kxc_plus["v3sigma2tau"][0][8]
        v3sigma2tau_cba_minus = kxc_minus["v3sigma2tau"][0][8]
        v4rhosigma2tau_acba_fd = (v3sigma2tau_cba_plus -
                                  v3sigma2tau_cba_minus) / (delta_h * 2)
        self.check_term(tol, "v4rhosigma2tau_acba", v4rhosigma2tau_acba,
                        v4rhosigma2tau_acba_fd)

        v4rhosigma2tau_acbb = lxc["v4rhosigma2tau"][0][9]
        v3sigma2tau_cbb_plus = kxc_plus["v3sigma2tau"][0][9]
        v3sigma2tau_cbb_minus = kxc_minus["v3sigma2tau"][0][9]
        v4rhosigma2tau_acbb_fd = (v3sigma2tau_cbb_plus -
                                  v3sigma2tau_cbb_minus) / (delta_h * 2)
        self.check_term(tol, "v4rhosigma2tau_acbb", v4rhosigma2tau_acbb,
                        v4rhosigma2tau_acbb_fd)

        v4rhosigma2tau_abba = lxc["v4rhosigma2tau"][0][10]
        v3sigma2tau_bba_plus = kxc_plus["v3sigma2tau"][0][10]
        v3sigma2tau_bba_minus = kxc_minus["v3sigma2tau"][0][10]
        v4rhosigma2tau_abba_fd = (v3sigma2tau_bba_plus -
                                  v3sigma2tau_bba_minus) / (delta_h * 2)
        self.check_term(tol, "v4rhosigma2tau_abba", v4rhosigma2tau_abba,
                        v4rhosigma2tau_abba_fd)

        v4rhosigma2tau_abbb = lxc["v4rhosigma2tau"][0][11]
        v3sigma2tau_bbb_plus = kxc_plus["v3sigma2tau"][0][11]
        v3sigma2tau_bbb_minus = kxc_minus["v3sigma2tau"][0][11]
        v4rhosigma2tau_abbb_fd = (v3sigma2tau_bbb_plus -
                                  v3sigma2tau_bbb_minus) / (delta_h * 2)
        self.check_term(tol, "v4rhosigma2tau_abbb", v4rhosigma2tau_abbb,
                        v4rhosigma2tau_abbb_fd)

        v4rhosigmatau2_aaaa = lxc["v4rhosigmatau2"][0][0]
        v3sigmatau2_aaa_plus = kxc_plus["v3sigmatau2"][0][0]
        v3sigmatau2_aaa_minus = kxc_minus["v3sigmatau2"][0][0]
        v4rhosigmatau2_aaaa_fd = (v3sigmatau2_aaa_plus -
                                  v3sigmatau2_aaa_minus) / (delta_h * 2)
        self.check_term(tol, "v4rhosigmatau2_aaaa", v4rhosigmatau2_aaaa,
                        v4rhosigmatau2_aaaa_fd)

        v4rhosigmatau2_aaab = lxc["v4rhosigmatau2"][0][1]
        v3sigmatau2_aab_plus = kxc_plus["v3sigmatau2"][0][1]
        v3sigmatau2_aab_minus = kxc_minus["v3sigmatau2"][0][1]
        v4rhosigmatau2_aaab_fd = (v3sigmatau2_aab_plus -
                                  v3sigmatau2_aab_minus) / (delta_h * 2)
        self.check_term(tol, "v4rhosigmatau2_aaab", v4rhosigmatau2_aaab,
                        v4rhosigmatau2_aaab_fd)

        v4rhosigmatau2_aabb = lxc["v4rhosigmatau2"][0][2]
        v3sigmatau2_abb_plus = kxc_plus["v3sigmatau2"][0][2]
        v3sigmatau2_abb_minus = kxc_minus["v3sigmatau2"][0][2]
        v4rhosigmatau2_aabb_fd = (v3sigmatau2_abb_plus -
                                  v3sigmatau2_abb_minus) / (delta_h * 2)
        self.check_term(tol, "v4rhosigmatau2_aabb", v4rhosigmatau2_aabb,
                        v4rhosigmatau2_aabb_fd)

        v4rhosigmatau2_acaa = lxc["v4rhosigmatau2"][0][3]
        v3sigmatau2_caa_plus = kxc_plus["v3sigmatau2"][0][3]
        v3sigmatau2_caa_minus = kxc_minus["v3sigmatau2"][0][3]
        v4rhosigmatau2_acaa_fd = (v3sigmatau2_caa_plus -
                                  v3sigmatau2_caa_minus) / (delta_h * 2)
        self.check_term(tol, "v4rhosigmatau2_acaa", v4rhosigmatau2_acaa,
                        v4rhosigmatau2_acaa_fd)

        v4rhosigmatau2_acab = lxc["v4rhosigmatau2"][0][4]
        v3sigmatau2_cab_plus = kxc_plus["v3sigmatau2"][0][4]
        v3sigmatau2_cab_minus = kxc_minus["v3sigmatau2"][0][4]
        v4rhosigmatau2_acab_fd = (v3sigmatau2_cab_plus -
                                  v3sigmatau2_cab_minus) / (delta_h * 2)
        self.check_term(tol, "v4rhosigmatau2_acab", v4rhosigmatau2_acab,
                        v4rhosigmatau2_acab_fd)

        v4rhosigmatau2_acbb = lxc["v4rhosigmatau2"][0][5]
        v3sigmatau2_cbb_plus = kxc_plus["v3sigmatau2"][0][5]
        v3sigmatau2_cbb_minus = kxc_minus["v3sigmatau2"][0][5]
        v4rhosigmatau2_acbb_fd = (v3sigmatau2_cbb_plus -
                                  v3sigmatau2_cbb_minus) / (delta_h * 2)
        self.check_term(tol, "v4rhosigmatau2_acbb", v4rhosigmatau2_acbb,
                        v4rhosigmatau2_acbb_fd)

        v4rhosigmatau2_abaa = lxc["v4rhosigmatau2"][0][6]
        v3sigmatau2_baa_plus = kxc_plus["v3sigmatau2"][0][6]
        v3sigmatau2_baa_minus = kxc_minus["v3sigmatau2"][0][6]
        v4rhosigmatau2_abaa_fd = (v3sigmatau2_baa_plus -
                                  v3sigmatau2_baa_minus) / (delta_h * 2)
        self.check_term(tol, "v4rhosigmatau2_abaa", v4rhosigmatau2_abaa,
                        v4rhosigmatau2_abaa_fd)

        v4rhosigmatau2_abab = lxc["v4rhosigmatau2"][0][7]
        v3sigmatau2_bab_plus = kxc_plus["v3sigmatau2"][0][7]
        v3sigmatau2_bab_minus = kxc_minus["v3sigmatau2"][0][7]
        v4rhosigmatau2_abab_fd = (v3sigmatau2_bab_plus -
                                  v3sigmatau2_bab_minus) / (delta_h * 2)
        self.check_term(tol, "v4rhosigmatau2_abab", v4rhosigmatau2_abab,
                        v4rhosigmatau2_abab_fd)

        v4rhosigmatau2_abbb = lxc["v4rhosigmatau2"][0][8]
        v3sigmatau2_bbb_plus = kxc_plus["v3sigmatau2"][0][8]
        v3sigmatau2_bbb_minus = kxc_minus["v3sigmatau2"][0][8]
        v4rhosigmatau2_abbb_fd = (v3sigmatau2_bbb_plus -
                                  v3sigmatau2_bbb_minus) / (delta_h * 2)
        self.check_term(tol, "v4rhosigmatau2_abbb", v4rhosigmatau2_abbb,
                        v4rhosigmatau2_abbb_fd)

        v4rhotau3_aaaa = lxc["v4rhotau3"][0][0]
        v3tau3_aaa_plus = kxc_plus["v3tau3"][0][0]
        v3tau3_aaa_minus = kxc_minus["v3tau3"][0][0]
        v4rhotau3_aaaa_fd = (v3tau3_aaa_plus - v3tau3_aaa_minus) / (delta_h * 2)
        self.check_term(tol, "v4rhotau3_aaaa", v4rhotau3_aaaa,
                        v4rhotau3_aaaa_fd)

        v4rhotau3_aaab = lxc["v4rhotau3"][0][1]
        v3tau3_aab_plus = kxc_plus["v3tau3"][0][1]
        v3tau3_aab_minus = kxc_minus["v3tau3"][0][1]
        v4rhotau3_aaab_fd = (v3tau3_aab_plus - v3tau3_aab_minus) / (delta_h * 2)
        self.check_term(tol, "v4rhotau3_aaab", v4rhotau3_aaab,
                        v4rhotau3_aaab_fd)

        v4rhotau3_aabb = lxc["v4rhotau3"][0][2]
        v3tau3_abb_plus = kxc_plus["v3tau3"][0][2]
        v3tau3_abb_minus = kxc_minus["v3tau3"][0][2]
        v4rhotau3_aabb_fd = (v3tau3_abb_plus - v3tau3_abb_minus) / (delta_h * 2)
        self.check_term(tol, "v4rhotau3_aabb", v4rhotau3_aabb,
                        v4rhotau3_aabb_fd)

        v4rhotau3_abbb = lxc["v4rhotau3"][0][3]
        v3tau3_bbb_plus = kxc_plus["v3tau3"][0][3]
        v3tau3_bbb_minus = kxc_minus["v3tau3"][0][3]
        v4rhotau3_abbb_fd = (v3tau3_bbb_plus - v3tau3_bbb_minus) / (delta_h * 2)
        self.check_term(tol, "v4rhotau3_abbb", v4rhotau3_abbb,
                        v4rhotau3_abbb_fd)

    def test_rho_a_fd_pkzb(self):

        self.run_rho_a_fd('pkzb', 1.0e-7)

    def test_rho_a_fd_tpssh(self):

        self.run_rho_a_fd('tpssh', 1.0e-7)

    def test_rho_a_fd_m06(self):

        self.run_rho_a_fd('m06', 1.0e-7)

    def test_rho_a_fd_scan(self):

        self.run_rho_a_fd('scan', 1.0e-7)
