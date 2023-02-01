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
        nabla_b_plus = nabla_b + delta_h
        inp_plus['sigma'][1] = nabla_b_plus * nabla_a
        inp_plus['sigma'][2] = nabla_b_plus * nabla_b_plus

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
        nabla_b_minus = nabla_b - delta_h
        inp_minus['sigma'][1] = nabla_b_minus * nabla_a
        inp_minus['sigma'][2] = nabla_b_minus * nabla_b_minus

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

        vsigma_b = vxc['vsigma'][0][2]
        vsigma_c = vxc['vsigma'][0][1]
        e_a_plus = exc_plus['zk'][0][0] * (rho_a + rho_b)
        e_a_minus = exc_minus['zk'][0][0] * (rho_a + rho_b)
        vsigma_b_fd = (e_a_plus - e_a_minus) / (delta_h * 2)
        self.check_term(tol, 'vnabla_b',
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

        v2sigmatau_ba = fxc["v2sigmatau"][0][4]
        v2sigmatau_ca = fxc["v2sigmatau"][0][2]
        vtau_a_plus = vxc_plus["vtau"][0][0]
        vtau_a_minus = vxc_minus["vtau"][0][0]
        v2sigmatau_ba_fd = (vtau_a_plus - vtau_a_minus) / (delta_h * 2)
        self.check_term(tol, "v2sigmatau_ba",
                        2 * v2sigmatau_ba * nabla_b + v2sigmatau_ca * nabla_a,
                        v2sigmatau_ba_fd)

        v2sigmatau_bb = fxc["v2sigmatau"][0][5]
        v2sigmatau_cb = fxc["v2sigmatau"][0][3]
        vtau_b_plus = vxc_plus["vtau"][0][1]
        vtau_b_minus = vxc_minus["vtau"][0][1]
        v2sigmatau_bb_fd = (vtau_b_plus - vtau_b_minus) / (delta_h * 2)
        self.check_term(tol, "v2sigmatau_bb",
                        2 * v2sigmatau_bb * nabla_b + v2sigmatau_cb * nabla_a,
                        v2sigmatau_bb_fd)

        v3sigma3_bbb = kxc["v3sigma3"][0][9]
        v3sigma3_cbb = kxc["v3sigma3"][0][8]
        v2sigma2_bb_plus = fxc_plus["v2sigma2"][0][5]
        v2sigma2_bb_minus = fxc_minus["v2sigma2"][0][5]
        v3sigma3_bbb_fd = (v2sigma2_bb_plus - v2sigma2_bb_minus) / (delta_h * 2)
        self.check_term(tol, "v3sigma3_bbb",
                        2 * v3sigma3_bbb * nabla_b + v3sigma3_cbb * nabla_a,
                        v3sigma3_bbb_fd)

        v3sigma2tau_bba = kxc["v3sigma2tau"][0][10]
        v3sigma2tau_cba = kxc["v3sigma2tau"][0][8]
        v2sigmatau_ba_plus = fxc_plus["v2sigmatau"][0][4]
        v2sigmatau_ba_minus = fxc_minus["v2sigmatau"][0][4]
        v3sigma2tau_bba_fd = (v2sigmatau_ba_plus -
                              v2sigmatau_ba_minus) / (delta_h * 2)
        self.check_term(
            tol, "v3sigma2tau_bba",
            2 * v3sigma2tau_bba * nabla_b + v3sigma2tau_cba * nabla_a,
            v3sigma2tau_bba_fd)

        v3sigma2tau_bbb = kxc["v3sigma2tau"][0][11]
        v3sigma2tau_cbb = kxc["v3sigma2tau"][0][9]
        v2sigmatau_bb_plus = fxc_plus["v2sigmatau"][0][5]
        v2sigmatau_bb_minus = fxc_minus["v2sigmatau"][0][5]
        v3sigma2tau_bbb_fd = (v2sigmatau_bb_plus -
                              v2sigmatau_bb_minus) / (delta_h * 2)
        self.check_term(
            tol, "v3sigma2tau_bbb",
            2 * v3sigma2tau_bbb * nabla_b + v3sigma2tau_cbb * nabla_a,
            v3sigma2tau_bbb_fd)

        v3sigmatau2_baa = kxc["v3sigmatau2"][0][6]
        v3sigmatau2_caa = kxc["v3sigmatau2"][0][3]
        v2tau2_aa_plus = fxc_plus["v2tau2"][0][0]
        v2tau2_aa_minus = fxc_minus["v2tau2"][0][0]
        v3sigmatau2_baa_fd = (v2tau2_aa_plus - v2tau2_aa_minus) / (delta_h * 2)
        self.check_term(
            tol, "v3sigmatau2_baa",
            2 * v3sigmatau2_baa * nabla_b + v3sigmatau2_caa * nabla_a,
            v3sigmatau2_baa_fd)

        v3sigmatau2_bab = kxc["v3sigmatau2"][0][7]
        v3sigmatau2_cab = kxc["v3sigmatau2"][0][4]
        v2tau2_ab_plus = fxc_plus["v2tau2"][0][1]
        v2tau2_ab_minus = fxc_minus["v2tau2"][0][1]
        v3sigmatau2_bab_fd = (v2tau2_ab_plus - v2tau2_ab_minus) / (delta_h * 2)
        self.check_term(
            tol, "v3sigmatau2_bab",
            2 * v3sigmatau2_bab * nabla_b + v3sigmatau2_cab * nabla_a,
            v3sigmatau2_bab_fd)

        v3sigmatau2_bbb = kxc["v3sigmatau2"][0][8]
        v3sigmatau2_cbb = kxc["v3sigmatau2"][0][5]
        v2tau2_bb_plus = fxc_plus["v2tau2"][0][2]
        v2tau2_bb_minus = fxc_minus["v2tau2"][0][2]
        v3sigmatau2_bbb_fd = (v2tau2_bb_plus - v2tau2_bb_minus) / (delta_h * 2)
        self.check_term(
            tol, "v3sigmatau2_bbb",
            2 * v3sigmatau2_bbb * nabla_b + v3sigmatau2_cbb * nabla_a,
            v3sigmatau2_bbb_fd)

        v4sigma4_bbbb = lxc["v4sigma4"][0][14]
        v4sigma4_cbbb = lxc["v4sigma4"][0][13]
        v3sigma3_bbb_plus = kxc_plus["v3sigma3"][0][9]
        v3sigma3_bbb_minus = kxc_minus["v3sigma3"][0][9]
        v4sigma4_bbbb_fd = (v3sigma3_bbb_plus - v3sigma3_bbb_minus) / (delta_h *
                                                                       2)
        self.check_term(tol, "v4sigma4_bbbb",
                        2 * v4sigma4_bbbb * nabla_b + v4sigma4_cbbb * nabla_a,
                        v4sigma4_bbbb_fd)

        v4sigma3tau_bbba = lxc["v4sigma3tau"][0][18]
        v4sigma3tau_cbba = lxc["v4sigma3tau"][0][16]
        v3sigma2tau_bba_plus = kxc_plus["v3sigma2tau"][0][10]
        v3sigma2tau_bba_minus = kxc_minus["v3sigma2tau"][0][10]
        v4sigma3tau_bbba_fd = (v3sigma2tau_bba_plus -
                               v3sigma2tau_bba_minus) / (delta_h * 2)
        self.check_term(
            tol, "v4sigma3tau_bbba",
            2 * v4sigma3tau_bbba * nabla_b + v4sigma3tau_cbba * nabla_a,
            v4sigma3tau_bbba_fd)

        v4sigma3tau_bbbb = lxc["v4sigma3tau"][0][19]
        v4sigma3tau_cbbb = lxc["v4sigma3tau"][0][17]
        v3sigma2tau_bbb_plus = kxc_plus["v3sigma2tau"][0][11]
        v3sigma2tau_bbb_minus = kxc_minus["v3sigma2tau"][0][11]
        v4sigma3tau_bbbb_fd = (v3sigma2tau_bbb_plus -
                               v3sigma2tau_bbb_minus) / (delta_h * 2)
        self.check_term(
            tol, "v4sigma3tau_bbbb",
            2 * v4sigma3tau_bbbb * nabla_b + v4sigma3tau_cbbb * nabla_a,
            v4sigma3tau_bbbb_fd)

        v4sigma2tau2_bbaa = lxc["v4sigma2tau2"][0][15]
        v4sigma2tau2_cbaa = lxc["v4sigma2tau2"][0][12]
        v3sigmatau2_baa_plus = kxc_plus["v3sigmatau2"][0][6]
        v3sigmatau2_baa_minus = kxc_minus["v3sigmatau2"][0][6]
        v4sigma2tau2_bbaa_fd = (v3sigmatau2_baa_plus -
                                v3sigmatau2_baa_minus) / (delta_h * 2)
        self.check_term(
            tol, "v4sigma2tau2_bbaa",
            2 * v4sigma2tau2_bbaa * nabla_b + v4sigma2tau2_cbaa * nabla_a,
            v4sigma2tau2_bbaa_fd)

        v4sigma2tau2_bbab = lxc["v4sigma2tau2"][0][16]
        v4sigma2tau2_cbab = lxc["v4sigma2tau2"][0][13]
        v3sigmatau2_bab_plus = kxc_plus["v3sigmatau2"][0][7]
        v3sigmatau2_bab_minus = kxc_minus["v3sigmatau2"][0][7]
        v4sigma2tau2_bbab_fd = (v3sigmatau2_bab_plus -
                                v3sigmatau2_bab_minus) / (delta_h * 2)
        self.check_term(
            tol, "v4sigma2tau2_bbab",
            2 * v4sigma2tau2_bbab * nabla_b + v4sigma2tau2_cbab * nabla_a,
            v4sigma2tau2_bbab_fd)

        v4sigma2tau2_bbbb = lxc["v4sigma2tau2"][0][17]
        v4sigma2tau2_cbbb = lxc["v4sigma2tau2"][0][14]
        v3sigmatau2_bbb_plus = kxc_plus["v3sigmatau2"][0][8]
        v3sigmatau2_bbb_minus = kxc_minus["v3sigmatau2"][0][8]
        v4sigma2tau2_bbbb_fd = (v3sigmatau2_bbb_plus -
                                v3sigmatau2_bbb_minus) / (delta_h * 2)
        self.check_term(
            tol, "v4sigma2tau2_bbbb",
            2 * v4sigma2tau2_bbbb * nabla_b + v4sigma2tau2_cbbb * nabla_a,
            v4sigma2tau2_bbbb_fd)

        v4sigmatau3_baaa = lxc["v4sigmatau3"][0][8]
        v4sigmatau3_caaa = lxc["v4sigmatau3"][0][4]
        v3tau3_aaa_plus = kxc_plus["v3tau3"][0][0]
        v3tau3_aaa_minus = kxc_minus["v3tau3"][0][0]
        v4sigmatau3_baaa_fd = (v3tau3_aaa_plus - v3tau3_aaa_minus) / (delta_h *
                                                                      2)
        self.check_term(
            tol, "v4sigmatau3_baaa",
            2 * v4sigmatau3_baaa * nabla_b + v4sigmatau3_caaa * nabla_a,
            v4sigmatau3_baaa_fd)

        v4sigmatau3_baab = lxc["v4sigmatau3"][0][9]
        v4sigmatau3_caab = lxc["v4sigmatau3"][0][5]
        v3tau3_aab_plus = kxc_plus["v3tau3"][0][1]
        v3tau3_aab_minus = kxc_minus["v3tau3"][0][1]
        v4sigmatau3_baab_fd = (v3tau3_aab_plus - v3tau3_aab_minus) / (delta_h *
                                                                      2)
        self.check_term(
            tol, "v4sigmatau3_baab",
            2 * v4sigmatau3_baab * nabla_b + v4sigmatau3_caab * nabla_a,
            v4sigmatau3_baab_fd)

        v4sigmatau3_babb = lxc["v4sigmatau3"][0][10]
        v4sigmatau3_cabb = lxc["v4sigmatau3"][0][6]
        v3tau3_abb_plus = kxc_plus["v3tau3"][0][2]
        v3tau3_abb_minus = kxc_minus["v3tau3"][0][2]
        v4sigmatau3_babb_fd = (v3tau3_abb_plus - v3tau3_abb_minus) / (delta_h *
                                                                      2)
        self.check_term(
            tol, "v4sigmatau3_babb",
            2 * v4sigmatau3_babb * nabla_b + v4sigmatau3_cabb * nabla_a,
            v4sigmatau3_babb_fd)

        v4sigmatau3_bbbb = lxc["v4sigmatau3"][0][11]
        v4sigmatau3_cbbb = lxc["v4sigmatau3"][0][7]
        v3tau3_bbb_plus = kxc_plus["v3tau3"][0][3]
        v3tau3_bbb_minus = kxc_minus["v3tau3"][0][3]
        v4sigmatau3_bbbb_fd = (v3tau3_bbb_plus - v3tau3_bbb_minus) / (delta_h *
                                                                      2)
        self.check_term(
            tol, "v4sigmatau3_bbbb",
            2 * v4sigmatau3_bbbb * nabla_b + v4sigmatau3_cbbb * nabla_a,
            v4sigmatau3_bbbb_fd)

    def test_sigma_b_fd_pkzb(self):

        self.run_sigma_b_fd('pkzb', 1.0e-6)

    def test_sigma_b_fd_tpssh(self):

        self.run_sigma_b_fd('tpssh', 1.0e-7)

    def test_sigma_b_fd_m06(self):

        self.run_sigma_b_fd('m06', 1.0e-6)

    def test_sigma_b_fd_scan(self):

        self.run_sigma_b_fd('scan', 1.0e-7)
