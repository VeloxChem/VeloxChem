import numpy as np

from veloxchem.veloxchemlib import XCNewIntegrator


class TestLibxcTauA:

    def check_term(self, tol, label, val, val_fd):

        if abs(val_fd) < tol:
            assert abs(val - val_fd) < tol
        else:
            assert abs(val - val_fd) / abs(val_fd) < tol

    def run_tau_a_fd(self, xcfun_label, tol):

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
        inp_plus['tau'][0] = tau_a + delta_h

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
        inp_minus['tau'][0] = tau_a - delta_h

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

        vtau_a = vxc['vtau'][0][0]
        e_a_plus = exc_plus['zk'][0][0] * (rho_a + rho_b)
        e_a_minus = exc_minus['zk'][0][0] * (rho_a + rho_b)
        vtau_a_fd = (e_a_plus - e_a_minus) / (delta_h * 2)
        self.check_term(tol, 'vtau_a', vtau_a, vtau_a_fd)

        v2tau2_aa = fxc["v2tau2"][0][0]
        vtau_a_plus = vxc_plus["vtau"][0][0]
        vtau_a_minus = vxc_minus["vtau"][0][0]
        v2tau2_aa_fd = (vtau_a_plus - vtau_a_minus) / (delta_h * 2)
        self.check_term(tol, "v2tau2_aa", v2tau2_aa, v2tau2_aa_fd)

        v2tau2_ab = fxc["v2tau2"][0][1]
        vtau_b_plus = vxc_plus["vtau"][0][1]
        vtau_b_minus = vxc_minus["vtau"][0][1]
        v2tau2_ab_fd = (vtau_b_plus - vtau_b_minus) / (delta_h * 2)
        self.check_term(tol, "v2tau2_ab", v2tau2_ab, v2tau2_ab_fd)

        v3tau3_aaa = kxc["v3tau3"][0][0]
        v2tau2_aa_plus = fxc_plus["v2tau2"][0][0]
        v2tau2_aa_minus = fxc_minus["v2tau2"][0][0]
        v3tau3_aaa_fd = (v2tau2_aa_plus - v2tau2_aa_minus) / (delta_h * 2)
        self.check_term(tol, "v3tau3_aaa", v3tau3_aaa, v3tau3_aaa_fd)

        v3tau3_aab = kxc["v3tau3"][0][1]
        v2tau2_ab_plus = fxc_plus["v2tau2"][0][1]
        v2tau2_ab_minus = fxc_minus["v2tau2"][0][1]
        v3tau3_aab_fd = (v2tau2_ab_plus - v2tau2_ab_minus) / (delta_h * 2)
        self.check_term(tol, "v3tau3_aab", v3tau3_aab, v3tau3_aab_fd)

        v3tau3_abb = kxc["v3tau3"][0][2]
        v2tau2_bb_plus = fxc_plus["v2tau2"][0][2]
        v2tau2_bb_minus = fxc_minus["v2tau2"][0][2]
        v3tau3_abb_fd = (v2tau2_bb_plus - v2tau2_bb_minus) / (delta_h * 2)
        self.check_term(tol, "v3tau3_abb", v3tau3_abb, v3tau3_abb_fd)

        v4tau4_aaaa = lxc["v4tau4"][0][0]
        v3tau3_aaa_plus = kxc_plus["v3tau3"][0][0]
        v3tau3_aaa_minus = kxc_minus["v3tau3"][0][0]
        v4tau4_aaaa_fd = (v3tau3_aaa_plus - v3tau3_aaa_minus) / (delta_h * 2)
        self.check_term(tol, "v4tau4_aaaa", v4tau4_aaaa, v4tau4_aaaa_fd)

        v4tau4_aaab = lxc["v4tau4"][0][1]
        v3tau3_aab_plus = kxc_plus["v3tau3"][0][1]
        v3tau3_aab_minus = kxc_minus["v3tau3"][0][1]
        v4tau4_aaab_fd = (v3tau3_aab_plus - v3tau3_aab_minus) / (delta_h * 2)
        self.check_term(tol, "v4tau4_aaab", v4tau4_aaab, v4tau4_aaab_fd)

        v4tau4_aabb = lxc["v4tau4"][0][2]
        v3tau3_abb_plus = kxc_plus["v3tau3"][0][2]
        v3tau3_abb_minus = kxc_minus["v3tau3"][0][2]
        v4tau4_aabb_fd = (v3tau3_abb_plus - v3tau3_abb_minus) / (delta_h * 2)
        self.check_term(tol, "v4tau4_aabb", v4tau4_aabb, v4tau4_aabb_fd)

        v4tau4_abbb = lxc["v4tau4"][0][3]
        v3tau3_bbb_plus = kxc_plus["v3tau3"][0][3]
        v3tau3_bbb_minus = kxc_minus["v3tau3"][0][3]
        v4tau4_abbb_fd = (v3tau3_bbb_plus - v3tau3_bbb_minus) / (delta_h * 2)
        self.check_term(tol, "v4tau4_abbb", v4tau4_abbb, v4tau4_abbb_fd)

    def test_tau_a_fd_pkzb(self):

        self.run_tau_a_fd('pkzb', 1.0e-7)

    def test_tau_a_fd_tpssh(self):

        self.run_tau_a_fd('tpssh', 1.0e-7)

    def test_tau_a_fd_m06(self):

        self.run_tau_a_fd('m06', 1.0e-7)

    def test_tau_a_fd_scan(self):

        self.run_tau_a_fd('scan', 1.0e-7)
