import numpy as np

from veloxchem.veloxchemlib import XCIntegrator


class TestLibxcRhoB:

    def check_term(self, tol, label, val, val_fd):

        if abs(val_fd) < tol:
            assert abs(val - val_fd) < tol
        else:
            assert abs(val - val_fd) / abs(val_fd) < tol

    def run_rho_b_fd(self, xcfun_label, tol):

        delta_h = 1.0e-6

        xc_drv = XCIntegrator()

        rho_a, rho_b = 0.17, 0.15

        inp = {
            'rho': np.array([rho_a, rho_b]),
        }

        lxc = xc_drv.compute_lxc_for_lda(xcfun_label, inp['rho'])
        kxc = xc_drv.compute_kxc_for_lda(xcfun_label, inp['rho'])
        fxc = xc_drv.compute_fxc_for_lda(xcfun_label, inp['rho'])
        vxc = xc_drv.compute_exc_vxc_for_lda(xcfun_label, inp['rho'])

        inp_plus = dict(inp)
        inp_plus['rho'][1] = rho_b + delta_h

        kxc_plus = xc_drv.compute_kxc_for_lda(xcfun_label, inp_plus['rho'])
        fxc_plus = xc_drv.compute_fxc_for_lda(xcfun_label, inp_plus['rho'])
        vxc_plus = xc_drv.compute_exc_vxc_for_lda(xcfun_label, inp_plus['rho'])
        exc_plus = {'zk': vxc_plus['exc']}

        inp_minus = dict(inp)
        inp_minus['rho'][1] = rho_b - delta_h

        kxc_minus = xc_drv.compute_kxc_for_lda(xcfun_label, inp_minus['rho'])
        fxc_minus = xc_drv.compute_fxc_for_lda(xcfun_label, inp_minus['rho'])
        vxc_minus = xc_drv.compute_exc_vxc_for_lda(xcfun_label,
                                                   inp_minus['rho'])
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

        v3rho3_bbb = kxc["v3rho3"][0][3]
        v2rho2_bb_plus = fxc_plus["v2rho2"][0][2]
        v2rho2_bb_minus = fxc_minus["v2rho2"][0][2]
        v3rho3_bbb_fd = (v2rho2_bb_plus - v2rho2_bb_minus) / (delta_h * 2)
        self.check_term(tol, "v3rho3_bbb", v3rho3_bbb, v3rho3_bbb_fd)

        v4rho4_bbbb = lxc["v4rho4"][0][4]
        v3rho3_bbb_plus = kxc_plus["v3rho3"][0][3]
        v3rho3_bbb_minus = kxc_minus["v3rho3"][0][3]
        v4rho4_bbbb_fd = (v3rho3_bbb_plus - v3rho3_bbb_minus) / (delta_h * 2)
        self.check_term(tol, "v4rho4_bbbb", v4rho4_bbbb, v4rho4_bbbb_fd)

    def test_rho_b_fd_slater(self):

        self.run_rho_b_fd('slater', 1.0e-7)

    def test_rho_b_fd_slda(self):

        self.run_rho_b_fd('slda', 1.0e-7)
