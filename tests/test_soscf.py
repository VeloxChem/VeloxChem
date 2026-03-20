import pytest

from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.soscf import Soscf


@pytest.mark.solvers
class TestSoscfSwitching:

    def test_should_activate_soscf_with_persistent_stall(self):

        scf_drv = ScfRestrictedDriver()
        scf_drv.acc_type = 'DIIS_SOSCF'
        scf_drv.soscf_switch_thresh = 1.0e-2
        scf_drv.soscf_switch_ratio = 0.7
        scf_drv.soscf_switch_window = 4
        scf_drv.soscf_switch_persistence = 2

        history = [
            {
                'energy': -2335.166073343963,
                'diff_energy': -0.0000978576,
                'gradient_norm': 0.00281164,
                'max_gradient': 0.00016019,
                'diff_density': 0.04845477,
            },
            {
                'energy': -2335.166097637520,
                'diff_energy': -0.0000242936,
                'gradient_norm': 0.00225906,
                'max_gradient': 0.00013435,
                'diff_density': 0.01142867,
            },
            {
                'energy': -2335.166153195724,
                'diff_energy': -0.0000555582,
                'gradient_norm': 0.00253967,
                'max_gradient': 0.00012442,
                'diff_density': 0.02814867,
            },
            {
                'energy': -2335.166218287050,
                'diff_energy': -0.0000650913,
                'gradient_norm': 0.00361708,
                'max_gradient': 0.00012870,
                'diff_density': 0.03515947,
            },
            {
                'energy': -2335.166263580340,
                'diff_energy': -0.0000452933,
                'gradient_norm': 0.00280575,
                'max_gradient': 0.00011244,
                'diff_density': 0.02633746,
            },
            {
                'energy': -2335.166347099875,
                'diff_energy': -0.0000835195,
                'gradient_norm': 0.00175920,
                'max_gradient': 0.00012595,
                'diff_density': 0.05348030,
            },
        ]

        scf_drv._history = history[:4]
        assert scf_drv._should_activate_soscf() is False
        assert scf_drv._soscf_switch_counter == 1

        scf_drv._history = history[:5]
        assert scf_drv._should_activate_soscf() is True
        assert scf_drv._soscf_switch_counter == 2

    def test_lbfgs_history_is_bounded_and_updates(self):

        soscf = Soscf()
        options = {
            'use_bfgs': True,
            'max_step_norm': 10.0,
            'damping': 1.0,
            'min_damping': 0.1,
            'damping_decay': 0.5,
            'damping_growth': 1.2,
            'reset_ratio': 10.0,
            'improve_ratio': 0.8,
            'curvature_tol': 1.0e-12,
            'diag_blend': 0.0,
            'history_size': 2,
            'reset_persistence': 2,
            'success_recovery_steps': 3,
        }

        diag_inv = [1.0, 0.5]

        step_1, info_1 = soscf.compute_step([1.0, -0.5], diag_inv, options)
        assert info_1['reset'] is True
        assert len(soscf.history) == 0

        step_2, info_2 = soscf.compute_step([0.6, -0.2], diag_inv, options)
        assert info_2['bfgs_updated'] is True
        assert len(soscf.history) == 1

        step_3, info_3 = soscf.compute_step([0.3, -0.1], diag_inv, options)
        assert info_3['bfgs_updated'] is True
        assert len(soscf.history) == 2

        step_4, info_4 = soscf.compute_step([0.1, -0.05], diag_inv, options)
        assert info_4['bfgs_updated'] is True
        assert len(soscf.history) == 2

        assert step_1.shape == (2,)
        assert step_2.shape == (2,)
        assert step_3.shape == (2,)
        assert step_4.shape == (2,)

    def test_single_bad_step_backs_off_without_reset(self):

        soscf = Soscf()
        options = {
            'use_bfgs': True,
            'max_step_norm': 10.0,
            'damping': 1.0,
            'min_damping': 0.1,
            'damping_decay': 0.5,
            'damping_growth': 1.2,
            'reset_ratio': 1.25,
            'improve_ratio': 0.8,
            'curvature_tol': 1.0e-12,
            'diag_blend': 0.0,
            'history_size': 4,
            'reset_persistence': 2,
            'success_recovery_steps': 3,
        }

        diag_inv = [1.0, 1.0]

        soscf.compute_step([1.0, 0.0], diag_inv, options)
        soscf.compute_step([0.6, 0.0], diag_inv, options)
        assert len(soscf.history) == 1

        _, info = soscf.compute_step([1.0, 0.0], diag_inv, options)
        assert info['backoff'] is True
        assert info['reset'] is False
        assert len(soscf.history) == 1
        assert soscf.damping == pytest.approx(0.5)

    def test_repeated_bad_steps_trigger_reset(self):

        soscf = Soscf()
        options = {
            'use_bfgs': True,
            'max_step_norm': 10.0,
            'damping': 1.0,
            'min_damping': 0.1,
            'damping_decay': 0.5,
            'damping_growth': 1.2,
            'reset_ratio': 1.25,
            'improve_ratio': 0.8,
            'curvature_tol': 1.0e-12,
            'diag_blend': 0.0,
            'history_size': 4,
            'reset_persistence': 2,
            'success_recovery_steps': 3,
        }

        diag_inv = [1.0, 1.0]

        soscf.compute_step([1.0, 0.0], diag_inv, options)
        soscf.compute_step([0.6, 0.0], diag_inv, options)

        soscf.compute_step([1.0, 0.0], diag_inv, options)
        _, info = soscf.compute_step([1.3, 0.0], diag_inv, options)

        assert info['backoff'] is True
        assert info['reset'] is True
        assert len(soscf.history) == 0

    def test_successful_steps_recover_damping(self):

        soscf = Soscf()
        options = {
            'use_bfgs': True,
            'max_step_norm': 10.0,
            'damping': 1.0,
            'min_damping': 0.1,
            'damping_decay': 0.5,
            'damping_growth': 1.2,
            'reset_ratio': 1.25,
            'improve_ratio': 0.95,
            'curvature_tol': 1.0e-12,
            'diag_blend': 0.0,
            'history_size': 4,
            'reset_persistence': 2,
            'success_recovery_steps': 2,
        }

        diag_inv = [1.0, 1.0]

        soscf.compute_step([1.0, 0.0], diag_inv, options)
        soscf.compute_step([0.6, 0.0], diag_inv, options)
        soscf.compute_step([1.0, 0.0], diag_inv, options)
        assert soscf.damping == pytest.approx(0.5)

        soscf.compute_step([0.8, 0.0], diag_inv, options)
        soscf.compute_step([0.6, 0.0], diag_inv, options)

        assert soscf.damping > 0.5
