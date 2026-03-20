import pytest

from veloxchem.scfrestdriver import ScfRestrictedDriver


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
