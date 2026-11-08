import numpy as np
import pytest

from veloxchem.outputstream import OutputStream
from veloxchem.solvationfepdriver import SolvationFepDriver

mm = pytest.importorskip("openmm")
pytest.importorskip("pymbar")


class TestSolvationFepOptions:

    def driver(self):

        return SolvationFepDriver(ostream=OutputStream(None))

    def test_defaults_keep_the_four_stage_protocol(self):

        driver = self.driver()
        assert driver.skip_vacuum is False
        assert driver.coulomb_global_parameters == []

    def test_coulomb_global_parameters_follow_the_charge_stages(self):

        system = mm.System()
        force = mm.CustomNonbondedForce('-(1-c4_lambda)*C41*C42/r^4')
        force.addGlobalParameter('c4_lambda', 0.0)
        force.addGlobalParameter('other', 0.25)
        force.addPerParticleParameter('C4')
        system.addForce(force)
        system.addForce(mm.HarmonicBondForce())

        driver = self.driver()
        driver.coulomb_global_parameters = ['c4_lambda']
        for stage, expected in [(1, 0.4), (2, 1.0), (3, 0.4), (4, 1.0)]:
            driver.stage = stage
            driver._scale_coulomb_global_parameters(system, 0.4)
            assert force.getGlobalParameterDefaultValue(0) == expected
            assert force.getGlobalParameterDefaultValue(1) == 0.25

    @pytest.mark.parametrize('skip_vacuum, stages, total',
                             [(False, [1, 2, 3, 4], 1.0 + 2.0 - (3.0 + 4.0)),
                              (True, [1, 2], 1.0 + 2.0)])
    def test_skip_vacuum_runs_only_the_solvated_stages(self, skip_vacuum,
                                                       stages, total):

        driver = self.driver()
        driver.skip_vacuum = skip_vacuum
        run = []

        def fake_stage(stage, vacuum=False):
            run.append(stage)
            delta_f = {
                'Delta_f': np.zeros((2, 2)),
                'dDelta_f': np.zeros((2, 2))
            }
            return delta_f, float(stage)

        driver._run_lambda_simulations = fake_stage
        result = driver._run_stages()
        assert run == stages
        assert result['free_energy'] == pytest.approx(total)
        assert result['Stage 3']['Delta_f'] == 0.0
        assert result['Stage 4']['Delta_f'] == 0.0
