from pathlib import Path
import pytest

from veloxchem.veloxchemlib import mpi_master
from veloxchem.mpitask import MpiTask
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.rspshg import SHG


@pytest.mark.solvers
class TestSHG:

    def run_scf(self, task, *, ri):

        scf_drv = ScfRestrictedDriver(task.mpi_comm, task.ostream)
        scf_drv.ri = ri
        scf_drv.update_settings(task.input_dict['scf'],
                                task.input_dict['method_settings'])
        scf_results = scf_drv.compute(task.molecule, task.ao_basis)

        if scf_drv.rank == mpi_master():
            if ri:
                assert scf_results['ri_jk'] or scf_results['ri_coulomb']
                assert not scf_results['ri_fallback']

        return scf_results

    def run_shg(self, inpfile, xcfun_label, shg_type, ref_result, *, ri):

        task = MpiTask([inpfile, None])
        task.input_dict['scf']['checkpoint_file'] = None

        if xcfun_label is not None:
            task.input_dict['method_settings']['xcfun'] = xcfun_label

        if shg_type is not None:
            task.input_dict['response']['shg_type'] = shg_type

        scf_results = self.run_scf(task, ri=ri)

        shg_prop = SHG(task.input_dict['response'],
                       task.input_dict['method_settings'])
        shg_prop.init_driver(task.mpi_comm, task.ostream)
        shg_prop.compute(task.molecule, task.ao_basis, scf_results)
        shg_result = shg_prop.rsp_property

        if task.mpi_rank == mpi_master():
            freq = list(shg_result['beta'].keys())[0]

            for ind, comp in enumerate('xyz'):
                assert abs(shg_result['beta'][freq][ind].real -
                           ref_result[comp].real) < 1.0e-5
                assert abs(shg_result['beta'][freq][ind].imag -
                           ref_result[comp].imag) < 1.0e-5

            assert abs(shg_result['beta_bar'][freq].real -
                       ref_result['beta_bar'].real) < 1.0e-5
            assert abs(shg_result['beta_bar'][freq].imag -
                       ref_result['beta_bar'].imag) < 1.0e-5

    @pytest.mark.parametrize('ri', [False, True])
    def test_shg_hf(self, ri):

        # vlxtag: RHF, SHG, QR

        here = Path(__file__).parent
        inpfile = str(here / 'data' / 'methanol_shg.inp')

        xcfun_label = None
        shg_type = 'full'

        ref_result = {
            False: {
                'x': -40.52766331 - 11.31636048j,
                'y': -0.11777265 - 0.02104016j,
                'z': 24.95268906 + 7.53241250j,
                'beta_bar': -39.79479881 - 11.62336253j,
            },
            True: {
                'x': -40.52652907 - 11.31660904j,
                'y': -0.11775838 - 0.02103991j,
                'z': 24.95256399 + 7.53263538j,
                'beta_bar': -39.79471569 - 11.62380796j
            }
        }[ri]

        self.run_shg(inpfile, xcfun_label, shg_type, ref_result, ri=ri)

    @pytest.mark.parametrize('ri', [False, True])
    def test_shg_b3lyp(self, ri):

        # vlxtag: RKS, SHG, QR

        here = Path(__file__).parent
        inpfile = str(here / 'data' / 'methanol_shg.inp')

        xcfun_label = 'b3lyp'
        shg_type = 'full'

        ref_result = {
            False: {
                'x': -48.55320751 - 20.31312232j,
                'y': -0.09597950 - 0.00334644j,
                'z': 32.09728459 + 15.20019960j,
                'beta_bar': -47.13862975 - 21.37978710j,
            },
            True: {
                'x': -48.56128288 - 20.32003954j,
                'y': -0.09608696 - 0.00340773j,
                'z': 32.09824726 + 15.20248473j,
                'beta_bar': -47.14396310 - 21.38499733j
            }
        }[ri]

        self.run_shg(inpfile, xcfun_label, shg_type, ref_result, ri=ri)

    @pytest.mark.parametrize('ri', [False, True])
    def test_shg_reduced_hf(self, ri):

        # vlxtag: RHF, SHG, QR

        here = Path(__file__).parent
        inpfile = str(here / 'data' / 'methanol_shg.inp')

        xcfun_label = None
        shg_type = 'reduced'

        ref_result = {
            False: {
                'x': -40.53256501 - 11.43101669j,
                'y': -0.11772170 - 0.02074979j,
                'z': 24.95878533 + 7.63353776j,
                'beta_bar': -39.80239953 - 11.76359430j,
            },
            True: {
                'x': -40.53139865 - 11.43131522j,
                'y': -0.11770724 - 0.02074993j,
                'z': 24.95864531 + 7.63377984j,
                'beta_bar': -39.80228929 - 11.76407952j
            }
        }[ri]

        self.run_shg(inpfile, xcfun_label, shg_type, ref_result, ri=ri)

    @pytest.mark.parametrize('ri', [False, True])
    def test_shg_reduced_b3lyp(self, ri):

        # vlxtag: RKS, SHG, QR

        here = Path(__file__).parent
        inpfile = str(here / 'data' / 'methanol_shg.inp')

        xcfun_label = 'b3lyp'
        shg_type = 'reduced'

        ref_result = {
            False: {
                'x': -48.45584560 - 20.92686505j,
                'y': -0.09586086 - 0.00356429j,
                'z': 32.02926487 + 15.65493314j,
                'beta_bar': -47.04068190 - 22.02152032j,
            },
            True: {
                'x': -48.46382386 - 20.93407136j,
                'y': -0.09596771 - 0.00362744j,
                'z': 32.03018306 + 15.65734724j,
                'beta_bar': -47.04593630 - 22.02697234j
            }
        }[ri]

        self.run_shg(inpfile, xcfun_label, shg_type, ref_result, ri=ri)
