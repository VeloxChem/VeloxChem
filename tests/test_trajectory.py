from pathlib import Path
import pytest
import sys

from veloxchem.mpitask import MpiTask
from veloxchem.trajectorydriver import TrajectoryDriver
from veloxchem.veloxchemlib import mpi_master

try:
    import pyframe
except ImportError:
    pass


@pytest.mark.solvers
class TestTrajectoryDriver:

    def run_trajectory(self, filename, ref_exc_energies, ref_osc_strengths):

        inpfile = f'{filename}.inp'
        xtcfile = f'{filename}.xtc'
        tprfile = f'{filename}.tpr'

        task = MpiTask([inpfile, None])

        task.input_dict['trajectory']['trajectory_file'] = xtcfile
        task.input_dict['trajectory']['topology_file'] = tprfile
        task.input_dict['trajectory']['filename'] = filename
        task.input_dict['trajectory']['charges'] = task.input_dict['charges']
        task.input_dict['trajectory']['polarizabilities'] = task.input_dict[
            'polarizabilities']

        traj_drv = TrajectoryDriver(task.mpi_comm, task.ostream)
        traj_drv.update_settings(task.input_dict['trajectory'],
                                 task.input_dict['spectrum_settings'],
                                 task.input_dict['response'],
                                 task.input_dict['method_settings'])
        traj_drv.compute(task.molecule, task.ao_basis, task.min_basis)

        if task.mpi_rank == mpi_master():
            exc_energies = []
            osc_strengths = []

            dir_path = Path(filename + '_files')
            files = sorted([x for x in dir_path.iterdir() if x.is_file()])
            for x in files:
                if (x.name.startswith(Path(filename).name) and
                        x.suffix == '.out'):
                    with x.open() as fh:
                        for line in fh:
                            if 'Osc.Str.' in line:
                                content = line.split()
                                exc_energies.append(float(content[5]))
                                osc_strengths.append(float(content[8]))
                x.unlink()

            npzfile = dir_path.parent / '.bithio.xtc_offsets.npz'
            if npzfile.is_file():
                npzfile.unlink()

            for (e, f, ref_e, ref_f) in zip(exc_energies, osc_strengths,
                                            ref_exc_energies,
                                            ref_osc_strengths):
                assert abs(e - ref_e) < 0.05
                assert abs(f - ref_f) < 0.01

    @pytest.mark.skipif('pyframe' not in sys.modules,
                        reason='pyframe not available')
    def test_trajectory_bithio(self):

        here = Path(__file__).parent
        filename = str(here / 'data' / 'bithio')

        ref_exc_energies = [5.790027, 7.172051, 5.628371, 7.199163]
        ref_osc_strengths = [0.5414, 0.0181, 0.5529, 0.1210]

        self.run_trajectory(filename, ref_exc_energies, ref_osc_strengths)
