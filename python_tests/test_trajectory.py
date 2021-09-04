import sys
from pathlib import Path

import pytest

try:
    import MDAnalysis
except ImportError:
    pass

try:
    import cppe
except ImportError:
    pass

from veloxchem.mpitask import MpiTask
from veloxchem.trajectorydriver import TrajectoryDriver
from veloxchem.veloxchemlib import is_mpi_master


@pytest.mark.skipif('cppe' not in sys.modules, reason='cppe not available')
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

        if is_mpi_master(task.mpi_comm):
            exc_energies = []
            osc_strengths = []

            dir_path = Path(filename + '_files')
            files = sorted([x for x in dir_path.iterdir() if x.is_file()])
            for x in files:
                if x.name.startswith(filename) and (x.suffix == '.out'):
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
                assert abs(e - ref_e) < 1.0e-6
                assert abs(f - ref_f) < 1.0e-4

    @pytest.mark.skipif('MDAnalysis' not in sys.modules,
                        reason='MDAnalysis not available')
    def test_trajectory_bithio(self):

        here = Path(__file__).parent
        filename = str(here / 'inputs' / 'bithio')

        ref_exc_energies = [5.790027, 7.172051, 5.628371, 7.199163]
        ref_osc_strengths = [0.5414, 0.0181, 0.5529, 0.1210]

        self.run_trajectory(filename, ref_exc_energies, ref_osc_strengths)
