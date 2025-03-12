from pathlib import Path
import pytest
import sys

from veloxchem.veloxchemlib import mpi_master
from veloxchem.mpitask import MpiTask
from veloxchem.mmforcefieldgenerator import MMForceFieldGenerator

try:
    import scipy
except ImportError:
    pass


class TestForceField:

    @pytest.mark.skipif('scipy' not in sys.modules,
                        reason='scipy not available')
    def test_force_field(self):

        # vlxtag: RKS, MM_Force_Field_Generation

        here = Path(__file__).parent
        inpfile = str(here / 'data' / 'butane.inp')

        task = MpiTask([inpfile, None])

        ff_dict = (task.input_dict['force_field']
                   if 'force_field' in task.input_dict else {})
        resp_dict = (task.input_dict['resp_charges']
                     if 'resp_charges' in task.input_dict else {})

        ff_dict['filename'] = task.input_dict['filename']
        resp_dict['filename'] = task.input_dict['filename']

        ff_gen = MMForceFieldGenerator(task.mpi_comm, task.ostream)
        ff_gen.update_settings(ff_dict, resp_dict)
        ff_gen.compute(task.molecule, task.ao_basis)

        if ff_gen.rank == mpi_master():

            assert ff_gen.fitting_summary['maximum_difference'] < 3.0
            assert ff_gen.fitting_summary['standard_deviation'] < 1.2

            scf_h5_file = Path(inpfile).with_suffix('.scf.h5')
            if scf_h5_file.is_file():
                scf_h5_file.unlink()

            scf_final_h5_file = Path(inpfile).with_suffix('.scf.results.h5')
            if scf_final_h5_file.is_file():
                scf_final_h5_file.unlink()

            mol_name = Path(inpfile).stem
            ff_dir = Path(inpfile).parent / (mol_name + '_files')
            for ftype in ['top', 'itp']:
                ff_file = ff_dir / (mol_name + f'.{ftype}')
                if ff_file.is_file():
                    ff_file.unlink()
