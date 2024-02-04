from pathlib import Path
import numpy as np
import pytest
import sys

try:
    import openmm
except ImportError:
    pass

from veloxchem.veloxchemlib import is_mpi_master
from veloxchem.mpitask import MpiTask
from veloxchem.forcefieldgenerator import ForceFieldGenerator


@pytest.mark.filterwarnings(
    'ignore:.*tostring.*tobytes:DeprecationWarning:geometric')
class TestForceField:

    @pytest.mark.skipif('openmm' not in sys.modules,
                        reason='openmm not available')
    def test_force_field(self):

        # vlxtag: RKS, Force_Field_Generation

        here = Path(__file__).parent
        inpfile = str(here / 'inputs' / 'butane.inp')

        task = MpiTask([inpfile, None])

        ff_dict = (task.input_dict['force_field']
                   if 'force_field' in task.input_dict else {})
        resp_dict = (task.input_dict['resp_charges']
                     if 'resp_charges' in task.input_dict else {})

        ff_dict['filename'] = task.input_dict['filename']
        resp_dict['filename'] = task.input_dict['filename']

        ff_drv = ForceFieldGenerator(task.mpi_comm, task.ostream)
        ff_drv.update_settings(ff_dict, resp_dict)
        ff_drv.ostream.mute()
        ff_drv.compute(task.molecule, task.ao_basis)

        if is_mpi_master(task.mpi_comm):
            ref_atom_inds = [1, 5, 8, 11]
            ref_funct_type = 3
            ref_parameters = np.array(
                [4.82579, 10.74424, 0.17548, -15.10881, 1.85225, -2.59463])

            dih_parameters = None
            new_itp_file = here / 'inputs' / 'butane_files' / 'butane_01.itp'

            with new_itp_file.open('r') as f_itp:
                for line in f_itp:
                    content = line.split()

                    try:
                        dih_atom_inds = [int(x) for x in content[:4]]
                        dih_funct_type = int(content[4])
                        if (dih_atom_inds == ref_atom_inds and
                                dih_funct_type == ref_funct_type):
                            dih_parameters = np.array(
                                [float(x) for x in content[5:11]])
                    except (ValueError, IndexError):
                        continue

            assert dih_parameters is not None
            assert np.max(np.abs(dih_parameters - ref_parameters)) < 1.0e-3

            pdb_file = Path(inpfile).with_suffix('.pdb')
            if pdb_file.is_file():
                pdb_file.unlink()

            scf_h5_file = Path(inpfile).with_suffix('.scf.h5')
            if scf_h5_file.is_file():
                scf_h5_file.unlink()

            scf_final_h5_file = Path(inpfile).with_suffix('.scf.tensors.h5')
            if scf_final_h5_file.is_file():
                scf_final_h5_file.unlink()

            mol_name = Path(inpfile).stem
            ff_dir = Path(inpfile).parent / (mol_name + '_files')
            for ffver in [0, 1]:
                for ftype in ['top', 'itp']:
                    ff_file = ff_dir / (mol_name + f'_{ffver:02d}.{ftype}')
                    if ff_file.is_file():
                        ff_file.unlink()
