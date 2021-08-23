from pathlib import Path
from unittest.mock import patch
import sys

from veloxchem.veloxchemlib import is_mpi_master
from veloxchem.mpitask import MpiTask
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.scfunrestdriver import ScfUnrestrictedDriver
from veloxchem.rspabsorption import Absorption
from veloxchem.rspc6 import C6
from veloxchem.rspcdspec import CircularDichroismSpectrum
from veloxchem.rsplinabscross import LinearAbsorptionCrossSection
from veloxchem.rsppolarizability import Polarizability
from veloxchem.rsptpa import TPA
from veloxchem.main import select_scf_driver, select_rsp_property, main


class TestMain:

    def test_select_scf_driver(self):

        # restricted

        here = Path(__file__).parent
        inpfile = here / 'inputs' / 'water.inp'

        task = MpiTask([str(inpfile), None])
        scf_type = 'restricted'

        scf_drv = select_scf_driver(task, scf_type)
        assert isinstance(scf_drv, ScfRestrictedDriver)

        # unrestricted

        here = Path(__file__).parent
        inpfile = here / 'inputs' / 'water_unrest.inp'

        task = MpiTask([str(inpfile), None])
        scf_type = 'unrestricted'

        scf_drv = select_scf_driver(task, scf_type)
        assert isinstance(scf_drv, ScfUnrestrictedDriver)

    def test_select_rsp_driver(self):

        here = Path(__file__).parent
        inpfile = here / 'inputs' / 'water.inp'

        task = MpiTask([str(inpfile), None])
        task.input_dict['scf']['checkpoint_file'] = None

        scf_drv = ScfRestrictedDriver(task.mpi_comm, task.ostream)
        scf_drv.update_settings(task.input_dict['scf'],
                                task.input_dict['method_settings'])
        scf_drv.compute(task.molecule, task.ao_basis, task.min_basis)

        method_dict = {}

        rsp_dict = {'property': 'polarizability'}
        rsp_prop = select_rsp_property(task, scf_drv.mol_orbs, rsp_dict,
                                       method_dict)
        assert isinstance(rsp_prop, Polarizability)

        rsp_dict = {'property': 'absorption'}
        rsp_prop = select_rsp_property(task, scf_drv.mol_orbs, rsp_dict,
                                       method_dict)
        assert isinstance(rsp_prop, Absorption)

        rsp_dict = {'property': 'absorption (cpp)'}
        rsp_prop = select_rsp_property(task, scf_drv.mol_orbs, rsp_dict,
                                       method_dict)
        assert isinstance(rsp_prop, LinearAbsorptionCrossSection)

        rsp_dict = {'property': 'ecd (cpp)'}
        rsp_prop = select_rsp_property(task, scf_drv.mol_orbs, rsp_dict,
                                       method_dict)
        assert isinstance(rsp_prop, CircularDichroismSpectrum)

        rsp_dict = {'property': 'c6'}
        rsp_prop = select_rsp_property(task, scf_drv.mol_orbs, rsp_dict,
                                       method_dict)
        assert isinstance(rsp_prop, C6)

        rsp_dict = {'property': 'tpa'}
        rsp_prop = select_rsp_property(task, scf_drv.mol_orbs, rsp_dict,
                                       method_dict)
        assert isinstance(rsp_prop, TPA)

    def test_main(self, capsys):

        here = Path(__file__).parent
        inpfile = here / 'inputs' / 'water_maxh.inp'

        with patch.object(sys, 'argv', ['vlx', str(inpfile)]):
            main()
            captured = capsys.readouterr()

            if is_mpi_master():
                assert 'VeloxChem execution started' in captured.out
                assert 'SCF converged' in captured.out
                assert 'Total Energy                       :' in captured.out
                assert 'Electronic Energy                  :' in captured.out
                assert 'Nuclear Repulsion Energy           :' in captured.out
                assert 'VeloxChem execution completed' in captured.out

                scf_h5 = inpfile.with_suffix('.scf.h5')
                if scf_h5.is_file():
                    scf_h5.unlink()
                scf_final_h5 = scf_h5.with_suffix('.tensors.h5')
                if scf_final_h5.is_file():
                    scf_final_h5.unlink()
