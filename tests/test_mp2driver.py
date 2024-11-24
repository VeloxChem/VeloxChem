from pathlib import Path
import pytest

from veloxchem.veloxchemlib import mpi_master
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.scfunrestdriver import ScfUnrestrictedDriver
from veloxchem.mp2driver import Mp2Driver
from veloxchem.mpitask import MpiTask


@pytest.mark.solvers
class TestMp2Driver:

    def test_h2se_mp2(self):

        # scf
        here = Path(__file__).parent
        inpfile = str(here / 'data' / 'h2se.inp')

        task = MpiTask([inpfile, None])

        scf_drv = ScfRestrictedDriver(task.mpi_comm, task.ostream)
        scf_results = scf_drv.compute(task.molecule, task.ao_basis,
                                      task.min_basis)

        # mp2
        e_ref = -0.28529088

        mp2_drv = Mp2Driver(task.mpi_comm, task.ostream)

        mp2_drv.conventional = True
        mp2_result = mp2_drv.compute(task.molecule, task.ao_basis, scf_results)
        if task.mpi_rank == mpi_master():
            assert abs(e_ref - mp2_result['mp2_energy']) < 1.0e-8

        mp2_drv.conventional = False
        mp2_result = mp2_drv.compute(task.molecule, task.ao_basis, scf_results)
        if task.mpi_rank == mpi_master():
            assert abs(e_ref - mp2_result['mp2_energy']) < 1.0e-8

        task.finish()

    def test_h2se_ump2(self):

        # scf
        here = Path(__file__).parent
        inpfile = str(here / 'data' / 'h2se.inp')

        task = MpiTask([inpfile, None])
        task.molecule.set_multiplicity(3)
        task.molecule.check_multiplicity()

        scf_drv = ScfUnrestrictedDriver(task.mpi_comm, task.ostream)
        scf_results = scf_drv.compute(task.molecule, task.ao_basis,
                                      task.min_basis)

        # mp2
        e_ref = -0.26775296

        mp2_drv = Mp2Driver(task.mpi_comm, task.ostream)

        mp2_drv.conventional = True
        mp2_result = mp2_drv.compute(task.molecule, task.ao_basis, scf_results)
        if task.mpi_rank == mpi_master():
            assert abs(e_ref - mp2_result['mp2_energy']) < 1.0e-7

        mp2_drv.conventional = False
        mp2_result = mp2_drv.compute(task.molecule, task.ao_basis, scf_results)
        if task.mpi_rank == mpi_master():
            assert abs(e_ref - mp2_result['mp2_energy']) < 1.0e-7

    def test_mp2_update_settings(self):

        mp2_dict = {
            'eri_thresh': 1e-13,
        }

        mp2_drv = Mp2Driver()

        for key, val in mp2_dict.items():
            assert getattr(mp2_drv, key) != val

        mp2_drv.update_settings(mp2_dict)

        for key, val in mp2_dict.items():
            assert getattr(mp2_drv, key) == val
