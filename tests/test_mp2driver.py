from pathlib import Path
import pytest
from mpi4py import MPI

from veloxchem.veloxchemlib import mpi_master
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.scfunrestdriver import ScfUnrestrictedDriver
from veloxchem.mp2driver import Mp2Driver
from veloxchem.mpitask import MpiTask
from veloxchem.molecularorbitals import MolecularOrbitals, molorb


class TestMp2Driver:

    def run_mp2(self, task, scf_type, e_ref, mp2_type, tol):

        if scf_type.lower() == 'restricted':
            scf_drv = ScfRestrictedDriver(task.mpi_comm, task.ostream)
        else:
            scf_drv = ScfUnrestrictedDriver(task.mpi_comm, task.ostream)
        scf_results = scf_drv.compute(task.molecule, task.ao_basis)

        mp2_drv = Mp2Driver(task.mpi_comm, task.ostream)
        if mp2_type.lower() == 'conventional':
            mp2_drv.conventional = True
        else:
            mp2_drv.conventional = False
        mp2_result = mp2_drv.compute(task.molecule, task.ao_basis, scf_results)

        if task.mpi_rank == mpi_master():
            assert abs(e_ref - mp2_result['mp2_energy']) < tol

    @pytest.mark.solvers
    def test_mp2_conventional(self):

        here = Path(__file__).parent
        inpfile = str(here / 'data' / 'h2se.inp')

        task = MpiTask([inpfile, None])

        e_ref = -0.28529088

        self.run_mp2(task, 'restricted', e_ref, 'conventional', 1.0e-8)

    @pytest.mark.solvers
    def test_mp2_fockdriven(self):

        here = Path(__file__).parent
        inpfile = str(here / 'data' / 'h2se.inp')

        task = MpiTask([inpfile, None])

        e_ref = -0.28529088

        self.run_mp2(task, 'restricted', e_ref, 'fockdriven', 1.0e-8)

    @pytest.mark.solvers
    def test_ump2_conventional(self):

        here = Path(__file__).parent
        inpfile = str(here / 'data' / 'h2se.inp')

        task = MpiTask([inpfile, None])
        task.molecule.set_multiplicity(3)

        e_ref = -0.26775296

        self.run_mp2(task, 'unrestricted', e_ref, 'conventional', 1.0e-7)

    @pytest.mark.timeconsuming
    def test_ump2_fockdriven(self):

        here = Path(__file__).parent
        inpfile = str(here / 'data' / 'h2se.inp')

        task = MpiTask([inpfile, None])
        task.molecule.set_multiplicity(3)

        e_ref = -0.26775296

        self.run_mp2(task, 'unrestricted', e_ref, 'fockdriven', 1.0e-7)

    @pytest.mark.solvers
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

    @pytest.mark.solvers
    @pytest.mark.skipif(
        MPI.COMM_WORLD.Get_size() > 1,
        reason='pytest.raises assertions require a single MPI rank')
    def test_mp2_rejects_non_aufbau_occupations(self):

        here = Path(__file__).parent
        inpfile = str(here / 'data' / 'h2se.inp')

        task = MpiTask([inpfile, None])
        scf_drv = ScfRestrictedDriver(task.mpi_comm, task.ostream)
        scf_results = scf_drv.compute(task.molecule, task.ao_basis)

        occ_alpha = task.molecule.get_aufbau_alpha_occupation(
            scf_results['E_alpha'].shape[0], task.ao_basis)
        occ_alpha[task.molecule.number_of_alpha_occupied_orbitals(
            task.ao_basis) - 1] = 0.0
        occ_alpha[task.molecule.number_of_alpha_occupied_orbitals(
            task.ao_basis)] = 1.0

        mol_orbs = MolecularOrbitals([scf_results['C_alpha']],
                                     [scf_results['E_alpha']], [occ_alpha],
                                     molorb.rest)

        mp2_drv = Mp2Driver(task.mpi_comm, task.ostream)

        with pytest.raises(
                AssertionError,
                match='Alpha occupation numbers must follow Aufbau occupations'):
            mp2_drv.compute(task.molecule, task.ao_basis, mol_orbs)
