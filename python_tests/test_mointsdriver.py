from pathlib import Path
import numpy as np

from veloxchem.veloxchemlib import is_mpi_master
from veloxchem.veloxchemlib import moints
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.mointsdriver import MOIntegralsDriver
from veloxchem.mp2driver import Mp2Driver
from veloxchem.mpitask import MpiTask


class TestMOIntegralsDriver:

    def test_moints_type(self):

        moints_drv = MOIntegralsDriver()

        assert moints_drv.get_moints_type("OOOO") == moints.oooo
        assert moints_drv.get_moints_type("OOOV") == moints.ooov
        assert moints_drv.get_moints_type("OVOV") == moints.ovov
        assert moints_drv.get_moints_type("OOVV") == moints.oovv
        assert moints_drv.get_moints_type("OVVV") == moints.ovvv
        assert moints_drv.get_moints_type("VVVV") == moints.vvvv

    def test_h2se_mp2(self):

        # scf
        here = Path(__file__).parent
        inpfile = str(here / 'inputs' / 'h2se.inp')

        task = MpiTask([inpfile, None])

        scf_drv = ScfRestrictedDriver(task.mpi_comm, task.ostream)
        scf_drv.compute(task.molecule, task.ao_basis, task.min_basis)
        mol_orbs = scf_drv.mol_orbs

        # mp2
        mp2_drv = Mp2Driver(task.mpi_comm, task.ostream)
        mp2_drv.compute(task.molecule, task.ao_basis, mol_orbs)

        if is_mpi_master(task.mpi_comm):
            e_ref = -0.28529088
            e_mp2 = mp2_drv.e_mp2
            assert abs(e_ref - e_mp2) < 1.0e-8

        mp2_drv.update_settings({'conventional': 'yes'})
        mp2_drv.compute(task.molecule, task.ao_basis, mol_orbs)

        if is_mpi_master(task.mpi_comm):
            assert abs(e_ref - mp2_drv.e_mp2) < 1.0e-8

        # extra test: collect moints batches to master node
        moints_drv = MOIntegralsDriver(task.mpi_comm, task.ostream)
        grps = [p for p in range(task.mpi_size)]
        oovv = moints_drv.compute(task.molecule, task.ao_basis, mol_orbs,
                                  'OOVV', grps)
        moints_drv.collect_moints_batches(oovv, grps)

        if is_mpi_master(task.mpi_comm):
            orb_ene = mol_orbs.ea_to_numpy()
            nocc = task.molecule.number_of_alpha_electrons()
            eocc = orb_ene[:nocc]
            evir = orb_ene[nocc:]
            e_vv = evir.reshape(-1, 1) + evir

            # loop over generator pairs for mp2 energy
            master_e_mp2 = 0.0
            for pair in oovv.get_gen_pairs():
                ab = oovv.to_numpy(pair)
                denom = e_vv - eocc[pair.first()] - eocc[pair.second()]
                master_e_mp2 -= np.sum(ab * (2.0 * ab - ab.T) / denom)

            assert abs(e_ref - master_e_mp2) < 1.0e-8

        # in-mem test (will be deprecated)
        phys_oovv = moints_drv.compute_in_mem(task.molecule, task.ao_basis,
                                              mol_orbs, 'OOVV')

        if is_mpi_master(task.mpi_comm):
            in_mem_e_mp2 = 0.0
            for i in range(phys_oovv.shape[0]):
                for j in range(phys_oovv.shape[1]):
                    ab = phys_oovv[i, j, :, :]
                    denom = e_vv - eocc[i] - eocc[j]
                    in_mem_e_mp2 -= np.sum(ab * (2.0 * ab - ab.T) / denom)

            assert abs(e_ref - in_mem_e_mp2) < 1.0e-8

        task.finish()

    def test_mp2_update_settings(self):

        mp2_dict = {
            'qq_type': 'QQ',
            'eri_thresh': 1e-13,
            'batch_size': 99,
            'conventional': True,
        }

        mp2_drv = Mp2Driver()

        for key, val in mp2_dict.items():
            assert getattr(mp2_drv, key) != val

        mp2_drv.update_settings(mp2_dict)

        for key, val in mp2_dict.items():
            assert getattr(mp2_drv, key) == val
