from mpi4py import MPI
import numpy as np
import unittest
import os

from veloxchem.veloxchemlib import mpi_master
from veloxchem.veloxchemlib import hartree_in_ev
from veloxchem.mpitask import MpiTask
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.tdaexcidriver import TDAExciDriver


class TestPolEmbed(unittest.TestCase):

    def test_pe_dft(self):

        try:
            import cppe
        except ImportError:
            return

        inpfile = os.path.join('inputs', 'pe_water.inp')
        if not os.path.isfile(inpfile):
            inpfile = os.path.join('python_tests', inpfile)

        potfile = os.path.join('inputs', 'pe_water.pot')
        if not os.path.isfile(potfile):
            potfile = os.path.join('python_tests', potfile)

        task = MpiTask([inpfile, None], MPI.COMM_WORLD)
        task.input_dict['method_settings']['potfile'] = potfile

        scf_drv = ScfRestrictedDriver(task.mpi_comm, task.ostream)
        scf_drv.update_settings(task.input_dict['scf'],
                                task.input_dict['method_settings'])
        scf_drv.compute(task.molecule, task.ao_basis, task.min_basis)

        #    Final DFT energy:            -76.426499129387
        #    Nuclear repulsion:             9.194963912878
        #    Electronic energy:           -85.593468025478
        #    Embedding energy:             -0.027995016788

        ref_scf_energy = -76.426499129387
        #ref_pe_energy = -0.027995016788

        if task.mpi_rank == mpi_master():
            scf_energy = scf_drv.get_scf_energy()
            self.assertTrue(abs(scf_energy - ref_scf_energy) < 1.0e-5)

    def test_pe_tda(self):

        try:
            import cppe
        except ImportError:
            return

        inpfile = os.path.join('inputs', 'pe_water.inp')
        if not os.path.isfile(inpfile):
            inpfile = os.path.join('python_tests', inpfile)

        potfile = os.path.join('inputs', 'pe_water.pot')
        if not os.path.isfile(potfile):
            potfile = os.path.join('python_tests', potfile)

        task = MpiTask([inpfile, None], MPI.COMM_WORLD)
        task.input_dict['method_settings']['potfile'] = potfile

        scf_drv = ScfRestrictedDriver(task.mpi_comm, task.ostream)
        scf_drv.update_settings(task.input_dict['scf'],
                                task.input_dict['method_settings'])
        scf_drv.compute(task.molecule, task.ao_basis, task.min_basis)

        tda_exci = TDAExciDriver(task.mpi_comm, task.ostream)
        tda_exci.update_settings({'nstates': 5},
                                 task.input_dict['method_settings'])
        tda_results = tda_exci.compute(task.molecule, task.ao_basis,
                                       scf_drv.scf_tensors)

        #  State Frequency   Oscillator Strength    Rotatory  Strength
        #          (eV)      Velocity     Length    Velocity    Length
        #  -----------------------------------------------------------
        #     1     8.1131     0.0667     0.0595     0.1212     0.3229
        #     2     9.9888     0.0085     0.0092    -0.1655     0.1105
        #     3    10.3548     0.0864     0.1098    -0.2905     2.5348
        #     4    10.9487     0.0001     0.0003    -0.0258    -0.0721
        #     5    11.8707     0.0026     0.0008     0.1443    -0.4440

        data = np.array([
            [8.1131, 0.0595, 0.1212],
            [9.9888, 0.0092, -0.1655],
            [10.3548, 0.1098, -0.2905],
            [10.9487, 0.0003, -0.0258],
            [11.8707, 0.0008, 0.1443],
        ])

        ref_exc_ene = data[:, 0]
        ref_osc_str = data[:, 1]
        ref_rot_str = data[:, 2]

        if task.mpi_rank == mpi_master():
            exc_ene = tda_results['eigenvalues'] * hartree_in_ev()
            osc_str = tda_results['oscillator_strengths']
            rot_str = tda_results['rotatory_strengths']

            self.assertTrue(np.max(np.abs(exc_ene - ref_exc_ene)) < 5.0e-4)
            self.assertTrue(np.max(np.abs(osc_str - ref_osc_str)) < 1.0e-4)
            self.assertTrue(np.max(np.abs(rot_str - ref_rot_str)) < 1.0e-2)


if __name__ == "__main__":
    unittest.main()
