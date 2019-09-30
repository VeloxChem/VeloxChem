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

        inpfile = os.path.join('inputs', 'pe_water.inp')
        if not os.path.isfile(inpfile):
            inpfile = os.path.join('python_tests', inpfile)

        potfile = os.path.join('inputs', 'pe_water.pot')
        if not os.path.isfile(potfile):
            potfile = os.path.join('python_tests', potfile)

        task = MpiTask([inpfile, None], MPI.COMM_WORLD)
        task.input_dict['method_settings']['potfile'] = potfile
        task.input_dict['method_settings']['xcfun'] = 'slater'

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
        #     1     6.2079     0.0536     0.0554     0.1250     0.5057
        #     2     7.5234     0.0049     0.0056    -0.0673     0.2430
        #     3     8.3137     0.0594     0.0841    -0.2521     1.5850
        #     4     9.1766     0.0030     0.0029    -0.1330    -0.5048
        #     5     9.6380     0.0059     0.0093     0.1047    -0.7910

        data = np.array([
            [6.2079, 0.0554, 0.1250],
            [7.5234, 0.0056, -0.0673],
            [8.3137, 0.0841, -0.2521],
            [9.1766, 0.0029, -0.1330],
            [9.6380, 0.0093, 0.1047],
        ])

        ref_exc_ene = data[:, 0]
        ref_osc_str = data[:, 1]
        ref_rot_str = data[:, 2]

        if task.mpi_rank == mpi_master():
            exc_ene = tda_results['eigenvalues'] * hartree_in_ev()
            osc_str = tda_results['oscillator_strengths']
            rot_str = tda_results['rotatory_strengths']

            self.assertTrue(np.max(np.abs(exc_ene - ref_exc_ene)) < 1.0e-4)
            self.assertTrue(np.max(np.abs(osc_str - ref_osc_str)) < 1.0e-4)
            self.assertTrue(np.max(np.abs(rot_str - ref_rot_str)) < 1.0e-2)


if __name__ == "__main__":
    unittest.main()
