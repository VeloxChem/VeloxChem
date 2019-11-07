from mpi4py import MPI
import numpy as np
import unittest
import os

from veloxchem.veloxchemlib import mpi_master
from veloxchem.veloxchemlib import hartree_in_ev
from veloxchem.mpitask import MpiTask
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.tdaexcidriver import TDAExciDriver


class TestTDA(unittest.TestCase):

    def test_tda_hf(self):

        inpfile = os.path.join('inputs', 'water.inp')
        if not os.path.isfile(inpfile):
            inpfile = os.path.join('python_tests', inpfile)

        task = MpiTask([inpfile, None], MPI.COMM_WORLD)

        scf_drv = ScfRestrictedDriver(task.mpi_comm, task.ostream)
        scf_drv.update_settings(task.input_dict['scf'],
                                task.input_dict['method_settings'])
        scf_drv.compute(task.molecule, task.ao_basis, task.min_basis)

        #  State Frequency   Oscillator Strength    Rotatory  Strength
        #          (eV)      Velocity     Length    Velocity    Length
        #  -----------------------------------------------------------
        raw_data = """
            1     8.7847     0.0695     0.0530     0.0000     0.0000
            2    10.4740     0.0000     0.0000    -0.0000    -0.0000
            3    11.1634     0.0869     0.1070    -0.0000    -0.0000
            4    12.1889     0.0014     0.0045    -0.0000    -0.0000
            5    12.8274     0.0117     0.0250     0.0000     0.0000
        """
        lines = raw_data.split(os.linesep)[1:-1]

        ref_exc_ene = [float(line.split()[1]) for line in lines]
        ref_osc_str = [float(line.split()[3]) for line in lines]

        tda_solver = TDAExciDriver(task.mpi_comm, task.ostream)
        tda_solver.update_settings({'nstates': 5},
                                   task.input_dict['method_settings'])
        tda_results = tda_solver.compute(task.molecule, task.ao_basis,
                                         scf_drv.scf_tensors)

        if task.mpi_rank == mpi_master():
            exc_ene = tda_results['eigenvalues'] * hartree_in_ev()
            osc_str = tda_results['oscillator_strengths']

            self.assertTrue(np.max(np.abs(exc_ene - ref_exc_ene)) < 5.0e-4)
            self.assertTrue(np.max(np.abs(osc_str - ref_osc_str)) < 5.0e-4)

    def test_tda_dft(self):

        inpfile = os.path.join('inputs', 'water.inp')
        if not os.path.isfile(inpfile):
            inpfile = os.path.join('python_tests', inpfile)

        task = MpiTask([inpfile, None], MPI.COMM_WORLD)
        task.input_dict['method_settings']['xcfun'] = 'b3lyp'

        scf_drv = ScfRestrictedDriver(task.mpi_comm, task.ostream)
        scf_drv.update_settings(task.input_dict['scf'],
                                task.input_dict['method_settings'])
        scf_drv.compute(task.molecule, task.ao_basis, task.min_basis)

        #  State Frequency   Oscillator Strength    Rotatory  Strength
        #          (eV)      Velocity     Length    Velocity    Length
        #  -----------------------------------------------------------
        raw_data = """
            1     6.9987     0.0571     0.0537     0.0000     0.0000
            2     8.4291     0.0000     0.0000    -0.0000    -0.0000
            3     9.2135     0.0611     0.0906     0.0000     0.0000
            4    10.3003     0.0005     0.0000    -0.0000    -0.0000
            5    10.6270     0.0044     0.0127    -0.0000    -0.0000
        """
        lines = raw_data.split(os.linesep)[1:-1]

        ref_exc_ene = [float(line.split()[1]) for line in lines]
        ref_osc_str = [float(line.split()[3]) for line in lines]

        tda_solver = TDAExciDriver(task.mpi_comm, task.ostream)
        tda_solver.update_settings({'nstates': 5},
                                   task.input_dict['method_settings'])
        tda_results = tda_solver.compute(task.molecule, task.ao_basis,
                                         scf_drv.scf_tensors)

        if task.mpi_rank == mpi_master():
            exc_ene = tda_results['eigenvalues'] * hartree_in_ev()
            osc_str = tda_results['oscillator_strengths']

            self.assertTrue(np.max(np.abs(exc_ene - ref_exc_ene)) < 5.0e-4)
            self.assertTrue(np.max(np.abs(osc_str - ref_osc_str)) < 5.0e-4)

    def test_tda_dft_slda(self):

        inpfile = os.path.join('inputs', 'water.inp')
        if not os.path.isfile(inpfile):
            inpfile = os.path.join('python_tests', inpfile)

        task = MpiTask([inpfile, None], MPI.COMM_WORLD)
        task.input_dict['method_settings']['xcfun'] = 'slda'

        scf_drv = ScfRestrictedDriver(task.mpi_comm, task.ostream)
        scf_drv.update_settings(task.input_dict['scf'],
                                task.input_dict['method_settings'])
        scf_drv.compute(task.molecule, task.ao_basis, task.min_basis)

        #  State Frequency   Oscillator Strength    Rotatory  Strength
        #          (eV)      Velocity     Length    Velocity    Length
        #  -----------------------------------------------------------
        raw_data = """
            1     6.7828     0.0588     0.0561    -0.0000     0.0000
            2     8.2221     0.0000     0.0000     0.0000     0.0000
            3     8.9101     0.0603     0.0901    -0.0000    -0.0000
            4    10.1323     0.0014     0.0003     0.0000     0.0000
            5    10.3444     0.0036     0.0115    -0.0000    -0.0000
        """
        lines = raw_data.split(os.linesep)[1:-1]

        ref_exc_ene = [float(line.split()[1]) for line in lines]
        ref_osc_str = [float(line.split()[3]) for line in lines]

        tda_solver = TDAExciDriver(task.mpi_comm, task.ostream)
        tda_solver.update_settings({'nstates': 5},
                                   task.input_dict['method_settings'])
        tda_results = tda_solver.compute(task.molecule, task.ao_basis,
                                         scf_drv.scf_tensors)

        if task.mpi_rank == mpi_master():
            exc_ene = tda_results['eigenvalues'] * hartree_in_ev()
            osc_str = tda_results['oscillator_strengths']

            self.assertTrue(np.max(np.abs(exc_ene - ref_exc_ene)) < 5.0e-4)
            self.assertTrue(np.max(np.abs(osc_str - ref_osc_str)) < 5.0e-4)

    def test_tda_hf_pe(self):

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

        #  State Frequency   Oscillator Strength    Rotatory  Strength
        #          (eV)      Velocity     Length    Velocity    Length
        #  -----------------------------------------------------------
        raw_data = """
            1     9.1467     0.0758     0.0602     0.1440     0.3887
            2    11.1635     0.0113     0.0121    -0.4940    -0.2295
            3    11.4450     0.1062     0.1234    -0.1259     2.8161
            4    11.9792     0.0004     0.0013     0.0916     0.1963
            5    12.8941     0.0007     0.0006     0.0132    -0.4366
        """
        lines = raw_data.split(os.linesep)[1:-1]

        ref_exc_ene = [float(line.split()[1]) for line in lines]
        ref_osc_str = [float(line.split()[3]) for line in lines]
        ref_rot_str = [float(line.split()[4]) for line in lines]

        tda_solver = TDAExciDriver(task.mpi_comm, task.ostream)
        tda_solver.update_settings({'nstates': 5},
                                   task.input_dict['method_settings'])
        tda_results = tda_solver.compute(task.molecule, task.ao_basis,
                                         scf_drv.scf_tensors)

        if task.mpi_rank == mpi_master():
            exc_ene = tda_results['eigenvalues'] * hartree_in_ev()
            osc_str = tda_results['oscillator_strengths']
            rot_str = tda_results['rotatory_strengths']

            self.assertTrue(np.max(np.abs(exc_ene - ref_exc_ene)) < 5.0e-4)
            self.assertTrue(np.max(np.abs(osc_str - ref_osc_str)) < 5.0e-4)
            self.assertTrue(np.max(np.abs(rot_str - ref_rot_str)) < 1.0e-2)

    def test_tda_dft_pe(self):

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
        task.input_dict['method_settings']['xcfun'] = 'b3lyp'
        task.input_dict['method_settings']['potfile'] = potfile

        scf_drv = ScfRestrictedDriver(task.mpi_comm, task.ostream)
        scf_drv.update_settings(task.input_dict['scf'],
                                task.input_dict['method_settings'])
        scf_drv.compute(task.molecule, task.ao_basis, task.min_basis)

        #  State Frequency   Oscillator Strength    Rotatory  Strength
        #          (eV)      Velocity     Length    Velocity    Length
        #  -----------------------------------------------------------
        raw_data = """
            1     7.3245     0.0609     0.0582     0.1219     0.3536
            2     9.0048     0.0067     0.0073    -0.0937     0.1897
            3     9.5204     0.0742     0.1000    -0.2837     2.3062
            4    10.1845     0.0003     0.0003    -0.1075    -0.2451
            5    11.0440     0.0076     0.0029     0.1461    -0.5130
        """
        lines = raw_data.split(os.linesep)[1:-1]

        ref_exc_ene = [float(line.split()[1]) for line in lines]
        ref_osc_str = [float(line.split()[3]) for line in lines]
        ref_rot_str = [float(line.split()[4]) for line in lines]

        tda_solver = TDAExciDriver(task.mpi_comm, task.ostream)
        tda_solver.update_settings({'nstates': 5},
                                   task.input_dict['method_settings'])
        tda_results = tda_solver.compute(task.molecule, task.ao_basis,
                                         scf_drv.scf_tensors)

        if task.mpi_rank == mpi_master():
            exc_ene = tda_results['eigenvalues'] * hartree_in_ev()
            osc_str = tda_results['oscillator_strengths']
            rot_str = tda_results['rotatory_strengths']

            self.assertTrue(np.max(np.abs(exc_ene - ref_exc_ene)) < 5.0e-4)
            self.assertTrue(np.max(np.abs(osc_str - ref_osc_str)) < 5.0e-4)
            self.assertTrue(np.max(np.abs(rot_str - ref_rot_str)) < 1.0e-2)


if __name__ == "__main__":
    unittest.main()
