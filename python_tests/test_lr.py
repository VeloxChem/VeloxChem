from mpi4py import MPI
import numpy as np
import unittest
import os

from veloxchem.veloxchemlib import mpi_master
from veloxchem.mpitask import MpiTask
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.lrsolver import LinearResponseSolver


class TestLR(unittest.TestCase):

    def test_lr_hf(self):

        inpfile = os.path.join('inputs', 'water.inp')
        if not os.path.isfile(inpfile):
            inpfile = os.path.join('python_tests', inpfile)

        task = MpiTask([inpfile, None], MPI.COMM_WORLD)
        task.input_dict['scf']['checkpoint_file'] = None

        scf_drv = ScfRestrictedDriver(task.mpi_comm, task.ostream)
        scf_drv.update_settings(task.input_dict['scf'],
                                task.input_dict['method_settings'])
        scf_drv.compute(task.molecule, task.ao_basis, task.min_basis)

        raw_data = """
              7.251835     -5.5859894e-17 -2.0733213e-17
            -5.5087606e-17   8.724516      2.6475728e-13
            -4.1144107e-17  2.0165388e-13   7.880586
              7.311119     -5.7629617e-17 -2.1419305e-17
            -5.6685597e-17   8.770615      2.6779629e-13
            -4.3186854e-17  2.0291932e-13   7.930012
              7.502485     -6.3579042e-17 -2.3830877e-17
            -6.2009897e-17   8.912738      2.7738668e-13
            -5.0236787e-17  2.0685122e-13   8.084815
        """
        ref_prop = np.array([float(x) for x in raw_data.split()])

        lr_solver = LinearResponseSolver(task.mpi_comm, task.ostream)
        lr_solver.update_settings({'frequencies': '0-0.15(0.05)'},
                                  task.input_dict['method_settings'])
        lr_results = lr_solver.compute(task.molecule, task.ao_basis,
                                       scf_drv.scf_tensors)

        if task.mpi_rank == mpi_master():
            prop = np.array([
                -lr_results[(a, b, w)] for w in [0.0, 0.05, 0.1] for a in 'xyz'
                for b in 'xyz'
            ])
            self.assertTrue(np.max(np.abs(prop - ref_prop)) < 1.0e-4)

    def test_lr_dft(self):

        inpfile = os.path.join('inputs', 'water.inp')
        if not os.path.isfile(inpfile):
            inpfile = os.path.join('python_tests', inpfile)

        task = MpiTask([inpfile, None], MPI.COMM_WORLD)
        task.input_dict['scf']['checkpoint_file'] = None
        task.input_dict['method_settings']['xcfun'] = 'b3lyp'

        scf_drv = ScfRestrictedDriver(task.mpi_comm, task.ostream)
        scf_drv.update_settings(task.input_dict['scf'],
                                task.input_dict['method_settings'])
        scf_drv.compute(task.molecule, task.ao_basis, task.min_basis)

        raw_data = """
              8.769176      1.4902375e-17 -2.2723900e-16
            -8.3225895e-17   9.696176      9.5276746e-14
            -1.2184939e-16  1.8169161e-13   9.066348
              8.893308      2.0111098e-17 -2.4887460e-16
            -8.1721460e-17   9.757660      8.5361466e-14
            -1.3749816e-16  1.8120910e-13   9.147369
              9.318430      4.0775332e-17 -3.3100383e-16
            -7.4493096e-17   9.948384      5.1184715e-14
            -1.9911308e-16  1.7908630e-13   9.407202
        """
        ref_prop = np.array([float(x) for x in raw_data.split()])

        lr_solver = LinearResponseSolver(task.mpi_comm, task.ostream)
        lr_solver.update_settings({'frequencies': '0-0.15(0.05)'},
                                  task.input_dict['method_settings'])
        lr_results = lr_solver.compute(task.molecule, task.ao_basis,
                                       scf_drv.scf_tensors)

        if task.mpi_rank == mpi_master():
            prop = np.array([
                -lr_results[(a, b, w)] for w in [0.0, 0.05, 0.1] for a in 'xyz'
                for b in 'xyz'
            ])
            self.assertTrue(np.max(np.abs(prop - ref_prop)) < 1.0e-4)

    def test_lr_dft_slda(self):

        inpfile = os.path.join('inputs', 'water.inp')
        if not os.path.isfile(inpfile):
            inpfile = os.path.join('python_tests', inpfile)

        task = MpiTask([inpfile, None], MPI.COMM_WORLD)
        task.input_dict['scf']['checkpoint_file'] = None
        task.input_dict['method_settings']['xcfun'] = 'slda'

        scf_drv = ScfRestrictedDriver(task.mpi_comm, task.ostream)
        scf_drv.update_settings(task.input_dict['scf'],
                                task.input_dict['method_settings'])
        scf_drv.compute(task.molecule, task.ao_basis, task.min_basis)

        raw_data = """
              9.287658     -1.1997346E-16 -9.5824875E-16
            -2.9828750E-16   9.980757      1.3461724E-12
            -7.1713889E-16  2.2342396E-13   9.449695
              9.433758     -1.2002012E-16 -1.0006654E-15
            -3.0817669E-16   10.04701      1.3973664E-12
            -7.3819715E-16  2.2624766E-13   9.541279
              9.940174     -1.1900461E-16 -1.1523585E-15
            -3.4263886E-16   10.25285      1.5697014E-12
            -8.0798267E-16  2.3530467E-13   9.836690
        """
        ref_prop = np.array([float(x) for x in raw_data.split()])

        lr_solver = LinearResponseSolver(task.mpi_comm, task.ostream)
        lr_solver.update_settings({'frequencies': '0-0.15(0.05)'},
                                  task.input_dict['method_settings'])
        lr_results = lr_solver.compute(task.molecule, task.ao_basis,
                                       scf_drv.scf_tensors)

        if task.mpi_rank == mpi_master():
            prop = np.array([
                -lr_results[(a, b, w)] for w in [0.0, 0.05, 0.1] for a in 'xyz'
                for b in 'xyz'
            ])
            self.assertTrue(np.max(np.abs(prop - ref_prop)) < 1.0e-4)

    def test_lr_hf_pe(self):

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
        task.input_dict['scf']['checkpoint_file'] = None
        task.input_dict['method_settings']['potfile'] = potfile

        scf_drv = ScfRestrictedDriver(task.mpi_comm, task.ostream)
        scf_drv.update_settings(task.input_dict['scf'],
                                task.input_dict['method_settings'])
        scf_drv.compute(task.molecule, task.ao_basis, task.min_basis)

        raw_data = """
             8.747257      0.3740598      0.2502176
            0.3740668       7.609508      0.2731078
            0.2502164      0.2731005       8.107264
             8.793038      0.3721059      0.2547393
            0.3721132       7.669906      0.2760281
            0.2547378      0.2760205       8.158091
             8.934352      0.3645516      0.2691748
            0.3645603       7.863375      0.2849957
            0.2691720      0.2849872       8.317053
        """
        ref_prop = np.array([float(x) for x in raw_data.split()])

        lr_solver = LinearResponseSolver(task.mpi_comm, task.ostream)
        lr_solver.update_settings({'frequencies': '0-0.15(0.05)'},
                                  task.input_dict['method_settings'])
        lr_results = lr_solver.compute(task.molecule, task.ao_basis,
                                       scf_drv.scf_tensors)

        if task.mpi_rank == mpi_master():
            prop = np.array([
                -lr_results[(a, b, w)] for w in [0.0, 0.05, 0.1] for a in 'xyz'
                for b in 'xyz'
            ])
            self.assertTrue(np.max(np.abs(prop - ref_prop)) < 1.0e-4)

    def test_lr_dft_pe(self):

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
        task.input_dict['scf']['checkpoint_file'] = None
        task.input_dict['method_settings']['xcfun'] = 'b3lyp'
        task.input_dict['method_settings']['potfile'] = potfile

        scf_drv = ScfRestrictedDriver(task.mpi_comm, task.ostream)
        scf_drv.update_settings(task.input_dict['scf'],
                                task.input_dict['method_settings'])
        scf_drv.compute(task.molecule, task.ao_basis, task.min_basis)

        raw_data = """
             9.715185      0.2695542      0.3272248
            0.2695358       9.122653      0.2927027
            0.3272155      0.2927145       9.298530
             9.777315      0.2585090      0.3352403
            0.2584891       9.241745      0.2951616
            0.3352328      0.2951721       9.379825
             9.971221      0.2175020      0.3616347
            0.2174766       9.642433      0.3015763
            0.3616332      0.3015831       9.639443
        """
        ref_prop = np.array([float(x) for x in raw_data.split()])

        lr_solver = LinearResponseSolver(task.mpi_comm, task.ostream)
        lr_solver.update_settings({'frequencies': '0-0.15(0.05)'},
                                  task.input_dict['method_settings'])
        lr_results = lr_solver.compute(task.molecule, task.ao_basis,
                                       scf_drv.scf_tensors)

        if task.mpi_rank == mpi_master():
            prop = np.array([
                -lr_results[(a, b, w)] for w in [0.0, 0.05, 0.1] for a in 'xyz'
                for b in 'xyz'
            ])
            self.assertTrue(np.max(np.abs(prop - ref_prop)) < 1.0e-4)


if __name__ == "__main__":
    unittest.main()
