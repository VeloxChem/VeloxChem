from mpi4py import MPI
import numpy as np
import unittest
from pathlib import Path

from veloxchem.veloxchemlib import AODensityMatrix
from veloxchem.veloxchemlib import denmat
from veloxchem.veloxchemlib import mpi_master
from veloxchem.mpitask import MpiTask
from veloxchem.cubicgrid import CubicGrid
from veloxchem.visualizationdriver import VisualizationDriver
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.tdaexcidriver import TDAExciDriver
from veloxchem.lreigensolver import LinearResponseEigenSolver


class TestDetachAttach(unittest.TestCase):

    def run_detach_attach(self, inpfile, xcfun_label, nstates, ref_cube_vals,
                          flag):

        task = MpiTask([inpfile, None], MPI.COMM_WORLD)
        task.input_dict['scf']['checkpoint_file'] = None

        if xcfun_label is not None:
            task.input_dict['method_settings']['xcfun'] = xcfun_label

        scf_drv = ScfRestrictedDriver(task.mpi_comm, task.ostream)
        scf_drv.update_settings(task.input_dict['scf'],
                                task.input_dict['method_settings'])
        scf_drv.compute(task.molecule, task.ao_basis, task.min_basis)

        if flag == 'tda':
            rsp_drv = TDAExciDriver(task.mpi_comm, task.ostream)
        elif flag == 'rpa':
            rsp_drv = LinearResponseEigenSolver(task.mpi_comm, task.ostream)
        rsp_drv.update_settings({'nstates': nstates},
                                task.input_dict['method_settings'])
        rsp_results = rsp_drv.compute(task.molecule, task.ao_basis,
                                      scf_drv.scf_tensors)

        nocc = task.molecule.number_of_alpha_electrons()

        if task.mpi_rank == mpi_master():
            mo = scf_drv.scf_tensors['C']
            mo_occ = mo[:, :nocc]
            mo_vir = mo[:, nocc:]

            eig_vecs = rsp_results['eigenvectors']

            cube_vals = []

        for s in range(nstates):
            if task.mpi_rank == mpi_master():
                eig_vec = eig_vecs[:, s].copy()

                if flag == 'tda':
                    t_mat = eig_vec.reshape(mo_occ.shape[1], mo_vir.shape[1])
                    dens_D, dens_A = rsp_drv.get_detach_attach_densities(
                        t_mat, mo_occ, mo_vir)
                    dens_DA = AODensityMatrix([dens_D, dens_A], denmat.rest)

                elif flag == 'rpa':
                    z_mat = eig_vec[:eig_vec.shape[0] // 2].reshape(
                        mo_occ.shape[1], mo_vir.shape[1]) * np.sqrt(2.0)
                    y_mat = eig_vec[eig_vec.shape[0] // 2:].reshape(
                        mo_occ.shape[1], mo_vir.shape[1]) * np.sqrt(2.0)

                    dens_Dz, dens_Az = rsp_drv.get_detach_attach_densities(
                        z_mat, mo_occ, mo_vir)
                    dens_Dy, dens_Ay = rsp_drv.get_detach_attach_densities(
                        y_mat, mo_occ, mo_vir)

                    dens_DA = AODensityMatrix(
                        [dens_Dz + dens_Dy, dens_Az + dens_Ay], denmat.rest)
            else:
                dens_DA = AODensityMatrix()
            dens_DA.broadcast(task.mpi_rank, task.mpi_comm)

            vis_drv = VisualizationDriver(task.mpi_comm)
            grid = CubicGrid([-0.1, -0.1, -0.1], [0.3, 0.3, 0.2], [2, 2, 3])

            vis_drv.compute(grid, task.molecule, task.ao_basis, dens_DA, 0,
                            'alpha')
            if task.mpi_rank == mpi_master():
                cube_vals.append(grid.values_to_numpy().flatten())

            vis_drv.compute(grid, task.molecule, task.ao_basis, dens_DA, 1,
                            'alpha')
            if task.mpi_rank == mpi_master():
                cube_vals.append(grid.values_to_numpy().flatten())

        if task.mpi_rank == mpi_master():
            cube_vals = np.array(cube_vals)

            for s in range(cube_vals.shape[0]):
                if np.vdot(cube_vals[s, :], ref_cube_vals[s, :]) < 0.0:
                    cube_vals[s, :] *= -1.0
            self.assertTrue(np.max(np.abs(cube_vals - ref_cube_vals)) < 1.0e-4)

    def test_detach_attach_tda(self):

        here = Path(__file__).parent
        inpfile = str(here / 'inputs' / 'c2h4.inp')

        xcfun_label = None

        nstates = 5

        raw_vals = """
           -3.81621e-03  -3.81621e-03  -4.09509e-03
           -3.66394e-03  -3.66394e-03  -3.91786e-03
           -6.23259e-03  -6.23259e-03  -7.04599e-03
           -5.98682e-03  -5.98682e-03  -6.73525e-03
            1.51129e-04   1.51129e-04   1.20872e-03
            1.46478e-04   1.46478e-04   1.14244e-03
            1.88781e-04   1.88781e-04   1.60112e-03
            1.81843e-04   1.81843e-04   1.50992e-03
           -3.39078e-04  -3.39078e-04  -4.73343e-04
           -3.45303e-04  -3.45303e-04  -6.38456e-04
           -1.25629e-03  -1.25629e-03  -1.58539e-03
           -1.22529e-03  -1.22529e-03  -1.68720e-03
            2.57158e-04   2.57158e-04   4.24425e-04
            9.43715e-04   9.43715e-04   1.21507e-03
            2.86828e-04   2.86828e-04   7.91893e-04
            9.43174e-04   9.43174e-04   1.52054e-03
           -9.71119e-04  -9.71119e-04  -1.18837e-03
           -9.32354e-04  -9.32354e-04  -1.13312e-03
           -3.66594e-03  -3.66594e-03  -4.46763e-03
           -3.52138e-03  -3.52138e-03  -4.26279e-03
            6.01344e-05   6.01344e-05   1.45026e-04
            5.44916e-05   5.44916e-05   1.22518e-04
            4.57288e-05   4.57288e-05   1.08057e-04
            4.16791e-05   4.16791e-05   9.10495e-05
           -6.31689e-04  -6.31689e-04  -8.06054e-04
           -6.16541e-04  -6.16541e-04  -8.64765e-04
           -2.41362e-03  -2.41362e-03  -2.97229e-03
           -2.32789e-03  -2.32789e-03  -2.92570e-03
            4.61014e-04   4.61014e-04   6.05328e-04
            1.74446e-03   1.74446e-03   2.10302e-03
            4.62356e-04   4.62356e-04   7.70773e-04
            1.69191e-03   1.69191e-03   2.18564e-03
           -7.10253e-02  -7.10253e-02  -7.33318e-02
           -6.77306e-02  -6.77306e-02  -6.95517e-02
           -6.90265e-02  -6.90265e-02  -7.11074e-02
           -6.58643e-02  -6.58643e-02  -6.74921e-02
            2.76857e-05   2.76857e-05   2.41246e-04
            2.77213e-05   2.77213e-05   2.41045e-04
            7.21540e-05   7.21540e-05   6.80185e-04
            6.94721e-05   6.94721e-05   6.51147e-04
        """
        cube_vals = np.array([float(x) for x in raw_vals.split()])
        cube_vals = cube_vals.reshape(nstates * 2, -1)

        self.run_detach_attach(inpfile, xcfun_label, nstates, cube_vals, 'tda')

    def test_detach_attach_rpa(self):

        here = Path(__file__).parent
        inpfile = str(here / 'inputs' / 'c2h4.inp')

        xcfun_label = 'b3lyp'

        nstates = 3

        raw_vals = """
           -3.27476e-05  -3.27476e-05  -1.32379e-04
           -5.99181e-05  -5.99181e-05  -4.18994e-04
           -8.17569e-05  -8.17569e-05  -1.86632e-04
           -1.05149e-04  -1.05149e-04  -4.50826e-04
            3.41319e-05   3.41319e-05   2.42345e-04
            6.59041e-05   6.59041e-05   2.65393e-04
            9.52284e-05   9.52284e-05   8.65482e-04
            1.21892e-04   1.21892e-04   8.45480e-04
           -3.73966e-03  -3.73967e-03  -4.05024e-03
           -3.59327e-03  -3.59329e-03  -3.87652e-03
           -6.30145e-03  -6.30147e-03  -7.23666e-03
           -6.05426e-03  -6.05429e-03  -6.91549e-03
            1.31340e-04   1.31439e-04   1.11607e-03
            1.36796e-04   1.36866e-04   1.06198e-03
            1.78907e-04   1.79094e-04   1.60605e-03
            1.81161e-04   1.81312e-04   1.52018e-03
           -9.51467e-04  -9.51467e-04  -1.18572e-03
           -9.13013e-04  -9.13013e-04  -1.12955e-03
           -3.58911e-03  -3.58911e-03  -4.45325e-03
           -3.44578e-03  -3.44578e-03  -4.24514e-03
            3.21810e-04   3.21810e-04   5.24184e-04
            2.94778e-04   2.94778e-04   4.61132e-04
            2.73616e-04   2.73616e-04   4.34842e-04
            2.51156e-04   2.51156e-04   3.82716e-04
        """
        cube_vals = np.array([float(x) for x in raw_vals.split()])
        cube_vals = cube_vals.reshape(nstates * 2, -1)

        self.run_detach_attach(inpfile, xcfun_label, nstates, cube_vals, 'rpa')


if __name__ == "__main__":
    unittest.main()
