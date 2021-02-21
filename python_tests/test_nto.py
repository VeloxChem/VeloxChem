from mpi4py import MPI
from pathlib import Path
import numpy as np
import unittest
import tempfile

from veloxchem.veloxchemlib import mpi_master
from veloxchem.veloxchemlib import denmat
from veloxchem.mpitask import MpiTask
from veloxchem.cubicgrid import CubicGrid
from veloxchem.aodensitymatrix import AODensityMatrix
from veloxchem.molecularorbitals import MolecularOrbitals
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.tdaexcidriver import TDAExciDriver
from veloxchem.lreigensolver import LinearResponseEigenSolver


class TestNTO(unittest.TestCase):

    def run_nto(self, inpfile, xcfun_label, ref_eig_vals, ref_nto_lambdas,
                ref_nto_cube_vals, ref_dens_cube_vals, flag):

        task = MpiTask([inpfile, None], MPI.COMM_WORLD)
        task.input_dict['scf']['checkpoint_file'] = None

        if xcfun_label is not None:
            task.input_dict['method_settings']['xcfun'] = xcfun_label

        # run SCF

        scf_drv = ScfRestrictedDriver(task.mpi_comm, task.ostream)
        scf_drv.update_settings(task.input_dict['scf'],
                                task.input_dict['method_settings'])
        scf_drv.compute(task.molecule, task.ao_basis, task.min_basis)

        # run TDA/RPA

        if flag == 'tda':
            rsp_drv = TDAExciDriver(task.mpi_comm, task.ostream)
        elif flag == 'rpa':
            rsp_drv = LinearResponseEigenSolver(task.mpi_comm, task.ostream)
        rsp_drv.update_settings({'nstates': ref_eig_vals.shape[0]},
                                task.input_dict['method_settings'])
        rsp_results = rsp_drv.compute(task.molecule, task.ao_basis,
                                      scf_drv.scf_tensors)

        # get eigenvalues and eigenvectors

        nocc = task.molecule.number_of_alpha_electrons()

        if task.mpi_rank == mpi_master():
            mo = scf_drv.scf_tensors['C']
            mo_occ = mo[:, :nocc]
            mo_vir = mo[:, nocc:]

            nocc = mo_occ.shape[1]
            nvir = mo_vir.shape[1]

            eig_vals = rsp_results['eigenvalues']
            eig_vecs = rsp_results['eigenvectors']

            nto_lambdas = []
            nto_cube_vals = []
            dens_cube_vals = []

        for s in range(ref_eig_vals.shape[0]):

            # calculate NTOs and densities

            if task.mpi_rank == mpi_master():
                eig_vec = eig_vecs[:, s].copy()

                if flag == 'tda':
                    z_mat = eig_vec.reshape(nocc, nvir)
                    dens_D, dens_A = rsp_drv.get_detach_attach_densities(
                        z_mat, None, mo_occ, mo_vir)

                elif flag == 'rpa':
                    z_mat = eig_vec[:eig_vec.shape[0] // 2].reshape(
                        nocc, nvir) * np.sqrt(2.0)
                    y_mat = eig_vec[eig_vec.shape[0] // 2:].reshape(
                        nocc, nvir) * np.sqrt(2.0)
                    dens_D, dens_A = rsp_drv.get_detach_attach_densities(
                        z_mat, y_mat, mo_occ, mo_vir)

                lam_diag, nto_mo = rsp_drv.get_nto(z_mat, mo_occ, mo_vir)
                dens_DA = AODensityMatrix([dens_D, dens_A], denmat.rest)

            else:
                lam_diag = None
                nto_mo = MolecularOrbitals()
                dens_DA = AODensityMatrix()

            lam_diag = task.mpi_comm.bcast(lam_diag, root=mpi_master())
            nto_mo.broadcast(task.mpi_rank, task.mpi_comm)
            dens_DA.broadcast(task.mpi_rank, task.mpi_comm)

            # check sum of lambdas and number of electrons

            if task.mpi_rank == mpi_master():
                nto_lambdas.append(lam_diag[0])

                if flag == 'tda':
                    self.assertAlmostEqual(np.sum(lam_diag), 1.0, 8)

                self.assertAlmostEqual(
                    np.sum(dens_D * scf_drv.scf_tensors['S']), -1.0, 8)
                self.assertAlmostEqual(
                    np.sum(dens_A * scf_drv.scf_tensors['S']), 1.0, 8)

            # set up grid points

            grid = CubicGrid([-0.1, -0.1, -0.1], [0.3, 0.3, 0.2], [2, 2, 3])

            # compute NTOs on grid points

            with tempfile.TemporaryDirectory() as temp_dir:
                rsp_drv.filename = str(Path(temp_dir, 'test'))

                nto_cube_fnames = rsp_drv.write_nto_cubes(
                    grid, task.molecule, task.ao_basis, s, lam_diag, nto_mo, 1)

                dens_cube_fnames = rsp_drv.write_detach_attach_cubes(
                    grid, task.molecule, task.ao_basis, s, dens_DA)

                if task.mpi_rank == mpi_master():

                    for fname in nto_cube_fnames:
                        read_grid = CubicGrid.read_cube(fname)
                        nto_cube_vals.append(
                            read_grid.values_to_numpy().flatten())

                    for fname in dens_cube_fnames:
                        read_grid = CubicGrid.read_cube(fname)
                        dens_cube_vals.append(
                            read_grid.values_to_numpy().flatten())

        # compare NTO and densities on grid points with reference

        if task.mpi_rank == mpi_master():
            nto_lambdas = np.array(nto_lambdas)
            nto_cube_vals = np.array(nto_cube_vals)
            dens_cube_vals = np.array(dens_cube_vals)

            thresh = 1.0e-6 if xcfun_label is None else 1.0e-5
            self.assertTrue(np.max(np.abs(eig_vals - ref_eig_vals)) < thresh)

            thresh = 1.0e-4 if xcfun_label is None else 1.0e-3
            self.assertTrue(
                np.max(np.abs(nto_lambdas - ref_nto_lambdas)) < thresh)

            for s in range(nto_cube_vals.shape[0]):
                if np.vdot(nto_cube_vals[s, :], ref_nto_cube_vals[s, :]) < 0.0:
                    nto_cube_vals[s, :] *= -1.0
                if np.vdot(dens_cube_vals[s, :],
                           ref_dens_cube_vals[s, :]) < 0.0:
                    dens_cube_vals[s, :] *= -1.0
            self.assertTrue(
                np.max(np.abs(nto_cube_vals - ref_nto_cube_vals)) < 1.0e-4)
            self.assertTrue(
                np.max(np.abs(dens_cube_vals - ref_dens_cube_vals)) < 1.0e-4)

    def test_nto_tda(self):

        here = Path(__file__).parent
        inpfile = str(here / 'inputs' / 'c2h4.inp')

        xcfun_label = None
        nstates = 5

        raw_vals = """
            0.31643416  0.939996
            0.34674141  0.658068
            0.34866292  0.998586
            0.35831265  0.660147
            0.38623776  0.611341
        """
        vals = np.array([float(x) for x in raw_vals.split()])
        vals = vals.reshape(nstates, -1)

        eig_vals = vals[:, 0]
        nto_lambdas = vals[:, 1]

        raw_vals = """
            3.08300E-02  3.08300E-02  3.41697E-02
            3.02083E-02  3.02083E-02  3.33643E-02
           -6.04152E-02 -6.04152E-02 -6.67266E-02
           -5.92120E-02 -5.92120E-02 -6.51781E-02
            4.19373E-03 -4.19373E-03 -1.33649E-02
            4.06161E-03 -4.06161E-03 -1.29179E-02
           -8.12234E-03  8.12234E-03  2.58328E-02
           -7.87085E-03  7.87085E-03  2.49837E-02
            3.11244E-03 -3.11244E-03 -1.00451E-02
           -6.02484E-03  6.02484E-03  1.94061E-02
            3.00892E-03 -3.00892E-03 -9.69226E-03
           -5.82780E-03  5.82780E-03  1.87359E-02
            4.54510E-03 -4.54510E-03 -1.45426E-02
            4.39760E-03 -4.39760E-03 -1.40423E-02
           -8.79436E-03  8.79436E-03  2.80818E-02
           -8.51372E-03  8.51372E-03  2.71320E-02
           -3.08300E-02 -3.08300E-02 -3.41697E-02
           -3.02083E-02 -3.02083E-02 -3.33643E-02
            6.04152E-02  6.04152E-02  6.67266E-02
            5.92120E-02  5.92120E-02  6.51781E-02
           -7.71797E-03 -7.71797E-03 -1.20185E-02
           -7.34438E-03 -7.34438E-03 -1.10426E-02
           -6.58723E-03 -6.58723E-03 -1.02655E-02
           -6.27868E-03 -6.27868E-03 -9.40533E-03
           -3.08300E-02 -3.08300E-02 -3.41697E-02
           -3.02083E-02 -3.02083E-02 -3.33643E-02
            6.04152E-02  6.04152E-02  6.67266E-02
            5.92120E-02  5.92120E-02  6.51781E-02
            2.62157E-02  2.62157E-02  2.84788E-02
           -5.13050E-02 -5.13050E-02 -5.55604E-02
            2.56648E-02  2.56648E-02  2.77928E-02
           -5.02375E-02 -5.02375E-02 -5.42392E-02
           -3.39884E-01 -3.39884E-01 -3.45196E-01
           -3.31899E-01 -3.31899E-01 -3.36173E-01
           -3.32479E-01 -3.32479E-01 -3.36811E-01
           -3.24750E-01 -3.24750E-01 -3.28112E-01
            4.99796E-03 -4.99796E-03 -1.59590E-02
            4.83409E-03 -4.83409E-03 -1.54049E-02
           -9.66732E-03  9.66732E-03  3.08069E-02
           -9.35542E-03  9.35542E-03  2.97545E-02
        """
        nto_cube_vals = np.array([float(x) for x in raw_vals.split()])
        nto_cube_vals = nto_cube_vals.reshape(nstates * 2, -1)

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
        dens_cube_vals = np.array([float(x) for x in raw_vals.split()])
        dens_cube_vals = dens_cube_vals.reshape(nstates * 2, -1)

        self.run_nto(inpfile, xcfun_label, eig_vals, nto_lambdas, nto_cube_vals,
                     dens_cube_vals, 'tda')

    def test_nto_rpa(self):

        here = Path(__file__).parent
        inpfile = str(here / 'inputs' / 'c2h4.inp')

        xcfun_label = 'b3lyp'
        nstates = 5

        raw_vals = """
            0.30465317  9.81738000e-01
            0.30537184  9.90163214e-01
            0.31770840  9.97655472e-01
            0.33176863  9.81756714e-01
            0.35264869  8.96190203e-01
        """
        vals = np.array([float(x) for x in raw_vals.split()])
        vals = vals.reshape(nstates, -1)

        eig_vals = vals[:, 0]
        nto_lambdas = vals[:, 1]

        raw_vals = """
            3.21430E-03 -3.21430E-03 -1.03578E-02
           -6.22010E-03  6.22010E-03  2.00043E-02
            3.10657E-03 -3.10657E-03 -9.99144E-03
           -6.01501E-03  6.01501E-03  1.93083E-02
            4.77717E-03 -4.77717E-03 -1.52606E-02
            4.62157E-03 -4.62157E-03 -1.47339E-02
           -9.24233E-03  9.24233E-03  2.94651E-02
           -8.94621E-03  8.94621E-03  2.84649E-02
            3.05252E-02  3.05252E-02  3.41449E-02
            2.99013E-02  2.99013E-02  3.33239E-02
           -5.98010E-02 -5.98010E-02 -6.66455E-02
           -5.85944E-02 -5.85944E-02 -6.50684E-02
            4.36168E-03 -4.36168E-03 -1.39084E-02
            4.22290E-03 -4.22290E-03 -1.34388E-02
           -8.44485E-03  8.44485E-03  2.68744E-02
           -8.18069E-03  8.18069E-03  2.59826E-02
            3.05252E-02  3.05252E-02  3.41449E-02
            2.99013E-02  2.99013E-02  3.33239E-02
           -5.98010E-02 -5.98010E-02 -6.66455E-02
           -5.85944E-02 -5.85944E-02 -6.50684E-02
            1.78746E-02  1.78746E-02  2.28482E-02
            1.71020E-02  1.71020E-02  2.14232E-02
            1.64037E-02  1.64037E-02  2.07118E-02
            1.57078E-02  1.57078E-02  1.94172E-02
            3.05252E-02  3.05252E-02  3.41449E-02
            2.99013E-02  2.99013E-02  3.33239E-02
           -5.98010E-02 -5.98010E-02 -6.66455E-02
           -5.85944E-02 -5.85944E-02 -6.50684E-02
           -2.38075E-02 -2.38075E-02 -2.61538E-02
            4.65993E-02  4.65993E-02  5.10177E-02
           -2.33117E-02 -2.33117E-02 -2.55212E-02
            4.56395E-02  4.56395E-02  4.98006E-02
           -3.05252E-02 -3.05252E-02 -3.41449E-02
           -2.99013E-02 -2.99013E-02 -3.33239E-02
            5.98010E-02  5.98010E-02  6.66455E-02
            5.85944E-02  5.85944E-02  6.50684E-02
            3.06410E-03 -3.06410E-03 -7.40479E-03
            3.37437E-03 -3.37437E-03 -8.51925E-03
            3.39094E-03 -3.39094E-03 -8.56887E-03
            3.67531E-03 -3.67531E-03 -9.58616E-03
        """
        nto_cube_vals = np.array([float(x) for x in raw_vals.split()])
        nto_cube_vals = nto_cube_vals.reshape(nstates * 2, -1)

        raw_vals = """
           -2.90190e-05  -2.90190e-05  -1.28491e-04
           -5.62718e-05  -5.62718e-05  -4.15115e-04
           -7.67176e-05  -7.67176e-05  -1.81100e-04
           -1.00246e-04  -1.00246e-04  -4.45388e-04
            3.38531e-05   3.38531e-05   2.42054e-04
            6.48336e-05   6.48336e-05   2.64405e-04
            9.49380e-05   9.49380e-05   8.65055e-04
            1.20815e-04   1.20815e-04   8.44361e-04
           -2.55384e-03  -2.55384e-03  -2.82839e-03
           -2.45412e-03  -2.45412e-03  -2.70607e-03
           -5.04742e-03  -5.04742e-03  -5.92456e-03
           -4.84917e-03  -4.84917e-03  -5.65868e-03
            7.37806e-05   7.37806e-05   6.50915e-04
            7.71850e-05   7.71850e-05   6.20781e-04
            1.20455e-04   1.20455e-04   1.12571e-03
            1.20853e-04   1.20853e-04   1.06480e-03
           -9.35000e-04  -9.35000e-04  -1.16864e-03
           -8.97194e-04  -8.97194e-04  -1.11320e-03
           -3.57282e-03  -3.57282e-03  -4.43632e-03
           -3.43012e-03  -3.43012e-03  -4.22891e-03
            3.16751e-04   3.16751e-04   5.19480e-04
            2.89827e-04   2.89827e-04   4.56550e-04
            2.67603e-04   2.67603e-04   4.29115e-04
            2.45298e-04   2.45298e-04   3.77157e-04
           -9.15314e-04  -9.15314e-04  -1.14688e-03
           -8.78782e-04  -8.78782e-04  -1.09781e-03
           -3.51135e-03  -3.51135e-03  -4.36264e-03
           -3.37155e-03  -3.37155e-03  -4.16365e-03
            5.56773e-04   5.56773e-04   6.75673e-04
            2.13204e-03   2.13204e-03   2.55907e-03
            5.34839e-04   5.34839e-04   6.54947e-04
            2.04606e-03   2.04606e-03   2.44922e-03
           -1.00391e-02  -1.00391e-02  -1.06531e-02
           -9.55977e-03  -9.55977e-03  -1.00809e-02
           -1.19952e-02  -1.19952e-02  -1.31019e-02
           -1.14461e-02  -1.14461e-02  -1.24258e-02
            1.06293e-05   1.06293e-05   7.21569e-05
            1.22848e-05   1.22848e-05   8.65397e-05
            1.85298e-05   1.85298e-05   1.51190e-04
            1.98236e-05   1.98236e-05   1.62052e-04
        """
        dens_cube_vals = np.array([float(x) for x in raw_vals.split()])
        dens_cube_vals = dens_cube_vals.reshape(nstates * 2, -1)

        self.run_nto(inpfile, xcfun_label, eig_vals, nto_lambdas, nto_cube_vals,
                     dens_cube_vals, 'rpa')


if __name__ == "__main__":
    unittest.main()
