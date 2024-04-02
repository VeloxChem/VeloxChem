from pathlib import Path
import numpy as np

from veloxchem.veloxchemlib import is_mpi_master, mpi_barrier
from veloxchem.mpitask import MpiTask
from veloxchem.cubicgrid import CubicGrid
from veloxchem.molecularorbitals import MolecularOrbitals
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.tdaeigensolver import TdaEigenSolver
from veloxchem.lreigensolver import LinearResponseEigenSolver


class TestNTO:

    def run_nto(self, inpfile, xcfun_label, ref_eig_vals, ref_nto_lambdas,
                ref_nto_cube_vals, ref_dens_cube_vals, flag):

        task = MpiTask([str(inpfile), None])
        task.input_dict['scf']['checkpoint_file'] = None

        if xcfun_label is not None:
            task.input_dict['method_settings']['xcfun'] = xcfun_label

        # run SCF

        scf_drv = ScfRestrictedDriver(task.mpi_comm, task.ostream)
        scf_drv.update_settings(task.input_dict['scf'],
                                task.input_dict['method_settings'])
        scf_results = scf_drv.compute(task.molecule, task.ao_basis,
                                      task.min_basis)

        # run TDA/RPA

        rsp_dict = {
            'filename': str(inpfile.parent / inpfile.stem),
            'nstates': ref_eig_vals.shape[0],
            'nto': 'yes',
            'detach_attach': 'yes',
            'cube_origin': '-0.1, -0.1, -0.1',
            'cube_stepsize': '0.3, 0.3, 0.2',
            'cube_points': '2, 2, 3',
        }

        if flag == 'tda':
            rsp_drv = TdaEigenSolver(task.mpi_comm, task.ostream)
        elif flag == 'rpa':
            rsp_drv = LinearResponseEigenSolver(task.mpi_comm, task.ostream)
        rsp_drv.update_settings(rsp_dict, task.input_dict['method_settings'])
        rsp_results = rsp_drv.compute(task.molecule, task.ao_basis, scf_results)

        if is_mpi_master(task.mpi_comm):
            eig_vals = rsp_results['eigenvalues']

            nto_lambdas = []
            nto_cube_vals = []
            dens_cube_vals = []

            for s in range(ref_eig_vals.shape[0]):
                lam_diag = rsp_results['nto_lambdas'][s]
                nto_lambdas.append(lam_diag[0])

                nto_h5_fname = rsp_results['nto_h5_files'][s]
                nto_cube_fnames = rsp_results['nto_cubes'][s]
                assert Path(
                    nto_cube_fnames[0]).stem == Path(nto_h5_fname).stem + '_H1'
                assert Path(
                    nto_cube_fnames[1]).stem == Path(nto_h5_fname).stem + '_P1'

                nto_mo = MolecularOrbitals.read_hdf5(nto_h5_fname)
                assert nto_mo.is_nto()

                for fname in nto_cube_fnames[:2]:
                    read_grid = CubicGrid.read_cube(fname)
                    nto_cube_vals.append(read_grid.values_to_numpy().flatten())

                dens_cube_fnames = rsp_results['density_cubes'][s]
                for fname in dens_cube_fnames:
                    read_grid = CubicGrid.read_cube(fname)
                    dens_cube_vals.append(read_grid.values_to_numpy().flatten())

        # compare with reference

        if is_mpi_master(task.mpi_comm):
            nto_lambdas = np.array(nto_lambdas)
            nto_cube_vals = np.array(nto_cube_vals)
            dens_cube_vals = np.array(dens_cube_vals)

            thresh = 1.0e-6 if xcfun_label is None else 1.0e-5
            assert np.max(np.abs(eig_vals - ref_eig_vals)) < thresh

            thresh = 1.0e-4 if xcfun_label is None else 1.0e-3
            assert np.max(np.abs(nto_lambdas - ref_nto_lambdas)) < thresh

            for s in range(nto_cube_vals.shape[0]):
                if np.vdot(nto_cube_vals[s, :], ref_nto_cube_vals[s, :]) < 0.0:
                    nto_cube_vals[s, :] *= -1.0
                if np.vdot(dens_cube_vals[s, :],
                           ref_dens_cube_vals[s, :]) < 0.0:
                    dens_cube_vals[s, :] *= -1.0
            assert np.max(np.abs(nto_cube_vals - ref_nto_cube_vals)) < 1.0e-4
            assert np.max(np.abs(dens_cube_vals - ref_dens_cube_vals)) < 1.0e-4

            # clean up

            scf_h5 = inpfile.with_suffix('.scf.h5')
            if scf_h5.is_file():
                scf_h5.unlink()

            scf_final_h5 = scf_h5.with_suffix('.results.h5')
            if scf_final_h5.is_file():
                scf_final_h5.unlink()

            rsp_h5 = inpfile.with_suffix('.rsp.h5')
            if rsp_h5.is_file():
                rsp_h5.unlink()

            rsp_solutions_h5 = rsp_h5.with_suffix('.solutions.h5')
            if rsp_solutions_h5.is_file():
                rsp_solutions_h5.unlink()

            for s in range(ref_eig_vals.shape[0]):
                nto_h5_fname = rsp_results['nto_h5_files'][s]
                nto_h5 = Path(nto_h5_fname)
                if nto_h5.is_file():
                    nto_h5.unlink()

                nto_cube_fnames = rsp_results['nto_cubes'][s]
                dens_cube_fnames = rsp_results['density_cubes'][s]
                for fname in (nto_cube_fnames + dens_cube_fnames):
                    fpath = Path(fname)
                    if fpath.is_file():
                        fpath.unlink()

        mpi_barrier()

    def test_nto_tda(self):

        here = Path(__file__).parent
        inpfile = here / 'inputs' / 'c2h4.inp'

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
        inpfile = here / 'inputs' / 'c2h4.inp'

        xcfun_label = 'b3lyp'
        nstates = 5

        raw_vals = """
            0.30465317  0.95542402
            0.30537184  0.72744400
            0.31770840  0.98575971
            0.33176863  0.97458205
            0.35264869  0.88407337
        """
        vals = np.array([float(x) for x in raw_vals.split()])
        vals = vals.reshape(nstates, -1)

        eig_vals = vals[:, 0]
        nto_lambdas = vals[:, 1]

        raw_vals = """
           -3.21430e-03   3.21430e-03   1.03578e-02
            6.22010e-03  -6.22010e-03  -2.00043e-02
           -3.10657e-03   3.10657e-03   9.99145e-03
            6.01501e-03  -6.01501e-03  -1.93083e-02
           -4.77930e-03   4.77930e-03   1.52684e-02
           -4.62358e-03   4.62358e-03   1.47413e-02
            9.24647e-03  -9.24647e-03  -2.94803e-02
            8.95013e-03  -8.95013e-03  -2.84793e-02
           -3.05252e-02  -3.05252e-02  -3.41449e-02
           -2.99013e-02  -2.99013e-02  -3.33239e-02
            5.98010e-02   5.98010e-02   6.66455e-02
            5.85944e-02   5.85944e-02   6.50683e-02
           -4.07426e-03   4.07426e-03   1.29847e-02
           -3.94665e-03   3.94665e-03   1.25527e-02
            7.89235e-03  -7.89235e-03  -2.51022e-02
            7.64945e-03  -7.64945e-03  -2.42819e-02
            3.05252e-02   3.05252e-02   3.41449e-02
            2.99013e-02   2.99013e-02   3.33239e-02
           -5.98010e-02  -5.98010e-02  -6.66455e-02
           -5.85944e-02  -5.85944e-02  -6.50683e-02
            1.64671e-02   1.64671e-02   2.15365e-02
            1.57028e-02   1.57028e-02   2.01180e-02
            1.50995e-02   1.50995e-02   1.95102e-02
            1.44099e-02   1.44099e-02   1.82198e-02
           -3.05252e-02  -3.05252e-02  -3.41449e-02
           -2.99013e-02  -2.99013e-02  -3.33239e-02
            5.98010e-02   5.98010e-02   6.66455e-02
            5.85944e-02   5.85944e-02   6.50683e-02
           -2.36979e-02  -2.36979e-02  -2.60558e-02
            4.63834e-02   4.63834e-02   5.08241e-02
           -2.32036e-02  -2.32036e-02  -2.54242e-02
            4.54265e-02   4.54265e-02   4.96091e-02
           -3.05252e-02  -3.05252e-02  -3.41449e-02
           -2.99013e-02  -2.99013e-02  -3.33239e-02
            5.98010e-02   5.98010e-02   6.66455e-02
            5.85944e-02   5.85944e-02   6.50683e-02
           -2.95965e-03   2.95965e-03   7.09633e-03
           -3.27570e-03   3.27570e-03   8.22654e-03
           -3.27604e-03   3.27604e-03   8.22798e-03
           -3.56641e-03   3.56641e-03   9.26185e-03
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
