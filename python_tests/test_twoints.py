from pathlib import Path
import numpy as np
import unittest
import math

from veloxchem.veloxchemlib import KineticEnergyMatrix
from veloxchem.veloxchemlib import NuclearPotentialMatrix
from veloxchem.veloxchemlib import ElectronRepulsionIntegralsDriver
from veloxchem.veloxchemlib import denmat
from veloxchem.veloxchemlib import fockmat
from veloxchem.veloxchemlib import ericut
from veloxchem.veloxchemlib import is_mpi_master
from veloxchem.mpitask import MpiTask
from veloxchem.aodensitymatrix import AODensityMatrix
from veloxchem.aofockmatrix import AOFockMatrix
from veloxchem.qqscheme import get_qq_scheme


class TestTwoInts(unittest.TestCase):

    def test_fock_matrix(self):

        data_j = [[1., .2], [.2, 1.]]
        data_k = [[.9, .5], [.5, .9]]

        arr_j = np.array(data_j)
        arr_k = np.array(data_k)
        arr_jk = arr_j + arr_k

        x = -0.5
        arr_kx = x * arr_k
        arr_jkx = arr_j + arr_kx

        fock = AOFockMatrix(
            [arr_jk, arr_jkx, arr_j, arr_k, arr_kx],
            [
                fockmat.restjk, fockmat.restjkx, fockmat.restj, fockmat.restk,
                fockmat.restkx
            ],
            [1.0, x, 1.0, 1.0, x],
            [0, 0, 0, 0, 0],
        )

        np_jkx = fock.to_numpy(1)
        np_j = fock.to_numpy(2)
        np_k = fock.to_numpy(3)

        self.assertTrue((arr_j == np_j).all())
        self.assertTrue((arr_k == np_k).all())
        self.assertTrue((arr_jkx == np_jkx).all())

        self.assertEqual(fockmat.restjk, fock.get_fock_type(0))
        self.assertEqual(x, fock.get_scale_factor(1))
        self.assertEqual(0, fock.get_density_identifier(2))

        self.assertEqual(5, fock.number_of_fock_matrices())

    def test_unrestricted(self):

        data_a = [[1., .2], [.2, 1.]]
        data_b = [[.9, .5], [.5, .9]]

        arr_a = np.array(data_a)
        arr_b = np.array(data_b)

        fock = AOFockMatrix(
            [arr_a, arr_b],
            [fockmat.unrestjk, fockmat.unrestjk],
            [0.7, 0.6],
            [1, 2],
        )

        self.assertEqual(1, fock.number_of_fock_matrices())

        self.assertTrue((fock.alpha_to_numpy(0) == arr_a).all())
        self.assertTrue((fock.beta_to_numpy(0) == arr_b).all())

        self.assertEqual(fockmat.unrestjk, fock.get_fock_type(0))
        self.assertEqual(fockmat.unrestjk, fock.get_fock_type(0, True))

        self.assertEqual(0.7, fock.get_scale_factor(0))
        self.assertEqual(0.6, fock.get_scale_factor(0, True))

        self.assertEqual(1, fock.get_density_identifier(0))
        self.assertEqual(2, fock.get_density_identifier(0, True))

    def test_add_hcore(self):

        arr_t = np.array([[3., .2], [.2, 3.]])
        arr_v = np.array([[-9., .5], [.5, -9.]])
        arr_jk = np.array([[5., .1], [.1, 5.]])
        arr_fock = arr_jk + arr_t - arr_v

        kin = KineticEnergyMatrix(arr_t)
        npot = NuclearPotentialMatrix(arr_v)
        fock = AOFockMatrix([arr_jk], [fockmat.restjk], [1.0], [0])
        fock.add_hcore(kin, npot, 0)

        diff = np.max(np.abs(fock.to_numpy(0) - arr_fock))
        self.assertAlmostEqual(0., diff, 13)

    def test_add_fock(self):

        arr_1 = np.array([[1., .2], [.2, 1.]])
        arr_2 = np.array([[.9, .5], [.5, .9]])

        fock_1 = AOFockMatrix([arr_1], [fockmat.restjk], [1.0], [0])
        fock_2 = AOFockMatrix([arr_2], [fockmat.restjk], [1.0], [0])

        fock_sum = AOFockMatrix(fock_1)
        fock_sum.add(fock_2)

        diff = np.max(np.abs(fock_sum.to_numpy(0) - (arr_1 + arr_2)))
        self.assertAlmostEqual(0., diff, 13)

    def test_fock_density(self):

        data_a = [[1., .2], [.2, 1.]]

        d_rest = AODensityMatrix([data_a], denmat.rest)

        f_rest = AOFockMatrix(d_rest)

        self.assertEqual(fockmat.restjk, f_rest.get_fock_type(0))
        self.assertEqual(1.0, f_rest.get_scale_factor(0))
        self.assertEqual(0, f_rest.get_density_identifier(0))

    def test_fock_hdf5(self):

        data_a = [[1., .2], [.2, 1.]]
        data_b = [[.9, .5], [.5, .9]]
        data_c = [[.8, .6], [.6, .8]]
        data_d = [[.7, .5], [.5, .7]]

        types = [fockmat.restk, fockmat.restjkx, fockmat.restk, fockmat.restjk]

        factors = [1.0, 0.2, 1.0, 1.0]

        indices = [0, 0, 1, 0]

        f_rest = AOFockMatrix([data_a, data_b, data_c, data_d], types, factors,
                              indices)

        # hdf5 read/write tests

        if is_mpi_master():

            here = Path(__file__).parent
            h5file = here / 'inputs' / 'dummy.h5'

            f_rest.write_hdf5(h5file)
            f2 = AOFockMatrix.read_hdf5(h5file)
            self.assertEqual(f_rest, f2)

    def test_fock_build(self):

        here = Path(__file__).parent
        inpfile = here / 'inputs' / 'h2se.inp'
        outfile = inpfile.with_suffix('.out')

        task = MpiTask([str(inpfile), str(outfile)])

        molecule = task.molecule
        ao_basis = task.ao_basis

        molecule.check_proximity(0.1)

        molecule.check_multiplicity()

        enuc = molecule.nuclear_repulsion_energy()

        ref_enuc = 34.0 / 2.8 + 34.0 / 2.8 + 1.0 / (2.8 * math.sqrt(2.0))

        self.assertAlmostEqual(enuc, ref_enuc, 13)

        comm = task.mpi_comm
        rank = task.mpi_rank
        size = task.mpi_size

        # read density

        if is_mpi_master(comm):
            densfile = str(here / 'inputs' / 'h2se.dens.h5')

            dmat = AODensityMatrix.read_hdf5(densfile)
        else:
            dmat = AODensityMatrix()

        dmat.broadcast(rank, comm)

        # compute Fock

        eridrv = ElectronRepulsionIntegralsDriver(comm)

        qqdata = eridrv.compute(ericut.qqden, 1.0e-12, molecule, ao_basis)

        self.assertEqual(get_qq_scheme("QQ"), ericut.qq)
        self.assertEqual(get_qq_scheme("QQR"), ericut.qqr)
        self.assertEqual(get_qq_scheme("QQ_DEN"), ericut.qqden)
        self.assertEqual(get_qq_scheme("QQR_DEN"), ericut.qqrden)

        num_screeners = qqdata.number_of_screeners()
        self.assertTrue(num_screeners > 0)
        self.assertFalse(qqdata.is_empty())

        for screener_index in range(num_screeners):
            screener = qqdata.get_screener(screener_index)
            self.assertEqual(1.0e-12, screener.get_threshold())
            self.assertEqual(ericut.qqden, screener.get_screening_scheme())

        fock = AOFockMatrix(dmat)

        eridrv.compute(fock, dmat, molecule, ao_basis, qqdata)

        fock.reduce_sum(rank, size, comm)

        # compare with unrestricted Fock

        if is_mpi_master(comm):
            da = dmat.alpha_to_numpy(0)
            db = dmat.beta_to_numpy(0)
            d2 = AODensityMatrix([da, db], denmat.unrest)
        else:
            d2 = AODensityMatrix()
        d2.broadcast(rank, comm)

        f2 = AOFockMatrix(d2)
        eridrv.compute(f2, d2, molecule, ao_basis, qqdata)
        f2.reduce_sum(rank, size, comm)

        self.assertEqual(1, f2.number_of_fock_matrices())
        self.assertEqual(fockmat.unrestjk, f2.get_fock_type(0))
        self.assertEqual(fockmat.unrestjk, f2.get_fock_type(0, True))

        self.assertTrue(
            np.max(np.abs(fock.alpha_to_numpy(0) -
                          f2.alpha_to_numpy(0))) < 1.0e-13)
        self.assertTrue(
            np.max(np.abs(fock.beta_to_numpy(0) -
                          f2.beta_to_numpy(0))) < 1.0e-13)

        # compare with reference

        if is_mpi_master(comm):

            twoefile = str(here / 'inputs' / 'h2se.twoe.h5')

            fock_ref = AOFockMatrix.read_hdf5(twoefile)

            F1 = fock.to_numpy(0)
            F2 = fock_ref.to_numpy(0)
            dF = np.max(np.abs(F1 - F2))

            self.assertTrue(dF < 1.0e-11)


if __name__ == "__main__":
    unittest.main()
