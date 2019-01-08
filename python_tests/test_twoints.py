from mpi4py import MPI
from veloxchem.mpitask import MpiTask
from veloxchem.veloxchemlib import MolecularBasis
from veloxchem.veloxchemlib import OverlapIntegralsDriver
from veloxchem.veloxchemlib import SADGuessDriver
from veloxchem.veloxchemlib import ElectronRepulsionIntegralsDriver
from veloxchem.veloxchemlib import ScreeningContainer
from veloxchem.veloxchemlib import denmat
from veloxchem.veloxchemlib import fockmat
from veloxchem.veloxchemlib import ericut
from veloxchem.veloxchemlib import mpi_master

from veloxchem.aodensitymatrix import AODensityMatrix
from veloxchem.aofockmatrix import AOFockMatrix

import numpy as np
import math
import unittest


class TestTwoInts(unittest.TestCase):

    def test_fock_matrix(self):

        data_j = [[ 1., .2, ], [ .2, 1., ]]
        data_k = [[ .9, .5, ], [ .5, .9, ]]

        arr_j = np.array(data_j)
        arr_k = np.array(data_k)
        arr_jk = arr_j + arr_k

        x = -0.5
        arr_kx = x * arr_k
        arr_jkx = arr_j + arr_kx

        fock = AOFockMatrix(
            [arr_jk, arr_jkx, arr_j, arr_k, arr_kx],
            [fockmat.restjk, fockmat.restjkx,
                fockmat.restj, fockmat.restk, fockmat.restkx],
            [1.0, x, 1.0, 1.0, x],
            [0, 0, 0, 0, 0])

        np_jkx = fock.to_numpy(1)
        np_j = fock.to_numpy(2)
        np_k = fock.to_numpy(3)

        self.assertEqual(0, np.max(np.abs(arr_j - np_j)))
        self.assertEqual(0, np.max(np.abs(arr_k - np_k)))
        self.assertEqual(0, np.max(np.abs(arr_jkx - np_jkx)))

        self.assertEqual(fockmat.restjk, fock.get_fock_type(0))
        self.assertEqual(x, fock.get_scale_factor(1))
        self.assertEqual(0, fock.get_density_identifier(2))

    def test_fock_density(self):

        data_a = [[ 1., .2, ], [ .2, 1., ]]

        d_rest = AODensityMatrix([data_a], denmat.rest)

        f_rest = AOFockMatrix(d_rest)

        self.assertEqual(fockmat.restjk, f_rest.get_fock_type(0))
        self.assertEqual(1.0, f_rest.get_scale_factor(0))
        self.assertEqual(0, f_rest.get_density_identifier(0))

    def test_fock_hdf5(self):

        data_a = [[ 1., .2, ], [ .2, 1., ]]
        data_b = [[ .9, .5, ], [ .5, .9, ]]
        data_c = [[ .8, .6, ], [ .6, .8, ]]
        data_d = [[ .7, .5, ], [ .5, .7, ]]

        types = [fockmat.restk, fockmat.restjkx, fockmat.restk, fockmat.restjk]

        factors = [1.0, 0.2, 1.0, 1.0]

        indices = [0, 0, 1, 0]

        f_rest = AOFockMatrix([data_a, data_b, data_c, data_d],
                              types, factors, indices)

        # hdf5 read/write tests

        if MPI.COMM_WORLD.Get_rank() == mpi_master():

            f_rest.write_hdf5("inputs/dummy.h5")

            f2 = AOFockMatrix.read_hdf5("inputs/dummy.h5")

            self.assertEqual(f_rest, f2)

    def test_fock_build(self):

        task = MpiTask(["inputs/h2se.inp", "inputs/h2se.out"], MPI.COMM_WORLD)

        molecule = task.molecule
        ao_basis = task.ao_basis
        min_basis = task.min_basis

        molecule.check_proximity(0.1)

        molecule.check_multiplicity()

        enuc = molecule.nuclear_repulsion_energy()

        ref_enuc = 34.0 / 2.8 + 34.0 / 2.8 + 1.0 / (2.8 * math.sqrt(2.0))

        self.assertEqual(enuc, ref_enuc)

        comm = task.mpi_comm
        rank = task.mpi_rank
        size = task.mpi_size

        # read density

        dmat = AODensityMatrix.read_hdf5("inputs/h2se.dens.h5")

        # compute Fock

        eridrv = ElectronRepulsionIntegralsDriver(rank, size, comm)

        qqdata = eridrv.compute(ericut.qq, 1.0e-12, molecule, ao_basis)

        fock = AOFockMatrix(dmat)

        eridrv.compute(fock, dmat, molecule, ao_basis, qqdata, comm)

        F1 = fock.to_numpy(0)

        # compare with reference

        fock_ref = AOFockMatrix.read_hdf5("inputs/h2se.twoe.h5")

        F2 = fock_ref.to_numpy(0)

        if rank == mpi_master():

            dF = np.max(np.abs(F1 - F2))
            self.assertTrue(dF < 1.0e-11)


if __name__ == "__main__":
    unittest.main()
