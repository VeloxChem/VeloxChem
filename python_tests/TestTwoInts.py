from mpi4py import MPI
from HelperClass import Task
from veloxchem.VeloxChemLib import ElectronRepulsionIntegralsDriver
from veloxchem.VeloxChemLib import OverlapIntegralsDriver
from veloxchem.VeloxChemLib import SADGuessDriver
from veloxchem.VeloxChemLib import AODensityMatrix
from veloxchem.VeloxChemLib import AOFockMatrix
from veloxchem.VeloxChemLib import denmat
from veloxchem.VeloxChemLib import fockmat
from veloxchem.VeloxChemLib import ericut
from veloxchem.VeloxChemLib import mpi_master

import numpy as np
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

        fock = AOFockMatrix.from_numpy_list(
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

    def test_fock_rest_density(self):

        data_a = [[ 1., .2, ], [ .2, 1., ]]

        d_rest = AODensityMatrix.from_numpy_list([data_a], denmat.rest)

        f_rest = AOFockMatrix(d_rest)

        self.assertEqual(fockmat.restjk, f_rest.get_fock_type(0))
        self.assertEqual(1.0, f_rest.get_scale_factor(0))
        self.assertEqual(0, f_rest.get_density_identifier(0))

    def test_fock_rest_build(self):

        # mpi settings

        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        size = comm.Get_size()

        # process input file on master node

        if (rank == mpi_master()):

            task = Task("inputs/water.inp", "inputs/water.out")
            molecule = task.molecule
            ao_basis = task.ao_basis
            min_basis = task.min_basis
            ostream = task.ostream

        else:

            molecule = Molecule()
            ao_basis = MolecularBasis()
            min_basis = MolecularBasis()
            ostream = OutputStream("")

        # broadcast molecule and basis

        molecule.broadcast(rank, comm)
        ao_basis.broadcast(rank, comm)
        min_basis.broadcast(rank, comm)

        # compute overlap

        ovldrv = OverlapIntegralsDriver.create(rank, size, comm)
        S12 = ovldrv.compute(molecule, min_basis, ao_basis, ostream, comm)
        S22 = ovldrv.compute(molecule, ao_basis, ostream, comm)

        # compute initial guess

        saddrv = SADGuessDriver.create(rank, size, comm)
        dsad = saddrv.compute(molecule, min_basis, ao_basis, S12, S22, ostream,
                              comm)

        # compute Fock

        eridrv = ElectronRepulsionIntegralsDriver.create(rank, size, comm)

        qqdata = eridrv.compute(ericut.qq, 1.0e-12, molecule, ao_basis, ostream,
                                comm)

        fock = AOFockMatrix(dsad)

        fock.zero()

        eridrv.compute(fock, dsad, molecule, ao_basis, qqdata, ostream, comm)

        self.assertEqual(fockmat.restjk, fock.get_fock_type(0))
        self.assertEqual(1.0, fock.get_scale_factor(0))
        self.assertEqual(0, fock.get_density_identifier(0))


if __name__ == "__main__":
    unittest.main()
