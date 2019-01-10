from mpi4py import MPI
from veloxchem.veloxchemlib import OverlapMatrix
from veloxchem.veloxchemlib import KineticEnergyMatrix
from veloxchem.veloxchemlib import NuclearPotentialMatrix
from veloxchem.veloxchemlib import OverlapIntegralsDriver
from veloxchem.veloxchemlib import KineticEnergyIntegralsDriver
from veloxchem.veloxchemlib import NuclearPotentialIntegralsDriver
from veloxchem.veloxchemlib import SADGuessDriver
from veloxchem.veloxchemlib import ElectronRepulsionIntegralsDriver
from veloxchem.veloxchemlib import denmat
from veloxchem.veloxchemlib import fockmat
from veloxchem.veloxchemlib import ericut
from veloxchem.veloxchemlib import mpi_master

from veloxchem.mpitask import MpiTask
from veloxchem.aodensitymatrix import AODensityMatrix
from veloxchem.aofockmatrix import AOFockMatrix

import h5py
import numpy as np
import unittest


class TestPen(unittest.TestCase):

    def test_penicillin(self):

        task = MpiTask(["inputs/pen.inp", "inputs/pen.out"], MPI.COMM_WORLD)

        molecule = task.molecule
        ao_basis = task.ao_basis
        min_basis = task.min_basis

        comm = task.mpi_comm
        rank = task.mpi_rank
        size = task.mpi_size

        # compute 1e integrals

        ovldrv = OverlapIntegralsDriver(rank, size, comm)
        S = ovldrv.compute(molecule, ao_basis, comm)
        S1 = S.to_numpy()

        kindrv = KineticEnergyIntegralsDriver(rank, size, comm)
        T = kindrv.compute(molecule, ao_basis, comm)
        T1 = T.to_numpy()

        npotdrv = NuclearPotentialIntegralsDriver(rank, size, comm)
        V = npotdrv.compute(molecule, ao_basis, comm)
        V1 = V.to_numpy()

        # compare with reference

        hf = h5py.File("inputs/pen.onee.h5", 'r')
        S2 = np.array(hf.get("overlap"))
        T2 = np.array(hf.get("kinetic_energy"))
        V2 = np.array(hf.get("nuclear_potential"))
        hf.close()

        if rank == mpi_master():

            dS = np.max(np.abs(S1 - S2))
            dT = np.max(np.abs(T1 - T2))
            dV = np.max(np.abs(V1 - V2))

            self.assertTrue(dS < 1.0e-13)
            self.assertTrue(dT < 1.0e-11)
            self.assertTrue(dV < 1.0e-11)

        # read density

        if rank == mpi_master():
            dmat = AODensityMatrix.read_hdf5("inputs/pen.dens.h5")
        else:
            dmat = AODensityMatrix()
        dmat.broadcast(rank, comm)

        # compute Fock

        eridrv = ElectronRepulsionIntegralsDriver(rank, size, comm)

        qqdata = eridrv.compute(ericut.qq, 1.0e-12, molecule, ao_basis)

        fock = AOFockMatrix(dmat)

        eridrv.compute(fock, dmat, molecule, ao_basis, qqdata, comm)

        fock.reduce_sum(rank, size, comm)

        # compare with reference

        if rank == mpi_master():

            fock_ref = AOFockMatrix.read_hdf5("inputs/pen.twoe.h5")

            F1 = fock.to_numpy(0)
            F2 = fock_ref.to_numpy(0)

            maxelem = max(np.max(np.abs(F1)), np.max(np.abs(F2)))
            maxdiff = np.max(np.abs(F1 - F2))

            print(maxelem, maxdiff, maxdiff/maxelem)

            self.assertTrue(maxdiff / maxelem < 1.0e-11)


if __name__ == "__main__":
    unittest.main()
