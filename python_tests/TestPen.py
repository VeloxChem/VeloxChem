from mpi4py import MPI
from HelperClass import Task
from veloxchem.VeloxChemLib import OverlapMatrix
from veloxchem.VeloxChemLib import KineticEnergyMatrix
from veloxchem.VeloxChemLib import NuclearPotentialMatrix
from veloxchem.VeloxChemLib import OverlapIntegralsDriver
from veloxchem.VeloxChemLib import KineticEnergyIntegralsDriver
from veloxchem.VeloxChemLib import NuclearPotentialIntegralsDriver
from veloxchem.VeloxChemLib import SADGuessDriver
from veloxchem.VeloxChemLib import ElectronRepulsionIntegralsDriver
from veloxchem.VeloxChemLib import ScreeningContainer
from veloxchem.VeloxChemLib import denmat
from veloxchem.VeloxChemLib import fockmat
from veloxchem.VeloxChemLib import ericut
from veloxchem.VeloxChemLib import mpi_master

from veloxchem.aodensitymatrix import AODensityMatrix
from veloxchem.aofockmatrix import AOFockMatrix

import h5py
import numpy as np
import unittest


class TestPen(unittest.TestCase):

    def test_penicillin(self):

        # mpi settings

        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        size = comm.Get_size()

        # process input file on master node

        if (rank == mpi_master()):

            task = Task("inputs/pen.inp", "inputs/pen.out")
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

        # compute 1e integrals

        ovldrv = OverlapIntegralsDriver.create(rank, size, comm)
        S = ovldrv.compute(molecule, ao_basis, ostream, comm)
        S1 = S.to_numpy()

        kindrv = KineticEnergyIntegralsDriver.create(rank, size, comm)
        T = kindrv.compute(molecule, ao_basis, ostream, comm)
        T1 = T.to_numpy()

        npotdrv = NuclearPotentialIntegralsDriver.create(rank, size, comm)
        V = npotdrv.compute(molecule, ao_basis, ostream, comm)
        V1 = V.to_numpy()

        # compare with reference

        hf = h5py.File("inputs/pen.onee.h5", 'r')
        S2 = np.array(hf.get("overlap"))
        T2 = np.array(hf.get("kinetic_energy"))
        V2 = np.array(hf.get("nuclear_potential"))
        hf.close()

        S_diff = np.max(np.abs(S1 - S2))
        T_diff = np.max(np.abs(T1 - T2))
        V_diff = np.max(np.abs(V1 - V2))

        self.assertTrue(S_diff < 1.0e-13)
        self.assertTrue(T_diff < 1.0e-11)
        self.assertTrue(V_diff < 1.0e-11)

        # read density

        dmat = AODensityMatrix.read_hdf5("inputs/pen.dens.h5")

        # compute Fock

        eridrv = ElectronRepulsionIntegralsDriver.create(rank, size, comm)

        """
        qqdata = eridrv.compute(ericut.qq, 1.0e-12, molecule, ao_basis, ostream,
                                comm)
        """

        qqdata = ScreeningContainer()

        fock = AOFockMatrix(dmat)

        eridrv.compute(fock, dmat, molecule, ao_basis, qqdata, ostream, comm)

        F1 = fock.to_numpy(0)

        # compare with reference

        F2 = AOFockMatrix.read_hdf5("inputs/pen.twoe.h5").to_numpy(0)

        maxelem = max(np.max(np.abs(F1)), np.max(np.abs(F2)))

        maxdiff = np.max(np.abs(F1 - F2))

        self.assertTrue(maxdiff / maxelem < 1.0e-11)


if __name__ == "__main__":
    unittest.main()
