from mpi4py import MPI
from HelperClass import Task
from VeloxChemLib import Molecule
from VeloxChemLib import MolecularBasis
from VeloxChemLib import OverlapIntegralsDriver
from VeloxChemLib import KineticEnergyIntegralsDriver
from VeloxChemLib import NuclearPotentialIntegralsDriver
from VeloxChemLib import mpi_master
from VeloxChemLib import assemble_overlap_matrices
from VeloxChemLib import assemble_kinetic_energy_matrices
from VeloxChemLib import assemble_nuclear_potential_matrices

import numpy as np
import unittest


class TestExciton(unittest.TestCase):

    def test_assemble_matrices(self):

        # set up MPI

        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        size = comm.Get_size()

        # process input

        if (rank == mpi_master()):

            task = Task("inputs/dimer.inp", "inputs/dimer.out")
            molecule = task.molecule
            basis = task.ao_basis
            ostream = task.ostream

        else:

            molecule = Molecule()
            basis = MolecularBasis()
            ostream = OutputStream("")

        # broadcast molecule and basis

        molecule.broadcast(rank, comm)
        basis.broadcast(rank, comm)

        # build sub molecules

        mol_1 = molecule.get_sub_molecule(0, 4)
        mol_2 = molecule.get_sub_molecule(4, 5)

        # compute overlap

        ovldrv = OverlapIntegralsDriver.create(rank, size, comm)
        S = ovldrv.compute(molecule, basis, ostream, comm)
        S11 = ovldrv.compute(mol_1, basis, ostream, comm)
        S22 = ovldrv.compute(mol_2, basis, ostream, comm)
        S12 = ovldrv.compute(mol_1, mol_2, basis, ostream, comm)
        S21 = ovldrv.compute(mol_2, mol_1, basis, ostream, comm)

        if (rank == mpi_master()):

            S_exmod = assemble_overlap_matrices(
                mol_1, mol_2, basis, basis, S11, S22, S12, S21)

            self.assertEqual(S, S_exmod)

        # compute kinetic energy

        kindrv = KineticEnergyIntegralsDriver.create(rank, size, comm)
        T = kindrv.compute(molecule, basis, ostream, comm)
        T11 = kindrv.compute(mol_1, basis, ostream, comm)
        T22 = kindrv.compute(mol_2, basis, ostream, comm)
        T12 = kindrv.compute(mol_1, mol_2, basis, ostream, comm)
        T21 = kindrv.compute(mol_2, mol_1, basis, ostream, comm)

        if (rank == mpi_master()):

            T_exmod = assemble_kinetic_energy_matrices(
                mol_1, mol_2, basis, basis, T11, T22, T12, T21)

            self.assertEqual(T, T_exmod)

        # compute nuclear potential

        npotdrv = NuclearPotentialIntegralsDriver.create(rank, size, comm)
        V = npotdrv.compute(molecule, basis, ostream, comm)
        V11 = npotdrv.compute(mol_1, basis, molecule, ostream, comm)
        V22 = npotdrv.compute(mol_2, basis, molecule, ostream, comm)
        V12 = npotdrv.compute(mol_1, mol_2, basis, molecule, ostream, comm)
        V21 = npotdrv.compute(mol_2, mol_1, basis, molecule, ostream, comm)

        if (rank == mpi_master()):

            V_exmod = assemble_nuclear_potential_matrices(
                mol_1, mol_2, basis, basis, V11, V22, V12, V21)

            self.assertEqual(V, V_exmod)


if __name__ == "__main__":
    unittest.main()
