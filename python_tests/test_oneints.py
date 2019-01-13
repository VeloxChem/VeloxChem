from mpi4py import MPI
from veloxchem.mpitask import MpiTask
from veloxchem.inputparser import InputParser
from veloxchem.veloxchemlib import Molecule
from veloxchem.veloxchemlib import OverlapMatrix
from veloxchem.veloxchemlib import KineticEnergyMatrix
from veloxchem.veloxchemlib import NuclearPotentialMatrix
from veloxchem.veloxchemlib import OverlapIntegralsDriver
from veloxchem.veloxchemlib import KineticEnergyIntegralsDriver
from veloxchem.veloxchemlib import NuclearPotentialIntegralsDriver
from veloxchem.veloxchemlib import mpi_master

import h5py
import math
import numpy as np
import unittest


class TestOneInts(unittest.TestCase):

    def test_overlap_matrix(self):

        array = np.array([[1.0, 0.2], [0.2, 1.0]])
        matrix = OverlapMatrix(array)
        array2 = matrix.to_numpy()
        matrix2 = OverlapMatrix(array2)

        self.assertTrue((array == array2).all())
        self.assertEqual(matrix, matrix2)

    def test_get_ortho_matrix(self):

        arr = np.array([[1.0, 0.2, 0.1], [0.2, 2.0, 0.3], [0.1, 0.3, 3.0]])

        mat = OverlapMatrix(arr)
        ortho_1 = mat.get_ortho_matrix(1.0e-12).to_numpy()

        evals, evecs = np.linalg.eigh(arr)
        evals_sqrt_inv = np.diag([1.0 / math.sqrt(x) for x in evals])
        ortho_2 = np.dot(evecs, np.dot(evals_sqrt_inv, evecs.T))

        diff = np.max(np.abs(ortho_1 - ortho_2))
        self.assertAlmostEqual(0., diff, 13)

    def test_kinetic_energy_matrix(self):

        array = np.array([[1.0, 0.2], [0.2, 1.0]])
        matrix = KineticEnergyMatrix(array)
        array2 = matrix.to_numpy()
        matrix2 = KineticEnergyMatrix(array2)

        self.assertTrue((array == array2).all())
        self.assertEqual(matrix, matrix2)

    def test_nuclear_potential_matrix(self):

        array = np.array([[1.0, 0.2], [0.2, 1.0]])
        matrix = NuclearPotentialMatrix(array)
        array2 = matrix.to_numpy()
        matrix2 = NuclearPotentialMatrix(array2)

        self.assertTrue((array == array2).all())
        self.assertEqual(matrix, matrix2)

    def test_1e_integrals(self):

        task = MpiTask(["inputs/h2se.inp", "inputs/h2se.out"], MPI.COMM_WORLD)

        molecule = task.molecule
        basis = task.ao_basis

        comm = task.mpi_comm
        rank = task.mpi_rank
        size = task.mpi_size

        # compute 1e integrals

        ovldrv = OverlapIntegralsDriver(rank, size, comm)
        S = ovldrv.compute(molecule, basis, comm)
        S1 = S.to_numpy()

        kindrv = KineticEnergyIntegralsDriver(rank, size, comm)
        T = kindrv.compute(molecule, basis, comm)
        T1 = T.to_numpy()

        npotdrv = NuclearPotentialIntegralsDriver(rank, size, comm)
        V = npotdrv.compute(molecule, basis, comm)
        V1 = V.to_numpy()

        # compare with reference

        if rank == mpi_master():

            hf = h5py.File("inputs/h2se.onee.h5", 'r')
            S2 = np.array(hf.get("overlap"))
            T2 = np.array(hf.get("kinetic_energy"))
            V2 = np.array(hf.get("nuclear_potential"))
            hf.close()

            dS = np.max(np.abs(S1 - S2))
            dT = np.max(np.abs(T1 - T2))
            dV = np.max(np.abs(V1 - V2))

            self.assertTrue(dS < 1.0e-13)
            self.assertTrue(dT < 1.0e-11)
            self.assertTrue(dV < 1.0e-11)

    def test_mixed_basis_1e(self):

        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        size = comm.Get_size()

        ovldrv = OverlapIntegralsDriver(rank, size, comm)
        kindrv = KineticEnergyIntegralsDriver(rank, size, comm)
        npotdrv = NuclearPotentialIntegralsDriver(rank, size, comm)

        bas_path = '../basis'

        # one molecule, one basis set

        mol_1 = InputParser.create_molecule_from_xyz('inputs/h2o.xyz')
        bas_1 = InputParser(bas_path + '/DEF2-SVP').create_basis_set(mol_1)

        S11 = ovldrv.compute(mol_1, bas_1, comm)
        T11 = kindrv.compute(mol_1, bas_1, comm)
        V11 = npotdrv.compute(mol_1, bas_1, comm)

        V11p = npotdrv.compute(mol_1, bas_1, mol_1, comm)
        self.assertTrue((V11.to_numpy() == V11p.to_numpy()).all())

        if rank == mpi_master():

            hf = h5py.File("inputs/mix_basis_1e.h5", 'r')
            ref_S11 = np.array(hf.get("S_h2o_def2-svp"))
            ref_T11 = np.array(hf.get("T_h2o_def2-svp"))
            ref_V11 = np.array(hf.get("V_h2o_def2-svp"))
            hf.close()

            dS = np.max(np.abs(S11.to_numpy() - ref_S11))
            dT = np.max(np.abs(T11.to_numpy() - ref_T11))
            dV = np.max(np.abs(V11.to_numpy() - ref_V11))

            self.assertTrue(dS < 1.0e-13)
            self.assertTrue(dT < 1.0e-13)
            self.assertTrue(dV < 1.0e-12)

        # one molecule, two basis sets

        mol_1 = InputParser.create_molecule_from_xyz('inputs/h2o.xyz')
        bas_1 = InputParser(bas_path + '/DEF2-SVP').create_basis_set(mol_1)
        bas_2 = InputParser(bas_path + '/CC-PVDZ').create_basis_set(mol_1)

        S12 = ovldrv.compute(mol_1, bas_1, bas_2, comm)
        T12 = kindrv.compute(mol_1, bas_1, bas_2, comm)
        V12 = npotdrv.compute(mol_1, bas_1, bas_2, mol_1, comm)

        if rank == mpi_master():

            hf = h5py.File("inputs/mix_basis_1e.h5", 'r')
            ref_S12 = np.array(hf.get("S_h2o_def2-svp_cc-pvdz"))
            ref_T12 = np.array(hf.get("T_h2o_def2-svp_cc-pvdz"))
            ref_V12 = np.array(hf.get("V_h2o_def2-svp_cc-pvdz"))
            hf.close()

            dS = np.max(np.abs(S12.to_numpy() - ref_S12))
            dT = np.max(np.abs(T12.to_numpy() - ref_T12))
            dV = np.max(np.abs(V12.to_numpy() - ref_V12))

            self.assertTrue(dS < 1.0e-13)
            self.assertTrue(dT < 1.0e-13)
            self.assertTrue(dV < 1.0e-12)

        # two molecules, one basis set

        mol_1 = InputParser.create_molecule_from_xyz('inputs/h2o.xyz')
        mol_2 = InputParser.create_molecule_from_xyz('inputs/nh3.xyz')
        mol = Molecule(mol_1, mol_2)
        bas = InputParser(bas_path + '/DEF2-SVP').create_basis_set(mol)

        S12 = ovldrv.compute(mol_1, mol_2, bas, comm)
        T12 = kindrv.compute(mol_1, mol_2, bas, comm)
        V12 = npotdrv.compute(mol_1, mol_2, bas, mol, comm)

        if rank == mpi_master():

            hf = h5py.File("inputs/mix_basis_1e.h5", 'r')
            ref_S12 = np.array(hf.get("S_h2o_nh3_def2-svp"))
            ref_T12 = np.array(hf.get("T_h2o_nh3_def2-svp"))
            ref_V12 = np.array(hf.get("V_h2o_nh3_def2-svp"))
            hf.close()

            dS = np.max(np.abs(S12.to_numpy() - ref_S12))
            dT = np.max(np.abs(T12.to_numpy() - ref_T12))
            dV = np.max(np.abs(V12.to_numpy() - ref_V12))

            self.assertTrue(dS < 1.0e-13)
            self.assertTrue(dT < 1.0e-11)
            self.assertTrue(dV < 1.0e-11)

        # two molecules, two basis sets

        mol_1 = InputParser.create_molecule_from_xyz('inputs/h2o.xyz')
        mol_2 = InputParser.create_molecule_from_xyz('inputs/nh3.xyz')
        mol = Molecule(mol_1, mol_2)
        bas_1 = InputParser(bas_path + '/DEF2-SVP').create_basis_set(mol_1)
        bas_2 = InputParser(bas_path + '/CC-PVDZ').create_basis_set(mol_2)

        S12 = ovldrv.compute(mol_1, mol_2, bas_1, bas_2, comm)
        T12 = kindrv.compute(mol_1, mol_2, bas_1, bas_2, comm)
        V12 = npotdrv.compute(mol_1, mol_2, bas_1, bas_2, mol, comm)

        if rank == mpi_master():

            hf = h5py.File("inputs/mix_basis_1e.h5", 'r')
            ref_S12 = np.array(hf.get("S_h2o_def2-svp_nh3_cc-pvdz"))
            ref_T12 = np.array(hf.get("T_h2o_def2-svp_nh3_cc-pvdz"))
            ref_V12 = np.array(hf.get("V_h2o_def2-svp_nh3_cc-pvdz"))
            hf.close()

            dS = np.max(np.abs(S12.to_numpy() - ref_S12))
            dT = np.max(np.abs(T12.to_numpy() - ref_T12))
            dV = np.max(np.abs(V12.to_numpy() - ref_V12))

            self.assertTrue(dS < 1.0e-13)
            self.assertTrue(dT < 1.0e-11)
            self.assertTrue(dV < 1.0e-11)


if __name__ == "__main__":
    unittest.main()
