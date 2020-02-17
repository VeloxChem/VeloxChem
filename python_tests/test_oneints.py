from mpi4py import MPI
import numpy as np
import unittest
import h5py
import os

from veloxchem.veloxchemlib import OverlapMatrix
from veloxchem.veloxchemlib import KineticEnergyMatrix
from veloxchem.veloxchemlib import NuclearPotentialMatrix
from veloxchem.veloxchemlib import OverlapIntegralsDriver
from veloxchem.veloxchemlib import KineticEnergyIntegralsDriver
from veloxchem.veloxchemlib import NuclearPotentialIntegralsDriver
from veloxchem.veloxchemlib import mpi_master
from veloxchem.mpitask import MpiTask
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis


class TestOneInts(unittest.TestCase):

    def test_overlap_matrix(self):

        array = np.array([[1.0, 0.2], [0.2, 1.0]])
        matrix = OverlapMatrix(array)
        array2 = matrix.to_numpy()
        matrix2 = OverlapMatrix(array2)

        self.assertTrue((array == array2).all())
        self.assertEqual(matrix, matrix2)

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

        inpfile = os.path.join('inputs', 'h2se.inp')
        if not os.path.isfile(inpfile):
            inpfile = os.path.join('python_tests', inpfile)
        outfile = inpfile.replace('.inp', '.out')

        task = MpiTask([inpfile, outfile], MPI.COMM_WORLD)
        breakpoint()

        molecule = task.molecule
        basis = task.ao_basis

        comm = task.mpi_comm
        rank = task.mpi_rank

        # compute 1e integrals

        ovldrv = OverlapIntegralsDriver(comm)
        S = ovldrv.compute(molecule, basis)
        S1 = S.to_numpy()

        kindrv = KineticEnergyIntegralsDriver(comm)
        T = kindrv.compute(molecule, basis)
        T1 = T.to_numpy()

        npotdrv = NuclearPotentialIntegralsDriver(comm)
        V = npotdrv.compute(molecule, basis)
        V1 = V.to_numpy()

        # compare with reference

        if rank == mpi_master():

            h5file = os.path.join('inputs', 'h2se.onee.h5')
            if not os.path.isfile(h5file):
                h5file = os.path.join('python_tests', h5file)

            hf = h5py.File(h5file, 'r')
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

        ovldrv = OverlapIntegralsDriver(comm)
        kindrv = KineticEnergyIntegralsDriver(comm)
        npotdrv = NuclearPotentialIntegralsDriver(comm)

        h2ofile = os.path.join('inputs', 'h2o.xyz')
        if not os.path.isfile(h2ofile):
            h2ofile = os.path.join('python_tests', h2ofile)

        nh3file = os.path.join('inputs', 'nh3.xyz')
        if not os.path.isfile(nh3file):
            nh3file = os.path.join('python_tests', nh3file)

        h5file = os.path.join('inputs', 'mix_basis_1e.h5')
        if not os.path.isfile(h5file):
            h5file = os.path.join('python_tests', h5file)

        # one molecule, one basis set

        mol_1 = Molecule.read_xyz(h2ofile)
        bas_1 = MolecularBasis.read(mol_1, 'def2-svp')

        S11 = ovldrv.compute(mol_1, bas_1)
        T11 = kindrv.compute(mol_1, bas_1)
        V11 = npotdrv.compute(mol_1, bas_1)

        V11p = npotdrv.compute(mol_1, bas_1, mol_1)
        self.assertTrue((V11.to_numpy() == V11p.to_numpy()).all())

        if rank == mpi_master():

            hf = h5py.File(h5file, 'r')
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

        mol_1 = Molecule.read_xyz(h2ofile)
        bas_1 = MolecularBasis.read(mol_1, 'def2-svp')
        bas_2 = MolecularBasis.read(mol_1, 'cc-pvdz')

        S12 = ovldrv.compute(mol_1, bas_1, bas_2)
        T12 = kindrv.compute(mol_1, bas_1, bas_2)
        V12 = npotdrv.compute(mol_1, bas_1, bas_2, mol_1)

        if rank == mpi_master():

            hf = h5py.File(h5file, 'r')
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

        mol_1 = Molecule.read_xyz(h2ofile)
        mol_2 = Molecule.read_xyz(nh3file)
        mol = Molecule(mol_1, mol_2)
        bas = MolecularBasis.read(mol, 'def2-svp')

        S12 = ovldrv.compute(mol_1, mol_2, bas)
        T12 = kindrv.compute(mol_1, mol_2, bas)
        V12 = npotdrv.compute(mol_1, mol_2, bas, mol)

        if rank == mpi_master():

            hf = h5py.File(h5file, 'r')
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

        mol_1 = Molecule.read_xyz(h2ofile)
        mol_2 = Molecule.read_xyz(nh3file)
        mol = Molecule(mol_1, mol_2)
        bas_1 = MolecularBasis.read(mol_1, 'def2-svp')
        bas_2 = MolecularBasis.read(mol_2, 'cc-pvdz')

        S12 = ovldrv.compute(mol_1, mol_2, bas_1, bas_2)
        T12 = kindrv.compute(mol_1, mol_2, bas_1, bas_2)
        V12 = npotdrv.compute(mol_1, mol_2, bas_1, bas_2, mol)

        if rank == mpi_master():

            hf = h5py.File(h5file, 'r')
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
