from pathlib import Path
import numpy as np
import h5py

from veloxchem.veloxchemlib import DenseMatrix
from veloxchem.veloxchemlib import OverlapMatrix
from veloxchem.veloxchemlib import KineticEnergyMatrix
from veloxchem.veloxchemlib import NuclearPotentialMatrix
from veloxchem.veloxchemlib import OverlapIntegralsDriver
from veloxchem.veloxchemlib import KineticEnergyIntegralsDriver
from veloxchem.veloxchemlib import NuclearPotentialIntegralsDriver
from veloxchem.veloxchemlib import ElectricDipoleMatrix
from veloxchem.veloxchemlib import LinearMomentumMatrix
from veloxchem.veloxchemlib import AngularMomentumMatrix
from veloxchem.veloxchemlib import ElectricFieldMatrix
from veloxchem.veloxchemlib import is_mpi_master
from veloxchem.mpitask import MpiTask
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis


class TestOneInts:

    def test_overlap_matrix(self):

        array = np.array([[1.0, 0.2], [0.2, 1.0]])
        matrix = OverlapMatrix(array)
        array2 = matrix.to_numpy()
        matrix2 = OverlapMatrix(array2)

        assert (array == array2).all()
        assert matrix == matrix2

    def test_kinetic_energy_matrix(self):

        array = np.array([[1.0, 0.2], [0.2, 1.0]])
        matrix = KineticEnergyMatrix(array)
        array2 = matrix.to_numpy()
        matrix2 = KineticEnergyMatrix(array2)

        assert (array == array2).all()
        assert matrix == matrix2

    def test_nuclear_potential_matrix(self):

        array = np.array([[1.0, 0.2], [0.2, 1.0]])
        matrix = NuclearPotentialMatrix(array)
        array2 = matrix.to_numpy()
        matrix2 = NuclearPotentialMatrix(array2)

        assert (array == array2).all()
        assert matrix == matrix2

    def test_electric_dipole_matrix(self):

        x_array = DenseMatrix(np.random.rand(2, 2))
        y_array = DenseMatrix(np.random.rand(2, 2))
        z_array = DenseMatrix(np.random.rand(2, 2))
        origin = list(np.random.rand(3))

        matrix = ElectricDipoleMatrix(x_array, y_array, z_array, *origin)
        matrix2 = ElectricDipoleMatrix([x_array, y_array, z_array], origin)

        assert matrix == matrix2
        assert matrix.origin == origin
        assert (matrix.x_to_numpy() == x_array.to_numpy()).all()
        assert (matrix.y_to_numpy() == y_array.to_numpy()).all()
        assert (matrix.z_to_numpy() == z_array.to_numpy()).all()

    def test_linear_momentum_matrix(self):

        x_array = DenseMatrix(np.random.rand(2, 2))
        y_array = DenseMatrix(np.random.rand(2, 2))
        z_array = DenseMatrix(np.random.rand(2, 2))

        matrix = LinearMomentumMatrix(x_array, y_array, z_array)
        matrix2 = LinearMomentumMatrix([x_array, y_array, z_array])

        assert matrix == matrix2
        assert (matrix.x_to_numpy() == x_array.to_numpy()).all()
        assert (matrix.y_to_numpy() == y_array.to_numpy()).all()
        assert (matrix.z_to_numpy() == z_array.to_numpy()).all()

    def test_angular_momentum_matrix(self):

        x_array = DenseMatrix(np.random.rand(2, 2))
        y_array = DenseMatrix(np.random.rand(2, 2))
        z_array = DenseMatrix(np.random.rand(2, 2))
        origin = list(np.random.rand(3))

        matrix = AngularMomentumMatrix(x_array, y_array, z_array, *origin)
        matrix2 = AngularMomentumMatrix([x_array, y_array, z_array], origin)

        assert matrix == matrix2
        assert matrix.origin == origin
        assert (matrix.x_to_numpy() == x_array.to_numpy()).all()
        assert (matrix.y_to_numpy() == y_array.to_numpy()).all()
        assert (matrix.z_to_numpy() == z_array.to_numpy()).all()

    def test_electric_field_matrix(self):

        x_array = DenseMatrix(np.random.rand(2, 2))
        y_array = DenseMatrix(np.random.rand(2, 2))
        z_array = DenseMatrix(np.random.rand(2, 2))

        matrix = ElectricFieldMatrix(x_array, y_array, z_array)
        matrix2 = ElectricFieldMatrix([x_array, y_array, z_array])

        assert matrix == matrix2
        assert (matrix.x_to_numpy() == x_array.to_numpy()).all()
        assert (matrix.y_to_numpy() == y_array.to_numpy()).all()
        assert (matrix.z_to_numpy() == z_array.to_numpy()).all()

    def test_1e_integrals(self):

        here = Path(__file__).parent
        inpfile = here / 'inputs' / 'h2se.inp'
        outfile = inpfile.with_suffix('.out')

        task = MpiTask([str(inpfile), str(outfile)])

        molecule = task.molecule
        basis = task.ao_basis

        # compute 1e integrals

        ovldrv = OverlapIntegralsDriver(task.mpi_comm)
        S = ovldrv.compute(molecule, basis)
        S1 = S.to_numpy()

        kindrv = KineticEnergyIntegralsDriver(task.mpi_comm)
        T = kindrv.compute(molecule, basis)
        T1 = T.to_numpy()

        npotdrv = NuclearPotentialIntegralsDriver(task.mpi_comm)
        V = npotdrv.compute(molecule, basis)
        V1 = V.to_numpy()

        # compare with reference

        if is_mpi_master(task.mpi_comm):

            h5file = here / 'inputs' / 'h2se.onee.h5'

            hf = h5py.File(h5file, 'r')
            S2 = np.array(hf.get('overlap'))
            T2 = np.array(hf.get('kinetic_energy'))
            V2 = np.array(hf.get('nuclear_potential'))
            hf.close()

            dS = np.max(np.abs(S1 - S2))
            dT = np.max(np.abs(T1 - T2))
            dV = np.max(np.abs(V1 - V2))

            assert dS < 1.0e-13
            assert dT < 1.0e-11
            assert dV < 1.0e-11

    def test_mixed_basis_1e(self):

        ovldrv = OverlapIntegralsDriver()
        kindrv = KineticEnergyIntegralsDriver()
        npotdrv = NuclearPotentialIntegralsDriver()

        here = Path(__file__).parent
        h2ofile = here / 'inputs' / 'h2o.xyz'

        nh3file = here / 'inputs' / 'nh3.xyz'

        h5file = here / 'inputs' / 'mix_basis_1e.h5'

        # one molecule, one basis set

        mol_1 = Molecule.read_xyz(h2ofile)
        bas_1 = MolecularBasis.read(mol_1, 'def2-svp')

        S11 = ovldrv.compute(mol_1, bas_1)
        T11 = kindrv.compute(mol_1, bas_1)
        V11 = npotdrv.compute(mol_1, bas_1)

        V11p = npotdrv.compute(mol_1, bas_1, mol_1)
        assert (V11.to_numpy() == V11p.to_numpy()).all()

        if is_mpi_master():

            hf = h5py.File(h5file, 'r')
            ref_S11 = np.array(hf.get('S_h2o_def2-svp'))
            ref_T11 = np.array(hf.get('T_h2o_def2-svp'))
            ref_V11 = np.array(hf.get('V_h2o_def2-svp'))
            hf.close()

            dS = np.max(np.abs(S11.to_numpy() - ref_S11))
            dT = np.max(np.abs(T11.to_numpy() - ref_T11))
            dV = np.max(np.abs(V11.to_numpy() - ref_V11))

            assert dS < 1.0e-13
            assert dT < 1.0e-13
            assert dV < 1.0e-12

        # one molecule, two basis sets

        mol_1 = Molecule.read_xyz(h2ofile)
        bas_1 = MolecularBasis.read(mol_1, 'def2-svp')
        bas_2 = MolecularBasis.read(mol_1, 'cc-pvdz')

        S12 = ovldrv.compute(mol_1, bas_1, bas_2)
        T12 = kindrv.compute(mol_1, bas_1, bas_2)
        V12 = npotdrv.compute(mol_1, bas_1, bas_2, mol_1)

        if is_mpi_master():

            hf = h5py.File(h5file, 'r')
            ref_S12 = np.array(hf.get('S_h2o_def2-svp_cc-pvdz'))
            ref_T12 = np.array(hf.get('T_h2o_def2-svp_cc-pvdz'))
            ref_V12 = np.array(hf.get('V_h2o_def2-svp_cc-pvdz'))
            hf.close()

            dS = np.max(np.abs(S12.to_numpy() - ref_S12))
            dT = np.max(np.abs(T12.to_numpy() - ref_T12))
            dV = np.max(np.abs(V12.to_numpy() - ref_V12))

            assert dS < 1.0e-13
            assert dT < 1.0e-13
            assert dV < 1.0e-12

        # two molecules, one basis set

        mol_1 = Molecule.read_xyz(h2ofile)
        mol_2 = Molecule.read_xyz(nh3file)
        mol = Molecule(mol_1, mol_2)
        bas = MolecularBasis.read(mol, 'def2-svp')

        S12 = ovldrv.compute(mol_1, mol_2, bas)
        T12 = kindrv.compute(mol_1, mol_2, bas)
        V12 = npotdrv.compute(mol_1, mol_2, bas, mol)

        if is_mpi_master():

            hf = h5py.File(h5file, 'r')
            ref_S12 = np.array(hf.get('S_h2o_nh3_def2-svp'))
            ref_T12 = np.array(hf.get('T_h2o_nh3_def2-svp'))
            ref_V12 = np.array(hf.get('V_h2o_nh3_def2-svp'))
            hf.close()

            dS = np.max(np.abs(S12.to_numpy() - ref_S12))
            dT = np.max(np.abs(T12.to_numpy() - ref_T12))
            dV = np.max(np.abs(V12.to_numpy() - ref_V12))

            assert dS < 1.0e-13
            assert dT < 1.0e-11
            assert dV < 1.0e-11

        # two molecules, two basis sets

        mol_1 = Molecule.read_xyz(h2ofile)
        mol_2 = Molecule.read_xyz(nh3file)
        mol = Molecule(mol_1, mol_2)
        bas_1 = MolecularBasis.read(mol_1, 'def2-svp')
        bas_2 = MolecularBasis.read(mol_2, 'cc-pvdz')

        S12 = ovldrv.compute(mol_1, mol_2, bas_1, bas_2)
        T12 = kindrv.compute(mol_1, mol_2, bas_1, bas_2)
        V12 = npotdrv.compute(mol_1, mol_2, bas_1, bas_2, mol)

        if is_mpi_master():

            hf = h5py.File(h5file, 'r')
            ref_S12 = np.array(hf.get('S_h2o_def2-svp_nh3_cc-pvdz'))
            ref_T12 = np.array(hf.get('T_h2o_def2-svp_nh3_cc-pvdz'))
            ref_V12 = np.array(hf.get('V_h2o_def2-svp_nh3_cc-pvdz'))
            hf.close()

            dS = np.max(np.abs(S12.to_numpy() - ref_S12))
            dT = np.max(np.abs(T12.to_numpy() - ref_T12))
            dV = np.max(np.abs(V12.to_numpy() - ref_V12))

            assert dS < 1.0e-13
            assert dT < 1.0e-11
            assert dV < 1.0e-11
