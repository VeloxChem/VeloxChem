from mpi4py import MPI
from veloxchem.mpitask import MpiTask
from veloxchem.veloxchemlib import Molecule
from veloxchem.veloxchemlib import MolecularBasis
from veloxchem.veloxchemlib import denmat
from veloxchem.veloxchemlib import molorb
from veloxchem.veloxchemlib import mpi_master

from veloxchem.aodensitymatrix import AODensityMatrix
from veloxchem.molecularorbitals import MolecularOrbitals
from veloxchem.inputparser import InputParser

import numpy as np
import unittest


class TestOrbData(unittest.TestCase):

    def test_get_label(self):

        task = MpiTask(["inputs/dimer.inp", "inputs/dimer.out"], MPI.COMM_WORLD)
        basis = task.ao_basis

        self.assertEqual(basis.get_label(), "DEF2-SVP")

    def test_density_matrix(self):

        data_a = [[ 1., .2, ], [ .2, 1., ]]
        data_b = [[ .9, .5, ], [ .5, .9, ]]

        d_rest = AODensityMatrix([data_a], denmat.rest)
        d_unrest = AODensityMatrix([data_a, data_b], denmat.unrest)

        den_a1 = d_rest.total_to_numpy(0)
        den_a2 = d_unrest.alpha_to_numpy(0)
        den_b2 = d_unrest.beta_to_numpy(0)

        self.assertTrue((data_a == den_a1).all())
        self.assertTrue((data_a == den_a2).all())
        self.assertTrue((data_b == den_b2).all())

        self.assertEqual(denmat.rest, d_rest.get_density_type())
        self.assertEqual(denmat.unrest, d_unrest.get_density_type())

        self.assertEqual(1, d_rest.number_of_density_matrices())
        self.assertEqual(1, d_unrest.number_of_density_matrices())

    def test_density_sub(self):

        arr_1 = np.array([[ 1., .2, ], [ .2,  1., ]])
        arr_2 = np.array([[ .9, .5, ], [ .5,  .9, ]])

        den_1 = AODensityMatrix([arr_1], denmat.rest)
        den_2 = AODensityMatrix([arr_2], denmat.rest)
        den_diff = den_1.sub(den_2)

        diff = np.max(np.abs(den_diff.total_to_numpy(0) - (arr_1 - arr_2)))
        self.assertAlmostEqual(0., diff, 13)

    def test_density_hdf5(self):

        data_a = [[ 1., .2, ], [ .2, 1., ]]
        data_b = [[ .9, .5, ], [ .5, .9, ]]
        data_c = [[ .8, .6, ], [ .6, .8, ]]
        data_d = [[ .7, .5, ], [ .5, .7, ]]

        d_rest = AODensityMatrix([data_a, data_b], denmat.rest)

        d_unrest = AODensityMatrix(
            [data_a, data_b, data_c, data_d], denmat.unrest)

        # hdf5 read/write tests

        if MPI.COMM_WORLD.Get_rank() == mpi_master():

            d_rest.write_hdf5("inputs/dummy.h5")
            dummy = AODensityMatrix.read_hdf5("inputs/dummy.h5")
            self.assertEqual(d_rest, dummy)

            d_unrest.write_hdf5("inputs/dummy.h5")
            dummy = AODensityMatrix.read_hdf5("inputs/dummy.h5")
            self.assertEqual(d_unrest, dummy)

    def test_orbitals_matrix(self):

        data_a = [[ .9, .2, ], [ .1, .3, ], [ .4, .9, ]]
        data_b = [[ .8, .3, ], [ .2, .4, ], [ .5, .8, ]]

        ener_a = [0.5, 0.9]
        ener_b = [0.3, 0.6]

        orb_rest = MolecularOrbitals([data_a], [ener_a], molorb.rest)

        orb_unrest = MolecularOrbitals(
            [data_a, data_b], [ener_a, ener_b], molorb.unrest)

        orb_a1 = orb_rest.alpha_to_numpy()
        orb_a2 = orb_unrest.alpha_to_numpy()
        orb_b2 = orb_unrest.beta_to_numpy()

        self.assertTrue((data_a == orb_a1).all())
        self.assertTrue((data_a == orb_a2).all())
        self.assertTrue((data_b == orb_b2).all())

        ene_a1 = orb_rest.ea_to_numpy()
        ene_a2 = orb_unrest.ea_to_numpy()
        ene_b2 = orb_unrest.eb_to_numpy()

        self.assertTrue((ener_a == ene_a1).all())
        self.assertTrue((ener_a == ene_a2).all())
        self.assertTrue((ener_b == ene_b2).all())

        self.assertEqual(molorb.rest, orb_rest.get_orbitals_type())
        self.assertEqual(molorb.unrest, orb_unrest.get_orbitals_type())
        
        self.assertEqual(2, orb_rest.get_number_mos())
        self.assertEqual(2, orb_unrest.get_number_mos())
        
        self.assertEqual(3, orb_rest.get_number_aos())
        self.assertEqual(3, orb_unrest.get_number_aos())

        # hdf5 read/write tests

        if MPI.COMM_WORLD.Get_rank() == mpi_master():

            orb_rest.write_hdf5("inputs/dummy.h5")
            dummy = MolecularOrbitals.read_hdf5("inputs/dummy.h5")
            self.assertEqual(orb_rest, dummy)

            orb_unrest.write_hdf5("inputs/dummy.h5")
            dummy = MolecularOrbitals.read_hdf5("inputs/dummy.h5")
            self.assertEqual(orb_unrest, dummy)

    def test_rest_density(self):

        arr = np.array([[.9, .2, .3], [.3, .8, .6], [.1, .5, .7]])
        ene = [.7, .8, .9]

        orb_rest = MolecularOrbitals([arr], [ene], molorb.rest)

        mol = Molecule(["H", "H"], [0.0, 0.0], [0.0, 0.0], [0.0, 1.4])
        den_rest = orb_rest.get_density(mol).total_to_numpy(0)

        arr_occ = arr[:, :1]
        den_ref = np.dot(arr_occ, arr_occ.T)

        self.assertTrue((den_ref == den_rest).all())


if __name__ == "__main__":
    unittest.main()
