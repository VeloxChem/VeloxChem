from pathlib import Path
import numpy as np

from veloxchem.veloxchemlib import OverlapIntegralsDriver
from veloxchem.veloxchemlib import DenseMatrix
from veloxchem.veloxchemlib import denmat
from veloxchem.veloxchemlib import molorb
from veloxchem.veloxchemlib import is_mpi_master
from veloxchem.veloxchemlib import ao_matrix_to_dalton
from veloxchem.veloxchemlib import get_basis_function_indices_for_atom
from veloxchem.mpitask import MpiTask
from veloxchem.molecule import Molecule
from veloxchem.aodensitymatrix import AODensityMatrix
from veloxchem.molecularorbitals import MolecularOrbitals


class TestOrbData:

    def test_get_label(self):

        here = Path(__file__).parent
        inpfile = here / 'inputs' / 'dimer.inp'
        outfile = inpfile.with_suffix('.out')

        task = MpiTask([str(inpfile), str(outfile)])
        assert task.ao_basis.get_label() == "DEF2-SVP"

    def test_density_matrix(self):

        data_a = [[1., .2], [.2, 1.]]
        data_b = [[.9, .5], [.5, .9]]

        d_rest = AODensityMatrix([data_a], denmat.rest)
        d_unrest = AODensityMatrix([data_a, data_b], denmat.unrest)

        den_a1 = d_rest.alpha_to_numpy(0)
        den_b1 = d_rest.beta_to_numpy(0)
        den_a2 = d_unrest.alpha_to_numpy(0)
        den_b2 = d_unrest.beta_to_numpy(0)

        assert (data_a == den_a1).all()
        assert (data_a == den_b1).all()
        assert (data_a == den_a2).all()
        assert (data_b == den_b2).all()

        assert d_rest.get_density_type() == denmat.rest
        assert d_unrest.get_density_type() == denmat.unrest

        assert d_rest.number_of_density_matrices() == 1
        assert d_unrest.number_of_density_matrices() == 1

        den_empty_1 = d_rest.alpha_to_numpy(1)
        den_empty_2 = d_rest.beta_to_numpy(3)
        den_empty_3 = d_unrest.alpha_to_numpy(2)
        den_empty_4 = d_unrest.beta_to_numpy(4)

        assert den_empty_1.size == 0
        assert den_empty_2.size == 0
        assert den_empty_3.size == 0
        assert den_empty_4.size == 0

    def test_density_sub(self):

        arr_1 = np.array([[1., .2], [.2, 1.]])
        arr_2 = np.array([[.9, .5], [.5, .9]])

        den_1 = AODensityMatrix([arr_1], denmat.rest)
        den_2 = AODensityMatrix([arr_2], denmat.rest)
        den_diff = den_1.sub(den_2)

        diff = np.max(np.abs(den_diff.alpha_to_numpy(0) - (arr_1 - arr_2)))
        assert abs(diff) < 1.0e-13

    def test_orbitals_matrix(self):

        data_a = [[.9, .2], [.1, .3], [.4, .9]]
        data_b = [[.8, .3], [.2, .4], [.5, .8]]

        ener_a = [0.5, 0.9]
        ener_b = [0.3, 0.6]

        orb_rest = MolecularOrbitals([data_a], [ener_a], molorb.rest)

        orb_unrest = MolecularOrbitals([data_a, data_b], [ener_a, ener_b],
                                       molorb.unrest)

        orb_a1 = orb_rest.alpha_to_numpy()
        orb_a2 = orb_unrest.alpha_to_numpy()
        orb_b2 = orb_unrest.beta_to_numpy()

        assert (data_a == orb_a1).all()
        assert (data_a == orb_a2).all()
        assert (data_b == orb_b2).all()

        ene_a1 = orb_rest.ea_to_numpy()
        ene_a2 = orb_unrest.ea_to_numpy()
        ene_b2 = orb_unrest.eb_to_numpy()

        assert (ener_a == ene_a1).all()
        assert (ener_a == ene_a2).all()
        assert (ener_b == ene_b2).all()

        assert orb_rest.get_orbitals_type() == molorb.rest
        assert orb_unrest.get_orbitals_type() == molorb.unrest

        assert orb_rest.number_mos() == 2
        assert orb_unrest.number_mos() == 2

        assert orb_rest.number_aos() == 3
        assert orb_unrest.number_aos() == 3

        # hdf5 read/write tests

        if is_mpi_master():

            here = Path(__file__).parent
            h5file = str(here / 'inputs' / 'dummy.h5')

            nuc_chg = np.array([1, 8, 1], dtype=np.int32)
            orb_rest.write_hdf5(h5file, nuc_chg, 'sto-3g')
            dummy = MolecularOrbitals.read_hdf5(h5file)
            assert orb_rest == dummy
            assert MolecularOrbitals.match_hdf5(h5file, nuc_chg, 'sto-3g',
                                                'restricted')

            nuc_chg = np.array([1, 1, 8], dtype=np.int32)
            orb_unrest.write_hdf5(h5file, nuc_chg, 'cc-pvdz')
            dummy = MolecularOrbitals.read_hdf5(h5file)
            assert orb_unrest == dummy
            assert MolecularOrbitals.match_hdf5(h5file, nuc_chg, 'cc-pvdz',
                                                'unrestricted')

    def test_rest_density(self):

        mol = Molecule(['H', 'H'], [[0.0, 0.0, 0.0], [0.0, 0.0, 1.4]])

        arr = np.array([[.9, .2, .3], [.3, .8, .6], [.1, .5, .7]])
        ene = np.array([.7, .8, .9])

        orb_rest = MolecularOrbitals([arr], [ene], molorb.rest)
        den_rest = orb_rest.get_density(mol, 'restricted')
        den_a = den_rest.alpha_to_numpy(0)
        den_b = den_rest.beta_to_numpy(0)

        arr_occ = arr[:, :1]
        den_ref = np.dot(arr_occ, arr_occ.T)

        assert (den_ref == den_a).all()
        assert (den_ref == den_b).all()

    def test_unrest_density(self):

        mol = Molecule(['H', 'H'], [[0.0, 0.0, 0.0], [0.0, 0.0, 1.4]])

        arr_a = np.array([[.9, .2, .3], [.3, .8, .6], [.1, .5, .7]])
        ene_a = np.array([.7, .8, .9])

        arr_b = arr_a * 0.9
        ene_b = ene_a * 0.8

        orb_unrest = MolecularOrbitals([arr_a, arr_b], [ene_a, ene_b],
                                       molorb.unrest)
        den_unrest = orb_unrest.get_density(mol, 'unrestricted')
        den_a = den_unrest.alpha_to_numpy(0)
        den_b = den_unrest.beta_to_numpy(0)

        arr_occ_a = arr_a[:, :1]
        den_ref_a = np.dot(arr_occ_a, arr_occ_a.T)

        arr_occ_b = arr_b[:, :1]
        den_ref_b = np.dot(arr_occ_b, arr_occ_b.T)

        assert (den_ref_a == den_a).all()
        assert (den_ref_b == den_b).all()

    def test_basis_function_indices(self):

        here = Path(__file__).parent
        inpfile = here / 'inputs' / 'h2se.inp'

        task = MpiTask([str(inpfile), None])

        bf_indices = []
        bf_angmoms = []

        for i in range(task.molecule.number_of_atoms()):
            indices, angmoms = get_basis_function_indices_for_atom(
                task.molecule, task.ao_basis, i)
            bf_indices += indices
            bf_angmoms += angmoms

        ref_angmoms = [
            0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2,
            2, 2, 2, 2, 2, 2, 2, 2, 2, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1
        ]
        assert bf_angmoms == ref_angmoms

        ovldrv = OverlapIntegralsDriver()
        S = ovldrv.compute(task.molecule, task.ao_basis)
        smat = S.to_numpy()

        if is_mpi_master():
            sdal = ao_matrix_to_dalton(DenseMatrix(smat), task.ao_basis,
                                       task.molecule).to_numpy()
            assert np.max(
                np.abs(sdal - smat[bf_indices, :][:, bf_indices])) < 1.0e-12
