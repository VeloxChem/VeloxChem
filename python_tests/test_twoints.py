from pathlib import Path
import numpy as np
import math

from veloxchem.veloxchemlib import DenseMatrix
from veloxchem.veloxchemlib import KineticEnergyMatrix
from veloxchem.veloxchemlib import NuclearPotentialMatrix
from veloxchem.veloxchemlib import ElectronRepulsionIntegralsDriver
from veloxchem.veloxchemlib import denmat
from veloxchem.veloxchemlib import fockmat
from veloxchem.veloxchemlib import ericut
from veloxchem.veloxchemlib import is_mpi_master
from veloxchem.mpitask import MpiTask
from veloxchem.aodensitymatrix import AODensityMatrix
from veloxchem.aofockmatrix import AOFockMatrix
from veloxchem.qqscheme import get_qq_scheme


class TestTwoInts:

    def test_restricted_fock_matrix(self):

        data_j = [[1., .2], [.2, 1.]]
        data_k = [[.9, .5], [.5, .9]]

        arr_j = np.array(data_j)
        arr_k = np.array(data_k)
        arr_jk = 2.0 * arr_j - arr_k

        x = 0.5
        arr_kx = x * arr_k
        arr_jkx = 2.0 * arr_j - arr_kx

        fock = AOFockMatrix(
            [arr_jk, arr_jkx, arr_j, arr_k, arr_kx],
            [
                fockmat.restjk, fockmat.restjkx, fockmat.restj, fockmat.restk,
                fockmat.restkx
            ],
            [1.0, x, 1.0, 1.0, x],
            [0, 0, 0, 0, 0],
        )

        np_jkx = fock.to_numpy(1)
        np_j = fock.to_numpy(2)
        np_k = fock.to_numpy(3)

        assert (arr_j == np_j).all()
        assert (arr_k == np_k).all()
        assert (arr_jkx == np_jkx).all()

        assert fock.get_fock_type(0) == fockmat.restjk
        assert fock.get_scale_factor(1) == x
        assert fock.get_density_identifier(2) == 0
        assert fock.number_of_fock_matrices() == 5
        assert fock.is_closed_shell()

        fock.set_scale_factor(2.0 * x, 1)
        assert fock.get_scale_factor(1) == 2.0 * x

    def test_unrestricted_fock_matrix(self):

        data_a = [[1., .2], [.2, 1.]]
        data_b = [[.9, .5], [.5, .9]]

        arr_a = np.array(data_a)
        arr_b = np.array(data_b)

        fock = AOFockMatrix(
            [arr_a, arr_b],
            [fockmat.unrestjk, fockmat.unrestjk],
            [0.7, 0.6],
            [1, 2],
        )

        assert fock.number_of_fock_matrices() == 1

        assert (fock.alpha_to_numpy(0) == arr_a).all()
        assert (fock.beta_to_numpy(0) == arr_b).all()

        assert fock.get_fock_type(0, 'alpha') == fockmat.unrestjk
        assert fock.get_fock_type(0, 'beta') == fockmat.unrestjk

        assert fock.get_scale_factor(0, 'alpha') == 0.7
        assert fock.get_scale_factor(0, 'beta') == 0.6

        assert fock.get_density_identifier(0, 'alpha') == 1
        assert fock.get_density_identifier(0, 'beta') == 2

        assert fock.number_of_fock_matrices() == 1
        assert not fock.is_closed_shell()

    def test_add_hcore(self):

        arr_t = np.array([[3., .2], [.2, 3.]])
        arr_v = np.array([[-9., .5], [.5, -9.]])
        arr_jk = np.array([[5., .1], [.1, 5.]])
        arr_fock = arr_jk + arr_t - arr_v

        kin = KineticEnergyMatrix(arr_t)
        npot = NuclearPotentialMatrix(arr_v)
        fock = AOFockMatrix([arr_jk], [fockmat.restjk], [1.0], [0])
        fock.add_hcore(kin, npot, 0)

        diff = np.max(np.abs(fock.to_numpy(0) - arr_fock))
        assert abs(diff) < 1.0e-13

        fock.add_matrix(DenseMatrix(arr_t), 0)
        diff = np.max(np.abs(fock.to_numpy(0) - (arr_fock + arr_t)))
        assert abs(diff) < 1.0e-13

        fock.scale(0.5, 0)
        diff = np.max(np.abs(fock.to_numpy(0) - 0.5 * (arr_fock + arr_t)))
        assert abs(diff) < 1.0e-13

    def test_add_fock(self):

        arr_1 = np.array([[1., .2], [.2, 1.]])
        arr_2 = np.array([[.9, .5], [.5, .9]])

        fock_1 = AOFockMatrix([arr_1], [fockmat.restjk], [1.0], [0])
        fock_2 = AOFockMatrix([arr_2], [fockmat.restjk], [1.0], [0])

        fock_sum = AOFockMatrix(fock_1)
        fock_sum.add(fock_2)

        diff = np.max(np.abs(fock_sum.to_numpy(0) - (arr_1 + arr_2)))
        assert abs(diff) < 1.0e-13

    def test_get_energy(self):

        arr_a = np.array([[1., .2], [.2, 1.]])
        arr_b = np.array([[.9, .5], [.5, .9]])

        fock_rest = AOFockMatrix([arr_a], [fockmat.restjk], [1.0], [0])
        fock_unrest = AOFockMatrix([arr_a, arr_b],
                                   [fockmat.unrestjk, fockmat.unrestjk],
                                   [1.0, 1.0], [0, 0])

        dens_rest = AODensityMatrix([arr_a], denmat.rest)
        dens_unrest = AODensityMatrix([arr_a, arr_b], denmat.unrest)

        energy_rest = fock_rest.get_energy(0, dens_rest, 0)
        energy_unrest = fock_unrest.get_energy(0, dens_unrest, 0)

        energy_a = np.trace(np.matmul(arr_a, arr_a))
        energy_b = np.trace(np.matmul(arr_b, arr_b))

        assert abs(energy_rest - energy_a) < 1.0e-12
        assert abs(energy_unrest - 0.5 * (energy_a + energy_b)) < 1.0e-12

    def test_fock_density(self):

        data_a = [[1., .2], [.2, 1.]]

        d_rest = AODensityMatrix([data_a], denmat.rest)

        f_rest = AOFockMatrix(d_rest)

        assert f_rest.get_fock_type(0) == fockmat.restjk
        assert f_rest.get_scale_factor(0) == 1.0
        assert f_rest.get_density_identifier(0) == 0

    def test_fock_build(self):

        here = Path(__file__).parent
        inpfile = here / 'inputs' / 'h2se.inp'
        outfile = inpfile.with_suffix('.out')

        task = MpiTask([str(inpfile), str(outfile)])

        molecule = task.molecule
        ao_basis = task.ao_basis

        molecule.check_proximity(0.1)

        molecule.check_multiplicity()

        enuc = molecule.nuclear_repulsion_energy()

        ref_enuc = 34.0 / 2.8 + 34.0 / 2.8 + 1.0 / (2.8 * math.sqrt(2.0))

        assert abs(enuc - ref_enuc) < 1.0e-13

        comm = task.mpi_comm
        rank = task.mpi_rank
        size = task.mpi_size

        # read density

        if is_mpi_master(comm):
            densfile = str(here / 'inputs' / 'h2se.dens.h5')

            dmat = AODensityMatrix.read_hdf5(densfile)
        else:
            dmat = AODensityMatrix()

        dmat.broadcast(rank, comm)

        # compute Fock

        eridrv = ElectronRepulsionIntegralsDriver(comm)

        qqdata = eridrv.compute(ericut.qqden, 1.0e-12, molecule, ao_basis)

        assert get_qq_scheme('QQ') == ericut.qq
        assert get_qq_scheme('QQR') == ericut.qqr
        assert get_qq_scheme('QQ_DEN') == ericut.qqden
        assert get_qq_scheme('QQR_DEN') == ericut.qqrden

        fock = AOFockMatrix(dmat)

        eridrv.compute(fock, dmat, molecule, ao_basis, qqdata)

        fock.reduce_sum(rank, size, comm)

        # compare with unrestricted Fock

        if is_mpi_master(comm):
            da = dmat.alpha_to_numpy(0)
            db = dmat.beta_to_numpy(0)
            d2 = AODensityMatrix([da, db], denmat.unrest)
        else:
            d2 = AODensityMatrix()
        d2.broadcast(rank, comm)

        f2 = AOFockMatrix(d2)
        eridrv.compute(f2, d2, molecule, ao_basis, qqdata)
        f2.reduce_sum(rank, size, comm)

        assert f2.number_of_fock_matrices() == 1
        assert f2.get_fock_type(0, 'alpha') == fockmat.unrestjk
        assert f2.get_fock_type(0, 'beta') == fockmat.unrestjk

        assert np.max(
            np.abs(fock.alpha_to_numpy(0) - f2.alpha_to_numpy(0))) < 1.0e-13
        assert np.max(
            np.abs(fock.beta_to_numpy(0) - f2.beta_to_numpy(0))) < 1.0e-13

        # compare with reference

        if is_mpi_master(comm):

            twoefile = str(here / 'inputs' / 'h2se.twoe.h5')

            fock_ref = AOFockMatrix.read_hdf5(twoefile)

            F1 = fock.to_numpy(0)
            F2 = fock_ref.to_numpy(0)
            dF = np.max(np.abs(F1 - F2))

            assert dF < 1.0e-11
