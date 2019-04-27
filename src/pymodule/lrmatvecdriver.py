import numpy as np
import itertools

from .veloxchemlib import ElectronRepulsionIntegralsDriver
from .veloxchemlib import ElectricDipoleIntegralsDriver
from .veloxchemlib import NuclearPotentialIntegralsDriver
from .veloxchemlib import KineticEnergyIntegralsDriver
from .veloxchemlib import OverlapIntegralsDriver
from .veloxchemlib import AODensityMatrix
from .veloxchemlib import AOFockMatrix
from .veloxchemlib import ExcitationVector
from .veloxchemlib import mpi_master
from .veloxchemlib import denmat
from .veloxchemlib import fockmat
from .veloxchemlib import szblock
from .errorhandler import assert_msg_critical


class LinearResponseMatrixVectorDriver:
    """Implements linear response solver"""

    def __init__(self, comm):
        """Initializes linear response matrix vector driver.

        Initializes linear response matrix vector driver to default setup.

        Parameters
        ----------
        comm
            The MPI communicator.
        """

        # mpi information
        self.comm = comm
        self.rank = self.comm.Get_rank()
        self.nodes = self.comm.Get_size()

    def comp_1e_ints(self, molecule, basis):
        """Computes 1e integrals"""

        overlap_drv = OverlapIntegralsDriver(self.comm)
        kinetic_drv = KineticEnergyIntegralsDriver(self.comm)
        potential_drv = NuclearPotentialIntegralsDriver(self.comm)
        dipole_drv = ElectricDipoleIntegralsDriver(self.comm)

        S = overlap_drv.compute(molecule, basis)
        T = kinetic_drv.compute(molecule, basis)
        V = potential_drv.compute(molecule, basis)
        Dpl = dipole_drv.compute(molecule, basis)

        if self.rank == mpi_master():
            overlap = S.to_numpy()
            hcore = T.to_numpy() - V.to_numpy()
            dipoles = (Dpl.x_to_numpy(), Dpl.y_to_numpy(), Dpl.z_to_numpy())
            return overlap, hcore, dipoles
        else:
            return None, None, None

    def comp_fock(self, hcore, dens, screening, molecule, basis):
        """Computes Fock matrices"""

        if self.rank == mpi_master():
            dabs = (dens,)
        else:
            dabs = None

        fabs = self.get_two_el_fock(dabs, screening, molecule, basis)

        if self.rank == mpi_master():
            fa, fb = fabs[0]
            fa += hcore
            fb += hcore
            return fa, fb
        else:
            return None, None

    def e2n(self, vecs, tensors, screening, molecule, basis):

        if self.rank == mpi_master():
            assert_msg_critical(
                len(vecs.shape) == 2,
                'LinearResponseSolver.e2n: invalid shape of vecs')

            mo = tensors['C']
            S = tensors['S']
            da, db = tensors['D']
            fa, fb = tensors['F']

            nocc = molecule.number_of_alpha_electrons()
            norb = mo.shape[1]

            dks = []
            kns = []

            for col in range(vecs.shape[1]):
                vec = vecs[:, col]

                kN = self.vec2mat(vec, nocc, norb).T
                kn = mo @ kN @ mo.T

                dak = kn.T @ S @ da - da @ S @ kn.T
                dbk = kn.T @ S @ db - db @ S @ kn.T

                dks.append((dak, dbk))
                kns.append(kn)

            dks = tuple(dks)
        else:
            dks = None

        fks = self.get_two_el_fock(dks, screening, molecule, basis)

        if self.rank == mpi_master():
            gv = np.zeros(vecs.shape)

            for col, (kn, (fak, fbk)) in enumerate(zip(kns, fks)):

                kfa = S @ kn @ fa - fa @ kn @ S
                kfb = S @ kn @ fb - fb @ kn @ S

                fat = fak + kfa
                fbt = fbk + kfb

                gao = S @ (da @ fat.T + db @ fbt.T) - (fat.T @ da +
                                                       fbt.T @ db) @ S
                gmo = mo.T @ gao @ mo

                gv[:, col] = -self.mat2vec(gmo, nocc, norb)
            return gv
        else:
            return None

    def get_two_el_fock(self, dabs, screening, molecule, basis):

        # TODO: make this routine more general (for both rest and unrest)

        if self.rank == mpi_master():
            dts = []
            for dab in dabs:
                da, db = dab
                dt = da + db
                dts.append(dt)

                # Note: skip spin density for restricted case
                # ds = da - db
                # dts.append(ds)

            dens = AODensityMatrix(dts, denmat.rest)
        else:
            dens = AODensityMatrix()
        dens.broadcast(self.rank, self.comm)

        fock = AOFockMatrix(dens)

        # for i in range(0, 2 * len(dabs), 2):
        #    fock.set_fock_type(fockmat.rgenjk, i)
        #    fock.set_fock_type(fockmat.rgenk, i + 1)

        # Note: skip spin density for restricted case
        # for i in range(len(dabs)):
        #    fock.set_fock_type(fockmat.rgenjk, i)
        for i in range(fock.number_of_fock_matrices()):
            fock.set_fock_type(fockmat.rgenjk, i)

        eri_drv = ElectronRepulsionIntegralsDriver(self.comm)
        eri_drv.compute(fock, dens, molecule, basis, screening)
        fock.reduce_sum(self.rank, self.nodes, self.comm)

        fabs = []
        if self.rank == mpi_master():

            # for i in range(0, 2 * len(dabs), 2):
            #    ft = fock.to_numpy(i).T
            #    fs = -fock.to_numpy(i + 1).T
            #    fa = (ft + fs) / 2
            #    fb = (ft - fs) / 2
            #    fabs.append((fa, fb))

            # Note: skip spin density for restricted case
            for i in range(len(dabs)):
                ft = fock.to_numpy(i).T
                fa = ft * 0.5
                fb = ft * 0.5
                fabs.append((fa, fb))

            return tuple(fabs)
        else:
            return None

    def s2n(self, vecs, tensors, shapes):

        assert_msg_critical(
            len(vecs.shape) == 2,
            'LinearResponseSolver.s2n: invalid shape of vecs')

        mo = tensors['C']
        S = tensors['S']
        D = tensors['D'][0] + tensors['D'][1]

        nocc = shapes['nocc']
        norb = shapes['norb']

        s2n_vecs = np.ndarray(vecs.shape)
        rows, columns = vecs.shape

        for c in range(columns):
            kappa = self.vec2mat(vecs[:, c], nocc, norb).T
            kappa_ao = mo @ kappa @ mo.T

            s2n_ao = kappa_ao.T @ S @ D - D @ S @ kappa_ao.T
            s2n_mo = mo.T @ S @ s2n_ao @ S @ mo
            s2n_vecs[:, c] = -self.mat2vec(s2n_mo, nocc, norb)

        return s2n_vecs

    def vec2mat(self, vec, nocc, norb):

        zlen = len(vec) // 2
        z, y = vec[:zlen], vec[zlen:]

        xv = ExcitationVector(szblock.aa, 0, nocc, nocc, norb, True)
        xv.set_yzcoefficients(z, y)

        kz = xv.get_zmatrix()
        ky = xv.get_ymatrix()

        rows = kz.number_of_rows() + ky.number_of_rows()
        cols = kz.number_of_columns() + ky.number_of_columns()

        kzy = np.zeros((rows, cols))
        kzy[:kz.number_of_rows(), ky.number_of_columns():] = kz.to_numpy()
        kzy[kz.number_of_rows():, :ky.number_of_columns()] = ky.to_numpy()

        return kzy

    def mat2vec(self, mat, nocc, norb):

        xv = ExcitationVector(szblock.aa, 0, nocc, nocc, norb, True)
        excitations = list(
            itertools.product(xv.bra_unique_indexes(), xv.ket_unique_indexes()))

        z = [mat[i, j] for i, j in excitations]
        y = [mat[j, i] for i, j in excitations]
        return np.array(z + y)
