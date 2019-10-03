from os.path import isfile
import numpy as np
import time as tm
import itertools
import h5py

from .veloxchemlib import ElectronRepulsionIntegralsDriver
from .veloxchemlib import ElectricDipoleIntegralsDriver
from .veloxchemlib import LinearMomentumIntegralsDriver
from .veloxchemlib import AngularMomentumIntegralsDriver
from .veloxchemlib import AODensityMatrix
from .veloxchemlib import AOFockMatrix
from .veloxchemlib import DenseMatrix
from .veloxchemlib import ExcitationVector
from .veloxchemlib import XCIntegrator
from .veloxchemlib import MolecularGrid
from .veloxchemlib import mpi_master
from .veloxchemlib import denmat
from .veloxchemlib import fockmat
from .veloxchemlib import szblock
from .errorhandler import assert_msg_critical
from .subcommunicators import SubCommunicators


class LinearResponseMatrixVectorDriver:
    """
    Implements linear response matrix-vector driver.

    :param comm:
        The MPI communicator.
    :param rank:
        The MPI rank.
    :param nodes:
        Number of MPI processes.
    """

    def __init__(self, comm):
        """
        Initializes linear response matrix vector driver to default setup.

        :param comm:
            The MPI communicator.
        """

        # mpi information
        self.comm = comm
        self.rank = self.comm.Get_rank()
        self.nodes = self.comm.Get_size()

        self.split_comm_ratio = None

    def e2n_half_size(self,
                      vecs_ger,
                      vecs_ung,
                      tensors,
                      screening,
                      molecule,
                      basis,
                      dft=False,
                      xcfun=None,
                      molgrid=None,
                      gs_density=None,
                      pe=False,
                      V_es=None,
                      potfile=None):
        """
        Computes the E2 b matrix vector product.

        :param vecs_ger:
            The gerade trial vectors in half-size.
        :param vecs_ung:
            The ungerade trial vectors in half-size.
        :param tensors:
            The dictionary of tensors from converged SCF wavefunction.
        :param screening:
            The electron repulsion integrals screening pattern.
        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        :param dft:
            The DFT flag if true compute XC contribution, if false otherwise.
        :param xcfun:
            The exchange correlation functional.
        :param molgrid:
            The molecular grid for XC contributtion computation.
        :param gs_density:
            The ground state density matrix.
        :param pe:
            The flag for running polarizable embedding calculation.
        :param V_es:
            The polarizable embedding matrix.
        :param potfile:
            The name of the potential file for polarizable embedding.

        :return:
            The gerade and ungerade E2 b matrix vector product in half-size.
        """

        if self.rank == mpi_master():
            assert_msg_critical(
                vecs_ger.ndim == 2 and vecs_ung.ndim == 2,
                'LinearResponse.e2n_half_size: invalid shape of trial vectors')

            mo = tensors['C']
            S = tensors['S']
            da, db = tensors['D']
            fa, fb = tensors['F']

            nocc = molecule.number_of_alpha_electrons()
            norb = mo.shape[1]

            dks = []
            kns = []

            n_ger = vecs_ger.shape[1] if vecs_ger is not None else 0
            n_ung = vecs_ung.shape[1] if vecs_ung is not None else 0

            for col in range(n_ger + n_ung):
                if col < n_ger:
                    # full-size gerade trial vector
                    vec = np.hstack((vecs_ger[:, col], vecs_ger[:, col]))
                else:
                    # full-size ungerade trial vector
                    vec = np.hstack(
                        (vecs_ung[:, col - n_ger], -vecs_ung[:, col - n_ger]))

                kN = lrvec2mat(vec, nocc, norb).T
                kn = np.linalg.multi_dot([mo, kN, mo.T])

                dak = np.linalg.multi_dot([kn.T, S, da])
                dak -= np.linalg.multi_dot([da, S, kn.T])
                # dbk = kn.T @ S @ db - db @ S @ kn.T

                dks.append(dak)
                kns.append(kn)
            dens = AODensityMatrix(dks, denmat.rest)
        else:
            dens = AODensityMatrix()
        dens.broadcast(self.rank, self.comm)

        fock = AOFockMatrix(dens)

        self.comp_lr_fock(fock, dens, molecule, basis, screening, dft, xcfun,
                          molgrid, gs_density, pe, V_es, potfile)

        if self.rank == mpi_master():
            if vecs_ger is not None:
                half_size = vecs_ger.shape[0]
            else:
                half_size = vecs_ung.shape[0]
            gv = np.zeros((half_size, n_ger + n_ung))

            for col, kn in enumerate(kns):

                fak = fock.alpha_to_numpy(col).T
                # fbk = fock.beta_to_numpy(col).T

                kfa = np.linalg.multi_dot([S, kn, fa])
                kfa -= np.linalg.multi_dot([fa, kn, S])
                # kfb = S @ kn @ fb - fb @ kn @ S

                fat = fak + kfa
                fbt = fat
                # fbt = fbk + kfb

                gao = np.matmul(S, np.matmul(da, fat.T) + np.matmul(db, fbt.T))
                gao -= np.matmul(np.matmul(fat.T, da) + np.matmul(fbt.T, db), S)

                gmo = np.linalg.multi_dot([mo.T, gao, mo])

                gv[:, col] = -lrmat2vec(gmo, nocc, norb)[:half_size]
            return gv[:, :n_ger], gv[:, n_ger:]
        else:
            return None, None

    def comp_lr_fock(self, fock, dens, molecule, basis, screening, dft, xcfun,
                     molgrid, gs_density, pe, V_es, potfile):
        """
        Computes Fock/Fxc matrix (2e part) for linear response calculation.

        :param fock:
            The Fock matrix (2e part).
        :param dens:
            The density matrix.
        :param molecule:
            The molecule.
        :param basis:
            The basis set.
        :param screening:
            The screening container object.
        :param dft:
            The DFT flag if true compute XC contribution, if false otherwise.
        :param xcfun:
            The exchange correlation functional.
        :param molgrid:
            The molecular grid for XC contributtion computation.
        :param gs_density:
            The ground state density matrix.
        :param pe:
            The flag for running polarizable embedding calculation.
        :param V_es:
            The polarizable embedding matrix.
        :param potfile:
            The name of the potential file for polarizable embedding.
        """

        # set flags for Fock matrices

        fock_flag = fockmat.rgenjk
        if dft:
            if xcfun.is_hybrid():
                fock_flag = fockmat.rgenjkx
                fact_xc = xcfun.get_frac_exact_exchange()
                for i in range(fock.number_of_fock_matrices()):
                    fock.set_scale_factor(fact_xc, i)
            else:
                fock_flag = fockmat.rgenj
        for i in range(fock.number_of_fock_matrices()):
            fock.set_fock_type(fock_flag, i)

        # calculate Fock on subcommunicators

        if (dft or pe) and self.nodes >= 4:
            self.comp_lr_fock_split_comm(fock, dens, molecule, basis, screening,
                                         dft, xcfun, molgrid, gs_density, pe,
                                         V_es, potfile)

        else:
            eri_drv = ElectronRepulsionIntegralsDriver(self.comm)
            eri_drv.compute(fock, dens, molecule, basis, screening)

            if dft:
                if not xcfun.is_hybrid():
                    for ifock in range(fock.number_of_fock_matrices()):
                        fock.scale(2.0, ifock)
                xc_drv = XCIntegrator(self.comm)
                molgrid.distribute(self.rank, self.nodes, self.comm)
                xc_drv.integrate(fock, dens, gs_density, molecule, basis,
                                 molgrid, xcfun.get_func_label())

            if pe:
                from .polembed import PolEmbed
                pe_drv = PolEmbed(molecule, basis, self.comm, potfile)
                pe_drv.V_es = V_es.copy()
                for ifock in range(fock.number_of_fock_matrices()):
                    dm = dens.alpha_to_numpy(ifock) + dens.beta_to_numpy(ifock)
                    e_pe, V_pe = pe_drv.get_pe_contribution(dm, elec_only=True)
                    if self.rank == mpi_master():
                        fock.add_matrix(DenseMatrix(V_pe), ifock)

            fock.reduce_sum(self.rank, self.nodes, self.comm)

    def comp_lr_fock_split_comm(self, fock, dens, molecule, basis, screening,
                                dft, xcfun, molgrid, gs_density, pe, V_es,
                                potfile):
        """
        Computes linear response Fock/Fxc matrix on split communicators.

        :param fock:
            The Fock matrix (2e part).
        :param dens:
            The density matrix.
        :param molecule:
            The molecule.
        :param basis:
            The basis set.
        :param screening:
            The screening container object.
        :param dft:
            The DFT flag if true compute XC contribution, if false otherwise.
        :param xcfun:
            The exchange correlation functional.
        :param molgrid:
            The molecular grid for XC contributtion computation.
        :param gs_density:
            The ground state density matrix.
        :param pe:
            The flag for running polarizable embedding calculation.
        :param V_es:
            The polarizable embedding matrix.
        :param potfile:
            The name of the potential file for polarizable embedding.
        """

        if self.split_comm_ratio is None:
            self.split_comm_ratio = [0.5, 0.25, 0.25]
            if not dft:
                self.split_comm_ratio = [0.5, 0.0, 0.5]
            if not pe:
                self.split_comm_ratio = [0.5, 0.5, 0.0]

        if dft:
            dft_nodes = int(float(self.nodes) * self.split_comm_ratio[1] + 0.5)
            dft_nodes = max(1, dft_nodes)
        else:
            dft_nodes = 0

        if pe:
            pe_nodes = int(float(self.nodes) * self.split_comm_ratio[2] + 0.5)
            pe_nodes = max(1, pe_nodes)
        else:
            pe_nodes = 0

        eri_nodes = max(1, self.nodes - dft_nodes - pe_nodes)

        if eri_nodes >= dft_nodes and eri_nodes >= pe_nodes:
            eri_nodes = self.nodes - dft_nodes - pe_nodes
        elif dft_nodes >= eri_nodes and dft_nodes >= pe_nodes:
            dft_nodes = self.nodes - eri_nodes - pe_nodes
        else:
            pe_nodes = self.nodes - eri_nodes - dft_nodes

        node_grps = [0] * eri_nodes + [1] * dft_nodes + [2] * pe_nodes
        eri_comm = (node_grps[self.rank] == 0)
        dft_comm = (node_grps[self.rank] == 1)
        pe_comm = (node_grps[self.rank] == 2)

        subcomms = SubCommunicators(self.comm, node_grps)
        local_comm = subcomms.local_comm
        cross_comm = subcomms.cross_comm

        # calculate Fock on ERI nodes
        if eri_comm:
            t0 = tm.time()
            eri_drv = ElectronRepulsionIntegralsDriver(local_comm)
            eri_drv.compute(fock, dens, molecule, basis, screening)
            if dft:
                if not xcfun.is_hybrid():
                    for ifock in range(fock.number_of_fock_matrices()):
                        fock.scale(2.0, ifock)
            dt = tm.time() - t0

        # reset molecular grid for DFT and V_es for PE
        if self.rank != mpi_master():
            molgrid = MolecularGrid()
            V_es = np.zeros(0)
        if dft:
            if local_comm.Get_rank() == mpi_master():
                molgrid.broadcast(cross_comm.Get_rank(), cross_comm)
        if pe:
            if local_comm.Get_rank() == mpi_master():
                V_es = cross_comm.bcast(V_es, root=mpi_master())

        # calculate Fxc on DFT nodes
        if dft_comm:
            t0 = tm.time()
            xc_drv = XCIntegrator(local_comm)
            molgrid.distribute(local_comm.Get_rank(), local_comm.Get_size(),
                               local_comm)
            xc_drv.integrate(fock, dens, gs_density, molecule, basis, molgrid,
                             xcfun.get_func_label())
        if dft_comm:
            dt = tm.time() - t0

        # calculate e_pe and V_pe on PE nodes
        if pe_comm:
            from .polembed import PolEmbed
            t0 = tm.time()
            pe_drv = PolEmbed(molecule, basis, local_comm, potfile)
            pe_drv.V_es = V_es.copy()
            for ifock in range(fock.number_of_fock_matrices()):
                dm = dens.alpha_to_numpy(ifock) + dens.beta_to_numpy(ifock)
                e_pe, V_pe = pe_drv.get_pe_contribution(dm, elec_only=True)
                if local_comm.Get_rank() == mpi_master():
                    fock.add_matrix(DenseMatrix(V_pe), ifock)
            dt = tm.time() - t0

        # collect Fock on master node
        fock.reduce_sum(self.rank, self.nodes, self.comm)

        dt = self.comm.gather(dt, root=mpi_master())
        if self.rank == mpi_master():
            time_eri = sum(dt[:eri_nodes])
            time_dft = sum(dt[eri_nodes:eri_nodes + dft_nodes])
            time_pe = sum(dt[eri_nodes + dft_nodes:])
            time_sum = time_eri + time_dft + time_pe
            self.split_comm_ratio = [
                time_eri / time_sum,
                time_dft / time_sum,
                time_pe / time_sum,
            ]
        self.split_comm_ratio = self.comm.bcast(self.split_comm_ratio,
                                                root=mpi_master())


def get_rhs(operator, components, molecule, basis, scf_tensors, rank, comm):
    """
    Creates right-hand side of linear response equations.

    :param operator:
        The string for the operator.
    :param components:
        The string for Cartesian components.
    :param molecule:
        The molecule.
    :param basis:
        The AO basis set.
    :param scf_tensors:
        The dictionary of tensors from converged SCF wavefunction.
    :param rank:
        Rank of the MPI process.
    :param comm:
        The MPI communicator.

    :return:
        The right-hand sides (gradients).
    """

    # compute 1e integral

    assert_msg_critical(
        operator in [
            'dipole', 'linear_momentum', 'linear momentum', 'angular_momentum',
            'angular momentum'
        ], 'get_rhs: unsupported operator {}'.format(operator))

    if operator == 'dipole':
        dipole_drv = ElectricDipoleIntegralsDriver(comm)
        dipole_mats = dipole_drv.compute(molecule, basis)

        if rank == mpi_master():
            integrals = (dipole_mats.x_to_numpy(), dipole_mats.y_to_numpy(),
                         dipole_mats.z_to_numpy())
        else:
            integrals = tuple()

    elif operator in ['linear_momentum', 'linear momentum']:
        linmom_drv = LinearMomentumIntegralsDriver(comm)
        linmom_mats = linmom_drv.compute(molecule, basis)

        if rank == mpi_master():
            integrals = (linmom_mats.x_to_numpy(), linmom_mats.y_to_numpy(),
                         linmom_mats.z_to_numpy())
        else:
            integrals = tuple()

    elif operator in ['angular_momentum', 'angular momentum']:
        angmom_drv = AngularMomentumIntegralsDriver(comm)
        angmom_mats = angmom_drv.compute(molecule, basis)

        if rank == mpi_master():
            integrals = (angmom_mats.x_to_numpy(), angmom_mats.y_to_numpy(),
                         angmom_mats.z_to_numpy())
        else:
            integrals = tuple()

    # compute right-hand side

    if rank == mpi_master():
        indices = {'x': 0, 'y': 1, 'z': 2}
        integral_comps = [integrals[indices[p]] for p in components]

        mo = scf_tensors['C']
        S = scf_tensors['S']
        D = scf_tensors['D'][0] + scf_tensors['D'][1]

        nocc = molecule.number_of_alpha_electrons()
        norb = mo.shape[1]

        matrices = tuple(
            np.linalg.multi_dot([
                mo.T,
                (np.linalg.multi_dot([S, D, P.T]) -
                 np.linalg.multi_dot([P.T, D, S])), mo
            ]) for P in integral_comps)

        gradients = tuple(lrmat2vec(m, nocc, norb) for m in matrices)
        return gradients

    else:
        return tuple()


def get_complex_rhs(operator, components, molecule, basis, scf_tensors, rank,
                    comm):
    """
    Creates right-hand side of linear response equations.

    :param operator:
        The string for the operator.
    :param components:
        The string for Cartesian components.
    :param molecule:
        The molecule.
    :param basis:
        The AO basis set.
    :param scf_tensors:
        The dictionary of tensors from converged SCF wavefunction.
    :param rank:
        Rank of the MPI process.
    :param comm:
        The MPI communicator.

    :return:
        The right-hand sides (gradients).
    """

    # compute 1e integral

    assert_msg_critical(
        operator in [
            'dipole', 'linear_momentum', 'linear momentum', 'angular_momentum',
            'angular momentum'
        ], 'get_rhs: unsupported operator {}'.format(operator))

    if operator == 'dipole':
        dipole_drv = ElectricDipoleIntegralsDriver(comm)
        dipole_mats = dipole_drv.compute(molecule, basis)

        if rank == mpi_master():
            integrals = (dipole_mats.x_to_numpy(), dipole_mats.y_to_numpy(),
                         dipole_mats.z_to_numpy())
        else:
            integrals = tuple()

    elif operator in ['linear_momentum', 'linear momentum']:
        linmom_drv = LinearMomentumIntegralsDriver(comm)
        linmom_mats = linmom_drv.compute(molecule, basis)

        if rank == mpi_master():
            integrals = (-1j * linmom_mats.x_to_numpy(),
                         -1j * linmom_mats.y_to_numpy(),
                         -1j * linmom_mats.z_to_numpy())
        else:
            integrals = tuple()

    elif operator in ['angular_momentum', 'angular momentum']:
        angmom_drv = AngularMomentumIntegralsDriver(comm)
        angmom_mats = angmom_drv.compute(molecule, basis)

        if rank == mpi_master():
            integrals = (-1j * angmom_mats.x_to_numpy(),
                         -1j * angmom_mats.y_to_numpy(),
                         -1j * angmom_mats.z_to_numpy())
        else:
            integrals = tuple()

    # compute right-hand side

    if rank == mpi_master():
        indices = {'x': 0, 'y': 1, 'z': 2}
        integral_comps = [integrals[indices[p]] for p in components]

        mo = scf_tensors['C']
        S = scf_tensors['S']
        D = scf_tensors['D'][0] + scf_tensors['D'][1]

        nocc = molecule.number_of_alpha_electrons()
        norb = mo.shape[1]

        matrices = tuple(
            np.linalg.multi_dot([
                mo.T,
                (np.linalg.multi_dot([S, D, P.conj().T]) -
                 np.linalg.multi_dot([P.conj().T, D, S])), mo
            ]) for P in integral_comps)

        gradients = tuple(lrmat2vec(m, nocc, norb) for m in matrices)
        return gradients

    else:
        return tuple()


def lrvec2mat(vec, nocc, norb):
    """
    Converts vectors to matrices.

    :param vec:
        The vectors.
    :param nocc:
        Number of occupied orbitals.
    :param norb:
        Number of orbitals.

    :return:
        The matrices.
    """

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


def lrmat2vec(mat, nocc, norb):
    """
    Converts matrices to vectors.

    :param mat:
        The matrices.
    :param nocc:
        Number of occupied orbitals.
    :param norb:
        Number of orbitals.

    :return:
        The vectors.
    """

    xv = ExcitationVector(szblock.aa, 0, nocc, nocc, norb, True)
    excitations = list(
        itertools.product(xv.bra_unique_indexes(), xv.ket_unique_indexes()))

    z = [mat[i, j] for i, j in excitations]
    y = [mat[j, i] for i, j in excitations]
    return np.array(z + y)


def remove_linear_dependence(basis, threshold):
    """
    Removes linear dependence in a set of vectors.

    :param basis:
        The set of vectors.
    :param threshold:
        The threshold for removing linear dependence.

    :return:
        The new set of vectors.
    """

    Sb = np.matmul(basis.T, basis)
    l, T = np.linalg.eigh(Sb)
    b_norm = np.sqrt(Sb.diagonal())
    mask = l > b_norm * threshold
    return np.matmul(basis, T[:, mask])


def remove_linear_dependence_half(basis, threshold):
    """
    Removes linear dependence in a set of symmetrized vectors.

    :param basis:
        The set of upper parts of symmetrized vectors.
    :param threshold:
        The threshold for removing linear dependence.

    :return:
        The new set of vectors.
    """

    Sb = 2 * np.matmul(basis.T, basis)
    l, T = np.linalg.eigh(Sb)
    b_norm = np.sqrt(Sb.diagonal())
    mask = l > b_norm * threshold
    return np.matmul(basis, T[:, mask])


def orthogonalize_gram_schmidt(tvecs):
    """
    Applies modified Gram Schmidt orthogonalization to trial vectors.

    :param tvecs:
        The trial vectors.

    :return:
        The orthogonalized trial vectors.
    """

    if tvecs.shape[1] > 0:

        f = 1.0 / np.linalg.norm(tvecs[:, 0])
        tvecs[:, 0] *= f

        for i in range(1, tvecs.shape[1]):
            for j in range(i):
                f = np.dot(tvecs[:, i], tvecs[:, j]) / np.dot(
                    tvecs[:, j], tvecs[:, j])
                tvecs[:, i] -= f * tvecs[:, j]
            f = 1.0 / np.linalg.norm(tvecs[:, i])
            tvecs[:, i] *= f

    return tvecs


def orthogonalize_gram_schmidt_half(tvecs):
    """
    Applies modified Gram Schmidt orthogonalization to trial vectors.

    :param tvecs:
        The trial vectors.

    :return:
        The orthogonalized trial vectors.
    """

    invsqrt2 = np.sqrt(2.0)

    if tvecs.shape[1] > 0:

        f = invsqrt2 / np.linalg.norm(tvecs[:, 0])
        tvecs[:, 0] *= f

        for i in range(1, tvecs.shape[1]):
            for j in range(i):
                f = np.dot(tvecs[:, i], tvecs[:, j]) / np.dot(
                    tvecs[:, j], tvecs[:, j])
                tvecs[:, i] -= f * tvecs[:, j]
            f = invsqrt2 / np.linalg.norm(tvecs[:, i])
            tvecs[:, i] *= f

    return tvecs


def normalize(vecs):
    """
    Normalizes vectors by dividing by vector norm.

    :param vecs:
        The vectors.

    :param Retruns:
        The normalized vectors.
    """

    if len(vecs.shape) != 1:
        for vec in range(vecs.shape[1]):
            invnorm = 1.0 / np.linalg.norm(vecs[:, vec])
            vecs[:, vec] *= invnorm
    else:
        invnorm = 1.0 / np.linalg.norm(vecs)
        vecs *= invnorm

    return vecs


def normalize_half(vecs):
    """
    Normalizes half-sized vectors by dividing by vector norm.

    :param vecs:
        The half-sized vectors.

    :param Retruns:
        The normalized vectors.
    """

    invsqrt2 = np.sqrt(2.0)

    if len(vecs.shape) != 1:
        for vec in range(vecs.shape[1]):
            invnorm = invsqrt2 / np.linalg.norm(vecs[:, vec])
            vecs[:, vec] *= invnorm
    else:
        invnorm = invsqrt2 / np.linalg.norm(vecs)
        vecs *= invnorm

    return vecs


def construct_ed_sd(orb_ene, nocc, norb):
    """
    Gets the E0 and S0 diagonal elements as arrays.

    :param orb_ene:
        Orbital energies.
    :param nocc:
        Number of occupied orbitals.
    :param norb:
        Number of orbitals.

    :return:
        The E0 and S0 diagonal elements as numpy arrays.
    """

    xv = ExcitationVector(szblock.aa, 0, nocc, nocc, norb, True)
    excitations = list(
        itertools.product(xv.bra_unique_indexes(), xv.ket_unique_indexes()))

    z = [2.0 * (orb_ene[j] - orb_ene[i]) for i, j in excitations]
    ediag = np.array(z + z)

    lz = len(excitations)
    sdiag = 2.0 * np.ones(2 * lz)
    sdiag[lz:] = -2.0

    return ediag, sdiag


def construct_ed_sd_half(orb_ene, nocc, norb):
    """
    Gets the upper half of E0 and S0 diagonal elements as arrays.

    :param orb_ene:
        Orbital energies.
    :param nocc:
        Number of occupied orbitals.
    :param norb:
        Number of orbitals.

    :return:
        The upper half of E0 and S0 diagonal elements as numpy arrays.
    """

    xv = ExcitationVector(szblock.aa, 0, nocc, nocc, norb, True)
    excitations = list(
        itertools.product(xv.bra_unique_indexes(), xv.ket_unique_indexes()))

    z = [2.0 * (orb_ene[j] - orb_ene[i]) for i, j in excitations]
    ediag = np.array(z)

    lz = len(excitations)
    sdiag = 2.0 * np.ones(lz)

    return ediag, sdiag


def write_rsp_hdf5(fname, arrays, labels, e_nuc, nuclear_charges, basis_set,
                   dft_func_label, ostream):
    """
    Writes response vectors to checkpoint file. Nuclear charges and basis
    set can also be written to the checkpoint file.

    :param fname:
        Name of the checkpoint file.
    :param arrays:
        The response vectors.
    :param labels:
        The list of labels for trial vecotrs and transformed vectors.
    :param e_nuc:
        Nuclear repulsion energy.
    :param nuclear_charges:
        Nuclear charges of the molecule.
    :param basis_set:
        Name of the AO basis set.
    :param dft_func_label:
        Name of the density functional.
    :param ostream:
        The output stream.
    """

    valid_checkpoint = (fname and isinstance(fname, str))

    if not valid_checkpoint:
        return

    hf = h5py.File(fname, 'w')

    for label, array in zip(labels, arrays):
        hf.create_dataset(label, data=array, compression="gzip")

    hf.create_dataset('nuclear_repulsion',
                      data=np.array([e_nuc]),
                      compression='gzip')

    hf.create_dataset('nuclear_charges',
                      data=nuclear_charges,
                      compression='gzip')

    hf.create_dataset('basis_set',
                      data=np.string_([basis_set]),
                      compression='gzip')

    hf.create_dataset('dft_func_label',
                      data=np.string_([dft_func_label]),
                      compression='gzip')

    hf.close()

    checkpoint_text = 'Checkpoint written to file: '
    checkpoint_text += fname
    ostream.print_info(checkpoint_text)
    ostream.print_blank()


def read_rsp_hdf5(fname, labels, e_nuc, nuclear_charges, basis_set,
                  dft_func_label, ostream):
    """
    Reads response vectors from checkpoint file. Nuclear charges and basis
    set will be used to validate the checkpoint file.

    :param fname:
        Name of the checkpoint file.
    :param labels:
        The list of labels for trial vecotrs and transformed vectors.
    :param e_nuc:
        Nuclear repulsion energy.
    :param nuclear_charges:
        Nuclear charges of the molecule.
    :param basis_set:
        Name of the AO basis set.
    :param dft_func_label:
        Name of the density functional.
    :param ostream:
        The output stream.

    :return:
        The trial vectors and transformed vectors.
    """

    valid_checkpoint = (fname and isinstance(fname, str) and isfile(fname))

    if not valid_checkpoint:
        return tuple([None] * len(labels))

    hf = h5py.File(fname, 'r')

    match_nuclear_repulsion = False
    if 'nuclear_repulsion' in hf:
        h5_e_nuc = np.array(hf.get('nuclear_repulsion'))[0]
        if h5_e_nuc > 0.0 and e_nuc > 0.0:
            match_nuclear_repulsion = (abs(1.0 - h5_e_nuc / e_nuc) < 1.0e-13)
        else:
            match_nuclear_repulsion = (h5_e_nuc == e_nuc)

    match_nuclear_charges = False
    if 'nuclear_charges' in hf:
        h5_nuclear_charges = np.array(hf.get('nuclear_charges'))
        if h5_nuclear_charges.shape == nuclear_charges.shape:
            match_nuclear_charges = (
                h5_nuclear_charges == nuclear_charges).all()

    match_basis_set = False
    if 'basis_set' in hf:
        h5_basis_set = hf.get('basis_set')[0].decode('utf-8')
        match_basis_set = (h5_basis_set.upper() == basis_set.upper())

    match_dft_func = False
    if 'dft_func_label' in hf:
        h5_func_label = hf.get('dft_func_label')[0].decode('utf-8')
        match_dft_func = (h5_func_label.upper() == dft_func_label.upper())

    arrays = [None] * len(labels)

    if (match_nuclear_repulsion and match_nuclear_charges and
            match_basis_set and match_dft_func):
        for i in range(len(labels)):
            if labels[i] in hf.keys():
                arrays[i] = np.array(hf.get(labels[i]))

    hf.close()

    is_empty = [a is None for a in arrays]
    if True not in is_empty:
        checkpoint_text = 'Restarting from checkpoint file: '
        checkpoint_text += fname
        ostream.print_info(checkpoint_text)
        ostream.print_blank()

    return tuple(arrays)
