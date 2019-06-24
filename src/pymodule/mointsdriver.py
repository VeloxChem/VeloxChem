import time as tm

from .veloxchemlib import ElectronRepulsionIntegralsDriver
from .veloxchemlib import mpi_master
from .veloxchemlib import fockmat
from .veloxchemlib import moints
from .veloxchemlib import MOIntsBatch
from .veloxchemlib import TwoIndexes
from .aofockmatrix import AOFockMatrix
from .aodensitymatrix import AODensityMatrix
from .subcommunicators import SubCommunicators
from .qqscheme import get_qq_type
from .qqscheme import get_qq_scheme


class MOIntegralsDriver:
    """Implements MO integrals driver.

    Implements MO integrals driver.

    Attributes
    ----------
    qq_type
        The electron repulsion integrals screening scheme.
    eri_thresh
        The electron repulsion integrals screening threshold.
    num_matrices
        Number of Fock matrices to be computed.
    batch_size
        Batch size for computing Fock matrices.
    comm
        The MPI communicator.
    rank
        The MPI rank.
    nodes
        Number of MPI processes.
    ostream
        The output stream.
    """

    def __init__(self, comm, ostream):
        """Initializes MO integrals driver.

        Initializes MO integrals driver  to default setup.

        Parameters
        ----------
        comm
            The MPI communicator.
        ostream
            The output stream.
        """

        # screening scheme
        self.qq_type = "QQ_DEN"
        self.eri_thresh = 1.0e-12

        # Fock matrices
        self.num_matrices = 0
        self.batch_size = 3000

        # mpi information
        self.comm = comm
        self.rank = self.comm.Get_rank()
        self.nodes = self.comm.Get_size()

        # output stream
        self.ostream = ostream

    def compute(self, molecule, ao_basis, mol_orbs, mints_type, grps):
        """Performs MO integrals calculation.

        Performs MO integrals calculation for a molecule and a basis set.

        Parameters
        ----------
        molecule
            The molecule.
        ao_basis
            The AO basis set.
        mol_orbs
            The molecular orbitals.
        mints_type
            The type of MO integrals to be calculated.
        grps
            The color group for creating MPI subcommunicators.

        Returns
        -------
        MOIntegralsBatch
            The computed MO integrals batch.
        """

        # start timer
        start_time = tm.time()

        # split communicators
        subcomm = SubCommunicators(self.comm, grps)
        local_comm = subcomm.local_comm
        cross_comm = subcomm.cross_comm

        global_master = (self.rank == mpi_master())

        local_rank = local_comm.Get_rank()
        local_nodes = local_comm.Get_size()
        local_master = (local_rank == mpi_master())

        cross_rank = cross_comm.Get_rank()
        cross_nodes = cross_comm.Get_size()

        # set up screening data
        eri_drv = ElectronRepulsionIntegralsDriver(local_comm)

        qq_data = eri_drv.compute(get_qq_scheme(self.qq_type), self.eri_thresh,
                                  molecule, ao_basis)

        # initialize MO integrals batch
        moints_batch = MOIntsBatch()

        # set up properties of MO integrals batch
        moints_batch.set_ext_indexes(self.get_external_indexes(mints_type))
        moints_batch.set_batch_type(self.get_moints_type(mints_type))

        # broadcast molecular orbitals to local master nodes
        if local_master:
            mol_orbs.broadcast(cross_rank, cross_comm)

            # transformation matrices
            xmat, ymat = self.get_transformation_vectors(
                mints_type, mol_orbs, molecule)

            # set up bra and ket orbitals list
            bra_ids, ket_ids = self.set_orbital_pairs(mints_type, mol_orbs,
                                                      molecule)

            # create MO integrals batch
            batch_pos, batch_dim = self.get_batch_dimensions(bra_ids)

            loc_bra_ids = [
                bra_ids[pos + cross_rank:pos + dim:cross_nodes]
                for pos, dim in zip(batch_pos, batch_dim)
            ]
            loc_ket_ids = [
                ket_ids[pos + cross_rank:pos + dim:cross_nodes]
                for pos, dim in zip(batch_pos, batch_dim)
            ]

        else:
            loc_bra_ids = []
            loc_ket_ids = []

        loc_bra_ids = local_comm.bcast(loc_bra_ids, root=mpi_master())
        loc_ket_ids = local_comm.bcast(loc_ket_ids, root=mpi_master())

        # print title
        if global_master:
            self.print_header()

        # loop over batches
        for cur_bra_ids, cur_ket_ids in zip(loc_bra_ids, loc_ket_ids):

            # compute pair densities on local master nodes
            if local_master:
                pair_den = mol_orbs.get_pair_density(cur_bra_ids, cur_ket_ids)
            else:
                pair_den = AODensityMatrix()

            # broadcast pair densities via local communicators
            pair_den.broadcast(local_rank, local_comm)

            # compute multiple Fock matrices
            fock_mat = AOFockMatrix(pair_den)

            self.set_fock_matrices_type(mints_type, fock_mat)
            eri_drv.compute(fock_mat, pair_den, molecule, ao_basis, qq_data)

            fock_mat.reduce_sum(local_rank, local_nodes, local_comm)

            if local_master:
                # transform Fock matrices and add MO integrals to batch
                moints_batch.append(fock_mat, xmat, ymat, cur_bra_ids,
                                    cur_ket_ids)

        if global_master:
            self.print_finish(start_time)

        return moints_batch

    def collect_moints_batches(self, moints_batch, grps):
        """Collects MO integrals.

        Collects MO integrals to the master node.

        Parameters
        ----------
        moints_batch
            The MO integrals batch.
        grps
            The color group for creating MPI subcommunicators.

        Returns
        -------
        MOIntegralsBatch
            The collected MO integrals batch.
        """

        subcomm = SubCommunicators(self.comm, grps)
        local_comm = subcomm.local_comm
        cross_comm = subcomm.cross_comm

        global_master = (self.rank == mpi_master())
        local_master = (local_comm.Get_rank() == mpi_master())

        cross_rank = cross_comm.Get_rank()
        cross_nodes = cross_comm.Get_size()

        if local_master:
            moints_batch.collect_batches(cross_rank, cross_nodes, cross_comm)

        if not global_master:
            moints_batch = MOIntsBatch()

        return moints_batch

    def print_header(self):
        """Prints header for the MO integrals driver.

        Prints header for the MO integrals driver.
        """

        self.ostream.print_blank()
        self.ostream.print_header(
            "Integrals Transformation (AO->MO) Driver Setup")
        self.ostream.print_header(46 * "=")
        self.ostream.print_blank()

        str_width = 80
        cur_str = "Number of Fock matrices      : " + str(self.num_matrices)
        self.ostream.print_header(cur_str.ljust(str_width))
        cur_str = "Size of Fock matrices batch  : " + str(self.batch_size)
        self.ostream.print_header(cur_str.ljust(str_width))
        cur_str = "ERI screening scheme         : " + get_qq_type(self.qq_type)
        self.ostream.print_header(cur_str.ljust(str_width))
        cur_str = "ERI Screening Threshold      : " + \
            "{:.1e}".format(self.eri_thresh)
        self.ostream.print_header(cur_str.ljust(str_width))
        self.ostream.print_blank()
        self.ostream.flush()

    def print_finish(self, start_time):
        """Prints header for the MO integrals driver.

        Prints header for the MO integrals driver.

        Parameters
        ----------
        start_time
            The starting time of the calculation.
        """

        valstr = "*** Integrals transformation (AO->MO) done in "
        valstr += "{:.2f}".format(tm.time() - start_time) + " sec."
        self.ostream.print_blank()
        self.ostream.print_header(valstr.ljust(92))
        self.ostream.print_blank()

    def get_batch_dimensions(self, vec_ids):
        """Gets the dimension of the batches.

        Gets the dimension of the batches.

        Parameters
        ----------
        vec_ids
            The list of vector id.

        Returns
        -------
        tuple
            The tuple containing the position and the dimension of the batch.
        """

        vdim = len(vec_ids)

        # determine number of batches
        nbatch = vdim // self.batch_size
        if (vdim % self.batch_size) != 0:
            nbatch = nbatch + 1

        # compute dimensions of batches
        bdim = vdim // nbatch
        brem = vdim % nbatch

        # populate dimensions list
        vec_dim = []
        for i in range(nbatch):
            if i < brem:
                vec_dim.append(bdim + 1)
            else:
                vec_dim.append(bdim)

        # populate positions list
        vec_pos = []
        cur_pos = 0
        for i in range(nbatch):
            vec_pos.append(cur_pos)
            cur_pos = cur_pos + vec_dim[i]

        return (vec_pos, vec_dim)

    def get_num_orbitals(self, mol_orbs, molecule):
        """Gets number of orbitals.

        Gets number of occupied and virtual orbitals.

        Parameters
        ----------
        mol_orbs
            The molecular orbitals.
        molecule
            The molecule.

        Returns
        -------
        tuple
            The number of occupied and virtual orbitals.
        """

        nocc = molecule.number_of_alpha_electrons()
        nvirt = mol_orbs.number_mos() - nocc

        return (nocc, nvirt)

    def set_orbital_pairs(self, mints_type, mol_orbs, molecule):
        """Sets up bra and ket orbitals list.

        Sets up bra and ket orbitals list.

        Parameters
        ----------
        mints_type
            Type of MO integrals.
        mol_orbs
            The molecular orbitals.
        molecule
            The molecule.

        Returns
        -------
        tuple
            The tuple containing the lists of bra and ket ids.
        """

        nocc, nvirt = self.get_num_orbitals(mol_orbs, molecule)

        # cases: oooo, ooov, oovv, ovov, asym_oooo, asym_ooov, asym_oovv
        bra_dim = (0, nocc)
        ket_dim = (0, nocc)

        # case: ovvv
        if mints_type == "OVVV":
            ket_dim = (nocc, nocc + nvirt)

        # case: vvvv
        if mints_type == "VVVV":
            bra_dim = (nocc, nocc + nvirt)
            ket_dim = (nocc, nocc + nvirt)

        # case: asym_ovov
        if mints_type == "ASYM_OVOV":
            ket_dim = (nocc, nocc + nvirt)

        # case: asym_ovvv
        if mints_type == "ASYM_OVVV":
            ket_dim = (nocc, nocc + nvirt)

        # case: asym_vvvv
        if mints_type == "ASYM_VVVV":
            bra_dim = (nocc, nocc + nvirt)
            ket_dim = (nocc, nocc + nvirt)

        usesym = self.use_symmetry(mints_type)

        # set up list of orbital pairs
        bra_ids = []
        ket_ids = []
        for i in range(bra_dim[0], bra_dim[1]):
            for j in range(ket_dim[0], ket_dim[1]):
                if (i > j) and usesym:
                    continue
                bra_ids.append(i)
                ket_ids.append(j)

        # set number of Fock matrices
        self.num_matrices = len(bra_ids)

        return (bra_ids, ket_ids)

    def get_transformation_vectors(self, mints_type, mol_orbs, molecule):
        """Gets transformation matrices.

        Gets transformation matrices.

        Parameters
        ----------
        mints_type
            Type of MO integrals.
        mol_orbs
            The molecular orbitals.
        molecule
            The molecule.

        Returns
        -------
        tuple
            The tuple containing transformation matrices.
        """

        nocc, nvirt = self.get_num_orbitals(mol_orbs, molecule)

        # cases: oovv, ovov, ovvv, vvvv, asym_oovv, asym_ovvv, asym_vvvv
        xmat = mol_orbs.alpha_orbitals(nocc, nvirt)
        ymat = mol_orbs.alpha_orbitals(nocc, nvirt)

        # case: oooo
        if mints_type == "OOOO":
            xmat = mol_orbs.alpha_orbitals(0, nocc)
            ymat = mol_orbs.alpha_orbitals(0, nocc)

        # case: ooov
        if mints_type == "OOOV":
            xmat = mol_orbs.alpha_orbitals(0, nocc)

        # case: asym_oooo
        if mints_type == "ASYM_OOOO":
            xmat = mol_orbs.alpha_orbitals(0, nocc)
            ymat = mol_orbs.alpha_orbitals(0, nocc)

        # case: asym_ooov
        if mints_type == "ASYM_OOOV":
            xmat = mol_orbs.alpha_orbitals(0, nocc)

        # case: asym_ovov
        if mints_type == "ASYM_OVOV":
            xmat = mol_orbs.alpha_orbitals(0, nocc)

        return (xmat, ymat)

    def set_fock_matrices_type(self, mints_type, fock_matrices):
        """Sets Fock matrix types.

        Sets Fock matrix types.

        Parameters
        ----------
        mints_type
            Type of MO integrals.
        fock_matrices
            The Fock matrices.
        """

        if mints_type.startswith("ASYM"):
            nfock = fock_matrices.number_of_fock_matrices()
            for i in range(nfock):
                fock_matrices.set_fock_type(fockmat.rgenk, i)
            return

        if mints_type != "OOVV":
            nfock = fock_matrices.number_of_fock_matrices()
            for i in range(nfock):
                fock_matrices.set_fock_type(fockmat.rgenj, i)

    def get_external_indexes(self, mints_type):
        """Gets external indices.

        Gets external indices.

        Parameters
        ----------
        mints_type
            Type of MO integrals.

        Returns
        -------
        TwoIndexes
            The external indices.
        """

        if mints_type == "OOVV":
            return TwoIndexes(2, 3)

        if mints_type == "ASYM_OOOO":
            return TwoIndexes(2, 3)

        if mints_type == "ASYM_OOOV":
            return TwoIndexes(2, 3)

        if mints_type == "ASYM_OVOV":
            return TwoIndexes(2, 3)

        if mints_type == "ASYM_OOVV":
            return TwoIndexes(2, 3)

        if mints_type == "ASYM_OVVV":
            return TwoIndexes(2, 3)

        if mints_type == "ASYM_VVVV":
            return TwoIndexes(2, 3)

        return TwoIndexes(1, 3)

    def get_moints_type(self, mints_type):
        """Gets MO integrals type.

        Gets MO integrals type.

        Parameters
        ----------
        mints_type
            Type of MO integrals.

        Returns
        -------
        moints
            The MO integrals type.
        """

        if mints_type == "OOOO":
            return moints.oooo

        if mints_type == "OOOV":
            return moints.ooov

        if mints_type == "OOVV":
            return moints.oovv

        if mints_type == "OVOV":
            return moints.ovov

        if mints_type == "OVVV":
            return moints.ovvv

        if mints_type == "VVVV":
            return moints.vvvv

        if mints_type == "ASYM_OOOO":
            return moints.asym_oooo

        if mints_type == "ASYM_OOOV":
            return moints.asym_ooov

        if mints_type == "ASYM_OVOV":
            return moints.asym_ovov

        if mints_type == "ASYM_OOVV":
            return moints.asym_oovv

        if mints_type == "ASYM_OVVV":
            return moints.asym_ovvv

        if mints_type == "ASYM_VVVV":
            return moints.asym_vvvv

        return None

    def use_symmetry(self, mints_type):
        """Checks if symmetry is used.

        Checks if symmetry is used.

        Parameters
        ----------
        mints_type
            Type of MO integrals.

        Returns
        -------
        bool
            Whether symmetry is used.
        """

        if mints_type == "ASYM_OOOO":
            return True

        if mints_type == "ASYM_OOOV":
            return True

        if mints_type == "ASYM_OOVV":
            return True

        if mints_type == "ASYM_VVVV":
            return True

        return False
