from pathlib import Path
import numpy as np
import pytest

from mpi4py import MPI

from veloxchem.veloxchemlib import mpi_master
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.molecule import Molecule
from veloxchem.fockdriver import FockDriver
from veloxchem.veloxchemlib import T4CScreener
from veloxchem.matrices import Matrices
from veloxchem.subcommunicators import SubCommunicators
from veloxchem.veloxchemlib import make_matrix
from veloxchem.veloxchemlib import mat_t


@pytest.mark.solvers
class TestFockDriverOnSubcomm:

    def get_data_h2o_dimer(self):

        h2o_dimer_str = """
            O -1.464  0.099  0.300
            H -1.956  0.624 -0.340
            H -1.797 -0.799  0.206
            O  1.369  0.146 -0.395
            H  1.894  0.486  0.335
            H  0.451  0.163 -0.083
        """
        mol = Molecule.read_str(h2o_dimer_str, 'angstrom')
        bas = MolecularBasis.read(mol, 'def2-svpd')

        return mol, bas

    def test_h2o_dimer_focks_svpd_with_subcomm(self):

        comm = MPI.COMM_WORLD

        t4c_drv = None
        num_densities = None
        den_mats = None

        if comm.Get_rank() == mpi_master():
            mol_h2o_dimer, bas_svpd = self.get_data_h2o_dimer()

            # load density matrix
            here = Path(__file__).parent
            D = np.load(str(here / 'data' / 'h2o.dimer.svpd.density.gen.npy'))

            dmat_0 = make_matrix(bas_svpd, mat_t.general)
            dmat_1 = make_matrix(bas_svpd, mat_t.general)

            dmat_0.set_values(D * 1.2)
            dmat_1.set_values(D * 0.8)

            den_mats = Matrices()
            den_mats.add(dmat_0, '0')
            den_mats.add(dmat_1, '1')

            num_densities = 2

            # screen basis function pairs
            t4c_drv = T4CScreener()
            t4c_drv.partition(bas_svpd, mol_h2o_dimer, 'eri')

        # broadcast data
        t4c_drv = comm.bcast(t4c_drv, root=mpi_master())
        num_densities = comm.bcast(num_densities, root=mpi_master())
        den_mats = Matrices.bcast(den_mats, comm, mpi_master())

        # compute Fock matrix
        fock_drv = FockDriver(comm)
        fock_mats = fock_drv.compute(t4c_drv, den_mats, ['j', 'k'], 0.0, 0.0,
                                     15)
        fock_mats = fock_mats.reduce(comm, mpi_master())

        node_grps = [p % 2 for p in range(comm.Get_size())]
        subcomm = SubCommunicators(comm, node_grps)
        local_comm = subcomm.local_comm
        # cross_comm = subcomm.cross_comm

        fock_mats_subcomm = fock_drv.compute_on_subcomm(local_comm, t4c_drv,
                                                        den_mats, ['j', 'k'],
                                                        0.0, 0.0, 15)
        fock_mats_subcomm = fock_mats_subcomm.reduce(local_comm, mpi_master())

        if comm.Get_rank() == mpi_master():
            ref_fock_0 = fock_mats.matrix('0').full_matrix().to_numpy()
            ref_fock_1 = fock_mats.matrix('1').full_matrix().to_numpy()
        else:
            ref_fock_0 = None
            ref_fock_1 = None

        ref_fock_0 = comm.bcast(ref_fock_0, root=mpi_master())
        ref_fock_1 = comm.bcast(ref_fock_1, root=mpi_master())

        for target_rank in [0, 1]:
            if comm.Get_rank() == target_rank:
                fock_subcomm_0 = fock_mats_subcomm.matrix(
                    '0').full_matrix().to_numpy()
                fock_subcomm_1 = fock_mats_subcomm.matrix(
                    '1').full_matrix().to_numpy()

                assert np.max(np.abs(fock_subcomm_0 - ref_fock_0)) < 1.0e-12
                assert np.max(np.abs(fock_subcomm_1 - ref_fock_1)) < 1.0e-12
