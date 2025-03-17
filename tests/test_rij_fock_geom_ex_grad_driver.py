from mpi4py import MPI
import numpy as np

from veloxchem import Molecule, MolecularBasis
from veloxchem import FockGeom1000Driver
from veloxchem import RIFockDriver, TwoCenterElectronRepulsionDriver
from veloxchem import T4CScreener, SubMatrix
from veloxchem import make_matrix, mat_t, mpi_master
from veloxchem import ScfRestrictedDriver
from veloxchem.rigradientdriver import RIFockGradDriver


class TestRIJFockGeomExGradDriver:

    def get_data_h2o(self):

        h2ostr = """
           O   1.206898105612   1.227024188072  -0.113458607566
           H   0.490554304233   1.138628097468   0.508036765414
           H   1.992568601678   1.153431301285   0.406410355109
        """

        mol = Molecule.read_str(h2ostr, 'angstrom')
        bas = MolecularBasis.read(mol, 'def2-svp')
        aux_bas = MolecularBasis.read(mol, 'def2-universal-jfit')

        return mol, bas, aux_bas

    def test_h2o_ri_j_grad_screened(self):

        mol, bas, bas_aux = self.get_data_h2o()

        scf_drv = ScfRestrictedDriver()
        scf_drv.ostream.mute()
        scf_res = scf_drv.compute(mol, bas)

        comm = scf_drv.comm

        if scf_drv.rank == mpi_master():
            density = scf_res['D_alpha']

            mo = scf_res['C_alpha']
            nocc = mol.number_of_alpha_electrons()
            mo_occ = mo[:, :nocc]
            mo_vir = mo[:, nocc:]
            nvir = mo_vir.shape[1]

            exc_mat = np.zeros((nocc, nvir))
            exc_mat[nocc - 1, 0] = 1.0
            exc_density = np.linalg.multi_dot([mo_occ, exc_mat, mo_vir.T])
        else:
            density = None
            exc_density = None
        density = comm.bcast(density, root=mpi_master())
        exc_density = comm.bcast(exc_density, root=mpi_master())
        sym_exc_density = 0.5 * (exc_density + exc_density.T)

        gs_den_mat = make_matrix(bas, mat_t.symmetric)
        gs_den_mat.set_values(density)

        rw_den_sym_mat = make_matrix(bas, mat_t.symmetric)
        rw_den_sym_mat.set_values(sym_exc_density)

        rw_den_gen_mat = make_matrix(bas, mat_t.general)
        rw_den_gen_mat.set_values(exc_density)

        fock_grad_drv = FockGeom1000Driver()

        ket_t4c = T4CScreener()
        ket_t4c.partition(bas, mol, 'eri')

        # compute J metric
        t2c_drv = TwoCenterElectronRepulsionDriver()
        matj = t2c_drv.compute(mol, bas_aux)
        rmatj = np.linalg.inv(matj.full_matrix().to_numpy())
        invmatj = SubMatrix([0, 0, rmatj.shape[0], rmatj.shape[0]])
        invmatj.set_values(rmatj)

        local_atoms = mol.partition_atoms(comm)

        ri_fock_drv = RIFockDriver(invmatj)
        ri_fock_drv.prepare_buffers(mol, bas, bas_aux, local_atoms)

        local_gs_bq = np.array(ri_fock_drv.compute_local_bq_vector(gs_den_mat))
        gs_bq = np.zeros(local_gs_bq.shape)
        comm.Allreduce(local_gs_bq, gs_bq, op=MPI.SUM)

        local_rw_bq = np.array(
            ri_fock_drv.compute_local_bq_vector(rw_den_sym_mat))
        rw_bq = np.zeros(local_rw_bq.shape)
        comm.Allreduce(local_rw_bq, rw_bq, op=MPI.SUM)

        ri_grad_drv = RIFockGradDriver()

        natoms = mol.number_of_atoms()
        ref_grad = np.zeros((natoms, 3))
        calc_grad = np.zeros((natoms, 3))

        for iatom in local_atoms:

            bra_t4c = T4CScreener()
            bra_t4c.partition_atom(bas, mol, 'eri', iatom)

            rgrad = fock_grad_drv.compute(bas, bra_t4c, ket_t4c, gs_den_mat,
                                          rw_den_gen_mat, iatom, 'j', 0.0, 0.0,
                                          12)
            ref_grad[iatom, :] += np.array(rgrad)

            cgrad = ri_grad_drv.direct_compute(ket_t4c, bas, bas_aux, mol,
                                               rw_bq, gs_bq, rw_den_sym_mat,
                                               gs_den_mat, iatom, 12)
            calc_grad[iatom, :] += np.array(cgrad.coordinates())

        ref_grad = comm.allreduce(ref_grad)
        calc_grad = comm.allreduce(calc_grad)

        assert np.max(np.abs(ref_grad - calc_grad)) < 1.0e-4
