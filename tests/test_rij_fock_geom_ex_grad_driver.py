import numpy as np
import pytest

from mpi4py import MPI

from veloxchem import MolecularBasis, Molecule
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.veloxchemlib import (PackedMatrix, mat_t, mpi_master, rimode,
                                    make_matrix, FockGeom1000Driver,
                                    T4CScreener, SimdRIJFockDriver,
                                    SimdRIJGradientDriver)


@pytest.mark.solvers
class TestRIJFockGeomExGradDriver:

    def test_h2o_excited_state_ri_j_gradient(self):

        h2ostr = """
            O   1.206898105612   1.227024188072  -0.113458607566
            H   0.490554304233   1.138628097468   0.508036765414
            H   1.992568601678   1.153431301285   0.406410355109
        """

        mol = Molecule.read_str(h2ostr, 'angstrom')
        bas = MolecularBasis.read(mol, 'def2-svp')
        aux_bas = MolecularBasis.read(mol, 'def2-universal-jfit')

        scf_drv = ScfRestrictedDriver()
        scf_drv.ri = False
        scf_drv.ostream.mute()
        scf_results = scf_drv.compute(mol, bas)

        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        nodes = comm.Get_size()

        natoms = mol.number_of_atoms()
        nao = bas.get_dimensions_of_basis()

        if rank == mpi_master():
            density = scf_results['D_alpha']

            mo = scf_results['C_alpha']
            nocc = mol.number_of_alpha_occupied_orbitals(bas)
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

        rw_den_gen_mat = make_matrix(bas, mat_t.general)
        rw_den_gen_mat.set_values(exc_density)

        # four-centre reference for the ground state/excited state cross term

        fock_grad_drv = FockGeom1000Driver()

        ket_t4c = T4CScreener()
        ket_t4c.partition(bas, mol, 'eri')

        local_atoms = mol.partition_atoms(comm)

        ref_grad = np.zeros((natoms, 3))

        for iatom in local_atoms:
            bra_t4c = T4CScreener()
            bra_t4c.partition_atom(bas, mol, 'eri', iatom)

            rgrad = fock_grad_drv.compute(bas, bra_t4c, ket_t4c, gs_den_mat,
                                          rw_den_gen_mat, iatom, 'j', 0.0,
                                          0.0, 12)
            ref_grad[iatom, :] = np.array(rgrad)

        ref_grad = comm.allreduce(ref_grad)

        # fitting coefficients of the ground and excited state densities

        ri_drv = SimdRIJFockDriver()

        metric = PackedMatrix()
        if rank == mpi_master():
            metric = ri_drv.make_metric(mol, aux_bas, 1.0e-12)
        metric = metric.broadcast(comm, root=mpi_master())

        ri_drv.prepare(mol, bas, aux_bas, 1.0e-12, 1 << 30, 1.0e-12,
                       rimode.in_memory, metric, rank, nodes)

        def _fitting(density_np):
            packed = PackedMatrix(nao, nao, mat_t.symmetric)
            packed.from_numpy(np.ascontiguousarray(density_np))

            local = np.array(ri_drv.compute_gamma(packed, ri_drv.owned_parts()),
                             dtype=np.float64)

            total = np.zeros_like(local)
            comm.Allreduce(local, total, op=MPI.SUM)

            return np.array(ri_drv.solve_fitting(total), dtype=np.float64)

        gamma_gs = _fitting(density)
        gamma_ex = _fitting(sym_exc_density)

        # the cross term from the self terms the SIMD driver forms:
        # [G(gex + ggs, Dex + Dgs) - G(gex - ggs, Dex - Dgs)] / 8

        grad_drv = SimdRIJGradientDriver(1.0e-12)

        all_atoms = list(range(natoms))

        den_plus = PackedMatrix(nao, nao, mat_t.symmetric)
        den_plus.from_numpy(np.ascontiguousarray(sym_exc_density + density))

        den_minus = PackedMatrix(nao, nao, mat_t.symmetric)
        den_minus.from_numpy(np.ascontiguousarray(sym_exc_density - density))

        grad_plus = grad_drv.compute(mol, bas, aux_bas, gamma_ex + gamma_gs,
                                     den_plus, all_atoms, list(local_atoms),
                                     rank == mpi_master())

        grad_minus = grad_drv.compute(mol, bas, aux_bas, gamma_ex - gamma_gs,
                                      den_minus, all_atoms,
                                      list(local_atoms),
                                      rank == mpi_master())

        calc_grad = (grad_plus.to_numpy() - grad_minus.to_numpy()) / 8.0
        calc_grad = comm.allreduce(calc_grad)

        assert np.allclose(ref_grad, calc_grad,
                           rtol=1.0e-4, atol=1.0e-4)
