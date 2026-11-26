from pathlib import Path

import numpy as np
import pytest

from mpi4py import MPI

from veloxchem import MolecularBasis, Molecule
from veloxchem.veloxchemlib import (PackedMatrix, mat_t, mpi_master, rimode,
                                    make_matrix, FockGeom1000Driver,
                                    SimdRIJFockDriver, SimdRIJGradientDriver)


@pytest.mark.solvers
class TestRIJFockGeomGradDriver:

    def get_data_h2o_trimer(self):

        h2ostr = """
            O    1.206898105612   1.227024188072  -0.113458607566
            H    0.490554304233   1.138628097468   0.508036765414
            H    1.992568601678   1.153431301285   0.406410355109
            O   -0.995720200000   0.016041500000   1.242255600000
            H   -1.454270300000  -0.566974100000   1.847281700000
            H   -0.937795000000  -0.481791200000   0.426756200000
            O   -0.243234300000  -1.019856600000  -1.195380800000
            H    0.436753600000  -0.375943300000  -0.997329700000
            H   -0.503183500000  -0.825149200000  -2.095795900000
        """

        mol = Molecule.read_str(h2ostr, 'angstrom')
        bas = MolecularBasis.read(mol, 'def2-svp')
        aux_bas = MolecularBasis.read(mol, 'def2-universal-jfit')

        return mol, bas, aux_bas

    def test_h2o_trimer_ri_j_gradient(self):

        mol, bas, aux_bas = self.get_data_h2o_trimer()

        natoms = mol.number_of_atoms()
        nao = bas.get_dimensions_of_basis()

        here = Path(__file__).parent
        density = np.load(here / 'data' / 'h2o.trimer.svp.density.npy')

        den_mat = make_matrix(bas, mat_t.symmetric)
        den_mat.set_values(density)

        # four-centre Coulomb gradient for reference

        fock_drv = FockGeom1000Driver()
        ref_grad = np.zeros((natoms, 3))

        for iatom in range(natoms):
            fmats = fock_drv.compute(bas, mol, den_mat, iatom, '2jkx', 0.0, 0.0)
            for xyz, label in enumerate('XYZ'):
                fmat = fmats.matrix(label).full_matrix().to_numpy()
                ref_grad[iatom, xyz] = np.trace(np.matmul(fmat, density))

        # fitting coefficients of the density

        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        nodes = comm.Get_size()

        packed_density = PackedMatrix(nao, nao, mat_t.symmetric)
        packed_density.from_numpy(np.ascontiguousarray(density))

        ri_drv = SimdRIJFockDriver()

        metric = PackedMatrix()
        if rank == mpi_master():
            metric = ri_drv.make_metric(mol, aux_bas, 1.0e-12)
        metric = metric.broadcast(comm, root=mpi_master())

        ri_drv.prepare(mol, bas, aux_bas, 1.0e-12, 1 << 30, 1.0e-12,
                       rimode.in_memory, metric, rank, nodes)

        mine = ri_drv.owned_parts()

        local = np.array(ri_drv.compute_gamma(packed_density, mine),
                         dtype=np.float64)

        total = np.zeros_like(local)
        comm.Allreduce(local, total, op=MPI.SUM)

        gamma = ri_drv.solve_fitting(total)

        # SIMD RI-J Coulomb gradient of every atom

        grad_drv = SimdRIJGradientDriver(1.0e-12)

        all_atoms = list(range(natoms))
        local_atoms = mol.partition_atoms(comm)

        atomgrad = grad_drv.compute(mol, bas, aux_bas, gamma, packed_density,
                                    all_atoms, list(local_atoms),
                                    rank == mpi_master())

        calc_grad = comm.allreduce(atomgrad.to_numpy())

        assert np.allclose(ref_grad, calc_grad,
                           rtol=1.0e-4, atol=1.0e-4)
