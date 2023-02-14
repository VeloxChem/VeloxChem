import numpy as np
import pytest

from veloxchem.veloxchemlib import GridDriver, XCNewIntegrator
from veloxchem.veloxchemlib import is_single_node, mpi_master, denmat
from veloxchem.aodensitymatrix import AODensityMatrix
from veloxchem.aofockmatrix import AOFockMatrix
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver


class TestExcVxc:

    @pytest.mark.skipif(not is_single_node(), reason='single node only')
    def test_slater(self):

        mol_str = """
            O  0.0000000000   0.0000000000  -0.0254395383
            H  0.0000000000   0.7695699584   0.5948147012
            H  0.0000000000  -0.7695699584   0.5948147012
            O  4.5000000000   0.0000000000  -0.0254395383
            H  4.5000000000   0.7695699584   0.5948147012
            H  4.5000000000  -0.7695699584   0.5948147012
        """
        basis_label = 'def2-svp'
        xcfun_label = 'slater'
        grid_level = 1
        tol = 1.0e-10

        molecule = Molecule.read_str(mol_str, units='angstrom')
        basis = MolecularBasis.read(molecule, basis_label, ostream=None)

        scf_drv = ScfRestrictedDriver()
        scf_drv.xcfun = xcfun_label
        scf_drv.grid_level = grid_level
        scf_drv.ostream.mute()
        scf_drv.compute(molecule, basis)
        gs_density = scf_drv.density

        if scf_drv.rank == 0:
            mo = scf_drv.scf_tensors['C_alpha']
            nocc = molecule.number_of_alpha_electrons()
            nvir = mo.shape[1] - nocc
            mo_occ = mo[:, :nocc]
            mo_vir = mo[:, nocc:]
            matrices = [np.zeros((nocc, nvir)) for i in range(4)]
            matrices[0][nocc - 1, 1] = 1.0
            matrices[1][nocc - 2, 0] = 1.0
            matrices[2][nocc - 2, 1] = 1.0
            matrices[3][nocc - 1, 0] = 1.0
            ao_matrices = [
                np.linalg.multi_dot([mo_occ, m, mo_vir.T]) for m in matrices
            ]
            rhow_density = AODensityMatrix(ao_matrices, denmat.rest)
        else:
            rhow_density = AODensityMatrix()
        rhow_density.broadcast(scf_drv.rank, scf_drv.comm)

        grid_drv = GridDriver()
        grid_drv.set_level(grid_level)
        molgrid = grid_drv.generate(molecule)

        xc_drv = XCNewIntegrator()
        vxc = xc_drv.integrate_vxc_fock(molecule, basis, gs_density, molgrid,
                                        xcfun_label)
        vxc.reduce_sum(scf_drv.rank, scf_drv.nodes, scf_drv.comm)

        fock = AOFockMatrix(rhow_density)
        xc_drv.integrate_fxc_fock(fock, molecule, basis, rhow_density,
                                  gs_density, molgrid, xcfun_label)
        fock.reduce_sum(scf_drv.rank, scf_drv.nodes, scf_drv.comm)

        if scf_drv.rank == mpi_master():
            gto = xc_drv.compute_gto_values(molecule, basis, molgrid)
            Fmat = np.matmul(gs_density.alpha_to_numpy(0), gto)
            rho_a = np.diag(np.matmul(gto.T, Fmat))
            rho_b = rho_a.copy()

            npoints = gto.shape[1]
            rho = np.zeros(npoints * 2)
            rho[0::2] = rho_a[:]
            rho[1::2] = rho_b[:]

            exc_vxc = xc_drv.compute_exc_vxc_for_lda(xcfun_label, rho)
            exc = exc_vxc['exc'][:, 0]
            vrho_a = exc_vxc['vrho'][:, 0]

            gw = molgrid.w_to_numpy()
            rho_total = rho_a + rho_b

            n_ene = np.sum(gw * rho_total)
            assert np.abs(n_ene - vxc.get_electrons()) < tol

            xc_ene = np.sum(gw * exc * rho_total)
            assert np.abs(xc_ene - vxc.get_energy()) < tol

            Gmat = gto * gw * vrho_a

            Vxcmat = np.matmul(gto, Gmat.T)
            assert np.max(np.abs(Vxcmat - vxc.get_matrix().to_numpy())) < tol

            fxc = xc_drv.compute_fxc_for_lda(xcfun_label, rho)
            v2rho2_a_a = fxc['v2rho2'][:, 0]
            v2rho2_a_b = fxc['v2rho2'][:, 1]

            for i in range(rhow_density.number_of_density_matrices()):
                Fmat = np.matmul(rhow_density.alpha_to_numpy(i), gto)
                rhow_a = np.diag(np.matmul(gto.T, Fmat))
                rhow_b = rhow_a.copy()

                Gmat = gto * gw * (v2rho2_a_a * rhow_a + v2rho2_a_b * rhow_b)

                Fxcmat = np.matmul(gto, Gmat.T)
                assert np.max(np.abs(Fxcmat - fock.alpha_to_numpy(i))) < tol

    @pytest.mark.skipif(not is_single_node(), reason='single node only')
    def test_b3lyp(self):

        mol_str = """
            O  0.0000000000   0.0000000000  -0.0254395383
            H  0.0000000000   0.7695699584   0.5948147012
            H  0.0000000000  -0.7695699584   0.5948147012
            O  4.5000000000   0.0000000000  -0.0254395383
            H  4.5000000000   0.7695699584   0.5948147012
            H  4.5000000000  -0.7695699584   0.5948147012
        """
        basis_label = 'def2-svp'
        xcfun_label = 'b3lyp'
        grid_level = 1
        tol = 1.0e-10

        molecule = Molecule.read_str(mol_str, units='angstrom')
        basis = MolecularBasis.read(molecule, basis_label, ostream=None)

        scf_drv = ScfRestrictedDriver()
        scf_drv.xcfun = xcfun_label
        scf_drv.grid_level = grid_level
        scf_drv.ostream.mute()
        scf_drv.compute(molecule, basis)
        gs_density = scf_drv.density

        if scf_drv.rank == 0:
            mo = scf_drv.scf_tensors['C_alpha']
            nocc = molecule.number_of_alpha_electrons()
            nvir = mo.shape[1] - nocc
            mo_occ = mo[:, :nocc]
            mo_vir = mo[:, nocc:]
            matrices = [np.zeros((nocc, nvir)) for i in range(4)]
            matrices[0][nocc - 1, 1] = 1.0
            matrices[1][nocc - 2, 0] = 1.0
            matrices[2][nocc - 2, 1] = 1.0
            matrices[3][nocc - 1, 0] = 1.0
            ao_matrices = [
                np.linalg.multi_dot([mo_occ, m, mo_vir.T]) for m in matrices
            ]
            rhow_density = AODensityMatrix(ao_matrices, denmat.rest)
        else:
            rhow_density = AODensityMatrix()
        rhow_density.broadcast(scf_drv.rank, scf_drv.comm)

        grid_drv = GridDriver()
        grid_drv.set_level(grid_level)
        molgrid = grid_drv.generate(molecule)

        xc_drv = XCNewIntegrator()
        vxc = xc_drv.integrate_vxc_fock(molecule, basis, gs_density, molgrid,
                                        xcfun_label)
        vxc.reduce_sum(scf_drv.rank, scf_drv.nodes, scf_drv.comm)

        fock = AOFockMatrix(rhow_density)
        xc_drv.integrate_fxc_fock(fock, molecule, basis, rhow_density,
                                  gs_density, molgrid, xcfun_label)
        fock.reduce_sum(scf_drv.rank, scf_drv.nodes, scf_drv.comm)

        if scf_drv.rank == mpi_master():
            gto, gto_x, gto_y, gto_z = xc_drv.compute_gto_values_and_derivatives(
                molecule, basis, molgrid)
            npoints = gto.shape[1]

            Dmat = gs_density.alpha_to_numpy(0)
            sym_Dmat = 0.5 * (Dmat + Dmat.T)
            Fmat = np.matmul(sym_Dmat, gto)
            rho_a = np.diag(np.matmul(gto.T, Fmat))
            rho_b = rho_a.copy()

            rho_a_x = 2.0 * np.diag(np.matmul(gto_x.T, Fmat))
            rho_a_y = 2.0 * np.diag(np.matmul(gto_y.T, Fmat))
            rho_a_z = 2.0 * np.diag(np.matmul(gto_z.T, Fmat))
            rho_b_x = rho_a_x.copy()
            rho_b_y = rho_a_y.copy()
            rho_b_z = rho_a_z.copy()

            sigma_a_a = rho_a_x * rho_a_x + rho_a_y * rho_a_y + rho_a_z * rho_a_z
            sigma_a_b = rho_a_x * rho_b_x + rho_a_y * rho_b_y + rho_a_z * rho_b_z
            sigma_b_b = rho_b_x * rho_b_x + rho_b_y * rho_b_y + rho_b_z * rho_b_z

            rho = np.zeros(npoints * 2)
            rho[0::2] = rho_a[:]
            rho[1::2] = rho_b[:]

            sigma = np.zeros(npoints * 3)
            sigma[0::3] = sigma_a_a[:]
            sigma[1::3] = sigma_a_b[:]
            sigma[2::3] = sigma_b_b[:]

            exc_vxc = xc_drv.compute_exc_vxc_for_gga(xcfun_label, rho, sigma)
            exc = exc_vxc['exc'][:, 0]
            vrho_a = exc_vxc['vrho'][:, 0]
            vsigma_aa = exc_vxc['vsigma'][:, 0]
            vsigma_ab = exc_vxc['vsigma'][:, 1]

            gw = molgrid.w_to_numpy()
            rho_total = rho_a + rho_b

            n_ene = np.sum(gw * rho_total)
            assert np.abs(n_ene - vxc.get_electrons()) < tol

            xc_ene = np.sum(gw * exc * rho_total)
            assert np.abs(xc_ene - vxc.get_energy()) < tol

            Gmat = gto * gw * vrho_a

            v_x = 2 * vsigma_aa * rho_a_x + vsigma_ab * rho_b_x
            v_y = 2 * vsigma_aa * rho_a_y + vsigma_ab * rho_b_y
            v_z = 2 * vsigma_aa * rho_a_z + vsigma_ab * rho_b_z

            Gmat_gga = (gto_x * v_x + gto_y * v_y + gto_z * v_z) * gw

            Vxcmat_gga = np.matmul(gto, Gmat_gga.T)
            Vxcmat = np.matmul(gto, Gmat.T) + Vxcmat_gga + Vxcmat_gga.T
            assert np.max(np.abs(Vxcmat - vxc.get_matrix().to_numpy())) < tol

            fxc = xc_drv.compute_fxc_for_gga(xcfun_label, rho, sigma)
            v2rho2_a_a = fxc['v2rho2'][:, 0]
            v2rho2_a_b = fxc['v2rho2'][:, 1]
            v2rhosigma_a_aa = fxc['v2rhosigma'][:, 0]
            v2rhosigma_a_ab = fxc['v2rhosigma'][:, 1]
            v2rhosigma_a_bb = fxc['v2rhosigma'][:, 2]
            v2rhosigma_b_aa = fxc['v2rhosigma'][:, 3]
            v2rhosigma_b_ab = fxc['v2rhosigma'][:, 4]
            v2sigma2_aa_aa = fxc['v2sigma2'][:, 0]
            v2sigma2_aa_ab = fxc['v2sigma2'][:, 1]
            v2sigma2_aa_bb = fxc['v2sigma2'][:, 2]
            v2sigma2_ab_ab = fxc['v2sigma2'][:, 3]
            v2sigma2_ab_bb = fxc['v2sigma2'][:, 4]

            for i in range(rhow_density.number_of_density_matrices()):
                Dmat = rhow_density.alpha_to_numpy(i)
                sym_Dmat = 0.5 * (Dmat + Dmat.T)
                Fmat = np.matmul(sym_Dmat, gto)
                rhow_a = np.diag(np.matmul(gto.T, Fmat))
                rhow_b = rhow_a.copy()

                rhow_a_x = 2.0 * np.diag(np.matmul(gto_x.T, Fmat))
                rhow_a_y = 2.0 * np.diag(np.matmul(gto_y.T, Fmat))
                rhow_a_z = 2.0 * np.diag(np.matmul(gto_z.T, Fmat))
                rhow_b_x = rhow_a_x.copy()
                rhow_b_y = rhow_a_y.copy()
                rhow_b_z = rhow_a_z.copy()

                grhow_grho_aa = 2 * (rhow_a_x * rho_a_x + rhow_a_y * rho_a_y +
                                     rhow_a_z * rho_a_z)

                grhow_grho_ab = (rhow_a_x * rho_b_x + rhow_a_y * rho_b_y +
                                 rhow_a_z * rho_b_z + rhow_b_x * rho_a_x +
                                 rhow_b_y * rho_a_y + rhow_b_z * rho_a_z)

                grhow_grho_bb = 2 * (rhow_b_x * rho_b_x + rhow_b_y * rho_b_y +
                                     rhow_b_z * rho_b_z)

                Gmat = gto * gw * (v2rho2_a_a * rhow_a + v2rho2_a_b * rhow_b +
                                   v2rhosigma_a_aa * grhow_grho_aa +
                                   v2rhosigma_a_ab * grhow_grho_ab +
                                   v2rhosigma_a_bb * grhow_grho_bb)

                f_aa = (v2rhosigma_a_aa * rhow_a + v2rhosigma_b_aa * rhow_b +
                        v2sigma2_aa_aa * grhow_grho_aa +
                        v2sigma2_aa_ab * grhow_grho_ab +
                        v2sigma2_aa_bb * grhow_grho_bb)

                f_ab = (v2rhosigma_a_ab * rhow_a + v2rhosigma_b_ab * rhow_b +
                        v2sigma2_aa_ab * grhow_grho_aa +
                        v2sigma2_ab_ab * grhow_grho_ab +
                        v2sigma2_ab_bb * grhow_grho_bb)

                f_x = 2 * f_aa * rho_a_x + f_ab * rho_b_x
                f_y = 2 * f_aa * rho_a_y + f_ab * rho_b_y
                f_z = 2 * f_aa * rho_a_z + f_ab * rho_b_z

                f_x += 2 * vsigma_aa * rhow_a_x + vsigma_ab * rhow_b_x
                f_y += 2 * vsigma_aa * rhow_a_y + vsigma_ab * rhow_b_y
                f_z += 2 * vsigma_aa * rhow_a_z + vsigma_ab * rhow_b_z

                Gmat_gga = (gto_x * f_x + gto_y * f_y + gto_z * f_z) * gw

                Fxcmat_gga = np.matmul(gto, Gmat_gga.T)
                Fxcmat = np.matmul(gto, Gmat.T) + Fxcmat_gga + Fxcmat_gga.T
                assert np.max(np.abs(Fxcmat - fock.alpha_to_numpy(i))) < tol
