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
        scf_drv.ostream.state = False
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
        molgrid.partition_grid_points()
        molgrid.distribute_counts_and_displacements(scf_drv.rank, scf_drv.nodes,
                                                    scf_drv.comm)

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
            rhoa = np.diag(np.matmul(gto.T, Fmat))
            rhob = rhoa.copy()

            npoints = gto.shape[1]
            rho = np.zeros(npoints * 2)
            rho[0::2] = rhoa[:]
            rho[1::2] = rhob[:]

            exc, vrho = xc_drv.compute_exc_vxc_for_lda(xcfun_label, rho)
            vrhoa = vrho[0::2]

            gw = molgrid.w_to_numpy()
            n_ene = np.sum(gw * (rhoa + rhob))
            assert np.abs(n_ene - vxc.get_electrons()) < tol

            xc_ene = np.sum(gw * exc * (rhoa + rhob))
            assert np.abs(xc_ene - vxc.get_energy()) < tol

            Gmat = gto * gw * vrhoa
            Vxcmat = np.matmul(gto, Gmat.T)
            assert np.max(np.abs(Vxcmat - vxc.get_matrix().to_numpy())) < tol

            v2rho2 = xc_drv.compute_fxc_for_lda(xcfun_label, rho)
            v2rho2_aa = v2rho2[0::3]
            v2rho2_ab = v2rho2[1::3]

            for i in range(rhow_density.number_of_density_matrices()):
                Fmat = np.matmul(rhow_density.alpha_to_numpy(i), gto)
                rhowa = np.diag(np.matmul(gto.T, Fmat))
                rhowb = rhowa.copy()

                Gmat = gto * gw * (v2rho2_aa * rhowa + v2rho2_ab * rhowb)
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
        scf_drv.ostream.state = False
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
        molgrid.partition_grid_points()
        molgrid.distribute_counts_and_displacements(scf_drv.rank, scf_drv.nodes,
                                                    scf_drv.comm)

        xc_drv = XCNewIntegrator()
        vxc = xc_drv.integrate_vxc_fock(molecule, basis, gs_density, molgrid,
                                        xcfun_label)
        vxc.reduce_sum(scf_drv.rank, scf_drv.nodes, scf_drv.comm)

        fock = AOFockMatrix(rhow_density)
        xc_drv.integrate_fxc_fock(fock, molecule, basis, rhow_density,
                                  gs_density, molgrid, xcfun_label)
        fock.reduce_sum(scf_drv.rank, scf_drv.nodes, scf_drv.comm)

        if scf_drv.rank == mpi_master():
            gto, gtox, gtoy, gtoz = xc_drv.compute_gto_values_and_derivatives(
                molecule, basis, molgrid)
            npoints = gto.shape[1]

            Dmat = gs_density.alpha_to_numpy(0)
            sym_Dmat = 0.5 * (Dmat + Dmat.T)
            Fmat = np.matmul(sym_Dmat, gto)
            rhoa = np.diag(np.matmul(gto.T, Fmat))
            rhob = rhoa.copy()

            rhoax = 2.0 * np.diag(np.matmul(gtox.T, Fmat))
            rhoay = 2.0 * np.diag(np.matmul(gtoy.T, Fmat))
            rhoaz = 2.0 * np.diag(np.matmul(gtoz.T, Fmat))
            rhobx = rhoax.copy()
            rhoby = rhoay.copy()
            rhobz = rhoaz.copy()

            sigma_aa = rhoax * rhoax + rhoay * rhoay + rhoaz * rhoaz
            sigma_ab = rhoax * rhobx + rhoay * rhoby + rhoaz * rhobz
            sigma_bb = rhobx * rhobx + rhoby * rhoby + rhobz * rhobz

            rho = np.zeros(npoints * 2)
            rho[0::2] = rhoa[:]
            rho[1::2] = rhob[:]

            sigma = np.zeros(npoints * 3)
            sigma[0::3] = sigma_aa[:]
            sigma[1::3] = sigma_ab[:]
            sigma[2::3] = sigma_bb[:]

            exc, vrho, vsigma = xc_drv.compute_exc_vxc_for_gga(
                xcfun_label, rho, sigma)
            vrhoa = vrho[0::2]
            vsigma_aa = vsigma[0::3]
            vsigma_ab = vsigma[1::3]

            gw = molgrid.w_to_numpy()
            n_ene = np.sum(gw * (rhoa + rhob))
            assert np.abs(n_ene - vxc.get_electrons()) < tol

            xc_ene = np.sum(gw * exc * (rhoa + rhob))
            assert np.abs(xc_ene - vxc.get_energy()) < tol

            Gmat = gto * gw * vrhoa
            Gmat_gga = (gtox * rhoax + gtoy * rhoay +
                        gtoz * rhoaz) * gw * (2.0 * vsigma_aa + vsigma_ab)
            Vxcmat_gga = np.matmul(gto, Gmat_gga.T)
            Vxcmat = np.matmul(gto, Gmat.T) + Vxcmat_gga + Vxcmat_gga.T
            assert np.max(np.abs(Vxcmat - vxc.get_matrix().to_numpy())) < tol

            v2rho2, v2rhosigma, v2sigma2 = xc_drv.compute_fxc_for_gga(
                xcfun_label, rho, sigma)
            v2rho2_aa = v2rho2[0::3]
            v2rho2_ab = v2rho2[1::3]
            v2rhosigma_aa = v2rhosigma[0::6]
            v2rhosigma_ac = v2rhosigma[1::6]
            v2rhosigma_ab = v2rhosigma[2::6]
            v2rhosigma_bc = v2rhosigma[4::6]
            v2sigma2_aa = v2sigma2[0::6]
            v2sigma2_ac = v2sigma2[1::6]
            v2sigma2_ab = v2sigma2[2::6]
            v2sigma2_cc = v2sigma2[3::6]
            v2sigma2_cb = v2sigma2[4::6]

            for i in range(rhow_density.number_of_density_matrices()):
                Dmat = rhow_density.alpha_to_numpy(i)
                sym_Dmat = 0.5 * (Dmat + Dmat.T)
                Fmat = np.matmul(sym_Dmat, gto)
                rhowa = np.diag(np.matmul(gto.T, Fmat))
                rhowb = rhowa.copy()

                rhowax = 2.0 * np.diag(np.matmul(gtox.T, Fmat))
                rhoway = 2.0 * np.diag(np.matmul(gtoy.T, Fmat))
                rhowaz = 2.0 * np.diag(np.matmul(gtoz.T, Fmat))
                rhowbx = rhowax.copy()
                rhowby = rhoway.copy()
                rhowbz = rhowaz.copy()

                grhow_grho_aa = rhowax * rhoax + rhoway * rhoay + rhowaz * rhoaz
                grhow_grho_ab = rhowax * rhobx + rhoway * rhoby + rhowaz * rhobz
                grhow_grho_ba = rhowbx * rhoax + rhowby * rhoay + rhowbz * rhoaz
                grhow_grho_bb = rhowbx * rhobx + rhowby * rhoby + rhowbz * rhobz

                grhow_grho_ab_ba = grhow_grho_ab + grhow_grho_ba

                Gmat = gto * gw * (v2rho2_aa * rhowa + v2rho2_ab * rhowb +
                                   2 * v2rhosigma_aa * grhow_grho_aa +
                                   v2rhosigma_ac * grhow_grho_ab_ba +
                                   2 * v2rhosigma_ab * grhow_grho_bb)
                Fxcmat = np.matmul(gto, Gmat.T)

                fac_1 = ((2 * v2rhosigma_aa + v2rhosigma_ac) * rhowa +
                         (2 * v2rhosigma_ab + v2rhosigma_bc) * rhowb +
                         (4 * v2sigma2_aa + 2 * v2sigma2_ac) * grhow_grho_aa +
                         (4 * v2sigma2_ab + 2 * v2sigma2_cb) * grhow_grho_bb +
                         (2 * v2sigma2_ac + v2sigma2_cc) * grhow_grho_ab_ba)

                fac_2 = 2 * vsigma_aa + vsigma_ab

                xcomp = fac_1 * rhoax + fac_2 * rhowax
                ycomp = fac_1 * rhoay + fac_2 * rhoway
                zcomp = fac_1 * rhoaz + fac_2 * rhowaz

                Gmat_gga = (gtox * xcomp + gtoy * ycomp + gtoz * zcomp) * gw

                Fxcmat_gga = np.matmul(gto, Gmat_gga.T)
                Fxcmat = np.matmul(gto, Gmat.T) + Fxcmat_gga + Fxcmat_gga.T
                assert np.max(np.abs(Fxcmat - fock.alpha_to_numpy(i))) < tol
