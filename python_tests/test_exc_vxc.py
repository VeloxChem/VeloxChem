import numpy as np
import pytest

from veloxchem.veloxchemlib import GridDriver, XCNewIntegrator
from veloxchem.veloxchemlib import is_single_node, mpi_master
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
        basis = MolecularBasis.read(molecule, basis_label)

        scf_drv = ScfRestrictedDriver()
        scf_drv.xcfun = xcfun_label
        scf_drv.grid_level = grid_level
        scf_drv.compute(molecule, basis)
        density = scf_drv.density

        grid_drv = GridDriver()
        grid_drv.set_level(grid_level)
        molgrid = grid_drv.generate(molecule)
        molgrid.partition_grid_points()
        molgrid.distribute_counts_and_displacements(scf_drv.rank, scf_drv.nodes,
                                                    scf_drv.comm)

        xc_drv = XCNewIntegrator()
        vxc = xc_drv.integrate_vxc_fock(molecule, basis, density, molgrid,
                                        xcfun_label)
        vxc.reduce_sum(scf_drv.rank, scf_drv.nodes, scf_drv.comm)

        if scf_drv.rank == mpi_master():
            gto = xc_drv.compute_gto_values(molecule, basis, molgrid)
            Fmat = np.matmul(density.alpha_to_numpy(0), gto)
            rhoa = np.diag(np.matmul(gto.T, Fmat))
            rhob = rhoa.copy()

            npoints = gto.shape[1]
            rho = np.zeros(npoints * 2)
            rho[0::2] = rhoa[:]
            rho[1::2] = rhob[:]

            exc = np.zeros(npoints)
            vrho = np.zeros(npoints * 2)
            exc, vrho = xc_drv.compute_exc_vxc_for_lda(xcfun_label, rho)

            gw = molgrid.w_to_numpy()
            n_ene = np.sum(gw * (rhoa + rhob))
            assert np.abs(n_ene - vxc.get_electrons()) < tol

            xc_ene = np.sum(gw * exc * (rhoa + rhob))
            assert np.abs(xc_ene - vxc.get_energy()) < tol

            vrhoa = vrho[0::2]
            Gmat = gto * gw * vrhoa
            Vxcmat = np.matmul(gto, Gmat.T)
            assert np.max(np.abs(Vxcmat - vxc.get_matrix().to_numpy())) < tol

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
        basis = MolecularBasis.read(molecule, basis_label)

        scf_drv = ScfRestrictedDriver()
        scf_drv.xcfun = xcfun_label
        scf_drv.grid_level = grid_level
        scf_drv.compute(molecule, basis)
        density = scf_drv.density

        grid_drv = GridDriver()
        grid_drv.set_level(grid_level)
        molgrid = grid_drv.generate(molecule)
        molgrid.partition_grid_points()
        molgrid.distribute_counts_and_displacements(scf_drv.rank, scf_drv.nodes,
                                                    scf_drv.comm)

        xc_drv = XCNewIntegrator()
        vxc = xc_drv.integrate_vxc_fock(molecule, basis, density, molgrid,
                                        xcfun_label)
        vxc.reduce_sum(scf_drv.rank, scf_drv.nodes, scf_drv.comm)

        if scf_drv.rank == mpi_master():
            gto, gtox, gtoy, gtoz = xc_drv.compute_gto_values_and_derivatives(
                molecule, basis, molgrid)
            npoints = gto.shape[1]

            Dmat = density.alpha_to_numpy(0)
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

            exc = np.zeros(npoints)
            vrho = np.zeros(npoints * 2)
            vsigma = np.zeros(npoints * 3)
            exc, vrho, vsigma = xc_drv.compute_exc_vxc_for_gga(
                xcfun_label, rho, sigma)

            gw = molgrid.w_to_numpy()
            n_ene = np.sum(gw * (rhoa + rhob))
            assert np.abs(n_ene - vxc.get_electrons()) < tol

            xc_ene = np.sum(gw * exc * (rhoa + rhob))
            assert np.abs(xc_ene - vxc.get_energy()) < tol

            vrhoa = vrho[0::2]
            vsigma_aa = vsigma[0::3]
            vsigma_ab = vsigma[1::3]
            Gmat = gto * gw * vrhoa
            Gmat_gga = (gtox * rhoax + gtoy * rhoay +
                        gtoz * rhoaz) * gw * (2.0 * vsigma_aa + vsigma_ab)
            Vxcmat_gga = np.matmul(gto, Gmat_gga.T)
            Vxcmat = np.matmul(gto, Gmat.T) + Vxcmat_gga + Vxcmat_gga.T
            assert np.max(np.abs(Vxcmat - vxc.get_matrix().to_numpy())) < tol
