import numpy as np
import pytest

from veloxchem.veloxchemlib import GridDriver, XCNewIntegrator, XCNewFunctional
from veloxchem.veloxchemlib import is_single_node, mpi_master, new_parse_xc_func
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver


class TestFunctionalExcVxc:

    @pytest.mark.skipif(not is_single_node(), reason="single node only")
    def test_lda(self):

        mol_str = """
            O  0.0000000000   0.0000000000  -0.0254395383
            H  0.0000000000   0.7695699584   0.5948147012
            H  0.0000000000  -0.7695699584   0.5948147012
            O  4.5000000000   0.0000000000  -0.0254395383
            H  4.5000000000   0.7695699584   0.5948147012
            H  4.5000000000  -0.7695699584   0.5948147012
        """
        basis_label = "def2-svp"
        xcfun_label = "slda"
        grid_level = 1
        tol = 1.0e-10

        molecule = Molecule.read_str(mol_str, units="angstrom")
        basis = MolecularBasis.read(molecule, basis_label, ostream=None)

        scf_drv = ScfRestrictedDriver()
        scf_drv.xcfun = xcfun_label
        scf_drv.grid_level = grid_level
        scf_drv.ostream.state = False
        scf_drv.compute(molecule, basis)
        gs_density = scf_drv.density

        grid_drv = GridDriver()
        grid_drv.set_level(grid_level)
        molgrid = grid_drv.generate(molecule)

        xc_drv = XCNewIntegrator()
        vxc = xc_drv.integrate_vxc_fock(molecule, basis, gs_density, molgrid,
                                        xcfun_label)
        vxc.reduce_sum(scf_drv.rank, scf_drv.nodes, scf_drv.comm)

        func = new_parse_xc_func('slda')

        func_ref = XCNewFunctional(
            'slda',
            ['LDA_X', 'LDA_C_VWN_RPA'],
            [1.0, 1.0],
        )
        assert func == func_ref

        if scf_drv.rank == mpi_master():
            gto = xc_drv.compute_gto_values(molecule, basis, molgrid)
            Fmat = np.matmul(gs_density.alpha_to_numpy(0), gto)
            rho_a = np.diag(np.matmul(gto.T, Fmat))
            rho_b = rho_a.copy()

            npoints = gto.shape[1]
            rho = np.zeros(npoints * 2)
            rho[0::2] = rho_a[:]
            rho[1::2] = rho_b[:]

            exc, vrho = func.compute_exc_vxc_for_lda(rho)
            vrho_a = vrho[0::2]

            gw = molgrid.w_to_numpy()
            rho_total = rho_a + rho_b

            n_ene = np.sum(gw * rho_total)
            assert np.abs(n_ene - vxc.get_electrons()) < tol

            xc_ene = np.sum(gw * exc * rho_total)
            assert np.abs(xc_ene - vxc.get_energy()) < tol

            Gmat = gto * gw * vrho_a

            Vxcmat = np.matmul(gto, Gmat.T)
            assert np.max(np.abs(Vxcmat - vxc.get_matrix().to_numpy())) < tol

    @pytest.mark.skipif(not is_single_node(), reason="single node only")
    def test_gga(self):

        mol_str = """
            O  0.0000000000   0.0000000000  -0.0254395383
            H  0.0000000000   0.7695699584   0.5948147012
            H  0.0000000000  -0.7695699584   0.5948147012
            O  4.5000000000   0.0000000000  -0.0254395383
            H  4.5000000000   0.7695699584   0.5948147012
            H  4.5000000000  -0.7695699584   0.5948147012
        """
        basis_label = "def2-svp"
        xcfun_label = "b3lyp"
        grid_level = 1
        tol = 1.0e-10

        molecule = Molecule.read_str(mol_str, units="angstrom")
        basis = MolecularBasis.read(molecule, basis_label, ostream=None)

        scf_drv = ScfRestrictedDriver()
        scf_drv.xcfun = xcfun_label
        scf_drv.grid_level = grid_level
        scf_drv.ostream.state = False
        scf_drv.compute(molecule, basis)
        gs_density = scf_drv.density

        grid_drv = GridDriver()
        grid_drv.set_level(grid_level)
        molgrid = grid_drv.generate(molecule)

        xc_drv = XCNewIntegrator()
        vxc = xc_drv.integrate_vxc_fock(molecule, basis, gs_density, molgrid,
                                        xcfun_label)
        vxc.reduce_sum(scf_drv.rank, scf_drv.nodes, scf_drv.comm)

        func = new_parse_xc_func('b3lyp')

        func_ref = XCNewFunctional(
            'b3lyp',
            ['LDA_X', 'GGA_X_B88', 'LDA_C_VWN_RPA', 'GGA_C_LYP'],
            [0.08, 0.72, 0.19, 0.81],
            0.2,
        )
        assert func == func_ref

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

            exc, vrho, vsigma = func.compute_exc_vxc_for_gga(rho, sigma)
            vrho_a = vrho[0::2]
            vsigma_aa = vsigma[0::3]
            vsigma_ab = vsigma[1::3]

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
