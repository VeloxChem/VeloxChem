from mpi4py import MPI
import numpy as np
import pytest

from veloxchem.veloxchemlib import XCIntegrator, XCPairDensityFunctional
from veloxchem.veloxchemlib import mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.griddriver import GridDriver
from veloxchem.scfrestopendriver import ScfRestrictedOpenDriver


@pytest.mark.solvers
class TestPDFT:

    def test_pfunc(self):

        pfunc = XCPairDensityFunctional('TSLATER', ['TSLATER'], [1.0])
        assert pfunc.get_func_label() == 'TSLATER'
        assert pfunc.get_func_type() == 'PLDA'

        pfunc = XCPairDensityFunctional('TLDA', ['TSLATER', 'TVWN'], [1.0, 1.0])
        assert pfunc.get_func_label() == 'TLDA'
        assert pfunc.get_func_type() == 'PLDA'

        pfunc = XCPairDensityFunctional('TPBE', ['TPBE_X', 'TPBE_C'],
                                        [1.0, 1.0])
        assert pfunc.get_func_label() == 'TPBE'
        assert pfunc.get_func_type() == 'PGGA'

        pfunc = XCPairDensityFunctional('TBLYP', ['TSLATER', 'TB88', 'TLYP'],
                                        [1.0, 1.0, 1.0])
        assert pfunc.get_func_label() == 'TBLYP'
        assert pfunc.get_func_type() == 'PGGA'

        pfunc2 = XCPairDensityFunctional('TSLATER', ['TSLATER'], [1.0])
        assert pfunc != pfunc2

        pfunc2 = XCPairDensityFunctional('TBLYP', ['TSLATER', 'TB88', 'TLYP'],
                                         [1.0, 1.0, 1.0])
        assert pfunc == pfunc2

    def run_RODFT(self, func, pfunc):

        O2_xyz = """
            O 0.0 0.0 -0.6
            O 0.0 0.0  0.6
        """
        molecule = Molecule.read_str(O2_xyz)
        molecule.set_multiplicity(3)
        basis = MolecularBasis.read(molecule, 'cc-pvdz', ostream=None)

        # Optimize ROHF wavefunction
        scfdrv = ScfRestrictedOpenDriver()
        scfdrv.xcfun = func
        scfdrv.ostream.mute()
        scf_results = scfdrv.compute(molecule, basis)

        # Compute SLDA correction
        grid_drv = GridDriver()
        molgrid = grid_drv.generate(molecule)

        xc_drv = XCIntegrator()
        if scfdrv.rank == mpi_master():
            Da, Db = scf_results['D_alpha'], scf_results['D_beta']
        else:
            Da, Db = None, None
        Da = scfdrv.comm.bcast(Da, root=mpi_master())
        Db = scfdrv.comm.bcast(Db, root=mpi_master())

        vxc_mat = xc_drv.integrate_vxc_fock(molecule, basis, [Da, Db], molgrid,
                                            func)

        xc_energy = scfdrv.comm.reduce(vxc_mat.get_energy(), root=mpi_master())

        np_xcmat_a = scfdrv.comm.reduce(vxc_mat.alpha_to_numpy(),
                                        root=mpi_master())
        np_xcmat_b = scfdrv.comm.reduce(vxc_mat.beta_to_numpy(),
                                        root=mpi_master())

        # Compute total and on-top pair densities
        if scfdrv.rank == mpi_master():
            total_density = scf_results['D_alpha'] + scf_results['D_beta']
        else:
            total_density = None
        total_density = scfdrv.comm.bcast(total_density, root=mpi_master())

        # Only in the 2 singly occupied orbitals
        Dact = np.identity(2)
        # True density minus "closed shell"
        D2act = -0.5 * np.einsum('mt,np->mntp', Dact, Dact)

        if scfdrv.rank == mpi_master():
            mo_act = scfdrv.mol_orbs.alpha_to_numpy()[:, [7, 8]]
        else:
            mo_act = None
        mo_act = scfdrv.comm.bcast(mo_act, root=mpi_master())

        pdft_vxc, pdft_wxc = xc_drv.integrate_vxc_pdft(total_density, D2act,
                                                       mo_act.T.copy(),
                                                       molecule, basis, molgrid,
                                                       pfunc)

        pdft_xc_energy = scfdrv.comm.reduce(pdft_vxc.get_energy(),
                                            root=mpi_master())

        pdft_np_xcmat_a = scfdrv.comm.reduce(pdft_vxc.alpha_to_numpy(),
                                             root=mpi_master())

        pdft_wxc = scfdrv.comm.reduce(pdft_wxc, op=MPI.SUM, root=mpi_master())

        if scfdrv.rank == mpi_master():
            C = scf_results['C_alpha']
            nAO = C.shape[0]

            # Compute gradients
            nIn = molecule.number_of_beta_electrons()
            nInAct = molecule.number_of_alpha_electrons()
            nAct = nInAct - nIn
            wxc = 2.0 * pdft_wxc.reshape(nAO, nAct, nAct, nAct)

            FA = np_xcmat_a
            FB = np_xcmat_b

            FAMO = np.linalg.multi_dot([C.T, FA, C])
            FBMO = np.linalg.multi_dot([C.T, FB, C])
            grad_ks = FAMO[:nInAct, :]
            grad_ks[:nIn, :] += FBMO[:nIn, :]

            Fock = pdft_np_xcmat_a
            Qpt = np.tensordot(wxc, D2act, axes=([1, 2, 3], [1, 2, 3]))

            grad_pdft = np.zeros_like(grad_ks)
            grad_pdft[:nIn, :] = np.linalg.multi_dot(
                [C[:, :nIn].T, 2.0 * Fock, C])
            grad_pdft[nIn:nInAct, :] = np.linalg.multi_dot(
                [Dact, mo_act.T, Fock, C])
            grad_pdft[nIn:nInAct, :] += np.matmul(C.T, Qpt).T

            return (xc_energy, pdft_xc_energy, grad_ks, grad_pdft)
        else:
            return None, None, None, None

    def test_O2_ROSlater(self):
        ksdft, pdft, ks_grad, pdft_grad = self.run_RODFT('slater', 'tslater')
        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            assert abs(ksdft - pdft) < 1.0e-6
            assert np.allclose(ks_grad, pdft_grad)

    def test_O2_ROLDA(self):
        ksdft, pdft, ks_grad, pdft_grad = self.run_RODFT('slda', 'tlda')
        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            assert abs(ksdft - pdft) < 1.0e-6
            assert np.allclose(ks_grad, pdft_grad)

    def test_O2_ROPBE(self):
        ksdft, pdft, ks_grad, pdft_grad = self.run_RODFT('pbe', 'tpbe')
        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            assert abs(-16.911864099412625 - pdft) < 1.0e-6

    def test_O2_ROBLYP(self):
        ksdft, pdft, ks_grad, pdft_grad = self.run_RODFT('blyp', 'tblyp')
        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            assert abs(-17.056873017749865 - pdft) < 1.0e-6
