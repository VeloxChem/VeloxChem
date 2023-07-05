from mpi4py import MPI
import numpy as np

from veloxchem.veloxchemlib import GridDriver, XCIntegrator
from veloxchem.veloxchemlib import is_mpi_master, mpi_master
from veloxchem.veloxchemlib import denmat
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestopendriver import ScfRestrictedOpenDriver
from veloxchem.aodensitymatrix import AODensityMatrix


class TestPDFT:

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
        vxc_mat = xc_drv.integrate_vxc_fock(molecule, basis, scfdrv.density,
                                            molgrid, func)
        vxc_mat.reduce_sum(scfdrv.rank, scfdrv.nodes, scfdrv.comm)

        # Compute total and on-top pair densities
        if scfdrv.rank == mpi_master():
            total_density = scf_results['D_alpha'] + scf_results['D_beta']
            den_mat = AODensityMatrix([total_density], denmat.rest)
        else:
            den_mat = AODensityMatrix()
        den_mat.broadcast(scfdrv.rank, scfdrv.comm)

        # Only in the 2 singly occupied orbitals
        Dact = np.identity(2)
        # True density minus "closed shell"
        D2act = -0.5 * np.einsum('mt,np->mntp', Dact, Dact)

        if scfdrv.rank == mpi_master():
            mo_act = scfdrv.mol_orbs.alpha_to_numpy()[:, [7, 8]]
        else:
            mo_act = None
        mo_act = scfdrv.comm.bcast(mo_act, root=mpi_master())

        pdft_vxc, pdft_wxc = xc_drv.integrate_vxc_pdft(den_mat, D2act,
                                                       mo_act.T.copy(),
                                                       molecule, basis, molgrid,
                                                       pfunc)
        pdft_vxc.reduce_sum(scfdrv.rank, scfdrv.nodes, scfdrv.comm)
        pdft_wxc = scfdrv.comm.reduce(pdft_wxc, op=MPI.SUM, root=mpi_master())

        if scfdrv.rank == mpi_master():
            C = scf_results['C_alpha']
            nAO = C.shape[0]

            # Compute gradients
            nIn = molecule.number_of_beta_electrons()
            nInAct = molecule.number_of_alpha_electrons()
            nAct = nInAct - nIn
            wxc = 2.0 * pdft_wxc.reshape(nAO, nAct, nAct, nAct)

            FA = vxc_mat.get_alpha_matrix().to_numpy()
            FB = vxc_mat.get_beta_matrix().to_numpy()

            FAMO = np.linalg.multi_dot([C.T, FA, C])
            FBMO = np.linalg.multi_dot([C.T, FB, C])
            grad_ks = FAMO[:nInAct, :]
            grad_ks[:nIn, :] += FBMO[:nIn, :]

            Fock = pdft_vxc.get_matrix().to_numpy()
            Qpt = np.tensordot(wxc, D2act, axes=([1, 2, 3], [1, 2, 3]))

            grad_pdft = np.zeros_like(grad_ks)
            grad_pdft[:nIn, :] = np.linalg.multi_dot(
                [C[:, :nIn].T, 2.0 * Fock, C])
            grad_pdft[nIn:nInAct, :] = np.linalg.multi_dot(
                [Dact, mo_act.T, Fock, C])
            grad_pdft[nIn:nInAct, :] += np.matmul(C.T, Qpt).T

            return (vxc_mat.get_energy(), pdft_vxc.get_energy(), grad_ks,
                    grad_pdft)
        else:
            return None, None, None, None

    def test_O2_ROSlater(self):
        ksdft, pdft, ks_grad, pdft_grad = self.run_RODFT('slater', 'tslater')
        if is_mpi_master():
            assert abs(ksdft - pdft) < 1.0e-6
            assert np.allclose(ks_grad, pdft_grad)

    def test_O2_ROLDA(self):
        ksdft, pdft, ks_grad, pdft_grad = self.run_RODFT('slda', 'tlda')
        if is_mpi_master():
            assert abs(ksdft - pdft) < 1.0e-6
            assert np.allclose(ks_grad, pdft_grad)

    def test_O2_ROGGA(self):
        ksdft, pdft, ks_grad, pdft_grad = self.run_RODFT('pbe', 'tpbe')
        if is_mpi_master():
            assert abs(-16.911864099412625 - pdft) < 1.0e-6
