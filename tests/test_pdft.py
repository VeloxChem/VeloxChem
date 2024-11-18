from mpi4py import MPI
import numpy as np
import pytest

from veloxchem.veloxchemlib import XCIntegrator, XCPairDensityFunctional, XCMolecularGradient
from veloxchem.veloxchemlib import available_pdft_functionals
from veloxchem.veloxchemlib import mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.griddriver import GridDriver
from veloxchem.scfrestopendriver import ScfRestrictedOpenDriver
from veloxchem.scfrestdriver import ScfRestrictedDriver


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

    def test_H(self):
        """
        Test all functional components for numerical stability for nbeta = 0
        """

        xyz = "H 0.0 0.0 0.0"

        molecule = Molecule.read_str(xyz)
        basis = MolecularBasis.read(molecule, "aug-cc-pVDZ")
        molecule.set_multiplicity(2)

        # Optimize ROHF wavefunction
        scfdrv = ScfRestrictedOpenDriver()
        scfdrv.ostream.mute()
        scf_results = scfdrv.compute(molecule, basis)

        if scfdrv.rank == mpi_master():
            total_density = scf_results['D_alpha']
        else:
            total_density = None
        total_density = scfdrv.comm.bcast(total_density, root=mpi_master())

        # True density minus "closed shell"
        D2act = np.array([[[[-0.5]]]])

        if scfdrv.rank == mpi_master():
            mo_act = scfdrv.mol_orbs.alpha_to_numpy()[:, [1]]
        else:
            mo_act = None
        mo_act = scfdrv.comm.bcast(mo_act, root=mpi_master())

        # Compute PDFT
        grid_drv = GridDriver()
        molgrid = grid_drv.generate(molecule)

        xc_drv = XCIntegrator()

        # References
        references = {
            "tSlater": -0.21286739663259885,
            "tSlater_erf": -0.0780214343193473,
            "tVWN_RPA": -0.05792602070111518,
            "tVWN5": -0.040929846047812316,
            "tPMGB06": -0.023054878539261928,
            "tP86": -0.018538415905553742,
            "tPBE_X": -0.2543064204128253,
            "tPBEX_erf": -0.09113476160995143,
            "tPBE_C": -0.014861790644530826,
            "tPBEC_erf": -0.01496306351838878,
            "tB88": -0.045962079987423625,
            "tB88_erf": -0.011016744566851068,
            "tLYP": -0.013596771624186594,
            "tLYP_erf": -0.01415140810845559,
            "HPG20": -0.0139680482003766
        }

        for func in available_pdft_functionals():
            pdft_vxc, pdft_wxc = xc_drv.integrate_vxc_pdft(
                total_density, D2act, mo_act.T.copy(), molecule, basis, molgrid,
                func, {func: 1.0}, 0.4)

            pdft_xc_energy = scfdrv.comm.reduce(pdft_vxc.get_energy(),
                                                root=mpi_master())

            pdft_vxc_mat = scfdrv.comm.reduce(pdft_vxc.alpha_to_numpy(),
                                              root=mpi_master())

            pdft_wxc = scfdrv.comm.reduce(pdft_wxc,
                                          op=MPI.SUM,
                                          root=mpi_master())

            if scfdrv.comm.rank == mpi_master():
                if func in references:
                    assert (pdft_xc_energy - references[func]) < 1.0e-8
                else:
                    assert not np.isnan(pdft_xc_energy)
                assert not np.isnan(pdft_vxc_mat).any()
                assert not np.isnan(pdft_vxc_mat).any()

    def test_O2_ROLDA(self):

        O2_xyz = """
            O 0.0 0.0 -0.6
            O 0.0 0.0  0.6
        """
        molecule = Molecule.read_str(O2_xyz)
        molecule.set_multiplicity(3)
        basis = MolecularBasis.read(molecule, 'cc-pvdz', ostream=None)

        func = "SLDA"
        pfunc = {
            "name": "tLDA",
            "components": {
                "TSLATER": 1.0,
                'TVWN_RPA': 1.0
            }
        }

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
                                                       pfunc["name"],
                                                       pfunc["components"], 0.0)

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

            assert abs(xc_energy - pdft_xc_energy) < 1.0e-6
            assert np.allclose(grad_ks, grad_pdft)

    def run_RDFT(self, func, pfunc):

        mol_str = """
        H    0.7493682    0.0000000    0.4424329
        O    0.0000000    0.0000000   -0.1653507
        H   -0.7493682    0.0000000    0.4424329
        """
        molecule = Molecule.read_str(mol_str, units='angstrom')
        basis = MolecularBasis.read(molecule, "def2-svp")

        # Optimize ROHF wavefunction
        scfdrv = ScfRestrictedDriver()
        scfdrv.ostream.mute()
        scf_results = scfdrv.compute(molecule, basis)

        # Compute DFT correction
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

        # Compute total and on-top pair densities
        if scfdrv.rank == mpi_master():
            total_density = scf_results['D_alpha'] + scf_results['D_beta']
        else:
            total_density = None
        total_density = scfdrv.comm.bcast(total_density, root=mpi_master())

        # Arrays of size 0
        D2act = np.empty((0, 0, 0, 0))
        mo_act = np.empty((Da.shape[0], 0))

        pdft_vxc, pdft_wxc = xc_drv.integrate_vxc_pdft(total_density, D2act,
                                                       mo_act.T.copy(),
                                                       molecule, basis, molgrid,
                                                       pfunc["name"],
                                                       pfunc["components"], 0.0)

        pdft_xc_energy = scfdrv.comm.reduce(pdft_vxc.get_energy(),
                                            root=mpi_master())

        pdft_np_xcmat_a = scfdrv.comm.reduce(pdft_vxc.alpha_to_numpy(),
                                             root=mpi_master())

        return xc_energy, pdft_xc_energy, np_xcmat_a, pdft_np_xcmat_a

    def test_H20_PBE(self):
        pfunc = {"name": "tPBE", "components": {"TPBE_X": 1.0, 'TPBE_C': 1.0}}
        xc_energy, pdft_xc_energy, np_xcmat_a, pdft_np_xcmat_a = self.run_RDFT(
            'pbe', pfunc)
        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            assert abs(xc_energy - pdft_xc_energy) < 3.0e-6
            assert np.allclose(np_xcmat_a, pdft_np_xcmat_a)

    def test_H2O_BLYP(self):
        pfunc = {
            "name": "tBLYP",
            "components": {
                "TSLATER": 1.0,
                'TB88': 1.0,
                'TLYP': 1.0
            }
        }
        xc_energy, pdft_xc_energy, np_xcmat_a, pdft_np_xcmat_a = self.run_RDFT(
            'blyp', pfunc)
        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            assert abs(xc_energy - pdft_xc_energy) < 1.0e-6
            assert np.allclose(np_xcmat_a, pdft_np_xcmat_a)

    def test_H2O_BP86(self):
        pfunc = {
            "name": "tBLYP",
            "components": {
                "TSLATER": 1.0,
                'TB88': 1.0,
                'TP86': 1.0
            }
        }
        xc_energy, pdft_xc_energy, np_xcmat_a, pdft_np_xcmat_a = self.run_RDFT(
            'bp86', pfunc)
        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            assert abs(xc_energy - pdft_xc_energy) < 1.0e-6
            assert np.allclose(np_xcmat_a, pdft_np_xcmat_a)

    def test_grad_tlda(self):
        H2O_xyz = """3

O   0.0   0.0   0.0
H   0.0   1.4   1.1
H   0.0  -1.4   1.1
        """
        molecule = Molecule.from_xyz_string(H2O_xyz)
        basis = MolecularBasis.read(molecule, "cc-pvdz")

        grid_drv = GridDriver()
        molgrid = grid_drv.generate(molecule)

        func = "tLDA"
        func_comp = {"TSLATER": 1.0, 'TVWN_RPA': 1.0}

        # densities and orbitals from a converged cas(4,4) MC-tBLYP calculation
        dmat_ao = np.load("data/h2o.ccpvdz.cas44.dao.npy")
        ns_d2 = np.load("data/h2o.ccpvdz.cas44.d2.npy")
        act_mos = np.load("data/h2o.ccpvdz.cas44.actorbs.npy")
        xc_drv = XCMolecularGradient()
        dft_gradients = xc_drv.integrate_vxc_pdft_gradient(molecule, basis, dmat_ao, ns_d2,
                                             act_mos, molgrid, func, func_comp, 0.0)

        ref_grad = np.array([[0, 0, -0.128566279], [0, 0.0894509471, 0.0642872016], [0, -0.0894509471, 0.0642872016]])
        assert np.allclose(dft_gradients, ref_grad)

    def test_grad_tblyp(self):
        H2O_xyz = """3

O   0.0   0.0   0.0
H   0.0   1.4   1.1
H   0.0  -1.4   1.1
        """
        molecule = Molecule.from_xyz_string(H2O_xyz)
        basis = MolecularBasis.read(molecule, "cc-pvdz")

        grid_drv = GridDriver()
        molgrid = grid_drv.generate(molecule)

        func = "tBLYP"
        func_comp = {"TSLATER": 1.0, 'TB88': 1.0, 'TLYP': 1.0}

        # densities and orbitals from a converged cas(4,4) MC-tBLYP calculation
        dmat_ao = np.load("data/h2o.ccpvdz.cas44.dao.npy")
        ns_d2 = np.load("data/h2o.ccpvdz.cas44.d2.npy")
        act_mos = np.load("data/h2o.ccpvdz.cas44.actorbs.npy")
        xc_drv = XCMolecularGradient()
        dft_gradients = xc_drv.integrate_vxc_pdft_gradient(molecule, basis, dmat_ao, ns_d2,
                                             act_mos, molgrid, func, func_comp, 0.0)

        ref_grad = np.array([[0, 0, -0.116818315], [0, 0.0808942190, 0.0584138119], [0, -0.0808942190, 0.0584138119]])
        assert np.allclose(dft_gradients, ref_grad)
