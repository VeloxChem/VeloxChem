import numpy as np

from .veloxchemlib import ElectronRepulsionIntegralsDriver
from .mointsdriver import MOIntegralsDriver
from .veloxchemlib import AODensityMatrix
from .veloxchemlib import AOFockMatrix
from .veloxchemlib import mpi_master
from .veloxchemlib import denmat
from .veloxchemlib import fockmat
from .veloxchemlib import XCIntegrator
from .orbitalresponse import OrbitalResponse
from .qqscheme import get_qq_scheme


class Mp2OrbitalResponse(OrbitalResponse):
    """
    Implements orbital response Lagrange multipliers computation using a
    conjugate gradient scheme for the MP2 level of theory.

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.
    """

    def __init__(self, comm, ostream):
        """
        Initializes orbital response computation driver to default setup.
        """

        super().__init__(comm, ostream)

    def update_settings(self, mp2_dict, method_dict=None):
        """
        Updates response and method settings in orbital response computation
        driver.

        :param mp2_dict:
            The dictionary of mp2 settings
            (essentially if conventional or distributed;
             only conventional for now)
        :param method_dict:
            The dictionary of method settings.
        """

        super().update_settings(method_dict=method_dict, mp2_dict=mp2_dict)

    def compute_rhs(self, molecule, basis, mol_orbs, dft_dict, profiler):
        """
        Computes the right-hand side (RHS) of the MP2 orbital response equation
        including the necessary density matrices using molecular data.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        :param mol_orbs:
            Molecular orbitals object from scfdriver.
        :param dft_dict:
            The dictionary containing DFT information.
            Kept for consistency with the parent class Orbital Response;
            will be an empty dictionary.
        :param profiler:
            The profiler.

        :return:
            A dictionary containing the orbital-response RHS and
            unrelaxed one-particle density.
        """

        profiler.start_timer(0, 'RHS')

        #TODO: consider removing scf_tensors; if everything we need is
        # contained in mol_orbs

        # Workflow:
        # 1) Construct the necessary density matrices
        # 2) Construct the RHS
        # 3) Construct the initial guess => in parent class
        # 4) Write the linear operator for matrix-vector product
        #    => in parent class
        # 5) Run the conjugate gradient => in parent class

        if self.rank == mpi_master():

            # 0) Preparing necessary quantities (MO coefficients, orbital energies)
            ### ovlp = scf_tensors['S'] #TODO
            mo = mol_orbs.alpha_to_numpy()

            nocc = molecule.number_of_alpha_electrons()
            mo_occ = mo[:, :nocc].copy()
            mo_vir = mo[:, nocc:].copy()
            nvir = mo_vir.shape[1]
            orb_ene = mol_orbs.ea_to_numpy()
            evir = orb_ene[nocc:]
            eocc = orb_ene[:nocc]
            evv = evir.reshape(-1,1) + evir
            eoo = eocc.reshape(-1,1) + eocc

            # Creating the 4D tensor of orbital energy differences
            eoovv = eoo.reshape((nocc,nocc,1,1)) - evv.reshape((1,1,nvir,nvir))


            # TODO: consider transforming the t-amplitudes to AO basis
            # t2_ao = np.linalg.multi_dot([mo_occ, t2_mo, mo_vir.T])
            # t2_ao = np.einsum('mi,nj,ijab,ta,pb->mntp', mo_occ, mo_occ, t2_mo, mo_vir, mo_vir)

            # Calculate the oovv integrals and anti-symmetrize them
            moints_drv = MOIntegralsDriver(self.comm, self.ostream)
            oovv = moints_drv.compute_in_mem(molecule, basis, mol_orbs, "OOVV")
            oovv_antisym = oovv - oovv.transpose(0,1,3,2)

            # TODO: check what is more efficient (memory and time)
            # the way it is implemented now, or using transpose directly inside
            # the np.einsum

            # 1) Calculate unrelaxed one-particle density matrix in MO basis
            dm_oo = - ( np.einsum('ikab,jkab->ij',
                                         oovv / eoovv,
                                        (oovv + oovv_antisym) / eoovv,
                                         optimize=True)
                            + np.einsum('ikab,jkab->ij',
                                        (oovv + oovv_antisym) / eoovv,
                                         oovv / eoovv, optimize=True)
                             )

            dm_vv =  ( np.einsum('ijac,ijbc->ab',
                                        oovv / eoovv,
                                       (oovv + oovv_antisym) / eoovv,
                                        optimize=True)
                           + np.einsum('ijac,ijbc->ab',
                                       (oovv + oovv_antisym) / eoovv,
                                        oovv / eoovv)
                            )

            # print("DM_OO")
            # print(dm_oo)

            # Transform unrelaxed one-particle density matrix to the AO basis
            unrel_dm_ao = (np.linalg.multi_dot([mo_occ, dm_oo, mo_occ.T]) +
                           np.linalg.multi_dot([mo_vir, dm_vv, mo_vir.T]))

            # 2) Construct the right-hand side
            dm_ao_rhs = AODensityMatrix([unrel_dm_ao], denmat.rest)
        else:
            dm_ao_rhs = AODensityMatrix()

        dm_ao_rhs.broadcast(self.rank, self.comm)

        molgrid = dft_dict['molgrid']
        gs_density = dft_dict['gs_density']

        # Fock matrices with corresponding type
        fock_ao_rhs = AOFockMatrix(dm_ao_rhs)

        eri_drv = ElectronRepulsionIntegralsDriver(self.comm)
        screening = eri_drv.compute(get_qq_scheme(self.qq_type),
                                    self.eri_thresh, molecule, basis)

        eri_drv.compute(fock_ao_rhs, dm_ao_rhs, molecule, basis, screening)

        fock_ao_rhs.reduce_sum(self.rank, self.nodes, self.comm)

        # Calculate the RHS and transform it to the MO basis
        if self.rank == mpi_master():
            fock_ao_rhs_0 = fock_ao_rhs.alpha_to_numpy(0)
            ### fock_ao_rhs_1 = fock_ao_rhs.alpha_to_numpy(1)

            fmo_rhs_0 = np.linalg.multi_dot(
                [mo_occ.T, 0.5 * fock_ao_rhs_0, mo_vir])

            # Compute the 2PDM contributions to the RHS:
            # TODO: is it worth doing this in AO basis, rather than MO?
            ooov = moints_drv.compute_in_mem(molecule, basis, mol_orbs, "OOOV")
            ovvv = moints_drv.compute_in_mem(molecule, basis, mol_orbs, "OVVV")
            ooov_antisym = ooov - ooov.transpose(1,0,2,3)
            ovvv_antisym = ovvv - ovvv.transpose(0,1,3,2)

            # Not sure about the "-" sign...
            rhs_2pdm_mo = -0.5*( np.einsum('jkab,jkib->ia',
                                           oovv / eoovv,
                                           ooov + ooov_antisym,
                                           optimize=True)
                               + np.einsum('jkba,kjib->ia',
                                           oovv / eoovv,
                                           ooov + ooov_antisym,
                                           optimize=True)
                               - np.einsum('ijbc,jacb->ia',
                                            oovv / eoovv,
                                            ovvv + ovvv_antisym,
                                            optimize=True)
                               - np.einsum('ijcb,jabc->ia',
                                            oovv / eoovv,
                                            ovvv + ovvv_antisym,
                                            optimize=True)
                               )

            rhs_mo = fmo_rhs_0 + rhs_2pdm_mo

        profiler.stop_timer(0, 'RHS')

        # TODO: delete print statements
        # print("1PDM contribution to RHS:\n")
        # print(fmo_rhs_0)
        # print()
        # print("2PDM contribution to RHS:\n")
        # print(rhs_2pdm_mo)
        # print("Right Hand SIDE:\n")
        # print(rhs_mo)

        if self.rank == mpi_master():
            return {
                'rhs_mo': rhs_mo,
                'dm_oo': dm_oo,
                'dm_vv': dm_vv,
                'unrel_dm_ao': unrel_dm_ao,
                'fock_ao_rhs': fock_ao_rhs,
            }
        else:
            return {}

    # TODO: fix omega; now it is from TDA;
    def compute_omega(self, ovlp, mo_occ, mo_vir, epsilon_dm_ao, mol_orbs,
                      fock_ao_rhs, fock_lambda):
        """
        Calculates the TDA Lagrange multipliers for the overlap matrix.

        :param ovlp:
            The overlap matrix.
        :param mo_occ:
            The occupied MO coefficients.
        :param mo_vir:
            The virtual MO coefficients.
        :param epsilon_dm_ao:
            The energy-weighted relaxed density matrix.
        :param mol_orbs:
            The molecular orbitals object from scfdriver.
        :param fock_ao_rhs:
            The AOFockMatrix from the right-hand side of the orbital response eq.
        :param fock_lambda:
            The Fock matrix from Lagrange multipliers.

        :return:
            a numpy array containing the Lagrange multipliers in AO basis.
        """

        # Get the excitation vector of interest and transform it to AO
        nocc = mo_occ.shape[1]
        nvir = mo_vir.shape[1]
        ###exc_vec = tda_results['eigenvectors'][:, self.n_state_deriv]
        ###exc_vec = exc_vec.reshape(nocc, nvir).copy()
        ###exc_vec_ao = np.linalg.multi_dot([mo_occ, exc_vec, mo_vir.T])

        # The density matrix; only alpha block;
        # Only works for the restricted case
        D_occ = np.matmul(mo_occ, mo_occ.T)
        D_vir = np.matmul(mo_vir, mo_vir.T)

        # Because the excitation vector is not symmetric,
        # we need both the matrix (OO block in omega, and probably VO)
        # and its transpose (VV, OV blocks)
        # this comes from the transformation of the 2PDM contribution
        # from MO to AO basis
        ### fock_ao_rhs_1 = fock_ao_rhs.alpha_to_numpy(1)
        ### Ft = np.linalg.multi_dot([0.5 * fock_ao_rhs_1.T, exc_vec_ao, ovlp])
        ### F = np.linalg.multi_dot([0.5 * fock_ao_rhs_1, exc_vec_ao.T, ovlp.T])

        # Compute the contributions from the 2PDM and the relaxed 1PDM
        # to the omega Lagrange multipliers:
        fmat = (fock_lambda.alpha_to_numpy(0) +
                fock_lambda.alpha_to_numpy(0).T +
                0.5 * fock_ao_rhs.alpha_to_numpy(0))

        omega_1pdm_2pdm_contribs = ( #np.linalg.multi_dot([D_occ, F, D_occ]) +
        #                            np.linalg.multi_dot([D_occ, Ft, D_vir]) +
        #                            np.linalg.multi_dot([D_occ, Ft, D_vir]).T +
        #                            np.linalg.multi_dot([D_vir, Ft, D_vir]) +
                                    np.linalg.multi_dot([D_occ, fmat, D_occ]))

        omega = -epsilon_dm_ao - omega_1pdm_2pdm_contribs

        return omega
