#
#                                   VELOXCHEM
#              ----------------------------------------------------
#                          An Electronic Structure Code
#
#  SPDX-License-Identifier: BSD-3-Clause
#
#  Copyright 2018-2025 VeloxChem developers
#
#  Redistribution and use in source and binary forms, with or without modification,
#  are permitted provided that the following conditions are met:
#
#  1. Redistributions of source code must retain the above copyright notice, this
#     list of conditions and the following disclaimer.
#  2. Redistributions in binary form must reproduce the above copyright notice,
#     this list of conditions and the following disclaimer in the documentation
#     and/or other materials provided with the distribution.
#  3. Neither the name of the copyright holder nor the names of its contributors
#     may be used to endorse or promote products derived from this software without
#     specific prior written permission.
#
#  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
#  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
#  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
#  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
#  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
#  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
#  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
#  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
#  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
#  OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

from collections import Counter
import numpy as np
import time as tm
import math

from .veloxchemlib import T4CScreener
from .veloxchemlib import XCMolecularHessian
from .veloxchemlib import mpi_master
from .veloxchemlib import make_matrix, mat_t, partition_atoms
from .veloxchemlib import (OverlapGeom100Driver, KineticEnergyGeom100Driver,
                           NuclearPotentialGeom100Driver,
                           NuclearPotentialGeom010Driver, FockGeom1000Driver)
from .matrices import Matrices
from .matrix import Matrix
from .profiler import Profiler
from .distributedarray import DistributedArray
from .cphfsolver import CphfSolver
from .errorhandler import assert_msg_critical
from .dftutils import get_default_grid_level
from .batchsize import get_batch_size


class UnrestrictedHessianOrbitalResponse(CphfSolver):
    """
    Implements solver for the coupled-perturbed Hartree-Fock (CPHF) equations
    for the unrestricted Hessian orbital response

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.

    Instance variables
        - use_subspace_solver: flag to use subspace solver
          instead of conjugate gradient.
    """

    def __init__(self, comm=None, ostream=None):
        """
        Initializes orbital response driver to default setup.

        """

        super().__init__(comm, ostream)

        self.orbrsp_type = 'hessian'
        self._embedding_hess_drv = None

    def update_settings(self, cphf_dict, method_dict=None):
        """
        Updates response and method settings in CPHF solver.

        :param cphf_dict:
            The dictionary of CPHF (orbital response) settings.
        :param method_dict:
            The dictionary of method settings.
        """

        super().update_settings(cphf_dict, method_dict)

    def compute_rhs(self,
                    molecule,
                    basis,
                    scf_tensors,
                    eri_dict,
                    dft_dict,
                    pe_dict,
                    atom_pairs=None):
        """
        Computes the right hand side for the CPHF equations for
        the analytical Hessian, all atomic coordinates.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        :param scf_tensors:
            The tensors from the converged SCF calculation.

        :returns:
            The RHS of the CPHF equations.
        """

        profiler = Profiler({
            'timing': self.timing,
            'profiling': self.profiling,
            'memory_profiling': self.memory_profiling,
            'memory_tracing': self.memory_tracing,
        })

        profiler.set_timing_key('RHS')

        natm = molecule.number_of_atoms()
        nocc_a = molecule.number_of_alpha_electrons()
        nocc_b = molecule.number_of_beta_electrons()

        if self.rank == mpi_master():
            density_a = scf_tensors['D_alpha']
            density_b = scf_tensors['D_beta']
            eocc_a = scf_tensors['E_alpha'][:nocc_a]
            eocc_b = scf_tensors['E_beta'][:nocc_b]
            mo_a = scf_tensors['C_alpha']
            mo_b = scf_tensors['C_beta']
            point_charges = scf_tensors.get('point_charges', None)
            qm_vdw_params = scf_tensors.get('qm_vdw_params', None)
        else:
            density_a = None
            density_b = None
            eocc_a = None
            eocc_b = None
            mo_a = None
            mo_b = None
            point_charges = None
            qm_vdw_params = None

        density_a = self.comm.bcast(density_a, root=mpi_master())
        eocc_a = self.comm.bcast(eocc_a, root=mpi_master())
        mo_a = self.comm.bcast(mo_a, root=mpi_master())
        density_b = self.comm.bcast(density_b, root=mpi_master())
        eocc_b = self.comm.bcast(eocc_b, root=mpi_master())
        mo_b = self.comm.bcast(mo_b, root=mpi_master())
        point_charges = self.comm.bcast(point_charges, root=mpi_master())
        qm_vdw_params = self.comm.bcast(qm_vdw_params, root=mpi_master())

        # TODO: double check dft_dict['gs_density']
        mol_grid = dft_dict['molgrid']
        gs_density = [density_a, density_b]

        # MO coefficients
        mo_occ_a = mo_a[:, :nocc_a].copy()
        mo_vir_a = mo_a[:, nocc_a:].copy()
        nvir_a = mo_vir_a.shape[1]
        mo_occ_b = mo_b[:, :nocc_b].copy()
        mo_vir_b = mo_b[:, nocc_b:].copy()
        nvir_b = mo_vir_b.shape[1]

        omega_ao_a = -1.0 * np.linalg.multi_dot([mo_occ_a, np.diag(eocc_a), mo_occ_a.T])
        omega_ao_b = -1.0 * np.linalg.multi_dot([mo_occ_b, np.diag(eocc_b), mo_occ_b.T])

        # partition atoms for parallellisation
        # TODO: use partition_atoms in e.g. scfgradientdriver

        if atom_pairs is None:
            local_atoms = partition_atoms(natm, self.rank, self.nodes)
        else:
            atoms_in_pairs = []
            for i, j in atom_pairs:
                if i not in atoms_in_pairs:
                    atoms_in_pairs.append(i)
                if j not in atoms_in_pairs:
                    atoms_in_pairs.append(j)
            natm_in_pairs = len(atoms_in_pairs)
            local_atoms = []
            # todo atompairs, is this the best way to go about the ordering here?
            for at in atoms_in_pairs:
                if at % self.nodes == self.rank:
                    local_atoms.append(at)

        # Gathers information of which rank has which atom,
        # and then broadcasts this to all ranks
        atom_idx_rank = [(iatom, self.rank) for iatom in local_atoms]
        gathered_atom_idx_rank = self.comm.gather(atom_idx_rank)
        if self.rank == mpi_master():
            all_atom_idx_rank = []
            for local_list in gathered_atom_idx_rank:
                for elem in local_list:
                    all_atom_idx_rank.append(elem)
            all_atom_idx_rank = sorted(all_atom_idx_rank)
        else:
            all_atom_idx_rank = None
        all_atom_idx_rank = self.comm.bcast(all_atom_idx_rank,
                                            root=mpi_master())

        # preparing the CPHF RHS

        t0 = tm.time()

        # overlap derivatives

        profiler.start_timer('dOvlp')

        ovlp_deriv_ao_dict = {
            (iatom, x): None for iatom in range(natm) for x in range(3)
        }

        ovlp_grad_drv = OverlapGeom100Driver()

        for iatom in local_atoms:

            gmats = ovlp_grad_drv.compute(molecule, basis, iatom)

            for x, label in enumerate(['X', 'Y', 'Z']):
                gmat = gmats.matrix_to_numpy(label)
                ovlp_deriv_ao_dict[(iatom, x)] = gmat + gmat.T

            gmats = Matrices()

        dist_ovlp_deriv_ao = {}

        for iatom, root_rank in all_atom_idx_rank:
            for x in range(3):
                key_ix = (iatom, x)

                dist_ovlp_deriv_ao[key_ix] = DistributedArray(
                    ovlp_deriv_ao_dict[key_ix], self.comm, root=root_rank)

                ovlp_deriv_ao_dict[key_ix] = None

        ovlp_deriv_ao_dict.clear()

        profiler.stop_timer('dOvlp')

        # Fock derivatives

        profiler.start_timer('dFock')

        fock_deriv_ao_dict_a = {
            (iatom, x): None for iatom in range(natm) for x in range(3) 
        }
        fock_deriv_ao_dict_b = {
            (iatom, x): None for iatom in range(natm) for x in range(3) 
        }

        for iatom in local_atoms:
            # TODO: the derivatives are saved in a tuple. Is this a good idea?
            fock_deriv_ao_i = self._compute_fmat_deriv_unrestricted(molecule, basis,
                                                       density_a, density_b,
                                                       iatom, eri_dict,
                                                       point_charges,
                                                       qm_vdw_params)

            fock_deriv_ao_dict_a[(iatom, 0)] = fock_deriv_ao_i[0][0]
            fock_deriv_ao_dict_a[(iatom, 1)] = fock_deriv_ao_i[0][1]
            fock_deriv_ao_dict_a[(iatom, 2)] = fock_deriv_ao_i[0][2]

            fock_deriv_ao_dict_b[(iatom, 0)] = fock_deriv_ao_i[1][0]
            fock_deriv_ao_dict_b[(iatom, 1)] = fock_deriv_ao_i[1][1]
            fock_deriv_ao_dict_b[(iatom, 2)] = fock_deriv_ao_i[1][2]

        dist_fock_deriv_ao_a = {}
        dist_fock_deriv_ao_b = {}

        for iatom, root_rank in all_atom_idx_rank:
            for x in range(3):
                key_ix = (iatom, x)

                dist_fock_deriv_ao_a[key_ix] = DistributedArray(
                    fock_deriv_ao_dict_a[key_ix], self.comm, root=root_rank)

                fock_deriv_ao_dict_a[key_ix] = None

                dist_fock_deriv_ao_b[key_ix] = DistributedArray(
                    fock_deriv_ao_dict_b[key_ix], self.comm, root=root_rank)

                fock_deriv_ao_dict_b[key_ix] = None

        fock_deriv_ao_dict_a.clear()
        fock_deriv_ao_dict_b.clear()

        profiler.stop_timer('dFock')

        # XC derivatives

        # Note: Parallelization of DFT integration is done over grid points
        # instead of atoms.

        assert_msg_critical(
            not self._dft, 'UnrestrictedHessianOrbitalResponse: ' +
            'DFT not yet implemented.')

        #profiler.start_timer('dXC')

        #if self._dft:
        #    xc_mol_hess = XCMolecularHessian()

        #    if atom_pairs is None:
        #        naos = basis.get_dimensions_of_basis()

        #        batch_size = get_batch_size(None, natm * 3, naos, self.comm)
        #        batch_size = batch_size // 3

        #        num_batches = natm // batch_size
        #        if natm % batch_size != 0:
        #            num_batches += 1
        #    else:
        #        ao_map = basis.get_ao_basis_map(molecule)
        #        # Create a dictionary of the amount of AOs per atom
        #        aos_per_atom = dict(
        #            Counter(int(ao.strip().split()[0]) - 1 for ao in ao_map))
        #        # Calculate total number of AOs for atoms in atoms_in_pairs
        #        naos_in_pairs = sum(
        #            aos_per_atom.get(atom, 0) for atom in atoms_in_pairs)
        #        batch_size = get_batch_size(None, natm_in_pairs * 3,
        #                                    naos_in_pairs, self.comm)
        #        batch_size = batch_size // 3

        #        num_batches = natm_in_pairs // batch_size
        #        if natm % batch_size != 0:
        #            num_batches += 1

        #    for batch_ind in range(num_batches):

        #        if atom_pairs is None:
        #            batch_start = batch_ind * batch_size
        #            batch_end = min(batch_start + batch_size, natm)

        #            atom_list = list(range(batch_start, batch_end))
        #        else:
        #            batch_start = batch_ind * batch_size
        #            batch_end = min(batch_start + batch_size, natm_in_pairs)

        #            atom_list = atoms_in_pairs[batch_start:batch_end]

        #        vxc_deriv_batch = xc_mol_hess.integrate_vxc_fock_gradient(
        #            molecule, basis, gs_density, mol_grid,
        #            self.xcfun.get_func_label(), atom_list)

        #        for vecind, (iatom, root_rank) in enumerate(
        #                all_atom_idx_rank[batch_start:batch_end]):

        #            for x in range(3):
        #                vxc_deriv_ix = self.comm.reduce(
        #                    vxc_deriv_batch[vecind * 3 + x])

        #                dist_vxc_deriv_ix = DistributedArray(
        #                    vxc_deriv_ix, self.comm)

        #                key_ix = (iatom, x)
        #                dist_fock_deriv_ao[
        #                    key_ix].data += dist_vxc_deriv_ix.data

        #profiler.stop_timer('dXC')

        t1 = tm.time()

        if self.rank == mpi_master():
            self.ostream.print_info("CPHF/CPKS integral derivatives computed" +
                                    ' in {:.2f} sec.'.format(t1 - t0))
            self.ostream.print_blank()
            self.ostream.flush()

        hessian_first_integral_derivatives = np.zeros((natm, 3, natm, 3))

        profiler.start_timer('1stHess')

        for iatom, root_rank_i in all_atom_idx_rank:
            for x in range(3):
                key_ix = (iatom, x)

                ovlp_deriv_ao_ix = dist_ovlp_deriv_ao[key_ix].get_full_matrix()
                fock_deriv_ao_ix_a = dist_fock_deriv_ao_a[key_ix].get_full_matrix()
                fock_deriv_ao_ix_b = dist_fock_deriv_ao_b[key_ix].get_full_matrix()

                if self.rank == mpi_master():
                    DFD_ix_a = np.linalg.multi_dot(
                        [density_a, fock_deriv_ao_ix_a, density_a])
                    DFD_ix_b = np.linalg.multi_dot(
                        [density_b, fock_deriv_ao_ix_b, density_b])
                    DSD_ix_a = np.linalg.multi_dot(
                        [density_a, ovlp_deriv_ao_ix, density_a])
                    DSD_ix_b = np.linalg.multi_dot(
                        [density_b, ovlp_deriv_ao_ix, density_b])
                    OmegaSD_ix_a = np.linalg.multi_dot(
                        [omega_ao_a, ovlp_deriv_ao_ix, density_a])
                    OmegaSD_ix_b = np.linalg.multi_dot(
                        [omega_ao_b, ovlp_deriv_ao_ix, density_b])
                else:
                    DFD_ix_a = None
                    DSD_ix_a = None
                    OmegaSD_ix_a = None
                    DFD_ix_b = None
                    DSD_ix_b = None
                    OmegaSD_ix_b = None

                dist_DFD_ix_a = DistributedArray(DFD_ix_a, self.comm)
                dist_DSD_ix_a = DistributedArray(DSD_ix_a, self.comm)
                dist_OmegaSD_ix_a = DistributedArray(OmegaSD_ix_a, self.comm)
                dist_DFD_ix_b = DistributedArray(DFD_ix_b, self.comm)
                dist_DSD_ix_b = DistributedArray(DSD_ix_b, self.comm)
                dist_OmegaSD_ix_b = DistributedArray(OmegaSD_ix_b, self.comm)

                for jatom, root_rank_j in all_atom_idx_rank:
                    if jatom < iatom:
                        continue

                    if atom_pairs is not None:
                        if ((iatom, jatom) not in atom_pairs and
                            (jatom, iatom) not in atom_pairs and
                                iatom != jatom):
                            continue
                    # These are the contributions from the OO block of the
                    # CPHF coefficients/orbital rsp. multipliers to the Hessian.
                    # uij = -1/2 Sij
                    # TODO: Check if these terms are correct!
                    for y in range(3):
                        key_jy = (jatom, y)

                        Fix_Sjy_a = np.dot(
                            dist_DFD_ix_a.data.reshape(-1),
                            dist_ovlp_deriv_ao[key_jy].data.reshape(-1))

                        Fjy_Six_a = np.dot(
                            dist_DSD_ix_a.data.reshape(-1),
                            dist_fock_deriv_ao_a[key_jy].data.reshape(-1))

                        Six_Sjy_a =  np.dot(
                            dist_OmegaSD_ix_a.data.reshape(-1),
                            dist_ovlp_deriv_ao[key_jy].data.reshape(-1))

                        Fix_Sjy_b = np.dot(
                            dist_DFD_ix_b.data.reshape(-1),
                            dist_ovlp_deriv_ao[key_jy].data.reshape(-1))

                        Fjy_Six_b = np.dot(
                            dist_DSD_ix_b.data.reshape(-1),
                            dist_fock_deriv_ao_b[key_jy].data.reshape(-1))

                        Six_Sjy_b =  2.0 * np.dot(
                            dist_OmegaSD_ix_b.data.reshape(-1),
                            dist_ovlp_deriv_ao[key_jy].data.reshape(-1))

                        hess_ijxy = -1.0 * (Fix_Sjy_a + Fjy_Six_a + Six_Sjy_a
                                          + Fix_Sjy_b + Fjy_Six_b + Six_Sjy_b)
                        hessian_first_integral_derivatives[iatom, x, jatom,
                                                           y] += hess_ijxy
                        if iatom != jatom:
                            hessian_first_integral_derivatives[jatom, y, iatom,
                                                               x] += hess_ijxy

        hessian_first_integral_derivatives = self.comm.reduce(
            hessian_first_integral_derivatives, root=mpi_master())

        profiler.stop_timer('1stHess')

        dist_cphf_rhs = []

        hessian_eri_overlap = np.zeros((natm, 3, natm, 3))

        for iatom, root_rank in all_atom_idx_rank:

            uij_t0 = tm.time()

            # the oo part of the CPHF coefficients in AO basis
            if self.rank == mpi_master():
                uij_ao_list_a = []
                uij_ao_list_b = []
            else:
                uij_ao_list_a = None
                uij_ao_list_b = None

            for x in range(3):
                key_ix = (iatom, x)
                ovlp_deriv_ao_ix = dist_ovlp_deriv_ao[key_ix].get_full_matrix()

                if self.rank == mpi_master():
                    uij_ao_list_a.append(
                        np.linalg.multi_dot(
                            [density_a, -0.5 * ovlp_deriv_ao_ix, density_a]))
                    uij_ao_list_b.append(
                        np.linalg.multi_dot(
                            [density_b, -0.5 * ovlp_deriv_ao_ix, density_b]))

            profiler.add_timing_info('UijAO', tm.time() - uij_t0)

            # create AODensity and Fock matrix objects, contract with ERI
            # This term contributes to the RHS of the CPHF/CPKS equations.
            # fock_uij is an array [(ax, bx), (ay, by), (az, bz)] where
            # a is an alpha component and b is beta.
            fock_uij = self._comp_lr_fock_unrestricted(
                    (uij_ao_list_a, uij_ao_list_b), molecule, basis,
                     eri_dict, dft_dict, pe_dict, profiler)

            if self.rank == mpi_master():
                uij_ao_list_a.clear()
                uij_ao_list_b.clear()

            uij_t0 = tm.time()

            for x in range(3):
                if self.rank == mpi_master():
                    DFD_ix_a = np.linalg.multi_dot(
                        [density_a, fock_uij[x * 2 + 0], density_a])
                    DFD_ix_b = np.linalg.multi_dot(
                        [density_b, fock_uij[x * 2 + 1], density_b])
                else:
                    DFD_ix_a = None
                    DFD_ix_b = None

                dist_DFD_ix_a = DistributedArray(DFD_ix_a, self.comm)
                dist_DFD_ix_b = DistributedArray(DFD_ix_b, self.comm)

                # Note: hessian_eri_overlap, i.e. inner product of P_P_Six_fock_ao and
                # P_P_Sjy, is obtained from fock_uij and uij_ao_jy in hessian orbital
                # response
                # TODO: Check if these are correct!
                for jatom, root_rank_j in all_atom_idx_rank:
                    if jatom < iatom:
                        continue

                    if atom_pairs is not None:
                        if ((iatom, jatom) not in atom_pairs and
                            (jatom, iatom) not in atom_pairs and
                                iatom != jatom):
                            continue

                    for y in range(3):
                        key_jy = (jatom, y)

                        # 2.0 * (np.dot(DFD_ix_list[x], (-0.5 * ovlp_deriv_ao_jy)))
                        hess_ijxy_a = -np.dot(
                            dist_DFD_ix_a.data.reshape(-1),
                            dist_ovlp_deriv_ao[key_jy].data.reshape(-1))
                        hess_ijxy_b = -np.dot(
                            dist_DFD_ix_b.data.reshape(-1),
                            dist_ovlp_deriv_ao[key_jy].data.reshape(-1))

                        hess_ijxy_a *= 2.0
                        hess_ijxy_b *= 2.0

                        hessian_eri_overlap[iatom, x, jatom, y] += hess_ijxy_a + hess_ijxy_b
                        if iatom != jatom:
                            hessian_eri_overlap[jatom, y, iatom, x] += hess_ijxy_a + hess_ijxy_b

            profiler.add_timing_info('UijDot', tm.time() - uij_t0)

            uij_t0 = tm.time()

            for x in range(3):

                # form dist_fock_uij_ov_2 from mpi_master

                if self.rank == mpi_master():
                    # transform to MO basis
                    fock_uij_ov_2_a =  np.linalg.multi_dot(
                        [mo_occ_a.T, fock_uij[x * 2 + 0], mo_vir_a])
                    fock_uij_ov_2_a = fock_uij_ov_2_a.reshape(nocc_a * nvir_a)
                    fock_uij_ov_2_b = np.linalg.multi_dot(
                        [mo_occ_b.T, fock_uij[x * 2 + 1], mo_vir_b])
                    fock_uij_ov_2_b = fock_uij_ov_2_b.reshape(nocc_b * nvir_b)
                else:
                    fock_uij_ov_2_a = None
                    fock_uij_ov_2_b = None

                dist_fock_uij_ov_2_a = DistributedArray(fock_uij_ov_2_a, self.comm)
                dist_fock_uij_ov_2_b = DistributedArray(fock_uij_ov_2_b, self.comm)

                # form dist_fock_deriv_ov_ix and dist_orben_ovlp_deriv_ov_ix
                # from root_rank

                key_ix = (iatom, x)
                ovlp_deriv_ao_ix = dist_ovlp_deriv_ao[key_ix].get_full_matrix(
                    root=root_rank)
                fock_deriv_ao_ix_a = dist_fock_deriv_ao_a[key_ix].get_full_matrix(
                    root=root_rank)
                fock_deriv_ao_ix_b = dist_fock_deriv_ao_b[key_ix].get_full_matrix(
                    root=root_rank)

                # TODO: do this on mpi_master
                if self.rank == root_rank:

                    fock_deriv_ov_dict_ix_a = np.linalg.multi_dot(
                        [mo_occ_a.T, fock_deriv_ao_ix_a,
                         mo_vir_a]).reshape(nocc_a * nvir_a)

                    ovlp_deriv_ov_ix_a = np.linalg.multi_dot(
                        [mo_occ_a.T, ovlp_deriv_ao_ix, mo_vir_a])

                    orben_ovlp_deriv_ov_ix_a = (eocc_a.reshape(-1, 1) *
                                              ovlp_deriv_ov_ix_a).reshape(nocc_a *
                                                                          nvir_a)
                    fock_deriv_ov_dict_ix_b = np.linalg.multi_dot(
                        [mo_occ_b.T, fock_deriv_ao_ix_b,
                         mo_vir_b]).reshape(nocc_b * nvir_b)

                    ovlp_deriv_ov_ix_b = np.linalg.multi_dot(
                        [mo_occ_b.T, ovlp_deriv_ao_ix, mo_vir_b])

                    orben_ovlp_deriv_ov_ix_b = (eocc_b.reshape(-1, 1) *
                                              ovlp_deriv_ov_ix_b).reshape(nocc_b *
                                                                          nvir_b)

                else:
                    fock_deriv_ov_dict_ix_a = None
                    orben_ovlp_deriv_ov_ix_a = None
                    fock_deriv_ov_dict_ix_b = None
                    orben_ovlp_deriv_ov_ix_b = None

                dist_fock_deriv_ov_ix_a = DistributedArray(fock_deriv_ov_dict_ix_a,
                                                           self.comm,
                                                           root=root_rank)

                dist_orben_ovlp_deriv_ov_ix_a = DistributedArray(
                    orben_ovlp_deriv_ov_ix_a, self.comm, root=root_rank)

                dist_cphf_rhs_ix_data_a = (dist_fock_uij_ov_2_a.data +
                                         dist_fock_deriv_ov_ix_a.data -
                                         dist_orben_ovlp_deriv_ov_ix_a.data)

                dist_fock_deriv_ov_ix_b = DistributedArray(fock_deriv_ov_dict_ix_b,
                                                           self.comm,
                                                           root=root_rank)

                dist_orben_ovlp_deriv_ov_ix_b = DistributedArray(
                    orben_ovlp_deriv_ov_ix_b, self.comm, root=root_rank)

                dist_cphf_rhs_ix_data_b = (dist_fock_uij_ov_2_b.data +
                                         dist_fock_deriv_ov_ix_b.data -
                                         dist_orben_ovlp_deriv_ov_ix_b.data)

                dist_cphf_rhs_ix_a = DistributedArray(dist_cphf_rhs_ix_data_a,
                                                    self.comm,
                                                    distribute=False)
                dist_cphf_rhs_ix_b = DistributedArray(dist_cphf_rhs_ix_data_b,
                                                    self.comm,
                                                    distribute=False)

                dist_cphf_rhs.append(dist_cphf_rhs_ix_a)
                dist_cphf_rhs.append(dist_cphf_rhs_ix_b)

            profiler.add_timing_info('distRHS', tm.time() - uij_t0)

        # fill up the missing values in dist_cphf_rhs with zeros
        # so later indexing does not fail
        # TODO: check if this works correctly.
        if atom_pairs is not None:
            for i in range(molecule.number_of_atoms()):
                if i not in atoms_in_pairs:
                    zer_a = np.zeros(len(dist_cphf_rhs[0].data))
                    zer_b = np.zeros(len(dist_cphf_rhs[1].data))
                    for j in range(3):
                        empty_dist_arr_a = DistributedArray(zer_a,
                                                          self.comm,
                                                          distribute=False)
                        empty_dist_arr_b = DistributedArray(zer_b,
                                                          self.comm,
                                                          distribute=False)
                        # 2 * index + 0 corresponds to alpha
                        dist_cphf_rhs.insert(2 * (i * 3 + j), empty_dist_arr_a)
                        # 2 * index + 1 corresponds to beta
                        dist_cphf_rhs.insert(2 * (i * 3 + j) + 1, empty_dist_arr_b)

        hessian_eri_overlap = self.comm.reduce(hessian_eri_overlap,
                                               root=mpi_master())

        t2 = tm.time()

        if self.rank == mpi_master():
            self.ostream.print_info('CPHF/CPKS right-hand side computed' +
                                    ' in {:.2f} sec.'.format(t2 - t1))
            self.ostream.print_blank()
            self.ostream.flush()

        profiler.print_timing(self.ostream)
        profiler.print_profiling_summary(self.ostream)

        return {
            'dist_cphf_rhs': dist_cphf_rhs,
            # 'ovlp_deriv_oo': ovlp_deriv_oo,
            # 'fock_deriv_ao': fock_deriv_ao,
            # 'fock_uij': fock_uij,
            'hessian_first_integral_derivatives': hessian_first_integral_derivatives,
            'hessian_eri_overlap': hessian_eri_overlap,
        }

    def _compute_fmat_deriv_unrestricted(self,
                            molecule,
                            basis,
                            density_a,
                            density_b,
                            i,
                            eri_dict,
                            point_charges=None,
                            qm_vdw_params=None):
        """
        Computes the derivative of the Fock matrix with respect
        to the coordinates of atom i.

        :param molecule:
            The molecule.
        :param basis:
            The basis set.
        :param density_a:
            The density matrix in AO basis corresponding to spin alpha.
        :param density_b:
            The density matrix in AO basis corresponding to spin beta.
        :param i:
            The atom index.
        :param eri_dict:
            The dictionary containing ERI information.

        :return fmat_deriv:
            The derivative of the Fock matrix wrt. atom i.
        """

        # number of atomic orbitals
        nao = basis.get_dimensions_of_basis()

        # initialize fmat variable
        fmat_deriv_a = np.zeros((3, nao, nao))
        fmat_deriv_b = np.zeros((3, nao, nao))

        # kinetic integral gadient
        kin_grad_drv = KineticEnergyGeom100Driver()

        gmats_kin = kin_grad_drv.compute(molecule, basis, i)

        for x, label in enumerate(['X', 'Y', 'Z']):
            gmat_kin = gmats_kin.matrix_to_numpy(label)
            fmat_deriv_a[x] += gmat_kin + gmat_kin.T
            fmat_deriv_b[x] += gmat_kin + gmat_kin.T

        gmats_kin = Matrices()

        # nuclear potential integral gradients
        npot_grad_100_drv = NuclearPotentialGeom100Driver()
        npot_grad_010_drv = NuclearPotentialGeom010Driver()

        gmats_npot_100 = npot_grad_100_drv.compute(molecule, basis, i)
        gmats_npot_010 = npot_grad_010_drv.compute(molecule, basis, i)

        for x, label in enumerate(['X', 'Y', 'Z']):
            gmat_npot_100 = gmats_npot_100.matrix_to_numpy(label)
            gmat_npot_010 = gmats_npot_010.matrix_to_numpy(label)
            fmat_deriv_a[x] -= gmat_npot_100 + gmat_npot_100.T + gmat_npot_010
            fmat_deriv_b[x] -= gmat_npot_100 + gmat_npot_100.T + gmat_npot_010

        gmats_npot_100 = Matrices()
        gmats_npot_010 = Matrices()

        # point charges contribution
        if point_charges is not None:
            npoints = point_charges.shape[1]

            mm_coords = []
            mm_charges = []
            for p in range(npoints):
                xyz_p = point_charges[:3, p]
                chg_p = point_charges[3, p]
                mm_coords.append(xyz_p.copy())
                mm_charges.append(chg_p)

            gmats_100 = npot_grad_100_drv.compute(molecule, basis, i, mm_coords,
                                                  mm_charges)

            for x, label in enumerate(['X', 'Y', 'Z']):
                gmat_100 = gmats_100.matrix_to_numpy(label)
                fmat_deriv_a[x] -= gmat_100 + gmat_100.T
                fmat_deriv_b[x] -= gmat_100 + gmat_100.T

            gmats_100 = Matrices()

        # TODO: check if this actually works in the unrestricted case
        if self._embedding_hess_drv is not None:
            pe_fock_grad_contr = (
                self._embedding_hess_drv.compute_pe_fock_gradient_contributions(
                    i=i))
            for x in range(3):
                fmat_deriv_a[x] += pe_fock_grad_contr[x]
                fmat_deriv_b[x] += pe_fock_grad_contr[x]

        if self._dft:
            if self.xcfun.is_hybrid():
                fock_type = '2jkx'
                exchange_scaling_factor = self.xcfun.get_frac_exact_exchange()
            else:
                fock_type = 'j'
                exchange_scaling_factor = 0.0
        else:
            exchange_scaling_factor = 1.0
            fock_type = "2jk"

        need_omega = (self._dft and self.xcfun.is_range_separated())
        if need_omega:
            exchange_scaling_factor = (self.xcfun.get_rs_alpha() +
                                       self.xcfun.get_rs_beta())
            erf_k_coef = -self.xcfun.get_rs_beta()
            omega = self.xcfun.get_rs_omega()
        else:
            erf_k_coef, omega = None, None

        Da_for_fock = make_matrix(basis, mat_t.symmetric)
        Da_for_fock.set_values(density_a)

        Db_for_fock = make_matrix(basis, mat_t.symmetric)
        Db_for_fock.set_values(density_b)

        Dab_for_fock = make_matrix(basis, mat_t.symmetric)
        Dab_for_fock.set_values(density_a + density_b)

        # ERI threshold
        thresh_int = int(-math.log10(self.eri_thresh))

        # screening
        screener = eri_dict['screening']

        screener_atom = T4CScreener()
        screener_atom.partition_atom(basis, molecule, 'eri', i)

        # Fock gradient
        fock_grad_drv = FockGeom1000Driver()
        # TODO: should we set_block_size_factor?

        gmats_eri_Jab = fock_grad_drv.compute(basis, screener_atom, screener,
                                          Dab_for_fock, i, 'j',
                                          0.0, 0.0, thresh_int)

        if fock_type != 'j':
            gmats_eri_Ka = fock_grad_drv.compute(basis, screener_atom, screener,
                                          Da_for_fock, i, 'kx', exchange_scaling_factor,
                                          0.0, thresh_int)
            gmats_eri_Kb = fock_grad_drv.compute(basis, screener_atom, screener,
                                          Db_for_fock, i, 'kx', exchange_scaling_factor,
                                          0.0, thresh_int)

        if need_omega:
            # for range-separated functional
            gmats_eri_rs_a = fock_grad_drv.compute(basis, screener_atom, screener,
                                                  Da_for_fock, i, 'kx_rs',
                                                  erf_k_coef, omega, thresh_int)
            gmats_eri_rs_b = fock_grad_drv.compute(basis, screener_atom, screener,
                                                  Db_for_fock, i, 'kx_rs',
                                                  erf_k_coef, omega, thresh_int)

        # calculate gradient contributions
        for x, label in enumerate(['X', 'Y', 'Z']):
            gmat_eri_Jab = gmats_eri_Jab.matrix_to_numpy(label)
            gmat_eri_Ka = gmats_eri_Ka.matrix_to_numpy(label)
            gmat_eri_Kb = gmats_eri_Kb.matrix_to_numpy(label)
            fmat_deriv_a[x] += gmat_eri_Jab - gmat_eri_Ka
            fmat_deriv_b[x] += gmat_eri_Jab - gmat_eri_Kb

            if need_omega:
                # range-separated functional contribution
                gmat_eri_rs_a = gmats_eri_rs_a.matrix_to_numpy(label)
                gmat_eri_rs_b = gmats_eri_rs_b.matrix_to_numpy(label)
                fmat_deriv_a[x] -= gmat_eri_rs_a
                fmat_deriv_b[x] -= gmat_eri_rs_b

        gmats_eri = Matrices()
        Da_for_fock = Matrix()
        Db_for_fock = Matrix()
        Dab_for_fock = Matrix()

        return (fmat_deriv_a, fmat_deriv_b)

    # TODO: move to parent class?
    def print_cphf_header(self, title):
        """
        Prints information on the solver setup
        """

        self.ostream.print_blank()
        self.ostream.print_header('{:s} Setup'.format(title))
        self.ostream.print_header('=' * (len(title) + 8))
        self.ostream.print_blank()

        str_width = 70

        # print general info
        cur_str = 'Solver Type                     : '
        if self.use_subspace_solver:
            cur_str += 'Iterative Subspace Algorithm'
        else:
            cur_str += 'Conjugate Gradient'
        self.ostream.print_header(cur_str.ljust(str_width))

        cur_str = 'Max. Number of Iterations       : ' + str(self.max_iter)
        self.ostream.print_header(cur_str.ljust(str_width))
        cur_str = 'Convergence Threshold           : {:.1e}'.format(
            self.conv_thresh)
        self.ostream.print_header(cur_str.ljust(str_width))

        if self._dft:
            cur_str = 'Exchange-Correlation Functional : '
            cur_str += self.xcfun.get_func_label().upper()
            self.ostream.print_header(cur_str.ljust(str_width))
            grid_level = (get_default_grid_level(self.xcfun)
                          if self.grid_level is None else self.grid_level)
            cur_str = 'Molecular Grid Level            : ' + str(grid_level)
            self.ostream.print_header(cur_str.ljust(str_width))

        self.ostream.print_blank()
        self.ostream.flush()
