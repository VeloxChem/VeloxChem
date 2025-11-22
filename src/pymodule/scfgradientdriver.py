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

from mpi4py import MPI
from copy import deepcopy
import numpy as np
import time
import math

from .veloxchemlib import (OverlapGeom100Driver, KineticEnergyGeom100Driver,
                           NuclearPotentialGeom100Driver,
                           NuclearPotentialGeom010Driver, FockGeom1000Driver)
from .veloxchemlib import XCFunctional, MolecularGrid, XCMolecularGradient
from .veloxchemlib import T4CScreener
from .veloxchemlib import RIFockGradDriver
from .veloxchemlib import mpi_master, mat_t
from .veloxchemlib import make_matrix
from .veloxchemlib import parse_xc_func
from .veloxchemlib import bohr_in_angstrom, hartree_in_kjpermol
from .molecularbasis import MolecularBasis
from .matrices import Matrices
from .outputstream import OutputStream
from .dispersionmodel import DispersionModel
from .gradientdriver import GradientDriver
from .errorhandler import assert_msg_critical
from .sanitychecks import (molecule_sanity_check, scf_results_sanity_check,
                           dft_sanity_check)


class ScfGradientDriver(GradientDriver):
    """
    Implements SCF gradient driver.

    :param scf_drv:
        The SCF driver.

    Instance variables
        - flag: The driver flag.
        - delta_h: The displacement for finite difference.
        - dispersion: The flag for calculating D4 dispersion correction.
    """

    def __init__(self, scf_drv):
        """
        Initializes gradient driver.
        """

        super().__init__(scf_drv.comm, scf_drv.ostream)

        self.scf_driver = scf_drv
        self.flag = 'SCF Gradient Driver'

        self.eri_thresh = scf_drv.eri_thresh
        self.timing = scf_drv.timing
        self._debug = scf_drv._debug

        self._block_size_factor = 4

        self._xcfun_ldstaging = scf_drv._xcfun_ldstaging

        # D4 dispersion correction
        self.dispersion = scf_drv.dispersion

    def compute(self, molecule, basis, scf_results=None):
        """
        Performs calculation of gradient.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        :param scf_results:
            The dictionary containing converged SCF results.
        """

        if scf_results is None:
            scf_results = self.scf_driver.scf_tensors

        start_time = time.time()
        self.print_header()

        if self.numerical:
            self.ostream.mute()
            self.compute_numerical(molecule, basis, scf_results)
            self.ostream.unmute()

        else:

            # sanity checks
            molecule_sanity_check(molecule)
            scf_results_sanity_check(self, self.scf_driver.scf_tensors)
            dft_sanity_check(self, 'compute')

            if self.rank == mpi_master():
                scf_type = scf_results['scf_type']
            else:
                scf_type = None
            scf_type = self.comm.bcast(scf_type, root=mpi_master())

            if scf_type == 'restricted':
                self.compute_analytical_restricted(molecule, basis, scf_results)

            elif scf_type == 'unrestricted':
                self.compute_analytical_unrestricted(molecule, basis,
                                                     scf_results)

            else:
                assert_msg_critical(
                    False,
                    'ScfGradientDriver: Not implemented for restricted open-shell'
                )

        # print gradient
        self.print_geometry(molecule)
        self.print_gradient(molecule)

        valstr = '*** Time spent in gradient calculation: '
        valstr += '{:.2f} sec ***'.format(time.time() - start_time)
        self.ostream.print_header(valstr)
        self.ostream.print_blank()
        self.ostream.flush()

    def compute_analytical_restricted(self, molecule, basis, scf_results):
        """
        Performs calculation of gradient for restricted SCF.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        :param scf_results:
            The dictionary containing converged SCF results.
        """

        grad_timing = {
            'Screening': 0.0,
            'Overlap_grad': 0.0,
            'Kinetic_energy_grad': 0.0,
            'Nuclear_potential_grad': 0.0,
            'Point_charges_grad': 0.0,
            'Fock_grad': 0.0,
            'XC_grad': 0.0,
            'PE_grad': 0.0,
            'CPCM_grad': 0.0,
            'D4_grad': 0.0,
            'Classical': 0.0,
        }

        if self.rank == mpi_master():
            D = scf_results['D_alpha']
            nocc = molecule.number_of_alpha_electrons()
            ene_occ = scf_results['E_alpha'][:nocc]
            mo_occ = scf_results['C_alpha'][:, :nocc].copy()
            W = np.linalg.multi_dot([mo_occ, np.diag(ene_occ), mo_occ.T])
        else:
            D = None
            W = None

        D = self.comm.bcast(D, root=mpi_master())
        W = self.comm.bcast(W, root=mpi_master())

        natoms = molecule.number_of_atoms()

        self.gradient = np.zeros((natoms, 3))

        local_atoms = molecule.partition_atoms(self.comm)

        # kinetic energy contribution to gradient

        t0 = time.time()

        kin_grad_drv = KineticEnergyGeom100Driver()

        for iatom in local_atoms:
            gmats = kin_grad_drv.compute(molecule, basis, iatom)

            for i, label in enumerate(['X', 'Y', 'Z']):
                gmat = gmats.matrix_to_numpy(label)
                self.gradient[iatom, i] += 2.0 * np.sum((gmat + gmat.T) * D)

            gmats = Matrices()

        grad_timing['Kinetic_energy_grad'] += time.time() - t0

        # nuclear potential contribution to gradient

        t0 = time.time()

        npot_grad_100_drv = NuclearPotentialGeom100Driver()
        npot_grad_010_drv = NuclearPotentialGeom010Driver()

        for iatom in local_atoms:
            gmats_100 = npot_grad_100_drv.compute(molecule, basis, iatom)
            gmats_010 = npot_grad_010_drv.compute(molecule, basis, iatom)

            for i, label in enumerate(['X', 'Y', 'Z']):
                gmat_100 = gmats_100.matrix_to_numpy(label)
                gmat_010 = gmats_010.matrix_to_numpy(label)

                # TODO: move minus sign into function call (such as in oneints)
                self.gradient[iatom, i] -= 2.0 * np.sum(
                    (gmat_100 + gmat_100.T) * D)
                self.gradient[iatom, i] -= 2.0 * np.sum(gmat_010 * D)

            gmats_100 = Matrices()
            gmats_010 = Matrices()

        grad_timing['Nuclear_potential_grad'] += time.time() - t0

        t0 = time.time()

        # point charges contribution
        if self.scf_driver.point_charges is not None:
            npoints = self.scf_driver.point_charges.shape[1]

            mm_coords = []
            mm_charges = []
            for p in range(npoints):
                xyz_p = self.scf_driver.point_charges[:3, p]
                chg_p = self.scf_driver.point_charges[3, p]
                mm_coords.append(xyz_p.copy())
                mm_charges.append(chg_p)

            for iatom in local_atoms:
                gmats_100 = npot_grad_100_drv.compute(molecule, basis, iatom,
                                                      mm_coords, mm_charges)

                for i, label in enumerate(['X', 'Y', 'Z']):
                    gmat_100 = gmats_100.matrix_to_numpy(label)

                    # TODO: move minus sign into function call (such as in oneints)
                    self.gradient[iatom, i] -= 2.0 * np.sum(
                        (gmat_100 + gmat_100.T) * D)

                gmats_100 = Matrices()

            grad_timing['Point_charges_grad'] += time.time() - t0

        # orbital contribution to gradient
    
        t0 = time.time()
    
        ovl_grad_drv = OverlapGeom100Driver()
    
        for iatom in local_atoms:
            gmats = ovl_grad_drv.compute(molecule, basis, iatom)
    
            for i, label in enumerate(['X', 'Y', 'Z']):
                gmat = gmats.matrix_to_numpy(label)
                # Note: minus sign for energy weighted density
                self.gradient[iatom, i] -= 2.0 * np.sum((gmat + gmat.T) * W)
    
            gmats = Matrices()
    
        grad_timing['Overlap_grad'] += time.time() - t0

        # ERI contribution to gradient

        if self.rank == mpi_master():
            use_dft = 'xcfun' in scf_results
        else:
            use_dft = False
        use_dft = self.comm.bcast(use_dft, root=mpi_master())

        # determine fock_type and exchange_scaling_factor
        if use_dft:
            if self.rank == mpi_master():
                xcfun = scf_results['xcfun']
            else:
                xcfun = None
            xcfun = self.comm.bcast(xcfun, root=mpi_master())
            xcfun = parse_xc_func(xcfun)

            if xcfun.is_hybrid():
                fock_type = '2jkx'
                exchange_scaling_factor = xcfun.get_frac_exact_exchange()
            else:
                fock_type = 'j'
                exchange_scaling_factor = 0.0
        else:
            fock_type = '2jk'
            exchange_scaling_factor = 1.0

        # further determine exchange_scaling_factor, erf_k_coef and omega
        need_omega = (use_dft and xcfun.is_range_separated())
        if need_omega:
            exchange_scaling_factor = (xcfun.get_rs_alpha() +
                                       xcfun.get_rs_beta())
            erf_k_coef = -xcfun.get_rs_beta()
            omega = xcfun.get_rs_omega()
        else:
            erf_k_coef, omega = None, None

        den_mat_for_fock = make_matrix(basis, mat_t.symmetric)
        den_mat_for_fock.set_values(D)

        den_mat_for_fock2 = make_matrix(basis, mat_t.general)
        den_mat_for_fock2.set_values(D)

        fock_grad_drv = FockGeom1000Driver()
        fock_grad_drv._set_block_size_factor(self._block_size_factor)

        t0 = time.time()

        screener = T4CScreener()
        screener.partition(basis, molecule, 'eri')

        grad_timing['Screening'] += time.time() - t0

        thresh_int = int(-math.log10(self.eri_thresh))

        if self.scf_driver.ri_coulomb:
            assert_msg_critical(
                basis.get_label().lower().startswith('def2-'),
                'ScfGradientDriver: Invalid basis set for RI-J')

            self.ostream.print_info(
                'Using the resolution of the identity (RI) approximation.')
            self.ostream.print_blank()
            self.ostream.flush()

            if self.rank == mpi_master():
                basis_ri_j = MolecularBasis.read(
                    molecule, self.scf_driver.ri_auxiliary_basis)
            else:
                basis_ri_j = None
            basis_ri_j = self.comm.bcast(basis_ri_j, root=mpi_master())

            ri_gvec = self.scf_driver._ri_drv.compute_bq_vector(
                den_mat_for_fock)

            ri_grad_drv = RIFockGradDriver()

            t0 = time.time()

            for iatom in local_atoms:
                atomgrad = ri_grad_drv.direct_compute(screener, basis,
                                                      basis_ri_j, molecule,
                                                      ri_gvec, den_mat_for_fock,
                                                      iatom, thresh_int)

                # Note: RI gradient already contains factor of 2 for
                # closed-shell
                self.gradient[iatom, :] += np.array(atomgrad.coordinates())

            grad_timing['Fock_grad'] += time.time() - t0

        else:

            for iatom in local_atoms:

                t0 = time.time()

                screener_atom = T4CScreener()
                screener_atom.partition_atom(basis, molecule, 'eri', iatom)

                grad_timing['Screening'] += time.time() - t0

                t0 = time.time()

                atomgrad = fock_grad_drv.compute(basis, screener_atom, screener,
                                                 den_mat_for_fock,
                                                 den_mat_for_fock2, iatom,
                                                 fock_type,
                                                 exchange_scaling_factor, 0.0,
                                                 thresh_int)

                factor = 2.0 if fock_type == 'j' else 1.0

                self.gradient[iatom, :] += np.array(atomgrad) * factor

                if need_omega:
                    # for range-separated functional
                    atomgrad_rs = fock_grad_drv.compute(
                        basis, screener_atom, screener, den_mat_for_fock,
                        den_mat_for_fock2, iatom, 'kx_rs', erf_k_coef, omega,
                        thresh_int)

                    self.gradient[iatom, :] -= np.array(atomgrad_rs)

                grad_timing['Fock_grad'] += time.time() - t0

        # XC contribution to gradient

        t0 = time.time()

        if use_dft:
            if self.rank == mpi_master():
                xcfun_label = scf_results['xcfun']
            else:
                xcfun_label = None
            xcfun_label = self.comm.bcast(xcfun_label, root=mpi_master())

            grad_drv = XCMolecularGradient()
            self.gradient += grad_drv.integrate_vxc_gradient(
                molecule, basis, [D], self.scf_driver._mol_grid, xcfun_label)

            grad_timing['XC_grad'] += time.time() - t0

        else:
            xcfun_label = 'hf'

        # Embedding contribution to the gradient

        t0 = time.time()

        if self.scf_driver._pe:
            from .embedding import PolarizableEmbeddingGrad

            # pass along emb object from scf, or make a new one? -> for ind dipoles.
            self._embedding_drv = PolarizableEmbeddingGrad(
                molecule=molecule,
                ao_basis=basis,
                options=self.scf_driver.embedding,
                comm=self.comm)

            pe_grad = self._embedding_drv.compute_pe_contributions(
                density_matrix=2.0 * D)

            if self.rank == mpi_master():
                self.gradient += pe_grad

            grad_timing['PE_grad'] += time.time() - t0

        # CPCM contribution to gradient

        if self.scf_driver._cpcm:
            assert_msg_critical(not self.scf_driver._smd,
                                'Cannot use SMD in gradient calculation')
            self.gradient += self.scf_driver.cpcm_drv.compute_gradient(
                molecule, basis, 2.0 * D)
            grad_timing['CPCM_grad'] += time.time() - t0

        # nuclear contribution to gradient
        # and D4 dispersion correction if requested
        # (only added on master rank)

        t0 = time.time()

        if self.rank == mpi_master():
            self.gradient += self.grad_nuc_contrib(molecule)

            grad_timing['Classical'] += time.time() - t0

            t0 = time.time()

            if self.dispersion or (
                    self.scf_driver._dft and
                    'D4' in self.scf_driver.xcfun.get_func_label().upper()):
                disp = DispersionModel()
                disp.compute(molecule, xcfun_label)
                self.gradient += disp.get_gradient()

                grad_timing['D4_grad'] += time.time() - t0

            t0 = time.time()

        # nuclei-point charges contribution to gradient

        t0 = time.time()

        if self.scf_driver.point_charges is not None:
            coords = molecule.get_coordinates_in_bohr()
            nuclear_charges = molecule.get_element_ids()
            npoints = self.scf_driver.point_charges.shape[1]

            for a in range(self.rank, natoms, self.nodes):
                z_a = nuclear_charges[a]
                r_a = coords[a]

                for p in range(npoints):
                    r_p = self.scf_driver.point_charges[:3, p]
                    q_p = self.scf_driver.point_charges[3, p]
                    r = np.linalg.norm(r_a - r_p)
                    f_ij = z_a * q_p * (r_p - r_a) / r**3

                    self.gradient[a] += f_ij

            if self.scf_driver.qm_vdw_params is not None:
                vdw_grad = np.zeros((natoms, 3))

                for a in range(self.rank, natoms, self.nodes):
                    xyz_i = coords[a]
                    sigma_i = self.scf_driver.qm_vdw_params[a, 0]
                    epsilon_i = self.scf_driver.qm_vdw_params[a, 1]

                    for p in range(npoints):
                        xyz_j = self.scf_driver.point_charges[:3, p]
                        sigma_j = self.scf_driver.point_charges[4, p]
                        epsilon_j = self.scf_driver.point_charges[5, p]

                        r_ij = xyz_j - xyz_i
                        distance_ij = np.linalg.norm(r_ij)
                        n_ij = r_ij / distance_ij

                        # bohr to nm
                        distance_ij *= bohr_in_angstrom() * 0.1

                        epsilon_ij = np.sqrt(epsilon_i * epsilon_j)
                        sigma_ij = 0.5 * (sigma_i + sigma_j)

                        sigma_r_6 = (sigma_ij / distance_ij)**6
                        sigma_r_12 = sigma_r_6**2

                        g = -24.0 * epsilon_ij * (2.0 * sigma_r_12 / distance_ij
                                                  - sigma_r_6 / distance_ij)

                        vdw_grad[a] += -g * n_ij

                # convert gradient to atomic unit
                vdw_grad /= (hartree_in_kjpermol() * 10.0 / bohr_in_angstrom())

                self.gradient += vdw_grad

        grad_timing['Classical'] += time.time() - t0

        # collect gradient

        self.gradient = self.comm.allreduce(self.gradient, op=MPI.SUM)

        if self.timing and self.rank == mpi_master():
            self.ostream.print_info('Gradient timing decomposition')
            for key, val in grad_timing.items():
                if val > 0.0:
                    self.ostream.print_info(f'    {key:<25}:  {val:.2f} sec')
            self.ostream.print_blank()

    def compute_analytical_unrestricted(self, molecule, basis, scf_results):
        """
        Performs calculation of gradient for unrestricted SCF.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        :param scf_results:
            The dictionary containing converged SCF results.
        """

        if self.rank == mpi_master():
            Da = scf_results['D_alpha']
            Db = scf_results['D_beta']

            nocc_a = molecule.number_of_alpha_electrons()
            nocc_b = molecule.number_of_beta_electrons()

            ene_occ_a = scf_results['E_alpha'][:nocc_a]
            ene_occ_b = scf_results['E_beta'][:nocc_b]

            mo_occ_a = scf_results['C_alpha'][:, :nocc_a].copy()
            mo_occ_b = scf_results['C_beta'][:, :nocc_b].copy()

            Wa = np.linalg.multi_dot([mo_occ_a, np.diag(ene_occ_a), mo_occ_a.T])
            Wb = np.linalg.multi_dot([mo_occ_b, np.diag(ene_occ_b), mo_occ_b.T])
        else:
            Da = None
            Db = None
            Wa = None
            Wb = None

        Da = self.comm.bcast(Da, root=mpi_master())
        Db = self.comm.bcast(Db, root=mpi_master())
        Wa = self.comm.bcast(Wa, root=mpi_master())
        Wb = self.comm.bcast(Wb, root=mpi_master())

        natoms = molecule.number_of_atoms()

        self.gradient = np.zeros((natoms, 3))

        local_atoms = molecule.partition_atoms(self.comm)

        # kinetic energy contribution to gradient

        kin_grad_drv = KineticEnergyGeom100Driver()

        for iatom in local_atoms:
            gmats = kin_grad_drv.compute(molecule, basis, iatom)

            for i, label in enumerate(['X', 'Y', 'Z']):
                gmat = gmats.matrix_to_numpy(label)
                self.gradient[iatom, i] += np.sum((gmat + gmat.T) * (Da + Db))

            gmats = Matrices()

        # nuclear potential contribution to gradient

        npot_grad_100_drv = NuclearPotentialGeom100Driver()
        npot_grad_010_drv = NuclearPotentialGeom010Driver()

        for iatom in local_atoms:
            gmats_100 = npot_grad_100_drv.compute(molecule, basis, iatom)
            gmats_010 = npot_grad_010_drv.compute(molecule, basis, iatom)

            for i, label in enumerate(['X', 'Y', 'Z']):
                gmat_100 = gmats_100.matrix_to_numpy(label)
                gmat_010 = gmats_010.matrix_to_numpy(label)

                # TODO: move minus sign into function call (such as in oneints)
                self.gradient[iatom, i] -= np.sum(
                    (gmat_100 + gmat_100.T) * (Da + Db))
                self.gradient[iatom, i] -= np.sum(gmat_010 * (Da + Db))

            gmats_100 = Matrices()
            gmats_010 = Matrices()

        # point charges contribution
        if self.scf_driver.point_charges is not None:
            npoints = self.scf_driver.point_charges.shape[1]

            mm_coords = []
            mm_charges = []
            for p in range(npoints):
                xyz_p = self.scf_driver.point_charges[:3, p]
                chg_p = self.scf_driver.point_charges[3, p]
                mm_coords.append(xyz_p.copy())
                mm_charges.append(chg_p)

            for iatom in local_atoms:
                gmats_100 = npot_grad_100_drv.compute(molecule, basis, iatom,
                                                      mm_coords, mm_charges)

                for i, label in enumerate(['X', 'Y', 'Z']):
                    gmat_100 = gmats_100.matrix_to_numpy(label)

                    # TODO: move minus sign into function call (such as in oneints)
                    self.gradient[iatom, i] -= np.sum(
                        (gmat_100 + gmat_100.T) * (Da + Db))

                gmats_100 = Matrices()

        # orbital contribution to gradient

        ovl_grad_drv = OverlapGeom100Driver()

        for iatom in local_atoms:
            gmats = ovl_grad_drv.compute(molecule, basis, iatom)

            for i, label in enumerate(['X', 'Y', 'Z']):
                gmat = gmats.matrix_to_numpy(label)
                # Note: minus sign for energy weighted density
                self.gradient[iatom, i] -= np.sum((gmat + gmat.T) * (Wa + Wb))

            gmats = Matrices()

        # ERI contribution to gradient

        if self.rank == mpi_master():
            use_dft = 'xcfun' in scf_results
        else:
            use_dft = False
        use_dft = self.comm.bcast(use_dft, root=mpi_master())

        # determine fock_type and exchange_scaling_factor
        if use_dft:
            if self.rank == mpi_master():
                xcfun = scf_results['xcfun']
            else:
                xcfun = None
            xcfun = self.comm.bcast(xcfun, root=mpi_master())
            xcfun = parse_xc_func(xcfun)

            if xcfun.is_hybrid():
                fock_type = '2jkx'
                exchange_scaling_factor = xcfun.get_frac_exact_exchange()
            else:
                fock_type = 'j'
                exchange_scaling_factor = 0.0
        else:
            fock_type = '2jk'
            exchange_scaling_factor = 1.0

        # further determine exchange_scaling_factor, erf_k_coef and omega
        need_omega = (use_dft and xcfun.is_range_separated())
        if need_omega:
            exchange_scaling_factor = (xcfun.get_rs_alpha() +
                                       xcfun.get_rs_beta())
            erf_k_coef = -xcfun.get_rs_beta()
            omega = xcfun.get_rs_omega()
        else:
            erf_k_coef, omega = None, None

        Da_for_fock = make_matrix(basis, mat_t.symmetric)
        Da_for_fock.set_values(Da)

        Db_for_fock = make_matrix(basis, mat_t.symmetric)
        Db_for_fock.set_values(Db)

        Dab_for_fock = make_matrix(basis, mat_t.symmetric)
        Dab_for_fock.set_values(Da + Db)

        Da_for_fock_2 = make_matrix(basis, mat_t.general)
        Da_for_fock_2.set_values(Da)

        Db_for_fock_2 = make_matrix(basis, mat_t.general)
        Db_for_fock_2.set_values(Db)

        Dab_for_fock_2 = make_matrix(basis, mat_t.general)
        Dab_for_fock_2.set_values(Da + Db)

        fock_grad_drv = FockGeom1000Driver()
        fock_grad_drv._set_block_size_factor(self._block_size_factor)

        screener = T4CScreener()
        screener.partition(basis, molecule, 'eri')

        thresh_int = int(-math.log10(self.eri_thresh))

        if self.scf_driver.ri_coulomb:
            assert_msg_critical(
                basis.get_label().lower().startswith('def2-'),
                'ScfGradientDriver: Invalid basis set for RI-J')

            self.ostream.print_info(
                'Using the resolution of the identity (RI) approximation.')
            self.ostream.print_blank()
            self.ostream.flush()

            if self.rank == mpi_master():
                basis_ri_j = MolecularBasis.read(
                    molecule, self.scf_driver.ri_auxiliary_basis)
            else:
                basis_ri_j = None
            basis_ri_j = self.comm.bcast(basis_ri_j, root=mpi_master())

            ri_gvec = self.scf_driver._ri_drv.compute_bq_vector(Dab_for_fock)

            ri_grad_drv = RIFockGradDriver()

            for iatom in local_atoms:
                atomgrad = ri_grad_drv.direct_compute(screener, basis,
                                                      basis_ri_j, molecule,
                                                      ri_gvec, Dab_for_fock,
                                                      iatom, thresh_int)
                # Note: RI gradient already contains factor of 2 for
                # closed-shell, so for unrestricted factor becomes 0.25
                self.gradient[iatom, :] += 0.25 * np.array(
                    atomgrad.coordinates())

        else:

            for iatom in local_atoms:

                screener_atom = T4CScreener()
                screener_atom.partition_atom(basis, molecule, 'eri', iatom)

                atomgrad_Jab = fock_grad_drv.compute(basis, screener_atom,
                                                     screener, Dab_for_fock,
                                                     Dab_for_fock_2, iatom, 'j',
                                                     0.0, 0.0, thresh_int)

                self.gradient[iatom, :] += 0.5 * np.array(atomgrad_Jab)

                if fock_type != 'j':
                    atomgrad_Ka = fock_grad_drv.compute(
                        basis, screener_atom, screener, Da_for_fock,
                        Da_for_fock_2, iatom, 'kx', exchange_scaling_factor,
                        0.0, thresh_int)
                    atomgrad_Kb = fock_grad_drv.compute(
                        basis, screener_atom, screener, Db_for_fock,
                        Db_for_fock_2, iatom, 'kx', exchange_scaling_factor,
                        0.0, thresh_int)

                    self.gradient[iatom, :] -= 0.5 * np.array(atomgrad_Ka)
                    self.gradient[iatom, :] -= 0.5 * np.array(atomgrad_Kb)

                    if need_omega:
                        # for range-separated functional
                        atomgrad_Ka_rs = fock_grad_drv.compute(
                            basis, screener_atom, screener, Da_for_fock,
                            Da_for_fock_2, iatom, 'kx_rs', erf_k_coef, omega,
                            thresh_int)
                        atomgrad_Kb_rs = fock_grad_drv.compute(
                            basis, screener_atom, screener, Db_for_fock,
                            Db_for_fock_2, iatom, 'kx_rs', erf_k_coef, omega,
                            thresh_int)

                        self.gradient[
                            iatom, :] -= 0.5 * np.array(atomgrad_Ka_rs)
                        self.gradient[
                            iatom, :] -= 0.5 * np.array(atomgrad_Kb_rs)

        # XC contribution to gradient

        if use_dft:
            if self.rank == mpi_master():
                xcfun_label = scf_results['xcfun']
            else:
                xcfun_label = None
            xcfun_label = self.comm.bcast(xcfun_label, root=mpi_master())

            grad_drv = XCMolecularGradient()
            self.gradient += grad_drv.integrate_vxc_gradient(
                molecule, basis, [Da, Db], self.scf_driver._mol_grid,
                xcfun_label)

        else:
            xcfun_label = 'hf'

        # Embedding contribution to the gradient
        if self.scf_driver._pe:
            from .embedding import PolarizableEmbeddingGrad

            # pass along emb object from scf, or make a new one? -> for ind dipoles.
            self._embedding_drv = PolarizableEmbeddingGrad(
                molecule=molecule,
                ao_basis=basis,
                options=self.scf_driver.embedding,
                comm=self.comm)

            pe_grad = self._embedding_drv.compute_pe_contributions(
                density_matrix=Da + Db)

            if self.rank == mpi_master():
                self.gradient += pe_grad

        # CPCM contribution to gradient

        if self.scf_driver._cpcm:
            assert_msg_critical(not self.scf_driver._smd,
                                'Cannot use SMD in gradient calculation')
            self.gradient += self.scf_driver.cpcm_drv.compute_gradient(
                molecule, basis, Da + Db)

        # nuclear contribution to gradient
        # and D4 dispersion correction if requested
        # (only added on master rank)

        if self.rank == mpi_master():
            self.gradient += self.grad_nuc_contrib(molecule)

            if self.dispersion or (
                    self.scf_driver._dft and
                    'D4' in self.scf_driver.xcfun.get_func_label().upper()):
                disp = DispersionModel()
                disp.compute(molecule, xcfun_label)
                self.gradient += disp.get_gradient()

        # nuclei-point charges contribution to gradient

        if self.scf_driver.point_charges is not None:
            coords = molecule.get_coordinates_in_bohr()
            nuclear_charges = molecule.get_element_ids()
            npoints = self.scf_driver.point_charges.shape[1]

            for a in range(self.rank, natoms, self.nodes):
                z_a = nuclear_charges[a]
                r_a = coords[a]

                for p in range(npoints):
                    r_p = self.scf_driver.point_charges[:3, p]
                    q_p = self.scf_driver.point_charges[3, p]
                    r = np.linalg.norm(r_a - r_p)
                    f_ij = z_a * q_p * (r_p - r_a) / r**3

                    self.gradient[a] += f_ij

            if self.scf_driver.qm_vdw_params is not None:
                vdw_grad = np.zeros((natoms, 3))

                for a in range(self.rank, natoms, self.nodes):
                    xyz_i = coords[a]
                    sigma_i = self.scf_driver.qm_vdw_params[a, 0]
                    epsilon_i = self.scf_driver.qm_vdw_params[a, 1]

                    for p in range(npoints):
                        xyz_j = self.scf_driver.point_charges[:3, p]
                        sigma_j = self.scf_driver.point_charges[4, p]
                        epsilon_j = self.scf_driver.point_charges[5, p]

                        r_ij = xyz_j - xyz_i
                        distance_ij = np.linalg.norm(r_ij)
                        n_ij = r_ij / distance_ij

                        # bohr to nm
                        distance_ij *= bohr_in_angstrom() * 0.1

                        epsilon_ij = np.sqrt(epsilon_i * epsilon_j)
                        sigma_ij = 0.5 * (sigma_i + sigma_j)

                        sigma_r_6 = (sigma_ij / distance_ij)**6
                        sigma_r_12 = sigma_r_6**2

                        g = -24.0 * epsilon_ij * (2.0 * sigma_r_12 / distance_ij
                                                  - sigma_r_6 / distance_ij)

                        vdw_grad[a] += -g * n_ij

                # convert gradient to atomic unit
                vdw_grad /= (hartree_in_kjpermol() * 10.0 / bohr_in_angstrom())

                self.gradient += vdw_grad

        # collect gradient

        self.gradient = self.comm.allreduce(self.gradient, op=MPI.SUM)

    def compute_energy(self, molecule, ao_basis, scf_results=None):
        """
        Computes the energy at current geometry.

        :param molecule:
            The molecule.
        :param ao_basis:
            The AO basis set.
        :param scf_results:
            The dictionary containing converged SCF results.

        :return:
            The energy.
        """

        # if not self._debug:
        #     self.ostream.mute()

        if self.numerical:
            # disable restarting scf for numerical gradient
            self.scf_driver.restart = False
        else:
            # always try restarting scf for analytical gradient
            self.scf_driver.restart = True
        new_scf_results = self.scf_driver.compute(molecule, ao_basis)
        assert_msg_critical(self.scf_driver.is_converged,
                            'ScfGradientDriver: SCF did not converge')

        # if not self._debug:
        #     self.ostream.unmute()

        if (self.rank == mpi_master()) and (scf_results is not None):
            scf_results.update(new_scf_results)

        return self.scf_driver.get_scf_energy()

    def __deepcopy__(self, memo):
        """
        Implements deepcopy.

        :param memo:
            The memo dictionary for deepcopy.

        :return:
            A deepcopy of self.
        """

        new_grad_drv = ScfGradientDriver(self.comm, self.ostream)

        for key, val in vars(self).items():
            if isinstance(val, (MPI.Intracomm, OutputStream)):
                pass
            elif isinstance(val, XCFunctional):
                new_grad_drv.key = XCFunctional(val)
            elif isinstance(val, MolecularGrid):
                new_grad_drv.key = MolecularGrid(val)
            else:
                new_grad_drv.key = deepcopy(val)

        return new_grad_drv
