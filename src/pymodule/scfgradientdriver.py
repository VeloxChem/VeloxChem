#
#                              VELOXCHEM
#         ----------------------------------------------------
#                     An Electronic Structure Code
#
#  Copyright © 2018-2024 by VeloxChem developers. All rights reserved.
#
#  SPDX-License-Identifier: LGPL-3.0-or-later
#
#  This file is part of VeloxChem.
#
#  VeloxChem is free software: you can redistribute it and/or modify it under
#  the terms of the GNU Lesser General Public License as published by the Free
#  Software Foundation, either version 3 of the License, or (at your option)
#  any later version.
#
#  VeloxChem is distributed in the hope that it will be useful, but WITHOUT
#  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
#  License for more details.
#
#  You should have received a copy of the GNU Lesser General Public License
#  along with VeloxChem. If not, see <https://www.gnu.org/licenses/>.

from mpi4py import MPI
from copy import deepcopy
import numpy as np
import time
import math

from .veloxchemlib import (OverlapGeom100Driver, KineticEnergyGeom100Driver,
                           NuclearPotentialGeom100Driver,
                           NuclearPotentialGeom010Driver, FockGeom1000Driver)
from .veloxchemlib import XCFunctional, MolecularGrid, XCMolecularGradient
from .veloxchemlib import DispersionModel
from .veloxchemlib import T4CScreener
from .veloxchemlib import mpi_master, mat_t
from .veloxchemlib import make_matrix
from .veloxchemlib import parse_xc_func
from .veloxchemlib import bohr_in_angstrom, hartree_in_kcalpermol
from .matrices import Matrices
from .profiler import Profiler
from .outputstream import OutputStream
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

    def _print_debug_info(self, label):
        """
        Prints debug information.

        :param label:
            The label of debug information.
        """

        if self._debug:
            profiler = Profiler()
            self.ostream.print_info(f'==DEBUG==   available memory {label}: ' +
                                    profiler.get_available_memory())
            self.ostream.flush()

    def partition_atoms(self, molecule):
        """
        Partition atoms for parallel computation of gradient.

        :param molecule:
            The molecule.

        :return:
            The list of atom indices for the current MPI rank.
        """

        if self.rank == mpi_master():
            elem_ids = molecule.get_identifiers()
            coords = molecule.get_coordinates_in_bohr()
            mol_com = molecule.center_of_mass_in_bohr()

            r2_array = np.sum((coords - mol_com)**2, axis=1)
            sorted_r2_list = sorted([
                (r2, nchg, i)
                for i, (r2, nchg) in enumerate(zip(r2_array, elem_ids))
            ])

            dict_atoms = {}
            for r2, nchg, i in sorted_r2_list:
                if nchg not in dict_atoms:
                    dict_atoms[nchg] = []
                dict_atoms[nchg].append(i)

            list_atoms = []
            for nchg in sorted(dict_atoms.keys(), reverse=True):
                list_atoms += dict_atoms[nchg]

        else:
            list_atoms = None

        list_atoms = self.comm.bcast(list_atoms, root=mpi_master())

        return list_atoms[self.rank::self.nodes]

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
            'classical': 0.0,
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

        local_atoms = self.partition_atoms(molecule)

        # kinetic energy contribution to gradient

        t0 = time.time()

        kin_grad_drv = KineticEnergyGeom100Driver()

        self._print_debug_info('before kin_grad')

        for iatom in local_atoms:
            gmats = kin_grad_drv.compute(molecule, basis, iatom)

            for i, label in enumerate(['X', 'Y', 'Z']):
                gmat = gmats.matrix_to_numpy(label)
                self.gradient[iatom, i] += 2.0 * np.sum((gmat + gmat.T) * D)

            gmats = Matrices()

        self._print_debug_info('after  kin_grad')

        grad_timing['Kinetic_energy_grad'] += time.time() - t0

        # nuclear potential contribution to gradient

        t0 = time.time()

        self._print_debug_info('before npot_grad')

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

        self._print_debug_info('after  npot_grad')

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

        self._print_debug_info('before ovl_grad')

        ovl_grad_drv = OverlapGeom100Driver()

        for iatom in local_atoms:
            gmats = ovl_grad_drv.compute(molecule, basis, iatom)

            for i, label in enumerate(['X', 'Y', 'Z']):
                gmat = gmats.matrix_to_numpy(label)
                # Note: minus sign for energy weighted density
                self.gradient[iatom, i] -= 2.0 * np.sum((gmat + gmat.T) * W)

            gmats = Matrices()

        self._print_debug_info('after  ovl_grad')

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

        # TODO: range-separated Fock
        need_omega = (use_dft and xcfun.is_range_separated())
        if need_omega:
            assert_msg_critical(
                False, 'ScfGradientDriver: Not implemented for' +
                ' range-separated functional')

        den_mat_for_fock = make_matrix(basis, mat_t.symmetric)
        den_mat_for_fock.set_values(D)

        den_mat_for_fock2 = make_matrix(basis, mat_t.general)
        den_mat_for_fock2.set_values(D)

        self._print_debug_info('before fock_grad')

        fock_grad_drv = FockGeom1000Driver()
        fock_grad_drv._set_block_size_factor(self._block_size_factor)

        t0 = time.time()

        screener = T4CScreener()
        screener.partition(basis, molecule, 'eri')

        grad_timing['Screening'] += time.time() - t0

        thresh_int = int(-math.log10(self.eri_thresh))

        for iatom in local_atoms:

            t0 = time.time()

            screener_atom = T4CScreener()
            screener_atom.partition_atom(basis, molecule, 'eri', iatom)

            grad_timing['Screening'] += time.time() - t0

            t0 = time.time()

            atomgrad = fock_grad_drv.compute(basis, screener_atom, screener,
                                             den_mat_for_fock,
                                             den_mat_for_fock2, iatom,
                                             fock_type, exchange_scaling_factor,
                                             0.0, thresh_int)

            grad_timing['Fock_grad'] += time.time() - t0

            factor = 2.0 if fock_type == 'j' else 1.0

            self.gradient[iatom, :] += np.array(atomgrad) * factor

        self._print_debug_info('after  fock_grad')

        # XC contribution to gradient

        t0 = time.time()

        self._print_debug_info('before xc_grad')

        if use_dft:
            if self.rank == mpi_master():
                xcfun_label = scf_results['xcfun']
            else:
                xcfun_label = None
            xcfun_label = self.comm.bcast(xcfun_label, root=mpi_master())

            grad_drv = XCMolecularGradient()
            self.gradient += grad_drv.integrate_vxc_gradient(
                molecule, basis, [D], self.scf_driver._mol_grid, xcfun_label)

        else:
            xcfun_label = 'hf'

        self._print_debug_info('after  xc_grad')

        grad_timing['XC_grad'] += time.time() - t0

        # nuclear contribution to gradient
        # and D4 dispersion correction if requested
        # (only added on master rank)

        t0 = time.time()

        if self.rank == mpi_master():
            self.gradient += self.grad_nuc_contrib(molecule)

            if self.dispersion:
                disp = DispersionModel()
                disp.compute(molecule, xcfun_label)
                self.gradient += disp.get_gradient().to_numpy()

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
                vdw_grad /= (4.184 * hartree_in_kcalpermol() * 10.0 /
                             bohr_in_angstrom())

                self.gradient += vdw_grad

        grad_timing['classical'] += time.time() - t0

        # collect gradient

        self.gradient = self.comm.allreduce(self.gradient, op=MPI.SUM)

        if self.timing and self.rank == mpi_master():
            self.ostream.print_info('Gradient timing decomposition')
            for key, val in grad_timing.items():
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

        local_atoms = self.partition_atoms(molecule)

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

        # TODO: range-separated Fock
        need_omega = (use_dft and xcfun.is_range_separated())
        if need_omega:
            assert_msg_critical(
                False, 'ScfGradientDriver: Not implemented for' +
                ' range-separated functional')

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

        for iatom in local_atoms:

            screener_atom = T4CScreener()
            screener_atom.partition_atom(basis, molecule, 'eri', iatom)

            atomgrad_Jab = fock_grad_drv.compute(basis, screener_atom, screener,
                                                 Dab_for_fock, Dab_for_fock_2,
                                                 iatom, 'j', 0.0, 0.0,
                                                 thresh_int)

            self.gradient[iatom, :] += 0.5 * np.array(atomgrad_Jab)

            if fock_type != 'j':
                atomgrad_Ka = fock_grad_drv.compute(basis, screener_atom,
                                                    screener, Da_for_fock,
                                                    Da_for_fock_2, iatom, 'kx',
                                                    exchange_scaling_factor,
                                                    0.0, thresh_int)
                atomgrad_Kb = fock_grad_drv.compute(basis, screener_atom,
                                                    screener, Db_for_fock,
                                                    Db_for_fock_2, iatom, 'kx',
                                                    exchange_scaling_factor,
                                                    0.0, thresh_int)

                self.gradient[iatom, :] -= 0.5 * np.array(atomgrad_Ka)
                self.gradient[iatom, :] -= 0.5 * np.array(atomgrad_Kb)

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

        # nuclear contribution to gradient
        # and D4 dispersion correction if requested
        # (only added on master rank)

        if self.rank == mpi_master():
            self.gradient += self.grad_nuc_contrib(molecule)

            if self.dispersion:
                disp = DispersionModel()
                disp.compute(molecule, xcfun_label)
                self.gradient += disp.get_gradient().to_numpy()

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
                vdw_grad /= (4.184 * hartree_in_kcalpermol() * 10.0 /
                             bohr_in_angstrom())

                self.gradient += vdw_grad

        # collect gradient

        self.gradient = self.comm.allreduce(self.gradient, op=MPI.SUM)

    def compute_energy(self, molecule, ao_basis, scf_results):
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

        if self.rank == mpi_master():
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
