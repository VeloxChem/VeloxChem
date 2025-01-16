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

import numpy as np
import time as tm
import math

from .veloxchemlib import OverlapGeom200Driver
from .veloxchemlib import OverlapGeom101Driver
from .veloxchemlib import OverlapGeom100Driver
from .veloxchemlib import KineticEnergyGeom200Driver
from .veloxchemlib import KineticEnergyGeom101Driver
from .veloxchemlib import NuclearPotentialGeom200Driver
from .veloxchemlib import NuclearPotentialGeom020Driver
from .veloxchemlib import NuclearPotentialGeom110Driver
from .veloxchemlib import NuclearPotentialGeom101Driver
from .veloxchemlib import ElectricDipoleMomentGeom100Driver
from .veloxchemlib import FockGeom2000Driver
from .veloxchemlib import FockGeom1100Driver
from .veloxchemlib import FockGeom1010Driver
from .veloxchemlib import XCMolecularHessian
from .veloxchemlib import T4CScreener
from .veloxchemlib import make_matrix, mat_t
from .veloxchemlib import mpi_master
from .molecule import Molecule
from .molecularbasis import MolecularBasis
from .griddriver import GridDriver
from .scfgradientdriver import ScfGradientDriver
from .hessiandriver import HessianDriver
from .firstorderprop import FirstOrderProperties
from .hessianorbitalresponse import HessianOrbitalResponse
from .profiler import Profiler
from .matrices import Matrices
from .dftutils import get_default_grid_level
from .errorhandler import assert_msg_critical
from .oneeints import compute_electric_dipole_integrals
from .sanitychecks import (molecule_sanity_check, scf_results_sanity_check,
                           dft_sanity_check)


class ScfHessianDriver(HessianDriver):
    """
    Implements SCF Hessian driver.

    :param scf_drv:
        The SCF driver.

    Instance variables
        - hessian: The Hessian in Hartree per Bohr**2.
        - flag: The type of Hessian driver.
        - perturbed_density: The perturbed density
    """

    def __init__(self, scf_drv):
        """
        Initializes SCF Hessian driver.
        """

        super().__init__(scf_drv.comm, scf_drv.ostream)

        self.flag = 'SCF Hessian Driver'
        self.scf_driver = scf_drv

        self.perturbed_density = None

        # TODO TEMPORARY FLAG
        # Only run orbital response for performance testing
        self.orbrsp_only = False

        # flag for printing the Hessian
        self.do_print_hessian = False

        # TODO: determine _block_size_factor for SCF Hessian driver
        # self._block_size_factor = 4

        self._xcfun_ldstaging = scf_drv._xcfun_ldstaging

        self.use_subcomms = False

        self._input_keywords['hessian'].update({
            'orbrsp_only':
                ('bool', 'whether to only run CPHF orbital response'),
            'use_subcomms':
                ('bool', 'whether to use subcommunicators in orbital response'),
        })

    def update_settings(self, method_dict, hess_dict=None, cphf_dict=None):
        """
        Updates settings in ScfHessianDriver.

        :param method_dict:
            The input dictionary of method settings group.
        :param hess_dict:
            The input dictionary of Hessian settings group.
        :param cphf_dict:
            The input dictionary of CPHF (orbital response) settings.
        :param rsp_dict:
            The input dictionary for linear response settings
            (needed to compute the polarizability gradient).
        """

        super().update_settings(method_dict, hess_dict)

        if cphf_dict is None:
            cphf_dict = {}

        self.cphf_dict = dict(cphf_dict)

    def compute(self, molecule, ao_basis):
        """
        Computes the analytical or numerical nuclear Hessian.

        :param molecule:
            The molecule.
        :param ao_basis:
            The AO basis set.
        """

        if self.rank == mpi_master():
            self.print_header()

        start_time = tm.time()

        profiler = Profiler({
            'timing': self.timing,
            'profiling': self.profiling,
            'memory_profiling': self.memory_profiling,
            'memory_tracing': self.memory_tracing,
        })

        # Save the electronic energy
        self.elec_energy = self.scf_driver.get_scf_energy()

        # TODO TEMPORARY
        if self.orbrsp_only:
            self.ostream.print_header(
                '*** WARNING only computing Hessian orbital response!')
            self.compute_orbital_response(molecule, ao_basis)
            self.ostream.print_header(
                '*** Hessian orbital response only: DONE  ***')
            self.ostream.flush()
            return

        if self.numerical:
            self.compute_numerical(molecule, ao_basis)
        else:
            self.compute_analytical(molecule, ao_basis, profiler)

        if self.rank == mpi_master():
            # print Hessian
            if self.do_print_hessian:
                self.print_geometry(molecule)
                self.ostream.print_blank()
                self.print_hessian(molecule)

            valstr = '*** Time spent in Hessian calculation: '
            valstr += '{:.2f} sec ***'.format(tm.time() - start_time)
            self.ostream.print_header(valstr)
            self.ostream.print_blank()
            self.ostream.print_blank()
            self.ostream.flush()

    def compute_numerical(self, molecule, ao_basis):
        """
        Performs the calculation of a numerical Hessian based only
        on the energy.

        :param molecule:
            The molecule.
        :param ao_basis:
            The AO basis set.
        """

        self.ostream.mute()

        # atom labels and atom basis labels
        labels = molecule.get_labels()
        atom_basis_labels = molecule.get_atom_basis_labels()

        # main basis label
        basis_label = ao_basis.get_label()

        # number of atoms
        natm = molecule.number_of_atoms()

        # atom coordinates (nx3)
        coords = molecule.get_coordinates_in_bohr()

        # charge and spin multiplicity
        charge = molecule.get_charge()
        multiplicity = molecule.get_multiplicity()

        # Hessian
        hessian = np.zeros((natm, 3, natm, 3))

        # First-order properties for gradient of dipole moment
        prop = FirstOrderProperties(self.comm, self.ostream)
        # numerical gradient (3 dipole components, no. atoms x 3 atom coords)
        if self.rank == mpi_master():
            self.dipole_gradient = np.zeros((3, 3 * natm))

        grad_drv = ScfGradientDriver(self.scf_driver)

        for i in range(natm):

            self.ostream.unmute()
            self.ostream.print_info(f'Processing atom {i + 1}/{natm}...')
            self.ostream.flush()
            self.ostream.mute()

            for x in range(3):

                coords[i, x] += self.delta_h
                new_mol = Molecule(labels, coords, 'au', atom_basis_labels)
                new_mol.set_charge(charge)
                new_mol.set_multiplicity(multiplicity)
                new_bas = MolecularBasis.read(new_mol, basis_label)
                scf_results = self.scf_driver.compute(new_mol, new_bas)
                grad_drv.compute(new_mol, new_bas, scf_results)
                grad_plus = grad_drv.get_gradient().copy()

                prop.compute_scf_prop(new_mol, new_bas, scf_results)
                if self.rank == mpi_master():
                    mu_plus = prop.get_property('dipole moment')

                coords[i, x] -= 2.0 * self.delta_h
                new_mol = Molecule(labels, coords, 'au', atom_basis_labels)
                new_mol.set_charge(charge)
                new_mol.set_multiplicity(multiplicity)
                new_bas = MolecularBasis.read(new_mol, basis_label)
                scf_results = self.scf_driver.compute(new_mol, new_bas)
                grad_drv.compute(new_mol, new_bas, scf_results)
                grad_minus = grad_drv.get_gradient().copy()

                prop.compute_scf_prop(new_mol, new_bas, scf_results)
                if self.rank == mpi_master():
                    mu_minus = prop.get_property('dipole moment')

                coords[i, x] += self.delta_h
                hessian[i, x, :, :] = ((grad_plus - grad_minus) /
                                       (2.0 * self.delta_h))

                if self.rank == mpi_master():
                    for c in range(3):
                        self.dipole_gradient[c, 3 * i +
                                             x] = ((mu_plus[c] - mu_minus[c]) /
                                                   (2.0 * self.delta_h))

        # reshaped Hessian as member variable
        self.hessian = hessian.reshape(3 * natm, 3 * natm)

        # restore scf_drv to initial state
        scf_results = self.scf_driver.compute(molecule, ao_basis)
        assert_msg_critical(self.scf_driver.is_converged,
                            'ScfHessianDriver: SCF did not converge')
        self.ostream.unmute()

    def compute_analytical(self, molecule, ao_basis, profiler):
        """
        Computes the analytical nuclear Hessian.

        :param molecule:
            The molecule.
        :param ao_basis:
            The AO basis set.
        :param profiler:
            The profiler.
        """

        # sanity checks
        molecule_sanity_check(molecule)
        scf_results_sanity_check(self, self.scf_driver.scf_tensors)
        dft_sanity_check(self, 'compute', 'hessian')

        self.ostream.print_info('Computing analytical Hessian...')
        self.ostream.print_blank()
        hess_ref = 'P. Deglmann, F. Furche, R. Ahlrichs,'
        hess_ref += ' Chem. Phys. Lett. 2002, 362, 511-518.'
        self.ostream.print_reference('Reference: ' + hess_ref)
        self.ostream.print_blank()
        self.ostream.flush()

        # Preparation

        natm = molecule.number_of_atoms()
        scf_tensors = self.scf_driver.scf_tensors

        if self.rank == mpi_master():
            density = scf_tensors['D_alpha']
            mo = scf_tensors['C_alpha']
            nocc = molecule.number_of_alpha_electrons()
            mo_occ = mo[:, :nocc]
            mo_energies = scf_tensors['E_alpha']
            eocc = mo_energies[:nocc]
            omega_ao = -np.linalg.multi_dot([mo_occ, np.diag(eocc), mo_occ.T])

        else:
            density = None
            omega_ao = None

        density = self.comm.bcast(density, root=mpi_master())
        omega_ao = self.comm.bcast(omega_ao, root=mpi_master())

        # CPHF equations

        cphf_solver = HessianOrbitalResponse(self.comm, self.ostream)

        # TODO: double check propagation of cphf settings
        profiler_keywords = {
            'timing', 'profiling', 'memory_profiling', 'memory_tracing',
            'use_subcomms'
        }
        for key in profiler_keywords:
            setattr(cphf_solver, key, getattr(self, key))

        cphf_solver.compute(molecule, ao_basis, scf_tensors)

        cphf_solution_dict = cphf_solver.cphf_results
        dist_cphf_ov = cphf_solution_dict['dist_cphf_ov']
        dist_cphf_rhs = cphf_solution_dict['dist_cphf_rhs']

        hessian_first_integral_derivatives = cphf_solution_dict[
            'hessian_first_integral_derivatives']
        hessian_eri_overlap = cphf_solution_dict['hessian_eri_overlap']

        # First-order contributions

        t1 = tm.time()

        # RHS contracted with CPHF coefficients (ov)
        hessian_cphf_coeff_rhs = np.zeros((natm, 3, natm, 3))

        for i in range(natm):
            for x in range(3):
                dist_cphf_ov_ix_data = dist_cphf_ov[i * 3 + x].data

                for j in range(i, natm):
                    for y in range(3):
                        hess_ijxy = 4.0 * (np.dot(
                            dist_cphf_ov_ix_data,
                            dist_cphf_rhs[j * 3 + y].data))

                        hessian_cphf_coeff_rhs[i, x, j, y] += hess_ijxy
                        if i != j:
                            hessian_cphf_coeff_rhs[j, y, i, x] += hess_ijxy

        hessian_cphf_coeff_rhs = self.comm.reduce(hessian_cphf_coeff_rhs,
                                                  root=mpi_master())

        if self.rank == mpi_master():
            hessian_first_order_derivatives = (
                hessian_cphf_coeff_rhs + hessian_first_integral_derivatives +
                hessian_eri_overlap)
        else:
            hessian_first_order_derivatives = None

        self.ostream.print_info('First order derivative contributions' +
                                ' to the Hessian computed in' +
                                ' {:.2f} sec.'.format(tm.time() - t1))
        self.ostream.print_blank()
        self.ostream.flush()

        # Second-order contributions

        t2 = tm.time()

        ovlp_hess_200_drv = OverlapGeom200Driver()
        ovlp_hess_101_drv = OverlapGeom101Driver()

        kin_hess_200_drv = KineticEnergyGeom200Driver()
        kin_hess_101_drv = KineticEnergyGeom101Driver()

        npot_hess_200_drv = NuclearPotentialGeom200Driver()
        npot_hess_020_drv = NuclearPotentialGeom020Driver()
        npot_hess_110_drv = NuclearPotentialGeom110Driver()
        npot_hess_101_drv = NuclearPotentialGeom101Driver()

        fock_hess_2000_drv = FockGeom2000Driver()
        fock_hess_1100_drv = FockGeom1100Driver()
        fock_hess_1010_drv = FockGeom1010Driver()

        # determine fock_type and exchange_scaling_factor
        if self._dft:
            xcfun = self.scf_driver.xcfun
            if xcfun.is_hybrid():
                fock_type = '2jkx'
                exchange_scaling_factor = xcfun.get_frac_exact_exchange()
                fock_factor = 1.0
            else:
                fock_type = 'j'
                exchange_scaling_factor = 0.0
                fock_factor = 2.0
        else:
            fock_type = '2jk'
            exchange_scaling_factor = 1.0
            fock_factor = 1.0

        # TODO: range-separated Fock
        need_omega = (self._dft and self.scf_driver.xcfun.is_range_separated())
        if need_omega:
            assert_msg_critical(
                False, 'ScfHessianDriver: Not implemented for' +
                ' range-separated functional')

        den_mat_for_fock = make_matrix(ao_basis, mat_t.symmetric)
        den_mat_for_fock.set_values(density)

        den_mat_for_fock2 = make_matrix(ao_basis, mat_t.general)
        den_mat_for_fock2.set_values(density)

        screener = T4CScreener()
        screener.partition(ao_basis, molecule, 'eri')

        thresh_int = int(-math.log10(self.scf_driver.eri_thresh))

        if self.scf_driver.point_charges is not None:
            mm_coords = []
            mm_charges = []
            npoints = self.scf_driver.point_charges.shape[1]
            for p in range(npoints):
                xyz_p = self.scf_driver.point_charges[:3, p]
                chg_p = self.scf_driver.point_charges[3, p]
                mm_coords.append(xyz_p.copy())
                mm_charges.append(chg_p)
        else:
            mm_coords = None
            mm_charges = None

        # Parts related to second-order integral derivatives
        hessian_2nd_order_derivatives = np.zeros((natm, natm, 3, 3))

        # TODO: use alternative way to partition atoms
        local_atoms = list(range(natm))[self.rank::self.nodes]

        for i in local_atoms:

            ovlp_hess_200_mats = ovlp_hess_200_drv.compute(
                molecule, ao_basis, i)

            for x, label_x in enumerate('XYZ'):
                for y, label_y in enumerate('XYZ'):
                    ovlp_label = label_x + label_y if x <= y else label_y + label_x
                    ovlp_iixy = ovlp_hess_200_mats.matrix_to_numpy(ovlp_label)
                    hessian_2nd_order_derivatives[i, i, x, y] += 2.0 * (np.sum(
                        omega_ao * (ovlp_iixy + ovlp_iixy.T)))

            ovlp_hess_200_mats = Matrices()

            kin_hess_200_mats = kin_hess_200_drv.compute(molecule, ao_basis, i)

            for x, label_x in enumerate('XYZ'):
                for y, label_y in enumerate('XYZ'):
                    kin_label = label_x + label_y if x <= y else label_y + label_x
                    kin_200_iixy = kin_hess_200_mats.matrix_to_numpy(kin_label)
                    hessian_2nd_order_derivatives[i, i, x, y] += 2.0 * (np.sum(
                        density * (kin_200_iixy + kin_200_iixy.T)))

            kin_hess_200_mats = Matrices()

            npot_hess_200_mats = npot_hess_200_drv.compute(
                molecule, ao_basis, i)
            npot_hess_020_mats = npot_hess_020_drv.compute(
                molecule, ao_basis, i)

            for x, label_x in enumerate('XYZ'):
                for y, label_y in enumerate('XYZ'):
                    npot_label = label_x + label_y if x <= y else label_y + label_x
                    npot_200_iixy = npot_hess_200_mats.matrix_to_numpy(
                        npot_label)
                    npot_020_iixy = npot_hess_020_mats.matrix_to_numpy(
                        npot_label)
                    # TODO: move minus sign into function call (such as in oneints)
                    hessian_2nd_order_derivatives[i, i, x, y] += -2.0 * (np.sum(
                        density *
                        (npot_200_iixy + npot_200_iixy.T + npot_020_iixy)))

            npot_hess_200_mats = Matrices()
            npot_hess_020_mats = Matrices()

            if self.scf_driver.point_charges is not None:
                hmats_200 = npot_hess_200_drv.compute(molecule, ao_basis, i,
                                                      mm_coords, mm_charges)

                for x, label_x in enumerate('XYZ'):
                    for y, label_y in enumerate('XYZ'):
                        npot_label = label_x + label_y if x <= y else label_y + label_x
                        npot_200_iixy = hmats_200.matrix_to_numpy(npot_label)
                        # TODO: move minus sign into function call (such as in oneints)
                        hessian_2nd_order_derivatives[i, i, x, y] += -2.0 * (
                            np.sum(density * (npot_200_iixy + npot_200_iixy.T)))

                hmats_200 = Matrices()

            screener_atom_i = T4CScreener()
            screener_atom_i.partition_atom(ao_basis, molecule, 'eri', i)

            fock_hess_2000 = fock_hess_2000_drv.compute(
                ao_basis, screener_atom_i, screener, den_mat_for_fock,
                den_mat_for_fock2, i, fock_type, exchange_scaling_factor, 0.0,
                thresh_int)

            # 'XX', 'XY', 'XZ', 'YY', 'YZ', 'ZZ'
            xy_pairs_upper_triang = [
                (x, y) for x in range(3) for y in range(x, 3)
            ]

            for idx, (x, y) in enumerate(xy_pairs_upper_triang):
                hess_val = fock_factor * fock_hess_2000[idx]
                hessian_2nd_order_derivatives[i, i, x, y] += hess_val
                if x != y:
                    hessian_2nd_order_derivatives[i, i, y, x] += hess_val

        # do only upper triangular matrix
        all_atom_pairs = [(i, j) for i in range(natm) for j in range(i, natm)]

        # TODO: use alternative way to partition atom pairs
        local_atom_pairs = all_atom_pairs[self.rank::self.nodes]

        for i, j in local_atom_pairs:

            ovlp_hess_101_mats = ovlp_hess_101_drv.compute(
                molecule, ao_basis, i, j)

            for x, label_x in enumerate('XYZ'):
                for y, label_y in enumerate('XYZ'):
                    ovlp_label = f'{label_x}_{label_y}'
                    ovlp_ijxy = ovlp_hess_101_mats.matrix_to_numpy(ovlp_label)
                    hessian_2nd_order_derivatives[i, j, x, y] += 2.0 * (np.sum(
                        omega_ao * (ovlp_ijxy + ovlp_ijxy.T)))

            ovlp_hess_101_mats = Matrices()

            kin_hess_101_mats = kin_hess_101_drv.compute(
                molecule, ao_basis, i, j)

            for x, label_x in enumerate('XYZ'):
                for y, label_y in enumerate('XYZ'):
                    kin_label = f'{label_x}_{label_y}'
                    kin_101_ijxy = kin_hess_101_mats.matrix_to_numpy(kin_label)
                    hessian_2nd_order_derivatives[i, j, x, y] += 2.0 * (np.sum(
                        density * (kin_101_ijxy + kin_101_ijxy.T)))

            kin_hess_101_mats = Matrices()

            npot_hess_110_mats_ij = npot_hess_110_drv.compute(
                molecule, ao_basis, i, j)
            npot_hess_110_mats_ji = npot_hess_110_drv.compute(
                molecule, ao_basis, j, i)
            npot_hess_101_mats = npot_hess_101_drv.compute(
                molecule, ao_basis, i, j)

            for x, label_x in enumerate('XYZ'):
                for y, label_y in enumerate('XYZ'):
                    npot_xy_label = f'{label_x}_{label_y}'
                    npot_yx_label = f'{label_y}_{label_x}'
                    npot_110_ijxy = (
                        npot_hess_110_mats_ij.matrix_to_numpy(npot_xy_label) +
                        npot_hess_110_mats_ji.matrix_to_numpy(npot_yx_label))
                    npot_101_ijxy = npot_hess_101_mats.matrix_to_numpy(
                        npot_xy_label)
                    # TODO: move minus sign into function call (such as in oneints)
                    hessian_2nd_order_derivatives[i, j, x, y] += -2.0 * (np.sum(
                        density * (npot_110_ijxy + npot_110_ijxy.T +
                                   npot_101_ijxy + npot_101_ijxy.T)))

            npot_hess_110_mats_ij = Matrices()
            npot_hess_110_mats_ji = Matrices()
            npot_hess_101_mats = Matrices()

            if self.scf_driver.point_charges is not None:
                hmats_101 = npot_hess_101_drv.compute(molecule, ao_basis, i, j,
                                                      mm_coords, mm_charges)

                for x, label_x in enumerate('XYZ'):
                    for y, label_y in enumerate('XYZ'):
                        npot_xy_label = f'{label_x}_{label_y}'
                        npot_101_ijxy = hmats_101.matrix_to_numpy(npot_xy_label)
                        # TODO: move minus sign into function call (such as in oneints)
                        hessian_2nd_order_derivatives[i, j, x, y] += -2.0 * (
                            np.sum(density * (npot_101_ijxy + npot_101_ijxy.T)))

                hmats_101 = Matrices()

            screener_atom_pair = T4CScreener()
            screener_atom_pair.partition_atom_pair(ao_basis, molecule, 'eri', i,
                                                   j)

            fock_hess_1100 = fock_hess_1100_drv.compute(
                ao_basis, screener_atom_pair, screener, den_mat_for_fock,
                den_mat_for_fock2, i, j, fock_type, exchange_scaling_factor,
                0.0, thresh_int)

            screener_atom_i = T4CScreener()
            screener_atom_i.partition_atom(ao_basis, molecule, 'eri', i)

            screener_atom_j = T4CScreener()
            screener_atom_j.partition_atom(ao_basis, molecule, 'eri', j)

            # Note: use general matrix on both sides
            fock_hess_1010 = fock_hess_1010_drv.compute(
                ao_basis, screener_atom_i, screener_atom_j, den_mat_for_fock2,
                den_mat_for_fock2, i, j, fock_type, exchange_scaling_factor,
                0.0, thresh_int)

            # 'X_X', 'X_Y', 'X_Z', 'Y_X', 'Y_Y', 'Y_Z', 'Z_X', 'Z_Y', 'Z_Z'
            xy_pairs = [(x, y) for x in range(3) for y in range(3)]

            for idx, (x, y) in enumerate(xy_pairs):
                hessian_2nd_order_derivatives[i, j, x, y] += fock_factor * (
                    fock_hess_1100[idx] + fock_hess_1010[idx])

            # lower triangle is transpose of the upper part
            if i != j:
                hessian_2nd_order_derivatives[j, i] += (
                    hessian_2nd_order_derivatives[i, j].T)

        hessian_2nd_order_derivatives = self.comm.reduce(
            hessian_2nd_order_derivatives, root=mpi_master())

        # DFT:
        if self._dft:
            grid_drv = GridDriver(self.comm)
            grid_level = (get_default_grid_level(self.scf_driver.xcfun)
                          if self.scf_driver.grid_level is None else
                          self.scf_driver.grid_level)
            # make sure to use high grid level for DFT Hessian
            if molecule.get_charge() == 0:
                grid_level = max(6, grid_level)
            else:
                grid_level = max(7, grid_level)
            grid_drv.set_level(grid_level)
            mol_grid = grid_drv.generate(molecule)

            xc_mol_hess = XCMolecularHessian()
            hessian_dft_xc = xc_mol_hess.integrate_exc_hessian(
                molecule, ao_basis, [density], mol_grid,
                self.scf_driver.xcfun.get_func_label())
            hessian_dft_xc = self.comm.reduce(hessian_dft_xc, root=mpi_master())

        if self.rank == mpi_master():

            # Nuclear-nuclear repulsion contribution
            hessian_nuclear_nuclear = self.hess_nuc_contrib(molecule)

            # nuclei-point charges contribution
            if self.scf_driver.point_charges is not None:

                qm_coords = molecule.get_coordinates_in_bohr()
                nuclear_charges = molecule.get_element_ids()

                for i in range(natm):
                    z_a = nuclear_charges[i]
                    r_a = qm_coords[i]

                    for j in range(len(mm_charges)):
                        z_b = mm_charges[j]
                        r_b = mm_coords[j]

                        vec_ab = r_b - r_a
                        rab = np.linalg.norm(vec_ab)

                        for k in range(3):
                            for l in range(3):

                                hess_iikl = (3.0 * z_a * z_b * vec_ab[k] *
                                             vec_ab[l] / rab**5)

                                if k == l:
                                    hess_iikl -= z_a * z_b / rab**3

                                hessian_nuclear_nuclear[i, i, k, l] += hess_iikl

            # Sum up the terms and reshape for final Hessian
            self.hessian = (
                hessian_first_order_derivatives +
                hessian_2nd_order_derivatives.transpose(0, 2, 1, 3) +
                hessian_nuclear_nuclear.transpose(0, 2, 1, 3))

            self.hessian = self.hessian.reshape(natm * 3, natm * 3)

            if self._dft:
                self.hessian += hessian_dft_xc

        self.ostream.print_info('Second order derivative contributions' +
                                ' to the Hessian computed in' +
                                ' {:.2f} sec.'.format(tm.time() - t2))
        self.ostream.print_blank()
        self.ostream.flush()

        # Calculate the gradient of the dipole moment for IR intensities
        if self.do_dipole_gradient:
            self.compute_dipole_gradient(molecule, ao_basis, dist_cphf_ov)

    def compute_dipole_gradient(self, molecule, ao_basis, dist_cphf_ov):
        """
        Computes the analytical gradient of the dipole moment.

        :param molecule:
            The molecule.
        :param ao_basis:
            The AO basis set.
        """

        # Number of atoms and atomic charges
        natm = molecule.number_of_atoms()
        nuclear_charges = molecule.get_element_ids()

        # Dipole integrals
        dipole_mats = compute_electric_dipole_integrals(molecule, ao_basis,
                                                        [0.0, 0.0, 0.0])
        # Note: compute_dipole_gradient uses r instead of mu for dipole operator
        dipole_ints = (
            -1.0 * dipole_mats[0],
            -1.0 * dipole_mats[1],
            -1.0 * dipole_mats[2],
        )

        if self.rank == mpi_master():
            # Initialize a local dipole gradient to zero
            dipole_gradient = np.zeros((3, natm, 3))

            # Put the nuclear contributions to the right place
            natm_zeros = np.zeros(natm)
            dipole_gradient[0, :, :] = np.vstack(
                (nuclear_charges, natm_zeros, natm_zeros)).T
            dipole_gradient[1, :, :] = np.vstack(
                (natm_zeros, nuclear_charges, natm_zeros)).T
            dipole_gradient[2, :, :] = np.vstack(
                (natm_zeros, natm_zeros, nuclear_charges)).T

            scf_tensors = self.scf_driver.scf_tensors
            density = scf_tensors['D_alpha']
            mo = scf_tensors['C_alpha']
            nocc = molecule.number_of_alpha_electrons()
            mo_occ = mo[:, :nocc]
            mo_vir = mo[:, nocc:]
            nvir = mo_vir.shape[1]

        ovlp_grad_drv = OverlapGeom100Driver()
        dip_grad_drv = ElectricDipoleMomentGeom100Driver()

        for a in range(natm):

            # overlap gradient

            gmats = ovlp_grad_drv.compute(molecule, ao_basis, a)

            ovlp_deriv_ao_a = []
            for x, label in enumerate(['X', 'Y', 'Z']):
                gmat = gmats.matrix_to_numpy(label)
                ovlp_deriv_ao_a.append(gmat + gmat.T)

            gmats = Matrices()

            # perturbed density

            perturbed_density_a = []

            for x in range(3):

                cphf_ov_ax = dist_cphf_ov[a * 3 + x].get_full_vector(0)

                if self.rank == mpi_master():

                    cphf_ov_ax = cphf_ov_ax.reshape(nocc, nvir)

                    perturbed_density_a.append(
                        # mj,xyij,ni->xymn
                        -np.linalg.multi_dot(
                            [density, ovlp_deriv_ao_a[x], density])
                        # ma,xyia,ni->xymn
                        + np.linalg.multi_dot([mo_vir, cphf_ov_ax.T, mo_occ.T])
                        # mi,xyia,na->xymn
                        + np.linalg.multi_dot([mo_occ, cphf_ov_ax, mo_vir.T]))

            # dipole gradient

            gmats_dip = dip_grad_drv.compute(molecule, ao_basis,
                                             [0.0, 0.0, 0.0], a)

            gmats_dip_components = [
                'X_X', 'X_Y', 'X_Z', 'Y_X', 'Y_Y', 'Y_Z', 'Z_X', 'Z_Y', 'Z_Z'
            ]

            comp_to_idx = {'X': 0, 'Y': 1, 'Z': 2}

            for i, label in enumerate(gmats_dip_components):
                gmat_dip = gmats_dip.matrix_to_numpy(label)

                icoord = comp_to_idx[label[0]]  # atom coordinate component
                icomp = comp_to_idx[label[-1]]  # dipole operator component

                if self.rank == mpi_master():

                    c, x = icomp, icoord

                    dipole_gradient[c, a, x] += -2.0 * (
                        np.sum(density * (gmat_dip + gmat_dip.T)) +
                        np.sum(perturbed_density_a[x] * dipole_ints[c]))

        if self.rank == mpi_master():
            self.dipole_gradient = dipole_gradient.reshape(3, 3 * natm)

    def compute_orbital_response(self, molecule, ao_basis):
        """
        TEMPORARY FUNCTION FOR PERFORMANCE TESTING
        Computes the CPHF orbital response.

        :param molecule:
            The molecule.
        :param ao_basis:
            Tha AO basis.
        """

        # get SCF tensors
        scf_tensors = self.scf_driver.scf_tensors

        # Set up a CPHF solver
        cphf_solver = HessianOrbitalResponse(self.comm, self.ostream)
        cphf_solver.update_settings(self.cphf_dict, self.method_dict)

        # Solve the CPHF equations
        cphf_solver.compute(molecule, ao_basis, scf_tensors, self.scf_driver)
