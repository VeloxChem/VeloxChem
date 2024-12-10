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
from .griddriver import GridDriver
from .hessiandriver import HessianDriver
from .scfgradientdriver import ScfGradientDriver
from .firstorderprop import FirstOrderProperties
from .hessianorbitalresponse import HessianOrbitalResponse
from .profiler import Profiler
from .matrices import Matrices
from .dftutils import get_default_grid_level
from .errorhandler import assert_msg_critical
from .oneeints import compute_electric_dipole_integrals


class ScfHessianDriver(HessianDriver):
    """
    Implements SCF Hessian driver.

    :param scf_drv:
        The SCF driver.

    Instance variables
        - hessian: The Hessian in Hartree per Bohr**2.
        - flag: The type of Hessian driver.
        - do_pople_hessian: Evaluate the Hessian the Pople or
                the Ahlrichs/Furche way.
        - numerical_grad: Perform numerical gradient calculation.
        - perturbed_density: The perturbed density
    """

    def __init__(self, scf_drv, comm=None, ostream=None):
        """
        Initializes SCF Hessian driver.
        """

        super().__init__(scf_drv.comm, scf_drv.ostream)

        self.flag = 'SCF Hessian Driver'
        self.scf_driver = scf_drv

        self.numerical_grad = False

        self.do_pople_hessian = False

        self.perturbed_density = None

        # TODO TEMPORARY FLAG
        # Only run orbital response for performance testing
        self.orbrsp_only = False

        # flag for printing the Hessian
        self.do_print_hessian = False

        # TODO: determine _block_size_factor for SCF Hessian driver
        # self._block_size_factor = 4

        self._xcfun_ldstaging = scf_drv._xcfun_ldstaging

        self._input_keywords['hessian'].update({
                'do_pople_hessian': ('bool', 'whether to compute Pople Hessian'),
                'numerical_grad': ('bool', 'whether the gradient is numerical'),
                'orbrsp_only': ('bool', 'whether to only run CPHF orbital response'),
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

        # TODO: make use of profiler
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
            self.ostream.print_header('*** WARNING only computing Hessian orbital response!')
            self.compute_orbital_response(molecule, ao_basis)
            self.ostream.print_header('*** Hessian orbital response only: DONE  ***')
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

    # TODO: check if compute numerical works with MPI
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

        # atom labels
        labels = molecule.get_labels()

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

        self.scf_driver.restart = False
        scf_results = self.scf_driver.compute(molecule, ao_basis)
        assert_msg_critical(self.scf_driver.is_converged,
                            'ScfHessianDriver: SCF did not converge')
        energy_0 = self.scf_driver.get_scf_energy()

        for i in range(natm):
            for x in range(3):
                # Plus x
                coords[i, x] += self.delta_h
                new_mol = Molecule(labels, coords, units='au')
                new_mol.set_charge(charge)
                new_mol.set_multiplicity(multiplicity)
                scf_results = self.scf_driver.compute(new_mol, ao_basis)
                assert_msg_critical(self.scf_driver.is_converged,
                                    'ScfHessianDriver: SCF did not converge')
                energy_ixp = self.scf_driver.get_scf_energy()

                prop.compute_scf_prop(new_mol, ao_basis, scf_results)
                if self.rank == mpi_master():
                    mu_plus = prop.get_property('dipole moment')

                # Minus x
                coords[i, x] -= 2.0 * self.delta_h
                new_mol = Molecule(labels, coords, units='au')
                new_mol.set_charge(charge)
                new_mol.set_multiplicity(multiplicity)
                scf_results = self.scf_driver.compute(new_mol, ao_basis)
                assert_msg_critical(self.scf_driver.is_converged,
                                    'ScfHessianDriver: SCF did not converge')
                energy_ixm = self.scf_driver.get_scf_energy()

                prop.compute_scf_prop(new_mol, ao_basis, scf_results)
                if self.rank == mpi_master():
                    mu_minus = prop.get_property('dipole moment')

                if self.rank == mpi_master():
                    for c in range(3):
                        self.dipole_gradient[c, 3 * i + x] = (
                                (mu_plus[c] - mu_minus[c]) / (2.0 * self.delta_h))

                hessian[i, x, i, x] = ((energy_ixp - 2 * energy_0 + energy_ixm) /
                              self.delta_h**2)
                coords[i, x] += self.delta_h

                for j in range(i, natm):
                    for y in range(3):
                        if (j == i and x != y) or (j != i):
                            # Plus y
                            coords[j, y] += self.delta_h
                            new_mol = Molecule(labels, coords, units='au')
                            new_mol.set_charge(charge)
                            new_mol.set_multiplicity(multiplicity)
                            scf_results = self.scf_driver.compute(
                                new_mol, ao_basis)
                            assert_msg_critical(
                                self.scf_driver.is_converged,
                                'ScfHessianDriver: SCF did not converge')
                            energy_jyp = self.scf_driver.get_scf_energy()

                            # Plus x, plus y
                            coords[i, x] += self.delta_h
                            new_mol = Molecule(labels, coords, units='au')
                            new_mol.set_charge(charge)
                            new_mol.set_multiplicity(multiplicity)
                            scf_results = self.scf_driver.compute(
                                new_mol, ao_basis)
                            assert_msg_critical(
                                self.scf_driver.is_converged,
                                'ScfHessianDriver: SCF did not converge')
                            energy_ixp_jyp = self.scf_driver.get_scf_energy()
                            coords[i, x] -= self.delta_h

                            # Minus y
                            coords[j, y] -= 2.0 * self.delta_h
                            new_mol = Molecule(labels, coords, units='au')
                            new_mol.set_charge(charge)
                            new_mol.set_multiplicity(multiplicity)
                            scf_results = self.scf_driver.compute(
                                new_mol, ao_basis)
                            assert_msg_critical(
                                self.scf_driver.is_converged,
                                'ScfHessianDriver: SCF did not converge')
                            energy_jym = self.scf_driver.get_scf_energy()

                            # Minus x, minus y:
                            coords[i, x] -= self.delta_h
                            new_mol = Molecule(labels, coords, units='au')
                            new_mol.set_charge(charge)
                            new_mol.set_multiplicity(multiplicity)
                            scf_results = self.scf_driver.compute(
                                new_mol, ao_basis)
                            assert_msg_critical(
                                self.scf_driver.is_converged,
                                'ScfHessianDriver: SCF did not converge')
                            energy_ixm_jym = self.scf_driver.get_scf_energy()

                            coords[i, x] += self.delta_h
                            coords[j, y] += self.delta_h

                            hessian[i, x, j, y] = (
                                (energy_ixp_jyp - energy_ixp - energy_jyp +
                                 2 * energy_0 - energy_ixm - energy_jym +
                                 energy_ixm_jym) / (2 * self.delta_h**2))
                            hessian[j, y, i, x] = hessian[i, x, j, y]
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
        So far only for restricted Hartree-Fock with PySCF integral derivatives...

        :param molecule:
            The molecule.
        :param ao_basis:
            The AO basis set.
        :param profiler:
            The profiler.
        """

        natm = molecule.number_of_atoms()
        scf_tensors = self.scf_driver.scf_tensors

        if self.rank == mpi_master(): 
            density = scf_tensors['D_alpha']
            mo = scf_tensors['C_alpha']
            nao = mo.shape[0]
            nocc = molecule.number_of_alpha_electrons()
            mo_occ = mo[:, :nocc]
            mo_vir = mo[:, nocc:]
            nvir = mo_vir.shape[1]
            mo_energies = scf_tensors['E_alpha']
            eocc = mo_energies[:nocc]
            eoo = eocc.reshape(-1, 1) + eocc  # ei+ej
            omega_ao = - np.linalg.multi_dot([mo_occ, np.diag(eocc), mo_occ.T])

        else:
            density = None
        density = self.comm.bcast(density, root=mpi_master())
        nao = density.shape[0]

        # Solve CPHF equations
        cphf_solver = HessianOrbitalResponse(self.comm, self.ostream)
        cphf_solver.update_settings(self.cphf_dict, self.method_dict)
        cphf_solver.compute(molecule, ao_basis, scf_tensors)

        cphf_solution_dict = cphf_solver.cphf_results
        dist_cphf_ov = cphf_solution_dict['dist_cphf_ov']
        dist_cphf_rhs = cphf_solution_dict['dist_cphf_rhs']

        hessian_first_integral_derivatives = cphf_solution_dict['hessian_first_integral_derivatives']
        hessian_eri_overlap = cphf_solution_dict['hessian_eri_overlap']
        
        if self.rank == mpi_master():

            t1 = tm.time()
        
            # Parts related to first-order integral derivatives
            if self.do_pople_hessian:
                fock_uij = cphf_solution_dict['fock_uij']

                fock_deriv_ao = cphf_solution_dict['fock_deriv_ao']

                fock_deriv_oo = np.zeros((natm, 3, nocc, nocc))
                orben_ovlp_deriv_oo = np.zeros((natm, 3, nocc, nocc))
                for x in range(natm):
                    for y in range(3):
                        # mi,xymn,nj->xyij
                        fock_deriv_oo[x, y] = np.linalg.multi_dot([
                            mo_occ.T, fock_deriv_ao[x, y], mo_occ])
                        # ij,xyij->xyij (element-wise multiplication)
                        orben_ovlp_deriv_oo[x, y] = np.multiply(eoo, ovlp_deriv_oo[x,y])

        else:
            if self.do_pople_hessian:
                fock_uij = None
                fock_deriv_ao = None
                fock_deriv_oo = None
                orben_ovlp_deriv_oo = None

        if self.do_pople_hessian:
            fock_uij = self.comm.bcast(fock_uij, root=mpi_master())
            fock_deriv_ao = self.comm.bcast(fock_deriv_ao, root=mpi_master())
            fock_deriv_oo = self.comm.bcast(fock_deriv_oo, root=mpi_master())
            orben_ovlp_deriv_oo = self.comm.bcast(orben_ovlp_deriv_oo,
                                                 root=mpi_master())
            hessian_first_order_derivatives = self.compute_pople(molecule,
                                    ao_basis, -0.5 * ovlp_deriv_oo, cphf_ov,
                                    fock_uij, fock_deriv_oo,
                                    orben_ovlp_deriv_oo,
                                    self.perturbed_density, profiler)
        else:
            hessian_first_order_derivatives = self.compute_furche(molecule,
                                    ao_basis, dist_cphf_rhs, dist_cphf_ov,
                                    hessian_first_integral_derivatives,
                                    hessian_eri_overlap, profiler)

        # DFT:
        if self._dft:
            xc_mol_hess = XCMolecularHessian()
            hessian_dft_xc = xc_mol_hess.integrate_exc_hessian(molecule,
                                                ao_basis,
                                                [density], self.scf_driver._mol_grid,
                                                self.scf_driver.xcfun.get_func_label())
            hessian_dft_xc = self.comm.reduce(hessian_dft_xc, root=mpi_master())

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
            if self.scf_driver.xcfun.is_hybrid():
                fock_type = '2jkx'
                exchange_scaling_factor = self.scf_driver.xcfun.get_frac_exact_exchange()
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

        # TODO: parallelize over atoms
        if self.rank == mpi_master():
            t2 = tm.time()
            self.ostream.print_info('First order derivative contributions'
                                    + ' to the Hessian computed in' +
                                     ' {:.2f} sec.'.format(t2 - t1))
            self.ostream.print_blank()
            self.ostream.flush()

            # Parts related to second-order integral derivatives
            hessian_2nd_order_derivatives = np.zeros((natm, natm, 3, 3))

            for i in range(natm):

                ovlp_hess_200_mats = ovlp_hess_200_drv.compute(molecule, ao_basis, i)

                for x, label_x in enumerate('XYZ'):
                    for y, label_y in enumerate('XYZ'):
                        ovlp_label = label_x + label_y if x <= y else label_y + label_x
                        ovlp_iixy = ovlp_hess_200_mats.matrix_to_numpy(ovlp_label)
                        hessian_2nd_order_derivatives[i, i, x, y] += 2.0 * (
                                np.sum(omega_ao * (ovlp_iixy + ovlp_iixy.T)))

                ovlp_hess_200_mats = Matrices()

                kin_hess_200_mats = kin_hess_200_drv.compute(molecule, ao_basis, i)

                for x, label_x in enumerate('XYZ'):
                    for y, label_y in enumerate('XYZ'):
                        kin_label = label_x + label_y if x <= y else label_y + label_x
                        kin_200_iixy = kin_hess_200_mats.matrix_to_numpy(kin_label)
                        hessian_2nd_order_derivatives[i, i, x, y] += 2.0 * (
                                np.sum(density * (kin_200_iixy  + kin_200_iixy.T)))

                kin_hess_200_mats = Matrices()

                npot_hess_200_mats = npot_hess_200_drv.compute(molecule, ao_basis, i)
                npot_hess_020_mats = npot_hess_020_drv.compute(molecule, ao_basis, i)

                for x, label_x in enumerate('XYZ'):
                    for y, label_y in enumerate('XYZ'):
                        npot_label = label_x + label_y if x <= y else label_y + label_x
                        npot_200_iixy = npot_hess_200_mats.matrix_to_numpy(npot_label)
                        npot_020_iixy = npot_hess_020_mats.matrix_to_numpy(npot_label)
                        # TODO: move minus sign into function call (such as in oneints)
                        hessian_2nd_order_derivatives[i, i, x, y] += -2.0 * (
                                np.sum(density * (npot_200_iixy + npot_200_iixy.T + npot_020_iixy)))

                npot_hess_200_mats = Matrices()
                npot_hess_020_mats = Matrices()

                screener_atom = T4CScreener()
                screener_atom.partition_atom(ao_basis, molecule, 'eri', i)

                fock_hess_2000 = fock_hess_2000_drv.compute(
                        ao_basis, screener_atom, screener,
                        den_mat_for_fock, den_mat_for_fock2, i,
                        fock_type, exchange_scaling_factor,
                        0.0, thresh_int)

                # 'XX', 'XY', 'XZ', 'YY', 'YZ', 'ZZ'
                xy_pairs = [(x, y) for x in range(3) for y in range(x, 3)]

                for idx, (x, y) in enumerate(xy_pairs):
                    hess_val = fock_factor * fock_hess_2000[idx]
                    hessian_2nd_order_derivatives[i, i, x, y] += hess_val
                    if x != y:
                        hessian_2nd_order_derivatives[i, i, y, x] += hess_val

                # do only upper triangular matrix
                for j in range(i, natm):

                    ovlp_hess_101_mats = ovlp_hess_101_drv.compute(molecule, ao_basis, i, j)

                    for x, label_x in enumerate('XYZ'):
                        for y, label_y in enumerate('XYZ'):
                            ovlp_label = f'{label_x}_{label_y}'
                            ovlp_ijxy = ovlp_hess_101_mats.matrix_to_numpy(ovlp_label)
                            hessian_2nd_order_derivatives[i, j, x, y] += 2.0 * (
                                    np.sum(omega_ao * (ovlp_ijxy + ovlp_ijxy.T)))

                    ovlp_hess_101_mats = Matrices()

                    kin_hess_101_mats = kin_hess_101_drv.compute(molecule, ao_basis, i, j)

                    for x, label_x in enumerate('XYZ'):
                        for y, label_y in enumerate('XYZ'):
                            kin_label = f'{label_x}_{label_y}'
                            kin_101_ijxy = kin_hess_101_mats.matrix_to_numpy(kin_label)
                            hessian_2nd_order_derivatives[i, j, x, y] += 2.0 * (
                                    np.sum(density * (kin_101_ijxy  + kin_101_ijxy.T)))

                    kin_hess_101_mats = Matrices()

                    npot_hess_110_mats_ij = npot_hess_110_drv.compute(molecule, ao_basis, i, j)
                    npot_hess_110_mats_ji = npot_hess_110_drv.compute(molecule, ao_basis, j, i)
                    npot_hess_101_mats = npot_hess_101_drv.compute(molecule, ao_basis, i, j)

                    for x, label_x in enumerate('XYZ'):
                        for y, label_y in enumerate('XYZ'):
                            npot_xy_label = f'{label_x}_{label_y}'
                            npot_yx_label = f'{label_y}_{label_x}'
                            npot_110_ijxy = (npot_hess_110_mats_ij.matrix_to_numpy(npot_xy_label) +
                                             npot_hess_110_mats_ji.matrix_to_numpy(npot_yx_label))
                            npot_101_ijxy = npot_hess_101_mats.matrix_to_numpy(npot_xy_label)
                            # TODO: move minus sign into function call (such as in oneints)
                            hessian_2nd_order_derivatives[i, j, x, y] += -2.0 * (
                                np.sum(density * (npot_110_ijxy + npot_110_ijxy.T + npot_101_ijxy + npot_101_ijxy.T)))

                    npot_hess_110_mats_ij = Matrices()
                    npot_hess_110_mats_ji = Matrices()
                    npot_hess_101_mats = Matrices()

                    fock_hess_1100_mats = fock_hess_1100_drv.compute(ao_basis, molecule, den_mat_for_fock, i, j, fock_type, exchange_scaling_factor, 0.0)
                    fock_hess_1010_mats = fock_hess_1010_drv.compute(ao_basis, molecule, den_mat_for_fock, i, j, fock_type, exchange_scaling_factor, 0.0)

                    for x, label_x in enumerate('XYZ'):
                        for y, label_y in enumerate('XYZ'):
                            fock_label = f'{label_x}_{label_y}'
                            fock_hess_1100_mats_xy = fock_hess_1100_mats.matrix_to_numpy(fock_label)
                            fock_hess_1010_mats_xy = fock_hess_1010_mats.matrix_to_numpy(fock_label)
                            hessian_2nd_order_derivatives[i, j, x, y] += fock_factor * np.sum(
                                density * (fock_hess_1100_mats_xy + fock_hess_1010_mats_xy))

                    # TODO: atom-based screening
                    # TODO: in-place accumulation with two densities

                    fock_hess_1100_mats = Matrices()
                    fock_hess_1010_mats = Matrices()

                # lower triangle is transpose of the upper part
                for j in range(i):
                    hessian_2nd_order_derivatives[i,j] += (
                                hessian_2nd_order_derivatives[j,i].T )

            ## Nuclear-nuclear repulsion contribution
            hessian_nuclear_nuclear = self.hess_nuc_contrib(molecule)

            ## Sum up the terms and reshape for final Hessian
            self.hessian = ( hessian_first_order_derivatives
                            + hessian_2nd_order_derivatives
                            + hessian_nuclear_nuclear
                            ).transpose(0,2,1,3).reshape(3*natm, 3*natm)

            if self._dft:
                self.hessian += hessian_dft_xc

            t3 = tm.time()
            self.ostream.print_info('Second order derivative contributions'
                                    + ' to the Hessian computed in' +
                                     ' {:.2f} sec.'.format(t3 - t2))
            self.ostream.print_blank()
            self.ostream.flush()

        # Calculate the gradient of the dipole moment for IR intensities
        if self.do_dipole_gradient:
            self.compute_dipole_gradient(molecule, ao_basis, dist_cphf_ov)

    def compute_pople(self, molecule, ao_basis, cphf_oo, cphf_ov, fock_uij,
                      fock_deriv_oo, orben_ovlp_deriv_oo, perturbed_density,
                      profiler):
        """
        Computes the analytical nuclear Hessian the Pople way.
        Int. J. Quantum Chem. Quantum Chem. Symp. 13, 225-241 (1979).
        DOI: 10.1002/qua.560160825

        :param molecule:
            The molecule.
        :param ao_basis:
            The AO basis set.
        :param cphf_oo:
            The oo block of the CPHF coefficients.
        :param cphf_ov:
            The ov block of the CPHF coefficients.
        :param fock_uij:
            The auxiliary Fock matrix constructed
            using the oo block of the CPHF coefficients
            and two electron integrals.
        :param fock_deriv_oo:
            The oo block of the derivative of the Fock matrix
            with respect to nuclear coordinates
        :param orben_ovlp_deriv_oo:
            The oo block of the derivative of the overlap matrix
            with respect to nuclear coordinates, multiplied with
            orbital energies (ei+ej)S^chi_ij
        :param perturbed_density:
            The perturbed density matrix.
        :param profiler:
            The profiler.
        """

        natm = molecule.number_of_atoms()
        nocc = molecule.number_of_alpha_electrons()

        if self.rank == mpi_master():
            scf_tensors = self.scf_driver.scf_tensors
            mo = scf_tensors['C_alpha']
            mo_occ = mo[:, :nocc].copy()
            mo_vir = mo[:, nocc:].copy()
            nao = mo.shape[0]
            density = scf_tensors['D_alpha']
            mo_energies = scf_tensors['E_alpha']
            eocc = mo_energies[:nocc]
            eo_diag = np.diag(eocc)
            epsilon_dm_ao = - np.linalg.multi_dot([mo_occ, eo_diag, mo_occ.T])

            dof = 3

            # Construct the perturbed density matrix, and perturbed omega
            # TODO: consider if using the transpose makes the
            # computation faster; consider using cphf coefficients in AO
            # to compute the perturbed density matrix.

            orben_perturbed_density = np.zeros((natm, dof, nao, nao))
            mo_e_occ = np.multiply(mo_occ, eocc)
            for x in range(natm):
                for y in range(dof):
                    orben_perturbed_density[x, y] = (
                            # i,mj,xyij,ni->xymn
                            np.linalg.multi_dot([
                                mo_occ, cphf_oo[x, y].T, mo_e_occ.T])
                            # i,mi,xyij,nj->xymn
                            + np.linalg.multi_dot([
                                mo_e_occ, cphf_oo[x, y], mo_occ.T])
                            # i,ma,xyia,ni->xymn
                            + np.linalg.multi_dot([
                                mo_vir, cphf_ov[x, y].T, mo_e_occ.T])
                            # i,mi,xyia,na->xymn
                            + np.linalg.multi_dot([
                                mo_e_occ, cphf_ov[x, y], mo_vir.T])
                            )
        else:
            density = None

        if self._dft:
            grid_drv = GridDriver()
            grid_level = (get_default_grid_level(self.xcfun)
                          if self.scf_driver.grid_level is None else self.scf_driver.grid_level)
            grid_drv.set_level(grid_level)
            mol_grid = grid_drv.generate(molecule)

        density = self.comm.bcast(density, root=mpi_master())

        fock_uia_numpy = self.construct_fock_matrix_cphf(molecule, ao_basis,
                                                         cphf_ov)

        if self.rank == mpi_master():
            fock_cphf_oo = np.zeros((natm, dof, nocc, nocc))
            fock_cphf_ov = np.zeros((natm, dof, nocc, nocc))
            perturbed_omega_ao = np.zeros((natm, dof, nao, nao))
            for x in range(natm):
                for y in range(dof):
                    fock_cphf_oo[x, y] = (
                            # mi,xymn,nj->xyij
                            np.linalg.multi_dot([
                                mo_occ.T, fock_uij[x, y], mo_occ])
                            )
                    fock_cphf_ov[x, y] = (
                            # mi,xymn,nj->xyij
                            np.linalg.multi_dot([
                                mo_occ.T, fock_uia_numpy[x, y], mo_occ])
                            # mj,xymn,ni->xyij
                            + np.linalg.multi_dot([
                                mo_occ.T, fock_uia_numpy[x, y], mo_occ]).T
                            )

                    # Construct the derivative of the omega multipliers:
                    perturbed_omega_ao[x, y] = -1.0 * (
                            orben_perturbed_density[x,y]
                            # mi,xyij,nj->xymn
                            + np.linalg.multi_dot([
                                mo_occ, fock_deriv_oo[x, y], mo_occ.T])
                            # mi,xyij,nj->xymn
                            - 0.5 * np.linalg.multi_dot([
                                mo_occ, orben_ovlp_deriv_oo[x,y], mo_occ.T])
                            # mi,xyij,nj->xymn
                            + 2.0 * np.linalg.multi_dot([
                                mo_occ, fock_cphf_oo[x, y], mo_occ.T])
                            # mi,xyij,nj->xymn
                            + np.linalg.multi_dot([
                                mo_occ, fock_cphf_ov[x, y], mo_occ.T])
                            )

            # First integral derivatives: partial Fock and overlap matrix
            # derivatives
            hessian_first_integral_derivatives = np.zeros((natm, natm, 3, 3))

        if self._dft: 
            xc_mol_hess = XCMolecularHessian()
        for i in range(natm):
            # upper triangular part
            for j in range(i, natm):
                if self._dft:
                    # First derivative of the Vxc matrix elements
                    vxc_deriv_j = xc_mol_hess.integrate_vxc_fock_gradient(
                                    molecule, ao_basis, [density], mol_grid,
                                    self.scf_driver.xcfun.get_func_label(), j)
                    vxc_deriv_j = self.comm.reduce(vxc_deriv_j,
                                                   root=mpi_master())
                if self.rank == mpi_master():
                    # First derivative of the Fock matrix
                    fock_deriv_j =  fock_deriv(molecule, ao_basis, density, j,
                                               self.scf_driver)
                    if self._dft:
                        fock_deriv_j += vxc_deriv_j

                    # First derivative of overlap matrix
                    ovlp_deriv_j = overlap_deriv(molecule, ao_basis, j)
                    # Add the contribution of the perturbed density matrix
                    for x in range(dof):
                        for y in range(dof):
                            hessian_first_integral_derivatives[i, j, x, y] += (
                                    # xmn,ymn->xy
                                    np.linalg.multi_dot([2.0 * perturbed_density[i, x].reshape(nao**2),
                                                         fock_deriv_j[y].reshape(nao**2)])
                                    # xmn,ymn->xy
                                    + np.linalg.multi_dot([2.0 * perturbed_omega_ao[i, x].reshape(nao**2),
                                                           ovlp_deriv_j[y].reshape(nao**2)])
                                    )

            if self.rank == mpi_master():
                # lower triangular part
                for j in range(i):
                    hessian_first_integral_derivatives[i,j] += (
                                     hessian_first_integral_derivatives[j,i].T )

        if self.rank == mpi_master():
            return hessian_first_integral_derivatives
        else:
            return None


    def compute_furche(self, molecule, ao_basis, dist_cphf_rhs, dist_cphf_ov,
                       hessian_first_integral_derivatives,
                       hessian_eri_overlap, profiler):
        """
        Computes the analytical nuclear Hessian the Furche/Ahlrichs way.
        Chem. Phys. Lett. 362, 511–518 (2002).
        DOI: 10.1016/S0009-2614(02)01084-9

        :param molecule:
            The molecule.
        :param ao_basis:
            The AO basis set.
        :param dist_cphf_rhs:
            The distributed RHS of the CPHF equations in MO basis.
        :param dist_cphf_ov:
            The distributed ov block of the CPHF coefficients.
        :param profiler:
            The profiler.
        """

        natm = molecule.number_of_atoms()

        # RHS contracted with CPHF coefficients (ov)
        if self.rank == mpi_master():
            hessian_cphf_coeff_rhs = np.zeros((natm, natm, 3, 3))

        for i in range(natm):
            for x in range(3):
                cphf_ov_ix = dist_cphf_ov[i*3+x].get_full_vector(0)

                for j in range(i, natm):
                    for y in range(3):
                        cphf_rhs_jy = dist_cphf_rhs[j*3+y].get_full_vector(0)

                        if self.rank == mpi_master():
                            hess_ijxy = 4.0 * (np.dot(cphf_ov_ix, cphf_rhs_jy))
                            hessian_cphf_coeff_rhs[i,j,x,y] += hess_ijxy
                            if i != j:
                                hessian_cphf_coeff_rhs[j,i,y,x] += hess_ijxy

        # return the sum of the three contributions
        if self.rank == mpi_master():
            return ( hessian_cphf_coeff_rhs
                   + hessian_first_integral_derivatives
                   + hessian_eri_overlap)
        else:
            return None


    def construct_fock_matrix_cphf(self, molecule, ao_basis, cphf_ov):
        """
        Contracts the CPHF coefficients with the two-electron
        integrals and returns an auxiliary fock matrix as a
        numpy array.

        :param molecule:
            The Molecule.
        :param ao_basis:
            The AO Basis.
        :param chpf_ov:
            The occupied-virtual block of the CPHF coefficients.
        """
        natm = molecule.number_of_atoms()
        nocc = molecule.number_of_alpha_electrons()

        if self.rank == mpi_master():
            scf_tensors = self.scf_driver.scf_tensors
            mo = scf_tensors['C_alpha']
            mo_occ = mo[:, :nocc].copy()
            mo_vir = mo[:, nocc:].copy()
            nao = mo.shape[0]
            density = scf_tensors['D_alpha']
            mo_energies = scf_tensors['E_alpha']
            eocc = mo_energies[:nocc]
            eo_diag = np.diag(eocc)
            epsilon_dm_ao = - np.linalg.multi_dot([mo_occ, eo_diag, mo_occ.T])
            # Transform the CPHF coefficients to AO:
            uia_ao = np.zeros((natm, 3, nao, nao))
            for x in range(natm):
                for y in range(3):
                    uia_ao[x,y] = np.linalg.multi_dot([
                        mo_occ, cphf_ov[x,y], mo_vir.T])
            uia_ao = uia_ao.reshape(3*natm, nao, nao)
            
            # create AODensity and Fock matrix objects, contract with ERI
            uia_ao_list = [uia_ao[x] for x in range(natm * 3)]
        else:
            scf_tensors = None
            uia_ao_list = None

        uia_ao_list = self.comm.bcast(uia_ao_list, root=mpi_master())

        # We use comp_lr_fock from CphfSolver to compute the eri
        # and xc contributions
        cphf_solver = HessianOrbitalResponse(self.comm, self.ostream)
        cphf_solver.update_settings(self.cphf_dict, self.method_dict)
        # ERI information
        eri_dict = cphf_solver._init_eri(molecule, ao_basis)
        # DFT information
        dft_dict = cphf_solver._init_dft(molecule, scf_tensors)
        # PE information
        pe_dict = cphf_solver._init_pe(molecule, ao_basis)

        profiler = Profiler({
            'timing': self.timing,
            'profiling': self.profiling,
            'memory_profiling': self.memory_profiling,
            'memory_tracing': self.memory_tracing,
        })

        fock_uia = cphf_solver._comp_lr_fock(uia_ao_list, molecule, ao_basis,
                                 eri_dict, dft_dict, pe_dict, profiler)

        if self.rank == mpi_master():
            # TODO: can this be done in a different way?
            fock_uia_numpy = np.zeros((natm,3,nao,nao))
            for i in range(natm):
                for x in range(3):
                    fock_uia_numpy[i,x] = fock_uia[3*i + x]

            return fock_uia_numpy
        else:
            return None

    # FIXME the function below is unused
    def compute_perturbed_energy_weighted_density_matrix(self, molecule,
                                                         ao_basis):
        """
        Calculates the perturbed energy weighted density matrix
        and returns it as a numpy array.

        :param molecule:
            The molecule.
        :param ao_basis:
            The AO basis set.
        """
        natm = molecule.number_of_atoms()
        nocc = molecule.number_of_alpha_electrons()
        scf_tensors = self.scf_driver.scf_tensors
        mo = scf_tensors['C']
        mo_occ = mo[:, :nocc].copy()
        mo_vir = mo[:, nocc:].copy()
        nao = mo.shape[0]
        nmo = mo.shape[1]
        nvir = nmo - nocc
        mo_energies = scf_tensors['E']
        eocc = mo_energies[:nocc]
        eoo = eocc.reshape(-1, 1) + eocc #ei+ej
        # Set up a CPHF solver
        cphf_solver = HessianOrbitalResponse(self.comm, self.ostream)
        cphf_solver.update_settings(self.cphf_dict, self.method_dict)

        # Solve the CPHF equations
        cphf_solver.compute(molecule, ao_basis, scf_tensors, self.scf_driver)

        # Extract the relevant results
        cphf_solution_dict = cphf_solver.cphf_results
        cphf_ov = cphf_solution_dict['cphf_ov'].reshape(natm, 3, nocc, nvir)
        ovlp_deriv_oo = cphf_solution_dict['ovlp_deriv_oo']
        cphf_oo = -0.5 * ovlp_deriv_oo

        fock_uij = cphf_solution_dict['fock_uij']
        fock_deriv_ao = cphf_solution_dict['fock_deriv_ao']
        fock_deriv_oo = np.einsum('mi,xymn,nj->xyij', mo_occ,
                                   fock_deriv_ao, mo_occ)
        fock_uia_numpy = self.construct_fock_matrix_cphf(molecule, ao_basis,
                                                         cphf_ov)

        # (ei+ej)S^chi_ij
        orben_ovlp_deriv_oo = np.einsum('ij,xyij->xyij', eoo, ovlp_deriv_oo)

        fock_cphf_oo = np.einsum('mi,xymn,nj->xyij', mo_occ, fock_uij, mo_occ)

        fock_cphf_ov = ( np.einsum('mi,xymn,nj->xyij', mo_occ,
                                   fock_uia_numpy, mo_occ)
                        +np.einsum('mj,xymn,ni->xyij', mo_occ,
                                   fock_uia_numpy, mo_occ)
                        )

        orben_perturbed_density = ( np.einsum('i,mj,xyij,ni->xymn',
                                            eocc, mo_occ, cphf_oo, mo_occ)
                                  + np.einsum('i,mi,xyij,nj->xymn',
                                            eocc, mo_occ, cphf_oo, mo_occ)
                                  + np.einsum('i,ma,xyia,ni->xymn',
                                            eocc, mo_vir, cphf_ov, mo_occ)
                                  +np.einsum('i,mi,xyia,na->xymn',
                                            eocc, mo_occ, cphf_ov, mo_vir)
                                 )

        # Construct the derivative of the omega multipliers:
        perturbed_omega_ao = - ( orben_perturbed_density
                                + np.einsum('mi,xyij,nj->xymn', mo_occ,
                                            fock_deriv_oo, mo_occ)
                                -0.5*np.einsum('mi,xyij,nj->xymn', mo_occ,
                                                orben_ovlp_deriv_oo, mo_occ)
                                + 2*np.einsum('mi,xyij,nj->xymn', mo_occ,
                                            fock_cphf_oo, mo_occ)
                                + np.einsum('mi,xyij,nj->xymn', mo_occ,
                                            fock_cphf_ov, mo_occ)
                                )
        return perturbed_omega_ao

    # TODO: make this option available
    def compute_numerical_with_analytical_gradient(self, molecule, ao_basis,
                                                   profiler=None):

        """
        Performs calculation of numerical Hessian.

        :param molecule:
            The molecule.
        :param ao_basis:
            The AO basis set.
        :param profiler:
            The profiler.
        """

        if profiler is None:
            profiler = Profiler({
                'timing': self.timing,
                'profiling': self.profiling,
                'memory_profiling': self.memory_profiling,
                'memory_tracing': self.memory_tracing,
            })


        # settings dictionary for gradient driver
        grad_dict = dict(self.hess_dict)
        if self.numerical_grad:
            grad_dict['numerical'] = 'yes'
            warn_msg = '*** Warning: Numerical Hessian will be calculated '
            warn_msg += 'based on numerical gradient.'
            self.ostream.print_header(warn_msg.ljust(56))
            warn_msg = '  This takes a long time and has limited accuracy.'
            self.ostream.print_header(warn_msg.ljust(56))
            self.ostream.print_blank()
            self.ostream.flush()
        else:
            grad_dict['numerical'] = 'no'

        self.ostream.mute()

        # atom labels
        labels = molecule.get_labels()

        # number of atoms
        natm = molecule.number_of_atoms()

        # atom coordinates (nx3)
        coords = molecule.get_coordinates()

        # gradient driver
        grad_drv = ScfGradientDriver(self.scf_driver)
        grad_drv.update_settings(grad_dict, self.method_dict)

        # number of atomic orbitals
        scf_tensors = self.scf_driver.scf_tensors
        nao = scf_tensors['D_alpha'].shape[0]

        # Hessian
        hessian = np.zeros((natm, 3, natm, 3))

        # First-order properties for gradient of dipole moment
        prop = FirstOrderProperties(self.comm, self.ostream)
        # numerical gradient (3 dipole components, no. atoms x 3 atom coords)
        self.dipole_gradient = np.zeros((3, 3 * natm))

        if not self.do_four_point:
            for i in range(natm):
                for d in range(3):
                    coords[i, d] += self.delta_h
                    new_mol = Molecule(labels, coords, units='au')
                    self.scf_driver.compute(new_mol, ao_basis)
                    grad_drv.compute(new_mol, ao_basis)
                    grad_plus = grad_drv.get_gradient()

                    prop.compute_scf_prop(new_mol, ao_basis,
                                          self.scf_driver.scf_tensors)
                    mu_plus = prop.get_property('dipole moment')

                    coords[i, d] -= 2.0 * self.delta_h
                    new_mol = Molecule(labels, coords, units='au')
                    self.scf_driver.compute(new_mol, ao_basis)
                    grad_drv.compute(new_mol, ao_basis)
                    grad_minus = grad_drv.get_gradient()

                    prop.compute_scf_prop(new_mol, ao_basis,
                                          self.scf_driver.scf_tensors)
                    mu_minus = prop.get_property('dipole moment')

                    for c in range(3):
                        self.dipole_gradient[c, 3*i + d] = (
                            (mu_plus[c] - mu_minus[c]) / (2.0 * self.delta_h) )
                    coords[i, d] += self.delta_h
                    hessian[i, d, :, :] = (
                           (grad_plus - grad_minus) / (2.0 * self.delta_h) )
        else:
            # Four-point numerical derivative approximation
            # for debugging of analytical Hessian:
            # [ f(x - 2h) - 8 f(x - h) + 8 f(x + h) - f(x + 2h) ] / ( 12h )
            for i in range(natm):
                for d in range(3):
                    coords[i, d] += self.delta_h
                    new_mol = Molecule(labels, coords, units='au')
                    self.scf_driver.compute(new_mol, ao_basis)
                    grad_drv.compute(new_mol, ao_basis)
                    grad_plus1 = grad_drv.get_gradient()

                    prop.compute_scf_prop(new_mol, ao_basis,
                                          self.scf_driver.scf_tensors)
                    mu_plus1 = prop.get_property('dipole moment')

                    coords[i, d] += self.delta_h
                    new_mol = Molecule(labels, coords, units='au')
                    self.scf_driver.compute(new_mol, ao_basis)
                    grad_drv.compute(new_mol, ao_basis)
                    grad_plus2 = grad_drv.get_gradient()

                    prop.compute_scf_prop(new_mol, ao_basis,
                                          self.scf_driver.scf_tensors)
                    mu_plus2 = prop.get_property('dipole moment')

                    coords[i, d] -= 3.0 * self.delta_h
                    new_mol = Molecule(labels, coords, units='au')
                    self.scf_driver.compute(new_mol, ao_basis)
                    grad_drv.compute(new_mol, ao_basis)
                    grad_minus1 = grad_drv.get_gradient()

                    prop.compute_scf_prop(new_mol, ao_basis,
                                          self.scf_driver.scf_tensors)
                    mu_minus1 = prop.get_property('dipole moment')

                    coords[i, d] -= self.delta_h
                    new_mol = Molecule(labels, coords, units='au')
                    self.scf_driver.compute(new_mol, ao_basis)
                    grad_drv.compute(new_mol, ao_basis)
                    grad_minus2 = grad_drv.get_gradient()

                    prop.compute_scf_prop(new_mol, ao_basis,
                                          self.scf_driver.scf_tensors)
                    mu_minus2 = prop.get_property('dipole moment')

                    for c in range(3):
                        self.dipole_gradient[c, 3*i + d] = (
                          ( mu_minus2[c] - 8.0 * mu_minus1[c] 
                        + 8.0 * mu_plus1[c] - mu_plus2[c])
                        / (12.0 * self.delta_h) )
                    coords[i, d] += 2.0 * self.delta_h
                    # f'(x) ~ [ f(x - 2h) - 8 f(x - h) + 8 f(x + h)
                    # - f(x + 2h) ] / ( 12h )
                    hessian[i, d] = (grad_minus2 - 8.0 * grad_minus1
                       + 8.0 * grad_plus1 - grad_plus2) / (12.0 * self.delta_h)

        # reshaped Hessian as member variable
        self.hessian = hessian.reshape(3*natm, 3*natm)

        self.scf_driver.compute(molecule, ao_basis)
        self.ostream.unmute()

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
            dipole_gradient[0,:,:] = np.vstack((nuclear_charges,
                                                natm_zeros,
                                                natm_zeros)).T
            dipole_gradient[1,:,:] = np.vstack((natm_zeros,
                                                nuclear_charges,
                                                natm_zeros)).T
            dipole_gradient[2,:,:] = np.vstack((natm_zeros,
                                                natm_zeros,
                                                nuclear_charges)).T

            scf_tensors = self.scf_driver.scf_tensors 
            density = scf_tensors['D_alpha']
            nao = density.shape[0]
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

                cphf_ov_ax = dist_cphf_ov[a*3+x].get_full_vector(0)

                if self.rank == mpi_master():

                    cphf_ov_ax = cphf_ov_ax.reshape(nocc, nvir)

                    perturbed_density_a.append(
                        # mj,xyij,ni->xymn
                        - np.linalg.multi_dot([density, ovlp_deriv_ao_a[x], density])
                        # ma,xyia,ni->xymn
                        + np.linalg.multi_dot([mo_vir, cphf_ov_ax.T, mo_occ.T]) 
                        # mi,xyia,na->xymn
                        + np.linalg.multi_dot([mo_occ, cphf_ov_ax, mo_vir.T]) 
                        )

            # dipole gradient

            gmats_dip = dip_grad_drv.compute(molecule, ao_basis, [0.0, 0.0, 0.0], a)

            gmats_dip_components = ['X_X', 'X_Y', 'X_Z', 'Y_X', 'Y_Y', 'Y_Z',
                                    'Z_X', 'Z_Y', 'Z_Z']

            comp_to_idx = {'X': 0, 'Y': 1, 'Z': 2}

            for i, label in enumerate(gmats_dip_components):
                gmat_dip = gmats_dip.matrix_to_numpy(label)

                icoord = comp_to_idx[label[0]] # atom coordinate component
                icomp = comp_to_idx[label[-1]] # dipole operator component

                if self.rank == mpi_master():

                    c, x = icomp, icoord

                    dipole_gradient[c,a,x] += -2.0 * (
                            np.sum(density * (gmat_dip + gmat_dip.T)) +
                            np.sum(perturbed_density_a[x] * dipole_ints[c]))

        if self.rank == mpi_master():
            self.dipole_gradient = dipole_gradient.reshape(3, 3 * natm)

    def compute_dipole_integral_derivatives(self, molecule, ao_basis):
        """
        Imports the analytical derivatives of dipole integrals.

        :param molecule:
            The molecule.
        :param ao_basis:
            The AO basis set.

        :return:
            The dipole integral derivatives.
        """

        # number of atoms
        natm = molecule.number_of_atoms()

        # number of atomic orbitals
        nao = ao_basis.get_dimensions_of_basis()
        
        dip_grad_drv = ElectricDipoleMomentGeom100Driver()
        dipole_integrals_gradient = np.zeros((3, natm, 3, nao, nao))

        for iatom in range(natm):
            gmats_dip = dip_grad_drv.compute(molecule, ao_basis, [0.0, 0.0, 0.0], iatom)

            # the keys of the dipole gmat
            gmats_dip_components = (['X_X', 'X_Y', 'X_Z', 'Y_X' ]
                                  + [ 'Y_Y', 'Y_Z', 'Z_X', 'Z_Y', 'Z_Z'])

            # dictionary to convert from string idx to integer idx
            comp_to_idx = {'X': 0, 'Y': 1, 'Z': 2}

            for i, label in enumerate(gmats_dip_components):
                gmat_dip = gmats_dip.matrix_to_numpy(label)
                gmat_dip += gmat_dip.T

                icoord = comp_to_idx[label[0]] # atom coordinate component
                icomp = comp_to_idx[label[-1]] # dipole operator component

                # reorder indices to first is operator comp, second is coord
                dipole_integrals_gradient[icomp, iatom, icoord] += gmat_dip

        return dipole_integrals_gradient

    def compute_orbital_response(self, molecule, ao_basis):
        """
        TEMPORARY FUNCTION FOR PERFORMACE TESTING
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

