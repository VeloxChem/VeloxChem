#
#                           VELOXCHEM 1.0-RC3
#         ----------------------------------------------------
#                     An Electronic Structure Code
#
#  Copyright © 2018-2022 by VeloxChem developers. All rights reserved.
#  Contact: https://veloxchem.org/contact
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
import numpy as np
import time as tm
import sys

from .molecule import Molecule
from .gradientdriver import GradientDriver
from .hessiandriver import HessianDriver
from .scfgradientdriver import ScfGradientDriver
from .outputstream import OutputStream
from .firstorderprop import FirstOrderProperties
from .cphfsolver import CphfSolver
from .hessianorbitalresponse import HessianOrbitalResponse
from .lrsolver import LinearResponseSolver
from .cppsolver import ComplexResponse
from .polarizabilitygradient import PolarizabilityGradient
from .profiler import Profiler
from .qqscheme import get_qq_scheme
from .dftutils import get_default_grid_level
from .veloxchemlib import mpi_master
from .veloxchemlib import ElectronRepulsionIntegralsDriver
from .veloxchemlib import denmat
from .veloxchemlib import fockmat
from .veloxchemlib import AODensityMatrix
from .veloxchemlib import AOFockMatrix
from .veloxchemlib import ElectricDipoleIntegralsDriver
from .veloxchemlib import xcfun
from .veloxchemlib import XCMolecularHessian
from .veloxchemlib import GridDriver
from .errorhandler import assert_msg_critical
from .dftutils import get_default_grid_level
from .inputparser import parse_input
from .sanitychecks import dft_sanity_check

# For PySCF integral derivatives
from .import_from_pyscf import overlap_deriv
from .import_from_pyscf import fock_deriv
from .import_from_pyscf import eri_deriv
from .import_from_pyscf import overlap_second_deriv
from .import_from_pyscf import hcore_second_deriv
from .import_from_pyscf import eri_second_deriv
from .import_from_pyscf import dipole_deriv

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
#        - numerical: Perform numerical Hessian calculation.
        - numerical_grad: Perform numerical gradient calculation.
#        - do_four_point: Perform four-point numerical approximation.
#        - delta_h: Nuclear displacement for finite differences.
#        - do_print_hessian: Flag for printing the Hessian.
#        - do_dipole_gradient: The gradient of the dipole moment.
#        - dipole_gradient: The gradient of the dipole moment.
        - perturbed_density: The perturbed density
    """

    def __init__(self, scf_drv, comm=None, ostream=None):
        """
        Initializes SCF Hessian driver.
        """

        super().__init__(scf_drv.comm, scf_drv.ostream)

        #if comm is None:
        #    comm = MPI.COMM_WORLD

        #if ostream is None:
        #    if comm.Get_rank() == mpi_master():
        #        ostream = OutputStream(sys.stdout)
        #    else:
        #        ostream = OutputStream(None)

        ## MPI information
        #self.comm = comm
        #self.rank = self.comm.Get_rank()
        #self.nodes = self.comm.Get_size()

        #self.ostream = ostream

        self.flag = 'SCF Hessian Driver'
        self.scf_driver = scf_drv

#        self.hessian = None

#        self.numerical = False
#        self.do_four_point = False
#        self.delta_h = 0.001

        self.numerical_grad = False

        self.do_pople_hessian = False

        self.perturbed_density = None

#        self.do_dipole_gradient = False
#        self.dipole_gradient = None

        # flag for printing the Hessian
        self.do_print_hessian = False

#        self._dft = False
#        self.grid_level = None
#        self.xcfun = None

        # Timing and profiling
#        self.timing = False
#        self.profiling = False
#        self.memory_profiling = False
#        self.memory_tracing = False

        self._input_keywords['hessian'].update({
        #self._input_keywords = {
            #'hessian': {
                'do_pople_hessian': ('bool', 'whether to compute Pople Hessian'),
                'numerical_grad': ('bool', 'whether the gradient is numerical'),
#                'do_dipole_gradient': ('bool', 'whether to compute the dipole gradient'),
#                'do_print_hessian': ('bool', 'whether to print the Hessian'),
#                'timing': ('bool', 'whether timing is needed'),
#                'profiling': ('bool', 'whether profiling is needed'),
#                'memory_profiling': ('bool', 'whether to profile memory'),
#                'memory_tracing': ('bool', 'whether to trace memory'),
            #    },
#             'method_settings': {
#                'xcfun': ('str_upper', 'exchange-correlation functional'),
#                'grid_level': ('int', 'accuracy level of DFT grid'),
            #    }
            })

    def update_settings(self, method_dict, hess_dict=None, cphf_dict=None):#,
                        #rsp_dict=None, polgrad_dict=None):
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

        #if method_dict is None:
        #    method_dict = {}
        #if hess_dict is None:
        #    hess_dict = {}
        if cphf_dict is None:
            cphf_dict = {}

        #hess_keywords = {
        #    key: val[0] for key, val in
        #    self._input_keywords['hessian'].items()
        #}

        #parse_input(self, hess_keywords, hess_dict)

        #method_keywords = {
        #    key: val[0]
        #    for key, val in self._input_keywords['method_settings'].items()
        #}

        #parse_input(self, method_keywords, method_dict)

        #dft_sanity_check(self, 'update_settings')

        #self.method_dict = dict(method_dict)
        #self.hess_dict = dict(hess_dict)
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

        # TODO: find a better place for the profiler?
        profiler = Profiler({
            'timing': self.timing,
            'profiling': self.profiling,
            'memory_profiling': self.memory_profiling,
            'memory_tracing': self.memory_tracing,
        })

        # Sanity check
        if (self.do_pople_hessian) and (not self._dft):
            raise ValueError('Pople Hessian only valid for DFT')

        # Save the electronic energy
        #self.elec_energy = self.scf_driver.get_scf_energy()

        if self.numerical:
            self.compute_numerical(molecule, ao_basis)
        else:
            self.compute_analytical(molecule, ao_basis, profiler)

        # Calculate the gradient of the dipole moment for IR intensities
        if self.do_dipole_gradient and (self.perturbed_density is not None):
            self.compute_dipole_gradient(molecule, ao_basis)

        # Calculate the analytical polarizability gradient for Raman intensities
        #if self.do_raman:
        #    self.compute_polarizability_gradient(molecule, ao_basis)

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
            mo_energies = scf_tensors['E']
            eocc = mo_energies[:nocc]
            eoo = eocc.reshape(-1, 1) + eocc #ei+ej
            omega_ao = - np.linalg.multi_dot([mo_occ, np.diag(eocc), mo_occ.T])

            gs_density = AODensityMatrix([density], denmat.rest)
        else:
            gs_density = AODensityMatrix()

        gs_density.broadcast(self.rank, self.comm)

        # Set up a CPHF solver
        cphf_solver = HessianOrbitalResponse(self.comm, self.ostream)
        cphf_solver.update_settings(self.cphf_dict, self.method_dict)

        # Solve the CPHF equations
        cphf_solver.compute(molecule, ao_basis, scf_tensors, self.scf_driver)
        
        if self.rank == mpi_master():
            hessian_first_order_derivatives = np.zeros((natm, natm, 3, 3))
            cphf_solution_dict = cphf_solver.cphf_results
            cphf_ov = cphf_solution_dict['cphf_ov'].reshape(natm, 3, nocc, nvir)
            ovlp_deriv_oo = cphf_solution_dict['ovlp_deriv_oo']

            dof = 3
            perturbed_density = np.zeros((natm, dof, nao, nao))
            for x in range(natm):
                for y in range(dof):
                    perturbed_density[x, y] = (
                            # mj,xyij,ni->xymn
                            - np.linalg.multi_dot([mo_occ, ovlp_deriv_oo[x, y].T, mo_occ.T])
                            # ma,xyia,ni->xymn
                            + np.linalg.multi_dot([mo_vir, cphf_ov[x, y].T, mo_occ.T]) 
                            # mi,xyia,na->xymn
                            + np.linalg.multi_dot([mo_occ, cphf_ov[x, y], mo_vir.T]) 
                            )

            t1 = tm.time()
        
            # Parts related to first-order integral derivatives
            if self.do_pople_hessian:
                fock_uij = cphf_solution_dict['fock_uij']

                fock_deriv_ao = cphf_solution_dict['fock_deriv_ao']

                fock_deriv_oo = np.zeros((natm, dof, nocc, nocc))
                orben_ovlp_deriv_oo = np.zeros((natm, dof, nocc, nocc))
                for x in range(natm):
                    for y in range(dof):
                        # mi,xymn,nj->xyij
                        fock_deriv_oo[x, y] = np.linalg.multi_dot([
                            mo_occ.T, fock_deriv_ao[x, y], mo_occ])
                        # ij,xyij->xyij (element-wise multiplication)
                        orben_ovlp_deriv_oo[x, y] = np.multiply(eoo, ovlp_deriv_oo[x,y])
            else:
                cphf_rhs = cphf_solution_dict['cphf_rhs'].reshape(natm, 3,
                                                                  nocc, nvir)
        else:
            ovlp_deriv_oo = None
            cphf_ov = None
            perturbed_density = None

            if self.do_pople_hessian:
                fock_uij = None
                fock_deriv_ao = None
                fock_deriv_oo = None
                orben_ovlp_deriv_oo = None
            else:
                cphf_rhs = None

        ovlp_deriv_oo = self.comm.bcast(ovlp_deriv_oo, root=mpi_master())
        cphf_ov = self.comm.bcast(cphf_ov, root=mpi_master())
        perturbed_density = self.comm.bcast(perturbed_density,
                                            root=mpi_master())

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
                                    perturbed_density, profiler)
        else:
            cphf_rhs = self.comm.bcast(cphf_rhs, root=mpi_master())
            hessian_first_order_derivatives = self.compute_furche(molecule,
                                    ao_basis, cphf_rhs, -0.5 * ovlp_deriv_oo,
                                    cphf_ov, profiler)
        # amount of exact exchange
        frac_K = 1.0

        # DFT:
        if self._dft:
            if self.scf_driver.xcfun.is_hybrid():
                frac_K = self.scf_driver.xcfun.get_frac_exact_exchange()
            else:
                frac_K = 0.0

            grid_drv = GridDriver()
            grid_level = (get_default_grid_level(self.xcfun)
                          if self.scf_driver.grid_level is None else self.scf_driver.grid_level)
            grid_drv.set_level(grid_level)

            xc_mol_hess = XCMolecularHessian()
            mol_grid = grid_drv.generate(molecule)

            hessian_dft_xc = xc_mol_hess.integrate_exc_hessian(molecule,
                                                ao_basis,
                                                gs_density, mol_grid,
                                                self.scf_driver.xcfun.get_func_label())
            hessian_dft_xc = self.comm.reduce(hessian_dft_xc, root=mpi_master())

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
                # do only upper triangular matrix
                for j in range(i, natm):
                    # Get integral second-order derivatives
                    ovlp_2nd_deriv_ii, ovlp_2nd_deriv_ij = overlap_second_deriv(
                                                    molecule, ao_basis, i, j)
                    hcore_2nd_deriv_ij = hcore_second_deriv(
                                            molecule, ao_basis, i, j)
                    eri_2nd_deriv_ii, eri_2nd_deriv_ij = eri_second_deriv(
                                            molecule, ao_basis, i, j)
                    if i == j:
                        # Add diagonal (same atom) contributions, 2S + 2J - K

                        # Build Fock matrices by contractig the density matrix
                        # with the second order derivatives of the two-electron
                        # integrals 

                        aux_ii_Fock_2nd_deriv_j = np.zeros((3, 3, nao, nao))
                        aux_ii_Fock_2nd_deriv_k = np.zeros((3, 3, nao, nao))
                        for x in range(dof):
                            for y in range(dof):
                                # mn,xymn->xy
                                hessian_2nd_order_derivatives[i, i, x, y] += 2.0 * (
                                        np.linalg.multi_dot([
                                            omega_ao.reshape(nao**2), 
                                            ovlp_2nd_deriv_ii[x, y].reshape(nao**2)])
                                        )
                                # kl,xymnkl->xymn
                                aux_ii_Fock_2nd_deriv_j[x, y] = np.linalg.multi_dot([
                                    density.reshape(nao**2), 
                                    (eri_2nd_deriv_ii[x, y].transpose(2,3,0,1)).reshape(nao**2,nao**2)
                                    ]).reshape(nao, nao)
                                # kl,xymknl->xymn
                                aux_ii_Fock_2nd_deriv_k[x, y] = np.linalg.multi_dot([
                                    density.reshape(nao**2), 
                                    (eri_2nd_deriv_ii[x, y].transpose(1,3,0,2)).reshape(nao**2,nao**2)
                                    ]).reshape(nao, nao)
                                # mn,xymn->xy
                                hessian_2nd_order_derivatives[i, i, x, y] += 2.0 * (
                                        np.linalg.multi_dot([
                                            density.reshape(nao**2), 
                                            aux_ii_Fock_2nd_deriv_j[x, y].reshape(nao**2)
                                        ]))
                                # mn,xymn->xy
                                hessian_2nd_order_derivatives[i, i, x, y] -= frac_K * (
                                        np.linalg.multi_dot([
                                            density.reshape(nao**2), 
                                            aux_ii_Fock_2nd_deriv_k[x, y].reshape(nao**2)
                                        ]))

                    # Add non-diagonal contributions, 2S + 2J - K + 2h

                    aux_ij_Fock_2nd_deriv_j = np.zeros((3, 3, nao, nao))
                    aux_ij_Fock_2nd_deriv_k = np.zeros((3, 3, nao, nao))
                    for x in range(dof):
                        for y in range(dof):
                            # kl,xymnkl->xymn
                            aux_ij_Fock_2nd_deriv_j[x, y] = np.linalg.multi_dot([
                                density.reshape(nao**2), 
                                (eri_2nd_deriv_ij[x, y].transpose(2, 3, 0, 1)).reshape(nao**2, nao**2)
                                ]).reshape(nao, nao)
                            # kl,xymknl->xymn
                            aux_ij_Fock_2nd_deriv_k[x, y] = np.linalg.multi_dot([
                                density.reshape(nao**2), 
                                (eri_2nd_deriv_ij[x, y].transpose(1, 3, 0, 2)).reshape(nao**2, nao**2)
                                ]).reshape(nao, nao)
                            # mn, xymn->xy
                            hessian_2nd_order_derivatives[i, j, x, y] += 2.0 * (
                                    np.linalg.multi_dot([
                                        density.reshape(nao**2),
                                        aux_ij_Fock_2nd_deriv_j[x, y].reshape(nao**2)]))
                            # mn,xymn->xy
                            hessian_2nd_order_derivatives[i, j, x, y] -= frac_K * (
                                    np.linalg.multi_dot([
                                        density.reshape(nao**2),
                                        aux_ij_Fock_2nd_deriv_k[x, y].reshape(nao**2)]))
                            # mn,xymn->xy
                            hessian_2nd_order_derivatives[i, j, x, y] += 2.0 * (
                                    np.linalg.multi_dot([
                                        omega_ao.reshape(nao**2),
                                        ovlp_2nd_deriv_ij[x, y].reshape(nao**2)]))
                            # mn,xymn->xy
                            hessian_2nd_order_derivatives[i, j, x, y] += 2.0 * (
                                    np.linalg.multi_dot([
                                        density.reshape(nao**2),
                                        hcore_2nd_deriv_ij[x, y].reshape(nao**2)]))

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

            # save perturbed density as instance variable: needed for dipole
            # gradient
            self.perturbed_density = perturbed_density

            t3 = tm.time()
            self.ostream.print_info('Second order derivative contributions'
                                    + ' to the Hessian computed in' +
                                     ' {:.2f} sec.'.format(t3 - t2))
            self.ostream.print_blank()
            self.ostream.flush()

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
            orbital energies (ei+ej)S^\chi_ij
        :param perturbed_density:
            The perturbed density matrix.
        :param profiler:
            The profiler.
        """

        natm = molecule.number_of_atoms()
        nocc = molecule.number_of_alpha_electrons()

        if self.rank == mpi_master():
            scf_tensors = self.scf_driver.scf_tensors
            mo = scf_tensors['C']
            mo_occ = mo[:, :nocc].copy()
            mo_vir = mo[:, nocc:].copy()
            nao = mo.shape[0]
            density = scf_tensors['D_alpha']
            mo_energies = scf_tensors['E']
            eocc = mo_energies[:nocc]
            eo_diag = np.diag(eocc)
            epsilon_dm_ao = - np.linalg.multi_dot([mo_occ, eo_diag, mo_occ.T])

            dof = 3

            # Construct the perturbed density matrix, and perturbed omega
            # TODO: consider if using the transpose makes the
            # computation faster; consider using cphf coefficients in AO
            # to compute the perturbed density matrix.
            gs_density = AODensityMatrix([density], denmat.rest)

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
            gs_density = AODensityMatrix()

        grid_drv = GridDriver()
        grid_level = (get_default_grid_level(self.xcfun)
                      if self.scf_driver.grid_level is None else self.scf_driver.grid_level)
        grid_drv.set_level(grid_level)

        mol_grid = grid_drv.generate(molecule)
        gs_density.broadcast(self.rank, self.comm)

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

        if self.scf_driver._dft: 
            xc_mol_hess = XCMolecularHessian()
        for i in range(natm):
            # upper triangular part
            for j in range(i, natm):
                if self.scf_driver._dft:
                    # First derivative of the Vxc matrix elements
                    vxc_deriv_j = xc_mol_hess.integrate_vxc_fock_gradient(
                                    molecule, ao_basis, gs_density, mol_grid,
                                    self.scf_driver.xcfun.get_func_label(), j)
                    vxc_deriv_j = self.comm.reduce(vxc_deriv_j,
                                                   root=mpi_master())
                if self.rank == mpi_master():
                    # First derivative of the Fock matrix
                    fock_deriv_j =  fock_deriv(molecule, ao_basis, density, j,
                                               self.scf_driver)
                    if self.scf_driver._dft:
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


    def compute_furche(self, molecule, ao_basis, cphf_rhs, cphf_oo, cphf_ov,
                       profiler):
        """
        Computes the analytical nuclear Hessian the Furche/Ahlrichs way.
        Chem. Phys. Lett. 362, 511–518 (2002).
        DOI: 10.1016/S0009-2614(02)01084-9

        :param molecule:
            The molecule.
        :param ao_basis:
            The AO basis set.
        :param cphf_rhs:
            The RHS of the CPHF equations in MO basis.
        :param cphf_oo:
            The oo block of the CPHF coefficients.
        :param cphf_ov:
            The ov block of the CPHF coefficients.
        :param profiler:
            The profiler.
        """

        natm = molecule.number_of_atoms()
        nocc = molecule.number_of_alpha_electrons()
        if self.rank == mpi_master():
            scf_tensors = self.scf_driver.scf_tensors
            mo = scf_tensors['C']
            mo_occ = mo[:, :nocc].copy()
            mo_vir = mo[:, nocc:].copy()
            nvir = mo_vir.shape[1]
            nao = mo.shape[0]
            density = scf_tensors['D_alpha']
            mo_energies = scf_tensors['E']
            eocc = mo_energies[:nocc]
            omega_ao = - np.linalg.multi_dot([mo_occ, np.diag(eocc), mo_occ.T])
            gs_density = AODensityMatrix([density], denmat.rest)
        else:
            density = None
            gs_density = AODensityMatrix()
            scf_tensors = None

        if self.scf_driver._dft:
            grid_drv = GridDriver()
            grid_level = (get_default_grid_level(self.xcfun)
                          if self.scf_driver.grid_level is None else self.scf_driver.grid_level)
            grid_drv.set_level(grid_level)

            mol_grid = grid_drv.generate(molecule)
        gs_density.broadcast(self.rank, self.comm)

        density = self.comm.bcast(density, root=mpi_master())
        scf_tensors = self.comm.bcast(scf_tensors, root=mpi_master())

        # We use comp_lr_fock from CphfSolver to compute the eri
        # and xc contributions
        #cphf_solver = CphfSolver(self.comm, self.ostream)
        cphf_solver = HessianOrbitalResponse(self.comm, self.ostream)
        cphf_solver.update_settings(self.cphf_dict, self.method_dict)
        # ERI information
        eri_dict = cphf_solver._init_eri(molecule, ao_basis)
        # DFT information
        dft_dict = cphf_solver._init_dft(molecule, scf_tensors)
        # PE information
        pe_dict = cphf_solver._init_pe(molecule, ao_basis)

        # TODO: the MPI is not done properly here,
        # fix once the new integral code is ready.
        if self.rank == mpi_master():
            # RHS contracted with CPHF coefficients (ov)

            hessian_cphf_coeff_rhs = np.zeros((natm, natm, 3, 3))
            for i in range(natm):
                for j in range(natm):
                    for x in range(3):
                        for y in range(3):
                            hessian_cphf_coeff_rhs[i,j,x,y] = 4.0 * (
                            np.linalg.multi_dot([cphf_ov[i,x].reshape(nocc*nvir),
                                                 cphf_rhs[j,y].reshape(nocc*nvir)])
                            )

            # First integral derivatives: partial Fock and overlap
            # matrix derivatives
            # TODO why is this initiated here?
            hessian_first_integral_derivatives = np.zeros((natm, natm, 3, 3))

        if self.scf_driver._dft: 
            xc_mol_hess = XCMolecularHessian()
        for i in range(natm):
            # First derivative of the Vxc matrix elements
            if self.scf_driver._dft:
                vxc_deriv_i = xc_mol_hess.integrate_vxc_fock_gradient(
                                molecule, ao_basis, gs_density, mol_grid,
                                self.scf_driver.xcfun.get_func_label(), i)
                vxc_deriv_i = self.comm.reduce(vxc_deriv_i, root=mpi_master())
            if self.rank == mpi_master():
                fock_deriv_i =  fock_deriv(molecule, ao_basis, density, i,
                                           self.scf_driver)
                if self.scf_driver._dft:
                    fock_deriv_i += vxc_deriv_i
                ovlp_deriv_i = overlap_deriv(molecule, ao_basis, i)

            # upper triangular part
            for j in range(i, natm):
                # First derivative of the Vxc matrix elements
                if self.scf_driver._dft:
                    vxc_deriv_j = xc_mol_hess.integrate_vxc_fock_gradient(
                                    molecule, ao_basis, gs_density, mol_grid,
                                    self.scf_driver.xcfun.get_func_label(), j)
                    vxc_deriv_j = self.comm.reduce(vxc_deriv_j,
                                                    root=mpi_master())
                if self.rank == mpi_master():
                    fock_deriv_j =  fock_deriv(molecule, ao_basis, density, j,
                                               self.scf_driver)
                    if self.scf_driver._dft:
                        fock_deriv_j += vxc_deriv_j
                    ovlp_deriv_j = overlap_deriv(molecule, ao_basis, j)
        
                    Fix_Sjy = np.zeros((3,3))
                    Fjy_Six = np.zeros((3,3))
                    Six_Sjy = np.zeros((3,3))
                    for x in range(3):
                        for y in range(3):
                            Fix_Sjy[x,y] = np.trace(np.linalg.multi_dot([
                                    density, fock_deriv_i[x], density,
                                    ovlp_deriv_j[y]]))
                            Fjy_Six[x,y] = np.trace(np.linalg.multi_dot([
                                    density, fock_deriv_j[y],
                                    density, ovlp_deriv_i[x]]))
                            Six_Sjy[x,y] = ( 2*np.trace(np.linalg.multi_dot([
                                    omega_ao, ovlp_deriv_i[x], density,
                                    ovlp_deriv_j[y]])) )
                    hessian_first_integral_derivatives[i,j] += -2 * (Fix_Sjy 
                                                            + Fjy_Six + Six_Sjy)
   
            if self.rank == mpi_master(): 
                # lower triangular part
                for j in range(i):
                    hessian_first_integral_derivatives[i,j] += (
                            hessian_first_integral_derivatives[j,i].T )

        if self.rank == mpi_master():
            # Overlap derivative with ERIs
            hessian_eri_overlap = np.zeros((natm, natm, 3, 3))
        for i in range(natm):
            if self.rank == mpi_master():
                ovlp_deriv_i = overlap_deriv(molecule, ao_basis, i)
            # upper triangular part
            for j in range(i, natm):
                if self.rank == mpi_master():
                    ovlp_deriv_j = overlap_deriv(molecule, ao_basis, j)

                    # Overlap derivative contracted with two density matrices
                    P_P_Six = np.zeros((3, nao, nao))
                    P_P_Sjy = np.zeros((3, nao, nao))
                    for x in range(3):
                        # mn,xnk,kl->xml
                        P_P_Six[x] = np.linalg.multi_dot([
                        density, ovlp_deriv_i[x], density])
                        # mn,xnk,kl->xml
                        P_P_Sjy[x] = np.linalg.multi_dot([
                            density, ovlp_deriv_j[x], density])

                    # Create a list of 2D numpy arrays to create
                    # AODensityMatrix objects
                    # and calculate auxiliary Fock matrix with ERI driver
                    P_P_Six_ao_list = list([P_P_Six[x] for x in range(3)])
                    P_P_Six_dm_ao = AODensityMatrix(P_P_Six_ao_list,
                                                    denmat.rest)
                else:
                    P_P_Six_dm_ao = AODensityMatrix()

                P_P_Six_dm_ao.broadcast(self.rank, self.comm)

                P_P_Six_fock_ao = AOFockMatrix(P_P_Six_dm_ao)
                # MPI issue
                cphf_solver._comp_lr_fock(P_P_Six_fock_ao, P_P_Six_dm_ao,
                                          molecule, ao_basis, eri_dict,
                                          dft_dict, pe_dict, profiler)

                if self.rank == mpi_master():
                    # Convert the auxiliary Fock matrices to numpy arrays 
                    # for further use
                    np_P_P_Six_fock = np.zeros((3,nao,nao))
                    for k in range(3):
                        np_P_P_Six_fock[k] = P_P_Six_fock_ao.to_numpy(k)

                    for x in range(3):
                        for y in range(3):
                            # xmn,ymn->xy
                            hessian_eri_overlap[i,j] += 2.0 * (np.linalg.multi_dot([
                                np_P_P_Six_fock[x].reshape(nao**2), 
                                P_P_Sjy[y].reshape(nao**2)]))

            # lower triangular part
            for j in range(i):
                if self.rank == mpi_master():
                    hessian_eri_overlap[i,j] += hessian_eri_overlap[j,i].T

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
            mo = scf_tensors['C']
            mo_occ = mo[:, :nocc].copy()
            mo_vir = mo[:, nocc:].copy()
            nao = mo.shape[0]
            density = scf_tensors['D_alpha']
            mo_energies = scf_tensors['E']
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
            uia_ao_list = list([uia_ao[x] for x in range(natm * 3)])
            ao_density_uia = AODensityMatrix(uia_ao_list, denmat.rest)
        else:
            ao_density_uia = AODensityMatrix()

        ao_density_uia.broadcast(self.rank, self.comm)

        fock_uia = AOFockMatrix(ao_density_uia)

        fock_flag = fockmat.rgenjk
        for i in range(natm*3):
            fock_uia.set_fock_type(fock_flag, i)

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

        cphf_solver._comp_lr_fock(fock_uia, ao_density_uia, molecule, ao_basis,
                                 eri_dict, dft_dict, pe_dict, profiler)

        if self.rank == mpi_master():
            # TODO: can this be done in a different way?
            fock_uia_numpy = np.zeros((natm,3,nao,nao))
            for i in range(natm):
                for x in range(3):
                    fock_uia_numpy[i,x] = fock_uia.to_numpy(3*i + x)

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

        # (ei+ej)S^\chi_ij
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

    def compute_dipole_gradient(self, molecule, ao_basis):
        """
        Computes the analytical gradient of the dipole moment.

        :param molecule:
            The molecule.
        :param ao_basis:
            The AO basis set.
        :param scf_drv:
            The SCF driver.
        :param perturbed_density: #TODO remove this if setup approved
            The perturbed density matrix.
        """

        # Number of atoms and atomic charges
        natm = molecule.number_of_atoms()
        nuclear_charges = molecule.elem_ids_to_numpy()

        scf_tensors = self.scf_driver.scf_tensors 
        perturbed_density = self.perturbed_density

        density = scf_tensors['D_alpha']
        nao = perturbed_density.shape[-1] # need for reshaping

        # Dipole integrals
        dipole_drv = ElectricDipoleIntegralsDriver(self.comm)
        dipole_mats = dipole_drv.compute(molecule, ao_basis)
        dipole_ints = np.array((dipole_mats.x_to_numpy(),
                                dipole_mats.y_to_numpy(),
                                dipole_mats.z_to_numpy()))

        # Initialize a local dipole gradient to zero
        dipole_gradient = np.zeros((3, natm, 3))

        # Put the nuclear contributions to the right place
        natm_zeros = np.zeros((natm))
        dipole_gradient[0] = np.vstack((nuclear_charges,
                                        natm_zeros, natm_zeros)).T
        dipole_gradient[1] = np.vstack((natm_zeros,
                                        nuclear_charges, natm_zeros)).T
        dipole_gradient[2] = np.vstack((natm_zeros,
                                        natm_zeros, nuclear_charges)).T

        # TODO: replace once veloxchem analytical integral derivatives 
        # are available
        # TODO: replace the full array of dipole integrals derivatives with
        # the derivative with respect to each atom and contract with 
        # the densities in the same way as done for the energy gradient.
        dipole_integrals_deriv = self.compute_dipole_integral_derivatives(
                                                        molecule, ao_basis)

        # Add the electronic contributions
        for a in range(natm):
            for c in range(3):
                for x in range(3):
                    dipole_gradient[c,a,x] += -2.0 * (
                            np.linalg.multi_dot([
                                density.reshape(nao**2),
                                dipole_integrals_deriv[c,a,x].reshape(nao**2)])
                            + np.linalg.multi_dot([
                                perturbed_density[a,x].reshape(nao**2),
                                dipole_ints[c].reshape(nao**2)]))

        self.dipole_gradient = dipole_gradient.reshape(3, 3 * natm)


    # TODO can stay here
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

        scf_tensors = self.scf_driver.scf_tensors
        # number of atomic orbitals
        nao = scf_tensors['D_alpha'].shape[0]

        # 3 dipole components x No. atoms x 3 atomic coordinates
        # x No. basis x No. basis
        dipole_integrals_gradient = np.zeros((3, natm, 3, nao, nao))

        for i in range(natm):
            dipole_integrals_gradient[:,i,:,:,:] = (
                                dipole_deriv(molecule, ao_basis, i)
                                )

        return dipole_integrals_gradient

