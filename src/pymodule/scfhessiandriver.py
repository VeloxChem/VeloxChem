#
#                           VELOXCHEM 1.0-RC2
#         ----------------------------------------------------
#                     An Electronic Structure Code
#
#  Copyright © 2018-2021 by VeloxChem developers. All rights reserved.
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
from .scfgradientdriver import ScfGradientDriver
from .hessiandriver import HessianDriver
from .outputstream import OutputStream
from .scffirstorderprop import ScfFirstOrderProperties
from .lrsolver import LinearResponseSolver
from .profiler import Profiler
from .qqscheme import get_qq_scheme
from .veloxchemlib import mpi_master
from .veloxchemlib import ElectronRepulsionIntegralsDriver
from .veloxchemlib import denmat
from .veloxchemlib import fockmat
from .veloxchemlib import AODensityMatrix
from .veloxchemlib import AOFockMatrix
from .veloxchemlib import ElectricDipoleIntegralsDriver

class ScfHessianDriver(HessianDriver):
    """
    Implements SCF Hessian driver.

    :param scf_drv:
        The SCF driver.
    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.

    Instance variables
        - scf_drv: The SCF driver.
        - do_raman: Additionally calculate Raman intensities
                (at significantly higher computational cost).
    """

    def __init__(self, scf_drv, comm=None, ostream=None):
        """
        Initializes SCF Hessian driver.
        """

        super().__init__(comm, ostream)

        self.flag = 'SCF Hessian Driver'
        self.scf_drv = scf_drv
        self.do_raman = False
        self.pople = False

        # Solver setup
        self.conv_thresh = 1.0e-4
        self.max_iter = 50
        self.iter_count = 0
        self.is_converged = False


    def update_settings(self, method_dict, freq_dict=None, cphf_dict=None):
        """
        Updates settings in ScfHessianDriver.

        :param method_dict:
            The input dictionary of method settings group.
        :param freq_dict:
            The input dictionary of Hessian/frequency settings group.
        :param cphf_dict:
            The input dictionary of CPHF (orbital response) settings.
        """

        super().update_settings(method_dict, freq_dict)

        if freq_dict is None:
            freq_dict = {}

        # Settings for orbital response module
        if cphf_dict is None:
            cphf_dict = {}
        self.cphf_dict = dict(cphf_dict)

        if 'conv_thresh' in cphf_dict:
            self.conv_thresh = float(cphf_dict['conv_thresh'])

        if 'max_iter' in cphf_dict:
            self.max_iter = int(cphf_dict['max_iter']) 

        # check if Raman intensities are to be calculated (numerically)
        if 'do_raman' in freq_dict:
            key = freq_dict['do_raman'].lower()
            self.do_raman = True if key in ['yes', 'y'] else False

        # The electronic energy
        self.elec_energy = self.scf_drv.get_scf_energy()


    def compute(self, molecule, ao_basis, min_basis=None):
        """
        Computes the analytical or numerical nuclear Hessian.

        :param molecule:
            The molecule.
        :param ao_basis:
            The AO basis set.
        :param min_basis:
            The minimal AO basis set.
        """
        self.print_header()

        start_time = tm.time()

        profiler = Profiler({
            'timing': self.timing,
            'profiling': self.profiling,
            'memory_profiling': self.memory_profiling,
            'memory_tracing': self.memory_tracing,
        })

        #if self.numerical:
        self.compute_numerical(molecule, ao_basis, min_basis, profiler)
        #else:
        #    raise ValueError("Only numerical Hessian currently implemented.")

        # print Hessian
        if self.do_print_hessian is True:
            self.print_geometry(molecule)
            self.ostream.print_blank()
            self.print_hessian(molecule)

        valstr = '*** Time spent in Hessian calculation: '
        valstr += '{:.2f} sec ***'.format(tm.time() - start_time)
        self.ostream.print_header(valstr)
        self.ostream.print_blank()
        self.ostream.print_blank()
        self.ostream.flush()


    def compute_numerical(self, molecule, ao_basis, min_basis, profiler):
        """
        Performs the calculation of a numerical Hessian based only
        on the energy.

        :param molecule:
            The molecule.
        :param ao_basis:
            The AO basis set.
        :param min_basis:
            The minimal AO basis set.
        :param profiler:
            The profiler.
        """

        scf_ostream_state = self.scf_drv.ostream.state
        self.scf_drv.ostream.state = False

        # atom labels
        labels = molecule.get_labels()

        # number of atoms
        natm = molecule.number_of_atoms()

        # atom coordinates (nx3)
        coords = molecule.get_coordinates()

        # number of atomic orbitals
        nao = self.scf_drv.scf_tensors['D_alpha'].shape[0]

        # Hessian
        hessian = np.zeros((natm, 3, natm, 3))

        # First-order properties for gradient of dipole moment
        prop = ScfFirstOrderProperties(self.comm, self.ostream)
        # numerical gradient (3 dipole components, no. atoms x 3 atom coords)
        #self.dipole_gradient = np.zeros((3, natm, 3))
        self.dipole_gradient = np.zeros((3, 3 * natm))

        # If Raman intensities are calculated,
        # set up LR solver and member variable
        if self.do_raman:
            # linear response driver for polarizability calculation
            lr_drv = LinearResponseSolver(self.comm, self.scf_drv.ostream)
            # polarizability: 3 coordinates x 3 coordinates
            # (ignoring frequencies)
            # polarizability gradient: dictionary goes through
            # 3 coordinates x 3 coordinates
            # each entry having values for no. atoms x 3 coordinates
            self.pol_gradient = np.zeros((3, 3, 3 * natm))
            # dictionary to translate from numbers to operator components 'xyz'
            component_dict = {0: 'x', 1: 'y', 2: 'z'}

        self.scf_drv.compute(molecule, ao_basis, min_basis)
        energy_0 = self.scf_drv.get_scf_energy()

        for i in range(natm):
            for x in range(3):
                # Plus x
                coords[i, x] += self.delta_h
                new_mol = Molecule(labels, coords, units='au')
                self.scf_drv.compute(new_mol, ao_basis, min_basis)
                energy_ixp = self.scf_drv.get_scf_energy()

                prop.compute(new_mol, ao_basis,
                             self.scf_drv.scf_tensors)
                mu_plus = prop.get_property('dipole moment')

                if self.do_raman:
                    lr_drv.is_converged = False
                    lr_results_p = lr_drv.compute(new_mol, ao_basis,
                                                  self.scf_drv.scf_tensors)

                # Minus x
                coords[i, x] -= 2.0 * self.delta_h
                new_mol = Molecule(labels, coords, units='au')
                self.scf_drv.compute(new_mol, ao_basis, min_basis)
                energy_ixm = self.scf_drv.get_scf_energy()

                prop.compute(new_mol, ao_basis,
                             self.scf_drv.scf_tensors)
                mu_minus = prop.get_property('dipole moment')

                if self.do_raman:
                    lr_drv.is_converged = False
                    lr_results_m = lr_drv.compute(new_mol, ao_basis,
                                                    self.scf_drv.scf_tensors)
                    for aop in range(3):
                        for bop in range(3):
                            # TODO: here the freq. is hard-coded to 0.0!
                            comp_plus = ( lr_results_p['response_functions'][
                                component_dict[aop], component_dict[bop], 0.0] )
                            comp_minus = ( lr_results_m['response_functions'][
                                component_dict[aop], component_dict[bop], 0.0] )
                            self.pol_gradient[aop, bop, 3*i + x] = (
                                ( comp_plus - comp_minus) /
                                    (2.0 * self.delta_h) )

                for c in range(3):
                    self.dipole_gradient[c, 3*i + x] = (
                        (mu_plus[c] - mu_minus[c]) / (2.0 * self.delta_h) )

                hessian[i, x, i, x] = (
                  ( energy_ixp - 2 * energy_0 + energy_ixm ) / self.delta_h**2
                                        )
                coords[i, x] += self.delta_h

                for j in range(i, natm):
                    for y in range(3):
                        if ( j == i and x != y ) or ( j != i ):
                            # Plus y
                            coords[j, y] += self.delta_h
                            new_mol = Molecule(labels, coords, units='au')
                            self.scf_drv.compute(new_mol, ao_basis, min_basis)
                            energy_jyp = self.scf_drv.get_scf_energy()

                            # Plus x, plus y
                            coords[i, x] += self.delta_h
                            new_mol = Molecule(labels, coords, units='au')
                            self.scf_drv.compute(new_mol, ao_basis, min_basis)
                            energy_ixp_jyp = self.scf_drv.get_scf_energy()
                            coords[i, x] -= self.delta_h
 
                            # Minus y
                            coords[j, y] -= 2.0 * self.delta_h
                            new_mol = Molecule(labels, coords, units='au')
                            self.scf_drv.compute(new_mol, ao_basis, min_basis)
                            energy_jym = self.scf_drv.get_scf_energy()

                            # Minus x, minus y:
                            coords[i, x] -= self.delta_h
                            new_mol = Molecule(labels, coords, units='au')
                            self.scf_drv.compute(new_mol, ao_basis, min_basis)
                            energy_ixm_jym = self.scf_drv.get_scf_energy()

                            coords[i, x] += self.delta_h
                            coords[j, y] += self.delta_h

                            hessian[i, x, j, y] = (
                        ( energy_ixp_jyp - energy_ixp - energy_jyp 
                        + 2 * energy_0 
                        - energy_ixm - energy_jym + energy_ixm_jym ) 
                                / (2 * self.delta_h**2 )
                                                )
                            hessian[j, y, i, x] = hessian[i, x, j, y]
        # reshaped Hessian as member variable
        self.hessian = hessian.reshape(3*natm, 3*natm)

        # restore scf_drv to initial state
        self.scf_drv.compute(molecule, ao_basis, min_basis)
        self.scf_drv.ostream.state = scf_ostream_state


    def compute_dipole_gradient(self, molecule, ao_basis, perturbed_density):
        """
        Computes the analytical gradient of the dipole moment.

        :param molecule:
            The molecule.
        :param ao_basis:
            The AO basis set.
        :param perturbed_density:
            The perturbed density matrix.
        """

        # Number of atoms and atomic charges
        natm = molecule.number_of_atoms()
        nuclear_charges = molecule.elem_ids_to_numpy()

        density = self.scf_drv.scf_tensors['D_alpha']

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

        # TODO: replace once analytical integral derivatives are available
        dipole_integrals_deriv = (
            self.compute_dipole_integral_derivatives(molecule, ao_basis) )

        # Add the electronic contributions
        dipole_gradient += -2 * (
              np.einsum('mn,caxmn->cax', density, dipole_integrals_deriv)
            + np.einsum('axmn,cmn->cax', perturbed_density, dipole_ints)
                           )

        self.dipole_gradient = dipole_gradient.reshape(3, 3 * natm)


    def compute_dipole_integral_derivatives(self, molecule, ao_basis):
        """
        Computes numerical derivatives of dipole integrals.

        :param molecule:
            The molecule.
        :param ao_basis:
            The AO basis set.

        :return:
            The dipole integral derivatives.
        """

        # atom labels
        labels = molecule.get_labels()

        # number of atoms
        natm = molecule.number_of_atoms()

        # atom coordinates (nx3)
        coords = molecule.get_coordinates()

        # number of atomic orbitals
        nao = self.scf_drv.scf_tensors['D_alpha'].shape[0]

        # Dipole integrals driver
        dipole_drv = ElectricDipoleIntegralsDriver(self.comm)

        # 3 dipole components x No. atoms x 3 atomic coordinates
        # x No. basis x No. basis
        dipole_integrals_gradient = np.zeros((3, natm, 3, nao, nao))

        # smaller delta_h values can be used here
        local_delta_h = 0.01 * self.delta_h

        for i in range(natm):
            for d in range(3):
                coords[i, d] += local_delta_h
                new_mol = Molecule(labels, coords, units='au')

                dipole_mats_p = dipole_drv.compute(new_mol, ao_basis)
                dipole_ints_p = ( dipole_mats_p.x_to_numpy(),
                                  dipole_mats_p.y_to_numpy(),
                                  dipole_mats_p.z_to_numpy())

                coords[i, d] -= 2.0 * local_delta_h
                new_mol = Molecule(labels, coords, units='au')

                dipole_mats_m = dipole_drv.compute(new_mol, ao_basis)
                dipole_ints_m = ( dipole_mats_m.x_to_numpy(),
                                  dipole_mats_m.y_to_numpy(),
                                  dipole_mats_m.z_to_numpy())

                for c in range(3):
                    dipole_integrals_gradient[c, i, d] = (
                        ( dipole_ints_p[c] - dipole_ints_m[c] ) 
                            / (2.0 * local_delta_h) )

        return dipole_integrals_gradient


