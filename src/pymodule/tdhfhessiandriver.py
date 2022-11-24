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
import h5py

from .molecule import Molecule
from .tdhfgradientdriver import TdhfGradientDriver
from .hessiandriver import HessianDriver
from .outputstream import OutputStream
#from .firstorderprop import FirstOrderProperties
from .veloxchemlib import mpi_master

# For PySCF integral derivatives
from .import_from_pyscf import overlap_deriv
from .import_from_pyscf import fock_deriv
from .import_from_pyscf import eri_deriv


class TdhfHessianDriver(HessianDriver):
    """
    Implements the Hessian at the TDHF and CIS levels of theory.

    :param scf_drv:
        The SCF driver.
    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.

    Instance variables
        - scf_drv: The SCF driver.
    """

    def __init__(self, scf_drv, comm=None, ostream=None):
        """
        Initializes gradient driver.
        """

        super().__init__(comm, ostream)

        self.flag = 'RPA Hessian Driver'
        self.scf_drv = scf_drv
        self.state_deriv_index = 0 # calculate 1st excited state by default
        self.tamm_dancoff = False

        self.checkpoint_file = f'{scf_drv._filename}.tdhessian.h5'
        self.save_checkpoint = False

    def update_settings(self, method_dict, rsp_dict, freq_dict=None,
                        orbrsp_dict=None):
        """
        Updates settings in ScfHessianDriver.

        :param method_dict:
            The input dictionary of method settings group.
        :param rsp_dict:
            The input dictionary of response settings group.
        :param freq_dict:
            The input dictionary of Hessian/frequency settings group.
        :param orbrsp_dict:
            The input dictionary of orbital response settings group.
        """

        if freq_dict is None:
            freq_dict = {}
        if orbrsp_dict is None:
            orbrsp_dict = {}

        super().update_settings(method_dict, freq_dict)

        # TODO: set default values in __init__ instead
        if 'state_deriv_index' in freq_dict:
            # user gives '1' for first excited state, but internal index is 0
            self.state_deriv_index = int(freq_dict['state_deriv_index']) - 1

        if 'tamm_dancoff' in rsp_dict:
            key = rsp_dict['tamm_dancoff'].lower()
            if key in ['yes', 'y']:
                self.tamm_dancoff = True
            else:
                self.tamm_dancoff = False

        if self.tamm_dancoff:
            self.flag = 'TDA Hessian Driver'

        if 'save_checkpoint' in freq_dict:
            key = freq_dict['save_checkpoint'].lower()
            self.save_checkpoint = True if key in ['yes', 'y'] else False

        self.rsp_dict = dict(rsp_dict)
        self.orbrsp_dict = dict(orbrsp_dict)


    def compute(self, molecule, ao_basis, rsp_drv, min_basis=None):
        """
        Computes the numerical nuclear Hessian at the RPA or TDA level of theory.

        :param molecule:
            The molecule.
        :param ao_basis:
            The AO basis set.
        :param rsp_drv:
            The RPA or TDA driver.
        :param min_basis:
            The minimal AO basis set.
        """
        self.print_header()

        start_time = tm.time()

        if True: #self.numerical:
            self.compute_numerical(molecule, ao_basis, rsp_drv, min_basis)
        #else:
        #    self.compute_analytical(molecule, ao_basis)


        # save Hessian in a checkpoint file
        if self.save_checkpoint:
            self.write_checkpoint()

        # print Hessian
        if self.do_print_hessian is True:
            self.print_geometry(molecule)
            self.ostream.print_blank()
            self.print_hessian(molecule)

        valstr = '*** Time spent in Hessian calculation: '
        valstr += '{:.2f} sec ***'.format(tm.time() - start_time)
        self.ostream.print_header(valstr)
        self.ostream.print_blank()
        self.ostream.flush()

###    def compute_analytical(self, molecule, ao_basis):
###        """
###        Computes the analytical nuclear Hessian.
###
###        :param molecule:
###            The molecule.
###        :param ao_basis:
###            The AO basis set.
###        """
###        return

    def compute_numerical(self, molecule, ao_basis, rsp_drv, min_basis=None):
        """
        Performs calculation of numerical Hessian.

        :param molecule:
            The molecule.
        :param ao_basis:
            The AO basis set.
        :param rsp_drv:
            The RPA or TDA driver.
        :param min_basis:
            The minimal AO basis set.
        """

        # save the ostream state and set to False
        rsp_ostream_state = rsp_drv.ostream.state
        rsp_drv.ostream.state = False

        # settings dictionary for gradient driver
        grad_dict = dict(self.freq_dict)
        grad_dict['save_checkpoint'] = 'no'

        # required to compute the relaxed dipole moment
        grad_dict['do_first_order_prop'] = 'yes'

        if self.numerical_grad:
            grad_dict['numerical'] = 'yes'
            warn_msg = '*** Warning: Numerical Hessian will be calculated based on numerical gradient.'
            self.ostream.print_header(warn_msg.ljust(56))
            warn_msg = '  This takes a long time and has limited accuracy.'
            self.ostream.print_header(warn_msg.ljust(56))
            self.ostream.print_blank()
            self.ostream.flush()
        else:
            grad_dict['numerical'] = 'no'

        scf_ostream_state = self.scf_drv.ostream.state
        self.scf_drv.ostream.state = False

        # atom labels
        labels = molecule.get_labels()

        # number of atoms
        natm = molecule.number_of_atoms()

        # atom coordinates (nx3)
        coords = molecule.get_coordinates()

        # gradient driver
        grad_drv = TdhfGradientDriver(self.scf_drv, self.scf_drv.comm, self.scf_drv.ostream)
        grad_drv.update_settings(grad_dict, self.rsp_dict, self.orbrsp_dict, self.method_dict)

        # numerical dipole moment gradient (3 dipole components, no. atoms x 3 atom coords)
        self.dipole_gradient = np.zeros((3, 3 * molecule.number_of_atoms()))

        # Hessian in temporary variable
        hessian = np.zeros((natm, 3, natm, 3))

        if not self.do_four_point:
            for i in range(molecule.number_of_atoms()):
                for d in range(3):
                    coords[i, d] += self.delta_h
                    new_mol = Molecule(labels, coords, units='au')
                    self.scf_drv.compute(new_mol, ao_basis, min_basis)
                    scf_tensors = self.scf_drv.scf_tensors
                    rsp_drv._is_converged = False  # only needed for RPA
                    rsp_results = rsp_drv.compute(new_mol, ao_basis,
                                                       scf_tensors)
                    grad_drv.compute(new_mol, ao_basis, rsp_drv, rsp_results)
                    grad_plus = grad_drv.get_gradient()
                    mu_plus = grad_drv.relaxed_dipole_moment


                    coords[i, d] -= 2.0 * self.delta_h
                    new_mol = Molecule(labels, coords, units='au')
                    self.scf_drv.compute(new_mol, ao_basis, min_basis)
                    scf_tensors = self.scf_drv.scf_tensors
                    rsp_drv._is_converged = False  # only needed for RPA
                    rsp_results = rsp_drv.compute(new_mol, ao_basis,
                                                       scf_tensors)
                    grad_drv.compute(new_mol, ao_basis, rsp_drv, rsp_results)
                    grad_minus = grad_drv.get_gradient()
                    mu_minus = grad_drv.relaxed_dipole_moment

                    for c in range(3):
                        self.dipole_gradient[c, 3*i + d] = (mu_plus[c] - mu_minus[c]) / (2.0 * self.delta_h)

                    coords[i, d] += self.delta_h
                    hessian[i, d, :, :] = (grad_plus - grad_minus) / (2.0 * self.delta_h)


        else:
            # Four-point numerical derivative approximation
            # for debugging of analytical Hessian:
            # [ f(x - 2h) - 8 f(x - h) + 8 f(x + h) - f(x + 2h) ] / ( 12h )
            for i in range(molecule.number_of_atoms()):
                for d in range(3):
                    coords[i, d] += self.delta_h
                    new_mol = Molecule(labels, coords, units='au')
                    self.scf_drv.compute(new_mol, ao_basis, min_basis)
                    scf_tensors = self.scf_drv.scf_tensors
                    rsp_drv._is_converged = False  # only needed for RPA
                    rsp_results = rsp_drv.compute(new_mol, ao_basis,
                                                       scf_tensors)
                    grad_drv.compute(new_mol, ao_basis, rsp_drv, rsp_results)
                    grad_plus1 = grad_drv.get_gradient()
                    mu_plus1 = grad_drv.relaxed_dipole_moment

                    coords[i, d] += self.delta_h
                    new_mol = Molecule(labels, coords, units='au')
                    self.scf_drv.compute(new_mol, ao_basis, min_basis)
                    scf_tensors = self.scf_drv.scf_tensors
                    rsp_drv._is_converged = False  # only needed for RPA
                    rsp_results = rsp_drv.compute(new_mol, ao_basis,
                                                       scf_tensors)
                    grad_drv.compute(new_mol, ao_basis, rsp_drv, rsp_results)
                    grad_plus2 = grad_drv.get_gradient()
                    mu_plus2 = grad_drv.relaxed_dipole_moment

                    coords[i, d] -= 3.0 * self.delta_h
                    new_mol = Molecule(labels, coords, units='au')
                    self.scf_drv.compute(new_mol, ao_basis, min_basis)
                    scf_tensors = self.scf_drv.scf_tensors
                    rsp_drv._is_converged = False  # only needed for RPA
                    rsp_results = rsp_drv.compute(new_mol, ao_basis,
                                                       scf_tensors)
                    grad_drv.compute(new_mol, ao_basis, rsp_drv, rsp_results)
                    grad_minus1 = grad_drv.get_gradient()
                    mu_minus1 = grad_drv.relaxed_dipole_moment

                    coords[i, d] -= self.delta_h
                    new_mol = Molecule(labels, coords, units='au')
                    self.scf_drv.compute(new_mol, ao_basis, min_basis)
                    scf_tensors = self.scf_drv.scf_tensors
                    rsp_drv._is_converged = False  # only needed for RPA
                    rsp_results = rsp_drv.compute(new_mol, ao_basis,
                                                       scf_tensors)
                    grad_drv.compute(new_mol, ao_basis, rsp_drv, rsp_results)
                    grad_minus2 = grad_drv.get_gradient()
                    mu_minus2 = grad_drv.relaxed_dipole_moment

                    for c in range(3):
                        self.dipole_gradient[c, 3*i + d] = (mu_minus2[c] - 8.0 * mu_minus1[c]
                                                         + 8.0 * mu_plus1[c] - mu_plus2[c]) / (12.0 * self.delta_h)
                    coords[i, d] += 2.0 * self.delta_h
                    # f'(x) ~ [ f(x - 2h) - 8 f(x - h) + 8 f(x + h) - f(x + 2h) ] / ( 12h )
                    hessian[i, d] = (grad_minus2 - 8.0 * grad_minus1
                                           + 8.0 * grad_plus1 - grad_plus2) / (12.0 * self.delta_h)

        # reshaped Hessian as member variable
        self.hessian = hessian.reshape(3*natm, 3*natm)

        #self.ostream.print_blank()

        self.scf_drv.compute(molecule, ao_basis, min_basis)
        rsp_drv._is_converged = False
        rsp_results = rsp_drv.compute(molecule, ao_basis, self.scf_drv.scf_tensors)
        self.elec_energy = (self.scf_drv.get_scf_energy()
                           + rsp_results['eigenvalues'][self.state_deriv_index])
        self.scf_drv.ostream.state = scf_ostream_state
        rsp_drv.ostream.state = rsp_ostream_state


    def write_checkpoint(self):
        """
        Writes the TDDFT/TDHF Hessian to a checkpoint file.

        """

        if self.rank == mpi_master():
            if self.checkpoint_file is None:
                self.checkpoint_file = f'{self.scf_drv._filename}.tdhessian.h5'

            h5f = h5py.File(self.checkpoint_file, 'a')

            label = "hessian_state_%d" % (self.state_deriv_index + 1)
            h5f.create_dataset(label,
                               data=self.hessian,
                               compression='gzip')

            h5f.close()
