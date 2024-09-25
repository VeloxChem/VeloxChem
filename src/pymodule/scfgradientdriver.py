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
from copy import deepcopy
import numpy as np
import time as tm

from .veloxchemlib import XCFunctional, MolecularGrid
from .veloxchemlib import mpi_master
from .outputstream import OutputStream
from .gradientdriver import GradientDriver
from .moleculargradientdriver import MolecularGradientDriver
from .errorhandler import assert_msg_critical


class ScfGradientDriver(GradientDriver):
    """
    Implements SCF gradient driver.

    :param scf_drv:
        The SCF driver.

    Instance variables
        - flag: The driver flag.
        - delta_h: The displacement for finite difference.
    """

    def __init__(self, scf_drv):
        """
        Initializes gradient driver.
        """

        super().__init__(scf_drv.comm, scf_drv.ostream)

        self.scf_driver = scf_drv
        self.flag = 'SCF Gradient Driver'

        self.numerical = True
        self.delta_h = 0.001

    def compute(self, molecule, basis, scf_results):
        """
        Performs calculation of gradient.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        :param scf_results:
            The dictionary containing converged SCF results.
        """

        start_time = tm.time()
        self.print_header()

        if self.rank == mpi_master():
            nocc = molecule.number_of_alpha_electrons()

            D = scf_results['D_alpha']

            ene_occ = scf_results['E_alpha'][:nocc]
            mo_occ = scf_results['C_alpha'][:, :nocc].copy()
            # TODO: move minus sign into MolecularGradientDriver
            W = -1.0 * np.linalg.multi_dot([mo_occ, np.diag(ene_occ), mo_occ.T])
        else:
            D = None
            W = None

        D = self.comm.bcast(D, root=mpi_master())
        W = self.comm.bcast(W, root=mpi_master())

        molgrad_drv = MolecularGradientDriver()
        # TODO: add exchange_factor and omega
        # TODO: add DFT contribution
        self.gradient = molgrad_drv.mpi_comp_electronic_grad(
            self.comm, molecule, basis, D, W, 0.0, 0.0)

        if self.rank == mpi_master():
            # add nuclear contribution
            self.gradient += self.grad_nuc_contrib(molecule)

        self.gradient = self.comm.allreduce(self.gradient, op=MPI.SUM)

        # print gradient
        self.print_geometry(molecule)
        self.print_gradient(molecule)

        valstr = '*** Time spent in gradient calculation: '
        valstr += '{:.2f} sec ***'.format(tm.time() - start_time)
        self.ostream.print_header(valstr)
        self.ostream.print_blank()
        self.ostream.flush()

    def compute_numerical_gradient(self, molecule, ao_basis, scf_results):
        """
        Performs calculation of gradient.

        :param molecule:
            The molecule.
        :param ao_basis:
            The AO basis set.
        :param scf_results:
            The dictionary containing converged SCF results.
        """

        start_time = tm.time()
        self.print_header()

        self.ostream.mute()
        # Currently, only numerical gradients activated
        self.compute_numerical(molecule, ao_basis, scf_results)
        self.ostream.unmute()

        # print gradient
        self.print_geometry(molecule)
        self.print_gradient(molecule)

        valstr = '*** Time spent in gradient calculation: '
        valstr += '{:.2f} sec ***'.format(tm.time() - start_time)
        self.ostream.print_header(valstr)
        self.ostream.print_blank()
        self.ostream.flush()

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

        self.ostream.mute()
        self.scf_driver.restart = False
        new_scf_results = self.scf_driver.compute(molecule, ao_basis)
        assert_msg_critical(self.scf_driver.is_converged,
                            'ScfGradientDriver: SCF did not converge')
        self.ostream.unmute()

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
