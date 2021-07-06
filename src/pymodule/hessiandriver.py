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

import numpy as np

class HessianDriver:
    """
    Implements the Hessian driver.

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.

    Instance variables
        - hessian: The Hessian.
        - flag: The type of Hessian driver.
        - numerical: Perform numerical Hessian calculation.
    """

    def __init__(self, comm, ostream):
        """
        Initializes Hessian driver.
        """

        self.comm = comm
        self.ostream = ostream

        self.hessian = None
        self.flag = None
        self.numerical = True
        # flag for two-point or four-point approximation
        self.do_four_point = False

    def compute(self, molecule, ao_basis=None, min_basis=None):
        """
        Performs calculation of numerical gradient.

        :param molecule:
            The molecule.
        :param ao_basis:
            The AO basis set.
        :param min_basis:
            The minimal AO basis set.
        """

        return

    def hess_nuc_contrib(self, molecule):
        """
        Calculates the contribution of the nuclear-nuclear repulsion
        to the analytical nuclear Hessian.

        :param molecule:
            The molecule.

        :return:
            The nuclear contribution to the Hessian.
        """
        # number of atoms
        natm = molecule.number_of_atoms()

        # nuclear repulsion energy contribution to Hessian
        nuc_contrib = np.zeros((natm, 3, natm, 3))

        # atom coordinates (nx3)
        coords = molecule.get_coordinates()

        # atomic charges
        nuclear_charges = molecule.elem_ids_to_numpy()

        # loop over all distinct atom pairs and add energy contribution
        for i in range(natm):
            z_a = nuclear_charges[i]
            r_a = coords[i]
            for j in range(natm):
                if i != j:
                    z_b = nuclear_charges[j]
                    r_b = coords[j]
                    r = np.sqrt(np.dot(r_a - r_b, r_a - r_b))
                    for k in range(3):
                        for l in range(3):
                    # off-diagonal parts
                            nuc_contrib[i, k, j, l] = - 3*z_a*z_b*(r_b[k] - r_a[k])*(r_b[l] - r_a[l]) / r**5
                            if k == l:
                                nuc_contrib[i, k, j, l] += z_a * z_b / r**3

                    # add the diagonal contribution
                            nuc_contrib[i, k, i, l] += 3*z_a*z_b*(r_b[k] - r_a[k])*(r_b[l] - r_a[l]) / r**5
                            if k == l:
                                nuc_contrib[i, k, i, l] -= z_a * z_b / r**3


        return nuc_contrib.reshape(3*natm, 3*natm)


    def get_gradient(self):
        """
        Gets the gradient.

        :return:
            The gradient.
        """

        return self.gradient

    def print_geometry(self, molecule):
        """
        Prints the geometry.

        :param molecule:
            The molecule.
        """

        self.ostream.print_block(molecule.get_string())

    def diagonalize_hessian(self, molecule, basis, rsp_drv=None):
        """
        Diagonalizes the Hessian matrix to obtain force constants as eigenvalues
        (and hence vibrational frequencies) and normal modes as eigenvectors.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        :param rsp_drv:
            The RPA or TDA driver (for TdhfHessianDriver).

        :return:
            The vibrational frequencies in atomic units.
        """
        # compute the Hessian if not done already
        if self.hessian is None:
            if rsp_drv is None:
                self.compute(molecule, basis)
            else:
                # TdhfHessianDriver inherits this function
                self.compute(molecule, basis, rsp_drv)

        natm = molecule.number_of_atoms()

        # take the square root of the atomic masses, repeat each three times (xyz components)
        # then form the direct product matrix
        masses_sqrt = np.sqrt(molecule.masses_to_numpy())
        masses_repeat = np.repeat(masses_sqrt, 3)
        masses_matrix = masses_repeat.reshape(-1, 1) * masses_repeat

        # reshape the Hessian as 3Nx3N and mass-weight it
        reshaped_hessian = self.hessian.reshape(3*natm, 3*natm)
        self.mass_weighted_hessian = reshaped_hessian / masses_matrix

        # diagonalize the mass-weighted Hessian
        hessian_eigvals, hessian_eigvecs = np.linalg.eigh(self.mass_weighted_hessian)

        # the first 6 elements should be close to zero (translation & rotation)
        return hessian_eigvals[6:]


    def get_hessian(self):
        """
        Gets the Hessian.

        :return:
            The Hessian.
        """

        return self.hessian

    def print_geometry(self, molecule):
        """
        Prints the geometry.

        :param molecule:
            The molecule.
        """

        self.ostream.print_block(molecule.get_string())

    def print_hessian(self, molecule):
        """
        Prints the Hessian.

        :param molecule:
            The molecule.
        """

        # atom labels
        labels = molecule.get_labels()

        if self.numerical:
            title = 'Numerical '
        else:
            title = 'Analytical '

        title += 'Hessian (Hartree/Bohr**2)'
        self.ostream.print_header(title)
        self.ostream.print_header('-' * (len(title) + 2))
        self.ostream.print_blank()

        natm = molecule.number_of_atoms()

        for k in range(0, natm, 2):

            valstr = 'Coord. '

            coord_dict = {0:'(x)', 1:'(y)', 2:'(z)'}
            end = k + 2
            if k + 2 > natm:
                 end = natm
            for i in range(k,end):
                for di in range(3):
                    valstr += '{:>20s}  '.format('     '+str(i+1)+' '+labels[i]+coord_dict[di]+'     ')

            self.ostream.print_line(valstr)
            self.ostream.print_blank()


            for i in range(natm):
                for di in range(3):
                    valstr = '  {:<4s}'.format(str(i+1)+' '+labels[i]+coord_dict[di])
                    for j in range(k,end):
                        #print("j = %d" % j)
                        for dj in range(3):
                            valstr += '{:22.12f}'.format(self.hessian[i*3+di,j*3+dj])
                    self.ostream.print_line(valstr)
            self.ostream.print_blank()
            self.ostream.print_blank()
        self.ostream.print_blank()
        self.ostream.flush()


    def print_header(self):
        """
        Prints Hessian calculation setup details to output stream.
        """

        self.ostream.print_blank()
        self.ostream.print_header(self.flag)
        self.ostream.print_header((len(self.flag) + 2) * '=')
        self.ostream.print_blank()
        self.ostream.flush()



