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
from .outputstream import OutputStream
import numpy as np
import sys
import h5py
import geometric

# TODO: add routines to read energy, gradient, Hessian in internal coordinates.
class ImpesCoordinates:
    """
    Implements routines for the coordinates required for potential
    energy surface interpolation.

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.

    Instance variables
        - hessian: The Cartesian Hessian in Hartree per Bohr**2.
        - gradient: The Cartesian gradient in Hartree per Bohr
        - energy: the total energy of the groound or excited state in Hartree
        - z_matrix: The Z-matrix (a list of coordinates which define bonds
          bond angles and dihedral angles.
        - internal_hessian: The Hessian in internal coordinates.
        - internal_gradient: The gradient in internal coordinates.
        - b_matrix: The Wilson B-matrix.
        - b2_matrix: The matrix of second order deriatives of 
          internal coordinates with respect to Cartesian coordinates
          with respect to Cartesian coordinates.
        - r_inverse: Flag for using the inverse bond length instead of the 
          bond length.
        - internal_coordinates: a list of geomeTRIC internal coordinates.
    """

    def __init__(self, comm=None, ostream=None):
        """
        Initializes the ImpesCoordinates object.
        """

        if comm is None:
            comm = MPI.COMM_WORLD

        if ostream is None:
            ostream = OutputStream(sys.stdout)

        # MPI information
        self.comm = comm
        self.rank = self.comm.Get_rank()
        self.nodes = self.comm.Get_size()
        self.ostream = ostream


        self.hessian = None
        self.gradient = None
        self.z_matrix = None
        self.internal_hessian = None
        self.internal_gradient = None
        self.b_matrix = None
        self.b2_matrix = None
        self.r_inverse = True
        self.internal_coordinates = None
        self.internal_coordinates_array = None
        self.cartesian_coordinates = None
        # TODO: add timing and profiling

    def update_settings(self, impes_dict=None):
        """
        Updates settings in the ImpesDriver.

        :param impes_dict:
            The input dictionary of settings for IMPES.
        """ 

        if impes_dict is None:
            impes_dict = {}

        if 'r_inverse' in impes_dict:
            key = impes_dict['r_inverse'].lower()
            self.r_inverse = True if key in ['yes', 'y'] else False

        if 'z_matrix' in impes_dict:
            self.z_matrix = self.parse_z_matrix(impes_dict['z_matrix'])

        if 'energy' in impes_dict:
            self.energy = float(impes_dict['energy'])

    @staticmethod
    def parse_z_matrix(input_z_matrix):
        """
        Parses the input Z-matrix to create a list of indices which define
        bonds, bond angles and dihedrals.
        
        :param input_z_matrix:
            The string containing the Z-matrix.

        :return:
            a list of bond, bond angle, and dihedral atomic indices.
        """

        z_matrix = []
        for q_str in input_z_matrix.split(';'):
            q = []
            for ind in q_str.split(','):
                q.append(int(ind))
            z_matrix.append(q)

        return z_matrix

    def define_internal_coordinates(self, z_matrix=None):
        """
        Constructs the list of internal coordinates based on the Z-matrix.

        :param z_matrix:
            The list of bond, angle, and dihedral indices which define the 
            Z-matrix.
        """

        if z_matrix is None and self.z_matrix is None:
            raise ValueError("No Z-matrix defined.")
        elif z_matrix is not None:
            self.z_matrix = z_matrix

        if self.internal_coordinates is None:
            self.internal_coordinates = []

        for z in self.z_matrix:
            if len(z) == 2:
                a = z[0]
                b = z[1]
                q = geometric.internal.Distance(a, b)
            elif len(z) == 3:
                a = z[0]
                b = z[1]
                c = z[2]
                q = geometric.internal.Angle(a, b, c)
            else:
                a = z[0]
                b = z[1]
                c = z[2]
                d = z[3]
                q = geometric.internal.Dihedral(a, b, c, d)
            self.internal_coordinates.append(q)

    def calculate_b_matrix(self):
        """
        Calculates the Wilson B-matrix of first order derivatives of the 
        internal coordinates with respect to Cartesian coordinates.
        """

        if self.internal_coordinates is None:
            raise ValueError("No internal coordinates are defined.")

        if self.cartesian_coordinates is None:
            raise ValueError("No cartesian coordinates are defined.")

        n_atoms = self.cartesian_coordinates.shape[0]
        coords = self.cartesian_coordinates.reshape((3*n_atoms))

        if self.b_matrix is None:
            self.b_matrix = np.zeros((len(self.z_matrix), 3*n_atoms))

        for i in range(len(self.z_matrix)):
            z = self.z_matrix[i]
            q = self.internal_coordinates[i]
            derivative = q.derivative(coords)
            for a in z:
                self.b_matrix[i,3*a:3*a+3] = derivative[a]

    def calculate_b2_matrix(self):
        """
        Calculates the B2-matrix of second order derivatives of the 
        internal coordinates with respect to Cartesian coordinates.
        """

        if self.internal_coordinates is None:
            raise ValueError("No internal coordinates are defined.")
        if self.cartesian_coordinates is None:
            raise ValueError("No cartesian coordinates are defined.")

        n_atoms = self.cartesian_coordinates.shape[0]
        coords = self.cartesian_coordinates.reshape((3*n_atoms))

        if self.b2_matrix is None:
            self.b2_matrix = np.zeros((len(self.z_matrix), 3*n_atoms, 3*n_atoms))

        for i in range(len(self.z_matrix)):
            z = self.z_matrix[i]
            q = self.internal_coordinates[i]
            second_derivative = q.second_derivative(coords)
            for m in z:
                for n in z:
                    self.b2_matrix[i,3*m:3*m+3,3*n:3*n+3] = (
                                        second_derivative[m, :, n, :])

    def transform_gradient_to_internal_coordinates(self, tol=1e-7):
        """
        Transforms the gradient from Cartesian to internal coordinates.

        :param tol:
            The tolerance of comparison for inversion.
        """
        if self.b_matrix is None:
            self.calculate_b_matrix()

        n_atoms = self.cartesian_coordinates.shape[0]
        g_matrix = np.matmul(self.b_matrix, self.b_matrix.T)
        gval, gvec = np.linalg.eig(g_matrix)

        # TODO: is there a better way of selecting out the non-zero eigenvalues
        # and computing g_minus?
        gval_inverse = []
        for g in gval:
            if abs(g) > tol:
                gval_inverse.append(1.0/g)
            else:
                gval_inverse.append(0.0)

        gval_inverse_matrix = np.diag(np.array(gval_inverse))
        self.g_minus_matrix = (
                np.linalg.multi_dot([gvec, gval_inverse_matrix, gvec.T]) )
        gradient = self.gradient.reshape((3*n_atoms))
        self.internal_gradient = (np.linalg.multi_dot(
                    [self.g_minus_matrix, self.b_matrix, gradient]) )

    def transform_hessian_to_internal_coordinates(self):
        """
        Transforms the Hessian from Cartesian to internal coordinates.

        """

        if self.internal_gradient is None:
            self.transform_gradient_to_internal_coordinates()
        if self.b2_matrix is None:
            self.calculate_b2_matrix()

        b2_gradient = np.einsum("qxy,q->xy", self.b2_matrix,
                                self.internal_gradient)
        self.internal_hessian = (
                np.linalg.multi_dot([self.g_minus_matrix, self.b_matrix,
                                     self.hessian - b2_gradient,
                                     self.b_matrix.T, self.g_minus_matrix.T])
                                )

    def symmetrize_hessian(self):
        """
        Symmetrizes the Hessian in internal coordinates.
        """
        n = self.internal_hessian.shape[0]
        for i in range(n):
            for j in range(i+1, n):
                avg = 0.5 * (  self.internal_hessian[i,j]
                             + self.internal_hessian[j,i]
                                )
                self.internal_hessian[i,j] = avg
                self.internal_hessian[j,i] = avg

    def transform_to_r_inverse(self):
        """
        Transforms the gradient and Hessian from using the bond length R, to 
        using the inverse bond length 1/R.

        """
        if self.internal_coordinates is None:
            raise ValueError("No internal coordinates are defined.")
        if self.cartesian_coordinates is None:
            raise ValueError("No cartesian coordinates are defined.")

        n_atoms = self.cartesian_coordinates.shape[0]
        coords = self.cartesian_coordinates.reshape((3*n_atoms))

        i = 0
        for z in self.z_matrix:
            j = 0
            if len(z) == 2:
                r = self.internal_coordinates[i].value(coords)
                self.internal_gradient[i] *= -r**2
            for t in self.z_matrix:
                if len(z) == 2:
                    r = self.internal_coordinates[i].value(coords)
                    self.internal_hessian[i, j] *= -r**2
                if len(t) == 2:
                    r = self.internal_coordinates[j].value(coords)
                    self.internal_hessian[i, j] *= -r**2
                if i == j and len(z) == 2:
                    r = self.internal_coordinates[j].value(coords)
                    self.internal_hessian[i, j] += (
                                        - 2.0 * r * self.internal_gradient[i]
                                        )
                j += 1
            i += 1

    def transform_gradient_and_hessian(self, z_matrix=None):
        """
        Performs all steps required to transform from Cartesian coordinates
        to the internal coordinates defined in the Z-matrix.

        :param z_matrix:
            The Z-matrix - a list of indices which define the bonds, angles,
            and dihedrals of interest.
        """

        # Create the internal coordinates through geomeTRIC
        self.define_internal_coordinates(z_matrix)

        # Transform the gradient and Hessian to internal coordinates
        self.transform_gradient_to_internal_coordinates()
        self.transform_hessian_to_internal_coordinates()

        # Symmetrize the Hessian matrix
        self.symmetrize_hessian()

        # Transform the gradient and Hessian to use 1/R
        if self.r_inverse:
            self.transform_to_r_inverse()

        # Save the values of the internal coordinates as a numpy array
        self.get_internal_coordinates()

    def get_internal_coordinates(self):
        """
        Returns an array with the internal coordinates.

        """
        if self.internal_coordinates is None:
            raise ValueError("No internal coordinates are defined.")
        if self.cartesian_coordinates is None:
            raise ValueError("No cartesian coordinates are defined.")

        n_atoms = self.cartesian_coordinates.shape[0]
        coords = self.cartesian_coordinates.reshape((3*n_atoms))

        int_coords = []
        for q in self.internal_coordinates:
            if isinstance(q, geometric.internal.Distance) and self.r_inverse:
                int_coords.append(1.0/q.value(coords))
            else:
                int_coords.append(q.value(coords))

        self.internal_coordinates_array = np.array(int_coords)

    def get_z_matrix_as_np_array(self):
        """
        Returns an array with the Z-matrix
        """
        if self.z_matrix is None:
            raise ValueError("No Z-matrix is defined.")

        return np.array(self.z_matrix)

    def write_hdf5(self, fname, label):
        """
        Writes the energy, internal coordinates, internal gradient, and
        internal Hessian to a checkpoint file.

        :param fname:
            Name of the checkpoint file to be written.
        :param label:
            A string describing the coordinate.
        """
        valid_checkpoint = (fname and isinstance(fname, str))

        if not valid_checkpoint:
            return False
        try:
            h5f = h5py.File(fname, 'a')
        except IOError:
            return False
 
        if self.r_inverse:
            label += "_rinv"

        if self.energy is None:
            return False
        else:
            # energy
            full_label = label + "_energy"
            h5f.create_dataset(full_label, data=self.energy)
        
        if self.internal_gradient is None:
            return False
        else:
            # gradient in internal coordinates
            full_label = label + "_gradient"
            h5f.create_dataset(full_label, data=self.internal_gradient,
                               compression='gzip')

        if self.internal_hessian is None:
            return False
        else:
            # Hessian in internal coordinates
            full_label = label + "_hessian"
            h5f.create_dataset(full_label, data=self.internal_hessian,
                               compression='gzip')

        if self.internal_coordinates is None:
            return False
        else:
            # internal coordinates
            full_label = label + "_internal_coordinates"
            self.get_internal_coordinates()
            h5f.create_dataset(full_label, data=self.internal_coordinates_array,
                               compression='gzip')

        ##if self.z_matrix is None:
        ##    return False
        ##else:
        ##    # internal coordinates
        ##    full_label = label + "_z_matrix"
        ##    zmat = self.get_z_matrix_as_np_array()
        ##    h5f.create_dataset(full_label, zmat, compression='gzip')

        h5f.close()
        return True


    # TODO: finish code/remove?
    def read_hdf5(self, fname, label):
        """
        Reads the energy, internal coordinates, gradient, and Hessian from
        a checkpoint file,

        :param fname:
            Name of the checkpoint file.
            The file must exist before calling this routine.
        :param label:
            The label of the selected coordinates.
        """
        pass

    def read_cartesian_coords_energy_gradient_hessian(self, fname, label):
        """
        Reads the energy, Cartesian gradient, and Hessian from checkpoint file.

        :param fname:
            Name of the checkpoint file.
            The file must exist before calling this routine.
        :param label:
            The label of the selected coordinates.
        """

        valid_checkpoint = (fname and isinstance(fname, str))

        if not valid_checkpoint:
            raise ValueError("Not a valid checkpoint file.")

        try:
            h5f = h5py.File(fname, 'r')
        except IOError:
            raise ValueError("Could not open checkpoint file.")

        energy_label = label + "_energy"
        gradient_label = label + "_gradient"
        hessian_label = label + "_hessian"
        coords_label = label + "_coordinates"
        self.energy = np.array(h5f.get(energy_label))
        self.gradient = np.array(h5f.get(gradient_label))
        self.hessian = np.array(h5f.get(hessian_label))
        self.cartesian_coordinates = np.array(h5f.get(coords_label)) 

class ImpesDriver():
    """
    Implements the potential energy surface interpolation driver

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.

    Instance variables
        - list_of_points: a list of impes coordinates.
    """

    def __init__(self, comm=None, ostream=None):
        """
        Initializes the IMPES driver.
        """

        if comm is None:
            comm = MPI.COMM_WORLD

        if ostream is None:
            ostream = OutputStream(sys.stdout)

        # MPI information
        self.comm = comm
        self.rank = self.comm.Get_rank()
        self.nodes = self.comm.Get_size()
        self.ostream = ostream

        # TODO: add timing and profiling
