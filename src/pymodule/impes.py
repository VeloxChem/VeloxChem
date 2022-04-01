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
          bond angles and dihedral angles).
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
            self.z_matrix = parse_z_matrix(impes_dict['z_matrix'])

        if 'energy' in impes_dict:
            self.energy = float(impes_dict['energy'])

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
        if self.internal_coordinates is None:
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
        self.compute_internal_coordinates()

    # TODO: is it a good idea to do it this way??
    def reset(self):
        """
        Resets some variables to be able to reuse the object.
        """
        self.internal_gradient = None
        self.internal_hessian = None
        self.internal_coordinates_array = None
        self.b_matrix = None
        self.b2_matrix = None
        self.g_minus = None

    def compute_internal_coordinates(self):
        """
        Creats an array with the internal coordinates.

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

    # TODO: return new variable instead of shifting the
    # cartesian coordinates themselves?
    def translate_to_zero(self):
        """Translates the molecule (Cartesian coordinates)
           to (0, 0, 0)
         """

        if self.cartesian_coordinates is None:
            raise ValueError("No Cartesian coordinates found.")
        else:
            natoms = self.cartesian_coordinates.shape[1]
            sum_x = 0.0
            sum_y = 0.0
            sum_z = 0.0

            for i in range(natoms):
                sum_x += self.cartesian_coordinates[i,0]
                sum_y += self.cartesian_coordinates[i,1]
                sum_z += self.cartesian_coordinates[i,2]
            sum_x /= natoms
            sum_y /= natoms
            sum_z /= natoms

            for i in range(natoms):
                self.cartesian_coordinates[i,0] -= sum_x
                self.cartesian_coordinates[i,1] -= sum_y
                self.cartesian_coordinates[i,2] -= sum_z

    # TODO: change to Eckart routine 
    # is it possible to make use of geometric here?
    # save in self or return?
    def rotate_to(self, reference):
        """Rotates the coordinates to align to a reference.

            :param reference: numpy array of Cartesian coordinates
                              to be used as a reference.
        """
        S = np.einsum("ix,iy->xy", reference, self.cartesian_coordinates)
        # TODO inverse and square root of S S.T 
        # (diagonal?, otherwise diagonalize) 
        sqS = np.sqrt(np.linalg.inv(np.matmul(S, S.T)))
        rotation_matrix = np.matmul(sqS, S)

        rotated_coords = np.matmul(rotation_matrix, self.cartesian_coordinates)

        return rotated_coords 

    def write_hdf5(self, fname, label, write_zmat=False):
        """
        Writes the energy, internal coordinates, internal gradient, and
        internal Hessian to a checkpoint file.

        :param fname:
            Name of the checkpoint file to be written.
        :param label:
            A string describing the coordinate.
        :param write_zmat:
            flag to determine if the Z-matrix should be writen to the
            checkpoint file. False by default.
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
            self.compute_internal_coordinates()
            h5f.create_dataset(full_label, data=self.internal_coordinates_array,
                               compression='gzip')

        if write_zmat:
            if self.z_matrix is None:
                return False
            else:
                # internal coordinates
                full_label = label + "_z_matrix"
                zmat = self.get_z_matrix_as_np_array()
                h5f.create_dataset(full_label, zmat, compression='gzip')

        h5f.close()
        return True


    # TODO: complpete
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
        valid_checkpoint = (fname and isinstance(fname, str))

        if not valid_checkpoint:
            raise ValueError("Not a valid checkpoint file.")

        try:
            h5f = h5py.File(fname, 'r')
        except IOError:
            raise ValueError("Could not open checkpoint file.")

        if self.r_inverse:
            label += "_rinv"
        energy_label = label + "_energy"
        gradient_label = label + "_gradient"
        hessian_label = label + "_hessian"
        coords_label = label + "_internal_coordinates"
        self.energy = np.array(h5f.get(energy_label))
        self.internal_gradient = np.array(h5f.get(gradient_label))
        self.internal_hessian = np.array(h5f.get(hessian_label))
        self.internal_coordinates_array = np.array(h5f.get(coords_label))

        h5f.close()


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

        h5f.close() 

class ImpesDriver():
    """
    Implements the potential energy surface interpolation driver

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.

    Instance variables
        - list_of_points: a list of impes coordinates?
        - z_matrix: the Z-matrix of indices defining the internal
          coordinates.
        - r_inv: flag which determines if to use 1/R.
        - impes_coordinate: the new molecular geometry at which the energy
          is to be calculated by interpolation.
        - energy: the energy determined by interpolation.
        - exponent_p: exponent required in the interpolation procedure (both
          simple and Shepard interpolation).
        - exponent_q: exponent required in the Shepard interpolation procedure.
        - confidence_radius: confidence radius required by the Shepard 
          interpolation procedure.
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

        # Z-matrix, energy and coordinate
        self.z_matrix = None
        self.energy = None
        self.impes_coordinate = None

        # r_inv: if to use 1/R instead of R
        self.r_inv = True
        
        # simple or Shepard interpolation
        self.interpolation_type = 'shepard'
        self.exponent_p = None
        self.exponent_q = None
        self.confidence_radius = None

    def update_settings(self, impes_dict=None):
        """
        Updates settings in the ImpesDriver.

        :param impes_dict:
            The input dictionary of settings for IMPES.
        """ 

        if impes_dict is None:
            impes_dict = {}

        if 'z_matrix' in impes_dict:
            self.z_matrix = parse_z_matrix(impes_dict['z_matrix'])

        if 'r_inverse' in impes_dict:
            key = impes_dict['r_inverse'].lower()
            self.r_inverse = True if key in ['yes', 'y'] else False

        if 'interpolation_type' in impes_dict:
            key = impes_dict['interpolation_type'].lower()
            if key in ['simple', 'shepard']:
                self.interpolation_type = key

        if 'exponent_p' in impes_dict:
            # 2p, JCTC 12, 5235 (2016)
            self.exponent_p = 2*float(impes_dict['exponent_p'])
        if 'exponent_q' in impes_dict:
            # 2q, JCTC 12, 5235 (2016)
            self.exponent_q = 2*float(impes_dict['exponent_q'])

        if self.interpolation_type == 'shepard' and self.exponent_q is None:
            self.exponent_q = self.exponent_p / 2.0

        if 'confidence_radius' in impes_dict:
            self.confidence_radius = float(impes_dict['confidence_radius'])

        self.impes_dict = dict(impes_dict)

    def define_impes_coordinate(self, coordinates):
        """Defines the current coordinate based on the molecule object.
           The energy of this molecular geometry is to  be determined by
           interpolation.

            :param coordinates:
                a numpy array of Cartesian coordinates.
        """
        self.impes_coordinate = ImpesCoordinates(self.comm, self.ostream)
        self.impes_coordinate.update_settings(self.impes_dict)
        self.impes_coordinate.cartesian_coordinates = coordinates
        self.impes_coordinate.define_internal_coordinates()
        self.impes_coordinate.compute_internal_coordinates()

    def compute(self, coordinates, fname, labels):
        """Computes the energy by interpolation between pre-defined points.

            :param coordinates:
                a numpy array of Cartesian coordinates.
            :param fname:
                the name of the checkpoint file  containing the pre-defined
                data points.
            :param labels:
                a list of data point labels.
        """

        self.define_impes_coordinate(coordinates)

        if self.interpolation_type == 'simple':
            self.simple_interpolation(fname, labels)
        elif self.interpolation_type == 'shepard':
            self.shepard_interpolation(fname, labels)
        else:
            errtxt = "Unrecognized interpolation type: "
            errtxt += self.interpolation_type
            raise ValueError(errtxt)

    def simple_interpolation(self, fname, labels):
        """Performs a simple interpolation.

        :param fname:
            The name of the checkpoint file  containing the pre-defined
            data points.
        :param labels:
            A list of data point labels.
        """
        n_points = len(labels)
        data_point = ImpesCoordinates(self.comm, self.ostream)
        data_point.update_settings(self.impes_dict)
        sum_weights = 0.0
        potentials = []
        weights = []
        # Determine weights and energy values for interpolation
        for label in labels:
            data_point.read_hdf5(fname, label)
            distance = np.linalg.norm( 
                            self.impes_coordinate.internal_coordinates_array
                          - data_point.internal_coordinates_array )
            weight = 1.0 / (distance**self.exponent_p)
            sum_weights += weight
            potential = self.compute_potential(data_point)
            potentials.append(potential)
            weights.append(weight)

        # Perform the interpolation and save the result in self.energy
        self.energy = 0
        weights /= sum_weights
        for i in range(n_points):
            self.energy += weights[i] * potentials[i]

    def shepard_interpolation(self, fname, labels):
        """Performs a simple interpolation.

        :param fname:
            The name of the checkpoint file  containing the pre-defined
            data points.
        :param labels:
            A list of data point labels.
        """
        n_points = len(labels)
        data_point = ImpesCoordinates(self.comm, self.ostream)
        data_point.update_settings(self.impes_dict)
        sum_weights = 0.0
        potentials = []
        weights = []
        # Determine weights and potential energy surface values for interpolation
        for label in labels:
            data_point.read_hdf5(fname, label)
            distance = np.linalg.norm(
                            self.impes_coordinate.internal_coordinates_array 
                          - data_point.internal_coordinates_array )
            denominator = ( 
                    ( distance / self.confidence_radius )**self.exponent_p
                  + ( distance / self.confidence_radius )**self.exponent_q )
            weight = 1.0 / denominator
            sum_weights += weight
            potential = self.compute_potential(data_point)
            potentials.append(potential)
            weights.append(weight)

        # Perform the interpolation and save the result in self.energy
        self.energy = 0
        weights /= sum_weights
        for i in range(n_points):
            self.energy += weights[i] * potentials[i]

    def compute_potential(self, data_point):
        """Calculates the potential energy surface at self.impes_coordinate
           based on the energy, gradient and Hessian of data_point.

           :param data_point:
                ImpesCoordinates object.
        """

        energy = data_point.energy
        grad = data_point.internal_gradient
        hessian = data_point.internal_hessian
        dist = (   self.impes_coordinate.internal_coordinates_array 
                 - data_point.internal_coordinates_array )

        pes = (   energy + np.matmul(dist.T, grad) 
                + 0.5 * np.linalg.multi_dot([dist.T, hessian, dist]) 
                )

        return pes

    def cartesian_distance(self, data_point):
        """Calculates and returns the cartesian distance between 
           self.coordinates and data_point coordinates.

           :param data_point:
                ImpesCoordinates object
        """

        # First, translate the data point to zero
        # self should also be shifted
        data_point.translate_to_zero()
        self.impes_coordinate.translate_to_zero()

        # Then, determine the rotation matrix which
        # aligns data_point and self.impes_coordinate
        # TODO: change to Eckart routine
        reference_coordinates = self.impes_coordinate.cartesian_coordinates
        data_point.rotate_to(reference_coordinates)

        distance = np.linalg.norm(  reference_coordinates
                                  - data_point.cartesian_coordinates)

        return distance

