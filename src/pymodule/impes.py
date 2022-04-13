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

au_in_fs = 2.418884326585747e-2
fs_in_au = 1.0 / au_in_fs
amu_in_au = 1822.888486217198 

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

def parse_labels(input_labels):
    """
    Parses the input for labels of interpolation data points
    and returns them as a list.
    
    :param input_z_matrix:
        The string containing the Z-matrix.

    :return:
        a list of data point labels.
    """

    labels = []
    for label_str in input_labels.split('\n'):
        label = label_str.replace(" ","")
        labels.append(label)

    return labels


# TODO: MPI-parallelization!
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
            self.b2_matrix = np.zeros(
                    (len(self.z_matrix), 3*n_atoms, 3*n_atoms))

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

    def transform_to_r(self):
        """
        Returns the internal gradient transformed from 
        using the inverse bond length 1/R to using R.
        Also transforms distance * Hessian from 1/R to R.
        (required to transform back to the Cartesian gradient).

        :param dist_hessian:
            result of distance.T Hessian matrix-matrix multiplication 

        """
        if self.internal_coordinates is None:
            self.define_internal_coordinates()
        if self.cartesian_coordinates is None:
            raise ValueError("No cartesian coordinates are defined.")
        if not self.r_inverse:
            raise ValueError("Internal Gradient is already in R.")

        n_atoms = self.cartesian_coordinates.shape[0]
        coords = self.cartesian_coordinates.reshape((3*n_atoms))
        gradient_in_r = self.internal_gradient.copy()
        hessian_in_r = self.internal_hessian.copy()
        i = 0
        for z in self.z_matrix:
            j = 0
            if len(z) == 2:
                r = self.internal_coordinates[i].value(coords)
                gradient_in_r[i] /= -r**2
            for t in self.z_matrix:
                # TODO: *** TEST FOR BIG MOLECULES ***
                if i == j and len(z) == 2:
                    r = self.internal_coordinates[j].value(coords)
                    hessian_in_r[i, j] -= (
                                        - 2.0 * r * self.internal_gradient[i]
                                        )
                if len(z) == 2:
                    r = self.internal_coordinates[i].value(coords)
                    hessian_in_r[i, j] /= -r**2
                if len(t) == 2:
                    r = self.internal_coordinates[j].value(coords)
                    hessian_in_r[i, j] /= -r**2

                j += 1
            i += 1
        return gradient_in_r, hessian_in_r

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

    def return_internal_coordinates(self, r_inverse=True):
        """
        Returns a numpy array with the internal coordinates.

        :param r_inverse:
            flag to decide if R or 1/R is used, where R is the bond length.
        """
        if self.internal_coordinates is None:
            self.define_internal_coordinates()
        if self.cartesian_coordinates is None:
            raise ValueError("No cartesian coordinates are defined.")

        n_atoms = self.cartesian_coordinates.shape[0]
        coords = self.cartesian_coordinates.reshape((3*n_atoms))

        int_coords = []
        for q in self.internal_coordinates:
            if isinstance(q, geometric.internal.Distance) and self.r_inverse:
                if r_inverse:
                    int_coords.append(1.0/q.value(coords))
                else:
                    int_coords.append(q.value(coords))
            else:
                int_coords.append(q.value(coords))

        return np.array(int_coords)

    def get_z_matrix_as_np_array(self):
        """
        Returns an array with the Z-matrix
        """
        if self.z_matrix is None:
            raise ValueError("No Z-matrix is defined.")

        return np.array(self.z_matrix)

    def translate_to_zero(self):
        """Retruns an array with the Cartesian coordinates
           translated to (0, 0, 0)
         """

        if self.cartesian_coordinates is None:
            raise ValueError("No Cartesian coordinates found.")
        else:
            natoms = self.cartesian_coordinates.shape[0]
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

            translated_coordinates = self.cartesian_coordinates.copy()
            for i in range(natoms):
                translated_coordinates[i,0] -= sum_x
                translated_coordinates[i,1] -= sum_y
                translated_coordinates[i,2] -= sum_z
            return translated_coordinates

    def rotate_to(self, reference):
        """Rotates the coordinates to align to a reference.

            :param reference: numpy array of target Cartesian coordinates
                              to be used as a reference.
        """
        N = reference.shape[0]
        rotation_matrix = geometric.rotate.get_rot(self.cartesian_coordinates,
                                                   reference)

        rotated_coords = np.dot(rotation_matrix, self.cartesian_coordinates.T).T
        self.cartesian_coordinates = rotated_coords
        rmsd = np.sqrt(np.sum((rotated_coords-reference)**2)/N)
        return rotation_matrix, rotated_coords, rmsd

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

        if self.cartesian_coordinates is None:
            return False
        else:
            # internal coordinates
            full_label = label + "_cartesian_coordinates"
            h5f.create_dataset(full_label, data=self.cartesian_coordinates,
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


    # TODO: read more variables in?
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
        cart_coords_label = label + "_cartesian_coordinates"
        # TODO
        self.energy = np.array(h5f.get(energy_label))
        self.internal_gradient = np.array(h5f.get(gradient_label))
        self.internal_hessian = np.array(h5f.get(hessian_label))
        self.internal_coordinates_array = np.array(h5f.get(coords_label))
        self.cartesian_coordinates = np.array(h5f.get(cart_coords_label))

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

# TODO: MPI-parallelization!
class ImpesDriver():
    """
    Implements the potential energy surface interpolation driver

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.

    Instance variables
        - z_matrix: the Z-matrix of indices defining the internal
          coordinates.
        - r_inverse: flag which determines if to use 1/R.
        - impes_coordinate: the new molecular geometry at which the energy
          is to be calculated by interpolation.
        - energy: the energy determined by interpolation.
        - gradient: the Cartesian gradient determined by interpolation.
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
        self.gradient = None
        self.impes_coordinate = None

        # r_inv: if to use 1/R instead of R
        self.r_inverse = True
        
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
        self.impes_coordinate.calculate_b_matrix() # required for gradient
        self.impes_coordinate.compute_internal_coordinates()

    def compute(self, coordinates, fname, labels=None):
        """Computes the energy by interpolation between pre-defined points.

            :param coordinates:
                a numpy array of Cartesian coordinates.
            :param fname:
                the name of the checkpoint file  containing the pre-defined
                data points.
            :param labels:
                a list of data point labels (optional).
                if this parameter is None, all datapoints from fname will be
                used.
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

    # TODO: double check weight gradients!
    def simple_interpolation(self, fname, labels):
        """Performs a simple interpolation.

        :param fname:
            The name of the checkpoint file  containing the pre-defined
            data points.
        :param labels:
            A list of data point labels. If labels is None, it will be read
            from fname.
        """

        natms = self.impes_coordinate.cartesian_coordinates.shape[0]
        data_point = ImpesCoordinates(self.comm, self.ostream)
        data_point.update_settings(self.impes_dict)
        sum_weights = 0.0
        potentials = []
        gradients = []
        weights = []
        weight_gradients = []
        sum_weight_gradients = np.zeros((natms, 3))

        if labels is None:
            labels = self.read_labels(fname)
        n_points = len(labels)

        # Determine weights, weight gradients,
        # energy values, and energy gradients for interpolation
        for label in labels:
            data_point.read_hdf5(fname, label)
            distance, weight_gradient = self.cartesian_distance(data_point) 
            #dist =     np.linalg.norm( 
            #                self.impes_coordinate.internal_coordinates_array
            #              - data_point.internal_coordinates_array )
            weight = 1.0 / (distance**self.exponent_p)
            sum_weights += weight
            sum_weight_gradients += weight_gradient
            potential = self.compute_potential(data_point)
            gradient = self.compute_gradient(data_point)
            gradients.append(gradient)
            potentials.append(potential)
            weights.append(weight)
            weight_gradients.append(weight_gradient)

        # Perform the interpolation and save the result in self.energy
        # and self.gradient
        self.energy = 0
        self.gradient = np.zeros((natms, 3))
        weights /= sum_weights
        for i in range(n_points):
            self.energy += weights[i] * potentials[i]
            self.gradient += ( weights[i] * gradients[i]
             + potentials[i] * weight_gradients[i] / sum_weights
             - potentials[i] * weights[i] * sum_weight_gradients / sum_weights 
                                )

    def shepard_interpolation(self, fname, labels):
        """Performs a simple interpolation.

        :param fname:
            The name of the checkpoint file  containing the pre-defined
            data points.
        :param labels:
            A list of data point labels.
        """
        natms = self.impes_coordinate.cartesian_coordinates.shape[0]
        data_point = ImpesCoordinates(self.comm, self.ostream)
        data_point.update_settings(self.impes_dict)
        sum_weights = 0.0
        potentials = []
        gradients = []
        weights = []
        weight_gradients = []
        sum_weight_gradients = np.zeros((natms, 3))
        if labels is None:
            labels = self.read_labels(fname)
        n_points = len(labels)

        # Determine weights, weight gradients,
        # energy values, and energy gradients for interpolation
        for label in labels:
            data_point.read_hdf5(fname, label)
            distance, weight_gradient = self.cartesian_distance(data_point)
                       # np.linalg.norm(
                       #     self.impes_coordinate.internal_coordinates_array 
                       #   - data_point.internal_coordinates_array )
            denominator = ( 
                    ( distance / self.confidence_radius )**self.exponent_p
                  + ( distance / self.confidence_radius )**self.exponent_q )
            weight = 1.0 / denominator
            sum_weights += weight
            sum_weight_gradients += weight_gradient
            potential = self.compute_potential(data_point)
            gradient = self.compute_gradient(data_point)
            gradients.append(gradient)
            potentials.append(potential)
            weights.append(weight)
            weight_gradients.append(weight_gradient)

        # Perform the interpolation and save the result
        # in self.energy, self.gradient
        self.energy = 0
        self.gradient = np.zeros((natms, 3))
        weights /= sum_weights
        for i in range(n_points):
            self.energy += weights[i] * potentials[i]
            self.gradient += ( weights[i] * gradients[i]
             + potentials[i] * weight_gradients[i] / sum_weights
             - potentials[i] * weights[i] * sum_weight_gradients / sum_weights 
                                )

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

    def compute_gradient(self, data_point):
        """Calculates part of the cartesian gradient 
           for self.impes_coordinate based on the gradient
           and Hessian of data_point.
           J ( g_a + delta_a H_a )
           eq. (4) JCTC 12, 5235-5246

           :param data_point:
                ImpesCoordinates object.
        """

        natm = data_point.cartesian_coordinates.shape[0]
        energy = data_point.energy
        grad = data_point.internal_gradient
        hessian = data_point.internal_hessian
        dist = (
            self.impes_coordinate.return_internal_coordinates(r_inverse=False) 
          - data_point.return_internal_coordinates(r_inverse=False)
                )

        if data_point.b_matrix is None:
            data_point.calculate_b_matrix()
        b_matrix = data_point.b_matrix #self.impes_coordinate.b_matrix

        if self.r_inverse:
            # trnasform gradient and hessian back to r.
            # this is required to be able to then transform
            # to the Cartesian gradient.

            grad, hessian = data_point.transform_to_r()

        dist_hessian = np.matmul(dist.T, hessian)

        gradient = (  np.matmul(b_matrix.T, grad)
                    + np.matmul(b_matrix.T, dist_hessian)
                    ).reshape(natm, 3)


        return gradient

    def simple_weight_gradient(self, distance_vector, distance):
        """ Returns the derivative of an unormalized simple interpolation
            weight with respect to the Cartesian coordinates in
            self.impes_coordinate

            :param distance_vector:
                The Cartesian distance vector between
                the current data_point and self.impes_coordinate.
            :param distance:
                The norm of the distance vector * sqrt(N), N number of atoms.
        """
        weight_gradient = ( - self.exponent_p * distance_vector /
                            distance**( self.exponent_p + 1 ) )
        # TODO double check +1 / +2 

        return weight_gradient
    
    def shepard_weight_gradient(self, distance_vector, distance):
        """ Returns the derivative of an unormalized Shepard interpolation
            weight with respect to the Cartesian coordinates in
            self.impes_coordinate

            :param distance_vector:
                The Cartesian distance vector between
                the current data_point and self.impes_coordinate.
            :param distance:
                The norm of the distance vector * sqrt(N), N number of atoms.
        """

        denominator = ( 
                    ( distance / self.confidence_radius )**self.exponent_p
                  + ( distance / self.confidence_radius )**self.exponent_q )
        derivative_p = ( self.exponent_p 
                * ( distance / self.confidence_radius )**( self.exponent_p - 1 )
                * distance_vector / self.confidence_radius
                        )
        derivative_q = ( self.exponent_q 
                * ( distance / self.confidence_radius )**( self.exponent_q - 1 )
                * distance_vector / self.confidence_radius
                        )

        weight_gradient = ( - 1.0 / denominator**2 * 
                                ( derivative_p + derivative_q )
                            )

        return weight_gradient

    def cartesian_distance(self, data_point):
        """Calculates and returns the cartesian distance between 
           self.coordinates and data_point coordinates.
           Besides the distance, it also returns the weight gradient,
           which requires the distance vector to be computed.

           :param data_point:
                ImpesCoordinates object
        """

        # First, translate the cartesian coordinates to zero
        target_coordinates = data_point.translate_to_zero()
        reference_coordinates = self.impes_coordinate.translate_to_zero()

        # Then, determine the rotation matrix which
        # aligns data_point (target_coordinates)
        # to self.impes_coordinate (reference_coordinates)

        rotation_matrix = geometric.rotate.get_rot(target_coordinates,
                                                   reference_coordinates)

        # Rotate the data point
        rotated_coordinates = np.dot(rotation_matrix, target_coordinates.T).T

        # Calculate the Cartesian distance
        N = reference_coordinates.shape[0] # number of atoms
        distance = np.sqrt(N) * np.linalg.norm(  reference_coordinates
                                  - rotated_coordinates)

        # Calculate the gradient of the interpolation weights
        # (required for energy gradient interpolation)
        distance_vector = reference_coordinates - rotated_coordinates

        if self.interpolation_type == 'shepard':
            weight_gradient = self.shepard_weight_gradient(distance_vector,
                                                           distance)
        elif self.interpolation_type == 'simple':
            weight_gradient = self.simple_weight_gradient(distance_vector,
                                                          distance)
        else:
            errtxt = "Unrecognized interpolation type: "
            errtxt += self.interpolation_type
            raise ValueError(errtxt)

        return distance, weight_gradient

    def read_labels(self, fname):
        """
        Read data point labels from checkpoint file.

        :param fname:
            Name of the checkpoint file.
        """

        h5f = h5py.File(fname, 'r')
        key = 'labels'
        labels = []
        for label_bytes in list(h5f.get(key)):
            label = label_bytes.decode()
            labels.append(label)
        h5f.close()

        return labels

class ImpesDynamicsDriver():
    """
    Simple interpolated dynamics driver based on Verlet integration.

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.

    Instance variables
        - impes_driver: an ImpesDriver object that provides the gradient.
        - coordinates: array of coordinates.
        - velocities: array of velocities
        - time_step: the time step to be used (in fs).
        - duration: the total time to run the dynamics (in fs).
        - current_step: time step at which dynamics has arrived-
        - name: name of the checkpoint file where the data points for
                interpolation are stored.
        - labels: labels of the data points to be used for interpolation.
    """

    def __init__(self, comm=None, ostream=None):
        """
            Initializes the IMPES dynamics driver.
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

        # ImpesDriver, coordinates, velocities
        self.impes_driver = ImpesDriver()
        self.coordinates = None
        self.velocities = None

        # time step and total duration
        self.time_step = 1 # fs = 41 a.u.  
        self.duration = 10 # fs = 410 a.u.
        self.current_step = 0

        # name of the checkpoint file with data points and labels
        self.name = None
        self.labels = None


    def update_settings(self, impes_dict=None):
        """
        Updates settings in the ImpesDynamicsDriver.

        :param impes_dict:
            The input dictionary of settings for IMPES.
        """ 
        if impes_dict is None:
            impes_dict = {}

        self.impes_driver.update_settings(impes_dict)

        if 'time_step' in impes_dict:
            self.time_step = float(impes_dict['time_step'])
        if 'duration' in impes_dict:
            self.duration = float(impes_dict['duration'])
        if 'name' in impes_dict:
            self.name = impes_dict['name']
        if 'labels' in impes_dict:
            self.labels = parse_labels(impes_dict['labels'])

    def run_dynamics_step(self, molecule, iteration=0):
        """
        Runs one step of the interpolated dynamics to update
        the atomic coordinates and velocities based on 
        Verlet integration.

        :param molecule:
            The molecule.
        """

        # atomic masses
        masses = molecule.masses_to_numpy() * amu_in_au # transform to au

        # Determine the atomic forces and accelerations
        forces = - self.impes_driver.gradient
        acc = np.einsum('ix,i->ix', forces, 1/masses) 

        # Calculate velocities at t + 1/2 dt
        dt = self.time_step * fs_in_au # transform fs to a.u. 
        velocities = self.velocities[iteration] + 0.5 * acc * dt
        coordinates = self.coordinates[iteration] + velocities * dt

        # calculate energy and gradient for the new coordinates
        # TODO: save potential and kinetic energies too? 
        self.impes_driver.compute(coordinates, self.name, self.labels)

        # recalculate the accelerations
        forces = - self.impes_driver.gradient
        acc_dt = np.einsum('ix,i->ix', forces, 1/masses)

        # calculate the velocity at t + dt
        velocities = self.velocities[iteration] + 0.5 * ( acc + acc_dt ) * dt

        # calculate kinetic energy:
        velocities_squared = np.einsum('ix,ix->ix', velocities, velocities)
        kinetic_energy = np.einsum('i,ix->', masses, 0.5 * velocities_squared)
        total_energy = kinetic_energy + self.impes_driver.energy

        # save results
        self.velocities.append(velocities)
        self.coordinates.append(coordinates)
        self.potential_energies.append(self.impes_driver.energy)
        self.kinetic_energies.append(kinetic_energy)
        self.total_energies.append(total_energy)

        self.print_iteration(iteration)

        # update current step
        self.current_step += self.time_step


    def compute(self, molecule):
        """
        Runs the interpolation dynamics.

        :param molecule:
            The molecule.
        """

        # Initialize coordinates and velocities
        starting_coordinates = molecule.get_coordinates()
        self.coordinates = [ starting_coordinates ]
        self.velocities = [np.zeros_like(starting_coordinates)]

        # Calculate energy and gradient for starting coordinates
        self.impes_driver.compute(starting_coordinates, self.name, self.labels)
        self.potential_energies = [self.impes_driver.energy]
        self.kinetic_energies = [0.0]
        self.total_energies = [self.impes_driver.energy]

        self.print_header()

        # Start dynamics
        iteration = 0
        while self.duration >= self.current_step:
            self.run_dynamics_step(molecule, iteration)
            iteration += 1

    def print_header(self):
        """
        Prints information about the ImpesDynamicsDriver.
        """
        self.ostream.print_blank()

        title = 'Interpolated Dynamics Driver'
        self.ostream.print_header(title)
        self.ostream.print_header('=' * (len(title) + 2))
        self.ostream.print_blank()

        width = 50

        cur_str = 'Time Step                : {:.3f} fs'.format(
            self.time_step)
        self.ostream.print_header(cur_str.ljust(width))
        cur_str = 'Total Duration           : {:.2f} fs'.format(
            self.duration)
        self.ostream.print_header(cur_str.ljust(width))

        self.ostream.print_blank()
        self.ostream.flush()

    def print_iteration(self, iteration=0):
        """
        Prints information about the current iteration.

        :param iteration:
            The current iteration.
        """
        width = 92
        output_header = '*** Iteration:   {:2d} '.format(iteration)
        output_header += '  Kinetic energy (H): '
        output_header += '{:7.2f}'.format(self.kinetic_energies[iteration])
        output_header += '  Potential energy (H): '
        output_header += '{:7.2f}'.format(self.potential_energies[iteration])
        output_header += '  Total energy (H): '
        output_header += '{:7.2f}'.format(self.total_energies[iteration])
        self.ostream.print_header(output_header.ljust(width))
        self.ostream.flush()
