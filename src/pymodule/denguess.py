#
#                           VELOXCHEM 1.0-RC
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

from os.path import isfile
import time as tm

from .veloxchemlib import OverlapIntegralsDriver
from .veloxchemlib import SADGuessDriver
from .veloxchemlib import mpi_master
from .molecularorbitals import MolecularOrbitals
from .aodensitymatrix import AODensityMatrix


class DensityGuess:
    """
    Implements initial density guess generator with set of different methods
    for generation of initial density.

    :param guess_type:
        The type of initial guess.
    :param checkpoint_file:
        The name of checkpoint file.
    """

    def __init__(self, guess_type='SAD', checkpoint_file=None):
        """
        Initializes initial guess generator by setting initial guess type.
        Default for initial guess type is a superposition of atomic densities.
        """

        self._guess_type = guess_type
        self._checkpoint_file = checkpoint_file

    def __str__(self):
        """
        Converts object to it's informal string.
        """

        return self._guess_type

    def __repr__(self):
        """
        Converts object to it's formal string.
        """

        return self._guess_type

    @property
    def guess_type(self):
        """
        Getter function for protected guess_type attribute.
        """

        return self._guess_type

    @guess_type.setter
    def guess_type(self, value):
        """
        Setter function for protected guess_type attribute.
        """

        self._guess_type = value

    def validate_checkpoint(self, rank, comm, nuclear_charges, basis_set,
                            restricted):
        """
        Validates the checkpoint file by checking nuclear charges and basis set.

        :param rank:
            The rank of the MPI process.
        :param comm:
            The MPI communicator.
        :param nuclear_charges:
            Numpy array of the nuclear charges.
        :param basis_set:
            Name of the AO basis.
        :param restricted:
            The flag for restricted molecular orbitals.

        :return:
            Validity of the checkpoint file.
        """

        if (self._checkpoint_file and isinstance(self._checkpoint_file, str) and
                rank == mpi_master() and isfile(self._checkpoint_file)):
            valid = MolecularOrbitals.match_hdf5(self._checkpoint_file,
                                                 nuclear_charges, basis_set,
                                                 restricted)
        else:
            valid = False

        valid = comm.bcast(valid, root=mpi_master())

        return valid

    def restart_density(self, molecule, rank, ostream):
        """
        Reads initial molecular orbitals and AO density from checkpoint file.

        :param molecule:
            The molecule.
        :param rank:
            The rank of the MPI process.
        :param ostream:
            The output stream.

        :return:
            The AO density matrix to restart from.
        """

        if rank == mpi_master():
            mol_orbs = MolecularOrbitals.read_hdf5(self._checkpoint_file)
            den_mat = mol_orbs.get_density(molecule)

            restart_text = 'Restarting from checkpoint file: '
            restart_text += self._checkpoint_file
            ostream.print_info(restart_text)
            ostream.print_blank()

        else:
            den_mat = AODensityMatrix()

        return den_mat

    def sad_density(self, molecule, ao_basis, min_basis, overlap_matrix,
                    restricted, comm, ostream):
        """
        Computes initial AO density using superposition of atomic densities
        scheme.

        :param molecule:
            The molecule.
        :param ao_basis:
            The AO basis.
        :param min_basis:
            The minimal AO basis for generation of atomic densities.
        :param overlap_matrix:
            The AO overlap matrix.
        :param restricted:
            The flag for generating restricted initial guess.
        :param comm:
            The local MPI communicator.
        :param ostream:
            The output stream.

        :return:
            The AO density matrix from SAD initial guess.
        """

        if self._guess_type == 'SAD':

            ovl_drv = OverlapIntegralsDriver(comm)

            ovl_mat_sb = ovl_drv.compute(molecule, min_basis, ao_basis)

            t0 = tm.time()

            sad_drv = SADGuessDriver(comm)

            den_mat = sad_drv.compute(molecule, min_basis, ao_basis, ovl_mat_sb,
                                      overlap_matrix, restricted)

            if comm.Get_rank() == mpi_master():

                ostream.print_info('SAD initial guess computed in %.2f sec.' %
                                   (tm.time() - t0))

                ostream.print_blank()

                ostream.flush()

            return den_mat

        return AODensityMatrix()

    def prcmo_density(self, molecule, ao_basis, red_basis, red_orbs):
        """
        Computes initial AO density from molecular orbitals obtained by
        inserting molecular orbitals from reduced basis into molecular
        orbitals in full AO basis.

        :param molecule:
            The molecule.
        :param ao_basis:
            The AO basis.
        :param red_basis:
            The reduced AO basis for generation of molecular orbitals.
        :param red_orbs:
            The molecular orbitals in reduced AO basis.

        :return:
            The AO density matrix from PRCMO.
        """

        if self._guess_type == 'PRCMO':

            proj_orbs = red_orbs.insert(molecule, ao_basis, red_basis)

            return proj_orbs.get_density(molecule)

        return AODensityMatrix()
