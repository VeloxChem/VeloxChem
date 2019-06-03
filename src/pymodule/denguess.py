from os.path import isfile
import time as tm

from .veloxchemlib import OverlapIntegralsDriver
from .veloxchemlib import SADGuessDriver
from .veloxchemlib import mpi_master
from .molecularorbitals import MolecularOrbitals
from .aodensitymatrix import AODensityMatrix


class DensityGuess:
    """Implements initial density guess generator.

    Implements initial density guess generator with set of different methods
    for generation of initial density.

    Attributes
    ----------
    _guess_type
        The type of initial guess.
    """

    def __init__(self, guess_type='SAD', checkpoint_file=None):
        """Initializes initial guess generator.

        Initializes initial guess generator by setting initial guess type.
        Default for initial guess type is a superposition of atomic densities.
        """

        self._guess_type = guess_type
        self._checkpoint_file = checkpoint_file

    def __str__(self):
        """Converts object to it's informal string"""
        return self._guess_type

    def __repr__(self):
        """Converts object to it's formal string"""
        return self._guess_type

    @property
    def guess_type(self):
        """Getter function for protected guess_type attribute."""
        return self._guess_type

    @guess_type.setter
    def guess_type(self, value):
        """Setter function for protected guess_type attribute."""
        self._guess_type = value

    def validate_checkpoint(self, rank, comm, nuclear_charges, basis_set):
        """Validates the checkpoint file.

        Validates the checkpoint file by checking nuclear charges and basis set.

        Parameters
        ----------
        rank
            The rank of the MPI process.
        comm
            The MPI communicator.
        nuclear_charges
            Numpy array of the nuclear charges.
        basis_set
            Name of the AO basis.
        """

        if (self._checkpoint_file and isinstance(self._checkpoint_file, str) and
                rank == mpi_master() and isfile(self._checkpoint_file)):
            valid = MolecularOrbitals.match_hdf5(self._checkpoint_file,
                                                 nuclear_charges, basis_set)
        else:
            valid = False

        valid = comm.bcast(valid, root=mpi_master())

        return valid

    def restart_density(self, molecule, rank, ostream):
        """Reads initial AO density from checkpoint file.

        Reads initial molecular orbitals and AO density from checkpoint file.

        Parameters
        ----------
        molecule
            The molecule.
        rank
            The rank of the MPI process.
        ostream
            The output stream.
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

    def sad_density(self, molecule, ao_basis, min_basis, overlap_matrix, comm,
                    ostream):
        """Generates initial AO density using SAD scheme.

        Computes initial AO density using superposition of atomic densities
        scheme.

        Parameters
        ----------
        molecule
            The molecule.
        ao_basis
            The AO basis.
        min_basis
            The minimal AO basis for generation of atomic densities.
        overlap_matrix
            The AO overlap matrix.
        comm
            The local MPI communicator.
        ostream
            The output stream.
        """

        if self._guess_type == 'SAD':

            ovl_drv = OverlapIntegralsDriver(comm)

            ovl_mat_sb = ovl_drv.compute(molecule, min_basis, ao_basis)

            t0 = tm.time()

            sad_drv = SADGuessDriver(comm)

            den_mat = sad_drv.compute(molecule, min_basis, ao_basis, ovl_mat_sb,
                                      overlap_matrix)

            if comm.Get_rank() == mpi_master():

                ostream.print_info('SAD initial guess computed in %.2f sec.' %
                                   (tm.time() - t0))

                ostream.print_blank()

                ostream.flush()

            return den_mat

        return AODensityMatrix()

    def prcmo_density(self, molecule, ao_basis, red_basis, red_orbs):
        """Generates initial AO density using PRCMO scheme.

        Computes initial AO density from molecular orbitals obtained by
        inserting molecular orbitals from reduced basis into molecular
        orbitals in full AO basis.

        Parameters
        ----------
        molecule
            The molecule.
        ao_basis
            The AO basis.
        red_basis
            The reduced AO basis for generation of molecular orbitals.
        red_orbs
            The molecular orbitals in reduced AO basis.
        """

        if self._guess_type == 'PRCMO':

            proj_orbs = red_orbs.insert(molecule, ao_basis, red_basis)

            return proj_orbs.get_density(molecule)

        return AODensityMatrix()
