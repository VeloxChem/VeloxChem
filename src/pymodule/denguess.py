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

    def validate_checkpoint(self, molecule, ao_basis, comm, ovl_thresh):
        """Validates the checkpoint file.

        Validates the checkpoint file by checking number of AOs and MOs.

        Parameters
        ----------
        molecule
            The molecule.
        ao_basis
            The AO basis.
        comm
            The MPI communicator.
        ovl_thresh
            The atomic orbitals linear dependency threshold.
        """

        rank = comm.Get_rank()
        size = comm.Get_size()
        ovl_drv = OverlapIntegralsDriver(rank, size, comm)
        ovl_mat = ovl_drv.compute(molecule, ao_basis, comm)

        if rank == mpi_master():
            oao_mat = ovl_mat.get_ortho_matrix(ovl_thresh)
            nao = oao_mat.number_of_rows()
            nmo = oao_mat.number_of_columns()

        valid = False
        if self._checkpoint_file and isinstance(self._checkpoint_file, str):
            if rank == mpi_master() and isfile(self._checkpoint_file):
                mol_orbs = MolecularOrbitals.read_hdf5(self._checkpoint_file)
                if (mol_orbs.number_aos() == nao and
                        mol_orbs.number_mos() == nmo):
                    valid = True

        valid = comm.bcast(valid, root=mpi_master())
        return valid

    def restart_density(self, molecule, comm, ostream):
        """Reads initial AO density from checkpoint file.

        Reads initial molecular orbitals and AO density from checkpoint file.

        Parameters
        ----------
        molecule
            The molecule.
        comm
            The MPI communicator.
        ostream
            The output stream.
        """

        if comm.Get_rank() == mpi_master():
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
                    loc_rank, loc_nodes, comm, ostream):
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
        loc_rank
            The local MPI rank.
        loc_nodes
            The local number of MPI processes.
        comm
            The local MPI communicator.
        ostream
            The output stream.
        """

        if self._guess_type == 'SAD':

            ovl_drv = OverlapIntegralsDriver(loc_rank, loc_nodes, comm)

            ovl_mat_sb = ovl_drv.compute(molecule, min_basis, ao_basis, comm)

            t0 = tm.time()

            sad_drv = SADGuessDriver(comm)

            den_mat = sad_drv.compute(molecule, min_basis, ao_basis, ovl_mat_sb,
                                      overlap_matrix)

            if loc_rank == mpi_master():

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
