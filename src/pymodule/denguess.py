from .veloxchemlib import OverlapMatrix
from .veloxchemlib import OverlapIntegralsDriver
from .veloxchemlib import SADGuessDriver

from .aodensitymatrix import AODensityMatrix

import numpy as np
import time as tm

class DensityGuess:
    """Implements initial density guess generator.
        
    Implements initial density guess generator with set of different methods
    for generation of initial density.
        
    Attributes
    ----------
    _guess_type
        The type of initial guess.
    """

    def __init__(self, guess_type = "SAD"):
        """Initializes initial guess generator.
            
        Initializes initial guess generator by setting initial guess type.
        Default for initial guess type is a superposition of atomic densities.
        """
    
        self._guess_type = guess_type
    
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
        
        if self._guess_type == "SAD":
            
            ovl_drv = OverlapIntegralsDriver.create(loc_rank, loc_nodes, comm)
            
            ovl_mat_sb = ovl_drv.compute(molecule, min_basis, ao_basis, comm)

            t0 = tm.time()
                                            
            sad_drv = SADGuessDriver.create(loc_rank, loc_nodes, comm)
            
            den_mat = sad_drv.compute(molecule, min_basis, ao_basis, ovl_mat_sb,
                                      overlap_matrix, comm)

            ostream.print_info("SAD initial guess computed in %.2f sec." %
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
        
        if self._guess_type == "PRCMO":

            proj_orbs = red_orbs.insert(molecule, ao_basis, red_basis)
        
            # add handling of unrestricted case
            return proj_orbs.get_rest_density(molecule)

        return AODensityMatrix()

    
