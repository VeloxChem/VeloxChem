import numpy as np
import time as tm

from .molecule import Molecule
from .outputstream import OutputStream
from .scfrestdriver import ScfDriver    #Restricted ?
from .veloxchemlib import ElectricDipoleIntegralsDriver
from .veloxchemlib import mpi_master

class ScfProperties:
    """
    Implements the calculation of first-order properties
    for the SCF ground state, such as the dipole moment.

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.

    Instance variables:
        - density: The AO density matrix.
        - dipole_moment: The electric dipole moment.
        - au2debye: Conversion factor from atomic units to Debye.
    """

    def __init__(self, comm, ostream):
        """
        Initializes the one-electron properties driver.
        """

        self.comm = comm
        self.ostream = ostream

        # Density matrix in the AO basis
        self.density = None

        self.dipole_moment = None
        self.au2debye = 2.541746229 # keep it for the moment

    def comp_nuclear_dipole_moment(self, molecule):
        """
        Computes the nuclear part of the electric dipole moment.

        :param molecule:
            The molecule.

        :return:
            The Cartesian components of the nuclear contribution
            to the electric dipole moment.
        """
        # Atomic coordinates (nx3)
        coords = molecule.get_coordinates()

        # Nuclear charges
        nuc_charges = molecule.elem_ids_to_numpy()

        # Calculate nuclear dipole moment components
        nuc_dipmom = np.einsum('i,ix->x', nuc_charges, coords)
        
        if self.comm.Get_rank() == mpi_master():
            return nuc_dipmom
        else:
            return None

    def comp_dipole_ints(self, molecule, ao_basis):
        """
        Computes one-electron dipole integrals.

        :param molecule:
            The molecule.
        :param ao_basis:
            The AO basis set.

        :return:
            The Cartesian components of one-electron dipole integrals.
        """

        dipole_drv = ElectricDipoleIntegralsDriver(self.comm)
        dipole_mats = dipole_drv.compute(molecule, ao_basis)

        if self.comm.Get_rank() == mpi_master():
            return (dipole_mats.x_to_numpy(), dipole_mats.y_to_numpy(),
                    dipole_mats.z_to_numpy())
        else:
            return None

    def compute(self, molecule, basis, scf_tensors):
        """
        Computes the SCF ground-state dipole moment

        :param molecule:
            The molecule
        :param basis:
            The AO basis set.
        :param scf_tensors:
            The tensors from the converged SCF calculation.

        :return:
            Cartesian components of the ground-state dipole moment
        """
        # Calculate the dipole integrals in the AO basis on the master node
        if self.comm.Get_rank() == mpi_master:
            dipole_ints = self.comp_dipole_ints(molecule, basis)

        # The total electron density is the alpha plus the beta part
        self.density = scf_tensors['D'][0] + scf_tensors['D'][1]

        # Calculate the electronic contribution
        electronic_dipole = []
        
        for d in range(3):
            elec_dipole = np.sum(dipole_ints[d] * self.density)
            electronic_dipole.append(elec_dipole)

        # Calculate the nuclear contribution
        nuclear_dipole = self.comp_nuclear_dipole_moment(molecule)

        # Element-wise add (or subtract because of the electrons' negative charge)
        # the nuclear and electronic contributions
        self.dipole_moment = [a - b for a, b in zip(nuclear_dipole, electronic_dipole)]
        #dipole_moment = [nuclear_dipole[i] - electronic_dipole[i] for i in range(3)]

        if self.comm.Get_rank() == mpi_master():
            return np.array(self.dipole_moment)
        else:
            return None

    def get_dipole_moment(self):
        """
        Gets the electric dipole moment (in atomic units).

        :return:
            The electric dipole moment.
        """

        return self.dipole_moment

    def print_scf_properties(self, molecule): #, basis):
        """
        Print first-order SCF ground-state properties.
        So far this includes only the electric dipole moment.

        :param molecule:
            The molecule
        """
        
        self.ostream.print_blank()

        self.ostream.print_header("Ground-State Properties".ljust(92))
        self.ostream.print_header("------------------------".ljust(92))
       
        # Check if the system is charged. If so, print a warning
        if molecule.get_charge() != 0:
            self.ostream.print_header("""
                *** Warning: The dipole moment for charged molecules is not well defined.""".ljust(92))
            self.ostream.print_blank()
            self.ostream.flush()

        # Format the electric dipole moment to have only four decimal digits
        fmtd_dipmom = [float('{:.4f}'.format(d)) for d in self.dipole_moment]

        # Print the results
        valstr = "Dipole Moment [a.u.]      : {}".format(fmtd_dipmom)
        self.ostream.print_header(valstr.ljust(92))
        total_dipole = np.linalg.norm(self.dipole_moment)
        valstr = "Total Dipole [Debye]      :{:7.4f}".format(total_dipole * self.au2debye)
        self.ostream.print_header(valstr.ljust(92))
        
        self.ostream.print_blank()
        self.ostream.flush()

