import numpy as np
import time as tm

from .molecule import Molecule
from .outputstream import OutputStream
from .scfrestdriver import ScfDriver    #Restricted
from .veloxchemlib import ElectricDipoleIntegralsDriver #MH
from .veloxchemlib import mpi_master

class OneElectronProperties:
    """
    Implements the calculation of one-electron properties
    for the SCF ground state, such as the dipole moment.

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.
    :param scf_tensors:
        The tensors from the converged SCF.

    Instance variables:
        - scf_tensors: The tensors from the converged SCF result.
        - density: The AO density matrix.
        - dipole_moment: The electric dipole moment.
        - total_dipole: The total electric dipole moment.
        - au2debye: Conversion factor from atomic units to Debye.

        - TODO: what else?
    """

    def __init__(self, comm, ostream, scf_tensors):
        """
        Initializes the one-electron properties driver.
        """

        self.comm = comm
        self.ostream = ostream

        self.scf_tensors = scf_tensors

        # Get the SCF density matrix in the AO basis
        # Total density = alpha + beta part
        self.density = self.scf_tensors['D'][0] + self.scf_tensors['D'][1]

        self.dipole_moment = None
        self.total_dipole = 0.0
        self.au2debye = 2.541746229
        #TODO more things to add

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
            return np.array([])

    #MH (from tdaexcidriver.py)
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
            return ()

    def compute_dipole_moment(self, molecule, basis):
        """
        Computes the SCF ground-state dipole moment

        :param molecule:
            The molecule
        :param basis:
            The AO basis set.
        #TODO: tbc

        :return:
            Cartesian components of the ground-state dipole moment
        """
        # Calculate the dipole integrals in the AO basis
        dipole_ints = self.comp_dipole_ints(molecule, basis)

        # Calculate the electronic contribution
        electronic_dipole = []
        
        for d in range(3):
            elec_dipole = np.sum(dipole_ints[d] * self.density)
            electronic_dipole.append(elec_dipole)

        nuclear_dipole = self.comp_nuclear_dipole_moment(molecule)

        # Element-wise add (or subtract because of the electrons' negative charge)
        # the nuclear and electronic contributions
        self.dipole_moment = [a - b for a, b in zip(nuclear_dipole, electronic_dipole)]
        #dipole_moment = [nuclear_dipole[i] - electronic_dipole[i] for i in range(3)]

        self.total_dipole = np.linalg.norm(self.dipole_moment)

        if self.comm.Get_rank() == mpi_master():
            return np.array(self.dipole_moment)
        else:
            return np.array([])

    def get_dipole_moment(self):
        """
        Gets the electric dipole moment (in atomic units).

        :return:
            The electric dipole moment.
        """

        return self.dipole_moment

    def get_total_dipole(self):
        """
        Gets the total electric dipole moment (in atomic units).

        :return:
            The total electric dipole moment.
        """

        return self.total_dipole

    def compute(self, molecule, basis):
        """
        Computes SCF ground-state properties.
        So far this includes only the electric dipole moment.

        :param molecule:
            The molecule
        :param basis:
            The AO basis set.
        """
        
        self.ostream.print_blank()

        self.ostream.print_header("Ground-State Properties".ljust(92))
        self.ostream.print_header("------------------------".ljust(92))
       
        # Check if the system is charged. If so, print a warning
        chg = molecule.get_charge()
        if chg != 0:
            self.ostream.print_header("""
                *** Warning: The dipole moment for charged molecules is ill-defined.""".ljust(92))
            self.ostream.flush()

        # Compute the electric dipole moment
        dipmom = self.compute_dipole_moment(molecule, basis)
        fmtd_dipmom = [float('{:.4f}'.format(d)) for d in dipmom]

        # Print the results
        valstr = "Dipole Moment [a.u.]      : {}".format(fmtd_dipmom)
        self.ostream.print_header(valstr.ljust(92))
        valstr = "Total Dipole [Debye]      :{:7.4f}".format(self.total_dipole * self.au2debye)
        self.ostream.print_header(valstr.ljust(92))
        
        self.ostream.print_blank()
        self.ostream.flush()

