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
        - dipole_moment: The electric dipole moment.
        - TODO: what else?
    """

    def __init__(self, comm, ostream, scf_tensors):
        """
        Initializes the one-electron properties driver.
        """

        self.comm = comm
        self.ostream = ostream

        #self.scf_drv = ScfDriver(self.comm, OutputStream()) #Restricted
        self.scf_tensors = scf_tensors

        self.dipole_moment = None
        #TODO more things to add

    def update_settings(self, scf_dict, method_dict=None):
        """
        Updates settings in gradient driver.

        :param scf_dict:
            The input dictionary of scf group.
        :param method_dict:
            The input dicitonary of method settings group.
        """

        self.scf_drv.update_settings(scf_dict, method_dict)

    def comp_nuclear_dipole_moment(self, molecule):
        """
        Computes the nuclear part of the electric dipole moment.

        :param molecule:
            The molecule.

        :return:
            The Cartesian components of the nuclear contribution
            to the electric dipole moment.
        """
        #nuc_dipmom = []
        
        # Get the Cartesian coordinates and nuclear charges
        #x_coords = molecule.x_to_numpy()
        #y_coords = molecule.y_to_numpy()
        #z_coords = molecule.z_to_numpy()
        #coords = np.vstack(
        #            (molecule.x_to_numpy(),
        #            molecule.y_to_numpy(),
        #            molecule.z_to_numpy())).T
        coords = molecule.get_coordinates()
        nuc_charges = molecule.elem_ids_to_numpy()

        # Calculate nuclear dipole moment components
        nuc_dipmom = np.einsum('i,ix->x', nuc_charges, coords)
        #for coords in x_coords, y_coords, z_coords:
        #    nuc_dipmoms = np.dot(nuc_charges, coords)
        #    nuc_dipmom.append(nuc_dipmoms)

        print("Nuclear dipole moment contribution: ", nuc_dipmom)
        return nuc_dipmom

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

    #MH Compute the ground-state dipole moment
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
        #print("dipole_ints = ", dipole_ints)

        # Get the AO density matrix
        # TODO: perhaps there is a factor of 2 missing, since it only takes
        # the alpha part of the density!?
        D = self.scf_tensors['D'][0]

        # Calculate the electronic contribution
        electronic_dipole = []
        
        for d in range(3):
            elec_dipole = np.sum(dipole_ints[d] * D)
            #elec_dipole = np.trace(D.dot(dipole_ints[d]))
            electronic_dipole.append(elec_dipole)

        print("Electronic contribution to dipole moment: ", electronic_dipole)

        nuclear_dipole = self.comp_nuclear_dipole_moment(molecule)

        # Element-wise add nuclear and electronic contributions
        #dipole_moment = [a + b for a, b in zip(nuclear_dipole, elec_dipole)]
        dipole_moment = [electronic_dipole[i] + nuclear_dipole[i] for i in range(3)]

        au2debye = 2.541746229
        total_dipole = np.linalg.norm(dipole_moment)
        print("Total dipole =", total_dipole, "a.u. =", total_dipole * au2debye, "D")

        return dipole_moment

