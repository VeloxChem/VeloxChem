import numpy as np
import time as tm

from .molecule import Molecule
from .outputstream import OutputStream
from .scfrestdriver import ScfRestrictedDriver
from .veloxchemlib import ElectricDipoleIntegralsDriver #MH

class OneElProperties:
    """
    Implements the calculation of one-electron properties
    for the SCF ground state, such as the dipole moment.

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.

    Instance variables:
        - scf_drv: The SCF driver.
        - dipole_moment: The electric dipole moment.
        - TODO: what else?
    """

    def __init__(self, comm, ostream):
        """
        Initializes the one-electron properties driver.
        """

        self.comm = comm
        self.ostream = ostream

        self.scf_drv = ScfRestrictedDriver(self.comm, OutputStream())

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

        print("scfdriver.py: comp_dipole_ints")
        dipole_drv = ElectricDipoleIntegralsDriver(self.comm)
        dipole_mats = dipole_drv.compute(molecule, ao_basis)

        if self.rank == mpi_master():
            return (dipole_mats.x_to_numpy(), dipole_mats.y_to_numpy(),
                    dipole_mats.z_to_numpy())
        else:
            return ()

