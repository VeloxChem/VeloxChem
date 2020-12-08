import numpy as np

from .veloxchemlib import ElectricDipoleIntegralsDriver
from .veloxchemlib import mpi_master


class ScfProperties:
    """
    Implements the calculation of first-order properties for the SCF ground
    state, such as the dipole moment.

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.

    Instance variables:
        - dipole_moment: The electric dipole moment.
        - au2debye: Conversion factor from atomic units to Debye.
    """

    def __init__(self, comm, ostream):
        """
        Initializes the SCF properties.
        """

        self.comm = comm
        self.rank = comm.Get_rank()

        self.ostream = ostream

        self.dipole_moment = None

        self.au2debye = 2.541746229  # keep it for the moment

    def compute(self, molecule, basis, scf_tensors):
        """
        Computes the SCF ground-state dipole moment.

        :param molecule:
            The molecule
        :param basis:
            The AO basis set.
        :param scf_tensors:
            The tensors from the converged SCF calculation.
        """

        # dipole integrals
        dipole_drv = ElectricDipoleIntegralsDriver(self.comm)
        dipole_mats = dipole_drv.compute(molecule, basis)

        if self.rank == mpi_master():
            dipole_ints = (dipole_mats.x_to_numpy(), dipole_mats.y_to_numpy(),
                           dipole_mats.z_to_numpy())

            # electronic contribution
            total_density = scf_tensors['D'][0] + scf_tensors['D'][1]
            electronic_dipole = np.array(
                [np.sum(dipole_ints[d] * total_density) for d in range(3)])

            # nuclear contribution
            coords = molecule.get_coordinates()
            nuclear_charges = molecule.elem_ids_to_numpy()
            nuclear_dipole = np.sum(coords.T * nuclear_charges, axis=1)

            self.dipole_moment = nuclear_dipole - electronic_dipole

            return {'dipole moment': self.dipole_moment}

        return {}

    def print_scf_properties(self, molecule):
        """
        Prints SCF ground-state properties.

        :param molecule:
            The molecule.
        """

        self.ostream.print_blank()

        title = 'Ground-State Dipole Moment'
        self.ostream.print_header(title)
        self.ostream.print_header('-' * (len(title) + 2))

        # prints warning if the molecule is charged
        if molecule.get_charge() != 0:
            warn_msg = '*** Warning: Molecule has non-zero charge.'
            self.ostream.print_header(warn_msg)
            warn_msg = '*** Dipole moment will be gauge-dependent.'
            self.ostream.print_header(warn_msg)

        self.ostream.print_blank()

        dip_au = list(self.dipole_moment) + [np.linalg.norm(self.dipole_moment)]
        dip_debye = [m * self.au2debye for m in dip_au]

        for i, a in enumerate(['  X', '  Y', '  Z', 'Total']):
            valstr = '{:<5s} :'.format(a)
            valstr += '{:17.6f} a.u.'.format(dip_au[i])
            valstr += '{:17.6f} Debye   '.format(dip_debye[i])
            self.ostream.print_header(valstr)

        self.ostream.print_blank()
        self.ostream.flush()
