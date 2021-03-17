import numpy as np

from .gradientdriver import GradientDriver
from .tdaorbitalresponse import TdaOrbitalResponse
from .rpaorbitalresponse import RpaOrbitalResponse
from .veloxchemlib import ElectricDipoleIntegralsDriver
from .veloxchemlib import mpi_master
from .veloxchemlib import dipole_in_debye
from .errorhandler import assert_msg_critical


class TdhfGradientDriver(GradientDriver):
    """
    Implements the analytic gradient driver for excited states at the
    Tamm-Dancoff approximation (TDA) and random phase approximation (RPA)
    level based on a Hartree-Fock ground state.
	DFT references will be implemented in the future.

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.

    Instance variables
        - scf_tensors: The results from the converged SCF calculation.
		- is_tda: Flag if Tamm-Dancoff approximation is employed.
        - n_state_deriv: The excited state of interest.
    """

    def __init__(self, comm, ostream):
        """
        Initializes gradient driver.
        """

        super().__init__(comm, ostream)
        self.rank = self.comm.Get_rank()

        self.flag = 'RPA Gradient Driver'

		# flag on whether RPA or TDA is calculated
        self.is_tda = False

        # excited state information, default to first excited state
        self.n_state_deriv = 0

    def update_settings(self, rsp_dict, method_dict=None):
        """
        Updates settings in gradient driver.

        :param rsp_dict:
            The input dictionary of response settings  group.
        :param method_dict:
            The input dicitonary of method settings group.
        """

        if method_dict is None:
            method_dict = {}

        if 'tamm_dancoff' in rsp_dict:
            if rsp_dict['tamm_dancoff'] == 'yes':
                self.is_tda = True
                self.flag = 'TDA Gradient Driver'

        # TODO: 'n_state_deriv' can be interpreted as the number of states, not
        # the index of the state. Perhaps 'state_deriv' is a better keyword.

        if 'n_state_deriv' in rsp_dict:
            # user gives '1' for first excited state, but internal index is 0
            self.n_state_deriv = int(rsp_dict['n_state_deriv']) - 1

        self.rsp_dict = rsp_dict
        self.method_dict = method_dict

    def compute(self, molecule, basis, scf_tensors, rsp_results):
        """
        Performs calculation of analytical gradient.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
		:param scf_tensors:
			The tensors from the converged SCF calculation.
		:param rsp_results:
			The results of the RPA or TDA calculation.
        """

        # sanity check for number of state
        if self.rank == mpi_master():
            assert_msg_critical(
                self.n_state_deriv < rsp_results['eigenvalues'].size,
                'TdhfGradientDriver: not enough states calculated')

        self.print_header()

        # orbital response driver
        if self.is_tda:
            orbrsp_drv = TdaOrbitalResponse(self.comm, self.ostream)
        else:
            orbrsp_drv = RpaOrbitalResponse(self.comm, self.ostream)

        orbrsp_drv.update_settings(self.rsp_dict, self.method_dict)
        orbrsp_results = orbrsp_drv.compute(molecule, basis, scf_tensors,
                                            rsp_results)

        # calculate the relaxed and unrelaxed excited-state dipole moment
        dipole_moments = self.compute_properties(molecule, basis, scf_tensors,
                                                 orbrsp_results)

        if self.rank == mpi_master():
            self.print_properties(molecule, dipole_moments)

    def compute_properties(self, molecule, basis, scf_tensors, orbrsp_results):
        """
        Calculates first-order properties of RPA or TDA excited states
        using the results of the orbital response calculation.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        :param scf_tensors:
            The tensors from the converged SCF calculation.
        :param orbrsp_results:
            The orbital response results.

        :return:
            A dictionary containing the properties.
        """

        if molecule.get_charge() != 0:
            coords = molecule.get_coordinates()
            nuclear_charges = molecule.elem_ids_to_numpy()
            origin = np.sum(coords.T * nuclear_charges,
                            axis=1) / np.sum(nuclear_charges)
        else:
            origin = np.zeros(3)

        # dipole integrals
        dipole_drv = ElectricDipoleIntegralsDriver(self.comm)
        dipole_drv.set_origin(*list(origin))
        dipole_mats = dipole_drv.compute(molecule, basis)

        if self.rank == mpi_master():
            dipole_ints = (dipole_mats.x_to_numpy(), dipole_mats.y_to_numpy(),
                           dipole_mats.z_to_numpy())

            # electronic contribution
            unrel_density = (scf_tensors['D'][0] + scf_tensors['D'][1] +
                             orbrsp_results['unrel_dm_ao'])
            unrel_electronic_dipole = -1.0 * np.array(
                [np.sum(dipole_ints[d] * unrel_density) for d in range(3)])

			# relaxed density only available if orbital response converged
            if 'rel_dm_ao' in orbrsp_results:
                rel_density = (scf_tensors['D'][0] + scf_tensors['D'][1] +
                               orbrsp_results['rel_dm_ao'])
                rel_electronic_dipole = -1.0 * np.array(
                    [np.sum(dipole_ints[d] * rel_density) for d in range(3)])

            # nuclear contribution
            coords = molecule.get_coordinates()
            nuclear_charges = molecule.elem_ids_to_numpy()
            nuclear_dipole = np.sum((coords - origin).T * nuclear_charges,
                                    axis=1)

            if 'rel_dm_ao' in orbrsp_results:
                return {
                    'unrelaxed_dipole_moment':
                        (nuclear_dipole + unrel_electronic_dipole),
                    'relaxed_dipole_moment':
                        (nuclear_dipole + rel_electronic_dipole),
                }
            else:
                return {
                    'unrelaxed_dipole_moment':
                        (nuclear_dipole + unrel_electronic_dipole),
                }
        else:
            return {}

    def print_properties(self, molecule, properties):
        """
        Prints RPA/TDA excited-state properties.

        :param molecule:
            The molecule.
        :param properties:
            The dictionary of properties.
        """

        self.ostream.print_blank()

        # prints warning if the molecule is charged
        if molecule.get_charge() != 0:
            warn_msg = '*** Warning: Molecule has non-zero charge. Dipole'
            self.ostream.print_header(warn_msg.ljust(56))
            warn_msg = '    moment will be dependent on the choice of origin.'
            self.ostream.print_header(warn_msg.ljust(56))
            warn_msg = '    Center of nuclear charge is chosen as the origin.'
            self.ostream.print_header(warn_msg.ljust(56))

        # Remove warning once DFT orbital response is implemented
        if 'xcfun' in self.method_dict:
            if self.method_dict['xcfun'] is not None:
                warn_msg = '*** Warning: Orbital response for TDDFT is not yet fully'
                self.ostream.print_header(warn_msg.ljust(56))
                warn_msg = '    implemented. Relaxed dipole moment will be wrong.'
                self.ostream.print_header(warn_msg.ljust(56))

        self.ostream.print_blank()

        prop_header = 'First-Order Properties for Excited State S'
        prop_header += str(self.n_state_deriv + 1)
        self.ostream.print_header(prop_header)
        self.ostream.print_header('-' * (len(prop_header) + 2))

        self.ostream.print_blank()

        if self.is_tda:
            title = 'TDA '
        else:
            title = 'RPA '

        title += 'Unrelaxed Dipole Moment'
        self.ostream.print_header(title)
        self.ostream.print_header('-' * (len(title) + 2))

        unrel_dip = properties['unrelaxed_dipole_moment']
        unrel_dip_au = list(unrel_dip) + [np.linalg.norm(unrel_dip)]
        unrel_dip_debye = [m * dipole_in_debye() for m in unrel_dip_au]

        for i, a in enumerate(['  X', '  Y', '  Z', 'Total']):
            valstr = '{:<5s} :'.format(a)
            valstr += '{:17.6f} a.u.'.format(unrel_dip_au[i])
            valstr += '{:17.6f} Debye   '.format(unrel_dip_debye[i])
            self.ostream.print_header(valstr)

        self.ostream.print_blank()
        self.ostream.flush()

        if 'relaxed_dipole_moment' in properties:
            if self.is_tda:
                title = 'TDA '
            else:
                title = 'RPA '

            title += 'Relaxed Dipole Moment'
            self.ostream.print_header(title)
            self.ostream.print_header('-' * (len(title) + 2))

            rel_dip = properties['relaxed_dipole_moment']
            rel_dip_au = list(rel_dip) + [np.linalg.norm(rel_dip)]
            rel_dip_debye = [m * dipole_in_debye() for m in rel_dip_au]

            for i, a in enumerate(['  X', '  Y', '  Z', 'Total']):
                valstr = '{:<5s} :'.format(a)
                valstr += '{:17.6f} a.u.'.format(rel_dip_au[i])
                valstr += '{:17.6f} Debye   '.format(rel_dip_debye[i])
                self.ostream.print_header(valstr)

            self.ostream.print_blank()
            self.ostream.flush()
        else:
            warn_msg = '*** Warning: Orbital response did not converge.'
            self.ostream.print_header(warn_msg.ljust(56))
            warn_msg = '    Relaxed dipole moment could not be calculated.'
            self.ostream.print_header(warn_msg.ljust(56))
            self.ostream.print_blank()

