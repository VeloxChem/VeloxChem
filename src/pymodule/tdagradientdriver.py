import numpy as np
import time as tm

from .molecule import Molecule
from .gradientdriver import GradientDriver
from .orbitalresponse import OrbitalResponse
from .veloxchemlib import ElectricDipoleIntegralsDriver
from .veloxchemlib import mpi_master
from .veloxchemlib import dipole_in_debye


class TdaGradientDriver(GradientDriver):
    """
    Implements the gradient driver for excited states at the
    Tamm-Dancoff approximation (TDA) level based on a Hartree-Fock
    ground state. DFT references will be implemented in the future.

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.

    Instance variables
        - scf_drv: The SCF driver.
        - tda_results: Results from the TDA driver.
        - n_state_deriv: The excited state of interest.
    """

    def __init__(self, comm, ostream, scf_drv, tda_results):
        """
        Initializes gradient driver.
        """

        super().__init__(comm, ostream)
        self.rank = comm.Get_rank()

        # excited state information, default to first excited state
        self.n_state_deriv = 0

        self.flag = 'TDA Gradient Driver'
        self.scf_drv = scf_drv
        self.tda_results = tda_results

    def update_settings(self, rsp_dict, method_dict=None):
        """
        Updates settings in gradient driver.

        :param rsp_dict:
            The input dictionary of response settings  group.
        :param method_dict:
            The input dicitonary of method settings group.
        """
        # TODO: what needs to be updated here?
        # n_state_deriv?
        if 'n_state_deriv' in rsp_dict:
            # user gives '1' for first excited state, but internal index is 0
            self.n_state_deriv = int(rsp_dict['n_state_deriv']) - 1

        return

    def compute(self, molecule, ao_basis, min_basis=None):
        """
        Performs calculation of analytical gradient.

        :param molecule:
            The molecule.
        :param ao_basis:
            The AO basis set.
        :param min_basis:
            The minimal AO basis set.
        """

        self.print_header()
        start_time = tm.time()

        scf_ostream_state = self.scf_drv.ostream.state
        self.scf_drv.ostream.state = False
        scf_tensors = self.scf_drv.scf_tensors

        # excitation energies and vectors
        exc_energies = self.tda_results["eigenvalues"]
        exc_vectors = self.tda_results["eigenvectors"]

        # orbital response driver
        orbrsp_drv = OrbitalResponse(self.comm, self.ostream)
        orbrsp_drv.update_settings(rsp_dict, method_dict) # where does he know those from?
        
        print("Calculating orbital response.")
        orbrsp_results = orbrsp_drv.compute(molecule, basis,
                                    scf_tensors, exc_vectors)
        print("Orbital response calculation done.")

        # Calculate the omega multipliers (only if nuclear gradient is requested?)
        #omega = self.compute_omega_mult(molecule, scf_tensors, exc_vectors)

        # Calculate the relaxed and unrelaxed excited-state dipole moment
        dipole_moments = self.compute_properties(molecule, scf_tensors, orbrsp_results)

        # atom labels
        #labels = molecule.get_labels()

        # atom coordinates (nx3)
        #coords = molecule.get_coordinates()

        # analytical gradient
        #self.gradient = np.zeros((molecule.number_of_atoms(), 3))

        # TODO: Contract densities and Lagrange multipliers
        #   with the integral derivatives ...

        ##self.ostream.print_blank()

        ##self.scf_drv.compute(molecule, ao_basis, min_basis)
        ##self.scf_drv.ostream.state = scf_ostream_state

        ### print gradient
        ##self.print_geometry(molecule)
        ##self.print_gradient(molecule, labels)

        ##valstr = '*** Time spent in gradient calculation: '
        ##valstr += '{:.2f} sec ***'.format(tm.time() - start_time)
        ##self.ostream.print_header(valstr)
        ##self.ostream.print_blank()
        ##self.ostream.flush()

    def compute_properties(self, molecule, basis, scf_tensors, orbrsp_results):
        """
        Calculates first-order properties of TDA excited states
        using the results of the orbital response calculation.

        :param molecule:
            The molecule.
		:param basis:
			The AO basis set.
        :param scf_tensors:
            The tensors from the converged SCF calculation.
        :param orbrsp_results:
            The results from orbital response.

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
            unrel_density = (scf_tensors['D'][0] + scf_tensors['D'][1]
                                + orbrsp_results['unrelaxed_density'])
            rel_density = (scf_tensors['D'][0] + scf_tensors['D'][1]
                                + orbrsp_results['relaxed_density'])
            unrel_electronic_dipole = -1.0 * np.array(
                [np.sum(dipole_ints[d] * unrel_density) for d in range(3)])
            rel_electronic_dipole = -1.0 * np.array(
                [np.sum(dipole_ints[d] * rel_density) for d in range(3)])

            # nuclear contribution
            coords = molecule.get_coordinates()
            nuclear_charges = molecule.elem_ids_to_numpy()
            nuclear_dipole = np.sum((coords - origin).T * nuclear_charges,
                                    axis=1)
    
            return {
                'unrelaxed_dipole_moment': (nuclear_dipole +
                                                unrel_electronic_dipole),
                'relaxed_dipole_moment': (nuclear_dipole +
                                                rel_electronic_dipole),
            }

    def compute_omega_mult(self, molecule, scf_tensors, exc_vectors):
        """
        Calculates the Lagrange multipliers for the overlap matrix.

        :param molecule:
            The molecule.
        :param scf_tensors:
            The tensors from the converged SCF calculation.
        :param exc_vectors:
            The excitation vectors from the TDA calculation.

        :return:
            A dictionary containing the Lagrange multipliers.
        """
        # Orbital information 
        nocc = molecule.number_of_alpha_electrons() 
        mo_occ = scf_tensors['C'][:, :nocc] # occupied MO coefficients
        mo_vir = scf_tensors['C'][:, nocc:] # virtual MO coefficients
        nvir = mo_vir.shape[1]
        ovlp = scf_tensors['S'] # overlap matrix
        eocc = scf_tensors['E'][:nocc]
        evir = scf_tensors['E'][nocc:]

        # Take vector of interest and convert to matrix form
        exc_vec = exc_vectors[:, self.n_state_deriv].copy().reshape(nocc, nvir)

        #TODO: to be continued...



    def print_properties(self, molecule, properties):
        """
        Prints TDA excited-state properties.

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

        self.ostream.print_blank()

        title = 'Unrelaxed Dipole Moment'
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

        title = 'Relaxed Dipole Moment'
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
