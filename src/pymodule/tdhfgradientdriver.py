import numpy as np
import time as tm

from .molecule import Molecule
from .gradientdriver import GradientDriver
from .tdaorbitalresponse import TdaOrbitalResponse
from .rpaorbitalresponse import RpaOrbitalResponse
from .veloxchemlib import ElectricDipoleIntegralsDriver
from .veloxchemlib import mpi_master
from .veloxchemlib import dipole_in_debye
from .errorhandler import assert_msg_critical
from .firstorderprop import FirstOrderProperties


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
        - gradient: The gradient.
        - scf_tensors: The results from the converged SCF calculation.
        - is_tda: Flag if Tamm-Dancoff approximation is employed.
        - n_state_deriv: The excited state of interest.
        - do_first_order_prop: Controls the calculation of first-order properties.
        - delta_h: The displacement for finite difference.
        - do_four_point: Flag for four-point finite difference.
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

        # flag on whether to calculate excited-state properties
        self.do_first_order_prop = False

        # for numerical gradient
        self.gradient = None
        self.delta_h = 0.001
        #self.scf_drv = None #scf_drv
        #self.rsp_drv = None
        # flag for two-point or four-point approximation
        self.do_four_point = False


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
            key = rsp_dict['tamm_dancoff'].lower()
            self.is_tda = True if key in ['yes', 'y'] else False

        if self.is_tda:
            self.flag = 'TDA Gradient Driver'

        if 'do_first_order_prop' in rsp_dict:
            key = rsp_dict['do_first_order_prop'].lower()
            self.do_first_order_prop = True if key in ['yes', 'y'] else False

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
            method = 'TDA'
            orbrsp_drv = TdaOrbitalResponse(self.comm, self.ostream)
        else:
            orbrsp_drv = RpaOrbitalResponse(self.comm, self.ostream)
            method = 'RPA'

        orbrsp_drv.update_settings(self.rsp_dict, self.method_dict)
        orbrsp_results = orbrsp_drv.compute(molecule, basis, scf_tensors,
                                            rsp_results)

        # If desired, calculate the relaxed and unrelaxed excited-state dipole moment
        if self.do_first_order_prop:
            firstorderprop = FirstOrderProperties(self.comm, self.ostream)

            unrel_density = (scf_tensors['D'][0] + scf_tensors['D'][1] +
                             orbrsp_results['unrel_dm_ao'])
            firstorderprop.compute(molecule, basis, unrel_density)
            if self.rank == mpi_master():
                title = method + ' Unrelaxed Dipole Moment for Excited State ' + str(self.n_state_deriv + 1)
                firstorderprop.print_properties(molecule, title)

            if 'rel_dm_ao' in orbrsp_results:
                rel_density = (scf_tensors['D'][0] + scf_tensors['D'][1] +
                               orbrsp_results['rel_dm_ao'])
                firstorderprop.compute(molecule, basis, rel_density)
                if self.rank == mpi_master():
                    # TODO: Remove warning once TDDFT orbital response is fully implemented
                    if 'xcfun' in self.method_dict and self.method_dict['xcfun'] is not None:
                        warn_msg = '*** Warning: Orbital response for TDDFT is not yet fully'
                        self.ostream.print_header(warn_msg.ljust(56))
                        warn_msg = '    implemented. Relaxed dipole moment will be wrong.'
                        self.ostream.print_header(warn_msg.ljust(56))

                    title = method + ' Relaxed Dipole Moment for Excited State ' + str(self.n_state_deriv + 1)
                    firstorderprop.print_properties(molecule, title)


    def compute_numerical_gradient(self, molecule, ao_basis, scf_drv,
                                   rsp_drv, min_basis=None):
        """
        Performs calculation of numerical gradient at RPA or TDA level.

        :param molecule:
            The molecule.
        :param ao_basis:
            The AO basis set.
        :param scf_drv:
            The SCF driver.
        :param rsp_drv:
            The response (RPA or TDA) driver.
        :param min_basis:
            The minimal AO basis set.
        """

        # self.print_header()
        start_time = tm.time()

        self.scf_drv = scf_drv
        scf_ostream_state = self.scf_drv.ostream.state
        self.scf_drv.ostream.state = False

        # This does not have any influence, but depends
        # on the ostream state of the scf driver
        self.rsp_drv = rsp_drv
        # rsp_ostream_state = self.rsp_drv.ostream.state
        # self.rsp_drv.ostream.state = False

        # atom labels
        labels = molecule.get_labels()

        # atom coordinates (nx3)
        coords = molecule.get_coordinates()

        # numerical gradient
        self.gradient = np.zeros((molecule.number_of_atoms(), 3))

        if not self.do_four_point:
            for i in range(molecule.number_of_atoms()):
                for d in range(3):
                    coords[i, d] += self.delta_h
                    new_mol = Molecule(labels, coords, units='au')
                    self.scf_drv.compute(new_mol, ao_basis, min_basis)
                    scf_tensors = self.scf_drv.scf_tensors
                    self.rsp_drv.is_converged = False  # only needed for RPA
                    rsp_results = self.rsp_drv.compute(new_mol, ao_basis,
                                                       scf_tensors)
                    exc_en_plus = rsp_results['eigenvalues'][self.n_state_deriv]
                    e_plus = self.scf_drv.get_scf_energy() + exc_en_plus

                    coords[i, d] -= 2.0 * self.delta_h
                    new_mol = Molecule(labels, coords, units='au')
                    self.scf_drv.compute(new_mol, ao_basis, min_basis)
                    self.rsp_drv.is_converged = False
                    rsp_results = self.rsp_drv.compute(new_mol, ao_basis,
                                                       self.scf_drv.scf_tensors)
                    exc_en_minus = rsp_results['eigenvalues'][self.n_state_deriv]
                    e_minus = self.scf_drv.get_scf_energy() + exc_en_minus

                    coords[i, d] += self.delta_h
                    self.gradient[i, d] = (e_plus - e_minus) / (2.0 * self.delta_h)

        else:
            # Four-point numerical derivative approximation
            # for debugging of analytical gradient:
            # [ f(x - 2h) - 8 f(x - h) + 8 f(x + h) - f(x + 2h) ] / ( 12h )
            for i in range(molecule.number_of_atoms()):
                for d in range(3):
                    coords[i, d] += self.delta_h
                    new_mol = Molecule(labels, coords, units='au')
                    self.scf_drv.compute(new_mol, ao_basis, min_basis)
                    scf_tensors = self.scf_drv.scf_tensors
                    self.rsp_drv.is_converged = False  # only needed for RPA
                    rsp_results = self.rsp_drv.compute(new_mol, ao_basis,
                                                       scf_tensors)
                    exc_en = rsp_results['eigenvalues'][self.n_state_deriv]
                    e_plus1 = self.scf_drv.get_scf_energy() + exc_en

                    coords[i, d] += self.delta_h
                    new_mol = Molecule(labels, coords, units='au')
                    self.scf_drv.compute(new_mol, ao_basis, min_basis)
                    scf_tensors = self.scf_drv.scf_tensors
                    self.rsp_drv.is_converged = False  # only needed for RPA
                    rsp_results = self.rsp_drv.compute(new_mol, ao_basis,
                                                       scf_tensors)
                    exc_en = rsp_results['eigenvalues'][self.n_state_deriv]
                    e_plus2 = self.scf_drv.get_scf_energy() + exc_en

                    coords[i, d] -= 3.0 * self.delta_h
                    new_mol = Molecule(labels, coords, units='au')
                    self.scf_drv.compute(new_mol, ao_basis, min_basis)
                    self.rsp_drv.is_converged = False
                    rsp_results = self.rsp_drv.compute(new_mol, ao_basis,
                                                       self.scf_drv.scf_tensors)
                    exc_en = rsp_results['eigenvalues'][self.n_state_deriv]
                    e_minus1 = self.scf_drv.get_scf_energy() + exc_en

                    coords[i, d] -= self.delta_h
                    new_mol = Molecule(labels, coords, units='au')
                    self.scf_drv.compute(new_mol, ao_basis, min_basis)
                    self.rsp_drv.is_converged = False
                    rsp_results = self.rsp_drv.compute(new_mol, ao_basis,
                                                       self.scf_drv.scf_tensors)
                    exc_en = rsp_results['eigenvalues'][self.n_state_deriv]
                    e_minus2 = self.scf_drv.get_scf_energy() + exc_en

                    coords[i, d] += 2.0 * self.delta_h
                    # f'(x) ~ [ f(x - 2h) - 8 f(x - h) + 8 f(x + h) - f(x + 2h) ] / ( 12h )
                    self.gradient[i, d] = (e_minus2 - 8.0 * e_minus1
                                           + 8.0 * e_plus1 - e_plus2) / (12.0 * self.delta_h)


        self.ostream.print_blank()

        self.scf_drv.compute(molecule, ao_basis, min_basis)
        self.scf_drv.ostream.state = scf_ostream_state

        # print gradient
        self.print_geometry(molecule)
        self.print_gradient(molecule, labels)

        valstr = '*** Time spent in gradient calculation: '
        valstr += '{:.2f} sec ***'.format(tm.time() - start_time)
        self.ostream.print_header(valstr)
        self.ostream.print_blank()
        self.ostream.flush()

    def compute_numerical_dipole(self, molecule, ao_basis, scf_drv,
                                 rsp_drv, field_strength=1e-5, min_basis=None):
        """
        Performs calculation of numerical dipole moment at RPA or TDA level.

        :param molecule:
            The molecule.
        :param ao_basis:
            The AO basis set.
        :param scf_drv:
            The SCF driver.
        :param rsp_drv:
            The response (RPA or TDA) driver.
        :param field_strength:
            The strength of the external electric field.
        :param min_basis:
            The minimal AO basis set.
        """

        # self.print_header()
        start_time = tm.time()

        self.scf_drv = scf_drv
        scf_ostream_state = self.scf_drv.ostream.state
        self.scf_drv.ostream.state = False

        # This does not have any influence, but depends
        # on the ostream state of the scf driver
        self.rsp_drv = rsp_drv
        # rsp_ostream_state = self.rsp_drv.ostream.state
        # self.rsp_drv.ostream.state = False

        # numerical gradient
        dipole_moment = np.zeros((3))
        field = [0.0, 0.0, 0.0]

        for i in range(3):
            field[i] = field_strength
            self.scf_drv.electric_field = field
            self.scf_drv.compute(molecule, ao_basis, min_basis)
            scf_tensors = self.scf_drv.scf_tensors
            self.rsp_drv.is_converged = False  # only needed for RPA
            rsp_results = self.rsp_drv.compute(molecule, ao_basis,
                                               scf_tensors)
            exc_en_plus = rsp_results['eigenvalues'][self.n_state_deriv]
            e_plus = self.scf_drv.get_scf_energy() + exc_en_plus

            field[i] = -field_strength
            self.scf_drv.compute(molecule, ao_basis, min_basis)
            self.rsp_drv.is_converged = False
            rsp_results = self.rsp_drv.compute(molecule, ao_basis,
                                               self.scf_drv.scf_tensors)
            exc_en_minus = rsp_results['eigenvalues'][self.n_state_deriv]
            e_minus = self.scf_drv.get_scf_energy() + exc_en_minus

            field[i] = 0.0
            dipole_moment[i] = - (e_plus - e_minus) / (2.0 * field_strength)

        return dipole_moment

    def print_geometry(self, molecule):
        """
        Prints the geometry.

        :param molecule:
            The molecule.
        """

        self.ostream.print_block(molecule.get_string())

    def print_gradient(self, molecule, labels):
        """
        Prints the gradient.

        :param molecule:
            The molecule.
        :param labels:
            The atom labels.
        """

        title = 'Gradient (Hartree/Bohr)'
        self.ostream.print_header(title)
        self.ostream.print_header('-' * (len(title) + 2))
        self.ostream.print_blank()

        valstr = '  Atom '
        valstr += '{:>20s}  '.format('Gradient X')
        valstr += '{:>20s}  '.format('Gradient Y')
        valstr += '{:>20s}  '.format('Gradient Z')
        self.ostream.print_header(valstr)
        self.ostream.print_blank()

        for i in range(molecule.number_of_atoms()):
            valstr = '  {:<4s}'.format(labels[i])
            for d in range(3):
                valstr += '{:22.12f}'.format(self.gradient[i, d])
            self.ostream.print_header(valstr)

        self.ostream.print_blank()
        self.ostream.flush()
