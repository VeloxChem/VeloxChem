from io import StringIO
from contextlib import redirect_stdout
import numpy as np

from .veloxchemlib import (compute_electric_field_integrals,
                           compute_electric_field_values)
from .oneeints import compute_nuclear_potential_integrals

try:
    from pyframe.embedding import (read_input, electrostatic_interactions,
                                   induction_interactions,
                                   repulsion_interactions,
                                   dispersion_interactions)
    from pyframe.embedding.logging_util import log_manager
except ImportError:
    raise ImportError('Unable to import PyFraME. Please install PyFraME.')


class EmbeddingIntegralDriver:

    def __init__(self, molecule, ao_basis):

        self.molecule = molecule
        self.basis = ao_basis

    def multipole_potential_integrals(
            self, multipole_coordinates: np.ndarray,
            multipole_orders: np.ndarray,
            multipoles: list[np.ndarray]) -> np.ndarray:
        """Calculate the electronic potential integrals and multiply with the
        multipoles.

        Args:
            multipole_coordinates: Coordinates of the Multipoles.
                Shape: (number of atoms, 3)
                Dtype: np.float64
            multipole_orders: Multipole orders of all multipoles.
                Shape: (number of atoms)
                Dtype: np.int64
            multipoles: Multipoles multiplied with degeneracy coefficients and
                        taylor coefficients.
                Shape: (number of atoms, number of multipole elements)
                Dtype: np.float64

        Returns:
            Product of electronic potential integrals and multipoles.
                Shape: (number of ao functions, number of ao functions)
                Dtype: np.float64
        """

        # if np.any(multipole_orders > 2):
        #     raise NotImplementedError("""Multipole potential integrals not
        #                                      implemented for order > 2.""")

        op = 0

        # 0 order
        idx = np.where(multipole_orders >= 0)[0]
        charge_coordinates = multipole_coordinates[idx]
        charges = np.array([multipoles[i][0] for i in idx])
        op += compute_nuclear_potential_integrals(
            molecule=self.molecule,
            basis=self.basis,
            coordinates=charge_coordinates,
            charges=charges)

        # 1 order
        if np.any(multipole_orders >= 1):
            idx = np.where(multipole_orders >= 1)[0]
            dipole_coordinates = multipole_coordinates[idx]
            dipoles = np.array([multipoles[i][1:4] for i in idx])
            op += compute_electric_field_integrals(self.molecule, self.basis,
                                                   dipole_coordinates, dipoles)

        return op

    def electronic_fields(self, coordinates: np.ndarray,
                          density_matrix: np.ndarray) -> np.ndarray:
        """Calculate the electronic fields on coordinates.

        Args:
            coordinates: Coordinates on which the fields are to be evaluated.
                Shape: (number of atoms, 3)
                Dtype: np.float64
            density_matrix: Density Matrix that is the source of the electronic field.
                Shape: (number of ao functions, number of ao functions)
                Dtype: np.float64

        Returns:
            Electronic fields. Shape: (number of atoms, 3) Dtype: np.float64.
        """

        return -1.0 * compute_electric_field_values(self.molecule, self.basis,
                                                    coordinates, density_matrix)

    def induced_dipoles_potential_integrals(
            self, induced_dipoles: np.ndarray,
            coordinates: np.ndarray) -> np.ndarray:
        """Calculate the electronic potential integrals and contract with the
        induced dipoles of Atoms.

        Args:
            induced_dipoles: Induced dipoles
                Shape (number of induced dipoles, 3)
                Dtype: np.float64
            coordinates: Coordinates of the induced dipoles on which the
                         integrals are to be evaluated.
                Shape (number of induced dipoles, 3)
                Dtype: np.float64

        Returns:
            Product of the electronic potential integrals and the induced dipoles.
                Shape: (number of ao functions, number of ao functions)
                Dtype: np.float64
        """

        return compute_electric_field_integrals(self.molecule, self.basis,
                                                coordinates, induced_dipoles)


class PolarizableEmbedding:
    """
    Implements the fragment-based polarizable embedding (PE) method.

    :param molecule: The molecule to be embedded.
    :param ao_basis: The AO basis set.
    :param options: Dictionary with settings and inputs.
        settings:
            - embedding_method: string to set embedding method.
            - vdw: dictionary with keys: method and combination_rule.
                - method
                - combination_rule
            - induced_dipoles: dictionary with keys: threshold, max_iterations,
                               and solver.
                - solver
                - threshold
                - max_iterations
            - environment_energy: bool which decides if the environment energy
                                  will be calculated.
        inputs:
            - json_file: string that is the path to the json file that contains
                         the embedding potentials.
            - objects:
                - quantum_subsystem: PyFraME class subsystem.QuantumSubsystem.
                - classical_subsystem: PyFraME class subsystem.ClassicalSubsystem.
    :param comm: The MPI communicator.
    :param log_level: Sets the logging level to be used.

    Instance variables:
        - comm: The MPI communicator.
        - molecule: The molecule.
        - basis: The AO basis set.
        - options: Dictionary with the options and inputs.
        - quantum_subsystem: Instance of PyFraME QuantumSubsystem.
        - classical_subsystem: Instance of PyFraME ClassicalSubsystem.
        - _integral_driver: Instance of EmbeddingIntegralDriver to calculate
                            all relevant integrals.
        - _e_induction: Induction energy.
        - _threshold: Threshold for solving induced dipoles.
        - _max_iterations: Maximum iterations for solving induced dipoles.
        - _solver: Solver method for induced dipoles.
    """

    def __init__(self, molecule, ao_basis, options, comm=None, log_level=20):

        self.log_manager = log_manager
        self.log_manager.set_level(log_level)

        self.options = options
        if self.options['settings']['embedding_method'] != 'PE':
            raise NotImplementedError(
                "Method set is not Polarizable Embedding.")

        self.comm = comm
        self.molecule = molecule
        self.basis = ao_basis

        self._integral_driver = EmbeddingIntegralDriver(molecule=self.molecule,
                                                        ao_basis=self.basis)
        self._create_pyframe_objects()
        self._e_induction = None

        induced_dipoles_options = self.options['settings'].get(
            'induced_dipoles', {})

        self._threshold = induced_dipoles_options.get('threshold', 1e-8)
        self._max_iterations = induced_dipoles_options.get(
            'max_iterations', 100)
        self._solver = induced_dipoles_options.get('solver', 'jacobi')

    def _create_pyframe_objects(self):

        if "json_file" in self.options['inputs']:
            with redirect_stdout(StringIO()) as fg_out:
                (self.quantum_subsystem,
                 self.classical_subsystem) = read_input.reader(
                     input_data=self.options['inputs']['json_file'],
                     comm=self.comm)

        elif "objects" in self.options['inputs']:
            self.quantum_subsystem, self.classical_subsystem = (
                self.options['inputs']['objects']['quantum_subsystem'],
                self.options['inputs']['objects']['classical_subsystem'])


class PolarizableEmbeddingSCF(PolarizableEmbedding):
    """
    Subclass for Self-Consistent Field (SCF) calculations with polarizable embedding.

    Instance variables:
        - e: Sum of PE contributions (E_rep + E_disp + E_ind + E_es).
        - v: Sum of PE Fock matrix contributions.
        - _e_es: Electrostatic energy.
        - vdw_method: The VDW method to be used.
        - vdw_combination_rule: The VDW combination rule to be used.
        - _f_elec_es: Electrostatic Fock matrix.
        - _e_nuc_es: Electrostatic nuclear energy.
        - _e_rep: VDW repulsion energy.
        - _e_disp: VDW dispersion energy.
        - _environment_energy: Interaction energy between multipoles in the
                               classical_subsystem.
    """

    def __init__(self, molecule, ao_basis, options, comm=None, log_level=20):
        super().__init__(molecule, ao_basis, options, comm, log_level)

        self._f_elec_es = electrostatic_interactions.es_fock_matrix_contributions(
            classical_subsystem=self.classical_subsystem,
            integral_driver=self._integral_driver)

        self._e_nuc_es = electrostatic_interactions.compute_electrostatic_nuclear_energy(
            quantum_subsystem=self.quantum_subsystem,
            classical_subsystem=self.classical_subsystem)

        self._e_es = None
        self.e = None
        self.v = None

        vdw_options = self.options['settings'].get('vdw', {})
        self.vdw_method = vdw_options.get('method', 'LJ')
        self.vdw_combination_rule = vdw_options.get('combination_rule',
                                                    'Lorentz-Berthelot')

        self._e_rep = repulsion_interactions.compute_repulsion_interactions(
            quantum_subsystem=self.quantum_subsystem,
            classical_subsystem=self.classical_subsystem,
            method=self.vdw_method,
            combination_rule=self.vdw_combination_rule
        ) if 'vdw' in self.options['settings'] else 0.0

        self._e_disp = dispersion_interactions.compute_dispersion_interactions(
            quantum_subsystem=self.quantum_subsystem,
            classical_subsystem=self.classical_subsystem,
            method=self.vdw_method,
            combination_rule=self.vdw_combination_rule
        ) if 'vdw' in self.options['settings'] else 0.0

        self._environment_energy = (self.classical_subsystem.environment_energy
                                    if self.options['settings'].get(
                                        'environment_energy', True) else 0.0)

    def compute_pe_contributions(self, density_matrix):
        self.log_manager.reset()

        if self._e_nuc_es is None:
            self._e_nuc_es = electrostatic_interactions.compute_electrostatic_nuclear_energy(
                quantum_subsystem=self.quantum_subsystem,
                classical_subsystem=self.classical_subsystem)

        if self._f_elec_es is None:
            self._f_elec_es = electrostatic_interactions.es_fock_matrix_contributions(
                classical_subsystem=self.classical_subsystem,
                integral_driver=self._integral_driver)

        e_elec_es = np.sum(self._f_elec_es * density_matrix)

        el_fields = self.quantum_subsystem.compute_electronic_fields(
            coordinates=self.classical_subsystem.coordinates,
            density_matrix=density_matrix,
            integral_driver=self._integral_driver)

        nuc_fields = self.quantum_subsystem.compute_nuclear_fields(
            self.classical_subsystem.coordinates)

        self.classical_subsystem.solve_induced_dipoles(
            external_fields=(el_fields + nuc_fields),
            threshold=self._threshold,
            max_iterations=self._max_iterations,
            solver=self._solver)

        self._e_induction = induction_interactions.compute_induction_energy(
            induced_dipoles=self.classical_subsystem.induced_dipoles.
            induced_dipoles,
            total_fields=el_fields + nuc_fields +
            self.classical_subsystem.multipole_fields)

        f_elec_induction = induction_interactions.ind_fock_matrix_contributions(
            classical_subsystem=self.classical_subsystem,
            integral_driver=self._integral_driver)

        self._e_es = self._e_nuc_es + e_elec_es
        self.e = self._e_induction + self._e_es + self._e_disp + self._e_rep
        self.v = self._f_elec_es - f_elec_induction

        self.log_manager.logger.info(
            'Polarizable Embedding (PE) Energy Contributions'.ljust(92))
        self.log_manager.logger.info(
            '------------------------------------'.ljust(92))
        self.log_manager.logger.info(
            f'Electrostatic Contributions (E_es) :{self._e_es:20.10f} a.u.')
        self.log_manager.logger.info(
            f'Induced Contributions (E_ind)      :{self._e_induction:20.10f} a.u.'
        )

        if self._environment_energy != 0.0:
            self.log_manager.logger.info(
                'Environment Contributions (E_mul)   :' +
                f'{self._environment_energy:20.10f} a.u.')

        if self._e_disp != 0.0:
            self.log_manager.logger.info(
                f'Dispersion Contributions (E_disp)   :{self._e_disp:20.10f} a.u.'
            )

        if self._e_rep != 0.0:
            self.log_manager.logger.info(
                f'Repulsion Contributions (E_disp)   :{self._e_rep:20.10f} a.u.'
            )

        self.log_manager.logger.info(
            '------------------------------------'.ljust(92))

        return self.e, self.v


class PolarizableEmbeddingLRS(PolarizableEmbedding):
    """
    Subclass for Linear Response Solver (LRS) calculations with polarizable embedding.
    """

    def __init__(self, molecule, ao_basis, options, comm=None, log_level=20):

        super().__init__(molecule, ao_basis, options, comm, log_level)

    def compute_pe_contributions(self, density_matrix):

        el_fields = self.quantum_subsystem.compute_electronic_fields(
            coordinates=self.classical_subsystem.coordinates,
            density_matrix=density_matrix,
            integral_driver=self._integral_driver)

        self.classical_subsystem.solve_induced_dipoles(
            external_fields=el_fields,
            threshold=self._threshold,
            max_iterations=self._max_iterations,
            solver=self._solver,
            exclude_static_internal_fields=True)

        f_elec_induction = induction_interactions.ind_fock_matrix_contributions(
            classical_subsystem=self.classical_subsystem,
            integral_driver=self._integral_driver)

        return -f_elec_induction
