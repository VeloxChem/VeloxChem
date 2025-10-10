#
#                                   VELOXCHEM
#              ----------------------------------------------------
#                          An Electronic Structure Code
#
#  SPDX-License-Identifier: BSD-3-Clause
#
#  Copyright 2018-2025 VeloxChem developers
#
#  Redistribution and use in source and binary forms, with or without modification,
#  are permitted provided that the following conditions are met:
#
#  1. Redistributions of source code must retain the above copyright notice, this
#     list of conditions and the following disclaimer.
#  2. Redistributions in binary form must reproduce the above copyright notice,
#     this list of conditions and the following disclaimer in the documentation
#     and/or other materials provided with the distribution.
#  3. Neither the name of the copyright holder nor the names of its contributors
#     may be used to endorse or promote products derived from this software without
#     specific prior written permission.
#
#  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
#  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
#  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
#  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
#  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
#  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
#  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
#  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
#  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
#  OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

from pathlib import Path
from io import StringIO
from contextlib import redirect_stdout
import numpy as np
import json

from .veloxchemlib import bohr_in_angstrom
from .oneeints import (compute_nuclear_potential_integrals,
                       compute_nuclear_potential_gradient_bfs,
                       compute_electrostatic_potential_hessian,
                       compute_electrostatic_integrals_gradient,
                       compute_electric_field_integrals,
                       compute_electric_field_values,
                       compute_electric_field_potential_gradient,
                       compute_electric_field_fock_gradient,
                       compute_electric_field_potential_gradient_for_mm,
                       compute_electric_field_potential_hessian)
from .errorhandler import assert_msg_critical
from .sanitychecks import write_pe_jsonfile

try:
    from pyframe.embedding import (read_input, electrostatic_interactions,
                                   induction_interactions,
                                   repulsion_interactions,
                                   dispersion_interactions)
    from pyframe.embedding.subsystem import ClassicalSubsystem, QuantumSubsystem
    from pyframe.simulation_box import SimulationBox
except ImportError:
    raise ImportError('Unable to import PyFraME. Please install PyFraME.')


def show_embedding(molecule,
                   potfile=None,
                   embedding=None,
                   width=600,
                   height=500,
                   mm_opacity=0.8):

    # Note: this method can involve writting to file so it should only be used
    # on single MPI rank.

    assert_msg_critical(
        not ((embedding is None) and (potfile is None)),
        'show_embedding: expecting embedding or potfile')

    assert_msg_critical(
        not ((embedding is not None) and (potfile is not None)),
        'show_embedding: expecting either embedding or potfile, not both')

    if embedding is not None:
        json_potfile = embedding['inputs']['json_file']

    elif potfile is not None:
        if Path(potfile).suffix == '.json':
            json_potfile = potfile
        else:
            json_potfile = write_pe_jsonfile(molecule, potfile)

    with open(json_potfile, 'r') as fh:
        emb_dict = json.load(fh)

    mm_atoms = []
    for subsys in emb_dict['classical_subsystems']:
        for frag in subsys['classical_fragments']:
            for atom in frag['atoms']:
                mm_atoms.append([atom['element']] + [str(x * bohr_in_angstrom()) for x in atom['coordinate']])
    mm_atoms_xyz = f'{len(mm_atoms)}\n\n'
    for a in mm_atoms:
        mm_atoms_xyz += ' '.join(a) + '\n'

    try:
        import py3Dmol
        viewer = py3Dmol.view(width=width, height=height)
        viewer.setViewStyle({"style": "outline", "width": 0.02})

        viewer.addModel(molecule.get_xyz_string())
        viewer.setStyle({'model': 0}, {"stick": {}, "sphere": {"scale": 0.25}})

        viewer.addModel(mm_atoms_xyz)
        viewer.setStyle(
            {
                'model': 1,
            },
            {
                "stick": {"radius": 0.1, "opacity": mm_opacity},
                "sphere": {"scale": 0.08, "opacity": mm_opacity},
            },
        )

        viewer.zoomTo()
        viewer.show()

    except ImportError:
        raise ImportError('Unable to import py3Dmol')


class EmbeddingIntegralDriver:
    """
    Implements a driver for computing embedding-related integrals.

    :param molecule:
        The molecule.
    :param ao_basis:
        The AO basis set.

    Instance variables:
        - molecule: The molecule.
        - basis: The AO basis set.
    """

    def __init__(self, molecule, ao_basis):
        """
        Initializes the embedding integral driver.

        :param molecule:
            The molecule.
        :param ao_basis:
            The AO basis set.
        """

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
            Electronic fields.
                Shape: (number of atoms, 3)
                Dtype: np.float64.
        """

        return compute_electric_field_values(molecule=self.molecule,
                                             basis=self.basis,
                                             dipole_coords=coordinates,
                                             density=density_matrix)

    def electronic_electrostatic_energy_gradients(self,
            multipole_coordinates: np.ndarray,
            multipole_orders: np.ndarray,
            multipoles: list[np.ndarray],
            density_matrix: np.ndarray) -> np.ndarray:
        """Calculate the electronic electrostatic energy gradients.

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
            density_matrix: Density Matrix that is the source of the electronic field.
                Shape: (number of ao functions, number of ao functions)
                Dtype: np.float64

        Returns:
            Electronic electrostatic energy gradients.
                Shape: (number of nuclei, 3)
                Dtype: np.float64
        """
        op = 0
        # 0 order
        idx = np.where(multipole_orders >= 0)[0]
        charge_coordinates = multipole_coordinates[idx]
        charges = np.array([multipoles[i][0] for i in idx])

        op += compute_nuclear_potential_gradient_bfs(
            molecule=self.molecule,
            basis=self.basis,
            coordinates=charge_coordinates,
            charges=charges,
            density=density_matrix)
        return op


    def electronic_induction_energy_gradients(self,
                                             induced_dipoles:np.ndarray,
                                             coordinates: np.ndarray,
                                             density_matrix: np.ndarray) -> np.ndarray:
        """Calculate the electronic induction energy gradients.

        Args:
            induced_dipoles: Induced dipoles
                Shape (number of induced dipoles, 3)
                Dtype: np.float64
            coordinates: Coordinates on which the fields are to be evaluated.
                Shape: (number of atoms, 3)
                Dtype: np.float64
            density_matrix: Density Matrix that is the source of the electronic field.
                Shape: (number of ao functions, number of ao functions)
                Dtype: np.float64

        Returns:
            Electronic induction energy gradients.
                Shape: (number of nuclei, 3)
                Dtype: np.float64
        """

        return compute_electric_field_potential_gradient(molecule=self.molecule,
                                                         basis=self.basis,
                                                         dipole_coords=coordinates,
                                                         dipole_moments=induced_dipoles,
                                                         density=density_matrix)

    def electronic_induction_fock_gradient(self,
                                           induced_dipoles:np.ndarray,
                                           coordinates: np.ndarray,
                                           i: int) -> np.ndarray:
        """Calculate the electronic induction Fock gradient.

        Args:
            induced_dipoles: Induced dipoles
                Shape (number of induced dipoles, 3)
                Dtype: np.float64
            coordinates: Coordinates on which the fields are to be evaluated.
                Shape: (number of atoms, 3)
                Dtype: np.float64
            i: Index of nucleus "i".
                Shape: (1)
                Dtype: np.int64

        Returns:
            Electronic induction energy gradients.
                Shape: (number of nuclei, 3)
                Dtype: np.float64
        """
        return compute_electric_field_fock_gradient(molecule=self.molecule,
                                                    basis=self.basis,
                                                    dipole_coords=coordinates,
                                                    dipole_moments=induced_dipoles,
                                                    qm_atom_index=i)

    def induced_dipoles_potential_integrals(self,
                                            induced_dipoles: np.ndarray,
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

        return compute_electric_field_integrals(self.molecule,
                                                self.basis,
                                                coordinates,
                                                induced_dipoles)

    def electronic_electrostatic_energy_hessian(self,
                                                multipole_coordinates: np.ndarray,
                                                multipole_orders: np.ndarray,
                                                multipoles: list[np.ndarray],
                                                density_matrix: np.ndarray,
                                                nuc_i: int,
                                                nuc_j: int):
        """Calculate the electronic electrostatic energy Hessian.

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
            density_matrix: Density Matrix that is the source of the electronic field.
                Shape: (number of ao functions, number of ao functions)
                Dtype: np.float64
            nuc_i: Index of nucleus "i" corresponding to the Hessian indexing H_ij.
                Shape: (1)
                Dtype: np.int64
            nuc_j: Index of nucleus "j" corresponding to the Hessian indexing H_ij.
                Shape: (1)
                Dtype: np.int64
        Returns:
            Electronic electrostatic energy Hessian.
                Shape: (3,3)
                Dtype: np.float64
        """
        op = 0
        # 0 order
        idx = np.where(multipole_orders >= 0)[0]
        charge_coordinates = multipole_coordinates[idx]
        charges = np.array([multipoles[i][0] for i in idx])

        op += compute_electrostatic_potential_hessian(molecule=self.molecule,
                                                      basis=self.basis,
                                                      mm_coordinates=charge_coordinates,
                                                      mm_charges=charges,
                                                      density=density_matrix,
                                                      qm_atom_index_i=nuc_i,
                                                      qm_atom_index_j=nuc_j)
        return op

    def electronic_electrostatic_fock_gradient(self,
                                                multipole_coordinates: np.ndarray,
                                                multipole_orders: np.ndarray,
                                                multipoles: list[np.ndarray],
                                                i: int):
        """Calculate the electronic electrostatic fock gradient.

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
            i: Index of nucleus "i".
                Shape: (1)
                Dtype: np.int64

        Returns:
            Electronic induction energy
        """
        op = 0
        # 0 order
        idx = np.where(multipole_orders >= 0)[0]
        charge_coordinates = multipole_coordinates[idx]
        charges = np.array([multipoles[i][0] for i in idx])
        op += compute_electrostatic_integrals_gradient(molecule=self.molecule,
                                                       basis=self.basis,
                                                       mm_charges=charges,
                                                       mm_coordinates=charge_coordinates,
                                                       qm_atom_index_i=i)
        return op

    def compute_electronic_field_gradients(self,
                                           coordinates: np.ndarray,
                                           density_matrix: np.ndarray,
                                           i: int):
        """Calculate the electronic field gradients on coordinates.

        Args:
            coordinates: Coordinates on which the field gradients are to be evaluated.
                Shape: (number of atoms, 3)
                Dtype: np.float64
            density_matrix: Density Matrix that is the source of the electronic field.
                Shape: (number of ao functions, number of ao functions)
                Dtype: np.float64
            i: Index of nucleus "i".
                Shape: (1)
                Dtype: np.int64

        Returns:
            Electronic field gradients.
                Shape: (3, number of atoms, 3)
                Dtype: np.float64
        """
        return compute_electric_field_potential_gradient_for_mm(
            molecule=self.molecule,
            basis=self.basis,
            dipole_coords=coordinates,
            density=density_matrix,
            qm_atom_index=i
        )

    def compute_electronic_field_hessian(self,
                                         coordinates,
                                         induced_dipoles,
                                         density_matrix):
        """Calculate the electronic field Hessian on coordinates.

        Args:
            coordinates: Coordinates on which the field Hessian is to be evaluated.
                Shape: (number of atoms, 3)
                Dtype: np.float64
            induced_dipoles: Induced dipoles at coordinates.
                Shape (number of induced dipoles, 3)
                Dtype: np.float64
            density_matrix: Density Matrix that is the source of the electronic field.
                Shape: (number of ao functions, number of ao functions)
                Dtype: np.float64
        Returns:
            Electronic field Hessian.
                Shape: (3 * number of atoms, 3 * number of atoms)
                Dtype: np.float64
        """
        return compute_electric_field_potential_hessian(molecule=self.molecule,
                                                        basis=self.basis,
                                                        dipole_coords=coordinates,
                                                        dipole_moments=induced_dipoles,
                                                        density=density_matrix)



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
                - mic
                - simulation_box
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
        - simulation_box: Instance of PyFraME SimulationBox.
        - _integral_driver: Instance of EmbeddingIntegralDriver to calculate
                            all relevant integrals.
        - _e_induction: Induction energy.
        - _threshold: Threshold for solving induced dipoles.
        - _max_iterations: Maximum iterations for solving induced dipoles.
        - _solver: Solver method for induced dipoles.
        - _mic: Bool which decides if minimum-image convention (MIC) is used.
    """

    def __init__(self, molecule, ao_basis, options, comm=None, log_level=20):
        """
        Initializes the polarizable embedding.

        :param molecule:
            The molecule to be embedded.
        :param ao_basis:
            The AO basis set.
        :param options:
            Dictionary with settings and inputs.
        :param comm:
            The MPI communicator.
        :param log_level:
            Sets the logging level to be used.
        """

        self.options = options

        assert_msg_critical(
            self.options['settings']['embedding_method'] == 'PE',
            'PolarizableEmbedding: Invalid embedding_method. Only PE is supported.'
        )

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
        self._mic = induced_dipoles_options.get('mic', False)
        # FIXME temporary update coords
        self.quantum_subsystem.coordinates = molecule.get_coordinates_in_bohr()

    def _create_pyframe_objects(self):
        """
        Creates PyFraME objects from input options.

        Reads the input file or uses provided objects to create the simulation box,
        quantum subsystem, and classical subsystem instances.
        """
        def categorize_subsystems(reader_output):
            """Categorize components from the reader output tuple."""
            simulation_box = SimulationBox()
            quantum_subsystems = []
            classical_subsystems = []

            for item in reader_output:
                if isinstance(item, SimulationBox):
                    simulation_box = item
                elif isinstance(item, QuantumSubsystem):
                    quantum_subsystems.append(item)
                elif isinstance(item, ClassicalSubsystem):
                    classical_subsystems.append(item)
                else:
                    raise ValueError(f"Unexpected object type returned by reader: {type(item)}")

            return simulation_box, quantum_subsystems, classical_subsystems

        if "json_file" in self.options['inputs']:
            with redirect_stdout(StringIO()) as fg_out:
                reader_output = read_input.reader(
                    input_data=self.options['inputs']['json_file'],
                    comm=self.comm
                )
                simulation_box, quantum_subsystems, classical_subsystems = categorize_subsystems(reader_output)

                # Handle single or multiple quantum/classical subsystems
                self.simulation_box = simulation_box
                self.quantum_subsystem = quantum_subsystems if len(quantum_subsystems) > 1 else \
                quantum_subsystems[0] if quantum_subsystems else None
                self.classical_subsystem = classical_subsystems if len(classical_subsystems) > 1 else \
                classical_subsystems[0] if classical_subsystems else None

        elif "objects" in self.options['inputs']:
            # Directly use the provided objects
            self.simulation_box = self.options['inputs']['objects']['simulation_box']
            self.quantum_subsystem = self.options['inputs']['objects']['quantum_subsystem']
            self.classical_subsystem = self.options['inputs']['objects']['classical_subsystem']


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
        - _e_vdw: VDW energy.
        - _environment_energy: Interaction energy between multipoles in the
                               classical_subsystem.
    """

    def __init__(self, molecule, ao_basis, options, comm=None, log_level=20):
        """
        Initializes the SCF driver with polarizable embedding.

        :param molecule:
            The molecule to be embedded.
        :param ao_basis:
            The AO basis set.
        :param options:
            Dictionary with settings and inputs.
        :param comm:
            The MPI communicator.
        :param log_level:
            Sets the logging level to be used.
        """

        super().__init__(molecule, ao_basis, options, comm, log_level)

        self._f_elec_es = electrostatic_interactions.es_fock_matrix_contributions(
            classical_subsystem=self.classical_subsystem,
            integral_driver=self._integral_driver)

        self._e_nuc_es = electrostatic_interactions.compute_electrostatic_nuclear_energy(
            quantum_subsystem=self.quantum_subsystem,
            classical_subsystem=self.classical_subsystem)

        vdw_options = self.options['settings'].get('vdw', {})
        self.vdw_method = vdw_options.get('method', 'LJ')
        self.vdw_combination_rule = vdw_options.get('combination_rule',
                                                    'Lorentz-Berthelot')

        if 'vdw' in self.options['settings']:
            self._e_vdw = repulsion_interactions.compute_repulsion_interactions(
                quantum_subsystem=self.quantum_subsystem,
                classical_subsystem=self.classical_subsystem,
                method=self.vdw_method,
                combination_rule=self.vdw_combination_rule)

            self._e_vdw += dispersion_interactions.compute_dispersion_interactions(
                quantum_subsystem=self.quantum_subsystem,
                classical_subsystem=self.classical_subsystem,
                method=self.vdw_method,
                combination_rule=self.vdw_combination_rule)

        else:
            self._e_vdw = 0.0

        self._environment_energy = (self.classical_subsystem.
                                    environment_energy(vdw_method=self.vdw_method,
                                                       vdw_combination_rule=self.vdw_combination_rule)
                                    if self.options['settings'].get(
                                        'environment_energy', False) else 0.0)

        self.pe_summary = None

    def compute_pe_contributions(self, density_matrix):
        """
        Computes polarizable embedding contributions for SCF.

        :param density_matrix:
            The density matrix.

        :return:
            A tuple containing the embedding energy and the embedding Fock matrix.
        """

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
            solver=self._solver,
            mic=self._mic,
            box=self.simulation_box.box)

        self._e_induction = induction_interactions.compute_induction_energy(
            induced_dipoles=self.classical_subsystem.induced_dipoles.
            induced_dipoles,
            total_fields=el_fields + nuc_fields +
            self.classical_subsystem.multipole_fields)

        f_elec_ind = induction_interactions.ind_fock_matrix_contributions(
            classical_subsystem=self.classical_subsystem,
            integral_driver=self._integral_driver)

        e_emb = self._e_induction + self._e_nuc_es + e_elec_es + self._e_vdw

        V_emb = self._f_elec_es - f_elec_ind

        self.pe_summary = []

        self.pe_summary.append(
            'Polarizable Embedding (PE) Energy Contributions')
        self.pe_summary.append('------------------------------------')

        self.pe_summary.append('Electrostatic Contribution         :' +
                               f'{self._e_nuc_es + e_elec_es:20.10f} a.u.')

        self.pe_summary.append('Induced Contribution               :' +
                               f'{self._e_induction:20.10f} a.u.')

        if self._e_vdw != 0.0:
            self.pe_summary.append('Van der Waals Contribution         :' +
                                   f'{self._e_vdw:20.10f} a.u.')

        if self._environment_energy != 0.0:
            self.pe_summary.append('Environment Contribution           :' +
                                   f'{self._environment_energy:20.10f} a.u.')

        self.pe_summary.append('------------------------------------')

        return e_emb, V_emb

    def get_pe_summary(self):
        """
        Gets the summary of polarizable embedding energy contributions.

        :return:
            A list containing the summary information.
        """

        return self.pe_summary


class PolarizableEmbeddingLRS(PolarizableEmbedding):
    """
    Subclass for Linear Response Solver (LRS) calculations with polarizable embedding.
    """

    def __init__(self, molecule, ao_basis, options, comm=None, log_level=20):
        """
        Initializes the LRS driver with polarizable embedding.

        :param molecule:
            The molecule to be embedded.
        :param ao_basis:
            The AO basis set.
        :param options:
            Dictionary with settings and inputs.
        :param comm:
            The MPI communicator.
        :param log_level:
            Sets the logging level to be used.
        """

        super().__init__(molecule, ao_basis, options, comm, log_level)

    def compute_pe_contributions(self, density_matrix):
        """
        Computes polarizable embedding contributions for linear response calculations.

        :param density_matrix:
            The density matrix.

        :return:
            The embedding Fock matrix contribution.
        """

        el_fields = self.quantum_subsystem.compute_electronic_fields(
            coordinates=self.classical_subsystem.coordinates,
            density_matrix=density_matrix,
            integral_driver=self._integral_driver)

        self.classical_subsystem.solve_induced_dipoles(
            external_fields=el_fields,
            threshold=self._threshold,
            max_iterations=self._max_iterations,
            solver=self._solver,
            exclude_static_internal_fields=True,
            mic=self._mic,
            box=self.simulation_box.box)

        f_elec_ind = induction_interactions.ind_fock_matrix_contributions(
            classical_subsystem=self.classical_subsystem,
            integral_driver=self._integral_driver)

        return -1.0 * f_elec_ind

class PolarizableEmbeddingGrad(PolarizableEmbedding):
    """
    Subclass for Gradient (Grad) calculations utilizing polarizable embedding.
    """

    def __init__(self, molecule, ao_basis, options, comm=None, log_level=20):
        """
        Initializes the gradient driver with polarizable embedding.

        :param molecule:
            The molecule to be embedded.
        :param ao_basis:
            The AO basis set.
        :param options:
            Dictionary with settings and inputs.
        :param comm:
            The MPI communicator.
        :param log_level:
            Sets the logging level to be used.
        """
        super().__init__(molecule, ao_basis, options, comm, log_level)
        # TODO read in from options gradient specific options?
        self._e_es_nuc_grad = electrostatic_interactions.compute_electrostatic_nuclear_gradients(
            quantum_subsystem=self.quantum_subsystem,
            classical_subsystem=self.classical_subsystem)
        self._nuc_field_grad = self.quantum_subsystem.compute_nuclear_field_gradients(
            coordinates=self.classical_subsystem.coordinates)
        vdw_options = self.options['settings'].get('vdw', {})
        self.vdw_method = vdw_options.get('method', 'LJ')
        self.vdw_combination_rule = vdw_options.get('combination_rule',
                                                    'Lorentz-Berthelot')

        if 'vdw' in self.options['settings']:
            self._e_vdw_grad = repulsion_interactions.compute_repulsion_interactions_gradient(
                quantum_subsystem=self.quantum_subsystem,
                classical_subsystem=self.classical_subsystem,
                method=self.vdw_method,
                combination_rule=self.vdw_combination_rule)

            self._e_vdw_grad += dispersion_interactions.compute_dispersion_interactions_gradient(
                quantum_subsystem=self.quantum_subsystem,
                classical_subsystem=self.classical_subsystem,
                method=self.vdw_method,
                combination_rule=self.vdw_combination_rule)
        else:
            self._e_vdw_grad = 0.0

    def compute_pe_contributions(self, density_matrix):
        """
        Computes polarizable embedding contributions to the gradient.

        :param density_matrix:
            The density matrix.

        :return:
            The gradient contribution from polarizable embedding.
        """
        # FIXME -> ind dipoles only necessary if not passed from scf results
        # usecase call gradients several times? -> yes different geometries though -> guess of recalculation of ind dipoles from previous set of ind dipoles -> DIIS scheme
        if np.all(self.classical_subsystem.induced_dipoles.induced_dipoles == 0):
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
                solver=self._solver,
                mic=self._mic,
                box=self.simulation_box.box)
        e_es_elec_grad = electrostatic_interactions.compute_electronic_electrostatic_energy_gradients(
            density_matrix=density_matrix,
            classical_subsystem=self.classical_subsystem,
            integral_driver=self._integral_driver)
        e_ind_el_grad = induction_interactions.compute_electronic_induction_energy_gradients(
            density_matrix=density_matrix,
            classical_subsystem=self.classical_subsystem,
            integral_driver=self._integral_driver)
        e_ind_nuc_grad = induction_interactions.compute_induction_energy_gradient(
            induced_dipoles=self.classical_subsystem.induced_dipoles.induced_dipoles,
            total_field_gradients=self._nuc_field_grad)

        return self._e_es_nuc_grad + e_ind_nuc_grad + e_ind_el_grad + self._e_vdw_grad + e_es_elec_grad


class PolarizableEmbeddingHess(PolarizableEmbedding):
    """
    Subclass for Hessian (Hess) calculations utilizing polarizable embedding.
    """

    def __init__(self, molecule, ao_basis, options, density, comm=None, log_level=20):
        """
        Initializes the Hessian driver with polarizable embedding.

        :param molecule:
            The molecule to be embedded.
        :param ao_basis:
            The AO basis set.
        :param options:
            Dictionary with settings and inputs.
        :param density:
            The density matrix.
        :param comm:
            The MPI communicator.
        :param log_level:
            Sets the logging level to be used.
        """
        super().__init__(molecule, ao_basis, options, comm, log_level)
        self._e_es_nuc_hess = electrostatic_interactions.compute_electrostatic_nuclear_hessian(
            quantum_subsystem=self.quantum_subsystem,
            classical_subsystem=self.classical_subsystem)
        self._nuc_field_hess = self.quantum_subsystem.compute_nuclear_field_gradients(
            coordinates=self.classical_subsystem.coordinates)
        vdw_options = self.options['settings'].get('vdw', {})
        self.vdw_method = vdw_options.get('method', 'LJ')
        self.vdw_combination_rule = vdw_options.get('combination_rule',
                                                    'Lorentz-Berthelot')

        if 'vdw' in self.options['settings']:
            self._e_vdw_hess = repulsion_interactions.compute_repulsion_interactions_hessian(
                quantum_subsystem=self.quantum_subsystem,
                classical_subsystem=self.classical_subsystem,
                method=self.vdw_method,
                combination_rule=self.vdw_combination_rule)

            self._e_vdw_hess += dispersion_interactions.compute_dispersion_interactions_hessian(
                quantum_subsystem=self.quantum_subsystem,
                classical_subsystem=self.classical_subsystem,
                method=self.vdw_method,
                combination_rule=self.vdw_combination_rule)
        else:
            self._e_vdw_hess = 0.0
        if np.all(self.classical_subsystem.induced_dipoles.induced_dipoles == 0):
            el_fields = self.quantum_subsystem.compute_electronic_fields(
                coordinates=self.classical_subsystem.coordinates,
                density_matrix=density,
                integral_driver=self._integral_driver)

            nuc_fields = self.quantum_subsystem.compute_nuclear_fields(
                self.classical_subsystem.coordinates)

            self.classical_subsystem.solve_induced_dipoles(
                external_fields=(el_fields + nuc_fields),
                threshold=self._threshold,
                max_iterations=self._max_iterations,
                solver=self._solver,
                mic=self._mic,
                box=self.simulation_box.box)
    def compute_pe_fock_gradient_contributions(self, i):
        """
        Computes polarizable embedding Fock gradient contributions.

        :param i:
            The index of the nucleus.

        :return:
            The Fock gradient contribution from polarizable embedding.
        """

        es_fock_grad = electrostatic_interactions.compute_electronic_electrostatic_fock_gradient(
            i=i,
            classical_subsystem=self.classical_subsystem,
            integral_driver=self._integral_driver
        )

        ind_fock_grad = induction_interactions.compute_electronic_induction_fock_gradient(
            i=i,
            classical_subsystem=self.classical_subsystem,
            integral_driver=self._integral_driver
        )
        return es_fock_grad + ind_fock_grad

    def compute_pe_energy_hess_contributions(self, density_matrix):
        """
        Computes polarizable embedding energy Hessian contributions.

        :param density_matrix:
            The density matrix.

        :return:
            The energy Hessian contribution from polarizable embedding.
        """
        nuc_list = np.arange(self.quantum_subsystem.num_nuclei, dtype=np.int64)
        # TODO: double check density_matrix
        e_es_elec_hess = electrostatic_interactions.compute_electronic_electrostatic_energy_hessian(
            nuc_list=nuc_list,
            density_matrix=density_matrix,
            classical_subsystem=self.classical_subsystem,
            integral_driver=self._integral_driver
            )
        e_es_nuc_hess = electrostatic_interactions.compute_electrostatic_nuclear_hessian(
            quantum_subsystem=self.quantum_subsystem,
            classical_subsystem= self.classical_subsystem)
        # TODO: double check density_matrix
        e_ind_hess = induction_interactions.compute_induction_energy_hessian(
            density_matrix=(-1.0) * density_matrix,
            classical_subsystem=self.classical_subsystem,
            quantum_subsystem=self.quantum_subsystem,
            integral_driver=self._integral_driver,
            threshold=self._threshold,
            max_iterations=self._max_iterations,
            mic=self._mic,
            box=self.simulation_box.box,
            solver=self._solver
            )
        return e_es_elec_hess + e_es_nuc_hess + e_ind_hess + self._e_vdw_hess
