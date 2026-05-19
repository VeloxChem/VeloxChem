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

from mpi4py import MPI
from pathlib import Path
import numpy as np
import sys
import time

from .veloxchemlib import mpi_master, bohr_in_angstrom
from .mmforcefieldgenerator import MMForceFieldGenerator
from .molecule import Molecule
from .outputstream import OutputStream
from .waterparameters import get_water_parameters
from .errorhandler import assert_msg_critical


class SolvationBuilder:
    """
    Builds systems for molecular dynamics simulations.

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.

    Instance variables:
        - threshold: The threshold for overlap between molecules in angstrom.
        - acceleration_threshold: The threshold of number of molecules to increase the batch size to speed up the process.
        - number_of_attempts: The number of attempts to insert a molecule.
        - random_rotation: Boolean flag to indicate if the solvent molecules will be randomly rotated. False for linear molecules.
        - solute: The VeloxChem molecule object of the solute.
        - solute_labels: The list of atom labels of the solute.
        - solute_ff: The ForceField object of the solute.
        - solvents: The list of VeloxChem molecule objects of the solvent molecules.
        - solvent_labels: The list of atom labels of the solvent molecules.
        - quantities: The list of quantities of each solvent molecule.
        - added_solvent_counts: The list of the number of solvent molecules added to the system.
        - solvent_name: The name of the solvent molecule.
        - system_molecule: The VeloxChem molecule object of the solvated system.
        - box: The dimensions of the box (x, y, z).
        - centered_solute: The list of tuples with the molecule id, atom labels and coordinates of the centered solute.
        - system: The list of tuples with the molecule id, atom labels and coordinates of the system.
        - equilibration_flag: Boolean flag to indicate if an has been requested.
        - temperature: The temperature of the for the equilibration in Kelvin.
        - pressure: The pressure for the equilibration in bar.
        - steps: The number of steps for the equilibration.
        - pcharge: The positive charge for the counterion. Default is 'Na'.
        - ncharge: The negative charge for the counterion. Default is 'Cl'.
        - counterion: The VeloxChem molecule object of the counterion.
        - parent_forcefield: The name of the parent forcefield. Default is 'amber03'.
    """

    def __init__(self, comm=None, ostream=None):
        '''
        Initialize the System Builder class.
        '''

        if comm is None:
            comm = MPI.COMM_WORLD

        assert_msg_critical(
            comm.Get_size() == 1,
            'SolvationBuilder: Only single-rank interactive use is supported. '
            'Run the solvation builder with one MPI rank.')

        if ostream is None:
            if comm.Get_rank() == mpi_master():
                ostream = OutputStream(sys.stdout)
            else:
                ostream = OutputStream(None)

        # Output stream
        self.ostream = ostream
        self.workdir = Path('.')

        # Packing configuration
        self.threshold = 1.8
        self.acceleration_threshold = 1000
        self.number_of_attempts = 100
        self.random_rotation = True
        self.failures_factor = 0.7

        # Molecules
        # Solute
        self.solute = None
        self.solute_labels = []
        self.solute_ff = None
        # Solvents
        self.solvents = []
        self.solvent_labels = []
        self.quantities = []
        self.added_solvent_counts = []
        self.solvent_name = None
        self.water_parameters = get_water_parameters()

        # System
        self.system_molecule = None
        self.box = None
        self.centered_solute = None
        # This is a list of tuples with the molecule id, atom labels and coordinates of the system
        self.system = []
        self.write_pdb_only = False

        # NPT Equilibration options
        self.equilibration_flag = False
        self.temperature = 300
        self.pressure = 1
        self.steps = 5000

        # Neutralization options
        self.pcharge = 'Na'
        self.ncharge = 'Cl'
        self.counterion = None
        self.added_counterions = 0

        # Standard forcefield
        self.parent_forcefield = 'amber03'

    def _path(self, filename):
        """
        Resolve a generated filename relative to the builder work directory.

        :param filename:
            The generated filename.
        """

        return Path(self.workdir) / filename

    def solvate(self,
                solute,
                solvent='cspce',
                solvent_molecule=None,
                padding=1.0,
                target_density=None,
                neutralize=True,
                equilibrate=False,
                equilibration_steps=None,
                box=None):
        """
        Create a solvated system with the most typical solvent molecules.

        :param solute:
            The VeloxChem molecule object of the solute.
        :param solvent:
            The name of the solvent molecule. The default is 'water'.
            Available options: 'cspce', 'ctip3p', 'spce', 'tip3p', 'ethanol', 'methanol', 'acetone',
            'chloroform', 'hexane', 'toluene', 'dcm', 'benzene', 'dmso', 'thf', 
            'acetonitrile', 'dmf', 'other' or 'itself'.
                * 'other': The solvent molecule must be provided.
                * 'itself': The solute molecule is used as the solvent as well.
        :param solvent_molecule:
            The VeloxChem molecule object of the solvent molecule. Required if solvent is 'other'.
        :param padding:
            The padding to be added to the bounding box of the molecule in nm.
        :param target_density:
            The target density of the solvent in kg/m^3. If None, the experimental density is used.
            If the solvent is 'other' the target density must be provided or will be estimated from the solvent molecule.
        :param neutralize:
            If True, neutralizes the total charge of the molecule with a counterion.
        :param equilibrate:
            If True, perform an equilibration of the system.
        :param equilibration_steps:
            Optional number of MD steps to use for the equilibration run.
            If None, the builder default stored in `self.steps` is used.
        """

        from scipy.spatial import cKDTree

        # Save the solvent name
        self.solvent_name = solvent
        self.equilibration_flag = False

        self._print_solvation_header(solvent, padding, equilibrate)

        # Add the solute to the system
        self._load_solute_molecule(solute)
        solute_volume = self._get_volume(solute) * 1e-3

        volume_nm3, box_center = self._prepare_solvation_box(
            solute, solute_volume, padding, box)
        self._reset_solvation_state()
        solvent_molecule, number_of_solvents = self._register_typical_solvent(
            solute, solvent, solvent_molecule, target_density, volume_nm3)

        tree = self._center_solute_and_build_tree(box_center, cKDTree)

        if neutralize:

            charge = int(self.solute.get_charge())
            self.added_counterions = abs(charge)

            if charge == 0:
                self.ostream.print_info(
                    "The solute is neutral, no counterions will be added")
                self.ostream.flush()

            elif charge > 0:
                self.ostream.print_info(
                    f"The solute has a charge of {charge}, adding {abs(charge)}{self.ncharge} counterions"
                )
                self.ostream.flush()

                self.ion_name = self.ncharge
                self.counterion = self._counterion_molecules()
                # Update the required number of solvents
                number_of_solvents = max(number_of_solvents - abs(charge), 0)
                self._set_last_solvent_quantity(number_of_solvents)

            elif charge < 0:
                self.ostream.print_info(
                    f"The solute has a charge of {charge}, adding {abs(charge)} {self.pcharge} counterions"
                )
                self.ostream.flush()

                self.ion_name = self.pcharge
                self.counterion = self._counterion_molecules()
                number_of_solvents = max(number_of_solvents - abs(charge), 0)
                self._set_last_solvent_quantity(number_of_solvents)

        # Solvate the solute with the solvent molecules
        # This dynamic batch size is used to avoid building too often the KDTree.

        max_batch_size = int(number_of_solvents / 2)
        min_batch_size = 10

        # Check linearity of the solvent molecule
        if self.solvent_name == 'itself':
            solvent_molecule = solute
        else:
            solvent_molecule = solvent_molecule

        if self._molecule_linearity(solvent_molecule):
            self.ostream.print_info(
                f"The solvent molecule is linear, random rotation is disabled")
            self.ostream.flush()
            self.random_rotation = False

        # Solvate the solute with the solvent molecules
        start = time.time()
        tree = self._pack_solvent_molecules(cKDTree,
                                            tree,
                                            use_dynamic_batch=True,
                                            max_batch_size=max_batch_size,
                                            min_batch_size=min_batch_size)

        # Insert the counterions
        if self.counterion:
            self._insert_counterions(cKDTree, tree, charge)

        end = time.time()
        self.ostream.print_info(
            f"Time to solvate the system: {end - start:.2f} s")
        self.ostream.flush()
        # Print results
        box_volume_nm3 = self.box[0] * self.box[1] * self.box[2] * 1e-3
        self.ostream.print_info(
            'The density of the system after packing is: ' +
            f'{self._compute_density(box_volume_nm3)} kg/m^3'
        )
        self.ostream.print_blank()
        self.ostream.flush()

        self.system_molecule = self._save_molecule()

        if equilibrate:
            run_steps = (self.steps if equilibration_steps is None else
                         equilibration_steps)
            # TODO: run perform_equilibration using openmm files
            try:
                start = time.time()
                self.perform_equilibration(steps=equilibration_steps)
                end = time.time()
                self.ostream.print_info("Equilibrating the system")
                self.ostream.print_blank()
                self.ostream.print_info(f"Duration: {run_steps/1000} ps")
                self.ostream.print_info(f"Temperature: {self.temperature} K")
                self.ostream.print_info(f"Pressure: {self.pressure} bar")
                self.ostream.print_blank()
                self.ostream.print_info(
                    f"Elapsed time to equilibrate the system: {end - start:.2f} s"
                )
                self.equilibration_flag = True
            except ValueError:
                # ValueError: Could not locate #include file: amber03.ff/forcefield.itp
                self.ostream.print_info(
                    "Equilibration skipped due to missing files")
                self.ostream.print_blank()
            self.ostream.flush()

    def custom_solvate(self, solute, solvents, proportion, box_size):
        """
        Solvate the solute

        :param solute:
            The VeloxChem molecule object of the solute
        :param solvents:
            The list of VeloxChem molecule objects of the solvent molecules
        :param proportion:
            The proportion of solvent molecules in term of number of molecules.
        :param box_size:
            The array with the dimensions of the box (x, y, z)
        """

        from scipy.spatial import cKDTree

        # Reset solvent name
        self.solvent_name = None
        self.equilibration_flag = False

        # Add the solute to the system
        self._load_solute_molecule(solute)
        solute_volume = self._get_volume(solute) * 1e-3

        solute_nonh_count = 0
        for atom_label in solute.get_labels():
            if atom_label != 'H':
                solute_nonh_count += 1

        # Define the box
        self._define_box(*box_size)

        # Estimate quantities
        # normalize proportion
        assert_msg_critical(
            len(solvents) == len(proportion),
            'SolvationBuilder: The number of solvents and proportions must match.'
        )
        assert_msg_critical(
            len(proportion) > 0,
            'SolvationBuilder: At least one solvent proportion must be provided.'
        )
        assert_msg_critical(
            all(p > 0 for p in proportion),
            'SolvationBuilder: All custom solvent proportions must be positive.'
        )
        sum_proportion = sum(proportion)
        normalized_proportion = [p / sum_proportion for p in proportion]
        # Make rough estimation of quantities based on the empirical
        # observation that the density of nonhydrogen atoms is around 30 per
        # nm^3. To give the solvation builder a bit more room for insertion we
        # use an approximate density of 28 nonhydrogen atoms per nm^3.
        norm_nonh_count = 0
        for solvent, norm_prop in zip(solvents, normalized_proportion):
            for atom_label in solvent.get_labels():
                if atom_label != 'H':
                    norm_nonh_count += norm_prop
        assert_msg_critical(
            norm_nonh_count > 0.0,
            'SolvationBuilder: Custom solvents must contain at least one non-hydrogen atom '
            'to estimate packing.')
        box_volume_nm3 = box_size[0] * box_size[1] * box_size[2] * 1e-3
        available_volume_nm3 = box_volume_nm3 - solute_volume
        assert_msg_critical(
            available_volume_nm3 > 0.0,
            'SolvationBuilder: The available solvent volume must be positive. '
            'Increase the box size or padding.')
        max_nonh_count = available_volume_nm3 * 28
        assert_msg_critical(
            max_nonh_count > 0.0,
            'SolvationBuilder: The custom solvation box is too small to fit any solvent molecules. '
            'Increase the box size.')
        total_solvent_count = int(max_nonh_count / norm_nonh_count)
        assert_msg_critical(
            total_solvent_count > 0,
            'SolvationBuilder: The custom solvation box is too small to fit any solvent molecules. '
            'Increase the box size.')
        raw_solvent_counts = np.array(
            normalized_proportion) * total_solvent_count
        solvent_counts = np.floor(raw_solvent_counts).astype(int)
        remaining_count = total_solvent_count - int(solvent_counts.sum())
        if remaining_count > 0:
            fractional_parts = raw_solvent_counts - solvent_counts
            order = np.argsort(-fractional_parts)
            for solvent_idx in order[:remaining_count]:
                solvent_counts[solvent_idx] += 1
        quantities = list(solvent_counts)

        # Load the solvent molecules
        self._clear_system()
        self._clear_solvent_molecules()
        self._clear_counterions()
        for solvent, quantity in zip(solvents, quantities):
            self._load_solvent_molecule(solvent, quantity)

        box_center = 0.5 * np.array(self.box)
        tree = self._center_solute_and_build_tree(box_center, cKDTree)

        # Solvate the solute with the solvent molecules. For mixed solvents we
        # insert them in an interleaved order so the achieved composition stays
        # closer to the requested ratio when the box becomes saturated.
        self._pack_solvent_mixture(cKDTree, tree)
        self.ostream.print_info(
            'The density of the system after packing is: ' +
            f'{self._compute_density(box_volume_nm3)} kg/m^3'
        )
        self.ostream.print_blank()
        self.ostream.flush()

        self.system_molecule = self._save_molecule()

    def _print_solvation_header(self, solvent, padding, equilibrate):
        """
        Print the standard header for the solvation workflow.
        """

        header_msg = "VeloxChem Solvation Builder"
        self.ostream.print_header(header_msg)
        self.ostream.print_header("=" * (len(header_msg) + 2))
        self.ostream.print_blank()
        self.ostream.flush()

        self.ostream.print_info(
            f"Solvating the solute with {solvent} molecules")
        self.ostream.print_info(f"Padding: {padding} nm")
        self.ostream.print_blank()
        self.ostream.flush()

        if equilibrate:
            self.ostream.print_info("NPT Equilibration of the box requested")
            self.ostream.print_blank()
            self.ostream.flush()

    def _prepare_solvation_box(self, solute, solute_volume, padding, box):
        """
        Define the simulation box and report the available solvent volume.
        """

        if box is None:
            box_size = self._determine_cubic_box_size(
                solute.get_coordinates_in_angstrom(), padding)
            volume_nm3 = box_size**3 * 0.001 - solute_volume
            self._define_box(box_size, box_size, box_size)
            self.ostream.print_info(
                "The box size is: {:.2f} x {:.2f} x {:.2f} nm^3".format(
                    box_size * 0.1, box_size * 0.1, box_size * 0.1))
            self.ostream.flush()
            box_center = [box_size / 2] * 3
        else:
            volume_nm3 = box[0] * box[1] * box[2] * 0.001 - solute_volume
            self._define_box(*box)
            box_center = [box[0] / 2, box[1] / 2, box[2] / 2]

        assert_msg_critical(
            volume_nm3 > 0.0,
            'SolvationBuilder: The available solvent volume must be positive. '
            'Increase the box size or padding.')

        self.ostream.print_info(
            "The volume of the solute is: {:.2f} nm^3".format(solute_volume))
        self.ostream.flush()
        self.ostream.print_info(
            "The volume available for the solvent is: {:.2f} nm^3".format(
                volume_nm3))
        self.ostream.flush()

        return volume_nm3, box_center

    def _reset_solvation_state(self):
        """
        Clear the current solvated system and registered solvent metadata.
        """

        self._clear_system()
        self._clear_solvent_molecules()
        self._clear_counterions()

    def _register_typical_solvent(self, solute, solvent, solvent_molecule,
                                  target_density, volume_nm3):
        """
        Register the requested solvent and estimate the solvent count.
        """

        if solvent == 'other':
            if solvent_molecule is None:
                raise ValueError(
                    "The solvent molecule must be provided if the solvent is 'other'"
                )
            if target_density is None:
                raise ValueError(
                    "The target density must be provided if the solvent is 'other'"
                )
            self.ostream.print_info(
                f"The target density of the solvent is: {target_density} kg/m^3"
            )
            self.ostream.flush()
            mols_per_nm3 = self._density_to_mols_per_nm3(
                solvent_molecule, target_density)

        elif solvent == 'itself':
            if target_density is None:
                raise ValueError(
                    "The target density must be provided if the solvent is 'itself'"
                )
            self.ostream.print_info(
                f"The target density of the solvent is: {target_density} kg/m^3"
            )
            self.ostream.flush()
            solvent_molecule = solute
            mols_per_nm3 = self._density_to_mols_per_nm3(
                solvent_molecule, target_density)

        else:
            mols_per_nm3, density, smiles_code = self._solvent_properties(
                solvent)
            self.ostream.print_info(
                f"The experimental density of {solvent} is: {density} kg/m^3")
            self.ostream.flush()
            solvent_molecule = Molecule.read_smiles(smiles_code)

        number_of_solvents = int(mols_per_nm3 * volume_nm3)
        self.ostream.print_info(
            f"The number of solvent molecules to be added to match the target density is: {number_of_solvents}"
        )
        self.ostream.flush()
        self._load_solvent_molecule(solvent_molecule, number_of_solvents)

        return solvent_molecule, number_of_solvents

    def _center_solute_and_build_tree(self, box_center, cKDTree):
        """
        Center the solute in the box and initialize the packing KD-tree.
        """

        solute_xyz = self.solute.get_coordinates_in_angstrom()
        centroid = np.mean(solute_xyz, axis=0)
        translation = np.array(box_center) - np.array(centroid)
        self.centered_solute = [
            (0, label, coord + translation)
            for label, coord in zip(self.solute_labels, solute_xyz)
        ]
        self.system.extend(self.centered_solute)

        existing_coords = np.array([atom[-1] for atom in self.system])
        return cKDTree(existing_coords) if existing_coords.size > 0 else None

    def _pack_solvent_molecules(self,
                                cKDTree,
                                tree,
                                use_dynamic_batch=False,
                                max_batch_size=None,
                                min_batch_size=10):
        """
        Insert the registered solvent molecules into the current system.
        """

        for solvent, quantity in zip(self.solvents, self.quantities):
            added_count = 0
            failure_count = 0
            partial_fill = False

            if use_dynamic_batch and quantity >= self.acceleration_threshold:
                batch_size = self._compute_batch_size(
                    added_count=added_count,
                    max_batch_size=max_batch_size,
                    min_batch_size=min_batch_size,
                    total_quantity=quantity)
                max_failures = self.failures_factor * quantity
            else:
                batch_size = 1
                max_failures = max(1, int(np.ceil(self.failures_factor *
                                                  quantity)))

            while added_count < quantity:
                attempts = min(batch_size, quantity - added_count)
                new_molecules = []

                for _ in range(attempts):
                    result = self._insert_molecule(solvent, tree)
                    if result:
                        new_molecules.extend(result)
                        added_count += 1
                        failure_count = 0
                    else:
                        failure_count += 1

                if failure_count >= max_failures:
                    self._report_partial_fill(quantity, added_count,
                                              failure_count)
                    partial_fill = True
                    break

                if new_molecules:
                    self.system.extend(new_molecules)
                    tree = cKDTree(np.array([atom[-1] for atom in self.system]))

            self.added_solvent_counts.append(added_count)
            if not partial_fill:
                self._report_solvation_result(added_count, quantity)

        return tree

    def _pack_solvent_mixture(self, cKDTree, tree):
        """
        Insert multiple solvent species in an interleaved order.

        Packing one solvent species to completion before attempting the next can
        distort the achieved mixture when the box gets crowded. This routine
        always tries to insert the most underfilled solvent first.
        """

        solvent_count = len(self.solvents)
        added_counts = [0] * solvent_count
        failure_counts = [0] * solvent_count
        partial_fill = [False] * solvent_count
        max_failures = [
            max(1, int(np.ceil(self.failures_factor * quantity)))
            for quantity in self.quantities
        ]
        active_indices = {
            idx for idx, quantity in enumerate(self.quantities) if quantity > 0
        }

        stop_reason = None
        stop_failure_idx = None

        while active_indices:
            ordered_indices = sorted(
                active_indices,
                key=lambda idx: (
                    added_counts[idx] / self.quantities[idx],
                    -self.quantities[idx],
                    idx,
                ))
            progress_made = False

            idx = ordered_indices[0]
            result = self._insert_molecule(self.solvents[idx], tree)
            if result:
                self.system.extend(result)
                tree = cKDTree(np.array([atom[-1] for atom in self.system]))
                added_counts[idx] += 1
                failure_counts[idx] = 0
                progress_made = True

                if added_counts[idx] >= self.quantities[idx]:
                    active_indices.remove(idx)
            else:
                failure_counts[idx] += 1
                if failure_counts[idx] >= max_failures[idx]:
                    stop_reason = (
                        "mixed-solvent packing stopped to preserve the "
                        "requested composition as the box became saturated")
                    stop_failure_idx = idx
                    for active_idx in active_indices:
                        if added_counts[active_idx] < self.quantities[active_idx]:
                            partial_fill[active_idx] = True
                    active_indices.clear()

            if not progress_made and not active_indices:
                break

        self.added_solvent_counts.extend(added_counts)
        for idx, (added_count,
                  quantity) in enumerate(zip(added_counts, self.quantities)):
            if partial_fill[idx]:
                failure_count = (
                    failure_counts[idx] if idx == stop_failure_idx else None)
                self._report_partial_fill(quantity,
                                          added_count,
                                          failure_count=failure_count,
                                          reason=stop_reason)
            else:
                self._report_solvation_result(added_count, quantity)

        return tree

    def _insert_counterions(self, cKDTree, tree, charge):
        """
        Insert counterions and update the recorded count if packing is partial.
        """

        inserted_counterions = 0
        for _ in range(abs(charge)):
            result = self._insert_molecule(self.counterion, tree)
            if result:
                self.system.extend(result)
                tree = cKDTree(np.array([atom[-1] for atom in self.system]))
                inserted_counterions += 1

        if inserted_counterions != self.added_counterions:
            self._report_partial_fill(self.added_counterions,
                                      inserted_counterions,
                                      self.added_counterions -
                                      inserted_counterions)
            self.added_counterions = inserted_counterions

    def write_gromacs_files(self,
                            solute_ff=None,
                            solvent_ffs=None,
                            equilibration=False):
        """
        Write the GROMACS topology and coordinate files for the system.

        :param solute_ff:
            The force-field object of the solute.
        :param solvent_ffs:
            The force-field objects of the solvent molecules.
        :param equilibration:
            If True, suppress informational output for internal equilibration use.
        """

        self._generate_forcefields(solute_ff, solvent_ffs, equilibration)

        # Special case for 'itself' solvent
        if self.solvent_name == 'itself':
            self._write_self_solvent_gromacs_files(equilibration)
            self._write_system_gro(filename='liquid.gro')

        else:
            self._write_component_gromacs_files(equilibration)
            atomtypes = self._collect_component_atomtypes()
            solute_atomtypes_lines = self._strip_component_atomtypes()
            self._inject_solute_atomtypes_into_top(solute_atomtypes_lines)
            self._write_system_top_file(atomtypes)

            if not equilibration:
                self.ostream.print_info("system.top file written")
                self.ostream.flush()

            # Write the system GRO file
            self._write_system_gro()

            if not equilibration:
                self.ostream.print_info("system.gro file written")
                self.ostream.flush()

    def _write_self_solvent_gromacs_files(self, equilibration):
        """
        Write GROMACS files for the pure-liquid special case.
        """

        liquid_itp = self._path('liquid.itp')
        liquid_top = self._path('liquid.top')
        self.solute_ff.write_itp(str(liquid_itp), 'MOL')
        self.solute_ff.write_top(str(liquid_top), str(liquid_itp), 'MOL')

        with open(liquid_top, 'r') as f:
            lines = f.readlines()
        with open(liquid_top, 'w') as f:
            for line in lines:
                if line.startswith('; Compound '):
                    f.write(line)
                    f.write(f'MOL               {self.added_solvent_counts[0] + 1}\n')
                    break
                f.write(line)

        if not equilibration:
            self.ostream.print_info(
                "liquid.itp, liquid.top, and solute.top files written")
            self.ostream.flush()

    def _write_component_gromacs_files(self, equilibration):
        """
        Write the component GROMACS files used to assemble the system topology.
        """

        solute_itp = self._path('solute.itp')
        solute_top = self._path('solute.top')
        solute_gro = self._path('solute.gro')
        self.solute_ff.write_itp(str(solute_itp), 'MOL')
        self.solute_ff.write_top(str(solute_top), str(solute_itp), 'MOL')
        self.solute_ff.write_gro(str(solute_gro), 'MOL')

        if not equilibration:
            self.ostream.print_info("solute.itp file written")
            self.ostream.flush()

        if self.solvent_ffs:
            for i, solvent_ff in enumerate(self.solvent_ffs):
                solvent_ff.write_itp(str(self._path(f'solvent_{i+1}.itp')),
                                     f'SOL{i+1}')
                if not equilibration:
                    self.ostream.print_info(f"solvent_{i+1}.itp file written")
                    self.ostream.flush()

    def _collect_component_atomtypes(self):
        """
        Collect unique atomtype entries from the component ITP files.
        """

        atomtypes = []
        atomtypes.extend(self._extract_atomtypes(str(self._path('solute.itp'))))

        if self.solvent_ffs:
            for i in range(len(self.solvent_ffs)):
                atomtypes.extend(
                    self._extract_atomtypes(
                        str(self._path(f'solvent_{i+1}.itp'))))

        return list(set(atomtypes))

    def _strip_component_atomtypes(self):
        """
        Remove atomtype sections from the component ITP files after extraction.
        """

        solute_atomtypes_lines = self._remove_atomtypes_section(
            str(self._path('solute.itp')))
        if self.solvent_ffs:
            for i in range(len(self.solvent_ffs)):
                self._remove_atomtypes_section(
                    str(self._path(f'solvent_{i+1}.itp')))

        return solute_atomtypes_lines

    def _inject_solute_atomtypes_into_top(self, solute_atomtypes_lines):
        """
        Insert the extracted solute atomtype block into the solute topology file.
        """

        solute_top = self._path('solute.top')

        with open(solute_top, 'r') as f:
            solute_lines = f.readlines()
        with open(solute_top, 'w') as f:
            for line in solute_lines:
                if '"solute.itp"' in line:
                    f.write(''.join(solute_atomtypes_lines) + '\n')
                f.write(line)

    def _write_system_top_file(self, atomtypes):
        """
        Assemble the combined system topology from the generated component files.
        """

        with open(self._path('system.top'), 'w') as f:
            f.write('[ defaults ]\n')
            f.write(
                '; nbfunc        comb-rule       gen-pairs        fudgeLJ   fudgeQQ\n'
            )
            f.write(
                '1               2               yes             0.500000  0.833333\n\n'
            )
            f.write('[ atomtypes ]\n')
            for atomtype in atomtypes:
                f.write(atomtype + '\n')
            f.write('\n')
            f.write(';Residue topologies\n')
            f.write('#include "solute.itp"\n')
            for i in range(len(self.solvent_ffs)):
                f.write(f'#include "solvent_{i+1}.itp"\n')
            if self.counterion:
                f.write(
                    f'#include "{self.parent_forcefield}.ff/forcefield.itp"\n'
                )
                f.write(f'#include "{self.parent_forcefield}.ff/ions.itp"\n')
            f.write('\n[ system ]\n')
            f.write('System\n\n')
            f.write('[ molecules ]\n')
            f.write('MOL 1\n')
            for i, count in enumerate(self.added_solvent_counts):
                f.write(f'SOL{i+1} {count}\n')
            if self.counterion:
                residue_name = self.ion_name.upper()
                f.write(f'{residue_name} {abs(self.added_counterions)}\n')

    def write_openmm_files(self, solute_ff=None, solvent_ffs=None):
        '''
        Generates the PDB and XML files for OpenMM.

        :param solute_ff:
            The ForceField object of the solute (optional).
        :param solvent_ffs:
            The list of ForceField objects of the solvent molecules (optional).
        '''

        # Generate the forcefields
        self._generate_forcefields(solute_ff, solvent_ffs)

        # Special case for 'itself' solvent
        if self.solvent_name == 'itself':
            # Write the XML and PDB files for the solute
            if not self.write_pdb_only:
                self.solute_ff.write_openmm_files(str(self._path('liquid')),
                                                  'MOL')
                self.ostream.print_info("liquid.xml file written")
                self.ostream.flush()
            # Write the system PDB file
            filename = 'liquid.pdb'

        else:
            # Solute
            if not self.write_pdb_only:
                self.solute_ff.write_openmm_files(str(self._path('solute')),
                                                  'MOL')
                self.ostream.print_info(
                    "solute.pdb, and solute.xml files written")
                self.ostream.flush()

            for i, solvent_ff in enumerate(self.solvent_ffs):
                solvent_ff.generate_residue_xml(
                                                str(self._path(
                                                    f'solvent_{i+1}.xml')),
                                                f'S{i+1:02d}')
                self.ostream.print_info(f"solvent_{i+1}.xml file written")
                self.ostream.flush()

            filename = 'system.pdb'

        # Write the system PDB file
        if self.equilibration_flag:
            # If the system was equilibrated, remove pdb to align with new resnames etc.
            self._unlink_if_exists(self._path('equilibrated_system.pdb'))

        self._write_system_pdb(filename=filename)
        # Print information
        self.ostream.print_info(f"{filename} file written")
        self.ostream.flush()

    def perform_equilibration(self, water_model=None, steps=None):
        """
        Performs an equilibration using OpenMM.

        :param water_model:
            Optional water model name to use when equilibrating a pure-water
            `solvent='itself'` system.
        :param steps:
            Optional number of MD steps to run. If None, `self.steps` is used.
        """

        try:
            import openmm as mm
            import openmm.app as app
            import openmm.unit as unit

        except ImportError:
            raise ImportError("OpenMM is required for this functionality")

        run_steps = self.steps if steps is None else steps

        # Generate the forcefields using semiempirical charges.

        solute_ff = MMForceFieldGenerator()
        solute_ff.ostream.mute()
        solute_ff.partial_charges = self.solute.get_partial_charges(
            self.solute.get_charge())

        if self.solvent_name == 'itself' and self.solute.is_water_molecule():
            assert_msg_critical(
                water_model is not None,
                'SolvationBuilder: water_model must be provided when equilibrating a pure-water solvent=\'itself\' system.'
            )
            solute_ff.create_topology(self.solute, water_model=water_model)
        else:
            solute_ff.create_topology(self.solute)

        use_water_model = (self.solvent_name in self.water_parameters)

        if self.solvent_name in ['itself']:
            solvent_ffs = None

        else:
            solvent_ffs = []
            for i in range(len(self.solvents)):
                solvent_ff = MMForceFieldGenerator()
                solvent_ff.ostream.mute()
                solvent_ff.partial_charges = self.solvents[
                    i].get_partial_charges(self.solvents[i].get_charge())
                if use_water_model:
                    solvent_ff.create_topology(self.solvents[i],
                                               water_model=self.solvent_name)
                else:
                    if self.solvents[i].is_water_molecule():
                        # auto-detect water molecule and use ctip3p as default
                        solvent_ff.create_topology(self.solvents[i],
                                                   water_model='ctip3p')
                    else:
                        solvent_ff.create_topology(self.solvents[i])
                solvent_ffs.append(solvent_ff)

        self.write_gromacs_files(solute_ff, solvent_ffs, equilibration=True)

        # Load the system
        if self.solvent_name == 'itself':
            gro = app.GromacsGroFile(str(self._path('liquid.gro')))
        else:
            gro = app.GromacsGroFile(str(self._path('system.gro')))

        # # Create the force field
        if self.solvent_name == 'itself':
            forcefield = app.GromacsTopFile(
                str(self._path('liquid.top')),
                periodicBoxVectors=gro.getPeriodicBoxVectors())
        else:
            forcefield = app.GromacsTopFile(
                str(self._path('system.top')),
                periodicBoxVectors=gro.getPeriodicBoxVectors())

        topology = forcefield.topology
        positions = gro.positions

        # Create the OpenMM system
        system = forcefield.createSystem(nonbondedMethod=app.PME,
                                         nonbondedCutoff=1.0 * unit.nanometers,
                                         constraints=app.HBonds)

        # Set the temperature and pressure
        integrator = mm.LangevinIntegrator(self.temperature * unit.kelvin,
                                           1.0 / unit.picosecond,
                                           2 * unit.femtosecond)
        barostat = mm.MonteCarloBarostat(self.pressure * unit.bar,
                                         self.temperature * unit.kelvin, 25)
        system.addForce(barostat)

        # Create the simulation
        simulation = app.Simulation(topology, system, integrator)
        simulation.context.setPositions(positions)

        # Minimize the energy
        simulation.minimizeEnergy()

        # Equilibrate
        simulation.reporters.append(
            app.StateDataReporter(str(self._path('equilibration.log')),
                                  1000,
                                  step=True,
                                  potentialEnergy=True,
                                  temperature=True,
                                  volume=True))
        simulation.step(run_steps)

        # Get the final positions
        positions = simulation.context.getState(
            getPositions=True).getPositions()

        # Exctract the new box size from the simulation context
        box_vectors = simulation.context.getState().getPeriodicBoxVectors()
        # Update the box size in angstrom
        self.box = [
            box_vectors[0][0].value_in_unit(unit.angstroms),
            box_vectors[1][1].value_in_unit(unit.angstroms),
            box_vectors[2][2].value_in_unit(unit.angstroms)
        ]
        self.ostream.print_info(
            f'The box size after equilibration is: {self.box[0] * 0.1:.2f} x {self.box[1] * 0.1:.2f} x {self.box[2] * 0.1:.2f} nm^3'
        )
        self.ostream.flush()
        box_volume_nm3 = self.box[0] * self.box[1] * self.box[2] * 1e-3
        # Recalculate the density of the system
        self.ostream.print_info(
            f'The density of the system after equilibration is: {self._compute_density(box_volume_nm3)} kg/m^3'
        )
        self.ostream.flush()
        # Write the PDB file
        with open(self._path('equilibrated_system.pdb'), 'w') as f:
            app.PDBFile.writeFile(simulation.topology, positions, f)

        # Delete the produced gro and top files with Path
        if self.solvent_name == 'itself':
            self._unlink_if_exists(self._path('liquid.gro'))
            self._unlink_if_exists(self._path('liquid.top'))
            self._unlink_if_exists(self._path('liquid.itp'))
        else:
            self._unlink_if_exists(self._path('system.gro'))
            self._unlink_if_exists(self._path('system.top'))
            self._unlink_if_exists(self._path('solute.itp'))
            for i in range(len(self.solvent_ffs)):
                self._unlink_if_exists(self._path(f'solvent_{i+1}.itp'))

        # Update the system molecule
        self.system_molecule = Molecule.read_pdb_file(
            str(self._path('equilibrated_system.pdb')))

    # Auxiliary functions

    def _determine_cubic_box_size(self, coordinates, padding):
        """
        Determines the size of a cubic box based on the molecular geometry and a padding parameter.
        
        :param coordinates: 
            The array of shape (n, 3) where n is the number of atoms, and each row contains the (x, y, z) coordinates of an atom.
        :param padding: 
            The padding to be added to the bounding box of the molecule in nm.
        :return: 
            The size of the cubic box in the same units as the input coordinates.
        """
        # Calculate the minimum and maximum coordinates along each axis
        min_coords = np.min(coordinates, axis=0)
        max_coords = np.max(coordinates, axis=0)

        # Determine the extent of the bounding box
        extent = max_coords - min_coords

        # Add the padding to the extent
        box_size = np.max(extent) + 2 * padding * 10

        return box_size

    def _unlink_if_exists(self, filename):
        """
        Remove a generated file if it exists.

        :param filename:
            The file to remove.
        """

        path = Path(filename)
        if path.exists():
            path.unlink()

    def _load_solute_molecule(self, solute):
        '''
        Register the solute molecule to be added to the system
        '''

        self.solute = solute
        self.solute_labels = solute.get_labels()

    def _define_box(self, dim_x, dim_y, dim_z):
        '''
        Create a box around the solute
        '''

        self.box = [dim_x, dim_y, dim_z]

    def _clear_system(self):
        """
        Clear the registered system.
        """

        self.system.clear()

    def _clear_solvent_molecules(self):
        """
        Clear the registered solvent molecules.
        """

        self.solvents.clear()
        self.solvent_labels.clear()
        self.quantities.clear()
        self.added_solvent_counts.clear()
        self.random_rotation = True

    def _clear_counterions(self):
        """
        Clear the registered counter ions.
        """

        self.counterion = None
        self.added_counterions = 0

    def _load_solvent_molecule(self, solvent, quantity):
        """
        Register the solvent molecule and its quantity to be added to the system

        :param solvent:
            The VeloxChem molecule object of the solvent
        :param quantity:
            The quantity of the solvent molecules to be added
        """

        self.solvents.append(solvent)
        self.quantities.append(quantity)
        self.solvent_labels.append(solvent.get_labels())

    def _set_last_solvent_quantity(self, quantity):
        """
        Update the quantity of the most recently registered solvent molecule.

        :param quantity:
            The quantity of the solvent molecules to be added.
        """

        if self.quantities:
            self.quantities[-1] = quantity

    def _report_partial_fill(self,
                             quantity,
                             added_count,
                             failure_count=None,
                             reason=None):
        """
        Reports a partial-fill outcome after repeated insertion failures.

        :param quantity:
            The requested number of molecules.
        :param added_count:
            The number of molecules that were added successfully.
        :param failure_count:
            The number of consecutive failed insertion attempts.
        :param reason:
            Optional explanation for why packing stopped early.
        """

        msg = (f"Solvated system with {added_count} solvent molecules out of "
               f"{quantity} requested")
        self.ostream.print_info(msg)
        self.ostream.flush()

    def _report_solvation_result(self, added_count, quantity):
        """
        Reports the final number of packed solvent molecules.

        :param added_count:
            The number of molecules that were added successfully.
        :param quantity:
            The requested number of molecules.
        """

        msg = f"Solvated system with {added_count} solvent molecules out of {quantity} requested"
        self.ostream.print_info(msg)
        self.ostream.flush()

    def _insert_molecule(self, new_molecule, tree):
        """
        Insert a molecule with a random rotation without rebuilding the KDTree each time.
        """

        from scipy.spatial.transform import Rotation

        new_molecule_xyz = new_molecule.get_coordinates_in_angstrom()
        new_molecule_labels = new_molecule.get_labels()

        total_attempts = self.number_of_attempts

        for _ in range(total_attempts):

            if self.random_rotation:
                # Generate a random rotation
                rotation_matrix = Rotation.random().as_matrix()
                # Rotate the molecule coordinates
                rotated_coords = new_molecule_xyz @ rotation_matrix.T

            else:
                rotated_coords = new_molecule_xyz

            # Generate a random position for the solvent molecule
            position = np.random.rand(3) * self.box

            # Translate the rotated solvent molecule to the new position
            translated_solvent_coords = rotated_coords + position

            # Check if the solvent molecule is within the box
            if not np.all((translated_solvent_coords >= 0) &
                          (translated_solvent_coords <= self.box)):
                continue

            # Check for overlap using KD-Tree
            if tree:
                # Query the KD-Tree for atoms within the threshold distance
                indices = tree.query_ball_point(translated_solvent_coords,
                                                self.threshold)
                if any(indices):
                    continue

            # No overlap detected; return the translated molecule
            molecule_id = len(self.system)
            translated_solvent = [
                (molecule_id, label, coord) for label, coord in zip(
                    new_molecule_labels, translated_solvent_coords)
            ]
            return translated_solvent

        return None

    def _check_overlap(self, new_coords, existing_tree):
        '''
        Check for overlap using a KD-Tree.
        '''
        indices = existing_tree.query_ball_point(new_coords, self.threshold)
        # If any indices are returned, overlap exists
        return any(indices)

    def _save_molecule(self):
        '''
        Creates a VeloxChem molecule object from the system
        '''

        # Extract atom labels and coordinates from the system
        labels = [atom[1] for atom in self.system]
        coordinates = np.array([atom[-1] for atom in self.system])

        # Create the VeloxChem molecule object
        molecule = Molecule(symbols=labels,
                            coordinates=coordinates,
                            units='angstrom')

        return molecule

    def _extract_atomtypes(self, itp_filename):
        """
        Extracts atom types from an ITP file, removes duplicates, and returns a list of unique atom types.
        
        :param itp_filename: The name of the ITP file from which to extract atom types.
        :return: A list of unique atom types.
        """
        atomtypes = []
        inside_atomtypes_block = False

        with open(itp_filename, 'r') as file:
            for line in file:
                if line.strip().startswith('[ atomtypes ]'):
                    inside_atomtypes_block = True
                    continue  # Skip the [ atomtypes ] header line

                if inside_atomtypes_block:
                    if line.strip() == '' or line.startswith('['):
                        # End of the atomtypes block
                        inside_atomtypes_block = False
                    elif not line.strip().startswith(';'):
                        # Add the atom type line to the list if it's not a comment
                        atomtypes.append(line.strip())

        # Remove duplicates by converting to a set and back to a list
        unique_atomtypes = list(set(atomtypes))

        return unique_atomtypes

    def _remove_atomtypes_section(self, itp_filename):
        """
        Removes the [ atomtypes ] section from an ITP file.
        
        :param itp_filename:
            The name of the ITP file from which to remove the atom types
            section.

        :return:
            The removed lines.
        """

        new_lines = []
        removed_lines = []
        inside_atomtypes_block = False

        with open(itp_filename, 'r') as file:
            for line in file:
                if line.strip().startswith('[ atomtypes ]'):
                    inside_atomtypes_block = True
                    removed_lines.append(line)
                    continue  # Skip the [ atomtypes ] header line

                if inside_atomtypes_block:
                    if line.strip() == '' or line.startswith('['):
                        # End of the atomtypes block, resume normal line processing
                        inside_atomtypes_block = False
                        new_lines.append(line)
                    else:
                        removed_lines.append(line)
                else:
                    new_lines.append(line)

        # Rewrite the itp file without the atomtypes block
        with open(itp_filename, 'w') as file:
            file.writelines(new_lines)

        return removed_lines

    def _solvent_properties(self, solvent):
        '''
        Contains the density of the solvent

        :param solvent:
            The name of the solvent
        '''

        if solvent.lower() in self.water_parameters:
            mols_per_nm3 = 33.3
            density = 1000
            smiles_code = 'O'

        elif solvent == 'ethanol':
            mols_per_nm3 = 10.3
            density = 789
            smiles_code = 'CCO'

        elif solvent == 'methanol':
            mols_per_nm3 = 14.9
            density = 791
            smiles_code = 'CO'

        elif solvent == 'acetone':
            mols_per_nm3 = 8.1
            density = 784
            smiles_code = 'CC(=O)C'

        elif solvent == 'chloroform':
            mols_per_nm3 = 7.6
            density = 1494
            smiles_code = 'C(Cl)(Cl)Cl'

        elif solvent == 'hexane':
            mols_per_nm3 = 4.7
            density = 655
            smiles_code = 'CCCCCC'

        elif solvent == 'toluene':
            mols_per_nm3 = 5.7
            density = 866
            smiles_code = 'CC1=CC=CC=C1'

        elif solvent == 'dcm':
            mols_per_nm3 = 8.3
            density = 1335
            smiles_code = 'ClCCl'

        elif solvent == 'benzene':
            mols_per_nm3 = 6.8
            density = 876
            smiles_code = 'c1ccccc1'

        elif solvent == 'dmso':
            mols_per_nm3 = 8.5
            density = 1101
            smiles_code = 'CS(=O)C'

        elif solvent == 'thf':
            mols_per_nm3 = 7.4
            density = 888
            smiles_code = 'C1CCOC1'

        elif solvent == 'acetonitrile':
            mols_per_nm3 = 11.53
            density = 786
            smiles_code = 'CC#N'

        elif solvent == 'dmf':
            mols_per_nm3 = 7.78
            density = 944
            smiles_code = 'CN(C)C=O'

        else:
            raise ValueError(
                f"Unsupported solvent '{solvent}'. "
                "Use 'other' with solvent_molecule for a custom solvent.")

        return mols_per_nm3, density, smiles_code

    def _counterion_molecules(self):
        """
        Generates a molecule object for the counterion.
        """
        coords = np.array([[0.0, 0.0, 0.0]])

        ion_definitions = {
            'Na': (['Na'], 1.0),
            'K': (['K'], 1.0),
            'Li': (['Li'], 1.0),
            'Cl': (['Cl'], -1.0),
        }
        assert_msg_critical(
            self.ion_name in ion_definitions,
            f"Unsupported counterion '{self.ion_name}'. Supported counterions are: "
            f"{', '.join(ion_definitions)}.")
        labels, charge = ion_definitions[self.ion_name]

        ion_mol = Molecule(labels, coords, 'angstrom')
        ion_mol.set_charge(charge)

        return ion_mol

    def _compute_density(self, volume):
        """
        Return the density of the whole periodic system in kg/m^3.

        This includes the solute, all solvent species, and any inserted
        counterions, divided by the full simulation box volume.

        :param volume:
            The volume in cubic nm.
        """

        volume_m3 = volume * 1e-27
        total_mass = 0.0

        solute_count = 1
        if self.solvent_name == 'itself' and self.added_solvent_counts:
            solute_count += self.added_solvent_counts[0]
            solvent_start = 1
        else:
            solvent_start = 0

        if self.solute is not None:
            solute_mass = sum(self.solute.get_masses()) * 1e-3
            total_mass += solute_mass * solute_count / 6.022e23

        for molecule, number_of_molecules in zip(self.solvents[solvent_start:],
                                                 self.added_solvent_counts[solvent_start:]):
            mass = sum(molecule.get_masses()) * 1e-3
            moles = number_of_molecules / 6.022e23
            total_mass += mass * moles

        if self.counterion is not None and self.added_counterions > 0:
            ion_mass = sum(self.counterion.get_masses()) * 1e-3
            total_mass += ion_mass * self.added_counterions / 6.022e23

        return int(total_mass / volume_m3)

    def _density_to_mols_per_nm3(self, molecule, density):
        """
        Given the density in kg/m^3, return the number of moles per nm^3.

        :param density:
            The density in kg/m^3.

        :return:
            The number of moles per nm^3.
        """
        assert_msg_critical(
            density > 0.0,
            'SolvationBuilder: The target density must be positive.')

        # Get the mass in kg/mol
        mass = sum(molecule.get_masses()) * 1e-3 / 6.022e23

        # Get the mols per m^3
        mols_per_m3 = density / mass

        # Convert to mols per nm^3
        mols_per_nm3 = mols_per_m3 * 1e-27

        return mols_per_nm3

    def _generate_forcefields(self,
                              solute_ff,
                              solvent_ffs,
                              equilibration=False):
        """
        Generate the force fields for the solute and solvent molecules.
        The forcefields get stored in the solute_ff and solvent_ffs attributes (lists).

        :param solute_ff:
            The ForceField object of the solute.
        :param solvent_ffs:
            The list of ForceField objects of the solvent molecules.
        :param equilibration:
            Boolean flag to indicate if the gromacs files will be used for equilibration.
            If True, printouts will not be displayed.
        """

        # Solute
        if not solute_ff:
            self.solute_ff = MMForceFieldGenerator()
            self.solute_ff.ostream.mute()
            if not equilibration:
                self.ostream.print_info(
                    'Generating the ForceField for the solute')
                self.ostream.flush()
            self.solute_ff.create_topology(self.solute)

        else:
            self.solute_ff = solute_ff

        # Solvents
        # Special case for itself
        if not solvent_ffs:
            use_water_model = (self.solvent_name in self.water_parameters)

            if self.solvent_name == 'itself':
                self.solvent_ffs = []
                self.solvent_ffs.append(self.solute_ff)

            else:
                self.solvent_ffs = []
                for solvent in self.solvents:
                    solvent_ff = MMForceFieldGenerator()
                    solvent_ff.ostream.mute()

                    if not equilibration:
                        self.ostream.print_info(
                            f'Generating the ForceField for the solvent')
                        self.ostream.flush()

                    if use_water_model:
                        solvent_ff.create_topology(
                            solvent, water_model=self.solvent_name)
                    else:
                        if solvent.is_water_molecule():
                            # auto-detect water molecule and use ctip3p as default
                            solvent_ff.create_topology(solvent,
                                                       water_model='ctip3p')
                        else:
                            solvent_ff.create_topology(solvent)

                    self.solvent_ffs.append(solvent_ff)

        else:
            self.solvent_ffs = solvent_ffs

    def _write_system_gro(self, filename='system.gro'):
        """
        Write the solvated system to a GRO file.
        """

        filename = self._path(filename)
        coords_in_nm = self.system_molecule.get_coordinates_in_angstrom() * 0.1

        with open(filename, 'w') as f:
            self._write_gro_header(f)

            if self.solvent_name == 'itself':
                self._write_self_solvent_gro_atoms(f, coords_in_nm)
            else:
                res_id_counter, atom_id_counter = self._write_solute_gro_atoms(
                    f, coords_in_nm)
                res_id_counter, atom_id_counter = self._write_solvent_gro_atoms(
                    f, coords_in_nm, res_id_counter, atom_id_counter)
                self._write_counterion_gro_atoms(f, coords_in_nm, res_id_counter,
                                                 atom_id_counter)

            self._write_gro_box(f)

    def _write_system_pdb(self, filename='system.pdb'):
        """
        Write the solvated system to a PDB file.
        """

        filename = self._path(filename)
        # Write the system PDB file
        with open(filename, 'w') as f:
            f.write("HEADER    Generated by VeloxChem\n")
            f.write("TITLE     Solvated system\n")
            # Write the box size to the PDB file
            # Box size is in nm, convert to Angstrom
            f.write(
                f"CRYST1{self.box[0]:9.3f}{self.box[1]:9.3f}{self.box[2]:9.3f}  90.00  90.00  90.00 P 1           1\n"
            )

            atom_counter = 1
            residue_counter = 1
            chain_ids = ['A', 'B', 'C']
            pdb_atom_numbers = {}
            coordinates = self.system_molecule.get_coordinates_in_angstrom()

            assert_msg_critical(
                self.system_molecule.number_of_atoms() <= 99999,
                "The total number of atoms exceeds 99999. " +
                "The PDB format does not support more than 99999 atoms.")

            atom_counter, residue_counter, coordinate_counter = self._write_solute_pdb_atoms(
                f, coordinates, chain_ids, pdb_atom_numbers, atom_counter,
                residue_counter)
            atom_counter, residue_counter, coordinate_counter = self._write_solvent_pdb_atoms(
                f, coordinates, chain_ids, pdb_atom_numbers, atom_counter,
                residue_counter, coordinate_counter)
            self._write_counterion_pdb_atoms(f, coordinates, chain_ids,
                                             pdb_atom_numbers, atom_counter,
                                             residue_counter, coordinate_counter)
            self._write_pdb_connect_records(f, pdb_atom_numbers)

            # Write END line
            f.write("END\n")

    def _write_gro_header(self, stream):
        """
        Write the standard GRO file header.
        """

        stream.write('Generated by VeloxChem\n')
        stream.write(f'{self.system_molecule.number_of_atoms()}\n')

    def _write_self_solvent_gro_atoms(self, stream, coords_in_nm):
        """
        Write GRO atom records for the pure-liquid special case.
        """

        atom_id_counter = 0
        for mol_id in range(self.added_solvent_counts[0] + 1):
            for atom in self.solute_ff.atoms.values():
                atom_name = atom['name']
                line_str = f'{(mol_id + 1) % 100000:>5d}'
                line_str += f'{"MOL":<5s}{atom_name:>5s}'
                line_str += f'{(atom_id_counter + 1) % 100000:>5d}'
                for d in range(3):
                    line_str += f'{coords_in_nm[atom_id_counter][d]:8.3f}'
                line_str += '\n'
                stream.write(line_str)
                atom_id_counter += 1

    def _write_solute_gro_atoms(self, stream, coords_in_nm):
        """
        Write the solute atom records in GRO format.
        """

        atom_id_counter = 0
        for i, atom in self.solute_ff.atoms.items():
            atom_name = atom['name']
            line_str = f'{1:>5d}{"MOL":<5s}{atom_name:>5s}{i + 1:>5d}'
            for d in range(3):
                line_str += f'{coords_in_nm[i][d]:8.3f}'
            line_str += '\n'
            stream.write(line_str)
            atom_id_counter += 1

        return 1, atom_id_counter

    def _write_solvent_gro_atoms(self, stream, coords_in_nm, res_id_counter,
                                 atom_id_counter):
        """
        Write the solvent atom records in GRO format.
        """

        for i, solvent in enumerate(self.solvent_ffs):
            for _ in range(self.added_solvent_counts[i]):
                resid_k = (res_id_counter + 1) % 100000
                for atom in solvent.atoms.values():
                    atom_name = atom['name']
                    atomid_j = (atom_id_counter + 1) % 100000
                    line_str = f'{resid_k:>5d}{f"SOL{i+1}":<5s}{atom_name:>5s}{atomid_j:>5d}'
                    for d in range(3):
                        line_str += f'{coords_in_nm[atom_id_counter][d]:8.3f}'
                    line_str += '\n'
                    stream.write(line_str)
                    atom_id_counter += 1
                res_id_counter += 1

        return res_id_counter, atom_id_counter

    def _write_counterion_gro_atoms(self, stream, coords_in_nm, res_id_counter,
                                    atom_id_counter):
        """
        Write the counterion atom records in GRO format.
        """

        if self.counterion:
            atom_name = self.ion_name.upper()
            for _ in range(self.added_counterions):
                resid_i = (res_id_counter + 1) % 100000
                atomid_i = (atom_id_counter + 1) % 100000
                line_str = f'{resid_i:>5d}{atom_name:<5s}{atom_name:>5s}{atomid_i:>5d}'
                for d in range(3):
                    line_str += f'{coords_in_nm[atom_id_counter][d]:8.3f}'
                line_str += '\n'
                stream.write(line_str)
                atom_id_counter += 1
                res_id_counter += 1

    def _write_gro_box(self, stream):
        """
        Write the GRO box line.
        """

        for d in range(3):
            stream.write(f'{0.1 * self.box[d]:10.5f}')
        stream.write('\n')

    def _write_solute_pdb_atoms(self, stream, coordinates, chain_ids,
                                pdb_atom_numbers, atom_counter,
                                residue_counter):
        """
        Write the solute atom records in PDB format.
        """

        residue_name = 'MOL'
        for i, atom in self.solute_ff.atoms.items():
            atom_name = atom['name']
            element = self.solute.get_labels()[i]
            x, y, z = coordinates[i]
            stream.write(
                "{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}\n"
                .format('HETATM', atom_counter, atom_name, '', residue_name,
                        chain_ids[0], residue_counter, '', x, y, z, 1.00,
                        0.00, element))
            pdb_atom_numbers[('solute', i)] = atom_counter
            atom_counter += 1

        residue_counter += 1
        coordinate_counter = len(self.solute_ff.atoms)
        return atom_counter, residue_counter, coordinate_counter

    def _write_solvent_pdb_atoms(self, stream, coordinates, chain_ids,
                                 pdb_atom_numbers, atom_counter,
                                 residue_counter, coordinate_counter):
        """
        Write the solvent atom records in PDB format for the current solvent mode.
        """

        if self.solvent_name == 'itself':
            residue_name = 'MOL'
            for mols in range(self.added_solvent_counts[0]):
                for i, atom in self.solute_ff.atoms.items():
                    atom_name = atom['name']
                    element = self.solute.get_labels()[i]
                    x, y, z = coordinates[coordinate_counter]
                    stream.write(
                        "{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}\n"
                        .format('HETATM', atom_counter, atom_name, '',
                                residue_name, chain_ids[1],
                                residue_counter % 10000, '', x, y, z, 1.00,
                                0.00, element))
                    pdb_atom_numbers[('solvent', mols, i)] = atom_counter
                    atom_counter += 1
                    coordinate_counter += 1
                residue_counter += 1
        else:
            for i, solvent_ff in enumerate(self.solvent_ffs):
                elements = self.solvents[i].get_labels()
                for j in range(self.added_solvent_counts[i]):
                    for k in range(len(elements)):
                        residue_name = f'S{i+1:02d}'
                        atom_name = solvent_ff.atoms[k]['name']
                        element = elements[k]
                        x, y, z = coordinates[coordinate_counter]
                        stream.write(
                            "{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}\n"
                            .format('HETATM', atom_counter, atom_name, '',
                                    residue_name, chain_ids[1],
                                    residue_counter % 10000, '', x, y, z, 1.00,
                                    0.00, element))
                        pdb_atom_numbers[('solvent', i, j, k)] = atom_counter
                        atom_counter += 1
                        coordinate_counter += 1
                    residue_counter += 1

        return atom_counter, residue_counter, coordinate_counter

    def _write_counterion_pdb_atoms(self, stream, coordinates, chain_ids,
                                    pdb_atom_numbers, atom_counter,
                                    residue_counter, coordinate_counter):
        """
        Write the counterion atom records in PDB format.
        """

        if self.counterion:
            for i in range(self.added_counterions):
                atom_name = self.ion_name.upper()
                residue_name = self.ion_name.upper()
                element = self.counterion.get_labels()[0]
                x, y, z = coordinates[coordinate_counter]
                stream.write(
                    "{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}\n"
                    .format('ATOM', atom_counter, atom_name, '', residue_name,
                            chain_ids[2], residue_counter % 10000, '', x, y, z,
                            1.00, 0.00, element))
                pdb_atom_numbers[('counterion', i)] = atom_counter
                atom_counter += 1
                coordinate_counter += 1

    def _write_pdb_connect_records(self, stream, pdb_atom_numbers):
        """
        Write the PDB bond connectivity records.
        """

        for (i_atom, j_atom) in self.solute_ff.bonds:
            pdb_i = pdb_atom_numbers[('solute', i_atom)]
            pdb_j = pdb_atom_numbers[('solute', j_atom)]
            stream.write(f"CONECT{pdb_i:>5d}{pdb_j:>5d}\n")

        if self.solvent_ffs:
            if self.solvent_name == 'itself':
                for mols in range(self.added_solvent_counts[0]):
                    for (i_atom, j_atom) in self.solute_ff.bonds:
                        pdb_i = pdb_atom_numbers[('solvent', mols, i_atom)]
                        pdb_j = pdb_atom_numbers[('solvent', mols, j_atom)]
                        stream.write(f"CONECT{pdb_i:>5d}{pdb_j:>5d}\n")
            else:
                for i, solvent_ff in enumerate(self.solvent_ffs):
                    for j in range(self.added_solvent_counts[i]):
                        for (i_atom, j_atom) in solvent_ff.bonds:
                            pdb_i = pdb_atom_numbers[('solvent', i, j, i_atom)]
                            pdb_j = pdb_atom_numbers[('solvent', i, j, j_atom)]
                            stream.write(f"CONECT{pdb_i:>5d}{pdb_j:>5d}\n")

    def _get_volume(self, molecule):
        """
        Determine the volume of the molecule based on the vdW radii, with careful intersection correction.
        
        :return:
            The volume of the molecule in cubic Angstrom.
        """
        # Get the atomic van der Waals radii in Angstrom
        atomic_radii = molecule.vdw_radii_to_numpy() * bohr_in_angstrom()

        # Get the coordinates
        coords = molecule.get_coordinates_in_angstrom()

        # Get the number of atoms
        natoms = molecule.number_of_atoms()

        # Initialize the volume to be the sum of individual atomic volumes
        volume = 0.0

        # The initial volume is the sum of the volumes of the vdW spheres
        for i in range(natoms):
            volume += (4.0 / 3.0) * np.pi * atomic_radii[i]**3

        # Correct the volume by subtracting the intersection volumes
        for i in range(natoms):
            for j in range(i + 1, natoms):
                rij = np.linalg.norm(coords[i] - coords[j])
                if rij < (atomic_radii[i] + atomic_radii[j]):
                    # Calculate the volume of the intersection of the spheres
                    r1 = atomic_radii[i]
                    r2 = atomic_radii[j]
                    d = rij

                    # Calculate the intersection volume using a more precise formula
                    # Formula from: https://mathworld.wolfram.com/Sphere-SphereIntersection.html
                    if d > 0:
                        h1 = (r1 - r2 + d) * (r1 + r2 - d) / (2 * d)
                        h2 = (r2 - r1 + d) * (r2 + r1 - d) / (2 * d)

                        intersection_volume = (np.pi * h1**2 *
                                               (3 * r1 - h1) + np.pi * h2**2 *
                                               (3 * r2 - h2)) / 6
                        volume -= intersection_volume

        # Round the volume to 2 decimal places
        volume = round(volume, 2)

        return volume

    def _molecule_linearity(self, molecule):
        """
        Determines if a molecule is linear based on its moment of inertia tensor.

        Parameters:
        - molecule: An object with methods `get_coordinates_in_angstrom()` and
        `number_of_atoms()`. Each atom should have a mass accessible.

        Returns:
        - bool: True if the molecule is linear, False otherwise.
        """

        # Get the coordinates and masses
        coords = molecule.get_coordinates_in_angstrom()
        masses = molecule.get_masses()

        # Calculate the center of mass (COM)
        com = np.average(coords, axis=0, weights=masses)

        # Calculate the moment of inertia tensor as:
        # I = sum(m_i * (ri^2 * 1 - r_i x r_i))
        inertia_tensor = np.zeros((3, 3))
        for i in range(molecule.number_of_atoms()):
            r = coords[i] - com
            inertia_tensor += masses[i] * np.outer(r, r)

        # Diagonalize the inertia tensor
        eigvals, _ = np.linalg.eigh(inertia_tensor)

        # Sort eigenvalues in ascending order
        eigvals = np.sort(eigvals)

        # In a linear molecule one of the eigenvalues is much larger than the other two
        l_1 = eigvals[0]
        l_2 = eigvals[1]
        l_3 = eigvals[2]

        # Linearity coefficient is the ratio of the sum of the two smallest eigenvalues
        # to the largest eigenvalue
        linearity_coefficient = 0.5 * (l_1 + l_2) / l_3

        if linearity_coefficient < 0.05:
            is_linear = True
        else:
            is_linear = False

        return is_linear

    def _compute_batch_size(self, added_count, max_batch_size, min_batch_size,
                            total_quantity):
        """
        Compute the batch size dynamically based on the number of molecules added.
        
        :param added_count: Number of molecules already added to the system.
        :param max_batch_size: Maximum batch size to start with.
        :param min_batch_size: Minimum batch size to avoid being too small.
        :param total_quantity: Total number of molecules to be added.
        :return: Adjusted batch size.
        """
        fraction_completed = added_count / total_quantity
        batch_size = int(max_batch_size * (1 - fraction_completed))
        batch_size = max(batch_size, min_batch_size)
        return batch_size

    def show_solvation_box(self, width=600, height=500, solvent_opacity=0.8):

        n_solute_atoms = self.solute.number_of_atoms()

        solute_atoms = []
        for atom in self.system[:n_solute_atoms]:
            solute_atoms.append([atom[1]] + [str(x) for x in atom[2]])
        solute_atoms_xyz = f'{len(solute_atoms)}\n\n'
        for a in solute_atoms:
            solute_atoms_xyz += ' '.join(a) + '\n'

        solvent_atoms = []
        for atom in self.system[n_solute_atoms:]:
            solvent_atoms.append([atom[1]] + [str(x) for x in atom[2]])
        solvent_atoms_xyz = f'{len(solvent_atoms)}\n\n'
        for a in solvent_atoms:
            solvent_atoms_xyz += ' '.join(a) + '\n'

        try:
            import py3Dmol
            viewer = py3Dmol.view(width=width, height=height)
            viewer.setViewStyle({"style": "outline", "width": 0.02})

            viewer.addModel(solute_atoms_xyz)
            viewer.setStyle({'model': 0}, {
                "stick": {},
                "sphere": {
                    "scale": 0.25
                }
            })

            viewer.addModel(solvent_atoms_xyz)
            viewer.setStyle(
                {
                    'model': 1,
                },
                {
                    "stick": {
                        "radius": 0.1,
                        "opacity": solvent_opacity
                    },
                    "sphere": {
                        "scale": 0.08,
                        "opacity": solvent_opacity
                    },
                },
            )

            viewer.zoomTo()
            viewer.show()

        except ImportError:
            raise ImportError('Unable to import py3Dmol')
