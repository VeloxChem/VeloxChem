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
import numpy as np
import sys
from pathlib import Path
from .veloxchemlib import mpi_master
from .sanitychecks import molecule_sanity_check
from .molecule import Molecule
from .molecularbasis import MolecularBasis
from .scfunrestdriver import ScfUnrestrictedDriver
from .scfrestdriver import ScfRestrictedDriver
from .scfhessiandriver import ScfHessianDriver
from .respchargesdriver import RespChargesDriver
from .xtbdriver import XtbDriver
from .optimizationdriver import OptimizationDriver
from .mmforcefieldgenerator import MMForceFieldGenerator
from .reactionmatcher import ReactionMatcher
from .outputstream import OutputStream
from .veloxchemlib import Point
from .errorhandler import assert_msg_critical, print_exception_if_debug
from .openmmdynamics import OpenMMDynamics
from .reactionsystembuilder import ReactionSystemBuilder

try:
    import openmm as mm
    import openmm.app as mmapp
    import openmm.unit as mmunit
except ImportError:
    pass


class ReactionForceFieldBuilder:

    def __init__(self, comm=None, ostream=None):
        if comm is None:
            comm = MPI.COMM_WORLD

        if ostream is None:
            if comm.Get_rank() == mpi_master():
                ostream = OutputStream(sys.stdout)
            else:
                ostream = OutputStream(None)

        # output stream
        self.ostream = ostream

        # mpi information
        self.comm = comm
        self.rank = self.comm.Get_rank()
        self.nodes = self.comm.Get_size()
        self.calculate_resp: bool = True
        self.optimize_mol: bool = False
        self.reparameterize_bonds: bool = False
        self.reparameterize_angles: bool = False
        self.optimize_ff: bool = True
        self.optimize_steps: int = 1000
        self.optimize_conformer_snapshots: int = 10
        self.optimize_temp: int = 600
        self.optimize_dist_restraint_offset = 0.5  # Angstrom
        self.optimize_dist_restraint_k = 1000.0  # kJ mol^-1 nm^-2
        self.mm_opt_constrain_bonds: bool = True
        self.water_model: str = 'cspce'
        self.product_mapping: dict[int, int] | None = None  # one-indexed
        self.mute_scf: bool = True
        self.skip_reaction_matching: bool = False

        self.hessian_xc_fun: str = 'B3LYP'
        self.hessian_basis = 'def2-SV_P_'

        # Relative spread, (max - min) / mean, above which the averaged
        # reference force constant of a reactive pair is reported as
        # unrepresentative. See _build_reactive_bond_registry.
        self.fc_ref_spread_warning: float = 0.5

    def _sub_ostream(self, silent):
        # Output stream to hand to a sub-driver. Never leave it to the driver
        # default: that opens a stream on whatever sys.stdout happens to be at
        # construction time, which then outlives the call.
        return OutputStream(None) if silent else self.ostream

    def build_force_fields(
        self,
        reactant: Molecule | list[Molecule],
        product: Molecule | list[Molecule],
        reactant_partial_charges: list[float] | list[list[float]] | None = None,
        product_partial_charges: list[float] | list[list[float]] | None = None,
        reactant_hessians: np.ndarray | list[np.ndarray | None] | None = None,
        product_hessians: np.ndarray | list[np.ndarray | None] | None = None,
        reactant_total_multiplicity: int = -1,
        product_total_multiplicity: int = -1,
        forced_breaking_bonds: set[tuple[int, int]] | tuple = (),
        forced_forming_bonds: set[tuple[int, int]] | tuple = (),
        product_mapping: dict[int, int] | None = None,
    ):
        """Create matched forcefields for the reactant and product

        Args:
            reactant (Molecule | list[Molecule]): The reactant molecule or a list of reactant molecules.
            product (Molecule | list[Molecule]): The product molecule or a list of product molecules.
            reactant_partial_charges (list[float], list[list[float]]): Partial charges for the reactant. Will be calculated if not provided. Defaults to None.
            product_partial_charges (list[float], list[list[float]]): Partial charges for the product. Will be calculated if not provided. Defaults to None.
            reactant_hessians (np.ndarray, list[np.ndarray]): Hessians for the reactant for the Seminario method. Will be calculated if not provided. Defaults to None.
            product_hessians (np.ndarray, list[np.ndarray]): Hessians for the product for the Seminario method. Will be calculated if not provided. Defaults to None.
            reactant_total_multiplicity (int): Total multiplicity for the reactant to override calculated value. Defaults to -1.
            product_total_multiplicity (int): Total multiplicity for the product to override calculated value. Defaults to -1.
            forced_breaking_bonds (set[tuple[int, int]] | tuple, optional): A (set of) 1-indexed tuple(s) specifying bonds that are forced to break. Defaults to ().
            forced_forming_bonds (set[tuple[int, int]] | tuple, optional): A (set of) 1-indexed tuple(s) specifying bonds that are forced to form. Defaults to ().
            product_mapping (dict[int, int] | None, optional): A dictionary specifying the 1-indexed mapping from reactant to product atoms. If provided, reaction matching is skipped. Defaults to None.

        Returns:
            tuple: reactant forcefield, product forcefield, forming bonds, breaking bonds, list of reactant forcefields, list of product forcefields, product mapping (0-indexed)
        """

        # Two-state wrapper over the general N-state builder, so that there is
        # exactly one implementation of matching / index-frame composition /
        # reaction summarising. Only the argument marshalling lives here: the
        # per-state arguments become two-element lists and the per-step ones
        # single-element lists.
        results = self.build_many_force_fields(
            states=[reactant, product],
            partial_charges=[reactant_partial_charges, product_partial_charges],
            hessians=[reactant_hessians, product_hessians],
            total_multiplicities=[
                reactant_total_multiplicity, product_total_multiplicity
            ],
            forced_breaking_bonds=[forced_breaking_bonds],
            forced_forming_bonds=[forced_forming_bonds],
            mappings=[product_mapping],
            _state_labels=["REA", "PRO"],
        )

        states = results['states']
        return (
            states[0],
            states[1],
            results['forming_bonds'][0],
            results['breaking_bonds'][0],
            results['fragments'][0],
            results['fragments'][1],
            results['mappings'][0],
        )

    def build_many_force_fields(
        self,
        states: list,
        partial_charges: list | None = None,
        hessians: list | None = None,
        total_multiplicities: list[int] | None = None,
        forced_breaking_bonds: list | None = None,
        forced_forming_bonds: list | None = None,
        mappings: list | None = None,
        _state_labels: list[str] | None = None,
    ) -> dict:
        """Create matched force fields for an arbitrary number of states along a
        single reaction path (state 0 -> state 1 -> ... -> state N-1).

        Every state is brought into the atom ordering of state 0, so all states,
        every step's systems and the shared topology are indexed identically.

        Args:
            states: N entries; each is a Molecule or a list of Molecule fragments.
            partial_charges: N entries; each is None, a flat list of charges, or
                one list of charges per fragment of that state.
            hessians: N entries; each is None, an ndarray, or one per fragment.
            total_multiplicities: N entries; -1 to use the calculated value.
            forced_breaking_bonds: N-1 entries, one per step; each a set of
                1-indexed atom pairs forced to break in that step.
            forced_forming_bonds: N-1 entries, one per step.
            mappings: N-1 entries; a 1-indexed mapping from state k to state k+1
                to skip reaction matching for that step, or None.

        Returns:
            dict with keys 'states' (list of combined force fields), 'fragments'
            (per state, the list of per-fragment force fields), 'forming_bonds'
            and 'breaking_bonds' (per step), 'reactive_bonds' (the global union),
            'registry' (see _build_reactive_bond_registry) and 'mappings'.
        """

        assert_msg_critical(
            len(states) >= 2,
            f"At least two states are required, got {len(states)}")

        n_states = len(states)
        n_steps = n_states - 1

        partial_charges = self._per_state_argument(partial_charges, n_states,
                                                   "partial_charges")
        hessians = self._per_state_argument(hessians, n_states, "hessians")
        total_multiplicities = self._per_state_argument(total_multiplicities,
                                                        n_states,
                                                        "total_multiplicities",
                                                        default=-1)
        forced_breaking_bonds = self._per_step_argument(forced_breaking_bonds,
                                                        n_steps,
                                                        "forced_breaking_bonds",
                                                        default=())
        forced_forming_bonds = self._per_step_argument(forced_forming_bonds,
                                                       n_steps,
                                                       "forced_forming_bonds",
                                                       default=())
        mappings = self._per_step_argument(mappings, n_steps, "mappings",
                                           default=None)

        if _state_labels is None:
            _state_labels = [f"STATE{k + 1}" for k in range(n_states)]

        state_ffs = []
        fragment_ffs = []
        total_charges = []
        for k, state in enumerate(states):
            ff, fragments, total_charge = self._create_combined_forcefield(
                state,
                partial_charges[k],
                hessians[k],
                _state_labels[k],
                total_multiplicities[k],
            )
            state_ffs.append(ff)
            fragment_ffs.append(fragments)
            total_charges.append(total_charge)

        for k, total_charge in enumerate(total_charges[1:], start=1):
            assert_msg_critical(
                total_charge == total_charges[0],
                f"Total charge of state {k + 1} ({total_charge}) must match "
                f"that of state 1 ({total_charges[0]})")

        # Bring every state into state 0's atom ordering. The left operand of
        # each match is already expressed in state 0's frame, so the mappings
        # compose implicitly and no separate composition step is needed.
        resolved_mappings = []
        for step in range(n_steps):
            mapping = self._match_step(
                state_ffs[step],
                state_ffs[step + 1],
                forced_breaking_bonds[step],
                forced_forming_bonds[step],
                mappings[step],
                step,
                n_steps,
            )
            resolved_mappings.append(mapping)
            if mapping is not None:
                state_ffs[step + 1] = self._apply_mapping_to_forcefield(
                    state_ffs[step + 1], mapping)
                state_ffs[step + 1].molecule = self._apply_mapping_to_molecule(
                    state_ffs[step + 1].molecule, mapping)

        ReactionSystemBuilder.assert_consistent_atom_ordering(
            state_ffs, "force field building")

        forming_bonds = []
        breaking_bonds = []
        for step in range(n_steps):
            step_label = None if n_steps == 1 else f"step {step + 1}"
            formed, broken = self._summarise_reaction(state_ffs[step],
                                                      state_ffs[step + 1],
                                                      step_label=step_label)
            forming_bonds.append(formed)
            breaking_bonds.append(broken)

            for bond in broken:
                state_ffs[step].bonds[bond]['comment'] += ', broken in reaction'
            for bond in formed:
                state_ffs[step +
                          1].bonds[bond]['comment'] += ', formed in reaction'

        self.ostream.flush()

        reactive_bonds = set()
        for step in range(n_steps):
            reactive_bonds |= set(forming_bonds[step])
            reactive_bonds |= set(breaking_bonds[step])

        registry = self._build_reactive_bond_registry(state_ffs, reactive_bonds)

        self._optimize_states(state_ffs, forming_bonds, breaking_bonds)

        return {
            'states': state_ffs,
            'fragments': fragment_ffs,
            'forming_bonds': forming_bonds,
            'breaking_bonds': breaking_bonds,
            'reactive_bonds': reactive_bonds,
            'registry': registry,
            'mappings': resolved_mappings,
        }

    @staticmethod
    def _per_state_argument(value, n_states, name, default=None):
        """Normalise a per-state argument to a list of length n_states."""
        if value is None:
            return [default] * n_states
        assert_msg_critical(
            isinstance(value, list),
            f"{name} must be a list with one entry per state")
        assert_msg_critical(
            len(value) == n_states,
            f"{name} must have one entry per state: expected {n_states} "
            f"entries, got {len(value)}")
        return list(value)

    @staticmethod
    def _per_step_argument(value, n_steps, name, default=None):
        """Normalise a per-step argument to a list of length n_steps."""
        if value is None:
            return [default] * n_steps
        assert_msg_critical(
            isinstance(value, list),
            f"{name} must be a list with one entry per step")
        assert_msg_critical(
            len(value) == n_steps,
            f"{name} must have one entry per step (one fewer than the number "
            f"of states): expected {n_steps} entries, got {len(value)}")
        return list(value)

    def _match_step(self, state_ff, next_state_ff, forced_breaking_bonds,
                    forced_forming_bonds, mapping, step, n_steps):
        """Resolve the 0-indexed mapping that takes the next state into the
        current state's (i.e. state 0's) atom ordering."""

        step_insert = "" if n_steps == 1 else f" for step {step + 1}"

        if not self.skip_reaction_matching and mapping is None:
            breaking_bonds_insert = "no forced breaking bonds"
            if len(forced_breaking_bonds) > 0:
                breaking_bonds_insert = f"forced breaking bonds: {forced_breaking_bonds}"

            forming_bonds_insert = "no forced forming bonds"
            if len(forced_forming_bonds) > 0:
                forming_bonds_insert = f"forced forming bonds: {forced_forming_bonds}"

            msg = f"Matching reactant and product force fields{step_insert} with "
            msg += breaking_bonds_insert + " and " + forming_bonds_insert + "."
            self.ostream.print_info(msg)
            self.ostream.flush()

            # adjust for 1-indexed input of breaking bonds
            forced_breaking_bonds = {(bond[0] - 1, bond[1] - 1)
                                     for bond in forced_breaking_bonds}
            forced_forming_bonds = {(bond[0] - 1, bond[1] - 1)
                                    for bond in forced_forming_bonds}

            mapping = self._match_reactant_and_product(
                state_ff.molecule,
                next_state_ff.molecule,
                forced_breaking_bonds,
                forced_forming_bonds,
            )
        elif mapping is not None:
            self.ostream.print_info(
                f"Skipping reaction matching{step_insert} because the mapping {mapping} is already provided"
            )
            mapping = {k - 1: v - 1 for k, v in mapping.items()}
        else:
            self.ostream.print_info(f"Skipping reaction matching{step_insert}")
        self.ostream.flush()
        return mapping

    def _build_reactive_bond_registry(self, state_ffs, reactive_bonds):
        """Collect, for every pair that reacts anywhere along the path, which
        states it is bonded in and one path-global reference force constant.

        The reference force constant is the average over the states in which the
        pair is bonded. It is only ever used for the surrogate bonded term that
        stands in for a dissociated pair during integration - never for a
        reported energy - but it has to be identical in every step, otherwise
        the integration Hamiltonian jumps where two steps meet.
        """

        registry = {}
        for pair in sorted(reactive_bonds):
            states_bonded = [
                k for k, ff in enumerate(state_ffs) if pair in ff.bonds
            ]
            fc_per_state = {
                k: state_ffs[k].bonds[pair]['force_constant']
                for k in states_bonded
            }
            assert_msg_critical(
                len(states_bonded) > 0,
                f"Reactive pair {pair} is not bonded in any state")

            force_constants = list(fc_per_state.values())
            fc_ref = float(np.mean(force_constants))
            registry[pair] = {
                'states_bonded': states_bonded,
                'fc_per_state': fc_per_state,
                'fc_ref': fc_ref,
            }

            if len(force_constants) > 1 and fc_ref > 0:
                spread = (max(force_constants) - min(force_constants)) / fc_ref
                if spread > self.fc_ref_spread_warning:
                    per_state = ", ".join(
                        f"state {k + 1}: {fc:.1f}"
                        for k, fc in sorted(fc_per_state.items()))
                    self.ostream.print_warning(
                        f"Reactive bond {(pair[0] + 1, pair[1] + 1)} has force "
                        f"constants that differ by {100 * spread:.0f}% between "
                        f"the states it is bonded in ({per_state}). Using the "
                        f"average {fc_ref:.1f} as the reference force constant "
                        "for the dissociated-pair restraint; no single value "
                        "represents these well.")
        self.ostream.flush()
        return registry

    def _optimize_states(self, state_ffs, forming_bonds, breaking_bonds):
        """Run the MM conformational optimisation for every state.

        The pairs that are restrained for state k are those that react in an
        adjacent step but are not bonded in state k itself, i.e. exactly the
        forming bonds of the step starting at k and the breaking bonds of the
        step ending at k. For two states this is the reactant's forming bonds
        and the product's breaking bonds, as before.
        """

        if not self.optimize_ff:
            return

        n_states = len(state_ffs)
        any_change = any(len(bonds) > 0 for bonds in forming_bonds) or any(
            len(bonds) > 0 for bonds in breaking_bonds)
        if not any_change:
            self.ostream.print_info(
                "Skipping optimization of the force fields because no bonds are breaking or forming."
            )
            return

        for k, ff in enumerate(state_ffs):
            changing_bonds = set()
            if k < n_states - 1:
                changing_bonds |= set(forming_bonds[k])
            if k > 0:
                changing_bonds |= set(breaking_bonds[k - 1])
            changing_bonds = {
                pair for pair in changing_bonds if pair not in ff.bonds
            }

            if n_states == 2:
                note = 'reactant' if k == 0 else 'product'
            else:
                note = f'state {k + 1}'

            # TODO this optimisation can likely be taken care of by the openmmdynamics class
            ff.molecule = self._optimize_molecule(
                ff.molecule.get_element_ids(),
                ff,
                changing_bonds,
                note=note,
            )

    def _create_combined_forcefield(
        self,
        molecules: list[Molecule],
        partial_charges: list[list[float]]
        | list[float] | None,
        hessians: np.ndarray | list[np.ndarray | None] | None,
        name: str,
        total_multiplicity: int,
    ):
        if isinstance(molecules, Molecule):
            molecules = [molecules]

        if partial_charges is None:
            partial_charges = [None] * len(molecules)  # type: ignore
        elif isinstance(partial_charges[0], float) or isinstance(
                partial_charges[0], int):
            partial_charges = [partial_charges]  # type: ignore
        assert isinstance(partial_charges, list)

        if isinstance(hessians, np.ndarray) or hessians is None:
            hessians = [hessians] * len(molecules)

        assert_msg_critical(
            len(molecules) == len(partial_charges),
            "Amount of input molecules and lists of partial charges must match")
        assert_msg_critical(
            len(molecules) == len(hessians),
            "Amount of input molecules and lists of hessians must match")

        single_ffs: list[MMForceFieldGenerator] = []
        for i, (molecule, partial_charge,
                hessian) in enumerate(zip(molecules, partial_charges,
                                          hessians)):
            molecule_sanity_check(molecule)
            if partial_charge is not None:
                # Casting to float is necessary for json serialization
                assert isinstance(partial_charge, list) or isinstance(
                    partial_charge, np.ndarray)
                partial_charge = [float(x) for x in partial_charge]
                cond = len(partial_charge) == molecule.number_of_atoms()
                msg = f"Number of partial charges {len(partial_charge)} must match the number of atoms {molecule.number_of_atoms()} in the molecule."
                assert_msg_critical(cond, msg)

            self.ostream.print_blank()
            self.ostream.print_header(f"Building force field for {name}_{i+1}")
            self.ostream.print_blank()
            self.ostream.flush()
            ff = self._create_single_forcefield(
                molecule,
                partial_charge,
                hessian,
            )
            single_ffs.append(ff)

        total_charge = sum([mol.get_charge() for mol in molecules])

        self.ostream.print_info("Creating combined reactant force field")
        self.ostream.flush()
        combined_mol = self._combine_molecule(
            molecules,
            total_multiplicity,
            name,
        )
        self.ostream.print_info(
            f"Combined reactant with total charge {combined_mol.get_charge()} and multiplicity {combined_mol.get_multiplicity()}"
        )
        self.ostream.flush()
        combined_ff = self._combine_forcefield(single_ffs)
        combined_ff.molecule = combined_mol

        return combined_ff, single_ffs, total_charge

    def _create_single_forcefield(
        self,
        molecule,
        partial_charges: list[float] | None,
        hessian: np.ndarray | None = None,
    ) -> MMForceFieldGenerator:

        # Transforms input difctionary with keys 'molecule', 'charges', 'forcefield', 'optimize' into a forcefield generator
        # If the forcefield is not provided, it will be created from the molecule and charges
        # Calculates resp charges with some fallbacks, does optimization if specified, reparameterizes the forcefield if necessary

        # # If charges exist, load them into the forcefield object before creating the topology

        # The topology creation will calculate charges if they're not already set
        if self.optimize_mol:
            scf_drv = XtbDriver(ostream=self._sub_ostream(self.mute_scf))
            opt_drv = OptimizationDriver(scf_drv)
            if self.mute_scf:
                self.ostream.print_info("Optimising the geometry with xtb.")
                self.ostream.flush()
            opt_results = opt_drv.compute(molecule)

            charge = molecule.get_charge()
            multiplicity = molecule.get_multiplicity()
            molecule = Molecule.from_xyz_string(opt_results["final_geometry"])
            molecule.set_charge(charge)
            molecule.set_multiplicity(multiplicity)

        forcefield = MMForceFieldGenerator(ostream=self.ostream)

        forcefield.eq_param = False

        # Load or calculate the charges
        if partial_charges is not None:
            assert len(partial_charges) == molecule.number_of_atoms(
            ), "The number of provided charges does not match the number of atoms in the molecule"
            charge_sum = sum(partial_charges)
            if charge_sum - round(charge_sum) > 0.001:
                self.ostream.print_warning(
                    f"Sum of charges is {charge_sum} is not close to an integer. Confirm that the input is correct."
                )
            forcefield.partial_charges = partial_charges
            self.ostream.print_info("Creating topology")
            self.ostream.flush()
            forcefield.create_topology(molecule, water_model=self.water_model)
        else:
            if self.calculate_resp:
                if max(molecule.get_masses()) > 84:
                    basis = MolecularBasis.read(molecule,
                                                "STO-6G",
                                                ostream=None)
                    self.ostream.print_info(
                        f"Heavy ({max(molecule.get_masses())}) atom found. Using STO-6G basis (only comes in for RESP calculation)."
                    )
                else:
                    basis = MolecularBasis.read(molecule,
                                                "6-31G*",
                                                ostream=None)
                sub_ostream = self._sub_ostream(self.mute_scf)
                if molecule.get_multiplicity() == 1:
                    scf_drv = ScfRestrictedDriver(ostream=sub_ostream)
                else:
                    scf_drv = ScfUnrestrictedDriver(ostream=sub_ostream)

                self.ostream.print_info("Calculating SCF for RESP charges")
                self.ostream.flush()
                scf_results = scf_drv.compute(molecule, basis)
                if not scf_drv.is_converged:
                    self.ostream.print_warning(
                        "SCF did not converge, increasing convergence threshold to 1.0e-4 and maximum itterations to 200."
                    )
                    self.ostream.flush()
                    scf_drv.conv_thresh = 1.0e-4
                    scf_drv.max_iter = 200
                    scf_results = scf_drv.compute(molecule, basis)
                # self.ostream.unmute()
                assert scf_drv.is_converged, "SCF calculation for RESP charges did not converge, aborting"
                resp_drv = RespChargesDriver(
                    ostream=self._sub_ostream(self.mute_scf))
                self.ostream.flush()
                if self.mute_scf:
                    self.ostream.print_info("Calculating RESP charges")
                    self.ostream.flush()
                forcefield.partial_charges = resp_drv.compute(
                    molecule, basis, scf_results, 'resp')
                self.ostream.print_info(
                    f"RESP charges: {forcefield.partial_charges}")
                self.ostream.flush()
                self.ostream.print_info("Creating topology")
                forcefield.create_topology(
                    molecule,
                    basis,
                    scf_results=scf_results,
                    water_model=self.water_model,
                )
            else:
                self.ostream.print_info(
                    "Assigning partial charges based on electronegativity")
                partial_charges = molecule.get_partial_charges(
                    molecule.get_charge())
                forcefield.partial_charges = partial_charges
                self.ostream.print_info("Creating topology")
                self.ostream.flush()
                forcefield.create_topology(molecule,
                                           water_model=self.water_model)

        # Reparameterize the forcefield if necessary and requested
        # todo let this be handled by the MMforcefieldgenerator
        unknown_pairs = set()
        unknown_params = set()
        if self.reparameterize_bonds:
            for key, bond in forcefield.bonds.items():
                if bond['comment'] == 'Guessed':
                    sorted_tuple = tuple(sorted(key))
                    unknown_pairs.add(sorted_tuple)
                    unknown_params.add(key)
        if self.reparameterize_angles:
            for key, angle in forcefield.angles.items():
                if angle['comment'] == 'Guessed':
                    sorted_tuple = tuple(sorted((key[0], key[1])))
                    unknown_pairs.add(sorted_tuple)
                    sorted_tuple = tuple(sorted((key[1], key[2])))
                    unknown_pairs.add(sorted_tuple)
                    unknown_params.add(key)

        if len(unknown_pairs) > 0:
            self.ostream.print_info("Reparameterising force field.")
            self.ostream.flush()
            if hessian is None:
                self.ostream.print_info(
                    f"Calculating hessian submatrices for atom pairs {unknown_pairs} to reparameterise the force field."
                )
                sub_ostream = self._sub_ostream(self.mute_scf)
                if molecule.get_multiplicity() == 1:
                    scf_drv = ScfRestrictedDriver(ostream=sub_ostream)
                else:
                    scf_drv = ScfUnrestrictedDriver(ostream=sub_ostream)
                self.ostream.flush()
                basis = MolecularBasis.read(molecule, self.hessian_basis)
                scf_drv.xcfun = self.hessian_xc_fun
                scf_drv.dispersion = True
                scf_drv.compute(molecule, basis)
                if scf_drv.is_converged is False:
                    self.ostream.print_warning(
                        "SCF did not converge, increasing convergence threshold to 1.0e-4 and maximum itterations to 200."
                    )
                    self.ostream.flush()
                    scf_drv.conv_thresh = 1.0e-4
                    scf_drv.max_iter = 200
                    scf_drv.compute(molecule, basis)
                assert scf_drv.is_converged, "SCF calculation for Hessian did not converge, aborting"

                # inherits scf_drv's (already silent) output stream
                hess_drv = ScfHessianDriver(scf_drv)
                hess_drv.atom_pairs = list(unknown_pairs)
                hess_drv.compute(molecule, basis)
                hessian = np.copy(hess_drv.hessian)  # type: ignore
            else:
                cond = np.shape(hessian) == (molecule.number_of_atoms() * 3,
                                             molecule.number_of_atoms() * 3)
                msg = f"Hessian shape {np.shape(hessian)} should be square with width 3 times the number"
                msg += f" of atoms 3*{molecule.number_of_atoms()}={3*molecule.number_of_atoms()} in the molecule."
                assert_msg_critical(cond, msg)
            self.ostream.flush()
            forcefield.reparameterize(hessian,
                                      reparameterize_keys=unknown_params)
        return forcefield

    def _combine_molecule(self, molecules, total_multiplicity, name='MOL'):

        # Guesses how to combine molecular structures without overlapping them

        combined_molecule = Molecule()
        # pos = []
        charge = 0
        Sm1 = 0
        for i, mol in enumerate(molecules):
            charge += mol.get_charge()
            Sm1 += mol.get_multiplicity() - 1
            if combined_molecule.number_of_atoms() > 0:
                coords = combined_molecule.get_coordinates_in_angstrom()
                max_x = max(coords[:, 0])
                min_x = min(coords[:, 0])
                shift = max_x - min_x + 2
                self.ostream.print_info(
                    f"max_x: {max_x}, min_x: {min_x}, shifting next molecule by {shift} angstrom"
                )
            else:
                shift = 0
            coords = mol.get_coordinates_in_angstrom()
            min_x = min(coords[:, 0])
            coords[:, 0] -= min_x
            for elem, coord in zip(mol.get_element_ids(), coords):
                coord[0] += shift
                # pos.append(coord)
                combined_molecule.add_atom(int(elem), Point(coord), 'angstrom')
        combined_molecule.set_charge(charge)
        if total_multiplicity > -1:
            combined_molecule.set_multiplicity(total_multiplicity)
        else:
            self.ostream.print_info(
                f"Setting multiplicity of the combined molecule to {Sm1 + 1} based on the multiplicities of the provided molecules."
            )
            self.ostream.flush()

            combined_molecule.set_multiplicity(Sm1 + 1)

        molecule_sanity_check(combined_molecule)
        return combined_molecule

    def _match_reactant_and_product(
        self,
        reactant: Molecule,
        product: Molecule,
        breaking_bonds: set[tuple[int, int]],
        forming_bonds: set[tuple[int, int]],
    ):

        # Match the indices of the reactant and product forcefield generators

        assert reactant.number_of_atoms() == product.number_of_atoms(
        ), "The number of atoms in the reactant and product do not match"
        # Turn the reactand and product into graphs

        rm = ReactionMatcher(ostream=self.ostream)
        total_mapping, breaking_bonds, forming_bonds = rm.get_mapping(
            reactant,
            product,
            breaking_bonds,
            forming_bonds,
        )  # type: ignore
        if total_mapping is None:
            raise ValueError(
                "Could not find a mapping between the reactant and product force fields."
            )
        total_mapping = {v: k for k, v in total_mapping.items()}
        print_mapping = {k + 1: v + 1 for k, v in total_mapping.items()}
        self.ostream.print_info(f"Mapping: {print_mapping}")
        self.ostream.flush()
        return total_mapping

    def _combine_forcefield(
            self,
            forcefields: list[MMForceFieldGenerator]) -> MMForceFieldGenerator:

        # Merge a list of forcefield generators into a single forcefield generator while taking care of the atom indices

        forcefield = MMForceFieldGenerator(ostream=self.ostream)
        forcefield.atoms = {}
        forcefield.bonds = {}
        forcefield.angles = {}
        forcefield.dihedrals = {}
        forcefield.impropers = {}
        atom_count = 0
        forcefield.unique_atom_types = []
        forcefield.pairs = {}
        forcefield.atom_info_dict = {}
        for l, ff in enumerate(forcefields):
            # Shift all atom keys by the current atom count so that every atom has a unique ID
            shift = atom_count
            mapping = {atom_key: atom_key + shift for atom_key in ff.atoms}
            ReactionForceFieldBuilder._apply_mapping_to_forcefield(ff, mapping)
            atom_count += len(ff.atoms)
            for atom in ff.atoms.values():
                atom['name'] = f"{atom['name']}{l+1}"
            forcefield.atom_info_dict.update(ff.atom_info_dict)
            forcefield.atoms.update(ff.atoms)
            forcefield.bonds.update(ff.bonds)
            forcefield.angles.update(ff.angles)
            forcefield.dihedrals.update(ff.dihedrals)
            forcefield.impropers.update(ff.impropers)
            if hasattr(ff, 'unique_atom_types'):
                forcefield.unique_atom_types += ff.unique_atom_types
            if hasattr(ff, 'pairs'):
                forcefield.pairs.update(ff.pairs)

        return forcefield

    @staticmethod
    def _apply_mapping_to_forcefield(
            forcefield: MMForceFieldGenerator,
            mapping: dict[int, int]) -> MMForceFieldGenerator:

        # Remap indices in the forcefield to the new indices

        new_ff_atoms = {}
        for atom_key in forcefield.atoms:
            key = mapping[atom_key]
            val = forcefield.atoms[atom_key]
            new_ff_atoms.update({key: val})

        new_atom_info = {}
        for atom_info_key in forcefield.atom_info_dict.keys():
            key = mapping[atom_info_key - 1] + 1
            val = forcefield.atom_info_dict[atom_info_key]
            val['ConnectedAtomsNumbers'] = [
                mapping[key - 1] + 1 for key in val['ConnectedAtomsNumbers']
            ]
            val['AtomNumber'] = mapping[val['AtomNumber'] - 1] + 1
            new_atom_info.update({key: val})

        # Sort the atoms by index
        forcefield.atoms = dict(sorted(new_ff_atoms.items()))
        forcefield.atom_info_dict = dict(sorted(new_atom_info.items()))

        forcefield.bonds = ReactionForceFieldBuilder._apply_mapping_to_parameters(
            forcefield.bonds, mapping)
        forcefield.angles = ReactionForceFieldBuilder._apply_mapping_to_parameters(
            forcefield.angles, mapping)
        forcefield.dihedrals = ReactionForceFieldBuilder._apply_mapping_to_parameters(
            forcefield.dihedrals, mapping)
        forcefield.impropers = ReactionForceFieldBuilder._apply_mapping_to_parameters(
            forcefield.impropers, mapping)
        return forcefield

    @staticmethod
    def _apply_mapping_to_molecule(molecule, mapping):
        new_molecule = Molecule()
        positions = molecule.get_coordinates_in_angstrom()
        element_ids = molecule.get_element_ids()
        sorted_ids = dict(sorted(mapping.items(),
                                 key=lambda item: item[1])).keys()
        for id in sorted_ids:
            new_molecule.add_atom(int(element_ids[id]), Point(positions[id]),
                                  'angstrom')
        return new_molecule

    @staticmethod
    def _apply_mapping_to_parameters(
            old_parameters: dict[tuple, dict],
            mapping: dict[int, int]) -> dict[tuple, dict]:

        # Remap the indices in a specific set of parameters

        new_parameters = {}
        for old_key in old_parameters:
            new_key = tuple([mapping[atom_key] for atom_key in old_key])
            val = old_parameters[old_key]

            # Make sure that the new key is still properly ordered§
            if new_key[-1] < new_key[0]:
                new_key = new_key[::-1]
            new_parameters.update({new_key: val})
        new_parameters = dict(sorted(new_parameters.items()))
        return new_parameters

    def _summarise_reaction(self, reactant, product, step_label=None):
        """
        Summarises the reaction by printing the bonds that are being broken and formed.

        Returns:
            None
        """
        reactant_bonds = set(reactant.bonds)
        product_bonds = set(product.bonds)
        formed_bonds = product_bonds - reactant_bonds
        broken_bonds = reactant_bonds - product_bonds
        if step_label is None:
            self.ostream.print_header("Reaction summary")
        else:
            self.ostream.print_header(f"Reaction summary for {step_label}")
        self.ostream.print_header(f"{len(broken_bonds)} breaking bonds:")

        if len(broken_bonds) > 0:
            self.ostream.print_header(
                "ReaType  ProType  ID - ReaType  ProType  ID")
        for bond_key in broken_bonds:
            reactant_type0 = reactant.atoms[bond_key[0]]["type"]
            product_type0 = product.atoms[bond_key[0]]["type"]
            id0 = bond_key[0] + 1
            reactant_type1 = reactant.atoms[bond_key[1]]["type"]
            product_type1 = product.atoms[bond_key[1]]["type"]
            id1 = bond_key[1] + 1
            self.ostream.print_header(
                f"{reactant_type0:^9}{product_type0:^9}{id0:^2} - {reactant_type1:^9}{product_type1:^9}{id1:^2}"
            )
        self.ostream.print_blank()
        self.ostream.print_header(f"{len(formed_bonds)} forming bonds:")
        if len(formed_bonds) > 0:
            self.ostream.print_header(
                "ReaType  ProType  ID - ReaType  ProType  ID")
        for bond_key in formed_bonds:
            reactant_type0 = reactant.atoms[bond_key[0]]["type"]
            product_type0 = product.atoms[bond_key[0]]["type"]
            id0 = bond_key[0] + 1
            reactant_type1 = reactant.atoms[bond_key[1]]["type"]
            product_type1 = product.atoms[bond_key[1]]["type"]
            id1 = bond_key[1] + 1
            self.ostream.print_header(
                f"{reactant_type0:^9}{product_type0:^9}{id0:^2} - {reactant_type1:^9}{product_type1:^9}{id1:^2}"
            )
        self.ostream.print_blank()
        self.ostream.flush()
        return formed_bonds, broken_bonds

    def _optimize_molecule(
        self,
        elemental_ids,
        forcefield,
        changing_bonds,
        name='MOL',
        note=None,
    ):

        # Does an FF optimization of the molecule.

        # merging of systems through openmm files is shaky, as it depends on the atom naming working. See atom renaming in combine_forcefield
        # todo find a better way to do this
        forcefield.write_openmm_files(name, name)

        # in case one wants to inspect reactant.itp and product.itp
        # forcefield.write_gromacs_files(note, name)

        # for bond in changing_bonds:
        #     forcefield.bonds.pop(bond)

        pdb = mmapp.PDBFile(f'{name}.pdb')
        ff = mmapp.ForceField(f'{name}.xml')
        sys_xml_path = Path(f'{name}_sys.xml')

        modeller = mmapp.Modeller(pdb.topology, pdb.positions)

        top = modeller.getTopology()
        # pos = modeller.getPositions()

        exiting = False
        while not exiting:
            try:
                mmsys = ff.createSystem(
                    top,
                    nonbondedMethod=mmapp.CutoffNonPeriodic,
                    nonbondedCutoff=1.0 * mmunit.nanometers,
                )

                if self.mm_opt_constrain_bonds:
                    mmsys = self._add_reaction_bonds(forcefield, mmsys,
                                                     changing_bonds, note)

                with sys_xml_path.open('w') as f:
                    f.write(mm.XmlSerializer.serialize(mmsys))

                opm_dyn = OpenMMDynamics(ostream=self._sub_ostream(True))
                opm_dyn.openmm_platform = "CPU"

                opm_dyn.pdb = pdb
                opm_dyn.system = mmsys

                self.ostream.print_info(
                    f"Running conformational sampling with {self.optimize_steps*self.optimize_conformer_snapshots} steps and {self.optimize_conformer_snapshots} snapshots at {self.optimize_temp} K for {note} molecule."
                )
                self.ostream.print_info(
                    "This can be turned off with optimize_ff = False.")
                self.ostream.flush()
                conformers_dict = opm_dyn.conformational_sampling(
                    ensemble='NVT',
                    nsteps=self.optimize_steps *
                    self.optimize_conformer_snapshots,
                    snapshots=self.optimize_conformer_snapshots,
                    temperature=self.optimize_temp,
                )

                min_arg = np.argmin(conformers_dict['energies'])
                new_molecule = conformers_dict['molecules'][min_arg]
                self.ostream.print_info(
                    f"Found {len(conformers_dict['molecules'])} conformers with energies {conformers_dict['energies']} during optimization of the {note} molecule."
                )
                self.ostream.flush()
                exiting = True

            except Exception as e:
                print_exception_if_debug()
                if self.optimize_dist_restraint_offset < 2.5:
                    self.optimize_dist_restraint_offset += 0.5
                    self.ostream.print_warning(
                        f"OpenMM optimization of the {note} molecule failed with error: {e}. Increasing distance restraint offset to {self.optimize_dist_restraint_offset} angstrom and retrying."
                    )
                    self.ostream.flush()
                else:
                    self.ostream.print_warning(
                        f"OpenMM optimization of the {note} molecule failed with error: {e}. Reverting to original geometry."
                    )
                    self.ostream.flush()
                    new_molecule = forcefield.molecule
                    exiting = True
        self.optimize_dist_restraint_offset = 0.5  # Reset for next use
        new_molecule.set_charge(forcefield.molecule.get_charge())
        new_molecule.set_multiplicity(forcefield.molecule.get_multiplicity())
        Path(f'{name}.xml').unlink()
        Path(f'{name}.pdb').unlink()
        sys_xml_path.unlink()
        return new_molecule

    def _add_reaction_bonds(self, forcefield, mmsys, changing_bonds, note):
        self.ostream.print_info(
            "Guessing intra-molecular constraints based on breaking bonds. This option can be turned off with mm_opt_constrain_bonds."
        )

        dist_restraint_expr = "0.5 * k * (r - r0)^2*step(r-r0)"
        dist_restraint_force = mm.CustomBondForce(dist_restraint_expr)
        dist_restraint_force.addPerBondParameter("r0")
        dist_restraint_force.addGlobalParameter("k",
                                                self.optimize_dist_restraint_k)

        for bond in changing_bonds:
            s1 = forcefield.atoms[bond[0]]['sigma']
            s2 = forcefield.atoms[bond[1]]['sigma']
            s = (s1 + s2) / 2

            r0 = s * (
                2**(1 / 6)
            ) + self.optimize_dist_restraint_offset * 0.1  # rmin = sigma * 2^(1/6)
            dist_restraint_force.addBond(bond[0], bond[1], [r0])
            self.ostream.print_info(
                f"Adding distance restraint for atoms {tuple(b+1 for b in bond)} with r0 {r0:.3f} and k {self.optimize_dist_restraint_k} for MM equilibration of the {note}"
            )
            self.ostream.flush()

        mmsys.addForce(dist_restraint_force)
        return mmsys
