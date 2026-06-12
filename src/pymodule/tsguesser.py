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
import sys
import numpy as np
import math
import copy
import h5py
import time
from random import getrandbits
from pathlib import Path

from .veloxchemlib import mpi_master, hartree_in_kjpermol, bohr_in_angstrom
from .outputstream import OutputStream
from .openmmdynamics import OpenMMDynamics
from .errorhandler import assert_msg_critical
from .mmforcefieldgenerator import MMForceFieldGenerator
from .molecule import Molecule
from .scfrestdriver import ScfRestrictedDriver
from .molecularbasis import MolecularBasis
from .reaffbuilder import ReactionForceFieldBuilder
from .evbsystembuilder import EvbSystemBuilder

try:
    import openmm as mm
    import openmm.app as mmapp
    import openmm.unit as mmunit
except ImportError:
    pass

# All positions are in Angsrom unless otherwise stated


class TransitionStateGuesser():

    def __init__(self, comm=None, ostream=None):
        '''
        Initialize the Transition State Guesser class.
        '''
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

        self.molecule = None
        self.results = {}
        self.lambda_vector = list(np.round(np.linspace(0, 1, 21), 3))
        self.mute_ff_build = True

        # File I/O controls
        # Set save_intermediates=True to write reactant/product JSON and
        # system XML files to disk — useful for debugging but off by default
        # to avoid unnecessary file I/O on HPC shared filesystems.
        self.save_intermediates = False
        # Set save_results_file=True (default) to write the HDF5 results file
        # after each scan. This enables crash recovery and load_results().
        self.save_results_file = True

        timing_str = str(int(time.time()))
        self.folder_name = f'ts_data_{timing_str}'
        self.results_file = f'ts_results_{timing_str}.h5'

        self.mm_temperature = 600
        self.mm_steps = 1000
        self.conformer_snapshots = 10
        self.mm_step_size = 0.001
        self.save_mm_traj = False
        self.force_conformer_search = False
        self.discont_conformer_search = False
        self.peak_conformer_search = False
        self.peak_conformer_search_range = 2
        self.mm_conformer_equivalence_threshold = 1e-1  # kJ/mol
        self.mm_scan_backward = False

        # Force constant (kJ/mol) for the periodic dihedral restraint used in
        # the integration systems when a conformational TS is detected (no
        # forming or breaking bonds).  The expression k*(1-cos(theta-theta0))
        # is used, which is periodic and harmonic-like near the minimum.
        self.conformer_k: float = 200.0

        # Optional explicit active torsion (1-indexed atom indices) for the
        # conformational TS scan.  When set, auto-detection via
        # _detect_active_dihedral is skipped and this dihedral is used instead.
        # Must be a 4-element sequence [i, j, k, l] with 1-based atom indices.
        self.active_torsion: tuple | None = None

        # Set by build_systems when a conformational TS is detected; used for
        # descriptive output during the scan.
        self._conformer_active_torsion: tuple | None = None
        self._conformer_phi_reactant: float | None = None
        self._conformer_phi_product: float | None = None

        # Implicit solvation during conformational sampling.
        # Set implicit_solvent_model to one of 'gbn', 'gbn2', 'obc1', 'obc2', 'hct' to enable GB solvation; None runs in vacuum (default).
        self.implicit_solvent_model: str | None = None
        self.solute_dielectric: float = 1.0
        self.solvent_dielectric: float = 78.39

        # Named solvent for the QM SMD calculation, used automatically when
        # implicit_solvent_model is set and do_qm_scan is True.
        # Must be a solvent name recognised by VeloxChem's SMD driver
        # (e.g. 'water', 'ethanol', 'acetonitrile'). Defaults to 'water'.
        # Note: XTB does not support SMD; a warning is emitted in that case.
        self.smd_solvent: str = 'water'

        self.scf_drv = None
        self.qm_xcfun = "PBE0"
        self.qm_basis = 'def2-svp'
        self.do_qm_scan = False
        self.max_qm_conformers = 5
        self.mute_scf = True

        self.sys_builder_configuration = {
            "name": "vacuum",
            "bonded_integration": True,
            "soft_core_coulomb_pes": True,
            "soft_core_lj_pes": True,
            "soft_core_coulomb_int": False,
            "soft_core_lj_int": False,
        }

        self._reaction_matcher_assist_min_depth = None

        self.ffbuilder = ReactionForceFieldBuilder(ostream=self.ostream)
        self.ffbuilder.calculate_resp = False

    def find_transition_state(
        self,
        reactant: Molecule | list[Molecule],
        product: Molecule | list[Molecule],
        constraints=None,
        **build_forcefields_kwargs,
    ):
        """Find a guess for the transition state using a force field scan.

        Args:
            reactant (Molecule | list[Molecule]): The reactant molecule or a list of reactant molecules.
            product (Molecule | list[Molecule]): The product molecule or a list of product molecules.
            constraints: list of constraints to be applied during the scan.
            **build_forcefields_kwargs: Additional keyword arguments to be passed to the force field builder.

        Returns:
            dict: molecule object of the guessed transition state and a dictionary with the results of the scan.
        """
        self.results = {}
        # Build forcefields and systems
        self.ffbuilder.calculate_resp = self.implicit_solvent_model is not None
        self.build_forcefields(reactant, product, **build_forcefields_kwargs)
        self.build_systems(constraints)

        # Scan MM
        self.scan_mm()

        # Scan QM
        if self.do_qm_scan:
            self.scan_qm(self.results)

        return self.results

    def build_forcefields(self, reactant, product, **build_forcefields_kwargs):
        if self.mute_ff_build:
            self.ostream.print_info(
                "Building forcefields. Disable mute_ff_build to see detailed output."
            )
            self.ostream.flush()
            self.ostream.mute()

        if self._reaction_matcher_assist_min_depth is not None:
            self.ffbuilder._reaction_matcher_assist_min_depth = int(self._reaction_matcher_assist_min_depth)

        self.reactant, self.product, self.forming_bonds, self.breaking_bonds, reactants, products, product_mapping = self.ffbuilder.build_forcefields(
            reactant=reactant,
            product=product,
            **build_forcefields_kwargs,
        )

        if self.mute_ff_build:
            self.ostream.unmute()
            self.ffbuilder._summarise_reaction(self.reactant, self.product)

        self.molecule = Molecule.read_xyz_string(
            self.reactant.molecule.get_xyz_string())
        self.molecule.set_charge(self.reactant.molecule.get_charge())
        self.molecule.set_multiplicity(
            self.reactant.molecule.get_multiplicity())

        self.ostream.print_info(
            f"System has charge {self.molecule.get_charge()} and multiplicity {self.molecule.get_multiplicity()}. Provide correct values if this is wrong."
        )
        self.mol_charge = self.molecule.get_charge()
        self.mol_multiplicity = self.molecule.get_multiplicity()

        if self.save_intermediates:
            Path(self.folder_name).mkdir(parents=True, exist_ok=True)
            self.ostream.print_info(
                f"Saving reactant and product forcefield as json to {self.folder_name}"
            )
            self.ostream.flush()
            self.reactant.save_forcefield_as_json(
                self.reactant, str(Path(self.folder_name) / "reactant.json"))
            self.product.save_forcefield_as_json(
                self.product, str(Path(self.folder_name) / "product.json"))

        rea_bonds = set(self.reactant.bonds.keys())
        pro_bonds = set(self.product.bonds.keys())
        static_bonds = rea_bonds & pro_bonds
        self.results.update({
            'breaking_bonds': self.breaking_bonds,
            'forming_bonds': self.forming_bonds,
            'static_bonds': static_bonds,
            'reactant': self.reactant,
            'product': self.product,
        })

        return self.results

    def build_systems(self, constraints=None):

        self.lambda_vector = [round(l, 3) for l in self.lambda_vector]
        self.ostream.print_info(
            f"Rounding lambda vector to 3 decimal places: {self.lambda_vector}")

        sysbuilder = EvbSystemBuilder()
        if self.mute_ff_build:
            sysbuilder.ostream.mute()
            self.ostream.print_info(
                "Building MM systems for the transition state guess. Disable mute_ff_build to see detailed output."
            )
            self.ostream.flush()

        configuration = dict(self.sys_builder_configuration)
        if self.implicit_solvent_model is not None:
            configuration[
                'implicit_solvent_model'] = self.implicit_solvent_model
            configuration['solute_dielectric'] = self.solute_dielectric
            configuration['solvent_dielectric'] = self.solvent_dielectric

        # Conformational TS: no bonds forming or breaking — replace the
        # active dihedral in the integration systems with a periodic restraint
        # whose minimum shifts from phi_reactant (λ=0) to phi_product (λ=1).
        if len(self.forming_bonds) == 0 and len(self.breaking_bonds) == 0:
            self.ostream.print_info(
                "No forming or breaking bonds detected. "
                "Treating as a conformational transition state: "
                "using a shifting dihedral restraint in the integration systems."
            )
            self.ostream.flush()
            if self.active_torsion is not None:
                assert_msg_critical(
                    len(self.active_torsion) == 4,
                    'TransitionStateGuesser: active_torsion must contain '
                    'exactly four 1-based atom indices')
                reactant_natoms = self.reactant.molecule.number_of_atoms()
                assert_msg_critical(
                    all(1 <= a <= reactant_natoms for a in self.active_torsion),
                    'TransitionStateGuesser: active_torsion indices must be '
                    f'between 1 and {reactant_natoms}')

                # Convert user-supplied 1-indexed tuple to 0-indexed
                torsion_0idx = tuple(a - 1 for a in self.active_torsion)
                one_based = list(self.active_torsion)
                phi_reactant = self.reactant.molecule.get_dihedral_in_degrees(
                    one_based)
                phi_product = self.product.molecule.get_dihedral_in_degrees(
                    one_based)
                # Wrap product angle to shortest path from reactant
                delta = phi_product - phi_reactant
                if delta > 180.0:
                    delta -= 360.0
                elif delta <= -180.0:
                    delta += 360.0
                self.ostream.print_info(
                    f"Using explicitly set active torsion (1-indexed): "
                    f"{tuple(self.active_torsion)}, "
                    f"phi_reactant = {phi_reactant:.1f}°, "
                    f"phi_product = {phi_product:.1f}°, "
                    f"|Δφ| = {abs(delta):.1f}°.")
                self.ostream.flush()
                active_torsion = torsion_0idx
            else:
                active_torsion, phi_reactant, phi_product = self._detect_active_dihedral(
                )
            self._conformer_active_torsion = active_torsion
            self._conformer_phi_reactant = float(phi_reactant)
            self._conformer_phi_product = float(phi_product)
            configuration['conformer_active_torsion'] = active_torsion
            configuration['conformer_phi_reactant'] = float(phi_reactant)
            configuration['conformer_phi_product'] = float(phi_product)
            configuration['conformer_k'] = float(self.conformer_k)
            self.ostream.print_info(
                f"Conformer restraint force constant: {self.conformer_k:.1f} kJ/mol."
            )
            self.ostream.flush()

        self.systems, self.topology, _ = sysbuilder.build_systems(
            self.reactant,
            self.product,
            list(self.lambda_vector),
            configuration,
            constraints,
        )
        if self.mute_ff_build:
            self.ostream.print_blank()
            self.ostream.flush()
            sysbuilder.ostream.unmute()
        if self.save_intermediates:
            systems_dir_path = Path(self.folder_name) / "systems"
            systems_dir_path.mkdir(parents=True, exist_ok=True)
            systems_dir = str(systems_dir_path)
            self.ostream.print_info(f"Saving systems as xml to {systems_dir}")
            self.ostream.flush()
            sysbuilder.save_systems_as_xml(self.systems, systems_dir)
        self.results.update({'lambda_vec': self.lambda_vector})

    def scan_mm(self):

        self.folder = Path.cwd() / self.folder_name

        # pdbs are saved in angstrom

        rea_int = mm.VerletIntegrator(1)
        rea_sim = mmapp.Simulation(
            self.topology,
            self.systems['reactant'],
            rea_int,
        )
        pro_int = mm.VerletIntegrator(1)
        pro_sim = mmapp.Simulation(
            self.topology,
            self.systems['product'],
            pro_int,
        )

        # pos in angstrom
        # pos = self.initial_positions
        rea_init_pos = self.reactant.molecule.get_coordinates_in_angstrom()
        pro_init_pos = self.product.molecule.get_coordinates_in_angstrom()
        scan_dict = {}
        try:
            if self.force_conformer_search:
                self.ostream.print_info(
                    "force_conformer_search true. Doing conformer search at every lambda."
                )
                self.ostream.flush()
                scan_dict = self._run_mm_scan(
                    self.lambda_vector,
                    rea_sim,
                    pro_sim,
                    conformer_search=True,
                    forward_init_pos=rea_init_pos,
                    backward_init_pos=pro_init_pos,
                )
            else:
                scan_dict = self._run_mm_scan(
                    self.lambda_vector,
                    rea_sim,
                    pro_sim,
                    conformer_search=False,
                    forward_init_pos=rea_init_pos,
                    backward_init_pos=pro_init_pos,
                )

                #     # Find peak
                searched_conformers_indices = []
                V, E1, E2, conf_indices = self._get_best_mm_E_from_scan_dict(
                    scan_dict)
                if self.peak_conformer_search:
                    peak_iteration = 0

                    while True:
                        peak_iteration += 1
                        peak_index = int(np.argmax(V))
                        peak_lambda = self.lambda_vector[peak_index]

                        # Stop as soon as the peak itself has already been
                        # conformer-searched. Checking neighbours is not
                        # sufficient: the peak can fall in a gap between two
                        # disjoint search windows whose edges are both marked.
                        if (peak_index in searched_conformers_indices
                                and (max(0, peak_index - 1)
                                     in searched_conformers_indices)
                                and (min(peak_index + 1,
                                         len(self.lambda_vector) - 1)
                                     in searched_conformers_indices)):
                            break

                        min_index = max(
                            0,
                            peak_index - self.peak_conformer_search_range,
                        )
                        max_index = min(
                            len(self.lambda_vector) - 1,
                            peak_index + self.peak_conformer_search_range,
                        )

                        # Only scan indices not yet covered by a prior
                        # iteration to avoid redundant MD runs.
                        new_indices = [
                            i for i in range(min_index, max_index + 1)
                            if i not in searched_conformers_indices
                        ]
                        if not new_indices:
                            # Entire window already searched but peak not
                            # recorded — guard against infinite loop.
                            break

                        new_lambdas = [
                            self.lambda_vector[i] for i in new_indices
                        ]

                        self.ostream.print_info(
                            f"Found peak MM E: {V[peak_index]:.3f} at Lambda: {peak_lambda}"
                            f" (iteration {peak_iteration}).")
                        self.ostream.print_info(
                            f"Doing conformer search from Lambda: {new_lambdas[0]} to Lambda: {new_lambdas[-1]}."
                        )
                        self.ostream.flush()

                        forward_init_pos = scan_dict[new_lambdas[0]][0]['pos']
                        backward_init_pos = scan_dict[new_lambdas[-1]][0]['pos']

                        scan_dict_peak_conf = self._run_mm_scan(
                            new_lambdas,
                            rea_sim,
                            pro_sim,
                            conformer_search=True,
                            forward_init_pos=forward_init_pos,
                            backward_init_pos=backward_init_pos,
                            skip_backward=True,
                        )
                        for l in scan_dict_peak_conf.keys():
                            scan_dict[l] += scan_dict_peak_conf[l]

                        searched_conformers_indices.extend(new_indices)
                        searched_conformers_indices = sorted(
                            list(set(searched_conformers_indices)))

                        # Re-evaluate so the next iteration and the
                        # discontinuity check both see up-to-date energies.
                        V, E1, E2, conf_indices = (
                            self._get_best_mm_E_from_scan_dict(scan_dict))

                if self.discont_conformer_search:
                    discont_indices = self._check_discontinuities(E1, E2)

                    while len(discont_indices) > 0 and len(
                            searched_conformers_indices) < len(
                                self.lambda_vector):
                        self.ostream.flush()
                        self.ostream.print_info(
                            f"Found discontinuities at indices: {discont_indices}."
                        )
                        to_search_indices = sorted([
                            i for i in discont_indices
                            if i not in searched_conformers_indices
                        ])
                        searched_conformers_indices.extend(to_search_indices)
                        searched_conformers_indices = sorted(
                            list(set(searched_conformers_indices)))
                        to_search_lambda = [
                            self.lambda_vector[i] for i in to_search_indices
                        ]
                        self.ostream.print_info(
                            f"Performing conformer search at lambda values: {to_search_lambda}."
                        )
                        self.ostream.flush()
                        forward_init_pos = scan_dict[self.lambda_vector[
                            to_search_indices[0]]][0]['pos']
                        backward_init_pos = scan_dict[self.lambda_vector[
                            to_search_indices[-1]]][0]['pos']
                        scan_dict_discont_conf = self._run_mm_scan(
                            to_search_lambda,
                            rea_sim,
                            pro_sim,
                            conformer_search=True,
                            forward_init_pos=forward_init_pos,
                            backward_init_pos=backward_init_pos,
                            skip_backward=True,
                        )
                        for l in scan_dict_discont_conf.keys():
                            scan_dict[l] += scan_dict_discont_conf[l]
                        V, E1, E2, conf_indices = self._get_best_mm_E_from_scan_dict(
                            scan_dict)
                        discont_indices = self._check_discontinuities(E1, E2)
        except Exception:
            err_str = "The MM scan crashed. Saving results in self.results and raising exception"
            self.ostream.print_warning(err_str)
            self.ostream.flush()
            self.results.update({'scan': scan_dict})
            raise

        max_mm_energy = None
        for i, (l, mm_result) in enumerate(scan_dict.items()):
            min_local_E = None
            for j, conformer in enumerate(mm_result):
                if min_local_E is None or conformer['v'] < min_local_E:
                    min_local_E = conformer['v']
                    min_local_conformer_index = j

            if max_mm_energy is None or min_local_E > max_mm_energy:
                max_mm_xyz = mm_result[min_local_conformer_index]['xyz']
                max_mm_energy = min_local_E
                max_mm_lambda = l
                min_mm_conformer_index = min_local_conformer_index

        self.ostream.print_info(
            f"Found highest MM E: {max_mm_energy:.3f} at Lambda: {max_mm_lambda} and conformer index: {min_mm_conformer_index}."
        )
        self.ostream.print_blank()
        self.results.update({
            'scan': scan_dict,
            'max_mm_xyz': max_mm_xyz,
            'max_mm_lambda': max_mm_lambda,
            'min_mm_conformer_index': min_mm_conformer_index,
        })
        self.molecule = Molecule.read_xyz_string(max_mm_xyz)
        self.molecule.set_multiplicity(self.mol_multiplicity)
        self.molecule.set_charge(self.mol_charge)
        if self.save_results_file:
            self.save_results(self.results_file, self.results)
        return self.results

    def _check_discontinuities(self, E1, E2):
        discont_indices = []
        for i, e1 in enumerate(E1[:-1]):
            if e1 > E1[i + 1]:
                discont_indices.append(i)
                discont_indices.append(i + 1)

        for i, e2 in enumerate(E2[:-1]):
            if e2 < E2[i + 1]:
                discont_indices.append(i)
                discont_indices.append(i + 1)
        discont_indices = set(sorted(discont_indices))
        return discont_indices

    def _run_mm_scan(self,
                     lambda_vals,
                     rea_sim,
                     pro_sim,
                     conformer_search,
                     forward_init_pos,
                     backward_init_pos,
                     skip_backward=False):
        pos = copy.copy(forward_init_pos)
        results = {}
        self._print_mm_header(conformer_search=conformer_search)
        for l in lambda_vals:

            result = self._get_mm_energy(
                self.topology,
                self.systems[l],
                l,
                pos,
                rea_sim,
                pro_sim,
                conformer_search,
            )
            results[l] = result
            arg = np.argmin([res['v'] for res in result])
            e1 = result[arg]['e1']
            e2 = result[arg]['e2']
            v = result[arg]['v']
            # e_int = result[arg]['e_int']
            pos = result[arg]['pos']
            n_conf = len(result)

            self._print_mm_iter(l, e1, e2, v, n_conf)

        if self.mm_scan_backward and not skip_backward:
            self.ostream.print_info(
                "mm_scan_backward turned on. Scanning in reverse direction.")
            self.ostream.flush()
            lambda_vals_rev = list(reversed(lambda_vals))
            self._print_mm_header(conformer_search=conformer_search)
            pos = copy.copy(backward_init_pos)
            for l in lambda_vals_rev:
                result = self._get_mm_energy(
                    self.topology,
                    self.systems[l],
                    l,
                    pos,
                    rea_sim,
                    pro_sim,
                    conformer_search,
                )
                results[l] += result
                arg = np.argmin([res['v'] for res in result])
                e1 = result[arg]['e1']
                e2 = result[arg]['e2']
                v = result[arg]['v']
                # e_int = result[arg]['e_int']
                pos = result[arg]['pos']
                n_conf = len(result)
                self._print_mm_iter(l, e1, e2, v, n_conf)
                self.ostream.flush()

        self.ostream.print_blank()

        for l, result in results.items():
            # remove duplicate conformers with identical MM energy 'v'
            unique_confs = []
            seen_v = set()
            for conf in result:
                v_val = conf['v']
                seen = False
                for v in seen_v:
                    if abs(v - v_val) < self.mm_conformer_equivalence_threshold:
                        seen = True
                if not seen:
                    seen_v.add(v_val)
                    unique_confs.append(conf)
            original_n_conf = len(result)
            new_n_conf = len(unique_confs)
            if original_n_conf != new_n_conf:
                self.ostream.print_info(
                    f"Lambda {l}: Reduced {original_n_conf} conformers to {new_n_conf} unique conformers using equivalence threshold of {self.mm_conformer_equivalence_threshold} kJ/mol."
                )
                self.ostream.flush()
            results[l] = unique_confs

        return results

    def _get_mm_energy(
        self,
        topology,
        system,
        l,
        init_pos,
        reasim,
        prosim,
        conformer_search,
    ):
        result = {}
        # else:
        opm_dyn = OpenMMDynamics()
        opm_dyn.ostream.mute()
        # platform settings for small molecule
        opm_dyn.openmm_platform = "CPU"
        # opm_dyn.create_system_from_molecule(mol, ff_gen)
        if self.save_intermediates:
            pdb_name = str(Path(self.folder_name) / f'conf_top_{l}.pdb')
        else:
            pdb_name = f'topology_{getrandbits(32):08x}.pdb'

        pdb_not_used = mmapp.PDBFile.writeFile(
            topology,
            init_pos * mmunit.angstrom,
            pdb_name,
        )
        opm_dyn.pdb = mmapp.PDBFile(pdb_name)
        opm_dyn.system = system

        if conformer_search:
            snapshots = self.conformer_snapshots
        else:
            snapshots = 1

        conformers_dict = opm_dyn.conformational_sampling(
            ensemble='NVT',
            nsteps=self.mm_steps * snapshots,
            snapshots=snapshots,
            temperature=self.mm_temperature,
        )
        if not self.save_intermediates:
            Path(pdb_name).unlink()

        result = []
        for e_int, temp_mol in zip(conformers_dict['energies'],
                                   conformers_dict['molecules']):
            pos = temp_mol.get_coordinates_in_angstrom()
            v, e1, e2 = self._recalc_mm_energy(pos, l, reasim, prosim)
            avg_x = np.mean(pos[:, 0])
            avg_y = np.mean(pos[:, 1])
            avg_z = np.mean(pos[:, 2])
            pos -= [avg_x, avg_y, avg_z]
            xyz = temp_mol.get_xyz_string()
            temp_result = {
                'v': v,
                'e1': e1,
                'e2': e2,
                'e_int': e_int,
                'pos': pos,
                'xyz': xyz
            }
            result.append(temp_result)

        return result

    def _recalc_mm_energy(self, pos, l, rea_sim, pro_sim):
        # unit conversion from angstrom to nm
        pos_nm = pos * 0.1
        rea_sim.context.setPositions(pos_nm)
        pro_sim.context.setPositions(pos_nm)
        e1 = rea_sim.context.getState(
            getEnergy=True).getPotentialEnergy().value_in_unit(
                mmunit.kilojoules_per_mole)
        e2 = pro_sim.context.getState(
            getEnergy=True).getPotentialEnergy().value_in_unit(
                mmunit.kilojoules_per_mole)
        em = e1 * (1 - l) + e2 * l

        return em, e1, e2

    def scan_qm(self, results=None):
        if results is None:
            results = self.results
        assert_msg_critical(
            'scan' in results.keys(),
            f'Could not find "scan" in results. Total keys: {results.keys()}',
        )

        self._print_qm_header()
        ref = None
        max_qm_energy = None
        min_qm_conf_index = 0
        try:
            for l in results['scan'].keys():

                min_qm_conf_E = None
                # min_conf_index = 0
                scan = sorted(results['scan'][l], key=lambda x: x['v'])
                results['scan'][l] = scan
                # Pick out lowest 5 conformers from scan

                scan_range = min(self.max_qm_conformers, len(scan))
                for i, conformer in enumerate(scan[:scan_range]):
                    qm_E = self._get_qm_energy(conformer['xyz'])

                    results['scan'][l][i]['qm_energy'] = qm_E
                    if math.isnan(qm_E):
                        continue
                    if min_qm_conf_E is None or qm_E < min_qm_conf_E:
                        min_qm_conf_E = qm_E
                        # min_conf_index = i

                    if ref is None:
                        ref = qm_E
                    dif = qm_E - ref
                    mm_E = results['scan'][l][i]['v']

                    self._print_qm_iter(l, qm_E, mm_E, dif, i)

                if max_qm_energy is None or min_qm_conf_E > max_qm_energy:
                    max_qm_energy = min_qm_conf_E
                    max_qm_lambda = l
                    # min_qm_conf_index = min_conf_index
                    # max_qm_xyz = scan[min_conf_index]['xyz']
                qm_scan = results['scan'][l][:scan_range]
                rest = results['scan'][l][scan_range:]
                sorted_qm_scan = sorted(qm_scan, key=lambda x: x['qm_energy'])
                results['scan'][l] = sorted_qm_scan + rest

            self.ostream.print_blank()
            min_qm_conf_index, _ = min(
                enumerate(results['scan'][max_qm_lambda]),
                key=lambda x: x[1].get('qm_energy', math.inf))
            max_qm_xyz = results['scan'][max_qm_lambda][min_qm_conf_index][
                'xyz']
            results.update({
                'max_qm_xyz': max_qm_xyz,
                'max_qm_lambda': max_qm_lambda,
                'min_qm_conformer_index': min_qm_conf_index,
            })
            self.ostream.print_info(
                f"Found highest QM E: {max_qm_energy:.3f} at Lambda: {max_qm_lambda} and conformer index: {min_qm_conf_index}."
            )
            self.ostream.flush()
            self.results.update(results)

            self.molecule = Molecule.read_xyz_string(max_qm_xyz)
            self.molecule.set_multiplicity(self.mol_multiplicity)
            self.molecule.set_charge(self.mol_charge)
        except Exception:
            err_str = "The QM scan crashed. Saving results in self.results and raising exception"
            self.ostream.print_warning(err_str)
            self.ostream.flush()
            self.results.update(results)
            raise

        if self.save_results_file:
            self.save_results(self.results_file, self.results)

        return self.results

    def _get_qm_energy(self, xyz):
        self.molecule = Molecule.read_xyz_string(xyz)
        self.molecule.set_multiplicity(self.mol_multiplicity)
        self.molecule.set_charge(self.mol_charge)
        if self.scf_drv is None:
            scf_drv = ScfRestrictedDriver()
            scf_drv.xcfun = self.qm_xcfun
            self.scf_drv = scf_drv

            # When MM implicit solvation is active, mirror it on the QM driver via
            # SMD.  Duck-type on solvation_model: SCF drivers expose it, XTB does
            # not.  If the attribute is absent we warn and skip rather than error.
            if self.implicit_solvent_model is not None:
                if hasattr(self.scf_drv, 'solvation_model'):
                    self.scf_drv.solvation_model = 'smd'
                    self.scf_drv.smd_solvent = self.smd_solvent
                else:
                    self.ostream.print_warning(
                        'QM driver does not support SMD solvation. '
                        'Implicit solvation will not be applied to the QM scan. '
                        'Use an SCF driver instead of XTB to enable SMD.')
                    self.ostream.flush()
        if self.implicit_solvent_model is not None:
            if hasattr(
                    self.scf_drv,
                    'solvation_model') and self.scf_drv.solvation_model is None:
                self.ostream.print_warning(
                    'Implicit solvation turned on, but explicitly provided SCF '
                    'driver has no solvation model activated. Continuing without QM solvation. '
                    'Provide an SCF driver with a solvation model activated to enable QM solvation.'
                )

        if self.mute_scf:
            self.scf_drv.ostream.mute()
        basis = MolecularBasis.read(self.molecule, self.qm_basis)
        scf_results = self.scf_drv.compute(self.molecule, basis)

        if not self.scf_drv.is_converged:
            self.ostream.print_warning(
                "SCF did not converge, increasing convergence threshold to 1.0e-4 and maximum itterations to 200."
            )
            self.ostream.flush()
            self.scf_drv.conv_thresh = 1.0e-4
            self.scf_drv.max_iter = 200
            scf_results = self.scf_drv.compute(self.molecule, basis)

            if not self.scf_drv.is_converged:
                self.ostream.print_warning(
                    "SCF still did not converge. Returning NaN for current calculation"
                )
                self.ostream.flush()
                return math.nan
        return scf_results['scf_energy'] * hartree_in_kjpermol()

    def show_results(self, ts_results=None, filename=None, **mol_show_kwargs):
        """Show the results of the transition state guesser.
        This function uses ipywidgets to create an interactive plot of the MM and SCF energies as a function of lambda.

        Args:
            ts_results (dict, optional): The results of the transition state guesser. If none (default), uses the current results.
            atom_indices (bool, optional): If true, shows the atom-indices of the molecule. Defaults to False.

        Raises:
            ImportError: If ipywidgets is not installed, an ImportError is raised.
        """

        try:
            import ipywidgets
            from IPython.display import display as ipy_display
        except ImportError:
            raise ImportError('ipywidgets is required for this functionality.')

        ostream = OutputStream(sys.stdout)

        if ts_results is None:
            if filename is not None:
                ostream.print_info(f"Loading results from {filename}")
                ts_results = TransitionStateGuesser.load_results(
                    filename,
                    ostream,
                )
            else:
                if self.results is None or 'scan' not in self.results.keys():
                    raise ValueError(
                        "No results provided. Provide either ts_results or filename."
                    )
                ts_results = self.results

        scan = ts_results['scan']

        if TransitionStateGuesser._has_qm_results(scan,
                                                  ts_results['lambda_vec']):
            final_lambda = round(float(ts_results.get('max_qm_lambda', 0)), 3)
        else:
            final_lambda = round(float(ts_results['max_mm_lambda']), 3)

        forming_bonds = set(ts_results.get('forming_bonds', set()))
        breaking_bonds = set(ts_results.get('breaking_bonds', set()))
        bonds = set(ts_results.get('static_bonds', set()))

        def _best_conformer_index(step):
            best_index = 0
            min_energy = None
            if TransitionStateGuesser._has_qm_results(scan,
                                                      ts_results['lambda_vec']):
                for i, conf in enumerate(scan[step]):
                    qm_E = conf.get('qm_energy', None)
                    if min_energy is None or (qm_E is not None
                                              and qm_E < min_energy):
                        min_energy = qm_E
                        best_index = i
            else:
                for i, conf in enumerate(scan[step]):
                    if min_energy is None or conf['v'] < min_energy:
                        min_energy = conf['v']
                        best_index = i
            return best_index

        initial_best = _best_conformer_index(final_lambda)

        rounded_lambda_vec = [
            round(float(l), 3) for l in ts_results['lambda_vec']
        ]
        step_selector = ipywidgets.SelectionSlider(
            options=rounded_lambda_vec,
            description='Lambda',
            value=final_lambda,
        )
        conformer_selector = ipywidgets.Dropdown(
            options=list(range(1,
                               len(scan[final_lambda]) + 1)),
            description='Conf. ID',
            value=initial_best + 1,
        )

        def _render_interact(step, conformer_id):
            # Clamp conformer_id in case it exceeds available conformers for
            # the current lambda.
            n_conf = len(scan[step])
            TransitionStateGuesser._show_conformer_iteration(
                step,
                min(conformer_id, n_conf),
                rounded_lambda_vec,
                scan,
                bonds=bonds,
                forming_bonds=forming_bonds,
                breaking_bonds=breaking_bonds,
                **mol_show_kwargs,
            )

        # Hide the output until the delayed update has populated it.  Keeping
        # the output in the layout avoids a visible resize flicker when it is
        # revealed.
        out_widget = ipywidgets.Output()
        out_widget.layout.visibility = 'hidden'
        window = ipywidgets.VBox([
            step_selector,
            conformer_selector,
            out_widget,
        ])

        lambda_change_in_progress = [False]

        def _render_current_selection():
            with out_widget:
                out_widget.clear_output(wait=True)
                _render_interact(step_selector.value, conformer_selector.value)

        def _on_lambda_change(change):
            lambda_change_in_progress[0] = True
            try:
                step = change['new']
                options = list(range(1, len(scan[step]) + 1))
                best = _best_conformer_index(step)
                conformer_selector.options = options
                conformer_selector.value = best + 1
            finally:
                lambda_change_in_progress[0] = False
            _render_current_selection()

        def _on_conformer_change(change):
            if not lambda_change_in_progress[0]:
                _render_current_selection()

        step_selector.observe(_on_lambda_change, names='value')
        conformer_selector.observe(_on_conformer_change, names='value')

        ipy_display(window)

        def _finish_render():
            out_widget.layout.visibility = 'visible'

        def _do_render():
            try:
                # Render once the widget output has reached the frontend, so
                # py3Dmol's script runs against an existing output area.
                _render_current_selection()
            except Exception:
                pass
            finally:
                _schedule_callback(_finish_render, delay=0.5)

        def _schedule_callback(callback, delay=0.0):
            try:
                loop = _ip.kernel.io_loop
                loop.call_later(delay, callback)
                return
            except Exception:
                pass

            if delay > 0.0:
                time.sleep(delay)
            callback()

        # Schedule _do_render from a post_execute hook so the delay begins when
        # the cell finishes and the widget output area has been sent to the
        # frontend.
        try:
            from IPython import get_ipython as _get_ipython
            _ip = _get_ipython()
        except ImportError:
            _ip = None

        if _ip is not None:

            def _on_post_execute():
                _ip.events.unregister('post_execute', _on_post_execute)
                _schedule_callback(_do_render, delay=0.3)

            _ip.events.register('post_execute', _on_post_execute)
        else:
            _schedule_callback(_do_render, delay=0.3)

    @staticmethod
    def _show_conformer_iteration(
        step,
        conformer_id,
        lambda_vec,
        scan,
        bonds=None,
        forming_bonds=None,
        breaking_bonds=None,
        **mol_show_kwargs,
    ):
        """
        Show the geometry at a specific iteration.
        """

        try:
            import matplotlib.pyplot as plt
        except ImportError:
            raise ImportError('matplotlib is required for this functionality.')

        mm_energies, _, _, _ = TransitionStateGuesser._get_best_mm_E_from_scan_dict(
            scan)
        mm_min = np.min(mm_energies)
        rel_mm_energies = np.asarray(mm_energies) - mm_min

        lam_index = lambda_vec.index(step)
        xyz_i = scan[step][conformer_id - 1]['xyz']
        total_steps = len(rel_mm_energies) - 1
        x = np.linspace(0, lambda_vec[-1], 100)
        y = np.interp(x, lambda_vec, rel_mm_energies)
        fig, ax1 = plt.subplots(figsize=(6.5, 4))

        if TransitionStateGuesser._has_qm_results(scan, lambda_vec):
            qm_energies, _ = TransitionStateGuesser._get_best_qm_E_from_scan_dict(
                scan)
            qm_min = np.min(qm_energies)
            rel_qm_energies = np.asarray(qm_energies) - qm_min
            print("  {:>9} {:>18} {:>19}  ".format("conformer",
                                                   "MM energy [kJ/mol]",
                                                   "QM energy [kJ/mol]"))

            for i, conf in enumerate(scan[step]):
                conf_str = f"{i+1}"
                mm_e = f"{conf['v'] - mm_min:.3f}"
                qm_e_raw = conf.get('qm_energy', None)
                if qm_e_raw is None:
                    qm_e = ""
                else:
                    qm_e = f"{conf['qm_energy'] - qm_min:.3f}"

                if i + 1 == conformer_id and len(scan[step]) > 1:
                    conf_str_formatted = f"{'-' * (9 - (1 + len(str(conf_str))))} {conf_str}"
                    mm_e_formatted = f"{'-' * (18 - (1 + len(str(mm_e))))} {mm_e}"
                    qm_e_formatted = f"{'-' * (19 - (1 + len(str(qm_e))))} {qm_e}"
                    print_str = "{:>1} {:>9} {:>18} {:>19} {:>1}".format(
                        ">",
                        conf_str_formatted,
                        mm_e_formatted,
                        qm_e_formatted,
                        "<",
                    )
                else:
                    print_str = "  {:>9} {:>18} {:>19}  ".format(
                        conf_str,
                        mm_e,
                        qm_e,
                    )
                print(print_str)

        else:
            print("  {:>9} {:>19}  ".format("conformer", "MM energy [kJ/mol]"))
            rel_qm_energies = None
            for i, conf in enumerate(scan[step]):
                conf_str = f"{i+1}"
                mm_e = f"{conf['v'] - mm_min:.3f}"
                if i + 1 == conformer_id and len(scan[step]) > 1:
                    conf_str_formatted = f"{'-' * (9 - (1 + len(str(conf_str))))} {conf_str}"
                    mm_e_formatted = f"{'-' * (18 - (1 + len(str(mm_e))))} {mm_e}"
                    print_str = "{:>1} {:>9} {:>19} {:>1}".format(
                        ">",
                        conf_str_formatted,
                        mm_e_formatted,
                        "<",
                    )
                else:
                    print_str = "  {:>9} {:>19}  ".format(conf_str, mm_e)
                print(print_str)

        # Collect all conformer MM energies for the stripe markers.
        conf_x_mm, conf_y_mm = [], []
        for lv in lambda_vec:
            for conf in scan[lv]:
                conf_x_mm.append(lv)
                conf_y_mm.append(conf['v'] - mm_min)

        # Collect all conformer QM energies for the stripe markers.
        conf_x_qm, conf_y_qm = [], []
        if rel_qm_energies is not None:
            for lv in lambda_vec:
                for conf in scan[lv]:
                    qm_e = conf.get('qm_energy', None)
                    if qm_e is not None and not math.isnan(qm_e):
                        conf_x_qm.append(lv)
                        conf_y_qm.append(qm_e - qm_min)

        ax1.plot(
            x,
            y,
            color='black',
            alpha=0.9,
            linewidth=2.5,
            linestyle='-',
            zorder=0,
            label='MM energy',
        )
        ax1.scatter(
            conf_x_mm,
            conf_y_mm,
            marker='_',
            color='darkcyan',
            alpha=0.4,
            s=30 / math.log(min(2, total_steps), 10),
            linewidths=2.0,
            zorder=0.5,
        )
        ax1.scatter(
            lambda_vec,
            rel_mm_energies,
            color='black',
            alpha=0.7,
            s=50 / math.log(min(2, total_steps), 10),
            facecolors="none",
            edgecolor="darkcyan",
            zorder=1,
        )
        selected_mm_e = scan[step][conformer_id - 1]['v'] - mm_min
        ax1.scatter(
            lambda_vec[lam_index],
            selected_mm_e,
            marker='o',
            color='darkcyan',
            alpha=1.0,
            s=50 / math.log(min(2, total_steps), 10),
            zorder=2,
        )
        ax1.set_xlabel(r'$\lambda$')
        ax1.set_ylabel('Relative MM energy [kJ/mol]')

        if rel_qm_energies is not None:
            ax1.plot(
                x,
                np.interp(x, lambda_vec, rel_qm_energies),
                color='darkred',
                alpha=0.9,
                linewidth=2.5,
                linestyle='--',
                zorder=0,
                label='QM energy',
            )
            ax1.scatter(
                conf_x_qm,
                conf_y_qm,
                marker='_',
                color='darkorange',
                alpha=0.4,
                s=30 / math.log(min(2, total_steps), 10),
                linewidths=2.0,
                zorder=0.5,
            )
            ax1.scatter(
                lambda_vec,
                rel_qm_energies,
                alpha=0.7,
                s=50 / math.log(min(2, total_steps), 10),
                facecolors="none",
                edgecolor="darkorange",
                zorder=1,
            )
            selected_qm_e_raw = scan[step][conformer_id - 1].get(
                'qm_energy', None)
            if selected_qm_e_raw is not None and not math.isnan(
                    selected_qm_e_raw):
                selected_qm_e = selected_qm_e_raw - qm_min
            else:
                selected_qm_e = rel_qm_energies[lam_index]
            ax1.scatter(
                lambda_vec[lam_index],
                selected_qm_e,
                marker='o',
                color='darkorange',
                alpha=1.0,
                s=120 / math.log(min(2, total_steps), 10),
                zorder=2,
            )
            ax1.set_ylabel('Relative QM energy [kJ/mol]')

        fig.legend(loc='upper right',
                   bbox_to_anchor=(1, 1),
                   bbox_transform=ax1.transAxes)
        ax1.set_title("Transition state finder")
        ax1.set_xticks(lambda_vec[::2])
        fig.tight_layout()
        plt.show()
        mol = Molecule.read_xyz_string(xyz_i)
        offset = 0.07
        add = 0.1
        breaking_width = offset + add * (1 - step)
        forming_width = offset + add * (step)

        if bonds is not None and (forming_bonds is not None
                                  or breaking_bonds is not None):
            mol.show(
                bonds=bonds,
                forming_bonds=list(forming_bonds),
                breaking_bonds=list(breaking_bonds),
                forming_width=forming_width,
                breaking_width=breaking_width,
                width=640,
                height=360,
                **mol_show_kwargs,
            )
        else:
            mol.show(width=640, height=360, **mol_show_kwargs)

    def _mm_to_xyz_str(self, positions, molecule=None):
        if molecule is None:
            molecule = self.molecule
        new_mol = TransitionStateGuesser._set_molecule_positions(
            molecule, positions)
        return new_mol.get_xyz_string()

    @staticmethod
    def _get_best_mm_E_from_scan_dict(scan):
        V = []
        E1 = []
        E2 = []
        conf_indices = []
        for conf_scan in scan.values():
            lowest_v = None
            index = 0
            for i, conf in enumerate(conf_scan):
                if lowest_v is None or conf['v'] < lowest_v:
                    lowest_v = conf['v']
                    lowest_e1 = conf['e1']
                    lowest_e2 = conf['e2']
                    index = i
            V.append(lowest_v)
            E1.append(lowest_e1)
            E2.append(lowest_e2)
            conf_indices.append(index)

        return V, E1, E2, conf_indices

    @staticmethod
    def _get_best_qm_E_from_scan_dict(scan):
        qm_energies = []
        conf_indices = []
        for conf_scan in scan.values():
            lowest_qm_e = None
            index = 0
            for i, conf in enumerate(conf_scan):
                qm_e = conf.get('qm_energy', None)
                if qm_e is not None and not math.isnan(qm_e):
                    if lowest_qm_e is None or qm_e < lowest_qm_e:
                        lowest_qm_e = qm_e
                        index = i
            qm_energies.append(lowest_qm_e)
            conf_indices.append(index)

        return qm_energies, conf_indices

    @staticmethod
    def _has_qm_results(scan, lambda_vec):
        """Return True if qm_energy values are present in the scan dict."""
        return scan[lambda_vec[0]][0].get('qm_energy', None) is not None

    @staticmethod
    def _set_molecule_positions(molecule, positions):
        positions_au = positions / bohr_in_angstrom()
        assert molecule.number_of_atoms() == len(positions_au)
        for i in range(molecule.number_of_atoms()):
            molecule.set_atom_coordinates(i, positions_au[i])
        return molecule

    def save_results(self, fname, results):
        required_keys = {
            'max_mm_xyz', 'max_mm_lambda', 'min_mm_conformer_index'
        }
        missing = required_keys - results.keys()
        if missing:
            raise ValueError("Cannot save results: the MM scan is incomplete. "
                             f"Missing keys: {sorted(missing)}. "
                             "Call scan_mm() before saving.")
        self.ostream.print_info(f"Saving results to {fname}")
        self.ostream.flush()

        def _bonds_to_array(bonds):
            # Always store as shape (N, 2) — even for empty sets — so that
            # load_results can iterate rows uniformly regardless of set size.
            if len(bonds) == 0:
                return np.empty((0, 2), dtype='i')
            return np.array(list(bonds), dtype='i')

        with h5py.File(fname, 'w') as hf:
            # breaking forming static bonds
            hf.create_dataset('breaking_bonds',
                              data=_bonds_to_array(results['breaking_bonds']))
            hf.create_dataset('forming_bonds',
                              data=_bonds_to_array(results['forming_bonds']))
            hf.create_dataset('static_bonds',
                              data=_bonds_to_array(results['static_bonds']))

            # lambda vec
            lambda_vec = results.get('lambda_vec', None)
            hf.create_dataset('lambda_vec', data=lambda_vec, dtype='f')

            # reactant and product
            rea_ff = results['reactant'].get_forcefield_as_json(
                results['reactant'])
            pro_ff = results['product'].get_forcefield_as_json(
                results['product'])
            rea_xyz = results['reactant'].molecule.get_xyz_string()
            pro_xyz = results['product'].molecule.get_xyz_string()
            hf.create_dataset('reactant_ff', data=rea_ff)
            hf.create_dataset('product_ff', data=pro_ff)
            hf.create_dataset('reactant_xyz', data=rea_xyz)
            hf.create_dataset('product_xyz', data=pro_xyz)

            # max_mm_xyz max_mm_lambda min_mm_conformer_index
            max_mm_xyz = results['max_mm_xyz']
            max_mm_lambda = results['max_mm_lambda']
            min_mm_conformer_index = results['min_mm_conformer_index']
            hf.create_dataset('max_mm_xyz', data=[max_mm_xyz])
            hf.create_dataset('max_mm_lambda', data=max_mm_lambda, dtype='f')
            hf.create_dataset('min_mm_conformer_index',
                              data=min_mm_conformer_index,
                              dtype='i')

            max_qm_xyz = results.get('max_qm_xyz', None)
            max_qm_lambda = results.get('max_qm_lambda', None)
            min_qm_conformer_index = results.get('min_qm_conformer_index', None)
            if max_qm_xyz is not None:
                hf.create_dataset('max_qm_xyz', data=[max_qm_xyz])
                hf.create_dataset('max_qm_lambda',
                                  data=max_qm_lambda,
                                  dtype='f')
                hf.create_dataset('min_qm_conformer_index',
                                  data=min_qm_conformer_index,
                                  dtype='i')

            scan_grp = hf.create_group('scan')
            for l, conf_scan in results['scan'].items():
                l_grp = scan_grp.create_group(f'{round(l,3)}')
                for i, conf in enumerate(conf_scan):
                    conf_grp = l_grp.create_group(str(i))
                    # v e1 e2 e_int xyz qm
                    conf_grp.create_dataset('v', data=conf['v'], dtype='f')
                    conf_grp.create_dataset('e1', data=conf['e1'], dtype='f')
                    conf_grp.create_dataset('e2', data=conf['e2'], dtype='f')
                    conf_grp.create_dataset('e_int',
                                            data=conf['e_int'],
                                            dtype='f')
                    conf_grp.create_dataset('xyz', data=[conf['xyz']])
                    qm_e = conf.get('qm_energy', None)
                    if qm_e is not None:
                        conf_grp.create_dataset('qm_energy',
                                                data=qm_e,
                                                dtype='f')

    @staticmethod
    def load_results(fname, ostream=None):
        if ostream is None:
            ostream = OutputStream(sys.stdout)
        ostream.print_info(f"Loading results from {fname}")
        ostream.flush()

        results = {}
        with h5py.File(fname, 'r') as hf:
            # Bonds
            results['breaking_bonds'] = {
                tuple(map(int, b))
                for b in hf['breaking_bonds'][()]
            }
            results['forming_bonds'] = {
                tuple(map(int, b))
                for b in hf['forming_bonds'][()]
            }
            results['static_bonds'] = {
                tuple(map(int, b))
                for b in hf['static_bonds'][()]
            }

            # Lambda vector
            results['lambda_vec'] = [
                round(float(l), 3) for l in hf['lambda_vec'][()]
            ]

            # Reactant and product forcefields and xyz
            reactant_ff = MMForceFieldGenerator.load_forcefield_from_json_string(
                hf['reactant_ff'][()])
            product_ff = MMForceFieldGenerator.load_forcefield_from_json_string(
                hf['product_ff'][()])
            reactant_xyz = hf['reactant_xyz'][()].decode('utf-8')
            product_xyz = hf['product_xyz'][()].decode('utf-8')
            reactant_mol = Molecule.read_xyz_string(reactant_xyz)
            product_mol = Molecule.read_xyz_string(product_xyz)
            reactant_ff.molecule = reactant_mol
            product_ff.molecule = product_mol
            results['reactant'] = reactant_ff
            results['product'] = product_ff

            # Max/min xyzs and lambdas
            results['max_mm_xyz'] = hf['max_mm_xyz'][()][0].decode('utf-8')
            results['max_mm_lambda'] = hf['max_mm_lambda'][()]
            results['min_mm_conformer_index'] = hf['min_mm_conformer_index'][()]

            # Optional QM results
            if 'max_qm_xyz' in hf:
                results['max_qm_xyz'] = hf['max_qm_xyz'][()][0].decode('utf-8')
                results['max_qm_lambda'] = hf['max_qm_lambda'][()]
                results['min_qm_conformer_index'] = hf[
                    'min_qm_conformer_index'][()]

            # Scan results
            results['scan'] = {}
            scan_grp = hf['scan']
            for l in scan_grp:
                l_grp = scan_grp[l]
                conf_scan = []
                for i in sorted(l_grp, key=lambda x: int(x)):
                    conf_grp = l_grp[i]
                    conf = {
                        'v': conf_grp['v'][()],
                        'e1': conf_grp['e1'][()],
                        'e2': conf_grp['e2'][()],
                        'e_int': conf_grp['e_int'][()],
                        'xyz': conf_grp['xyz'][()][0].decode('utf-8')
                    }
                    if 'qm_energy' in conf_grp:
                        conf['qm_energy'] = conf_grp['qm_energy'][()]
                    conf_scan.append(conf)
                results['scan'][round(float(l), 3)] = conf_scan

        return results

    @staticmethod
    def _param(label, value, lw=24, vw=20):
        """Format one parameter line with fixed label and value widths.

        Because print_header centers text, all lines must be the same total
        length to appear left-aligned relative to each other.
        Total width = lw + len(' : ') + vw = 47 chars (default).
        """
        return f"{label:<{lw}} : {str(value):>{vw}}"

    def _print_mm_header(self, conformer_search=False):
        self.ostream.print_blank()
        if conformer_search:
            self.ostream.print_header(
                "Starting MM Scan  (with conformer search)")
        else:
            self.ostream.print_header("Starting MM Scan")
        self.ostream.print_blank()

        conf_snapshots = self.conformer_snapshots if conformer_search else 1
        self.ostream.print_header(self._param("MD steps", self.mm_steps))
        self.ostream.print_header(self._param("conf. snapshots",
                                              conf_snapshots))
        self.ostream.print_header(
            self._param("MD temperature", f"{self.mm_temperature} K"))
        self.ostream.print_header(
            self._param("MD step size", f"{self.mm_step_size} ps"))
        self.ostream.print_header(self._param("folder name", self.folder_name))
        self.ostream.print_header(
            self._param("saving MD traj", self.save_mm_traj))
        if self.implicit_solvent_model is not None:
            self.ostream.print_header(
                self._param("implicit solvent", self.implicit_solvent_model))
            self.ostream.print_header(
                self._param("solvent dielectric",
                            f"{self.solvent_dielectric:.2f}"))
            self.ostream.print_header(
                self._param("solute dielectric",
                            f"{self.solute_dielectric:.2f}"))

        if self._conformer_active_torsion is not None:
            one_based = tuple(a + 1 for a in self._conformer_active_torsion)
            self.ostream.print_blank()
            self.ostream.print_header(self._param("conformational TS", "yes"))
            self.ostream.print_header(
                self._param("active torsion (1-idx)", str(one_based)))
            self.ostream.print_header(
                self._param("phi reactant",
                            f"{self._conformer_phi_reactant:.1f} deg"))
            self.ostream.print_header(
                self._param("phi product",
                            f"{self._conformer_phi_product:.1f} deg"))
            self.ostream.print_header(
                self._param("restraint k", f"{self.conformer_k:.1f} kJ/mol"))

        self.ostream.print_blank()
        valstr = '{} | {} | {} | {} | {}'.format(
            'Lambda',
            'E1 (kj/mol)',
            'E2 (kj/mol)',
            'V (kj/mol)',
            'n_conf',
        )
        self.ostream.print_header(valstr)
        self.ostream.print_header(60 * '-')
        self.ostream.flush()

    def _print_mm_iter(self, l, e1, e2, v, n_conf=None):

        valstr = "{:6.2f}   {:10.2f}   {:10.2f}   {:9.2f}   {:6}".format(
            l, e1, e2, v, n_conf)
        self.ostream.print_header(valstr)
        self.ostream.flush()

    def _print_qm_header(self):
        if self.mute_scf:
            self.ostream.print_info("Disable mute_scf to see detailed output.")
        self.ostream.print_blank()
        self.ostream.print_header("Starting QM Scan")
        self.ostream.print_blank()

        self.ostream.print_header(self._param("basis set", self.qm_basis))
        self.ostream.print_header(
            self._param("DFT xc functional", self.qm_xcfun))
        self.ostream.print_header(
            self._param("max conformers", self.max_qm_conformers))
        if self.implicit_solvent_model is not None:
            self.ostream.print_header(self._param("solvation", "SMD"))
            self.ostream.print_header(
                self._param("SMD solvent", self.smd_solvent))

        self.ostream.print_blank()
        self.ostream.flush()
        valsltr = '{} | {} | {} | {}'.format(
            'Lambda',
            'Conf. i',
            'Rel. E (kJ/mol)',
            'MM V (kJ/mol)',
        )
        self.ostream.print_header(valsltr)
        self.ostream.print_header(60 * '-')
        self.ostream.flush()

    def _print_qm_iter(self, l, qm_E, mm_E, dif, conf_index):

        valstr = "{:6.2f}   {:7d}   {:15.3f}   {:14.3f}".format(
            l, conf_index + 1, dif, mm_E)
        self.ostream.print_header(valstr)
        self.ostream.flush()

    def _detect_active_dihedral(self):
        """Identify the dihedral that differs most between reactant and product.

        Called when no bonds are forming or breaking (conformational TS).
        Iterates over every dihedral in the reactant force field, one
        representative per unique central bond, and picks the one with the
        largest absolute angle difference between the two geometries.

        Returns a tuple (active_torsion, phi_reactant, phi_product) where
        active_torsion is a 0-indexed (i, j, k, l) tuple and the angles are
        in degrees.

        Raises:
            ValueError: if no dihedral with |Δφ| > 20° is found.
        """
        threshold = 20.0  # degrees

        max_delta = 0.0
        active_torsion = None
        best_phi_rea = None
        best_phi_pro = None

        # Deduplicate by central bond so each bond is measured exactly once.
        seen_central = set()
        for key in self.reactant.dihedrals:
            central = tuple(sorted([key[1], key[2]]))
            if central in seen_central:
                continue
            seen_central.add(central)

            # get_dihedral_in_degrees expects 1-based indices
            one_based = [idx + 1 for idx in key]
            phi_rea = self.reactant.molecule.get_dihedral_in_degrees(one_based)
            phi_pro = self.product.molecule.get_dihedral_in_degrees(one_based)

            # Shortest-path difference, wrapped to (-180, 180]
            delta = phi_pro - phi_rea
            if delta > 180.0:
                delta -= 360.0
            elif delta <= -180.0:
                delta += 360.0

            if abs(delta) > max_delta:
                max_delta = abs(delta)
                active_torsion = key
                best_phi_rea = phi_rea
                best_phi_pro = phi_pro

        if active_torsion is None or max_delta < threshold:
            raise ValueError(
                f"No dihedral with |Δφ| > {threshold:.0f}° found "
                f"between reactant and product geometries "
                f"(largest difference: {max_delta:.1f}°). "
                "Verify that the two input geometries differ by a dihedral rotation."
            )

        assert best_phi_rea is not None
        assert best_phi_pro is not None
        self.ostream.print_info(f"Active torsion (1-indexed): "
                                f"{tuple(a + 1 for a in active_torsion)}, "
                                f"phi_reactant = {best_phi_rea:.1f}°, "
                                f"phi_product = {best_phi_pro:.1f}°, "
                                f"|Δφ| = {max_delta:.1f}°.")
        self.ostream.flush()
        return active_torsion, best_phi_rea, best_phi_pro

    # todo add option for reading geometry (bond distances, angles, etc.) from transition state instead of averaging them
    # todo add option for recalculating charges from ts_mol
    def get_ts_ffgen(self,
                     reaffgen=None,
                     proffgen=None,
                     l=0.5,
                     ts_mol=None,
                     recalculate=True):
        if reaffgen is None:
            reaffgen = self.reactant
        if proffgen is None:
            proffgen = self.product
        if ts_mol is None:
            ts_mol = self.molecule

        if recalculate:
            assert ts_mol is not None, "Please provide a molecule, or turn of recalculation"
        ts_ffgen = MMForceFieldGenerator()

        ts_ffgen.bonds = self._average_params(
            reaffgen.bonds,
            proffgen.bonds,
            l,
            ts_mol,
        )
        ts_ffgen.angles = self._average_params(
            reaffgen.angles,
            proffgen.angles,
            l,
            ts_mol,
        )
        ts_ffgen.dihedrals = self._mix_dihedrals(
            reaffgen.dihedrals,
            proffgen.dihedrals,
            l,
            ts_mol,
        )
        ts_ffgen.impropers = self._mix_dihedrals(
            reaffgen.impropers,
            proffgen.impropers,
            l,
            ts_mol,
        )
        ts_ffgen.atoms = self._merge_atoms(
            reaffgen.atoms,
            proffgen.atoms,
            l,
            ts_mol,
        )
        return ts_ffgen

    @staticmethod
    def _average_params(reaparams, proparams, l, mol):
        ts_params = {}
        for id in reaparams.keys() | proparams.keys():
            param = {}
            if id in reaparams.keys() and id in proparams.keys():
                reaparam = reaparams[id]
                proparam = proparams[id]
                # todo change this to measurement from molecule
                eq = (1 -
                      l) * reaparam['equilibrium'] + l * proparam['equilibrium']
                fc = (1 - l) * reaparam['force_constant'] + l * proparam[
                    'force_constant']
                comment = f"averaged {reaparam['comment']} {proparam['comment']}"
                param = {
                    'force_constant': fc,
                    'equilibrium': eq,
                    'comment': comment,
                }
            elif id in reaparams.keys():
                param = reaparams[id]
                param['force_constant'] *= (1 - l)
            elif id in proparams.keys():
                param = proparams[id]
                param['force_constant'] *= l
            else:
                continue
            ts_params.update({id: copy.copy(param)})
        if not ts_params:
            return None
        return ts_params

    @staticmethod
    def _mix_dihedrals(rea_dihedrals, pro_dihedrals, l, mol):
        ts_params = {}
        for dict, scaling in zip([rea_dihedrals, pro_dihedrals], [1 - l, l]):
            for id, param in dict.items():
                # todo reassign value of phase?
                new_param = copy.copy(param)
                if new_param.get('multiple'):
                    new_param['barrier'] = [
                        bar * scaling for bar in new_param['barrier']
                    ]
                else:
                    new_param['barrier'] *= scaling
                ts_params.update({id: new_param})
        return ts_params

    @staticmethod
    def _merge_atoms(rea_atoms, pro_atoms, l, mol):
        atoms = {}
        for i, (rea_atom, pro_atom) in enumerate(
                zip(rea_atoms.values(), pro_atoms.values())):
            assert rea_atom["mass"] == pro_atom[
                "mass"], "Atoms in reactont and product are not the same"

            name = ''.join(x for x in rea_atom['name'] if x.isalpha())

            if rea_atom['type'] == pro_atom['type']:
                typ = rea_atom['type']
            else:
                typ = name.lower() + "rx"

            if rea_atom.get('comment') == pro_atom.get(
                    'comment') and rea_atom.get('comment') is not None:
                comment = rea_atom['comment']
            else:
                comment = "Merged TS atom"

            charge = (1 - l) * rea_atom['charge'] + l * pro_atom['charge']
            eps = (1 - l) * rea_atom['epsilon'] + l * pro_atom['epsilon']
            sigma = (1 - l) * rea_atom['sigma'] + l * pro_atom['sigma']

            atom = {
                'name': f"{name}{i}",
                'sigma': sigma,
                'epsilon': eps,
                'charge': charge,
                'mass': rea_atom["mass"],
                'type': typ,
                'comment': comment
            }

            if rea_atom.get('equivalent_atom') == pro_atom.get(
                    'equivalent_atom'):
                equi_atom = rea_atom['equivalent_atom']
                atom.update({'equivalent_atom': equi_atom})
            atoms.update({i: atom})
        return atoms
