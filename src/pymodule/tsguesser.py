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
import os
import math
import copy
import h5py
import time
from pathlib import Path

from .veloxchemlib import mpi_master, hartree_in_kjpermol, bohr_in_angstrom
from .outputstream import OutputStream
from .openmmdynamics import OpenMMDynamics
from .errorhandler import assert_msg_critical
from .sanitychecks import molecule_sanity_check
from .mmforcefieldgenerator import MMForceFieldGenerator
from .evbdriver import EvbDriver
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

        self.lambda_vec = list(np.round(np.linspace(0, 1, 21), 3))
        self.scf_xcfun = "b3lyp"
        self.scf_basis = 'def2-svp'
        self.mute_scf = True
        self.mute_ff_build = True
        self.mm_temperature = 600
        self.mm_steps = 1000
        self.mm_step_size = 0.001
        self.save_mm_traj = False
        self.scf_drv = None
        self.folder_name = 'ts_data_' + str(int(time.time()))
        self.force_conformer_search = False
        self.skip_discont_conformer_search = True
        self.peak_conformer_search = False
        self.peak_conformer_search_range = 0
        self.conformer_steps = 10000
        self.conformer_snapshots = 10
        self.results_file = 'ts_results.h5'

        self.ffbuilder = ReactionForceFieldBuilder()

    def find_TS(
        self,
        reactant: Molecule | list[Molecule],
        product: Molecule | list[Molecule],
        scf=True,
        **ff_kwargs,
    ):
        """Find a guess for the transition state using a force field scan.

        Args:
            reactant (Molecule | list[Molecule]): The reactant molecule or a list of reactant molecules.
            product (Molecule | list[Molecule]): The product molecule or a list of product molecules.
            scf (bool, optional): If True, performs an SCF scan after the MM scan. Defaults to True.
            reactant_partial_charges (list[float], list[list[float]]): Partial charges for the reactant. Will be calculated if not provided. Defaults to None.
            product_partial_charges (list[float], list[list[float]]): Partial charges for the product. Will be calculated if not provided. Defaults to None.
            reparameterize (bool): If True, reparameterizes unknown force constants with the Seminario method. Defaults to True
            reactant_hessians (np.ndarray, list[np.ndarray]): Hessians for the reactant for the Seminario method. Will be calculated if not provided. Defaults to None.
            product_hessians (np.ndarray, list[np.ndarray]): Hessians for the product for the Seminario method. Will be calculated if not provided. Defaults to None.
            mm_opt_constrain_bonds (list[tuple[int, int]]): Bonds to constrain during MM optimization.
            reactant_total_multiplicity (int): Total multiplicity for the reactant to override calculated value. Defaults to -1.
            product_total_multiplicity (int): Total multiplicity for the product to override calculated value. Defaults to -1.
            breaking_bonds (list[tuple[int, int]]): (List of) Bond(s) that is forced to break and is not allowed to recombine over the reaction. Defaults to None.
            mute_ff_scf (bool): If True, mutes SCF output from RESP calculations. Has no effect if mute_ff_build is True. Defaults to True.
            optimize_mol (bool): If True, does an xtb optimization of every provided molecule object before reparameterisation. Defaults to False.
            optimize_ff (bool): If True, does an mm optimization of the combined reactant and product after reparameterisation. Defaults to True.
            water_model (str): The water model used by the ffbuilder. Only has effect if there is a water molecule involved in the reaction. Defaults to "spce".
            
        Raises:
            ff_exception: If for whatever reason the force field scan crashes, an exception is raised.

        Returns:
            molecule, dict: molecule object of the guessed transition state and a dictionary with the results of the scan.
        """

        # Build forcefields and systems
        self.build_forcefields(reactant, product, **ff_kwargs)

        # Scan MM
        self.scan_mm()

        # Scan SCF
        if scf:
            # assert False, 'Not implemented yet'
            self.scan_scf(self.results)

        return self.results

    def build_forcefields(self, reactant, product, constraints=[], **ts_kwargs):
        if self.mute_ff_build:
            self.ffbuilder.ostream.mute()
            self.ostream.print_info(
                "Building forcefields. Disable mute_ff_build to see detailed output."
            )
            self.ostream.flush()
        self.ffbuilder.read_keywords(**ts_kwargs)

        self.reactant, self.product, self.forming_bonds, self.breaking_bonds, reactants, products, product_mapping = self.ffbuilder.build_forcefields(
            reactant=reactant,
            product=product,
        )

        self.molecule = Molecule.read_xyz_string(
            self.reactant.molecule.get_xyz_string())
        self.molecule.set_charge(self.reactant.molecule.get_charge())
        self.molecule.set_multiplicity(
            self.reactant.molecule.get_multiplicity())

        if self.mute_ff_build:
            self.ostream.print_blank()
            self.ostream.flush()
            self.ffbuilder.ostream.unmute()
            self.ffbuilder._summarise_reaction(self.reactant, self.product)
            self.ffbuilder.ostream.mute()
        self.ostream.print_info(
            f"System has charge {self.molecule.get_charge()} and multiplicity {self.molecule.get_multiplicity()}. Provide correct values if this is wrong."
        )
        self.mol_charge = self.molecule.get_charge()
        self.mol_multiplicity = self.molecule.get_multiplicity()

        conf = {
            "name": "vacuum",
            "bonded_integration": True,
            "soft_core_coulomb_pes": True,
            "soft_core_lj_pes": True,
            "soft_core_coulomb_int": False,
            "soft_core_lj_int": False,
        }
        sysbuilder = EvbSystemBuilder()
        if self.mute_ff_build:
            sysbuilder.ostream.mute()
            self.ostream.print_info(
                "Building MM systems for the transition state guess. Disable mute_ff_build to see detailed output."
            )
            self.ostream.flush()

        self.ostream.print_info(
            f"Saving reactant and product forcefield as json to {self.folder_name}"
        )
        self.ostream.flush()
        self.reactant.save_forcefield_as_json(
            self.reactant, f"{self.folder_name}/reactant.json")
        self.product.save_forcefield_as_json(
            self.product, f"{self.folder_name}/product.json")

        # self.initial_positions in angstrom
        self.systems, self.topology, self.initial_positions = sysbuilder.build_systems(
            self.reactant,
            self.product,
            list(self.lambda_vec),
            conf,
            constraints,
        )
        self.ostream.print_info(
            f"Saving systems as xml to {self.folder_name}/systems")
        self.ostream.flush()
        sysbuilder.save_systems_as_xml(self.systems,
                                       self.folder_name + "/systems")
        self.results.update({
            'breaking_bonds': self.breaking_bonds,
            'forming_bonds': self.forming_bonds,
            'lambda_vec': self.lambda_vec
        })

        return

    def scan_mm(self):

        self.folder = Path().cwd() / self.folder_name

        # pdbs are saved in angstrom
        if self.save_mm_traj:
            mmapp.PDBFile.writeFile(self.topology, self.initial_positions,
                                    str(self.folder / 'topology.pdb'))
        exception = None

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
        pos = self.initial_positions

        try:
            if self.force_conformer_search:
                self.ostream.print_info(
                    "force_conformer_search true. Doing conformer search at every lambda."
                )
                self.ostream.flush()
                positions, V, E1, E2, E_int, N_conf = self._run_mm_scan(
                    self.lambda_vec,
                    rea_sim,
                    pro_sim,
                    conformer_search=True,
                    init_pos=pos,
                )
            else:
                positions, V, E1, E2, E_int, N_conf = self._run_mm_scan(
                    self.lambda_vec,
                    rea_sim,
                    pro_sim,
                    conformer_search=False,
                    init_pos=pos,
                )

                # Find peak
                searched_conformers_indices = []
                if self.peak_conformer_search:
                    peak_index = np.argmax(V)
                    self.ostream.print_info(
                        f"Found peak MM E: {V[peak_index]:.3f} at Lambda: {self.lambda_vec[peak_index]}."
                    )
                    max_index = min(
                        len(self.lambda_vec),
                        peak_index + self.peak_conformer_search_range,
                    )
                    min_index = max(
                        0,
                        peak_index - self.peak_conformer_search_range,
                    )
                    self.ostream.print_info(
                        f"Doing conformer search from Lambda: {self.lambda_vec[min_index]} to Lambda: {self.lambda_vec[max_index]}."
                    )
                    self.ostream.flush()
                    searched_conformers_indices.extend(
                        range(min_index, max_index + 1))
                    positions_cs, Em_cs, E1_cs, E2_cs, E_int_cs, N_conf_cs = self._run_mm_scan(
                        self.lambda_vec[min_index:max_index + 1],
                        rea_sim,
                        pro_sim,
                        conformer_search=True,
                        init_pos=positions[min_index],
                    )

                    for i, l in enumerate(self.lambda_vec[min_index:max_index +
                                                          1]):
                        if Em_cs[i] < V[min_index + i]:
                            positions[min_index + i] = positions_cs[i]
                            V[min_index + i] = Em_cs[i]
                            E1[min_index + i] = E1_cs[i]
                            E2[min_index + i] = E2_cs[i]
                            E_int[min_index + i] = E_int_cs[i]
                            N_conf[min_index + i] = N_conf_cs[i]

                if not self.skip_discont_conformer_search:
                    discont_indices = self._check_discontinuities(E1, E2)

                    while len(discont_indices) > 0 and len(
                            searched_conformers_indices) < len(self.lambda_vec):
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
                            self.lambda_vec[i] for i in to_search_indices
                        ]
                        self.ostream.print_info(
                            f"Performing conformer search at lambda values: {to_search_lambda}."
                        )
                        self.ostream.flush()
                        positions_cs, Em_cs, E1_cs, E2_cs, E_int_cs, N_conf_cs = self._run_mm_scan(
                            to_search_lambda,
                            rea_sim,
                            pro_sim,
                            conformer_search=True,
                            init_pos=positions[to_search_indices[
                                0]],  # TODO improve initial guess when searching multiple regions
                        )

                        for i, i_l in enumerate(to_search_indices):
                            if Em_cs[i] < V[i_l]:
                                positions[i_l] = positions_cs[i]
                                V[i_l] = Em_cs[i]
                                E1[i_l] = E1_cs[i]
                                E2[i_l] = E2_cs[i]
                                E_int[i_l] = E_int_cs[i]
                                N_conf[i_l] = N_conf_cs[i]

                        discont_indices = self._check_discontinuities(E1, E2)

            pass
        except Exception as e:
            self.ostream.print_warning(f"Error in the ff scan: {e}")
            exception = e

        # Save the results and return them or raise the exception

        structures = []
        for mm_pos in positions:
            xyz = self._mm_to_xyz_str(mm_pos, self.molecule)
            structures.append(xyz)
        results = {
            'mm_energies': V,
            'mm_energies_reactant': E1,
            'mm_energies_product': E2,
            'int_energies': E_int,
            'structures': structures,
        }

        if exception is None:
            max_mm_index = np.argmax(V)
            max_mm_structure = structures[max_mm_index]
            max_mm_energy = V[max_mm_index]
            max_mm_lambda = self.lambda_vec[max_mm_index]
            self.ostream.print_info(
                f"Found highest MM E: {max_mm_energy:.3f} at Lammba: {max_mm_lambda}."
            )
            self.ostream.print_blank()
            results.update({
                'max_mm_structure': max_mm_structure,
                'max_mm_lambda': max_mm_lambda,
            })
            self.results.update(results)

            self.molecule = Molecule.read_xyz_string(max_mm_structure)
            self.molecule.set_multiplicity(self.mol_multiplicity)
            self.molecule.set_charge(self.mol_charge)
            self._save_results(self.results_file, self.results)
            return self.results
        else:
            self.ostream.flush()
            self.results.update(results)
            # self._save_results(self.results_file, self.results)
            self.ostream.print_warning(
                "The force field scan crashed. Saving results in self.results and raising exception"
            )
            raise exception

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

    def _run_mm_scan(self, lambda_vals, rea_sim, pro_sim, conformer_search,
                     init_pos):

        pos = copy.copy(init_pos)
        positions = []
        V = []
        E1 = []
        E2 = []
        E_int = []
        N_conf = []
        self._print_mm_header(conformer_search=conformer_search,
                              lambda_vals=lambda_vals)
        for l in lambda_vals:

            v, e1, e2, e_int, return_pos, n_conf = self._get_mm_energy(
                self.topology,
                self.systems[l],
                l,
                pos,
                rea_sim,
                pro_sim,
                conformer_search,
            )
            positions.append(return_pos)
            V.append(v)
            E1.append(e1)
            E2.append(e2)
            E_int.append(e_int)
            N_conf.append(n_conf)
            if conformer_search:
                self._print_mm_iter(l, e1, e2, v, e_int, n_conf)
            else:
                self._print_mm_iter(l, e1, e2, v, e_int)
            pos = return_pos

        return positions, V, E1, E2, E_int, N_conf

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
        if not conformer_search:
            integrator = mm.VerletIntegrator(self.mm_step_size *
                                             mmunit.picoseconds)
            # TODO add more flexible options if we need them, coordinate with conformergenerator as well
            # platform settings for small molecule
            platform = mm.Platform.getPlatformByName("CPU")
            platform.setPropertyDefaultValue("Threads", "1")
            simulation = mmapp.Simulation(
                topology,
                system,
                integrator,
                platform,
            )

            # units converted from angstrom to nm
            init_pos_nm = init_pos * 0.1
            simulation.context.setPositions(init_pos_nm)

            simulation.minimizeEnergy()
            simulation.context.setVelocitiesToTemperature(self.mm_temperature *
                                                          mmunit.kelvin)
            if self.save_mm_traj:
                state_pos = simulation.context.getState(
                    getPositions=True).getPositions(asNumpy=True)
                mmapp.PDBFile.writeFile(
                    topology,
                    state_pos,
                    str(self.folder / f'{l}_begin_minim.pdb'),
                )
                simulation.reporters.append(
                    mmapp.XTCReporter(str(self.folder / f'{l}_traj.xtc'), 1))
            simulation.step(self.mm_steps)
            simulation.minimizeEnergy()
            if self.save_mm_traj:
                state_pos = simulation.context.getState(
                    getPositions=True).getPositions(asNumpy=True)
                mmapp.PDBFile.writeFile(
                    topology,
                    state_pos,
                    str(self.folder / f'{l}_end_minim.pdb'),
                )

            state = simulation.context.getState(getEnergy=True,
                                                getPositions=True)
            e_int = state.getPotentialEnergy().value_in_unit(
                mmunit.kilojoules_per_mole)
            pos = state.getPositions(asNumpy=True).value_in_unit(
                mmunit.angstrom)

            n_conf = -1
        else:
            opm_dyn = OpenMMDynamics()
            opm_dyn.ostream.mute()
            # platform settings for small molecule
            opm_dyn.openmm_platform = "CPU"
            # opm_dyn.create_system_from_molecule(mol, ff_gen)
            pdb_name = self.folder_name + f'/conf_top_{l}.pdb'

            pdb = mmapp.PDBFile.writeFile(
                topology,
                init_pos * mmunit.angstrom,
                pdb_name,
            )
            opm_dyn.pdb = mmapp.PDBFile(pdb_name)
            opm_dyn.system = system
            conformers_dict = opm_dyn.conformational_sampling(
                ensemble='NVT',
                nsteps=self.conformer_steps,
                snapshots=self.conformer_snapshots,
            )
            n_conf = len(conformers_dict['energies'])
            arg = np.argmin(conformers_dict['energies'])
            e_int = conformers_dict['energies'][arg]
            temp_mol = conformers_dict['molecules'][arg]
            pos = temp_mol.get_coordinates_in_angstrom()

        v, e1, e2 = self._recalc_mm_energy(pos, l, reasim, prosim)
        avg_x = np.mean(pos[:, 0])
        avg_y = np.mean(pos[:, 1])
        avg_z = np.mean(pos[:, 2])
        pos -= [avg_x, avg_y, avg_z]

        return v, e1, e2, e_int, pos, n_conf

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

    def scan_scf(self, results=None):
        if results is None:
            results = self.results
        assert_msg_critical(
            'structures' in results.keys(),
            'Could not find "structures" in results. Total keys: {results.keys()}',
        )
        structures = results['structures']
        self._print_scf_header()
        scf_energies = []
        ref = 0
        for i, l in enumerate(self.lambda_vec):
            scf_E = self._get_scf_energy(structures[i])
            if i == 0:
                ref = scf_E
                dif = 0
            else:
                dif = scf_E - ref
            scf_energies.append(scf_E)
            self._print_scf_iter(l, scf_E, dif)
        self.ostream.print_blank()

        max_scf_index = np.argmax(scf_energies)
        max_scf_structure = structures[max_scf_index]
        max_scf_energy = scf_energies[max_scf_index]
        max_scf_lambda = self.lambda_vec[max_scf_index]
        results = {
            'scf_energies': scf_energies,
            'max_scf_structure': max_scf_structure,
            'max_scf_lambda': max_scf_lambda
        }
        self.ostream.print_info(
            f"Found highest SCF E: {max_scf_energy:.3f} at Lambda: {max_scf_lambda}."
        )
        self.ostream.flush()
        self.results.update(results)

        self.molecule = Molecule.read_xyz_string(max_scf_structure)
        self.molecule.set_multiplicity(self.mol_multiplicity)
        self.molecule.set_charge(self.mol_charge)
        self._save_results(self.results_file, self.results)
        return self.results

    def _get_scf_energy(self, structure):
        self.molecule = Molecule.read_xyz_string(structure)
        self.molecule.set_multiplicity(self.mol_multiplicity)
        self.molecule.set_charge(self.mol_charge)
        if self.scf_drv is None:
            scf_drv = ScfRestrictedDriver()
            scf_drv.xcfun = self.scf_xcfun
            self.scf_drv = scf_drv

        if self.mute_scf:
            self.scf_drv.ostream.mute()
        basis = MolecularBasis.read(self.molecule, self.scf_basis)
        scf_results = self.scf_drv.compute(self.molecule, basis)

        if not self.scf_drv.is_converged:
            self.ostream.print_warning(
                "SCF did not converge, increasing convergence threshold to 1.0e-4 and maximum itterations to 200."
            )
            self.ostream.flush()
            self.scf_drv.conv_thresh = 1.0e-4
            self.scf_drv.max_iter = 200
            scf_results = self.scf_drv.compute(self.molecule, basis)
        return scf_results['scf_energy'] * hartree_in_kjpermol()

    def show_results(self, ts_results=None, atom_indices=False, filename=None):
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
        except ImportError:
            raise ImportError('ipywidgets is required for this functionality.')

        if ts_results is None:
            if self.results is not None and self.results != {}:
                self.ostream.print_info("Loading self.results")
                ts_results = self.results
            else:
                try:
                    if filename is not None:
                        self.results_file = filename
                    self.ostream.print_info(
                        f"Loading results from {self.results_file}")
                    ts_results = self.load_results(self.results_file)
                except Exception as e:
                    raise e

        mm_energies = ts_results.get('mm_energies', None)
        structures = ts_results.get('structures', None)
        lambda_vec = ts_results.get('lambda_vec', None)
        scf_energies = ts_results.get('scf_energies', None)
        if scf_energies is not None:
            final_lambda = ts_results.get('max_scf_lambda', None)
        else:
            final_lambda = ts_results.get('max_mm_lambda', None)
        ipywidgets.interact(
            self._show_iteration,
            mm_energies=ipywidgets.fixed(mm_energies),
            scf_energies=ipywidgets.fixed(scf_energies),
            structures=ipywidgets.fixed(structures),
            lambda_vec=ipywidgets.fixed(lambda_vec),
            step=ipywidgets.SelectionSlider(
                options=lambda_vec,
                description='Lambda',
                value=final_lambda,
            ),
            atom_indices=ipywidgets.fixed(atom_indices),
        )

    def _show_iteration(
        self,
        mm_energies,
        structures,
        lambda_vec,
        step,
        scf_energies=None,
        atom_indices=False,
    ):
        """
        Show the geometry at a specific iteration.
        """

        try:
            import matplotlib.pyplot as plt
        except ImportError:
            raise ImportError('matplotlib is required for this functionality.')

        rel_mm_energies = mm_energies - np.min(mm_energies)

        lam_index = np.where(lambda_vec == np.array(step))[0][0]
        structure_i = structures[lam_index]
        total_steps = len(rel_mm_energies) - 1
        x = np.linspace(0, lambda_vec[-1], 100)
        y = np.interp(x, lambda_vec, rel_mm_energies)
        fig, ax1 = plt.subplots(figsize=(6.5, 4))
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
            lambda_vec,
            rel_mm_energies,
            color='black',
            alpha=0.7,
            s=120 / math.log(total_steps, 10),
            facecolors="none",
            edgecolor="darkcyan",
            zorder=1,
        )
        ax1.scatter(
            lambda_vec[lam_index],
            rel_mm_energies[lam_index],
            marker='o',
            color='darkcyan',
            alpha=1.0,
            s=120 / math.log(total_steps, 10),
            zorder=2,
        )
        ax1.set_xlabel(r'$\lambda$')
        ax1.set_ylabel('Relative MM energy [kJ/mol]')

        if scf_energies is not None:
            ax2 = ax1.twinx()
            rel_scf_energies = scf_energies - np.min(scf_energies)
            ax2.plot(
                x,
                np.interp(x, lambda_vec, rel_scf_energies),
                color='darkred',
                alpha=0.9,
                linewidth=2.5,
                linestyle='--',
                zorder=0,
                label='SCF energy',
            )
            ax2.scatter(
                lambda_vec,
                rel_scf_energies,
                alpha=0.7,
                s=120 / math.log(total_steps, 10),
                facecolors="none",
                edgecolor="darkorange",
                zorder=1,
            )
            ax2.scatter(
                lambda_vec[lam_index],
                rel_scf_energies[lam_index],
                marker='o',
                color='darkorange',
                alpha=1.0,
                s=120 / math.log(total_steps, 10),
                zorder=2,
            )
            ax2.set_ylabel('Relative SCF energy [kJ/mol]')

        fig.legend(loc='upper right',
                   bbox_to_anchor=(1, 1),
                   bbox_transform=ax1.transAxes)
        ax1.set_title("Transition state finder")
        # Ensure x-axis displays as integers
        ax1.set_xticks(lambda_vec[::2])
        fig.tight_layout()
        plt.show()

        mol = Molecule.read_xyz_string(structure_i)

        if self.reactant.bonds is not None and self.product.bonds is not None:
            reactant_bonds = set(self.reactant.bonds.keys())
            product_bonds = set(self.product.bonds.keys())
            changing_bonds = list(self.forming_bonds) + list(
                self.breaking_bonds)
            mol.show(
                bonds=reactant_bonds | product_bonds,
                dashed_bonds=changing_bonds,
                atom_indices=atom_indices,
                width=640,
                height=360,
            )
        else:
            mol.show(atom_indices=atom_indices, width=640, height=360)

    def _mm_to_xyz_str(self, positions, molecule=None):
        if molecule is None:
            molecule = self.molecule
        new_mol = TransitionStateGuesser._set_molecule_positions(
            molecule, positions)
        return new_mol.get_xyz_string()

    @staticmethod
    def _set_molecule_positions(molecule, positions):
        positions_au = positions / bohr_in_angstrom()
        assert molecule.number_of_atoms() == len(positions_au)
        for i in range(molecule.number_of_atoms()):
            molecule.set_atom_coordinates(i, positions_au[i])
        return molecule

    def _save_results(self, fname, results):
        self.ostream.print_info(f"Saving results to {self.results_file}")
        self.ostream.flush()
        if os.path.exists(fname):
            self.ostream.print_warning(
                f"File {fname} already exists. Overwriting it.")
            self.ostream.flush()
            os.remove(fname)
        hf = h5py.File(fname, 'w')
        lambda_vec = results.get('lambda_vec', None)
        mm_energies = results.get('mm_energies', None)
        max_mm_structure = results.get('max_mm_structure', None)
        max_mm_lambda = results.get('max_mm_lambda', None)
        scf_energies = results.get('scf_energies', None)
        max_scf_structure = results.get('max_scf_structure', None)
        max_scf_lambda = results.get('max_scf_lambda', None)

        structures = results.get('structures', None)

        reactant = self.reactant.get_forcefield_as_json(self.reactant)
        product = self.product.get_forcefield_as_json(self.product)
        rea_xyz = self.reactant.molecule.get_xyz_string()
        pro_xyz = self.product.molecule.get_xyz_string()

        forming_bonds = results.get('forming_bonds', None)
        breaking_bonds = results.get('breaking_bonds', None)

        hf.create_dataset('reactant_ff', data=reactant)
        hf.create_dataset('product_ff', data=product)
        hf.create_dataset('reactant_xyz', data=rea_xyz)
        hf.create_dataset('product_xyz', data=pro_xyz)

        hf.create_dataset('mm_energies', data=mm_energies, dtype='f')
        hf.create_dataset('lambda_vec', data=lambda_vec, dtype='f')
        hf.create_dataset('max_mm_structure', data=[max_mm_structure])
        hf.create_dataset('max_mm_lambda', data=max_mm_lambda, dtype='f')

        hf.create_dataset('forming_bonds',
                          data=np.array(np.array(list(forming_bonds)),
                                        dtype='i'))
        hf.create_dataset('breaking_bonds',
                          data=np.array(np.array(list(breaking_bonds)),
                                        dtype='i'))

        dt = h5py.string_dtype(encoding='utf-8')
        hf.create_dataset('structures', data=np.array(structures, dtype=dt))

        if scf_energies is not None:
            hf.create_dataset('scf_energies', data=scf_energies, dtype='f')
            hf.create_dataset('max_scf_structure', data=[max_scf_structure])
            hf.create_dataset('max_scf_lambda', data=max_scf_lambda, dtype='f')

    def load_results(self, fname):
        self.ostream.print_info(f"Loading results from {fname}")
        self.ostream.flush()

        hf = h5py.File(fname, 'r')
        results = {}
        results['mm_energies'] = hf['mm_energies'][:]
        results['lambda_vec'] = hf['lambda_vec'][:]
        results['max_mm_structure'] = hf['max_mm_structure'][0].decode('utf-8')
        results['max_mm_lambda'] = hf['max_mm_lambda'][()]
        results['structures'] = [s.decode('utf-8') for s in hf['structures'][:]]

        reactant_ff = hf['reactant_ff'][()]
        product_ff = hf['product_ff'][()]
        reactant_xyz = hf['reactant_xyz'][()].decode('utf-8')
        product_xyz = hf['product_xyz'][()].decode('utf-8')

        forming_bonds = hf['forming_bonds'][()]
        breaking_bonds = hf['breaking_bonds'][()]

        self.reactant = MMForceFieldGenerator.load_forcefield_from_json_string(
            reactant_ff)
        self.reactant.molecule = Molecule.read_xyz_string(reactant_xyz)
        self.product = MMForceFieldGenerator.load_forcefield_from_json_string(
            product_ff)
        self.product.molecule = Molecule.read_xyz_string(product_xyz)

        self.forming_bonds = set([tuple(bond) for bond in forming_bonds])
        self.breaking_bonds = set([tuple(bond) for bond in breaking_bonds])

        if 'scf_energies' in hf:
            results['scf_energies'] = hf['scf_energies'][:]
            results['max_scf_structure'] = hf['max_scf_structure'][0].decode(
                'utf-8')
            results['max_scf_lambda'] = hf['max_scf_lambda'][()]
        return results

    def _print_mm_header(self, conformer_search=False, lambda_vals=None):
        self.ostream.print_blank()
        if conformer_search is False:
            if lambda_vals is None:
                self.ostream.print_header("Starting MM scan")
            else:
                self.ostream.print_header(
                    f"Starting MM scan for lambda values {lambda_vals}")
            self.ostream.print_blank()
            self.ostream.print_header("MM parameters:")
            self.ostream.print_header(
                f"MD steps:              {self.mm_steps:>10}")
            self.ostream.print_header(
                f"MD temperature:        {self.mm_temperature:>8} K")
            self.ostream.print_header(
                f"MD step size:          {self.mm_step_size:>7} ps")
            self.ostream.print_header(f"folder name: {self.folder_name:>20}")
            self.ostream.print_header(
                f"saving MD traj:        {str(self.save_mm_traj):>10}")
            # self.ostream.print_header(
            #     f"conf. search:   {str(conformer_search):>10}")
            valstr = '{} | {} | {} | {}'.format(
                'Lambda',
                '    E1',
                '    E2',
                '     V',
                ' E_int',
            )
        else:
            if lambda_vals is None:
                self.ostream.print_header(
                    "Starting MM scan with conformer search")
            else:
                self.ostream.print_header(
                    f"Starting MM scan with conformer search for lambda values {lambda_vals}"
                )
            self.ostream.print_header(
                f"conf. steps:           {self.conformer_steps:>10}")
            self.ostream.print_header(
                f"conf. snapshots:       {self.conformer_snapshots:>10}")
            self.ostream.print_header(
                f"MD temperature:        {self.mm_temperature:>8} K")
            self.ostream.print_header(
                f"MD step size:          {self.mm_step_size:>7} ps")
            self.ostream.print_header(f"folder name: {self.folder_name:>20}")
            self.ostream.print_header(
                f"saving MD traj:        {str(self.save_mm_traj):>10}")
            valstr = '{} | {} | {} | {} | {}'.format(
                'Lambda',
                '    E1',
                '    E2',
                '     V',
                ' E_int',
                'n_conf',
            )
        self.ostream.print_header(valstr)
        self.ostream.print_header(45 * '-')
        self.ostream.flush()

    def _print_mm_iter(self, l, e1, e2, v, e_int, n_conf=None):
        if n_conf is None:
            valstr = "{:8.2f}  {:7.1f}  {:7.1f}  {:7.1f}".format(l, e1, e2, v)
        else:
            valstr = "{:8.2f}  {:7.1f}  {:7.1f}  {:8}  {:8}".format(
                l, e1, e2, v, n_conf)
        self.ostream.print_header(valstr)
        self.ostream.flush()

    def _print_scf_header(self):
        if self.mute_scf:
            self.ostream.print_info("Disable mute_scf to see detailed output.")
        self.ostream.print_blank()
        self.ostream.print_header(f"Starting SCF scan")
        self.ostream.print_blank()
        self.ostream.print_header("SCF parameters:")
        self.ostream.print_header(f"basis:       {self.scf_basis:>10}")
        self.ostream.print_header(f"DFT xc fun:  {self.scf_xcfun:>10}")
        self.ostream.print_blank()
        self.ostream.flush()
        valsltr = '{} | {} | {}'.format(
            'Lambda',
            'SCF Energy',
            'Difference',
        )
        self.ostream.print_header(valsltr)
        self.ostream.print_header(40 * '-')
        self.ostream.flush()

    def _print_scf_iter(self, l, scf_E, dif):
        valstr = "{:8.2f}  {:15.5f}  {:8.3f}".format(l, scf_E, dif)
        self.ostream.print_header(valstr)
        self.ostream.flush()

    #todo add option for reading geometry (bond distances, angles, etc.) from transition state instead of averaging them
    #todo add option for recalculating charges from ts_mol
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
                #todo change this to measurement from molecule
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
        if ts_params is {}:
            return None
        return ts_params

    @staticmethod
    def _mix_dihedrals(rea_dihedrals, pro_dihedrals, l, mol):
        ts_params = {}
        for dict, scaling in zip([rea_dihedrals, pro_dihedrals], [1 - l, l]):
            for id, param in dict.items():
                #todo reassign value of phase?
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
