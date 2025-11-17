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


        self.lambda_vec = list(np.round(np.linspace(0, 1, 21), 3))
        self.scf_xcfun = "b3lyp"
        self.scf_basis = 'def2-svp'
        self.mute_scf = True
        self.mute_ff_build = True
        self.mm_temperature = 600
        self.mm_steps = 1000
        self.conformer_snapshots = 10
        self.mm_step_size = 0.001
        self.save_mm_traj = False
        self.scf_drv = None
        self.folder_name = 'ts_data_' + str(int(time.time()))
        self.force_conformer_search = False
        self.discont_conformer_search = False
        self.peak_conformer_search = False
        self.peak_conformer_search_range = 1
        self.mm_scan_backward = False
        self.scf_scan = False
        self.mm_conformer_equivalence_threshold = 1e-1 # kJ/mol
        self.calculate_resp = False

        self.results_file = 'ts_results.h5'

        self.sys_builder_configuration = conf = {
            "name": "vacuum",
            "bonded_integration": True,
            "soft_core_coulomb_pes": True,
            "soft_core_lj_pes": True,
            "soft_core_coulomb_int": False,
            "soft_core_lj_int": False,
        }

        self.ffbuilder = ReactionForceFieldBuilder()
        # override default options in the ffbuilder
        self.ffbuilder.calculate_resp = False

    def find_transition_state(
        self,
        reactant: Molecule | list[Molecule],
        product: Molecule | list[Molecule],
        **build_forcefields_kwargs,
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
            water_model (str): The water model used by the ffbuilder. Only has effect if there is a water molecule involved in the reaction. Defaults to "cspce".
            
        Raises:
            ff_exception: If for whatever reason the force field scan crashes, an exception is raised.

        Returns:
            molecule, dict: molecule object of the guessed transition state and a dictionary with the results of the scan.
        """
        self.results = {}
        # Build forcefields and systems

        self.build_forcefields(reactant, product, **build_forcefields_kwargs)

        # Scan MM
        self.scan_mm()

        # Scan SCF
        if self.scf_scan:
            # assert False, 'Not implemented yet'
            self.scan_scf(self.results)

        return self.results

    def build_forcefields(self,
                          product,
                          reactant,
                          constraints=None,
                          **build_forcefields_kwargs):
        if self.mute_ff_build:
            self.ostream.print_info(
                "Building forcefields. Disable mute_ff_build to see detailed output."
            )
            self.ostream.flush()
            self.ffbuilder.ostream.mute()

        self.reactant, self.product, self.forming_bonds, self.breaking_bonds, reactants, products, product_mapping = self.ffbuilder.build_forcefields(
            reactant=product,
            product=reactant,
            **build_forcefields_kwargs,
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
            self.ffbuilder._summarise_reaction(self.reactant, self.product,
                                               self.ostream)
            self.ffbuilder.ostream.mute()
        self.ostream.print_info(
            f"System has charge {self.molecule.get_charge()} and multiplicity {self.molecule.get_multiplicity()}. Provide correct values if this is wrong."
        )
        self.mol_charge = self.molecule.get_charge()
        self.mol_multiplicity = self.molecule.get_multiplicity()

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

        
        self.systems, self.topology, _ = sysbuilder.build_systems(
            self.reactant,
            self.product,
            list(self.lambda_vec),
            self.sys_builder_configuration,
            constraints,
        )
        self.ostream.print_info(
            f"Saving systems as xml to {self.folder_name}/systems")
        self.ostream.flush()
        sysbuilder.save_systems_as_xml(self.systems,
                                       self.folder_name + "/systems")
        rea_bonds = set(self.reactant.bonds.keys())
        pro_bonds = set(self.product.bonds.keys())
        static_bonds = rea_bonds & pro_bonds
        self.results.update({
            'breaking_bonds': self.breaking_bonds,
            'forming_bonds': self.forming_bonds,
            'static_bonds': static_bonds,
            'lambda_vec': self.lambda_vec,
            'reactant': self.reactant,
            'product': self.product,
        })

        return self.results

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
                # positions, V, E1, E2, E_int, N_conf = self._run_mm_scan(
                scan_dict = self._run_mm_scan(
                    self.lambda_vec,
                    rea_sim,
                    pro_sim,
                    conformer_search=True,
                    forward_init_pos=rea_init_pos,
                    backward_init_pos=pro_init_pos,
                )
            else:
                scan_dict = self._run_mm_scan(
                    self.lambda_vec,
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
                    peak_index = np.argmax(V)
                    peak_lambda = self.lambda_vec[peak_index]

                    min_index = max(
                        0,
                        peak_index - self.peak_conformer_search_range,
                    )
                    max_index = min(
                        len(self.lambda_vec),
                        peak_index + self.peak_conformer_search_range,
                    )

                    self.ostream.print_info(
                        f"Found peak MM E: {V[peak_index]:.3f} at Lambda: {peak_lambda}."
                    )
                    self.ostream.print_info(
                        f"Doing conformer search from Lambda: {self.lambda_vec[min_index]} to Lambda: {self.lambda_vec[max_index]}."
                    )
                    self.ostream.flush()

                    searched_conformers_indices.extend(
                        range(min_index, max_index + 1))
                    forward_init_pos = scan_dict[self.lambda_vec[min_index]][0]['pos']
                    backward_init_pos = scan_dict[self.lambda_vec[max_index]][0]['pos']

                    scan_dict_peak_conf = self._run_mm_scan(
                        self.lambda_vec[min_index:max_index + 1],
                        rea_sim,
                        pro_sim,
                        conformer_search=True,
                        forward_init_pos=forward_init_pos,
                        backward_init_pos=backward_init_pos,
                        skip_backward=True,
                    )
                    for l in scan_dict_peak_conf.keys():
                        scan_dict[l] += scan_dict_peak_conf[l]
                    pass

                if self.discont_conformer_search:
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
                        forward_init_pos = scan_dict[self.lambda_vec[
                            to_search_indices[0]]][0]['pos']
                        backward_init_pos = scan_dict[self.lambda_vec[
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

        except Exception as e:
            self.ostream.print_warning(f"Error in the ff scan: {e}")
            self.ostream.flush()
            exception = e

        if exception is None:
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
                f"Found highest MM E: {max_mm_energy:.3f} at Lammba: {max_mm_lambda} and conformer index: {min_mm_conformer_index}."
            )
            self.ostream.print_blank()
            self.results.update({
                'scan':
                scan_dict,
                'max_mm_xyz':
                max_mm_xyz,
                'max_mm_lambda':
                max_mm_lambda,
                'min_mm_conformer_index':
                min_mm_conformer_index,
            })
            self.molecule = Molecule.read_xyz_string(max_mm_xyz)
            self.molecule.set_multiplicity(self.mol_multiplicity)
            self.molecule.set_charge(self.mol_charge)
            self.save_results(self.results_file, self.results)
            return self.results
        else:
            self.ostream.flush()
            self.results.update({'scan': scan_dict})
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
                     forward_init_pos, backward_init_pos,skip_backward=False):
        pos = copy.copy(forward_init_pos)
        results = {}
        self._print_mm_header(lambda_vals=lambda_vals,
                              conformer_search=conformer_search)
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
            e_int = result[arg]['e_int']
            pos = result[arg]['pos']
            n_conf = len(result)

            self._print_mm_iter(l, e1, e2, v, e_int, n_conf)
            
        if self.mm_scan_backward and not skip_backward:
            self.ostream.print_info("mm_scan_backward turned on. Scanning in reverse direction.")
            self.ostream.flush()
            lambda_vals_rev = list(reversed(lambda_vals))
            self._print_mm_header(lambda_vals=lambda_vals,
                              conformer_search=conformer_search)
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
                e_int = result[arg]['e_int']
                pos = result[arg]['pos']
                n_conf = len(result)
                self._print_mm_iter(l, e1, e2, v, e_int, n_conf)
                self.ostream.flush()
                
        self.ostream.print_blank()
                
        
        for l, result in results.items():
            # remove duplicate conformers with identical MM energy 'v'
            unique_confs = []
            seen_v = set()
            for conf in result:
                v_val  = conf['v']
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
                self.ostream.print_info(f"Lambda {l}: Reduced {original_n_conf} conformers to {new_n_conf} unique conformers using equivalence threshold of {self.mm_conformer_equivalence_threshold} kJ/mol.")
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
        pdb_name = self.folder_name + f'/conf_top_{l}.pdb'

        pdb = mmapp.PDBFile.writeFile(
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

    def scan_scf(self, results=None):
        if results is None:
            results = self.results
        assert_msg_critical(
            'scan' in results.keys(),
            f'Could not find "scan" in results. Total keys: {results.keys()}',
        )

        self._print_scf_header()
        scf_energies = []
        ref = None
        max_scf_energy = None
        min_scf_conf_index = 0
        try:
            for l, scan in results['scan'].items():
                min_scf_conf_E = None
                min_conf_index = 0
                for i, conformer in enumerate(scan):
                    scf_E = self._get_scf_energy(conformer['xyz'])
                    
                    results['scan'][l][i]['scf_energy'] = scf_E
                    if math.isnan(scf_E):
                        continue
                    if min_scf_conf_E is None or scf_E < min_scf_conf_E:
                        min_scf_conf_E = scf_E
                        min_conf_index = i

                    if ref is None:
                        ref = scf_E
                    dif = scf_E - ref
                    mm_E = results['scan'][l][i]['v']

                    self._print_scf_iter(l, scf_E, mm_E, dif, i)

                if max_scf_energy is None or min_scf_conf_E > max_scf_energy:
                    max_scf_energy = min_scf_conf_E
                    max_scf_lambda = l
                    min_scf_conf_index = min_conf_index
                    max_scf_xyz = scan[min_conf_index]['xyz']

            self.ostream.print_blank()

            results = {
                'max_scf_xyz': max_scf_xyz,
                'max_scf_lambda': max_scf_lambda,
                'min_scf_conformer_index': min_scf_conf_index,
            }
            self.ostream.print_info(
                f"Found highest SCF E: {max_scf_energy:.3f} at Lambda: {max_scf_lambda} and conformer index: {min_scf_conf_index}."
            )
            self.ostream.flush()
            self.results.update(results)

            self.molecule = Molecule.read_xyz_string(max_scf_xyz)
            self.molecule.set_multiplicity(self.mol_multiplicity)
            self.molecule.set_charge(self.mol_charge)
        except Exception as e:
            self.ostream.print_warning(f"Error in the SCF scan: {e}")
            self.ostream.flush()
            self.results.update(results)
        self.save_results(self.results_file, self.results)
        return self.results

    def _get_scf_energy(self, xyz):
        self.molecule = Molecule.read_xyz_string(xyz)
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

            if not self.scf_drv.is_converged:
                self.ostream.print_warning("SCF still did not converge. Returning NaN for current calculation")
                self.ostream.flush()
                return math.nan
        return scf_results['scf_energy'] * hartree_in_kjpermol()

    @staticmethod
    def show_results(ts_results=None, filename=None, **mol_show_kwargs):
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

        ostream = OutputStream(sys.stdout)

        if ts_results is None:
            if filename is not None:
                try:
                    ostream.print_info(f"Loading results from {filename}")
                    ts_results = TransitionStateGuesser.load_results(
                        filename,
                        ostream,
                    )
                except Exception as e:
                    raise e
            else:
                raise ValueError(
                    "No results provided. Provide either ts_results or filename."
                )
        lambda_vec = [round(float(l),3) for l in ts_results['lambda_vec']]

        # if there are scf energies, get the best scf energies and everything corresponding to that
        # otherwise, get the best mm energies

        if ts_results['scan'][0][0].get('scf_energy', None) is not None:
            final_lambda = round(float(ts_results.get('max_scf_lambda', 0)),3)
        else:
            final_lambda = round(float(ts_results['max_mm_lambda']),3)
                

        forming_bonds = set(ts_results.get('forming_bonds', None))
        breaking_bonds = set(ts_results.get('breaking_bonds', None))
        bonds = set(ts_results.get('static_bonds', None))
        dashed_bonds = forming_bonds | breaking_bonds

        # lambda_slider = step=
        ipywidgets.interact(
            TransitionStateGuesser._show_iteration_inter,
            step = ipywidgets.SelectionSlider(
                options=lambda_vec,
                description='Lambda',
                value=final_lambda,
            ),
            lambda_vec=ipywidgets.fixed(lambda_vec),
            scan = ipywidgets.fixed(ts_results['scan']),
            bonds=ipywidgets.fixed(bonds),
            dashed_bonds=ipywidgets.fixed(dashed_bonds),
            **mol_show_kwargs,
        )

    @staticmethod
    def _show_iteration_inter(
        step,
        lambda_vec,
        scan,
        bonds=None,
        dashed_bonds=None,
        **mol_show_kwargs,):
        
        try:
            import ipywidgets
        except ImportError:
            raise ImportError('ipywidgets is required for this functionality.')
        
        # conformer_dropdown = 
        options = list(range(len(scan[step])))
        best_index = 0
        min_energy = None
        if scan[0][0].get('scf_energy', None) is not None:
            for i, conf in enumerate(scan[step]):
                if min_energy is None or conf['scf_energy'] < min_energy:
                    min_energy = conf['scf_energy']
                    best_index = i
        else:
            for i, conf in enumerate(scan[step]):
                if min_energy is None or conf['v'] < min_energy:
                    min_energy = conf['v']
                    best_index = i
        ipywidgets.interact(
            TransitionStateGuesser._show_iteration,
            step=ipywidgets.fixed(step),
            conformer_id = ipywidgets.Dropdown(
                options=options, 
                description='Conf. ID', 
                value=best_index,
            ),
            lambda_vec=ipywidgets.fixed(lambda_vec),
            scan = ipywidgets.fixed(scan),
            bonds=ipywidgets.fixed(bonds),
            dashed_bonds=ipywidgets.fixed(dashed_bonds),
            **mol_show_kwargs,
        )


    @staticmethod
    def _show_iteration(
        step,
        conformer_id,
        lambda_vec,
        scan,
        bonds=None,
        dashed_bonds=None,
        **mol_show_kwargs,
    ):
        """
        Show the geometry at a specific iteration.
        """

        try:
            import matplotlib.pyplot as plt
        except ImportError:
            raise ImportError('matplotlib is required for this functionality.')

        # todo add a nicer visualisation to this
        mm_energies,_,_,_ = TransitionStateGuesser._get_best_mm_E_from_scan_dict(
            scan)
        mm_min = np.min(mm_energies)
        rel_mm_energies = np.asarray(mm_energies) - mm_min
        
        lam_index = np.where(lambda_vec == np.array(step))[0][0]
        xyz_i = scan[step][conformer_id]['xyz']
        total_steps = len(rel_mm_energies) - 1
        x = np.linspace(0, lambda_vec[-1], 100)
        y = np.interp(x, lambda_vec, rel_mm_energies)
        fig, ax1 = plt.subplots(figsize=(6.5, 4))
        
        print(f"Energies for conformers:")
        if scan[0][0].get('scf_energy', None) is not None:
            scf_energies,_ = TransitionStateGuesser._get_best_scf_E_from_scan_dict(
                scan)
            scf_min = np.min(scf_energies)
            rel_scf_energies = np.asarray(scf_energies) - scf_min
            for i, conf in enumerate(scan[step]):
                print(f"  Conf. ID {i}: Relative MM energy = {conf['v']-mm_min:.3f} kJ/mol, Relative SCF energy = {conf['scf_energy']-scf_min:.3f} kJ/mol")
        else:
            rel_scf_energies = None
            for i, conf in enumerate(scan[step]):
                print(f"  Conf. ID {i}: Relative MM energy = {conf['v']-mm_min:.3f} kJ/mol")
            
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

        if rel_scf_energies is not None:
            ax2 = ax1.twinx()
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

        mol = Molecule.read_xyz_string(xyz_i)

        if bonds is not None and dashed_bonds is not None:
            mol.show(
                bonds=bonds,
                dashed_bonds=dashed_bonds,
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
    def _get_best_scf_E_from_scan_dict(scan):
        scf_energies = []
        conf_indices = []
        for conf_scan in scan.values():
            lowest_scf_e = None
            index = 0
            for i, conf in enumerate(conf_scan):
                scf_e = conf.get('scf_energy', None)
                if scf_e is not None and not math.isnan(scf_e):
                    if lowest_scf_e is None or scf_e < lowest_scf_e:
                        lowest_scf_e = scf_e
                        index = i
            scf_energies.append(lowest_scf_e)
            conf_indices.append(index)

        return scf_energies, conf_indices

    @staticmethod
    def _set_molecule_positions(molecule, positions):
        positions_au = positions / bohr_in_angstrom()
        assert molecule.number_of_atoms() == len(positions_au)
        for i in range(molecule.number_of_atoms()):
            molecule.set_atom_coordinates(i, positions_au[i])
        return molecule

    def save_results(self, fname, results):
        self.ostream.print_info(f"Saving results to {fname}")
        self.ostream.flush()
        with h5py.File(fname, 'w') as hf:
            # breaking forming static bonds
            breaking_bonds = np.array(list(results['breaking_bonds']),
                                      dtype='i')
            forming_bonds = np.array(list(results['forming_bonds']), dtype='i')
            static_bonds = np.array(list(results['static_bonds']), dtype='i')
            hf.create_dataset('breaking_bonds', data=breaking_bonds)
            hf.create_dataset('forming_bonds', data=forming_bonds)
            hf.create_dataset('static_bonds', data=static_bonds)

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

            max_scf_xyz = results.get('max_scf_xyz', None)
            max_scf_lambda = results.get('max_scf_lambda', None)
            min_scf_conformer_index = results.get('min_scf_conformer_index',
                                                  None)
            if max_scf_xyz is not None:
                hf.create_dataset('max_scf_xyz', data=[max_scf_xyz])
                hf.create_dataset('max_scf_lambda',
                                  data=max_scf_lambda,
                                  dtype='f')
                hf.create_dataset('min_scf_conformer_index',
                                  data=min_scf_conformer_index,
                                  dtype='i')

            scan_grp = hf.create_group('scan')
            for l, conf_scan in results['scan'].items():
                l_grp = scan_grp.create_group(f'{l}')
                for i, conf in enumerate(conf_scan):
                    conf_grp = l_grp.create_group(str(i))
                    # v e1 e2 e_int xyz scf
                    conf_grp.create_dataset('v', data=conf['v'], dtype='f')
                    conf_grp.create_dataset('e1', data=conf['e1'], dtype='f')
                    conf_grp.create_dataset('e2', data=conf['e2'], dtype='f')
                    conf_grp.create_dataset('e_int',
                                            data=conf['e_int'],
                                            dtype='f')
                    conf_grp.create_dataset('xyz', data=[conf['xyz']])
                    scf_e = conf.get('scf_energy', None)
                    if scf_e is not None:
                        conf_grp.create_dataset('scf_energy',
                                                data=scf_e,
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
            results['lambda_vec'] = [float(l) for l in hf['lambda_vec'][()]]

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

            # Optional SCF results
            if 'max_scf_xyz' in hf:
                results['max_scf_xyz'] = hf['max_scf_xyz'][(
                )][0].decode('utf-8')
                results['max_scf_lambda'] = hf['max_scf_lambda'][()]
                results['min_scf_conformer_index'] = hf[
                    'min_scf_conformer_index'][()]

            # Scan results
            results['scan'] = {}
            scan_grp = hf['scan']
            for l in scan_grp:
                l_grp = scan_grp[l]
                conf_scan = []
                for i in l_grp:
                    conf_grp = l_grp[i]
                    conf = {
                        'v': conf_grp['v'][()],
                        'e1': conf_grp['e1'][()],
                        'e2': conf_grp['e2'][()],
                        'e_int': conf_grp['e_int'][()],
                        'xyz': conf_grp['xyz'][()][0].decode('utf-8')
                    }
                    if 'scf_energy' in conf_grp:
                        conf['scf_energy'] = conf_grp['scf_energy'][()]
                    conf_scan.append(conf)
                results['scan'][float(l)] = conf_scan

        return results

    def _print_mm_header(self, lambda_vals=None, conformer_search=False):
        self.ostream.print_blank()
        if conformer_search:
            if lambda_vals is None:
                self.ostream.print_header(
                    "Starting MM scan with conformer search")
            else:
                self.ostream.print_header(
                    f"Starting MM scan with conformer search for lambda values {lambda_vals}"
                )
        else:
            self.ostream.print_header("Starting MM scan")
        self.ostream.print_header(f"MD steps:              {self.mm_steps:>10}")
        if conformer_search:
            conf_snapshots = self.conformer_snapshots
        else:
            conf_snapshots = 1
        self.ostream.print_header(
            f"conf. snapshots:       {conf_snapshots:>10}")
        self.ostream.print_header(
            f"MD temperature:        {self.mm_temperature:>8} K")
        self.ostream.print_header(
            f"MD step size:          {self.mm_step_size:>7} ps")
        self.ostream.print_header(f"folder name: {self.folder_name:>20}")
        self.ostream.print_header(
            f"saving MD traj:        {str(self.save_mm_traj):>10}")
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

    def _print_mm_iter(self, l, e1, e2, v, e_int, n_conf=None):

        valstr = "{:6.2f}   {:10.2f}   {:10.2f}   {:9.2f}   {:6}".format(
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
        valsltr = '{} | {} | {} | {}'.format(
            'Lambda',
            'Conf. i',
            'Rel. E (kJ/mol)',
            'MM V (kJ/mol)',
        )
        self.ostream.print_header(valsltr)
        self.ostream.print_header(60 * '-')
        self.ostream.flush()

    def _print_scf_iter(self, l, scf_E, mm_E, dif, conf_index):

        valstr = "{:6.2f}   {:7d}   {:15.3f}   {:14.3f}".format(
            l, conf_index+1, dif, mm_E)
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
        if not ts_params:
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
