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
from pathlib import Path

from .veloxchemlib import mpi_master, hartree_in_kjpermol
from .outputstream import OutputStream
from .errorhandler import assert_msg_critical
from .sanitychecks import molecule_sanity_check
from .mmforcefieldgenerator import MMForceFieldGenerator
from .evbdriver import EvbDriver
from .molecule import Molecule
from .scfrestdriver import ScfRestrictedDriver
from .molecularbasis import MolecularBasis

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
        self.lambda_vec = np.round(np.linspace(0, 1, 21), 3)
        self.scf_xcfun = "b3lyp"
        self.scf_basis = 'def2-svp'
        self.mute_scf = True
        self.mm_temperature = 600
        self.mm_steps = 1000
        self.mm_step_size = 0.001 * mmunit.picoseconds
        self.save_mm_traj = False
        self.scf_drv = None
        self.evb_drv: None | EvbDriver = None
        self.molecule = None

    def find_TS(
        self,
        evb_drv,
        scf=True,
        scf_drv=None,
        constraints=None,
        charge=None,
        multiplicity=None,
    ):
        """Find a guess for the transition state using a force field scan.

        Args:
            evb (evbDriver): An EVB driver object with built forcefields
            scf (bool, optional): If an scf energy scan should be performed. Defaults to True.
            scf_drv (scfDriver, optional): The scf driver to be used for the scf scan.
            constraints (dict, optional): Dictionary of constraints, see the EVB documentation for details. Defaults to None.
            charge (int, optional): The charge of the total system. If None (default), the charge of the evb reactant molecule will be used.
            multiplicity (int, optional): The multiplicity of the total system. If None (default), the multiplicity of the evb reactant molecule will be used.

        Raises:
            ff_exception: If for whatever reason the force field scan crashes, an exception is raised.

        Returns:
            molecule, dict: molecule object of the guessed transition state and a dictionary with the results of the scan.
        """
        self.evb_drv = evb_drv
        self.scf_drv = scf_drv
        if scf:
            molecule_sanity_check(self.evb_drv.reactant.molecule)

        self.molecule = self.evb_drv.reactant.molecule
        if charge is not None:
            self.molecule.set_charge(charge)
        if multiplicity is not None:
            self.molecule.set_multiplicity(multiplicity)
        print(
            f"System has charge {self.molecule.get_charge()} and multiplicity {self.molecule.get_multiplicity()}. Provide correct values if this is wrong."
        )

        self.evb_drv.temperature = self.mm_temperature
        self.evb_drv.build_systems(
            ['ts_guesser'],
            self.lambda_vec,
            save_output=False,
            constraints=constraints,
        )
        self.ostream.print_blank()
        self.ostream.print_header("Starting MM scan")
        self.ostream.print_info(f"Lambda vector: {self.lambda_vec}")
        self.ostream.flush()
        mm_energies, mm_geometries, ff_exception = self._scan_ff(self.evb_drv)
        xyz_geometries = []
        for mm_geom in mm_geometries:
            xyz_geom = self._mm_to_xyz_geom(mm_geom,
                                            self.evb_drv.reactant.molecule)
            xyz_geometries.append(xyz_geom)
        self.results = {
            'mm_energies': mm_energies,
            'mm_geometries': mm_geometries,
            'xyz_geometries': xyz_geometries,
            'lambda_vec': self.lambda_vec,
            'broken_bonds': self.evb_drv.broken_bonds,
            'formed_bonds': self.evb_drv.formed_bonds,
        }
        # for i, (E, P) in enumerate(zip(mm_energies, mm_geometries)):
        #     self.results[self.lambda_vec[i]] = {
        #         'mm_energy': E,
        #         'mm_geometry': P,
        #     }
        self.results['broken_bonds'] = self.evb_drv.broken_bonds
        self.results['formed_bonds'] = self.evb_drv.formed_bonds
        if ff_exception:
            self.ostream.print_warning(
                "The force field scan crashed. Saving results in self.results and raising exception"
            )
            raise ff_exception

        max_mm_index = np.argmax(mm_energies)
        max_xyz_geom = xyz_geometries[max_mm_index]
        max_mm_energy = mm_energies[max_mm_index]
        max_mm_lambda = self.lambda_vec[max_mm_index]
        self.ostream.print_info(
            f"Found highest MM E: {max_mm_energy:.3f} at Lammba: {max_mm_lambda}."
        )
        self.ostream.print_blank()
        if not scf:
            self.results.update({'final_geometry': max_xyz_geom})
            self.results.update({'final_lambda': max_mm_lambda})
        else:
            # assert False, 'Not implemented yet'

            self.ostream.print_header(f"Starting SCF scan")
            self.ostream.flush()
            scf_energies = self._scan_scf()
            max_scf_index = np.argmax(scf_energies)
            max_scf_geom = self.results['xyz_geometries'][max_scf_index]
            max_scf_energy = scf_energies[max_scf_index]
            max_scf_lambda = self.lambda_vec[max_scf_index]
            self.results.update({'final_geometry': max_scf_geom})
            self.results.update({'final_lambda': max_scf_lambda})
            self.ostream.print_info(
                f"Found highest SCF E: {max_scf_energy:.3f} at Lammba: {max_scf_lambda}."
            )
            self.ostream.flush()

        self.molecule = Molecule.read_xyz_string(self.results['final_geometry'])
        return self.molecule, self.results

    #todo add option for reading geometry (bond distances, angles, etc.) from transition state instead of averaging them
    #todo add option for recalculating charges from ts_mol
    def get_ts_ffgen(self, reaffgen=None, proffgen=None, l=0.5, ts_mol=None, recalculate = True):
        if reaffgen is None:
            reaffgen = self.evb_drv.reactant
        if proffgen is None:
            proffgen = self.evb_drv.product
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
    def _mix_dihedrals(rea_dihedrals, pro_dihedrals, l,mol):
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
    def _merge_atoms(rea_atoms, pro_atoms, l,mol):
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

    def _scan_ff(self, EVB):
        energies = []
        positions = []
        initial_positions = EVB.system_confs[0]['initial_positions']  #in nm

        topology = EVB.system_confs[0]['topology']
        if self.save_mm_traj:
            if not os.path.exists('ts_data'):
                os.makedirs('ts_data')
            else:
                for file in os.listdir('ts_data'):
                    file_path = os.path.join('ts_data', file)
                    if os.path.isfile(file_path):
                        os.unlink(file_path)
            ts_folder = Path().cwd() / 'ts_data'

            mmapp.PDBFile.writeFile(topology, initial_positions,
                                    str(ts_folder / 'topology.pdb'))
        exception = None
        mm_P = initial_positions
        try:
            for l in self.lambda_vec:
                integrator = mm.VerletIntegrator(self.mm_step_size)
                simulation = mmapp.Simulation(topology,
                                              EVB.system_confs[0]['systems'][l],
                                              integrator)
                simulation.context.setPositions(mm_P)

                simulation.minimizeEnergy()
                if self.save_mm_traj:
                    mmapp.PDBFile.writeFile(
                        topology,
                        simulation.context.getState(
                            getPositions=True).getPositions(),
                        str(ts_folder / f'{l}_begin_minim.pdb'),
                    )
                    simulation.reporters.append(
                        mmapp.XTCReporter(str(ts_folder / f'{l}_traj.xtc'), 1))
                simulation.step(self.mm_steps)
                simulation.minimizeEnergy()
                if self.save_mm_traj:
                    mmapp.PDBFile.writeFile(
                        topology,
                        simulation.context.getState(
                            getPositions=True).getPositions(),
                        str(ts_folder / f'{l}_end_minim.pdb'),
                    )

                state = simulation.context.getState(getEnergy=True,
                                                    getPositions=True)
                mm_E = state.getPotentialEnergy().value_in_unit(
                    mmunit.kilojoules_per_mole)
                mm_P = state.getPositions(asNumpy=True).value_in_unit(
                    mmunit.bohr)
                avg_x = np.mean(mm_P[:, 0])
                avg_y = np.mean(mm_P[:, 1])
                avg_z = np.mean(mm_P[:, 2])
                mm_P -= [avg_x, avg_y, avg_z]

                self.ostream.print_info(
                    f"Lambda: {l}, MM Energy: {mm_E:.3f} kJ/mol")
                self.ostream.flush()
                energies.append(mm_E)
                positions.append(mm_P)
                bohr_to_nm = 0.0529177249
                mm_P = mm_P * bohr_to_nm
        except Exception as e:
            self.ostream.print_warning(f"Error in the ff scan: {e}")
            exception = e
        return energies, positions, exception

    def _scan_scf(self):
        scf_energies = []
        for i, l in enumerate(self.lambda_vec):
            geom = self.results['mm_geometries'][i]
            scf_E = self._get_scf_energy(geom)
            scf_energies.append(scf_E)
            self.ostream.print_info(
                f"Lambda: {l}, SCF Energy: {scf_E:.3f} kJ/mol")
            self.ostream.flush()

        self.results['scf_energies'] = scf_energies

        return scf_energies

    def _get_scf_energy(self, positions):
        self.molecule = self._set_molecule_positions(self.molecule, positions)
        if self.scf_drv is None:
            scf_drv = ScfRestrictedDriver()
            scf_drv.xcfun = self.scf_xcfun
            self.scf_drv = scf_drv

        if self.mute_scf:
            self.scf_drv.ostream.mute()
        basis = MolecularBasis.read(self.molecule, self.scf_basis)
        scf_results = self.scf_drv.compute(self.molecule, basis)
        return scf_results['scf_energy'] * hartree_in_kjpermol()

    def show_results(self, ts_results=None, atom_indices=False):
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
            ts_results = self.results

        mm_energies = ts_results['mm_energies']
        geometries = ts_results['xyz_geometries']
        lambda_vec = ts_results['lambda_vec']
        scf_energies = ts_results.get('scf_energies', None)
        final_lambda = ts_results['final_lambda']
        ipywidgets.interact(
            self._show_iteration,
            mm_energies=ipywidgets.fixed(mm_energies),
            scf_energies=ipywidgets.fixed(scf_energies),
            geometries=ipywidgets.fixed(geometries),
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
        geometries,
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

        lam_index = np.where(lambda_vec == step)[0][0]
        xyz_data_i = geometries[lam_index]
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

        mol = Molecule.read_xyz_string(xyz_data_i)
        mol.show(atom_indices=atom_indices, width=640, height=360)

    @staticmethod
    def _mm_to_xyz_geom(geom, molecule):
        new_mol = TransitionStateGuesser._set_molecule_positions(molecule, geom)
        return new_mol.get_xyz_string()

    @staticmethod
    def _set_molecule_positions(molecule, positions):
        assert molecule.number_of_atoms() == len(positions)
        for i in range(molecule.number_of_atoms()):
            molecule.set_atom_coordinates(i, positions[i])
        return molecule
