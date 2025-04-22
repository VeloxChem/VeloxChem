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
from pathlib import Path

from .veloxchemlib import mpi_master, hartree_in_kjpermol
from .outputstream import OutputStream
from .errorhandler import assert_msg_critical
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

    def find_TS(self, evb, scf=True, constraints=None):
        #todo ideally find a way to efficiently deal with the input to the EVB object
        evb.temperature = self.mm_temperature
        config = evb.default_system_configurations("ts_guesser")

        evb.build_systems([config], self.lambda_vec, save_output=False, constraints=constraints)
        self.ostream.print_blank()
        self.ostream.print_header("Starting MM scan")
        self.ostream.print_info(f"Lambda vector: {self.lambda_vec}")
        self.ostream.flush()
        mm_energies, mm_geometries, ff_exception = self.scan_ff(evb)
        xyz_geometries = []
        for mm_geom in mm_geometries:
            xyz_geom = self.mm_to_xyz_geom(mm_geom, evb.reactant.molecule)
            xyz_geometries.append(xyz_geom)
        self.results = {
            'mm_energies':mm_energies,
            'mm_geometries':mm_geometries,
            'xyz_geometries':xyz_geometries,
            'lambda_vec':self.lambda_vec,
            'broken_bonds':evb.broken_bonds,
            'formed_bonds':evb.formed_bonds,
        }
        # for i, (E, P) in enumerate(zip(mm_energies, mm_geometries)):
        #     self.results[self.lambda_vec[i]] = {
        #         'mm_energy': E,
        #         'mm_geometry': P,
        #     }
        self.results['broken_bonds'] = evb.broken_bonds
        self.results['formed_bonds'] = evb.formed_bonds
        if ff_exception:
            self.ostream.print_warning(
                "The force field scan crashed. Saving results in self.results and raising exception"
            )
            raise ff_exception

        self.molecule = evb.reactant.molecule

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
            
            self.ostream.print_header(
                f"Starting SCF scan")
            self.ostream.flush()
            scf_energies = self.scan_scf()
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

            # max_scf_lambda = self.scan_scf(max_mm_index)
            # if max_scf_lambda is not None:
            #     self.ostream.print_info(
            #         f"Found highest SCF E at Lambda: {max_scf_lambda}. Returning"
            #     )
            #     self.ostream.flush()
            #     self.molecule = self.set_molecule_positions(
            #         self.molecule, self.results[max_scf_lambda]['mm_geometry'])
            #     return self.molecule, self.results
            # else:
            #     self.ostream.print_info(
            #         f"Could not find a maximum in the SCF energy. Returning mm results"
            #     )
            #     self.ostream.flush()
            #     self.molecule = self.set_molecule_positions(
            #         self.molecule, max_mm_geom)
            #     return self.molecule, self.results
        self.molecule = Molecule.read_xyz_string(self.results['final_geometry'])
        return self.molecule, self.results

    def get_constraint_string(self, geom=None):
        #todo return a string that can be directly inserted in any geometry optimisation
        pass

    def scan_ff(self, EVB):
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

    def scan_scf(self):
        scf_energies = []
        for i,l in enumerate(self.lambda_vec):
            geom = self.results['mm_geometries'][i]
            scf_E = self._get_scf_energy(geom)
            scf_energies.append(scf_E)
            self.ostream.print_info(
                f"Lambda: {l}, SCF Energy: {scf_E:.3f} kJ/mol")
            self.ostream.flush()

        self.results['scf_energies'] = scf_energies

        return scf_energies
    # def scan_scf(self, starting_index=None, starting_lambda=None):
        # if starting_index is not None:
        #     l = self.lambda_vec[starting_index]
        # elif starting_lambda is not None:
        #     assert starting_lambda in self.lambda_vec
        #     l = starting_lambda
        #     starting_index = np.where(self.lambda_vec == l)[0][0]
        # else:
        #     assert_msg_critical(
        #         False,
        #         "Either starting_index or starting_lambda must be provided")

        # if l == 1:
        #     l = self.lambda_vec[-1]
        #     lp1 = 1
        # else:
        #     lp1 = self.lambda_vec[starting_index + 1]
        # scf_E_1 = self._get_scf_energy(self.results[l]['mm_geometry'])
        # self.ostream.print_info(
        #     f"Lambda: {l}, SCF Energy: {scf_E_1:.3f} Hartree")
        # self.ostream.flush()
        # scf_E_2 = self._get_scf_energy(self.results[lp1]['mm_geometry'])
        # self.ostream.print_info(
        #     f"Lambda: {lp1}, SCF Energy: {scf_E_2:.3f} Hartree")
        # self.ostream.flush()
        # self.results[l]['scf_energy'] = scf_E_1
        # self.results[lp1]['scf_energy'] = scf_E_2

        # if scf_E_2 > scf_E_1:
        #     start = starting_index + 2
        #     direction = 1
        #     end = len(self.lambda_vec)
        #     scf_E = scf_E_2
        # else:
        #     start = starting_index - 1
        #     direction = -1
        #     end = 0
        #     scf_E = scf_E_1

        # for i in range(start, end, direction):
        #     last_scf_E = scf_E
        #     l = self.lambda_vec[i]
        #     geom = self.results[l]['mm_geometry']
        #     scf_E = self._get_scf_energy(geom)
        #     self.results[l]['scf_energy'] = scf_E

        #     self.ostream.print_info(
        #         f"Lambda: {l}, SCF Energy: {scf_E:.3f} Hartree")
        #     self.ostream.flush()
        #     self.molecule = self.set_molecule_positions(self.molecule, geom)
        #     if last_scf_E > scf_E:

        #         return self.lambda_vec[i - direction]

    def _get_scf_energy(self, positions):
        self.molecule = self.set_molecule_positions(self.molecule, positions)
        scf_drv = ScfRestrictedDriver()
        if self.mute_scf:
            scf_drv.ostream.mute()
        basis = MolecularBasis.read(self.molecule, self.scf_basis)
        scf_drv.xcfun = self.scf_xcfun
        scf_results = scf_drv.compute(self.molecule, basis)
        return scf_results['scf_energy'] * hartree_in_kjpermol()


    def show_results(self,ts_results=None, atom_indices=False):
        """
            Plot the convergence of the optimization

            :param ts_results:
                The dictionary of transition finder results.
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
            step=ipywidgets.SelectionSlider(options=lambda_vec,
                                            description='Lambda',
                                            value=final_lambda),
            atom_indices=ipywidgets.fixed(atom_indices),
        )

    def _show_iteration(self,
                    mm_energies,
                    geometries,
                    lambda_vec,
                    step,
                    scf_energies=None,
                    atom_indices=False):
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
    def mm_to_xyz_geom(geom,molecule):
        new_mol = TransitionStateGuesser.set_molecule_positions(molecule, geom)
        return new_mol.get_xyz_string()

    @staticmethod
    def set_molecule_positions(molecule, positions):
        assert molecule.number_of_atoms() == len(positions)
        for i in range(molecule.number_of_atoms()):
            molecule.set_atom_coordinates(i, positions[i])
        return molecule


