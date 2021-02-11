from mpi4py import MPI
import numpy as np
import time as tm
import contextlib
import json
import sys
import os
from pathlib import Path
import matplotlib.pyplot as plt
import MDAnalysis as mda
from MDAnalysis.topology.guessers import guess_atom_element

from .veloxchemlib import mpi_master
from .veloxchemlib import bohr_in_angstroms
from .veloxchemlib import hartree_in_ev
from .veloxchemlib import hartree_in_inverse_nm
from .molecule import Molecule
from .molecularbasis import MolecularBasis
from .scfrestdriver import ScfRestrictedDriver
from .inputparser import InputParser
from .outputstream import OutputStream
from .rspabsorption import Absorption
from .subcommunicators import SubCommunicators


class QMMMDriver:
    """
    Implements QMMM driver.

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.

    Instance variables
        - comm: The MPI communicator.
        - rank: The rank of MPI process.
        - nodes: The number of MPI processes.
        - ostream: The output stream.
        - topology_file: The topology filename.
        - trajectory_file: The trajectory filename.
        - sampling_time: The time for extracting frames from trajectory
          (default unit: picosecond).
        - qm_region: The string for selecting QM region.
        - mm_pol_region: The cutoff radius for selecting polarizable MM region.
        - mm_nonpol_region: The cutoff radius for selecting non-polarizable MM
          region.
        - filename: The filename for the calculation.
        - nstates: The number of excited states.
        - diff_tol: if the first excitation energy is not too far away from
          averaged spectrum that np.abs(excitation energy - average) >=
          diff_tol * average
        - line_profile: Gaussian or Lorentzian
        - line_param: line broadening parameter
    """

    def __init__(self, comm=None, ostream=None):
        """
        Initializes QMMM driver.
        """

        if comm is None:
            comm = MPI.COMM_WORLD

        if ostream is None:
            ostream = OutputStream(sys.stdout)

        self.comm = comm
        self.rank = self.comm.Get_rank()
        self.nodes = self.comm.Get_size()

        self.ostream = ostream

        self.topology_file = None
        self.trajectory_file = None
        self.sampling_time = np.zeros(1)

        self.qm_region = None
        self.mm_pol_region = None
        self.mm_nonpol_region = None

        self.filename = None
        self.method_dict = None
        self.description = 'N/A'

        self.nstates = 3

        self.diff_tol = 0.4

        self.line_profile = 'Gaussian'
        self.line_param = 0.015

        self.charges = None
        self.polarizabilities = None

    def update_settings(self,
                        traj_dict,
                        spect_dict,
                        rsp_dict,
                        method_dict=None):
        """
        Updates settings in qmmm driver.

        :param traj_dict:
            The input dictionary of trajectory group.
        :param spect_dict:
            The input dictionary of spectrum settings group.
        :param rsp_dict:
            The input dictionary of response settings group.
        :param method_dict:
            The input dicitonary of method settings group.
        """

        time_factor = 1.0
        length_factor = 1.0

        if 'units' in traj_dict:
            units = traj_dict['units'].replace(',', ' ').split()
            # unit for time
            if 'ps' in units:
                time_factor = 1.0
            elif 'fs' in units:
                time_factor = 1.0e-3
            elif 'ns' in units:
                time_factor = 1.0e+3
            # unit for length
            if 'angstrom' in units:
                length_factor = 1.0
            elif 'bohr' in units:
                length_factor = bohr_in_angstroms()
            elif 'nm' in units:
                length_factor = 10.0

        if 'topology_file' in traj_dict:
            self.topology_file = traj_dict['topology_file']
        if 'trajectory_file' in traj_dict:
            self.trajectory_file = traj_dict['trajectory_file']

        if 'sampling_time' in traj_dict:
            self.sampling_time = np.array(
                InputParser.parse_frequencies(traj_dict['sampling_time']))
            self.sampling_time *= time_factor

        if 'quantum_region' in traj_dict:
            self.qm_region = traj_dict['quantum_region']
        if 'classical_polarizable_region' in traj_dict:
            self.mm_pol_region = float(
                traj_dict['classical_polarizable_region'])
            self.mm_pol_region *= length_factor
        if 'classical_non-polarizable_region' in traj_dict:
            self.mm_nonpol_region = float(
                traj_dict['classical_non-polarizable_region'])
            self.mm_nonpol_region *= length_factor

        if 'filename' in traj_dict:
            self.filename = traj_dict['filename']
        if 'description' in traj_dict:
            self.description = traj_dict['description']

        if 'charges' in traj_dict:
            self.charges = traj_dict['charges']
        if 'polarizabilities' in traj_dict:
            self.polarizabilities = traj_dict['polarizabilities']

        if 'line_profile' in spect_dict:
            self.line_profile = spect_dict['line_profile'].capitalize()
        if 'broadening_parameter' in spect_dict:
            self.line_param = float(spect_dict['broadening_parameter'])
            if 'units' in spect_dict and spect_dict['units'].lower() == 'ev':
                self.line_param /= hartree_in_ev()

        if 'nstates' in rsp_dict:
            self.nstates = int(rsp_dict['nstates'])

        if method_dict is not None:
            self.method_dict = dict(method_dict)

    def compute(self, molecule, basis, min_basis):
        """
        Performs QMMM calculation.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        :param min_basis:
            The minimal AO basis set.
        """

        self.print_header()

        start_time = tm.time()

        qm_charge = molecule.get_charge()
        qm_multiplicity = molecule.get_multiplicity()

        output_dir = Path(self.filename + '_files')
        if self.rank == mpi_master():
            output_dir.mkdir(parents=True, exist_ok=True)
        self.comm.barrier()

        u = mda.Universe(self.topology_file,
                         self.trajectory_file,
                         refresh_offsets=True)

        grps = [p for p in range(self.nodes)]
        subcomm = SubCommunicators(self.comm, grps)
        local_comm = subcomm.local_comm
        cross_comm = subcomm.cross_comm
        local_rank = local_comm.Get_rank()

        local_sampling_time = self.sampling_time[self.rank::self.nodes]

        local_result = []

        # go through frames in trajectory
        for ts in u.trajectory:
            # skip frames that are not in sampling_time
            if np.min(np.abs(local_sampling_time - u.trajectory.time)) > 1e-6:
                continue

            info_text = f'Processing frame {ts.frame} '
            info_text += f'({u.trajectory.time:.3f} ps) on master node...'
            self.ostream.print_info(info_text)
            self.ostream.flush()

            if local_rank == mpi_master():

                # select QM, MM_pol and MM_nonpol regions
                qm = u.select_atoms(self.qm_region)

                mm_pol_select = f'byres around {self.mm_pol_region} group qm'
                mm_pol = u.select_atoms(mm_pol_select, qm=qm)

                mm_nonpol_select = f'byres around {self.mm_nonpol_region} group qm'
                mm_nonpol = u.select_atoms(mm_nonpol_select, qm=qm)
                mm_nonpol = mm_nonpol - mm_pol

                # crate pdb files
                with open(os.devnull, 'w') as f_devnull:
                    with contextlib.redirect_stderr(f_devnull):
                        qm.write(output_dir /
                                 f'{self.filename}_frame_{ts.frame}.pdb')

                # make QM molecule whole
                qm.unwrap()
                # shift QM, MM_pol and MM_nonpol to center of the box
                box_center = 0.5 * np.sum(ts.triclinic_dimensions, axis=0)
                shift = box_center - qm.center_of_geometry()
                qm.positions += shift
                mm_pol.positions += shift
                mm_nonpol.positions += shift
                # apply periodic boundary condition for MM
                mm_pol.pack_into_box()
                mm_nonpol.pack_into_box()
                # shift QM, MM_pol and MM_nonpol to (0,0,0)
                qm.positions -= box_center
                mm_pol.positions -= box_center
                mm_nonpol.positions -= box_center

                # create molecule
                qm_elems = [guess_atom_element(name) for name in qm.names]
                qm_mol = Molecule(qm_elems, qm.positions, 'angstrom')
                qm_mol.set_charge(qm_charge)
                qm_mol.set_multiplicity(qm_multiplicity)

                # create basis set
                basis_path = '.'
                if 'basis_path' in self.method_dict:
                    basis_path = self.method_dict['basis_path']
                basis_name = self.method_dict['basis'].upper()
                qm_basis = MolecularBasis.read(qm_mol, basis_name, basis_path)

            else:

                qm_mol = Molecule()
                qm_basis = MolecularBasis()

            # create potential file
            potfile = output_dir / f'{self.filename}_frame_{ts.frame}.pot'

            if local_rank == mpi_master():

                with open(str(potfile), 'w') as f_pot:
                    f_pot.write('@environment' + os.linesep)
                    f_pot.write('units: angstrom' + os.linesep)
                    f_pot.write('xyz:' + os.linesep)

                    res_count = 0

                    # write coordinates for the polarizable region
                    res_name = 'water'
                    for res in mm_pol.residues:
                        res_count += 1
                        for atom in res.atoms:
                            atom_label = guess_atom_element(atom.name)
                            f_pot.write(f'{atom_label:<5s}' +
                                        f'{atom.position[0]:15.8f}' +
                                        f'{atom.position[1]:15.8f}' +
                                        f'{atom.position[2]:15.8f}' +
                                        f'{res_name:>10s}' + f'{res_count:8d}' +
                                        os.linesep)

                    # write coordinate for the non-polarizable region
                    res_name = 'water-n'
                    for res in mm_nonpol.residues:
                        res_count += 1
                        for atom in res.atoms:
                            atom_label = guess_atom_element(atom.name)
                            f_pot.write(f'{atom_label:<5s}' +
                                        f'{atom.position[0]:15.8f}' +
                                        f'{atom.position[1]:15.8f}' +
                                        f'{atom.position[2]:15.8f}' +
                                        f'{res_name:>10s}' + f'{res_count:8d}' +
                                        os.linesep)

                    f_pot.write('@end' + os.linesep + os.linesep)

                    # charges
                    if self.charges:
                        f_pot.write('@charges' + os.linesep)
                        for line in self.charges:
                            f_pot.write(line + os.linesep)
                        f_pot.write('@end' + os.linesep + os.linesep)

                    # polarizabilities
                    if self.polarizabilities:
                        f_pot.write('@polarizabilities' + os.linesep)
                        for line in self.polarizabilities:
                            f_pot.write(line + os.linesep)
                        f_pot.write('@end' + os.linesep + os.linesep)

            # update method_dict with potential file
            if 'pe_options' in self.method_dict:
                self.method_dict['pe_options']['potfile'] = str(potfile)
            else:
                self.method_dict['potfile'] = str(potfile)

            # broadcast molecule and basis set
            qm_mol.broadcast(local_rank, local_comm)
            qm_basis.broadcast(local_rank, local_comm)

            # setup output stream
            output = output_dir / '{}_frame_{}.out'.format(
                self.filename, ts.frame)
            ostream = OutputStream(str(output))
            self.print_molecule_and_basis(qm_mol, qm_basis, ostream)

            # run SCF
            scf_drv = ScfRestrictedDriver(local_comm, ostream)
            scf_drv.update_settings({}, self.method_dict)
            scf_drv.compute(qm_mol, qm_basis)

            if local_rank == mpi_master():
                scf_energy = scf_drv.get_scf_energy()

            # run response for spectrum
            abs_spec = Absorption({'nstates': self.nstates}, self.method_dict)
            abs_spec.init_driver(local_comm, ostream)
            abs_spec.compute(qm_mol, qm_basis, scf_drv.scf_tensors)

            if local_rank == mpi_master():
                abs_spec.print_property(ostream)

                excitation_energies = abs_spec.get_property('eigenvalues')
                oscillator_strengths = abs_spec.get_property(
                    'oscillator_strengths')

                local_result.append((
                    ts.frame,
                    {
                        'trajectory_time': u.trajectory.time,
                        'scf_energy': scf_energy,
                        'excitation_energies': excitation_energies,
                        'oscillator_strengths': oscillator_strengths,
                    },
                ))

        self.ostream.print_blank()
        self.ostream.flush()

        if local_rank == mpi_master():
            results = cross_comm.gather(local_result, root=mpi_master())

        frame_numbers = []
        list_ex_energy = []
        list_osci_strength = []

        if self.rank == mpi_master():
            flatten_results = [tup for result in results for tup in result]
            for tup in sorted(flatten_results):
                frame = tup[0]
                time = tup[1]['trajectory_time']
                scf_energy = tup[1]['scf_energy']
                excitation_energies = tup[1]['excitation_energies']
                oscillator_strengths = tup[1]['oscillator_strengths']

                frame_numbers.append(frame)
                list_ex_energy.append(excitation_energies)
                list_osci_strength.append(oscillator_strengths)

                self.print_frame_and_time(frame, time)
                self.print_scf_energy(scf_energy)
                self.print_excited_states(excitation_energies,
                                          oscillator_strengths)

        if self.rank == mpi_master():

            # create json file for spectra data
            json_dict = {self.qm_region: {}}

            key_1 = 'description'
            json_dict[self.qm_region][key_1] = self.description

            key_2 = 'excitation energies and ocillator strengths'
            json_dict[self.qm_region][key_2] = []

            for item_1, item_2 in zip(list_ex_energy, list_osci_strength):
                json_dict[self.qm_region][key_2].append(
                    [list(item_1), list(item_2)])

            with open(str(output_dir / 'spectrum.json'), 'w') as f_json:
                f_json.write(json.dumps(json_dict, indent=4))

            # run spectrum broadening
            self.spectrum_broadening(list_ex_energy, list_osci_strength,
                                     frame_numbers, output_dir)

        # print time spent in QMMM
        valstr = '*** Time spent in QMMM calculation: '
        valstr += '{:.2f} sec ***'.format(tm.time() - start_time)
        self.ostream.print_header(valstr)
        self.ostream.print_blank()
        self.ostream.flush()

    def print_header(self):
        """
        Prints header for QMMM driver.
        """

        self.ostream.print_blank()
        self.ostream.print_header('QMMM Driver Setup')
        self.ostream.print_header(19 * '=')
        self.ostream.print_blank()

        lines = []
        lines.append(f'Topology file            :    {self.topology_file}')
        lines.append(f'Trajectory file          :    {self.trajectory_file}')
        lines.append(f'QM region                :    {self.qm_region}')
        lines.append(
            f'Pol. MM Region           :    {self.mm_pol_region:.3f} Angstrom')
        lines.append(
            f'Non-Pol. MM Region       :    {self.mm_nonpol_region:.3f} Angstrom'
        )

        lines.append(f'Spect. Line Profile      :    {self.line_profile}')
        lines.append(f'Broadening parameter     :    {self.line_param:.5f} au')

        maxlen = max([len(line) for line in lines])
        for line in lines:
            self.ostream.print_header(line.ljust(maxlen))

        self.ostream.print_blank()
        self.ostream.flush()

    def print_frame_and_time(self, frame, time):
        """
        Prints frame number and simulation time to output stream.

        :param frame:
            The frame number.
        :param time:
            The simulation time.
        """

        info_text = 'Frame: {:d}  (Time: {:.3f} ps)'.format(frame, time)
        self.ostream.print_info(info_text)
        self.ostream.flush()

    def print_molecule_and_basis(self, molecule, basis, ostream):
        """
        Prints molecule and basis set to output stream.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        :param ostream:
            The output stream.
        """

        ostream.print_block(molecule.get_string())
        ostream.print_block(molecule.more_info())
        ostream.print_blank()
        ostream.print_block(basis.get_string('Atomic Basis', molecule))
        ostream.flush()

    def print_scf_energy(self, scf_energy):
        """
        Prints SCF energy to output stream.

        :param scf_energy:
            The SCF energy.
        """

        self.ostream.print_info('SCF energy: {:.10f} a.u.'.format(scf_energy))
        self.ostream.flush()

    def print_excited_states(self, excitation_energies, oscillator_strengths):
        """
        Prints excited states to output stream.

        :param excitation_energies:
            The excitation energies.
        :param oscillator_strengths:
            The oscillator strengths.
        """

        for s in range(self.nstates):
            e = excitation_energies[s] * hartree_in_ev()
            f = oscillator_strengths[s]
            valstr = 'Excited state S{:d}:'.format(s + 1)
            valstr += '{:12.6f} eV     Osc.Str. {:10.4f}'.format(e, f)
            self.ostream.print_info(valstr)
        self.ostream.print_blank()
        self.ostream.flush()

    def spectrum_broadening(self, list_ex_energy, list_osci_strength,
                            frame_numbers, output_dir):
        """
        Calculate the spectrum with either a Gaussian or Lourenzian line profile

        :param list_ex_energy:
            The list of excitation energies in each frame.
        :param list_osci_strength:
            The list of oscillator strengths in each frame.
        :param frame_numbers:
            The list of frame indices.
        :param output_dir:
            Path of the output directory.
        """

        x_min = np.amin(list_ex_energy) - 0.02
        x_max = np.amax(list_ex_energy) + 0.02
        x = np.arange(x_min, x_max, 0.001)
        y = np.zeros((len(list_ex_energy), len(x)))

        # go through the frames and calculate the spectrum for each frame
        for i in range(len(list_ex_energy)):
            for xp in range(len(x)):
                for e, f in zip(list_ex_energy[i], list_osci_strength[i]):
                    if self.line_profile == 'Gaussian':
                        y[i][xp] += f * np.exp(-(
                            (e - x[xp]) / self.line_param)**2)
                    elif self.line_profile == 'Lorentzian':
                        y[i][xp] += 0.5 * self.line_param * f / (np.pi * (
                            (x[xp] - e)**2 + 0.25 * self.line_param**2))

        # x_max: the corresponding excitation energy
        # absorption_max: the largest peak
        x_max, absorption_max = self.plot_spectra(x, y, output_dir, 'initial')

        # Decide which frames to be deleted
        # Depending on whether the first absorption lies too far away from average
        del_indicies = []
        for i in range(len(y)):
            if np.abs(list_ex_energy[i][0] - x_max) >= self.diff_tol * x_max:
                del_indicies.append(i)

        if del_indicies:
            self.ostream.print_blank()
            for i in del_indicies:
                self.ostream.print_info(
                    f'Frame {frame_numbers[i]} is very distorted, check geometry'
                )
            self.ostream.print_blank()
            self.ostream.flush()

            y = np.delete(y, del_indicies, axis=0)

        self.plot_spectra(x, y, output_dir, 'final')

    def plot_spectra(self, x, y, output_dir, label):
        """
        calculate & plot averaged spectra, return peak value

        :param x:
            The array for the horizontal axis.
        :param y:
            The multi-line array for the vertical axis.

        :return:
            The excitation energy and absorption at the peak.
        """

        # calculate the average spectrum
        y_averaged = np.sum(y, axis=0) / y.shape[0]

        # obtain absorption_max and x_max
        x_max, absorption_max = x[0], y_averaged[0]
        for a, b in zip(x, y_averaged):
            if b > absorption_max:
                absorption_max = b
                x_max = a

        hartee_to_nm = 1.0 / hartree_in_inverse_nm()

        # plot spectra in output_dir folder
        plt.figure(figsize=(8, 4))
        plt.plot(hartee_to_nm / x, y_averaged)
        plt.ylabel('Absorption')
        plt.xlabel('Wavelength (nm)')
        plt.savefig(output_dir / 'average-spec_{}nm-{}.pdf'.format(
            round(hartee_to_nm / x_max, 2), label))

        return x_max, absorption_max
