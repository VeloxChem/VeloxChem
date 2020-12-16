import numpy as np
import time as tm
import MDAnalysis as mda
from pathlib import Path

from .molecule import Molecule
from .scfrestdriver import ScfRestrictedDriver
from .inputparser import InputParser
from .outputstream import OutputStream
from .rspabsorption import Absorption
from .veloxchemlib import hartree_in_ev


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
        - ostream: The output stream.
        - tpr_file: The tpr filename.
        - xtc_file: The xtc filename.
        - sampling_time: The time (in picosecond) for extracting frames from
          trajectory.
        - qm_selection: The string for selecting QM region.
        - mm_pol_selection: The string for selecting polarizable MM region.
        - mm_nonpol_selection: The string for selecting non-polarizable MM
          region.
        - filename: The filename for the calculation.
        - nstates: The number of excited states.
    """

    def __init__(self, comm, ostream):
        """
        Initializes QMMM driver.
        """

        self.comm = comm
        self.rank = comm.Get_rank()
        self.ostream = ostream

        self.tpr_file = None
        self.xtc_file = None
        self.sampling_time = np.zeros(1)

        self.qm_selection = None
        self.mm_pol_selection = None
        self.mm_nonpol_selection = None

        self.filename = None
        self.method_dict = None
        self.nstates = 3

    def update_settings(self, qmmm_dict, method_dict=None):
        """
        Updates settings in qmmm driver.

        :param qmmm_dict:
            The input dictionary of qmmm group.
        :param method_dict:
            The input dicitonary of method settings group.
        """

        if 'tpr_file' in qmmm_dict:
            self.tpr_file = qmmm_dict['tpr_file']
        if 'xtc_file' in qmmm_dict:
            self.xtc_file = qmmm_dict['xtc_file']
        if 'sampling_time' in qmmm_dict:
            self.sampling_time = np.array(
                InputParser.parse_frequencies(qmmm_dict['sampling_time']))

        if 'qm_selection' in qmmm_dict:
            self.qm_selection = qmmm_dict['qm_selection']
        if 'mm_pol_selection' in qmmm_dict:
            self.mm_pol_selection = qmmm_dict['mm_pol_selection']
        if 'mm_nonpol_selection' in qmmm_dict:
            self.mm_nonpol_selection = qmmm_dict['mm_nonpol_selection']

        if 'filename' in qmmm_dict:
            self.filename = qmmm_dict['filename']
        if method_dict is not None:
            self.method_dict = dict(method_dict)
        if 'nstates' in qmmm_dict:
            self.nstates = int(qmmm_dict['nstates'])

    def compute(self, molecule, basis, min_basis):
        """
        Performs QMMM calculation.
        """

        self.print_header()
        start_time = tm.time()

        qm_elems = molecule.get_labels()
        qm_charge = molecule.get_charge()
        qm_multiplicity = molecule.get_multiplicity()

        output_dir = Path(self.filename + '_files')
        output_dir.mkdir(parents=True, exist_ok=True)

        u = mda.Universe(self.tpr_file, self.xtc_file, refresh_offsets=True)

        # go through frames in trajectory
        for ts in u.trajectory:

            # skip frames that are not in sampling_time
            if np.min(np.abs(self.sampling_time - u.trajectory.time)) > 1e-6:
                continue

            self.print_frame_and_time(ts.frame, u.trajectory.time)

            # select QM, MM_pol and MM_nonpol regions
            qm = u.select_atoms(self.qm_selection)
            if 'group qm' in self.mm_pol_selection:
                mm_pol = u.select_atoms(self.mm_pol_selection, qm=qm)
            else:
                mm_pol = u.select_atoms(self.mm_pol_selection)
            if 'group qm' in self.mm_nonpol_selection:
                mm_nonpol = u.select_atoms(self.mm_nonpol_selection, qm=qm)
            else:
                mm_nonpol = u.select_atoms(self.mm_nonpol_selection)
            mm_nonpol = mm_nonpol - mm_pol

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
            qm_mol = Molecule(qm_elems, qm.positions, 'angstrom')
            qm_mol.set_charge(qm_charge)
            qm_mol.set_multiplicity(qm_multiplicity)

            # create potential file
            potfile = output_dir / '{}_frame_{}.pot'.format(
                self.filename, ts.frame)
            # TODO: write potential file for MM region
            # with open(str(potfile), 'w') as f_pot:
            #     ...

            # update method_dict with potential file
            if Path(potfile).is_file():
                self.method_dict['potfile'] = str(potfile)
            else:
                self.method_dict.pop('potfile', None)

            # setup output stream
            output = output_dir / '{}_frame_{}.out'.format(
                self.filename, ts.frame)
            ostream = OutputStream(str(output))
            self.print_molecule_and_basis(qm_mol, basis, ostream)

            # run SCF
            scf_drv = ScfRestrictedDriver(self.comm, ostream)
            scf_drv.update_settings({}, self.method_dict)
            scf_drv.compute(qm_mol, basis)

            scf_energy = scf_drv.get_scf_energy()
            self.print_scf_energy(scf_energy)

            # run response for spectrum
            abs_spec = Absorption({'nstates': self.nstates}, self.method_dict)
            abs_spec.init_driver(self.comm, ostream)
            abs_spec.compute(qm_mol, basis, scf_drv.scf_tensors)
            abs_spec.print_property(ostream)

            excitation_energies = abs_spec.get_property('eigenvalues')
            oscillator_strengths = abs_spec.get_property('oscillator_strengths')
            self.print_excited_states(excitation_energies, oscillator_strengths)

        # TODO: compute average spectrum

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
        lines.append('TPR file                 :    ' + self.tpr_file)
        lines.append('XTC file                 :    ' + self.xtc_file)
        lines.append('QM selection             :    ' + self.qm_selection)
        lines.append('Pol. MM selection        :    ' + self.mm_pol_selection)
        lines.append('Non-Pol. MM selection    :    ' +
                     self.mm_nonpol_selection)

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
        self.ostream.print_blank()
        self.ostream.print_info(info_text)
        self.ostream.print_blank()
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
        self.ostream.print_blank()
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
