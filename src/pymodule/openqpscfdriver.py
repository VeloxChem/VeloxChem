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

from contextlib import contextmanager
from copy import deepcopy
from mpi4py import MPI
import hashlib
import os
import sys
import time
import tempfile
import uuid
import numpy as np

from .veloxchemlib import mpi_master
from .outputstream import OutputStream
from .errorhandler import assert_msg_critical
from .openqpstatetracker import native_overlap_to_previous_current

try:
    import oqp
    from oqp.molecule import Molecule as OQPMolecule
    from oqp.library.single_point import (SinglePoint, Gradient, LastStep,
                                          BasisOverlap, NACME)
    from oqp.utils.input_checker import check_input_values
except ImportError:
    pass


class OpenQPScfDriver:
    """
    Implements OpenQP (PyOQP) SCF driver.

    OpenQP performs the entire electronic-structure pipeline internally
    (SCF reference, and -- for the excited-state drivers -- MRSF/SF-TDDFT and
    analytic gradients).  This driver is a thin VeloxChem wrapper around the
    ``oqp`` Python API: it translates a VeloxChem molecule and the driver
    settings into an OpenQP configuration dictionary, runs the requested stage
    through ``oqp.library.single_point.SinglePoint`` / ``Gradient``, and maps
    the resulting energies and gradients back to VeloxChem.

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.

    Instance variables
        - method: Electronic-structure method (`hf` or `dft`).
        - scf_type: OpenQP SCF reference type (`rhf`, `uhf`, or `rohf`).
        - basis: Basis set label for OpenQP.
        - functional: DFT functional label (empty string => Hartree-Fock).
        - multiplicity: SCF spin multiplicity override (None => from molecule).
        - max_cycles: Maximum number of SCF iterations.
        - scratch_dir: Base scratch directory for OpenQP files.
        - openqp_verbose: Print OpenQP log to the OpenQP log file verbosely.
    """

    def __init__(self, comm=None, ostream=None):
        """
        Initializes OpenQP SCF driver.
        """

        if comm is None:
            comm = MPI.COMM_WORLD

        if ostream is None:
            if comm.Get_rank() == mpi_master():
                ostream = OutputStream(sys.stdout)
            else:
                ostream = OutputStream(None)

        self.comm = comm
        self.rank = self.comm.Get_rank()
        self.nodes = self.comm.Get_size()
        self.ostream = ostream

        # method / model settings
        self.method = 'hf'
        self.scf_type = 'rhf'
        self.scf_mode = 'auto'
        self.basis = '6-31g*'
        self.functional = ''
        self.multiplicity = None

        # SCF convergence controls
        self.max_cycles = 100
        self.conv_thresh = None
        self.d4 = False

        # OpenMP threads for the OpenQP native kernels
        self.omp_threads = None

        self.scratch_dir = None
        self.openqp_verbose = False
        self._title_printed = False

        # cached OpenQP molecule and calculation state
        self._oqp_mol = None
        self._calc_signature = None
        self._active_geom_signature = None
        self._scf_geom_signature = None

        self._energy = None
        self._energies = None
        self._gradient = None
        self._scf_results = None

        self.print_openqp_log = True
        self.save_openqp_json = True
        self._openqp_log_file = None
        self._openqp_log_offset = 0

    @staticmethod
    def is_available():
        """
        Returns whether the OpenQP (oqp) Python module is available.
        """

        return 'oqp' in sys.modules

    # -- settings ---------------------------------------------------------

    def set_method(self, method_label):
        """
        Sets OpenQP method or DFT functional.

        :param method_label:
            'hf'/'rhf'/'uhf'/'rohf' for Hartree-Fock, 'dft' or a functional
            name (e.g. 'bhhlyp') for DFT.
        """

        if self.rank != mpi_master():
            return

        label = str(method_label).strip().lower()
        hf_aliases = {'hf', 'rhf', 'uhf', 'rohf'}
        dft_aliases = {'dft', 'rks', 'uks', 'roks'}

        if label in hf_aliases:
            self.method = 'hf'
        elif label in dft_aliases:
            self.method = 'dft'
        else:
            # Interpret unknown labels as DFT functional names.
            self.method = 'dft'
            self.functional = label

        self._invalidate_cache()

    def set_basis(self, basis_label):
        """
        Sets OpenQP basis label.
        """

        if self.rank == mpi_master():
            self.basis = str(basis_label).strip()
            self._invalidate_cache()

    def set_functional(self, functional_label):
        """
        Sets OpenQP DFT functional label.
        """

        if self.rank == mpi_master():
            self.method = 'dft'
            self.functional = str(functional_label).strip().lower()
            self._invalidate_cache()

    def set_scf_type(self, scf_type):
        """
        Sets OpenQP SCF reference type (`rhf`, `uhf`, or `rohf`).
        """

        if self.rank != mpi_master():
            return

        stype = str(scf_type).strip().lower()
        assert_msg_critical(
            stype in ('rhf', 'uhf', 'rohf'),
            f'OpenQPScfDriver: invalid scf_type "{scf_type}"')
        self.scf_type = stype
        self._invalidate_cache()

    def set_scf_mode(self, scf_mode):
        """
        Sets OpenQP SCF mode (`restricted`, `unrestricted`, or `auto`).

        This is a convenience alias mapping onto :func:`set_scf_type`.
        """

        if self.rank != mpi_master():
            return

        mode = str(scf_mode).strip().lower()
        if mode == 'auto':
            self.scf_mode = 'auto'
        elif mode in ('restricted', 'rhf', 'rks'):
            self.scf_mode = 'restricted'
            self.scf_type = 'rhf'
        elif mode in ('unrestricted', 'uhf', 'uks'):
            self.scf_mode = 'unrestricted'
            self.scf_type = 'uhf'
        else:
            assert_msg_critical(
                False, f'OpenQPScfDriver: invalid scf_mode "{scf_mode}"')
        self._invalidate_cache()

    def set_multiplicity(self, multiplicity):
        """
        Sets an explicit SCF spin multiplicity (overriding the molecule).
        """

        if self.rank == mpi_master():
            self.multiplicity = int(multiplicity)
            self._invalidate_cache()

    def set_max_cycles(self, max_cycles):
        """
        Sets the maximum number of SCF iterations.
        """

        if self.rank == mpi_master():
            self.max_cycles = int(max_cycles)
            self._invalidate_cache()

    def set_scratch_dir(self, scratch_dir):
        """
        Sets the base scratch directory for OpenQP files.
        """

        if self.rank == mpi_master():
            path = os.path.abspath(str(scratch_dir))
            self.scratch_dir = path + os.sep if not path.endswith(
                os.sep) else path
            self._invalidate_cache()

    def get_method(self):
        """
        Gets the OpenQP method label.
        """

        if self.rank == mpi_master():
            if self.method == 'dft':
                return f'dft:{self.functional}'
            return self.method
        return None

    def get_energy(self):
        """
        Gets the OpenQP SCF (reference) energy.
        """

        if self.rank == mpi_master() and self._energy is not None:
            return self._energy
        return None

    def get_gradient(self):
        """
        Gets the last computed OpenQP gradient.
        """

        if self.rank == mpi_master() and self._gradient is not None:
            return self._gradient
        return None

    # -- public compute API ----------------------------------------------

    def compute(self, molecule, basis=None):
        """
        Performs an OpenQP SCF (reference) calculation.

        :param molecule:
            The molecule.
        :param basis:
            Unused; present for interface compatibility with other drivers.

        :return:
            The OpenQP SCF results dictionary on the master rank, otherwise
            None.
        """

        errmsg = 'OpenQPScfDriver: oqp is not available. '
        errmsg += 'Please install OpenQP (pip install openqp) in this '
        errmsg += 'environment.'
        assert_msg_critical(self.is_available(), errmsg)

        if self.rank == mpi_master():
            self.print_title()
            mol = self._run_reference(molecule)
            self._energy = float(mol.energies[0])
            self._scf_results = {
                'scf_energy': self._energy,
                'scf_type': self.scf_type,
                'method': self.method,
                'functional': self.functional if self.method == 'dft' else 'HF',
                'basis_label_openqp': self.basis,
            }
            energy = self._energy
            results = dict(self._scf_results)
        else:
            energy = None
            results = None

        self._energy = self.comm.bcast(energy, root=mpi_master())

        if self.rank == mpi_master():
            self.ostream.print_blank()
            self.ostream.flush()
            return results
        return None

    def compute_gradient(self, molecule):
        """
        Performs an OpenQP ground-state (SCF) gradient calculation.

        :param molecule:
            The molecule.

        :return:
            A results dictionary with 'energy' and 'gradient' on the master
            rank, otherwise None.
        """

        errmsg = 'OpenQPScfDriver: oqp is not available. '
        errmsg += 'Please install OpenQP (pip install openqp) in this '
        errmsg += 'environment.'
        assert_msg_critical(self.is_available(), errmsg)

        if self.rank == mpi_master():
            mol = self._run_reference(molecule, need_gradient=True)
            self._energy = float(mol.energies[0])
            self._gradient = np.array(mol.grads[0], dtype=float).reshape(
                (molecule.number_of_atoms(), 3))
            results = {
                'energy': self._energy,
                'gradient': self._gradient.copy(),
            }
            energy = self._energy
            gradient = self._gradient
        else:
            results = None
            energy = None
            gradient = None

        self._energy = self.comm.bcast(energy, root=mpi_master())
        self._gradient = self.comm.bcast(gradient, root=mpi_master())

        if self.rank == mpi_master():
            return results
        return None

    # -- OpenQP engine ----------------------------------------------------

    def _run_reference(self, molecule, need_gradient=False):
        """
        Ground-state (method='hf'/'dft') driver: runs the OpenQP reference SCF
        (and, optionally, its analytic gradient) and returns the OpenQP
        molecule with ``energies``/``grads`` populated.
        """

        method = 'hf'  # OpenQP 'method' key: hf covers both HF and KS-DFT
        functional = self.functional if self.method == 'dft' else ''
        scf_type = self._effective_scf_type(molecule)
        multiplicity = self._effective_multiplicity(molecule)

        return self.run_oqp(
            molecule,
            method=method,
            functional=functional,
            scf_type=scf_type,
            multiplicity=multiplicity,
            need_gradient=need_gradient,
        )

    def run_oqp(self, molecule, *, method, functional, scf_type, multiplicity, exc_multiplicity=1,
                tdhf_type=None, nstate=1, grad_state=None, need_gradient=False):
        """
        Central OpenQP engine shared by all OpenQP drivers.

        Builds (or reuses) an ``oqp`` molecule for the requested calculation,
        synchronizes its geometry, runs the reference (and excitation, for
        method='tdhf') energy, and optionally the requested gradient.

        Returns the ``oqp`` molecule with ``energies`` and, if requested,
        ``grads`` populated.
        """

        calc_signature = self._get_calc_signature(
            molecule, method, functional, scf_type, multiplicity, tdhf_type,
            nstate)
        geom_signature = self._get_geometry_signature(molecule)

        rebuild = (self._oqp_mol is None or
                   self._calc_signature != calc_signature)

        if rebuild:
            config = self._build_config(
                molecule, method=method, functional=functional,
                scf_type=scf_type, multiplicity=multiplicity, exc_multiplicity=exc_multiplicity,
                tdhf_type=tdhf_type, nstate=nstate, grad_state=grad_state)
            self._oqp_mol = self._new_oqp_molecule(config)
            self._calc_signature = calc_signature
            self._active_geom_signature = geom_signature
            self._scf_geom_signature = None
        else:
            # Reuse the OpenQP molecule; update the gradient target and sync
            # the geometry if it changed (e.g. an MD or optimization step).
            if grad_state is not None:
                self._oqp_mol.config['properties']['grad'] = [int(grad_state)]
            if geom_signature != self._active_geom_signature:
                self._oqp_mol.update_system(
                    np.ascontiguousarray(
                        molecule.get_coordinates_in_bohr(), dtype=np.float64))
                self._active_geom_signature = geom_signature
                self._scf_geom_signature = None

        mol = self._oqp_mol

        energy_recomputed = self._scf_geom_signature != geom_signature

        if energy_recomputed:
            with self._openqp_output_context(mol):
                try:
                    SinglePoint(mol).energy()
                except SystemExit as exc:
                    assert_msg_critical(
                        False,
                        'OpenQPScfDriver: OpenQP calculation did not converge.')
                    raise exc

        if need_gradient:
            with self._openqp_output_context(mol):
                Gradient(mol).gradient()

        # OpenQP applies D4 in LastStep, not in SinglePoint or Gradient.
        if self.d4 and (energy_recomputed or need_gradient):
            with self._openqp_output_context(mol):
                d4_step = LastStep(mol)

                if need_gradient:
                    grad_list = mol.config['properties']['grad']
                else:
                    grad_list = None

                d4_energy, d4_gradient = d4_step.get_dispersion(
                    mol, grad_list)

                # Do not add the energy twice when the SCF result was cached.
                if energy_recomputed:
                    d4_step.final_energy(d4_energy)

                if need_gradient:
                    d4_step.final_grad(d4_gradient, grad_list)

        if energy_recomputed:
            self._energies = list(np.asarray(mol.energies, dtype=float))
            self._scf_geom_signature = geom_signature

        return mol
                #mol = self._oqp_mol

        ## Reference (+ excitation) energy: recompute only when the geometry
        ## changed since the last SCF on this molecule.
        #if self._scf_geom_signature != geom_signature:
        #    with self._openqp_output_context(mol):
        #        try:
        #            SinglePoint(mol).energy()
        #        except SystemExit as exc:  # OpenQP exit() on non-convergence
        #            assert_msg_critical(
        #                False, 'OpenQPScfDriver: OpenQP calculation did not '
        #                'converge.')
        #            raise exc
        #    self._energies = list(np.array(mol.energies, dtype=float))
        #    self._scf_geom_signature = geom_signature

        #if need_gradient:
        #    with self._openqp_output_context(mol):
        #        Gradient(mol).gradient()
        #
        #self._save_openqp_results(mol)
        #self._print_new_openqp_log()

        #return mol

    def run_oqp_tracked_mrsf(self,
                             molecule,
                             *,
                             functional,
                             multiplicity,
                             nstate,
                             tracker,
                             exc_multiplicity=1):
        """Runs one atomic, state-tracked MRSF energy/gradient calculation.

        The OpenQP energy workflow is split deliberately.  ``BasisOverlap``
        must align the current and previous MO spaces after the reference SCF
        but before the current MRSF response calculation.  The native state
        overlap is then assigned to persistent identities before the gradient
        root is selected.  The reference is committed only after the selected
        gradient succeeds.

        :return:
            ``(openqp_molecule, tracking_result)``.
        """

        if str(getattr(tracker, 'method_name', '')).lower() != \
                'openqp_mrsf_native_overlap':
            raise TypeError(
                'OpenQPScfDriver.run_oqp_tracked_mrsf: invalid state tracker.')
        if int(tracker.nstates) != int(nstate):
            raise ValueError(
                'OpenQPScfDriver.run_oqp_tracked_mrsf: tracker/root-count '
                'mismatch.')

        scf_type = 'rohf'
        calc_signature = self._get_calc_signature(
            molecule, 'tdhf', functional, scf_type, multiplicity, 'mrsf',
            nstate)
        geom_signature = self._get_geometry_signature(molecule)

        # A fresh OpenQP molecule is required because the native overlap
        # pipeline injects an immutable previous-geometry snapshot through the
        # back-door interface between the reference and response stages.
        #
        # No gradient root is configured yet.  The root is only known after
        # the response and the assignment, and passing a provisional one here
        # would put a stale target into the OpenQP configuration.
        config = self._build_config(
            molecule,
            method='tdhf',
            functional=functional,
            scf_type=scf_type,
            multiplicity=multiplicity,
            exc_multiplicity=exc_multiplicity,
            tdhf_type='mrsf',
            nstate=nstate,
            grad_state=None)
        mol = self._new_oqp_molecule(config)
        self._oqp_mol = mol
        self._calc_signature = calc_signature
        self._active_geom_signature = geom_signature
        self._scf_geom_signature = None

        # Seed the triplet ROHF from the last optimizer-accepted orbitals.
        # Without this the SCF restarts from an extended-Huckel guess at every
        # geometry and can converge to a different ROHF solution, which makes
        # the whole MRSF spectrum jump discontinuously between geometries that
        # differ by less than 1e-3 Angstrom.  Only the *committed* payload is
        # used, so a rejected nuclear trial never steers these orbitals, and
        # the geometry is deliberately not taken from the payload.
        had_committed_reference = bool(tracker.has_reference)
        guess_payload = tracker.committed_guess_payload()
        seeded_from_reference = guess_payload is not None
        if seeded_from_reference:
            # put_data() reports the config blocks that a tag-only snapshot
            # does not carry; that chatter belongs in the OpenQP output
            # channel, not on the VeloxChem console.
            with self._openqp_output_context(mol):
                mol.put_data(guess_payload)
            # check_input_values() rejects guess.type=json without a restart
            # file, so the type is switched after validation.  oqp.guess_json
            # builds the density from the OQP::VEC_MO_* tags just injected,
            # and guess.continue_geom stays False so the previous geometry is
            # not reused.
            mol.config['guess']['type'] = 'json'
            mol.config['guess']['continue_geom'] = False

        try:
            with self._openqp_output_context(mol):
                single_point = SinglePoint(mol)
                reference_energy = single_point.reference()

                if tracker.has_reference:
                    mol.config['properties']['back_door'] = True
                    mol.back_door = tracker.get_reference_payload()
                    BasisOverlap(mol).overlap()

                single_point.excitation(reference_energy)

                state_energies = np.asarray(
                    mol.energies[1:int(nstate) + 1], dtype=float)
                if tracker.has_reference:
                    NACME(mol).nacme()
                    # The tagarray hands back Fortran memory read in C order,
                    # so the delivered matrix is the index transpose of
                    # s_st(previous, current).  Restore the documented
                    # (previous, current) orientation before assigning.
                    state_overlap = native_overlap_to_previous_current(
                        mol.data['OQP::td_states_overlap'])
                    tracking_result = tracker.evaluate(
                        state_energies, state_overlap)
                else:
                    tracking_result = tracker.initialize_result(state_energies)

                if (not tracking_result.is_reliable and
                        tracker.failure_mode == 'error'):
                    details = '; '.join(tracking_result.warnings)
                    raise RuntimeError(
                        'OpenQP MRSF state tracking is unreliable; the '
                        f'gradient was not evaluated: {details}')

                # The gradient target is written only now, after the response
                # and the assignment, so energy and gradient always refer to
                # the same newly assigned raw root.
                selected_root = int(tracking_result.selected_raw_root)
                mol.config['properties']['grad'] = [selected_root]

                # Capture the aligned response state before Gradient runs.
                # mol.get_data() materializes plain Python lists, so this is a
                # value copy and later OpenQP mutations cannot reach it.
                snapshot_xyz = np.asarray(mol.get_system(),
                                          dtype=float).copy()
                snapshot_data = deepcopy(mol.get_data())

                Gradient(mol).gradient()

                if self.d4:
                    d4_step = LastStep(mol)
                    grad_list = mol.config['properties']['grad']
                    d4_energy, d4_gradient = d4_step.get_dispersion(
                        mol, grad_list)
                    d4_step.final_energy(d4_energy)
                    d4_step.final_grad(d4_gradient, grad_list)

        except SystemExit as exc:
            assert_msg_critical(
                False,
                'OpenQPScfDriver: tracked OpenQP calculation did not converge.')
            raise exc

        tracking_result.state_energies = np.asarray(
            mol.energies[1:int(nstate) + 1], dtype=float)
        tracking_result.seeded_from_reference = seeded_from_reference
        # Staged, not committed: the caller commits once the consumer (for a
        # relaxed optimization, geomeTRIC) has accepted this nuclear step.
        tracker.propose(tracking_result, snapshot_xyz, snapshot_data)

        if not had_committed_reference:
            # Bootstrap, committed exactly once.  The first geometry defines
            # the persistent identity and is accepted by construction: there
            # is no earlier reference to roll back to, and leaving it staged
            # would make the next geometry re-initialize the identity from raw
            # energy order instead of following the requested state.
            tracker.commit()

        self._energies = list(np.asarray(mol.energies, dtype=float))
        self._scf_geom_signature = geom_signature

        return mol, tracking_result

    def _new_oqp_molecule(self, config):
        """
        Creates and configures a fresh ``oqp`` molecule from a config dict.
        """

        self._ensure_scratch_dir()
        project = f'vlx_openqp_{uuid.uuid4().hex[:10]}'
        log = os.path.join(self.scratch_dir, project + '.log')

        mol = OQPMolecule(project, None, log, silent=1)
        mol.usempi = False
        mol.load_config(config)
        check_input_values(mol.config, input_dir=None)

        # Raise (rather than call exit()) inside OpenQP on non-convergence, so
        # VeloxChem stays in control.  The parser stores [tests] verbatim, so
        # the Python bool must be injected after validation.
        mol.config['tests']['exception'] = True

        mol.start_time = time.time()
        mol.data['OQP::log_filename'] = log

        # Initialize the OpenQP Fortran log unit.  This is mandatory: the native
        # integral/basis setup (apply_basis) writes to this unit, and without
        # the banner call it aborts trying to open an empty filename.
        with self._openqp_output_context(mol):
            oqp.oqp_banner(mol)

        self._openqp_log_file = log
        self._openqp_log_offset = 0

        return mol

    def _build_config(self, molecule, *, method, functional, scf_type,
                      multiplicity, exc_multiplicity=1, tdhf_type=None, nstate=1, grad_state=None):
        """
        Builds an OpenQP configuration dictionary (ini-style, string values)
        from the VeloxChem molecule and the requested calculation settings.
        """

        runtype = 'grad' if grad_state is not None else 'energy'

        config = {
            'input': {
                'charge': str(int(molecule.get_charge())),
                'basis': self.basis,
                'functional': functional,
                'method': method,
                'runtype': runtype,
                'system': self._system_block(molecule),
                'd4': str(bool(self.d4)),
            },
            'guess': {
                'type': 'huckel',
            },
            'scf': {
                'type': scf_type,
                'multiplicity': str(int(multiplicity)),
                'maxit': str(int(self.max_cycles)),
            },
            'properties': {
                'grad': str(int(grad_state)) if grad_state is not None else '0',
            },
        }

        if self.conv_thresh is not None:
            config['scf']['conv'] = repr(float(self.conv_thresh))

        if self.omp_threads is not None:
            config['input']['omp_threads'] = str(int(self.omp_threads))

        if method == 'tdhf':
            config['tdhf'] = {
                'type': str(tdhf_type),
                'nstate': str(int(nstate)),
                'multiplicity':str(int(exc_multiplicity))
            }

        return config

    @staticmethod
    def _system_block(molecule):
        """
        Builds the OpenQP inline ``system`` coordinate block.

        OpenQP's ``read_system`` treats a non-empty first line as an xyz file
        name, so an inline block must start with an empty line.  Coordinates
        are given in Angstrom (OpenQP converts them to Bohr internally).
        """

        labels = list(molecule.get_labels())
        coords_ang = np.array(molecule.get_coordinates_in_angstrom(),
                              dtype=float).reshape((-1, 3))
        lines = ['']
        for label, (x, y, z) in zip(labels, coords_ang):
            lines.append(f'{label} {x:.12f} {y:.12f} {z:.12f}')
        return '\n'.join(lines)

    # -- helpers ----------------------------------------------------------

    def _effective_scf_type(self, molecule):
        if self.scf_type in ('rhf', 'uhf', 'rohf'):
            # honour an explicit closed-/open-shell request but fall back to
            # unrestricted for an open-shell molecule requested as rhf.
            if self.scf_type == 'rhf' and self._effective_multiplicity(
                    molecule) != 1:
                return 'uhf'
            return self.scf_type
        return 'rhf'

    def _effective_multiplicity(self, molecule):
        if self.multiplicity is not None:
            return int(self.multiplicity)
        return int(molecule.get_multiplicity())

    def _get_calc_signature(self, molecule, method, functional, scf_type,
                            multiplicity, tdhf_type, nstate):
        return (
            method,
            functional.lower(),
            scf_type,
            int(multiplicity),
            str(tdhf_type),
            int(nstate),
            self.basis.lower(),
            int(molecule.get_charge()),
            tuple(molecule.get_labels()),
            bool(self.d4),
        )
        #return (
        #    method,
        #    functional.lower(),
        #    scf_type,
        #    int(multiplicity),
        #    str(tdhf_type),
        #    int(nstate),
        #    self.basis.lower(),
        #    int(molecule.get_charge()),
        #    tuple(molecule.get_labels()),
        #)


    def get_openqp_log_file(self):
        """Returns the current OpenQP log filename."""

        if self.rank != mpi_master():
            return None

        if self._oqp_mol is None:
            return None

        return os.path.abspath(self._oqp_mol.log)

    @staticmethod
    def _get_geometry_signature(molecule):
        coords = np.ascontiguousarray(molecule.get_coordinates_in_bohr(),
                                      dtype=np.float64)
        hasher = hashlib.sha1()
        hasher.update(str(coords.shape).encode('ascii'))
        hasher.update(coords.tobytes())
        return hasher.hexdigest()

    def _ensure_scratch_dir(self):
        if self.scratch_dir is None:
            self.scratch_dir = tempfile.mkdtemp(prefix='vlx_openqp_') + os.sep
        elif not self.scratch_dir.endswith(os.sep):
            self.scratch_dir += os.sep
        os.makedirs(self.scratch_dir, exist_ok=True)

    @contextmanager
    def _openqp_output_context(self, mol):
        """
        Silences OpenQP's stdout chatter unless verbose output is requested.

        OpenQP writes its detailed log to ``mol.log`` regardless; this only
        controls console noise.
        """

        if self.openqp_verbose and not self.ostream.is_muted:
            yield
            return

        devnull = open(os.devnull, 'w')
        old_stdout = sys.stdout
        try:
            sys.stdout = devnull
            yield
        finally:
            sys.stdout = old_stdout
            devnull.close()

    def _invalidate_cache(self):
        """
        Resets the cached OpenQP molecule and calculation state.
        """

        self._oqp_mol = None
        self._calc_signature = None
        self._active_geom_signature = None
        self._scf_geom_signature = None
        self._energy = None
        self._energies = None
        self._gradient = None
        self._scf_results = None
        self._openqp_log_file = None
        self._openqp_log_offset = 0

    def print_title(self):
        """
        Prints the title for the OpenQP calculation (once per driver instance).
        """

        if self._title_printed:
            return
        self._title_printed = True

        self.ostream.print_blank()
        self.ostream.print_header('OpenQP SCF Driver')
        self.ostream.print_header(17 * '=')
        self.ostream.print_blank()
        self.ostream.print_reference('Reference:')
        self.ostream.print_reference(self.get_reference())
        self.ostream.flush()

    def set_print_openqp_log(self, enabled=True):
        """Echo detailed OpenQP output into the VeloxChem output."""
        if self.rank == mpi_master():
            self.print_openqp_log = bool(enabled)


    def set_save_openqp_json(self, enabled=True):
        """Write the complete OpenQP result bundle."""
        if self.rank == mpi_master():
            self.save_openqp_json = bool(enabled)


    def get_openqp_log_file(self):
        """Returns the current OpenQP log filename."""
        if self.rank == mpi_master():
            return self._openqp_log_file
        return None


    def _print_new_openqp_log(self):
        """Print only new log content, avoiding duplication during optimization."""

        if not self.print_openqp_log:
            return

        filename = self._openqp_log_file
        if filename is None or not os.path.isfile(filename):
            return

        filesize = os.path.getsize(filename)
        if filesize < self._openqp_log_offset:
            self._openqp_log_offset = 0

        with open(filename, 'rb') as logfile:
            logfile.seek(self._openqp_log_offset)
            text = logfile.read().decode('utf-8', errors='replace')
            self._openqp_log_offset = logfile.tell()

        if text:
            self.ostream.print_header('OpenQP Detailed Output')
            for line in text.splitlines():
                self.ostream.print_info(line)
            self.ostream.flush()


    def _save_openqp_results(self, mol):
        """Save OpenQP numerical results and its fully resolved configuration."""

        if not self.save_openqp_json:
            return None

        # OpenQP's full numerical bundle, including internal arrays.
        try:
            mol.save_data(lean=False)
        except TypeError:
            # Compatibility with older OpenQP versions.
            mol.save_data()

        results_file = os.path.splitext(mol.log)[0] + '.json'
        settings_file = os.path.splitext(mol.log)[0] + '.settings.json'

        # mol.config contains the validated OpenQP settings, including defaults.
        import json

        def json_default(value):
            if isinstance(value, np.ndarray):
                return value.tolist()
            if isinstance(value, np.generic):
                return value.item()
            return str(value)

        with open(settings_file, 'w', encoding='utf-8') as outfile:
            json.dump(mol.config, outfile, indent=2, default=json_default)

        return {
            'openqp_log_file': mol.log,
            'openqp_results_file': results_file,
            'openqp_settings_file': settings_file,
        }

    @staticmethod
    def get_reference():
        """
        Gets the reference string for OpenQP.
        """

        ref_str = 'S. Lee, M. Filatov, C. H. Choi et al., '
        ref_str += 'OpenQP: Open Quantum Platform for MRSF-TDDFT. '
        ref_str += 'J. Chem. Theory Comput. / J. Chem. Phys. (2024).'
        return ref_str
