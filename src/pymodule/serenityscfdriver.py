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

from contextlib import nullcontext
from mpi4py import MPI
import atexit
import ctypes
import hashlib
import os
import re
import shutil
import sys
import tempfile
import uuid
import weakref
import numpy as np

from .veloxchemlib import mpi_master
from .outputstream import OutputStream
from .errorhandler import assert_msg_critical

from pathlib import Path
from .molecularbasis import MolecularBasis
from .resultsio import create_hdf5, write_scf_results_to_hdf5
from .inputparser import write_unparsed_input_to_hdf5

from .firstorderprop import FirstOrderProperties
from .oneeints import compute_overlap_integrals
from .resultsio import create_hdf5, write_scf_results_to_hdf5, write_results_to_hdf5

try:
    from qcserenity import serenipy as spy
    import qcserenity as qc
except ImportError:
    pass


class SerenityCalculationError(RuntimeError):
    """
    A Serenity step that must not be used as a valid electronic point.

    Serenity reports a non-converged Davidson solver (and, with
    ``allowNotConverged``, a non-converged SCF) only as a printed warning and
    carries on with the unconverged numbers.  The VeloxChem interface raises
    this error instead, so an optimizer can never receive an energy or
    gradient from such a step.

    :param message:
        The error message.
    :param stage:
        ``'scf'``, ``'response'``, ``'gradient'``, ``'state_selection'`` or
        ``'root_identity'``.
    :param details:
        Optional diagnostics dictionary.
    """

    def __init__(self, message, stage=None, details=None):
        super().__init__(message)
        self.stage = stage
        self.details = {} if details is None else dict(details)


class AdiabaticStateSelectionError(SerenityCalculationError):
    """The computed root window does not contain the requested manifold state."""

    def __init__(self, message, details=None):
        super().__init__(message, stage='state_selection', details=details)


# Serenity's own warning texts (src/postHF/LRSCF/Tools/IterativeSolver.h,
# src/scf/Scf.cpp, src/postHF/LRSCF/Tools/LRSCFRestart.cpp).  They are the only
# convergence/restart information Serenity exposes to Python.
SERENITY_LR_NOT_CONVERGED = 'Convergence criterion not reached'
SERENITY_SCF_NOT_CONVERGED = 'SCF did NOT converge'
SERENITY_LR_RESTART_LOADED = 'Successfully loaded'
SERENITY_LR_RESTART_SCRATCH = 'Will continue from scratch'
_SERENITY_LR_ITERATIONS = re.compile(
    r'Iterative solver converged in\s+(\d+)\s+iterations')


def parse_serenity_lr_output(text):
    """
    Extracts convergence and restart information from captured LRSCF output.

    The non-convergence warning is printed at every print level, but
    "Iterative solver converged in N iterations" and the restart messages
    are silenced at print level MINIMUM.  Convergence is therefore only
    accepted on positive evidence: without either message ``converged`` is
    None, and callers must treat that as "not verified".

    :param text:
        Everything Serenity wrote to stdout during one LRSCF solve (an
        ``LRSCFTask`` or the LRSCF step inside a ``GradientTask``).

    :return:
        Dictionary with ``converged`` (True, False or None),
        ``restart_loaded`` (True, False or None if no restart message was
        printed), ``davidson_iterations`` (last reported count),
        ``n_converged_solves`` and ``warnings``.
    """

    text = '' if text is None else str(text)
    warnings = [
        line.strip() for line in text.splitlines()
        if SERENITY_LR_NOT_CONVERGED in line
    ]
    iterations = [int(match) for match in
                  _SERENITY_LR_ITERATIONS.findall(text)]
    if warnings:
        converged = False
    elif iterations:
        converged = True
    else:
        converged = None
    if SERENITY_LR_RESTART_LOADED in text:
        restart_loaded = True
    elif SERENITY_LR_RESTART_SCRATCH in text:
        restart_loaded = False
    else:
        restart_loaded = None
    return {
        'converged': converged,
        'restart_loaded': restart_loaded,
        'davidson_iterations': iterations[-1] if iterations else None,
        'n_converged_solves': len(iterations),
        'warnings': warnings,
    }


# Serenity ends an SCF as soon as any two of |dE|, rmsd(P) and ||[F,P]|| are
# below their thresholds (ConvergenceController, _nNecessaryToConverge = 2).
# The ROHF branch of Serenity's excited-state gradient (setupROHFReference)
# rejects a reference whose alpha and beta occupied spaces are not nested to
# SERENITY_ROHF_NESTING_TOLERANCE.  The nesting error of a CUHF solution
# follows the orbital gradient: CUHF-BHHLYP/6-31G* SCFs of a 38-atom azo
# compound, stopped by the default ||[F,P]|| < 5e-7, were nested to 1.3-2.6
# times ||[F,P]|| (and up to 15 times rmsd(P)), i.e. 4e-7 ... 1.3e-6, so the
# default thresholds fail the gradient check at random.  Near convergence
# |dE| < 5e-8 always holds, so whichever of the two tightened criteria below
# ends the SCF, the nesting error stays below 1e-7.  |dE| keeps its default:
# it is second order in the orbital error and scattered by 1e-8 between
# converged iterations.  The same SCF stagnated at rmsd(P) ~ 2e-9 and
# ||[F,P]|| ~ 9e-9 (1e-10/1e-10/1e-9 did not converge in 52 cycles).
SERENITY_ROHF_NESTING_TOLERANCE = 1.0e-6
CUHF_SCF_THRESHOLDS = {'rmsd': 5.0e-9, 'diis': 2.0e-8}


def occupied_space_nesting(C_alpha, C_beta, S, occ_alpha, occ_beta):
    """
    Measures how far the minority-spin occupied space lies outside the
    majority-spin occupied space.

    The singular values of C_minority,virt^T S C_majority,occ are 1 for the
    n_open open-shell orbitals, followed by sin(theta_k) for the principal
    angles theta_k between the two occupied spaces.  Serenity's ROHF gradient
    prints sigma_{n_open} and sigma_{n_open + 1} as its "open-shell overlap
    eigenvalues" and requires them to be 1 and 0 within
    SERENITY_ROHF_NESTING_TOLERANCE.  The nesting error is first order in the
    angle, whereas <S^2> - S(S + 1) = sum_k sin^2(theta_k) is second order and
    cannot detect it.

    :return:
        Dictionary with ``n_open``, ``open_shell_singular_value``
        (sigma_{n_open}), ``nesting_error`` (sigma_{n_open + 1}) and
        ``within_serenity_rohf_tolerance``.
    """

    occ_alpha = np.asarray(occ_alpha, dtype=float)
    occ_beta = np.asarray(occ_beta, dtype=float)
    if np.sum(occ_alpha) >= np.sum(occ_beta):
        c_major, occ_major, c_minor, occ_minor = (C_alpha, occ_alpha, C_beta,
                                                  occ_beta)
    else:
        c_major, occ_major, c_minor, occ_minor = (C_beta, occ_beta, C_alpha,
                                                  occ_alpha)
    major_occ = np.flatnonzero(occ_major > 0.5)
    minor_virt = np.flatnonzero(occ_minor <= 0.5)
    n_open = int(major_occ.size - np.count_nonzero(occ_minor > 0.5))

    overlap = (np.asarray(c_minor, dtype=float)[:, minor_virt].T @
               np.asarray(S, dtype=float) @
               np.asarray(c_major, dtype=float)[:, major_occ])
    sigma = (np.linalg.svd(overlap, compute_uv=False)
             if overlap.size else np.zeros(0))

    if n_open == 0:
        open_value = 1.0
    elif n_open <= sigma.size:
        open_value = float(sigma[n_open - 1])
    else:
        open_value = float('nan')
    error = float(sigma[n_open]) if sigma.size > n_open else 0.0
    return {
        'n_open': n_open,
        'open_shell_singular_value': open_value,
        'nesting_error': error,
        'within_serenity_rohf_tolerance': bool(
            abs(1.0 - open_value) <= SERENITY_ROHF_NESTING_TOLERANCE and
            error <= SERENITY_ROHF_NESTING_TOLERANCE),
    }


def _flush_c_stdio():
    """Flushes C stdio buffers; Serenity mixes printf and std::cout."""

    try:
        ctypes.CDLL(None).fflush(None)
    except Exception:
        pass


class SerenityOutputCapture:
    """
    Captures everything written to file descriptor 1 during a Serenity task.

    ``qcserenity.redirectOutputToFile`` does not flush C stdio before it
    restores the descriptor, so ``printf`` output such as the Davidson
    iteration count can be lost or emitted after the task.  This capture
    flushes both Python and C buffers on entry and exit, keeps the text for
    inspection and, when ``echo`` is set, replays it to the real stdout.

    :param directory:
        Directory for the temporary capture file; ``None`` uses TMPDIR.
    :param echo:
        Write the captured text to stdout after the task (verbose mode).
    """

    def __init__(self, directory=None, echo=False):
        self.directory = directory
        self.echo = bool(echo)
        self.text = ''
        self._file = None
        self._saved_fd = None

    def __enter__(self):
        try:
            sys.stdout.flush()
        except Exception:
            pass
        _flush_c_stdio()
        directory = self.directory
        if directory is not None and not os.path.isdir(directory):
            directory = None
        self._file = tempfile.TemporaryFile(mode='w+b', dir=directory)
        self._saved_fd = os.dup(1)
        os.dup2(self._file.fileno(), 1)
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        try:
            sys.stdout.flush()
        except Exception:
            pass
        _flush_c_stdio()
        os.dup2(self._saved_fd, 1)
        os.close(self._saved_fd)
        self._saved_fd = None
        try:
            self._file.seek(0)
            self.text = self._file.read().decode('utf-8', errors='replace')
        finally:
            self._file.close()
            self._file = None
        if self.echo and self.text:
            try:
                sys.stdout.write(self.text)
                sys.stdout.flush()
            except Exception:
                pass
        return False


class SerenityScfDriver:
    """
    Implements Serenity SCF driver.

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.

    Instance variables
        - method: Electronic-structure method (`hf` or `dft`).
        - scf_mode: SCF mode (`restricted`, `unrestricted`, or `auto`).
        - basis: Basis set label for Serenity.
        - dft_functional: DFT functional label for Serenity.
        - scratch_dir: Base scratch directory for Serenity files.
        - serenity_verbose: Print Serenity output directly to stdout.
        - scf_energy_threshold, scf_rmsd_threshold, scf_diis_threshold:
          Serenity SCF convergence thresholds; None keeps Serenity's default
          (CUHF references default to CUHF_SCF_THRESHOLDS).
    """

    def __init__(self, comm=None, ostream=None):
        """
        Initializes Serenity SCF driver.
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

        self.method = 'hf'
        self.rohf_type = None
        self.scf_mode = 'auto'
        self.basis = '6-31GS'
        self.dft_functional = 'bp86'
        
        self.dispersion = None

        #custom functional part
        self.basic_functional = None
        self.mixing_factors = None
        self.hfexchange_ratio = None
        self.lrexchange_ratio = None
        self.mu = None
        self.impl = 'xcfun'

        self.densfit_j = 'none'
        self.grid_accuracy = 7
        self.small_grid_accuracy = 7

        self.scratch_dir = None
        self.serenity_verbose = False

        # A scratch root created by this driver is owned by it and removed at
        # interpreter shutdown; one supplied through set_scratch_dir() belongs
        # to the caller and is left alone.
        self.keep_scratch_dir = False
        self._owns_scratch_dir = False
        self._scratch_cleanup_hook = None

        self._system = None # This is the systemController object from Serenity
        self._scf_task = None # This is the Task object that is executing the given Task here: SCFTask
        self._fat_task = None
        self._gradient_task = None # Uses gradient task object to perform the scf calculation

        # Variables that store information about the SCF 
        # !!! Currently needs to be reseted within the dyanmics as BasisFunctionOnGridController is giving an Error !!!
        self._system_signature = None 
        self._active_geom_signature = None
        self._last_scf_geom_signature = None
        self._last_grad_geom_signature = None
        
        self.max_cycles = 1000

        # SCF convergence thresholds (Serenity settings.scf.energyThreshold,
        # rmsdThreshold and diisThreshold).  None keeps Serenity's default,
        # except for CUHF references, which use CUHF_SCF_THRESHOLDS.
        self.scf_energy_threshold = None
        self.scf_rmsd_threshold = None
        self.scf_diis_threshold = None

        self._energy = None
        self._gradient = None
        self._scf_results = None
        self._current_scf_mode = None

        # SCF provenance.  ``_scf_revision`` increases with every SCF that
        # actually runs, so an LR solution can be tied to the exact MO
        # coefficients it was solved in (same geometry is not enough: a new
        # SCF at the same geometry may return orbitals with different phases).
        # ``_system_name`` identifies the Serenity System, whose scratch
        # directory holds the LRSCF restart file.
        self._scf_revision = 0
        self._system_name = None
        self._system_generation = 0
        self._system_scf_count = 0
        self._scf_provenance = None
        self.last_serenity_output = {}

        # h5 part of the file
        self.filename = None
        self._skip_writing_h5 = False

    @staticmethod
    def is_available():
        """
        Returns if Serenity python driver is available.
        """

        return ('qcserenity' in sys.modules and 'spy' in globals() and
                'qc' in globals())

    def set_method(self, method_label):
        """
        Sets Serenity method or DFT functional.
        """

        if self.rank != mpi_master():
            return

        label = str(method_label).strip().lower()

        if label in ('rohf', 'roks'):
            assert_msg_critical(
                False,
                f'SerenityScfDriver: "{label}" does not select an open-shell '
                'reference constraint; use set_rohf_type("CUHF") together '
                f'with set_method("{"hf" if label == "rohf" else "dft"}").')

        hf_aliases = {'hf', 'rhf', 'uhf'}
        dft_aliases = {'dft', 'rks', 'uks'}

        if label in hf_aliases:
            self.method = 'hf'
        elif label in dft_aliases:
            self.method = 'dft'
        else:
            # Interpret unknown labels as DFT functionals.
            self.method = 'dft'
            self.dft_functional = label

        self._invalidate_cache()

    def set_basis(self, basis_label):
        """
        Sets Serenity basis label.
        """

        if self.rank == mpi_master():
            self.basis = str(basis_label)
            self._invalidate_cache()

    def set_functional(self, functional_label=None, custom_functional=None):
        """
        Sets Serenity DFT functional label.
        """

        if self.rank == mpi_master():
            self.method = 'dft'
            if functional_label:
                self.dft_functional = str(functional_label).strip().lower()
            elif custom_functional:
                if 'basicfunctional' in custom_functional:
                    self.basic_functional = custom_functional.get('basicfunctional')
                else:
                    raise ValueError('Basicfunctionals missing')
                if 'mixingfactors' in custom_functional:
                    self.mixing_factors = custom_functional.get('mixingfactors')
                else:
                    raise ValueError('Missing mixingfactors for basic functionals')
                if 'hfexchangeratio' in custom_functional:
                    self.hfexchange_ratio = custom_functional.get('hfexchangeratio')
                else:
                    raise ValueError('Missing hfexchangeratio for custom functional construction')
                if 'lrexchangeratio' in custom_functional:
                    self.lrexchange_ratio = custom_functional.get('lrexchangeratio')
                else:
                    raise ValueError('Missing lrexchangeratio for custom functional construction')
                if 'mu' in custom_functional:
                    self.mu = custom_functional.get('mu')
                else:
                    raise ValueError('Missing mu for custom functional construction')

            self._invalidate_cache()

    def set_rohf_type(self, rohf_type):
        """
        Sets the open-shell reference constraint (NONE, CUHF or SUHF).

        CUHF constrains the unrestricted high-spin reference to <S^2> = S(S+1)
        exactly, which is the ROHF-equivalent reference used by OpenQP.
        """

        if self.rank != mpi_master():
            return

        if rohf_type is None:
            self.rohf_type = None
        else:
            label = str(rohf_type).strip().upper()
            assert_msg_critical(
                label in ('NONE', 'CUHF', 'SUHF'),
                f'SerenityScfDriver: Invalid rohf_type "{rohf_type}"')
            self.rohf_type = label

        self._invalidate_cache()

    def get_scf_thresholds(self):
        """
        Returns the SCF convergence thresholds passed to Serenity.

        Explicitly set thresholds win.  A CUHF reference takes the remaining
        ones from CUHF_SCF_THRESHOLDS, because the ROHF branch of Serenity's
        excited-state gradient needs the alpha and beta occupied spaces
        nested to SERENITY_ROHF_NESTING_TOLERANCE, which Serenity's default
        thresholds do not guarantee.  None keeps Serenity's default.

        :return:
            Dictionary with the ``energy``, ``rmsd`` and ``diis`` thresholds.
        """

        defaults = CUHF_SCF_THRESHOLDS if self.rohf_type == 'CUHF' else {}
        explicit = {
            'energy': self.scf_energy_threshold,
            'rmsd': self.scf_rmsd_threshold,
            'diis': self.scf_diis_threshold,
        }
        return {
            key: (float(value) if value is not None else defaults.get(key))
            for key, value in explicit.items()
        }

    def set_scf_mode(self, scf_mode):
        """
        Sets Serenity SCF mode.
        """

        if self.rank != mpi_master():
            return

        mode = str(scf_mode).strip().lower()
        if mode in ('auto',):
            self.scf_mode = 'auto'
        elif mode in ('restricted', 'rhf', 'rks'):
            self.scf_mode = 'restricted'
        elif mode in ('unrestricted', 'uhf', 'uks'):
            self.scf_mode = 'unrestricted'
        else:
            assert_msg_critical(False,
                                f'SerenityScfDriver: Invalid scf_mode "{mode}"')

        self._invalidate_cache()

    def set_scratch_dir(self, scratch_dir):
        """
        Sets base scratch directory for Serenity files.
        """

        if self.rank == mpi_master():
            path = os.path.abspath(str(scratch_dir))
            self.scratch_dir = path + os.sep if not path.endswith(
                os.sep) else path
            self._owns_scratch_dir = False
            self._invalidate_cache()

    def get_method(self):
        """
        Gets Serenity method.
        """

        if self.rank == mpi_master():
            if self.method == 'dft':
                return f'dft:{self.dft_functional}'
            return self.method
        return None

    def get_energy(self):
        """
        Gets Serenity energy.
        """

        if self.rank == mpi_master() and self._energy is not None:
            return self._energy
        return None

    def get_gradient(self):
        """
        Gets Serenity gradient.
        """

        if self.rank == mpi_master() and self._gradient is not None:
            return self._gradient
        return None

    def compute(self, molecule, basis=None):
        """
        Performs Serenity SCF calculation.

        :param molecule:
            The molecule.

        :return:
            The Serenity SCF results on master rank, otherwise None.
        """

        errmsg = 'SerenityScfDriver: qcserenity is not available. '
        errmsg += 'Please install/build Serenity python bindings.'
        assert_msg_critical(self.is_available(), errmsg)
        if self.rank == mpi_master():
            results = self._compute_energy_master(molecule)
            energy = float(results['scf_energy'])
        else:
            results = None
            energy = None

        self._energy = self.comm.bcast(energy, root=mpi_master())

        self.write_final_hdf5(self.get_final_hdf5_file(), molecule, basis)
        
        if self.rank == mpi_master():
            return dict(results)
        return None

    def get_final_hdf5_file(self):
        if self.filename is None:
            return None
        return f'{self.filename}.h5'

    def compute_gradient(self, molecule):
        """
        Performs Serenity gradient calculation.

        :param molecule:
            The molecule.

        :return:
            A results dictionary on master rank, otherwise None.
        """

        errmsg = 'SerenityScfDriver: qcserenity is not available. '
        errmsg += 'Please install/build Serenity python bindings.'
        assert_msg_critical(self.is_available(), errmsg)

        if self.rank == mpi_master():
            self._compute_gradient_master(molecule)
            results = {
                'energy': float(self._energy),
                'gradient': self._gradient.copy(),
            }
            gradient = self._gradient
        else:
            results = None
            gradient = None

        self._gradient = self.comm.bcast(gradient, root=mpi_master())

        if self.rank == mpi_master():
            return results
        return None

    def print_title(self):
        """
        Prints title for Serenity calculation.
        """

        self.ostream.print_blank()
        self.ostream.print_header('Serenity SCF Driver')
        self.ostream.print_header(19 * '=')
        self.ostream.print_blank()
        self.ostream.print_reference('Reference:')
        self.ostream.print_reference(self.get_reference())
        self.ostream.flush()

    def get_reference(self):
        """
        Gets reference string for Serenity.
        """

        ref_str = 'J. P. Unsleber, T. Dresselhaus, K. Klahr, D. Schnieders, '
        ref_str += 'M. Bockers, D. Barton, J. Neugebauer, '
        ref_str += 'J. Comput. Chem. 2018, 39, 788-798'
        return ref_str

    def _invalidate_cache(self):
        """
        Forces a new Serenity System at the next calculation.

        ``_system_signature`` is reset, so ``_ensure_system`` builds a new
        System (new name, new scratch directory) instead of moving the atoms
        of the current one.  A new System means a fresh SCF guess and no LRSCF
        restart file.  Surface-hopping measurements found that moving the
        atoms of a System on which a gradient had been evaluated shifted the
        SF reference energy (basis-function-on-grid data stayed tied to the
        old geometry), so callers that evaluate gradients reset the System
        after each evaluation unless same-System reuse has been validated.
        """

        # self._system = None
        self._scf_task = None
        self._gradient_task = None
        self._system_signature = None
        self._active_geom_signature = None
        self._last_scf_geom_signature = None
        self._last_grad_geom_signature = None
        self._energy = None
        self._gradient = None
        self._scf_results = None
        self._current_scf_mode = None

    def _compute_energy_master(self, molecule):
        self._ensure_system(molecule)

        geom_signature = self._sync_geometry_if_needed(molecule)

        if self._last_scf_geom_signature != geom_signature:

            # SCF warm start: on a reused System Serenity starts the SCF from
            # the electronic structure it already holds, i.e. the orbitals of
            # the previous geometry.  This is independent of (and never
            # implies) an LRSCF restart; see SerenityLinearResponseSolver.
            warm_start = self._system_scf_count > 0
            capture = self.capture_serenity_output('scf')
            try:
                with capture:
                    self.print_title()
                    self._scf_task.run()
            except Exception as error:
                raise SerenityCalculationError(
                    f'Serenity SCF failed: {error}', stage='scf',
                    details={'serenity_output_tail':
                             capture.text[-4000:]}) from error

            energy = float(self._system.getEnergy())
            converged = (SERENITY_SCF_NOT_CONVERGED not in capture.text and
                         np.isfinite(energy))
            if not converged:
                raise SerenityCalculationError(
                    'Serenity SCF did not converge (or returned a nonfinite '
                    'energy); the point is not usable.', stage='scf',
                    details={'energy': energy,
                             'serenity_output_tail': capture.text[-4000:]})

            self._energy = energy
            with self._serenity_output_context():
                ao_basis = self._veloxchem_basis(molecule)
                self._scf_results = self._collect_scf_results(molecule, ao_basis)
            self._last_scf_geom_signature = geom_signature
            self._last_grad_geom_signature = None
            self._gradient = None

            self._scf_revision += 1
            self._system_scf_count += 1
            self._scf_provenance = {
                'scf_revision': int(self._scf_revision),
                'system_name': self._system_name,
                'system_generation': int(self._system_generation),
                'system_reused': bool(warm_start),
                'scf_warm_start': bool(warm_start),
                'geometry_signature': geom_signature,
                'scf_converged': True,
                'reference_energy': float(energy),
                'scf_mode': self._current_scf_mode,
                'rohf_type': self.rohf_type or 'NONE',
                'reference_type': self._reference_type_label(),
                'scf_thresholds': self.get_scf_thresholds(),
                'occupied_space_nesting':
                    self._scf_results.get('occupied_space_nesting'),
            }

        return self._scf_results

    def _reference_type_label(self):
        """Human-readable reference label, e.g. ``CUHF-DFT(bhlyp)``."""

        method = ('HF' if self.method == 'hf' else
                  f'DFT({self.dft_functional})')
        if self._current_scf_mode == 'restricted':
            prefix = 'R'
        elif self.rohf_type not in (None, 'NONE'):
            prefix = self.rohf_type + '-'
        else:
            prefix = 'U'
        return prefix + method

    def get_scf_provenance(self):
        """
        Returns the provenance of the SCF solution currently held.

        :return:
            Dictionary with the SCF revision, System identity, whether the
            SCF started from a previous geometry's orbitals (warm start), the
            reference energy and the reference type, or ``None``.
        """

        if self._scf_provenance is None:
            return None
        return dict(self._scf_provenance)

    def _compute_gradient_master(self, molecule):
        
        self._compute_energy_master(molecule)

        geom_signature = self._active_geom_signature
        if self._last_grad_geom_signature != geom_signature:
            with self._serenity_output_context():
                self._gradient_task.run()
            self._gradient = np.array(self._system.getGeometry().getGradients(),
                                      dtype=float)
            self._last_grad_geom_signature = geom_signature

    # def _collect_scf_results(self):
    #     results = {'energy': float(self._energy)}

    #     mode = self._current_scf_mode
    #     try:
    #         if mode == 'restricted':
    #             es = self._system.getElectronicStructure_R()
    #             results['orbital_energies'] = np.array(es.orbEn(), dtype=float)
    #         else:
    #             es = self._system.getElectronicStructure_U()
    #             results['orbital_energies_alpha'] = np.array(es.alphaOrbEn(),
    #                                                          dtype=float)
    #             results['orbital_energies_beta'] = np.array(es.betaOrbEn(),
    #                                                         dtype=float)
    #     except Exception:
    #         # Keep SCF results minimal if optional details are unavailable.
    #         pass

    #     return results

    def _ensure_system(self, molecule):

        """
        In Sereentiy the systemController 
        """

        mode = self._get_effective_scf_mode(molecule)
        signature = self._get_system_signature(molecule, mode)

        if (self._system is not None and self._system_signature == signature):
            # Reused system: geometry may still have changed (e.g., MD steps).
            self._sync_geometry_if_needed(molecule)
            return

        self._ensure_scratch_dir()

        settings = spy.Settings()
        settings.name = f'vlx_serenity_{uuid.uuid4().hex[:10]}'
        settings.path = self.scratch_dir
        settings.charge = int(molecule.get_charge())
        settings.spin = int(molecule.get_multiplicity() - 1)
        settings.basis.label = self.basis
        settings.basis.densFitJ = self.densfit_j
        settings.grid.accuracy = self.grid_accuracy
        settings.grid.smallGridAccuracy = self.small_grid_accuracy
        settings.scf.maxCycles = self.max_cycles
        thresholds = self.get_scf_thresholds()
        if thresholds['energy'] is not None:
            settings.scf.energyThreshold = thresholds['energy']
        if thresholds['rmsd'] is not None:
            settings.scf.rmsdThreshold = thresholds['rmsd']
        if thresholds['diis'] is not None:
            settings.scf.diisThreshold = thresholds['diis']

        # if mode == 'restricted':
        #     settings.scfMode = spy.SCF_MODES.RESTRICTED
        # else:
        #     settings.scfMode = spy.SCF_MODES.UNRESTRICTED
        if mode == 'restricted':
            settings.scfMode = spy.SCF_MODES.RESTRICTED
            assert_msg_critical(
                self.rohf_type in (None, 'NONE'),
                'SerenityScfDriver: rohf_type requires an unrestricted '
                'open-shell reference.')
        else:
            settings.scfMode = spy.SCF_MODES.UNRESTRICTED
            if self.rohf_type is not None and self.rohf_type != 'NONE':
                assert_msg_critical(
                    settings.spin != 0,
                    'SerenityScfDriver: rohf_type requires multiplicity > 1.')
                settings.scf.rohf = self.rohf_type

        if self.method == 'hf':
            settings.method = spy.ELECTRONIC_STRUCTURE_THEORIES.HF
        else:
            settings.method = spy.ELECTRONIC_STRUCTURE_THEORIES.DFT
            settings.dft.functional = self.dft_functional
        if self.dispersion is not None:
            settings.dft.dispersion = self.dispersion

        labels = list(molecule.get_labels())
        coords = np.array(molecule.get_coordinates_in_bohr(), dtype=float)
        geometry = spy.Geometry(labels, coords)

        with self._serenity_output_context():
            self._system = spy.System(geometry, settings)
            if mode == 'restricted':
                self._scf_task = spy.ScfTask_R(self._system)
                self._gradient_task = spy.GradientTask_R([self._system], [])
            else:
                self._scf_task = spy.ScfTask_U(self._system)
                self._gradient_task = spy.GradientTask_U([self._system], [])

        self._configure_tasks()

        self._system_signature = signature
        self._active_geom_signature = self._get_geometry_signature(molecule)
        self._last_scf_geom_signature = None
        self._last_grad_geom_signature = None
        self._energy = None
        self._gradient = None
        self._scf_results = None
        self._current_scf_mode = mode
        self._system_name = settings.name
        self._system_generation += 1
        self._system_scf_count = 0
        self._scf_provenance = None

    def _sync_geometry_if_needed(self, molecule):
        """
        Synchronizes Serenity system geometry with the current molecule.

        Returns:
            The geometry signature of the input molecule.
        """

        geom_signature = self._get_geometry_signature(molecule)
        
        if self._system is None:
            return geom_signature

        if geom_signature != self._active_geom_signature:
            coords = np.array(molecule.get_coordinates_in_bohr(), dtype=float)
            with self._serenity_output_context():
                self._system.getGeometry().setCoordinates(coords)
            self._active_geom_signature = geom_signature
            self._last_scf_geom_signature = None
            self._last_grad_geom_signature = None
            self._energy = None
            self._gradient = None
            self._scf_results = None

        return geom_signature

    def _configure_tasks(self):
        if hasattr(self._scf_task, 'generalSettings'):
            self._scf_task.generalSettings.printLevel = (
                spy.GLOBAL_PRINT_LEVELS.MINIMUM)

        if hasattr(self._gradient_task, 'generalSettings'):
            self._gradient_task.generalSettings.printLevel = (
                spy.GLOBAL_PRINT_LEVELS.MINIMUM)

        self._gradient_task.settings.gradType = spy.GRADIENT_TYPES.ANALYTICAL

    def _ensure_scratch_dir(self):
        if self.scratch_dir is None:
            # mkdtemp() has no finalizer; without this hook every driver
            # instance leaves a Serenity scratch tree in TMPDIR forever.
            self.scratch_dir = tempfile.mkdtemp(prefix='vlx_serenity_') + os.sep
            self._owns_scratch_dir = True
            self._register_scratch_cleanup()
        elif not self.scratch_dir.endswith(os.sep):
            self.scratch_dir += os.sep

    def _register_scratch_cleanup(self):
        """Arranges for a driver-created scratch root to be removed at exit."""

        if self._scratch_cleanup_hook is not None:
            return

        root = self.scratch_dir
        driver_ref = weakref.ref(self)

        def remove_scratch_root():
            driver = driver_ref()
            if driver is not None and driver.keep_scratch_dir:
                return
            shutil.rmtree(root, ignore_errors=True)

        atexit.register(remove_scratch_root)
        self._scratch_cleanup_hook = remove_scratch_root

    def cleanup_scratch_dir(self):
        """Removes a scratch root this driver created."""

        if self.rank != mpi_master() or not self._owns_scratch_dir:
            return

        shutil.rmtree(self.scratch_dir, ignore_errors=True)
        if self._scratch_cleanup_hook is not None:
            atexit.unregister(self._scratch_cleanup_hook)
            self._scratch_cleanup_hook = None
        self.scratch_dir = None
        self._owns_scratch_dir = False

    def _get_effective_scf_mode(self, molecule=None):
        if self.scf_mode in ('restricted', 'unrestricted'):
            return self.scf_mode
        if molecule is not None and molecule.get_multiplicity() == 1:
            return 'restricted'
        return 'unrestricted'

    def _get_system_signature(self, molecule, mode):
        labels = tuple(molecule.get_labels())
        return (
            mode,
            self.method,
            self.basis.upper(),
            self.dft_functional.lower(),
            self.rohf_type or 'NONE',
            tuple(sorted(self.get_scf_thresholds().items())),
            int(molecule.get_charge()),
            int(molecule.get_multiplicity()),
            labels,
        )

    @staticmethod
    def _get_geometry_signature(molecule):
        coords = np.ascontiguousarray(molecule.get_coordinates_in_bohr(),
                                      dtype=np.float64)
        hasher = hashlib.sha1()
        hasher.update(str(coords.shape).encode('ascii'))
        hasher.update(coords.tobytes())
        return hasher.hexdigest()

    def _serenity_output_context(self):
        if self.serenity_verbose and not self.ostream.is_muted:
            return nullcontext()
        return qc.redirectOutputToFile(os.devnull)

    def capture_serenity_output(self, stage):
        """
        Returns a context manager that captures Serenity's stdout for a task.

        The captured text is kept in ``last_serenity_output[stage]`` so the
        caller can inspect Serenity's convergence warnings.  In verbose mode
        it is replayed to stdout once the task has finished.

        :param stage:
            Label of the task, e.g. ``'scf'``, ``'response'``, ``'gradient'``.
        """

        echo = bool(self.serenity_verbose and not self.ostream.is_muted)
        capture = SerenityOutputCapture(self.scratch_dir, echo=echo)
        driver = self

        class _StageCapture:

            text = ''

            def __enter__(self):
                capture.__enter__()
                return self

            def __exit__(self, exc_type, exc_value, traceback):
                result = capture.__exit__(exc_type, exc_value, traceback)
                self.text = capture.text
                driver.last_serenity_output[stage] = capture.text
                return result

        return _StageCapture()
    
    def get_final_h5py_file(self):
        
        if self.filename is None:
            return None
        else:
            return f"{self.filename}.h5"
        
    def _veloxchem_basis(self, molecule, basis_label=None):

        if basis_label is None:
            vlx_label = self.basis
            if self.basis == "6-31GS":
                vlx_label = "6-31g*"
            return MolecularBasis.read(molecule, vlx_label)
        else:
            vlx_label = basis_label
            if self.basis == "6-31GS":
                vlx_label = "6-31g*"

            return MolecularBasis.read(molecule, vlx_label)
    
    def create_final_hdf5(self, fname, molecule, basis=None):
        if self.rank != mpi_master() or self._skip_writing_h5:
            return False
        if Path(fname).is_file():
            return True

        ao_basis = self._veloxchem_basis(molecule, basis)
        xc_label = self.dft_functional if self.method == 'dft' else 'HF'
        create_hdf5(fname, molecule, ao_basis, xc_label, '')
        return Path(fname).is_file()

    def write_final_hdf5(self, fname, molecule, basis=None):
        if fname is None:
            return False
        if not self.create_final_hdf5(fname, molecule, basis):
            return False

        if self._scf_results is not None:
            write_scf_results_to_hdf5(fname, self._scf_results)

            write_results_to_hdf5(fname, 'serenity', {
                'program': 'Serenity',
                'method': self.method,
                'scf_mode': self._current_scf_mode,
                'basis_label_serenity': self.basis,
            })

        return True

    def _collect_scf_results(self, molecule, ao_basis):
        mode = self._current_scf_mode

        if mode == 'restricted':
            es = self._system.getElectronicStructure_R()
            C_alpha = np.array(es.coeff(), dtype=float)
            C_beta = C_alpha.copy()
            E_alpha = np.array(es.orbEn(), dtype=float)
            E_beta = E_alpha.copy()
            F_alpha = np.array(es.fock(), dtype=float)
            F_beta = F_alpha.copy()
            S = np.array(es.overlap(), dtype=float)

            # Serenity totalDens is total alpha+beta density.
            D_total = np.array(es.totalDens(), dtype=float)
            D_alpha = 0.5 * D_total
            D_beta = 0.5 * D_total

            scf_type = 'restricted'

        else:
            es = self._system.getElectronicStructure_U()
            C_alpha = np.array(es.alphaCoeff(), dtype=float)
            C_beta = np.array(es.betaCoeff(), dtype=float)
            E_alpha = np.array(es.alphaOrbEn(), dtype=float)
            E_beta = np.array(es.betaOrbEn(), dtype=float)
            F_alpha = np.array(es.alphaFock(), dtype=float)
            F_beta = np.array(es.betaFock(), dtype=float)
            D_alpha = np.array(es.alphaDens(), dtype=float)
            D_beta = np.array(es.betaDens(), dtype=float)
            S = np.array(es.overlap(), dtype=float)

            scf_type = 'unrestricted'
        
        perm = self._serenity_to_veloxchem_ao_indices(molecule, ao_basis)

        S = self._reorder_serenity_ao_matrix_to_veloxchem(S, perm)
        D_alpha = self._reorder_serenity_ao_matrix_to_veloxchem(D_alpha, perm)
        D_beta = self._reorder_serenity_ao_matrix_to_veloxchem(D_beta, perm)
        F_alpha = self._reorder_serenity_ao_matrix_to_veloxchem(F_alpha, perm)
        F_beta = self._reorder_serenity_ao_matrix_to_veloxchem(F_beta, perm)

        C_alpha = self._reorder_serenity_mo_coefficients_to_veloxchem(C_alpha, perm)
        C_beta = self._reorder_serenity_mo_coefficients_to_veloxchem(C_beta, perm)

        S_vlx = compute_overlap_integrals(molecule, ao_basis)
        max_overlap_error = np.max(np.abs(S - S_vlx))

        assert_msg_critical(
            np.allclose(S, S_vlx, atol=1.0e-10, rtol=1.0e-10),
            'Serenity AO reordering failed: reordered overlap matrix does not match '
            f'VeloxChem ordering. Max error: {max_overlap_error:.3e}'
        )

        n_mo = C_alpha.shape[1]
        occ_alpha = molecule.get_aufbau_alpha_occupation(n_mo, ao_basis)
        occ_beta = molecule.get_aufbau_beta_occupation(n_mo, ao_basis)

        scf_results = {
            'eri_thresh': 1.0e-12,
            'scf_type': scf_type,
            'rohf_type': self.rohf_type or 'NONE',
            'scf_energy': float(self._energy),
            'restart': False,
            'filename': self.filename,
            'scf_history': [],
            'S': S,
            'C_alpha': C_alpha,
            'C_beta': C_beta,
            'E_alpha': E_alpha,
            'E_beta': E_beta,
            'occ_alpha': occ_alpha,
            'occ_beta': occ_beta,
            'D_alpha': D_alpha,
            'D_beta': D_beta,
            'F_alpha': F_alpha,
            'F_beta': F_beta,
            'F': (F_alpha, F_beta),
        }

        if scf_type == 'unrestricted':
            scf_results['occupied_space_nesting'] = occupied_space_nesting(
                C_alpha, C_beta, S, occ_alpha, occ_beta)

        if self.method == 'dft':
            scf_results['xcfun'] = self.dft_functional.upper()

        try:
            prop = FirstOrderProperties(self.comm, self.ostream)
            prop.compute_scf_prop(molecule, ao_basis, scf_results)
            scf_results['dipole_moment'] = np.array(
                prop.get_property('dipole_moment'))
        except Exception:
            pass

        return scf_results
    
    def _serenity_to_veloxchem_ao_indices(self, molecule, ao_basis):
        serenity_order = []
    
        for atomidx in range(molecule.number_of_atoms()):
            for ang in range(ao_basis.max_angular_momentum() + 1):
                nrad = ao_basis.number_of_basis_functions([atomidx], ang)
                ncomp = 2 * ang + 1

                for radial in range(nrad):
                    for comp in range(ncomp):
                        serenity_order.append((atomidx, ang, radial, comp))

        serenity_index = {
            label: idx for idx, label in enumerate(serenity_order)
        }

        perm = []
        for ang in range(ao_basis.max_angular_momentum() + 1):
            ncomp = 2 * ang + 1

            for comp in range(ncomp):
                for atomidx in range(molecule.number_of_atoms()):
                    nrad = ao_basis.number_of_basis_functions([atomidx], ang)

                    for radial in range(nrad):
                        perm.append(serenity_index[(atomidx, ang, radial, comp)])

        nao = ao_basis.get_dimensions_of_basis()
        assert_msg_critical(
            len(perm) == nao,
            'Serenity AO reordering: inconsistent number of AO indices.'
        )

        return np.array(perm, dtype=int)


    @staticmethod
    def _reorder_serenity_ao_matrix_to_veloxchem(matrix, perm):
        matrix = np.array(matrix, dtype=float)
        return matrix[np.ix_(perm, perm)]


    @staticmethod
    def _reorder_serenity_mo_coefficients_to_veloxchem(coefficients, perm):
        coefficients = np.array(coefficients, dtype=float)
        return coefficients[perm, :]
                