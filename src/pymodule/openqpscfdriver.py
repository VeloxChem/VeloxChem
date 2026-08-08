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
    from oqp.library import set_basis as oqp_set_basis
    from oqp.library import ints_1e as oqp_ints_1e
    from oqp.library.single_point import (SinglePoint, Gradient, LastStep,
                                          BasisOverlap, NACME,
                                          SCFnotConverged)
    from oqp.utils.input_checker import check_input_values
except ImportError:
    pass

# OpenQP signals a failed SCF either by raising (tests.exception = True, which
# every molecule built here sets) or by calling exit().  RuntimeError covers
# the message OpenQP raises when the exception flag is off.
_OPENQP_CONVERGENCE_FAILURES = tuple(
    failure for failure in (
        globals().get('SCFnotConverged'), SystemExit, RuntimeError)
    if failure is not None)


# These arrays are outputs of an overlap/NAC evaluation for one ordered pair
# of geometries.  They are not electronic reference data for the newer
# geometry.  Re-injecting them through ``put_data`` or ``back_door`` on the
# following step makes OpenQP consume a stale pair descriptor as if it
# belonged to the new pair.  Near the ethene S1/S2 near-degeneracy this made
# even an identical-geometry overlap almost rank deficient.
_TRANSIENT_OPENQP_PAIR_KEYS = frozenset({
    'OQP::td_states_overlap',
    'OQP::dc_matrix',
    'OQP::nac_matrix',
})


def _reference_safe_openqp_data(data):
    """Returns a defensive OpenQP snapshot without pair-specific outputs."""

    return {
        key: deepcopy(value)
        for key, value in data.items()
        if key not in _TRANSIENT_OPENQP_PAIR_KEYS
    }


# The converged electronic reference of one geometry.  These are the only
# arrays an SCF initial guess needs: OpenQP rebuilds the one-electron
# integrals, the Fock matrices and the response vectors from the current
# geometry before they are used, so seeding those would inject stale data for
# no benefit.  Alpha/beta orbitals are paired with the density they generate.
_SCF_GUESS_ORBITAL_KEYS = ('OQP::VEC_MO_A', 'OQP::VEC_MO_B')
_SCF_GUESS_DENSITY_KEYS = ('OQP::DM_A', 'OQP::DM_B')
_SCF_GUESS_ELECTRON_TAGS = ('nelec_A', 'nelec_B')
_SCF_GUESS_KEYS = (_SCF_GUESS_ORBITAL_KEYS + _SCF_GUESS_DENSITY_KEYS +
                   ('OQP::E_MO_A', 'OQP::E_MO_B'))


def scf_guess_payload(data):
    """Extracts the geometry-independent electronic reference of a snapshot.

    Accepts either a plain snapshot dictionary or a live ``mol.data`` object;
    tags that the calculation did not produce are skipped.
    """
    
    payload = {}
    for key in _SCF_GUESS_KEYS:
        try:
            value = data[key]
        except (KeyError, AttributeError, TypeError):
            continue
        if value is None:
            continue
        payload[key] = np.array(value, dtype=float)
    return payload


def unpack_triangular(packed, nbf):
    """Expands an OpenQP packed lower triangle into a symmetric matrix."""

    values = np.asarray(packed, dtype=float).ravel()
    if values.size != nbf * (nbf + 1) // 2:
        raise ValueError(
            'unpack_triangular: packed length does not match the requested '
            'dimension.')
    matrix = np.zeros((nbf, nbf))
    rows, cols = np.tril_indices(nbf)
    matrix[rows, cols] = values
    matrix[cols, rows] = values
    return matrix


def pack_triangular(matrix):
    """Packs a symmetric matrix into the OpenQP lower-triangle layout."""

    rows, cols = np.tril_indices(np.asarray(matrix).shape[0])
    return np.asarray(matrix, dtype=float)[rows, cols]


def align_scf_guess_payload(payload, overlap, occupations=None):
    """Re-orthonormalizes seeded orbitals in the AO metric of a new geometry.

    MO coefficients converged at one geometry are orthonormal under *that*
    geometry's overlap matrix only.  Injected verbatim into the next geometry
    they define a density that is not idempotent, so the first SCF energy of
    the new geometry drops *below* its own variational minimum and DIIS has to
    spend several iterations repairing the guess before it can converge.

    A symmetric (Loewdin) transformation ``C <- C (C^T S C)^{-1/2}`` in the new
    metric restores orthonormality while perturbing the orbitals as little as
    possible, and it is applied identically to the alpha and beta sets so an
    ROHF reference keeps sharing one orbital set.  The seeded densities are
    rebuilt from the repaired orbitals, which matters for the restricted
    guess: OpenQP rebuilds the density from the orbitals for ROHF/UHF but uses
    the injected ``OQP::DM_A`` as given for RHF.

    :param payload:
        Guess payload from :func:`scf_guess_payload`; modified in place.
    :param overlap:
        Square AO overlap matrix of the *new* geometry.
    :param occupations:
        Optional ``{orbital key: number of occupied orbitals}`` mapping used to
        rebuild the densities.

    :return:
        The payload and the largest orthonormality violation that was
        repaired, or ``None`` when the payload could not be aligned.
    """

    overlap = np.asarray(overlap, dtype=float)
    nbf = overlap.shape[0]
    occupations = occupations or {}
    deviation = 0.0

    for orbital_key, density_key in zip(_SCF_GUESS_ORBITAL_KEYS,
                                        _SCF_GUESS_DENSITY_KEYS):
        if orbital_key not in payload:
            continue

        stored = np.asarray(payload[orbital_key], dtype=float)
        if stored.size != nbf * nbf:
            # A different basis dimension: the payload cannot describe this
            # calculation at all, so it is left untouched for the caller to
            # reject.
            return payload, None

        # Fortran stores MOs column-wise, so the C-order view has them as rows.
        coefficients = stored.reshape(nbf, nbf).T
        metric = coefficients.T @ overlap @ coefficients
        deviation = max(deviation,
                        float(np.abs(metric - np.eye(nbf)).max()))

        eigenvalues, eigenvectors = np.linalg.eigh(metric)
        if not np.all(np.isfinite(eigenvalues)) or eigenvalues.min() <= 0.0:
            # A singular metric means the stored orbitals do not span the
            # current AO space; a repaired guess would be meaningless.
            return payload, None
        inverse_sqrt = (eigenvectors * eigenvalues**-0.5) @ eigenvectors.T
        coefficients = coefficients @ inverse_sqrt

        payload[orbital_key] = np.ascontiguousarray(coefficients.T)

        n_occupied = occupations.get(orbital_key)
        if (density_key in payload and n_occupied is not None and
                0 <= int(n_occupied) <= nbf):
            occupied = coefficients[:, :int(n_occupied)]
            payload[density_key] = pack_triangular(occupied @ occupied.T)

    return payload, deviation


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

        # Previous-geometry SCF guess.  A relaxed scan or an optimization
        # evaluates a sequence of neighbouring geometries, and restarting every
        # SCF from an extended-Huckel guess throws that continuity away.  The
        # payload is keyed by the electronic problem (not the geometry), so a
        # changed basis, charge, functional or reference type retires it.
        #
        # Unlike the tracking reference of run_oqp_tracked_mrsf, this payload
        # is *not* gated on an accepted nuclear step.  An initial guess cannot
        # change a converged result, only how the SCF gets there, and the
        # orbitals of a rejected trial geometry are just as close to the next
        # trial as the accepted ones.  The tracking reference is gated because
        # it decides which state carries the target identity, which is a
        # correctness question rather than a convergence one.
        self.reuse_scf_guess = True
        self._guess_payload = None
        self._guess_signature = None
        self.last_scf_guess_info = None
        self.scf_guess_reuse_count = 0

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

    def set_reuse_scf_guess(self, enabled=True):
        """Enables or disables the previous-geometry SCF initial guess.

        Disabling it makes every geometry start from the neutral
        extended-Huckel guess again, which is slower but removes any
        dependence of the SCF solution on the order in which geometries were
        visited.
        """

        if self.rank == mpi_master():
            self.reuse_scf_guess = bool(enabled)
            if not self.reuse_scf_guess:
                self.clear_scf_guess()

    def clear_scf_guess(self):
        """Discards the stored previous-geometry SCF guess."""

        self._guess_payload = None
        self._guess_signature = None
        self.last_scf_guess_info = None

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
        guess_signature = self._get_guess_signature(
            molecule, functional, scf_type, multiplicity)
        geom_signature = self._get_geometry_signature(molecule)

        def build_molecule():
            config = self._build_config(
                molecule, method=method, functional=functional,
                scf_type=scf_type, multiplicity=multiplicity,
                exc_multiplicity=exc_multiplicity, tdhf_type=tdhf_type,
                nstate=nstate, grad_state=grad_state)
            return self._new_oqp_molecule(config)

        rebuild = (self._oqp_mol is None or
                   self._calc_signature != calc_signature)

        if rebuild:
            self._oqp_mol = build_molecule()
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
            seeded = self._seed_scf_guess(mol, guess_signature)

            try:
                with self._openqp_output_context(mol):
                    SinglePoint(mol).energy()
            except _OPENQP_CONVERGENCE_FAILURES as exc:
                if not seeded:
                    assert_msg_critical(
                        False,
                        'OpenQPScfDriver: OpenQP calculation did not converge.')
                    raise exc

                # A previous-geometry guess is an accelerator, never a
                # requirement.  If it steers the SCF into a solution it cannot
                # converge, the geometry is retried once from the neutral
                # Huckel guess on a fresh molecule, and the payload is dropped
                # so the failure cannot propagate to later geometries.
                self.clear_scf_guess()
                self.ostream.print_info(
                    'OpenQPScfDriver: the previous-geometry SCF guess did not '
                    'converge; retrying this geometry from the Huckel guess.')
                self.ostream.flush()

                mol = build_molecule()
                self._oqp_mol = mol
                self._calc_signature = calc_signature
                self._active_geom_signature = geom_signature
                self._scf_geom_signature = None
                if grad_state is not None:
                    mol.config['properties']['grad'] = [int(grad_state)]

                with self._openqp_output_context(mol):
                    try:
                        SinglePoint(mol).energy()
                    except SystemExit as retry_exc:
                        assert_msg_critical(
                            False, 'OpenQPScfDriver: OpenQP calculation did '
                            'not converge.')
                        raise retry_exc

            self._capture_scf_guess(mol, guess_signature)

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

    # -- previous-geometry SCF guess --------------------------------------

    def _seed_scf_guess(self, mol, guess_signature, payload=None):
        """Initializes the SCF of this geometry from a previous one.

        The payload is aligned to the AO metric of the *current* geometry
        before it is injected, see :func:`align_scf_guess_payload`.  Seeding is
        skipped -- never fatal -- when no payload is available, when it
        describes a different electronic problem, or when it cannot be aligned.

        :param mol:
            The OpenQP molecule about to run its reference SCF.
        :param guess_signature:
            Signature of the electronic problem, or ``None`` to accept the
            supplied payload without a signature check.
        :param payload:
            Explicit payload; defaults to the driver's stored one.

        :return:
            ``True`` when the guess was injected.
        """

        self.last_scf_guess_info = None

        if not self.reuse_scf_guess:
            return False

        if payload is None:
            if (self._guess_payload is None or guess_signature is None or
                    self._guess_signature != guess_signature):
                return False
            payload = self._guess_payload

        payload = scf_guess_payload(payload)
        if 'OQP::VEC_MO_A' not in payload:
            return False

        # The basis and the one-electron integrals of the current geometry
        # define the metric the orbitals have to be aligned to.  SinglePoint
        # recomputes both in _prep_guess(); repeating them here costs one
        # O(N^2) integral pass and is negligible next to the SCF itself.
        with self._openqp_output_context(mol):
            oqp_set_basis(mol)
            oqp_ints_1e(mol)

        nbf = int(round(np.sqrt(payload['OQP::VEC_MO_A'].size)))
        try:
            overlap = unpack_triangular(mol.data['OQP::SM'], nbf)
        except (AttributeError, KeyError, ValueError):
            # A payload from a different basis dimension cannot seed this SCF.
            return False

        occupations = {}
        for orbital_key, electron_tag in zip(_SCF_GUESS_ORBITAL_KEYS,
                                             _SCF_GUESS_ELECTRON_TAGS):
            try:
                occupations[orbital_key] = int(mol.data[electron_tag])
            except (AttributeError, KeyError, TypeError, ValueError):
                continue

        payload, deviation = align_scf_guess_payload(
            payload, overlap, occupations)
        if deviation is None:
            return False

        # put_data() reports the config blocks that a tag-only snapshot does
        # not carry; that chatter belongs in the OpenQP output channel, not on
        # the VeloxChem console.
        with self._openqp_output_context(mol):
            mol.put_data(payload)

        # check_input_values() rejects guess.type=json without a restart file,
        # so the type is switched after validation.  oqp.guess_json builds the
        # density from the OQP::VEC_MO_* tags just injected, and
        # guess.continue_geom stays False so the previous geometry is not
        # reused as the nuclear configuration.
        mol.config['guess']['type'] = 'json'
        mol.config['guess']['continue_geom'] = False

        self.scf_guess_reuse_count += 1
        self.last_scf_guess_info = {
            'seeded': True,
            'basis_functions': nbf,
            'orthonormality_error': float(deviation),
        }
        return True

    def _capture_scf_guess(self, mol, guess_signature):
        """Stores the converged orbitals of this geometry for the next one.

        Only the MO/density tags are copied out of OpenQP; ``mol.get_data()``
        would materialize every array the calculation holds, including the
        response vectors, which is far more expensive and of no use to an
        initial guess.
        """

        if not self.reuse_scf_guess or guess_signature is None:
            return

        payload = scf_guess_payload(mol.data)
        if 'OQP::VEC_MO_A' not in payload:
            return

        self._guess_payload = payload
        self._guess_signature = guess_signature

    def _get_guess_signature(self, molecule, functional, scf_type,
                             multiplicity):
        """Identity of the electronic problem an SCF guess belongs to.

        The geometry is deliberately absent: carrying orbitals across
        geometries is the whole point.  Everything that changes the meaning or
        the dimension of an MO coefficient is present, so a payload is retired
        automatically when the calculation changes.
        """

        return (
            str(functional).lower(),
            str(scf_type),
            int(multiplicity),
            self.basis.lower(),
            int(molecule.get_charge()),
            tuple(molecule.get_labels()),
        )

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
        # the geometry is deliberately not taken from the payload.  The
        # signature check is bypassed because the tracker payload is by
        # construction the same calculation at the previous geometry.
        had_committed_reference = bool(tracker.has_reference)
        guess_payload = tracker.committed_guess_payload()
        seeded_from_reference = (
            guess_payload is not None and
            self._seed_scf_guess(mol, None, payload=guess_payload))

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

                # A ground-state collision is survivable -- the guard has
                # already kept the identity on the lowest excited root -- but
                # it means the path reached the S1/S0 seam, which a strict run
                # must be told about before a gradient is spent on it.
                if ((not tracking_result.is_reliable or
                     tracking_result.ground_state_collision) and
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
                snapshot_data = _reference_safe_openqp_data(mol.get_data())

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

    def run_oqp_mrsf_states(self,
                            molecule,
                            *,
                            functional,
                            nstate,
                            previous_payload=None,
                            exc_multiplicity=1,
                            multiplicity=3):
        """Runs one MRSF energy/response calculation without a gradient.

        This is the multi-state entry point used by surface hopping.  It
        differs from :meth:`run_oqp_tracked_mrsf` in exactly one way: it makes
        no persistent-state decision.  It returns the *complete* MRSF energy
        array and, when a previous payload is supplied, OpenQP's native state
        overlap in ``(previous, current)`` orientation.  Which raw target the
        dynamics wants is decided afterwards, by the surface-hopping state
        tracker, and only then is a gradient requested through
        :meth:`compute_oqp_state_gradient`.

        The split matters physically: ``BasisOverlap`` must align the current
        and previous MO spaces after the reference SCF but *before* the
        current MRSF response, so the overlap cannot be recovered from two
        finished calculations.

        :param molecule:
            The molecule at the current geometry.
        :param functional:
            Exchange-correlation functional label, empty for Hartree-Fock.
        :param nstate:
            Number of physical MRSF target states.
        :param previous_payload:
            The payload returned by a previous call, or ``None`` at the first
            geometry.  Only a previously *accepted* payload may be passed.
        :param exc_multiplicity:
            Multiplicity of the requested MRSF target manifold.
        :param multiplicity:
            Multiplicity of the high-spin ROHF working reference.

        :return:
            A dictionary with the OpenQP molecule, the full ``energies``
            array, the response energies, the geometry and data snapshot, and
            the native ``overlap`` or ``None``.
        """

        n_states = int(nstate)
        assert_msg_critical(
            n_states >= 1,
            'OpenQPScfDriver.run_oqp_mrsf_states: nstate must be positive.')

        scf_type = 'rohf'
        geom_signature = self._get_geometry_signature(molecule)

        # A fresh OpenQP molecule per geometry.  The native overlap pipeline
        # injects an immutable previous-geometry snapshot through the
        # back-door interface between the reference and the response, and a
        # reused molecule would carry the wrong MO space into that step.  It
        # also isolates the scratch/log identity of every calculation, so a
        # restarted or concurrent run cannot parse another run's output.
        config = self._build_config(
            molecule,
            method='tdhf',
            functional=functional,
            scf_type=scf_type,
            multiplicity=multiplicity,
            exc_multiplicity=exc_multiplicity,
            tdhf_type='mrsf',
            nstate=n_states,
            grad_state=None)
        mol = self._new_oqp_molecule(config)

        self._oqp_mol = mol
        self._calc_signature = self._get_calc_signature(
            molecule, 'tdhf', functional, scf_type, multiplicity, 'mrsf',
            n_states)
        self._active_geom_signature = geom_signature
        self._scf_geom_signature = None

        # Seed the triplet ROHF from the previous accepted orbitals.  Without
        # it the SCF restarts from an extended-Huckel guess at every geometry
        # and can converge to a different ROHF solution, which makes the whole
        # MRSF spectrum jump discontinuously between neighbouring geometries.
        if previous_payload is not None:
            self._seed_scf_guess(mol, None,
                                 payload=previous_payload['data'])

        try:
            with self._openqp_output_context(mol):
                single_point = SinglePoint(mol)
                reference_energy = single_point.reference()

                if previous_payload is not None:
                    mol.config['properties']['back_door'] = True
                    mol.back_door = (deepcopy(previous_payload['xyz']),
                                     deepcopy(previous_payload['data']))
                    BasisOverlap(mol).overlap()

                single_point.excitation(reference_energy)

                overlap = None
                if previous_payload is not None:
                    NACME(mol).nacme()
                    # The tagarray hands back Fortran memory read in C order,
                    # so the delivered matrix is the index transpose of
                    # s_st(previous, current).
                    overlap = native_overlap_to_previous_current(
                        mol.data['OQP::td_states_overlap'])

                # Captured before any gradient runs, exactly as the tracked
                # optimization path does, so the payload describes the aligned
                # response state and not a post-gradient one.
                snapshot_xyz = np.asarray(mol.get_system(), dtype=float).copy()
                snapshot_data = _reference_safe_openqp_data(mol.get_data())

        except SystemExit as exc:
            assert_msg_critical(
                False,
                'OpenQPScfDriver: the MRSF energy/response calculation did '
                'not converge.')
            raise exc

        energies = np.asarray(mol.energies, dtype=float).reshape(-1)
        self._energies = list(energies)
        self._scf_geom_signature = geom_signature

        response_energies = None
        if 'OQP::td_energies' in snapshot_data:
            response_energies = np.asarray(
                snapshot_data['OQP::td_energies'],
                dtype=float).reshape(-1)[:n_states]

        return {
            'molecule': mol,
            'energies': energies,
            'response_energies': response_energies,
            'overlap': overlap,
            'xyz': snapshot_xyz,
            'data': snapshot_data,
            'scf_converged': True,
            'response_converged': True,
        }

    def compute_oqp_state_gradient(self, mol, selector, number_of_atoms):
        """Evaluates the gradient of one MRSF target on an existing molecule.

        The gradient is taken from the *same* OpenQP molecule that produced
        the energies, so a gradient can never belong to a different response
        calculation, geometry or state than its active energy.

        :param mol:
            The OpenQP molecule returned by :meth:`run_oqp_mrsf_states`.
        :param selector:
            OpenQP ``properties.grad`` selector; raw target ``r`` uses
            ``r + 1``.
        :param number_of_atoms:
            Number of atoms, used to reshape the flat OpenQP gradient.

        :return:
            Array of shape (number_of_atoms, 3) in Hartree/bohr.
        """

        state = int(selector)
        assert_msg_critical(
            state >= 1,
            'OpenQPScfDriver.compute_oqp_state_gradient: the OpenQP gradient '
            'selector is one-based; selector 1 is the lowest MRSF target '
            'state, not the working reference.')

        mol.config['properties']['grad'] = [state]

        try:
            with self._openqp_output_context(mol):
                Gradient(mol).gradient()

                if self.d4:
                    d4_step = LastStep(mol)
                    grad_list = mol.config['properties']['grad']
                    d4_energy, d4_gradient = d4_step.get_dispersion(
                        mol, grad_list)
                    d4_step.final_grad(d4_gradient, grad_list)
        except SystemExit as exc:
            assert_msg_critical(
                False,
                'OpenQPScfDriver: the MRSF state gradient did not converge.')
            raise exc

        gradient = np.array(mol.grads[state], dtype=float).reshape(
            (int(number_of_atoms), 3))

        assert_msg_critical(
            np.all(np.isfinite(gradient)),
            'OpenQPScfDriver: OpenQP returned a nonfinite MRSF state '
            'gradient.')

        return gradient

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

        The previous-geometry SCF guess deliberately survives: the
        excited-state gradient driver invalidates this cache after *every*
        optimization step, and dropping the orbitals there would restart each
        geometry from the Huckel guess again.  A settings change cannot leak
        through, because the guess carries its own signature (see
        :meth:`_get_guess_signature`) and is retired when it stops matching.
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
