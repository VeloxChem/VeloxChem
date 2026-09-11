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

from copy import deepcopy

import numpy as np

from .veloxchemlib import mpi_master, hartree_in_ev

from .errorhandler import assert_msg_critical
from .gradientdriver import GradientDriver
from .serenityscfdriver import (SerenityScfDriver, SerenityCalculationError,
                                AdiabaticStateSelectionError,
                                parse_serenity_lr_output)
from .serenitylrrspeigensolver import SerenityLinearResponseSolver
from .transitiondensitytracker import (StateTrackingError,
                                       StateTrackingResult,
                                       StateTrackingStatus,
                                       TransitionDensityTracker)

try:
    from qcserenity import serenipy as spy
except ImportError:
    pass


_MULTIPLICITY_LETTERS = {1: 'S', 2: 'D', 3: 'T', 4: 'Q'}


def manifold_state_label(multiplicity, manifold_state_index):
    """
    Physical label of the n-th state of a spin manifold.

    Singlets count from the ground state (index 1 -> S0, 2 -> S1); other
    manifolds count from 1 (triplet index 1 -> T1).
    """

    multiplicity = int(multiplicity)
    index = int(manifold_state_index)
    if multiplicity == 1:
        return f'S{index - 1}'
    letter = _MULTIPLICITY_LETTERS.get(multiplicity, f'M{multiplicity}_')
    return f'{letter}{index}'


def spin_classification_boundaries(multiplicities):
    """<S^2> values halfway between neighbouring ideal multiplicities."""

    mults = np.asarray(multiplicities, dtype=int).reshape(-1)
    first = 2 if mults.size and int(mults[0]) % 2 == 0 else 1
    candidates = np.arange(first, first + 12, 2)
    spins = 0.5 * (candidates - 1)
    ideal = spins * (spins + 1.0)
    return 0.5 * (ideal[:-1] + ideal[1:])


def classify_spin_quality(state_s2, s2_deviation, multiplicities,
                          ambiguity_margin=0.25, contamination_threshold=0.3):
    """
    Per-root spin-quality diagnostics for approximate SF <S^2> values.

    ``spin_ambiguous`` marks roots whose <S^2> lies within
    ``ambiguity_margin`` of a classification boundary (<S^2> ~ 1 between
    singlet and triplet), where the nearest-multiplicity assignment is not
    meaningful.  ``spin_contaminated`` marks roots further than
    ``contamination_threshold`` from their nearest ideal <S^2>.
    """

    s2 = np.asarray(state_s2, dtype=float).reshape(-1)
    deviation = np.asarray(s2_deviation, dtype=float).reshape(-1)
    boundaries = spin_classification_boundaries(multiplicities)
    distance = np.min(np.abs(s2[:, None] - boundaries[None, :]), axis=1)
    ambiguous = distance < float(ambiguity_margin)
    contaminated = deviation > float(contamination_threshold)
    quality = np.where(ambiguous, 'ambiguous',
                       np.where(contaminated, 'contaminated', 'clean'))
    return {
        'spin_ambiguous': ambiguous,
        'spin_contaminated': contaminated,
        'distance_to_spin_boundary': distance,
        'spin_quality': [str(value) for value in quality],
    }


def select_adiabatic_manifold_root(excitation_energies_ev,
                                   multiplicities,
                                   state_s2,
                                   s2_deviation,
                                   target_multiplicity=1,
                                   manifold_state_index=2,
                                   manifold_filter='nearest',
                                   s2_tolerance=0.5,
                                   spin_ambiguity_margin=0.25,
                                   spin_contamination_threshold=0.3,
                                   near_crossing_threshold_ev=0.15):
    """
    Selects the adiabatic n-th state of a spin manifold from one spectrum.

    A pure function of the CURRENT spectrum -- no previous geometry enters:
    the raw roots are ordered by current excitation energy, the roots of
    ``target_multiplicity`` are kept, and the ``manifold_state_index``-th is
    returned.  With the defaults this is S1, the second-lowest singlet,
    whatever its raw Serenity root number.

    Roots are assigned to a multiplicity by the nearest ideal <S^2>
    (``manifold_filter='nearest'``).  The legacy hard tolerance
    (``'strict'``: |<S^2> - ideal| <= s2_tolerance) is reported for
    comparison only: it moved S0 in and out of the counted manifold when
    its <S^2> crossed 0.5 and so switched the optimization between S1 and
    S2.  Any root at or below the selected energy whose <S^2> lies near a
    classification boundary makes the selection non-robust
    (``selection_robust = False``); it is flagged, never silently trusted.

    :param excitation_energies_ev:
        Excitation energies of the raw roots (eV, any common zero).
    :param multiplicities:
        Nearest-multiplicity classification of every raw root.
    :param state_s2:
        Approximate <S^2> of every raw root.
    :param s2_deviation:
        |<S^2> - ideal <S^2>| of every raw root.
    :param target_multiplicity:
        Multiplicity of the manifold (1 = singlets).
    :param manifold_state_index:
        One-based position in the manifold (2 = S1 for singlets).
    :param manifold_filter:
        ``'nearest'`` or ``'strict'``.
    :param s2_tolerance:
        Tolerance of the strict filter.
    :param spin_ambiguity_margin:
        <S^2> margin around classification boundaries flagged as ambiguous.
    :param spin_contamination_threshold:
        Deviation above which a root is flagged as spin contaminated.
    :param near_crossing_threshold_ev:
        Gap below which a near crossing is flagged.

    :return:
        ``(raw_root, info)`` with a one-based raw root.

    :raises AdiabaticStateSelectionError:
        When the window holds fewer than ``manifold_state_index`` roots of
        the target multiplicity.
    """

    energies = np.asarray(excitation_energies_ev, dtype=float).reshape(-1)
    mults = np.asarray(multiplicities, dtype=int).reshape(-1)
    s2 = np.asarray(state_s2, dtype=float).reshape(-1)
    deviation = np.asarray(s2_deviation, dtype=float).reshape(-1)
    nroots = int(energies.size)
    filter_label = str(manifold_filter).strip().lower()
    target = int(target_multiplicity)
    index = int(manifold_state_index)

    if not (mults.size == s2.size == deviation.size == nroots):
        raise AdiabaticStateSelectionError(
            'energy, multiplicity and <S^2> arrays have different sizes')
    if filter_label not in ('nearest', 'strict'):
        raise ValueError(
            f"manifold_filter must be 'nearest' or 'strict', not "
            f'{manifold_filter!r}')
    if index < 1:
        raise ValueError('manifold_state_index must be >= 1')
    if not (np.all(np.isfinite(energies)) and np.all(np.isfinite(s2))):
        raise AdiabaticStateSelectionError(
            'nonfinite excitation energy or <S^2> in the spectrum')

    order = [int(i) for i in np.argsort(energies, kind='stable')]
    quality = classify_spin_quality(s2, deviation, mults,
                                    spin_ambiguity_margin,
                                    spin_contamination_threshold)
    nearest = [i for i in order if mults[i] == target]
    strict = [i for i in nearest if deviation[i] <= float(s2_tolerance)]
    manifold = nearest if filter_label == 'nearest' else strict
    label = manifold_state_label(target, index)

    info = {
        'target_multiplicity': target,
        'manifold_state_index': index,
        'state_label': label,
        'manifold_filter': filter_label,
        's2_tolerance': float(s2_tolerance),
        'n_roots': nroots,
        'energy_ordered_roots': [i + 1 for i in order],
        'manifold_roots': [i + 1 for i in manifold],
        'manifold_labels': [manifold_state_label(target, k + 1)
                            for k in range(len(manifold))],
        'nearest_multiplicity_roots': [i + 1 for i in nearest],
        'strict_s2_roots': [i + 1 for i in strict],
        'spin_quality': quality['spin_quality'],
        'distance_to_spin_boundary':
            quality['distance_to_spin_boundary'].tolist(),
        'spin_ambiguous_roots': [
            i + 1 for i in range(nroots) if quality['spin_ambiguous'][i]],
        'spin_contaminated_roots': [
            i + 1 for i in range(nroots) if quality['spin_contaminated'][i]],
        'warnings': [],
    }

    if len(manifold) < index:
        info['selected_raw_root'] = None
        raise AdiabaticStateSelectionError(
            f'the computed window of {nroots} root(s) contains only '
            f'{len(manifold)} root(s) of multiplicity {target} '
            f'({filter_label} classification: raw roots '
            f'{[i + 1 for i in manifold]}), but {label} (manifold state '
            f'{index}) was requested; enlarge the response root window',
            details=info)

    chosen = manifold[index - 1]
    position = order.index(chosen)
    ambiguous_below = [i + 1 for i in order[:position + 1]
                       if quality['spin_ambiguous'][i]]
    lower = manifold[index - 2] if index >= 2 else None
    upper = manifold[index] if len(manifold) > index else None
    gap_lower = (None if lower is None else
                 float(energies[chosen] - energies[lower]))
    gap_upper = (None if upper is None else
                 float(energies[upper] - energies[chosen]))
    others = [i for i in range(nroots) if i != chosen]
    closest = (min(others, key=lambda i: abs(energies[i] - energies[chosen]))
               if others else None)
    closest_gap = (None if closest is None else
                   float(abs(energies[closest] - energies[chosen])))
    strict_choice = strict[index - 1] + 1 if len(strict) >= index else None
    threshold = float(near_crossing_threshold_ev)

    info.update({
        'selected_raw_root': chosen + 1,
        'selected_excitation_energy_ev': float(energies[chosen]),
        'selected_multiplicity': int(mults[chosen]),
        'selected_s2': float(s2[chosen]),
        'selected_s2_deviation': float(deviation[chosen]),
        'selected_spin_quality': quality['spin_quality'][chosen],
        'manifold_ground_root': manifold[0] + 1,
        'lower_manifold_root': None if lower is None else lower + 1,
        'upper_manifold_root': None if upper is None else upper + 1,
        'gap_to_lower_manifold_state_ev': gap_lower,
        'gap_to_upper_manifold_state_ev': gap_upper,
        'closest_root': None if closest is None else closest + 1,
        'gap_to_closest_root_ev': closest_gap,
        'near_crossing_threshold_ev': threshold,
        'near_crossing_lower': bool(gap_lower is not None and
                                    gap_lower < threshold),
        'near_crossing_upper': bool(gap_upper is not None and
                                    gap_upper < threshold),
        'near_degenerate_any_root': bool(closest_gap is not None and
                                         closest_gap < threshold),
        'ambiguous_roots_at_or_below_selected': ambiguous_below,
        'selection_robust': not ambiguous_below,
        'strict_filter_selected_root': strict_choice,
        'strict_filter_disagrees': strict_choice != chosen + 1,
    })
    if target == 1 and index == 2:
        info['gap_s0_s1_ev'] = gap_lower
        info['gap_s1_s2_ev'] = gap_upper

    if ambiguous_below:
        info['warnings'].append(
            f'raw root(s) {ambiguous_below} at or below the selected energy '
            f'have <S^2> within {float(spin_ambiguity_margin):.2f} of a '
            f'multiplicity boundary; the {label} assignment depends on '
            'their classification')
    if quality['spin_contaminated'][chosen]:
        info['warnings'].append(
            f'selected raw root {chosen + 1} is spin contaminated '
            f'(<S^2> = {s2[chosen]:.4f}, deviation {deviation[chosen]:.4f})')
    if info['near_crossing_lower'] or info['near_crossing_upper']:
        info['warnings'].append(
            f'near crossing: gap below {gap_lower} eV / above {gap_upper} eV '
            f'(threshold {threshold} eV)')
    return chosen + 1, info


class _PresolvedResponse:
    """
    A solved LR spectrum presented through the ``getLRSCFController()``
    interface of a Serenity GradientTask, so the state-selection hooks can be
    applied before any gradient is requested.
    """

    def __init__(self, controller):
        self._controller = controller

    def getLRSCFController(self):
        return self._controller


class SerenityExcitedStateGradientDriver(GradientDriver):
    """
    Implements Serenity excited-state gradient driver.

    :param serenity_scf_drv:
        The Serenity SCF driver.
    :param serenity_rsp_drv:
        Optional Serenity LR response solver.

    Instance variables
        - state_deriv_index: Excited-state index of interest (1-based).
        - exc_method: Excited-state method (`tda` or `tddft`).
    """

    def __init__(self, serenity_scf_drv, serenity_rsp_drv=None):
        """
        Initializes Serenity excited-state gradient driver.
        """

        errmsg = 'SerenityExcitedStateGradientDriver: invalid Serenity SCF '
        errmsg += 'driver.'
        assert_msg_critical(isinstance(serenity_scf_drv, SerenityScfDriver),
                            errmsg)

        super().__init__(serenity_scf_drv.comm, serenity_scf_drv.ostream)

        self.serenity_driver = serenity_scf_drv
        self.rsp_driver = (serenity_rsp_drv if serenity_rsp_drv is not None else
                           SerenityLinearResponseSolver(serenity_scf_drv))

        self.state_deriv_index = 1
        self.exc_method = 'tda'
        self.enforce_same_multiplicity = True
        self.target_multiplicity = None
        self.s2_tolerance = 0.5

        # State selection.  'adiabatic' optimizes the manifold_state_index-th
        # root of target_multiplicity ordered by CURRENT energy (S1 = second
        # singlet) and uses transition-density tracking only as a diagnostic
        # guard.  None keeps the historical behaviour: overlap tracking when
        # enforce_same_multiplicity is set or a tracker is attached, the raw
        # root otherwise.  See set_adiabatic_state().
        self.state_selection_mode = None
        self.manifold_state_index = None
        self.manifold_filter = 'nearest'
        self.spin_ambiguity_margin = 0.25
        self.spin_contamination_threshold = 0.3
        self.near_crossing_threshold_ev = 0.15
        self.max_root_window_expansions = 2
        self.root_window_increment = 5
        self.expand_root_window_on_low_overlap = True
        self.root_identity_tolerance_ev = 1.0e-3
        self.root_identity_min_overlap = 0.99

        # SCF warm start: keep the Serenity System (and hence the previous
        # orbitals as SCF guess) across geometries.  Off by default; see
        # SerenityScfDriver._invalidate_cache.  Never implies an LR restart.
        self.reuse_scf_system = False

        self.state_selection_info = None
        self.evaluation_record = None
        self.evaluation_history = []
        self.last_gradient_task_provenance = None
        self._guard_initialized_this_evaluation = False

        self.excited_state_energy = None
        self.selected_excitation_energy = None
        self.total_energy = None
        self.computes_energy_with_gradient = True

        # The high-spin working-reference energy, taken directly from
        # Serenity's ``system.getEnergy()``.  It is recorded explicitly rather
        # than back-solved as ``total_energy - selected_excitation_energy``,
        # because the driver invalidates the Serenity system cache after every
        # ``compute`` and the reconstruction would then be the only remaining
        # source.  The high-spin determinant is a working reference: it is
        # never an SF dynamics surface.
        self.reference_energy = None

        self.reference_s2 = None
        self.delta_s2 = None
        self.state_s2 = None
        self.state_multiplicities = None
        self.ideal_state_s2 = None
        self.s2_deviation = None
        self.multiplicity_valid_mask = None
        self.selected_s2 = None
        self.selected_multiplicity = None
        self.selected_s2_deviation = None

        self.flag = 'Serenity Excited-State Gradient Driver'

        self._grad_task = None
        self.state_tracker = None
        self.tracking_info = None
        self.last_lrscf_controller = None
        self._last_tracking_system = None
        self._last_tracking_mode = None
        self._last_tracking_molecule = None
        self._tracking_applied_in_compute = False
        # OptimizationEngine enables this while bounded LOW_OVERLAP recovery
        # is still available.  Non-strict fallbacks are applied only after
        # those retries have been exhausted.
        self._tracking_recovery_active = False
        self._current_gradient_task_roots = []
        self.tracking_history = []
        self._tracking_evaluation_counter = 0

        self._input_keywords['gradient'].update({
            'enforce_same_multiplicity':
                ('bool', 'restrict spin-flip root tracking by multiplicity'),
            'target_multiplicity':
                ('int', 'target spin multiplicity; omit to infer initially'),
            's2_tolerance':
                ('float', 'maximum deviation from ideal spin squared'),
            'state_selection_mode':
                ('str_lower', 'raw, tracked or adiabatic state selection'),
            'manifold_state_index':
                ('int', 'one-based state index in the target spin manifold'),
            'manifold_filter':
                ('str_lower', 'nearest or strict multiplicity classification'),
            'spin_ambiguity_margin':
                ('float', '<S^2> margin flagged ambiguous at boundaries'),
            'near_crossing_threshold_ev':
                ('float', 'gap below which a near crossing is flagged'),
            'reuse_scf_system':
                ('bool', 'reuse the Serenity System across geometries'),
        })


    def set_state_deriv_index(self, state_deriv_index):
        """
        Sets the excited-state index of interest (1-based).
        """

        state = int(state_deriv_index)
        assert_msg_critical(
            state > 0,
            'SerenityExcitedStateGradientDriver: state index must be > 0')
        self.state_deriv_index = state

    def set_exc_method(self, exc_method):
        """
        Sets excited-state method (`tda` or `tddft`).
        """

        method = str(exc_method).strip().lower()
        if method in ('tda', 'cis'):
            self.exc_method = 'tda'
        elif method in ('tddft', 'rpa'):
            self.exc_method = 'tddft'
        else:
            errmsg = 'SerenityExcitedStateGradientDriver: invalid exc_method '
            errmsg += f'"{exc_method}"'
            assert_msg_critical(False, errmsg)

        self.rsp_driver.set_exc_method(self.exc_method)

    def update_settings(self,
                        grad_dict,
                        rsp_dict=None,
                        method_dict=None):
        """
        Updates settings in excited-state gradient driver.

        :param grad_dict:
            The input dictionary of gradient settings group.
        :param rsp_dict:
            The input dictionary of response settings group.
        :param method_dict:
            The input dictionary of method settings group.
        """

        if grad_dict is None:
            grad_dict = {}

        if rsp_dict is None:
            rsp_dict = {}

        if method_dict is None:
            method_dict = {}

        # parse common gradient settings: numerical, do_four_point, delta_h
        super().update_settings(grad_dict, method_dict={})

        if 'state_deriv_index' in grad_dict:
            state_val = grad_dict['state_deriv_index']
            if isinstance(state_val, (list, tuple, np.ndarray)):
                assert_msg_critical(
                    len(state_val) > 0,
                    'SerenityExcitedStateGradientDriver: empty state list')
                state_val = state_val[0]
            self.set_state_deriv_index(state_val)

        if 'exc_method' in grad_dict:
            self.set_exc_method(grad_dict['exc_method'])

        if 'method' in rsp_dict:
            self.set_exc_method(rsp_dict['method'])

        if 'tamm_dancoff' in rsp_dict:
            self.set_exc_method('tda' if bool(rsp_dict['tamm_dancoff']) else
                                'tddft')

        self.rsp_driver.update_settings(rsp_dict, method_dict)

        assert_msg_critical(
            self.target_multiplicity is None or self.target_multiplicity > 0,
            'SerenityExcitedStateGradientDriver: target_multiplicity must be '
            'positive.')
        assert_msg_critical(
            self.s2_tolerance >= 0.0,
            'SerenityExcitedStateGradientDriver: s2_tolerance must be '
            'non-negative.')
        self._validate_selection_settings()

    def compute(self, molecule):
        """
        Performs calculation of Serenity excited-state gradient.

        :param molecule:
            The molecule.
        """

        self._tracking_applied_in_compute = False
        self._current_gradient_task_roots = []
        self.evaluation_record = None
        self._guard_initialized_this_evaluation = False

        tracking_error = None
        calculation_error = None
        try:
            if self.numerical:
                if self.rank == mpi_master():
                    self.compute_numerical(molecule)
                else:
                    self.gradient = None
            else:
                if self.rank == mpi_master():
                    self.gradient = self._compute_analytical_master(molecule)
                else:
                    self.gradient = None
        except StateTrackingError as error:
            tracking_error = (str(error), error.result)
        except SerenityCalculationError as error:
            # Fail closed: a failed SCF, LR solve, state selection or
            # gradient never yields an energy or gradient.
            self._reset_serenity_state_after_failure()
            if self.comm.Get_size() == 1:
                raise
            calculation_error = (str(error), error.stage, error.details)

        tracking_error = self.comm.bcast(
            tracking_error, root=mpi_master())
        if tracking_error is not None:
            message, result = tracking_error
            raise StateTrackingError(message, result)
        if self.comm.Get_size() > 1:
            calculation_error = self.comm.bcast(
                calculation_error, root=mpi_master())
            if calculation_error is not None:
                raise SerenityCalculationError(*calculation_error)

        self.gradient = self.comm.bcast(self.gradient, root=mpi_master())
        if self.rank == mpi_master():
            state_payload = {
                'state_deriv_index': int(self.state_deriv_index),
                'target_multiplicity': self.target_multiplicity,
                'excited_state_energy': self.excited_state_energy,
                'selected_excitation_energy':
                    self.selected_excitation_energy,
                'total_energy': self.total_energy,
                'reference_energy': self.reference_energy,
                'reference_s2': self.reference_s2,
                'delta_s2': self.delta_s2,
                'state_s2': self.state_s2,
                'state_multiplicities': self.state_multiplicities,
                'ideal_state_s2': self.ideal_state_s2,
                's2_deviation': self.s2_deviation,
                'multiplicity_valid_mask':
                    self.multiplicity_valid_mask,
                'selected_s2': self.selected_s2,
                'selected_multiplicity': self.selected_multiplicity,
                'selected_s2_deviation': self.selected_s2_deviation,
                'tracking_info': self.tracking_info,
                'tracking_applied': self._tracking_applied_in_compute,
                'state_selection_info': self.state_selection_info,
                'evaluation_record': self.evaluation_record,
            }
        else:
            state_payload = None

        state_payload = self.comm.bcast(
            state_payload, root=mpi_master())
        self.state_deriv_index = int(state_payload['state_deriv_index'])
        self.target_multiplicity = state_payload['target_multiplicity']
        self.excited_state_energy = state_payload['excited_state_energy']
        self.selected_excitation_energy = (
            state_payload['selected_excitation_energy'])
        self.total_energy = state_payload['total_energy']
        self.reference_energy = state_payload['reference_energy']
        self.reference_s2 = state_payload['reference_s2']
        self.delta_s2 = state_payload['delta_s2']
        self.state_s2 = state_payload['state_s2']
        self.state_multiplicities = (
            state_payload['state_multiplicities'])
        self.ideal_state_s2 = state_payload['ideal_state_s2']
        self.s2_deviation = state_payload['s2_deviation']
        self.multiplicity_valid_mask = (
            state_payload['multiplicity_valid_mask'])
        self.selected_s2 = state_payload['selected_s2']
        self.selected_multiplicity = state_payload['selected_multiplicity']
        self.selected_s2_deviation = (
            state_payload['selected_s2_deviation'])
        self.tracking_info = state_payload['tracking_info']
        self._tracking_applied_in_compute = bool(
            state_payload['tracking_applied'])
        self.state_selection_info = state_payload['state_selection_info']
        self.evaluation_record = state_payload['evaluation_record']
        if (self._tracking_applied_in_compute and
                self.tracking_info is not None):
            self.tracking_history.append(deepcopy(self.tracking_info))
        if self.evaluation_record is not None:
            self.evaluation_history.append(
                self._evaluation_summary(self.evaluation_record))

        self.print_geometry(molecule)
        self.print_gradient(molecule, [self.state_deriv_index])

        if not self.reuse_scf_system:
            # Fresh Serenity System for the next geometry: fresh SCF guess,
            # no LR restart file.  With reuse_scf_system the System (and its
            # orbitals, i.e. the SCF warm start) is kept; an LR restart is
            # still refused across geometries by decide_lr_restart().
            self.serenity_driver._invalidate_cache()
            self.rsp_driver._invalidate_rsp_cache()

        self.ostream.print_blank()
        self.ostream.flush()

    def compute_energy(self, molecule):
        """
        Computes excited-state total energy at current geometry.

        :param molecule:
            The molecule.

        :return:
            The excited-state total energy.
        """

        if self.rank != mpi_master():
            return None

        adiabatic = self._effective_selection_mode() == 'adiabatic'
        if adiabatic:
            # Energy-only evaluations (numerical gradients) follow the same
            # adiabatic definition; the tracking guard is not run for them.
            spectrum, selection, _ = (
                self._adiabatic_spectrum_selection_and_guard(
                    molecule, None, with_guard=False))
            rsp_results = spectrum['results']
            self.state_deriv_index = int(selection['selected_raw_root'])
            self.state_selection_info = selection
        else:
            rsp_results = self.rsp_driver.compute(molecule, broadcast=False)
        eigenvalues = np.asarray(
            rsp_results['eigenvalues'], dtype=float).reshape(-1)

        if self.rsp_driver.spinflip:
            self._set_spin_metadata(rsp_results)
            if self.enforce_same_multiplicity and not adiabatic:
                self.state_deriv_index = (
                    self._select_energy_root_by_multiplicity(
                        self.state_deriv_index))
            self._update_selected_spin_diagnostics()

        self.excited_state_energy = eigenvalues
        self.selected_excitation_energy = self._extract_excitation_energy(
            rsp_results)
        self.reference_energy = float(self.serenity_driver.get_energy())
        self.total_energy = float(self.reference_energy +
                                  self.selected_excitation_energy)

        return self.total_energy

    def _compute_analytical_master(self, molecule):
        """
        Energy and gradient of the selected state at the current geometry.

        The response spectrum is solved and inspected BEFORE any gradient is
        requested, so the gradient task runs once, for the selected root:

          1. SCF at the current geometry (fail-closed);
          2. LR spectrum at the current geometry, solved in the current MO
             basis -- a stored LR solution seeds the Davidson solver only if
             it belongs to this geometry, System and SCF solution;
          3. convergence/completeness checks and spin diagnostics;
          4. selection: ``adiabatic`` (n-th root of a spin manifold by current
             energy), overlap ``tracked``, or ``raw``;
          5. comparison with the last accepted state (a diagnostic guard in
             adiabatic mode);
          6. excited-state gradient for the selected raw root; its LRSCF step
             may restart from step 2, which was solved at this geometry;
          7. verification that the gradient task's root is the selected one.
        """

        # 1. SCF (raises SerenityCalculationError when not converged).
        self.serenity_driver._compute_energy_master(molecule)

        mode = self.serenity_driver._current_scf_mode
        requested_state = int(self.state_deriv_index)
        selection_mode = self._effective_selection_mode()
        self._current_gradient_task_roots = []
        self.state_selection_info = None
        tracking_applied = False
        selection = None
        guard = None

        if selection_mode == 'adiabatic':
            # 2.-5. Spectrum, adiabatic selection, tracking guard.
            spectrum, selection, guard = (
                self._adiabatic_spectrum_selection_and_guard(molecule, mode))
            selected_state = int(selection['selected_raw_root'])
            tracking_applied = guard is not None
        else:
            # 2. Spectrum at the current geometry; the historical selection
            # hooks receive it through the GradientTask-like adapter.
            spectrum = self._solve_response_spectrum(
                molecule, minimum_roots=requested_state)
            presolved = _PresolvedResponse(spectrum['controller'])
            if self.rsp_driver.spinflip:
                self._update_spin_metadata(presolved.getLRSCFController())

            if selection_mode == 'raw':
                selected_state = requested_state
            elif self.rsp_driver.spinflip:
                if self.enforce_same_multiplicity:
                    selected_state = self._select_same_multiplicity_state(
                        molecule, presolved, mode, requested_state)
                    tracking_applied = self.state_tracker is not None
                elif self.state_tracker is not None:
                    selected_state = self._select_tracked_state(
                        presolved, mode, requested_state)
                    tracking_applied = True
                else:
                    selected_state = requested_state
            elif self.state_tracker is not None:
                # Ordinary TDA/TDDFT has no interleaved spin manifolds, but its
                # raw roots can still exchange character. Select the
                # continuation from the transition-density overlap before the
                # gradient is requested, without any spin filter.
                selected_state = self._select_tracked_state(
                    presolved, mode, requested_state)
                tracking_applied = True
            else:
                selected_state = requested_state

        # 6. One gradient task, for the selected raw root only.
        self.state_deriv_index = int(selected_state)
        grad_task = self._run_excited_gradient_task(mode)

        # 7. The gradient task re-solves the LR problem; check that its root
        # with this index is the root that was selected.
        verification = self._verify_gradient_root_identity(
            grad_task, spectrum, int(selected_state))

        self.state_deriv_index = int(selected_state)
        if self.rsp_driver.spinflip:
            self._update_selected_spin_diagnostics()
            self.ostream.print_info(
                'Spin-flip state selection: '
                f'root {self.state_deriv_index}, '
                f'{self._format_selected_spin_label()}, '
                f'<S^2> = {self.selected_s2:.6f}, '
                f'deviation = {self.selected_s2_deviation:.6f}.')
            self.ostream.flush()

        controller = grad_task.getLRSCFController()
        excitation_energies = np.asarray(
            controller.getExcitationEnergies('isolated'),
            dtype=float).reshape(-1) / hartree_in_ev()

        assert_msg_critical(
            self.state_deriv_index <= excitation_energies.size,
            'SerenityExcitedStateGradientDriver: selected state index '
            f'{self.state_deriv_index} but only {excitation_energies.size} '
            'state(s) are available.')

        self.excited_state_energy = excitation_energies
        self.selected_excitation_energy = float(
            excitation_energies[self.state_deriv_index - 1])
        self.reference_energy = float(self.serenity_driver.get_energy())
        self.total_energy = float(self.reference_energy +
                                  self.selected_excitation_energy)
        gradient = np.array(
            self.serenity_driver._system.getGeometry().getGradients(),
            dtype=float)
        if not (np.isfinite(self.total_energy) and
                np.all(np.isfinite(gradient))):
            raise SerenityCalculationError(
                'Serenity returned a nonfinite excited-state energy or '
                'gradient.', stage='gradient',
                details={'root': int(self.state_deriv_index)})

        self._grad_task = grad_task
        self.last_lrscf_controller = controller
        self._last_tracking_system = self.serenity_driver._system
        self._last_tracking_mode = mode
        self._last_tracking_molecule = molecule

        if selection_mode == 'adiabatic':
            self._finalize_adiabatic_evaluation(
                spectrum, selection, guard, verification, gradient, mode)
            return gradient

        if tracking_applied:
            if (self.tracking_info is not None and
                    self.tracking_info.get(
                        'reference_update_recommended', False) and
                    not self.tracking_info.get('initialized', False)):
                # The selected-root gradient and its controller snapshot are
                # staged together. geomeTRIC's accepted-step hook commits both;
                # a rejected trial discards the proposal.
                self.state_tracker.propose_reference(
                    self._last_tracking_system,
                    self.last_lrscf_controller,
                    self._last_tracking_mode,
                    self.state_deriv_index,
                )
                self.tracking_info['reference_staged'] = True
            self._tracking_applied_in_compute = True

        if self.tracking_info is not None:
            tracking_update = {
                'step': int(len(self.tracking_history)),
                'evaluation_id': int(self._tracking_evaluation_counter),
                'gradient_root': int(self.state_deriv_index),
                'gradient_task_roots':
                    list(self._current_gradient_task_roots),
                'gradient_recomputed':
                    len(self._current_gradient_task_roots) > 1,
                'selected_excitation_energy':
                    float(self.selected_excitation_energy),
                'total_energy': float(self.total_energy),
                'gradient_rms': float(np.sqrt(np.mean(gradient**2))),
                'state_excitation_energies':
                    np.asarray(self.excited_state_energy, dtype=float).copy(),
                'state_excitation_energies_hartree':
                    np.asarray(self.excited_state_energy, dtype=float).copy(),
            }
            if self.rsp_driver.spinflip:
                tracking_update.update({
                    'state_multiplicities': np.asarray(
                        self.state_multiplicities, dtype=int).copy(),
                    'state_s2': np.asarray(
                        self.state_s2, dtype=float).copy(),
                    's2_deviation': np.asarray(
                        self.s2_deviation, dtype=float).copy(),
                })
            self.tracking_info.update(tracking_update)
            self._tracking_evaluation_counter += 1
            self._print_tracking_diagnostics(self.tracking_info)

        self.evaluation_record = self._build_evaluation_record(
            selection_mode, spectrum, None, None, verification, gradient)
        return gradient

    # ------------------------------------------------------------------
    # State-selection configuration
    # ------------------------------------------------------------------

    def set_adiabatic_state(self, target_multiplicity=1,
                            manifold_state_index=2,
                            manifold_filter='nearest'):
        """
        Optimizes the n-th state of a spin manifold defined by current energy.

        With the defaults the driver follows the adiabatic S1 surface: at
        every geometry all roots are classified by nearest multiplicity, the
        singlets are ordered by their current energy and the gradient of the
        second one is computed.  Raw Serenity root numbers may change freely;
        an attached transition-density tracker only reports on changes of
        state character and never overrides the energy ordering.

        :param target_multiplicity:
            Spin multiplicity of the manifold (1 = singlets).
        :param manifold_state_index:
            One-based position in the manifold (1 = S0, 2 = S1, ...).
        :param manifold_filter:
            ``'nearest'`` (recommended) or ``'strict'`` (legacy hard
            <S^2> tolerance, for comparisons only).
        """

        self.state_selection_mode = 'adiabatic'
        self.target_multiplicity = int(target_multiplicity)
        self.manifold_state_index = int(manifold_state_index)
        self.manifold_filter = str(manifold_filter).strip().lower()
        self._validate_selection_settings()

    def _validate_selection_settings(self):
        mode = getattr(self, 'state_selection_mode', None)
        assert_msg_critical(
            mode in (None, 'raw', 'tracked', 'adiabatic'),
            'SerenityExcitedStateGradientDriver: state_selection_mode must '
            'be raw, tracked or adiabatic.')
        assert_msg_critical(
            getattr(self, 'manifold_filter', 'nearest') in ('nearest',
                                                            'strict'),
            'SerenityExcitedStateGradientDriver: manifold_filter must be '
            'nearest or strict.')
        if mode == 'adiabatic':
            # Adiabatic S1 unless configured otherwise.
            if self.target_multiplicity is None:
                self.target_multiplicity = 1
            if self.manifold_state_index is None:
                self.manifold_state_index = 2
            assert_msg_critical(
                int(self.manifold_state_index) >= 1 and
                int(self.target_multiplicity) >= 1,
                'SerenityExcitedStateGradientDriver: adiabatic selection '
                'needs target_multiplicity >= 1 and manifold_state_index '
                '>= 1.')

    def _effective_selection_mode(self):
        """``'adiabatic'``, ``'tracked'`` or ``'raw'``."""

        mode = getattr(self, 'state_selection_mode', None)
        if mode is not None:
            return str(mode)
        spinflip = bool(getattr(self.rsp_driver, 'spinflip', False))
        if (getattr(self, 'state_tracker', None) is not None or
                (spinflip and getattr(self, 'enforce_same_multiplicity',
                                      False))):
            return 'tracked'
        return 'raw'

    # ------------------------------------------------------------------
    # Response spectrum and adiabatic selection
    # ------------------------------------------------------------------

    def _solve_response_spectrum(self, molecule, minimum_roots=None):
        """
        Solves (or reuses) the LR spectrum at the current geometry.

        Serenity's spin-flip gradient is always SF-TDA, so the spectrum used
        for the selection is solved with SF-TDA too; otherwise the raw root
        numbers of the two solves need not refer to the same states.

        :return:
            Dictionary with ``results``, ``controller``, ``energies_ev``,
            ``vectors`` (first component, copied before the gradient task
            runs) and ``provenance``.
        """

        rsp = self.rsp_driver
        if rsp.spinflip and rsp.exc_method != 'tda':
            self.ostream.print_info(
                'Serenity spin-flip gradients are SF-TDA; the spectrum used '
                'for state selection is solved with SF-TDA as well.')
            self.ostream.flush()
            rsp.set_exc_method('tda')
            self.exc_method = 'tda'
        if minimum_roots is not None and int(minimum_roots) > int(rsp.nstates):
            rsp.set_nstates(int(minimum_roots))

        results = rsp.compute(molecule, broadcast=False)
        controller = rsp.get_lr_controller()
        energies_ev = (np.asarray(results['eigenvalues'],
                                  dtype=float).reshape(-1) * hartree_in_ev())
        return {
            'results': results,
            'controller': controller,
            'energies_ev': energies_ev,
            'vectors': self._first_vector_component(controller),
            'provenance': results.get('lr_provenance'),
        }

    @staticmethod
    def _first_vector_component(controller):
        """Copy of the first excitation-vector component, or ``None``."""

        try:
            vectors = controller.getExcitationVectors('isolated')
            first = np.array(vectors[0], dtype=float, copy=True)
        except Exception:
            return None
        return first if first.ndim == 2 else None

    def _adiabatic_spectrum_selection_and_guard(self, molecule, mode,
                                                with_guard=True):
        """
        Solves the spectrum, selects the adiabatic state and runs the guard.

        The root window is enlarged (at most ``max_root_window_expansions``
        times) when the target manifold is incomplete, and once when the
        previous state's character overlaps with no current root, which can
        mean that it has left the window.

        :return:
            ``(spectrum, selection, guard)``; ``guard`` is ``None`` without a
            tracker or with ``with_guard=False``.
        """

        assert_msg_critical(
            bool(self.rsp_driver.spinflip),
            'SerenityExcitedStateGradientDriver: adiabatic manifold selection '
            'requires a spin-flip response, whose spectrum contains the '
            'ground state and a multiplicity for every root.')
        self._validate_selection_settings()

        expansions = []
        low_overlap_retry = False
        while True:
            spectrum = self._solve_response_spectrum(molecule)
            self._set_spin_metadata(spectrum['results'])
            try:
                selection = self._select_adiabatic_root(spectrum['energies_ev'])
            except AdiabaticStateSelectionError as error:
                if len(expansions) >= int(self.max_root_window_expansions):
                    error.details['root_window_expansions'] = expansions
                    raise
                expansions.append(self._expand_root_window(
                    'target manifold incomplete'))
                continue

            guard = (self._adiabatic_tracking_guard(spectrum, mode, selection)
                     if with_guard else None)
            if (guard is not None and guard.get('low_overlap') and
                    self.expand_root_window_on_low_overlap and
                    not low_overlap_retry and
                    len(expansions) < int(self.max_root_window_expansions)):
                low_overlap_retry = True
                expansions.append(self._expand_root_window(
                    'the previous state character overlaps with no root in '
                    'the window'))
                continue
            break

        selection['root_window_expansions'] = expansions
        selection['nstates'] = int(self.rsp_driver.nstates)
        selection['reference_s2'] = (None if self.reference_s2 is None else
                                     float(self.reference_s2))
        if guard is not None and guard.get('low_overlap') and low_overlap_retry:
            guard['warnings'].append(
                'low overlap persists after enlarging the root window; the '
                'adiabatic selection is kept')
        self.state_selection_info = selection
        return spectrum, selection, guard

    def _select_adiabatic_root(self, energies_ev):
        """Applies select_adiabatic_manifold_root to the current metadata."""

        _, info = select_adiabatic_manifold_root(
            energies_ev,
            self.state_multiplicities,
            self.state_s2,
            self.s2_deviation,
            target_multiplicity=int(self.target_multiplicity),
            manifold_state_index=int(self.manifold_state_index),
            manifold_filter=self.manifold_filter,
            s2_tolerance=float(self.s2_tolerance),
            spin_ambiguity_margin=float(self.spin_ambiguity_margin),
            spin_contamination_threshold=float(
                self.spin_contamination_threshold),
            near_crossing_threshold_ev=float(self.near_crossing_threshold_ev))
        return info

    def _expand_root_window(self, reason):
        """Enlarges the response root window by ``root_window_increment``."""

        old = int(self.rsp_driver.nstates)
        new = old + int(self.root_window_increment)
        self.rsp_driver.set_nstates(new)
        self.ostream.print_info(
            f'Serenity adiabatic selection: enlarging the response window '
            f'from {old} to {new} roots ({reason}).')
        self.ostream.flush()
        return {'from_nstates': old, 'to_nstates': new, 'reason': str(reason)}

    def _adiabatic_tracking_guard(self, spectrum, mode, selection):
        """
        Compares the current spectrum with the last accepted adiabatic state.

        Diagnostic only; the result never changes the selected root.  It
        reports where the previous state's character went (largest
        transition-density overlap, ground state included), how much of it
        the current adiabatic root carries, whether a character crossing is
        suggested and whether continuity is lost (low overlap).
        """

        tracker = getattr(self, 'state_tracker', None)
        if tracker is None:
            return None

        controller = spectrum['controller']
        selected = int(selection['selected_raw_root'])
        nroots = int(np.asarray(spectrum['energies_ev']).size)
        in_manifold = np.zeros(nroots, dtype=bool)
        for root in selection['manifold_roots']:
            in_manifold[int(root) - 1] = True
        metadata = self._tracking_candidate_metadata(controller, in_manifold)
        previous_root = (int(tracker.reference_state)
                         if tracker.has_reference() else None)
        label = selection['state_label']

        guard = {
            'state_label': label,
            'previous_selected_root': previous_root,
            'current_adiabatic_root': selected,
            'reference_staged': False,
            'reference_updated': False,
            'warnings': [],
        }
        try:
            result = tracker.track(
                self.serenity_driver._system, controller, mode,
                active_reference_state=selected, allowed_states=None,
                candidate_metadata=metadata)
        except Exception as error:
            guard.update({
                'status': 'UNAVAILABLE', 'initialized': False,
                'guard_available': False, 'low_overlap': False,
                'character_crossing_suspected': False,
                'tracked_matches_adiabatic': None, 'overlap_matrix': None,
                'continuity': 'unavailable',
            })
            guard['warnings'].append(
                f'transition-density tracking failed: {error}')
            return guard

        info = result.to_dict()
        initialized = bool(info.get('initialized', False))
        if initialized:
            self._guard_initialized_this_evaluation = True
        overlap = info.get('overlap_matrix')
        guard.update({
            'status': info.get('status'),
            'initialized': initialized,
            'reference_state': info.get('reference_state'),
            'tracked_root': info.get('new_state'),
            'max_overlap': info.get('max_overlap'),
            'second_state': info.get('second_state'),
            'second_overlap': info.get('second_overlap'),
            'overlap_ratio': info.get('overlap_ratio'),
            'ground_state_root': info.get('ground_state_root'),
            'ground_state_collision': bool(
                info.get('ground_state_collision', False)),
            'global_state': info.get('global_state'),
            'overlap_matrix': overlap,
            'candidate_table': info.get('candidate_table'),
        })
        guard['warnings'].extend(info.get('warnings') or [])

        if initialized:
            guard.update({
                'guard_available': True, 'tracked_root_any': selected,
                'tracked_overlap_any': 1.0,
                'overlap_with_adiabatic_root': 1.0,
                'tracked_matches_adiabatic': True,
                'character_crossing_suspected': False, 'low_overlap': False,
                'continuity': 'reference_initialized',
            })
            return guard

        reference = guard.get('reference_state')
        matrix = None if overlap is None else np.asarray(overlap)
        if (matrix is None or matrix.ndim != 2 or reference is None or
                not 1 <= int(reference) <= matrix.shape[1] or
                matrix.shape[0] != nroots):
            guard.update({
                'guard_available': False, 'low_overlap': False,
                'character_crossing_suspected': False,
                'tracked_matches_adiabatic': None,
                'continuity': 'unavailable',
            })
            guard['warnings'].append(
                'no transition-density overlap is available for this '
                'evaluation')
            return guard

        column = np.abs(np.real(matrix[:, int(reference) - 1])).astype(float)
        scores = np.where(np.isfinite(column), column, -np.inf)
        any_root = int(np.argmax(scores)) + 1
        any_overlap = float(column[any_root - 1])
        adiabatic_overlap = float(column[selected - 1])
        minimum = float(getattr(tracker, 'min_overlap', 0.5))
        low = (not np.isfinite(any_overlap)) or any_overlap < minimum
        matches = any_root == selected
        crossing = (not low) and (not matches)
        guard.update({
            'guard_available': True,
            'tracked_root_any': any_root,
            'tracked_overlap_any': any_overlap,
            'overlap_with_adiabatic_root': adiabatic_overlap,
            'tracked_matches_adiabatic': bool(matches),
            'character_crossing_suspected': bool(crossing),
            'low_overlap': bool(low),
            'continuity': ('continuous' if matches and not low else
                           'character_crossing' if crossing else
                           'low_overlap'),
        })
        if crossing:
            guard['warnings'].append(
                f'the character of the previous {label} (raw root '
                f'{reference}) now overlaps most ({any_overlap:.3f}) with raw '
                f'root {any_root}, while the adiabatic {label} is raw root '
                f'{selected} (overlap {adiabatic_overlap:.3f}); the adiabatic '
                'definition is kept (possible physical crossing)')
        if low:
            guard['warnings'].append(
                f'largest overlap with the previous {label} character is '
                f'{any_overlap:.3f} < {minimum:.3f}; continuity is not '
                'established')
        return guard

    def _finalize_adiabatic_evaluation(self, spectrum, selection, guard,
                                       verification, gradient, mode):
        """Stages the tracking reference and writes the evaluation records."""

        selected = int(selection['selected_raw_root'])
        if guard is not None:
            if not guard.get('initialized', False):
                # In adiabatic mode the accepted-step reference is always the
                # current adiabatic state; the overlap classification does not
                # decide what is staged.  geomeTRIC's accepted-step hook
                # commits it and a rejected trial rolls it back, so every
                # comparison is against the last ACCEPTED geometry.
                self.state_tracker.propose_reference(
                    self._last_tracking_system, spectrum['controller'], mode,
                    selected)
                guard['reference_staged'] = True
            self._tracking_applied_in_compute = True

        record = self._build_evaluation_record(
            'adiabatic', spectrum, selection, guard, verification, gradient)
        self.evaluation_record = record
        self.state_selection_info = selection

        g = guard or {}
        warnings = list(selection.get('warnings') or [])
        warnings.extend(g.get('warnings') or [])
        self.tracking_info = {
            'tracking_framework': 'serenity_spinflip',
            'state_selection_mode': 'adiabatic',
            'overlap_source': ('serenity_transition_density_overlap'
                               if guard is not None else None),
            'status': 'ADIABATIC' if guard is None else g.get('status'),
            'initialized': bool(g.get('initialized', False)),
            'state_label': selection['state_label'],
            'selected_raw_root': selected,
            'selected_state': selected,
            'new_state': selected,
            'old_state': g.get('previous_selected_root'),
            'assignment_confident': bool(selection['selection_robust'] and
                                         verification['verified']),
            'selection_robust': bool(selection['selection_robust']),
            'tracked_root': g.get('tracked_root'),
            'tracked_root_any': g.get('tracked_root_any'),
            'tracked_overlap_any': g.get('tracked_overlap_any'),
            'overlap_with_adiabatic_root': g.get('overlap_with_adiabatic_root'),
            'tracked_matches_adiabatic': g.get('tracked_matches_adiabatic'),
            'character_crossing_suspected':
                g.get('character_crossing_suspected'),
            'low_overlap': g.get('low_overlap'),
            'max_overlap': g.get('max_overlap'),
            'second_overlap': g.get('second_overlap'),
            'reference_staged': bool(g.get('reference_staged', False)),
            'reference_updated': False,
            'overlap_matrix': g.get('overlap_matrix'),
            'gradient_root': selected,
            'gradient_task_roots': list(self._current_gradient_task_roots),
            'gradient_recomputed': False,
            'selected_excitation_energy': float(self.selected_excitation_energy),
            'total_energy': float(self.total_energy),
            'gradient_rms': float(np.sqrt(np.mean(gradient**2))),
            'state_excitation_energies_hartree': np.asarray(
                self.excited_state_energy, dtype=float).copy(),
            'state_multiplicities': np.asarray(
                self.state_multiplicities, dtype=int).copy(),
            'state_s2': np.asarray(self.state_s2, dtype=float).copy(),
            's2_deviation': np.asarray(self.s2_deviation, dtype=float).copy(),
            'gap_to_lower_manifold_state_ev':
                selection.get('gap_to_lower_manifold_state_ev'),
            'gap_to_upper_manifold_state_ev':
                selection.get('gap_to_upper_manifold_state_ev'),
            'near_crossing_lower': selection.get('near_crossing_lower'),
            'near_crossing_upper': selection.get('near_crossing_upper'),
            'warnings': warnings,
        }
        self._tracking_evaluation_counter += 1
        self._print_adiabatic_selection(record)

    # ------------------------------------------------------------------
    # Gradient task
    # ------------------------------------------------------------------

    def _gradient_lr_signature(self):
        """Response-settings signature of the gradient task's LRSCF step."""

        rsp = self.rsp_driver
        nstates = max(int(self.state_deriv_index), int(rsp.nstates))
        method = 'tda' if rsp.spinflip else self.exc_method
        return rsp._get_rsp_signature(nstates=nstates, exc_method=method)

    def _run_excited_gradient_task(self, mode, lr_restart=None):
        """
        Runs one Serenity gradient task for the current raw response root.

        The task re-solves the LR problem.  Its Serenity restart flag is set
        only when the System's stored LR solution was obtained at this
        geometry, with this SCF solution and the same response settings
        (normally the spectrum just solved for the state selection); after
        any geometry change it is False.

        :param mode:
            ``'restricted'`` or ``'unrestricted'``.
        :param lr_restart:
            Optional explicit restart flag; ``None`` applies the policy of
            ``SerenityLinearResponseSolver.decide_lr_restart``.

        :return:
            The finished Serenity gradient task.
        """

        self._current_gradient_task_roots.append(
            int(self.state_deriv_index))

        signature = self._gradient_lr_signature()
        if lr_restart is None:
            decision = self.rsp_driver.decide_lr_restart(signature)
        else:
            decision = {
                'restart': bool(lr_restart),
                'reason': 'explicitly requested by the caller',
                'policy': getattr(self.rsp_driver, 'lr_restart_policy', None),
            }

        if mode == 'restricted':
            grad_task = spy.GradientTask_R(self.serenity_driver._system)
        else:
            grad_task = spy.GradientTask_U(self.serenity_driver._system)

        self._configure_excited_gradient_task(
            grad_task, lr_restart=decision['restart'])

        capture = self.serenity_driver.capture_serenity_output('gradient')
        try:
            with capture:
                grad_task.run()
        except Exception as error:
            self.rsp_driver.invalidate_lr_restart_ledger()
            raise SerenityCalculationError(
                f'Serenity excited-state gradient task failed: {error}',
                stage='gradient',
                details={'root': int(self.state_deriv_index),
                         'scf_provenance': getattr(
                             self.serenity_driver, 'get_scf_provenance',
                             lambda: None)(),
                         'serenity_output_tail': capture.text[-4000:]}
            ) from error

        # The captured text covers every iterative solve of the task (the
        # LRSCF step and the Z-vector equations).
        parsed = parse_serenity_lr_output(capture.text)
        provenance = {
            'root': int(self.state_deriv_index),
            'lr_restart_requested': bool(decision['restart']),
            # None if Serenity printed no restart message.
            'lr_restart_used': (parsed['restart_loaded']
                                if decision['restart'] else False),
            'lr_restart_reason': decision['reason'],
            'lr_converged': parsed['converged'],
            'davidson_iterations': parsed['davidson_iterations'],
            'n_converged_solves': parsed['n_converged_solves'],
            'warnings': parsed['warnings'],
        }
        self.last_gradient_task_provenance = provenance
        if parsed['converged'] is not True:
            self.rsp_driver.invalidate_lr_restart_ledger()
            if parsed['converged'] is None:
                reason = ('Serenity printed no convergence status for the '
                          'iterative solves of the gradient task, so '
                          'convergence cannot be verified')
            else:
                reason = ('An iterative solve of the Serenity gradient task '
                          '(LRSCF step or Z-vector) did not converge '
                          '("Convergence criterion not reached")')
            raise SerenityCalculationError(
                f'{reason}; the gradient is not usable.', stage='gradient',
                details=provenance)

        # The gradient task left its converged LR solution in the System.
        self.rsp_driver.record_lr_solution(signature)
        return grad_task

    def _verify_gradient_root_identity(self, grad_task, spectrum,
                                       selected_state):
        """
        Checks that the gradient task's root is the selected root.

        Compares the excitation energy (and, when available, the excitation
        vector) of ``selected_state`` between the spectrum used for the
        selection and the LRSCF step of the gradient task.

        :raises SerenityCalculationError:
            When the two roots differ.
        """

        controller = grad_task.getLRSCFController()
        grad_energies = np.asarray(
            controller.getExcitationEnergies('isolated'),
            dtype=float).reshape(-1)
        spectrum_energies = np.asarray(spectrum['energies_ev'],
                                       dtype=float).reshape(-1)
        index = int(selected_state) - 1
        info = {
            'requested_root': int(selected_state),
            'n_roots_spectrum': int(spectrum_energies.size),
            'n_roots_gradient_task': int(grad_energies.size),
            'selected_energy_difference_ev': None,
            'max_energy_difference_ev': None,
            'selected_vector_overlap': None,
        }
        problems = []
        if not (0 <= index < grad_energies.size and
                index < spectrum_energies.size):
            problems.append('the selected root is outside the gradient-task '
                            'spectrum')
        else:
            difference = float(abs(grad_energies[index] -
                                   spectrum_energies[index]))
            count = min(grad_energies.size, spectrum_energies.size)
            info['selected_energy_difference_ev'] = difference
            info['max_energy_difference_ev'] = float(np.max(np.abs(
                grad_energies[:count] - spectrum_energies[:count])))
            tolerance = float(getattr(self, 'root_identity_tolerance_ev',
                                      1.0e-3))
            minimum_overlap = float(getattr(self, 'root_identity_min_overlap',
                                            0.99))
            if not np.isfinite(difference) or difference > tolerance:
                problems.append(
                    f'the excitation energy of raw root {selected_state} '
                    f'differs by {difference:.3e} eV between the selection '
                    'spectrum and the gradient task')

            before = spectrum.get('vectors')
            after = self._first_vector_component(controller)
            if (before is not None and after is not None and
                    before.shape[0] == after.shape[0] and
                    index < before.shape[1] and index < after.shape[1]):
                a = before[:, index]
                b = after[:, index]
                norm = float(np.linalg.norm(a) * np.linalg.norm(b))
                overlap = float(abs(a @ b) / norm) if norm > 0.0 else 0.0
                info['selected_vector_overlap'] = overlap
                if overlap < minimum_overlap:
                    problems.append(
                        f'the excitation vector of raw root {selected_state} '
                        f'has overlap {overlap:.6f} with the selected root')
        info['verified'] = not problems
        info['problems'] = problems
        if problems:
            raise SerenityCalculationError(
                'The Serenity gradient task did not solve for the selected '
                'root: ' + '; '.join(problems), stage='root_identity',
                details=info)
        return info

    def _reset_serenity_state_after_failure(self):
        """Drops the System, the LR ledger and a just-created reference."""

        for action in (getattr(self.serenity_driver, '_invalidate_cache',
                               None),
                       getattr(self.rsp_driver, '_invalidate_rsp_cache',
                               None)):
            if action is not None:
                try:
                    action()
                except Exception:
                    pass
        tracker = getattr(self, 'state_tracker', None)
        if (tracker is not None and
                getattr(self, '_guard_initialized_this_evaluation', False)):
            clear = getattr(tracker, 'clear_reference', None)
            if clear is not None:
                clear()
        self._guard_initialized_this_evaluation = False

    # ------------------------------------------------------------------
    # Records
    # ------------------------------------------------------------------

    def _build_evaluation_record(self, selection_mode, spectrum, selection,
                                 guard, verification, gradient):
        """
        Everything needed to reconstruct the state selection of one
        evaluation: SCF provenance, LR provenance (restart requested/used),
        full spectrum with spin diagnostics, the selection, the tracking
        guard and the gradient task with its root verification.
        """

        results = spectrum.get('results') or {}
        energies = np.asarray(results.get('eigenvalues', []),
                              dtype=float).reshape(-1)
        reference_energy = self.reference_energy
        gradient = np.asarray(gradient, dtype=float)
        scf_provenance = getattr(self.serenity_driver, 'get_scf_provenance',
                                 None)

        def as_list(value, dtype=float):
            if value is None:
                return None
            return np.asarray(value, dtype=dtype).reshape(-1).tolist()

        tracking = None
        if guard is not None:
            tracking = {key: value for key, value in guard.items()
                        if key not in ('overlap_matrix', 'candidate_table')}

        return {
            'evaluation_id': len(getattr(self, 'evaluation_history', []) or
                                 []),
            'state_selection_mode': selection_mode,
            'state_label': (None if selection is None else
                            selection.get('state_label')),
            'geometry_signature': getattr(self.serenity_driver,
                                          '_active_geom_signature', None),
            'scf': scf_provenance() if callable(scf_provenance) else None,
            'reference_s2': getattr(self, 'reference_s2', None),
            'response': spectrum.get('provenance'),
            'spectrum': {
                'excitation_energies_hartree': energies.tolist(),
                'excitation_energies_ev':
                    (energies * hartree_in_ev()).tolist(),
                'total_energies_hartree': (
                    None if reference_energy is None else
                    (float(reference_energy) + energies).tolist()),
                'delta_s2': as_list(getattr(self, 'delta_s2', None)),
                'state_s2': as_list(getattr(self, 'state_s2', None)),
                'multiplicities': as_list(
                    getattr(self, 'state_multiplicities', None), int),
                's2_deviation': as_list(getattr(self, 's2_deviation', None)),
            },
            'selection': (selection if selection is not None else {
                'selected_raw_root': int(self.state_deriv_index),
                'state_selection_mode': selection_mode,
            }),
            'tracking': tracking,
            'gradient': {
                'requested_root': int(self.state_deriv_index),
                'gradient_task_roots':
                    list(self._current_gradient_task_roots),
                'task': getattr(self, 'last_gradient_task_provenance', None),
                'root_identity': verification,
                'gradient_rms': float(np.sqrt(np.mean(gradient**2))),
                'gradient_max': (float(np.max(np.linalg.norm(
                    gradient.reshape(-1, 3), axis=1))) if gradient.size
                                 else None),
                'finite': bool(np.all(np.isfinite(gradient))),
            },
            'reference_energy': reference_energy,
            'selected_excitation_energy_hartree':
                self.selected_excitation_energy,
            'total_energy': self.total_energy,
        }

    @staticmethod
    def _evaluation_summary(record):
        """Compact per-evaluation summary kept in ``evaluation_history``."""

        selection = record.get('selection') or {}
        tracking = record.get('tracking') or {}
        response = record.get('response') or {}
        scf = record.get('scf') or {}
        task = (record.get('gradient') or {}).get('task') or {}
        return {
            'evaluation_id': record.get('evaluation_id'),
            'geometry_signature': record.get('geometry_signature'),
            'state_label': record.get('state_label'),
            'selected_raw_root': selection.get('selected_raw_root'),
            'manifold_roots': selection.get('manifold_roots'),
            'total_energy': record.get('total_energy'),
            'selected_excitation_energy_ev':
                selection.get('selected_excitation_energy_ev'),
            'selected_s2': selection.get('selected_s2'),
            'selection_robust': selection.get('selection_robust'),
            'strict_filter_selected_root':
                selection.get('strict_filter_selected_root'),
            'gap_to_lower_manifold_state_ev':
                selection.get('gap_to_lower_manifold_state_ev'),
            'gap_to_upper_manifold_state_ev':
                selection.get('gap_to_upper_manifold_state_ev'),
            'near_crossing_lower': selection.get('near_crossing_lower'),
            'near_crossing_upper': selection.get('near_crossing_upper'),
            'tracked_root_any': tracking.get('tracked_root_any'),
            'tracked_overlap_any': tracking.get('tracked_overlap_any'),
            'overlap_with_adiabatic_root':
                tracking.get('overlap_with_adiabatic_root'),
            'character_crossing_suspected':
                tracking.get('character_crossing_suspected'),
            'low_overlap': tracking.get('low_overlap'),
            'tracking_status': tracking.get('status'),
            'lr_restart_requested': response.get('restart_requested'),
            'lr_restart_requested': response.get('restart_requested'),
            'lr_restart_used': response.get('restart_used'),
            'lr_converged': response.get('converged'),
            'gradient_lr_restart_requested': task.get('lr_restart_requested'),
            'gradient_lr_restart_used': task.get('lr_restart_used'),
            'gradient_lr_converged': task.get('lr_converged'),
            'scf_warm_start': scf.get('scf_warm_start'),
            'reference_s2': record.get('reference_s2'),
            'reference_nesting_error':
                (scf.get('occupied_space_nesting') or {}).get('nesting_error'),
        }

    def _print_adiabatic_selection(self, record):
        """Prints the adiabatic selection table of one evaluation."""

        selection = record['selection']
        tracking = record.get('tracking') or {}
        response = record.get('response') or {}
        gradient = record.get('gradient') or {}
        task = gradient.get('task') or {}
        identity = gradient.get('root_identity') or {}
        scf = record.get('scf') or {}
        spectrum = record['spectrum']
        label = selection['state_label']
        info = self.ostream.print_info

        def ev(value):
            return 'n/a' if value is None else f'{value:.4f} eV'

        self.ostream.print_header(
            'Serenity Spin-Flip Adiabatic State Selection')
        info(f"Rule                 : {label} = state "
             f"{selection['manifold_state_index']} of multiplicity "
             f"{selection['target_multiplicity']} by current energy "
             f"({selection['manifold_filter']} classification)")
        info(f"Selected raw root    : {selection['selected_raw_root']} "
             f"({selection['selected_excitation_energy_ev']:.6f} eV, "
             f"<S^2> = {selection['selected_s2']:.4f}, "
             f"{selection['selected_spin_quality']})")
        info('Manifold             : ' + ', '.join(
            f'{name}=raw {root}' for name, root in
            zip(selection['manifold_labels'], selection['manifold_roots'])))
        info(f"Gap lower / upper    : "
             f"{ev(selection['gap_to_lower_manifold_state_ev'])} / "
             f"{ev(selection['gap_to_upper_manifold_state_ev'])}; closest "
             f"root {selection['closest_root']} at "
             f"{ev(selection['gap_to_closest_root_ev'])}")
        info(f"Selection robust     : {selection['selection_robust']}; "
             f"strict <S^2> filter would pick raw root "
             f"{selection['strict_filter_selected_root']}")
        info(f"Reference            : {scf.get('reference_type')}, "
             f"<S^2>_ref = {record.get('reference_s2')}, SCF warm start = "
             f"{scf.get('scf_warm_start')}")
        nesting = scf.get('occupied_space_nesting') or {}
        if nesting.get('nesting_error') is not None:
            info(f"Nesting              : alpha/beta occupied spaces nested to "
                 f"{nesting['nesting_error']:.2e} (within ROHF gradient "
                 f"tolerance: {nesting['within_serenity_rohf_tolerance']}); "
                 f"SCF thresholds {scf.get('scf_thresholds')}")
        info(f"LR spectrum          : restart requested "
             f"{response.get('restart_requested')}, used "
             f"{response.get('restart_used')} "
             f"({response.get('restart_reason')}); converged "
             f"{response.get('converged')} in "
             f"{response.get('davidson_iterations')} iterations")
        info(f"Gradient task        : root {task.get('root')}, LR restart "
             f"requested {task.get('lr_restart_requested')}, used "
             f"{task.get('lr_restart_used')} "
             f"({task.get('lr_restart_reason')}); converged "
             f"{task.get('lr_converged')}; root verified "
             f"{identity.get('verified')} (dE = "
             f"{identity.get('selected_energy_difference_ev')} eV, overlap "
             f"{identity.get('selected_vector_overlap')})")
        if tracking:
            info(f"Tracking guard       : {tracking.get('status')} / "
                 f"{tracking.get('continuity')}; previous {label} raw root "
                 f"{tracking.get('previous_selected_root')} -> max overlap "
                 f"raw root {tracking.get('tracked_root_any')} "
                 f"({tracking.get('tracked_overlap_any')}), adiabatic root "
                 f"overlap {tracking.get('overlap_with_adiabatic_root')}")
        for message in (list(selection.get('warnings') or []) +
                        list(tracking.get('warnings') or [])):
            info(f'WARNING: {message}')

        energies = spectrum['excitation_energies_ev']
        spins = spectrum['state_s2'] or [float('nan')] * len(energies)
        mults = spectrum['multiplicities'] or [0] * len(energies)
        names = dict(zip(selection['manifold_roots'],
                         selection['manifold_labels']))
        info(' raw   excitation/eV     <S^2>  mult  spin quality  manifold')
        for root in selection['energy_ordered_roots']:
            k = root - 1
            marker = (f'  <== {label}' if root == selection['selected_raw_root']
                      else '')
            info(f'{root:4d} {energies[k]:16.6f} {spins[k]:9.4f} '
                 f'{mults[k]:5d}  {selection["spin_quality"][k]:>12s}  '
                 f'{names.get(root, "-"):>8s}{marker}')
        self.ostream.print_blank()
        self.ostream.flush()

    def _update_spin_metadata(self, controller):
        metadata = self.rsp_driver.get_spinflip_metadata(controller)
        self._set_spin_metadata(metadata)

    def _set_spin_metadata(self, metadata):
        self.reference_s2 = metadata['reference_s2']
        self.delta_s2 = np.asarray(metadata['delta_s2'], dtype=float)
        self.state_s2 = np.asarray(metadata['state_s2'], dtype=float)
        self.state_multiplicities = np.asarray(
            metadata['state_multiplicities'], dtype=int)
        self.ideal_state_s2 = np.asarray(
            metadata['ideal_state_s2'], dtype=float)
        self.s2_deviation = np.asarray(
            metadata['s2_deviation'], dtype=float)

    def _update_selected_spin_diagnostics(self):
        selected_index = int(self.state_deriv_index) - 1
        assert_msg_critical(
            0 <= selected_index < self.state_s2.size,
            'SerenityExcitedStateGradientDriver: selected spin-flip root is '
            'outside the spin metadata range.')

        self.selected_s2 = float(self.state_s2[selected_index])
        self.selected_multiplicity = int(
            self.state_multiplicities[selected_index])
        self.selected_s2_deviation = float(
            self.s2_deviation[selected_index])

    def _format_selected_spin_label(self):
        """
        Formats the spin assignment of the selected spin-flip root.

        Nearest-multiplicity classification always returns some multiplicity,
        so an SF root with ``<S^2>`` near 1 -- equidistant from the singlet
        (0) and triplet (2) values -- is reported as a clean singlet or
        triplet depending only on which side of the tie it falls.  Those roots
        are spin-incomplete rather than spin-pure, and printing a bare
        "multiplicity 3" for ``<S^2> = 1.02`` claims a triplet that is not
        there.  Qualify the label whenever the root sits further from its
        nearest physical ``<S^2>`` than the configured tolerance allows.
        """

        if self.selected_s2_deviation > self.s2_tolerance:
            return ('spin-incomplete (nearest multiplicity '
                    f'{self.selected_multiplicity})')

        return f'multiplicity {self.selected_multiplicity}'

    def _select_energy_root_by_multiplicity(self, requested_state):
        """
        Validates a raw root and records spin diagnostics for energy-only use.

        Serenity SF-TDA roots are not guaranteed spin-pure, so the approximate
        ``<S^2>`` classification is diagnostic here rather than a hard mask.
        """

        requested_state = int(requested_state)
        nstates = int(self.state_s2.size)
        assert_msg_critical(
            1 <= requested_state <= nstates,
            'SerenityExcitedStateGradientDriver: requested state index is '
            'outside the spin-flip spectrum.')

        if self.target_multiplicity is None:
            self.target_multiplicity = int(
                self.state_multiplicities[requested_state - 1])

        self.multiplicity_valid_mask = (
            (self.state_multiplicities == int(self.target_multiplicity)) &
            (self.s2_deviation <= float(self.s2_tolerance))
        )
        return requested_state

    def _select_same_multiplicity_state(self, molecule, grad_task, mode,
                                        requested_state):
        """
        Tracks all valid roots, then applies the configured failure policy.

        Approximate ``<S^2>`` never removes a finite response root before the
        overlap comparison.  A multiplicity mismatch is reported separately as
        ``SPIN_CONFLICT``.
        """

        del molecule

        nstates = int(np.asarray(self.state_s2).size)
        assert_msg_critical(
            1 <= requested_state <= nstates,
            'SerenityExcitedStateGradientDriver: requested state index is '
            'outside the spin-flip spectrum.')

        if self.target_multiplicity is None:
            self.target_multiplicity = int(
                self.state_multiplicities[requested_state - 1])

        same_inferred_multiplicity = (
            (np.asarray(self.state_multiplicities, dtype=int) ==
             int(self.target_multiplicity))
        )
        within_s2_tolerance = (
            np.asarray(self.s2_deviation, dtype=float) <=
            float(self.s2_tolerance)
        )
        self.multiplicity_valid_mask = (
            same_inferred_multiplicity & within_s2_tolerance)

        if self.state_tracker is None:
            self.state_tracker = TransitionDensityTracker(
                self.rsp_driver,
                target_state=requested_state,
            )

        controller = grad_task.getLRSCFController()
        # Spin compatibility must combine the inferred multiplicity with the
        # <S^2> tolerance.  Nearest-multiplicity classification alone labels a
        # spin-incomplete SF root (<S^2> ~ 1, equidistant from 0 and 2) as a
        # singlet, so following it would be reported as CONFIDENT instead of
        # SPIN_CONFLICT.
        candidate_metadata = self._tracking_candidate_metadata(
            controller, self.multiplicity_valid_mask)
        result = self.state_tracker.track(
            self.serenity_driver._system,
            controller,
            mode,
            active_reference_state=requested_state,
            # Only invalid response vectors are hard-excluded by the tracker.
            allowed_states=None,
            candidate_metadata=candidate_metadata,
        )

        if (result['initialized'] and
                not bool(self.multiplicity_valid_mask[requested_state - 1])):
            diagnostics = result.to_dict()
            diagnostics.pop('status', None)
            diagnostics['assignment_confident'] = False
            diagnostics['warnings'] = list(diagnostics['warnings']) + [
                'initial raw root conflicts with target multiplicity'
            ]
            result = StateTrackingResult(
                StateTrackingStatus.SPIN_CONFLICT, diagnostics)

        selected_state, result = self._apply_tracking_failure_policy(result)
        selected_index = selected_state - 1
        info = result.to_dict()
        info.update({
            'tracking_framework': 'serenity_spinflip',
            'overlap_source': 'serenity_transition_density_overlap',
            'multiplicity_filter_applied': True,
            'target_multiplicity': int(self.target_multiplicity),
            'selected_state': int(selected_state),
            'selected_multiplicity':
                int(self.state_multiplicities[selected_index]),
            'selected_s2': float(self.state_s2[selected_index]),
            'selected_s2_deviation':
                float(self.s2_deviation[selected_index]),
            'spin_diagnostic_roots': [
                int(state) for state in
                np.flatnonzero(same_inferred_multiplicity) + 1
            ],
            's2_tolerance_roots': [
                int(state) for state in
                np.flatnonzero(self.multiplicity_valid_mask) + 1
            ],
            'reference_updated': False,
            'reference_staged': False,
        })
        self.tracking_info = info

        return int(selected_state)

    def _select_tracked_state(self, grad_task, mode, requested_state):
        """Selects a response root without multiplicity filtering.

        Ordinary TDA/TDDFT always uses this path. Spin-flip can also use it when
        ``enforce_same_multiplicity`` is disabled. No approximate ``<S^2>``
        classification or multiplicity mask is constructed here.
        """

        controller = grad_task.getLRSCFController()
        energies = np.asarray(
            controller.getExcitationEnergies('isolated'),
            dtype=float).reshape(-1)
        assert_msg_critical(
            1 <= int(requested_state) <= energies.size,
            'SerenityExcitedStateGradientDriver: requested state index is '
            'outside the ordinary TDA/TDDFT spectrum.')

        candidate_metadata = self._tracking_candidate_metadata(
            controller, nstates=energies.size)
        result = self.state_tracker.track(
            self.serenity_driver._system,
            controller,
            mode,
            active_reference_state=int(requested_state),
            allowed_states=None,
            candidate_metadata=candidate_metadata,
        )
        selected_state, result = self._apply_tracking_failure_policy(result)
        info = result.to_dict()
        framework = ('serenity_spinflip' if self.rsp_driver.spinflip else
                     'serenity_tddft')
        info.update({
            'tracking_framework': framework,
            'overlap_source': 'serenity_transition_density_overlap',
            'multiplicity_filter_applied': False,
            'selected_state': int(selected_state),
            'reference_updated': False,
            'reference_staged': False,
        })
        self.tracking_info = info
        return int(selected_state)

    def _tracking_candidate_metadata(self, controller, spin_compatible=None,
                                     nstates=None):
        """Builds per-root diagnostics, including optional solver metadata."""

        if nstates is None:
            if spin_compatible is not None:
                nstates = int(np.asarray(spin_compatible).size)
            else:
                nstates = int(np.asarray(
                    controller.getExcitationEnergies('isolated')).size)
        residuals = self._optional_controller_vector(
            controller,
            ('getResidualNorms', 'getResiduals', 'getEigenpairResiduals'),
            nstates,
        )
        converged = self._optional_controller_vector(
            controller,
            ('getConvergedRoots', 'getRootConvergence'),
            nstates,
        )
        metadata = {
            'solver_residual': residuals,
            'solver_converged': converged,
        }
        if spin_compatible is not None:
            metadata.update({
                's2': np.asarray(self.state_s2, dtype=float),
                'inferred_multiplicity': np.asarray(
                    self.state_multiplicities, dtype=int),
                's2_deviation': np.asarray(self.s2_deviation, dtype=float),
                'spin_compatible': np.asarray(spin_compatible, dtype=bool),
            })
        return metadata

    @staticmethod
    def _optional_controller_vector(controller, method_names, nstates):
        for method_name in method_names:
            method = getattr(controller, method_name, None)
            if method is None:
                continue
            try:
                values = np.asarray(method(), dtype=float).reshape(-1)
            except Exception:
                continue
            if values.size == nstates:
                return values
        return None

    def _apply_tracking_failure_policy(self, result):
        status = result.status
        selected_state = result.get('new_state')
        if status is StateTrackingStatus.CONFIDENT:
            return int(selected_state), result

        if status in (StateTrackingStatus.NO_ELIGIBLE_ROOT,
                      StateTrackingStatus.INVALID_RESPONSE):
            raise self._tracking_error(result)

        policy = self.state_tracker.failure_policy
        if (policy == 'strict' or
                (status is StateTrackingStatus.LOW_OVERLAP and
                 getattr(self, '_tracking_recovery_active', False))):
            raise self._tracking_error(result)

        diagnostics = result.to_dict()
        diagnostics.pop('status', None)
        diagnostics['provisional'] = True
        diagnostics['reference_update_recommended'] = False
        diagnostics['policy_action'] = policy

        if policy == 'best_effort':
            if status is StateTrackingStatus.AMBIGUOUS:
                raise self._tracking_error(result)
            diagnostics['warnings'] = list(diagnostics['warnings']) + [
                'BEST_EFFORT: using the unique argmax provisionally; the last '
                'valid reference is preserved'
            ]
            provisional = StateTrackingResult(status, diagnostics)
            self._print_provisional_warning(provisional)
            return int(selected_state), provisional

        # Adiabatic fallback: energy-order the spin-compatible roots for SF, or
        # all eligible response roots for ordinary TDA/TDDFT. No rotated
        # diabatic state is constructed; a real Serenity raw-root gradient is
        # always requested.
        candidates = [
            row for row in diagnostics['candidate_table']
            if row['eligible'] and row.get('spin_compatible') is not False and
            row['excitation_energy_ev'] is not None
        ]
        if not candidates:
            raise self._tracking_error(result)
        adiabatic = min(
            candidates, key=lambda row: row['excitation_energy_ev'])
        selected_state = int(adiabatic['raw_root'])
        diagnostics.update({
            'new_state': selected_state,
            'selected_by_overlap': int(result['new_state']),
            'reference_update_recommended': True,
        })
        spinflip = bool(getattr(getattr(self, 'rsp_driver', None),
                               'spinflip', False))
        manifold = 'inferred target-spin ' if spinflip else ''
        diagnostics['warnings'] = list(diagnostics['warnings']) + [
            'ADIABATIC: following the lowest-energy ' + manifold +
            f'raw root {selected_state}'
        ]
        provisional = StateTrackingResult(status, diagnostics)
        self._print_provisional_warning(provisional)
        return selected_state, provisional

    def _tracking_error(self, result):
        message = (
            'Serenity state tracking failed with status '
            f'{result.status.value}: reference root {result.get("old_state")}, '
            f'candidate root {result.get("new_state")}, normalized overlap '
            f'{result.get("max_overlap")}, second/best ratio '
            f'{result.get("overlap_ratio")}.')
        self.tracking_info = result.to_dict()
        self._print_tracking_diagnostics(self.tracking_info)
        return StateTrackingError(message, result)

    def _print_provisional_warning(self, result):
        self.ostream.print_info(
            'WARNING: Serenity state assignment is provisional '
            f'({result.status.value}, policy='
            f'{self.state_tracker.failure_policy}).')
        self.ostream.flush()

    def _configure_excited_gradient_task(self, grad_task, lr_restart=False):
        """
        Configures an excited-state gradient task for ``state_deriv_index``.

        :param lr_restart:
            Serenity restart flag of the task's LRSCF step; only True when the
            stored LR solution was obtained at this geometry and SCF solution.
        """

        # Print level NORMAL keeps Serenity's convergence and restart
        # messages in the captured output (see parse_serenity_lr_output).
        if hasattr(grad_task, 'generalSettings'):
            grad_task.generalSettings.printLevel = (
                spy.GLOBAL_PRINT_LEVELS.NORMAL)

        grad_task.settings.gradType = 'analytical'
        grad_task.settings.excMethod = self.exc_method
        grad_task.settings.excGradList = [int(self.state_deriv_index)]

        nstates_req = max(int(self.state_deriv_index),
                          int(self.rsp_driver.nstates))

        grad_task.settings.lrscfSettings.method = self.exc_method
        grad_task.settings.lrscfSettings.nEigen = int(nstates_req)
        grad_task.settings.lrscfSettings.restart = bool(lr_restart)
        if self.rsp_driver.conv_thresh is not None:
            grad_task.settings.lrscfSettings.conv = float(
                self.rsp_driver.conv_thresh)
            

        if self.rsp_driver.max_cycles is not None:
            grad_task.settings.lrscfSettings.maxCycles = int(
                self.rsp_driver.max_cycles)

        if self.rsp_driver.max_subspace_dimension is not None:
            grad_task.settings.lrscfSettings.maxSubspaceDimension = int(
                self.rsp_driver.max_subspace_dimension)

        if self.rsp_driver.densfit_j is not None:
            grad_task.settings.lrscfSettings.densFitJ = self.rsp_driver.densfit_j

        if self.rsp_driver.grid_accuracy is not None:
            grad_task.settings.lrscfSettings.grid.accuracy = int(
                self.rsp_driver.grid_accuracy)

        if self.rsp_driver.small_grid_accuracy is not None:
            grad_task.settings.lrscfSettings.grid.smallGridAccuracy = int(
                self.rsp_driver.small_grid_accuracy)
        if self.rsp_driver.spinflip:
            grad_task.settings.excMethod = 'sftda'
            if hasattr(grad_task.settings.lrscfSettings, 'scfstab'):
                grad_task.settings.lrscfSettings.scfstab = 'spinflip'
            elif hasattr(grad_task.settings, 'scfstab'):
                grad_task.settings.scfstab = 'spinflip'

    def _extract_excitation_energy(self, rsp_results):
        eig = rsp_results['eigenvalues']

        state = int(self.state_deriv_index)
        nst = len(eig)

        errmsg = 'SerenityExcitedStateGradientDriver: requested state index '
        errmsg += f'{state} but only {nst} state(s) are available.'
        assert_msg_critical(state <= nst, errmsg)

        return float(eig[state - 1])
    
    def set_state_tracker(self, tracker):
        assert_msg_critical(
            tracker is None or isinstance(tracker, TransitionDensityTracker),
            "SerenityExcitedStateGradientDriver: invalid transition-density tracker.")
        self.state_tracker = tracker
        self.tracking_info = None
        self.tracking_history = []

    def _print_tracking_diagnostics(self, info):
        """Prints one Serenity transition-density tracking decision."""

        spinflip = bool(getattr(getattr(self, 'rsp_driver', None),
                               'spinflip', False))
        title = ('Serenity Spin-Flip State Tracking'
                 if spinflip else
                 'Serenity TDA/TDDFT State Tracking')
        self.ostream.print_header(title)
        self.ostream.print_info(
            f"Status                     : {info.get('status', 'CONFIDENT')}")
        self.ostream.print_info(
            f"Previous raw root          : {info['old_state']}")
        self.ostream.print_info(
            f"Selected candidate root    : {info.get('new_state')}")
        self.ostream.print_info(
            'Hard-eligible roots         : ' +
            ' '.join(str(state)
                     for state in info.get('allowed_states', [])))
        if info.get('spin_diagnostic_roots') is not None:
            self.ostream.print_info(
                'Target-spin diagnostics : ' +
                ' '.join(str(state)
                         for state in info['spin_diagnostic_roots']))
        if info.get('s2_tolerance_roots') is not None:
            # The roots that are actually spin-compatible: target multiplicity
            # and within the <S^2> tolerance.  Roots listed above but missing
            # here are spin-incomplete SF roots.
            self.ostream.print_info(
                'Spin-compatible roots   : ' +
                ' '.join(str(state)
                         for state in info['s2_tolerance_roots']))
        if info.get('initialized', False):
            self.ostream.print_info(
                'Tracking reference initialized at the current geometry.')
        else:
            if info.get('max_overlap') is not None:
                self.ostream.print_info(
                    f"Target overlap             : "
                    f"{info['max_overlap']:.8f}")
            if info.get('second_overlap') is not None:
                self.ostream.print_info(
                    f"Second eligible overlap    : "
                    f"{info['second_overlap']:.8f}")
            if info.get('overlap_ratio') is not None:
                self.ostream.print_info(
                    f"Ambiguity ratio            : "
                    f"{info['overlap_ratio']:.8f}")
            self.ostream.print_info(
                f"Assignment confident       : "
                f"{info.get('assignment_confident', False)}")
            if info.get('global_state') is not None:
                self.ostream.print_info(
                    f"Hungarian target root      : {info['global_state']}")
        self.ostream.print_info(
            f"Raw root changed           : {info.get('swapped', False)}")
        if info.get('gradient_task_roots') is not None:
            self.ostream.print_info(
                'Gradient task root sequence: ' +
                ' '.join(str(state)
                         for state in info['gradient_task_roots']))

        if info.get('candidate_table'):
            has_spin = any(
                row.get('s2') is not None or
                row.get('inferred_multiplicity') is not None
                for row in info['candidate_table'])
            if has_spin:
                self.ostream.print_info(
                    'Candidates: raw  energy/eV residual conv       <S^2> mult '
                    'raw-overlap norm-overlap eligible exclusion')
            else:
                self.ostream.print_info(
                    'Candidates: raw  energy/eV residual conv '
                    'raw-overlap norm-overlap eligible exclusion')
            for row in info['candidate_table']:
                def value(key, fmt):
                    item = row.get(key)
                    return 'n/a' if item is None else format(item, fmt)

                prefix = (
                    f"  {row['raw_root']:3d} "
                    f"{value('excitation_energy_ev', '.8f'):>12s} "
                    f"{value('solver_residual', '.3e'):>8s} "
                    f"{str(row.get('solver_converged')):>5s} ")
                if has_spin:
                    prefix += (
                        f"{value('s2', '.6f'):>10s} "
                        f"{str(row.get('inferred_multiplicity')):>4s} ")
                self.ostream.print_info(
                    prefix +
                    f"{value('raw_overlap', '.8f'):>11s} "
                    f"{value('normalized_overlap', '.8f'):>12s} "
                    f"{str(row['eligible']):>8s} "
                    f"{row.get('exclusion_reason') or '-'}")
        self.ostream.print_blank()
        self.ostream.flush()

    def commit_tracking_reference(self):
        """Commits the staged reference after an accepted nuclear step."""

        if self.state_tracker is None:
            return False
        if self.rank == mpi_master():
            committed = self.state_tracker.commit()
            reference_state = int(self.state_tracker.reference_state)
            if committed and self.tracking_info is not None:
                self.tracking_info['reference_updated'] = True
                self.tracking_info['reference_staged'] = False
                if self.tracking_history:
                    self.tracking_history[-1]['reference_updated'] = True
                    self.tracking_history[-1]['reference_staged'] = False
        else:
            committed = None
            reference_state = None
        committed, reference_state = self.comm.bcast(
            (committed, reference_state), root=mpi_master())
        if committed:
            self.state_deriv_index = int(reference_state)
        return committed

    def rollback_tracking_reference(self):
        """Rolls back a rejected trial and its selected raw-root identity."""

        if self.state_tracker is None:
            return False
        if self.rank == mpi_master():
            rolled_back = self.state_tracker.rollback()
            reference_state = int(self.state_tracker.reference_state)
        else:
            rolled_back = None
            reference_state = None
        rolled_back, reference_state = self.comm.bcast(
            (rolled_back, reference_state), root=mpi_master())
        self.state_deriv_index = int(reference_state)
        self._tracking_applied_in_compute = False
        return rolled_back

    def begin_tracking_retry(self):
        """Restores the committed root and starts midpoint subdivision."""

        if self.state_tracker is None:
            return False
        if self.rank == mpi_master():
            self.state_tracker.begin_retry_chain()
            reference_state = int(self.state_tracker.reference_state)
        else:
            reference_state = None
        reference_state = self.comm.bcast(
            reference_state, root=mpi_master())
        self.state_deriv_index = int(reference_state)
        return True

    def promote_tracking_reference_for_retry(self):
        """Promotes a successful midpoint only within the retry transaction."""

        if self.state_tracker is None:
            return False
        if self.rank == mpi_master():
            promoted = self.state_tracker.promote_pending_for_retry()
            reference_state = int(self.state_tracker.reference_state)
        else:
            promoted = None
            reference_state = None
        promoted, reference_state = self.comm.bcast(
            (promoted, reference_state), root=mpi_master())
        if promoted:
            self.state_deriv_index = int(reference_state)
        return promoted

    def track_state(self, molecule, recompute_on_switch=True,
                update_reference=True):
        if self.state_tracker is None:
            self.tracking_info = None
            return None

        # Analytical calculations now track as part of compute(), so energy,
        # gradient and root are selected atomically for both ordinary TDDFT and
        # spin-flip. OpenMM may still call track_state() afterwards; return the
        # completed decision instead of tracking twice.
        if self._tracking_applied_in_compute:
            if update_reference:
                self.commit_tracking_reference()
            return self.tracking_info

        assert_msg_critical(
            self.last_lrscf_controller is not None,
            "track_state: call compute() before track_state()."
        )
        old_state = int(self.state_deriv_index)
        result = self.state_tracker.track(
            self._last_tracking_system,
            self.last_lrscf_controller,
            self._last_tracking_mode,
            active_reference_state=old_state,
        )
        info = result.to_dict()
        new_state = int(info["new_state"])
        info["gradient_recomputed"] = False
        info["reference_updated"] = False
        if new_state != old_state:
            self.set_state_deriv_index(new_state)
            if recompute_on_switch:
                self.gradient = self._compute_analytical_master(molecule)
                info["gradient_recomputed"] = True
        if update_reference:
            self.state_tracker.accept_reference(
                self._last_tracking_system,
                self.last_lrscf_controller,
                self._last_tracking_mode,
                int(self.state_deriv_index),
            )
            info["reference_updated"] = True
        self.tracking_info = info

        return self.tracking_info
