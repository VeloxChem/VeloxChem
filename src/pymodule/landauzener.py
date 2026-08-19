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

"""
Pure numerical core of the adiabatic Landau-Zener surface-hopping model.

This module deliberately contains **no** VeloxChem or OpenMM imports so that
the mathematics can be unit tested in isolation with plain NumPy arrays.

All quantities are in atomic units:

    * energies and energy gaps      : Hartree
    * time and electronic timestep  : atomic units of time
    * gap curvature                 : Hartree / (atomic time unit)**2

With hbar = 1 the adiabatic Landau-Zener hopping probability evaluated at a
local minimum of the raw adiabatic energy gap Z_ab(t) = |E_b - E_a| reads

.. math::

    P_{a \\rightarrow b} =
        \\exp\\left[ -\\frac{\\pi}{2}
        \\sqrt{ \\frac{Z_c^{3}}{\\ddot{Z}_c} } \\right]

where :math:`Z_c` is the interpolated minimum gap and :math:`\\ddot{Z}_c` its
second time derivative at the minimum.
"""

from dataclasses import dataclass, asdict
from collections import deque
import numpy as np

# Default ceiling for the Landau-Zener exponent.  exp(-700) underflows to
# exactly 0.0 in IEEE double precision, so clipping here loses no information
# while protecting np.exp from spurious overflow warnings.
DEFAULT_EXPONENT_CLIP = 700.0


@dataclass
class GapFitResult:
    """
    Result of a local fit of a raw adiabatic energy gap around a candidate
    minimum.

    :param minimum_gap:
        Interpolated minimum gap Z_c, in Hartree.
    :param curvature:
        Second time derivative of the gap at the minimum, in Hartree per
        atomic time unit squared.
    :param event_time_offset:
        Fitted event time tau_c relative to the central frame, in atomic
        time units.
    :param fit_quality:
        Quality of the local fit in [0, 1].  1.0 means the quadratic model
        reproduces the sampled gap values exactly.
    :param n_points:
        Number of gap samples used by the fit.
    """

    minimum_gap: float
    curvature: float
    event_time_offset: float
    fit_quality: float
    n_points: int


@dataclass
class LandauZenerEvent:
    """
    A validated Landau-Zener hopping candidate.

    All energies are in Hartree and all times in atomic time units.
    """

    source_state: int
    target_state: int
    central_step: int
    event_time: float
    minimum_gap: float
    gap_curvature: float
    exponent: float
    probability: float
    fit_quality: float
    tracking_quality: float
    event_time_offset: float = 0.0

    def as_dict(self):
        """
        :return:
            A plain dictionary representation for machine-readable logging.
        """

        return asdict(self)


def landau_zener_probability(minimum_gap,
                             gap_curvature,
                             exponent_clip=DEFAULT_EXPONENT_CLIP):
    """
    Evaluates the adiabatic Landau-Zener hopping probability.

    .. math::

        P = \\exp\\left[-\\frac{\\pi}{2}
            \\sqrt{\\frac{Z_c^3}{\\ddot Z_c}}\\right]

    :param minimum_gap:
        The interpolated minimum gap Z_c, in Hartree.  Must be positive.
    :param gap_curvature:
        The gap curvature at the minimum, in Hartree per atomic time unit
        squared.  Must be positive.
    :param exponent_clip:
        Upper bound applied to the exponent before exponentiation.

    :return:
        The tuple (exponent, probability).  Returns (inf, 0.0) for invalid
        input rather than raising, so that callers can reject the event.
    """

    gap = float(minimum_gap)
    curv = float(gap_curvature)

    if not (np.isfinite(gap) and np.isfinite(curv)):
        return float('inf'), 0.0
    if gap <= 0.0 or curv <= 0.0:
        return float('inf'), 0.0

    exponent = 0.5 * np.pi * np.sqrt(gap**3 / curv)

    if not np.isfinite(exponent):
        return float('inf'), 0.0

    probability = float(np.exp(-min(exponent, float(exponent_clip))))

    return float(exponent), probability


class GapMinimumFitter:
    """
    Base class for local fits of a raw adiabatic energy gap around a minimum.

    Subclasses implement :func:`fit` and declare how many uniformly spaced
    samples they consume via :attr:`n_points`.  The dynamics controller only
    ever talks to this interface, so a five-point or Savitzky-Golay fitter can
    replace the three-point fitter without touching the controller.
    """

    n_points = 3

    def fit(self, gaps, timestep):
        """
        Fits a local quadratic through uniformly spaced gap samples.

        :param gaps:
            Sequence of :attr:`n_points` gap values, in Hartree.
        :param timestep:
            The uniform electronic timestep h, in atomic time units.

        :return:
            A :class:`GapFitResult`, or None when the samples do not describe
            an acceptable minimum.
        """

        raise NotImplementedError


class ThreePointGapFitter(GapMinimumFitter):
    """
    Three-point quadratic fit of a gap minimum.

    Given Z_{n-1}, Z_n, Z_{n+1} at uniform spacing h the exact quadratic
    through the three points yields

    .. math::

        \\ddot Z_c &= \\frac{Z_{n-1} - 2 Z_n + Z_{n+1}}{h^2} \\\\
        \\tau_c &= \\frac{h (Z_{n-1} - Z_{n+1})}
                        {2 (Z_{n-1} - 2 Z_n + Z_{n+1})} \\\\
        Z_c &= Z_n - \\frac{(Z_{n+1} - Z_{n-1})^2}
                           {8 (Z_{n-1} - 2 Z_n + Z_{n+1})}

    The quadratic is exact through three points, so :attr:`fit_quality` is
    always 1.0 here; use :class:`FivePointGapFitter` when a genuine residual
    based quality measure is wanted.
    """

    n_points = 3

    def fit(self, gaps, timestep):
        """
        See :func:`GapMinimumFitter.fit`.
        """

        gaps = np.asarray(gaps, dtype=float)

        if gaps.size != 3:
            return None
        if not np.all(np.isfinite(gaps)):
            return None

        h = float(timestep)
        if not np.isfinite(h) or h <= 0.0:
            return None

        z_prev, z_cent, z_next = (float(gaps[0]), float(gaps[1]),
                                  float(gaps[2]))

        denominator = z_prev - 2.0 * z_cent + z_next

        # A nonpositive denominator means the samples do not bracket a
        # convex minimum; reject rather than extrapolate.
        if not np.isfinite(denominator) or denominator <= 0.0:
            return None

        curvature = denominator / (h * h)
        tau_c = h * (z_prev - z_next) / (2.0 * denominator)
        z_c = z_cent - (z_next - z_prev)**2 / (8.0 * denominator)

        if not (np.isfinite(curvature) and np.isfinite(tau_c)
                and np.isfinite(z_c)):
            return None
        if curvature <= 0.0 or z_c <= 0.0:
            return None
        # The fitted minimum must lie inside the sampled interval.
        if abs(tau_c) > h:
            return None

        return GapFitResult(minimum_gap=z_c,
                            curvature=curvature,
                            event_time_offset=tau_c,
                            fit_quality=1.0,
                            n_points=3)


class FivePointGapFitter(GapMinimumFitter):
    """
    Five-point least-squares quadratic fit of a gap minimum.

    The gap is modelled as Z(t) = c0 + c1 t + c2 t^2 with t measured from the
    central frame in units of the timestep.  The fit is overdetermined, so the
    normalized residual provides a genuine fit-quality measure that can reject
    minima whose local shape is not quadratic.
    """

    n_points = 5

    def fit(self, gaps, timestep):
        """
        See :func:`GapMinimumFitter.fit`.
        """

        gaps = np.asarray(gaps, dtype=float)

        if gaps.size != 5:
            return None
        if not np.all(np.isfinite(gaps)):
            return None

        h = float(timestep)
        if not np.isfinite(h) or h <= 0.0:
            return None

        # Offsets in units of h, centred on the middle frame.
        s = np.array([-2.0, -1.0, 0.0, 1.0, 2.0])
        design = np.vstack([np.ones_like(s), s, s * s]).T

        coeffs, *_ = np.linalg.lstsq(design, gaps, rcond=None)
        c0, c1, c2 = (float(coeffs[0]), float(coeffs[1]), float(coeffs[2]))

        if not np.all(np.isfinite(coeffs)):
            return None
        if c2 <= 0.0:
            return None

        # d^2 Z / dt^2 with t in atomic time units.
        curvature = 2.0 * c2 / (h * h)
        # Vertex of the parabola, in units of h, then converted to au.
        tau_c = -0.5 * c1 / c2 * h
        z_c = c0 - 0.25 * c1 * c1 / c2

        if not (np.isfinite(curvature) and np.isfinite(tau_c)
                and np.isfinite(z_c)):
            return None
        if curvature <= 0.0 or z_c <= 0.0:
            return None
        if abs(tau_c) > h:
            return None

        residual = gaps - design @ coeffs
        scale = float(np.max(np.abs(gaps)))
        if scale <= 0.0:
            return None
        rms = float(np.sqrt(np.mean(residual**2))) / scale
        fit_quality = float(np.clip(1.0 - rms, 0.0, 1.0))

        return GapFitResult(minimum_gap=z_c,
                            curvature=curvature,
                            event_time_offset=tau_c,
                            fit_quality=fit_quality,
                            n_points=5)


GAP_FITTERS = {
    'three_point': ThreePointGapFitter,
    'five_point': FivePointGapFitter,
}


def create_gap_fitter(name):
    """
    Instantiates a gap fitter by name.

    :param name:
        One of the keys of :data:`GAP_FITTERS`.

    :return:
        A :class:`GapMinimumFitter` instance.
    """

    key = str(name).lower()

    if key not in GAP_FITTERS:
        valid = ', '.join(sorted(GAP_FITTERS))
        raise ValueError(
            f'create_gap_fitter: unknown gap fitter {name}; valid: {valid}')

    return GAP_FITTERS[key]()


class LandauZenerEventDetector:
    """
    Detects local minima of raw adiabatic energy gaps and converts them
    into validated :class:`LandauZenerEvent` objects.

    The detector is purely numerical: it consumes state energies in raw
    adiabatic-root order, keyed by raw source/target indices, and knows nothing
    about character labels, wavefunctions, OpenMM or VeloxChem drivers.  State
    tracking may gate an event through ``tracking_quality`` but does not
    permute these energies or indices.

    :param timestep:
        The uniform electronic timestep, in atomic time units.
    :param fitter:
        A :class:`GapMinimumFitter` instance, or a name accepted by
        :func:`create_gap_fitter`.
    :param gap_screening_threshold:
        Gaps above this value (Hartree) never produce events.
    :param minimum_curvature:
        Curvatures at or below this value are rejected.
    :param exponent_clip:
        Ceiling applied to the Landau-Zener exponent.
    :param minimum_fit_quality:
        Fits below this quality are rejected.
    :param multistate_warning_gap:
        When more than one candidate gap is below this value at the same
        central frame the event is flagged as non-isolated.
    """

    def __init__(self,
                 timestep,
                 fitter='three_point',
                 gap_screening_threshold=None,
                 minimum_curvature=0.0,
                 exponent_clip=DEFAULT_EXPONENT_CLIP,
                 minimum_fit_quality=0.0,
                 multistate_warning_gap=None):

        self.timestep = float(timestep)
        if not np.isfinite(self.timestep) or self.timestep <= 0.0:
            raise ValueError(
                'LandauZenerEventDetector: timestep must be positive.')

        if isinstance(fitter, GapMinimumFitter):
            self.fitter = fitter
        else:
            self.fitter = create_gap_fitter(fitter)

        self.gap_screening_threshold = gap_screening_threshold
        self.minimum_curvature = float(minimum_curvature)
        self.exponent_clip = float(exponent_clip)
        self.minimum_fit_quality = float(minimum_fit_quality)
        self.multistate_warning_gap = multistate_warning_gap

        n_hist = self.fitter.n_points

        # (source, target) -> deque of gap values
        self._gap_history = {}
        # (source, target) -> deque of step indices
        self._step_history = {}
        # (source, target) -> last consumed central step
        self._last_event_step = {}

        self._history_length = n_hist
        self.last_multistate_warning = False

    @property
    def history_length(self):
        """
        :return:
            Number of gap samples required before an event can be detected.
        """

        return self._history_length

    def reset(self):
        """
        Clears all stored gap history and event bookkeeping.
        """

        self._gap_history.clear()
        self._step_history.clear()
        self._last_event_step.clear()
        self.last_multistate_warning = False

    def reset_pairs_involving(self, state):
        """
        Clears gap history for every pair that involves ``state``.

        Called after an accepted hop, where gaps referenced to the old active
        state are no longer meaningful.

        :param state:
            Raw adiabatic state index.
        """

        state = int(state)
        stale = [key for key in self._gap_history if state in key]

        for key in stale:
            self._gap_history.pop(key, None)
            self._step_history.pop(key, None)
            self._last_event_step.pop(key, None)

    def push(self, step, active_state, state_energies, candidate_states):
        """
        Appends one frame of raw adiabatic energies to the gap history.

        :param step:
            Integer index of the electronic frame.
        :param active_state:
            Raw adiabatic index of the currently active state.
        :param state_energies:
            Mapping or sequence of raw adiabatic state energies, in Hartree.
        :param candidate_states:
            Iterable of raw adiabatic indices considered as hopping targets.
        """

        active_state = int(active_state)
        e_active = float(_lookup_energy(state_energies, active_state))

        for target in candidate_states:
            target = int(target)
            if target == active_state:
                continue

            e_target = float(_lookup_energy(state_energies, target))
            gap = abs(e_target - e_active)

            key = (active_state, target)

            if key not in self._gap_history:
                self._gap_history[key] = deque(maxlen=self._history_length)
                self._step_history[key] = deque(maxlen=self._history_length)

            self._gap_history[key].append(gap)
            self._step_history[key].append(int(step))

    def detect(self, tracking_quality=1.0, minimum_tracking_quality=0.0):
        """
        Scans the stored histories for validated Landau-Zener events.

        A candidate exists when the central sample of a full history window is
        strictly smaller than both of its neighbours.  Candidates are rejected
        when the local fit fails, the curvature or gap is nonpositive, the
        fitted event time leaves the sampled interval, the fit quality is too
        low, or the state tracking confidence is too low.

        :param tracking_quality:
            Confidence of the state tracking for the current window, in
            [0, 1].
        :param minimum_tracking_quality:
            Events below this tracking confidence are rejected.

        :return:
            A list of :class:`LandauZenerEvent`, ordered by the deterministic
            selection policy: earliest fitted event time first, ties broken by
            smallest minimum gap.
        """

        events = []
        small_gap_count = 0

        for key, history in self._gap_history.items():
            if len(history) < self._history_length:
                continue

            source, target = key
            gaps = list(history)
            steps = list(self._step_history[key])

            centre = self._history_length // 2
            z_centre = gaps[centre]
            central_step = steps[centre]

            if (self.multistate_warning_gap is not None
                    and z_centre < float(self.multistate_warning_gap)):
                small_gap_count += 1

            # Do not re-process a minimum that has already been consumed.
            if self._last_event_step.get(key) == central_step:
                continue

            if not np.isfinite(z_centre):
                continue

            if (self.gap_screening_threshold is not None
                    and z_centre > float(self.gap_screening_threshold)):
                continue

            # Strict local minimum with respect to the immediate neighbours.
            if not (z_centre < gaps[centre - 1] and z_centre < gaps[centre + 1]):
                continue

            fit = self.fitter.fit(gaps, self.timestep)

            if fit is None:
                continue
            if fit.curvature <= self.minimum_curvature:
                continue
            if fit.fit_quality < self.minimum_fit_quality:
                continue
            if float(tracking_quality) < float(minimum_tracking_quality):
                continue

            exponent, probability = landau_zener_probability(
                fit.minimum_gap, fit.curvature, self.exponent_clip)

            if not np.isfinite(exponent):
                continue

            events.append(
                LandauZenerEvent(
                    source_state=int(source),
                    target_state=int(target),
                    central_step=int(central_step),
                    event_time=float(central_step) * self.timestep +
                    fit.event_time_offset,
                    minimum_gap=float(fit.minimum_gap),
                    gap_curvature=float(fit.curvature),
                    exponent=float(exponent),
                    probability=float(probability),
                    fit_quality=float(fit.fit_quality),
                    tracking_quality=float(tracking_quality),
                    event_time_offset=float(fit.event_time_offset)))

        self.last_multistate_warning = (small_gap_count > 1)

        # Deterministic policy: earliest fitted event first, then tightest gap.
        events.sort(key=lambda ev: (round(ev.event_time, 12), ev.minimum_gap))

        return events

    def mark_processed(self, event):
        """
        Records that an event has been consumed so that the same minimum is
        not detected again on the following frame.

        :param event:
            The :class:`LandauZenerEvent` that has been processed.
        """

        key = (int(event.source_state), int(event.target_state))
        self._last_event_step[key] = int(event.central_step)

    def get_history_snapshot(self):
        """
        :return:
            A picklable snapshot of the detector history for checkpointing.
        """

        return {
            'gap_history': {
                key: list(val) for key, val in self._gap_history.items()
            },
            'step_history': {
                key: list(val) for key, val in self._step_history.items()
            },
            'last_event_step': dict(self._last_event_step),
        }

    def set_history_snapshot(self, snapshot):
        """
        Restores detector history from a checkpoint snapshot.

        :param snapshot:
            A dictionary produced by :func:`get_history_snapshot`.
        """

        self.reset()

        for key, val in snapshot.get('gap_history', {}).items():
            self._gap_history[tuple(key)] = deque(
                val, maxlen=self._history_length)

        for key, val in snapshot.get('step_history', {}).items():
            self._step_history[tuple(key)] = deque(
                val, maxlen=self._history_length)

        self._last_event_step.update(
            {tuple(k): v for k, v in snapshot.get('last_event_step', {}).items()})


def _lookup_energy(state_energies, index):
    """
    Reads one state energy from either a mapping or a sequence.

    :param state_energies:
        Mapping from tracked index to energy, or an indexable sequence.
    :param index:
        Tracked state index.

    :return:
        The energy in Hartree.
    """

    if isinstance(state_energies, dict):
        if index not in state_energies:
            raise KeyError(
                f'Landau-Zener detector: missing energy for state {index}.')
        return state_energies[index]

    energies = np.asarray(state_energies, dtype=float)

    if index < 0 or index >= energies.size:
        raise IndexError(
            f'Landau-Zener detector: state index {index} out of range.')

    return energies[index]
