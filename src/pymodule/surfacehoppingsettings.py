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
Configuration and unit handling for Landau-Zener surface-hopping dynamics.

All unit conversions used by the surface-hopping layer are defined here once
and imported elsewhere, so that no approximate physical constant is duplicated
across modules.
"""

from dataclasses import dataclass, asdict
import hashlib
import json
import math

from .veloxchemlib import (bohr_in_angstrom, amu_in_kg, amu_in_electron_mass,
                           hartree_in_kjpermol)
from .errorhandler import assert_msg_critical


def _atomic_time_unit_in_seconds():
    """
    Derives the atomic unit of time from the CODATA constants already exposed
    by VeloxChem, so that no new numerical constant is introduced.

    Using ``t_au = a_0 sqrt(m_e / E_h)`` with

        * ``a_0``  from :func:`bohr_in_angstrom`
        * ``m_e``  from :func:`amu_in_kg` and :func:`amu_in_electron_mass`
        * ``E_h``  from :func:`hartree_in_kjpermol` together with the Avogadro
          constant recovered as ``N_A = 1e-3 / amu_in_kg()``

    :return:
        The atomic unit of time in seconds.
    """

    # E_h in joule: kJ/mol -> J/mol -> J, with N_A = 1e-3 kg/mol / amu_in_kg.
    hartree_in_joule = hartree_in_kjpermol() * 1.0e6 * amu_in_kg()
    bohr_in_meter = bohr_in_angstrom() * 1.0e-10
    electron_mass_in_kg = amu_in_kg() / amu_in_electron_mass()

    return bohr_in_meter * math.sqrt(electron_mass_in_kg / hartree_in_joule)


#: Conversion factor from femtoseconds to atomic time units.
FS_IN_AU = 1.0e-15 / _atomic_time_unit_in_seconds()

#: Conversion factor from atomic time units to femtoseconds.
AU_IN_FS = 1.0 / FS_IN_AU

#: Conversion factor from nanometre (OpenMM length unit) to bohr.
NM_IN_BOHR = 10.0 / bohr_in_angstrom()

#: Conversion factor from bohr to nanometre.
BOHR_IN_NM = 1.0 / NM_IN_BOHR

#: Conversion factor from nm/ps (OpenMM velocity unit) to bohr per atomic
#: time unit.
NM_PER_PS_IN_AU_VELOCITY = NM_IN_BOHR / (1000.0 * FS_IN_AU)

#: Conversion factor from bohr per atomic time unit to nm/ps.
AU_VELOCITY_IN_NM_PER_PS = 1.0 / NM_PER_PS_IN_AU_VELOCITY

#: Conversion factor from unified atomic mass units to electron masses.
AMU_IN_ELECTRON_MASS = amu_in_electron_mass()


@dataclass
class SurfaceHoppingSettings:
    """
    Validated configuration of a Landau-Zener surface-hopping trajectory.

    Attribute names follow the VeloxChem input-keyword convention (lower case
    with underscores) and map one-to-one onto user input keywords.

    :param surface_hopping_method:
        Hopping model. Only ``'landau_zener'`` is currently implemented.
    :param number_of_states:
        Total number of native adiabatic states, including the ground state
        for conventional ground-plus-excitation providers.
    :param initial_active_state:
        Native adiabatic-root index of the initially active state. 0 is the
        ground state for conventional providers or the lowest target for a
        spin-flip provider.
    :param allowed_transitions:
        ``'adjacent'`` restricts hopping to neighbouring adiabatic roots;
        ``'all'`` permits any pair involving the active state; a list of
        ``(a, b)`` pairs defines an explicit set.
    :param electronic_timestep:
        Electronic/nuclear timestep in femtoseconds.
    :param gap_fit_points:
        Number of samples in the local gap fit, 3 or 5.
    :param gap_screening_threshold:
        Gaps above this value, in Hartree, never generate events.
    :param minimum_curvature:
        Minimum accepted gap curvature, in Hartree per atomic time unit
        squared.
    :param maximum_lz_exponent:
        Ceiling applied to the Landau-Zener exponent.
    :param minimum_fit_quality:
        Minimum accepted local-fit quality, in [0, 1].
    :param electronic_backend:
        Explicit electronic backend: ``'openqp_mrsf'``, ``'serenity_sf'`` or
        ``'veloxchem_tda'``.  Production runs must name a spin-flip backend;
        it is never inferred and never silently replaced.
    :param state_tracking_method:
        ``'backend_overlap'``, ``'transition_density'``, ``'descriptor'`` or
        ``'none'``.  The last two cannot see a change of electronic
        character and are refused in production; see
        ``allow_unsafe_state_tracking``.
    :param minimum_tracking_overlap:
        Minimum accepted state-tracking confidence, in [0, 1].
    :param maximum_tracking_ambiguity:
        Maximum accepted ratio of the best rejected similarity to the
        assigned one, per tracked state, in [0, 1].
    :param allow_unsafe_state_tracking:
        Opt-in switch that permits ``'descriptor'`` or ``'none'`` tracking
        with a production backend.  Such a run establishes state identity
        from excitation-energy order alone and is diagnostic only; it must
        never be used for production data.
    :param momentum_rescaling:
        ``'global_qm'`` or ``'gap_gradient'``.
    :param rescaling_atoms:
        ``'qm'`` or an explicit list of atom indices.
    :param remove_translation:
        Exclude and restore center-of-mass translation when rescaling.
    :param remove_rotation:
        Exclude and restore rigid-body rotation when rescaling.
    :param frustrated_hop_policy:
        Only ``'reject'`` is implemented: the active state is unchanged and
        velocities are untouched.
    :param random_seed:
        Seed of the surface-hopping random number generator.
    :param multistate_warning_gap:
        Diagnostic threshold, in Hartree, for non-isolated multistate events.
    :param hop_energy_tolerance:
        Maximum accepted energy-conservation residual on a hop, in Hartree.
    :param checkpoint_interval:
        Number of steps between checkpoint writes; 0 disables checkpointing.
    :param diagnostics_file:
        Path of the machine-readable per-step diagnostics file.
    :param checkpoint_file:
        Path of the restart checkpoint file.
    :param trajectory_id:
        Identifier written into every diagnostics record.
    """

    surface_hopping_method: str = 'landau_zener'
    number_of_states: int = 4
    initial_active_state: int = 1
    allowed_transitions: object = 'adjacent'
    electronic_timestep: float = 0.5
    gap_fit_points: int = 3
    gap_screening_threshold: float = 0.15
    minimum_curvature: float = 0.0
    maximum_lz_exponent: float = 700.0
    minimum_fit_quality: float = 0.0
    electronic_backend: str = 'veloxchem_tda'
    state_tracking_method: str = 'descriptor'
    minimum_tracking_overlap: float = 0.5
    maximum_tracking_ambiguity: float = 0.8
    allow_unsafe_state_tracking: bool = False
    momentum_rescaling: str = 'global_qm'
    rescaling_atoms: object = 'qm'
    remove_translation: bool = True
    remove_rotation: bool = False
    frustrated_hop_policy: str = 'reject'
    random_seed: int = 42
    multistate_warning_gap: float = 0.005
    hop_energy_tolerance: float = 1.0e-6
    checkpoint_interval: int = 10
    diagnostics_file: str = 'surface_hopping.jsonl'
    checkpoint_file: str = 'surface_hopping_checkpoint.npz'
    trajectory_id: str = 'traj0'

    #: Names of settings whose change invalidates a checkpoint.
    FINGERPRINT_KEYS = ('surface_hopping_method', 'number_of_states',
                        'allowed_transitions', 'electronic_timestep',
                        'gap_fit_points', 'gap_screening_threshold',
                        'minimum_curvature', 'maximum_lz_exponent',
                        'minimum_fit_quality', 'electronic_backend',
                        'state_tracking_method',
                        'minimum_tracking_overlap',
                        'maximum_tracking_ambiguity',
                        'allow_unsafe_state_tracking', 'momentum_rescaling',
                        'rescaling_atoms', 'remove_translation',
                        'remove_rotation', 'frustrated_hop_policy',
                        'hop_energy_tolerance')

    #: Electronic backends the production driver may select.
    PRODUCTION_BACKENDS = ('openqp_mrsf', 'serenity_sf')

    #: Every accepted backend identifier.
    VALID_BACKENDS = PRODUCTION_BACKENDS + ('veloxchem_tda',)

    @property
    def timestep_in_au(self):
        """
        :return:
            The electronic timestep in atomic time units.
        """

        return float(self.electronic_timestep) * FS_IN_AU

    @classmethod
    def from_dict(cls, sh_dict):
        """
        Builds settings from a user input dictionary, ignoring unknown keys
        only after reporting them.

        :param sh_dict:
            Dictionary of surface-hopping keywords.

        :return:
            A validated :class:`SurfaceHoppingSettings`.
        """

        sh_dict = dict(sh_dict or {})
        valid_keys = set(cls.__dataclass_fields__)
        unknown = sorted(set(sh_dict) - valid_keys)

        assert_msg_critical(
            not unknown,
            'SurfaceHoppingSettings: unknown keyword(s): ' + ', '.join(unknown))

        settings = cls(**sh_dict)
        settings.validate()

        return settings

    def validate(self):
        """
        Validates every setting before dynamics starts.

        Raises a critical error through :func:`assert_msg_critical` on the
        first invalid setting encountered.
        """

        assert_msg_critical(
            str(self.surface_hopping_method).lower() == 'landau_zener',
            'SurfaceHoppingSettings: only the landau_zener method is '
            'implemented.')

        assert_msg_critical(
            int(self.number_of_states) >= 2,
            'SurfaceHoppingSettings: number_of_states must be >= 2.')

        assert_msg_critical(
            0 <= int(self.initial_active_state) < int(self.number_of_states),
            'SurfaceHoppingSettings: initial_active_state must lie in '
            f'[0, {int(self.number_of_states) - 1}].')

        assert_msg_critical(
            float(self.electronic_timestep) > 0.0,
            'SurfaceHoppingSettings: electronic_timestep must be positive.')

        assert_msg_critical(
            int(self.gap_fit_points) in (3, 5),
            'SurfaceHoppingSettings: gap_fit_points must be 3 or 5.')

        assert_msg_critical(
            float(self.gap_screening_threshold) > 0.0,
            'SurfaceHoppingSettings: gap_screening_threshold must be '
            'positive.')

        assert_msg_critical(
            float(self.minimum_curvature) >= 0.0,
            'SurfaceHoppingSettings: minimum_curvature must be '
            'non-negative.')

        assert_msg_critical(
            float(self.maximum_lz_exponent) > 0.0,
            'SurfaceHoppingSettings: maximum_lz_exponent must be positive.')

        assert_msg_critical(
            0.0 <= float(self.minimum_fit_quality) <= 1.0,
            'SurfaceHoppingSettings: minimum_fit_quality must lie in [0, 1].')

        assert_msg_critical(
            0.0 <= float(self.minimum_tracking_overlap) <= 1.0,
            'SurfaceHoppingSettings: minimum_tracking_overlap must lie in '
            '[0, 1].')

        assert_msg_critical(
            0.0 <= float(self.maximum_tracking_ambiguity) <= 1.0,
            'SurfaceHoppingSettings: maximum_tracking_ambiguity must lie in '
            '[0, 1].')

        assert_msg_critical(
            str(self.state_tracking_method).lower()
            in ('backend_overlap', 'transition_density', 'descriptor', 'none'),
            'SurfaceHoppingSettings: state_tracking_method must be one of '
            'backend_overlap, transition_density, descriptor, none.')

        assert_msg_critical(
            str(self.electronic_backend).lower() in self.VALID_BACKENDS,
            'SurfaceHoppingSettings: electronic_backend must be one of ' +
            ', '.join(self.VALID_BACKENDS) + '.')

        if str(self.electronic_backend).lower() in self.PRODUCTION_BACKENDS:
            assert_msg_critical(
                str(self.state_tracking_method).lower() == 'backend_overlap'
                or bool(self.allow_unsafe_state_tracking),
                'SurfaceHoppingSettings: electronic_backend '
                f'{self.electronic_backend} requires '
                "state_tracking_method='backend_overlap'. The descriptor and "
                'identity trackers establish state identity from excitation-'
                'energy order alone, because neither spin-flip backend '
                'exposes an oscillator strength or a transition dipole; that '
                'is the silent energy-order fallback production surface '
                'hopping must not perform. Set '
                'allow_unsafe_state_tracking=True only for diagnostics.')

        assert_msg_critical(
            str(self.momentum_rescaling).lower()
            in ('global_qm', 'gap_gradient'),
            'SurfaceHoppingSettings: momentum_rescaling must be global_qm or '
            'gap_gradient.')

        assert_msg_critical(
            str(self.frustrated_hop_policy).lower() == 'reject',
            'SurfaceHoppingSettings: only the reject frustrated_hop_policy '
            'is implemented.')

        assert_msg_critical(
            float(self.hop_energy_tolerance) > 0.0,
            'SurfaceHoppingSettings: hop_energy_tolerance must be positive.')

        assert_msg_critical(
            int(self.checkpoint_interval) >= 0,
            'SurfaceHoppingSettings: checkpoint_interval must be '
            'non-negative.')

        self._validate_allowed_transitions()

        if self.gap_fit_points == 5 and self.minimum_fit_quality == 0.0:
            # Not an error, but the five-point fitter is only useful with an
            # active quality threshold.
            pass

    def _validate_allowed_transitions(self):
        """
        Validates the ``allowed_transitions`` specification.
        """

        spec = self.allowed_transitions

        if isinstance(spec, str):
            assert_msg_critical(
                spec.lower() in ('adjacent', 'all'),
                'SurfaceHoppingSettings: allowed_transitions must be '
                'adjacent, all, or a list of state pairs.')
            return

        try:
            pairs = [tuple(int(x) for x in pair) for pair in spec]
        except (TypeError, ValueError):
            assert_msg_critical(
                False, 'SurfaceHoppingSettings: allowed_transitions list must '
                'contain pairs of integers.')
            return

        n_states = int(self.number_of_states)

        for pair in pairs:
            assert_msg_critical(
                len(pair) == 2,
                'SurfaceHoppingSettings: allowed_transitions entries must be '
                'pairs.')
            assert_msg_critical(
                all(0 <= idx < n_states for idx in pair),
                'SurfaceHoppingSettings: allowed_transitions index out of '
                f'range [0, {n_states - 1}].')
            assert_msg_critical(
                pair[0] != pair[1],
                'SurfaceHoppingSettings: allowed_transitions pair must '
                'involve two different states.')

    def candidate_states(self, active_state):
        """
        Returns the adiabatic roots considered as hopping targets for a given
        active raw root.

        :param active_state:
            Native adiabatic-root index of the active state.

        :return:
            A sorted list of candidate target-state indices.
        """

        active_state = int(active_state)
        n_states = int(self.number_of_states)
        spec = self.allowed_transitions

        if isinstance(spec, str):
            if spec.lower() == 'all':
                candidates = [i for i in range(n_states) if i != active_state]
            else:
                candidates = [
                    i for i in (active_state - 1, active_state + 1)
                    if 0 <= i < n_states
                ]
        else:
            candidates = []
            for pair in spec:
                a, b = int(pair[0]), int(pair[1])
                if a == active_state:
                    candidates.append(b)
                elif b == active_state:
                    candidates.append(a)

        return sorted(set(candidates))

    def fingerprint(self):
        """
        Builds a hash of the settings that must not change across a restart.

        :return:
            A hexadecimal digest string.
        """

        payload = {}

        for key in self.FINGERPRINT_KEYS:
            value = getattr(self, key)
            if not isinstance(value, (str, int, float, bool, type(None))):
                value = repr(value)
            payload[key] = value

        blob = json.dumps(payload, sort_keys=True).encode('utf-8')

        return hashlib.sha256(blob).hexdigest()

    def as_dict(self):
        """
        :return:
            A plain dictionary of all settings.
        """

        return asdict(self)
