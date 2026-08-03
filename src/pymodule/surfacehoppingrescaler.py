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
Momentum rescaling for surface hopping.

This module is pure NumPy: it operates on masses, velocities and gradients in
atomic units and has no VeloxChem or OpenMM dependencies.

Units throughout:

    * masses      : electron masses
    * velocities  : bohr / atomic time unit
    * energies    : Hartree
"""

from dataclasses import dataclass
import numpy as np


@dataclass
class RescalingResult:
    """
    Outcome of a momentum-rescaling attempt.

    :param accepted:
        True when the hop is energetically allowed and the velocities were
        rescaled successfully.
    :param frustrated:
        True when the hop was rejected because insufficient kinetic energy is
        available along the rescaling direction.
    :param velocities:
        The rescaled velocities, or the unchanged input velocities when the
        hop was frustrated.
    :param kinetic_before:
        Total kinetic energy of the rescaling atom set before the hop.
    :param kinetic_after:
        Total kinetic energy of the rescaling atom set after the hop.
    :param scaling_factor:
        The uniform scaling factor k for ``global_qm`` mode, else None.
    :param gamma:
        The solved coefficient gamma for ``gap_gradient`` mode, else None.
    :param energy_residual:
        Residual of the energy-conservation check
        (K_after + E_b) - (K_before + E_a), in Hartree.
    :param message:
        Human-readable explanation, primarily for frustrated hops.
    """

    accepted: bool
    frustrated: bool
    velocities: np.ndarray
    kinetic_before: float
    kinetic_after: float
    scaling_factor: float = None
    gamma: float = None
    energy_residual: float = 0.0
    message: str = ''


def kinetic_energy(masses, velocities):
    """
    Evaluates the kinetic energy of an atom set.

    :param masses:
        Array of shape (n,) with masses in electron masses.
    :param velocities:
        Array of shape (n, 3) with velocities in bohr per atomic time unit.

    :return:
        The kinetic energy in Hartree.
    """

    masses = np.asarray(masses, dtype=float)
    velocities = np.asarray(velocities, dtype=float)

    return float(0.5 * np.sum(masses[:, None] * velocities**2))


def center_of_mass_velocity(masses, velocities):
    """
    :param masses:
        Array of shape (n,) with masses in electron masses.
    :param velocities:
        Array of shape (n, 3) with velocities.

    :return:
        The center-of-mass velocity, shape (3,).
    """

    masses = np.asarray(masses, dtype=float)
    velocities = np.asarray(velocities, dtype=float)
    total_mass = float(np.sum(masses))

    if total_mass <= 0.0:
        return np.zeros(3)

    return np.sum(masses[:, None] * velocities, axis=0) / total_mass


def angular_velocity(masses, positions, velocities):
    """
    Evaluates the rigid-body angular velocity of an atom set about its center
    of mass.

    Solves ``I omega = L`` with the inertia tensor I and angular momentum L,
    using a pseudo-inverse so that linear and single-atom geometries, whose
    inertia tensors are singular, are handled gracefully.

    :param masses:
        Array of shape (n,) with masses in electron masses.
    :param positions:
        Array of shape (n, 3) with positions in bohr.
    :param velocities:
        Array of shape (n, 3) with velocities in bohr per atomic time unit.

    :return:
        The angular velocity vector, shape (3,).
    """

    masses = np.asarray(masses, dtype=float)
    positions = np.asarray(positions, dtype=float)
    velocities = np.asarray(velocities, dtype=float)

    total_mass = float(np.sum(masses))
    if total_mass <= 0.0 or masses.size < 2:
        return np.zeros(3)

    com = np.sum(masses[:, None] * positions, axis=0) / total_mass
    rel = positions - com

    angular_momentum = np.sum(masses[:, None] * np.cross(rel, velocities),
                              axis=0)

    inertia = np.zeros((3, 3))
    for m, r in zip(masses, rel):
        inertia += m * (np.dot(r, r) * np.eye(3) - np.outer(r, r))

    # Pseudo-inverse tolerates linear/atomic geometries with singular inertia.
    omega = np.linalg.pinv(inertia, rcond=1.0e-10) @ angular_momentum

    if not np.all(np.isfinite(omega)):
        return np.zeros(3)

    return omega


def decompose_velocities(masses,
                         positions,
                         velocities,
                         remove_translation=True,
                         remove_rotation=False):
    """
    Splits velocities into rigid-body translation, rigid-body rotation and an
    internal remainder.

    :param masses:
        Array of shape (n,) with masses in electron masses.
    :param positions:
        Array of shape (n, 3) with positions in bohr.
    :param velocities:
        Array of shape (n, 3) with velocities.
    :param remove_translation:
        Whether to split off the center-of-mass translation.
    :param remove_rotation:
        Whether to split off the rigid-body rotation.

    :return:
        The tuple (v_translation, v_rotation, v_internal), each of shape
        (n, 3) and summing to ``velocities``.
    """

    masses = np.asarray(masses, dtype=float)
    positions = np.asarray(positions, dtype=float)
    velocities = np.asarray(velocities, dtype=float)

    v_trans = np.zeros_like(velocities)
    v_rot = np.zeros_like(velocities)

    if remove_translation:
        v_trans = np.tile(center_of_mass_velocity(masses, velocities),
                          (velocities.shape[0], 1))

    residual = velocities - v_trans

    if remove_rotation:
        total_mass = float(np.sum(masses))
        com = np.sum(masses[:, None] * positions, axis=0) / total_mass
        rel = positions - com
        omega = angular_velocity(masses, positions, residual)
        v_rot = np.cross(np.tile(omega, (velocities.shape[0], 1)), rel)

    v_internal = velocities - v_trans - v_rot

    return v_trans, v_rot, v_internal


class SurfaceHopMomentumRescaler:
    """
    Adjusts nuclear velocities so that a surface hop conserves total energy.

    Two modes are supported.

    ``global_qm``
        Uniform scaling of the internal velocity component of the rescaling
        atom set.  With ``K'= K - Delta E`` the scaling factor is

        .. math:: k = \\sqrt{1 - \\Delta E / K}

        Upward hops with ``K < Delta E`` are frustrated.

    ``gap_gradient``
        Scaling along the gap-gradient direction
        ``u = grad(E_b) - grad(E_a)``.  With ``P' = P + gamma u`` energy
        conservation gives ``A gamma^2 + B gamma + Delta E = 0`` with

        .. math::

            A = \\frac{1}{2}\\sum_I \\frac{u_I \\cdot u_I}{M_I}, \\qquad
            B = \\sum_I \\frac{P_I \\cdot u_I}{M_I}

        A negative discriminant frustrates the hop; otherwise the root of
        smallest magnitude is taken.

    :param mode:
        ``'global_qm'`` or ``'gap_gradient'``.
    :param remove_translation:
        Whether center-of-mass translation is excluded from rescaling and
        restored afterwards (``global_qm`` only).
    :param remove_rotation:
        Whether rigid-body rotation is excluded from rescaling and restored
        afterwards (``global_qm`` only).
    :param energy_tolerance:
        Maximum acceptable magnitude of the energy-conservation residual, in
        Hartree.
    """

    VALID_MODES = ('global_qm', 'gap_gradient')

    def __init__(self,
                 mode='global_qm',
                 remove_translation=True,
                 remove_rotation=False,
                 energy_tolerance=1.0e-8):

        mode = str(mode).lower()

        if mode not in self.VALID_MODES:
            valid = ', '.join(self.VALID_MODES)
            raise ValueError(
                f'SurfaceHopMomentumRescaler: unknown mode {mode}; '
                f'valid modes are {valid}.')

        self.mode = mode
        self.remove_translation = bool(remove_translation)
        self.remove_rotation = bool(remove_rotation)
        self.energy_tolerance = float(energy_tolerance)

    def requires_target_gradient(self):
        """
        :return:
            True when the configured mode needs the target-state gradient, so
            that the controller can avoid the extra gradient evaluation
            whenever it is not required.
        """

        return self.mode == 'gap_gradient'

    def rescale(self,
                masses,
                positions,
                velocities,
                delta_energy,
                direction=None):
        """
        Attempts to rescale velocities for a hop with energy change
        ``delta_energy = E_target - E_source``.

        :param masses:
            Array of shape (n,) with masses in electron masses.
        :param positions:
            Array of shape (n, 3) with positions in bohr.
        :param velocities:
            Array of shape (n, 3) with velocities in bohr per atomic time
            unit.
        :param delta_energy:
            ``E_target - E_source`` in Hartree.  Positive for upward hops.
        :param direction:
            Array of shape (n, 3) with the gap-gradient direction
            ``grad(E_b) - grad(E_a)``, in Hartree per bohr.  Required for
            ``gap_gradient`` mode and ignored otherwise.

        :return:
            A :class:`RescalingResult`.
        """

        masses = np.asarray(masses, dtype=float)
        positions = np.asarray(positions, dtype=float)
        velocities = np.asarray(velocities, dtype=float)
        delta_energy = float(delta_energy)

        if masses.ndim != 1 or velocities.shape != (masses.size, 3):
            raise ValueError(
                'SurfaceHopMomentumRescaler: inconsistent masses/velocities '
                'shapes.')
        if not np.all(np.isfinite(velocities)):
            raise ValueError(
                'SurfaceHopMomentumRescaler: nonfinite input velocities.')
        if not np.isfinite(delta_energy):
            raise ValueError(
                'SurfaceHopMomentumRescaler: nonfinite energy difference.')
        if np.any(masses <= 0.0):
            raise ValueError(
                'SurfaceHopMomentumRescaler: nonpositive atomic mass.')

        if self.mode == 'global_qm':
            return self._rescale_global(masses, positions, velocities,
                                        delta_energy)

        return self._rescale_directional(masses, velocities, delta_energy,
                                         direction)

    def _rescale_global(self, masses, positions, velocities, delta_energy):
        """
        Uniform scaling of the internal velocity component.  See
        :class:`SurfaceHopMomentumRescaler`.
        """

        k_before = kinetic_energy(masses, velocities)

        v_trans, v_rot, v_internal = decompose_velocities(
            masses,
            positions,
            velocities,
            remove_translation=self.remove_translation,
            remove_rotation=self.remove_rotation)

        k_internal = kinetic_energy(masses, v_internal)

        if k_internal <= 0.0:
            return RescalingResult(
                accepted=False,
                frustrated=True,
                velocities=velocities,
                kinetic_before=k_before,
                kinetic_after=k_before,
                message='no internal kinetic energy available for rescaling')

        if k_internal < delta_energy:
            return RescalingResult(
                accepted=False,
                frustrated=True,
                velocities=velocities,
                kinetic_before=k_before,
                kinetic_after=k_before,
                message=(f'frustrated hop: internal kinetic energy '
                         f'{k_internal:.6e} < energy demand '
                         f'{delta_energy:.6e} Hartree'))

        factor = np.sqrt(1.0 - delta_energy / k_internal)

        if not np.isfinite(factor):
            return RescalingResult(
                accepted=False,
                frustrated=True,
                velocities=velocities,
                kinetic_before=k_before,
                kinetic_after=k_before,
                message='nonfinite scaling factor')

        # Translation and rotation are restored untouched.
        new_velocities = v_trans + v_rot + factor * v_internal

        k_after = kinetic_energy(masses, new_velocities)

        # Energy check is done on the rescaled (internal) component, which is
        # the part that absorbs the electronic energy change.
        residual = (kinetic_energy(masses, factor * v_internal) +
                    delta_energy) - k_internal

        if abs(residual) > self.energy_tolerance:
            return RescalingResult(
                accepted=False,
                frustrated=False,
                velocities=velocities,
                kinetic_before=k_before,
                kinetic_after=k_before,
                scaling_factor=float(factor),
                energy_residual=float(residual),
                message=(f'energy conservation residual {residual:.3e} '
                         f'exceeds tolerance {self.energy_tolerance:.3e}'))

        return RescalingResult(accepted=True,
                               frustrated=False,
                               velocities=new_velocities,
                               kinetic_before=k_before,
                               kinetic_after=k_after,
                               scaling_factor=float(factor),
                               energy_residual=float(residual),
                               message='')

    def _rescale_directional(self, masses, velocities, delta_energy,
                             direction):
        """
        Scaling along the gap-gradient direction.  See
        :class:`SurfaceHopMomentumRescaler`.
        """

        if direction is None:
            raise ValueError(
                'SurfaceHopMomentumRescaler: gap_gradient mode requires a '
                'direction vector.')

        direction = np.asarray(direction, dtype=float)

        if direction.shape != velocities.shape:
            raise ValueError(
                'SurfaceHopMomentumRescaler: direction shape does not match '
                'velocities.')
        if not np.all(np.isfinite(direction)):
            raise ValueError(
                'SurfaceHopMomentumRescaler: nonfinite direction vector.')

        k_before = kinetic_energy(masses, velocities)
        momenta = masses[:, None] * velocities

        a_coef = 0.5 * float(np.sum(np.sum(direction * direction, axis=1) /
                                    masses))
        b_coef = float(np.sum(np.sum(momenta * direction, axis=1) / masses))

        if a_coef <= 0.0:
            return RescalingResult(
                accepted=False,
                frustrated=True,
                velocities=velocities,
                kinetic_before=k_before,
                kinetic_after=k_before,
                message='vanishing gap-gradient direction')

        discriminant = b_coef * b_coef - 4.0 * a_coef * delta_energy

        if discriminant < 0.0:
            return RescalingResult(
                accepted=False,
                frustrated=True,
                velocities=velocities,
                kinetic_before=k_before,
                kinetic_after=k_before,
                message=(f'frustrated hop: negative discriminant '
                         f'{discriminant:.6e}'))

        sqrt_disc = np.sqrt(discriminant)
        roots = ((-b_coef + sqrt_disc) / (2.0 * a_coef),
                 (-b_coef - sqrt_disc) / (2.0 * a_coef))
        gamma = float(min(roots, key=abs))

        if not np.isfinite(gamma):
            return RescalingResult(
                accepted=False,
                frustrated=True,
                velocities=velocities,
                kinetic_before=k_before,
                kinetic_after=k_before,
                message='nonfinite gamma')

        new_momenta = momenta + gamma * direction
        new_velocities = new_momenta / masses[:, None]

        k_after = kinetic_energy(masses, new_velocities)
        residual = (k_after + delta_energy) - k_before

        if abs(residual) > self.energy_tolerance:
            return RescalingResult(
                accepted=False,
                frustrated=False,
                velocities=velocities,
                kinetic_before=k_before,
                kinetic_after=k_before,
                gamma=gamma,
                energy_residual=float(residual),
                message=(f'energy conservation residual {residual:.3e} '
                         f'exceeds tolerance {self.energy_tolerance:.3e}'))

        return RescalingResult(accepted=True,
                               frustrated=False,
                               velocities=new_velocities,
                               kinetic_before=k_before,
                               kinetic_after=k_after,
                               gamma=gamma,
                               energy_residual=float(residual),
                               message='')
