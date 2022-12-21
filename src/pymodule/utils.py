#
#                           VELOXCHEM 1.0-RC3
#         ----------------------------------------------------
#                     An Electronic Structure Code
#
#  Copyright © 2018-2022 by VeloxChem developers. All rights reserved.
#  Contact: https://veloxchem.org/contact
#
#  SPDX-License-Identifier: LGPL-3.0-or-later
#
#  This file is part of VeloxChem.
#
#  VeloxChem is free software: you can redistribute it and/or modify it under
#  the terms of the GNU Lesser General Public License as published by the Free
#  Software Foundation, either version 3 of the License, or (at your option)
#  any later version.
#
#  VeloxChem is distributed in the hope that it will be useful, but WITHOUT
#  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
#  License for more details.
#
#  You should have received a copy of the GNU Lesser General Public License
#  along with VeloxChem. If not, see <https://www.gnu.org/licenses/>.

import numpy as np
import math

from .veloxchemlib import (fine_structure_constant, hartree_in_ev,
                           extinction_coefficient_from_beta,
                           rotatory_strength_in_cgs)
from .errorhandler import assert_msg_critical


def lorentzian_absorption_spectrum(energy_unit, exc_ene, osc_str, e_min, e_max,
                                   e_step, gamma):
    """
    Broadens absorption stick spectrum.

    :param energy_unit:
        The unit of excitation energies.
    :param exc_ene:
        Excitation energies in energy_unit.
    :param osc_str:
        Oscillator strengths.
    :param e_min:
        Minimal excitation energy in energy_unit in broadened spectrum.
    :param e_max:
        Maximum excitation energy in energy_unit in broadened spectrum.
    :param e_step:
        Step size of excitation energy in energy_unit in broadened spectrum.
    :param gamma:
        The broadening parameter in energy_unit.

    :return:
        The excitation energies in energy_unit and sigma(w) in a.u.
    """

    assert_msg_critical(energy_unit.lower() in ['ev', 'au'],
                        'lorentzian_absorption_spectrum: Invalid energy_unit')

    exc_ene_copy = exc_ene.copy()

    x_i = np.arange(e_min, e_max + e_step / 100.0, e_step, dtype='float64')
    y_i = np.zeros_like(x_i)

    if energy_unit.lower() == 'ev':
        x_i /= hartree_in_ev()
        gamma /= hartree_in_ev()
        exc_ene_copy /= hartree_in_ev()

    factor = 2.0 * math.pi * fine_structure_constant()

    for i in range(x_i.size):
        for s in range(exc_ene_copy.size):
            y_i[i] += factor * gamma / (
                (x_i[i] - exc_ene_copy[s])**2 + gamma**2) * osc_str[s]

    if energy_unit.lower() == 'ev':
        x_i *= hartree_in_ev()

    return x_i, y_i


def lorentzian_ecd_spectrum(energy_unit, exc_ene, rot_str, e_min, e_max, e_step,
                            gamma):
    """
    Broadens ECD stick spectrum.

    :param energy_unit:
        The unit of excitation energies.
    :param exc_ene:
        Excitation energies in energy_unit
    :param rot_str:
        Rotatory strengths in 10**(-40) cgs unit.
    :param e_min:
        Minimal excitation energy in energy_unit in broadened spectrum.
    :param e_max:
        Maximum excitation energy in energy_unit in broadened spectrum.
    :param e_step:
        Step size of excitation energy in energy_unit in broadened spectrum.
    :param gamma:
        The broadening parameter in energy_unit.

    :return:
        The excitation energies in energy_unit and Delta_epsilon in L mol^-1 cm^-1
    """

    assert_msg_critical(energy_unit.lower() in ['ev', 'au'],
                        'lorentzian_ecd_spectrum: Invalid energy_unit')

    exc_ene_copy = exc_ene.copy()

    x_i = np.arange(e_min, e_max + e_step / 100.0, e_step, dtype='float64')
    y_i = np.zeros_like(x_i)

    if energy_unit.lower() == 'ev':
        x_i /= hartree_in_ev()
        gamma /= hartree_in_ev()
        exc_ene_copy /= hartree_in_ev()

    f = 1.0 / rotatory_strength_in_cgs()  # convert rot_str to a.u.
    f *= extinction_coefficient_from_beta() / 3.0

    for i in range(x_i.size):
        for s in range(exc_ene_copy.size):
            y_i[i] += f * gamma / ((x_i[i] - exc_ene_copy[s])**2 +
                                   gamma**2) * exc_ene_copy[s] * rot_str[s]

    if energy_unit.lower() == 'ev':
        x_i *= hartree_in_ev()

    return x_i, y_i
