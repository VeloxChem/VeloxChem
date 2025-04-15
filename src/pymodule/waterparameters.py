#                              VELOXCHEM
#         ----------------------------------------------------
#                     An Electronic Structure Code
#
#  Copyright © 2018-2024 by VeloxChem developers. All rights reserved.
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


def get_water_parameters():
    """
    Initialize water model parameters. Reference: T. Luchko, S. Gusarov, D. R. Roe, C. Simmerling, 
    D. A. Case, J. Tuszynski, A. Kovalenko. J. Chem. Theory Comput. 2010 6 (3), 607-624.

    :return:
        A dictionary containing water model parameters.
    """


    water_parameters = {}


    # TODO: Add more water models
    water_parameters['spce'] = {
        'bonds' : {
                'type': 'harmonic',
                'force_constant': 345000,
                'equilibrium': 0.1,
                'comment': 'SPC/E water'
                },
        'angles' : {
                'type': 'harmonic',
                'force_constant': 383,
                'equilibrium': 109.47,
                'comment': 'SPC/E water'
                },
        'ow' : {
                'type': 'ow',
                'name': 'O',
                'mass': 15.994915,
                'charge': -0.8476,
                'sigma': 3.1658e-01,
                'epsilon': 6.49775e-01,  
                'equivalent_atom': 'SPC/E water'
                },
        'hw' : {
                'type': 'hw',
                'name': 'H1',
                'mass': 1.007825,
                'charge': 0.4238,
                'sigma': 1.1658e-01,
                'epsilon': 0.64978e-01,
                'equivalent_atom': 'SPC/E water'
                }
    }

    water_parameters['tip3p'] = {

        'bonds' : {
                'type': 'harmonic',
                'force_constant': 462750.4,
                'equilibrium': 0.09572,
                'comment': 'TIP-3P water'
                },
        'angles' : {
                'type': 'harmonic',
                'force_constant': 836.8,
                'equilibrium': 104.52,
                'comment': 'TIP-3P water'
                },
        'ow' : {
                'type': 'ow',
                'name': 'O',
                'mass': 15.994915,
                'charge': -0.834,
                'sigma': 3.1507e-01, 
                'epsilon': 0.635968,        
                'equivalent_atom': 'TIP-3P water'
                },
        'hw' : {    
                'type': 'hw',
                'name': 'H1',
                'mass': 1.007825,
                'charge': 0.417,
                'sigma': 1.2363e-01,
                'epsilon': 0.63536e-01,
                'equivalent_atom': 'TIP-3P water'
                }
    }


    return water_parameters