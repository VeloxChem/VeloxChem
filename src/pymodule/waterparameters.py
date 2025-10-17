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


def get_water_parameters():
    """
    Initialize water model parameters.

    :return:
        A dictionary containing water model parameters.
    """

    water_parameters = {}

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
        'O' : {
                'type': 'ow',
                'name': 'O',
                'mass': 15.994915,
                'charge': -0.8476,
                'sigma': 3.1658e-01,
                'epsilon': 6.49775e-01,  
                'equivalent_atom': 'ow_00',
                'comment': 'SPC/E water'
                },
        'H' : {
                'type': 'hw',
                'name': 'H',
                'mass': 1.007825,
                'charge': 0.4238,
                'sigma': 1.0,
                'epsilon': 0.0,
                'equivalent_atom': 'hw_00',
                'comment': 'SPC/E water'
                },
        'ref' : 'Berendsen H, Grigera J, Straatsma T. J Phys Chem. 1987;91:6269–6271.'
    }

    water_parameters['cspce'] = {
        'bonds' : {
                'type': 'harmonic',
                'force_constant': 345000,
                'equilibrium': 0.1,
                'comment': 'cSPC/E water'
                },
        'angles' : {
                'type': 'harmonic',
                'force_constant': 383,
                'equilibrium': 109.47,
                'comment': 'cSPC/E water'
                },
        'O' : {
                'type': 'ow',
                'name': 'O',
                'mass': 15.994915,
                'charge': -0.8476,
                'sigma': 3.1658e-01,
                'epsilon': 6.49775e-01,
                'equivalent_atom': 'ow_00',  
                'comment': 'cSPC/E water'
                },
        'H' : {
                'type': 'hw',
                'name': 'H',
                'mass': 1.007825,
                'charge': 0.4238,
                'sigma': 1.1658e-01,
                'epsilon': 0.64978e-01,
                'equivalent_atom': 'hw_00',
                'comment': 'cSPC/E water'
                },
        'ref' : 'T. Luchko, S. Gusarov, D. R. Roe, C. Simmerling, D. A. Case, J. Tuszynski, A. Kovalenko. J. Chem. Theory Comput. 2010 6 (3), 607-624.'
    }
    
    water_parameters['tip3p'] = {

        'bonds' : {
                'type': 'harmonic',
                'force_constant': 502416.0,
                'equilibrium': 0.09572,
                'comment': 'TIP-3P water'
                },
        'angles' : {
                'type': 'harmonic',
                'force_constant': 628.02,
                'equilibrium': 104.52,
                'comment': 'TIP-3P water'
                },
        'O' : {
                'type': 'ow',
                'name': 'O',
                'mass': 15.994915,
                'charge': -0.834,
                'sigma': 3.1507e-01, 
                'epsilon': 0.635968,
                'equivalent_atom': 'ow_00',        
                'comment': 'TIP-3P water'
                },
        'H' : {    
                'type': 'hw',
                'name': 'H',
                'mass': 1.007825,
                'charge': 0.417,
                'sigma': 1.0,
                'epsilon': 0.0,
                'equivalent_atom': 'hw_00',
                'comment': 'TIP-3P water'
                },
        'ref' : 'Jorgensen W, Chandrasekhar J, Madura J, Impey R, Klein M. J Chem Phys. 1983;79:926–935.'
    }

    water_parameters['ctip3p'] = {

        'bonds' : {
                'type': 'harmonic',
                'force_constant': 502416.0,
                'equilibrium': 0.09572,
                'comment': 'cTIP-3P water'
                },
        'angles' : {
                'type': 'harmonic',
                'force_constant': 628.02,
                'equilibrium': 104.52,
                'comment': 'cTIP-3P water'
                },
        'O' : {
                'type': 'ow',
                'name': 'O',
                'mass': 15.994915,
                'charge': -0.834,
                'sigma': 3.1507e-01, 
                'epsilon': 0.635968,
                'equivalent_atom': 'ow_00',        
                'comment': 'cTIP-3P water'
                },
        'H' : {    
                'type': 'hw',
                'name': 'H',
                'mass': 1.007825,
                'charge': 0.417,
                'sigma': 1.2363e-01,
                'epsilon': 0.63536e-01,
                'equivalent_atom': 'hw_00',
                'comment': 'cTIP-3P water'
                },
        'ref' : 'T. Luchko, S. Gusarov, D. R. Roe, C. Simmerling, D. A. Case, J. Tuszynski, A. Kovalenko. J. Chem. Theory Comput. 2010 6 (3), 607-624.'
    }

    return water_parameters
