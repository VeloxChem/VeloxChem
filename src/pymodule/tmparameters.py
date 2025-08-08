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


def get_tm_parameters():
    """
    Initialize transition metal (TM) parameters. Reference: F. Šebesta, V. Sláma, J. Melcr, 
    Z. Futera, and J. V. Burda. J. Chem. Theory Comput. 2016 12 (8), 3681-3688.

    :return:
        A dictionary containing TM parameters. When multiple oxidation states available in the reference, 
        parameters are chosen for the lower oxidation state. 
    """

    tm_elements = [
        'Sc','Ti','V','Cr','Mn',
        'Fe','Co','Ni','Cu','Zn',
        'Ru','Rh','Pt','Hg'
    ]   
    
    sigma_values = [
        3.330, 2.869, 2.767, 2.734, 2.904,
        2.969, 2.805, 2.720, 2.668, 2.847, 
        2.939, 2.715, 2.670, 2.801
    ]
    
    epsilon_values = [
        0.110, 0.917, 1.904, 1.518, 1.089,
        0.698, 1.196, 2.650, 2.148, 1.095, 
        0.418, 3.946, 5.073, 1.955
    ]

    tm_parameters = {}

    for idx, (s, e) in enumerate(zip(sigma_values, epsilon_values)):
        key = tm_elements[idx]

        sigma = 0.1 * s * 2.0
        epsilon = e * 4.184  

        tm_parameters[key] = {'sigma': sigma, 'epsilon': epsilon}


    return tm_parameters
