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


def get_smd_solvent_properties():
    """
    From https://comp.chem.umn.edu/solvation/mnsddb.pdf
    """

    solvents =[
        '1,1,1-trichloroethane', '1,1,2-trichloroethane', '1,2,4-trimethylbenzene', '1,2-dibromoethane', '1,2-dichloroethane', '1,2-ethanediol', 
        '1,4-dioxane', '1-bromo-2-methylpropane', '1-bromooctane', '1-bromopentane', '1-bromopropane', '1-butanol', '1-chlorohexane', '1-chloropentane', 
        '1-chloropropane', '1-decanol', '1-fluorooctane', '1-heptanol', '1-hexanol', '1-hexene', '1-hexyne', '1-iodobutane', '1-iodohexadecane', 
        '1-iodopentane', '1-iodopropane', '1-nitropropane', '1-nonanol', '1-octanol', '1-pentanol', '1-pentene', '1-propanol', '2,2,2-trifluoroethanol', 
        '2,2,4-trimethylpentane', '2,4-dimethylpentane', '2,4-dimethylpyridine', '2,6-dimethylpyridine', '2-bromopropane', '2-butanol', '2-chlorobutane',
        '2-heptanone', '2-hexanone', '2-methoxyethanol', '2-methyl-1-propanol', '2-methyl-2-propanol', '2-methylpentane', '2-methylpyridine', 
        '2-nitropropane', '2-octanone', '2-pentanone', '2-propanol', '2-propen-1-ol', 'E-2-pentene', '3-methylpyridine', '3-pentanone', '4-heptanone', 
        '4-methyl-2-pentanone', '4-methylpyridine', '5-nonanone', 'acetic acid', 'acetone', 'acetonitrile', 'acetophenone', 'aniline', 'anisole', 
        'benzaldehyde', 'benzene', 'benzonitrile', 'benzyl alcohol', 'bromobenzene', 'bromoethane', 'bromoform', 'butanal', 'butanoic acid', 'butanone', 
        'butanonitrile', 'butyl ethanoate', 'butylamine', 'n-butylbenzene', 'sec-butylbenzene', 'tert-butylbenzene', 'carbon disulfide', 'carbon tetrachloride', 
        'chlorobenzene', 'chloroform', 'a-chlorotoluene', 'o-chlorotoluene', 'm-cresol', 'o-cresol', 'cyclohexane', 'cyclohexanone', 'cyclopentane', 'cyclopentanol', 
        'cyclopentanone', 'decalin (cis/trans mixture)', 'cis-decalin', 'n-decane', 'dibromomethane', 'dibutyl ether', 'o-dichlorobenzene', 'E-1,2-dichloroethene', 
        'Z-1,2-dichloroethene', 'dichloromethane', 'diethyl ether', 'diethyl sulfide', 'diethylamine', 'diiodomethane', 'diisopropyl ether', 'cis-1,2-dimethylcyclohexane', 
        'dimethyl disulfide', 'N,N-dimethylacetamide', 'N,N-dimethylformamide', 'dimethyl sulfoxide (DMSO)', 'diphenyl ether', 'dipropylamine', 'n-dodecane', 
        'ethanethiol', 'ethanol', 'ethyl ethanoate', 'ethyl methanoate', 'ethyl phenyl ether', 'ethylbenzene', 'fluorobenzene', 'formamide', 'formic acid', 
        'n-heptane', 'n-hexadecane', 'n-hexane', 'hexanoic acid', 'iodobenzene', 'iodoethane', 'iodomethane', 'isopropylbenzene', 'p-isopropyltoluene', 'mesitylene',
        'methanol', 'methyl benzoate', 'methyl butanoate', 'methyl ethanoate', 'methyl methanoate', 'methyl propanoate', 'N-methylaniline', 'methylcyclohexane', 
        'N-methylformamide (E/Z mixture)', 'nitrobenzene', 'nitroethane', 'nitromethane', 'o-nitrotoluene', 'n-nonane', 'n-octane', 'n-pentadecane', 'pentanal', 
        'n-pentane','pentanoic acid', 'pentyl ethanoate', 'pentylamine', 'perfluorobenzene', 'propanal', 'propanoic acid', 'propanonitrile', 'propyl ethanoate', 
        'propylamine', 'pyridine', 'tetrachloroethene', 'tetrahydrofuran', 'tetrahydrothiophene-S,S-dioxide', 'tetralin', 'thiophene', 'thiophenol', 'toluene', 
        'trans-decalin', 'tributylphosphate', 'trichloroethene', 'triethylamine', 'n-undecane', 'water', 'xylene (mixture)', 'm-xylene', 'o-xylene', 'p-xylene'
    ]
    
    epsilon_values = [
        7.0826, 7.1937, 2.3653, 4.9313, 10.125, 40.245, 2.2099, 7.7792, 5.0244, 6.269, 8.0496, 17.332, 5.9491, 6.5022, 8.3548, 7.5305, 3.89, 11.321, 12.51, 
        2.0717, 2.615, 6.173, 3.5338, 5.6973, 6.9626, 23.73, 8.5991, 9.8629, 15.13, 1.9905, 20.524, 26.726, 1.9358, 1.8939, 9.4176, 7.1735, 9.361, 15.944, 
        8.393, 11.658, 14.136, 17.2, 16.777, 12.47, 1.89, 9.9533, 25.654, 9.4678, 15.2, 19.264, 19.011, 2.051, 11.645, 16.78, 12.257, 12.887, 11.957, 10.6, 
        6.2528, 20.493, 35.688, 17.44, 6.8882, 4.2247, 18.22, 2.2706, 25.592, 12.457, 5.3954, 9.01, 4.2488, 13.45, 2.9931, 18.246, 24.291, 4.9941, 4.6178, 
        2.36, 2.3446, 2.3447, 2.6105, 2.228, 5.6968, 4.7113, 6.7175, 4.6331, 12.44, 6.76, 2.0165, 15.619, 1.9608, 16.989, 13.58, 2.196, 2.2139, 1.9846, 7.2273,
        3.0473, 9.9949, 2.14, 9.2, 8.93, 4.24, 5.723, 3.5766, 5.32, 3.38, 2.06, 9.6, 37.781, 37.219, 46.826, 3.73, 2.9112, 2.006, 6.667, 24.852, 5.9867, 8.331, 
        4.1797, 2.4339, 5.42, 108.94, 51.1, 1.9113, 2.0402, 1.8819, 2.6, 4.547, 7.6177, 6.865, 2.3712, 2.2322, 2.265, 32.613, 6.7367, 5.5607, 6.8615, 8.8377, 
        6.0777, 5.96, 2.024, 181.56, 34.809, 28.29, 36.562, 25.669, 1.9605, 1.9406, 2.0333, 10.0, 1.8371, 2.6924, 4.7297, 4.201, 2.029, 18.5, 3.44, 29.324, 
        5.5205, 4.9912, 12.978, 2.268, 7.4257, 43.962, 2.771, 2.727, 4.2728, 2.3741, 2.1781, 8.1781, 3.422, 2.3832, 1.991, 78.355, 2.3879, 2.3478, 2.5454, 2.2705
    ]
    
    alpha_values = [
        0.0, 0.13, 0.0, 0.1, 0.1, 0.58, 0.0, 0.0, 0.0, 0.0, 0.0, 0.37, 0.0, 0.0, 0.0, 0.37, 0.0, 0.37, 0.37, 0.0, 0.12, 0.0, 0.0, 0.0, 0.0, 0.0, 
        0.37, 0.37, 0.37, 0.0, 0.37, 0.57, 0.0, 0.0, 0.0, 0.0, 0.0, 0.33, 0.0, 0.0, 0.0, 0.3, 0.37, 0.31, 0.0, 0.0, 0.0, 0.0, 0.0, 0.33, 0.38, 0.0, 0.0, 0.0, 
        0.0, 0.0, 0.0, 0.0, 0.61, 0.04, 0.07, 0.0, 0.26, 0.0, 0.0, 0.0, 0.0, 0.33, 0.0, 0.0, 0.15, 0.0, 0.6, 0.0, 0.0, 0.0, 0.16, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
        0.15, 0.0, 0.0, 0.57, 0.52, 0.0, 0.0, 0.0, 0.32, 0.0, 0.0, 0.0, 0.0, 0.1, 0.0, 0.0, 0.09, 0.11, 0.1, 0.0, 0.0, 0.08, 0.05, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
        0.0, 0.08, 0.0, 0.0, 0.37, 0.0, 0.0, 0.0, 0.0, 0.0, 0.62, 0.75, 0.0, 0.0, 0.0, 0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.43, 0.0, 0.0, 0.0, 0.0, 0.0, 0.17, 
        0.0, 0.4, 0.0, 0.02, 0.06, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.6, 0.0, 0.16, 0.0, 0.0, 0.6, 0.02, 0.0, 0.16, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.09, 0.0, 0.0, 
        0.0, 0.08, 0.0, 0.0, 0.82, 0.0, 0.0, 0.0, 0.0
    ]
    
    
    beta_values = [
        0.09, 0.13, 0.19, 0.17, 0.11, 0.78, 0.64, 0.12, 0.12, 0.12, 0.12, 0.48, 0.1, 0.1, 0.1, 0.48, 0.1, 0.48, 0.48, 0.07, 0.1, 0.15, 0.15, 0.15, 0.15, 0.31, 
        0.48, 0.48, 0.48, 0.07, 0.48, 0.25, 0.0, 0.0, 0.63, 0.63, 0.14, 0.56, 0.12, 0.51, 0.51, 0.84, 0.48, 0.6, 0.0, 0.58, 0.33, 0.51, 0.51, 0.56, 0.48, 0.07, 
        0.54, 0.51, 0.51, 0.51, 0.54, 0.51, 0.44, 0.49, 0.32, 0.48, 0.41, 0.29, 0.39, 0.14, 0.33, 0.56, 0.09, 0.12, 0.06, 0.45, 0.45, 0.51, 0.36, 0.45, 0.61, 
        0.15, 0.16, 0.16, 0.07, 0.0, 0.07, 0.02, 0.33, 0.07, 0.34, 0.3, 0.0, 0.56, 0.0, 0.56, 0.52, 0.0, 0.0, 0.0, 0.1, 0.45, 0.04, 0.05, 0.05, 0.05, 0.41, 0.32, 
        0.69, 0.23, 0.41, 0.0, 0.28, 0.78, 0.74, 0.88, 0.2, 0.69, 0.0, 0.24, 0.48, 0.45, 0.38, 0.32, 0.15, 0.1, 0.6, 0.38, 0.0, 0.0, 0.0, 0.45, 0.12, 0.15, 0.13, 
        0.16, 0.19, 0.19, 0.47, 0.46, 0.45, 0.45, 0.38, 0.45, 0.43, 0.0, 0.55, 0.28, 0.33, 0.31, 0.27, 0.0, 0.0, 0.0, 0.45, 0.0, 0.45, 0.45, 0.61, 0.0, 0.45, 
        0.45, 0.36, 0.45, 0.61, 0.52, 0.0, 0.48, 0.88, 0.19, 0.15, 0.16, 0.14, 0.0, 1.21, 0.03, 0.79, 0.0, 0.35, 0.16, 0.16, 0.16, 0.16
    ]
    
    # the macroscopic surface tension of the solvent at air/solvent interface at 298.15 K
    # cal mol-1 Å-2 (1 dyn/cm = 1.439 32 cal mol-1Å-2)
    gamma_values = [
        36.24, 48.97, 42.03, 56.93, 45.86, 69.07, 47.14, 34.69, 41.28, 38.7, 36.36, 35.88, 37.03, 35.12, 30.66, 41.04, 33.92, 38.5, 37.15, 25.76, 28.79, 40.65, 
        46.48, 41.56, 41.45, 43.32, 40.14, 39.01, 36.5, 22.24, 33.57, 42.02, 26.38, 25.42, 46.86, 44.64, 33.46, 32.44, 31.1, 37.6, 36.63, 44.39, 32.38, 28.73, 
        24.3, 47.5, 42.16, 37.29, 33.46, 30.13, 36.39, 23.62, 49.61, 35.61, 35.98, 33.83, 50.17, 37.83, 39.01, 33.77, 41.25, 56.19, 60.62, 50.52, 54.69, 40.62, 
        55.83, 52.96, 50.72, 34.0, 64.58, 35.06, 37.49, 34.5, 38.75, 35.81, 33.74, 41.33, 40.35, 39.78, 45.45, 38.04, 47.48, 38.39, 53.04, 47.43, 51.37, 53.11, 
        35.48, 49.76, 31.49, 46.8, 47.21, 43.82, 45.45, 33.64, 56.21, 35.98, 52.72, 37.13, 39.8, 39.15, 23.96, 35.36, 28.57, 95.25, 24.86, 36.28, 48.06, 47.62, 
        49.56, 61.78, 38.5, 32.11, 35.85, 33.22, 31.62, 33.67, 33.36, 46.65, 41.38, 38.37, 82.08, 53.44, 28.28, 38.93, 25.75, 39.65, 55.72, 40.96, 43.67, 39.85, 
        38.34, 39.65, 31.77, 53.5, 35.44, 35.59, 35.06, 35.18, 53.11, 33.52, 55.44, 57.54, 46.25, 52.58, 59.12, 32.21, 30.43, 38.34, 36.62, 22.3, 38.4, 36.23, 
        35.54, 31.74, 32.48, 37.71, 38.5, 34.26, 31.31, 52.62, 45.19, 39.44, 87.49, 47.74, 44.16, 55.24, 40.2, 42.19, 27.55, 41.45, 29.1, 34.85, 0.0, 41.38, 
        40.98, 42.83, 40.32
    ]
    
    
    phi_values = [
        0.0, 0.0, 0.667, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.625, 0.625, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.714, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.714, 0.0, 0.0, 0.0, 0.714, 
        0.0, 0.0, 0.0, 0.0, 0.667, 0.857, 0.75, 0.857, 1.0, 0.75, 0.75, 0.857, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.6, 0.6, 0.6, 0.0, 0.0, 0.857, 0.0, 0.75, 
        0.75, 0.75, 0.75, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.75, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.923, 0.0, 0.0, 
        0.0, 0.0, 0.0, 0.0, 0.667, 0.75, 0.857, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.857, 0.0, 0.0, 0.667, 0.6, 0.667, 0.0, 0.6, 0.0, 0.0, 0.0, 0.0, 0.75, 0.0, 0.0, 
        0.667, 0.0, 0.0, 0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.833, 0.0, 0.0, 0.0, 0.6, 0.8, 0.857, 0.857, 0.0, 0.0, 0.0, 
        0.0, 0.0, 0.0, 0.75, 0.75, 0.75, 0.75
    ]
    
    
    psi_values = [
        0.6, 0.6, 0.0, 0.5, 0.5, 0.0, 0.0, 0.2, 0.111, 0.167, 0.25, 0.0, 0.143, 0.167, 0.25, 0.0, 0.111, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
        0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.25, 0.0, 0.2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.143, 0.333, 0.75, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.8, 0.143, 0.75, 0.125, 
        0.125, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.667, 0.0, 0.25, 0.5, 0.5, 0.667, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.143, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.667, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.6, 0.0, 0.0, 
        0.0, 0.0, 0.0, 0.0, 0.0
    ]
    
    # At 293 K
    refractive_idx = [
        1.4379, 1.4714, 1.5048, 1.5387, 1.4448, 1.4318, 1.4224, 1.4348, 1.4524, 1.4447, 1.4343, 1.3993, 1.4199, 1.4127, 1.3879, 1.4372, 1.3935, 
        1.4249, 1.4178, 1.3837, 1.3989, 1.5001, 1.4806, 1.4959, 1.5058, 1.4018, 1.4333, 1.4295, 1.4101, 1.3715, 1.385, 1.2907, 1.3915, 1.3815, 1.501, 1.4953, 
        1.4251, 1.3978, 1.3971, 1.4088, 1.4007, 1.4024, 1.3955, 1.3878, 1.3715, 1.4957, 1.3944, 1.4151, 1.3895, 1.3776, 1.4135, 1.3793, 1.504, 1.3924, 1.4069, 
        1.3962, 1.5037, 1.4195, 1.372, 1.3588, 1.3442, 1.5372, 1.5863, 1.5174, 1.5463, 1.5011, 1.5289, 1.5396, 1.5597, 1.4239, 1.6005, 1.3843, 1.398, 1.3788, 
        1.3842, 1.3941, 1.4031, 1.4898, 1.4895, 1.4927, 1.6319, 1.4601, 1.5241, 1.4459, 1.5391, 1.5268, 1.5438, 1.5361, 1.4266, 1.4507, 1.4065, 1.453, 1.4366, 
        1.4753, 1.481, 1.4102, 1.542, 1.3992, 1.5515, 1.4454, 1.449, 1.4242, 1.3526, 1.443, 1.3864, 1.7425, 1.3679, 1.436, 1.5289, 1.438, 1.4305, 1.4783, 1.5787, 
        1.405, 1.4216, 1.431, 1.3611, 1.3723, 1.3599, 1.5076, 1.4959, 1.4684, 1.4472, 1.3714, 1.3878, 1.4345, 1.3749, 1.4163, 1.62, 1.5133, 1.538, 1.4915, 1.4909, 
        1.4994, 1.3288, 1.5164, 1.3878, 1.3614, 1.3433, 1.3775, 1.5684, 1.4231, 1.4319, 1.5562, 1.3917, 1.3817, 1.545, 1.4054, 1.3974, 1.4315, 1.3944, 1.3575, 
        1.4085, 1.4023, 1.448, 1.3777, 1.3636, 1.3869, 1.3655, 1.3842, 1.387, 1.5095, 1.5053, 1.405, 1.4833, 1.5413, 1.5289, 1.5893, 1.4961, 1.4695, 1.4224, 
        1.4773, 1.401, 1.4398, 1.3328, 1.4995, 1.4972, 1.5055, 1.4958
    ]


    smd_solvent_parameters = {}

    for idx, (epsilon,alpha,beta,gamma,phi,psi,n) in enumerate(zip(epsilon_values,alpha_values,beta_values,gamma_values,phi_values,psi_values,refractive_idx)):
        key = solvents[idx]

        smd_solvent_parameters[key] = {
            'epsilon': epsilon, 
            'alpha':alpha, 
            'beta':beta, 
            'gamma':gamma,
            'phi':phi,
            'psi':psi,
            'n':n, 
        }


    return smd_solvent_parameters

def get_sigma_properties():

    # sigma_[n], sigma_[alpha], sigma_[beta]
    atom_dep_param = {
        ('C',): [58.10, 48.10, 32.87],
        ('H', 'C'): [-36.37, 0, 0],
        ('C', 'C'): [-62.05, 0, 0],
        ('O',): [-17.56, 193.06, -43.79],
        ('H', 'O'): [-19.39, 0, 0],
        ('O', 'C'): [-15.70, 95.99, 0],
        ('O', 'O'): [0, 0, -128.16],
        ('N',): [32.62, 0, 0],
        ('C', 'N'): [-99.76, 152.20, 0], 
        ('N', 'C'): [0, -41.00, 0],
        ('O', 'N'): [0, 0, 79.13],
        ('Cl',): [-24.31, 0, 0],
        ('Br',): [-35.42, 0, 0],
        ('S',): [-33.17, 0, 0],
        ('Si',): [-18.04, 0, 0]
    }

    # sigma_[water]
    water_param = {
        ('H',): 48.69,
        ('C',): 129.74,
        ('H', 'C'): -60.77,
        ('C', 'C'): -72.95,
        ('O', 'C'): 68.69,
        ('N', 'C'): -48.22,
        ('N', 'c3'): 84.10,  
        ('O', 'N'): 121.98,
        ('F',): 38.18,
        ('Cl',): 9.82,
        ('Br',): -8.72,
        ('S',): -9.10,
        ('O', 'P'): 68.85
    }

    return atom_dep_param, water_param

def get_rzz_parameters():
    # r_ZZ', dr_ZZ'
    r_ZZ_param = {
        ('H', 'C'): [1.55, 0.3],
        ('H', 'O'): [1.55, 0.3],
        ('C', 'H'): [1.55, 0.3],
        ('C', 'C'): [1.84, 0.3],
        ('C', 'N'): [1.84, 0.3], 
        ('C', 'O'): [1.84, 0.3],
        ('C', 'F'): [1.84, 0.3],
        ('C', 'P'): [2.2, 0.3],
        ('C', 'S'): [2.2, 0.3],
        ('C', 'Cl'): [2.1, 0.3],
        ('C', 'Br'): [2.3, 0.3],
        ('C', 'I'): [2.6, 0.3],
        ('N', 'C'): [1.84, 0.3],
        ('N', 'c3'): [1.225, 0.065],  
        ('O', 'C'): [1.33, 0.1],
        ('O', 'N'): [1.5, 0.3],
        ('O', 'O'): [1.8, 0.3],
        ('O', 'P'): [2.1, 0.3]
    }

    return r_ZZ_param