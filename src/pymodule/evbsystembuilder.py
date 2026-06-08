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

from mpi4py import MPI
import numpy as np
import networkx as nx
from networkx.algorithms.isomorphism import GraphMatcher
from networkx.algorithms.isomorphism import categorical_node_match
from pathlib import Path
import copy
import math
import sys
from enum import Enum, auto

from .veloxchemlib import hartree_in_kcalpermol, bohr_in_angstrom, Point
from .veloxchemlib import mpi_master
from .outputstream import OutputStream
from .mmforcefieldgenerator import MMForceFieldGenerator
from .solvationbuilder import SolvationBuilder
from .molecule import Molecule
from .errorhandler import assert_msg_critical
from .mathutils import safe_arccos
from .waterparameters import get_water_parameters

try:
    import openmm as mm
    import openmm.app as mmapp
    import openmm.unit as mmunit
except ImportError:
    pass

# ---------------------------------------------------------------------------
# GB implicit-solvent constants — ported from openmm.app.internal.customgbforces
# so the public CustomGBForce API is the only OpenMM dependency.
# ---------------------------------------------------------------------------

# Neck lookup tables (GBn / GBn2).  Original values are in Angstrom; rescaled
# here to nanometers (d0 /= 10, m0 *= 10) to match OpenMM's internal units.
_GB_D0 = [
    2.26685, 2.32548, 2.38397, 2.44235, 2.50057, 2.55867, 2.61663, 2.67444,
    2.73212, 2.78965, 2.84705, 2.90430, 2.96141, 3.01840, 3.07524, 3.13196,
    3.18854, 3.24498, 3.30132, 3.35752, 3.41360,
    2.31191, 2.37017, 2.42830, 2.48632, 2.54420, 2.60197, 2.65961, 2.71711,
    2.77449, 2.83175, 2.88887, 2.94586, 3.00273, 3.05948, 3.11610, 3.17260,
    3.22897, 3.28522, 3.34136, 3.39738, 3.45072,
    2.35759, 2.41549, 2.47329, 2.53097, 2.58854, 2.64600, 2.70333, 2.76056,
    2.81766, 2.87465, 2.93152, 2.98827, 3.04490, 3.10142, 3.15782, 3.21411,
    3.27028, 3.32634, 3.38230, 3.43813, 3.49387,
    2.40380, 2.46138, 2.51885, 2.57623, 2.63351, 2.69067, 2.74773, 2.80469,
    2.86152, 2.91826, 2.97489, 3.03140, 3.08781, 3.14410, 3.20031, 3.25638,
    3.31237, 3.36825, 3.42402, 3.47970, 3.53527,
    2.45045, 2.50773, 2.56492, 2.62201, 2.67900, 2.73590, 2.79270, 2.84940,
    2.90599, 2.96250, 3.01890, 3.07518, 3.13138, 3.18748, 3.24347, 3.29937,
    3.35515, 3.41085, 3.46646, 3.52196, 3.57738,
    2.49750, 2.55450, 2.61143, 2.66825, 2.72499, 2.78163, 2.83818, 2.89464,
    2.95101, 3.00729, 3.06346, 3.11954, 3.17554, 3.23143, 3.28723, 3.34294,
    3.39856, 3.45409, 3.50952, 3.56488, 3.62014,
    2.54489, 2.60164, 2.65830, 2.71488, 2.77134, 2.82780, 2.88412, 2.94034,
    2.99650, 3.05256, 3.10853, 3.16442, 3.22021, 3.27592, 3.33154, 3.38707,
    3.44253, 3.49789, 3.55316, 3.60836, 3.66348,
    2.59259, 2.64910, 2.70553, 2.76188, 2.81815, 2.87434, 2.93044, 2.98646,
    3.04241, 3.09827, 3.15404, 3.20974, 3.26536, 3.32089, 3.37633, 3.43170,
    3.48699, 3.54219, 3.59731, 3.65237, 3.70734,
    2.64054, 2.69684, 2.75305, 2.80918, 2.86523, 2.92122, 2.97712, 3.03295,
    3.08870, 3.14437, 3.19996, 3.25548, 3.31091, 3.36627, 3.42156, 3.47677,
    3.53190, 3.58695, 3.64193, 3.69684, 3.75167,
    2.68873, 2.74482, 2.80083, 2.85676, 2.91262, 2.96841, 3.02412, 3.07976,
    3.13533, 3.19082, 3.24623, 3.30157, 3.35685, 3.41205, 3.46718, 3.52223,
    3.57721, 3.63213, 3.68696, 3.74174, 3.79644,
    2.73713, 2.79302, 2.84884, 2.90459, 2.96027, 3.01587, 3.07140, 3.12686,
    3.18225, 3.23757, 3.29282, 3.34801, 3.40313, 3.45815, 3.51315, 3.56805,
    3.62290, 3.67767, 3.73237, 3.78701, 3.84159,
    2.78572, 2.84143, 2.89707, 2.95264, 3.00813, 3.06356, 3.11892, 3.17422,
    3.22946, 3.28462, 3.33971, 3.39474, 3.44971, 3.50460, 3.55944, 3.61421,
    3.66891, 3.72356, 3.77814, 3.83264, 3.88710,
    2.83446, 2.89000, 2.94547, 3.00088, 3.05621, 3.11147, 3.16669, 3.22183,
    3.27689, 3.33191, 3.38685, 3.44174, 3.49656, 3.55132, 3.60602, 3.66066,
    3.71523, 3.76975, 3.82421, 3.87860, 3.93293,
    2.88335, 2.93873, 2.99404, 3.04929, 3.10447, 3.15959, 3.21464, 3.26963,
    3.32456, 3.37943, 3.43424, 3.48898, 3.54366, 3.59830, 3.65287, 3.70737,
    3.76183, 3.81622, 3.87056, 3.92484, 3.97905,
    2.93234, 2.98760, 3.04277, 3.09786, 3.15291, 3.20787, 3.26278, 3.31764,
    3.37242, 3.42716, 3.48184, 3.53662, 3.59100, 3.64551, 3.69995, 3.75435,
    3.80867, 3.86295, 3.91718, 3.97134, 4.02545,
    2.98151, 3.03660, 3.09163, 3.14659, 3.20149, 3.25632, 3.31110, 3.36581,
    3.42047, 3.47507, 3.52963, 3.58411, 3.63855, 3.69293, 3.74725, 3.80153,
    3.85575, 3.90991, 3.96403, 4.01809, 4.07211,
    3.03074, 3.08571, 3.14061, 3.19543, 3.25021, 3.30491, 3.35956, 3.41415,
    3.46869, 3.52317, 3.57759, 3.63196, 3.68628, 3.74054, 3.79476, 3.84893,
    3.90303, 3.95709, 4.01111, 4.06506, 4.11897,
    3.08008, 3.13492, 3.18970, 3.24440, 3.29905, 3.35363, 3.40815, 3.46263,
    3.51704, 3.57141, 3.62572, 3.67998, 3.73418, 3.78834, 3.84244, 3.89650,
    3.95051, 4.00447, 4.05837, 4.11224, 4.16605,
    3.12949, 3.18422, 3.23888, 3.29347, 3.34800, 3.40247, 3.45688, 3.51124,
    3.56554, 3.61980, 3.67400, 3.72815, 3.78225, 3.83629, 3.89030, 3.94425,
    3.99816, 4.05203, 4.10583, 4.15961, 4.21333,
    3.17899, 3.23361, 3.28815, 3.34264, 3.39706, 3.45142, 3.50571, 3.55997,
    3.61416, 3.66831, 3.72241, 3.77645, 3.83046, 3.88440, 3.93831, 3.99216,
    4.04598, 4.09974, 4.15347, 4.20715, 4.26078,
    3.22855, 3.28307, 3.33751, 3.39188, 3.44620, 3.50046, 3.55466, 3.60880,
    3.66290, 3.71694, 3.77095, 3.82489, 3.87880, 3.93265, 3.98646, 4.04022,
    4.09395, 4.14762, 4.20126, 4.25485, 4.30840,
]
_GB_D0 = [d / 10.0 for d in _GB_D0]

_GB_M0 = [
    0.0381511, 0.0338587, 0.0301776, 0.0270030, 0.0242506, 0.0218529,
    0.0197547, 0.0179109, 0.0162844, 0.0148442, 0.0135647, 0.0124243,
    0.0114047, 0.0104906, 0.00966876, 0.008928, 0.0082587, 0.00765255,
    0.00710237, 0.00660196, 0.00614589,
    0.0396198, 0.0351837, 0.0313767, 0.0280911, 0.0252409, 0.0227563,
    0.0205808, 0.0186681, 0.0169799, 0.0154843, 0.014155, 0.0129696,
    0.0119094, 0.0109584, 0.0101031, 0.00933189, 0.0086348, 0.00800326,
    0.00742986, 0.00690814, 0.00643255,
    0.041048, 0.0364738, 0.0325456, 0.0291532, 0.0262084, 0.0236399,
    0.0213897, 0.0194102, 0.0176622, 0.0161129, 0.0147351, 0.0135059,
    0.0124061, 0.0114192, 0.0105312, 0.00973027, 0.00900602, 0.00834965,
    0.0077535, 0.00721091, 0.00671609,
    0.0424365, 0.0377295, 0.0336846, 0.0301893, 0.0271533, 0.0245038,
    0.0221813, 0.0201371, 0.018331, 0.0167295, 0.0153047, 0.014033,
    0.0128946, 0.0118727, 0.0109529, 0.0101229, 0.00937212, 0.00869147,
    0.00807306, 0.00751003, 0.00699641,
    0.0437861, 0.0389516, 0.0347944, 0.0311998, 0.0280758, 0.0253479,
    0.0229555, 0.0208487, 0.0189864, 0.0173343, 0.0158637, 0.0145507,
    0.0133748, 0.0123188, 0.0113679, 0.0105096, 0.0097329, 0.00902853,
    0.00838835, 0.00780533, 0.0072733,
    0.0450979, 0.0401406, 0.0358753, 0.0321851, 0.0289761, 0.0261726,
    0.0237125, 0.0215451, 0.0196282, 0.017927, 0.0164121, 0.0150588,
    0.0138465, 0.0127573, 0.0117761, 0.0108902, 0.0100882, 0.00936068,
    0.00869923, 0.00809665, 0.00754661,
    0.0463729, 0.0412976, 0.0369281, 0.0331456, 0.0298547, 0.026978,
    0.0244525, 0.0222264, 0.0202567, 0.0185078, 0.0169498, 0.0155575,
    0.0143096, 0.0131881, 0.0121775, 0.0112646, 0.010438, 0.00968781,
    0.00900559, 0.00838388, 0.00781622,
    0.0476123, 0.0424233, 0.0379534, 0.034082, 0.0307118, 0.0277645,
    0.0251757, 0.0228927, 0.0208718, 0.0190767, 0.0174768, 0.0160466,
    0.0147642, 0.0136112, 0.0125719, 0.0116328, 0.0107821, 0.0100099,
    0.00930735, 0.00866695, 0.00808206,
    0.0488171, 0.0435186, 0.038952, 0.0349947, 0.0315481, 0.0285324,
    0.0258824, 0.0235443, 0.0214738, 0.0196339, 0.0179934, 0.0165262,
    0.0152103, 0.0140267, 0.0129595, 0.0119947, 0.0111206, 0.0103268,
    0.00960445, 0.00894579, 0.00834405,
    0.0499883, 0.0445845, 0.0399246, 0.0358844, 0.032364, 0.0292822,
    0.0265729, 0.0241815, 0.0220629, 0.0201794, 0.0184994, 0.0169964,
    0.0156479, 0.0144345, 0.0133401, 0.0123504, 0.0114534, 0.0106386,
    0.00989687, 0.00922037, 0.00860216,
    0.0511272, 0.0456219, 0.040872, 0.0367518, 0.0331599, 0.0300142,
    0.0272475, 0.0248045, 0.0226392, 0.0207135, 0.0189952, 0.0174574,
    0.0160771, 0.0148348, 0.0137138, 0.0126998, 0.0117805, 0.0109452,
    0.0101846, 0.00949067, 0.00885636,
    0.0522348, 0.0466315, 0.0417948, 0.0375973, 0.0339365, 0.030729,
    0.0279067, 0.0254136, 0.023203, 0.0212363, 0.0194809, 0.0179092,
    0.016498, 0.0152275, 0.0140807, 0.013043, 0.012102, 0.0112466,
    0.0104676, 0.00975668, 0.00910664,
    0.0533123, 0.0476145, 0.042694, 0.0384218, 0.0346942, 0.0314268,
    0.0285507, 0.026009, 0.0237547, 0.0217482, 0.0199566, 0.018352,
    0.0169108, 0.0156128, 0.0144408, 0.0133801, 0.0124179, 0.011543,
    0.010746, 0.0100184, 0.00935302,
    0.0543606, 0.0485716, 0.04357, 0.0392257, 0.0354335, 0.0321082,
    0.02918, 0.0265913, 0.0242943, 0.0222492, 0.0204225, 0.0187859,
    0.0173155, 0.0159908, 0.0147943, 0.0137111, 0.0127282, 0.0118343,
    0.0110197, 0.0102759, 0.00959549,
    0.0553807, 0.0495037, 0.0444239, 0.0400097, 0.0361551, 0.0327736,
    0.0297949, 0.0271605, 0.0248222, 0.0227396, 0.0208788, 0.0192111,
    0.0177122, 0.0163615, 0.0151413, 0.0140361, 0.013033, 0.0121206,
    0.0112888, 0.0105292, 0.00983409,
    0.0563738, 0.0504116, 0.0452562, 0.0407745, 0.0368593, 0.0334235,
    0.0303958, 0.0277171, 0.0253387, 0.0232197, 0.0213257, 0.0196277,
    0.0181013, 0.0167252, 0.0154817, 0.0143552, 0.0133325, 0.0124019,
    0.0115534, 0.0107783, 0.0100688,
    0.0573406, 0.0512963, 0.0460676, 0.0415206, 0.0375468, 0.0340583,
    0.030983, 0.0282614, 0.0258441, 0.0236896, 0.0217634, 0.020036,
    0.0184826, 0.017082, 0.0158158, 0.0146685, 0.0136266, 0.0126783,
    0.0118135, 0.0110232, 0.0102998,
    0.0582822, 0.0521584, 0.0468589, 0.0422486, 0.038218, 0.0346784,
    0.0315571, 0.0287938, 0.0263386, 0.0241497, 0.0221922, 0.0204362,
    0.0188566, 0.0174319, 0.0161437, 0.0149761, 0.0139154, 0.0129499,
    0.0120691, 0.0112641, 0.0105269,
    0.0591994, 0.0529987, 0.0476307, 0.042959, 0.0388734, 0.0352843,
    0.0321182, 0.0293144, 0.0268225, 0.0246002, 0.0226121, 0.0208283,
    0.0192232, 0.0177751, 0.0164654, 0.015278, 0.0141991, 0.0132167,
    0.0123204, 0.0115009, 0.0107504,
    0.0600932, 0.053818, 0.0483836, 0.0436525, 0.0395136, 0.0358764,
    0.0326669, 0.0298237, 0.0272961, 0.0250413, 0.0230236, 0.0212126,
    0.0195826, 0.0181118, 0.0167811, 0.0155744, 0.0144778, 0.0134789,
    0.0125673, 0.0117338, 0.0109702,
    0.0609642, 0.0546169, 0.0491183, 0.0443295, 0.0401388, 0.036455,
    0.0332033, 0.030322, 0.0277596, 0.0254732, 0.0234266, 0.0215892,
    0.0199351, 0.018442, 0.0170909, 0.0158654, 0.0147514, 0.0137365,
    0.0128101, 0.0119627, 0.0111863,
]
_GB_M0 = [m * 10.0 for m in _GB_M0]

# Screen (scaling) parameters: element -> (normal, GBn, GBn2, GBn2_nucleic)
# Populated lazily after openmm.app is available; resolved in _add_implicit_solvent.
_GB_SCREEN_BY_ELEMENT_SYMBOL = {
    'H':  (0.85,   1.09085413633,  1.425952,  1.696538),
    'D':  (0.85,   1.09085413633,  1.425952,  1.696538),
    'C':  (0.72,   0.48435382330,  1.058554,  1.268902),
    'N':  (0.79,   0.700147318409, 0.733599,  1.4259728),
    'O':  (0.85,   1.06557401132,  1.061039,  0.1840098),
    'F':  (0.88,   0.5,            0.5,       0.5),
    'P':  (0.86,   0.5,            0.5,       1.5450597),
    'S':  (0.96,   0.602256336067, -0.703469, 0.05),
}
_GB_SCREEN_DEFAULT = (0.8, 0.5, 0.5, 0.5)

_GB_NUCLEIC_RESIDUES = frozenset(['A', 'C', 'G', 'U', 'DA', 'DC', 'DG', 'DT'])

# GBn2 per-element alpha/beta/gamma (protein and nucleic variants)
_GBN2_ATOM_PARAMS = {
    'H':  [0.788440, 0.798699, 0.437334],
    'D':  [0.788440, 0.798699, 0.437334],
    'C':  [0.733756, 0.506378, 0.205844],
    'N':  [0.503364, 0.316828, 0.192915],
    'O':  [0.867814, 0.876635, 0.387882],
    'S':  [0.867814, 0.876635, 0.387882],
}
_GBN2_ATOM_PARAMS_NUCLEIC = {
    'H':  [0.537050, 0.362861, 0.116704],
    'D':  [0.537050, 0.362861, 0.116704],
    'C':  [0.331670, 0.196842, 0.093422],
    'N':  [0.686311, 0.463189, 0.138722],
    'O':  [0.606344, 0.463006, 0.142262],
    'S':  [0.606344, 0.463006, 0.142262],
    'P':  [0.418365, 0.290054, 0.1064245],
}
_GBN2_DEFAULT_PARAMS = [1.0, 0.8, 4.851]


class EvbSystemBuilder():

    def __init__(self, comm=None, ostream=None):
        '''
        Initialize the EVB driver class.
        '''

        if comm is None:
            comm = MPI.COMM_WORLD

        if ostream is None:
            if comm.Get_rank() == mpi_master():
                ostream = OutputStream(sys.stdout)
            else:
                ostream = OutputStream(None)

        # output stream
        self.ostream = ostream

        # MPI information
        self.comm = comm
        self.rank = comm.Get_rank()
        self.size = comm.Get_size()

        self.temperature: float = 300
        self.Lambda: list[float]

        self.sc_alpha_lj: float = 0.7
        self.sc_alpha_q: float = 0.3
        self.sc_sigma_q: float = 1.0
        self.sc_power: float = 1 / 6  # The exponential power in the soft core expression
        self.morse_D_default: float = 500  # kj/mol, default dissociation energy if none is given
        self.morse_couple: float = 1  # kj/mol, scaling for the morse potential to emulate a coupling between two overlapping bonded states

        # self.centroid_k: float = 50000  # kj/mol nm, force constant for the position restraints
        # self.centroid_offset: float = -0.2 #A
        # self.centroid_complete_ligand: bool = True # If false, the centroid force will only be applied to the reacting atoms, if true, the complete ligand will be used for the centroid force

        self.posres_residue_radius: float = -1  # A, cutoff measured from the com of the ligand for which residues to add to the position restraints, -1 includes the whole molecule
        self.posres_k = 1000  # kj/mol nm^2, default in gromacs

        # self.restraint_r_default: float = 0.5  # nm, default position restraint distance if none is given
        # self.restraint_r_offset: float = 0.1  # nm, distance added to the measured distance in a structure to set the position restraint distance
        self.coul14_scale: float = 0.833
        self.lj14_scale: float = 0.5
        self.nb_cutoff: float = 1.  # nm, minimal cutoff for the nonbonded force

        self.pressure: float = -1.
        self.solvent: str = None  # type: ignore
        self.padding: float = 1.5
        self.no_reactant: bool = False
        self.E_field: list[float] = [0, 0, 0]
        self.neutralize: bool = False

        self.soft_core_coulomb_pes_static = False
        self.soft_core_lj_pes_static = False

        self.soft_core_coulomb_pes_dynamic = True
        self.soft_core_lj_pes_dynamic = True

        self.soft_core_coulomb_int = False
        self.soft_core_lj_int = False

        # self.bonded_integration: bool = True  # If the integration potential should use bonded (harmonic/morse) forces for forming/breaking bonds, instead of replacing them with nonbonded potentials
        self.bonded_integration_bond_fac: float = 0.2  # Scaling factor for the bonded integration forces.
        self.bonded_integration_angle_fac: float = 0.1  # Scaling factor for the bonded integration forces.
        self.torsion_lambda_switch: float = 0.25  # The minimum (1-maximum) lambda value at which to start turning on (have turned of) the proper torsion for the product (reactant)

        self.int_nb_const_exceptions = True  # If the exceptions for the integration nonbonded force should be kept constant over the entire simulation

        self.dynamic_bond_tightening = 3.25  # Power for tigthening breaking and forming bonds during the integration
        self.dynamic_bond_fc_factor = 0.8  # Force constant factor for the dynamic bond tightening

        self.verbose = False

        self.constraints: list[dict] = []

        self.k = 4.184 * hartree_in_kcalpermol() * 0.1 * bohr_in_angstrom(
        )  # Coulombic pre-factor

        self.deg_to_rad: float = np.pi / 180

        self.data_folder: str | None = None
        self.run_folder: str | None = None
        self.pdb: str | None = None
        self.pdb_active_res: list[
            dict] | None = None  # Residue ids of the residues active in the reaction in the PDB file

        self.no_force_groups: bool = False
        self.nb_switching_function: bool = True

        self.graphene = False
        self.graphene_size_nm = 4
        self.CNT = False
        self.CNT_radius_nm = 0.5

        self.begin_index = 0
        self.water_model: str
        self.decompose_bonded = True
        self.decompose_nb: list | None = None

        self.implicit_solvent_model: str | None = None
        self.solute_dielectric: float = 1.0
        self.solvent_dielectric: float = 78.39

        # Conformational TS parameters.  Set by build_systems when
        # 'conformer_active_torsion' is present in the configuration dict.
        # conformer_active_torsion is handled outside the keywords loop
        # because it is a tuple rather than a plain scalar.
        self.conformer_active_torsion: tuple | None = None
        self.conformer_phi_reactant: float = 0.0
        self.conformer_phi_product: float = 0.0
        self.conformer_k: float = 200.0

        self.keywords = {
            "temperature": float,  # -> system dependent
            "nb_cutoff": float,  # ->
            # "bonded_integration": bool,
            "bonded_integration_bond_fac": float,
            "bonded_integration_angle_fac": float,
            "torsion_lambda_switch": float,
            "soft_core_coulomb_pes_static": bool,
            "soft_core_lj_pes_static": bool,
            "soft_core_coulomb_pes_dynamic": bool,
            "soft_core_lj_pes_dynamic": bool,
            "soft_core_coulomb_int": bool,
            "soft_core_lj_int": bool,
            "int_nb_const_exceptions": bool,
            "dynamic_bond_tightening": float,
            "dynamic_bond_fc_factor": float,
            "pressure": float,
            "solvent": str,
            "padding": float,
            "no_reactant": bool,
            "E_field": list,
            "neutralize": bool,
            "morse_D_default": float,
            "morse_couple": float,
            "posres_k": float,
            "posres_residue_radius": float,
            "coul14_scale": float,
            "lj14_scale": float,
            "pdb": str,
            "pdb_active_res": list,
            "no_force_groups": bool,
            "nb_switching_function": bool,
            'graphene': bool,
            'graphene_size_nm': float,
            "CNT": bool,
            "CNT_radius_nm": float,
            "decompose_nb": list,
            "decompose_bonded": bool,
            "implicit_solvent_model": str,
            "solute_dielectric": float,
            "solvent_dielectric": float,
            "conformer_phi_reactant": float,
            "conformer_phi_product": float,
            "conformer_k": float,
        }

    def build_systems(
        self,
        reactant: MMForceFieldGenerator,
        product: MMForceFieldGenerator,
        Lambda: list,
        configuration: dict,
        constraints: list | None = None,
    ):

        assert_msg_critical('openmm' in sys.modules,
                            'openmm is required for EvbSystemBuilder.')

        for keyword, val in self.keywords.items():
            if keyword in configuration:
                if (not isinstance(configuration[keyword], val)
                        and not (isinstance(configuration[keyword], int)
                                 and val == float)):
                    raise TypeError(
                        f"Configuration option {keyword} should be of type {val}, but got {type(configuration[keyword])} instead."
                    )
                else:
                    setattr(self, keyword, configuration[keyword])
                    self.ostream.print_info(
                        f"{keyword}: {getattr(self, keyword)}")

            else:
                self.ostream.print_info(
                    f"{keyword}: {getattr(self, keyword)} (default)")
        # conformer_active_torsion is a tuple, handled separately
        if 'conformer_active_torsion' in configuration:
            self.conformer_active_torsion = tuple(
                configuration['conformer_active_torsion'])
            self.ostream.print_info(
                f"conformer_active_torsion: {self.conformer_active_torsion}")
        self.ostream.flush()
        self.reactant = reactant
        self.product = product

        if constraints is not None:
            self.constraints = constraints

        if self.pdb is None:
            system = mm.System()
            topology = mmapp.Topology()

            nb_force = mm.NonbondedForce()

            cmm_remover = mm.CMMotionRemover()
            system.addForce(nb_force)
            system.addForce(cmm_remover)
            vlx_mol = Molecule(reactant.molecule)
            self.positions = vlx_mol.get_coordinates_in_angstrom()
            self.reaction_atoms = dict()
        else:
            system, topology, vlx_mol, self.reaction_atoms = self._system_from_pdb(
            )
            nb_force = [
                force for force in system.getForces()
                if isinstance(force, mm.NonbondedForce)
            ][0]
            cmm_remover = [
                force for force in system.getForces()
                if isinstance(force, mm.CMMotionRemover)
            ][0]

        nb_force.setName("NonbondedForce")
        if self.solvent:
            nb_force.setNonbondedMethod(mm.NonbondedForce.PME)
        else:
            nb_force.setNonbondedMethod(mm.NonbondedForce.NoCutoff)
        cmm_remover.setName("CMMotionRemover")
        if not self.no_force_groups:
            nb_force.setForceGroup(EvbForceGroup.NB_FORCE_INT.value)
            cmm_remover.setForceGroup(EvbForceGroup.CMM_REMOVER.value)

        self._add_reactant(system, topology, nb_force)

        # Set the positions and make a box for it

        box = None
        if self.CNT or self.graphene:
            box = self._add_CNT_graphene(system, nb_force, topology, vlx_mol)

        # assert False, "Rethink CNT/graphene input"

        if self.solvent:
            box = self._configure_pbc(system, topology, nb_force, box)  # A
            # todo what about the box, and especially giving it to the solvator
            box = self._add_solvent(system, vlx_mol, self.solvent, topology,
                                    nb_force, self.neutralize, self.padding,
                                    box)

        if self.pressure > 0:
            barostat_not_used = self._add_barostat(system)

        if self.implicit_solvent_model is not None:
            assert_msg_critical(
                not self.solvent,
                'EvbSystemBuilder: implicit_solvent_model and solvent (explicit) '
                'cannot be used simultaneously.')
            assert_msg_critical(
                not self.CNT and not self.graphene,
                'EvbSystemBuilder: implicit_solvent_model is not compatible '
                'with CNT or graphene environments.')
            self._add_implicit_solvent(system, topology, nb_force)

        E_field = None
        if np.any(np.array(self.E_field) > 0.001):
            E_field = self._add_E_field(system, self.E_field)

        self.topology: mmapp.Topology = topology
        self.systems = self._interpolate_system(
            system,
            Lambda,
            nb_force,
            E_field,
        )
        # Add all lambda dependent parameters

        self.ostream.flush()
        return self.systems, self.topology, self.positions

    def _system_from_pdb(self):
        pdb_file = mmapp.PDBFile(self.pdb)
        topology = pdb_file.getTopology()
        system_mol = Molecule.read_pdb_file(self.pdb)
        posres_atoms = [atom for atom in topology.atoms()]

        forcefield = mmapp.ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
        templates, residues = forcefield.generateTemplatesForUnmatchedResidues(
            topology)

        for t, template in enumerate(templates):
            for i, atom in enumerate(template.atoms):
                forcefield.registerAtomType({
                    'name':
                    f'evb_{atom.name}_{t}',
                    'class':
                    'evb_placeholder',
                    'mass':
                    atom.element.mass.value_in_unit(mmunit.dalton),
                    'element':
                    atom.element
                })
                template.atoms[i].type = f'evb_{atom.name}_{t}'
            forcefield.registerResidueTemplate(template)

        nbforce = forcefield.getGenerators()[2]
        nbforce.registerAtom({
            'class': 'evb_placeholder',
            'charge': 0.0,
            'epsilon': 0.0,
            'sigma': 1.0,
        })

        system = forcefield.createSystem(
            topology,
            removeCMMotion=True,
        )

        self.positions = system_mol.get_coordinates_in_angstrom()
        self._add_posres(system, posres_atoms, self.positions)

        reaction_atoms = {}
        if self.no_reactant:
            return system, topology, system_mol, reaction_atoms

        for res_dict in self.pdb_active_res:
            # all atoms already exist in the topology, pdb: top, matching chain id with removed_chain
            # the pdb residue should be checked against the reactant, and an index mapping should be found out #todo
            # these should be added to the reaction_atoms
            chain = [c for c in topology.chains()
                     if c.id == res_dict['chain']][0]
            residue = [
                r for r in chain.residues() if int(r.id) == res_dict['residue']
            ][0]

            mapping = self._get_mapped_atom_ids_from_residue(residue)
            res_atoms = [atom for atom in residue.atoms()]

            for id in mapping.keys():
                self.reactant.atoms[id]['pdb'] = 'sys'
            for vlx_id, res_id in mapping.items():
                reaction_atoms.update({
                    vlx_id:
                    [atom for atom in res_atoms if atom.index == res_id][0]
                })

            system = self._delete_pdb_forces(system,
                                             [atom.index for atom in res_atoms])

            self.ostream.print_info(
                f"Reacting residue {residue.name} {residue.id} with {len(res_atoms)} atoms added to reaction_atoms"
            )
            self.ostream.flush()

        return system, topology, system_mol, reaction_atoms

    def _get_mapped_atom_ids_from_residue(self, residue):
        # create a graph of the residue
        # if the bonds are not available, create them based on proximity from the positions

        # figure out a mapping from self.reactant.atoms to the atoms in the residue with networkx
        vlx_elements = self.reactant.molecule.get_element_ids()
        vlx_ids = list(self.reactant.atoms.keys())
        vlx_bonds = list(self.reactant.bonds.keys())
        residue_elements = [
            atom.element.atomic_number for atom in residue.atoms()
        ]
        residue_ids = [atom.index for atom in residue.atoms()]

        # Depending on how the residue is stored in the PDB, the bonds might not be available.
        # In that case, the bonds are derived from the connectivity-matrix of a newly created molecule
        residue_bonds = [bond for bond in residue.bonds()]
        if len(residue_bonds) == 0:
            residue_bonds = []
            mol = Molecule()
            for id, element in zip(residue_ids, residue_elements):
                coord = self.positions[id]
                mol.add_atom(int(element), Point(coord), 'angstrom')
            connectivity_matrix = mol.get_connectivity_matrix()
            for i, id in enumerate(residue_ids):
                for j, jd in enumerate(residue_ids):
                    if i >= j:
                        continue
                    if connectivity_matrix[i, j] == 1:
                        residue_bonds.append((id, jd))
        else:
            residue_bonds = [(bond.atom1.index, bond.atom2.index)
                             for bond in residue.bonds()
                             if bond.atom1.index in residue_ids
                             and bond.atom2.index in residue_ids]

        vlx_graph = nx.Graph()
        vlx_graph.add_nodes_from(vlx_ids)
        vlx_graph.add_edges_from(vlx_bonds)
        for i, elem in zip(vlx_ids, vlx_elements):
            vlx_graph.nodes[i]['elem'] = elem

        res_graph = nx.Graph()
        res_graph.add_nodes_from(residue_ids)
        res_graph.add_edges_from(residue_bonds)
        for i, elem in zip(residue_ids, residue_elements):
            res_graph.nodes[i]['elem'] = elem

        GM = GraphMatcher(vlx_graph, res_graph,
                          categorical_node_match('elem', ''))
        if not GM.subgraph_is_isomorphic():
            raise ValueError(
                f"Could not find subgraph isomorphism between the residue {residue.name} {residue.index} and the reactant molecule"
            )
        mapping = next(GM.subgraph_isomorphisms_iter())
        return mapping

    def _delete_pdb_forces(self, system, del_indices):
        # set the right force groups, give descriptive names to the pdb forces, and delete all contributions that solely have the given indices
        pdb_forces = system.getForces()
        bond_force = [
            force for force in pdb_forces
            if isinstance(force, mm.HarmonicBondForce)
        ][0]
        angle_force = [
            force for force in pdb_forces
            if isinstance(force, mm.HarmonicAngleForce)
        ][0]
        torsion_force = [
            force for force in pdb_forces
            if isinstance(force, mm.PeriodicTorsionForce)
        ][0]
        bond_force.setForceGroup(EvbForceGroup.PDB.value)
        angle_force.setForceGroup(EvbForceGroup.PDB.value)
        torsion_force.setForceGroup(EvbForceGroup.PDB.value)
        bond_force.setName("PDB Bond Force")
        angle_force.setName("PDB Angle Force")
        torsion_force.setName("PDB Torsion Force")

        for i in range(bond_force.getNumBonds()):
            params = bond_force.getBondParameters(i)
            if params[0] in del_indices and params[1] in del_indices:
                bond_force.setBondParameters(i, params[0], params[1], 1, 0)

        for i in range(angle_force.getNumAngles()):
            params = angle_force.getAngleParameters(i)
            if params[0] in del_indices and params[1] in del_indices and params[
                    2] in del_indices:
                angle_force.setAngleParameters(i, params[0], params[1],
                                               params[2], 1, 0)

        for i in range(torsion_force.getNumTorsions()):
            params = torsion_force.getTorsionParameters(i)
            if params[0] in del_indices and params[1] in del_indices and params[
                    2] in del_indices and params[3] in del_indices:
                torsion_force.setTorsionParameters(i, params[0], params[1],
                                                   params[2], params[3], 1, 0,
                                                   0)
        return system

    def _configure_pbc(self, system, topology, nb_force, box=None):
        if box is None:
            box = [-1., -1., -1.]
        minim = [
            2 * self.padding + 0.1 *
            (max(self.positions[:, 0]) - min(self.positions[:, 0])),
            2 * self.padding + 0.1 *
            (max(self.positions[:, 1]) - min(self.positions[:, 1])),
            2 * self.padding + 0.1 *
            (max(self.positions[:, 2]) - min(self.positions[:, 2]))
        ]
        dims = ['x', 'y', 'z']
        for i in range(3):
            if box[i] == -1:
                box[i] = minim[i]
            else:
                self.ostream.print_info(
                    f"Box size calculation for {dims[i]}-component is being overridden by graphene or CNT with value of {box[i]}"
                )
                # if box[i]<minim[i]:
                #     self.ostream.print_warning(f"Provided {dims[i]}-component of box size {box[i]} is smaller then {minim[i]}. Consider increasing the box size (or the size of the graphene or CNT)")
        self.ostream.print_info(
            f"Size of the system: {minim[0]:.3f} x {minim[1]:.3f} x {minim[2]:.3f}, padding: {self.padding:.3f} nm."
        )
        self.ostream.print_info(
            f"Building system in box with dimensions {box[0]:.3f} x {box[1]:.3f} x {box[2]:.3f} nm"
        )
        vector_box: list[mm.Vec3] = [
            mm.Vec3(box[0], 0, 0),
            mm.Vec3(0, box[1], 0),
            mm.Vec3(0, 0, box[2])
        ]
        list_box = [[box[0], 0, 0], [0, box[1], 0], [0, 0, box[2]]]
        system.setDefaultPeriodicBoxVectors(*vector_box)
        topology.setPeriodicBoxVectors(list_box)

        nb_force.setCutoffDistance(self.nb_cutoff)
        self.ostream.print_info(
            f"Setting nonbonded cutoff to {self.nb_cutoff:.3f} nm")
        if self.nb_switching_function:
            nb_force.setUseSwitchingFunction(True)
            nb_force.setSwitchingDistance(0.9 * self.nb_cutoff)
        else:
            nb_force.setUseSwitchingFunction(False)
            nb_force.setSwitchingDistance(-1)
        box[0] *= 10
        box[1] *= 10
        box[2] *= 10
        self.ostream.flush()
        return box

    def _add_reactant(self, system, topology, nb_force):
        reaction_chain = topology.addChain()
        reaction_residue = topology.addResidue(name="REA", chain=reaction_chain)
        if not self.no_reactant:
            # add atoms of the solute to the topology and the system
            elements = self.reactant.molecule.get_labels()

            for id, atom in self.reactant.atoms.items():
                # if the pdb field is sys, the atom is already in the system and added to the nbforce
                if atom.get('pdb') == 'sys':
                    continue
                mm_element = mmapp.Element.getBySymbol(elements[id])
                name = f"{elements[id]}{id}"
                system.addParticle(mm_element.mass)
                nb_force.addParticle(
                    0, 1, 0
                )  # Placeholder values, actual values depend on lambda and will be set later

                # If a pdb field is defined, the atom is already in the topoolgy and hence also in the reaction_atoms
                if atom.get('pdb') is not None:
                    continue
                reaction_atom = topology.addAtom(
                    name,
                    mm_element,
                    reaction_residue,
                )
                self.reaction_atoms[id] = reaction_atom

            # Make sure the solute does not interact with itself through the default nonbonded force, as there will be another nonbonded force to take care of this
            # exception_params = [{nb_force.getExceptionParameters(i)[0],nb_force.getExceptionParameters(i)[1]} for i in range(nb_force.getNumExceptions())]
            atom_indices = [atom.index for atom in self.reaction_atoms.values()]
            set_exceptions = []
            for i in range(nb_force.getNumExceptions()):
                [atom_i, atom_j, charge, sigma,
                 epsilon] = nb_force.getExceptionParameters(i)
                if atom_i in atom_indices and atom_j in atom_indices:
                    nb_force.setExceptionParameters(i, atom_i, atom_j, 0.0, 1.0,
                                                    0.0)
                    set_exceptions.append((atom_i, atom_j))

            for atom_i in atom_indices:
                for atom_j in atom_indices:
                    # Skip any capping hydrogens, because they are not in the topology
                    if atom_i >= atom_j:
                        continue

                    if (atom_i, atom_j) in set_exceptions or (
                            atom_j, atom_i) in set_exceptions:

                        continue
                    nb_force.addException(
                        atom_i,
                        atom_j,
                        0.0,
                        1.0,
                        0.0,
                    )

            self.ostream.flush()

            for bond in self.reactant.bonds.keys():
                if not (self.reactant.atoms[bond[0]].get('pdb') is None
                        and self.reactant.atoms[bond[1]].get('pdb') is None):
                    continue
                topology.addBond(
                    self.reaction_atoms[bond[0]],
                    self.reaction_atoms[bond[1]],
                )

        return

    def _add_posres(self, system, atoms, positions):
        posres_expr = "posres_k*periodicdistance(x, y, z, x0, y0, z0)^2"
        posres_force = mm.CustomExternalForce(posres_expr)
        posres_force.setName("protein_ligand_posres")
        posres_force.setForceGroup(EvbForceGroup.POSRES.value)
        posres_force.addGlobalParameter('posres_k', self.posres_k)
        posres_force.addPerParticleParameter('x0')
        posres_force.addPerParticleParameter('y0')
        posres_force.addPerParticleParameter('z0')
        count = 0
        for atom in atoms:
            if atom.element is not mmapp.element.hydrogen:
                index = atom.index
                position = self.positions[index]
                posres_force.addParticle(index, position * 0.1)
                count += 1

        self.ostream.print_info(f"Adding {count} particles to posres force")
        self.ostream.flush()
        system.addForce(posres_force)

    def _add_CNT_graphene(self, system, nb_force, topology, system_mol):

        assert_msg_critical('openmm' in sys.modules,
                            'openmm is required for EvbSystemBuilder.')

        if self.CNT and self.graphene:
            raise ValueError(
                "CNT and Graphene cannot be used simultaneously, pick one please"
            )

        # cc_eq = 0.14080860183951827  # nm
        cc_eq = 0.1418  # nm, distance between two carbon atoms in graphene
        x_disp = 3 * cc_eq * 10  # A, size of the graphene unit cell in the x direction
        y_disp = math.sqrt(
            3
        ) * cc_eq * 10  # A, size of the graphene unit cell in the y direction
        box = [-1., -1., -1.]
        if self.graphene:
            x_minim = self.graphene_size_nm * 10  # A
            y_minim = self.graphene_size_nm * 10  # A
            M = math.ceil(x_minim / x_disp)
            N = math.ceil(y_minim / y_disp)
            X = M * x_disp  # A
            Y = N * y_disp  # A
            self.ostream.print_info(
                f"Box size x: {box[0]:.3f} y: {box[1]:.3f} nm, minimum graphene size: {self.graphene_size_nm:.3f} nm. Using largest to determine the size of the graphene sheet and size of the box"
            )
            self.ostream.print_info(
                f"Building graphene sheet with X: {X/10:.3f} nm, Y: {Y/10:.3f} nm"
            )
            box[0] = X / 10  # nm
            box[1] = Y / 10  # nm

        else:
            x_minim = 2 * self.padding * 10 + (
                max(self.positions[:, 0]) - min(self.positions[:, 0]))  # A
            y_minim = self.CNT_radius_nm * 2 * np.pi * 10  # A
            M = math.ceil(x_minim / x_disp)
            N = math.ceil(y_minim / y_disp)
            X = M * x_disp  # A
            Y = N * y_disp  # A
            R = Y / (2 * np.pi)  # A

            box[0] = X / 10
            self.ostream.print_info(
                f"Box size x: {box[0]:.3f}, minimum CNT radius: {self.CNT_radius_nm:.3f} nm."
            )
            self.ostream.print_info(
                f"Building CNT with X: {X/10:.3f} nm, R: {R/10.:3f} nm")
            self.ostream.flush()
        graphene_xyz = """4

        C        0.00000        1.21944        0.00000
        C        0.70404        0.00000        0.00000
        C        2.11212        0.00000        0.00000
        C        2.81616        1.21944        0.00000"""
        graphene_cell = Molecule.read_xyz_string(graphene_xyz)

        C_chain = topology.addChain()
        C_residue = topology.addResidue("CCC", chain=C_chain)
        CC_atoms = []

        # https://pubs.acs.org/doi/10.1021/jp011344u
        graphene_bond: dict = {
            'type': 'harmonic',
            'force_constant': 4579.87 * 100,
            'equilibrium': cc_eq,
            'comment': 'Graphene bond, cg-cg GAFF parameters'
        }

        graphene_angle: dict = {
            'type': 'harmonic',
            'force_constant': 504.75776,
            'equilibrium': 120,
            'comment': 'Graphene angle, ce-ce-ce GAFF parameters'
        }

        graphene_dihedral: dict = {
            "type": "Fourier",
            "barrier": 25,
            "phase": 180,
            "periodicity": 2,
            "comment": 'Graphene torsion, X -X -c -o GAFF parameters'
        }

        sigma = 0.3851
        epsilon = 0.4396

        bonds = {}
        angles = {}
        dihedrals = {}
        impropers = {}
        carbon_atoms = []

        def cell(_m, _n):
            return 4 * (_m % M * N + _n % N)

        index = 0
        mol_positions = system_mol.get_coordinates_in_angstrom()
        middle_x = (max(mol_positions[:, 0]) + min(mol_positions[:, 0])) / 2
        middle_y = (max(mol_positions[:, 1]) + min(mol_positions[:, 1])) / 2
        min_z = min(mol_positions[:, 2])
        min_y = min(mol_positions[:, 1])
        if self.graphene:
            x_offset = middle_x - X / 2
            y_offset = middle_y - Y / 2
            z_offset = min_z - 5
        else:
            x_offset = middle_x - X / 2
            y_offset = min_y - 0.7 * (R + 1)
            z_offset = min_z - 0.7 * (R + 1)

        for m in range(0, M):
            for n in range(0, N):
                for i, (atomic_number, coord) in enumerate(
                        zip(graphene_cell.get_element_ids(),
                            graphene_cell.get_coordinates_in_angstrom())):

                    coord[0] += m * x_disp
                    coord[1] += n * y_disp
                    if self.CNT:
                        phi = 2 * np.pi * coord[1] / Y
                        coord = np.array(
                            [coord[0], R * math.sin(phi), R * math.cos(phi)]
                        )  # A, X/2 factor is to center the CNT in the middle of the box. The X length is used to get the proper box size
                    coord += np.array([x_offset, y_offset, z_offset])

                    mm_element = mmapp.Element.getByAtomicNumber(atomic_number)
                    name = index
                    index += 1
                    atom = topology.addAtom(f"{name}", mm_element, C_residue)
                    carbon_atoms.append(atom)
                    CC_atoms.append(atom)

                    system_mol.add_atom(int(atomic_number), Point(coord),
                                        'angstrom')
                    self.positions = np.vstack([self.positions, coord])
                    system.addParticle(mm_element.mass)
                    nb_force.addParticle(0, sigma, epsilon)

                # Obtain the starting indices for the current unit cell and all adjecent unit cells
                center = cell(m, n)
                left = cell(m - 1, n)
                right = cell(m + 1, n)
                up = cell(m, n + 1)
                upright = cell(m + 1, n + 1)
                down = cell(m, n - 1)

                # Add bonds and angles
                bonds.update({(center + 0, center + 1): graphene_bond})
                bonds.update({(center + 1, center + 2): graphene_bond})
                bonds.update({(center + 2, center + 3): graphene_bond})
                bonds.update({(center + 3, right + 0): graphene_bond})
                bonds.update({(center + 0, up + 1): graphene_bond})
                bonds.update({(center + 3, up + 2): graphene_bond})

                angles.update({
                    (center + 0, center + 1, center + 2):
                    graphene_angle
                })
                angles.update({
                    (center + 1, center + 2, center + 3):
                    graphene_angle
                })
                angles.update({
                    (up + 1, center + 0, center + 1): graphene_angle
                })
                angles.update({
                    (center + 2, center + 3, up + 2): graphene_angle
                })
                angles.update({(left + 3, center + 0, up + 1): graphene_angle})
                angles.update({
                    (left + 3, center + 0, center + 1): graphene_angle
                })
                angles.update({
                    (center + 0, center + 1, down + 0): graphene_angle
                })
                angles.update({
                    (down + 0, center + 1, center + 2): graphene_angle
                })
                angles.update({
                    (center + 1, center + 2, down + 3): graphene_angle
                })
                angles.update({
                    (center + 3, center + 2, down + 3): graphene_angle
                })
                angles.update({
                    (center + 2, center + 3, right + 0):
                    graphene_angle
                })
                angles.update({(right + 0, center + 3, up + 2): graphene_angle})

                dihedrals.update({
                    (left + 3, center + 0, center + 1, down + 0):
                    graphene_dihedral
                })
                dihedrals.update({
                    (up + 1, center + 0, center + 1, down + 0):
                    graphene_dihedral
                })
                dihedrals.update({
                    (left + 3, center + 0, center + 1, center + 2):
                    graphene_dihedral
                })
                dihedrals.update({
                    (left + 3, center + 0, center + 1, center + 2):
                    graphene_dihedral
                })

                dihedrals.update({
                    (center + 0, center + 1, center + 2, center + 3):
                    graphene_dihedral
                })
                dihedrals.update({
                    (down + 0, center + 1, center + 2, center + 3):
                    graphene_dihedral
                })
                dihedrals.update({
                    (center + 0, center + 1, center + 2, down + 3):
                    graphene_dihedral
                })
                dihedrals.update({
                    (down + 0, center + 1, center + 2, down + 3):
                    graphene_dihedral
                })

                dihedrals.update({
                    (center + 1, center + 2, center + 3, up + 2):
                    graphene_dihedral
                })
                dihedrals.update({
                    (down + 3, center + 2, center + 3, up + 2):
                    graphene_dihedral
                })
                dihedrals.update({
                    (center + 1, center + 2, center + 3, right + 0):
                    graphene_dihedral
                })
                dihedrals.update({
                    (down + 3, center + 2, center + 3, right + 0):
                    graphene_dihedral
                })

                dihedrals.update({
                    (center + 2, center + 3, right + 0, right + 1):
                    graphene_dihedral
                })
                dihedrals.update({
                    (up + 2, center + 3, right + 0, right + 1):
                    graphene_dihedral
                })
                dihedrals.update({
                    (center + 2, center + 3, right + 0, upright + 1):
                    graphene_dihedral
                })
                dihedrals.update({
                    (up + 2, center + 3, right + 0, upright + 1):
                    graphene_dihedral
                })

                dihedrals.update({
                    (left + 3, center + 0, up + 1, up + 0):
                    graphene_dihedral
                })
                dihedrals.update({
                    (center + 1, center + 0, up + 1, up + 0):
                    graphene_dihedral
                })
                dihedrals.update({
                    (left + 3, center + 0, up + 1, up + 2):
                    graphene_dihedral
                })
                dihedrals.update({
                    (center + 1, center + 0, up + 1, up + 2):
                    graphene_dihedral
                })

                dihedrals.update({
                    (center + 2, center + 3, up + 2, up + 1):
                    graphene_dihedral
                })
                dihedrals.update({
                    (center + 2, center + 3, up + 2, up + 3):
                    graphene_dihedral
                })
                dihedrals.update({
                    (right + 0, center + 3, up + 2, up + 1):
                    graphene_dihedral
                })
                dihedrals.update({
                    (right + 0, center + 3, up + 2, up + 3):
                    graphene_dihedral
                })

                impropers.update({
                    (center + 0, left + 3, center + 1, up + 2):
                    graphene_dihedral
                })
                impropers.update({
                    (center + 1, center + 0, down + 0, center + 2):
                    graphene_dihedral
                })
                impropers.update({
                    (center + 2, center + 1, down + 3, center + 3):
                    graphene_dihedral
                })
                impropers.update({
                    (center + 3, center + 2, right + 0, up + 2):
                    graphene_dihedral
                })

        carbon_harmonic_bond_force = mm.HarmonicBondForce()
        carbon_harmonic_bond_force.setName("Carbon harmonic bond")
        carbon_harmonic_bond_force.setUsesPeriodicBoundaryConditions(True)

        carbon_harmonic_angle_force = mm.HarmonicAngleForce()
        carbon_harmonic_angle_force.setName("Carbon harmonic angle")
        carbon_harmonic_angle_force.setUsesPeriodicBoundaryConditions(True)

        carbon_fourier_dihedral_force = mm.PeriodicTorsionForce()
        carbon_fourier_dihedral_force.setName("Carbon improper fourier torsion")
        carbon_fourier_dihedral_force.setUsesPeriodicBoundaryConditions(True)

        for key, bond in bonds.items():
            atom_ids = self._key_to_id(key, CC_atoms)
            self._add_bond(
                carbon_harmonic_bond_force,
                atom_ids,
                bond['equilibrium'],
                bond['force_constant'],
            )

        for key, angle in angles.items():
            atom_ids = self._key_to_id(key, CC_atoms)
            self._add_angle(
                carbon_harmonic_angle_force,
                atom_ids,
                angle['equilibrium'],
                angle['force_constant'],
            )

        for key, dihedral in dihedrals.items():
            atom_ids = self._key_to_id(key, CC_atoms)
            self._add_torsion(carbon_fourier_dihedral_force, dihedral, atom_ids)

        for key, improper in impropers.items():
            atom_ids = self._key_to_id(key, CC_atoms)
            self._add_torsion(carbon_fourier_dihedral_force,
                              improper,
                              atom_ids,
                              improper=True)

        for i in range(len(carbon_atoms)):
            for j in range(len(carbon_atoms)):
                if j > i:
                    nb_force.addException(carbon_atoms[i].index,
                                          carbon_atoms[j].index, 0, 1, 0)
        if not self.no_force_groups:
            carbon_harmonic_bond_force.setForceGroup(EvbForceGroup.CARBON.value)
            carbon_harmonic_angle_force.setForceGroup(
                EvbForceGroup.CARBON.value)
            carbon_fourier_dihedral_force.setForceGroup(
                EvbForceGroup.CARBON.value)
        system.addForce(carbon_harmonic_bond_force)
        system.addForce(carbon_harmonic_angle_force)
        system.addForce(carbon_fourier_dihedral_force)

        return box

    def _add_solvent(self, system, system_mol, solvent, topology, nb_force,
                     neutralize, padding, box):

        assert_msg_critical('openmm' in sys.modules,
                            'openmm is required for EvbSystemBuilder.')

        vlxsysbuilder = SolvationBuilder(ostream=self.ostream)

        mols_per_nm3, density, smiles_code = vlxsysbuilder._solvent_properties(
            solvent)
        vlxsysbuilder.solvate(
            solute=system_mol,
            solvent=solvent,
            padding=padding,
            neutralize=neutralize,
            target_density=density,
            box=box  # A
        )

        self.positions = vlxsysbuilder.system_molecule.get_coordinates_in_angstrom(
        )
        box = [side * 0.1 for side in vlxsysbuilder.box]

        solvents = vlxsysbuilder.solvents
        solvent_counts = vlxsysbuilder.added_solvent_counts
        resnames = ["SOL"] * len(vlxsysbuilder.added_solvent_counts)
        if hasattr(vlxsysbuilder, "added_counterions"):
            if vlxsysbuilder.added_counterions > 0:
                solvents = solvents + [vlxsysbuilder._counterion_molecules()
                                       ] * vlxsysbuilder.added_counterions
                solvent_counts = solvent_counts + [
                    vlxsysbuilder.added_counterions
                ]
                resnames = resnames + ["ION"]

        for i, (solvent_count, vlx_solvent_molecule,
                resname) in enumerate(zip(solvent_counts, solvents, resnames)):

            num_solvent_atoms_per_molecule = vlx_solvent_molecule.number_of_atoms(
            )
            elements = vlx_solvent_molecule.get_labels()

            harmonic_bond_force = mm.HarmonicBondForce()
            harmonic_bond_force.setName(f"Solvent {i} harmonic bond")
            harmonic_angle_force = mm.HarmonicAngleForce()
            harmonic_angle_force.setName(f"Solvent {i} harmonic angle")
            fourier_force = mm.PeriodicTorsionForce()
            fourier_force.setName(f"Solvent {i} proper fourier torsion")
            fourier_imp_force = mm.PeriodicTorsionForce()
            fourier_imp_force.setName(f"Solvent {i} improper fourier torsion")

            # Loop over all solvent molecules
            solvent_nb_atom_count = 0
            solvent_system_atom_count = 0
            solvent_chain = topology.addChain()
            solvent_ff = MMForceFieldGenerator(ostream=self.ostream)

            solvent_ff.create_topology(vlx_solvent_molecule,
                                       water_model=self.water_model)

            atom_types = [atom['type'] for atom in solvent_ff.atoms.values()]
            if 'ow' in atom_types and 'hw' in atom_types and len(
                    atom_types) == 3:
                water_model = get_water_parameters()[self.water_model]
                for atom_id, atom in solvent_ff.atoms.items():
                    solvent_ff.atoms[atom_id] = copy.copy(
                        water_model[atom['type']])
                    solvent_ff.atoms[atom_id]['name'] = solvent_ff.atoms[
                        atom_id]['name'][0] + str(atom_id)
                for bond_id in solvent_ff.bonds.keys():
                    solvent_ff.bonds[bond_id] = water_model['bonds']
                for ang_id in solvent_ff.angles.keys():
                    solvent_ff.angles[ang_id] = water_model['angles']

            self.solvent_atom_ids = []
            for i in range(solvent_count):

                solvent_residue = topology.addResidue(name=resname,
                                                      chain=solvent_chain)
                local_solvent_atoms = []
                # Loop over all atoms in the solvent molecule
                for j in range(num_solvent_atoms_per_molecule):
                    # Figure out the element of the atom
                    mm_element = mmapp.Element.getBySymbol(elements[j])
                    name = f"{elements[j]}{j}"
                    # add the atom to the topology
                    solvent_atom = topology.addAtom(name, mm_element,
                                                    solvent_residue)
                    local_solvent_atoms.append(solvent_atom)
                    self.solvent_atom_ids.append(solvent_atom.index)
                    # Add the atom as a particle to the system
                    system.addParticle(mm_element.mass)
                    solvent_system_atom_count += 1
                    # Add this atom to the nonbonded force

                # add all bonded interactions in this molecule
                for key, bond in solvent_ff.bonds.items():
                    atom_ids = self._key_to_id(key, local_solvent_atoms)
                    self._add_bond(
                        harmonic_bond_force,
                        atom_ids,
                        bond['equilibrium'],
                        bond['force_constant'],
                    )

                for key, angle in solvent_ff.angles.items():
                    atom_ids = self._key_to_id(key, local_solvent_atoms)
                    self._add_angle(
                        harmonic_angle_force,
                        atom_ids,
                        angle['equilibrium'],
                        angle['force_constant'],
                    )

                for key, dihedral in solvent_ff.dihedrals.items():
                    atom_ids = self._key_to_id(key, local_solvent_atoms)
                    self._add_torsion(fourier_force, dihedral, atom_ids)
                for key, dihedral in solvent_ff.impropers.items():
                    atom_ids = self._key_to_id(key, local_solvent_atoms)
                    self._add_torsion(
                        fourier_imp_force,
                        dihedral,
                        atom_ids,
                        improper=True,
                    )

                exceptions = self._create_exceptions_from_bonds(
                    solvent_ff.atoms, solvent_ff.bonds)

                for key, atom in solvent_ff.atoms.items():

                    sigma = atom["sigma"]
                    if sigma == 0:
                        if atom["epsilon"] == 0:
                            sigma = 1
                        else:
                            raise ValueError(
                                "Sigma is 0 while epsilon is not, which will cause division by 0"
                            )
                    solvent_nb_atom_count += 1
                    nb_force.addParticle(atom["charge"], atom["sigma"],
                                         atom["epsilon"])
                # #Loop over all atoms, and check if their id's are part of any exceptions
                for i in solvent_ff.atoms.keys():
                    for j in solvent_ff.atoms.keys():
                        if i < j:
                            key = (i, j)
                            if key in exceptions.keys():
                                epsilon = exceptions[key]["epsilon"]
                                sigma = exceptions[key]["sigma"]
                                qq = exceptions[key]["qq"]

                                nb_force.addException(
                                    local_solvent_atoms[i].index,
                                    local_solvent_atoms[j].index,
                                    qq,
                                    sigma,
                                    epsilon,
                                )

            if harmonic_bond_force.getNumBonds() > 0:
                if not self.no_force_groups:
                    harmonic_bond_force.setForceGroup(
                        EvbForceGroup.SOLVENT.value)
                system.addForce(harmonic_bond_force)
            if harmonic_angle_force.getNumAngles() > 0:
                if not self.no_force_groups:
                    harmonic_angle_force.setForceGroup(
                        EvbForceGroup.SOLVENT.value)
                system.addForce(harmonic_angle_force)
            if fourier_force.getNumTorsions() > 0:
                if not self.no_force_groups:
                    fourier_force.setForceGroup(EvbForceGroup.SOLVENT.value)
                system.addForce(fourier_force)
            if fourier_imp_force.getNumTorsions() > 0:
                if not self.no_force_groups:
                    fourier_imp_force.setForceGroup(EvbForceGroup.SOLVENT.value)
                system.addForce(fourier_imp_force)
            self.ostream.print_info(
                f"Added {solvent_nb_atom_count} atoms to the nonbonded force and {solvent_system_atom_count} atoms to the system"
            )
            self.ostream.flush()
        return box

    def _add_barostat(self, system):
        if not (self.CNT or self.graphene):
            barostat = mm.MonteCarloBarostat(
                self.pressure * mmunit.bar,  # type: ignore
                self.temperature * mmunit.kelvin,  # type: ignore
            )
        else:
            barostat = mm.MonteCarloFlexibleBarostat(
                self.pressure * mmunit.bar,  # type: ignore
                self.temperature * mmunit.kelvin,  # type: ignore
            )
        if not self.no_force_groups:
            barostat.setForceGroup(EvbForceGroup.BAROSTAT.value)
        system.addForce(barostat)
        return barostat

    # ------------------------------------------------------------------
    # GB implicit-solvent helpers (ported from openmm.app.internal so
    # only the public CustomGBForce API is required).
    # ------------------------------------------------------------------

    @staticmethod
    def _gb_bonded_atoms(topology):
        from collections import defaultdict
        bonds = defaultdict(list)
        for a1, a2 in topology.bonds():
            bonds[a1].append(a2)
            bonds[a2].append(a1)
        return bonds

    @staticmethod
    def _gb_is_carboxylate_O(atom, bonds):
        if atom.element is None or atom.element.symbol != 'O':
            return False
        nbrs = bonds[atom]
        if len(nbrs) != 1 or nbrs[0].element is None or nbrs[0].element.symbol != 'C':
            return False
        c_nbrs = bonds[nbrs[0]]
        if len(c_nbrs) != 3:
            return False
        for a3 in c_nbrs:
            if a3 is not atom and a3.element is not None and a3.element.symbol == 'O':
                return True
        return False

    @staticmethod
    def _gb_bondi_radii(topology):
        radii_map = {'C': 1.7, 'H': 1.2, 'D': 1.2, 'N': 1.55, 'O': 1.5,
                     'F': 1.5, 'Si': 2.1, 'P': 1.85, 'S': 1.8, 'Cl': 1.7}
        radii = []
        for atom in topology.atoms():
            sym = atom.element.symbol if atom.element is not None else ''
            radii.append(radii_map.get(sym, 1.5))
        return radii  # in Angstrom

    @staticmethod
    def _gb_mbondi_radii(topology):
        bonds = EvbSystemBuilder._gb_bonded_atoms(topology)
        const_map = {'N': 1.55, 'O': 1.5, 'F': 1.5, 'Si': 2.1,
                     'P': 1.85, 'S': 1.8, 'Cl': 1.7}
        radii = []
        for atom in topology.atoms():
            sym = atom.element.symbol if atom.element is not None else ''
            if sym in ('H', 'D'):
                nbrs = bonds[atom]
                nsym = nbrs[0].element.symbol if nbrs and nbrs[0].element else ''
                if nsym in ('C', 'N'):
                    r = 1.3
                elif nsym in ('O', 'S'):
                    r = 0.8
                else:
                    r = 1.2
            elif sym == 'C':
                r = 1.7
            else:
                r = const_map.get(sym, 1.5)
            radii.append(r)
        return radii  # in Angstrom

    @staticmethod
    def _gb_mbondi2_radii(topology):
        bonds = EvbSystemBuilder._gb_bonded_atoms(topology)
        const_map = {'N': 1.55, 'O': 1.5, 'F': 1.5, 'Si': 2.1,
                     'P': 1.85, 'S': 1.8, 'Cl': 1.7}
        radii = []
        for atom in topology.atoms():
            sym = atom.element.symbol if atom.element is not None else ''
            if sym in ('H', 'D'):
                nbrs = bonds[atom]
                nsym = nbrs[0].element.symbol if nbrs and nbrs[0].element else ''
                r = 1.3 if nsym == 'N' else 1.2
            elif sym == 'C':
                r = 1.7
            else:
                r = const_map.get(sym, 1.5)
            radii.append(r)
        return radii  # in Angstrom

    @staticmethod
    def _gb_mbondi3_radii(topology):
        bonds = EvbSystemBuilder._gb_bonded_atoms(topology)
        radii = EvbSystemBuilder._gb_mbondi2_radii(topology)
        for i, atom in enumerate(topology.atoms()):
            if EvbSystemBuilder._gb_is_carboxylate_O(atom, bonds):
                radii[i] = 1.4
            elif atom.residue.name == 'ARG':
                name = atom.name
                if name.startswith('HH') or name.startswith('HE'):
                    radii[i] = 1.17
        return radii  # in Angstrom

    @staticmethod
    def _gb_neck_tables(unique_ors, offset):
        """Return (d0_table, m0_table) reduced to the unique-radius subspace."""
        n = len(unique_ors)
        positions = [(r + offset - 0.1) * 200 for r in unique_ors]
        idx1 = [0] * n
        idx2 = [0] * n
        w1 = [0.0] * n
        w2 = [0.0] * n
        for i, p in enumerate(positions):
            if p <= 0:
                w1[i] = 1.0
            elif p >= 20:
                idx1[i] = 20
                w1[i] = 1.0
            else:
                idx1[i] = int(math.floor(p))
                idx2[i] = idx1[i] + 1
                w1[i] = idx2[i] - p
                w2[i] = 1.0 - w1[i]

        def _reduce(full):
            table = []
            for i in range(n):
                for j in range(n):
                    v = (w1[i] * w1[j] * full[idx1[i] * 21 + idx1[j]] +
                         w1[i] * w2[j] * full[idx1[i] * 21 + idx2[j]] +
                         w2[i] * w1[j] * full[idx2[i] * 21 + idx1[j]] +
                         w2[i] * w2[j] * full[idx2[i] * 21 + idx2[j]])
                    table.append(v)
            return table

        return _reduce(_GB_D0), _reduce(_GB_M0)

    @staticmethod
    def _gb_screen(atom, idx):
        sym = atom.element.symbol if atom.element is not None else None
        return _GB_SCREEN_BY_ELEMENT_SYMBOL.get(sym, _GB_SCREEN_DEFAULT)[idx]

    @staticmethod
    def _gb_energy_terms(force, sv_diel, sl_diel, sa, offset):
        params = (f'; solventDielectric={sv_diel:.16g}'
                  f'; soluteDielectric={sl_diel:.16g}'
                  f'; kappa=0; offset={offset:.16g}')
        force.addEnergyTerm(
            '-0.5*138.935485*(1/soluteDielectric-1/solventDielectric)*charge^2/B' + params,
            mm.CustomGBForce.SingleParticle)
        if sa == 'ACE':
            force.addEnergyTerm(
                '28.3919551*(radius+0.14)^2*(radius/B)^6; radius=or+offset' + params,
                mm.CustomGBForce.SingleParticle)
        force.addEnergyTerm(
            '-138.935485*(1/soluteDielectric-1/solventDielectric)*charge1*charge2/f;'
            'f=sqrt(r^2+B1*B2*exp(-r^2/(4*B1*B2)))' + params,
            mm.CustomGBForce.ParticlePairNoExclusions)

    def _add_implicit_solvent(self, system, topology, nb_force):
        """Add a GB implicit solvent force to the system.

        Builds the CustomGBForce directly using the public OpenMM API,
        without depending on openmm.app.internal. Radii tables and formulas
        mirror those in the internal module exactly.

        :param system: the OpenMM System to add the force to.
        :param topology: the OpenMM Topology (all atoms must already be added).
        :param nb_force: the NonbondedForce already present in the system.
        """
        assert_msg_critical('openmm' in sys.modules,
                            'openmm is required for EvbSystemBuilder.')

        valid_models = ('hct', 'obc1', 'obc2', 'gbn', 'gbn2')
        model_key = self.implicit_solvent_model.lower()
        assert_msg_critical(
            model_key in valid_models,
            f'EvbSystemBuilder: unknown implicit_solvent_model '
            f'"{self.implicit_solvent_model}". '
            f'Valid options are: {list(valid_models)}')

        sv = self.solvent_dielectric
        sl = self.solute_dielectric

        _I_simple = (
            'select(step(r+sr2-or1), 0.5*(1/L-1/U+0.25*(r-sr2^2/r)*(1/(U^2)-1/(L^2))+0.5*log(L/U)/r), 0);'
            'U=r+sr2;'
            'L=max(or1, D);'
            'D=abs(r-sr2)')

        gb_force = mm.CustomGBForce()

        if model_key == 'hct':
            offset = 0.009
            radii_ang = self._gb_mbondi_radii(topology)
            gb_force.addPerParticleParameter('charge')
            gb_force.addPerParticleParameter('or')
            gb_force.addPerParticleParameter('sr')
            gb_force.addComputedValue('I', _I_simple,
                                      mm.CustomGBForce.ParticlePairNoExclusions)
            gb_force.addComputedValue('B', '1/(1/or-I)',
                                      mm.CustomGBForce.SingleParticle)
            self._gb_energy_terms(gb_force, sv, sl, 'ACE', offset)
            for i, atom in enumerate(topology.atoms()):
                charge = nb_force.getParticleParameters(i)[0]
                or_ = radii_ang[i] / 10.0 - offset
                sr = self._gb_screen(atom, 0) * or_
                gb_force.addParticle([charge, or_, sr])

        elif model_key == 'obc1':
            offset = 0.009
            radii_ang = self._gb_mbondi2_radii(topology)
            gb_force.addPerParticleParameter('charge')
            gb_force.addPerParticleParameter('or')
            gb_force.addPerParticleParameter('sr')
            gb_force.addComputedValue('I', _I_simple,
                                      mm.CustomGBForce.ParticlePairNoExclusions)
            gb_force.addComputedValue(
                'B',
                '1/(1/or-tanh(0.8*psi+2.909125*psi^3)/radius);'
                'psi=I*or; radius=or+offset; offset=0.009',
                mm.CustomGBForce.SingleParticle)
            self._gb_energy_terms(gb_force, sv, sl, 'ACE', offset)
            for i, atom in enumerate(topology.atoms()):
                charge = nb_force.getParticleParameters(i)[0]
                or_ = radii_ang[i] / 10.0 - offset
                sr = self._gb_screen(atom, 0) * or_
                gb_force.addParticle([charge, or_, sr])

        elif model_key == 'obc2':
            offset = 0.009
            radii_ang = self._gb_mbondi2_radii(topology)
            gb_force.addPerParticleParameter('charge')
            gb_force.addPerParticleParameter('or')
            gb_force.addPerParticleParameter('sr')
            gb_force.addComputedValue('I', _I_simple,
                                      mm.CustomGBForce.ParticlePairNoExclusions)
            gb_force.addComputedValue(
                'B',
                '1/(1/or-tanh(psi-0.8*psi^2+4.85*psi^3)/radius);'
                'psi=I*or; radius=or+offset; offset=0.009',
                mm.CustomGBForce.SingleParticle)
            self._gb_energy_terms(gb_force, sv, sl, 'ACE', offset)
            for i, atom in enumerate(topology.atoms()):
                charge = nb_force.getParticleParameters(i)[0]
                or_ = radii_ang[i] / 10.0 - offset
                sr = self._gb_screen(atom, 0) * or_
                gb_force.addParticle([charge, or_, sr])

        elif model_key == 'gbn':
            offset = 0.009
            radii_ang = self._gb_bondi_radii(topology)
            atoms = list(topology.atoms())
            ors = [radii_ang[i] / 10.0 - offset for i in range(len(atoms))]
            unique_ors = sorted(set(ors))
            or_to_idx = {r: i for i, r in enumerate(unique_ors)}
            n = len(unique_ors)
            d0_tab, m0_tab = self._gb_neck_tables(unique_ors, offset)

            gb_force.addPerParticleParameter('charge')
            gb_force.addPerParticleParameter('or')
            gb_force.addPerParticleParameter('sr')
            gb_force.addPerParticleParameter('radindex')
            gb_force.addTabulatedFunction(
                'getd0', mm.Discrete2DFunction(n, n, d0_tab))
            gb_force.addTabulatedFunction(
                'getm0', mm.Discrete2DFunction(n, n, m0_tab))
            gb_force.addComputedValue(
                'I',
                'Ivdw+neckScale*Ineck;'
                'Ineck=step(radius1+radius2+neckCut-r)*getm0(radindex1,radindex2)/'
                '(1+100*(r-getd0(radindex1,radindex2))^2+'
                '0.3*1000000*(r-getd0(radindex1,radindex2))^6);'
                'Ivdw=select(step(r+sr2-or1), 0.5*(1/L-1/U+0.25*(r-sr2^2/r)*(1/(U^2)-1/(L^2))+0.5*log(L/U)/r), 0);'
                'U=r+sr2;'
                'L=max(or1, D);'
                'D=abs(r-sr2);'
                'radius1=or1+offset; radius2=or2+offset;'
                'neckScale=0.361825; neckCut=0.68; offset=0.009',
                mm.CustomGBForce.ParticlePairNoExclusions)
            gb_force.addComputedValue(
                'B',
                '1/(1/or-tanh(1.09511284*psi-1.907992938*psi^2+2.50798245*psi^3)/radius);'
                'psi=I*or; radius=or+offset; offset=0.009',
                mm.CustomGBForce.SingleParticle)
            self._gb_energy_terms(gb_force, sv, sl, 'ACE', offset)
            for i, atom in enumerate(atoms):
                charge = nb_force.getParticleParameters(i)[0]
                or_ = ors[i]
                sr = self._gb_screen(atom, 1) * or_
                gb_force.addParticle([charge, or_, sr, or_to_idx[or_]])

        elif model_key == 'gbn2':
            offset = 0.0195141
            radii_ang = self._gb_mbondi3_radii(topology)
            atoms = list(topology.atoms())
            ors = [radii_ang[i] / 10.0 - offset for i in range(len(atoms))]
            unique_ors = sorted(set(ors))
            or_to_idx = {r: i for i, r in enumerate(unique_ors)}
            n = len(unique_ors)
            d0_tab, m0_tab = self._gb_neck_tables(unique_ors, offset)

            gb_force.addPerParticleParameter('charge')
            gb_force.addPerParticleParameter('or')
            gb_force.addPerParticleParameter('sr')
            gb_force.addPerParticleParameter('alpha')
            gb_force.addPerParticleParameter('beta')
            gb_force.addPerParticleParameter('gamma')
            gb_force.addPerParticleParameter('radindex')
            gb_force.addTabulatedFunction(
                'getd0', mm.Discrete2DFunction(n, n, d0_tab))
            gb_force.addTabulatedFunction(
                'getm0', mm.Discrete2DFunction(n, n, m0_tab))
            gb_force.addComputedValue(
                'I',
                'Ivdw+neckScale*Ineck;'
                'Ineck=step(radius1+radius2+neckCut-r)*getm0(radindex1,radindex2)/'
                '(1+100*(r-getd0(radindex1,radindex2))^2+'
                '0.3*1000000*(r-getd0(radindex1,radindex2))^6);'
                'Ivdw=select(step(r+sr2-or1), 0.5*(1/L-1/U+0.25*(r-sr2^2/r)*(1/(U^2)-1/(L^2))+0.5*log(L/U)/r), 0);'
                'U=r+sr2;'
                'L=max(or1, D);'
                'D=abs(r-sr2);'
                'radius1=or1+offset; radius2=or2+offset;'
                'neckScale=0.826836; neckCut=0.68; offset=0.0195141',
                mm.CustomGBForce.ParticlePairNoExclusions)
            gb_force.addComputedValue(
                'B',
                '1/(1/or-tanh(alpha*psi-beta*psi^2+gamma*psi^3)/radius);'
                'psi=I*or; radius=or+offset; offset=0.0195141',
                mm.CustomGBForce.SingleParticle)
            self._gb_energy_terms(gb_force, sv, sl, 'ACE', offset)
            for i, atom in enumerate(atoms):
                charge = nb_force.getParticleParameters(i)[0]
                or_ = ors[i]
                is_nucleic = atom.residue.name in _GB_NUCLEIC_RESIDUES
                sym = atom.element.symbol if atom.element is not None else ''
                if is_nucleic:
                    screen_idx = 3
                    abg = _GBN2_ATOM_PARAMS_NUCLEIC.get(sym, _GBN2_DEFAULT_PARAMS)
                else:
                    screen_idx = 2
                    abg = _GBN2_ATOM_PARAMS.get(sym, _GBN2_DEFAULT_PARAMS)
                sr = self._gb_screen(atom, screen_idx) * or_
                gb_force.addParticle(
                    [charge, or_, sr, abg[0], abg[1], abg[2], or_to_idx[or_]])

        gb_force.setNonbondedMethod(mm.CustomGBForce.NoCutoff)

        # Required by OpenMM when GB is active: set reaction-field dielectric
        # of the NonbondedForce to 1 so it does not double-count dielectric
        # screening.
        nb_force.setReactionFieldDielectric(1)

        system.addForce(gb_force)
        self.ostream.print_info(
            f'Added implicit solvent ({self.implicit_solvent_model}) with '
            f'solvent_dielectric={self.solvent_dielectric}, '
            f'solute_dielectric={self.solute_dielectric}')
        self.ostream.flush()

    def _interpolate_system(self, system, lambda_vec, nb_force, E_field_force):
        systems = {}

        # Locate the GB force once so we can keep its per-particle charges in
        # sync with the interpolated NB charges before each deepcopy.
        gb_force = None
        if self.implicit_solvent_model is not None:
            gb_candidates = [
                f for f in system.getForces()
                if isinstance(f, mm.CustomGBForce)
            ]
            if gb_candidates:
                gb_force = gb_candidates[0]

        for lam in lambda_vec:
            total_charge = 0
            if not self.no_reactant:
                # Set interpolated nonbonded parameters for the system-environment interaction
                for i, (reactant_atom, product_atom) in enumerate(
                        zip(self.reactant.atoms.values(),
                            self.product.atoms.values())):
                    if reactant_atom.get('pdb') == 'cap':
                        continue
                    charge = (1 - lam) * reactant_atom[
                        "charge"] + lam * product_atom["charge"]
                    total_charge += charge
                    sigma = (1 - lam) * reactant_atom[
                        "sigma"] + lam * product_atom["sigma"]
                    epsilon = (1 - lam) * reactant_atom[
                        "epsilon"] + lam * product_atom["epsilon"]
                    if sigma == 0:
                        if epsilon == 0:
                            sigma = 1
                        else:
                            raise ValueError(
                                "Sigma is 0 while epsilon is not, which will cause division by 0"
                            )

                    nb_force.setParticleParameters(self.reaction_atoms[i].index,
                                                   charge, sigma, epsilon)

                    # Mirror the interpolated charge into the GB force so each
                    # deepcopy carries the correct lambda-dependent charge.
                    if gb_force is not None:
                        atom_index = self.reaction_atoms[i].index
                        gb_params = list(
                            gb_force.getParticleParameters(atom_index))
                        gb_params[0] = charge
                        gb_force.setParticleParameters(atom_index, gb_params)

            # Add the interpolated charge to the E_field for both the system and environment
            for i in range(system.getNumParticles()):
                charge = nb_force.getParticleParameters(i)[0]
                if E_field_force is not None:
                    E_field_force.setParticleParameters(i, i, [charge])

            if not round(total_charge, 5).is_integer():
                self.ostream.print_warning(
                    f"Total charge for lambda {lam} is {total_charge} and is not a whole number"
                )

            new_system = copy.deepcopy(system)
            if lam == 0:
                rea_system = copy.deepcopy(system)
            if lam == 1:
                pro_system = copy.deepcopy(system)
            # Add the bonded forces for the reaction system
            if not self.no_reactant:
                self._add_reaction_forces(new_system, lam)
            systems[lam] = new_system

        self._add_reaction_forces(rea_system, 0, pes=True)
        self._add_reaction_forces(pro_system, 1, pes=True)

        if self.decompose_nb is not None:
            if self.solvent:
                self.ostream.print_info(
                    f"Adding nonbonded force decomposition reporting for particles {self.decompose_nb}"
                )
                systems.update(self._add_nb_decompositions(rea_system, 'rea'))
                systems.update(self._add_nb_decompositions(pro_system, 'pro'))
            else:
                self.ostream.print_info(
                    "Skipping nonbonded force decompositions")
            self.ostream.flush()

        if self.decompose_bonded:
            rea_bond_decomp = self._add_bonded_decompositions(rea_system)
            pro_bond_decomp = self._add_bonded_decompositions(pro_system)
            systems['reactant_bonded_decomp'] = rea_bond_decomp
            systems['product_bonded_decomp'] = pro_bond_decomp

        # rea_system = self._split_nb_force(rea_system)
        # pro_system = self._split_nb_force(pro_system)

        systems['reactant'] = rea_system
        systems['product'] = pro_system
        self.ostream.flush()
        return systems

    def _add_reaction_forces(self, system, lam, pes=False):

        assert_msg_critical('openmm' in sys.modules,
                            'openmm is required for EvbSystemBuilder.')
        static_bonded_harmonic = self._create_static_harmonic_bond_forces(lam)
        dynamic_bonded_harmonic = self._create_dynamic_harmonic_bond_forces(lam)
        morse = self._create_morse_force(lam, model_broken=False)
        angle = self._create_harmonic_angle_forces(lam, model_broken=False)

        # angle, angle_integration = self._create_angle_forces(lam)
        torsion_forces = self._create_proper_torsion_forces(lam, pes=pes)
        improper = self._create_improper_torsion_forces(lam)

        if not pes:
            system.addForce(dynamic_bonded_harmonic)
            syslj, syscoul = self._create_nonbonded_forces(
                lam,
                lj_soft_core=self.soft_core_lj_int,
                coul_soft_core=self.soft_core_coulomb_int,
                broken_exceptions=True,
            )
            syslj.setForceGroup(EvbForceGroup.SYSLJ_STATIC.value)
            syscoul.setForceGroup(EvbForceGroup.SYSCOUL_STATIC.value)
            system.addForce(syslj)
            system.addForce(syscoul)
        else:
            system.addForce(morse)
            syslj_static, syscoul_static = self._create_nonbonded_forces(
                lam,
                lj_soft_core=self.soft_core_lj_pes_static,
                coul_soft_core=self.soft_core_coulomb_pes_static,
                exclude_changing_bonds=True,
            )
            syslj_static.setForceGroup(EvbForceGroup.SYSLJ_STATIC.value)
            syscoul_static.setForceGroup(EvbForceGroup.SYSCOUL_STATIC.value)
            syslj_static.setName('Reaction internal LJ static')
            syscoul_static.setName('Reaction internal Coul static')
            system.addForce(syslj_static)
            system.addForce(syscoul_static)
            syslj_dynamic, syscoul_dynamic = self._create_nonbonded_forces(
                lam,
                lj_soft_core=self.soft_core_lj_pes_dynamic,
                coul_soft_core=self.soft_core_coulomb_pes_dynamic,
                only_changing_bonds=True,
            )
            syslj_dynamic.setForceGroup(EvbForceGroup.SYSLJ_DYNAMIC.value)
            syscoul_dynamic.setForceGroup(EvbForceGroup.SYSCOUL_DYNAMIC.value)
            syslj_dynamic.setName('Reaction internal LJ dynamic')
            syscoul_dynamic.setName('Reaction internal Coul dynamic')
            system.addForce(syslj_dynamic)
            system.addForce(syscoul_dynamic)

        bond_constraint, constant_force, angle_constraint, torsion_constraint = self._create_constraint_forces(
            lam)

        system.addForce(static_bonded_harmonic)
        system.addForce(angle)
        for torsion_force in torsion_forces:
            system.addForce(torsion_force)
        system.addForce(improper)
        system.addForce(bond_constraint)
        system.addForce(constant_force)
        system.addForce(angle_constraint)
        system.addForce(torsion_constraint)

        return system

    def _add_bonded_decompositions(self, system):
        system = copy.deepcopy(system)
        bonded_fgs = [
            EvbForceGroup.REA_HARM_BOND_STATIC.
            value,  # Bonded forces for the reaction atoms
            EvbForceGroup.REA_HARM_BOND_DYNAMIC.
            value,  # Static bonded forces for the reaction atoms
            EvbForceGroup.REA_MORSE_BOND.value,
            EvbForceGroup.REA_ANGLE.value,
            EvbForceGroup.REA_TORSION.value,
            EvbForceGroup.REA_IMP.value,
        ]
        to_remove = []
        for i, force in enumerate(system.getForces()):
            if force.getForceGroup() not in bonded_fgs:
                to_remove.append(i)
        for i in reversed(to_remove):
            system.removeForce(i)

        for force in system.getForces():
            if force.getForceGroup() == EvbForceGroup.REA_MORSE_BOND.value:
                for i in range(force.getNumBonds()):
                    p1, p2, (D, a, r) = force.getBondParameters(i)
                    force.setBondParameters(i, p1, p2, [0, 0, r])
            if force.getForceGroup(
            ) == EvbForceGroup.REA_HARM_BOND_STATIC.value or force.getForceGroup(
            ) == EvbForceGroup.REA_HARM_BOND_DYNAMIC.value:
                for i in range(force.getNumBonds()):
                    p1, p2, r, k = force.getBondParameters(i)
                    force.setBondParameters(i, p1, p2, r, 0)
            if force.getForceGroup() == EvbForceGroup.REA_ANGLE.value:
                for i in range(force.getNumAngles()):
                    p1, p2, p3, theta, k = force.getAngleParameters(i)
                    force.setAngleParameters(i, p1, p2, p3, theta, 0)
            if force.getForceGroup(
            ) == EvbForceGroup.REA_TORSION.value or force.getForceGroup(
            ) == EvbForceGroup.REA_IMP.value:
                for i in range(force.getNumTorsions()):
                    p1, p2, p3, p4, periodicity, phase, barrier = force.getTorsionParameters(
                        i)
                    force.setTorsionParameters(i, p1, p2, p3, p4, periodicity,
                                               phase, 0)

        return system

    def _add_nb_decompositions(self, system, state_name):
        systems = {}
        nbforce = copy.deepcopy([
            force for force in system.getForces()
            if isinstance(force, mm.NonbondedForce)
        ][0])
        nbforce.setExceptionsUsePeriodicBoundaryConditions(True)

        # forces for solvent solvent interactions
        solcoul = copy.deepcopy(nbforce)
        sollj = copy.deepcopy(nbforce)

        # remove all nonbonded parameters for the compound
        for i, _ in enumerate(self.reactant.atoms.values()):
            atom_id = self.reaction_atoms[i].index
            solcoul.setParticleParameters(atom_id, 0, 1, 0)
            sollj.setParticleParameters(atom_id, 0, 1, 0)

        # set the charges or epsilons of all solvent atoms to 0
        for atom_id in self.solvent_atom_ids:
            charge, sigma, epsilon = nbforce.getParticleParameters(atom_id)
            sollj.setParticleParameters(atom_id, 0, sigma, epsilon)
            solcoul.setParticleParameters(atom_id, charge, 1, 0)

        solcoul.setName('Coul solvent')
        coul_system = copy.deepcopy(system)
        self._remove_forces(coul_system)
        coul_system.addForce(solcoul)

        sollj.setName('LJ solvent')
        lj_system = copy.deepcopy(system)
        self._remove_forces(lj_system)
        lj_system.addForce(sollj)

        systems.update({
            f"decomp_{state_name}_solvent_Coul": coul_system,
            f"decomp_{state_name}_solvent_LJ": lj_system
        })

        # Remove all solvent solvent interactions in the other forces

        for to_decompose in self.decompose_nb:

            lj_dec = copy.deepcopy(nbforce)
            coul_dec = copy.deepcopy(nbforce)

            for i, _ in enumerate(self.reactant.atoms.values()):
                atom_id = self.reaction_atoms[i].index
                charge, sigma, epsilon = nbforce.getParticleParameters(atom_id)

                if i in to_decompose:
                    # Setting all system-solvent interactions as exceptions and defaulting everything to 0 is easier
                    # than setting all solvent-solvent interactions as exceptions with original parameters
                    for solvent_atom in self.solvent_atom_ids:
                        solv_charge, solv_sigma, solv_epsilon = nbforce.getParticleParameters(
                            solvent_atom)
                        qq = charge * solv_charge
                        sig = 0.5 * (solv_sigma + sigma)
                        eps = math.sqrt((epsilon * solv_epsilon).value_in_unit(
                            mmunit.kilojoule_per_mole**2))
                        lj_dec.addException(atom_id, solvent_atom, 0, sig, eps)
                        coul_dec.addException(atom_id, solvent_atom, qq, 1, 0)

                # Everything is handeled by the exceptions, so the default parameters can be set to 0
                lj_dec.setParticleParameters(atom_id, 0, 1, 0)
                coul_dec.setParticleParameters(atom_id, 0, 1, 0)

            # Since all system-solvent interactions are added as exceptions, all the parameters can be set to 0
            for solvent_id in self.solvent_atom_ids:
                lj_dec.setParticleParameters(solvent_id, 0, 1, 0)
                coul_dec.setParticleParameters(solvent_id, 0, 1, 0)
            name = "-".join(str(x) for x in to_decompose)

            coul_system = copy.deepcopy(system)
            self._remove_forces(coul_system)
            coul_system.addForce(coul_dec)

            lj_system = copy.deepcopy(system)
            self._remove_forces(lj_system)
            lj_system.addForce(lj_dec)

            systems.update({
                f"decomp_{state_name}_{name}_Coul": coul_system,
                f"decomp_{state_name}_{name}_LJ": lj_system
            })

        return systems

    @staticmethod
    def _split_nb_force(system):
        nb_force = [
            force for force in system.getForces()
            if isinstance(force, mm.NonbondedForce)
        ][0]
        # nb_force.setNonbondedMethod(mm.NonbondedForce.CutoffNonPeriodic)
        coul_force = nb_force
        lj_force = copy.copy(nb_force)
        coul_force.setName('Solvent coul')
        coul_force.setForceGroup(EvbForceGroup.SOL_COUL.value)
        lj_force.setName('Solvent lj')
        lj_force.setForceGroup(EvbForceGroup.SOL_LJ.value)

        for i in range(nb_force.getNumParticles()):
            charge, sigma, epsilon = nb_force.getParticleParameters(i)
            coul_force.setParticleParameters(i, charge, 1, 0)
            lj_force.setParticleParameters(i, 0, sigma, epsilon)

        for i in range(nb_force.getNumExceptions()):
            p1, p2, charge, sigma, epsilon = nb_force.getExceptionParameters(i)
            coul_force.setExceptionParameters(i, p1, p2, charge, 1, 0)
            lj_force.setExceptionParameters(i, p1, p2, 0, sigma, epsilon)

        system.addForce(lj_force)

    @staticmethod
    def _remove_forces(system):
        for i in range(system.getNumForces()):
            system.removeForce(0)

    def _add_E_field(self, system, E_field):

        assert_msg_critical('openmm' in sys.modules,
                            'openmm is required for EvbSystemBuilder.')

        E_field_force = mm.CustomExternalForce("-q*(Ex*x+Ey*y+Ez*z)")
        E_field_force.addGlobalParameter("Ex", E_field[0])
        E_field_force.addGlobalParameter("Ey", E_field[1])
        E_field_force.addGlobalParameter("Ez", E_field[2])
        E_field_force.addPerParticleParameter("q")
        E_field_force.setName("Electric field")
        system.addForce(E_field_force)
        for i in range(system.getNumParticles()):
            E_field_force.addParticle(
                i,
                [0])  # Actual charges are lambda dependent, and are set later
        if not self.no_force_groups:
            E_field_force.setForceGroup(EvbForceGroup.E_FIELD.value)
        return E_field_force

    def _create_static_harmonic_bond_forces(self, lam):

        assert_msg_critical('openmm' in sys.modules,
                            'openmm is required for EvbSystemBuilder.')

        harmonic_force = mm.HarmonicBondForce()
        harmonic_force.setName("Static reaction harmonic bond")
        if not self.no_force_groups:
            harmonic_force.setForceGroup(
                EvbForceGroup.REA_HARM_BOND_STATIC.value)

        bond_keys = list(set(self.reactant.bonds) | set(self.product.bonds))
        for key in bond_keys:
            atom_ids = self._key_to_id(key, self.reaction_atoms)
            if atom_ids is None:
                continue
            if key in self.reactant.bonds and key in self.product.bonds:
                bondA = self.reactant.bonds[key]
                fcA = bondA['force_constant']
                eqA = bondA['equilibrium']
                bondB = self.product.bonds[key]
                fcB = bondB['force_constant']
                eqB = bondB['equilibrium']
                eq = eqA * (1 - lam) + eqB * lam
                fc = fcA * (1 - lam) + fcB * lam
                self._add_bond(harmonic_force, atom_ids, eq, fc)
        return harmonic_force

    def _create_dynamic_harmonic_bond_forces(self, lam):

        assert_msg_critical('openmm' in sys.modules,
                            'openmm is required for EvbSystemBuilder.')

        harmonic_force = mm.HarmonicBondForce()
        harmonic_force.setName("Dynamic reaction harmonic bond")
        bond_keys = list(set(self.reactant.bonds) | set(self.product.bonds))
        if not self.no_force_groups:
            harmonic_force.setForceGroup(
                EvbForceGroup.REA_HARM_BOND_DYNAMIC.value)
        for key in bond_keys:
            atom_ids = self._key_to_id(key, self.reaction_atoms)
            fcA = fcB = 0
            eqA = eqB = 1
            if key in self.reactant.bonds and key not in self.product.bonds:
                bondA = self.reactant.bonds[key]
                fcA = bondA['force_constant']
                eqA = bondA['equilibrium']

                s1 = self.product.atoms[key[0]]['sigma']
                s2 = self.product.atoms[key[1]]['sigma']
                eqB = 0.5 * (s1 + s2)
                fcB = fcA * self.bonded_integration_bond_fac
                # if model_broken:
            elif key not in self.reactant.bonds and key in self.product.bonds:
                bondB = self.product.bonds[key]
                fcB = bondB['force_constant']
                eqB = bondB['equilibrium']

                s1 = self.product.atoms[key[0]]['sigma']
                s2 = self.product.atoms[key[1]]['sigma']
                eqA = 0.5 * (s1 + s2)

                fcA = fcB * self.bonded_integration_bond_fac

            p = self.dynamic_bond_tightening
            eq = eqA * (1 - lam) + eqB * lam

            # Apply bond tightening. This causes breaking and forming bonds to be more tightly bonded for intermediate lambda values
            min_eq = min(eqA, eqB)
            dif_eq = abs(eqA - eqB)
            if dif_eq > 0:
                eq = (((eq - min_eq) / dif_eq)**p * dif_eq) + min_eq

            gamma = self.dynamic_bond_fc_factor
            fc_scaling = 4 * (lam - 0.5)**2 * (1 - gamma) + gamma

            fc = fcA * (1 - lam) + fcB * lam
            fc = fc * fc_scaling
            self._add_bond(harmonic_force, atom_ids, eq, fc)

        return harmonic_force

    def _create_morse_force(self, lam, model_broken=False):
        morse_expr = "D*(1-exp(-a*(r-re)))^2;"
        morse_force = mm.CustomBondForce(morse_expr)
        morse_force.setName("Reaction morse bond")
        morse_force.addPerBondParameter("D")
        morse_force.addPerBondParameter("a")
        morse_force.addPerBondParameter("re")
        if not self.no_force_groups:
            morse_force.setForceGroup(EvbForceGroup.REA_MORSE_BOND.value)

        bond_keys = list(set(self.reactant.bonds) | set(self.product.bonds))

        for key in bond_keys:

            breaking = key in self.reactant.bonds and key not in self.product.bonds
            forming = key not in self.reactant.bonds and key in self.product.bonds
            if not (breaking or forming):
                continue

            if breaking:
                bond = self.reactant.bonds[key]
                scaling = 1 - lam
            elif forming:
                bond = self.product.bonds[key]
                scaling = lam

            atom_ids = self._key_to_id(key, self.reaction_atoms)
            if not model_broken:
                eq = bond['equilibrium']
                fc = bond['force_constant'] * scaling
            else:
                if breaking:
                    eqA = bond['equilibrium']
                    fcA = bond['force_constant']
                    s1 = self.product.atoms[key[0]]['sigma']
                    s2 = self.product.atoms[key[1]]['sigma']
                    eqB = 0.5 * (s1 + s2)
                    fcB = fcA

                elif forming:
                    s1 = self.product.atoms[key[0]]['sigma']
                    s2 = self.product.atoms[key[1]]['sigma']
                    eqA = 0.5 * (s1 + s2)
                    eqB = bond['equilibrium']
                    fcB = bond['force_constant']
                    fcA = fcB

                eq = eqA * (1 - lam) + eqB * lam
                fc = fcA * (1 - lam) + fcB * lam
            self.ostream.flush()
            self._add_morse_bond(morse_force, atom_ids, eq, fc, bond)

        self.ostream.flush()
        return morse_force

    def _create_harmonic_angle_forces(self, lam: float, model_broken):

        assert_msg_critical('openmm' in sys.modules,
                            'openmm is required for EvbSystemBuilder.')

        harmonic_force = mm.HarmonicAngleForce()
        harmonic_force.setName("Reaction angle")
        if not self.no_force_groups:
            harmonic_force.setForceGroup(EvbForceGroup.REA_ANGLE.value)

        angle_keys = list(set(self.reactant.angles) | set(self.product.angles))
        for key in angle_keys:
            atom_ids = self._key_to_id(key, self.reaction_atoms)
            if atom_ids is None:
                continue
            if (key in self.reactant.angles.keys()
                    and key in self.product.angles.keys()):

                angleA = self.reactant.angles[key]
                angleB = self.product.angles[key]

                fcA = angleA['force_constant']
                eqA = angleA['equilibrium']

                fcB = angleB['force_constant']
                eqB = angleB['equilibrium']

            elif key in self.reactant.angles.keys():
                # take angle from reactant, and from product structure
                angleA = self.reactant.angles[key]

                fcA = angleA['force_constant']
                eqA = angleA['equilibrium']
                if model_broken:
                    coords = self.product.molecule.get_coordinates_in_angstrom()
                    eqB = self.measure_angle(coords[key[0]],
                                             coords[key[1]],
                                             coords[key[2]],
                                             angle_unit='degree')
                    fcB = fcA * self.bonded_integration_angle_fac
                    self.ostream.print_warning(
                        f"Using product structure to determine angle {key} in reactant, with equilibrium {eqB} and force constant {fcB} for lambda {lam}. Check if this is indeed correct."
                    )
                else:
                    eqB = eqA
                    fcB = 0
            else:
                angleB = self.product.angles[key]
                eqB = angleB['equilibrium']
                fcB = angleB['force_constant']
                fcA = fcB
                if model_broken:
                    coords = self.reactant.molecule.get_coordinates_in_angstrom(
                    )
                    eqA = self.measure_angle(coords[key[0]],
                                             coords[key[1]],
                                             coords[key[2]],
                                             angle_unit='degree')
                    fcA = fcB * self.bonded_integration_angle_fac
                    self.ostream.print_warning(
                        f"Using reactant structure to determine angle {key} in reactant, with equilibrium {eqB} and force constant {fcB} for lambda {lam}. Check if this is indeed correct."
                    )
                else:
                    eqA = eqB
                    fcA = 0
            eq = eqA * (1 - lam) + eqB * lam
            fc = fcA * (1 - lam) + fcB * lam
            self._add_angle(harmonic_force, atom_ids, eq, fc)
        return harmonic_force

    def _create_proper_torsion_forces(self, lam, pes=False):
        """Build proper torsion forces for the reaction system at the given lambda.

        For conformational TS integration systems (self.conformer_active_torsion
        is set and pes=False), the Fourier terms on the active central bond are
        omitted and a periodic restraint  k*(1 - cos(theta - theta0))  is added
        instead, with theta0 shifting from phi_reactant (lam=0) to phi_product
        (lam=1) along the shortest arc.

        Returns:
            list[mm.Force]: always a list — one PeriodicTorsionForce for the
            normal case, plus a CustomTorsionForce appended for conformational
            integration systems.
        """
        assert_msg_critical('openmm' in sys.modules,
                            'openmm is required for EvbSystemBuilder.')

        fourier_force = mm.PeriodicTorsionForce()
        fourier_force.setName("Reaction proper fourier torsion")
        if not self.no_force_groups:
            fourier_force.setForceGroup(EvbForceGroup.REA_TORSION.value)

        # When modelling a conformational TS integration system, identify the
        # central bond to skip so the restraint below takes its place.
        active_torsion = self.conformer_active_torsion
        if active_torsion is not None and not pes:
            active_torsion_central = tuple(
                sorted([active_torsion[1], active_torsion[2]]))
        else:
            active_torsion_central = None

        dihedral_keys = list(
            set(self.reactant.dihedrals) | set(self.product.dihedrals))
        for key in dihedral_keys:
            # Skip Fourier terms on the active central bond for integration
            # systems; the periodic restraint below replaces them.
            if (active_torsion_central is not None and tuple(
                    sorted([key[1], key[2]])) == active_torsion_central):
                continue

            atom_ids = self._key_to_id(key, self.reaction_atoms)
            if atom_ids is None:
                continue
            if (key in self.reactant.dihedrals.keys()
                    and key in self.product.dihedrals.keys()):
                dihedA = self.reactant.dihedrals[key]
                dihedB = self.product.dihedrals[key]
                # if the torsion on one side of the reaction completely disapears, then turn of the torsion much earlier / turn it on later
                if dihedA['barrier'] == 0 or dihedB['barrier'] == 0:
                    reascale, proscale = self._get_lambda_scaling(
                        lam, self.torsion_lambda_switch)
                else:
                    reascale = 1 - lam
                    proscale = lam
                self._add_torsion(fourier_force, dihedA, atom_ids, reascale)
                self._add_torsion(fourier_force, dihedB, atom_ids, proscale)
            else:
                # if the torsion on one side of the reaction completely disapears, then turn of the torsion much earlier / turn it on later
                reascale, proscale = self._get_lambda_scaling(
                    lam, self.torsion_lambda_switch)
                if key in self.reactant.dihedrals.keys():
                    scale = reascale
                    dihed = self.reactant.dihedrals[key]
                else:
                    scale = proscale
                    dihed = self.product.dihedrals[key]
                if scale > 0:
                    self._add_torsion(fourier_force, dihed, atom_ids, scale)

        forces = [fourier_force]

        # For conformational TS integration systems, append the periodic
        # restraint whose minimum shifts from phi_reactant to phi_product.
        if active_torsion is not None and not pes:
            phi_rea_rad = self.conformer_phi_reactant * self.deg_to_rad
            phi_pro_rad = self.conformer_phi_product * self.deg_to_rad
            delta_phi_rad = phi_pro_rad - phi_rea_rad
            if delta_phi_rad > math.pi:
                delta_phi_rad -= 2 * math.pi
            elif delta_phi_rad <= -math.pi:
                delta_phi_rad += 2 * math.pi
            theta0 = phi_rea_rad + lam * delta_phi_rad

            restraint_force = mm.CustomTorsionForce(
                "k_dih * (1 - cos(theta - theta0))")
            restraint_force.setName("Conformer dihedral restraint")
            restraint_force.addGlobalParameter("k_dih", self.conformer_k)
            restraint_force.addPerTorsionParameter("theta0")
            if not self.no_force_groups:
                restraint_force.setForceGroup(EvbForceGroup.REA_TORSION.value)

            active_atom_ids = self._key_to_id(active_torsion,
                                              self.reaction_atoms)
            if active_atom_ids is not None:
                restraint_force.addTorsion(
                    active_atom_ids[0],
                    active_atom_ids[1],
                    active_atom_ids[2],
                    active_atom_ids[3],
                    [theta0],
                )
            forces.append(restraint_force)

        return forces

    # Create a linear switching function that turns the force on only past the lambda-switch
    def _get_lambda_scaling(self, lam, lambda_switch):
        if lambda_switch == 0:
            return 1 - lam, lam
        if lambda_switch == 1:
            return 0, 0

        a = -1 / lambda_switch
        b = 1
        reascale = max(0, a * lam + b)
        a = 1 / (1 - lambda_switch)
        b = 1 - a
        proscale = max(0, a * lam + b)
        return reascale, proscale

    def _create_improper_torsion_forces(self, lam):

        assert_msg_critical('openmm' in sys.modules,
                            'openmm is required for EvbSystemBuilder.')

        fourier_force = mm.PeriodicTorsionForce()
        fourier_force.setName("Reaction improper fourier torsion")
        if not self.no_force_groups:
            fourier_force.setForceGroup(EvbForceGroup.REA_IMP.value)

        dihedral_keys = list(
            set(self.reactant.impropers) | set(self.product.impropers))
        for key in dihedral_keys:
            atom_ids = self._key_to_id(key, self.reaction_atoms)
            if atom_ids is None:
                continue
            if (key in self.reactant.impropers.keys()
                    and key in self.product.impropers.keys()):
                dihedA = self.reactant.impropers[key]
                dihedB = self.product.impropers[key]

                if dihedA['barrier'] == 0 or dihedB['barrier'] == 0:
                    reascale, proscale = self._get_lambda_scaling(
                        lam, self.torsion_lambda_switch)
                else:
                    reascale = 1 - lam
                    proscale = lam

                self._add_torsion(fourier_force,
                                  dihedA,
                                  atom_ids,
                                  reascale,
                                  improper=True)
                self._add_torsion(fourier_force,
                                  dihedB,
                                  atom_ids,
                                  proscale,
                                  improper=True)
            else:
                reascale, proscale = self._get_lambda_scaling(
                    lam, self.torsion_lambda_switch)
                if key in self.reactant.impropers.keys():
                    scale = reascale
                    dihed = self.reactant.impropers[key]
                else:
                    scale = proscale
                    dihed = self.product.impropers[key]
                if scale > 0:
                    self._add_torsion(fourier_force,
                                      dihed,
                                      atom_ids,
                                      scale,
                                      improper=True)
        return fourier_force

    def _add_bond(self, bond_force, atom_id, equil, fc):
        if fc > 0:
            bond_force.addBond(
                atom_id[0],
                atom_id[1],
                equil,
                fc,
            )

    def _add_morse_bond(self, bond_force, atom_id, equil, fc, bond_dict):
        if "D" not in bond_dict.keys():
            D = self.morse_D_default
            if self.verbose:
                self.ostream.print_info(
                    f"No D value associated with bond {atom_id[0]} {atom_id[1]}. Setting to default value {self.morse_D_default}"
                )
        else:
            D = bond_dict["D"]
        a = math.sqrt(fc / (2 * D))

        if fc > 0:
            bond_force.addBond(
                atom_id[0],
                atom_id[1],
                [D, a, equil],
            )

    def _add_angle(self, angle_force, atom_id, equil, fc):
        if fc > 0:
            angle_force.addAngle(
                atom_id[0],
                atom_id[1],
                atom_id[2],
                equil * self.deg_to_rad,
                fc,
            )

    def _add_torsion(self,
                     fourier_force,
                     torsion_dict,
                     atom_id,
                     barrier_scaling=1.,
                     improper=False):
        assert torsion_dict["type"] == "Fourier", "Unknown dihedral type"
        if improper:
            atom_id = [atom_id[1], atom_id[0], atom_id[3], atom_id[2]]
        if barrier_scaling > 0:
            if torsion_dict.get("multiple", False):
                for periodicity, phase, barrier in zip(
                        torsion_dict["periodicity"], torsion_dict["phase"],
                        torsion_dict["barrier"]):
                    if barrier > 0:
                        fourier_force.addTorsion(
                            atom_id[0],
                            atom_id[1],
                            atom_id[2],
                            atom_id[3],
                            abs(periodicity),
                            phase * self.deg_to_rad,
                            barrier_scaling * barrier,
                        )
            else:
                if torsion_dict["barrier"] > 0:
                    fourier_force.addTorsion(
                        atom_id[0],
                        atom_id[1],
                        atom_id[2],
                        atom_id[3],
                        torsion_dict["periodicity"],
                        torsion_dict["phase"] * self.deg_to_rad,
                        barrier_scaling * torsion_dict["barrier"],
                    )

    def _get_couloumb_expression(self, soft_core: bool):
        soft_core_expression = (
            " (1-l) * CoultotA + l  * CoultotB; "
            "CoultotA   = step(r-rqA) * CoulA + step(rqA-r) * CoullinA;"
            "CoultotB   = step(r-rqB) * CoulB + step(rqB-r) * CoullinB;"
            "CoullinA   = k * ( ( qqA / rqAdiv^3 ) * r^2 - 3 * ( qqA / rqAdiv^2 ) * r + 3 * ( qqA / rqAdiv ) );"
            "CoullinB   = k * ( ( qqB / rqBdiv^3 ) * r^2 - 3 * ( qqB / rqBdiv^2 ) * r + 3 * ( qqB / rqBdiv ) );"
            "rqAdiv     = select(rqA,rqA,1);"
            "rqBdiv     = select(rqB,rqB,1);"
            "rqA        = alphaq * ( 1 + sigmaq * qqA );"
            "rqB        = alphaq * ( 1 + sigmaq * qqB );"
            "CoulA      = k*qqA/r; "
            "CoulB      = k*qqB/r; "
            ""
            f"k         = {self.k};"
            ""
            f"alphalj   = {self.sc_alpha_lj};"
            f"alphaq    = {self.sc_alpha_q};"
            f"sigmaq    = {self.sc_sigma_q};"
            f"pow       = {self.sc_power};")

        hard_core_expression = (" (1-l) * k*qqA/r + l * k*qqB/r; "
                                ""
                                f"k         = {self.k};")

        if soft_core:
            return soft_core_expression
        else:
            return hard_core_expression

    def _get_lj_expression(self, soft_core: bool):
        soft_core_expression = (
            " (1-l) * LjtotA + l * LjtotB; "
            ""
            "LjtotA     = (step(r - rljA) * LjA + step(rljA - r) * LjlinA);"
            "LjtotB     = (step(r - rljB) * LjB + step(rljB - r) * LjlinB);"
            "LjlinA     = ( (78*A12) / (rljAdiv^14) - (21*A6) / (rljAdiv^8) )*r^2 - ( (168*A12) / (rljAdiv^13) - (48*A6) / (rljAdiv^7) )*r + ( (91*A12) / (rljAdiv^12) - (28*A6) / (rljAdiv^6) );"
            "LjlinB     = ( (78*B12) / (rljBdiv^14) - (21*B6) / (rljBdiv^8) )*r^2 - ( (168*B12) / (rljBdiv^13) - (48*B6) / (rljBdiv^7) )*r + ( (91*B12) / (rljBdiv^12) - (28*B6) / (rljBdiv^6) );"
            # if rljA = 0, returns 1, otherwise returns rljA. Prevents division by 0 while the step factor is already 0
            "rljAdiv    = select(rljA,rljA,1);"
            "rljBdiv    = select(rljB,rljB,1);"
            "rljA       = alphalj * ( (26/7 ) * A6 ) ^ pow;"
            "rljB       = alphalj * ( (26/7 ) * B6 ) ^ pow;"
            "LjA        = A12 / r^12 - A6 / r^6;"
            "LjB        = B12 / r^12 - B6 / r^6;"
            "A12        = 4 * epsilonA * sigmaA ^ 12; "
            "B12        = 4 * epsilonB * sigmaB ^ 12; "
            "A6         = 4 * epsilonA * sigmaA ^ 6; "
            "B6         = 4 * epsilonB * sigmaB ^ 6; "
            ""
            f"alphalj   = {self.sc_alpha_lj};"
            f"alphaq    = {self.sc_alpha_q};"
            f"sigmaq    = {self.sc_sigma_q};"
            f"pow       = {self.sc_power};")

        hard_core_expression = (" (1-l) * LjtotA "
                                "  + l  * LjtotB; "
                                ""
                                "LjtotA     = A12 / r^12 - A6 / r^6;"
                                "LjtotB     = B12 / r^12 - B6 / r^6;"
                                "A12        = 4 * epsilonA * sigmaA ^ 12; "
                                "B12        = 4 * epsilonB * sigmaB ^ 12; "
                                "A6         = 4 * epsilonA * sigmaA ^ 6; "
                                "B6         = 4 * epsilonB * sigmaB ^ 6; ")

        if soft_core:
            return soft_core_expression
        else:
            return hard_core_expression

    # Merge exceptions takes the union of the reactant and product bonds, and bases all exceptions on that.
    # This causes all exceptdions to be constant. This will include constant 1-4 interactions between atoms that bonded in either the reactant or the product, but not in the other.

    # Broken exceptions will add exceptions for bonds that are changing between the reactant and product.
    # This is weaker then merge_exceptions, as it will not add 1-4 interactions between atoms that are not bonded, even if they are bonded on the other side of the reaction.
    # This is used to remove the nonbonded interactions that are modelled by other forces

    # Only changing bonds will only add nonbonded interactions for atom pairs that are part of changing bonds.
    # Exclude changing bonds will do the opposite

    # Turning both off gives the proper reactant or product for lam=0 or 1
    def _create_nonbonded_forces(
        self,
        lam,
        merge_exceptions=False,
        broken_exceptions=False,
        lj_soft_core=False,
        coul_soft_core=False,
        only_changing_bonds=False,
        exclude_changing_bonds=False,
    ):

        assert_msg_critical('openmm' in sys.modules,
                            'openmm is required for EvbSystemBuilder.')

        coulomb_force = mm.CustomBondForce(
            self._get_couloumb_expression(coul_soft_core))
        name = "Reaction internal coulomb"
        if coul_soft_core:
            name += " soft-core"

        coulomb_force.setName(name)
        coulomb_force.addPerBondParameter("qqA")
        coulomb_force.addPerBondParameter("qqB")
        coulomb_force.addGlobalParameter("l", lam)

        lj_force = mm.CustomBondForce(self._get_lj_expression(lj_soft_core))
        name = "Reaction internal lennard-jones"
        if lj_soft_core:
            name += " soft-core"
        lj_force.setName(name)
        lj_force.addPerBondParameter("sigmaA")
        lj_force.addPerBondParameter("sigmaB")
        lj_force.addPerBondParameter("epsilonA")
        lj_force.addPerBondParameter("epsilonB")
        lj_force.addGlobalParameter("l", lam)

        # The bonded parameter forces the same exceptions across the entire reaction
        if merge_exceptions:
            bonds = list(set(self.reactant.bonds) | set(self.product.bonds))
            reactant_bonds = bonds
            product_bonds = bonds
        else:
            reactant_bonds = self.reactant.bonds
            product_bonds = self.product.bonds
        reactant_exceptions = self._create_exceptions_from_bonds(
            self.reactant.atoms, reactant_bonds)
        product_exceptions = self._create_exceptions_from_bonds(
            self.product.atoms, product_bonds)

        changing_bonds = list(
            (set(self.reactant.bonds) ^ set(self.product.bonds)))

        atom_keys = self.reactant.atoms.keys()

        # Loop over all atoms, and check if their id's are part of any exceptions
        for i in atom_keys:
            for j in atom_keys:
                if i < j:
                    key = (i, j)
                    atom_ids = self._key_to_id(key, self.reaction_atoms)
                    if atom_ids is None:
                        continue

                    if only_changing_bonds and key not in changing_bonds:
                        continue
                    if exclude_changing_bonds and key in changing_bonds:
                        continue

                    if broken_exceptions:
                        if key in changing_bonds:
                            continue
                    # Remove any exception from the nonbondedforce
                    # and add it instead to the exception bond force
                    if key in reactant_exceptions.keys():
                        epsilonA = reactant_exceptions[key]["epsilon"]
                        sigmaA = reactant_exceptions[key]["sigma"]
                        qqA = reactant_exceptions[key]["qq"]

                    else:
                        atomA1 = self.reactant.atoms[key[0]]
                        atomA2 = self.reactant.atoms[key[1]]
                        epsilonA = math.sqrt(atomA1["epsilon"] *
                                             atomA2["epsilon"])
                        sigmaA = 0.5 * (atomA1["sigma"] + atomA2["sigma"])
                        qqA = atomA1["charge"] * atomA2["charge"]

                    if key in product_exceptions.keys():
                        epsilonB = product_exceptions[key]["epsilon"]
                        sigmaB = product_exceptions[key]["sigma"]
                        qqB = product_exceptions[key]["qq"]
                    else:
                        atomB1 = self.product.atoms[key[0]]
                        atomB2 = self.product.atoms[key[1]]
                        epsilonB = math.sqrt(atomB1["epsilon"] *
                                             atomB2["epsilon"])
                        sigmaB = 0.5 * (atomB1["sigma"] + atomB2["sigma"])
                        qqB = atomB1["charge"] * atomB2["charge"]

                    if sigmaA == 1.0:
                        sigmaA = sigmaB
                    elif sigmaB == 1.0:
                        sigmaB = sigmaA

                    if not (qqA == 0.0 and qqB == 0.0):
                        coulomb_force.addBond(
                            atom_ids[0],
                            atom_ids[1],
                            [qqA, qqB],
                        )
                    if not (epsilonA == 0.0 and epsilonB == 0.0):
                        lj_force.addBond(
                            atom_ids[0],
                            atom_ids[1],
                            [sigmaA, sigmaB, epsilonA, epsilonB],
                        )

        return lj_force, coulomb_force

    def _create_exceptions_from_bonds(
            self, particles, bonds) -> dict[tuple[int, int], dict[str, float]]:

        exclusions = [set() for _ in range(len(particles))]
        bonded12 = [set() for _ in range(len(particles))]
        exceptions = {}

        # Populate bonded12 with bonds.
        for bond in bonds:
            bonded12[bond[0]].add(bond[1])
            bonded12[bond[1]].add(bond[0])

        # Find particles separated by 1, 2, or 3 bonds.
        for i in range(len(particles)):
            self._add_exclusions_to_set(bonded12, exclusions[i], i, i, 2)

        # Find particles separated by 1 or 2 bonds and create exceptions.
        for i in range(len(exclusions)):
            bonded13 = set()
            self._add_exclusions_to_set(bonded12, bonded13, i, i, 1)
            for j in exclusions[i]:
                if j < i:
                    if j not in bonded13:
                        # This is a 1-4 interaction.
                        particle1 = particles[j]
                        particle2 = particles[i]
                        charge_prod = (self.coul14_scale * particle1["charge"] *
                                       particle2["charge"])
                        sigma = 0.5 * (particle1["sigma"] + particle2["sigma"])
                        epsilon = (
                            self.lj14_scale *
                            (particle1["epsilon"] * particle2["epsilon"])**0.5)
                        exceptions[tuple(sorted((i, j)))] = {
                            "qq": charge_prod,
                            "sigma": sigma,
                            "epsilon": epsilon,
                        }
                    else:
                        # This interaction should be completely excluded.
                        exceptions[tuple(sorted((i, j)))] = {
                            "qq": 0.0,
                            "sigma": 1.0,
                            "epsilon": 0.0,
                        }
        return exceptions

    @staticmethod
    def _add_exclusions_to_set(bonded12, exclusions, base_particle,
                               from_particle, current_level):
        for i in bonded12[from_particle]:
            if i != base_particle:
                exclusions.add(i)
            if current_level > 0:
                EvbSystemBuilder._add_exclusions_to_set(bonded12, exclusions,
                                                        base_particle, i,
                                                        current_level - 1)

    def _create_constraint_forces(self, lam):

        assert_msg_critical('openmm' in sys.modules,
                            'openmm is required for EvbSystemBuilder.')

        bond_constraint = mm.HarmonicBondForce()
        bond_constraint.setName("Bond constraint")

        constant_force = mm.CustomBondForce("-k*r")
        constant_force.setName("Constant bond force constraint")
        constant_force.addPerBondParameter("k")

        angle_constraint = mm.HarmonicAngleForce()
        angle_constraint.setName("Angle constraint")
        torsion_constraint = mm.CustomTorsionForce("0.5*k*(theta-theta0)^2")
        torsion_constraint.setName("Harmonic torsion constraint")
        torsion_constraint.addPerTorsionParameter("theta0")
        torsion_constraint.addPerTorsionParameter("k")
        if not self.no_force_groups:
            bond_constraint.setForceGroup(EvbForceGroup.CONSTRAINT.value)
            constant_force.setForceGroup(EvbForceGroup.CONSTRAINT.value)
            angle_constraint.setForceGroup(EvbForceGroup.CONSTRAINT.value)
            torsion_constraint.setForceGroup(EvbForceGroup.CONSTRAINT.value)

        if len(self.constraints) > 0:
            if self.verbose:
                self.ostream.print_info(
                    f"Adding constraints: {self.constraints}")

        for constraint in self.constraints:
            key = list(constraint.keys())[0]
            if "lambda_couple" in constraint[key].keys():
                if constraint[key]["lambda_couple"] == 1:
                    scale = lam
                elif constraint[key]["lambda_couple"] == -1:
                    scale = 1 - lam
                else:
                    scale = 1
            else:
                scale = 1

            atomids = self._key_to_id(key, self.reaction_atoms)
            if len(key) == 2:
                if constraint[key]['type'] == 'harmonic':
                    self._add_bond(
                        bond_constraint,
                        atomids,
                        constraint[key]["equilibrium"],
                        constraint[key]["force_constant"] * scale,
                    )
                elif constraint[key]['type'] == 'linear':
                    constant_force.addBond(
                        atomids[0],
                        atomids[1],
                        [constraint[key]["force_constant"] * scale],
                    )
                else:
                    raise ValueError(
                        f"Unknown constraint bond type {constraint[key]['type']}"
                    )
            if len(key) == 3:
                self._add_angle(
                    angle_constraint,
                    atomids,
                    constraint[key]["equilibrium"],
                    constraint[key]["force_constant"] * scale,
                )

            if len(key) == 4:
                torsion_constraint.addTorsion(
                    atomids[0],
                    atomids[1],
                    atomids[2],
                    atomids[3],
                    [
                        constraint[key]["equilibrium"] * self.deg_to_rad,
                        constraint[key]["force_constant"] * scale,
                    ],
                )
        return bond_constraint, constant_force, angle_constraint, torsion_constraint

    def _key_to_id(self, key: tuple[int, ...],
                   atoms: dict | list) -> list[int] | None:
        for i in key:
            if isinstance(atoms, dict):
                if i not in atoms.keys():
                    return None
            elif isinstance(atoms, list):
                if i >= len(atoms):
                    return None

        return [atoms[key[i]].index for i in range(len(key))]

    @staticmethod
    def measure_length(v1, v2):
        """
        Calculates the distance between v1 and v2
        """

        return np.linalg.norm(np.array(v1) - np.array(v2))

    @staticmethod
    def measure_angle(v1, v2, v3, angle_unit='radian'):
        """
        Calculates the angle between v1 and v2 and v3
        """

        vector_u = (v1 - v2) / np.linalg.norm(v1 - v2)
        vector_v = (v3 - v2) / np.linalg.norm(v3 - v2)

        # Calculate the angle using the dot product
        dot_product = np.dot(vector_u, vector_v)
        angle = np.arccos(dot_product)
        angle = angle * 180 / np.pi  # Convert to degrees

        return angle

    @staticmethod
    def measure_dihedral(v1, v2, v3, v4, angle_unit='radian'):
        """
        Calculates the dihedral angle between v1, v2, v3 and v4
        """

        v21 = v2 - v1
        v32 = v3 - v2
        v43 = v4 - v3

        u21 = v21 / np.linalg.norm(v21)
        u32 = v32 / np.linalg.norm(v32)
        u43 = v43 / np.linalg.norm(v43)

        cos_theta_123 = -np.vdot(u21, u32)
        cos_theta_234 = -np.vdot(u32, u43)

        sin_theta_123 = math.sqrt(1.0 - cos_theta_123**2)
        sin_theta_234 = math.sqrt(1.0 - cos_theta_234**2)

        cos_phi = ((cos_theta_123 * cos_theta_234 - np.vdot(u21, u43)) /
                   (sin_theta_123 * sin_theta_234))
        sin_phi = -(np.vdot(u43, np.cross(u21, u32)) /
                    (sin_theta_123 * sin_theta_234))

        phi_in_radian = safe_arccos(cos_phi)
        if sin_phi < 0.0:
            phi_in_radian *= -1.0

        assert_msg_critical(angle_unit.lower() in ['degree', 'radian'],
                            'Molecule.get_dihedral: Invalid angle unit')

        if angle_unit.lower() == 'degree':
            return 180.0 * phi_in_radian / math.pi
        else:
            return phi_in_radian

    def save_systems_as_xml(self, systems: dict, folder: str):
        """Save the systems as xml files to the given folder.

        Args:
            systems (dict): The systems to save
            folder (str): The folder relative to the current working directory to save the systems to.
        """

        assert_msg_critical('openmm' in sys.modules,
                            'openmm is required for EvbDriver.')

        path = Path().cwd() / folder
        self.ostream.print_info(f"Saving systems to {path}")
        self.ostream.flush()
        if not path.exists():
            path.mkdir(parents=True, exist_ok=True)
        for name, system in systems.items():
            if isinstance(name, float) or isinstance(name, int):
                filename = f"{name:.3f}_sys.xml"
            else:
                filename = f"{name}_sys.xml"
            with open(path / filename, mode="w", encoding="utf-8") as output:
                output.write(mm.XmlSerializer.serialize(system))

    def load_systems_from_xml(self, folder: str):
        """Load the systems from xml files in the given folder.

        Args:
            folder (str): The folder relative to the current working directory to load the systems from.
        Returns:
            dict: The loaded systems
        """

        assert_msg_critical('openmm' in sys.modules,
                            'openmm is required for EvbDriver.')

        systems = {}
        path = Path().cwd() / folder
        for lam in self.Lambda:
            with open(path / f"{lam:.3f}_sys.xml", mode="r",
                      encoding="utf-8") as input:
                systems[lam] = mm.XmlSerializer.deserialize(input.read())
        return systems


# The EVB procedure uses two different potentials, one for the integration of the EOMs to explore phase space, and another one for the calculation of the PES
# The integration potential is optimized to explore all relevant areas of phase space in an efficient manner. This includes constraints and distance restraints.
# The PES potential is optimized to calculate the potential energy surface accurate. Constraints and restraints are ommitted for this.
# Soft core long range potentials provide faster convergence for averages and are thus included in the PES potential. They cause unstable integration though.
class EvbForceGroup(Enum):
    DEFAULT = auto(
    )  # Default force group, included for both integration and energy calculations
    # THERMOSTAT = auto()  # Thermostat
    SYSLJ_STATIC = auto()  # Integration lennard-jones potential
    SYSCOUL_STATIC = auto()  # Integration coulombic potential
    SYSLJ_DYNAMIC = auto()  # Dynamic lennard-jones potential
    SYSCOUL_DYNAMIC = auto()  # Dynamic coulombic potential
    CONSTRAINT = auto()

    CMM_REMOVER = auto()  # Center of mass motion remover
    NB_FORCE_INT = auto()  # Solvent-solvent and solvent-solute nb force
    BAROSTAT = auto()  # Barostat
    E_FIELD = auto()  # Electric field force
    REA_HARM_BOND_STATIC = auto()  # Bonded forces for the reaction atoms
    REA_HARM_BOND_DYNAMIC = auto(
    )  # Static bonded forces for the reaction atoms
    REA_MORSE_BOND = auto()
    REA_ANGLE = auto()
    REA_TORSION = auto()
    REA_IMP = auto()

    # Constraints that also should be included in the PES calculations. Currently only used for the linear bond constraint
    PES_CONSTRAINT = auto()
    RESTRAINT = auto()
    SOLVENT = auto(
    )  # All solvent-solvent interactions. Does not include the solute-solvent long range interaction
    CARBON = auto()  # Graphene and CNTs
    PDB = auto()  # Bonded forces added from the PDB
    SOL_COUL = auto()
    SOL_LJ = auto()
    POSRES = auto()

    @classmethod
    def pes_forcegroups(cls):
        max_ind = cls.POSRES.value
        return set(range(1, max_ind))

    @classmethod
    def get_header(cls):
        # Simple method for printing a descrpitive header to be used in force group logging files
        header = ""
        # pes_forcegroups = cls.pes_force_groups()
        for fg in cls:
            header += f"{fg.name}({fg.value}), "
        header = header[:-2] + '\n'
        return header
