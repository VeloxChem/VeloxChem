#
# This file contains derivative work of the OpenMM molecular simulation toolkit
# and is licensed under the MIT License. Portions of this file are adapted from
# openmm.app.internal.customgbforces.
#
# SPDX-License-Identifier: MIT
#
# Portions copyright (c) 2012-2022 University of Virginia and the Authors.
# Authors: Christoph Klein, Michael R. Shirts
# Contributors: Jason M. Swails, Peter Eastman, Justin L. MacCallum
#
# Modifications made by VeloxChem developers, 2026.
#
# Permission is hereby granted, free of charge, to any person obtaining a
# copy of this software and associated documentation files (the "Software"),
# to deal in the Software without restriction, including without limitation
# the rights to use, copy, modify, merge, publish, distribute, sublicense,
# and/or sell copies of the Software, and to permit persons to whom the
# Software is furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
# THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
# DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
# OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
# USE OR OTHER DEALINGS IN THE SOFTWARE.

import math
import sys
from collections import defaultdict

try:
    import openmm as mm
except ImportError:
    pass


class GBImplicitSolvent:
    """Builds a CustomGBForce for any of the five Amber GB variants.

    Usage::

        gb_force = GBImplicitSolvent.build(
            model='obc2',
            topology=topology,
            nb_force=nb_force,
            sv_dielectric=78.39,
            sl_dielectric=1.0,
        )
        system.addForce(gb_force)
    """

    @staticmethod
    def get_valid_models():
        return ('hct', 'obc1', 'obc2', 'gbn', 'gbn2')

    @staticmethod
    def _get_d0():
        return [
            d / 10.0 for d in [
                2.26685, 2.32548, 2.38397, 2.44235, 2.50057, 2.55867, 2.61663,
                2.67444, 2.73212, 2.78965, 2.84705, 2.90430, 2.96141, 3.01840,
                3.07524, 3.13196, 3.18854, 3.24498, 3.30132, 3.35752, 3.41360,
                2.31191, 2.37017, 2.42830, 2.48632, 2.54420, 2.60197, 2.65961,
                2.71711, 2.77449, 2.83175, 2.88887, 2.94586, 3.00273, 3.05948,
                3.11610, 3.17260, 3.22897, 3.28522, 3.34136, 3.39738, 3.45072,
                2.35759, 2.41549, 2.47329, 2.53097, 2.58854, 2.64600, 2.70333,
                2.76056, 2.81766, 2.87465, 2.93152, 2.98827, 3.04490, 3.10142,
                3.15782, 3.21411, 3.27028, 3.32634, 3.38230, 3.43813, 3.49387,
                2.40380, 2.46138, 2.51885, 2.57623, 2.63351, 2.69067, 2.74773,
                2.80469, 2.86152, 2.91826, 2.97489, 3.03140, 3.08781, 3.14410,
                3.20031, 3.25638, 3.31237, 3.36825, 3.42402, 3.47970, 3.53527,
                2.45045, 2.50773, 2.56492, 2.62201, 2.67900, 2.73590, 2.79270,
                2.84940, 2.90599, 2.96250, 3.01890, 3.07518, 3.13138, 3.18748,
                3.24347, 3.29937, 3.35515, 3.41085, 3.46646, 3.52196, 3.57738,
                2.49750, 2.55450, 2.61143, 2.66825, 2.72499, 2.78163, 2.83818,
                2.89464, 2.95101, 3.00729, 3.06346, 3.11954, 3.17554, 3.23143,
                3.28723, 3.34294, 3.39856, 3.45409, 3.50952, 3.56488, 3.62014,
                2.54489, 2.60164, 2.65830, 2.71488, 2.77134, 2.82780, 2.88412,
                2.94034, 2.99650, 3.05256, 3.10853, 3.16442, 3.22021, 3.27592,
                3.33154, 3.38707, 3.44253, 3.49789, 3.55316, 3.60836, 3.66348,
                2.59259, 2.64910, 2.70553, 2.76188, 2.81815, 2.87434, 2.93044,
                2.98646, 3.04241, 3.09827, 3.15404, 3.20974, 3.26536, 3.32089,
                3.37633, 3.43170, 3.48699, 3.54219, 3.59731, 3.65237, 3.70734,
                2.64054, 2.69684, 2.75305, 2.80918, 2.86523, 2.92122, 2.97712,
                3.03295, 3.08870, 3.14437, 3.19996, 3.25548, 3.31091, 3.36627,
                3.42156, 3.47677, 3.53190, 3.58695, 3.64193, 3.69684, 3.75167,
                2.68873, 2.74482, 2.80083, 2.85676, 2.91262, 2.96841, 3.02412,
                3.07976, 3.13533, 3.19082, 3.24623, 3.30157, 3.35685, 3.41205,
                3.46718, 3.52223, 3.57721, 3.63213, 3.68696, 3.74174, 3.79644,
                2.73713, 2.79302, 2.84884, 2.90459, 2.96027, 3.01587, 3.07140,
                3.12686, 3.18225, 3.23757, 3.29282, 3.34801, 3.40313, 3.45815,
                3.51315, 3.56805, 3.62290, 3.67767, 3.73237, 3.78701, 3.84159,
                2.78572, 2.84143, 2.89707, 2.95264, 3.00813, 3.06356, 3.11892,
                3.17422, 3.22946, 3.28462, 3.33971, 3.39474, 3.44971, 3.50460,
                3.55944, 3.61421, 3.66891, 3.72356, 3.77814, 3.83264, 3.88710,
                2.83446, 2.89000, 2.94547, 3.00088, 3.05621, 3.11147, 3.16669,
                3.22183, 3.27689, 3.33191, 3.38685, 3.44174, 3.49656, 3.55132,
                3.60602, 3.66066, 3.71523, 3.76975, 3.82421, 3.87860, 3.93293,
                2.88335, 2.93873, 2.99404, 3.04929, 3.10447, 3.15959, 3.21464,
                3.26963, 3.32456, 3.37943, 3.43424, 3.48898, 3.54366, 3.59830,
                3.65287, 3.70737, 3.76183, 3.81622, 3.87056, 3.92484, 3.97905,
                2.93234, 2.98760, 3.04277, 3.09786, 3.15291, 3.20787, 3.26278,
                3.31764, 3.37242, 3.42716, 3.48184, 3.53662, 3.59100, 3.64551,
                3.69995, 3.75435, 3.80867, 3.86295, 3.91718, 3.97134, 4.02545,
                2.98151, 3.03660, 3.09163, 3.14659, 3.20149, 3.25632, 3.31110,
                3.36581, 3.42047, 3.47507, 3.52963, 3.58411, 3.63855, 3.69293,
                3.74725, 3.80153, 3.85575, 3.90991, 3.96403, 4.01809, 4.07211,
                3.03074, 3.08571, 3.14061, 3.19543, 3.25021, 3.30491, 3.35956,
                3.41415, 3.46869, 3.52317, 3.57759, 3.63196, 3.68628, 3.74054,
                3.79476, 3.84893, 3.90303, 3.95709, 4.01111, 4.06506, 4.11897,
                3.08008, 3.13492, 3.18970, 3.24440, 3.29905, 3.35363, 3.40815,
                3.46263, 3.51704, 3.57141, 3.62572, 3.67998, 3.73418, 3.78834,
                3.84244, 3.89650, 3.95051, 4.00447, 4.05837, 4.11224, 4.16605,
                3.12949, 3.18422, 3.23888, 3.29347, 3.34800, 3.40247, 3.45688,
                3.51124, 3.56554, 3.61980, 3.67400, 3.72815, 3.78225, 3.83629,
                3.89030, 3.94425, 3.99816, 4.05203, 4.10583, 4.15961, 4.21333,
                3.17899, 3.23361, 3.28815, 3.34264, 3.39706, 3.45142, 3.50571,
                3.55997, 3.61416, 3.66831, 3.72241, 3.77645, 3.83046, 3.88440,
                3.93831, 3.99216, 4.04598, 4.09974, 4.15347, 4.20715, 4.26078,
                3.22855, 3.28307, 3.33751, 3.39188, 3.44620, 3.50046, 3.55466,
                3.60880, 3.66290, 3.71694, 3.77095, 3.82489, 3.87880, 3.93265,
                3.98646, 4.04022, 4.09395, 4.14762, 4.20126, 4.25485, 4.30840
            ]
        ]

    @staticmethod
    def _get_m0():
        return [
            m * 10.0 for m in [
                0.0381511, 0.0338587, 0.0301776, 0.0270030, 0.0242506,
                0.0218529, 0.0197547, 0.0179109, 0.0162844, 0.0148442,
                0.0135647, 0.0124243, 0.0114047, 0.0104906, 0.00966876,
                0.008928, 0.0082587, 0.00765255, 0.00710237, 0.00660196,
                0.00614589, 0.0396198, 0.0351837, 0.0313767, 0.0280911,
                0.0252409, 0.0227563, 0.0205808, 0.0186681, 0.0169799,
                0.0154843, 0.014155, 0.0129696, 0.0119094, 0.0109584, 0.0101031,
                0.00933189, 0.0086348, 0.00800326, 0.00742986, 0.00690814,
                0.00643255, 0.041048, 0.0364738, 0.0325456, 0.0291532,
                0.0262084, 0.0236399, 0.0213897, 0.0194102, 0.0176622,
                0.0161129, 0.0147351, 0.0135059, 0.0124061, 0.0114192,
                0.0105312, 0.00973027, 0.00900602, 0.00834965, 0.0077535,
                0.00721091, 0.00671609, 0.0424365, 0.0377295, 0.0336846,
                0.0301893, 0.0271533, 0.0245038, 0.0221813, 0.0201371, 0.018331,
                0.0167295, 0.0153047, 0.014033, 0.0128946, 0.0118727, 0.0109529,
                0.0101229, 0.00937212, 0.00869147, 0.00807306, 0.00751003,
                0.00699641, 0.0437861, 0.0389516, 0.0347944, 0.0311998,
                0.0280758, 0.0253479, 0.0229555, 0.0208487, 0.0189864,
                0.0173343, 0.0158637, 0.0145507, 0.0133748, 0.0123188,
                0.0113679, 0.0105096, 0.0097329, 0.00902853, 0.00838835,
                0.00780533, 0.0072733, 0.0450979, 0.0401406, 0.0358753,
                0.0321851, 0.0289761, 0.0261726, 0.0237125, 0.0215451,
                0.0196282, 0.017927, 0.0164121, 0.0150588, 0.0138465, 0.0127573,
                0.0117761, 0.0108902, 0.0100882, 0.00936068, 0.00869923,
                0.00809665, 0.00754661, 0.0463729, 0.0412976, 0.0369281,
                0.0331456, 0.0298547, 0.026978, 0.0244525, 0.0222264, 0.0202567,
                0.0185078, 0.0169498, 0.0155575, 0.0143096, 0.0131881,
                0.0121775, 0.0112646, 0.010438, 0.00968781, 0.00900559,
                0.00838388, 0.00781622, 0.0476123, 0.0424233, 0.0379534,
                0.034082, 0.0307118, 0.0277645, 0.0251757, 0.0228927, 0.0208718,
                0.0190767, 0.0174768, 0.0160466, 0.0147642, 0.0136112,
                0.0125719, 0.0116328, 0.0107821, 0.0100099, 0.00930735,
                0.00866695, 0.00808206, 0.0488171, 0.0435186, 0.038952,
                0.0349947, 0.0315481, 0.0285324, 0.0258824, 0.0235443,
                0.0214738, 0.0196339, 0.0179934, 0.0165262, 0.0152103,
                0.0140267, 0.0129595, 0.0119947, 0.0111206, 0.0103268,
                0.00960445, 0.00894579, 0.00834405, 0.0499883, 0.0445845,
                0.0399246, 0.0358844, 0.032364, 0.0292822, 0.0265729, 0.0241815,
                0.0220629, 0.0201794, 0.0184994, 0.0169964, 0.0156479,
                0.0144345, 0.0133401, 0.0123504, 0.0114534, 0.0106386,
                0.00989687, 0.00922037, 0.00860216, 0.0511272, 0.0456219,
                0.040872, 0.0367518, 0.0331599, 0.0300142, 0.0272475, 0.0248045,
                0.0226392, 0.0207135, 0.0189952, 0.0174574, 0.0160771,
                0.0148348, 0.0137138, 0.0126998, 0.0117805, 0.0109452,
                0.0101846, 0.00949067, 0.00885636, 0.0522348, 0.0466315,
                0.0417948, 0.0375973, 0.0339365, 0.030729, 0.0279067, 0.0254136,
                0.023203, 0.0212363, 0.0194809, 0.0179092, 0.016498, 0.0152275,
                0.0140807, 0.013043, 0.012102, 0.0112466, 0.0104676, 0.00975668,
                0.00910664, 0.0533123, 0.0476145, 0.042694, 0.0384218,
                0.0346942, 0.0314268, 0.0285507, 0.026009, 0.0237547, 0.0217482,
                0.0199566, 0.018352, 0.0169108, 0.0156128, 0.0144408, 0.0133801,
                0.0124179, 0.011543, 0.010746, 0.0100184, 0.00935302, 0.0543606,
                0.0485716, 0.04357, 0.0392257, 0.0354335, 0.0321082, 0.02918,
                0.0265913, 0.0242943, 0.0222492, 0.0204225, 0.0187859,
                0.0173155, 0.0159908, 0.0147943, 0.0137111, 0.0127282,
                0.0118343, 0.0110197, 0.0102759, 0.00959549, 0.0553807,
                0.0495037, 0.0444239, 0.0400097, 0.0361551, 0.0327736,
                0.0297949, 0.0271605, 0.0248222, 0.0227396, 0.0208788,
                0.0192111, 0.0177122, 0.0163615, 0.0151413, 0.0140361, 0.013033,
                0.0121206, 0.0112888, 0.0105292, 0.00983409, 0.0563738,
                0.0504116, 0.0452562, 0.0407745, 0.0368593, 0.0334235,
                0.0303958, 0.0277171, 0.0253387, 0.0232197, 0.0213257,
                0.0196277, 0.0181013, 0.0167252, 0.0154817, 0.0143552,
                0.0133325, 0.0124019, 0.0115534, 0.0107783, 0.0100688,
                0.0573406, 0.0512963, 0.0460676, 0.0415206, 0.0375468,
                0.0340583, 0.030983, 0.0282614, 0.0258441, 0.0236896, 0.0217634,
                0.020036, 0.0184826, 0.017082, 0.0158158, 0.0146685, 0.0136266,
                0.0126783, 0.0118135, 0.0110232, 0.0102998, 0.0582822,
                0.0521584, 0.0468589, 0.0422486, 0.038218, 0.0346784, 0.0315571,
                0.0287938, 0.0263386, 0.0241497, 0.0221922, 0.0204362,
                0.0188566, 0.0174319, 0.0161437, 0.0149761, 0.0139154,
                0.0129499, 0.0120691, 0.0112641, 0.0105269, 0.0591994,
                0.0529987, 0.0476307, 0.042959, 0.0388734, 0.0352843, 0.0321182,
                0.0293144, 0.0268225, 0.0246002, 0.0226121, 0.0208283,
                0.0192232, 0.0177751, 0.0164654, 0.015278, 0.0141991, 0.0132167,
                0.0123204, 0.0115009, 0.0107504, 0.0600932, 0.053818, 0.0483836,
                0.0436525, 0.0395136, 0.0358764, 0.0326669, 0.0298237,
                0.0272961, 0.0250413, 0.0230236, 0.0212126, 0.0195826,
                0.0181118, 0.0167811, 0.0155744, 0.0144778, 0.0134789,
                0.0125673, 0.0117338, 0.0109702, 0.0609642, 0.0546169,
                0.0491183, 0.0443295, 0.0401388, 0.036455, 0.0332033, 0.030322,
                0.0277596, 0.0254732, 0.0234266, 0.0215892, 0.0199351, 0.018442,
                0.0170909, 0.0158654, 0.0147514, 0.0137365, 0.0128101,
                0.0119627, 0.0111863
            ]
        ]

    @staticmethod
    def _get_screen():
        return {
            'H': (0.85, 1.09085413633, 1.425952, 1.696538),
            'D': (0.85, 1.09085413633, 1.425952, 1.696538),
            'C': (0.72, 0.48435382330, 1.058554, 1.268902),
            'N': (0.79, 0.700147318409, 0.733599, 1.4259728),
            'O': (0.85, 1.06557401132, 1.061039, 0.1840098),
            'F': (0.88, 0.5, 0.5, 0.5),
            'P': (0.86, 0.5, 0.5, 1.5450597),
            'S': (0.96, 0.602256336067, -0.703469, 0.05),
        }

    @staticmethod
    def _get_screen_default():
        return (0.8, 0.5, 0.5, 0.5)

    @staticmethod
    def _get_nucleic_residues():
        return frozenset(['A', 'C', 'G', 'U', 'DA', 'DC', 'DG', 'DT'])

    @staticmethod
    def _get_gbn2_params():
        return {
            'H': [0.788440, 0.798699, 0.437334],
            'D': [0.788440, 0.798699, 0.437334],
            'C': [0.733756, 0.506378, 0.205844],
            'N': [0.503364, 0.316828, 0.192915],
            'O': [0.867814, 0.876635, 0.387882],
            'S': [0.867814, 0.876635, 0.387882],
        }

    @staticmethod
    def _get_gbn2_params_nucleic():
        return {
            'H': [0.537050, 0.362861, 0.116704],
            'D': [0.537050, 0.362861, 0.116704],
            'C': [0.331670, 0.196842, 0.093422],
            'N': [0.686311, 0.463189, 0.138722],
            'O': [0.606344, 0.463006, 0.142262],
            'S': [0.606344, 0.463006, 0.142262],
            'P': [0.418365, 0.290054, 0.1064245],
        }

    @staticmethod
    def _get_gbn2_default():
        return [1.0, 0.8, 4.851]

    @staticmethod
    def _get_i_simple():
        return (
            'select(step(r+sr2-or1), 0.5*(1/L-1/U+0.25*(r-sr2^2/r)*(1/(U^2)-1/(L^2))+0.5*log(L/U)/r), 0);'
            'U=r+sr2;'
            'L=max(or1, D);'
            'D=abs(r-sr2)')

    @staticmethod
    def build(model,
              topology,
              nb_force,
              sv_dielectric=78.39,
              sl_dielectric=1.0):
        """Return a configured CustomGBForce for *model*.

        :param model: one of 'hct', 'obc1', 'obc2', 'gbn', 'gbn2'
        :param topology: openmm.app.Topology with all atoms present
        :param nb_force: the NonbondedForce in the system (charges are read from it)
        :param sv_dielectric: solvent dielectric constant
        :param sl_dielectric: solute dielectric constant
        :returns: configured CustomGBForce (not yet added to any System)
        """
        if 'openmm' not in sys.modules:
            raise ImportError(
                'openmm is required to build implicit solvent forces.')

        key = model.lower()
        if key not in GBImplicitSolvent.get_valid_models():
            raise ValueError(
                f'Unknown GB model "{model}". '
                f'Valid options: {list(GBImplicitSolvent.get_valid_models())}')

        sv, sl = sv_dielectric, sl_dielectric
        force = mm.CustomGBForce()

        if key == 'hct':
            GBImplicitSolvent._build_hct(force, topology, nb_force, sv, sl)
        elif key == 'obc1':
            GBImplicitSolvent._build_obc1(force, topology, nb_force, sv, sl)
        elif key == 'obc2':
            GBImplicitSolvent._build_obc2(force, topology, nb_force, sv, sl)
        elif key == 'gbn':
            GBImplicitSolvent._build_gbn(force, topology, nb_force, sv, sl)
        elif key == 'gbn2':
            GBImplicitSolvent._build_gbn2(force, topology, nb_force, sv, sl)

        force.setNonbondedMethod(mm.CustomGBForce.NoCutoff)
        return force

    # ------------------------------------------------------------------
    # Per-model builders
    # ------------------------------------------------------------------

    @staticmethod
    def _build_hct(force, topology, nb_force, sv, sl):
        offset = 0.009
        radii = GBImplicitSolvent._mbondi_radii(topology)
        force.addPerParticleParameter('charge')
        force.addPerParticleParameter('or')
        force.addPerParticleParameter('sr')
        force.addComputedValue('I', GBImplicitSolvent._get_i_simple(),
                               mm.CustomGBForce.ParticlePairNoExclusions)
        force.addComputedValue('B', '1/(1/or-I)',
                               mm.CustomGBForce.SingleParticle)
        GBImplicitSolvent._add_energy_terms(force, sv, sl, offset)
        for i, atom in enumerate(topology.atoms()):
            charge = nb_force.getParticleParameters(i)[0]
            or_ = radii[i] / 10.0 - offset
            sr = GBImplicitSolvent._screen(atom, 0) * or_
            force.addParticle([charge, or_, sr])

    @staticmethod
    def _build_obc1(force, topology, nb_force, sv, sl):
        offset = 0.009
        radii = GBImplicitSolvent._mbondi2_radii(topology)
        force.addPerParticleParameter('charge')
        force.addPerParticleParameter('or')
        force.addPerParticleParameter('sr')
        force.addComputedValue('I', GBImplicitSolvent._get_i_simple(),
                               mm.CustomGBForce.ParticlePairNoExclusions)
        force.addComputedValue(
            'B', '1/(1/or-tanh(0.8*psi+2.909125*psi^3)/radius);'
            'psi=I*or; radius=or+offset; offset=0.009',
            mm.CustomGBForce.SingleParticle)
        GBImplicitSolvent._add_energy_terms(force, sv, sl, offset)
        for i, atom in enumerate(topology.atoms()):
            charge = nb_force.getParticleParameters(i)[0]
            or_ = radii[i] / 10.0 - offset
            sr = GBImplicitSolvent._screen(atom, 0) * or_
            force.addParticle([charge, or_, sr])

    @staticmethod
    def _build_obc2(force, topology, nb_force, sv, sl):
        offset = 0.009
        radii = GBImplicitSolvent._mbondi2_radii(topology)
        force.addPerParticleParameter('charge')
        force.addPerParticleParameter('or')
        force.addPerParticleParameter('sr')
        force.addComputedValue('I', GBImplicitSolvent._get_i_simple(),
                               mm.CustomGBForce.ParticlePairNoExclusions)
        force.addComputedValue(
            'B', '1/(1/or-tanh(psi-0.8*psi^2+4.85*psi^3)/radius);'
            'psi=I*or; radius=or+offset; offset=0.009',
            mm.CustomGBForce.SingleParticle)
        GBImplicitSolvent._add_energy_terms(force, sv, sl, offset)
        for i, atom in enumerate(topology.atoms()):
            charge = nb_force.getParticleParameters(i)[0]
            or_ = radii[i] / 10.0 - offset
            sr = GBImplicitSolvent._screen(atom, 0) * or_
            force.addParticle([charge, or_, sr])

    @staticmethod
    def _build_gbn(force, topology, nb_force, sv, sl):
        offset = 0.009
        radii = GBImplicitSolvent._bondi_radii(topology)
        atoms = list(topology.atoms())
        ors = [radii[i] / 10.0 - offset for i in range(len(atoms))]
        unique_ors = sorted(set(ors))
        or_to_idx = {r: i for i, r in enumerate(unique_ors)}
        n = len(unique_ors)
        d0_tab, m0_tab = GBImplicitSolvent._neck_tables(unique_ors, offset)

        force.addPerParticleParameter('charge')
        force.addPerParticleParameter('or')
        force.addPerParticleParameter('sr')
        force.addPerParticleParameter('radindex')
        force.addTabulatedFunction('getd0', mm.Discrete2DFunction(n, n, d0_tab))
        force.addTabulatedFunction('getm0', mm.Discrete2DFunction(n, n, m0_tab))
        force.addComputedValue(
            'I', 'Ivdw+neckScale*Ineck;'
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
        force.addComputedValue(
            'B',
            '1/(1/or-tanh(1.09511284*psi-1.907992938*psi^2+2.50798245*psi^3)/radius);'
            'psi=I*or; radius=or+offset; offset=0.009',
            mm.CustomGBForce.SingleParticle)
        GBImplicitSolvent._add_energy_terms(force, sv, sl, offset)
        for i, atom in enumerate(atoms):
            charge = nb_force.getParticleParameters(i)[0]
            or_ = ors[i]
            sr = GBImplicitSolvent._screen(atom, 1) * or_
            force.addParticle([charge, or_, sr, or_to_idx[or_]])

    @staticmethod
    def _build_gbn2(force, topology, nb_force, sv, sl):
        offset = 0.0195141
        radii = GBImplicitSolvent._mbondi3_radii(topology)
        atoms = list(topology.atoms())
        ors = [radii[i] / 10.0 - offset for i in range(len(atoms))]
        unique_ors = sorted(set(ors))
        or_to_idx = {r: i for i, r in enumerate(unique_ors)}
        n = len(unique_ors)
        d0_tab, m0_tab = GBImplicitSolvent._neck_tables(unique_ors, offset)

        force.addPerParticleParameter('charge')
        force.addPerParticleParameter('or')
        force.addPerParticleParameter('sr')
        force.addPerParticleParameter('alpha')
        force.addPerParticleParameter('beta')
        force.addPerParticleParameter('gamma')
        force.addPerParticleParameter('radindex')
        force.addTabulatedFunction('getd0', mm.Discrete2DFunction(n, n, d0_tab))
        force.addTabulatedFunction('getm0', mm.Discrete2DFunction(n, n, m0_tab))
        force.addComputedValue(
            'I', 'Ivdw+neckScale*Ineck;'
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
        force.addComputedValue(
            'B', '1/(1/or-tanh(alpha*psi-beta*psi^2+gamma*psi^3)/radius);'
            'psi=I*or; radius=or+offset; offset=0.0195141',
            mm.CustomGBForce.SingleParticle)
        GBImplicitSolvent._add_energy_terms(force, sv, sl, offset)
        nucleic_residues = GBImplicitSolvent._get_nucleic_residues()
        gbn2_params = GBImplicitSolvent._get_gbn2_params()
        gbn2_params_nucleic = GBImplicitSolvent._get_gbn2_params_nucleic()
        gbn2_default = GBImplicitSolvent._get_gbn2_default()
        for i, atom in enumerate(atoms):
            charge = nb_force.getParticleParameters(i)[0]
            or_ = ors[i]
            sym = atom.element.symbol if atom.element is not None else ''
            if atom.residue.name in nucleic_residues:
                abg = gbn2_params_nucleic.get(sym, gbn2_default)
                sr = GBImplicitSolvent._screen(atom, 3) * or_
            else:
                abg = gbn2_params.get(sym, gbn2_default)
                sr = GBImplicitSolvent._screen(atom, 2) * or_
            force.addParticle(
                [charge, or_, sr, abg[0], abg[1], abg[2], or_to_idx[or_]])

    # ------------------------------------------------------------------
    # Radius helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _bonded_atoms(topology):
        bonds = defaultdict(list)
        for a1, a2 in topology.bonds():
            bonds[a1].append(a2)
            bonds[a2].append(a1)
        return bonds

    @staticmethod
    def _is_carboxylate_O(atom, bonds):
        if atom.element is None or atom.element.symbol != 'O':
            return False
        nbrs = bonds[atom]
        if len(nbrs) != 1 or nbrs[0].element is None or nbrs[
                0].element.symbol != 'C':
            return False
        c_nbrs = bonds[nbrs[0]]
        if len(c_nbrs) != 3:
            return False
        for a3 in c_nbrs:
            if a3 is not atom and a3.element is not None and a3.element.symbol == 'O':
                return True
        return False

    @staticmethod
    def _bondi_radii(topology):
        table = {
            'C': 1.7,
            'H': 1.2,
            'D': 1.2,
            'N': 1.55,
            'O': 1.5,
            'F': 1.5,
            'Si': 2.1,
            'P': 1.85,
            'S': 1.8,
            'Cl': 1.7
        }
        return [
            table.get(a.element.symbol if a.element else '', 1.5)
            for a in topology.atoms()
        ]

    @staticmethod
    def _mbondi_radii(topology):
        bonds = GBImplicitSolvent._bonded_atoms(topology)
        const = {
            'N': 1.55,
            'O': 1.5,
            'F': 1.5,
            'Si': 2.1,
            'P': 1.85,
            'S': 1.8,
            'Cl': 1.7
        }
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
                r = const.get(sym, 1.5)
            radii.append(r)
        return radii

    @staticmethod
    def _mbondi2_radii(topology):
        bonds = GBImplicitSolvent._bonded_atoms(topology)
        const = {
            'N': 1.55,
            'O': 1.5,
            'F': 1.5,
            'Si': 2.1,
            'P': 1.85,
            'S': 1.8,
            'Cl': 1.7
        }
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
                r = const.get(sym, 1.5)
            radii.append(r)
        return radii

    @staticmethod
    def _mbondi3_radii(topology):
        bonds = GBImplicitSolvent._bonded_atoms(topology)
        radii = GBImplicitSolvent._mbondi2_radii(topology)
        for i, atom in enumerate(topology.atoms()):
            if GBImplicitSolvent._is_carboxylate_O(atom, bonds):
                radii[i] = 1.4
            elif atom.residue.name == 'ARG':
                name = atom.name
                if name.startswith('HH') or name.startswith('HE'):
                    radii[i] = 1.17
        return radii

    # ------------------------------------------------------------------
    # Neck table and screen helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _neck_tables(unique_ors, offset):
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

        return _reduce(GBImplicitSolvent._get_d0()), _reduce(
            GBImplicitSolvent._get_m0())

    @staticmethod
    def _screen(atom, idx):
        sym = atom.element.symbol if atom.element is not None else ''
        return GBImplicitSolvent._get_screen().get(
            sym, GBImplicitSolvent._get_screen_default())[idx]

    @staticmethod
    def _add_energy_terms(force, sv, sl, offset):
        params = (f'; solventDielectric={sv:.16g}'
                  f'; soluteDielectric={sl:.16g}'
                  f'; kappa=0; offset={offset:.16g}')
        force.addEnergyTerm(
            '-0.5*138.935485*(1/soluteDielectric-1/solventDielectric)*charge^2/B'
            + params, mm.CustomGBForce.SingleParticle)
        force.addEnergyTerm(
            '28.3919551*(radius+0.14)^2*(radius/B)^6; radius=or+offset' +
            params, mm.CustomGBForce.SingleParticle)
        force.addEnergyTerm(
            '-138.935485*(1/soluteDielectric-1/solventDielectric)*charge1*charge2/f;'
            'f=sqrt(r^2+B1*B2*exp(-r^2/(4*B1*B2)))' + params,
            mm.CustomGBForce.ParticlePairNoExclusions)
