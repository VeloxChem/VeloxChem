//
// This file contains derivative work of multicharge 0.3.0,
// which is licensed under the Apache License, Version 2.0.
//
// SPDX-License-Identifier: Apache-2.0
//
// Copyright 2021 Sebastian Ehlert
//
// Modifications made by VeloxChem developer, 2024, for use in C++.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include "AtomicPartialChargesModel.hpp"

#include <cmath>

#include "Codata.hpp"
#include "DenseLinearAlgebra.hpp"
#include "DenseMatrix.hpp"
#include "Molecule.hpp"

namespace atmparchg {  // atmparchg namespace

auto
getPartialCharges(const CMolecule& molecule, const double netcharge) -> std::vector<double>
{
    // Electronegativity equilibration charge model
    // E. Caldeweyher, S. Ehlert, A. Hansen, H. Neugebauer, S. Spicher, C. Bannwarth, S. Grimme,
    // J. Chem. Phys., 2019, 150, 154122.

    // Element-specific electronegativity for the electronegativity equilibration charges

    std::vector<double> eeq_chi({0.00000000, 1.23695041, 1.26590957, 0.54341808, 0.99666991, 1.26691604, 1.40028282, 1.55819364, 1.56866440, 1.57540015,
                                 1.15056627, 0.55936220, 0.72373742, 1.12910844, 1.12306840, 1.52672442, 1.40768172, 1.48154584, 1.31062963, 0.40374140,
                                 0.75442607, 0.76482096, 0.98457281, 0.96702598, 1.05266584, 0.93274875, 1.04025281, 0.92738624, 1.07419210, 1.07900668,
                                 1.04712861, 1.15018618, 1.15388455, 1.36313743, 1.36485106, 1.39801837, 1.18695346, 0.36273870, 0.58797255, 0.71961946,
                                 0.96158233, 0.89585296, 0.81360499, 1.00794665, 0.92613682, 1.09152285, 1.14907070, 1.13508911, 1.08853785, 1.11005982,
                                 1.12452195, 1.21642129, 1.36507125, 1.40340000, 1.16653482, 0.34125098, 0.58884173, 0.68441115, 0.56999999, 0.56999999,
                                 0.56999999, 0.56999999, 0.56999999, 0.56999999, 0.56999999, 0.56999999, 0.56999999, 0.56999999, 0.56999999, 0.56999999,
                                 0.56999999, 0.56999999, 0.87936784, 1.02761808, 0.93297476, 1.10172128, 0.97350071, 1.16695666, 1.23997927, 1.18464453,
                                 1.14191734, 1.12334192, 1.01485321, 1.12950808, 1.30804834, 1.33689961, 1.27465977});

    // Element-specific chemical hardnesses for the electronegativity equilibration charges

    std::vector<double> eeq_eta(
        {0.00000000,  -0.35015861, 1.04121227, 0.09281243,  0.09412380,  0.26629137,  0.19408787,  0.05317918,  0.03151644,  0.32275132,  1.30996037,
         0.24206510,  0.04147733,  0.11634126, 0.13155266,  0.15350650,  0.15250997,  0.17523529,  0.28774450,  0.42937314,  0.01896455,  0.07179178,
         -0.01121381, -0.03093370, 0.02716319, -0.01843812, -0.15270393, -0.09192645, -0.13418723, -0.09861139, 0.18338109,  0.08299615,  0.11370033,
         0.19005278,  0.10980677,  0.12327841, 0.25345554,  0.58615231,  0.16093861,  0.04548530,  -0.02478645, 0.01909943,  0.01402541,  -0.03595279,
         0.01137752,  -0.03697213, 0.08009416, 0.02274892,  0.12801822,  -0.02078702, 0.05284319,  0.07581190,  0.09663758,  0.09547417,  0.07803344,
         0.64913257,  0.15348654,  0.05054344, 0.11000000,  0.11000000,  0.11000000,  0.11000000,  0.11000000,  0.11000000,  0.11000000,  0.11000000,
         0.11000000,  0.11000000,  0.11000000, 0.11000000,  0.11000000,  0.11000000,  -0.02786741, 0.01057858,  -0.03892226, -0.04574364, -0.03874080,
         -0.03782372, -0.07046855, 0.09546597, 0.21953269,  0.02522348,  0.15263050,  0.08042611,  0.01878626,  0.08715453,  0.10500484});

    // Element-specific CN scaling constant for the electronegativity equilibration charges

    std::vector<double> eeq_kcn({0.00000000,  0.04916110,  0.10937243,  -0.12349591, -0.02665108, -0.02631658, 0.06005196,  0.09279548,  0.11689703,
                                 0.15704746,  0.07987901,  -0.10002962, -0.07712863, -0.02170561, -0.04964052, 0.14250599,  0.07126660,  0.13682750,
                                 0.14877121,  -0.10219289, -0.08979338, -0.08273597, -0.01754829, -0.02765460, -0.02558926, -0.08010286, -0.04163215,
                                 -0.09369631, -0.03774117, -0.05759708, 0.02431998,  -0.01056270, -0.02692862, 0.07657769,  0.06561608,  0.08006749,
                                 0.14139200,  -0.05351029, -0.06701705, -0.07377246, -0.02927768, -0.03867291, -0.06929825, -0.04485293, -0.04800824,
                                 -0.01484022, 0.07917502,  0.06619243,  0.02434095,  -0.01505548, -0.03030768, 0.01418235,  0.08953411,  0.08967527,
                                 0.07277771,  -0.02129476, -0.06188828, -0.06568203, -0.11000000, -0.11000000, -0.11000000, -0.11000000, -0.11000000,
                                 -0.11000000, -0.11000000, -0.11000000, -0.11000000, -0.11000000, -0.11000000, -0.11000000, -0.11000000, -0.11000000,
                                 -0.03585873, -0.03132400, -0.05902379, -0.02827592, -0.07606260, -0.02123839, 0.03814822,  0.02146834,  0.01580538,
                                 -0.00894298, -0.05864876, -0.01817842, 0.07721851,  0.07936083,  0.05849285});

    // Element-specific charge widths for the electronegativity equilibration charges

    std::vector<double> eeq_rad({0.00000000, 0.55159092, 0.66205886, 0.90529132, 1.51710827, 2.86070364, 1.88862966, 1.32250290, 1.23166285, 1.77503721,
                                 1.11955204, 1.28263182, 1.22344336, 1.70936266, 1.54075036, 1.38200579, 2.18849322, 1.36779065, 1.27039703, 1.64466502,
                                 1.58859404, 1.65357953, 1.50021521, 1.30104175, 1.46301827, 1.32928147, 1.02766713, 1.02291377, 0.94343886, 1.14881311,
                                 1.47080755, 1.76901636, 1.98724061, 2.41244711, 2.26739524, 2.95378999, 1.20807752, 1.65941046, 1.62733880, 1.61344972,
                                 1.63220728, 1.60899928, 1.43501286, 1.54559205, 1.32663678, 1.37644152, 1.36051851, 1.23395526, 1.65734544, 1.53895240,
                                 1.97542736, 1.97636542, 2.05432381, 3.80138135, 1.43893803, 1.75505957, 1.59815118, 1.76401732, 1.63999999, 1.63999999,
                                 1.63999999, 1.63999999, 1.63999999, 1.63999999, 1.63999999, 1.63999999, 1.63999999, 1.63999999, 1.63999999, 1.63999999,
                                 1.63999999, 1.63999999, 1.47055223, 1.81127084, 1.40189963, 1.54015481, 1.33721475, 1.57165422, 1.04815857, 1.78342098,
                                 2.79106396, 1.78160840, 2.47588882, 2.37670734, 1.76613217, 2.66172302, 2.82773085});

    auto natoms = molecule.number_of_atoms();

    auto ids_elem = molecule.identifiers();

    std::vector<double> chg_chi(natoms), chg_eta(natoms), chg_kcn(natoms), chg_rad(natoms);

    for (int i = 0; i < natoms; i++)
    {
        chg_chi[i] = eeq_chi[ids_elem[i]];
        chg_eta[i] = eeq_eta[ids_elem[i]];
        chg_kcn[i] = eeq_kcn[ids_elem[i]];
        chg_rad[i] = eeq_rad[ids_elem[i]];
    }

    // form right-hand side vector
    // source: get_vrhs

    std::vector<double> vec_rhs(natoms + 1);

    auto cn = atmparchg::getCoordinationNumber(molecule);

    for (int i = 0; i < natoms; i++)
    {
        double val = chg_kcn[i] / (std::sqrt(cn[i]) + 1.0e-14);

        vec_rhs[i] = -chg_chi[i] + val * cn[i];
    }

    vec_rhs[natoms] = netcharge;

    // form left-hand side matrix
    // source: get_amat_0d

    CDenseMatrix Amat(natoms + 1, natoms + 1);

    const double sqrtpi = std::sqrt(mathconst::pi_value());

    const double sqrt2pi = std::sqrt(2.0 / mathconst::pi_value());

    auto xyzcoord = molecule.coordinates("bohr");

    for (int i = 0; i < natoms; i++)
    {
        const auto rixyz = xyzcoord[i].coordinates();

        auto Amat_row_i = Amat.row(i);

        for (int j = 0; j < i; j++)
        {
            const auto rjxyz = xyzcoord[j].coordinates();

            auto Amat_row_j = Amat.row(j);

            std::vector<double> rij({rixyz[0] - rjxyz[0], rixyz[1] - rjxyz[1], rixyz[2] - rjxyz[2]});

            double r2 = rij[0] * rij[0] + rij[1] * rij[1] + rij[2] * rij[2];

            double gam = 1.0 / (chg_rad[i] * chg_rad[i] + chg_rad[j] * chg_rad[j]);

            double val = std::erf(std::sqrt(gam * r2)) / std::sqrt(r2);

            Amat_row_j[i] = val;

            Amat_row_i[j] = val;
        }

        Amat_row_i[i] = chg_eta[i] + sqrt2pi / chg_rad[i];
    }

    for (int i = 0; i < natoms; i++)
    {
        Amat.row(i)[natoms] = 1.0;

        Amat.row(natoms)[i] = 1.0;
    }

    Amat.row(natoms)[natoms] = 0.0;

    // get partial charges

    auto solutions = denblas::serialSolve(Amat, vec_rhs);

    auto solptr = solutions.data();

    std::vector<double> partialcharges(solptr, solptr + natoms);

    return partialcharges;
}

auto
getCoordinationNumber(const CMolecule& molecule) -> std::vector<double>
{
    // Covalent radii (taken from Pyykko and Atsumi, Chem. Eur. J. 15, 2009,
    // 188-197), values for metals decreased by 10 %

    std::vector<double> covalent_radii({
        0.0,                                             // dummy atom
        0.32, 0.46,                                      // H,He
        1.20, 0.94, 0.77, 0.75, 0.71, 0.63, 0.64, 0.67,  // Li-Ne
        1.40, 1.25, 1.13, 1.04, 1.10, 1.02, 0.99, 0.96,  // Na-Ar
        1.76, 1.54,                                      // K,Ca
        1.33, 1.22, 1.21, 1.10, 1.07,                    // Sc-
        1.04, 1.00, 0.99, 1.01, 1.09,                    // -Zn
        1.12, 1.09, 1.15, 1.10, 1.14, 1.17,              // Ga-Kr
        1.89, 1.67,                                      // Rb,Sr
        1.47, 1.39, 1.32, 1.24, 1.15,                    // Y-
        1.13, 1.13, 1.08, 1.15, 1.23,                    // -Cd
        1.28, 1.26, 1.26, 1.23, 1.32, 1.31,              // In-Xe
        2.09, 1.76,                                      // Cs,Ba
        1.62, 1.47, 1.58, 1.57, 1.56, 1.55, 1.51,        // La-Eu
        1.52, 1.51, 1.50, 1.49, 1.49, 1.48, 1.53,        // Gd-Yb
        1.46, 1.37, 1.31, 1.23, 1.18,                    // Lu-
        1.16, 1.11, 1.12, 1.13, 1.32,                    // -Hg
        1.30, 1.30, 1.36, 1.31, 1.38, 1.42,              // Tl-Rn
        2.01, 1.81,                                      // Fr,Ra
        1.67, 1.58, 1.52, 1.53, 1.54, 1.55, 1.49,        // Ac-Am
        1.49, 1.51, 1.51, 1.48, 1.50, 1.56, 1.58,        // Cm-No
        1.45, 1.41, 1.34, 1.29, 1.27,                    // Lr-
        1.21, 1.16, 1.15, 1.09, 1.22,                    // -Cn
        1.36, 1.43, 1.46, 1.58, 1.48, 1.57               // Nh-Og
    });

    for (size_t i = 0; i < covalent_radii.size(); i++)
    {
        covalent_radii[i] /= units::bohr_in_angstrom();

        // D3 covalent radii used to construct the coordination number
        covalent_radii[i] *= 4.0 / 3.0;
    }

    // prepare parameters

    const double kcn = 7.50;

    const double cn_max = 8.0;

    const double cutoff2 = 25.0 * 25.0;

    const double sqrtpi = std::sqrt(mathconst::pi_value());

    // get molecular information

    auto natoms = molecule.number_of_atoms();

    auto ids_elem = molecule.identifiers();

    auto xyzcoord = molecule.coordinates("bohr");

    // compute coordination numbers with error function
    // source: ncoord_erf

    std::vector<double> cn(natoms, 0.0);

    for (int i = 0; i < natoms; i++)
    {
        auto rcov_i = covalent_radii[ids_elem[i]];

        const auto rixyz = xyzcoord[i].coordinates();

        for (int j = 0; j < i; j++)
        {
            auto rcov_j = covalent_radii[ids_elem[j]];

            const auto rjxyz = xyzcoord[j].coordinates();

            std::vector<double> rij({rjxyz[0] - rixyz[0], rjxyz[1] - rixyz[1], rjxyz[2] - rixyz[2]});

            double r2 = rij[0] * rij[0] + rij[1] * rij[1] + rij[2] * rij[2];

            if (r2 > cutoff2) continue;

            double r = std::sqrt(r2);

            double rc = rcov_i + rcov_j;

            double cn_val = 0.5 * (1.0 + std::erf(-kcn * (r - rc) / rc));

            cn[i] += cn_val;

            cn[j] += cn_val;
        }
    }

    // apply cutoff function for large coordination numbers
    // source: log_cn_cut

    for (int i = 0; i < natoms; i++)
    {
        cn[i] = std::log(1.0 + std::exp(cn_max)) - std::log(1.0 + std::exp(cn_max - cn[i]));
    }

    return cn;
}

}  // namespace atmparchg
