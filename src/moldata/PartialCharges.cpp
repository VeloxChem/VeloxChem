//
//                           VELOXCHEM 1.0-RC
//         ----------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2021 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact
//
//  SPDX-License-Identifier: LGPL-3.0-or-later
//
//  This file is part of VeloxChem.
//
//  VeloxChem is free software: you can redistribute it and/or modify it under
//  the terms of the GNU Lesser General Public License as published by the Free
//  Software Foundation, either version 3 of the License, or (at your option)
//  any later version.
//
//  VeloxChem is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
//  License for more details.
//
//  You should have received a copy of the GNU Lesser General Public License
//  along with VeloxChem. If not, see <https://www.gnu.org/licenses/>.
//
//  This file contains derivative work of dftd4 (v2.4.0):
//  Copyright © 2017-2019 Stefan Grimme, Sebastian Ehlert, Eike Caldeweyher

#include "PartialCharges.hpp"

#include "CoordinationNumber.hpp"
#include "DenseDiagonalizer.hpp"
#include "DenseLinearAlgebra.hpp"
#include "DenseMatrix.hpp"
#include "ErrorHandler.hpp"
#include "MathConst.hpp"
#include "MathFunc.hpp"
#include "Molecule.hpp"

namespace parchg {  // parchg namespace

std::vector<double>
getPartialCharges(const CMolecule& molecule, const double netcharge)
{
    CDenseMatrix dqdr;

    return getPartialCharges(molecule, netcharge, dqdr);
}

std::vector<double>
getPartialCharges(const CMolecule& molecule, const double netcharge, CDenseMatrix& dqdr)
{
    // Reference: dftd4 (v2.4.0)

    // prepare parameters

    std::vector<double> en({0.00000000, 1.23695041, 1.26590957, 0.54341808, 0.99666991, 1.26691604, 1.40028282, 1.55819364, 1.56866440, 1.57540015,
                            1.15056627, 0.55936220, 0.72373742, 1.12910844, 1.12306840, 1.52672442, 1.40768172, 1.48154584, 1.31062963, 0.40374140,
                            0.75442607, 0.76482096, 0.98457281, 0.96702598, 1.05266584, 0.93274875, 1.04025281, 0.92738624, 1.07419210, 1.07900668,
                            1.04712861, 1.15018618, 1.15388455, 1.36313743, 1.36485106, 1.39801837, 1.18695346, 0.36273870, 0.58797255, 0.71961946,
                            0.96158233, 0.89585296, 0.81360499, 1.00794665, 0.92613682, 1.09152285, 1.14907070, 1.13508911, 1.08853785, 1.11005982,
                            1.12452195, 1.21642129, 1.36507125, 1.40340000, 1.16653482, 0.34125098, 0.58884173, 0.68441115, 0.56999999, 0.56999999,
                            0.56999999, 0.56999999, 0.56999999, 0.56999999, 0.56999999, 0.56999999, 0.56999999, 0.56999999, 0.56999999, 0.56999999,
                            0.56999999, 0.56999999, 0.87936784, 1.02761808, 0.93297476, 1.10172128, 0.97350071, 1.16695666, 1.23997927, 1.18464453,
                            1.14191734, 1.12334192, 1.01485321, 1.12950808, 1.30804834, 1.33689961, 1.27465977});

    std::vector<double> gamm(
        {0.00000000,  -0.35015861, 1.04121227, 0.09281243,  0.09412380,  0.26629137,  0.19408787,  0.05317918,  0.03151644,  0.32275132,  1.30996037,
         0.24206510,  0.04147733,  0.11634126, 0.13155266,  0.15350650,  0.15250997,  0.17523529,  0.28774450,  0.42937314,  0.01896455,  0.07179178,
         -0.01121381, -0.03093370, 0.02716319, -0.01843812, -0.15270393, -0.09192645, -0.13418723, -0.09861139, 0.18338109,  0.08299615,  0.11370033,
         0.19005278,  0.10980677,  0.12327841, 0.25345554,  0.58615231,  0.16093861,  0.04548530,  -0.02478645, 0.01909943,  0.01402541,  -0.03595279,
         0.01137752,  -0.03697213, 0.08009416, 0.02274892,  0.12801822,  -0.02078702, 0.05284319,  0.07581190,  0.09663758,  0.09547417,  0.07803344,
         0.64913257,  0.15348654,  0.05054344, 0.11000000,  0.11000000,  0.11000000,  0.11000000,  0.11000000,  0.11000000,  0.11000000,  0.11000000,
         0.11000000,  0.11000000,  0.11000000, 0.11000000,  0.11000000,  0.11000000,  -0.02786741, 0.01057858,  -0.03892226, -0.04574364, -0.03874080,
         -0.03782372, -0.07046855, 0.09546597, 0.21953269,  0.02522348,  0.15263050,  0.08042611,  0.01878626,  0.08715453,  0.10500484});

    std::vector<double> cnfak({0.00000000,  0.04916110,  0.10937243,  -0.12349591, -0.02665108, -0.02631658, 0.06005196,  0.09279548,  0.11689703,
                               0.15704746,  0.07987901,  -0.10002962, -0.07712863, -0.02170561, -0.04964052, 0.14250599,  0.07126660,  0.13682750,
                               0.14877121,  -0.10219289, -0.08979338, -0.08273597, -0.01754829, -0.02765460, -0.02558926, -0.08010286, -0.04163215,
                               -0.09369631, -0.03774117, -0.05759708, 0.02431998,  -0.01056270, -0.02692862, 0.07657769,  0.06561608,  0.08006749,
                               0.14139200,  -0.05351029, -0.06701705, -0.07377246, -0.02927768, -0.03867291, -0.06929825, -0.04485293, -0.04800824,
                               -0.01484022, 0.07917502,  0.06619243,  0.02434095,  -0.01505548, -0.03030768, 0.01418235,  0.08953411,  0.08967527,
                               0.07277771,  -0.02129476, -0.06188828, -0.06568203, -0.11000000, -0.11000000, -0.11000000, -0.11000000, -0.11000000,
                               -0.11000000, -0.11000000, -0.11000000, -0.11000000, -0.11000000, -0.11000000, -0.11000000, -0.11000000, -0.11000000,
                               -0.03585873, -0.03132400, -0.05902379, -0.02827592, -0.07606260, -0.02123839, 0.03814822,  0.02146834,  0.01580538,
                               -0.00894298, -0.05864876, -0.01817842, 0.07721851,  0.07936083,  0.05849285});

    std::vector<double> alp({0.00000000, 0.55159092, 0.66205886, 0.90529132, 1.51710827, 2.86070364, 1.88862966, 1.32250290, 1.23166285, 1.77503721,
                             1.11955204, 1.28263182, 1.22344336, 1.70936266, 1.54075036, 1.38200579, 2.18849322, 1.36779065, 1.27039703, 1.64466502,
                             1.58859404, 1.65357953, 1.50021521, 1.30104175, 1.46301827, 1.32928147, 1.02766713, 1.02291377, 0.94343886, 1.14881311,
                             1.47080755, 1.76901636, 1.98724061, 2.41244711, 2.26739524, 2.95378999, 1.20807752, 1.65941046, 1.62733880, 1.61344972,
                             1.63220728, 1.60899928, 1.43501286, 1.54559205, 1.32663678, 1.37644152, 1.36051851, 1.23395526, 1.65734544, 1.53895240,
                             1.97542736, 1.97636542, 2.05432381, 3.80138135, 1.43893803, 1.75505957, 1.59815118, 1.76401732, 1.63999999, 1.63999999,
                             1.63999999, 1.63999999, 1.63999999, 1.63999999, 1.63999999, 1.63999999, 1.63999999, 1.63999999, 1.63999999, 1.63999999,
                             1.63999999, 1.63999999, 1.47055223, 1.81127084, 1.40189963, 1.54015481, 1.33721475, 1.57165422, 1.04815857, 1.78342098,
                             2.79106396, 1.78160840, 2.47588882, 2.37670734, 1.76613217, 2.66172302, 2.82773085});

    auto natoms = molecule.getNumberOfAtoms();

    auto ids_elem = molecule.getIdsElemental();

    std::vector<double> chg_xi, chg_gam, chg_kappa, chg_alpha;

    for (int32_t i = 0; i < natoms; i++)
    {
        chg_xi.push_back(en[ids_elem[i]]);

        chg_gam.push_back(gamm[ids_elem[i]]);

        chg_kappa.push_back(cnfak[ids_elem[i]]);

        chg_alpha.push_back(std::pow(alp[ids_elem[i]], 2));
    }

    // form right-hand side vector

    int32_t natoms_1 = natoms + 1;

    CDenseMatrix Xvec(natoms_1, 1);

    CDenseMatrix Xfac(natoms, 1);

    CDenseMatrix dcndr(3 * natoms, natoms);

    auto cn = coordnum::getCoordinationNumber(molecule, dcndr);

    for (int32_t i = 0; i < natoms; i++)
    {
        double val = chg_kappa[i] / (std::sqrt(cn[i]) + 1.0e-14);

        Xvec.values()[i] = -chg_xi[i] + val * cn[i];

        Xfac.values()[i] = 0.5 * val;
    }

    Xvec.values()[natoms] = netcharge;

    // form left-hand side matrix

    CDenseMatrix Amat(natoms_1, natoms_1);

    const double sqrtpi = std::sqrt(mathconst::getPiValue());

    const double sqrt2pi = std::sqrt(2.0 / mathconst::getPiValue());

    auto xcoord = molecule.getCoordinatesX();

    auto ycoord = molecule.getCoordinatesY();

    auto zcoord = molecule.getCoordinatesZ();

    for (int32_t i = 0; i < natoms; i++)
    {
        for (int32_t j = 0; j < i; j++)
        {
            auto r = mathfunc::distance(xcoord[j], ycoord[j], zcoord[j], xcoord[i], ycoord[i], zcoord[i]);

            double gamij = 1.0 / std::sqrt(chg_alpha[i] + chg_alpha[j]);

            double val = std::erf(gamij * r) / r;

            Amat.values()[j * natoms_1 + i] = val;

            Amat.values()[i * natoms_1 + j] = val;
        }

        Amat.values()[i * natoms_1 + i] = chg_gam[i] + sqrt2pi / std::sqrt(chg_alpha[i]);
    }

    for (int32_t i = 0; i < natoms; i++)
    {
        Amat.values()[i * natoms_1 + natoms] = 1.0;

        Amat.values()[natoms * natoms_1 + i] = 1.0;
    }

    Amat.values()[natoms * natoms_1 + natoms] = 0.0;

    // get partial charges

    CDenseDiagonalizer diagdrv;

    diagdrv.diagonalize(Amat);

    std::string err_diag("parchg::getPartialCharges: Matrix diagonalization failed");

    errors::assertMsgCritical(diagdrv.getState(), err_diag);

    auto Ainv = diagdrv.getInvertedMatrix();

    auto solution = denblas::multAB(Ainv, Xvec);

    std::vector<double> partialcharges(natoms);

    for (int32_t i = 0; i < natoms; i++)
    {
        partialcharges[i] = solution.values()[i];
    }

    if (dqdr.getNumberOfElements() == 0)
    {
        return partialcharges;
    }

    // compute derivative of partial charges

    std::string err_size("PartialCharges - Mismatch in dqdr matrix size");

    errors::assertMsgCritical(dqdr.getNumberOfRows() == 3 * natoms, err_size);

    errors::assertMsgCritical(dqdr.getNumberOfColumns() == natoms, err_size);

    CDenseMatrix dAmatdr(3 * natoms, natoms_1);

    CDenseMatrix dXvecdr(3 * natoms, natoms_1);

    CDenseMatrix Afac(3, natoms);

    for (int32_t i = 0; i < natoms; i++)
    {
        for (int32_t dj = 0; dj < 3 * natoms; dj++)
        {
            dXvecdr.values()[dj * natoms_1 + i] = dcndr.values()[dj * natoms + i] * Xfac.values()[i];
        }

        for (int32_t j = 0; j < i; j++)
        {
            std::vector<double> rij({xcoord[i] - xcoord[j], ycoord[i] - ycoord[j], zcoord[i] - zcoord[j]});

            double r2 = rij[0] * rij[0] + rij[1] * rij[1] + rij[2] * rij[2];

            double gamij = 1.0 / std::sqrt(chg_alpha[i] + chg_alpha[j]);

            double arg = gamij * gamij * r2;

            double dval = 2.0 * gamij * std::exp(-arg) / (sqrtpi * r2) - Amat.values()[j * natoms_1 + i] / r2;

            for (int32_t d = 0; d < 3; d++)
            {
                int32_t di = d * natoms + i;

                int32_t dj = d * natoms + j;

                Afac.values()[di] += dval * rij[d] * partialcharges[j];

                Afac.values()[dj] -= dval * rij[d] * partialcharges[i];

                dAmatdr.values()[di * natoms_1 + j] = dval * rij[d] * partialcharges[i];

                dAmatdr.values()[dj * natoms_1 + i] = -dval * rij[d] * partialcharges[j];
            }
        }
    }

    for (int32_t d = 0; d < 3; d++)
    {
        for (int32_t i = 0; i < natoms; i++)
        {
            int32_t di = d * natoms + i;

            dAmatdr.values()[di * natoms_1 + i] += Afac.values()[di];
        }
    }

    auto prod_dmat = denblas::multAB(dAmatdr, Ainv);

    auto prod_dvec = denblas::multAB(dXvecdr, Ainv);

    auto prod_sum = denblas::addAB(prod_dmat, prod_dvec, 1.0);

    for (int32_t dj = 0; dj < 3 * natoms; dj++)
    {
        for (int32_t i = 0; i < natoms; i++)
        {
            dqdr.values()[dj * natoms + i] = -prod_sum.values()[dj * natoms_1 + i];
        }
    }

    return partialcharges;
}

}  // namespace parchg
