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

#include "AtomBasisTest.hpp"

#include "AtomBasis.hpp"
#include "AtomBasisSetter.hpp"
#include "CheckFunctions.hpp"

TEST_F(CAtomBasisTest, DefaultConstructor)
{
    CAtomBasis abas;

    CAtomBasis bbas = vlxbas::getAtomBasisEmpty();

    ASSERT_EQ(abas, bbas);
}

TEST_F(CAtomBasisTest, CopyConstructor)
{
    CAtomBasis abas = vlxbas::getAtomBasisForLi();

    CAtomBasis bbas(abas);

    ASSERT_EQ(abas, bbas);
}

TEST_F(CAtomBasisTest, MoveConstructor)
{
    CAtomBasis abas = vlxbas::getAtomBasisForLi();

    CAtomBasis bbas(vlxbas::getAtomBasisForLi());

    ASSERT_EQ(abas, bbas);
}

TEST_F(CAtomBasisTest, CopyAssignment)
{
    CAtomBasis abas = vlxbas::getAtomBasisForLi();

    CAtomBasis bbas = abas;

    ASSERT_EQ(abas, bbas);
}

TEST_F(CAtomBasisTest, MoveAssignment)
{
    CAtomBasis abas = vlxbas::getAtomBasisForLi();

    ASSERT_EQ(abas, vlxbas::getAtomBasisForLi());
}

TEST_F(CAtomBasisTest, SetAndGetIdElemental)
{
    CAtomBasis abas = vlxbas::getAtomBasisForLi();

    ASSERT_EQ(3, abas.getIdElemental());

    abas.setIdElemental(2);

    ASSERT_EQ(2, abas.getIdElemental());
}

TEST_F(CAtomBasisTest, SetAndGetMaxAngularMomentum)
{
    CAtomBasis abas = vlxbas::getAtomBasisForLi();

    ASSERT_EQ(1, abas.getMaxAngularMomentum());

    abas.setMaxAngularMomentum(3);

    ASSERT_EQ(3, abas.getMaxAngularMomentum());
}

TEST_F(CAtomBasisTest, AddBasisFunction)
{
    CAtomBasis abas;

    abas.setIdElemental(1);

    abas.addBasisFunction(CBasisFunction(
        {1.301070100000e+01, 1.962257200000e+00, 4.445379600000e-01}, {1.968215800000e-02, 1.379652400000e-01, 4.783193500000e-01}, 0));

    abas.addBasisFunction(CBasisFunction({1.219496200000e-01}, {1.000000000000e+00}, 0));

    abas.addBasisFunction(CBasisFunction({8.000000000000e-01}, {1.000000000000e+00}, 1));

    ASSERT_EQ(abas, vlxbas::getAtomBasisForH());
}

TEST_F(CAtomBasisTest, GetNumberOfBasisFunctions)
{
    CAtomBasis abas = vlxbas::getAtomBasisForLi();

    ASSERT_EQ(0, abas.getNumberOfBasisFunctions(3));

    ASSERT_EQ(3, abas.getNumberOfBasisFunctions(0));

    ASSERT_EQ(2, abas.getNumberOfBasisFunctions(1));
}

TEST_F(CAtomBasisTest, GetNumberOfPrimitiveFunctions)
{
    CAtomBasis abas = vlxbas::getAtomBasisForLi();

    ASSERT_EQ(0, abas.getNumberOfPrimitiveFunctions(3));

    ASSERT_EQ(7, abas.getNumberOfPrimitiveFunctions(0));

    ASSERT_EQ(3, abas.getNumberOfPrimitiveFunctions(1));
}

TEST_F(CAtomBasisTest, GetContractionString)
{
    CAtomBasis abas = vlxbas::getAtomBasisForLi();

    ASSERT_EQ(std::string("(3S,2P)"), abas.getContractionString());
}

TEST_F(CAtomBasisTest, GetPrimitivesString)
{
    CAtomBasis abas = vlxbas::getAtomBasisForLi();

    ASSERT_EQ(std::string("(7S,3P)"), abas.getPrimitivesString());
}

TEST_F(CAtomBasisTest, GetBasisFunctions)
{
    CAtomBasis abas = vlxbas::getAtomBasisForLi();

    auto sbfs = abas.getBasisFunctions(0);

    ASSERT_EQ(3u, sbfs.size());

    ASSERT_EQ(sbfs[0],
              CBasisFunction({2.662778551600e+02, 4.006978344700e+01, 9.055994438900e+00, 2.450300905100e+00, 7.220957185500e-01},
                             {6.492015032500e-03, 4.774786321500e-02, 2.026879611100e-01, 4.860657481700e-01, 4.362697795500e-01},
                             0));

    ASSERT_EQ(sbfs[1], CBasisFunction({5.281088472100e-02}, {1.000000000000e+00}, 0));

    ASSERT_EQ(sbfs[2], CBasisFunction({2.096094879800e-02}, {1.000000000000e+00}, 0));

    auto pbfs = abas.getBasisFunctions(1);

    ASSERT_EQ(2u, pbfs.size());

    ASSERT_EQ(pbfs[0], CBasisFunction({1.450000000000e+00, 3.000000000000e-01}, {2.586000000000e-01, 1.000000000000e+00}, 1));

    ASSERT_EQ(pbfs[1], CBasisFunction({8.200000000000e-02}, {1.000000000000e+00}, 1));

    auto dbfs = abas.getBasisFunctions(2);

    ASSERT_EQ(0u, dbfs.size());
}

TEST_F(CAtomBasisTest, ReduceToValenceBasis)
{
    CAtomBasis abas = vlxbas::getAtomBasisForLi();

    CAtomBasis bbas = abas.reduceToValenceBasis();

    CAtomBasis cbas;

    cbas.setIdElemental(3);

    cbas.addBasisFunction(CBasisFunction({2.662778551600e+02, 4.006978344700e+01, 9.055994438900e+00, 2.450300905100e+00, 7.220957185500e-01},
                                         {6.492015032500e-03, 4.774786321500e-02, 2.026879611100e-01, 4.860657481700e-01, 4.362697795500e-01},
                                         0));

    cbas.addBasisFunction(CBasisFunction({5.281088472100e-02}, {1.000000000000e+00}, 0));

    cbas.addBasisFunction(CBasisFunction({2.096094879800e-02}, {1.000000000000e+00}, 0));

    ASSERT_EQ(bbas, cbas);
}
