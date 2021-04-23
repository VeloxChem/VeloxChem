//
//                           VELOXCHEM 1.0-RC2
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

#include "ScreeningContainerTest.hpp"

#include "ScreeningContainer.hpp"
#include "MoleculeSetter.hpp"
#include "MolecularBasisSetter.hpp"

TEST_F(CScreeningContainerTest, DefaultConstructor)
{
    CGtoPairsContainer bppcont;
    
    CVecMemBlock<double> bqvals;
    
    CScreeningContainer qconta;
    
    CScreeningContainer qcontb(bqvals, bqvals, bppcont, bppcont, ericut::qq,
                               1.0e-13);
    
    ASSERT_EQ(qconta, qcontb);
}

TEST_F(CScreeningContainerTest, CopyConstructor)
{
    auto mlih = vlxmol::getTestLiH();
    
    auto mbas = vlxbas::getTestBasisForLiH();
    
    CGtoBlock bgtos(mlih, mbas, 0);
    
    CGtoBlock kgtos(mlih, mbas, 1);
    
    CGtoPairsBlock bpairs(bgtos, kgtos, 1.0e-13);
    
    CGtoPairsBlock kpairs(bgtos, 1.0e-13);
    
    CMemBlock<double> bqvals({1.0, 2.0, 3.0, 4.0});
    
    CMemBlock<double> kqvals({0.0, 1.0, 2.0, 7.0, 7.0, 2.0});
    
    CVecMemBlock<double> bqvec{bqvals, kqvals};
    
    CVecMemBlock<double> kqvec{kqvals, kqvals};
    
    CGtoPairsContainer bppcont({bpairs, kpairs});
    
    CGtoPairsContainer kppcont({kpairs, kpairs});
    
    CScreeningContainer qscra(bqvec, kqvec, bppcont, kppcont, ericut::qq,
                              1.0e-13);
    
    CScreeningContainer qscrb(qscra);
    
    ASSERT_EQ(qscra, qscrb);
}

TEST_F(CScreeningContainerTest, MoveConstructor)
{
    auto mlih = vlxmol::getTestLiH();
    
    auto mbas = vlxbas::getTestBasisForLiH();
    
    CGtoBlock bgtos(mlih, mbas, 0);
    
    CGtoBlock kgtos(mlih, mbas, 1);
    
    CGtoPairsBlock bpairs(bgtos, kgtos, 1.0e-13);
    
    CGtoPairsBlock kpairs(bgtos, 1.0e-13);
    
    CMemBlock<double> bqvals({1.0, 2.0, 3.0, 4.0});
    
    CMemBlock<double> kqvals({0.0, 1.0, 2.0, 7.0, 7.0, 2.0});
    
    CVecMemBlock<double> bqvec{bqvals, kqvals};
    
    CVecMemBlock<double> kqvec{kqvals, kqvals};
    
    CGtoPairsContainer bppcont({bpairs, kpairs});
    
    CGtoPairsContainer kppcont({kpairs, kpairs});
    
    CScreeningContainer qscra(bqvec, kqvec, bppcont, kppcont, ericut::qq,
                              1.0e-13);
    
    CScreeningContainer qscrb(CScreeningContainer(bqvec, kqvec, bppcont,
                                                  kppcont, ericut::qq,
                                                  1.0e-13));
    
    ASSERT_EQ(qscra, qscrb);
}

TEST_F(CScreeningContainerTest, CopyAssignment)
{
    auto mlih = vlxmol::getTestLiH();
    
    auto mbas = vlxbas::getTestBasisForLiH();
    
    CGtoBlock bgtos(mlih, mbas, 0);
    
    CGtoBlock kgtos(mlih, mbas, 1);
    
    CGtoPairsBlock bpairs(bgtos, kgtos, 1.0e-13);
    
    CGtoPairsBlock kpairs(bgtos, 1.0e-13);
    
    CMemBlock<double> bqvals({1.0, 2.0, 3.0, 4.0});
    
    CMemBlock<double> kqvals({0.0, 1.0, 2.0, 7.0, 7.0, 2.0});
    
    CVecMemBlock<double> bqvec{bqvals, kqvals};
    
    CVecMemBlock<double> kqvec{kqvals, kqvals};
    
    CGtoPairsContainer bppcont({bpairs, kpairs});
    
    CGtoPairsContainer kppcont({kpairs, kpairs});
    
    CScreeningContainer qscra(bqvec, kqvec, bppcont, kppcont, ericut::qq,
                              1.0e-13);
    
    CScreeningContainer qscrb = qscra;
    
    ASSERT_EQ(qscra, qscrb);
}

TEST_F(CScreeningContainerTest, MoveAssignment)
{
    auto mlih = vlxmol::getTestLiH();
    
    auto mbas = vlxbas::getTestBasisForLiH();
    
    CGtoBlock bgtos(mlih, mbas, 0);
    
    CGtoBlock kgtos(mlih, mbas, 1);
    
    CGtoPairsBlock bpairs(bgtos, kgtos, 1.0e-13);
    
    CGtoPairsBlock kpairs(bgtos, 1.0e-13);
    
    CMemBlock<double> bqvals({1.0, 2.0, 3.0, 4.0});
    
    CMemBlock<double> kqvals({0.0, 1.0, 2.0, 7.0, 7.0, 2.0});
    
    CVecMemBlock<double> bqvec{bqvals, kqvals};
    
    CVecMemBlock<double> kqvec{kqvals, kqvals};
    
    CGtoPairsContainer bppcont({bpairs, kpairs});
    
    CGtoPairsContainer kppcont({kpairs, kpairs});
    
    CScreeningContainer qscra(bqvec, kqvec, bppcont, kppcont, ericut::qq,
                              1.0e-13);
    
    CScreeningContainer qscrb = CScreeningContainer(bqvec, kqvec, bppcont,
                                                    kppcont, ericut::qq,
                                                    1.0e-13);
    
    ASSERT_EQ(qscra, qscrb);
}

TEST_F(CScreeningContainerTest, IsEmpty)
{
    auto mlih = vlxmol::getTestLiH();
    
    auto mbas = vlxbas::getTestBasisForLiH();
    
    CGtoBlock bgtos(mlih, mbas, 0);
    
    CGtoBlock kgtos(mlih, mbas, 1);
    
    CGtoPairsBlock bpairs(bgtos, kgtos, 1.0e-13);
    
    CGtoPairsBlock kpairs(bgtos, 1.0e-13);
    
    CMemBlock<double> bqvals({1.0, 2.0, 3.0, 4.0});
    
    CMemBlock<double> kqvals({0.0, 1.0, 2.0, 7.0, 7.0, 2.0});
    
    CVecMemBlock<double> bqvec{bqvals, kqvals};
    
    CVecMemBlock<double> kqvec{kqvals, kqvals};
    
    CGtoPairsContainer bppcont({bpairs, kpairs});
    
    CGtoPairsContainer kppcont({kpairs, kpairs});
    
    CScreeningContainer qscra(bqvec, kqvec, bppcont, kppcont, ericut::qq,
                              1.0e-13);
    
    CScreeningContainer qscrb;
    
    ASSERT_FALSE(qscra.isEmpty());
    
    ASSERT_TRUE(qscrb.isEmpty());
}

TEST_F(CScreeningContainerTest, GetNumberOfScreeners)
{
    auto mlih = vlxmol::getTestLiH();
    
    auto mbas = vlxbas::getTestBasisForLiH();
    
    CGtoBlock bgtos(mlih, mbas, 0);
    
    CGtoBlock kgtos(mlih, mbas, 1);
    
    CGtoPairsBlock bpairs(bgtos, kgtos, 1.0e-13);
    
    CGtoPairsBlock kpairs(bgtos, 1.0e-13);
    
    CMemBlock<double> bqvals({1.0, 2.0, 3.0, 4.0});
    
    CMemBlock<double> kqvals({0.0, 1.0, 2.0, 7.0, 7.0, 2.0});
    
    CVecMemBlock<double> bqvec{bqvals, kqvals};
    
    CVecMemBlock<double> kqvec{kqvals, kqvals};
    
    CGtoPairsContainer bppcont({bpairs, kpairs});
    
    CGtoPairsContainer kppcont({kpairs, kpairs});
    
    CScreeningContainer qscra(bqvec, kqvec, bppcont, kppcont, ericut::qq,
                              1.0e-13);
    
    CScreeningContainer qscrb(bqvec, bqvec, bppcont, bppcont, ericut::qq,
                              1.0e-13);
    
    ASSERT_EQ(qscra.getNumberOfScreeners(), 4);
    
    ASSERT_EQ(qscrb.getNumberOfScreeners(), 3);
}

TEST_F(CScreeningContainerTest, GetScreener)
{
    auto mlih = vlxmol::getTestLiH();
    
    auto mbas = vlxbas::getTestBasisForLiH();
    
    CGtoBlock bgtos(mlih, mbas, 0);
    
    CGtoBlock kgtos(mlih, mbas, 1);
    
    CGtoPairsBlock bpairs(bgtos, kgtos, 1.0e-13);
    
    CGtoPairsBlock kpairs(bgtos, 1.0e-13);
    
    CMemBlock<double> bqvals({1.0, 2.0, 3.0, 4.0});
    
    CMemBlock<double> kqvals({0.0, 1.0, 2.0, 7.0, 7.0, 2.0});
    
    CVecMemBlock<double> bqvec{bqvals, kqvals};
    
    CVecMemBlock<double> kqvec{kqvals, kqvals};
    
    CGtoPairsContainer bppcont({bpairs, kpairs});
    
    CGtoPairsContainer kppcont({kpairs, kpairs});
    
    CScreeningContainer qscra(bqvec, kqvec, bppcont, kppcont, ericut::qq,
                              1.0e-13);
    
    CCauchySchwarzScreener qq0dat(bqvals, kqvals, bpairs, kpairs, ericut::qq,
                                  1.0e-13);
    
    CCauchySchwarzScreener qq2dat(kqvals, kqvals, kpairs, kpairs, ericut::qq,
                                  1.0e-13);
    
    ASSERT_EQ(qscra.getScreener(0), qq2dat);
    
    ASSERT_EQ(qscra.getScreener(1), qq2dat);
    
    ASSERT_EQ(qscra.getScreener(2), qq0dat);
    
    ASSERT_EQ(qscra.getScreener(3), qq0dat);
    
    ASSERT_EQ(qscra.getScreener(4), CCauchySchwarzScreener());
}
