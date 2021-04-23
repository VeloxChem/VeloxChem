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

#include "TwoIntsDistributorTest.hpp"

#include "TwoIntsDistributor.hpp"
#include "MoleculeSetter.hpp"
#include "MolecularBasisSetter.hpp"
#include "AODensityMatrix.hpp"
#include "AOFockMatrix.hpp"

TEST_F(CTwoIntsDistributionTest, DefaultConstructor)
{
    CTwoIntsDistribution dista;
    
    CTwoIntsDistribution distb(nullptr, 0, 0, dist2e::batch);
    
    ASSERT_EQ(dista, distb);
}

TEST_F(CTwoIntsDistributionTest, CopyConstructor)
{
    double val = 1.0;
    
    CTwoIntsDistribution dista(&val, 10, 20, 2, dist2e::qvalues);
    
    CTwoIntsDistribution distb(dista);
    
    ASSERT_EQ(dista, distb);
}

TEST_F(CTwoIntsDistributionTest, CopyAssignment)
{
    double val = 1.0;
    
    CTwoIntsDistribution dista(&val, 10, 20, 2, dist2e::qvalues);
    
    CTwoIntsDistribution distb = dista;
    
    ASSERT_EQ(dista, distb);
}

TEST_F(CTwoIntsDistributionTest, DistributeWithBatchNoSym)
{
    auto mlih = vlxmol::getTestLiH();
    
    auto mbas = vlxbas::getTestBasisForLiH();
    
    CGtoBlock bgtos(mlih, mbas, 0);
    
    CGtoBlock kgtos(mlih, mbas, 1);
    
    CGtoPairsBlock bpairs(bgtos, kgtos, 1.0e-13);
    
    CGtoPairsBlock kpairs(bgtos, 1.0e-13);
   
    CMemBlock<double> fints(3);
    
    fints.zero();
   
    CTwoIntsDistribution dista(fints.data(), 1, 1, 0, dist2e::batch);
    
    CMemBlock2D<double> refvals({1.0, 3.0, 4.0, 5.0, 7.0, 3.0}, 2, 3);
    
    dista.distribute(refvals, bpairs.pick(2), kpairs.pick(3), false, 0, 0);
    
    ASSERT_EQ(fints, CMemBlock<double>({1.0, 4.0, 7.0}));
    
    CMemBlock<double> tints(9);
    
    tints.zero();
    
    CTwoIntsDistribution distb(tints.data(), 1, 1, 0, dist2e::batch);
    
    distb.distribute(refvals, bpairs.pick(2), kpairs.pick(3), false, 0, 1);
    
    ASSERT_EQ(tints, CMemBlock<double>({0.0, 0.0, 0.0, 1.0, 4.0, 7.0,
                                        0.0, 0.0, 0.0}));
    
    distb.distribute(refvals, bpairs.pick(2), kpairs.pick(3), false, 0, 2);
    
    ASSERT_EQ(tints, CMemBlock<double>({0.0, 0.0, 0.0, 1.0, 4.0, 7.0,
                                        1.0, 4.0, 7.0}));
    
    distb.distribute(refvals, bpairs.pick(2), kpairs.pick(3), false, 0, 0);
    
    ASSERT_EQ(tints, CMemBlock<double>({1.0, 4.0, 7.0, 1.0, 4.0, 7.0,
                                        1.0, 4.0, 7.0}));
    
    CMemBlock<double> gints(12);
    
    gints.zero();
    
    CTwoIntsDistribution distc(gints.data(), 1, 2, 0, dist2e::batch);
    
    distc.distribute(refvals, bpairs.pick(2), kpairs.pick(3), false, 0, 1);
    
    ASSERT_EQ(gints, CMemBlock<double>({0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                        1.0, 4.0, 7.0, 0.0, 0.0, 0.0}));
    
    distc.distribute(refvals, bpairs.pick(2), kpairs.pick(3), false, 0, 0);
    
    ASSERT_EQ(gints, CMemBlock<double>({1.0, 4.0, 7.0, 0.0, 0.0, 0.0,
                                        1.0, 4.0, 7.0, 0.0, 0.0, 0.0}));
}

TEST_F(CTwoIntsDistributionTest, DistributeWithBatchSym)
{
    auto mlih = vlxmol::getTestLiH();
    
    auto mbas = vlxbas::getTestBasisForLiH();
    
    CGtoBlock bgtos(mlih, mbas, 0);
    
    CGtoPairsBlock bpairs(bgtos, 1.0e-13);
    
    CMemBlock<double> fints(6);
    
    fints.zero(); 
    
    CTwoIntsDistribution dista(fints.data(), 1, 1, 0, dist2e::batch);
    
    CMemBlock2D<double> ref0vals({2.0}, 1, 1);
    
    dista.distribute(ref0vals, bpairs.pick(2), bpairs.pick(1), true, 0, 0);
    
    ASSERT_EQ(fints, CMemBlock<double>({2.0, 0.0, 0.0, 0.0, 0.0, 0.0}));
    
    CMemBlock2D<double> ref1vals({4.0, 3.0}, 2, 1);
    
    dista.distribute(ref1vals, bpairs.pick(2), bpairs.pick(1), true, 0, 1);
    
    ASSERT_EQ(fints, CMemBlock<double>({2.0, 4.0, 3.0, 0.0, 0.0, 0.0}));
    
    CMemBlock2D<double> ref2vals({1.0, 2.0, 7.0}, 3, 1);
    
    dista.distribute(ref2vals, bpairs.pick(2), bpairs.pick(1), true, 0, 2);
    
    ASSERT_EQ(fints, CMemBlock<double>({2.0, 4.0, 3.0, 1.0, 2.0, 7.0}));
}

TEST_F(CTwoIntsDistributionTest, DistributeWithQValues)
{
    auto mlih = vlxmol::getTestLiH();
    
    auto mbas = vlxbas::getTestBasisForLiH();
    
    CGtoBlock bgtos(mlih, mbas, 0);
    
    CGtoBlock kgtos(mlih, mbas, 1);
    
    CGtoPairsBlock bpairs(bgtos, kgtos, 1.0e-13);
    
    CGtoPairsBlock kpairs(bgtos, 1.0e-13);
    
    CMemBlock<double> fints(3);
    
    fints.zero();
    
    CTwoIntsDistribution dista(fints.data(), 1, 1, 2, dist2e::qvalues);
    
    CMemBlock2D<double> ref0vals({1.0, 3.0, 4.0}, 1, 3);
    
    dista.distribute(ref0vals, bpairs.pick(0), kpairs.pick(1), false, 0, 5);
    
    ASSERT_EQ(fints, CMemBlock<double>({0.0, 0.0, 2.0}));
    
    CTwoIntsDistribution distb(fints.data(), 1, 1, 0, dist2e::qvalues);
    
    CMemBlock2D<double> ref1vals({1.0, 9.0, 4.0}, 1, 3);
    
    distb.distribute(ref1vals, bpairs.pick(0), kpairs.pick(1), false, 0, 5);
    
    ASSERT_EQ(fints, CMemBlock<double>({3.0, 0.0, 2.0}));
    
    CTwoIntsDistribution distc(fints.data(), 1, 2, 1, dist2e::qvalues);
    
    CMemBlock2D<double> ref2vals({16.0, 9.0, 4.0}, 1, 3);
    
    distc.distribute(ref2vals, bpairs.pick(0), kpairs.pick(1), true, 0, 2);
    
    ASSERT_EQ(fints, CMemBlock<double>({3.0, 4.0, 2.0}));
}
