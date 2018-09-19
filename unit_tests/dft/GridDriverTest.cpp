//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "GridDriverTest.hpp"

#include <cmath>

#include "GridDriver.hpp"
#include "OutputStream.hpp"
#include "MoleculeSetter.hpp"
#include "MathConst.hpp"
#include "MpiFunc.hpp"

TEST_F(CGridDriverTest, DefaultConstructor)
{
    CGridDriver gdrv(mpi::master(), mpi::nodes(MPI_COMM_WORLD), execmode::cpu,
                     MPI_COMM_WORLD);
    
    gdrv.setLevel(6, MPI_COMM_WORLD);
    
    COutputStream ost(std::string("dummy.out"));
    
    auto mlih = vlxmol::getMoleculeLiH();
    
    auto mgrid = gdrv.generate(mlih, ost, MPI_COMM_WORLD);
    
    auto npnt = mgrid.getNumberOfGridPoints();
    
    EXPECT_EQ(741140, npnt);
    
    auto rx = mgrid.getCoordinatesX();
    
    auto ry = mgrid.getCoordinatesY();
    
    auto rz = mgrid.getCoordinatesZ();
    
    auto w = mgrid.getWeights();
    
    double fa  = 0.0;
    
    double fb  = 0.0;
    
    double fab = 0.0;
    
    for (int32_t i = 0; i < npnt; i++)
    {
        auto r2a = rx[i] * rx[i] + ry[i] * ry[i] + rz[i] * rz[i];
        
        fa += w[i] * std::exp(-2.3 * r2a);
        
        auto r2b = rx[i] * rx[i] + ry[i] * ry[i]
        
                 + (rz[i] - 1.20) * (rz[i] - 1.20);
        
       fb += w[i] * std::exp(-0.5 * r2b);
        
       fab += w[i] * std::exp(-2.3 * r2a) * std::exp(-0.5 * r2b);
    }
    
    EXPECT_NEAR(fa, 1.5963681525241624, 1.0e-10);
    
    EXPECT_NEAR(fb, 15.749609945385632, 1.0e-10);
    
    EXPECT_NEAR(fab, 0.65786017622805693, 1.0e-10);
}
