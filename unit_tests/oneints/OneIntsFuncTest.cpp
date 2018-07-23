//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include "OneIntsFuncTest.hpp"

#include "OneIntsFunc.hpp"
#include "CheckFunctions.hpp"

TEST_F(COneIntsFuncTest, CompXiAndZeta)
{
    double kexp[5] __attribute__ ((aligned(VLX_ALIGN))) = {5.0, 8.0, 1.0, 3.0, 4.0};
    
    double fxi[5] __attribute__ ((aligned(VLX_ALIGN)));
    
    double fzeta[5] __attribute__ ((aligned(VLX_ALIGN)));
    
    intsfunc::compXiAndZeta(fxi, fzeta, 2.0, kexp, 5);
    
    double refxi[5] __attribute__ ((aligned(VLX_ALIGN))) = {1.0/7.0, 0.1, 1.0/3.0,
                                                            0.2, 1.0/6.0};
    
    vlxtest::compare(fxi, refxi, 5);
    
    double refzeta[5] __attribute__ ((aligned(VLX_ALIGN))) = {10.0/7.0, 1.6, 2.0/3.0,
                                                              1.2, 8.0/6.0};
    
    vlxtest::compare(fzeta, refzeta, 5);
}

TEST_F(COneIntsFuncTest, CompDistancesPA)
{
    double kexp[5] __attribute__ ((aligned(VLX_ALIGN))) = {5.0, 8.0, 1.0, 3.0, 4.0};
    
    double fxi[5] __attribute__ ((aligned(VLX_ALIGN))) = {0.2, 0.3, 0.4, 0.5, 0.6};
    
    double abx[5] __attribute__ ((aligned(VLX_ALIGN))) = {1.0, 2.0, 3.0, 4.0, 5.0};
    
    double aby[5] __attribute__ ((aligned(VLX_ALIGN))) = {2.0, 3.0, 4.0, 1.0, 0.0};
    
    double abz[5] __attribute__ ((aligned(VLX_ALIGN))) = {1.0, 1.0, 3.0, 2.0, 1.0};
    
    double pax[5] __attribute__ ((aligned(VLX_ALIGN)));
    
    double pay[5] __attribute__ ((aligned(VLX_ALIGN)));
    
    double paz[5] __attribute__ ((aligned(VLX_ALIGN)));
    
    intsfunc::compDistancesPA(pax, pay, paz, kexp, fxi, abx, aby, abz, 5);
    
    double refpax[5] __attribute__ ((aligned(VLX_ALIGN))) = {-1.0, -4.8, -1.2,
                                                             -6.0, -12.0};
    
    vlxtest::compare(pax, refpax, 5);
    
    double refpay[5] __attribute__ ((aligned(VLX_ALIGN))) = {-2.0, -7.2, -1.6,
                                                             -1.5, 0.0};
    
    vlxtest::compare(pay, refpay, 5);
    
    double refpaz[5] __attribute__ ((aligned(VLX_ALIGN))) = {-1.0, -2.4, -1.2,
                                                             -3.0, -2.4};
    
    vlxtest::compare(paz, refpaz, 5);
}

TEST_F(COneIntsFuncTest, CompDistancesPB)
{
    double fxi[5] __attribute__ ((aligned(VLX_ALIGN))) = {0.2, 0.3, 0.4, 0.5, 0.6};
    
    double abx[5] __attribute__ ((aligned(VLX_ALIGN))) = {1.0, 2.0, 3.0, 4.0, 5.0};
    
    double aby[5] __attribute__ ((aligned(VLX_ALIGN))) = {2.0, 3.0, 4.0, 1.0, 0.0};
    
    double abz[5] __attribute__ ((aligned(VLX_ALIGN))) = {1.0, 1.0, 3.0, 2.0, 1.0};
    
    double pbx[5] __attribute__ ((aligned(VLX_ALIGN)));
    
    double pby[5] __attribute__ ((aligned(VLX_ALIGN)));
    
    double pbz[5] __attribute__ ((aligned(VLX_ALIGN)));
    
    intsfunc::compDistancesPB(pbx, pby, pbz, 2.0, fxi, abx, aby, abz, 5);
    
    double refpbx[5] __attribute__ ((aligned(VLX_ALIGN))) = {0.4, 1.2, 2.4, 4.0, 6.0};
    
    vlxtest::compare(pbx, refpbx, 5);
    
    double refpby[5] __attribute__ ((aligned(VLX_ALIGN))) = {0.8, 1.8, 3.2, 1.0, 0.0};
    
    vlxtest::compare(pby, refpby, 5);
    
    double refpbz[5] __attribute__ ((aligned(VLX_ALIGN))) = {0.4, 0.6, 2.4, 2.0, 1.2};
    
    vlxtest::compare(pbz, refpbz, 5);
}

