//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include "AtomBasisSetter.hpp"

namespace vlxbas { // vlxbas namespace
    
CAtomBasis
getAtomBasisEmpty()
{
    CAtomBasis bas;

    bas.setMaxAngularMomentum(-1);

    bas.setIdElemental(-1);

    return bas;
}

CAtomBasis
getAtomBasisForH()
{
    CAtomBasis bas;

    bas.setIdElemental(1);

    bas.addBasisFunction(CBasisFunction({1.301070100000e+01, 1.962257200000e+00,
                                         4.445379600000e-01},
                                        {1.968215800000e-02, 1.379652400000e-01,
                                         4.783193500000e-01}, 0));

    bas.addBasisFunction(CBasisFunction({1.219496200000e-01},
                                        {1.000000000000e+00}, 0));

    bas.addBasisFunction(CBasisFunction({8.000000000000e-01},
                                        {1.000000000000e+00}, 1));

    return bas;
}

CAtomBasis
getAtomBasisForLi()
{
    CAtomBasis bas;

    bas.setIdElemental(3);

    bas.addBasisFunction(CBasisFunction({2.662778551600e+02, 4.006978344700e+01,
                                         9.055994438900e+00, 2.450300905100e+00,
                                         7.220957185500e-01},
                                        {6.492015032500e-03, 4.774786321500e-02,
                                         2.026879611100e-01, 4.860657481700e-01,
                                         4.362697795500e-01}, 0));

    bas.addBasisFunction(CBasisFunction({1.450000000000e+00,
                                         3.000000000000e-01},
                                        {2.586000000000e-01,
                                         1.000000000000e+00}, 1));

    bas.addBasisFunction(CBasisFunction({5.281088472100e-02},
                                        {1.000000000000e+00}, 0));

    bas.addBasisFunction(CBasisFunction({2.096094879800e-02},
                                        {1.000000000000e+00}, 0));

    bas.addBasisFunction(CBasisFunction({8.200000000000e-02},
                                        {1.000000000000e+00}, 1));

    return bas;
}
  
CAtomBasis
getTestBasisForH(const int32_t angularMomentum)
{
    CAtomBasis bas;
    
    bas.setIdElemental(1);
    
    CBasisFunction bfone({3.0, 2.0}, {0.5, 0.5}, angularMomentum);
    
    bfone.normalize();
    
    CBasisFunction bftwo({0.8}, {1.0}, angularMomentum);
    
    bftwo.normalize();
    
    bas.addBasisFunction(bfone);
    
    bas.addBasisFunction(bftwo);
    
    return bas; 
}
    
CAtomBasis
getTestBasisForLi(const int32_t angularMomentum)
{
    CAtomBasis bas;
    
    bas.setIdElemental(3);
        
    CBasisFunction bfone({2.8, 1.5}, {0.5, 0.5}, angularMomentum);
        
    bfone.normalize();
        
    CBasisFunction bftwo({1.2}, {1.0}, angularMomentum);
        
    bftwo.normalize();
        
    bas.addBasisFunction(bfone);
        
    bas.addBasisFunction(bftwo);
        
    return bas;
}
    
} // vlxbas namespace
