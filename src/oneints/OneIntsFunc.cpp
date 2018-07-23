//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include "OneIntsFunc.hpp"


namespace intsfunc { // intsfunc namespace

    void
    compXiAndZeta(      double* factorsXi,
                        double* factorsZeta,
                  const double  braExponent,
                  const double* ketExponents,
                  const int32_t nElements)
    {
        #pragma omp simd aligned(factorsXi, factorsZeta, ketExponents: VLX_ALIGN)
        for (int32_t i = 0; i < nElements; i++)
        {
            factorsXi[i] = 1.0 / (braExponent + ketExponents[i]);
            
            factorsZeta[i] = braExponent * ketExponents[i] * factorsXi[i];
        }
    }
    
    void
    compDistancesPA(      double* xDistancesPA,
                          double* yDistancesPA,
                          double* zDistancesPA,
                    const double* ketExponents,
                    const double* factorsXi,
                    const double* xDistancesAB,
                    const double* yDistancesAB,
                    const double* zDistancesAB,
                    const int32_t nElements)
    {
        #pragma omp simd aligned(xDistancesPA, yDistancesPA, zDistancesPA,\
                                 ketExponents, factorsXi, xDistancesAB, \
                                 yDistancesAB, zDistancesAB: VLX_ALIGN)
        for (int32_t i = 0; i < nElements; i++)
        {
            double fact = - ketExponents[i] * factorsXi[i];
            
            xDistancesPA[i] = fact * xDistancesAB[i];
            
            yDistancesPA[i] = fact * yDistancesAB[i];
            
            zDistancesPA[i] = fact * zDistancesAB[i];
        }
    }
    
    void
    compDistancesPB(      double* xDistancesPB,
                          double* yDistancesPB,
                          double* zDistancesPB,
                    const double  braExponent,
                    const double* factorsXi,
                    const double* xDistancesAB,
                    const double* yDistancesAB,
                    const double* zDistancesAB,
                    const int32_t nElements)
    {
       #pragma omp simd aligned(xDistancesPB, yDistancesPB, zDistancesPB,\
                                factorsXi, xDistancesAB, yDistancesAB,\
                                zDistancesAB: VLX_ALIGN)
        for (int32_t i = 0; i < nElements; i++)
        {
            double fact = braExponent * factorsXi[i];
            
            xDistancesPB[i] = fact * xDistancesAB[i];
            
            yDistancesPB[i] = fact * yDistancesAB[i];
            
            zDistancesPB[i] = fact * zDistancesAB[i];
        }
    }
    
} // intsfunc namespace
