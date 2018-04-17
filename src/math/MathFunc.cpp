//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include "MathFunc.hpp"

namespace mathfunc { // mathfunc namespace

    void zero(double* vector, const int32_t nElements)
    {
        #pragma omp simd aligned (vector: VLX_ALIGN)
        for (int32_t i = 0; i < nElements; i++) vector[i] = 0.0;
    }

    void set_to(double* vector, const double value, const int32_t nElements)
    {
        #pragma omp simd aligned (vector: VLX_ALIGN)
        for (int32_t i = 0; i < nElements; i++) vector[i] = value;
    }

    double sum(const double* vector, const int32_t nElements)
    {
        double fsum = 0.0;

        #pragma omp simd aligned (vector: VLX_ALIGN)
        for (int32_t i = 0; i < nElements; i++) fsum += vector[i];

        return fsum;
    }

    int32_t sum(const int32_t* vector, const int32_t nElements)
    {
        int32_t isum = 0;

        #pragma omp simd aligned (vector: VLX_ALIGN)
        for (int32_t i = 0; i < nElements; i++) isum += vector[i];

        return isum;
    }

    void normalize(double* vector, const int32_t nElements)
    {
        auto factor = 1.0 / mathfunc::sum(vector, nElements);

        #pragma omp simd aligned (vector: VLX_ALIGN)
        for (int32_t i = 0; i < nElements; i++) vector[i] *= factor;
    }

    void indexes(int32_t* aVector, const int32_t* bVector,
                 const int32_t nElements)
    {
        int32_t index = 0;

        for (int32_t i = 0; i < nElements; i++)
        {
            aVector[i] = index;

            index += bVector[i];
        }
    }

} // mathfunc namespace


