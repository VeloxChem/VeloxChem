//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include "Log3Quadrature.hpp"

#include <cmath>

#include "MathConst.hpp"
#include "MathFunc.hpp"

CLog3Quadrature::CLog3Quadrature(const int32_t nRadialPoints,
                                 const int32_t idElemental)

    : _nRadialPoints(nRadialPoints)

    , _idElemental(idElemental)
{

}

CLog3Quadrature::~CLog3Quadrature()
{

}

CMemBlock2D<double> CLog3Quadrature::generate() const
{
    if (_nRadialPoints > 0)
    {
        CMemBlock2D<double> qpoints(_nRadialPoints, 2);
        
        // set up aliases to grid points data
        
        auto rcoords  = qpoints.data(0);
        
        auto rweights = qpoints.data(1);
        
        // compute Chebyshev quadrature of second kind

        mathfunc::quadChebyshevOfKindTwo(rcoords, rweights, _nRadialPoints);

        // apply M4 mapping

        auto fxi = _getXiFactor() / std::log(2.0);

        #pragma omp simd aligned(rcoords:VLX_ALIGN, rweights:VLX_ALIGN)
        for (int32_t i = 0; i < _nRadialPoints; i++)
        {
            // set up parameters

            auto x  = rcoords[i];

            auto mx = 1.0 - x;

            auto px = 1.0 + x;

            // compute auxilary functions 

            auto fone = std::pow(px, 0.6) ;

            auto ftwo = std::log(2.0 / mx);

            // transform coordinates using M4 mapping

            rcoords[i] = fxi * fone * ftwo;

            // accumulate tranformation factors into weights

            rweights[i] *= 1.0 / std::sqrt(1.0 - x * x);

            rweights[i] *= fxi * (fone / mx + 0.6 * ftwo * std::pow(px, -0.4));
        }

        return qpoints;
    }

    return CMemBlock2D<double>();
}

double CLog3Quadrature::_getXiFactor() const
{
    // H atom

    if (_idElemental == 1) return 0.8;

    // He atom

    if (_idElemental == 2) return 0.9;

    // Li atom

    if (_idElemental == 3) return 1.8;

    // Be atom

    if (_idElemental == 4) return 1.4;

    // B atom

    if (_idElemental == 5) return 1.3;

    // C atom

    if (_idElemental == 6) return 1.1;

    // N-Ne atoms

    if ((_idElemental > 6) && (_idElemental < 11)) return 0.9;

    // Na atom

    if (_idElemental == 11) return 1.4;

    // Mg, Al atoms

    if ((_idElemental == 12) || (_idElemental == 13)) return 1.3;

    // Si atom

    if (_idElemental == 14) return 1.2;

    // P atom

    if (_idElemental == 15) return 1.1;

    // K atom

    if (_idElemental == 19) return 1.5;

    // Ca atom

    if (_idElemental == 20) return 1.4;

    // Sc atom

    if (_idElemental == 21) return 1.3;

    // Ti-Co atoms

    if ((_idElemental > 21) && (_idElemental < 28)) return 1.2;

    // Ni-Ga atoms

    if ((_idElemental > 27) && (_idElemental < 32)) return 1.1;

    // As-Kr atoms

    if ((_idElemental > 32) && (_idElemental < 37)) return 0.9;

    // default value

    return 1.0;
}
