//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#ifndef Log3Quadrature_hpp
#define Log3Quadrature_hpp

#include <cstdint>

#include "MemBlock2D.hpp"

/**
 Class CLog3Quadrature class generates M4T2 aka "Log3" quadrature.
 Reference: O. Treutler and R. Ahlrichs, J. Chem. Phys. 102, 346 (1995).
 
 @author Z. Rinkevicius
 */
class CLog3Quadrature
{
    /**
     The number of radial points.
     */
    int32_t _nRadialPoints;

    /**
     The identifier of chemical element.
     */
    int32_t _idElemental;

    /**
     Gets Xi factor of radial quadrature for specific chemical element.

     @return the Xi factor of radial quadrature.
     */
    double _getXiFactor() const;

public:

    /**
     Creates a Log3 quadrature object.
     
     @param nRadialPoints the number of radial points.
     @param idElemental the identifier of chemical element.
     */
    CLog3Quadrature(const int32_t nRadialPoints, const int32_t idElemental);

    /**
     Destroys a Log3 quadrature object.
     */
    ~CLog3Quadrature();

    /**
     Generates quadrature points for radial Log3 quadrature.
     
     @return the quadrature points (coordinates, weights).
     */
    CMemBlock2D<double> generate() const;
};

#endif /* Log3Quadrature_hpp */
