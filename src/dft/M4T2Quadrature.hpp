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

#ifndef M4T2Quadrature_hpp
#define M4T2Quadrature_hpp

#include <cstdint>

#include "MemBlock2D.hpp"

/**
 Class CM4T2Quadrature class generates M4T2 quadrature.
 Reference: O. Treutler and R. Ahlrichs, J. Chem. Phys. 102, 346 (1995).

 @author Z. Rinkevicius
 */
class CM4T2Quadrature
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
     Creates a M4T2 quadrature object.

     @param nRadialPoints the number of radial points.
     @param idElemental the identifier of chemical element.
     */
    CM4T2Quadrature(const int32_t nRadialPoints, const int32_t idElemental);

    /**
     Destroys a M4T2 quadrature object.
     */
    ~CM4T2Quadrature();

    /**
     Generates quadrature points for radial M4T2 quadrature.

     @return the quadrature points (coordinates, weights).
     */
    CMemBlock2D<double> generate() const;
};

#endif /* M4T2Quadrature_hpp */
