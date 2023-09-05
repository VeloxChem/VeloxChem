//
//                              VELOXCHEM
//         ----------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2023 by VeloxChem developers. All rights reserved.
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

#ifndef Log3Quadrature_hpp
#define Log3Quadrature_hpp

#include <cstdint>

#include "DenseMatrix.hpp"

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
    auto _getXiFactor() const -> double;

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
    auto generate() const -> CDenseMatrix;
};

#endif /* Log3Quadrature_hpp */
