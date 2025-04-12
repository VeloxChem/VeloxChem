//
//                                   VELOXCHEM
//              ----------------------------------------------------
//                          An Electronic Structure Code
//
//  SPDX-License-Identifier: BSD-3-Clause
//
//  Copyright 2018-2025 VeloxChem developers
//
//  Redistribution and use in source and binary forms, with or without modification,
//  are permitted provided that the following conditions are met:
//
//  1. Redistributions of source code must retain the above copyright notice, this
//     list of conditions and the following disclaimer.
//  2. Redistributions in binary form must reproduce the above copyright notice,
//     this list of conditions and the following disclaimer in the documentation
//     and/or other materials provided with the distribution.
//  3. Neither the name of the copyright holder nor the names of its contributors
//     may be used to endorse or promote products derived from this software without
//     specific prior written permission.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
//  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
//  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
//  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
//  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
//  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
//  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
//  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
//  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
//  OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#ifndef M4T2Quadrature_hpp
#define M4T2Quadrature_hpp

#include <cstdint>

#include "DenseMatrix.hpp"

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
    int64_t _nRadialPoints;

    /**
     The identifier of chemical element.
     */
    int64_t _idElemental;

    /**
     Gets Xi factor of radial quadrature for specific chemical element.

     @return the Xi factor of radial quadrature.
     */
    auto _getXiFactor() const -> double;

   public:
    /**
     Creates a M4T2 quadrature object.

     @param nRadialPoints the number of radial points.
     @param idElemental the identifier of chemical element.
     */
    CM4T2Quadrature(const int64_t nRadialPoints, const int64_t idElemental);

    /**
     Destroys a M4T2 quadrature object.
     */
    ~CM4T2Quadrature();

    /**
     Generates quadrature points for radial M4T2 quadrature.

     @return the quadrature points (coordinates, weights).
     */
    auto generate() const -> CDenseMatrix;
};

#endif /* M4T2Quadrature_hpp */
