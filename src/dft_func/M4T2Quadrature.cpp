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

#include "M4T2Quadrature.hpp"

#include <cmath>

#include "MathConst.hpp"
#include "MathFunc.hpp"

CM4T2Quadrature::CM4T2Quadrature(const int nRadialPoints, const int idElemental)

    : _nRadialPoints(nRadialPoints)

    , _idElemental(idElemental)
{
}

CM4T2Quadrature::~CM4T2Quadrature()
{
}

auto
CM4T2Quadrature::generate() const -> CDenseMatrix
{
    if (_nRadialPoints > 0)
    {
        CDenseMatrix qpoints(2, _nRadialPoints);

        // set up aliases to grid points data

        auto rcoords = qpoints.row(0);

        auto rweights = qpoints.row(1);

        // compute Chebyshev quadrature of second kind

        mathfunc::quadChebyshevOfKindTwo(rcoords, rweights, _nRadialPoints);

        // apply M4 mapping

        auto fxi = _getXiFactor() / std::log(2.0);

#pragma omp simd
        for (int i = 0; i < _nRadialPoints; i++)
        {
            // set up parameters

            auto x = rcoords[i];

            auto mx = 1.0 - x;

            auto px = 1.0 + x;

            // compute auxilary functions

            auto fone = std::pow(px, 0.6);

            auto ftwo = std::log(2.0 / mx);

            // transform coordinates using M4 mapping

            rcoords[i] = fxi * fone * ftwo;

            // accumulate tranformation factors into weights

            rweights[i] *= 1.0 / std::sqrt(1.0 - x * x);

            rweights[i] *= fxi * (fone / mx + 0.6 * ftwo * std::pow(px, -0.4));
        }

        return qpoints;
    }

    return CDenseMatrix();
}

auto
CM4T2Quadrature::_getXiFactor() const -> double
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
