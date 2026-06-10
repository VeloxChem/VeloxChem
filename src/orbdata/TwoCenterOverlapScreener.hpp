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

#ifndef TwoCenterOverlapScreener_hpp
#define TwoCenterOverlapScreener_hpp

#include "BasisFunction.hpp"
#include "Point.hpp"

/// @brief Class CTwoCenterOverlapScreener is a light screening predicate for a
/// pair of basis functions placed on two centers. It is constructed from the
/// two basis functions: the exponent dependent data (which is the same for all
/// atom pairs sharing these functions) is precomputed once from the most
/// diffuse (smallest exponent) primitive on each side, and the center-only
/// call operator keeps an atom pair when the estimated overlap is greater than
/// or equal to a threshold.
class CTwoCenterOverlapScreener
{
   public:
    /// @brief The constructor with bra and ket basis functions and threshold.
    /// Precomputes the screening factor, decay exponent and total angular
    /// momentum from the smallest exponent (and its coefficient) on each side.
    /// @param bra_function The basis function on bra side.
    /// @param ket_function The basis function on ket side.
    /// @param threshold The screening threshold; pairs with an estimated
    /// overlap below this value are discarded.
    CTwoCenterOverlapScreener(const CBasisFunction &bra_function, const CBasisFunction &ket_function, const double threshold = 1.0e-14);

    /// @brief Evaluates the screening predicate for the basis function pair on
    /// two centers. The estimate is
    ///   est = (pi/s)^{3/2} * |c_a| * |c_b| * max(1, R^{la+lb}) * exp(-a_min*b_min/s * R^2),
    /// with s = a_min + b_min and R the distance between the centers.
    /// @param bra_center The Cartesian coordinates of the bra atom.
    /// @param ket_center The Cartesian coordinates of the ket atom.
    /// @return True if the estimated overlap is greater than or equal to the
    /// threshold, False otherwise.
    auto operator()(const TPoint<double> &bra_center, const TPoint<double> &ket_center) const -> bool;

    /// @brief Gets the screening threshold.
    /// @return The screening threshold.
    auto threshold() const -> double;

   private:
    /// @brief The precomputed screening factor (pi/s)^{3/2} * |c_a| * |c_b|.
    double _factor;

    /// @brief The precomputed decay exponent a_min * b_min / (a_min + b_min).
    double _exponent;

    /// @brief The total angular momentum l_a + l_b of the basis function pair.
    int _total_angular_momentum;

    /// @brief The screening threshold.
    double _threshold;
};

#endif /* TwoCenterOverlapScreener_hpp */
