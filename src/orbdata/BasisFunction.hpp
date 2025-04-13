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

#ifndef BasisFunction_hpp
#define BasisFunction_hpp

#include <cstddef>
#include <cstdint>
#include <vector>

/**
 Class CBasisFunction stores data about single contracted GTO and
 provides set of methods for handling of GTO data.
 */
class CBasisFunction
{
    /**
     The vector of exponents of primitive Gaussian functions.
     */
    std::vector<double> _exponents;

    /**
     The vector of normalization factors of primitive Gaussian functions.
     */
    std::vector<double> _norms;

    /**
     The angular momentum of basis function.
     */
    int64_t _angmom;

    /**
     Rescales normalization factors to match normalization of spherical (l,0)
     component of basis function.
     */
    auto _rescale() -> void;

    /**
     Computes overlap between two primitive Gaussian functions.

     @param i the index of first primitve Gaussain function.
     @param j the index of second primitve Gaussain function.
     @return the overlap between two primitive Gaussian functions.
     */
    auto _overlap(const size_t i, const size_t j) const -> double;

   public:
    /**
     Creates an empty basis function object.
     */
    CBasisFunction() = default;

    /**
     Creates a basis function object.

     @param exponents the vector of exponents of primitive Gaussian functions.
     @param norms the vector of normalization factors of primitive Gaussian functions.
     @param angmom the angular momentum of basis function.
     */
    CBasisFunction(const std::vector<double>& exponents, const std::vector<double>& norms, const int64_t angmom);

    /**
     Sets exponents of primittive Gaussian functions with specific vector of
     exponents.

     @param exponents the vector of exponents.
     */
    auto setExponents(const std::vector<double>& exponents) -> void;

    /**
     Sets normalization factors of primitive Gaussian functions with specific
     vector of normalization factors.

     @param norms the vector of normalization factors.
     */
    auto setNormalizationFactors(const std::vector<double>& norms) -> void;

    /**
     Set angular momentum of basis function.

     @param angmom the angular momentum.
     */
    auto setAngularMomentum(const int64_t angmom) -> void;

    /**
     Adds primittive Gaussian function to basis function.

     @param exponent the exponent of primitive Gaussian function.
     @param norm the normalization factor of primitive Gaussian function.
     */
    auto add(const double exponent, const double norm) -> void;

    /**
     Normalizes basis function.
     */
    auto normalize() -> void;

    /**
     Gets vector of exponents of primitive Gaussian functions.

     @return the vector of exponents.
     */
    auto getExponents() const -> std::vector<double>;

    /**
     Gets vector of normalization factors of primitive Gaussian functions.

     @return the vector of normalization factors.
     */
    auto getNormalizationFactors() const -> std::vector<double>;

    /**
     Gets angular momentum of basis function.

     @return the angular momentum.
     */
    auto getAngularMomentum() const -> int64_t;

    /**
     Gets number of primitive Gaussian functions in basis function.

     @return the number of primitive Gaussian functions.
     */
    auto getNumberOfPrimitiveFunctions() const -> int64_t;
};

#endif /* BasisFunction_hpp */
