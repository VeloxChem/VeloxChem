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

#ifndef GtoPairBlock_hpp
#define GtoPairBlock_hpp

#include <array>
#include <cstdint>
#include <string>
#include <vector>

#include "GtoBlock.hpp"
#include "Point.hpp"
#include "T2Index.hpp"

/**
 Class CGtoPairBlock stores data about contarcted GTO pairs block and provides set of methods
 for manipulating with contarcted GTO pairs block.

 @author Z. Rinkevicius
 */
class CGtoPairBlock
{
    /**
     The vector of Cartesian coordinates of contracted GTO pairs.
     */
    std::vector<TPairOfPoints3D> _coordinates;

    /**
     The vector of exponents of primitive GTOs pairs.
     */
    std::vector<TPoint2D> _exponents;

    /**
     The vector of normalization factors of primitive GTO pairs.
     */
    std::vector<TPoint2D> _norms;

    /**
     The vector of AO indexes of contracted GTO pairs.
     */
    std::vector<T2Index> _orb_indexes;

    /**
     The vector of atomic indexes of contracted GTO pairs.
     */
    std::vector<T2Index> _atm_indexes;

    /**
     The angular momentums of contracted GTO pair.
     */
    T2Index _angmoms;

    /**
     The number of primitive GTOs in contracted GTO pair.
     */
    int64_t _nppairs;

   public:
    /**
     Creates an empty contarcted GTO pairs block.
     */
    CGtoPairBlock() = default;

    /**
     Creates a contarcted GTO pairs block.

     @param coordinates the vector of basis functions pair coordinates.
     @param exponents the vector of exponents of primitive GTO pairs.
     @param norms the vector of normalization factors of primitive GTO pairs.
     @param orb_indexes the vector of  AO indexes.
     @param atm_indexes the vector of  atomic indexes.
     @param angmoms the angular momentums of GTO pair.
     @param nppairs the number of primitive GTO pairs in contracted GTO pairs.
     */
    CGtoPairBlock(const std::vector<TPairOfPoints3D>& coordinates,
                  const std::vector<TPoint2D>&        exponents,
                  const std::vector<TPoint2D>&        norms,
                  const std::vector<T2Index>&         orb_indexes,
                  const std::vector<T2Index>&         atm_indexes,
                  const T2Index&                      angmoms,
                  const int64_t                       nppairs);

    /**
     Creates a contarcted GTO pairs block.

     @param gto_block the GTOs block.
     */
    CGtoPairBlock(const CGtoBlock& gto_block);

    /**
     Creates a contarcted GTO pairs block.

     @param bra_gto_block the GTOs block on bra side.
     @param ket_gto_block the GTOs block on ket side.
     */
    CGtoPairBlock(const CGtoBlock& bra_gto_block, const CGtoBlock& ket_gto_block);

    /**
     Gets vector of GTO pair coordinates.

     @return the vector of GTO pair coordinates.
     */
    auto getCoordinates() const -> std::vector<TPairOfPoints3D>;

    /**
     Gets vector of GTO pair exponents.

     @return the vector of GTO pair exponents.
     */
    auto getExponents() const -> std::vector<TPoint2D>;

    /**
     Gets vector of GTO pair normalization factors.

     @return the vector of GTO pair normalization factors.
     */
    auto getNormalizationFactors() const -> std::vector<TPoint2D>;

    /**
     Gets vector of orbital indexes of contracted GTO pairs.

     @return the vector of orbital indexes of GTO pairs.
     */
    auto getOrbitalIndexes() const -> std::vector<T2Index>;

    /**
     Gets vector of atomic indexes of contracted GTO pairs.

     @return the vector of atomic indexes of GTO pairs.
     */
    auto getAtomicIndexes() const -> std::vector<T2Index>;

    /**
     Gets angular momentums of GTO pair.

     @return the angular momentums of GTO pair.
     */
    auto getAngularMomentums() const -> T2Index;

    /**
     Gets number of primitive GTO pairs in contracted GTO pair.

     @return the number of primtive GTO pairs in contracted GTO pair.
     */
    auto getNumberOfPrimitivePairs() const -> int64_t;

    /**
     Gets number of GTO pairs in contracted GTO pairs block.

     @return the number of GTO pairs in contracted GTO pairs block.
     */
    auto getNumberOfContractedPairs() const -> int64_t;
};

#endif /* GtoPairBlock_hpp */
