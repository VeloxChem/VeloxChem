//
//                              VELOXCHEM
//         ----------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2024 by VeloxChem developers. All rights reserved.
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
