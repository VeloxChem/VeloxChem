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

#include "GtoPairBlock.hpp"

#include "MatrixIndex.hpp"

CGtoPairBlock::CGtoPairBlock(const std::vector<TPairOfPoints3D>& coordinates,
                             const std::vector<TPoint2D>&        exponents,
                             const std::vector<TPoint2D>&        norms,
                             const std::vector<T2Index>&         orb_indexes,
                             const std::vector<T2Index>&         atm_indexes,
                             const T2Index&                      angmoms,
                             const int64_t                       nppairs)

    : _coordinates(coordinates)

    , _exponents(exponents)

    , _norms(norms)

    , _orb_indexes(orb_indexes)

    , _atm_indexes(atm_indexes)

    , _angmoms(angmoms)

    , _nppairs(nppairs)
{
}

CGtoPairBlock::CGtoPairBlock(const CGtoBlock& gto_block)

    : _coordinates(std::vector<TPairOfPoints3D>())

    , _exponents(std::vector<TPoint2D>())

    , _norms(std::vector<TPoint2D>())

    , _orb_indexes(std::vector<T2Index>())

    , _atm_indexes(std::vector<T2Index>())

    , _angmoms(T2Index({-1, -1}))

    , _nppairs(-1)
{
    // fetch GTOs block data

    const auto coords = gto_block.getCoordinates();

    const auto gexps = gto_block.getExponents();

    const auto gnorms = gto_block.getNormalizationFactors();

    const auto orbidx = gto_block.getOrbitalIndexes();

    const auto atmidx = gto_block.getAtomicIndexes();

    // set up dimensions

    auto ncgtos = gto_block.getNumberOfBasisFunctions();

    auto npgtos = gto_block.getNumberOfPrimitives();

    // allocate memory

    const auto cdim = ncgtos * (ncgtos + 1) / 2;

    _coordinates = std::vector<TPairOfPoints3D>(cdim, TPairOfPoints3D());

    _orb_indexes = std::vector<T2Index>(cdim + 1, T2Index());

    _atm_indexes = std::vector<T2Index>(cdim, T2Index());

    const auto pdim = cdim * npgtos * npgtos;

    _exponents = std::vector<TPoint2D>(pdim, TPoint2D({0.0, 0.0}));

    _norms = std::vector<TPoint2D>(pdim, TPoint2D({0.0, 0.0}));

    // set up GTO pairs data

    const auto angmom = gto_block.getAngularMomentum();

    _orb_indexes[0] = T2Index({orbidx[0], orbidx[0]});

    _angmoms = {angmom, angmom};

    _nppairs = npgtos * npgtos;

    for (int64_t i = 0; i < ncgtos; i++)
    {
        for (int64_t j = i; j < ncgtos; j++)
        {
            const auto ij_idx = mathfunc::uplo_index(i, j);

            _coordinates[ij_idx] = TPairOfPoints3D({coords[i], coords[j]});

            _orb_indexes[ij_idx + 1] = T2Index({orbidx[i + 1], orbidx[j + 1]});

            _atm_indexes[ij_idx] = T2Index({atmidx[i], atmidx[j]});

            for (int64_t k = 0; k < npgtos; k++)
            {
                const auto koff = k * ncgtos + i;

                for (int64_t l = 0; l < npgtos; l++)
                {
                    const auto loff = l * ncgtos + j;

                    const auto ijoff = (k * npgtos + l) * cdim + ij_idx;

                    _exponents[ijoff] = TPoint2D({gexps[koff], gexps[loff]});

                    _norms[ijoff] = TPoint2D({gnorms[koff], gnorms[loff]});
                }
            }
        }
    }
}

CGtoPairBlock::CGtoPairBlock(const CGtoBlock& bra_gto_block, const CGtoBlock& ket_gto_block)

    : _coordinates(std::vector<TPairOfPoints3D>())

    , _exponents(std::vector<TPoint2D>())

    , _norms(std::vector<TPoint2D>())

    , _orb_indexes(std::vector<T2Index>())

    , _atm_indexes(std::vector<T2Index>())

    , _angmoms(T2Index({-1, -1}))

    , _nppairs(-1)
{
    // fetch GTOs block data on bra side

    const auto bcoords = bra_gto_block.getCoordinates();

    const auto bexps = bra_gto_block.getExponents();

    const auto bnorms = bra_gto_block.getNormalizationFactors();

    const auto borbidx = bra_gto_block.getOrbitalIndexes();

    const auto batmidx = bra_gto_block.getAtomicIndexes();

    // fetch GTOs block data on ket side

    const auto kcoords = ket_gto_block.getCoordinates();

    const auto kexps = ket_gto_block.getExponents();

    const auto knorms = ket_gto_block.getNormalizationFactors();

    const auto korbidx = ket_gto_block.getOrbitalIndexes();

    const auto katmidx = ket_gto_block.getAtomicIndexes();

    // set up dimensions

    auto bcgtos = bra_gto_block.getNumberOfBasisFunctions();

    auto bpgtos = bra_gto_block.getNumberOfPrimitives();

    auto kcgtos = ket_gto_block.getNumberOfBasisFunctions();

    auto kpgtos = ket_gto_block.getNumberOfPrimitives();

    // reserve memory

    const auto cdim = bcgtos * kcgtos;

    _coordinates = std::vector<TPairOfPoints3D>(cdim, TPairOfPoints3D());

    _orb_indexes = std::vector<T2Index>(cdim + 1, T2Index());

    _atm_indexes = std::vector<T2Index>(cdim, T2Index());

    const auto pdim = cdim * bpgtos * kpgtos;

    _exponents = std::vector<TPoint2D>(pdim, TPoint2D({0.0, 0.0}));

    _norms = std::vector<TPoint2D>(pdim, TPoint2D({0.0, 0.0}));

    // set up GTO pairs data

    const auto bangmom = bra_gto_block.getAngularMomentum();

    const auto kangmom = ket_gto_block.getAngularMomentum();

    _orb_indexes[0] = T2Index({borbidx[0], korbidx[0]});

    _angmoms = {bangmom, kangmom};

    _nppairs = bpgtos * kpgtos;

    for (int64_t i = 0; i < bcgtos; i++)
    {
        for (int64_t j = 0; j < kcgtos; j++)
        {
            const auto ij_idx = i * kcgtos + j;

            _coordinates[ij_idx] = TPairOfPoints3D({bcoords[i], kcoords[j]});

            _orb_indexes[ij_idx + 1] = T2Index({borbidx[i + 1], korbidx[j + 1]});

            _atm_indexes[ij_idx] = T2Index({batmidx[i], katmidx[j]});

            for (int64_t k = 0; k < bpgtos; k++)
            {
                const auto koff = k * bcgtos + i;

                for (int64_t l = 0; l < kpgtos; l++)
                {
                    const auto loff = l * kcgtos + j;

                    const auto ijoff = (k * kpgtos + l) * cdim + ij_idx;

                    _exponents[ijoff] = TPoint2D({bexps[koff], kexps[loff]});

                    _norms[ijoff] = TPoint2D({bnorms[koff], knorms[loff]});
                }
            }
        }
    }
}

auto
CGtoPairBlock::getCoordinates() const -> std::vector<TPairOfPoints3D>
{
    return _coordinates;
}

auto
CGtoPairBlock::getExponents() const -> std::vector<TPoint2D>
{
    return _exponents;
}

auto
CGtoPairBlock::getNormalizationFactors() const -> std::vector<TPoint2D>
{
    return _norms;
}

auto
CGtoPairBlock::getOrbitalIndexes() const -> std::vector<T2Index>
{
    return _orb_indexes;
}

auto
CGtoPairBlock::getAtomicIndexes() const -> std::vector<T2Index>
{
    return _atm_indexes;
}

auto
CGtoPairBlock::getAngularMomentums() const -> T2Index
{
    return _angmoms;
}

auto
CGtoPairBlock::getNumberOfPrimitivePairs() const -> int64_t
{
    return _nppairs;
}

auto
CGtoPairBlock::getNumberOfContractedPairs() const -> int64_t
{
    return static_cast<int64_t>(_coordinates.size());
}
