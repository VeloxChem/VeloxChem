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

#include "BlockedGtoPairBlock.hpp"

#include <cmath>
#include <ranges>

#include "MathFunc.hpp"

CBlockedGtoPairBlock::CBlockedGtoPairBlock(const std::array<CGtoPairBlock, 16>& gto_pair_blocks)

    : _gto_pair_blocks(gto_pair_blocks)
{
    
}

CBlockedGtoPairBlock::CBlockedGtoPairBlock(const std::vector<CGtoPairBlock> &gto_pair_blocks, const std::vector<int> &block_indices)

    : _gto_pair_blocks{}
{
    std::ranges::for_each(block_indices, [&](const int i) { _gto_pair_blocks[i] = gto_pair_blocks[i]; });
}

CBlockedGtoPairBlock::CBlockedGtoPairBlock(const CGtoPairBlock &gto_pair_block, const std::vector<double> &integrals)

    : _gto_pair_blocks(std::array<CGtoPairBlock, 16>())
{
    // set up block map from integrals

    std::vector<int> block_map;

    block_map.reserve(integrals.size());

    std::ranges::transform(integrals, std::back_inserter(block_map), [](auto tint) -> int {
        const auto fact = std::sqrt(tint);
        if (fact > 1.0e-15)
        {
            return (fact < 1.0) ? int(-std::floor(std::log10(fact))) : 0;
        }
        else
        {
            return -1;
        }
    });

    // set up dimensions of GTOs pair block

    const auto ncpairs = gto_pair_block.number_of_contracted_pairs();

    const auto nppairs = gto_pair_block.number_of_primitive_pairs();

    // set up data of GTOs pair block

    const auto gp_ang_pair = gto_pair_block.angular_momentums();

    const auto gp_bra_coords = gto_pair_block.bra_coordinates();

    const auto gp_ket_coords = gto_pair_block.ket_coordinates();

    const auto gp_norms = gto_pair_block.normalization_factors();

    const auto gp_ovls = gto_pair_block.overlap_factors();

    const auto gp_bra_exponents = gto_pair_block.bra_exponents();

    const auto gp_ket_exponents = gto_pair_block.ket_exponents();

    const auto gp_bra_orb_indices = gto_pair_block.bra_orbital_indices();

    const auto gp_ket_orb_indices = gto_pair_block.ket_orbital_indices();

    const auto gp_bra_atm_indices = gto_pair_block.bra_atomic_indices();

    const auto gp_ket_atm_indices = gto_pair_block.ket_atomic_indices();

    // partition GTOs pair block
    // TODO: replace with algorithms based loops

    for (int i = 0; i < 16; i++)
    {
        if (const auto icpairs = mathfunc::count_elements_by_values(block_map, i); icpairs > 0)
        {
            auto igp_bra_coords = std::vector<TPoint<double>>(icpairs, TPoint<double>({0.0, 0.0, 0.0}));

            auto igp_ket_coords = std::vector<TPoint<double>>(icpairs, TPoint<double>({0.0, 0.0, 0.0}));

            auto igp_bra_orb_indices = std::vector<size_t>(icpairs + 1, 0);

            auto igp_ket_orb_indices = std::vector<size_t>(icpairs + 1, 0);

            auto igp_bra_atm_indices = std::vector<int>(icpairs, 0);

            auto igp_ket_atm_indices = std::vector<int>(icpairs, 0);

            auto igp_bra_exponents = std::vector<double>(icpairs * nppairs, 0.0);

            auto igp_ket_exponents = std::vector<double>(icpairs * nppairs, 0.0);

            auto igp_norms = std::vector<double>(icpairs * nppairs, 0.0);

            auto igp_ovls = std::vector<double>(icpairs * nppairs, 0.0);

            igp_bra_orb_indices[0] = gp_bra_orb_indices[0];

            igp_ket_orb_indices[0] = gp_ket_orb_indices[0];

            int idx = 0;

            for (int j = 0; j < ncpairs; j++)
            {
                if (block_map[j] == i)
                {
                    igp_bra_coords[idx] = gp_bra_coords[j];

                    igp_ket_coords[idx] = gp_ket_coords[j];

                    igp_bra_orb_indices[idx + 1] = gp_bra_orb_indices[j + 1];

                    igp_ket_orb_indices[idx + 1] = gp_ket_orb_indices[j + 1];

                    igp_bra_atm_indices[idx] = gp_bra_atm_indices[j];

                    igp_ket_atm_indices[idx] = gp_ket_atm_indices[j];

                    for (int k = 0; k < nppairs; k++)
                    {
                        igp_norms[k * icpairs + idx] = gp_norms[k * ncpairs + j];

                        igp_ovls[k * icpairs + idx] = gp_ovls[k * ncpairs + j];

                        igp_bra_exponents[k * icpairs + idx] = gp_bra_exponents[k * ncpairs + j];

                        igp_ket_exponents[k * icpairs + idx] = gp_ket_exponents[k * ncpairs + j];
                    }

                    idx++;
                }
            }

            _gto_pair_blocks[i] = CGtoPairBlock(igp_bra_coords,
                                                igp_ket_coords,
                                                igp_bra_exponents,
                                                igp_ket_exponents,
                                                igp_norms,
                                                igp_ovls,
                                                igp_bra_orb_indices,
                                                igp_ket_orb_indices,
                                                igp_bra_atm_indices,
                                                igp_ket_atm_indices,
                                                gp_ang_pair,
                                                nppairs);
        }
    }
}

CBlockedGtoPairBlock::CBlockedGtoPairBlock(const CBlockedGtoPairBlock &other)

    : _gto_pair_blocks(other._gto_pair_blocks)
{
}

CBlockedGtoPairBlock::CBlockedGtoPairBlock(CBlockedGtoPairBlock &&other) noexcept

    : _gto_pair_blocks(std::move(other._gto_pair_blocks))
{
}

auto
CBlockedGtoPairBlock::operator=(const CBlockedGtoPairBlock &other) -> CBlockedGtoPairBlock &
{
    _gto_pair_blocks = other._gto_pair_blocks;

    return *this;
}

auto
CBlockedGtoPairBlock::operator=(CBlockedGtoPairBlock &&other) noexcept -> CBlockedGtoPairBlock &
{
    if (this != &other)
    {
        _gto_pair_blocks = std::move(other._gto_pair_blocks);
    }

    return *this;
}

auto
CBlockedGtoPairBlock::operator==(const CBlockedGtoPairBlock &other) const -> bool
{
    return std::ranges::equal(_gto_pair_blocks, other._gto_pair_blocks);
}

auto
CBlockedGtoPairBlock::operator!=(const CBlockedGtoPairBlock &other) const -> bool
{
    return !(*this == other);
}

auto
CBlockedGtoPairBlock::gto_pair_block(const int index) const -> CGtoPairBlock
{
    return _gto_pair_blocks.at(index);
}

auto
CBlockedGtoPairBlock::is_empty_gto_pair_block(const int index) const -> bool
{
    return (_gto_pair_blocks[index].number_of_contracted_pairs() == 0);
}

auto
CBlockedGtoPairBlock::gto_pair_blocks() const -> std::array<CGtoPairBlock, 16>
{
    return _gto_pair_blocks;
}
