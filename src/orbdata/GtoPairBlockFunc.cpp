#include "GtoPairBlockFunc.hpp"

#include <ranges>

#include "CustomViews.hpp"
#include "GtoFunc.hpp"

namespace gtofunc {  // gtofunc namespace

auto
make_gto_pair_blocks(const CMolecularBasis& basis, const CMolecule& molecule) -> std::vector<CGtoPairBlock>
{
    const auto gblocks = gtofunc::make_gto_blocks(basis, molecule);

    if (const auto nblocks = gblocks.size(); nblocks > 0)
    {
        std::vector<CGtoPairBlock> gpblocks;

        gpblocks.reserve(nblocks * (nblocks + 1) / 2);

        std::ranges::for_each(views::triangular(nblocks), [&](const auto& index) {
            const auto [i, j] = index;
            if (i == j)
            {
                gpblocks.push_back(CGtoPairBlock(gblocks[i]));
            }
            else
            {
                gpblocks.push_back(CGtoPairBlock(gblocks[i], gblocks[j]));
            }
        });

        return gpblocks;
    }
    else
    {
        return std::vector<CGtoPairBlock>();
    }
}

auto
make_gto_pair_blocks(const std::vector<CGtoBlock>& bra_gto_blocks,
                     const std::vector<CGtoBlock>& ket_gto_blocks) -> std::vector<CGtoPairBlock>
{
    const auto nbra_blocks = bra_gto_blocks.size();
    
    const auto nket_blocks = ket_gto_blocks.size();
    
    if  ((nbra_blocks > 0) && (nket_blocks > 0))
    {
        std::vector<CGtoPairBlock> gto_pair_blocks;
            
        std::ranges::for_each(views::rectangular(nbra_blocks, nket_blocks), [&](const auto& index) {
            const auto [i, j] = index;
            gto_pair_blocks.push_back(CGtoPairBlock(bra_gto_blocks[i], ket_gto_blocks[j]));
        });
            
        return gto_pair_blocks;
    }
    else
    {
        return std::vector<CGtoPairBlock>();
    }
}

}  // namespace gtofunc
