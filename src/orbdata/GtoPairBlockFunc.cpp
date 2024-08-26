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

}  // namespace gtofunc
