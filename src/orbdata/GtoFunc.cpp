#include "GtoFunc.hpp"

#include <algorithm>
#include <ranges>

namespace gtofunc {  // gtofunc namespace

auto
make_gto_blocks(const CMolecularBasis &basis, const CMolecule &molecule) -> std::vector<CGtoBlock>
{
    std::vector<CGtoBlock> gto_blocks;

    std::ranges::for_each(std::views::iota(0, basis.max_angular_momentum() + 1), [&](const int i) {
        std::ranges::for_each(basis.contraction_depths(i),
                              [&](const int j) { gto_blocks.push_back(CGtoBlock(basis, molecule, i, j)); });
    });

    return gto_blocks;
}

auto
make_gto_blocks(const CMolecularBasis &basis, const CMolecule &molecule, const std::vector<int> &atoms)
    -> std::vector<CGtoBlock>
{
    std::vector<CGtoBlock> gto_blocks;

    std::ranges::for_each(std::views::iota(0, basis.max_angular_momentum(atoms) + 1), [&](const int i) {
        std::ranges::for_each(basis.contraction_depths(atoms, i),
                              [&](const int j) { gto_blocks.push_back(CGtoBlock(basis, molecule, atoms, i, j)); });
    });

    return gto_blocks;
}

}  // namespace gtofunc
