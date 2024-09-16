#include "GtoFunc.hpp"

#include <algorithm>
#include <ranges>

namespace gtofunc {  // gtofunc namespace

auto
make_gto_blocks(const CMolecularBasis &basis, const CMolecule &molecule) -> std::vector<CGtoBlock>
{
    std::vector<CGtoBlock> gto_blocks;

    std::ranges::for_each(std::views::iota(0, basis.max_angular_momentum() + 1), [&](const int i) {
        std::ranges::for_each(basis.contraction_depths(i), [&](const int j) { gto_blocks.push_back(CGtoBlock(basis, molecule, i, j)); });
    });

    return gto_blocks;
}

auto
make_gto_blocks(const CMolecularBasis &basis, const CMolecule &molecule, const std::vector<int> &atoms) -> std::vector<CGtoBlock>
{
    std::vector<CGtoBlock> gto_blocks;

    std::ranges::for_each(std::views::iota(0, basis.max_angular_momentum(atoms) + 1), [&](const int i) {
        std::ranges::for_each(basis.contraction_depths(atoms, i),
                              [&](const int j) { gto_blocks.push_back(CGtoBlock(basis, molecule, atoms, i, j)); });
    });

    return gto_blocks;
}

auto
getNumberOfAtomicOrbitals(const std::vector<CGtoBlock>& gto_blocks) -> int
{
    int naos = 0;

    for (const auto& gto_block : gto_blocks)
    {
        const auto ncgtos = gto_block.number_of_basis_functions();

        const auto ang = gto_block.angular_momentum();

        naos += ncgtos * (ang * 2 + 1);
    }

    return naos;
}

}  // namespace gtofunc
