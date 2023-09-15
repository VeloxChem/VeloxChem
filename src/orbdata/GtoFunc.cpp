#include "GtoFunc.hpp"

namespace gtofunc {  // gtofunc namespace

auto
makeGtoBlocks(const CMolecularBasis& basis, const CMolecule& molecule) -> std::vector<CGtoBlock>
{
    if (const auto mang = basis.getMaxAngularMomentum(); mang >= 0)
    {
        std::vector<CGtoBlock> gto_blocks;

        for (int64_t i = 0; i <= mang; i++)
        {
            for (const auto npgtos : basis.getContractionDepths(i))
            {
                gto_blocks.push_back(CGtoBlock(basis, molecule, i, npgtos));
            }
        }

        return gto_blocks;
    }
    else
    {
        return std::vector<CGtoBlock>();
    }
}

auto
makeGtoBlocks(const CMolecularBasis& basis, const CMolecule& molecule, const std::vector<int64_t>& atoms) -> std::vector<CGtoBlock>
{
    if (const auto mang = basis.getMaxAngularMomentum(atoms); mang >= 0)
    {
        std::vector<CGtoBlock> gto_blocks;

        for (int64_t i = 0; i <= mang; i++)
        {
            for (const auto npgtos : basis.getContractionDepths(atoms, i))
            {
                gto_blocks.push_back(CGtoBlock(basis, molecule, atoms, i, npgtos));
            }
        }

        return gto_blocks;
    }
    else
    {
        return std::vector<CGtoBlock>();
    }
}

auto
getNumberOfAtomicOrbitals(const std::vector<CGtoBlock>& gto_blocks) -> int64_t
{
    int64_t naos = 0;

    for (const auto& gto_block : gto_blocks)
    {
        const auto ncgtos = gto_block.getNumberOfBasisFunctions();

        const auto ang = gto_block.getAngularMomentum();

        naos += ncgtos * (ang * 2 + 1);
    }

    return naos;
}

}  // namespace gtofunc
