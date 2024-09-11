#include "T4CScreener.hpp"

#include <algorithm>
#include <ranges>

#include "GtoPairBlockFunc.hpp"
#include "OpenMPFunc.hpp"
#include "T4CDiagonalDistributor.hpp"
#include "DiagonalElectronRepulsionFunc.hpp"

auto
CT4CScreener::partition(const CMolecularBasis& basis,
                        const CMolecule&       molecule,
                        const std::string&     label) -> void
{
    const auto gto_pair_blocks = gtofunc::make_gto_pair_blocks(basis, molecule);
    
    if (const auto nblocks = gto_pair_blocks.size(); nblocks > 0)
    {
        // set up max. integral values
        
        std::vector<std::vector<double>> max_values;
        
        max_values.reserve(nblocks);

        std::ranges::transform(gto_pair_blocks, std::back_inserter(max_values), [](const auto& gpblock) {
            return std::vector<double>(gpblock.number_of_contracted_pairs(), 0.0);
        });
        
        // prepare pointers for OMP parallel region

        auto ptr_gto_pair_blocks = &gto_pair_blocks;

        auto ptr_max_values = max_values.data();

        // execute OMP tasks with static scheduling

        #pragma omp parallel shared(ptr_gto_pair_blocks, ptr_max_values, label)
        {
            #pragma omp single nowait
            {
                const auto work_tasks = omp::make_diag_work_group(*ptr_gto_pair_blocks);
                
                std::ranges::for_each(std::views::reverse(work_tasks), [&](const auto& task) {
                    const auto index = task[0];
                    const auto gto_range = std::pair<size_t, size_t>{task[1], task[2]};
                    #pragma omp task firstprivate(index, gto_range)
                    {
                        const auto gto_pair_block = ptr_gto_pair_blocks->at(index);
                        CT4CDiagonalDistributor distributor(ptr_max_values[index].data());
                        if (label == "eri")
                        {
                            erifunc::diag_compute<CT4CDiagonalDistributor>(distributor, gto_pair_block, gto_range);
                        }
                    }
                });
            }
        }
        
        // apply Cauchy–Schwarz partitioning
        
        _gto_pair_blocks.clear();
        
        _gto_pair_blocks.reserve(nblocks);
        
        std::ranges::for_each(std::views::iota(size_t{0}, nblocks), [&](const auto index) {
            _gto_pair_blocks.push_back(CBlockedGtoPairBlock(gto_pair_blocks[index], max_values[index]));
        } );
    }
}

auto
CT4CScreener::gto_pair_blocks() const -> std::vector<CBlockedGtoPairBlock>
{
    return _gto_pair_blocks;
}
