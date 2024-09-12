#include "OpenMPFunc.hpp"

#include <algorithm>
#include <ranges>
#include <utility>

#include "BatchFunc.hpp"
#include "CustomViews.hpp"
#include "SimdArray.hpp"

namespace omp {  // omp namespace

auto
make_work_tasks(const std::vector<CGtoBlock>& gto_blocks) -> std::vector<std::array<size_t, 6>>
{
    auto tasks = std::vector<std::array<size_t, 6>>{};

    if (const auto nblocks = gto_blocks.size(); nblocks > 0)
    {
        const auto nthreads = omp::get_number_of_threads();

        std::ranges::for_each(views::triangular(nblocks), [&](const auto& index) {
            const auto [i, j]    = index;
            const auto bra_size  = static_cast<size_t>(gto_blocks[i].number_of_basis_functions());
            const auto ket_size  = static_cast<size_t>(gto_blocks[j].number_of_basis_functions());
            auto       bra_bsize = bra_size / nthreads;
            auto       ket_bsize = ket_size / nthreads;
            if (bra_bsize < simd::width<double>()) bra_bsize = simd::width<double>();
            if (ket_bsize < simd::width<double>()) ket_bsize = simd::width<double>();
            const auto bra_blocks = batch::number_of_batches(bra_size, bra_bsize);
            const auto ket_blocks = batch::number_of_batches(ket_size, ket_bsize);
            if (i == j)
            {
                std::ranges::for_each(views::triangular(bra_blocks), [&](const auto& bkpair) {
                    const auto [k, l]   = bkpair;
                    const auto bindices = batch::batch_range(k, bra_size, bra_bsize);
                    const auto kindices = batch::batch_range(l, bra_size, bra_bsize);
                    tasks.push_back({index.first, index.second, bindices.first, bindices.second, kindices.first, kindices.second});
                });
            }
            else
            {
                std::ranges::for_each(views::rectangular(bra_blocks, ket_blocks), [&](const auto& bkpair) {
                    const auto [k, l]   = bkpair;
                    const auto bindices = batch::batch_range(k, bra_size, bra_bsize);
                    const auto kindices = batch::batch_range(l, ket_size, ket_bsize);
                    tasks.push_back({index.first, index.second, bindices.first, bindices.second, kindices.first, kindices.second});
                });
            }
        });
    }

    return tasks;
}

auto
make_work_tasks(const std::vector<CGtoBlock>& bra_gto_blocks, const std::vector<CGtoBlock>& ket_gto_blocks) -> std::vector<std::array<size_t, 6>>
{
    auto tasks = std::vector<std::array<size_t, 6>>{};

    if (const auto bra_nblocks = bra_gto_blocks.size(), ket_nblocks = ket_gto_blocks.size(); (bra_nblocks > 0) && (ket_nblocks > 0))
    {
        const auto nthreads = omp::get_number_of_threads();

        std::ranges::for_each(views::rectangular(bra_nblocks, ket_nblocks), [&](const auto& index) {
            const auto [i, j]    = index;
            const auto bra_size  = static_cast<size_t>(bra_gto_blocks[i].number_of_basis_functions());
            const auto ket_size  = static_cast<size_t>(ket_gto_blocks[j].number_of_basis_functions());
            auto       bra_bsize = bra_size / nthreads;
            auto       ket_bsize = ket_size / nthreads;
            if (bra_bsize < simd::width<double>()) bra_bsize = simd::width<double>();
            if (ket_bsize < simd::width<double>()) ket_bsize = simd::width<double>();
            const auto bra_blocks = batch::number_of_batches(bra_size, bra_bsize);
            const auto ket_blocks = batch::number_of_batches(ket_size, ket_bsize);
            std::ranges::for_each(views::rectangular(bra_blocks, ket_blocks), [&](const auto& bkpair) {
                const auto [k, l]   = bkpair;
                const auto bindices = batch::batch_range(k, bra_size, bra_bsize);
                const auto kindices = batch::batch_range(l, ket_size, ket_bsize);
                tasks.push_back({index.first, index.second, bindices.first, bindices.second, kindices.first, kindices.second});
            });
        });
    }

    return tasks;
}

auto
make_diag_work_group(const std::vector<CGtoPairBlock>& gto_pair_blocks) -> std::vector<std::array<size_t, 3>>
{
    auto wtasks = std::vector<std::array<size_t, 3>>();

    if (const auto nblocks = gto_pair_blocks.size(); nblocks > 0)
    {
        const auto nthreads = omp::get_number_of_threads();
        
        for (size_t i = 0; i < nblocks; i++)
        {
            const auto gp_size = gto_pair_blocks[i].number_of_contracted_pairs();
            
            auto bsize = gp_size / nthreads;
            
            if (bsize < simd::width<double>()) bsize = simd::width<double>();
            
            const auto bblocks = batch::number_of_batches(gp_size, bsize);
            
            for (size_t j = 0; j < bblocks; j++)
            {
                const auto bindices = batch::batch_range(j, gp_size, bsize);
                
                wtasks.push_back(std::array<size_t, 3>({i, bindices.first, bindices.second}));
            }
        }
    }
   
    return wtasks;
}

auto
make_work_group(const std::vector<CBlockedGtoPairBlock>& gto_pair_blocks, const int ithreshold) -> std::vector<std::array<size_t, 8>>
{
    auto wtasks = std::vector<std::array<size_t, 8>>();
   
    if (const auto nblocks = gto_pair_blocks.size(); nblocks > 0)
    {
        const auto nthreads = omp::get_number_of_threads();
        
        for (size_t i = 0; i < nblocks; i++)
        {
            // apply threshold to blocked GTOs pair blocks on bra side

            for (int bra_idx = 0; bra_idx < 16; bra_idx++)
            {
                if (gto_pair_blocks[i].is_empty_gto_pair_block(bra_idx)) continue;

                const auto bra_block = gto_pair_blocks[i].gto_pair_block(bra_idx);
                
                const auto bra_angpair = bra_block.angular_momentums();

                const auto bra_size = bra_block.number_of_contracted_pairs();
                
                auto bra_bsize = omp::angular_momentum_scale(bra_angpair) * simd::width<double>();
                
                const auto bra_blocks = batch::number_of_batches(bra_size, bra_bsize);
                
                for (int ket_idx = bra_idx; ket_idx < 16; ket_idx++)
                {
                    if (gto_pair_blocks[i].is_empty_gto_pair_block(ket_idx)) continue;
                    
                    const auto ket_block = gto_pair_blocks[i].gto_pair_block(ket_idx);
                    
                    const auto ket_angpair = ket_block.angular_momentums();

                    const auto ket_size = ket_block.number_of_contracted_pairs();
                    
                    auto ket_bsize = omp::angular_momentum_scale(ket_angpair)  * simd::width<double>();
                    
                    const auto ket_blocks = batch::number_of_batches(ket_size, ket_bsize);
                    
                    // create task graph

                    if ((bra_idx + ket_idx) <= ithreshold)
                    {
                        for (size_t itask = 0; itask < bra_blocks; itask++)
                        {
                            const auto bindices = batch::batch_range(itask, bra_size, bra_bsize);
                            
                            const auto jstart = (bra_idx == ket_idx) ? itask : size_t{0};
                                              
                            for (size_t jtask = jstart; jtask < ket_blocks; jtask++)
                            {
                                const auto kindices = batch::batch_range(jtask, ket_size, ket_bsize);

                                wtasks.push_back({i, i,
                                                  static_cast<size_t>(bra_idx),
                                                  static_cast<size_t>(ket_idx),
                                                  bindices.first, bindices.second,
                                                  kindices.first, kindices.second});
                            }
                        }
                    }
                }
            }

            for (size_t j = i + 1; j < nblocks; j++)
            {
                // apply threshold to blocked GTOs pair blocks on bra and ket side

                for (int bra_idx = 0; bra_idx < 16; bra_idx++)
                {
                    if (gto_pair_blocks[i].is_empty_gto_pair_block(bra_idx)) continue;
                    
                    const auto bra_block = gto_pair_blocks[i].gto_pair_block(bra_idx);
                    
                    const auto bra_angpair = bra_block.angular_momentums();

                    const auto bra_size = bra_block.number_of_contracted_pairs();
                
                    auto bra_bsize = omp::angular_momentum_scale(bra_angpair) * simd::width<double>();
                    
                    const auto bra_blocks = batch::number_of_batches(bra_size, bra_bsize);

                    for (int ket_idx = 0; ket_idx < 16; ket_idx++)
                    {
                        if (gto_pair_blocks[j].is_empty_gto_pair_block(ket_idx)) continue;
                        
                        const auto ket_block = gto_pair_blocks[j].gto_pair_block(ket_idx);
                        
                        const auto ket_angpair = ket_block.angular_momentums();

                        const auto ket_size = ket_block.number_of_contracted_pairs();
                        
                        auto ket_bsize = omp::angular_momentum_scale(ket_angpair) * simd::width<double>();
                        
                        const auto ket_blocks = batch::number_of_batches(ket_size, ket_bsize);
                    
                        // create task graph

                        if ((bra_idx + ket_idx) <= ithreshold)
                        {
                            for (size_t itask = 0; itask < bra_blocks; itask++)
                            {
                                const auto bindices = batch::batch_range(itask, bra_size, bra_bsize);
                                
                                for (size_t jtask = 0; jtask < ket_blocks; jtask++)
                                {
                                    const auto kindices = batch::batch_range(jtask, ket_size, ket_bsize);

                                    wtasks.push_back({i, j,
                                                      static_cast<size_t>(bra_idx),
                                                      static_cast<size_t>(ket_idx),
                                                      bindices.first, bindices.second,
                                                      kindices.first, kindices.second});
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    return wtasks;
}

auto
angular_momentum_scale(const std::pair<int, int>& ang_pair) -> size_t
{
    const auto angmom = ang_pair.first + ang_pair.second;
    
    if (angmom > 8) return 16;
    
    if (angmom > 4) return 32;
    
    return 64;
}

}  // namespace omp
