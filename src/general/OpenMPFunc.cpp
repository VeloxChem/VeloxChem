#include "OpenMPFunc.hpp"

#include <algorithm>
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

}  // namespace omp
