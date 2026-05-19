//
//  Tabula — custom-recursion molecular-integral machinery.
//  Two-center overlap driver.
//

#include "TabulaOverlapDriver.hpp"

#include <algorithm>
#include <chrono>
#include <cstddef>
#include <vector>

#include "omp.h"

#include "GtoBlock.hpp"
#include "GtoFunc.hpp"
#include "Point.hpp"
#include "TabulaOverlapKernel.hpp"

namespace tabula {  // tabula namespace

namespace {  // unnamed namespace

/// @brief Wall seconds spanned by a steady-clock interval.
inline auto
seconds(const std::chrono::steady_clock::duration &interval) -> double
{
    return std::chrono::duration<double>(interval).count();
}

// per-thread load balance of the most recent compute's task loop
ThreadBalance g_balance;

/// @brief One basis-function block's primitive data, fetched once from a
/// `CGtoBlock` — the fused kernel reads it through an `OverlapBlockData` view.
struct BlockArrays
{
    /// @brief The primitive exponents, primitive-major.
    std::vector<double> exponents;
    /// @brief The primitive normalization factors.
    std::vector<double> norms;
    /// @brief The x / y / z coordinate of each contracted GTO.
    std::vector<double> x, y, z;
    /// @brief The AO indices — `[0]` the component stride, `[c+1]` contracted
    /// GTO `c`'s global AO offset.
    std::vector<std::size_t> orb_indices;
    /// @brief The number of contracted GTOs.
    int ncgtos{0};
    /// @brief The number of primitives per contracted GTO.
    int nprims{0};
    /// @brief The block angular momentum.
    int angular_momentum{0};

    /// @brief A pointer view of this block, as the kernel consumes it.
    auto view() const -> OverlapBlockData
    {
        return OverlapBlockData{exponents.data(), norms.data(), x.data(), y.data(), z.data(), ncgtos, nprims};
    }
};

/// @brief Fetches one basis-function block's primitive data.
auto
fetch_block(const CGtoBlock &block) -> BlockArrays
{
    BlockArrays arrays;
    arrays.exponents        = block.exponents();
    arrays.norms            = block.normalization_factors();
    arrays.orb_indices      = block.orbital_indices();
    arrays.ncgtos           = block.number_of_basis_functions();
    arrays.nprims           = block.number_of_primitives();
    arrays.angular_momentum = block.angular_momentum();

    const auto coords = block.coordinates();
    arrays.x.resize(static_cast<std::size_t>(arrays.ncgtos));
    arrays.y.resize(static_cast<std::size_t>(arrays.ncgtos));
    arrays.z.resize(static_cast<std::size_t>(arrays.ncgtos));
    for (int c = 0; c < arrays.ncgtos; c++)
    {
        const auto xyz = coords[static_cast<std::size_t>(c)].coordinates();
        arrays.x[static_cast<std::size_t>(c)] = xyz[0];
        arrays.y[static_cast<std::size_t>(c)] = xyz[1];
        arrays.z[static_cast<std::size_t>(c)] = xyz[2];
    }
    return arrays;
}

/// @brief Evaluates one task — one block pair restricted to the bra
/// contracted-GTO range `[bra_begin, bra_end)` — into the overlap matrix: the
/// fused overlap kernel followed by the scatter.
///
/// Each element is scattered once into the matrix's upper triangle; the lower
/// triangle is filled by a single `symmetrize` once every task is done. For a
/// `diagonal` task (a block paired with itself) a contracted pair and its
/// mirror both occur, so only the `r ≤ c` half is written. When `profile` is
/// non-null, the per-phase wall times are accumulated into it.
auto
evaluate_task(DenseMatrix       &matrix,
              const BlockArrays &bra,
              const BlockArrays &ket,
              const int          bra_begin,
              const int          bra_end,
              const bool         diagonal,
              OverlapProfile    *profile) -> void
{
    const int l_a = bra.angular_momentum;
    const int l_c = ket.angular_momentum;

    const std::size_t cdim =
        static_cast<std::size_t>(bra_end - bra_begin) * static_cast<std::size_t>(ket.ncgtos);
    if (cdim == 0) return;

    const auto t_start = std::chrono::steady_clock::now();

    // the fused overlap kernel — the primitive-pair weight, the seed ladder,
    // the primitive contraction, the single-centre MD recursion, and the
    // Cartesian-to-spherical assembly, in one tiled pass
    const auto stride         = ((cdim + 7) / 8) * 8;
    const auto bra_components = 2 * l_a + 1;
    const auto ket_components = 2 * l_c + 1;

    std::vector<double> spherical(static_cast<std::size_t>(bra_components * ket_components) * stride, 0.0);
    overlap_kernel(l_a, l_c, bra.view(), bra_begin, bra_end, ket.view(), spherical.data());

    const auto t_kernel = std::chrono::steady_clock::now();

    // scatter the spherical block into the matrix; the orbital indices carry
    // the AO component stride ([0]) and the per-GTO offset ([·+1]). The
    // contracted pair (i, j) is column (i-bra_begin)·kcgtos + j
    const int kcgtos = ket.ncgtos;
    for (int ca = 0; ca < bra_components; ca++)
    {
        for (int cc = 0; cc < ket_components; cc++)
        {
            const auto *row = spherical.data() + static_cast<std::size_t>(ca * ket_components + cc) * stride;

            for (int i = bra_begin; i < bra_end; i++)
            {
                const auto r = static_cast<std::size_t>(ca) * bra.orb_indices[0] +
                               bra.orb_indices[static_cast<std::size_t>(i) + 1];
                const std::size_t cp_base = static_cast<std::size_t>(i - bra_begin) * static_cast<std::size_t>(kcgtos);

                for (int j = 0; j < kcgtos; j++)
                {
                    const auto c = static_cast<std::size_t>(cc) * ket.orb_indices[0] +
                                   ket.orb_indices[static_cast<std::size_t>(j) + 1];

                    // a diagonal task also produces every mirror (j,i), so the
                    // r > c half is skipped — its mirror writes that slot
                    if (diagonal && r > c) continue;

                    matrix(std::min(r, c), std::max(r, c)) = row[cp_base + static_cast<std::size_t>(j)];
                }
            }
        }
    }

    if (profile != nullptr)
    {
        const auto t_scatter = std::chrono::steady_clock::now();
        profile->kernel  += seconds(t_kernel - t_start);
        profile->scatter += seconds(t_scatter - t_kernel);
    }
}

}  // namespace

auto
overlap_thread_balance() -> ThreadBalance
{
    return g_balance;
}

auto
OverlapDriver::compute(const CMolecule       &molecule,
                       const CMolecularBasis &basis,
                       const double           threshold,
                       OverlapProfile        *profile) const -> DenseMatrix
{
    // screening is not yet ported to the fused-kernel path — every pair is kept
    static_cast<void>(threshold);

    const auto t_blocks = std::chrono::steady_clock::now();

    const auto gto_blocks = gtofunc::make_gto_blocks(basis, molecule);

    if (profile != nullptr)
    {
        profile->make_blocks += seconds(std::chrono::steady_clock::now() - t_blocks);
    }

    const auto dimension = static_cast<std::size_t>(gtofunc::getNumberOfAtomicOrbitals(gto_blocks));

    DenseMatrix matrix(dimension, dimension, Symmetry::symmetric);

    // fetch each basis-function block's primitive data once — the kernel reads
    // it directly, so there is no per-pair block to construct
    const auto t_setup = std::chrono::steady_clock::now();

    std::vector<BlockArrays> blocks;
    blocks.reserve(gto_blocks.size());
    for (const auto &block : gto_blocks)
    {
        blocks.push_back(fetch_block(block));
    }

    if (profile != nullptr)
    {
        profile->pair_setup += seconds(std::chrono::steady_clock::now() - t_setup);
    }

    // the unit of parallel work is a task — one block pair restricted to a
    // sub-range of the bra block's contracted GTOs. The triangular block pairs
    // span orders of magnitude in cost; splitting the bra dimension so each
    // task holds roughly `target_pairs` contracted pairs turns dozens of
    // coarse units into thousands of even ones. Each task writes a disjoint AO
    // region, so the concurrent writes do not race.
    constexpr int target_pairs = 1 << 16;

    struct Task
    {
        std::size_t i;
        std::size_t j;
        int         bra_begin;
        int         bra_end;
        bool        diagonal;
        double      cost;
    };

    std::vector<Task> tasks;
    for (std::size_t i = 0; i < blocks.size(); i++)
    {
        for (std::size_t j = 0; j <= i; j++)
        {
            const int bcgtos = blocks[i].ncgtos;
            const int kcgtos = blocks[j].ncgtos;
            const int order  = blocks[i].angular_momentum + blocks[j].angular_momentum;

            const int chunk = std::max(1, target_pairs / std::max(1, kcgtos));

            for (int a = 0; a < bcgtos; a += chunk)
            {
                const int    b    = std::min(a + chunk, bcgtos);
                const double cost = static_cast<double>(b - a) * kcgtos * (order + 1) * (order + 2);
                tasks.push_back({i, j, a, b, i == j, cost});
            }
        }
    }

    // longest task first — dynamic scheduling then packs the small tasks into
    // the tail rather than stalling on a late big one
    std::sort(tasks.begin(), tasks.end(), [](const Task &lhs, const Task &rhs) { return lhs.cost > rhs.cost; });

    const auto nthreads = static_cast<std::size_t>(omp_get_max_threads());

    // per-thread profile accumulators, summed after the parallel region
    std::vector<OverlapProfile> thread_profiles;
    if (profile != nullptr)
    {
        thread_profiles.resize(nthreads);
    }

    // per-thread load-balance capture of the task loop
    std::vector<double> thread_busy(nthreads, 0.0);
    std::vector<long>   thread_pairs(nthreads, 0);

    const auto t_region = std::chrono::steady_clock::now();

#pragma omp parallel for schedule(dynamic)
    for (int p = 0; p < static_cast<int>(tasks.size()); p++)
    {
        const auto tid    = static_cast<std::size_t>(omp_get_thread_num());
        const auto t_body = std::chrono::steady_clock::now();

        OverlapProfile *slot = (profile != nullptr) ? &thread_profiles[tid] : nullptr;

        const auto &task = tasks[static_cast<std::size_t>(p)];

        evaluate_task(matrix, blocks[task.i], blocks[task.j], task.bra_begin, task.bra_end, task.diagonal, slot);

        thread_busy[tid] += seconds(std::chrono::steady_clock::now() - t_body);
        thread_pairs[tid] += 1;
    }

    g_balance.wall  = seconds(std::chrono::steady_clock::now() - t_region);
    g_balance.busy  = thread_busy;
    g_balance.pairs = thread_pairs;

    // the scatter filled the upper triangle only — mirror it into the lower
    const auto t_symmetrize = std::chrono::steady_clock::now();
    matrix.symmetrize();

    if (profile != nullptr)
    {
        profile->symmetrize += seconds(std::chrono::steady_clock::now() - t_symmetrize);

        for (const auto &thread_profile : thread_profiles)
        {
            profile->kernel  += thread_profile.kernel;
            profile->scatter += thread_profile.scatter;
        }
    }

    return matrix;
}

}  // namespace tabula
