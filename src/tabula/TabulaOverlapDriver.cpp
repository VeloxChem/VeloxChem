//
//  Tabula — custom-recursion molecular-integral machinery.
//  Two-center overlap driver.
//

#include "TabulaOverlapDriver.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

#include "omp.h"

#include "GtoBlock.hpp"
#include "GtoFunc.hpp"
#include "MathConst.hpp"
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

/// @brief A maximal contiguous run of one atom's contracted GTOs within a
/// block — the unit a task's bra side covers, and the unit ket screening
/// keeps or drops. Carries the atom centre and a conservative screening
/// datum folded over the run.
struct AtomSpan
{
    /// @brief The molecular atom index.
    int atom{0};
    /// @brief The first / one-past-last contracted-GTO index of the run.
    int cgto_begin{0}, cgto_end{0};
    /// @brief The atom centre.
    double x{0.0}, y{0.0}, z{0.0};
    /// @brief The largest contraction-coefficient magnitude over the run.
    double cmax{0.0};
    /// @brief The smallest / largest primitive exponent over the run.
    double exp_min{0.0}, exp_max{0.0};
};

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
    /// @brief The maximal contiguous same-atom runs of contracted GTOs.
    std::vector<AtomSpan> spans;
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

/// @brief Fetches one basis-function block's primitive data and atom spans.
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

    // group the contracted GTOs into maximal contiguous same-atom runs and
    // fold a conservative screening datum over each run — the largest
    // coefficient and the widest exponent span
    const auto atoms = block.atomic_indices();
    for (int c = 0; c < arrays.ncgtos; c++)
    {
        if (c == 0 || atoms[static_cast<std::size_t>(c)] != atoms[static_cast<std::size_t>(c - 1)])
        {
            AtomSpan span;
            span.atom       = atoms[static_cast<std::size_t>(c)];
            span.cgto_begin = c;
            span.cgto_end   = c;
            span.x          = arrays.x[static_cast<std::size_t>(c)];
            span.y          = arrays.y[static_cast<std::size_t>(c)];
            span.z          = arrays.z[static_cast<std::size_t>(c)];
            span.cmax       = 0.0;
            span.exp_min    = std::numeric_limits<double>::max();
            span.exp_max    = 0.0;
            arrays.spans.push_back(span);
        }

        auto      &span = arrays.spans.back();
        const auto data = block.screening_data(static_cast<std::size_t>(c));
        span.cgto_end   = c + 1;
        span.cmax       = std::max(span.cmax, data.max_coefficient);
        span.exp_min    = std::min(span.exp_min, data.min_exponent);
        span.exp_max    = std::max(span.exp_max, data.max_exponent);
    }

    return arrays;
}

/// @brief A conservative upper bound on the overlap magnitude of the
/// shell-pair sub-block formed by bra span `A` and ket span `B`, separated by
/// `r² = r2`. See `docs/screening-predicates.md`; using the slowest-decaying
/// (smallest) exponents keeps it an upper bound, so a sub-block whose estimate
/// is below the threshold is provably negligible.
auto
overlap_estimate(const int l_a, const int l_c, const AtomSpan &A, const AtomSpan &B, const double r2) -> double
{
    const int    m       = l_a + l_c;
    const double r       = std::sqrt(r2);
    const double rho_min = A.exp_min * B.exp_min / (A.exp_min + B.exp_min);
    const double rho_max = A.exp_max * B.exp_max / (A.exp_max + B.exp_max);

    double estimate = A.cmax * B.cmax;

    // r^m bounds the R-monomial growth; (2 rho_max)^m the seed ladder
    const double two_rho = 2.0 * rho_max;
    for (int t = 0; t < m; t++) estimate *= r * two_rho;

    // kernelBound — the Gaussian decay
    estimate *= std::exp(-rho_min * r2);

    // the Gaussian-product prefactor magnitude at the slowest-decaying exponents
    const double pe = mathconst::pi_value() / (A.exp_min + B.exp_min);
    estimate *= pe * std::sqrt(pe);

    const double ia = 0.5 / A.exp_min;
    for (int t = 0; t < l_a; t++) estimate *= ia;

    const double ib = 0.5 / B.exp_min;
    for (int t = 0; t < l_c; t++) estimate *= ib;

    return estimate;
}

/// @brief Evaluates one task — bra span `A` of `bra` against the ket block —
/// into the overlap matrix: the ket atom spans are screened against `A`, the
/// survivors are gathered into a contiguous block, the fused kernel runs, and
/// the spherical block is scattered.
///
/// Each element is scattered once into the matrix's upper triangle. For a
/// `diagonal` task (a block paired with itself) a contracted pair and its
/// mirror both occur, so only the `r ≤ c` half is written. When `profile` is
/// non-null, the per-phase wall times are accumulated into it.
auto
evaluate_task(DenseMatrix       &matrix,
              const BlockArrays &bra,
              const AtomSpan    &A,
              const BlockArrays &ket,
              const bool         diagonal,
              const double       threshold,
              OverlapProfile    *profile) -> void
{
    const int l_a   = bra.angular_momentum;
    const int l_c   = ket.angular_momentum;
    const int k_bra = A.cgto_end - A.cgto_begin;

    const auto t_start = std::chrono::steady_clock::now();

    // screen the ket atom spans against the bra span — a same-atom span and
    // any span whose conservative estimate clears the threshold are kept; a
    // threshold at or below 0 keeps every span and skips the estimate
    thread_local std::vector<int> orig;
    orig.clear();
    const bool screen = threshold > 0.0;
    for (const auto &B : ket.spans)
    {
        bool keep = !screen || (A.atom == B.atom);
        if (!keep)
        {
            const double dx = A.x - B.x;
            const double dy = A.y - B.y;
            const double dz = A.z - B.z;
            keep = overlap_estimate(l_a, l_c, A, B, dx * dx + dy * dy + dz * dz) >= threshold;
        }
        if (keep)
        {
            for (int c = B.cgto_begin; c < B.cgto_end; c++) orig.push_back(c);
        }
    }

    const int ket_kept = static_cast<int>(orig.size());
    if (ket_kept == 0) return;

    // gather the surviving ket contracted GTOs into a contiguous block; when
    // every span survives the kernel reads the block directly
    OverlapBlockData ket_view = ket.view();
    if (ket_kept != ket.ncgtos)
    {
        thread_local std::vector<double> cmp_exp, cmp_norm, cmp_x, cmp_y, cmp_z;
        cmp_exp.resize(static_cast<std::size_t>(ket.nprims) * ket_kept);
        cmp_norm.resize(static_cast<std::size_t>(ket.nprims) * ket_kept);
        cmp_x.resize(static_cast<std::size_t>(ket_kept));
        cmp_y.resize(static_cast<std::size_t>(ket_kept));
        cmp_z.resize(static_cast<std::size_t>(ket_kept));

        for (int c = 0; c < ket_kept; c++)
        {
            const int kc = orig[static_cast<std::size_t>(c)];
            cmp_x[static_cast<std::size_t>(c)] = ket.x[static_cast<std::size_t>(kc)];
            cmp_y[static_cast<std::size_t>(c)] = ket.y[static_cast<std::size_t>(kc)];
            cmp_z[static_cast<std::size_t>(c)] = ket.z[static_cast<std::size_t>(kc)];
        }
        for (int lk = 0; lk < ket.nprims; lk++)
        {
            for (int c = 0; c < ket_kept; c++)
            {
                const int         kc  = orig[static_cast<std::size_t>(c)];
                const std::size_t dst = static_cast<std::size_t>(lk) * ket_kept + c;
                cmp_exp[dst]  = ket.exponents[static_cast<std::size_t>(lk) * ket.ncgtos + kc];
                cmp_norm[dst] = ket.norms[static_cast<std::size_t>(lk) * ket.ncgtos + kc];
            }
        }
        ket_view = OverlapBlockData{cmp_exp.data(), cmp_norm.data(), cmp_x.data(), cmp_y.data(), cmp_z.data(),
                                    ket_kept, ket.nprims};
    }

    const auto t_screen = std::chrono::steady_clock::now();

    // the fused overlap kernel — bra span A against the kept ket block. The
    // kernel writes every component's [0, cdim) so the (grow-only,
    // thread-local) spherical buffer needs no zero-fill
    const std::size_t cdim   = static_cast<std::size_t>(k_bra) * static_cast<std::size_t>(ket_kept);
    const auto        stride = ((cdim + 7) / 8) * 8;
    const auto        bra_components = 2 * l_a + 1;
    const auto        ket_components = 2 * l_c + 1;

    thread_local std::vector<double> spherical;
    const std::size_t sph_size = static_cast<std::size_t>(bra_components * ket_components) * stride;
    if (spherical.size() < sph_size) spherical.resize(sph_size);

    overlap_kernel(l_a, l_c, bra.view(), A.cgto_begin, A.cgto_end, ket_view, spherical.data());

    const auto t_kernel = std::chrono::steady_clock::now();

    // scatter — the contracted pair (bra-local, c) is column bra-local·ket_kept + c;
    // its ket contracted GTO is orig[c]
    for (int ca = 0; ca < bra_components; ca++)
    {
        for (int cc = 0; cc < ket_components; cc++)
        {
            const auto *row = spherical.data() + static_cast<std::size_t>(ca * ket_components + cc) * stride;

            for (int bl = 0; bl < k_bra; bl++)
            {
                const auto r = static_cast<std::size_t>(ca) * bra.orb_indices[0] +
                               bra.orb_indices[static_cast<std::size_t>(A.cgto_begin + bl) + 1];
                const std::size_t cp_base = static_cast<std::size_t>(bl) * static_cast<std::size_t>(ket_kept);

                for (int c = 0; c < ket_kept; c++)
                {
                    const int  kc = orig[static_cast<std::size_t>(c)];
                    const auto cao = static_cast<std::size_t>(cc) * ket.orb_indices[0] +
                                     ket.orb_indices[static_cast<std::size_t>(kc) + 1];

                    if (diagonal && r > cao) continue;

                    matrix(std::min(r, cao), std::max(r, cao)) = row[cp_base + static_cast<std::size_t>(c)];
                }
            }
        }
    }

    if (profile != nullptr)
    {
        const auto t_scatter = std::chrono::steady_clock::now();
        profile->screen  += seconds(t_screen - t_start);
        profile->kernel  += seconds(t_kernel - t_screen);
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
    const auto t_blocks = std::chrono::steady_clock::now();

    const auto gto_blocks = gtofunc::make_gto_blocks(basis, molecule);

    if (profile != nullptr)
    {
        profile->make_blocks += seconds(std::chrono::steady_clock::now() - t_blocks);
    }

    const auto dimension = static_cast<std::size_t>(gtofunc::getNumberOfAtomicOrbitals(gto_blocks));

    DenseMatrix matrix(dimension, dimension, Symmetry::symmetric);

    // fetch each basis-function block's primitive data and atom spans once
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

    // the unit of parallel work is a task — one block pair restricted to one
    // bra atom span. Each task screens the ket atom spans against its bra
    // span, so far atoms never reach the kernel. Each task writes a disjoint
    // AO region, so the concurrent writes do not race.
    struct Task
    {
        std::size_t i;
        std::size_t j;
        std::size_t span;
        bool        diagonal;
        double      cost;
    };

    std::vector<Task> tasks;
    for (std::size_t i = 0; i < blocks.size(); i++)
    {
        for (std::size_t j = 0; j <= i; j++)
        {
            const int kcgtos = blocks[j].ncgtos;
            const int order  = blocks[i].angular_momentum + blocks[j].angular_momentum;

            for (std::size_t span = 0; span < blocks[i].spans.size(); span++)
            {
                const auto  &A    = blocks[i].spans[span];
                const double cost = static_cast<double>(A.cgto_end - A.cgto_begin) * kcgtos * (order + 1) * (order + 2);
                tasks.push_back({i, j, span, i == j, cost});
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

        evaluate_task(matrix,
                      blocks[task.i],
                      blocks[task.i].spans[task.span],
                      blocks[task.j],
                      task.diagonal,
                      threshold,
                      slot);

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
            profile->screen  += thread_profile.screen;
            profile->kernel  += thread_profile.kernel;
            profile->scatter += thread_profile.scatter;
        }
    }

    return matrix;
}

}  // namespace tabula
