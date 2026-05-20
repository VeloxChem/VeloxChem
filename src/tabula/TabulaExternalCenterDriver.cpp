//
//  Tabula — custom-recursion molecular-integral machinery.
//  Shared external-center driver core — block fetch, atom-span screening, the
//  task list, the parallel loop, and the dense / block-sparse scatter. The
//  operator supplies the fused kernel (which sums over the point charges) and
//  the screening estimate (`ExternalCenterOperator`). Sibling of the
//  two-center core.
//

#include "TabulaExternalCenterDriver.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <limits>
#include <utility>
#include <vector>

#include "omp.h"

#include "GtoBlock.hpp"
#include "GtoFunc.hpp"
#include "Point.hpp"

namespace tabula {  // tabula namespace

namespace {  // unnamed namespace

/// @brief Wall seconds spanned by a steady-clock interval.
inline auto
seconds(const std::chrono::steady_clock::duration &interval) -> double
{
    return std::chrono::duration<double>(interval).count();
}

/// @brief One basis-function block's primitive data, fetched once from a
/// `CGtoBlock` — the fused kernel reads it through a `KernelBlockData` view.
struct BlockArrays
{
    std::vector<double>      exponents;
    std::vector<double>      norms;
    std::vector<double>      x, y, z;
    std::vector<std::size_t> orb_indices;
    std::vector<int>         atom_of_cgto;
    std::vector<AtomSpan>    spans;
    int                      ncgtos{0};
    int                      nprims{0};
    int                      angular_momentum{0};

    auto view() const -> KernelBlockData
    {
        return KernelBlockData{exponents.data(), norms.data(), x.data(), y.data(), z.data(), ncgtos, nprims};
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
    arrays.atom_of_cgto     = block.atomic_indices();
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

    for (int c = 0; c < arrays.ncgtos; c++)
    {
        if (c == 0 || arrays.atom_of_cgto[static_cast<std::size_t>(c)] != arrays.atom_of_cgto[static_cast<std::size_t>(c - 1)])
        {
            AtomSpan span;
            span.atom       = arrays.atom_of_cgto[static_cast<std::size_t>(c)];
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

/// @brief Whether ket span `B` is kept against bra span `A`.
inline auto
ket_kept(const int                     l_a,
         const int                     l_c,
         const AtomSpan               &A,
         const AtomSpan               &B,
         const double                  threshold,
         const bool                    restrict_ket_le_bra,
         const ExternalCenterOperator &op) -> bool
{
    if (restrict_ket_le_bra && B.atom > A.atom) return false;
    if (threshold <= 0.0 || A.atom == B.atom) return true;

    const double dx = A.x - B.x;
    const double dy = A.y - B.y;
    const double dz = A.z - B.z;
    return op.estimate(l_a, l_c, A, B, dx * dx + dy * dy + dz * dz, op.charge_factor) >= threshold;
}

/// @brief Screens the ket atom spans against bra span `A`, gathers the
/// surviving ket contracted GTOs into a contiguous block, and runs the fused
/// kernel (which sums over its bound point sources). Returns the kept ket
/// contracted-GTO count.
auto
screen_and_kernel(const BlockArrays            &bra,
                  const AtomSpan               &A,
                  const BlockArrays            &ket,
                  const double                  threshold,
                  const bool                    restrict_ket_le_bra,
                  const ExternalCenterOperator &op,
                  std::vector<int>             &orig,
                  std::vector<double>          &spherical,
                  KernelProfile                *profile) -> int
{
    const int l_a   = bra.angular_momentum;
    const int l_c   = ket.angular_momentum;
    const int k_bra = A.cgto_end - A.cgto_begin;

    const auto t_start = std::chrono::steady_clock::now();

    orig.clear();
    for (const auto &B : ket.spans)
    {
        if (ket_kept(l_a, l_c, A, B, threshold, restrict_ket_le_bra, op))
        {
            for (int c = B.cgto_begin; c < B.cgto_end; c++) orig.push_back(c);
        }
    }

    const int kept = static_cast<int>(orig.size());
    if (kept == 0) return 0;

    KernelBlockData ket_view = ket.view();
    if (kept != ket.ncgtos)
    {
        thread_local std::vector<double> cmp_exp, cmp_norm, cmp_x, cmp_y, cmp_z;
        cmp_exp.resize(static_cast<std::size_t>(ket.nprims) * kept);
        cmp_norm.resize(static_cast<std::size_t>(ket.nprims) * kept);
        cmp_x.resize(static_cast<std::size_t>(kept));
        cmp_y.resize(static_cast<std::size_t>(kept));
        cmp_z.resize(static_cast<std::size_t>(kept));

        for (int c = 0; c < kept; c++)
        {
            const int kc = orig[static_cast<std::size_t>(c)];
            cmp_x[static_cast<std::size_t>(c)] = ket.x[static_cast<std::size_t>(kc)];
            cmp_y[static_cast<std::size_t>(c)] = ket.y[static_cast<std::size_t>(kc)];
            cmp_z[static_cast<std::size_t>(c)] = ket.z[static_cast<std::size_t>(kc)];
        }
        for (int lk = 0; lk < ket.nprims; lk++)
        {
            for (int c = 0; c < kept; c++)
            {
                const int         kc  = orig[static_cast<std::size_t>(c)];
                const std::size_t dst = static_cast<std::size_t>(lk) * kept + c;
                cmp_exp[dst]  = ket.exponents[static_cast<std::size_t>(lk) * ket.ncgtos + kc];
                cmp_norm[dst] = ket.norms[static_cast<std::size_t>(lk) * ket.ncgtos + kc];
            }
        }
        ket_view = KernelBlockData{cmp_exp.data(), cmp_norm.data(), cmp_x.data(), cmp_y.data(), cmp_z.data(),
                                   kept, ket.nprims};
    }

    const auto t_screen = std::chrono::steady_clock::now();

    const std::size_t cdim   = static_cast<std::size_t>(k_bra) * static_cast<std::size_t>(kept);
    const auto        stride = ((cdim + 7) / 8) * 8;
    const auto        components = static_cast<std::size_t>((2 * l_a + 1) * (2 * l_c + 1));

    if (spherical.size() < components * stride) spherical.resize(components * stride);

    op.kernel(l_a, l_c, bra.view(), A.cgto_begin, A.cgto_end, ket_view, spherical.data());

    if (profile != nullptr)
    {
        const auto t_kernel = std::chrono::steady_clock::now();
        profile->screen += seconds(t_screen - t_start);
        profile->kernel += seconds(t_kernel - t_screen);
    }

    return kept;
}

/// @brief Evaluates one task into the dense matrix.
auto
evaluate_task(DenseMatrix                  &matrix,
              const BlockArrays            &bra,
              const AtomSpan               &A,
              const BlockArrays            &ket,
              const bool                    diagonal,
              const double                  threshold,
              const ExternalCenterOperator &op,
              KernelProfile                *profile) -> void
{
    thread_local std::vector<int>    orig;
    thread_local std::vector<double> spherical;

    const int kept = screen_and_kernel(bra, A, ket, threshold, false, op, orig, spherical, profile);
    if (kept == 0) return;

    const auto t_kernel = std::chrono::steady_clock::now();

    const int         l_a    = bra.angular_momentum;
    const int         l_c    = ket.angular_momentum;
    const int         k_bra  = A.cgto_end - A.cgto_begin;
    const std::size_t cdim   = static_cast<std::size_t>(k_bra) * static_cast<std::size_t>(kept);
    const auto        stride = ((cdim + 7) / 8) * 8;
    const auto        bra_components = 2 * l_a + 1;
    const auto        ket_components = 2 * l_c + 1;

    double *const     m   = matrix.values();
    const std::size_t dim = matrix.columns();
    for (int ca = 0; ca < bra_components; ca++)
    {
        for (int cc = 0; cc < ket_components; cc++)
        {
            const auto *row = spherical.data() + static_cast<std::size_t>(ca * ket_components + cc) * stride;

            for (int bl = 0; bl < k_bra; bl++)
            {
                const auto r = static_cast<std::size_t>(ca) * bra.orb_indices[0] +
                               bra.orb_indices[static_cast<std::size_t>(A.cgto_begin + bl) + 1];
                const std::size_t cp_base = static_cast<std::size_t>(bl) * static_cast<std::size_t>(kept);

                for (int c = 0; c < kept; c++)
                {
                    const int  kc  = orig[static_cast<std::size_t>(c)];
                    const auto cao = static_cast<std::size_t>(cc) * ket.orb_indices[0] +
                                     ket.orb_indices[static_cast<std::size_t>(kc) + 1];

                    if (diagonal && r > cao) continue;

                    const double v = row[cp_base + static_cast<std::size_t>(c)];
                    m[r * dim + cao] = v;
                    if (r != cao) m[cao * dim + r] = v;
                }
            }
        }
    }

    if (profile != nullptr)
    {
        profile->scatter += seconds(std::chrono::steady_clock::now() - t_kernel);
    }
}

/// @brief Evaluates one task into the block-sparse matrix.
auto
evaluate_task_sparse(BlockSparseMatrix              &matrix,
                     const BlockArrays              &bra,
                     const AtomSpan                 &A,
                     const BlockArrays              &ket,
                     const bool                      diagonal_block_pair,
                     const double                    threshold,
                     const ExternalCenterOperator   &op,
                     const std::vector<std::size_t> &ao_local,
                     const std::vector<long>        &block_of,
                     const std::size_t               n_atoms,
                     KernelProfile                  *profile) -> void
{
    thread_local std::vector<int>    orig;
    thread_local std::vector<double> spherical;

    const int kept = screen_and_kernel(bra, A, ket, threshold, diagonal_block_pair, op, orig, spherical, profile);
    if (kept == 0) return;

    const auto t_kernel = std::chrono::steady_clock::now();

    const int         l_a    = bra.angular_momentum;
    const int         l_c    = ket.angular_momentum;
    const int         k_bra  = A.cgto_end - A.cgto_begin;
    const std::size_t cdim   = static_cast<std::size_t>(k_bra) * static_cast<std::size_t>(kept);
    const auto        stride = ((cdim + 7) / 8) * 8;
    const auto        bra_components = 2 * l_a + 1;
    const auto        ket_components = 2 * l_c + 1;

    double *const values = matrix.values();
    for (int c = 0; c < kept; c++)
    {
        const int  kc       = orig[static_cast<std::size_t>(c)];
        const int  ket_atom = ket.atom_of_cgto[static_cast<std::size_t>(kc)];
        const auto group_a  = static_cast<std::size_t>(std::max(A.atom, ket_atom));
        const auto group_b  = static_cast<std::size_t>(std::min(A.atom, ket_atom));
        const long block    = block_of[group_a * n_atoms + group_b];
        if (block < 0) continue;

        const bool transpose = A.atom < ket_atom;
        const bool diagonal  = A.atom == ket_atom;

        const auto        blk  = matrix.block(static_cast<std::size_t>(block));
        double *const     base = values + blk.offset;
        const std::size_t ncol = blk.columnCount;

        for (int ca = 0; ca < bra_components; ca++)
        {
            for (int bl = 0; bl < k_bra; bl++)
            {
                const auto r = static_cast<std::size_t>(ca) * bra.orb_indices[0] +
                               bra.orb_indices[static_cast<std::size_t>(A.cgto_begin + bl) + 1];
                const auto lr = ao_local[r];

                for (int cc = 0; cc < ket_components; cc++)
                {
                    const auto cao = static_cast<std::size_t>(cc) * ket.orb_indices[0] +
                                     ket.orb_indices[static_cast<std::size_t>(kc) + 1];
                    const auto lc = ao_local[cao];

                    const auto v = spherical[static_cast<std::size_t>(ca * ket_components + cc) * stride +
                                             static_cast<std::size_t>(bl) * static_cast<std::size_t>(kept) +
                                             static_cast<std::size_t>(c)];

                    if (transpose)
                    {
                        base[lc * ncol + lr] = v;
                    }
                    else
                    {
                        base[lr * ncol + lc] = v;
                        if (diagonal) base[lc * ncol + lr] = v;
                    }
                }
            }
        }
    }

    if (profile != nullptr)
    {
        profile->scatter += seconds(std::chrono::steady_clock::now() - t_kernel);
    }
}

/// @brief The triangular block-pair / bra-span task list, longest first.
struct Task
{
    std::size_t i;
    std::size_t j;
    std::size_t span;
    bool        diagonal;
    double      cost;
};

auto
build_tasks(const std::vector<BlockArrays> &blocks) -> std::vector<Task>
{
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
    std::sort(tasks.begin(), tasks.end(), [](const Task &lhs, const Task &rhs) { return lhs.cost > rhs.cost; });
    return tasks;
}

}  // namespace

auto
external_center_compute(const CMolecule&              molecule,
                        const CMolecularBasis&        basis,
                        const double                  threshold,
                        const ExternalCenterOperator& op,
                        KernelProfile*                profile,
                        ThreadBalance&                balance) -> DenseMatrix
{
    const auto t_blocks = std::chrono::steady_clock::now();

    const auto gto_blocks = gtofunc::make_gto_blocks(basis, molecule);

    if (profile != nullptr)
    {
        profile->make_blocks += seconds(std::chrono::steady_clock::now() - t_blocks);
    }

    const auto dimension = static_cast<std::size_t>(gtofunc::getNumberOfAtomicOrbitals(gto_blocks));

    const bool initialize = threshold > 0.0;

    DenseMatrix matrix(dimension, dimension, Symmetry::symmetric, initialize);

    const auto t_setup = std::chrono::steady_clock::now();

    std::vector<BlockArrays> blocks;
    blocks.reserve(gto_blocks.size());
    for (const auto &block : gto_blocks)
    {
        blocks.push_back(fetch_block(block));
    }

    const auto tasks = build_tasks(blocks);

    if (profile != nullptr)
    {
        profile->pair_setup += seconds(std::chrono::steady_clock::now() - t_setup);
    }

    const auto nthreads = static_cast<std::size_t>(omp_get_max_threads());

    std::vector<KernelProfile> thread_profiles;
    if (profile != nullptr)
    {
        thread_profiles.resize(nthreads);
    }

    std::vector<double> thread_busy(nthreads, 0.0);
    std::vector<long>   thread_pairs(nthreads, 0);

    const auto t_region = std::chrono::steady_clock::now();

#pragma omp parallel for schedule(dynamic)
    for (int p = 0; p < static_cast<int>(tasks.size()); p++)
    {
        const auto tid    = static_cast<std::size_t>(omp_get_thread_num());
        const auto t_body = std::chrono::steady_clock::now();

        KernelProfile *slot = (profile != nullptr) ? &thread_profiles[tid] : nullptr;

        const auto &task = tasks[static_cast<std::size_t>(p)];

        evaluate_task(matrix, blocks[task.i], blocks[task.i].spans[task.span], blocks[task.j], task.diagonal, threshold, op, slot);

        thread_busy[tid] += seconds(std::chrono::steady_clock::now() - t_body);
        thread_pairs[tid] += 1;
    }

    balance.wall  = seconds(std::chrono::steady_clock::now() - t_region);
    balance.busy  = thread_busy;
    balance.pairs = thread_pairs;

    if (profile != nullptr)
    {
        for (const auto &thread_profile : thread_profiles)
        {
            profile->screen  += thread_profile.screen;
            profile->kernel  += thread_profile.kernel;
            profile->scatter += thread_profile.scatter;
        }
    }

    return matrix;
}

auto
external_center_compute_sparse(const CMolecule&              molecule,
                               const CMolecularBasis&        basis,
                               const double                  threshold,
                               const ExternalCenterOperator& op,
                               KernelProfile*                profile,
                               ThreadBalance&                balance) -> BlockSparseMatrix
{
    const auto t_blocks = std::chrono::steady_clock::now();

    const auto gto_blocks = gtofunc::make_gto_blocks(basis, molecule);

    if (profile != nullptr)
    {
        profile->make_blocks += seconds(std::chrono::steady_clock::now() - t_blocks);
    }

    const auto dimension = static_cast<std::size_t>(gtofunc::getNumberOfAtomicOrbitals(gto_blocks));

    const auto t_setup = std::chrono::steady_clock::now();

    std::vector<BlockArrays> blocks;
    blocks.reserve(gto_blocks.size());
    for (const auto &block : gto_blocks)
    {
        blocks.push_back(fetch_block(block));
    }

    std::size_t n_atoms = 0;
    for (const auto &block : blocks)
    {
        for (const auto &span : block.spans)
        {
            n_atoms = std::max(n_atoms, static_cast<std::size_t>(span.atom) + 1);
        }
    }

    std::vector<std::vector<std::size_t>> group_global_ao(n_atoms);
    std::vector<std::size_t>              ao_local(dimension, 0);
    for (const auto &block : blocks)
    {
        const int components = 2 * block.angular_momentum + 1;
        for (const auto &span : block.spans)
        {
            auto &group = group_global_ao[static_cast<std::size_t>(span.atom)];
            for (int cgto = span.cgto_begin; cgto < span.cgto_end; cgto++)
            {
                for (int ca = 0; ca < components; ca++)
                {
                    const auto ao = static_cast<std::size_t>(ca) * block.orb_indices[0] +
                                    block.orb_indices[static_cast<std::size_t>(cgto) + 1];
                    ao_local[ao] = group.size();
                    group.push_back(ao);
                }
            }
        }
    }

    const auto tasks = build_tasks(blocks);

    std::vector<char> significant(n_atoms * n_atoms, 0);
#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < static_cast<int>(blocks.size()); i++)
    {
        for (std::size_t j = 0; j <= static_cast<std::size_t>(i); j++)
        {
            const auto &bra = blocks[static_cast<std::size_t>(i)];
            const auto &ket = blocks[j];
            const bool  diagonal_block_pair = (static_cast<std::size_t>(i) == j);

            for (const auto &A : bra.spans)
            {
                for (const auto &B : ket.spans)
                {
                    if (ket_kept(bra.angular_momentum, ket.angular_momentum, A, B, threshold, diagonal_block_pair, op))
                    {
                        const auto ga = static_cast<std::size_t>(std::max(A.atom, B.atom));
                        const auto gb = static_cast<std::size_t>(std::min(A.atom, B.atom));
                        significant[ga * n_atoms + gb] = 1;
                    }
                }
            }
        }
    }

    std::vector<std::pair<std::size_t, std::size_t>> group_pairs;
    std::vector<long>                                block_of(n_atoms * n_atoms, -1);
    for (std::size_t ga = 0; ga < n_atoms; ga++)
    {
        for (std::size_t gb = 0; gb <= ga; gb++)
        {
            if (significant[ga * n_atoms + gb] != 0)
            {
                block_of[ga * n_atoms + gb] = static_cast<long>(group_pairs.size());
                group_pairs.push_back({ga, gb});
            }
        }
    }

    BlockSparseMatrix matrix(dimension, group_global_ao, group_pairs);

    if (profile != nullptr)
    {
        profile->pair_setup += seconds(std::chrono::steady_clock::now() - t_setup);
    }

    const auto nthreads = static_cast<std::size_t>(omp_get_max_threads());

    std::vector<KernelProfile> thread_profiles;
    if (profile != nullptr)
    {
        thread_profiles.resize(nthreads);
    }

    std::vector<double> thread_busy(nthreads, 0.0);
    std::vector<long>   thread_pairs(nthreads, 0);

    const auto t_region = std::chrono::steady_clock::now();

#pragma omp parallel for schedule(dynamic)
    for (int p = 0; p < static_cast<int>(tasks.size()); p++)
    {
        const auto tid    = static_cast<std::size_t>(omp_get_thread_num());
        const auto t_body = std::chrono::steady_clock::now();

        KernelProfile *slot = (profile != nullptr) ? &thread_profiles[tid] : nullptr;

        const auto &task = tasks[static_cast<std::size_t>(p)];

        evaluate_task_sparse(matrix,
                             blocks[task.i],
                             blocks[task.i].spans[task.span],
                             blocks[task.j],
                             task.diagonal,
                             threshold,
                             op,
                             ao_local,
                             block_of,
                             n_atoms,
                             slot);

        thread_busy[tid] += seconds(std::chrono::steady_clock::now() - t_body);
        thread_pairs[tid] += 1;
    }

    balance.wall  = seconds(std::chrono::steady_clock::now() - t_region);
    balance.busy  = thread_busy;
    balance.pairs = thread_pairs;

    if (profile != nullptr)
    {
        for (const auto &thread_profile : thread_profiles)
        {
            profile->screen  += thread_profile.screen;
            profile->kernel  += thread_profile.kernel;
            profile->scatter += thread_profile.scatter;
        }
    }

    return matrix;
}

auto
external_center_field_compute(const CMolecule&         molecule,
                              const CMolecularBasis&   basis,
                              const DenseMatrix&       density,
                              const FieldKernelFn      kernel,
                              const ExternalEstimateFn estimate,
                              const double*            point_x,
                              const double*            point_y,
                              const double*            point_z,
                              const int                n_points,
                              const double             threshold,
                              ThreadBalance&           balance) -> std::vector<std::array<double, 3>>
{
    const auto gto_blocks = gtofunc::make_gto_blocks(basis, molecule);

    std::vector<BlockArrays> blocks;
    blocks.reserve(gto_blocks.size());
    for (const auto& block : gto_blocks)
    {
        blocks.push_back(fetch_block(block));
    }

    const auto tasks = build_tasks(blocks);

    const auto        nthreads = static_cast<std::size_t>(omp_get_max_threads());
    const std::size_t np       = static_cast<std::size_t>(n_points);
    const bool        screen   = threshold > 0.0;

    // each thread accumulates the full per-point field (axis-major), reduced after
    std::vector<std::vector<double>> thread_field(nthreads, std::vector<double>(3 * np, 0.0));
    std::vector<double>              thread_busy(nthreads, 0.0);
    std::vector<long>                thread_pairs(nthreads, 0);

    const double* const d    = density.values();
    const std::size_t   n_ao = density.columns();

    const auto t_region = std::chrono::steady_clock::now();

#pragma omp parallel
    {
        std::vector<double> dblk;

#pragma omp for schedule(dynamic)
        for (int p = 0; p < static_cast<int>(tasks.size()); p++)
        {
            const auto tid    = static_cast<std::size_t>(omp_get_thread_num());
            const auto t_body = std::chrono::steady_clock::now();

            const auto &task = tasks[static_cast<std::size_t>(p)];
            const auto &bra  = blocks[task.i];
            const auto &A    = bra.spans[task.span];
            const auto &ket  = blocks[task.j];

            const int la        = bra.angular_momentum;
            const int lc        = ket.angular_momentum;
            const int k_bra     = A.cgto_end - A.cgto_begin;
            const int bra_comps = 2 * la + 1;
            const int ket_comps = 2 * lc + 1;

            const double  weight = task.diagonal ? 1.0 : 2.0;
            double *const E      = thread_field[tid].data();

            for (const auto &B : ket.spans)
            {
                const int         k_ket = B.cgto_end - B.cgto_begin;
                const std::size_t cdim  = static_cast<std::size_t>(k_bra) * static_cast<std::size_t>(k_ket);
                dblk.resize(static_cast<std::size_t>(bra_comps) * ket_comps * cdim);

                // gather the density sub-block (out_row = ca·ket_comps + cc) and Σ|D|
                double dsum = 0.0;
                for (int ca = 0; ca < bra_comps; ca++)
                {
                    for (int cc = 0; cc < ket_comps; cc++)
                    {
                        const std::size_t out_row = static_cast<std::size_t>(ca) * ket_comps + cc;
                        for (int il = 0; il < k_bra; il++)
                        {
                            const std::size_t ao_a = static_cast<std::size_t>(ca) * bra.orb_indices[0] +
                                                     bra.orb_indices[static_cast<std::size_t>(A.cgto_begin + il) + 1];
                            for (int jl = 0; jl < k_ket; jl++)
                            {
                                const std::size_t ao_c = static_cast<std::size_t>(cc) * ket.orb_indices[0] +
                                                         ket.orb_indices[static_cast<std::size_t>(B.cgto_begin + jl) + 1];
                                const double v = d[ao_a * n_ao + ao_c];
                                dblk[out_row * cdim + static_cast<std::size_t>(il) * k_ket + jl] = v;
                                dsum += std::abs(v);
                            }
                        }
                    }
                }

                // density-weighted shell-pair screen: skip (A,B) whose field bound
                // (largest field element scaled by Σ|D_ac|) is negligible for any
                // point. Same-atom pairs are exempt — the estimate's R_AB^m factor
                // degenerates at R_AB = 0, and they carry the dominant near-field.
                if (screen && A.atom != B.atom)
                {
                    const double dx = A.x - B.x, dy = A.y - B.y, dz = A.z - B.z;
                    if (dsum * estimate(la, lc, A, B, dx * dx + dy * dy + dz * dz, 1.0) < threshold) continue;
                }

                kernel(la, lc, bra.view(), A.cgto_begin, A.cgto_end, ket.view(), B.cgto_begin, B.cgto_end, point_x,
                       point_y, point_z, n_points, dblk.data(), weight, E);
            }

            thread_busy[tid] += seconds(std::chrono::steady_clock::now() - t_body);
            thread_pairs[tid] += 1;
        }
    }

    balance.wall  = seconds(std::chrono::steady_clock::now() - t_region);
    balance.busy  = thread_busy;
    balance.pairs = thread_pairs;

    std::vector<std::array<double, 3>> result(np, {0.0, 0.0, 0.0});
    for (const auto& tf : thread_field)
    {
        for (std::size_t pt = 0; pt < np; pt++)
        {
            for (int a = 0; a < 3; a++) result[pt][static_cast<std::size_t>(a)] += tf[static_cast<std::size_t>(a) * np + pt];
        }
    }

    return result;
}

auto
external_center_field_compute_sparse(const CMolecule&         molecule,
                                     const CMolecularBasis&   basis,
                                     const BlockSparseMatrix& density,
                                     const FieldKernelFn      kernel,
                                     const ExternalEstimateFn estimate,
                                     const double*            point_x,
                                     const double*            point_y,
                                     const double*            point_z,
                                     const int                n_points,
                                     const double             threshold,
                                     ThreadBalance&           balance) -> std::vector<std::array<double, 3>>
{
    const auto gto_blocks = gtofunc::make_gto_blocks(basis, molecule);

    std::vector<BlockArrays> blocks;
    blocks.reserve(gto_blocks.size());
    for (const auto& block : gto_blocks)
    {
        blocks.push_back(fetch_block(block));
    }

    const auto tasks = build_tasks(blocks);

    // index the sparse density: the stored block per atom-pair (ga ≥ gb), and
    // each global AO's local position within its atom group
    const std::size_t n_groups = density.number_of_groups();
    std::vector<long> block_of(n_groups * n_groups, -1);
    for (std::size_t b = 0; b < density.number_of_blocks(); b++)
    {
        const auto desc = density.block(b);
        block_of[desc.groupA * n_groups + desc.groupB] = static_cast<long>(b);
    }
    std::vector<std::size_t> ao_local(density.dimension(), 0);
    for (std::size_t g = 0; g < n_groups; g++)
    {
        const auto aos = density.group_global_ao(g);
        for (std::size_t l = 0; l < aos.size(); l++) ao_local[aos[l]] = l;
    }

    const auto        nthreads = static_cast<std::size_t>(omp_get_max_threads());
    const std::size_t np       = static_cast<std::size_t>(n_points);
    const bool        screen   = threshold > 0.0;

    std::vector<std::vector<double>> thread_field(nthreads, std::vector<double>(3 * np, 0.0));
    std::vector<double>              thread_busy(nthreads, 0.0);
    std::vector<long>                thread_pairs(nthreads, 0);

    const double* const dvals = density.values();

    const auto t_region = std::chrono::steady_clock::now();

#pragma omp parallel
    {
        std::vector<double> dblk;

#pragma omp for schedule(dynamic)
        for (int p = 0; p < static_cast<int>(tasks.size()); p++)
        {
            const auto tid    = static_cast<std::size_t>(omp_get_thread_num());
            const auto t_body = std::chrono::steady_clock::now();

            const auto &task = tasks[static_cast<std::size_t>(p)];
            const auto &bra  = blocks[task.i];
            const auto &A    = bra.spans[task.span];
            const auto &ket  = blocks[task.j];

            const int la        = bra.angular_momentum;
            const int lc        = ket.angular_momentum;
            const int k_bra     = A.cgto_end - A.cgto_begin;
            const int bra_comps = 2 * la + 1;
            const int ket_comps = 2 * lc + 1;

            const double  weight = task.diagonal ? 1.0 : 2.0;
            double *const E      = thread_field[tid].data();

            for (const auto &B : ket.spans)
            {
                const auto        atom_a = static_cast<std::size_t>(A.atom);
                const auto        atom_b = static_cast<std::size_t>(B.atom);
                const std::size_t ga     = std::max(atom_a, atom_b);
                const std::size_t gb     = std::min(atom_a, atom_b);
                const long        blk    = block_of[ga * n_groups + gb];
                if (blk < 0) continue;  // density block not stored — atom-pair negligible

                const auto          desc  = density.block(static_cast<std::size_t>(blk));
                const double* const dblock = dvals + desc.offset;
                const std::size_t   ncol  = desc.columnCount;
                const bool          ahi   = atom_a >= atom_b;  // A is the row group (ga)

                const int         k_ket = B.cgto_end - B.cgto_begin;
                const std::size_t cdim  = static_cast<std::size_t>(k_bra) * static_cast<std::size_t>(k_ket);
                dblk.resize(static_cast<std::size_t>(bra_comps) * ket_comps * cdim);

                // gather the density sub-block from the sparse block (transpose
                // local indices when A is the lower-numbered atom), and Σ|D|
                double dsum = 0.0;
                for (int ca = 0; ca < bra_comps; ca++)
                {
                    for (int cc = 0; cc < ket_comps; cc++)
                    {
                        const std::size_t out_row = static_cast<std::size_t>(ca) * ket_comps + cc;
                        for (int il = 0; il < k_bra; il++)
                        {
                            const std::size_t ao_a = static_cast<std::size_t>(ca) * bra.orb_indices[0] +
                                                     bra.orb_indices[static_cast<std::size_t>(A.cgto_begin + il) + 1];
                            const std::size_t la_loc = ao_local[ao_a];
                            for (int jl = 0; jl < k_ket; jl++)
                            {
                                const std::size_t ao_c = static_cast<std::size_t>(cc) * ket.orb_indices[0] +
                                                         ket.orb_indices[static_cast<std::size_t>(B.cgto_begin + jl) + 1];
                                const std::size_t lc_loc = ao_local[ao_c];
                                const double      v = ahi ? dblock[la_loc * ncol + lc_loc] : dblock[lc_loc * ncol + la_loc];
                                dblk[out_row * cdim + static_cast<std::size_t>(il) * k_ket + jl] = v;
                                dsum += std::abs(v);
                            }
                        }
                    }
                }

                if (screen && A.atom != B.atom)
                {
                    const double dx = A.x - B.x, dy = A.y - B.y, dz = A.z - B.z;
                    if (dsum * estimate(la, lc, A, B, dx * dx + dy * dy + dz * dz, 1.0) < threshold) continue;
                }

                kernel(la, lc, bra.view(), A.cgto_begin, A.cgto_end, ket.view(), B.cgto_begin, B.cgto_end, point_x,
                       point_y, point_z, n_points, dblk.data(), weight, E);
            }

            thread_busy[tid] += seconds(std::chrono::steady_clock::now() - t_body);
            thread_pairs[tid] += 1;
        }
    }

    balance.wall  = seconds(std::chrono::steady_clock::now() - t_region);
    balance.busy  = thread_busy;
    balance.pairs = thread_pairs;

    std::vector<std::array<double, 3>> result(np, {0.0, 0.0, 0.0});
    for (const auto& tf : thread_field)
    {
        for (std::size_t pt = 0; pt < np; pt++)
        {
            for (int a = 0; a < 3; a++) result[pt][static_cast<std::size_t>(a)] += tf[static_cast<std::size_t>(a) * np + pt];
        }
    }

    return result;
}

}  // namespace tabula
