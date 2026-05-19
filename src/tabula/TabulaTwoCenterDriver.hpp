//
//  Tabula — custom-recursion molecular-integral machinery.
//  Shared two-center driver core — operator-agnostic scaffolding.
//

#ifndef TabulaTwoCenterDriver_hpp
#define TabulaTwoCenterDriver_hpp

#include "MolecularBasis.hpp"
#include "Molecule.hpp"
#include "TabulaBlockSparseMatrix.hpp"
#include "TabulaDenseMatrix.hpp"
#include "TabulaKernelBlockData.hpp"
#include "TabulaThreadBalance.hpp"

namespace tabula {  // tabula namespace

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

/// @brief Per-phase wall-time breakdown of a two-center compute run, in
/// seconds, summed over the OpenMP threads.
struct KernelProfile
{
    /// @brief `gtofunc::make_gto_blocks`.
    double make_blocks{0.0};
    /// @brief Per-block primitive-data fetch.
    double pair_setup{0.0};
    /// @brief Atom-span screening and the gather of the surviving ket block.
    double screen{0.0};
    /// @brief The fused kernel.
    double kernel{0.0};
    /// @brief The scatter into the matrix.
    double scatter{0.0};
    /// @brief Unused — kept for profiler-dict layout parity.
    double symmetrize{0.0};
};

/// @brief A fused two-center kernel dispatch — `op_kernel(l_a, l_c, …)`.
using KernelFn = void (*)(int, int, const KernelBlockData&, int, int, const KernelBlockData&, double*);

/// @brief A conservative upper bound on an operator's integral magnitude for
/// the shell-pair sub-block of bra span `A` and ket span `B` separated by
/// `r²` — the screening predicate.
using EstimateFn = double (*)(int, int, const AtomSpan&, const AtomSpan&, double);

/// @brief The operator-specific halves of a two-center driver — the fused
/// kernel and the screening estimate. Everything else is shared.
struct TwoCenterOperator
{
    /// @brief The fused `(l_a, l_c)` kernel dispatch.
    KernelFn kernel;
    /// @brief The screening upper-bound estimate.
    EstimateFn estimate;
};

/// @brief Computes a two-center integral matrix for `op` — the shared dense
/// driver behind `OverlapDriver` / `KineticDriver`.
/// @param molecule The molecule.
/// @param basis The molecular basis.
/// @param threshold The screening threshold; `0` keeps every atom span.
/// @param op The operator's fused kernel and screening estimate.
/// @param profile Optional — receives the per-phase wall-time breakdown.
/// @param balance Receives the task loop's per-thread load balance.
/// @return The dense integral matrix.
auto two_center_compute(const CMolecule&         molecule,
                        const CMolecularBasis&   basis,
                        const double             threshold,
                        const TwoCenterOperator& op,
                        KernelProfile*           profile,
                        ThreadBalance&           balance) -> DenseMatrix;

/// @brief Computes a two-center integral matrix for `op` in block-sparse
/// storage — the shared sparse driver.
/// @param molecule The molecule.
/// @param basis The molecular basis.
/// @param threshold The screening threshold; `0` keeps every atom pair.
/// @param op The operator's fused kernel and screening estimate.
/// @param profile Optional — receives the per-phase wall-time breakdown.
/// @param balance Receives the task loop's per-thread load balance.
/// @return The block-sparse integral matrix.
auto two_center_compute_sparse(const CMolecule&         molecule,
                               const CMolecularBasis&   basis,
                               const double             threshold,
                               const TwoCenterOperator& op,
                               KernelProfile*           profile,
                               ThreadBalance&           balance) -> BlockSparseMatrix;

}  // namespace tabula

#endif /* TabulaTwoCenterDriver_hpp */
