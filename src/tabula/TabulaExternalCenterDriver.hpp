//
//  Tabula — custom-recursion molecular-integral machinery.
//  Shared external-center driver core — operator-agnostic scaffolding for
//  integrals over two Gaussian centers and a set of external point sources
//  (nuclear attraction's charges; charge–dipole's dipoles). Sibling of the
//  two-center core; the kernel closure carries its own source set, so the
//  scaffolding never names it.
//

#ifndef TabulaExternalCenterDriver_hpp
#define TabulaExternalCenterDriver_hpp

#include <array>
#include <functional>
#include <vector>

#include "MolecularBasis.hpp"
#include "Molecule.hpp"
#include "TabulaBlockSparseMatrix.hpp"
#include "TabulaDenseMatrix.hpp"
#include "TabulaKernelBlockData.hpp"
#include "TabulaThreadBalance.hpp"
#include "TabulaTwoCenterDriver.hpp"  // AtomSpan, KernelProfile, EstimateFn

namespace tabula {  // tabula namespace

/// @brief A fused external-center kernel call — the `(l_a, l_c)` dispatch with
/// its point-source set already bound (the closure captures the charges or
/// dipoles), so the shared scaffolding stays source-agnostic.
using ExternalKernelFn =
    std::function<void(int, int, const KernelBlockData&, int, int, const KernelBlockData&, double*)>;

/// @brief A conservative upper bound on an external-center operator's block
/// magnitude — like `EstimateFn` but with the charge factor `Σ|Z_N|` (the
/// charge sum that scales the bound; the point-source kernel is bounded
/// charge-position-independently).
using ExternalEstimateFn = double (*)(int, int, const AtomSpan&, const AtomSpan&, double, double);

/// @brief The operator-specific halves of an external-center driver — the
/// fused kernel, the screening estimate, and the charge factor the estimate
/// scales by. Everything else is shared.
struct ExternalCenterOperator
{
    /// @brief The fused `(l_a, l_c)` kernel dispatch.
    ExternalKernelFn kernel;
    /// @brief The screening upper-bound estimate.
    ExternalEstimateFn estimate;
    /// @brief `Σ_N |Z_N|` — the charge factor the estimate scales by.
    double charge_factor;
};

/// @brief Computes an external-center integral matrix for `op` — the shared
/// dense driver behind `NuclearAttractionDriver`.
/// @param molecule The molecule.
/// @param basis The molecular basis.
/// @param threshold The screening threshold; `0` keeps every atom span.
/// @param op The operator's fused kernel (source-bound) and screening estimate.
/// @param profile Optional — receives the per-phase wall-time breakdown.
/// @param balance Receives the task loop's per-thread load balance.
/// @return The dense integral matrix.
auto external_center_compute(const CMolecule&              molecule,
                             const CMolecularBasis&        basis,
                             const double                  threshold,
                             const ExternalCenterOperator& op,
                             KernelProfile*                profile,
                             ThreadBalance&                balance) -> DenseMatrix;

/// @brief Computes an external-center integral matrix for `op` in block-sparse
/// storage — the shared sparse driver.
auto external_center_compute_sparse(const CMolecule&              molecule,
                                    const CMolecularBasis&        basis,
                                    const double                  threshold,
                                    const ExternalCenterOperator& op,
                                    KernelProfile*                profile,
                                    ThreadBalance&                balance) -> BlockSparseMatrix;

/// @brief A field-contraction kernel call — one (bra-block, ket-block) pair's
/// contribution to the per-point field, given the gathered density block, the
/// block-pair multiplicity, and the field accumulator. The external points are
/// bound into the closure.
using FieldKernelFn =
    std::function<void(int, int, const KernelBlockData&, int, int, const KernelBlockData&, const double*, double, double*)>;

/// @brief Contracts an external-center field operator against a density matrix
/// to the field at each external point — the transpose of the matrix path. It
/// sweeps the atom-block-pair triangle (off-diagonal doubled, since the field
/// and density are symmetric), gathers each pair's density sub-block, runs the
/// field `kernel`, and reduces the per-thread per-point accumulators.
/// @param molecule The molecule.
/// @param basis The molecular basis.
/// @param density The AO density matrix.
/// @param kernel The field-contraction kernel (its points bound).
/// @param n_points The number of external points.
/// @param balance Receives the task loop's per-thread load balance.
/// @return The field at each point — `n_points` entries of `{Ex, Ey, Ez}`.
auto external_center_field_compute(const CMolecule&       molecule,
                                   const CMolecularBasis& basis,
                                   const DenseMatrix&     density,
                                   const FieldKernelFn&   kernel,
                                   int                    n_points,
                                   ThreadBalance&         balance) -> std::vector<std::array<double, 3>>;

}  // namespace tabula

#endif /* TabulaExternalCenterDriver_hpp */
