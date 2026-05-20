//
//  Tabula — custom-recursion molecular-integral machinery.
//  Shared external-center driver core — operator-agnostic scaffolding for
//  integrals over two Gaussian centers and a set of external point charges
//  (nuclear attraction; later, charge–dipole). Sibling of the two-center
//  core; differs only in that the kernel call carries the charge set.
//

#ifndef TabulaExternalCenterDriver_hpp
#define TabulaExternalCenterDriver_hpp

#include "MolecularBasis.hpp"
#include "Molecule.hpp"
#include "TabulaBlockSparseMatrix.hpp"
#include "TabulaChargeSet.hpp"
#include "TabulaDenseMatrix.hpp"
#include "TabulaKernelBlockData.hpp"
#include "TabulaThreadBalance.hpp"
#include "TabulaTwoCenterDriver.hpp"  // AtomSpan, KernelProfile, EstimateFn

namespace tabula {  // tabula namespace

/// @brief A fused external-center kernel dispatch — like `KernelFn` but with
/// the point-charge set the kernel sums over.
using ExternalKernelFn =
    void (*)(int, int, const KernelBlockData&, int, int, const KernelBlockData&, const ChargeSet&, double*);

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
/// @param op The operator's fused kernel and screening estimate.
/// @param charges The point charges to sum over.
/// @param profile Optional — receives the per-phase wall-time breakdown.
/// @param balance Receives the task loop's per-thread load balance.
/// @return The dense integral matrix.
auto external_center_compute(const CMolecule&              molecule,
                             const CMolecularBasis&        basis,
                             const double                  threshold,
                             const ExternalCenterOperator& op,
                             const ChargeSet&              charges,
                             KernelProfile*                profile,
                             ThreadBalance&                balance) -> DenseMatrix;

/// @brief Computes an external-center integral matrix for `op` in block-sparse
/// storage — the shared sparse driver.
auto external_center_compute_sparse(const CMolecule&              molecule,
                                    const CMolecularBasis&        basis,
                                    const double                  threshold,
                                    const ExternalCenterOperator& op,
                                    const ChargeSet&              charges,
                                    KernelProfile*                profile,
                                    ThreadBalance&                balance) -> BlockSparseMatrix;

}  // namespace tabula

#endif /* TabulaExternalCenterDriver_hpp */
