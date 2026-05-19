//
//  Tabula — custom-recursion molecular-integral machinery.
//  Two-center overlap driver.
//

#ifndef TabulaOverlapDriver_hpp
#define TabulaOverlapDriver_hpp

#include "MolecularBasis.hpp"
#include "Molecule.hpp"
#include "TabulaBlockSparseMatrix.hpp"
#include "TabulaDenseMatrix.hpp"
#include "TabulaThreadBalance.hpp"
#include "TabulaTwoCenterDriver.hpp"

namespace tabula {  // tabula namespace

/// @brief Driver for the two-center overlap integral `∫ χ_a χ_c`. A thin
/// wrapper over the shared two-center core — supplies the overlap kernel and
/// the overlap screening estimate.
class OverlapDriver
{
   public:
    /// @brief Creates an overlap driver.
    OverlapDriver() = default;

    /// @brief Computes the overlap matrix.
    /// @param molecule The molecule.
    /// @param basis The molecular basis.
    /// @param threshold The screening threshold; `0` (the default) keeps
    /// every atom span.
    /// @param profile Optional — when non-null, receives the per-phase
    /// wall-time breakdown of the run.
    /// @return The overlap matrix.
    auto compute(const CMolecule&       molecule,
                 const CMolecularBasis& basis,
                 const double           threshold = 0.0,
                 KernelProfile         *profile   = nullptr) const -> DenseMatrix;

    /// @brief Computes the overlap matrix in block-sparse storage.
    /// @param molecule The molecule.
    /// @param basis The molecular basis.
    /// @param threshold The screening threshold; `0` keeps every atom pair.
    /// @param profile Optional — when non-null, receives the per-phase
    /// wall-time breakdown of the run.
    /// @return The block-sparse overlap matrix.
    auto computeSparse(const CMolecule&       molecule,
                       const CMolecularBasis& basis,
                       const double           threshold = 0.0,
                       KernelProfile         *profile   = nullptr) const -> BlockSparseMatrix;
};

/// @brief Gets the per-thread balance captured by the most recent
/// `OverlapDriver` compute run.
/// @return The thread-balance capture.
auto overlap_thread_balance() -> ThreadBalance;

}  // namespace tabula

#endif /* TabulaOverlapDriver_hpp */
