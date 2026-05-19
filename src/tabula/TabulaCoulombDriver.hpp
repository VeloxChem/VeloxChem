//
//  Tabula — custom-recursion molecular-integral machinery.
//  Two-center Coulomb driver.
//

#ifndef TabulaCoulombDriver_hpp
#define TabulaCoulombDriver_hpp

#include "MolecularBasis.hpp"
#include "Molecule.hpp"
#include "TabulaBlockSparseMatrix.hpp"
#include "TabulaDenseMatrix.hpp"
#include "TabulaThreadBalance.hpp"
#include "TabulaTwoCenterDriver.hpp"

namespace tabula {  // tabula namespace

/// @brief Driver for the two-center Coulomb integral
/// `(a|b) = ∫∫ χ_a(r₁)·r₁₂⁻¹·χ_c(r₂)`. A thin wrapper over the shared
/// two-center core — supplies the Coulomb kernel and the Coulomb screening
/// estimate.
class CoulombDriver
{
   public:
    /// @brief Creates a Coulomb driver.
    CoulombDriver() = default;

    /// @brief Computes the Coulomb matrix.
    /// @param molecule The molecule.
    /// @param basis The molecular basis.
    /// @param threshold The screening threshold; `0` (the default) keeps
    /// every atom span.
    /// @param profile Optional — when non-null, receives the per-phase
    /// wall-time breakdown of the run.
    /// @return The Coulomb matrix.
    auto compute(const CMolecule&       molecule,
                 const CMolecularBasis& basis,
                 const double           threshold = 0.0,
                 KernelProfile         *profile   = nullptr) const -> DenseMatrix;

    /// @brief Computes the Coulomb matrix in block-sparse storage.
    /// @param molecule The molecule.
    /// @param basis The molecular basis.
    /// @param threshold The screening threshold; `0` keeps every atom pair.
    /// @param profile Optional — when non-null, receives the per-phase
    /// wall-time breakdown of the run.
    /// @return The block-sparse Coulomb matrix.
    auto computeSparse(const CMolecule&       molecule,
                       const CMolecularBasis& basis,
                       const double           threshold = 0.0,
                       KernelProfile         *profile   = nullptr) const -> BlockSparseMatrix;
};

/// @brief Gets the per-thread balance captured by the most recent
/// `CoulombDriver` compute run.
/// @return The thread-balance capture.
auto coulomb_thread_balance() -> ThreadBalance;

}  // namespace tabula

#endif /* TabulaCoulombDriver_hpp */
