//
//  Tabula — custom-recursion molecular-integral machinery.
//  Two-center kinetic-energy driver.
//

#ifndef TabulaKineticDriver_hpp
#define TabulaKineticDriver_hpp

#include "MolecularBasis.hpp"
#include "Molecule.hpp"
#include "TabulaBlockSparseMatrix.hpp"
#include "TabulaDenseMatrix.hpp"
#include "TabulaThreadBalance.hpp"
#include "TabulaTwoCenterDriver.hpp"

namespace tabula {  // tabula namespace

/// @brief Driver for the two-center kinetic-energy integral `−½∫χ_a∇²χ_c`. A
/// thin wrapper over the shared two-center core — supplies the kinetic kernel
/// and the kinetic screening estimate.
class KineticDriver
{
   public:
    /// @brief Creates a kinetic-energy driver.
    KineticDriver() = default;

    /// @brief Computes the kinetic-energy matrix.
    /// @param molecule The molecule.
    /// @param basis The molecular basis.
    /// @param threshold The screening threshold; `0` (the default) keeps
    /// every atom span.
    /// @param profile Optional — when non-null, receives the per-phase
    /// wall-time breakdown of the run.
    /// @return The kinetic-energy matrix.
    auto compute(const CMolecule&       molecule,
                 const CMolecularBasis& basis,
                 const double           threshold = 0.0,
                 KernelProfile         *profile   = nullptr) const -> DenseMatrix;

    /// @brief Computes the kinetic-energy matrix in block-sparse storage.
    /// @param molecule The molecule.
    /// @param basis The molecular basis.
    /// @param threshold The screening threshold; `0` keeps every atom pair.
    /// @param profile Optional — when non-null, receives the per-phase
    /// wall-time breakdown of the run.
    /// @return The block-sparse kinetic-energy matrix.
    auto computeSparse(const CMolecule&       molecule,
                       const CMolecularBasis& basis,
                       const double           threshold = 0.0,
                       KernelProfile         *profile   = nullptr) const -> BlockSparseMatrix;
};

/// @brief Gets the per-thread balance captured by the most recent
/// `KineticDriver` compute run.
/// @return The thread-balance capture.
auto kinetic_thread_balance() -> ThreadBalance;

}  // namespace tabula

#endif /* TabulaKineticDriver_hpp */
