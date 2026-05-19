//
//  Tabula — custom-recursion molecular-integral machinery.
//  Two-center kinetic-energy driver.
//

#ifndef TabulaKineticDriver_hpp
#define TabulaKineticDriver_hpp

#include <vector>

#include "MolecularBasis.hpp"
#include "Molecule.hpp"
#include "TabulaBlockSparseMatrix.hpp"
#include "TabulaDenseMatrix.hpp"
#include "TabulaThreadBalance.hpp"

namespace tabula {  // tabula namespace

/// @brief Per-phase wall-time breakdown of a `KineticDriver::compute` run,
/// in seconds, summed over the OpenMP threads.
struct KineticProfile
{
    /// @brief `gtofunc::make_gto_blocks`.
    double make_blocks{0.0};
    /// @brief Per-block primitive-data fetch.
    double pair_setup{0.0};
    /// @brief Atom-span screening and the gather of the surviving ket block.
    double screen{0.0};
    /// @brief The fused kinetic kernel — the primitive-pair weight, the seed
    /// ladder, the primitive contraction, the single-centre MD recursion, and
    /// the Cartesian-to-spherical assembly.
    double kernel{0.0};
    /// @brief The scatter into the matrix.
    double scatter{0.0};
    /// @brief Unused — kept for layout parity with the overlap profile.
    double symmetrize{0.0};
};

/// @brief Driver for the two-center kinetic-energy integral `−½ ∫ χ_a ∇² χ_c`.
///
/// Mirrors `OverlapDriver` — intakes a VeloxChem molecule and molecular basis,
/// builds the basis-function blocks, and for each task (one block pair
/// restricted to a bra contracted-GTO range) runs the fused kinetic kernel
/// directly from the block data and scatters the result into the matrix.
class KineticDriver
{
   public:
    /// @brief Creates a kinetic-energy driver.
    KineticDriver() = default;

    /// @brief Computes the kinetic-energy matrix.
    /// @param molecule The molecule.
    /// @param basis The molecular basis.
    /// @param threshold The screening threshold — a ket atom span is kept
    /// when its conservative kinetic estimate against the bra span is at or
    /// above this value. `0` (the default) keeps every span.
    /// @param profile Optional — when non-null, receives the per-phase
    /// wall-time breakdown of the run.
    /// @return The kinetic-energy matrix.
    auto compute(const CMolecule&       molecule,
                 const CMolecularBasis& basis,
                 const double           threshold = 0.0,
                 KineticProfile        *profile   = nullptr) const -> DenseMatrix;

    /// @brief Computes the kinetic-energy matrix in block-sparse storage.
    ///
    /// Only the atom-pair blocks whose conservative kinetic estimate clears
    /// the threshold are stored; the rest of the matrix is dropped.
    /// @param molecule The molecule.
    /// @param basis The molecular basis.
    /// @param threshold The screening threshold.
    /// @param profile Optional — when non-null, receives the per-phase
    /// wall-time breakdown of the run.
    /// @return The block-sparse kinetic-energy matrix.
    auto computeSparse(const CMolecule&       molecule,
                       const CMolecularBasis& basis,
                       const double           threshold = 0.0,
                       KineticProfile        *profile   = nullptr) const -> BlockSparseMatrix;
};

/// @brief Gets the per-thread balance captured by the most recent
/// `KineticDriver` compute run.
/// @return The thread-balance capture.
auto kinetic_thread_balance() -> ThreadBalance;

}  // namespace tabula

#endif /* TabulaKineticDriver_hpp */
