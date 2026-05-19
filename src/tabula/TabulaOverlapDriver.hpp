//
//  Tabula — custom-recursion molecular-integral machinery.
//  Two-center overlap driver.
//

#ifndef TabulaOverlapDriver_hpp
#define TabulaOverlapDriver_hpp

#include <vector>

#include "MolecularBasis.hpp"
#include "Molecule.hpp"
#include "TabulaBlockSparseMatrix.hpp"
#include "TabulaDenseMatrix.hpp"

namespace tabula {  // tabula namespace

/// @brief Per-phase wall-time breakdown of an `OverlapDriver::compute` run,
/// in seconds, summed over the OpenMP threads.
struct OverlapProfile
{
    /// @brief `gtofunc::make_gto_blocks`.
    double make_blocks{0.0};
    /// @brief Per-block primitive-data fetch.
    double pair_setup{0.0};
    /// @brief Atom-span screening and the gather of the surviving ket block.
    double screen{0.0};
    /// @brief The fused overlap kernel — the primitive-pair weight, the seed
    /// ladder, the primitive contraction, the single-centre MD recursion, and
    /// the Cartesian-to-spherical assembly.
    double kernel{0.0};
    /// @brief The scatter into the matrix's upper triangle.
    double scatter{0.0};
    /// @brief Mirroring the upper triangle into the lower — one `symmetrize`
    /// after every block pair (wall time, not thread-summed).
    double symmetrize{0.0};
};

/// @brief Driver for the two-center overlap integral.
///
/// Intakes a VeloxChem molecule and molecular basis, builds the basis-function
/// blocks with `gtofunc::make_gto_blocks`, and for each task — one block pair
/// restricted to a bra contracted-GTO range — runs the fused overlap kernel
/// directly from the block data and scatters the result into the matrix.
class OverlapDriver
{
   public:
    /// @brief Creates an overlap driver.
    OverlapDriver() = default;

    /// @brief Computes the overlap matrix.
    /// @param molecule The molecule.
    /// @param basis The molecular basis.
    /// @param threshold The screening threshold — a ket atom span is kept
    /// when its conservative overlap estimate against the bra span is at or
    /// above this value. `0` (the default) keeps every span.
    /// @param profile Optional — when non-null, receives the per-phase
    /// wall-time breakdown of the run.
    /// @return The overlap matrix.
    auto compute(const CMolecule&     molecule,
                 const CMolecularBasis& basis,
                 const double         threshold = 0.0,
                 OverlapProfile      *profile   = nullptr) const -> DenseMatrix;

    /// @brief Computes the overlap matrix in block-sparse storage.
    ///
    /// Only the atom-pair blocks whose conservative overlap estimate clears
    /// the threshold are stored; the rest of the matrix is dropped, giving an
    /// O(N) footprint for an extended molecule. `0` keeps every atom pair.
    /// @param molecule The molecule.
    /// @param basis The molecular basis.
    /// @param threshold The screening threshold.
    /// @param profile Optional — when non-null, receives the per-phase
    /// wall-time breakdown of the run.
    /// @return The block-sparse overlap matrix.
    auto computeSparse(const CMolecule&     molecule,
                       const CMolecularBasis& basis,
                       const double         threshold = 0.0,
                       OverlapProfile      *profile   = nullptr) const -> BlockSparseMatrix;
};

/// @brief Per-thread load balance of the task parallel loop in the most
/// recent `OverlapDriver::compute` run.
struct ThreadBalance
{
    /// @brief Wall time of the whole parallel region.
    double wall{0.0};
    /// @brief Busy wall seconds of each OpenMP thread.
    std::vector<double> busy;
    /// @brief Number of tasks handled by each OpenMP thread.
    std::vector<long> pairs;
};

/// @brief Gets the per-thread balance captured by the most recent
/// `OverlapDriver::compute` run.
/// @return The thread-balance capture.
auto overlap_thread_balance() -> ThreadBalance;

}  // namespace tabula

#endif /* TabulaOverlapDriver_hpp */
