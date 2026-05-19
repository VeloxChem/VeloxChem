//
//  Tabula — custom-recursion molecular-integral machinery.
//  Two-center overlap driver.
//

#ifndef TabulaOverlapDriver_hpp
#define TabulaOverlapDriver_hpp

#include "MolecularBasis.hpp"
#include "Molecule.hpp"
#include "TabulaDenseMatrix.hpp"

namespace tabula {  // tabula namespace

/// @brief Per-phase wall-time breakdown of an `OverlapDriver::compute` run,
/// in seconds, summed over the OpenMP threads.
struct OverlapProfile
{
    /// @brief `gtofunc::make_gto_blocks`.
    double make_blocks{0.0};
    /// @brief Screened `CGtoPairBlock` construction.
    double pair_setup{0.0};
    /// @brief Step (a) — the seed ladder `[0]^m`.
    double seed{0.0};
    /// @brief Step (b) — the primitive-pair contraction.
    double contract{0.0};
    /// @brief Step (c) — the single-centre MD recursion.
    double md{0.0};
    /// @brief Step (d) — the Cartesian-to-spherical assembly.
    double transform{0.0};
    /// @brief Step (e) — the scatter into the matrix.
    double scatter{0.0};
};

/// @brief Driver for the two-center overlap integral.
///
/// Intakes a VeloxChem molecule and molecular basis, builds the basis-function
/// blocks with `gtofunc::make_gto_blocks`, and for each block pair constructs
/// a screened `CGtoPairBlock` (the overlap screening estimator drops the
/// negligible contracted-GTO pairs). The overlap matrix is then evaluated on
/// the late-contraction recursion path — seed ladder, primitive-pair
/// contraction, single-centre MD recursion, Cartesian-to-spherical assembly,
/// and the scatter into the matrix.
class OverlapDriver
{
   public:
    /// @brief Creates an overlap driver.
    OverlapDriver() = default;

    /// @brief Computes the overlap matrix.
    /// @param molecule The molecule.
    /// @param basis The molecular basis.
    /// @param threshold The screening threshold — a contracted-GTO pair is
    /// kept when its overlap screening estimate is at or above this value.
    /// `0` (the default) keeps every pair.
    /// @param profile Optional — when non-null, receives the per-phase
    /// wall-time breakdown of the run.
    /// @return The overlap matrix.
    auto compute(const CMolecule&     molecule,
                 const CMolecularBasis& basis,
                 const double         threshold = 0.0,
                 OverlapProfile      *profile   = nullptr) const -> DenseMatrix;
};

/// @brief Accumulated wall time of the overlap scatter, in seconds, summed
/// over the OpenMP threads — split into the direct `matrix(r,c)` writes and
/// the transposed `matrix(c,r)` writes an off-diagonal block pair adds.
struct ScatterProfile
{
    /// @brief The direct `matrix(r,c)` writes.
    double direct{0.0};
    /// @brief The transposed `matrix(c,r)` writes.
    double transpose{0.0};
};

/// @brief Gets the accumulated overlap-scatter profile.
/// @return The scatter profile.
auto scatter_profile() -> ScatterProfile;

/// @brief Resets the accumulated overlap-scatter profile.
auto reset_scatter_profile() -> void;

}  // namespace tabula

#endif /* TabulaOverlapDriver_hpp */
