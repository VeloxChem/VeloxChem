//
//  Tabula — custom-recursion molecular-integral machinery.
//  Two-center overlap recursion.
//

#ifndef TabulaOverlapRecursion_hpp
#define TabulaOverlapRecursion_hpp

#include <vector>

#include "TabulaGtoPairBlock.hpp"

namespace tabula {  // tabula namespace

/// @brief Per-section wall time of `compute_overlap_seed`, in seconds,
/// accumulated across calls and summed over the OpenMP threads.
struct SeedProfile
{
    /// @brief Allocating and zero-filling the seed buffer.
    double allocate{0.0};
    /// @brief Building row 0 — `[0]^0` and its `(1/2α)^l_a · (−1/2γ)^l_c`
    /// power passes.
    double row0{0.0};
    /// @brief Building the seed ladder rows `1 … L`.
    double ladder{0.0};
};

/// @brief Gets the accumulated `compute_overlap_seed` profile.
/// @return The seed profile.
auto seed_profile() -> SeedProfile;

/// @brief Resets the accumulated `compute_overlap_seed` profile.
auto reset_seed_profile() -> void;

/// @brief Computes the overlap seed ladder `[0]^m`, `m = 0 … l_a+l_c`, for
/// every primitive pair of a basis-function-pair block — step (a) of the
/// late-contraction overlap recursion.
///
/// The result is a row-major buffer of `L+1` rows (`L = l_a + l_c`), one per
/// seed order `m`. A row holds one value per primitive pair, in the pair
/// block's primitive-pair-major order (`ijoff = pp·cdim + ij`); the row
/// stride is the primitive-pair count padded up to a multiple of 8, so each
/// row is aligned for SIMD.
///
/// `[0]^0 = weight · (1/2α)^l_a · (−1/2γ)^l_c`, reusing the pair block's
/// precomputed weights (the primitive normalization and overlap factors
/// folded into one, the contraction weight `c_i·c_j` premultiplied in) and
/// primitive exponents — so the contraction step is a plain sum. The ladder
/// is `[0]^m = (−2ρ)·[0]^(m−1)`, `ρ = αγ/(α+γ)`.
///
/// @param pair_block The basis-function-pair block.
/// @return The seed-ladder buffer.
auto compute_overlap_seed(const GtoPairBlock& pair_block) -> std::vector<double>;

}  // namespace tabula

#endif /* TabulaOverlapRecursion_hpp */
