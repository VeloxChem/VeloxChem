//
//  Tabula — custom-recursion molecular-integral machinery.
//  Two-center overlap screening predicate.
//

#ifndef TabulaOverlapScreener_hpp
#define TabulaOverlapScreener_hpp

#include <functional>

#include "GtoBlockScreeningData.hpp"

namespace tabula {  // tabula namespace

/// @brief Builds the overlap shell-pair screening estimator for an
/// `(l_a, l_c)` basis-function-block pair.
///
/// The returned callable is the `estimator` argument of the screened
/// `CGtoPairBlock` constructor: given the screening data of a bra and a ket
/// contracted GTO and the distance `|R|` between their centers, it returns
/// the overlap screening estimate (see `docs/screening-predicates.md`). A
/// contracted-GTO pair is kept when the estimate is at or above the
/// screening threshold.
///
/// Same-center pairs (`|R| = 0`) are always kept — the estimate's
/// `r^(l_a+l_c)` factor is an asymptotic, large-`R` bound and is not valid
/// at `R = 0`.
///
/// @param l_a The angular momentum of the bra block.
/// @param l_c The angular momentum of the ket block.
/// @return The overlap screening estimator.
auto make_overlap_screening_estimator(const int l_a, const int l_c)
    -> std::function<double(const CGtoBlockScreeningData &, const CGtoBlockScreeningData &, const double)>;

}  // namespace tabula

#endif /* TabulaOverlapScreener_hpp */
