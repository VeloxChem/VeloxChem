//
//  Tabula — custom-recursion molecular-integral machinery.
//  Boys function — the Coulomb / nuclear-attraction integral seed.
//

#ifndef TabulaBoys_hpp
#define TabulaBoys_hpp

namespace tabula {  // tabula namespace

/// @brief The largest Boys-function order the evaluator covers.
inline constexpr int boys_max_order = 32;

/// @brief Evaluates the Boys function `F_m(x) = ∫₀¹ t^(2m)·e^(−x·t²) dt` for
/// every order `m = 0 … order`, writing `F_m` to `results[m]`.
///
/// A three-region rational-minimax approximation (Vikhamar-Sandberg &
/// Repisky, arXiv:2512.10059v3): a per-order minimax then a downward
/// recursion on `[0, 11.90)`, an `F_0` minimax then an upward recursion on
/// `[11.90, 28.99)`, and the asymptotic `F_0` then an upward recursion above.
/// Guaranteed ~5e-14 absolute error over `m = 0 … 32`.
///
/// @param order The highest order; `0 … boys_max_order`.
/// @param x The argument; `x >= 0`.
/// @param results The output buffer — at least `order + 1` doubles.
auto boys(const int order, const double x, double *results) -> void;

}  // namespace tabula

#endif /* TabulaBoys_hpp */
