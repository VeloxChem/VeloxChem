//
//  Tabula — custom-recursion molecular-integral machinery.
//  Single-centre McMurchie–Davidson recursion.
//

#ifndef TabulaMDRecursion_hpp
#define TabulaMDRecursion_hpp

#include <array>
#include <cstddef>
#include <vector>

namespace tabula {  // tabula namespace

/// @brief Per-section wall time of `compute_one_center_md`, in seconds,
/// accumulated across calls and summed over the OpenMP threads.
struct MDRecursionProfile
{
    /// @brief Zero-filling the per-degree buffers.
    double allocate{0.0};
    /// @brief Copying the contracted seed ladder into level 0.
    double seed_copy{0.0};
    /// @brief Building each degree level; `build_by_degree[d]` is the wall
    /// time spent on level `d` (the recursion fills degrees `1 … L`).
    std::array<double, 9> build_by_degree{};
};

/// @brief Gets the accumulated `compute_one_center_md` profile.
/// @return The recursion profile.
auto md_recursion_profile() -> MDRecursionProfile;

/// @brief Resets the accumulated `compute_one_center_md` profile.
auto reset_md_recursion_profile() -> void;

/// @brief The single-centre McMurchie–Davidson recursion — step (c) of the
/// late-contraction recursion, generic to any integral.
///
/// From the contracted seed ladder `[0]^m`, `m = 0…L`, it builds the
/// `R_AC`-monomial integrals `[r]^0` for `|r| = L` via
///
///     `[r]^(m) = R_i·[r−1_i]^(m+1) + (r_i−1)·[r−2_i]^(m+1)`
///
/// with `i` the first Cartesian axis (x, then y, then z) for which `r_i ≥ 1`,
/// and `R = AC`. The `[r]^m` are built level by level in the monomial degree.
///
/// @param contracted_seed The contracted seed ladder `[0]^m` — `order+1` rows
/// of `cdim` values, row stride padded to a multiple of 8.
/// @param order The total angular momentum `L = l_a + l_c`.
/// @param cdim The number of contracted pairs.
/// @param ac_x The x components of `AC = A − C`, one per contracted pair.
/// @param ac_y The y components of `AC`.
/// @param ac_z The z components of `AC`.
/// @return `[r]^0` for the degree-`L` Cartesian monomials `r` in canonical
/// order (x descending, then y descending) — `(L+1)(L+2)/2` rows of `cdim`
/// values, row stride padded to a multiple of 8.
auto compute_one_center_md(const std::vector<double> &contracted_seed,
                           const std::size_t          order,
                           const std::size_t          cdim,
                           const std::vector<double> &ac_x,
                           const std::vector<double> &ac_y,
                           const std::vector<double> &ac_z) -> std::vector<double>;

}  // namespace tabula

#endif /* TabulaMDRecursion_hpp */
