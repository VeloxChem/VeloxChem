//
//  Tabula — custom-recursion molecular-integral machinery.
//  Tabula-owned basis-function-pair block.
//

#ifndef TabulaGtoPairBlock_hpp
#define TabulaGtoPairBlock_hpp

#include <cstddef>
#include <functional>
#include <memory>
#include <utility>
#include <vector>

#include "GtoBlock.hpp"
#include "GtoBlockScreeningData.hpp"
#include "Point.hpp"

namespace tabula {  // tabula namespace

/// @brief The contracted-GTO-pair screening estimator: given the screening
/// data of a bra and a ket contracted GTO and the distance `|R|` between
/// their centers, returns the screening estimate. A contracted-GTO pair is
/// kept when its estimate is at or above the screening threshold.
using ScreeningEstimator =
    std::function<double(const CGtoBlockScreeningData &, const CGtoBlockScreeningData &, const double)>;

/// @brief A Tabula-owned basis-function-pair block.
///
/// Built from two `CGtoBlock`s, this carries exactly the structure-of-arrays
/// the late-contraction overlap recursion consumes — and nothing else. It
/// differs from VeloxChem's `CGtoPairBlock` in three ways that cut the
/// construction cost:
///  - the primitive-pair arrays are allocated uninitialized (the fill writes
///    every element, so the value-initialization is pure waste);
///  - the primitive normalization and overlap factors are folded into a
///    single `weights` array (the seed step only ever uses their product);
///  - the accessors hand back pointers / const references, never by-value
///    copies of the per-primitive arrays.
///
/// The primitive-pair arrays are primitive-pair-major: the value for
/// primitive pair `pp` of contracted pair `ij` sits at `pp·cdim + ij`.
class GtoPairBlock
{
   public:
    /// @brief Creates an empty pair block.
    GtoPairBlock() = default;

    /// @brief Creates an unscreened pair block — every contracted-GTO pair of
    /// the two blocks is kept.
    /// @param bra_gto_block The basis functions block on the bra side.
    /// @param ket_gto_block The basis functions block on the ket side.
    GtoPairBlock(const CGtoBlock &bra_gto_block, const CGtoBlock &ket_gto_block);

    /// @brief Creates a screened pair block — keeps the contracted-GTO pairs
    /// whose screening estimate is at or above the threshold.
    /// @param bra_gto_block The basis functions block on the bra side.
    /// @param ket_gto_block The basis functions block on the ket side.
    /// @param estimator The contracted-GTO-pair screening estimator.
    /// @param threshold The screening threshold; a threshold at or below `0`
    /// keeps every pair and skips the estimator entirely.
    GtoPairBlock(const CGtoBlock          &bra_gto_block,
                 const CGtoBlock          &ket_gto_block,
                 const ScreeningEstimator &estimator,
                 const double              threshold);

    /// @brief Creates an unscreened pair block over a sub-range of the bra
    /// block's contracted GTOs — for splitting one block pair into finer
    /// units of parallel work. The pair block holds the contracted pairs
    /// `(i, j)` with `i` in `[bra_begin, bra_end)` and `j` over the whole ket
    /// block; the orbital indices still address the global AO layout.
    /// @param bra_gto_block The basis functions block on the bra side.
    /// @param ket_gto_block The basis functions block on the ket side.
    /// @param bra_begin The first bra contracted-GTO index of the range.
    /// @param bra_end The one-past-last bra contracted-GTO index of the range.
    GtoPairBlock(const CGtoBlock &bra_gto_block, const CGtoBlock &ket_gto_block, const int bra_begin, const int bra_end);

    /// @brief Gets the angular momentums of the basis-function pair.
    /// @return The `(l_a, l_c)` angular-momentum pair.
    auto angular_momentums() const -> std::pair<int, int>;

    /// @brief Gets the number of contracted-GTO pairs.
    /// @return The contracted-pair count.
    auto number_of_contracted_pairs() const -> std::size_t;

    /// @brief Gets the number of primitive pairs per contracted pair.
    /// @return The primitive-pair count.
    auto number_of_primitive_pairs() const -> int;

    /// @brief Gets the Cartesian coordinates of the bra centers.
    /// @return The per-contracted-pair bra coordinates.
    auto bra_coordinates() const -> const std::vector<TPoint<double>> &;

    /// @brief Gets the Cartesian coordinates of the ket centers.
    /// @return The per-contracted-pair ket coordinates.
    auto ket_coordinates() const -> const std::vector<TPoint<double>> &;

    /// @brief Gets the AO indices on the bra side; element `0` is the AO
    /// component stride, element `ij+1` the offset of contracted pair `ij`.
    /// @return The bra orbital indices.
    auto bra_orbital_indices() const -> const std::vector<std::size_t> &;

    /// @brief Gets the AO indices on the ket side; element `0` is the AO
    /// component stride, element `ij+1` the offset of contracted pair `ij`.
    /// @return The ket orbital indices.
    auto ket_orbital_indices() const -> const std::vector<std::size_t> &;

    /// @brief Gets the primitive-pair bra exponents (primitive-pair-major).
    /// @return The bra exponents.
    auto bra_exponents() const -> const double *;

    /// @brief Gets the primitive-pair ket exponents (primitive-pair-major).
    /// @return The ket exponents.
    auto ket_exponents() const -> const double *;

    /// @brief Gets the primitive-pair weights — the normalization factor and
    /// the overlap factor folded into one (primitive-pair-major).
    /// @return The weights.
    auto weights() const -> const double *;

   private:
    /// @brief The Cartesian coordinates of the bra centers.
    std::vector<TPoint<double>> _bra_coordinates;

    /// @brief The Cartesian coordinates of the ket centers.
    std::vector<TPoint<double>> _ket_coordinates;

    /// @brief The AO indices on the bra side.
    std::vector<std::size_t> _bra_orb_indices;

    /// @brief The AO indices on the ket side.
    std::vector<std::size_t> _ket_orb_indices;

    /// @brief The primitive-pair bra exponents.
    std::unique_ptr<double[]> _bra_exponents;

    /// @brief The primitive-pair ket exponents.
    std::unique_ptr<double[]> _ket_exponents;

    /// @brief The primitive-pair weights (normalization × overlap factor).
    std::unique_ptr<double[]> _weights;

    /// @brief The angular momentums of the basis-function pair.
    std::pair<int, int> _angular_momentums{-1, -1};

    /// @brief The number of primitive pairs per contracted pair.
    int _nppairs{0};

    /// @brief Builds the pair block from two basis-function blocks, over the
    /// bra contracted-GTO range `[bra_begin, bra_end)`.
    /// @param bra_gto_block The basis functions block on the bra side.
    /// @param ket_gto_block The basis functions block on the ket side.
    /// @param estimator The screening estimator, or `nullptr` to keep every
    /// contracted-GTO pair.
    /// @param threshold The screening threshold.
    /// @param bra_begin The first bra contracted-GTO index of the range.
    /// @param bra_end The one-past-last bra contracted-GTO index of the range.
    auto _build(const CGtoBlock          &bra_gto_block,
                const CGtoBlock          &ket_gto_block,
                const ScreeningEstimator *estimator,
                const double              threshold,
                const int                 bra_begin,
                const int                 bra_end) -> void;
};

}  // namespace tabula

#endif /* TabulaGtoPairBlock_hpp */
