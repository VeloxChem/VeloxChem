//
//                                   VELOXCHEM
//              ----------------------------------------------------
//                          An Electronic Structure Code
//
//  SPDX-License-Identifier: BSD-3-Clause
//
//  Copyright 2018-2025 VeloxChem developers
//
//  Redistribution and use in source and binary forms, with or without modification,
//  are permitted provided that the following conditions are met:
//
//  1. Redistributions of source code must retain the above copyright notice, this
//     list of conditions and the following disclaimer.
//  2. Redistributions in binary form must reproduce the above copyright notice,
//     this list of conditions and the following disclaimer in the documentation
//     and/or other materials provided with the distribution.
//  3. Neither the name of the copyright holder nor the names of its contributors
//     may be used to endorse or promote products derived from this software without
//     specific prior written permission.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
//  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
//  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
//  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
//  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
//  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
//  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
//  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
//  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
//  OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


#ifndef ScreeningFunc_hpp
#define ScreeningFunc_hpp

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <ranges>
#include <vector>

#include "BasisFunction.hpp"
#include "MathConst.hpp"

/// @brief Defines supported integral bounds for screening of atom pairs:
/// screener::overlap           - the two-center overlap integral bound
/// screener::kinetic_energy    - the two-center kinetic energy integral bound
/// screener::nuclear_potential - the two-center nuclear potential integral bound
/// screener::electron_repulsion - the three-center electron repulsion integral bound
enum class screener
{
    overlap,
    kinetic_energy,
    nuclear_potential,
    electron_repulsion
};

namespace screenfunc {  // screenfunc namespace

/// @brief Counts leading elements of range satisfying predicate.
/// @param values The range of values to count leading elements in.
/// @param predicate The predicate to partition range with.
/// @return The number of leading elements satisfying predicate.
/// @note The range must be partitioned with respect to predicate, i.e. all
/// elements satisfying predicate must precede those which do not. Bisection is
/// used to locate the partition point, so the predicate is evaluated a
/// logarithmic number of times. A range which is not partitioned yields a wrong
/// count instead of an error, as detecting it requires the linear scan avoided
/// here.
template <typename T, typename P>
inline auto
number_of_leading(const std::vector<T> &values, P predicate) -> size_t
{
    return static_cast<size_t>(std::ranges::distance(values.begin(), std::ranges::partition_point(values, predicate)));
}

/// @brief Counts leading indices satisfying predicate.
/// @param nvalues The number of indices to count leading indices in.
/// @param predicate The predicate to partition indices with.
/// @return The number of leading indices satisfying predicate.
/// @note The indices must be partitioned with respect to predicate, as for the
/// range taking overload. This overload is used where the values are computed
/// from an index rather than stored, so that a linear pass over them is avoided.
template <typename P>
inline auto
number_of_leading(const size_t nvalues, P predicate) -> size_t
{
    const auto indices = std::views::iota(size_t{0}, nvalues);

    return static_cast<size_t>(std::ranges::distance(indices.begin(), std::ranges::partition_point(indices, predicate)));
}

/// @brief Computes the distance prefactor of the two-center integral bounds.
/// @param r The distance between the atoms carrying the basis functions.
/// @param lsum The sum of angular momenta of the basis functions.
/// @return The distance prefactor.
/// @note The prefactor is clamped to one below unit distance, where it would
/// lower the bound instead of raising it. Repeated multiplication is used, as
/// the sum of angular momenta is small.
inline auto
distance_factor(const double r, const int lsum) -> double
{
    auto fact = 1.0;

    if (r > 1.0)
    {
        for (int i = 0; i < lsum; i++)
        {
            fact *= r;
        }
    }

    return fact;
}

/// @brief Computes upper bound of two-center overlap integral between basis
/// functions on bra and ket sides.
/// @param bra The basis function on bra side.
/// @param ket The basis function on ket side.
/// @param r The distance between the atoms carrying the basis functions.
/// @return The upper bound of two-center overlap integral.
/// @note The bound rises with distance beyond unit distance before the Gaussian
/// decay takes over, so it is not monotone. It is nonetheless partitioned by any
/// threshold below its value at unit distance, which is the regime of screening
/// thresholds, and may thus be bisected with number_of_leading.
inline auto
two_center_overlap_bound(const CBasisFunction &bra, const CBasisFunction &ket, const double r) -> double
{
    const auto fexp = bra.smallest_exponent() + ket.smallest_exponent();

    const auto fmu = bra.smallest_exponent() * ket.smallest_exponent() / fexp;

    const auto fpi = mathconst::pi_value() / fexp;

    const auto fact = bra.sum_of_absolute_norms() * ket.sum_of_absolute_norms() *
                      distance_factor(r, bra.get_angular_momentum() + ket.get_angular_momentum());

    return fact * fpi * std::sqrt(fpi) * std::exp(-fmu * r * r);
}

/// @brief Computes upper bound of two-center kinetic energy integral between
/// basis functions on bra and ket sides.
/// @param bra The basis function on bra side.
/// @param ket The basis function on ket side.
/// @param r The distance between the atoms carrying the basis functions.
/// @return The upper bound of two-center kinetic energy integral.
/// @note The overlap part is computed with the smallest exponents, as for the
/// overlap bound, while the kinetic energy part is computed with the largest
/// exponents. The bound is partitioned by any threshold below its value at unit
/// distance, and may thus be bisected with number_of_leading.
inline auto
two_center_kinetic_energy_bound(const CBasisFunction &bra, const CBasisFunction &ket, const double r) -> double
{
    const auto fexp = bra.smallest_exponent() + ket.smallest_exponent();

    const auto fmu = bra.smallest_exponent() * ket.smallest_exponent() / fexp;

    const auto fpi = mathconst::pi_value() / fexp;

    const auto gexp = bra.largest_exponent() + ket.largest_exponent();

    const auto gmu = bra.largest_exponent() * ket.largest_exponent() / gexp;

    const auto fact = bra.sum_of_absolute_norms() * ket.sum_of_absolute_norms() *
                      distance_factor(r, bra.get_angular_momentum() + ket.get_angular_momentum());

    return fact * gmu * (3.0 + 2.0 * gmu * r * r) * fpi * std::sqrt(fpi) * std::exp(-fmu * r * r);
}

/// @brief Computes upper bound of two-center nuclear potential integral between
/// basis functions on bra and ket sides.
/// @param bra The basis function on bra side.
/// @param ket The basis function on ket side.
/// @param r The distance between the atoms carrying the basis functions.
/// @return The upper bound of two-center nuclear potential integral.
/// @note The Boys function is set to its maximum value of one. The bound is
/// partitioned by any threshold below its value at unit distance, and may thus
/// be bisected with number_of_leading.
inline auto
two_center_nuclear_potential_bound(const CBasisFunction &bra, const CBasisFunction &ket, const double r) -> double
{
    const auto fexp = bra.smallest_exponent() + ket.smallest_exponent();

    const auto fmu = bra.smallest_exponent() * ket.smallest_exponent() / fexp;

    const auto lsum = bra.get_angular_momentum() + ket.get_angular_momentum();

    const auto fact = bra.sum_of_absolute_norms() * ket.sum_of_absolute_norms() * distance_factor(r, lsum);

    return fact * static_cast<double>(lsum + 1) * (2.0 * mathconst::pi_value() / fexp) * std::exp(-fmu * r * r);
}

/// @brief Computes upper bound of three-center electron repulsion integral
/// between basis functions on a, b and c sides.
/// @param a The basis function on a side.
/// @param b The basis function on b side.
/// @param c The basis function on c side.
/// @param r The distance between the atoms carrying the basis functions on a and
/// b sides.
/// @return The upper bound of three-center electron repulsion integral.
/// @note The Boys function is set to its maximum value of one, so the distance
/// to the atom on c side does not enter and the c side contributes through the
/// exponents and contraction of its basis function only. The sum of angular
/// momenta runs over all three sides. The bound is partitioned by any threshold
/// below its value at unit distance, and may thus be bisected with
/// number_of_leading.
inline auto
three_center_electron_repulsion_bound(const CBasisFunction &a, const CBasisFunction &b, const CBasisFunction &c, const double r) -> double
{
    const auto fexp = a.smallest_exponent() + b.smallest_exponent();

    const auto fmu = a.smallest_exponent() * b.smallest_exponent() / fexp;

    const auto cexp = c.smallest_exponent();

    const auto qexp = fexp + cexp;

    const auto lsum = a.get_angular_momentum() + b.get_angular_momentum() + c.get_angular_momentum();

    constexpr auto fpi = mathconst::pi_value();

    const auto fact = a.sum_of_absolute_norms() * b.sum_of_absolute_norms() * c.sum_of_absolute_norms() * distance_factor(r, lsum) *
                      static_cast<double>(lsum + 1);

    return fact * (2.0 * fpi * fpi * std::sqrt(fpi)) / (fexp * cexp * std::sqrt(qexp)) * std::exp(-fmu * r * r);
}

/// @brief Computes upper bound of two-center overlap integral between primitive
/// basis functions on bra and ket sides.
/// @param bra The basis function on bra side.
/// @param iprim The index of primitive of basis function on bra side.
/// @param ket The basis function on ket side.
/// @param jprim The index of primitive of basis function on ket side.
/// @param r The distance between the atoms carrying the basis functions.
/// @return The upper bound of two-center overlap integral.
/// @note The bound is the contracted overlap bound with the exponent and the
/// absolute normalization factor of a single primitive in place of the smallest
/// exponent and the sum of absolute normalization factors.
inline auto
two_center_overlap_primitive_bound(const CBasisFunction &bra, const size_t iprim, const CBasisFunction &ket, const size_t jprim,
                                   const double r) -> double
{
    const auto aexp = bra.exponents()[iprim];

    const auto bexp = ket.exponents()[jprim];

    const auto fexp = aexp + bexp;

    const auto fmu = aexp * bexp / fexp;

    const auto fpi = mathconst::pi_value() / fexp;

    const auto fact = std::fabs(bra.normalization_factors()[iprim]) * std::fabs(ket.normalization_factors()[jprim]) *
                      distance_factor(r, bra.get_angular_momentum() + ket.get_angular_momentum());

    return fact * fpi * std::sqrt(fpi) * std::exp(-fmu * r * r);
}

}  // namespace screenfunc

#endif /* ScreeningFunc_hpp */
