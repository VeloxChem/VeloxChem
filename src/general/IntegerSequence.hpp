//
//                           VELOXCHEM 1.0-RC2
//         ----------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2021 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact
//
//  SPDX-License-Identifier: LGPL-3.0-or-later
//
//  This file is part of VeloxChem.
//
//  VeloxChem is free software: you can redistribute it and/or modify it under
//  the terms of the GNU Lesser General Public License as published by the Free
//  Software Foundation, either version 3 of the License, or (at your option)
//  any later version.
//
//  VeloxChem is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
//  License for more details.
//
//  You should have received a copy of the GNU Lesser General Public License
//  along with VeloxChem. If not, see <https://www.gnu.org/licenses/>.

#ifndef IntegerSequence_hpp
#define IntegerSequence_hpp

#include <cstddef>
#include <type_traits>

namespace buffer {
// forward declarations
template <typename First, typename Size, typename Increment>
class IntegerSequence;
template <typename First, typename Size, typename Increment>
[[nodiscard]] auto seqN(First first, Size size, Increment inc) -> IntegerSequence<First, Size, Increment>;

template <typename First, typename Size, typename Increment>
class IntegerSequence final
{
    static_assert(std::is_integral_v<First>, "First element in sequence must be of integer type!");
    static_assert(std::is_integral_v<Size>, "Size of sequence must be of integer type!");
    static_assert(std::is_integral_v<Increment>, "Increment of sequence must be of integer type!");

   private:
    First     _first;
    Size      _size;
    Increment _increment;

   public:
    IntegerSequence(First f, Size s, Increment i) : _first(f), _size(s), _increment(i)
    {
    }

    std::ptrdiff_t
    size() const
    {
        return _size;
    }

    std::ptrdiff_t
    operator[](std::ptrdiff_t i) const
    {
        return _first + i * _increment;
    }

    auto
    reverse() const -> decltype(seqN(_first + (_size - 1) * _increment, _size, -_increment))
    {
        return seqN(_first + (_size - 1) * _increment, _size, -_increment);
    }
};

/** A sequence of integers of given length, starting point, and increment.
 *
 * @tparam First type of first element in sequence.
 * @tparam Size type of sequence size.
 * @tparam Increment type of increment.
 * @param[in] f first element in sequence.
 * @param[in] s size of sequence.
 * @param[in] inc increment between adjacent elements in sequence.
 *
 * For example: `seqN(2, 3, 3)` is a sequence of 3 integers, starting from 2 and
 * with increment of 3: `{2, 5, 8}`.
 */
template <typename First, typename Size, typename Increment>
[[nodiscard]] auto
seqN(First f, Size s, Increment inc) -> IntegerSequence<First, Size, Increment>
{
    return IntegerSequence(f, s, inc);
}

/** A sequence of integers of given length and starting point, with unit increment.
 *
 * @tparam First type of first element in sequence.
 * @tparam Size type of sequence size.
 * @param[in] f first element in sequence.
 * @param[in] s size of sequence.
 *
 * For example: `seqN(2, 5)` is a sequence of 5 integers, starting from 2 and
 * with increment of 1: `{2, 3, 4, 5, 6}`.
 */
template <typename First, typename Size>
[[nodiscard]] auto
seqN(First f, Size s) -> IntegerSequence<First, Size, decltype(1)>
{
    return IntegerSequence(f, s, 1);
}

/** A sequence of integers of given starting point, end point, and increment.
 *
 * @tparam First type of first element in sequence.
 * @tparam Last type of last element in sequence.
 * @tparam Increment type of increment.
 * @param[in] f first element in sequence.
 * @param[in] l last element in sequence.
 * @param[in] inc increment between adjacent elements in sequence.
 *
 * For example: `seq(2, 8, 2)` is the sequence starting at 2, ending at 8, with
 * an increment of 2: `{2, 4, 6, 8}`.
 */
template <typename First, typename Last, typename Increment>
[[nodiscard]] auto
seq(First f, Last l, Increment inc) -> IntegerSequence<First, decltype((l - f + inc) / inc), Increment>
{
    return seqN(f, (l - f + inc) / inc, inc);
}

/** A sequence of integers of given starting point and end point, with increment of 1.
 *
 * @tparam First type of first element in sequence.
 * @tparam Last type of last element in sequence.
 * @param[in] f first element in sequence.
 * @param[in] l last element in sequence.
 *
 * For example: `seq(2, 5)` is the sequence starting at 2, ending at 5, with
 * an increment of 1: `{2, 3, 4, 5}`.
 */
template <typename First, typename Last>
[[nodiscard]] auto
seq(First f, Last l) -> IntegerSequence<First, decltype((l - f + 1)), decltype(1)>
{
    return seq(f, l, 1);
}
}  // namespace buffer

#endif
