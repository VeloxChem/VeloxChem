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

#ifndef MetaUtils_hpp
#define MetaUtils_hpp

#include <mpi.h>

#include <climits>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numeric>
#include <string>
#include <type_traits>

namespace metautils {
/* \@{ From type to MPI_Datatype enum variant and string representation */
template <typename T>
struct type_to_mpi_datatype;

template <>
struct type_to_mpi_datatype<bool>
{
    inline const static MPI_Datatype value{MPI_C_BOOL};
};

template <>
struct type_to_mpi_datatype<char>
{
    inline const static MPI_Datatype value{MPI_CHAR};
};

template <>
struct type_to_mpi_datatype<int32_t>
{
    inline const static MPI_Datatype value{MPI_INT32_T};
};

template <>
struct type_to_mpi_datatype<int64_t>
{
    inline const static MPI_Datatype value{MPI_INT64_T};
};

template <>
struct type_to_mpi_datatype<std::size_t>
{
#if SIZE_MAX == UCHAR_MAX
    inline const static MPI_Datatype value{MPI_UNSIGNED_CHAR};
#elif SIZE_MAX == USHRT_MAX
    inline const static MPI_Datatype value{MPI_UNSIGNED_SHORT};
#elif SIZE_MAX == UINT_MAX
    inline const static MPI_Datatype value{MPI_UNSIGNED};
#elif SIZE_MAX == ULONG_MAX
    inline const static MPI_Datatype value{MPI_UNSIGNED_LONG};
#elif SIZE_MAX == ULLONG_MAX
    inline const static MPI_Datatype value{MPI_UNSIGNED_LONG_LONG};
#else
#error "Your compiler does not know what SIZE_MAX is?"
#endif
};

template <>
struct type_to_mpi_datatype<double>
{
    inline const static MPI_Datatype value{MPI_DOUBLE};
};
/* \@} */

/* \@{ From type to string representation */
template <typename T>
struct type_to_string;

template <>
struct type_to_string<bool>
{
    inline const static std::string name{"bool"};
};

template <>
struct type_to_string<char>
{
    inline const static std::string name{"char"};
};

template <>
struct type_to_string<int32_t>
{
    inline const static std::string name{"int32_t"};
};

template <>
struct type_to_string<int64_t>
{
    inline const static std::string name{"int64_t"};
};

template <>
struct type_to_string<std::size_t>
{
    inline const static std::string name{"std::size_t"};
};

template <>
struct type_to_string<float>
{
    inline const static std::string name{"float"};
};

template <>
struct type_to_string<double>
{
    inline const static std::string name{"double"};
};
/* \@} */

/** @{ Compile-time condition for types. */
/** Compile-time condition.
 *
 * @tparam B condition.
 * @tparam T type.
 */
template <bool B, typename T>
struct condition
{
    static constexpr bool value = B;
    using type                  = T;
};
/* \@} */

namespace detail {
/** @{ Compile-time switch-statement for types.
 *
 * See: https://stackoverflow.com/a/32785263
 */
template <typename Head, typename... Tail>
struct select : std::conditional_t<Head::value, Head, select<Tail...>>
{
};

/** Recursion bottom. */
template <typename T>
struct select<T>
{
    using type = T;
};

/** Recursion bottom. */
template <bool B, typename T>
struct select<condition<B, T>>
{
    static_assert(B, "None of the conditions were met!");
    using type = T;
};
/* \@} */
}  // namespace detail

/** @{ Compile-time switch-statement for types. */
/** Compile-time switch-statement for types.
 *
 * Usage:
 *
 * \code{.cpp}
 * using T = select_t<condition<..., T1>, condition<..., T2>, ..., Tf>;
 * \endcode
 *
 * `Tf` is the fallback type, selected if none of the conditions are met.
 *
 * See: https://stackoverflow.com/a/32785263
 */
template <typename Head, typename... Tail>
using select_t = typename detail::select<Head, Tail...>::type;
/**@}*/
}  // namespace metautils

#endif /* MetaUtils_hpp */
