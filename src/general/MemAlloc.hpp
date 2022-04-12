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

#ifndef MemAlloc_hpp
#define MemAlloc_hpp

#ifdef _MSC_VER
#include <malloc.h>
#else
#include <cstdlib>
#endif
#include <cstddef>
#include <numeric>
#include <string>
#include <tuple>
#include <type_traits>

#ifdef ENABLE_MKL
#include <mkl_service.h>  // needed for mkl_malloc and mkl_free
#endif

#include "Device.hpp"
#include "ErrorHandler.hpp"

namespace mem {
/**
 * Tag struct for host backend.
 */
struct Host
{
    inline const static std::string name{"Host"};
};

template <typename B>
inline constexpr bool is_on_host_v = std::is_same_v<B, Host>;

/** Compute pitch for aligned allocation.
 *
 * @tparam T type of elements.
 * @param[in] alignment requested alignment, in bytes.
 * @param[in] count number of elements.
 * @return pitch of the allocation: _i.e._ the `count` plus the padding, if any.
 */
template <typename T>
constexpr auto
get_pitch(size_t alignment, size_t count) -> size_t
{
    auto tile        = alignment / sizeof(T);
    auto [quot, rem] = std::ldiv(count, tile);

    if (rem) quot += 1;

    return (quot * tile);
}

namespace detail {
/** Aligned malloc to be used.
 *
 * - When linking against MKL, use `mkl_malloc`. Otherwise:
 *   * On Windows, use the intrinsic `_aligned_malloc`. See:
 * https://developercommunity.visualstudio.com/t/c17-stdaligned-alloc%E7%BC%BA%E5%A4%B1/468021
 *   * On Linux, use `::posix_memalign`
 */
constexpr inline auto __vlx_malloc =
#ifdef ENABLE_MKL
    mkl_malloc;
#else  // ENABLE_MKL
#ifdef _MSC_VER
    _aligned_alloc;
#else  // _MSC_VER
    [](size_t count_bytes, size_t alignment) -> void * {
    void* ptr = nullptr;

    auto ierr = ::posix_memalign(&ptr, alignment, count_bytes);

    errors::assertMsgCritical(ierr == 0, "posix_memalign failed");

    return ptr;
};
#endif
#endif

/** Aligned free to be used.
 *
 * - When linking against MKL, use `mkl_free`. Otherwise:
 *   * On Windows, use the intrinsic `_aligned_free`. See:
 * https://developercommunity.visualstudio.com/t/c17-stdaligned-alloc%E7%BC%BA%E5%A4%B1/468021
 *   * On Linux, use `::free`
 */
constexpr inline auto __vlx_free =
#ifdef ENABLE_MKL
    mkl_free;
#else  // ENABLE_MKL
#ifdef _MSC_VER
    _aligned_free;
#else  // _MSC_VER
    ::free;
#endif
#endif

/** Aligned allocation of given type on the host.
 *
 * @tparam T scalar type of the allocation.
 * @param[in] alignment desired alignment of allocation.
 * @param[in] count number of element in allocation, including padding.
 * @return pointer to the allocation
 *
 * @note Alignment must be a power of 2.
 */
template <typename T>
auto
host_allocate(size_t alignment, size_t count) -> T *
{
    if (count == 0)
    {
        return nullptr;
    }

    // check that alignment is a power of 2
    // http://www.graphics.stanford.edu/~seander/bithacks.html#DetermineIfPowerOf2
    if (alignment && !(alignment & (alignment - 1))) errors::msgCritical(std::string(__func__) + ": alignment must be a power of 2");

    // check that we are not trying to allocate too big a chunk
    if (count > std::numeric_limits<size_t>::max() / sizeof(T))
    {
        // equivalent of: throw std::bad_array_new_length();
        errors::msgCritical(std::string(__func__) + ": you cannot allocate a memory block with " + std::to_string(count) + " elements");

        // the useless return statement is to avoid warnings from the compiler
        return nullptr;
    }

    if (auto p = static_cast<T *>(__vlx_malloc(count * sizeof(T), alignment)))
    {
        return p;
    }
    else
    {
        // equivalent of: throw std::bad_alloc();
        errors::msgCritical(std::string(__func__) + ": aligned allocation failed");

        // the useless return statement is to avoid warnings from the compiler
        return nullptr;
    }
}

/** Deallocate chunk of memory of given type on host.
 *
 * @tparam T type of chunk.
 * @param[in] p the pointer to the memory chunk.
 */
template <typename T>
auto
host_deallocate(T *p) noexcept -> void
{
    if (p) __vlx_free(p);
}
}  // namespace detail

template <typename T, typename B = Host>
auto
malloc_1d(size_t alignment, size_t count) -> T *
{
    if constexpr (is_on_host_v<B>)
    {
        return detail::host_allocate<T>(alignment, count);
    }
    else
    {
        return detail::device_allocate<T>(count);
    }
}

template <typename T, typename B = Host>
auto
malloc(size_t count) -> T *
{
    return malloc_1d<T, B>(VLX_ALIGN, count);
}

template <typename T, typename B = Host>
auto
malloc_2d(size_t alignment, size_t height, size_t width) -> std::tuple<size_t, T *>
{
    if constexpr (is_on_host_v<B>)
    {
        auto pitch = get_pitch<T>(alignment, width);
        return {pitch, detail::host_allocate<T>(alignment, pitch * height)};
    }
    else
    {
        return detail::device_allocate<T>(height, width);
    }
}

template <typename T, typename B = Host>
auto
malloc(size_t height, size_t width) -> std::tuple<size_t, T *>
{
    return malloc_2d<T, B>(VLX_ALIGN, height, width);
}

/** Deallocate memory chunk pointed to by pointer.
 *
 * @tparam T type of memory chunk.
 * @tparam B memory allocation backend.
 * @param pointer the pointer to memory chunk.
 */
template <typename T, typename B>
auto
free(T *pointer) -> void
{
    if constexpr (is_on_host_v<B>)
    {
        detail::host_deallocate(pointer);
    }
    else
    {
        detail::device_deallocate(pointer);
    }
}

/** Deallocate memory chunk pointed to by host pointer.
 *
 * @tparam T type of memory chunk.
 * @param pointer the pointer to memory chunk.
 */
template <typename T>
auto
free(T *pointer) -> void
{
    detail::host_deallocate(pointer);
}
}  // namespace mem

#endif /* MemAlloc_hpp */
