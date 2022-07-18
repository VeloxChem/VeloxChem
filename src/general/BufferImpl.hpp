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

#ifndef BufferImpl_hpp
#define BufferImpl_hpp

#include <array>
#include <cstddef>
#include <type_traits>
#include <vector>

#include "Device.hpp"
#include "ErrorHandler.hpp"
#include "MathFunc.hpp"
#include "MemAlloc.hpp"
#include "MetaUtils.hpp"
#include "mdspan.hpp"

/** Tag for size known at run-time.
 *
 * This is of `size_t` type, its value is not important.
 */
inline constexpr auto Dynamic = std::experimental::dynamic_extent;

namespace buffer {
namespace stdex = std::experimental;

enum class Kind
{
    X,
    N,
    XY,
    MY,
    XN,
    MN
};

template <decltype(Dynamic) NRows, decltype(Dynamic) NCols>
constexpr auto
getKind() -> Kind
{
    if constexpr (NRows == 1 && NCols == Dynamic)
    {
        return Kind::X;
    }
    else if constexpr (NRows == 1 && NCols != Dynamic)
    {
        return Kind::N;
    }
    else if constexpr (NRows == Dynamic && NCols == Dynamic)
    {
        return Kind::XY;
    }
    else if constexpr (NRows != Dynamic && NCols == Dynamic)
    {
        return Kind::MY;
    }
    else if constexpr (NRows == Dynamic && NCols != Dynamic)
    {
        return Kind::XN;
    }
    else
    {
        return Kind::MN;
    }
    // default case
    return Kind::XY;
}

/** Class to manage memory buffer allocation, manipulation, and deallocation.
 *
 * @tparam T scalar type of buffer elements. Must be an arithmetic type.
 * @tparam B backend of buffer allocation.
 * @tparam NRows number of rows at compile-time.
 * @tparam NCols number of columns at compile-time.
 * @author Roberto Di Remigio
 *
 * This class is implicitly convertible to `mdspan`.
 * @note The buffer is aligned to the byte-boundary appropriate for the backend.
 */
template <typename T, typename B = mem::Host, auto NRows = Dynamic, auto NCols = Dynamic>
class CBuffer
{
    static_assert(std::is_arithmetic_v<T>, "CBuffer can only be instantiated with arithmetic types.");

    // Rank-1 objects are all stored and handled as row vectors!
    static_assert(NCols != 1, "CBuffer does not allow for storing column vectors.");

    static constexpr bool Layout1D                  = (NRows == 1);
    static constexpr bool CompileTimeRows           = (NRows != Dynamic);
    static constexpr bool CompileTimeColumns        = (NCols != Dynamic);
    static constexpr bool CompileTimeRowsAndColumns = (CompileTimeRows && CompileTimeColumns);
    static constexpr auto NResizableExtents         = CompileTimeRowsAndColumns ? 0 : ((CompileTimeRows || CompileTimeColumns) ? 1 : 2);

   public:
    using size_type       = std::decay_t<decltype(Dynamic)>;
    using value_type      = T;
    using backend_type    = B;
    using reference       = value_type &;
    using const_reference = const value_type &;
    using pointer         = value_type *;
    using const_pointer   = const value_type *;

    static constexpr size_type Alignment = [] {
        if constexpr (mem::is_on_host_v<B>)
        {
            return VLX_ALIGN;
        }
        else
        {
            // if we are on the device, the device allocation functions will
            // decide alignment and padding.  Thus, we set the alignment to the
            // size (in bytes) of the type, such that the computation of the
            // pitch does *not* add spurious padding.
            return sizeof(T);
        }
    }();

    static constexpr Kind kind = getKind<NRows, NCols>();

   private:
    /** The number of rows in memory buffer. */
    size_type _nRows{NRows};

    /** The number of columns in memory buffer. */
    size_type _nColumns{NCols};

    /** The number of columns plus padding in memory buffer. This is also called the pitch. */
    size_type _nPaddedColumns{mem::get_pitch<value_type>(Alignment, NCols)};

    /** Total number of elements in the allocation: _nRows * _nPaddedColumns */
    size_type _nElements{0};

    /** The contiguous and aligned memory buffer. */
    pointer _data{nullptr};

    /** Copy aligned data from `src` to `_data`.
     *
     * @param src the pointer to data source.
     *
     * @note Only available for buffers on the host. The padding elements are copied as well.
     */
    template <typename B_ = backend_type, typename = std::enable_if_t<mem::is_on_host_v<B_>>>
    auto
    _copy_aligned(const_pointer src) -> void
    {
        auto pdata = _data;

#pragma omp simd aligned(pdata, src : Alignment)
        for (size_type i = 0; i < _nElements; ++i)
            pdata[i] = src[i];
    }

    /** Copy 1-dimensional unaligned data to buffer.
     *
     * @param[in] src source memory.
     *
     * The `src` pointer is the start address of unaligned 1D data of `_nColumns` elements.
     * The `_data` pointer in `CBuffer` is aligned, so `_nElements == _nPaddedColumns`.
     *
     * For example, a 1D buffer holding 9 `double`s aligned to 64-byte boundary,
     * will be padded with 7 additional `double`s:
     *
     * [ o o o o o o o o | o x x x x x x x ]
     *
     * However, we only need to copy the first `_nColumns` elements!
     */
    template <auto L_ = Layout1D, std::enable_if_t<L_, bool> = true>
    auto
    _copy_unaligned(const_pointer src) -> void
    {
        if constexpr (mem::is_on_host_v<backend_type>)
        {
            auto pdata = _data;

#pragma omp simd
            for (size_type i = 0; i < _nColumns; ++i)
                pdata[i] = src[i];
        }
        else
        {
            DEVICE_CHECK(deviceMemcpy(_data, const_cast<value_type *>(src), _nColumns * sizeof(value_type), H2D));
        }
    }

    /** Fill 1-dimensional buffer with given value.
     *
     * @param[in] fill_value the filler for the buffer.
     */
    template <auto L_ = Layout1D, std::enable_if_t<L_, bool> = true>
    auto
    _fill(value_type fill_value) -> void
    {
        if constexpr (mem::is_on_host_v<backend_type>)
        {
            auto pdata = _data;
#pragma omp simd aligned(pdata : Alignment)
            for (size_type i = 0; i < _nColumns; ++i)
                pdata[i] = fill_value;
        }
        else
        {
#ifdef VLX_USE_DEVICE
            const auto block = dim3(256);
            const auto grid  = dim3((_nColumns + block.x - 1) / block.x);
            deviceLaunch((device::full1D<value_type>), grid, block, 0, 0, _data, fill_value, _nColumns);
            DEVICE_CHECK(deviceStreamSynchronize(0));
#endif
        }
    }

    /** Copy 2-dimensional unaligned data to buffer.
     *
     * @param[in] src source memory.
     */
    template <auto L_ = Layout1D, std::enable_if_t<!L_, bool> = true>
    auto
    _copy_unaligned(const_pointer src) -> void
    {
        if constexpr (mem::is_on_host_v<backend_type>)
        {
            for (size_type i = 0; i < _nRows; i++)
            {
                auto doff = i * _nPaddedColumns;

                auto soff = i * _nColumns;
                for (size_type j = 0; j < _nColumns; j++)
                {
                    _data[doff + j] = src[soff + j];
                }
            }
        }
        else
        {
            DEVICE_CHECK(deviceMemcpy2D(_data,
                                        _nPaddedColumns * sizeof(value_type), /* dpitch (bytes) */
                                        src,
                                        _nColumns * sizeof(value_type), /* spitch (bytes) */
                                        _nColumns * sizeof(value_type), /* width (bytes) */
                                        _nRows,
                                        H2D));
        }
    }

    /** Fill 2-dimensional buffer with given value.
     *
     * @param[in] fill_value the filler for the buffer.
     */
    template <auto L_ = Layout1D, std::enable_if_t<!L_, bool> = true>
    auto
    _fill(value_type fill_value) -> void
    {
        if constexpr (mem::is_on_host_v<backend_type>)
        {
            for (size_type i = 0; i < _nRows; ++i)
            {
                auto pdata = _data;
#pragma omp simd aligned(pdata : Alignment)
                for (size_type j = 0; j < _nColumns; ++j)
                {
                    pdata[i * _nPaddedColumns + j] = fill_value;
                }
            }
        }
        else
        {
#ifdef VLX_USE_DEVICE
            const auto block = dim3(16, 16);
            const auto grid  = dim3((_nRows + block.x - 1) / block.x, (_nColumns + block.y - 1) / block.y);
            deviceLaunch((device::full2D<value_type>), grid, block, 0, 0, _data, fill_value, _nRows, _nColumns, _nPaddedColumns);
            DEVICE_CHECK(deviceStreamSynchronize(0));
#endif
        }
    }

    /** Fill buffer with random values in interval.
     *
     * @param[in] lower lower bound of interval.
     * @param[in] upper upper bound of interval.
     *
     * @note The random numbers are generated on the CPU, then copied into the buffer.
     * This method uses the C++ default random engine with random seeding.
     * If you need more control, generate the random sequence (as std::vector or
     * raw array) and then use any of the available constructors.
     */
    auto
    _random_fill(T lower, T upper) -> void
    {
        // generate
        auto tmp = new T[_nElements];
        mathfunc::fill_random(tmp, lower, upper, _nElements);

        // copy into buffer
        _copy_unaligned(tmp);

        // delete temporary
        delete[] tmp;
    }

    /** Allocation function.
     *
     * @tparam Extents variadic pack of extents dimensions.
     *
     * @param[in] extents dimensions of the extents.
     */
    template <typename... Extents>
    auto
    _allocate(Extents... extents) -> void
    {
        static_assert(std::conjunction_v<std::is_integral<Extents>...>, "Extents must be of integral type.");

        if (_data)
        {
            errors::msgCritical(std::string(__func__) + ": buffer already allocated!");
        }

        // unpack parameter
        std::array<size_type, NResizableExtents> tmp{static_cast<size_type>(extents)...};

        if constexpr (kind == Kind::X)
        {
            _nColumns       = tmp[0];
            _nPaddedColumns = mem::get_pitch<value_type>(Alignment, _nColumns);
            _nElements      = _nRows * _nPaddedColumns;
            _data           = mem::malloc<value_type, backend_type>(_nElements);
        }
        else if constexpr (kind == Kind::XY)
        {
            _nRows                           = tmp[0];
            _nColumns                        = tmp[1];
            std::tie(_nPaddedColumns, _data) = mem::malloc<value_type, backend_type>(_nRows, _nColumns);
            _nElements                       = _nRows * _nPaddedColumns;
        }
        else if constexpr (kind == Kind::XN)
        {
            _nRows                           = tmp[0];
            std::tie(_nPaddedColumns, _data) = mem::malloc<value_type, backend_type>(_nRows, _nColumns);
            _nElements                       = _nRows * _nPaddedColumns;
        }
        else if constexpr (kind == Kind::MY)
        {
            _nColumns                        = tmp[0];
            std::tie(_nPaddedColumns, _data) = mem::malloc<value_type, backend_type>(_nRows, _nColumns);
            _nElements                       = _nRows * _nPaddedColumns;
        }
        else
        {
            _nElements = _nRows * _nPaddedColumns;
            _data      = mem::malloc<value_type, backend_type>(_nElements);
        }
    }

   public:
    /** Default CTOR.
     *
     * @warning This CTOR **only** allocates memory when *all* the extents of the
     * buffer are known at compile-time.  In all other cases, you will have to
     * call `resize` to allocate or, alternatively one of
     * `setConstant`/`setZero`/`setRandom` to allocate *and* initialize.
     */
    CBuffer()
    {
        if constexpr (kind == Kind::N || kind == Kind::MN)
        {
            _allocate();
        }
    }

    /** CTOR from dimensions.
     *
     * @tparam Extents variadic pack of extents dimensions.
     *
     * @param[in] extents dimensions of the extents.
     *
     * @note Only defined when some of the buffer extents are not known at
     * compile-time.
     * @note These CTOR **allocates** memory, but leaves it uninitialized.
     * You can call one of `setConstant`/`setZero`/`setRandom` to do so.
     */
    template <typename... Extents,
              auto CT_                             = CompileTimeRowsAndColumns,
              std::enable_if_t<!CT_, bool>         = true,
              auto AllIntegral_                    = std::conjunction_v<std::is_integral<Extents>...>,
              std::enable_if_t<AllIntegral_, bool> = true>
    explicit CBuffer(Extents... extents)
    {
        _allocate(extents...);
    }

    /** @{ CTORs from unaligned data buffers.
     *
     * @note These CTORs **allocate** memory *and* initialize it with the data
     * in the unaligned buffer, _i.e._ a `std::vector` or a raw pointer.
     */
    /** CTOR for buffer with run-time number of elements.
     *
     * @tparam Extents variadic pack of extents dimensions.
     *
     * @param[in] v raw array of data to copy into 1D buffer.
     * @param[in] extents dimensions of the extents.
     *
     * @warning the input array **must** have at least the same dimensions in each extent!
     */
    template <typename... Extents>
    explicit CBuffer(const_pointer v, Extents... extents) : CBuffer<value_type, backend_type, NRows, NCols>{static_cast<size_type>(extents)...}
    {
        _copy_unaligned(v);
    }

    /** CTOR for buffer with run-time number of elements.
     *
     * @tparam Extents variadic pack of extents dimensions.
     *
     * @param[in] v vector of data to copy into 1D buffer.
     * @param[in] extents dimensions of the extents.
     *
     * @warning the input array **must** have at least the same dimensions in each extent!
     */
    template <typename... Extents>
    explicit CBuffer(const std::vector<value_type> &v, Extents... extents)
        : CBuffer<value_type, backend_type, NRows, NCols>{v.data(), static_cast<size_type>(extents)...}
    {
        errors::assertMsgCritical(v.size() >= _nColumns,
                                  std::string(__func__) + ": input vector must have at least " + std::to_string(_nColumns) + " elements");
    }

    /** CTOR for 1D buffer with run-time number of elements.
     *
     * @param[in] v vector of data to copy into 1D buffer
     *
     * @note This CTOR takes the dimension from `std::vector`.
     */
    template <auto K_ = kind, std::enable_if_t<(K_ == Kind::X), bool> = true>
    explicit CBuffer(const std::vector<value_type> &v) : CBuffer<value_type, backend_type, NRows, NCols>{v.data(), v.size()}
    {
    }
    /** @} */

    /** @{ Copy CTOR. */
    /** Create a buffer object by copying other buffer object.
     *
     * @param src the source buffer object.
     * @note Takes care of host-to-host and device-to-device copy-constructions.
     */
    CBuffer(const CBuffer &src)
        : _nRows{src._nRows}
        , _nColumns{src._nColumns}
        , _nPaddedColumns{src._nPaddedColumns}
        , _nElements{src._nElements}
        , _data{mem::malloc<value_type, backend_type>(_nElements)}
    {
        if constexpr (mem::is_on_host_v<B>)
        {
            _copy_aligned(src._data);
        }
        else
        {
            if constexpr (Layout1D)
            {
                DEVICE_CHECK(deviceMemcpy(_data, src._data, _nElements * sizeof(T), D2D));
            }
            else
            {
                DEVICE_CHECK(deviceMemcpy2D(_data,
                                            _nPaddedColumns * sizeof(value_type), /* dpitch (bytes) */
                                            src._data,
                                            src._nPaddedColumns * sizeof(value_type), /* spitch (bytes) */
                                            _nColumns * sizeof(value_type),           /* width (bytes) */
                                            _nRows,
                                            D2D));
            }
        }
    }
    /**@}*/

    /** @{ Copy-assignment operator. */
    /** Create a buffer object by copy-assignment of other buffer object.
     *
     * @param src the source buffer object.
     * @note This overload takes care of host-to-host and device-to-device
     * copy-assignments.
     */
    auto
    operator=(const CBuffer &src) -> CBuffer &
    {
        if (&src != this)
        {
            _nRows          = src._nRows;
            _nColumns       = src._nColumns;
            _nPaddedColumns = src._nPaddedColumns;
            _nElements      = src._nElements;

            _data = mem::malloc<value_type, B>(_nElements);

            if constexpr (mem::is_on_host_v<B>)
            {
                _copy_aligned(src._data);
            }
            else
            {
                if constexpr (Layout1D)
                {
                    DEVICE_CHECK(deviceMemcpy(_data, src._data, _nElements * sizeof(T), D2D));
                }
                else
                {
                    DEVICE_CHECK(deviceMemcpy2D(_data,
                                                _nPaddedColumns * sizeof(value_type), /* dpitch (bytes) */
                                                src._data,
                                                src._nPaddedColumns * sizeof(value_type), /* spitch (bytes) */
                                                _nColumns * sizeof(value_type),           /* width (bytes) */
                                                _nRows,
                                                D2D));
                }
            }
        }

        return *this;
    }
    /**@}*/

    /** @{ Move CTORs. */
    /** Create a buffer object by moving other 1D buffer object.
     *
     * @param src the source buffer object.
     * @note Move-CTOR only applicable to host-side buffers.
     */
    template <typename B_ = B, typename = std::enable_if_t<!mem::is_on_device_v<B_>>>
    CBuffer(CBuffer &&src) noexcept
        : _nRows{src._nRows}, _nColumns{src._nColumns}, _nPaddedColumns{src._nPaddedColumns}, _nElements{src._nElements}, _data{src._data}
    {
        src._nRows          = NRows;
        src._nColumns       = NCols;
        src._nPaddedColumns = mem::get_pitch<value_type>(Alignment, NCols);
        src._nElements      = 0;
        src._data           = nullptr;
    }
    /**@}*/

    /** @{ Move-assignment operator. */
    /** Create a buffer object by move-assignment of other buffer object.
     *
     * @param src the source buffer object.
     * @note Move-assignment only applicable to host-side buffers.
     */
    template <typename B_ = B, typename = std::enable_if_t<!mem::is_on_device_v<B_>>>
    auto
    operator=(CBuffer &&src) noexcept -> CBuffer &
    {
        if (&src != this)
        {
            _nRows          = src._nRows;
            _nColumns       = src._nColumns;
            _nPaddedColumns = src._nPaddedColumns;
            _nElements      = src._nElements;

            mem::free<value_type, backend_type>(_data);

            _data = src._data;

            src._data = nullptr;
        }

        return *this;
    }
    /**@}*/

    /** Destroys a memory buffer object. */
    ~CBuffer() noexcept
    {
        _nElements = 0;
        mem::free<value_type, backend_type>(_data);
    }

    /** Resize existing buffer.
     *
     * @tparam Extents variadic pack of extents dimensions.
     *
     * @param[in] extents new dimension of the extents.
     *
     * @note This is a **no-op** when all dimensions are known at compile-time.
     * @warning Resizing is **always** non-conservative, _i.e._ the current
     * contents of the buffer  will be discarded.
     */
    template <typename... Extents>
    auto
    resize(Extents... extents) -> void
    {
        static_assert(std::conjunction_v<std::is_integral<Extents>...>, "New extents must be of integral type.");

        if constexpr ((sizeof...(extents) == NResizableExtents) && (!CompileTimeRowsAndColumns))
        {
            clear();

            // resize and reallocate
            _allocate(extents...);
        }
    }

    friend class CBuffer<T, mem::Host, NRows, NCols>;
    friend class CBuffer<T, mem::Device, NRows, NCols>;

    /** @{ Host-to-device and device-to-host converting CTORs */
    /** Converts a 1D buffer object between host and device backends, by creating a copy.
     *
     * @tparam BSource backend of source buffer object.
     * @param src the source buffer object.
     * @note This function takes care of host-to-device and device-to-host
     * copy-constructions.
     */
    template <auto L_ = Layout1D,
              typename BSource,
              std::enable_if_t<L_, bool> = true,
              typename                   = std::enable_if_t<!std::is_same_v<backend_type, BSource>>>
    explicit CBuffer(const CBuffer<T, BSource, NRows, NCols> &src)
        : _nRows{src._nRows}
        , _nColumns{src._nColumns}
        , _nPaddedColumns{src._nPaddedColumns}
        , _nElements{src._nElements}
        , _data{mem::malloc<value_type, backend_type>(_nElements)}
    {
        if constexpr (mem::is_on_host_v<BSource>)
        {
            DEVICE_CHECK(deviceMemcpy(_data, src._data, _nElements * sizeof(T), H2D));
        }
        else
        {
            DEVICE_CHECK(deviceMemcpy(_data, src._data, _nElements * sizeof(T), D2H));
        }
    }

    /** Converts a 2D buffer object between host and device backends, by creating a copy.
     *
     * @tparam BSource backend of source buffer object.
     * @param src the source buffer object.
     * @note This function takes care of host-to-device and device-to-host
     * copy-constructions.
     */
    template <auto L_ = Layout1D,
              typename BSource,
              std::enable_if_t<!L_, bool> = true,
              typename                    = std::enable_if_t<!std::is_same_v<backend_type, BSource>>>
    explicit CBuffer(const CBuffer<T, BSource, NRows, NCols> &src) : _nRows{src._nRows}, _nColumns{src._nColumns}
    {
        // pitched allocation
        std::tie(_nPaddedColumns, _data) = mem::malloc<value_type, backend_type>(_nRows, _nColumns);
        _nElements                       = _nRows * _nPaddedColumns;

        if constexpr (mem::is_on_host_v<BSource>)
        {
            DEVICE_CHECK(deviceMemcpy2D(_data,
                                        _nPaddedColumns * sizeof(value_type), /* dpitch (bytes) */
                                        src._data,
                                        src._nPaddedColumns * sizeof(value_type), /* spitch (bytes) */
                                        _nColumns * sizeof(value_type),           /* width (bytes) */
                                        _nRows,
                                        H2D));
        }
        else
        {
            DEVICE_CHECK(deviceMemcpy2D(_data,
                                        _nPaddedColumns * sizeof(value_type), /* dpitch (bytes) */
                                        src._data,
                                        src._nPaddedColumns * sizeof(value_type), /* spitch (bytes) */
                                        _nColumns * sizeof(value_type),           /* width (bytes) */
                                        _nRows,
                                        D2H));
        }
    }
    /**@}*/

    /** @{ Host-to-device and device-to-host converting assignment operators */
    /** Creates a 1D buffer object by copy-assignment of other 1D buffer object.
     *
     * @tparam BSource backend of source buffer object.
     * @param src the source buffer object.
     * @note This overload takes care of host-to-device and device-to-host
     * copy-assignments.
     */
    template <auto L_ = Layout1D,
              typename BSource,
              std::enable_if_t<L_, bool> = true,
              typename                   = std::enable_if_t<!std::is_same_v<backend_type, BSource>>>
    auto
    operator=(const CBuffer<T, BSource, NRows, NCols> &src) -> CBuffer<T, B, NRows, NCols> &
    {
        if (&src != this)
        {
            _nRows          = src._nRows;
            _nColumns       = src._nColumns;
            _nPaddedColumns = src._nPaddedColumns;
            _nElements      = src._nElements;

            _data = mem::malloc<value_type, B>(_nElements);

            if constexpr (mem::is_on_host_v<BSource>)
            {
                DEVICE_CHECK(deviceMemcpy(_data, src._data, _nElements * sizeof(T), H2D));
            }
            else
            {
                DEVICE_CHECK(deviceMemcpy(_data, src._data, _nElements * sizeof(T), D2H));
            }
        }

        return *this;
    }

    /** Creates a 2D buffer object by copy-assignment of other 2D buffer object.
     *
     * @tparam BSource backend of source buffer object.
     * @param src the source buffer object.
     * @note This overload takes care of host-to-device and device-to-host
     * copy-assignments.
     */
    template <auto L_ = Layout1D,
              typename BSource,
              std::enable_if_t<!L_, bool> = true,
              typename                    = std::enable_if_t<!std::is_same_v<backend_type, BSource>>>
    auto
    operator=(const CBuffer<T, BSource, NRows, NCols> &src) -> CBuffer<T, B, NRows, NCols> &
    {
        if (&src != this)
        {
            _nRows    = src._nRows;
            _nColumns = src._nColumns;

            // pitched allocation
            std::tie(_nPaddedColumns, _data) = mem::malloc<value_type, B>(_nRows, _nColumns);
            _nElements                       = _nRows * _nPaddedColumns;

            if constexpr (mem::is_on_host_v<BSource>)
            {
                DEVICE_CHECK(deviceMemcpy2D(_data,
                                            _nPaddedColumns * sizeof(value_type), /* dpitch (bytes) */
                                            src._data,
                                            src._nPaddedColumns * sizeof(value_type), /* spitch (bytes) */
                                            _nColumns * sizeof(value_type),           /* width (bytes) */
                                            _nRows,
                                            H2D));
            }
            else
            {
                DEVICE_CHECK(deviceMemcpy2D(_data,
                                            _nPaddedColumns * sizeof(value_type), /* dpitch (bytes) */
                                            src._data,
                                            src._nPaddedColumns * sizeof(value_type), /* spitch (bytes) */
                                            _nColumns * sizeof(value_type),           /* width (bytes) */
                                            _nRows,
                                            D2H));
            }
        }

        return *this;
    }
    /**@}*/

    /** @{ Accessors */
    /** Return number of rows (height) of the buffer. */
    __host__ __device__ [[nodiscard]] auto
             nRows() const -> size_type
    {
        return _nRows;
    }

    /** Return number of columns (width) of the buffer. */
    __host__ __device__ [[nodiscard]] auto
             nColumns() const -> size_type
    {
        return _nColumns;
    }

    /** Return number of padded columns (pitch) of the buffer. */
    __host__ __device__ [[nodiscard]] auto
             nPaddedColumns() const -> size_type
    {
        return _nPaddedColumns;
    }

    /** Return unpadded number of elements in buffer. */
    __host__ __device__ [[nodiscard]] auto
             size() const -> size_type
    {
        return _nRows * _nColumns;
    }

    /** Return padded number of elements in buffer. */
    __host__ __device__ [[nodiscard]] auto
             paddedSize() const -> size_type
    {
        return _nElements;
    }

    /** Get element in 1D buffer
     *
     * @param[in] i index.
     * @note Access is not bounds-checked at run-time!
     */
    template <auto L_ = Layout1D, std::enable_if_t<L_, bool> = true>
    __host__ __device__ auto
    operator()(size_type i) -> reference
    {
        return _data[i];
    }

    /** Get element in 1D buffer
     *
     * @param[in] i index.
     * @note Access is not bounds-checked at run-time!
     */
    template <auto L_ = Layout1D, std::enable_if_t<L_, bool> = true>
    __host__ __device__ auto
    operator()(size_type i) const -> const_reference
    {
        return _data[i];
    }

    /** Get element in 2D buffer
     *
     * @param[in] i row index.
     * @param[in] j column index.
     * @note Access is not bounds-checked at run-time!
     */
    template <auto L_ = Layout1D, std::enable_if_t<!L_, bool> = true>
    __host__ __device__ auto
    operator()(size_type i, size_type j) -> reference
    {
        return _data[i * _nPaddedColumns + j];
    }

    /** Get element in 2D buffer
     *
     * @param[in] i row index.
     * @param[in] j column index.
     * @note Access is not bounds-checked at run-time!
     */
    template <auto L_ = Layout1D, std::enable_if_t<!L_, bool> = true>
    __host__ __device__ auto
    operator()(size_type i, size_type j) const -> const_reference
    {
        return _data[i * _nPaddedColumns + j];
    }

    /** Get element in 1D buffer
     *
     * @param[in] i index.
     * @note Access is bounds-checked at run-time!
     */
    template <auto L_ = Layout1D, std::enable_if_t<L_, bool> = true>
    __host__ __device__ auto
    at(size_type i) -> reference
    {
        errors::assertMsgCritical(i < _nColumns, std::string(__func__) + ": you cannot access an element beyond " + std::to_string(_nColumns));
        return _data[i];
    }

    /** Get element in 1D buffer
     *
     * @param[in] i index.
     * @note Access is bounds-checked at run-time!
     */
    template <auto L_ = Layout1D, std::enable_if_t<L_, bool> = true>
    __host__ __device__ auto
    at(size_type i) const -> const_reference
    {
        errors::assertMsgCritical(i < _nColumns, std::string(__func__) + ": you cannot access an element beyond " + std::to_string(_nColumns));
        return _data[i];
    }

    /** Get element in 2D buffer
     *
     * @param[in] i row index.
     * @param[in] j column index.
     * @note Access is bounds-checked at run-time!
     */
    template <auto L_ = Layout1D, std::enable_if_t<!L_, bool> = true>
    __host__ __device__ auto
    at(size_type i, size_type j) -> reference
    {
        errors::assertMsgCritical(i < _nRows, std::string(__func__) + ": you cannot access an element beyond row " + std::to_string(_nRows));
        errors::assertMsgCritical(i < _nColumns, std::string(__func__) + ": you cannot access an element beyond column " + std::to_string(_nColumns));
        return _data[i * _nPaddedColumns + j];
    }

    /** Get element in 2D buffer
     *
     * @param[in] i row index.
     * @param[in] j column index.
     * @note Access is bounds-checked at run-time!
     */
    template <auto L_ = Layout1D, std::enable_if_t<!L_, bool> = true>
    __host__ __device__ auto
    at(size_type i, size_type j) const -> const_reference
    {
        errors::assertMsgCritical(i < _nRows, std::string(__func__) + ": you cannot access an element beyond row " + std::to_string(_nRows));
        errors::assertMsgCritical(i < _nColumns, std::string(__func__) + ": you cannot access an element beyond column " + std::to_string(_nColumns));
        return _data[i * _nPaddedColumns + j];
    }
    /**@}*/

    /** @{ Pointer accesssors. */
    /** Return pointer to buffer data. */
    __host__ __device__ auto
    data() -> pointer
    {
        return _data;
    }

    /** @param irow the index of row. */
    /** Return pointer to buffer data. */
    __host__ __device__ [[nodiscard]] auto
             data(const size_type irow) -> pointer
    {
        return &_data[irow * _nPaddedColumns];
    }

    /** Return pointer to buffer data. */
    __host__ __device__ [[nodiscard]] auto
             data() const -> const_pointer
    {
        return _data;
    }

    /** @param irow the index of row. */
    /** Return pointer to buffer data. */
    __host__ __device__ [[nodiscard]] auto
             data(const size_type irow) const -> const_pointer
    {
        return &_data[irow * _nPaddedColumns];
    }

    /**@}*/

    /** @{ Iterators */
    // TODO
    /**@}*/

    /** @{ mdspan interface */
    template <typename W>
    struct restrict_accessor
    {
        using element_type     = W;
        using data_handle_type = W *__restrict;
        using reference        = W &;
        auto
        access(data_handle_type p, ptrdiff_t i) const noexcept -> reference
        {
            return p[i];
        }
        auto
        offset(data_handle_type p, ptrdiff_t i) const noexcept -> data_handle_type
        {
            return p + i;
        }
    };

    using extents_type = metautils::select_t<
        // 1D buffer with run-time number of elements.
        metautils::condition<(kind == Kind::X), stdex::extents<size_type, Dynamic>>,
        // 1D buffer with compile-time number of elements.
        metautils::condition<(kind == Kind::N), stdex::extents<size_type, NCols>>,
        // 2D buffer with run-time number of rows and columns.
        metautils::condition<(kind == Kind::XY), stdex::extents<size_type, Dynamic, Dynamic>>,
        // 2D buffer with compile-time number of rows.
        metautils::condition<(kind == Kind::MY), stdex::extents<size_type, NRows, Dynamic>>,
        // 2D buffer with compile-time number of columns.
        metautils::condition<(kind == Kind::XN), stdex::extents<size_type, Dynamic, NCols>>,
        // 2D buffer with compile-time number of rows and columns.
        metautils::condition<(kind == Kind::MN), stdex::extents<size_type, NRows, NCols>>,
        // fallback type
        stdex::mdspan<T, stdex::extents<size_type, Dynamic, Dynamic>>>;

    using mdspan_type = metautils::select_t<
        // 1D buffer with run-time number of elements.
        metautils::condition<(kind == Kind::X), stdex::mdspan<value_type, extents_type, stdex::layout_right, restrict_accessor<value_type>>>,
        // 1D buffer with compile-time number of elements.
        metautils::condition<(kind == Kind::N), stdex::mdspan<value_type, extents_type, stdex::layout_right, restrict_accessor<value_type>>>,
        // 2D buffer with run-time number of rows and columns.
        metautils::condition<(kind == Kind::XY), stdex::mdspan<value_type, extents_type, stdex::layout_stride, restrict_accessor<value_type>>>,
        // 2D buffer with compile-time number of rows.
        metautils::condition<(kind == Kind::MY), stdex::mdspan<value_type, extents_type, stdex::layout_stride, restrict_accessor<value_type>>>,
        // 2D buffer with compile-time number of columns.
        metautils::condition<(kind == Kind::XN), stdex::mdspan<value_type, extents_type, stdex::layout_stride, restrict_accessor<value_type>>>,
        // 2D buffer with compile-time number of rows and columns.
        metautils::condition<(kind == Kind::MN), stdex::mdspan<value_type, extents_type, stdex::layout_stride, restrict_accessor<value_type>>>,
        // fallback type
        stdex::mdspan<value_type, extents_type, stdex::layout_stride, restrict_accessor<value_type>>>;

    using const_mdspan_type = metautils::select_t<
        // 1D buffer with run-time number of elements.
        metautils::condition<(kind == Kind::X),
                             stdex::mdspan<const value_type, extents_type, stdex::layout_right, restrict_accessor<const value_type>>>,
        // 1D buffer with compile-time number of elements.
        metautils::condition<(kind == Kind::N),
                             stdex::mdspan<const value_type, extents_type, stdex::layout_right, restrict_accessor<const value_type>>>,
        // 2D buffer with run-time number of rows and columns.
        metautils::condition<(kind == Kind::XY),
                             stdex::mdspan<const value_type, extents_type, stdex::layout_stride, restrict_accessor<const value_type>>>,
        // 2D buffer with compile-time number of rows.
        metautils::condition<(kind == Kind::MY),
                             stdex::mdspan<const value_type, extents_type, stdex::layout_stride, restrict_accessor<const value_type>>>,
        // 2D buffer with compile-time number of columns.
        metautils::condition<(kind == Kind::XN),
                             stdex::mdspan<const value_type, extents_type, stdex::layout_stride, restrict_accessor<const value_type>>>,
        // 2D buffer with compile-time number of rows and columns.
        metautils::condition<(kind == Kind::MN),
                             stdex::mdspan<const value_type, extents_type, stdex::layout_stride, restrict_accessor<const value_type>>>,
        // fallback type
        stdex::mdspan<const value_type, extents_type, stdex::layout_stride, restrict_accessor<const value_type>>>;

    /** @{ Implicit conversions */
    /** Conversion to mdspan. */
    operator mdspan_type()
    {
        if constexpr (kind == buffer::Kind::X)
        {
            return mdspan_type(_data, _nColumns);
        }
        else if constexpr (kind == buffer::Kind::MY)
        {
            stdex::layout_stride::mapping layout{extents_type{_nColumns}, std::array<size_type, 2>{_nPaddedColumns, 1}};
            return mdspan_type(_data, layout);
        }
        else if constexpr (kind == buffer::Kind::XN)
        {
            stdex::layout_stride::mapping layout{extents_type{_nRows}, std::array<size_type, 2>{_nPaddedColumns, 1}};
            return mdspan_type(_data, layout);
        }
        else if constexpr (kind == buffer::Kind::XY)
        {
            stdex::layout_stride::mapping layout{extents_type{_nRows, _nColumns}, std::array<size_type, 2>{_nPaddedColumns, 1}};
            return mdspan_type(_data, layout);
        }
        else if constexpr (kind == buffer::Kind::MN)
        {
            stdex::layout_stride::mapping layout{extents_type{}, std::array<size_type, 2>{_nPaddedColumns, 1}};
            return mdspan_type(_data, layout);
        }
        else
        {
            return mdspan_type(_data);
        }
    }

    /** Conversion to const mdspan. */
    operator const_mdspan_type() const
    {
        if constexpr (kind == buffer::Kind::X)
        {
            return const_mdspan_type(_data, _nColumns);
        }
        else if constexpr (kind == buffer::Kind::MY)
        {
            stdex::layout_stride::mapping layout{extents_type{_nColumns}, std::array<size_type, 2>{_nPaddedColumns, 1}};
            return const_mdspan_type(_data, layout);
        }
        else if constexpr (kind == buffer::Kind::XN)
        {
            stdex::layout_stride::mapping layout{extents_type{_nRows}, std::array<size_type, 2>{_nPaddedColumns, 1}};
            return const_mdspan_type(_data, layout);
        }
        else if constexpr (kind == buffer::Kind::XY)
        {
            stdex::layout_stride::mapping layout{extents_type{_nRows, _nColumns}, std::array<size_type, 2>{_nPaddedColumns, 1}};
            return const_mdspan_type(_data, layout);
        }
        else if constexpr (kind == buffer::Kind::MN)
        {
            stdex::layout_stride::mapping layout{extents_type{}, std::array<size_type, 2>{_nPaddedColumns, 1}};
            return const_mdspan_type(_data, layout);
        }
        else
        {
            return const_mdspan_type(_data);
        }
    }
    /**@}*/

    auto
    getMDSpan() -> mdspan_type
    {
        return this->operator mdspan_type();
    }

    auto
    getMDSpan() const -> const_mdspan_type
    {
        return this->operator const_mdspan_type();
    }
    /**@}*/

    /** @{ Equality/inequality operators */
    /** Check approximate equality.
     *
     * @tparam Bother backend of right-hand side in comparison.
     * @param[in] lhs left-hand side of comparison.
     * @param[in] rhs right-hand side of comparison.
     *
     * @note We don't compare the pitch, since the padding of the rows might be
     * different between LHS and RHS, especially if the two are allocated on
     * different backends!  The number of elements might also be different, due
     * to the padding.
     * We compare the actual elements in `_data`, *i.e.* excluding the padding elements
     * which are, most likely, garbage!
     */
    template <typename Bother>
    friend auto
    operator==(const CBuffer<T, B, NRows, NCols> &lhs, const CBuffer<T, Bother, NRows, NCols> &rhs) -> bool
    {
        if (lhs._nElements != rhs._nElements)
        {
            return false;
        }

        if (lhs._nElements == 0)
        {
            return true;
        }

        if (lhs._nRows != rhs._nRows)
        {
            return false;
        }

        if (lhs._nColumns != rhs._nColumns)
        {
            return false;
        }

        for (size_type i = 0; i < lhs._nRows; ++i)
        {
            for (size_type j = 0; j < lhs._nColumns; ++j)
            {
                auto lhs_v = lhs._data[i * lhs._nPaddedColumns + j];
                auto rhs_v = rhs._data[i * rhs._nPaddedColumns + j];

                if (std::abs(lhs_v - rhs_v) > 1.0e-13)
                {
                    return false;
                }
            }
        }

        return true;
    }

    /** Check approximate inequality.
     *
     * @tparam Bother backend of right-hand side in comparison.
     * @param[in] lhs left-hand side of comparison.
     * @param[in] rhs right-hand side of comparison.
     */
    template <typename Bother>
    friend auto
    operator!=(const CBuffer<T, B, NRows, NCols> &lhs, const CBuffer<T, Bother, NRows, NCols> &rhs) -> bool
    {
        return !(lhs == rhs);
    }
    /**@}*/

    /** Set all elements in buffer to a given value.
     *
     * @tparam Extents variadic pack of extents dimensions.
     *
     * @param[in] fill_value value of fill element.
     * @param[in] extents new dimension of the extents.
     *
     * @note The padding elements are left uninitialized!
     */
    template <typename... Extents>
    auto
    setConstant(value_type fill_value, Extents... extents) -> void
    {
        resize(extents...);
        _fill(fill_value);
    }

    /** Creates a buffer with all elements equal to a given value.
     *
     * @tparam Extents variadic pack of extents dimensions.
     *
     * @param[in] fill_value value of fill element.
     * @param[in] extents new dimension of the extents.
     *
     * @note The padding elements are left uninitialized!
     */
    template <typename... Extents>
    static auto
    Constant(value_type fill_value, Extents... extents) -> CBuffer<T, B, NRows, NCols>
    {
        auto buf = CBuffer<T, B, NRows, NCols>{extents...};
        buf.setConstant(fill_value);
        return buf;
    }

    /** Set all elements in buffer to zero.
     *
     * @tparam Extents variadic pack of extents dimensions.
     *
     * @param[in] extents new dimension of the extents.
     *
     * @note The padding elements are left uninitialized!
     */
    template <typename... Extents>
    auto
    setZero(Extents... extents) -> void
    {
        resize(extents...);
        _fill(value_type{0});
    }

    /** Creates a buffer with zero elements.
     *
     * @tparam Extents variadic pack of extents dimensions.
     *
     * @param[in] extents new dimension of the extents.
     *
     * @note The padding elements are left uninitialized!
     */
    template <typename... Extents>
    static auto
    Zero(Extents... extents) -> CBuffer<T, B, NRows, NCols>
    {
        auto buf = CBuffer<T, B, NRows, NCols>{extents...};
        buf.setZero();
        return buf;
    }

    /** Sets buffer elements to uniform random values in given inteval.
     *
     * @tparam Extents variadic pack of extents dimensions.
     *
     * @param[in] lower lower bound of interval.
     * @param[in] upper upper bound of interval.
     * @param[in] extents new dimension of the extents.
     *
     *
     * @note The padding elements are left uninitialized!
     *
     * The random numbers are generated on the CPU, then copied into the buffer.
     * This method uses the C++ default random engine with random seeding.
     * If you need more control, generate the random sequence (as std::vector or
     * raw array) and then use any of the available constructors.
     */
    template <typename... Extents>
    auto
    setRandom(value_type lower, value_type upper, Extents... extents) -> void
    {
        resize(extents...);
        _random_fill(lower, upper);
    }

    /** Creates a buffer with uniform random values in given inteval.
     *
     * @tparam Extents variadic pack of extents dimensions.
     *
     * @param[in] lower lower bound of interval.
     * @param[in] upper upper bound of interval.
     * @param[in] extents new dimension of the extents.
     *
     * @note The padding elements are left uninitialized!
     *
     * The random numbers are generated on the CPU, then copied into the buffer.
     * This method uses the C++ default random engine with random seeding.
     * If you need more control, generate the random sequence (as std::vector or
     * raw array) and then use any of the available constructors.
     */
    template <typename... Extents>
    static auto
    Random(value_type lower, value_type upper, Extents... extents) -> CBuffer<T, B, NRows, NCols>
    {
        auto buf = CBuffer<T, B, NRows, NCols>{extents...};
        buf.setRandom(lower, upper);
        return buf;
    }

    /** Clears contents of buffer.
     */
    auto
    clear() -> void
    {
        // deallocate existing _data
        if (_data)
        {
            mem::free<value_type, backend_type>(_data);
        }
    }

    /** Checks if buffer is empty.
     *
     * @return True if buffer is empty, false otherwise.
     */
    [[nodiscard]] auto
    empty() const -> bool
    {
        return (_nElements == 0);
    }
};
}  // namespace buffer

#endif /* BufferImpl_hpp */
