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
#include <tuple>
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

/** Class to manage memory buffer allocation, manipulation, and deallocation.
 *
 * @tparam T scalar type of buffer elements. Must be an arithmetic type.
 * @tparam B backend of buffer allocation.
 * @tparam NRows number of rows at compile-time.
 * @tparam NColumns number of columns at compile-time.
 * @author R. Di Remigio
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

    static constexpr bool Layout1D           = (NRows == 1);
    static constexpr bool CompileTimeRows    = (NRows != Dynamic);
    static constexpr bool CompileTimeColumns = (NCols != Dynamic);

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

    static constexpr Kind kind = [] {
        if constexpr (Layout1D && !CompileTimeColumns)
        {
            return Kind::X;
        }
        else if constexpr (Layout1D && CompileTimeColumns)
        {
            return Kind::N;
        }
        else if constexpr (!CompileTimeRows && !CompileTimeColumns)
        {
            return Kind::XY;
        }
        else if constexpr (CompileTimeRows && !CompileTimeColumns)
        {
            return Kind::MY;
        }
        else if constexpr (!CompileTimeRows && CompileTimeColumns)
        {
            return Kind::XN;
        }
        else
        {
            return Kind::MN;
        }
    }();

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

    /** TODO Write Doxygen string
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

    /** TODO Write Doxygen string
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

    /** TODO Write Doxygen string
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

    /** TODO Write Doxygen string
     */
    auto
    _allocate() -> void
    {
        if (_data)
        {
            errors::msgCritical(std::string(__func__) + ": buffer already allocated!");
        }

        if constexpr (Layout1D)
        {
            _nPaddedColumns = mem::get_pitch<value_type>(Alignment, _nColumns);
            _nElements      = _nRows * _nPaddedColumns;
            _data           = mem::malloc<value_type, backend_type>(_nElements);
        }
        else
        {
            std::tie(_nPaddedColumns, _data) = mem::malloc<value_type, backend_type>(_nRows, _nColumns);
            _nElements                       = _nRows * _nPaddedColumns;
        }
    }

    /** TODO Write Doxygen string */
    auto
    _resize(size_type nRows, size_type nCols) -> void
    {
        // deallocate existing _data
        if (_data)
        {
            mem::free<value_type, backend_type>(_data);
        }

        // reset number of rows
        _nRows = nRows;

        // reset number of columns
        _nColumns = nCols;

        // reallocate
        _allocate();
    }

   public:
    /** @{ Default CTORs. */
    /** Default CTOR when some dimensions are not known at compile-time.
     *
     * @warning These CTORs **do not** allocate memory. You will have to call
     * `resize` to do so or, alternatively one of
     * `setConstant`/`setZero`/`setRandom` to allocate *and* initialize.
     */
    template <auto K_ = kind, std::enable_if_t<!(K_ == Kind::N || K_ == Kind::MN), bool> = true>
    CBuffer()
    {
    }

    /** Default CTOR for 1D and 2D buffers with compile-time number of elements in each dimension. */
    template <auto K_ = kind, std::enable_if_t<(K_ == Kind::N || K_ == Kind::MN), bool> = true>
    CBuffer() : _nElements{_nRows * _nPaddedColumns}, _data{mem::malloc<value_type, backend_type>(_nElements)}
    {
    }
    /** @} */

    /** @{ CTORs from dimensions.
     *
     * @note Only defined when some of the buffer dimensions are not known at
     * compile-time.
     * @note These CTORs **allocate** memory, but leave uninitialized.
     * You can call one of `setConstant`/`setZero`/`setRandom` to do so.
     */
    /** CTOR for 1D buffer with run-time number of elements.
     *
     * @param[in] nCols number of elements in 1D buffer
     */
    template <auto K_ = kind, std::enable_if_t<(K_ == Kind::X), bool> = true>
    explicit CBuffer(size_type nCols)
        : _nColumns{nCols}
        , _nPaddedColumns{mem::get_pitch<value_type>(Alignment, _nColumns)}
        , _nElements{_nRows * _nPaddedColumns}
        , _data{mem::malloc<value_type, backend_type>(_nElements)}
    {
    }

    /** CTOR for 2D buffer with run-time number of rows and columns.
     *
     * @param[in] nRows number of rows (height).
     * @param[in] nCols number of columns (width).
     */
    template <auto K_ = kind, std::enable_if_t<(K_ == Kind::XY), bool> = true>
    CBuffer(size_type nRows, size_type nCols) : _nRows{nRows}, _nColumns{nCols}
    {
        _allocate();
    }

    /** CTOR for 2D buffer with compile-time number of rows.
     *
     * @param[in] nCols number of columns (width).
     */
    template <auto K_ = kind, std::enable_if_t<(K_ == Kind::MY), bool> = true>
    explicit CBuffer(size_type nCols) : _nColumns{nCols}
    {
        _allocate();
    }

    /** CTOR for 2D buffer with compile-time number of columns.
     *
     * @param[in] nRows number of rows (height).
     */
    template <auto K_ = kind, std::enable_if_t<(K_ == Kind::XN), bool> = true>
    explicit CBuffer(size_type nRows) : _nRows{nRows}
    {
        _allocate();
    }
    /** @} */

    /** @{ CTORs from unaligned data buffers.
     *
     * @note These CTORs **allocate** memory *and* initialize it with the data
     * in the unaligned buffer, _i.e._ a `std::vector` or a raw pointer.
     */
    /** CTOR for 1D buffer with run-time number of elements.
     *
     * @param[in] nCols number of elements in 1D buffer.
     * @param[in] v raw array of data to copy into 1D buffer.
     * @warning the input array **must** have at least `nCols` elements!
     */
    template <auto K_ = kind, std::enable_if_t<(K_ == Kind::X), bool> = true>
    CBuffer(size_type nCols, const_pointer v) : CBuffer<value_type, backend_type, 1, Dynamic>{nCols}
    {
        _copy_unaligned(v);
    }

    /** CTOR for 1D buffer with run-time number of elements.
     *
     * @param[in] nCols number of elements in 1D buffer.
     * @param[in] v vector of data to copy into 1D buffer.
     * @warning the input vector **must** have at least `nCols` elements!
     */
    template <auto K_ = kind, std::enable_if_t<(K_ == Kind::X), bool> = true>
    CBuffer(size_type nCols, const std::vector<value_type> &v) : CBuffer<value_type, backend_type, 1, Dynamic>{nCols, v.data()}
    {
        errors::assertMsgCritical(v.size() <= _nColumns,
                                  std::string(__func__) + ": input vector must have at least " + std::to_string(_nColumns) + " elements");
    }

    /** CTOR for 1D buffer with run-time number of elements.
     *
     * @param[in] v vector of data to copy into 1D buffer
     * @note This CTOR takes the dimension from `std::vector`.
     */
    template <auto K_ = kind, std::enable_if_t<(K_ == Kind::X), bool> = true>
    explicit CBuffer(const std::vector<value_type> &v) : CBuffer<value_type, backend_type, 1, Dynamic>{v.size(), v.data()}
    {
    }

    /** CTORs for 1D buffer with compile-time number of elements.
     *
     * @param[in] v raw array of data to copy into 1D buffer.
     * @warning the input array **must** have at least `NCols` elements!
     */
    template <auto K_ = kind, std::enable_if_t<(K_ == Kind::N), bool> = true>
    explicit CBuffer(const_pointer v) : CBuffer<value_type, backend_type, 1, NCols>{}
    {
        _copy_unaligned(v);
    }

    /** CTORs for 1D buffer with compile-time number of elements.
     *
     * @param[in] v vector of data to copy into 1D buffer
     * @warning the input vector **must** have at least `NCols` elements!
     */
    template <auto K_ = kind, std::enable_if_t<(K_ == Kind::N), bool> = true>
    explicit CBuffer(const std::vector<value_type> &v) : CBuffer<value_type, backend_type, 1, NCols>{v.data()}
    {
    }

    /** CTOR for 2D buffer with run-time number of rows and columns.
     *
     * @param[in] nRows number of rows (height).
     * @param[in] nCols number of columns (width).
     * @param[in] v raw array of data to copy into buffer.
     */
    template <auto K_ = kind, std::enable_if_t<(K_ == Kind::XY), bool> = true>
    CBuffer(size_type nRows, size_type nCols, const_pointer v) : CBuffer<value_type, backend_type, Dynamic, Dynamic>{nRows, nCols}
    {
        _copy_unaligned(v);
    }

    /** CTOR for 2D buffer with run-time number of rows and columns.
     *
     * @param[in] nRows number of rows (height).
     * @param[in] nCols number of columns (width).
     * @param[in] v vector of data to copy into buffer.
     * @warning the input vector **must** have at least `nRows*nCols` elements!
     */
    template <auto K_ = kind, std::enable_if_t<(K_ == Kind::XY), bool> = true>
    CBuffer(size_type nRows, size_type nCols, const std::vector<value_type> &v)
        : CBuffer<value_type, backend_type, Dynamic, Dynamic>{nRows, nCols, v.data()}
    {
        errors::assertMsgCritical(v.size() <= _nRows * _nColumns,
                                  std::string(__func__) + ": input vector must have at least " + std::to_string(_nRows * _nColumns) + " elements");
    }

    /** CTOR for 2D buffer with compile-time number of rows.
     *
     * @param[in] nCols number of columns (width).
     * @param[in] v raw array of data to copy into buffer.
     * @warning the input array **must** have at least `nRows*nCols` elements!
     */
    template <auto K_ = kind, std::enable_if_t<(K_ == Kind::MY), bool> = true>
    CBuffer(size_type nCols, const_pointer v) : CBuffer<value_type, backend_type, NRows, Dynamic>{nCols}
    {
        _copy_unaligned(v);
    }

    /** CTOR for 2D buffer with compile-time number of rows.
     *
     * @param[in] nCols number of columns (width).
     * @param[in] v vector of data to copy into buffer.
     * @warning the input vector **must** have at least `nRows*nCols` elements!
     */
    template <auto K_ = kind, std::enable_if_t<(K_ == Kind::MY), bool> = true>
    CBuffer(size_type nCols, const std::vector<value_type> &v) : CBuffer<value_type, backend_type, NRows, Dynamic>{nCols, v.data()}
    {
        errors::assertMsgCritical(v.size() <= _nRows * _nColumns,
                                  std::string(__func__) + ": input vector must have at least " + std::to_string(_nRows * _nColumns) + " elements");
    }

    /** CTOR for 2D buffer with compile-time number of columns.
     *
     * @param[in] nRows number of rows (height).
     * @param[in] v raw array of data to copy into buffer.
     * @warning the input array **must** have at least `nRows*nCols` elements!
     */
    template <auto K_ = kind, std::enable_if_t<(K_ == Kind::XN), bool> = true>
    CBuffer(size_type nRows, const_pointer v) : CBuffer<value_type, backend_type, Dynamic, NCols>{nRows}
    {
        _copy_unaligned(v);
    }

    /** CTOR for 2D buffer with compile-time number of columns.
     *
     * @param[in] nRows number of rows (height).
     * @param[in] v vector of data to copy into buffer.
     * @warning the input vector **must** have at least `nRows*nCols` elements!
     */
    template <auto K_ = kind, std::enable_if_t<(K_ == Kind::XN), bool> = true>
    CBuffer(size_type nRows, const std::vector<T> &v) : CBuffer<value_type, backend_type, Dynamic, NCols>{nRows, v.data()}
    {
        errors::assertMsgCritical(v.size() <= _nRows * _nColumns,
                                  std::string(__func__) + ": input vector must have at least " + std::to_string(_nRows * _nColumns) + " elements");
    }

    /** CTOR for 2D buffer with compile-time number of rows and columns.
     *
     * @param[in] v raw array of data to copy into buffer.
     * @warning the input array **must** have at least `_nRows*_nColumns` elements!
     */
    template <auto K_ = kind, std::enable_if_t<(K_ == Kind::MN), bool> = true>
    explicit CBuffer(const_pointer v) : CBuffer<value_type, backend_type, NRows, NCols>{}
    {
        _copy_unaligned(v);
    }

    /** CTOR for 2D buffer with compile-time number of rows and columns.
     *
     * @param[in] v vector of data to copy into buffer.
     * @warning the input vector **must** have at least `_nRows*_nColumns` elements!
     */
    template <auto K_ = kind, std::enable_if_t<(K_ == Kind::MN), bool> = true>
    explicit CBuffer(const std::vector<T> &v) : CBuffer<value_type, backend_type, NRows, NCols>{v.data()}
    {
        errors::assertMsgCritical(v.size() <= _nRows * _nColumns,
                                  std::string(__func__) + ": input vector must have at least " + std::to_string(_nRows * _nColumns) + " elements");
    }
    /** @} */

    /** @{ Resize functions
     *
     * @note Available when some dimensions are not known at compile-time.
     * @note If the requested dimensions are equal to the current ones, it is a
     * no-op.
     * @warning Resizing is **always** non-conservative, _i.e._ the current
     * contents of the buffer  will be discarded.
     */
    /** Resize existing 1-dimensional buffer.
     *
     * @param[in] nCols number of elements in 1D buffer.
     *
     */
    template <auto K_ = kind, std::enable_if_t<(K_ == Kind::X), bool> = true>
    auto
    resize(size_type nCols) -> void
    {
        _resize(1, nCols);
    }

    /** Resize existing 2-dimensional buffer.
     *
     * @param[in] nRows number of rows (height).
     * @param[in] nCols number of columns (width).
     */
    template <auto K_ = kind, std::enable_if_t<(K_ == Kind::XY), bool> = true>
    auto
    resize(size_type nRows, size_type nCols) -> void
    {
        _resize(nRows, nCols);
    }

    /** Resize existing 2-dimensional buffer with compile-time number of rows.
     *
     * @param[in] nCols number of columns (width).
     */
    template <auto K_ = kind, std::enable_if_t<(K_ == Kind::MY), bool> = true>
    auto
    resize(size_type nCols) -> void
    {
        _resize(_nRows, nCols);
    }

    /** Resize existing 2-dimensional buffer with compile-time number of columns.
     *
     * @param[in] nRows number of rows (height).
     */
    template <auto K_ = kind, std::enable_if_t<(K_ == Kind::XN), bool> = true>
    auto
    resize(size_type nRows) -> void
    {
        _resize(nRows, _nColumns);
    }
    /** @} */

    friend class CBuffer<T, mem::Host, NRows, NCols>;
    friend class CBuffer<T, mem::Device, NRows, NCols>;

    /** @{ Copy CTOR. */
    /** Create a buffer object by copying other buffer object.
     *
     * @param source the source buffer object.
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
                                            source._data,
                                            source._nPaddedColumns * sizeof(value_type), /* spitch (bytes) */
                                            _nColumns * sizeof(value_type),              /* width (bytes) */
                                            _nRows,
                                            D2D));
            }
        }
    }
    /**@}*/

    /** @{ Host-to-device and device-to-host copies */
    /** Converts a 1D buffer object between host and device backends, by creating a copy.
     *
     * @tparam BSource backend of source buffer object.
     * @param source the source buffer object.
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
     * @param source the source buffer object.
     * @note This function takes care of host-to-device and device-to-host
     * copy-constructions.
     */
    template <auto L_ = Layout1D,
              typename BSource,
              std::enable_if_t<!L_, bool> = true,
              typename                    = std::enable_if_t<!std::is_same_v<backend_type, BSource>>>
    explicit CBuffer(const CBuffer<T, BSource, NRows, NCols> &source) : _nRows{source._nRows}, _nColumns{source._nColumns}
    {
        // pitched allocation
        std::tie(_nPaddedColumns, _data) = mem::malloc<value_type, backend_type>(_nRows, _nColumns);
        _nElements                       = _nRows * _nPaddedColumns;

        if constexpr (mem::is_on_host_v<BSource>)
        {
            DEVICE_CHECK(deviceMemcpy2D(_data,
                                        _nPaddedColumns * sizeof(value_type), /* dpitch (bytes) */
                                        source._data,
                                        source._nPaddedColumns * sizeof(value_type), /* spitch (bytes) */
                                        _nColumns * sizeof(value_type),              /* width (bytes) */
                                        _nRows,
                                        H2D));
        }
        else
        {
            DEVICE_CHECK(deviceMemcpy2D(_data,
                                        _nPaddedColumns * sizeof(value_type), /* dpitch (bytes) */
                                        source._data,
                                        source._nPaddedColumns * sizeof(value_type), /* spitch (bytes) */
                                        _nColumns * sizeof(value_type),              /* width (bytes) */
                                        _nRows,
                                        D2H));
        }
    }
    /**@}*/

    /** @{ Move CTORs. */
    /** Create a buffer object by moving other 1D buffer object.
     *
     * @param source the source buffer object.
     * @note Move-CTOR only applicable to host-side buffers.
     */
    template <typename BSource = B, typename = std::enable_if_t<!mem::is_on_device_v<BSource>>>
    CBuffer(CBuffer<T, B, NRows, NCols> &&source) noexcept
        : _nRows{source._nRows}
        , _nColumns{source._nColumns}
        , _nPaddedColumns{source._nPaddedColumns}
        , _nElements{source._nElements}
        , _data{source._data}
    {
        source._data = nullptr;
    }
    /**@}*/

    /** @{ Copy-assignment operators. */
    /** Create a 1D buffer object by copy-assignment of other 1D buffer object.
     *
     * @param source the source buffer object.
     * @note This overload takes care of host-to-host and device-to-device
     * copy-assignments.
     */
    template <auto L_ = Layout1D, std::enable_if_t<L_, bool> = true>
    auto
    operator=(const CBuffer<T, B, NRows, NCols> &source) -> CBuffer<T, B, NRows, NCols> &
    {
        if (this == &source) return *this;

        _nRows          = source._nRows;
        _nColumns       = source._nColumns;
        _nPaddedColumns = source._nPaddedColumns;
        _nElements      = source._nElements;

        _data = mem::malloc<value_type, B>(_nElements);

        if constexpr (mem::is_on_host_v<B>)
        {
            _copy_aligned(source._data);
        }
        else
        {
            DEVICE_CHECK(deviceMemcpy(_data, source._data, _nElements * sizeof(T), D2D));
        }

        return *this;
    }

    /** Creates a 1D buffer object by copy-assignment of other 1D buffer object.
     *
     * @tparam BSource backend of source buffer object.
     * @param source the source buffer object.
     * @note This overload takes care of host-to-device and device-to-host
     * copy-assignments.
     */
    template <auto L_ = Layout1D,
              typename BSource,
              std::enable_if_t<L_, bool> = true,
              typename                   = std::enable_if_t<!std::is_same_v<backend_type, BSource>>>
    auto
    operator=(const CBuffer<T, BSource, NRows, NCols> &source) -> CBuffer<T, B, NRows, NCols> &
    {
        if (this == &source) return *this;

        _nRows          = source._nRows;
        _nColumns       = source._nColumns;
        _nPaddedColumns = source._nPaddedColumns;
        _nElements      = source._nElements;

        _data = mem::malloc<value_type, B>(_nElements);

        if constexpr (mem::is_on_host_v<BSource>)
        {
            DEVICE_CHECK(deviceMemcpy(_data, source._data, _nElements * sizeof(T), H2D));
        }
        else
        {
            DEVICE_CHECK(deviceMemcpy(_data, source._data, _nElements * sizeof(T), D2H));
        }

        return *this;
    }

    /** Create a 2D buffer object by copy-assignment of other 2D buffer object.
     *
     * @param source the source buffer object.
     * @note This overload takes care of host-to-host and device-to-device
     * copy-assignments.
     */
    template <auto L_ = Layout1D, std::enable_if_t<!L_, bool> = true>
    auto
    operator=(const CBuffer<T, B, NRows, NCols> &source) -> CBuffer<T, B, NRows, NCols> &
    {
        if (this == &source) return *this;

        _nRows          = source._nRows;
        _nColumns       = source._nColumns;
        _nPaddedColumns = source._nPaddedColumns;
        _nElements      = source._nElements;

        _data = mem::malloc<value_type, backend_type>(_nElements);

        if constexpr (mem::is_on_host_v<B>)
        {
            _copy_aligned(source._data);
        }
        else
        {
            DEVICE_CHECK(deviceMemcpy2D(_data,
                                        _nPaddedColumns * sizeof(value_type), /* dpitch (bytes) */
                                        source._data,
                                        source._nPaddedColumns * sizeof(value_type), /* spitch (bytes) */
                                        _nColumns * sizeof(value_type),              /* width (bytes) */
                                        _nRows,
                                        D2D));
        }

        return *this;
    }

    /** Creates a 2D buffer object by copy-assignment of other 2D buffer object.
     *
     * @tparam BSource backend of source buffer object.
     * @param source the source buffer object.
     * @note This overload takes care of host-to-device and device-to-host
     * copy-assignments.
     */
    template <auto L_ = Layout1D,
              typename BSource,
              std::enable_if_t<!L_, bool> = true,
              typename                    = std::enable_if_t<!std::is_same_v<backend_type, BSource>>>
    auto
    operator=(const CBuffer<T, BSource, NRows, NCols> &source) -> CBuffer<T, B, NRows, NCols> &
    {
        if (this == &source) return *this;

        _nRows    = source._nRows;
        _nColumns = source._nColumns;

        // pitched allocation
        std::tie(_nPaddedColumns, _data) = mem::malloc<value_type, B>(_nRows, _nColumns);
        _nElements                       = _nRows * _nPaddedColumns;

        if constexpr (mem::is_on_host_v<BSource>)
        {
            DEVICE_CHECK(deviceMemcpy2D(_data,
                                        _nPaddedColumns * sizeof(value_type), /* dpitch (bytes) */
                                        source._data,
                                        source._nPaddedColumns * sizeof(value_type), /* spitch (bytes) */
                                        _nColumns * sizeof(value_type),              /* width (bytes) */
                                        _nRows,
                                        H2D));
        }
        else
        {
            DEVICE_CHECK(deviceMemcpy2D(_data,
                                        _nPaddedColumns * sizeof(value_type), /* dpitch (bytes) */
                                        source._data,
                                        source._nPaddedColumns * sizeof(value_type), /* spitch (bytes) */
                                        _nColumns * sizeof(value_type),              /* width (bytes) */
                                        _nRows,
                                        D2H));
        }

        return *this;
    }
    /**@}*/

    /** @{ Move-assignment operators. */
    /** Create a buffer object by move-assignment of other buffer object.
     *
     * @param source the source buffer object.
     * @note Move-assignment only applicable to host-side buffers.
     */
    template <typename BSource = B, typename = std::enable_if_t<!mem::is_on_device_v<BSource>>>
    auto
    operator=(CBuffer<T, B, NRows, NCols> &&source) noexcept -> CBuffer<T, B, NRows, NCols> &
    {
        if (this == &source) return *this;

        _nRows          = source._nRows;
        _nColumns       = source._nColumns;
        _nPaddedColumns = source._nPaddedColumns;
        _nElements      = source._nElements;

        mem::free<value_type, backend_type>(_data);

        _data = source._data;

        source._data = nullptr;

        return *this;
    }
    /**@}*/

    /** Destroys a memory buffer object. */
    ~CBuffer()
    {
        _nElements = 0;
        mem::free<value_type, backend_type>(_data);
    }

    /** @{ Accessors */
    /** Return number of rows (height) of the buffer. */
    __host__ __device__ auto
    nRows() const -> size_type
    {
        return _nRows;
    }

    /** Return number of columns (width) of the buffer. */
    __host__ __device__ auto
    nColumns() const -> size_type
    {
        return _nColumns;
    }

    /** Return number of padded columns (pitch) of the buffer. */
    __host__ __device__ auto
    nPaddedColumns() const -> size_type
    {
        return _nPaddedColumns;
    }

    /** Return unpadded number of elements in buffer. */
    __host__ __device__ auto
    size() const -> size_type
    {
        return _nRows * _nColumns;
    }

    /** Return padded number of elements in buffer. */
    __host__ __device__ auto
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
    __host__ __device__ auto
    data(const size_type irow) -> pointer
    {
        return &_data[irow * _nPaddedColumns];
    }

    /** Return pointer to buffer data. */
    __host__ __device__ auto
    data() const -> const_pointer
    {
        return _data;
    }
    
    /** @param irow the index of row. */
    /** Return pointer to buffer data. */
    __host__ __device__ auto
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
        using element_type = W;
        using pointer      = W *__restrict;
        using reference    = W &;
        reference
        access(pointer p, ptrdiff_t i) const noexcept
        {
            return p[i];
        }
        pointer
        offset(pointer p, ptrdiff_t i) const noexcept
        {
            return p + i;
        }
    };

    using extents_type = metautils::select_t<
        // 1D buffer with run-time number of elements.
        metautils::condition<(kind == Kind::X), stdex::extents<Dynamic>>,
        // 1D buffer with compile-time number of elements.
        metautils::condition<(kind == Kind::N), stdex::extents<NCols>>,
        // 2D buffer with run-time number of rows and columns.
        metautils::condition<(kind == Kind::XY), stdex::extents<Dynamic, Dynamic>>,
        // 2D buffer with compile-time number of rows.
        metautils::condition<(kind == Kind::MY), stdex::extents<NRows, Dynamic>>,
        // 2D buffer with compile-time number of columns.
        metautils::condition<(kind == Kind::XN), stdex::extents<Dynamic, NCols>>,
        // 2D buffer with compile-time number of rows and columns.
        metautils::condition<(kind == Kind::MN), stdex::extents<NRows, NCols>>,
        // fallback type
        stdex::mdspan<T, stdex::extents<Dynamic, Dynamic>>>;

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
    friend inline auto
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
    friend inline auto
    operator!=(const CBuffer<T, B, NRows, NCols> &lhs, const CBuffer<T, Bother, NRows, NCols> &rhs) -> bool
    {
        return !(lhs == rhs);
    }
    /**@}*/

    /** @{ Set all elements in buffer to a given value.
     *
     * @note The padding elements are left uninitialized!
     */
    /** 1D/2D buffer.
     *
     * @param v value of fill element.
     * @note This is always valid when all dimensions are fixed at compile-time.
     * When one or more dimensions are known at run-time, this method is only
     * valid after allocation of the buffer, _e.g._ after calling `resize` or by
     * creating the object with any of the allocating CTORs first.
     */
    auto
    setConstant(value_type fill_value) -> void
    {
        _fill(fill_value);
    }

    template <auto K_ = kind, std::enable_if_t<(K_ == Kind::N || K_ == Kind::MN), bool> = true>
    inline static auto
    Constant(value_type fill_value) -> CBuffer<T, B, NRows, NCols>
    {
        auto buf = CBuffer<T, B, NRows, NCols>{};
        buf.setConstant(fill_value);
        return buf;
    }

    /** 1D buffer.
     *
     * @param v value of fill element.
     * @param nCols number of columns.
     */
    template <auto K_ = kind, std::enable_if_t<(K_ == Kind::X), bool> = true>
    auto
    setConstant(value_type fill_value, size_type nCols) -> void
    {
        resize(nCols);
        _fill(fill_value);
    }

    template <auto K_ = kind, std::enable_if_t<(K_ == Kind::X), bool> = true>
    inline static auto
    Constant(value_type fill_value, size_type nCols) -> CBuffer<T, B, 1, NCols>
    {
        auto buf = CBuffer<T, B, 1, NCols>{nCols};
        buf.setConstant(fill_value);
        return buf;
    }

    /** 2D buffer of XY kind.
     *
     * @param v value of fill element.
     * @param nRows number of rows.
     * @param nCols number of columns.
     */
    template <auto K_ = kind, std::enable_if_t<(K_ == Kind::XY), bool> = true>
    auto
    setConstant(value_type fill_value, size_type nRows, size_type nCols) -> void
    {
        resize(nRows, nCols);
        _fill(fill_value);
    }

    template <auto K_ = kind, std::enable_if_t<(K_ == Kind::XY), bool> = true>
    inline static auto
    Constant(value_type fill_value, size_type nRows, size_type nCols) -> CBuffer<T, B, Dynamic, Dynamic>
    {
        auto buf = CBuffer<T, B, Dynamic, Dynamic>{nRows, nCols};
        buf.setConstant(fill_value);
        return buf;
    }

    /** 2D buffer of MY kind.
     *
     * @param v value of fill element.
     * @param nCols number of columns.
     */
    template <auto K_ = kind, std::enable_if_t<(K_ == Kind::MY), bool> = true>
    auto
    setConstant(value_type fill_value, size_type nCols) -> void
    {
        resize(nCols);
        _fill(fill_value);
    }

    template <auto K_ = kind, std::enable_if_t<(K_ == Kind::MY), bool> = true>
    inline static auto
    Constant(value_type fill_value, size_type nCols) -> CBuffer<T, B, NRows, Dynamic>
    {
        auto buf = CBuffer<T, B, NRows, Dynamic>{nCols};
        buf.setConstant(fill_value);
        return buf;
    }

    /** 2D buffer of XN kind.
     *
     * @param v value of fill element.
     * @param nRows number of rows.
     */
    template <auto K_ = kind, std::enable_if_t<(K_ == Kind::XN), bool> = true>
    auto
    setConstant(value_type fill_value, size_type nRows) -> void
    {
        resize(nRows);
        _fill(fill_value);
    }

    template <auto K_ = kind, std::enable_if_t<(K_ == Kind::XN), bool> = true>
    inline static auto
    Constant(value_type fill_value, size_type nRows) -> CBuffer<T, B, Dynamic, NCols>
    {
        auto buf = CBuffer<T, B, Dynamic, NCols>{nRows};
        buf.setConstant(fill_value);
        return buf;
    }
    /**@}*/

    /** @{ Set all elements in buffer to zero.
     *
     * @note The padding elements are left uninitialized!
     */
    /** 1D/2D buffer.
     *
     * @param v value of fill element.
     * @note This is always valid when all dimensions are fixed at compile-time.
     * When one or more dimensions are known at run-time, this method is only
     * valid after allocation of the buffer, _e.g._ after calling `resize` or by
     * creating the object with any of the allocating CTORs first.
     */
    auto
    setZero() -> void
    {
        setConstant(value_type{0});
    }

    template <auto K_ = kind, std::enable_if_t<(K_ == Kind::N || K_ == Kind::MN), bool> = true>
    inline static auto
    Zero() -> CBuffer<T, B, NRows, NCols>
    {
        auto buf = CBuffer<T, B, NRows, NCols>{};
        buf.setZero();
        return buf;
    }

    /** 1D buffer.
     *
     * @param nCols number of columns.
     */
    template <auto K_ = kind, std::enable_if_t<(K_ == Kind::X), bool> = true>
    auto
    setZero(size_type nCols) -> void
    {
        setConstant(value_type{0}, nCols);
    }

    template <auto K_ = kind, std::enable_if_t<(K_ == Kind::X), bool> = true>
    inline static auto
    Zero(size_type nCols) -> CBuffer<T, B, 1, NCols>
    {
        auto buf = CBuffer<T, B, 1, NCols>{nCols};
        buf.setZero();
        return buf;
    }

    /** 2D buffer of XY kind.
     *
     * @param nRows number of rows.
     * @param nCols number of columns.
     */
    template <auto K_ = kind, std::enable_if_t<(K_ == Kind::XY), bool> = true>
    auto
    setZero(size_type nRows, size_type nCols) -> void
    {
        setConstant(value_type{0}, nRows, nCols);
    }

    template <auto K_ = kind, std::enable_if_t<(K_ == Kind::XY), bool> = true>
    inline static auto
    Zero(size_type nRows, size_type nCols) -> CBuffer<T, B, Dynamic, Dynamic>
    {
        auto buf = CBuffer<T, B, Dynamic, Dynamic>{nRows, nCols};
        buf.setZero();
        return buf;
    }

    /** 2D buffer of MY kind.
     *
     * @param nCols number of columns.
     */
    template <auto K_ = kind, std::enable_if_t<(K_ == Kind::MY), bool> = true>
    auto
    setZero(size_type nCols) -> void
    {
        setConstant(value_type{0}, nCols);
    }

    template <auto K_ = kind, std::enable_if_t<(K_ == Kind::MY), bool> = true>
    inline static auto
    Zero(size_type nCols) -> CBuffer<T, B, NRows, Dynamic>
    {
        auto buf = CBuffer<T, B, NRows, Dynamic>{nCols};
        buf.setZero();
        return buf;
    }

    /** 2D buffer of XN kind.
     *
     * @param v value of fill element.
     * @param nRows number of rows.
     */
    template <auto K_ = kind, std::enable_if_t<(K_ == Kind::XN), bool> = true>
    auto
    setZero(size_type nRows) -> void
    {
        setConstant(value_type{0}, nRows);
    }

    template <auto K_ = kind, std::enable_if_t<(K_ == Kind::XN), bool> = true>
    inline static auto
    Zero(size_type nRows) -> CBuffer<T, B, Dynamic, NCols>
    {
        auto buf = CBuffer<T, B, Dynamic, NCols>{nRows};
        buf.setZero();
        return buf;
    }
    /**@}*/

    /** @{ Sets buffer elements to uniform random values in given inteval.
     *
     * @note The padding elements are left uninitialized!
     *
     * The random numbers are generated on the CPU, then copied into the buffer.
     * This method uses the C++ default random engine with random seeding.
     * If you need more control, generate the random sequence (as std::vector or
     * raw array) and then use any of the available constructors.
     */
    /** 1D/2D buffer.
     *
     * @param[in] lower lower bound of interval.
     * @param[in] upper upper bound of interval.
     * @note This is always valid when all dimensions are fixed at compile-time.
     * When one or more dimensions are known at run-time, this method is only
     * valid after allocation of the buffer, _e.g._ after calling `resize` or by
     * creating the object with any of the allocating CTORs first.
     */
    auto
    setRandom(value_type lower, value_type upper) -> void
    {
        _random_fill(lower, upper);
    }

    template <auto K_ = kind, std::enable_if_t<(K_ == Kind::N || K_ == Kind::MN), bool> = true>
    inline static auto
    Random(value_type upper, value_type lower) -> CBuffer<T, B, NRows, NCols>
    {
        auto buf = CBuffer<T, B, NRows, NCols>{};
        buf.setRandom(lower, upper);
        return buf;
    }

    /** 1D buffer
     *
     * @param[in] lower lower bound of interval.
     * @param[in] upper upper bound of interval.
     * @param nCols number of columns.
     */
    template <auto K_ = kind, std::enable_if_t<(K_ == Kind::X), bool> = true>
    auto
    setRandom(value_type lower, value_type upper, size_type nCols) -> void
    {
        resize(nCols);

        _random_fill(lower, upper);
    }

    template <auto K_ = kind, std::enable_if_t<(K_ == Kind::X), bool> = true>
    inline static auto
    Random(value_type upper, value_type lower, size_type nCols) -> CBuffer<T, B, 1, NCols>
    {
        auto buf = CBuffer<T, B, 1, NCols>{nCols};
        buf.setRandom(lower, upper);
        return buf;
    }

    /** 2D buffer of XY kind.
     *
     * @param[in] lower lower bound of interval.
     * @param[in] upper upper bound of interval.
     * @param nRows number of rows.
     * @param nCols number of columns.
     */
    template <auto K_ = kind, std::enable_if_t<(K_ == Kind::XY), bool> = true>
    auto
    setRandom(value_type lower, value_type upper, size_type nRows, size_type nCols) -> void
    {
        resize(nRows, nCols);

        _random_fill(lower, upper);
    }

    template <auto K_ = kind, std::enable_if_t<(K_ == Kind::XY), bool> = true>
    inline static auto
    Random(value_type upper, value_type lower, size_type nRows, size_type nCols) -> CBuffer<T, B, Dynamic, Dynamic>
    {
        auto buf = CBuffer<T, B, Dynamic, Dynamic>{nRows, nCols};
        buf.setRandom(lower, upper);
        return buf;
    }

    /** 2D buffer of MY kind.
     *
     * @param[in] lower lower bound of interval.
     * @param[in] upper upper bound of interval.
     * @param nCols number of columns.
     */
    template <auto K_ = kind, std::enable_if_t<(K_ == Kind::MY), bool> = true>
    auto
    setRandom(value_type lower, value_type upper, size_type nCols) -> void
    {
        resize(nCols);

        _random_fill(lower, upper);
    }

    template <auto K_ = kind, std::enable_if_t<(K_ == Kind::MY), bool> = true>
    inline static auto
    Random(value_type upper, value_type lower, size_type nCols) -> CBuffer<T, B, NRows, Dynamic>
    {
        auto buf = CBuffer<T, B, NRows, Dynamic>{nCols};
        buf.setRandom(lower, upper);
        return buf;
    }

    /** 2D buffer of XN kind.
     *
     * @param[in] lower lower bound of interval.
     * @param[in] upper upper bound of interval.
     * @param nRows number of rows.
     */
    template <auto K_ = kind, std::enable_if_t<(K_ == Kind::XN), bool> = true>
    auto
    setRandom(value_type lower, value_type upper, size_type nRows) -> void
    {
        resize(nRows);

        _random_fill(lower, upper);
    }

    template <auto K_ = kind, std::enable_if_t<(K_ == Kind::XN), bool> = true>
    inline static auto
    Random(value_type upper, value_type lower, size_type nRows) -> CBuffer<T, B, Dynamic, NCols>
    {
        auto buf = CBuffer<T, B, Dynamic, NCols>{nRows};
        buf.setRandom(lower, upper);
        return buf;
    }
    /** @} */
    
    /** Checks if buffer is empty.
     *
     * @return True if buffer is empty, false otherwise.
     */
    inline auto
    empty() const -> bool
    {
        return _nElements == 0;
    }
};
}  // namespace buffer

#endif /* BufferImpl_hpp */
