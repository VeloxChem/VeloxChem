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

   private:
    static constexpr bool Layout1D           = (NRows == 1);
    static constexpr bool CompileTimeRows    = (NRows != Dynamic);
    static constexpr bool CompileTimeColumns = (NCols != Dynamic);

    enum class Kind
    {
        X,
        N,
        XY,
        MY,
        XN,
        MN
    };

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
        for (auto i = 0; i < _nElements; ++i)
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
    template <auto K_ = kind, std::enable_if_t<(K_ == Kind::X || K_ == Kind::N), bool> = true>
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

    /** Copy 2-dimensional unaligned data to buffer.
     *
     * @param[in] src source memory.
     */
    template <auto K_ = kind, std::enable_if_t<!(K_ == Kind::X || K_ == Kind::N), bool> = true>
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

   public:
    /** Default CTOR when some dimensions are not known at compile-time.
     *
     * @warning This CTOR **does not** allocate memory.  You will have to call
     * one of `allocate`/`setConstant`/`setZero`/`setRandom` to do so.
     */
    template <auto K_ = kind, std::enable_if_t<!(K_ == Kind::N) || !(K_ == Kind::MN), bool> = true>
    explicit CBuffer()
    {
    }

    /** @{ CTORs 1D buffer with run-time number of elements. */
    /** CTOR
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

    /** CTOR
     *
     * @param[in] v vector of data to copy into 1D buffer
     */
    template <auto K_ = kind, std::enable_if_t<(K_ == Kind::X), bool> = true>
    explicit CBuffer(const std::vector<value_type> &v) : CBuffer<value_type, backend_type, 1, Dynamic>{v.size()}
    {
        _copy_unaligned(v.data());
    }

    /** CTOR
     *
     * @param[in] nCols number of elements in 1D buffer.
     * @param[in] v vector of data to copy into 1D buffer.
     * @warning the input vector **must** have at least `nCols` elements!
     */
    template <auto K_ = kind, std::enable_if_t<(K_ == Kind::X), bool> = true>
    explicit CBuffer(size_type nCols, const std::vector<value_type> &v) : CBuffer<value_type, backend_type, 1, Dynamic>{nCols}
    {
        errors::assertMsgCritical(v.size() <= _nColumns,
                                  std::string(__func__) + ": input vector must have at least " + std::to_string(_nColumns) + " elements");
        _copy_unaligned(v.data());
    }

    /** Resize existing buffer.
     *
     * @param[in] nCols number of elements in 1D buffer.
     *
     * @warning The resize is **always** non-conservative, _i.e._ the current
     * contents of the buffer  will be discarded.
     */
    template <auto K_ = kind, std::enable_if_t<(K_ == Kind::X), bool> = true>
    auto
    resize(size_type nCols) -> void
    {
        // deallocate existing _data
        if (_data)
        {
            mem::free<value_type, backend_type>(_data);
        }

        // reset number of columns
        _nColumns = nCols;
        // reallocate
        _allocate();
    }
    /** @} */

    /** @{ CTORs for 1D buffer with compile-time number of elements. */
    /** Default CTOR. */
    template <auto K_ = kind, std::enable_if_t<(K_ == Kind::N), bool> = true>
    CBuffer() : _nElements{_nRows * _nPaddedColumns}, _data{mem::malloc<value_type, backend_type>(_nElements)}
    {
    }

    /** CTOR
     *
     * @param[in] v vector of data to copy into 1D buffer
     * @warning the input vector **must** have at least `NCols` elements!
     */
    template <auto K_ = kind, std::enable_if_t<(K_ == Kind::N), bool> = true>
    explicit CBuffer(const std::vector<value_type> &v) : CBuffer<value_type, backend_type, 1, NCols>{}
    {
        errors::assertMsgCritical(v.size() <= _nColumns,
                                  std::string(__func__) + ": input vector must have at least " + std::to_string(_nColumns) + " elements");
        _copy_unaligned(v.data());
    }
    /** @} */

    /** @{ CTORs for 2D buffer with run-time number of rows and columns. */
    /** CTOR
     *
     * @param[in] nRows number of rows (height).
     * @param[in] nCols number of columns (width).
     */
    template <auto K_ = kind, std::enable_if_t<(K_ == Kind::XY), bool> = true>
    explicit CBuffer(size_type nRows, size_type nCols) : _nRows{nRows}, _nColumns{nCols}
    {
        std::tie(_nPaddedColumns, _data) = mem::malloc<value_type, backend_type>(_nRows, _nColumns);
        _nElements                       = _nRows * _nPaddedColumns;
    }

    /** CTOR
     *
     * @param[in] nRows number of rows (height).
     * @param[in] nCols number of columns (width).
     * @param[in] v vector of data to copy into buffer.
     * @warning the input vector **must** have at least `nRows*nCols` elements!
     */
    template <auto K_ = kind, std::enable_if_t<(K_ == Kind::XY), bool> = true>
    CBuffer(size_type nRows, size_type nCols, const std::vector<value_type> &v) : CBuffer<value_type, backend_type, Dynamic, Dynamic>{nRows, nCols}
    {
        errors::assertMsgCritical(v.size() <= _nRows * _nColumns,
                                  std::string(__func__) + ": input vector must have at least " + std::to_string(_nRows * _nColumns) + " elements");
        _copy_unaligned(v.data());
    }

    /** Allocate buffer.
     *
     * @param[in] nRows number of rows (height).
     * @param[in] nCols number of columns (width).
     * @warning it is invalid to call this method on an already allocated buffer!
     */
    template <auto K_ = kind, std::enable_if_t<(K_ == Kind::XY), bool> = true>
    auto
    allocate(size_type nRows, size_type nCols) -> void
    {
        if (!_data)
        {
            _nRows                           = nRows;
            _nColumns                        = nCols;
            std::tie(_nPaddedColumns, _data) = mem::malloc<value_type, backend_type>(_nRows, _nColumns);
            _nElements                       = _nRows * _nPaddedColumns;
        }
        else
        {
            errors::msgCritical(std::string(__func__) + ": 2D buffer already allocated!");
        }
    }
    /** @} */

    /** @{ CTORs for 2D buffer with compile-time number of rows. */
    /** CTOR
     *
     * @param[in] nCols number of columns (width).
     */
    template <auto K_ = kind, std::enable_if_t<(K_ == Kind::MY), bool> = true>
    explicit CBuffer(size_type nCols) : _nColumns{nCols}
    {
        std::tie(_nPaddedColumns, _data) = mem::malloc<value_type, backend_type>(_nRows, _nColumns);
        _nElements                       = _nRows * _nPaddedColumns;
    }

    /** CTOR
     *
     * @param[in] nCols number of columns (width).
     * @param[in] v vector of data to copy into buffer.
     * @warning the input vector **must** have at least `nRows*nCols` elements!
     */
    template <auto K_ = kind, std::enable_if_t<(K_ == Kind::MY), bool> = true>
    CBuffer(size_type nCols, const std::vector<value_type> &v) : CBuffer<value_type, backend_type, NRows, Dynamic>{nCols}
    {
        errors::assertMsgCritical(v.size() <= _nRows * _nColumns,
                                  std::string(__func__) + ": input vector must have at least " + std::to_string(_nRows * _nColumns) + " elements");
        _copy_unaligned(v.data());
    }

    /** Allocate buffer.
     *
     * @param[in] nCols number of columns (width).
     * @warning it is invalid to call this method on an already allocated buffer!
     */
    template <auto K_ = kind, std::enable_if_t<(K_ == Kind::MY), bool> = true>
    auto
    allocate(size_type nCols) -> void
    {
        if (!_data)
        {
            _nColumns                        = nCols;
            std::tie(_nPaddedColumns, _data) = mem::malloc<value_type, backend_type>(_nRows, _nColumns);
            _nElements                       = _nRows * _nPaddedColumns;
        }
        else
        {
            errors::msgCritical(std::string(__func__) + ": 2D buffer already allocated!");
        }
    }
    /** @} */

    /** @{ CTORs for 2D buffer with compile-time number of columns. */
    /** CTOR
     *
     * @param[in] nRows number of rows (height).
     */
    template <auto K_ = kind, std::enable_if_t<(K_ == Kind::XN), bool> = true>
    explicit CBuffer(size_type nRows) : _nRows{nRows}
    {
        std::tie(_nPaddedColumns, _data) = mem::malloc<value_type, backend_type>(_nRows, _nColumns);
        _nElements                       = _nRows * _nPaddedColumns;
    }

    /** CTOR
     *
     * @param[in] nRows number of rows (height).
     * @param[in] v vector of data to copy into buffer.
     * @warning the input vector **must** have at least `nRows*nCols` elements!
     */
    template <auto K_ = kind, std::enable_if_t<(K_ == Kind::XN), bool> = true>
    CBuffer(size_type nRows, const std::vector<T> &v) : CBuffer<value_type, backend_type, Dynamic, NCols>{nRows}
    {
        errors::assertMsgCritical(v.size() <= _nRows * _nColumns,
                                  std::string(__func__) + ": input vector must have at least " + std::to_string(_nRows * _nColumns) + " elements");
        _copy_unaligned(v.data());
    }

    /** Allocate buffer.
     *
     * @param[in] nRows number of rows (height).
     * @warning it is invalid to call this method on an already allocated buffer!
     */
    template <auto K_ = kind, std::enable_if_t<(K_ == Kind::XN), bool> = true>
    auto
    allocate(size_type nRows) -> void
    {
        if (!_data)
        {
            _nRows                           = nRows;
            std::tie(_nPaddedColumns, _data) = mem::malloc<value_type, backend_type>(_nRows, _nColumns);
            _nElements                       = _nRows * _nPaddedColumns;
        }
        else
        {
            errors::msgCritical(std::string(__func__) + ": 2D buffer already allocated!");
        }
    }
    /** @} */

    /** @{ CTORs for 2D buffer with compile-time number of rows and columns*/
    /** Default CTOR. */
    template <auto K_ = kind, std::enable_if_t<(K_ == Kind::MN), bool> = true>
    CBuffer()
    {
        std::tie(_nPaddedColumns, _data) = mem::malloc<value_type, backend_type>(_nRows, _nColumns);
        _nElements                       = _nRows * _nPaddedColumns;
    }

    /** CTOR
     *
     * @param[in] v vector of data to copy into buffer.
     * @warning the input vector **must** have at least `_nRows*_nColumns` elements!
     */
    template <auto K_ = kind, std::enable_if_t<(K_ == Kind::MN), bool> = true>
    explicit CBuffer(const std::vector<T> &v) : CBuffer<value_type, backend_type, NRows, NCols>{}
    {
        errors::assertMsgCritical(v.size() <= _nRows * _nColumns,
                                  std::string(__func__) + ": input vector must have at least " + std::to_string(_nRows * _nColumns) + " elements");
        _copy_unaligned(v.data());
    }
    /**@}*/

    friend class CBuffer<T, mem::Host, NRows, NCols>;
    friend class CBuffer<T, mem::Device, NRows, NCols>;

    /** @{ Copy CTORs. */
    /** Create a 1D buffer object by copying other 1D buffer object.
     *
     * @param source the source buffer object.
     * @note This overload takes care of host-to-host and device-to-device
     * copy-constructions.
     */
    template <auto K_ = kind, std::enable_if_t<(K_ == Kind::X || K_ == Kind::N), bool> = true>
    CBuffer(const CBuffer<T, B, NRows, NCols> &source)
        : _nRows{source._nRows}
        , _nColumns{source._nColumns}
        , _nPaddedColumns{source._nPaddedColumns}
        , _nElements{source._nElements}
        , _data{mem::malloc<value_type, backend_type>(_nElements)}
    {
        if constexpr (mem::is_on_host_v<B>)
        {
            _copy_aligned(source._data);
        }
        else
        {
            DEVICE_CHECK(deviceMemcpy(_data, source._data, _nElements * sizeof(T), D2D));
        }
    }

    /** Creates a 1D buffer object by copying other 1D buffer object.
     *
     * @tparam BSource backend of source buffer object.
     * @param source the source buffer object.
     * @note This overload takes care of host-to-device and device-to-host
     * copy-constructions.
     */
    template <auto K_ = kind,
              typename BSource,
              std::enable_if_t<(K_ == Kind::X || K_ == Kind::N), bool> = true,
              typename                                                 = std::enable_if_t<!std::is_same_v<backend_type, BSource>>>
    CBuffer(const CBuffer<T, BSource, NRows, NCols> &source)
        : _nRows{source._nRows}
        , _nColumns{source._nColumns}
        , _nPaddedColumns{source._nPaddedColumns}
        , _nElements{source._nElements}
        , _data{mem::malloc<value_type, backend_type>(_nElements)}
    {
        if constexpr (mem::is_on_host_v<BSource>)
        {
            DEVICE_CHECK(deviceMemcpy(_data, source._data, _nElements * sizeof(T), H2D));
        }
        else
        {
            DEVICE_CHECK(deviceMemcpy(_data, source._data, _nElements * sizeof(T), D2H));
        }
    }

    /** Create a 2D buffer object by copying other 2D buffer object.
     *
     * @param source the source buffer object.
     * @note This overload takes care of host-to-host and device-to-device
     * copy-constructions.
     */
    template <auto K_ = kind, std::enable_if_t<!(K_ == Kind::X || K_ == Kind::N), bool> = true>
    CBuffer(const CBuffer<T, B, NRows, NCols> &source)
        : _nRows{source._nRows}
        , _nColumns{source._nColumns}
        , _nPaddedColumns{source._nPaddedColumns}
        , _nElements{source._nElements}
        , _data{mem::malloc<value_type, backend_type>(_nElements)}
    {
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
    }

    /** Creates a 2D buffer object by copying other 2D buffer object.
     *
     * @tparam BSource backend of source buffer object.
     * @param source the source buffer object.
     * @note This overload takes care of host-to-device and device-to-host
     * copy-constructions.
     */
    template <auto K_ = kind,
              typename BSource,
              std::enable_if_t<!(K_ == Kind::X || K_ == Kind::N), bool> = true,
              typename                                                  = std::enable_if_t<!std::is_same_v<backend_type, BSource>>>
    CBuffer(const CBuffer<T, BSource, NRows, NCols> &source) : _nRows{source._nRows}, _nColumns{source._nColumns}
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
    CBuffer(const CBuffer<T, B, NRows, NCols> &&source) noexcept
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
    template <auto K_ = kind, std::enable_if_t<(K_ == Kind::X || K_ == Kind::N), bool> = true>
    CBuffer<T, B, NRows, NCols> &
    operator=(const CBuffer<T, B, NRows, NCols> &source)
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
    template <auto K_ = kind,
              typename BSource,
              std::enable_if_t<(K_ == Kind::X || K_ == Kind::N), bool> = true,
              typename                                                 = std::enable_if_t<!std::is_same_v<backend_type, BSource>>>
    CBuffer<T, B, NRows, NCols> &
    operator=(const CBuffer<T, BSource, NRows, NCols> &source)
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
    template <auto K_ = kind, std::enable_if_t<!(K_ == Kind::X || K_ == Kind::N), bool> = true>
    CBuffer<T, B, NRows, NCols> &
    operator=(const CBuffer<T, B, NRows, NCols> &source)
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
    template <auto K_ = kind,
              typename BSource,
              std::enable_if_t<!(K_ == Kind::X || K_ == Kind::N), bool> = true,
              typename                                                  = std::enable_if_t<!std::is_same_v<backend_type, BSource>>>
    CBuffer<T, B, NRows, NCols> &
    operator=(const CBuffer<T, BSource, NRows, NCols> &source)
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
    CBuffer<T, B, NRows, NCols> &
    operator=(const CBuffer<T, B, NRows, NCols> &&source) noexcept
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
    template <auto K_ = kind, std::enable_if_t<(K_ == Kind::X || K_ == Kind::N), bool> = true>
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
    template <auto K_ = kind, std::enable_if_t<(K_ == Kind::X || K_ == Kind::N), bool> = true>
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
    template <auto K_ = kind, std::enable_if_t<!(K_ == Kind::X || K_ == Kind::N), bool> = true>
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
    template <auto K_ = kind, std::enable_if_t<!(K_ == Kind::X || K_ == Kind::N), bool> = true>
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
    template <auto K_ = kind, std::enable_if_t<(K_ == Kind::X || K_ == Kind::N), bool> = true>
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
    template <auto K_ = kind, std::enable_if_t<(K_ == Kind::X || K_ == Kind::N), bool> = true>
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
    template <auto K_ = kind, std::enable_if_t<!(K_ == Kind::X || K_ == Kind::N), bool> = true>
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
    template <auto K_ = kind, std::enable_if_t<!(K_ == Kind::X || K_ == Kind::N), bool> = true>
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

    /** Return pointer to buffer data. */
    __host__ __device__ auto
    data() const -> const_pointer
    {
        return _data;
    }
    /**@}*/

    /** @{ Iterators */
    // FIXME
    /**@}*/

    /** @{ mdspan interface */
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
        metautils::condition<(kind == Kind::X), stdex::mdspan<value_type, extents_type>>,
        // 1D buffer with compile-time number of elements.
        metautils::condition<(kind == Kind::N), stdex::mdspan<value_type, extents_type>>,
        // 2D buffer with run-time number of rows and columns.
        metautils::condition<(kind == Kind::XY), stdex::mdspan<value_type, extents_type, stdex::layout_stride>>,
        // 2D buffer with compile-time number of rows.
        metautils::condition<(kind == Kind::MY), stdex::mdspan<value_type, extents_type, stdex::layout_stride>>,
        // 2D buffer with compile-time number of columns.
        metautils::condition<(kind == Kind::XN), stdex::mdspan<value_type, extents_type, stdex::layout_stride>>,
        // 2D buffer with compile-time number of rows and columns.
        metautils::condition<(kind == Kind::MN), stdex::mdspan<value_type, extents_type, stdex::layout_stride>>,
        // fallback type
        stdex::mdspan<value_type, extents_type, stdex::layout_stride>>;

    using const_mdspan_type = metautils::select_t<
        // 1D buffer with run-time number of elements.
        metautils::condition<(kind == Kind::X), stdex::mdspan<const value_type, extents_type>>,
        // 1D buffer with compile-time number of elements.
        metautils::condition<(kind == Kind::N), stdex::mdspan<const value_type, extents_type>>,
        // 2D buffer with run-time number of rows and columns.
        metautils::condition<(kind == Kind::XY), stdex::mdspan<const value_type, extents_type, stdex::layout_stride>>,
        // 2D buffer with compile-time number of rows.
        metautils::condition<(kind == Kind::MY), stdex::mdspan<const value_type, extents_type, stdex::layout_stride>>,
        // 2D buffer with compile-time number of columns.
        metautils::condition<(kind == Kind::XN), stdex::mdspan<const value_type, extents_type, stdex::layout_stride>>,
        // 2D buffer with compile-time number of rows and columns.
        metautils::condition<(kind == Kind::MN), stdex::mdspan<const value_type, extents_type, stdex::layout_stride>>,
        // fallback type
        stdex::mdspan<const value_type, extents_type, stdex::layout_stride>>;

    /** @{ Implicit conversions */
    /** Conversion to mdspan. */
    // operator mdspan_type()
    //{
    //     if constexpr (kind == Kind::X || kind == Kind::N)
    //     {
    //         return mdspan_type(_data, _nColumns);
    //     }
    //     else
    //     {
    //         stdex::layout_stride::mapping layout{stdex::extents{_nRows, _nColumns}, std::array<size_t, 2>{_nPaddedColumns, 1}};
    //         return mdspan_type(_data, layout);
    //     }
    // }

    ///** Conversion to const mdspan. */
    // operator const_mdspan_type() const
    //{
    //     if constexpr (kind == Kind::X || kind == Kind::N)
    //     {
    //         return const_mdspan_type(_data, _nColumns);
    //     }
    //     else
    //     {
    //         stdex::layout_stride::mapping layout{stdex::extents{_nRows, _nColumns}, std::array<size_t, 2>{_nPaddedColumns, 1}};
    //         return const_mdspan_type(_data, layout);
    //     }
    // }
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
    /** Check approximate equality of 1-dimensional buffers.
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
    template <typename Bother, auto K_ = kind, std::enable_if_t<(K_ == Kind::X || K_ == Kind::N), bool> = true>
    friend inline auto
    operator==(const CBuffer<T, B, 1, NCols> &lhs, const CBuffer<T, Bother, 1, NCols> &rhs) -> bool
    {
        if (lhs._nColumns != rhs._nColumns)
        {
            return false;
        }

        for (auto i = 0; i < lhs._nColumns; ++i)
        {
            if (std::abs(lhs(i) - rhs(i)) > 1.0e-13)
            {
                return false;
            }
        }

        return true;
    }

    /** Check approximate inequality of 1-dimensional buffers.
     *
     * @tparam Bother backend of right-hand side in comparison.
     * @param[in] lhs left-hand side of comparison.
     * @param[in] rhs right-hand side of comparison.
     */
    template <typename Bother, auto K_ = kind, std::enable_if_t<(K_ == Kind::X || K_ == Kind::N), bool> = true>
    friend inline auto
    operator!=(const CBuffer<T, B, 1, NCols> &lhs, const CBuffer<T, Bother, 1, NCols> &rhs) -> bool
    {
        return !(lhs == rhs);
    }

    /** Check approximate equality of 2-dimensional buffers.
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
    template <typename Bother, auto K_ = kind, std::enable_if_t<!(K_ == Kind::X || K_ == Kind::N), bool> = true>
    friend inline auto
    operator==(const CBuffer<T, B, NRows, NCols> &lhs, const CBuffer<T, Bother, NRows, NCols> &rhs) -> bool
    {
        if (lhs._nRows != rhs._nRows)
        {
            return false;
        }

        if (lhs._nColumns != rhs._nColumns)
        {
            return false;
        }

        for (auto i = 0; i < lhs._nRows; ++i)
        {
            for (auto j = 0; j < lhs._nColumns; ++j)
            {
                if (std::abs(lhs(i, j) - rhs(i, j)) > 1.0e-13)
                {
                    return false;
                }
            }
        }

        return true;
    }

    /** Check approximate inequality of 2-dimensional buffers.
     *
     * @tparam Bother backend of right-hand side in comparison.
     * @param[in] lhs left-hand side of comparison.
     * @param[in] rhs right-hand side of comparison.
     */
    template <typename Bother, auto K_ = kind, std::enable_if_t<!(K_ == Kind::X || K_ == Kind::N), bool> = true>
    friend inline auto
    operator!=(const CBuffer<T, B, NRows, NCols> &lhs, const CBuffer<T, Bother, NRows, NCols> &rhs) -> bool
    {
        return !(lhs == rhs);
    }
    /**@}*/

    /** @{ Fillers/generators */
    /** Sets all elements in 1D buffer to given value.
     *
     * @param v value of fill element.
     * @param nCols number of columns. Defaults to the NCols template parameter.
     * @note The padding elements are left uninitialized!
     */
    template <auto K_ = kind, std::enable_if_t<(K_ == Kind::X || K_ == Kind::N), bool> = true>
    auto
    setConstant(value_type v, size_type nCols = NCols) -> void
    {
        // allocate
        if constexpr (kind == Kind::X)
        {
            allocate(nCols);
        }
        else
        { /* nothing to do if the buffer is of MN kind! */
        }

        // set all of _data, except padding elements, to the specified constant
        if constexpr (mem::is_on_host_v<backend_type>)
        {
            auto pdata = _data;
#pragma omp simd aligned(pdata : Alignment)
            for (size_type i = 0; i < _nColumns; ++i)
                pdata[i] = v;
        }
        else
        {
#ifdef VLX_USE_DEVICE
            const auto block = dim3(256);
            const auto grid  = dim3((_nColumns + block.x - 1) / block.x);
            deviceLaunch((device::full1D<value_type>), grid, block, 0, 0, _data, v, _nColumns);
            DEVICE_CHECK(deviceStreamSynchronize(0));
#endif
        }
    }

    /** Sets all elements in 2D buffer to given value.
     *
     * @param v value of fill element.
     * @param nRows number of rows. Defaults to the NRows template parameter.
     * @param nCols number of columns. Defaults to the NCols template parameter.
     * @note The padding elements are left uninitialized!
     */
    template <auto K_ = kind, std::enable_if_t<!(K_ == Kind::X || K_ == Kind::N), bool> = true>
    auto
    setConstant(value_type fill_value, size_type nRows = NRows, size_type nCols = NCols) -> void
    {
        // allocate
        if constexpr (kind == Kind::MY)
        {
            allocate(nCols);
        }
        else if constexpr (kind == Kind::XN)
        {
            allocate(nRows);
        }
        else if constexpr (kind == Kind::XY)
        {
            allocate(nRows, nCols);
        }
        else
        { /* nothing to do if the buffer is of MN kind! */
        }

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

    /** Zero out 1D buffer.
     *
     * @param nCols number of columns. Defaults to the NCols template parameter.
     * @note The padding elements are left uninitialized!
     */
    template <auto K_ = kind, std::enable_if_t<(K_ == Kind::X || K_ == Kind::N), bool> = true>
    auto
    setZero(size_type nCols = NCols) -> void
    {
        setConstant(T{0}, nCols);
    }

    /** Zero out 2D buffer.
     *
     * @param nRows number of rows. Defaults to the NRows template parameter.
     * @param nCols number of columns. Defaults to the NCols template parameter.
     * @note The padding elements are left uninitialized!
     */
    template <auto K_ = kind, std::enable_if_t<!(K_ == Kind::X || K_ == Kind::N), bool> = true>
    auto
    setZero(size_type nRows = NRows, size_type nCols = NCols) -> void
    {
        setConstant(T{0}, nRows, nCols);
    }

    /** Sets 1D buffer elements to uniform random values in given inteval.
     *
     * @param[in] lower lower bound of interval.
     * @param[in] upper upper bound of interval.
     * @param nCols number of columns. Defaults to the NCols template parameter.
     * @note The padding elements are left uninitialized!
     *
     * The random numbers are generated on the CPU, then copied into the buffer.
     * This method uses the C++ default random engine with random seeding.
     * If you need more control, generate the random sequence (as std::vector or
     * raw array) and then use any of the available constructors.
     */
    template <auto K_ = kind, std::enable_if_t<(K_ == Kind::X || K_ == Kind::N), bool> = true>
    auto
    setRandom(T lower, T upper, size_type nCols = NCols) -> void
    {
        // allocate
        if constexpr (kind == Kind::X)
        {
            allocate(nCols);
        }
        else
        { /* nothing to do if the buffer is of MN kind! */
        }

        // generate
        auto tmp = new T[_nElements];
        mathfunc::fill_random(tmp, lower, upper, _nElements);

        // copy into buffer
        _copy_unaligned(tmp);

        // delete temporary
        delete[] tmp;
    }

    /** Sets 2D buffer elements to uniform random values in given inteval.
     *
     * @param[in] lower lower bound of interval.
     * @param[in] upper upper bound of interval.
     * @param nRows number of rows. Defaults to the NRows template parameter.
     * @param nCols number of columns. Defaults to the NCols template parameter.
     * @note The padding elements are left uninitialized!
     *
     * The random numbers are generated on the CPU, then copied into the buffer.
     * This method uses the C++ default random engine with random seeding.
     * If you need more control, generate the random sequence (as std::vector or
     * raw array) and then use any of the available constructors.
     */
    template <auto K_ = kind, std::enable_if_t<!(K_ == Kind::X || K_ == Kind::N), bool> = true>
    auto
    setRandom(T lower, T upper, size_type nRows = NRows, size_type nCols = NCols) -> void
    {
        // allocate
        if constexpr (kind == Kind::MY)
        {
            allocate(nCols);
        }
        else if constexpr (kind == Kind::XN)
        {
            allocate(nRows);
        }
        else if constexpr (kind == Kind::XY)
        {
            allocate(nRows, nCols);
        }
        else
        { /* nothing to do if the buffer is of MN kind! */
        }

        // generate
        auto tmp = new T[_nElements];
        mathfunc::fill_random(tmp, lower, upper, _nElements);

        // copy into buffer
        _copy_unaligned(tmp);

        // delete temporary
        delete[] tmp;
    }
    /**@}*/
};
}  // namespace buffer

#endif /* BufferImpl_hpp */
