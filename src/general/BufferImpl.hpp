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

#include <mpi.h>

#include <array>
#include <cstddef>
#include <cstdint>
#include <sstream>
#include <string>
#include <type_traits>
#include <vector>

#include "Device.hpp"
#include "ErrorHandler.hpp"
#include "IntegerSequence.hpp"
#include "MathFunc.hpp"
#include "MemAlloc.hpp"
#include "MetaUtils.hpp"
#include "MpiFunc.hpp"
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

    static constexpr bool Layout1D           = (NRows == 1);
    static constexpr bool CompileTimeRows    = (NRows != Dynamic);
    static constexpr bool CompileTimeColumns = (NCols != Dynamic);
    /** Whether the buffer has both number of rows and columns known at compile-time.
     *
     * @note This is `true` for `Kind::N` and `Kind::MN` buffers.
     */
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

    /** Generate a linearly-spaced 1-dimensional buffer within interval.
     *
     * @param[in] start lower bound of interval.
     * @param[in] end upper bound of interval.
     * @note the data is generated in [start, end)
     */
    template <auto L_ = Layout1D, std::enable_if_t<L_, bool> = true>
    auto
    _linspace(value_type start, value_type end) -> void
    {
        auto delta = (end - start) / static_cast<value_type>(_nColumns);

        if constexpr (mem::is_on_host_v<backend_type>)
        {
            auto pdata = _data;
#pragma omp simd aligned(pdata : Alignment)
            for (size_type i = 0; i < _nColumns; ++i)
                pdata[i] = start + delta * i;
        }
        else
        {
#ifdef VLX_USE_DEVICE
            const auto block = dim3(256);
            const auto grid  = dim3((_nColumns + block.x - 1) / block.x);
            // TODO write kernel!
            // deviceLaunch((device::full1D<value_type>), grid, block, 0, 0, _data, fill_value, _nColumns);
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
            auto pdata = _data;
            for (size_type i = 0; i < _nRows; ++i)
            {
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

    /** Generate a linearly-spaced 2-dimensional buffer within interval.
     *
     * @param[in] start lower bound of interval.
     * @param[in] end upper bound of interval.
     * @note the data is generated in [start, end)
     */
    template <auto L_ = Layout1D, std::enable_if_t<!L_, bool> = true>
    auto
    _linspace(value_type start, value_type end) -> void
    {
        auto delta = (end - start) / static_cast<value_type>(_nRows * _nColumns);

        if constexpr (mem::is_on_host_v<backend_type>)
        {
            auto pdata = _data;
            for (size_type i = 0; i < _nRows; ++i)
            {
#pragma omp simd aligned(pdata : Alignment)
                for (size_type j = 0; j < _nColumns; ++j)
                {
                    pdata[i * _nPaddedColumns + j] = start + delta * (i * _nColumns + j);
                }
            }
        }
        else
        {
#ifdef VLX_USE_DEVICE
            const auto block = dim3(16, 16);
            const auto grid  = dim3((_nRows + block.x - 1) / block.x, (_nColumns + block.y - 1) / block.y);
            // TODO
            // deviceLaunch((device::full2D<value_type>), grid, block, 0, 0, _data, fill_value, _nRows, _nColumns, _nPaddedColumns);
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

    /** Set buffer dimensions.
     *
     * @tparam Extents variadic pack of extents dimensions.
     *
     * @param[in] extents dimensions of the extents.
     */
    template <typename... Extents>
    auto
    _setDimensions(Extents... extents) -> void
    {
        static_assert(std::conjunction_v<std::is_integral<Extents>...>, "Extents must be of integral type.");

        auto tup = std::tuple<Extents...>(static_cast<size_type>(extents)...);

        constexpr auto tup_sz = std::tuple_size_v<decltype(tup)>;

        static_assert(tup_sz <= 2, "Dimensioning tuple must be at most of size 2!");

        if constexpr (kind == Kind::XY)
        {
            static_assert(tup_sz == 2, "For Kind::XY buffer the dimensioning tuple must be exactly of size 2!");
            std::tie(_nRows, _nColumns) = tup;
        }
        else if constexpr (kind == Kind::X)
        {
            static_assert(tup_sz > 0, "For Kind::X buffer the dimensioning tuple must be >0!");

            if constexpr (tup_sz > 1)
            {
                std::tie(std::ignore, _nColumns) = tup;
            }
            else
            {
                std::tie(_nColumns) = tup;
            }
        }
        else if constexpr (kind == Kind::XN)
        {
            static_assert(tup_sz > 0, "For Kind::XN buffer the dimensioning tuple must be >0!");

            if constexpr (tup_sz > 1)
            {
                std::tie(_nRows, std::ignore) = tup;
            }
            else
            {
                std::tie(_nRows) = tup;
            }
        }
        else if constexpr (kind == Kind::MY)
        {
            static_assert(tup_sz > 0, "For Kind::MY buffer the dimensioning tuple must be >0!");

            if constexpr (tup_sz > 1)
            {
                std::tie(std::ignore, _nColumns) = tup;
            }
            else
            {
                std::tie(_nColumns) = tup;
            }
        }
        else
        {
            // for Kind::N and Kind::MN there are no resizable extents: ignore whatever was passed in!
        }
    }

    /** Allocation function.
     *
     * @tparam Extents variadic pack of extents dimensions.
     *
     * @param[in] extents dimensions of the extents.
     *
     * @note This function performs a re-allocation if `_data` already points to
     * some memory, *i.e.* is not `nullptr`.
     */
    template <typename... Extents>
    auto
    _allocate(Extents... extents) -> void
    {
        static_assert(std::conjunction_v<std::is_integral<Extents>...>, "Extents must be of integral type.");

        if (_data)
        {
            mem::free<value_type, backend_type>(_data);
        }

        // unpack variadic parameter and set _nRows and _nColumns
        _setDimensions(extents...);

        if constexpr (kind == Kind::X)
        {
            _nPaddedColumns = mem::get_pitch<value_type>(Alignment, _nColumns);
            _nElements      = _nRows * _nPaddedColumns;
            _data           = mem::malloc<value_type, backend_type>(_nElements);
        }
        else if constexpr (kind == Kind::XY)
        {
            std::tie(_nPaddedColumns, _data) = mem::malloc<value_type, backend_type>(_nRows, _nColumns);
            _nElements                       = _nRows * _nPaddedColumns;
        }
        else if constexpr (kind == Kind::XN)
        {
            std::tie(_nPaddedColumns, _data) = mem::malloc<value_type, backend_type>(_nRows, _nColumns);
            _nElements                       = _nRows * _nPaddedColumns;
        }
        else if constexpr (kind == Kind::MY)
        {
            std::tie(_nPaddedColumns, _data) = mem::malloc<value_type, backend_type>(_nRows, _nColumns);
            _nElements                       = _nRows * _nPaddedColumns;
        }
        else
        {
            _nElements = _nRows * _nPaddedColumns;
            _data      = mem::malloc<value_type, backend_type>(_nElements);
        }
    }

    /** String representation of this bufffer's type. */
    auto
    _type_to_string() const -> std::string
    {
        auto dim_str = [](size_type dim) {
            if (dim == Dynamic)
            {
                return std::string{"Dynamic"};
            }
            else
            {
                return std::to_string(NRows);
            }
        };

        std::ostringstream os;

        os << "[CBuffer<" << metautils::type_to_string<value_type>::name << ", " << backend_type::name << ", " << dim_str(NRows) << ", "
           << dim_str(NCols) << "> (Object): " << this << "]";

        return os.str();
    }

    /** @{ Slicing functions */
    /** Obtain a new 1D buffer as slice of current 1D buffer.
     *
     * @tparam Indices the type of the array of indices describing the slice.
     * @param[in] idxs the indices describing the slice.
     * @return a 1D slice of the current 1D buffer.
     *
     * The `idxs` argument can be any object that implements:
     *
     * \code
     * <integral type> operator[](<integral type>) const;
     * <integral type> size() const;
     * \endcode
     *
     * for example, `std::array`, `std::vector`, and `buffer::IndexSequence`.
     */
    template <typename Indices,
              typename ReturnType        = CBuffer<value_type, backend_type, size_type{1}, Dynamic>,
              auto L_                    = Layout1D,
              std::enable_if_t<L_, bool> = true>
    [[nodiscard]] auto
    _slice(const Indices &idxs) const -> ReturnType
    {
        // allocate new buffer of appropriate dimensions
        auto buf = ReturnType(idxs.size());

        if constexpr (mem::is_on_host_v<backend_type>)
        {
            // loop over given row and column indices
            // and copy elements from `this` buffer to destination
            auto src = _data;
            auto dst = buf.data();
#pragma omp simd aligned(src, dst : Alignment)
            for (auto i = 0; i < idxs.size(); ++i)
            {
                dst[i] = src[idxs[i]];
            }
        }
        else
        {
            // TODO GPU implementation
        }

        return buf;
    }

    /** Obtain a new 2D buffer as slice of current 2D buffer.
     *
     * @tparam RowIndices the type of the array of indices describing the rows in the slice.
     * @tparam ColIndices the type of the array of indices describing the columns in the slice.
     * @param[in] rows the indices describing the rows in the slice.
     * @param[in] cols the indices describing the columns in the slice.
     * @return a 1D slice of the current 1D buffer.
     *
     * The `rows` and `cols` arguments can be any object that implements:
     *
     * \code
     * <integral type> operator[](<integral type>) const;
     * <integral type> size() const;
     * \endcode
     *
     * for example, `std::array`, `std::vector`, and `buffer::IndexSequence`.
     */
    template <typename RowIndices,
              typename ColIndices,
              typename ReturnType         = CBuffer<value_type, backend_type, Dynamic, Dynamic>,
              auto L_                     = Layout1D,
              std::enable_if_t<!L_, bool> = true>
    [[nodiscard]] auto
    _slice(const RowIndices &rows, const ColIndices &cols) const -> ReturnType
    {
        // allocate new buffer of appropriate dimensions
        auto buf = ReturnType(rows.size(), cols.size());

        if constexpr (mem::is_on_host_v<backend_type>)
        {
            // loop over given row and column indices
            // and copy elements from `this` buffer to destination
            auto src           = _data;
            auto dst           = buf.data();
            auto n_padded_cols = buf.nPaddedColumns();

            for (auto i = 0; i < rows.size(); ++i)
            {
#pragma omp simd aligned(src, dst : Alignment)
                for (auto j = 0; j < cols.size(); ++j)
                {
                    dst[i * n_padded_cols + j] = src[rows[i] * _nPaddedColumns + cols[j]];
                }
            }
        }
        else
        {
            // TODO GPU implementation
        }

        return buf;
    }
    /**}@*/

    /** @{ MPI functions */
    /** Compute pattern for MPI scattering/gathering of 1D buffer.
     *
     * @param[in] rank the rank of the calling process.
     * @param[in] nodes number of nodes to scatter to.
     * @return triple of types `int32_t`, `int32_t*` and `int32_t*`.
     * - In a **scattering operation**, these correspond to `recvcount`,
     *   `sendcounts`, and `displs`, respectively, and describe the buffer on
     *   the receiving and sending ends of the scatter operation.
     * - In a **gathering operation**, these correspond to `sendcount`,
     *   `recvcounts`, and `displs`, respectively, and describe the buffer on
     *   the sending and receiving ends of the gather operation.
     *
     * @note For 1D arrays, we partition the data in the row, excluding the
     * padding elements.
     *
     * For example, a 1D buffer holding 9 `double`s aligned to 64-byte boundary,
     * will be padded with 7 additional `double`s:
     *
     * [ o o o o o o o o | o x x x x x x x ]
     *
     * When scattering over 4 processes, each rank will "allocate" a similar
     * buffer of 9 `double`s aligned to 64-byte boundary. This means 9 columns
     * plus 7 padding elements.
     * After scattering, the ranks will contain 3, 2, 2, and 2 elements,
     * respectively, of the original buffer, **at their head**:
     *
     * rank 0: [ o o o * * * * * | x x x x x x x x ]
     *
     * rank 1: [ o o * * * * * * | x x x x x x x x ]
     *
     * rank 2: [ o o * * * * * * | x x x x x x x x ]
     *
     * rank 3: [ o o * * * * * * | x x x x x x x x ]
     *
     * Here, the "x" are padding elements, which are inaccessible when accessing
     * elements of the buffer through `operator()` and `at()`.
     * The "*" are placeholder elements: they are rubbish, but still accessible!
     *
     * When gathering over 4 processes, we pick 3, 2, 2, and 2 elements,
     * respectively, from the **head** of the buffer on each process.
     * So starting from:
     *
     * rank 0: [ o o o * * * * * | x x x x x x x x ]
     *
     * rank 1: [ o o * * * * * * | x x x x x x x x ]
     *
     * rank 2: [ o o * * * * * * | x x x x x x x x ]
     *
     * rank 3: [ o o * * * * * * | x x x x x x x x ]
     *
     * we get:
     *
     * rank 0: [ o o o o o o o o | o x x x x x x x ]
     */
    template <auto L_ = Layout1D, std::enable_if_t<L_, bool> = true>
    [[nodiscard]] auto
    _pattern_scatter_gather(int32_t rank, int32_t nodes) -> std::tuple<int32_t, int32_t *, int32_t *>
    {
        // how many elements should "rank" receive from/send to the root?
        auto count = mpi::batch_size(_nColumns, rank, nodes);

        // how many elements should the root send to/receive from "rank"?
        int32_t *counts = (rank == mpi::master()) ? mem::malloc<int32_t>(nodes) : nullptr;

        // where to read the elements sent to/received from "rank"?
        int32_t *displs = (rank == mpi::master()) ? mem::malloc<int32_t>(nodes) : nullptr;

        // fill counts and displs
        if (rank == mpi::master())
        {
            mpi::batches_pattern(counts, _nColumns, nodes);

            mathfunc::indexes(displs, counts, nodes);
        }

        return {count, counts, displs};
    }

    /** Compute pattern for MPI scattering/gathering of 2D buffer.
     *
     * @param[in] rank the rank of the calling process.
     * @param[in] nodes number of nodes to scatter to.
     * @return triple of types `int32_t`, `int32_t*` and `int32_t*`.
     * - In a **scattering operation**, these correspond to `recvcount`,
     *   `sendcounts`, and `displs`, respectively, and describe the buffer on
     *   the receiving and sending ends of the scatter operation.
     * - In a **gathering operation**, these correspond to `sendcount`,
     *   `recvcounts`, and `displs`, respectively, and describe the buffer on
     *   the sending and receiving ends of the gather operation.
     *
     * @note For 2D arrays, we partition the data by rows.
     *
     * For example, a 2D buffer with 3 rows and 9 columns holds 27 `double`s
     * aligned to 64-byte boundary.  Each row will be padded with 7 additional
     * `double`s, for an allocation of 48 `double`-s:
     *
     * / o o o o o o o o | o x x x x x x x \
     * | o o o o o o o o | o x x x x x x x |
     * \ o o o o o o o o | o x x x x x x x /
     *
     * When scattering over 4 processes, each rank will "allocate" a similar
     * buffer of 3x9 `double`s aligned to 64-byte boundary. This means 3 rows
     * each with 9 columns plus 7 padding elements.
     * After scattering, the ranks will contain (or not) complete rows of the
     * original buffer as **first row**:
     *
     * rank 0: / o o o o o o o o | o x x x x x x x \
     *         | * * * * * * * * | * x x x x x x x |
     *         \ * * * * * * * * | * x x x x x x x /
     *
     * rank 1: / o o o o o o o o | o x x x x x x x \
     *         | * * * * * * * * | * x x x x x x x |
     *         \ * * * * * * * * | * x x x x x x x /
     *
     * rank 2: / o o o o o o o o | o x x x x x x x \
     *         | * * * * * * * * | * x x x x x x x |
     *         \ * * * * * * * * | * x x x x x x x /
     *
     * rank 3: / * * * * * * * * | * x x x x x x x \
     *         | * * * * * * * * | * x x x x x x x |
     *         \ * * * * * * * * | * x x x x x x x /
     *
     * Here, the "x" are padding elements, which are inaccessible when accessing
     * elements of the buffer through `operator()` and `at()`.
     * The "*" are placeholder elements: they are rubbish, but still accessible!
     * Note how rank 3 (the 4th process) gets no actual values after the scatter!
     *
     * When gathering over 4 processes, we pick the first 1, 1, 1, and 0 rows,
     * respectively, from the buffer on each rank.
     * So starting from:
     *
     * rank 0: / o o o o o o o o | o x x x x x x x \
     *         | * * * * * * * * | * x x x x x x x |
     *         \ * * * * * * * * | * x x x x x x x /
     *
     * rank 1: / o o o o o o o o | o x x x x x x x \
     *         | * * * * * * * * | * x x x x x x x |
     *         \ * * * * * * * * | * x x x x x x x /
     *
     * rank 2: / o o o o o o o o | o x x x x x x x \
     *         | * * * * * * * * | * x x x x x x x |
     *         \ * * * * * * * * | * x x x x x x x /
     *
     * rank 3: / * * * * * * * * | * x x x x x x x \
     *         | * * * * * * * * | * x x x x x x x |
     *         \ * * * * * * * * | * x x x x x x x /
     *
     * we get:
     *
     * rank 0: / o o o o o o o o | o x x x x x x x \
     *         | o o o o o o o o | o x x x x x x x |
     *         \ o o o o o o o o | o x x x x x x x /
     */
    template <auto L_ = Layout1D, std::enable_if_t<!L_, bool> = true>
    [[nodiscard]] auto
    _pattern_scatter_gather(int32_t rank, int32_t nodes) -> std::tuple<int32_t, int32_t *, int32_t *>
    {
        // how many rows should "rank" receive from/send to the root?
        int32_t *rows_in_batch = mem::malloc<int32_t>(nodes);
        auto [quot, rem]       = std::ldiv(_nRows, nodes);
        for (auto i = 0; i < nodes; ++i)
        {
            rows_in_batch[i] = ((rem) && (i < rem)) ? quot + 1 : quot;
        }

        // how many elements should "rank" receive from/send to the root?
        // compute as:
        // \f[
        //   R_{i} * N + (R_{i} - 1) * p
        // \f]
        // where: \f$R_{i}\f$ is the number of rows in the rank \f$i\f$ batch, \f$N\f$ the number of columns, and \f$p\f$ the number of padding
        // elements.
        // Set to zero if there are no rows to send to the given rank!
        auto count = (rows_in_batch[rank] > 0) ? (rows_in_batch[rank] - 1) * _nPaddedColumns + _nColumns : 0;

        // how many elements should the root send to/receive from "rank"?
        int32_t *counts = (rank == mpi::master()) ? mem::malloc<int32_t>(nodes) : nullptr;

        // where to read the elements sent to/received from "rank"?
        int32_t *displs = (rank == mpi::master()) ? mem::malloc<int32_t>(nodes) : nullptr;

        // fill counts and displs, taking the padding into account
        if (rank == mpi::master())
        {
            for (auto i = 0; i < nodes; ++i)
            {
                counts[i] = (rows_in_batch[i] > 0) ? (rows_in_batch[i] - 1) * _nPaddedColumns + _nColumns : 0;
            }

            // the displs for rank i accumulates the counts for rank i-1 plus
            // the padding to the displacement for rank-1
            displs[0] = 0;
            for (auto i = 1; i < nodes; ++i)
            {
                displs[i] = displs[i - 1] + rows_in_batch[i - 1] * _nPaddedColumns;
            }
        }

        mem::free(rows_in_batch);

        return {count, counts, displs};
    }
    /**}@*/

   public:
    /** Default CTOR.
     *
     * @warning This CTOR **only** allocates memory when *all* the extents of the
     * buffer are known at compile-time.
     * In all other cases, this CTOR is **non-allocating**. Thus you will have
     * to:
     *
     * - Allocate with either `resize` or `reserve`. The former makes sense
     *   when not all extents are known at compile-time.
     * - Allocate *and* initialize with one of
     *   `setConstant`/`setZero`/`setRandom`/`setLinspace`.
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
     * @note This CTOR **allocates** memory, but leaves it uninitialized.
     * You can call one of `setConstant`/`setZero`/`setRandom` to do so.
     */
    template <typename... Extents,
              auto CT_                                                                 = CompileTimeRowsAndColumns,
              std::enable_if_t<!CT_, bool>                                             = true,
              std::enable_if_t<std::conjunction_v<std::is_integral<Extents>...>, bool> = true>
    explicit CBuffer(Extents... extents)
    {
        static_assert(std::conjunction_v<std::is_integral<Extents>...>, "Extents must be of integral type.");

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

        if constexpr ((sizeof...(extents) == NResizableExtents))
        {
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
     * @tparam B_RHS backend of right-hand side in comparison.
     * @tparam NRows_RHS compile-time number of rows of right-hand side in comparison.
     * @tparam NCols_RHS compile-time number of columns of right-hand side in comparison.
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
    template <typename B_RHS, auto NRows_RHS, auto NCols_RHS>
    friend auto
    operator==(const CBuffer<T, B, NRows, NCols> &lhs, const CBuffer<T, B_RHS, NRows_RHS, NCols_RHS> &rhs) -> bool
    {
        if (lhs.size() != rhs.size())
        {
            return false;
        }

        if (lhs._nRows != rhs.nRows())
        {
            return false;
        }

        if (lhs._nColumns != rhs.nColumns())
        {
            return false;
        }

        auto rhs_padded_cols = rhs.nPaddedColumns();
        auto rhs_data        = rhs.data();
        for (size_type i = 0; i < lhs._nRows; ++i)
        {
            for (size_type j = 0; j < lhs._nColumns; ++j)
            {
                auto lhs_v = lhs._data[i * lhs._nPaddedColumns + j];
                auto rhs_v = rhs_data[i * rhs_padded_cols + j];

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
    [[nodiscard]] static auto
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
    [[nodiscard]] static auto
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
    [[nodiscard]] static auto
    Random(value_type lower, value_type upper, Extents... extents) -> CBuffer<T, B, NRows, NCols>
    {
        auto buf = CBuffer<T, B, NRows, NCols>{extents...};
        buf.setRandom(lower, upper);
        return buf;
    }

    /** Set all elements in buffer to a given value.
     *
     * @tparam Extents variadic pack of extents dimensions.
     *
     * @param[in] start lower bound of interval.
     * @param[in] end upper bound of interval.
     * @param[in] extents new dimension of the extents.
     *
     * @note The data is generated in [start, end)
     */
    template <typename... Extents>
    auto
    setLinspace(value_type start, value_type end, Extents... extents) -> void
    {
        resize(extents...);
        _linspace(start, end);
    }

    /** Creates a buffer with all elements equal to a given value.
     *
     * @tparam Extents variadic pack of extents dimensions.
     *
     * @param[in] start lower bound of interval.
     * @param[in] end upper bound of interval.
     * @param[in] extents new dimension of the extents.
     *
     * @note The data is generated in [start, end)
     */
    template <typename... Extents>
    [[nodiscard]] static auto
    Linspace(value_type start, value_type end, Extents... extents) -> CBuffer<T, B, NRows, NCols>
    {
        auto buf = CBuffer<T, B, NRows, NCols>{extents...};
        buf.setLinspace(start, end);
        return buf;
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

    /** @{ Slicing functions */
    /** Obtain a new buffer as slice of current buffer.
     *
     * @tparam SliceIndices the type(s) of the array(s) of indices describing the slice.
     * @param[in] idxs the array(s) of indices describing the slice.
     * @return a slice of the current buffer.
     *
     * The `idxs` arguments can be any object that implements:
     *
     * \code
     * <integral type> operator[](<integral type>) const;
     * <integral type> size() const;
     * \endcode
     *
     * for example, `std::array`, `std::vector`, and `buffer::IndexSequence`.
     */
    template <typename... SliceIndices, auto _NRows = (sizeof...(SliceIndices) == 1) ? 1 : Dynamic>
    [[nodiscard]] auto
    slice(SliceIndices &&...idxs) -> CBuffer<value_type, backend_type, _NRows, Dynamic>
    {
        static_assert((sizeof...(SliceIndices) == 1) || (sizeof...(SliceIndices) == 2), "");

        return _slice(std::forward<SliceIndices>(idxs)...);
    }
    /**}@*/

    /** Converts buffer object to text output. */
    [[nodiscard]] auto
    repr() const -> std::string
    {
        std::ostringstream os;

        os << std::endl;

        os << _type_to_string() << std::endl;

        os << "_nRows: " << _nRows << std::endl;

        os << "_nColumns: " << _nColumns << std::endl;

        os << "_nPaddedColumns: " << _nPaddedColumns << std::endl;

        os << "_nElements: " << _nElements << std::endl;

        os << "_data (" << &(_data) << "):" << std::endl;

        os << "[";
        auto first_row = true;
        for (size_type i = 0; i < _nRows; ++i)
        {
            if (!first_row) os << "," << std::endl;
            os << "[";
            auto first_element = true;
            for (size_type j = 0; j < _nColumns; ++j)
            {
                if (!first_element) os << ", ";
                os << _data[i * _nPaddedColumns + j];
                first_element = false;
            }
            os << "]";
            first_row = false;
        }
        os << "]";

        os << std::endl;

        return os.str();
    }

    /** @{ MPI functions */
    /** Broadcasts buffer within domain of MPI communicator.
     *
     * @param[in] rank the rank of MPI process.
     * @param[in] comm the MPI communicator.
     */
    template <typename B_ = B, typename = std::enable_if_t<!mem::is_on_device_v<B_>>>
    auto
    broadcast(int32_t rank, MPI_Comm comm) -> void
    {
        if constexpr (ENABLE_MPI)
        {
            // broadcast _nRows
            mpi::bcast(_nRows, comm);

            // broadcast _nColumns
            mpi::bcast(_nColumns, comm);

            // allocate _data array
            if (rank != mpi::master()) _allocate(_nRows, _nColumns);

            // wait for allocation to finish on all processes
            MPI_Barrier(comm);

            // broadcast _data
            // clang-format off
            auto merror = MPI_Bcast(_data
                                  , _nElements
                                  , mpi::type_v<T>
                                  , mpi::master()
                                  , comm);
            // clang-format on

            if (merror != MPI_SUCCESS) mpi::abort(merror, "broadcast(" + _type_to_string() + ")");
        }
    }

    /** Reasigns buffer on all MPI process within MPI communicator by scattering buffer from master MPI process.
     *
     * @param[in] rank the rank of MPI process.
     * @param[in] nodes the number of MPI processes in MPI communicator.
     * @param[in] comm the MPI communicator.
     *
     * Each receiving buffer will have the same dimensionality.
     */
    template <typename B_ = B, typename = std::enable_if_t<!mem::is_on_device_v<B_>>>
    auto
    scatter(int32_t rank, int32_t nodes, MPI_Comm comm) -> void
    {
        if constexpr (ENABLE_MPI)
        {
            if (nodes == 1) return;

            // set up scattering pattern
            auto [recvcount, sendcounts, displs] = _pattern_scatter_gather(rank, nodes);

            // allocate data chunk on the receiving side of the scatter operation
            auto recvbuf = mem::malloc<value_type>(_nElements);

            // scatter data chunks from master node
            // clang-format off
            auto merror = MPI_Scatterv(_data
                                     , sendcounts
                                     , displs
                                     , mpi::type_v<value_type>
                                     , recvbuf
                                     , recvcount
                                     , mpi::type_v<value_type>
                                     , mpi::master()
                                     , comm);
            // clang-format on

            if (merror != MPI_SUCCESS) mpi::abort(merror, "scatter(" + _type_to_string() + ")");

            // deallocate scattering pattern data
            mem::free(sendcounts);

            mem::free(displs);

            // swap data chunk pointers
            if (_data) mem::free(_data);

            _data = recvbuf;
        }
    }

    /** Creates buffer on master MPI process by gathering buffers from all MPI processes within MPI communicator.
     *
     * @param[in] rank the rank of MPI process.
     * @param[in] nodes the number of MPI processes in MPI communicator.
     * @param[in] comm the MPI communicator.
     * @return the buffer object:
     *   - on master node with gathered data;
     *   - on worker nodes empty.
     */
    template <typename B_ = B, typename = std::enable_if_t<!mem::is_on_device_v<B_>>>
    auto
    gather(int32_t rank, int32_t nodes, MPI_Comm comm) -> CBuffer
    {
        if constexpr (ENABLE_MPI)
        {
            if (nodes == 1) return CBuffer(*this);

            // set up gathering pattern
            auto [sendcount, recvcounts, displs] = _pattern_scatter_gather(rank, nodes);

            // declare buffer for gathered data
            CBuffer<value_type, backend_type, NRows, NCols> buf;

            if (rank == mpi::master())
            {
                // allocate buffer for gathering
                buf._allocate(_nRows, _nColumns);
            }

            // gather data chunks to master node
            // clang-format off
            auto merror = MPI_Gatherv(_data
                                    , sendcount
                                    , mpi::type_v<value_type>
                                    , buf._data
                                    , recvcounts
                                    , displs
                                    , mpi::type_v<value_type>
                                    , mpi::master()
                                    , comm);
            // clang-format on

            if (merror != MPI_SUCCESS) mpi::abort(merror, "gather(" + _type_to_string() + ")");

            // deallocate gathering pattern data
            mem::free(recvcounts);

            mem::free(displs);

            return buf;
        }

        return CBuffer(*this);
    }

    /** Sum-Reduce buffers from all MPI process within MPI communicator into buffer on master node.
     *
     * @param[in] rank the rank of MPI process.
     * @param[in] nodes the number of MPI processes in MPI communicator.
     * @param[in] comm the MPI communicator.
     */
    template <typename B_ = B, typename = std::enable_if_t<!mem::is_on_device_v<B_>>>
    auto
    reduce_sum(int32_t rank, int32_t nodes, MPI_Comm comm) -> void
    {
        if constexpr (ENABLE_MPI)
        {
            if (nodes == 1) return;

            value_type *bdata = (rank == mpi::master()) ? mem::malloc<value_type>(_nElements) : nullptr;

            mpi::reduce_sum(_data, bdata, _nElements, comm);

            if (rank == mpi::master())
            {
                mem::free(_data);

                _data = bdata;
            }
        }
    }
    /**}@*/
};
}  // namespace buffer

template <typename T, typename B, auto NRows = Dynamic, auto NCols = Dynamic>
std::ostream &
operator<<(std::ostream &output, const buffer::CBuffer<T, B, NRows, NCols> &source)
{
    return (output << source.repr());
}

#endif /* BufferImpl_hpp */
