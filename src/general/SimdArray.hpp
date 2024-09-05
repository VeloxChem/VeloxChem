#ifndef SimdArray_hpp
#define SimdArray_hpp

#include <array>
#include <cstddef>
#include <cstdlib>
#include <vector>

#include "Point.hpp"
#include "omp.h"

namespace simd {  // simd namespace

/// @brief Gets IPO optimization size for SIMD vector size.
/// @return The IPO multiplier for cashline size.
constexpr auto
ipo_size() -> size_t
{
    return 2;
}

/// @brief Gets width of SIMD vector.
/// @return The size of SIMD vector.
template <typename T>
constexpr auto
width() -> size_t
{
    return 64 * ipo_size() / sizeof(T);
}

}  // namespace simd

/// @brief Class CSimdArray stores 2D array data as plain vector in form suitable for SIMD operations.
template <typename T>
class CSimdArray
{
   public:
    /// @brief The default constructor.
    /// @param rows The number of rows in array.
    /// @param columns The number of columns in SIMD widths.
    CSimdArray(const size_t rows, const size_t columns)
        : _data{nullptr}

        , _rows{rows}

        , _columns{columns}

        , _active_width{simd::width<T>()}
    {
        if (auto nelems = _rows * _columns; nelems > 0)
        {
            void *ptr = nullptr;

            ::posix_memalign(&ptr, 64, nelems * simd::width<T>() * sizeof(T));

            _data = (T*)ptr;
        }
    }

    /// @brief The default copy constructor.
    /// @param other The SIMD array to be copied.
    CSimdArray(const CSimdArray &other) = delete;

    /// @brief The default move constructor.
    /// @param other The SIMD array to be moved.
    CSimdArray(CSimdArray &&other) noexcept = delete;

    /// @brief The custom destructor.
    ~CSimdArray()
    {
        if (_data != nullptr) ::free((void*)_data);
    };

    /// @brief The default copy assignment operator.
    /// @param other The SIMD array to be copy assigned.
    /// @return The assigned SIMD array.
    auto operator=(const CSimdArray &other) -> CSimdArray & = delete;

    /// @brief The default move assignment operator.
    /// @param other The SIMD array to be move assigned.
    /// @return The assigned SIMD array.
    auto operator=(CSimdArray &&other) noexcept -> CSimdArray & = delete;

    /// @brief The equality operator.
    /// @param other The SIMD array to be compared.
    /// @return True if SIMD arrays are equal, False otherwise.
    auto operator==(const CSimdArray &other) const -> bool = delete;

    /// @brief The non-equality operator.
    /// @param other The SIMD array to be compared.
    /// @return True if SIMD arrays are not equal, False otherwise.
    auto operator!=(const CSimdArray &other) const -> bool = delete;

    /// @brief Sets active width.
    /// @param width The active width of SIMD width.
    auto
    set_active_width(const size_t width) -> void
    {
        _active_width = width;
    };

    /// @brief Gets pointer to SIMD array row.
    /// @param irow The requested row.
    /// @return The pointer to SIMD array row.
    inline auto
    data(const size_t irow) -> T *
    {
        return &_data[simd::width<T>() * _columns * irow];
    }

    /// @brief Gets constant pointer to SIMD array row.
    /// @param irow The requested row.
    /// @return The constant pointer to SIMD array row.
    inline auto
    data(const size_t irow) const -> const T *
    {
        return &_data[simd::width<T>() * _columns * irow];
    }

    /// @brief Gets number of rows in SIMD array.
    /// @return The number of rows.
    inline auto
    number_of_rows() const -> size_t
    {
        return _rows;
    }

    /// @brief Gets number of columns in SIMD array measured in SIMD blocks.
    /// @return The number of columns.
    inline auto
    number_of_columns() const -> size_t
    {
        return _columns;
    }

    /// @brief Gets number of active elements in column.
    /// @return The number of active elements.
    inline auto
    number_of_active_elements() const -> size_t
    {
        return _columns * _active_width;
    }

    /// @brief Sets all elements of SIMD array  to zero.
    inline auto
    zero() -> void
    {
        const auto nelems = _rows * _columns * simd::width<T>();

        const auto fact = T{0.0};

        auto ptr_data = _data;

#pragma omp simd aligned(ptr_data : 64)
        for (size_t i = 0; i < nelems; i++)
        {
            ptr_data[i] = fact;
        }
    }

    /// @brief Replicates data from vector into SIMD array.
    /// @param source  The vector to load data from.
    /// @param indices  The range [first, last) of data to replicate.
    /// @param position The row position of replicated data in SIMD array.
    /// @param nreps  The number of repetitions of source data.
    inline auto
    replicate(const std::vector<T> &source, std::pair<size_t, size_t> indices, const size_t position, const size_t nreps) -> void
    {
        if (const size_t ndim = indices.second - indices.first; ndim <= simd::width<T>())
        {
            _active_width = ndim;

            const size_t row_off = simd::width<T>() * _columns * position;

            for (size_t i = 0; i < nreps; i++)
            {
                for (size_t j = indices.first; j < indices.second; j++)
                {
                    _data[row_off + i * ndim + j - indices.first] = source[j];
                }
            }
        }
    }

    /// @brief Replicates data from vector into SIMD array.
    /// @param source  The vector to load data from.
    /// @param indices  The range [first, last) of data to replicate.
    /// @param position The row position of replicated data in SIMD array.
    /// @param nreps  The number of repetitions of source data.
    inline auto
    replicate_points(const std::vector<TPoint<T>> &source, std::pair<size_t, size_t> indices, const size_t position, const size_t nreps) -> void
    {
        if (const size_t ndim = indices.second - indices.first; ndim <= simd::width<T>())
        {
            _active_width = ndim;

            const size_t row_off_x = simd::width<T>() * _columns * position;

            const size_t row_off_y = row_off_x + simd::width<T>() * _columns;

            const size_t row_off_z = row_off_y + simd::width<T>() * _columns;

            for (size_t i = 0; i < nreps; i++)
            {
                for (size_t j = indices.first; j < indices.second; j++)
                {
                    const auto [x, y, z] = source[j].coordinates();

                    _data[row_off_x + i * ndim + j - indices.first] = x;

                    _data[row_off_y + i * ndim + j - indices.first] = y;

                    _data[row_off_z + i * ndim + j - indices.first] = z;
                }
            }
        }
    }

    /// @brief Load data from vector into SIMD array.
    /// @param source  The vector to load data from.
    /// @param indices  The range [first, last) of data to replicate.
    /// @param position The row position of replicated data in SIMD array.
    /// @param nblocks : the number of blocks in source vector.
    inline auto
    load(const std::vector<T> &source, std::pair<size_t, size_t> indices, const size_t position, size_t nblocks) -> void
    {
        if (const size_t ndim = indices.second - indices.first; ndim <= simd::width<T>())
        {
            _active_width = ndim;

            const size_t row_off = simd::width<T>() * _columns * position;

            const size_t sdim = source.size() / nblocks;

            for (size_t i = 0; i < nblocks; i++)
            {
                for (size_t j = indices.first; j < indices.second; j++)
                {
                    _data[row_off + i * ndim + j - indices.first] = source[i * sdim + j];
                }
            }
        }
    }

    /// @brief Load data from vector of Cartesian points into SIMD array.
    /// @param source  The vector of Cartesian points to load data from.
    /// @param indices  The range [first, last) of data to replicate.
    /// @param position The row position of replicated data in SIMD array.
    /// @param nblocks : the number of blocks in source vector.
    inline auto
    load_points(const std::vector<TPoint<T>> &source, std::pair<size_t, size_t> indices, const size_t position, size_t nblocks) -> void
    {
        if (const size_t ndim = indices.second - indices.first; ndim <= simd::width<T>())
        {
            _active_width = ndim;

            const size_t row_off_x = simd::width<T>() * _columns * position;

            const size_t row_off_y = row_off_x + simd::width<T>() * _columns;

            const size_t row_off_z = row_off_y + simd::width<T>() * _columns;

            const size_t sdim = source.size() / nblocks;

            for (size_t i = 0; i < nblocks; i++)
            {
                for (size_t j = indices.first; j < indices.second; j++)
                {
                    const auto [x, y, z] = source[i * sdim + j];

                    _data[row_off_x + i * ndim + j - indices.first] = x;

                    _data[row_off_y + i * ndim + j - indices.first] = y;

                    _data[row_off_z + i * ndim + j - indices.first] = z;
                }
            }
        }
    }

   private:
    /// @brief Memory block for data storage.
    T *_data;

    /// @brief Number of rows in 2D array.
    size_t _rows;

    /// @brief Number of columns in 2D array in SIMD blocks.
    size_t _columns;

    /// @brief Number of avtive elements in default SIMD block.
    size_t _active_width;
};

#endif /* SimdArray_hpp */
