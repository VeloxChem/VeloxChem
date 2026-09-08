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


#ifndef SimdMatrix_hpp
#define SimdMatrix_hpp

#include <algorithm>
#include <cstddef>
#include <new>
#include <string>
#include <utility>

#include "ErrorHandler.hpp"
#include "SimdAlign.hpp"

/// @brief Class CSimdMatrix stores a two-dimensional array of values in a form
/// suitable for SIMD operations. The rows are padded, so that every row starts
/// at a cache line boundary and can be loaded with aligned SIMD instructions.
class CSimdMatrix
{
   public:
    /// @brief The default constructor.
    CSimdMatrix()

        : _data(nullptr)

        , _rows(0)

        , _columns(0)

        , _pitch(0)

        , _capacity(0)

        , _owned(true)
    {
    }

    /// @brief The constructor with number of rows and columns.
    /// @param rows The number of rows in matrix.
    /// @param columns The number of columns in matrix.
    CSimdMatrix(const size_t rows, const size_t columns)

        : _data(nullptr)

        , _rows(rows)

        , _columns(columns)

        , _pitch(simd::pitch_of(columns))

        , _capacity(rows * simd::pitch_of(columns))

        , _owned(true)
    {
        _allocate();
    }

    /// @brief The constructor with borrowed values.
    /// @param values The values the matrix takes the shape over, which it does not
    /// own and does not free.
    /// @param capacity The number of values at values, which the shape of the matrix
    /// may not exceed.
    /// @note The matrix has no shape until it is reshaped. A driver which forms one
    /// arena for a block of atom pairs and reshapes a view of it for every
    /// combination of basis functions keeps the rows of a combination as close
    /// together as an owned matrix would, while asking the allocator once.
    CSimdMatrix(double *values, const size_t capacity)

        : _data(values)

        , _rows(0)

        , _columns(0)

        , _pitch(0)

        , _capacity(capacity)

        , _owned(false)
    {
    }

    /// @brief The copy constructor.
    /// @param other The matrix to be copied.
    CSimdMatrix(const CSimdMatrix &other)

        : _data(nullptr)

        , _rows(other._rows)

        , _columns(other._columns)

        , _pitch(other._pitch)

        , _capacity(other._rows * other._pitch)

        , _owned(true)
    {
        _allocate();

        if (_data != nullptr) std::copy(other._data, other._data + number_of_elements(), _data);
    }

    /// @brief The move constructor.
    /// @param other The matrix to be moved.
    CSimdMatrix(CSimdMatrix &&other) noexcept

        : _data(other._data)

        , _rows(other._rows)

        , _columns(other._columns)

        , _pitch(other._pitch)

        , _capacity(other._capacity)

        , _owned(other._owned)
    {
        other._data = nullptr;

        other._rows = 0;

        other._columns = 0;

        other._pitch = 0;

        other._capacity = 0;

        other._owned = true;
    }

    /// @brief The destructor.
    ~CSimdMatrix()
    {
        _deallocate();
    }

    /// @brief The copy assignment operator.
    /// @param other The matrix to be copy assigned.
    /// @return The assigned matrix.
    auto
    operator=(const CSimdMatrix &other) -> CSimdMatrix &
    {
        if (this != &other)
        {
            _deallocate();

            _rows = other._rows;

            _columns = other._columns;

            _pitch = other._pitch;

            _capacity = other._rows * other._pitch;

            _owned = true;

            _allocate();

            if (_data != nullptr) std::copy(other._data, other._data + number_of_elements(), _data);
        }

        return *this;
    }

    /// @brief The move assignment operator.
    /// @param other The matrix to be move assigned.
    /// @return The assigned matrix.
    auto
    operator=(CSimdMatrix &&other) noexcept -> CSimdMatrix &
    {
        if (this != &other)
        {
            _deallocate();

            _data = other._data;

            _rows = other._rows;

            _columns = other._columns;

            _pitch = other._pitch;

            _capacity = other._capacity;

            _owned = other._owned;

            other._data = nullptr;

            other._rows = 0;

            other._columns = 0;

            other._pitch = 0;

            other._capacity = 0;

            other._owned = true;
        }

        return *this;
    }

    /// @brief Sets the shape of matrix over the values it borrows.
    /// @param rows The number of rows in matrix.
    /// @param columns The number of columns in matrix.
    /// @note Only a matrix which borrows its values may be reshaped, as an owned one
    /// would have to be reallocated. The values are left as they are.
    auto
    reshape(const size_t rows, const size_t columns) -> void
    {
        if (_owned) errors::assertMsgCritical(false, std::string("SimdMatrix.reshape: Matrix owns its values"));

        const auto pitch = simd::pitch_of(columns);

        if (rows * pitch > _capacity)
        {
            errors::assertMsgCritical(false, std::string("SimdMatrix.reshape: Shape exceeds capacity of borrowed values"));
        }

        _rows = rows;

        _columns = columns;

        _pitch = pitch;
    }

    /// @brief Gets number of values matrix may hold.
    /// @return The number of values matrix owns or borrows.
    auto
    capacity() const -> size_t
    {
        return _capacity;
    }

    /// @brief Sets all values of matrix, padding included, to zero.
    auto
    zero() -> void
    {
        if (_data != nullptr) std::fill(_data, _data + number_of_elements(), 0.0);
    }

    /// @brief Gets values of matrix.
    /// @return The pointer to the values of matrix.
    auto
    data() -> double *
    {
        return _data;
    }

    /// @brief Gets values of matrix.
    /// @return The constant pointer to the values of matrix.
    auto
    data() const -> const double *
    {
        return _data;
    }

    /// @brief Gets values of specific row of matrix.
    /// @param row The index of row.
    /// @return The pointer to the values of row, aligned to a cache line boundary.
    auto
    data(const size_t row) -> double *
    {
        // NOTE: the message is built only when the check fails, as constructing
        // it eagerly allocates on every access of a row.

        if (row >= _rows) errors::assertMsgCritical(false, std::string("SimdMatrix.data: Index of row is out of range"));

        return _data + row * _pitch;
    }

    /// @brief Gets values of specific row of matrix.
    /// @param row The index of row.
    /// @return The constant pointer to the values of row, aligned to a cache line
    /// boundary.
    auto
    data(const size_t row) const -> const double *
    {
        if (row >= _rows) errors::assertMsgCritical(false, std::string("SimdMatrix.data: Index of row is out of range"));

        return _data + row * _pitch;
    }

    /// @brief Gets number of rows in matrix.
    /// @return The number of rows.
    auto
    number_of_rows() const -> size_t
    {
        return _rows;
    }

    /// @brief Gets number of columns in matrix.
    /// @return The number of columns.
    auto
    number_of_columns() const -> size_t
    {
        return _columns;
    }

    /// @brief Gets padded number of columns in a row of matrix.
    /// @return The padded number of columns.
    auto
    pitch() const -> size_t
    {
        return _pitch;
    }

    /// @brief Gets number of values in matrix, padding included.
    /// @return The number of values.
    auto
    number_of_elements() const -> size_t
    {
        return _rows * _pitch;
    }

    /// @brief Gets memory required to store the values of matrix.
    /// @return The memory in bytes.
    auto
    memory_size() const -> size_t
    {
        return number_of_elements() * sizeof(double);
    }

   public:
    /// @brief Class CBlockReuse turns on the reuse of the freed blocks of values
    /// for the calling thread while it is alive.
    /// @note The reuse is off unless a scope asks for it, and a driver which does
    /// not construct this guard allocates and frees exactly as it did before the
    /// cache existed. That is deliberate. A driver which forms a few matrices per
    /// block of atom pairs gains nothing from the reuse, as the allocations are
    /// already rare against the work of the block, and loses five to eight percent
    /// on the largest cases of the overlap and the kinetic energy to the memory
    /// the cache holds back from the allocator. A driver which forms them once per
    /// atom on c side, as the three-center electron repulsion does, gains a factor
    /// of three and a half on the same measurement. The two live in the same range
    /// of sizes, so no floor and no budget separates them and the scope has to say
    /// which of the two it is.
    /// @note The guard nests. An inner guard leaves the reuse on when it ends and
    /// the outermost one frees what the thread still holds, so a scope never keeps
    /// memory beyond itself.
    class CBlockReuse
    {
       public:
        /// @brief The constructor, which turns the reuse on for the thread.
        CBlockReuse()

            : _previous(_reusing())
        {
            _reusing() = true;
        }

        /// @brief The deleted copy constructor, as the guard owns a state of its
        /// thread.
        CBlockReuse(const CBlockReuse &other) = delete;

        /// @brief The deleted copy assignment operator.
        auto operator=(const CBlockReuse &other) -> CBlockReuse & = delete;

        /// @brief The destructor, which restores the reuse and frees the blocks
        /// the thread holds if no guard is left.
        ~CBlockReuse()
        {
            _reusing() = _previous;

            if (!_reusing()) _cache().clear();
        }

       private:
        /// @brief The state of the reuse before the guard was constructed.
        bool _previous;
    };

   private:
    /// @brief Class CBlockCache keeps the blocks of values a thread has freed, so
    /// that a matrix of a shape the thread has just freed takes its values back
    /// instead of asking the allocator of the system for them.
    /// @note The matrices of the integral drivers are formed and destroyed inside
    /// a parallel region, in a handful of shapes which repeat over the atoms and
    /// the blocks. The allocator of the system serializes the aligned requests of
    /// the threads, so those requests cost more wall time on all threads together
    /// than on one alone. The cache removes them from the parallel region: the
    /// blocks it holds belong to the thread and no other thread reaches them.
    /// @note The cache is a fallback and never a requirement. A shape it does not
    /// hold is allocated as it was before, and one too large for its budget is
    /// freed as it was before, so the matrices are correct whatever the cache
    /// holds and however much of it the thread has used.
    class CBlockCache
    {
       public:
        /// @brief The default constructor.
        CBlockCache() = default;

        /// @brief The deleted copy constructor, as a cache belongs to its thread.
        CBlockCache(const CBlockCache &other) = delete;

        /// @brief The deleted copy assignment operator.
        auto operator=(const CBlockCache &other) -> CBlockCache & = delete;

        /// @brief The destructor, which frees the blocks the thread still holds.
        ~CBlockCache()
        {
            for (size_t i = 0; i < _nentries; i++)
            {
                for (size_t j = 0; j < _entries[i].count; j++) _free(_entries[i].blocks[j]);
            }
        }

        /// @brief Frees every block the cache holds, leaving it empty.
        auto
        clear() -> void
        {
            for (size_t i = 0; i < _nentries; i++)
            {
                for (size_t j = 0; j < _entries[i].count; j++) _free(_entries[i].blocks[j]);

                _entries[i] = CEntry{};
            }

            _nentries = 0;

            _bytes = 0;
        }

        /// @brief Takes a block of the given size from the cache.
        /// @param nbytes The size of the block in bytes.
        /// @return The block, whose content is undefined, or null if the cache
        /// holds no block of that size.
        auto
        take(const size_t nbytes) -> double *
        {
            for (size_t i = 0; i < _nentries; i++)
            {
                if ((_entries[i].nbytes != nbytes) || (_entries[i].count == 0)) continue;

                _entries[i].tick = ++_tick;

                _bytes -= nbytes;

                return _entries[i].blocks[--_entries[i].count];
            }

            return nullptr;
        }

        /// @brief Gives a block of the given size to the cache.
        /// @param values The block to give.
        /// @param nbytes The size of the block in bytes.
        /// @return True if the cache took the block, false if the caller must
        /// free it.
        auto
        give(double *values, const size_t nbytes) -> bool
        {
            // NOTE: a block below the smallest size is left to the allocator of
            // the system, which keeps the small blocks of a thread in a cache of
            // its own and does not serialize their requests. A block above the
            // budget of the cache never fits and is left to it as well.

            if ((nbytes < _min_bytes) || (nbytes > _max_bytes)) return false;

            // NOTE: room is made by evicting the sizes which have gone unused for
            // the longest, rather than by refusing the block. A thread which moves
            // from one block of atom pairs to another of a different size would
            // otherwise hold the sizes of the block it has left until it ends, as
            // a block leaves the cache only when a matrix of its own size is
            // formed, and would never cache the sizes it has moved to.

            while ((_bytes + nbytes > _max_bytes) && _evict_oldest(true))
            {
            }

            if (_bytes + nbytes > _max_bytes) return false;

            for (size_t i = 0; i < _nentries; i++)
            {
                if (_entries[i].nbytes != nbytes) continue;

                if (_entries[i].count == _max_blocks) return false;

                _entries[i].blocks[_entries[i].count++] = values;

                _entries[i].tick = ++_tick;

                _bytes += nbytes;

                return true;
            }

            if ((_nentries == _max_entries) && !_evict_oldest(false)) return false;

            auto &entry = _entries[_nentries++];

            entry.nbytes    = nbytes;
            entry.blocks[0] = values;
            entry.count     = 1;
            entry.tick      = ++_tick;

            _bytes += nbytes;

            return true;
        }

       private:
        /// @brief The number of blocks of one size the cache holds. Two are kept
        /// rather than one, as a matrix is copied while the matrix it is copied
        /// from is alive.
        static constexpr size_t _max_blocks = 2;

        /// @brief The number of sizes the cache tracks. The shapes of a block of
        /// atom pairs are the coordinates and the solid harmonics of the angular
        /// momenta below the highest, so a handful of sizes covers a block.
        static constexpr size_t _max_entries = 16;

        /// @brief The smallest block the cache holds, in bytes.
        static constexpr size_t _min_bytes = 4096;

        /// @brief The memory the cache holds for its thread, in bytes.
        static constexpr size_t _max_bytes = 16 * 1024 * 1024;

        /// @brief Struct CEntry holds the blocks of one size.
        struct CEntry
        {
            /// @brief The size of the blocks in bytes.
            size_t nbytes = 0;

            /// @brief The number of blocks held.
            size_t count = 0;

            /// @brief The value of the counter when the entry was last used.
            size_t tick = 0;

            /// @brief The blocks held.
            double *blocks[_max_blocks] = {};
        };

        /// @brief Frees a block of the cache.
        /// @param values The block to free.
        static auto
        _free(double *values) -> void
        {
            ::operator delete[](values, std::align_val_t{simd::cache_line_size()});
        }

        /// @brief Evicts the size which has gone unused for the longest, freeing
        /// the blocks it holds.
        /// @param held True to evict only a size which holds blocks, as is wanted
        /// when the budget is what runs out, false to evict any size, as is
        /// wanted when the table is what runs out.
        /// @return True if a size was evicted.
        auto
        _evict_oldest(const bool held) -> bool
        {
            auto slot = _max_entries;

            for (size_t i = 0; i < _nentries; i++)
            {
                if (held && (_entries[i].count == 0)) continue;

                if ((slot == _max_entries) || (_entries[i].tick < _entries[slot].tick)) slot = i;
            }

            if (slot == _max_entries) return false;

            for (size_t j = 0; j < _entries[slot].count; j++) _free(_entries[slot].blocks[j]);

            _bytes -= _entries[slot].count * _entries[slot].nbytes;

            // NOTE: the last entry takes the slot of the evicted one, so that the
            // entries in use stay at the front of the table.

            _entries[slot] = _entries[--_nentries];

            _entries[_nentries] = CEntry{};

            return true;
        }

        /// @brief The entries of the cache.
        CEntry _entries[_max_entries] = {};

        /// @brief The number of entries in use.
        size_t _nentries = 0;

        /// @brief The counter which orders the entries by their last use.
        size_t _tick = 0;

        /// @brief The memory the cache holds, in bytes.
        size_t _bytes = 0;
    };

    /// @brief Gets the cache of blocks of the calling thread.
    /// @return The cache of the thread.
    static auto
    _cache() -> CBlockCache &
    {
        thread_local CBlockCache cache;

        return cache;
    }

    /// @brief Gets whether the calling thread reuses the blocks it frees.
    /// @return The reference to the state of the thread, which CBlockReuse sets.
    static auto
    _reusing() -> bool &
    {
        thread_local bool reusing = false;

        return reusing;
    }

    /// @brief Allocates the values of matrix, leaving their content undefined.
    auto
    _allocate() -> void
    {
        if (_capacity > 0)
        {
            const auto nbytes = _capacity * sizeof(double);

            if (_reusing())
            {
                if (auto *values = _cache().take(nbytes); values != nullptr)
                {
                    _data = values;

                    return;
                }
            }

            _data = static_cast<double *>(::operator new[](nbytes, std::align_val_t{simd::cache_line_size()}));
        }
    }

    /// @brief Deallocates the values of matrix.
    auto
    _deallocate() -> void
    {
        if (_data != nullptr)
        {
            if (_owned && (!_reusing() || !_cache().give(_data, _capacity * sizeof(double))))
            {
                ::operator delete[](_data, std::align_val_t{simd::cache_line_size()});
            }

            _data = nullptr;
        }
    }

    /// @brief The values of matrix, stored row wise with padded rows.
    double *_data;

    /// @brief The number of rows in matrix.
    size_t _rows;

    /// @brief The number of columns in matrix.
    size_t _columns;

    /// @brief The number of values matrix owns or borrows, which its shape may not
    /// exceed.
    size_t _capacity;

    /// @brief Whether matrix owns the values it holds and frees them.
    bool _owned;

    /// @brief The padded number of columns in a row of matrix.
    size_t _pitch;
};

#endif /* SimdMatrix_hpp */
