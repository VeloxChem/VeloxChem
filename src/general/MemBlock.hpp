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

#ifndef MemBlock_hpp
#define MemBlock_hpp

#include <mpi.h>

#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <sstream>
#include <string>
#include <vector>

#include "MathFunc.hpp"
#include "MemAlloc.hpp"
#include "MpiFunc.hpp"
#include "NumaPolicy.hpp"

/**
 Templated class CMemBlock manages contiguous memory block allocation,
 manipulation, and deallocation.

 @author Z. Rinkevicius
 */
template <class T>
class CMemBlock
{
    /**
     The pointer to contiguous memory block.
     */
    T* _data{nullptr};

    /**
     The number of elements in memory block.
     */
    int32_t _nElements{0};

    /**
     The NUMA policy for memory block handling.
     */
    numa _numaPolicy{numa::serial};

    /**
     Allocates memory block.
     */
    void _allocate();

    /**
     Copies aligned data to memory block from data source.

     @param source the pointer to data source.
     */
    void _copy_aligned(const T* source);

    /**
     Copies unaligned data to memory block from data source.

     @param source the pointer to data source.
     */
    void _copy_unaligned(const T* source);

   public:
    /**
     Creates an empty memory block object.
     */
    CMemBlock();

    /**
     Creates a memory block object.

     @param nElements - the number of elements.
     */
    explicit CMemBlock(const int32_t nElements);

    /**
     Creates a memory block object with specific NUMA policy.

     @param nElements - the number of elements.
     @param numaPolicy - the numa policy for data initialization.
     */
    CMemBlock(const int32_t nElements, const numa numaPolicy);

    /**
     Creates a memory block object.

     @param dataVector - the vector with data elements.
     */
    CMemBlock(const std::vector<T>& dataVector);

    /**
     Create a memory block object.

     @param source - a raw pointer
     @param size - size of data
     @param numaPolicy - the numa policy for data initialization.
     */
    CMemBlock(const T* source, int32_t size, numa numaPolicy);

    /**
     Creates a memory block object by copying other memory block object.

     @param source the memory block object.
     */
    CMemBlock(const CMemBlock<T>& source);

    /**
     Creates a memory block object by moving other memory block object.

     @param source the memory block object.
     */
    CMemBlock(CMemBlock<T>&& source) noexcept;

    /**
     Destroys a memory block object.
     */
    ~CMemBlock();

    /**
     Assigns a memory block object by copying other memory block object.

     @param source the memory block object.
     */
    CMemBlock<T>& operator=(const CMemBlock<T>& source);

    /**
     Assigns a memory block object by moving other memory block object.

     @param source the memory block object.
     */
    CMemBlock<T>& operator=(CMemBlock<T>&& source) noexcept;

    /**
     Compares memory block object with other memory block object.

     @param other the memory block object.
     @return true if memory block objects are equal, false otherwise.
     */
    bool operator==(const CMemBlock<T>& other) const;

    /**
     Compares memory block object with other memory block object.

     @param other the memory block object.
     @return true if memory block objects are not equal, false otherwise.
     */
    bool operator!=(const CMemBlock<T>& other) const;

    /**
     Gets pointer to memory block.

     @return the pointer to contiguous memory block.
     */
    T* data();

    /**
     Gets pointer to constant memory block.

     @return the pointer to constant contiguous memory block.
     */
    const T* data() const;

    /**
     Gets pointer to specific element in memory block.

     @param iElement the index of element in memory block.
     @return the pointer to element in memory block.
     */
    T* data(const int32_t iElement);

    /**
     Gets constant pointer to specific element in memory block.

     @param iElement the index of element in memory block.
     @return the constant pointer to element in memory block.
     */
    const T* data(const int32_t iElement) const;

    /**
     Gets number of elements in memory block.

     @return the number of elements.
     */
    int32_t size() const;

    /**
     Sets memory block elements to zero.
     */
    void zero();

    /**
     Gets reference to specific element of memory block.

     @param iElement the index of element in memory block.
     @return the reference to element.
     */
    T& at(const int32_t iElement);

    /**
     Gets constant reference to specific element of memory block.

     @param iElement the index of element in memory block.
     @return the reference to element.
     */
    const T& at(const int32_t iElement) const;

    /**
     Creates a memory block object by subdividing data in memory block object
     into specific number of subblocks and summing up all subblocks.

     @param nElements the number elements in subblock.
     @param nBlocks the number of subblocks.
     @return the memory block object.
     */
    CMemBlock<T> pack(const int32_t nElements, const int32_t nBlocks) const;

    /**
     Creates a memory block object by subdividing data in memory block object
     into specific number of subblocks and picking specific element from each
     subblock.

     @param nElements the number elements in subblock.
     @param nBlocks the number of subblocks.
     @param iPosition the position of picked element in subblock.
     @return the memory block object.
     */
    CMemBlock<T> pick(const int32_t nElements, const int32_t nBlocks, const int32_t iPosition) const;

    /**
     Shrinks memory block object by discarding all elements beyond given elements
     number. NOTE: No action taken if given elements number is equal or larger
     than number of elements in memory block object.

     @param nElements the number of elements.
     */
    void shrink(const int32_t nElements);

    /**
     Creates memory block object by slicing given number of elements at specific
     position from memory block object.

     @param iPosition the position of first sliced element.
     @param nElements the number of sliced elements.
     @return the memory block object with sliced elements.
     */
    CMemBlock<T> slice(const int32_t iPosition, const int32_t nElements) const;

    /**
     Broadcasts memory block object within domain of MPI communicator.

     @param rank the rank of MPI process.
     @param comm the MPI communicator.
     */
    void broadcast(int32_t rank, MPI_Comm comm);

    /**
     Creates memory block object on master MPI process by gathering memory
     block object from all MPI processes within domain of MPI communicator.

     @param rank the rank of MPI process.
     @param nodes the number of MPI processes in MPI communicator.
     @param comm the MPI communicator.
     @return the memory block object: (a) on master node with gathered data;
             (b) on worker nodes empty.
     */
    CMemBlock<T> gather(int32_t rank, int32_t nodes, MPI_Comm comm);

    /**
     Reasigns memory block object on all MPI process within domain of MPI
     communicator by scattering memory object from master MPI process.

     @param rank the rank of MPI process.
     @param nodes the number of MPI processes in MPI communicator.
     @param comm the MPI communicator.
     */
    void scatter(int32_t rank, int32_t nodes, MPI_Comm comm);

    /**
     Reduces memory block objects from all MPI process within domain of MPI
     communicator into memory block object on master node by summing them.

     @param rank the rank of MPI process.
     @param nodes the number of MPI processes in MPI communicator.
     @param comm the MPI communicator.
     */
    void reduce_sum(int32_t rank, int32_t nodes, MPI_Comm comm);

    /**
     Converts memory block object to text output.
     */
    std::string repr() const;
};

template <class T>
CMemBlock<T>::CMemBlock()

    : _data(nullptr)

    , _nElements(0)

    , _numaPolicy(numa::serial)
{
}

template <class T>
CMemBlock<T>::CMemBlock(const int32_t nElements)

    : _data(nullptr)

    , _nElements(nElements)

    , _numaPolicy(numa::serial)
{
    _allocate();

    zero();
}

template <class T>
CMemBlock<T>::CMemBlock(const int32_t nElements, const numa numaPolicy)

    : _data(nullptr)

    , _nElements(nElements)

    , _numaPolicy(numaPolicy)
{
    _allocate();

    zero();
}

template <class T>
CMemBlock<T>::CMemBlock(const std::vector<T>& dataVector)

    : _data(nullptr)

    , _nElements(static_cast<int32_t>(dataVector.size()))

    , _numaPolicy(numa::serial)
{
    _allocate();

    _copy_unaligned(dataVector.data());
}

template <typename T>
CMemBlock<T>::CMemBlock(const T* source, int32_t nElements, numa numaPolicy)

    : _data(nullptr)

    , _nElements(nElements)

    , _numaPolicy(numaPolicy)
{
    _allocate();

    _copy_unaligned(source);
}

template <class T>
CMemBlock<T>::CMemBlock(const CMemBlock<T>& source)

    : _data(nullptr)

    , _nElements(source._nElements)

    , _numaPolicy(source._numaPolicy)
{
    _allocate();

    _copy_aligned(source.data());
}

template <class T>
CMemBlock<T>::CMemBlock(CMemBlock<T>&& source) noexcept

    : _data(nullptr)

    , _nElements(std::move(source._nElements))

    , _numaPolicy(std::move(source._numaPolicy))
{
    _data = source._data;

    source._data = nullptr;
}

template <class T>
CMemBlock<T>::~CMemBlock()
{
    mem::free(_data);
}

template <class T>
CMemBlock<T>&
CMemBlock<T>::operator=(const CMemBlock<T>& source)
{
    if (this == &source) return *this;

    _nElements = source._nElements;

    _numaPolicy = source._numaPolicy;

    _allocate();

    _copy_aligned(source.data());

    return *this;
}

template <class T>
CMemBlock<T>&
CMemBlock<T>::operator=(CMemBlock<T>&& source) noexcept
{
    if (this == &source) return *this;

    _nElements = std::move(source._nElements);

    _numaPolicy = std::move(source._numaPolicy);

    mem::free(_data);

    _data = source._data;

    source._data = nullptr;

    return *this;
}

template <class T>
bool
CMemBlock<T>::operator==(const CMemBlock<T>& other) const
{
    if (_nElements != other._nElements) return false;

    if (_numaPolicy != other._numaPolicy) return false;

    for (int32_t i = 0; i < _nElements; i++)
    {
        if (std::fabs(_data[i] - other._data[i]) > 1.0e-13) return false;
    }

    return true;
}

template <class T>
bool
CMemBlock<T>::operator!=(const CMemBlock<T>& other) const
{
    return !(*this == other);
}

template <class T>
T*
CMemBlock<T>::data()
{
    return _data;
}

template <class T>
const T*
CMemBlock<T>::data() const
{
    return _data;
}

template <class T>
T*
CMemBlock<T>::data(const int32_t iElement)
{
    return &(_data[iElement]);
}

template <class T>
const T*
CMemBlock<T>::data(const int32_t iElement) const
{
    return &(_data[iElement]);
}

template <class T>
int32_t
CMemBlock<T>::size() const
{
    return _nElements;
}

template <class T>
void
CMemBlock<T>::zero()
{
    auto t0value = static_cast<T>(0);

    auto pdata = _data;

    if (_numaPolicy == numa::parallel)
    {
#pragma omp      parallel for schedule(static)
        for (int32_t i = 0; i < _nElements; i++)
            pdata[i] = t0value;
    }
    else
    {
#pragma omp simd aligned(pdata : VLX_ALIGN)
        for (int32_t i = 0; i < _nElements; i++)
            pdata[i] = t0value;
    }
}

template <class T>
inline T&
CMemBlock<T>::at(const int32_t iElement)
{
    return _data[iElement];
}

template <class T>
inline const T&
CMemBlock<T>::at(const int32_t iElement) const
{
    return _data[iElement];
}

template <class T>
CMemBlock<T>
CMemBlock<T>::pack(const int32_t nElements, const int32_t nBlocks) const
{
    auto nelem = nElements * nBlocks;

    if (nelem <= _nElements)
    {
        CMemBlock<T> mblock(nElements);

        mblock.zero();

        for (int32_t i = 0; i < nBlocks; i++)
        {
            for (int32_t j = 0; j < nElements; j++)
            {
                mblock.at(j) += _data[i * nElements + j];
            }
        }

        return mblock;
    }

    return CMemBlock<T>();
}

template <class T>
CMemBlock<T>
CMemBlock<T>::pick(const int32_t nElements, const int32_t nBlocks, const int32_t iPosition) const
{
    auto nelem = nElements * nBlocks;

    if ((nelem <= _nElements) && (iPosition < nElements))
    {
        CMemBlock<T> mblock(nBlocks);

        for (int32_t i = 0; i < nBlocks; i++)
        {
            mblock.at(i) = _data[i * nElements + iPosition];
        }

        return mblock;
    }

    return CMemBlock<T>();
}

template <class T>
void
CMemBlock<T>::shrink(const int32_t nElements)
{
    if (nElements < _nElements)
    {
        auto tvals = mem::malloc<T>(nElements);

        auto pdata = _data;

#pragma omp simd aligned(pdata, tvals : VLX_ALIGN)
        for (int32_t i = 0; i < nElements; i++)
            tvals[i] = pdata[i];

        mem::free(_data);

        _nElements = nElements;

        _data = tvals;
    }
}

template <class T>
CMemBlock<T>
CMemBlock<T>::slice(const int32_t iPosition, const int32_t nElements) const
{
    if ((iPosition + nElements) <= _nElements)
    {
        CMemBlock<T> mblock(nElements);

        auto pdat = mblock.data();

        for (int32_t i = 0; i < nElements; i++)
        {
            pdat[i] = _data[iPosition + i];
        }

        return mblock;
    }

    return CMemBlock<T>();
}

template <>
inline void
CMemBlock<int32_t>::broadcast(int32_t rank, MPI_Comm comm)
{
    if constexpr (ENABLE_MPI)
    {
        // broadcast numa policy

        int32_t nmpol = 0;

        if (rank == mpi::master()) nmpol = to_int(_numaPolicy);

        mpi::bcast(nmpol, comm);

        if (rank != mpi::master()) _numaPolicy = to_numa(nmpol);

        // broadcast memory block data

        mpi::bcast(_nElements, comm);

        if (rank != mpi::master()) _allocate();

        MPI_Barrier(comm);

        auto merror = MPI_Bcast(_data, _nElements, MPI_INT32_T, mpi::master(), comm);

        if (merror != MPI_SUCCESS) mpi::abort(merror, "broadcast(CMemBlock)");
    }
}

template <>
inline void
CMemBlock<double>::broadcast(int32_t rank, MPI_Comm comm)
{
    if constexpr (ENABLE_MPI)
    {
        // broadcast numa policy

        int32_t nmpol = 0;

        if (rank == mpi::master()) nmpol = to_int(_numaPolicy);

        mpi::bcast(nmpol, comm);

        if (rank != mpi::master()) _numaPolicy = to_numa(nmpol);

        // broadcast memory block data

        mpi::bcast(_nElements, comm);

        if (rank != mpi::master()) _allocate();

        MPI_Barrier(comm);

        auto merror = MPI_Bcast(_data, _nElements, MPI_DOUBLE, mpi::master(), comm);

        if (merror != MPI_SUCCESS) mpi::abort(merror, "broadcast(CMemBlock)");
    }
}

template <>
inline CMemBlock<int32_t>
CMemBlock<int32_t>::gather(int32_t rank, int32_t nodes, MPI_Comm comm)
{
    if constexpr (ENABLE_MPI)
    {
        if (nodes == 1) return CMemBlock<int32_t>(*this);

        CMemBlock<int32_t> mblock;

        // setup gathering pattern

        int32_t* bsizes = nullptr;

        if (rank == mpi::master()) bsizes = mem::malloc<int32_t>(nodes);

        mpi::gather(bsizes, _nElements, comm);

        // set indexes of data blocks on master node

        int32_t* bindexes = nullptr;

        if (rank == mpi::master())
        {
            bindexes = mem::malloc<int32_t>(nodes);

            mathfunc::indexes(bindexes, bsizes, nodes);

            mblock = CMemBlock<int32_t>(mathfunc::sum(bsizes, nodes));
        }

        // gather data on master node

        auto merror = MPI_Gatherv(_data, _nElements, MPI_INT32_T, mblock.data(), bsizes, bindexes, MPI_INT32_T, mpi::master(), comm);

        if (merror != MPI_SUCCESS) mpi::abort(merror, "gather(CMemBlock)");

        // deallocate gathering pattern data

        mem::free(bsizes);

        mem::free(bindexes);

        return mblock;
    }

    return CMemBlock<int32_t>(*this);
}

template <>
inline CMemBlock<double>
CMemBlock<double>::gather(int32_t rank, int32_t nodes, MPI_Comm comm)
{
    if constexpr (ENABLE_MPI)
    {
        if (nodes == 1) return CMemBlock<double>(*this);

        CMemBlock<double> mblock;

        // setup for gathering pattern

        int32_t* bsizes = nullptr;

        if (rank == mpi::master()) bsizes = mem::malloc<int32_t>(nodes);

        mpi::gather(bsizes, _nElements, comm);

        // set indexes of data blocks on master node

        int32_t* bindexes = nullptr;

        if (rank == mpi::master())
        {
            bindexes = mem::malloc<int32_t>(nodes);

            mathfunc::indexes(bindexes, bsizes, nodes);

            mblock = CMemBlock<double>(mathfunc::sum(bsizes, nodes));
        }

        // gather data on master node

        auto merror = MPI_Gatherv(_data, _nElements, MPI_DOUBLE, mblock.data(), bsizes, bindexes, MPI_DOUBLE, mpi::master(), comm);

        if (merror != MPI_SUCCESS) mpi::abort(merror, "gather(CMemBlock)");

        // deallocate gathering pattern data

        mem::free(bsizes);

        mem::free(bindexes);

        return mblock;
    }

    return CMemBlock<double>(*this);
}

template <>
inline void
CMemBlock<int32_t>::scatter(int32_t rank, int32_t nodes, MPI_Comm comm)
{
    if constexpr (ENABLE_MPI)
    {
        if (nodes == 1) return;

        // broadcast total number of elements

        int32_t totelem = _nElements;

        mpi::bcast(totelem, comm);

        // update number of elements

        _nElements = mpi::batch_size(totelem, rank, nodes);

        // allocate local data chunk

        auto bdata = mem::malloc<int32_t>(_nElements);

        // setup for scaterring pattern

        int32_t* bsizes = nullptr;

        int32_t* bindexes = nullptr;

        if (rank == mpi::master())
        {
            bsizes = mem::malloc<int32_t>(nodes);

            mpi::batches_pattern(bsizes, totelem, nodes);

            bindexes = mem::malloc<int32_t>(nodes);

            mathfunc::indexes(bindexes, bsizes, nodes);
        }

        // scatter data chunk

        auto merror = MPI_Scatterv(_data, bsizes, bindexes, MPI_INT32_T, bdata, _nElements, MPI_INT32_T, mpi::master(), comm);

        if (merror != MPI_SUCCESS) mpi::abort(merror, "scatter(CMemBlock)");

        // deallocate scattering pattern data

        mem::free(bsizes);

        mem::free(bindexes);

        // swap data chuck pointers

        if (_data != nullptr) mem::free(_data);

        _data = bdata;
    }
}

template <>
inline void
CMemBlock<double>::scatter(int32_t rank, int32_t nodes, MPI_Comm comm)
{
    if constexpr (ENABLE_MPI)
    {
        if (nodes == 1) return;

        // broadcast total number of elements

        int32_t totelem = _nElements;

        mpi::bcast(totelem, comm);

        // update number of elements

        _nElements = mpi::batch_size(totelem, rank, nodes);

        // allocate local data chunk

        auto bdata = mem::malloc<double>(_nElements);

        // setup for scaterring pattern

        int32_t* bsizes = nullptr;

        int32_t* bindexes = nullptr;

        if (rank == mpi::master())
        {
            bsizes = mem::malloc<int32_t>(nodes);

            mpi::batches_pattern(bsizes, totelem, nodes);

            bindexes = mem::malloc<int32_t>(nodes);

            mathfunc::indexes(bindexes, bsizes, nodes);
        }

        // scatter data chunk

        auto merror = MPI_Scatterv(_data, bsizes, bindexes, MPI_DOUBLE, bdata, _nElements, MPI_DOUBLE, mpi::master(), comm);

        if (merror != MPI_SUCCESS) mpi::abort(merror, "scatter(CMemBlock)");

        // deallocate scattering pattern data

        mem::free(bsizes);

        mem::free(bindexes);

        // swap data chuck pointers

        if (_data != nullptr) mem::free(_data);

        _data = bdata;
    }
}

template <>
inline void
CMemBlock<double>::reduce_sum(int32_t rank, int32_t nodes, MPI_Comm comm)
{
    if constexpr (ENABLE_MPI)
    {
        if (nodes == 1) return;

        double* bdata = nullptr;

        if (rank == mpi::master()) bdata = mem::malloc<double>(_nElements);

        mpi::reduce_sum(_data, bdata, _nElements, comm);

        if (rank == mpi::master())
        {
            mem::free(_data);

            _data = bdata;
        }
    }
}

template <class T>
void
CMemBlock<T>::_allocate()
{
    if (_nElements > 0)
    {
        if (_data != nullptr) mem::free(_data);

        auto numelem = static_cast<size_t>(_nElements);

        _data = mem::malloc<T>(numelem);
    }
}

template <class T>
void
CMemBlock<T>::_copy_unaligned(const T* source)
{
    auto pdata = _data;

    if (_numaPolicy == numa::parallel)
    {
#pragma omp      parallel for schedule(static)
        for (int32_t i = 0; i < _nElements; i++)
            pdata[i] = source[i];
    }
    else
    {
#pragma omp simd aligned(pdata : VLX_ALIGN)
        for (int32_t i = 0; i < _nElements; i++)
            pdata[i] = source[i];
    }
}

template <class T>
void
CMemBlock<T>::_copy_aligned(const T* source)
{
    auto pdata = _data;

    if (_numaPolicy == numa::parallel)
    {
#pragma omp      parallel for schedule(static)
        for (int32_t i = 0; i < _nElements; i++)
            pdata[i] = source[i];
    }
    else
    {
#pragma omp simd aligned(pdata, source : VLX_ALIGN)
        for (int32_t i = 0; i < _nElements; i++)
            pdata[i] = source[i];
    }
}

template <typename T>
std::string
CMemBlock<T>::repr() const
{
    std::ostringstream os;

    os << std::endl;

    os << "[CMemBlock (Object):" << this << "]" << std::endl;

    os << "_numaPolicy: " << to_string(_numaPolicy) << std::endl;

    os << "_nElements: " << _nElements << std::endl;

    os << "_data (" << &(_data) << "):" << std::endl;

    for (auto i = 0; i < _nElements; ++i)
    {
        os << " " << _data[i];
    }

    os << std::endl;

    return os.str();
}

template <typename T>
std::ostream&
operator<<(std::ostream& output, const CMemBlock<T>& source)
{
    return (output << source.repr());
}

#endif /* MemBlock_hpp */
