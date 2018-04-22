//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#ifndef MemBlock_hpp
#define MemBlock_hpp

#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <cmath>

#include "MemAlloc.hpp"
#include "MpiFunc.hpp"
#include "MathFunc.hpp"

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
    T* _data;

    /**
     The number of elements in memory block.
     */
    int32_t _nElements;

    /**
     Allocates memory block.
     */
    void _allocate();

    /**
     Copies data to memory block from data source.

     @param source the pointer to data source.
     */
    void _copy(const T* source);

public:

    /**
     Creates an empty memory block object.
     */
    CMemBlock();

    /**
     Creates a memory block object.
     
     @param nElements - the number of elements.
     */
    CMemBlock(const int32_t nElements);

    /**
     Creates a memory block object.
     
     @param dataVector - the vector with data elements.
     */
    CMemBlock(const std::vector<T>& dataVector);

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

     @param nSubBlocks the number of subblocks.
     @return the memory block object.
     */
    CMemBlock<T> pack(const int32_t nSubBlocks) const;
    
    /**
     Creates a memory block object by subdividing data in memory block object
     into specific number of subblocks and picking specific element from each
     subblock.

     @param nSubBlocks the number of subblocks.
     @param iPosition the starting position of subblock.
     @return the memory block object.
     */
    CMemBlock<T> pick(const int32_t nSubBlocks, const int32_t iPosition) const;
    
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
     Converts memory block object to text output and insert it into output text
     stream.
     
     @param output the output text stream.
     @param source the memory block object.
     */
    template <class U>
    friend std::ostream& operator<<(std::ostream& output,
                                    const CMemBlock<U>& source);
};

template <class T>
CMemBlock<T>::CMemBlock()

    : _data(nullptr)

    , _nElements(0)
{

}

template <class T>
CMemBlock<T>::CMemBlock(const int32_t nElements)

    : _data(nullptr)

    , _nElements(nElements)
{
    _allocate();
}

template<class T>
CMemBlock<T>::CMemBlock(const std::vector<T>& dataVector)

    : _data(nullptr)

    , _nElements(static_cast<int32_t>(dataVector.size()))
{
    _allocate();

    for (int32_t i = 0; i < _nElements; i++) _data[i] = dataVector[i];
}

template <class T>
CMemBlock<T>::CMemBlock(const CMemBlock<T>& source)

    : _data(nullptr)

    , _nElements(source._nElements)
{
    _allocate();

    _copy(source.data());
}

template <class T>
CMemBlock<T>::CMemBlock(CMemBlock<T>&& source) noexcept

    : _data(nullptr)

    , _nElements(std::move(source._nElements))
{
    _data = source._data;

    source._data = nullptr;
}

template<class T>
CMemBlock<T>::~CMemBlock()
{
    mem::free(_data);
}

template<class T>
CMemBlock<T>& CMemBlock<T>::operator=(const CMemBlock<T>& source)
{
    if (this == &source) return *this;

    _nElements = source._nElements;

    _allocate();

    _copy(source.data());

    return *this;
}

template<class T>
CMemBlock<T>& CMemBlock<T>::operator=(CMemBlock<T>&& source) noexcept
{
    if (this == &source) return *this;

    _nElements = std::move(source._nElements);

    mem::free(_data);

    _data = source._data;

    source._data = nullptr;

    return *this;
}

template<class T>
bool CMemBlock<T>::operator==(const CMemBlock<T>& other) const
{
    if (_nElements != other._nElements) return false;

    for (int32_t i = 0; i < _nElements; i++)
    {
        if (std::fabs(_data[i] - other._data[i]) > 1.0e-13) return false;
    }

    return true;
}

template<class T>
bool CMemBlock<T>::operator!=(const CMemBlock<T>& other) const
{
    return !(*this == other);
}

template<class T>
T* CMemBlock<T>::data()
{
    return _data;
}

template<class T>
const T* CMemBlock<T>::data() const
{
    return _data;
}

template<class T>
T* CMemBlock<T>::data(const int32_t iElement)
{
    return &(_data[iElement]);
}

template<class T>
const T* CMemBlock<T>::data(const int32_t iElement) const
{
    return &(_data[iElement]);
}

template <class T>
int32_t CMemBlock<T>::size() const
{
    return _nElements;
}

template <class T>
void CMemBlock<T>::zero()
{
    auto t0value = static_cast<T>(0);

    auto pdata = _data;

    #pragma omp simd aligned(pdata:VLX_ALIGN)
    for (int32_t i = 0; i < _nElements; i++) pdata[i] = t0value;
}

template <class T>
inline T& CMemBlock<T>::at(const int32_t iElement)
{
    return _data[iElement];
}

template <class T>
inline const T& CMemBlock<T>::at(const int32_t iElement) const 
{
    return _data[iElement];
}

template <class T>
CMemBlock<T> CMemBlock<T>::pack(const int32_t nSubBlocks) const
{
    auto numelem = _nElements / nSubBlocks;
    
    if (numelem > 0)
    {
        CMemBlock<T> mblock(numelem);
        
        mblock.zero();
        
        for (int32_t i = 0; i < nSubBlocks; i++)
        {
            for (int32_t j = 0; j < numelem; j++)
            {
                mblock.at(j) += _data[i * numelem + j];
            }
        }
        
        auto nremind = _nElements % nSubBlocks;
        
        if (nremind != 0)
        {
            for (int32_t i = 0; i < nremind; i++)
            {
                mblock.at(i) += _data[nSubBlocks * numelem + i]; 
            }
        }
        
        return mblock;
    }
    
    return CMemBlock<T>(*this);
}

template <class T>
CMemBlock<T> CMemBlock<T>::pick(const int32_t nSubBlocks,
                                const int32_t iPosition) const
{
    auto numelem = _nElements / nSubBlocks;
    
    if ((numelem > 0) && (nSubBlocks > 1))
    {
        CMemBlock<T> mblock(numelem);
        
        for (int32_t i = 0; i < numelem; i++)
        {
            mblock.at(i) = _data[i * numelem + iPosition];
        }
        
        return mblock;
    }
    
    return CMemBlock<T>(*this);
}

template <>
inline void CMemBlock<int32_t>::broadcast(int32_t rank, MPI_Comm comm)
{
    if (ENABLE_MPI)
    {
        mpi::bcast(_nElements, comm);
        
        if (rank != mpi::master()) _allocate();
        
        auto merror = MPI_Bcast(_data, _nElements, MPI_INT32_T, mpi::master(),
                                comm);
        
        if (merror != MPI_SUCCESS) mpi::abort(merror, "broadcast(CMemBlock)");
    }
}

template <>
inline void CMemBlock<double>::broadcast(int32_t rank, MPI_Comm comm)
{
    if (ENABLE_MPI)
    {
        mpi::bcast(_nElements, comm);
        
        if (rank != mpi::master()) _allocate();
        
        auto merror = MPI_Bcast(_data, _nElements, MPI_DOUBLE, mpi::master(),
                                comm);
        
        if (merror != MPI_SUCCESS) mpi::abort(merror, "broadcast(CMemBlock)");
    }
}

template <>
inline CMemBlock<int32_t>
CMemBlock<int32_t>::gather(int32_t rank, int32_t nodes, MPI_Comm comm)
{
    if (ENABLE_MPI)
    {
        if (nodes == 1) return CMemBlock<int32_t>(*this);
        
        CMemBlock<int32_t> mblock;
        
        // setup gathering pattern
        
        size_t nsizes = static_cast<size_t>(nodes) * sizeof(int32_t);

        int32_t* bsizes = nullptr;

        if (rank == mpi::master())
        {
            bsizes = (int32_t*) mem::malloc(nsizes);
        }

        mpi::gather(bsizes, _nElements, rank, comm);
        
        // set indexes of data blocks on master node
        
        int32_t* bindexes = nullptr;

        if (rank == mpi::master())
        {
            bindexes = (int32_t*) mem::malloc(nsizes);

            mathfunc::indexes(bindexes, bsizes, nodes);

            mblock = CMemBlock<int32_t>(mathfunc::sum(bsizes, nodes));
        }

        // gather data on master node

        auto merror = MPI_Gatherv(_data, _nElements, MPI_INT32_T, mblock.data(),
                                  bsizes, bindexes, MPI_INT32_T, mpi::master(),
                                  comm);

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
    if (ENABLE_MPI)
    {
        if (nodes == 1) return CMemBlock<double>(*this);
        
        CMemBlock<double> mblock;

        // setup for gathering pattern
        
        size_t nsizes = static_cast<size_t>(nodes) * sizeof(int32_t);

        int32_t* bsizes = nullptr;

        if (rank == mpi::master()) bsizes = (int32_t*) mem::malloc(nsizes);

        mpi::gather(bsizes, _nElements, rank, comm);

        // set indexes of data blocks on master node

        int32_t* bindexes = nullptr;

        if (rank == mpi::master())
        {
            bindexes = (int32_t*) mem::malloc(nsizes);

            mathfunc::indexes(bindexes, bsizes, nodes);

            mblock = CMemBlock<double>(mathfunc::sum(bsizes, nodes));
        }

        // gather data on master node

        auto merror = MPI_Gatherv(_data, _nElements, MPI_DOUBLE, mblock.data(),
                                  bsizes, bindexes, MPI_DOUBLE, mpi::master(),
                                  comm);

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
    if (ENABLE_MPI)
    {
        if (nodes == 1) return;
        
        // broadcast total number of elements
        
        int32_t totelem = _nElements;

        mpi::bcast(totelem, comm);

        // update number of elements

        _nElements = mpi::batch_size(totelem, rank, nodes);

        // allocate local data chunk

        auto nsizes = static_cast<size_t>(_nElements) * sizeof(int32_t);

        int32_t* bdata = (int32_t*) mem::malloc(nsizes);

        // setup for scaterring pattern

        int32_t* bsizes = nullptr;

        int32_t* bindexes = nullptr;

        if (rank == mpi::master())
        {
            nsizes = static_cast<size_t>(nodes) * sizeof(int32_t);

            bsizes = (int32_t*) mem::malloc(nsizes);

            mpi::batches_pattern(bsizes, totelem, nodes);

            bindexes = (int32_t*) mem::malloc(nsizes);

            mathfunc::indexes(bindexes, bsizes, nodes);
        }

        // scatter data chunk

        auto merror = MPI_Scatterv(_data, bsizes, bindexes, MPI_INT32_T, bdata,
                                   _nElements, MPI_INT32_T, mpi::master(), comm);

        if (merror != MPI_SUCCESS) mpi::abort(merror, "scatter(CMemBlock)");
        
        // deallocate scattering pattern data

        mem::free(bsizes);

        mem::free(bindexes);

        // swap data chuck pointers

        if (_data != nullptr)  mem::free(_data);

        _data = bdata;
    }
}

template <>
inline void
CMemBlock<double>::scatter(int32_t rank, int32_t nodes, MPI_Comm comm)
{
    if (ENABLE_MPI)
    {
        if (nodes == 1) return;
        
        // broadcast total number of elements
        
        int32_t totelem = _nElements;

        mpi::bcast(totelem, comm);

        // update number of elements

        _nElements = mpi::batch_size(totelem, rank, nodes);

        // allocate local data chunk

        auto nsizes =  static_cast<size_t>(_nElements) * sizeof(double);

        double* bdata = (double*) mem::malloc(nsizes);

        // setup for scaterring pattern

        int32_t* bsizes = nullptr;

        int32_t* bindexes = nullptr;

        if (rank == mpi::master())
        {
            nsizes = static_cast<size_t>(nodes) * sizeof(int32_t);

            bsizes = (int32_t*) mem::malloc(nsizes);

            mpi::batches_pattern(bsizes, totelem, nodes);

            bindexes = (int32_t*) mem::malloc(nsizes);

            mathfunc::indexes(bindexes, bsizes, nodes);
        }

        // scatter data chunk

        auto merror = MPI_Scatterv(_data, bsizes, bindexes, MPI_DOUBLE, bdata,
                                   _nElements, MPI_DOUBLE, mpi::master(), comm);

        if (merror != MPI_SUCCESS) mpi::abort(merror, "scatter(CMemBlock)");
        
        // deallocate scattering pattern data
        
        mem::free(bsizes);
        
        mem::free(bindexes);
        
        // swap data chuck pointers

        if (_data != nullptr)  mem::free(_data);

        _data = bdata;
    }
}

template <class T>
void CMemBlock<T>::_allocate()
{
    if (_nElements > 0)
    {
        if (_data != nullptr) mem::free(_data); 

        auto numelem = static_cast<size_t>(_nElements);

        _data = (T*) mem::malloc(numelem * sizeof(T));
    }
}

template <class T>
void CMemBlock<T>::_copy(const T* source)
{
    auto pdata = _data;

    #pragma omp simd aligned(pdata:VLX_ALIGN, source:VLX_ALIGN)
    for (int32_t i = 0; i < _nElements; i++) pdata[i] = source[i];
}


template<class U>
std::ostream& operator<<(std::ostream& output, const CMemBlock<U>& source)
{
    output << std::endl;

    output << "[CMemBlock (Object):" << &source << "]" << std::endl;

    output << "_nElements: " << source._nElements << std::endl;

    output << "_data (" << &(source._data) << "):" <<  std::endl;

    for (int32_t i = 0; i < source._nElements; i++)
    {
        output << " " << source._data[i];
    }

    output << std::endl;

    return output;
}

#endif /* MemBlock_hpp */
