//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#ifndef MemBlock2D_hpp
#define MemBlock2D_hpp

#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <cmath>

#include "mpi.h"

#include "MemAlloc.hpp"
#include "MpiFunc.hpp"
#include "MemBlock.hpp"
#include "MatOrder.hpp"

/**
 Templated class CMemBlock2D manages 2D memory block allocation, manipulation,
 and deallocation.
 
 @author Z. Rinkevicius
 */
template <class T>
class CMemBlock2D
{
    /**
     The contiguous memory block.
     */
    CMemBlock<T> _data;
    
    /**
     The original sizes of data chunks in memory block.
     */
    CMemBlock<int32_t> _originalSizes;
    
    /**
     The padded sizes of data chunks in memory block.
     */
    CMemBlock<int32_t> _paddedSizes;
    
    /**
     The start positions of data chunks in memory block.
     */
    CMemBlock<int32_t> _positions;
    
    /**
     The number of elements in memory block.
     */
    int32_t _nElements;
    
    /**
     Sets uniform original size for all data chunks.

     @param nElements the number of elements in data chunk.
     @param nBlocks the number of data chunks.
     */
    void _setOriginalSizes(const int32_t nElements,
                           const int32_t nBlocks);
    
    /**
     Sets dimensions, padding for all data chunks.
     */
    void _setDimensions();
    
    /**
     Copies data elements from vector to memory block.

     @param dataVector the vector with data elements.
     @param dataOrder the order for copying data elements from vector to memory.
     */
    void _copy(const std::vector<T>& dataVector,
               const matorder        dataOrder);
    
public:

    /**
     Creates an empty 2D memory block object.
     */
    CMemBlock2D();

    /**
     Creates an 2D memory block object.
     
     @param nElements the number of elements in data chunk.
     @param nBlocks the number of data chunks.
     */
    CMemBlock2D(const int32_t nElements,
                const int32_t nBlocks);
    
    /**
     Creates an 2D memory block object.
     
     @param dimVector the vector of data chunk sizes.
     */
    CMemBlock2D(const std::vector<int32_t>& dimVector);
    
    /**
     Creates an 2D memory block object.
     
     @param originalSizes the memory block object with data chunk sizes.
     */
    CMemBlock2D(const CMemBlock<int32_t> originalSizes);

    /**
     Creates an 2D memory block object.
     
     @param dataVector the vector with data elements.
     @param nElements the number of elements in data chunk.
     @param nBlocks the number of data chunks.
     */
    CMemBlock2D(const std::vector<T>& dataVector,
                const int32_t         nElements,
                const int32_t         nBlocks);
    
    /**
     Creates an 2D memory block object.
     
     @param dataVector the vector with data elements.
     @param dataOrder the order for copying data elements from vector to memory.
     @param nElements the number of elements in data chunk.
     @param nBlocks the number of data chunks.
     */
    CMemBlock2D(const std::vector<T>& dataVector,
                const matorder        dataOrder,
                const int32_t         nElements,
                const int32_t         nBlocks);
    
    /**
     Creates an 2D memory block object.
     
     @param dataVector the vector with data elements.
     @param dimVector the vector of data chunk sizes.
     */
    CMemBlock2D(const std::vector<T>&       dataVector,
                const std::vector<int32_t>& dimVector);
    
    /**
     Creates an 2D memory block object.
     
     @param dataVector the vector with data elements.
     @param dataOrder the order for copying data elements from vector to memory.
     @param dimVector the vector of data chunk sizes.
     */
    CMemBlock2D(const std::vector<T>&       dataVector,
                const matorder              dataOrder,
                const std::vector<int32_t>& dimVector);

    /**
     Creates a 2D memory block object by copying other 2D memory block object.
     
     @param source the 2D memory block object.
     */
    CMemBlock2D(const CMemBlock2D<T>& source);

    /**
     Creates a 2D memory block object by moving other 2D memory block object.
     
     @param source the 2D memory block object.
     */
    CMemBlock2D(CMemBlock2D<T>&& source) noexcept;

    /**
     Destroys a 2D memory block object.
     */
    ~CMemBlock2D();

    /**
     Assigns a 2D memory block object by copying other 2D memory block object.
     
     @param source the 2D memory block object.
     */
    CMemBlock2D<T>& operator=(const CMemBlock2D<T>& source);

    /**
     Assigns a 2D memory block object by moving other 2D memory block object.
     
     @param source the 2D memory block object.
     */
    CMemBlock2D<T>& operator=(CMemBlock2D<T>&& source) noexcept;

    /**
     Compares 2D memory block object with other 2D memory block object.
     
     @param other the 2D memory block object.
     @return true if 2D memory block objects are equal, false otherwise.
     */
    bool operator==(const CMemBlock2D<T>& other) const;

    /**
     Compares 2D memory block object with other 2D memory block object.
     
     @param other the 2D memory block object.
     @return true if 2D memory block objects are not equal, false otherwise.
     */
    bool operator!=(const CMemBlock2D<T>& other) const;

    /**
     Creates 2D memory block object by slicing part of this 2D memory block
     object.

     @param iPosition the position of first element in data chunks.
     @param nElements the number of elements in data chunks.
     @return the 2D memory block object.
     */
    CMemBlock2D<T> slice(const int32_t iPosition,
                         const int32_t nElements) const;
    
    /**
     Sets all elements of contiguous memory block to zero.
     */
    void zero();
    
    /**
     Checks if 2D memory block object is empty.

     @return true if 2D memory block object is empty, false otherwise.
     */
    bool isEmpty() const;
    
    /**
     Gets pointer to specific data chunk in 2D memory block.
     
     @param iBlock the index of data chunk.
     @return the pointer to data chunk.
     */
    T* data(const int32_t iBlock);
    
    /**
     Gets constant pointer to specific data chunk in 2D memory block.
     
     @param iBlock the index of data chunk.
     @return the pointer to data chunk.
     */
    const T* data(const int32_t iBlock) const;
    
    /**
     Gets pointer to specific element of selected data chunk in 2D memory block.
     
     @param iBlock the index of data chunk.
     @param iElement the index of element in memory block.
     @return the pointer to data chunk.
     */
    T* data(const int32_t iBlock,
            const int32_t iElement);
    
    /**
     Gets constant pointer to specific elment of selected data chunk in 2D
     memory block.
     
     @param iBlock the index of data chunk.
     @param iElement the index of element in memory block.
     @return the pointer to data chunk.
     */
    const T* data(const int32_t iBlock,
                  const int32_t iElement) const;
    
    /**
     Gets number of elements in specific data chunk.

     @param iBlock the index of data chunk.
     @return the number of elements in data chunk.
     */
    int32_t size(const int32_t iBlock) const;
    
    /**
     Gets number of data chunk.
     
     @return the number of data chunk.
     */
    int32_t blocks() const;
    
    /**
     Broadcasts 2D memory block within domain of MPI communicator.
     
     @param rank the rank of MPI process.
     @param comm the MPI communicator.
     */
    void broadcast(int32_t  rank,
                   MPI_Comm comm);
    
    /**
     Creates 2D memory block object on master MPI process by gathering 2D memory
     block object from all MPI processes within domain of MPI communicator.
     
     @param rank the rank of MPI process.
     @param nodes the number of MPI processes in MPI communicator.
     @param comm the MPI communicator.
     @return the 2D memory block object: (a) on master node with gathered data;
             (b) on worker nodes empty.
     */
    CMemBlock2D<T> gather(int32_t  rank,
                          int32_t  nodes,
                          MPI_Comm comm);
    
    /**
     Reasigns 2D memory block object on all MPI process within domain of MPI
     communicator by scattering 2D memory object from master MPI process.
     
     @param rank the rank of MPI process.
     @param nodes the number of MPI processes in MPI communicator.
     @param comm the MPI communicator.
     */
    void scatter(int32_t  rank,
                 int32_t  nodes,
                 MPI_Comm comm);
    
    /**
     Converts 2D memory block object to text output and insert it into output
     text stream.
     
     @param output the output text stream.
     @param source the 2D memory block object.
     */
    template <class U>
    friend std::ostream& operator<<(      std::ostream&   output,
                                    const CMemBlock2D<U>& source);
};

template <class T>
CMemBlock2D<T>::CMemBlock2D()

    : _nElements(0)
{

}

template <class T>
CMemBlock2D<T>::CMemBlock2D(const int32_t nElements,
                            const int32_t nBlocks)

    : _nElements(0)
{
    _setOriginalSizes(nElements, nBlocks);
    
    _setDimensions();
    
    _data = CMemBlock<T>(_nElements);
    
    _data.zero();
}

template<class T>
CMemBlock2D<T>::CMemBlock2D(const std::vector<int32_t>& dimVector)
{
    _originalSizes = CMemBlock<int32_t>(dimVector);
    
    _setDimensions();
    
    _data = CMemBlock<T>(_nElements);
    
    _data.zero();
}

template<class T>
CMemBlock2D<T>::CMemBlock2D(const CMemBlock<int32_t> originalSizes)
{
    _originalSizes = originalSizes;
    
    _setDimensions();
    
    _data = CMemBlock<T>(_nElements);
    
    _data.zero();
}

template<class T>
CMemBlock2D<T>::CMemBlock2D(const std::vector<T>& dataVector,
                            const int32_t         nElements,
                            const int32_t         nBlocks)

    : CMemBlock2D<T>(dataVector, matorder::row_major, nElements, nBlocks)
{
    
}

template<class T>
CMemBlock2D<T>::CMemBlock2D(const std::vector<T>& dataVector,
                            const matorder        dataOrder,
                            const int32_t         nElements,
                            const int32_t         nBlocks)

    : _nElements(0)
{
    _setOriginalSizes(nElements, nBlocks);
    
    _setDimensions();
    
    _data = CMemBlock<T>(_nElements);
    
    _data.zero();
    
    _copy(dataVector, dataOrder);
}

template<class T>
CMemBlock2D<T>::CMemBlock2D(const std::vector<T>&       dataVector,
                            const std::vector<int32_t>& dimVector)

    : CMemBlock2D<T>(dataVector, matorder::row_major, dimVector)
{
    
}

template<class T>
CMemBlock2D<T>::CMemBlock2D(const std::vector<T>&       dataVector,
                            const matorder              dataOrder,
                            const std::vector<int32_t>& dimVector)
{
    _originalSizes = CMemBlock<int32_t>(dimVector);
    
    _setDimensions();
    
    _data = CMemBlock<T>(_nElements);
    
    _data.zero();
    
    _copy(dataVector, dataOrder);
}

template <class T>
CMemBlock2D<T>::CMemBlock2D(const CMemBlock2D<T>& source)

    : _data(source._data)

    , _originalSizes(source._originalSizes)

    , _paddedSizes(source._paddedSizes)

    , _positions(source._positions)

    , _nElements(source._nElements)
{
}

template <class T>
CMemBlock2D<T>::CMemBlock2D(CMemBlock2D<T>&& source) noexcept

    : _data(std::move(source._data))

    , _originalSizes(std::move(source._originalSizes))

    , _paddedSizes(std::move(source._paddedSizes))

    , _positions(std::move(source._positions))

    , _nElements(std::move(source._nElements))
{

}

template<class T>
CMemBlock2D<T>::~CMemBlock2D()
{
    
}

template<class T>
CMemBlock2D<T>&
CMemBlock2D<T>::operator=(const CMemBlock2D<T>& source)
{
    if (this == &source) return *this;

    _nElements = source._nElements;
    
    _positions = source._positions;
    
    _paddedSizes = source._paddedSizes;
    
    _originalSizes = source._originalSizes;
    
    _data = source._data;
   
    return *this;
}

template<class T>
CMemBlock2D<T>&
CMemBlock2D<T>::operator=(CMemBlock2D<T>&& source) noexcept
{
    if (this == &source) return *this;

    _nElements = std::move(source._nElements);
    
    _positions = std::move(source._positions);
    
    _paddedSizes = std::move(source._paddedSizes);
    
    _originalSizes = std::move(source._originalSizes);
    
    _data = std::move(source._data);

    return *this;
}

template<class T>
bool
CMemBlock2D<T>::operator==(const CMemBlock2D<T>& other) const
{
    if (_nElements != other._nElements) return false;
    
    if (_positions != other._positions) return false;
    
    if (_paddedSizes != other._paddedSizes) return false;
    
    if (_originalSizes != other._originalSizes) return false;
    
    if (_data != other._data) return false;

    return true;
}

template<class T>
bool
CMemBlock2D<T>::operator!=(const CMemBlock2D<T>& other) const
{
    return !(*this == other);
}

template <class T>
CMemBlock2D<T>
CMemBlock2D<T>::slice(const int32_t iPosition,
                      const int32_t nElements) const
{
    CMemBlock2D<T> mblock(nElements, blocks());
        
    for (int32_t i = 0; i < blocks(); i++)
    {
        // set up pointers to data chunks
        
        auto idata = data(i, iPosition);
        
        auto odata = mblock.data(i);
        
        // copy elements of data chunks 
        
        for (int32_t j = 0; j < nElements; j++) odata[j] = idata[j];
    }
    
    return mblock;
}

template <class T>
void
CMemBlock2D<T>::zero()
{
    _data.zero(); 
}

template <class T>
bool
CMemBlock2D<T>::isEmpty() const
{
    if (_nElements == 0) return true;
    
    return false;
}

template<class T>
T*
CMemBlock2D<T>::data(const int32_t iBlock)
{
    if (_originalSizes.size() > 0)
    {
        if (iBlock < 0) return nullptr;
        
        if (iBlock >= blocks()) return nullptr;
        
        return _data.data(_positions.at(iBlock));
    }
    
    return nullptr;
}

template<class T>
const T*
CMemBlock2D<T>::data(const int32_t iBlock) const
{
    if (_originalSizes.size() > 0)
    {
        if (iBlock < 0) return nullptr;
        
        if (iBlock >= blocks()) return nullptr;
        
        return _data.data(_positions.at(iBlock));
    }
    
    return nullptr;
}

template<class T>
T*
CMemBlock2D<T>::data(const int32_t iBlock,
                     const int32_t iElement)
{
    if (_originalSizes.size() > 0)
    {
        auto pdata = _data.data(_positions.at(iBlock));
        
        if (iElement < _originalSizes.at(iBlock))
        {
            return &(pdata[iElement]);
        }
        
        return nullptr;
    }
    
    return nullptr;
}

template<class T>
const T*
CMemBlock2D<T>::data(const int32_t iBlock,
                     const int32_t iElement) const
{
    if (_originalSizes.size() > 0)
    {
        auto pdata = _data.data(_positions.at(iBlock));
        
        if (iElement < _originalSizes.at(iBlock))
        {
            return &(pdata[iElement]);
        }
        
        return nullptr;
    }
    
    return nullptr;
}

template <class T>
int32_t
CMemBlock2D<T>::size(const int32_t iBlock) const
{
    if (iBlock < _originalSizes.size()) return _originalSizes.at(iBlock);
    
    return 0;
}

template <class T>
int32_t
CMemBlock2D<T>::blocks() const
{
    return _originalSizes.size();
}

template <>
inline void
CMemBlock2D<int32_t>::broadcast(int32_t  rank,
                                MPI_Comm comm)
{
    if (ENABLE_MPI)
    {
        mpi::bcast(_nElements, comm);
        
        _positions.broadcast(rank, comm);
        
        _paddedSizes.broadcast(rank, comm);
        
        _originalSizes.broadcast(rank, comm);
        
        _data.broadcast(rank, comm);
    }
}

template <>
inline void
CMemBlock2D<double>::broadcast(int32_t  rank,
                               MPI_Comm comm)
{
    if (ENABLE_MPI)
    {
        mpi::bcast(_nElements, comm);
        
        _positions.broadcast(rank, comm);
        
        _paddedSizes.broadcast(rank, comm);
        
        _originalSizes.broadcast(rank, comm);
        
        _data.broadcast(rank, comm);
    }
}

template <>
inline CMemBlock2D<int32_t>
CMemBlock2D<int32_t>::gather(int32_t  rank,
                             int32_t  nodes,
                             MPI_Comm comm)
{
    if (ENABLE_MPI)
    {
        if (nodes == 1) return CMemBlock2D<int32_t>(*this);
        
        // get dimensions from nodes
        
        auto orgsizes = _originalSizes.gather(rank, nodes, comm);
        
        // create 2D memory block
        
        CMemBlock2D<int32_t> mblock;
        
        if (rank == mpi::master())
        {
            auto nsizes = orgsizes.pack(blocks(), nodes);
            
            mblock = CMemBlock2D<int32_t>(nsizes);
        }
        
        // set up gathering pattern
        
        CMemBlock<int32_t> bsizes;
        
        CMemBlock<int32_t> bindexes;
        
        if (rank == mpi::master()) bindexes = CMemBlock<int32_t>(nodes);
        
        // gather data chunks
        
        for (int32_t i = 0; i < blocks(); i++)
        {
            if (rank == mpi::master())
            {
                bsizes = orgsizes.pick(blocks(), nodes, i);
                
                mathfunc::indexes(bindexes.data(), bsizes.data(), nodes);
            }
           
            // gather data on master node
            
            auto merror = MPI_Gatherv(data(i), size(i), MPI_INT32_T,
                                      mblock.data(i), bsizes.data(),
                                      bindexes.data(), MPI_INT32_T,
                                      mpi::master(), comm);
            
            if (merror != MPI_SUCCESS) mpi::abort(merror, "gather(CMemBlock2D)");
          
            MPI_Barrier(comm);
        }
        
        return mblock;
    }

    return CMemBlock2D<int32_t>(*this);
}

template <>
inline CMemBlock2D<double>
CMemBlock2D<double>::gather(int32_t  rank,
                            int32_t  nodes,
                            MPI_Comm comm)
{
    if (ENABLE_MPI)
    {
        if (nodes == 1) return CMemBlock2D<double>(*this);
        
        // get dimensions from nodes
        
        auto orgsizes = _originalSizes.gather(rank, nodes, comm);
        
        // create 2D memory block
        
        CMemBlock2D<double> mblock;
        
        if (rank == mpi::master())
        {
            auto nsizes = orgsizes.pack(blocks(), nodes);
            
            mblock = CMemBlock2D<double>(nsizes);
        }
        
        // set up gathering pattern
        
        CMemBlock<int32_t> bsizes;
        
        CMemBlock<int32_t> bindexes;
        
        if (rank == mpi::master()) bindexes = CMemBlock<int32_t>(nodes);
        
        // gather data chunks
        
        for (int32_t i = 0; i < blocks(); i++)
        {
            if (rank == mpi::master())
            {
                bsizes = orgsizes.pick(blocks(), nodes, i);
                
                mathfunc::indexes(bindexes.data(), bsizes.data(), nodes);
            }
            
            // gather data on master node
            
            auto merror = MPI_Gatherv(data(i), size(i), MPI_DOUBLE,
                                      mblock.data(i), bsizes.data(),
                                      bindexes.data(), MPI_DOUBLE,
                                      mpi::master(), comm);
            
            if (merror != MPI_SUCCESS) mpi::abort(merror, "gather(CMemBlock2D)");
            
            MPI_Barrier(comm);
        }
        
        return mblock;
    }
    
    return CMemBlock2D<double>(*this);
}

template <>
inline void
CMemBlock2D<int32_t>::scatter(int32_t  rank,
                              int32_t  nodes,
                              MPI_Comm comm)
{
    if (ENABLE_MPI)
    {
        if (nodes == 1) return;
        
        // save old dimensions
        
        CMemBlock<int32_t> oldsizes;
        
        CMemBlock<int32_t> oldpositions;
        
        if (rank == mpi::master())
        {
            oldsizes = _originalSizes;
            
            oldpositions = _positions;
        }
        
        // update dimensions on all nodes
        
        _originalSizes.broadcast(rank, comm);
        
        for (int32_t i = 0; i < _originalSizes.size(); i++)
        {
            _originalSizes.at(i) = mpi::batch_size(_originalSizes.at(i),
                                                   rank, nodes);
        }
        
        _setDimensions();
        
        // set up scattered data storage
        
        CMemBlock<int32_t> mblock(_nElements);
        
        mblock.zero(); 
        
        // allocate scaterring pattern data
        
        CMemBlock<int32_t> bsizes;
        
        CMemBlock<int32_t> bindexes;
        
        if (rank == mpi::master())
        {
            bsizes = CMemBlock<int32_t>(nodes);
            
            bindexes = CMemBlock<int32_t>(nodes);
        }
        
        // scatter data chunks
        
        for (int32_t i = 0; i < _originalSizes.size(); i++)
        {
            // set up scattering pattern for data chunk
            
            int32_t* pdata = nullptr;
            
            if (rank == mpi::master())
            {
                mpi::batches_pattern(bsizes.data(), oldsizes.at(i), nodes);
                
                mathfunc::indexes(bindexes.data(), bsizes.data(), nodes);
                
                pdata = _data.data(oldpositions.at(i));
            }
            
            // scatter data chunk

            int32_t* bdata = mblock.data(_positions.at(i));
            
            auto merror = MPI_Scatterv(pdata, bsizes.data(), bindexes.data(),
                                       MPI_INT32_T, bdata, _originalSizes.at(i),
                                       MPI_INT32_T, mpi::master(), comm);
            
            if (merror != MPI_SUCCESS) mpi::abort(merror, "scatter(CMemBlock2D)");
            
            MPI_Barrier(comm);
        }
        
        // update data
        
        _data = mblock;
    }
}

template <>
inline void
CMemBlock2D<double>::scatter(int32_t  rank,
                             int32_t  nodes,
                             MPI_Comm comm)
{
    if (ENABLE_MPI)
    {
        if (nodes == 1) return;
        
        // save old dimensions
        
        CMemBlock<int32_t> oldsizes;
        
        CMemBlock<int32_t> oldpositions;
        
        if (rank == mpi::master())
        {
            oldsizes = _originalSizes;
            
            oldpositions = _positions;
        }
        
        // update dimensions on all nodes
        
        _originalSizes.broadcast(rank, comm);
        
        for (int32_t i = 0; i < _originalSizes.size(); i++)
        {
            _originalSizes.at(i) = mpi::batch_size(_originalSizes.at(i),
                                                   rank, nodes);
        }
        
        _setDimensions();
        
        // set up scattered data storage
        
        CMemBlock<double> mblock(_nElements);
        
        mblock.zero();
        
        // allocate scaterring pattern data
        
        CMemBlock<int32_t> bsizes;
        
        CMemBlock<int32_t> bindexes;
        
        if (rank == mpi::master())
        {
            bsizes = CMemBlock<int32_t>(nodes);
            
            bindexes = CMemBlock<int32_t>(nodes);
        }
        
        // scatter data chunks
        
        for (int32_t i = 0; i < _originalSizes.size(); i++)
        {
            // set up scattering pattern for data chunk
            
            double* pdata = nullptr;
            
            if (rank == mpi::master())
            {
                mpi::batches_pattern(bsizes.data(), oldsizes.at(i), nodes);
                
                mathfunc::indexes(bindexes.data(), bsizes.data(), nodes);
                
                pdata = _data.data(oldpositions.at(i));
            }
            
            // scatter data chunk

            double* bdata = mblock.data(_positions.at(i));
            
            auto merror = MPI_Scatterv(pdata, bsizes.data(), bindexes.data(),
                                       MPI_DOUBLE, bdata, _originalSizes.at(i),
                                       MPI_DOUBLE, mpi::master(), comm);
            
            if (merror != MPI_SUCCESS) mpi::abort(merror, "scatter(CMemBlock2D)");
            
            MPI_Barrier(comm);
        }
        
        // update data
        
        _data = mblock;
    }
}

template <class T>
void
CMemBlock2D<T>::_setOriginalSizes(const int32_t nElements,
                                  const int32_t nBlocks)
{
    _originalSizes = CMemBlock<int32_t>(nBlocks);
    
    for (int32_t i = 0; i < nBlocks; i++) _originalSizes.at(i) = nElements;
}

template <class T>
void
CMemBlock2D<T>::_setDimensions()
{
    auto numblocks = _originalSizes.size();
    
    _paddedSizes = CMemBlock<int32_t>(numblocks);
    
    _positions = CMemBlock<int32_t>(numblocks);
    
    // loop over data chunks
    
    int32_t primdim = VLX_ALIGN / sizeof(T);
    
    _nElements = 0;
    
    for (int32_t i = 0; i < numblocks; i++)
    {
        // compute padded size of data chunk
        
        auto pblocks = _originalSizes.at(i) / primdim;
        
        if ((_originalSizes.at(i) % primdim) != 0) pblocks++;
        
        _paddedSizes.at(i) = pblocks * primdim;
        
        // determine start position of data chunk
        
        _positions.at(i) = _nElements;
        
        _nElements += _paddedSizes.at(i);
    }
}

template <class T>
void
CMemBlock2D<T>::_copy(const std::vector<T>& dataVector,
                      const matorder        dataOrder)
{
    // row major ordering of data in vector
    
    if (dataOrder == matorder::row_major)
    {
        int32_t pindex = 0;
        
        for (int32_t i = 0; i < _originalSizes.size(); i++)
        {
            auto currpos = _positions.at(i);
        
            for (int32_t j = 0; j < _originalSizes.at(i); j++)
            {
                _data.at(currpos + j) = dataVector[pindex];
                
                pindex++;
            }
        }
        
        return;
    }
    
    // column major ordering of data in vector
    
    if (dataOrder == matorder::col_major)
    {
        int32_t pindex = 0;
        
        auto ndim = mathfunc::max(_originalSizes.data(), _originalSizes.size());
        
        for (int32_t i = 0; i < ndim; i++)
        {
            for (int32_t j = 0; j < _originalSizes.size(); j++)
            {
                if (i < _originalSizes.at(j))
                {
                    auto currpos = _positions.at(j);
                    
                    _data.at(currpos + i) = dataVector[pindex];
                    
                    pindex++;
                }
            }
        }
    }
}

template<class U>
std::ostream&
operator<<(      std::ostream&   output, 
           const CMemBlock2D<U>& source)
{
    output << std::endl;

    output << "[CMemBlock2D (Object):" << &source << "]" << std::endl;

    output << "_nElements: " << source._nElements << std::endl;
    
    output << "_originalSizes: " << source._originalSizes << std::endl;
    
    output << "_paddedSizes: " << source._paddedSizes << std::endl;

    output << "_positions: " << source._positions << std::endl;

    output << "_data: " << source._data <<  std::endl;

    output << std::endl;

    return output;
}

#endif /* MemBlock2D_hpp */
