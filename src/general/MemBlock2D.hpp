//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

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
    void _setOriginalSizes(const int32_t nElements, const int32_t nBlocks);
    
    /**
     Sets dimensions, padding for all data chunks.
     */
    void _setDimensions();
    
    /**
     Copies data elements from vector to memory block.

     @param dataVector the vector with data elements.
     */
    void _copy(const std::vector<T>& dataVector);
    
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
    CMemBlock2D(const int32_t nElements, const int32_t nBlocks);
    
    /**
     Creates an 2D memory block object.
     
     @param dimVector the vector of data chunk sizes.
     */
    CMemBlock2D(const std::vector<int32_t>& dimVector);

    /**
     Creates an 2D memory block object.
     
     @param dataVector the vector with data elements.
     @param nElements the number of elements in data chunk.
     @param nBlocks the number of data chunks.
     */
    CMemBlock2D(const std::vector<T>& dataVector, const int32_t nElements,
                const int32_t nBlocks);
    
    /**
     Creates an 2D memory block object.
     
     @param dataVector the vector with data elements.
     @param dimVector the vector of data chunk sizes.
     */
    CMemBlock2D(const std::vector<T>& dataVector,
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
     Sets all elements of contiguous memory block to zero.
     */
    void zero();
    
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
     Converts 2D memory block object to text output and insert it into output
     text stream.
     
     @param output the output text stream.
     @param source the 2D memory block object.
     */
    template <class U>
    friend std::ostream& operator<<(std::ostream& output,
                                    const CMemBlock2D<U>& source);
};

template <class T>
CMemBlock2D<T>::CMemBlock2D()

    : _nElements(0)
{

}

template <class T>
CMemBlock2D<T>::CMemBlock2D(const int32_t nElements, const int32_t nBlocks)

    : _nElements(0)
{
    _setOriginalSizes(nElements, nBlocks);
    
    _setDimensions();
    
    _data = CMemBlock<T>(_nElements);
}

template<class T>
CMemBlock2D<T>::CMemBlock2D(const std::vector<int32_t>& dimVector)
{
    _originalSizes = CMemBlock<int32_t>(dimVector);
    
    _setDimensions();
    
    _data = CMemBlock<T>(_nElements);
}

template<class T>
CMemBlock2D<T>::CMemBlock2D(const std::vector<T>& dataVector,
                            const int32_t nElements, const int32_t nBlocks)

    : _nElements(0)
{
    _setOriginalSizes(nElements, nBlocks);
    
    _setDimensions();
    
    _data = CMemBlock<T>(_nElements);
    
    _data.zero();
    
    _copy(dataVector);
 }

template<class T>
CMemBlock2D<T>::CMemBlock2D(const std::vector<T>& dataVector,
                            const std::vector<int32_t>& dimVector)
{
    _originalSizes = CMemBlock<int32_t>(dimVector);
    
    _setDimensions();
    
    _data = CMemBlock<T>(_nElements);
    
    _data.zero();
    
    _copy(dataVector);
}

template <class T>
CMemBlock2D<T>::CMemBlock2D(const CMemBlock2D<T>& source)

    : _nElements(source._nElements)

    , _positions(source._positions)

    , _paddedSizes(source._paddedSizes)

    , _originalSizes(source._originalSizes)

    , _data(source._data)
{
}

template <class T>
CMemBlock2D<T>::CMemBlock2D(CMemBlock2D<T>&& source) noexcept

    : _nElements(std::move(source._nElements))

    , _positions(std::move(source._positions))

    , _paddedSizes(std::move(source._paddedSizes))

    , _originalSizes(std::move(source._originalSizes))

    , _data(std::move(source._data))
{

}

template<class T>
CMemBlock2D<T>::~CMemBlock2D()
{
    
}

template<class T>
CMemBlock2D<T>& CMemBlock2D<T>::operator=(const CMemBlock2D<T>& source)
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
CMemBlock2D<T>& CMemBlock2D<T>::operator=(CMemBlock2D<T>&& source) noexcept
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
bool CMemBlock2D<T>::operator==(const CMemBlock2D<T>& other) const
{
    if (_nElements != other._nElements) return false;
    
    if (_positions != other._positions) return false;
    
    if (_paddedSizes != other._paddedSizes) return false;
    
    if (_originalSizes != other._originalSizes) return false;
    
    if (_data != other._data) return false;

    return true;
}

template<class T>
bool CMemBlock2D<T>::operator!=(const CMemBlock2D<T>& other) const
{
    return !(*this == other);
}

template <class T>
void CMemBlock2D<T>::zero()
{
    _data.zero(); 
}

template<class T>
T* CMemBlock2D<T>::data(const int32_t iBlock)
{
    auto pdata = _data.data();
    
    return &(pdata[_positions.at(iBlock)]);
}

template<class T>
const T* CMemBlock2D<T>::data(const int32_t iBlock) const 
{
    auto pdata = _data.data();
    
    return &(pdata[_positions.at(iBlock)]);
}

template <class T>
int32_t CMemBlock2D<T>::size(const int32_t iBlock) const
{
    return _originalSizes.at(iBlock);
}

template <class T>
int32_t CMemBlock2D<T>::blocks() const
{
    return _originalSizes.size();
}

template <class T>
void CMemBlock2D<T>::_setOriginalSizes(const int32_t nElements,
                                       const int32_t nBlocks)
{
    _originalSizes = CMemBlock<int32_t>(nBlocks);
    
    for (int32_t i = 0; i < nBlocks; i++) _originalSizes.at(i) = nElements;
}

template <class T>
void CMemBlock2D<T>::_setDimensions()
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
void CMemBlock2D<T>::_copy(const std::vector<T>& dataVector)
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
