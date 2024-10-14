//
//                              VELOXCHEM
//         ----------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2024 by VeloxChem developers. All rights reserved.
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

#ifndef MemBlock2D_hpp
#define MemBlock2D_hpp

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <sstream>
#include <string>

#include "MemBlock.hpp"

/**
 Templated class CMemBlock2D manages 2D memory block allocation, manipulation,
 and deallocation.
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
    CMemBlock<int> _originalSizes;

    /**
     The padded sizes of data chunks in memory block.
     */
    CMemBlock<int> _paddedSizes;

    /**
     The start positions of data chunks in memory block.
     */
    CMemBlock<int> _positions;

    /**
     The number of elements in memory block.
     */
    int _nElements;

    /**
     Sets uniform original size for all data chunks.

     @param nElements the number of elements in data chunk.
     @param nBlocks the number of data chunks.
     */
    void _setOriginalSizes(const int nElements, const int nBlocks);

    /**
     Sets dimensions, padding for all data chunks.
     */
    void _setDimensions();

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
    CMemBlock2D(const int nElements, const int nBlocks);

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
    CMemBlock2D<T> slice(const int iPosition, const int nElements) const;

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
    T* data(const int iBlock);

    /**
     Gets constant pointer to specific data chunk in 2D memory block.

     @param iBlock the index of data chunk.
     @return the pointer to data chunk.
     */
    const T* data(const int iBlock) const;

    /**
     Gets pointer to specific element of selected data chunk in 2D memory block.

     @param iBlock the index of data chunk.
     @param iElement the index of element in memory block.
     @return the pointer to data chunk.
     */
    T* data(const int iBlock, const int iElement);

    /**
     Gets constant pointer to specific elment of selected data chunk in 2D
     memory block.

     @param iBlock the index of data chunk.
     @param iElement the index of element in memory block.
     @return the pointer to data chunk.
     */
    const T* data(const int iBlock, const int iElement) const;

    /**
     Gets number of elements in specific data chunk.

     @param iBlock the index of data chunk.
     @return the number of elements in data chunk.
     */
    int size(const int iBlock) const;

    /**
     Gets number of elements in specific pitched data chunk.

     @param iBlock the index of data chunk.
     @return the number of elements in pitched data chunk.
     */
    int pitched_size(const int iBlock) const;

    /**
     Gets number of data chunk.

     @return the number of data chunk.
     */
    int blocks() const;
};

template <class T>
CMemBlock2D<T>::CMemBlock2D()

    : _nElements(0)
{
}

template <class T>
CMemBlock2D<T>::CMemBlock2D(const int nElements, const int nBlocks)

    : _nElements(0)
{
    _setOriginalSizes(nElements, nBlocks);

    _setDimensions();

    _data = CMemBlock<T>(_nElements);

    _data.zero();
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

template <class T>
CMemBlock2D<T>::~CMemBlock2D()
{
}

template <class T>
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

template <class T>
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

template <class T>
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

template <class T>
bool
CMemBlock2D<T>::operator!=(const CMemBlock2D<T>& other) const
{
    return !(*this == other);
}

template <class T>
CMemBlock2D<T>
CMemBlock2D<T>::slice(const int iPosition, const int nElements) const
{
    CMemBlock2D<T> mblock(nElements, blocks());

    for (int i = 0; i < blocks(); i++)
    {
        // set up pointers to data chunks

        auto idata = data(i, iPosition);

        auto odata = mblock.data(i);

        // copy elements of data chunks

        for (int j = 0; j < nElements; j++)
            odata[j] = idata[j];
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

template <class T>
T*
CMemBlock2D<T>::data(const int iBlock)
{
    if (_originalSizes.size() > 0)
    {
        if (iBlock < 0) return nullptr;

        if (iBlock >= blocks()) return nullptr;

        return _data.data(_positions.at(iBlock));
    }

    return nullptr;
}

template <class T>
const T*
CMemBlock2D<T>::data(const int iBlock) const
{
    if (_originalSizes.size() > 0)
    {
        if (iBlock < 0) return nullptr;

        if (iBlock >= blocks()) return nullptr;

        return _data.data(_positions.at(iBlock));
    }

    return nullptr;
}

template <class T>
T*
CMemBlock2D<T>::data(const int iBlock, const int iElement)
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
const T*
CMemBlock2D<T>::data(const int iBlock, const int iElement) const
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
int
CMemBlock2D<T>::size(const int iBlock) const
{
    if (iBlock < _originalSizes.size()) return _originalSizes.at(iBlock);

    return 0;
}

template <class T>
int
CMemBlock2D<T>::pitched_size(const int iBlock) const
{
    if (iBlock < _paddedSizes.size()) return _paddedSizes.at(iBlock);

    return 0;
}

template <class T>
int
CMemBlock2D<T>::blocks() const
{
    return _originalSizes.size();
}

template <class T>
void
CMemBlock2D<T>::_setOriginalSizes(const int nElements, const int nBlocks)
{
    _originalSizes = CMemBlock<int>(nBlocks);

    for (int i = 0; i < nBlocks; i++)
        _originalSizes.at(i) = nElements;
}

template <class T>
void
CMemBlock2D<T>::_setDimensions()
{
    auto numblocks = _originalSizes.size();

    _paddedSizes = CMemBlock<int>(numblocks);

    _positions = CMemBlock<int>(numblocks);

    // loop over data chunks

    int primdim = VLX_ALIGN / sizeof(T);

    _nElements = 0;

    for (int i = 0; i < numblocks; i++)
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

#endif /* MemBlock2D_hpp */
