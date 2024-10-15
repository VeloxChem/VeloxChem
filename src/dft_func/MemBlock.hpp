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

#ifndef MemBlock_hpp
#define MemBlock_hpp

#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <sstream>
#include <string>
#include <vector>

#include "MathFunc.hpp"
#include "MemAlloc.hpp"

/**
 Templated class CMemBlock manages contiguous memory block allocation,
 manipulation, and deallocation.
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
    int _nElements{0};

    /**
     Allocates memory block.
     */
    void _allocate();

    /**
     Copies aligned data to memory block from data source.

     @param source the pointer to data source.
     */
    void _copy_aligned(const T* source);

   public:
    /**
     Creates an empty memory block object.
     */
    CMemBlock();

    /**
     Creates a memory block object.

     @param nElements - the number of elements.
     */
    explicit CMemBlock(const int nElements);

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
    T* data(const int iElement);

    /**
     Gets constant pointer to specific element in memory block.

     @param iElement the index of element in memory block.
     @return the constant pointer to element in memory block.
     */
    const T* data(const int iElement) const;

    /**
     Gets number of elements in memory block.

     @return the number of elements.
     */
    int size() const;

    /**
     Sets memory block elements to zero.
     */
    void zero();

    /**
     Gets reference to specific element of memory block.

     @param iElement the index of element in memory block.
     @return the reference to element.
     */
    T& at(const int iElement);

    /**
     Gets constant reference to specific element of memory block.

     @param iElement the index of element in memory block.
     @return the reference to element.
     */
    const T& at(const int iElement) const;
};

template <class T>
CMemBlock<T>::CMemBlock()

    : _data(nullptr)

    , _nElements(0)
{
}

template <class T>
CMemBlock<T>::CMemBlock(const int nElements)

    : _data(nullptr)

    , _nElements(nElements)
{
    _allocate();

    zero();
}

template <class T>
CMemBlock<T>::CMemBlock(const CMemBlock<T>& source)

    : _data(nullptr)

    , _nElements(source._nElements)
{
    _allocate();

    _copy_aligned(source.data());
}

template <class T>
CMemBlock<T>::CMemBlock(CMemBlock<T>&& source) noexcept

    : _data(nullptr)

    , _nElements(std::move(source._nElements))
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

    for (int i = 0; i < _nElements; i++)
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
CMemBlock<T>::data(const int iElement)
{
    return &(_data[iElement]);
}

template <class T>
const T*
CMemBlock<T>::data(const int iElement) const
{
    return &(_data[iElement]);
}

template <class T>
int
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

#pragma omp simd aligned(pdata : VLX_ALIGN)
    for (int i = 0; i < _nElements; i++)
    {
        pdata[i] = t0value;
    }
}

template <class T>
inline T&
CMemBlock<T>::at(const int iElement)
{
    return _data[iElement];
}

template <class T>
inline const T&
CMemBlock<T>::at(const int iElement) const
{
    return _data[iElement];
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
CMemBlock<T>::_copy_aligned(const T* source)
{
    auto pdata = _data;

#pragma omp simd aligned(pdata, source : VLX_ALIGN)
    for (int i = 0; i < _nElements; i++)
    {
        pdata[i] = source[i];
    }
}

#endif /* MemBlock_hpp */
