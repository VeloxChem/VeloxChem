//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include "SparseVector.hpp"

CSparseVector::CSparseVector()

    : _nMaxElements(0)

    , _nElements(0)

    , _threshold(1.0e-13)
{
    
}

CSparseVector::CSparseVector(const std::vector<double>&  values,
                             const std::vector<int32_t>& indexes,
                             const double                threshold)

    : _values(CMemBlock<double>(values))

    , _indexes(CMemBlock<int32_t>(indexes))

    , _nElements(_values.size())

    , _nMaxElements(_values.size())

    , _threshold(threshold)
{
    
}

CSparseVector::CSparseVector(const int32_t nElements,
                             const double  threshold)

    : _values(CMemBlock<double>(nElements))

    , _indexes(CMemBlock<int32_t>(nElements))

    , _nMaxElements(nElements)

    , _nElements(0)

    , _threshold(threshold)
{
    
}

CSparseVector::CSparseVector(const CSparseVector& source)

    : _values(source._values)

    , _indexes(source._indexes)

    , _nMaxElements(source._nMaxElements)

    , _nElements(source._nElements)

    , _threshold(source._threshold)
{
    
}

CSparseVector::CSparseVector(CSparseVector&& source) noexcept

    : _values(std::move(source._values))

    , _indexes(std::move(source._indexes))

    , _nMaxElements(std::move(source._nMaxElements))

    , _nElements(std::move(source._nElements))

    , _threshold(std::move(source._threshold))
{
    
}

CSparseVector::~CSparseVector()
{
    
}

CSparseVector&
CSparseVector::operator=(const CSparseVector& source)
{
    if (this == &source) return *this;
    
    _values = source._values;
    
    _indexes = source._indexes;
    
    _nMaxElements = source._nMaxElements;
    
    _nElements = source._nElements;
    
    _threshold= source._threshold;
    
    return *this;
}

CSparseVector&
CSparseVector::operator=(CSparseVector&& source) noexcept
{
    if (this == &source) return *this;
    
    _values = std::move(source._values);
    
    _indexes = std::move(source._indexes);
    
    _nMaxElements = std::move(source._nMaxElements);
    
    _nElements = std::move(source._nElements);
    
    _threshold = std::move(source._threshold);
    
    return *this;
}

bool
CSparseVector::operator==(const CSparseVector& other) const
{
    // NOTE: max size of buffer is not uniquely defined and depends on
    // constructor used to initialize sparse vector.
    
    if (_nElements != other._nElements) return false;
    
    // compare relevant elements
    
    for (int32_t i = 0; i < _nElements; i++)
    {
        if (std::fabs(_values.at(i) - other._values.at(i)) > 1.0e-13) return false;
        
        if (_indexes.at(i) != other._indexes.at(i)) return false;
    }
    
    if (std::fabs(_threshold - other._threshold) > 1.0e-13) return false;
    
    return true;
}

bool
CSparseVector::operator!=(const CSparseVector& other) const
{
    return !(*this == other);
}

std::ostream&
operator<<(      std::ostream&  output,
           const CSparseVector& source)
{
    output << std::endl;
    
    output << "[CSparseVector (Object):" << &source << "]" << std::endl;
    
    output << "_values: " << source._values <<  std::endl;
    
    output << "_indexes: " << source._indexes <<  std::endl;
    
    output << "_nMaxElements: " << source._nMaxElements  <<  std::endl;
    
    output << "_nElements: " << source._nElements  <<  std::endl;
    
    output << "_threshold: " << source._threshold  <<  std::endl;
    
    return output;
}
