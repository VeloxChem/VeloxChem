//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#include "SparseVector.hpp"

CSparseVector::CSparseVector()

    : _nLength(0)

    , _nMaxElements(0)

    , _nElements(0)

    , _threshold(1.0e-15)
{
}

CSparseVector::CSparseVector(const std::vector<double>& values, const std::vector<int32_t>& indexes, const int32_t nLength, const double threshold)

    : _nLength(nLength)

    , _values(CMemBlock<double>(values))

    , _indexes(CMemBlock<int32_t>(indexes))

    , _nMaxElements(_values.size())

    , _nElements(_values.size())

    , _threshold(threshold)
{
}

CSparseVector::CSparseVector(const int32_t nLength, const double threshold)

    : _nLength(nLength)

    , _nElements(0)

    , _threshold(threshold)
{
    _nMaxElements = _setMaxNumberOfElements();

    _values = CMemBlock<double>(_nMaxElements);

    _indexes = CMemBlock<int32_t>(_nMaxElements);
}

CSparseVector::CSparseVector(const CSparseVector& source)

    : _nLength(source._nLength)

    , _values(source._values)

    , _indexes(source._indexes)

    , _nMaxElements(source._nMaxElements)

    , _nElements(source._nElements)

    , _threshold(source._threshold)
{
}

CSparseVector::CSparseVector(CSparseVector&& source) noexcept

    : _nLength(std::move(source._nLength))

    , _values(std::move(source._values))

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

    _nLength = source._nLength;

    _values = source._values;

    _indexes = source._indexes;

    _nMaxElements = source._nMaxElements;

    _nElements = source._nElements;

    _threshold = source._threshold;

    return *this;
}

CSparseVector&
CSparseVector::operator=(CSparseVector&& source) noexcept
{
    if (this == &source) return *this;

    _nLength = std::move(source._nLength);

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
    if (_nLength != other._nLength) return false;

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

const double*
CSparseVector::values() const
{
    return _values.data();
}

double*
CSparseVector::values()
{
    return _values.data();
}

const int32_t*
CSparseVector::indexes() const
{
    return _indexes.data();
}

double
CSparseVector::getThreshold() const
{
    return _threshold;
}

double
CSparseVector::getSparsity() const
{
    if (_nLength > 0)
    {
        return static_cast<double>(_nElements) / static_cast<double>(_nLength);
    }

    return 0.0;
}

int32_t
CSparseVector::_setMaxNumberOfElements() const
{
    if (_nLength > 100000) return _nLength / 5;

    if (_nLength > 20000) return _nLength / 3;

    return _nLength / 2;
}

std::ostream&
operator<<(std::ostream& output, const CSparseVector& source)
{
    output << std::endl;

    output << "[CSparseVector (Object):" << &source << "]" << std::endl;

    output << "_nLength: " << source._nLength << std::endl;

    output << "_values: " << source._values << std::endl;

    output << "_indexes: " << source._indexes << std::endl;

    output << "_nMaxElements: " << source._nMaxElements << std::endl;

    output << "_nElements: " << source._nElements << std::endl;

    output << "_threshold: " << source._threshold << std::endl;

    return output;
}
