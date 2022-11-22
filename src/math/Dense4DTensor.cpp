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

#include "Dense4DTensor.hpp"

#include <cmath>
#include <sstream>
#include <utility>

#include <mpi.h>

#include "StringFormat.hpp"

CDense4DTensor::CDense4DTensor()

    : _iIndex(0)

    , _jIndex(0)

    , _kIndex(0)

    , _lIndex(0)
{
}

CDense4DTensor::CDense4DTensor(const std::vector<double>& values,
                           const int32_t              iIndex,
                           const int32_t              jIndex,
                           const int32_t              kIndex,
                           const int32_t              lIndex)

    : _iIndex(iIndex)

    , _jIndex(jIndex)

    , _kIndex(kIndex)

    , _lIndex(lIndex)

    , _values(CMemBlock<double>(values))
{
}

CDense4DTensor::CDense4DTensor(const int32_t              iIndex,
                           const int32_t              jIndex,
                           const int32_t              kIndex,
                           const int32_t              lIndex)

    : _iIndex(iIndex)

    , _jIndex(jIndex)

    , _kIndex(kIndex)

    , _lIndex(lIndex)

    , _values(CMemBlock<double>(iIndex * jIndex * kIndex * lIndex))
{
}

CDense4DTensor::CDense4DTensor(const int32_t nRows)

    : _iIndex(nRows)

    , _jIndex(nRows)

    , _kIndex(nRows)

    , _lIndex(nRows)

    , _values(CMemBlock<double>(nRows * nRows * nRows * nRows))
{
}

CDense4DTensor::CDense4DTensor(const CDense4DTensor& source)

    : _iIndex(source._iIndex)

    , _jIndex(source._jIndex)

    , _kIndex(source._lIndex)

    , _lIndex(source._lIndex)

    , _values(source._values)
{
}

CDense4DTensor::CDense4DTensor(CDense4DTensor&& source) noexcept

    : _iIndex(std::move(source._iIndex))

    , _jIndex(std::move(source._jIndex))

    , _kIndex(std::move(source._kIndex))

    , _lIndex(std::move(source._lIndex))

    , _values(std::move(source._values))
{
}

CDense4DTensor::~CDense4DTensor()
{
}

CDense4DTensor&
CDense4DTensor::operator=(const CDense4DTensor& source)
{
    if (this == &source) return *this;

    _iIndex = source._iIndex;

    _jIndex = source._jIndex;

    _lIndex = source._kIndex;

    _lIndex = source._lIndex;

    _values = source._values;

    return *this;
}

CDense4DTensor&
CDense4DTensor::operator=(CDense4DTensor&& source) noexcept
{
    if (this == &source) return *this;

    _iIndex = std::move(source._iIndex);

    _jIndex = std::move(source._jIndex);

    _kIndex = std::move(source._kIndex);

    _lIndex = std::move(source._lIndex);

    _values = std::move(source._values);

    return *this;
}

bool
CDense4DTensor::operator==(const CDense4DTensor& other) const
{
    if (_iIndex != other._iIndex) return false;

    if (_jIndex != other._jIndex) return false;

    if (_kIndex != other._kIndex) return false;

    if (_lIndex != other._lIndex) return false;

    if (_values != other._values) return false;

    return true;
}

bool
CDense4DTensor::operator!=(const CDense4DTensor& other) const
{
    return !(*this == other);
}

void
CDense4DTensor::zero()
{
    mathfunc::zero(_values.data(), _iIndex * _jIndex * _kIndex * _lIndex);
}

int32_t
CDense4DTensor::getiIndex() const
{
    return _iIndex;
}

int32_t
CDense4DTensor::getjIndex() const
{
    return _jIndex;
}

int32_t
CDense4DTensor::getkIndex() const
{
    return _kIndex;
}

int32_t
CDense4DTensor::getlIndex() const
{
    return _lIndex;
}

int32_t
CDense4DTensor::getNumberOfElements() const
{
    return _iIndex * _jIndex * _kIndex * _lIndex;
}

const double*
CDense4DTensor::values() const
{
    return _values.data();
}

double*
CDense4DTensor::values()
{
    return _values.data();
}


std::string
CDense4DTensor::getString() const
{
    std::stringstream sst("");

    auto vals = _values.data();

    sst << "[Dimension " << _iIndex << " x " << _jIndex << " x " << _kIndex << " x " << _lIndex  <<"]\n";

    for (int32_t i = 0; i < _iIndex; i++)
    {
        for (int32_t j = 0; j < _jIndex; j++)
        {
            sst << i <<" " <<j <<"\n";
            int32_t ij = (i * _jIndex +j) * _kIndex * _lIndex;
            for (int32_t k = 0; k < _kIndex; k++)
            {
                for (int32_t l = 0; l < _lIndex; l++)
                {
                    sst << fstr::to_string(vals[ij + k * _lIndex + l], 8, 15, fmt::right);
                }
                sst << "\n";
            }
            sst << "\n";
        }

        sst << "\n";
    }

    return sst.str();
}

void
CDense4DTensor::broadcast(int32_t rank, MPI_Comm comm)
{
    if constexpr (ENABLE_MPI)
    {
        mpi::bcast(_iIndex, comm);

        mpi::bcast(_jIndex, comm);

        mpi::bcast(_kIndex, comm);

        mpi::bcast(_lIndex, comm);

        _values.broadcast(rank, comm);
    }
}

std::ostream&
operator<<(std::ostream& output, const CDense4DTensor& source)
{
    output << std::endl;

    output << "[CDense4DTensor (Object):" << &source << "]" << std::endl;

    output << "_iIndex: " << source._iIndex << std::endl;

    output << "_jIndex: " << source._jIndex << std::endl;

    output << "_kIndex: " << source._kIndex << std::endl;

    output << "_lIndex: " << source._lIndex << std::endl;

    output << "_values: " << source._values << std::endl;

    return output;
}
