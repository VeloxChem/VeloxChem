//
//                           VELOXCHEM 1.0-RC
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

#include "RecursionMap.hpp"

CRecursionMap::CRecursionMap()
{
    _angularForm = recblock::cc;

    _repeatUnits = 0;
}

CRecursionMap::CRecursionMap(const recblock angularForm, const int32_t repeatUnits)
{
    _angularForm = angularForm;

    _repeatUnits = repeatUnits;
}

CRecursionMap::CRecursionMap(const std::vector<CRecursionTerm>& recursionTerms,
                             const std::vector<int32_t>&        recursionIndexes,
                             const recblock                     angularForm,
                             const int32_t                      repeatUnits)

    : _recursionTerms(recursionTerms)

    , _recursionIndexes(recursionIndexes)

    , _angularForm(angularForm)

    , _repeatUnits(repeatUnits)
{
}

CRecursionMap::CRecursionMap(const CRecursionMap& source)

    : _recursionTerms(source._recursionTerms)

    , _recursionIndexes(source._recursionIndexes)

    , _angularForm(source._angularForm)

    , _repeatUnits(source._repeatUnits)
{
}

CRecursionMap::CRecursionMap(CRecursionMap&& source) noexcept

    : _recursionTerms(std::move(source._recursionTerms))

    , _recursionIndexes(std::move(source._recursionIndexes))

    , _angularForm(std::move(source._angularForm))

    , _repeatUnits(std::move(source._repeatUnits))
{
}

CRecursionMap::~CRecursionMap()
{
}

CRecursionMap&
CRecursionMap::operator=(const CRecursionMap& source)
{
    if (this == &source) return *this;

    _recursionTerms = source._recursionTerms;

    _recursionIndexes = source._recursionIndexes;

    _angularForm = source._angularForm;

    _repeatUnits = source._repeatUnits;

    return *this;
}

CRecursionMap&
CRecursionMap::operator=(CRecursionMap&& source) noexcept
{
    if (this == &source) return *this;

    _recursionTerms = std::move(source._recursionTerms);

    _recursionIndexes = std::move(source._recursionIndexes);

    _angularForm = std::move(source._angularForm);

    _repeatUnits = std::move(source._repeatUnits);

    return *this;
}

bool
CRecursionMap::operator==(const CRecursionMap& other) const
{
    if (_recursionTerms.size() != other._recursionTerms.size()) return false;

    for (size_t i = 0; i < _recursionTerms.size(); i++)
    {
        if (_recursionTerms[i] != other._recursionTerms[i]) return false;
    }

    if (_recursionIndexes.size() != other._recursionIndexes.size()) return false;

    for (size_t i = 0; i < _recursionIndexes.size(); i++)
    {
        if (_recursionIndexes[i] != other._recursionIndexes[i]) return false;
    }

    if (_angularForm != other._angularForm) return false;

    if (_repeatUnits != other._repeatUnits) return false;

    return true;
}

bool
CRecursionMap::operator!=(const CRecursionMap& other) const
{
    return !(*this == other);
}

void
CRecursionMap::add(const CRecursionTerm& recursionTerm)
{
    if (recursionTerm.isValid())
    {
        if (!find(recursionTerm))
        {
            auto ncomps = getNumberOfComponents();

            _recursionTerms.push_back(recursionTerm);

            _recursionIndexes.push_back(ncomps);
        }
    }
}

void
CRecursionMap::append(const CRecursionMap& source)
{
    if (_angularForm == source._angularForm)
    {
        for (size_t i = 0; i < source._recursionTerms.size(); i++)
        {
            add(source._recursionTerms[i]);
        }
    }
}

void
CRecursionMap::append(const std::vector<CRecursionTerm>& recursionTerms)
{
    for (size_t i = 0; i < recursionTerms.size(); i++)
    {
        add(recursionTerms[i]);
    }
}

int32_t
CRecursionMap::index(const CRecursionTerm& recursionTerm) const
{
    for (size_t i = 0; i < _recursionTerms.size(); i++)
    {
        if (recursionTerm == _recursionTerms[i])
        {
            return static_cast<int32_t>(i);
        }
    }
    
    return -1;
}

int32_t
CRecursionMap::getNumberOfComponents() const
{
    int32_t ncomps = 0;

    for (size_t i = 0; i < _recursionTerms.size(); i++)
    {
        ncomps += _recursionTerms[i].getNumberOfComponents(_angularForm);
    }

    return _repeatUnits * ncomps;
}

int32_t
CRecursionMap::getNumberOfTerms() const
{
    return static_cast<int32_t>(_recursionTerms.size());
}

CRecursionTerm
CRecursionMap::getTerm(const int32_t iRecursionTerm) const
{
    if (iRecursionTerm < getNumberOfTerms())
    {
        return _recursionTerms[iRecursionTerm];
    }

    return CRecursionTerm();
}

int32_t
CRecursionMap::getIndexOfTerm(const CRecursionTerm& recursionTerm) const
{
    for (size_t i = 0; i < _recursionTerms.size(); i++)
    {
        if (recursionTerm == _recursionTerms[i])
        {
            return _recursionIndexes[i];
        }
    }

    return -1;
}

int32_t
CRecursionMap::getMaxOrder(const std::string&  label,
                           const CFourIndexes& braAngularMomentum,
                           const CFourIndexes& ketAngularMomentum,
                           const int32_t       braCenters,
                           const int32_t       ketCenters) const
{
    int32_t mord = -1;

    for (size_t i = 0; i < _recursionTerms.size(); i++)
    {
        if (_recursionTerms[i].isIntegral(label, braAngularMomentum, ketAngularMomentum, braCenters, ketCenters))
        {
            auto cord = _recursionTerms[i].getOrder();

            if (cord > mord) mord = cord;
        }
    }

    return mord;
}

bool
CRecursionMap::find(const CRecursionTerm& recursionTerm) const
{
    for (size_t i = 0; i < _recursionTerms.size(); i++)
    {
        if (recursionTerm == _recursionTerms[i]) return true;
    }

    return false;
}

CMemBlock2D<double>*
CRecursionMap::createBuffer(const int32_t length) const
{
    CMemBlock2D<double>* buffer = new CMemBlock2D<double>[_recursionTerms.size()];
    
    for (size_t i = 0; i < _recursionTerms.size(); i++)
    {
        buffer[i] = CMemBlock2D<double>(length, _recursionTerms[i].getNumberOfComponents(_angularForm) * _repeatUnits);
    }
    
    return buffer;
}

void
CRecursionMap::destroyBuffer(CMemBlock2D<double>* buffer) const
{
    delete [] buffer;
}

std::ostream&
operator<<(std::ostream& output, const CRecursionMap& source)
{
    output << std::endl;

    output << "[CRecursionMap (Object):" << &source << "]" << std::endl;

    output << "_recursionTerms: " << std::endl;

    for (size_t i = 0; i < source._recursionTerms.size(); i++)
    {
        output << "_recursionTerms[" << i << "]: " << std::endl;

        output << source._recursionTerms[i] << std::endl;
    }

    output << "_recursionIndexes: " << std::endl;

    for (size_t i = 0; i < source._recursionIndexes.size(); i++)
    {
        output << "_recursionIndexes[" << i << "]: " << std::endl;

        output << source._recursionIndexes[i] << std::endl;
    }

    output << "_angularForm: " << to_string(source._angularForm) << std::endl;

    return output;
}
