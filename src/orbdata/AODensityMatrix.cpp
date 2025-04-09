//
//                                   VELOXCHEM
//              ----------------------------------------------------
//                          An Electronic Structure Code
//
//  SPDX-License-Identifier: BSD-3-Clause
//
//  Copyright 2018-2025 VeloxChem developers
//
//  Redistribution and use in source and binary forms, with or without modification,
//  are permitted provided that the following conditions are met:
//
//  1. Redistributions of source code must retain the above copyright notice, this
//     list of conditions and the following disclaimer.
//  2. Redistributions in binary form must reproduce the above copyright notice,
//     this list of conditions and the following disclaimer in the documentation
//     and/or other materials provided with the distribution.
//  3. Neither the name of the copyright holder nor the names of its contributors
//     may be used to endorse or promote products derived from this software without
//     specific prior written permission.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
//  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
//  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
//  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
//  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
//  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
//  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
//  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
//  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
//  OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include "AODensityMatrix.hpp"

#include "ErrorHandler.hpp"

CAODensityMatrix::CAODensityMatrix()

    : _denType(denmat::rest)
{
}

CAODensityMatrix::CAODensityMatrix(const std::vector<CDenseMatrix>& denMatrices, const denmat denType)

    : _denMatrices(denMatrices)

    , _denType(denType)
{
    if (!isClosedShell())
    {
        errors::assertMsgCritical(denMatrices.size() % 2 == 0,
                                  "AODensityMatrix: Odd number of matrices for unrestricted or restricted open-shell density");
    }
}

CAODensityMatrix::CAODensityMatrix(const CAODensityMatrix& source)

    : _denMatrices(source._denMatrices)

    , _denType(source._denType)
{
}

CAODensityMatrix::CAODensityMatrix(CAODensityMatrix&& source) noexcept

    : _denMatrices(std::move(source._denMatrices))

    , _denType(std::move(source._denType))
{
}

CAODensityMatrix::~CAODensityMatrix()
{
}

CAODensityMatrix&
CAODensityMatrix::operator=(const CAODensityMatrix& source)
{
    if (this == &source) return *this;

    _denMatrices = source._denMatrices;

    _denType = source._denType;

    return *this;
}

CAODensityMatrix&
CAODensityMatrix::operator=(CAODensityMatrix&& source) noexcept
{
    if (this == &source) return *this;

    _denMatrices = std::move(source._denMatrices);

    _denType = std::move(source._denType);

    return *this;
}

bool
CAODensityMatrix::operator==(const CAODensityMatrix& other) const
{
    if (_denType != other._denType) return false;

    if (_denMatrices.size() != other._denMatrices.size()) return false;

    for (size_t i = 0; i < _denMatrices.size(); i++)
    {
        if (_denMatrices[i] != other._denMatrices[i]) return false;
    }

    return true;
}

bool
CAODensityMatrix::operator!=(const CAODensityMatrix& other) const
{
    return !(*this == other);
}

bool
CAODensityMatrix::isClosedShell() const
{
    return (_denType == denmat::rest);
}

int
CAODensityMatrix::getNumberOfDensityMatrices() const
{
    if (isClosedShell())
    {
        return static_cast<int>(_denMatrices.size());
    }
    else
    {
        return static_cast<int>(_denMatrices.size()) / 2;
    }
}

denmat
CAODensityMatrix::getDensityType() const
{
    return _denType;
}

int
CAODensityMatrix::getNumberOfRows(const int iDensityMatrix) const
{
    if (iDensityMatrix < getNumberOfDensityMatrices())
    {
        if (isClosedShell())
        {
            return _denMatrices[iDensityMatrix].getNumberOfRows();
        }
        else
        {
            return _denMatrices[2 * iDensityMatrix].getNumberOfRows();
        }
    }

    return 0;
}

int
CAODensityMatrix::getNumberOfColumns(const int iDensityMatrix) const
{
    if (iDensityMatrix < getNumberOfDensityMatrices())
    {
        if (isClosedShell())
        {
            return _denMatrices[iDensityMatrix].getNumberOfColumns();
        }
        else
        {
            return _denMatrices[2 * iDensityMatrix].getNumberOfColumns();
        }
    }

    return 0;
}

int
CAODensityMatrix::getNumberOfElements(const int iDensityMatrix) const
{
    if (iDensityMatrix < getNumberOfDensityMatrices())
    {
        if (isClosedShell())
        {
            return _denMatrices[iDensityMatrix].getNumberOfElements();
        }
        else
        {
            return _denMatrices[2 * iDensityMatrix].getNumberOfElements();
        }
    }

    return 0;
}

const double*
CAODensityMatrix::alphaDensity(const int iDensityMatrix) const
{
    if (iDensityMatrix < getNumberOfDensityMatrices())
    {
        if (isClosedShell())
        {
            return _denMatrices[iDensityMatrix].values();
        }
        else
        {
            return _denMatrices[2 * iDensityMatrix].values();
        }
    }

    return nullptr;
}

const double*
CAODensityMatrix::betaDensity(const int iDensityMatrix) const
{
    if (iDensityMatrix < getNumberOfDensityMatrices())
    {
        if (isClosedShell())
        {
            return _denMatrices[iDensityMatrix].values();
        }
        else
        {
            return _denMatrices[2 * iDensityMatrix + 1].values();
        }
    }

    return nullptr;
}

const CDenseMatrix&
CAODensityMatrix::getReferenceToDensity(const int iDensityMatrix) const
{
    return _denMatrices[iDensityMatrix];
}
