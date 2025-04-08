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

#include "AOKohnShamMatrix.hpp"

#include <cmath>
#include <string>

#include "ErrorHandler.hpp"
#include "StringFormat.hpp"

CAOKohnShamMatrix::CAOKohnShamMatrix()

    : _xcRestricted(true)

    , _xcElectrons(0.0)

    , _xcEnergy(0.0)
{
    _xcMatrices = std::vector<CDenseMatrix>({CDenseMatrix()});
}

CAOKohnShamMatrix::CAOKohnShamMatrix(const int nRows, const int nColumns, const bool xcRestricted)

    : _xcRestricted(xcRestricted)

    , _xcElectrons(0.0)

    , _xcEnergy(0.0)
{
    _xcMatrices = std::vector<CDenseMatrix>({CDenseMatrix(nRows, nColumns)});

    if (!_xcRestricted) _xcMatrices.push_back(CDenseMatrix(nRows, nColumns));
}

CAOKohnShamMatrix::CAOKohnShamMatrix(const int nRows, const int nColumns, const std::string& flag)

    : _xcElectrons(0.0)

    , _xcEnergy(0.0)
{
    _xcRestricted = ((format::lower_case(flag) == std::string("closedshell")) ||
                     (format::lower_case(flag) == std::string("closed-shell")) ||
                     (format::lower_case(flag) == std::string("closed_shell")));

    _xcMatrices = std::vector<CDenseMatrix>({CDenseMatrix(nRows, nColumns)});

    if (!_xcRestricted) _xcMatrices.push_back(CDenseMatrix(nRows, nColumns));
}

CAOKohnShamMatrix::CAOKohnShamMatrix(const CAOKohnShamMatrix& source)

    : _xcMatrices(source._xcMatrices)

    , _xcRestricted(source._xcRestricted)

    , _xcElectrons(source._xcElectrons)

    , _xcEnergy(source._xcEnergy)
{
}

CAOKohnShamMatrix::CAOKohnShamMatrix(CAOKohnShamMatrix&& source) noexcept

    : _xcMatrices(std::move(source._xcMatrices))

    , _xcRestricted(std::move(source._xcRestricted))

    , _xcElectrons(std::move(source._xcElectrons))

    , _xcEnergy(std::move(source._xcEnergy))
{
}

CAOKohnShamMatrix::~CAOKohnShamMatrix()
{
}

auto
CAOKohnShamMatrix::operator=(const CAOKohnShamMatrix& source) -> CAOKohnShamMatrix&
{
    if (this == &source) return *this;

    _xcMatrices = source._xcMatrices;

    _xcRestricted = source._xcRestricted;

    _xcElectrons = source._xcElectrons;

    _xcEnergy = source._xcEnergy;

    return *this;
}

auto
CAOKohnShamMatrix::operator=(CAOKohnShamMatrix&& source) noexcept -> CAOKohnShamMatrix&
{
    if (this == &source) return *this;

    _xcMatrices = std::move(source._xcMatrices);

    _xcRestricted = std::move(source._xcRestricted);

    _xcElectrons = std::move(source._xcElectrons);

    _xcEnergy = std::move(source._xcEnergy);

    return *this;
}

auto
CAOKohnShamMatrix::operator==(const CAOKohnShamMatrix& other) const -> bool
{
    if (_xcMatrices.size() != other._xcMatrices.size()) return false;

    for (size_t i = 0; i < _xcMatrices.size(); i++)
    {
        if (_xcMatrices[i] != other._xcMatrices[i]) return false;
    }

    if (_xcRestricted != other._xcRestricted) return false;

    if (std::fabs(_xcElectrons - other._xcElectrons) > 1.0e-13) return false;

    if (std::fabs(_xcEnergy - other._xcEnergy) > 1.0e-13) return false;

    return true;
}

auto
CAOKohnShamMatrix::operator!=(const CAOKohnShamMatrix& other) const -> bool
{
    return !(*this == other);
}

auto
CAOKohnShamMatrix::zero() -> void
{
    for (size_t i = 0; i < _xcMatrices.size(); i++)
    {
        _xcMatrices[i].zero();
    }
}

auto
CAOKohnShamMatrix::setNumberOfElectrons(const double xcElectrons) -> void
{
    _xcElectrons = xcElectrons;
}

auto
CAOKohnShamMatrix::setExchangeCorrelationEnergy(const double xcEnergy) -> void
{
    _xcEnergy = xcEnergy;
}

auto
CAOKohnShamMatrix::isRestricted() const -> bool
{
    return _xcRestricted;
}

auto
CAOKohnShamMatrix::getNumberOfElectrons() const -> double
{
    return _xcElectrons;
}

auto
CAOKohnShamMatrix::getExchangeCorrelationEnergy() const -> double
{
    return _xcEnergy;
}

auto
CAOKohnShamMatrix::getNumberOfRows() const -> int
{
    return _xcMatrices[0].getNumberOfRows();
}

auto
CAOKohnShamMatrix::getNumberOfColumns() const -> int
{
    return _xcMatrices[0].getNumberOfColumns();
}

auto
CAOKohnShamMatrix::getNumberOfElements() const -> int
{
    return _xcMatrices[0].getNumberOfElements();
}

auto
CAOKohnShamMatrix::alphaValues() const -> const double*
{
    return _xcMatrices[0].values();
}

auto
CAOKohnShamMatrix::alphaValues() -> double*
{
    return _xcMatrices[0].values();
}

auto
CAOKohnShamMatrix::betaValues() const -> const double*
{
    if (_xcRestricted) return _xcMatrices[0].values();

    return _xcMatrices[1].values();
}

auto
CAOKohnShamMatrix::betaValues() -> double*
{
    if (_xcRestricted) return _xcMatrices[0].values();

    return _xcMatrices[1].values();
}
