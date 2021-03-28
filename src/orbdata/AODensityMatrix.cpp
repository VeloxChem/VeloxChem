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

#include "AODensityMatrix.hpp"

#include <mpi.h>

#include "DenseLinearAlgebra.hpp"
#include "ErrorHandler.hpp"
#include "MpiFunc.hpp"
#include "StringFormat.hpp"

CAODensityMatrix::CAODensityMatrix()

    : _denType(denmat::rest)
{
}

CAODensityMatrix::CAODensityMatrix(const std::vector<CDenseMatrix>& denMatrices, const denmat denType)

    : _denMatrices(denMatrices)

    , _denType(denType)
{
    if ((denType == denmat::unrest) || (denType == denmat::umoij) || (denType == denmat::osrest))
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

int32_t
CAODensityMatrix::_getNumberOfMatricesPerDensity() const
{
    if ((_denType == denmat::rest) || (_denType == denmat::rmoij) || (_denType == denmat::rgen))
    {
        return 1;
    }

    if ((_denType == denmat::unrest) || (_denType == denmat::umoij) || (_denType == denmat::osrest))
    {
        return 2;
    }

    return 1;
}

bool
CAODensityMatrix::isClosedShell() const
{
    return (_getNumberOfMatricesPerDensity() == 1);
}

int32_t
CAODensityMatrix::_getMatrixID(const int32_t iDensityMatrix, const std::string& spin) const
{
    if (isClosedShell())
    {
        return iDensityMatrix;
    }
    else
    {
        if (fstr::upcase(spin) == std::string("ALPHA"))
        {
            return 2 * iDensityMatrix;
        }
        else
        {
            return 2 * iDensityMatrix + 1;
        }
    }
}

void
CAODensityMatrix::setDensityType(const denmat denType)
{
    _denType = denType;
}

void
CAODensityMatrix::append(const CAODensityMatrix& other)
{
    if (_denType == other._denType)
    {
        for (size_t i = 0; i < other._denMatrices.size(); i++)
        {
            _denMatrices.push_back(other._denMatrices[i]);
        }
    }
}

CAODensityMatrix
CAODensityMatrix::sub(const CAODensityMatrix& other) const
{
    std::vector<CDenseMatrix> dmats;

    for (size_t i = 0; i < _denMatrices.size(); i++)
    {
        dmats.push_back(denblas::subAB(_denMatrices[i], other._denMatrices[i]));
    }

    return CAODensityMatrix(dmats, _denType);
}

int32_t
CAODensityMatrix::getNumberOfDensityMatrices() const
{
    if (isClosedShell())
    {
        return static_cast<int32_t>(_denMatrices.size());
    }
    else
    {
        return static_cast<int32_t>(_denMatrices.size()) / 2;
    }

    return 0;
}

int32_t
CAODensityMatrix::getNumberOfMatrices() const
{
    return static_cast<int32_t>(_denMatrices.size());
}

denmat
CAODensityMatrix::getDensityType() const
{
    return _denType;
}

int32_t
CAODensityMatrix::getNumberOfRows(const int32_t iDensityMatrix) const
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

int32_t
CAODensityMatrix::getNumberOfColumns(const int32_t iDensityMatrix) const
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

int32_t
CAODensityMatrix::getNumberOfElements(const int32_t iDensityMatrix) const
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
CAODensityMatrix::alphaDensity(const int32_t iDensityMatrix) const
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
CAODensityMatrix::betaDensity(const int32_t iDensityMatrix) const
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

const double*
CAODensityMatrix::getDensity(const int32_t iDensityMatrix) const
{
    if (iDensityMatrix < static_cast<int32_t>(_denMatrices.size()))
    {
        return _denMatrices[iDensityMatrix].values();
    }

    return nullptr;
}

const CDenseMatrix&
CAODensityMatrix::getReferenceToDensity(const int32_t iDensityMatrix) const
{
    return _denMatrices[iDensityMatrix];
}

std::string
CAODensityMatrix::getString() const
{
    std::string dmat_str;

    dmat_str += "Density Type: " + to_string(_denType) + "\n";

    for (size_t i = 0; i < _denMatrices.size(); i++)
    {
        dmat_str += _denMatrices[i].getString();
    }

    return dmat_str;
}

void
CAODensityMatrix::broadcast(int32_t rank, MPI_Comm comm)
{
    if (ENABLE_MPI)
    {
        // broadcast density matrix type

        int32_t dmtyp = 0;

        if (rank == mpi::master()) dmtyp = to_int(_denType);

        mpi::bcast(dmtyp, comm);

        if (rank != mpi::master()) _denType = to_denmat(dmtyp);

        // broadcast vector of AO density matrices

        int32_t ndmat = static_cast<int32_t>(_denMatrices.size());

        mpi::bcast(ndmat, comm);

        for (int32_t i = 0; i < ndmat; i++)
        {
            CDenseMatrix dmat;

            if (rank == mpi::master()) dmat = _denMatrices[i];

            dmat.broadcast(rank, comm);

            if (rank != mpi::master()) _denMatrices.push_back(dmat);

            MPI_Barrier(comm);
        }
    }
}

std::ostream&
operator<<(std::ostream& output, const CAODensityMatrix& source)
{
    output << std::endl;

    output << "[CAODensityMatrix (Object):" << &source << "]" << std::endl;

    output << "_denType: " << to_string(source._denType);

    output << "_denMatrices: " << std::endl;

    for (size_t i = 0; i < source._denMatrices.size(); i++)
    {
        output << "_denMatrices[" << i << "]: " << std::endl;

        output << source._denMatrices[i] << std::endl;
    }

    return output;
}
