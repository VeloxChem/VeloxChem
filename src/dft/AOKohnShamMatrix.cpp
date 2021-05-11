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

#include <string>

#include <mpi.h>

#include "AOKohnShamMatrix.hpp"
#include "ErrorHandler.hpp"
#include "MpiFunc.hpp"

CAOKohnShamMatrix::CAOKohnShamMatrix()
    : _xcMatrices(std::vector<CDenseMatrix>())

    , _xcRestricted(true)

    , _xcElectrons(0.0)

    , _xcEnergy(0.0)
{
    
}

CAOKohnShamMatrix::CAOKohnShamMatrix(const std::vector<CDenseMatrix>& xcMatrices,
                                     const bool                       xcRestricted,
                                     const double                     xcElectrons,
                                     const double                     xcEnergy)

    : _xcMatrices(xcMatrices)

    , _xcRestricted(xcRestricted)

    , _xcElectrons(xcElectrons)

    , _xcEnergy(xcEnergy)
{
    
}

CAOKohnShamMatrix::CAOKohnShamMatrix(const int32_t nRows,
                                     const int32_t nColumns,
                                     const int32_t nMatrices, 
                                     const bool    xcRestricted)
{
    _xcRestricted = xcRestricted;
    
    _xcElectrons = 0.0;
    
    _xcEnergy = 0.0;
    
    for (int32_t i = 0; i < nMatrices; i++)
    {
        _xcMatrices.push_back(CDenseMatrix(nRows, nColumns));
    
        if (!_xcRestricted) _xcMatrices.push_back(CDenseMatrix(nRows, nColumns));
    }
}

CAOKohnShamMatrix::CAOKohnShamMatrix(const int32_t nRows,
                                     const int32_t nColumns,
                                     const bool    xcRestricted)

    : CAOKohnShamMatrix(nRows, nColumns, 1, xcRestricted)
{
    
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

CAOKohnShamMatrix&
CAOKohnShamMatrix::operator=(const CAOKohnShamMatrix& source)
{
    if (this == &source) return *this;
    
    _xcMatrices = source._xcMatrices;
    
    _xcRestricted = source._xcRestricted;
    
    _xcElectrons = source._xcElectrons;
    
    _xcEnergy = source._xcEnergy;
    
    return *this;
}

CAOKohnShamMatrix&
CAOKohnShamMatrix::operator=(CAOKohnShamMatrix&& source) noexcept
{
    if (this == &source) return *this;
    
    _xcMatrices = std::move(source._xcMatrices);
    
    _xcRestricted = std::move(source._xcRestricted);
    
    _xcElectrons = std::move(source._xcElectrons);
    
    _xcEnergy = std::move(source._xcEnergy);
    
    return *this;
}

bool
CAOKohnShamMatrix::operator==(const CAOKohnShamMatrix& other) const
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

bool
CAOKohnShamMatrix::operator!=(const CAOKohnShamMatrix& other) const
{
    return !(*this == other);
}

void
CAOKohnShamMatrix::zero()
{
    for (size_t i = 0; i < _xcMatrices.size(); i++)
    {
        _xcMatrices[i].zero();
    }
}

void
CAOKohnShamMatrix::symmetrize()
{
    for (size_t i = 0; i < _xcMatrices.size(); i++)
    {
        _xcMatrices[i].symmetrize();
    }
}

void
CAOKohnShamMatrix::setNumberOfElectrons(const double xcElectrons)
{
    _xcElectrons = xcElectrons;
}

void
CAOKohnShamMatrix::setExchangeCorrelationEnergy(const double xcEnergy)
{
    _xcEnergy = xcEnergy;
}

void
CAOKohnShamMatrix::reduce_sum(int32_t  rank,
                              int32_t  nodes,
                              MPI_Comm comm)
{
    for (size_t i = 0; i < _xcMatrices.size(); i++)
    {
        _xcMatrices[i].reduce_sum(rank, nodes, comm);
        
        MPI_Barrier(comm);
    }
    
    auto fsum = mpi::reduce_sum(_xcElectrons, comm);
    
    _xcElectrons = fsum;
    
    fsum = mpi::reduce_sum(_xcEnergy, comm);
    
    _xcEnergy = fsum; 
}

void
CAOKohnShamMatrix::collect(int32_t rank, int32_t nodes, MPI_Comm comm, int32_t source)
{
    if (ENABLE_MPI)
    {
        std::string errsource("AOKohnShamMatrix.collect: Invalid rank for the source");

        errors::assertMsgCritical(0 <= source && source < nodes, errsource);

        // master: receive data

        if (rank == mpi::master())
        {
            int32_t tag_id = source;

            MPI_Status mstat;

            // xc_integers: nxcmats, nrows, ncols, rest

            std::vector<int32_t> xc_integers(4);

            auto merror = MPI_Recv(xc_integers.data(), xc_integers.size(), MPI_INT, source, tag_id++, comm, &mstat);

            if (merror != MPI_SUCCESS) mpi::abort(merror, "collectKohnShamMatrix: xc_integers");

            int32_t nxcmats = xc_integers[0];

            int32_t nrows = xc_integers[1];

            int32_t ncols = xc_integers[2];

            _xcRestricted = (xc_integers[3] == 1) ? true : false;

            // xc_doubles: xc_electrons, xc_energy

            std::vector<double> xc_doubles(2);

            merror = MPI_Recv(xc_doubles.data(), xc_doubles.size(), MPI_DOUBLE, source, tag_id++, comm, &mstat);

            if (merror != MPI_SUCCESS) mpi::abort(merror, "collectKohnShamMatrix: xc_doubles");

            _xcElectrons = xc_doubles[0];

            _xcEnergy = xc_doubles[1];

            // _xcMatrices

            std::vector<double> data(nrows * ncols);

            for (int32_t imat = 0; imat < nxcmats; imat++)
            {
                merror = MPI_Recv(data.data(), nrows * ncols, MPI_DOUBLE, source, tag_id++, comm, &mstat);

                if (merror != MPI_SUCCESS) mpi::abort(merror, "collectKohnShamMatrix: xc_matrices");

                _xcMatrices.push_back(CDenseMatrix(data, nrows, ncols));
            }
        }

        // source: send data

        else if (rank == source)
        {
            int32_t tag_id = rank;

            // xc_integers: nxcmats, nrows, ncols, rest

            int32_t nxcmats = static_cast<int32_t>(_xcMatrices.size());

            int32_t nrows = getNumberOfRows();

            int32_t ncols = getNumberOfColumns();

            int32_t rest = _xcRestricted ? 1 : 0;

            std::vector<int32_t> xc_integers({nxcmats, nrows, ncols, rest});

            auto merror = MPI_Send(xc_integers.data(), xc_integers.size(), MPI_INT, mpi::master(), tag_id++, comm);

            if (merror != MPI_SUCCESS) mpi::abort(merror, "collectKohnShamMatrix: xc_integers");

            // xc_doubles: xc_electrons, xc_energy

            std::vector<double> xc_doubles({_xcElectrons, _xcEnergy});

            merror = MPI_Send(xc_doubles.data(), xc_doubles.size(), MPI_DOUBLE, mpi::master(), tag_id++, comm);

            if (merror != MPI_SUCCESS) mpi::abort(merror, "collectKohnShamMatrix: xc_doubles");

            // _xcMatrices

            for (int32_t imat = 0; imat < nxcmats; imat++)
            {
                merror = MPI_Send(_xcMatrices[imat].values(), nrows * ncols, MPI_DOUBLE, mpi::master(), tag_id++, comm);

                if (merror != MPI_SUCCESS) mpi::abort(merror, "collectKohnShamMatrix: xc_matrices");
            }
        }
    }
}

bool
CAOKohnShamMatrix::isRestricted() const
{
    return _xcRestricted;
}

double
CAOKohnShamMatrix::getNumberOfElectrons() const
{
    return _xcElectrons;
}

double
CAOKohnShamMatrix::getExchangeCorrelationEnergy() const
{
    return _xcEnergy; 
}

int32_t
CAOKohnShamMatrix::getNumberOfRows() const
{
    if (!_xcMatrices.empty())
    {
        return _xcMatrices[0].getNumberOfRows();
    }
    
    return 0;
}

int32_t
CAOKohnShamMatrix::getNumberOfColumns() const
{
    if (!_xcMatrices.empty())
    {
        return _xcMatrices[0].getNumberOfColumns();
    }
    
    return 0;
}

int32_t
CAOKohnShamMatrix::getNumberOfElements() const
{
    if (!_xcMatrices.empty())
    {
        return _xcMatrices[0].getNumberOfElements();
    }
    
    return 0;
}

int32_t
CAOKohnShamMatrix::getNumberOfMatrices() const
{
    return static_cast<int32_t>(_xcMatrices.size());
}

const double*
CAOKohnShamMatrix::getKohnSham(const bool beta) const
{
    if (!_xcMatrices.empty())
    {
        if (_xcRestricted)
        {
            return _xcMatrices[0].values();
        }
        else if (!beta)
        {
            return _xcMatrices[0].values();
        }
        else
        {
            return _xcMatrices[1].values();
        }
    }
    
    return nullptr;
}

double*
CAOKohnShamMatrix::getKohnSham(const bool beta)
{
    if (!_xcMatrices.empty())
    {
        if (_xcRestricted)
        {
            return _xcMatrices[0].values();
        }
        else if (!beta)
        {
            return _xcMatrices[0].values();
        }
        else
        {
            return _xcMatrices[1].values();
        }
    }
    
    return nullptr;
}

const CDenseMatrix&
CAOKohnShamMatrix::getReferenceToKohnSham(const bool beta) const
{
    if (_xcRestricted)
    {
        return _xcMatrices[0];
    }
    else if (!beta)
    {
        return _xcMatrices[0];
    }
    else
    {
        return _xcMatrices[1];
    }
}

const double*
CAOKohnShamMatrix::getMatrix(const int32_t iMatrix) const
{
    if (iMatrix < getNumberOfMatrices())
    {
        return _xcMatrices[iMatrix].values();
    }
    
    return nullptr;
}

double*
CAOKohnShamMatrix::getMatrix(const int32_t iMatrix)
{
    if (iMatrix < getNumberOfMatrices())
    {
        return _xcMatrices[iMatrix].values();
    }
    
    return nullptr;
}

const CDenseMatrix&
CAOKohnShamMatrix::getReferenceToMatrix(const int32_t iMatrix) const
{
    return _xcMatrices[iMatrix]; 
}

std::string
CAOKohnShamMatrix::getString() const
{
    if (_xcMatrices.empty()) return std::string();
 
    std::string ksmat_str;
    
    ksmat_str += "Is restricted: ";
    
    ksmat_str += (_xcRestricted) ? "Yes" : "No";

    ksmat_str += "Exchange-correlation energy: " + std::to_string(_xcEnergy) + "\n";
    
    ksmat_str += "Number of electrons: " + std::to_string(_xcElectrons) + "\n";
    
    if (_xcRestricted)
    {
        ksmat_str += "Total Kohn-Sham matrix: \n";
        
        ksmat_str += _xcMatrices[0].getString();
    }
    else
    {
        ksmat_str += "Alpha Kohn-Sham matrix: \n";
        
        ksmat_str += _xcMatrices[0].getString();
        
        ksmat_str += "Beta Kohn-Sham matrix: \n";
        
        ksmat_str += _xcMatrices[1].getString();
    }
    
    return ksmat_str;
}

std::ostream&
operator<<(      std::ostream&     output,
           const CAOKohnShamMatrix& source)
{
    output << std::endl;
    
    output << "[CAOKohnShamMatrix (Object):" << &source << "]" << std::endl;
    
    output << "_xcMatrices: " << std::endl;
    
    for (size_t i = 0; i < source._xcMatrices.size(); i++)
    {
        output << "_xcMatrices[" << i << "]: " << std::endl;
        
        output << source._xcMatrices[i] << std::endl;
    }
    
    output << "_xcRestricted: " << source._xcRestricted << std::endl;
    
    output << "_xcElectrons: " << source._xcElectrons << std::endl;
    
    output << "_xcEnergy" << source._xcEnergy;
    
    return output;
}
