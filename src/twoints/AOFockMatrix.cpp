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

#include "AOFockMatrix.hpp"

#include <mpi.h>

#include <cmath>

#include "DenseLinearAlgebra.hpp"
#include "NumaPolicy.hpp"
#include "StringFormat.hpp"

CAOFockMatrix::CAOFockMatrix()
{
}

CAOFockMatrix::CAOFockMatrix(const std::vector<CDenseMatrix>& fockMatrices,
                             const std::vector<fockmat>&      fockTypes,
                             const std::vector<double>&       scaleFactors,
                             const std::vector<int32_t>&      idDensityMatrices)

    : _fockMatrices(fockMatrices)

    , _fockTypes(fockTypes)

    , _scaleFactors(scaleFactors)

    , _idDensityMatrices(idDensityMatrices)
{
}

CAOFockMatrix::CAOFockMatrix(const CAODensityMatrix& aoDensityMatrix)
{
    auto dmtyp = aoDensityMatrix.getDensityType();

    for (int32_t i = 0; i < aoDensityMatrix.getNumberOfDensityMatrices(); i++)
    {
        // set up dimensions of Fock matrix

        auto nrow = aoDensityMatrix.getNumberOfRows(i);

        auto ncol = aoDensityMatrix.getNumberOfColumns(i);

        // spin restricted closed-shell Hartree-Fock
        // spin restricted AO to MO transformation
        // spin restricted general Fock matrix

        if ((dmtyp == denmat::rest) || (dmtyp == denmat::rmoij) || (dmtyp == denmat::rgen))
        {
            _fockMatrices.push_back(CDenseMatrix(nrow, ncol));

            if (dmtyp == denmat::rest)
            {
                _fockTypes.push_back(fockmat::restjk);
            }
            else if (dmtyp == denmat::rmoij)
            {
                _fockTypes.push_back(fockmat::rgenk);
            }
            else if (dmtyp == denmat::rgen)
            {
                _fockTypes.push_back(fockmat::rgenjk);
            }

            _scaleFactors.push_back(1.0);

            _idDensityMatrices.push_back(i);
        }

        // spin unrestricted open-shell Hartree-Fock

        if (dmtyp == denmat::unrest)
        {
            _fockMatrices.push_back(CDenseMatrix(nrow, ncol));

            _fockMatrices.push_back(CDenseMatrix(nrow, ncol));

            _fockTypes.push_back(fockmat::unrestjk);

            _fockTypes.push_back(fockmat::unrestjk);

            _scaleFactors.push_back(1.0);

            _scaleFactors.push_back(1.0);

            _idDensityMatrices.push_back(i);

            _idDensityMatrices.push_back(i);
        }

        // FIX ME: Add restricted open-shell Hartree-Fock
    }
}

CAOFockMatrix::CAOFockMatrix(const CAOFockMatrix& source)

    : _fockMatrices(source._fockMatrices)

    , _fockTypes(source._fockTypes)

    , _scaleFactors(source._scaleFactors)

    , _idDensityMatrices(source._idDensityMatrices)
{
}

CAOFockMatrix::CAOFockMatrix(CAOFockMatrix&& source) noexcept

    : _fockMatrices(std::move(source._fockMatrices))

    , _fockTypes(std::move(source._fockTypes))

    , _scaleFactors(std::move(source._scaleFactors))

    , _idDensityMatrices(std::move(source._idDensityMatrices))
{
}

CAOFockMatrix::~CAOFockMatrix()
{
}

CAOFockMatrix&
CAOFockMatrix::operator=(const CAOFockMatrix& source)
{
    if (this == &source) return *this;

    _fockMatrices = source._fockMatrices;

    _fockTypes = source._fockTypes;

    _scaleFactors = source._scaleFactors;

    _idDensityMatrices = source._idDensityMatrices;

    return *this;
}

CAOFockMatrix&
CAOFockMatrix::operator=(CAOFockMatrix&& source) noexcept
{
    if (this == &source) return *this;

    _fockMatrices = std::move(source._fockMatrices);

    _fockTypes = std::move(source._fockTypes);

    _scaleFactors = std::move(source._scaleFactors);

    _idDensityMatrices = std::move(source._idDensityMatrices);

    return *this;
}

bool
CAOFockMatrix::operator==(const CAOFockMatrix& other) const
{
    if (_fockMatrices.size() != other._fockMatrices.size()) return false;

    for (size_t i = 0; i < _fockMatrices.size(); i++)
    {
        if (_fockMatrices[i] != other._fockMatrices[i]) return false;
    }

    if (_fockTypes.size() != other._fockTypes.size()) return false;

    for (size_t i = 0; i < _fockTypes.size(); i++)
    {
        if (_fockTypes[i] != other._fockTypes[i]) return false;
    }

    if (_scaleFactors.size() != other._scaleFactors.size()) return false;

    for (size_t i = 0; i < _scaleFactors.size(); i++)
    {
        if (std::fabs(_scaleFactors[i] - other._scaleFactors[i]) > 1.0e-13) return false;
    }

    if (_idDensityMatrices.size() != other._idDensityMatrices.size()) return false;

    for (size_t i = 0; i < _idDensityMatrices.size(); i++)
    {
        if (_idDensityMatrices[i] != other._idDensityMatrices[i]) return false;
    }

    return true;
}

bool
CAOFockMatrix::operator!=(const CAOFockMatrix& other) const
{
    return !(*this == other);
}

int32_t
CAOFockMatrix::_getNumberOfMatricesPerFock() const
{
    for (size_t i = 0; i < _fockMatrices.size(); i++)
    {
        if ((_fockTypes[i] == fockmat::restjk) ||

            (_fockTypes[i] == fockmat::restjkx) ||

            (_fockTypes[i] == fockmat::restj) ||

            (_fockTypes[i] == fockmat::restk) ||

            (_fockTypes[i] == fockmat::restkx) ||

            (_fockTypes[i] == fockmat::rgenjk) ||

            (_fockTypes[i] == fockmat::rgenjkx) ||

            (_fockTypes[i] == fockmat::rgenj) ||

            (_fockTypes[i] == fockmat::rgenk) ||

            (_fockTypes[i] == fockmat::rgenkx))
        {
            return 1;
        }

        if ((_fockTypes[i] == fockmat::unrestjk) ||

            (_fockTypes[i] == fockmat::unrestj) ||

            (_fockTypes[i] == fockmat::unrestjkx))
        {
            return 2;
        }
    }

    return 1;
}

bool
CAOFockMatrix::isClosedShell() const
{
    return (_getNumberOfMatricesPerFock() == 1);
}

int32_t
CAOFockMatrix::_getMatrixID(const int32_t iFockMatrix, const std::string& spin) const
{
    if (isClosedShell())
    {
        return iFockMatrix;
    }
    else
    {
        if (fstr::upcase(spin) == std::string("ALPHA"))
        {
            return 2 * iFockMatrix;
        }
        else
        {
            return 2 * iFockMatrix + 1;
        }
    }
}

void
CAOFockMatrix::setFockType(const fockmat& fockType, const int32_t iFockMatrix, const std::string& spin)
{
    auto fockMatrixID = _getMatrixID(iFockMatrix, spin);

    _fockTypes[fockMatrixID] = fockType;
}

void
CAOFockMatrix::setFockScaleFactor(const double factor, const int32_t iFockMatrix, const std::string& spin)
{
    auto fockMatrixID = _getMatrixID(iFockMatrix, spin);

    _scaleFactors[fockMatrixID] = factor;
}

void
CAOFockMatrix::zero()
{
    for (size_t i = 0; i < _fockMatrices.size(); i++)
    {
        _fockMatrices[i].zero();
    }
}

void
CAOFockMatrix::symmetrize()
{
    if (isClosedShell())
    {
        for (int32_t i = 0; i < getNumberOfFockMatrices(); i++)
        {
            if (isSymmetric(i))
            {
                _fockMatrices[i].symmetrize();
            }
        }
    }
    else
    {
        for (int32_t i = 0; i < getNumberOfFockMatrices(); i++)
        {
            if (isSymmetric(i))
            {
                _fockMatrices[2 * i].symmetrize();

                _fockMatrices[2 * i + 1].symmetrize();
            }
        }
    }

    // FIX ME: Add antisymmetric matrices
}

void
CAOFockMatrix::add(const CAOFockMatrix& source)
{
    if (isClosedShell())
    {
        for (int32_t i = 0; i < getNumberOfFockMatrices(); i++)
        {
            _fockMatrices[i] = denblas::addAB(_fockMatrices[i], source._fockMatrices[i], 1.0);
        }
    }
    else
    {
        for (int32_t i = 0; i < getNumberOfFockMatrices(); i++)
        {
            _fockMatrices[2 * i] = denblas::addAB(_fockMatrices[2 * i], source._fockMatrices[2 * i], 1.0);

            _fockMatrices[2 * i + 1] = denblas::addAB(_fockMatrices[2 * i + 1], source._fockMatrices[2 * i + 1], 1.0);
        }
    }
}

void
CAOFockMatrix::addCoreHamiltonian(const CKineticEnergyMatrix&    kineticEnergyMatrix,
                                  const CNuclearPotentialMatrix& nuclearPotentialMatrix,
                                  const int32_t                  iFockMatrix)
{
    // set up pointer to kinetic energy matrix

    auto pkin = kineticEnergyMatrix.values();

    // set up pointer to nuclear potential matrix

    auto pnucpot = nuclearPotentialMatrix.values();

    if (isClosedShell())
    {
        // set up pointer to Fock matrix

        auto pfock = _fockMatrices[iFockMatrix].values();

        // add core Hamiltonian contributions

        auto ndim = _fockMatrices[iFockMatrix].getNumberOfElements();

#pragma omp simd aligned(pfock, pkin, pnucpot : VLX_ALIGN)
        for (int32_t i = 0; i < ndim; i++)
        {
            pfock[i] += pkin[i] - pnucpot[i];
        }
    }
    else
    {
        // set up pointer to Fock matrix

        auto pfock_a = _fockMatrices[2 * iFockMatrix].values();

        auto pfock_b = _fockMatrices[2 * iFockMatrix + 1].values();

        // add core Hamiltonian contributions

        auto ndim = _fockMatrices[2 * iFockMatrix].getNumberOfElements();

#pragma omp simd aligned(pfock_a, pfock_b, pkin, pnucpot : VLX_ALIGN)
        for (int32_t i = 0; i < ndim; i++)
        {
            pfock_a[i] += pkin[i] - pnucpot[i];

            pfock_b[i] += pkin[i] - pnucpot[i];
        }
    }
}

void
CAOFockMatrix::addOneElectronMatrix(const CDenseMatrix& oneElectronMatrix, const int32_t iFockMatrix, const std::string& spin)
{
    auto fockMatrixID = _getMatrixID(iFockMatrix, spin);

    // set up pointer to one electron matrix

    auto pone = oneElectronMatrix.values();

    // set up pointer to Fock matrix

    auto pfock = _fockMatrices[fockMatrixID].values();

    // add one electron operator contribution contributions

    auto ndim = _fockMatrices[fockMatrixID].getNumberOfElements();

#pragma omp simd aligned(pfock, pone : VLX_ALIGN)
    for (int32_t i = 0; i < ndim; i++)
    {
        pfock[i] += pone[i];
    }
}

void
CAOFockMatrix::scale(const double factor, const int32_t iFockMatrix, const std::string& spin)
{
    auto fockMatrixID = _getMatrixID(iFockMatrix, spin);

    auto pfock = _fockMatrices[fockMatrixID].values();

    // add one electron operator contribution contributions

    auto ndim = _fockMatrices[fockMatrixID].getNumberOfElements();

#pragma omp simd aligned(pfock : VLX_ALIGN)
    for (int32_t i = 0; i < ndim; i++)
    {
        pfock[i] *= factor;
    }
}

void
CAOFockMatrix::reduce_sum(int32_t rank, int32_t nodes, MPI_Comm comm)
{
    for (size_t i = 0; i < _fockMatrices.size(); i++)
    {
        _fockMatrices[i].reduce_sum(rank, nodes, comm);

        MPI_Barrier(comm);
    }
}

int32_t
CAOFockMatrix::getNumberOfFockMatrices() const
{
    if (isClosedShell())
    {
        return static_cast<int32_t>(_fockMatrices.size());
    }
    else
    {
        return static_cast<int32_t>(_fockMatrices.size()) / 2;
    }
}

int32_t
CAOFockMatrix::getNumberOfRows(const int32_t iFockMatrix) const
{
    if (iFockMatrix < getNumberOfFockMatrices())
    {
        if (isClosedShell())
        {
            return _fockMatrices[iFockMatrix].getNumberOfRows();
        }
        else
        {
            return _fockMatrices[2 * iFockMatrix].getNumberOfRows();
        }
    }

    return 0;
}

int32_t
CAOFockMatrix::getNumberOfColumns(const int32_t iFockMatrix) const
{
    if (iFockMatrix < getNumberOfFockMatrices())
    {
        if (isClosedShell())
        {
            return _fockMatrices[iFockMatrix].getNumberOfColumns();
        }
        else
        {
            return _fockMatrices[2 * iFockMatrix].getNumberOfColumns();
        }
    }

    return 0;
}

int32_t
CAOFockMatrix::getNumberOfElements(const int32_t iFockMatrix) const
{
    if (iFockMatrix < getNumberOfFockMatrices())
    {
        if (isClosedShell())
        {
            return _fockMatrices[iFockMatrix].getNumberOfElements();
        }
        else
        {
            return _fockMatrices[2 * iFockMatrix].getNumberOfElements();
        }
    }

    return 0;
}

const double*
CAOFockMatrix::getFock(const int32_t iFockMatrix, const std::string& spin) const
{
    if (iFockMatrix < getNumberOfFockMatrices())
    {
        auto fockMatrixID = _getMatrixID(iFockMatrix, spin);

        return _fockMatrices[fockMatrixID].values();
    }

    return nullptr;
}

double*
CAOFockMatrix::getFock(const int32_t iFockMatrix, const std::string& spin)
{
    if (iFockMatrix < getNumberOfFockMatrices())
    {
        auto fockMatrixID = _getMatrixID(iFockMatrix, spin);

        return _fockMatrices[fockMatrixID].values();
    }

    return nullptr;
}

const CDenseMatrix&
CAOFockMatrix::getReferenceToFock(const int32_t iFockMatrix, const std::string& spin) const
{
    auto fockMatrixID = _getMatrixID(iFockMatrix, spin);

    return _fockMatrices[fockMatrixID];
}

fockmat
CAOFockMatrix::getFockType(const int32_t iFockMatrix, const std::string& spin) const
{
    auto fockMatrixID = _getMatrixID(iFockMatrix, spin);

    return _fockTypes[fockMatrixID];
}

double
CAOFockMatrix::getScaleFactor(const int32_t iFockMatrix, const std::string& spin) const
{
    if (iFockMatrix < getNumberOfFockMatrices())
    {
        auto fockMatrixID = _getMatrixID(iFockMatrix, spin);

        return _scaleFactors[fockMatrixID];
    }

    return 0.0;
}

int32_t
CAOFockMatrix::getDensityIdentifier(const int32_t iFockMatrix, const std::string& spin) const
{
    if (iFockMatrix < getNumberOfFockMatrices())
    {
        auto fockMatrixID = _getMatrixID(iFockMatrix, spin);

        return _idDensityMatrices[fockMatrixID];
    }

    return -1;
}

bool
CAOFockMatrix::isSymmetric(const int32_t iFockMatrix) const
{
    if (iFockMatrix < getNumberOfFockMatrices())
    {
        auto fockindex = iFockMatrix;

        if (!isClosedShell())
        {
            fockindex = 2 * iFockMatrix;
        }

        // check if Fock matrix is square

        if (_fockMatrices[fockindex].getNumberOfRows() != _fockMatrices[fockindex].getNumberOfColumns())
        {
            return false;
        }

        // determine symmetry by Fock matrix type

        auto fcktyp = _fockTypes[fockindex];

        if (fcktyp == fockmat::restjk) return true;

        if (fcktyp == fockmat::restjkx) return true;

        if (fcktyp == fockmat::restj) return true;

        if (fcktyp == fockmat::restk) return true;

        if (fcktyp == fockmat::restkx) return true;

        if (fcktyp == fockmat::unrestjk) return true;

        if (fcktyp == fockmat::unrestj) return true;

        if (fcktyp == fockmat::unrestjkx) return true;
    }

    return false;
}

double
CAOFockMatrix::getElectronicEnergy(const int32_t iFockMatrix, const CAODensityMatrix& aoDensityMatrix, const int32_t iDensityMatrix) const
{
    if ((iFockMatrix < getNumberOfFockMatrices()) && (iDensityMatrix < aoDensityMatrix.getNumberOfMatrices()))
    {
        if (isClosedShell())
        {
            return denblas::trace(_fockMatrices[iFockMatrix], aoDensityMatrix.getReferenceToDensity(iDensityMatrix));
        }
        else
        {
            auto ifock_a = 2 * iFockMatrix;

            auto ifock_b = 2 * iFockMatrix + 1;

            auto idensity_a = 2 * iDensityMatrix;

            auto idensity_b = 2 * iDensityMatrix + 1;

            auto e_a = 0.5 * denblas::trace(_fockMatrices[ifock_a], aoDensityMatrix.getReferenceToDensity(idensity_a));

            auto e_b = 0.5 * denblas::trace(_fockMatrices[ifock_b], aoDensityMatrix.getReferenceToDensity(idensity_b));

            return e_a + e_b;
        }
    }

    return 0.0;
}

std::string
CAOFockMatrix::getString() const
{
    std::string dmat_str;

    for (size_t i = 0; i < _fockMatrices.size(); i++)
    {
        dmat_str += "Fock Type: " + to_string(_fockTypes[i]) + "\n";

        dmat_str += "Density Identifier: " + std::to_string(_idDensityMatrices[i]) + "\n";

        dmat_str += _fockMatrices[i].getString();
    }

    return dmat_str;
}

std::ostream&
operator<<(std::ostream& output, const CAOFockMatrix& source)
{
    output << std::endl;

    output << "[CAOFockMatrix (Object):" << &source << "]" << std::endl;

    output << "_fockMatrices: " << std::endl;

    for (size_t i = 0; i < source._fockMatrices.size(); i++)
    {
        output << "_fockMatrices[" << i << "]: " << std::endl;

        output << source._fockMatrices[i] << std::endl;
    }

    output << "_fockTypes: " << std::endl;

    for (size_t i = 0; i < source._fockTypes.size(); i++)
    {
        output << "_fockTypes[" << i << "]: ";

        output << to_string(source._fockTypes[i]) << std::endl;
    }

    output << "_scaleFactors: " << std::endl;

    for (size_t i = 0; i < source._fockTypes.size(); i++)
    {
        output << "_scaleFactors[" << i << "]: ";

        output << source._scaleFactors[i] << std::endl;
    }

    output << "_idDensityMatrices: " << std::endl;

    for (size_t i = 0; i < source._idDensityMatrices.size(); i++)
    {
        output << "_idDensityMatrices[" << i << "]: ";

        output << source._idDensityMatrices[i] << std::endl;
    }

    return output;
}
