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

#include "AtomBasis.hpp"

#include <sstream>

#include <mpi.h>

#include "ChemicalElement.hpp"
#include "MpiFunc.hpp"
#include "StringFormat.hpp"

CAtomBasis::CAtomBasis()

    : _idElemental(-1)

    , _maxAngularMomentum(-1)
{
}

CAtomBasis::CAtomBasis(const CAtomBasis& source)

    : _basisFunctions(source._basisFunctions)

    , _idElemental(source._idElemental)

    , _maxAngularMomentum(source._maxAngularMomentum)
{
}

CAtomBasis::CAtomBasis(CAtomBasis&& source) noexcept

    : _basisFunctions(std::move(source._basisFunctions))

    , _idElemental(std::move(source._idElemental))

    , _maxAngularMomentum(std::move(source._maxAngularMomentum))
{
}

CAtomBasis::~CAtomBasis()
{
}

CAtomBasis&
CAtomBasis::operator=(const CAtomBasis& source)
{
    if (this == &source) return *this;

    _basisFunctions = source._basisFunctions;

    _idElemental = source._idElemental;

    _maxAngularMomentum = source._maxAngularMomentum;

    return *this;
}

CAtomBasis&
CAtomBasis::operator=(CAtomBasis&& source) noexcept
{
    if (this == &source) return *this;

    _basisFunctions = std::move(source._basisFunctions);

    _idElemental = std::move(source._idElemental);

    _maxAngularMomentum = std::move(source._maxAngularMomentum);

    return *this;
}

bool
CAtomBasis::operator==(const CAtomBasis& other) const
{
    if (_basisFunctions.size() != other._basisFunctions.size()) return false;

    for (size_t i = 0; i < _basisFunctions.size(); i++)
    {
        if (_basisFunctions[i] != other._basisFunctions[i]) return false;
    }

    if (_idElemental != other._idElemental) return false;

    if (_maxAngularMomentum != other._maxAngularMomentum) return false;

    return true;
}

bool
CAtomBasis::operator!=(const CAtomBasis& other) const
{
    return !(*this == other);
}

void
CAtomBasis::setIdElemental(const int32_t idElemental)
{
    _idElemental = idElemental;
}

void
CAtomBasis::setMaxAngularMomentum(const int32_t maxAngularMomentum)
{
    _maxAngularMomentum = maxAngularMomentum;
}

void
CAtomBasis::addBasisFunction(const CBasisFunction& basisFunction)
{
    _basisFunctions.push_back(basisFunction);

    auto bAngularMomentum = basisFunction.getAngularMomentum();

    if (bAngularMomentum > _maxAngularMomentum)
    {
        _maxAngularMomentum = bAngularMomentum;
    }
}

int32_t
CAtomBasis::getIdElemental() const
{
    return _idElemental;
}

int32_t
CAtomBasis::getMaxAngularMomentum() const
{
    return _maxAngularMomentum;
}

int32_t
CAtomBasis::getNumberOfBasisFunctions(const int32_t angularMomentum) const
{
    if (angularMomentum > _maxAngularMomentum) return 0;

    int32_t nbfuncs = 0;

    for (size_t i = 0; i < _basisFunctions.size(); i++)
    {
        if (_basisFunctions[i].getAngularMomentum() == angularMomentum)
        {
            nbfuncs++;
        }
    }

    return nbfuncs;
}

int32_t
CAtomBasis::getNumberOfBasisFunctions(const int32_t angularMomentum,
                                      const int32_t nPrimitiveGtos) const
{
    if (angularMomentum > _maxAngularMomentum) return 0;

    int32_t nbfuncs = 0;

    for (size_t i = 0; i < _basisFunctions.size(); i++)
    {
        if ((_basisFunctions[i].getAngularMomentum() == angularMomentum) &&
            (_basisFunctions[i].getNumberOfPrimitiveFunctions() == nPrimitiveGtos))
        {
            nbfuncs++;
        }
    }

    return nbfuncs;
}

int32_t
CAtomBasis::getNumberOfPrimitiveFunctions(const int32_t angularMomentum) const
{
    if (angularMomentum > _maxAngularMomentum) return 0;

    int32_t npfuncs = 0;

    for (size_t i = 0; i < _basisFunctions.size(); i++)
    {
        if (_basisFunctions[i].getAngularMomentum() == angularMomentum)
        {
            npfuncs += _basisFunctions[i].getNumberOfPrimitiveFunctions();
        }
    }

    return npfuncs;
}

std::set<int32_t>
CAtomBasis::getContractionDepths(const int32_t angularMomentum) const
{
    std::set<int32_t> cnums;
    
    for (const auto& bfunc : _basisFunctions)
    {
        if (bfunc.getAngularMomentum() == angularMomentum)
        {
            cnums.insert(bfunc.getNumberOfPrimitiveFunctions());
        }
    }
    
    return cnums;
}

std::string
CAtomBasis::getContractionString() const
{
    std::string str("(");

    for (int32_t i = 0; i <= _maxAngularMomentum; i++)
    {
        str.append(std::to_string(getNumberOfBasisFunctions(i)));

        str.append(fstr::to_AngularMomentum(i));

        if (i != _maxAngularMomentum) str.append(",");
    }

    str.append(")");

    return str;
}

std::string
CAtomBasis::getPrimitivesString() const
{
    std::string str("(");

    for (int32_t i = 0; i <= _maxAngularMomentum; i++)
    {
        str.append(std::to_string(getNumberOfPrimitiveFunctions(i)));

        str.append(fstr::to_AngularMomentum(i));

        if (i != _maxAngularMomentum) str.append(",");
    }

    str.append(")");

    return str;
}

std::vector<CBasisFunction>
CAtomBasis::getBasisFunctions(const int32_t angularMomentum) const
{
    std::vector<CBasisFunction> basvector;

    for (size_t i = 0; i < _basisFunctions.size(); i++)
    {
        if (_basisFunctions[i].getAngularMomentum() == angularMomentum)
        {
            basvector.push_back(_basisFunctions[i]);
        }
    }

    return basvector;
}

CAtomBasis
CAtomBasis::reduceToValenceBasis() const
{
    // set atomic shell max. angular momentum

    CChemicalElement chemele;

    chemele.setAtomType(_idElemental);

    auto mang = chemele.getMaxAngularMomentum();

    // generate valence basis

    CAtomBasis valbas;

    valbas.setIdElemental(_idElemental);

    for (size_t i = 0; i < _basisFunctions.size(); i++)
    {
        if (_basisFunctions[i].getAngularMomentum() <= mang)
        {
            valbas.addBasisFunction(_basisFunctions[i]);
        }
    }

    return valbas;
}

void
CAtomBasis::broadcast(int32_t rank, MPI_Comm comm)
{
    if constexpr (ENABLE_MPI)
    {
        mpi::bcast(_idElemental, comm);

        mpi::bcast(_maxAngularMomentum, comm);

        int32_t nbasfuncs = static_cast<int32_t>(_basisFunctions.size());

        mpi::bcast(nbasfuncs, comm);

        for (int32_t i = 0; i < nbasfuncs; i++)
        {
            CBasisFunction bfunc;

            if (rank == mpi::master()) bfunc = _basisFunctions[i];

            bfunc.broadcast(rank, comm);

            if (rank != mpi::master()) addBasisFunction(bfunc);

            MPI_Barrier(comm);
        }
    }
}

std::string CAtomBasis::repr() const {
    std::ostringstream os;

    os << std::endl;

    os << "[CAtomBasis (Object):" << this << "]" << std::endl;

    os << "_idElemental: " << _idElemental << std::endl;

    os << "_maxAngularMomentum: " << _maxAngularMomentum;

    os << std::endl;

    os << "_basisFunctions: " << std::endl;

    for (size_t i = 0; i < _basisFunctions.size(); i++)
    {
        os << "_basisFunctions[" << i << "]: " << std::endl;

        os << _basisFunctions[i] << std::endl;
    }

    return os.str();
}

std::ostream&
operator<<(std::ostream& output, const CAtomBasis& source)
{
    return (output << source.repr());
}
