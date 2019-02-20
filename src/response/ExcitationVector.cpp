//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "ExcitationVector.hpp"

#include <utility>

#include "DenseLinearAlgebra.hpp"

CExcitationVector::CExcitationVector()

    : _excitationType(szblock::aa)
{
    
}

CExcitationVector::CExcitationVector(const szblock               excitationType,
                                     const std::vector<int32_t>& braIndexes,
                                     const std::vector<int32_t>& ketIndexes,
                                     const std::vector<double>&  zCoefficients,
                                     const std::vector<double>&  yCoefficients)

    : _excitationType(excitationType)

    , _braIndexes(CMemBlock<int32_t>(braIndexes))

    , _ketIndexes(CMemBlock<int32_t>(ketIndexes))

    , _zCoefficents(CMemBlock<double>(zCoefficients))

    , _yCoefficents(CMemBlock<double>(yCoefficients))
{
    
}

CExcitationVector::CExcitationVector(const szblock excitationType,
                                     const int32_t braStartPosition,
                                     const int32_t braEndPosition,
                                     const int32_t ketStartPosition,
                                     const int32_t ketEndPosition)

    : _excitationType(excitationType)
{
    auto nbra = braEndPosition - braStartPosition;
    
    auto nket = ketEndPosition - ketStartPosition;
    
    _braIndexes = CMemBlock<int32_t>(nket * nbra);
    
    _ketIndexes = CMemBlock<int32_t>(nket * nbra);
    
    _zCoefficents = CMemBlock<double>(nket * nbra);
    
    _yCoefficents = CMemBlock<double>(nket * nbra);
    
    _zCoefficents.zero();
    
    _yCoefficents.zero();
    
    int32_t idx = 0;
    
    for (int32_t i = braStartPosition; i < braEndPosition; i++)
    {
        for (int32_t j = ketStartPosition; j < ketEndPosition; j++)
        {
            _braIndexes.at(idx) = i;
            
            _ketIndexes.at(idx) = j; 
            
            idx++;
        }
    }
}

CExcitationVector::CExcitationVector(const CExcitationVector& source)

    : _excitationType(source._excitationType)

    , _braIndexes(source._braIndexes)

    , _ketIndexes(source._ketIndexes)

    , _zCoefficents(source._zCoefficents)

    , _yCoefficents(source._yCoefficents)
{
    
}

CExcitationVector::CExcitationVector(CExcitationVector&& source) noexcept

    : _excitationType(std::move(source._excitationType))

    , _braIndexes(std::move(source._braIndexes))

    , _ketIndexes(std::move(source._ketIndexes))

    , _zCoefficents(std::move(source._zCoefficents))

    , _yCoefficents(std::move(source._yCoefficents))
{
    
}

CExcitationVector::~CExcitationVector()
{
    
}

CExcitationVector&
CExcitationVector::operator=(const CExcitationVector& source)
{
    if (this == &source) return *this;
    
    _excitationType = source._excitationType;
    
    _braIndexes = source._braIndexes;
    
    _ketIndexes = source._ketIndexes;
    
    _zCoefficents = source._zCoefficents;
    
    _yCoefficents = source._yCoefficents;
    
    return *this;
}

CExcitationVector&
CExcitationVector::operator=(CExcitationVector&& source) noexcept
{
    if (this == &source) return *this;
    
    _excitationType = std::move(source._excitationType);
    
    _braIndexes = std::move(source._braIndexes);
    
    _ketIndexes = std::move(source._ketIndexes);
    
    _zCoefficents = std::move(source._zCoefficents);
    
    _yCoefficents = std::move(source._yCoefficents);
    
    return *this;
}

bool
CExcitationVector::operator==(const CExcitationVector& other) const
{
    if (_excitationType != other._excitationType) return false;
    
    if (_braIndexes != other._braIndexes) return false;
    
    if (_ketIndexes != other._ketIndexes) return false;
    
    if (_zCoefficents != other._zCoefficents) return false;
    
    if (_yCoefficents != other._yCoefficents) return false;
    
    return true;
}

bool
CExcitationVector::operator!=(const CExcitationVector& other) const
{
    return !(*this == other);
}

void
CExcitationVector::setCoefficientsZY(const CMemBlock<double>& zCoefficients,
                                     const CMemBlock<double>& yCoefficients)
{
    _zCoefficents = zCoefficients;
    
    _yCoefficents = yCoefficients; 
}

double*
CExcitationVector::getCoefficientsZ()
{
    return _zCoefficents.data();
}

const double*
CExcitationVector::getCoefficientsZ() const
{
    return _zCoefficents.data();
}

double*
CExcitationVector::getCoefficientsY()
{
    return _yCoefficents.data();
}

const double*
CExcitationVector::getCoefficientsY() const
{
    return _yCoefficents.data();
}

int32_t
CExcitationVector::getNumberOfExcitations() const
{
    return _zCoefficents.size();
}

CAODensityMatrix
CExcitationVector::getTransformedDensity(const CMolecularOrbitals& molecularOrbitals) const
{
    auto ndim = molecularOrbitals.getNumberOfRows();
    
    CDenseMatrix mden(ndim, ndim);
    
    mden.zero();
    
    return CAODensityMatrix({mden}, denmat::rgen);
}


std::ostream&
operator<<(      std::ostream&  output,
           const CExcitationVector& source)
{
    output << std::endl;
    
    output << "[CExcitationVector (Object):" << &source << "]" << std::endl;
    
    output << "_excitationType: " << to_string(source._excitationType) << std::endl;
    
    output << "_braIndexes: " << source._braIndexes << std::endl;
    
    output << "_ketIndexes: " << source._ketIndexes <<  std::endl;
    
    output << "_zCoefficents: " << source._zCoefficents <<  std::endl;
    
    output << "_yCoefficents: " << source._yCoefficents <<  std::endl;
    
    return output;
}
